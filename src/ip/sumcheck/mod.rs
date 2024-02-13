use ark_ff::PrimeField;
use ark_poly::{
    multivariate::{SparsePolynomial, SparseTerm},
    Polynomial,
};

/// Utility types
pub type HyperCube<F> = Vec<Vec<F>>;

/// Utility to get all points on the hypercube
pub fn get_hypercube_points<F: PrimeField>(v: usize) -> HyperCube<F> {
    let mut points = Vec::with_capacity(1 << v);
    for i in 0..(1 << v) {
        let mut point = Vec::with_capacity(v);
        for j in 0..v {
            point.push(F::from(i >> j & 1u64));
        }
        points.push(point);
    }
    points
}

pub fn sample_random_vector<F: PrimeField>(v: usize) -> Vec<F> {
    let mut rng = ark_std::test_rng();
    let mut x = Vec::with_capacity(v);
    for _ in 0..v {
        x.push(F::rand(&mut rng));
    }
    x
}

pub fn get_evaluations_f_over_hypercube<F: PrimeField>(
    f: &SparsePolynomial<F, SparseTerm>,
    h: &HyperCube<F>,
) -> Vec<F> {
    let mut evaluations = Vec::with_capacity(h.len());
    for point in h {
        evaluations.push(f.evaluate(&point));
    }
    evaluations
}

pub fn compute_chi_w<F: PrimeField>(w: &Vec<F>, x: &Vec<F>) -> F {
    let mut chi_w = F::one();
    for (w_i, x_i) in w.iter().zip(x.iter()) {
        chi_w *= x_i.mul(w_i) + (F::one() - x_i) * (F::one() - w_i);
    }
    chi_w
}

/// Naive M.L.E. evaluations
/// Follows Thaler's notation in Proofs, Args and zk (lemma 3.6.) f, w, Chi, x
pub fn naive_mle_evaluation<F: PrimeField>(poly_evals: &Vec<F>, h: &HyperCube<F>, x: Vec<F>) -> F {
    let mut sum = F::zero();
    for (point, coeff) in h.iter().zip(poly_evals.iter()) {
        let chi_w = compute_chi_w::<F>(&point, &x);
        sum += *coeff * chi_w;
    }
    sum
}

pub fn binary_vec_to_usize<F: PrimeField>(binary_vec: &Vec<F>) -> usize {
    let mut result = 0;
    for (i, bit) in binary_vec.iter().enumerate() {
        if bit == &F::one() {
            result += 1 << i;
        }
    }
    result
}

/// Implementing lemma 3.8
/// At stage j, build table A^{j} of size 2^{j}
/// Follows Thaler's notation in Proofs, Args and zk (lemma 3.8.)
pub fn build_memoized_chi_table<F: PrimeField>(
    stage: usize,
    prev_table: Vec<F>,
    h: &HyperCube<F>,
    r: &Vec<F>,
) -> Vec<F> {
    let mut table = Vec::<F>::with_capacity(h.len());
    for w in h {
        let table_idx = binary_vec_to_usize::<F>(&w);
        let a_j = prev_table[table_idx]
            * (w[stage] * r[stage] + (F::one() - w[stage]) * (F::one() - r[stage]));
        table.push(a_j);
    }
    if stage == r.len() - 1 {
        return table;
    }
    build_memoized_chi_table(stage + 1, table, h, r)
}

pub fn memoized_mle_evaluation<F: PrimeField>(
    poly_evals: &Vec<F>,
    memoized_chi_table: &Vec<F>,
) -> F {
    let mut sum = F::zero();
    for (coeff, a_j) in poly_evals.iter().zip(memoized_chi_table.iter()) {
        sum += *coeff * a_j;
    }
    sum
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::Field;
    use ark_pallas::Fr;
    use ark_poly::{
        multivariate::{SparsePolynomial, SparseTerm},
        DenseMVPolynomial, DenseMultilinearExtension, MultilinearExtension,
    };
    use ark_std::test_rng;

    #[test]
    fn test_naive_mle() {
        let mut rng = test_rng();
        let n_vars = 5;
        let poly: SparsePolynomial<Fr, SparseTerm> = SparsePolynomial::rand(2, n_vars, &mut rng);
        let hypercube = get_hypercube_points::<Fr>(n_vars);
        let evaluations = get_evaluations_f_over_hypercube::<Fr>(&poly, &hypercube);
        let mle = DenseMultilinearExtension::from_evaluations_vec(n_vars, evaluations.clone());
        let x = sample_random_vector::<Fr>(n_vars);
        let naive_eval = naive_mle_evaluation::<Fr>(&evaluations, &hypercube, x.clone());
        let mle_eval = mle.evaluate(&x).unwrap();
        assert_eq!(naive_eval, mle_eval);
    }

    #[test]
    fn test_build_memoized_chi_table() {
        let n_vars = 5;
        let r = sample_random_vector::<Fr>(n_vars);
        let hypercube = get_hypercube_points::<Fr>(n_vars);
        let prev_table = vec![Fr::ONE; hypercube.len()];
        let table = build_memoized_chi_table::<Fr>(0, prev_table, &hypercube, &r);
        for (i, w) in hypercube.iter().enumerate() {
            let chi_w = compute_chi_w::<Fr>(&w, &r);
            assert_eq!(chi_w, table[i]);
        }
    }

    #[test]
    fn test_memoized_mle_evaluation() {
        let mut rng = test_rng();
        let n_vars = 5;
        let poly: SparsePolynomial<Fr, SparseTerm> = SparsePolynomial::rand(2, n_vars, &mut rng);
        let hypercube = get_hypercube_points::<Fr>(n_vars);
        let evaluations = get_evaluations_f_over_hypercube::<Fr>(&poly, &hypercube);
        let mle = DenseMultilinearExtension::from_evaluations_vec(n_vars, evaluations.clone());
        let x = sample_random_vector::<Fr>(n_vars);
        let chi_table =
            build_memoized_chi_table::<Fr>(0, vec![Fr::ONE; hypercube.len()], &hypercube, &x);
        let memoized_eval = memoized_mle_evaluation::<Fr>(&evaluations, &chi_table);
        let mle_eval = mle.evaluate(&x).unwrap();
        assert_eq!(memoized_eval, mle_eval);
    }
}
