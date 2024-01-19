use ark_ff::PrimeField;
use ark_poly::{
    univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Evaluations,
    GeneralEvaluationDomain,
};

/// Computes the lagrange interpolation for the set of points:
/// (\omega^{0}, y_0), (\omega^{1}, y_1), ..., (\omega^{n}, y_n)
/// where \omega is a primitive n-th root of unity.
pub fn compute_lagrange_interpolation_on_roots_of_unity<F: PrimeField>(
    evals: &Vec<F>,
) -> DensePolynomial<F> {
    let k = evals.len() as usize;
    let omegas = GeneralEvaluationDomain::<F>::new(k).unwrap();
    let lagrange: DensePolynomial<F> =
        Evaluations::<F>::from_vec_and_domain(evals.to_vec(), omegas).interpolate();
    lagrange
}

/// Computes the lagrange interpolated polynomial from the given points `p_i`
/// From https://github.com/privacy-scaling-explorations/folding-schemes
pub fn compute_lagrange_interpolation<F: PrimeField>(p_i: &[F]) -> DensePolynomial<F> {
    // domain is 0..p_i.len(), to fit `interpolate_uni_poly` from hyperplonk
    let domain: Vec<usize> = (0..p_i.len()).collect();

    // compute l(x), common to every basis polynomial
    let mut l_x = DensePolynomial::from_coefficients_vec(vec![F::ONE]);
    for x_m in domain.clone() {
        let prod_m = DensePolynomial::from_coefficients_vec(vec![-F::from(x_m as u64), F::ONE]);
        l_x = &l_x * &prod_m;
    }

    // compute each w_j - barycentric weights
    let mut w_j_vector: Vec<F> = vec![];
    for x_j in domain.clone() {
        let mut w_j = F::ONE;
        for x_m in domain.clone() {
            if x_m != x_j {
                let prod = (F::from(x_j as u64) - F::from(x_m as u64))
                    .inverse()
                    .unwrap(); // an inverse always exists since x_j != x_m (!=0)
                               // hence, we call unwrap() here without checking the Option's content
                w_j *= prod;
            }
        }
        w_j_vector.push(w_j);
    }

    // compute each polynomial within the sum L(x)
    let mut lagrange_poly = DensePolynomial::from_coefficients_vec(vec![F::ZERO]);
    for (j, w_j) in w_j_vector.iter().enumerate() {
        let x_j = domain[j];
        let y_j = p_i[j];
        // we multiply by l(x) here, otherwise the below division will not work - deg(0)/deg(d)
        let poly_numerator = &(&l_x * (*w_j)) * (y_j);
        let poly_denominator =
            DensePolynomial::from_coefficients_vec(vec![-F::from(x_j as u64), F::ONE]);
        let poly = &poly_numerator / &poly_denominator;
        lagrange_poly = &lagrange_poly + &poly;
    }

    lagrange_poly
}
