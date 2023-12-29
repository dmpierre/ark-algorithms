// How to turn an R1CS into a QAP and verify its satisfiability.
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};

use crate::utils::lagrange::compute_lagrange_interpolation;
use crate::utils::linear_algebra::Matrix;

pub fn get_omega_domain<F: PrimeField>(n: usize) -> (GeneralEvaluationDomain<F>, Vec<F>) {
    // Builds the domain consisting of n roots of unity in F
    let omegas = GeneralEvaluationDomain::<F>::new(n).unwrap();
    let mut domain_elements: Vec<F> = vec![];
    for element in omegas.elements() {
        domain_elements.push(element);
    }
    (omegas, domain_elements)
}

pub fn compute_lagrange_polynomial_from_matrix<F: PrimeField>(
    mat: &Matrix<F>,
) -> Vec<DensePolynomial<F>> {
    let mut lagrange_polys: Vec<DensePolynomial<F>> = Vec::with_capacity(mat.num_cols);
    let n_cols = mat.num_cols;
    for i in 0..n_cols {
        let mut evals: Vec<F> = Vec::with_capacity(mat.num_rows);
        for j in 0..mat.num_rows {
            evals.push(mat.rows[j].elements[i]);
        }
        // lagrange polynomial for the i-th column
        let lagrange_poly = compute_lagrange_interpolation(&evals);
        lagrange_polys.push(lagrange_poly);
    }
    lagrange_polys
}

#[cfg(test)]
pub mod tests {

    use crate::qap::get_omega_domain;
    use crate::utils::linear_algebra::{Matrix, Vector};
    use ark_ff::One;
    use ark_ff::PrimeField;
    use ark_ff::Zero;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::Polynomial;
    use ark_test_curves::bls12_381::Fr;

    use super::compute_lagrange_polynomial_from_matrix;

    pub fn get_test_r1cs<F: PrimeField>() -> (Matrix<F>, Matrix<F>, Matrix<F>) {
        // Taken from vb: https://medium.com/@VitalikButerin/quadratic-arithmetic-programs-from-zero-to-hero-f6d558cea649
        let a: Vec<Vec<F>> = vec![
            vec![
                F::from(0 as u8),
                F::from(1 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
            ],
            vec![
                F::from(0 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
                F::from(1 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
            ],
            vec![
                F::from(0 as u8),
                F::from(1 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
                F::from(1 as u8),
                F::from(0 as u8),
            ],
            vec![
                F::from(5 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
                F::from(1 as u8),
            ],
        ];
        let b: Vec<Vec<F>> = vec![
            vec![
                F::from(0 as u8),
                F::from(1 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
            ],
            vec![
                F::from(0 as u8),
                F::from(1 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
            ],
            vec![
                F::from(1 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
            ],
            vec![
                F::from(1 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
            ],
        ];
        let c: Vec<Vec<F>> = vec![
            vec![
                F::from(0 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
                F::from(1 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
            ],
            vec![
                F::from(0 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
                F::from(1 as u8),
                F::from(0 as u8),
            ],
            vec![
                F::from(0 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
                F::from(1 as u8),
            ],
            vec![
                F::from(0 as u8),
                F::from(0 as u8),
                F::from(1 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
                F::from(0 as u8),
            ],
        ];
        (
            Matrix::new_from_vecs(&a),
            Matrix::new_from_vecs(&b),
            Matrix::new_from_vecs(&c),
        )
    }
    pub fn get_test_satisfying_witness<F: PrimeField>() -> Vector<F> {
        // Taken from vb: https://medium.com/@VitalikButerin/quadratic-arithmetic-programs-from-zero-to-hero-f6d558cea649
        Vector::new(&vec![
            F::from(1 as u8),
            F::from(3 as u8),
            F::from(35 as u8),
            F::from(9 as u8),
            F::from(27 as u8),
            F::from(30 as u8),
        ])
    }

    #[test]
    pub fn test_r1cs_is_satisfied() {
        let (a, b, c): (Matrix<Fr>, Matrix<Fr>, Matrix<Fr>) = get_test_r1cs();
        let witness: Vector<Fr> = get_test_satisfying_witness();
        let a_dot_w = a.dot_vector(&witness);
        let b_dot_w = b.dot_vector(&witness);
        let c_dot_w = c.dot_vector(&witness);
        let a_times_b = a_dot_w * b_dot_w;
        let a_times_b_minus_c = a_times_b - c_dot_w;
        assert!(a_times_b_minus_c.is_zero_vector());
    }

    #[test]
    pub fn test_qap_is_satisfied() {
        let (a, b, c): (Matrix<Fr>, Matrix<Fr>, Matrix<Fr>) = get_test_r1cs();
        let witness: Vector<Fr> = get_test_satisfying_witness();
        
        // we lagrange-interpolate polynomials over an n-roots of unity domain
        // i.e.: f(\omega^{i}) == vec[i]
        let a_polys = compute_lagrange_polynomial_from_matrix(&a);
        let b_polys = compute_lagrange_polynomial_from_matrix(&b);
        let c_polys = compute_lagrange_polynomial_from_matrix(&c);

        // retrieving the domain over which the polynomials have been interpolated (the roots of unity)
        let (domain, omegas) = get_omega_domain::<Fr>(a_polys[0].coeffs.len());
        
        // to illustrate, we can retrieve the last row of A:
        assert_eq!(a_polys[0].evaluate(&omegas[3]), a.rows[3].elements[0]);
        assert_eq!(a_polys[1].evaluate(&omegas[3]), a.rows[3].elements[1]);
        assert_eq!(a_polys[2].evaluate(&omegas[3]), a.rows[3].elements[2]);
        assert_eq!(a_polys[3].evaluate(&omegas[3]), a.rows[3].elements[3]);
        assert_eq!(a_polys[4].evaluate(&omegas[3]), a.rows[3].elements[4]);
        assert_eq!(a_polys[5].evaluate(&omegas[3]), a.rows[3].elements[5]);

        // compute a*s, b*s, c*s
        let mut a_final_poly: DensePolynomial<Fr> = DensePolynomial::zero();
        let mut b_final_poly: DensePolynomial<Fr> = DensePolynomial::zero();
        let mut c_final_poly: DensePolynomial<Fr> = DensePolynomial::zero();
        for i in 0..a_polys.len() {
            a_final_poly = &a_final_poly + &(&a_polys[i] * witness.elements[i]);
            b_final_poly = &b_final_poly + &(&b_polys[i] * witness.elements[i]);
            c_final_poly = &c_final_poly + &(&c_polys[i] * witness.elements[i]);
        }

        // compute (a*s) * (b*s) - (c*s)
        let final_poly: DensePolynomial<Fr> = &(&a_final_poly * &b_final_poly) - &c_final_poly;

        // check that the division by the vanishing polynomial leaves no remainder
        let (_, remainder) = final_poly.divide_by_vanishing_poly(domain).unwrap();
        assert!(remainder.is_zero());

        // tamper with the witness
        // and check that the division by the vanishing polynomial leaves a remainder
        let mut a_final_tampered_poly: DensePolynomial<Fr> = DensePolynomial::zero();
        for i in 0..a_polys.len() {
            a_final_tampered_poly = &a_final_tampered_poly + &(&a_polys[i] * (witness.elements[i] + Fr::one()));
        }
        let final_tampered_poly: DensePolynomial<Fr> = &(&a_final_tampered_poly * &b_final_poly) - &c_final_poly;
        let (_, remainder) = final_tampered_poly.divide_by_vanishing_poly(domain).unwrap();
        assert!(remainder.is_zero() == false);

    }
}
