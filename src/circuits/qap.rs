// How to turn an R1CS into a QAP and verify its satisfiability.
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;

use crate::utils::lagrange::compute_lagrange_interpolation_on_roots_of_unity;
use crate::utils::linear_algebra::Matrix;

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
        let lagrange_poly = compute_lagrange_interpolation_on_roots_of_unity(&evals);
        lagrange_polys.push(lagrange_poly);
    }
    lagrange_polys
}

#[cfg(test)]
pub mod tests {

    use crate::circuits::r1cs::{get_test_r1cs, get_test_satisfying_witness};
    use crate::utils::get_omega_domain;
    use crate::utils::linear_algebra::{Matrix, Vector};
    use ark_ff::One;
    use ark_ff::Zero;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::Polynomial;
    use ark_test_curves::bls12_381::Fr;

    use super::compute_lagrange_polynomial_from_matrix;

    #[test]
    pub fn test_qap_is_satisfied() {
        let (a, b, c): (Matrix<Fr>, Matrix<Fr>, Matrix<Fr>) = get_test_r1cs();
        let witness: Vector<Fr> = get_test_satisfying_witness(3);

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
            a_final_tampered_poly =
                &a_final_tampered_poly + &(&a_polys[i] * (witness.elements[i] + Fr::one()));
        }
        let final_tampered_poly: DensePolynomial<Fr> =
            &(&a_final_tampered_poly * &b_final_poly) - &c_final_poly;
        let (_, remainder) = final_tampered_poly
            .divide_by_vanishing_poly(domain)
            .unwrap();
        assert!(remainder.is_zero() == false);
    }
}
