use ark_ec::pairing::Pairing;
use ark_ff::{Field, One};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, Polynomial};
use ark_std::Zero;

use crate::utils::{build_zero_polynomial, lagrange::compute_lagrange_interpolation};

pub struct KZG<E: Pairing> {
    pub g1: E::G1,
    pub g2: E::G2,
    pub degree: usize,
    pub crs: Vec<E::G1>,
    pub crs_2: Vec<E::G2>,
    pub vk: E::G2,
}

impl<E: Pairing> KZG<E> {
    pub fn new(g1: E::G1, g2: E::G2, degree: usize) -> Self {
        Self {
            g1,
            g2,
            degree,
            crs: vec![],
            crs_2: vec![],
            vk: g2,
        }
    }

    pub fn setup(&mut self, tau: E::ScalarField) {
        let vk = self.g2 * tau;
        for pow in 0..self.degree + 1 {
            let tau_i: E::ScalarField = tau.pow([pow as u64]);
            let crs_point_g1 = self.g1 * tau_i;
            let crs_point_g2 = self.g2 * tau_i;
            self.crs.push(crs_point_g1);
            self.crs_2.push(crs_point_g2);
        }
        self.vk = vk;
    }

    pub fn commit(&mut self, polynomial: &DensePolynomial<E::ScalarField>) -> E::G1 {
        let mut commitment = E::G1::zero();
        for i in 0..self.degree + 1 {
            let value = self.crs[i as usize] * polynomial.coeffs[i as usize];
            commitment += value;
        }
        commitment
    }

    /// Single point kzg opening
    pub fn open(
        &self,
        polynomial: &DensePolynomial<E::ScalarField>,
        z: E::ScalarField,
        y: E::ScalarField,
    ) -> E::G1 {
        // Opening at y = p(z). Notation from here: https://hackmd.io/@gnark/kzg-bls24
        let y_polynomial = DensePolynomial::from_coefficients_vec(vec![y]);
        let numerator = polynomial - &y_polynomial;
        let denominator = DensePolynomial::from_coefficients_vec(vec![-z, E::ScalarField::ONE]);
        let q_x = &numerator / &denominator;
        let mut pi = E::G1::zero();
        for (i, coeff) in q_x.coeffs.iter().enumerate() {
            pi += self.crs[i] * coeff;
        }
        pi
    }

    /// Multi-point kzg opening, also referred as "batch opening"
    pub fn multi_open(
        &self,
        polynomial: &DensePolynomial<E::ScalarField>,
        z_values: &Vec<E::ScalarField>, // z_values = {0, 1, 2, 3, ...}
    ) -> (
        E::G2,
        DensePolynomial<E::ScalarField>,
        DensePolynomial<E::ScalarField>,
    ) {
        let mut y_values = Vec::new();
        for z in z_values.iter() {
            let y = polynomial.evaluate(z);
            y_values.push(y);
        }
        let lagrange_polynomial = compute_lagrange_interpolation::<E::ScalarField>(&y_values);
        let zero_polynomial = build_zero_polynomial::<E::ScalarField>(&z_values);
        let q = &(polynomial - &lagrange_polynomial) / &zero_polynomial;
        let mut pi = self.g2 * E::ScalarField::ZERO;
        for (i, coeff) in q.coeffs.iter().enumerate() {
            pi += self.crs_2[i] * coeff;
        }
        (pi, lagrange_polynomial, zero_polynomial)
    }

    /// Single point kzg verification
    pub fn verify(
        &self,
        y: E::ScalarField,
        z: E::ScalarField,
        commitment: E::G1,
        pi: E::G1,
    ) -> bool {
        let py = self.g1 * y;
        let pz = self.g2 * z;
        let lhs = E::pairing(pi, self.vk - pz);
        let rhs = E::pairing(commitment - py, self.g2);
        lhs == rhs
    }

    /// This is the same as `verify` but re-wrote as to avoid any operations in G2
    /// This is useful for testing the EVM implementation.
    pub fn verify_no_g2_ops(
        &self,
        y: E::ScalarField,
        z: E::ScalarField,
        commitment: E::G1,
        pi: E::G1,
    ) -> bool {
        let py = self.g1 * y;
        let g2_neg = -self.g2;
        let lhs_1 = E::pairing(pi, self.vk);
        let lhs_2 = E::pairing(pi * z, g2_neg);
        let lhs = lhs_1.0 * lhs_2.0;
        let rhs = E::pairing(commitment - py, self.g2);
        lhs == rhs.0
    }

    /// This is the same as `verify_no_g2_ops` but with the pairing written as an EVM opcode.
    /// This is useful for testing the EVM implementation.
    pub fn verify_no_g2_ops_evm_opcode(
        &self,
        y: E::ScalarField,
        z: E::ScalarField,
        commitment: E::G1,
        pi: E::G1,
    ) -> bool {
        let py = self.g1 * y;
        let lhs = E::pairing(pi, self.vk);
        let rhs = E::pairing(pi * -z - commitment + py, self.g2);
        (lhs.0 * rhs.0).is_one()
    }

    pub fn verify_from_encrypted_y(
        &self,
        py: E::G1,
        z: E::ScalarField,
        commitment: E::G1,
        pi: E::G1,
    ) -> bool {
        let pz = self.g2 * z;
        let lhs = E::pairing(pi, self.vk - pz);
        let rhs = E::pairing(commitment - py, self.g2);
        lhs == rhs
    }

    /// Verify a multi-open proof for a polynomial `p` at points `z_values`.
    /// Verifier can compute the zero polynomial on its own.
    /// The zero poly is provided here since we emulate what would be done within the EVM.
    /// This is because computing the zero poly would be quite expensive on the EVM.
    pub fn verify_multi_open_no_g2_ops(
        &self,
        commitment: &E::G1,
        z_values: &Vec<E::ScalarField>,
        y_values: &Vec<E::ScalarField>, // evaluations of \phi(z)
        lagrange_polynomial: &DensePolynomial<E::ScalarField>,
        zero_polynomial: &DensePolynomial<E::ScalarField>,
        pi: &E::G2,
    ) -> bool {
        // 1. check that lagrange interpolated poly is correct
        let _ = z_values
            .iter()
            .zip(y_values)
            .map(|(z, y)| assert_eq!(lagrange_polynomial.evaluate(&z), *y));

        // 2. check that the zero polynomial is zero at all z_values
        let _ = z_values
            .iter()
            .map(|z| assert_eq!(zero_polynomial.evaluate(&z), E::ScalarField::ZERO));

        // 3. Compute input values to pairing
        let z_tau = zero_polynomial
            .coeffs
            .iter()
            .zip(&self.crs)
            .fold(E::G1::zero(), |acc, (coeff, tau)| acc + *tau * coeff);

        let i_tau = lagrange_polynomial
            .coeffs
            .iter()
            .zip(&self.crs)
            .fold(E::G1::zero(), |acc, (coeff, tau)| acc + *tau * coeff);

        (E::pairing(z_tau, pi).0 * E::pairing(-*commitment + i_tau, self.g2).0).is_one()
    }
}

#[cfg(test)]
mod tests {
    use crate::cs::pcs::kzg::KZG;
    use ark_bn254::{Bn254, Fr, G1Projective, G2Projective};
    use ark_ff::UniformRand;
    use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, Polynomial};
    use ark_std::test_rng;

    #[test]
    pub fn test_full_kzg() {
        let mut rng = test_rng();
        let degree = 10;
        let tau = Fr::rand(&mut rng);
        let g1 = G1Projective::rand(&mut rng);
        let g2 = G2Projective::rand(&mut rng);
        let mut kzg = KZG::<Bn254>::new(g1, g2, degree);
        let polynomial: DensePolynomial<Fr> = DensePolynomial::rand(degree, &mut rng);
        let _ = kzg.setup(tau);
        let commitment = kzg.commit(&polynomial);
        let z = Fr::rand(&mut rng);
        let y = polynomial.evaluate(&z);
        let pi = kzg.open(&polynomial, z, y);
        assert!(kzg.verify(y, z, commitment, pi));
        assert!(kzg.verify_no_g2_ops(y, z, commitment, pi));
        assert!(kzg.verify_no_g2_ops_evm_opcode(y, z, commitment, pi));
    }

    #[test]
    pub fn test_multi_open_kzg_with_no_g2_ops() {
        let mut rng = test_rng();
        let degree = 10;
        let tau = Fr::rand(&mut rng);
        let g1 = G1Projective::rand(&mut rng);
        let g2 = G2Projective::rand(&mut rng);
        let mut kzg = KZG::<Bn254>::new(g1, g2, degree);
        let polynomial: DensePolynomial<Fr> = DensePolynomial::rand(degree, &mut rng);
        let _ = kzg.setup(tau);
        let commitment = kzg.commit(&polynomial);
        let z_values = (0..5).map(|i| Fr::from(i as u64)).collect::<Vec<_>>();
        let y_values = z_values
            .iter()
            .map(|z| polynomial.evaluate(z))
            .collect::<Vec<_>>();
        let (pi, lagrange_polynomial, zero_polynomial) = kzg.multi_open(&polynomial, &z_values);
        let result = kzg.verify_multi_open_no_g2_ops(
            &commitment,
            &z_values,
            &y_values,
            &lagrange_polynomial,
            &zero_polynomial,
            &pi,
        );
        assert!(result);

        // make it fail with wrong commitment
        let wrong_commitment = commitment + kzg.g1;
        let result = kzg.verify_multi_open_no_g2_ops(
            &wrong_commitment,
            &z_values,
            &y_values,
            &lagrange_polynomial,
            &zero_polynomial,
            &pi,
        );
        assert!(!result);
    }
}
