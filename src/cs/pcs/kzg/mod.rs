use ark_ec::pairing::Pairing;
use ark_ff::{Field, One};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial};
use ark_std::Zero;

use crate::utils::build_zero_polynomial;

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

    pub fn multi_open(
        &self,
        polynomial: &DensePolynomial<E::ScalarField>,
        lagrange_polynomial: &DensePolynomial<E::ScalarField>,
        z_values: Vec<E::ScalarField>,
    ) -> E::G1 {
        let zero_polynomial = build_zero_polynomial::<E::ScalarField>(&z_values);
        let q = &(polynomial - lagrange_polynomial) / &zero_polynomial;
        let mut pi = self.g1 * E::ScalarField::ZERO;
        for (i, coeff) in q.coeffs.iter().enumerate() {
            pi += self.crs[i] * coeff;
        }
        pi
    }

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

    pub fn verify_multi_open(
        &self,
        commitment: E::G1,
        pi: E::G1,
        zero_polynomial: &DensePolynomial<E::ScalarField>,
        lagrange_polynomial: &DensePolynomial<E::ScalarField>,
    ) -> bool {
        let mut pz = self.g2 * E::ScalarField::ZERO;
        for (i, coeff) in zero_polynomial.coeffs.iter().enumerate() {
            pz += self.crs_2[i] * coeff;
        }

        let mut py = self.g1 * E::ScalarField::ZERO;
        for (i, coeff) in lagrange_polynomial.coeffs.iter().enumerate() {
            py += self.crs[i] * coeff;
        }

        let lhs = E::pairing(pi, pz);
        let rhs = E::pairing(commitment - py, self.g2);

        lhs == rhs
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bn254::{Bn254, Fr, G1Projective, G2Projective};
    use ark_ff::UniformRand;
    use ark_poly::{univariate::DensePolynomial, Polynomial};
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
}
