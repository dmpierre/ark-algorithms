// Pedersen commitments
#[cfg(test)]
mod test {
    use ark_ff::UniformRand;
    use ark_pallas::Affine;
    use ark_pallas::Fr;
    use ark_std::test_rng;

    #[test]
    pub fn test_pedersen_commitment() {
        let mut rng = test_rng();
        let g = Affine::rand(&mut rng);
        let h = Affine::rand(&mut rng);

        let m_1 = Fr::rand(&mut rng);
        let m_2 = Fr::rand(&mut rng);
        let r_1 = Fr::rand(&mut rng);
        let r_2 = Fr::rand(&mut rng);

        let c_1 = g * m_1 + h * r_1;
        let c_2 = g * m_2 + h * r_2;

        // check the homomorphism
        let c1_plus_c2 = c_1 + c_2;
        let homomorphic_sum = g * (m_1 + m_2) + h * (r_1 + r_2);
        assert!(c1_plus_c2.eq(&homomorphic_sum));
    }
}
