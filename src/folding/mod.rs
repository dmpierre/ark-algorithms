#[cfg(test)]
mod tests {
    use crate::circuits::r1cs::utils::{get_test_r1cs, get_test_satisfying_witness};
    use crate::utils::linear_algebra::Vector;
    use ark_ff::Field;
    use ark_ff::UniformRand;
    use ark_std::test_rng;
    use ark_test_curves::bls12_381::Fr;

    #[test]
    pub fn test_nolib_folding_relaxed_r1cs() {
        // this is "raw" folding, without any external libraries
        // intent is to show how to fold relaxed r1cs
        let mut rng = test_rng();

        // check the satisfiability of folded relaxed r1cs using folded witness
        let (a, b, c) = get_test_r1cs::<Fr>();

        let w_1 = get_test_satisfying_witness::<Fr>(3);
        let w_2 = get_test_satisfying_witness::<Fr>(5);
        let (u_1, u_2) = (Fr::ONE, Fr::ONE);
        let e_1 = Vector::<Fr>::new_zero_vector(a.num_rows); // a.dot(w) has dimensions a.num_rows. see below.
        let e_2 = e_1.clone();

        // u <-- u_1 + r * u_2
        // u is a scalar
        let r = Fr::rand(&mut rng);
        let u = u_1 + r * u_2;

        // E <-- E_1 + r * (AZ_1 o BZ_2 + AZ_2 o BZ_1 - u_1CZ_2 - u_2CZ_1) + r^2 * E_2
        // E is the F^m error vector
        let e = (a.dot_vector(&w_1) * b.dot_vector(&w_2) + a.dot_vector(&w_2) * b.dot_vector(&w_1)
            - c.dot_vector(&w_2).scalar_mul(&u_1)
            - c.dot_vector(&w_1).scalar_mul(&u_2))
        .scalar_mul(&r)
            + e_1
            + e_2.scalar_mul(&(r.square()));

        // AZ o BZ
        let az_bz = a.dot_vector(&w_1) * b.dot_vector(&w_1)
            + (a.dot_vector(&w_1) * b.dot_vector(&w_2) + a.dot_vector(&w_2) * b.dot_vector(&w_1))
                .scalar_mul(&r)
            + (a.dot_vector(&w_2) * b.dot_vector(&w_2)).scalar_mul(&(r.square()));

        // uCZ + E
        let cz = c.dot_vector(&w_1) + (c.dot_vector(&w_2)).scalar_mul(&r);
        let u_cz_plus_e = cz.scalar_mul(&u) + e;

        // checks that relaxed r1cs is satisfied
        assert!((az_bz - u_cz_plus_e).is_zero_vector());
    }
}
