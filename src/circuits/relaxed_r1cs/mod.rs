use std::ops::Add;

use ark_ff::PrimeField;

use crate::utils::linear_algebra::{Matrix, Vector};

use super::r1cs::R1CS;

pub type R1CSRelaxedInstanceWitness<F> = Vector<F>;
pub type R1CSRelaxedErrorTerm<F> = Vector<F>;

/// A relaxed R1CS equation
pub struct R1CSRelaxed<F: PrimeField> {
    pub n_constraints: usize,
    pub n_witness: usize,
    pub n_instance: usize,
    pub a: Matrix<F>,
    pub b: Matrix<F>,
    pub c: Matrix<F>,
    pub e: R1CSRelaxedErrorTerm<F>,
    pub u: F,
}

/// An instance to a relaxed R1CS equation
pub struct R1CSRelaxedInstance<F: PrimeField> {
    pub e: R1CSRelaxedErrorTerm<F>,
    pub u: F,
    pub x: Vector<F>,
}

impl<F: PrimeField> From<R1CS<F>> for R1CSRelaxed<F> {
    fn from(value: R1CS<F>) -> Self {
        Self {
            n_constraints: value.n_constraints,
            n_witness: value.n_witness,
            n_instance: value.n_instance,
            a: value.a,
            b: value.b,
            c: value.c,
            e: Vector::new_zero_vector(value.n_constraints),
            u: F::ONE,
        }
    }
}

impl<F: PrimeField> R1CSRelaxed<F> {
    /// Creates a relaxed r1cs by providing all necessary r1cs components with error term and u
    pub fn from_relaxed_r1cs(
        a: Matrix<F>,
        b: Matrix<F>,
        c: Matrix<F>,
        u: F,
        e: R1CSRelaxedErrorTerm<F>,
    ) -> Self {
        Self {
            n_constraints: a.num_rows,
            n_witness: a.num_cols,
            n_instance: b.num_cols,
            a,
            b,
            c,
            e,
            u,
        }
    }

    /// Checks if the relaxed r1cs is satisfied
    pub fn is_satisfied(&self, z: &R1CSRelaxedInstanceWitness<F>) -> bool {
        let az = self.a.dot_vector(z);
        let bz = self.b.dot_vector(z);
        let cz = self.c.dot_vector(z);
        ((az * bz) - cz.scalar_mul(&self.u).add(self.e.clone())).is_zero_vector()
    }

    /// Computes the T term, where:
    /// T = AZ_1 o BZ_2 + AZ_2 o BZ_1 - u_1CZ_2 - u_2CZ_1
    /// It is required for computing the error vector
    /// T is the cross term that pops up when taking linear combinations with naive r1cs
    pub fn compute_t(
        &self,
        rhs: &R1CSRelaxed<F>,
        z1: &R1CSRelaxedInstanceWitness<F>,
        z2: &R1CSRelaxedInstanceWitness<F>,
    ) -> Vector<F> {
        let (u1, u2) = (self.u, rhs.u);
        (self.a.dot_vector(z1)) * (self.b.dot_vector(z2))
            + (self.a.dot_vector(z2)) * (self.b.dot_vector(z1))
            - (self.c.dot_vector(z2)).scalar_mul(&u1)
            - (self.c.dot_vector(z1)).scalar_mul(&u2)
    }

    /// Computes the u term, where:
    /// u = u_1 + r * u_2
    /// This linear combination will be the "new" u for our updated relaxed r1cs
    pub fn compute_u(&self, rhs: &R1CSRelaxed<F>, r: &F) -> F {
        self.u + rhs.u * r
    }

    /// Computes the E term, where:
    /// E = E_1 + r * T + r^2 * E_2
    /// This linear combination will be the "new" E for our updated relaxed r1cs
    pub fn compute_e(
        &self,
        rhs: &R1CSRelaxed<F>,
        r: &F,
        z1: &R1CSRelaxedInstanceWitness<F>,
        z2: &R1CSRelaxedInstanceWitness<F>,
    ) -> Vector<F> {
        let e1 = &self.e;
        let e2 = &rhs.e;
        let r_square = r.square();
        let t = &self.compute_t(rhs, z1, z2);
        (t.scalar_mul(r) + e1.clone()) + e2.scalar_mul(&r_square)
    }

    /// Computes the Z term, where:
    /// Z = Z_1 + r * Z_2
    /// This linear combination will be the "new" Z for our updated relaxed r1cs
    pub fn compute_z(
        &self,
        r: &F,
        z1: &R1CSRelaxedInstanceWitness<F>,
        z2: &R1CSRelaxedInstanceWitness<F>,
    ) -> R1CSRelaxedInstanceWitness<F> {
        z1.clone() + z2.scalar_mul(r)
    }
}

#[cfg(test)]
mod tests {
    use ark_bn254::Fr;
    use ark_std::test_rng;

    use crate::circuits::{
        r1cs::{
            utils::{get_r1cs_from_cs, get_z_from_cs, TestPythagoreCircuit},
            R1CS,
        },
        relaxed_r1cs::R1CSRelaxedInstanceWitness,
    };

    use super::R1CSRelaxed;
    use ark_ff::UniformRand;

    #[test]
    pub fn test_valid_rlc_of_two_relaxed_r1cs() {
        let circuit = TestPythagoreCircuit::new(Fr::from(2), Fr::from(3), Fr::from(13));
        let r1cs: R1CS<Fr> = get_r1cs_from_cs(circuit.clone()).unwrap();
        let relaxed_r1cs_1 = R1CSRelaxed::from(r1cs.clone());
        let z_1: R1CSRelaxedInstanceWitness<Fr> = get_z_from_cs(circuit.clone()).unwrap();
        assert!(relaxed_r1cs_1.is_satisfied(&z_1));

        let circuit = TestPythagoreCircuit::new(Fr::from(5), Fr::from(10), Fr::from(125));
        let relaxed_r1cs_2 = R1CSRelaxed::from(r1cs.clone()); // same relaxed r1cs as before
        let z_2: R1CSRelaxedInstanceWitness<Fr> = get_z_from_cs(circuit.clone()).unwrap();
        assert!(relaxed_r1cs_2.is_satisfied(&z_2));

        let mut rng = test_rng();
        let r = Fr::rand(&mut rng);

        // Compute the new relaxed r1cs
        let e_3 = relaxed_r1cs_1.compute_e(&relaxed_r1cs_2, &r, &z_1, &z_2);
        let u_3 = relaxed_r1cs_1.compute_u(&relaxed_r1cs_2, &r);
        let relaxed_r1cs_3 = R1CSRelaxed::from_relaxed_r1cs(r1cs.a, r1cs.b, r1cs.c, u_3, e_3);

        // Compute satisfying instance-witness for the new relaxed r1cs
        let z_3 = relaxed_r1cs_1.compute_z(&r, &z_1, &z_2);
        assert!(relaxed_r1cs_3.is_satisfied(&z_3));
    }
}
