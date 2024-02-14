pub mod utils;
use crate::utils::linear_algebra::{Matrix, Vector};
/// A lot of code has been forked from https://github.com/privacy-scaling-explorations/folding-schemes
/// It includes things such as how r1cs matrices or the z vector are extracted
/// It has been adapted here and there, in minor ways.
/// Thanks Arnau! :)
use ark_ff::PrimeField;

use self::utils::R1CSInstanceWitness;

/// A "regular" R1CS equation
#[derive(Clone, Debug)]
pub struct R1CS<F: PrimeField> {
    pub n_constraints: usize,
    pub n_witness: usize,
    pub n_instance: usize,
    pub a: Matrix<F>,
    pub b: Matrix<F>,
    pub c: Matrix<F>,
}

impl<F: PrimeField> R1CS<F> {
    pub fn is_satisfied(&self, z: &R1CSInstanceWitness<F>) -> bool {
        let az = self.a.dot_vector(z);
        let bz = self.b.dot_vector(z);
        let cz = self.c.dot_vector(z);
        ((az * bz) - cz).is_zero_vector()
    }
}

#[cfg(test)]
mod test {
    use ark_pallas::Fr;
    use ark_relations::r1cs::{ConstraintSynthesizer, ConstraintSystem};

    use crate::{
        circuits::r1cs::{
            utils::{
                get_r1cs_from_cs, get_test_r1cs, get_test_satisfying_witness, get_z_from_cs,
                TestPythagoreCircuit,
            },
            R1CS,
        },
        utils::linear_algebra::{Matrix, Vector},
    };

    #[test]
    pub fn test_raw_r1cs_is_satisfied() {
        let (a, b, c): (Matrix<Fr>, Matrix<Fr>, Matrix<Fr>) = get_test_r1cs();
        let witness: Vector<Fr> = get_test_satisfying_witness(5);
        let a_dot_w = a.dot_vector(&witness);
        let b_dot_w = b.dot_vector(&witness);
        let c_dot_w = c.dot_vector(&witness);
        let a_times_b = a_dot_w * b_dot_w;
        let a_times_b_minus_c = a_times_b - c_dot_w;
        assert!(a_times_b_minus_c.is_zero_vector());
    }

    #[test]
    pub fn test_pythagore_circuit() {
        let cs = ConstraintSystem::<Fr>::new_ref();
        let circuit = TestPythagoreCircuit::new(Fr::from(5), Fr::from(10), Fr::from(125));
        circuit.generate_constraints(cs.clone()).unwrap();
        assert!(cs.is_satisfied().unwrap());

        let circuit = TestPythagoreCircuit::new(Fr::from(5), Fr::from(10), Fr::from(10));
        circuit.generate_constraints(cs.clone()).unwrap();
        assert!(!cs.is_satisfied().unwrap());
    }

    #[test]
    pub fn test_valid_r1cs_is_satisfying() {
        let circuit = TestPythagoreCircuit::new(Fr::from(5), Fr::from(10), Fr::from(125));
        let r1cs: R1CS<Fr> = get_r1cs_from_cs(circuit.clone()).unwrap();
        let z = get_z_from_cs(circuit.clone()).unwrap();
        assert!(r1cs.is_satisfied(&z));
    }

    #[test]
    pub fn test_invalid_r1cs_is_not_satisfying() {
        let circuit = TestPythagoreCircuit::new(Fr::from(1), Fr::from(1), Fr::from(100));
        let r1cs: R1CS<Fr> = get_r1cs_from_cs(circuit.clone()).unwrap();
        let z = get_z_from_cs(circuit.clone()).unwrap();
        assert!(!r1cs.is_satisfied(&z));
    }
}
