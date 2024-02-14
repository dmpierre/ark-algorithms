use ark_ff::PrimeField;
use ark_r1cs_std::{alloc::AllocVar, eq::EqGadget, fields::fp::FpVar};
use ark_relations::r1cs::{ConstraintSynthesizer, ConstraintSystem, ConstraintSystemRef};

use crate::utils::linear_algebra::{Matrix, Vector};

use super::R1CS;

/// Type for an instance-witness tuple
pub type R1CSInstanceWitness<F> = Vector<F>;

/// a^2 + b^2 = c^2 where a, b are witnesses, and c is public io
#[derive(Clone, Debug)]
pub struct TestPythagoreCircuit<F: PrimeField> {
    a: F,
    b: F,
    c_square: F,
}

impl<F: PrimeField> TestPythagoreCircuit<F> {
    pub fn new(a: F, b: F, c_square: F) -> Self {
        Self { a, b, c_square }
    }
}

pub fn generate_constraint_system<F: PrimeField>(
    circuit: impl ConstraintSynthesizer<F>,
) -> Result<ConstraintSystem<F>, String> {
    let cs = ConstraintSystem::<F>::new_ref();
    circuit
        .generate_constraints(cs.clone())
        .map_err(|_| "Error generating constraints")?;
    cs.finalize();
    let cs = cs.into_inner().ok_or("Error returning inner cs")?;
    Ok(cs)
}

pub fn get_r1cs_from_cs<F: PrimeField>(
    circuit: impl ConstraintSynthesizer<F>,
) -> Result<R1CS<F>, String> {
    let cs = generate_constraint_system(circuit)?;
    let r1cs = extract_r1cs::<F>(&cs);
    Ok(r1cs)
}

/// Returns the instance-witness vector in (W, 1, x) format
pub fn get_z_from_cs<F: PrimeField>(
    circuit: impl ConstraintSynthesizer<F>,
) -> Result<R1CSInstanceWitness<F>, String> {
    let cs = generate_constraint_system(circuit)?;
    let z = extract_z::<F>(&cs);
    Ok(z)
}

/// forked from https://github.com/privacy-scaling-explorations/folding-schemes
/// adapted here and there to fit my types
/// thanks Arnau! :)
pub fn extract_r1cs<F: PrimeField>(cs: &ConstraintSystem<F>) -> R1CS<F> {
    let m = cs.to_matrices().unwrap();
    let n_constraints = cs.num_constraints;
    let n_witness = cs.num_witness_variables;
    let n_instance = cs.num_instance_variables;
    let n_rows = cs.num_constraints;
    let n_cols = cs.num_instance_variables + cs.num_witness_variables; // cs.num_instance_variables already counts the 1
    let a = Matrix::new_from_ark_matrix(&m.a, n_rows, n_cols);
    let b = Matrix::new_from_ark_matrix(&m.b, n_rows, n_cols);
    let c = Matrix::new_from_ark_matrix(&m.c, n_rows, n_cols);
    R1CS {
        n_constraints,
        n_witness,
        n_instance,
        a,
        b,
        c,
    }
}

pub fn extract_z<F: PrimeField>(cs: &ConstraintSystem<F>) -> R1CSInstanceWitness<F> {
    let mut z = cs.instance_assignment.clone(); // starts with pub io
    let mut witness = cs.witness_assignment.clone();
    _ = z.append(&mut witness);
    Vector::new(&z)
}

impl<F: PrimeField> ConstraintSynthesizer<F> for TestPythagoreCircuit<F> {
    fn generate_constraints(
        self,
        cs: ConstraintSystemRef<F>,
    ) -> Result<(), ark_relations::r1cs::SynthesisError> {
        let a = FpVar::new_witness(cs.clone(), || Ok(self.a)).unwrap();
        let b = FpVar::new_witness(cs.clone(), || Ok(self.b)).unwrap();
        let c_square = FpVar::new_input(cs.clone(), || Ok(self.c_square)).unwrap();
        let a_square = a.clone() * a.clone();
        let b_square = b.clone() * b.clone();
        let c_square_computed = a_square.clone() + b_square.clone();
        c_square_computed.enforce_equal(&c_square)
    }
}

/// Returns a "raw" r1cs, composed of three matrices A, B and C
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

/// From: https://github.com/privacy-scaling-explorations/folding-schemes/blob/05f49918ac35fba62bd43389943f6f5d33e78cd7/src/ccs/r1cs.rs#L105
pub fn get_test_satisfying_witness<F: PrimeField>(input: usize) -> Vector<F> {
    // z = (1, io, w)
    let input = F::from(input as u64);
    Vector::new(&vec![
        F::ONE,
        input,                                            // io
        input * input * input + input + F::from(5 as u8), // x^3 + x + 5
        input * input,                                    // x^2
        input * input * input,                            // x^2 * x
        input * input * input + input,                    // x^3 + x
    ])
}
