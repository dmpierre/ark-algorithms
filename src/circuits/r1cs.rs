use ark_ff::PrimeField;

use crate::utils::linear_algebra::{Matrix, Vector};

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
