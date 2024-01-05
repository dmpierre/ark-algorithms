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

#[cfg(test)]
mod test {
    use ark_pallas::Fr;

    use crate::{
        circuits::r1cs::{get_test_r1cs, get_test_satisfying_witness},
        utils::linear_algebra::{Matrix, Vector},
    };

    #[test]
    pub fn test_r1cs_is_satisfied() {
        let (a, b, c): (Matrix<Fr>, Matrix<Fr>, Matrix<Fr>) = get_test_r1cs();
        let witness: Vector<Fr> = get_test_satisfying_witness(5);
        let a_dot_w = a.dot_vector(&witness);
        let b_dot_w = b.dot_vector(&witness);
        let c_dot_w = c.dot_vector(&witness);
        let a_times_b = a_dot_w * b_dot_w;
        let a_times_b_minus_c = a_times_b - c_dot_w;
        assert!(a_times_b_minus_c.is_zero_vector());
    }
}
