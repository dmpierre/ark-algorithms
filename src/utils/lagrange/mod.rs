use ark_ff::PrimeField;
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, Evaluations, GeneralEvaluationDomain,
};

/// Computes the lagrange interpolation for the set of points:
/// (\omega^{0}, y_0), (\omega^{1}, y_1), ..., (\omega^{n}, y_n)
/// where \omega is a primitive n-th root of unity.
pub fn compute_lagrange_interpolation<F: PrimeField>(evals: &Vec<F>) -> DensePolynomial<F> {
    let k = evals.len() as usize;
    let omegas = GeneralEvaluationDomain::<F>::new(k).unwrap();
    let lagrange: DensePolynomial<F> =
        Evaluations::<F>::from_vec_and_domain(evals.to_vec(), omegas).interpolate();
    lagrange
}