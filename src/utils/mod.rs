use ark_ff::PrimeField;
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial};

pub mod linear_algebra;
pub mod lagrange;

pub fn build_zero_polynomial<F: PrimeField>(roots: &Vec<F>) -> DensePolynomial<F> {
    // roots are the values at which the polynomial will be zero
    // (X - roots[0]) * (X - roots[1]) * ... * (X - roots[n])
    let mut polys = vec![];
    for root in roots {
        let poly = DensePolynomial::from_coefficients_vec(vec![*root * (-F::ONE), F::ONE]);
        polys.push(poly);
    }
    // multiply all the different polys together to get one single polynomial
    let mut zero_poly = polys[0].clone();
    for i in 1..polys.len() {
        zero_poly = &zero_poly * &polys[i];
    }
    zero_poly
}
