use ark_ec::pairing::Pairing;
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial};
use ark_ff::Field;

pub fn build_zero_polynomial<E: Pairing>(roots: &Vec<E::ScalarField>) -> DensePolynomial<<E as Pairing>::ScalarField> {
    // roots are the values at which the polynomial will be zero
    // (X - roots[0]) * (X - roots[1]) * ... * (X - roots[n])
    let mut polys = vec![];
    for root in roots {
        let poly = DensePolynomial::from_coefficients_vec(vec![*root * (-E::ScalarField::ONE), E::ScalarField::ONE]);
        polys.push(poly);
    }
    // multiply all the different polys together to get one single polynomial
    let mut zero_poly = polys[0].clone();
    for i in 1..polys.len() {
        zero_poly = &zero_poly * &polys[i];
    }
    zero_poly
}