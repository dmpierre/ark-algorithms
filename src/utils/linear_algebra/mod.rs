use ark_ff::PrimeField;
use std::ops::{Add, Mul, Sub};

#[derive(Clone, Debug)]
pub struct Matrix<F: PrimeField> {
    pub rows: Vec<Vector<F>>,
    pub num_rows: usize,
    pub num_cols: usize,
}

#[derive(Clone, Debug)]
pub struct Vector<F: PrimeField> {
    pub elements: Vec<F>,
    pub size: usize,
}

impl<F: PrimeField> Matrix<F> {
    pub fn new(rows: &Vec<Vector<F>>) -> Self {
        Self {
            rows: rows.clone(),
            num_rows: rows.len(),
            num_cols: rows[0].elements.len(),
        }
    }

    pub fn new_from_vecs(rows: &Vec<Vec<F>>) -> Self {
        let mut vec_rows = vec![];
        for row in rows {
            vec_rows.push(Vector::new(row));
        }
        Self::new(&vec_rows)
    }
}

impl<F: PrimeField> Vector<F> {
    pub fn new(elements: &Vec<F>) -> Self {
        Self {
            elements: elements.clone(),
            size: elements.len(),
        }
    }

    pub fn new_zero_vector(size: usize) -> Self {
        Self {
            elements: vec![F::zero(); size],
            size,
        }
    }
}

impl<F: PrimeField> Matrix<F> {
    pub fn dot(&self, rhs: &Matrix<F>) -> Matrix<F> {
        assert_eq!(self.num_cols, rhs.num_rows);
        let mut res = vec![];
        for i in 0..self.num_rows {
            let mut row = vec![];
            for j in 0..rhs.num_cols {
                let mut sum = F::zero();
                for k in 0..self.num_cols {
                    sum += self.rows[i].elements[k] * rhs.rows[k].elements[j];
                }
                row.push(sum);
            }
            res.push(Vector::new(&row));
        }
        Matrix::new(&res)
    }

    pub fn dot_vector(&self, rhs: &Vector<F>) -> Vector<F> {
        assert_eq!(self.num_cols, rhs.size);
        let mut res = vec![F::zero(); self.num_rows];
        for i in 0..self.num_rows {
            for j in 0..self.num_cols {
                res[i] += self.rows[i].elements[j] * rhs.elements[j];
            }
        }
        Vector::new(&res)
    }
}

impl<F: PrimeField> Sub for Vector<F> {
    type Output = Vector<F>;

    fn sub(self, rhs: Self) -> Self::Output {
        assert_eq!(self.elements.len(), rhs.elements.len());
        let mut res = vec![F::zero(); self.elements.len()];
        for i in 0..self.elements.len() {
            res[i] = self.elements[i] - rhs.elements[i];
        }
        Vector::new(&res)
    }
}

impl<F: PrimeField> Mul for Vector<F> {
    type Output = Vector<F>;

    fn mul(self, rhs: Self) -> Self::Output {
        assert_eq!(self.elements.len(), rhs.elements.len());
        let mut res = vec![F::zero(); self.elements.len()];
        for i in 0..self.elements.len() {
            res[i] = self.elements[i] * rhs.elements[i];
        }
        Vector::new(&res)
    }
}

impl<F: PrimeField> Add for Vector<F> {
    type Output = Vector<F>;

    fn add(self, rhs: Self) -> Self::Output {
        assert_eq!(self.elements.len(), rhs.elements.len());
        let mut res = vec![F::zero(); self.elements.len()];
        for i in 0..self.elements.len() {
            res[i] = self.elements[i] + rhs.elements[i];
        }
        Vector::new(&res)
    }
}

impl<F: PrimeField> Vector<F> {
    pub fn is_zero_vector(&self) -> bool {
        for i in 0..self.elements.len() {
            if self.elements[i] != F::zero() {
                return false;
            }
        }
        true
    }

    pub fn scalar_mul(&self, scalar: &F) -> Vector<F> {
        let mut res = vec![F::zero(); self.elements.len()];
        for i in 0..self.elements.len() {
            res[i] = self.elements[i] * scalar;
        }
        Vector::new(&res)
    }
}
