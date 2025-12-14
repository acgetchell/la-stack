//! LU decomposition and solves.

use crate::LaError;
use crate::matrix::Matrix;
use crate::vector::Vector;

/// LU decomposition (PA = LU) with partial pivoting.
#[must_use]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Lu<const D: usize> {
    factors: Matrix<D>,
    piv: [usize; D],
    piv_sign: f64,
    tol: f64,
}

impl<const D: usize> Lu<D> {
    #[inline]
    pub(crate) fn factor(a: Matrix<D>, tol: f64) -> Result<Self, LaError> {
        let mut lu = a;

        let mut piv = [0usize; D];
        for (i, p) in piv.iter_mut().enumerate() {
            *p = i;
        }

        let mut piv_sign = 1.0;

        for k in 0..D {
            // Choose pivot row.
            let mut pivot_row = k;
            let mut pivot_abs = lu.rows[k][k].abs();
            if !pivot_abs.is_finite() {
                return Err(LaError::NonFinite { pivot_col: k });
            }

            for r in (k + 1)..D {
                let v = lu.rows[r][k].abs();
                if !v.is_finite() {
                    return Err(LaError::NonFinite { pivot_col: k });
                }
                if v > pivot_abs {
                    pivot_abs = v;
                    pivot_row = r;
                }
            }

            if pivot_abs <= tol {
                return Err(LaError::Singular { pivot_col: k });
            }

            if pivot_row != k {
                lu.rows.swap(k, pivot_row);
                piv.swap(k, pivot_row);
                piv_sign = -piv_sign;
            }

            let pivot = lu.rows[k][k];
            if !pivot.is_finite() {
                return Err(LaError::NonFinite { pivot_col: k });
            }

            // Eliminate below pivot.
            for r in (k + 1)..D {
                let mult = lu.rows[r][k] / pivot;
                if !mult.is_finite() {
                    return Err(LaError::NonFinite { pivot_col: k });
                }
                lu.rows[r][k] = mult;

                for c in (k + 1)..D {
                    lu.rows[r][c] = (-mult).mul_add(lu.rows[k][c], lu.rows[r][c]);
                }
            }
        }

        Ok(Self {
            factors: lu,
            piv,
            piv_sign,
            tol,
        })
    }

    /// Solve `A x = b` using this LU factorization.
    ///
    /// # Errors
    /// Returns [`LaError::Singular`] if a diagonal of `U` is (numerically) zero.
    /// Returns [`LaError::NonFinite`] if NaN/∞ is detected.
    #[inline]
    pub fn solve_vec(&self, b: Vector<D>) -> Result<Vector<D>, LaError> {
        let mut x = [0.0; D];
        for (i, x_i) in x.iter_mut().enumerate() {
            *x_i = b.data[self.piv[i]];
        }

        // Forward substitution for L (unit diagonal).
        for i in 0..D {
            let mut sum = x[i];
            let row = self.factors.rows[i];
            for (j, x_j) in x.iter().enumerate().take(i) {
                sum = (-row[j]).mul_add(*x_j, sum);
            }
            if !sum.is_finite() {
                return Err(LaError::NonFinite { pivot_col: i });
            }
            x[i] = sum;
        }

        // Back substitution for U.
        for ii in 0..D {
            let i = D - 1 - ii;
            let mut sum = x[i];
            let row = self.factors.rows[i];
            for (j, x_j) in x.iter().enumerate().skip(i + 1) {
                sum = (-row[j]).mul_add(*x_j, sum);
            }

            let diag = row[i];
            if !diag.is_finite() || !sum.is_finite() {
                return Err(LaError::NonFinite { pivot_col: i });
            }
            if diag.abs() <= self.tol {
                return Err(LaError::Singular { pivot_col: i });
            }

            x[i] = sum / diag;
        }

        Ok(Vector::new(x))
    }

    /// Determinant of the original matrix.
    #[inline]
    #[must_use]
    pub fn det(&self) -> f64 {
        let mut det = self.piv_sign;
        for i in 0..D {
            det *= self.factors.rows[i][i];
        }
        det
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::DEFAULT_PIVOT_TOL;

    fn assert_approx(a: f64, b: f64, eps: f64) {
        assert!((a - b).abs() <= eps, "{a} !~= {b} (eps={eps})");
    }

    #[test]
    fn solve_2x2_basic() {
        let a = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let lu = Lu::factor(a, DEFAULT_PIVOT_TOL).unwrap();
        let b = Vector::<2>::new([5.0, 11.0]);
        let x = lu.solve_vec(b).unwrap().into_array();
        assert_approx(x[0], 1.0, 1e-12);
        assert_approx(x[1], 2.0, 1e-12);
    }

    #[test]
    fn solve_requires_pivoting() {
        let a = Matrix::<2>::from_rows([[0.0, 1.0], [1.0, 0.0]]);
        let lu = a.lu(DEFAULT_PIVOT_TOL).unwrap();
        let b = Vector::<2>::new([1.0, 2.0]);
        let x = lu.solve_vec(b).unwrap().into_array();
        assert_approx(x[0], 2.0, 1e-12);
        assert_approx(x[1], 1.0, 1e-12);
    }

    #[test]
    fn singular_detected() {
        let a = Matrix::<2>::from_rows([[1.0, 2.0], [2.0, 4.0]]);
        let err = a.lu(DEFAULT_PIVOT_TOL).unwrap_err();
        assert_eq!(err, LaError::Singular { pivot_col: 1 });
    }
}
