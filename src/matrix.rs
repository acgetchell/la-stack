//! Fixed-size, stack-allocated square matrices.

use crate::LaError;
use crate::lu::Lu;

/// Fixed-size square matrix `D×D`, stored inline.
#[must_use]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Matrix<const D: usize> {
    pub(crate) rows: [[f64; D]; D],
}

impl<const D: usize> Matrix<D> {
    /// Construct from row-major storage.
    #[inline]
    pub const fn from_rows(rows: [[f64; D]; D]) -> Self {
        Self { rows }
    }

    /// All-zeros matrix.
    #[inline]
    pub const fn zero() -> Self {
        Self {
            rows: [[0.0; D]; D],
        }
    }

    /// Identity matrix.
    #[inline]
    pub fn identity() -> Self {
        let mut m = Self::zero();
        for i in 0..D {
            m.rows[i][i] = 1.0;
        }
        m
    }

    /// Get an element with bounds checking.
    #[inline]
    #[must_use]
    pub fn get(&self, r: usize, c: usize) -> Option<f64> {
        if r < D && c < D {
            Some(self.rows[r][c])
        } else {
            None
        }
    }

    /// Set an element with bounds checking.
    ///
    /// Returns `true` if the index was in-bounds.
    #[inline]
    pub fn set(&mut self, r: usize, c: usize, value: f64) -> bool {
        if r < D && c < D {
            self.rows[r][c] = value;
            true
        } else {
            false
        }
    }

    /// Infinity norm (maximum absolute row sum).
    #[inline]
    #[must_use]
    pub fn inf_norm(&self) -> f64 {
        let mut max_row_sum = 0.0;
        for r in 0..D {
            let mut row_sum = 0.0;
            for c in 0..D {
                row_sum += self.rows[r][c].abs();
            }
            if row_sum > max_row_sum {
                max_row_sum = row_sum;
            }
        }
        max_row_sum
    }

    /// Compute an LU decomposition with partial pivoting.
    ///
    /// # Errors
    /// Returns [`LaError::Singular`] if no suitable pivot (|pivot| > `tol`) exists for a column.
    /// Returns [`LaError::NonFinite`] if NaN/∞ is detected during factorization.
    #[inline]
    pub fn lu(self, tol: f64) -> Result<Lu<D>, LaError> {
        Lu::factor(self, tol)
    }

    /// Determinant computed via LU decomposition.
    ///
    /// # Errors
    /// Propagates LU factorization errors (e.g. singular matrices).
    #[inline]
    pub fn det(self, tol: f64) -> Result<f64, LaError> {
        self.lu(tol).map(|lu| lu.det())
    }
}

impl<const D: usize> Default for Matrix<D> {
    #[inline]
    fn default() -> Self {
        Self::zero()
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
    fn get_set_bounds_checked() {
        let mut m = Matrix::<2>::zero();
        assert!(m.set(0, 0, 1.0));
        assert_eq!(m.get(0, 0), Some(1.0));

        assert!(!m.set(2, 0, 1.0));
        assert_eq!(m.get(2, 0), None);
    }

    #[test]
    fn inf_norm_max_row_sum() {
        let m = Matrix::<2>::from_rows([[1.0, -2.0], [3.0, 4.0]]);
        assert_approx(m.inf_norm(), 7.0, 0.0);
    }

    #[test]
    fn det_identity_is_one() {
        let det = Matrix::<3>::identity().det(DEFAULT_PIVOT_TOL).unwrap();
        assert_approx(det, 1.0, 1e-12);
    }
}
