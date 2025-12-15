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
    ///
    /// # Examples
    /// ```
    /// #![allow(unused_imports)]
    /// use la_stack::prelude::*;
    ///
    /// let m = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
    /// assert_eq!(m.get(0, 1), Some(2.0));
    /// ```
    #[inline]
    pub const fn from_rows(rows: [[f64; D]; D]) -> Self {
        Self { rows }
    }

    /// All-zeros matrix.
    ///
    /// # Examples
    /// ```
    /// #![allow(unused_imports)]
    /// use la_stack::prelude::*;
    ///
    /// let z = Matrix::<2>::zero();
    /// assert_eq!(z.get(1, 1), Some(0.0));
    /// ```
    #[inline]
    pub const fn zero() -> Self {
        Self {
            rows: [[0.0; D]; D],
        }
    }

    /// Identity matrix.
    ///
    /// # Examples
    /// ```
    /// #![allow(unused_imports)]
    /// use la_stack::prelude::*;
    ///
    /// let i = Matrix::<3>::identity();
    /// assert_eq!(i.get(0, 0), Some(1.0));
    /// assert_eq!(i.get(0, 1), Some(0.0));
    /// assert_eq!(i.get(2, 2), Some(1.0));
    /// ```
    #[inline]
    pub fn identity() -> Self {
        let mut m = Self::zero();

        let mut i = 0;
        while i < D {
            m.rows[i][i] = 1.0;
            i += 1;
        }

        m
    }

    /// Get an element with bounds checking.
    ///
    /// # Examples
    /// ```
    /// #![allow(unused_imports)]
    /// use la_stack::prelude::*;
    ///
    /// let m = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
    /// assert_eq!(m.get(1, 0), Some(3.0));
    /// assert_eq!(m.get(2, 0), None);
    /// ```
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
    ///
    /// # Examples
    /// ```
    /// #![allow(unused_imports)]
    /// use la_stack::prelude::*;
    ///
    /// let mut m = Matrix::<2>::zero();
    /// assert!(m.set(0, 1, 2.5));
    /// assert_eq!(m.get(0, 1), Some(2.5));
    /// assert!(!m.set(10, 0, 1.0));
    /// ```
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
    ///
    /// # Examples
    /// ```
    /// #![allow(unused_imports)]
    /// use la_stack::prelude::*;
    ///
    /// let m = Matrix::<2>::from_rows([[1.0, -2.0], [3.0, 4.0]]);
    /// assert!((m.inf_norm() - 7.0).abs() <= 1e-12);
    /// ```
    #[inline]
    #[must_use]
    pub fn inf_norm(&self) -> f64 {
        let mut max_row_sum: f64 = 0.0;

        for row in &self.rows {
            let row_sum: f64 = row.iter().map(|&x| x.abs()).sum();
            if row_sum > max_row_sum {
                max_row_sum = row_sum;
            }
        }

        max_row_sum
    }

    /// Compute an LU decomposition with partial pivoting.
    ///
    /// # Examples
    /// ```
    /// #![allow(unused_imports)]
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let a = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
    /// let lu = a.lu(DEFAULT_PIVOT_TOL)?;
    ///
    /// let b = Vector::<2>::new([5.0, 11.0]);
    /// let x = lu.solve_vec(b)?.into_array();
    ///
    /// assert!((x[0] - 1.0).abs() <= 1e-12);
    /// assert!((x[1] - 2.0).abs() <= 1e-12);
    /// # Ok(())
    /// # }
    /// ```
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
    /// # Examples
    /// ```
    /// #![allow(unused_imports)]
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let det = Matrix::<3>::identity().det(DEFAULT_PIVOT_TOL)?;
    /// assert!((det - 1.0).abs() <= 1e-12);
    /// # Ok(())
    /// # }
    /// ```
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

    #[test]
    fn identity_has_ones_on_diag_and_zeros_off_diag() {
        let m = Matrix::<3>::identity();

        for r in 0..3 {
            for c in 0..3 {
                let expected = if r == c { 1.0 } else { 0.0 };
                assert_approx(m.get(r, c).unwrap(), expected, 0.0);
            }
        }
    }

    #[test]
    fn default_is_zero() {
        let m = Matrix::<3>::default();
        assert_approx(m.inf_norm(), 0.0, 0.0);
    }
}
