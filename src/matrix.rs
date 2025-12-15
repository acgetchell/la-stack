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
    pub const fn identity() -> Self {
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
    pub const fn get(&self, r: usize, c: usize) -> Option<f64> {
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
    pub const fn set(&mut self, r: usize, c: usize, value: f64) -> bool {
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

    use approx::assert_abs_diff_eq;
    use pastey::paste;

    macro_rules! gen_public_api_matrix_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<public_api_matrix_from_rows_get_set_bounds_checked_ $d d>]() {
                    let mut rows = [[0.0f64; $d]; $d];
                    rows[0][0] = 1.0;
                    rows[$d - 1][$d - 1] = -2.0;

                    let mut m = Matrix::<$d>::from_rows(rows);

                    assert_eq!(m.get(0, 0), Some(1.0));
                    assert_eq!(m.get($d - 1, $d - 1), Some(-2.0));

                    // Out-of-bounds is None.
                    assert_eq!(m.get($d, 0), None);

                    // Out-of-bounds set fails.
                    assert!(!m.set($d, 0, 3.0));

                    // In-bounds set works.
                    assert!(m.set(0, $d - 1, 3.0));
                    assert_eq!(m.get(0, $d - 1), Some(3.0));
                }

                #[test]
                fn [<public_api_matrix_zero_and_default_are_zero_ $d d>]() {
                    let z = Matrix::<$d>::zero();
                    assert_abs_diff_eq!(z.inf_norm(), 0.0, epsilon = 0.0);

                    let d = Matrix::<$d>::default();
                    assert_abs_diff_eq!(d.inf_norm(), 0.0, epsilon = 0.0);
                }

                #[test]
                fn [<public_api_matrix_inf_norm_max_row_sum_ $d d>]() {
                    let mut rows = [[0.0f64; $d]; $d];

                    // Row 0 has absolute row sum = D.
                    for c in 0..$d {
                        rows[0][c] = -1.0;
                    }

                    // Row 1 has smaller absolute row sum.
                    for c in 0..$d {
                        rows[1][c] = 0.5;
                    }

                    let m = Matrix::<$d>::from_rows(rows);
                    assert_abs_diff_eq!(m.inf_norm(), f64::from($d), epsilon = 0.0);
                }

                #[test]
                fn [<public_api_matrix_identity_lu_det_solve_vec_ $d d>]() {
                    let m = Matrix::<$d>::identity();

                    // Identity has ones on diag and zeros off diag.
                    for r in 0..$d {
                        for c in 0..$d {
                            let expected = if r == c { 1.0 } else { 0.0 };
                            assert_abs_diff_eq!(m.get(r, c).unwrap(), expected, epsilon = 0.0);
                        }
                    }

                    // Determinant is 1.
                    let det = m.det(DEFAULT_PIVOT_TOL).unwrap();
                    assert_abs_diff_eq!(det, 1.0, epsilon = 1e-12);

                    // LU solve on identity returns the RHS.
                    let lu = m.lu(DEFAULT_PIVOT_TOL).unwrap();

                    let b_arr = {
                        let mut arr = [0.0f64; $d];
                        let values = [1.0f64, 2.0, 3.0, 4.0, 5.0];
                        for (dst, src) in arr.iter_mut().zip(values.iter()) {
                            *dst = *src;
                        }
                        arr
                    };

                    let b = crate::Vector::<$d>::new(b_arr);
                    let x = lu.solve_vec(b).unwrap().into_array();

                    for (x_i, b_i) in x.iter().zip(b_arr.iter()) {
                        assert_abs_diff_eq!(*x_i, *b_i, epsilon = 1e-12);
                    }
                }
            }
        };
    }

    // Mirror delaunay-style multi-dimension tests.
    gen_public_api_matrix_tests!(2);
    gen_public_api_matrix_tests!(3);
    gen_public_api_matrix_tests!(4);
    gen_public_api_matrix_tests!(5);
}
