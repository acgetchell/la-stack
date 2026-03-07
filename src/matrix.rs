//! Fixed-size, stack-allocated square matrices.

use crate::LaError;
use crate::ldlt::Ldlt;
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
    /// Returns [`LaError::Singular`] if, for some column `k`, the largest-magnitude candidate pivot
    /// in that column satisfies `|pivot| <= tol` (so no numerically usable pivot exists).
    /// Returns [`LaError::NonFinite`] if NaN/∞ is detected during factorization.
    #[inline]
    pub fn lu(self, tol: f64) -> Result<Lu<D>, LaError> {
        Lu::factor(self, tol)
    }

    /// Compute an LDLT factorization (`A = L D Lᵀ`) without pivoting.
    ///
    /// This is intended for symmetric positive definite (SPD) and positive semi-definite (PSD)
    /// matrices such as Gram matrices.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let a = Matrix::<2>::from_rows([[4.0, 2.0], [2.0, 3.0]]);
    /// let ldlt = a.ldlt(DEFAULT_SINGULAR_TOL)?;
    ///
    /// // det(A) = 8
    /// assert!((ldlt.det() - 8.0).abs() <= 1e-12);
    ///
    /// // Solve A x = b
    /// let b = Vector::<2>::new([1.0, 2.0]);
    /// let x = ldlt.solve_vec(b)?.into_array();
    /// assert!((x[0] - (-0.125)).abs() <= 1e-12);
    /// assert!((x[1] - 0.75).abs() <= 1e-12);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::Singular`] if, for some step `k`, the required diagonal entry `d = D[k,k]`
    /// is `<= tol` (non-positive or too small). This treats PSD degeneracy (and indefinite inputs)
    /// as singular/degenerate.
    /// Returns [`LaError::NonFinite`] if NaN/∞ is detected during factorization.
    #[inline]
    pub fn ldlt(self, tol: f64) -> Result<Ldlt<D>, LaError> {
        Ldlt::factor(self, tol)
    }

    /// Closed-form determinant for dimensions 0–4, bypassing LU factorization.
    ///
    /// Returns `Some(det)` for `D` ∈ {0, 1, 2, 3, 4}, `None` for D ≥ 5.
    /// `D = 0` returns `Some(1.0)` (empty product).
    /// This is a `const fn` (Rust 1.94+) and uses fused multiply-add (`mul_add`)
    /// for improved accuracy and performance.
    ///
    /// For a determinant that works for any dimension (falling back to LU for D ≥ 5),
    /// use [`det`](Self::det).
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// let m = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
    /// assert!((m.det_direct().unwrap() - (-2.0)).abs() <= 1e-12);
    ///
    /// // D = 0 is the empty product.
    /// assert_eq!(Matrix::<0>::zero().det_direct(), Some(1.0));
    ///
    /// // D ≥ 5 returns None.
    /// assert!(Matrix::<5>::identity().det_direct().is_none());
    /// ```
    #[inline]
    #[must_use]
    pub const fn det_direct(&self) -> Option<f64> {
        match D {
            0 => Some(1.0),
            1 => Some(self.rows[0][0]),
            2 => {
                // ad - bc
                Some(self.rows[0][0].mul_add(self.rows[1][1], -(self.rows[0][1] * self.rows[1][0])))
            }
            3 => {
                // Cofactor expansion on first row.
                let m00 =
                    self.rows[1][1].mul_add(self.rows[2][2], -(self.rows[1][2] * self.rows[2][1]));
                let m01 =
                    self.rows[1][0].mul_add(self.rows[2][2], -(self.rows[1][2] * self.rows[2][0]));
                let m02 =
                    self.rows[1][0].mul_add(self.rows[2][1], -(self.rows[1][1] * self.rows[2][0]));
                Some(
                    self.rows[0][0]
                        .mul_add(m00, (-self.rows[0][1]).mul_add(m01, self.rows[0][2] * m02)),
                )
            }
            4 => {
                // Cofactor expansion on first row → four 3×3 sub-determinants.
                // Hoist the 6 unique 2×2 minors from rows 2–3 (each used twice).
                let r = &self.rows;

                // 2×2 minors: s_ij = r[2][i]*r[3][j] - r[2][j]*r[3][i]
                let s23 = r[2][2].mul_add(r[3][3], -(r[2][3] * r[3][2])); // cols 2,3
                let s13 = r[2][1].mul_add(r[3][3], -(r[2][3] * r[3][1])); // cols 1,3
                let s12 = r[2][1].mul_add(r[3][2], -(r[2][2] * r[3][1])); // cols 1,2
                let s03 = r[2][0].mul_add(r[3][3], -(r[2][3] * r[3][0])); // cols 0,3
                let s02 = r[2][0].mul_add(r[3][2], -(r[2][2] * r[3][0])); // cols 0,2
                let s01 = r[2][0].mul_add(r[3][1], -(r[2][1] * r[3][0])); // cols 0,1

                // 3×3 cofactors via row 1 expansion using hoisted minors.
                let c00 = r[1][1].mul_add(s23, (-r[1][2]).mul_add(s13, r[1][3] * s12));
                let c01 = r[1][0].mul_add(s23, (-r[1][2]).mul_add(s03, r[1][3] * s02));
                let c02 = r[1][0].mul_add(s13, (-r[1][1]).mul_add(s03, r[1][3] * s01));
                let c03 = r[1][0].mul_add(s12, (-r[1][1]).mul_add(s02, r[1][2] * s01));

                Some(r[0][0].mul_add(
                    c00,
                    (-r[0][1]).mul_add(c01, r[0][2].mul_add(c02, -(r[0][3] * c03))),
                ))
            }
            _ => None,
        }
    }

    /// Determinant, using closed-form formulas for D ≤ 4 and LU decomposition for D ≥ 5.
    ///
    /// For D ∈ {1, 2, 3, 4}, this bypasses LU factorization entirely for a significant
    /// speedup (see [`det_direct`](Self::det_direct)). The `tol` parameter is only used
    /// by the LU fallback path for D ≥ 5.
    ///
    /// # Examples
    /// ```
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
    /// Returns [`LaError::NonFinite`] if the result contains NaN or infinity.
    /// For D ≥ 5, propagates LU factorization errors (e.g. [`LaError::Singular`]).
    #[inline]
    pub fn det(self, tol: f64) -> Result<f64, LaError> {
        if let Some(d) = self.det_direct() {
            return if d.is_finite() {
                Ok(d)
            } else {
                Err(LaError::NonFinite { pivot_col: 0 })
            };
        }
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
    use std::hint::black_box;

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

    // === det_direct tests ===

    #[test]
    fn det_direct_d0_is_one() {
        assert_eq!(Matrix::<0>::zero().det_direct(), Some(1.0));
    }

    #[test]
    fn det_direct_d1_returns_element() {
        let m = Matrix::<1>::from_rows([[42.0]]);
        assert_eq!(m.det_direct(), Some(42.0));
    }

    #[test]
    fn det_direct_d2_known_value() {
        // [[1,2],[3,4]] → det = 1*4 - 2*3 = -2
        // black_box prevents compile-time constant folding of the const fn.
        let m = black_box(Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]));
        assert_abs_diff_eq!(m.det_direct().unwrap(), -2.0, epsilon = 1e-15);
    }

    #[test]
    fn det_direct_d3_known_value() {
        // Classic 3×3: det = 0
        let m = black_box(Matrix::<3>::from_rows([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
        ]));
        assert_abs_diff_eq!(m.det_direct().unwrap(), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn det_direct_d3_nonsingular() {
        // [[2,1,0],[0,3,1],[1,0,2]] → det = 2*(6-0) - 1*(0-1) + 0 = 13
        let m = black_box(Matrix::<3>::from_rows([
            [2.0, 1.0, 0.0],
            [0.0, 3.0, 1.0],
            [1.0, 0.0, 2.0],
        ]));
        assert_abs_diff_eq!(m.det_direct().unwrap(), 13.0, epsilon = 1e-12);
    }

    #[test]
    fn det_direct_d4_identity() {
        let m = black_box(Matrix::<4>::identity());
        assert_abs_diff_eq!(m.det_direct().unwrap(), 1.0, epsilon = 1e-15);
    }

    #[test]
    fn det_direct_d4_known_value() {
        // Diagonal matrix: det = product of diagonal entries.
        let mut rows = [[0.0f64; 4]; 4];
        rows[0][0] = 2.0;
        rows[1][1] = 3.0;
        rows[2][2] = 5.0;
        rows[3][3] = 7.0;
        let m = black_box(Matrix::<4>::from_rows(rows));
        assert_abs_diff_eq!(m.det_direct().unwrap(), 210.0, epsilon = 1e-12);
    }

    #[test]
    fn det_direct_d5_returns_none() {
        assert_eq!(Matrix::<5>::identity().det_direct(), None);
    }

    #[test]
    fn det_direct_d8_returns_none() {
        assert_eq!(Matrix::<8>::zero().det_direct(), None);
    }

    macro_rules! gen_det_direct_agrees_with_lu {
        ($d:literal) => {
            paste! {
                #[test]
                #[allow(clippy::cast_precision_loss)] // r, c, D are tiny integers
                fn [<det_direct_agrees_with_lu_ $d d>]() {
                    // Well-conditioned matrix: diagonally dominant.
                    let mut rows = [[0.0f64; $d]; $d];
                    for r in 0..$d {
                        for c in 0..$d {
                            rows[r][c] = if r == c {
                                (r as f64) + f64::from($d) + 1.0
                            } else {
                                0.1 / ((r + c + 1) as f64)
                            };
                        }
                    }
                    let m = Matrix::<$d>::from_rows(rows);
                    let direct = m.det_direct().unwrap();
                    let lu_det = m.lu(DEFAULT_PIVOT_TOL).unwrap().det();
                    let eps = lu_det.abs().mul_add(1e-12, 1e-12);
                    assert_abs_diff_eq!(direct, lu_det, epsilon = eps);
                }
            }
        };
    }

    gen_det_direct_agrees_with_lu!(1);
    gen_det_direct_agrees_with_lu!(2);
    gen_det_direct_agrees_with_lu!(3);
    gen_det_direct_agrees_with_lu!(4);

    #[test]
    fn det_direct_identity_all_dims() {
        assert_abs_diff_eq!(
            Matrix::<1>::identity().det_direct().unwrap(),
            1.0,
            epsilon = 0.0
        );
        assert_abs_diff_eq!(
            Matrix::<2>::identity().det_direct().unwrap(),
            1.0,
            epsilon = 0.0
        );
        assert_abs_diff_eq!(
            Matrix::<3>::identity().det_direct().unwrap(),
            1.0,
            epsilon = 0.0
        );
        assert_abs_diff_eq!(
            Matrix::<4>::identity().det_direct().unwrap(),
            1.0,
            epsilon = 0.0
        );
    }

    #[test]
    fn det_direct_zero_matrix() {
        assert_abs_diff_eq!(
            Matrix::<2>::zero().det_direct().unwrap(),
            0.0,
            epsilon = 0.0
        );
        assert_abs_diff_eq!(
            Matrix::<3>::zero().det_direct().unwrap(),
            0.0,
            epsilon = 0.0
        );
        assert_abs_diff_eq!(
            Matrix::<4>::zero().det_direct().unwrap(),
            0.0,
            epsilon = 0.0
        );
    }

    #[test]
    fn det_returns_nonfinite_error_for_nan_d2() {
        let m = Matrix::<2>::from_rows([[f64::NAN, 1.0], [1.0, 1.0]]);
        assert_eq!(
            m.det(DEFAULT_PIVOT_TOL),
            Err(LaError::NonFinite { pivot_col: 0 })
        );
    }

    #[test]
    fn det_returns_nonfinite_error_for_inf_d3() {
        let m =
            Matrix::<3>::from_rows([[f64::INFINITY, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        assert_eq!(
            m.det(DEFAULT_PIVOT_TOL),
            Err(LaError::NonFinite { pivot_col: 0 })
        );
    }

    #[test]
    fn det_direct_is_const_evaluable_d2() {
        // Const evaluation proves the function is truly const fn.
        const DET: Option<f64> = {
            let m = Matrix::<2>::from_rows([[1.0, 0.0], [0.0, 1.0]]);
            m.det_direct()
        };
        assert_eq!(DET, Some(1.0));
    }

    #[test]
    fn det_direct_is_const_evaluable_d3() {
        const DET: Option<f64> = {
            let m = Matrix::<3>::from_rows([[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 5.0]]);
            m.det_direct()
        };
        assert_eq!(DET, Some(30.0));
    }
}
