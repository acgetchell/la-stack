//! Fixed-size, stack-allocated square matrices.

use core::hint::cold_path;

use crate::LaError;
use crate::ldlt::Ldlt;
use crate::lu::Lu;
use crate::{ERR_COEFF_2, ERR_COEFF_3, ERR_COEFF_4};

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
    /// # Non-finite handling
    /// If any entry is NaN, the result is NaN.  NaN is detected explicitly
    /// because a naive `row_sum > max_row_sum` comparison silently skips NaN
    /// rows (every ordered comparison against NaN is `false`).  If any entry
    /// is infinite (and no entry is NaN), the result is `+∞`.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// let m = Matrix::<2>::from_rows([[1.0, -2.0], [3.0, 4.0]]);
    /// assert!((m.inf_norm() - 7.0).abs() <= 1e-12);
    ///
    /// // NaN entries propagate to the norm.
    /// let nan = Matrix::<2>::from_rows([[f64::NAN, 1.0], [2.0, 3.0]]);
    /// assert!(nan.inf_norm().is_nan());
    /// ```
    #[inline]
    #[must_use]
    pub const fn inf_norm(&self) -> f64 {
        let mut max_row_sum: f64 = 0.0;

        let mut r = 0;
        while r < D {
            // Iterator chains like `row.iter().map(|x| x.abs()).sum()` are
            // not yet const-stable, so accumulate the absolute row sum with
            // a manual `while` loop.
            let row = &self.rows[r];
            let mut row_sum: f64 = 0.0;
            let mut c = 0;
            while c < D {
                row_sum += row[c].abs();
                c += 1;
            }
            // Propagate NaN explicitly: `f64::max` drops NaN (IEEE 754 `maxNum`)
            // and `f64::maximum` (IEEE 754-2019 `maximum`) is still unstable,
            // so we short-circuit on NaN instead.
            if row_sum.is_nan() {
                cold_path();
                return f64::NAN;
            }
            if row_sum > max_row_sum {
                max_row_sum = row_sum;
            }
            r += 1;
        }

        max_row_sum
    }

    /// Returns `true` if the matrix is symmetric within a relative tolerance.
    ///
    /// Two entries `self[r][c]` and `self[c][r]` are considered equal (for the
    /// purposes of symmetry) when
    /// `|self[r][c] - self[c][r]| <= rel_tol * max(1.0, self.inf_norm())`.
    /// This mirrors the predicate used internally by the debug-build symmetry
    /// check inside [`ldlt`](Self::ldlt), so callers can pre-validate matrices
    /// that may come from untrusted sources without relying on a debug-only
    /// panic.
    ///
    /// Use [`first_asymmetry`](Self::first_asymmetry) to locate the first
    /// offending pair when this returns `false`.
    ///
    /// # NaN / infinity handling
    /// Any non-finite `|self[r][c] - self[c][r]|` (NaN or ±∞) causes this
    /// predicate to return `false`.  This catches both NaN off-diagonals and
    /// asymmetric pairs where one side is infinite and the other is finite
    /// (which would otherwise slip through when `inf_norm()` blows `eps` up
    /// to `+∞` and makes `diff > eps` trivially false).  A matrix whose
    /// [`inf_norm`](Self::inf_norm) is `+∞` can still tolerate *finite*
    /// asymmetries under an infinite `eps` — callers who need strict equality
    /// on large-magnitude finite entries should validate finiteness
    /// separately.
    ///
    /// # Panics
    /// In debug builds, panics if `rel_tol` is negative or NaN; in release
    /// builds these are silently treated as garbage-in garbage-out, matching
    /// the convention of [`lu`](Self::lu) and [`ldlt`](Self::ldlt).
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// let a = Matrix::<2>::from_rows([[4.0, 2.0], [2.0, 3.0]]);
    /// assert!(a.is_symmetric(1e-12));
    ///
    /// let b = Matrix::<2>::from_rows([[4.0, 2.0], [3.0, 3.0]]);
    /// assert!(!b.is_symmetric(1e-12));
    /// ```
    #[inline]
    #[must_use]
    pub fn is_symmetric(&self, rel_tol: f64) -> bool {
        self.first_asymmetry(rel_tol).is_none()
    }

    /// Returns the indices `(r, c)` (with `r < c`) of the first off-diagonal
    /// pair that violates symmetry, or `None` if the matrix is symmetric
    /// within `rel_tol`.
    ///
    /// Iteration order is row-major over the strict upper triangle, so the
    /// returned indices are the lexicographically smallest such pair.  The
    /// predicate is the same as [`is_symmetric`](Self::is_symmetric):
    /// `|self[r][c] - self[c][r]| <= rel_tol * max(1.0, self.inf_norm())`.
    ///
    /// # Panics
    /// In debug builds, panics if `rel_tol` is negative or NaN.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// let a = Matrix::<3>::from_rows([
    ///     [1.0, 2.0, 0.0],
    ///     [2.0, 4.0, 5.0],
    ///     [0.0, 6.0, 9.0], // 6.0 breaks symmetry with a[1][2] = 5.0
    /// ]);
    /// assert_eq!(a.first_asymmetry(1e-12), Some((1, 2)));
    /// assert_eq!(Matrix::<3>::identity().first_asymmetry(1e-12), None);
    /// ```
    #[inline]
    #[must_use]
    pub fn first_asymmetry(&self, rel_tol: f64) -> Option<(usize, usize)> {
        debug_assert!(
            rel_tol >= 0.0,
            "rel_tol must be non-negative (got {rel_tol})"
        );
        let eps = rel_tol * self.inf_norm().max(1.0);
        for r in 0..D {
            for c in (r + 1)..D {
                let diff = (self.rows[r][c] - self.rows[c][r]).abs();
                // Any non-finite `diff` is reported as asymmetric:
                //  * NaN contaminates one side only, and `diff > eps` would
                //    silently skip it because ordered comparisons against NaN
                //    are always `false`.
                //  * ±∞ arises when exactly one of `self[r][c]` / `self[c][r]`
                //    is infinite; a naive `diff > eps` misses this when the
                //    matrix's `inf_norm()` pushes `eps` to `+∞` (because
                //    `∞ > ∞` is `false`).
                if !diff.is_finite() || diff > eps {
                    cold_path();
                    return Some((r, c));
                }
            }
        }
        None
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
    /// # Preconditions
    /// **The input matrix `self` must be symmetric** — that is, `self[i][j] == self[j][i]`
    /// (within rounding) for all `i`, `j`.  This is a *correctness* precondition, not merely
    /// a performance hint.
    ///
    /// - In **debug builds** a `debug_assert!` verifies symmetry via
    ///   [`is_symmetric`](Self::is_symmetric) (relative tolerance scaled by the matrix's
    ///   infinity norm) and panics if it fails.
    /// - In **release builds** the check is compiled out for performance.  An asymmetric
    ///   input will be accepted silently and produce a mathematically meaningless
    ///   factorization — subsequent calls to [`Ldlt::det`] and [`Ldlt::solve_vec`] will
    ///   return wrong results with no error.
    ///
    /// Callers who cannot statically guarantee symmetry should pre-validate with
    /// [`is_symmetric`](Self::is_symmetric) (or locate the offending pair with
    /// [`first_asymmetry`](Self::first_asymmetry)) before calling `ldlt`.  If you need a
    /// general-purpose factorization that tolerates non-symmetric inputs, use
    /// [`lu`](Self::lu) instead.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// // Note the symmetric layout: a[0][1] == a[1][0] == 2.0.
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
    ///
    /// Note that an *asymmetric* input is **not** reported as an error in release builds —
    /// see the [Preconditions](#preconditions) section above.
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
            _ => {
                // Cold in the common D ≤ 4 case; callers fall back to LU for D ≥ 5.
                cold_path();
                None
            }
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
                cold_path();
                // Scan for the first non-finite entry to preserve coordinates.
                for r in 0..D {
                    for c in 0..D {
                        if !self.rows[r][c].is_finite() {
                            return Err(LaError::non_finite_cell(r, c));
                        }
                    }
                }
                // All entries are finite but the determinant overflowed.
                Err(LaError::non_finite_at(0))
            };
        }
        self.lu(tol).map(|lu| lu.det())
    }

    /// Conservative absolute error bound for `det_direct()`.
    ///
    /// Returns `Some(bound)` such that `|det_direct() - det_exact| ≤ bound`,
    /// or `None` for D ≥ 5 where no fast bound is available.
    ///
    /// For D ≤ 4, the bound is derived from the absolute Leibniz sum using
    /// Shewchuk-style error analysis (see `REFERENCES.md` \[8\] and the
    /// per-constant docs on [`ERR_COEFF_2`](crate::ERR_COEFF_2),
    /// [`ERR_COEFF_3`](crate::ERR_COEFF_3), and
    /// [`ERR_COEFF_4`](crate::ERR_COEFF_4)). For D = 0 or 1, returns
    /// `Some(0.0)` since the determinant computation is exact (no
    /// arithmetic).
    ///
    /// This method does NOT require the `exact` feature — the bounds use
    /// pure f64 arithmetic and are useful for custom adaptive-precision logic.
    ///
    /// # When to use
    ///
    /// Use this to build adaptive-precision logic: if `|det_direct()| > bound`,
    /// the f64 sign is provably correct. Otherwise fall back to exact arithmetic.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// let m = Matrix::<3>::from_rows([
    ///     [1.0, 2.0, 3.0],
    ///     [4.0, 5.0, 6.0],
    ///     [7.0, 8.0, 9.0],
    /// ]);
    /// let bound = m.det_errbound().unwrap();
    /// let det_approx = m.det_direct().unwrap();
    /// // If |det_approx| > bound, the sign is guaranteed correct.
    /// ```
    ///
    /// # Adaptive precision pattern (requires `exact` feature)
    /// ```ignore
    /// use la_stack::prelude::*;
    ///
    /// let m = Matrix::<3>::identity();
    /// if let Some(bound) = m.det_errbound() {
    ///     let det = m.det_direct().unwrap();
    ///     if det.abs() > bound {
    ///         // f64 sign is guaranteed correct
    ///         let sign = det.signum() as i8;
    ///     } else {
    ///         // Fall back to exact arithmetic (requires `exact` feature)
    ///         let sign = m.det_sign_exact().unwrap();
    ///     }
    /// } else {
    ///     // D ≥ 5: no fast filter, use exact directly
    ///     let sign = m.det_sign_exact().unwrap();
    /// }
    /// ```
    #[must_use]
    #[inline]
    pub const fn det_errbound(&self) -> Option<f64> {
        match D {
            0 | 1 => Some(0.0), // No arithmetic — result is exact.
            2 => {
                let r = &self.rows;
                let permanent = (r[0][0] * r[1][1]).abs() + (r[0][1] * r[1][0]).abs();
                Some(ERR_COEFF_2 * permanent)
            }
            3 => {
                let r = &self.rows;
                let pm00 = (r[1][1] * r[2][2]).abs() + (r[1][2] * r[2][1]).abs();
                let pm01 = (r[1][0] * r[2][2]).abs() + (r[1][2] * r[2][0]).abs();
                let pm02 = (r[1][0] * r[2][1]).abs() + (r[1][1] * r[2][0]).abs();
                let permanent = r[0][2]
                    .abs()
                    .mul_add(pm02, r[0][1].abs().mul_add(pm01, r[0][0].abs() * pm00));
                Some(ERR_COEFF_3 * permanent)
            }
            4 => {
                let r = &self.rows;
                // 2×2 minor permanents from rows 2–3.
                let sp23 = (r[2][2] * r[3][3]).abs() + (r[2][3] * r[3][2]).abs();
                let sp13 = (r[2][1] * r[3][3]).abs() + (r[2][3] * r[3][1]).abs();
                let sp12 = (r[2][1] * r[3][2]).abs() + (r[2][2] * r[3][1]).abs();
                let sp03 = (r[2][0] * r[3][3]).abs() + (r[2][3] * r[3][0]).abs();
                let sp02 = (r[2][0] * r[3][2]).abs() + (r[2][2] * r[3][0]).abs();
                let sp01 = (r[2][0] * r[3][1]).abs() + (r[2][1] * r[3][0]).abs();
                // 3×3 cofactor permanents from row 1.
                let pc0 = r[1][3]
                    .abs()
                    .mul_add(sp12, r[1][2].abs().mul_add(sp13, r[1][1].abs() * sp23));
                let pc1 = r[1][3]
                    .abs()
                    .mul_add(sp02, r[1][2].abs().mul_add(sp03, r[1][0].abs() * sp23));
                let pc2 = r[1][3]
                    .abs()
                    .mul_add(sp01, r[1][1].abs().mul_add(sp03, r[1][0].abs() * sp13));
                let pc3 = r[1][2]
                    .abs()
                    .mul_add(sp01, r[1][1].abs().mul_add(sp02, r[1][0].abs() * sp12));
                let permanent = r[0][3].abs().mul_add(
                    pc3,
                    r[0][2]
                        .abs()
                        .mul_add(pc2, r[0][1].abs().mul_add(pc1, r[0][0].abs() * pc0)),
                );
                Some(ERR_COEFF_4 * permanent)
            }
            _ => None,
        }
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
            Err(LaError::NonFinite {
                row: Some(0),
                col: 0
            })
        );
    }

    #[test]
    fn det_returns_nonfinite_error_for_inf_d3() {
        let m =
            Matrix::<3>::from_rows([[f64::INFINITY, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        assert_eq!(
            m.det(DEFAULT_PIVOT_TOL),
            Err(LaError::NonFinite {
                row: Some(0),
                col: 0
            })
        );
    }

    #[test]
    fn det_returns_nonfinite_error_for_overflow_with_finite_entries() {
        // det_direct produces an overflowing f64 (1e300 * 1e300 = ∞) even
        // though every matrix entry is finite.  The entry scan in `det`
        // falls through and returns NonFinite { row: None, col: 0 } to signal
        // a computed overflow rather than a NaN/∞ input.
        let m = Matrix::<2>::from_rows([[1e300, 0.0], [0.0, 1e300]]);
        assert_eq!(
            m.det(DEFAULT_PIVOT_TOL),
            Err(LaError::NonFinite { row: None, col: 0 })
        );
    }

    // === det_direct const-evaluability tests (D = 2..=5) ===
    //
    // Every dimension hits a distinct arm of the `match D { … }` body inside
    // `det_direct`, so exercising each at compile time is the tightest
    // const-fn proof available.

    macro_rules! gen_det_direct_const_eval_tests {
        ($d:literal) => {
            paste! {
                /// `Matrix::<D>::det_direct()` on the identity must const-evaluate
                /// to `Some(1.0)` for every closed-form dimension `D ∈ {1, 2, 3, 4}`.
                #[test]
                fn [<det_direct_const_eval_ $d d>]() {
                    const DET: Option<f64> = Matrix::<$d>::identity().det_direct();
                    assert_eq!(DET, Some(1.0));
                }
            }
        };
    }

    gen_det_direct_const_eval_tests!(2);
    gen_det_direct_const_eval_tests!(3);
    gen_det_direct_const_eval_tests!(4);

    #[test]
    fn det_direct_const_eval_d5_is_none() {
        // D ≥ 5 has no closed-form arm; `det_direct` returns `None`.  Verify
        // that the wildcard arm is reachable in a `const { … }` context.
        const DET: Option<f64> = Matrix::<5>::identity().det_direct();
        assert_eq!(DET, None);
    }

    // === det_errbound tests (no `exact` feature required) ===

    #[test]
    fn det_errbound_available_without_exact_feature() {
        // Verify det_errbound is accessible without exact feature
        let m = Matrix::<3>::identity();
        let bound = m.det_errbound();
        assert!(bound.is_some());
        assert!(bound.unwrap() > 0.0);
    }

    #[test]
    fn det_errbound_d5_returns_none() {
        // D=5 has no fast filter
        assert_eq!(Matrix::<5>::identity().det_errbound(), None);
    }

    // === det_errbound const-evaluability tests (D = 2..=5) ===

    macro_rules! gen_det_errbound_const_eval_tests {
        ($d:literal) => {
            paste! {
                /// `Matrix::<D>::det_errbound()` on the identity must const-evaluate
                /// to `Some(bound)` with `bound > 0` for every closed-form dimension
                /// `D ∈ {2, 3, 4}`.  Each dimension hits a distinct arm of
                /// `det_errbound` with a dimension-specific permanent computation.
                #[test]
                fn [<det_errbound_const_eval_ $d d>]() {
                    const BOUND: Option<f64> = Matrix::<$d>::identity().det_errbound();
                    assert!(BOUND.is_some());
                    assert!(BOUND.unwrap() > 0.0);
                }
            }
        };
    }

    gen_det_errbound_const_eval_tests!(2);
    gen_det_errbound_const_eval_tests!(3);
    gen_det_errbound_const_eval_tests!(4);

    #[test]
    fn det_errbound_const_eval_d5_is_none() {
        // D ≥ 5 has no fast-filter bound; `det_errbound` returns `None`.
        const BOUND: Option<f64> = Matrix::<5>::identity().det_errbound();
        assert_eq!(BOUND, None);
    }

    // === inf_norm const-evaluability tests (D = 2..=5) ===

    macro_rules! gen_inf_norm_const_eval_tests {
        ($d:literal) => {
            paste! {
                /// `Matrix::<D>::inf_norm()` on the identity must const-evaluate
                /// to `1.0` for every `D ≥ 1` — each row has a single `1.0`
                /// entry, so the max absolute row sum is exactly `1.0`.
                #[test]
                fn [<inf_norm_const_eval_ $d d>]() {
                    const NORM: f64 = Matrix::<$d>::identity().inf_norm();
                    assert!((NORM - 1.0).abs() <= 1e-12);
                }
            }
        };
    }

    gen_inf_norm_const_eval_tests!(2);
    gen_inf_norm_const_eval_tests!(3);
    gen_inf_norm_const_eval_tests!(4);
    gen_inf_norm_const_eval_tests!(5);

    // === inf_norm NaN / Inf propagation (regression tests for #85) ===

    macro_rules! gen_inf_norm_nonfinite_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<inf_norm_all_nan_returns_nan_ $d d>]() {
                    // Before the fix, `NaN > max_row_sum` was always false, so a
                    // matrix full of NaN silently produced inf_norm == 0.0.
                    let m = Matrix::<$d>::from_rows([[f64::NAN; $d]; $d]);
                    assert!(m.inf_norm().is_nan());
                }

                #[test]
                fn [<inf_norm_single_nan_entry_returns_nan_ $d d>]() {
                    // A single NaN entry must contaminate its row sum and
                    // propagate through `f64::maximum` to the final result.
                    let mut rows = [[0.0f64; $d]; $d];
                    rows[0][0] = f64::NAN;
                    rows[$d - 1][$d - 1] = 1.0;
                    let m = Matrix::<$d>::from_rows(rows);
                    assert!(m.inf_norm().is_nan());
                }

                #[test]
                fn [<inf_norm_infinity_entry_propagates_ $d d>]() {
                    // Infinity entries should propagate to +∞ via the row sum,
                    // not be silently dropped.  The norm is a sum of absolute
                    // values, so any infinite result is necessarily +∞.
                    let mut rows = [[0.0f64; $d]; $d];
                    rows[0][0] = f64::INFINITY;
                    let m = Matrix::<$d>::from_rows(rows);
                    let norm = m.inf_norm();
                    assert!(norm.is_infinite() && norm.is_sign_positive());
                }
            }
        };
    }

    gen_inf_norm_nonfinite_tests!(2);
    gen_inf_norm_nonfinite_tests!(3);
    gen_inf_norm_nonfinite_tests!(4);
    gen_inf_norm_nonfinite_tests!(5);

    // === is_symmetric / first_asymmetry (public LDLT preconditions helpers) ===

    macro_rules! gen_is_symmetric_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<is_symmetric_true_for_identity_ $d d>]() {
                    let m = Matrix::<$d>::identity();
                    assert!(m.is_symmetric(1e-12));
                    assert_eq!(m.first_asymmetry(1e-12), None);
                }

                #[test]
                fn [<is_symmetric_true_for_zero_ $d d>]() {
                    let m = Matrix::<$d>::zero();
                    assert!(m.is_symmetric(1e-12));
                    assert_eq!(m.first_asymmetry(1e-12), None);
                }

                #[test]
                fn [<is_symmetric_true_for_constructed_symmetric_ $d d>]() {
                    // Construct A = M + Mᵀ so A is provably symmetric.
                    let mut m = [[0.0f64; $d]; $d];
                    for r in 0..$d {
                        for c in 0..$d {
                            #[allow(clippy::cast_precision_loss)]
                            {
                                m[r][c] = (r * $d + c) as f64;
                            }
                        }
                    }
                    let mut sym = [[0.0f64; $d]; $d];
                    for r in 0..$d {
                        for c in 0..$d {
                            sym[r][c] = m[r][c] + m[c][r];
                        }
                    }
                    let a = Matrix::<$d>::from_rows(sym);
                    assert!(a.is_symmetric(1e-12));
                    assert_eq!(a.first_asymmetry(1e-12), None);
                }

                #[test]
                fn [<is_symmetric_false_for_asymmetric_offdiagonal_ $d d>]() {
                    // Perturb a single off-diagonal entry so symmetry fails.
                    let mut rows = [[0.0f64; $d]; $d];
                    for i in 0..$d {
                        rows[i][i] = 1.0;
                    }
                    rows[0][$d - 1] = 1.0;
                    rows[$d - 1][0] = -1.0; // breaks symmetry
                    let a = Matrix::<$d>::from_rows(rows);
                    assert!(!a.is_symmetric(1e-12));
                    assert_eq!(a.first_asymmetry(1e-12), Some((0, $d - 1)));
                }

                #[test]
                fn [<is_symmetric_false_for_nan_offdiagonal_ $d d>]() {
                    // A NaN off-diagonal must be detected as asymmetric.
                    let mut rows = [[0.0f64; $d]; $d];
                    for i in 0..$d {
                        rows[i][i] = 1.0;
                    }
                    rows[0][1] = f64::NAN;
                    rows[1][0] = f64::NAN;
                    let a = Matrix::<$d>::from_rows(rows);
                    assert!(!a.is_symmetric(1e-12));
                    // (0, 1) is the first upper-triangular pair involving the NaN.
                    assert_eq!(a.first_asymmetry(1e-12), Some((0, 1)));
                }
            }
        };
    }

    gen_is_symmetric_tests!(2);
    gen_is_symmetric_tests!(3);
    gen_is_symmetric_tests!(4);
    gen_is_symmetric_tests!(5);

    #[test]
    fn is_symmetric_tolerance_scales_with_inf_norm() {
        // Off-diagonal entries differ by 1e-6.  With inf_norm ≈ 2e6, the
        // relative tolerance 1e-12 yields eps ≈ 2e-6, which accepts the gap;
        // a stricter tol of 1e-15 rejects it.
        let a = Matrix::<2>::from_rows([[1.0e6, 1.0e6 + 1.0e-6], [1.0e6, 1.0e6]]);
        assert!(a.is_symmetric(1e-12));
        assert!(!a.is_symmetric(1e-15));
    }

    #[test]
    fn first_asymmetry_returns_lexicographically_first_pair() {
        // Two asymmetric pairs: (0, 2) and (1, 2).  We must get (0, 2) first.
        let a = Matrix::<3>::from_rows([[1.0, 0.0, 2.0], [0.0, 1.0, 3.0], [-2.0, -3.0, 1.0]]);
        assert_eq!(a.first_asymmetry(1e-12), Some((0, 2)));
    }

    /// Regression: a single infinite off-diagonal paired with a finite entry
    /// used to slip through as "symmetric" because `inf_norm()` blew `eps` up
    /// to `+∞` and `∞ > ∞` evaluates to `false`.  After the fix, any
    /// non-finite `|a[r][c] - a[c][r]|` is reported as an asymmetry regardless
    /// of `eps`.
    #[test]
    fn first_asymmetry_flags_infinite_offdiagonal_against_finite() {
        let a = Matrix::<2>::from_rows([[1.0, f64::INFINITY], [0.0, 1.0]]);
        assert_eq!(a.first_asymmetry(1e-12), Some((0, 1)));
        assert!(!a.is_symmetric(1e-12));
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "rel_tol must be non-negative")]
    fn first_asymmetry_debug_panics_on_negative_tol() {
        // Mirrors the `debug_assert!(tol >= 0.0)` convention used by
        // `Matrix::lu` / `Matrix::ldlt`.
        let _ = Matrix::<2>::identity().first_asymmetry(-1.0);
    }
}
