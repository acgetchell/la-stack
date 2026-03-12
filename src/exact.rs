//! Exact determinant sign via arbitrary-precision rational arithmetic.
//!
//! This module is only compiled when the `"exact"` Cargo feature is enabled.
//!
//! # Architecture
//!
//! `det_sign_exact` uses a two-stage adaptive-precision approach inspired by
//! Shewchuk's robust geometric predicates:
//!
//! 1. **Fast filter (D ≤ 4)**: compute `det_direct()` and a conservative error
//!    bound. If `|det| > bound`, the f64 sign is provably correct — return
//!    immediately without allocating.
//! 2. **Exact fallback**: convert entries to `BigRational` and run the Bareiss
//!    algorithm (fraction-free Gaussian elimination) for a guaranteed-correct
//!    sign.

use num_bigint::{BigInt, Sign};
use num_rational::BigRational;
use num_traits::ToPrimitive;

use crate::LaError;
use crate::matrix::Matrix;

/// Validate that all entries in a `D×D` matrix are finite (not NaN or infinite).
///
/// Returns `Ok(())` if all entries are finite, or `Err(LaError::NonFinite)` with
/// the column of the first non-finite entry found.
fn validate_finite<const D: usize>(m: &Matrix<D>) -> Result<(), LaError> {
    for r in 0..D {
        for c in 0..D {
            if !m.rows[r][c].is_finite() {
                return Err(LaError::NonFinite { col: c });
            }
        }
    }
    Ok(())
}

/// Convert an `f64` to an exact `BigRational`.
///
/// Every finite `f64` is exactly representable as a rational number (`m × 2^e`),
/// so this conversion is lossless.
///
/// # Panics
/// Panics if `x` is NaN or infinite.
fn f64_to_bigrational(x: f64) -> BigRational {
    BigRational::from_float(x).expect("non-finite matrix entry in det_sign_exact")
}

/// Compute the exact determinant of a `D×D` matrix using the Bareiss algorithm
/// (fraction-free Gaussian elimination) in `BigRational` arithmetic.
///
/// Returns the determinant as an exact `BigRational`.
fn bareiss_det<const D: usize>(m: &Matrix<D>) -> BigRational {
    if D == 0 {
        return BigRational::from_integer(BigInt::from(1));
    }
    if D == 1 {
        return f64_to_bigrational(m.rows[0][0]);
    }

    // Convert f64 entries to exact BigRational.
    let mut a: Vec<Vec<BigRational>> = Vec::with_capacity(D);
    for r in 0..D {
        let mut row = Vec::with_capacity(D);
        for c in 0..D {
            row.push(f64_to_bigrational(m.rows[r][c]));
        }
        a.push(row);
    }

    let zero = BigRational::from_integer(BigInt::from(0));
    let mut prev_pivot = BigRational::from_integer(BigInt::from(1));
    let mut sign: i8 = 1;

    for k in 0..D {
        // Partial pivoting: find non-zero entry in column k at or below row k.
        if a[k][k] == zero {
            let mut found = false;
            for i in (k + 1)..D {
                if a[i][k] != zero {
                    a.swap(k, i);
                    sign = -sign;
                    found = true;
                    break;
                }
            }
            if !found {
                // Entire column below (and including) diagonal is zero → singular.
                return zero;
            }
        }

        // Bareiss elimination for rows below k.
        for i in (k + 1)..D {
            for j in (k + 1)..D {
                // a[i][j] = (a[k][k] * a[i][j] - a[i][k] * a[k][j]) / prev_pivot
                a[i][j] = (&a[k][k] * &a[i][j] - &a[i][k] * &a[k][j]) / &prev_pivot;
            }
            // Zero out the eliminated column entry (not needed for det, but keeps
            // the matrix consistent for potential debugging).
            a[i][k] = zero.clone();
        }

        prev_pivot = a[k][k].clone();
    }

    let det = &a[D - 1][D - 1];
    if sign < 0 { -det } else { det.clone() }
}

impl<const D: usize> Matrix<D> {
    /// Exact determinant using arbitrary-precision rational arithmetic.
    ///
    /// Returns the determinant as an exact [`BigRational`] value. Every finite
    /// `f64` is exactly representable as a rational, so the conversion is
    /// lossless and the result is provably correct.
    ///
    /// # When to use
    ///
    /// Use this when you need the exact determinant *value* — for example,
    /// exact volume computation or distinguishing truly-degenerate simplices
    /// from near-degenerate ones.  If you only need the *sign*, prefer
    /// [`det_sign_exact`](Self::det_sign_exact) which has a fast f64 filter.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// let m = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
    /// let det = m.det_exact().unwrap();
    /// // det = 1*4 - 2*3 = -2  (exact)
    /// assert_eq!(det, BigRational::from_integer((-2).into()));
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] if any matrix entry is NaN or infinite.
    #[inline]
    pub fn det_exact(&self) -> Result<BigRational, LaError> {
        validate_finite(self)?;
        Ok(bareiss_det(self))
    }

    /// Exact determinant converted to `f64`.
    ///
    /// Computes the exact [`BigRational`] determinant via [`det_exact`](Self::det_exact)
    /// and converts it to the nearest `f64`.  This is useful when you want the
    /// most accurate f64 determinant possible without committing to `BigRational`
    /// in your downstream code.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// let m = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
    /// let det = m.det_exact_f64().unwrap();
    /// assert!((det - (-2.0)).abs() <= f64::EPSILON);
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] if any matrix entry is NaN or infinite,
    /// or if the exact determinant is too large to represent as a finite `f64`.
    #[inline]
    pub fn det_exact_f64(&self) -> Result<f64, LaError> {
        let exact = self.det_exact()?;
        let val = exact.to_f64().unwrap_or(f64::INFINITY);
        if val.is_finite() {
            Ok(val)
        } else {
            Err(LaError::NonFinite { col: 0 })
        }
    }

    /// Exact determinant sign using adaptive-precision arithmetic.
    ///
    /// Returns `1` if `det > 0`, `-1` if `det < 0`, and `0` if `det == 0` (singular).
    ///
    /// For D ≤ 4, a fast f64 filter is tried first: `det_direct()` is compared
    /// against a conservative error bound derived from the matrix permanent.
    /// If the f64 result clearly exceeds the bound, the sign is returned
    /// immediately without allocating.  Otherwise (and always for D ≥ 5), the
    /// Bareiss algorithm runs in exact [`BigRational`] arithmetic.
    ///
    /// # When to use
    ///
    /// Use this when the sign of the determinant must be correct regardless of
    /// floating-point conditioning (e.g. geometric predicates on near-degenerate
    /// configurations).  For well-conditioned matrices the fast filter resolves
    /// the sign without touching `BigRational`, so the overhead is minimal.
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
    /// // This matrix is singular (row 3 = row 1 + row 2 in exact arithmetic).
    /// assert_eq!(m.det_sign_exact().unwrap(), 0);
    ///
    /// assert_eq!(Matrix::<3>::identity().det_sign_exact().unwrap(), 1);
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] if any matrix entry is NaN or infinite.
    #[inline]
    pub fn det_sign_exact(&self) -> Result<i8, LaError> {
        validate_finite(self)?;

        // Stage 1: f64 fast filter for D ≤ 4.
        if let (Some(det_f64), Some(err)) = (self.det_direct(), self.det_errbound()) {
            // When entries are large (e.g. near f64::MAX) the determinant can
            // overflow to infinity even though every individual entry is finite.
            // In that case the fast filter is inconclusive; fall through to the
            // exact Bareiss path.
            if det_f64.is_finite() {
                if det_f64 > err {
                    return Ok(1);
                }
                if det_f64 < -err {
                    return Ok(-1);
                }
            }
        }

        // Stage 2: exact Bareiss fallback.
        let det = bareiss_det(self);
        Ok(match det.numer().sign() {
            Sign::Plus => 1,
            Sign::Minus => -1,
            Sign::NoSign => 0,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::DEFAULT_PIVOT_TOL;

    #[test]
    fn det_sign_exact_d0_is_positive() {
        assert_eq!(Matrix::<0>::zero().det_sign_exact().unwrap(), 1);
    }

    #[test]
    fn det_sign_exact_d1_positive() {
        let m = Matrix::<1>::from_rows([[42.0]]);
        assert_eq!(m.det_sign_exact().unwrap(), 1);
    }

    #[test]
    fn det_sign_exact_d1_negative() {
        let m = Matrix::<1>::from_rows([[-3.5]]);
        assert_eq!(m.det_sign_exact().unwrap(), -1);
    }

    #[test]
    fn det_sign_exact_d1_zero() {
        let m = Matrix::<1>::from_rows([[0.0]]);
        assert_eq!(m.det_sign_exact().unwrap(), 0);
    }

    #[test]
    fn det_sign_exact_identity_2d() {
        assert_eq!(Matrix::<2>::identity().det_sign_exact().unwrap(), 1);
    }

    #[test]
    fn det_sign_exact_identity_3d() {
        assert_eq!(Matrix::<3>::identity().det_sign_exact().unwrap(), 1);
    }

    #[test]
    fn det_sign_exact_identity_4d() {
        assert_eq!(Matrix::<4>::identity().det_sign_exact().unwrap(), 1);
    }

    #[test]
    fn det_sign_exact_identity_5d() {
        assert_eq!(Matrix::<5>::identity().det_sign_exact().unwrap(), 1);
    }

    #[test]
    fn det_sign_exact_singular_duplicate_rows() {
        let m = Matrix::<3>::from_rows([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [1.0, 2.0, 3.0], // duplicate of row 0
        ]);
        assert_eq!(m.det_sign_exact().unwrap(), 0);
    }

    #[test]
    fn det_sign_exact_singular_linear_combination() {
        // Row 2 = row 0 + row 1 in exact arithmetic.
        let m = Matrix::<3>::from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [5.0, 7.0, 9.0]]);
        assert_eq!(m.det_sign_exact().unwrap(), 0);
    }

    #[test]
    fn det_sign_exact_negative_det_row_swap() {
        // Swapping two rows of the identity negates the determinant.
        let m = Matrix::<3>::from_rows([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
        assert_eq!(m.det_sign_exact().unwrap(), -1);
    }

    #[test]
    fn det_sign_exact_negative_det_known() {
        // det([[1,2],[3,4]]) = 1*4 - 2*3 = -2
        let m = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        assert_eq!(m.det_sign_exact().unwrap(), -1);
    }

    #[test]
    fn det_sign_exact_agrees_with_det_for_spd() {
        // SPD matrix → positive determinant.
        let m = Matrix::<3>::from_rows([[4.0, 2.0, 0.0], [2.0, 5.0, 1.0], [0.0, 1.0, 3.0]]);
        assert_eq!(m.det_sign_exact().unwrap(), 1);
        assert!(m.det(DEFAULT_PIVOT_TOL).unwrap() > 0.0);
    }

    /// Near-singular matrix with an exact perturbation.
    ///
    /// The base matrix `[[1,2,3],[4,5,6],[7,8,9]]` is exactly singular (rows in
    /// arithmetic progression).  Adding `2^-50` to entry (0,0) makes
    /// `det = 2^-50 × cofactor(0,0) = 2^-50 × (5×9 − 6×8) = −3 × 2^-50 < 0`.
    /// Both f64 `det_direct()` and `det_sign_exact()` should agree here.
    #[test]
    fn det_sign_exact_near_singular_perturbation() {
        let perturbation = f64::from_bits(0x3CD0_0000_0000_0000); // 2^-50
        let m = Matrix::<3>::from_rows([
            [1.0 + perturbation, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
        ]);
        // Exact: det = perturbation × (5×9 − 6×8) = perturbation × (−3) < 0.
        assert_eq!(m.det_sign_exact().unwrap(), -1);
    }

    /// For D ≤ 4, well-conditioned matrices should hit the fast filter
    /// and never allocate `BigRational`.  We can't directly observe this,
    /// but we verify correctness for a range of known signs.
    #[test]
    fn det_sign_exact_fast_filter_positive_4x4() {
        let m = Matrix::<4>::from_rows([
            [2.0, 1.0, 0.0, 0.0],
            [1.0, 3.0, 1.0, 0.0],
            [0.0, 1.0, 4.0, 1.0],
            [0.0, 0.0, 1.0, 5.0],
        ]);
        // SPD tridiagonal → positive det.
        assert_eq!(m.det_sign_exact().unwrap(), 1);
    }

    #[test]
    fn det_sign_exact_fast_filter_negative_4x4() {
        // Swap rows 0 and 1 of the above → negate det.
        let m = Matrix::<4>::from_rows([
            [1.0, 3.0, 1.0, 0.0],
            [2.0, 1.0, 0.0, 0.0],
            [0.0, 1.0, 4.0, 1.0],
            [0.0, 0.0, 1.0, 5.0],
        ]);
        assert_eq!(m.det_sign_exact().unwrap(), -1);
    }

    #[test]
    fn det_sign_exact_subnormal_entries() {
        // Subnormal f64 values should convert losslessly.
        let tiny = 5e-324_f64; // smallest positive subnormal
        assert!(tiny.is_subnormal());

        let m = Matrix::<2>::from_rows([[tiny, 0.0], [0.0, tiny]]);
        // det = tiny^2 > 0
        assert_eq!(m.det_sign_exact().unwrap(), 1);
    }

    #[test]
    fn det_sign_exact_returns_err_on_nan() {
        let m = Matrix::<2>::from_rows([[f64::NAN, 0.0], [0.0, 1.0]]);
        assert_eq!(m.det_sign_exact(), Err(LaError::NonFinite { col: 0 }));
    }

    #[test]
    fn det_sign_exact_returns_err_on_infinity() {
        let m = Matrix::<2>::from_rows([[f64::INFINITY, 0.0], [0.0, 1.0]]);
        assert_eq!(m.det_sign_exact(), Err(LaError::NonFinite { col: 0 }));
    }

    #[test]
    fn det_sign_exact_returns_err_on_nan_5x5() {
        // D ≥ 5 bypasses the fast filter, exercising the bareiss_det path.
        let mut m = Matrix::<5>::identity();
        m.set(2, 3, f64::NAN);
        assert_eq!(m.det_sign_exact(), Err(LaError::NonFinite { col: 3 }));
    }

    #[test]
    fn det_sign_exact_returns_err_on_infinity_5x5() {
        let mut m = Matrix::<5>::identity();
        m.set(0, 0, f64::INFINITY);
        assert_eq!(m.det_sign_exact(), Err(LaError::NonFinite { col: 0 }));
    }

    #[test]
    fn det_sign_exact_pivot_needed_5x5() {
        // D ≥ 5 skips the fast filter → exercises Bareiss pivoting.
        // Permutation matrix with a single swap (rows 0↔1) → det = −1.
        let m = Matrix::<5>::from_rows([
            [0.0, 1.0, 0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0],
        ]);
        assert_eq!(m.det_sign_exact().unwrap(), -1);
    }

    #[test]
    fn det_sign_exact_5x5_known() {
        // det of a permutation matrix with two swaps = +1 (even permutation).
        let m = Matrix::<5>::from_rows([
            [0.0, 1.0, 0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0],
        ]);
        // Two transpositions → even permutation → det = +1
        assert_eq!(m.det_sign_exact().unwrap(), 1);
    }

    // -----------------------------------------------------------------------
    // Direct tests for internal helpers (coverage of private functions)
    // -----------------------------------------------------------------------

    #[test]
    fn det_errbound_d0_is_zero() {
        assert_eq!(Matrix::<0>::zero().det_errbound(), Some(0.0));
    }

    #[test]
    fn det_errbound_d1_is_zero() {
        assert_eq!(Matrix::<1>::from_rows([[42.0]]).det_errbound(), Some(0.0));
    }

    #[test]
    fn det_errbound_d2_positive() {
        let m = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let bound = m.det_errbound().unwrap();
        assert!(bound > 0.0);
        // bound = ERR_COEFF_2 * (|1*4| + |2*3|) = ERR_COEFF_2 * 10
        assert!(crate::ERR_COEFF_2.mul_add(-10.0, bound).abs() < 1e-30);
    }

    #[test]
    fn det_errbound_d3_positive() {
        let m = Matrix::<3>::identity();
        let bound = m.det_errbound().unwrap();
        assert!(bound > 0.0);
    }

    #[test]
    fn det_errbound_d3_non_identity() {
        // Non-identity matrix to exercise all code paths in D=3 case
        let m = Matrix::<3>::from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 10.0]]);
        let bound = m.det_errbound().unwrap();
        assert!(bound > 0.0);
    }

    #[test]
    fn det_errbound_d4_positive() {
        let m = Matrix::<4>::identity();
        let bound = m.det_errbound().unwrap();
        assert!(bound > 0.0);
    }

    #[test]
    fn det_errbound_d4_non_identity() {
        // Non-identity matrix to exercise all code paths in D=4 case
        let m = Matrix::<4>::from_rows([
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 2.0, 0.0, 0.0],
            [0.0, 0.0, 3.0, 0.0],
            [0.0, 0.0, 0.0, 4.0],
        ]);
        let bound = m.det_errbound().unwrap();
        assert!(bound > 0.0);
    }

    #[test]
    fn det_errbound_d5_is_none() {
        assert_eq!(Matrix::<5>::identity().det_errbound(), None);
    }

    #[test]
    fn bareiss_det_d0_is_one() {
        let det = bareiss_det(&Matrix::<0>::zero());
        assert_eq!(det, BigRational::from_integer(BigInt::from(1)));
    }

    #[test]
    fn bareiss_det_d1_returns_entry() {
        let det = bareiss_det(&Matrix::<1>::from_rows([[7.0]]));
        assert_eq!(det, f64_to_bigrational(7.0));
    }

    #[test]
    fn bareiss_det_d3_with_pivoting() {
        // First column has zero on diagonal → exercises pivot swap + break.
        let m = Matrix::<3>::from_rows([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
        let det = bareiss_det(&m);
        // det of this permutation matrix = -1
        assert_eq!(det, BigRational::from_integer(BigInt::from(-1)));
    }

    #[test]
    fn bareiss_det_singular_all_zeros_in_column() {
        // Column 1 is all zeros below diagonal after elimination → singular.
        let m = Matrix::<3>::from_rows([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
        let det = bareiss_det(&m);
        assert_eq!(det, BigRational::from_integer(BigInt::from(0)));
    }

    #[test]
    fn det_sign_exact_overflow_determinant_finite_entries() {
        // Entries near f64::MAX are finite, but the f64 determinant overflows
        // to infinity.  The fast filter should be skipped and Bareiss should
        // compute the correct positive sign.
        let big = f64::MAX / 2.0;
        assert!(big.is_finite());
        let m = Matrix::<3>::from_rows([[0.0, 0.0, 1.0], [big, 0.0, 1.0], [0.0, big, 1.0]]);
        // det = big^2 > 0
        assert_eq!(m.det_sign_exact().unwrap(), 1);
    }

    // -----------------------------------------------------------------------
    // det_exact tests
    // -----------------------------------------------------------------------

    #[test]
    fn det_exact_d0_is_one() {
        let det = Matrix::<0>::zero().det_exact().unwrap();
        assert_eq!(det, BigRational::from_integer(BigInt::from(1)));
    }

    #[test]
    fn det_exact_identity_3d() {
        let det = Matrix::<3>::identity().det_exact().unwrap();
        assert_eq!(det, BigRational::from_integer(BigInt::from(1)));
    }

    #[test]
    fn det_exact_known_2x2() {
        // det([[1,2],[3,4]]) = 1*4 - 2*3 = -2
        let m = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let det = m.det_exact().unwrap();
        assert_eq!(det, BigRational::from_integer(BigInt::from(-2)));
    }

    #[test]
    fn det_exact_singular_returns_zero() {
        // Rows in arithmetic progression → exactly singular.
        let m = Matrix::<3>::from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]);
        let det = m.det_exact().unwrap();
        assert_eq!(det, BigRational::from_integer(BigInt::from(0)));
    }

    #[test]
    fn det_exact_near_singular_perturbation() {
        // Same 2^-50 perturbation case: exact det = -3 × 2^-50.
        let perturbation = f64::from_bits(0x3CD0_0000_0000_0000); // 2^-50
        let m = Matrix::<3>::from_rows([
            [1.0 + perturbation, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
        ]);
        let det = m.det_exact().unwrap();
        // det should be exactly -3 × 2^-50.
        let expected = BigRational::new(BigInt::from(-3), BigInt::from(1u64 << 50));
        assert_eq!(det, expected);
    }

    #[test]
    fn det_exact_5x5_permutation() {
        // Single swap (rows 0↔1) → det = -1.
        let m = Matrix::<5>::from_rows([
            [0.0, 1.0, 0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0],
        ]);
        let det = m.det_exact().unwrap();
        assert_eq!(det, BigRational::from_integer(BigInt::from(-1)));
    }

    #[test]
    fn det_exact_err_on_nan() {
        let m = Matrix::<2>::from_rows([[f64::NAN, 0.0], [0.0, 1.0]]);
        assert_eq!(m.det_exact(), Err(LaError::NonFinite { col: 0 }));
    }

    #[test]
    fn det_exact_err_on_inf() {
        let m = Matrix::<2>::from_rows([[f64::INFINITY, 0.0], [0.0, 1.0]]);
        assert_eq!(m.det_exact(), Err(LaError::NonFinite { col: 0 }));
    }

    // -----------------------------------------------------------------------
    // det_exact_f64 tests
    // -----------------------------------------------------------------------

    #[test]
    fn det_exact_f64_identity_3d() {
        let det = Matrix::<3>::identity().det_exact_f64().unwrap();
        assert!((det - 1.0).abs() <= f64::EPSILON);
    }

    #[test]
    fn det_exact_f64_known_2x2() {
        let m = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let det = m.det_exact_f64().unwrap();
        assert!((det - (-2.0)).abs() <= f64::EPSILON);
    }

    #[test]
    fn det_exact_f64_agrees_with_det_direct_2d() {
        let m = Matrix::<2>::from_rows([[4.0, 2.0], [2.0, 3.0]]);
        let exact = m.det_exact_f64().unwrap();
        let direct = m.det_direct().unwrap();
        assert!((exact - direct).abs() <= f64::EPSILON);
    }

    #[test]
    fn det_exact_f64_agrees_with_det_direct_3d() {
        let m = Matrix::<3>::from_rows([[2.0, 1.0, 0.0], [0.0, 3.0, 1.0], [1.0, 0.0, 2.0]]);
        let exact = m.det_exact_f64().unwrap();
        let direct = m.det_direct().unwrap();
        assert!((exact - direct).abs() <= 1e-12);
    }

    #[test]
    fn det_exact_f64_agrees_with_det_direct_4d() {
        let m = Matrix::<4>::from_rows([
            [2.0, 1.0, 0.0, 0.0],
            [1.0, 3.0, 1.0, 0.0],
            [0.0, 1.0, 4.0, 1.0],
            [0.0, 0.0, 1.0, 5.0],
        ]);
        let exact = m.det_exact_f64().unwrap();
        let direct = m.det_direct().unwrap();
        assert!((exact - direct).abs() <= 1e-12);
    }

    #[test]
    fn det_exact_f64_err_on_nan() {
        let m = Matrix::<2>::from_rows([[f64::NAN, 0.0], [0.0, 1.0]]);
        assert_eq!(m.det_exact_f64(), Err(LaError::NonFinite { col: 0 }));
    }

    #[test]
    fn det_exact_f64_overflow_returns_err() {
        // Entries near f64::MAX produce a determinant too large for f64.
        let big = f64::MAX / 2.0;
        let m = Matrix::<3>::from_rows([[0.0, 0.0, 1.0], [big, 0.0, 1.0], [0.0, big, 1.0]]);
        // det = big^2, which overflows f64.
        assert_eq!(m.det_exact_f64(), Err(LaError::NonFinite { col: 0 }));
    }

    // -----------------------------------------------------------------------
    // validate_finite tests
    // -----------------------------------------------------------------------

    #[test]
    fn validate_finite_ok_for_finite() {
        assert!(validate_finite(&Matrix::<3>::identity()).is_ok());
    }

    #[test]
    fn validate_finite_err_on_nan() {
        let mut m = Matrix::<2>::identity();
        m.set(1, 0, f64::NAN);
        assert_eq!(validate_finite(&m), Err(LaError::NonFinite { col: 0 }));
    }

    #[test]
    fn validate_finite_err_on_inf() {
        let mut m = Matrix::<2>::identity();
        m.set(0, 1, f64::NEG_INFINITY);
        assert_eq!(validate_finite(&m), Err(LaError::NonFinite { col: 1 }));
    }
}
