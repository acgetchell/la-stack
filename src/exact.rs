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

use crate::matrix::Matrix;

// ---------------------------------------------------------------------------
// Error-bound constants for the f64 fast filter.
//
// Each constant bounds the absolute error of `det_direct()` relative to the
// *permanent* (sum of absolute products in the Leibniz expansion).  The
// constants are conservative over-estimates following Shewchuk's methodology;
// over-estimating just means we fall through to Bareiss more often.
// ---------------------------------------------------------------------------

const EPS: f64 = f64::EPSILON; // 2^-52

/// D=2: one f64 multiply + one FMA → 2 rounding events.
const ERR_COEFF_2: f64 = 3.0 * EPS + 16.0 * EPS * EPS;

/// D=3: three 2×2 FMA minors + nested FMA combination.
const ERR_COEFF_3: f64 = 8.0 * EPS + 64.0 * EPS * EPS;

/// D=4: six hoisted 2×2 minors → four 3×3 cofactors → FMA row combination.
const ERR_COEFF_4: f64 = 12.0 * EPS + 128.0 * EPS * EPS;

/// Conservative absolute error bound for `det_direct()`, or `None` for D ≥ 5.
///
/// Returns `Some(bound)` such that `|det_direct() - det_exact| ≤ bound`.
fn det_errbound<const D: usize>(m: &Matrix<D>) -> Option<f64> {
    match D {
        0 | 1 => Some(0.0), // No arithmetic — result is exact.
        2 => {
            let r = &m.rows;
            let permanent = (r[0][0] * r[1][1]).abs() + (r[0][1] * r[1][0]).abs();
            Some(ERR_COEFF_2 * permanent)
        }
        3 => {
            let r = &m.rows;
            let pm00 = (r[1][1] * r[2][2]).abs() + (r[1][2] * r[2][1]).abs();
            let pm01 = (r[1][0] * r[2][2]).abs() + (r[1][2] * r[2][0]).abs();
            let pm02 = (r[1][0] * r[2][1]).abs() + (r[1][1] * r[2][0]).abs();
            let permanent = r[0][2]
                .abs()
                .mul_add(pm02, r[0][1].abs().mul_add(pm01, r[0][0].abs() * pm00));
            Some(ERR_COEFF_3 * permanent)
        }
        4 => {
            let r = &m.rows;
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
    /// assert_eq!(m.det_sign_exact(), 0);
    ///
    /// assert_eq!(Matrix::<3>::identity().det_sign_exact(), 1);
    /// ```
    ///
    /// # Panics
    /// Panics if any matrix entry is NaN or infinite.
    #[must_use]
    pub fn det_sign_exact(&self) -> i8 {
        // Stage 1: f64 fast filter for D ≤ 4.
        if let (Some(det_f64), Some(err)) = (self.det_direct(), det_errbound(self)) {
            if det_f64 > err {
                return 1;
            }
            if det_f64 < -err {
                return -1;
            }
        }

        // Stage 2: exact Bareiss fallback.
        let det = bareiss_det(self);
        match det.numer().sign() {
            Sign::Plus => 1,
            Sign::Minus => -1,
            Sign::NoSign => 0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::DEFAULT_PIVOT_TOL;

    #[test]
    fn det_sign_exact_d0_is_positive() {
        assert_eq!(Matrix::<0>::zero().det_sign_exact(), 1);
    }

    #[test]
    fn det_sign_exact_d1_positive() {
        let m = Matrix::<1>::from_rows([[42.0]]);
        assert_eq!(m.det_sign_exact(), 1);
    }

    #[test]
    fn det_sign_exact_d1_negative() {
        let m = Matrix::<1>::from_rows([[-3.5]]);
        assert_eq!(m.det_sign_exact(), -1);
    }

    #[test]
    fn det_sign_exact_d1_zero() {
        let m = Matrix::<1>::from_rows([[0.0]]);
        assert_eq!(m.det_sign_exact(), 0);
    }

    #[test]
    fn det_sign_exact_identity_2d() {
        assert_eq!(Matrix::<2>::identity().det_sign_exact(), 1);
    }

    #[test]
    fn det_sign_exact_identity_3d() {
        assert_eq!(Matrix::<3>::identity().det_sign_exact(), 1);
    }

    #[test]
    fn det_sign_exact_identity_4d() {
        assert_eq!(Matrix::<4>::identity().det_sign_exact(), 1);
    }

    #[test]
    fn det_sign_exact_identity_5d() {
        assert_eq!(Matrix::<5>::identity().det_sign_exact(), 1);
    }

    #[test]
    fn det_sign_exact_singular_duplicate_rows() {
        let m = Matrix::<3>::from_rows([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [1.0, 2.0, 3.0], // duplicate of row 0
        ]);
        assert_eq!(m.det_sign_exact(), 0);
    }

    #[test]
    fn det_sign_exact_singular_linear_combination() {
        // Row 2 = row 0 + row 1 in exact arithmetic.
        let m = Matrix::<3>::from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [5.0, 7.0, 9.0]]);
        assert_eq!(m.det_sign_exact(), 0);
    }

    #[test]
    fn det_sign_exact_negative_det_row_swap() {
        // Swapping two rows of the identity negates the determinant.
        let m = Matrix::<3>::from_rows([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
        assert_eq!(m.det_sign_exact(), -1);
    }

    #[test]
    fn det_sign_exact_negative_det_known() {
        // det([[1,2],[3,4]]) = 1*4 - 2*3 = -2
        let m = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        assert_eq!(m.det_sign_exact(), -1);
    }

    #[test]
    fn det_sign_exact_agrees_with_det_for_spd() {
        // SPD matrix → positive determinant.
        let m = Matrix::<3>::from_rows([[4.0, 2.0, 0.0], [2.0, 5.0, 1.0], [0.0, 1.0, 3.0]]);
        assert_eq!(m.det_sign_exact(), 1);
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
        assert_eq!(m.det_sign_exact(), -1);
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
        assert_eq!(m.det_sign_exact(), 1);
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
        assert_eq!(m.det_sign_exact(), -1);
    }

    #[test]
    fn det_sign_exact_subnormal_entries() {
        // Subnormal f64 values should convert losslessly.
        let tiny = 5e-324_f64; // smallest positive subnormal
        assert!(tiny.is_subnormal());

        let m = Matrix::<2>::from_rows([[tiny, 0.0], [0.0, tiny]]);
        // det = tiny^2 > 0
        assert_eq!(m.det_sign_exact(), 1);
    }

    #[test]
    #[should_panic(expected = "non-finite matrix entry")]
    fn det_sign_exact_panics_on_nan() {
        let m = Matrix::<2>::from_rows([[f64::NAN, 0.0], [0.0, 1.0]]);
        let _ = m.det_sign_exact();
    }

    #[test]
    #[should_panic(expected = "non-finite matrix entry")]
    fn det_sign_exact_panics_on_infinity() {
        let m = Matrix::<2>::from_rows([[f64::INFINITY, 0.0], [0.0, 1.0]]);
        let _ = m.det_sign_exact();
    }

    #[test]
    fn det_sign_exact_pivot_needed_in_first_column() {
        // First column has zero on diagonal but non-zero below → needs pivoting.
        let m = Matrix::<3>::from_rows([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
        assert_eq!(m.det_sign_exact(), -1);
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
        assert_eq!(m.det_sign_exact(), 1);
    }
}
