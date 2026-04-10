//! Exact arithmetic operations via arbitrary-precision rational numbers.
//!
//! This module is only compiled when the `"exact"` Cargo feature is enabled.
//!
//! # Architecture
//!
//! ## Determinants
//!
//! All three determinant methods (`det_exact`, `det_exact_f64`, `det_sign_exact`)
//! share the same core: entries are converted losslessly to `BigRational` and
//! the Bareiss algorithm (fraction-free Gaussian elimination) computes the
//! exact determinant.
//!
//! `det_sign_exact` adds a two-stage adaptive-precision optimisation inspired
//! by Shewchuk's robust geometric predicates:
//!
//! 1. **Fast filter (D ≤ 4)**: compute `det_direct()` and a conservative error
//!    bound. If `|det| > bound`, the f64 sign is provably correct — return
//!    immediately without allocating.
//! 2. **Exact fallback**: run the Bareiss algorithm for a guaranteed-correct
//!    sign.
//!
//! ## Linear system solve
//!
//! `solve_exact` and `solve_exact_f64` solve `A x = b` using Gaussian
//! elimination with first-non-zero pivoting in `BigRational` arithmetic.
//! Since all arithmetic is exact, any non-zero pivot gives the correct result
//! (there is no numerical stability concern).  Every finite `f64` is exactly
//! representable as a rational, so the result is provably correct.
//!
//! ## f64 → `BigRational` conversion
//!
//! All entry conversions use `f64_to_bigrational`, which decomposes the
//! IEEE 754 binary64 bit representation (\[9\]) into sign, exponent, and
//! significand and constructs a `BigRational` directly — avoiding the GCD
//! normalization that `BigRational::from_float` performs.  See Goldberg
//! \[10\] for background on floating-point representation and exact
//! rational reconstruction.  Reference numbers refer to `REFERENCES.md`.

use num_bigint::{BigInt, Sign};
use num_rational::BigRational;
use num_traits::ToPrimitive;

use crate::LaError;
use crate::matrix::Matrix;
use crate::vector::Vector;

/// Validate that all entries in a `D×D` matrix are finite (not NaN or infinite).
///
/// Returns `Ok(())` if all entries are finite, or `Err(LaError::NonFinite)` with
/// the column of the first non-finite entry found.
fn validate_finite<const D: usize>(m: &Matrix<D>) -> Result<(), LaError> {
    for r in 0..D {
        for c in 0..D {
            if !m.rows[r][c].is_finite() {
                return Err(LaError::NonFinite {
                    row: Some(r),
                    col: c,
                });
            }
        }
    }
    Ok(())
}

/// Validate that all entries in a length-`D` vector are finite.
///
/// Returns `Ok(())` if all entries are finite, or `Err(LaError::NonFinite)` with
/// the index of the first non-finite entry found.
fn validate_finite_vec<const D: usize>(v: &Vector<D>) -> Result<(), LaError> {
    for (i, &x) in v.data.iter().enumerate() {
        if !x.is_finite() {
            return Err(LaError::NonFinite { row: None, col: i });
        }
    }
    Ok(())
}

/// Convert an `f64` to an exact `BigRational` via IEEE 754 bit decomposition.
///
/// Every finite `f64` is exactly representable as `±m × 2^e` where `m` is a
/// non-negative integer and `e` is an integer.  This function extracts `(m, e)`
/// directly from the IEEE 754 binary64 bit layout \[9\], strips trailing zeros
/// from `m` so the resulting fraction is already in lowest terms, then
/// constructs a `BigRational` via `new_raw` — bypassing the GCD reduction
/// that `BigRational::from_float` performs internally.
///
/// See `REFERENCES.md` \[9-10\] for the IEEE 754 standard and Goldberg's
/// survey of floating-point representation.
///
/// # Panics
/// Panics if `x` is NaN or infinite.
fn f64_to_bigrational(x: f64) -> BigRational {
    let bits = x.to_bits();
    let biased_exp = ((bits >> 52) & 0x7FF) as i32;
    let fraction = bits & 0x000F_FFFF_FFFF_FFFF;

    // ±0.0
    if biased_exp == 0 && fraction == 0 {
        return BigRational::from_integer(BigInt::from(0));
    }

    // NaN / Inf — callers must validate finiteness before reaching here.
    assert!(biased_exp != 0x7FF, "non-finite f64 in exact conversion");

    let (mantissa, raw_exp) = if biased_exp == 0 {
        // Subnormal: (-1)^s × 0.fraction × 2^(-1022)
        //          = (-1)^s × fraction × 2^(-1074)
        (fraction, -1074_i32)
    } else {
        // Normal: (-1)^s × 1.fraction × 2^(biased_exp - 1023)
        //       = (-1)^s × (2^52 | fraction) × 2^(biased_exp - 1075)
        ((1u64 << 52) | fraction, biased_exp - 1075)
    };

    // Strip trailing zeros so the fraction is already in lowest terms:
    // after stripping, mantissa is odd and the denominator (if any) is a
    // power of 2, so gcd(mantissa, denom) = 1.
    let tz = mantissa.trailing_zeros();
    let mantissa = mantissa >> tz;
    let exponent = raw_exp + tz.cast_signed();

    let is_negative = bits >> 63 != 0;
    let numer = if is_negative {
        -BigInt::from(mantissa)
    } else {
        BigInt::from(mantissa)
    };

    if exponent >= 0 {
        BigRational::new_raw(numer << exponent.cast_unsigned(), BigInt::from(1u32))
    } else {
        BigRational::new_raw(numer, BigInt::from(1u32) << (-exponent).cast_unsigned())
    }
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
    let mut a: [[BigRational; D]; D] =
        std::array::from_fn(|r| std::array::from_fn(|c| f64_to_bigrational(m.rows[r][c])));

    let zero = BigRational::from_integer(BigInt::from(0));
    let mut prev_pivot = BigRational::from_integer(BigInt::from(1));
    let mut sign: i8 = 1;

    for k in 0..D {
        // Find first non-zero entry in column k at or below row k.
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

/// Solve `A x = b` using Gaussian elimination with first-non-zero pivoting
/// in `BigRational` arithmetic.
///
/// Since all arithmetic is exact, any non-zero pivot gives the correct result
/// (no numerical stability concern).  This matches the pivoting strategy used
/// by `bareiss_det`.
///
/// Returns the exact solution as `[BigRational; D]`.
/// Returns `Err(LaError::Singular)` if the matrix is exactly singular.
fn gauss_solve<const D: usize>(m: &Matrix<D>, b: &Vector<D>) -> Result<[BigRational; D], LaError> {
    let zero = BigRational::from_integer(BigInt::from(0));

    // Build matrix and RHS separately (cannot use [BigRational; D+1] augmented
    // columns because const-generic expressions are unstable).
    let mut mat: [[BigRational; D]; D] =
        std::array::from_fn(|r| std::array::from_fn(|c| f64_to_bigrational(m.rows[r][c])));
    let mut rhs: [BigRational; D] = std::array::from_fn(|r| f64_to_bigrational(b.data[r]));

    // Forward elimination with first-non-zero pivoting.
    for k in 0..D {
        // Find first non-zero pivot in column k at or below row k.
        if mat[k][k] == zero {
            if let Some(swap_row) = ((k + 1)..D).find(|&i| mat[i][k] != zero) {
                mat.swap(k, swap_row);
                rhs.swap(k, swap_row);
            } else {
                return Err(LaError::Singular { pivot_col: k });
            }
        }

        // Eliminate below pivot.
        let pivot = mat[k][k].clone();
        for i in (k + 1)..D {
            if mat[i][k] != zero {
                let factor = &mat[i][k] / &pivot;
                // We need index `j` to read mat[k][j] and write mat[i][j]
                // (two distinct rows) — iterators can't borrow both.
                #[allow(clippy::needless_range_loop)]
                for j in (k + 1)..D {
                    let term = &factor * &mat[k][j];
                    mat[i][j] -= term;
                }
                let rhs_term = &factor * &rhs[k];
                rhs[i] -= rhs_term;
                mat[i][k] = zero.clone();
            }
        }
    }

    // Back-substitution.
    let mut x: [BigRational; D] = std::array::from_fn(|_| zero.clone());
    for i in (0..D).rev() {
        let mut sum = rhs[i].clone();
        for j in (i + 1)..D {
            sum -= &mat[i][j] * &x[j];
        }
        x[i] = sum / &mat[i][i];
    }

    Ok(x)
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
    /// Returns [`LaError::NonFinite`] if any matrix entry is NaN or infinite.
    /// Returns [`LaError::Overflow`] if the exact determinant is too large to
    /// represent as a finite `f64`.
    #[inline]
    pub fn det_exact_f64(&self) -> Result<f64, LaError> {
        let exact = self.det_exact()?;
        let val = exact.to_f64().unwrap_or(f64::INFINITY);
        if val.is_finite() {
            Ok(val)
        } else {
            Err(LaError::Overflow { index: None })
        }
    }

    /// Exact linear system solve using arbitrary-precision rational arithmetic.
    ///
    /// Solves `A x = b` where `A` is `self` and `b` is the given vector.
    /// Returns the exact solution as `[BigRational; D]`.  Every finite `f64`
    /// is exactly representable as a rational, so the conversion is lossless
    /// and the result is provably correct.
    ///
    /// # When to use
    ///
    /// Use this when you need a provably correct solution — for example,
    /// exact circumcenter computation for near-degenerate simplices where
    /// f64 arithmetic may produce wildly wrong results.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// // A x = b  where A = [[1,2],[3,4]], b = [5, 11]  →  x = [1, 2]
    /// let a = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
    /// let b = Vector::<2>::new([5.0, 11.0]);
    /// let x = a.solve_exact(b).unwrap();
    /// assert_eq!(x[0], BigRational::from_integer(1.into()));
    /// assert_eq!(x[1], BigRational::from_integer(2.into()));
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] if any matrix or vector entry is NaN or
    /// infinite.
    /// Returns [`LaError::Singular`] if the matrix is exactly singular.
    #[inline]
    pub fn solve_exact(&self, b: Vector<D>) -> Result<[BigRational; D], LaError> {
        validate_finite(self)?;
        validate_finite_vec(&b)?;
        gauss_solve(self, &b)
    }

    /// Exact linear system solve converted to `f64`.
    ///
    /// Computes the exact [`BigRational`] solution via
    /// [`solve_exact`](Self::solve_exact) and converts each component to the
    /// nearest `f64`.  This is useful when you want the most accurate f64
    /// solution possible without committing to `BigRational` downstream.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// let a = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
    /// let b = Vector::<2>::new([5.0, 11.0]);
    /// let x = a.solve_exact_f64(b).unwrap().into_array();
    /// assert!((x[0] - 1.0).abs() <= f64::EPSILON);
    /// assert!((x[1] - 2.0).abs() <= f64::EPSILON);
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] if any matrix or vector entry is NaN or
    /// infinite.
    /// Returns [`LaError::Singular`] if the matrix is exactly singular.
    /// Returns [`LaError::Overflow`] if any component of the exact solution is
    /// too large to represent as a finite `f64`.
    #[inline]
    pub fn solve_exact_f64(&self, b: Vector<D>) -> Result<Vector<D>, LaError> {
        let exact = self.solve_exact(b)?;
        let mut result = [0.0f64; D];
        for (i, val) in exact.iter().enumerate() {
            let f = val.to_f64().unwrap_or(f64::INFINITY);
            if !f.is_finite() {
                return Err(LaError::Overflow { index: Some(i) });
            }
            result[i] = f;
        }
        Ok(Vector::new(result))
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

    use pastey::paste;

    // -----------------------------------------------------------------------
    // Macro-generated per-dimension tests (D=2..5)
    // -----------------------------------------------------------------------

    macro_rules! gen_det_exact_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<det_exact_identity_ $d d>]() {
                    let det = Matrix::<$d>::identity().det_exact().unwrap();
                    assert_eq!(det, BigRational::from_integer(BigInt::from(1)));
                }

                #[test]
                fn [<det_exact_err_on_nan_ $d d>]() {
                    let mut m = Matrix::<$d>::identity();
                    m.set(0, 0, f64::NAN);
                    assert_eq!(m.det_exact(), Err(LaError::NonFinite { row: Some(0), col: 0 }));
                }

                #[test]
                fn [<det_exact_err_on_inf_ $d d>]() {
                    let mut m = Matrix::<$d>::identity();
                    m.set(0, 0, f64::INFINITY);
                    assert_eq!(m.det_exact(), Err(LaError::NonFinite { row: Some(0), col: 0 }));
                }
            }
        };
    }

    gen_det_exact_tests!(2);
    gen_det_exact_tests!(3);
    gen_det_exact_tests!(4);
    gen_det_exact_tests!(5);

    macro_rules! gen_det_exact_f64_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<det_exact_f64_identity_ $d d>]() {
                    let det = Matrix::<$d>::identity().det_exact_f64().unwrap();
                    assert!((det - 1.0).abs() <= f64::EPSILON);
                }

                #[test]
                fn [<det_exact_f64_err_on_nan_ $d d>]() {
                    let mut m = Matrix::<$d>::identity();
                    m.set(0, 0, f64::NAN);
                    assert_eq!(m.det_exact_f64(), Err(LaError::NonFinite { row: Some(0), col: 0 }));
                }
            }
        };
    }

    gen_det_exact_f64_tests!(2);
    gen_det_exact_f64_tests!(3);
    gen_det_exact_f64_tests!(4);
    gen_det_exact_f64_tests!(5);

    /// For D ≤ 4, `det_exact_f64` should agree with `det_direct` on
    /// well-conditioned matrices.
    macro_rules! gen_det_exact_f64_agrees_with_det_direct {
        ($d:literal) => {
            paste! {
                #[test]
                #[allow(clippy::cast_precision_loss)]
                fn [<det_exact_f64_agrees_with_det_direct_ $d d>]() {
                    // Diagonally dominant → well-conditioned.
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
                    let exact = m.det_exact_f64().unwrap();
                    let direct = m.det_direct().unwrap();
                    let eps = direct.abs().mul_add(1e-12, 1e-12);
                    assert!((exact - direct).abs() <= eps);
                }
            }
        };
    }

    gen_det_exact_f64_agrees_with_det_direct!(2);
    gen_det_exact_f64_agrees_with_det_direct!(3);
    gen_det_exact_f64_agrees_with_det_direct!(4);

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
        assert_eq!(
            m.det_sign_exact(),
            Err(LaError::NonFinite {
                row: Some(0),
                col: 0
            })
        );
    }

    #[test]
    fn det_sign_exact_returns_err_on_infinity() {
        let m = Matrix::<2>::from_rows([[f64::INFINITY, 0.0], [0.0, 1.0]]);
        assert_eq!(
            m.det_sign_exact(),
            Err(LaError::NonFinite {
                row: Some(0),
                col: 0
            })
        );
    }

    #[test]
    fn det_sign_exact_returns_err_on_nan_5x5() {
        // D ≥ 5 bypasses the fast filter, exercising the bareiss_det path.
        let mut m = Matrix::<5>::identity();
        m.set(2, 3, f64::NAN);
        assert_eq!(
            m.det_sign_exact(),
            Err(LaError::NonFinite {
                row: Some(2),
                col: 3
            })
        );
    }

    #[test]
    fn det_sign_exact_returns_err_on_infinity_5x5() {
        let mut m = Matrix::<5>::identity();
        m.set(0, 0, f64::INFINITY);
        assert_eq!(
            m.det_sign_exact(),
            Err(LaError::NonFinite {
                row: Some(0),
                col: 0
            })
        );
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
    // det_exact: dimension-specific tests
    // -----------------------------------------------------------------------

    #[test]
    fn det_exact_d0_is_one() {
        let det = Matrix::<0>::zero().det_exact().unwrap();
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

    // -----------------------------------------------------------------------
    // det_exact_f64: dimension-specific tests
    // -----------------------------------------------------------------------

    #[test]
    fn det_exact_f64_known_2x2() {
        let m = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let det = m.det_exact_f64().unwrap();
        assert!((det - (-2.0)).abs() <= f64::EPSILON);
    }

    #[test]
    fn det_exact_f64_overflow_returns_err() {
        // Entries near f64::MAX produce a determinant too large for f64.
        let big = f64::MAX / 2.0;
        let m = Matrix::<3>::from_rows([[0.0, 0.0, 1.0], [big, 0.0, 1.0], [0.0, big, 1.0]]);
        // det = big^2, which overflows f64.
        assert_eq!(m.det_exact_f64(), Err(LaError::Overflow { index: None }));
    }

    // -----------------------------------------------------------------------
    // solve_exact: macro-generated per-dimension tests (D=2..5)
    // -----------------------------------------------------------------------

    /// Helper: build an arbitrary RHS vector for dimension `$d`.
    fn arbitrary_rhs<const D: usize>() -> Vector<D> {
        let values = [1.0, -2.5, 3.0, 0.25, -4.0];
        let mut arr = [0.0f64; D];
        for (dst, src) in arr.iter_mut().zip(values.iter()) {
            *dst = *src;
        }
        Vector::<D>::new(arr)
    }

    macro_rules! gen_solve_exact_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<solve_exact_identity_ $d d>]() {
                    let a = Matrix::<$d>::identity();
                    let b = arbitrary_rhs::<$d>();
                    let x = a.solve_exact(b).unwrap();
                    for (i, xi) in x.iter().enumerate() {
                        assert_eq!(*xi, f64_to_bigrational(b.data[i]));
                    }
                }

                #[test]
                fn [<solve_exact_err_on_nan_matrix_ $d d>]() {
                    let mut a = Matrix::<$d>::identity();
                    a.set(0, 0, f64::NAN);
                    let b = arbitrary_rhs::<$d>();
                    assert_eq!(a.solve_exact(b), Err(LaError::NonFinite { row: Some(0), col: 0 }));
                }

                #[test]
                fn [<solve_exact_err_on_inf_matrix_ $d d>]() {
                    let mut a = Matrix::<$d>::identity();
                    a.set(0, 0, f64::INFINITY);
                    let b = arbitrary_rhs::<$d>();
                    assert_eq!(a.solve_exact(b), Err(LaError::NonFinite { row: Some(0), col: 0 }));
                }

                #[test]
                fn [<solve_exact_err_on_nan_vector_ $d d>]() {
                    let a = Matrix::<$d>::identity();
                    let mut b_arr = [1.0f64; $d];
                    b_arr[0] = f64::NAN;
                    let b = Vector::<$d>::new(b_arr);
                    assert_eq!(a.solve_exact(b), Err(LaError::NonFinite { row: None, col: 0 }));
                }

                #[test]
                fn [<solve_exact_err_on_inf_vector_ $d d>]() {
                    let a = Matrix::<$d>::identity();
                    let mut b_arr = [1.0f64; $d];
                    b_arr[$d - 1] = f64::INFINITY;
                    let b = Vector::<$d>::new(b_arr);
                    assert_eq!(a.solve_exact(b), Err(LaError::NonFinite { row: None, col: $d - 1 }));
                }

                #[test]
                fn [<solve_exact_singular_ $d d>]() {
                    // Zero matrix is singular.
                    let a = Matrix::<$d>::zero();
                    let b = arbitrary_rhs::<$d>();
                    assert_eq!(a.solve_exact(b), Err(LaError::Singular { pivot_col: 0 }));
                }
            }
        };
    }

    gen_solve_exact_tests!(2);
    gen_solve_exact_tests!(3);
    gen_solve_exact_tests!(4);
    gen_solve_exact_tests!(5);

    macro_rules! gen_solve_exact_f64_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<solve_exact_f64_identity_ $d d>]() {
                    let a = Matrix::<$d>::identity();
                    let b = arbitrary_rhs::<$d>();
                    let x = a.solve_exact_f64(b).unwrap().into_array();
                    for i in 0..$d {
                        assert!((x[i] - b.data[i]).abs() <= f64::EPSILON);
                    }
                }

                #[test]
                fn [<solve_exact_f64_err_on_nan_ $d d>]() {
                    let mut a = Matrix::<$d>::identity();
                    a.set(0, 0, f64::NAN);
                    let b = arbitrary_rhs::<$d>();
                    assert_eq!(a.solve_exact_f64(b), Err(LaError::NonFinite { row: Some(0), col: 0 }));
                }
            }
        };
    }

    gen_solve_exact_f64_tests!(2);
    gen_solve_exact_f64_tests!(3);
    gen_solve_exact_f64_tests!(4);
    gen_solve_exact_f64_tests!(5);

    /// For D ≤ 4, `solve_exact_f64` should agree with `Lu::solve_vec` on
    /// well-conditioned matrices.
    macro_rules! gen_solve_exact_f64_agrees_with_lu {
        ($d:literal) => {
            paste! {
                #[test]
                #[allow(clippy::cast_precision_loss)]
                fn [<solve_exact_f64_agrees_with_lu_ $d d>]() {
                    // Diagonally dominant → well-conditioned.
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
                    let a = Matrix::<$d>::from_rows(rows);
                    let b = arbitrary_rhs::<$d>();
                    let exact = a.solve_exact_f64(b).unwrap().into_array();
                    let lu_sol = a.lu(DEFAULT_PIVOT_TOL).unwrap()
                        .solve_vec(b).unwrap().into_array();
                    for i in 0..$d {
                        let eps = lu_sol[i].abs().mul_add(1e-12, 1e-12);
                        assert!((exact[i] - lu_sol[i]).abs() <= eps);
                    }
                }
            }
        };
    }

    gen_solve_exact_f64_agrees_with_lu!(2);
    gen_solve_exact_f64_agrees_with_lu!(3);
    gen_solve_exact_f64_agrees_with_lu!(4);
    gen_solve_exact_f64_agrees_with_lu!(5);

    // -----------------------------------------------------------------------
    // solve_exact: dimension-specific tests
    // -----------------------------------------------------------------------

    #[test]
    fn solve_exact_d0_returns_empty() {
        let a = Matrix::<0>::zero();
        let b = Vector::<0>::zero();
        let x = a.solve_exact(b).unwrap();
        assert!(x.is_empty());
    }

    #[test]
    fn solve_exact_known_2x2() {
        // [[1,2],[3,4]] x = [5, 11] → x = [1, 2]
        let a = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let b = Vector::<2>::new([5.0, 11.0]);
        let x = a.solve_exact(b).unwrap();
        assert_eq!(x[0], BigRational::from_integer(BigInt::from(1)));
        assert_eq!(x[1], BigRational::from_integer(BigInt::from(2)));
    }

    #[test]
    fn solve_exact_pivoting_needed() {
        // First column has zero on diagonal → pivot swap required.
        let a = Matrix::<3>::from_rows([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
        let b = Vector::<3>::new([2.0, 3.0, 4.0]);
        let x = a.solve_exact(b).unwrap();
        // x = [3, 2, 4]
        assert_eq!(x[0], f64_to_bigrational(3.0));
        assert_eq!(x[1], f64_to_bigrational(2.0));
        assert_eq!(x[2], f64_to_bigrational(4.0));
    }

    #[test]
    fn solve_exact_fractional_result() {
        // [[2, 1], [1, 3]] x = [1, 1] → x = [2/5, 1/5]
        let a = Matrix::<2>::from_rows([[2.0, 1.0], [1.0, 3.0]]);
        let b = Vector::<2>::new([1.0, 1.0]);
        let x = a.solve_exact(b).unwrap();
        assert_eq!(x[0], BigRational::new(BigInt::from(2), BigInt::from(5)));
        assert_eq!(x[1], BigRational::new(BigInt::from(1), BigInt::from(5)));
    }

    #[test]
    fn solve_exact_singular_duplicate_rows() {
        let a = Matrix::<3>::from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [1.0, 2.0, 3.0]]);
        let b = Vector::<3>::new([1.0, 2.0, 3.0]);
        assert!(matches!(a.solve_exact(b), Err(LaError::Singular { .. })));
    }

    #[test]
    fn solve_exact_5x5_permutation() {
        // Permutation matrix (swap rows 0↔1): P x = b → x = P^T b.
        let a = Matrix::<5>::from_rows([
            [0.0, 1.0, 0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0],
        ]);
        let b = Vector::<5>::new([10.0, 20.0, 30.0, 40.0, 50.0]);
        let x = a.solve_exact(b).unwrap();
        assert_eq!(x[0], f64_to_bigrational(20.0));
        assert_eq!(x[1], f64_to_bigrational(10.0));
        assert_eq!(x[2], f64_to_bigrational(30.0));
        assert_eq!(x[3], f64_to_bigrational(40.0));
        assert_eq!(x[4], f64_to_bigrational(50.0));
    }

    // -----------------------------------------------------------------------
    // solve_exact_f64: dimension-specific tests
    // -----------------------------------------------------------------------

    #[test]
    fn solve_exact_f64_known_2x2() {
        let a = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let b = Vector::<2>::new([5.0, 11.0]);
        let x = a.solve_exact_f64(b).unwrap().into_array();
        assert!((x[0] - 1.0).abs() <= f64::EPSILON);
        assert!((x[1] - 2.0).abs() <= f64::EPSILON);
    }

    #[test]
    fn solve_exact_f64_overflow_returns_err() {
        // [[1/big, 0], [0, 1/big]] x = [big, big] → x = [big², big²],
        // which overflows f64.
        let big = f64::MAX / 2.0;
        let a = Matrix::<2>::from_rows([[1.0 / big, 0.0], [0.0, 1.0 / big]]);
        let b = Vector::<2>::new([big, big]);
        assert_eq!(
            a.solve_exact_f64(b),
            Err(LaError::Overflow { index: Some(0) })
        );
    }

    // -----------------------------------------------------------------------
    // gauss_solve: internal helper tests
    // -----------------------------------------------------------------------

    #[test]
    fn gauss_solve_d0_returns_empty() {
        let a = Matrix::<0>::zero();
        let b = Vector::<0>::zero();
        assert_eq!(gauss_solve(&a, &b).unwrap().len(), 0);
    }

    #[test]
    fn gauss_solve_d1() {
        let a = Matrix::<1>::from_rows([[2.0]]);
        let b = Vector::<1>::new([6.0]);
        let x = gauss_solve(&a, &b).unwrap();
        assert_eq!(x[0], f64_to_bigrational(3.0));
    }

    #[test]
    fn gauss_solve_singular_column_all_zero() {
        let a = Matrix::<3>::from_rows([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
        let b = Vector::<3>::new([1.0, 2.0, 3.0]);
        assert_eq!(gauss_solve(&a, &b), Err(LaError::Singular { pivot_col: 1 }));
    }

    // -----------------------------------------------------------------------
    // f64_to_bigrational tests
    // -----------------------------------------------------------------------

    #[test]
    fn f64_to_bigrational_positive_zero() {
        let r = f64_to_bigrational(0.0);
        assert_eq!(r, BigRational::from_integer(BigInt::from(0)));
    }

    #[test]
    fn f64_to_bigrational_negative_zero() {
        let r = f64_to_bigrational(-0.0);
        assert_eq!(r, BigRational::from_integer(BigInt::from(0)));
    }

    #[test]
    fn f64_to_bigrational_one() {
        let r = f64_to_bigrational(1.0);
        assert_eq!(r, BigRational::from_integer(BigInt::from(1)));
    }

    #[test]
    fn f64_to_bigrational_negative_one() {
        let r = f64_to_bigrational(-1.0);
        assert_eq!(r, BigRational::from_integer(BigInt::from(-1)));
    }

    #[test]
    fn f64_to_bigrational_half() {
        let r = f64_to_bigrational(0.5);
        assert_eq!(r, BigRational::new(BigInt::from(1), BigInt::from(2)));
    }

    #[test]
    fn f64_to_bigrational_quarter() {
        let r = f64_to_bigrational(0.25);
        assert_eq!(r, BigRational::new(BigInt::from(1), BigInt::from(4)));
    }

    #[test]
    fn f64_to_bigrational_negative_three_and_a_half() {
        // -3.5 = -7/2
        let r = f64_to_bigrational(-3.5);
        assert_eq!(r, BigRational::new(BigInt::from(-7), BigInt::from(2)));
    }

    #[test]
    fn f64_to_bigrational_integer() {
        let r = f64_to_bigrational(42.0);
        assert_eq!(r, BigRational::from_integer(BigInt::from(42)));
    }

    #[test]
    fn f64_to_bigrational_power_of_two() {
        let r = f64_to_bigrational(1024.0);
        assert_eq!(r, BigRational::from_integer(BigInt::from(1024)));
    }

    #[test]
    fn f64_to_bigrational_subnormal() {
        let tiny = 5e-324_f64; // smallest positive subnormal
        assert!(tiny.is_subnormal());
        let r = f64_to_bigrational(tiny);
        // 5e-324 = 1 × 2^(-1074)
        assert_eq!(
            r,
            BigRational::new(BigInt::from(1), BigInt::from(1u32) << 1074u32)
        );
    }

    #[test]
    fn f64_to_bigrational_already_lowest_terms() {
        // 0.5 should produce numer=1, denom=2 (already reduced).
        let r = f64_to_bigrational(0.5);
        assert_eq!(*r.numer(), BigInt::from(1));
        assert_eq!(*r.denom(), BigInt::from(2));
    }

    #[test]
    fn f64_to_bigrational_round_trip() {
        // -0.0 is excluded: it maps to BigRational(0) which round-trips
        // to +0.0 (correct; tested separately in f64_to_bigrational_negative_zero).
        let values = [
            0.0,
            1.0,
            -1.0,
            0.5,
            0.25,
            0.1,
            42.0,
            -3.5,
            1e10,
            1e-10,
            f64::MAX / 2.0,
            f64::MIN_POSITIVE,
            5e-324,
        ];
        for &v in &values {
            let r = f64_to_bigrational(v);
            let back = r.to_f64().expect("round-trip to_f64 failed");
            assert!(
                v.to_bits() == back.to_bits(),
                "round-trip failed for {v}: got {back}"
            );
        }
    }

    #[test]
    #[should_panic(expected = "non-finite f64 in exact conversion")]
    fn f64_to_bigrational_panics_on_nan() {
        f64_to_bigrational(f64::NAN);
    }

    #[test]
    #[should_panic(expected = "non-finite f64 in exact conversion")]
    fn f64_to_bigrational_panics_on_inf() {
        f64_to_bigrational(f64::INFINITY);
    }

    #[test]
    #[should_panic(expected = "non-finite f64 in exact conversion")]
    fn f64_to_bigrational_panics_on_neg_inf() {
        f64_to_bigrational(f64::NEG_INFINITY);
    }

    // -----------------------------------------------------------------------
    // validate_finite_vec tests
    // -----------------------------------------------------------------------

    #[test]
    fn validate_finite_vec_ok() {
        assert!(validate_finite_vec(&Vector::<3>::new([1.0, 2.0, 3.0])).is_ok());
    }

    #[test]
    fn validate_finite_vec_err_on_nan() {
        assert_eq!(
            validate_finite_vec(&Vector::<2>::new([f64::NAN, 1.0])),
            Err(LaError::NonFinite { row: None, col: 0 })
        );
    }

    #[test]
    fn validate_finite_vec_err_on_inf() {
        assert_eq!(
            validate_finite_vec(&Vector::<2>::new([1.0, f64::NEG_INFINITY])),
            Err(LaError::NonFinite { row: None, col: 1 })
        );
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
        assert_eq!(
            validate_finite(&m),
            Err(LaError::NonFinite {
                row: Some(1),
                col: 0
            })
        );
    }

    #[test]
    fn validate_finite_err_on_inf() {
        let mut m = Matrix::<2>::identity();
        m.set(0, 1, f64::NEG_INFINITY);
        assert_eq!(
            validate_finite(&m),
            Err(LaError::NonFinite {
                row: Some(0),
                col: 1
            })
        );
    }
}
