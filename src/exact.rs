//! Exact arithmetic operations via arbitrary-precision rational numbers.
//!
//! This module is only compiled when the `"exact"` Cargo feature is enabled.
//!
//! # Architecture
//!
//! ## Determinants
//!
//! All three determinant methods (`det_exact`, `det_exact_f64`, `det_sign_exact`)
//! share the same integer-only Bareiss core (`bareiss_det_int`).  Each f64
//! entry is decomposed via `f64_decompose` into `mantissa × 2^exponent`,
//! all entries are scaled to a common `BigInt` matrix (shifting by
//! `e - e_min`), and Bareiss elimination runs entirely in `BigInt`
//! arithmetic — no `BigRational`, no GCD, no denominator tracking.
//! The result is `(det_int, total_exp)` where `det = det_int × 2^(D × e_min)`.
//! `bareiss_det` wraps this with `bigint_exp_to_bigrational` to reconstruct
//! a reduced `BigRational`; `det_sign_exact` reads the sign directly from
//! `det_int` (the scale factor is always positive).
//!
//! `det_sign_exact` adds a two-stage adaptive-precision optimisation inspired
//! by Shewchuk's robust geometric predicates:
//!
//! 1. **Fast filter (D ≤ 4)**: compute `det_direct()` and a conservative error
//!    bound. If `|det| > bound`, the f64 sign is provably correct — return
//!    immediately without allocating.
//! 2. **Exact fallback**: run integer-only Bareiss for a guaranteed-correct
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

use std::array::from_fn;

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

/// Decompose a finite `f64` into its IEEE 754 components.
///
/// Returns `None` for ±0.0, or `Some((mantissa, exponent, is_negative))` where
/// the value is exactly `(-1)^is_negative × mantissa × 2^exponent` and
/// `mantissa` is odd (trailing zeros stripped).  See `REFERENCES.md` \[9-10\].
///
/// # Panics
/// Panics if `x` is NaN or infinite.
fn f64_decompose(x: f64) -> Option<(u64, i32, bool)> {
    let bits = x.to_bits();
    let biased_exp = ((bits >> 52) & 0x7FF) as i32;
    let fraction = bits & 0x000F_FFFF_FFFF_FFFF;

    // ±0.0
    if biased_exp == 0 && fraction == 0 {
        return None;
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

    // Strip trailing zeros so the mantissa is odd.
    let tz = mantissa.trailing_zeros();
    let mantissa = mantissa >> tz;
    let exponent = raw_exp + tz.cast_signed();
    let is_negative = bits >> 63 != 0;

    Some((mantissa, exponent, is_negative))
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
    let Some((mantissa, exponent, is_negative)) = f64_decompose(x) else {
        return BigRational::from_integer(BigInt::from(0));
    };

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

/// Convert a `BigInt × 2^exp` pair to a reduced `BigRational`.
///
/// When `exp < 0` (denominator is `2^(-exp)`), shared factors of 2 are
/// stripped from `value` to keep the fraction in lowest terms without a
/// full GCD computation.
fn bigint_exp_to_bigrational(mut value: BigInt, mut exp: i32) -> BigRational {
    if value == BigInt::from(0) {
        return BigRational::from_integer(BigInt::from(0));
    }

    // Strip shared powers of 2 between value and the 2^(-exp) denominator.
    if exp < 0
        && let Some(tz) = value.trailing_zeros()
    {
        let reduce = tz.min(u64::from((-exp).cast_unsigned()));
        value >>= reduce;
        exp += i32::try_from(reduce).expect("reduce ≤ -exp which fits in i32");
    }

    if exp >= 0 {
        BigRational::new_raw(value << exp.cast_unsigned(), BigInt::from(1u32))
    } else {
        BigRational::new_raw(value, BigInt::from(1u32) << (-exp).cast_unsigned())
    }
}

/// Compute the exact determinant using integer-only Bareiss elimination.
///
/// Returns `(det_int, scale_exp)` where the true determinant is
/// `det_int × 2^scale_exp`.  Since the scale factor `2^scale_exp` is always
/// positive, `det_int.sign()` gives the sign of the determinant directly.
///
/// All arithmetic is in `BigInt` — no `BigRational`, no GCD, no denominator
/// tracking.  Each f64 entry is decomposed into `mantissa × 2^exponent` and
/// scaled to a common base `2^e_min` so every entry becomes an integer.
/// The Bareiss inner-loop division is exact (guaranteed by the algorithm).
fn bareiss_det_int<const D: usize>(m: &Matrix<D>) -> (BigInt, i32) {
    if D == 0 {
        return (BigInt::from(1), 0);
    }
    if D == 1 {
        return match f64_decompose(m.rows[0][0]) {
            None => (BigInt::from(0), 0),
            Some((mant, exp, neg)) => {
                let v = if neg {
                    -BigInt::from(mant)
                } else {
                    BigInt::from(mant)
                };
                (v, exp)
            }
        };
    }

    // Decompose all entries and find the minimum exponent.
    let mut components = [[(0u64, 0i32, false); D]; D];
    let mut e_min = i32::MAX;

    for (r, row) in m.rows.iter().enumerate() {
        for (c, &entry) in row.iter().enumerate() {
            if let Some((mant, exp, neg)) = f64_decompose(entry) {
                components[r][c] = (mant, exp, neg);
                e_min = e_min.min(exp);
            }
            // Zero entries keep the default (0, 0, false); their exponent is
            // excluded from e_min.
        }
    }

    // All entries are zero → singular.
    if e_min == i32::MAX {
        return (BigInt::from(0), 0);
    }

    // Build the integer matrix: a[r][c] = (±mantissa) × 2^(exp − e_min).
    let mut a: [[BigInt; D]; D] = from_fn(|r| {
        from_fn(|c| {
            let (mant, exp, neg) = components[r][c];
            if mant == 0 {
                BigInt::from(0)
            } else {
                let shift = (exp - e_min).cast_unsigned();
                let v = BigInt::from(mant) << shift;
                if neg { -v } else { v }
            }
        })
    });

    // Bareiss elimination in BigInt.
    let zero = BigInt::from(0);
    let mut prev_pivot = BigInt::from(1);
    let mut sign: i8 = 1;

    for k in 0..D {
        // Pivot search.
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
                return (BigInt::from(0), 0);
            }
        }

        // Elimination.
        for i in (k + 1)..D {
            for j in (k + 1)..D {
                a[i][j] = (&a[k][k] * &a[i][j] - &a[i][k] * &a[k][j]) / &prev_pivot;
            }
            a[i][k].clone_from(&zero);
        }

        prev_pivot.clone_from(&a[k][k]);
    }

    let det_int = if sign < 0 {
        -&a[D - 1][D - 1]
    } else {
        a[D - 1][D - 1].clone()
    };

    // det(original) = det_int × 2^(D × e_min)
    let d_i32 = i32::try_from(D).expect("dimension exceeds i32");
    let total_exp = e_min
        .checked_mul(d_i32)
        .expect("exponent overflow in bareiss_det_int");

    (det_int, total_exp)
}

/// Compute the exact determinant of a `D×D` matrix using integer-only Bareiss
/// elimination and return the result as a `BigRational`.
fn bareiss_det<const D: usize>(m: &Matrix<D>) -> BigRational {
    let (det_int, total_exp) = bareiss_det_int(m);
    bigint_exp_to_bigrational(det_int, total_exp)
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
    let mut mat: [[BigRational; D]; D] = from_fn(|r| from_fn(|c| f64_to_bigrational(m.rows[r][c])));
    let mut rhs: [BigRational; D] = from_fn(|r| f64_to_bigrational(b.data[r]));

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
    let mut x: [BigRational; D] = from_fn(|_| zero.clone());
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
    /// immediately without allocating.  Otherwise (and always for D ≥ 5),
    /// integer-only Bareiss elimination (`bareiss_det_int`) computes the exact
    /// sign without constructing any `BigRational` values.
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

        // Stage 2: integer Bareiss fallback — the 2^(D×e_min) scale factor
        // is always positive, so det_int.sign() == det(A).sign().
        let (det_int, _) = bareiss_det_int(self);
        Ok(match det_int.sign() {
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

    // -----------------------------------------------------------------------
    // f64_decompose tests
    // -----------------------------------------------------------------------

    #[test]
    fn f64_decompose_zero() {
        assert!(f64_decompose(0.0).is_none());
        assert!(f64_decompose(-0.0).is_none());
    }

    #[test]
    fn f64_decompose_one() {
        let (mant, exp, neg) = f64_decompose(1.0).unwrap();
        assert_eq!(mant, 1);
        assert_eq!(exp, 0);
        assert!(!neg);
    }

    #[test]
    fn f64_decompose_negative() {
        let (mant, exp, neg) = f64_decompose(-3.5).unwrap();
        // -3.5 = -7 × 2^(-1), mantissa is 7 (odd after stripping)
        assert_eq!(mant, 7);
        assert_eq!(exp, -1);
        assert!(neg);
    }

    #[test]
    fn f64_decompose_subnormal() {
        let tiny = 5e-324_f64;
        assert!(tiny.is_subnormal());
        let (mant, exp, neg) = f64_decompose(tiny).unwrap();
        assert_eq!(mant, 1);
        assert_eq!(exp, -1074);
        assert!(!neg);
    }

    #[test]
    fn f64_decompose_power_of_two() {
        let (mant, exp, neg) = f64_decompose(1024.0).unwrap();
        assert_eq!(mant, 1);
        assert_eq!(exp, 10); // 1024 = 2^10
        assert!(!neg);
    }

    #[test]
    #[should_panic(expected = "non-finite f64 in exact conversion")]
    fn f64_decompose_panics_on_nan() {
        f64_decompose(f64::NAN);
    }

    // -----------------------------------------------------------------------
    // bareiss_det_int tests
    // -----------------------------------------------------------------------

    #[test]
    fn bareiss_det_int_d0() {
        let (det, exp) = bareiss_det_int(&Matrix::<0>::zero());
        assert_eq!(det, BigInt::from(1));
        assert_eq!(exp, 0);
    }

    #[test]
    fn bareiss_det_int_d1_value() {
        // 7.0 = 7 × 2^0
        let (det, exp) = bareiss_det_int(&Matrix::<1>::from_rows([[7.0]]));
        assert_eq!(det, BigInt::from(7));
        assert_eq!(exp, 0);
    }

    #[test]
    fn bareiss_det_int_d1_zero() {
        let (det, _) = bareiss_det_int(&Matrix::<1>::from_rows([[0.0]]));
        assert_eq!(det, BigInt::from(0));
    }

    #[test]
    fn bareiss_det_int_d2_known() {
        // det([[1,2],[3,4]]) = -2
        let m = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let (det_int, total_exp) = bareiss_det_int(&m);
        // Reconstruct and verify.
        let det = bigint_exp_to_bigrational(det_int, total_exp);
        assert_eq!(det, BigRational::from_integer(BigInt::from(-2)));
    }

    #[test]
    fn bareiss_det_int_all_zeros() {
        let (det, _) = bareiss_det_int(&Matrix::<3>::zero());
        assert_eq!(det, BigInt::from(0));
    }

    #[test]
    fn bareiss_det_int_sign_matches_det_sign_exact() {
        // The sign of det_int should match det_sign_exact for various matrices.
        let m = Matrix::<3>::from_rows([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
        let (det_int, _) = bareiss_det_int(&m);
        assert_eq!(det_int.sign(), Sign::Minus); // det = -1
    }

    #[test]
    fn bareiss_det_int_fractional_entries() {
        // Entries with negative exponents: 0.5 = 1×2^(-1), 0.25 = 1×2^(-2).
        // det([[0.5, 0.25], [1.0, 1.0]]) = 0.5×1.0 − 0.25×1.0 = 0.25
        let m = Matrix::<2>::from_rows([[0.5, 0.25], [1.0, 1.0]]);
        let (det_int, total_exp) = bareiss_det_int(&m);
        let det = bigint_exp_to_bigrational(det_int, total_exp);
        assert_eq!(det, BigRational::new(BigInt::from(1), BigInt::from(4)));
    }

    #[test]
    fn bareiss_det_int_d1_negative() {
        // -3.5 = -7 × 2^(-1)
        let (det, exp) = bareiss_det_int(&Matrix::<1>::from_rows([[-3.5]]));
        assert_eq!(det, BigInt::from(-7));
        assert_eq!(exp, -1);
    }

    #[test]
    fn bareiss_det_int_d1_fractional() {
        // 0.5 = 1 × 2^(-1)
        let (det, exp) = bareiss_det_int(&Matrix::<1>::from_rows([[0.5]]));
        assert_eq!(det, BigInt::from(1));
        assert_eq!(exp, -1);
    }

    #[test]
    fn bareiss_det_int_d3_with_pivoting() {
        // Zero on diagonal → exercises pivot swap inside bareiss_det_int.
        let m = Matrix::<3>::from_rows([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
        let (det_int, total_exp) = bareiss_det_int(&m);
        let det = bigint_exp_to_bigrational(det_int, total_exp);
        assert_eq!(det, BigRational::from_integer(BigInt::from(-1)));
    }

    /// Per AGENTS.md: dimension-generic tests must cover D=2–5.
    macro_rules! gen_bareiss_det_int_identity_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<bareiss_det_int_identity_ $d d>]() {
                    let (det_int, total_exp) = bareiss_det_int(&Matrix::<$d>::identity());
                    let det = bigint_exp_to_bigrational(det_int, total_exp);
                    assert_eq!(det, BigRational::from_integer(BigInt::from(1)));
                }
            }
        };
    }

    gen_bareiss_det_int_identity_tests!(2);
    gen_bareiss_det_int_identity_tests!(3);
    gen_bareiss_det_int_identity_tests!(4);
    gen_bareiss_det_int_identity_tests!(5);

    // -----------------------------------------------------------------------
    // bigint_exp_to_bigrational tests
    // -----------------------------------------------------------------------

    #[test]
    fn bigint_exp_to_bigrational_zero() {
        let r = bigint_exp_to_bigrational(BigInt::from(0), -50);
        assert_eq!(r, BigRational::from_integer(BigInt::from(0)));
    }

    #[test]
    fn bigint_exp_to_bigrational_positive_exp() {
        // 3 × 2^2 = 12
        let r = bigint_exp_to_bigrational(BigInt::from(3), 2);
        assert_eq!(r, BigRational::from_integer(BigInt::from(12)));
    }

    #[test]
    fn bigint_exp_to_bigrational_negative_exp_reduced() {
        // 6 × 2^(-2) = 6/4 → reduced to 3/2 (strip one shared factor of 2)
        let r = bigint_exp_to_bigrational(BigInt::from(6), -2);
        assert_eq!(*r.numer(), BigInt::from(3));
        assert_eq!(*r.denom(), BigInt::from(2));
    }

    #[test]
    fn bigint_exp_to_bigrational_negative_exp_already_odd() {
        // 3 × 2^(-2) = 3/4 (already in lowest terms since 3 is odd)
        let r = bigint_exp_to_bigrational(BigInt::from(3), -2);
        assert_eq!(*r.numer(), BigInt::from(3));
        assert_eq!(*r.denom(), BigInt::from(4));
    }

    #[test]
    fn bigint_exp_to_bigrational_negative_value() {
        // -5 × 2^1 = -10
        let r = bigint_exp_to_bigrational(BigInt::from(-5), 1);
        assert_eq!(r, BigRational::from_integer(BigInt::from(-10)));
    }

    #[test]
    fn bigint_exp_to_bigrational_negative_value_with_denominator() {
        // -3 × 2^(-2) = -3/4
        let r = bigint_exp_to_bigrational(BigInt::from(-3), -2);
        assert_eq!(*r.numer(), BigInt::from(-3));
        assert_eq!(*r.denom(), BigInt::from(4));
    }

    // -----------------------------------------------------------------------
    // bareiss_det (wrapper) tests
    // -----------------------------------------------------------------------

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
