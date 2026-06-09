#![forbid(unsafe_code)]
#![deny(missing_docs)]
#![doc = include_str!("../README.md")]

#[cfg(doc)]
mod readme_doctests {
    //! Executable versions of README examples.
    /// ```rust
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// // This system requires pivoting (a[0][0] = 0), so it's a good LU demo.
    /// let a = Matrix::<5>::try_from_rows([
    ///     [0.0, 1.0, 1.0, 1.0, 1.0],
    ///     [1.0, 0.0, 1.0, 1.0, 1.0],
    ///     [1.0, 1.0, 0.0, 1.0, 1.0],
    ///     [1.0, 1.0, 1.0, 0.0, 1.0],
    ///     [1.0, 1.0, 1.0, 1.0, 0.0],
    /// ])?;
    ///
    /// let b = Vector::<5>::try_new([14.0, 13.0, 12.0, 11.0, 10.0])?;
    ///
    /// let lu = a.lu(DEFAULT_SINGULAR_TOL)?;
    /// let x = lu.solve(b)?.into_array();
    ///
    /// // Floating-point rounding is expected; compare with a tolerance.
    /// let expected = [1.0, 2.0, 3.0, 4.0, 5.0];
    /// for (x_i, e_i) in x.iter().zip(expected.iter()) {
    ///     assert!((*x_i - *e_i).abs() <= 1e-12);
    /// }
    /// # Ok(())
    /// # }
    /// ```
    fn solve_5x5_example() {}

    /// ```rust
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// // This matrix is symmetric positive-definite (A = L*L^T) so LDLT works without pivoting.
    /// let a = Matrix::<5>::try_from_rows([
    ///     [1.0, 1.0, 0.0, 0.0, 0.0],
    ///     [1.0, 2.0, 1.0, 0.0, 0.0],
    ///     [0.0, 1.0, 2.0, 1.0, 0.0],
    ///     [0.0, 0.0, 1.0, 2.0, 1.0],
    ///     [0.0, 0.0, 0.0, 1.0, 2.0],
    /// ])?;
    ///
    /// let ldlt = match a.ldlt(DEFAULT_SINGULAR_TOL) {
    ///     Ok(ldlt) => ldlt,
    ///     Err(err @ LaError::Asymmetric { row, col, .. }) => {
    ///         eprintln!("LDLT requires symmetry; first mismatch at ({row}, {col})");
    ///         return Err(err);
    ///     }
    ///     Err(err) => return Err(err),
    /// };
    ///
    /// let det = ldlt.det()?;
    /// assert!((det - 1.0).abs() <= 1e-12);
    /// # Ok(())
    /// # }
    /// ```
    fn det_5x5_ldlt_example() {}

    /// ```rust
    /// use la_stack::prelude::*;
    ///
    /// // Evaluated entirely at compile time — no runtime cost.
    /// const DET: Result<Option<f64>, LaError> = match Matrix::<4>::try_from_rows([
    ///     [2.0, 0.0, 0.0, 0.0],
    ///     [0.0, 3.0, 0.0, 0.0],
    ///     [0.0, 0.0, 5.0, 0.0],
    ///     [0.0, 0.0, 0.0, 7.0],
    /// ]) {
    ///     Ok(matrix) => matrix.det_direct(),
    ///     Err(err) => Err(err),
    /// };
    ///
    /// # fn main() -> Result<(), LaError> {
    /// assert_eq!(DET?, Some(210.0));
    /// # Ok(())
    /// # }
    /// ```
    fn det_direct_4x4_const_example() {}

    #[cfg(feature = "exact")]
    /// ```rust
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// // Exact determinant
    /// let m = Matrix::<3>::try_from_rows([
    ///     [1.0, 2.0, 3.0],
    ///     [4.0, 5.0, 6.0],
    ///     [7.0, 8.0, 9.0],
    /// ])?;
    /// assert_eq!(m.det_sign_exact()?, 0); // exactly singular
    ///
    /// let det = m.det_exact()?;
    /// assert_eq!(det, BigRational::from_integer(0.into())); // exact zero
    /// let det_f64 = m.det_exact_f64()?;
    /// assert_eq!(det_f64, 0.0);
    ///
    /// // If strict exact-to-f64 conversion would require rounding, opt in
    /// // explicitly with the rounded API.
    /// let inexact = Matrix::<2>::try_from_rows([
    ///     [1.0 + f64::EPSILON, 0.0],
    ///     [0.0, 1.0 - f64::EPSILON],
    /// ])?;
    /// let rounded_det = match inexact.det_exact_f64() {
    ///     Ok(det) => det,
    ///     Err(err) if err.requires_rounding() => inexact.det_exact_rounded_f64()?,
    ///     Err(err) => return Err(err),
    /// };
    /// assert_eq!(rounded_det.to_bits(), 1.0f64.to_bits());
    ///
    /// // If the exact determinant cannot fit in f64, keep the BigRational value.
    /// let big = f64::MAX / 2.0;
    /// let huge = Matrix::<3>::try_from_rows([
    ///     [0.0, 0.0, 1.0],
    ///     [big, 0.0, 1.0],
    ///     [0.0, big, 1.0],
    /// ])?;
    /// let huge_det = huge.det_exact()?;
    /// assert_eq!(
    ///     huge.det_exact_f64()
    ///         .err()
    ///         .and_then(|err| err.unrepresentable_reason()),
    ///     Some(UnrepresentableReason::NotFinite)
    /// );
    /// println!("exact determinant = {huge_det}");
    ///
    /// // Exact linear system solve
    /// let a = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]])?;
    /// let b = Vector::<2>::try_new([5.0, 11.0])?;
    /// let x = a.solve_exact_f64(b)?.into_array();
    /// assert!((x[0] - 1.0).abs() <= f64::EPSILON);
    /// assert!((x[1] - 2.0).abs() <= f64::EPSILON);
    /// # Ok(())
    /// # }
    /// ```
    fn exact_arithmetic_example() {}
}

mod error;
#[cfg(feature = "exact")]
mod exact;
mod ldlt;
mod lu;
mod matrix;
mod tolerance;
mod vector;

#[cfg(feature = "exact")]
pub use num_bigint::BigInt;
#[cfg(feature = "exact")]
pub use num_rational::BigRational;
#[cfg(feature = "exact")]
pub use num_traits::{FromPrimitive, Signed, ToPrimitive};

// ---------------------------------------------------------------------------
// Error-bound constants for `Matrix::det_errbound()`.
//
// For `D ∈ {2, 3, 4}`, `Matrix::det_direct()` evaluates the Leibniz expansion
// of the determinant as a tree of f64 multiplies and fused multiply-adds
// (FMAs).  Following Shewchuk's error-analysis methodology (REFERENCES.md
// [8]), the absolute error of that computation is bounded by
//
//     |det_direct(A) - det_exact(A)|  ≤  ERR_COEFF_D · p(|A|)
//
// where `p(|A|)` is the **absolute Leibniz sum**
//
//     p(|A|) = Σ_σ ∏ᵢ |A[i, σ(i)]|,
//
// i.e. the same cofactor-expansion tree as `det_direct` but with each
// entry replaced by its magnitude.  Note that `p(|A|)` is *not* the
// combinatorial matrix permanent — the name "permanent" appears in the
// source for brevity and to match the cited literature.
//
// Each constant has the shape `a · EPS + b · EPS²`: the linear term bounds
// the first-order rounding and the quadratic term absorbs the interaction
// of errors in nested FMAs.  The coefficients `a` and `b` are conservative
// over-estimates derived from the longest dependency chain of `det_direct`
// at that dimension.
//
// These constants are NOT feature-gated — they rely only on f64 arithmetic
// and are useful for adaptive-precision logic even without the `exact`
// feature.  Most callers should prefer `Matrix::det_errbound()`, which
// applies these constants to the actual matrix; the raw constants are
// exposed for advanced use cases (composing the bound with a pre-reduced
// permanent, rolling a custom adaptive filter, etc.).  See
// `Matrix::det_sign_exact()` (behind the `exact` feature) for the
// reference adaptive-filter that consumes these internally.
// ---------------------------------------------------------------------------

const EPS: f64 = f64::EPSILON; // 2^-52

/// Absolute error coefficient for [`Matrix::<2>::det_direct`](crate::Matrix::det_direct).
///
/// This constant is not a caller-tuned tolerance. It is the dimension-specific
/// multiplier that turns the matrix's absolute Leibniz sum into a conservative
/// bound on floating-point roundoff in the closed-form 2×2 determinant formula.
///
/// For any 2×2 matrix `A = [[a, b], [c, d]]` with finite f64 entries,
///
/// ```text
/// |A.det_direct() - det_exact(A)|  ≤  ERR_COEFF_2 · (|a·d| + |b·c|)
/// ```
///
/// `det_direct` evaluates `a·d - b·c` as one multiply followed by one FMA
/// (2 rounding events); the linear `3·EPS` term bounds those roundings
/// and the quadratic `16·EPS²` term is a conservative cushion for their
/// interaction.  Derivation follows Shewchuk's framework; see
/// `REFERENCES.md` \[8\].
///
/// Prefer [`Matrix::det_errbound`](crate::Matrix::det_errbound) unless
/// you already have the absolute-Leibniz sum available; see
/// `Matrix::det_sign_exact` (under the `exact` feature) for the reference
/// adaptive-precision filter.
///
/// # Example
/// ```
/// use la_stack::prelude::*;
///
/// # fn main() -> Result<(), LaError> {
/// let m = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]])?;
/// let Some(det) = m.det_direct()? else {
///     return Ok(());
/// };
/// assert_eq!(det, -2.0);
/// // Compute the bound from the raw constant for illustration; most
/// // callers would match on `m.det_errbound()?` instead.
/// let p = (1.0_f64 * 4.0).abs() + (2.0_f64 * 3.0).abs();
/// let bound = ERR_COEFF_2 * p;
/// if det.abs() > bound {
///     // The f64 sign is provably correct without exact arithmetic.
/// }
/// # Ok(())
/// # }
/// ```
pub const ERR_COEFF_2: f64 = 3.0 * EPS + 16.0 * EPS * EPS;

/// Absolute error coefficient for [`Matrix::<3>::det_direct`](crate::Matrix::det_direct).
///
/// This constant is not a caller-tuned tolerance. It is the dimension-specific
/// multiplier that turns the matrix's absolute Leibniz sum into a conservative
/// bound on floating-point roundoff in the closed-form 3×3 determinant formula.
///
/// For any 3×3 matrix `A` with finite f64 entries,
///
/// ```text
/// |A.det_direct() - det_exact(A)|  ≤  ERR_COEFF_3 · p(|A|)
/// ```
///
/// where `p(|A|)` is the absolute Leibniz sum (the same cofactor
/// expansion as `det_direct` but with `|·|` at every leaf).
/// `det_direct` for D=3 uses three 2×2 FMA minors combined by a nested
/// FMA, yielding the `8·EPS + 64·EPS²` bound.  See `REFERENCES.md`
/// \[8\] for the Shewchuk framework these bounds follow.
///
/// Prefer [`Matrix::det_errbound`](crate::Matrix::det_errbound) over this
/// constant for typical use; see [`ERR_COEFF_2`] for a worked example.
pub const ERR_COEFF_3: f64 = 8.0 * EPS + 64.0 * EPS * EPS;

/// Absolute error coefficient for [`Matrix::<4>::det_direct`](crate::Matrix::det_direct).
///
/// This constant is not a caller-tuned tolerance. It is the dimension-specific
/// multiplier that turns the matrix's absolute Leibniz sum into a conservative
/// bound on floating-point roundoff in the closed-form 4×4 determinant formula.
///
/// For any 4×4 matrix `A` with finite f64 entries,
///
/// ```text
/// |A.det_direct() - det_exact(A)|  ≤  ERR_COEFF_4 · p(|A|)
/// ```
///
/// where `p(|A|)` is the absolute Leibniz sum.  `det_direct` for D=4
/// hoists six 2×2 minors, combines them into four 3×3 cofactors, then
/// reduces those with an FMA row combination, yielding the
/// `12·EPS + 128·EPS²` bound.  See `REFERENCES.md` \[8\] for the
/// Shewchuk framework these bounds follow.
///
/// Prefer [`Matrix::det_errbound`](crate::Matrix::det_errbound) over this
/// constant for typical use; see [`ERR_COEFF_2`] for a worked example.
pub const ERR_COEFF_4: f64 = 12.0 * EPS + 128.0 * EPS * EPS;

/// Largest dimension supported by [`try_with_stack_matrix!`].
///
/// The crate can represent `Matrix<D>` for any compile-time `D`, but runtime
/// dispatch must enumerate a finite set of concrete stack types.  Dimensions
/// `0..=7` cover downstream geometric predicate matrices while keeping the
/// dispatch surface explicit.
pub const MAX_STACK_MATRIX_DISPATCH_DIM: usize = 7;

pub use error::{LaError, UnrepresentableReason};
pub use ldlt::Ldlt;
pub use lu::Lu;
pub use matrix::Matrix;
pub(crate) use tolerance::LDLT_SYMMETRY_REL_TOL;
pub use tolerance::{DEFAULT_SINGULAR_TOL, Tolerance};
pub use vector::Vector;

/// Fallibly dispatch a runtime dimension to a concrete stack-allocated matrix.
///
/// The macro creates a zero matrix with type `Matrix<N>` for the selected
/// runtime dimension `N`, then evaluates the supplied closure body.  Supported
/// runtime dimensions run from `0` through [`MAX_STACK_MATRIX_DISPATCH_DIM`].
/// Unsupported dimensions return
/// `Err(LaError::UnsupportedDimension { requested, max })` converted with
/// `From<LaError>`, so downstream crates can use their own public error type.
///
/// # Errors
/// Returns [`LaError::UnsupportedDimension`] (converted through `From<LaError>`)
/// when the requested runtime dimension is greater than
/// [`MAX_STACK_MATRIX_DISPATCH_DIM`].  The closure body may return any other
/// error representable by its declared `Result` type.
///
/// # Examples
/// ```
/// use la_stack::prelude::*;
///
/// # fn main() -> Result<(), LaError> {
/// let requested = 2usize;
/// let det = try_with_stack_matrix!(requested, |mut m| -> Result<f64, LaError> {
///     m.set_checked(0, 0, 1.0)?;
///     m.set_checked(1, 1, 1.0)?;
///     m.det()
/// })?;
///
/// assert_eq!(det, 1.0);
/// # Ok(())
/// # }
/// ```
#[macro_export]
macro_rules! try_with_stack_matrix {
    ($dim:expr, |$matrix:ident| -> $ret:ty $body:block $(,)?) => {{
        let __la_stack_requested_dim: usize = $dim;
        match __la_stack_requested_dim {
            0 => $crate::try_with_stack_matrix!(@arm 0, $matrix, $ret, $body),
            1 => $crate::try_with_stack_matrix!(@arm 1, $matrix, $ret, $body),
            2 => $crate::try_with_stack_matrix!(@arm 2, $matrix, $ret, $body),
            3 => $crate::try_with_stack_matrix!(@arm 3, $matrix, $ret, $body),
            4 => $crate::try_with_stack_matrix!(@arm 4, $matrix, $ret, $body),
            5 => $crate::try_with_stack_matrix!(@arm 5, $matrix, $ret, $body),
            6 => $crate::try_with_stack_matrix!(@arm 6, $matrix, $ret, $body),
            7 => $crate::try_with_stack_matrix!(@arm 7, $matrix, $ret, $body),
            requested => Err(::core::convert::From::from(
                $crate::LaError::unsupported_dimension(
                    requested,
                    $crate::MAX_STACK_MATRIX_DISPATCH_DIM,
                ),
            )),
        }
    }};
    ($dim:expr, |mut $matrix:ident| -> $ret:ty $body:block $(,)?) => {{
        let __la_stack_requested_dim: usize = $dim;
        match __la_stack_requested_dim {
            0 => $crate::try_with_stack_matrix!(@arm_mut 0, $matrix, $ret, $body),
            1 => $crate::try_with_stack_matrix!(@arm_mut 1, $matrix, $ret, $body),
            2 => $crate::try_with_stack_matrix!(@arm_mut 2, $matrix, $ret, $body),
            3 => $crate::try_with_stack_matrix!(@arm_mut 3, $matrix, $ret, $body),
            4 => $crate::try_with_stack_matrix!(@arm_mut 4, $matrix, $ret, $body),
            5 => $crate::try_with_stack_matrix!(@arm_mut 5, $matrix, $ret, $body),
            6 => $crate::try_with_stack_matrix!(@arm_mut 6, $matrix, $ret, $body),
            7 => $crate::try_with_stack_matrix!(@arm_mut 7, $matrix, $ret, $body),
            requested => Err(::core::convert::From::from(
                $crate::LaError::unsupported_dimension(
                    requested,
                    $crate::MAX_STACK_MATRIX_DISPATCH_DIM,
                ),
            )),
        }
    }};
    (@arm $d:literal, $matrix:ident, $ret:ty, $body:block) => {{
        let __la_stack_body = |$matrix: $crate::Matrix<$d>| -> $ret { $body };
        __la_stack_body($crate::Matrix::<$d>::zero())
    }};
    (@arm_mut $d:literal, $matrix:ident, $ret:ty, $body:block) => {{
        let __la_stack_body = |mut $matrix: $crate::Matrix<$d>| -> $ret { $body };
        __la_stack_body($crate::Matrix::<$d>::zero())
    }};
}

/// Common imports for ergonomic usage.
///
/// This prelude re-exports the primary types and constants: [`Matrix`],
/// [`Vector`], [`Lu`], [`Ldlt`], [`Tolerance`], [`LaError`],
/// [`UnrepresentableReason`], [`DEFAULT_SINGULAR_TOL`], and the determinant
/// error bound coefficients [`ERR_COEFF_2`], [`ERR_COEFF_3`], and
/// [`ERR_COEFF_4`]. It also re-exports [`MAX_STACK_MATRIX_DISPATCH_DIM`] and
/// [`try_with_stack_matrix!`] for runtime-to-const matrix dispatch.
///
/// When the `exact` feature is enabled, `BigInt` and `BigRational` are also
/// re-exported so callers can construct exact values (e.g. as the expected
/// result of `Matrix::det_exact`) without adding `num-bigint` / `num-rational`
/// to their own dependencies. The most commonly needed `num-traits` items are
/// re-exported alongside them: `FromPrimitive` for `BigRational::from_f64` /
/// `from_i64`, `ToPrimitive` for `BigRational::to_f64` / `to_i64`, and `Signed`
/// for `.is_positive()` / `.is_negative()` / `.abs()`.
pub mod prelude {
    pub use crate::{
        DEFAULT_SINGULAR_TOL, ERR_COEFF_2, ERR_COEFF_3, ERR_COEFF_4, LaError, Ldlt, Lu,
        MAX_STACK_MATRIX_DISPATCH_DIM, Matrix, Tolerance, UnrepresentableReason, Vector,
        try_with_stack_matrix,
    };

    #[cfg(feature = "exact")]
    pub use crate::{BigInt, BigRational, FromPrimitive, Signed, ToPrimitive};
}

#[cfg(test)]
mod tests {
    use super::*;

    use approx::assert_abs_diff_eq;

    mod prelude_tests {
        use approx::assert_abs_diff_eq;

        use crate::prelude::*;

        #[test]
        fn prelude_reexports_compile_and_work() -> Result<(), LaError> {
            // Use the items so we know they are in scope and usable.
            let m = Matrix::<2>::identity();
            let v = Vector::<2>::try_new([1.0, 2.0])?;
            let tol = Tolerance::new(0.0)?;
            assert_abs_diff_eq!(tol.get(), 0.0, epsilon = 0.0);
            assert_abs_diff_eq!(m.inf_norm()?, 1.0, epsilon = 0.0);
            assert_abs_diff_eq!(v.norm2_sq()?, 5.0, epsilon = 0.0);
            let _ = m.lu(DEFAULT_SINGULAR_TOL)?.solve(v)?;
            let _ = m.ldlt(DEFAULT_SINGULAR_TOL)?.solve(v)?;
            assert_eq!(
                LaError::unrepresentable(None, UnrepresentableReason::RequiresRounding),
                LaError::Unrepresentable {
                    index: None,
                    reason: UnrepresentableReason::RequiresRounding,
                }
            );
            assert_eq!(MAX_STACK_MATRIX_DISPATCH_DIM, 7);
            Ok(())
        }
    }

    macro_rules! gen_stack_matrix_dispatch_tests {
        ($d:literal) => {
            pastey::paste! {
                #[test]
                fn [<try_with_stack_matrix_dispatches_ $d d>]() {
                    let requested = $d;
                    let got = try_with_stack_matrix!(requested, |mut m| -> Result<usize, LaError> {
                        if $d > 0 {
                            m.set_checked($d - 1, $d - 1, f64::from($d))?;
                            assert_abs_diff_eq!(
                                m.get_checked($d - 1, $d - 1)?,
                                f64::from($d),
                                epsilon = 0.0
                            );
                        }
                        Ok($d)
                    });

                    assert_eq!(got, Ok($d));
                }
            }
        };
    }

    gen_stack_matrix_dispatch_tests!(2);
    gen_stack_matrix_dispatch_tests!(3);
    gen_stack_matrix_dispatch_tests!(4);
    gen_stack_matrix_dispatch_tests!(5);
    gen_stack_matrix_dispatch_tests!(6);
    gen_stack_matrix_dispatch_tests!(7);

    #[test]
    fn try_with_stack_matrix_supports_zero_dimension() {
        let got = try_with_stack_matrix!(0usize, |m| -> Result<Option<f64>, LaError> {
            m.det_direct()
        });

        assert_eq!(got, Ok(Some(1.0)));
    }

    #[test]
    fn try_with_stack_matrix_reports_unsupported_dimension() {
        let got = try_with_stack_matrix!(8usize, |m| -> Result<f64, LaError> { m.det() });

        assert_eq!(
            got,
            Err(LaError::UnsupportedDimension {
                requested: 8,
                max: MAX_STACK_MATRIX_DISPATCH_DIM,
            })
        );
    }

    #[derive(Debug, PartialEq)]
    struct DownstreamError(LaError);

    impl From<LaError> for DownstreamError {
        fn from(err: LaError) -> Self {
            Self(err)
        }
    }

    #[test]
    fn try_with_stack_matrix_converts_unsupported_dimension_error() {
        let got = try_with_stack_matrix!(9usize, |m| -> Result<usize, DownstreamError> {
            assert_abs_diff_eq!(m.inf_norm()?, 0.0, epsilon = 0.0);
            Ok(0)
        });

        assert_eq!(
            got,
            Err(DownstreamError(LaError::UnsupportedDimension {
                requested: 9,
                max: MAX_STACK_MATRIX_DISPATCH_DIM,
            }))
        );
    }

    /// Exercise every exact-feature re-export via the prelude so a future
    /// refactor that drops one (e.g. removing `Signed` from the prelude
    /// list) fails to compile rather than silently breaking downstream.
    #[cfg(feature = "exact")]
    #[test]
    fn prelude_exact_reexports_compile_and_work() {
        use crate::prelude::*;

        // `BigInt` and `BigRational` constructors.
        let n = BigInt::from(7);
        let r = BigRational::from_integer(n.clone());
        assert_eq!(*r.numer(), n);

        // `FromPrimitive::from_f64` / `from_i64` on `BigRational`.
        let half = BigRational::new(BigInt::from(1), BigInt::from(2));
        let two = BigRational::from_integer(BigInt::from(2));
        assert_eq!(BigRational::from_f64(0.5), Some(half.clone()));
        assert_eq!(BigRational::from_i64(2), Some(two.clone()));
        assert_eq!(
            half.clone() + half.clone(),
            BigRational::from_integer(BigInt::from(1))
        );

        // `Signed::is_positive` / `is_negative` / `abs`.
        assert!(half.is_positive());
        assert!(!half.is_negative());
        let neg = -half.clone();
        assert!(neg.is_negative());
        assert_eq!(neg.abs(), half);

        // `ToPrimitive::to_f64` / `to_i64`.
        assert_eq!(half.to_f64(), Some(0.5));
        assert_eq!(two.to_i64(), Some(2));
    }
}
