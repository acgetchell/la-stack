#![forbid(unsafe_code)]
#![deny(missing_docs)]
#![cfg_attr(docsrs, feature(doc_cfg))]
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

    /// ```rust
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let x = Interval::try_from_subtraction(0.1, 0.0)?;
    /// let y = Interval::try_from_subtraction(0.1, 0.0)?;
    /// let z = Interval::try_from_subtraction(0.1, 0.0)?;
    /// let lifted = x
    ///     .try_square()?
    ///     .try_add(&y.try_square()?)?
    ///     .try_add(&z.try_square()?)?;
    ///
    /// let matrix = IntervalMatrix::<4>::from_rows([
    ///     [Interval::ONE, Interval::ZERO, Interval::ZERO, Interval::ONE],
    ///     [Interval::ZERO, Interval::ONE, Interval::ZERO, Interval::ONE],
    ///     [Interval::ZERO, Interval::ZERO, Interval::ONE, Interval::ONE],
    ///     [x, y, z, lifted],
    /// ]);
    /// assert_eq!(
    ///     matrix.det_sign()?,
    ///     IntervalDeterminantSign::Negative,
    /// );
    /// # Ok(())
    /// # }
    /// ```
    fn interval_determinant_example() {}

    /// ```rust
    /// use la_stack::prelude::*;
    ///
    /// fn is_separated<const D: usize>(
    ///     axis: &Vector<D>,
    ///     left: &Vector<D>,
    ///     right: &Vector<D>,
    ///     threshold: f64,
    /// ) -> Result<Option<bool>, LaError> {
    ///     let Some(value) = axis.dot_difference_with_errbound(left, right)? else {
    ///         return Ok(None);
    ///     };
    ///     if value.lower_bound() > threshold {
    ///         Ok(Some(true))
    ///     } else if value.upper_bound() <= threshold {
    ///         Ok(Some(false))
    ///     } else {
    ///         Ok(None)
    ///     }
    /// }
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let axis = Vector::<2>::try_new([2.0, -1.0])?;
    /// let left = Vector::<2>::try_new([4.0, 1.0])?;
    /// let right = Vector::<2>::try_new([1.0, 3.0])?;
    /// assert_eq!(is_separated(&axis, &left, &right, 1.0)?, Some(true));
    /// # Ok(())
    /// # }
    /// ```
    fn certified_linear_form_example() {}

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
    /// assert_eq!(m.det_sign_exact(), DeterminantSign::Zero); // exactly singular
    ///
    /// let det = m.det_exact()?;
    /// assert_eq!(det, BigRational::from_integer(0.into())); // exact zero
    /// let det_f64 = det.try_to_f64()?;
    /// assert_eq!(det_f64, 0.0);
    ///
    /// // If strict exact-to-f64 conversion would require rounding, opt in
    /// // explicitly with the rounded API.
    /// let inexact = Matrix::<2>::try_from_rows([
    ///     [1.0 + f64::EPSILON, 0.0],
    ///     [0.0, 1.0 - f64::EPSILON],
    /// ])?;
    /// let exact_det = inexact.det_exact()?;
    /// let rounded_det = match exact_det.try_to_f64() {
    ///     Ok(det) => det,
    ///     Err(err) if err.requires_rounding() => exact_det.to_rounded_f64()?,
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
    ///     huge_det
    ///         .try_to_f64()
    ///         .err()
    ///         .and_then(|err| err.unrepresentable_reason()),
    ///     Some(UnrepresentableReason::NotFinite)
    /// );
    /// println!("exact determinant = {huge_det}");
    ///
    /// // Exact linear system solve
    /// let a = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]])?;
    /// let b = Vector::<2>::try_new([5.0, 11.0])?;
    /// let exact_x = a.solve_exact(b)?;
    /// let x = exact_x.try_to_f64()?.into_array();
    /// assert!((x[0] - 1.0).abs() <= f64::EPSILON);
    /// assert!((x[1] - 2.0).abs() <= f64::EPSILON);
    /// # Ok(())
    /// # }
    /// ```
    fn exact_arithmetic_example() {}

    #[cfg(feature = "exact")]
    /// ```rust
    /// use core::assert_matches;
    ///
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let epsilon = BigRational::new(1.into(), (1_u64 << 60).into());
    /// let one = BigRational::from_integer(1.into());
    /// let zero = BigRational::from_integer(0.into());
    ///
    /// let matrix = RationalMatrix::<5>::try_from_fn(|row, col| match (row, col) {
    ///     (0, 0 | 1) | (1, 0) => one.clone(),
    ///     (1, 1) => &one + &epsilon,
    ///     _ if row == col => one.clone(),
    ///     _ => zero.clone(),
    /// })?;
    /// assert_eq!(matrix.det_sign(), DeterminantSign::Positive);
    /// assert_eq!(matrix.det(), epsilon);
    ///
    /// let rhs = RationalVector::try_new([
    ///     zero,
    ///     -&epsilon,
    ///     BigRational::from_integer(2.into()),
    ///     BigRational::from_integer(3.into()),
    ///     BigRational::from_integer(4.into()),
    /// ])?;
    /// let exact_solution = matrix.solve(&rhs)?;
    /// assert_eq!(
    ///     exact_solution.as_array(),
    ///     &[
    ///         BigRational::from_integer(1.into()),
    ///         BigRational::from_integer((-1).into()),
    ///         BigRational::from_integer(2.into()),
    ///         BigRational::from_integer(3.into()),
    ///         BigRational::from_integer(4.into()),
    ///     ]
    /// );
    ///
    /// let epsilon_f64 = epsilon.try_to_f64()?;
    /// assert_eq!((1.0 + epsilon_f64).to_bits(), 1.0_f64.to_bits());
    /// let f64_matrix = Matrix::<5>::try_from_rows([
    ///     [1.0, 1.0, 0.0, 0.0, 0.0],
    ///     [1.0, 1.0 + epsilon_f64, 0.0, 0.0, 0.0],
    ///     [0.0, 0.0, 1.0, 0.0, 0.0],
    ///     [0.0, 0.0, 0.0, 1.0, 0.0],
    ///     [0.0, 0.0, 0.0, 0.0, 1.0],
    /// ])?;
    /// let f64_rhs = Vector::<5>::try_new([0.0, -epsilon_f64, 2.0, 3.0, 4.0])?;
    /// let f64_solve = f64_matrix
    ///     .lu(DEFAULT_SINGULAR_TOL)
    ///     .and_then(|lu| lu.solve(f64_rhs));
    /// assert_matches!(
    ///     f64_solve,
    ///     Err(LaError::Singular { .. })
    /// );
    /// # Ok(())
    /// # }
    /// ```
    fn rational_input_example() {}

    #[cfg(feature = "exact")]
    /// ```rust
    /// use la_stack::prelude::*;
    ///
    /// fn adaptive_det_sign<const D: usize>(
    ///     matrix: &Matrix<D>,
    /// ) -> DeterminantSign {
    ///     if let Ok(Some(estimate)) = matrix.det_direct_with_errbound() {
    ///         if estimate.determinant().abs() > estimate.absolute_error_bound() {
    ///             return if estimate.determinant() > 0.0 {
    ///                 DeterminantSign::Positive
    ///             } else {
    ///                 DeterminantSign::Negative
    ///             };
    ///         }
    ///     }
    ///
    ///     matrix.det_sign_exact()
    /// }
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let identity = Matrix::<3>::identity();
    /// assert_eq!(
    ///     adaptive_det_sign(&identity),
    ///     DeterminantSign::Positive
    /// );
    ///
    /// let singular = Matrix::<3>::try_from_rows([
    ///     [1.0, 2.0, 3.0],
    ///     [4.0, 5.0, 6.0],
    ///     [7.0, 8.0, 9.0],
    /// ])?;
    /// assert_eq!(adaptive_det_sign(&singular), DeterminantSign::Zero);
    ///
    /// let big = f64::MAX / 2.0;
    /// let overflowing = Matrix::<3>::try_from_rows([
    ///     [0.0, 0.0, 1.0],
    ///     [big, 0.0, 1.0],
    ///     [0.0, big, 1.0],
    /// ])?;
    /// assert_eq!(
    ///     adaptive_det_sign(&overflowing),
    ///     DeterminantSign::Positive
    /// );
    /// # Ok(())
    /// # }
    /// ```
    fn adaptive_precision_example() {}
}

mod error;
#[cfg(feature = "exact")]
mod exact;
mod interval;
mod ldlt;
mod lu;
mod matrix;
#[cfg(feature = "exact")]
mod rational;
mod scaled_product;
mod tolerance;
mod vector;

#[cfg(feature = "exact")]
#[cfg_attr(docsrs, doc(cfg(feature = "exact")))]
pub use exact::{DeterminantSign, ExactF64Conversion};
#[cfg(feature = "exact")]
#[cfg_attr(docsrs, doc(cfg(feature = "exact")))]
pub use num_bigint::BigInt;
#[cfg(feature = "exact")]
#[cfg_attr(docsrs, doc(cfg(feature = "exact")))]
pub use num_rational::BigRational;
#[cfg(feature = "exact")]
#[cfg_attr(docsrs, doc(cfg(feature = "exact")))]
pub use num_traits::{FromPrimitive, Signed, ToPrimitive};
#[cfg(feature = "exact")]
#[cfg_attr(docsrs, doc(cfg(feature = "exact")))]
pub use rational::{RationalMatrix, RationalVector};

// ---------------------------------------------------------------------------
// Error-bound constants for `Matrix::det_direct_with_errbound()` and
// `Matrix::det_errbound()`.
//
// For `D ∈ {2, 3, 4}`, `Matrix::det_direct()` evaluates the Leibniz expansion
// of the determinant as a tree of f64 multiplies and fused multiply-adds
// (FMAs).  When every rounded intermediate is normal or an exact structural
// zero, Shewchuk's error-analysis methodology (REFERENCES.md [8]) bounds the
// absolute error of that computation by
//
//     |det_direct(A) - det_exact(A)|  ≤  ERR_COEFF_D · p(|A|)
//
// where `p(|A|)` is the **absolute Leibniz sum**
//
//     p(|A|) = Σ_σ ∏ᵢ |A[i, σ(i)]|,
//
// i.e. exactly the combinatorial matrix permanent `perm(|A|)`. The
// implementation evaluates the corresponding fixed-size expansion in f64, so
// the computed `permanent` value used by the bound may itself be rounded even
// though the mathematical quantity above is exact.
//
// Each constant has the shape `a · EPS + b · EPS²`: the linear term bounds
// the first-order rounding and the quadratic term absorbs the interaction
// of errors in nested FMAs.  The coefficients `a` and `b` are conservative
// over-estimates derived from the longest dependency chain of `det_direct`
// at that dimension.
//
// These constants are NOT feature-gated — they rely only on f64 arithmetic
// and are useful for adaptive-precision logic even without the `exact`
// feature. Most callers should prefer `Matrix::det_direct_with_errbound()`
// when they need the approximation and bound together, or
// `Matrix::det_errbound()` when they need only the bound. Those methods apply
// these constants to the actual matrix; the raw constants are
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
/// For a 2×2 matrix `A = [[a, b], [c, d]]` whose closed-form determinant
/// intermediates do not undergo gradual underflow,
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
/// Prefer
/// [`Matrix::det_direct_with_errbound`](crate::Matrix::det_direct_with_errbound)
/// unless you need only the bound or already have the absolute-Leibniz sum;
/// see
/// `Matrix::det_sign_exact` (under the `exact` feature) for the reference
/// adaptive-precision filter.
///
/// # Example
/// ```
/// use la_stack::{prelude::*, ERR_COEFF_2};
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
/// For a 3×3 matrix `A` whose closed-form determinant intermediates do not
/// undergo gradual underflow,
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
/// Prefer
/// [`Matrix::det_direct_with_errbound`](crate::Matrix::det_direct_with_errbound)
/// over this constant for typical use; see [`ERR_COEFF_2`] for a worked
/// example.
pub const ERR_COEFF_3: f64 = 8.0 * EPS + 64.0 * EPS * EPS;

/// Absolute error coefficient for [`Matrix::<4>::det_direct`](crate::Matrix::det_direct).
///
/// This constant is not a caller-tuned tolerance. It is the dimension-specific
/// multiplier that turns the matrix's absolute Leibniz sum into a conservative
/// bound on floating-point roundoff in the closed-form 4×4 determinant formula.
///
/// For a 4×4 matrix `A` whose closed-form determinant intermediates do not
/// undergo gradual underflow,
///
/// ```text
/// |A.det_direct() - det_exact(A)|  ≤  ERR_COEFF_4 · p(|A|)
/// ```
///
/// where `p(|A|)` is the absolute Leibniz sum. `det_direct` for D=4
/// evaluates four nested 3×3 cofactors, sharing their six 2×2 minors when
/// every coefficient in the first two rows is non-zero, and reduces them with
/// an FMA row combination, yielding the
/// `12·EPS + 128·EPS²` bound.  See `REFERENCES.md` \[8\] for the
/// Shewchuk framework these bounds follow.
///
/// Prefer
/// [`Matrix::det_direct_with_errbound`](crate::Matrix::det_direct_with_errbound)
/// over this constant for typical use; see [`ERR_COEFF_2`] for a worked
/// example.
pub const ERR_COEFF_4: f64 = 12.0 * EPS + 128.0 * EPS * EPS;

/// Largest dimension supported by [`try_with_stack_matrix!`].
///
/// The crate can represent `Matrix<D>` for any compile-time `D`, but runtime
/// dispatch must enumerate a finite set of concrete stack types.  Dimensions
/// `0..=7` cover downstream geometric predicate matrices while keeping the
/// dispatch surface explicit.
pub const MAX_STACK_MATRIX_DISPATCH_DIM: usize = 7;

/// Largest dimension supported by [`try_with_rational_matrix!`].
///
/// This bound covers exact downstream geometric systems through dimension 8
/// while keeping runtime-to-const dispatch explicit and available on stable
/// Rust.
#[cfg(feature = "exact")]
#[cfg_attr(docsrs, doc(cfg(feature = "exact")))]
pub const MAX_RATIONAL_MATRIX_DISPATCH_DIM: usize = 8;

pub use error::{
    ArithmeticOperation, FactorizationKind, IntervalBound, IntervalOperand, InvalidToleranceReason,
    LaError, NonFiniteLocation, NonFiniteOrigin, PositiveSemidefiniteViolation, SingularityReason,
    UnrepresentableReason,
};
pub use interval::{Interval, IntervalDeterminantSign, IntervalMatrix, MAX_INTERVAL_MATRIX_DIM};
pub use ldlt::Ldlt;
pub use lu::Lu;
pub use matrix::{DeterminantWithErrorBound, Matrix};
pub use tolerance::{DEFAULT_SINGULAR_TOL, Tolerance};
pub use vector::{ScalarWithErrorBound, Vector};

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
///     m.set(0, 0, 1.0)?;
///     m.set(1, 1, 1.0)?;
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

/// Fallibly dispatch a runtime dimension to a concrete interval matrix.
///
/// The macro creates a zero [`IntervalMatrix`] with the selected const-generic
/// dimension, then evaluates the closure body. Supported dimensions run from
/// `0` through [`MAX_INTERVAL_MATRIX_DIM`]. Unsupported dimensions return
/// [`LaError::UnsupportedDimension`] converted through `From<LaError>`.
///
/// # Errors
/// Returns [`LaError::UnsupportedDimension`] (converted through
/// `From<LaError>`) when the requested dimension is greater than
/// [`MAX_INTERVAL_MATRIX_DIM`]. The closure body may return any other error
/// representable by its declared `Result` type.
///
/// # Examples
/// ```
/// use la_stack::prelude::*;
///
/// # fn main() -> Result<(), LaError> {
/// let requested = 3usize;
/// let sign = try_with_interval_matrix!(requested, |mut matrix| -> Result<
///     IntervalDeterminantSign,
///     LaError,
/// > {
///     for index in 0..requested {
///         matrix.set(index, index, Interval::ONE)?;
///     }
///     matrix.det_sign()
/// })?;
/// assert_eq!(sign, IntervalDeterminantSign::Positive);
/// # Ok(())
/// # }
/// ```
#[macro_export]
macro_rules! try_with_interval_matrix {
    ($dim:expr, |$matrix:ident| -> $ret:ty $body:block $(,)?) => {{
        let __la_stack_requested_dim: usize = $dim;
        match __la_stack_requested_dim {
            0 => $crate::try_with_interval_matrix!(@arm 0, $matrix, $ret, $body),
            1 => $crate::try_with_interval_matrix!(@arm 1, $matrix, $ret, $body),
            2 => $crate::try_with_interval_matrix!(@arm 2, $matrix, $ret, $body),
            3 => $crate::try_with_interval_matrix!(@arm 3, $matrix, $ret, $body),
            4 => $crate::try_with_interval_matrix!(@arm 4, $matrix, $ret, $body),
            5 => $crate::try_with_interval_matrix!(@arm 5, $matrix, $ret, $body),
            6 => $crate::try_with_interval_matrix!(@arm 6, $matrix, $ret, $body),
            7 => $crate::try_with_interval_matrix!(@arm 7, $matrix, $ret, $body),
            requested => Err(::core::convert::From::from(
                $crate::LaError::unsupported_dimension(
                    requested,
                    $crate::MAX_INTERVAL_MATRIX_DIM,
                ),
            )),
        }
    }};
    ($dim:expr, |mut $matrix:ident| -> $ret:ty $body:block $(,)?) => {{
        let __la_stack_requested_dim: usize = $dim;
        match __la_stack_requested_dim {
            0 => $crate::try_with_interval_matrix!(@arm_mut 0, $matrix, $ret, $body),
            1 => $crate::try_with_interval_matrix!(@arm_mut 1, $matrix, $ret, $body),
            2 => $crate::try_with_interval_matrix!(@arm_mut 2, $matrix, $ret, $body),
            3 => $crate::try_with_interval_matrix!(@arm_mut 3, $matrix, $ret, $body),
            4 => $crate::try_with_interval_matrix!(@arm_mut 4, $matrix, $ret, $body),
            5 => $crate::try_with_interval_matrix!(@arm_mut 5, $matrix, $ret, $body),
            6 => $crate::try_with_interval_matrix!(@arm_mut 6, $matrix, $ret, $body),
            7 => $crate::try_with_interval_matrix!(@arm_mut 7, $matrix, $ret, $body),
            requested => Err(::core::convert::From::from(
                $crate::LaError::unsupported_dimension(
                    requested,
                    $crate::MAX_INTERVAL_MATRIX_DIM,
                ),
            )),
        }
    }};
    (@arm $d:literal, $matrix:ident, $ret:ty, $body:block) => {{
        let __la_stack_body = |$matrix: $crate::IntervalMatrix<$d>| -> $ret { $body };
        __la_stack_body($crate::IntervalMatrix::<$d>::zero())
    }};
    (@arm_mut $d:literal, $matrix:ident, $ret:ty, $body:block) => {{
        let __la_stack_body = |mut $matrix: $crate::IntervalMatrix<$d>| -> $ret { $body };
        __la_stack_body($crate::IntervalMatrix::<$d>::zero())
    }};
}

/// Fallibly dispatch a runtime dimension to a concrete exact rational matrix.
///
/// The macro creates a zero [`RationalMatrix`] with the selected const-generic
/// dimension, then evaluates the supplied closure body. Dimensions `0..=8` are
/// supported on stable Rust. The closure may fill the matrix through
/// [`RationalMatrix::set`] or replace it with a value built by
/// [`RationalMatrix::try_from_fn`].
///
/// # Errors
/// Returns [`LaError::UnsupportedDimension`] (converted through
/// `From<LaError>`) when the requested dimension is greater than
/// [`MAX_RATIONAL_MATRIX_DISPATCH_DIM`]. The closure body may return any other
/// error representable by its declared `Result` type.
///
/// # Examples
/// ```
/// use la_stack::prelude::*;
///
/// # fn main() -> Result<(), LaError> {
/// let requested = 3usize;
/// let sign = try_with_rational_matrix!(requested, |mut matrix| -> Result<
///     DeterminantSign,
///     LaError,
/// > {
///     for index in 0..requested {
///         matrix.set(index, index, BigRational::from_integer(1.into()))?;
///     }
///     Ok(matrix.det_sign())
/// })?;
/// assert_eq!(sign, DeterminantSign::Positive);
/// # Ok(())
/// # }
/// ```
#[cfg(feature = "exact")]
#[cfg_attr(docsrs, doc(cfg(feature = "exact")))]
#[macro_export]
macro_rules! try_with_rational_matrix {
    ($dim:expr, |$matrix:ident| -> $ret:ty $body:block $(,)?) => {{
        let __la_stack_requested_dim: usize = $dim;
        match __la_stack_requested_dim {
            0 => $crate::try_with_rational_matrix!(@arm 0, $matrix, $ret, $body),
            1 => $crate::try_with_rational_matrix!(@arm 1, $matrix, $ret, $body),
            2 => $crate::try_with_rational_matrix!(@arm 2, $matrix, $ret, $body),
            3 => $crate::try_with_rational_matrix!(@arm 3, $matrix, $ret, $body),
            4 => $crate::try_with_rational_matrix!(@arm 4, $matrix, $ret, $body),
            5 => $crate::try_with_rational_matrix!(@arm 5, $matrix, $ret, $body),
            6 => $crate::try_with_rational_matrix!(@arm 6, $matrix, $ret, $body),
            7 => $crate::try_with_rational_matrix!(@arm 7, $matrix, $ret, $body),
            8 => $crate::try_with_rational_matrix!(@arm 8, $matrix, $ret, $body),
            requested => Err(::core::convert::From::from(
                $crate::LaError::unsupported_dimension(
                    requested,
                    $crate::MAX_RATIONAL_MATRIX_DISPATCH_DIM,
                ),
            )),
        }
    }};
    ($dim:expr, |mut $matrix:ident| -> $ret:ty $body:block $(,)?) => {{
        let __la_stack_requested_dim: usize = $dim;
        match __la_stack_requested_dim {
            0 => $crate::try_with_rational_matrix!(@arm_mut 0, $matrix, $ret, $body),
            1 => $crate::try_with_rational_matrix!(@arm_mut 1, $matrix, $ret, $body),
            2 => $crate::try_with_rational_matrix!(@arm_mut 2, $matrix, $ret, $body),
            3 => $crate::try_with_rational_matrix!(@arm_mut 3, $matrix, $ret, $body),
            4 => $crate::try_with_rational_matrix!(@arm_mut 4, $matrix, $ret, $body),
            5 => $crate::try_with_rational_matrix!(@arm_mut 5, $matrix, $ret, $body),
            6 => $crate::try_with_rational_matrix!(@arm_mut 6, $matrix, $ret, $body),
            7 => $crate::try_with_rational_matrix!(@arm_mut 7, $matrix, $ret, $body),
            8 => $crate::try_with_rational_matrix!(@arm_mut 8, $matrix, $ret, $body),
            requested => Err(::core::convert::From::from(
                $crate::LaError::unsupported_dimension(
                    requested,
                    $crate::MAX_RATIONAL_MATRIX_DISPATCH_DIM,
                ),
            )),
        }
    }};
    (@arm $d:literal, $matrix:ident, $ret:ty, $body:block) => {{
        let __la_stack_body = |$matrix: $crate::RationalMatrix<$d>| -> $ret { $body };
        __la_stack_body($crate::RationalMatrix::<$d>::zero())
    }};
    (@arm_mut $d:literal, $matrix:ident, $ret:ty, $body:block) => {{
        let __la_stack_body = |mut $matrix: $crate::RationalMatrix<$d>| -> $ret { $body };
        __la_stack_body($crate::RationalMatrix::<$d>::zero())
    }};
}

/// Common imports for ergonomic usage.
///
/// This prelude re-exports the primary types and common constants: [`Matrix`],
/// [`DeterminantWithErrorBound`], [`Interval`], [`IntervalMatrix`],
/// [`IntervalDeterminantSign`], [`ScalarWithErrorBound`], [`Vector`], [`Lu`],
/// [`Ldlt`], [`Tolerance`],
/// and [`LaError`]. Its typed
/// error categories include [`ArithmeticOperation`], [`FactorizationKind`],
/// [`IntervalBound`], [`IntervalOperand`], [`InvalidToleranceReason`],
/// [`NonFiniteLocation`], [`NonFiniteOrigin`], [`PositiveSemidefiniteViolation`],
/// [`SingularityReason`], and [`UnrepresentableReason`]. It also re-exports
/// [`DEFAULT_SINGULAR_TOL`],
/// [`MAX_STACK_MATRIX_DISPATCH_DIM`], [`MAX_INTERVAL_MATRIX_DIM`],
/// [`try_with_stack_matrix!`], and [`try_with_interval_matrix!`] for
/// runtime-to-const matrix dispatch. Advanced custom-filter code should import
/// [`ERR_COEFF_2`], [`ERR_COEFF_3`], and [`ERR_COEFF_4`] explicitly from the
/// crate root; those raw coefficients intentionally stay out of the prelude.
#[cfg_attr(feature = "exact", doc = "")]
#[cfg_attr(
    feature = "exact",
    doc = "When the `exact` feature is enabled, [`RationalMatrix`], [`RationalVector`],"
)]
#[cfg_attr(
    feature = "exact",
    doc = "[`DeterminantSign`], [`ExactF64Conversion`], [`BigInt`], and [`BigRational`]"
)]
#[cfg_attr(
    feature = "exact",
    doc = "are also re-exported, together with [`MAX_RATIONAL_MATRIX_DISPATCH_DIM`] and"
)]
#[cfg_attr(
    feature = "exact",
    doc = "[`try_with_rational_matrix!`] for runtime-to-const exact-matrix dispatch."
)]
#[cfg_attr(
    feature = "exact",
    doc = "[`ExactF64Conversion`] converts an already-computed exact determinant or solution"
)]
#[cfg_attr(
    feature = "exact",
    doc = "under either the strict or explicitly rounded binary64 contract, without repeating"
)]
#[cfg_attr(
    feature = "exact",
    doc = "exact elimination. The number types let callers construct expected exact values"
)]
#[cfg_attr(
    feature = "exact",
    doc = "without adding `num-bigint` / `num-rational` to their own dependencies. The most"
)]
#[cfg_attr(
    feature = "exact",
    doc = "commonly needed `num-traits` items are re-exported alongside them: [`FromPrimitive`]"
)]
#[cfg_attr(
    feature = "exact",
    doc = "for `BigRational::from_f64` / `from_i64`, [`ToPrimitive`] for"
)]
#[cfg_attr(
    feature = "exact",
    doc = "`BigRational::to_f64` / `to_i64`, and [`Signed`] for `is_positive` / `is_negative` /"
)]
#[cfg_attr(feature = "exact", doc = "`abs`.")]
pub mod prelude {
    pub use crate::{
        ArithmeticOperation, DEFAULT_SINGULAR_TOL, DeterminantWithErrorBound, FactorizationKind,
        Interval, IntervalBound, IntervalDeterminantSign, IntervalMatrix, IntervalOperand,
        InvalidToleranceReason, LaError, Ldlt, Lu, MAX_INTERVAL_MATRIX_DIM,
        MAX_STACK_MATRIX_DISPATCH_DIM, Matrix, NonFiniteLocation, NonFiniteOrigin,
        PositiveSemidefiniteViolation, ScalarWithErrorBound, SingularityReason, Tolerance,
        UnrepresentableReason, Vector, try_with_interval_matrix, try_with_stack_matrix,
    };

    #[cfg(feature = "exact")]
    #[cfg_attr(docsrs, doc(cfg(feature = "exact")))]
    pub use crate::{
        BigInt, BigRational, DeterminantSign, ExactF64Conversion, FromPrimitive,
        MAX_RATIONAL_MATRIX_DISPATCH_DIM, RationalMatrix, RationalVector, Signed, ToPrimitive,
        try_with_rational_matrix,
    };
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;
    use pastey::paste;

    use super::*;

    macro_rules! gen_stack_matrix_dispatch_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<try_with_stack_matrix_dispatches_ $d d>]() {
                    let requested = $d;
                    let got = try_with_stack_matrix!(requested, |mut m| -> Result<usize, LaError> {
                        if $d > 0 {
                            m.set($d - 1, $d - 1, f64::from($d))?;
                            assert_abs_diff_eq!(
                                m.try_get($d - 1, $d - 1)?,
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

    gen_stack_matrix_dispatch_tests!(1);
    gen_stack_matrix_dispatch_tests!(2);
    gen_stack_matrix_dispatch_tests!(3);
    gen_stack_matrix_dispatch_tests!(4);
    gen_stack_matrix_dispatch_tests!(5);
    gen_stack_matrix_dispatch_tests!(6);
    gen_stack_matrix_dispatch_tests!(7);

    macro_rules! gen_interval_matrix_dispatch_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<try_with_interval_matrix_dispatches_ $d d>]() {
                    let requested = $d;
                    let got = try_with_interval_matrix!(
                        requested,
                        |mut matrix| -> Result<IntervalDeterminantSign, LaError> {
                            let mut index = 0;
                            while index < $d {
                                matrix.set(index, index, Interval::ONE)?;
                                index += 1;
                            }
                            matrix.det_sign()
                        },
                    );

                    assert_eq!(got, Ok(IntervalDeterminantSign::Positive));
                }
            }
        };
    }

    gen_interval_matrix_dispatch_tests!(1);
    gen_interval_matrix_dispatch_tests!(2);
    gen_interval_matrix_dispatch_tests!(3);
    gen_interval_matrix_dispatch_tests!(4);
    gen_interval_matrix_dispatch_tests!(5);
    gen_interval_matrix_dispatch_tests!(6);
    gen_interval_matrix_dispatch_tests!(7);

    #[cfg(feature = "exact")]
    macro_rules! gen_rational_matrix_dispatch_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<try_with_rational_matrix_dispatches_ $d d>]() {
                    let requested = $d;
                    let got = try_with_rational_matrix!(
                        requested,
                        |mut matrix| -> Result<DeterminantSign, LaError> {
                            let mut index = 0;
                            while index < $d {
                                matrix.set(
                                    index,
                                    index,
                                    BigRational::from_integer(BigInt::from(1)),
                                )?;
                                index += 1;
                            }
                            Ok(matrix.det_sign())
                        },
                    );

                    assert_eq!(got, Ok(DeterminantSign::Positive));
                }
            }
        };
    }

    #[cfg(feature = "exact")]
    gen_rational_matrix_dispatch_tests!(1);
    #[cfg(feature = "exact")]
    gen_rational_matrix_dispatch_tests!(2);
    #[cfg(feature = "exact")]
    gen_rational_matrix_dispatch_tests!(3);
    #[cfg(feature = "exact")]
    gen_rational_matrix_dispatch_tests!(4);
    #[cfg(feature = "exact")]
    gen_rational_matrix_dispatch_tests!(5);
    #[cfg(feature = "exact")]
    gen_rational_matrix_dispatch_tests!(6);
    #[cfg(feature = "exact")]
    gen_rational_matrix_dispatch_tests!(7);
    #[cfg(feature = "exact")]
    gen_rational_matrix_dispatch_tests!(8);

    #[cfg(feature = "exact")]
    #[test]
    fn try_with_rational_matrix_dispatches_zero_dimension() {
        let got = try_with_rational_matrix!(0usize, |matrix| -> Result<DeterminantSign, LaError> {
            Ok(matrix.det_sign())
        });

        assert_eq!(got, Ok(DeterminantSign::Positive));
    }

    #[test]
    fn try_with_stack_matrix_supports_zero_dimension() {
        let got = try_with_stack_matrix!(0usize, |m| -> Result<Option<f64>, LaError> {
            m.det_direct()
        });

        assert_eq!(got, Ok(Some(1.0)));
    }

    #[test]
    fn try_with_interval_matrix_supports_zero_dimension() {
        let got = try_with_interval_matrix!(0usize, |matrix| -> Result<
            IntervalDeterminantSign,
            LaError,
        > { matrix.det_sign() },);

        assert_eq!(got, Ok(IntervalDeterminantSign::Positive));
    }

    #[test]
    fn try_with_stack_matrix_evaluates_dimension_once() {
        let mut evaluations = 0;
        let got = try_with_stack_matrix!(
            {
                evaluations += 1;
                2usize
            },
            |matrix| -> Result<f64, LaError> { matrix.try_get(1, 1) },
        );

        assert_eq!(evaluations, 1);
        assert_eq!(got, Ok(0.0));
    }

    #[test]
    fn try_with_interval_matrix_evaluates_dimension_once() {
        let mut evaluations = 0;
        let got = try_with_interval_matrix!(
            {
                evaluations += 1;
                2usize
            },
            |matrix| -> Result<Interval, LaError> { matrix.try_get(1, 1) },
        );

        assert_eq!(evaluations, 1);
        assert_eq!(got, Ok(Interval::ZERO));
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

    #[test]
    fn try_with_interval_matrix_reports_unsupported_dimension() {
        let got = try_with_interval_matrix!(8usize, |matrix| -> Result<
            IntervalDeterminantSign,
            LaError,
        > { matrix.det_sign() },);

        assert_eq!(
            got,
            Err(LaError::UnsupportedDimension {
                requested: 8,
                max: MAX_INTERVAL_MATRIX_DIM,
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

    #[test]
    fn try_with_interval_matrix_converts_unsupported_dimension_error() {
        let got = try_with_interval_matrix!(8usize, |matrix| -> Result<
            IntervalDeterminantSign,
            DownstreamError,
        > { Ok(matrix.det_sign()?) },);

        assert_eq!(
            got,
            Err(DownstreamError(LaError::UnsupportedDimension {
                requested: 8,
                max: MAX_INTERVAL_MATRIX_DIM,
            }))
        );
    }

    #[cfg(feature = "exact")]
    #[test]
    fn try_with_rational_matrix_reports_unsupported_dimension() {
        let got = try_with_rational_matrix!(9usize, |matrix| -> Result<BigRational, LaError> {
            Ok(matrix.det())
        });

        assert_eq!(
            got,
            Err(LaError::UnsupportedDimension {
                requested: 9,
                max: MAX_RATIONAL_MATRIX_DISPATCH_DIM,
            })
        );
    }

    #[cfg(feature = "exact")]
    #[test]
    fn try_with_rational_matrix_converts_unsupported_dimension_error() {
        let got =
            try_with_rational_matrix!(9usize, |matrix| -> Result<BigRational, DownstreamError> {
                Ok(matrix.det())
            });

        assert_eq!(
            got,
            Err(DownstreamError(LaError::UnsupportedDimension {
                requested: 9,
                max: MAX_RATIONAL_MATRIX_DISPATCH_DIM,
            }))
        );
    }
}
