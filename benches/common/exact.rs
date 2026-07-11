#![forbid(unsafe_code)]

//! Shared helpers for exact-arithmetic benchmark input generation and tests.

use core::array::{self, from_fn};
use core::cmp::Ordering;
use std::fmt::{self, Display};
use std::num::NonZeroU64;

use la_stack::{LaError, Matrix, UnrepresentableReason, Vector};
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{FromPrimitive, Signed, ToPrimitive, Zero};

/// Number of matrices in each deterministic random benchmark corpus.
pub const RANDOM_INPUT_ARRAY_LEN: usize = 50;
/// Stable global seed used to derive one random corpus per dimension.
pub const RANDOM_SEED: [u8; 32] = [0; 32];

/// Configuration errors for exact-arithmetic benchmark input generation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum ExactBenchConfigError {
    /// The inclusive lower bound was greater than the inclusive upper bound.
    #[non_exhaustive]
    UnorderedRange {
        /// Inclusive lower bound.
        min: i16,
        /// Inclusive upper bound.
        max: i16,
    },
}

impl Display for ExactBenchConfigError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::UnorderedRange { min, max } => {
                write!(f, "random integer range must be ordered: {min}..={max}")
            }
        }
    }
}

impl std::error::Error for ExactBenchConfigError {}

/// Inclusive integer range used by the fixed-seed exact benchmark generator.
#[derive(Clone, Copy)]
#[must_use]
pub struct I16Range {
    min: i16,
    width: NonZeroU64,
}

impl I16Range {
    /// Validate an inclusive `i16` range and cache its sampling width.
    ///
    /// # Errors
    ///
    /// Returns [`ExactBenchConfigError::UnorderedRange`] when `min > max`.
    pub fn try_new(min: i16, max: i16) -> Result<Self, ExactBenchConfigError> {
        if min > max {
            return Err(ExactBenchConfigError::UnorderedRange { min, max });
        }

        let raw_width = u64::from(max.abs_diff(min)) + 1;
        // For ordered i16 bounds, the inclusive width is always 1..=65,536.
        // The `None` arm is therefore an internal invariant violation, not a
        // caller-reachable configuration error.
        let Some(width) = NonZeroU64::new(raw_width) else {
            unreachable!("an ordered inclusive i16 range has positive width");
        };
        Ok(Self { min, width })
    }
}

/// Deterministic `SplitMix64` generator for reproducible benchmark corpora.
#[must_use]
pub struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    /// Initialize the generator with a fixed state.
    pub const fn new(state: u64) -> Self {
        Self { state }
    }

    /// Advance the generator and return the next 64 random bits.
    const fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Draw a random `i16` inside a validated inclusive range.
    #[must_use]
    pub fn next_i16(&mut self, range: I16Range) -> i16 {
        #[expect(
            clippy::cast_possible_truncation,
            reason = "an inclusive i16 range has width at most 65,536, so its modulo offset fits i32"
        )]
        let offset = (self.next_u64() % range.width.get()) as i32;
        let value_i32 = i32::from(range.min) + offset;
        #[expect(
            clippy::cast_possible_truncation,
            reason = "the validated inclusive range guarantees min plus its modulo offset stays within i16"
        )]
        let value = value_i32 as i16;
        value
    }
}

/// Matrix/RHS pair used by exact-arithmetic benchmarks and smoke tests.
#[derive(Clone, Copy)]
#[must_use]
pub struct ExactInput<const D: usize> {
    /// Finite matrix under test.
    pub matrix: Matrix<D>,
    /// Finite right-hand side under test.
    pub rhs: Vector<D>,
}

/// Exact-arithmetic benchmark input whose results have been checked against
/// independent mathematical oracles.
///
/// Values of this type can only be produced by [`validate_exact_fixture`].
/// Keeping its fields private prevents Criterion helpers from accidentally
/// accepting a raw fixture whose preconditions have not been checked.
#[derive(Clone, Copy)]
#[must_use]
pub struct ValidatedExactInput<const D: usize> {
    matrix: Matrix<D>,
    rhs: Vector<D>,
}

impl<const D: usize> ValidatedExactInput<D> {
    /// Borrow the independently validated benchmark matrix.
    pub const fn matrix(&self) -> &Matrix<D> {
        &self.matrix
    }

    /// Return the independently validated benchmark right-hand side.
    pub const fn rhs(&self) -> Vector<D> {
        self.rhs
    }
}

/// Return a successful fixture-construction result or panic with context.
fn require_ok<T, E: Display>(result: Result<T, E>, operation: &str) -> T {
    match result {
        Ok(value) => value,
        Err(err) => panic!("{operation} failed: {err}"),
    }
}

/// Return one matrix entry through the bounds-checked API shared by current and
/// v0.4.3 releases.
fn stored_matrix_entry<const D: usize>(matrix: &Matrix<D>, row: usize, col: usize) -> f64 {
    matrix
        .get(row, col)
        .unwrap_or_else(|| panic!("matrix entry ({row}, {col}) is outside dimension {D}"))
}

/// Normalize the exact determinant-sign API across the v0.4.3 compatibility
/// boundary used only by historical benchmark worktrees.
#[cfg(not(la_stack_v0_4_3_api))]
fn checked_det_sign<const D: usize>(matrix: &Matrix<D>) -> i8 {
    matrix.det_sign_exact().as_i8()
}

#[cfg(la_stack_v0_4_3_api)]
fn checked_det_sign<const D: usize>(matrix: &Matrix<D>) -> i8 {
    require_ok(
        matrix.det_sign_exact(),
        "exact determinant sign oracle check",
    )
}

/// Return a deterministic, strictly diagonally-dominant matrix entry.
#[inline]
#[expect(
    clippy::cast_precision_loss,
    reason = "benchmark dimensions and indices are small enough to be represented exactly as f64"
)]
const fn matrix_entry<const D: usize>(r: usize, c: usize) -> f64 {
    if r == c {
        (r as f64).mul_add(1.0e-3, (D as f64) + 1.0)
    } else {
        0.1 / ((r + c + 1) as f64)
    }
}

/// Build the deterministic baseline matrix rows for dimension `D`.
#[inline]
#[must_use]
pub const fn make_matrix_rows<const D: usize>() -> [[f64; D]; D] {
    let mut rows = [[0.0; D]; D];
    let mut r = 0;
    while r < D {
        let mut c = 0;
        while c < D {
            rows[r][c] = matrix_entry::<D>(r, c);
            c += 1;
        }
        r += 1;
    }
    rows
}

/// Build the deterministic baseline right-hand-side vector for dimension `D`.
#[inline]
#[expect(
    clippy::cast_precision_loss,
    reason = "benchmark vector indices are small enough to be represented exactly as f64"
)]
#[must_use]
pub fn make_vector_array<const D: usize>() -> [f64; D] {
    from_fn(|i| (i as f64) + 1.0)
}

/// Derive a stable per-dimension seed from the global random benchmark seed.
fn random_seed_for_dim<const D: usize>() -> u64 {
    let mut seed =
        0xC0DE_CAFE_D15C_A11Au64 ^ require_ok(u64::try_from(D), "dimension seed conversion");
    for (i, byte) in RANDOM_SEED.iter().copied().enumerate() {
        let shift = require_ok(u32::try_from((i % 8) * 8), "seed shift conversion");
        seed ^= u64::from(byte) << shift;
        seed = seed.rotate_left(7) ^ require_ok(u64::try_from(i), "seed index conversion");
    }
    seed
}

/// Build a fixed random corpus of finite, strictly diagonally-dominant inputs.
pub fn make_random_input_corpus<const D: usize>() -> [ExactInput<D>; RANDOM_INPUT_ARRAY_LEN] {
    let mut rng = SplitMix64::new(random_seed_for_dim::<D>());
    let entry_range = require_ok(I16Range::try_new(-10, 10), "random integer range");
    array::from_fn(|_| {
        let mut rows = [[0.0; D]; D];
        let mut diag = [0_i16; D];

        for (r, row) in rows.iter_mut().enumerate() {
            for (c, entry) in row.iter_mut().enumerate() {
                if r == c {
                    diag[r] = rng.next_i16(entry_range);
                } else {
                    *entry = f64::from(rng.next_i16(entry_range));
                }
            }
        }

        let shift =
            f64::from(require_ok(u8::try_from(D), "dimension shift conversion")).mul_add(10.0, 1.0);
        for (i, row) in rows.iter_mut().enumerate() {
            row[i] = if diag[i] >= 0 {
                f64::from(diag[i]) + shift
            } else {
                f64::from(diag[i]) - shift
            };
        }

        let rhs = from_fn(|_| f64::from(rng.next_i16(entry_range)));

        ExactInput {
            matrix: require_ok(
                Matrix::<D>::try_from_rows(rows),
                "random matrix construction",
            ),
            rhs: require_ok(Vector::<D>::try_new(rhs), "random RHS vector construction"),
        }
    })
}

/// Build the fixed near-singular 3×3 benchmark input.
pub fn near_singular_3x3_input() -> ExactInput<3> {
    let perturbation = f64::from_bits(0x3CD0_0000_0000_0000); // 2^-50
    ExactInput {
        matrix: require_ok(
            Matrix::<3>::try_from_rows([
                [1.0 + perturbation, 2.0, 3.0],
                [4.0, 5.0, 6.0],
                [7.0, 8.0, 9.0],
            ]),
            "near-singular matrix construction",
        ),
        rhs: require_ok(
            Vector::<3>::try_new([1.0, 2.0, 3.0]),
            "near-singular RHS vector construction",
        ),
    }
}

/// Build the fixed extreme-magnitude 3×3 benchmark input.
pub fn large_entries_3x3_input() -> ExactInput<3> {
    let big = f64::MAX / 2.0;
    ExactInput {
        matrix: require_ok(
            Matrix::<3>::try_from_rows([[big, 1.0, 1.0], [1.0, big, 1.0], [1.0, 1.0, big]]),
            "large-entry matrix construction",
        ),
        rhs: require_ok(
            Vector::<3>::try_new([1.0, 1.0, 1.0]),
            "large-entry RHS vector construction",
        ),
    }
}

/// Build a Hilbert-matrix benchmark input with an all-ones RHS.
#[expect(
    clippy::cast_precision_loss,
    reason = "Hilbert benchmark dimensions and indices are small enough to be represented exactly as f64"
)]
pub fn hilbert_input<const D: usize>() -> ExactInput<D> {
    let rows = from_fn(|r| from_fn(|c| 1.0 / ((r + c + 1) as f64)));
    ExactInput {
        matrix: require_ok(
            Matrix::<D>::try_from_rows(rows),
            "Hilbert matrix construction",
        ),
        rhs: require_ok(
            Vector::<D>::try_new([1.0; D]),
            "Hilbert RHS vector construction",
        ),
    }
}

/// Convert one finite binary64 value to its exact rational value independently.
fn rational_from_f64(value: f64) -> BigRational {
    let Some(exact) = BigRational::from_f64(value) else {
        panic!("finite binary64 fixture {value:?} must convert exactly");
    };
    exact
}

/// Return whether one permutation has even parity.
fn permutation_is_even(perm: &[usize]) -> bool {
    let mut inversions = 0usize;
    for i in 0..perm.len() {
        for j in (i + 1)..perm.len() {
            if perm[i] > perm[j] {
                inversions += 1;
            }
        }
    }
    inversions.is_multiple_of(2)
}

/// Advance a slice to its next lexicographic permutation.
fn next_permutation(values: &mut [usize]) -> bool {
    if values.len() < 2 {
        return false;
    }

    let mut pivot = values.len() - 2;
    loop {
        if values[pivot] < values[pivot + 1] {
            break;
        }
        if pivot == 0 {
            return false;
        }
        pivot -= 1;
    }

    let mut successor = values.len() - 1;
    while values[successor] <= values[pivot] {
        successor -= 1;
    }
    values.swap(pivot, successor);
    values[(pivot + 1)..].reverse();
    true
}

/// Compute a determinant with the independent factorial-time Leibniz formula.
fn determinant_leibniz<const D: usize>(matrix: &Matrix<D>) -> BigRational {
    let mut determinant = BigRational::zero();
    let mut permutation: [usize; D] = from_fn(|index| index);

    loop {
        let mut term = BigRational::from_integer(BigInt::from(1));
        for (row, &col) in permutation.iter().enumerate() {
            term *= rational_from_f64(stored_matrix_entry(matrix, row, col));
        }
        if permutation_is_even(&permutation) {
            determinant += term;
        } else {
            determinant -= term;
        }
        if !next_permutation(&mut permutation) {
            break;
        }
    }

    determinant
}

/// Return the exact determinant sign implied by an independent rational value.
fn determinant_sign(exact: &BigRational) -> i8 {
    match exact.cmp(&BigRational::zero()) {
        Ordering::Less => -1,
        Ordering::Equal => 0,
        Ordering::Greater => 1,
    }
}

/// Derive the strict finite-binary64 outcome from an independently checked rational.
fn expected_strict_f64(exact: &BigRational) -> Result<f64, UnrepresentableReason> {
    let Some(rounded) = exact.to_f64() else {
        return Err(UnrepresentableReason::NotFinite);
    };
    if !rounded.is_finite() {
        return Err(UnrepresentableReason::NotFinite);
    }
    if BigRational::from_f64(rounded).as_ref() == Some(exact) {
        Ok(rounded)
    } else {
        Err(UnrepresentableReason::RequiresRounding)
    }
}

/// Return the exact overflow midpoint for binary64 round-to-nearest.
///
/// Magnitudes strictly below this value round to a finite binary64 value. The
/// midpoint itself rounds to infinity because `f64::MAX` has an odd least
/// significand bit while the hypothetical `2^1024` endpoint is even.
fn finite_rounding_limit() -> BigRational {
    let half_max_ulp = BigRational::from_integer(BigInt::from(1) << 970usize);
    rational_from_f64(f64::MAX) + half_max_ulp
}

/// Verify directly that a finite binary64 result is the nearest-even value.
fn assert_nearest_even_f64(actual: f64, exact: &BigRational) {
    assert!(actual.is_finite());
    if actual == 0.0 {
        assert_eq!(actual.is_sign_negative(), exact.is_negative());
    }

    let actual_exact = rational_from_f64(actual);
    let actual_distance = (&actual_exact - exact).abs();
    for neighbor in [actual.next_down(), actual.next_up()] {
        if !neighbor.is_finite() {
            continue;
        }
        let neighbor_distance = (rational_from_f64(neighbor) - exact).abs();
        assert!(
            actual_distance <= neighbor_distance,
            "rounded value {actual:?} is farther from {exact} than adjacent value {neighbor:?}",
        );
        if actual_distance == neighbor_distance {
            assert_eq!(
                actual.to_bits() & 1,
                0,
                "halfway value must select the even binary64 significand"
            );
        }
    }
}

/// Check one strict scalar conversion against an independently derived outcome.
fn assert_strict_scalar(
    actual: Result<f64, LaError>,
    exact: &BigRational,
    expected_index: Option<usize>,
) {
    match (actual, expected_strict_f64(exact)) {
        (Ok(actual), Ok(expected)) => assert_eq!(actual.to_bits(), expected.to_bits()),
        (Err(LaError::Unrepresentable { index, reason, .. }), Err(expected_reason)) => {
            assert_eq!(index, expected_index);
            assert_eq!(reason, expected_reason);
        }
        (actual, expected) => {
            panic!("strict conversion mismatch: actual={actual:?}, expected={expected:?}")
        }
    }
}

/// Check one rounded scalar conversion against the exact rational oracle.
fn assert_rounded_scalar(actual: Result<f64, LaError>, exact: &BigRational) {
    if exact.abs() < finite_rounding_limit() {
        match actual {
            Ok(actual) => assert_nearest_even_f64(actual, exact),
            Err(error) => panic!("finite nearest-even conversion failed: {error}"),
        }
    } else {
        assert!(matches!(
            actual,
            Err(LaError::Unrepresentable {
                reason: UnrepresentableReason::NotFinite,
                ..
            })
        ));
    }
}

/// Verify `A · x = b` exactly using independently reconstructed binary64 inputs.
fn assert_exact_residual<const D: usize>(input: &ExactInput<D>, solution: &[BigRational; D]) {
    for row in 0..D {
        let mut observed = BigRational::zero();
        for (col, value) in solution.iter().enumerate() {
            observed += rational_from_f64(stored_matrix_entry(&input.matrix, row, col)) * value;
        }
        assert_eq!(observed, rational_from_f64(input.rhs.as_array()[row]));
    }
}

/// Validate every exact benchmark operation against independent mathematical evidence.
///
/// The returned proof-bearing fixture is the only input accepted by Criterion
/// registration and timed-operation helpers.
///
/// This runs only during benchmark setup and smoke tests, never inside a timed
/// Criterion closure.
///
/// # Panics
///
/// Panics when any benchmark operation disagrees with the independent exact
/// oracle or when a fixture unexpectedly violates an operation precondition.
pub fn validate_exact_fixture<const D: usize>(input: ExactInput<D>) -> ValidatedExactInput<D> {
    let determinant = determinant_leibniz(&input.matrix);
    assert_eq!(
        require_ok(input.matrix.det_exact(), "exact determinant oracle check"),
        determinant
    );
    assert_eq!(
        checked_det_sign(&input.matrix),
        determinant_sign(&determinant)
    );
    assert_strict_scalar(input.matrix.det_exact_f64(), &determinant, None);
    assert_rounded_scalar(input.matrix.det_exact_rounded_f64(), &determinant);

    let solution = require_ok(
        input.matrix.solve_exact(input.rhs),
        "exact solve oracle check",
    );
    assert_exact_residual(&input, &solution);

    let strict_solution = input.matrix.solve_exact_f64(input.rhs);
    let first_failure = solution.iter().enumerate().find_map(|(index, value)| {
        expected_strict_f64(value)
            .err()
            .map(|reason| (index, reason))
    });
    match (strict_solution, first_failure) {
        (Ok(actual), None) => {
            for (index, exact) in solution.iter().enumerate() {
                let Ok(expected) = expected_strict_f64(exact) else {
                    panic!("strict solution component {index} unexpectedly requires rounding");
                };
                assert_eq!(actual.as_array()[index].to_bits(), expected.to_bits());
            }
        }
        (
            Err(LaError::Unrepresentable {
                index: Some(index),
                reason,
                ..
            }),
            Some((expected_index, expected_reason)),
        ) => {
            assert_eq!(index, expected_index);
            assert_eq!(reason, expected_reason);
        }
        (actual, expected) => {
            panic!(
                "strict exact-solve conversion mismatch: actual={actual:?}, expected={expected:?}"
            )
        }
    }

    let rounding_limit = finite_rounding_limit();
    let first_rounded_failure = solution
        .iter()
        .position(|exact| exact.abs() >= rounding_limit);
    match (
        input.matrix.solve_exact_rounded_f64(input.rhs),
        first_rounded_failure,
    ) {
        (Ok(rounded), None) => {
            for (actual, exact) in rounded.as_array().iter().copied().zip(&solution) {
                assert_nearest_even_f64(actual, exact);
            }
        }
        (
            Err(LaError::Unrepresentable {
                index: Some(index),
                reason: UnrepresentableReason::NotFinite,
                ..
            }),
            Some(expected_index),
        ) => assert_eq!(index, expected_index),
        (actual, expected) => {
            panic!(
                "rounded exact-solve conversion mismatch: actual={actual:?}, expected failing index={expected:?}"
            )
        }
    }

    ValidatedExactInput {
        matrix: input.matrix,
        rhs: input.rhs,
    }
}
