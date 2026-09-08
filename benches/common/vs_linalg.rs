#![forbid(unsafe_code)]

//! Shared helpers for the `vs_linalg` benchmark and its smoke tests.

use core::array::from_fn;

use faer::linalg::solvers::{Ldlt as FaerLdlt, PartialPivLu};
use faer::perm::PermRef;
use nalgebra::SMatrix;

use la_stack::{LaError, Matrix, Tolerance, Vector};

#[cfg(not(la_stack_v0_4_3_api))]
use crate::bench_utils::OrAbort;

/// Evaluate la-stack's dot product through the ownership contract used by the
/// selected library revision.
///
/// # Errors
///
/// Returns the selected revision's typed error if finite inputs overflow during
/// dot-product accumulation.
#[cfg(not(la_stack_v0_4_3_api))]
#[inline]
pub fn la_stack_dot<const D: usize>(left: &Vector<D>, right: &Vector<D>) -> Result<f64, LaError> {
    left.dot(right)
}

/// Evaluate the v0.4.3 by-value dot-product API without changing benchmark
/// inputs or the mathematical operation.
///
/// # Errors
///
/// Returns v0.4.3's typed error if finite inputs overflow during dot-product
/// accumulation.
#[cfg(la_stack_v0_4_3_api)]
#[inline]
pub fn la_stack_dot<const D: usize>(left: &Vector<D>, right: &Vector<D>) -> Result<f64, LaError> {
    (*left).dot(*right)
}

/// Evaluate the squared norm through the current public API.
///
/// # Errors
///
/// Returns the library's typed error if squared-norm accumulation overflows.
#[cfg(not(any(la_stack_pre_rational_input_api, la_stack_v0_4_3_api)))]
#[inline]
pub const fn la_stack_norm_squared<const D: usize>(vector: &Vector<D>) -> Result<f64, LaError> {
    vector.norm_squared()
}

/// Evaluate the same squared norm through the pre-v0.4.6 method name.
///
/// This adapter preserves historical comparisons without adding library aliases.
///
/// # Errors
///
/// Returns the selected revision's typed error if squared-norm accumulation overflows.
#[cfg(any(la_stack_pre_rational_input_api, la_stack_v0_4_3_api))]
#[inline]
pub const fn la_stack_norm_squared<const D: usize>(vector: &Vector<D>) -> Result<f64, LaError> {
    vector.norm2_sq()
}

/// Evaluate the matrix infinity norm through the current public API.
///
/// # Errors
///
/// Returns the library's typed error if an absolute row sum overflows.
#[cfg(not(any(la_stack_pre_rational_input_api, la_stack_v0_4_3_api)))]
#[inline]
pub const fn la_stack_norm_inf<const D: usize>(matrix: &Matrix<D>) -> Result<f64, LaError> {
    matrix.norm_inf()
}

/// Evaluate the same matrix infinity norm through the pre-v0.4.6 method name.
///
/// This adapter preserves historical comparisons without adding library aliases.
///
/// # Errors
///
/// Returns the selected revision's typed error if an absolute row sum overflows.
#[cfg(any(la_stack_pre_rational_input_api, la_stack_v0_4_3_api))]
#[inline]
pub const fn la_stack_norm_inf<const D: usize>(matrix: &Matrix<D>) -> Result<f64, LaError> {
    matrix.inf_norm()
}

/// Parse a tolerance through the constructor exposed by the selected library
/// revision.
///
/// # Errors
///
/// Returns a typed error when `value` is negative or non-finite.
#[cfg(not(la_stack_v0_4_3_api))]
#[inline]
pub const fn la_stack_tolerance(value: f64) -> Result<Tolerance, LaError> {
    Tolerance::try_new(value)
}

/// Parse a tolerance through v0.4.3's pre-`try_` constructor name.
///
/// # Errors
///
/// Returns a typed error when `value` is negative or non-finite.
#[cfg(la_stack_v0_4_3_api)]
#[inline]
pub const fn la_stack_tolerance(value: f64) -> Result<Tolerance, LaError> {
    Tolerance::new(value)
}

/// Return `det(P)` for faer's permutation representation.
///
/// Sign(det(P)) is +1 for even permutations and -1 for odd. Parity is computed
/// from the number of cycles: `sign = (-1)^(n - cycles)`. The benchmark
/// dimensions use an allocation-free bit mask; larger permutations use an
/// allocation-free inversion count so the function remains total.
#[must_use]
pub fn faer_perm_sign(p: PermRef<'_, usize>) -> f64 {
    let (forward, _inverse) = p.arrays();
    let is_odd = if forward.len() <= u128::BITS as usize {
        permutation_is_odd_by_cycles(forward)
    } else {
        permutation_is_odd_by_inversions(forward)
    };

    if is_odd { -1.0 } else { 1.0 }
}

/// Return whether a permutation of at most 128 elements is odd.
fn permutation_is_odd_by_cycles(forward: &[usize]) -> bool {
    let mut seen = 0u128;
    let mut cycles = 0usize;

    for start in 0..forward.len() {
        if seen & (1u128 << start) != 0 {
            continue;
        }
        cycles += 1;

        let mut index = start;
        while seen & (1u128 << index) == 0 {
            seen |= 1u128 << index;
            index = forward[index];
        }
    }

    !(forward.len() - cycles).is_multiple_of(2)
}

/// Return whether a permutation is odd using an allocation-free fallback.
fn permutation_is_odd_by_inversions(forward: &[usize]) -> bool {
    let mut is_odd = false;

    for (index, &left) in forward.iter().enumerate() {
        for &right in &forward[index + 1..] {
            is_odd ^= left > right;
        }
    }

    is_odd
}

/// Borrow a faer partial-pivot LU factorization together with its precomputed
/// permutation sign.
///
/// The private sign field is derived from the borrowed factorization, so callers
/// cannot accidentally combine one LU decomposition with another permutation.
#[must_use]
pub struct PreparedFaerLuDet<'a> {
    lu: &'a PartialPivLu<f64>,
    permutation_sign: f64,
}

impl<'a> PreparedFaerLuDet<'a> {
    /// Prepare repeated determinant queries for one borrowed LU factorization.
    pub fn new(lu: &'a PartialPivLu<f64>) -> Self {
        Self {
            lu,
            permutation_sign: faer_perm_sign(lu.P()),
        }
    }

    /// Compute the determinant from the prepared factorization.
    #[must_use]
    pub fn det(&self) -> f64 {
        // For PA = LU with unit-lower L, det(A) = det(P) * det(U).
        let u = self.lu.U();
        let mut det = 1.0;
        for i in 0..u.nrows() {
            det *= u[(i, i)];
        }
        det * self.permutation_sign
    }
}

/// Compute a determinant from a faer LDLT factorization.
#[must_use]
pub fn faer_det_from_ldlt(ldlt: &FaerLdlt<f64>) -> f64 {
    let d = ldlt.D().column_vector();
    let mut det = 1.0;
    for i in 0..d.nrows() {
        det *= d[i];
    }
    det
}

/// Return a deterministic, strictly diagonally-dominant benchmark matrix entry.
#[inline]
#[expect(
    clippy::cast_precision_loss,
    reason = "benchmark dimensions and indices are small enough to be represented exactly as f64"
)]
#[must_use]
pub fn matrix_entry<const D: usize>(r: usize, c: usize) -> f64 {
    if r == c {
        // Strict diagonal dominance for stability.
        (r as f64).mul_add(1.0e-3, (D as f64) + 1.0)
    } else {
        // Small, varying off-diagonals.
        0.1 / ((r + c + 1) as f64)
    }
}

/// Build the shared matrix rows used by all crates for a dimension.
#[inline]
#[must_use]
pub fn make_matrix_rows<const D: usize>() -> [[f64; D]; D] {
    let mut rows = [[0.0; D]; D];

    for (r, row) in rows.iter_mut().enumerate() {
        for (c, entry) in row.iter_mut().enumerate() {
            *entry = matrix_entry::<D>(r, c);
        }
    }

    rows
}

/// Build a well-conditioned matrix whose first LU column requires a row swap.
#[inline]
#[must_use]
pub fn make_pivoting_matrix_rows<const D: usize>() -> [[f64; D]; D] {
    let mut rows = make_matrix_rows();
    if D > 1 {
        rows.swap(0, 1);
    }
    rows
}

/// Diagnostic solve families, separate from the headline comparison matrix.
#[derive(Clone, Copy, Debug)]
#[cfg(not(la_stack_v0_4_3_api))]
pub enum LuSolveScenario {
    /// Rotate the well-conditioned rows, requiring repeated LU row swaps.
    Pivoting,
    /// Dense SPD matrix J + 2^-16 I, with condition number 1 + D * 2^16.
    DenseIllConditioned,
}

#[cfg(not(la_stack_v0_4_3_api))]
impl LuSolveScenario {
    /// Stable diagnostic group label.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Pivoting => "pivoting",
            Self::DenseIllConditioned => "dense_ill_conditioned",
        }
    }
}

/// A diagnostic system checked against a known solution and a scaled residual.
#[must_use]
#[cfg(not(la_stack_v0_4_3_api))]
pub struct ValidatedLuSolveInput<const D: usize> {
    matrix: Matrix<D>,
    rhs: Vector<D>,
    expected: [f64; D],
}

#[cfg(not(la_stack_v0_4_3_api))]
impl<const D: usize> ValidatedLuSolveInput<D> {
    /// Borrow the finite matrix.
    pub const fn matrix(&self) -> &Matrix<D> {
        &self.matrix
    }

    /// Return the finite right-hand side.
    pub const fn rhs(&self) -> Vector<D> {
        self.rhs
    }

    /// Validate another implementation before timing it on this system.
    ///
    /// # Panics
    /// Panics if a non-finite value, excessive forward error, or excessive
    /// scaled residual is observed.
    pub fn validate_solution(&self, solution: &[f64; D]) {
        let mut residual = 0.0_f64;
        let mut matrix_norm = 0.0_f64;
        let mut solution_norm = 0.0_f64;
        let mut rhs_norm = 0.0_f64;
        for (&actual, &expected) in solution.iter().zip(&self.expected) {
            assert!(actual.is_finite());
            assert!((actual - expected).abs() <= 1e-7 * expected.abs().max(1.0));
            solution_norm = solution_norm.max(actual.abs());
        }
        for (row, &rhs) in self.matrix.as_rows().iter().zip(self.rhs.as_array()) {
            // Deliberately use an ordinary product/sum residual, independently
            // of the timed triangular FMA recurrence.
            let observed: f64 = row.iter().zip(solution).map(|(&a, &x)| a * x).sum();
            residual = residual.max((observed - rhs).abs());
            matrix_norm = matrix_norm.max(row.iter().map(|value| value.abs()).sum());
            rhs_norm = rhs_norm.max(rhs.abs());
        }
        assert!(residual <= 1e-12 * matrix_norm.mul_add(solution_norm, rhs_norm));
    }
}

/// Construct a diagnostic system and validate its complete LU solve.
///
/// The dense family uses an exactly representable RHS for x[i] = i+1;
/// it exercises cancellation without diagonal or sparse shortcuts.
///
/// # Panics
/// Panics for dimensions outside the diagnostic domain or a failed oracle check.
#[cfg(not(la_stack_v0_4_3_api))]
pub fn validated_lu_solve_input<const D: usize>(
    scenario: LuSolveScenario,
) -> ValidatedLuSolveInput<D> {
    assert!((2..=64).contains(&D));
    let expected = from_fn(|i| f64::from(u32::try_from(i + 1).or_abort("solution index")));
    let (rows, rhs) = match scenario {
        LuSolveScenario::Pivoting => {
            let mut rows = make_matrix_rows::<D>();
            rows.rotate_left(1);
            let rhs = rows.map(|row| row.iter().zip(&expected).map(|(&a, &x)| a * x).sum());
            (rows, rhs)
        }
        LuSolveScenario::DenseIllConditioned => {
            const DELTA: f64 = 1.0 / 65_536.0;
            let rows = from_fn(|i| from_fn(|j| if i == j { 1.0 + DELTA } else { 1.0 }));
            let total = f64::from(u32::try_from(D * (D + 1) / 2).or_abort("solution sum"));
            let rhs = expected.map(|value| DELTA.mul_add(value, total));
            (rows, rhs)
        }
    };
    let input = ValidatedLuSolveInput {
        matrix: Matrix::try_from_rows(rows).or_abort("diagnostic matrix"),
        rhs: Vector::try_new(rhs).or_abort("diagnostic RHS"),
        expected,
    };
    let solution = input
        .matrix
        .lu(la_stack_tolerance(0.0).or_abort("diagnostic tolerance"))
        .or_abort("diagnostic LU")
        .solve(input.rhs)
        .or_abort("diagnostic solve");
    input.validate_solution(solution.as_array());
    input
}

/// Build a positive-definite diagonal matrix spanning 112 binary exponents at D=8.
///
/// Each successive pivot is `2^-16` times the previous one. Benchmarks use a
/// zero tolerance so the complete, finite factorization remains in scope. The
/// fixed return dimension prevents extending the progression until a diagonal
/// entry underflows to zero and destroys positive-definiteness.
#[inline]
#[must_use]
pub fn make_ill_conditioned_matrix_rows() -> [[f64; 8]; 8] {
    let mut rows = [[0.0; 8]; 8];
    let mut diagonal = 1.0;
    for (index, row) in rows.iter_mut().enumerate() {
        row[index] = diagonal;
        diagonal *= 1.0 / 65_536.0;
    }
    rows
}

/// Build a positive diagonal matrix whose factors underflow in sequential order
/// but whose balanced exact product is one when `D` is a multiple of four.
#[inline]
#[must_use]
pub fn make_balanced_dynamic_range_rows<const D: usize>() -> [[f64; D]; D] {
    const TWO_NEG_800: f64 = f64::from_bits(223_u64 << 52);
    const TWO_POS_800: f64 = f64::from_bits(1823_u64 << 52);

    let mut rows = [[0.0; D]; D];
    for (index, row) in rows.iter_mut().enumerate() {
        row[index] = if index % 4 < 2 {
            TWO_NEG_800
        } else {
            TWO_POS_800
        };
    }
    rows
}

/// Return a deterministic benchmark vector entry.
#[inline]
#[expect(
    clippy::cast_precision_loss,
    reason = "benchmark vector indices are small enough to be represented exactly as f64"
)]
#[must_use]
pub fn vector_entry(i: usize, offset: f64) -> f64 {
    (i as f64) + 1.0 + offset
}

/// Build the shared vector input used by all crates for a dimension.
#[inline]
#[must_use]
pub fn make_vector_array<const D: usize>(offset: f64) -> [f64; D] {
    let mut data = [0.0; D];

    for (i, entry) in data.iter_mut().enumerate() {
        *entry = vector_entry(i, offset);
    }

    data
}

/// Build a norm input whose non-zero magnitudes decrease after the first entry.
#[inline]
#[must_use]
pub fn make_norm_descending_array<const D: usize>() -> [f64; D] {
    from_fn(|index| {
        let magnitude = vector_entry(D - index - 1, 0.0);
        if index % 2 == 0 {
            magnitude
        } else {
            -magnitude
        }
    })
}

/// Build a norm input whose entries repeatedly have the same non-zero magnitude.
#[inline]
#[must_use]
pub fn make_norm_repeated_scale_array<const D: usize>() -> [f64; D] {
    from_fn(|index| if index % 2 == 0 { 3.0 } else { -3.0 })
}

/// Build a norm input with one non-zero entry and otherwise skipped zeros.
#[inline]
#[must_use]
pub fn make_norm_sparse_array<const D: usize>() -> [f64; D] {
    let mut data = [0.0; D];
    if D > 0 {
        data[D / 2] = -3.0;
    }
    data
}

/// Build a finite norm input spanning normal, subnormal, and zero magnitudes.
#[inline]
#[must_use]
pub fn make_norm_wide_dynamic_range_array<const D: usize>() -> [f64; D] {
    const VALUES: [f64; 8] = [
        1.0e200,
        -1.0e-200,
        f64::from_bits(1),
        0.0,
        -1.0e100,
        1.0e-100,
        1.0,
        -1.0e200,
    ];

    from_fn(|index| VALUES[index % VALUES.len()])
}

/// Build the named Euclidean-norm scenario corpus in stable benchmark order.
#[inline]
#[must_use]
pub fn norm_scenarios<const D: usize>() -> [(&'static str, [f64; D]); 4] {
    [
        ("descending", make_norm_descending_array()),
        ("repeated_scale", make_norm_repeated_scale_array()),
        ("sparse", make_norm_sparse_array()),
        ("wide_dynamic_range", make_norm_wide_dynamic_range_array()),
    ]
}

/// Compute a Euclidean norm by chaining the standard library's binary `hypot`.
///
/// This is an independent overflow- and underflow-safe reference kernel for the
/// focused `norm2` benchmarks.
#[inline]
#[must_use]
pub fn iterative_hypot<const D: usize>(values: &[f64; D]) -> f64 {
    values.iter().fold(0.0, |norm, &value| norm.hypot(value))
}

/// Reproduce the generalized scaled norm currently owned by Delaunay.
///
/// D=2 delegates to binary `hypot`; larger dimensions make one pass for the
/// maximum magnitude and a second pass for the scaled sum of squares.
#[inline]
#[must_use]
pub fn delaunay_scaled_norm<const D: usize>(values: &[f64; D]) -> f64 {
    match D {
        0 => 0.0,
        1 => values[0].abs(),
        2 => values[0].hypot(values[1]),
        _ => {
            if values.iter().any(|value| value.is_nan()) {
                return f64::NAN;
            }
            if values.iter().any(|value| value.is_infinite()) {
                return f64::INFINITY;
            }

            let scale = values.iter().map(|value| value.abs()).fold(0.0, f64::max);
            if scale == 0.0 {
                return 0.0;
            }

            let scaled_sum = values.iter().fold(0.0, |sum, &value| {
                let ratio = value / scale;
                ratio.mul_add(ratio, sum)
            });
            scale * scaled_sum.sqrt()
        }
    }
}

/// Compute nalgebra's matrix infinity norm using la-stack's row-sum convention.
#[inline]
#[must_use]
pub fn nalgebra_inf_norm<const D: usize>(m: &SMatrix<f64, D, D>) -> f64 {
    // Infinity norm = max absolute row sum.
    let mut max_row_sum = 0.0;

    for r in 0..D {
        let mut row_sum = 0.0;
        for c in 0..D {
            row_sum += m[(r, c)].abs();
        }
        if row_sum > max_row_sum {
            max_row_sum = row_sum;
        }
    }

    max_row_sum
}
