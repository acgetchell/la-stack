#![forbid(unsafe_code)]

//! Shared helpers for the `vs_linalg` benchmark and its smoke tests.

use faer::linalg::solvers::{Ldlt as FaerLdlt, PartialPivLu};
use faer::perm::PermRef;
use nalgebra::SMatrix;

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

/// Compute a determinant from a faer partial-pivot LU factorization.
#[must_use]
pub fn faer_det_from_partial_piv_lu(lu: &PartialPivLu<f64>) -> f64 {
    // For PA = LU with unit-lower L, det(A) = det(P) * det(U).
    let u = lu.U();
    let mut det = 1.0;
    for i in 0..u.nrows() {
        det *= u[(i, i)];
    }
    det * faer_perm_sign(lu.P())
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

/// Build a positive-definite diagonal matrix spanning 112 binary exponents at D=8.
///
/// Each successive pivot is `2^-16` times the previous one. Benchmarks use a
/// zero tolerance so the complete, finite factorization remains in scope.
#[inline]
#[must_use]
pub fn make_ill_conditioned_matrix_rows<const D: usize>() -> [[f64; D]; D] {
    let mut rows = [[0.0; D]; D];
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
