//! Shared helpers for the `vs_linalg` benchmark and its smoke tests.

use faer::linalg::solvers::{Ldlt as FaerLdlt, PartialPivLu};
use faer::perm::PermRef;
use nalgebra::SMatrix;

/// Return `det(P)` for faer's permutation representation.
///
/// Sign(det(P)) is +1 for even permutations and -1 for odd. Parity is computed
/// from the number of cycles: `sign = (-1)^(n - cycles)`.
pub fn faer_perm_sign(p: PermRef<'_, usize>) -> f64 {
    let (forward, _inverse) = p.arrays();
    let n = forward.len();

    let mut seen = vec![false; n];
    let mut cycles = 0usize;

    for start in 0..n {
        if seen[start] {
            continue;
        }
        cycles += 1;

        let mut i = start;
        while !seen[i] {
            seen[i] = true;
            i = forward[i];
        }
    }

    if (n - cycles).is_multiple_of(2) {
        1.0
    } else {
        -1.0
    }
}

/// Compute a determinant from a faer partial-pivot LU factorization.
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
#[allow(clippy::cast_precision_loss)] // D, r, c are small integers, precision loss is not an issue.
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
pub fn make_matrix_rows<const D: usize>() -> [[f64; D]; D] {
    let mut rows = [[0.0; D]; D];

    for (r, row) in rows.iter_mut().enumerate() {
        for (c, entry) in row.iter_mut().enumerate() {
            *entry = matrix_entry::<D>(r, c);
        }
    }

    rows
}

/// Return a deterministic benchmark vector entry.
#[inline]
#[allow(clippy::cast_precision_loss)] // i is a small integer, precision loss is not an issue.
pub fn vector_entry(i: usize, offset: f64) -> f64 {
    (i as f64) + 1.0 + offset
}

/// Build the shared vector input used by all crates for a dimension.
#[inline]
pub fn make_vector_array<const D: usize>(offset: f64) -> [f64; D] {
    let mut data = [0.0; D];

    for (i, entry) in data.iter_mut().enumerate() {
        *entry = vector_entry(i, offset);
    }

    data
}

/// Compute nalgebra's matrix infinity norm using la-stack's row-sum convention.
#[inline]
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
