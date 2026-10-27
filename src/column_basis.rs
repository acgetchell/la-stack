#![forbid(unsafe_code)]

//! Provisional numerical column selection for fixed rectangular arrays.

use crate::{ArithmeticOperation, LaError, Tolerance};

/// Numerical column selection for an `R`-row, `C`-column input.
///
/// Created by [`column_basis`]. Both input dimensions remain part of the type.
/// The selected indices are unique, strictly increasing, and less than `C`;
/// their count is at most `min(R, C)`. Storage is stack-backed, with capacity
/// `R`. Neither the indices nor their count certify exact independence or rank.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[must_use]
pub struct ColumnBasis<const R: usize, const C: usize> {
    columns: [usize; R],
    rank: usize,
}

impl<const R: usize, const C: usize> ColumnBasis<R, C> {
    /// Borrow selected column indices in increasing input-column order.
    #[must_use]
    pub const fn columns(&self) -> &[usize] {
        self.columns.split_at(self.rank).0
    }

    /// Number of pivots accepted under the requested absolute tolerance.
    ///
    /// This is a provisional numerical rank, not an exact rank certificate,
    /// including when the requested tolerance is zero.
    #[must_use]
    pub const fn rank(&self) -> usize {
        self.rank
    }
}

/// Select columns of a rectangular array by partial-pivoted elimination.
///
/// Columns are visited from left to right. At each step the largest absolute
/// entry among remaining rows is the candidate pivot; ties select the first
/// row in the current row order. A pivot is accepted only when its magnitude
/// is **strictly greater** than `tolerance.get()`. Rejected columns are skipped;
/// selection ends after all columns or `R` accepted pivots. Empty and zero
/// matrices return rank zero; tall, square, and wide shapes are supported.
///
/// The tolerance is absolute, in input coefficient units. There is no automatic
/// scaling or normalization. To use a relative policy, compute a finite absolute
/// threshold from the original input scale and pass it through
/// [`Tolerance::try_new`]. Rescaling rows or columns can change the selection.
/// Elimination uses partial row pivoting and fused multiply-add updates
/// (Gaussian elimination, `REFERENCES.md` \[11, 12\]). No absolute rounding-error
/// bound is provided. Rounding, underflow, growth, and ill-conditioning can alter
/// the selected columns and rank even at zero tolerance. The result is provisional:
/// exact certification of a square selected basis remains the caller's job.
/// Failure to find a pivot never produces an exact singularity claim.
///
/// The input is borrowed and preserved. A local `R × C` copy and `R` indices
/// use stack storage; no allocation or optional feature is required. Dimensions
/// have no dispatch cap. Work is `O(R (1 + C (1 + min(R, C))))`, including
/// validation and initialization for empty shapes.
///
/// # Errors
/// All inputs are checked in row-major order before elimination, including
/// columns beyond a full basis. NaN/infinity returns [`LaError::NonFinite`]
/// with input origin and the original matrix cell. A non-finite update returns
/// the same error with [`ArithmeticOperation::ColumnBasis`] and its cell in the
/// row-permuted work array. Insufficient numerical rank is a successful partial
/// selection, not an error. No partially computed result is returned on error.
///
/// # Examples
/// ```
/// use la_stack::prelude::*;
/// # fn main() -> Result<(), LaError> {
/// let rows = [[0.0, 1.0, 2.0, 0.0], [0.0, 2.0, 4.0, 1.0]];
/// let basis: ColumnBasis<2, 4> = column_basis(&rows, Tolerance::try_new(1e-12)?)?;
/// assert_eq!(basis.columns(), &[1, 3]);
/// assert_eq!(basis.rank(), 2);
/// # Ok(())
/// # }
/// ```
pub const fn column_basis<const R: usize, const C: usize>(
    rows: &[[f64; C]; R],
    tolerance: Tolerance,
) -> Result<ColumnBasis<R, C>, LaError> {
    let mut row = 0;
    while row < R {
        let mut col = 0;
        while col < C {
            if !rows[row][col].is_finite() {
                return Err(LaError::non_finite_input_matrix(row, col));
            }
            col += 1;
        }
        row += 1;
    }
    let mut work = *rows;
    let mut basis = ColumnBasis {
        columns: [0; R],
        rank: 0,
    };
    let mut col = 0;
    while col < C && basis.rank < R {
        let pivot_row = basis.rank;
        let mut candidate = pivot_row;
        row = pivot_row + 1;
        while row < R {
            if work[row][col].abs() > work[candidate][col].abs() {
                candidate = row;
            }
            row += 1;
        }
        if work[candidate][col].abs() > tolerance.get() {
            work.swap(pivot_row, candidate);
            let pivot = work[pivot_row][col];
            row = pivot_row + 1;
            while row < R {
                // Partial pivoting bounds the multiplier by one. A nonzero
                // accepted finite pivot makes division finite, even if it underflows.
                let factor = work[row][col] / pivot;
                work[row][col] = 0.0;
                let mut trailing = col + 1;
                while trailing < C {
                    let value = (-factor).mul_add(work[pivot_row][trailing], work[row][trailing]);
                    if !value.is_finite() {
                        return Err(LaError::non_finite_computation_matrix(
                            ArithmeticOperation::ColumnBasis,
                            row,
                            trailing,
                        ));
                    }
                    work[row][trailing] = value;
                    trailing += 1;
                }
                row += 1;
            }
            basis.columns[basis.rank] = col;
            basis.rank += 1;
        }
        col += 1;
    }
    Ok(basis)
}
