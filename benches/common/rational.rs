#![forbid(unsafe_code)]

//! Independently validated rational inputs shared by timing and allocation probes.

use la_stack::{BigInt, BigRational, DeterminantSign, RationalMatrix, RationalVector};

use super::bench_utils::OrAbort;

/// Fixed rational component-size families.
#[derive(Clone, Copy, Debug)]
pub enum RationalInputKind {
    /// Original small-integer components from the rational-input benchmark.
    Small,
    /// Row factors with at least 256-bit numerator/denominator inputs.
    Wide256,
    /// Row factors with at least 1024-bit numerator/denominator inputs.
    Wide1024,
}

impl RationalInputKind {
    /// Stable fixture families, in registration order.
    pub const ALL: [Self; 3] = [Self::Small, Self::Wide256, Self::Wide1024];

    /// Stable family name used in benchmark and allocation output.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Small => "small",
            Self::Wide256 => "wide256",
            Self::Wide1024 => "wide1024",
        }
    }

    const fn component_bits(self) -> Option<usize> {
        match self {
            Self::Small => None,
            Self::Wide256 => Some(256),
            Self::Wide1024 => Some(1024),
        }
    }
}

/// Independently validated exact-input rational benchmark fixture.
#[must_use]
pub struct ValidatedRationalInput<const D: usize> {
    matrix: RationalMatrix<D>,
    rhs: RationalVector<D>,
}

impl<const D: usize> ValidatedRationalInput<D> {
    /// Borrow the independently checked matrix.
    pub const fn matrix(&self) -> &RationalMatrix<D> {
        &self.matrix
    }

    /// Borrow the independently checked right-hand side.
    pub const fn rhs(&self) -> &RationalVector<D> {
        &self.rhs
    }
}

/// Deterministic, strictly diagonally-dominant exact-input rational fixture.
///
/// # Panics
/// Panics if construction or the independent determinant/solve checks fail.
pub fn rational_input<const D: usize>(kind: RationalInputKind) -> ValidatedRationalInput<D> {
    let mut rows = std::array::from_fn(|row| {
        std::array::from_fn(|col| {
            if row == col {
                let diagonal = 2 * D + row + 1;
                BigRational::from_integer(BigInt::from(diagonal))
            } else {
                let raw_numerator = (row * 3 + col * 5) % 5;
                let numerator = i64::try_from(raw_numerator).or_abort("rational numerator") - 2;
                let denominator = (row + col) % 7 + 2;
                BigRational::new(BigInt::from(numerator), BigInt::from(denominator))
            }
        })
    });
    if let Some(bits) = kind.component_bits() {
        for (row, entries) in rows.iter_mut().enumerate() {
            let shift = bits + 16 * row;
            let factor = BigRational::new(
                (BigInt::from(1) << shift) + BigInt::from(1),
                (BigInt::from(1) << (shift + 1)) - BigInt::from(1),
            );
            for entry in entries {
                *entry *= &factor;
            }
        }
    }
    let expected_solution = std::array::from_fn(|index| {
        BigRational::new(BigInt::from(index + 1), BigInt::from(index + 2))
    });
    let rhs_data = rational_matvec(&rows, &expected_solution);
    let matrix = RationalMatrix::try_from_rows(rows.clone())
        .or_abort("rational benchmark matrix construction");
    let rhs =
        RationalVector::try_new(rhs_data.clone()).or_abort("rational benchmark RHS construction");

    let reference_determinant = rational_determinant_gaussian(rows.clone());
    assert_eq!(matrix.det(), reference_determinant);
    assert_eq!(matrix.det_sign(), determinant_sign(&reference_determinant));
    let reference_solution = rational_solve_gaussian(rows, rhs_data)
        .or_abort("rational Gaussian benchmark validation solve");
    assert_eq!(reference_solution, expected_solution);
    assert_eq!(
        matrix
            .solve(&rhs)
            .or_abort("row-cleared Bareiss benchmark validation solve")
            .into_array(),
        expected_solution
    );

    ValidatedRationalInput { matrix, rhs }
}

/// Exact rational matrix-vector multiplication used only to assemble and
/// validate benchmark fixtures.
fn rational_matvec<const D: usize>(
    rows: &[[BigRational; D]; D],
    vector: &[BigRational; D],
) -> [BigRational; D] {
    std::array::from_fn(|row| {
        rows[row]
            .iter()
            .zip(vector.iter())
            .map(|(coefficient, component)| coefficient * component)
            .sum()
    })
}

/// Straightforward cubic rational Gaussian determinant reference.
#[must_use]
pub fn rational_determinant_gaussian<const D: usize>(
    mut rows: [[BigRational; D]; D],
) -> BigRational {
    let zero = BigRational::from_integer(BigInt::from(0));
    let mut determinant = BigRational::from_integer(BigInt::from(1));
    let mut odd_swaps = false;

    for pivot_col in 0..D {
        let Some(pivot_row) = (pivot_col..D).find(|&row| rows[row][pivot_col] != zero) else {
            return zero;
        };
        if pivot_row != pivot_col {
            rows.swap(pivot_col, pivot_row);
            odd_swaps = !odd_swaps;
        }

        let (pivot_rows, rows_below) = rows.split_at_mut(pivot_col + 1);
        let pivot_entries = &pivot_rows[pivot_col];
        let pivot = &pivot_entries[pivot_col];
        determinant *= pivot;
        for row_entries in rows_below {
            let factor = &row_entries[pivot_col] / pivot;
            for (entry, pivot_entry) in row_entries
                .iter_mut()
                .zip(pivot_entries.iter())
                .skip(pivot_col + 1)
            {
                *entry -= &factor * pivot_entry;
            }
            row_entries[pivot_col] = zero.clone();
        }
    }

    if odd_swaps { -determinant } else { determinant }
}

/// Straightforward cubic rational Gaussian solve reference.
#[must_use]
pub fn rational_solve_gaussian<const D: usize>(
    mut rows: [[BigRational; D]; D],
    mut rhs: [BigRational; D],
) -> Option<[BigRational; D]> {
    let zero = BigRational::from_integer(BigInt::from(0));
    for pivot_col in 0..D {
        let pivot_row = (pivot_col..D).find(|&row| rows[row][pivot_col] != zero)?;
        if pivot_row != pivot_col {
            rows.swap(pivot_col, pivot_row);
            rhs.swap(pivot_col, pivot_row);
        }

        let (pivot_rows, rows_below) = rows.split_at_mut(pivot_col + 1);
        let pivot_entries = &pivot_rows[pivot_col];
        let pivot = &pivot_entries[pivot_col];
        let (pivot_rhs_entries, rhs_below) = rhs.split_at_mut(pivot_col + 1);
        let pivot_rhs = &pivot_rhs_entries[pivot_col];
        for (row_entries, rhs_entry) in rows_below.iter_mut().zip(rhs_below) {
            let factor = &row_entries[pivot_col] / pivot;
            for (entry, pivot_entry) in row_entries
                .iter_mut()
                .zip(pivot_entries.iter())
                .skip(pivot_col + 1)
            {
                *entry -= &factor * pivot_entry;
            }
            let rhs_update = &factor * pivot_rhs;
            *rhs_entry -= rhs_update;
            row_entries[pivot_col] = zero.clone();
        }
    }

    let mut solution = std::array::from_fn(|_| zero.clone());
    for row in (0..D).rev() {
        let mut value = rhs[row].clone();
        for (coefficient, component) in rows[row].iter().zip(solution.iter()).skip(row + 1) {
            value -= coefficient * component;
        }
        solution[row] = value / &rows[row][row];
    }
    Some(solution)
}

fn determinant_sign(value: &BigRational) -> DeterminantSign {
    let zero = BigRational::from_integer(BigInt::from(0));
    match value.cmp(&zero) {
        std::cmp::Ordering::Less => DeterminantSign::Negative,
        std::cmp::Ordering::Equal => DeterminantSign::Zero,
        std::cmp::Ordering::Greater => DeterminantSign::Positive,
    }
}
