#![forbid(unsafe_code)]

//! Structure-sensitive exact solve fixtures and the captured Delaunay control.

use std::array::from_fn;

use la_stack::{BigInt, BigRational, LaError, RationalMatrix, RationalVector};
use num_traits::Zero;

use super::bench_utils::OrAbort;

/// Matrix structures covered by the solve benchmark.
#[derive(Clone, Copy, Debug)]
pub enum Structure {
    /// Non-zero diagonal only.
    Diagonal,
    /// Dense upper triangle.
    Upper,
    /// Dense lower triangle.
    Lower,
    /// Upper triangle with rows rotated left once (issue #246).
    PermutedUpper,
    /// Lower triangle with rows rotated left once.
    PermutedLower,
    /// Strictly diagonally dominant tridiagonal matrix.
    Sparse,
    /// Hilbert matrix, or a dense dyadic diagonally dominant control.
    Dense,
}

impl Structure {
    /// Stable registration order.
    pub const ALL: [Self; 7] = [
        Self::Diagonal,
        Self::Upper,
        Self::Lower,
        Self::PermutedUpper,
        Self::PermutedLower,
        Self::Sparse,
        Self::Dense,
    ];

    /// Stable benchmark name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Diagonal => "diagonal",
            Self::Upper => "upper",
            Self::Lower => "lower",
            Self::PermutedUpper => "permuted_upper",
            Self::PermutedLower => "permuted_lower",
            Self::Sparse => "sparse",
            Self::Dense => "dense",
        }
    }
}

/// Inputs validated against a manufactured solution and independent elimination.
#[must_use]
pub struct ValidatedSolveInput<const D: usize> {
    matrix: RationalMatrix<D>,
    rhs: RationalVector<D>,
    rows: Vec<Vec<BigRational>>,
    values: Vec<BigRational>,
}

impl<const D: usize> ValidatedSolveInput<D> {
    /// Borrow the accepted matrix.
    #[must_use]
    pub const fn matrix(&self) -> &RationalMatrix<D> {
        &self.matrix
    }

    /// Borrow the accepted right-hand side.
    #[must_use]
    pub const fn rhs(&self) -> &RationalVector<D> {
        &self.rhs
    }

    /// Run the captured consuming solver, including its required input copies.
    #[must_use]
    pub fn legacy(&self) -> Option<Vec<BigRational>> {
        legacy_solve(self.rows.clone(), self.values.clone())
    }

    /// Construct fixed-size inputs using the downstream setter-based adapter.
    ///
    /// # Errors
    /// Propagates public rational construction errors.
    pub fn construct(&self) -> Result<(RationalMatrix<D>, RationalVector<D>), LaError> {
        let mut matrix = RationalMatrix::zero();
        for (row, entries) in self.rows.iter().enumerate() {
            for (col, value) in entries.iter().enumerate() {
                matrix.set(row, col, value.clone())?;
            }
        }
        let rhs = RationalVector::try_from_fn(|index| self.values[index].clone())?;
        Ok((matrix, rhs))
    }

    /// Full downstream runtime-dispatch adapter, including output collection.
    ///
    /// # Errors
    /// Propagates construction, dimension, and exact singularity errors.
    pub fn adapter(&self) -> Result<Vec<BigRational>, LaError> {
        la_stack::try_with_rational_matrix!(self.values.len(), |mut matrix| -> Result<
            Vec<BigRational>,
            LaError,
        > {
            for (row, entries) in self.rows.iter().enumerate() {
                for (col, value) in entries.iter().enumerate() {
                    matrix.set(row, col, value.clone())?;
                }
            }
            let rhs = RationalVector::try_from_fn(|index| self.values[index].clone())?;
            Ok(matrix.solve(&rhs)?.into_array().into_iter().collect())
        })
    }
}

/// Generate exact coefficients, with no intermediate binary64 reconstruction.
///
/// # Panics
/// Panics if any independent fixture check or public construction fails.
pub fn solve_input<const D: usize>(structure: Structure, dyadic: bool) -> ValidatedSolveInput<D> {
    let mut rows = from_fn(|row| {
        from_fn(|col| {
            let present = match structure {
                Structure::Diagonal => row == col,
                Structure::Upper | Structure::PermutedUpper => col >= row,
                Structure::Lower | Structure::PermutedLower => col <= row,
                Structure::Sparse => row.abs_diff(col) <= 1,
                Structure::Dense => true,
            };
            if !present {
                return BigRational::zero();
            }
            if matches!(structure, Structure::Dense) && !dyadic {
                return BigRational::new(BigInt::from(1), BigInt::from(row + col + 1));
            }
            let numerator =
                if row == col && matches!(structure, Structure::Sparse | Structure::Dense) {
                    4 * D * D + 1
                } else {
                    row + col + 1
                };
            let denominator = if dyadic {
                1_usize << (row % 4 + 1)
            } else {
                3 + 2 * row
            };
            BigRational::new(BigInt::from(numerator), BigInt::from(denominator))
        })
    });
    if D > 0
        && matches!(
            structure,
            Structure::PermutedUpper | Structure::PermutedLower | Structure::Dense
        )
    {
        rows.rotate_left(1);
    }
    let expected: [BigRational; D] = from_fn(|col| {
        BigRational::new(
            BigInt::from(col + 1),
            BigInt::from(if dyadic { 8 } else { 11 }),
        )
    });
    let rhs = from_fn(|row| rows[row].iter().zip(&expected).map(|(a, x)| a * x).sum());
    let input = ValidatedSolveInput {
        matrix: RationalMatrix::try_from_rows(rows.clone()).or_abort("solve fixture matrix"),
        rhs: RationalVector::try_new(rhs).or_abort("solve fixture RHS"),
        rows: rows.into_iter().map(Vec::from).collect(),
        values: Vec::new(),
    };
    let input = ValidatedSolveInput {
        values: input.rhs.as_array().to_vec(),
        ..input
    };
    assert_eq!(input.legacy().or_abort("legacy fixture solve"), expected);
    assert_eq!(
        input
            .matrix
            .solve(&input.rhs)
            .or_abort("prepared fixture solve")
            .as_array(),
        &expected
    );
    assert_eq!(input.adapter().or_abort("adapter fixture solve"), expected);
    let (matrix, rhs) = input.construct().or_abort("construction fixture check");
    assert_eq!(matrix, input.matrix);
    assert_eq!(rhs, input.rhs);
    input
}

/// Captured Delaunay Gaussian solver with zero-coefficient skips.
///
/// From commit 2b7017bb3981babed2aab64826282bb06ac0e312,
/// `src/geometry/matrix.rs`, lines 215–256. Arithmetic, cloning, and zero skips
/// are retained; only the function name and visibility differ.
#[must_use]
#[expect(
    clippy::needless_range_loop,
    reason = "preserve the captured legacy control"
)]
pub fn legacy_solve(
    mut matrix: Vec<Vec<BigRational>>,
    mut rhs: Vec<BigRational>,
) -> Option<Vec<BigRational>> {
    let dimension = rhs.len();
    if matrix.len() != dimension || matrix.iter().any(|row| row.len() != dimension) {
        return None;
    }
    for pivot_col in 0..dimension {
        let pivot_row = (pivot_col..dimension).find(|&row| !matrix[row][pivot_col].is_zero())?;
        if pivot_row != pivot_col {
            matrix.swap(pivot_col, pivot_row);
            rhs.swap(pivot_col, pivot_row);
        }
        let pivot = matrix[pivot_col][pivot_col].clone();
        for row in pivot_col + 1..dimension {
            if matrix[row][pivot_col].is_zero() {
                continue;
            }
            let factor = matrix[row][pivot_col].clone() / pivot.clone();
            matrix[row][pivot_col] = BigRational::from_integer(0.into());
            for column in pivot_col + 1..dimension {
                matrix[row][column] = matrix[row][column].clone()
                    - factor.clone() * matrix[pivot_col][column].clone();
            }
            rhs[row] = rhs[row].clone() - factor * rhs[pivot_col].clone();
        }
    }
    let zero = BigRational::from_integer(0.into());
    let mut solution = vec![zero; dimension];
    for row in (0..dimension).rev() {
        let mut value = rhs[row].clone();
        for column in row + 1..dimension {
            value -= matrix[row][column].clone() * solution[column].clone();
        }
        solution[row] = value / matrix[row][row].clone();
    }
    Some(solution)
}
