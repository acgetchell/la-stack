#![forbid(unsafe_code)]

//! Select a provisional column basis from a rectangular matrix.
//!
//! See the [exact certification workflow](https://docs.rs/la-stack/latest/la_stack/guide/exact/index.html#exact-reductions-and-provisional-column-selection)
//! for applications that need a proof about the selected minor.

use std::array::from_fn;

use la_stack::prelude::*;

fn main() -> Result<(), LaError> {
    // Column 0 is zero, column 2 is twice column 1, and column 3 adds a
    // nonzero second coordinate. The expected basis uses columns 1 and 3.
    let rows = [[0.0, 4.0, 8.0, 1.0], [0.0, 0.0, 0.0, 2.0]];
    let basis = column_basis(&rows, Tolerance::try_new(1e-12)?)?;

    // Check for full numerical row rank before extracting a square minor.
    // Even at zero tolerance the rank is provisional because elimination rounds.
    assert_eq!(basis.rank(), rows.len());
    let columns = basis.columns();
    assert_eq!(columns, &[1, 3]);
    let minor = Matrix::<2>::try_from_rows(rows.map(|row| from_fn(|c| row[columns[c]])))?;
    assert_eq!(minor.det()?.to_bits(), 8.0_f64.to_bits());
    println!("Selected columns: {columns:?}");

    // The absolute tolerance 2 rejects the second pivot (equal to 2).
    // The partial selection describes this threshold's decision. The same
    // input still has a nonzero 2×2 minor with determinant 8.
    let partial = column_basis(&rows, Tolerance::try_new(2.0)?)?;
    assert_eq!(partial.rank(), 1);
    assert_eq!(partial.columns(), &[1]);
    println!("Selected columns with tolerance 2: {:?}", partial.columns());

    Ok(())
}
