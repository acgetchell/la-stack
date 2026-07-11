#![forbid(unsafe_code)]

//! Compile-time 4×4 determinant via `det_direct()`.
//!
//! Because `det_direct` is a `const fn` (Rust 1.94+), the determinant is
//! evaluated entirely at compile time — zero runtime cost.

use la_stack::prelude::*;

/// An example 4×4 matrix with small integer entries.
const MAT: Result<Matrix<4>, LaError> = Matrix::<4>::try_from_rows([
    [1.0, 2.0, 3.0, 4.0],
    [5.0, 6.0, 7.0, 8.0],
    [2.0, 6.0, 1.0, 5.0],
    [3.0, 8.0, 2.0, 9.0],
]);

/// Determinant computed at compile time.
const DET: Result<Option<f64>, LaError> = match MAT {
    Ok(matrix) => matrix.det_direct(),
    Err(err) => Err(err),
};

fn main() -> Result<(), LaError> {
    let mat = MAT?;

    println!("4×4 matrix:");
    for row in mat.as_rows() {
        print!("  [");
        for (col, value) in row.iter().enumerate() {
            if col > 0 {
                print!(", ");
            }
            print!("{value:5.1}");
        }
        println!("]");
    }
    println!();
    match DET? {
        Some(det) => println!("det (computed at compile time) = {det}"),
        None => println!("det_direct is only available for D <= 4"),
    }

    Ok(())
}
