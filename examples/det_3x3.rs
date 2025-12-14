//! Compute the determinant of a 3×3 matrix via explicit LU factorization.

use la_stack::{DEFAULT_PIVOT_TOL, LaError, Matrix};

fn main() -> Result<(), LaError> {
    let a = Matrix::<3>::from_rows([[1.0, 2.0, 3.0], [0.0, 4.0, 5.0], [1.0, 0.0, 6.0]]);

    // Compute via explicit LU factorization.
    let lu = a.lu(DEFAULT_PIVOT_TOL)?;
    let det = lu.det();

    println!("det = {det}");
    Ok(())
}
