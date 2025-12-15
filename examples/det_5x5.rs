//! Compute the determinant of a 5×5 matrix via explicit LU factorization.

use la_stack::prelude::*;

fn main() -> Result<(), LaError> {
    // 5×5 matrix with zeros on diagonal and ones elsewhere (J - I).
    // det(J - I) = (D - 1) * (-1)^(D-1) = 4 for D=5.
    let a = Matrix::<5>::from_rows([
        [0.0, 1.0, 1.0, 1.0, 1.0],
        [1.0, 0.0, 1.0, 1.0, 1.0],
        [1.0, 1.0, 0.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0, 1.0, 0.0],
    ]);

    // Compute via explicit LU factorization.
    let lu = a.lu(DEFAULT_PIVOT_TOL)?;
    let det = lu.det();

    println!("det = {det}");
    Ok(())
}
