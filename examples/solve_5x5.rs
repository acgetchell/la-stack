//! Solve a 5×5 linear system via LU factorization (with pivoting).

use la_stack::prelude::*;

fn main() -> Result<(), LaError> {
    // This system requires pivoting (a[0][0] = 0), so it's a good LU demo.
    // A = J - I: zeros on diagonal, ones elsewhere.
    let a = Matrix::<5>::from_rows([
        [0.0, 1.0, 1.0, 1.0, 1.0],
        [1.0, 0.0, 1.0, 1.0, 1.0],
        [1.0, 1.0, 0.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0, 1.0, 0.0],
    ]);

    // Choose x = [1, 2, 3, 4, 5]. Then b = A x = [14, 13, 12, 11, 10].
    let b = Vector::<5>::new([14.0, 13.0, 12.0, 11.0, 10.0]);

    let lu = a.lu(DEFAULT_PIVOT_TOL)?;
    let x = lu.solve_vec(b)?.into_array();

    println!("x = {x:?}");
    Ok(())
}
