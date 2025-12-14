//! Solve a 3×3 linear system via LU factorization (with pivoting).

use la_stack::{DEFAULT_PIVOT_TOL, LaError, Matrix, Vector};

fn main() -> Result<(), LaError> {
    // This system requires pivoting (a[0][0] = 0), so it's a good LU demo.
    let a = Matrix::<3>::from_rows([[0.0, 1.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 0.0]]);
    let b = Vector::<3>::new([5.0, 4.0, 3.0]);

    let lu = a.lu(DEFAULT_PIVOT_TOL)?;
    let x = lu.solve_vec(b)?.into_array();

    println!("x = [{:.6}, {:.6}, {:.6}]", x[0], x[1], x[2]);
    Ok(())
}
