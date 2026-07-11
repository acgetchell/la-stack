#![forbid(unsafe_code)]

//! Solve a 3×3 symmetric positive definite system via LDLT factorization.
//!
//! LDLT is the natural choice for SPD matrices (e.g. Gram matrices, covariance
//! matrices, stiffness matrices).  It avoids the row swaps of LU and exploits
//! symmetry, using roughly half the work.
//!
//! Run with: `cargo run --example ldlt_solve_3x3`

use la_stack::prelude::*;

fn main() -> Result<(), LaError> {
    // Symmetric positive definite 3×3 matrix (classic SPD tridiagonal).
    let a = Matrix::<3>::try_from_rows([[4.0, -1.0, 0.0], [-1.0, 4.0, -1.0], [0.0, -1.0, 4.0]])?;

    // Choose x = [1, 2, 3].  Then b = A x = [2, 4, 10].
    let b = Vector::<3>::try_new([2.0, 4.0, 10.0])?;

    let ldlt = a.ldlt(DEFAULT_SINGULAR_TOL)?;
    let x = ldlt.solve(b)?.into_array();
    let det = ldlt.det()?;

    println!("A (3×3 SPD tridiagonal):");
    for row in a.as_rows() {
        print!("  [");
        for (col, value) in row.iter().enumerate() {
            if col > 0 {
                print!(", ");
            }
            print!("{value:5.1}");
        }
        println!("]");
    }
    println!(
        "b = [{}, {}, {}]",
        b.as_array()[0],
        b.as_array()[1],
        b.as_array()[2]
    );
    println!();
    println!("x   = [{:.6}, {:.6}, {:.6}]", x[0], x[1], x[2]);
    println!("det = {det}");

    Ok(())
}
