#![forbid(unsafe_code)]

//! Exact linear system solve for a near-singular 3×3 system.
//!
//! This example demonstrates `solve_exact()` and [`ExactF64Conversion`]. The exact
//! solve uses arbitrary-precision rational arithmetic to compute a provably
//! correct solution — even when the matrix is so close to singular that the f64
//! LU solve produces a wildly inaccurate result. `try_to_f64()` only succeeds
//! when every exact component is exactly representable as `f64`, while
//! `to_rounded_f64()` explicitly opts into nearest-even rounding.
//!
//! Run with: `cargo run --features exact --example exact_solve_3x3`

use la_stack::prelude::*;

fn main() -> Result<(), LaError> {
    // Near-singular 3×3 system.
    //
    // The base matrix [[1,2,3],[4,5,6],[7,8,9]] is exactly singular (rows in
    // arithmetic progression).  Perturbing entry (0,0) by 2^-50 ≈ 8.9e-16
    // makes it invertible but extremely ill-conditioned.
    let perturbation = f64::from_bits(0x3CD0_0000_0000_0000); // 2^-50
    let a = Matrix::<3>::try_from_rows([
        [1.0 + perturbation, 2.0, 3.0],
        [4.0, 5.0, 6.0],
        [7.0, 8.0, 9.0],
    ])?;

    let b = Vector::<3>::try_new([1.0, 2.0, 3.0])?;

    // f64 LU solve (using zero pivot tolerance since the matrix is nearly singular
    // and would be rejected by DEFAULT_SINGULAR_TOL).
    let lu_x = a.lu(Tolerance::try_new(0.0)?)?.solve(b)?.into_array();

    // Exact solve.
    let exact_x = a.solve_exact(b)?;
    println!("Near-singular 3×3 system (perturbation = 2^-50 ≈ {perturbation:.2e}):");
    for row in a.as_rows() {
        print!("  [");
        for (col, value) in row.iter().enumerate() {
            if col > 0 {
                print!(", ");
            }
            print!("{value:22.18}");
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
    println!(
        "f64 LU solve:       x = [{:+.6e}, {:+.6e}, {:+.6e}]",
        lu_x[0], lu_x[1], lu_x[2]
    );
    println!(
        "solve_exact():      x = [{}, {}, {}]",
        exact_x[0], exact_x[1], exact_x[2]
    );
    match exact_x.try_to_f64() {
        Ok(x) => {
            let x = x.into_array();
            println!(
                "exact try_to_f64(): x = [{:+.6e}, {:+.6e}, {:+.6e}]",
                x[0], x[1], x[2]
            );
        }
        Err(err) if err.requires_rounding() => {
            println!("exact try_to_f64(): {err}");
            let x = exact_x.to_rounded_f64()?.into_array();
            println!(
                "exact to_rounded_f64(): x = [{:+.6e}, {:+.6e}, {:+.6e}]",
                x[0], x[1], x[2]
            );
        }
        Err(err) => return Err(err),
    }
    Ok(())
}
