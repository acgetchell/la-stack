//! Exact linear system solve for a near-singular 3×3 system.
//!
//! This example demonstrates `solve_exact()` and `solve_exact_f64()`, which use
//! arbitrary-precision rational arithmetic to compute a provably correct
//! solution — even when the matrix is so close to singular that the f64 LU
//! solve produces a wildly inaccurate result.
//!
//! Run with: `cargo run --features exact --example exact_solve_3x3`

use la_stack::prelude::*;

fn main() {
    // Near-singular 3×3 system.
    //
    // The base matrix [[1,2,3],[4,5,6],[7,8,9]] is exactly singular (rows in
    // arithmetic progression).  Perturbing entry (0,0) by 2^-50 ≈ 8.9e-16
    // makes it invertible but extremely ill-conditioned.
    let perturbation = f64::from_bits(0x3CD0_0000_0000_0000); // 2^-50
    let a = Matrix::<3>::from_rows([
        [1.0 + perturbation, 2.0, 3.0],
        [4.0, 5.0, 6.0],
        [7.0, 8.0, 9.0],
    ]);

    let b = Vector::<3>::new([1.0, 2.0, 3.0]);

    // f64 LU solve (using zero pivot tolerance since the matrix is nearly singular
    // and would be rejected by DEFAULT_PIVOT_TOL).
    let lu_x = a.lu(0.0).unwrap().solve_vec(b).unwrap().into_array();

    // Exact solve.
    let exact_x = a.solve_exact(b).unwrap();
    let exact_x_f64 = a.solve_exact_f64(b).unwrap().into_array();

    println!("Near-singular 3×3 system (perturbation = 2^-50 ≈ {perturbation:.2e}):");
    for r in 0..3 {
        print!("  [");
        for c in 0..3 {
            if c > 0 {
                print!(", ");
            }
            print!("{:22.18}", a.get(r, c).unwrap());
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
    println!(
        "solve_exact_f64():  x = [{:+.6e}, {:+.6e}, {:+.6e}]",
        exact_x_f64[0], exact_x_f64[1], exact_x_f64[2]
    );
}
