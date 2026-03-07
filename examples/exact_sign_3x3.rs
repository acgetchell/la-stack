//! Exact determinant sign for a near-singular 3×3 matrix.
//!
//! This example demonstrates `det_sign_exact()`, which uses adaptive-precision
//! arithmetic to return the provably correct sign of the determinant — even when
//! the matrix is so close to singular that f64 rounding could give the wrong answer.
//!
//! Run with: `cargo run --features exact --example exact_sign_3x3`

use la_stack::prelude::*;

fn main() {
    // Base matrix: rows in arithmetic progression → exactly singular (det = 0).
    //   [[1, 2, 3],
    //    [4, 5, 6],
    //    [7, 8, 9]]
    //
    // Perturb entry (0,0) by 2^-50 ≈ 8.9e-16.
    // Exact det = 2^-50 × cofactor(0,0) = 2^-50 × (5×9 − 6×8) = −3 × 2^-50 < 0.
    let perturbation = f64::from_bits(0x3CD0_0000_0000_0000); // 2^-50
    let m = Matrix::<3>::from_rows([
        [1.0 + perturbation, 2.0, 3.0],
        [4.0, 5.0, 6.0],
        [7.0, 8.0, 9.0],
    ]);

    let sign = m.det_sign_exact();
    let det_f64 = m.det(DEFAULT_PIVOT_TOL).unwrap();

    println!("Near-singular 3×3 matrix (perturbation = 2^-50 ≈ {perturbation:.2e}):");
    for r in 0..3 {
        print!("  [");
        for c in 0..3 {
            if c > 0 {
                print!(", ");
            }
            print!("{:22.18}", m.get(r, c).unwrap());
        }
        println!("]");
    }
    println!();
    println!("f64 det()          = {det_f64:+.6e}");
    println!("det_sign_exact()   = {sign}");
    println!();
    println!("The exact sign is −1 (negative), matching the analytical result.");
}
