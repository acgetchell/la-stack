//! Compile-time 4×4 determinant via `det_direct()`.
//!
//! Because `det_direct` is a `const fn` (Rust 1.94+), the determinant is
//! evaluated entirely at compile time — zero runtime cost.

use la_stack::prelude::*;

/// A 4×4 Hilbert-like matrix with exact rational entries scaled to integers.
const MAT: Matrix<4> = Matrix::<4>::from_rows([
    [1.0, 2.0, 3.0, 4.0],
    [5.0, 6.0, 7.0, 8.0],
    [2.0, 6.0, 1.0, 5.0],
    [3.0, 8.0, 2.0, 9.0],
]);

/// Determinant computed at compile time.
const DET: f64 = match MAT.det_direct() {
    Some(d) => d,
    None => panic!("det_direct only supports D <= 4"),
};

fn main() {
    println!("4×4 matrix:");
    for r in 0..4 {
        print!("  [");
        for c in 0..4 {
            if c > 0 {
                print!(", ");
            }
            print!("{:5.1}", MAT.get(r, c).unwrap());
        }
        println!("]");
    }
    println!();
    println!("det (computed at compile time) = {DET}");
}
