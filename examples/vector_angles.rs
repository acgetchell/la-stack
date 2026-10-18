#![forbid(unsafe_code)]

//! Unsigned angles for fixed vectors and borrowed coordinate slices.
//!
//! Results are in radians in `[0, π]`, with no certified error bound. The
//! assertion tolerances below apply only to these known examples.
//!
//! Run with: `cargo run --locked --example vector_angles`

use std::f64::consts::FRAC_PI_2;

use approx::assert_abs_diff_eq;

use la_stack::prelude::*;

fn main() -> Result<(), LaError> {
    // Fixed vectors validate finite coordinates at construction. Neither
    // vector needs to be normalized before computing the angle.
    let x = Vector::try_new([1.0, 0.0, 0.0])?;
    let y = Vector::try_new([0.0, 2.0, 0.0])?;
    let angle = x.angle(&y)?;
    assert_abs_diff_eq!(angle, FRAC_PI_2, epsilon = f64::EPSILON);
    println!(
        "Fixed vectors: {angle} radians ({} degrees)",
        angle.to_degrees()
    );

    // The prelude also imports VectorAngle for borrowed slices. Each operand
    // is scaled independently, even when its Euclidean norm would overflow.
    let large = [f64::MAX, f64::MAX, 0.0, 0.0];
    let small = [-f64::MIN_POSITIVE, f64::MIN_POSITIVE, 0.0, 0.0];
    let angle = large.as_slice().angle(small.as_slice())?;
    assert_abs_diff_eq!(angle, FRAC_PI_2, epsilon = f64::EPSILON);
    println!("Borrowed slices with different magnitudes: {angle} radians");

    // A representable tiny separation survives without normalizing or using
    // acos. Shared zero coordinates do not change the angle.
    let tiny = f64::from_bits(1);
    let left = [1.0, 0.0, 0.0, 0.0];
    let right = [1.0, tiny, 0.0, 0.0];
    let angle = left.as_slice().angle(right.as_slice())?;
    assert_eq!(angle, tiny);
    println!("Nearly parallel slices: {angle:e} radians");

    Ok(())
}
