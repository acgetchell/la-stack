#![forbid(unsafe_code)]

//! Solve a rational 5×5 system that becomes singular after conversion to binary64.

use la_stack::prelude::*;

fn main() -> Result<(), LaError> {
    let epsilon = BigRational::new(1.into(), (1_u64 << 60).into());
    let one = BigRational::from_integer(1.into());
    let zero = BigRational::from_integer(0.into());

    let matrix = RationalMatrix::<5>::try_from_fn(|row, col| match (row, col) {
        (0, 0 | 1) | (1, 0) => one.clone(),
        (1, 1) => &one + &epsilon,
        _ if row == col => one.clone(),
        _ => zero.clone(),
    })?;
    let rhs = RationalVector::try_new([
        zero,
        -&epsilon,
        BigRational::from_integer(2.into()),
        BigRational::from_integer(3.into()),
        BigRational::from_integer(4.into()),
    ])?;

    let exact_solution = matrix.solve(&rhs)?;
    println!("exact determinant: {}", matrix.det());
    println!("exact determinant sign: {:?}", matrix.det_sign());
    println!("exact solution: {:?}", exact_solution.as_array());

    let epsilon_f64 = epsilon.try_to_f64()?;
    assert_eq!((1.0 + epsilon_f64).to_bits(), 1.0_f64.to_bits());
    let f64_matrix = Matrix::<5>::try_from_rows([
        [1.0, 1.0, 0.0, 0.0, 0.0],
        [1.0, 1.0 + epsilon_f64, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 1.0],
    ])?;
    let f64_rhs = Vector::<5>::try_new([0.0, -epsilon_f64, 2.0, 3.0, 4.0])?;
    let f64_solve = f64_matrix
        .lu(DEFAULT_SINGULAR_TOL)
        .and_then(|lu| lu.solve(f64_rhs));
    assert!(matches!(&f64_solve, Err(LaError::Singular { .. })));

    println!("1.0 + 2^-60 supplied as f64: {}", 1.0 + epsilon_f64);
    println!("solve from f64 inputs: {f64_solve:?}");

    Ok(())
}
