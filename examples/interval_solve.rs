#![forbid(unsafe_code)]

//! Divide by a zero-free interval and enclose the solutions of an interval system.

use la_stack::prelude::*;

fn main() -> Result<(), LaError> {
    let divisor = Interval::try_new(2.0, 4.0)?;
    assert_eq!(divisor.try_recip()?, Interval::try_new(0.25, 0.5)?);
    let quotient = Interval::try_new(4.0, 8.0)?.try_div(&divisor)?;
    assert_eq!(quotient, Interval::try_new(1.0, 4.0)?);

    // For every a ∈ [2, 4] and b ∈ [6, 10], the equations a*x + y = b
    // and 2*y = 4 have a unique solution with y = 2 and x = (b - 2)/a.
    let two = Interval::point(2.0)?;
    let matrix = IntervalMatrix::from_rows([[divisor, Interval::ONE], [Interval::ZERO, two]]);
    let rhs = [Interval::try_new(6.0, 10.0)?, Interval::point(4.0)?];
    let solution = matrix.solve(&rhs)?;
    assert_eq!(solution, [Interval::try_new(1.0, 4.0)?, two]);
    println!("Solution enclosure for every admitted system: {solution:?}");

    // Direct division rejects a supplied divisor containing zero. A solve
    // instead reports when elimination cannot certify a pivot.
    let crosses_zero = Interval::try_new(-1.0, 1.0)?;
    assert!(matches!(
        Interval::ONE.try_div(&crosses_zero),
        Err(LaError::IntervalZeroDivisor {
            lower: -1.0,
            upper: 1.0,
            ..
        })
    ));

    // The exact determinant is 3 × 2^-52. Rounding in elimination can still
    // leave the final pivot enclosing zero.
    let regular = IntervalMatrix::try_from_point_rows([[3.0, 3.0], [1.0, 1.0_f64.next_up()]])?;
    match regular.solve(&[Interval::ONE; 2]) {
        Err(LaError::IntervalSolveInconclusive { pivot_col, .. }) => {
            assert_eq!(pivot_col, 1);
            println!("Pivot {pivot_col} is inconclusive. Use a certified fallback.");
        }
        Err(error) => return Err(error),
        Ok(bounds) => panic!("expected an inconclusive pivot, got {bounds:?}"),
    }

    // Here the unique real solution is 2 × f64::MAX, outside the finite
    // endpoint range. A certified fallback needs a wider representation.
    let small_pivot = IntervalMatrix::try_from_point_rows([[0.5]])?;
    match small_pivot.solve(&[Interval::point(f64::MAX)?]) {
        Err(LaError::IntervalRangeExhausted { operation, .. }) => {
            assert_eq!(operation, ArithmeticOperation::IntervalSolve);
            println!("Finite endpoint range exhausted. Use a certified fallback.");
        }
        Err(error) => return Err(error),
        Ok(bounds) => panic!("expected finite endpoint range exhaustion, got {bounds:?}"),
    }

    Ok(())
}
