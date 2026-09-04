#![forbid(unsafe_code)]

//! Downstream-style contract tests for the public prelude and explicit root exports.

use core::assert_matches;

use approx::assert_abs_diff_eq;

use la_stack::prelude::*;
use la_stack::{ERR_COEFF_2, ERR_COEFF_3, ERR_COEFF_4};

const _: [f64; 3] = [ERR_COEFF_2, ERR_COEFF_3, ERR_COEFF_4];

#[test]
fn common_prelude_supports_downstream_composition() -> Result<(), LaError> {
    let matrix = Matrix::<2>::identity();
    let vector = Vector::<2>::try_new([1.0, 2.0])?;
    let tolerance = Tolerance::try_new(0.0)?;
    let Some(estimate): Option<DeterminantWithErrorBound> = matrix.det_direct_with_errbound()?
    else {
        panic!("the 2x2 identity must have a certified direct determinant bound");
    };
    assert_abs_diff_eq!(estimate.determinant(), 1.0, epsilon = 0.0);
    assert!(estimate.absolute_error_bound() >= 0.0);
    let dot_estimate: ScalarWithErrorBound = vector
        .dot_with_errbound(&vector)?
        .expect("ordinary vector inputs must have a certified dot bound");
    assert_abs_diff_eq!(dot_estimate.estimate(), 5.0, epsilon = 0.0);
    assert!(dot_estimate.lower_bound() <= 5.0);
    assert!(dot_estimate.upper_bound() >= 5.0);

    let lu: Lu<2> = matrix.lu(tolerance)?;
    let ldlt: Ldlt<2> = matrix.ldlt(tolerance)?;
    let lu_solution = lu.solve(vector)?.into_array();
    let ldlt_solution = ldlt.solve(vector)?.into_array();
    for (actual, expected) in lu_solution.into_iter().zip([1.0, 2.0]) {
        assert_abs_diff_eq!(actual, expected, epsilon = 1e-12);
    }
    for (actual, expected) in ldlt_solution.into_iter().zip([1.0, 2.0]) {
        assert_abs_diff_eq!(actual, expected, epsilon = 1e-12);
    }

    assert_abs_diff_eq!(DEFAULT_SINGULAR_TOL.get(), 1e-12, epsilon = 0.0);
    assert_eq!(ArithmeticOperation::LuSolve.to_string(), "LU solve");
    assert_eq!(FactorizationKind::Lu.to_string(), "LU");

    assert_matches!(
        LaError::invalid_tolerance(-1.0),
        LaError::InvalidTolerance {
            reason: InvalidToleranceReason::Negative,
            ..
        }
    );
    assert_matches!(
        LaError::non_finite_input_vector(1),
        LaError::NonFinite {
            location: NonFiniteLocation::VectorEntry { index: 1, .. },
            origin: NonFiniteOrigin::Input,
            ..
        }
    );
    assert_matches!(
        LaError::not_positive_semidefinite_negative(0, -1.0),
        LaError::NotPositiveSemidefinite {
            violation: PositiveSemidefiniteViolation::NegativePivot { value: -1.0, .. },
            ..
        }
    );
    assert_matches!(
        LaError::singular_exact(0),
        LaError::Singular {
            reason: SingularityReason::Exact,
            ..
        }
    );
    assert!(
        LaError::unrepresentable(None, UnrepresentableReason::RequiresRounding).requires_rounding()
    );

    let dispatched = try_with_stack_matrix!(2usize, |mut dynamic| -> Result<f64, LaError> {
        dynamic.set(0, 0, 1.0)?;
        dynamic.set(1, 1, 1.0)?;
        dynamic.det()
    })?;
    assert_abs_diff_eq!(dispatched, 1.0, epsilon = 0.0);
    assert_eq!(MAX_STACK_MATRIX_DISPATCH_DIM, 7);

    let interval = Interval::try_from_subtraction(1.0, 0.1)?;
    assert!(interval.lower() <= 1.0 - 0.1);
    assert!(interval.upper() >= 1.0 - 0.1);
    let interval_matrix = IntervalMatrix::<2>::identity();
    assert_eq!(
        interval_matrix.det_sign()?,
        IntervalDeterminantSign::Positive
    );
    let interval_sign =
        try_with_interval_matrix!(MAX_INTERVAL_MATRIX_DIM, |matrix| -> Result<
            IntervalDeterminantSign,
            LaError,
        > { matrix.det_sign() },)?;
    assert_eq!(interval_sign, IntervalDeterminantSign::Zero);
    assert_matches!(
        LaError::inverted_interval(2.0, 1.0),
        LaError::InvertedInterval {
            lower: 2.0,
            upper: 1.0,
            ..
        }
    );
    Ok(())
}

#[test]
fn interval_error_categories_are_available_from_the_prelude() {
    assert_matches!(
        LaError::non_finite_input_interval_bound(IntervalBound::Lower),
        LaError::NonFinite {
            location: NonFiniteLocation::IntervalBound {
                bound: IntervalBound::Lower,
                ..
            },
            ..
        }
    );
    assert_matches!(
        LaError::non_finite_input_interval_operand(IntervalOperand::Right),
        LaError::NonFinite {
            location: NonFiniteLocation::IntervalOperand {
                operand: IntervalOperand::Right,
                ..
            },
            ..
        }
    );
    assert_matches!(
        LaError::interval_range_exhausted(ArithmeticOperation::IntervalAddition),
        LaError::IntervalRangeExhausted {
            operation: ArithmeticOperation::IntervalAddition,
            ..
        }
    );
}

#[cfg(feature = "exact")]
#[test]
fn exact_prelude_supports_downstream_composition() {
    let half = BigRational::new(BigInt::from(1), BigInt::from(2));
    let two = BigRational::from_integer(BigInt::from(2));

    assert_eq!(BigRational::from_f64(0.5).as_ref(), Some(&half));
    assert!(half.is_positive());
    assert_eq!(half.to_f64(), Some(0.5));
    assert_eq!(two.to_i64(), Some(2));
    assert_eq!(DeterminantSign::Positive.as_i8(), 1);

    let rational_matrix = RationalMatrix::<1>::try_from_rows([[half]]).unwrap();
    let rational_rhs = RationalVector::<1>::try_new([two]).unwrap();
    assert_eq!(
        rational_matrix.solve(&rational_rhs).unwrap().into_array(),
        [BigRational::from_integer(BigInt::from(4))]
    );
    let dispatched =
        try_with_rational_matrix!(MAX_RATIONAL_MATRIX_DISPATCH_DIM, |matrix| -> Result<
            DeterminantSign,
            LaError,
        > {
            Ok(matrix.det_sign())
        },)
        .unwrap();
    assert_eq!(dispatched, DeterminantSign::Zero);
}
