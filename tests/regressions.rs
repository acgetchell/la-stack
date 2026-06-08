//! Regression tests for bugs caught in public API behavior.

use la_stack::prelude::*;

#[test]
#[cfg(feature = "exact")]
fn det_exact_f64_preserves_min_positive_subnormal() -> Result<(), LaError> {
    let tiny = f64::from_bits(1);
    let m = Matrix::<1>::try_from_rows([[tiny]])?;

    assert_eq!(m.det_exact_f64().unwrap().to_bits(), tiny.to_bits());
    Ok(())
}

#[test]
#[cfg(feature = "exact")]
fn det_exact_f64_strict_vs_rounded_inexact_det() -> Result<(), LaError> {
    let m = Matrix::<2>::try_from_rows([[1.0 + f64::EPSILON, 0.0], [0.0, 1.0 - f64::EPSILON]])?;

    assert_eq!(
        m.det_exact_f64(),
        Err(LaError::Unrepresentable {
            index: None,
            reason: UnrepresentableReason::RequiresRounding,
        })
    );
    assert_eq!(
        m.det_exact_rounded_f64().unwrap().to_bits(),
        1.0f64.to_bits()
    );
    Ok(())
}

#[test]
#[cfg(feature = "exact")]
fn solve_exact_f64_strict_vs_rounded_non_dyadic() -> Result<(), LaError> {
    let a = Matrix::<1>::try_from_rows([[3.0]])?;
    let b = Vector::<1>::try_new([1.0])?;

    assert_eq!(
        a.solve_exact_f64(b),
        Err(LaError::Unrepresentable {
            index: Some(0),
            reason: UnrepresentableReason::RequiresRounding,
        })
    );
    assert_eq!(
        a.solve_exact_rounded_f64(b).unwrap().into_array()[0].to_bits(),
        (1.0f64 / 3.0).to_bits()
    );
    Ok(())
}

#[test]
#[cfg(feature = "exact")]
fn requires_rounding_error_can_fall_back_to_rounded_solve() -> Result<(), LaError> {
    let a = Matrix::<1>::try_from_rows([[3.0]])?;
    let b = Vector::<1>::try_new([1.0])?;

    let rounded = match a.solve_exact_f64(b) {
        Ok(x) => x,
        Err(err) if err.requires_rounding() => a.solve_exact_rounded_f64(b)?,
        Err(err) => return Err(err),
    };

    assert_eq!(rounded.into_array()[0].to_bits(), (1.0f64 / 3.0).to_bits());
    Ok(())
}

#[test]
fn det_direct_skips_zero_coefficient_terms_that_would_overflow() -> Result<(), LaError> {
    let d3 = Matrix::<3>::try_from_rows([
        [1.0, 0.0, 0.0],
        [1.0e300, 1.0, 1.0e300],
        [1.0e300, 0.0, 1.0e300],
    ])?;
    assert_eq!(d3.det_direct(), Ok(Some(1.0e300)));

    let d4 = Matrix::<4>::try_from_rows([
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [1.0e300, 0.0, 1.0e300, 1.0e300],
        [1.0e300, 0.0, 1.0e300, -1.0e300],
    ])?;
    assert_eq!(d4.det_direct(), Ok(Some(0.0)));

    Ok(())
}

#[test]
fn det_errbound_skips_zero_coefficient_terms_that_would_overflow() -> Result<(), LaError> {
    let d3 = Matrix::<3>::try_from_rows([
        [1.0, 0.0, 0.0],
        [1.0e300, 1.0, 1.0e300],
        [1.0e300, 0.0, 1.0e300],
    ])?;
    assert_eq!(d3.det_errbound(), Ok(Some(ERR_COEFF_3 * 1.0e300)));

    let d4 = Matrix::<4>::try_from_rows([
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [1.0e300, 0.0, 1.0e300, 1.0e300],
        [1.0e300, 0.0, 1.0e300, -1.0e300],
    ])?;
    assert_eq!(d4.det_errbound(), Ok(Some(0.0)));

    Ok(())
}
