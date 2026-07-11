#![forbid(unsafe_code)]

//! Regression coverage for range-scaled public factor determinants.

use la_stack::{LaError, Matrix, Tolerance};

const SUBNORMAL_ROUNDING_LEFT: f64 = f64::from_bits(0x3cb2_e219_27ac_435a);
const SUBNORMAL_ROUNDING_RIGHT: f64 = f64::from_bits(0x0014_55e5_f80b_50eb);
const TWO_NEG_800: f64 = f64::from_bits(223_u64 << 52);
const TWO_POS_800: f64 = f64::from_bits(1823_u64 << 52);
const LEAST_SUBNORMAL_BITS: u64 = 1;

/// Build a finite diagonal matrix from the supplied entries.
fn diagonal_matrix<const D: usize>(diagonal: [f64; D]) -> Result<Matrix<D>, LaError> {
    let mut rows = [[0.0; D]; D];
    for (index, value) in diagonal.into_iter().enumerate() {
        rows[index][index] = value;
    }
    Matrix::try_from_rows(rows)
}

#[test]
fn public_lu_and_ldlt_determinants_round_final_subnormal_once() -> Result<(), LaError> {
    let expected = SUBNORMAL_ROUNDING_LEFT * SUBNORMAL_ROUNDING_RIGHT;
    assert_eq!(expected.to_bits(), LEAST_SUBNORMAL_BITS);

    let matrix = diagonal_matrix([SUBNORMAL_ROUNDING_LEFT, SUBNORMAL_ROUNDING_RIGHT])?;
    let zero_tolerance = Tolerance::try_new(0.0)?;
    let lu_det = matrix.lu(zero_tolerance)?.det()?;
    let ldlt_det = matrix.ldlt(zero_tolerance)?.det()?;

    assert_eq!(lu_det.to_bits(), LEAST_SUBNORMAL_BITS);
    assert_eq!(ldlt_det.to_bits(), LEAST_SUBNORMAL_BITS);
    Ok(())
}

#[test]
fn public_factor_determinants_round_subnormal_after_earlier_range_loss() -> Result<(), LaError> {
    let expected = SUBNORMAL_ROUNDING_LEFT * SUBNORMAL_ROUNDING_RIGHT;
    assert_eq!(expected.to_bits(), LEAST_SUBNORMAL_BITS);

    let matrix = diagonal_matrix([
        TWO_NEG_800,
        TWO_NEG_800,
        TWO_POS_800,
        TWO_POS_800,
        SUBNORMAL_ROUNDING_LEFT,
        SUBNORMAL_ROUNDING_RIGHT,
    ])?;
    let zero_tolerance = Tolerance::try_new(0.0)?;
    let lu_det = matrix.lu(zero_tolerance)?.det()?;
    let ldlt_det = matrix.ldlt(zero_tolerance)?.det()?;

    assert_eq!(lu_det.to_bits(), LEAST_SUBNORMAL_BITS);
    assert_eq!(ldlt_det.to_bits(), LEAST_SUBNORMAL_BITS);
    Ok(())
}
