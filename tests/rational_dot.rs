//! Exact identities and shape contracts for fixed and runtime rational products.

#![forbid(unsafe_code)]
#![cfg(feature = "exact")]

use core::assert_matches;

use la_stack::prelude::*;
use pastey::paste;

fn ratio(n: i64, d: i64) -> BigRational {
    BigRational::new(n.into(), d.into())
}

fn reciprocal_identity<const D: usize>() {
    let left = RationalVector::<D>::try_from_fn(|i| {
        BigRational::new(BigInt::from(i + 1), BigInt::from(2 * i + 3))
    })
    .unwrap();
    let right = RationalVector::<D>::try_from_fn(|i| {
        BigRational::new(BigInt::from(2 * i + 3), BigInt::from(i + 1))
    })
    .unwrap();
    let expected = BigRational::from_integer(BigInt::from(D));
    assert_eq!(left.dot(&right), expected);
    assert_eq!(
        rational_dot(left.as_array(), right.as_array()).unwrap(),
        expected
    );
}

macro_rules! dimensions {
    ($($d:literal),+ $(,)?) => { $(paste! {
        #[test]
        fn [<reciprocal_identity_ $d>]() { reciprocal_identity::<$d>(); }
    })+ };
}

// Basis sizes 2–8, original objectives 2(D+1), and Phase I sizes 3D+4,
// including Delaunay's six-dimensional 8-row / 14- and 22-variable programs.
dimensions!(0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 13, 14, 16, 19, 22, 64);

#[test]
fn cancellation_negative_and_non_dyadic_coefficients() {
    let left = RationalVector::try_new([ratio(1, 3), ratio(-2, 5), ratio(7, 11)]).unwrap();
    let right = RationalVector::try_new([ratio(3, 7), ratio(5, 14), ratio(-11, 13)]).unwrap();
    // 1/7 - 1/7 - 7/13.
    assert_eq!(left.dot(&right), ratio(-7, 13));
    assert_eq!(
        rational_dot(left.as_array(), right.as_array()).unwrap(),
        ratio(-7, 13)
    );
}

#[test]
fn large_coefficients_cancel_without_a_binary64_intermediate() {
    let huge = (BigInt::from(1) << 4096usize) + 1;
    let denominator = (BigInt::from(1) << 8192usize) + 3;
    let a = BigRational::new(huge, denominator);
    let b = a.recip();
    let left = RationalVector::try_new([a.clone(), -a, ratio(1, 3)]).unwrap();
    let right = RationalVector::try_new([b.clone(), b, ratio(1, 7)]).unwrap();
    assert_eq!(left.dot(&right), ratio(1, 21));
    assert_eq!(
        rational_dot(left.as_array(), right.as_array()).unwrap(),
        ratio(1, 21)
    );
}

#[test]
fn runtime_adapter_accepts_raw_quotients_and_preserves_inputs() {
    let left = [
        BigRational::new_raw(2.into(), (-6).into()),
        BigRational::new_raw(0.into(), (-7).into()),
    ];
    let right = [BigRational::new_raw((-9).into(), 3.into()), ratio(3, 11)];
    let result = rational_dot(&left, &right).unwrap();
    assert_eq!(result.numer(), &BigInt::from(1));
    assert_eq!(result.denom(), &BigInt::from(1));
    assert_eq!(left[0].numer(), &BigInt::from(2));
    assert_eq!(left[0].denom(), &BigInt::from(-6));
    let zero = rational_dot(&left[1..], &right[1..]).unwrap();
    assert_eq!(zero.numer(), &BigInt::from(0));
    assert_eq!(zero.denom(), &BigInt::from(1));
}

#[test]
fn shape_is_checked_before_denominators_in_both_directions() {
    let invalid = [BigRational::new_raw(1.into(), 0.into())];
    for (left, right) in [(invalid.as_slice(), &[][..]), (&[][..], invalid.as_slice())] {
        assert_matches!(rational_dot(left, right),
            Err(LaError::DimensionMismatch { left: l, right: r, .. })
                if l == left.len() && r == right.len());
    }
}

#[test]
fn denominator_errors_retain_operand_and_index_even_for_zero_products() {
    let valid = [ratio(0, 1), ratio(0, 1)];
    let invalid = [ratio(0, 1), BigRational::new_raw(0.into(), 0.into())];
    for (left, right, expected_operand) in [
        (&invalid, &valid, VectorOperand::Left),
        (&valid, &invalid, VectorOperand::Right),
        (&invalid, &invalid, VectorOperand::Left),
    ] {
        assert_matches!(rational_dot(left, right), Err(LaError::NonFinite {
            origin: NonFiniteOrigin::Input,
            location: NonFiniteLocation::VectorOperandEntry { operand, index: 1, .. }, ..
        }) if operand == expected_operand);
    }
}
