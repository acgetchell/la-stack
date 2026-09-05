#![forbid(unsafe_code)]

//! Independent integer matrix-product oracles for Gram construction.

use core::array::from_fn;

use approx::assert_abs_diff_eq;
use pastey::paste;
use proptest::prelude::*;

use la_stack::prelude::*;

#[path = "common/proptest_config.rs"]
mod proptest_config;
use proptest_config::with_default_cases;

fn check_product<const M: usize, const N: usize>(entries: &[i16]) {
    let integers: [[i16; N]; M] = from_fn(|i| from_fn(|k| entries[i * N + k]));
    let vectors = integers.map(|row| Vector::try_new(row.map(f64::from)).unwrap());
    let result = gram_matrix(&vectors).unwrap();
    for i in 0..M {
        for j in 0..M {
            // Separate integer multiplication and summation, not the FMA kernel.
            let expected: i32 = (0..N)
                .map(|k| i32::from(integers[i][k]) * i32::from(integers[j][k]))
                .sum();
            assert_abs_diff_eq!(result.as_rows()[i][j], f64::from(expected), epsilon = 0.0);
            assert_eq!(
                result.as_rows()[i][j].to_bits(),
                result.as_rows()[j][i].to_bits()
            );
        }
    }
}

fn check_orthogonal_and_dependent<const D: usize>() {
    let vectors: [Vector<D>; D] =
        from_fn(|i| Vector::try_new(from_fn(|j| if i == j { 2.0 } else { 0.0 })).unwrap());
    let gram = gram_matrix(&vectors).unwrap();
    for i in 0..D {
        for j in 0..D {
            assert_abs_diff_eq!(
                gram.as_rows()[i][j],
                if i == j { 4.0 } else { 0.0 },
                epsilon = 0.0
            );
        }
    }
    let expected = (0..D).fold(1.0, |product, _| product * 4.0);
    assert_abs_diff_eq!(
        gram.ldlt(Tolerance::try_new(0.0).unwrap())
            .unwrap()
            .det()
            .unwrap(),
        expected,
        epsilon = 0.0
    );
    let dependent = [vectors[0]; D];
    let dependent_gram = gram_matrix(&dependent).unwrap();
    assert_eq!(dependent_gram.as_rows(), &[[4.0; D]; D]);
    assert_eq!(
        dependent_gram
            .ldlt(Tolerance::try_new(0.0).unwrap())
            .unwrap_err(),
        LaError::singular_numerical(1, FactorizationKind::Ldlt, 0.0, 0.0)
    );
}

// A separately rounded product or reversed reduction loses the exact residual.
fn check_fused_reduction<const D: usize>() {
    let delta = 2.0_f64.powi(-27);
    let left: Vector<D> = Vector::try_new(from_fn(|k| match k {
        0 => 1.0,
        1 => 1.0 + delta,
        _ => 0.0,
    }))
    .unwrap();
    let right = Vector::try_new(from_fn(|k| match k {
        0 => -1.0,
        1 => 1.0 - delta,
        _ => 0.0,
    }))
    .unwrap();
    let gram = gram_matrix(&[left, right]).unwrap();
    // -1 + (1 + 2^-27)(1 - 2^-27) = -2^-54, exactly representable.
    let expected = (-2.0_f64.powi(-54)).to_bits();
    assert_eq!(gram.as_rows()[0][1].to_bits(), expected);
    assert_eq!(gram.as_rows()[1][0].to_bits(), expected);
}

macro_rules! gram_cases {
    ($d:literal) => {
        paste! {
            proptest! {
                #![proptest_config(with_default_cases(64))]
                #[test]
                fn [<gram_product_ $d d>](entries in prop::collection::vec(-100i16..=100, 64)) {
                    check_product::<$d, $d>(&entries);
                    if $d != 8 {
                        check_product::<$d, 8>(&entries);
                        check_product::<8, $d>(&entries);
                    }
                }
            }

            #[test]
            fn [<gram_orthogonal_and_dependent_ $d d>]() {
                check_orthogonal_and_dependent::<$d>();
            }

            #[test]
            fn [<gram_fused_reduction_ $d d>]() {
                check_fused_reduction::<$d>();
            }
        }
    };
}

gram_cases!(2);
gram_cases!(3);
gram_cases!(4);
gram_cases!(5);
gram_cases!(8);

#[test]
fn empty_dimensions_and_const_evaluation() {
    const EMPTY: Result<Matrix<0>, LaError> = gram_matrix::<0, 8>(&[]);
    const ZERO: Result<Matrix<3>, LaError> = gram_matrix(&[Vector::<0>::zero(); 3]);
    assert_eq!(EMPTY.unwrap().as_rows(), &[] as &[[f64; 0]; 0]);
    assert_eq!(ZERO.unwrap().as_rows(), &[[0.0; 3]; 3]);
    assert_eq!(
        gram_matrix::<0, 0>(&[]).unwrap().as_rows(),
        &[] as &[[f64; 0]; 0]
    );
}

#[test]
fn nonempty_const_evaluation() {
    const VECTOR: Vector<2> = match Vector::try_new([3.0, 4.0]) {
        Ok(vector) => vector,
        Err(_) => panic!("finite literal input"),
    };
    const GRAM: Result<Matrix<1>, LaError> = gram_matrix(&[VECTOR]);
    assert_eq!(GRAM.unwrap().as_rows(), &[[25.0]]);
}

#[test]
fn underflow_preserves_mirrored_negative_zero() {
    let vectors = [
        Vector::try_new([2.0_f64.powi(-537)]).unwrap(),
        Vector::try_new([-2.0_f64.powi(-538)]).unwrap(),
    ];
    let gram = gram_matrix(&vectors).unwrap();
    assert_eq!(gram.as_rows()[0][0].to_bits(), 1); // Smallest positive subnormal.
    assert_eq!(gram.as_rows()[0][1].to_bits(), (-0.0_f64).to_bits());
    assert_eq!(gram.as_rows()[1][0].to_bits(), (-0.0_f64).to_bits());
    assert_eq!(gram.as_rows()[1][1].to_bits(), 0);
}

#[test]
fn overflow_checks_off_diagonal_before_later_diagonal() {
    let vectors = [
        Vector::try_new([0.0, 2.0]).unwrap(),
        Vector::try_new([f64::MAX, f64::MAX]).unwrap(),
    ];
    // (0,0) succeeds; (0,1) fails at coordinate 1. Evaluating (1,1) first
    // would instead report coordinate 0, violating upper-triangle row order.
    assert_eq!(
        gram_matrix(&vectors).unwrap_err(),
        LaError::non_finite_computation_step(ArithmeticOperation::VectorDotProduct, 1)
    );

    let later_diagonal = [
        Vector::try_new([1.0, 0.0]).unwrap(),
        Vector::try_new([0.0, f64::MAX]).unwrap(),
    ];
    assert_eq!(
        gram_matrix(&later_diagonal).unwrap_err(),
        LaError::non_finite_computation_step(ArithmeticOperation::VectorDotProduct, 1)
    );
}

#[test]
fn mixed_scale_and_signed_zero() {
    let vectors = [
        Vector::try_new([2.0_f64.powi(500), 2.0_f64.powi(-500), -0.0]).unwrap(),
        Vector::try_new([2.0_f64.powi(-500), -2.0_f64.powi(500), 0.0]).unwrap(),
    ];
    let gram = gram_matrix(&vectors).unwrap();
    assert_eq!(
        gram.as_rows(),
        &[[2.0_f64.powi(1000), 0.0], [0.0, 2.0_f64.powi(1000)]]
    );
    assert_eq!(
        gram.as_rows()[0][1].to_bits(),
        gram.as_rows()[1][0].to_bits()
    );
}

#[test]
fn overflow_preserves_dot_diagnostics() {
    for (data, index) in [([f64::MAX, 0.0], 0), ([1.0e154, 1.0e154], 1)] {
        let vectors = [Vector::try_new(data).unwrap()];
        assert!(matches!(
            gram_matrix(&vectors),
            Err(LaError::NonFinite {
                origin: NonFiniteOrigin::Computation { operation: ArithmeticOperation::VectorDotProduct, .. },
                location: NonFiniteLocation::Step { index: actual, .. },
                ..
            }) if actual == index
        ));
    }
}
