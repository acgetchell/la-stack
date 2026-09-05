//! Allocation evidence for canonical rational row clearing (issue #233).
//!
//! The counting allocator is linked only into this test executable; Criterion
//! timing executables retain their normal allocator. All fixture construction,
//! row scales, and exact oracle checks precede the measured closures.

#![cfg(all(feature = "bench", feature = "exact"))]
#![forbid(unsafe_code)]

#[path = "../benches/common/bench_utils.rs"]
mod bench_utils;
#[path = "../benches/common/rational.rs"]
pub mod rational_bench;

use std::{array::from_fn, hint::black_box};

use allocation_counter::{AllocationInfo, measure};
use la_stack::{BigInt, BigRational};
use pastey::paste;
use rational_bench::{RationalInputKind, rational_input};

/// Reproduce the positive-canonical branch before a62e72a removed its clones.
fn scaled_cloned(value: &BigRational, scale: &BigInt) -> BigInt {
    let denominator = value.denom().clone();
    let multiplier = scale / denominator;
    let numerator = value.numer().clone();
    numerator * multiplier
}

/// Model the expression used by the private production `integer_at_scale`.
fn scaled_borrowed(value: &BigRational, scale: &BigInt) -> BigInt {
    value.numer() * (scale / value.denom())
}

/// Compute LCM via rational reduction, independently of production's GCD loop.
fn row_scale<'a>(values: impl Iterator<Item = &'a BigRational>) -> BigInt {
    values.fold(BigInt::from(1), |scale, value| {
        let reduced = BigRational::new(scale, value.denom().clone());
        reduced.numer() * value.denom()
    })
}

fn record<const D: usize>(kind: RationalInputKind, operation: &str, run: impl Fn()) {
    let counts = measure(&run);
    // Counts must be repeatable, and each complete call drops its outputs.
    assert_eq!(counts.count_current, 0);
    assert_eq!(counts.bytes_current, 0);
    for _ in 0..2 {
        assert_eq!(measure(&run), counts);
    }
    let AllocationInfo {
        count_total,
        bytes_total,
        ..
    } = counts;
    println!(
        "allocation,{},{D},{operation},{count_total},{bytes_total}",
        kind.name()
    );
}

fn rational_allocation_evidence<const D: usize>() {
    for kind in RationalInputKind::ALL {
        let input = rational_input::<D>(kind);
        let matrix = input.matrix();
        let rhs = input.rhs();
        record::<D>(kind, "det", || {
            black_box(black_box(matrix).det());
        });
        record::<D>(kind, "det_sign", || {
            let _ = black_box(black_box(matrix).det_sign());
        });
        record::<D>(kind, "solve", || {
            let _ = black_box(black_box(matrix).solve(black_box(rhs)).unwrap());
        });

        for augmented in [false, true] {
            let values_in_row = |row: usize| {
                matrix.as_rows()[row]
                    .iter()
                    .chain(augmented.then_some(&rhs.as_array()[row]))
            };
            let scales: [BigInt; D] = from_fn(|row| row_scale(values_in_row(row)));
            let entries: Vec<_> = scales
                .iter()
                .enumerate()
                .flat_map(|(row, scale)| values_in_row(row).map(move |value| (value, scale)))
                .collect();
            for &(value, scale) in &entries {
                let expected = value * BigRational::from_integer(scale.clone());
                assert_eq!(
                    BigRational::from_integer(scaled_cloned(value, scale)),
                    expected
                );
                assert_eq!(
                    BigRational::from_integer(scaled_borrowed(value, scale)),
                    expected
                );
            }

            let prefix = if augmented { "augmented" } else { "matrix" };
            record::<D>(kind, &format!("{prefix}_component_clones"), || {
                for &(value, _) in black_box(&entries) {
                    black_box((value.numer().clone(), value.denom().clone()));
                }
            });
            record::<D>(kind, &format!("{prefix}_scaled_cloned"), || {
                for &(value, scale) in black_box(&entries) {
                    black_box(scaled_cloned(value, scale));
                }
            });
            record::<D>(kind, &format!("{prefix}_scaled_borrowed"), || {
                for &(value, scale) in black_box(&entries) {
                    black_box(scaled_borrowed(value, scale));
                }
            });
        }
    }
}

macro_rules! gen_rational_allocation_tests {
    ($d:literal) => {
        paste! {
            #[test]
            fn [<rational_allocation_evidence_ $d d>]() {
                rational_allocation_evidence::<$d>();
            }
        }
    };
}

gen_rational_allocation_tests!(2);
gen_rational_allocation_tests!(3);
gen_rational_allocation_tests!(4);
gen_rational_allocation_tests!(5);
gen_rational_allocation_tests!(6);
gen_rational_allocation_tests!(7);
gen_rational_allocation_tests!(8);
