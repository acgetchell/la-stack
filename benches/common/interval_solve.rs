#![forbid(unsafe_code)]

//! Deterministic systems with manufactured exact solutions, shared by the
//! timing harness and independent rational fixture tests. D is in 1..=6.

use std::array::from_fn;

use super::bench_utils::OrAbort;
use la_stack::{ArithmeticOperation, Interval, IntervalMatrix, LaError};

pub(crate) const CASES: usize = 32;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Outcome {
    Enclosed,
    Inconclusive,
    RangeExhausted,
}

pub(crate) struct Input<const D: usize> {
    pub(crate) matrix: IntervalMatrix<D>,
    pub(crate) rhs: [Interval; D],
    pub(crate) known_solution: [f64; D],
    pub(crate) expected: Outcome,
}

pub(crate) const FAMILIES: [&str; 5] = [
    "dense",
    "uncertain_dense",
    "near_singular",
    "singular",
    "range_exhausted",
];

fn point(value: f64) -> Interval {
    Interval::point(value).or_abort("finite benchmark coefficient")
}

pub(crate) fn input<const D: usize>(family: &str, index: usize) -> Input<D> {
    assert!((1..=6).contains(&D) && index < CASES);
    let known_solution = from_fn(|i| if i.is_multiple_of(2) { 1.0 } else { -1.0 });
    let mut centers = [[0.0; D]; D];
    let mut radii = [[0.0; D]; D];
    let expected = match family {
        "dense" | "uncertain_dense" => {
            for (i, row) in centers.iter_mut().enumerate() {
                for (j, value) in row.iter_mut().enumerate() {
                    *value = if (i + 1) % D == j {
                        f64::from(u32::try_from(8 * D).or_abort("small dimension"))
                            + f64::from(u32::try_from(index).or_abort("corpus index")) / 32.0
                    } else {
                        f64::from(
                            i32::try_from((i * D + j + index) % 9).or_abort("bounded coefficient")
                                - 4,
                        ) / 8.0
                    };
                    if family == "uncertain_dense" {
                        radii[i][j] = 2.0_f64.powi(-30);
                    }
                }
            }
            Outcome::Enclosed
        }
        "near_singular" => {
            for (i, row) in centers.iter_mut().enumerate() {
                row[i] = 1.0;
            }
            let delta = 2.0_f64.powi(-20 - i32::try_from(index).or_abort("corpus exponent"));
            let radius = 2.0_f64.powi(-48);
            if D == 1 {
                centers[0][0] = delta;
            } else {
                centers[D - 2][D - 1] = 1.0;
                centers[D - 1][D - 2] = 1.0;
                centers[D - 1][D - 1] = 1.0 + delta;
            }
            radii[D - 1][D - 1] = radius;
            if delta > radius {
                Outcome::Enclosed
            } else {
                Outcome::Inconclusive
            }
        }
        "singular" | "range_exhausted" => {
            for (i, row) in centers.iter_mut().enumerate().take(D - 1) {
                row[i] = 1.0;
            }
            if family == "singular" {
                Outcome::Inconclusive
            } else {
                centers[D - 1][D - 1] = f64::from_bits(1);
                Outcome::RangeExhausted
            }
        }
        _ => panic!("unknown interval solve family"),
    };
    // These small dyadic sums are exact in binary64; the range fixture replaces
    // its last RHS and deliberately has no finite enclosure.
    let mut rhs =
        centers.map(|row| point(row.iter().zip(known_solution).map(|(a, x)| a * x).sum()));
    if expected == Outcome::RangeExhausted {
        rhs[D - 1] = point(f64::MAX);
    }
    let rows = from_fn(|i| {
        from_fn(|j| {
            Interval::try_new(centers[i][j] - radii[i][j], centers[i][j] + radii[i][j])
                .or_abort("finite coefficient bounds")
        })
    });
    Input {
        matrix: IntervalMatrix::from_rows(rows),
        rhs,
        known_solution,
        expected,
    }
}

/// Require the intended outcome before registering a timing row; no oracle work
/// or assertions run in the timed closure.
pub(crate) fn validate<const D: usize>(input: &Input<D>) {
    match (input.expected, input.matrix.solve(&input.rhs)) {
        (Outcome::Enclosed, Ok(bounds)) => {
            for (bound, value) in bounds.into_iter().zip(input.known_solution) {
                assert!(bound.contains(value));
            }
        }
        (Outcome::Inconclusive, Err(LaError::IntervalSolveInconclusive { pivot_col, .. })) => {
            assert_eq!(pivot_col, D - 1);
        }
        (
            Outcome::RangeExhausted,
            Err(LaError::IntervalRangeExhausted {
                operation: ArithmeticOperation::IntervalSolve,
                ..
            }),
        ) => {}
        (expected, actual) => panic!("expected {expected:?}, got {actual:?}"),
    }
}
