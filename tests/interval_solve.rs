#![forbid(unsafe_code)]

//! Solution enclosures, fallback outcomes, and independent rational oracles.

use core::assert_matches;
use std::array::from_fn;

use la_stack::prelude::*;
use pastey::paste;

fn point(value: f64) -> Interval {
    Interval::point(value).expect("finite fixture")
}

#[test]
fn division_domain_and_range_failures_are_distinct() {
    for divisor in [(-1.0, 1.0), (0.0, 0.0), (-0.0, 2.0), (-2.0, -0.0)] {
        let divisor = Interval::try_new(divisor.0, divisor.1).unwrap();
        for result in [
            divisor.try_recip(),
            Interval::ZERO.try_div(&divisor),
            Interval::ONE.try_div(&divisor),
        ] {
            match result {
                Err(LaError::IntervalZeroDivisor { lower, upper, .. }) => {
                    assert_eq!(lower, divisor.lower());
                    assert_eq!(upper, divisor.upper());
                }
                other => panic!("expected zero-divisor failure: {other:?}"),
            }
        }
    }
    let tiny = point(f64::from_bits(1));
    assert_eq!(tiny.try_div(&tiny), Ok(Interval::ONE));
    assert_eq!(
        tiny.try_recip(),
        Err(LaError::interval_range_exhausted(
            ArithmeticOperation::IntervalReciprocal
        ))
    );
    assert_eq!(
        point(f64::MAX).try_div(&point(0.5)),
        Err(LaError::interval_range_exhausted(
            ArithmeticOperation::IntervalDivision
        ))
    );
    assert_eq!(
        tiny.try_div(&point(2.0)),
        Interval::try_new(0.0, f64::from_bits(1))
    );
    assert_eq!(
        tiny.negate().try_div(&point(2.0)),
        Interval::try_new(-f64::from_bits(1), 0.0)
    );
    assert_eq!(point(-0.0).try_div(&point(-3.0)), Ok(Interval::ZERO));
}

#[test]
fn constant_evaluation_and_empty_solution() {
    const RECIPROCAL: Result<Interval, LaError> = Interval::ONE.try_recip();
    const SOLUTION: Result<[Interval; 2], LaError> =
        IntervalMatrix::<2>::identity().solve(&[Interval::ONE; 2]);
    assert_eq!(RECIPROCAL, Ok(Interval::ONE));
    assert_eq!(SOLUTION, Ok([Interval::ONE; 2]));
    assert_eq!(IntervalMatrix::<0>::zero().solve(&[]), Ok([]));
    // The determinant's subset-workspace limit does not apply to solves.
    assert_eq!(
        IntervalMatrix::<8>::identity().solve(&[Interval::ONE; 8]),
        Ok([Interval::ONE; 8])
    );
}

fn known_solution<const D: usize>() {
    let rows = from_fn(|row| {
        from_fn(|col| {
            if col == (row + 1) % D {
                point(-2.0)
            } else {
                Interval::ZERO
            }
        })
    });
    let matrix = IntervalMatrix::<D>::from_rows(rows);
    let rhs = [point(-6.0); D];
    assert_eq!(matrix.solve(&rhs), Ok([point(3.0); D]));
    assert_eq!(matrix.as_rows(), &rows);
    assert_eq!(rhs, [point(-6.0); D]);
    let uncertain_rhs = [Interval::try_new(-6.0, -4.0).unwrap(); D];
    assert_eq!(
        matrix.solve(&uncertain_rhs),
        Ok([Interval::try_new(2.0, 3.0).unwrap(); D])
    );
}

macro_rules! dimension_tests {
    ($d:literal) => {
        paste! {
            #[test]
            fn [<known_solution_d $d>]() { known_solution::<$d>(); }
        }
    };
}
dimension_tests!(1);
dimension_tests!(2);
dimension_tests!(3);
dimension_tests!(4);
dimension_tests!(5);
dimension_tests!(6);

#[test]
fn pivots_use_separation_and_fail_inconclusively_including_last_column() {
    let wide = Interval::try_new(-100.0, 100.0).unwrap();
    let matrix =
        IntervalMatrix::from_rows([[wide, Interval::ONE], [Interval::ONE, Interval::ZERO]]);
    let solution = matrix.solve(&[Interval::ONE, Interval::ZERO]).unwrap();
    assert_eq!(solution, [Interval::ZERO, Interval::ONE]);
    for rhs in [[Interval::ONE; 2], [Interval::ONE, Interval::ZERO]] {
        let singular = IntervalMatrix::try_from_point_rows([[1.0, 1.0], [1.0, 1.0]]).unwrap();
        assert_matches!(
            singular.solve(&rhs),
            Err(LaError::IntervalSolveInconclusive { pivot_col: 1, .. })
        );
    }
    let uncertain = IntervalMatrix::from_rows([[Interval::try_new(-1.0, 1.0).unwrap()]]);
    assert_matches!(
        uncertain.solve(&[Interval::ZERO]),
        Err(LaError::IntervalSolveInconclusive { pivot_col: 0, .. })
    );
    // Exact determinant = 3 * 2^-52 > 0, but the rounded elimination encloses
    // zero in its final pivot. Inconclusive must not be confused with singular.
    let regular =
        IntervalMatrix::try_from_point_rows([[3.0, 3.0], [1.0, 1.0_f64.next_up()]]).unwrap();
    assert_matches!(
        regular.solve(&[Interval::ONE; 2]),
        Err(LaError::IntervalSolveInconclusive { pivot_col: 1, .. })
    );
    // A point matrix with a tiny nonzero determinant must still succeed.
    let almost_singular =
        IntervalMatrix::try_from_point_rows([[1.0, 1.0], [1.0, 1.0_f64.next_up()]]).unwrap();
    let x = almost_singular.solve(&[point(2.0), point(2.0)]).unwrap();
    assert!(x[0].contains(2.0) && x[1].contains(0.0));
}

#[test]
fn solve_extreme_scale_and_typed_range_failures() {
    for scale in [f64::from_bits(1), f64::MIN_POSITIVE, f64::MAX] {
        let matrix = IntervalMatrix::try_from_point_rows([[scale]]).unwrap();
        assert_eq!(matrix.solve(&[point(scale)]), Ok([Interval::ONE]));
    }
    let matrix = IntervalMatrix::try_from_point_rows([[2.0]]).unwrap();
    assert_eq!(
        matrix.solve(&[point(f64::from_bits(1))]),
        Ok([Interval::try_new(0.0, f64::from_bits(1)).unwrap()])
    );
    for (matrix, rhs) in [
        ([[0.5, 0.0], [0.0, 1.0]], [f64::MAX, 0.0]),
        ([[1.0, f64::MAX], [1.0, -f64::MAX]], [0.0, 0.0]),
        ([[1.0, f64::MAX], [0.0, 1.0]], [0.0, 2.0]),
    ] {
        let matrix = IntervalMatrix::try_from_point_rows(matrix).unwrap();
        let before = matrix;
        assert_eq!(
            matrix.solve(&rhs.map(point)),
            Err(LaError::interval_range_exhausted(
                ArithmeticOperation::IntervalSolve
            ))
        );
        assert_eq!(matrix, before);
    }
}

#[test]
fn solve_is_allocation_free_and_errors_remain_descriptive() {
    let matrix = IntervalMatrix::<6>::identity();
    let info = allocation_counter::measure(|| {
        assert_eq!(matrix.solve(&[Interval::ONE; 6]), Ok([Interval::ONE; 6]));
    });
    assert_eq!(info.count_total, 0);
    assert_eq!(
        ArithmeticOperation::IntervalDivision.to_string(),
        "interval division"
    );
    assert_eq!(
        ArithmeticOperation::IntervalReciprocal.to_string(),
        "interval reciprocal"
    );
    assert_eq!(
        ArithmeticOperation::IntervalSolve.to_string(),
        "interval solve"
    );
    assert_eq!(
        Interval::try_new(-1.0, 2.0)
            .unwrap()
            .try_recip()
            .unwrap_err()
            .to_string(),
        "interval divisor [-1, 2] contains zero"
    );
    assert_eq!(
        IntervalMatrix::<1>::zero()
            .solve(&[Interval::ONE])
            .unwrap_err()
            .to_string(),
        "interval solve is inconclusive at pivot column 0; no pivot excludes zero"
    );
    assert_matches!(
        Interval::point(f64::NAN),
        Err(LaError::NonFinite {
            origin: NonFiniteOrigin::Input,
            location: NonFiniteLocation::Scalar,
            ..
        })
    );
}

#[cfg(feature = "exact")]
#[path = "common/interval_solve_oracle.rs"]
mod oracle;

#[cfg(feature = "exact")]
#[path = "common/proptest_config.rs"]
mod proptest_config;

#[path = "../benches/common/interval_solve.rs"]
mod solve_fixtures;

#[path = "../benches/common/bench_utils.rs"]
mod bench_utils;

fn benchmark_fixtures<const D: usize>() {
    for family in solve_fixtures::FAMILIES {
        for index in 0..solve_fixtures::CASES {
            let input = solve_fixtures::input::<D>(family, index);
            solve_fixtures::validate(&input);
            #[cfg(feature = "exact")]
            {
                use oracle::{assert_encloses, exact, exact_solve};
                // The midpoint system recovers the exact dyadic construction;
                // solve it independently before accepting benchmark evidence.
                let midpoint = |x: Interval| (exact(x.lower()) + exact(x.upper())) / exact(2.0);
                let rows = input.matrix.into_rows().map(|r| r.map(midpoint));
                let rhs = input.rhs.map(midpoint);
                let exact_solution = exact_solve(rows, rhs);
                if family == "singular" {
                    assert!(exact_solution.is_none());
                } else {
                    let exact_solution = exact_solution.unwrap();
                    if family == "range_exhausted" {
                        assert!(exact_solution[D - 1].abs() > exact(f64::MAX));
                    } else {
                        assert_eq!(exact_solution, input.known_solution.map(exact));
                        if input.expected == solve_fixtures::Outcome::Enclosed {
                            assert_encloses(
                                &input.matrix.solve(&input.rhs).unwrap(),
                                &exact_solution,
                            );
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn benchmark_corpus_has_independently_validated_outcomes() {
    benchmark_fixtures::<1>();
    benchmark_fixtures::<2>();
    benchmark_fixtures::<3>();
    benchmark_fixtures::<4>();
    benchmark_fixtures::<5>();
    benchmark_fixtures::<6>();
}

#[cfg(feature = "exact")]
mod exact_tests {
    use super::oracle::{assert_encloses, exact, exact_solve};
    use super::proptest_config::with_default_cases;
    use super::*;
    use proptest::prelude::*;

    fn check_point_system<const D: usize>(a: [[i8; D]; D], b: [i8; D]) {
        let a = a.map(|row| row.map(|x| f64::from(x) / 10.0));
        let b = b.map(|x| f64::from(x) / 10.0);
        let exact_solution = exact_solve(a.map(|row| row.map(exact)), b.map(exact));
        match IntervalMatrix::try_from_point_rows(a)
            .unwrap()
            .solve(&b.map(point))
        {
            Ok(solution) => assert_encloses(
                &solution,
                &exact_solution.expect("a successful enclosure proves regularity"),
            ),
            Err(LaError::IntervalSolveInconclusive { pivot_col, .. }) => assert!(pivot_col < D),
            other => panic!("unexpected result for bounded random system: {other:?}"),
        }
    }

    fn check_uncertain_system<const D: usize>(a: [[i8; D]; D], b: [i8; D]) {
        let radius = 1.0 / 1024.0;
        let rows = from_fn(|i| {
            from_fn(|j| {
                let center = if i == j {
                    f64::from(u32::try_from(4 * D).unwrap())
                } else {
                    f64::from(a[i][j]) / 128.0
                };
                Interval::try_new(center - radius, center + radius).unwrap()
            })
        });
        let rhs =
            b.map(|x| Interval::try_new(f64::from(x) - radius, f64::from(x) + radius).unwrap());
        let solution = IntervalMatrix::from_rows(rows)
            .solve(&rhs)
            .expect("small diagonally dominant intervals must certify");
        // Sample lower, upper, checkerboard, and strictly interior rational
        // coefficients. Each sample solves a different admitted exact system.
        for sample in 0..4 {
            let pick = |x: Interval, index: usize| match sample {
                0 => exact(x.lower()),
                2 if index.is_multiple_of(2) => exact(x.lower()),
                1 | 2 => exact(x.upper()),
                _ => (exact(x.lower()) + exact(x.upper()) * exact(2.0)) / exact(3.0),
            };
            assert_encloses(
                &solution,
                &exact_solve(
                    from_fn(|i| from_fn(|j| pick(rows[i][j], i * D + j))),
                    from_fn(|i| pick(rhs[i], D * D + i)),
                )
                .unwrap(),
            );
        }
    }

    fn circumsphere_system<const D: usize>() {
        // Consumer-style equation assembly stays in tests, outside the API.
        // Translation by non-dyadic stored coordinates makes subtraction and
        // squared-norm uncertainty part of the system being certified.
        for scale in [1.0, 2.0_f64.powi(-500), 2.0_f64.powi(500)] {
            let origin: [f64; D] =
                from_fn(|i| f64::from(u32::try_from(i + 1).unwrap()) * 0.1 * scale);
            let vertices: [[f64; D]; D] = from_fn(|i| {
                from_fn(|j| {
                    origin[j]
                        + if (i + 1) % D == j {
                            2.0 * scale
                        } else {
                            0.25 * scale
                        }
                })
            });
            let rows = from_fn(|i| {
                from_fn(|j| Interval::try_from_subtraction(vertices[i][j], origin[j]).unwrap())
            });
            let rhs = rows.map(|row| {
                row.into_iter()
                    .fold(Interval::ZERO, |sum, x| {
                        sum.try_add(&x.try_square().unwrap()).unwrap()
                    })
                    .try_mul(&point(0.5))
                    .unwrap()
            });
            let exact_rows: [[BigRational; D]; D] =
                from_fn(|i| from_fn(|j| exact(vertices[i][j]) - exact(origin[j])));
            let exact_rhs =
                from_fn(|i| exact_rows[i].iter().map(|x| x * x).sum::<BigRational>() / exact(2.0));
            let solution = IntervalMatrix::from_rows(rows).solve(&rhs).unwrap();
            assert_encloses(&solution, &exact_solve(exact_rows, exact_rhs).unwrap());
        }
    }

    macro_rules! geometry_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<circumsphere_equations_d $d>]() { circumsphere_system::<$d>(); }
                proptest! {
                    #![proptest_config(with_default_cases(24))]
                    #[test]
                    fn [<point_systems_match_rational_oracle_d $d>](a in any::<[[i8; $d]; $d]>(), b in any::<[i8; $d]>()) {
                        check_point_system(a, b);
                    }
                    #[test]
                    fn [<uncertain_systems_enclose_rational_samples_d $d>](a in any::<[[i8; $d]; $d]>(), b in any::<[i8; $d]>()) {
                        check_uncertain_system(a, b);
                    }
                }
            }
        };
    }
    geometry_tests!(1);
    geometry_tests!(2);
    geometry_tests!(3);
    geometry_tests!(4);
    geometry_tests!(5);
    geometry_tests!(6);

    #[test]
    fn uncertain_two_by_two_encloses_every_endpoint_system_and_interior() {
        let rows = [
            [
                Interval::try_new(3.0, 4.0).unwrap(),
                Interval::try_new(-0.25, 0.25).unwrap(),
            ],
            [
                Interval::try_new(-0.5, 0.5).unwrap(),
                Interval::try_new(-4.0, -3.0).unwrap(),
            ],
        ];
        let rhs = [
            Interval::try_new(-2.0, 1.0).unwrap(),
            Interval::try_new(1.0, 3.0).unwrap(),
        ];
        let solution = IntervalMatrix::from_rows(rows).solve(&rhs).unwrap();
        for mask in 0_u32..64 {
            let pick = |value: Interval, bit: usize| {
                if mask & (1 << bit) == 0 {
                    exact(value.lower())
                } else {
                    exact(value.upper())
                }
            };
            let a = from_fn(|i| from_fn(|j| pick(rows[i][j], i * 2 + j)));
            let b = from_fn(|i| pick(rhs[i], 4 + i));
            assert_encloses(&solution, &exact_solve(a, b).unwrap());
        }
        let interior =
            |x: Interval| (exact(x.lower()) + exact(x.upper()) * exact(2.0)) / exact(3.0);
        assert_encloses(
            &solution,
            &exact_solve(rows.map(|r| r.map(interior)), rhs.map(interior)).unwrap(),
        );
    }
}
