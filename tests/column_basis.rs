#![forbid(unsafe_code)]

//! Analytical rectangular ranks, numerical boundaries, and exact certification.

use core::assert_matches;
use std::array::from_fn;

use la_stack::prelude::*;
use pastey::paste;

fn tolerance(value: f64) -> Tolerance {
    Tolerance::try_new(value).unwrap()
}

#[cfg(feature = "exact")]
fn certify_full_row_rank<const R: usize, const C: usize>(
    rows: &[[f64; C]; R],
    basis: &ColumnBasis<R, C>,
) {
    assert_eq!(basis.rank(), R);
    let selected =
        Matrix::<R>::try_from_rows(from_fn(|r| from_fn(|c| rows[r][basis.columns()[c]]))).unwrap();
    // Independent fraction-free integer determinant of the represented inputs.
    assert_ne!(
        selected.det_exact().unwrap(),
        BigRational::from_integer(0.into())
    );
}

#[test]
fn empty_and_zero_shapes() {
    assert_eq!(
        column_basis::<0, 5>(&[], tolerance(0.0)).unwrap().columns(),
        &[]
    );
    assert_eq!(
        column_basis::<5, 0>(&[[]; 5], tolerance(0.0))
            .unwrap()
            .rank(),
        0
    );
    assert_eq!(column_basis::<0, 0>(&[], tolerance(0.0)).unwrap().rank(), 0);
    assert_eq!(
        column_basis(&[[-0.0; 4]; 3], tolerance(0.0))
            .unwrap()
            .columns(),
        &[]
    );
}

#[test]
fn tall_wide_and_deficient_shapes() {
    let tall = [[1.0, 2.0], [2.0, 4.0], [0.0, 1.0], [3.0, 7.0]];
    assert_eq!(
        column_basis(&tall, tolerance(0.0)).unwrap().columns(),
        &[0, 1]
    );
    let wide = [[0.0, 1.0, 2.0, 0.0], [0.0, 2.0, 4.0, 1.0]];
    let basis = column_basis(&wide, tolerance(0.0)).unwrap();
    assert_eq!(basis.columns(), &[1, 3]);
    #[cfg(feature = "exact")]
    certify_full_row_rank(&wide, &basis);
    let deficient = [[1.0, 2.0, -3.0], [2.0, 4.0, -6.0], [0.0, 0.0, 0.0]];
    assert_eq!(
        column_basis(&deficient, tolerance(0.0)).unwrap().columns(),
        &[0]
    );
}

#[test]
fn strict_threshold_and_near_dependence_are_provisional() {
    let delta = 2.0_f64.powi(-40);
    let rows = [[1.0, 1.0], [1.0, 1.0 + delta]];
    assert_eq!(
        column_basis(&rows, tolerance(delta)).unwrap().columns(),
        &[0]
    );
    let basis = column_basis(&rows, tolerance(delta.next_down())).unwrap();
    assert_eq!(basis.columns(), &[0, 1]);
    #[cfg(feature = "exact")]
    certify_full_row_rank(&rows, &basis);
    assert_eq!(column_basis(&rows, tolerance(f64::MAX)).unwrap().rank(), 0);
}

#[test]
fn severe_uniform_and_nonuniform_rescaling_obeys_absolute_policy() {
    for scale in [
        f64::from_bits(1),
        2.0_f64.powi(-900),
        1.0,
        2.0_f64.powi(900),
        f64::MAX,
    ] {
        let rows = [[scale, 0.0, scale], [0.0, -scale, -scale]];
        assert_eq!(
            column_basis(&rows, tolerance(0.0)).unwrap().columns(),
            &[0, 1]
        );
        assert_eq!(column_basis(&rows, tolerance(scale)).unwrap().rank(), 0);
    }
    let rows = [[2.0_f64.powi(-900), 0.0], [0.0, 2.0_f64.powi(900)]];
    assert_eq!(
        column_basis(&rows, tolerance(1e-12)).unwrap().columns(),
        &[1]
    );
    assert_eq!(
        column_basis(&rows, tolerance(0.0)).unwrap().columns(),
        &[0, 1]
    );
}

#[test]
fn permutations_and_ties_have_deterministic_order() {
    let rows = [[0.0, 1.0, 2.0, 1.0], [0.0, -1.0, -2.0, 1.0]];
    for row_order in [[0, 1], [1, 0]] {
        for columns in [[0, 1, 2, 3], [3, 2, 0, 1], [2, 0, 1, 3]] {
            let permuted = row_order.map(|r| columns.map(|c| rows[r][c]));
            let expected = if columns[0] == 0 {
                [1, 3]
            } else if columns[0] == 3 {
                [0, 1]
            } else {
                [0, 3]
            };
            let basis = column_basis(&permuted, tolerance(0.0)).unwrap();
            assert_eq!(basis.columns(), &expected);
            assert_eq!(column_basis(&permuted, tolerance(0.0)).unwrap(), basis);
            #[cfg(feature = "exact")]
            certify_full_row_rank(&permuted, &basis);
        }
    }
    // Equal magnitudes pick the first current row: its second coordinate is
    // the pivot-row term in FMA. The residual is exactly 2^-52, not its negative.
    let rows = [[1.0, 1.0], [-1.0, -1.0 + f64::EPSILON]];
    assert_eq!(column_basis(&rows, tolerance(0.0)).unwrap().rank(), 2);
}

#[test]
fn all_input_cells_are_validated_before_early_rank_completion() {
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_matches!(
            column_basis(&[[1.0, bad]], tolerance(0.0)),
            Err(LaError::NonFinite {
                origin: NonFiniteOrigin::Input,
                location: NonFiniteLocation::MatrixCell { row: 0, col: 1, .. },
                ..
            })
        );
        assert_matches!(
            column_basis(&[[0.0, bad], [bad, 1.0]], tolerance(0.0)),
            Err(LaError::NonFinite {
                origin: NonFiniteOrigin::Input,
                location: NonFiniteLocation::MatrixCell { row: 0, col: 1, .. },
                ..
            })
        );
    }
}

#[test]
fn nonfinite_updates_report_work_cell_and_preserve_input() {
    let rows = [[1.0, f64::MAX], [-1.0, f64::MAX]];
    assert_matches!(
        column_basis(&rows, tolerance(0.0)),
        Err(LaError::NonFinite {
            origin: NonFiniteOrigin::Computation {
                operation: ArithmeticOperation::ColumnBasis,
                ..
            },
            location: NonFiniteLocation::MatrixCell { row: 1, col: 1, .. },
            ..
        })
    );
    assert_eq!(rows, [[1.0, f64::MAX], [-1.0, f64::MAX]]);
    // All first-column magnitudes tie. Choosing the first current row must
    // encounter this overflow; choosing the last row would avoid it. Reordering
    // the zero second-coordinate row first therefore gives a finite selection.
    let tied = [[1.0, f64::MAX], [1.0, -f64::MAX], [1.0, 0.0]];
    assert_matches!(
        column_basis(&tied, tolerance(0.0)),
        Err(LaError::NonFinite {
            origin: NonFiniteOrigin::Computation {
                operation: ArithmeticOperation::ColumnBasis,
                ..
            },
            location: NonFiniteLocation::MatrixCell { row: 1, col: 1, .. },
            ..
        })
    );
    assert_eq!(
        column_basis(&[tied[2], tied[1], tied[0]], tolerance(0.0))
            .unwrap()
            .columns(),
        &[0, 1]
    );
    assert_eq!(
        ArithmeticOperation::ColumnBasis.to_string(),
        "column basis selection"
    );
}

fn affine_shared_face<const R: usize, const C: usize>() {
    // Rows are [point, 1] for {0, e_0, ..., e_(R-2)} in ambient C-1.
    // Their affine independence is analytical; the selected minor is ±1.
    let rows: [[f64; C]; R] =
        from_fn(|r| from_fn(|c| f64::from(c == C - 1 || (r > 0 && c == r - 1))));
    let expected: Vec<_> = (0..R - 1).chain([C - 1]).collect();
    let basis = column_basis(&rows, tolerance(1e-12)).unwrap();
    assert_eq!(basis.rank(), R);
    assert_eq!(basis.columns(), expected);
    #[cfg(feature = "exact")]
    certify_full_row_rank(&rows, &basis);
}

fn simplex_constraints<const R: usize, const C: usize>() {
    // R=D+2, C=2(D+1): two barycentric sum equations followed by
    // coordinate agreement for standard simplices on opposite sides of a facet.
    let vertices = C / 2;
    let rows: [[f64; C]; R] = from_fn(|r| {
        from_fn(|c| {
            if r < 2 {
                return f64::from(c / vertices == r);
            }
            let vertex = c % vertices;
            if vertex != r - 1 {
                return 0.0;
            }
            if c < vertices || r == R - 1 {
                1.0
            } else {
                -1.0
            }
        })
    });
    let basis = column_basis(&rows, tolerance(1e-12)).unwrap();
    assert_eq!(basis.columns(), (0..R).collect::<Vec<_>>());
    #[cfg(feature = "exact")]
    certify_full_row_rank(&rows, &basis);
}

macro_rules! geometry_dimensions {
    ($d:literal, $rows:literal, $cols:literal, $affine_cols:literal; $($face:literal),+) => { paste! {
        #[test]
        fn [<shared_faces_ $d d>]() { $(affine_shared_face::<$face, $affine_cols>();)+ }
        #[test]
        fn [<intersection_constraints_ $d d>]() { simplex_constraints::<$rows, $cols>(); }
    }};
}

geometry_dimensions!(2, 4, 6, 3; 1, 2);
geometry_dimensions!(3, 5, 8, 4; 1, 2, 3);
geometry_dimensions!(4, 6, 10, 5; 1, 2, 3, 4);
geometry_dimensions!(5, 7, 12, 6; 1, 2, 3, 4, 5);
geometry_dimensions!(6, 8, 14, 7; 1, 2, 3, 4, 5, 6);

#[test]
fn allocation_free_and_const_evaluable() {
    const BASIS: Result<ColumnBasis<2, 3>, LaError> = match Tolerance::try_new(0.0) {
        Ok(tol) => column_basis(&[[0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], tol),
        Err(error) => Err(error),
    };
    assert_eq!(BASIS.unwrap().columns(), &[1, 2]);
    let allocations = allocation_counter::measure(|| {
        let rows = std::hint::black_box([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        assert_eq!(
            std::hint::black_box(column_basis(&rows, tolerance(0.0)).unwrap()).rank(),
            2
        );
    });
    assert_eq!(allocations.count_total, 0);
}
