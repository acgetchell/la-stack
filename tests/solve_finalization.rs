#![forbid(unsafe_code)]

//! Regression coverage for finite solve results and substitution error order.

use la_stack::{ArithmeticOperation, DEFAULT_SINGULAR_TOL, LaError, Matrix, Vector};
use pastey::paste;

fn assert_diagonal_overflow<const D: usize>() {
    // Exercise every output coordinate, including the last one finalized by LU.
    for index in 0..D {
        let mut rows = Matrix::<D>::identity().into_rows();
        rows[index][index] = 1.0e-11;
        let a = Matrix::try_from_rows(rows).unwrap();
        let mut rhs = [0.0; D];
        rhs[index] = 1.0e300;
        let b = Vector::try_new(rhs).unwrap();
        assert_eq!(
            a.lu(DEFAULT_SINGULAR_TOL).unwrap().solve(b),
            Err(LaError::non_finite_computation_step(
                ArithmeticOperation::LuSolve,
                index,
            )),
        );
        assert_eq!(
            a.ldlt(DEFAULT_SINGULAR_TOL).unwrap().solve(b),
            Err(LaError::non_finite_computation_step(
                ArithmeticOperation::LdltSolve,
                index,
            )),
        );
    }
}

fn assert_forward_overflow<const D: usize>() {
    let mut rows = Matrix::<D>::identity().into_rows();
    rows[D - 1][0] = -1.0;
    let mut rhs = [0.0; D];
    rhs[0] = f64::MAX;
    rhs[D - 1] = f64::MAX;
    let b = Vector::try_new(rhs).unwrap();
    let lu = Matrix::try_from_rows(rows)
        .unwrap()
        .lu(DEFAULT_SINGULAR_TOL)
        .unwrap();
    assert_eq!(
        lu.solve(b),
        Err(LaError::non_finite_computation_step(
            ArithmeticOperation::LuSolve,
            D - 1,
        )),
    );

    // A = L Lᵀ with L[D-1, 0] = -1 and unit diagonal, so the same
    // forward substitution overflows before either solve reaches finalization.
    rows[0][D - 1] = -1.0;
    rows[D - 1][D - 1] = 2.0;
    let ldlt = Matrix::try_from_rows(rows)
        .unwrap()
        .ldlt(DEFAULT_SINGULAR_TOL)
        .unwrap();
    assert_eq!(
        ldlt.solve(b),
        Err(LaError::non_finite_computation_step(
            ArithmeticOperation::LdltSolve,
            D - 1,
        )),
    );
}

fn assert_back_substitution_overflow<const D: usize>() {
    // Test each possible failing row. In the all-rows case, a final ascending
    // scan would report 0 instead of the required descending solve step D-2.
    for failing_row in 0..D - 1 {
        for all_rows in [false, true] {
            let mut upper = Matrix::<D>::identity().into_rows();
            let mut spd = upper;
            let mut last_diagonal = 1.0;
            for i in 0..D - 1 {
                if all_rows || i == failing_row {
                    upper[i][D - 1] = 2.0;
                    spd[i][D - 1] = 2.0;
                    spd[D - 1][i] = 2.0;
                    last_diagonal += 4.0;
                }
            }
            // A = L Lᵀ, where L's final row contains the selected multipliers.
            spd[D - 1][D - 1] = last_diagonal;
            let mut rhs = [0.0; D];
            rhs[D - 1] = 1.0e308;
            let b = Vector::try_new(rhs).unwrap();
            let expected_index = if all_rows { D - 2 } else { failing_row };
            let lu = Matrix::try_from_rows(upper)
                .unwrap()
                .lu(DEFAULT_SINGULAR_TOL)
                .unwrap();
            assert_eq!(
                lu.solve(b),
                Err(LaError::non_finite_computation_step(
                    ArithmeticOperation::LuSolve,
                    expected_index,
                )),
            );
            let ldlt = Matrix::try_from_rows(spd)
                .unwrap()
                .ldlt(DEFAULT_SINGULAR_TOL)
                .unwrap();
            assert_eq!(
                ldlt.solve(b),
                Err(LaError::non_finite_computation_step(
                    ArithmeticOperation::LdltSolve,
                    expected_index,
                )),
            );
        }
    }
}

macro_rules! gen_tests {
    ($d:literal) => {
        paste! {
            #[test]
            fn [<diagonal_overflow_ $d d>]() {
                assert_diagonal_overflow::<$d>();
            }

            #[test]
            fn [<forward_overflow_ $d d>]() {
                assert_forward_overflow::<$d>();
            }

            #[test]
            fn [<back_substitution_overflow_ $d d>]() {
                assert_back_substitution_overflow::<$d>();
            }
        }
    };
}

gen_tests!(2);
gen_tests!(3);
gen_tests!(4);
gen_tests!(5);
gen_tests!(8);
gen_tests!(16);
gen_tests!(32);
gen_tests!(64);
