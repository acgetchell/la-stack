//! Property-based tests for the exact-arithmetic APIs
//! (requires `exact` feature).
//!
//! Covers:
//! - `det_sign_exact` on diagonal and full small-integer matrices
//! - `solve_exact` round-trip with integer inputs (`A · x0` in f64 is
//!   exact for small integers, so `solve(A, A · x0) == x0`)
//! - `solve_exact` residual property (`A · solve(A, b) == b` in
//!   `BigRational` arithmetic) on random RHS vectors

#![cfg(feature = "exact")]

use std::array::from_fn;

use pastey::paste;
use proptest::prelude::*;

use la_stack::prelude::*;

fn small_nonzero_f64() -> impl Strategy<Value = f64> {
    prop_oneof![(-1000i16..=-1i16), (1i16..=1000i16)].prop_map(|x| f64::from(x) / 10.0)
}

/// Small non-zero integers in `[-5, 5] \ {0}`, returned as `f64`.
///
/// Used to populate the off-diagonal of diagonally-dominant integer
/// matrices in the `solve_exact` proptests.  Every value is an exact
/// `f64` integer, so `A · x` never loses precision for small `x`.
fn small_nonzero_int_f64() -> impl Strategy<Value = f64> {
    prop_oneof![(-5i32..=-1i32), (1i32..=5i32)].prop_map(f64::from)
}

/// Small signed integers in `[-10, 10]` as `f64` (may be zero).
fn small_int_f64() -> impl Strategy<Value = f64> {
    (-10i32..=10i32).prop_map(f64::from)
}

/// Multiply `A · x` entirely in `BigRational`, lifting each f64 matrix
/// entry via `BigRational::from_f64`.  Used by residual assertions.
///
/// All f64 inputs in the proptests are small exact integers, so
/// `from_f64` always succeeds with an exact rational reconstruction.
fn bigrational_matvec<const D: usize>(a: &[[f64; D]; D], x: &[BigRational; D]) -> [BigRational; D] {
    from_fn(|i| {
        let mut sum = BigRational::from_integer(BigInt::from(0));
        for (aij, xj) in a[i].iter().zip(x.iter()) {
            let entry = BigRational::from_f64(*aij).expect("small int fits in BigRational");
            sum += entry * xj;
        }
        sum
    })
}

/// Build a strictly diagonally-dominant f64 matrix from:
/// - an off-diagonal matrix of small integers (entries in `[-10, 10]`
///   per `small_int_f64`), and
/// - diagonal entries shifted by `D · 10 + 1` so every row satisfies
///   `|A[i][i]| > Σ_{j≠i} |A[i][j]|`, which guarantees invertibility
///   (Levy–Desplanques).
///
/// The shift must match the off-diagonal strategy's maximum magnitude
/// (`max_off_diag = 10`): with `D - 1` off-diagonals of magnitude ≤ 10
/// the row sum is at most `10 (D - 1) < 10 D + 1`, so the shifted
/// diagonal strictly dominates.  The shift keeps every entry a small
/// exact `f64` integer, so matrix × small-integer-vector products are
/// exact in f64.
fn make_diagonally_dominant<const D: usize>(
    offdiag: [[f64; D]; D],
    diag: [f64; D],
) -> [[f64; D]; D] {
    let mut rows = offdiag;
    // Must track `small_int_f64`'s `max_off_diag = 10`: `D · 10 + 1`
    // strictly dominates the worst-case row sum of `10 (D - 1)`.
    let shift = f64::from(u8::try_from(D).unwrap_or(u8::MAX)).mul_add(10.0, 1.0);
    for i in 0..D {
        rows[i][i] = if diag[i] >= 0.0 {
            diag[i] + shift
        } else {
            diag[i] - shift
        };
    }
    rows
}

macro_rules! gen_det_sign_exact_proptests {
    ($d:literal) => {
        paste! {
            proptest! {
                #![proptest_config(ProptestConfig::with_cases(64))]

                #[test]
                fn [<det_sign_exact_agrees_with_det_for_diagonal_ $d d>](
                    diag in proptest::array::[<uniform $d>](small_nonzero_f64()),
                ) {
                    // Diagonal matrix: determinant sign = product of diagonal signs.
                    let mut rows = [[0.0f64; $d]; $d];
                    for i in 0..$d {
                        rows[i][i] = diag[i];
                    }
                    let m = Matrix::<$d>::from_rows(rows);

                    let exact_sign = m.det_sign_exact().unwrap();

                    // Expected sign from the product of diagonal entries.
                    let neg_count = diag.iter().filter(|&&x| x < 0.0).count();
                    let expected_sign: i8 = if neg_count % 2 == 0 { 1 } else { -1 };

                    prop_assert_eq!(exact_sign, expected_sign);
                }

                #[test]
                fn [<det_sign_exact_agrees_with_det_signum_ $d d>](
                    diag in proptest::array::[<uniform $d>](small_nonzero_f64()),
                ) {
                    // For well-conditioned diagonal matrices, det().signum()
                    // should agree with det_sign_exact().
                    let mut rows = [[0.0f64; $d]; $d];
                    for i in 0..$d {
                        rows[i][i] = diag[i];
                    }
                    let m = Matrix::<$d>::from_rows(rows);

                    let exact_sign = m.det_sign_exact().unwrap();
                    let fp_det = m.det(DEFAULT_PIVOT_TOL).unwrap();
                    let fp_sign: i8 = if fp_det > 0.0 {
                        1
                    } else if fp_det < 0.0 {
                        -1
                    } else {
                        0
                    };

                    prop_assert_eq!(exact_sign, fp_sign);
                }
            }
        }
    };
}

gen_det_sign_exact_proptests!(2);
gen_det_sign_exact_proptests!(3);
gen_det_sign_exact_proptests!(4);
gen_det_sign_exact_proptests!(5);

/// Round-trip property: for random small-integer `x0` and a random
/// diagonally-dominant integer matrix `A`, `A · x0` is exactly
/// representable in `f64` (small integer products stay well under the
/// 53-bit mantissa) and `solve_exact(A, A · x0)` must return `x0`
/// exactly as `BigRational`.  This exercises the full Bareiss forward
/// elimination + rational back-substitution pipeline on a different
/// input for every case.
macro_rules! gen_solve_exact_roundtrip_proptests {
    ($d:literal) => {
        paste! {
            proptest! {
                #![proptest_config(ProptestConfig::with_cases(64))]

                #[test]
                fn [<solve_exact_integer_roundtrip_ $d d>](
                    offdiag in proptest::array::[<uniform $d>](
                        proptest::array::[<uniform $d>](small_int_f64()),
                    ),
                    diag in proptest::array::[<uniform $d>](small_nonzero_int_f64()),
                    x0 in proptest::array::[<uniform $d>](small_int_f64()),
                ) {
                    let rows = make_diagonally_dominant::<$d>(offdiag, diag);
                    let a = Matrix::<$d>::from_rows(rows);

                    // b = A · x0, computed in f64.  Small integers keep
                    // every partial sum exact.
                    let mut b_arr = [0.0f64; $d];
                    for i in 0..$d {
                        let mut sum = 0.0f64;
                        for j in 0..$d {
                            sum += rows[i][j] * x0[j];
                        }
                        b_arr[i] = sum;
                    }
                    let b = Vector::<$d>::new(b_arr);
                    let x = a.solve_exact(b).expect("diagonally-dominant A is non-singular");

                    let expected: [BigRational; $d] = from_fn(|i| {
                        BigRational::from_f64(x0[i]).expect("small int fits in BigRational")
                    });
                    for i in 0..$d {
                        prop_assert_eq!(&x[i], &expected[i]);
                    }
                }
            }
        }
    };
}

gen_solve_exact_roundtrip_proptests!(2);
gen_solve_exact_roundtrip_proptests!(3);
gen_solve_exact_roundtrip_proptests!(4);
gen_solve_exact_roundtrip_proptests!(5);

/// Residual property: for a random diagonally-dominant integer matrix
/// `A` and a random integer RHS `b`, `solve_exact` must return an `x`
/// such that `A · x` equals `b` *exactly* in `BigRational`
/// arithmetic.  Unlike the round-trip test above, the exact solution
/// is generally fractional — this catches back-substitution bugs that
/// preserve integer inputs but mishandle denominators.
macro_rules! gen_solve_exact_residual_proptests {
    ($d:literal) => {
        paste! {
            proptest! {
                #![proptest_config(ProptestConfig::with_cases(32))]

                #[test]
                fn [<solve_exact_residual_ $d d>](
                    offdiag in proptest::array::[<uniform $d>](
                        proptest::array::[<uniform $d>](small_int_f64()),
                    ),
                    diag in proptest::array::[<uniform $d>](small_nonzero_int_f64()),
                    b_arr in proptest::array::[<uniform $d>](small_int_f64()),
                ) {
                    let rows = make_diagonally_dominant::<$d>(offdiag, diag);
                    let a = Matrix::<$d>::from_rows(rows);
                    let b = Vector::<$d>::new(b_arr);
                    let x = a.solve_exact(b).expect("diagonally-dominant A is non-singular");

                    let ax = bigrational_matvec::<$d>(&rows, &x);
                    for i in 0..$d {
                        let b_rat = BigRational::from_f64(b_arr[i])
                            .expect("small int fits in BigRational");
                        prop_assert_eq!(&ax[i], &b_rat);
                    }
                }
            }
        }
    };
}

gen_solve_exact_residual_proptests!(2);
gen_solve_exact_residual_proptests!(3);
gen_solve_exact_residual_proptests!(4);
gen_solve_exact_residual_proptests!(5);

/// On full (non-diagonal) random small-integer matrices,
/// `det_sign_exact()` must agree with `det_exact().signum()`.  This
/// exercises the adaptive fast-filter / Bareiss-fallback boundary on
/// inputs the existing diagonal-only proptests don't touch (e.g.
/// matrices where the f64 det is near its error bound and the filter
/// must defer to Bareiss).
macro_rules! gen_det_sign_agrees_with_det_exact_proptests {
    ($d:literal) => {
        paste! {
            proptest! {
                #![proptest_config(ProptestConfig::with_cases(64))]

                #[test]
                fn [<det_sign_exact_agrees_with_det_exact_ $d d>](
                    entries in proptest::array::[<uniform $d>](
                        proptest::array::[<uniform $d>](small_int_f64()),
                    ),
                ) {
                    let m = Matrix::<$d>::from_rows(entries);
                    let sign = m.det_sign_exact().unwrap();
                    let det = m.det_exact().unwrap();
                    let expected: i8 = if det.is_positive() {
                        1
                    } else if det.is_negative() {
                        -1
                    } else {
                        0
                    };
                    prop_assert_eq!(sign, expected);
                }
            }
        }
    };
}

gen_det_sign_agrees_with_det_exact_proptests!(2);
gen_det_sign_agrees_with_det_exact_proptests!(3);
gen_det_sign_agrees_with_det_exact_proptests!(4);
gen_det_sign_agrees_with_det_exact_proptests!(5);

/// Fast-filter invariant: whenever `|det_direct()| > det_errbound()`,
/// the f64 sign is provably correct — so
/// `det_direct().signum() == det_sign_exact()`.  This is the
/// correctness guarantee the Shewchuk-style filter inside
/// `det_sign_exact` relies on.  The proptest cross-checks that the
/// fast-filter boundary itself is honoured, independent of whether
/// `det_sign_exact` ended up using the filter or the Bareiss fallback
/// on any particular input.  Only D=2..=4 have a closed-form
/// `det_direct` / `det_errbound` pair.
macro_rules! gen_det_sign_fast_filter_boundary_proptests {
    ($d:literal) => {
        paste! {
            proptest! {
                #![proptest_config(ProptestConfig::with_cases(64))]

                #[test]
                fn [<det_sign_exact_agrees_with_det_direct_when_filter_conclusive_ $d d>](
                    entries in proptest::array::[<uniform $d>](
                        proptest::array::[<uniform $d>](small_int_f64()),
                    ),
                ) {
                    let m = Matrix::<$d>::from_rows(entries);
                    let det = m.det_direct().expect("D<=4 has closed-form det_direct");
                    let bound = m.det_errbound().expect("D<=4 has a det_errbound");
                    let sign = m.det_sign_exact().unwrap();

                    // Only assert when the filter is conclusive.  When
                    // `|det| <= bound` the f64 sign may disagree with the
                    // exact sign; that case is covered by the other
                    // proptests via the Bareiss fallback.
                    if det.abs() > bound {
                        let direct_sign: i8 = if det > 0.0 {
                            1
                        } else if det < 0.0 {
                            -1
                        } else {
                            0
                        };
                        prop_assert_eq!(direct_sign, sign);
                    }
                }
            }
        }
    };
}

gen_det_sign_fast_filter_boundary_proptests!(2);
gen_det_sign_fast_filter_boundary_proptests!(3);
gen_det_sign_fast_filter_boundary_proptests!(4);
