#![forbid(unsafe_code)]

//! Property-based tests for the exact-arithmetic APIs
//! (requires `exact` feature).
//!
//! Covers:
//! - `det_sign_exact` on diagonal and full small-integer matrices
//! - `det_exact` on full small-integer matrices against an independent
//!   `BigRational` Leibniz-expansion oracle
//! - determinant sign and error-bound filtering across independently mixed
//!   binary64 exponent regimes
//! - `solve_exact` round-trip with integer inputs (`A · x0` in f64 is
//!   exact for small integers, so `solve(A, A · x0) == x0`)
//! - `solve_exact` residual property (`A · solve(A, b) == b` in
//!   `BigRational` arithmetic) on random integer and mixed-exponent inputs

#![cfg(feature = "exact")]

use std::array::from_fn;

use pastey::paste;
use proptest::{array, prelude::*};

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

/// Construct the exactly representable finite binary64 value `2^exponent`.
///
/// Building subnormal powers from their bit representation avoids the
/// intermediate underflow that `f64::powi` can incur for exponents below
/// -1023 on some targets.
fn exact_power_of_two(exponent: i32) -> f64 {
    assert!(
        (-1074..=1023).contains(&exponent),
        "binary64 power-of-two exponent must be finite"
    );
    if exponent < -1022 {
        f64::from_bits(1_u64 << (exponent + 1074).cast_unsigned())
    } else {
        let biased_exponent =
            u64::try_from(exponent + 1023).expect("normal exponent bias is non-negative");
        f64::from_bits(biased_exponent << 52)
    }
}

fn mixed_scale_finite_f64() -> impl Strategy<Value = f64> {
    prop_oneof![
        Just(0.0),
        Just(-0.0),
        Just(f64::from_bits(1)),
        Just(-f64::from_bits(1)),
        Just(f64::MIN_POSITIVE),
        Just(-f64::MIN_POSITIVE),
        Just(1.0),
        Just(-1.0),
        Just(3.5),
        Just(-7.25),
        Just(f64::MAX / 4.0),
        Just(-f64::MAX / 4.0),
    ]
}

/// Finite binary64 values whose exponent regime is selected independently for
/// each generated entry.
///
/// Non-zero values are a small integer coefficient times an exact power of
/// two. The ranges deliberately cover subnormal, tiny normal, ordinary, and
/// large finite values without approaching overflow during value generation.
fn mixed_exponent_finite_f64() -> impl Strategy<Value = f64> {
    prop_oneof![
        1 => Just(0.0),
        1 => Just(-0.0),
        3 => (small_nonzero_int_f64(), -1074i32..=-1023i32)
            .prop_map(|(coefficient, exponent)| coefficient * exact_power_of_two(exponent)),
        3 => (small_nonzero_int_f64(), -1022i32..=-900i32)
            .prop_map(|(coefficient, exponent)| coefficient * exact_power_of_two(exponent)),
        4 => (small_nonzero_int_f64(), -20i32..=20i32)
            .prop_map(|(coefficient, exponent)| coefficient * exact_power_of_two(exponent)),
        3 => (small_nonzero_int_f64(), 900i32..=1018i32)
            .prop_map(|(coefficient, exponent)| coefficient * exact_power_of_two(exponent)),
    ]
}

/// Independent exact power-of-two row scales used by the dense solve corpus.
///
/// Exponents leave ample upper-range headroom for the small integer diagonal
/// while spanning subnormal, tiny normal, ordinary, and large finite rows.
fn solve_row_exponent() -> impl Strategy<Value = i32> {
    prop_oneof![
        2 => -1074i32..=-1022i32,
        3 => -900i32..=-300i32,
        3 => -20i32..=20i32,
        3 => 300i32..=900i32,
    ]
}

fn is_unrepresentable<T>(
    result: &Result<T, LaError>,
    expected_index: Option<usize>,
    expected_reason: UnrepresentableReason,
) -> bool {
    matches!(
        result,
        Err(LaError::Unrepresentable { index, reason, .. })
            if *index == expected_index && *reason == expected_reason
    )
}

/// Multiply `A · x` entirely in `BigRational`, lifting each f64 matrix
/// entry via `BigRational::from_f64`.  Used by residual assertions.
///
/// Every accepted finite f64 has an exact rational reconstruction, including
/// subnormal values and values with large binary exponents.
fn big_rational_matvec<const D: usize>(
    a: &[[f64; D]; D],
    x: &[BigRational; D],
) -> [BigRational; D] {
    from_fn(|i| {
        let mut sum = BigRational::from_integer(BigInt::from(0));
        for (aij, xj) in a[i].iter().zip(x.iter()) {
            let entry = BigRational::from_f64(*aij).expect("finite f64 converts exactly");
            sum += entry * xj;
        }
        sum
    })
}

/// Compute an exact determinant via the Leibniz permutation expansion.
///
/// This is intentionally independent from the production Bareiss core. It is
/// factorial-time, but the proptests only use D=2..=5, so it stays tiny while
/// giving `det_exact` a separate dense-matrix oracle.
fn big_rational_det_leibniz<const D: usize>(a: &[[f64; D]; D]) -> BigRational {
    let mut det = BigRational::from_integer(BigInt::from(0));
    let mut perm: [usize; D] = from_fn(|i| i);

    loop {
        let mut term = BigRational::from_integer(BigInt::from(1));
        for (row, &col) in perm.iter().enumerate() {
            let entry = BigRational::from_f64(a[row][col]).expect("finite f64 converts exactly");
            term *= entry;
        }

        if permutation_is_even(&perm) {
            det += term;
        } else {
            det -= term;
        }

        if !next_permutation(&mut perm) {
            break;
        }
    }

    det
}

fn determinant_sign(value: &BigRational) -> DeterminantSign {
    if value.is_positive() {
        DeterminantSign::Positive
    } else if value.is_negative() {
        DeterminantSign::Negative
    } else {
        DeterminantSign::Zero
    }
}

fn permutation_is_even(perm: &[usize]) -> bool {
    let mut inversions = 0usize;
    for i in 0..perm.len() {
        for j in (i + 1)..perm.len() {
            if perm[i] > perm[j] {
                inversions += 1;
            }
        }
    }
    inversions.is_multiple_of(2)
}

fn next_permutation(values: &mut [usize]) -> bool {
    if values.len() < 2 {
        return false;
    }

    let mut pivot = values.len() - 2;
    loop {
        if values[pivot] < values[pivot + 1] {
            break;
        }
        if pivot == 0 {
            return false;
        }
        pivot -= 1;
    }

    let mut successor = values.len() - 1;
    while values[successor] <= values[pivot] {
        successor -= 1;
    }
    values.swap(pivot, successor);
    values[(pivot + 1)..].reverse();
    true
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
    let dimension = u32::try_from(D).expect("proptest matrix dimension must fit in u32");
    let shift = f64::from(dimension).mul_add(10.0, 1.0);
    for i in 0..D {
        rows[i][i] = if diag[i] >= 0.0 {
            diag[i] + shift
        } else {
            diag[i] - shift
        };
    }
    rows
}

/// Build a dense, strictly diagonally-dominant matrix, then scale each row by
/// an independently generated exact power of two.
///
/// Every off-diagonal entry is a non-zero small integer. Before scaling, the
/// diagonal is one greater than the row's off-diagonal absolute sum, so
/// `|A[i][i]| > Σ_{j≠i} |A[i][j]|`. Positive row scaling preserves that
/// inequality and therefore invertibility (Levy–Desplanques). The selected
/// exponent ranges keep every scaled entry finite and exactly representable,
/// including at the subnormal boundary.
fn make_dense_mixed_exponent_matrix<const D: usize>(
    offdiag: [[f64; D]; D],
    row_exponents: [i32; D],
) -> [[f64; D]; D] {
    let mut rows = offdiag;
    for i in 0..D {
        let offdiag_abs_sum = rows[i]
            .iter()
            .enumerate()
            .filter(|&(j, _)| j != i)
            .map(|(_, value)| value.abs())
            .sum::<f64>();
        rows[i][i] = offdiag_abs_sum + 1.0;

        let row_scale = exact_power_of_two(row_exponents[i]);
        for value in &mut rows[i] {
            *value *= row_scale;
        }
    }
    rows
}

#[test]
fn solve_exact_handles_bit_exact_subnormal_row_scales() {
    let row_0_scale = exact_power_of_two(-1023);
    let row_1_scale = exact_power_of_two(-1024);
    assert_eq!(exact_power_of_two(-1074).to_bits(), 1);
    assert_eq!(row_1_scale.to_bits(), 1_u64 << 50);
    assert_eq!(row_0_scale.to_bits(), 1_u64 << 51);

    let rows = [
        [3.0 * row_0_scale, -2.0 * row_0_scale],
        [-2.0 * row_1_scale, 3.0 * row_1_scale],
    ];
    let determinant = big_rational_det_leibniz(&rows);
    let expected_determinant = BigRational::new(BigInt::from(5), BigInt::from(1_u8) << 2047_u32);
    assert_eq!(determinant, expected_determinant);
    assert!(determinant.is_positive());

    let matrix = Matrix::<2>::try_from_rows(rows).unwrap();
    let rhs = Vector::<2>::try_new([row_0_scale, row_1_scale]).unwrap();
    let solution = matrix.solve_exact(rhs).unwrap();
    let one = BigRational::from_integer(BigInt::from(1));
    assert_eq!(solution, [one.clone(), one]);

    let residual = big_rational_matvec(&rows, &solution);
    assert_eq!(
        residual,
        [
            BigRational::from_f64(row_0_scale).unwrap(),
            BigRational::from_f64(row_1_scale).unwrap(),
        ]
    );
}

macro_rules! gen_det_sign_exact_proptests {
    ($d:literal) => {
        paste! {
            proptest! {
                #![proptest_config(ProptestConfig::with_cases(64))]

                #[test]
                fn [<det_sign_exact_agrees_with_diagonal_product_and_det_ $d d>](
                    diag in array::[<uniform $d>](small_nonzero_f64()),
                ) {
                    // Diagonal matrix: determinant sign = product of diagonal signs.
                    let mut rows = [[0.0f64; $d]; $d];
                    for i in 0..$d {
                        rows[i][i] = diag[i];
                    }
                    let m = Matrix::<$d>::try_from_rows(rows).unwrap();

                    let exact_sign = m.det_sign_exact();

                    // Expected sign from the product of diagonal entries.
                    let neg_count = diag.iter().filter(|&&x| x < 0.0).count();
                    let expected_sign = if neg_count % 2 == 0 {
                        DeterminantSign::Positive
                    } else {
                        DeterminantSign::Negative
                    };

                    prop_assert_eq!(exact_sign, expected_sign);
                    let fp_det = m.det().unwrap();
                    let fp_sign = if fp_det > 0.0 {
                        DeterminantSign::Positive
                    } else if fp_det < 0.0 {
                        DeterminantSign::Negative
                    } else {
                        DeterminantSign::Zero
                    };

                    prop_assert_eq!(fp_sign, expected_sign);
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
                    offdiag in array::[<uniform $d>](
                        array::[<uniform $d>](small_int_f64()),
                    ),
                    diag in array::[<uniform $d>](small_nonzero_int_f64()),
                    x0 in array::[<uniform $d>](small_int_f64()),
                ) {
                    let rows = make_diagonally_dominant::<$d>(offdiag, diag);
                    let a = Matrix::<$d>::try_from_rows(rows).unwrap();

                    // b = A · x0, computed in f64.  Small integers keep
                    // every partial sum exact.
                    let mut b_arr = [0.0f64; $d];
                    for i in 0..$d {
                        let mut sum = 0.0f64;
                        for j in 0..$d {
                            sum = rows[i][j].mul_add(x0[j], sum);
                        }
                        b_arr[i] = sum;
                    }
                    let b = Vector::<$d>::try_new(b_arr).unwrap();
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
                    offdiag in array::[<uniform $d>](
                        array::[<uniform $d>](small_int_f64()),
                    ),
                    diag in array::[<uniform $d>](small_nonzero_int_f64()),
                    b_arr in array::[<uniform $d>](small_int_f64()),
                ) {
                    let rows = make_diagonally_dominant::<$d>(offdiag, diag);
                    let a = Matrix::<$d>::try_from_rows(rows).unwrap();
                    let b = Vector::<$d>::try_new(b_arr).unwrap();
                    let x = a.solve_exact(b).expect("diagonally-dominant A is non-singular");

                    let ax = big_rational_matvec::<$d>(&rows, &x);
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

/// Mixed-exponent residual property: dense matrices remain exactly solvable
/// when every row has an independent power-of-two scale and every RHS entry
/// independently ranges from zero/subnormal through large finite values.
///
/// The matrix constructor guarantees strict diagonal dominance. The residual
/// oracle independently reconstructs the original f64 inputs as rationals and
/// verifies `A · solve_exact(A, b) == b` without reusing the Bareiss core.
macro_rules! gen_solve_exact_mixed_exponent_residual_proptests {
    ($d:literal) => {
        paste! {
            proptest! {
                #![proptest_config(ProptestConfig::with_cases(24))]

                #[test]
                fn [<solve_exact_mixed_exponent_residual_ $d d>](
                    offdiag in array::[<uniform $d>](
                        array::[<uniform $d>](small_nonzero_int_f64()),
                    ),
                    row_exponents in array::[<uniform $d>](solve_row_exponent()),
                    b_arr in array::[<uniform $d>](mixed_exponent_finite_f64()),
                ) {
                    let rows = make_dense_mixed_exponent_matrix::<$d>(offdiag, row_exponents);
                    let a = Matrix::<$d>::try_from_rows(rows).unwrap();
                    let b = Vector::<$d>::try_new(b_arr).unwrap();
                    let x = a
                        .solve_exact(b)
                        .expect("strict diagonal dominance guarantees invertibility");

                    let ax = big_rational_matvec::<$d>(&rows, &x);
                    for i in 0..$d {
                        let b_rat = BigRational::from_f64(b_arr[i])
                            .expect("finite f64 converts exactly");
                        prop_assert_eq!(&ax[i], &b_rat);
                    }
                }
            }
        }
    };
}

gen_solve_exact_mixed_exponent_residual_proptests!(2);
gen_solve_exact_mixed_exponent_residual_proptests!(3);
gen_solve_exact_mixed_exponent_residual_proptests!(4);
gen_solve_exact_mixed_exponent_residual_proptests!(5);

/// Dense determinant oracle: both the exact value and adaptive sign must match
/// one independent `BigRational` Leibniz expansion.
macro_rules! gen_det_exact_and_sign_leibniz_oracle_proptests {
    ($d:literal) => {
        paste! {
            proptest! {
                #![proptest_config(ProptestConfig::with_cases(64))]

                #[test]
                fn [<det_exact_and_sign_agree_with_leibniz_oracle_ $d d>](
                    entries in array::[<uniform $d>](
                        array::[<uniform $d>](small_int_f64()),
                    ),
                ) {
                    let m = Matrix::<$d>::try_from_rows(entries).unwrap();
                    let expected = big_rational_det_leibniz::<$d>(&entries);
                    let expected_sign = determinant_sign(&expected);

                    prop_assert_eq!(m.det_exact().unwrap(), expected);
                    prop_assert_eq!(m.det_sign_exact(), expected_sign);
                }
            }
        }
    };
}

gen_det_exact_and_sign_leibniz_oracle_proptests!(2);
gen_det_exact_and_sign_leibniz_oracle_proptests!(3);
gen_det_exact_and_sign_leibniz_oracle_proptests!(4);
gen_det_exact_and_sign_leibniz_oracle_proptests!(5);

/// Fast-filter invariant: whenever `|det_direct()| > det_errbound()`,
/// the f64 sign is provably correct. The expected sign comes directly from an
/// independent `BigRational` Leibniz expansion rather than `det_sign_exact`,
/// avoiding a self-referential comparison with the filter's own consumer.
/// Only D=2..=4 have a closed-form `det_direct` / `det_errbound` pair.
macro_rules! gen_det_sign_fast_filter_boundary_proptests {
    ($d:literal) => {
        paste! {
            proptest! {
                #![proptest_config(ProptestConfig::with_cases(64))]

                #[test]
                fn [<det_direct_sign_agrees_with_leibniz_when_filter_conclusive_ $d d>](
                    entries in array::[<uniform $d>](
                        array::[<uniform $d>](small_int_f64()),
                    ),
                ) {
                    let m = Matrix::<$d>::try_from_rows(entries).unwrap();
                    let det = m
                        .det_direct()
                        .unwrap()
                        .expect("D<=4 has closed-form det_direct");
                    let exact = big_rational_det_leibniz::<$d>(&entries);
                    let exact_sign = determinant_sign(&exact);

                    // Only assert when the filter is conclusive.  When
                    // `det_errbound` is unavailable or `|det| <= bound`, the
                    // f64 sign may disagree with the exact sign; those cases
                    // fall through to direct exact-integer evaluation.
                    if let Some(bound) = m.det_errbound().unwrap() {
                        if det.abs() > bound {
                            let direct_sign = if det > 0.0 {
                                DeterminantSign::Positive
                            } else if det < 0.0 {
                                DeterminantSign::Negative
                            } else {
                                DeterminantSign::Zero
                            };
                            prop_assert_eq!(direct_sign, exact_sign);
                        }
                    }
                }
            }
        }
    };
}

gen_det_sign_fast_filter_boundary_proptests!(2);
gen_det_sign_fast_filter_boundary_proptests!(3);
gen_det_sign_fast_filter_boundary_proptests!(4);

/// Error-bound invariant: for every dense D≤4 matrix in this corpus,
/// `det_errbound()` must bound the absolute error of `det_direct()` against an
/// independent exact Leibniz expansion.  The entries include decimal fractions
/// that are not generally exactly representable in binary64, so this exercises
/// non-trivial rounding in the closed-form determinant path while keeping the
/// magnitudes far from overflow.
macro_rules! gen_det_errbound_leibniz_oracle_proptests {
    ($d:literal) => {
        paste! {
            proptest! {
                #![proptest_config(ProptestConfig::with_cases(64))]

                #[test]
                fn [<det_errbound_bounds_det_direct_error_ $d d>](
                    entries in array::[<uniform $d>](
                        array::[<uniform $d>](
                            (-50i16..=50i16).prop_map(|x| f64::from(x) / 10.0)
                        ),
                    ),
                ) {
                    let m = Matrix::<$d>::try_from_rows(entries).unwrap();
                    let det_direct = m
                        .det_direct()
                        .unwrap()
                        .expect("D<=4 has closed-form det_direct");
                    let exact = big_rational_det_leibniz::<$d>(&entries);
                    if let Some(bound) = m.det_errbound().unwrap() {
                        let direct_exact = BigRational::from_f64(det_direct)
                            .expect("det_direct returned finite f64");
                        let bound_exact = BigRational::from_f64(bound)
                            .expect("det_errbound returned finite f64");
                        let error = (direct_exact - exact).abs();

                        prop_assert!(
                            error <= bound_exact,
                            "det_direct error exceeded det_errbound for D={}: error={error}, bound={bound_exact}",
                            $d
                        );
                    }
                }
            }
        }
    };
}

gen_det_errbound_leibniz_oracle_proptests!(2);
gen_det_errbound_leibniz_oracle_proptests!(3);
gen_det_errbound_leibniz_oracle_proptests!(4);

/// Exercise the determinant certificate with independently mixed per-entry
/// exponents spanning zero, subnormal, tiny normal, ordinary, and large finite
/// regimes. `det_sign_exact` must always match the independent Leibniz oracle;
/// whenever the filter publishes a bound, the same oracle must confirm it.
/// Inconclusive or overflowed scalar filter arithmetic defers to direct
/// exact-integer evaluation for these D≤4 cases.
macro_rules! gen_extreme_exponent_det_filter_proptests {
    ($d:literal) => {
        paste! {
            proptest! {
                #![proptest_config(ProptestConfig::with_cases(32))]

                #[test]
                fn [<det_filter_is_conservative_across_mixed_exponents_ $d d>](
                    entries in array::[<uniform $d>](
                        array::[<uniform $d>](mixed_exponent_finite_f64()),
                    ),
                ) {
                    let matrix = Matrix::<$d>::try_from_rows(entries).unwrap();
                    let exact = big_rational_det_leibniz::<$d>(&entries);
                    let expected_sign = determinant_sign(&exact);

                    prop_assert_eq!(matrix.det_sign_exact(), expected_sign);

                    if let Ok(Some(bound)) = matrix.det_errbound() {
                        let direct = matrix
                            .det_direct()
                            .unwrap()
                            .expect("D<=4 has closed-form det_direct");
                        let direct_exact = BigRational::from_f64(direct)
                            .expect("det_direct returned finite f64");
                        let bound_exact = BigRational::from_f64(bound)
                            .expect("det_errbound returned finite f64");
                        let error = (direct_exact - exact).abs();

                        prop_assert!(
                            error <= bound_exact,
                            "mixed-exponent determinant error exceeded bound for D={}: error={error}, bound={bound_exact}",
                            $d,
                        );
                    }
                }
            }
        }
    };
}

gen_extreme_exponent_det_filter_proptests!(2);
gen_extreme_exponent_det_filter_proptests!(3);
gen_extreme_exponent_det_filter_proptests!(4);

/// Mixed-scale diagonal matrices stress the shared-exponent conversion path:
/// zeros, subnormals, tiny normal values, ordinary values, and very large
/// finite values can all appear in the same determinant.  The independent
/// expectation uses `BigRational::from_f64` on each diagonal value.
macro_rules! gen_mixed_scale_diagonal_exact_det_proptests {
    ($d:literal) => {
        paste! {
            proptest! {
                #![proptest_config(ProptestConfig::with_cases(32))]

                #[test]
                fn [<det_exact_handles_mixed_scale_diagonal_ $d d>](
                    diag in array::[<uniform $d>](mixed_scale_finite_f64()),
                ) {
                    let mut rows = [[0.0f64; $d]; $d];
                    let mut expected = BigRational::from_integer(BigInt::from(1));
                    for i in 0..$d {
                        rows[i][i] = diag[i];
                        expected *= BigRational::from_f64(diag[i])
                            .expect("strategy only emits finite f64 values");
                    }
                    let m = Matrix::<$d>::try_from_rows(rows).unwrap();

                    let expected_sign = determinant_sign(&expected);
                    let expected_f64 = expected.to_f64();

                    prop_assert_eq!(m.det_sign_exact(), expected_sign);

                    match expected_f64 {
                        Some(expected_f64)
                            if expected_f64.is_finite()
                                && BigRational::from_f64(expected_f64).as_ref()
                                    == Some(&expected) =>
                        {
                            prop_assert_eq!(
                                m.det_exact_f64().unwrap().to_bits(),
                                expected_f64.to_bits()
                            );
                        }
                        Some(expected_f64) if expected_f64.is_finite() => {
                            let result = m.det_exact_f64();
                            prop_assert!(is_unrepresentable(
                                &result,
                                None,
                                UnrepresentableReason::RequiresRounding
                            ));
                        }
                        _ => {
                            let result = m.det_exact_f64();
                            prop_assert!(is_unrepresentable(
                                &result,
                                None,
                                UnrepresentableReason::NotFinite
                            ));
                        }
                    }

                    prop_assert_eq!(m.det_exact().unwrap(), expected);
                }
            }
        }
    };
}

gen_mixed_scale_diagonal_exact_det_proptests!(2);
gen_mixed_scale_diagonal_exact_det_proptests!(3);
gen_mixed_scale_diagonal_exact_det_proptests!(4);
gen_mixed_scale_diagonal_exact_det_proptests!(5);
