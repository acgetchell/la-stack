#![forbid(unsafe_code)]

//! Benchmarks for exact arithmetic operations.
//!
//! These benchmarks measure the performance of the `exact` feature's
//! arbitrary-precision methods.  They are organised into four classes:
//!
//! 1. **General-case benches** (`exact_d{2..5}`) — a single
//!    well-conditioned diagonally-dominant matrix per dimension.  These
//!    measure typical-case performance and track regressions against a
//!    reproducible input. D=2..=4 also include the direct determinant and
//!    certified error-bound f64 baselines.
//! 2. **Adversarial / extreme-input benches** — matrices chosen to
//!    stress specific corners of the exact-arithmetic pipeline:
//!    near-singularity (forces the exact integer fallback), large f64 entries
//!    (stresses intermediate `BigInt` growth), and Hilbert-style
//!    ill-conditioning (wide range of `(mantissa, exponent)` pairs in
//!    the `decompose_f64 → BigInt` path).  These measure tail behaviour
//!    that fixed well-conditioned inputs miss and provide stronger
//!    empirical evidence for `docs/PERFORMANCE.md`.
//! 3. **Random corpus benches** (`exact_random_corpus_d{2..5}`) — a
//!    fixed-seed corpus of diagonally-dominant random matrices per dimension.
//!    Every measured iteration executes the full corpus in its stable order,
//!    so current and baseline revisions receive identical workloads.
//! 4. **Exact-input rational benches** (`rational_input_d{2..8}`) — compare
//!    row-denominator clearing plus integer Bareiss elimination with
//!    straightforward cubic `BigRational` Gaussian elimination on identical
//!    rational matrices and right-hand sides.
//!
//! Fallible exact-to-f64 conversions use a `_result` suffix. Those rows measure
//! the full `Result` path, including valid `Err(Unrepresentable)` outcomes for
//! inputs whose exact answer cannot be represented as finite binary64.

use std::hint::black_box;

use criterion::{BenchmarkGroup, Criterion, Throughput, measurement::WallTime};

#[cfg(not(any(la_stack_pre_rational_input_api, la_stack_v0_4_3_api)))]
use la_stack::{BigInt, BigRational, DeterminantSign, RationalMatrix, RationalVector};
use la_stack::{Matrix, Vector};

#[path = "common/bench_utils.rs"]
mod bench_utils;
#[path = "common/exact.rs"]
pub mod exact_bench;

use bench_utils::OrAbort;
use exact_bench::{
    ExactInput, RANDOM_INPUT_ARRAY_LEN, ValidatedExactInput, hilbert_input,
    large_entries_3x3_input, make_matrix_rows, make_random_input_corpus, make_vector_array,
    near_singular_3x3_input, validate_exact_fixture, validate_f64_determinant_benchmarks,
};

/// Exact operation measured by a benchmark group.
#[derive(Clone, Copy)]
enum ExactOperation {
    DetSignExact,
    DetExact,
    DetExactF64Result,
    DetExactRoundedF64,
    SolveExact,
    SolveExactF64Result,
    SolveExactRoundedF64,
}

impl ExactOperation {
    /// Return the benchmark-name stem for this exact operation.
    const fn name(self) -> &'static str {
        match self {
            Self::DetSignExact => "det_sign_exact",
            Self::DetExact => "det_exact",
            Self::DetExactF64Result => "det_exact_f64_result",
            Self::DetExactRoundedF64 => "det_exact_rounded_f64",
            Self::SolveExact => "solve_exact",
            Self::SolveExactF64Result => "solve_exact_f64_result",
            Self::SolveExactRoundedF64 => "solve_exact_rounded_f64",
        }
    }
}

const GENERAL_OPERATIONS: &[ExactOperation] = &[
    ExactOperation::DetExact,
    ExactOperation::DetExactF64Result,
    ExactOperation::DetExactRoundedF64,
    ExactOperation::DetSignExact,
    ExactOperation::SolveExact,
    ExactOperation::SolveExactF64Result,
    ExactOperation::SolveExactRoundedF64,
];

const CORPUS_AND_EXTREME_OPERATIONS: &[ExactOperation] = &[
    ExactOperation::DetSignExact,
    ExactOperation::DetExact,
    ExactOperation::SolveExact,
    ExactOperation::SolveExactF64Result,
    ExactOperation::SolveExactRoundedF64,
];

/// Independently validated exact-input rational benchmark fixture.
#[cfg(not(any(la_stack_pre_rational_input_api, la_stack_v0_4_3_api)))]
#[must_use]
struct ValidatedRationalInput<const D: usize> {
    matrix: RationalMatrix<D>,
    rhs: RationalVector<D>,
}

/// Deterministic, strictly diagonally-dominant exact-input rational fixture.
#[cfg(not(any(la_stack_pre_rational_input_api, la_stack_v0_4_3_api)))]
fn rational_input<const D: usize>() -> ValidatedRationalInput<D> {
    let rows = std::array::from_fn(|row| {
        std::array::from_fn(|col| {
            if row == col {
                let diagonal = 2 * D + row + 1;
                BigRational::from_integer(BigInt::from(diagonal))
            } else {
                let raw_numerator = (row * 3 + col * 5) % 5;
                let numerator = i64::try_from(raw_numerator).or_abort("rational numerator") - 2;
                let denominator = (row + col) % 7 + 2;
                BigRational::new(BigInt::from(numerator), BigInt::from(denominator))
            }
        })
    });
    let expected_solution = std::array::from_fn(|index| {
        BigRational::new(BigInt::from(index + 1), BigInt::from(index + 2))
    });
    let rhs_data = rational_matvec(&rows, &expected_solution);
    let matrix = RationalMatrix::try_from_rows(rows.clone())
        .or_abort("rational benchmark matrix construction");
    let rhs =
        RationalVector::try_new(rhs_data.clone()).or_abort("rational benchmark RHS construction");

    let reference_determinant = rational_determinant_gaussian(rows.clone());
    assert_eq!(matrix.det(), reference_determinant);
    assert_eq!(matrix.det_sign(), determinant_sign(&reference_determinant));
    let reference_solution = rational_solve_gaussian(rows, rhs_data)
        .or_abort("rational Gaussian benchmark validation solve");
    assert_eq!(reference_solution, expected_solution);
    assert_eq!(
        matrix
            .solve(&rhs)
            .or_abort("row-cleared Bareiss benchmark validation solve")
            .into_array(),
        expected_solution
    );

    ValidatedRationalInput { matrix, rhs }
}

/// Exact rational matrix-vector multiplication used only to assemble and
/// validate benchmark fixtures.
#[cfg(not(any(la_stack_pre_rational_input_api, la_stack_v0_4_3_api)))]
fn rational_matvec<const D: usize>(
    rows: &[[BigRational; D]; D],
    vector: &[BigRational; D],
) -> [BigRational; D] {
    std::array::from_fn(|row| {
        rows[row]
            .iter()
            .zip(vector.iter())
            .map(|(coefficient, component)| coefficient * component)
            .sum()
    })
}

/// Straightforward cubic rational Gaussian determinant reference.
#[cfg(not(any(la_stack_pre_rational_input_api, la_stack_v0_4_3_api)))]
fn rational_determinant_gaussian<const D: usize>(mut rows: [[BigRational; D]; D]) -> BigRational {
    let zero = BigRational::from_integer(BigInt::from(0));
    let mut determinant = BigRational::from_integer(BigInt::from(1));
    let mut odd_swaps = false;

    for pivot_col in 0..D {
        let Some(pivot_row) = (pivot_col..D).find(|&row| rows[row][pivot_col] != zero) else {
            return zero;
        };
        if pivot_row != pivot_col {
            rows.swap(pivot_col, pivot_row);
            odd_swaps = !odd_swaps;
        }

        let pivot = rows[pivot_col][pivot_col].clone();
        let pivot_entries = rows[pivot_col].clone();
        determinant *= &pivot;
        for row_entries in rows.iter_mut().skip(pivot_col + 1) {
            let factor = &row_entries[pivot_col] / &pivot;
            for (entry, pivot_entry) in row_entries
                .iter_mut()
                .zip(pivot_entries.iter())
                .skip(pivot_col + 1)
            {
                *entry -= &factor * pivot_entry;
            }
            row_entries[pivot_col] = zero.clone();
        }
    }

    if odd_swaps { -determinant } else { determinant }
}

/// Straightforward cubic rational Gaussian solve reference.
#[cfg(not(any(la_stack_pre_rational_input_api, la_stack_v0_4_3_api)))]
fn rational_solve_gaussian<const D: usize>(
    mut rows: [[BigRational; D]; D],
    mut rhs: [BigRational; D],
) -> Option<[BigRational; D]> {
    let zero = BigRational::from_integer(BigInt::from(0));
    for pivot_col in 0..D {
        let pivot_row = (pivot_col..D).find(|&row| rows[row][pivot_col] != zero)?;
        if pivot_row != pivot_col {
            rows.swap(pivot_col, pivot_row);
            rhs.swap(pivot_col, pivot_row);
        }

        let pivot = rows[pivot_col][pivot_col].clone();
        let pivot_entries = rows[pivot_col].clone();
        let pivot_rhs = rhs[pivot_col].clone();
        for (row_entries, rhs_entry) in rows.iter_mut().zip(rhs.iter_mut()).skip(pivot_col + 1) {
            let factor = &row_entries[pivot_col] / &pivot;
            for (entry, pivot_entry) in row_entries
                .iter_mut()
                .zip(pivot_entries.iter())
                .skip(pivot_col + 1)
            {
                *entry -= &factor * pivot_entry;
            }
            let rhs_update = &factor * &pivot_rhs;
            *rhs_entry -= rhs_update;
            row_entries[pivot_col] = zero.clone();
        }
    }

    let mut solution = std::array::from_fn(|_| zero.clone());
    for row in (0..D).rev() {
        let mut value = rhs[row].clone();
        for (coefficient, component) in rows[row].iter().zip(solution.iter()).skip(row + 1) {
            value -= coefficient * component;
        }
        solution[row] = value / &rows[row][row];
    }
    Some(solution)
}

#[cfg(not(any(la_stack_pre_rational_input_api, la_stack_v0_4_3_api)))]
fn determinant_sign(value: &BigRational) -> DeterminantSign {
    let zero = BigRational::from_integer(BigInt::from(0));
    match value.cmp(&zero) {
        std::cmp::Ordering::Less => DeterminantSign::Negative,
        std::cmp::Ordering::Equal => DeterminantSign::Zero,
        std::cmp::Ordering::Greater => DeterminantSign::Positive,
    }
}

/// Compare exact-input row clearing and Bareiss elimination with direct
/// `BigRational` Gaussian elimination.
#[cfg(not(any(la_stack_pre_rational_input_api, la_stack_v0_4_3_api)))]
fn bench_rational_input<const D: usize>(criterion: &mut Criterion) {
    let input = rational_input::<D>();
    let mut group = criterion.benchmark_group(format!("rational_input_d{D}"));

    group.bench_function("det_sign_row_cleared_bareiss", |bencher| {
        bencher.iter(|| {
            let sign = black_box(&input.matrix).det_sign();
            let _ = black_box(sign);
        });
    });
    group.bench_function("det_row_cleared_bareiss", |bencher| {
        bencher.iter(|| {
            let determinant = black_box(&input.matrix).det();
            black_box(determinant);
        });
    });
    group.bench_function("det_big_rational_gaussian", |bencher| {
        bencher.iter(|| {
            let rows = black_box(input.matrix.as_rows()).clone();
            let determinant = rational_determinant_gaussian(rows);
            black_box(determinant);
        });
    });
    group.bench_function("solve_row_cleared_bareiss", |bencher| {
        bencher.iter(|| {
            let solution = black_box(&input.matrix)
                .solve(black_box(&input.rhs))
                .or_abort("row-cleared Bareiss rational benchmark solve");
            let _ = black_box(solution);
        });
    });
    group.bench_function("solve_big_rational_gaussian", |bencher| {
        bencher.iter(|| {
            let rows = black_box(input.matrix.as_rows()).clone();
            let rhs = black_box(input.rhs.as_array()).clone();
            let solution =
                rational_solve_gaussian(rows, rhs).or_abort("BigRational Gaussian benchmark solve");
            black_box(solution);
        });
    });

    group.finish();
}

/// Execute one exact operation on a borrowed, independently validated input.
fn run_exact_operation<const D: usize>(operation: ExactOperation, input: &ValidatedExactInput<D>) {
    match operation {
        ExactOperation::DetSignExact => {
            let sign = black_box(input.matrix()).det_sign_exact();
            let _ = black_box(sign);
        }
        ExactOperation::DetExact => {
            let det = black_box(input.matrix())
                .det_exact()
                .or_abort("exact determinant");
            black_box(det);
        }
        ExactOperation::DetExactF64Result => {
            let det = black_box(input.matrix()).det_exact_f64();
            let _ = black_box(det);
        }
        ExactOperation::DetExactRoundedF64 => {
            let det = black_box(input.matrix())
                .det_exact_rounded_f64()
                .or_abort("exact determinant rounded to f64");
            let _ = black_box(det);
        }
        ExactOperation::SolveExact => {
            let x = black_box(input.matrix())
                .solve_exact(black_box(input.rhs()))
                .or_abort("exact linear solve");
            let _ = black_box(x);
        }
        ExactOperation::SolveExactF64Result => {
            let x = black_box(input.matrix()).solve_exact_f64(black_box(input.rhs()));
            let _ = black_box(x);
        }
        ExactOperation::SolveExactRoundedF64 => {
            let x = black_box(input.matrix())
                .solve_exact_rounded_f64(black_box(input.rhs()))
                .or_abort("exact linear solve rounded to f64");
            let _ = black_box(x);
        }
    }
}

/// Add one exact-arithmetic operation benchmark over a validated fixed input pair.
fn bench_exact_operation<const D: usize>(
    group: &mut BenchmarkGroup<'_, WallTime>,
    operation: ExactOperation,
    input: &ValidatedExactInput<D>,
) {
    group.bench_function(operation.name(), |bencher| {
        bencher.iter(|| {
            run_exact_operation(operation, input);
        });
    });
}

/// Add one Criterion benchmark that executes the complete validated random corpus.
fn bench_random_corpus_operation<const D: usize>(
    group: &mut BenchmarkGroup<'_, WallTime>,
    corpus: &[ValidatedExactInput<D>; RANDOM_INPUT_ARRAY_LEN],
    operation: ExactOperation,
) {
    group.bench_function(operation.name(), |bencher| {
        bencher.iter(|| {
            for input in corpus {
                run_exact_operation(operation, input);
            }
        });
    });
}

/// Populate a Criterion group with the five headline exact-arithmetic
/// benches on a single `(matrix, rhs)` pair: `det_sign_exact`,
/// `det_exact`, `solve_exact`, `solve_exact_f64_result`, and
/// `solve_exact_rounded_f64`.
///
/// Used by every adversarial-input group so each one measures the same
/// operations, making the resulting tables directly comparable.
fn bench_extreme_group<const D: usize>(
    group: &mut BenchmarkGroup<'_, WallTime>,
    input: &ValidatedExactInput<D>,
) {
    for &operation in CORPUS_AND_EXTREME_OPERATIONS {
        bench_exact_operation(group, operation, input);
    }
}

/// Add the direct-determinant baseline for a dimension that supports it.
fn bench_det_direct<const D: usize>(
    group: &mut BenchmarkGroup<'_, WallTime>,
    input: &ValidatedExactInput<D>,
) {
    let Some(_) = input
        .matrix()
        .det_direct()
        .or_abort("direct determinant setup")
    else {
        panic!("det_direct must support this benchmark dimension");
    };
    group.bench_function("det_direct", |bencher| {
        bencher.iter(|| {
            let det = black_box(input.matrix())
                .det_direct()
                .or_abort("direct f64 determinant");
            let Some(det) = det else {
                panic!("det_direct support changed after benchmark setup");
            };
            black_box(det);
        });
    });
}

/// Add the standalone certified determinant-bound baseline.
fn bench_det_errbound<const D: usize>(
    group: &mut BenchmarkGroup<'_, WallTime>,
    input: &ValidatedExactInput<D>,
) {
    let bound = input
        .matrix()
        .det_errbound()
        .or_abort("determinant error-bound setup")
        .or_abort("determinant error-bound setup");
    black_box(bound);
    group.bench_function("det_errbound", |bencher| {
        bencher.iter(|| {
            let bound = black_box(input.matrix())
                .det_errbound()
                .or_abort("f64 determinant error bound")
                .or_abort("f64 determinant error bound");
            black_box(bound);
        });
    });
}

/// Add the paired direct-determinant and certified-bound baseline.
#[cfg(not(la_stack_v0_4_3_api))]
fn bench_det_direct_with_errbound<const D: usize>(
    group: &mut BenchmarkGroup<'_, WallTime>,
    input: &ValidatedExactInput<D>,
) {
    let estimate = input
        .matrix()
        .det_direct_with_errbound()
        .or_abort("paired determinant-bound setup")
        .or_abort("paired determinant-bound setup");
    black_box((estimate.determinant(), estimate.absolute_error_bound()));
    group.bench_function("det_direct_with_errbound", |bencher| {
        bencher.iter(|| {
            let estimate = black_box(input.matrix())
                .det_direct_with_errbound()
                .or_abort("paired f64 determinant and error bound")
                .or_abort("paired f64 determinant and error bound");
            black_box((estimate.determinant(), estimate.absolute_error_bound()));
        });
    });
}

macro_rules! register_det_filter_benchmarks {
    ($group:expr, $input:expr, supported) => {{
        bench_det_direct(&mut $group, &$input);
        #[cfg(not(la_stack_v0_4_3_api))]
        bench_det_direct_with_errbound(&mut $group, &$input);
        bench_det_errbound(&mut $group, &$input);
    }};
    ($group:expr, $matrix:expr, unsupported) => {};
}

macro_rules! gen_exact_benches_for_dim {
    ($c:expr, $d:literal, $direct:ident) => {{
        let input = validate_exact_fixture(ExactInput {
            matrix: Matrix::<$d>::try_from_rows(make_matrix_rows::<$d>())
                .or_abort("benchmark matrix construction"),
            rhs: Vector::<$d>::try_new(make_vector_array::<$d>())
                .or_abort("benchmark RHS vector construction"),
        });
        validate_f64_determinant_benchmarks(&input);

        let mut group = ($c).benchmark_group(concat!("exact_d", stringify!($d)));

        // === f64 baselines ===
        group.bench_function("det", |bencher| {
            bencher.iter(|| {
                let det = black_box(input.matrix()).det().or_abort("f64 determinant");
                black_box(det);
            });
        });

        register_det_filter_benchmarks!(group, input, $direct);

        for &operation in GENERAL_OPERATIONS {
            bench_exact_operation(&mut group, operation, &input);
        }

        group.finish();
    }};
}

macro_rules! gen_random_corpus_benches_for_dim {
    ($c:expr, $d:literal) => {{
        let corpus = make_random_input_corpus::<$d>().map(validate_exact_fixture);

        let mut group = ($c).benchmark_group(concat!("exact_random_corpus_d", stringify!($d)));
        let input_count =
            u64::try_from(corpus.len()).or_abort("random corpus throughput conversion");
        group.throughput(Throughput::Elements(input_count));

        for &operation in CORPUS_AND_EXTREME_OPERATIONS {
            bench_random_corpus_operation(&mut group, &corpus, operation);
        }

        group.finish();
    }};
}

fn main() {
    let mut c = Criterion::default().configure_from_args();

    {
        gen_exact_benches_for_dim!(&mut c, 2, supported);
        gen_exact_benches_for_dim!(&mut c, 3, supported);
        gen_exact_benches_for_dim!(&mut c, 4, supported);
        gen_exact_benches_for_dim!(&mut c, 5, unsupported);
    }

    // === Fixed random-corpus groups ===
    //
    // Each measured iteration executes all 50 strictly diagonally-dominant
    // integer inputs in their fixed-seed order. Baseline and current revisions
    // therefore receive exactly the same workload.
    {
        gen_random_corpus_benches_for_dim!(&mut c, 2);
        gen_random_corpus_benches_for_dim!(&mut c, 3);
        gen_random_corpus_benches_for_dim!(&mut c, 4);
        gen_random_corpus_benches_for_dim!(&mut c, 5);
    }

    #[cfg(not(any(la_stack_pre_rational_input_api, la_stack_v0_4_3_api)))]
    {
        // === Already-exact rational-input comparisons ===
        //
        // These compare the production row-cleared integer Bareiss backend with
        // straightforward BigRational Gaussian elimination on dimensions needed
        // by downstream geometric predicates and runtime-selected basis systems.
        bench_rational_input::<2>(&mut c);
        bench_rational_input::<3>(&mut c);
        bench_rational_input::<4>(&mut c);
        bench_rational_input::<5>(&mut c);
        bench_rational_input::<6>(&mut c);
        bench_rational_input::<7>(&mut c);
        bench_rational_input::<8>(&mut c);
    }

    // === Adversarial / extreme-input groups ===
    //
    // Each group runs the same five exact-arithmetic benches
    // (`det_sign_exact`, `det_exact`, `solve_exact`, `solve_exact_f64_result`,
    // `solve_exact_rounded_f64`)
    // via `bench_extreme_group`, so the resulting tables are directly
    // comparable across input classes.

    // Near-singular 3×3: forces the direct BigInt fallback in det_sign_exact
    // and exercises an ill-conditioned exact solve.
    {
        let input = validate_exact_fixture(near_singular_3x3_input());
        let mut group = c.benchmark_group("exact_near_singular_3x3");
        bench_extreme_group(&mut group, &input);
        group.finish();
    }

    // Large-entry 3×3: diagonal entries near `f64::MAX / 2` stress
    // BigInt growth during Bareiss forward elimination.
    {
        let input = validate_exact_fixture(large_entries_3x3_input());
        let mut group = c.benchmark_group("exact_large_entries_3x3");
        bench_extreme_group(&mut group, &input);
        group.finish();
    }

    // Hilbert 4×4 and 5×5: classically ill-conditioned matrices whose
    // entries have varied binary mantissas and exponents, exercising the
    // f64 → BigInt scaling path.
    {
        let input = validate_exact_fixture(hilbert_input::<4>());
        let mut group = c.benchmark_group("exact_hilbert_4x4");
        bench_extreme_group(&mut group, &input);
        group.finish();
    }

    {
        let input = validate_exact_fixture(hilbert_input::<5>());
        let mut group = c.benchmark_group("exact_hilbert_5x5");
        bench_extreme_group(&mut group, &input);
        group.finish();
    }

    c.final_summary();
}
