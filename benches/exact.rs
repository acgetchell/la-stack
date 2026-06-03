//! Benchmarks for exact arithmetic operations.
//!
//! These benchmarks measure the performance of the `exact` feature's
//! arbitrary-precision methods.  They are organised into three classes:
//!
//! 1. **General-case benches** (`exact_d{2..5}`) — a single
//!    well-conditioned diagonally-dominant matrix per dimension.  These
//!    measure typical-case performance and track regressions against a
//!    reproducible input.
//! 2. **Adversarial / extreme-input benches** — matrices chosen to
//!    stress specific corners of the exact-arithmetic pipeline:
//!    near-singularity (forces the Bareiss fallback), large f64 entries
//!    (stresses intermediate `BigInt` growth), and Hilbert-style
//!    ill-conditioning (wide range of `(mantissa, exponent)` pairs in
//!    the `f64_decompose → BigInt` path).  These measure tail behaviour
//!    that fixed well-conditioned inputs miss and provide stronger
//!    empirical evidence for `docs/PERFORMANCE.md`.
//! 3. **Random percentile benches** (`exact_random_percentile_d{2..5}`) —
//!    a fixed-seed corpus of diagonally-dominant random matrices per
//!    dimension.  Each operation is pre-timed across the corpus to select
//!    p50/p95/p99 cumulative input subsets, then measured with Criterion.

use std::array;
use std::cell::Cell;
use std::fmt::{self, Display};
use std::hint::black_box;
use std::num::NonZeroUsize;
use std::time::Instant;

use criterion::{BatchSize, BenchmarkGroup, Criterion, measurement::WallTime};
use pastey::paste;

use la_stack::{Matrix, Vector};

const RANDOM_INPUTS_PER_DIM: SampleCount = SampleCount::new_unchecked(50);
const RANDOM_INPUT_ARRAY_LEN: usize = RANDOM_INPUTS_PER_DIM.get();
const RANDOM_TIMING_PASSES: SampleCount = SampleCount::new_unchecked(5);
const RANDOM_SEED: [u8; 32] = [0; 32];
const RANDOM_PERCENTILES: [RandomPercentile; 3] = [
    RandomPercentile::P50,
    RandomPercentile::P95,
    RandomPercentile::P99,
];

/// Return a successful benchmark operation result or panic with the named operation.
fn require_ok<T, E: Display>(result: Result<T, E>, operation: &str) -> T {
    match result {
        Ok(value) => value,
        Err(err) => panic!("{operation} failed: {err}"),
    }
}

/// Configuration errors for exact-arithmetic benchmark input generation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ExactBenchConfigError {
    EmptyCorpus,
    UnorderedRange { min: i16, max: i16 },
}

impl Display for ExactBenchConfigError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::EmptyCorpus => f.write_str("random input corpus must be nonempty"),
            Self::UnorderedRange { min, max } => {
                write!(f, "random integer range must be ordered: {min}..={max}")
            }
        }
    }
}

/// Non-zero sample count used when selecting percentile benchmark inputs.
#[derive(Clone, Copy)]
struct SampleCount {
    len: NonZeroUsize,
}

impl SampleCount {
    /// Construct a sample count for compile-time constants with visible nonzero values.
    const fn new_unchecked(len: usize) -> Self {
        match NonZeroUsize::new(len) {
            Some(len) => Self { len },
            None => panic!("random input corpus must be nonempty"),
        }
    }

    /// Validate a runtime sample count before percentile calculations use it.
    const fn new(len: usize) -> Result<Self, ExactBenchConfigError> {
        if let Some(len) = NonZeroUsize::new(len) {
            Ok(Self { len })
        } else {
            Err(ExactBenchConfigError::EmptyCorpus)
        }
    }

    /// Return the proven nonzero sample count as a raw `usize`.
    const fn get(self) -> usize {
        self.len.get()
    }
}

/// Inclusive integer range used by the fixed-seed exact benchmark generator.
#[derive(Clone, Copy)]
struct I16Range {
    min: i16,
    width: u64,
}

impl I16Range {
    /// Validate an inclusive `i16` range and cache its sampling width.
    fn new(min: i16, max: i16) -> Result<Self, ExactBenchConfigError> {
        if min > max {
            return Err(ExactBenchConfigError::UnorderedRange { min, max });
        }

        let width = i32::from(max) - i32::from(min) + 1;
        Ok(Self {
            min,
            width: u64::try_from(width)
                .map_err(|_| ExactBenchConfigError::UnorderedRange { min, max })?,
        })
    }
}

/// Percentiles selected from a pre-timed random-input corpus.
#[derive(Clone, Copy)]
enum RandomPercentile {
    P50,
    P95,
    P99,
}

impl RandomPercentile {
    /// Return the percentile value as an integer percentage.
    const fn value(self) -> usize {
        match self {
            Self::P50 => 50,
            Self::P95 => 95,
            Self::P99 => 99,
        }
    }

    /// Return the benchmark-name suffix for this percentile.
    const fn name(self) -> &'static str {
        match self {
            Self::P50 => "p50",
            Self::P95 => "p95",
            Self::P99 => "p99",
        }
    }
}

/// Return a deterministic, strictly diagonally-dominant benchmark matrix entry.
#[inline]
#[allow(clippy::cast_precision_loss)]
const fn matrix_entry<const D: usize>(r: usize, c: usize) -> f64 {
    if r == c {
        (r as f64).mul_add(1.0e-3, (D as f64) + 1.0)
    } else {
        0.1 / ((r + c + 1) as f64)
    }
}

/// Build the deterministic baseline matrix rows for dimension `D`.
#[inline]
const fn make_matrix_rows<const D: usize>() -> [[f64; D]; D] {
    let mut rows = [[0.0; D]; D];
    let mut r = 0;
    while r < D {
        let mut c = 0;
        while c < D {
            rows[r][c] = matrix_entry::<D>(r, c);
            c += 1;
        }
        r += 1;
    }
    rows
}

/// Build the deterministic baseline right-hand-side vector for dimension `D`.
#[inline]
#[allow(clippy::cast_precision_loss)]
fn make_vector_array<const D: usize>() -> [f64; D] {
    let mut data = [0.0; D];
    let mut i = 0;
    while i < D {
        data[i] = (i as f64) + 1.0;
        i += 1;
    }
    data
}

/// Matrix/RHS pair used by random percentile exact-arithmetic benchmarks.
#[derive(Clone, Copy)]
struct ExactRandomInput<const D: usize> {
    matrix: Matrix<D>,
    rhs: Vector<D>,
}

/// Exact operation timed when selecting representative random inputs.
#[derive(Clone, Copy)]
enum ExactRandomOperation {
    DetSignExact,
    DetExact,
    SolveExact,
    SolveExactF64,
}

impl ExactRandomOperation {
    /// Return the benchmark-name stem for this exact operation.
    const fn name(self) -> &'static str {
        match self {
            Self::DetSignExact => "det_sign_exact",
            Self::DetExact => "det_exact",
            Self::SolveExact => "solve_exact",
            Self::SolveExactF64 => "solve_exact_f64",
        }
    }
}

/// Deterministic `SplitMix64` generator for reproducible benchmark corpora.
struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    /// Initialize the generator with a fixed state.
    const fn new(state: u64) -> Self {
        Self { state }
    }

    /// Advance the generator and return the next 64 random bits.
    const fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    #[allow(clippy::cast_possible_truncation)]
    /// Draw a random `i16` inside a validated inclusive range.
    fn next_i16(&mut self, range: I16Range) -> i16 {
        let offset = (self.next_u64() % range.width) as i32;
        let value = i32::from(range.min) + offset;
        value as i16
    }
}

/// Derive a stable per-dimension seed from the global random benchmark seed.
#[allow(clippy::cast_possible_truncation)]
fn random_seed_for_dim<const D: usize>() -> u64 {
    let mut seed =
        0xC0DE_CAFE_D15C_A11Au64 ^ require_ok(u64::try_from(D), "dimension seed conversion");
    for (i, byte) in RANDOM_SEED.iter().copied().enumerate() {
        let shift = require_ok(u32::try_from((i % 8) * 8), "seed shift conversion");
        seed ^= u64::from(byte) << shift;
        seed = seed.rotate_left(7) ^ require_ok(u64::try_from(i), "seed index conversion");
    }
    seed
}

/// Build a fixed random corpus of finite, strictly diagonally-dominant inputs.
fn make_random_input_corpus<const D: usize>() -> [ExactRandomInput<D>; RANDOM_INPUT_ARRAY_LEN] {
    let mut rng = SplitMix64::new(random_seed_for_dim::<D>());
    let entry_range = require_ok(I16Range::new(-10, 10), "random integer range");
    array::from_fn(|_| {
        let mut rows = [[0.0; D]; D];
        let mut diag = [0_i16; D];

        for (r, row) in rows.iter_mut().enumerate() {
            for (c, entry) in row.iter_mut().enumerate() {
                if r == c {
                    diag[r] = rng.next_i16(entry_range);
                } else {
                    *entry = f64::from(rng.next_i16(entry_range));
                }
            }
        }

        let shift =
            f64::from(require_ok(u8::try_from(D), "dimension shift conversion")).mul_add(10.0, 1.0);
        for (i, row) in rows.iter_mut().enumerate() {
            row[i] = if diag[i] >= 0 {
                f64::from(diag[i]) + shift
            } else {
                f64::from(diag[i]) - shift
            };
        }

        let rhs = array::from_fn(|_| f64::from(rng.next_i16(entry_range)));

        ExactRandomInput {
            matrix: require_ok(
                Matrix::<D>::try_from_rows(rows),
                "random matrix construction",
            ),
            rhs: require_ok(Vector::<D>::try_new(rhs), "random RHS vector construction"),
        }
    })
}

/// Execute one exact operation on a random benchmark input.
fn run_random_operation<const D: usize>(
    operation: ExactRandomOperation,
    input: ExactRandomInput<D>,
) {
    match operation {
        ExactRandomOperation::DetSignExact => {
            let sign = require_ok(
                black_box(input.matrix).det_sign_exact(),
                "exact determinant sign",
            );
            black_box(sign);
        }
        ExactRandomOperation::DetExact => {
            let det = require_ok(black_box(input.matrix).det_exact(), "exact determinant");
            black_box(det);
        }
        ExactRandomOperation::SolveExact => {
            let x = require_ok(
                black_box(input.matrix).solve_exact(black_box(input.rhs)),
                "exact linear solve",
            );
            let _ = black_box(x);
        }
        ExactRandomOperation::SolveExactF64 => {
            let x = require_ok(
                black_box(input.matrix).solve_exact_f64(black_box(input.rhs)),
                "exact linear solve converted to f64",
            );
            let _ = black_box(x);
        }
    }
}

/// Time one exact operation on one random input in nanoseconds.
fn time_random_operation<const D: usize>(
    operation: ExactRandomOperation,
    input: ExactRandomInput<D>,
) -> u128 {
    let start = Instant::now();
    run_random_operation(operation, input);
    start.elapsed().as_nanos()
}

/// Time one exact operation repeatedly on one random input.
fn time_random_operation_repeated<const D: usize>(
    operation: ExactRandomOperation,
    input: ExactRandomInput<D>,
) -> u128 {
    let mut elapsed = 0;
    for _ in 0..RANDOM_TIMING_PASSES.get() {
        elapsed += time_random_operation(operation, input);
    }
    elapsed
}

/// Convert a percentile request into an index in a sorted timing corpus.
const fn percentile_index(count: SampleCount, percentile: RandomPercentile) -> usize {
    ((count.get() - 1) * percentile.value() + 50) / 100
}

/// Select cumulative corpus index sets by pre-timing every input for one operation.
fn percentile_input_indices<const D: usize>(
    corpus: &[ExactRandomInput<D>; RANDOM_INPUT_ARRAY_LEN],
    operation: ExactRandomOperation,
) -> [Vec<usize>; RANDOM_PERCENTILES.len()] {
    let input_count = require_ok(SampleCount::new(corpus.len()), "random input corpus size");
    let mut timings = [(0_u128, 0_usize); RANDOM_INPUT_ARRAY_LEN];
    for (i, input) in corpus.iter().copied().enumerate() {
        timings[i] = (time_random_operation_repeated(operation, input), i);
    }
    timings.sort_unstable();

    RANDOM_PERCENTILES.map(|percentile| {
        let timing_idx = percentile_index(input_count, percentile);
        let threshold = timings[timing_idx].0;
        let mut indices = Vec::new();
        for &(elapsed, input_idx) in &timings {
            if elapsed <= threshold {
                indices.push(input_idx);
            }
        }
        indices
    })
}

/// Add p50/p95/p99 Criterion benches over percentile input sets.
fn bench_random_percentile_operation<const D: usize>(
    group: &mut BenchmarkGroup<'_, WallTime>,
    corpus: &[ExactRandomInput<D>; RANDOM_INPUT_ARRAY_LEN],
    operation: ExactRandomOperation,
) {
    let index_sets = percentile_input_indices(corpus, operation);

    for (percentile, input_indices) in RANDOM_PERCENTILES.into_iter().zip(index_sets) {
        let input_count = require_ok(
            SampleCount::new(input_indices.len()),
            "percentile input set size",
        );
        let cursor = Cell::new(0);
        group.bench_function(
            format!("{}_{}", operation.name(), percentile.name()),
            move |bencher| {
                bencher.iter_batched(
                    || {
                        let cursor_pos = cursor.get();
                        cursor.set((cursor_pos + 1) % input_count.get());
                        corpus[input_indices[cursor_pos]]
                    },
                    |sample| run_random_operation(operation, sample),
                    BatchSize::SmallInput,
                );
            },
        );
    }
}

/// Near-singular matrix: base singular matrix + tiny perturbation.
///
/// The base `[[1,2,3],[4,5,6],[7,8,9]]` is exactly singular; adding
/// `2^-50` to entry (0,0) makes `det = -3 × 2^-50 ≠ 0`.  The f64 filter
/// in `det_sign_exact` cannot resolve this sign, so Bareiss is forced;
/// `solve_exact` is the primary use case for near-degenerate inputs
/// (exact circumcenter etc.) and exercises the largest intermediate
/// `BigInt` values in the hybrid solve.
#[inline]
fn near_singular_3x3() -> Matrix<3> {
    let perturbation = f64::from_bits(0x3CD0_0000_0000_0000); // 2^-50
    require_ok(
        Matrix::<3>::try_from_rows([
            [1.0 + perturbation, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
        ]),
        "near-singular matrix construction",
    )
}

/// Large-entry 3×3: strictly diagonally-dominant matrix with diagonal
/// entries near `f64::MAX / 2` and ones elsewhere.
///
/// Each big entry decomposes into a 53-bit mantissa with exponent `~970`;
/// the unit off-diagonals have exponent `0`, so the shared `e_min = 0`
/// shift in `component_to_bigint` produces `BigInt`s of `~1023` bits for
/// the diagonal and small integers elsewhere.  Bareiss fraction-free
/// updates then multiply these together, stressing the big-integer
/// multiply and allocator along the full `O(D³)` elimination phase.  The
/// matrix is non-singular (det ≈ `big³`) so both `det_*` and `solve_*`
/// paths complete.
#[inline]
fn large_entries_3x3() -> Matrix<3> {
    let big = f64::MAX / 2.0;
    require_ok(
        Matrix::<3>::try_from_rows([[big, 1.0, 1.0], [1.0, big, 1.0], [1.0, 1.0, big]]),
        "large-entry matrix construction",
    )
}

/// Hilbert matrix `H[i][j] = 1 / (i + j + 1)`.
///
/// Most entries (`1/3`, `1/5`, `1/6`, `1/7`, …) are non-terminating in
/// binary, so every cell has a distinct 53-bit mantissa and a small
/// negative exponent.  `f64_decompose` therefore produces a wide mix of
/// `(mantissa, exponent)` pairs with no shared power-of-two factors,
/// and the scaling shift to the common `e_min` yields `BigInt` values
/// of varied bit-lengths — a different kind of adversarial input from
/// the large-entries case.  Hilbert matrices are also classically
/// ill-conditioned (condition number grows exponentially with D), so
/// they are a realistic stand-in for the near-degenerate geometric
/// predicate inputs that motivate exact arithmetic.
#[inline]
#[allow(clippy::cast_precision_loss)]
fn hilbert<const D: usize>() -> Matrix<D> {
    let mut rows = [[0.0; D]; D];
    let mut r = 0;
    while r < D {
        let mut c = 0;
        while c < D {
            rows[r][c] = 1.0 / ((r + c + 1) as f64);
            c += 1;
        }
        r += 1;
    }
    require_ok(
        Matrix::<D>::try_from_rows(rows),
        "Hilbert matrix construction",
    )
}

/// Populate a Criterion group with the four headline exact-arithmetic
/// benches on a single `(matrix, rhs)` pair: `det_sign_exact`,
/// `det_exact`, `solve_exact`, `solve_exact_f64`.
///
/// Used by every adversarial-input group so each one measures the same
/// operations, making the resulting tables directly comparable.
fn bench_extreme_group<const D: usize>(
    group: &mut BenchmarkGroup<'_, WallTime>,
    m: Matrix<D>,
    rhs: Vector<D>,
) {
    group.bench_function("det_sign_exact", |bencher| {
        bencher.iter(|| {
            let sign = require_ok(black_box(m).det_sign_exact(), "exact determinant sign");
            black_box(sign);
        });
    });

    group.bench_function("det_exact", |bencher| {
        bencher.iter(|| {
            let det = require_ok(black_box(m).det_exact(), "exact determinant");
            black_box(det);
        });
    });

    group.bench_function("solve_exact", |bencher| {
        bencher.iter(|| {
            let x = require_ok(
                black_box(m).solve_exact(black_box(rhs)),
                "exact linear solve",
            );
            let _ = black_box(x);
        });
    });

    group.bench_function("solve_exact_f64", |bencher| {
        bencher.iter(|| {
            let x = require_ok(
                black_box(m).solve_exact_f64(black_box(rhs)),
                "exact linear solve converted to f64",
            );
            let _ = black_box(x);
        });
    });
}

macro_rules! gen_exact_benches_for_dim {
    ($c:expr, $d:literal) => {
        paste! {{
            let a = require_ok(
                Matrix::<$d>::try_from_rows(make_matrix_rows::<$d>()),
                "benchmark matrix construction",
            );
            let rhs = require_ok(
                Vector::<$d>::try_new(make_vector_array::<$d>()),
                "benchmark RHS vector construction",
            );

            let mut [<group_d $d>] = ($c).benchmark_group(concat!("exact_d", stringify!($d)));

            // === f64 baselines ===
            [<group_d $d>].bench_function("det", |bencher| {
                bencher.iter(|| {
                    let det = require_ok(black_box(a).det(), "f64 determinant");
                    black_box(det);
                });
            });

            [<group_d $d>].bench_function("det_direct", |bencher| {
                bencher.iter(|| {
                    let det = black_box(a).det_direct();
                    black_box(det);
                });
            });

            // === det_exact (BigRational result) ===
            [<group_d $d>].bench_function("det_exact", |bencher| {
                bencher.iter(|| {
                    let det = require_ok(black_box(a).det_exact(), "exact determinant");
                    black_box(det);
                });
            });

            // === det_exact_f64 (exact → f64) ===
            [<group_d $d>].bench_function("det_exact_f64", |bencher| {
                bencher.iter(|| {
                    let det = require_ok(
                        black_box(a).det_exact_f64(),
                        "exact determinant converted to f64",
                    );
                    black_box(det);
                });
            });

            // === det_sign_exact (adaptive: fast filter + exact fallback) ===
            [<group_d $d>].bench_function("det_sign_exact", |bencher| {
                bencher.iter(|| {
                    let sign = require_ok(black_box(a).det_sign_exact(), "exact determinant sign");
                    black_box(sign);
                });
            });

            // === solve_exact (BigRational result) ===
            [<group_d $d>].bench_function("solve_exact", |bencher| {
                bencher.iter(|| {
                    let x = require_ok(
                        black_box(a).solve_exact(black_box(rhs)),
                        "exact linear solve",
                    );
                    black_box(x);
                });
            });

            // === solve_exact_f64 (exact → f64) ===
            [<group_d $d>].bench_function("solve_exact_f64", |bencher| {
                bencher.iter(|| {
                    let x = require_ok(
                        black_box(a).solve_exact_f64(black_box(rhs)),
                        "exact linear solve converted to f64",
                    );
                    black_box(x);
                });
            });

            [<group_d $d>].finish();
        }};
    };
}

macro_rules! gen_random_percentile_benches_for_dim {
    ($c:expr, $d:literal) => {
        paste! {{
            let corpus = make_random_input_corpus::<$d>();
            let mut [<group_random_percentile_d $d>] =
                ($c).benchmark_group(concat!("exact_random_percentile_d", stringify!($d)));

            bench_random_percentile_operation(
                &mut [<group_random_percentile_d $d>],
                &corpus,
                ExactRandomOperation::DetSignExact,
            );
            bench_random_percentile_operation(
                &mut [<group_random_percentile_d $d>],
                &corpus,
                ExactRandomOperation::DetExact,
            );
            bench_random_percentile_operation(
                &mut [<group_random_percentile_d $d>],
                &corpus,
                ExactRandomOperation::SolveExact,
            );
            bench_random_percentile_operation(
                &mut [<group_random_percentile_d $d>],
                &corpus,
                ExactRandomOperation::SolveExactF64,
            );

            [<group_random_percentile_d $d>].finish();
        }};
    };
}

fn main() {
    let mut c = Criterion::default().configure_from_args();

    #[allow(unused_must_use)]
    {
        gen_exact_benches_for_dim!(&mut c, 2);
        gen_exact_benches_for_dim!(&mut c, 3);
        gen_exact_benches_for_dim!(&mut c, 4);
        gen_exact_benches_for_dim!(&mut c, 5);
    }

    // === Random percentile groups ===
    //
    // Each dimension uses a fixed-seed corpus of strictly
    // diagonally-dominant integer matrices.  For each operation, the corpus
    // is pre-timed repeatedly to select cumulative p50/p95/p99 input sets,
    // then Criterion cycles through each set with normal sampling.
    #[allow(unused_must_use)]
    {
        gen_random_percentile_benches_for_dim!(&mut c, 2);
        gen_random_percentile_benches_for_dim!(&mut c, 3);
        gen_random_percentile_benches_for_dim!(&mut c, 4);
        gen_random_percentile_benches_for_dim!(&mut c, 5);
    }

    // === Adversarial / extreme-input groups ===
    //
    // Each group runs the same four exact-arithmetic benches
    // (`det_sign_exact`, `det_exact`, `solve_exact`, `solve_exact_f64`)
    // via `bench_extreme_group`, so the resulting tables are directly
    // comparable across input classes.

    // Near-singular 3×3: forces Bareiss fallback in det_sign_exact and
    // exercises the largest intermediate BigInt values in solve_exact
    // (the primary motivating use case for exact solve).
    {
        let mut group = c.benchmark_group("exact_near_singular_3x3");
        bench_extreme_group(
            &mut group,
            near_singular_3x3(),
            require_ok(
                Vector::<3>::try_new([1.0, 2.0, 3.0]),
                "near-singular RHS vector construction",
            ),
        );
        group.finish();
    }

    // Large-entry 3×3: diagonal entries near `f64::MAX / 2` stress
    // BigInt growth during Bareiss forward elimination.
    {
        let mut group = c.benchmark_group("exact_large_entries_3x3");
        bench_extreme_group(
            &mut group,
            large_entries_3x3(),
            require_ok(
                Vector::<3>::try_new([1.0, 1.0, 1.0]),
                "large-entry RHS vector construction",
            ),
        );
        group.finish();
    }

    // Hilbert 4×4 and 5×5: classically ill-conditioned matrices whose
    // entries span many orders of magnitude in `(mantissa, exponent)`
    // space, exercising the f64 → BigInt scaling path.
    {
        let mut group = c.benchmark_group("exact_hilbert_4x4");
        bench_extreme_group(
            &mut group,
            hilbert::<4>(),
            require_ok(
                Vector::<4>::try_new([1.0; 4]),
                "Hilbert RHS vector construction",
            ),
        );
        group.finish();
    }

    {
        let mut group = c.benchmark_group("exact_hilbert_5x5");
        bench_extreme_group(
            &mut group,
            hilbert::<5>(),
            require_ok(
                Vector::<5>::try_new([1.0; 5]),
                "Hilbert RHS vector construction",
            ),
        );
        group.finish();
    }

    c.final_summary();
}
