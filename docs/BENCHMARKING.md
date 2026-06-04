# Benchmarking

This guide covers how to run, compare, and track performance for la-stack.

## Benchmark suites

la-stack has two Criterion benchmark suites:

- **`vs_linalg`** (`benches/vs_linalg.rs`) — compares la-stack against
  nalgebra and faer across D=2–64 for LU, solve, det, dot, norm, etc.
  Use this to answer "why choose la-stack over other crates?"
  The suite also includes SPD factorization rows for la-stack LDLT, faer
  LDLT, and nalgebra Cholesky. The nalgebra rows are labelled Cholesky
  because nalgebra does not expose a dense LDLT factorization in the
  dependency version used here.

- **`exact`** (`benches/exact.rs`) — measures exact-arithmetic methods
  (`det_exact`, `solve_exact`, `det_sign_exact`, etc.) alongside f64
  baselines (`det`, `det_direct`) across D=2–5. Use this to understand
  the cost of exact arithmetic and track optimization progress.
  In addition to the fixed per-dimension groups (`exact_d{2..5}`), the
  suite includes random percentile and adversarial-input groups designed
  to capture variance and stress specific corners of the pipeline:

  - `exact_random_percentile_d{2..5}` — fixed-seed corpora of 50
    strictly diagonally-dominant random matrices per dimension. Each
    operation is pre-timed across the corpus to select representative
    p50/p95/p99 inputs, then Criterion measures those inputs normally.
  - `exact_near_singular_3x3` — a 2^-50 perturbation of a singular base
    matrix; forces the Bareiss fallback in `det_sign_exact` and
    exercises the largest intermediate `BigInt` values in `solve_exact`.
  - `exact_large_entries_3x3` — diagonal entries near `f64::MAX / 2`
    stress `BigInt` growth during Bareiss forward elimination.
  - `exact_hilbert_4x4` / `exact_hilbert_5x5` — classically
    ill-conditioned matrices whose non-terminating-in-binary entries
    stress the `f64_decompose → BigInt` scaling path.

  Each random percentile and adversarial group runs the same four
  exact-arithmetic benches (`det_sign_exact`, `det_exact`, `solve_exact`,
  `solve_exact_f64`) so the resulting tables are directly comparable
  across input classes.

## `vs_linalg` methodology

`vs_linalg` is a per-kernel comparison, not a single aggregate score. Each
reported row compares one operation for one dimension `D`, using Criterion's
selected statistic from `target/criterion/d{D}/{benchmark}/{sample}/estimates.json`.
The default report and README table use Criterion's `median.point_estimate`
in nanoseconds. Lower is better.

All three crates receive equivalent deterministic inputs for a given
dimension:

- matrix entries come from the same strictly diagonally-dominant generator
  (`matrix_entry::<D>`)
- right-hand sides and vector inputs come from the same deterministic vector
  generator
- each benchmark uses `black_box` around inputs and outputs to keep the
  measured operation visible to the optimizer

The integration smoke test `tests/vs_linalg_inputs.rs` reuses the benchmark
input helpers and verifies that la-stack, nalgebra, and faer agree on the
determinant, solve, dot, and infinity-norm results for D=2..=5. Run it with
`cargo test --features bench --test vs_linalg_inputs` when changing benchmark
input construction, adding comparable kernels, or updating the `faer` or
`nalgebra` benchmark dependencies.

The main comparable metrics are:

- `det_via_lu` — factor the matrix and compute determinant from the LU factor
- `lu` — factorization only
- `lu_solve` — factor the matrix and solve one right-hand side
- `solve_from_lu` — solve one right-hand side using a precomputed LU factor
- `det_from_lu` — compute determinant using a precomputed LU factor
- `dot` — vector dot product
- `norm2_sq` — squared Euclidean vector norm
- `inf_norm` — matrix infinity norm, implemented as maximum absolute row sum

Additional SPD metrics compare la-stack LDLT against faer LDLT and nalgebra
Cholesky. These rows are labelled by algorithm (`ldlt` or `cholesky`) because
nalgebra does not expose a dense LDLT factorization in the dependency version
used here. They should be read as SPD factorization/solve/determinant
comparisons, not as identical algorithm comparisons across all three crates.

Release-signal reports compare latest la-stack measurements against a saved
la-stack baseline, and show saved nalgebra/faer baseline timings as context
where a matching peer benchmark exists. That keeps iteration cheap while still
making the release signal auditable. The full `vs_linalg` run remains the
source of README plots and crate-to-crate comparison tables.

## Quick reference

```bash
# Run vs_linalg benchmarks
just bench-vs-linalg

# Run only la-stack rows from vs_linalg
just bench-vs-linalg-la-stack

# Run exact-arithmetic benchmarks
just bench-exact

# Run the cheaper latest measurements used for latest-vs-last reports
just bench-latest

# Save a full baseline named "last"
just bench-save-last

# Compare latest measurements against the saved "last" baseline
just bench-compare

# Run latest measurements and compare against "last"
just bench-latest-vs-last
```

## Comparing performance across releases

Criterion baselines are saved into `target/criterion/` and persist across
`git checkout` but **not** across `cargo clean`. Published releases also attach
a compressed Criterion baseline to the GitHub Release so historical release
baselines can be restored later.

### Latest vs last

The default workflow is optimized for the common maintenance question:
"how does latest la-stack compare to the last release?"

At release time, save a full baseline:

```bash
just bench-save-last
```

During development, run the cheaper latest path:

```bash
just bench-latest-vs-last
```

`bench-latest` runs exact arithmetic plus only the la-stack rows from
`vs_linalg`. The comparison report still shows the last-release nalgebra
and faer timings for matching rows, so you can see whether a la-stack
change improves or weakens the release signal without rerunning third-party
benchmarks on every iteration.

### Workflow

```bash
# 1. Check out the old release and save its full baseline
git checkout v0.2.0
just bench-save-baseline v0.2.0

# 2. Switch to current code and run latest la-stack measurements
git checkout main   # or your feature branch
just bench-latest   # populates target/criterion/*/new/

# 3. Generate a local comparison report
just bench-compare v0.2.0
```

You can save multiple baselines and compare against any of them.

If the release baseline is already present in `target/criterion/`, skip the
checkout step and compare directly. For example, to compare current code against
the saved `v0.4.2` release baseline:

```bash
just bench-latest          # gather latest la-stack measurements
just bench-compare v0.4.2  # compare latest measurements against v0.4.2
```

If the release baseline is not present locally, download and restore the release
asset first:

```bash
gh release download v0.4.2 --pattern "la-stack-v0.4.2-criterion-baseline.tar.gz"  # fetch archived release baseline
mkdir -p target                                                                    # ensure Criterion parent directory exists
tar -C target -xzf la-stack-v0.4.2-criterion-baseline.tar.gz                       # restore target/criterion baseline data
just bench-latest                                                                  # gather latest la-stack measurements
just bench-compare v0.4.2                                                          # compare latest measurements against v0.4.2
```

### Output

`just bench-compare` writes `target/bench-reports/performance.md` by
default. The file contains machine-specific timings and is intentionally
local. The report includes per-dimension tables showing median times,
percent change, speedup, and last-release nalgebra/faer context where a
matching `vs_linalg` peer exists.

To generate a current snapshot without a saved baseline:

```bash
uv run bench-compare --snapshot
```

## vs\_linalg plotting

The `criterion_dim_plot.py` script generates CSV/SVG plots and updates the
README benchmark table from vs\_linalg results:

```bash
# Run benchmarks + update README table and SVG plot
just bench-vs-linalg
just plot-vs-linalg-readme
```

See `scripts/criterion_dim_plot.py --help` for options.

## Release workflow

At release time, save a local baseline so future work can compare against it:

```bash
just bench-save-baseline $TAG
just bench-save-last
```

When the GitHub Release is published, `.github/workflows/release-benchmarks.yml`
saves a full release baseline and attaches
`la-stack-$TAG-criterion-baseline.tar.gz` to the release as the durable archive.
See `docs/RELEASING.md` step 5 for where this fits in the release process.
