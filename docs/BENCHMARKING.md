# Benchmarking

This guide explains how to run, compare, and publish performance results for
`la-stack`. Start with the workflow table below; the later sections explain what
the commands measure and where their outputs go.

## Contents

- [Start Here](#start-here)
- [Benchmark Suites](#benchmark-suites)
- [Common Workflows](#common-workflows)
  - [Compare Current Code With The Latest Release](#compare-current-code-with-the-latest-release)
  - [Compare Current Code With A Specific Release](#compare-current-code-with-a-specific-release)
  - [Iterate Against A Local Saved Baseline](#iterate-against-a-local-saved-baseline)
  - [Update The README nalgebra/faer Table](#update-the-readme-nalgebrafaer-table)
  - [Create The Release Performance Report](#create-the-release-performance-report)
  - [Compare Published Release Artifacts](#compare-published-release-artifacts)
- [Output Locations](#output-locations)
- [`vs_linalg` Methodology](#vs_linalg-methodology)
- [Exact-Arithmetic Notes](#exact-arithmetic-notes)
- [Release Notes](#release-notes)

## Start Here

| Goal | Use | Output | Notes |
|------|-----|--------|-------|
| Clean local audit against the latest published release | `just performance-local` | `target/bench-reports/performance.md` | Self-contained; creates temporary worktrees and regenerates the release baseline locally. |
| Non-exact release-signal check against a specific release | `just performance-local-vs-linalg v0.4.3 v0.4.2` | `target/bench-reports/performance.md` | Narrower than `performance-local`; useful for LU/LDLT/dot/norm work. |
| Fast repeated comparisons while tuning one kernel | `just bench-save-baseline <name> <suite>` then `just bench-compare <name> <suite> all-benches` | `target/bench-reports/performance.md` | Uses local `target/criterion/`; fastest loop after the baseline exists. |
| Full current la-stack vs nalgebra/faer comparison | `just bench-vs-linalg` | `target/criterion/` | Measures current la-stack, nalgebra, and faer rows. |
| README benchmark table and SVG plot | `just plot-vs-linalg-readme` after `just bench-vs-linalg` | `README.md`, `docs/assets/bench/` | Uses current `target/criterion` data. |
| Release PR performance artifact | `just performance-release v0.4.3 v0.4.2` | `docs/PERFORMANCE.md`, `docs/archive/performance/` | Mutates committed docs. Run during release preparation. |
| Compare already-published release assets | `just performance-github-assets v0.4.3 v0.4.2` | `target/bench-reports/github-assets-performance.md` | Uses GitHub Release baseline assets instead of local cargo runs. |

Rule of thumb:

- Use `performance-local*` for clean, self-contained answers.
- Use `bench-save-*` plus `bench-compare` for tight local optimization loops.
- Use `bench-vs-linalg` plus plotting when updating README crate-to-crate
  comparisons.
- Use `performance-release` only when preparing committed release artifacts.

## Benchmark Suites

`la-stack` has two Criterion benchmark suites.

**`vs_linalg`** (`benches/vs_linalg.rs`) compares `la-stack` against
`nalgebra` and `faer` across D=2-64 for LU, solve, determinant, dot, norm, and
SPD factorization operations. Use this suite to answer "why choose la-stack over
other crates?"

The SPD rows compare la-stack LDLT, faer LDLT, and nalgebra Cholesky. They are
labelled by algorithm because nalgebra does not expose a dense LDLT
factorization in the dependency version used here.

**`exact`** (`benches/exact.rs`) measures exact-arithmetic methods
(`det_exact`, `solve_exact`, `det_sign_exact`, strict `*_result` conversions,
and lossy `*_rounded_f64` conversions) alongside f64 baselines (`det`,
`det_direct`) across D=2-5. Use this suite to understand exact-arithmetic cost
and track optimization progress.

## Common Workflows

### Compare Current Code With The Latest Release

Use this when you want a clean local answer to "how does this checkout compare
with the latest published release?"

```bash
just performance-local
```

This creates isolated temporary worktrees, generates the latest published
release baseline locally, benchmarks the current tree on the same machine, and
writes `target/bench-reports/performance.md`.

This command does not depend on existing local `target/criterion/` baselines.
It is slower than reusing a saved baseline, but less sensitive to stale local
benchmark state.

### Compare Current Code With A Specific Release

For a narrower non-exact check against a known release pair, run:

```bash
just performance-local-vs-linalg v0.4.3 v0.4.2
```

This generates a local `v0.4.2` `vs_linalg` baseline, measures the current
la-stack `vs_linalg` rows, and renders a `vs_linalg` report. The report includes
saved baseline nalgebra/faer timings as context where matching peer rows exist,
without rerunning current peer crates.

### Iterate Against A Local Saved Baseline

Use local saved baselines when tuning one kernel and comparing several edits
against the same starting point. These baselines are local scratch data, not
release artifacts.

For example, before optimizing `Matrix::inf_norm`, save a named baseline:

```bash
just bench-save-baseline inf-norm-before vs_linalg
```

Then make a change, rerun only the current measurements you care about, and
compare:

```bash
just bench-vs-linalg-la-stack
just bench-compare inf-norm-before vs_linalg all-benches
```

The `just bench-compare` recipe uses positional arguments:
`just bench-compare <baseline> <suite> <scope>`. The underlying
`uv run bench-compare` CLI accepts the explicit `--suite` and `--scope` flags.

`just bench-save-baseline <name>` writes Criterion samples under
`target/criterion/`. `just bench-save-last` saves the conventional local
baseline named `last`, which enables shortcuts such as:

```bash
just bench-latest-vs-last
just bench-vs-linalg-latest-vs
just bench-compare
```

Saved baselines persist across `git checkout` but not across `cargo clean`, and
they are not pushed to GitHub.

### Update The README nalgebra/faer Table

The README benchmark table and SVG plot are crate-to-crate comparisons from the
current checkout:

```bash
just bench-vs-linalg
just plot-vs-linalg-readme
```

`just bench-vs-linalg` measures current la-stack, nalgebra, and faer rows.
`just plot-vs-linalg-readme` reads those Criterion results and updates:

- `README.md`
- `docs/assets/bench/vs_linalg_lu_solve_median.csv`
- `docs/assets/bench/vs_linalg_lu_solve_median.svg`

See `scripts/criterion_dim_plot.py --help` for plotting options.

### Create The Release Performance Report

Release PRs promote one curated release-to-release comparison into committed
docs:

```bash
just performance-release v0.4.3 v0.4.2
```

With no arguments, `just performance-release` infers the current release tag
from `Cargo.toml` and discovers the previous stable published release. During
release preparation, passing both tags explicitly removes ambiguity.

This command creates temporary worktrees, generates the comparison, writes
`docs/PERFORMANCE.md`, and archives the previous committed report under
`docs/archive/performance/`. Archive filenames are release-pair names such as
`v0.4.2-vs-v0.4.1.md`.

### Compare Published Release Artifacts

After releases are published, the GitHub Release benchmark workflow attaches a
compressed Criterion baseline artifact. To compare those stored artifacts
without running cargo locally:

```bash
just performance-github-assets v0.4.3 v0.4.2
```

With no arguments, the recipe discovers the latest stable published GitHub
release and its previous stable release automatically.

## Output Locations

| Path | Committed? | Producer | Purpose |
|------|------------|----------|---------|
| `target/criterion/` | No | `cargo bench`, `bench-save-*` | Local Criterion measurements and named baselines. |
| `target/bench-reports/performance.md` | No | `bench-compare`, `performance-local*` | Local comparison report. |
| `target/bench-reports/github-assets-performance.md` | No | `performance-github-assets` | Local report from published release artifacts. |
| `docs/PERFORMANCE.md` | Yes | `performance-release` | Latest curated release-to-release comparison. |
| `docs/archive/performance/` | Yes | `performance-release` | Older curated release-to-release comparisons. |
| `docs/assets/bench/` | Yes | `plot-vs-linalg-readme` | README benchmark CSV/SVG assets. |
| GitHub Release asset `la-stack-$TAG-criterion-baseline.tar.gz` | Remote release artifact | `.github/workflows/release-benchmarks.yml` | Durable Criterion baseline archive for published releases. |

## `vs_linalg` Methodology

`vs_linalg` is a per-kernel comparison, not a single aggregate score. Each row
compares one operation for one dimension `D`, using Criterion's selected
statistic from `target/criterion/d{D}/{benchmark}/{sample}/estimates.json`.
The README table uses `median.point_estimate` in nanoseconds. Lower is better.

All three crates receive equivalent deterministic inputs for a given dimension:

- matrix entries come from the same strictly diagonally-dominant generator
  (`matrix_entry::<D>`)
- right-hand sides and vector inputs come from the same deterministic vector
  generator
- each benchmark uses `black_box` around inputs and outputs to keep the
  measured operation visible to the optimizer

The integration smoke test `tests/vs_linalg_inputs.rs` reuses the benchmark
input helpers and verifies that la-stack, nalgebra, and faer agree on the
determinant, solve, dot, and infinity-norm results for D=2..=5. Run it with:

```bash
cargo test --features bench --test vs_linalg_inputs
```

Run that test when changing benchmark input construction, adding comparable
kernels, or updating the `faer` or `nalgebra` benchmark dependencies.

The main comparable metrics are:

- `det_via_lu` — factor the matrix and compute determinant from the LU factor
- `lu` — LU factorization only
- `lu_solve` — factor the matrix and solve one right-hand side
- `solve_from_lu` — solve one right-hand side using a precomputed LU factor
- `det_from_lu` — compute determinant using a precomputed LU factor
- `dot` — vector dot product
- `norm2_sq` — squared Euclidean vector norm
- `inf_norm` — matrix infinity norm, implemented as maximum absolute row sum

Additional SPD metrics compare la-stack LDLT against faer LDLT and nalgebra
Cholesky:

- `ldlt` / `cholesky` — SPD factorization only
- `ldlt_solve` / `cholesky_solve` — factor and solve one right-hand side
- `solve_from_ldlt` / `solve_from_cholesky` — solve using a precomputed factor
- `det_from_ldlt` / `det_from_cholesky` — determinant from a precomputed factor

Read these as SPD factorization/solve/determinant comparisons, not as identical
algorithm comparisons across all three crates.

Release-signal reports compare latest la-stack measurements against a saved
la-stack baseline, and show saved nalgebra/faer baseline timings as context
where a matching peer benchmark exists. That keeps iteration cheap while still
making the release signal auditable. The full `vs_linalg` run remains the source
of README plots and crate-to-crate comparison tables.

## Exact-Arithmetic Notes

The exact suite includes fixed per-dimension groups (`exact_d{2..5}`), random
percentile groups, and adversarial-input groups:

- `exact_random_percentile_d{2..5}` — fixed-seed corpora of 50 strictly
  diagonally-dominant random matrices per dimension. Each operation is
  pre-timed across the corpus to select representative p50/p95/p99 inputs, then
  Criterion measures those inputs normally.
- `exact_near_singular_3x3` — a 2^-50 perturbation of a singular base matrix;
  forces the Bareiss fallback in `det_sign_exact` and exercises the largest
  intermediate `BigInt` values in `solve_exact`.
- `exact_large_entries_3x3` — diagonal entries near `f64::MAX / 2` stress
  `BigInt` growth during Bareiss forward elimination.
- `exact_hilbert_4x4` / `exact_hilbert_5x5` — classically ill-conditioned
  matrices whose non-terminating-in-binary entries stress the
  `f64_decompose -> BigInt` scaling path.

Each random percentile and adversarial group runs the same exact-arithmetic
benches (`det_sign_exact`, `det_exact`, `solve_exact`,
`solve_exact_f64_result`, `solve_exact_rounded_f64`) so tables are comparable
across input classes.

For exact-arithmetic comparisons against v0.4.2 or older baselines, rows such
as `det_exact_rounded_f64 (vs det_exact_f64)` mean the current rounded API is
being compared to the historical lossy `*_exact_f64` benchmark. Rows such as
`det_exact_f64_result (vs det_exact_f64)` intentionally show the overhead of the
new strict conversion contract against that same historical baseline.

The default `release-signal` scope reports exact-arithmetic rows whose inputs
are fixed across versions: deterministic D=2..=5 cases plus adversarial fixed
matrices. Random percentile groups are exploratory tail probes; each benchmark
run selects p50/p95/p99 input sets by timing the implementation under test, so
those rows can measure different corpus subsets across versions. Include them
when investigating tails with:

```bash
uv run bench-compare v0.4.2 --suite exact --scope all-benches
```

To generate a current snapshot without a saved baseline:

```bash
uv run bench-compare --snapshot
```

## Release Notes

Local Criterion baselines are optional during release. Save them only if you
want convenience baselines for follow-up development on the same machine:

```bash
just bench-save-baseline <tag>
just bench-save-last
```

The durable published baseline is the GitHub Release artifact created by
`.github/workflows/release-benchmarks.yml`. The committed release comparison is
`docs/PERFORMANCE.md`, created by `just performance-release`.
