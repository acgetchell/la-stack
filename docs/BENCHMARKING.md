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

| Goal | Recipe |
|------|--------|
| Latest-release local audit | `just performance-local` |
| Release-signal check against tags | `just performance-local-vs-linalg v0.4.3 v0.4.2` |
| Fast saved-baseline loop | `just bench-save-baseline <name> <suite>` then `just bench-compare <name> <suite> all-benches` |
| Full crate comparison | `just bench-vs-linalg` |
| README table and plot | `just plot-vs-linalg-readme` |
| Release report | `just performance-release v0.4.3 v0.4.2` |
| Published-asset comparison | `just performance-github-assets v0.4.3 v0.4.2` |

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
and lossy `*_rounded_f64` conversions) alongside the f64 `det` baseline across
D=2-5 and `det_direct` across its supported D=2-4 range. Use this suite to
understand exact-arithmetic cost and track optimization progress.

## Common Workflows

### Compare Current Code With The Latest Release

Use this when you want a clean local answer to "how does this checkout compare
with the latest published release?"

```bash
just performance-local
```

This creates isolated temporary worktrees and runs both library revisions on the
same machine with the current checkout's benchmark sources, manifests, lockfile,
benchmark-input tests, recipes, and Rust toolchain. Only the baseline library
implementation comes from the release tag. Before either timing run, the command
runs `just test-bench-inputs` against that revision under the shared current
fixture harness. It writes `target/bench-reports/performance.md` and records both
commits, CPU, operating system, Rust toolchain, lockfile and harness digests,
Criterion selection/commands, and both correctness-gate results. The report
reader rejects malformed or mismatched provenance and incomplete selected-suite
coverage.

The shared harness carries an explicit v0.4.3-only API adapter for renamed or
ownership-adjusted calls (`det_sign_exact`, `Tolerance`, and vector dot
products). The adapter changes only how the same operation is invoked; it does
not patch either library implementation. Comparison builds cap lint diagnostics
at warning for both revisions because the current manifest's lint policy may
reject historical source that predates a lint, even though that source remains
valid benchmark input.

The v0.4.3 LU/LDLT balanced-range determinant paths return an incorrect zero,
so their two D=8 stress rows are deliberately not timed as baselines. Reports
leave those baselines explicitly unavailable rather than presenting invalid
performance evidence. The other v0.4.3 D=8 pivoting and ill-conditioned rows
remain in the comparison.

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
`uv run --locked bench-compare` CLI accepts the explicit `--suite` and
`--scope` flags.

`just bench-save-baseline <name>` writes Criterion samples under
`target/criterion/`. `just bench-save-last` saves the conventional local
baseline named `last`, which enables shortcuts such as:

```bash
just bench-latest-vs-last
just bench-vs-linalg-latest-vs
just bench-compare
```

Saved baselines persist across `git checkout` but not across `cargo clean`, and
they are not pushed to GitHub. A manually saved baseline is comparable only when
its harness has not changed; use `performance-local` for a checked
revision-to-revision comparison.

### Update The README nalgebra/faer Table

The README benchmark table and SVG plot are crate-to-crate comparisons from the
current checkout:

```bash
just plot-vs-linalg-readme
```

This publication recipe validates the benchmark fixtures, runs a fresh complete
`vs_linalg` benchmark, and requires la-stack, nalgebra, and faer results for every
canonical dimension (D=2, 3, 4, 5, 8, 16, 32, and 64) before updating:

- `README.md`
- `docs/assets/bench/vs_linalg_lu_solve_median.csv`
- `docs/assets/bench/vs_linalg_lu_solve_median.svg`
- `docs/assets/bench/vs_linalg_lu_solve_median.provenance.json`

The provenance sidecar records the measured source state, CPU, operating system,
Rust toolchain, dependency lock and harness digests, Criterion dependency and
selection, dimensions, benchmark command, and correctness-gate result. Missing
coverage or provenance aborts publication. Use `--allow-partial` only for
exploratory CSV/SVG output; it cannot update README and its sidecar explicitly
marks measurement provenance unavailable.

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
`v0.4.2-vs-v0.4.1.md`. Publication fails rather than emitting a partial report
when a selected suite or required dimension is absent.

### Compare Published Release Artifacts

After releases are published, the GitHub Release benchmark workflow attaches a
compressed Criterion baseline artifact. To compare those stored artifacts
without running cargo locally:

```bash
just performance-github-assets v0.4.3 v0.4.2
```

With no arguments, the recipe discovers the latest and previous stable
published GitHub releases.

Published artifacts preserve each release's original benchmark harness. Their
historical timing environments may not have been recorded, so report provenance
labels those fields unavailable rather than reconstructing them. The workflow
still runs the current independent fixture gate against both source revisions
under the shared current fixture harness before reading the assets. Use a local
shared-harness workflow before attributing a difference solely to library code.

## Output Locations

| Path | Committed? | Producer | Purpose |
|------|------------|----------|---------|
| `target/criterion/` | No | `cargo bench`, `bench-save-*` | Local Criterion measurements and named baselines. |
| `target/bench-reports/performance.md` | No | `bench-compare`, `performance-local*` | Local comparison report. |
| `target/bench-reports/github-assets-performance.md` | No | `performance-github-assets` | Local report from published release artifacts. |
| `docs/PERFORMANCE.md` | Yes | `performance-release` | Latest curated release-to-release comparison. |
| `docs/archive/performance/` | Yes | `performance-release` | Older curated release-to-release comparisons. |
| `docs/assets/bench/` | Yes | `plot-vs-linalg-readme` | README benchmark CSV/SVG assets and JSON provenance. |
| GitHub Release | Remote | `.github/workflows/release-benchmarks.yml` | Criterion baseline archive. |

Published baseline assets use the filename
`la-stack-$TAG-criterion-baseline.tar.gz`.

## `vs_linalg` Methodology

`vs_linalg` is a per-kernel comparison, not a single aggregate score. Each row
compares one operation for one dimension `D`, using Criterion's selected
statistic from `target/criterion/d{D}/{benchmark}/{sample}/estimates.json`.
The README table uses `median.point_estimate` in nanoseconds. Lower is better,
but point-estimate ratios alone are descriptive and do not establish a
statistically supported performance difference. Preserve Criterion confidence
intervals or repeat controlled runs when making a stronger claim.

All three crates receive equivalent deterministic inputs for a given dimension:

- matrix entries come from the same strictly diagonally-dominant generator
  (`matrix_entry::<D>`)
- right-hand sides and vector inputs come from the same deterministic vector
  generator
- each benchmark uses `black_box` around inputs and outputs to keep the
  measured operation visible to the optimizer
- precomputed-factor benchmarks pass the factor itself through `black_box`
  before each solve or determinant query, preventing invariant captured factors
  from being hoisted out of the measured closure
- consuming stack-matrix inputs are copied in Criterion batch setup, outside
  the measured closure, so factorization rows measure kernels rather than the
  harness's need to reuse one input
- borrowed operations receive references through `black_box`; in particular,
  `inf_norm` does not copy the matrix inside the measured closure

The integration smoke test `tests/vs_linalg_inputs.rs` reuses the benchmark
input helpers and verifies that la-stack, nalgebra, and faer agree on the
determinant, solve, dot, squared-norm, and infinity-norm results for every
measured dimension: D=2, 3, 4, 5, 8, 16, 32, and 64. The same focused recipe
also tests exact-benchmark range and deterministic-generator configuration:

```bash
just test-bench-inputs
```

Run that test when changing benchmark input construction, adding comparable
kernels, or updating the `faer` or `nalgebra` benchmark dependencies.

The D=8 group also includes la-stack stress rows for a forced LU row swap, a
successful diagonal factorization spanning 112 binary exponents, and a balanced
dynamic-range determinant whose sequential factor product leaves the binary64
range even though the final result is one. These rows keep pivoting,
ill-conditioning, and scaled-product cold paths visible alongside the shared
well-conditioned peer fixture.

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

The exact suite includes fixed per-dimension groups (`exact_d{2..5}`), fixed
random-corpus groups, and adversarial-input groups:

- `exact_random_corpus_d{2..5}` — fixed-seed corpora of 50 strictly
  diagonally-dominant random matrices per dimension. Every Criterion iteration
  executes the full corpus in its stable order, so baseline and current
  revisions receive identical workloads. Criterion reports time per complete
  50-input corpus and records throughput in elements.
- `exact_near_singular_3x3` — a 2^-50 perturbation of a singular base matrix;
  forces the direct `BigInt` fallback in `det_sign_exact` and exercises an
  ill-conditioned exact solve.
- `exact_large_entries_3x3` — diagonal entries near `f64::MAX / 2` stress
  `BigInt` growth during Bareiss forward elimination.
- `exact_hilbert_4x4` / `exact_hilbert_5x5` — classically ill-conditioned
  matrices whose binary64 entries have varied mantissas and exponents, stressing
  the `decompose_f64 -> BigInt` scaling path.

Each random-corpus and adversarial group runs the same exact-arithmetic
benches (`det_sign_exact`, `det_exact`, `solve_exact`,
`solve_exact_f64_result`, `solve_exact_rounded_f64`) so tables are comparable
across input classes.

Before timing begins, every fixed, adversarial, and corpus input is consumed into
a private-field `ValidatedExactInput` after checks by an independent exact
oracle. Timed and registration helpers accept only that proof-bearing wrapper. A
factorial-time Leibniz determinant over exact rational reconstructions verifies
determinant values and signs; exact residuals verify `A x = b`; and
strict/rounded binary64 results are checked for their exact bits, typed reason,
and first failing component. These checks run outside timed Criterion closures.
Any disagreement or unexpected error fails setup instead of becoming an
artificially fast measurement.

For exact-arithmetic comparisons against v0.4.2 or older baselines, rows such
as `det_exact_rounded_f64 (vs det_exact_f64)` mean the current rounded API is
being compared to the historical lossy `*_exact_f64` benchmark. Rows such as
`det_exact_f64_result (vs det_exact_f64)` intentionally show the overhead of the
new strict conversion contract against that same historical baseline.

The default `release-signal` scope includes all exact-arithmetic groups because
their inputs and execution order are fixed across revisions. Historical
baselines created before the `exact_random_corpus_d*` names were introduced do
not have comparable full-corpus rows, so those rows appear once both sides of a
comparison provide the stable group.

To generate a current snapshot without a saved baseline:

```bash
uv run --locked bench-compare --snapshot
```

## Release Notes

Local Criterion baselines are optional during release. Save them only if you
want convenience baselines for follow-up development on the same machine:

```bash
just bench-save-baseline <tag>
just bench-save-last
```

The durable published baseline is the GitHub Release artifact created by
`.github/workflows/release-benchmarks.yml`. That workflow runs the benchmark-input
correctness gate before timing or packaging the artifact. The committed release
comparison is `docs/PERFORMANCE.md`, created by `just performance-release`.
