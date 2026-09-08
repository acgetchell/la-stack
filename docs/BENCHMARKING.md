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
  - [Hosted Release Runtime Budget](#hosted-release-runtime-budget)
  - [Validate The Release Workflow](#validate-the-release-workflow)

## Start Here

| Goal | Recipe |
|------|--------|
| Latest-release local audit | `just performance-local` |
| Non-exact release-signal check against tags | `just performance-local-non-exact v0.4.5 v0.4.4` |
| Fast saved-baseline loop | `just bench-save-baseline <name> <suite>` then `just bench-compare <name> <suite> all-benches` |
| Full crate comparison | `just bench-vs-linalg` |
| Interval determinant filter | `just bench-interval` |
| Certified dot/linear-form filter | `just bench-linear-form` |
| README table and plot | `just performance-release` then `just performance-readme` |
| Release report | `just performance-release v0.4.5 v0.4.4` |
| Build docs from retained release inputs | `just performance-doc` |
| Published-asset comparison | `just performance-github-assets v0.4.5 v0.4.4` |

Rule of thumb:

- Use `performance-local*` for clean, self-contained answers.
- Use `bench-save-*` plus `bench-compare` for tight local optimization loops.
- Use `bench-vs-linalg` plus `plot-vs-linalg` for exploratory crate-to-crate
  plots.
- Use `performance-release` only when preparing committed release artifacts.
- After `performance-release`, use `performance-readme` to publish its
  retained crate-to-crate measurements without benchmarking again.
- Use `performance-doc` for report-format changes after a valid, promotable
  comparison dataset has already been retained.

The three canonical workflows compose around one artifact schema, metric set,
and renderer:

| Recipe | Measure | Retain CSV/JSON | Promote release docs |
|--------|---------|-----------------|----------------------|
| `performance-local` | Yes | Yes | No |
| `performance-doc` | No | Consumes retained inputs | Yes |
| `performance-release` | Yes | Yes | Yes |

For a distinct release pair with no intervening source or configuration changes,
running `performance-local` followed by `performance-doc` produces the same
report and committed documentation as `performance-release`.
`performance-release` exists as the safer one-step release operation: it keeps
fresh measurement, validated artifact publication, and rollback-capable document
promotion in one command.

## Benchmark Suites

`la-stack` has five Criterion benchmark suites.

Newly rendered reports use one table per selected suite. Dimension and
adversarial-input group appear in a `Case` column instead of creating a separate
table for every group. The `vs_linalg` table is the one wider variant because it
adds nalgebra and faer context columns where matching peer measurements exist.

**`vs_linalg`** (`benches/vs_linalg.rs`) compares `la-stack` against
`nalgebra` and `faer` across D=2-64 for LU, solve, determinant, dot, norm, and
SPD factorization operations. Use this suite to answer "why choose la-stack over
other crates?"

The SPD rows compare la-stack LDLT, faer LDLT, and nalgebra Cholesky. They are
labelled by algorithm because nalgebra does not expose a dense LDLT
factorization in the dependency version used here.

**`exact`** (`benches/exact.rs`) measures exact-arithmetic methods over lifted
binary64 inputs (`det_exact`, `solve_exact`, `det_sign_exact`, strict `*_result`
conversions, and lossy `*_rounded_f64` conversions) alongside the f64 `det`
baseline across D=2-5. Its supported D=2-4 range also includes `det_direct`, the
paired `det_direct_with_errbound`, and the bound-only `det_errbound`. The same
suite compares row-cleared Bareiss operations with direct `BigRational` Gaussian
operations over already-exact rational inputs across D=2-8. Use it to understand
exact-arithmetic cost and track optimization progress.

**`gram`** (`benches/gram.rs`) compares `gram_matrix` with checked hand-written
assembly for square and embedded vector sets with coordinate dimensions 2-8.
The orthogonal, dependent, near-dependent, and mixed-scale fixtures are checked
against an independent integer matrix-product oracle before timing. Run it with
`cargo bench --locked --features bench --bench gram`. This focused construction
signal is not part of the release-to-release report schema.

**`interval`** (`benches/interval.rs`) measures the default-feature,
division-free interval determinant sign filter. Its fixtures cover a conclusive
4×4 relative-coordinate lifted predicate, the corresponding inconclusive
boundary regime, and a conclusive 7×7 lifted workload at the supported dimension
limit. Fixture construction and expected-sign validation occur outside the
timed closures. This suite is a focused kernel signal; it is not part of the
release-to-release `vs_linalg` or `exact` report schema.

**`linear_form`** (`benches/linear_form.rs`) measures the default-feature
certified dot-product and affine-difference filters at D=4. It compares the
well-separated bounded dot product with the same plain `Vector::dot` input and
also covers dot and affine-difference cases whose certified intervals overlap
zero. Fixture construction and exact small-integer expectations are validated
outside the timed closures. This focused kernel signal is not part of the
release-to-release report schema.

## Common Workflows

### Compare Current Code With The Latest Release

Use this when you want a clean local answer to "how does this checkout compare
with the latest published release?"

```bash
just performance-local
```

This creates isolated temporary worktrees and runs both library revisions on the
same machine with the current checkout's benchmark sources, manifests, lockfile,
benchmark-input tests, recipes, and Rust toolchain. Staged and unstaged changes
to tracked files are applied to the current worktree. Untracked files are
excluded; stage a new file before running the command if it must participate in
the comparison. Only the baseline library implementation comes from the release
tag. Before either timing run, the command runs `just test-bench-inputs` against
that revision under the shared current fixture harness. This is a prerequisite
correctness gate over the deterministic fixtures and operations, not validation
of each timed Criterion sample. It writes
`target/bench-reports/performance.md` plus retained `performance.csv` and
`performance.provenance.json` comparison inputs. The report and sidecar embed
both commits, CPU, operating system, Rust toolchain, lockfile and harness
digests, Criterion selection/commands, and both correctness-gate results. The
report reader rejects malformed or mismatched provenance and incomplete
selected-suite coverage.

The shared harness carries an explicit v0.4.3-only API adapter for renamed or
ownership-adjusted calls (`det_sign_exact`, `Tolerance`, and vector dot
products). The adapter changes only how the same operation is invoked; it does
not patch either library implementation. Comparison builds cap lint diagnostics
at warning for both revisions because the current manifest's lint policy may
reject historical source that predates a lint, even though that source remains
valid benchmark input.

The paired `det_direct_with_errbound` API postdates v0.4.3, so its D=2-4
baseline rows remain explicitly unavailable in that comparison. The older
bound-only `det_errbound` API is source-compatible and remains a required
shared-harness baseline.

The v0.4.3 LU/LDLT balanced-range determinant paths return an incorrect zero,
so their two D=8 stress rows are deliberately not timed as baselines. Reports
leave those baselines explicitly unavailable rather than presenting invalid
performance evidence. The other v0.4.3 D=8 pivoting and ill-conditioned rows
remain in the comparison.

This command does not depend on existing local `target/criterion/` baselines.
It is slower than reusing a saved baseline, but less sensitive to stale local
benchmark state.

The workflow streams its correctness-gate, Cargo, and Criterion output as each
phase runs. `[performance]` markers identify baseline validation, baseline
timing, current validation, and current timing, so a long comparison exposes
completed samples and its active phase instead of remaining silent until the
final report is rendered.

The local report may compare a checkout whose package version is identical to
the latest published release. Commit/ref and source-state provenance distinguish
the modified checkout from the tagged baseline even though both display the same
package version. Release artifact publication remains stricter and requires two
distinct release identifiers.

### Compare Current Code With A Specific Release

For a narrower non-exact check against a known release pair, run:

```bash
just performance-local-non-exact v0.4.5 v0.4.4
```

This generates a local `v0.4.4` `vs_linalg` baseline, measures the current
la-stack `vs_linalg` rows, and renders a `vs_linalg` report. The report includes
saved baseline nalgebra/faer timings as context where matching peer rows exist,
without rerunning current peer crates.

This narrowed peer-context view uses the same metrics and renderer but writes a
separate `performance-non-exact.*` scratch bundle so it cannot replace the
canonical full comparison inputs accidentally.

When tags are provided explicitly, the current tag must match the package
version in the `HEAD` checkout. A mismatch is rejected before tags are fetched,
worktrees are created, or benchmarks run.

### Iterate Against A Local Saved Baseline

Use local saved baselines when tuning one kernel and comparing several edits
against the same starting point. These baselines are local scratch data, not
release artifacts.

For example, before optimizing `Matrix::norm_inf`, save a named baseline:

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

The README benchmark table and SVG plot are derived from the retained release
comparison. Generate and validate that comparison first, then publish its
crate-to-crate measurements:

```bash
just performance-release
just performance-readme
```

`performance-readme` does not run benchmarks. It loads the canonical
`target/bench-reports/performance.csv` and adjacent provenance JSON retained by
`performance-release`, then uses the current la-stack result and the peer
nalgebra/faer results measured by that same shared current harness. It requires
all three timings for every canonical dimension (D=2, 3, 4, 5, 8, 16, 32, and
64) before updating:

- `README.md`
- `docs/assets/bench/vs_linalg_lu_solve_median.csv`
- `docs/assets/bench/vs_linalg_lu_solve_median.svg`
- `docs/assets/bench/vs_linalg_lu_solve_median.provenance.json`

The publisher verifies the retained artifact digest and schema, release version,
measured source state, commit, dependency lock, complete peer coverage, and
recorded measurement provenance before writing anything. When the retained
provenance includes a benchmark-contract digest, it also verifies the benchmark
code, inputs, dependencies, and toolchain against the current checkout. Legacy
retained artifacts without that field remain publishable but are labeled
`legacy-retained-artifact` rather than contract-matched. The derived provenance
sidecar preserves the measurement commands and environment and records the
retained CSV/JSON digests. Missing, stale, or inconsistent input aborts
publication atomically.

For exploratory CSV/SVG output, run `bench-vs-linalg` and `plot-vs-linalg`.
That path still reads raw Criterion output; `--allow-partial` remains
exploratory-only and cannot update README.

See `uv run --locked criterion-dim-plot --help` for plotting options.

### Create The Release Performance Report

Release PRs promote one curated release-to-release comparison into committed
docs:

```bash
just performance-release v0.4.5 v0.4.4
```

With no arguments, `just performance-release` infers the current release tag
from `Cargo.toml` and discovers the previous stable published release. During
release preparation, passing both tags explicitly removes ambiguity.

This command creates temporary worktrees, validates the complete comparison,
and writes the exact report inputs to
`target/bench-reports/performance.csv` with adjacent
`performance.provenance.json`. The CSV records deterministic benchmark keys,
coverage status and notes, baseline/current median estimates, and
same-current-harness nalgebra/faer peer estimates with complete confidence
intervals in nanoseconds. Those peer fields are the source for
`performance-readme`. The JSON sidecar binds the CSV digest and row count to
the release pair, source states, commands, toolchain, Criterion version,
harness/configuration digests, host, and schema version.

Before creating worktrees or running either benchmark revision, structured
local and release-report workflows require an identifiable CPU model. Raw
Criterion benchmark recipes may still record measurements when that metadata
is unavailable, but those measurements cannot be promoted as reproducible
release evidence.

The pair is validated and published before the temporary worktree is removed.
`docs/performance.md` is then rendered from a validated reload of that retained
pair, and the previous committed report is archived under
`docs/archive/performance/`. Archive filenames are release-pair names such as
`v0.4.2-vs-v0.4.1.md`. Serialization, validation, rendering, coverage, or
promotion failures preserve the previous valid report and artifact pair.

To reproduce and promote the report without running Cargo or creating Git
worktrees, use:

```bash
just performance-doc
```

This command fails closed on a missing, partial, malformed, mismatched, or
unsupported artifact pair. It consumes the default CSV/JSON pair retained by a
successful `performance-local` or `performance-release` run, rewrites the
scratch Markdown, promotes it to `docs/performance.md`, and archives the previous
committed report when the release pair changes. Promotion requires distinct
current and baseline package versions, so a same-version local comparison is
retained and reproducible but cannot become release documentation. Use promotion
for presentation-only report corrections; changes to benchmark inputs, code,
toolchains, or measurement configuration require a fresh local or release run.

### Compare Published Release Artifacts

After releases are published, the GitHub Release benchmark workflow attaches a
compressed Criterion baseline artifact. To compare those stored artifacts
without running cargo locally, install the GitHub CLI (`gh`) and authenticate it
with access to the repository (`gh auth login` or an equivalent token). The
requirement applies even when both release tags are supplied explicitly because
the recipe still downloads their GitHub Release assets:

```bash
just performance-github-assets v0.4.5 v0.4.4
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
| `target/bench-reports/performance.md` | No | `bench-compare`, `performance-local`, `performance-release`, `performance-doc` | Canonical local comparison report. |
| `target/bench-reports/performance.csv` | No | `performance-local`, `performance-release` | Validated tabular inputs for the canonical comparison and README publisher. |
| `target/bench-reports/performance.provenance.json` | No | `performance-local`, `performance-release` | Schema, package identifiers, source, command, toolchain, host, digest, and harness provenance consumed by the README publisher. |
| `target/bench-reports/performance-non-exact.*` | No | `performance-local-non-exact` | Narrowed non-exact report and retained peer-context comparison inputs. |
| `target/bench-reports/github-assets-performance.md` | No | `performance-github-assets` | Local report from published release artifacts. |
| `target/bench-reports/github-assets-performance.csv` | No | `performance-github-assets` | Tabular inputs derived from published native archives. |
| `target/bench-reports/github-assets-performance.provenance.json` | No | `performance-github-assets` | Provenance for the published-asset report inputs. |
| `docs/performance.md` | Yes | `performance-release`, `performance-doc` | Latest curated release-to-release comparison. |
| `docs/archive/performance/` | Yes | `performance-release`, `performance-doc` | Older curated release-to-release comparisons. |
| `docs/assets/bench/` | Yes | `performance-readme` | README benchmark CSV/SVG assets and JSON provenance. |
| GitHub Release | Remote | `.github/workflows/release-benchmarks.yml` | Criterion baseline archive. |

Published baseline assets use the filename
`la-stack-$TAG-criterion-baseline.tar.gz`.

Everything under `target/bench-reports/` is reproducible local scratch owned by
the performance-report workflows. It survives temporary-worktree cleanup but
may be removed by `just clean` or `cargo clean`; retain or copy the CSV/JSON pair
while reviewing or re-rendering a release PR. The compact CSV is the analysis and
Markdown-reproduction layer for both local and release comparisons. It does not
replace the full native Criterion `.tar.gz` archive attached to each GitHub
Release, which remains the durable raw
baseline for post-release comparisons.

## `vs_linalg` Methodology

`vs_linalg` is a per-kernel comparison, not a single aggregate score. Each row
compares one operation for one dimension `D`, using Criterion's selected
statistic from `target/criterion/d{D}/{benchmark}/{sample}/estimates.json`.
The README table uses `median.point_estimate` in nanoseconds. Lower is better,
but point-estimate ratios alone are descriptive and do not establish a
statistically supported performance difference. Preserve Criterion confidence
intervals or repeat controlled runs when making a stronger claim.
For experimental-design background on controlled repetitions and uncertainty,
see [REFERENCES.md](../REFERENCES.md) \[13\]; these workflows do not claim to
implement every recommendation in that study.

The harness calls native crate APIs where they expose the same operation. Where
a peer crate does not expose a matching convenience method, repository-owned
adapter code computes the agreed mathematical kernel inside the timed closure:

| Metric family | la-stack implementation | nalgebra implementation | faer implementation |
|---------------|-------------------------|-------------------------|---------------------|
| LU factorization and solve rows | Native `Lu` APIs | Native `LU` APIs | Native partial-pivoting LU APIs |
| LDLT/Cholesky factorization and solve rows | Native `Ldlt` APIs | Native `Cholesky` APIs | Native LDLT APIs |
| `det_via_lu`, `det_from_lu` | Native `Lu::det` | Native `LU::determinant` | Harness adapter: product of the U diagonal and permutation sign |
| `det_from_ldlt` / `det_from_cholesky` | Native `Ldlt::det` | Native `Cholesky::determinant` | Harness adapter: product of the D diagonal |
| `dot` | Native `Vector::dot` | Native `dot` | Harness adapter: left-to-right fused multiply-add loop |
| `norm2_sq` | Native `Vector::norm_squared` | Native `norm_squared` | Native `squared_norm_l2` |
| `norm2` | Native `Vector::norm` | Native `norm` | Native `norm_l2` |
| `inf_norm` | Native `Matrix::norm_inf` | Harness adapter: maximum absolute row sum | Harness adapter: maximum absolute row sum |

The `norm2`, `norm2_sq`, and `inf_norm` benchmark IDs, including scenario suffixes, are
retained for continuity with saved baselines. The current public methods are
`Vector::norm()`, `Vector::norm_squared()`, and `Matrix::norm_inf()`; the historical
benchmark adapters call `norm2_sq()` and `inf_norm()` only when building against
older library releases.

The `norm2` family also measures iterative `f64::hypot` and Delaunay's existing
dimension-specialized scaled implementation as labeled reference kernels. These
adapter timings are benchmark-kernel comparisons, not claims about the speed of
an identically named public convenience method in every crate. Adapter
implementations are versioned with the benchmark harness, included in the
benchmark-contract digest, and covered by the cross-crate input smoke tests.

The ordinary cross-crate norm row uses positive magnitudes in increasing order,
which changes the running scale at every entry in la-stack's one-pass recurrence.
Additional `norm2` rows cover decreasing magnitudes, repeated scales, sparse
vectors, and finite values spanning normal, subnormal, and zero magnitudes. Those
scenario rows compare only la-stack with the two overflow- and underflow-safe
reference kernels; peer methods whose behavior is not contract-equivalent on the
wide-range input are deliberately excluded. Fixture construction and agreement
checks occur before Criterion's measured closures. Run the focused corpus with
`just bench-vs-linalg-quick norm2_scenario`.

All three crates receive equivalent deterministic inputs for a given dimension:

- matrix entries come from the same strictly diagonally-dominant generator
  (`matrix_entry::<D>`)
- right-hand sides and vector inputs come from the same deterministic vector
  generator
- each benchmark uses `black_box` around inputs and outputs to keep the
  measured operation visible to the optimizer
- the `lu_solve` factor-plus-solve comparison uses direct Criterion
  `bencher.iter` measurement; batching overhead is large enough at D=2 to
  distort both the absolute timing and cross-crate ratio
- precomputed-factor benchmarks pass the factor itself through `black_box`
  before each solve or determinant query, preventing invariant captured factors
  from being hoisted out of the measured closure
- `lu_solve` receives owned fixed-size `black_box` inputs inside the measured
  closure, applying the same complete-operation protocol to la-stack and
  nalgebra
- borrowed operations receive references through `black_box`; in particular,
  `Matrix::norm_inf()` does not copy the matrix inside the measured closure

Use `iter_batched` only when fixture construction is explicitly outside the
scientific quantity being measured. The exclusion must be symmetric across the
compared implementations, documented beside the benchmark, and checked against
direct `iter` in the same binary to show that batching does not materially alter
the reported kernel time or cross-crate ratio.

The integration smoke test `tests/vs_linalg_inputs.rs` reuses the benchmark
input helpers and verifies that la-stack, nalgebra, and faer agree on the
determinant, solve, dot, Euclidean-norm, squared-norm, and infinity-norm results for every
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

Local `lu_solve_{pivoting,dense_ill_conditioned}_d{D}` groups cover the same
eight README dimensions. They measure both complete LU solves and solves from
precomputed factors for all three libraries. The pivoting input cyclically
rotates the baseline matrix rows; the dense input is `J + 2^-16 I`, where `J`
is the all-ones matrix, with 2-norm condition number `1 + D * 2^16`. Both use
manufactured solution references and scaled-residual checks before timing.
The pivoting RHS uses rounded dot products; the dense RHS is exactly
representable. These diagnostic calculations provide no rigorous absolute
rounding-error bound. The groups are
separate from the release-signal registry and README plots.

Run both diagnostic families across all eight dimensions with:

```bash
just bench-vs-linalg '^lu_solve_'
```

The main comparable metrics are:

- `det_via_lu` — factor the matrix and compute determinant from the LU factor
- `lu` — LU factorization only
- `lu_solve` — factor the matrix and solve one right-hand side
- `solve_from_lu` — solve one right-hand side using a precomputed LU factor
- `det_from_lu` — compute determinant using a precomputed LU factor
- `dot` — vector dot product
- `norm2` — overflow- and underflow-safe Euclidean vector norm
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

The [solve finalization decision](performance/solve-finalization.md) records
why #234 retained the existing result construction after testing all eight
README dimensions. It links the regression tests and current benchmark command.

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
- `rational_input_d{2..8}` — already-exact, diagonally-dominant rational
  systems. These compare public `RationalMatrix::det_sign`, `det`, and `solve`
  calls using row-denominator clearing plus integer Bareiss elimination with
  straightforward cubic `BigRational` Gaussian determinant and solve
  references on identical matrices and right-hand sides \[7, 11-12\].
- `rational_input_wide{256,1024}_d{2..8}` — diagnostic rational groups with
  large positive rational row factors. These retain exact solutions and strict
  diagonal dominance while exercising heap-backed numerator and denominator
  storage. See the [row-clearing study](performance/rational-row-clearing.md)
  for fixtures, allocation evidence, focused timing commands, and limitations.

Rational Gaussian determinant and solve references use direct `bencher.iter`,
including their required working copies inside the timed operation. Both
implementations therefore start from borrowed, accepted input and include
workspace preparation and result destruction. Earlier Gaussian measurements
used `iter_batched` and excluded input cloning; do not compare those timings
directly with the current complete-operation measurements. Remeasure both
revisions using the same current harness.

Two additional local diagnostic families are excluded from the release signal:

- `canonical_conversion_*_d{2..5}` isolates strict and rounded conversion of
  already-canonical vectors. The dyadic, non-dyadic, and wide-component inputs
  are joined by the smallest subnormal, negative half-subnormal, and integers
  immediately below and exactly at the overflow midpoint. Before timing, both
  canonical and raw-array conversions must return the known output bits
  (including negative zero) or the exact typed reason and component index.
- `det4_diagnostic_*` checks determinant values and signs against rational
  Gaussian elimination. Dense, sparse, and singular controls are joined by
  positive/negative nonzero near-singular determinants, mixed row exponents,
  and extreme diagonal entries. The singular and near-singular inputs must
  leave the floating-point sign filter inconclusive. Near-singular determinants
  are independently checked against ±2^-50. Mixed-exponent rows are scaled by
  `2^[900, -1074, 700, -526]`; these exact scales preserve the dense determinant
  while including subnormal entries. Large entries use `f64::MAX / 2` on the
  diagonal and one elsewhere, producing an exact result beyond binary64 range.

These complement the existing LU pivoting and dense ill-conditioned solve
groups at D=2, 3, 4, 5, 8, 16, 32, and 64, plus the exact near-singular,
large-entry, and Hilbert families described above. The local adversarial
diagnostics can be smoke-tested without collecting timing estimates:

```bash
cargo bench --locked --features bench,exact --bench exact -- \
  '^canonical_conversion_|^det4_diagnostic_' --test
```

Allocation evidence for canonical conversion uses a separate test executable:

```bash
cargo test --locked --release --features bench,exact \
  --test canonical_conversion_allocations -- --nocapture --test-threads=1
```

The counting allocator is not linked into Criterion timing executables.

The local performance audit on 2026-09-05 retained two exact-arithmetic
optimizations: strict conversion uses the existing proof that rational vectors
are canonical, and dense exact D=4 determinants share six lower-row minors.
Matching before/after runs supported both changes. LU loop-expansion and
combined-check prototypes were removed after failing performance acceptance
across the required dimensions.

Those experiments used Rust 1.98.0 on an Apple M4 Max / AArch64 with Criterion
0.8.2, 50 samples, one-second warm-up, and three-second measurement. Correctness
gates ran before timing. The benchmarks above remain runnable with the current
toolchain; historical timing estimates, source fingerprints, and discarded
prototypes are optional local analysis artifacts rather than prerequisites for
building, testing, or benchmarking the crate.

The f64-input random-corpus and adversarial groups run the same exact-arithmetic
benches (`det_sign_exact`, `det_exact`, `solve_exact`,
`solve_exact_f64_result`, `solve_exact_rounded_f64`) so tables are comparable
across input classes.

Before timing begins, every fixed, adversarial, and corpus input is consumed into
a private-field `ValidatedExactInput` after checks by an independent exact
oracle. Timed and registration helpers accept only that proof-bearing wrapper. A
factorial-time Leibniz determinant over exact rational reconstructions verifies
determinant values and signs, including each direct determinant's certified
absolute bound; exact residuals verify `A x = b`; and
strict/rounded binary64 results are checked for their exact bits, typed reason,
and first failing component. These checks run outside timed Criterion closures.
Any disagreement or unexpected error fails setup instead of becoming an
artificially fast measurement.

The original `rational_input_d{2..8}` groups are part of the canonical exact
release signal. The wide-component
diagnostics remain local Criterion evidence and are not included in the release
report's group registry.
Releases produced with the rational-input harness include their Criterion point
estimates and confidence intervals. When the comparison baseline predates the
rational-input API, the report retains current-only rows with an explicit
coverage note and does not calculate a cross-release ratio; historical reports
whose shared harness predates these groups omit them on both sides.

The proof-bearing fixture is therefore a prerequisite correctness gate, not a
claim that every timed sample is revalidated. Criterion closures remain free of
oracle work so their measurements cover only the named operation; the operation
is deterministic for the already-validated input.

The original v0.4.2-to-v0.4.3 report is not library-only evidence for its
headline `det_sign_exact` rows. The v0.4.3 harness routed those calls through an
operation enum and constructed a complete matrix/RHS input inside each timed
iteration, while the saved v0.4.2 samples used direct closures. The current
shared harness borrows a prevalidated input for both revisions. It also verifies
that the D=2–4 headline fixtures resolve through the floating-point filter, so
those rows continue to measure the intended common path. Use a current
shared-harness comparison before attributing the historical D=2 or D=3 changes
to the library implementation.

For exact-arithmetic comparisons against v0.4.2 or older baselines, rows such
as `det_exact_rounded_f64 (vs det_exact_f64)` mean the current rounded API is
being compared to the historical lossy `*_exact_f64` benchmark. Rows such as
`det_exact_f64_result (vs det_exact_f64)` intentionally show the overhead of the
new strict conversion contract against that same historical baseline.

The default `release-signal` scope includes the canonical exact-arithmetic
groups because their inputs and execution order are fixed across revisions. Historical
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
comparison is `docs/performance.md`, created by `just performance-release`.

Follow [Releasing](RELEASING.md#5-create-the-draft-github-release): create the
tagged stable release as a draft, dispatch the workflow with
`--ref "$TAG" -f tag="$TAG"`, and let the workflow upload and verify the archive before it
publishes the draft. Dispatch requires a stable `vX.Y.Z` tag and exactly one
mutable draft with that tag as its title. The producer checks out the resolved
tag commit only when it matches the workflow's own commit. The dispatch ref
must be that same tag, which keeps execution in the tag's cache scope; the
publisher rechecks the commit and captured release ID.
Missing releases, prereleases, and published releases are rejected before
benchmarking. Publication makes the attached evidence immutable.

### Hosted Release Runtime Budget

The producer runs full `vs_linalg` and `exact` suites sequentially on one
`ubuntu-latest` runner. Separate steps allow 150 and 90 minutes respectively;
the outer job allows 285 minutes. Checkout has a 2-minute timeout, followed
by a composite preparation step with one shared 28-minute timeout covering
all tool installation, input validation, and inventory. These two limits bound
setup execution to 30 minutes; exceeding either fails the producer before
benchmarking or publication. The remaining 15 minutes provide headroom for
dataset validation, packaging, upload, diagnostics, and runner overhead.
Discovery compiles both suites before measurement and shares the preparation
deadline with the preceding work. Dependency caches remain disabled.

The [v0.4.5 run](https://github.com/acgetchell/la-stack/actions/runs/32444040827)
measured 304 comparative benchmarks in about 55m 34s, plus 2m 9s compilation
(57m 43s total). Its exact suite was cancelled during `exact_d2/det_exact`;
that run provides no complete exact-suite runtime. The old 60-minute outer
limit therefore could not accommodate even that smaller harness.

Inventory on 2026-09-07 found the following current workloads. The planning
estimate uses 12 seconds per benchmark, allowing analysis/report overhead
above the configured 3-second warmup and 5-second measurement target. This is
a capacity estimate, not a measured runtime or a guaranteed upper bound;
Criterion can extend measurement for expensive iterations.

| Suite | Benchmarks | Nominal warmup + measurement | Planning estimate | Step limit | Headroom above estimate |
| --- | --- | --- | --- | --- | --- |
| `vs_linalg` | 533 | 71.1 min | 106.6 min | 150 min | 43.4 min (41%) |
| `exact` | 264 | 35.2 min | 52.8 min | 90 min | 37.2 min (70%) |

Both commands retain Criterion's release defaults: 100 samples, 100,000
resamples, 95% confidence, 3-second warmup, and 5-second measurement target.
There are no quick-mode flags, sampling reductions, or benchmark filters.
Diagnostic families and all peer rows remain included. The workflow logs
the current inventory counts and planning estimates before timing; its final
summary records each suite's outcome, elapsed seconds, and budget even after
a step failure. A runner-level termination can still prevent that summary.

`just bench-release-inventory` requires a fresh Criterion directory and obtains
the complete IDs from each compiled binary's `--list` output. It also checks
that the inventory includes the report registry and every canonical README
peer row. `just bench-release-check <tag>` then requires every discovered ID,
valid mean/median estimates with 95% intervals, and 100 finite positive samples.
Each named baseline's four raw JSON files must match its `new` measurement.
Missing diagnostics, failed Criterion writes, stale baselines, and malformed
measurements all stop publication. Only successful validation permits packaging
the single `criterion/` archive, including the inventory manifest, and uploading
the temporary Actions artifact. The separate publisher attaches that archive
as `la-stack-$TAG-criterion-baseline.tar.gz` to the draft, verifies its uploaded
state, size, and SHA-256 digest, and only then publishes the release.

### Validate The Release Workflow

Run the applicable local gates for workflow changes:

```bash
just lint-config
just python-ci
just markdown-ci
just doc-check
```

The Python suite executes the workflow's shell with simulated GitHub API
responses to test draft rejection, commit changes, upload failures, asset
verification, and safe reruns. Archive fixtures check the complete dataset.

Every hosted dispatch now requires a real release draft and authorizes its
publication; there is no producer-only manual mode. Use the
[release sequence](RELEASING.md#6-run-benchmarks-and-publish-the-draft) for a
planned release. Record the successful run URL and both elapsed suite times
when validating a budget change. The estimates above still require this
representative hosted run; local tests and archive fixtures do not establish
GitHub-runner runtime or upload success.

The temporary artifact is named
`bench-baseline-$TAG-<run-id>-<producer-attempt>` and retained for 30 days.
Retry a failed producer by rerunning all jobs or dispatching again, so the
draft is checked in the same attempt before setup or measurement starts.
Rerunning only failed publisher jobs reuses the successful producer's artifact.
An existing draft asset is reused only when its bytes match; conflicting assets
require inspection and manual removal while the release is still a draft.
Published releases are always rejected, and `--clobber` is never used. See
[failed-run recovery](RELEASING.md#recovering-a-failed-run) before retrying.
