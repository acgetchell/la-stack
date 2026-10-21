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
- [Output Locations](#output-locations)
- [`vs_linalg` Methodology](#vs_linalg-methodology)
- [Exact-Arithmetic Notes](#exact-arithmetic-notes)
- [Release Notes](#release-notes)
  - [Hosted Release Runtime Budget](#hosted-release-runtime-budget)
  - [Validate The Release Workflow](#validate-the-release-workflow)

## Start Here

Run `just performance-prepare` before measurements to fetch the release tags.
This explicit Git operation and the measurement recipes' temporary worktrees
require authorization. Both labels of new measurements must be v0.5.0 or newer.

| Goal | Recipe |
| --- | --- |
| Working tree versus latest release | `just performance-local` |
| Non-exact working tree versus latest release | `just performance-local-non-exact` |
| Two existing tags, full or non-exact | `just performance-measure "$CURRENT_TAG" "$PREVIOUS_TAG" full` |
| Fast saved-baseline loop | `just bench-save-baseline <name> <suite>`, then `just bench-compare <name>` |
| Native comparison suite | `just bench-vs-linalg` |
| Interval determinant filter | `just bench-interval` |
| Certified dot/linear-form filter | `just bench-linear-form` |
| Measure and promote a release report | `just performance-release` |
| Promote an explicit full scratch run | `just performance-promote` |
| Render the retained release report offline | `just performance-doc` |

`performance-local*` creates scratch evidence. `performance-promote` selects
that scratch pair explicitly and requires a distinct release pair.
`performance-doc` selects archived evidence explicitly. Missing or incomplete
scratch cannot silently select history. `performance-release` composes fresh
full measurement and promotion; a failed measurement stops before promotion.

README and plot publication require independent current/baseline tag and full
commit assertions. The commands are shown below. Keep those assertions from the
reviewed source checkouts or workflow inputs, rather than reading them from the
evidence being validated.

## Benchmark Suites

`la-stack` has seven Criterion benchmark suites.

**`angle`** (`benches/angle.rs`) measures unsigned angles in ambient dimensions
3–6 using dense, nearly parallel, nearly antipodal, mixed-scale, and subnormal
inputs. Structured sparse fixtures vary support size and position; dense and
sparse endpoint cases cover identical and opposite directions. Each fixture
checks independent analytical or integer Gram/exterior-product
references before timing. Run it with
`cargo bench --locked --features bench --bench angle`.
The ordinary max-scaled Kahan control computes the same angles on these vetted
fixtures, but lacks the public API's general input and underflow handling.
Its timings describe kernel cost only. Prepared `Vector::angle`, checked slices
through `VectorAngle::angle`, and construction followed by the method have separate
rows so validation and adapter costs remain visible. Every row uses `iter`
with borrowed, preconstructed inputs; only `construct_then_angle` includes
construction. This focused signal is outside the release-report schema.

Complete-run reports use a table for each statistic, with columns for the
selected series. Benchmark identifiers retain the suite, dimension, and input
case. Saved-baseline comparisons and historical reports retain their separate
layouts.

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

The `exact_exponent_spread_row_d{2..5}` diagnostics vary row scales over
0–2,000 binary exponent bits while preserving the exact determinant. The
`exact_exponent_spread_column_d4` control varies column scales, retaining
exponent differences within each row. Both measure exact values and signs;
all four determinant APIs are checked against rational Gaussian elimination
outside timing. These focused diagnostics are separate from release metrics.

**`rational_solve`** (`benches/rational_solve.rs`) measures diagonal,
upper/lower triangular, row-permuted triangular, general sparse, and dense
rational systems through D=8, with both dyadic and non-dyadic coefficients.
It separates setter-based construction, prepared solves, the full runtime
dispatch adapter, and the captured Delaunay skip-zero Gaussian reference
(including input copies). Fixtures verify manufactured exact solutions against
the independent reference before timing. Run it with:

```bash
cargo bench --locked --features bench,exact --bench rational_solve -- \
  --sample-size 30 --warm-up-time 1 --measurement-time 3 --noplot
```

The non-dyadic `permuted_upper` family and row-permuted Hilbert control
reproduce #246. The dyadic dense control is a row permutation of a strictly
diagonally dominant matrix; the sparse control is tridiagonal.
`tests/exact_bench_config.rs` validates every family through D=8 (including
empty and scalar cases). This focused suite stays outside the release-report
schema.

**`gram`** (`benches/gram.rs`) compares `gram_matrix` with checked hand-written
assembly for square and embedded vector sets with coordinate dimensions 2-8.
The orthogonal, dependent, near-dependent, and mixed-scale fixtures are checked
against an independent integer matrix-product oracle before timing. Run it with
`cargo bench --locked --features bench --bench gram`. This focused construction
signal is not part of the release-to-release report schema.

**`interval`** (`benches/interval.rs`) measures the default-feature,
division-free interval determinant sign filter. It separates scalar addition,
multiplication, and squaring from lifted-matrix assembly and prepared determinant
evaluation. Matrices cover 3×3 through 7×7, with sparse lifted, dense conclusive,
and dense inconclusive inputs. Scalar cases include exact values, cancellation,
zero, and underflow. Only explicitly named assembly measurements include matrix
construction; expected-result validation stays outside timing.
This suite is a focused kernel signal; it is not part of the
release-to-release `vs_linalg` or `exact` report schema.

**`linear_form`** (`benches/linear_form.rs`) measures the default-feature
certified dot-product and affine-difference filters through D=2–6. Prepared
operands cover dense, sparse, cancellation, underflow, and overflow cases.
Separate `*_endpoints` cases request both outward endpoints, accounting for
certificate construction and endpoint derivation together.
Repeated projection batches reuse a common vertex at the origin and after an
exact translation, with separate construction-only, prepared, and
construction-inclusive measurements. The original D=4 plain-dot control remains.
Fixture expectations are validated outside the timed closures.
This focused kernel signal is not part of the
release-to-release report schema.

## Common Workflows

### Compare Current Code With The Latest Release

```bash
just performance-prepare
just performance-local
```

This compares the current working tree with the latest stable GitHub release.
The shared snapshot includes tracked changes and nonignored untracked regular
files, and rechecks them before publishing measurement outputs. This replaces
the former tracked-only snapshot policy. Keep generated files under ignored
`target/`; symlinks, junctions and submodules are outside this snapshot contract.

Both sources use the captured current harness, dependency lock and toolchain.
The declarations in `tooling/measurements/full.toml` and
`tooling/performance.just` select the native host target and recorded Rust
code-generation settings. Inherited Rust flag, compiler-wrapper or target
overrides do not customize these recipes. Edit and review the declarations for
a different measurement environment.

The shared engine discovers complete native inventories on both sources, checks
independent required scientific IDs, and runs both correctness gates before any
timing. The gates cover `vs_linalg_inputs` and `exact_bench_config`.
Baseline measurements include nalgebra and faer; current measurements reuse
those baseline-phase peers and measure la-stack and the exact suite. Diagnostic
cases remain in the complete samples even when omitted from selected report rows.

Every case requires 100 samples, mean and median estimates, 95% confidence
intervals and nanoseconds. Host, CPU, architecture, tool versions, dependencies,
source/harness membership and hashes, commands, environments, discovery and
acquisition declarations are retained with the run. Runs with incompatible
required provenance are rejected.

Outputs are `target/bench-reports/full/{run.json,evidence.json,report.md}`.
They survive worktree cleanup. Same-version local comparisons require distinct
source identities and cannot be promoted as release documentation.

For a shorter comparison that omits the exact suite:

```bash
just performance-local-non-exact
```

It uses `tooling/measurements/non-exact.toml` and writes the separate
`target/bench-reports/non-exact/` bundle. It retains the same sampling,
correctness gates and original reference phases.

### Compare Current Code With A Specific Release

The installed explicit-tag mode measures two existing immutable tags:

```bash
just performance-measure "$CURRENT_TAG" "$PREVIOUS_TAG" full
just performance-measure "$CURRENT_TAG" "$PREVIOUS_TAG" non-exact
```

Use `performance-local*` for working-tree changes against the latest release.
The old positional labels on `performance-local-non-exact` and the local
`--suite`/`--scope` adapter flags are retired. Full and non-exact selection now
uses named declarations. To include an open-ended report selection, maintain a
named measurement declaration using the shared `series.select` contract;
selection never weakens complete sampling or the independent required inventory.

### Iterate Against A Local Saved Baseline

```bash
just bench-save-baseline before vs_linalg
# Make the change, then measure again.
just bench-vs-linalg
just bench-compare before
```

This is exploratory saved-Criterion comparison. It includes added and omitted
rows and does not establish common-harness release eligibility. Use
`performance-local` for that gate. Native Criterion reports remain useful for
individual benchmarks.

The hosted release workflow discovers the full inventory, measures each suite
once, saves its requested baseline name and verifies exact saved copies before
packaging. The installed verifier is also available directly:

```bash
just bench-release-check "$CURRENT_TAG" target/criterion
```

The second argument explicitly selects the root-relative Criterion home.
`tooling/criterion-saved.toml` independently declares required scientific IDs
and their canonical storage paths. Fresh and saved samples must contain the
same exact `benchmark.json`, `estimates.json`, `sample.json` and
`tukey.json` bytes. A shortened discovery inventory or matching samples at
the wrong path are errors.

### Update The README nalgebra/faer Table

After promoting a full release comparison, supply the source identities
captured independently before measurement:

```bash
just performance-readme "$CURRENT_TAG" "$CURRENT_COMMIT" "$PREVIOUS_TAG" "$PREVIOUS_COMMIT"
```

The canonical README view is LU solve, median, logarithmic time axis, declared
in `tooling/performance-readme.toml`. It selects archived evidence and checks
release eligibility, current source/harness membership and hashes, package and
report release references, and existing-tag artifact identity. The current
release may be prospective; tag-pinned links become available when its reviewed
tag is pushed.

The installed publisher formats and checks the candidate README before replacing
the document and CSV/SVG together. The complete-run payload and evidence envelope
supply provenance links; new publication does not create a redundant plot JSON
sidecar. Original sidecars remain frozen. The table retains recorded intervals,
a confidence legend, versioned originating series and explicit peer/library
ratio direction.

For plot/CSV-only output, select a named view and input source:

```bash
just plot-vs-linalg "$CURRENT_TAG" "$CURRENT_COMMIT" "$PREVIOUS_TAG" "$PREVIOUS_COMMIT" dot-mean-log retained
just plot-vs-linalg "$CURRENT_TAG" "$CURRENT_COMMIT" "$PREVIOUS_TAG" "$PREVIOUS_COMMIT" lu_solve-median-linear scratch
just plot-vs-linalg "$CURRENT_TAG" "$CURRENT_COMMIT" "$PREVIOUS_TAG" "$PREVIOUS_COMMIT" lu_solve-median-log non-exact
```

View names combine a metric, `mean` or `median`, and `linear` or `log`.
Metrics are `det_from_lu`, `det_via_lu`, `dot`, `inf_norm`, `lu`,
`lu_solve`, `norm2_sq` and `solve_from_lu`. All 32 views have static
declarations under each of `tooling/plots/non-exact/`,
`tooling/plots/retained/` and `tooling/plots/scratch/`. Scratch views consume
the full scratch pair; retained views work after scratch cleanup. Non-exact
views consume the non-exact scratch pair and use separate
`vs_linalg_non_exact_*` output paths. The old separate metric/statistic,
`sample=new` and boolean log arguments are replaced by this named view.
There is no implicit input fallback.

Plots preserve scientific coordinates, nanoseconds, the selected statistic
and original baseline-phase peers. They render evidence without measuring or
probing the current host. Independent identity and measurement-policy checks
still apply. README publication additionally requires live source freshness.

### Create The Release Performance Report

After updating the package version and fetching tags:

```bash
just performance-release
```

With no labels, the installed CLI infers the current package's release pair
using publication chronology. A prospective version compares the working tree
with the latest release; a published version compares existing tags. Explicit
labels select existing tags:

```bash
just performance-release "$CURRENT_TAG" "$PREVIOUS_TAG"
```

The shared release report configuration requires the full measurement policy,
independent eligibility and distinct release labels. Non-exact runs retain
their local report and plot workflow. Successful promotion writes the immutable complete
run, evidence envelope and report under `docs/performance-runs/runs/<id>/`,
and updates the indexes, latest pointer and `docs/performance.md` in a single
publication transaction. Repeated measurements coexist.

To promote an already measured full scratch pair or rerender retained history:

```bash
just performance-promote
just performance-doc
```

Those commands deliberately select different inputs. A partial scratch pair
fails `performance-promote`; it never triggers a retained-history fallback.

Historical reports, plots, CSV/JSON, indexes and immutable complete runs keep
their original bytes and links during this migration. The existing
v0.5.0/v0.4.6 comparison lacks the newer acquisition declaration and is now a
frozen record, with no regeneration requirement. Current publication checks
must not be relaxed and missing provenance must not be invented to replay it.
Future runs acquired with the current declarations support offline rendering.

The metadata-free `performance-github-assets` command and historical benchmark
execution compatibility remain retired. Published release assets remain
downloadable as original records.

## Output Locations

| Path | Commit? | Owner |
| --- | --- | --- |
| `target/bench-reports/full/` | No | Complete full scratch run, envelope and report |
| `target/bench-reports/non-exact/` | No | Separate non-exact scratch run, envelope and report |
| `target/bench-reports/performance.md` | No | Exploratory `bench-compare` report |
| `docs/performance.md` | Yes | Current promoted release report |
| `docs/performance-runs/` | Yes | Immutable complete runs, reports and retained-history indexes |
| `docs/performance/` | Yes | Frozen historical CSV/JSON evidence |
| `docs/archive/performance/` | Yes | Frozen historical reports and maintained study records |
| `docs/assets/bench/` | Yes | Named-view CSV/SVG outputs and frozen original sidecars |
| `target/criterion/` | No | Native raw Criterion output and saved baselines |

The Markdown recipes exclude immutable performance history and the generated
current report from formatting. Keep those exclusions aligned with
`tooling/performance-report.toml`. No historical file is normalized to the
current formatter or regenerated by an adoption check.

## `vs_linalg` Methodology

`vs_linalg` is a per-kernel comparison, not a single aggregate score. Each row
compares one operation for one dimension `D`. Exploratory Criterion reports use
the selected statistic from `target/criterion/d{D}/{benchmark}/{sample}/estimates.json`.
The README table reads median estimates and intervals in nanoseconds from validated
complete-run evidence. Lower is better,
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
| `dot` | Native `Vector::dot` | Native `dot` | Native `linalg::matmul::dot::inner_prod` |
| `norm2_sq` | Native `Vector::norm_squared` | Native `norm_squared` | Native `squared_norm_l2` |
| `norm2` | Native `Vector::norm` | Native `norm` | Native `norm_l2` |
| `inf_norm` | Native `Matrix::norm_inf` | Harness adapter: maximum absolute row sum | Harness adapter: maximum absolute row sum |

Native faer dot products use the `faer_dot_native` benchmark ID. Historical
`faer_dot` samples measured a repository-owned scalar fused multiply-add loop;
report readers never substitute those samples for the native operation. A saved
baseline without `faer_dot_native` has no native faer dot context. The README
chart and table compare LU factorization plus one solve; the detailed release
report also includes dot products and other kernels.

The `norm2`, `norm2_sq`, and `inf_norm` benchmark IDs, including scenario suffixes, are
retained for continuity with saved baselines. The current public methods are
`Vector::norm()`, `Vector::norm_squared()`, and `Matrix::norm_inf()`, which the
benchmark harness calls directly.

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

The integration smoke test `benches/comparison/tests/vs_linalg_inputs.rs` reuses the benchmark
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

The [solve finalization decision](archive/performance/studies/solve-finalization.md) records
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
  storage. See the [row-clearing study](archive/performance/studies/rational-row-clearing.md)
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
New exact comparisons require these groups, with Criterion point estimates and
confidence intervals, in both baseline and current phases. Frozen historical
reports may retain current-only rows with a coverage note and no cross-release
ratio when their baseline predates the rational-input API. Historical reports
whose shared harness predates these groups omit them on both sides; these
historical omissions do not qualify a run for current complete-run publication.

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
shared-harness comparison between supported releases before attributing timing
changes to the library implementation. v0.4.3 and older cannot be rerun with the
current workflow.

Archived reports may show rows such as `det_exact_rounded_f64 (vs det_exact_f64)`
or `det_exact_f64_result (vs det_exact_f64)`. Those used the historical lossy
`*_exact_f64` baseline. New comparisons require matching benchmark identifiers;
the pre-v0.4.3 name substitutions have been removed.

The default `release-signal` scope includes the canonical exact-arithmetic
groups because their inputs and execution order are fixed across revisions. Historical
baselines created before the `exact_random_corpus_d*` names were introduced do
not have comparable full-corpus rows, so those rows appear once both sides of a
comparison provide the stable group.

For a current snapshot without a saved baseline, use Criterion’s native reports
under `target/criterion/`. The retired `bench-compare --snapshot` command is no
longer provided.

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
`--ref "$TAG" -f tag="$TAG"`, and let the workflow upload and verify the archive.
Review the draft, then run `just release-publish "$TAG"` as the maintainer to
start the separate crates.io workflow. Dispatch requires a stable `vX.Y.Z` tag and exactly one
mutable draft with that tag as its title. The producer checks out the resolved
tag commit only when it matches the workflow's own commit. The dispatch ref
must be that same tag, which keeps execution in the tag's cache scope; the
publisher rechecks the commit and captured release ID.
Missing releases, prereleases, and published releases are rejected before
benchmarking. Publication makes the attached evidence immutable.

### Hosted Release Runtime Budget

The producer runs full `vs_linalg` and `exact` suites sequentially on one
`ubuntu-latest` runner. Separate steps allow 150 and 90 minutes respectively;
the outer job allows 350 minutes. Checkout has a 2-minute timeout, followed
by a composite preparation step with one shared 90-minute timeout covering
all tool installation, input validation, and inventory. Preparation reserves
35 minutes for uncached tools, the two declared 25-minute discovery deadlines,
and 5 minutes for input validation. Exceeding either setup deadline fails the
producer before benchmarking or publication. The remaining 18 minutes provide headroom for
dataset validation, packaging, upload, diagnostics, and runner overhead.
Discovery compiles both suites before measurement and shares the preparation
deadline with the preceding work. Dependency caches remain disabled.

The [first v0.5.0 run](https://github.com/acgetchell/la-stack/actions/runs/38998792985)
spent about 27 minutes installing uncached tools, then exhausted its old
28-minute preparation limit during inventory compilation. Neither measurement
suite started. The expanded setup budget accommodates compilation separately
from tool installation without reducing measurement time or sampling.

The [v0.4.5 run](https://github.com/acgetchell/la-stack/actions/runs/32444040827)
measured 304 comparative benchmarks in about 55m 34s, plus 2m 9s compilation
(57m 43s total). Its exact suite was cancelled during `exact_d2/det_exact`;
that run provides no complete exact-suite runtime. The old 60-minute outer
limit therefore could not accommodate even that smaller harness.

The comparative inventory on 2026-09-07 contained 533 cases. The current exact
suite contains 354 cases: 264 original cases plus 90 diagnostic cases retained
by the [row-scaling study](archive/performance/studies/exact-row-scaling.md).
The planning estimate uses 12 seconds per benchmark, allowing analysis/report overhead
above the configured 3-second warmup and 5-second measurement target. This is
a capacity estimate, not a measured runtime or a guaranteed upper bound;
Criterion can extend measurement for expensive iterations.

| Suite | Benchmarks | Nominal warmup + measurement | Planning estimate | Step limit | Headroom above estimate |
| --- | --- | --- | --- | --- | --- |
| `vs_linalg` | 533 | 71.1 min | 106.6 min | 150 min | 43.4 min (41%) |
| `exact` | 354 | 47.2 min | 70.8 min | 90 min | 19.2 min (27%) |

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
state, size, and SHA-256 digest, and leaves the draft ready for review.

### Validate The Release Workflow

Run the applicable local gates for workflow changes:

```bash
just lint-config
just python-ci
just markdown-ci
just doc-check
```

The Python suite executes the workflow's shell with the pinned installed CLI
and a simulated GitHub provider to prove captured-target handoff, commit-change
rejection and preservation of the draft for review. Generic provider failure matrices belong
upstream. Archive fixtures retain scientific coverage and exact raw bytes.

Every hosted dispatch requires a real release draft and attaches verified
evidence for review. Use the
[release sequence](RELEASING.md#6-run-benchmarks-and-review-the-draft) for a
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
