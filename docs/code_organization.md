# Code Organization

Module, feature, and file-placement guidance for the published Rust library and
its unpublished comparison-benchmark package.
Read [AGENTS.md](../AGENTS.md) for agent rules and use the
[contributor workflow](../CONTRIBUTING.md) for setup and validation commands.

## Contents

- [Library modules](#library-modules)
- [Feature boundaries](#feature-boundaries)
- [Tests and examples](#tests-and-examples)
- [Benchmarks and support tooling](#benchmarks-and-support-tooling)
- [Documentation owners](#documentation-owners)

## Library modules

`src/lib.rs` wires private implementation modules and public re-exports. It
also owns the documentation-only `guide` module, the private README doctest
mirrors, and the public prelude. There is no `src/main.rs`.

| Module | Owns |
|--------|------|
| [`src/angle.rs`](../src/angle.rs) | `VectorAngle` slice extension trait, direction scaling, and the shared unsigned-angle kernel |
| [`src/error.rs`](../src/error.rs) | `LaError` and typed singularity, non-finite, positive-semidefinite, tolerance, factorization, arithmetic-operation, and exact-conversion categories |
| [`src/exact.rs`](../src/exact.rs) | Exact determinants and solves, determinant-sign filtering, and exact-to-`f64` conversion |
| [`src/gram.rs`](../src/gram.rs) | Fixed-size Gram construction from vector dot products |
| [`src/interval.rs`](../src/interval.rs) | Interval scalars, matrices, and determinant-sign certification |
| [`src/ldlt.rs`](../src/ldlt.rs) | Unpivoted `Ldlt<D>` factorization for exactly symmetric positive-definite matrices, solves, and determinants |
| [`src/lu.rs`](../src/lu.rs) | Partially pivoted `Lu<D>` factorization, solves, and determinants |
| [`src/matrix.rs`](../src/matrix.rs) | `Matrix<D>`, finite storage, accessors, matrix operations, norms, and direct determinant APIs |
| [`src/norm.rs`](../src/norm.rs) | Exact range-boundary fallback for the approximate Euclidean norm |
| [`src/rational.rs`](../src/rational.rs) | `RationalMatrix<D>` and `RationalVector<D>`, canonical rational inputs, and denominator clearing for exact arithmetic |
| [`src/rounding.rs`](../src/rounding.rs) | Shared binary64 rounding primitives for certified arithmetic |
| [`src/scaled_product.rs`](../src/scaled_product.rs) | Allocation-free scaled products for floating-point factor diagonals |
| [`src/tolerance.rs`](../src/tolerance.rs) | Validated singular-tolerance policy |
| [`src/vector.rs`](../src/vector.rs) | `Vector<D>`, finite storage, reductions, norms, and certified scalar results |

Keep shared arithmetic invariants in the lowest owning module and expose public
items through `src/lib.rs`. Mathematical derivations and algorithm references
belong in [Mathematical basis](mathematical_basis.md) and
[REFERENCES.md](../REFERENCES.md).

## Feature boundaries

The [Cargo manifest](../Cargo.toml) owns feature and dependency declarations.

- **`exact`** gates `src/exact.rs`, `src/rational.rs`, their additional tests,
  exact guides, and exact-arithmetic examples. It enables `det_exact()`,
  strict `det_exact_f64()`, rounded `det_exact_rounded_f64()`,
  `det_sign_exact()`, `solve_exact()`, strict `solve_exact_f64()`, and
  rounded `solve_exact_rounded_f64()` for finite binary64 inputs, plus the
  rational-input types.
- `det_sign_exact()` is infallible for finite-by-construction `Matrix` values.
  Exact-value, conversion, and solve APIs retain their genuine scale,
  representation, and singularity failures.
- `ExactF64Conversion` converts an already-computed determinant or solution
  under the strict or rounded contract without rerunning exact elimination.
  Feature-gated re-exports include `DeterminantSign`, `ExactF64Conversion`,
  `RationalMatrix`, `RationalVector`, `BigInt`, `BigRational`, `FromPrimitive`,
  `ToPrimitive`, and `Signed`.
- Typed error categories, including `UnrepresentableReason`, remain available
  without `exact`; matching errors must not require optional arithmetic
  dependencies.
- The scaled `BigInt` determinant core uses direct expansions for D≤4 and
  fraction-free Bareiss elimination for D≥5. Exact signs add a floating-point
  filter for D≤4. Exact solves use fraction-free forward elimination with
  first-non-zero pivoting and integer back-substitution over a common denominator,
  reducing each final component into `BigRational`. Rational inputs
  use direct substitution for triangular systems (including row permutations),
  and otherwise clear denominators before reusing the integer backend.
- **`bench`** is a cfg-only gate for benchmark targets. Benchmark libraries
  remain dev-dependencies. Nalgebra and faer belong only to the unpublished
  `la-stack-comparison` package under `benches/comparison`, so exact benchmark
  builds do not compile either peer library.

## Tests and examples

- Rust unit tests live in inline `#[cfg(test)]` modules in `src/*.rs`.
- `tests/proptest_*.rs` covers angles, matrix, vector, factorization, exact, rational,
  interval, and Gram properties. Exact and rational suites require `exact`.
- Other `tests/*.rs` suites cover regressions, API contracts, allocation
  behavior, conversion boundaries, and benchmark inputs. Shared property-test
  configuration lives in `tests/common/`.
- `tests/vector_angles.rs` owns independent analytical angle regressions,
  typed operand/shape errors, and allocation evidence. `tests/angle_bench_inputs.rs`
  runs the independent fixture gates shared with the angle benchmark.
- `tests/rational_structured_solve.rs` checks triangular solves against the
  captured Gaussian reference, exact residuals, and typed singularity metadata
  through D=8, including arbitrary row permutations and wide rational inputs.
- `benches/comparison/tests/vs_linalg_inputs.rs` checks the shared comparison
  fixtures. `just test-bench-inputs` and the full CI test pass include both
  workspace packages.
- `tests/semgrep/` contains deliberate static-analysis fixtures, not ordinary
  runnable API examples.
- Python tests live under `tests/tooling/` and run with pytest through
  `just test-python` or the aggregate validation workflows.
- `examples/` contains complete runnable workflows. README mirrors and public
  guide doctests live in `src/lib.rs`.
- `tooling/examples.toml` owns the Cargo example package, feature overrides and
  time budgets. The shared runner discovers every binary example from metadata;
  executable numerical assertions stay in `examples/`.
- `examples/vector_angles.rs` demonstrates fixed-vector and borrowed-slice angles,
  independent magnitude scaling, and subnormal separations.

Use [Testing guidance](dev/testing.md) for dimension coverage and test design,
and [Documentation guidance](dev/docs.md) for executable-example ownership.
Standalone integration and support-tool tests belong under the owning package's
`tests/` directory; inline Rust unit tests and doctests stay with their code.
`tests/__init__.py` keeps the `tests.tooling` package distinct from the
production `tooling` adapters during pytest collection.

## Benchmarks and support tooling

`benches/public_api.rs` owns public construction, storage, mutation, symmetry,
direct determinant certificates, interval determinant and exact scalar conversion
workloads at D=2–5 where applicable. Kernel operations keep their existing owners.
`benches/` contains Criterion suites for angles, exact arithmetic, Gram construction,
intervals, linear forms, structured rational solves, and nalgebra/faer comparisons. Helpers under
`benches/common/` own fixture and oracle validation. Exact benchmark helpers
accept only `ValidatedExactInput`, after independent validation outside timing.
Adversarial groups include near-singular, large-entry, and Hilbert inputs.
`benches/common/angle.rs` owns the independently checked angle fixtures and
the stable comparison kernel; `benches/angle.rs` separates prepared vectors,
borrowed slices, and construction costs.
`benches/common/rational_solve.rs` owns the structured rational fixtures and
captured Delaunay Gaussian control. `benches/rational_solve.rs` separates
construction, prepared solves, and full runtime-dispatch adapter costs;
`tests/exact_bench_config.rs` includes its independent fixture gates.

The root package is the default workspace member. `benches/comparison/Cargo.toml`
owns the `vs_linalg` target at `benches/vs_linalg.rs` and its input test. Select
it with `just bench-vs-linalg` or Cargo's `-p la-stack-comparison`. Both packages
inherit the root workspace lints and share the lockfile and target directory.

The [Benchmarking guide](BENCHMARKING.md) owns benchmark commands, methodology,
baselines, output locations, and report promotion. The [Tooling guide](../tooling/README.md)
owns the performance declaration inventory, consumer-test scope, and recipe ownership.
The pinned published `research-repo-tools` dependency
owns changelog generation, normalization, minor-series archiving, note lookup,
and tag preparation through its CLI. Consumer identity and dependency-body policy
stay in `[tool.research-repo-tools.changelog]` in `pyproject.toml`, with
formatting in `changelog-rumdl.toml`; no common template is copied locally.
Focused integration checks live in `tests/tooling/test_changelog_integration.py`.
`tooling/release_changelog.py` temporarily resolves omitted prospective release
dates from the prepared citation metadata before delegating generation and
preview to the installed CLI. It does not update citation metadata or use a clock.
The same dependency owns CodeRabbit review orchestration through thin Just
wrappers. `tests/tooling/test_review_integration.py` owns consumer wiring checks
with local stubs; the [contributor review workflow](../CONTRIBUTING.md#coderabbit-review)
owns prerequisites and invocation policy.
The [justfile](../justfile) owns executable development workflows.
Shared complete-inventory Python and file runners own discovery and bounded
execution; repository lint policy and archive exclusions stay in the recipes
and native tool configuration. The [capability inventory](../tooling/README.md#shared-capability-adoption)
records the v0.2.3 adoption boundary.

`osv-scanner.toml` owns temporary advisory-specific dependency exceptions, and
`.gitleaks.toml` owns narrow secret-scan false-positive exceptions. Their rationale
and review policy belong in [Security Checks](../SECURITY.md#security-checks).

The shared package also owns managed tool installation, verification, and update
implementation. `.python-version`, `rust-toolchain.toml`, and `pyproject.toml`
own consumer declarations; `tests/tooling/test_toolchain_integration.py` checks
the exact published pin, registry lock, installed version and distribution hashes,
alongside actual recipe sequencing and native managed execution. The
`.github/actions/setup-tools/action.yml` composite synchronizes the locked PyPI
package and exports verified paths for CI. Release callers disable its caches.

`.github/actions-allowlist.json` owns approved action identities for local
validation and updates. Historical payloads under `.github/settings/` describe
live GitHub settings separately; consult the
[rollout guide](dev/MANAGING_CHANGES.md#dependabot-approval-rollout) before applying them.

Release metadata, version checks, Markdown line checks, and Semgrep fixture
validation also belong to the shared CLI. Consumer policy stays in
`pyproject.toml`; `tests/tooling/test_maintenance_integration.py` verifies the
actual release selectors and preservation of scientific evidence.
`tests/tooling/test_cargo_update_integration.py` exercises native dependency
upgrades and coupled exclusions against a disposable local registry.

The shared package also owns process discovery, captured/live execution, exact
byte transport, CPU metadata, diagnostics, and zizmor authentication. Consumers
call its public process API directly. Shared Criterion estimates replace local
timing wrappers. Tests use shared Just inspection and retain scientific
selection and consumer policy checks.

Performance consumers use shared Criterion parsing and estimate/comparison
validation, digest verification, archive extraction, byte-preserving document
sections, and multi-file transactions. Shared rendering produces complete-run
reports, coordinate SVG, timing CSV and document tables. Historical artifact
bytes and fingerprint framing remain frozen and readable directly. Published v0.2.2
owns common-harness orchestration, completeness, worktrees, complete run retention,
validated latest selection, and transactional report publication.
Its `RunTableStyle` API owns compact README table formatting; the consumer selects
dimension headings, intervals below timings, a common confidence legend and
ratio-only speedups without rewriting Markdown cells or prose.
Generic parsing, staging, and rollback tests belong upstream; local
tests verify the scientific and retained-artifact integration boundaries.

Hosted Dependabot approvals use the pinned shared
GitHub workflow; the [rollout guide](dev/MANAGING_CHANGES.md#dependabot-approval-rollout)
owns settings and deployment verification.

`tooling/criterion-inventory.toml` declares hosted la-stack-only list commands and
independent public API coverage. `tooling/criterion-saved.toml` declares
complete saved-sample requirements. Shared discovery and saved-sample APIs own
inventory parsing and exact raw-file verification. The release workflow packages only datasets that pass
that gate; its regression and archive tests live in
`tests/tooling/test_release_baseline.py`.

`tooling/performance.just` declares native Cargo discovery, correctness gates
and measurements. `tooling/measurements/{full,non-exact}.toml` declares the
minimum release, complete sampling, independent required cases, original
baseline-phase peers, sources, harness, probes, environment and compatibility.
`tooling/performance-policy.toml` remains the independent scientific catalog
used by consumer acceptance tests and the local budget display. Full acquisition
adds only 16 peer LU-and-solve cases in the baseline phase; the non-exact
diagnostic retains broader peers. Runtime commands consume the shared schemas.

`tooling/performance-report.toml` binds release eligibility to the full
measurement declaration and owns the current report and immutable history.
`tooling/performance-readme.toml` declares the canonical README view, formatter,
freshness checks and tagged evidence links. Static named configurations under
`tooling/plots/{retained,scratch}/` preserve LU solve with mean/median and
linear/log choices. `tooling/plots/non-exact/` retains eight diagnostic metrics.
All views require explicit input ownership.

`tooling/release_performance.py` temporarily captures independent release tags and
commits during `performance-prepare` and supplies them to the separate
`performance-readme` publication step, pending equivalent published upstream commands.
The former measurement/publication implementations and root `scripts/` directory are retired.
The installed CLI owns orchestration, validation and publication transactions.
Working-tree measurements now include nonignored untracked regular files.
`tests/tooling/test_performance_workflow.py` exercises actual recipes and
scientific policy. Shared failure matrices remain upstream. The
[adoption record](dev/performance-tooling-adoption.md) records replacement
ownership and native validation gates.

Original artifacts, including the v0.5.0/v0.4.6 complete run, remain frozen
without regeneration or invented acquisition records. New complete runs
support offline rendering through the declared eligibility policy.

The root Python 3.15 environment owns Semgrep alongside the other development
tools. Its exact pin, native-wheel requirement and compatibility overrides live
in `pyproject.toml` and `uv.lock`. Local scans, shared fixture validation and
hosted SARIF use that lock and retain the Linux stack-limit workaround for
[Semgrep #11979](https://github.com/semgrep/semgrep/issues/11979).
The root `.semgrepignore` keeps ordinary Rust and Python tests in those scans
while excluding deliberate `tests/semgrep/` fixtures, which `just semgrep-test`
validates separately. Git-ignored dependencies, environments, caches, and build
outputs remain outside production scans.

`.github/actions/prepare-release-benchmarks/action.yml` groups tool installation,
input validation, and inventory under the release workflow's shared setup timeout.

`.github/workflows/ci.yml` owns native checks, default/exact Cargo dry-runs and
signed main-push validation evidence. `tooling/release_evidence.py` supplies the
Cargo-specific policy around shared run/job validation and artifact retrieval;
`tests/tooling/test_release_evidence.py` owns its provider-boundary regressions.
`.github/workflows/publish-crates.yml` owns the reviewed-release event, reuse of
that exact CI attempt and artifact, deployment approval and OIDC upload. The publishing declaration
in `pyproject.toml` selects `la-stack`, native checks and the Criterion asset;
shared release gates own generic identity, ancestry and registry checks.
`tests/tooling/test_release_publishing.py` owns installed-package consumer checks
and upload-barrier tests. The [release guide](RELEASING.md) owns account setup,
manual publication, registry verification and recovery.

## Documentation owners

Use [Documentation guidance](dev/docs.md) for README, references, mathematical
background, rustdoc, and generated-file ownership. Focused operational rules
live in `docs/dev/`; human setup and validation guidance lives in
[CONTRIBUTING.md](../CONTRIBUTING.md). Release procedures belong in
[Releasing](RELEASING.md).

[Managing changes](dev/MANAGING_CHANGES.md) owns Git and GitHub procedures.
[Measuring coverage](MEASURING_COVERAGE.md) owns local and CI coverage execution.
The generated [performance report](performance.md) records measured results;
[Benchmarking](BENCHMARKING.md) owns the commands that produce and compare them.
[Local benchmark summaries](performance/README.md) owns complete versioned local
datasets. Completed optimization studies belong under
[archived performance studies](archive/performance/studies/README.md).
The [Rust 1.99 study](archive/performance/studies/rust-1.99.md) owns the compiler
migration baseline and its retained measurements; it is separate from generated
release-to-release reports.
The [interval and certified-reduction study](archive/performance/studies/interval-reductions.md)
owns the #247/#248 same-toolchain comparisons and downstream adoption evidence.
The [vector-angle study](archive/performance/studies/vector-angles.md) owns the
#249 formulation decision, comparative timings, and retained source provenance.
The [sparse vector-angle study](archive/performance/studies/sparse-vector-angles.md)
owns the #287 shared-zero optimization, dense controls, and retained evidence.
The [triangular rational solve study](archive/performance/studies/rational-triangular.md)
owns the #246 Rust 1.99 baseline, adapter/kernel comparisons, and retained evidence.
The [exact determinant row scaling study](archive/performance/studies/exact-row-scaling.md)
owns the row-local scaling comparison, exponent-spread diagnostics, and workload tradeoffs.
The [exact solve normalization study](archive/performance/studies/exact-solve-normalization.md)
owns the shared-denominator back-substitution comparison and diagonal controls.
The [algebraic floating-point study](archive/performance/studies/algebraic-float.md)
owns the #250 numerical/performance experiment and its raw evidence.

The [validation index](validation/README.md) maps arithmetic contracts
to independent tests and owns dated validation decisions and their limitations.
The [algebraic floating-point decision](validation/algebraic-floating-point.md)
records rejection and the maintainer's closure of unperformed prototype work.
Mathematical guarantees remain in `mathematical_basis.md`; performance archives
retain measurement data and experimental sources.

When adding, removing, renaming, or moving files, update the applicable ownership
rows here. Prefer links to the detailed owner over copying its procedure into
this map or `AGENTS.md`.
