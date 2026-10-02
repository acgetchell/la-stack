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
  first-non-zero pivoting and `BigRational` back-substitution. Rational inputs
  clear denominators before reusing the integer backend.
- **`bench`** is a cfg-only gate for benchmark targets. Benchmark libraries
  remain dev-dependencies. Nalgebra and faer belong only to the unpublished
  `la-stack-comparison` package under `benches/comparison`, so exact benchmark
  builds do not compile either peer library.

## Tests and examples

- Rust unit tests live in inline `#[cfg(test)]` modules in `src/*.rs`.
- `tests/proptest_*.rs` covers matrix, vector, factorization, exact, rational,
  interval, and Gram properties. Exact and rational suites require `exact`.
- Other `tests/*.rs` suites cover regressions, API contracts, allocation
  behavior, conversion boundaries, and benchmark inputs. Shared property-test
  configuration lives in `tests/common/`.
- `benches/comparison/tests/vs_linalg_inputs.rs` checks the shared comparison
  fixtures. `just test-bench-inputs` and the full CI test pass include both
  workspace packages.
- `tests/semgrep/` contains deliberate static-analysis fixtures, not ordinary
  runnable API examples.
- Python tests live under `scripts/tests/` and run with pytest through
  `just test-python` or the aggregate validation workflows.
- `examples/` contains complete runnable workflows. README mirrors and public
  guide doctests live in `src/lib.rs`.

Use [Testing guidance](dev/testing.md) for dimension coverage and test design,
and [Documentation guidance](dev/docs.md) for executable-example ownership.

## Benchmarks and support tooling

`benches/` contains Criterion suites for exact arithmetic, Gram construction,
intervals, linear forms, and nalgebra/faer comparisons. Helpers under
`benches/common/` own fixture and oracle validation. Exact benchmark helpers
accept only `ValidatedExactInput`, after independent validation outside timing.
Adversarial groups include near-singular, large-entry, and Hilbert inputs.

The root package is the default workspace member. `benches/comparison/Cargo.toml`
owns the `vs_linalg` target at `benches/vs_linalg.rs` and its input test. Select
it with `just bench-vs-linalg` or Cargo's `-p la-stack-comparison`. Both packages
inherit the root workspace lints and share the lockfile and target directory.

The [Benchmarking guide](BENCHMARKING.md) owns benchmark commands, methodology,
baselines, output locations, and report promotion. The [Scripts guide](../scripts/README.md)
owns the Python script inventory and entry points for comparisons and plotting.
The pinned published `research-repo-tools` dependency
owns changelog generation, normalization, minor-series archiving, note lookup,
and tag preparation through its CLI. Consumer policy stays in `cliff.toml`,
`changelog-rumdl.toml`, and `[tool.research-repo-tools]` in `pyproject.toml`;
focused integration checks live in `scripts/tests/test_changelog_integration.py`.
The same dependency owns CodeRabbit review orchestration through thin Just
wrappers. `scripts/tests/test_review_integration.py` owns consumer wiring checks
with local stubs; the [contributor review workflow](../CONTRIBUTING.md#coderabbit-review)
owns prerequisites and invocation policy.
The [justfile](../justfile) owns executable development workflows.

`osv-scanner.toml` owns temporary advisory-specific dependency exceptions, and
`.gitleaks.toml` owns narrow secret-scan false-positive exceptions. Their rationale
and review policy belong in [Security Checks](../SECURITY.md#security-checks).

The shared package also owns managed tool installation, verification, and update
implementation. `.python-version`, `rust-toolchain.toml`, and `pyproject.toml`
own consumer declarations; `scripts/tests/test_toolchain_integration.py` checks
actual recipe sequencing and native managed execution. The
`.github/actions/setup-tools/action.yml` composite synchronizes the locked PyPI
package and exports verified paths for CI. Release callers disable its caches.

Release metadata, version checks, Markdown line checks, and Semgrep fixture
validation also belong to the shared CLI. Consumer policy stays in
`pyproject.toml`; `scripts/tests/test_maintenance_integration.py` verifies the
actual release selectors and preservation of scientific evidence.
`scripts/tests/test_cargo_update_integration.py` exercises native dependency
upgrades and coupled exclusions against a disposable local registry.

The shared package also owns process discovery, captured/live execution, exact
byte transport, CPU metadata, diagnostics, and zizmor authentication. The thin
`scripts/benchmark_process.py` adapter retains benchmark phase signatures and
consumer root selection. Generic `subprocess_utils.py` and `run_zizmor.sh`
implementations and their duplicate tests are retired; scientific schemas and
benchmark policy remain here.

Performance consumers use shared Criterion parsing and estimate/comparison
validation, digest verification, archive extraction, byte-preserving document
sections, and multi-file transactions. Local rendering produces complete candidate
outputs before publication. Historical artifact schemas and fingerprint framing,
benchmark selection and eligibility, common-harness orchestration, and complete
run retention remain in the consumer pending the corresponding shared workflow
contract. Generic parsing, staging, and rollback tests belong upstream; local
tests verify the scientific and retained-artifact integration boundaries.

Hosted Dependabot approvals use the pinned shared
GitHub workflow; the [rollout guide](dev/MANAGING_CHANGES.md#dependabot-approval-rollout)
owns settings and deployment verification.

`scripts/release_baseline.py` owns release-suite inventory and complete raw
Criterion validation. The release workflow packages only datasets that pass
that gate; its regression and archive tests live in
`scripts/tests/test_release_baseline.py`.

`scripts/criterion_measurements.py` owns raw sample validation shared by hosted
archives and local summaries. `scripts/benchmark_summaries.py` owns complete
local summary serialization, snapshot identity, and lookup after cleanup.

`.github/actions/prepare-release-benchmarks/action.yml` groups tool installation,
input validation, and inventory under the release workflow's shared setup timeout.

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

When adding, removing, renaming, or moving files, update the applicable ownership
rows here. Prefer links to the detailed owner over copying its procedure into
this map or `AGENTS.md`.
