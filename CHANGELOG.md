# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.4.4] - 2026-07-13

### ⚠️ Breaking Changes

- Make numerical invariants and errors explicit
- Require Rust 1.97.0

### Added

- [**breaking**] Make numerical invariants and errors explicit
  [`668daed`](https://github.com/acgetchell/la-stack/commit/668daed6f5d83a03fc95a7872a0a302d097e3fc4)

  - add structured error reasons, locations, origins, and factorization context
  - make determinant, exact-conversion, LU, and LDLT paths range-safe and mathematically explicit
  - validate benchmark inputs independently and require reproducible, provenance-backed performance evidence
  - centralize tool versions, adopt nextest profiles, and replace Codacy with repository-owned checks and SARIF reporting

### Changed

- Consolidate setup abort handling [`b9e453e`](https://github.com/acgetchell/la-stack/commit/b9e453e540da3c63b0a14ca94b4bcfbc9c8ff6bd)

  - Share `OrAbort` across exact and comparison benchmarks while preserving contextual panic behavior.
  - Centralize Clippy policy in Cargo.toml for consistent local and CI checks.
  - Keep `unreachable_pub` enforced while allowing intentional crate-visible helpers.
  - Enable CodeRabbit's request-changes review workflow.

### Dependencies

- Bump python-multipart in the uv group across 1 directory [`6bb4fb7`](https://github.com/acgetchell/la-stack/commit/6bb4fb7ce48e49d4d2a5bfe4ff5d55b09e6300e2)
- Bump zizmorcore/zizmor-action from 0.5.6 to 0.5.7 [`4a999a6`](https://github.com/acgetchell/la-stack/commit/4a999a664211126e581da47713f6d2426d6db698)
- Bump github/codeql-action/upload-sarif [`85c0a7b`](https://github.com/acgetchell/la-stack/commit/85c0a7bd6a1b3d3e0e00f4e82c6d544c2ca0f22b)
- Bump actions/checkout from 6.0.2 to 7.0.0 [`0f0f956`](https://github.com/acgetchell/la-stack/commit/0f0f956c64ddbc5700f3c2f96b90242703e0a410)
- Bump astral-sh/setup-uv from 8.1.0 to 8.2.0 [`259acaa`](https://github.com/acgetchell/la-stack/commit/259acaa06cc2d6b50f5c024c60cbdf9fe3172f7f)
- Bump codecov/codecov-action from 6.0.1 to 7.0.0 [`035807c`](https://github.com/acgetchell/la-stack/commit/035807c9312ce1ff818dcab9dd986980dfcb37d8)
- Bump vulnerable Python tooling dependencies [`4990260`](https://github.com/acgetchell/la-stack/commit/49902604db70de7c6831f559aee967bf3c376bc1)

### Documentation

- Add mathematical basis and clarify numerical guarantees [`31b7b1e`](https://github.com/acgetchell/la-stack/commit/31b7b1e7c547fc69edbb601e68b2f29e83d5cd24)

  - define the binary64, factorization, determinant-filter, exact-arithmetic, tolerance, and error contracts
  - align README, API documentation, examples, and LDLT diagnostics with guarantees over stored values
  - audit academic references and extend Criterion coverage for determinant error-bound paths
- Clarify determinant and LDLT failure contracts [`e257343`](https://github.com/acgetchell/la-stack/commit/e2573434c541a8115750b3a965882f7b484b6be0)
  - State the determinant sign condition without ambiguous Markdown delimiters.
  - Distinguish singular zero pivots from coupled zero-pivot LDLT rejections.
  - Align citation and contributor guidance with the positive-definite domain.

### Fixed

- Make v0.4.3 comparisons correctness-aware [`4ac5af9`](https://github.com/acgetchell/la-stack/commit/4ac5af9ce104c5126bfd6217414977cca269b62e)

  - Adapt the shared benchmark harness across v0.4.3 API differences without changing measured operations
  - Exclude invalid balanced-range baselines while requiring current samples and reporting unavailable comparisons
  - Preserve benchmark provenance, suite-specific fallback commands, and publication rollback guarantees
  - Harden Windows Git input, changelog links, and version-reference parsing across platforms
- Harden exact arithmetic and benchmark publication [`db2fad5`](https://github.com/acgetchell/la-stack/commit/db2fad50087c773f7ae37f181f720da2a9ff4a23)
  - add `DeterminantWithErrorBound` for paired determinant estimates and certified bounds
  - scale exact systems independently and round exact values directly to IEEE-754
  - fail benchmark publication closed on invalid samples or mismatched provenance
  - make release and changelog tooling transactional, path-safe, and Windows-portable
  - align benchmark CI with pinned local tools and least-privilege publishing

### Maintenance

- Update setuptools requirement in the dependencies group [`6a6c339`](https://github.com/acgetchell/la-stack/commit/6a6c339e7fbec1c9ca4226005794e8ab87f0c2d6)

  Updates the requirements on [setuptools](https://github.com/pypa/setuptools) to permit the latest version.

  Updates `setuptools` to 83.0.0

  - [Release notes](https://github.com/pypa/setuptools/releases)
  - [Changelog](https://github.com/pypa/setuptools/blob/main/NEWS.rst)
  - [Commits](https://github.com/pypa/setuptools/compare/v65.0.0...v83.0.0)
- Avoid duplicate review checks [`5e15d92`](https://github.com/acgetchell/la-stack/commit/5e15d9250455fd1467e47a735412f65449f30b56)
  - Use the legacy required status while preserving automatic approvals.
  - Add explicit dense D3 determinant coverage and document the D4 error bound.
- [**breaking**] Require Rust 1.97.0 [`4033f1f`](https://github.com/acgetchell/la-stack/commit/4033f1f39e411ed066fec808805371979995964b)
  - Align Cargo, the pinned toolchain, Clippy, dependency-audit triggers, and maintainer guidance with the new baseline.
  - Audit Rust 1.97 integer bit helpers against exact-arithmetic hot paths and retain existing operations where alternatives do not preserve performance.
  - Correct v0.4.3 citation metadata and document exact decomposition and solve scaling accurately.

### Performance

- Restore det_direct throughput for D=2..4 [`f6d6885`](https://github.com/acgetchell/la-stack/commit/f6d68855a314fc4802c559080d9faca85db81bb9)

  - Reintroduce branch-free dense expansions for D=2 and D=3.
  - Share D=4 minors while retaining guarded sparse evaluation.
  - Preserve non-finite handling for mathematically inactive terms.
- Restore small-matrix det_sign_exact throughput [`b8bfa9f`](https://github.com/acgetchell/la-stack/commit/b8bfa9fa0c908ea034d294ceb123c8191eff3395)
  - Reuse proof-bearing shared minors for D4 determinant and permanent evaluation.
  - Restore the dense D3 filter while preserving sparse, overflow, and underflow fallbacks.
  - Require headline benchmarks to exercise the intended filter and document historical harness overhead.
- Improve inf_norm throughput [`dc28f98`](https://github.com/acgetchell/la-stack/commit/dc28f9828a1e90b5604e4a55284175e50f84d7e2)
  - Check row-sum finiteness once per completed row on the success path.
  - Replay only overflowing rows to preserve the first failing matrix coordinate.
  - Preserve left-to-right accumulation and const evaluation across dimensions.
- Improve dot and norm2_sq throughput [`ece54d7`](https://github.com/acgetchell/la-stack/commit/ece54d7b2a73111c2bc378a417189e2951a39cbc)
  - Check accumulator finiteness once after the success-path reduction.
  - Replay only non-finite reductions to preserve the first failing step.
  - Preserve left-to-right fused accumulation and const evaluation.

## [0.4.3] - 2026-06-09

### ⚠️ Breaking Changes

- Make exact f64 conversions strict

### Added

- Guard README dependency snippets [`7137fee`](https://github.com/acgetchell/la-stack/commit/7137fee16ab33e08f4dc6a60e02417e3e7c4e020)

  - Add a generic docs-version sync check that compares Markdown dependency
    snippets against the Cargo package name and version
  - Run the docs-version check from the repository Semgrep policy lane
  - Refresh README determinant examples with explicit fallible handling and
    hidden doctest mirrors
  - Update CI uv pins to 0.11.19
- Report determinant scale overflow precisely [`928f62b`](https://github.com/acgetchell/la-stack/commit/928f62bba0d837afe04cb8ccb3fbfed6b095d8f7)
  - Add a typed LaError::DeterminantScaleOverflow path for exact determinant scale exponent failures
  - Convert det_exact_f64 directly from the shared Bareiss integer/exponent pair while preserving Overflow for finite-f64 conversion failures
  - Reuse vector finiteness scanning across raw and proof-bearing constructors
  - Harden docs version sync checks for reordered inline-table dependency snippets and pruned Markdown traversal
- Add release performance comparison workflow [`53b5fde`](https://github.com/acgetchell/la-stack/commit/53b5fde13f3e4afbd3db80d324185a091203cb75)
  - Extend vs_linalg with LDLT/Cholesky benchmark rows and shared deterministic inputs.
  - Add smoke coverage that checks la-stack, nalgebra, and faer agree on benchmark inputs.
  - Expand bench-compare to support latest-vs-last reports, suite/scope selection, peer baseline context, and clearer malformed Criterion diagnostics.
  - Document the benchmark methodology, release baseline workflow, roadmap direction, and contributor guidance.
- Publish release benchmark baselines [`9497ca5`](https://github.com/acgetchell/la-stack/commit/9497ca5f88dd7800bdf3823123e2a295fcb5ced1)
  - Add a release-only benchmark workflow that saves full Criterion baselines for published releases and attaches the archived baseline to the GitHub Release.
  - Keep the regular benchmark workflow focused on PR and main-branch comparison runs.
  - Document how to restore archived release baselines for future performance comparisons.
- Feat!(api): make Matrix and Vector finite by construction [`1fa2f55`](https://github.com/acgetchell/la-stack/commit/1fa2f55cfac6f249a7e2bf30922901539e580dd8)
- [**breaking**] Make exact f64 conversions strict [`8e33f1a`](https://github.com/acgetchell/la-stack/commit/8e33f1a8ec291bfcb6312375969efce076421e96)
  - Add explicit rounded exact-to-f64 APIs for determinant and solve results
  - Report exact conversion failures with typed Unrepresentable reasons
  - Remove finite proof wrapper APIs now that Matrix and Vector carry finiteness directly
  - Move error and tolerance contracts into first-class modules with prelude exports
  - Update exact benchmarks to distinguish strict Result paths from rounded f64 paths
  - Document and exercise the rounded fallback pattern for RequiresRounding errors
- [**breaking**] Make exact f64 conversions strict [`89f3720`](https://github.com/acgetchell/la-stack/commit/89f3720ecde9f12d7a0f42e79394836615e8fd97)
  - Make Matrix and Vector the finite-by-construction public types for exact arithmetic.
  - Add rounded exact-to-f64 APIs for determinant and solve callers that want explicit lossy conversion.
  - Return typed Unrepresentable reasons when strict exact-to-f64 conversion would round or become non-finite.
  - Specialize D4 exact determinants and keep determinant/error-bound zero coefficients from evaluating overflowing absent terms.
  - Update exact benchmark comparison reporting to compare strict and rounded APIs against legacy v0.4.2 rows.
- Archive release performance reports [`2817d01`](https://github.com/acgetchell/la-stack/commit/2817d01374ad0aeab98d6f48a3dae9b30f878a8a)
  - Add an archive-performance utility that promotes curated benchmark reports into docs/PERFORMANCE.md while archiving prior release comparisons
  - Generate release comparisons in isolated temporary worktrees, including legacy command fallback for published tags
  - Wire release and historical archive recipes into just, Python packaging, and release documentation
- Automate published performance report archiving [`d31e26a`](https://github.com/acgetchell/la-stack/commit/d31e26a9d7a47a6c3089028630640bcff5afe7c0)
  - Track the latest curated release comparison in docs/PERFORMANCE.md and archive older comparisons under docs/archive/performance/
  - Let performance-archive-published discover the latest stable GitHub release and previous stable baseline automatically
  - Generate release comparisons in isolated temporary worktrees, with release-asset restore and local baseline fallback paths
  - Update benchmark and release docs to use the scripted workflow instead of manual checkout steps
- Split local and release performance comparisons [`7258525`](https://github.com/acgetchell/la-stack/commit/7258525590f2ed68d41879e71c833010e408e7f7)
  - Add default performance-local and performance-release workflows that infer the relevant release tags and run in temporary worktrees.
  - Add a performance-github-assets workflow for comparing stored GitHub Actions release benchmark assets without local cargo runs.
  - Normalize release tags before fetching, downloading assets, or checking out detached worktrees.
  - Update performance docs, release guidance, and generated report instructions to use the new benchmark workflows.
- Add vs_linalg-only performance checks [`d7c1487`](https://github.com/acgetchell/la-stack/commit/d7c1487115e1a8e5bb1ec4fcc7592786e300e2ce)
  - Add local workflows for comparing current non-exact la-stack kernels against a release baseline without rerunning current nalgebra/faer or exact benchmarks.
  - Route archive-performance baseline and current benchmark commands by suite, with legacy fallback support for older release worktrees.
  - Document the faster release-signal workflow and expand Semgrep fixtures for benchmark, example, doctest, and public panic-path rules.

### Changed

- Cover determinant scale overflow boundaries [`532093a`](https://github.com/acgetchell/la-stack/commit/532093a1ed9f65ace159f80aac290907d570ea8a)

  - Extract determinant scale exponent calculation into a private helper
  - Assert typed DeterminantScaleOverflow errors for dimension conversion and exponent product overflow
- Harden support script parsing [`87e1d00`](https://github.com/acgetchell/la-stack/commit/87e1d0042ad2de7888c1a065ab78524a63f4c045)
  - Require Python 3.13 for support-script tooling and align Ruff/Ty with that baseline.
  - Replace mypy with strict Ty checking in the Python workflow.
  - Parse TOML, JSON, argparse, and Semgrep inputs into typed boundary objects before downstream use.
  - Reject malformed Criterion estimates, non-finite timings, invalid confidence intervals, and malformed Semgrep result shapes.
- Harden Rust release hygiene [`8e12c93`](https://github.com/acgetchell/la-stack/commit/8e12c935fe54e265e8ceb640702267ec0e71b7b1)
  - Promote missing documentation and dead code lints to deny-level checks.
  - Forbid unsafe code explicitly across Rust modules and benchmark targets.
  - Document the LU/LDLT empty-matrix convention for D=0.
  - Move exact benchmark input generation into typed helpers and consolidate exact benchmark operation dispatch.

### Documentation

- Sync citation metadata for v0.4.2 [`f473ec5`](https://github.com/acgetchell/la-stack/commit/f473ec50946e5f668e9ad9a2d978e499dcb10f04)

  - Update CITATION.cff with the v0.4.2 version and release date.
  - Align the Python utility package metadata and lockfile with the crate release.
  - Add citation metadata validation to the release checklist and config lint flow.
  - Include CITATION.cff in YAML/CFF formatting checks.
- Make release documentation links portable [`dd4ad19`](https://github.com/acgetchell/la-stack/commit/dd4ad192a42e28d9cc72b336b40802fe65cea4f0)

### Fixed

- Escape path regex in benchmark parser test [`1222c93`](https://github.com/acgetchell/la-stack/commit/1222c9325c7b7ea0e1a174e4ee6f0f08e1f7e94b)

  Use a literal regex pattern for the malformed Criterion JSON diagnostic so
  Windows paths with backslashes do not break pytest's match expression.
- Align ty with Python 3.13 [`b9e0ba0`](https://github.com/acgetchell/la-stack/commit/b9e0ba08e54a15d8eddd5c5c53edc37bbc03939a)
- Preserve coordinates for overflowed accumulators [`1d976b3`](https://github.com/acgetchell/la-stack/commit/1d976b346172ad4eca37c68a3ec31817eeca8529)

  - Return matrix-cell metadata when inf-norm row sums or symmetry tolerance scaling overflow.
  - Avoid reparsing finite-by-construction RHS vectors in LU and LDLT solves.
- Re-raise unexpected archive failures [`7938386`](https://github.com/acgetchell/la-stack/commit/7938386166f1f3f5cf594c5def67458d48e19a98)
  - Limit archive-performance CLI error handling to expected validation, filesystem, subprocess, and runtime failures.
  - Let unexpected exceptions propagate so benchmark archiving bugs surface during development.

### Performance

- Improve factorization kernel [`8837df1`](https://github.com/acgetchell/la-stack/commit/8837df1f54a9fa2c20abc1487cfce4de8c8e09c5)

  - Preserve the tiny-dimension update shape for D2-D5 to avoid regressing the core fixed-size path
  - Fuse multiplier computation with trailing updates for larger dimensions to reduce extra column walks
  - Rely on the LDLT factorization proof instead of a redundant final finite-storage scan
- Optimize exact and factorized solve kernels [`1690355`](https://github.com/acgetchell/la-stack/commit/1690355bf27c2cbba685ba0cd70486275c7620b8)
  - Split LU and LDLT solve paths so tiny matrices keep the direct kernels while larger fixed dimensions avoid extra substitution work.
  - Convert dyadic exact solve results directly to finite f64 and preserve UnrepresentableReason recovery semantics on strict conversion failures.
  - Modernize release branch commands and keep just recipes sorted.

## [0.4.2] - 2026-06-04

### Added

- Feat!(matrix): enforce fallible matrix invariants [`e26c283`](https://github.com/acgetchell/la-stack/commit/e26c28358b2358100353b2895441b68892e92cd7)
- Feat!(api): enforce fallible numeric invariants [`adfc33b`](https://github.com/acgetchell/la-stack/commit/adfc33b945b259721bd1067e797ed2e7d4ec0e6e)
- Feat!(matrix): make determinant API tolerance-free [`11a355c`](https://github.com/acgetchell/la-stack/commit/11a355c099eaf366daec8c95af61b6934f914960)
- Feat!(api): hide finite and symmetry proofs behind matrix APIs
  [`7219336`](https://github.com/acgetchell/la-stack/commit/721933671c28eb71953f1386a201622d6171caf7)
- Guard public Rust examples against unwrap [`df1130a`](https://github.com/acgetchell/la-stack/commit/df1130a7ad0ba69a1072ef231e14f3efb7e4b8de)

  - Add repository-owned Semgrep rules for unwrap and expect usage in public doctests, examples, and benchmarks.
  - Add fixture-based Semgrep rule tests and include them in the lint workflow.
  - Update examples and benchmarks to model typed fallible flow or operation-labeled benchmark failures.
- Feat!(api): enforce finite Matrix and Vector construction [`92ba403`](https://github.com/acgetchell/la-stack/commit/92ba4034b194875c62a27f24dfbf6d43f380f54e)

### Changed

- Remove redundant cache restore keys for cargo-llvm-cov [`f75a01c`](https://github.com/acgetchell/la-stack/commit/f75a01c99c8dbcc8b6ffc36ae9f94ba968a2f111)

  Remove the broad restore-keys pattern to ensure only exact version matches
  are used from the cache, preventing potential version mismatches during
  CI runs.
- Revert "ci: modernize tooling checks and example execution"
  [`a4b9f64`](https://github.com/acgetchell/la-stack/commit/a4b9f64235f53274270ef9272634a61f653dad87)
- Reapply "ci: modernize tooling checks and example execution"
  [`758321a`](https://github.com/acgetchell/la-stack/commit/758321acf872b1f17286ff3bb7bee6a807e4b440)
- Encode nonzero mantissas in exact decomposition [`7a664ed`](https://github.com/acgetchell/la-stack/commit/7a664ede2f4add168c5813f8d24e16732fa03b30)

  - Replace the exact-arithmetic zero mantissa sentinel with `Option&lt;NonZeroU64&gt;`.
  - Carry nonzero mantissa proof through matrix/vector decomposition and BigInt scaling.
  - Clarify determinant documentation around uncertified `det()` bounds.
  - Keep SPD determinant proptests on the tolerance-aware LU path.
- Simplify finite proof wrappers [`54b603c`](https://github.com/acgetchell/la-stack/commit/54b603c21f0eb5b4d63ff334fac7f8cc325ebc2e)
  - Use the checked proof-wrapper constructors as the single internal path for finite matrices and vectors.
  - Remove exact-arithmetic tests that duplicated the matrix and vector non-finite boundary checks.

### Dependencies

- Bump taiki-e/install-action from 2.75.18 to 2.75.22 [`d6c944b`](https://github.com/acgetchell/la-stack/commit/d6c944bb7dd30bb00dfe820bc355c4351cb1f242)
- Bump actions/setup-node from 6.3.0 to 6.4.0 [`1fc57e1`](https://github.com/acgetchell/la-stack/commit/1fc57e1d5a2e9aed2cc67c90a6ce63b9eb469e2e)
- Bump pastey from 0.2.1 to 0.2.2 in the dependencies group [`5b215d5`](https://github.com/acgetchell/la-stack/commit/5b215d5a87a25744d313867471943a614e81abd0)
- Bump taiki-e/install-action from 2.75.22 to 2.77.3 [`3e5ce42`](https://github.com/acgetchell/la-stack/commit/3e5ce4273cbe02e0792747b8170efde5d079ae0c)
- Bump actions-rust-lang/setup-rust-toolchain [`8532b08`](https://github.com/acgetchell/la-stack/commit/8532b084168e07cfe6b48bb90fa2e533cb69d990)
- Bump taiki-e/install-action from 2.77.3 to 2.78.3 [`bef2caf`](https://github.com/acgetchell/la-stack/commit/bef2cafcedd072857742da53256dd472a5855655)
- Bump the dependencies group across 1 directory with 2 updates
  [`ba01a14`](https://github.com/acgetchell/la-stack/commit/ba01a14346b84fd6385ee8a6acff8b2180080df5)
- Bump taiki-e/install-action from 2.78.3 to 2.80.0 [`480e61b`](https://github.com/acgetchell/la-stack/commit/480e61b85638fb044edccd8116eec777f5588f7d)
- Bump codecov/codecov-action from 6.0.0 to 6.0.1 [`3ca8eed`](https://github.com/acgetchell/la-stack/commit/3ca8eed89997379e6cbad863ecc274bba1651749)

### Documentation

- Document feature requirement for exact APIs [`19b10d5`](https://github.com/acgetchell/la-stack/commit/19b10d552e83b6a7f9e91695b4850b8fab3f4550)
- Document scalar scope and release roadmap [`bfb0393`](https://github.com/acgetchell/la-stack/commit/bfb039386588f94b95561c610181ca6d486acd6e)

  - Clarify that la-stack intentionally supports f64 floating-point APIs plus optional exact rationals, not alternate scalar families.
  - Add a roadmap covering the v0.4.x stable-Rust issue sequence and the v0.5.0 generic_const_exprs anchor.
  - Refresh generated changelog entries and archived changelog grouping.
- Document finite RHS solve validation [`075aed7`](https://github.com/acgetchell/la-stack/commit/075aed78cf8264fc920258f1f1d977ddd589ffd7)
  - Document that LU and LDLT solve_vec reject non-finite RHS entries with LaError::NonFinite metadata.
  - Cite the Bareiss reference in the exact solve helper docs and describe exact-arithmetic growth and complexity.
  - Cover finite proof defaults and non-finite RHS solve boundaries in unit tests.
- Clarify finite solve and norm guarantees [`fb71485`](https://github.com/acgetchell/la-stack/commit/fb71485cac0b464b2fa9ee949140a558b7738781)
  - State that LU and LDLT solve_vec use floating-point substitution without a certified absolute rounding-error bound.
  - Clarify that inf_norm reports NonFinite for unchecked stored NaN/∞ as well as row-sum overflow.
  - Exercise the unchecked finite-proof fixture path directly in exact tests.
- Update v0.4.2 release notes [`7e11f93`](https://github.com/acgetchell/la-stack/commit/7e11f930b94bbba99c1c426e68a515bbefb8c489)

### Fixed

- Reject overflowed symmetry tolerance scaling [`a7b052a`](https://github.com/acgetchell/la-stack/commit/a7b052af5dc6361198bbfe1e17d6b1f0ba225ed7)

  - Enforce the tolerance contract around symmetry checks by surfacing scaled
    tolerance overflow as a typed non-finite intermediate error.
  - Document finite, non-negative tolerance requirements across tolerance-taking
    matrix APIs.
  - Add regression coverage for invalid tolerance construction and symmetry
    tolerance overflow.
  - Update exact examples to propagate typed crate errors instead of unwrapping.
- Harden Semgrep fixture parsing [`ac44c07`](https://github.com/acgetchell/la-stack/commit/ac44c078cc4435d5beca27f1890fbb4046cf5952)
  - Ignore non-canonical todoruleid annotations when counting expected rule hits.
  - Reject malformed Semgrep JSON results with clear stderr diagnostics instead of propagating KeyError.
- Revalidate finite proof conversions [`419a90f`](https://github.com/acgetchell/la-stack/commit/419a90f7267608051736498154ac5e6faf0909c5)

  Ensure internal finite proof conversions cannot accept raw Matrix or Vector storage without checking the invariant.
  - Revalidate TryFrom&lt;Matrix&lt;D&gt;&gt; and TryFrom&lt;Vector&lt;D&gt;&gt; before constructing finite wrappers.
  - Measure exact random percentile benchmarks over repeated corpus timings and cumulative input sets.
  - Tighten Codecov status thresholds and extend benchmark workflow timeout.
  - Keep Semgrep constructor fixtures aligned with public API guardrails.
- Revalidate public compute inputs [`ffca00e`](https://github.com/acgetchell/la-stack/commit/ffca00e9dd6fde5e57c8064f69807a76e45a469e)
  - Parse Matrix and Vector storage into private finite proof-bearing types at public compute boundaries.
  - Reject unchecked non-finite storage before LU, LDLT, determinant, norm, dot, and exact-arithmetic paths can proceed.
  - Keep unchecked proof-wrapper constructors crate-private for internal paths with local finiteness proofs.
  - Document the private proof-bearing invariant model in the API overview and roadmap.

### Maintenance

- Migrate coverage to cargo-llvm-cov [`66f2117`](https://github.com/acgetchell/la-stack/commit/66f21173344033b3003dd6e46b8400b08e6fdd40)
- Align shared security and release tooling [`b303509`](https://github.com/acgetchell/la-stack/commit/b3035095964d9491f0a3d5c4fb790ab5c22cd06f)

  - Run CI through pinned uv and cached Cargo tooling across all supported platforms.
  - Add repository Semgrep rules, action allowlist coverage, CodeQL, zizmor, and SARIF workflows.
  - Raise the Rust and Python tooling floors to Rust 1.96 and Python 3.12.
  - Add archive-aware changelog generation, post-processing, and tag-release tooling.
  - Preserve exact arithmetic overflow reporting without non-finite sentinel defaults.
- Modernize tooling checks and example execution [`da626bc`](https://github.com/acgetchell/la-stack/commit/da626bcca899aa91d58f728db433a53a46179e92)
  - Run examples from prebuilt binaries instead of invoking cargo run for each example.
  - Add check/fix recipe aliases and guard documented command ordering with Semgrep.
  - Document the Rust-native Markdown, YAML, TOML, spelling, and workflow action policy.
  - Tighten Python/tooling setup by pinning pytest, using .python-version in CI, and hardening setup-tool version checks.
  - Keep CI running the full just ci workflow across Ubuntu, macOS, and Windows.

### Performance

- Restore checked factorization throughput [`f2e6d56`](https://github.com/acgetchell/la-stack/commit/f2e6d560bd1cb327ff0616644ef704a865520e32)

  - Move LU and LDLT non-finite factor checks out of cubic update loops while preserving completed-storage validation before factors escape.
  - Borrow factor rows and finite vector arrays in solve, dot, and norm paths to avoid avoidable copies.
  - Refresh v0.4.2 release artifacts, LU solve benchmark docs, and release workflow docs.
  - Preserve README benchmark table spacing in the Criterion plot updater.

## [0.4.1] - 2026-04-21

### Merged Pull Requests

- Propagate NaN in Matrix::inf_norm [#85](https://github.com/acgetchell/la-stack/pull/85)
- Add adversarial-input coverage for exact arithmetic [#80](https://github.com/acgetchell/la-stack/pull/80)
- Integer-only forward elimination for gauss_solve [#72](https://github.com/acgetchell/la-stack/pull/72)

### Added

- Regression tests for solver and determinant overflow handling
  [`f763b11`](https://github.com/acgetchell/la-stack/commit/f763b119bcc57276b83f370b0bf7abce654c7eb8)

  Add test cases to verify that LDLT and LU solvers, as well as determinant
  calculations, correctly detect and return `LaError::NonFinite` when
  intermediate calculations overflow to infinity despite having finite
  inputs.
- Defensive-path test coverage for LU and LDLT solve_vec [`87d426f`](https://github.com/acgetchell/la-stack/commit/87d426fca1627445b804fd26b62fc7d9d4f0ae48)

  Add unit tests to exercise internal safety nets in the LU and LDLT
  diagonal solve routines. These tests manually construct factors with
  invalid diagonals (NaN or sub-tolerance) to verify that solve_vec
  correctly surfaces NonFinite and Singular errors, even though these
  states are unreachable via the standard factorization APIs.
- Const-ify Lu/Ldlt det + solve_vec and Matrix inf_norm + det_errbound
  [`81ecb35`](https://github.com/acgetchell/la-stack/commit/81ecb35bdaf159f1f44d1eb24274ecf82c6567d5)

  Make six public methods `const fn`, completing const-evaluation parity
  with `Matrix::det_direct` and `Vector::dot` now that MSRV 1.95 exposes
  `f64::mul_add` as const (1.94) and `core::hint::cold_path` as const
  (1.95):

  - `Lu::det`, `Lu::solve_vec`
  - `Ldlt::det`, `Ldlt::solve_vec`
  - `Matrix::inf_norm`, `Matrix::det_errbound`

  Iterator chains (`.iter().map().sum()`, `.enumerate().take(i)`,
  `.enumerate().skip(i + 1)`, `.iter_mut().enumerate().take(D)`) were
  rewritten as `while` loops since they are not const-stable.

  Fix error-variant correctness in both solve_vec paths:
  - A corrupt stored `U` / `D` diagonal at `(i, i)` now surfaces as
    `LaError::NonFinite { row: Some(i), col: i }`, matching the
    convention used by `Matrix::det`, `Lu::factor`, and `Ldlt::factor`.
  - A computed-intermediate overflow keeps `row: None, col: i`.
  - Previously both were conflated into `row: None, col: i`, defeating
    debuggability for callers who construct factorizations directly.

  Sharpen `LaError::NonFinite` variant docs and the `# Errors` sections of
  `Lu::solve_vec` / `Ldlt::solve_vec` to spell out the `(row, col)`
  contract for each failure mode.
- Fast-filter boundary proptests for exact determinant sign [`6357db3`](https://github.com/acgetchell/la-stack/commit/6357db35c70bca1b93e6bbf9a4fd231913631950)

  Introduce proptests to verify that the floating-point determinant sign
  matches the exact sign whenever the error bound is satisfied. This
  validates the correctness of the Shewchuk-style filter used in determinant
  calculations.

  Additionally, improve test reliability by correcting the diagonal dominance
  calculation to properly account for the maximum off-diagonal magnitude and
  refactor matrix-vector multiplication helpers to use idiomatic iterators.

### Changed

- Finalize documentation, benchmarks, and error handling [`0b98d3f`](https://github.com/acgetchell/la-stack/commit/0b98d3f6dbdd74699c318c4744a2b2f9a1b78481)

  Synchronize the CHANGELOG, README, and REFERENCES for the upcoming release.
  This includes documenting the new absolute error bound coefficients for
  determinants, adding a dedicated section for AI agent governance, and
  refining internal error helpers to use a more consistent coordinate
  convention for non-finite values.
- Consolidate and expand const-evaluability tests via macros [`f8d80a0`](https://github.com/acgetchell/la-stack/commit/f8d80a01e9913f87e1b19b2ad5ffbc0994e2bfdb)

  Refactor manual const-evaluation tests in `Ldlt` and `Matrix` into macros
  to standardize coverage across matrix dimensions 2 through 5. This
  ensures all determinant and norm variants are fully exercised at compile
  time.
- Refactor solve_exact to use hybrid Bareiss forward elimination
  [`ecbbe8a`](https://github.com/acgetchell/la-stack/commit/ecbbe8a571ccaeb9cfedbf0269b8db44d43a5773)

  Migrate the exact linear solver to a hybrid approach that performs O(D³)
  forward elimination using fraction-free Bareiss updates on integers,
  limiting BigRational arithmetic to the O(D²) back-substitution phase.
  This eliminates per-entry GCD overhead during elimination. Internal f64
  decomposition is unified via f64_decompose, and f64_to_bigrational is
  moved to test scope.
- Polish exact module (Component struct, errors, perf) [`53a5be6`](https://github.com/acgetchell/la-stack/commit/53a5be6abecc0af332398236ed6803ed75564b03)

  Major refactor of `src/exact.rs`:

  - Extract shared Bareiss primitives (`decompose_matrix`, `decompose_vec`,
    `component_to_bigint`, `build_bigint_matrix`, `build_bigint_vec`,
    `bareiss_forward_eliminate`) so `bareiss_det_int` and `gauss_solve`
    share a single integer-Bareiss core.
  - Hybrid BigInt/BigRational solve: forward elimination runs entirely in
    `BigInt` with fraction-free Bareiss updates on `(A | b)`; only the
    `O(D²)` back-substitution phase lifts into `BigRational`.
  - `mem::take` in back-substitution eliminates per-entry `clone()` calls
    on `rhs[i]`, `a[i][j]`, `a[i][i]`.

  Structured errors replace preconditions:
  - `decompose_matrix` / `decompose_vec` fold `is_finite()` into the same
    pass that decomposes each entry and return
    `Err(LaError::NonFinite { row, col })` on the first non-finite cell.
  - `bareiss_det_int`, `bareiss_det`, `gauss_solve` now return
    `Result&lt;_, LaError&gt;`; error propagates to every public entry point via
    `?`.
  - `validate_finite` and `validate_finite_vec` removed (dead after the
    refactor); `det_sign_exact` relies on IEEE 754 NaN/∞ propagation
    through `det_direct()` plus `bareiss_det_int` for Stage-2 validation.
- Add adversarial-input coverage for exact arithmetic [#80](https://github.com/acgetchell/la-stack/pull/80)
  [`5bf5815`](https://github.com/acgetchell/la-stack/commit/5bf5815cb165c3b6145c5592420a58085a66efaa)

  Expand benchmark, unit-test, and proptest coverage of the exact-arithmetic
  APIs to catch tail cases that fixed well-conditioned inputs miss.

  Benchmarks (benches/exact.rs):
  - Factor out `bench_extreme_group` helper running the same four benches
    (`det_sign_exact`, `det_exact`, `solve_exact`, `solve_exact_f64`) so
    adversarial groups are directly comparable.
  - Extend `exact_near_singular_3x3` with the two solve benches (the
    primary motivating use case for exact solve was previously unmeasured).
  - Add `exact_large_entries_3x3` (entries near `f64::MAX / 2`) to stress
    intermediate BigInt growth during Bareiss forward elimination.
  - Add `exact_hilbert_4x4` / `exact_hilbert_5x5` to stress the
    `f64_decompose → BigInt` scaling path on ill-conditioned inputs.
  - Tighten bench `expect(...)` messages to name the invariant each call
    relies on (e.g. "non-singular matrix with finite entries") so panics
    identify both where (Criterion bench name) and why.

  Unit tests (src/exact.rs):
  - `solve_exact_near_singular_3x3_integer_x0` — integer-x0 round-trip
    through the 2^-50-perturbed matrix.
  - `solve_exact_large_entries_3x3_unit_vector` — `A·[1,0,0] = [big,1,1]`
    round-trip with `f64::MAX/2` diagonal.
  - `det_sign_exact_large_entries_3x3_positive` — asserts the fast filter
    falls through (`det_direct` is non-finite) and `det_exact_f64` returns
    `Overflow { index: None }`.
  - `det_sign_exact_hilbert_positive_{3,4,5}d` — Hilbert is SPD, sign = 1.
  - `solve_exact_hilbert_residual_{3,4,5}d` — residual `A·x - b` is exactly
    zero in `BigRational`, stronger than integer round-trips since Hilbert
    entries are non-terminating in binary.

  Proptests (tests/proptest_exact.rs):
  - `solve_exact_integer_roundtrip_{2..5}d` — random diagonally-dominant
    integer `A` and small-integer `x0`, verify `solve(A, A·x0) == x0` exactly.
  - `solve_exact_residual_{2..5}d` — random `A` + small-integer `b`, verify
    `A · solve(A, b) == b` exactly (catches back-sub bugs on fractional
    solutions).
  - `det_sign_exact_agrees_with_det_exact_{2..5}d` — on full (non-diagonal)
    small-integer matrices, asserts `det_sign_exact() == det_exact().sign()`
    (exercises the filter/fallback boundary previously only diagonal-tested).

  Prelude (src/lib.rs):
  - Re-export `BigInt` alongside `BigRational` (crate root + prelude).
  - Re-export `FromPrimitive`, `Signed`, `ToPrimitive` from `num-traits` in
    the prelude so the re-exported `BigRational` is actually usable for
    construction (`from_f64`, `from_i64`) and sign queries without forcing
    downstream users to add `num-bigint` / `num-rational` / `num-traits` to
    their own Cargo.toml. Additive; no public API breaks.

  Tooling (scripts/bench_compare.py + test_bench_compare.py):
  - Register the three new adversarial groups in `EXACT_GROUPS`.
  - Extend `_group_heading` with human-readable titles
    ("Large entries 3x3", "Hilbert 4x4", "Hilbert 5x5").
  - Add group-heading unit tests for the new cases.

  Test results (`just ci`):
  - cargo test --features exact: 368 lib + 20 proptest_exact + 40 other
    proptests + 34 doc-tests — all pass
  - cargo test (no features): 175 lib + 40 proptests + 29 doc-tests — all pass
  - Python: 104 tests pass (ty, mypy, ruff clean)
  - Clippy (pedantic + nursery + cargo, `-D warnings`): clean
  - fmt, taplo, yamllint, shellcheck, spell-check, bench-compile, examples: clean
- Expand exact-arithmetic re-exports and adversarial benchmarks
  [`b1a491d`](https://github.com/acgetchell/la-stack/commit/b1a491d902ccdaba6f9cd2e6f8e05514b6dfa3de)

  Update the `exact` feature to re-export `BigInt` and essential `num-traits`
  (FromPrimitive, ToPrimitive, Signed) alongside `BigRational`. This
  simplifies downstream usage by removing the need for direct dependencies
  on the `num` ecosystem for common operations.

  Additionally, formalize adversarial-input benchmarking for near-singular,
  large-entry, and Hilbert matrices to better track performance on
  ill-conditioned data. CI benchmark comparisons are now more resilient to
  newly added benchmarks via `--baseline-lenient`.
- Update AGENTS.md [`1e0648d`](https://github.com/acgetchell/la-stack/commit/1e0648dad3147ef127c775d8969b7cd214a2a6ed)

### Dependencies

- Bump taiki-e/install-action from 2.73.0 to 2.75.7 [`a1d1c1e`](https://github.com/acgetchell/la-stack/commit/a1d1c1edba7d5bf6928de4e8cab86ba853685430)
- Bump actions/download-artifact from 4.3.0 to 8.0.1 [`8772ea2`](https://github.com/acgetchell/la-stack/commit/8772ea2f18cf1b94a21e269b4bc0c8e0a3b650eb)
- Bump rand from 0.9.2 to 0.9.4 [`b91670a`](https://github.com/acgetchell/la-stack/commit/b91670a46f82306dde5433ca7c406a8757e3fc59)
- Bump actions/upload-artifact from 7.0.0 to 7.0.1 [`7579274`](https://github.com/acgetchell/la-stack/commit/7579274d56ab8c46eae6906c07a3c0d0fd41c84f)
- Bump actions/github-script from 8.0.0 to 9.0.0 [`a522abd`](https://github.com/acgetchell/la-stack/commit/a522abdabf211d47f0f52cbfc2a349fa72ce87e2)
- Bump actions/cache from 5.0.4 to 5.0.5 [`b12d95c`](https://github.com/acgetchell/la-stack/commit/b12d95c8ce242ea6014c6a6aaa7845ff501d2e4d)
- Bump actions-rust-lang/setup-rust-toolchain [`c728774`](https://github.com/acgetchell/la-stack/commit/c72877460efba3d598f64ff253d2df8f3695acb1)
- Bump astral-sh/setup-uv from 8.0.0 to 8.1.0 [`7a18733`](https://github.com/acgetchell/la-stack/commit/7a18733ba0f39102678bbe30772fcde10f20e9d8)
- Bump taiki-e/install-action from 2.75.7 to 2.75.18 [`433cfc1`](https://github.com/acgetchell/la-stack/commit/433cfc1b01c9128e93c82cb553aa63d4091bace3)

### Documentation

- Strengthen LDLT symmetry precondition and add is_symmetric API
  [`1693307`](https://github.com/acgetchell/la-stack/commit/1693307c58e30a804e9b79abc783af99be8dc947)

  LDLT factorization assumes the input matrix is symmetric, but the contract
  was only implied by one sentence in the module and struct docs. Asymmetric
  inputs silently produce mathematically meaningless factorizations in
  release builds (the `debug_assert_symmetric` check is compiled out), and
  callers had no supported way to validate symmetry up front.

### Fixed

- Propagate NaN in Matrix::inf_norm [#85](https://github.com/acgetchell/la-stack/pull/85)
  [`16ffa45`](https://github.com/acgetchell/la-stack/commit/16ffa45ade11b21a179cad8fcecc51d802086a1d)
- Report infinite vs finite off-diagonal pairs as asymmetric [`1805779`](https://github.com/acgetchell/la-stack/commit/1805779dbca49183fbfa95c68ec00984966aa551)

  Update Matrix::first_asymmetry to flag any non-finite difference between
  off-diagonal pairs as an asymmetry. This prevents cases where a single
  infinite entry paired with a finite entry would incorrectly pass as
  symmetric because the matrix's infinite norm blew the tolerance up to
  infinity, making the comparison `diff &gt; eps` return false.

### Maintenance

- Bump MSRV to Rust 1.95 and adopt new stable features [`0ab3c33`](https://github.com/acgetchell/la-stack/commit/0ab3c336074f2b866256fbe5db8a8ec5306d580a)

  - Bump rust-version to 1.95 in Cargo.toml
  - Bump channel to 1.95.0 in rust-toolchain.toml
  - Add core::hint::cold_path() hints at cold/error branches:
    - src/exact.rs: validate_finite, validate_finite_vec, gauss_solve
      singular return, det_exact_f64 / solve_exact_f64 overflow returns,
      det_sign_exact Stage 2 Bareiss fallback
    - src/lu.rs: Lu::factor and Lu::solve_vec NonFinite / Singular returns
    - src/ldlt.rs: Ldlt::factor and Ldlt::solve_vec NonFinite / Singular
      returns
    - src/matrix.rs: det_direct D &gt;= 5 fallback arm (legal because
      cold_path is const fn in 1.95) and det NonFinite / overflow scan
  - Refactor det_sign_exact Stage 1 fast filter to use match + if let
    guard with let-chain, replacing the tuple destructure; semantics
    unchanged

  Test results (local `just ci`):
  - cargo fmt --all -- --check: clean
  - cargo clippy --workspace --all-targets --all-features
    -D warnings -W clippy::pedantic -W clippy::nursery -W clippy::cargo:
    clean
  - cargo test --features exact --lib: 258 passed, 0 failed
  - cargo test --features exact (doctests + examples): 31 passed, 0 failed
  - RUSTDOCFLAGS='-D warnings' cargo doc --no-deps --features exact: clean
  - python tests: 101 passed

### Performance

- Integer-only forward elimination for gauss_solve [#72](https://github.com/acgetchell/la-stack/pull/72)
  [`a1d3bdb`](https://github.com/acgetchell/la-stack/commit/a1d3bdbdb6fcb778a78a2a3d0cb66b79484e1472)

  Replace BigRational-only gauss_solve with a hybrid that runs
  Bareiss fraction-free forward elimination in BigInt on the
  augmented (A | b) system, then back-substitutes in BigRational.
  Eliminates GCD normalization from the O(D^3) phase while keeping
  rational overhead limited to the cheaper O(D^2) back-sub.

  Scope f64_to_bigrational to cfg(test); production code paths now
  use f64_decompose directly (shared with bareiss_det_int).

  Closes #72

  Co-Authored-By: Oz <oz-agent@warp.dev>

## [0.4.0] - 2026-04-11

### Merged Pull Requests

- Integer-only Bareiss determinant via BigInt [#64](https://github.com/acgetchell/la-stack/pull/64)
- Custom f64 → BigRational via IEEE 754 bit decomposition [#63](https://github.com/acgetchell/la-stack/pull/63)

### Added

- Add benchmark comparison tool and performance documentation
  [`d448fc3`](https://github.com/acgetchell/la-stack/commit/d448fc3e3c884d6145580e5a40eb265bf1f83cb8)

  Introduce a Python script to compare Criterion results across baselines
  and generate markdown reports. Integrate performance snapshots into the
  release workflow and add justfile recipes for baseline management to
  track exact arithmetic performance.

### Changed

- Update version to 0.3.0 and refine metadata keywords [`cc0e967`](https://github.com/acgetchell/la-stack/commit/cc0e967aedb4eabbcd98714b16ef67a8290c112e)

  Update CITATION.cff to reflect the v0.3.0 release. Keywords are updated to
  include exact arithmetic and robust predicates, with existing terms
  standardized to kebab-case.
- Refine benchmark comparison reporting and documentation [`e1b5955`](https://github.com/acgetchell/la-stack/commit/e1b5955fb5024232e34e9df9701bb24fb98efa15)

  Refactor the `bench_compare.py` script to group results by dimension and
  add a new test suite. Consolidate all benchmarking workflows into a new
  dedicated guide and move machine-specific performance snapshots to
  `.gitignore`.
- Expand test coverage for benchmark comparison edge cases [`bced7d9`](https://github.com/acgetchell/la-stack/commit/bced7d988bd2f42cb6bb5af9c54dabdcf787a5fc)
- Update documentation and tests for integer-only Bareiss [`2ee3f05`](https://github.com/acgetchell/la-stack/commit/2ee3f05caecfdf1a23b61257a7465b3bb6d63614)

  Update architecture descriptions in AGENTS.md and REFERENCES.md to
  reflect the move from BigRational to integer-only BigInt arithmetic.
  Add unit tests for bareiss_det_int covering negative inputs and pivot
  swapping. Refine gh CLI patterns for automated issue listing.
- Restrict benchmark baselines to main and improve reporting [`9a7caa2`](https://github.com/acgetchell/la-stack/commit/9a7caa241b5f476c2659772f3189468b10fcba2e)

  Update internal CI logic to ensure Criterion baselines are only saved
  and uploaded for the main branch. Add error handling to regression
  reporting to provide clear feedback when comparison data is missing.

### Dependencies

- Bump astral-sh/setup-uv from 7.3.1 to 7.5.0 [`0505061`](https://github.com/acgetchell/la-stack/commit/050506123b1da9bc79b04e7d8b44d8195d115199)
- Bump taiki-e/install-action from 2.68.25 to 2.68.31 [`dcff501`](https://github.com/acgetchell/la-stack/commit/dcff501354df2e6434fd6f817a7c15b86a40685f)
- Bump actions-rust-lang/setup-rust-toolchain [`4678d2e`](https://github.com/acgetchell/la-stack/commit/4678d2e8b121c3116d638b59b03eca95510888fe)
- Bump actions/cache from 5.0.3 to 5.0.4 [`633c6e2`](https://github.com/acgetchell/la-stack/commit/633c6e234dc4a94744fa433f41bceb84efe1c99c)
- Bump taiki-e/install-action from 2.68.34 to 2.69.6 [`89a306b`](https://github.com/acgetchell/la-stack/commit/89a306b64bf586ec201957953afac0a3d9c13eb7)
- Bump astral-sh/setup-uv from 7.5.0 to 7.6.0 [`676a6ff`](https://github.com/acgetchell/la-stack/commit/676a6ffe3df416ab316daa84119161893eeeadd0)
- Bump codecov/codecov-action from 5.5.2 to 5.5.3 [`2c15dc1`](https://github.com/acgetchell/la-stack/commit/2c15dc1a0c93f75dd23e0d540041bc01f2e3b71a)
- Bump the dependencies group with 2 updates [`23ade13`](https://github.com/acgetchell/la-stack/commit/23ade1343e575edead08ab49a4e361c9e37f5ca0)
- Bump taiki-e/install-action from 2.69.6 to 2.70.1 [`4073197`](https://github.com/acgetchell/la-stack/commit/40731976f465f69bacb0a6dc4b8520fb9bc9c8f4)
- Bump codecov/codecov-action from 5.5.3 to 6.0.0 [`b3e1380`](https://github.com/acgetchell/la-stack/commit/b3e1380e0b8df85478648036eeb1d9b2c79aaac5)
- Bump taiki-e/install-action from 2.70.1 to 2.73.0 [`7de720d`](https://github.com/acgetchell/la-stack/commit/7de720dfb8328d01843cfabd15d23086ee98832b)
- Bump astral-sh/setup-uv from 7.6.0 to 8.0.0 [`af40753`](https://github.com/acgetchell/la-stack/commit/af40753130fac56d48da7fce2f18b11dc391ebe6)

### Fixed

- Improve non-finite error reporting and coordinate accuracy [`d4b1452`](https://github.com/acgetchell/la-stack/commit/d4b14524a02b49a1a8e832818913ab875c4e5fe0)

  Refine error reporting by ensuring correct column indices in LDLT
  decomposition and LU back-substitution. Add missing finiteness checks to
  LU solve paths and update determinant calculations to provide specific
  row/column coordinates for non-finite inputs.

### Maintenance

- Add performance regression detection for exact-arithmetic benchmarks
  [`44bce99`](https://github.com/acgetchell/la-stack/commit/44bce99bcae8f852dbf7500e71fa182730e08bca)

  Add a GitHub Actions workflow that detects performance regressions in
  the exact-arithmetic benchmark suite on PRs and pushes to main.

  - On push to main: run exact benchmarks, save Criterion baseline as
    artifact (30-day retention)
  - On PRs: find latest main baseline artifact, download, run benchmarks
    with --baseline comparison, report regressions in job summary
  - Regressions are warning-only — the workflow never blocks PRs
  - Uses artifact-based baselines (not cache) for robustness: explicit
    provenance, no silent eviction, fork-safe
  - Validated locally: --save-baseline and --baseline flags confirmed
    working with Criterion

### Performance

- Stack-backed storage for exact arithmetic kernels [`b64ac60`](https://github.com/acgetchell/la-stack/commit/b64ac60e286bcdc1c21334ab271e6ffacfadc187)

  Replace heap-allocated Vec&lt;Vec&lt;BigRational&gt;&gt; with stack-backed
  [[BigRational; D]; D] arrays in bareiss_det and gauss_solve, and return
  [BigRational; D] from gauss_solve directly — eliminating per-call heap
  traffic and a redundant clone in solve_exact.

  Exact kernel changes (src/exact.rs):

  - bareiss_det: Vec&lt;Vec&lt;BigRational&gt;&gt; → [[BigRational; D]; D] via
    std::array::from_fn
  - gauss_solve: Vec&lt;Vec&lt;BigRational&gt;&gt; augmented matrix → separate
    [[BigRational; D]; D] + [BigRational; D] (avoids unstable
    const-generic [T; D+1])
  - gauss_solve return type: Vec&lt;BigRational&gt; → [BigRational; D]
  - solve_exact: drop intermediate clone, forward array directly

  Error enrichment (src/lib.rs, src/lu.rs, src/ldlt.rs, src/matrix.rs,
  src/exact.rs):
  - LaError::NonFinite gains row: Option&lt;usize&gt; for full (row, col)
    location on matrix inputs; None for vectors/intermediates
  - LaError::Overflow becomes a struct variant with index: Option&lt;usize&gt;
    to identify which solution component overflowed
  - Display messages updated to show positional context
- Custom f64 → BigRational via IEEE 754 bit decomposition [#63](https://github.com/acgetchell/la-stack/pull/63)
  [`0a8ce5b`](https://github.com/acgetchell/la-stack/commit/0a8ce5b45e21fa88e4e9832bc6ada38b3ba68eeb)

  Replace `BigRational::from_float(x)` in `f64_to_bigrational` with manual
  IEEE 754 binary64 bit decomposition and `BigRational::new_raw`, bypassing
  the unnecessary GCD normalization that `from_float` performs internally.
  - Decompose f64 into sign, biased exponent, and significand fields
  - Strip trailing zeros from the significand so the fraction is already
    in lowest terms (odd numerator over power-of-two denominator)
  - Handle zero (±0.0), subnormals, and normal values; panic on NaN/Inf
  - Add 15 dedicated unit tests: ±zero, integers, fractions, powers of
    two, subnormals, round-trip fidelity, lowest-terms, and panic cases
  - Add IEEE 754-2019 [9] and Goldberg [10] references to REFERENCES.md
  - Add module-level and function-level doc citations
- Integer-only Bareiss determinant via BigInt [#64](https://github.com/acgetchell/la-stack/pull/64)
  [`d422b25`](https://github.com/acgetchell/la-stack/commit/d422b251ca86a914522f80285964d4513bca1817)

  Replace the BigRational Bareiss determinant with a pure BigInt path:
  all f64 entries are decomposed into mantissa × 2^exponent, scaled to
  a common integer base, and eliminated without any rational arithmetic.
  The result is reconstructed as BigRational only at the end.
  - Add f64_decompose helper (extracted from f64_to_bigrational)
  - Add bareiss_det_int: integer-only Bareiss returning (BigInt, i32)
  - Add bigint_exp_to_bigrational: reconstruction with trailing-zero
    reduction (no full GCD needed)
  - Refactor bareiss_det as thin wrapper over bareiss_det_int
  - Optimize det_sign_exact to read sign from BigInt directly
    (skip BigRational reconstruction entirely)
  - Import std::array::from_fn for shorter call sites
  - Replace clippy cast suppression with exact try_from conversions
  - Update AGENTS.md with non-interactive gh CLI patterns
  - 24 new tests (256 total): f64_decompose, bareiss_det_int (D=0–5,
    fractional, all-zeros, sign), bigint_exp_to_bigrational (positive/
    negative exp, reduction, negative values)

  Performance (vs pre-bigint baseline on Apple M4 Max):
  det_exact: 16x (D=2) → 39x (D=5) faster
  det_exact_f64: 10x (D=2) → 38x (D=5) faster
  det_sign_exact: 40x at D=5 (bypasses BigRational entirely)
  Near-singular: 18x faster

## Archives

Older releases are archived by minor series:

- [0.3.x](docs/archive/changelog/0.3.md)
- [0.2.x](docs/archive/changelog/0.2.md)
- [0.1.x](docs/archive/changelog/0.1.md)

[0.4.4]: https://github.com/acgetchell/la-stack/compare/v0.4.3...v0.4.4
[0.4.3]: https://github.com/acgetchell/la-stack/compare/v0.4.2...v0.4.3
[0.4.2]: https://github.com/acgetchell/la-stack/compare/v0.4.1...v0.4.2
[0.4.1]: https://github.com/acgetchell/la-stack/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/acgetchell/la-stack/compare/v0.3.0...v0.4.0
