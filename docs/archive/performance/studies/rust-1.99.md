# Rust 1.99 Compiler Migration

This study owns the strict floating-point baseline for
[#251](https://github.com/acgetchell/la-stack/issues/251) and the subsequent
[#250](https://github.com/acgetchell/la-stack/issues/250) evaluation. It compares
compilers on the same library and benchmark sources, separately from the
generated [release performance report](../../../performance.md).

## Contents

- [Scope and compatibility decisions](#scope-and-compatibility-decisions)
- [Correctness and portability](#correctness-and-portability)
- [Measurement provenance](#measurement-provenance)
- [Measurements](#measurements)
- [Code generation](#code-generation)
- [Decision and future baseline](#decision-and-future-baseline)

## Scope and compatibility decisions

Reviewed the final [Rust 1.99 release notes][release], the
[Cargo changelog][cargo], the [LLVM 23 integration][llvm], and the
[inclusive-range optimization][range]. Repository-specific conclusions:

- **Cargo:** keep the existing release profile, resolver, and feature declarations.
  The new `debug` profile does not require renaming `dev`; there are no inherited
  workspace dependencies whose default features change meaning. CI's new default
  of disabling incremental compilation does not change release benchmark settings.
- **Compiler and language:** no affected extern statics, runtime symbols,
  `no_mangle` generics, inline-module path attributes, or `cfg_select!` branches
  were found. Macro-expression semicolons and inferred-pattern diagnostics are
  covered by compiling every workspace target with warnings denied.
- **Libraries:** no deprecated legacy floating-point or integral module constants
  were found in repository Rust. New boxed-array iteration, collection, string,
  filesystem, raw-pointer, allocation, and C-variadic APIs do not simplify a
  current library path. Unsafe remains forbidden.
- **Rustdoc:** retain warning denial and execute both doctest configurations to
  cover unused footnotes, unattached attributes, and changed doctest filtering.
- **Targets:** retain Ubuntu, macOS, and Windows MSVC CI. No new target promise is
  implied by Rust's platform promotions.

[release]: https://github.com/rust-lang/rust/releases/tag/1.99.0
[cargo]: https://doc.rust-lang.org/cargo/CHANGELOG.html#cargo-199-2026-10-01
[llvm]: https://github.com/rust-lang/rust/pull/158734
[range]: https://github.com/rust-lang/rust/pull/155114

The actionable diagnostics were Clippy's
[`manual_bit_width`](https://rust-lang.github.io/rust-clippy/rust-1.99.0/index.html#manual_bit_width),
[`assert_is_empty`](https://rust-lang.github.io/rust-clippy/rust-1.99.0/index.html#assert_is_empty),
and strengthened `must_use_candidate` checks. Four integer expressions now use
`bit_width` (three library sites and one test oracle), seven accessors/fixture
constructors gain `#[must_use]`, and eleven empty-array assertions show values on
failure. No lint suppression was added. `bit_width(x)` equals
`BITS - leading_zeros(x)`, including zero; the existing nonzero precondition of
the rounding helper still protects its subsequent subtraction by one. These
source changes are shared by both compiler measurements.

General matrix-vector and matrix-matrix multiplication are not implemented in
this source revision. Vector-angle issue #249 was open at measurement time.
They are not reported as measured APIs. The existing fixed-size determinants,
factorizations, vector reductions, scale-aware norms, interval signs, certified
linear forms, and rational/Bareiss controls are the applicable workloads.

## Correctness and portability

The unmodified 1.98.1 starting point passed all 893 runnable Rust tests through
`just test-rust-ci`. After the lint fixes, the same 893 tests passed independently
under both compilers with `PROPTEST_CASES=256` and `PROPTEST_RNG_SEED=251`, in
release mode with every workspace feature enabled. The seed fixes the generated
corpus; it is not a claim of exhaustive binary64 coverage.

The existing independent checks include:

- Integer oracles for vector dot products, squared norms, and Gram construction.
- Rational Gaussian elimination and permutation expansion for determinant,
  interval containment, exact-sign, and certified-bound properties.
- Exact squared-midpoint comparisons for the norm overflow boundary and
  explicit binary64 bit patterns for signed-zero/subnormal conversion ties.
- Singular, near-singular, pivoting, Hilbert, asymmetric, indefinite,
  cancellation, mixed-magnitude, NaN/infinity, and overflow/underflow cases with
  typed error fields checked.
- Deterministic comparison fixtures, including the Delaunay-style scaled-norm
  reference and lifted interval predicate cases.

`just check` and the final `just ci` passed on Rust 1.99. The latter includes
all 893 runnable Rust tests, both default/exact doctest configurations, Python
consumer tests, warning-denying Clippy/rustdoc, static-analysis fixtures,
benchmark compilation, and runnable examples. The source and benchmark digest
below remained unchanged through that pass and the measurement phase.

No correctness failure or documented-contract change was observed. Agreement
between compiler runs supplements these independent oracles; it does not replace
them. Rust 1.99 `cargo check --lib --no-default-features` also passed both with
and without `--features exact` for `x86_64-unknown-linux-gnu` and
`x86_64-pc-windows-msvc`. Native Linux/Windows execution remains the hosted CI
matrix's responsibility; cross-compilation checks do not claim execution coverage.

An isolated copy of downstream Delaunay commit
`6e13461453040eec7b505f1f27607ac954e11257` passed 721 geometry unit tests and
18 predicate property tests with Rust 1.99 and a Cargo patch to this local
la-stack. The sibling checkout was clean and remained unchanged. These tests
exercise filtered/exact orientation, in-sphere classification, dimensional
dispatch, and geometry behavior through the downstream API.

## Measurement provenance

Measurements use the same Apple M4 Max, `aarch64-apple-darwin`, macOS 27.0.1
(26A434), with Cargo's bench/release optimization, fat LTO, and one codegen unit.
No target-CPU override or relaxed floating-point flag is used.

| Compiler | Commit | LLVM |
|----------|--------|------|
| Rust 1.98.1 | `48a229ceaefd4985c50990b14116b6d856af0985` (2026-09-01) | 22.1.8 |
| Rust 1.99.0 | `b940084d7eb6a299eb4bfeb8e34901bc051e7ac4` (2026-09-28) | 23.1.1 |

The base commit is `1615b1bcad80a209ced1116e81ffd545e26d85b3`, plus this
migration's source changes and the maintainer's existing staged dependency
updates. Both compilers use that same working-tree lockfile. Historical release
numbers in generated performance artifacts are intentionally unchanged.

- Cargo.lock SHA-256:
  `c23b100c8f28ee06cbcf789186e556977ee2088024c205147cdedb68aa94cba3`.
- Source/benchmark inventory SHA-256:
  `1a51ad89102f6264ccfadd3021b3263b0dee436eda6064e868b3de14b40765c5`,
  computed by `rg --files src benches | LC_ALL=C sort | xargs shasum -a 256 |
  shasum -a 256` from the repository root.
- Code-generation probe SHA-256:
  `7577b1a08aff650b516d0c0fef0b82a7ba0f337dfe64414a37d3e60975a4966a`.

The [evidence archive](rust-1.99/evidence.tar.gz) retains both compilers' raw
Criterion samples and estimates, comparison logs (including significance and
outlier diagnostics), assembly/LLVM IR, and validation logs. Its SHA-256 is
`b6830437c51c4cb11cdae1fe14dd16a8c5c7bd324707dff3d3028de6121e7fa3`.
Each summarized measurement was checked with the repository's
`criterion_measurements.validate_measurement` before inclusion.

The old-compiler commands explicitly use `--ignore-rust-version` to run the
same final source after the manifest's MSRV bump. This is an experimental
comparison, not an extension of the new support promise. Separate target
directories prevent compiler artifacts from mixing. Tests and builds finish
before each measurement; timing jobs run serially, with no overlapping audit
builds or validators. Other desktop activity is not controlled.

The [retained measurement script](rust-1.99/measure.sh) runs from the repository
root after correctness validation. All five suites use `bench,exact`, 100
Criterion samples, a one-second warmup, a three-second measurement window,
95% confidence intervals, and the default 1% noise threshold. Inputs are the
existing validated fixtures; no construction or correctness assertions are
added to timed closures. The first pass saves `rust-1.98.1`; the second compares
Rust 1.99 against it in `target/toolchain-audit/criterion`.

The selected 54 cases cover all D=2-5 routine determinant, dot, norm,
squared-norm, LU, LDLT, factor-plus-solve paths; four D=4 norm regimes; six exact
or rational controls; three interval signs; five certified/plain linear forms;
and four D=4 Gram construction regimes. This is a representative compiler study,
not a replacement for the complete release-signal suite.

## Measurements

Both measurement passes were completed on 2026-10-05. Their purpose is attribution
of compiler effects, not a requirement that every benchmark improve before the
migration proceeds.

The complete [54-case comparison](rust-1.99/measurements.csv) retains unrounded
means, medians, absolute mean confidence intervals, and Criterion's relative
mean-change intervals. The table below shows representative first-pass means in
nanoseconds; negative change means a smaller runtime. These are throughput
measurements of the existing Criterion closures, not isolated-call latency.

| Case | Rust 1.98.1 ns | Rust 1.99.0 ns | Mean change | Criterion change 95% interval |
|------|---------------|---------------|-------------|------------------------------|
| `d4/la_stack_det` | 2.608 | 2.613 | +0.2% | -0.1% to +0.6% |
| `d5/la_stack_dot` | 0.850 | 0.846 | -0.5% | -1.2% to +0.3% |
| `d3/la_stack_norm2_sq` | 0.460 | 0.507 | +10.1% | +9.7% to +10.5% |
| `d4/la_stack_norm2` | 2.364 | 2.446 | +3.5% | +1.9% to +5.2% |
| `d4/la_stack_norm2_scenario_sparse` | 1.542 | 1.666 | +8.1% | +5.4% to +10.7% |
| `d4/la_stack_norm2_scenario_wide_dynamic_range` | 2.252 | 2.251 | -0.1% | -0.6% to +0.5% |
| `d4/la_stack_ldlt` | 24.388 | 18.280 | -25.0% | -27.5% to -22.7% |
| `d4/la_stack_ldlt_solve` | 24.374 | 18.321 | -24.8% | -25.1% to -24.5% |
| `d5/la_stack_ldlt` | 44.426 | 25.589 | -42.4% | -42.6% to -42.2% |
| `d5/la_stack_ldlt_solve` | 47.768 | 31.708 | -33.6% | -34.2% to -32.7% |
| `d5/la_stack_lu_solve` | 46.878 | 49.269 | +5.1% | +4.5% to +5.8% |
| `exact_d4/det_exact` | 1009.513 | 983.616 | -2.6% | -3.0% to -2.1% |
| `exact_d5/det_exact` | 3526.345 | 3605.950 | +2.3% | +1.3% to +3.5% |
| `exact_d5/solve_exact` | 159040.988 | 160320.874 | +0.8% | +0.5% to +1.1% |
| `rational_input_d5/solve_row_cleared_bareiss` | 9582.322 | 10061.106 | +5.0% | +4.7% to +5.4% |
| `interval_det_sign/d4_inconclusive_lifted` | 122.632 | 113.380 | -7.5% | -8.4% to -6.7% |
| `interval_det_sign/d7_conclusive_lifted` | 1386.706 | 1366.206 | -1.5% | -2.5% to -0.5% |
| `linear_form_d4/dot_bounded_well_separated` | 32.255 | 31.469 | -2.4% | -2.7% to -2.1% |
| `linear_form_d4/dot_difference_bounded_well_separated` | 65.078 | 63.844 | -1.9% | -2.3% to -1.4% |
| `gram/4x4/orthogonal/la_stack` | 7.847 | 7.684 | -2.1% | -2.5% to -1.6% |

The [downstream comparison](rust-1.99/downstream.csv) retains all eight Delaunay
cases. It uses that crate's `perf` profile (thin LTO, one codegen unit, line-table
debug information), default crate features plus la-stack's `exact` dependency,
and the same 100/1s/3s Criterion configuration. Its lockfile is resolved once with
`--config 'patch.crates-io.la-stack.path="/Users/adam/projects/la-stack"'` in the
isolated copy; both timing runs use that lock with `--locked --offline`.

Seventeen library cases were selected for longer reverse-order repeats:
absolute first-pass mean changes of at least 5%, regressions whose interval
crossed 5%, the noisy D=4 dot row, and the adjacent D=5 LU factorization. All
eight downstream cases were repeated. The [library repeat commands](rust-1.99/repeat.sh)
use 100 samples, a two-second warmup, and five seconds of measurement, running
1.99 before 1.98. This selection rule preserves unfavorable and noisy results.

The [repeat data](rust-1.99/repeat.csv) use the same canonical direction as the
first pass: positive means Rust 1.99 is slower. Criterion's reverse-order
relative interval is transformed with `1 / (1 + change) - 1`, swapping its
endpoints. Results are kept separate; they are not pooled across sessions.

| Case | First-pass change | Reverse-repeat change | Repeat 95% interval |
|------|-------------------|-----------------------|---------------------|
| `d3/la_stack_dot` | +53.5% | -0.3% | -0.6% to +0.1% |
| `d3/la_stack_ldlt` | +10.2% | -11.6% | -12.3% to -10.9% |
| `d3/la_stack_ldlt_solve` | -14.0% | -21.4% | -22.2% to -20.6% |
| `d3/la_stack_lu` | -28.4% | -4.1% | -4.9% to -3.3% |
| `d3/la_stack_norm2` | -9.6% | -1.2% | -3.8% to +1.4% |
| `d3/la_stack_norm2_sq` | +10.1% | +6.9% | +6.3% to +7.4% |
| `d4/la_stack_dot` | -1.4% | +3.4% | +2.4% to +4.4% |
| `d4/la_stack_ldlt` | -25.0% | -23.0% | -23.6% to -22.4% |
| `d4/la_stack_ldlt_solve` | -24.8% | -28.4% | -29.1% to -27.7% |
| `d4/la_stack_norm2` | +3.5% | -1.5% | -3.1% to +0.2% |
| `d4/la_stack_norm2_scenario_sparse` | +8.1% | +0.3% | -2.9% to +3.4% |
| `d5/la_stack_ldlt` | -42.4% | -42.2% | -42.5% to -41.8% |
| `d5/la_stack_ldlt_solve` | -33.6% | -32.7% | -33.3% to -32.2% |
| `d5/la_stack_lu` | +3.3% | +4.7% | +4.2% to +5.3% |
| `d5/la_stack_lu_solve` | +5.1% | +13.4% | +12.5% to +14.3% |
| `interval_det_sign/d4_inconclusive_lifted` | -7.5% | -4.8% | -6.4% to -3.5% |
| `rational_input_d5/solve_row_cleared_bareiss` | +5.0% | +0.9% | +0.3% to +1.4% |

The D=4–5 LDLT improvements persist in both orders: factorization is about
23–42% faster in the repeat and factor-plus-solve about 28–33% faster.
The D=3 LDLT first-pass regression reverses, so it is not a stable regression.
The large D=3 dot change disappears; D=3/D=4 norms and sparse D=4 norms also
fail to retain their first-pass changes. Those cases demonstrate session noise
and should not support a broad compiler-speed claim.

The unfavorable persistent library results are D=3 squared norm
(+10.1%, then +6.9%), D=5 LU (+3.3%, then +4.7%), and D=5 LU plus solve
(+5.1%, then +13.4%). Their direction repeats, but the LU-solve magnitude is
session-sensitive. The initial +5.0% rational solve change falls to +0.9%
(within the default 1% practical-noise threshold by its mean); it is not a
durable 5% regression. The inconclusive D=4 interval filter remains faster,
although the repeat gain is 4.8% instead of 7.5%.

The unchanged or modestly changed exact, interval, Gram, and certified-linear-
form controls in the full first-pass data guard against interpreting the LDLT
gains as a uniform machine-wide speedup. Only the selected cases were repeated;
small first-pass effects elsewhere are descriptive observations.

All eight [downstream repeats](rust-1.99/downstream-repeat.csv) are retained:

| Delaunay case | First-pass change | Reverse-repeat change | Repeat 95% interval |
|---------------|-------------------|-----------------------|---------------------|
| `geometry/axis_simplex/circumradius/3d` | -1.9% | -2.8% | -3.2% to -2.4% |
| `geometry/axis_simplex/volume/3d` | +14.5% | +17.0% | +16.5% to +17.5% |
| `orientation/exact_nonzero/3d` | -16.8% | -12.6% | -13.1% to -12.1% |
| `orientation/exact_nonzero/4d` | -4.7% | -0.2% | -1.0% to +0.5% |
| `orientation/exact_zero/3d` | +11.6% | +15.5% | +14.1% to +17.0% |
| `orientation/exact_zero/4d` | +13.5% | +13.3% | +12.7% to +13.9% |
| `orientation/well_conditioned/3d` | +35.8% | +36.3% | +20.2% to +53.7% |
| `orientation/well_conditioned/4d` | +4.9% | +4.4% | +3.2% to +5.6% |

The volume and exact-zero orientation slowdowns persist, as does the exact
nonzero D=3 improvement. D=3 volume uses Delaunay's own scalar triple-product
implementation, so that change is a downstream compiler effect rather than
a la-stack kernel regression. The orientation fixtures explicitly assert the
expected filter/fallback path and analytical sign before timing under both
compilers; a changed scientific result or unintended fallback does not explain
the timing differences.

The well-conditioned D=3 means are dominated by a slow tail: repeat medians
are 10.545 ns on 1.98 and 10.432 ns on 1.99, despite means of 10.474 and
14.272 ns. Both sessions show that discrepancy. Retain the mean regression,
but do not interpret it as a uniform 36% increase in ordinary-call cost.
The precise backend causes of the remaining downstream differences are
unresolved performance questions. They are accepted migration costs, with
correctness preserved; this study does not claim that every workload improves.

The [downstream command script](rust-1.99/downstream.sh) accepts an isolated
Delaunay `Cargo.toml` path and optional `repeat` mode. Reproduce with the
downstream commit above, resolve the local Cargo patch once, and build both
compiler/profile combinations before timing. The sibling checkout is not a
measurement output directory.

## Code generation

The retained [probe source](rust-1.99/codegen.rs) exposes thirteen public operations
without modifying their implementations or timed benchmark closures. Compile it
against each compiler's `bench,exact` library with `rustc --edition=2024
--crate-type=lib --crate-name codegen_audit -Dwarnings -Copt-level=3 -Clto=fat
-Ccodegen-units=1 --emit=asm,llvm-ir`, adding the matching dependency directory
with `-L dependency=...` and `--extern la_stack=...`.

Static instruction counts include cold error branches and exclude assembler
directives. They describe these wrappers, not dynamic instruction counts or a
proof of benchmark speed:

| Probe | 1.98.1 instructions | 1.99.0 instructions |
|-------|---------------------|---------------------|
| Certified dot D=4 | 91 | 91 |
| Determinant D=4 | 242 | 239 |
| Dot D=4 | 60 | 57 |
| Exact determinant D=5 | 789 | 790 |
| Interval determinant D=4 | 497 | 496 |
| LDLT + solve D=4 | 560 | 492 |
| LDLT + solve D=5 | 742 | 626 |
| LU + solve D=4 | 543 | 519 |
| LU + solve D=5 | 639 | 606 |
| Norm D=4 | 37 | 37 |
| Rational solve D=5 | 2363 | 2373 |
| Squared norm D=3 | 45 | 43 |
| Squared norm D=5 | 22 | 21 |

Both compilers retain the four dependent scalar `fmadd` instructions for dot
and five for squared norm, preserving reduction order. LLVM 23 uses `ubfx` and
an immediate exponent comparison for finite checks, replacing a full-width mask
and materialized constant. Similar changes shorten LU and direct-determinant
error checks. The scale-aware norm retains its divisions, FMA recurrence, square
root, and overflow-boundary fallback.

LDLT's changed inclusive-range lowering and unrolling produce smaller D=4 and
D=5 wrappers, consistent with the repeated gains and the upstream range
optimization. This comparison changes the whole Rust toolchain, so it does not
isolate LLVM as the sole cause. LLVM 23 also infers `nnan` on sixteen products
of already-checked finite multipliers and pivots across those two wrappers.
Those products can overflow to infinity but cannot produce NaN; later pivot,
multiplier, and solve-result checks remain present. This inferred fact does not authorize
reassociation or global non-finite assumptions. Neither probe module contains
`fast`, `reassoc`, `ninf`, `nsz`, `arcp`, `afn`, or `contract` arithmetic flags.
No source-level floating-point optimization was adopted.

The D=3 squared-norm probe retains the same three ordered FMAs; its finite
check changes from a mask/materialized constant to an exponent extraction.
The D=5 LU probe also gets shorter through finite-check lowering, yet runs
slower in the measured harness. Static code size is therefore not a speed
oracle. Scheduling, inlining, and layout in the timed closures are possible
contributors, but these wrappers do not establish a causal explanation for
either regression. The rational control changes little in static size and
loses its material timing regression on repeat. No speculative source
optimization is justified by this evidence.

## Decision and future baseline

Rust 1.99.0 is the new declared MSRV and shared contributor/CI/release toolchain.
Proceed with the migration: validation found no correctness or compatibility
blocker. The measured regressions above remain documented optimization
candidates, not reasons to weaken numerical guarantees or defer the upgrade.

The strict arithmetic policy remains in force. Future #250 measurements must use
this compiler, the same valid fixtures and profiles, and fresh paired runs on
the same machine when evaluating source changes. The Rust 1.98 release report
cannot substitute for that comparison.

The forthcoming v0.4.6 / Rust 1.98.1 versus v0.4.7 / Rust 1.99.0 release
comparison should report both toolchains and measure their combined
compiler-and-library effect. That is a valid release comparison. This separate
same-source study provides attribution when needed; it does not replace or
delay the normal release benchmarks.
