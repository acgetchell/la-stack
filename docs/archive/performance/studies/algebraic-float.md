# Algebraic Floating-Point Reduction Evaluation

The local evaluation for [#250](https://github.com/acgetchell/la-stack/issues/250)
records the decision to **retain strict arithmetic**. The selected algebraic reduction
can suppress a required overflow error on finite inputs, including D=4 and
D=5. This is a negative adoption decision; no library API, arithmetic policy,
or Semgrep exception changes are proposed.

The decision is supported by a native macOS counterexample and compiler-contract
analysis. On 2026-10-15 the maintainer rejected adoption and ended the evaluation,
including unperformed native Linux/Windows prototype runs. The
[validation decision](../../../validation/algebraic-floating-point.md) records
that scope closure; no claim is made that those runs passed.

## Contents

- [Candidate selection](#candidate-selection)
- [Numerical findings](#numerical-findings)
- [Measurement method and results](#measurement-method-and-results)
- [Provenance and reproduction](#provenance-and-reproduction)
- [Validation and scope closure](#validation-and-scope-closure)
- [Decision](#decision)

## Candidate selection

Select one existing kernel: the positive row-sum reduction in
`Matrix::norm_inf`. This fits #250's reductions category without changing an
FMA-based algorithm. The complete operation takes absolute values, sums each
row, checks overflow, and selects the largest sum. No certified rounding bound
is promised, but computed non-finite values must retain typed error provenance.

The [existing performance record](../../../performance.md) puts strict D=2–5
matrix infinity norms at roughly 0.6–3.5 ns, increasing to about 8.4 ns at D=8.
These historical measurements justify investigating the addition dependency
chain; they are **not** the baseline for this experiment. The public operation
already has a comparison fixture in `benches/vs_linalg.rs`, and its ordinary
additions offer a more direct candidate than an FMA substitution.

The disposable prototype changes only the success-path row additions to
`algebraic_add`. It accepts an already-validated `Matrix<D>`, retains absolute
values, finite-result checks and maximum selection, and delegates a non-finite
completed row to the original strict method for error recovery. That fallback
is deliberately conservative but cannot catch an incorrectly finite result.
Both versions return the complete `Result<f64, LaError>` in timed closures.

Other plausible candidates are excluded before implementation:

- **Dot products and squared norms:** the current APIs explicitly accumulate
  with ordered `mul_add`. The [compiler migration study](rust-1.99.md) observes
  the same scalar FMA chain in generated code. Breaking that chain by replacing
  FMA with algebraic multiplication/addition violates this issue's exclusion.
- **General matrix multiplication and matrix-vector multiplication:** these
  are not existing APIs at the captured source revision.
- **Predicates, certified/interval reductions, exact fallbacks, tolerance and
  branch decisions:** excluded by the issue. The candidate is not substituted
  into those paths or into a downstream classifier.

## Numerical findings

### A finite result hides required overflow

Let `a = 2^969` and `M = f64::MAX = 2^1024 - 2^971`. Use the row
`[a, a, M, 0]`, or append another zero for D=5. All entries are finite.

The strict recurrence first obtains `a + a = 2^970`, then adds `M`. The exact
sum is the rounding midpoint `2^1024 - 2^970`; round-to-nearest, ties-to-even
overflows. The original API returns:

```text
LaError::NonFinite
  origin: Computation { operation: MatrixInfinityNorm }
  location: MatrixCell { row: 0, col: 2 }
```

The optimized algebraic prototype instead returns the finite bit pattern
`0x7fefffffffffffff`, or `f64::MAX`, in both dimensions. Reassociation can absorb
each `a` separately when added to `M`, avoiding the overflowing strict prefix.
The retained inlined experiment demonstrates the changed classification; an
out-of-line probe need not choose the same association. A final `is_finite`
check and a strict fallback triggered by non-finite results cannot detect this
failure. No input validation or exact fallback was removed to produce it.

### Finite accuracy and boundary checks

An independent `BigRational` oracle sums the absolute values of the exact binary64
inputs and selects the largest rational row sum. The experiment covers D=0, 1,
2, 3, 4, 5, 8 and 16. There is no finite maximum supported dimension for the
const-generic API; D=0 and D=1 cover its lower boundaries, while D=8 and D=16
supplement the primary D=2–5 scope.

Explicit cases cover signed zero, minimum subnormals, minimum normal values,
mixed and extreme magnitudes, absorption, definite overflow, and near-overflow
permutations. NaN and both infinities are rejected at construction with exact
input origin and coordinates. Cancellation of signed terms is not applicable:
this kernel sums absolute values. Underflow from multiplication or division is
also absent; subnormal addition is exercised directly.

For `[2^53, 1, -1, 1]`, the exact absolute row sum is `2^53 + 3`. Strict
accumulation returns `2^53`, whereas the algebraic prototype returns `2^53 + 2`.
The latter is closer to the exact sum, but changes the ordered result already
protected by `matrix_norm_inf_preserves_left_to_right_row_sum_order`. The
overflow counterexample is stronger evidence against adoption than a bitwise
difference alone.

A deterministic xorshift corpus with seed 250 adds 1,000 matrices per dimension,
spanning finite exponent fields, signs and mantissas. It reports no random
classification changes, but 2 value changes at D=8 and 9 at D=16. Random success
does not erase the explicit D=4/D=5 failure. The archive retains all explicit
case outputs, random differences, exact-oracle relative errors, and the source
that reconstructs every input.

### Compiler guarantees

[Rust 1.99's algebraic-operator contract][algebraic-contract] allows
reassociation and unspecified precision, permits unexpected special-value
behavior, and does not promise repeatable results even within one run. Finite
inputs, an output range check, and successful observed tests therefore cannot
establish a portable accuracy guarantee for this candidate. The absence of a
certified bound on the existing API does not waive its ordered arithmetic or
typed overflow contract.

Generated LLVM IR for the prototype uses `reassoc nsz arcp contract` additions.
The native D=4 probe groups additions into pairs, whereas the strict probe
retains the ordered chain. These probes explain available transformations, not
the exact schedule of every inlined call or a causal timing proof.

[algebraic-contract]: https://doc.rust-lang.org/std/primitive.f32.html#algebraic-operators

## Measurement method and results

Measurements compare the captured strict public method and the isolated
prototype in the same optimized executable. The shared fixture reproduces
`benches/common/vs_linalg.rs::matrix_entry`, including its diagonal `mul_add`;
an additional dense signed-integer fixture has exactly representable sums.
Both use the same preconstructed matrix per pair.

Before timing, both results must be finite and agree with the independent exact
reference within the ordinary summation envelope `γ_D × exact_norm`, where
`γ_D = D u / (1 - D u)` and `u = 2^-53`. This is an empirical fixture gate,
**not** a guarantee supplied by algebraic operations. Adversarial failing cases
are never timed as successful workloads.

Criterion uses `iter`, 100 samples, 100,000 resamples, a one-second warmup,
three-second measurement target, 95% confidence intervals and a 1% noise
threshold. The first pass runs strict then algebraic; the second reverses that
order. All ten dimension/fixture pairs are repeated, with no selection based on
favorable results. Builds and validation finish before measurement; timing jobs
run sequentially. Other desktop activity is uncontrolled.

All 20 comparisons are retained in
[measurements.csv](algebraic-float/measurements.csv), including unrounded means,
mean confidence intervals and medians. The shared-fixture results below report
mean nanoseconds with marginal 95% confidence intervals. Negative changes mean
less time. Changes are descriptive point-estimate ratios, not paired confidence
intervals or application-level speedups.

| D | Forward strict ns [95% CI] | Forward algebraic ns [95% CI] | Forward change | Reverse change |
|---|---------------------------|------------------------------|----------------|----------------|
| 2 | 0.928 [0.923, 0.932] | 0.952 [0.949, 0.954] | +2.6% | +1.4% |
| 3 | 1.387 [1.386, 1.389] | 1.766 [1.757, 1.776] | +27.3% | +28.0% |
| 4 | 2.413 [2.409, 2.420] | 1.856 [1.849, 1.864] | -23.1% | -21.6% |
| 5 | 3.529 [3.524, 3.535] | 2.633 [2.630, 2.636] | -25.4% | -25.0% |
| 8 | 8.737 [8.713, 8.764] | 5.281 [5.272, 5.294] | -39.6% | -39.6% |

The D=4/D=5 reductions save roughly 0.5–0.9 ns per measured operation in both
orders; D=8 saves about 3.4 ns. The integer fixtures confirm the direction:
D=4 improves 23.5–23.6%, D=5 about 25%, and D=8 about 39–40%. The unfavorable
results also repeat: D=2 slows roughly 1–3%, and D=3 about 27–29%, across both
fixtures and orders. Thus algebraic arithmetic does unlock a material local
timing benefit at some dimensions, but not uniformly. These are throughput
measurements in Criterion's loop, not isolated-call latency.

These timings characterize an **ineligible prototype on the checked fixtures**.
They cannot support a claim of a correctness-preserving library optimization,
and are not promoted into the release performance record.

## Provenance and reproduction

The strict source is commit `d1f85bb606735d5213a2fa5bd2692e5685a9adc4`, which
includes PR #278 / #246. The worktree was clean before the study. This captures
the landed strict candidate rather than the older #251 source; subsequent
release changes would require a refreshed baseline before any future adoption.

- Compiler: Rust 1.99.0, commit
  `b940084d7eb6a299eb4bfeb8e34901bc051e7ac4`, LLVM 23.1.1.
- Machine: Apple M4 Max, `aarch64-apple-darwin`, macOS 27.0.1 (26A434).
- Profile: optimized release, fat LTO, one codegen unit; no target-CPU override
  or global floating-point flags.
- Library features: no default or optional features. Criterion 0.8.2 and the
  rational oracle are isolated harness dependencies, not library dependencies.
- Source/benchmark inventory SHA-256:
  `eec58161929f41d1a2d183aacca2b83f7af1fcfca9abe0c899fc9dd7a8961b30`.
- Repository Cargo.lock SHA-256:
  `0cd6cc2d847bca04695955c4907d0a61550d4148f5f0568d514d27bf0c8e2c13`.
- Harness Cargo.lock SHA-256:
  `28dea92b824d7856d83a6d68e40fbdb1ded51a947e54c8d8293a038563c7acdf`.

The [evidence archive](algebraic-float/evidence.tar.gz) retains the disposable
harness and its dependency lock, a captured strict-source build tree, numerical
results, complete raw Criterion samples/estimates, assembly/IR, validation logs,
and downstream lockfile. Its SHA-256 is
`5d03f40ecf6df8d73733e79b69d44dbbbb7bec599dac6a566637573942c5cce0`.
The [provenance record](algebraic-float/provenance.json) also records experimental
source hashes and exact validation commands.
Extract the archive **outside** repository-owned Rust. From the
extracted directory, `bash measure.sh` rebuilds, checks the fixtures, and runs
both measurement orders. Offline reproduction requires the locked dependencies
already in Cargo's cache. `summarize.py`, run through this repository's locked
`uv` environment, validates the raw samples and produces the CSV.

## Validation and scope closure

The unchanged strict baseline passed 18 focused unit tests and 5 matrix property
tests selected by `norm_inf`, with `PROPTEST_CASES=256` and
`PROPTEST_RNG_SEED=250`. The isolated numerical experiment ran in release mode
on the native macOS host.

An isolated copy of Delaunay commit
`bde38f4826cc248b425c4300cadd233417f7646c`, patched to the captured la-stack
implementation, passed 728 geometry unit tests, 26 geometry property tests and
18 predicate property tests under the same seed/case settings. These exercise
the **retained strict baseline**. They do not validate algebraic adoption:
the rejected candidate was never installed in a predicate or topology path.
The sibling checkout remained unchanged.

The strict library and algebraic library probes also cross-compiled for
`x86_64-unknown-linux-gnu` and `x86_64-pc-windows-msvc`. Assembly/IR is retained.
Cross-compilation is not native numerical validation or performance evidence.
Native Linux/Windows execution remains unperformed; no suitable native runner
was available in this local session.

Documentation validation passed `just markdown-fix` and `just markdown-ci`
(including spelling). Because ordinary archive content is excluded from the
active-document recipe, the new report and index also passed explicit rumdl
checks with exclusions disabled; the report/provenance passed explicit spelling
and JSON checks. The summary parser verified all 40 saved measurements have
100 positive finite samples and valid mean/median estimates.

The original remaining work was native Linux/Windows execution of the prototype.
The maintainer discontinued it when accepting rejection and requesting closure
of #250 on 2026-10-15. It is not required for an adoption that will not proceed,
and is not transferred to a follow-up issue. The original native-platform
criterion remains unperformed. The archived provenance retains its measurement-
time status; the [validation decision](../../../validation/algebraic-floating-point.md)
owns the final scope disposition.

## Decision

| Function | Performance observation | Numerical result | Disposition |
|----------|-------------------------|------------------|-------------|
| `Matrix::norm_inf` | D=4/5 improves about 22–25%; D=2/3 regresses | D=4/D=5 overflow becomes a finite success; portable accuracy and error guarantees are unsupported | **Retain strict** |

The experiment establishes that algebraic reassociation is not a scientifically
compatible replacement for this kernel, regardless of its fixture timing.
An opt-in name or feature flag would avoid changing old callers but would not
repair the missing numerical contract. No algebraic API or feature is planned.
Keep existing APIs and policy intact.
