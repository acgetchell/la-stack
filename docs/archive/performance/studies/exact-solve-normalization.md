# Exact Solve Normalization

This October 18, 2026 study retains integer back-substitution over a shared
denominator for general exact solves. Correctness and canonical output remain
the acceptance criteria. The [mathematical basis](../../../mathematical_basis.md#exact-arithmetic-over-binary64-inputs)
owns the invariant and denominator-growth proof.

## Contents

- [Decision and correctness](#decision-and-correctness)
- [Measurements](#measurements)
- [Diagonal control](#diagonal-control)
- [Provenance and replay](#provenance-and-replay)

## Decision and correctness

The old back-substitution reduces rational products and differences throughout
its quadratic inner sum. The replacement accumulates integer numerators over
one positive denominator, reduces once per solved row, and canonicalizes each
output. The denominator grows only when a newly solved component requires it.
This removes repeated GCD work without introducing rounding or changing pivot
selection, exact singularity metadata, or strict/rounded conversion behavior.
The rational triangular shortcut and exact determinant implementation are
unchanged.

Sixteen added tests cover manufactured systems at D=0–8 and independent
Cramer's-rule comparisons at D=2–8. They check exact residuals and canonical
numerators/denominators, including zero components, fractional solutions,
negative pivots, row permutations, and large common factors. Existing
subnormal, mixed-scale, Hilbert, singularity, and conversion-boundary tests
also pass.

The focused release-mode gate passed 728 tests before and 744 after, with
`bench,exact` enabled. Candidate property runs used `PROPTEST_CASES=64` and
`PROPTEST_RNG_SEED=285001`. Final `just ci` passed 1,020 Rust tests, 181 Python
tests, default/exact doctests, Clippy, static analysis, documentation builds,
benchmark compilation, and examples. The implementation and test patch was
verified byte for byte against the validated snapshot when preparing this PR.
The additional study files received focused documentation and data checks.
Linux and Windows execution remains with GitHub CI.

An earlier fixed-denominator prototype retained unnecessary determinant factors.
Although exact, it increased the 1024-bit D=4 solve time by about 9% and raised
wide-input allocation traffic. Its source and measurements are retained for
audit; it is not the implementation selected here.

## Measurements

Both revisions use Rust 1.99.1, Criterion 0.8.2, `bench,exact`, the same locked
dependencies and fixtures, fat LTO, and one codegen unit. The host is Apple M4
Max / aarch64 macOS 27.0.1 (`26A434`), with LLVM 23.1.1. Benchmark inputs are
validated independently outside timing, and timings include output destruction.
No compilation or correctness suite from this task overlaps timing. Normal host
activity, CPU placement, and frequency are uncontrolled.

The broad comparison covers 47 cases: 39 exercise the changed backend, and
eight are triangular controls. All 39 affected cases improve in the first pass.
A reverse-order repeat covers 21 cases, including the wide-input cases and all
16 structured sparse/dense/triangular cases. Each run uses 30 samples, a
one-second warm-up, and a two-second measurement target.

Arithmetic means per complete solve, in microseconds:

| Fixture | Before | After | Time reduction |
|---|---:|---:|---:|
| Binary64 input, D=2 | 7.189 | 3.270 | 54.5% |
| Binary64 input, D=3 | 30.780 | 9.477 | 69.2% |
| Binary64 input, D=4 | 80.040 | 19.546 | 75.6% |
| Binary64 input, D=5 | 158.095 | 31.966 | 79.8% |
| Small rational input, D=8 | 45.270 | 18.284 | 59.6% |
| 1024-bit rational input, D=4 | 46.057 | 35.854 | 22.2% |
| 1024-bit rational input, D=8 | 866.680 | 678.531 | 21.7% |
| Dense rational Hilbert input, D=8 | 30.915 | 11.342 | 63.3% |
| Sparse rational input, D=8 | 20.963 | 7.543 | 64.0% |

The [complete comparison](exact-solve-normalization/comparison.json) retains
every case, arithmetic means, marginal 95% Criterion intervals, and the
reverse-order repeats. Ratios are descriptive, not paired confidence intervals
or promises for other hardware or inputs. Random-corpus rows time batches of
50 inputs; the table above contains only individual solves. This is a focused
optimization study, not a release comparison or a workload-weighted result.

The [allocation comparison](exact-solve-normalization/allocations.json) uses
separate counting-allocator probes. Small rational D=8 falls from 461 to 412
allocations; 1024-bit D=8 falls from 4,113 to 3,816. Across 21 fixtures, the only
increase is small rational D=6: 78 to 79 allocations, or 1,368 to 1,384 cumulative
allocated bytes. Its time falls from 18.449 to 6.309 microseconds. Allocated
bytes describe cumulative traffic, not peak memory.

## Diagonal control

Initial rational diagonal D=8 means increased by 1–3%, despite bypassing the
changed backend. Disassembly confirms an indirect compiler change: the enclosing
solve's stack frame grows by 96 bytes and register allocation differs.
Unchanged source alone therefore cannot establish unchanged performance.

A subsequent predeclared experiment uses eight balanced ABBA/BAAB comparison
blocks (32 processes) and two same-binary placebo blocks (eight processes).
Each process uses 100 samples, two seconds of warm-up, and five seconds of
measurement. The aggregate means are 780.921 ns before and 779.730 ns after.
The geometric mean block ratio is 0.998461, with a block-bootstrap 95% interval
of [0.991337, 1.006374]: 0.87% faster to 0.64% slower. Five blocks favor the
candidate; three favor the baseline. The earlier slowdown does not reproduce
consistently, so no compiler attributes or additional algorithm changes were
made to chase it.

[All process means and block results](exact-solve-normalization/diagonal.json)
are retained. The resampling unit is an entire four-process block, not individual
Criterion samples. The sequence is balanced, not randomized; the interval is
conditional on this host and these binaries and does not exclude smaller or
platform-specific effects. Same-binary placebo label differences are −0.72%
and −0.38%. These controls are descriptive and are not used to adjust the main
comparison. Earlier unfavorable measurements remain in the evidence archive.

## Provenance and replay

The baseline is `99b060307476fb3df0f25fad5bab67dfdff0f5b2`. The
[evidence archive](exact-solve-normalization/evidence.tar.gz) contains raw
Criterion records, validation and timing logs, exact source snapshots, the
implementation/test patch, source/harness/binary hashes, and the original
experiment and analysis scripts. It excludes compiled binaries. The candidate
`src/exact.rs` SHA-256 is
`7226789c37158eb87b43260c28e32dccd93b7ab4f43aea9ccd6ececf3d51bd5a`.
The archived commands document the original uncommitted experiment; later
archival documentation changes do not affect the measured library or harness.
The `dprint` 0.64.0 tool-pin update was present during validation and is included
in this PR. It does not participate in Rust builds or timing.

Build each variant separately using
`cargo bench --locked --features bench,exact --bench exact --bench rational_solve --no-run`.
Preserve the resulting executables as `before-exact`, `before-rational-solve`,
`after-exact`, and `after-rational-solve` under
`target/exact-normalization-20261018/`. Extract the archive there before using
its `measure.sh`, `repeat.sh`, or `diagonal-investigation/measure.sh`; use a fresh
copy to avoid overwriting the retained measurements. The archived `commands.md`
contains all filters, gates, settings, and run order.

To verify archive integrity and replay the two timing summaries without running
benchmarks, execute from the repository root:

```bash
study_dir=docs/archive/performance/studies/exact-solve-normalization
(cd "$study_dir" && shasum -a 256 -c SHA256SUMS)
replay_dir="$(mktemp -d)"
tar -xzf "$study_dir/evidence.tar.gz" -C "$replay_dir"
rg --files --hidden --no-ignore "$replay_dir/criterion" \
  | rg '/(before|after|before-repeat|after-repeat)/estimates.json$' \
  | xargs jq -Sn -f "$replay_dir/compare.jq" > "$replay_dir/comparison-replay.json"
cmp "$study_dir/comparison.json" "$replay_dir/comparison-replay.json"
uv run --locked python "$replay_dir/diagonal-investigation/analyze.py" \
  > "$replay_dir/diagonal-replay.json"
cmp "$study_dir/diagonal.json" "$replay_dir/diagonal-replay.json"
jq -n --slurpfile before "$replay_dir/before-allocations.json" \
  --slurpfile after "$replay_dir/after-allocations.json" \
  -f "$replay_dir/compare-allocations.jq" > "$replay_dir/allocations-replay.json"
cmp "$study_dir/allocations.json" "$replay_dir/allocations-replay.json"
```

The diagonal analysis uses 100,000 whole-block bootstrap replicates and seed
285002. These commands reproduce the published JSON files byte for byte; the
original per-run records and marginal intervals remain available for other
analyses.
