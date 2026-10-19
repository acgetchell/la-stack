# Exact Determinant Row Scaling

This October 19, 2026 study retains row-local integer scaling for exact
determinants. The broader comparison supports substantial gains for near-singular
and mixed-exponent inputs, with costs in some simple controls. Correctness remains
mandatory; an individual control regression is not an automatic rejection.

## Contents

- [Decision and correctness](#decision-and-correctness)
- [Experiment](#experiment)
- [Results](#results)
- [Control rechecks](#control-rechecks)
- [Workload tradeoffs](#workload-tradeoffs)
- [Provenance and replay](#provenance-and-replay)

## Decision and correctness

Each row is lifted to integers using its own minimum nonzero binary exponent.
The determinant restores the sum of those row exponents. This reduces integer
growth when the rows admit different binary scales. Direct integer expansions,
Bareiss elimination, floating-point sign filtering, and strict/rounded conversion
algorithms remain unchanged. The [mathematical basis](../../../mathematical_basis.md#exact-arithmetic-over-binary64-inputs)
owns the scaling identity and exponent-range argument.

This revisits the final compact prototype from the October 18 investigation.
Its arithmetic implementation and forced-inlining annotation are preserved;
a targeted Clippy expectation documents the annotation. The earlier rejection
overweighted small controls without evidence about workload prevalence. This
study instead reports absolute costs, repeated comparisons, and hypothetical
mixtures. It does not establish a universal speedup or a typical user mix.

The baseline and candidate each passed 676 focused checks. These include
independent rational Leibniz tests across D=2–5, canonical numerators and
denominators, signs, zero rows, subnormal inputs, and conversion outcomes.
The expanded benchmark fixtures also compare all four determinant APIs against
rational Gaussian elimination and the determinant of an unscaled integer matrix.
Both benchmark executables passed all 354 fixture checks before timing.

Final `just ci` passed on the measured candidate: 1,031 Rust tests, 181 Python
tests, default/exact doctests, Clippy, Semgrep, documentation checks, benchmark
compilation, and examples. Source and binary hashes verified unchanged after
measurement. Linux and Windows execution remains with hosted CI.

## Experiment

The full `exact` suite contains 264 existing cases and 90 added diagnostics.
The new cases measure determinant values and signs at row-exponent spreads of
0, 8, 32, 64, 128, 256, 512, 1024, and 2000 bits across D=2–5. Symmetric
power-of-two scales multiply to one, preserving the exact output determinant.
The base matrices are diagonally dominant integer fixtures with some zero entries.
D=4 column-scaled controls retain large exponent differences within rows; zero
entries can still permit partial row-scale savings.

All 354 cases ran baseline then candidate. All 136 determinant cases ran again
in the opposite order, candidate then baseline. Unchanged solves, floating-point
methods, rational-input algorithms, and standalone conversions remain in the
full comparison. A further four-process recheck covers 17 unchanged controls
whose initial timings were suspicious. No compilation or correctness suite
overlapped timing.

Both binaries use the same Apple M4 Max host, macOS 27.0.1, Rust 1.99.1,
Criterion 0.8.2, `bench,exact`, lockfile, fixtures, fat LTO, and one codegen unit.
Every run uses 30 samples, one-second warm-up, two-second measurement targets,
and no plots. Host activity, CPU placement, and frequency are uncontrolled;
a snapshot during the repeat showed substantial CPU use by other processes.

The [complete comparison](exact-row-scaling/comparison.json) retains all 354
first-pass cases and all 136 repeats, including unfavorable results and marginal
95% intervals. Ratios compare mean point estimates; they are not paired
confidence intervals. Random-corpus rows measure batches of 50 inputs, not
individual operations. No aggregate across unlike operations is interpreted as
an application speedup.

## Results

Selected exact determinant results use reverse-pass means, in nanoseconds:

| Case | Before | After | Less time, first / reverse |
| --- | ---: | ---: | ---: |
| Near-singular D=3 | 339.22 | 195.27 | 46.4% / 42.4% |
| Near-singular positive D=4 | 1,170.60 | 450.98 | 58.1% / 61.5% |
| Original mixed-exponent D=4 | 1,826.73 | 427.28 | 73.9% / 76.6% |
| Row spread 32 bits, D=4 | 557.36 | 430.59 | 15.0% / 22.7% |
| Row spread 64 bits, D=4 | 836.56 | 422.44 | 49.5% / 49.5% |
| Row spread 64 bits, D=5 | 2,046.48 | 1,158.31 | 43.7% / 43.4% |
| Row spread 2000 bits, D=5 | 11,305.46 | 1,145.21 | 90.7% / 89.9% |

The whole row-spread sweep below reports reverse-pass after/before ratios for
`det_exact`; lower is better and values above one mean more time. The first
pass is retained in the comparison JSON. Benefits at 64 bits repeat across
all four dimensions, while 8 bits is insufficient for a consistent benefit
across dimensions. This is a fixture-dependent crossover, not a dispatch rule.

| Spread, bits | D=2 | D=3 | D=4 | D=5 |
| --- | ---: | ---: | ---: | ---: |
| 0 | 1.158 | 1.047 | 1.046 | 1.049 |
| 8 | 0.897 | 0.927 | 0.961 | 1.050 |
| 32 | 0.910 | 0.911 | 0.773 | 0.720 |
| 64 | 0.613 | 0.664 | 0.505 | 0.566 |
| 128 | 0.609 | 0.585 | 0.612 | 0.424 |
| 256 | 0.609 | 0.549 | 0.486 | 0.328 |
| 512 | 0.557 | 0.504 | 0.433 | 0.262 |
| 1024 | 0.513 | 0.441 | 0.364 | 0.179 |
| 2000 | 0.417 | 0.345 | 0.279 | 0.101 |

Simple controls have costs that should remain visible:

| Exact determinant case | Added ns, first / reverse | More time, first / reverse |
| --- | ---: | ---: |
| Equal-scale D=2 row fixture | 10.61 / 8.09 | 21.3% / 15.8% |
| Equal-scale D=3 row fixture | 12.66 / 7.76 | 7.8% / 4.7% |
| Original dense D=4 | 48.89 / 27.58 | 12.3% / 6.7% |
| Equal-scale D=5 row fixture | 156.83 / 49.85 | 16.0% / 4.9% |
| Row spread 8 bits, D=5 | 91.89 / 52.25 | 8.8% / 5.0% |
| Equal-scale D=4 column fixture | 54.62 / 68.12 | 13.8% / 17.5% |

The zero-spread D=4 row and column fixtures contain the same matrix, yet their
timing ratios differ. This cautions against universal fine-grained estimates
from these benchmark contexts. The ordinary D=3 determinant changes from 2.8%
slower to 2.0% faster between passes. D=2 strict and rounded determinant
conversions add 4.77 ns and 3.30 ns in the reverse pass (6.8% and 4.6%);
their first pass was substantially noisier. The D=2 random-corpus determinant
adds about 7.63 ns per input in the reverse pass (11.9% per 50-input batch).

The first-pass Hilbert determinant improvements do not repeat: D=4 changes
from 26% less time to 4.8% more, and D=5 changes from 32% less time to essentially
unchanged. No Hilbert improvement is claimed. Sign-only results also need
separate interpretation because D≤4 may return from the unchanged floating-point
filter before entering the optimized integer path.

## Control rechecks

Seventeen unchanged controls initially showed more than 10% additional time
and more than 5 ns of added cost, including apparent solve slowdowns of up to
76%. The candidate samples were often much more variable. The complete selected
set was rerun in baseline–candidate–candidate–baseline order with the same
measurement settings; the selection is retained in `control-filter.txt` inside
the evidence archive.

The large slowdowns did not reproduce. Forward recheck ratios range from
0.946 to 1.016; reverse ratios range from 0.946 to 1.043. The D=2 strict solve
conversion changes from an initial ratio of 1.762 to 1.014 and 1.005. The
[control results](exact-row-scaling/controls.json) retain every original and
recheck value. These observations support treating the initial large changes
as unconfirmed, rather than claiming a large solver regression or a precise
zero effect. Smaller or platform-specific effects remain possible.

## Workload tradeoffs

There is no evidence here that simple inputs are more prevalent than difficult
ones. Acceptance rests on the demonstrated benefit region and absolute costs,
not an assumed frequency. A workload dominated entirely by simple exact-value
calls can become slower; near-singular and row-separated inputs can save much
more time per call than simple cases add.

For an illustrative two-case mix, let `c` be the extra time for a simple case
and `s` the saved time for a beneficial case. The estimated break-even share
of beneficial calls is `c / (c + s)`. These are linear cost models, not measured
mixed workloads or frequency estimates. The [mixture calculations](exact-row-scaling/mixtures.json)
retain both passes and their inputs.

| Simple case → beneficial case | Break-even share, first / reverse |
| --- | ---: |
| Dense D=4 → near-singular D=4 | 7.1% / 3.7% |
| Dense D=4 → original mixed-exponent D=4 | 3.7% / 1.9% |
| Equal-scale D=2 → 64-bit row spread D=2 | 27.0% / 17.5% |
| Equal-scale D=3 → 64-bit row spread D=3 | 14.8% / 8.2% |
| Equal-scale D=4 → 64-bit row spread D=4 | 11.3% / 4.3% |
| Equal-scale D=5 → 64-bit row spread D=5 | 14.9% / 5.3% |

These conditional tradeoffs support retaining the optimization under the revised
acceptance criterion. They do not justify a single break-even percentage for
all callers or a claim that every small regression is merely noise.

## Provenance and replay

The baseline is `1f8492786e6a047ba7263788571f42409bf0c2b1`, plus the common
benchmark/test changes in `harness.patch`. The candidate additionally applies
`prototype.patch`. The [evidence archive](exact-row-scaling/evidence.tar.gz)
contains those patches, source and binary hashes, toolchain/host records,
build/validation logs, measurement and analysis scripts, and raw Criterion data.
Compiled executables are excluded. The [checksums](exact-row-scaling/SHA256SUMS)
cover the published archive and derived JSON.

Both versions were built with:

```sh
cargo bench --locked --features bench,exact --bench exact --no-run
```

Keep the two resulting executables as `bin/exact-before` and `bin/exact-after`
beside the extracted scripts. `measure.sh` records the full and determinant-only
filters; `control-repeat.sh` records the control sequence. Use fresh result names
and a separate extraction for new timing runs to preserve the historical samples.
The archive's `candidate.patch` includes the complete pre-measurement source diff
for inspection, while `harness.patch` plus `prototype.patch` separates the builds.

To regenerate every published JSON file from retained samples, from the repository
root after extracting the archive into `target/row-scaling-replay`:

```sh
study=target/row-scaling-replay
published=docs/archive/performance/studies/exact-row-scaling
sh "$study/analyze.sh"
cmp "$published/comparison.json" "$study/comparison.json"
cmp "$published/controls.json" "$study/controls.json"
cmp "$published/mixtures.json" "$study/mixtures.json"
```

All 1,048 named measurements retain 30 finite, positive samples each. Inventory,
sample, source-hash, and binary-hash checks passed. Later archival documentation
changes received focused Markdown, spelling, link, and JSON checks; they do not
change the measured implementation or harness.
