# Triangular Rational Solve Study (#246)

## Contents

- [Decision](#decision)
- [Environment and method](#environment-and-method)
- [Correctness](#correctness)
- [Measurements](#measurements)
- [Evidence and limitations](#evidence-and-limitations)

## Decision

Recognize row permutations of upper and lower triangular rational matrices and
solve them by direct exact substitution. Recognition maps distinct first or
last non-zero columns to their source rows. It requires no matrix copy and
proves every pivot non-zero before division. Residuals accumulate as exact
integer numerator/denominator pairs, with reduction once per solved component.
General and singular inputs keep
the existing denominator-clearing and Bareiss path. Public construction,
canonicalization, exactness, and error contracts are unchanged.

The [mathematical basis](../../../mathematical_basis.md#exact-arithmetic-over-rational-inputs)
owns the proof and arithmetic contract. The shared integer backend and
determinant paths are unchanged. This optimization is specific to
`RationalMatrix::solve`; it does not change `Matrix::solve_exact`.

## Environment and method

Measured October 9, 2026, on an Apple M4 Max running macOS 27.0.1 (26A434),
target `aarch64-apple-darwin`. Both phases use
`rustc 1.99.0 (b940084d7 2026-09-28)`, LLVM 23.1.1, Criterion 0.8.2,
the same Cargo lockfile, and the repository release profile (fat LTO, one
codegen unit, no debug information). No `RUSTFLAGS`,
`CARGO_ENCODED_RUSTFLAGS`, or `CARGO_BUILD_TARGET` override was set.

Both measurements used `Cargo.lock` from the baseline commit below. The later
lockfile refresh of `cc` and `syn` was not part of these timing runs; reproduce
the recorded measurements with the baseline lockfile and the recorded harness.

The baseline library is commit
`d12bdb499e0a7ee889c2a8cc7ad7871c7084d759`, with the new benchmark harness.
The optimized library is the same checkout with the triangular solve change.
The [provenance record](rational-triangular.provenance.json) records exact source
and harness hashes. The Rust 1.99 baseline reproduced the triangular regression
before implementation began. The Rust 1.98.1 issue measurements are historical
context, not the comparison baseline.

Both phases use this command, replacing `PHASE` with `before` or `final`:

```bash
filter='permuted_upper/rational_d[2-8]|'
filter+='(diagonal|upper|lower|permuted_lower|sparse|dense)/.*_d8|'
filter+='permuted_upper/dyadic_d8'
cargo bench --locked --features bench,exact --bench rational_solve -- \
  "$filter" \
  --sample-size 30 --warm-up-time 1 --measurement-time 3 --noplot \
  --save-baseline "rational-246-$PHASE"
```

These 80 measurements cover the issue's non-dyadic row-permuted upper family
at every size 2–8, plus all seven structural families with both coefficient
domains at size 8. The permanent harness registers all families at every
size 2–8; the narrower timing selection keeps this study focused.
Compilation precedes each phase. No other compilation, validation, or timing
job from this task overlaps measurement.

Each fixture first verifies its manufactured solution with the public solver
and the independent captured Gaussian solver. Criterion `iter` closures use
borrowed inputs through `black_box` and include output destruction:

- `adapter`: runtime `try_with_rational_matrix!` dispatch, zero matrix filled
  through checked setters, RHS construction, solve, and collection into a vector.
- `construct`: the same setter-based matrix and RHS construction, without a solve.
- `legacy`: the captured Delaunay Gaussian solver, including owned input copies.
- `prepared`: public `RationalMatrix::solve` on preconstructed validated inputs.

The [captured reference](../../../../benches/common/rational_solve.rs) retains
the arithmetic, cloning, and zero skips from Delaunay commit
`2b7017bb3981babed2aab64826282bb06ac0e312`, `src/geometry/matrix.rs`, lines 215–256.
It excludes the removed rational-to-binary64 probe, as required by the issue.
All coefficients and manufactured right-hand sides use exact rational arithmetic.

## Correctness

The [structured regressions](../../../../tests/rational_structured_solve.rs)
cover dimensions 0–8, arbitrary row permutations, upper and lower triangles,
mixed signs and denominators, zero RHS, and rational components wider than
256 bits. They require exact `A x == b`, the manufactured solution, and
agreement with independent Gaussian elimination. Dense controls force the
general fallback. Zero pivots, zero rows, nonzero rank-deficient rows, and
duplicate rows preserve the exact singularity reason and pivot column,
including after row permutations. Existing construction tests retain invalid
raw-denominator rejection and setter failure atomicity.

## Measurements

All values below are arithmetic means in microseconds. Brackets show marginal
95% Criterion bootstrap intervals. The [complete CSV](rational-triangular.csv)
retains mean, median, and slope estimates with their intervals for all 160
measurements; [raw samples and estimates](rational-triangular.measurements.json)
retain the underlying Criterion evidence.

For the reproduced row-permuted upper family:

| D | Legacy (final phase) | Prepared before | Prepared final | Adapter before | Adapter final |
|---|---:|---:|---:|---:|---:|
| 2 | 0.804 [0.801, 0.807] | 0.879 | 0.233 | 1.362 | 0.713 [0.711, 0.714] |
| 3 | 1.845 [1.838, 1.852] | 1.930 | 0.403 | 2.804 | 1.229 [1.227, 1.231] |
| 4 | 3.413 [3.397, 3.428] | 3.637 | 0.661 | 4.888 | 1.878 [1.873, 1.883] |
| 5 | 5.979 [5.871, 6.085] | 6.342 | 0.927 | 8.018 | 2.633 [2.628, 2.637] |
| 6 | 8.126 [8.098, 8.154] | 9.609 | 1.200 | 12.103 | 3.724 [3.680, 3.768] |
| 7 | 11.253 [11.215, 11.291] | 15.352 | 1.533 | 18.277 | 4.687 [4.676, 4.698] |
| 8 | 14.922 [14.856, 14.981] | 22.022 | 1.975 | 25.706 | 6.041 [6.018, 6.064] |

The final full adapter's upper interval endpoint is below the legacy control's
lower endpoint at every reproduced dimension. At D=8, prepared solve time falls
from 22.022 to 1.975 µs (11.15×), and the full adapter from 25.706 to 6.041 µs
(4.26×). The final adapter is 2.47× faster than the same-phase legacy mean.
These are descriptive ratios for this family. Construction is unchanged:
3.766 µs before and 3.770 µs after at D=8.

General-system controls retain the Bareiss path:

| D=8 control | Prepared before | Prepared final | Adapter before | Adapter final |
|---|---:|---:|---:|---:|
| dense, rational | 32.202 | 32.201 | 37.862 | 37.770 |
| dense, dyadic | 27.052 | 27.249 | 31.625 | 32.186 |
| sparse, rational | 22.024 | 21.330 | 25.507 | 25.085 |
| sparse, dyadic | 20.155 | 19.878 | 22.606 | 22.671 |

The dense Hilbert kernel is effectively unchanged. The dyadic dense kernel
increases 0.73% and its full adapter 1.77%; no zero-cost detection claim is made.
Sparse kernel means decrease 1.38–3.15%, while sparse adapter means move between
−1.66% and +0.28%. These controls show no material loss of the dense-system
advantage in the measured cases. All other triangular families, including lower
and row-permuted lower systems, retain their full before/after measurements in
the CSV.

## Evidence and limitations

These synthetic systems establish a focused solver improvement, not a measured
distribution of Delaunay LP bases or a whole-triangulation speedup. General
sparse inputs are controls; this change does not add a general sparse solver.
Arithmetic operation counts exclude arbitrary-precision bit complexity.

Point-estimate ratios are descriptive. Marginal Criterion confidence intervals
are not paired confidence intervals for a change. Construction, solve, and
adapter rows are independent measurements with their own destruction costs;
their times need not add. Repeated measurements of the unchanged legacy and
construction controls help expose environmental drift.
