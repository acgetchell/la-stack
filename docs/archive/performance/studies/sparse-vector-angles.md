# Sparse Vector Angle Costs

This study addresses [#287](https://github.com/acgetchell/la-stack/issues/287)
on 2026-10-18. It follows the [angle formulation study](vector-angles.md).
The numerical contract and derivation remain in
[Mathematical basis](../../../mathematical_basis.md#unsigned-vector-angles).

## Contents

- [Scope](#scope)
- [Implementation](#implementation)
- [Measurements](#measurements)
- [Validation](#validation)
- [Provenance](#provenance)

## Scope

The motivating inputs are lower-dimensional embeddings and axis-aligned
directions, where both operands have exact zeros on many of the same axes.
No application-wide frequency or workload-weighted speedup has been measured.
The dense controls represent inputs without shared zeros; these results should
not be generalized to all point clouds or platforms.

The common harness has 14 deterministic families in ambient dimensions 3–6:
dense sign changes, varied dense integer coordinates, dense identical/opposite
directions, adjacent and alternating sparse support, sparse identical/opposite
directions, orthogonal subnormal axes, maximum-magnitude small angles, nearly
parallel/antipodal directions, mixed scales, and subnormal angles. Support sizes
and positions vary. Analytical angles and integer Gram/exterior-product
references are checked outside timing, including the ordinary Kahan control.
The control is valid on these fixtures but does not carry the public API's
general numerical contract.

## Implementation

After complete validation, an axis with both raw coordinates equal to zero
contributes zero to the dot and every incident minor. Signed zeros satisfy this
test. Removing those terms preserves the order of all remaining arithmetic,
including power-of-two scaling, compensated products, compensated dot summation,
`hypot` reduction, and the small-angle quotient. A six-coordinate planar input
with two active axes needs one compensated minor instead of fifteen.

A once-per-call scan selects a const-generic reduction helper. The dense
specialization has no shared-zero checks in its inner loop. The sparse
specialization skips zero axes in both loops. Both use the same arithmetic
source, constant auxiliary space, and arbitrary borrowed-slice lengths.
Coordinate validation and typed error precedence remain unchanged.

The first prototype put the checks directly in both loops. Exploratory quick
measurements showed a dense regression, so it was rejected. These preliminary
runs overlapped tool installation and are diagnostic evidence only. The final
comparison uses the separate dense specialization. Scaling caches, endpoint
shortcuts, and changes to compensation or angle formulas were unnecessary.

A three-second macOS `sample` profile of the original six-coordinate nearly
parallel slice call placed 890 of 2,076 samples (43%) in `hypot` or its dynamic
call stub. Another 1,001 samples were in the shared angle kernel. This locates
the repeated reduction work; it is not a cycle-accurate attribution to individual
multiplications. Profiling was separate from the latency measurements.

## Measurements

Both revisions use Rust 1.99.1, the same expanded harness, locked dependencies,
`bench` feature, release profile with fat LTO and one codegen unit, and borrowed
preconstructed inputs. Binaries were built before timing with:

```bash
cargo +1.99.1 bench --locked --features bench --bench angle --no-run
```

The retained measurement script runs each binary with identical Criterion
arguments: 50 samples, 0.5-second warmup, 1.5-second measurement target,
95% confidence, and no plots. It selects all 56 borrowed-slice rows plus the
14 prepared-vector and 14 unchanged-control rows at ambient length six.
The reverse-order confirmation selects 47 rows: every dense slice row,
six-coordinate dense prepared vectors, all six-coordinate slice/control rows,
and nearly parallel slices at dimensions 3–5.

No compilation or validation from this task overlaps these measurement passes.
Separate Delaunay validation/build activity was observed during both passes;
the host was not exclusively reserved. Unchanged controls expose some temporal
variation, but are not used to normalize away candidate costs. Point-estimate
ratios are descriptive; separate estimate intervals are not a paired confidence
interval for the change.

The [full comparison](sparse-vector-angles/measurements.csv) and
[reverse-order confirmation](sparse-vector-angles/confirmation.csv) retain
median estimates and their 95% intervals. Six-coordinate borrowed-slice results:

| Fixture | Before (ns) | After (ns) | Change | Reverse-order change |
| --- | ---: | ---: | ---: | ---: |
| Dense | 65.83 | 65.62 | -0.3% | -2.7% |
| Dense varied | 98.36 | 102.03 | +3.7% | +1.2% |
| Dense identical | 48.00 | 49.95 | +4.1% | +2.4% |
| Dense antipodal | 49.09 | 51.25 | +4.4% | +2.4% |
| Sparse adjacent | 49.76 | 19.76 | -60.3% | -61.1% |
| Sparse alternating | 53.88 | 30.82 | -42.8% | -43.0% |
| Sparse identical | 49.80 | 13.11 | -73.7% | -73.4% |
| Sparse antipodal | 50.84 | 14.17 | -72.1% | -77.3% |
| Subnormal axes | 50.79 | 18.93 | -62.7% | -64.8% |
| Maximum small angle | 49.73 | 15.64 | -68.5% | -68.4% |
| Nearly parallel | 49.64 | 15.32 | -69.1% | -68.3% |
| Nearly antipodal | 52.57 | 21.73 | -58.7% | -67.2% |
| Mixed scale | 49.65 | 16.87 | -66.0% | -68.4% |
| Subnormal angle | 49.67 | 15.30 | -69.2% | -68.4% |

The nearly parallel slice gains increase with ambient length: the reverse-order
comparison measures reductions of 24.4%, 50.2%, 64.8%, and 68.3% at lengths
3, 4, 5, and 6. All dense slice rows have at most a 4.4% increase in the first
pass and 2.4% in the repeat. Dense prepared-vector increases peak at 5.8% and
2.8%, respectively. These are observed ranges, not regression bounds.

Unchanged controls range from -2.8% to +5.6% in the first pass and -19.9% to
+2.7% in the reverse pass. In particular, the reverse baseline's mixed-scale
and near-antipodal controls slowed substantially; do not interpret their larger
candidate gains as additional optimization. The nearly parallel control moves
by only -1.5% and +0.3% in the two passes, while its public slice call improves
by about 68–69% in both.

**Decision:** retain the shared-zero specialization. The sparse gain is large
and repeatable; the dense dispatch cost is small in these measurements, though
zero overhead is not promised. The evidence supports this local tradeoff, not a
workload-wide speedup, a precise dense-regression bound, or cross-platform timing
claims. The [rejected prototype](sparse-vector-angles/inner-branches.csv) raised
the two exploratory dense medians by roughly 22% and 8%; its inline checks are
not retained.

## Validation

Focused tests cover analytical angles, exact-rational Gram oracles, cancellation,
near-antipodal directions, signed zeros, and representable subnormal angles.
Random signed-zero padding must preserve the result bits. Every two-axis support
placement at length six exercises cancellation and subnormal references.
A 129-coordinate borrowed slice verifies zero allocations and late non-finite
input diagnostics; existing fixed-vector allocation checks remain in place.

`just check` and the final `just ci` passed on Rust 1.99.1. The full local run
passed 1,004 Rust tests, 181 Python tests, 80 default-feature doctests, and 104
exact-feature doctests, plus example execution and benchmark compilation.
Each doctest configuration also reports one ignored `Matrix::det_errbound` example.
This is local macOS validation; Linux and Windows execution remains with CI.

## Provenance

The baseline library is commit
`b2ad940bb49c77f2d56b9f8c2ff8a58cd0cdfeac`, compiled with the Rust 1.99.1 pins
and expanded harness used by the candidate. This compares the library change
on one compiler; it does not estimate the effect of upgrading Rust.
The host is an Apple M4 Max running macOS 27.0.1 (`26A434`),
`aarch64-apple-darwin`, with LLVM 23.1.1.
Rust's full revision is `58cf8f9c35b355199727627d07ae34a4e2846d67`.

The [evidence archive](sparse-vector-angles/evidence.tar.gz) contains raw Criterion
records, timing logs, scripts, source snapshots, the implementation/test patch,
and SHA-256 source/binary manifests. It contains no compiled binaries. Build
each variant with the command above and preserve its executable as `before` or
`after` under `target/angle-sparse-287/` before running the retained scripts.
The archived `after-profile.txt` belongs to the rejected inner-branch prototype,
not the final specialization. All 262 formal measurement records contain 50
finite, positive iteration counts and elapsed times; the exploratory quick
records are separate. The [provenance manifest](sparse-vector-angles/provenance.json)
records artifact hashes and measurement settings.
