# Interval and Certified-Reduction Optimization

This study records the Rust 1.99 investigation for
[#247](https://github.com/acgetchell/la-stack/issues/247) and
[#248](https://github.com/acgetchell/la-stack/issues/248). It separates library
changes from the downstream adapter needed to remove Delaunay's private kernels.
It supplements the [compiler baseline](rust-1.99.md); it does not replace the
generated [release performance report](../../../performance.md).

## Contents

- [Implementation and contracts](#implementation-and-contracts)
- [Baseline and profiling](#baseline-and-profiling)
- [Measurement provenance](#measurement-provenance)
- [Library measurements](#library-measurements)
- [Downstream measurements](#downstream-measurements)
- [Rejected intermediate candidate](#rejected-intermediate-candidate)
- [Correctness and decision](#correctness-and-decision)

## Implementation and contracts

- Certified reductions retain the ordered FMA estimate tree and product/result
  normality requirements. While proof is available, zero products preserve the
  normal or positive-zero estimate exactly and skip both FMA and magnitude work.
  After proof loss, their FMA still runs to preserve possible negative-zero
  behavior. After the first nonzero product, an upward-widened FMA accumulates
  the magnitude bound. The first
  product steps upward unconditionally unless it has an exact unit multiplier;
  this can slightly widen the magnitude bound compared with residual testing.
  A single normal product with a unit multiplier has a zero error bound, proved
  by exact unit scaling and the remaining zero-product FMAs. This stronger
  proof can return an exact certificate at `f64::MAX`, where a conservative
  positive bound previously made the endpoints unavailable.
- Determinants keep their existing subset expansion and operation order. Const
  dimension dispatch sizes the inline workspace to `max(2, 2^D)` intervals.
- Interval arithmetic reuses singleton squares/sums, recognizes exact
  subtraction by/from zero, and rounds only the needed endpoints for general
  addition. Finite endpoint and typed range-error checks remain in place.
- Product rounding uses an FMA residual only when the rounded product magnitude
  is at least `2^-968`. Below that threshold, the original integer-significand
  comparison remains authoritative, including normal products with a residual
  too small for binary64. Power-of-two scaling in the admitted range is exact.
- Scalar certificates store their estimate and absolute bound. Construction
  proves that both outward endpoints fit the finite range, using the dominating
  magnitude `|estimate| + bound`; accessors derive only the endpoints requested.
  This removes duplicate stored state while retaining infallible const accessors.
- The measured hot kernels and endpoint accessors explicitly inline to avoid
  per-coordinate or per-endpoint calls and large `Result` aggregates. The rare
  integer product comparison stays out of line.

The [mathematical basis](../../../mathematical_basis.md#certified-fixed-vector-reductions)
owns the reduction proof, and its
[interval section](../../../mathematical_basis.md#outward-rounded-interval-expressions)
owns the FMA threshold derivation. Public method signatures, const evaluation, checked
construction, typed errors, and exact-fallback meanings remain unchanged.

## Baseline and profiling

The prerequisite #251 was merged before this work. Fresh Rust 1.99 runs against
Delaunay's private implementation reproduced both regressions before arithmetic
was edited. The unoptimized adapter retained exact fallback on interval zero,
inconclusive signs, errors, and unavailable scalar certificates.

| Workload | Private kernel | Unoptimized public adapter |
|----------|---------------:|---------------------------:|
| In-sphere 2D, 10,000 queries | 1.0794 ms | 2.0188 ms |
| In-sphere 3D, 10,000 queries | 2.1143 ms | 3.3247 ms |
| In-sphere 4D, 10,000 queries | 9.6806 ms | 11.505 ms |
| In-sphere 5D, 10,000 queries | 17.803 ms | 20.147 ms |
| Single shared vertex 2D | 98.596 ns | 141.83 ns |
| Single shared vertex 3D | 163.32 ns | 226.65 ns |
| Single shared vertex 4D | 1.3572 µs | 1.5055 µs |
| Single shared vertex 5D | 2.1929 µs | 2.3664 µs |
| Whole Level 4, 2D / 500 vertices | 6.5003 ms | 7.3381 ms |

These are short diagnostic Criterion central estimates: 20 samples, 0.2-second
warmup, one-second measurement, and 1,000 resamples. They establish the blocker;
the final comparisons below use fresh private baselines and longer measurements.

Five-second sampling profiles located interval product rounding, outward
addition, and workspace clearing, plus certified magnitude accumulation and
certificate endpoint construction. A later profile and symbol inspection showed
per-coordinate `CertifiedReduction::add_product` calls copying reduction state.
Inlining limits attribution, so sample counts are not presented as a partition
of the full regression. Diagnostic logs retain unsuccessful intermediate trials.

## Measurement provenance

Measurements use Rust 1.99.0 (`b940084d7eb6a299eb4bfeb8e34901bc051e7ac4`),
LLVM 23.1.1, an Apple M4 Max, `aarch64-apple-darwin`, and macOS 27.0.1 (26A434).
The library base is `2095d18`; Delaunay's base is `2f2f01451`.
The sibling Delaunay checkout remains unchanged; builds and adapter trials use
an isolated source copy and a Cargo path patch.

Upstream binaries use `bench,exact`, optimization level 3, fat LTO, and one
codegen unit. Both revisions use the same final benchmark sources. Downstream
uses Delaunay's default features and its `exact` la-stack dependency, with the
`perf` profile: optimization level 3, thin LTO, one codegen unit, and line-table
debug information. No relaxed floating-point or target-CPU flags are used.

The [measurement script](interval-reductions/measure.sh) runs prebuilt binaries
serially after correctness validation. The maintained helper now resolves input
directories before changing working directory and saves every invocation's logs
and raw Criterion data under a fresh `<phase>.*` directory beside the binaries.
It sets `CRITERION_HOME` separately for each comparison, preserving `new` and
`change` when later phases run. The original script remains in the unchanged
evidence archive; these reproduction fixes did not produce the historical tables.
Upstream uses 100 samples, a 0.2-second
warmup, a one-second measurement window, and 10,000 resamples. Downstream uses
100 samples, a one-second warmup, three seconds of measurement, and 10,000
resamples. Every downstream case is repeated in reverse order. No build,
validator, or second timing job overlaps these measurements. Other desktop
activity is uncontrolled.

The downstream performance gate uses a 5% slowdown band for materiality and
checks both measurement orders, including whole-validation controls. This is a
local acceptance criterion, not a cross-platform guarantee. All individual
results and Criterion change confidence intervals are retained, including
unfavorable or inconclusive comparisons.

The [provenance record](interval-reductions/provenance.json) includes full
revision and lockfile hashes. The final source/benchmark inventory SHA-256 is
`3d963d404299248fe0edd93bdcae375c105daef5bbe77afbc16064482a815725`.
It remained unchanged through the recorded full CI and final measurement. Both
baseline benchmark files had the same SHA-256 as their candidate counterparts.
Subsequent review added regression tests, helper documentation, and stronger
scalar benchmark setup assertions. The archive retains the measured source and
harness; the arithmetic implementation and timed closures are unchanged.

The [evidence archive](interval-reductions/evidence.tar.gz) retains raw Criterion
samples and confidence intervals, profiles, validation/build logs, source
patches, lockfiles, and the CSV exporter. It includes rejected candidates as
separate datasets. Its SHA-256 is
`aebe0581fb5fdfc331e2f22f843553d65bbc143b3184cd47519999800ab4728e`.
The archive's README identifies baseline labels and export commands; the
exporter checks sample counts, finite positive measurements, and complete 95%
confidence intervals before producing tables.

Final downstream narrow-phase runs use Criterion's `Flat` sampling: every
sample performs the same number of iterations. Apply the retained sampling-mode
patch to both private and adopted source copies before building their
`realization_validation` binaries. Whole-validation cases retain their original
sampling configuration. The investigation below retains the earlier linear and
flat measurements that motivated the final zero-product optimization.

Build the upstream binaries with `cargo bench --locked --features bench,exact
--bench interval --bench linear_form --no-run --target-dir <revision-specific-output>`,
running in each source directory and using the expanded harness for both source
revisions. Build downstream with `cargo bench --profile perf
--bench cold_path_predicates --bench realization_validation --no-run`, adding
`--config 'patch.crates-io.la-stack.path="<library-source>"'` for the selected
library. Preserve separate upstream build directories and save each executable.
A pilot that shared an output directory was excluded after symbol inspection
found the candidate-only integer-fallback symbol in its baseline executable.
Final baseline binaries must contain no `compare_product_with_rounded_integer`
symbol; final candidate binaries retain that fallback. The script accepts the
saved-binary directory, the upstream directory, and the isolated downstream
directory; its optional fourth
argument `repeat` runs the downstream reverse-order comparison; `downstream`
reruns just the downstream suites. The `control` phase runs the extended 3D
whole-validation comparison with 200 samples. The `confirm`
phase repeats unfavorable upstream cases and measures the selected downstream
controls with longer windows. The historical `flat` phase uses separately saved
`realization-flat-private` and `realization-flat-after` binaries, each rebuilt
with the identical [sampling-mode patch](interval-reductions/flat-sampling.patch).
Their only harness change selects Criterion's `Flat` sampling for narrow-phase
cases. Fixture construction and correctness assertions remain outside timing.

The helper prints its output directory. Within it, export the preserved data
using `compare.py` from the evidence archive, running from the repository root:

```sh
uv run --locked <extracted>/compare.py <measure-output>/upstream-data issue-before
uv run --locked <extracted>/compare.py <measure-output>/downstream-data private-final
uv run --locked <extracted>/compare.py <repeat-output>/downstream-repeat-data adopted-repeat reverse
uv run --locked <extracted>/compare.py <control-output>/control-3d-data control-3d reverse 200
```

The `downstream` phase also writes `downstream-data`. For `confirm`, export
`upstream-confirm-data` with `adopted-confirm reverse`, `confirm-forward-data`
with `confirm-forward`, and `confirm-reverse-data` with `confirm-reverse reverse`.
The `flat` phase similarly writes `flat-forward-data` and `flat-reverse-data`
with their matching baseline labels; only the reverse export takes `reverse`.

## Library measurements

The [complete library table](interval-reductions/upstream.csv) contains all 183
cases: 43 interval and 140 vector measurements. Times below are sample means;
percentage changes compare those means. The CSV also retains medians, absolute
mean confidence intervals, and Criterion's bootstrapped change intervals.
These are focused comparisons, not a release-wide speed claim.

| Workload | Before (ns) | After (ns) | Mean change |
|----------|------------:|-----------:|------------:|
| Point interval square | 11.27 | 2.18 | -80.7% |
| 3×3 lifted assembly + sign | 115.15 | 39.90 | -65.3% |
| 7×7 lifted assembly + sign | 1570.60 | 1056.37 | -32.7% |
| 3×3 dense determinant sign | 98.29 | 42.64 | -56.6% |
| 7×7 dense determinant sign | 1932.45 | 1575.33 | -18.5% |
| D=2 dense dot + endpoints | 17.27 | 6.24 | -63.9% |
| D=3 dense dot + endpoints | 24.71 | 11.02 | -55.4% |
| D=6 dense dot + endpoints | 47.60 | 20.12 | -57.7% |
| D=3 sparse dot + endpoints | 17.13 | 2.78 | -83.8% |
| D=3 dense difference + endpoints | 50.09 | 28.35 | -43.4% |
| D=3 origin batch, prepared | 132.29 | 25.24 | -80.9% |
| D=3 origin batch, checked construction | 148.78 | 29.85 | -79.9% |
| D=3 translated batch, prepared | 193.74 | 84.67 | -56.3% |
| D=3 translated batch, checked construction | 205.89 | 88.36 | -57.1% |

Three interval shortcut microbenchmarks regress: zero addition by 0.15 ns,
zero multiplication by 0.92 ns, and multiplication in the cancellation fixture
by 0.73 ns in the final pass. The latter uses a unit interval multiplier.
These paths benefit little from arithmetic optimization; their absolute costs
remain below 2 ns. The tradeoff is explicit: the hot determinant and downstream
workflows improve, while these isolated shortcuts do not. Construction-only
D=2 controls fall by 0.17–0.21 ns even though `Vector::try_new` is unchanged;
earlier runs moved in the other direction. No constructor optimization is
claimed, and isolated timings are not an additive decomposition of downstream
runtime.

## Downstream measurements

The [retained adapter patch](interval-reductions/adoption.patch) keeps geometry, labels, orientation interpretation, and
exact fallback in Delaunay. It constructs checked upstream values and accepts
only certified separation. Beyond the original reuse of the common projection,
the final adapter:

1. Scales the axis into a fixed-size array, letting `Vector::try_new` own finite
   validation instead of allocating and then validating the same values twice.
2. Computes outward thresholds `shared.upper + 1` and `shared.lower - 1` once.
   Each vertex compares its certified estimate strictly against the rounded
   threshold plus its error (or threshold minus error on the negative side).
   For representable `x`, `x > RN(y)` implies `x > y`, because rounding is
   monotone and `RN(x) = x`; the symmetric argument holds for `<`. Thus these
   comparisons prove the same unit margin without constructing each vertex's
   endpoint. Ambiguity and overflow retain exact fallback.
3. Uses direct short-circuiting loops for per-vertex checks. Sampling the
   intermediate adapter found outlined `Iterator::try_fold` bodies in these
   checks; the direct loops expose the same proof and rejection paths together.

These downstream changes are part of the adoption candidate, not library-only
speedups. The upstream batch benchmarks still distinguish checked construction
from prepared arithmetic and measure the original per-vertex subtraction
workflow for both library revisions.

The [first pass](interval-reductions/downstream.csv) and
[reverse-order pass](interval-reductions/downstream-repeat.csv) retain all 15
cases. Times below are first-pass sample means. Negative changes are faster.
The reverse column consistently reports adopted/private, even though the
execution order is reversed.

| Workload | Private | Adopted | Change | Reverse-order change |
|----------|--------:|--------:|-------:|---------------------:|
| In-sphere boundary 2D | 1.438 µs | 1.304 µs | -9.4% | -3.1% |
| In-sphere boundary 3D | 2.566 µs | 2.620 µs | +2.1% | +2.3% |
| In-sphere boundary 4D | 4.943 µs | 4.824 µs | -2.4% | -1.2% |
| In-sphere 2D, 10,000 queries | 1.154 ms | 951.404 µs | -17.6% | -17.8% |
| In-sphere 3D, 10,000 queries | 2.273 ms | 1.797 ms | -21.0% | -20.1% |
| In-sphere 4D, 10,000 queries | 11.145 ms | 10.004 ms | -10.2% | -8.2% |
| In-sphere 5D, 10,000 queries | 18.992 ms | 16.275 ms | -14.3% | -10.1% |
| Single shared vertex 2D | 108.39 ns | 94.26 ns | -13.0% | -8.8% |
| Single shared vertex 3D | 170.49 ns | 165.07 ns | -3.2% | +0.2% |
| Single shared vertex 4D | 1.389 µs | 1.373 µs | -1.2% | -0.4% |
| Single shared vertex 5D | 2.199 µs | 2.212 µs | +0.6% | +0.3% |
| Whole Level 4 2D, 500 vertices | 6.786 ms | 6.769 ms | -0.3% | +1.3% |
| Whole Level 4 3D, 20 vertices | 111.670 ms | 115.509 ms | +3.4% | +4.4% |
| Whole Level 4 4D, 10 vertices | 115.405 ms | 118.668 ms | +2.8% | +3.6% |
| Whole Level 4 5D, 8 vertices | 17.415 ms | 17.725 ms | +1.8% | +3.6% |

All final mean changes stay below the study's 5% slowdown band. The 2D
shared-vertex workload improves by 13.0%/8.8%, and whole 2D validation changes
by −0.3%/+1.3%. The larger whole-validation controls retain small slowdowns;
this is not a claim that adopting the stronger contracts has zero cost.

The reverse whole-3D comparison has a 95% change interval of +3.83% to +5.12%,
which narrowly crosses the band. An
[extended control](interval-reductions/downstream-control.csv) used 200 samples,
a three-second warmup, a twenty-second requested measurement window, and the
same adopted-first order. It measured −1.32% (95% interval −1.56% to −1.09%),
so the material slowdown did not repeat. Keep the original result: desktop
timing variation limits precision, particularly for this control.

## Rejected intermediate candidate

The candidate before the final zero-product change failed its 2D shared-vertex
gate. Its [first pass](interval-reductions/initial-candidate/downstream.csv) and
[reverse pass](interval-reductions/initial-candidate/downstream-repeat.csv)
remain available. Those original linear-sampling
runs have early samples around 135 ns and later samples near 90 ns. The
candidate's first-pass mean is 113.42 ns versus 106.55 ns privately, despite
medians of 93.22 ns versus 102.06 ns. Thus the first mean comparison crosses the
5% band; the faster median alone does not establish acceptance.
[Longer forward](interval-reductions/initial-candidate/downstream-confirm-forward.csv) and
[longer reverse](interval-reductions/initial-candidate/downstream-confirm-reverse.csv) runs use a
three-second warmup and five-second measurement window, retaining 100 samples.
They reproduce this split (+5.5% and +3.9% in mean, with overlapping change
intervals), while whole 2D changes by +1.9%/+2.4% and whole 5D by +2.9%/+1.0%.
The initial candidate also regressed in fixed-work-per-sample controls:
[+17.9% forward](interval-reductions/initial-candidate/downstream-flat-forward.csv)
and [+9.5% reversed](interval-reductions/initial-candidate/downstream-flat-reverse.csv).
That candidate was rejected. The final library omits provably identical
zero-product FMAs while proof is available; it retains their signed-zero
behavior after proof loss. The
[next adapter repeat](interval-reductions/intermediate-adapter/downstream-repeat.csv)
still showed a +7.3% mean slowdown in 2D. Its profile located outlined
per-vertex iterator bodies. The final adapter's direct loops remove those call
layers without changing certificate arithmetic or geometric decisions. Final
results use fresh comparisons of that adapter with the same library code.

The [shortcut repeat](interval-reductions/initial-candidate/upstream-repeat.csv), run in reverse
order, confirms the three scalar interval regressions: approximately
0.20–0.97 ns, with final costs below 2 ns. Construction-only D=2 controls move
by +0.7% and −0.5% in that repeat. No vector-constructor optimization is claimed.

## Correctness and decision

Independent exact-rational tests check tight interval containment, determinant
signs, and scalar certificate containment. Added cases straddle the FMA
residual-underflow threshold and cover both signs, subnormals, maximum finite
products, signed-zero reduction state, and D=6 dot/difference certificates.
Existing tests retain non-finite provenance, range exhaustion, cancellation,
inconclusive intervals, and const-evaluation coverage.

`just check` and `just ci` passed. Full CI ran 904 Rust tests, 499 Python tests,
both doctest configurations, warning-denying Clippy/rustdoc, static-analysis
fixtures, benchmark compilation, and runnable examples. The final downstream
adapter passed Clippy and all 97 selected predicate/intersection tests, including
exact agreement and known shared-face/overlap cases through D=6. The additional
54 upstream exact-arithmetic properties passed with
`PROPTEST_CASES=512` and `PROPTEST_RNG_SEED=248`.

Release and downstream adoption remain separate maintainer actions. A published
version number alone does not establish downstream acceptance; the released
implementation and retained adapter must pass the same checks together.

The local evidence supports the #247/#248 implementation and the retained
adapter under the stated 5% comparison band. It does not establish performance
on other CPUs, compilers, profiles, or fixture distributions. The scalar shortcut
regressions and the whole-3D uncertainty remain visible above; the performance
claim is scoped to the measured predicate and realization workflows.
