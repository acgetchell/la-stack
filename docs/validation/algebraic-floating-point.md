# Algebraic Floating-Point Decision

**Decision, 2026-10-15:** reject algebraic floating-point operations and retain
strict arithmetic. The maintainer ended the investigation and requested closure
of [#250](https://github.com/acgetchell/la-stack/issues/250) with this negative
result. The [closing decision](https://github.com/acgetchell/la-stack/issues/250#issuecomment-6196794247)
records the evidence and scope disposition. No algebraic API, feature, or
adoption follow-up is planned.

## Contents

- [Reason for rejection](#reason-for-rejection)
- [Recorded validation](#recorded-validation)
- [Evidence and scope closure](#evidence-and-scope-closure)

## Reason for rejection

The isolated `Matrix::norm_inf` prototype changed positive row additions to
`algebraic_add`. For a D=4 row `[2^969, 2^969, f64::MAX, 0]`, strict accumulation
reports `LaError::NonFinite`, with operation `MatrixInfinityNorm` and location
`MatrixCell { row: 0, col: 2 }`. The algebraic version instead returns the finite
value `f64::MAX`. Appending a zero reproduces the failure in D=5. A final
finite-result check cannot catch an overflow hidden by reassociation.

The prototype reduced measured runtime about 22–25% at D=4/5 and 39–40% at D=8,
with both measurement orders confirming the direction. D=2/3 regressed. Those
timings describe an ineligible prototype on independently checked fixtures;
they do not establish a correctness-preserving optimization. The
[study](../archive/performance/studies/algebraic-float.md) retains absolute
timings, uncertainty, the mathematical counterexample, and compiler analysis.

Rust's [algebraic-operator semantics][algebraic-contract] also provide no
portable accuracy or repeatability guarantee for these operations. An opt-in
API would not repair the unsupported scientific contract. Existing ordinary
operators, deliberate single-rounding `mul_add`, certified bounds, and exact
fallbacks remain the supported arithmetic model; see the
[arithmetic contracts](README.md#arithmetic-contracts-and-evidence).

[algebraic-contract]: https://doc.rust-lang.org/std/primitive.f32.html#algebraic-operators

## Recorded validation

The la-stack source was
`d1f85bb606735d5213a2fa5bd2692e5685a9adc4`, including PR #278. Native execution
used Rust 1.99.0 / LLVM 23.1.1 on an Apple M4 Max with macOS 27.0.1,
`aarch64-apple-darwin`. Measurements used optimized release builds with fat LTO
and one codegen unit. No library behavior changed.

| Evidence | Observed outcome | Scope |
|----------|------------------|-------|
| Independent exact rational summation | Explicit D=4/5 overflow classification failures; finite-value differences also retained | Prototype compared with the strict method and exact binary64-value reference |
| Adversarial and random inputs | Signed zero, subnormals, mixed/extreme magnitudes, range boundaries, input rejection, and 1,000 random matrices per dimension | D=0, 1, 2, 3, 4, 5, 8, 16; cancellation is inapplicable to absolute row sums |
| Upstream norm checks | 18 unit tests and 5 matrix property tests passed | Unchanged strict implementation; 256 property cases, seed 250 |
| Downstream Delaunay checks | 728 geometry unit tests, 26 geometry properties, and 18 predicate properties passed | Unchanged strict la-stack patched into Delaunay commit `bde38f4826cc248b425c4300cadd233417f7646c`; prototype not installed |
| Criterion measurements | 40 saved measurements, each with 100 samples; all ten dimension/fixture pairs repeated in reverse order | D=2–5 and D=8; exact fixture checks outside timing |
| Linux and Windows targets | Strict library and prototype probes cross-compiled | No native numerical execution or performance measurement |

These are historical results for the recorded revisions, not a full validation
of every arithmetic API or a cross-platform reproducibility promise.

## Evidence and scope closure

The [full study](../archive/performance/studies/algebraic-float.md),
[measurements](../archive/performance/studies/algebraic-float/measurements.csv),
[provenance](../archive/performance/studies/algebraic-float/provenance.json), and
[evidence archive](../archive/performance/studies/algebraic-float/evidence.tar.gz)
retain the sources, locks, raw samples, numerical outputs, compiler probes, and
test logs. Experimental Rust must be extracted outside repository-owned source.

The original plan called for native Linux and Windows prototype execution.
Those runs were not performed. On 2026-10-15 the maintainer accepted rejection
and discontinued that remaining experimental work: there is no candidate to
ship and no adoption decision awaiting those runs. This is an explicit scope
closure, not a claim that the original native-platform criterion passed.
The archived provenance preserves the earlier execution status; this record
owns the final disposition. Normal cross-platform validation of the strict
library continues independently.
