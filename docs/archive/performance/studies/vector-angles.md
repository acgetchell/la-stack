# Unsigned Vector Angle Study (#249)

## Contents

- [Decision](#decision)
- [Environment and method](#environment-and-method)
- [Measurements](#measurements)
  - [Before and after the small-angle correction](#before-and-after-the-small-angle-correction)
- [Evidence and limitations](#evidence-and-limitations)

## Decision

Use power-of-two scaling and compensated exterior minors, followed by
`atan2(exterior_norm, dot)` except for sufficiently small positive angles,
where a bounded direct quotient avoids platform transcendental underflow. The
[mathematical basis](../../../mathematical_basis.md#unsigned-vector-angles)
owns the derivation, range analysis, and rounded contract.

The norm-weighted Kahan formulation from Delaunay was evaluated as the stable
control. Its normalized directions can coincide for inputs whose true
separation is representable: with `n = 2^51`, `[n,n-1]` and `[n+1,n]`
have angle `atan(2^-103)`. Compensated minors preserve this cancellation
residual and the least positive subnormal angle. The quadratic work is an
explicit tradeoff for the crate's small-dimension scope; storage stays constant
and no optional arithmetic dependency is required.

This records implementation evidence for #249 against the completed #251
baseline. Publication of v0.4.7 and downstream adoption remain separate.
The package version remains 0.4.6 in this change.

## Environment and method

Measurements were taken October 8, 2026, on an Apple M4 Max running macOS
27.0.1 (26A434), target `aarch64-apple-darwin`. The compiler was
`rustc 1.99.0 (b940084d7 2026-09-28)`, LLVM 23.1.1.
Cargo's release profile uses fat LTO and one codegen unit; no `RUSTFLAGS`,
`CARGO_ENCODED_RUSTFLAGS`, or `CARGO_BUILD_TARGET` override was set.

The starting la-stack revision was
`b121ea33884ce73d3c5f066c859e3d75675c023f`. The local Delaunay reference was
`55e2e260e62598cf41ff5ae75bc57d798ee592a9`, whose
`src/topology/spaces/spherical.rs` contains the corrected norm-weighted
distance calculation. Delaunay was inspected without modification.

After Windows CI exposed subnormal failures in the initial implementation
`09020aa18a2b233831a8975b343095d2de1670eb`, the corrected kernel was measured
in an isolated checkout with the original locked dependencies. This avoids
concurrent dependency edits in the development checkout. The Criterion 0.8.2
invocation repeated the original command:

```bash
cargo bench --locked --features bench --bench angle -- \
  --sample-size 50 --warm-up-time 1 --measurement-time 2 \
  --noplot --save-baseline angle-249-final
```

All 80 measurements used 50 samples, a one-second warmup, and a two-second
measurement target. No other compilation, validation, or benchmark job from
this task ran during timing. An earlier exploratory run was interrupted when the
normalization counterexample was identified; it is excluded from these results.

The benchmark uses the same fixture arrays for four paths:

- `stable_control`: ordinary max-scaled, norm-weighted Kahan kernel, without
  input validation. Its domain is restricted to independently checked fixtures;
  it is not a general-purpose replacement for the new API.
- `prepared_vector`: already-validated `Vector::angle`; construction excluded.
- `borrowed_slice`: `angle_between`, including shape and coordinate checks.
- `construct_then_angle`: fixed-array construction/validation and the method.

All rows use Criterion `iter` with borrowed inputs through `black_box`.
Fixture construction and independent analytical/integer correctness checks run
before timing. Dense fixtures use integer dot and exterior-norm references;
endpoint cases use `atan(2^-30)` and `π - atan(2^-30)`, mixed-scale inputs
use `π/4`, and the subnormal fixture's analytical answer rounds to
`f64::from_bits(1)`. The inaccurate acos path is not timed.

## Measurements

The table reports arithmetic mean point estimates in nanoseconds. The
[complete retained CSV](vector-angles.csv) includes all 80 mean estimates,
95% marginal bootstrap intervals, and medians.

| Ambient length | Fixture | Stable control (ns) | Prepared Vector (ns) | Borrowed slice (ns) | Construct + angle (ns) |
|---|---|---:|---:|---:|---:|
| 3 | dense | 26.27 | 22.84 | 24.04 | 23.29 |
| 3 | near parallel | 21.53 | 13.64 | 14.43 | 13.67 |
| 3 | near antipodal | 21.56 | 19.45 | 21.02 | 19.36 |
| 3 | mixed scale | 16.76 | 14.43 | 15.30 | 14.61 |
| 3 | subnormal | 19.97 | 13.71 | 14.74 | 13.66 |
| 4 | dense | 39.69 | 33.43 | 35.77 | 34.13 |
| 4 | near parallel | 27.81 | 22.94 | 25.07 | 23.36 |
| 4 | near antipodal | 27.68 | 28.44 | 30.20 | 28.94 |
| 4 | mixed scale | 22.63 | 24.59 | 26.47 | 25.04 |
| 4 | subnormal | 25.69 | 23.05 | 25.21 | 23.45 |
| 5 | dense | 52.06 | 47.91 | 50.21 | 47.71 |
| 5 | near parallel | 30.93 | 34.79 | 37.14 | 34.81 |
| 5 | near antipodal | 31.16 | 37.88 | 40.25 | 38.64 |
| 5 | mixed scale | 26.10 | 34.87 | 36.88 | 34.96 |
| 5 | subnormal | 28.92 | 34.17 | 36.55 | 34.48 |
| 6 | dense | 70.78 | 61.33 | 64.72 | 61.72 |
| 6 | near parallel | 34.39 | 45.91 | 48.67 | 46.42 |
| 6 | near antipodal | 34.50 | 49.02 | 51.43 | 49.48 |
| 6 | mixed scale | 29.64 | 46.39 | 49.44 | 46.70 |
| 6 | subnormal | 32.81 | 45.70 | 48.88 | 46.31 |

The prepared dense path averaged 22.8–61.3 ns across lengths 3–6. Sparse
length-six cases cost more than the linear Kahan control, consistent with
enumerating all 15 exterior minors. Borrowed-slice and construction-inclusive
rows expose caller costs separately. These estimates are descriptive:
subtracting them does not isolate a guaranteed adapter cost, and marginal
intervals do not provide a paired confidence interval for a performance change.

### Before and after the small-angle correction

The [initial measurements at `09020aa`](https://github.com/acgetchell/la-stack/blob/09020aa18a2b233831a8975b343095d2de1670eb/docs/archive/performance/studies/vector-angles.csv)
used the same machine, compiler, dependencies, harness, fixtures, command, and
sampling settings. The prepared-vector means below show representative costs
before and after the correction. Percentage changes describe point estimates;
they do not isolate the code change from run-to-run variation or establish a
paired statistical speedup claim.

| Ambient length | Fixture | Before (ns) | After (ns) | Change |
|---|---|---:|---:|---:|
| 3 | dense | 22.56 | 22.84 | +1.2% |
| 3 | near parallel | 18.53 | 13.64 | -26.4% |
| 3 | subnormal | 15.81 | 13.71 | -13.3% |
| 4 | dense | 32.57 | 33.43 | +2.6% |
| 4 | near parallel | 27.14 | 22.94 | -15.5% |
| 4 | subnormal | 24.91 | 23.05 | -7.5% |
| 5 | dense | 45.31 | 47.91 | +5.7% |
| 5 | near parallel | 38.68 | 34.79 | -10.0% |
| 5 | subnormal | 34.42 | 34.17 | -0.7% |
| 6 | dense | 59.52 | 61.33 | +3.0% |
| 6 | near parallel | 47.74 | 45.91 | -3.8% |
| 6 | subnormal | 45.98 | 45.70 | -0.6% |

## Evidence and limitations

Full `just ci` passed before timing (963 Rust and 502 Python tests, default/exact
doctests, Clippy, static checks, benchmark compilation, and examples). Tests include
small-angle transition and normal/subnormal boundary cases. Existing strict
subnormal assertions remain unchanged. The focused tests also passed: analytical boundary regressions,
integer Gram-determinant properties, full-binary64-range properties, and all
benchmark input gates. The regression corpus includes lengths 0, 1, 2–6, 8,
and 64, signed zero, invalid operands and lengths, overflowing norms, independent
power-of-two rescaling, normalization-collapse examples, and subnormal angles.
Allocation counting confirms both public paths allocate zero times.

The [provenance sidecar](vector-angles.provenance.json) records source and harness
SHA-256 digests. Raw Criterion estimates and 50-sample records remain locally
under `/private/tmp/la-stack-274-validation/target/criterion/angle_d*/<operation>/angle-249-final/`;
the retained CSV
does not include raw samples. Every final sample file was checked for 50 finite,
positive iteration counts and times.

This is one same-machine, same-toolchain comparison of equivalent results on
vetted inputs. It is not a cross-release comparison, a compiler speedup claim,
or a benchmark of Delaunay's radius validation and arc-length conversion.
No certified error bound or correct-rounding claim is inferred from the tests.
