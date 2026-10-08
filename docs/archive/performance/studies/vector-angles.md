# Unsigned Vector Angle Study (#249)

## Contents

- [Decision](#decision)
- [Environment and method](#environment-and-method)
- [Measurements](#measurements)
- [Evidence and limitations](#evidence-and-limitations)

## Decision

Use power-of-two scaling and compensated exterior minors, followed by
`atan2(exterior_norm, dot)` except for sufficiently small positive angles,
where a bounded direct quotient avoids platform transcendental underflow. The
[mathematical basis](../../../mathematical_basis.md#unsigned-vector-angles)
owns the derivation, range analysis, and rounded contract.

The public borrowed-coordinate API is the `VectorAngle` extension trait on
`[f64]`: import it directly or from the prelude, then call `left.angle(right)`.
`Vector<D>::angle` uses the same kernel with finite fixed-size storage and
avoids repeating coordinate validation.

The norm-weighted Kahan formulation from Delaunay is the stable control.
Its normalized directions can coincide for inputs whose true separation is
representable: with `n = 2^51`, `[n,n-1]` and `[n+1,n]` have angle
`atan(2^-103)`. Compensated minors preserve this cancellation residual and the
least positive subnormal angle. The quadratic work is an explicit tradeoff for
the crate's small-dimension scope; storage stays constant and no optional
arithmetic dependency is required.

This records implementation evidence for #249 against the completed #251
baseline. Publication of v0.4.7 and downstream adoption remain separate.
The package version remains 0.4.6 in this change.

## Environment and method

Measurements were taken October 8, 2026, on an Apple M4 Max running macOS
27.0.1 (26A434), target `aarch64-apple-darwin`. The compiler was
`rustc 1.99.0 (b940084d7 2026-09-28)`, LLVM 23.1.1.
Cargo's release profile uses fat LTO and one codegen unit; no `RUSTFLAGS`,
`CARGO_ENCODED_RUSTFLAGS`, or `CARGO_BUILD_TARGET` override was set.

The measured working tree starts from
`7ee4bdee0855aa7381ba6ace17f847feb37fa3ce` and includes the extension-trait
API and current locked dependencies. The
[provenance sidecar](vector-angles.provenance.json) identifies that exact source,
harness, and dependency state by SHA-256 digests, verified unchanged before and
after timing. The Delaunay formulation reference is commit
`55e2e260e62598cf41ff5ae75bc57d798ee592a9`,
`src/topology/spaces/spherical.rs`.

The pinned repository toolchain runner invoked Criterion 0.8.2 with:

```bash
cargo bench --locked --features bench --bench angle -- \
  --sample-size 50 --warm-up-time 1 --measurement-time 2 \
  --noplot --save-baseline angle-249-extension
```

All 80 measurements used 50 samples, a one-second warmup, and a two-second
measurement target. No other compilation, validation, or benchmark job from
this task ran during timing.

The benchmark uses the same fixture arrays for four paths:

- `stable_control`: ordinary max-scaled, norm-weighted Kahan kernel, without
  input validation. Its domain is restricted to independently checked fixtures;
  it is not a general-purpose replacement for the new API.
- `prepared_vector`: already-validated `Vector::angle`; construction excluded.
- `borrowed_slice`: `VectorAngle::angle` on slices, including shape and
  coordinate checks.
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
95% marginal bootstrap intervals, and medians from this run.

| Ambient length | Fixture | Stable control (ns) | Prepared Vector (ns) | Slice extension (ns) | Construct + angle (ns) |
|---|---|---:|---:|---:|---:|
| 3 | dense | 26.87 | 23.27 | 24.24 | 23.15 |
| 3 | near parallel | 21.49 | 13.46 | 14.25 | 13.44 |
| 3 | near antipodal | 22.20 | 19.30 | 20.81 | 19.31 |
| 3 | mixed scale | 16.86 | 14.51 | 15.36 | 14.50 |
| 3 | subnormal | 19.68 | 13.42 | 14.31 | 13.57 |
| 4 | dense | 40.30 | 33.58 | 35.03 | 33.32 |
| 4 | near parallel | 27.14 | 22.48 | 24.67 | 22.56 |
| 4 | near antipodal | 27.04 | 27.64 | 29.26 | 28.05 |
| 4 | mixed scale | 21.85 | 23.84 | 25.80 | 24.37 |
| 4 | subnormal | 25.00 | 22.40 | 24.61 | 22.38 |
| 5 | dense | 50.41 | 46.23 | 48.39 | 46.20 |
| 5 | near parallel | 29.98 | 33.44 | 35.54 | 33.74 |
| 5 | near antipodal | 30.14 | 36.71 | 38.80 | 37.12 |
| 5 | mixed scale | 24.95 | 34.30 | 36.62 | 34.06 |
| 5 | subnormal | 28.35 | 33.89 | 35.65 | 33.58 |
| 6 | dense | 68.95 | 60.34 | 62.93 | 61.44 |
| 6 | near parallel | 33.89 | 45.12 | 48.28 | 45.74 |
| 6 | near antipodal | 34.08 | 50.49 | 50.88 | 48.78 |
| 6 | mixed scale | 29.44 | 45.54 | 48.55 | 45.99 |
| 6 | subnormal | 32.06 | 45.16 | 48.11 | 45.57 |

The prepared dense path averaged 23.3–60.3 ns across lengths 3–6. Sparse
length-six cases cost more than the linear Kahan control, consistent with
enumerating all 15 exterior minors. Slice-extension and construction-inclusive
rows expose caller costs separately. These estimates are descriptive:
subtracting them does not isolate a guaranteed adapter cost, and marginal
intervals do not provide a paired confidence interval for a performance change.

## Evidence and limitations

Full `just ci` passed before timing: 970 Rust tests, 502 Python tests, 80 default
and 104 exact-feature runnable doctests, Clippy, static checks, benchmark
compilation, and examples. Each doctest configuration also has one intentional
ignored example. The numerical and error tests call the slice extension method
and compare it with the fixed-vector entry point. A borrowed-subslice test
checks both prelude method resolution and the crate-root trait export.

The regression corpus includes lengths 0, 1, 2–6, 8, and 64, signed zero,
invalid operands and lengths, overflowing norms, independent power-of-two
rescaling, normalization-collapse examples, small-angle transition and
normal/subnormal boundaries, and exact rational Gram references. Allocation
counting covers both the direct-quotient and general `atan2` paths through both
public entry points.

Raw Criterion estimates and sample records remain locally under
`target/criterion/angle_d*/<operation>/angle-249-extension/`.
Every sample file was checked for 50 finite, positive iteration counts and
times; every expected dimension/fixture/operation combination is present.
The retained CSV contains summary estimates, not raw samples.

This is one same-machine, same-toolchain comparison of equivalent results on
vetted inputs. It is not a cross-release comparison, a compiler speedup claim,
or a benchmark of Delaunay's radius validation and arc-length conversion.
No certified error bound or correct-rounding claim is inferred from the tests.
