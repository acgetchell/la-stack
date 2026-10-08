# Unsigned Vector Angle Study (#249)

## Contents

- [Decision](#decision)
- [Environment and method](#environment-and-method)
- [Measurements](#measurements)
- [Evidence and limitations](#evidence-and-limitations)

## Decision

Use power-of-two scaling and compensated exterior minors, followed by
`atan2(exterior_norm, dot)`. The
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

The single final Criterion 0.8.2 invocation used this command:

```bash
cargo bench --locked --features bench --bench angle -- \
  --sample-size 50 --warm-up-time 1 --measurement-time 2 \
  --noplot --save-baseline angle-249-final
```

All 80 measurements used 50 samples, a one-second warmup, and a two-second
measurement target. No other agent-owned compilation, validation, or benchmark
job ran during timing. An earlier exploratory run was interrupted when the
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
| 3 | dense | 25.84 | 22.56 | 23.80 | 22.63 |
| 3 | near parallel | 21.48 | 18.53 | 19.99 | 18.54 |
| 3 | near antipodal | 21.69 | 18.80 | 20.23 | 18.80 |
| 3 | mixed scale | 16.79 | 14.15 | 15.16 | 14.16 |
| 3 | subnormal | 19.46 | 15.81 | 17.07 | 15.92 |
| 4 | dense | 39.10 | 32.57 | 34.03 | 32.62 |
| 4 | near parallel | 26.42 | 27.14 | 28.48 | 27.46 |
| 4 | near antipodal | 26.59 | 27.24 | 28.92 | 27.78 |
| 4 | mixed scale | 21.58 | 23.37 | 25.18 | 23.72 |
| 4 | subnormal | 24.44 | 24.91 | 26.54 | 25.17 |
| 5 | dense | 48.96 | 45.31 | 50.12 | 48.63 |
| 5 | near parallel | 32.81 | 38.68 | 41.39 | 40.29 |
| 5 | near antipodal | 30.15 | 36.44 | 38.37 | 36.84 |
| 5 | mixed scale | 24.90 | 33.54 | 35.79 | 33.82 |
| 5 | subnormal | 27.80 | 34.42 | 36.70 | 34.82 |
| 6 | dense | 66.77 | 59.52 | 62.46 | 60.10 |
| 6 | near parallel | 33.26 | 47.74 | 50.37 | 48.25 |
| 6 | near antipodal | 33.72 | 47.85 | 50.31 | 48.41 |
| 6 | mixed scale | 29.04 | 44.93 | 48.04 | 45.56 |
| 6 | subnormal | 31.75 | 45.98 | 52.28 | 49.80 |

The prepared dense path averaged 22.6–59.5 ns across lengths 3–6. Sparse
length-six cases cost more than the linear Kahan control, consistent with
enumerating all 15 exterior minors. Borrowed-slice and construction-inclusive
rows expose caller costs separately. These estimates are descriptive:
subtracting them does not isolate a guaranteed adapter cost, and marginal
intervals do not provide a paired confidence interval for a performance change.

## Evidence and limitations

The focused tests passed before timing: analytical boundary regressions,
integer Gram-determinant properties, full-binary64-range properties, and all
benchmark input gates. The regression corpus includes lengths 0, 1, 2–6, 8,
and 64, signed zero, invalid operands and lengths, overflowing norms, independent
power-of-two rescaling, normalization-collapse examples, and subnormal angles.
Allocation counting confirms both public paths allocate zero times.

The [provenance sidecar](vector-angles.provenance.json) records source and harness
SHA-256 digests. Raw Criterion estimates and 50-sample records remain locally
under `target/criterion/angle_d*/<operation>/angle-249-final/`; the retained CSV
does not include raw samples. Every final sample file was checked for 50 finite,
positive iteration counts and times.

This is one same-machine, same-toolchain comparison of equivalent results on
vetted inputs. It is not a cross-release comparison, a compiler speedup claim,
or a benchmark of Delaunay's radius validation and arc-length conversion.
No certified error bound or correct-rounding claim is inferred from the tests.
