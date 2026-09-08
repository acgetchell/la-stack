# Benchmark Performance

**la-stack** v0.4.6 · `d36a9e9` (HEAD)
**Source revision timestamp**: 2026-09-08 01:59:08 UTC (deterministic report metadata; not the benchmark measurement time)
**Benchmark measurement timestamp**: not recorded by Criterion; use the provenance below to identify the measured revisions and environment.
**Statistic**: median
**Suite**: all
**Scope**: release-signal

## Benchmark Results

Comparison against baseline **v0.4.5**:

Negative point-estimate change means the current point estimate is smaller; a baseline/current point-estimate ratio above 1.00 has the same meaning.
The CI-relation column reports only whether the two marginal Criterion intervals overlap. These are not paired confidence intervals
for the change, so the report makes no statistical-significance or performance-improvement claim from interval separation.

### Reproducibility Provenance

**Measurement environment**: recorded for both samples under one shared current harness.

- CPU: `Apple M4 Max (arm64)`
- OS: `Darwin 25.6.0 arm64`
- rustc: `rustc 1.98.1 (48a229cea 2026-09-01)`
- Current commit: `d36a9e99ec223c4cbbf8a9191818f186159c27f7`
- Current Git clean: `false`
- Current source-state SHA-256: `b8fcf706fa6e287ca6682a1ecef5a69576413b9c88040184a9b1ecb8f6e92336`
- Baseline commit: `78ba1c25beaed7321e175113e0876aac6d763a59`
- Baseline Git clean: `false`
- Baseline source-state SHA-256: `7d72086128b12e2a89d5ca78f0e8cc0cb393a44f1062d027a069722ff9dae007`
- Cargo.lock SHA-256: `37c6f3339dc4eedcc5587e44133a7d0299fef2bd9da52fc0464e2ceeb0138edb`
- Benchmark harness SHA-256: `27acd15e854e46899ff1c1f8e1b399ff240b0b9732e56a36cf8edb39c4a58db3`

**Publication and validation environment**:

- Publication CPU: `Apple M4 Max (arm64)`
- Publication OS: `Darwin 25.6.0 arm64`
- Publication rustc: `rustc 1.98.1 (48a229cea 2026-09-01)`
- Publication commit: `d36a9e99ec223c4cbbf8a9191818f186159c27f7`
- Publication Git clean: `false`
- Publication source-state SHA-256: `b8fcf706fa6e287ca6682a1ecef5a69576413b9c88040184a9b1ecb8f6e92336`
- Publication Cargo.lock SHA-256: `37c6f3339dc4eedcc5587e44133a7d0299fef2bd9da52fc0464e2ceeb0138edb`
- Publication harness SHA-256: `27acd15e854e46899ff1c1f8e1b399ff240b0b9732e56a36cf8edb39c4a58db3`
- Criterion suite/scope: `all` / `release-signal`
- Criterion statistic/sample: `median` / `new`
- Criterion dependency version: `0.8.2`
- Baseline command: `just bench-save-baseline v0.4.5`
- Current command: `just bench-latest`
- Correctness gate: `just test-bench-inputs` passed against both the current and baseline revisions using the shared current fixture harness.
- Validated current revision: `d36a9e99ec223c4cbbf8a9191818f186159c27f7` (Git clean: `false`;
  source-state SHA-256: `b8fcf706fa6e287ca6682a1ecef5a69576413b9c88040184a9b1ecb8f6e92336`)
- Validated baseline revision: `78ba1c25beaed7321e175113e0876aac6d763a59` (Git clean: `false`;
  source-state SHA-256: `7d72086128b12e2a89d5ca78f0e8cc0cb393a44f1062d027a069722ff9dae007`)
- Baseline API compatibility: `la_stack_pre_rational_input_api` selects only source-compatible benchmark calls;
  one-sided rows outside the baseline's correctness domain are identified by the retained CSV coverage status and note.

## Exact arithmetic

| Case | Benchmark | v0.4.5 (point + CI) | Latest (point + CI) | Point-estimate change | CI relation | Point-estimate ratio |
|:-----|:----------|-------:|-------:|-------:|:-----------|--------:|
| D=2 | det | 0.4 ns [0.4 ns, 0.4 ns] | 0.4 ns [0.4 ns, 0.4 ns] | -0.2% | marginal CIs overlap | 1.00x |
| D=2 | det_direct | 0.4 ns [0.4 ns, 0.4 ns] | 0.4 ns [0.4 ns, 0.4 ns] | -0.3% | marginal CIs overlap | 1.00x |
| D=2 | det_direct_with_errbound | 1.7 ns [1.7 ns, 1.7 ns] | 1.7 ns [1.7 ns, 1.7 ns] | -0.0% | marginal CIs overlap | 1.00x |
| D=2 | det_errbound | 1.7 ns [1.6 ns, 1.7 ns] | 1.6 ns [1.6 ns, 1.6 ns] | -2.5% | faster point estimate; marginal CIs separated | 1.03x |
| D=2 | det_exact | 110.5 ns [110.4 ns, 110.7 ns] | 100.4 ns [100.3 ns, 100.5 ns] | -9.2% | faster point estimate; marginal CIs separated | 1.10x |
| D=2 | det_exact_f64_result | 71.3 ns [71.2 ns, 71.5 ns] | 73.1 ns [72.9 ns, 73.2 ns] | +2.5% | slower point estimate; marginal CIs separated | 0.98x |
| D=2 | det_exact_rounded_f64 | 72.7 ns [72.6 ns, 72.7 ns] | 73.3 ns [73.2 ns, 73.4 ns] | +0.8% | slower point estimate; marginal CIs separated | 0.99x |
| D=2 | det_sign_exact | 2.7 ns [2.7 ns, 2.7 ns] | 3.0 ns [3.0 ns, 3.0 ns] | +9.1% | slower point estimate; marginal CIs separated | 0.92x |
| D=2 | solve_exact | 7.06 µs [7.05 µs, 7.07 µs] | 7.10 µs [7.09 µs, 7.12 µs] | +0.7% | slower point estimate; marginal CIs separated | 0.99x |
| D=2 | solve_exact_f64_result | 8.28 µs [8.26 µs, 8.28 µs] | 7.33 µs [7.32 µs, 7.34 µs] | -11.5% | faster point estimate; marginal CIs separated | 1.13x |
| D=2 | solve_exact_rounded_f64 | 7.48 µs [7.47 µs, 7.49 µs] | 7.55 µs [7.54 µs, 7.56 µs] | +0.8% | slower point estimate; marginal CIs separated | 0.99x |
| D=3 | det | 0.8 ns [0.8 ns, 0.8 ns] | 0.8 ns [0.8 ns, 0.8 ns] | +0.1% | marginal CIs overlap | 1.00x |
| D=3 | det_direct | 0.8 ns [0.8 ns, 0.8 ns] | 0.8 ns [0.8 ns, 0.8 ns] | +2.9% | slower point estimate; marginal CIs separated | 0.97x |
| D=3 | det_direct_with_errbound | 3.4 ns [3.3 ns, 3.4 ns] | 3.3 ns [3.3 ns, 3.3 ns] | -0.5% | faster point estimate; marginal CIs separated | 1.01x |
| D=3 | det_errbound | 3.3 ns [3.3 ns, 3.3 ns] | 3.3 ns [3.3 ns, 3.3 ns] | +0.1% | slower point estimate; marginal CIs separated | 1.00x |
| D=3 | det_exact | 370.6 ns [369.9 ns, 371.5 ns] | 338.9 ns [338.5 ns, 339.4 ns] | -8.6% | faster point estimate; marginal CIs separated | 1.09x |
| D=3 | det_exact_f64_result | 338.5 ns [338.0 ns, 339.2 ns] | 313.9 ns [313.2 ns, 314.4 ns] | -7.3% | faster point estimate; marginal CIs separated | 1.08x |
| D=3 | det_exact_rounded_f64 | 340.1 ns [338.8 ns, 341.2 ns] | 317.2 ns [316.0 ns, 318.2 ns] | -6.7% | faster point estimate; marginal CIs separated | 1.07x |
| D=3 | det_sign_exact | 4.7 ns [4.7 ns, 4.7 ns] | 4.7 ns [4.7 ns, 4.7 ns] | -0.0% | marginal CIs overlap | 1.00x |
| D=3 | solve_exact | 30.84 µs [30.78 µs, 30.89 µs] | 30.60 µs [30.51 µs, 30.67 µs] | -0.8% | faster point estimate; marginal CIs separated | 1.01x |
| D=3 | solve_exact_f64_result | 32.80 µs [32.76 µs, 32.87 µs] | 30.92 µs [30.86 µs, 30.96 µs] | -5.7% | faster point estimate; marginal CIs separated | 1.06x |
| D=3 | solve_exact_rounded_f64 | 31.46 µs [31.35 µs, 31.60 µs] | 31.60 µs [31.55 µs, 31.72 µs] | +0.5% | marginal CIs overlap | 1.00x |
| D=4 | det | 3.2 ns [3.2 ns, 3.2 ns] | 3.2 ns [3.2 ns, 3.2 ns] | -0.1% | marginal CIs overlap | 1.00x |
| D=4 | det_direct | 2.4 ns [2.4 ns, 2.4 ns] | 2.4 ns [2.4 ns, 2.4 ns] | -0.0% | marginal CIs overlap | 1.00x |
| D=4 | det_direct_with_errbound | 6.8 ns [6.8 ns, 6.8 ns] | 6.8 ns [6.8 ns, 6.8 ns] | +0.2% | slower point estimate; marginal CIs separated | 1.00x |
| D=4 | det_errbound | 6.7 ns [6.7 ns, 6.7 ns] | 6.7 ns [6.7 ns, 6.7 ns] | +0.1% | slower point estimate; marginal CIs separated | 1.00x |
| D=4 | det_exact | 1.22 µs [1.21 µs, 1.22 µs] | 1.01 µs [1.01 µs, 1.01 µs] | -17.0% | faster point estimate; marginal CIs separated | 1.20x |
| D=4 | det_exact_f64_result | 1.17 µs [1.16 µs, 1.17 µs] | 1.06 µs [1.05 µs, 1.06 µs] | -9.6% | faster point estimate; marginal CIs separated | 1.11x |
| D=4 | det_exact_rounded_f64 | 1.17 µs [1.17 µs, 1.17 µs] | 1.01 µs [1.01 µs, 1.01 µs] | -13.8% | faster point estimate; marginal CIs separated | 1.16x |
| D=4 | det_sign_exact | 7.7 ns [7.7 ns, 7.7 ns] | 8.0 ns [8.0 ns, 8.0 ns] | +3.6% | slower point estimate; marginal CIs separated | 0.97x |
| D=4 | solve_exact | 80.39 µs [80.27 µs, 80.49 µs] | 79.18 µs [79.08 µs, 79.34 µs] | -1.5% | faster point estimate; marginal CIs separated | 1.02x |
| D=4 | solve_exact_f64_result | 83.39 µs [83.24 µs, 83.63 µs] | 79.34 µs [79.13 µs, 79.50 µs] | -4.9% | faster point estimate; marginal CIs separated | 1.05x |
| D=4 | solve_exact_rounded_f64 | 81.02 µs [80.74 µs, 81.52 µs] | 79.86 µs [79.66 µs, 80.05 µs] | -1.4% | faster point estimate; marginal CIs separated | 1.01x |
| D=5 | det | 24.7 ns [24.7 ns, 24.8 ns] | 25.0 ns [25.0 ns, 25.1 ns] | +1.4% | slower point estimate; marginal CIs separated | 0.99x |
| D=5 | det_exact | 3.54 µs [3.52 µs, 3.54 µs] | 3.74 µs [3.71 µs, 3.79 µs] | +5.7% | slower point estimate; marginal CIs separated | 0.95x |
| D=5 | det_exact_f64_result | 3.49 µs [3.48 µs, 3.50 µs] | 3.54 µs [3.54 µs, 3.55 µs] | +1.7% | slower point estimate; marginal CIs separated | 0.98x |
| D=5 | det_exact_rounded_f64 | 3.47 µs [3.46 µs, 3.48 µs] | 3.61 µs [3.60 µs, 3.61 µs] | +3.8% | slower point estimate; marginal CIs separated | 0.96x |
| D=5 | det_sign_exact | 3.58 µs [3.57 µs, 3.59 µs] | 3.70 µs [3.70 µs, 3.71 µs] | +3.5% | slower point estimate; marginal CIs separated | 0.97x |
| D=5 | solve_exact | 159.88 µs [159.63 µs, 160.24 µs] | 157.60 µs [157.28 µs, 157.89 µs] | -1.4% | faster point estimate; marginal CIs separated | 1.01x |
| D=5 | solve_exact_f64_result | 163.37 µs [162.91 µs, 163.71 µs] | 157.72 µs [157.38 µs, 157.93 µs] | -3.5% | faster point estimate; marginal CIs separated | 1.04x |
| D=5 | solve_exact_rounded_f64 | 161.50 µs [161.20 µs, 161.73 µs] | 164.26 µs [163.54 µs, 164.87 µs] | +1.7% | slower point estimate; marginal CIs separated | 0.98x |
| Hilbert 4x4 | det_exact | 1.34 µs [1.33 µs, 1.35 µs] | 1.04 µs [1.04 µs, 1.05 µs] | -22.3% | faster point estimate; marginal CIs separated | 1.29x |
| Hilbert 4x4 | det_sign_exact | 7.7 ns [7.7 ns, 7.7 ns] | 8.0 ns [8.0 ns, 8.0 ns] | +3.4% | slower point estimate; marginal CIs separated | 0.97x |
| Hilbert 4x4 | solve_exact | 60.21 µs [60.13 µs, 60.29 µs] | 59.56 µs [59.47 µs, 59.63 µs] | -1.1% | faster point estimate; marginal CIs separated | 1.01x |
| Hilbert 4x4 | solve_exact_f64_result | 62.75 µs [62.65 µs, 62.89 µs] | 59.74 µs [59.59 µs, 59.79 µs] | -4.8% | faster point estimate; marginal CIs separated | 1.05x |
| Hilbert 4x4 | solve_exact_rounded_f64 | 61.09 µs [60.95 µs, 61.18 µs] | 60.47 µs [60.26 µs, 60.62 µs] | -1.0% | faster point estimate; marginal CIs separated | 1.01x |
| Hilbert 5x5 | det_exact | 3.68 µs [3.66 µs, 3.69 µs] | 3.75 µs [3.74 µs, 3.76 µs] | +1.9% | slower point estimate; marginal CIs separated | 0.98x |
| Hilbert 5x5 | det_sign_exact | 3.87 µs [3.87 µs, 3.88 µs] | 3.92 µs [3.92 µs, 3.93 µs] | +1.1% | slower point estimate; marginal CIs separated | 0.99x |
| Hilbert 5x5 | solve_exact | 124.69 µs [124.44 µs, 124.97 µs] | 122.82 µs [122.53 µs, 122.95 µs] | -1.5% | faster point estimate; marginal CIs separated | 1.02x |
| Hilbert 5x5 | solve_exact_f64_result | 127.73 µs [127.44 µs, 128.19 µs] | 122.06 µs [121.91 µs, 122.40 µs] | -4.4% | faster point estimate; marginal CIs separated | 1.05x |
| Hilbert 5x5 | solve_exact_rounded_f64 | 126.57 µs [126.11 µs, 126.93 µs] | 124.52 µs [124.21 µs, 124.64 µs] | -1.6% | faster point estimate; marginal CIs separated | 1.02x |
| Large entries 3x3 | det_exact | 338.5 ns [337.7 ns, 339.5 ns] | 341.6 ns [340.6 ns, 343.9 ns] | +0.9% | slower point estimate; marginal CIs separated | 0.99x |
| Large entries 3x3 | det_sign_exact | 320.8 ns [319.6 ns, 323.0 ns] | 334.9 ns [334.1 ns, 337.0 ns] | +4.4% | slower point estimate; marginal CIs separated | 0.96x |
| Large entries 3x3 | solve_exact | 93.99 µs [93.81 µs, 94.26 µs] | 93.14 µs [92.93 µs, 93.30 µs] | -0.9% | faster point estimate; marginal CIs separated | 1.01x |
| Large entries 3x3 | solve_exact_f64_result | 94.71 µs [94.65 µs, 94.80 µs] | 93.37 µs [93.24 µs, 93.46 µs] | -1.4% | faster point estimate; marginal CIs separated | 1.01x |
| Large entries 3x3 | solve_exact_rounded_f64 | 94.49 µs [94.43 µs, 94.61 µs] | 93.32 µs [93.20 µs, 93.38 µs] | -1.2% | faster point estimate; marginal CIs separated | 1.01x |
| Near-singular 3x3 | det_exact | 342.0 ns [341.2 ns, 342.8 ns] | 334.5 ns [333.8 ns, 336.3 ns] | -2.2% | faster point estimate; marginal CIs separated | 1.02x |
| Near-singular 3x3 | det_sign_exact | 398.3 ns [394.8 ns, 401.2 ns] | 360.9 ns [360.5 ns, 361.6 ns] | -9.4% | faster point estimate; marginal CIs separated | 1.10x |
| Near-singular 3x3 | solve_exact | 2.60 µs [2.60 µs, 2.61 µs] | 2.81 µs [2.80 µs, 2.81 µs] | +7.9% | slower point estimate; marginal CIs separated | 0.93x |
| Near-singular 3x3 | solve_exact_f64_result | 2.70 µs [2.70 µs, 2.71 µs] | 2.81 µs [2.80 µs, 2.81 µs] | +3.8% | slower point estimate; marginal CIs separated | 0.96x |
| Near-singular 3x3 | solve_exact_rounded_f64 | 2.60 µs [2.60 µs, 2.61 µs] | 2.78 µs [2.77 µs, 2.80 µs] | +6.9% | slower point estimate; marginal CIs separated | 0.94x |
| Random corpus D=2 | det_exact | 3.12 µs [3.12 µs, 3.12 µs] | 3.16 µs [3.14 µs, 3.18 µs] | +1.2% | slower point estimate; marginal CIs separated | 0.99x |
| Random corpus D=2 | det_sign_exact | 136.4 ns [135.6 ns, 137.0 ns] | 149.9 ns [148.9 ns, 150.9 ns] | +9.9% | slower point estimate; marginal CIs separated | 0.91x |
| Random corpus D=2 | solve_exact | 67.47 µs [67.33 µs, 67.59 µs] | 78.46 µs [78.37 µs, 78.65 µs] | +16.3% | slower point estimate; marginal CIs separated | 0.86x |
| Random corpus D=2 | solve_exact_f64_result | 78.30 µs [77.84 µs, 78.74 µs] | 79.81 µs [79.67 µs, 79.96 µs] | +1.9% | slower point estimate; marginal CIs separated | 0.98x |
| Random corpus D=2 | solve_exact_rounded_f64 | 68.56 µs [68.41 µs, 68.68 µs] | 82.23 µs [82.02 µs, 82.37 µs] | +19.9% | slower point estimate; marginal CIs separated | 0.83x |
| Random corpus D=3 | det_exact | 8.35 µs [8.08 µs, 8.64 µs] | 8.47 µs [8.44 µs, 8.49 µs] | +1.5% | marginal CIs overlap | 0.99x |
| Random corpus D=3 | det_sign_exact | 225.1 ns [225.0 ns, 225.2 ns] | 227.9 ns [227.4 ns, 228.3 ns] | +1.2% | slower point estimate; marginal CIs separated | 0.99x |
| Random corpus D=3 | solve_exact | 214.24 µs [213.80 µs, 214.55 µs] | 239.62 µs [239.14 µs, 240.19 µs] | +11.9% | slower point estimate; marginal CIs separated | 0.89x |
| Random corpus D=3 | solve_exact_f64_result | 229.92 µs [229.52 µs, 230.51 µs] | 240.54 µs [240.00 µs, 241.28 µs] | +4.6% | slower point estimate; marginal CIs separated | 0.96x |
| Random corpus D=3 | solve_exact_rounded_f64 | 217.63 µs [217.29 µs, 217.94 µs] | 242.67 µs [242.06 µs, 243.32 µs] | +11.5% | slower point estimate; marginal CIs separated | 0.90x |
| Random corpus D=4 | det_exact | 25.26 µs [25.21 µs, 25.30 µs] | 21.20 µs [21.16 µs, 21.27 µs] | -16.1% | faster point estimate; marginal CIs separated | 1.19x |
| Random corpus D=4 | det_sign_exact | 424.5 ns [423.1 ns, 425.6 ns] | 454.6 ns [452.7 ns, 456.0 ns] | +7.1% | slower point estimate; marginal CIs separated | 0.93x |
| Random corpus D=4 | solve_exact | 489.45 µs [488.84 µs, 490.29 µs] | 535.37 µs [533.84 µs, 537.11 µs] | +9.4% | slower point estimate; marginal CIs separated | 0.91x |
| Random corpus D=4 | solve_exact_f64_result | 507.75 µs [507.05 µs, 508.51 µs] | 534.53 µs [533.89 µs, 535.33 µs] | +5.3% | slower point estimate; marginal CIs separated | 0.95x |
| Random corpus D=4 | solve_exact_rounded_f64 | 491.51 µs [490.45 µs, 491.98 µs] | 535.39 µs [535.09 µs, 535.96 µs] | +8.9% | slower point estimate; marginal CIs separated | 0.92x |
| Random corpus D=5 | det_exact | 47.41 µs [47.34 µs, 47.50 µs] | 55.99 µs [55.88 µs, 56.13 µs] | +18.1% | slower point estimate; marginal CIs separated | 0.85x |
| Random corpus D=5 | det_sign_exact | 49.54 µs [48.76 µs, 50.88 µs] | 61.93 µs [61.88 µs, 62.04 µs] | +25.0% | slower point estimate; marginal CIs separated | 0.80x |
| Random corpus D=5 | solve_exact | 965.82 µs [964.45 µs, 966.73 µs] | 1.04 ms [1.04 ms, 1.04 ms] | +8.1% | slower point estimate; marginal CIs separated | 0.92x |
| Random corpus D=5 | solve_exact_f64_result | 987.58 µs [986.29 µs, 989.22 µs] | 1.04 ms [1.04 ms, 1.04 ms] | +5.5% | slower point estimate; marginal CIs separated | 0.95x |
| Random corpus D=5 | solve_exact_rounded_f64 | 966.85 µs [965.47 µs, 968.31 µs] | 1.04 ms [1.04 ms, 1.05 ms] | +8.1% | slower point estimate; marginal CIs separated | 0.93x |

## vs_linalg

| Case | Benchmark | v0.4.5 (point + CI) | Latest (point + CI) | Point-estimate change | CI relation | Point-estimate ratio | v0.4.5 nalgebra | v0.4.5 faer |
|:-----|:----------|-------:|-------:|-------:|:-----------|--------:|-------:|-------:|
| D=16 | la_stack_det | 433.3 ns [431.0 ns, 435.4 ns] | 463.8 ns [460.4 ns, 472.3 ns] | +7.0% | slower point estimate; marginal CIs separated | 0.93x | — | — |
| D=16 | la_stack_det_from_ldlt | 2.9 ns [2.9 ns, 2.9 ns] | 2.9 ns [2.9 ns, 2.9 ns] | -0.1% | marginal CIs overlap | 1.00x | 1.7 ns [1.7 ns, 1.7 ns] | 4.3 ns [4.3 ns, 4.4 ns] |
| D=16 | la_stack_det_from_lu | 3.3 ns [3.3 ns, 3.3 ns] | 3.3 ns [3.3 ns, 3.3 ns] | +0.2% | marginal CIs overlap | 1.00x | 1.8 ns [1.8 ns, 1.8 ns] | 4.8 ns [4.7 ns, 4.8 ns] |
| D=16 | la_stack_det_via_lu | 403.1 ns [402.1 ns, 405.0 ns] | 405.2 ns [403.8 ns, 406.9 ns] | +0.5% | marginal CIs overlap | 0.99x | 449.5 ns [448.4 ns, 450.3 ns] | 671.6 ns [665.6 ns, 677.2 ns] |
| D=16 | la_stack_dot | 2.3 ns [2.3 ns, 2.3 ns] | 2.3 ns [2.3 ns, 2.3 ns] | +0.5% | slower point estimate; marginal CIs separated | 1.00x | 1.9 ns [1.9 ns, 1.9 ns] | 3.3 ns [3.3 ns, 3.3 ns] |
| D=16 | la_stack_inf_norm | 32.2 ns [32.2 ns, 32.3 ns] | 32.2 ns [32.1 ns, 32.2 ns] | -0.3% | marginal CIs overlap | 1.00x | 31.4 ns [31.3 ns, 31.4 ns] | 32.7 ns [32.5 ns, 32.7 ns] |
| D=16 | la_stack_ldlt | 383.4 ns [382.4 ns, 383.9 ns] | 403.8 ns [403.4 ns, 404.5 ns] | +5.3% | slower point estimate; marginal CIs separated | 0.95x | 397.0 ns [396.6 ns, 397.6 ns] | 428.2 ns [426.1 ns, 430.8 ns] |
| D=16 | la_stack_ldlt_solve | 425.4 ns [424.9 ns, 425.8 ns] | 447.6 ns [446.3 ns, 448.1 ns] | +5.2% | slower point estimate; marginal CIs separated | 0.95x | 607.4 ns [606.6 ns, 608.8 ns] | 626.1 ns [624.9 ns, 627.1 ns] |
| D=16 | la_stack_lu | 393.8 ns [392.8 ns, 395.0 ns] | 395.6 ns [394.4 ns, 397.1 ns] | +0.5% | marginal CIs overlap | 1.00x | 453.8 ns [452.9 ns, 454.1 ns] | 641.8 ns [639.4 ns, 645.7 ns] |
| D=16 | la_stack_lu_solve | 640.8 ns [639.4 ns, 641.5 ns] | 643.5 ns [641.5 ns, 644.9 ns] | +0.4% | slower point estimate; marginal CIs separated | 1.00x | 569.4 ns [568.9 ns, 569.9 ns] | 903.1 ns [900.3 ns, 905.3 ns] |
| D=16 | la_stack_norm2_sq | 2.1 ns [2.1 ns, 2.1 ns] | 2.1 ns [2.1 ns, 2.1 ns] | +0.0% | marginal CIs overlap | 1.00x | 1.5 ns [1.5 ns, 1.5 ns] | 4.1 ns [4.1 ns, 4.1 ns] |
| D=16 | la_stack_solve_from_ldlt | 27.7 ns [27.7 ns, 27.8 ns] | 27.7 ns [27.7 ns, 27.7 ns] | -0.1% | marginal CIs overlap | 1.00x | 116.1 ns [115.7 ns, 116.2 ns] | 181.7 ns [181.4 ns, 182.2 ns] |
| D=16 | la_stack_solve_from_lu | 188.0 ns [187.6 ns, 188.6 ns] | 188.8 ns [188.3 ns, 189.4 ns] | +0.4% | marginal CIs overlap | 1.00x | 92.2 ns [92.1 ns, 92.3 ns] | 258.8 ns [255.7 ns, 262.5 ns] |
| D=2 | la_stack_det | 0.6 ns [0.6 ns, 0.6 ns] | 0.6 ns [0.6 ns, 0.6 ns] | +0.0% | marginal CIs overlap | 1.00x | — | — |
| D=2 | la_stack_det_from_ldlt | 0.5 ns [0.5 ns, 0.5 ns] | 0.5 ns [0.5 ns, 0.5 ns] | +0.0% | marginal CIs overlap | 1.00x | 0.4 ns [0.4 ns, 0.4 ns] | 0.5 ns [0.5 ns, 0.6 ns] |
| D=2 | la_stack_det_from_lu | 0.5 ns [0.5 ns, 0.5 ns] | 0.5 ns [0.5 ns, 0.5 ns] | -0.1% | marginal CIs overlap | 1.00x | 0.4 ns [0.4 ns, 0.4 ns] | 0.7 ns [0.7 ns, 0.7 ns] |
| D=2 | la_stack_det_via_lu | 1.8 ns [1.8 ns, 1.8 ns] | 1.8 ns [1.8 ns, 1.8 ns] | +0.3% | slower point estimate; marginal CIs separated | 1.00x | 0.8 ns [0.8 ns, 0.8 ns] | 131.7 ns [131.1 ns, 132.0 ns] |
| D=2 | la_stack_dot | 0.6 ns [0.6 ns, 0.6 ns] | 0.6 ns [0.6 ns, 0.6 ns] | -0.2% | marginal CIs overlap | 1.00x | 0.6 ns [0.6 ns, 0.6 ns] | 2.5 ns [2.5 ns, 2.7 ns] |
| D=2 | la_stack_inf_norm | 0.6 ns [0.6 ns, 0.6 ns] | 0.6 ns [0.6 ns, 0.6 ns] | -0.0% | marginal CIs overlap | 1.00x | 0.5 ns [0.5 ns, 0.5 ns] | 0.7 ns [0.7 ns, 0.7 ns] |
| D=2 | la_stack_ldlt | 6.5 ns [6.5 ns, 6.5 ns] | 6.6 ns [6.6 ns, 6.6 ns] | +2.6% | slower point estimate; marginal CIs separated | 0.97x | 1.7 ns [1.7 ns, 1.7 ns] | 98.0 ns [97.6 ns, 98.4 ns] |
| D=2 | la_stack_ldlt_solve | 9.3 ns [9.3 ns, 9.3 ns] | 9.4 ns [9.4 ns, 9.5 ns] | +1.2% | slower point estimate; marginal CIs separated | 0.99x | 2.5 ns [2.5 ns, 2.5 ns] | 145.9 ns [143.5 ns, 147.3 ns] |
| D=2 | la_stack_lu | 1.6 ns [1.6 ns, 1.6 ns] | 1.6 ns [1.6 ns, 1.6 ns] | -0.3% | faster point estimate; marginal CIs separated | 1.00x | 1.6 ns [1.6 ns, 1.6 ns] | 120.9 ns [120.0 ns, 121.7 ns] |
| D=2 | la_stack_lu_solve | 2.0 ns [2.0 ns, 2.0 ns] | 2.0 ns [2.0 ns, 2.0 ns] | -0.5% | marginal CIs overlap | 1.01x | 4.5 ns [4.5 ns, 4.5 ns] | 182.5 ns [182.1 ns, 183.0 ns] |
| D=2 | la_stack_norm2_sq | 0.4 ns [0.4 ns, 0.4 ns] | 0.4 ns [0.4 ns, 0.4 ns] | -0.0% | marginal CIs overlap | 1.00x | 0.4 ns [0.4 ns, 0.4 ns] | 4.1 ns [4.1 ns, 4.1 ns] |
| D=2 | la_stack_solve_from_ldlt | 1.2 ns [1.2 ns, 1.2 ns] | 1.2 ns [1.2 ns, 1.2 ns] | +1.6% | slower point estimate; marginal CIs separated | 0.98x | 1.2 ns [1.2 ns, 1.2 ns] | 44.0 ns [43.7 ns, 44.2 ns] |
| D=2 | la_stack_solve_from_lu | 1.2 ns [1.2 ns, 1.2 ns] | 1.2 ns [1.2 ns, 1.2 ns] | -0.0% | marginal CIs overlap | 1.00x | 2.9 ns [2.9 ns, 2.9 ns] | 54.8 ns [54.6 ns, 54.9 ns] |
| D=3 | la_stack_det | 1.3 ns [1.3 ns, 1.3 ns] | 1.3 ns [1.3 ns, 1.3 ns] | -0.2% | marginal CIs overlap | 1.00x | — | — |
| D=3 | la_stack_det_from_ldlt | 0.6 ns [0.6 ns, 0.6 ns] | 0.6 ns [0.6 ns, 0.6 ns] | +0.0% | marginal CIs overlap | 1.00x | 0.5 ns [0.5 ns, 0.5 ns] | 0.7 ns [0.7 ns, 0.7 ns] |
| D=3 | la_stack_det_from_lu | 0.7 ns [0.7 ns, 0.7 ns] | 0.7 ns [0.7 ns, 0.7 ns] | +0.0% | marginal CIs overlap | 1.00x | 0.5 ns [0.5 ns, 0.5 ns] | 1.0 ns [1.0 ns, 1.0 ns] |
| D=3 | la_stack_det_via_lu | 8.7 ns [8.7 ns, 8.7 ns] | 8.7 ns [8.7 ns, 8.8 ns] | +0.4% | marginal CIs overlap | 1.00x | 17.8 ns [17.7 ns, 17.9 ns] | 162.9 ns [161.8 ns, 164.2 ns] |
| D=3 | la_stack_dot | 0.7 ns [0.7 ns, 0.7 ns] | 0.7 ns [0.7 ns, 0.7 ns] | +0.1% | marginal CIs overlap | 1.00x | 0.7 ns [0.7 ns, 0.7 ns] | 3.0 ns [3.0 ns, 3.0 ns] |
| D=3 | la_stack_inf_norm | 1.3 ns [1.3 ns, 1.3 ns] | 1.3 ns [1.3 ns, 1.3 ns] | +0.0% | marginal CIs overlap | 1.00x | 1.1 ns [1.1 ns, 1.1 ns] | 1.3 ns [1.2 ns, 1.3 ns] |
| D=3 | la_stack_ldlt | 13.4 ns [13.3 ns, 13.4 ns] | 13.6 ns [13.5 ns, 13.7 ns] | +1.6% | slower point estimate; marginal CIs separated | 0.98x | 3.9 ns [3.8 ns, 3.9 ns] | 108.2 ns [107.8 ns, 109.1 ns] |
| D=3 | la_stack_ldlt_solve | 12.5 ns [12.5 ns, 12.6 ns] | 18.3 ns [17.8 ns, 20.3 ns] | +45.7% | slower point estimate; marginal CIs separated | 0.69x | 5.6 ns [5.6 ns, 5.6 ns] | 155.4 ns [154.6 ns, 155.9 ns] |
| D=3 | la_stack_lu | 8.3 ns [8.3 ns, 8.3 ns] | 17.7 ns [16.7 ns, 18.4 ns] | +114.1% | slower point estimate; marginal CIs separated | 0.47x | 14.8 ns [14.7 ns, 14.8 ns] | 156.6 ns [155.4 ns, 158.0 ns] |
| D=3 | la_stack_lu_solve | 10.1 ns [10.0 ns, 10.2 ns] | 10.0 ns [9.9 ns, 10.0 ns] | -1.1% | faster point estimate; marginal CIs separated | 1.01x | 22.8 ns [22.7 ns, 23.0 ns] | 217.6 ns [216.2 ns, 219.0 ns] |
| D=3 | la_stack_norm2_sq | 0.4 ns [0.4 ns, 0.4 ns] | 0.4 ns [0.4 ns, 0.4 ns] | +0.1% | marginal CIs overlap | 1.00x | 0.4 ns [0.4 ns, 0.4 ns] | 4.1 ns [4.1 ns, 4.1 ns] |
| D=3 | la_stack_solve_from_ldlt | 1.8 ns [1.8 ns, 1.8 ns] | 1.8 ns [1.8 ns, 1.8 ns] | +0.0% | marginal CIs overlap | 1.00x | 2.7 ns [2.7 ns, 2.8 ns] | 45.8 ns [45.4 ns, 46.5 ns] |
| D=3 | la_stack_solve_from_lu | 2.1 ns [2.1 ns, 2.1 ns] | 2.1 ns [2.1 ns, 2.1 ns] | +0.1% | marginal CIs overlap | 1.00x | 4.5 ns [4.5 ns, 4.5 ns] | 56.7 ns [56.5 ns, 57.0 ns] |
| D=32 | la_stack_det | 2.18 µs [2.17 µs, 2.18 µs] | 2.17 µs [2.16 µs, 2.18 µs] | -0.5% | marginal CIs overlap | 1.01x | — | — |
| D=32 | la_stack_det_from_ldlt | 7.7 ns [7.7 ns, 7.7 ns] | 7.8 ns [7.7 ns, 7.8 ns] | +0.6% | slower point estimate; marginal CIs separated | 0.99x | 3.0 ns [3.0 ns, 3.0 ns] | 8.3 ns [8.2 ns, 8.3 ns] |
| D=32 | la_stack_det_from_lu | 9.1 ns [9.0 ns, 9.1 ns] | 9.0 ns [9.0 ns, 9.0 ns] | -0.1% | marginal CIs overlap | 1.00x | 3.1 ns [3.1 ns, 3.1 ns] | 8.6 ns [8.6 ns, 8.7 ns] |
| D=32 | la_stack_det_via_lu | 2.05 µs [2.04 µs, 2.05 µs] | 2.05 µs [2.04 µs, 2.06 µs] | +0.3% | marginal CIs overlap | 1.00x | 2.19 µs [2.19 µs, 2.20 µs] | 2.31 µs [2.30 µs, 2.32 µs] |
| D=32 | la_stack_dot | 4.0 ns [4.0 ns, 4.0 ns] | 4.0 ns [4.0 ns, 4.0 ns] | -0.1% | marginal CIs overlap | 1.00x | 4.7 ns [4.7 ns, 4.7 ns] | 4.9 ns [4.9 ns, 4.9 ns] |
| D=32 | la_stack_inf_norm | 127.0 ns [126.9 ns, 127.2 ns] | 127.6 ns [127.4 ns, 127.7 ns] | +0.5% | slower point estimate; marginal CIs separated | 0.99x | 158.5 ns [154.9 ns, 159.6 ns] | 163.6 ns [163.4 ns, 164.2 ns] |
| D=32 | la_stack_ldlt | 2.44 µs [2.43 µs, 2.44 µs] | 2.45 µs [2.45 µs, 2.45 µs] | +0.5% | slower point estimate; marginal CIs separated | 0.99x | 2.06 µs [2.05 µs, 2.06 µs] | 1.45 µs [1.45 µs, 1.46 µs] |
| D=32 | la_stack_ldlt_solve | 2.78 µs [2.77 µs, 2.79 µs] | 2.77 µs [2.76 µs, 2.79 µs] | -0.4% | marginal CIs overlap | 1.00x | 2.69 µs [2.69 µs, 2.69 µs] | 1.94 µs [1.94 µs, 1.94 µs] |
| D=32 | la_stack_lu | 2.11 µs [2.10 µs, 2.12 µs] | 2.12 µs [2.11 µs, 2.13 µs] | +0.4% | marginal CIs overlap | 1.00x | 2.14 µs [2.13 µs, 2.14 µs] | 2.26 µs [2.26 µs, 2.27 µs] |
| D=32 | la_stack_lu_solve | 2.78 µs [2.77 µs, 2.80 µs] | 2.76 µs [2.75 µs, 2.76 µs] | -1.0% | faster point estimate; marginal CIs separated | 1.01x | 2.71 µs [2.70 µs, 2.72 µs] | 2.92 µs [2.91 µs, 2.93 µs] |
| D=32 | la_stack_norm2_sq | 4.0 ns [4.0 ns, 4.0 ns] | 4.0 ns [4.0 ns, 4.0 ns] | +0.0% | marginal CIs overlap | 1.00x | 3.8 ns [3.8 ns, 3.8 ns] | 4.2 ns [4.2 ns, 4.2 ns] |
| D=32 | la_stack_solve_from_ldlt | 307.5 ns [307.0 ns, 307.8 ns] | 310.4 ns [309.5 ns, 311.1 ns] | +0.9% | slower point estimate; marginal CIs separated | 0.99x | 533.2 ns [531.2 ns, 534.3 ns] | 466.9 ns [466.3 ns, 467.6 ns] |
| D=32 | la_stack_solve_from_lu | 661.0 ns [660.0 ns, 662.2 ns] | 652.5 ns [650.7 ns, 653.7 ns] | -1.3% | faster point estimate; marginal CIs separated | 1.01x | 327.4 ns [327.0 ns, 327.8 ns] | 631.0 ns [629.5 ns, 632.5 ns] |
| D=4 | la_stack_det | 2.6 ns [2.6 ns, 2.6 ns] | 2.6 ns [2.6 ns, 2.6 ns] | +0.1% | slower point estimate; marginal CIs separated | 1.00x | — | — |
| D=4 | la_stack_det_from_ldlt | 0.8 ns [0.8 ns, 0.8 ns] | 0.8 ns [0.8 ns, 0.8 ns] | +0.0% | marginal CIs overlap | 1.00x | 0.5 ns [0.5 ns, 0.5 ns] | 1.0 ns [1.0 ns, 1.0 ns] |
| D=4 | la_stack_det_from_lu | 0.9 ns [0.9 ns, 0.9 ns] | 0.9 ns [0.9 ns, 0.9 ns] | +0.0% | marginal CIs overlap | 1.00x | 0.6 ns [0.6 ns, 0.6 ns] | 1.2 ns [1.2 ns, 1.2 ns] |
| D=4 | la_stack_det_via_lu | 14.4 ns [14.4 ns, 14.5 ns] | 14.5 ns [14.4 ns, 14.5 ns] | +0.1% | marginal CIs overlap | 1.00x | 30.0 ns [30.0 ns, 30.0 ns] | 179.5 ns [179.0 ns, 180.9 ns] |
| D=4 | la_stack_dot | 0.7 ns [0.7 ns, 0.7 ns] | 0.7 ns [0.7 ns, 0.7 ns] | +0.0% | marginal CIs overlap | 1.00x | 0.6 ns [0.6 ns, 0.6 ns] | 2.6 ns [2.6 ns, 2.7 ns] |
| D=4 | la_stack_inf_norm | 2.2 ns [2.2 ns, 2.2 ns] | 2.2 ns [2.2 ns, 2.2 ns] | +0.1% | marginal CIs overlap | 1.00x | 2.0 ns [2.0 ns, 2.0 ns] | 2.0 ns [2.0 ns, 2.0 ns] |
| D=4 | la_stack_ldlt | 21.1 ns [21.1 ns, 21.2 ns] | 21.3 ns [21.3 ns, 21.4 ns] | +0.8% | slower point estimate; marginal CIs separated | 0.99x | 7.5 ns [7.4 ns, 7.5 ns] | 126.4 ns [126.1 ns, 126.9 ns] |
| D=4 | la_stack_ldlt_solve | 22.9 ns [22.9 ns, 23.0 ns] | 38.6 ns [38.0 ns, 39.7 ns] | +68.2% | slower point estimate; marginal CIs separated | 0.59x | 10.4 ns [10.4 ns, 10.4 ns] | 172.6 ns [172.4 ns, 173.2 ns] |
| D=4 | la_stack_lu | 14.0 ns [13.9 ns, 14.0 ns] | 13.9 ns [13.8 ns, 13.9 ns] | -0.9% | faster point estimate; marginal CIs separated | 1.01x | 29.1 ns [29.1 ns, 29.2 ns] | 175.5 ns [174.3 ns, 176.3 ns] |
| D=4 | la_stack_lu_solve | 21.9 ns [21.9 ns, 22.0 ns] | 22.1 ns [22.1 ns, 22.1 ns] | +0.8% | slower point estimate; marginal CIs separated | 0.99x | 51.7 ns [51.5 ns, 51.9 ns] | 241.5 ns [239.3 ns, 243.2 ns] |
| D=4 | la_stack_norm2_sq | 0.5 ns [0.5 ns, 0.5 ns] | 0.5 ns [0.5 ns, 0.5 ns] | +0.1% | marginal CIs overlap | 1.00x | 0.5 ns [0.5 ns, 0.5 ns] | 4.1 ns [4.1 ns, 4.1 ns] |
| D=4 | la_stack_solve_from_ldlt | 2.5 ns [2.5 ns, 2.5 ns] | 2.5 ns [2.5 ns, 2.5 ns] | -0.0% | marginal CIs overlap | 1.00x | 5.3 ns [5.2 ns, 5.3 ns] | 44.8 ns [44.6 ns, 44.9 ns] |
| D=4 | la_stack_solve_from_lu | 4.0 ns [4.0 ns, 4.0 ns] | 3.9 ns [3.9 ns, 3.9 ns] | -1.8% | faster point estimate; marginal CIs separated | 1.02x | 5.7 ns [5.6 ns, 5.7 ns] | 59.9 ns [59.5 ns, 60.2 ns] |
| D=5 | la_stack_det | 39.4 ns [39.3 ns, 39.6 ns] | 40.6 ns [40.0 ns, 58.5 ns] | +2.9% | slower point estimate; marginal CIs separated | 0.97x | — | — |
| D=5 | la_stack_det_from_ldlt | 0.9 ns [0.9 ns, 0.9 ns] | 0.9 ns [0.9 ns, 0.9 ns] | +0.1% | marginal CIs overlap | 1.00x | 0.6 ns [0.6 ns, 0.6 ns] | 1.2 ns [1.2 ns, 1.2 ns] |
| D=5 | la_stack_det_from_lu | 1.1 ns [1.1 ns, 1.1 ns] | 1.1 ns [1.1 ns, 1.1 ns] | +0.1% | marginal CIs overlap | 1.00x | 0.7 ns [0.7 ns, 0.7 ns] | 1.5 ns [1.5 ns, 1.5 ns] |
| D=5 | la_stack_det_via_lu | 32.1 ns [32.0 ns, 32.3 ns] | 32.7 ns [32.6 ns, 32.7 ns] | +1.7% | slower point estimate; marginal CIs separated | 0.98x | 54.8 ns [54.7 ns, 54.9 ns] | 223.3 ns [221.9 ns, 224.9 ns] |
| D=5 | la_stack_dot | 0.8 ns [0.8 ns, 0.8 ns] | 0.8 ns [0.8 ns, 0.8 ns] | +5.2% | slower point estimate; marginal CIs separated | 0.95x | 0.7 ns [0.7 ns, 0.8 ns] | 2.8 ns [2.8 ns, 2.8 ns] |
| D=5 | la_stack_inf_norm | 3.4 ns [3.4 ns, 3.4 ns] | 3.5 ns [3.5 ns, 3.5 ns] | +1.9% | slower point estimate; marginal CIs separated | 0.98x | 3.2 ns [3.2 ns, 3.2 ns] | 3.2 ns [3.2 ns, 3.2 ns] |
| D=5 | la_stack_ldlt | 40.5 ns [40.4 ns, 40.6 ns] | 41.9 ns [41.6 ns, 42.3 ns] | +3.6% | slower point estimate; marginal CIs separated | 0.97x | 26.1 ns [25.1 ns, 26.5 ns] | 144.0 ns [143.6 ns, 144.6 ns] |
| D=5 | la_stack_ldlt_solve | 52.3 ns [51.9 ns, 52.5 ns] | 45.4 ns [45.3 ns, 45.5 ns] | -13.3% | faster point estimate; marginal CIs separated | 1.15x | 40.4 ns [40.3 ns, 40.4 ns] | 234.3 ns [231.2 ns, 246.0 ns] |
| D=5 | la_stack_lu | 32.0 ns [31.9 ns, 32.1 ns] | 32.3 ns [32.3 ns, 32.4 ns] | +1.0% | slower point estimate; marginal CIs separated | 0.99x | 55.1 ns [54.9 ns, 55.2 ns] | 217.7 ns [214.6 ns, 220.0 ns] |
| D=5 | la_stack_lu_solve | 44.3 ns [43.9 ns, 44.8 ns] | 46.1 ns [46.0 ns, 46.1 ns] | +4.0% | slower point estimate; marginal CIs separated | 0.96x | 69.0 ns [68.9 ns, 69.2 ns] | 326.8 ns [324.0 ns, 361.6 ns] |
| D=5 | la_stack_norm2_sq | 0.5 ns [0.5 ns, 0.5 ns] | 0.5 ns [0.5 ns, 0.5 ns] | +0.1% | slower point estimate; marginal CIs separated | 1.00x | 0.6 ns [0.6 ns, 0.6 ns] | 4.1 ns [4.1 ns, 4.1 ns] |
| D=5 | la_stack_solve_from_ldlt | 3.9 ns [3.9 ns, 3.9 ns] | 3.9 ns [3.9 ns, 3.9 ns] | +0.9% | slower point estimate; marginal CIs separated | 0.99x | 8.9 ns [8.8 ns, 8.9 ns] | 63.4 ns [63.3 ns, 63.6 ns] |
| D=5 | la_stack_solve_from_lu | 6.0 ns [6.0 ns, 6.0 ns] | 6.0 ns [6.0 ns, 6.0 ns] | -0.7% | faster point estimate; marginal CIs separated | 1.01x | 8.8 ns [8.8 ns, 8.8 ns] | 93.6 ns [91.6 ns, 95.0 ns] |
| D=64 | la_stack_det | 15.21 µs [15.17 µs, 15.25 µs] | 15.11 µs [15.05 µs, 15.16 µs] | -0.7% | faster point estimate; marginal CIs separated | 1.01x | — | — |
| D=64 | la_stack_det_from_ldlt | 22.6 ns [22.5 ns, 22.6 ns] | 22.8 ns [22.8 ns, 22.8 ns] | +1.0% | slower point estimate; marginal CIs separated | 0.99x | 8.2 ns [8.2 ns, 8.3 ns] | 20.2 ns [20.2 ns, 20.3 ns] |
| D=64 | la_stack_det_from_lu | 23.0 ns [22.9 ns, 23.0 ns] | 22.9 ns [22.9 ns, 23.0 ns] | -0.2% | marginal CIs overlap | 1.00x | 8.5 ns [8.5 ns, 8.6 ns] | 21.1 ns [21.1 ns, 21.2 ns] |
| D=64 | la_stack_det_via_lu | 15.21 µs [15.17 µs, 15.23 µs] | 15.03 µs [14.99 µs, 15.06 µs] | -1.1% | faster point estimate; marginal CIs separated | 1.01x | 13.41 µs [13.40 µs, 13.43 µs] | 10.55 µs [10.55 µs, 10.56 µs] |
| D=64 | la_stack_dot | 10.7 ns [10.6 ns, 10.7 ns] | 10.8 ns [10.8 ns, 10.8 ns] | +1.4% | slower point estimate; marginal CIs separated | 0.99x | 8.9 ns [8.9 ns, 8.9 ns] | 8.3 ns [8.3 ns, 8.3 ns] |
| D=64 | la_stack_inf_norm | 612.7 ns [612.2 ns, 614.1 ns] | 611.1 ns [610.2 ns, 613.2 ns] | -0.3% | marginal CIs overlap | 1.00x | 1.09 µs [1.08 µs, 1.09 µs] | 1.54 µs [1.54 µs, 1.54 µs] |
| D=64 | la_stack_ldlt | 19.64 µs [19.60 µs, 19.65 µs] | 20.87 µs [20.81 µs, 20.91 µs] | +6.3% | slower point estimate; marginal CIs separated | 0.94x | 11.48 µs [11.46 µs, 11.50 µs] | 8.81 µs [8.79 µs, 8.82 µs] |
| D=64 | la_stack_ldlt_solve | 21.85 µs [21.79 µs, 21.90 µs] | 23.16 µs [23.08 µs, 23.26 µs] | +6.0% | slower point estimate; marginal CIs separated | 0.94x | 12.12 µs [12.09 µs, 12.14 µs] | 10.15 µs [10.13 µs, 10.16 µs] |
| D=64 | la_stack_lu | 15.02 µs [14.99 µs, 15.03 µs] | 14.42 µs [14.40 µs, 14.44 µs] | -4.0% | faster point estimate; marginal CIs separated | 1.04x | 12.93 µs [12.91 µs, 12.94 µs] | 10.41 µs [10.40 µs, 10.43 µs] |
| D=64 | la_stack_lu_solve | 17.22 µs [17.16 µs, 17.44 µs] | 17.52 µs [17.36 µs, 17.56 µs] | +1.7% | marginal CIs overlap | 0.98x | 14.31 µs [14.30 µs, 14.35 µs] | 12.11 µs [12.09 µs, 12.12 µs] |
| D=64 | la_stack_norm2_sq | 10.6 ns [10.6 ns, 10.7 ns] | 10.8 ns [10.7 ns, 10.8 ns] | +1.3% | slower point estimate; marginal CIs separated | 0.99x | 7.5 ns [7.5 ns, 7.5 ns] | 6.2 ns [6.2 ns, 6.2 ns] |
| D=64 | la_stack_solve_from_ldlt | 1.46 µs [1.45 µs, 1.46 µs] | 1.08 µs [1.07 µs, 1.08 µs] | -26.1% | faster point estimate; marginal CIs separated | 1.35x | 1.25 µs [1.24 µs, 1.25 µs] | 1.26 µs [1.26 µs, 1.27 µs] |
| D=64 | la_stack_solve_from_lu | 2.52 µs [2.51 µs, 2.52 µs] | 2.47 µs [2.47 µs, 2.48 µs] | -1.7% | faster point estimate; marginal CIs separated | 1.02x | 793.3 ns [792.9 ns, 793.9 ns] | 1.67 µs [1.67 µs, 1.68 µs] |
| D=8 | la_stack_det | 90.7 ns [90.5 ns, 91.0 ns] | 90.5 ns [90.3 ns, 90.6 ns] | -0.3% | marginal CIs overlap | 1.00x | — | — |
| D=8 | la_stack_det_from_ldlt | 1.3 ns [1.3 ns, 1.3 ns] | 1.3 ns [1.3 ns, 1.3 ns] | -0.0% | marginal CIs overlap | 1.00x | 0.9 ns [0.9 ns, 0.9 ns] | 1.9 ns [1.9 ns, 1.9 ns] |
| D=8 | la_stack_det_from_ldlt_balanced_range | 8.9 ns [8.9 ns, 8.9 ns] | 8.9 ns [8.9 ns, 9.0 ns] | +0.4% | slower point estimate; marginal CIs separated | 1.00x | — | — |
| D=8 | la_stack_det_from_lu | 1.5 ns [1.5 ns, 1.5 ns] | 1.5 ns [1.5 ns, 1.5 ns] | +0.3% | slower point estimate; marginal CIs separated | 1.00x | 1.0 ns [1.0 ns, 1.0 ns] | 2.2 ns [2.2 ns, 2.2 ns] |
| D=8 | la_stack_det_from_lu_balanced_range | 9.7 ns [9.7 ns, 9.8 ns] | 9.1 ns [9.0 ns, 9.2 ns] | -6.4% | faster point estimate; marginal CIs separated | 1.07x | — | — |
| D=8 | la_stack_det_via_lu | 84.1 ns [83.9 ns, 84.4 ns] | 84.6 ns [84.3 ns, 85.0 ns] | +0.6% | marginal CIs overlap | 0.99x | 144.2 ns [143.9 ns, 144.6 ns] | 294.0 ns [292.1 ns, 294.7 ns] |
| D=8 | la_stack_dot | 0.9 ns [0.9 ns, 0.9 ns] | 1.0 ns [1.0 ns, 1.0 ns] | +3.2% | slower point estimate; marginal CIs separated | 0.97x | 1.1 ns [1.1 ns, 1.1 ns] | 2.6 ns [2.6 ns, 2.7 ns] |
| D=8 | la_stack_inf_norm | 8.4 ns [8.4 ns, 8.4 ns] | 8.4 ns [8.4 ns, 8.4 ns] | -0.1% | marginal CIs overlap | 1.00x | 8.0 ns [8.0 ns, 8.0 ns] | 8.0 ns [8.0 ns, 8.1 ns] |
| D=8 | la_stack_ldlt | 90.7 ns [89.7 ns, 92.3 ns] | 90.5 ns [88.6 ns, 91.4 ns] | -0.2% | marginal CIs overlap | 1.00x | 97.7 ns [97.3 ns, 98.1 ns] | 215.2 ns [214.5 ns, 216.1 ns] |
| D=8 | la_stack_ldlt_ill_conditioned | 91.2 ns [89.9 ns, 92.6 ns] | 90.3 ns [88.8 ns, 91.4 ns] | -1.0% | marginal CIs overlap | 1.01x | — | — |
| D=8 | la_stack_ldlt_solve | 100.9 ns [100.6 ns, 101.1 ns] | 99.3 ns [98.9 ns, 99.5 ns] | -1.6% | faster point estimate; marginal CIs separated | 1.02x | 142.0 ns [134.5 ns, 144.4 ns] | 296.7 ns [291.9 ns, 300.7 ns] |
| D=8 | la_stack_lu | 82.2 ns [82.0 ns, 82.4 ns] | 82.9 ns [82.6 ns, 83.3 ns] | +0.9% | slower point estimate; marginal CIs separated | 0.99x | 144.0 ns [143.7 ns, 144.5 ns] | 275.8 ns [273.4 ns, 277.3 ns] |
| D=8 | la_stack_lu_ill_conditioned | 81.9 ns [81.5 ns, 82.1 ns] | 82.3 ns [82.2 ns, 82.6 ns] | +0.6% | slower point estimate; marginal CIs separated | 0.99x | — | — |
| D=8 | la_stack_lu_pivoting | 90.8 ns [90.6 ns, 90.9 ns] | 90.4 ns [90.2 ns, 90.7 ns] | -0.4% | marginal CIs overlap | 1.00x | — | — |
| D=8 | la_stack_lu_solve | 137.4 ns [135.7 ns, 143.1 ns] | 128.4 ns [128.2 ns, 128.7 ns] | -6.5% | faster point estimate; marginal CIs separated | 1.07x | 165.9 ns [164.8 ns, 168.9 ns] | 401.4 ns [396.0 ns, 407.1 ns] |
| D=8 | la_stack_norm2_sq | 0.7 ns [0.7 ns, 0.7 ns] | 0.7 ns [0.7 ns, 0.7 ns] | +0.4% | marginal CIs overlap | 1.00x | 0.7 ns [0.7 ns, 0.7 ns] | 4.2 ns [4.2 ns, 4.2 ns] |
| D=8 | la_stack_solve_from_ldlt | 8.1 ns [8.1 ns, 8.1 ns] | 8.0 ns [8.0 ns, 8.1 ns] | -0.4% | faster point estimate; marginal CIs separated | 1.00x | 21.0 ns [21.0 ns, 21.0 ns] | 71.1 ns [70.8 ns, 71.1 ns] |
| D=8 | la_stack_solve_from_lu | 13.3 ns [13.3 ns, 13.3 ns] | 13.5 ns [13.5 ns, 13.5 ns] | +1.2% | slower point estimate; marginal CIs separated | 0.99x | 16.0 ns [15.9 ns, 16.0 ns] | 99.5 ns [99.0 ns, 99.6 ns] |

## Coverage Notes

One-sided rows retain the available measurement but are excluded from point-estimate change and ratio calculations.

| Benchmark | Coverage | v0.4.5 (point + CI) | Latest (point + CI) | Note |
|:----------|:---------|-----------------------------:|--------------------:|:-----|
| rational_input_d2/det_big_rational_gaussian | current-only | — | 1.18 µs [1.18 µs, 1.18 µs] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d2/det_row_cleared_bareiss | current-only | — | 191.3 ns [190.9 ns, 192.0 ns] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d2/det_sign_row_cleared_bareiss | current-only | — | 82.3 ns [82.2 ns, 82.4 ns] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d2/solve_big_rational_gaussian | current-only | — | 2.22 µs [2.22 µs, 2.22 µs] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d2/solve_row_cleared_bareiss | current-only | — | 967.5 ns [966.5 ns, 968.2 ns] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d3/det_big_rational_gaussian | current-only | — | 4.89 µs [4.88 µs, 4.90 µs] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d3/det_row_cleared_bareiss | current-only | — | 578.6 ns [578.0 ns, 579.1 ns] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d3/det_sign_row_cleared_bareiss | current-only | — | 368.6 ns [367.8 ns, 369.6 ns] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d3/solve_big_rational_gaussian | current-only | — | 8.44 µs [8.43 µs, 8.46 µs] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d3/solve_row_cleared_bareiss | current-only | — | 2.50 µs [2.49 µs, 2.51 µs] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d4/det_big_rational_gaussian | current-only | — | 14.16 µs [14.14 µs, 14.19 µs] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d4/det_row_cleared_bareiss | current-only | — | 1.17 µs [1.16 µs, 1.17 µs] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d4/det_sign_row_cleared_bareiss | current-only | — | 769.1 ns [767.2 ns, 770.7 ns] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d4/solve_big_rational_gaussian | current-only | — | 24.94 µs [24.90 µs, 24.98 µs] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d4/solve_row_cleared_bareiss | current-only | — | 5.51 µs [5.51 µs, 5.52 µs] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d5/det_big_rational_gaussian | current-only | — | 24.97 µs [24.55 µs, 25.75 µs] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d5/det_row_cleared_bareiss | current-only | — | 2.15 µs [2.15 µs, 2.16 µs] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d5/det_sign_row_cleared_bareiss | current-only | — | 1.58 µs [1.58 µs, 1.59 µs] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d5/solve_big_rational_gaussian | current-only | — | 41.11 µs [41.07 µs, 41.14 µs] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d5/solve_row_cleared_bareiss | current-only | — | 9.67 µs [9.64 µs, 9.72 µs] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d6/det_big_rational_gaussian | current-only | — | 57.91 µs [57.84 µs, 58.02 µs] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d6/det_row_cleared_bareiss | current-only | — | 3.75 µs [3.75 µs, 3.76 µs] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d6/det_sign_row_cleared_bareiss | current-only | — | 3.10 µs [3.09 µs, 3.10 µs] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d6/solve_big_rational_gaussian | current-only | — | 92.50 µs [92.35 µs, 92.60 µs] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d6/solve_row_cleared_bareiss | current-only | — | 18.03 µs [18.01 µs, 18.05 µs] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d7/det_big_rational_gaussian | current-only | — | 112.35 µs [112.25 µs, 112.50 µs] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d7/det_row_cleared_bareiss | current-only | — | 6.15 µs [6.14 µs, 6.16 µs] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d7/det_sign_row_cleared_bareiss | current-only | — | 5.30 µs [5.29 µs, 5.31 µs] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d7/solve_big_rational_gaussian | current-only | — | 165.07 µs [164.88 µs, 165.48 µs] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d7/solve_row_cleared_bareiss | current-only | — | 29.46 µs [29.44 µs, 29.50 µs] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d8/det_big_rational_gaussian | current-only | — | 194.98 µs [193.24 µs, 197.90 µs] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d8/det_row_cleared_bareiss | current-only | — | 9.31 µs [9.29 µs, 9.33 µs] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d8/det_sign_row_cleared_bareiss | current-only | — | 8.47 µs [8.46 µs, 8.48 µs] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d8/solve_big_rational_gaussian | current-only | — | 272.40 µs [272.07 µs, 272.65 µs] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |
| rational_input_d8/solve_row_cleared_bareiss | current-only | — | 43.37 µs [43.31 µs, 43.43 µs] | Baseline v0.4.5 has no correctness-compatible benchmark row under la_stack_pre_rational_input_api. |

## How to Update

Local performance reports are generated in isolated temporary worktrees:

```bash
# Local development: compare the current tree with the latest release
just performance-local

# Release PR: update docs/performance.md and archive the previous report
just performance-release

# Build release docs from retained CSV/JSON inputs (no benchmarks)
just performance-doc

# GitHub Actions release assets
just performance-github-assets

# Explicit repair
just performance-release <current-tag> <previous-tag>
```

`just performance-local` writes `performance.md` plus retained `performance.csv` and
`performance.provenance.json` comparison inputs under `target/bench-reports/` without promoting documentation.
It applies staged and unstaged tracked changes; untracked files are excluded.
`just performance-github-assets` writes `target/bench-reports/github-assets-performance.md`.
`just performance-release` also preserves every recorded local benchmark summary and the report inputs under
`docs/performance/<current>-vs-<baseline>/<run-digest>/`, then promotes distinct-release documentation.
`just performance-doc` consumes the retained pair from either workflow without benchmarking and promotes it when the package versions differ.
After `just clean`, `performance-doc` and `performance-readme` use the latest complete snapshot in `docs/performance/`.
For a distinct pair, `performance-local` followed by `performance-doc` is equivalent to the atomic `performance-release` workflow.

See the [local summary index](https://github.com/acgetchell/la-stack/tree/main/docs/performance)
for the saved-data schema, provenance, and historical availability.
Older curated release-to-release reports are archived in `docs/archive/performance/`.

See `docs/BENCHMARKING.md` for the full comparison workflow.
