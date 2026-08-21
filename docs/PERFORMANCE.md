# Benchmark Performance

**la-stack** v0.4.5 · `99d3392` (HEAD)
**Source revision timestamp**: 2026-08-20 21:16:55 UTC (deterministic report metadata; not the benchmark measurement time)
**Benchmark measurement timestamp**: not recorded by Criterion; use the provenance below to identify the measured revisions and environment.
**Statistic**: median
**Suite**: all
**Scope**: release-signal

## Benchmark Results

Comparison against baseline **v0.4.4**:

Negative point-estimate change means the current point estimate is smaller; a baseline/current point-estimate ratio above 1.00 has the same meaning.
The CI-relation column reports only whether the two marginal Criterion intervals overlap. These are not paired confidence intervals
for the change, so the report makes no statistical-significance or performance-improvement claim from interval separation.

### Reproducibility Provenance

**Measurement environment**: recorded for both samples under one shared current harness.

- CPU: `Apple M4 Max (arm64)`
- OS: `Darwin 25.6.0 arm64`
- rustc: `rustc 1.98.0 (88d9e12ae 2026-08-18)`
- Current commit: `99d33927e389018c599c009603a7eafffcc25089`
- Current Git clean: `false`
- Current source-state SHA-256: `7d72086128b12e2a89d5ca78f0e8cc0cb393a44f1062d027a069722ff9dae007`
- Baseline commit: `03a6dc751b7bf7c69b833aeb4e20e2acb6da2e4c`
- Baseline Git clean: `false`
- Baseline source-state SHA-256: `227a2780c989352be2eacf3f32ca67032cc69fc843e596acc386a7c4f7d459f0`
- Cargo.lock SHA-256: `72b99b4f7f3917d668bb71448198d5d496334684097eac282c859d08d8cf4492`
- Benchmark harness SHA-256: `164d6bbac2e6cb19b85c81677e4ab4428c3abb130dcfd38a2f5f162a5dd8a905`

**Publication and validation environment**:

- Publication CPU: `Apple M4 Max (arm64)`
- Publication OS: `Darwin 25.6.0 arm64`
- Publication rustc: `rustc 1.98.0 (88d9e12ae 2026-08-18)`
- Publication commit: `99d33927e389018c599c009603a7eafffcc25089`
- Publication Git clean: `false`
- Publication source-state SHA-256: `7d72086128b12e2a89d5ca78f0e8cc0cb393a44f1062d027a069722ff9dae007`
- Publication Cargo.lock SHA-256: `72b99b4f7f3917d668bb71448198d5d496334684097eac282c859d08d8cf4492`
- Publication harness SHA-256: `164d6bbac2e6cb19b85c81677e4ab4428c3abb130dcfd38a2f5f162a5dd8a905`
- Criterion suite/scope: `all` / `release-signal`
- Criterion statistic/sample: `median` / `new`
- Criterion dependency version: `0.8.2`
- Baseline command: `just bench-save-baseline v0.4.4`
- Current command: `just bench-latest`
- Correctness gate: `just test-bench-inputs` passed against both the current and baseline revisions using the shared current fixture harness.
- Validated current revision: `99d33927e389018c599c009603a7eafffcc25089` (Git clean: `false`;
  source-state SHA-256: `7d72086128b12e2a89d5ca78f0e8cc0cb393a44f1062d027a069722ff9dae007`)
- Validated baseline revision: `03a6dc751b7bf7c69b833aeb4e20e2acb6da2e4c` (Git clean: `false`;
  source-state SHA-256: `227a2780c989352be2eacf3f32ca67032cc69fc843e596acc386a7c4f7d459f0`)

## Exact arithmetic

| Case | Benchmark | v0.4.4 (point + CI) | Latest (point + CI) | Point-estimate change | CI relation | Point-estimate ratio |
|:-----|:----------|-------:|-------:|-------:|:-----------|--------:|
| D=2 | det | 0.4 ns [0.4 ns, 0.4 ns] | 0.4 ns [0.4 ns, 0.4 ns] | -1.4% | faster point estimate; marginal CIs separated | 1.01x |
| D=2 | det_direct | 0.4 ns [0.4 ns, 0.4 ns] | 0.4 ns [0.4 ns, 0.4 ns] | -1.5% | marginal CIs overlap | 1.02x |
| D=2 | det_direct_with_errbound | 1.7 ns [1.7 ns, 1.7 ns] | 1.7 ns [1.7 ns, 1.7 ns] | -0.7% | faster point estimate; marginal CIs separated | 1.01x |
| D=2 | det_errbound | 1.6 ns [1.6 ns, 1.6 ns] | 1.6 ns [1.6 ns, 1.6 ns] | -0.9% | faster point estimate; marginal CIs separated | 1.01x |
| D=2 | det_exact | 91.7 ns [91.6 ns, 91.9 ns] | 91.3 ns [91.2 ns, 91.5 ns] | -0.4% | faster point estimate; marginal CIs separated | 1.00x |
| D=2 | det_exact_f64_result | 65.6 ns [65.3 ns, 65.8 ns] | 64.7 ns [64.6 ns, 64.7 ns] | -1.3% | faster point estimate; marginal CIs separated | 1.01x |
| D=2 | det_exact_rounded_f64 | 66.7 ns [66.6 ns, 66.8 ns] | 66.6 ns [66.5 ns, 66.7 ns] | -0.1% | marginal CIs overlap | 1.00x |
| D=2 | det_sign_exact | 2.8 ns [2.7 ns, 2.8 ns] | 2.7 ns [2.7 ns, 2.7 ns] | -1.1% | faster point estimate; marginal CIs separated | 1.01x |
| D=2 | solve_exact | 6.90 µs [6.88 µs, 6.91 µs] | 6.78 µs [6.77 µs, 6.79 µs] | -1.8% | faster point estimate; marginal CIs separated | 1.02x |
| D=2 | solve_exact_f64_result | 8.14 µs [8.13 µs, 8.16 µs] | 7.91 µs [7.89 µs, 7.92 µs] | -2.9% | faster point estimate; marginal CIs separated | 1.03x |
| D=2 | solve_exact_rounded_f64 | 7.25 µs [7.23 µs, 7.27 µs] | 7.10 µs [7.10 µs, 7.11 µs] | -2.1% | faster point estimate; marginal CIs separated | 1.02x |
| D=3 | det | 0.8 ns [0.8 ns, 0.8 ns] | 0.8 ns [0.8 ns, 0.8 ns] | -1.1% | faster point estimate; marginal CIs separated | 1.01x |
| D=3 | det_direct | 0.8 ns [0.8 ns, 0.8 ns] | 0.8 ns [0.8 ns, 0.8 ns] | -0.9% | faster point estimate; marginal CIs separated | 1.01x |
| D=3 | det_direct_with_errbound | 3.3 ns [3.3 ns, 3.4 ns] | 3.3 ns [3.3 ns, 3.3 ns] | -0.0% | marginal CIs overlap | 1.00x |
| D=3 | det_errbound | 3.3 ns [3.3 ns, 3.3 ns] | 3.3 ns [3.3 ns, 3.3 ns] | -0.6% | faster point estimate; marginal CIs separated | 1.01x |
| D=3 | det_exact | 341.9 ns [340.5 ns, 343.5 ns] | 323.7 ns [323.2 ns, 324.3 ns] | -5.3% | faster point estimate; marginal CIs separated | 1.06x |
| D=3 | det_exact_f64_result | 300.9 ns [300.5 ns, 301.6 ns] | 293.7 ns [293.2 ns, 293.9 ns] | -2.4% | faster point estimate; marginal CIs separated | 1.02x |
| D=3 | det_exact_rounded_f64 | 300.5 ns [300.1 ns, 301.0 ns] | 293.7 ns [293.5 ns, 294.0 ns] | -2.2% | faster point estimate; marginal CIs separated | 1.02x |
| D=3 | det_sign_exact | 4.7 ns [4.7 ns, 4.7 ns] | 4.7 ns [4.7 ns, 4.7 ns] | -0.6% | faster point estimate; marginal CIs separated | 1.01x |
| D=3 | solve_exact | 30.23 µs [30.12 µs, 30.31 µs] | 29.93 µs [29.89 µs, 29.97 µs] | -1.0% | faster point estimate; marginal CIs separated | 1.01x |
| D=3 | solve_exact_f64_result | 32.20 µs [32.13 µs, 32.26 µs] | 31.88 µs [31.83 µs, 31.95 µs] | -1.0% | faster point estimate; marginal CIs separated | 1.01x |
| D=3 | solve_exact_rounded_f64 | 31.46 µs [31.19 µs, 31.53 µs] | 30.49 µs [30.44 µs, 30.52 µs] | -3.1% | faster point estimate; marginal CIs separated | 1.03x |
| D=4 | det | 3.3 ns [3.3 ns, 3.3 ns] | 3.2 ns [3.2 ns, 3.2 ns] | -2.1% | faster point estimate; marginal CIs separated | 1.02x |
| D=4 | det_direct | 2.4 ns [2.4 ns, 2.4 ns] | 2.4 ns [2.4 ns, 2.4 ns] | -1.4% | faster point estimate; marginal CIs separated | 1.01x |
| D=4 | det_direct_with_errbound | 6.9 ns [6.9 ns, 6.9 ns] | 6.8 ns [6.8 ns, 6.8 ns] | -1.4% | faster point estimate; marginal CIs separated | 1.01x |
| D=4 | det_errbound | 6.8 ns [6.8 ns, 6.8 ns] | 6.7 ns [6.7 ns, 6.7 ns] | -1.3% | faster point estimate; marginal CIs separated | 1.01x |
| D=4 | det_exact | 1.02 µs [1.02 µs, 1.03 µs] | 1.01 µs [1.01 µs, 1.01 µs] | -1.1% | faster point estimate; marginal CIs separated | 1.01x |
| D=4 | det_exact_f64_result | 993.0 ns [990.5 ns, 996.0 ns] | 976.7 ns [974.5 ns, 978.7 ns] | -1.6% | faster point estimate; marginal CIs separated | 1.02x |
| D=4 | det_exact_rounded_f64 | 991.3 ns [989.1 ns, 993.3 ns] | 973.8 ns [973.0 ns, 976.3 ns] | -1.8% | faster point estimate; marginal CIs separated | 1.02x |
| D=4 | det_sign_exact | 7.8 ns [7.8 ns, 7.8 ns] | 7.7 ns [7.7 ns, 7.7 ns] | -1.2% | faster point estimate; marginal CIs separated | 1.01x |
| D=4 | solve_exact | 78.81 µs [78.73 µs, 78.97 µs] | 77.91 µs [77.80 µs, 77.99 µs] | -1.1% | faster point estimate; marginal CIs separated | 1.01x |
| D=4 | solve_exact_f64_result | 81.82 µs [81.72 µs, 81.92 µs] | 80.71 µs [80.62 µs, 80.79 µs] | -1.4% | faster point estimate; marginal CIs separated | 1.01x |
| D=4 | solve_exact_rounded_f64 | 79.66 µs [79.57 µs, 79.83 µs] | 78.66 µs [78.57 µs, 78.76 µs] | -1.3% | faster point estimate; marginal CIs separated | 1.01x |
| D=5 | det | 25.1 ns [25.1 ns, 25.2 ns] | 25.0 ns [24.9 ns, 25.1 ns] | -0.6% | marginal CIs overlap | 1.01x |
| D=5 | det_exact | 3.08 µs [3.06 µs, 3.12 µs] | 2.94 µs [2.94 µs, 2.95 µs] | -4.4% | faster point estimate; marginal CIs separated | 1.05x |
| D=5 | det_exact_f64_result | 3.02 µs [3.01 µs, 3.03 µs] | 2.91 µs [2.90 µs, 2.91 µs] | -3.7% | faster point estimate; marginal CIs separated | 1.04x |
| D=5 | det_exact_rounded_f64 | 2.98 µs [2.98 µs, 2.99 µs] | 2.96 µs [2.95 µs, 2.97 µs] | -0.6% | faster point estimate; marginal CIs separated | 1.01x |
| D=5 | det_sign_exact | 3.15 µs [3.14 µs, 3.17 µs] | 3.05 µs [3.04 µs, 3.05 µs] | -3.4% | faster point estimate; marginal CIs separated | 1.04x |
| D=5 | solve_exact | 156.87 µs [156.58 µs, 157.10 µs] | 155.34 µs [155.10 µs, 155.58 µs] | -1.0% | faster point estimate; marginal CIs separated | 1.01x |
| D=5 | solve_exact_f64_result | 159.75 µs [159.55 µs, 160.17 µs] | 158.29 µs [158.09 µs, 158.43 µs] | -0.9% | faster point estimate; marginal CIs separated | 1.01x |
| D=5 | solve_exact_rounded_f64 | 159.08 µs [158.45 µs, 160.00 µs] | 156.37 µs [156.22 µs, 156.57 µs] | -1.7% | faster point estimate; marginal CIs separated | 1.02x |
| Hilbert 4x4 | det_exact | 1.13 µs [1.12 µs, 1.13 µs] | 1.10 µs [1.09 µs, 1.10 µs] | -2.8% | faster point estimate; marginal CIs separated | 1.03x |
| Hilbert 4x4 | det_sign_exact | 8.0 ns [7.9 ns, 8.0 ns] | 7.7 ns [7.7 ns, 7.7 ns] | -3.2% | faster point estimate; marginal CIs separated | 1.03x |
| Hilbert 4x4 | solve_exact | 59.49 µs [59.37 µs, 59.59 µs] | 58.24 µs [58.12 µs, 58.35 µs] | -2.1% | faster point estimate; marginal CIs separated | 1.02x |
| Hilbert 4x4 | solve_exact_f64_result | 61.76 µs [61.70 µs, 61.84 µs] | 60.89 µs [60.80 µs, 60.95 µs] | -1.4% | faster point estimate; marginal CIs separated | 1.01x |
| Hilbert 4x4 | solve_exact_rounded_f64 | 60.17 µs [60.12 µs, 60.28 µs] | 58.88 µs [58.75 µs, 59.03 µs] | -2.1% | faster point estimate; marginal CIs separated | 1.02x |
| Hilbert 5x5 | det_exact | 3.13 µs [3.12 µs, 3.13 µs] | 3.05 µs [3.03 µs, 3.05 µs] | -2.6% | faster point estimate; marginal CIs separated | 1.03x |
| Hilbert 5x5 | det_sign_exact | 3.33 µs [3.32 µs, 3.34 µs] | 3.27 µs [3.26 µs, 3.28 µs] | -1.9% | faster point estimate; marginal CIs separated | 1.02x |
| Hilbert 5x5 | solve_exact | 122.74 µs [122.59 µs, 122.86 µs] | 120.11 µs [119.87 µs, 120.28 µs] | -2.1% | faster point estimate; marginal CIs separated | 1.02x |
| Hilbert 5x5 | solve_exact_f64_result | 125.85 µs [125.61 µs, 126.01 µs] | 123.42 µs [123.24 µs, 123.55 µs] | -1.9% | faster point estimate; marginal CIs separated | 1.02x |
| Hilbert 5x5 | solve_exact_rounded_f64 | 124.06 µs [123.87 µs, 124.28 µs] | 122.12 µs [121.79 µs, 122.49 µs] | -1.6% | faster point estimate; marginal CIs separated | 1.02x |
| Large entries 3x3 | det_exact | 309.9 ns [308.3 ns, 310.7 ns] | 290.4 ns [290.0 ns, 291.4 ns] | -6.3% | faster point estimate; marginal CIs separated | 1.07x |
| Large entries 3x3 | det_sign_exact | 291.0 ns [289.7 ns, 292.2 ns] | 282.1 ns [281.0 ns, 283.3 ns] | -3.0% | faster point estimate; marginal CIs separated | 1.03x |
| Large entries 3x3 | solve_exact | 95.06 µs [94.85 µs, 95.29 µs] | 92.60 µs [92.46 µs, 92.76 µs] | -2.6% | faster point estimate; marginal CIs separated | 1.03x |
| Large entries 3x3 | solve_exact_f64_result | 96.49 µs [96.09 µs, 97.15 µs] | 93.61 µs [93.42 µs, 93.81 µs] | -3.0% | faster point estimate; marginal CIs separated | 1.03x |
| Large entries 3x3 | solve_exact_rounded_f64 | 95.16 µs [94.97 µs, 95.34 µs] | 93.36 µs [93.20 µs, 93.55 µs] | -1.9% | faster point estimate; marginal CIs separated | 1.02x |
| Near-singular 3x3 | det_exact | 322.6 ns [321.3 ns, 325.5 ns] | 303.0 ns [302.2 ns, 303.5 ns] | -6.1% | faster point estimate; marginal CIs separated | 1.06x |
| Near-singular 3x3 | det_sign_exact | 372.6 ns [370.5 ns, 374.2 ns] | 345.8 ns [344.2 ns, 347.4 ns] | -7.2% | faster point estimate; marginal CIs separated | 1.08x |
| Near-singular 3x3 | solve_exact | 2.42 µs [2.41 µs, 2.43 µs] | 2.34 µs [2.34 µs, 2.35 µs] | -3.0% | faster point estimate; marginal CIs separated | 1.03x |
| Near-singular 3x3 | solve_exact_f64_result | 2.51 µs [2.49 µs, 2.52 µs] | 2.46 µs [2.45 µs, 2.47 µs] | -1.9% | faster point estimate; marginal CIs separated | 1.02x |
| Near-singular 3x3 | solve_exact_rounded_f64 | 2.42 µs [2.42 µs, 2.43 µs] | 2.35 µs [2.35 µs, 2.36 µs] | -2.9% | faster point estimate; marginal CIs separated | 1.03x |
| Random corpus D=2 | det_exact | 3.20 µs [3.19 µs, 3.20 µs] | 3.13 µs [3.13 µs, 3.13 µs] | -2.0% | faster point estimate; marginal CIs separated | 1.02x |
| Random corpus D=2 | det_sign_exact | 141.6 ns [140.8 ns, 143.3 ns] | 136.5 ns [135.7 ns, 137.2 ns] | -3.6% | faster point estimate; marginal CIs separated | 1.04x |
| Random corpus D=2 | solve_exact | 69.25 µs [69.09 µs, 69.50 µs] | 68.00 µs [67.75 µs, 68.99 µs] | -1.8% | faster point estimate; marginal CIs separated | 1.02x |
| Random corpus D=2 | solve_exact_f64_result | 80.02 µs [79.86 µs, 80.24 µs] | 78.04 µs [77.88 µs, 78.24 µs] | -2.5% | faster point estimate; marginal CIs separated | 1.03x |
| Random corpus D=2 | solve_exact_rounded_f64 | 70.51 µs [70.05 µs, 70.90 µs] | 67.94 µs [67.89 µs, 68.00 µs] | -3.6% | faster point estimate; marginal CIs separated | 1.04x |
| Random corpus D=3 | det_exact | 9.03 µs [8.95 µs, 9.10 µs] | 8.81 µs [8.24 µs, 8.89 µs] | -2.5% | faster point estimate; marginal CIs separated | 1.03x |
| Random corpus D=3 | det_sign_exact | 228.6 ns [228.3 ns, 228.8 ns] | 225.0 ns [224.8 ns, 225.1 ns] | -1.6% | faster point estimate; marginal CIs separated | 1.02x |
| Random corpus D=3 | solve_exact | 219.99 µs [219.73 µs, 220.21 µs] | 216.21 µs [215.83 µs, 216.36 µs] | -1.7% | faster point estimate; marginal CIs separated | 1.02x |
| Random corpus D=3 | solve_exact_f64_result | 235.73 µs [234.74 µs, 236.39 µs] | 230.14 µs [229.96 µs, 230.31 µs] | -2.4% | faster point estimate; marginal CIs separated | 1.02x |
| Random corpus D=3 | solve_exact_rounded_f64 | 222.81 µs [222.46 µs, 223.31 µs] | 217.76 µs [217.35 µs, 218.01 µs] | -2.3% | faster point estimate; marginal CIs separated | 1.02x |
| Random corpus D=4 | det_exact | 25.74 µs [25.66 µs, 25.87 µs] | 24.87 µs [24.84 µs, 24.98 µs] | -3.4% | faster point estimate; marginal CIs separated | 1.03x |
| Random corpus D=4 | det_sign_exact | 434.0 ns [432.4 ns, 435.8 ns] | 424.5 ns [423.9 ns, 425.6 ns] | -2.2% | faster point estimate; marginal CIs separated | 1.02x |
| Random corpus D=4 | solve_exact | 511.24 µs [509.70 µs, 513.54 µs] | 493.20 µs [492.65 µs, 493.60 µs] | -3.5% | faster point estimate; marginal CIs separated | 1.04x |
| Random corpus D=4 | solve_exact_f64_result | 531.44 µs [528.39 µs, 535.08 µs] | 513.06 µs [512.52 µs, 513.31 µs] | -3.5% | faster point estimate; marginal CIs separated | 1.04x |
| Random corpus D=4 | solve_exact_rounded_f64 | 502.60 µs [502.00 µs, 503.33 µs] | 494.39 µs [493.41 µs, 495.23 µs] | -1.6% | faster point estimate; marginal CIs separated | 1.02x |
| Random corpus D=5 | det_exact | 47.95 µs [47.47 µs, 48.29 µs] | 46.58 µs [46.42 µs, 46.65 µs] | -2.9% | faster point estimate; marginal CIs separated | 1.03x |
| Random corpus D=5 | det_sign_exact | 48.52 µs [47.83 µs, 48.94 µs] | 51.72 µs [50.78 µs, 52.35 µs] | +6.6% | slower point estimate; marginal CIs separated | 0.94x |
| Random corpus D=5 | solve_exact | 994.36 µs [993.76 µs, 995.50 µs] | 972.74 µs [971.25 µs, 974.13 µs] | -2.2% | faster point estimate; marginal CIs separated | 1.02x |
| Random corpus D=5 | solve_exact_f64_result | 1.02 ms [1.02 ms, 1.02 ms] | 992.23 µs [990.24 µs, 993.69 µs] | -2.9% | faster point estimate; marginal CIs separated | 1.03x |
| Random corpus D=5 | solve_exact_rounded_f64 | 1.03 ms [1.03 ms, 1.04 ms] | 972.11 µs [970.46 µs, 973.74 µs] | -6.0% | faster point estimate; marginal CIs separated | 1.06x |

## vs_linalg

| Case | Benchmark | v0.4.4 (point + CI) | Latest (point + CI) | Point-estimate change | CI relation | Point-estimate ratio | v0.4.4 nalgebra | v0.4.4 faer |
|:-----|:----------|-------:|-------:|-------:|:-----------|--------:|-------:|-------:|
| D=16 | la_stack_det | 432.1 ns [428.1 ns, 433.9 ns] | 436.8 ns [435.3 ns, 438.7 ns] | +1.1% | slower point estimate; marginal CIs separated | 0.99x | — | — |
| D=16 | la_stack_det_from_ldlt | 3.5 ns [3.4 ns, 3.6 ns] | 2.9 ns [2.9 ns, 2.9 ns] | -17.7% | faster point estimate; marginal CIs separated | 1.21x | 1.7 ns [1.7 ns, 1.7 ns] | 4.5 ns [4.5 ns, 4.6 ns] |
| D=16 | la_stack_det_from_lu | 4.2 ns [4.1 ns, 4.3 ns] | 3.3 ns [3.3 ns, 3.3 ns] | -22.4% | faster point estimate; marginal CIs separated | 1.29x | 1.8 ns [1.8 ns, 1.8 ns] | 5.0 ns [4.9 ns, 5.0 ns] |
| D=16 | la_stack_det_via_lu | 423.2 ns [418.2 ns, 427.0 ns] | 422.8 ns [420.8 ns, 424.7 ns] | -0.1% | marginal CIs overlap | 1.00x | 458.2 ns [457.0 ns, 458.6 ns] | 679.3 ns [672.2 ns, 683.9 ns] |
| D=16 | la_stack_dot | 2.4 ns [2.4 ns, 2.4 ns] | 2.3 ns [2.3 ns, 2.3 ns] | -1.0% | faster point estimate; marginal CIs separated | 1.01x | 1.9 ns [1.9 ns, 1.9 ns] | 4.6 ns [4.6 ns, 4.6 ns] |
| D=16 | la_stack_inf_norm | 33.5 ns [33.3 ns, 33.6 ns] | 32.7 ns [32.7 ns, 32.7 ns] | -2.4% | faster point estimate; marginal CIs separated | 1.02x | 33.4 ns [33.2 ns, 33.6 ns] | 33.4 ns [33.3 ns, 33.5 ns] |
| D=16 | la_stack_ldlt | 394.2 ns [392.7 ns, 395.9 ns] | 390.6 ns [389.5 ns, 391.2 ns] | -0.9% | faster point estimate; marginal CIs separated | 1.01x | 408.4 ns [407.3 ns, 410.1 ns] | 446.8 ns [443.5 ns, 448.3 ns] |
| D=16 | la_stack_ldlt_solve | 433.6 ns [432.5 ns, 434.6 ns] | 436.1 ns [435.2 ns, 437.4 ns] | +0.6% | slower point estimate; marginal CIs separated | 0.99x | 638.5 ns [637.9 ns, 640.1 ns] | 624.6 ns [623.7 ns, 626.6 ns] |
| D=16 | la_stack_lu | 398.6 ns [397.3 ns, 400.9 ns] | 403.5 ns [401.3 ns, 405.9 ns] | +1.2% | slower point estimate; marginal CIs separated | 0.99x | 465.6 ns [464.7 ns, 466.9 ns] | 658.8 ns [654.2 ns, 665.1 ns] |
| D=16 | la_stack_lu_solve | 667.8 ns [665.2 ns, 672.5 ns] | 672.5 ns [670.9 ns, 674.0 ns] | +0.7% | marginal CIs overlap | 0.99x | 585.3 ns [584.6 ns, 586.3 ns] | 897.2 ns [895.5 ns, 900.9 ns] |
| D=16 | la_stack_norm2_sq | 2.1 ns [2.1 ns, 2.1 ns] | 2.1 ns [2.1 ns, 2.1 ns] | -3.5% | faster point estimate; marginal CIs separated | 1.04x | 1.5 ns [1.5 ns, 1.5 ns] | 4.2 ns [4.2 ns, 4.2 ns] |
| D=16 | la_stack_solve_from_ldlt | 28.0 ns [28.0 ns, 28.1 ns] | 27.7 ns [27.7 ns, 27.8 ns] | -1.0% | faster point estimate; marginal CIs separated | 1.01x | 126.5 ns [126.1 ns, 127.1 ns] | 179.6 ns [178.8 ns, 180.2 ns] |
| D=16 | la_stack_solve_from_lu | 204.7 ns [204.0 ns, 205.2 ns] | 203.1 ns [202.5 ns, 204.0 ns] | -0.8% | faster point estimate; marginal CIs separated | 1.01x | 95.3 ns [95.2 ns, 95.5 ns] | 245.0 ns [244.1 ns, 247.1 ns] |
| D=2 | la_stack_det | 0.6 ns [0.6 ns, 0.6 ns] | 0.6 ns [0.6 ns, 0.6 ns] | +2.3% | slower point estimate; marginal CIs separated | 0.98x | — | — |
| D=2 | la_stack_det_from_ldlt | 0.5 ns [0.5 ns, 0.5 ns] | 0.5 ns [0.5 ns, 0.5 ns] | -0.4% | faster point estimate; marginal CIs separated | 1.00x | 0.4 ns [0.4 ns, 0.4 ns] | 0.6 ns [0.6 ns, 0.6 ns] |
| D=2 | la_stack_det_from_lu | 0.5 ns [0.5 ns, 0.5 ns] | 0.5 ns [0.5 ns, 0.5 ns] | -0.1% | marginal CIs overlap | 1.00x | 0.5 ns [0.5 ns, 0.5 ns] | 0.8 ns [0.8 ns, 0.8 ns] |
| D=2 | la_stack_det_via_lu | 1.8 ns [1.8 ns, 1.8 ns] | 1.8 ns [1.8 ns, 1.8 ns] | +0.6% | slower point estimate; marginal CIs separated | 0.99x | 0.8 ns [0.8 ns, 0.8 ns] | 98.3 ns [97.9 ns, 99.0 ns] |
| D=2 | la_stack_dot | 0.6 ns [0.6 ns, 0.6 ns] | 0.6 ns [0.6 ns, 0.6 ns] | +5.3% | slower point estimate; marginal CIs separated | 0.95x | 0.6 ns [0.6 ns, 0.6 ns] | 0.7 ns [0.7 ns, 0.7 ns] |
| D=2 | la_stack_inf_norm | 0.7 ns [0.7 ns, 0.7 ns] | 0.6 ns [0.6 ns, 0.6 ns] | -8.9% | faster point estimate; marginal CIs separated | 1.10x | 0.5 ns [0.5 ns, 0.5 ns] | 0.7 ns [0.7 ns, 0.7 ns] |
| D=2 | la_stack_ldlt | 6.8 ns [6.8 ns, 6.9 ns] | 7.1 ns [7.1 ns, 7.1 ns] | +3.7% | slower point estimate; marginal CIs separated | 0.96x | 1.7 ns [1.7 ns, 1.7 ns] | 81.6 ns [81.4 ns, 81.8 ns] |
| D=2 | la_stack_ldlt_solve | 9.9 ns [9.9 ns, 9.9 ns] | 10.1 ns [10.0 ns, 10.1 ns] | +2.1% | slower point estimate; marginal CIs separated | 0.98x | 2.7 ns [2.7 ns, 2.7 ns] | 124.6 ns [124.2 ns, 125.1 ns] |
| D=2 | la_stack_lu | 1.6 ns [1.6 ns, 1.6 ns] | 1.6 ns [1.6 ns, 1.6 ns] | +2.4% | slower point estimate; marginal CIs separated | 0.98x | 1.6 ns [1.6 ns, 1.6 ns] | 95.9 ns [95.2 ns, 96.3 ns] |
| D=2 | la_stack_lu_solve | 2.1 ns [2.0 ns, 2.1 ns] | 2.0 ns [2.0 ns, 2.0 ns] | -0.6% | marginal CIs overlap | 1.01x | 4.6 ns [4.6 ns, 4.6 ns] | 151.9 ns [151.7 ns, 152.5 ns] |
| D=2 | la_stack_norm2_sq | 0.4 ns [0.4 ns, 0.4 ns] | 0.4 ns [0.4 ns, 0.4 ns] | +1.0% | slower point estimate; marginal CIs separated | 0.99x | 0.4 ns [0.4 ns, 0.4 ns] | 4.2 ns [4.2 ns, 4.2 ns] |
| D=2 | la_stack_solve_from_ldlt | 1.2 ns [1.2 ns, 1.2 ns] | 1.3 ns [1.3 ns, 1.3 ns] | +1.3% | slower point estimate; marginal CIs separated | 0.99x | 1.3 ns [1.3 ns, 1.3 ns] | 37.6 ns [37.4 ns, 37.9 ns] |
| D=2 | la_stack_solve_from_lu | 1.3 ns [1.3 ns, 1.3 ns] | 1.3 ns [1.3 ns, 1.3 ns] | +2.6% | slower point estimate; marginal CIs separated | 0.98x | 2.9 ns [2.9 ns, 2.9 ns] | 48.8 ns [48.6 ns, 49.1 ns] |
| D=3 | la_stack_det | 1.4 ns [1.3 ns, 1.4 ns] | 1.3 ns [1.3 ns, 1.3 ns] | -0.8% | faster point estimate; marginal CIs separated | 1.01x | — | — |
| D=3 | la_stack_det_from_ldlt | 0.7 ns [0.7 ns, 0.7 ns] | 0.6 ns [0.6 ns, 0.6 ns] | -19.6% | faster point estimate; marginal CIs separated | 1.24x | 0.5 ns [0.5 ns, 0.5 ns] | 0.8 ns [0.8 ns, 0.8 ns] |
| D=3 | la_stack_det_from_lu | 0.7 ns [0.7 ns, 0.7 ns] | 0.7 ns [0.7 ns, 0.7 ns] | -6.0% | faster point estimate; marginal CIs separated | 1.06x | 0.5 ns [0.5 ns, 0.5 ns] | 1.0 ns [1.0 ns, 1.0 ns] |
| D=3 | la_stack_det_via_lu | 8.8 ns [8.7 ns, 8.8 ns] | 8.8 ns [8.8 ns, 8.9 ns] | +0.7% | slower point estimate; marginal CIs separated | 0.99x | 16.9 ns [16.8 ns, 17.0 ns] | 141.6 ns [140.9 ns, 142.3 ns] |
| D=3 | la_stack_dot | 0.7 ns [0.7 ns, 0.7 ns] | 0.7 ns [0.7 ns, 0.7 ns] | +0.1% | marginal CIs overlap | 1.00x | 0.7 ns [0.7 ns, 0.7 ns] | 0.9 ns [0.9 ns, 0.9 ns] |
| D=3 | la_stack_inf_norm | 1.3 ns [1.3 ns, 1.3 ns] | 1.3 ns [1.3 ns, 1.3 ns] | -0.4% | faster point estimate; marginal CIs separated | 1.00x | 1.1 ns [1.1 ns, 1.1 ns] | 1.3 ns [1.3 ns, 1.3 ns] |
| D=3 | la_stack_ldlt | 14.8 ns [14.5 ns, 15.1 ns] | 14.1 ns [14.1 ns, 14.2 ns] | -4.5% | faster point estimate; marginal CIs separated | 1.05x | 4.0 ns [4.0 ns, 4.0 ns] | 97.8 ns [97.6 ns, 98.2 ns] |
| D=3 | la_stack_ldlt_solve | 11.6 ns [11.6 ns, 11.7 ns] | 11.7 ns [11.7 ns, 11.8 ns] | +0.9% | slower point estimate; marginal CIs separated | 0.99x | 6.0 ns [6.0 ns, 6.0 ns] | 142.4 ns [141.9 ns, 143.0 ns] |
| D=3 | la_stack_lu | 8.5 ns [8.4 ns, 8.5 ns] | 8.4 ns [8.4 ns, 8.4 ns] | -1.2% | faster point estimate; marginal CIs separated | 1.01x | 15.4 ns [15.3 ns, 15.5 ns] | 138.0 ns [137.4 ns, 138.6 ns] |
| D=3 | la_stack_lu_solve | 9.9 ns [9.9 ns, 10.0 ns] | 10.0 ns [10.0 ns, 10.0 ns] | +0.5% | marginal CIs overlap | 1.00x | 23.5 ns [23.4 ns, 23.6 ns] | 196.4 ns [194.9 ns, 198.8 ns] |
| D=3 | la_stack_norm2_sq | 0.5 ns [0.5 ns, 0.5 ns] | 0.4 ns [0.4 ns, 0.4 ns] | -1.5% | faster point estimate; marginal CIs separated | 1.01x | 0.4 ns [0.4 ns, 0.4 ns] | 4.2 ns [4.1 ns, 4.2 ns] |
| D=3 | la_stack_solve_from_ldlt | 1.8 ns [1.8 ns, 1.8 ns] | 1.8 ns [1.8 ns, 1.8 ns] | +0.8% | slower point estimate; marginal CIs separated | 0.99x | 2.9 ns [2.9 ns, 2.9 ns] | 40.6 ns [40.1 ns, 41.0 ns] |
| D=3 | la_stack_solve_from_lu | 2.1 ns [2.1 ns, 2.1 ns] | 2.1 ns [2.1 ns, 2.1 ns] | +0.5% | marginal CIs overlap | 1.00x | 4.6 ns [4.6 ns, 4.6 ns] | 48.6 ns [48.5 ns, 48.7 ns] |
| D=32 | la_stack_det | 2.22 µs [2.21 µs, 2.23 µs] | 2.17 µs [2.17 µs, 2.18 µs] | -2.3% | faster point estimate; marginal CIs separated | 1.02x | — | — |
| D=32 | la_stack_det_from_ldlt | 14.1 ns [14.1 ns, 14.1 ns] | 7.9 ns [7.9 ns, 7.9 ns] | -44.0% | faster point estimate; marginal CIs separated | 1.79x | 3.0 ns [3.0 ns, 3.0 ns] | 8.7 ns [8.7 ns, 8.8 ns] |
| D=32 | la_stack_det_from_lu | 14.1 ns [14.0 ns, 14.1 ns] | 9.1 ns [9.1 ns, 9.1 ns] | -35.2% | faster point estimate; marginal CIs separated | 1.54x | 3.1 ns [3.1 ns, 3.1 ns] | 9.0 ns [8.9 ns, 9.0 ns] |
| D=32 | la_stack_det_via_lu | 2.23 µs [2.21 µs, 2.23 µs] | 2.04 µs [2.03 µs, 2.04 µs] | -8.4% | faster point estimate; marginal CIs separated | 1.09x | 2.53 µs [2.52 µs, 2.54 µs] | 2.37 µs [2.36 µs, 2.38 µs] |
| D=32 | la_stack_dot | 4.1 ns [4.1 ns, 4.1 ns] | 4.0 ns [4.0 ns, 4.0 ns] | -1.2% | faster point estimate; marginal CIs separated | 1.01x | 4.8 ns [4.8 ns, 4.8 ns] | 15.2 ns [15.2 ns, 15.2 ns] |
| D=32 | la_stack_inf_norm | 131.6 ns [131.4 ns, 131.8 ns] | 130.4 ns [130.0 ns, 130.7 ns] | -0.9% | faster point estimate; marginal CIs separated | 1.01x | 157.7 ns [157.2 ns, 158.1 ns] | 167.1 ns [166.9 ns, 167.3 ns] |
| D=32 | la_stack_ldlt | 2.47 µs [2.46 µs, 2.47 µs] | 2.52 µs [2.51 µs, 2.52 µs] | +2.0% | slower point estimate; marginal CIs separated | 0.98x | 2.13 µs [2.12 µs, 2.13 µs] | 1.50 µs [1.49 µs, 1.51 µs] |
| D=32 | la_stack_ldlt_solve | 2.83 µs [2.82 µs, 2.84 µs] | 2.84 µs [2.83 µs, 2.86 µs] | +0.4% | marginal CIs overlap | 1.00x | 2.78 µs [2.77 µs, 2.79 µs] | 1.97 µs [1.97 µs, 1.98 µs] |
| D=32 | la_stack_lu | 2.15 µs [2.14 µs, 2.15 µs] | 2.07 µs [2.07 µs, 2.08 µs] | -3.6% | faster point estimate; marginal CIs separated | 1.04x | 2.21 µs [2.20 µs, 2.22 µs] | 2.30 µs [2.29 µs, 2.30 µs] |
| D=32 | la_stack_lu_solve | 2.87 µs [2.85 µs, 2.89 µs] | 2.78 µs [2.77 µs, 2.78 µs] | -3.3% | faster point estimate; marginal CIs separated | 1.03x | 2.50 µs [2.49 µs, 2.51 µs] | 2.95 µs [2.95 µs, 2.96 µs] |
| D=32 | la_stack_norm2_sq | 4.0 ns [4.0 ns, 4.1 ns] | 4.0 ns [4.0 ns, 4.0 ns] | -1.2% | faster point estimate; marginal CIs separated | 1.01x | 3.9 ns [3.9 ns, 3.9 ns] | 4.3 ns [4.3 ns, 4.3 ns] |
| D=32 | la_stack_solve_from_ldlt | 316.1 ns [315.9 ns, 316.5 ns] | 313.6 ns [313.2 ns, 314.1 ns] | -0.8% | faster point estimate; marginal CIs separated | 1.01x | 568.0 ns [566.1 ns, 569.7 ns] | 466.2 ns [464.8 ns, 467.1 ns] |
| D=32 | la_stack_solve_from_lu | 728.4 ns [724.4 ns, 732.7 ns] | 686.0 ns [684.4 ns, 687.3 ns] | -5.8% | faster point estimate; marginal CIs separated | 1.06x | 335.2 ns [334.7 ns, 336.0 ns] | 620.9 ns [619.3 ns, 621.8 ns] |
| D=4 | la_stack_det | 2.6 ns [2.5 ns, 2.6 ns] | 2.6 ns [2.6 ns, 2.6 ns] | +1.3% | slower point estimate; marginal CIs separated | 0.99x | — | — |
| D=4 | la_stack_det_from_ldlt | 1.0 ns [1.0 ns, 1.0 ns] | 0.8 ns [0.8 ns, 0.8 ns] | -20.3% | faster point estimate; marginal CIs separated | 1.25x | 0.5 ns [0.5 ns, 0.5 ns] | 1.0 ns [1.0 ns, 1.0 ns] |
| D=4 | la_stack_det_from_lu | 1.0 ns [1.0 ns, 1.0 ns] | 0.9 ns [0.9 ns, 0.9 ns] | -13.1% | faster point estimate; marginal CIs separated | 1.15x | 0.6 ns [0.6 ns, 0.6 ns] | 1.3 ns [1.3 ns, 1.3 ns] |
| D=4 | la_stack_det_via_lu | 14.3 ns [14.2 ns, 14.3 ns] | 14.4 ns [14.4 ns, 14.4 ns] | +0.8% | slower point estimate; marginal CIs separated | 0.99x | 31.3 ns [31.2 ns, 31.4 ns] | 172.0 ns [169.9 ns, 177.8 ns] |
| D=4 | la_stack_dot | 0.7 ns [0.7 ns, 0.7 ns] | 0.7 ns [0.7 ns, 0.7 ns] | -0.2% | marginal CIs overlap | 1.00x | 0.6 ns [0.6 ns, 0.6 ns] | 1.2 ns [1.2 ns, 1.2 ns] |
| D=4 | la_stack_inf_norm | 2.3 ns [2.2 ns, 2.3 ns] | 2.3 ns [2.3 ns, 2.3 ns] | +0.3% | slower point estimate; marginal CIs separated | 1.00x | 2.0 ns [2.0 ns, 2.0 ns] | 2.0 ns [2.0 ns, 2.0 ns] |
| D=4 | la_stack_ldlt | 22.2 ns [22.2 ns, 22.2 ns] | 22.0 ns [22.0 ns, 22.1 ns] | -0.7% | faster point estimate; marginal CIs separated | 1.01x | 8.0 ns [8.0 ns, 8.0 ns] | 121.3 ns [121.0 ns, 122.0 ns] |
| D=4 | la_stack_ldlt_solve | 24.1 ns [24.1 ns, 24.2 ns] | 24.0 ns [24.0 ns, 24.1 ns] | -0.4% | marginal CIs overlap | 1.00x | 11.6 ns [11.6 ns, 11.6 ns] | 169.0 ns [168.3 ns, 169.4 ns] |
| D=4 | la_stack_lu | 14.0 ns [13.9 ns, 14.0 ns] | 13.9 ns [13.8 ns, 13.9 ns] | -0.8% | faster point estimate; marginal CIs separated | 1.01x | 30.0 ns [29.9 ns, 30.1 ns] | 165.5 ns [164.6 ns, 168.1 ns] |
| D=4 | la_stack_lu_solve | 22.0 ns [21.9 ns, 22.0 ns] | 21.9 ns [21.8 ns, 21.9 ns] | -0.5% | faster point estimate; marginal CIs separated | 1.00x | 54.7 ns [54.6 ns, 54.9 ns] | 223.9 ns [222.6 ns, 225.8 ns] |
| D=4 | la_stack_norm2_sq | 0.5 ns [0.5 ns, 0.5 ns] | 0.5 ns [0.5 ns, 0.5 ns] | +0.5% | marginal CIs overlap | 1.00x | 0.5 ns [0.5 ns, 0.5 ns] | 4.2 ns [4.2 ns, 4.2 ns] |
| D=4 | la_stack_solve_from_ldlt | 2.6 ns [2.5 ns, 2.6 ns] | 2.5 ns [2.5 ns, 2.5 ns] | -2.5% | faster point estimate; marginal CIs separated | 1.03x | 5.9 ns [5.9 ns, 5.9 ns] | 41.5 ns [41.0 ns, 41.9 ns] |
| D=4 | la_stack_solve_from_lu | 4.0 ns [4.0 ns, 4.1 ns] | 4.0 ns [4.0 ns, 4.0 ns] | -1.9% | faster point estimate; marginal CIs separated | 1.02x | 5.9 ns [5.9 ns, 5.9 ns] | 59.5 ns [57.6 ns, 61.1 ns] |
| D=5 | la_stack_det | 41.0 ns [40.8 ns, 41.2 ns] | 62.2 ns [41.0 ns, 62.6 ns] | +51.6% | marginal CIs overlap | 0.66x | — | — |
| D=5 | la_stack_det_from_ldlt | 1.3 ns [1.3 ns, 1.3 ns] | 1.0 ns [1.0 ns, 1.0 ns] | -24.3% | faster point estimate; marginal CIs separated | 1.32x | 0.6 ns [0.6 ns, 0.6 ns] | 1.3 ns [1.3 ns, 1.3 ns] |
| D=5 | la_stack_det_from_lu | 2.0 ns [1.9 ns, 2.0 ns] | 1.1 ns [1.1 ns, 1.1 ns] | -42.9% | faster point estimate; marginal CIs separated | 1.75x | 0.7 ns [0.7 ns, 0.7 ns] | 1.5 ns [1.5 ns, 1.5 ns] |
| D=5 | la_stack_det_via_lu | 33.3 ns [33.2 ns, 33.5 ns] | 32.6 ns [32.6 ns, 32.8 ns] | -2.0% | faster point estimate; marginal CIs separated | 1.02x | 57.6 ns [57.5 ns, 57.8 ns] | 206.8 ns [206.4 ns, 207.7 ns] |
| D=5 | la_stack_dot | 0.9 ns [0.8 ns, 0.9 ns] | 0.8 ns [0.8 ns, 0.8 ns] | -0.9% | faster point estimate; marginal CIs separated | 1.01x | 0.7 ns [0.7 ns, 0.8 ns] | 1.5 ns [1.5 ns, 1.5 ns] |
| D=5 | la_stack_inf_norm | 3.5 ns [3.5 ns, 3.5 ns] | 3.5 ns [3.5 ns, 3.5 ns] | +1.1% | slower point estimate; marginal CIs separated | 0.99x | 3.3 ns [3.3 ns, 3.3 ns] | 3.2 ns [3.2 ns, 3.2 ns] |
| D=5 | la_stack_ldlt | 43.8 ns [43.7 ns, 43.9 ns] | 43.5 ns [43.4 ns, 43.5 ns] | -0.8% | faster point estimate; marginal CIs separated | 1.01x | 27.7 ns [26.9 ns, 28.7 ns] | 150.0 ns [149.1 ns, 150.8 ns] |
| D=5 | la_stack_ldlt_solve | 47.2 ns [47.1 ns, 47.3 ns] | 47.1 ns [47.0 ns, 47.1 ns] | -0.2% | marginal CIs overlap | 1.00x | 42.3 ns [42.1 ns, 42.4 ns] | 210.3 ns [209.8 ns, 210.9 ns] |
| D=5 | la_stack_lu | 32.8 ns [32.7 ns, 32.9 ns] | 33.0 ns [32.6 ns, 33.2 ns] | +0.6% | marginal CIs overlap | 0.99x | 57.7 ns [57.5 ns, 57.9 ns] | 199.0 ns [198.0 ns, 200.1 ns] |
| D=5 | la_stack_lu_solve | 44.9 ns [44.4 ns, 45.9 ns] | 44.5 ns [44.3 ns, 44.7 ns] | -0.8% | marginal CIs overlap | 1.01x | 71.2 ns [71.1 ns, 71.3 ns] | 293.4 ns [290.1 ns, 306.4 ns] |
| D=5 | la_stack_norm2_sq | 0.5 ns [0.5 ns, 0.5 ns] | 0.5 ns [0.5 ns, 0.5 ns] | -1.1% | faster point estimate; marginal CIs separated | 1.01x | 0.6 ns [0.6 ns, 0.6 ns] | 4.2 ns [4.2 ns, 4.2 ns] |
| D=5 | la_stack_solve_from_ldlt | 3.9 ns [3.9 ns, 3.9 ns] | 4.0 ns [3.9 ns, 4.0 ns] | +1.2% | slower point estimate; marginal CIs separated | 0.99x | 9.7 ns [9.7 ns, 9.8 ns] | 66.3 ns [63.8 ns, 67.5 ns] |
| D=5 | la_stack_solve_from_lu | 6.1 ns [6.1 ns, 6.1 ns] | 6.1 ns [6.1 ns, 6.1 ns] | -0.0% | marginal CIs overlap | 1.00x | 9.4 ns [9.4 ns, 9.4 ns] | 89.2 ns [88.9 ns, 89.5 ns] |
| D=64 | la_stack_det | 14.99 µs [14.95 µs, 15.03 µs] | 14.95 µs [14.88 µs, 14.98 µs] | -0.3% | marginal CIs overlap | 1.00x | — | — |
| D=64 | la_stack_det_from_ldlt | 32.5 ns [32.4 ns, 32.5 ns] | 23.4 ns [23.4 ns, 23.4 ns] | -27.9% | faster point estimate; marginal CIs separated | 1.39x | 8.7 ns [8.7 ns, 8.7 ns] | 21.6 ns [21.5 ns, 21.6 ns] |
| D=64 | la_stack_det_from_lu | 33.5 ns [33.1 ns, 33.5 ns] | 23.6 ns [23.6 ns, 23.7 ns] | -29.3% | faster point estimate; marginal CIs separated | 1.41x | 8.9 ns [8.9 ns, 9.0 ns] | 22.4 ns [22.4 ns, 22.5 ns] |
| D=64 | la_stack_det_via_lu | 15.03 µs [15.01 µs, 15.06 µs] | 14.86 µs [14.84 µs, 14.88 µs] | -1.1% | faster point estimate; marginal CIs separated | 1.01x | 13.17 µs [13.14 µs, 13.19 µs] | 10.70 µs [10.68 µs, 10.71 µs] |
| D=64 | la_stack_dot | 11.2 ns [11.2 ns, 11.3 ns] | 11.0 ns [11.0 ns, 11.0 ns] | -2.1% | faster point estimate; marginal CIs separated | 1.02x | 9.1 ns [9.1 ns, 9.1 ns] | 38.2 ns [38.1 ns, 38.3 ns] |
| D=64 | la_stack_inf_norm | 643.8 ns [642.5 ns, 645.5 ns] | 630.3 ns [628.7 ns, 631.3 ns] | -2.1% | faster point estimate; marginal CIs separated | 1.02x | 1.31 µs [1.30 µs, 1.31 µs] | 1.59 µs [1.58 µs, 1.59 µs] |
| D=64 | la_stack_ldlt | 21.88 µs [21.82 µs, 22.01 µs] | 19.83 µs [19.80 µs, 19.85 µs] | -9.3% | faster point estimate; marginal CIs separated | 1.10x | 10.93 µs [10.92 µs, 10.94 µs] | 9.03 µs [9.02 µs, 9.04 µs] |
| D=64 | la_stack_ldlt_solve | 24.52 µs [24.39 µs, 24.60 µs] | 22.45 µs [22.22 µs, 22.49 µs] | -8.5% | faster point estimate; marginal CIs separated | 1.09x | 12.91 µs [12.89 µs, 12.94 µs] | 10.37 µs [10.35 µs, 10.40 µs] |
| D=64 | la_stack_lu | 14.45 µs [14.42 µs, 14.52 µs] | 14.75 µs [14.74 µs, 14.76 µs] | +2.1% | slower point estimate; marginal CIs separated | 0.98x | 12.98 µs [12.94 µs, 13.01 µs] | 10.56 µs [10.55 µs, 10.58 µs] |
| D=64 | la_stack_lu_solve | 17.51 µs [17.47 µs, 17.55 µs] | 17.36 µs [17.26 µs, 17.41 µs] | -0.9% | faster point estimate; marginal CIs separated | 1.01x | 13.88 µs [13.85 µs, 13.90 µs] | 12.20 µs [12.18 µs, 12.22 µs] |
| D=64 | la_stack_norm2_sq | 11.2 ns [11.2 ns, 11.3 ns] | 11.0 ns [10.9 ns, 11.0 ns] | -2.2% | faster point estimate; marginal CIs separated | 1.02x | 7.4 ns [7.4 ns, 7.4 ns] | 6.3 ns [6.3 ns, 6.3 ns] |
| D=64 | la_stack_solve_from_ldlt | 1.11 µs [1.11 µs, 1.12 µs] | 1.09 µs [1.09 µs, 1.10 µs] | -1.8% | faster point estimate; marginal CIs separated | 1.02x | 1.32 µs [1.31 µs, 1.32 µs] | 1.24 µs [1.24 µs, 1.24 µs] |
| D=64 | la_stack_solve_from_lu | 2.72 µs [2.72 µs, 2.74 µs] | 2.76 µs [2.70 µs, 2.79 µs] | +1.3% | marginal CIs overlap | 0.99x | 796.6 ns [795.2 ns, 797.9 ns] | 1.65 µs [1.65 µs, 1.66 µs] |
| D=8 | la_stack_det | 94.9 ns [94.2 ns, 95.9 ns] | 94.7 ns [94.4 ns, 95.1 ns] | -0.2% | marginal CIs overlap | 1.00x | — | — |
| D=8 | la_stack_det_from_ldlt | 2.7 ns [2.7 ns, 2.7 ns] | 1.3 ns [1.3 ns, 1.3 ns] | -51.7% | faster point estimate; marginal CIs separated | 2.07x | 0.9 ns [0.9 ns, 0.9 ns] | 2.0 ns [2.0 ns, 2.0 ns] |
| D=8 | la_stack_det_from_ldlt_balanced_range | 9.3 ns [9.3 ns, 9.4 ns] | 8.8 ns [8.8 ns, 8.8 ns] | -6.0% | faster point estimate; marginal CIs separated | 1.06x | — | — |
| D=8 | la_stack_det_from_lu | 2.8 ns [2.8 ns, 2.9 ns] | 1.5 ns [1.5 ns, 1.5 ns] | -47.4% | faster point estimate; marginal CIs separated | 1.90x | 1.0 ns [1.0 ns, 1.0 ns] | 2.3 ns [2.3 ns, 2.3 ns] |
| D=8 | la_stack_det_from_lu_balanced_range | 9.3 ns [9.3 ns, 9.3 ns] | 9.3 ns [9.3 ns, 9.4 ns] | +0.3% | marginal CIs overlap | 1.00x | — | — |
| D=8 | la_stack_det_via_lu | 88.0 ns [87.4 ns, 88.7 ns] | 85.2 ns [84.9 ns, 85.3 ns] | -3.2% | faster point estimate; marginal CIs separated | 1.03x | 143.3 ns [143.0 ns, 143.7 ns] | 295.2 ns [293.3 ns, 297.7 ns] |
| D=8 | la_stack_dot | 1.0 ns [1.0 ns, 1.0 ns] | 0.9 ns [0.9 ns, 0.9 ns] | -0.9% | faster point estimate; marginal CIs separated | 1.01x | 1.1 ns [1.1 ns, 1.1 ns] | 2.3 ns [2.3 ns, 2.3 ns] |
| D=8 | la_stack_inf_norm | 8.5 ns [8.5 ns, 8.5 ns] | 8.5 ns [8.5 ns, 8.5 ns] | -0.3% | faster point estimate; marginal CIs separated | 1.00x | 8.1 ns [8.1 ns, 8.1 ns] | 8.2 ns [8.1 ns, 8.2 ns] |
| D=8 | la_stack_ldlt | 94.1 ns [93.1 ns, 94.7 ns] | 92.3 ns [91.4 ns, 93.2 ns] | -1.8% | marginal CIs overlap | 1.02x | 107.4 ns [106.7 ns, 107.9 ns] | 230.5 ns [226.0 ns, 233.4 ns] |
| D=8 | la_stack_ldlt_ill_conditioned | 93.3 ns [91.8 ns, 94.6 ns] | 91.9 ns [90.9 ns, 92.7 ns] | -1.5% | marginal CIs overlap | 1.02x | — | — |
| D=8 | la_stack_ldlt_solve | 103.4 ns [103.1 ns, 103.8 ns] | 100.8 ns [100.6 ns, 101.1 ns] | -2.4% | faster point estimate; marginal CIs separated | 1.02x | 147.7 ns [141.5 ns, 150.9 ns] | 289.2 ns [288.4 ns, 290.9 ns] |
| D=8 | la_stack_lu | 85.6 ns [85.3 ns, 86.0 ns] | 83.5 ns [83.1 ns, 83.7 ns] | -2.5% | faster point estimate; marginal CIs separated | 1.03x | 143.1 ns [141.4 ns, 145.1 ns] | 278.4 ns [277.5 ns, 280.4 ns] |
| D=8 | la_stack_lu_ill_conditioned | 86.4 ns [85.9 ns, 87.0 ns] | 83.1 ns [82.9 ns, 83.2 ns] | -3.9% | faster point estimate; marginal CIs separated | 1.04x | — | — |
| D=8 | la_stack_lu_pivoting | 94.6 ns [94.4 ns, 94.8 ns] | 90.5 ns [90.4 ns, 90.7 ns] | -4.3% | faster point estimate; marginal CIs separated | 1.05x | — | — |
| D=8 | la_stack_lu_solve | 149.5 ns [149.1 ns, 150.1 ns] | 145.4 ns [138.0 ns, 146.6 ns] | -2.7% | faster point estimate; marginal CIs separated | 1.03x | 188.4 ns [187.4 ns, 189.0 ns] | 381.9 ns [380.5 ns, 383.1 ns] |
| D=8 | la_stack_norm2_sq | 0.7 ns [0.7 ns, 0.7 ns] | 0.7 ns [0.7 ns, 0.7 ns] | +0.1% | marginal CIs overlap | 1.00x | 0.7 ns [0.7 ns, 0.7 ns] | 4.2 ns [4.2 ns, 4.2 ns] |
| D=8 | la_stack_solve_from_ldlt | 8.1 ns [8.1 ns, 8.2 ns] | 8.1 ns [8.1 ns, 8.2 ns] | +0.3% | marginal CIs overlap | 1.00x | 23.5 ns [23.4 ns, 23.7 ns] | 71.9 ns [71.4 ns, 72.2 ns] |
| D=8 | la_stack_solve_from_lu | 13.8 ns [13.7 ns, 13.8 ns] | 13.7 ns [13.6 ns, 13.7 ns] | -0.9% | faster point estimate; marginal CIs separated | 1.01x | 16.8 ns [16.8 ns, 16.8 ns] | 96.9 ns [96.7 ns, 97.4 ns] |

## How to Update

Local performance reports are generated in isolated temporary worktrees:

```bash
# Local development: compare the current tree with the latest release
just performance-local

# Release PR: update docs/PERFORMANCE.md and archive the previous report
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
`just performance-release` performs the same measurement and retention work, then promotes distinct-release documentation.
`just performance-doc` consumes the retained pair from either workflow without benchmarking and promotes it when the package versions differ.
For a distinct pair, `performance-local` followed by `performance-doc` is equivalent to the atomic `performance-release` workflow.

Older curated release-to-release reports are archived in `docs/archive/performance/`.

See `docs/BENCHMARKING.md` for the full comparison workflow.
