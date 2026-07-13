# Benchmark Performance

**la-stack** v0.4.4 · `e736c5f` (HEAD)
**Source revision timestamp**: 2026-07-12 17:11:09 UTC (deterministic report metadata; not the benchmark measurement time)
**Benchmark measurement timestamp**: not recorded by Criterion; use the provenance below to identify the measured revisions and environment.
**Statistic**: median
**Suite**: all
**Scope**: release-signal

## Benchmark Results

Comparison against baseline **v0.4.3**:

Negative point-estimate change means the current point estimate is smaller; a baseline/current point-estimate ratio above 1.00 has the same meaning.
The CI-relation column reports only whether the two marginal Criterion intervals overlap. These are not paired confidence intervals
for the change, so the report makes no statistical-significance or performance-improvement claim from interval separation.

### Reproducibility Provenance

**Measurement environment**: recorded for both samples under one shared current harness.

- CPU: `arm`
- OS: `Darwin 25.5.0 arm64`
- rustc: `rustc 1.97.0 (2d8144b78 2026-07-07)`
- Current commit: `e736c5fda155ef23c8712f89ae15bf5369ff3787`
- Current Git clean: `false`
- Current source-state SHA-256: `fb4e2675e4c50d859a2fe358f9dfd8fd7a4e3a31fcb9f6094bf1c43492138be2`
- Baseline commit: `dd4ad192a42e28d9cc72b336b40802fe65cea4f0`
- Baseline Git clean: `false`
- Baseline source-state SHA-256: `1d77b9640c3fe8dce45b486052177fe4d9a44b840605661da0717b8cb0cec9b8`
- Cargo.lock SHA-256: `0c275998d6fe18f8b4def36611598860e96c250303ba459da280ed64e2afd3cd`
- Benchmark harness SHA-256: `c2cc116cf2d77d415f12f30a60bcf3e1b9e0d553c884bc3c57ce6c953124bc08`

**Publication and validation environment**:

- Publication CPU: `arm`
- Publication OS: `Darwin 25.5.0 arm64`
- Publication rustc: `rustc 1.97.0 (2d8144b78 2026-07-07)`
- Publication commit: `e736c5fda155ef23c8712f89ae15bf5369ff3787`
- Publication Git clean: `false`
- Publication source-state SHA-256: `fb4e2675e4c50d859a2fe358f9dfd8fd7a4e3a31fcb9f6094bf1c43492138be2`
- Publication Cargo.lock SHA-256: `0c275998d6fe18f8b4def36611598860e96c250303ba459da280ed64e2afd3cd`
- Publication harness SHA-256: `c2cc116cf2d77d415f12f30a60bcf3e1b9e0d553c884bc3c57ce6c953124bc08`
- Criterion suite/scope: `all` / `release-signal`
- Criterion statistic/sample: `median` / `new`
- Criterion dependency version: `0.8.2`
- Baseline command: `just bench-save-baseline v0.4.3`
- Current command: `just bench-latest`
- Correctness gate: `just test-bench-inputs` passed against both the current and baseline revisions using the shared current fixture harness.
- Validated current revision: `e736c5fda155ef23c8712f89ae15bf5369ff3787` (Git clean: `false`;
  source-state SHA-256: `fb4e2675e4c50d859a2fe358f9dfd8fd7a4e3a31fcb9f6094bf1c43492138be2`)
- Validated baseline revision: `dd4ad192a42e28d9cc72b336b40802fe65cea4f0` (Git clean: `false`;
  source-state SHA-256: `1d77b9640c3fe8dce45b486052177fe4d9a44b840605661da0717b8cb0cec9b8`)
- Baseline API compatibility: `la_stack_v0_4_3_api` selects only source-compatible benchmark calls;
  rows outside the baseline's correctness domain remain explicitly unavailable.
- Baseline-unavailable rows: `d8/la_stack_det_from_lu_balanced_range` and
  `d8/la_stack_det_from_ldlt_balanced_range` were not timed because v0.4.3 returns zero for a
  fixture whose exact determinant is one; current samples remain required, but no speedup is claimed.
- Baseline-unavailable rows: `exact_d2/det_direct_with_errbound`,
  `exact_d3/det_direct_with_errbound`, and `exact_d4/det_direct_with_errbound` were not timed
  because v0.4.3 predates the paired API; the comparable `det_errbound` baselines remain required.

## Exact arithmetic

### D=2

| Benchmark | v0.4.3 (point + CI) | Latest (point + CI) | Point-estimate change | CI relation | Point-estimate ratio |
|-----------|-------:|-------:|-------:|:-----------|--------:|
| det | 0.5 ns [0.5 ns, 0.5 ns] | 0.4 ns [0.4 ns, 0.4 ns] | -11.0% | faster point estimate; marginal CIs separated | 1.12x |
| det_direct | 0.6 ns [0.6 ns, 0.6 ns] | 0.4 ns [0.4 ns, 0.4 ns] | -28.6% | faster point estimate; marginal CIs separated | 1.40x |
| det_errbound | 0.7 ns [0.7 ns, 0.7 ns] | 1.6 ns [1.6 ns, 1.6 ns] | +128.0% | slower point estimate; marginal CIs separated | 0.44x |
| det_exact | 95.7 ns [95.5 ns, 95.8 ns] | 89.3 ns [89.2 ns, 89.4 ns] | -6.7% | faster point estimate; marginal CIs separated | 1.07x |
| det_exact_f64_result | 77.0 ns [76.9 ns, 77.1 ns] | 65.5 ns [65.4 ns, 65.5 ns] | -15.0% | faster point estimate; marginal CIs separated | 1.18x |
| det_exact_rounded_f64 | 245.4 ns [244.9 ns, 245.9 ns] | 66.5 ns [66.5 ns, 66.6 ns] | -72.9% | faster point estimate; marginal CIs separated | 3.69x |
| det_sign_exact | 2.0 ns [2.0 ns, 2.0 ns] | 2.8 ns [2.8 ns, 2.8 ns] | +41.2% | slower point estimate; marginal CIs separated | 0.71x |
| solve_exact | 6.73 µs [6.72 µs, 6.74 µs] | 7.33 µs [7.32 µs, 7.34 µs] | +8.9% | slower point estimate; marginal CIs separated | 0.92x |
| solve_exact_f64_result | 6.89 µs [6.88 µs, 6.91 µs] | 8.38 µs [8.37 µs, 8.39 µs] | +21.6% | slower point estimate; marginal CIs separated | 0.82x |
| solve_exact_rounded_f64 | 7.08 µs [7.07 µs, 7.10 µs] | 7.66 µs [7.64 µs, 7.68 µs] | +8.1% | slower point estimate; marginal CIs separated | 0.93x |

### D=3

| Benchmark | v0.4.3 (point + CI) | Latest (point + CI) | Point-estimate change | CI relation | Point-estimate ratio |
|-----------|-------:|-------:|-------:|:-----------|--------:|
| det | 1.0 ns [1.0 ns, 1.0 ns] | 0.8 ns [0.8 ns, 0.8 ns] | -17.2% | faster point estimate; marginal CIs separated | 1.21x |
| det_direct | 1.0 ns [1.0 ns, 1.0 ns] | 0.8 ns [0.8 ns, 0.8 ns] | -16.6% | faster point estimate; marginal CIs separated | 1.20x |
| det_errbound | 1.6 ns [1.6 ns, 1.6 ns] | 3.4 ns [3.4 ns, 3.4 ns] | +112.5% | slower point estimate; marginal CIs separated | 0.47x |
| det_exact | 334.2 ns [333.8 ns, 335.3 ns] | 320.7 ns [320.2 ns, 321.0 ns] | -4.1% | faster point estimate; marginal CIs separated | 1.04x |
| det_exact_f64_result | 311.4 ns [310.4 ns, 313.1 ns] | 294.9 ns [294.0 ns, 295.3 ns] | -5.3% | faster point estimate; marginal CIs separated | 1.06x |
| det_exact_rounded_f64 | 505.9 ns [503.2 ns, 507.1 ns] | 294.9 ns [294.2 ns, 295.5 ns] | -41.7% | faster point estimate; marginal CIs separated | 1.72x |
| det_sign_exact | 3.7 ns [3.7 ns, 3.7 ns] | 4.7 ns [4.7 ns, 4.7 ns] | +26.9% | slower point estimate; marginal CIs separated | 0.79x |
| solve_exact | 29.70 µs [29.65 µs, 29.77 µs] | 31.46 µs [31.37 µs, 31.51 µs] | +5.9% | slower point estimate; marginal CIs separated | 0.94x |
| solve_exact_f64_result | 29.88 µs [29.80 µs, 30.01 µs] | 33.01 µs [32.96 µs, 33.11 µs] | +10.5% | slower point estimate; marginal CIs separated | 0.90x |
| solve_exact_rounded_f64 | 30.48 µs [30.41 µs, 30.58 µs] | 31.82 µs [31.76 µs, 31.84 µs] | +4.4% | slower point estimate; marginal CIs separated | 0.96x |

### D=4

| Benchmark | v0.4.3 (point + CI) | Latest (point + CI) | Point-estimate change | CI relation | Point-estimate ratio |
|-----------|-------:|-------:|-------:|:-----------|--------:|
| det | 4.5 ns [4.5 ns, 4.5 ns] | 2.3 ns [2.3 ns, 2.3 ns] | -47.7% | faster point estimate; marginal CIs separated | 1.91x |
| det_direct | 4.3 ns [4.3 ns, 4.3 ns] | 2.3 ns [2.3 ns, 2.3 ns] | -47.0% | faster point estimate; marginal CIs separated | 1.89x |
| det_errbound | 6.9 ns [6.9 ns, 6.9 ns] | 6.7 ns [6.7 ns, 6.7 ns] | -2.4% | faster point estimate; marginal CIs separated | 1.02x |
| det_exact | 1.10 µs [1.10 µs, 1.10 µs] | 994.3 ns [990.6 ns, 998.2 ns] | -9.5% | faster point estimate; marginal CIs separated | 1.11x |
| det_exact_f64_result | 1.09 µs [1.09 µs, 1.10 µs] | 959.7 ns [955.9 ns, 961.8 ns] | -12.1% | faster point estimate; marginal CIs separated | 1.14x |
| det_exact_rounded_f64 | 1.28 µs [1.27 µs, 1.29 µs] | 963.6 ns [962.2 ns, 964.4 ns] | -24.6% | faster point estimate; marginal CIs separated | 1.33x |
| det_sign_exact | 11.3 ns [11.3 ns, 11.3 ns] | 7.7 ns [7.7 ns, 7.7 ns] | -31.5% | faster point estimate; marginal CIs separated | 1.46x |
| solve_exact | 77.25 µs [77.10 µs, 77.40 µs] | 79.77 µs [79.57 µs, 79.95 µs] | +3.3% | slower point estimate; marginal CIs separated | 0.97x |
| solve_exact_f64_result | 76.96 µs [76.80 µs, 77.13 µs] | 83.27 µs [83.05 µs, 83.39 µs] | +8.2% | slower point estimate; marginal CIs separated | 0.92x |
| solve_exact_rounded_f64 | 77.73 µs [77.58 µs, 77.82 µs] | 81.50 µs [81.39 µs, 81.62 µs] | +4.8% | slower point estimate; marginal CIs separated | 0.95x |

### D=5

| Benchmark | v0.4.3 (point + CI) | Latest (point + CI) | Point-estimate change | CI relation | Point-estimate ratio |
|-----------|-------:|-------:|-------:|:-----------|--------:|
| det | 22.0 ns [22.0 ns, 22.1 ns] | 25.5 ns [25.4 ns, 25.6 ns] | +16.0% | slower point estimate; marginal CIs separated | 0.86x |
| det_exact | 3.08 µs [3.08 µs, 3.09 µs] | 2.93 µs [2.93 µs, 2.94 µs] | -5.0% | faster point estimate; marginal CIs separated | 1.05x |
| det_exact_f64_result | 3.07 µs [3.07 µs, 3.07 µs] | 2.88 µs [2.88 µs, 2.89 µs] | -6.2% | faster point estimate; marginal CIs separated | 1.07x |
| det_exact_rounded_f64 | 3.30 µs [3.29 µs, 3.30 µs] | 2.90 µs [2.89 µs, 2.90 µs] | -12.1% | faster point estimate; marginal CIs separated | 1.14x |
| det_sign_exact | 3.06 µs [3.05 µs, 3.07 µs] | 2.99 µs [2.98 µs, 3.00 µs] | -2.4% | faster point estimate; marginal CIs separated | 1.03x |
| solve_exact | 151.99 µs [151.61 µs, 152.39 µs] | 159.24 µs [158.84 µs, 159.40 µs] | +4.8% | slower point estimate; marginal CIs separated | 0.95x |
| solve_exact_f64_result | 151.95 µs [151.73 µs, 152.26 µs] | 162.98 µs [162.62 µs, 163.59 µs] | +7.3% | slower point estimate; marginal CIs separated | 0.93x |
| solve_exact_rounded_f64 | 153.08 µs [152.82 µs, 153.19 µs] | 159.84 µs [159.58 µs, 160.22 µs] | +4.4% | slower point estimate; marginal CIs separated | 0.96x |

### Random corpus D=2

| Benchmark | v0.4.3 (point + CI) | Latest (point + CI) | Point-estimate change | CI relation | Point-estimate ratio |
|-----------|-------:|-------:|-------:|:-----------|--------:|
| det_sign_exact | 99.8 ns [99.2 ns, 100.4 ns] | 136.1 ns [136.0 ns, 136.7 ns] | +36.4% | slower point estimate; marginal CIs separated | 0.73x |
| det_exact | 3.22 µs [3.22 µs, 3.23 µs] | 3.06 µs [3.05 µs, 3.06 µs] | -5.2% | faster point estimate; marginal CIs separated | 1.05x |
| solve_exact | 65.74 µs [65.66 µs, 65.83 µs] | 67.56 µs [67.42 µs, 67.79 µs] | +2.8% | slower point estimate; marginal CIs separated | 0.97x |
| solve_exact_f64_result | 66.67 µs [66.58 µs, 66.75 µs] | 77.80 µs [77.55 µs, 77.91 µs] | +16.7% | slower point estimate; marginal CIs separated | 0.86x |
| solve_exact_rounded_f64 | 66.67 µs [66.58 µs, 66.69 µs] | 69.23 µs [69.09 µs, 69.35 µs] | +3.8% | slower point estimate; marginal CIs separated | 0.96x |

### Random corpus D=3

| Benchmark | v0.4.3 (point + CI) | Latest (point + CI) | Point-estimate change | CI relation | Point-estimate ratio |
|-----------|-------:|-------:|-------:|:-----------|--------:|
| det_sign_exact | 174.3 ns [174.2 ns, 174.5 ns] | 225.1 ns [225.0 ns, 225.3 ns] | +29.1% | slower point estimate; marginal CIs separated | 0.77x |
| det_exact | 8.43 µs [8.39 µs, 8.47 µs] | 7.76 µs [7.72 µs, 7.81 µs] | -8.0% | faster point estimate; marginal CIs separated | 1.09x |
| solve_exact | 212.03 µs [211.85 µs, 212.21 µs] | 216.24 µs [216.02 µs, 216.47 µs] | +2.0% | slower point estimate; marginal CIs separated | 0.98x |
| solve_exact_f64_result | 213.65 µs [213.11 µs, 214.14 µs] | 229.74 µs [229.51 µs, 230.18 µs] | +7.5% | slower point estimate; marginal CIs separated | 0.93x |
| solve_exact_rounded_f64 | 213.73 µs [213.55 µs, 214.10 µs] | 218.02 µs [217.72 µs, 218.45 µs] | +2.0% | slower point estimate; marginal CIs separated | 0.98x |

### Random corpus D=4

| Benchmark | v0.4.3 (point + CI) | Latest (point + CI) | Point-estimate change | CI relation | Point-estimate ratio |
|-----------|-------:|-------:|-------:|:-----------|--------:|
| det_sign_exact | 534.8 ns [534.5 ns, 535.0 ns] | 426.5 ns [425.8 ns, 427.4 ns] | -20.3% | faster point estimate; marginal CIs separated | 1.25x |
| det_exact | 28.60 µs [28.34 µs, 28.91 µs] | 24.12 µs [24.09 µs, 24.16 µs] | -15.7% | faster point estimate; marginal CIs separated | 1.19x |
| solve_exact | 486.14 µs [485.47 µs, 487.12 µs] | 495.70 µs [494.68 µs, 496.46 µs] | +2.0% | slower point estimate; marginal CIs separated | 0.98x |
| solve_exact_f64_result | 487.15 µs [486.57 µs, 487.66 µs] | 512.02 µs [511.01 µs, 513.36 µs] | +5.1% | slower point estimate; marginal CIs separated | 0.95x |
| solve_exact_rounded_f64 | 489.02 µs [487.59 µs, 489.93 µs] | 497.20 µs [495.79 µs, 498.14 µs] | +1.7% | slower point estimate; marginal CIs separated | 0.98x |

### Random corpus D=5

| Benchmark | v0.4.3 (point + CI) | Latest (point + CI) | Point-estimate change | CI relation | Point-estimate ratio |
|-----------|-------:|-------:|-------:|:-----------|--------:|
| det_sign_exact | 51.48 µs [51.17 µs, 51.83 µs] | 48.07 µs [47.59 µs, 48.89 µs] | -6.6% | faster point estimate; marginal CIs separated | 1.07x |
| det_exact | 52.09 µs [51.90 µs, 52.27 µs] | 47.94 µs [47.67 µs, 48.20 µs] | -8.0% | faster point estimate; marginal CIs separated | 1.09x |
| solve_exact | 963.73 µs [961.57 µs, 964.57 µs] | 978.88 µs [977.33 µs, 979.71 µs] | +1.6% | slower point estimate; marginal CIs separated | 0.98x |
| solve_exact_f64_result | 960.79 µs [958.85 µs, 962.50 µs] | 995.01 µs [993.26 µs, 997.47 µs] | +3.6% | slower point estimate; marginal CIs separated | 0.97x |
| solve_exact_rounded_f64 | 964.92 µs [963.08 µs, 966.28 µs] | 982.14 µs [980.67 µs, 983.66 µs] | +1.8% | slower point estimate; marginal CIs separated | 0.98x |

### Near-singular 3x3

| Benchmark | v0.4.3 (point + CI) | Latest (point + CI) | Point-estimate change | CI relation | Point-estimate ratio |
|-----------|-------:|-------:|-------:|:-----------|--------:|
| det_sign_exact | 350.2 ns [348.8 ns, 352.5 ns] | 349.4 ns [348.4 ns, 350.5 ns] | -0.2% | marginal CIs overlap | 1.00x |
| det_exact | 306.5 ns [305.5 ns, 307.4 ns] | 304.5 ns [303.9 ns, 305.1 ns] | -0.7% | faster point estimate; marginal CIs separated | 1.01x |
| solve_exact | 2.30 µs [2.29 µs, 2.30 µs] | 2.32 µs [2.31 µs, 2.32 µs] | +0.9% | slower point estimate; marginal CIs separated | 0.99x |
| solve_exact_f64_result | 2.30 µs [2.30 µs, 2.31 µs] | 2.42 µs [2.41 µs, 2.42 µs] | +5.0% | slower point estimate; marginal CIs separated | 0.95x |
| solve_exact_rounded_f64 | 2.32 µs [2.31 µs, 2.32 µs] | 2.32 µs [2.32 µs, 2.33 µs] | +0.3% | marginal CIs overlap | 1.00x |

### Large entries 3x3

| Benchmark | v0.4.3 (point + CI) | Latest (point + CI) | Point-estimate change | CI relation | Point-estimate ratio |
|-----------|-------:|-------:|-------:|:-----------|--------:|
| det_sign_exact | 283.6 ns [283.2 ns, 286.8 ns] | 281.7 ns [280.8 ns, 282.9 ns] | -0.7% | faster point estimate; marginal CIs separated | 1.01x |
| det_exact | 286.3 ns [285.5 ns, 288.1 ns] | 295.7 ns [294.0 ns, 296.9 ns] | +3.3% | slower point estimate; marginal CIs separated | 0.97x |
| solve_exact | 90.64 µs [90.57 µs, 90.79 µs] | 94.07 µs [94.01 µs, 94.20 µs] | +3.8% | slower point estimate; marginal CIs separated | 0.96x |
| solve_exact_f64_result | 91.18 µs [91.03 µs, 91.25 µs] | 95.07 µs [94.88 µs, 95.21 µs] | +4.3% | slower point estimate; marginal CIs separated | 0.96x |
| solve_exact_rounded_f64 | 91.51 µs [91.41 µs, 91.63 µs] | 94.47 µs [94.34 µs, 94.60 µs] | +3.2% | slower point estimate; marginal CIs separated | 0.97x |

### Hilbert 4x4

| Benchmark | v0.4.3 (point + CI) | Latest (point + CI) | Point-estimate change | CI relation | Point-estimate ratio |
|-----------|-------:|-------:|-------:|:-----------|--------:|
| det_sign_exact | 11.3 ns [11.3 ns, 11.3 ns] | 7.7 ns [7.7 ns, 7.7 ns] | -31.5% | faster point estimate; marginal CIs separated | 1.46x |
| det_exact | 1.24 µs [1.24 µs, 1.25 µs] | 1.09 µs [1.09 µs, 1.09 µs] | -12.3% | faster point estimate; marginal CIs separated | 1.14x |
| solve_exact | 57.23 µs [57.13 µs, 57.41 µs] | 59.81 µs [59.71 µs, 59.94 µs] | +4.5% | slower point estimate; marginal CIs separated | 0.96x |
| solve_exact_f64_result | 57.31 µs [57.18 µs, 57.62 µs] | 61.72 µs [61.51 µs, 61.85 µs] | +7.7% | slower point estimate; marginal CIs separated | 0.93x |
| solve_exact_rounded_f64 | 58.02 µs [57.92 µs, 58.14 µs] | 60.30 µs [60.18 µs, 60.43 µs] | +3.9% | slower point estimate; marginal CIs separated | 0.96x |

### Hilbert 5x5

| Benchmark | v0.4.3 (point + CI) | Latest (point + CI) | Point-estimate change | CI relation | Point-estimate ratio |
|-----------|-------:|-------:|-------:|:-----------|--------:|
| det_sign_exact | 3.23 µs [3.22 µs, 3.24 µs] | 3.23 µs [3.23 µs, 3.24 µs] | +0.1% | marginal CIs overlap | 1.00x |
| det_exact | 3.26 µs [3.25 µs, 3.27 µs] | 3.00 µs [2.99 µs, 3.01 µs] | -8.0% | faster point estimate; marginal CIs separated | 1.09x |
| solve_exact | 119.87 µs [119.55 µs, 120.21 µs] | 124.37 µs [123.95 µs, 124.57 µs] | +3.8% | slower point estimate; marginal CIs separated | 0.96x |
| solve_exact_f64_result | 119.44 µs [119.30 µs, 119.59 µs] | 127.21 µs [126.96 µs, 127.34 µs] | +6.5% | slower point estimate; marginal CIs separated | 0.94x |
| solve_exact_rounded_f64 | 120.33 µs [120.13 µs, 120.48 µs] | 125.50 µs [125.24 µs, 125.78 µs] | +4.3% | slower point estimate; marginal CIs separated | 0.96x |

## vs_linalg

### D=2

| Benchmark | v0.4.3 (point + CI) | Latest (point + CI) | Point-estimate change | CI relation | Point-estimate ratio | v0.4.3 nalgebra | v0.4.3 faer |
|-----------|-------:|-------:|-------:|:-----------|--------:|-------:|-------:|
| la_stack_det_via_lu | 1.6 ns [1.6 ns, 1.6 ns] | 2.7 ns [2.7 ns, 2.7 ns] | +67.9% | slower point estimate; marginal CIs separated | 0.60x | 0.8 ns [0.8 ns, 0.8 ns] | 99.0 ns [97.7 ns, 99.4 ns] |
| la_stack_det | 0.7 ns [0.7 ns, 0.7 ns] | 0.6 ns [0.6 ns, 0.6 ns] | -16.0% | faster point estimate; marginal CIs separated | 1.19x | — | — |
| la_stack_lu | 1.3 ns [1.3 ns, 1.4 ns] | 1.9 ns [1.9 ns, 1.9 ns] | +40.5% | slower point estimate; marginal CIs separated | 0.71x | 1.6 ns [1.6 ns, 1.6 ns] | 87.9 ns [87.6 ns, 88.2 ns] |
| la_stack_ldlt | 2.1 ns [2.1 ns, 2.1 ns] | 6.7 ns [6.6 ns, 6.8 ns] | +225.3% | slower point estimate; marginal CIs separated | 0.31x | 1.7 ns [1.7 ns, 1.7 ns] | 79.3 ns [78.8 ns, 79.7 ns] |
| la_stack_lu_solve | 2.0 ns [2.0 ns, 2.0 ns] | 2.0 ns [2.0 ns, 2.0 ns] | -0.6% | faster point estimate; marginal CIs separated | 1.01x | 4.5 ns [4.5 ns, 4.5 ns] | 147.8 ns [147.2 ns, 148.8 ns] |
| la_stack_ldlt_solve | 4.0 ns [4.0 ns, 4.0 ns] | 10.0 ns [9.9 ns, 10.0 ns] | +146.2% | slower point estimate; marginal CIs separated | 0.41x | 2.8 ns [2.8 ns, 2.8 ns] | 122.8 ns [122.2 ns, 123.1 ns] |
| la_stack_solve_from_lu | 1.2 ns [1.2 ns, 1.2 ns] | 1.3 ns [1.2 ns, 1.3 ns] | +0.7% | slower point estimate; marginal CIs separated | 0.99x | 2.7 ns [2.7 ns, 2.7 ns] | 47.2 ns [47.1 ns, 47.4 ns] |
| la_stack_solve_from_ldlt | 1.2 ns [1.2 ns, 1.2 ns] | 1.2 ns [1.2 ns, 1.2 ns] | +3.4% | slower point estimate; marginal CIs separated | 0.97x | 1.3 ns [1.3 ns, 1.3 ns] | 38.6 ns [38.0 ns, 39.0 ns] |
| la_stack_det_from_lu | 0.5 ns [0.5 ns, 0.5 ns] | 0.5 ns [0.5 ns, 0.5 ns] | +4.7% | slower point estimate; marginal CIs separated | 0.95x | 0.5 ns [0.5 ns, 0.5 ns] | 0.7 ns [0.7 ns, 0.7 ns] |
| la_stack_det_from_ldlt | 0.5 ns [0.5 ns, 0.5 ns] | 0.5 ns [0.5 ns, 0.5 ns] | +0.4% | slower point estimate; marginal CIs separated | 1.00x | 0.4 ns [0.4 ns, 0.4 ns] | 0.6 ns [0.6 ns, 0.6 ns] |
| la_stack_dot | 0.7 ns [0.7 ns, 0.7 ns] | 0.6 ns [0.6 ns, 0.6 ns] | -14.0% | faster point estimate; marginal CIs separated | 1.16x | 0.6 ns [0.6 ns, 0.6 ns] | 0.7 ns [0.7 ns, 0.7 ns] |
| la_stack_norm2_sq | 0.5 ns [0.5 ns, 0.5 ns] | 0.4 ns [0.4 ns, 0.4 ns] | -20.2% | faster point estimate; marginal CIs separated | 1.25x | 0.4 ns [0.4 ns, 0.4 ns] | 4.2 ns [4.2 ns, 4.2 ns] |
| la_stack_inf_norm | 0.8 ns [0.8 ns, 0.8 ns] | 0.6 ns [0.6 ns, 0.6 ns] | -24.2% | faster point estimate; marginal CIs separated | 1.32x | 0.5 ns [0.5 ns, 0.5 ns] | 0.8 ns [0.8 ns, 0.8 ns] |

### D=3

| Benchmark | v0.4.3 (point + CI) | Latest (point + CI) | Point-estimate change | CI relation | Point-estimate ratio | v0.4.3 nalgebra | v0.4.3 faer |
|-----------|-------:|-------:|-------:|:-----------|--------:|-------:|-------:|
| la_stack_det_via_lu | 7.5 ns [7.5 ns, 7.5 ns] | 10.6 ns [10.6 ns, 10.6 ns] | +41.0% | slower point estimate; marginal CIs separated | 0.71x | 16.7 ns [16.6 ns, 16.8 ns] | 137.9 ns [137.1 ns, 138.9 ns] |
| la_stack_det | 1.5 ns [1.5 ns, 1.5 ns] | 2.2 ns [2.2 ns, 2.2 ns] | +46.4% | slower point estimate; marginal CIs separated | 0.68x | — | — |
| la_stack_lu | 9.1 ns [9.1 ns, 9.1 ns] | 8.5 ns [8.4 ns, 8.5 ns] | -7.0% | faster point estimate; marginal CIs separated | 1.08x | 15.2 ns [15.1 ns, 15.4 ns] | 127.3 ns [126.3 ns, 128.1 ns] |
| la_stack_ldlt | 7.1 ns [7.1 ns, 7.1 ns] | 14.0 ns [14.0 ns, 14.1 ns] | +97.2% | slower point estimate; marginal CIs separated | 0.51x | 4.1 ns [4.1 ns, 4.1 ns] | 93.6 ns [93.2 ns, 94.0 ns] |
| la_stack_lu_solve | 9.6 ns [9.6 ns, 9.7 ns] | 9.9 ns [9.9 ns, 9.9 ns] | +2.8% | slower point estimate; marginal CIs separated | 0.97x | 23.0 ns [22.9 ns, 23.2 ns] | 191.2 ns [190.5 ns, 191.8 ns] |
| la_stack_ldlt_solve | 8.6 ns [8.6 ns, 8.7 ns] | 15.4 ns [15.3 ns, 15.4 ns] | +78.5% | slower point estimate; marginal CIs separated | 0.56x | 8.2 ns [8.2 ns, 8.2 ns] | 137.6 ns [137.3 ns, 138.2 ns] |
| la_stack_solve_from_lu | 2.1 ns [2.1 ns, 2.1 ns] | 2.1 ns [2.1 ns, 2.1 ns] | -0.3% | faster point estimate; marginal CIs separated | 1.00x | 4.4 ns [4.4 ns, 4.4 ns] | 49.5 ns [49.3 ns, 49.7 ns] |
| la_stack_solve_from_ldlt | 1.8 ns [1.8 ns, 1.8 ns] | 1.8 ns [1.8 ns, 1.8 ns] | +0.0% | marginal CIs overlap | 1.00x | 2.9 ns [2.9 ns, 2.9 ns] | 38.2 ns [37.8 ns, 38.5 ns] |
| la_stack_det_from_lu | 0.6 ns [0.6 ns, 0.6 ns] | 0.7 ns [0.7 ns, 0.7 ns] | +28.1% | slower point estimate; marginal CIs separated | 0.78x | 0.5 ns [0.5 ns, 0.5 ns] | 1.0 ns [1.0 ns, 1.0 ns] |
| la_stack_det_from_ldlt | 0.5 ns [0.5 ns, 0.5 ns] | 0.7 ns [0.7 ns, 0.7 ns] | +40.3% | slower point estimate; marginal CIs separated | 0.71x | 0.5 ns [0.5 ns, 0.5 ns] | 0.7 ns [0.7 ns, 0.7 ns] |
| la_stack_dot | 0.8 ns [0.8 ns, 0.8 ns] | 0.7 ns [0.7 ns, 0.7 ns] | -12.8% | faster point estimate; marginal CIs separated | 1.15x | 0.7 ns [0.7 ns, 0.7 ns] | 0.9 ns [0.9 ns, 0.9 ns] |
| la_stack_norm2_sq | 0.5 ns [0.5 ns, 0.5 ns] | 0.4 ns [0.4 ns, 0.4 ns] | -18.4% | faster point estimate; marginal CIs separated | 1.22x | 0.4 ns [0.4 ns, 0.4 ns] | 4.2 ns [4.1 ns, 4.2 ns] |
| la_stack_inf_norm | 1.7 ns [1.7 ns, 1.7 ns] | 1.3 ns [1.3 ns, 1.3 ns] | -24.1% | faster point estimate; marginal CIs separated | 1.32x | 1.1 ns [1.1 ns, 1.1 ns] | 1.3 ns [1.3 ns, 1.3 ns] |

### D=4

| Benchmark | v0.4.3 (point + CI) | Latest (point + CI) | Point-estimate change | CI relation | Point-estimate ratio | v0.4.3 nalgebra | v0.4.3 faer |
|-----------|-------:|-------:|-------:|:-----------|--------:|-------:|-------:|
| la_stack_det_via_lu | 12.6 ns [12.6 ns, 12.7 ns] | 17.0 ns [16.9 ns, 17.0 ns] | +34.3% | slower point estimate; marginal CIs separated | 0.74x | 30.6 ns [30.5 ns, 30.6 ns] | 162.0 ns [161.2 ns, 164.0 ns] |
| la_stack_det | 6.1 ns [6.1 ns, 6.1 ns] | 4.9 ns [4.9 ns, 4.9 ns] | -20.1% | faster point estimate; marginal CIs separated | 1.25x | — | — |
| la_stack_lu | 13.7 ns [13.6 ns, 13.7 ns] | 13.8 ns [13.8 ns, 13.8 ns] | +0.9% | slower point estimate; marginal CIs separated | 0.99x | 29.6 ns [29.5 ns, 29.7 ns] | 148.2 ns [147.6 ns, 148.9 ns] |
| la_stack_ldlt | 12.6 ns [12.6 ns, 12.7 ns] | 21.5 ns [21.4 ns, 21.5 ns] | +69.8% | slower point estimate; marginal CIs separated | 0.59x | 10.0 ns [10.0 ns, 10.1 ns] | 115.8 ns [115.2 ns, 116.2 ns] |
| la_stack_lu_solve | 22.2 ns [22.1 ns, 22.2 ns] | 21.6 ns [21.6 ns, 21.6 ns] | -2.8% | faster point estimate; marginal CIs separated | 1.03x | 51.9 ns [51.8 ns, 52.0 ns] | 217.5 ns [216.4 ns, 218.8 ns] |
| la_stack_ldlt_solve | 17.5 ns [17.5 ns, 17.6 ns] | 26.8 ns [26.8 ns, 26.8 ns] | +52.7% | slower point estimate; marginal CIs separated | 0.65x | 14.9 ns [14.8 ns, 14.9 ns] | 159.4 ns [159.2 ns, 159.8 ns] |
| la_stack_solve_from_lu | 4.0 ns [4.0 ns, 4.0 ns] | 4.0 ns [4.0 ns, 4.0 ns] | -0.2% | faster point estimate; marginal CIs separated | 1.00x | 5.1 ns [5.1 ns, 5.1 ns] | 51.8 ns [51.7 ns, 51.9 ns] |
| la_stack_solve_from_ldlt | 2.5 ns [2.5 ns, 2.5 ns] | 2.5 ns [2.5 ns, 2.5 ns] | +0.1% | marginal CIs overlap | 1.00x | 5.5 ns [5.5 ns, 5.5 ns] | 39.0 ns [38.7 ns, 39.5 ns] |
| la_stack_det_from_lu | 0.7 ns [0.7 ns, 0.7 ns] | 0.9 ns [0.9 ns, 0.9 ns] | +21.1% | slower point estimate; marginal CIs separated | 0.83x | 0.6 ns [0.6 ns, 0.6 ns] | 1.2 ns [1.2 ns, 1.2 ns] |
| la_stack_det_from_ldlt | 0.7 ns [0.7 ns, 0.7 ns] | 1.0 ns [0.9 ns, 1.0 ns] | +47.9% | slower point estimate; marginal CIs separated | 0.68x | 0.5 ns [0.5 ns, 0.5 ns] | 1.0 ns [1.0 ns, 1.0 ns] |
| la_stack_dot | 1.0 ns [1.0 ns, 1.0 ns] | 0.7 ns [0.7 ns, 0.7 ns] | -29.3% | faster point estimate; marginal CIs separated | 1.41x | 0.7 ns [0.7 ns, 0.7 ns] | 1.2 ns [1.2 ns, 1.2 ns] |
| la_stack_norm2_sq | 0.7 ns [0.7 ns, 0.7 ns] | 0.5 ns [0.5 ns, 0.5 ns] | -32.7% | faster point estimate; marginal CIs separated | 1.49x | 0.5 ns [0.5 ns, 0.5 ns] | 4.2 ns [4.1 ns, 4.2 ns] |
| la_stack_inf_norm | 3.0 ns [3.0 ns, 3.0 ns] | 2.2 ns [2.2 ns, 2.2 ns] | -25.4% | faster point estimate; marginal CIs separated | 1.34x | 2.0 ns [2.0 ns, 2.0 ns] | 2.0 ns [2.0 ns, 2.0 ns] |

### D=5

| Benchmark | v0.4.3 (point + CI) | Latest (point + CI) | Point-estimate change | CI relation | Point-estimate ratio | v0.4.3 nalgebra | v0.4.3 faer |
|-----------|-------:|-------:|-------:|:-----------|--------:|-------:|-------:|
| la_stack_det_via_lu | 26.1 ns [26.0 ns, 26.1 ns] | 37.8 ns [37.8 ns, 38.0 ns] | +45.0% | slower point estimate; marginal CIs separated | 0.69x | 55.6 ns [55.4 ns, 55.7 ns] | 203.2 ns [201.2 ns, 206.1 ns] |
| la_stack_det | 26.5 ns [26.4 ns, 26.5 ns] | 39.6 ns [39.6 ns, 39.9 ns] | +49.6% | slower point estimate; marginal CIs separated | 0.67x | — | — |
| la_stack_lu | 26.8 ns [26.7 ns, 26.9 ns] | 32.3 ns [32.2 ns, 32.5 ns] | +20.6% | slower point estimate; marginal CIs separated | 0.83x | 54.7 ns [54.6 ns, 54.9 ns] | 194.4 ns [191.2 ns, 196.0 ns] |
| la_stack_ldlt | 19.8 ns [19.7 ns, 19.9 ns] | 42.2 ns [42.1 ns, 42.2 ns] | +113.0% | slower point estimate; marginal CIs separated | 0.47x | 14.7 ns [14.6 ns, 14.9 ns] | 137.3 ns [136.8 ns, 137.7 ns] |
| la_stack_lu_solve | 45.7 ns [45.6 ns, 46.1 ns] | 45.6 ns [44.9 ns, 46.1 ns] | -0.2% | marginal CIs overlap | 1.00x | 68.9 ns [68.8 ns, 69.0 ns] | 323.2 ns [302.8 ns, 329.6 ns] |
| la_stack_ldlt_solve | 24.7 ns [24.7 ns, 24.8 ns] | 54.9 ns [46.3 ns, 55.5 ns] | +121.7% | slower point estimate; marginal CIs separated | 0.45x | 61.5 ns [61.0 ns, 61.9 ns] | 219.5 ns [215.2 ns, 229.6 ns] |
| la_stack_solve_from_lu | 6.2 ns [6.2 ns, 6.3 ns] | 6.4 ns [6.4 ns, 6.4 ns] | +2.5% | slower point estimate; marginal CIs separated | 0.98x | 8.0 ns [8.0 ns, 8.1 ns] | 89.3 ns [88.2 ns, 89.9 ns] |
| la_stack_solve_from_ldlt | 3.9 ns [3.9 ns, 3.9 ns] | 3.9 ns [3.9 ns, 3.9 ns] | +1.4% | slower point estimate; marginal CIs separated | 0.99x | 9.2 ns [9.1 ns, 9.3 ns] | 66.4 ns [65.7 ns, 66.8 ns] |
| la_stack_det_from_lu | 0.9 ns [0.9 ns, 0.9 ns] | 1.7 ns [1.7 ns, 1.7 ns] | +98.3% | slower point estimate; marginal CIs separated | 0.50x | 0.7 ns [0.7 ns, 0.7 ns] | 1.5 ns [1.5 ns, 1.5 ns] |
| la_stack_det_from_ldlt | 0.8 ns [0.8 ns, 0.8 ns] | 1.2 ns [1.2 ns, 1.2 ns] | +51.4% | slower point estimate; marginal CIs separated | 0.66x | 0.6 ns [0.6 ns, 0.6 ns] | 1.2 ns [1.2 ns, 1.2 ns] |
| la_stack_dot | 1.2 ns [1.2 ns, 1.2 ns] | 0.8 ns [0.8 ns, 0.8 ns] | -32.9% | faster point estimate; marginal CIs separated | 1.49x | 0.7 ns [0.7 ns, 0.7 ns] | 1.4 ns [1.4 ns, 1.4 ns] |
| la_stack_norm2_sq | 0.9 ns [0.9 ns, 0.9 ns] | 0.5 ns [0.5 ns, 0.5 ns] | -38.6% | faster point estimate; marginal CIs separated | 1.63x | 0.6 ns [0.6 ns, 0.6 ns] | 4.3 ns [4.2 ns, 4.3 ns] |
| la_stack_inf_norm | 4.7 ns [4.7 ns, 4.7 ns] | 3.4 ns [3.4 ns, 3.4 ns] | -26.9% | faster point estimate; marginal CIs separated | 1.37x | 3.2 ns [3.2 ns, 3.2 ns] | 3.2 ns [3.2 ns, 3.2 ns] |

### D=8

| Benchmark | v0.4.3 (point + CI) | Latest (point + CI) | Point-estimate change | CI relation | Point-estimate ratio | v0.4.3 nalgebra | v0.4.3 faer |
|-----------|-------:|-------:|-------:|:-----------|--------:|-------:|-------:|
| la_stack_det_via_lu | 83.5 ns [83.1 ns, 83.7 ns] | 88.3 ns [88.1 ns, 88.7 ns] | +5.7% | slower point estimate; marginal CIs separated | 0.95x | 137.9 ns [137.5 ns, 138.5 ns] | 280.8 ns [279.8 ns, 281.8 ns] |
| la_stack_det | 81.9 ns [81.6 ns, 82.1 ns] | 90.3 ns [90.1 ns, 91.1 ns] | +10.3% | slower point estimate; marginal CIs separated | 0.91x | — | — |
| la_stack_lu | 83.7 ns [82.9 ns, 84.6 ns] | 82.5 ns [82.2 ns, 82.6 ns] | -1.4% | faster point estimate; marginal CIs separated | 1.01x | 123.5 ns [123.2 ns, 124.5 ns] | 265.6 ns [263.8 ns, 266.7 ns] |
| la_stack_ldlt | 100.9 ns [100.7 ns, 101.1 ns] | 92.0 ns [90.4 ns, 92.3 ns] | -8.8% | faster point estimate; marginal CIs separated | 1.10x | 97.7 ns [97.4 ns, 97.8 ns] | 207.9 ns [205.9 ns, 209.1 ns] |
| la_stack_lu_solve | 127.2 ns [126.9 ns, 127.5 ns] | 142.9 ns [141.7 ns, 143.7 ns] | +12.3% | slower point estimate; marginal CIs separated | 0.89x | 165.3 ns [164.7 ns, 166.2 ns] | 367.7 ns [367.0 ns, 368.7 ns] |
| la_stack_ldlt_solve | 119.7 ns [119.5 ns, 119.8 ns] | 105.1 ns [104.8 ns, 105.4 ns] | -12.2% | faster point estimate; marginal CIs separated | 1.14x | 148.2 ns [147.6 ns, 148.3 ns] | 272.8 ns [272.0 ns, 273.8 ns] |
| la_stack_solve_from_lu | 13.5 ns [13.5 ns, 13.5 ns] | 13.5 ns [13.5 ns, 13.5 ns] | -0.2% | faster point estimate; marginal CIs separated | 1.00x | 13.4 ns [13.3 ns, 13.4 ns] | 94.9 ns [94.7 ns, 95.6 ns] |
| la_stack_solve_from_ldlt | 8.3 ns [8.3 ns, 8.3 ns] | 8.2 ns [8.2 ns, 8.2 ns] | -0.5% | faster point estimate; marginal CIs separated | 1.00x | 21.9 ns [21.9 ns, 22.0 ns] | 69.5 ns [69.4 ns, 69.7 ns] |
| la_stack_det_from_lu | 1.3 ns [1.3 ns, 1.3 ns] | 2.2 ns [2.2 ns, 2.5 ns] | +69.6% | slower point estimate; marginal CIs separated | 0.59x | 1.0 ns [1.0 ns, 1.0 ns] | 2.2 ns [2.2 ns, 2.2 ns] |
| la_stack_det_from_ldlt | 1.3 ns [1.3 ns, 1.3 ns] | 2.5 ns [2.4 ns, 2.5 ns] | +94.9% | slower point estimate; marginal CIs separated | 0.51x | 0.9 ns [0.9 ns, 0.9 ns] | 2.1 ns [2.0 ns, 2.1 ns] |
| la_stack_dot | 1.5 ns [1.5 ns, 1.6 ns] | 1.0 ns [1.0 ns, 1.0 ns] | -33.3% | faster point estimate; marginal CIs separated | 1.50x | 1.1 ns [1.1 ns, 1.1 ns] | 2.3 ns [2.3 ns, 2.3 ns] |
| la_stack_norm2_sq | 1.3 ns [1.3 ns, 1.3 ns] | 0.7 ns [0.7 ns, 0.7 ns] | -48.8% | faster point estimate; marginal CIs separated | 1.95x | 0.7 ns [0.7 ns, 0.7 ns] | 4.1 ns [4.1 ns, 4.1 ns] |
| la_stack_inf_norm | 12.4 ns [12.4 ns, 12.5 ns] | 8.4 ns [8.4 ns, 8.4 ns] | -32.2% | faster point estimate; marginal CIs separated | 1.48x | 8.0 ns [8.0 ns, 8.1 ns] | 8.1 ns [8.1 ns, 8.1 ns] |
| la_stack_lu_pivoting | 103.7 ns [103.4 ns, 104.0 ns] | 91.3 ns [90.9 ns, 91.6 ns] | -12.0% | faster point estimate; marginal CIs separated | 1.14x | — | — |
| la_stack_lu_ill_conditioned | 96.7 ns [96.5 ns, 96.9 ns] | 82.5 ns [82.4 ns, 82.7 ns] | -14.7% | faster point estimate; marginal CIs separated | 1.17x | — | — |
| la_stack_ldlt_ill_conditioned | 100.4 ns [100.2 ns, 100.5 ns] | 92.0 ns [90.2 ns, 92.5 ns] | -8.4% | faster point estimate; marginal CIs separated | 1.09x | — | — |

### D=16

| Benchmark | v0.4.3 (point + CI) | Latest (point + CI) | Point-estimate change | CI relation | Point-estimate ratio | v0.4.3 nalgebra | v0.4.3 faer |
|-----------|-------:|-------:|-------:|:-----------|--------:|-------:|-------:|
| la_stack_det_via_lu | 398.9 ns [396.4 ns, 402.7 ns] | 427.7 ns [426.4 ns, 428.5 ns] | +7.2% | slower point estimate; marginal CIs separated | 0.93x | 460.8 ns [458.2 ns, 466.5 ns] | 662.5 ns [659.9 ns, 664.7 ns] |
| la_stack_det | 393.5 ns [392.8 ns, 395.5 ns] | 457.5 ns [444.9 ns, 461.4 ns] | +16.3% | slower point estimate; marginal CIs separated | 0.86x | — | — |
| la_stack_lu | 367.8 ns [367.3 ns, 368.7 ns] | 399.6 ns [398.4 ns, 404.0 ns] | +8.6% | slower point estimate; marginal CIs separated | 0.92x | 468.6 ns [468.4 ns, 469.5 ns] | 642.3 ns [637.6 ns, 646.4 ns] |
| la_stack_ldlt | 449.6 ns [448.6 ns, 450.4 ns] | 393.0 ns [392.5 ns, 393.5 ns] | -12.6% | faster point estimate; marginal CIs separated | 1.14x | 405.6 ns [405.2 ns, 405.8 ns] | 416.8 ns [415.5 ns, 418.1 ns] |
| la_stack_lu_solve | 638.4 ns [636.2 ns, 639.9 ns] | 654.7 ns [653.6 ns, 656.0 ns] | +2.6% | slower point estimate; marginal CIs separated | 0.98x | 577.2 ns [576.7 ns, 577.9 ns] | 894.6 ns [891.7 ns, 896.5 ns] |
| la_stack_ldlt_solve | 518.9 ns [517.5 ns, 520.2 ns] | 451.7 ns [451.4 ns, 452.3 ns] | -12.9% | faster point estimate; marginal CIs separated | 1.15x | 649.4 ns [647.5 ns, 651.8 ns] | 590.8 ns [588.9 ns, 593.7 ns] |
| la_stack_solve_from_lu | 196.9 ns [196.4 ns, 197.5 ns] | 197.2 ns [196.3 ns, 198.0 ns] | +0.2% | marginal CIs overlap | 1.00x | 93.4 ns [93.3 ns, 93.5 ns] | 239.7 ns [239.3 ns, 240.2 ns] |
| la_stack_solve_from_ldlt | 28.2 ns [28.1 ns, 28.3 ns] | 27.7 ns [27.7 ns, 27.7 ns] | -1.7% | faster point estimate; marginal CIs separated | 1.02x | 123.6 ns [123.3 ns, 123.8 ns] | 177.5 ns [177.1 ns, 177.8 ns] |
| la_stack_det_from_lu | 2.6 ns [2.6 ns, 2.6 ns] | 3.7 ns [3.7 ns, 3.7 ns] | +42.5% | slower point estimate; marginal CIs separated | 0.70x | 1.8 ns [1.8 ns, 1.8 ns] | 4.8 ns [4.7 ns, 4.9 ns] |
| la_stack_det_from_ldlt | 2.5 ns [2.5 ns, 2.6 ns] | 3.3 ns [3.2 ns, 3.5 ns] | +30.0% | slower point estimate; marginal CIs separated | 0.77x | 1.7 ns [1.7 ns, 1.7 ns] | 4.4 ns [4.4 ns, 4.4 ns] |
| la_stack_dot | 3.1 ns [3.0 ns, 3.3 ns] | 2.3 ns [2.3 ns, 2.3 ns] | -26.3% | faster point estimate; marginal CIs separated | 1.36x | 1.9 ns [1.9 ns, 1.9 ns] | 4.6 ns [4.5 ns, 4.6 ns] |
| la_stack_norm2_sq | 2.7 ns [2.6 ns, 2.7 ns] | 2.1 ns [2.0 ns, 2.1 ns] | -23.1% | faster point estimate; marginal CIs separated | 1.30x | 1.5 ns [1.5 ns, 1.5 ns] | 4.1 ns [4.1 ns, 4.1 ns] |
| la_stack_inf_norm | 49.8 ns [49.5 ns, 50.0 ns] | 32.7 ns [32.6 ns, 32.7 ns] | -34.3% | faster point estimate; marginal CIs separated | 1.52x | 31.7 ns [31.7 ns, 31.7 ns] | 32.6 ns [32.5 ns, 32.6 ns] |

### D=32

| Benchmark | v0.4.3 (point + CI) | Latest (point + CI) | Point-estimate change | CI relation | Point-estimate ratio | v0.4.3 nalgebra | v0.4.3 faer |
|-----------|-------:|-------:|-------:|:-----------|--------:|-------:|-------:|
| la_stack_det_via_lu | 2.11 µs [2.10 µs, 2.11 µs] | 2.21 µs [2.21 µs, 2.22 µs] | +5.0% | slower point estimate; marginal CIs separated | 0.95x | 2.51 µs [2.50 µs, 2.52 µs] | 2.30 µs [2.28 µs, 2.31 µs] |
| la_stack_det | 2.01 µs [2.00 µs, 2.02 µs] | 2.60 µs [2.60 µs, 2.61 µs] | +29.6% | slower point estimate; marginal CIs separated | 0.77x | — | — |
| la_stack_lu | 2.20 µs [2.18 µs, 2.21 µs] | 2.07 µs [2.06 µs, 2.07 µs] | -6.1% | faster point estimate; marginal CIs separated | 1.07x | 2.16 µs [2.16 µs, 2.16 µs] | 2.26 µs [2.25 µs, 2.27 µs] |
| la_stack_ldlt | 2.81 µs [2.81 µs, 2.82 µs] | 2.53 µs [2.52 µs, 2.54 µs] | -9.9% | faster point estimate; marginal CIs separated | 1.11x | 2.11 µs [2.10 µs, 2.11 µs] | 1.42 µs [1.42 µs, 1.42 µs] |
| la_stack_lu_solve | 2.68 µs [2.68 µs, 2.68 µs] | 2.79 µs [2.79 µs, 2.80 µs] | +4.3% | slower point estimate; marginal CIs separated | 0.96x | 2.45 µs [2.45 µs, 2.46 µs] | 2.86 µs [2.86 µs, 2.87 µs] |
| la_stack_ldlt_solve | 3.27 µs [3.26 µs, 3.27 µs] | 2.94 µs [2.92 µs, 2.95 µs] | -10.1% | faster point estimate; marginal CIs separated | 1.11x | 2.79 µs [2.79 µs, 2.79 µs] | 1.91 µs [1.91 µs, 1.91 µs] |
| la_stack_solve_from_lu | 688.9 ns [686.2 ns, 692.0 ns] | 674.0 ns [672.6 ns, 674.9 ns] | -2.2% | faster point estimate; marginal CIs separated | 1.02x | 331.4 ns [331.1 ns, 331.8 ns] | 619.2 ns [617.2 ns, 622.7 ns] |
| la_stack_solve_from_ldlt | 303.6 ns [302.9 ns, 304.0 ns] | 307.3 ns [306.8 ns, 307.7 ns] | +1.2% | slower point estimate; marginal CIs separated | 0.99x | 566.9 ns [564.8 ns, 570.2 ns] | 463.7 ns [463.0 ns, 464.2 ns] |
| la_stack_det_from_lu | 7.0 ns [7.0 ns, 7.0 ns] | 13.9 ns [13.4 ns, 14.0 ns] | +98.2% | slower point estimate; marginal CIs separated | 0.50x | 3.1 ns [3.1 ns, 3.1 ns] | 8.8 ns [8.8 ns, 8.8 ns] |
| la_stack_det_from_ldlt | 6.3 ns [6.3 ns, 6.4 ns] | 12.9 ns [12.9 ns, 13.0 ns] | +105.4% | slower point estimate; marginal CIs separated | 0.49x | 3.0 ns [3.0 ns, 3.0 ns] | 8.5 ns [8.5 ns, 8.5 ns] |
| la_stack_dot | 7.5 ns [7.5 ns, 7.6 ns] | 4.0 ns [4.0 ns, 4.0 ns] | -46.4% | faster point estimate; marginal CIs separated | 1.87x | 4.7 ns [4.7 ns, 4.7 ns] | 15.1 ns [15.1 ns, 15.1 ns] |
| la_stack_norm2_sq | 7.5 ns [7.5 ns, 7.5 ns] | 4.0 ns [4.0 ns, 4.0 ns] | -46.9% | faster point estimate; marginal CIs separated | 1.88x | 3.9 ns [3.8 ns, 3.9 ns] | 4.2 ns [4.2 ns, 4.2 ns] |
| la_stack_inf_norm | 202.4 ns [202.2 ns, 202.7 ns] | 128.9 ns [128.6 ns, 129.1 ns] | -36.3% | faster point estimate; marginal CIs separated | 1.57x | 159.6 ns [157.6 ns, 161.3 ns] | 163.8 ns [163.5 ns, 164.1 ns] |

### D=64

| Benchmark | v0.4.3 (point + CI) | Latest (point + CI) | Point-estimate change | CI relation | Point-estimate ratio | v0.4.3 nalgebra | v0.4.3 faer |
|-----------|-------:|-------:|-------:|:-----------|--------:|-------:|-------:|
| la_stack_det_via_lu | 15.18 µs [15.16 µs, 15.21 µs] | 15.60 µs [15.58 µs, 15.65 µs] | +2.7% | slower point estimate; marginal CIs separated | 0.97x | 15.01 µs [14.78 µs, 15.10 µs] | 10.70 µs [10.68 µs, 10.71 µs] |
| la_stack_det | 15.15 µs [15.13 µs, 15.17 µs] | 15.48 µs [15.44 µs, 15.51 µs] | +2.2% | slower point estimate; marginal CIs separated | 0.98x | — | — |
| la_stack_lu | 15.34 µs [15.31 µs, 15.37 µs] | 14.66 µs [14.64 µs, 14.69 µs] | -4.4% | faster point estimate; marginal CIs separated | 1.05x | 13.67 µs [13.65 µs, 13.70 µs] | 10.50 µs [10.48 µs, 10.56 µs] |
| la_stack_ldlt | 21.53 µs [21.49 µs, 21.58 µs] | 20.42 µs [20.39 µs, 20.49 µs] | -5.2% | faster point estimate; marginal CIs separated | 1.05x | 11.41 µs [10.97 µs, 11.50 µs] | 8.86 µs [8.85 µs, 8.87 µs] |
| la_stack_lu_solve | 17.48 µs [17.44 µs, 17.50 µs] | 17.24 µs [17.21 µs, 17.41 µs] | -1.3% | faster point estimate; marginal CIs separated | 1.01x | 14.84 µs [14.80 µs, 14.86 µs] | 12.20 µs [12.18 µs, 12.23 µs] |
| la_stack_ldlt_solve | 24.60 µs [24.49 µs, 24.76 µs] | 22.70 µs [22.63 µs, 22.84 µs] | -7.7% | faster point estimate; marginal CIs separated | 1.08x | 14.06 µs [14.02 µs, 14.10 µs] | 10.12 µs [10.11 µs, 10.12 µs] |
| la_stack_solve_from_lu | 2.64 µs [2.63 µs, 2.65 µs] | 2.55 µs [2.54 µs, 2.56 µs] | -3.4% | faster point estimate; marginal CIs separated | 1.04x | 786.0 ns [778.1 ns, 788.8 ns] | 1.66 µs [1.66 µs, 1.66 µs] |
| la_stack_solve_from_ldlt | 1.10 µs [1.10 µs, 1.10 µs] | 1.07 µs [1.07 µs, 1.07 µs] | -2.7% | faster point estimate; marginal CIs separated | 1.03x | 1.29 µs [1.29 µs, 1.29 µs] | 1.25 µs [1.24 µs, 1.25 µs] |
| la_stack_det_from_lu | 27.0 ns [27.0 ns, 27.1 ns] | 31.9 ns [31.8 ns, 32.5 ns] | +17.9% | slower point estimate; marginal CIs separated | 0.85x | 8.7 ns [8.7 ns, 8.7 ns] | 21.6 ns [21.6 ns, 21.7 ns] |
| la_stack_det_from_ldlt | 27.5 ns [27.0 ns, 28.7 ns] | 31.6 ns [31.6 ns, 31.6 ns] | +14.8% | slower point estimate; marginal CIs separated | 0.87x | 8.5 ns [8.5 ns, 8.5 ns] | 21.0 ns [20.9 ns, 21.1 ns] |
| la_stack_dot | 38.4 ns [38.4 ns, 38.5 ns] | 11.0 ns [10.9 ns, 11.0 ns] | -71.5% | faster point estimate; marginal CIs separated | 3.51x | 9.0 ns [9.0 ns, 9.0 ns] | 36.6 ns [36.6 ns, 36.7 ns] |
| la_stack_norm2_sq | 38.2 ns [38.2 ns, 38.3 ns] | 10.8 ns [10.8 ns, 10.8 ns] | -71.8% | faster point estimate; marginal CIs separated | 3.54x | 7.4 ns [7.4 ns, 7.4 ns] | 6.2 ns [6.2 ns, 6.2 ns] |
| la_stack_inf_norm | 1.92 µs [1.92 µs, 1.92 µs] | 621.3 ns [620.2 ns, 621.7 ns] | -67.6% | faster point estimate; marginal CIs separated | 3.09x | 1.14 µs [1.13 µs, 1.15 µs] | 1.57 µs [1.57 µs, 1.57 µs] |

## How to Update

Local performance reports are generated in isolated temporary worktrees:

```bash
# Local development: compare the current tree with the latest release
just performance-local

# Release PR: update docs/PERFORMANCE.md and archive the previous report
just performance-release

# GitHub Actions release assets
just performance-github-assets

# Explicit repair
just performance-release <current-tag> <previous-tag>
```

`just performance-local` writes `target/bench-reports/performance.md`.
`just performance-github-assets` writes `target/bench-reports/github-assets-performance.md`.

Older curated release-to-release reports are archived in `docs/archive/performance/`.

See `docs/BENCHMARKING.md` for the full comparison workflow.
