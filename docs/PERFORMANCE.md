# Exact Arithmetic Performance

**la-stack** v0.3.0 · `0a8ce5b` (perf/63-custom-f64-to-bigrational) · 2026-04-10 23:57:42 UTC
**Statistic**: median

## Benchmark Results

Current performance snapshot (no baseline comparison).

| Group | Benchmark | Median | 95% CI |
|-------|-----------|-------:|-------:|
| exact_d2 | det | 0.6 ns | [0.6 ns, 0.6 ns] |
| exact_d2 | det_direct | 0.7 ns | [0.7 ns, 0.7 ns] |
| exact_d2 | det_exact | 3.93 µs | [3.89 µs, 3.97 µs] |
| exact_d2 | det_exact_f64 | 4.19 µs | [4.17 µs, 4.21 µs] |
| exact_d2 | det_sign_exact | 1.1 ns | [1.1 ns, 1.1 ns] |
| exact_d2 | solve_exact | 10.72 µs | [10.72 µs, 10.72 µs] |
| exact_d2 | solve_exact_f64 | 11.22 µs | [11.19 µs, 11.25 µs] |
| exact_d3 | det | 1.4 ns | [1.4 ns, 1.4 ns] |
| exact_d3 | det_direct | 5.1 ns | [5.1 ns, 5.1 ns] |
| exact_d3 | det_exact | 21.32 µs | [21.20 µs, 21.45 µs] |
| exact_d3 | det_exact_f64 | 20.90 µs | [20.78 µs, 21.01 µs] |
| exact_d3 | det_sign_exact | 3.3 ns | [3.2 ns, 3.3 ns] |
| exact_d3 | solve_exact | 46.73 µs | [46.53 µs, 46.92 µs] |
| exact_d3 | solve_exact_f64 | 51.25 µs | [50.61 µs, 51.90 µs] |
| exact_d4 | det | 2.3 ns | [2.3 ns, 2.3 ns] |
| exact_d4 | det_direct | 2.4 ns | [2.4 ns, 2.4 ns] |
| exact_d4 | det_exact | 71.86 µs | [71.40 µs, 72.33 µs] |
| exact_d4 | det_exact_f64 | 72.30 µs | [72.20 µs, 72.40 µs] |
| exact_d4 | det_sign_exact | 6.3 ns | [6.3 ns, 6.3 ns] |
| exact_d4 | solve_exact | 140.23 µs | [137.78 µs, 142.68 µs] |
| exact_d4 | solve_exact_f64 | 144.25 µs | [142.93 µs, 145.58 µs] |
| exact_d5 | det | 26.6 ns | [26.1 ns, 27.2 ns] |
| exact_d5 | det_direct | 2.2 ns | [2.2 ns, 2.2 ns] |
| exact_d5 | det_exact | 163.15 µs | [162.53 µs, 163.77 µs] |
| exact_d5 | det_exact_f64 | 176.47 µs | [176.10 µs, 176.84 µs] |
| exact_d5 | det_sign_exact | 176.38 µs | [173.58 µs, 179.18 µs] |
| exact_d5 | solve_exact | 328.97 µs | [325.89 µs, 332.05 µs] |
| exact_d5 | solve_exact_f64 | 333.97 µs | [333.03 µs, 334.91 µs] |
| exact_near_singular_3x3 | det_sign_exact | 12.71 µs | [12.63 µs, 12.78 µs] |
| exact_near_singular_3x3 | det_exact | 12.83 µs | [12.75 µs, 12.91 µs] |

## How to Update

```bash
# Save a baseline at the current release
just bench-save-baseline v0.3.0

# Compare current code against a saved baseline
just bench-compare v0.3.0

# Generate a snapshot without comparison
just bench-compare
```

See `docs/RELEASING.md` for the full release workflow.
