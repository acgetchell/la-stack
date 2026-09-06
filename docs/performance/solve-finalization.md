# Keep existing solve finalization (#234)

**Decision:** retain the existing LU and LDLT result finalization. The safe
alternative provided no repeatable solve speedup at D=2, 3, 4, 5, 8, 16, 32,
or 64, so the prototype was removed.

The experiment measured solves using precomputed factors, with identical
inputs and Criterion settings for both implementations and a reverse-order
repeat. On Rust 1.98.0 / AArch64, the inspected compiler output already removed
the final finite-value scan through D=8. At larger dimensions where a scan
remained, the alternative still showed no repeatable overall improvement.
This supports keeping the simpler implementation; it does not prove that
finalization has zero cost in every compiler or calling context.

The retained tests in
[`tests/solve_finalization.rs`](../../tests/solve_finalization.rs) exercise
division, forward-substitution, and back-substitution overflow across all eight
dimensions, preserving the exact typed error and failing-step order.

Correctness tests can be rerun with:

```bash
cargo test --locked --test solve_finalization
```

The existing benchmark suite can measure current precomputed-factor solves:

```bash
just bench-vs-linalg 'la_stack_solve_from_(lu|ldlt)$'
```

Tests verify correctness; benchmarks measure performance on the current
toolchain and machine. The original prototype and measurement evidence are
local experimental artifacts rather than required repository files.
