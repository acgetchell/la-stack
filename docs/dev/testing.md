# Testing Guidance

Test-design details for [AGENTS.md](../../AGENTS.md). Command selection belongs
in the [contributor validation workflow](../../CONTRIBUTING.md#validation-workflow).

## Contents

- [Scientific coverage](#scientific-coverage)
- [Dimension coverage](#dimension-coverage)
- [Focused execution](#focused-execution)
- [Benchmark evidence](#benchmark-evidence)

## Scientific coverage

- Unit tests cover known values, error paths, and dimension-generic behavior.
- Match exact error variants, typed reasons, origins, locations, and structured
  fields. Never replace an unexpected error with a numeric sentinel or assert
  only `is_err()`.
- Property tests under `tests/proptest_*.rs` verify algebraic invariants such
  as round trips, residuals, and sign agreement, rather than merely checking
  that an operation does not panic.
- Near-singular, large-entry, and Hilbert-style ill-conditioned inputs accompany
  well-conditioned cases in both tests and benchmarks.
- When a public API has a fast filter and an exact fallback for the same
  question, a property test verifies agreement wherever both are defined.

## Dimension coverage

Dimension-generic code must cover D=2 through D=5 whenever possible. Use a macro
that accepts a dimension literal and generates the corresponding tests:

```rust
macro_rules! gen_tests {
    ($d:literal) => {
        paste! {
            #[test]
            fn [<test_foo_ $d d>]() {
                // Call a const-generic helper and assert the result.
            }
        }
    };
}

gen_tests!(2);
gen_tests!(3);
gen_tests!(4);
gen_tests!(5);
```

Keep the macro body thin. Shared setup and assertions belong in readable,
independently testable const-generic helper functions.

Existing patterns include:

- [`src/matrix.rs`](../../src/matrix.rs): `gen_matrix_tests!`.
- [`src/lu.rs`](../../src/lu.rs): `gen_pivoting_solve_and_det_tests!` and
  `gen_tridiagonal_smoke_solve_and_det_tests!`.
- [`src/ldlt.rs`](../../src/ldlt.rs): `gen_ldlt_identity_tests!` and
  `gen_ldlt_diagonal_tests!`.
- [`src/exact.rs`](../../src/exact.rs): `gen_exact_identity_tests!`,
  `gen_det_exact_f64_agrees_with_det_direct!`, `gen_solve_exact_tests!`, and
  `gen_solve_exact_f64_agrees_with_lu!`.

Single-dimension tests are appropriate for inherently dimension-specific known
values or errors requiring a particular layout; they need not be macro-generated.

## Focused execution

Select tests that exercise the changed behavior, then run the required final
validation for the affected surfaces. Useful targeted commands include:

```bash
cargo nextest run solve_2x2_basic
cargo nextest run -- --exact lu::tests::solve_2x2_basic
cargo nextest run --test proptest_matrix
```

Doctests run separately from nextest. README mirrors, guide examples, and
feature-dependent doctests follow [Documentation guidance](docs.md).

## Benchmark evidence

Performance remains subordinate to mathematical correctness, API stability,
composability, and clarity. When those conflict with speed, re-scope the problem
rather than weaken an invariant.

- Use comparable before-and-after measurements from the same representative
  command, inputs, features, and environment. Choose `bench-vs-linalg` for
  nalgebra/faer comparisons or `bench-exact` for exact arithmetic as appropriate.
- Prefer Criterion `bencher.iter` for nanosecond-scale fixed-size kernels so
  the complete public operation is measured symmetrically.
- Use `iter_batched` only when setup is explicitly outside the estimand, its
  exclusion is comparable across implementations, and a same-binary comparison
  shows batching does not materially distort the result.
- Preserve provenance. Point-estimate ratios are descriptive; marginal
  Criterion interval separation is not a paired confidence interval for a change.
- Runs that violate documented invariants are invalid performance evidence.
  Exact Criterion helpers accept only `ValidatedExactInput`; independent oracle
  validation occurs outside timed closures before that type is constructed.

The [Benchmarking guide](../BENCHMARKING.md) owns command matrices, comparison
methodology, input validation, report promotion, and release-artifact provenance.
