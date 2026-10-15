# Scientific Validation Records

This directory records scientific validation outcomes, decisions, and known
limits. The [mathematical basis](../mathematical_basis.md) owns the numerical
contracts and derivations; the records here link those contracts to retained
evidence. A passing test corpus supplements a contract or proof, rather than
establishing a guarantee for every representable input.

## Contents

- [Arithmetic contracts and evidence](#arithmetic-contracts-and-evidence)
- [Decision records](#decision-records)
- [Record scope](#record-scope)

## Arithmetic contracts and evidence

The distinction between rounded values, certified bounds, and exact results is
part of the public contract. The table maps those distinctions to existing
documentation and tests; it is not a claim that every suite was rerun for the
latest decision record.

| Arithmetic path | Contract and limitation | Representative evidence |
|-----------------|-------------------------|-------------------------|
| Certified floating reductions and determinants | [Bounds apply under documented assumptions](../mathematical_basis.md#certified-fixed-vector-reductions); an inconclusive filter requires fallback | [Independent rational bound/sign oracles](../../tests/proptest_exact.rs) |
| Exact arithmetic over binary64 inputs | [Exact for the stored values](../mathematical_basis.md#exact-arithmetic-over-binary64-inputs); earlier rounding cannot be recovered | [Exact determinants, residuals, and mixed-scale properties](../../tests/proptest_exact.rs) |
| Exact arithmetic over rational inputs | [Coefficients enter the exact domain directly](../mathematical_basis.md#exact-arithmetic-over-rational-inputs) | [Rational Gaussian and residual oracles](../../tests/proptest_rational.rs), [structured solves](../../tests/rational_structured_solve.rs) |
| Exact-to-binary64 conversion | [Strict conversion rejects required rounding](../mathematical_basis.md#exact-to-binary64-conversion); explicit rounded conversion opts into ties-to-even | [Rounding, signed-zero, subnormal, and overflow boundaries](../../tests/exact_conversion_boundaries.rs) |
| Outward-rounded interval arithmetic | [Encloses the expression assembled with interval operations](../mathematical_basis.md#outward-rounded-interval-expressions); not an exact point value | [Rational containment and sign checks](../../tests/proptest_interval.rs) |
| Strict `f64` arithmetic | [Defined binary64 operations and deliberate FMA](../mathematical_basis.md#represented-values-and-arithmetic-model); ordinary rounded APIs are not exact or automatically certified | [Matrix properties](../../tests/proptest_matrix.rs), [integer vector oracles](../../tests/proptest_vector.rs), [compiler migration study](../archive/performance/studies/rust-1.99.md) |

`strict f64` describes the arithmetic model. `strict exact-to-f64 conversion`
describes whether an exact value can be returned without rounding. These are
different contracts: a strict floating-point algorithm can accumulate rounding
error, while a successful strict exact conversion preserves the exact result.

## Decision records

- [Algebraic floating-point rejection (#250)](algebraic-floating-point.md):
  retained timing benefits do not justify a demonstrated overflow-contract
  violation; keep strict arithmetic and end the evaluation.

## Record scope

Each record identifies the source revision, environment, independent reference,
executed checks, observed failures, decision, and unperformed work. Raw
performance samples and experimental source remain with their
[archived study](../archive/performance/studies/README.md) rather than being
copied here. Cross-compilation is distinct from native execution, and downstream
tests of unchanged code do not validate a rejected prototype.

Test-design policy belongs in [testing guidance](../dev/testing.md); current
validation commands belong in the
[contributor workflow](../../CONTRIBUTING.md#validation-workflow).
