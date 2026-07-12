# Mathematical basis

`la-stack` provides fixed-dimension numerical linear algebra over finite IEEE
754 binary64 values. For a compile-time dimension `D`, `Matrix<D>` is a dense
`D × D` square matrix and `Vector<D>` is a length-`D` vector over the finite
binary64 set `F`. The default algorithms operate in binary64 and are therefore
approximate. The optional `exact` feature instead lifts each stored binary64
value to the exact rational number it represents.

This document separates three questions that are easy to conflate:

1. Which mathematical factorization or determinant identity is being used?
2. Which part of the computation is rounded in binary64?
3. What does a tolerance or typed error prove about the supplied values?

Reference numbers point to [REFERENCES.md](../REFERENCES.md).

## Represented values and arithmetic model

Public matrix and vector construction rejects NaN and infinity. Because the
backing arrays are private and mutation is validated, every stored entry is
finite. Arithmetic can still overflow, underflow, or round; computed non-finite
values are reported through `LaError::NonFinite` rather than stored in a public
matrix, vector, or factorization.

Every finite binary64 value is a dyadic rational with an exact representation

```text
x = (-1)^s m × 2^e,
```

for integer `m` and exponent `e` \[9-11\]. Both signed-zero bit patterns
represent rational zero. Exact APIs exploit this representation, but exactness
begins only after construction: they cannot recover information lost when a
decimal value or an earlier computation was rounded to `f64`.

`Matrix`, `Vector`, `Lu`, and `Ldlt` use inline fixed-size storage.
Arbitrary-precision `BigInt` and `BigRational` values allocate when the `exact` feature is
used. Const-generic `D` is not itself a mathematical dimension limit; “small,
fixed dimension” is the intended performance scope. `try_with_stack_matrix!` is
separately limited to `D = 0..=MAX_STACK_MATRIX_DISPATCH_DIM` (currently 7)
because it enumerates concrete stack types.

Except for the determinant filter described below, the floating-point APIs do
not provide certified forward, backward, or absolute error bounds. This includes
dot products, squared norms, matrix norms, factorizations, and solves. Some
kernels use FMA to reduce rounding steps, but that does not make them exact.

## Floating-point factorizations

### LU with partial pivoting

LU factorization targets

```text
P A = L U,
```

where `P` is a row permutation, `L` is unit lower triangular, and `U` is upper
triangular. At column `k`, the implementation selects the largest-magnitude
remaining entry in the active column. It rejects the factorization when that
magnitude is less than or equal to the caller's tolerance.

The tolerance is an absolute finite, non-negative threshold. It is not divided
by a matrix norm, so rescaling a system can change whether a pivot is accepted.
Partial pivoting is a practical stability strategy, not an unconditional
accuracy guarantee: worst-case element growth and typical behavior are distinct
questions \[1-3, 11-12\]. `Lu::solve` does not estimate conditioning or refine
the result.

`Lu::det` combines permutation parity with the product of the diagonal of `U`.
Scaled product accumulation avoids some premature overflow and underflow, but
the final binary64 determinant is still rounded. A returned zero or a numerical
`LaError::Singular` is therefore not proof that the represented matrix is
exactly singular.

### LDLT without pivoting

LDLT factorization targets

```text
A = L D Lᵀ,
```

with unit lower-triangular `L` and diagonal `D` \[4-5, 11-12\]. The input
must be exactly symmetric under binary64 comparison: every mirrored pair must
satisfy `A[i][j] == A[j][i]`. Signed zeros compare equal and are accepted.

A successful `Ldlt` requires every computed diagonal pivot to be positive and
greater than the caller's tolerance. An uncoupled computed zero or a positive
pivot at or below tolerance returns `LaError::Singular`; a negative pivot or
zero pivot with remaining coupling returns `LaError::NotPositiveSemidefinite`
with a typed violation.

This is not a pivoted symmetric-indefinite factorization such as Bunch-Kaufman
\[6, 11-12\]. Pivots are computed in binary64, so a singular represented matrix
can produce a small positive pivot above a low tolerance. Successful
factorization is therefore not an exact positive-definiteness certificate for
the stored matrix, much less for ideal values before binary64 input conversion.

## Determinants and certified sign filtering

`det_direct()` evaluates closed forms for `D = 0..=4`, with the empty-product
convention `det(Matrix::<0>) = 1`. `Matrix::det()` uses that path through `D = 4`
and zero-tolerance LU for `D ≥ 5`. The general `det()` result has no certified
roundoff bound, and its LU fallback can report numerical singularity even when
exact arithmetic over the stored entries would find a nonzero determinant.

For `D = 0` and `D = 1`, `det_direct_with_errbound()` returns the exact direct
determinant with a zero bound. For `D = 2..=4`, it can return a determinant
estimate and a conservative absolute bound. Let `ι(A)` denote the exact rational
lift of the stored entries and let

```text
p(|A|) = Σ_(σ ∈ S_D) Π_i |a_(i,σ(i))| = perm(|A|).
```

When every rounded intermediate in the implemented determinant and `p(|A|)`
evaluation trees is normal or an exact structural zero, the implementation uses

```text
|det_direct(A) - det(ι(A))| ≤ c_D × p(|A|),
```

with `ε = f64::EPSILON` and project-specific coefficients

```text
c_2 = 3ε + 16ε²
c_3 = 8ε + 64ε²
c_4 = 12ε + 128ε².
```

The analysis follows Shewchuk's adaptive-filter framework \[8\], while the
FMA evaluation trees and constants are derived for this crate rather than copied
from that source. IEEE 754 and standard floating-point error analysis provide
the arithmetic model \[9-11\].

If `|determinant| > absolute_error_bound`, the sign is certified. Otherwise the
filter is inconclusive; it does not prove that the matrix is singular. The API
returns `Ok(None)` for `D ≥ 5` or when gradual underflow could invalidate the
relative-error model. A non-finite computed determinant or bound returns
`LaError::NonFinite`.

## Exact arithmetic over binary64 inputs

The `exact` feature decomposes each stored entry into an integer mantissa and a
power of two. The entries are scaled to integer matrices without changing their
represented rational values \[9-10\].

Exact determinants use direct `BigInt` expansions for `D ≤ 4` and fraction-free
Bareiss elimination for `D ≥ 5` \[7\]. Exact solves apply Bareiss updates to an
integer augmented system, then use `BigRational` for back-substitution. Matrix
and right-hand-side scaling are tracked separately and reconciled with an exact
power-of-two factor. First-nonzero pivoting is sufficient for correctness in
exact arithmetic, although pivot choice can still affect computational cost.

`det_sign_exact()` first attempts the certified binary64 filter for `D ≤ 4` and
falls back to exact integer arithmetic when the filter is inconclusive. It
returns the exact determinant sign for every finite stored matrix. `det_exact()`
and `solve_exact()` return values exact for the stored `A` and `b`, subject to
their documented scale and singularity errors.

Strict `*_exact_f64` conversions succeed only when the exact result already has
an exactly representable finite binary64 value. `RequiresRounding` means a
finite output exists only after rounding; `NotFinite` means even the rounded
result cannot be finite. The explicit rounded conversions use round-to-nearest,
ties-to-even \[9-10\]. A nonzero exact value may consequently round to zero.

## Tolerances and typed errors

`Tolerance::try_new` accepts finite values greater than or equal to zero.
`DEFAULT_SINGULAR_TOL` is the absolute value `1e-12`. A tolerance is a rejection
policy for numerical pivots or diagnostics, not an error estimate and not a
condition-number threshold.

The error model keeps distinct mathematical conclusions separate:

- `LaError::Singular` distinguishes numerical rejection from exact singularity.
  Numerical context retains the factorization kind, observed pivot magnitude,
  and tolerance.
- `LaError::Asymmetric` reports the mirrored values that violate LDLT's exact
  symmetry precondition.
- `LaError::NotPositiveSemidefinite` records a negative pivot or a zero pivot
  with remaining coupling.
- `LaError::NonFinite` distinguishes invalid input locations from arithmetic
  operations that overflowed during computation.
- `LaError::Unrepresentable` distinguishes required rounding from the absence
  of any finite binary64 output.

`Matrix::is_symmetric` and `Matrix::first_asymmetry` are tolerance-based,
scale-aware diagnostics. They do not establish the exact mirrored equality
required by `Matrix::ldlt`.

## Choosing an API

| Need | API | Important boundary |
|------|-----|--------------------|
| General floating solve | `lu(tol)` then `solve` | Approximate; absolute pivot policy |
| Positive-definite floating solve | `ldlt(tol)` then `solve` | Exact symmetry; computed pivots must exceed tolerance; success is not a certificate |
| Floating determinant, any `D` | `det` | No certified bound; zero is not exact singularity |
| `D ≤ 4` error-bounded determinant/sign test | `det_direct_with_errbound` | Sign is certified when estimate magnitude exceeds bound; otherwise inconclusive |
| Exact determinant sign | `det_sign_exact` | Exact for stored binary64 entries |
| Exact determinant value or solve | `det_exact`, `solve_exact` | Exact for represented inputs |
| Binary64 output from an exact result | Strict or rounded conversions | Strict conversion forbids rounding |

## Geometry relationship and scope

Orientation, in-sphere, and related geometric predicates can be reduced to
determinant signs, which is why an adaptive exact sign is useful near degeneracy
\[8\]. `la-stack` supplies the determinant primitive; callers construct the
problem-specific predicate matrix and remain responsible for any rounding that
occurs during that construction. The crate originated to support
[`delaunay`](https://crates.io/crates/delaunay), but its matrix, factorization,
and exact-arithmetic APIs are general numerical infrastructure.

The deliberate anti-goals are dynamically sized or rectangular matrices,
sparse storage, broad decomposition coverage, alternate floating scalar
families, and GPU-, parallel-, BLAS-, or LAPACK-scale throughput. Those problems
need different storage models and numerical policies rather than extensions to
this small fixed-dimension design.
