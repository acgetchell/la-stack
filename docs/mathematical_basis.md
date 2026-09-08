# Mathematical basis

## Contents

- [Introduction](#introduction)
- [Choosing an API](#choosing-an-api)
- [Geometry relationship and scope](#geometry-relationship-and-scope)
- [Represented values and arithmetic model](#represented-values-and-arithmetic-model)
- [Tolerances and typed errors](#tolerances-and-typed-errors)
- [Linear algebra algorithms](#linear-algebra-algorithms)
  - [Certified fixed-vector reductions](#certified-fixed-vector-reductions)
  - [Determinants and certified sign filtering](#determinants-and-certified-sign-filtering)
    - [Derivation of the returned determinant bound](#derivation-of-the-returned-determinant-bound)
  - [Exact arithmetic over binary64 inputs](#exact-arithmetic-over-binary64-inputs)
  - [Exact arithmetic over rational inputs](#exact-arithmetic-over-rational-inputs)
  - [Exact-to-binary64 conversion](#exact-to-binary64-conversion)
  - [Gram matrices and geometric measures](#gram-matrices-and-geometric-measures)
  - [LDLT without pivoting](#ldlt-without-pivoting)
  - [LU with partial pivoting](#lu-with-partial-pivoting)
  - [Outward-rounded interval expressions](#outward-rounded-interval-expressions)
  - [Scaled determinant products](#scaled-determinant-products)
  - [Scaled Euclidean vector norm](#scaled-euclidean-vector-norm)

## Introduction

`la-stack` provides fixed-dimension numerical linear algebra over two deliberate
point-value input domains plus one bounded layer. `Matrix<D>` and `Vector<D>`
store finite IEEE 754 binary64 values; their default algorithms operate in
binary64 and are therefore approximate. `Interval` and `IntervalMatrix<D>`
enclose exact-real expression values between outward-rounded finite binary64
bounds when interval operations are used throughout expression assembly.
Lifting an already-rounded `Matrix` encloses only its stored values. The
optional `exact` feature can either lift stored binary64 values to
the exact rational numbers they represent or accept caller-supplied
`BigRational` values through `RationalMatrix<D>` and `RationalVector<D>`.

This document separates three questions that are easy to conflate:

1. Which mathematical factorization or determinant identity is being used?
2. Which part of the computation is rounded in binary64?
3. What does a tolerance or typed error prove about the supplied values?

Reference numbers point to [REFERENCES.md](../REFERENCES.md).

## Choosing an API

| Need | API | Important boundary |
|------|-----|--------------------|
| General floating solve | `lu(tol)` then `solve` | Approximate; absolute pivot policy |
| Positive-definite floating solve | `ldlt(tol)` then `solve` | Exact symmetry; computed pivots must exceed tolerance; success is not a certificate |
| Floating determinant, any `D` | `det` | No certified bound; zero is not exact singularity |
| `D ≤ 4` error-bounded determinant/sign test | `det_direct_with_errbound` | Sign is certified when estimate magnitude exceeds bound; otherwise inconclusive |
| Derived-expression determinant sign through `D ≤ 7` | `IntervalMatrix::det_sign` | Outward-rounded proof; overlap with zero is explicitly inconclusive |
| Exact determinant sign | `det_sign_exact` | Exact for stored binary64 entries |
| Exact determinant value or solve | `det_exact`, `solve_exact` | Exact for represented inputs |
| Exact operations over preassembled rationals | `RationalMatrix::det_sign`, `det`, `solve` | No intermediate binary64 reconstruction |
| Binary64 output from an exact result | Strict or rounded conversions | Strict conversion forbids rounding |

## Geometry relationship and scope

Orientation, in-sphere, and related geometric predicates can be reduced to
determinant signs, which is why an adaptive exact sign is useful near degeneracy
\[8\]. `la-stack` supplies both point-matrix and bounded-expression determinant
primitives; callers still own problem-specific matrix assembly and semantic
classification. The crate originated to support
[`delaunay`](https://crates.io/crates/delaunay), but its matrix, factorization,
and exact-arithmetic APIs are general numerical infrastructure.

The deliberate anti-goals are dynamically sized or rectangular matrices,
sparse storage, broad decomposition coverage, alternate floating scalar
families, and GPU-, parallel-, BLAS-, or LAPACK-scale throughput. Those problems
need different storage models and numerical policies rather than extensions to
this small fixed-dimension design.

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
represent rational zero. Exact methods on `Matrix` and `Vector` exploit this
representation, but their exactness begins only after binary64 construction:
they cannot recover information lost when a decimal value or an earlier
computation was rounded to `f64`. Caller-supplied `RationalMatrix` and
`RationalVector` values instead enter the exact domain directly and do not pass
through binary64.

`Matrix`, `Vector`, `IntervalMatrix`, `Lu`, and `Ldlt` use inline fixed-size storage.
Arbitrary-precision `BigInt` and `BigRational` values allocate when the `exact` feature is
used. Const-generic `D` is not itself a mathematical dimension limit; “small,
fixed dimension” is the intended performance scope. `try_with_stack_matrix!` is
separately limited to `D = 0..=MAX_STACK_MATRIX_DISPATCH_DIM` (currently 7)
because it enumerates concrete stack types.

Except for the fixed-vector reduction and determinant filters described below,
the floating-point APIs do not provide certified forward, backward, or absolute
error bounds. This includes plain `Vector::dot`, Euclidean and squared norms,
matrix norms, factorizations, and solves. Some kernels use FMA to reduce rounding
steps, but that does not make them exact.

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

## Linear algebra algorithms

### Certified fixed-vector reductions

`Vector::dot_with_errbound()` binds the ordinary left-to-right FMA estimate to
a certified absolute error bound for the exact-real dot product of the stored
binary64 inputs. Starting from `s₀ = 0`, its arithmetic tree is

```text
sᵢ₊₁ = fma(leftᵢ, rightᵢ, sᵢ).
```

`Vector::dot_difference_with_errbound()` targets the exact-real expression

```text
Σᵢ axisᵢ(leftᵢ - rightᵢ)
```

without rounding coordinate differences first. For each coordinate it applies
`fma(axisᵢ, leftᵢ, s)` followed by `fma(-axisᵢ, rightᵢ, s)`, giving a specified
`2D`-event tree.

Let `u = 2^-53` be binary64 unit roundoff and
`γₙ = n·u / (1 - n·u)`. When every estimate FMA result is normal or an exact
zero, standard floating-point reduction analysis gives [9-11]

```text
|estimate - exact value| ≤ γₙ Σⱼ |aⱼbⱼ|,
```

with `n = D` for a dot product and `n = 2D` for the affine difference. The
implementation constructs an upper bound on the magnitude sum: exact
integer-significand comparison determines whether each rounded product must move
to its next representable value, and every positive accumulation is rounded
upward. The division forming `γₙ` and its final multiplication are also rounded
upward. Magnitude-ordered `FastTwoSum` then selects finite outward endpoints for
`estimate ± bound` \[17\].

The relative-error argument is not used across gradual underflow. A nonzero
product or estimate FMA in the subnormal range, an invalid `γₙ`, or finite-range
exhaustion in proof-only arithmetic returns `Ok(None)`. This is inconclusive
evidence, not equality. A non-finite estimate FMA returns `LaError::NonFinite`
with the failing reduction index and operation. A returned
`ScalarWithErrorBound` can certify a sign or threshold comparison from its
lower/upper endpoints; overlap requires an exact fallback that reconstructs the
same expression over the original inputs. The certificate is a roundoff bound
for this arithmetic tree, not a numerical tolerance chosen by the caller.

### Determinants and certified sign filtering

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

#### Derivation of the returned determinant bound

The implementation uses a rounded permanent `p_hat`, so its actual returned
bound is `B = fl(c_D × p_hat)`. The following argument includes that rounding.
Assume IEEE 754 round-to-nearest, ties-to-even and gradual underflow, with every
rounded operation finite and normal or an exact zero \[9-11\].
For this proof use `ε = 2^-52`, twice the unit roundoff `u = 2^-53`.
Then each nonzero operation satisfies `fl(t) = t(1 + δ)`, `|δ| ≤ ε`.
This allowance also covers an exact `t` just below the smallest normal that
rounds up to a normal result: its relative error is at most `u / (1 - u) < ε`.
Exact zeros require no error term. Per-operation tracking conservatively accepts
only structural zeros; the grid shortcut below also proves cancellation zeros
exact. Subnormal rounded results are excluded by the filter.

Expanding the implemented arithmetic tree expresses each determinant monomial
as its exact signed value times at most `k_D` factors of the form `(1 + δ)`.
FMA contributes one rounding factor, including when its two terms cancel.
The same count applies to the unsigned monomials of the permanent:

| Dimension | Longest determinant path | Longest permanent path | `k_D` |
|-----------|--------------------------|------------------------|-------|
| 2 | One multiply, one FMA | One multiply, one addition | 2 |
| 3 | A 2×2 minor, then three outer operations | An absolute 2×2 sum, then three outer operations | 5 |
| 4 | A 3×3 cofactor, then four outer operations | An absolute 3×3 sum, then four outer operations | 9 |

Sharing the six 2×2 minors in the dense D=4 path changes reuse, not the number
of rounding factors per monomial. Skipping zero coefficients in sparse paths
can only remove terms and operations. The standard product-of-errors bound,
`|Π(1 + δ_j) - 1| ≤ γ_k`, with `γ_k = kε / (1 - kε)`, therefore gives \[11\]

```text
|d_hat - det(ι(A))| ≤ γ_k p
p_hat ≥ (1 - γ_k) p,
```

where `p = p(|A|)` is exact. No independence of the rounding errors is assumed.
Accounting also for rounding the final multiplication yields

```text
B ≥ (1 - ε) c_D p_hat ≥ (1 - ε) c_D (1 - γ_k) p.
```

Consequently, a sufficient condition for `B ≥ |d_hat - det(ι(A))|` is

```text
c_D ≥ γ_k / ((1 - ε)(1 - γ_k))
    = kε / ((1 - ε)(1 - 2kε)).
```

All three stored `c_D` values are exactly representable in binary64. For
`ε = 2^-52` and `k ≤ 9`, the denominator `(1 - ε)(1 - 2kε)` is greater than
`3/4`. The required coefficient is thus less than `(4k/3)ε`, which is at most
the respective linear term `3ε`, `8ε`, or `12ε`. The positive quadratic terms
provide additional margin. This proves the returned bound, including downward
rounding of the permanent and bound, for the stated domain. If `p = 0`, all
Leibniz terms vanish and the same inequalities give zero error and bound.

The fast underflow check also has a simple grid argument. When every nonzero
input has magnitude at least `2^-16`, every entry is an integer multiple of
`2^-68`. The homogeneous degree-`j` arithmetic nodes preserve the grid
`2^(-68j)` under rounding. Through degree four, a nonzero result is therefore
at least `2^-272`; since each coefficient is a multiple of `2^-100`, a nonzero
final bound is at least `2^-372`. These are far above the smallest normal
`2^-1022`. Smaller inputs use per-operation underflow tracking instead.
Overflow remains a typed failure in either path.

The independent rational Leibniz checks in `tests/proptest_exact.rs` exercise
the error inequality for D=2–4. They support regression detection; the argument
above supplies the bound for every input satisfying the arithmetic assumptions.

### Exact arithmetic over binary64 inputs

The `exact` feature decomposes each stored entry into an integer mantissa and a
power of two. The entries are scaled to integer matrices without changing their
represented rational values \[9-10\].

Exact determinants use direct `BigInt` expansions for `D ≤ 4` and fraction-free
Bareiss elimination for `D ≥ 5` \[7\]. Exact solves apply Bareiss updates to an
integer augmented system, then use `BigRational` for back-substitution. Matrix
and right-hand-side scales start independently. Writing the selected exponents
as `s_A` and `s_b`, the integer forms satisfy `A = 2ˢᴬ · A_int` and
`b = 2ˢᵇ · b_int` \[9\]. Gaps of at most 64 bits use the lower shared scale,
while larger gaps remain separate; multiplying the integer system's solution by
the exact factor `2^(s_b − s_A)` preserves `A x = b`. First-nonzero pivoting is
sufficient for correctness in exact arithmetic, although pivot choice can still
affect computational cost.

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

### Exact arithmetic over rational inputs

`RationalMatrix<D>` and `RationalVector<D>` are a separate input domain for
coefficients assembled exactly before linear algebra begins. Their constructors
reject zero denominators and store each accepted quotient in lowest terms with
a positive denominator, so equivalent raw representations have identical
storage and denominator-clearing cost. For each matrix row `i`, the
implementation selects a positive least common multiple `sᵢ` of
the rational denominators and forms the integer row `A_int[i, :] = sᵢ A[i, :]`.
Therefore

```text
sign(det(A_int)) = sign(det(A))
det(A) = det(A_int) / ∏ᵢ sᵢ.
```

The sign path reads the `BigInt` determinant sign directly and never constructs
the rational determinant value. Solves include the corresponding right-hand
side denominator in `sᵢ`, so each augmented row is multiplied by the same
positive factor and the solution set is unchanged. The integer determinant and
solve states use the same direct-expansion/Bareiss backend as exact operations
over binary64 inputs, followed by rational back-substitution for solves \[7\].

The rational types are const-generic and shape-safe after construction. The
`try_with_rational_matrix!` helper provides explicit runtime dispatch through
D=8 without unstable generic const expressions. Conversion to binary64 remains
a separate strict or explicitly rounded `ExactF64Conversion` operation.

### Exact-to-binary64 conversion

For strict conversion, a canonical nonzero rational must have a power-of-two
denominator. Removing powers of two from its numerator produces
`x = sign × m × 2^e`, with `m` positive and odd. Let `b` be the bit length of
`m`. Exact finite binary64 representation is possible precisely when \[9-10\]

```text
b ≤ 53,    e ≥ -1074,    e + b - 1 ≤ 1023.
```

These conditions bound significand precision, the least representable bit,
and the highest exponent. Subnormal output stores `m × 2^(e + 1074)` in the
fraction field; normal output aligns the significand to 53 bits and stores the
biased exponent. Zero is handled separately. Raw `BigRational` inputs with
common factors or negative denominators are normalized as needed before this
decision; a zero denominator is rejected. Canonical storage allows
`RationalVector` conversion to skip that repeated normalization.

The integer-and-exponent determinant path rounds directly from `BigInt` digits.
It retains 53 significand bits for normal results, or the bits on the `2^-1074`
grid for subnormals. A discarded guard bit triggers an increment only if some
lower discarded bit is set (the sticky condition) or the retained integer is
odd. This implements nearest-even ties without constructing a rational
denominator. A carry can promote a subnormal to normal, advance the exponent,
or overflow. At `f64::MAX + 2^970`, the tie rounds to infinity; at magnitude
`2^-1075`, the tie rounds to signed zero \[9-10\].

For an already-computed `BigRational` value, rounded conversion delegates to
`num-rational`'s `ToPrimitive::to_f64` and checks for finite output. The strict
path uses that rounded result only when needed to distinguish
`RequiresRounding` from `NotFinite`; it never returns a rounded value as an exact
conversion. Neither conversion reruns determinant or solve elimination.

### Gram matrices and geometric measures

For `M` vectors in `N` coordinates, let `V` contain the vectors as rows. Their
Gram matrix is `G = V Vᵀ`, with `G[i,j] = v_i · v_j`. In exact real arithmetic
`G` is positive semidefinite because `zᵀGz = ||Vᵀz||² ≥ 0`; it is positive
definite exactly when the vectors are linearly independent.

For `M ≤ N`, `det(G)` is the squared M-dimensional volume of the parallelotope
spanned by those vectors. When the vectors are edges from a common simplex
vertex, the simplex volume is `sqrt(det(G)) / M!` \[16\]. These identities
apply to lower-dimensional simplices embedded in a larger coordinate space.

`gram_matrix` evaluates each upper-triangle dot product once using
`Vector::dot`'s left-to-right FMA reduction and copies it to the other triangle.
The returned binary64 matrix is therefore exactly symmetric, but rounding can
destroy positive semidefiniteness or rank. Construction does not certify linear
independence or volume accuracy and provides no absolute rounding-error bound.
A dot-product overflow remains a typed `LaError::NonFinite` failure.

For nonempty `V` with full row rank, `κ₂(G) = κ₂(V)²`: forming a Gram matrix
squares the spectral condition number \[11-12\]. Nearly dependent vectors require care
when interpreting a computed determinant or passing the matrix to `ldlt`.
Its exact-symmetry and computed-pivot checks remain in force. The
[Gram references](../REFERENCES.md#gram-matrices-and-geometric-measures)
provide the volume identity and numerical background.

<!-- Preserve links to the former factorization grouping. -->
<a name="floating-point-factorizations"></a>

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

### Outward-rounded interval expressions

`Interval` owns the invariant `-∞ < lower ≤ upper < +∞`. Its public constructor
rejects non-finite or inverted bounds, its fields are private, and every
arithmetic operation either returns another valid enclosure or a typed range
failure. Both signed-zero inputs represent exact real zero and are canonicalized
to `+0.0`; finite subnormal endpoints remain valid.

Point construction introduces no width. Exact-real subtraction and interval
addition use an error-free, magnitude-ordered `FastTwoSum` residual to determine
whether the rounded result is exact or which adjacent binary64 value is required
for the outward endpoint. Ordering the operands prevents intermediate overflow
whenever the rounded sum is finite \[17\]. Multiplication decomposes each nonzero
binary64 operand into its exact integer significand and power of two, compares
the exact 106-bit significand
product with the rounded result, and widens only in the required direction.
This comparison also handles products that underflow to zero: a positive result
is enclosed by `[0, f64::from_bits(1)]`, and a negative result by the mirrored
interval. Squaring uses multiplication bounds but gives every interval spanning
zero the exact lower bound zero. These guarantees rely on IEEE-754 binary64
round-to-nearest, ties-to-even, and gradual underflow \[9-11\]. IEEE 1788
provides the broader standardized interval arithmetic model \[14\]; this crate's
deliberately smaller, undecorated surface does not claim conformance.

If the exact result lies outside `[-f64::MAX, f64::MAX]`, no interval with finite
binary64 endpoints can contain it. The operation then returns
`LaError::IntervalRangeExhausted` with the responsible interval
`ArithmeticOperation`; it never stores infinity as an interval bound.

For D≤7, `IntervalMatrix::det()` evaluates the Leibniz determinant with subset
dynamic programming. A state for each column subset stores the determinant
enclosure for the corresponding leading-row minor, requiring 128 inline states
at D=7 and `D × 2^(D-1)` interval products. The expansion performs no division,
so a pivot interval containing zero cannot make the algorithm unsound.
`det_sign()` classifies a strictly positive or negative enclosure accordingly,
returns `Zero` only for the singleton `[0, 0]`, and returns `Inconclusive` for
every other overlap with zero. Inconclusive evidence is not a singularity
classification.

Lifting a completed `Matrix` creates point intervals for its stored values. It
does not recover rounding from earlier subtraction, dot products, or lifted-norm
construction. Robust callers instead assemble those derived coefficients with
interval operations, use `IntervalMatrix::det_sign()` as a fast proof, and
rebuild the expression in `RationalMatrix` or another exact representation when
the filter is inconclusive or loses range.

### Scaled determinant products

Both factorizations compute their determinant from the stored diagonal factors,
with the row-permutation sign included for LU. Direct multiplication is retained
while every running product is finite and normal. Otherwise `ScaledProduct`
replays all factors as a sign, a normalized mantissa, and an integer exponent.
IEEE 754 bit decomposition gives each nonzero factor as `m × 2^e`, including
subnormals, with `1 ≤ m < 2` \[9-10\]. Multiplication of two such mantissas
rounds below 4, so at most one exact division by two restores the invariant.
The separate exponent records the corresponding scale without range loss.

The most recent factor stays pending until finalization. If the final exponent
is in the subnormal range, the two remaining operands are scaled by exact
powers of two before multiplication. This places their product directly on the
destination's subnormal grid, avoiding a normal-mantissa rounding followed by a
second rounding to that coarser grid. The final result can be signed zero;
overflow returns `LaError::NonFinite`.

This is range management for a product of already-rounded factor diagonals.
Earlier mantissa multiplications still round, and factorization error remains.
Scaling guarantees neither an exact product nor correct rounding of the exact
determinant. No certified absolute error bound is provided.

### Scaled Euclidean vector norm

`Vector::norm` computes `sqrt(Σᵢ xᵢ²)` with a left-to-right scaled
sum-of-squares recurrence. Its state represents the accumulated squared norm as
`scale² × scaled_sum`. For each nonzero `|xᵢ|`, either `|xᵢ| / scale` is
squared and accumulated, or a larger `|xᵢ|` becomes the new scale and the old
sum is rescaled. Every squared ratio is therefore at most one, which avoids raw
square overflow and scales all-subnormal vectors into a safe range \[15\].

Division, FMA, square root, and final rescaling still round in binary64. The
method is deterministic for a fixed coordinate order but does not generally
claim correct rounding or publish a certified error bound.

Near the upper range, a fixed-size integer accumulator avoids both false and
hidden overflow from the rounded recurrence. Every coordinate square is an
integer multiple of `2^-2148` and is below `2^2048`, so 4196 bits plus
`usize::BITS` carry bits suffice for every representable vector length. The
fallback sums these integer squares exactly, retaining even the smallest
subnormal square. It compares against squared binary64 rounding midpoints to
return the nearest norm, ties to even, without a floating square root. The
overflow midpoint is `f64::MAX + 2^970`; equality rounds to infinity \[9-10\].

The fallback is selected from the largest coordinate magnitude, independently
of the rounded norm. With `b = bit_length(D)`, a largest magnitude at most
`2^(1023-b)` gives the conservative L1 upper bound `D × scale < 2^1023`, so
the ordinary recurrence has ample overflow margin. Larger magnitudes use the
exact boundary calculation. Its fixed storage stays on the stack and requires
no `exact` feature.

A scalar `LaError::NonFinite` tagged with `ArithmeticOperation::VectorNorm`
therefore means the exact norm rounds to infinity. `Vector::norm_squared`
intentionally remains the direct FMA sum `Σᵢ xᵢ²` and may therefore fail even
when `norm` succeeds.
