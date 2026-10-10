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
  - [Unsigned vector angles](#unsigned-vector-angles)

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
| Unsigned vector angle in radians | `VectorAngle::angle`, `Vector::angle` | Finite nonzero vectors; rounded, without a certified bound |

## Geometry relationship and scope

Orientation, in-sphere, and related geometric predicates can be reduced to
determinant signs, which is why an adaptive exact sign is useful near degeneracy
\[8\]. `la-stack` supplies both point-matrix and bounded-expression determinant
primitives; callers still own problem-specific matrix assembly and semantic
classification. The crate originated to support
[`delaunay`](https://crates.io/crates/delaunay), but its matrix, factorization,
and exact-arithmetic APIs are general numerical infrastructure.
The borrowed-slice angle API also supports runtime ambient coordinate lengths
without changing the const-generic storage model. Spherical callers retain
radius validation and angle-to-arc-length conversion.

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
implementation constructs an upper bound on the magnitude sum. A first nonzero
product with a unit (`±1`) multiplier is exact; otherwise the rounded product
moves up one representable value. This conservative step avoids computing a
rounding residual for the magnitude bound. Subsequent products enter the bound
through an FMA followed by one upward step. This encloses the previous upper
bound plus the exact product; for the same incoming bound, monotonicity makes
this step no wider than separately rounding the product upward before adding
and widening. The first product's unconditional step can make the final bound
slightly wider than the former residual-based calculation. Each nonzero product
must still round to a normal value, preserving the relative-error model's
admission rule.
When the only nonzero product is normal and has a unit multiplier, the
implementation instead certifies a zero error bound: its FMA has a zero addend,
and all other FMAs add exact zero. Encountering another nonzero product clears
this special proof. This also permits an exact certificate at `f64::MAX` when
adding a conservative positive error bound would have exhausted the endpoint
range.
The division forming `γₙ` and its final multiplication are also rounded
upward. Magnitude-ordered `FastTwoSum` then selects finite outward endpoints for
`estimate ± bound` \[17\].

The certificate stores the estimate and bound; endpoints are derived when
requested. Construction first proves that both endpoints can remain finite:
the exact value `|estimate| + bound` dominates their magnitudes. If its rounded
value is below `f64::MAX`, even the next representable value is finite. If it
rounds to `f64::MAX`, a nonpositive `FastTwoSum` residual proves the exact sum
does not exceed that limit. An infinite rounded value or positive residual at
the limit makes the certificate unavailable. Endpoint accessors can therefore
remain infallible and use the same tight directed rounding as before.

Zero products leave the magnitude bound unchanged. While the proof is
available, the estimate is normal or positive zero: it starts at positive zero,
and exact cancellation of nonzero terms rounds to positive zero. Adding either
signed zero then preserves every estimate bit, so that FMA can be omitted.
After proof loss, the FMA still executes: underflow may have produced negative
zero, whose sign a subsequent zero product can change. Nonzero FMAs keep their
original order, and the returned estimate matches the specified FMA tree.

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
the rational determinant value. General solves include the corresponding right-hand
side denominator in `sᵢ`, so each augmented row is multiplied by the same
positive factor and the solution set is unchanged. The integer determinant and
solve states use the same direct-expansion/Bareiss backend as exact operations
over binary64 inputs, followed by rational back-substitution for solves \[7\].

Before clearing denominators, rational solves recognize row permutations of
upper or lower triangular matrices. For an upper triangular ordering, each
row must have a distinct first non-zero column. There are D rows and D columns,
so these columns form a permutation of `0..D`; assigning each row to its first
non-zero column gives a non-zero diagonal and zeros strictly below it. The
same argument using last non-zero columns gives a lower triangular ordering.
Only row indices are stored, and the right-hand side follows the same row map.

Back- or forward-substitution then computes
`x[k] = (b[row(k)] - Σⱼ A[row(k), j] x[j]) / A[row(k), k]`, with the sum over
already-solved columns \[11\]. The residual is accumulated as two arbitrary-precision
integers `N/Q`, initially the canonical RHS numerator and positive denominator.
For a non-zero product `a × x = P/R`, cross multiplication updates
`N ← N R - P Q` and `Q ← Q R`; equal denominators need only `N ← N - P`.
Both `Q` and `R` remain positive. Dividing by the non-zero pivot `p/q` constructs
and reduces `N q / (Q p)` once, so every stored solution component is canonical.
This avoids a GCD reduction after each product and subtraction without creating
unchecked rational values. Zero coefficients or solution components need no
multiplication. Recognition and substitution cost at most
`O(D²)` coefficient operations, with `O(D)` additional storage; bit complexity
still depends on integer component growth. Deferred reduction may increase
intermediate bit lengths relative to reducing every operation. Numerical conditioning introduces
no rounding error in this exact domain.

A zero row or repeated edge column rejects that ordering; if neither ordering
works, the unchanged Bareiss backend handles the general system and owns exact
singularity diagnostics. This includes singular triangular inputs. No Bareiss
update is skipped: even a zero elimination coefficient would still require its
pivot-ratio scaling. Empty systems retain their unique empty solution.

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
whenever the rounded sum is finite \[17\]. Multiplication uses the exact residual
`left.mul_add(right, -rounded)` of `TwoProductFMA` \[18\] when
`|rounded| ≥ 2^-968`. To justify this threshold, let `E` be the rounded product's
leading binary exponent. The exact product has at most 106 significant bits and
leading exponent at least `E-1`, allowing for a rounding carry. Its least bit
therefore has exponent at least `E-106 ≥ -1074`. The residual fits binary64
precision and cannot lose low bits below the subnormal range. Its sign selects
the same tight outward endpoint as exact integer comparison.
In this range a normal power-of-two operand makes the product exact by exponent
scaling, so no residual evaluation is needed.

Below that conservative threshold, multiplication decomposes each nonzero
operand into an integer significand and power of two and compares the exact
106-bit significand product with the rounded result. A zero FMA residual alone
would not prove exactness here, even for a normal rounded product.
This comparison also handles products that underflow to zero: a positive result
is enclosed by `[0, f64::from_bits(1)]`, and a negative result by the mirrored
interval. Squaring uses multiplication bounds but gives every interval spanning
zero the exact lower bound zero. Singleton squares and sums evaluate the shared
endpoint only once, preserving the same enclosure and error provenance.
These guarantees rely on IEEE-754 binary64
round-to-nearest, ties-to-even, and gradual underflow \[9-11\]. IEEE 1788
provides the broader standardized interval arithmetic model \[14\]; this crate's
deliberately smaller, undecorated surface does not claim conformance.

If the exact result lies outside `[-f64::MAX, f64::MAX]`, no interval with finite
binary64 endpoints can contain it. The operation then returns
`LaError::IntervalRangeExhausted` with the responsible interval
`ArithmeticOperation`; it never stores infinity as an interval bound.

For D≤7, `IntervalMatrix::det()` evaluates the Leibniz determinant with subset
dynamic programming. A state for each column subset stores the determinant
enclosure for the corresponding leading-row minor, requiring `max(2, 2^D)`
inline states and `D × 2^(D-1)` interval products for positive D. The const
dimension selects the workspace size without changing the expansion order.
The expansion performs no division,
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

### Unsigned vector angles

The `VectorAngle` extension trait provides `left.angle(right)` on borrowed,
equal-length coordinate slices;
`Vector<D>::angle(&other)` reuses the kernel with validated finite storage.
Both compute an unsigned angle in radians in `[0, π]`, in quadratic time and
constant auxiliary space, without allocation or optional dependencies.
Every positive ambient length is accepted, including Delaunay's lengths 3–6
and lengths beyond the matrix dispatch limit. Shape mismatch is checked
before coordinates; non-finite coordinates are checked left operand first,
then right, in index order. Two empty vectors return `EmptyVector`; an all-zero
operand returns `ZeroVector`, left before right. Non-finite input metadata
retains `VectorOperandEntry { operand, index }` and `NonFiniteOrigin::Input`.

For nonzero vectors `u, v`, Lagrange's identity gives

```text
w² = Σᵢ<ⱼ (uᵢvⱼ - uⱼvᵢ)² = ‖u‖² ‖v‖² - (u·v)²,
θ = atan2(w, u·v).
```

The exterior-product norm `w` is `‖u‖ ‖v‖ sin(θ)`; the dot is the same
positive norm factor times `cos(θ)`. `atan2` cancels that factor without forming
the norms and remains well-conditioned near both endpoints and `π/2`.
The implementation reduces the minors with `hypot`; it never subtracts rounded
Gram products or squares a tiny minor. For length one the exterior sum is empty
and the sign of the dot selects zero or `π`.

Each input has its own positive maximum magnitude `m`. Extract its binary
exponent `e = floor(log2(m))` and scale its coordinates by `2^(500-e)`, using
up to two power-of-two multiplications. The largest coordinate is then in
`[2^500, 2^501)`. This preserves significands without dividing by a rounded norm
or non-power-of-two maximum. Coordinates below roughly `2^-1574` relative to
their maximum may underflow during scaling; their total angular contribution
stays far below the output range even for addressable slice lengths.

Products of scaled coordinates stay below `2^1002`; differences stay below
`2^1003`. Each minor uses two FMA product residuals \[18\] plus a magnitude-ordered
`FastTwoSum` subtraction residual \[17\]. This retains small differences between
large, nearly equal products. All dot terms and completed minors are multiplied
by `2^-500` before reduction. Dot accumulation also tracks addition residuals.
For 64-bit slice lengths the accumulated dot and exterior norm remain below
roughly `2^567`, with ample overflow headroom. Minors contributing to a
representable tiny angle stay well above underflow at the enlarged product scale.

Kahan's norm-weighted half-angle formula \[19\], used in Delaunay's spherical
implementation, was considered and is retained as a benchmark control on valid
fixtures. General direction normalization can still erase representable tiny
angles: with `n = 2^51`, the planar vectors `[n,n-1]` and `[n+1,n]` have
determinant 1 and dot `2n²`, hence angle `atan(2^-103)`. Dividing each vector by
its maximum can make the directions identical in binary64. Compensated minors
avoid that information loss. They also preserve the least positive angle in
`[1,0]`, `[1,2^-1074]` and the separation of `[2,t]`, `[2,-t]` with
`t = 2^-1074`; separately dividing the latter small coordinates by 2 loses both.

For a positive dot and ratio `r = w / (u·v) ≤ 2^-27`, the implementation
returns the quotient directly. The alternating arctangent series gives
`0 ≤ r - atan(r) ≤ r³/3 ≤ r × 2^-55`; this approximation error is below
half an ulp for normal results and negligible for subnormal results. Division
avoids platform-dependent `atan2` underflow, observed in Windows CI for
`[3,t]`, `[3,-t]` and `[1,1,0]`, `[1,1,t]` with `t = 2^-1074`.
Both true angles round to the least positive subnormal. The bound here concerns
only replacing arctangent by its argument, not the preceding dot and exterior
reductions or their rounding. Larger ratios and nonpositive dots use `atan2`.

The exterior norm is nonnegative, so the result stays in `[0, π]` without clamping
or a half-angle intermediate. Identical inputs yield exactly positive zero;
opposite inputs yield the binary64 constant `π`, including signed-zero changes.
The price of retaining all minors is `N(N-1)/2` compensated differences and
`hypot` steps. This is intended for small dimensions, including Delaunay's 3–6,
while remaining allocation-free and without an artificial dimension limit.

This is approximate arithmetic over the stored inputs. There is no certified
absolute error bound, exact parallelism classification, or correct-rounding
promise. Products, compensated differences, reductions, and transcendental
evaluation still round; accumulated error depends on ambient length. Subnormal
outputs have coarse relative spacing, unrepresentable angles may round to
zero, and deviations from `π` smaller than its output spacing are lost.
Coordinate rescaling that itself rounds or underflows changes the input
direction. Rust's `hypot` and `atan2` can also vary across platforms.
The guarantees assume ordinary binary64 arithmetic with gradual underflow.
