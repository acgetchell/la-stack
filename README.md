# la-stack

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18158926.svg)](https://doi.org/10.5281/zenodo.18158926)
[![Crates.io](https://badgen.net/crates/v/la-stack)](https://crates.io/crates/la-stack)
[![Downloads](https://badgen.net/crates/d/la-stack)](https://crates.io/crates/la-stack)
[![License](https://badgen.net/github/license/acgetchell/la-stack)](https://github.com/acgetchell/la-stack/blob/v0.4.5/LICENSE)
[![Docs.rs](https://docs.rs/la-stack/badge.svg)](https://docs.rs/la-stack)
[![CI](https://github.com/acgetchell/la-stack/actions/workflows/ci.yml/badge.svg)](https://github.com/acgetchell/la-stack/actions/workflows/ci.yml)
[![rust-clippy analyze][clippy-badge]][clippy-workflow]
[![codecov](https://codecov.io/gh/acgetchell/la-stack/graph/badge.svg?token=4eKXa5QjuZ)](https://codecov.io/gh/acgetchell/la-stack)
[![Audit dependencies][audit-badge]][audit-workflow]

![la-stack](https://raw.githubusercontent.com/acgetchell/la-stack/main/docs/assets/la-stack.jpg)

Fast, stack-allocated linear algebra for fixed dimensions in Rust.

This crate grew from the need to support [`delaunay`](https://crates.io/crates/delaunay) with fast, stack-allocated linear algebra primitives and algorithms
while keeping the API intentionally small and explicit.

## 📐 Introduction

`la-stack` provides a handful of const-generic, stack-backed building blocks:

- `Vector<const D: usize>` for fixed-length `f64` vectors backed by `[f64; D]`
- `Matrix<const D: usize>` for fixed-size square `f64` matrices backed by `[[f64; D]; D]`
- `Interval` and `IntervalMatrix<const D: usize>` for outward-rounded,
  proof-bearing determinant filters through D=7
- `ScalarWithErrorBound` for proof-bearing fixed-vector dot products and
  affine differences over finite `f64` inputs
- `RationalVector<const D: usize>` and `RationalMatrix<const D: usize>` for
  exact rational inputs behind the optional `"exact"` feature
- `Lu<const D: usize>` for LU factorization with partial pivoting (solve + det)
- `Ldlt<const D: usize>` for no-pivot factorization intended for exactly
  symmetric positive-definite matrices (solve + det; typed pivot diagnostics)

## 🧮 Mathematical basis

`la-stack` operates on finite IEEE 754 binary64 values in small, fixed
dimensions. Its floating-point paths use LU with partial pivoting, LDLT without
pivoting for exactly symmetric positive-definite matrices, and closed-form
determinants through D=4. These results remain subject to conditioning and
binary64 rounding;
factorization tolerances are rejection thresholds, not accuracy guarantees. For
D≤4, direct determinants can be paired with a conservative absolute roundoff
bound when its range preconditions hold. Fixed-vector dot products and direct
affine differences can likewise return a paired estimate and certified absolute
roundoff bound without enabling arbitrary-precision dependencies.

Derived binary64 expressions can instead be assembled with `Interval`
subtraction, addition, multiplication, negation, and square. The resulting
`IntervalMatrix<D>` determinant sign is certified through D=7 when its enclosure
separates zero; the singleton `[0, 0]` also certifies exact zero. Every other
overlap with zero is explicitly inconclusive. This default-feature surface is
distinct from arbitrary-precision exact arithmetic.

With `features = ["exact"]`, callers can either lift stored binary64 inputs
losslessly or supply already-exact rational inputs for exact determinant signs,
determinant values, and solves. Exactness over binary64 input starts at the
stored values and cannot recover information rounded away before construction.
See the
[mathematical basis](https://github.com/acgetchell/la-stack/blob/v0.4.5/docs/mathematical_basis.md)
for the algorithms, validity boundaries, and supporting references.

## ✨ Design goals

- ✅ `const fn` where possible (compile-time evaluation of determinants, dot products, etc.)
- ✅ Const-generic storage (no dynamically sized matrix or vector representation)
- ✅ `Copy` types where possible
- ✅ Defined binary64 arithmetic semantics: Rust's `f64::algebraic_*`
  operations are forbidden because their unspecified reassociation, precision,
  and special-value behavior is incompatible with the crate's error bounds,
  non-finite classification, exact fallbacks, and reproducibility contract;
  deliberate `f64::mul_add` remains allowed for its defined single-rounding
  semantics
- ✅ Error-bounded f64 dot, affine-difference, and determinant filtering plus
  optional exact signs (`dot_with_errbound`, `dot_difference_with_errbound`,
  `det_errbound`, `det_sign_exact`)
- ✅ Outward-rounded interval expressions and division-free determinant signs
  through D=7, with explicit inconclusive evidence
- ✅ Exact determinant values and linear solves via optional arbitrary-precision
  arithmetic (`det_exact`, `solve_exact`, strict/rounded f64 conversions)
- ✅ Explicit algorithms (LU, solve, determinant)
- ✅ Inline, stack-backed storage for core types; optional arbitrary-precision
  exact values allocate as required
- ✅ No runtime dependencies by default (optional features may add deps)
- ✅ `unsafe` forbidden

See [CHANGELOG.md](https://github.com/acgetchell/la-stack/blob/v0.4.5/CHANGELOG.md)
for release history and
[docs/roadmap.md](https://github.com/acgetchell/la-stack/blob/v0.4.5/docs/roadmap.md)
for current release planning.

## 🚫 Anti-goals

- Alternate scalar families: `la-stack` deliberately supports finite `f64` and
  optional exact `BigRational` input domains, not `f32`, `f16`, complex, or
  generic scalar APIs
- Bare-metal performance: use [`blas`](https://crates.io/crates/blas) or
  [`lapack`](https://crates.io/crates/lapack) with a native backend selected
  through [`blas-src`](https://crates.io/crates/blas-src),
  [`lapack-src`](https://crates.io/crates/lapack-src), or
  [`openblas-src`](https://crates.io/crates/openblas-src)
- Broad general-purpose linear algebra: use [`nalgebra`](https://crates.io/crates/nalgebra)
- Large matrices/dimensions with parallelism: use [`faer`](https://crates.io/crates/faer)

## ✅ Use this crate when

- Your matrices and vectors have small, fixed dimensions known at compile time
- Stack allocation and `Copy` value semantics fit your data flow
- You want explicit LU / LDLT / determinant APIs rather than a broad algebra toolkit
- You need exact determinants, exact determinant signs, or exact linear solves
  for fixed-size systems
- You need a cheap, sound interval filter for determinant expressions assembled
  from rounded binary64 operations
- You need a certified sign or threshold comparison for a fixed-vector dot
  product or `axis · (left - right)` expression
- Robust predicates matter for geometry-style workloads near degeneracy
- You prefer a default build with no runtime dependencies

## 🔢 Scalar and bounded-value types

The public point-value scalar model deliberately has two input domains:

- finite `f64` through `Matrix<D>` and `Vector<D>` for floating-point work;
- arbitrary-precision `BigRational` through `RationalMatrix<D>` and
  `RationalVector<D>` behind the optional `"exact"` feature.

`Interval` is a separate bounded-value layer over finite `f64` endpoints. It
encloses exact-real values during a small set of outward-rounded operations and
feeds `IntervalMatrix<D>` determinant proofs; it does not make `Matrix` generic
over alternate scalars or provide a general interval package.

This is not a generic scalar-parameterized API. Exact support intentionally
covers the robustness-sensitive operations that require it: determinant sign,
determinant value, and linear solve, followed by explicit strict or rounded
conversion when an `f64` result is required. It does not promise a
`BigRational` counterpart for every floating-point helper or factorization.

Lower-precision `f32` / `f16` throughput-oriented workloads are outside the
crate's scope; they usually indicate large-matrix or accelerator-oriented use
cases better served by broader linear-algebra libraries.

## 🚀 Quickstart

The minimum supported Rust version (MSRV) is 1.98.0.

Add this to your `Cargo.toml`:

```toml
[dependencies]
la-stack = "0.4.5"
```

### Feature flags

- `default`: no runtime dependencies; includes outward-rounded `Interval` and
  `IntervalMatrix` APIs
- `exact`: exact determinant signs, determinant values, and solves over stored
  `f64` values or caller-supplied `BigRational` inputs
- `bench`: repository-development gate used only by benchmark targets and
  benchmark-input tests; application crates should not enable it

### LU solve

Solve a 5×5 system via LU:

```rust
use la_stack::prelude::*;

fn main() -> Result<(), LaError> {
    // This system requires pivoting (a[0][0] = 0), so it's a good LU demo.
    // A = J - I: zeros on diagonal, ones elsewhere.
    let a = Matrix::<5>::try_from_rows([
        [0.0, 1.0, 1.0, 1.0, 1.0],
        [1.0, 0.0, 1.0, 1.0, 1.0],
        [1.0, 1.0, 0.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0, 1.0, 0.0],
    ])?;

    let b = Vector::<5>::try_new([14.0, 13.0, 12.0, 11.0, 10.0])?;

    let lu = a.lu(DEFAULT_SINGULAR_TOL)?;
    let x = lu.solve(b)?.into_array();

    // Floating-point rounding is expected; compare with a tolerance.
    let expected = [1.0, 2.0, 3.0, 4.0, 5.0];
    for (x_i, e_i) in x.iter().zip(expected.iter()) {
        assert!((*x_i - *e_i).abs() <= 1e-12);
    }

    Ok(())
}
```

### LDLT determinant

Compute a determinant for a symmetric positive-definite matrix via LDLT (no
pivoting).

For these matrices, `LDLᵀ` is a square-root-free Cholesky form. Multiplying each
column of `L` by the square root of the corresponding diagonal entry yields a
Cholesky factor:

```rust
use la_stack::prelude::*;

fn main() -> Result<(), LaError> {
    // This matrix is symmetric positive-definite (A = L*L^T) so LDLT works without pivoting.
    let a = Matrix::<5>::try_from_rows([
        [1.0, 1.0, 0.0, 0.0, 0.0],
        [1.0, 2.0, 1.0, 0.0, 0.0],
        [0.0, 1.0, 2.0, 1.0, 0.0],
        [0.0, 0.0, 1.0, 2.0, 1.0],
        [0.0, 0.0, 0.0, 1.0, 2.0],
    ])?;

    let ldlt = match a.ldlt(DEFAULT_SINGULAR_TOL) {
        Ok(ldlt) => ldlt,
        Err(err @ LaError::Asymmetric {
            row,
            col,
            upper,
            lower,
            allowed_abs_diff,
            ..
        }) => {
            eprintln!(
                "LDLT mismatch at ({row}, {col}): {upper} vs {lower} (allowed {allowed_abs_diff})"
            );
            return Err(err);
        }
        Err(err) => return Err(err),
    };

    let det = ldlt.det()?;
    assert!((det - 1.0).abs() <= 1e-12);

    Ok(())
}
```

> ⚠️ **LDLT invariant:** The input matrix must be **exactly symmetric**: every
> mirrored pair must compare equal (`+0.0 == -0.0` is accepted). Asymmetric
> inputs passed to
> [`Matrix::ldlt`](https://docs.rs/la-stack/latest/la_stack/struct.Matrix.html#method.ldlt)
> return a typed `LaError::Asymmetric` containing both observed values and the
> required allowed difference of zero. The tolerance-based
> [`Matrix::first_asymmetry`](https://docs.rs/la-stack/latest/la_stack/struct.Matrix.html#method.first_asymmetry)
> and `Matrix::is_symmetric` methods remain useful diagnostics, but do not prove
> the exact precondition required by LDLT. Use `lu()` when exact symmetry or
> positive definiteness is not guaranteed. A negative LDLT diagonal or a zero
> diagonal with nonzero remaining coupling returns
> `LaError::NotPositiveSemidefinite` with a typed
> `PositiveSemidefiniteViolation`. An uncoupled zero or positive pivot
> at or below the caller's tolerance returns `LaError::Singular` with a
> numerical `SingularityReason`. Because these pivots are computed in binary64,
> success is not an exact positive-definiteness certificate for the stored
> matrix.

## ⚡ Compile-time determinants (D ≤ 4)

`det_direct()` is a `const fn` providing closed-form determinants for D=0–4,
using fused multiply-add where applicable. It returns `Ok(Some(det))` for those
dimensions and `Ok(None)` for D ≥ 5. `Matrix::<0>::zero().det_direct()` returns
`Ok(Some(1.0))` (the empty-product convention). For D=1–4, direct formulas
bypass LU factorization entirely. This enables compile-time evaluation when
inputs are known:

```rust
use la_stack::prelude::*;

// Evaluated entirely at compile time — no runtime cost.
const DET: Result<Option<f64>, LaError> = match Matrix::<4>::try_from_rows([
    [2.0, 0.0, 0.0, 0.0],
    [0.0, 3.0, 0.0, 0.0],
    [0.0, 0.0, 5.0, 0.0],
    [0.0, 0.0, 0.0, 7.0],
]) {
    Ok(matrix) => matrix.det_direct(),
    Err(err) => Err(err),
};

fn main() -> Result<(), LaError> {
    assert_eq!(DET?, Some(210.0));
    Ok(())
}
```

The public `det()` method automatically dispatches through the closed-form path
for D ≤ 4 and falls back to zero-tolerance LU for D ≥ 5. Tiny nonzero
determinants are not flattened by a configured pivot tolerance. The LU fallback
returns `LaError::Singular` when floating-point elimination cannot produce a
non-zero pivot; it does not misreport that numerical failure as an exact zero.
Use `lu()` directly when you need a different tolerance policy, and use the
exact determinant APIs when exact singularity classification matters.

## 📦 Outward-rounded interval determinants

`Interval` encloses expression construction that has not yet been reduced to a
single stored `f64`. Point intervals preserve finite binary64 values exactly;
`try_from_subtraction`, `try_add`, `try_mul`, `negate`, and `try_square` enclose
the corresponding exact-real operations. `IntervalMatrix<D>::det_sign()` then
uses a division-free subset expansion through D=7, returning positive,
negative, zero, or inconclusive evidence.

```rust
use la_stack::prelude::*;

fn main() -> Result<(), LaError> {
    // Relative coordinates and the lifted norm retain their construction error.
    let x = Interval::try_from_subtraction(0.1, 0.0)?;
    let y = Interval::try_from_subtraction(0.1, 0.0)?;
    let z = Interval::try_from_subtraction(0.1, 0.0)?;
    let lifted = x
        .try_square()?
        .try_add(&y.try_square()?)?
        .try_add(&z.try_square()?)?;

    let matrix = IntervalMatrix::<4>::from_rows([
        [Interval::ONE, Interval::ZERO, Interval::ZERO, Interval::ONE],
        [Interval::ZERO, Interval::ONE, Interval::ZERO, Interval::ONE],
        [Interval::ZERO, Interval::ZERO, Interval::ONE, Interval::ONE],
        [x, y, z, lifted],
    ]);
    assert_eq!(
        matrix.det_sign()?,
        IntervalDeterminantSign::Negative,
    );
    Ok(())
}
```

Every successful interval keeps finite ordered endpoints. Subnormal bounds are
preserved, both signed zeros are treated as real zero and canonicalized to
`+0.0`, and underflowed nonzero products widen toward the least subnormal value.
If an exact result range cannot fit between finite binary64 endpoints, the
operation returns `LaError::IntervalRangeExhausted` with its interval operation
recorded in `ArithmeticOperation`.

`Positive`, `Negative`, and `Zero` are proofs. `Inconclusive` only means that
the determinant enclosure overlaps zero; it must not be converted to equality
or singularity. A filtered-exact caller should rebuild the same derived
expression with `RationalMatrix` and call `det_sign()` when the interval result
is inconclusive or reports range failure. Lifting a finished `Matrix` with
`IntervalMatrix::from_matrix` encloses its stored entries, but cannot recover
rounding that occurred while those entries were assembled.

## 🔬 Exact arithmetic (`"exact"` feature)

The default build has **zero runtime dependencies**. Enable the optional
`exact` Cargo feature to add exact arithmetic methods using arbitrary-precision
rationals (this pulls in `num-bigint`, `num-rational`, and `num-traits` for
`BigRational`):

```toml
[dependencies]
la-stack = { version = "0.4.5", features = ["exact"] }
```

The feature exposes two deliberate input domains:

- `Matrix<D>` / `Vector<D>` store finite binary64 inputs. Their exact methods
  treat each stored bit pattern as its exact rational value, so the determinant
  or solve stage introduces no further roundoff. They cannot recover information
  already lost before construction.
- `RationalMatrix<D>` / `RationalVector<D>` accept coefficients already
  assembled as `BigRational`. They preserve derived differences, squared norms,
  affine coefficients, and other rational expressions without an intermediate
  `f64` conversion.

**Determinants:**

- **`det_exact()`** — returns the exact determinant as a `BigRational`
- **`det_exact_f64()`** — returns the exact determinant as `f64` only when
  it is exactly representable (or `LaError::Unrepresentable` otherwise)
- **`det_exact_rounded_f64()`** — returns the exact determinant rounded to a
  finite `f64` using IEEE 754 round-to-nearest, ties-to-even
- **`det_sign_exact()`** — infallibly returns the provably correct
  `DeterminantSign` variant (`Negative`, `Zero`, or `Positive`)

**Linear system solve:**

- **`solve_exact(b)`** — solves `Ax = b` exactly, returning a
  `RationalVector<D>`
- **`solve_exact_f64(b)`** — solves `Ax = b` exactly, returning `Vector<D>` only when
  every component is exactly representable as `f64`
- **`solve_exact_rounded_f64(b)`** — solves `Ax = b` exactly, returning each
  component rounded to finite `f64` using IEEE 754 round-to-nearest,
  ties-to-even
- **`ExactF64Conversion`** — converts an existing exact determinant or solution
  under the strict or rounded contract without repeating exact elimination

**Already-exact rational input:**

- **`RationalMatrix::det_sign()`** — returns the exact sign without constructing
  a rational determinant
- **`RationalMatrix::det()`** — returns the exact `BigRational` determinant
- **`RationalMatrix::solve(&rhs)`** — returns a `RationalVector<D>` exact
  solution
- **`try_with_rational_matrix!`** — dispatches a runtime-selected dimension
  through D=8 to a const-generic rational matrix on stable Rust

The `Matrix::det_exact*` value and conversion methods return
`LaError::DeterminantScaleOverflow` if their aggregate power-of-two scaling
exceeds the internal exponent representation. `RationalMatrix::det()` is
infallible because it clears rational row denominators without an exponent-scale
conversion. The exact solve methods for both input domains return
`LaError::Singular` with `SingularityReason::Exact` when the stored matrix is
exactly singular.

For exact-to-f64 output, strict conversions use
`UnrepresentableReason::RequiresRounding` when explicit rounding can produce a
finite value and `UnrepresentableReason::NotFinite` otherwise. Rounded
conversions opt into nearest-even rounding but still report `NotFinite` when no
finite `f64` exists.

The following 5×5 system has exact determinant 2^-60. Its exact rational inputs
therefore produce a unique solution through the general Bareiss path. Supplying
the same coefficients as `f64` inputs loses the `2^-60` perturbation at `1.0`,
making the leading rows identical and the binary64 system singular.

```rust,ignore
use core::assert_matches;

use la_stack::prelude::*;

fn main() -> Result<(), LaError> {
    // This is far below one binary64 ULP at 1.0, so 1.0 + 2^-60 rounds to 1.0.
    let epsilon = BigRational::new(1.into(), (1_u64 << 60).into());
    let one = BigRational::from_integer(1.into());
    let zero = BigRational::from_integer(0.into());

    // The leading block is [[1, 1], [1, 1 + 2^-60]]. The remaining diagonal
    // extends the example to D=5, where the general Bareiss path is used.
    let matrix = RationalMatrix::<5>::try_from_fn(|row, col| match (row, col) {
        (0, 0 | 1) | (1, 0) => one.clone(),
        (1, 1) => &one + &epsilon,
        _ if row == col => one.clone(),
        _ => zero.clone(),
    })?;
    assert_eq!(matrix.det_sign(), DeterminantSign::Positive);
    assert_eq!(matrix.det(), epsilon);

    let rhs = RationalVector::try_new([
        zero,
        -&epsilon,
        BigRational::from_integer(2.into()),
        BigRational::from_integer(3.into()),
        BigRational::from_integer(4.into()),
    ])?;
    let exact_solution = matrix.solve(&rhs)?;
    assert_eq!(
        exact_solution.as_array(),
        &[
            BigRational::from_integer(1.into()),
            BigRational::from_integer((-1).into()),
            BigRational::from_integer(2.into()),
            BigRational::from_integer(3.into()),
            BigRational::from_integer(4.into()),
        ]
    );

    // Supplying the same coefficients as f64 inputs destroys the perturbation
    // and makes the matrix singular, even though the exact solution is integral.
    let epsilon_f64 = epsilon.try_to_f64()?;
    assert_eq!((1.0 + epsilon_f64).to_bits(), 1.0_f64.to_bits());
    let f64_matrix = Matrix::<5>::try_from_rows([
        [1.0, 1.0, 0.0, 0.0, 0.0],
        [1.0, 1.0 + epsilon_f64, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 1.0],
    ])?;
    let f64_rhs = Vector::<5>::try_new([0.0, -epsilon_f64, 2.0, 3.0, 4.0])?;
    let f64_solve = f64_matrix
        .lu(DEFAULT_SINGULAR_TOL)
        .and_then(|lu| lu.solve(f64_rhs));
    assert_matches!(
        f64_solve,
        Err(LaError::Singular { .. })
    );
    Ok(())
}
```

```rust,ignore
use la_stack::prelude::*;

fn main() -> Result<(), LaError> {
    // Exact determinant
    let m = Matrix::<3>::try_from_rows([
        [1.0, 2.0, 3.0],
        [4.0, 5.0, 6.0],
        [7.0, 8.0, 9.0],
    ])?;
    assert_eq!(m.det_sign_exact(), DeterminantSign::Zero); // exactly singular

    let det = m.det_exact()?;
    assert_eq!(det, BigRational::from_integer(0.into())); // exact zero
    let det_f64 = det.try_to_f64()?;
    assert_eq!(det_f64, 0.0);

    // If strict exact-to-f64 conversion would require rounding, opt in
    // explicitly with the rounded API.
    let inexact = Matrix::<2>::try_from_rows([
        [1.0 + f64::EPSILON, 0.0],
        [0.0, 1.0 - f64::EPSILON],
    ])?;
    let exact_det = inexact.det_exact()?;
    let rounded_det = match exact_det.try_to_f64() {
        Ok(det) => det,
        Err(err) if err.requires_rounding() => exact_det.to_rounded_f64()?,
        Err(err) => return Err(err),
    };
    assert_eq!(rounded_det.to_bits(), 1.0f64.to_bits());

    // If the exact determinant cannot fit in f64, keep the BigRational value.
    let big = f64::MAX / 2.0;
    let huge = Matrix::<3>::try_from_rows([
        [0.0, 0.0, 1.0],
        [big, 0.0, 1.0],
        [0.0, big, 1.0],
    ])?;
    let huge_det = huge.det_exact()?;
    assert_eq!(
        huge_det
            .try_to_f64()
            .err()
            .and_then(|err| err.unrepresentable_reason()),
        Some(UnrepresentableReason::NotFinite)
    );
    println!("exact determinant = {huge_det}");

    // Exact linear system solve
    let a = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]])?;
    let b = Vector::<2>::try_new([5.0, 11.0])?;
    let exact_x = a.solve_exact(b)?;
    let x = exact_x.try_to_f64()?.into_array();
    assert!((x[0] - 1.0).abs() <= f64::EPSILON);
    assert!((x[1] - 2.0).abs() <= f64::EPSILON);

    Ok(())
}
```

With the `exact` feature enabled, `RationalMatrix`, `RationalVector`,
`DeterminantSign`, `ExactF64Conversion`, `BigInt`, and `BigRational` are
re-exported from the crate root and prelude,
alongside the most commonly needed `num-traits` items (`FromPrimitive`,
`ToPrimitive`, `Signed`). This lets consumers construct exact values
(`BigRational::from_f64`, `from_i64`), query sign (`is_positive` /
`is_negative`), and convert back to `f64` (`try_to_f64`, `to_rounded_f64`, or
the raw `to_f64`) with a single
`use la_stack::prelude::*;` — no need to add `num-bigint`, `num-rational`,
or `num-traits` to their own `Cargo.toml`. Use
`DeterminantSign::as_i8()` only when numeric −1/0/+1 interoperability is
required.

For `det_sign_exact()`, D ≤ 4 matrices first use a fast f64 filter
(error-bounded `det_direct()`) when its rounded intermediates stay in the normal
range or are exact structural zeros. An inconclusive filter falls back to the
same direct determinant expansion in `BigInt`. D ≥ 5 skips the closed-form
filter and uses fraction-free Bareiss elimination in `BigInt`.
Because `Matrix` stores only finite entries, arithmetic range failures in the
filter are inconclusive rather than errors and the exact fallback is total.

## 🎯 Certified dot products and affine differences

`Vector::dot_with_errbound()` evaluates the same left-to-right FMA tree as
`Vector::dot()` and returns its estimate together with a certified absolute
roundoff bound. `Vector::dot_difference_with_errbound()` directly evaluates

```text
Σᵢ axis[i] × (left[i] - right[i])
```

as two FMAs per coordinate. It does not first round `left - right` into a new
`Vector`, so the certificate covers the intended expression over the original
stored binary64 coordinates.

The opaque `ScalarWithErrorBound` exposes the estimate, absolute error bound,
and finite outward-rounded lower and upper bounds. Those endpoints support
positive, negative, and caller-selected threshold proofs:

```rust
use la_stack::prelude::*;

fn is_separated<const D: usize>(
    axis: &Vector<D>,
    left: &Vector<D>,
    right: &Vector<D>,
    threshold: f64,
) -> Result<Option<bool>, LaError> {
    let Some(value) = axis.dot_difference_with_errbound(left, right)? else {
        return Ok(None);
    };
    if value.lower_bound() > threshold {
        Ok(Some(true))
    } else if value.upper_bound() <= threshold {
        Ok(Some(false))
    } else {
        Ok(None)
    }
}

# fn main() -> Result<(), LaError> {
let axis = Vector::<2>::try_new([2.0, -1.0])?;
let left = Vector::<2>::try_new([4.0, 1.0])?;
let right = Vector::<2>::try_new([1.0, 3.0])?;
assert_eq!(is_separated(&axis, &left, &right, 1.0)?, Some(true));
# Ok(())
# }
```

An interval that overlaps the threshold is inconclusive, not equal. Likewise,
`Ok(None)` means gradual underflow or proof-only range exhaustion prevented a
certificate. A filtered-exact caller should rebuild the same dot or affine
expression in `BigRational` (available through the `exact` feature) or another
exact backend. A `LaError::NonFinite` instead reports that the specified FMA
estimate itself overflowed. These certified bounds describe roundoff in a fixed
arithmetic tree; they are distinct from user-selected numerical tolerances.

## 🛡️ Adaptive determinant filtering (D ≤ 4)

`det_direct_with_errbound()` returns a closed-form determinant together with
the conservative absolute error bound used by the fast filter, computed from
one call that evaluates the determinant once and computes its matching bound.
It returns `None` when a D ≤ 4 computation may be affected by gradual
underflow, as well as for unsupported D ≥ 5 dimensions.
It returns `LaError::NonFinite` if the determinant or bound computation
overflows to NaN or infinity.
This method does NOT require the `exact` feature — it uses pure f64 arithmetic
and is available by default. Use `det_errbound()` when only the bound is needed.
The paired API enables custom adaptive-precision logic for geometric predicates:

```rust,ignore
use la_stack::prelude::*;

fn adaptive_det_sign<const D: usize>(
    matrix: &Matrix<D>,
) -> DeterminantSign {
    if let Ok(Some(estimate)) = matrix.det_direct_with_errbound() {
        if estimate.determinant().abs() > estimate.absolute_error_bound() {
            return if estimate.determinant() > 0.0 {
                DeterminantSign::Positive
            } else {
                DeterminantSign::Negative
            };
        }
    }

    matrix.det_sign_exact()
}

fn main() -> Result<(), LaError> {
    let identity = Matrix::<3>::identity();
    assert_eq!(
        adaptive_det_sign(&identity),
        DeterminantSign::Positive
    );

    // A zero determinant cannot pass the f64 sign filter, so this exercises
    // the exact fallback.
    let singular = Matrix::<3>::try_from_rows([
        [1.0, 2.0, 3.0],
        [4.0, 5.0, 6.0],
        [7.0, 8.0, 9.0],
    ])?;
    assert_eq!(adaptive_det_sign(&singular), DeterminantSign::Zero);

    // The f64 filter overflows for this finite matrix, but the exact fallback
    // still resolves its positive determinant sign.
    let big = f64::MAX / 2.0;
    let overflowing = Matrix::<3>::try_from_rows([
        [0.0, 0.0, 1.0],
        [big, 0.0, 1.0],
        [0.0, big, 1.0],
    ])?;
    assert_eq!(
        adaptive_det_sign(&overflowing),
        DeterminantSign::Positive
    );

    Ok(())
}
```

The error coefficients (`ERR_COEFF_2`, `ERR_COEFF_3`, `ERR_COEFF_4`) are
conservative, dimension-specific constants, not caller-tunable tolerances. The
[mathematical basis](https://github.com/acgetchell/la-stack/blob/v0.4.5/docs/mathematical_basis.md#determinants-and-certified-sign-filtering)
documents the bound and states its range preconditions. The constants are explicit
crate-root exports for advanced users who want to compose the same bound:
`use la_stack::{ERR_COEFF_2, ERR_COEFF_3, ERR_COEFF_4};`. They intentionally stay
out of the common prelude.

## 🧩 API at a glance

| Type | Storage | Purpose | Key methods |
|---|---|---|---|
| `Vector<D>` | `[f64; D]` | Finite fixed-length vector for input and computation | `try_new`, `as_array`, `into_array`, `dot`, `dot_with_errbound`, `dot_difference_with_errbound`, `norm2_sq` |
| `Matrix<D>` | `[[f64; D]; D]` | Finite square matrix for input and computation | See below |
| `Interval` | Two finite ordered `f64` bounds | Outward-rounded exact-real enclosure | `try_new`, `point`, `try_from_subtraction`, `try_add`, `try_mul`, `negate`, `try_square` |
| `IntervalMatrix<D>` | `[[Interval; D]; D]` | Division-free determinant enclosure and sign proof through D=7 | `from_rows`, `try_from_point_rows`, `from_matrix`, `det`, `det_sign` |
| `IntervalDeterminantSign` | enum | Positive, negative, zero, or inconclusive determinant evidence | — |
| `RationalVector<D>`¹ | `[BigRational; D]` | Exact rational right-hand side and solution | `try_new`, `try_from_fn`, `as_array`, `into_array`, `get` |
| `RationalMatrix<D>`¹ | `[[BigRational; D]; D]` | Exact rational matrix for determinant and solve operations | `try_from_rows`, `try_from_fn`, `as_rows`, `det_sign`, `det`, `solve` |
| `DeterminantWithErrorBound` | Opaque validated pair | Paired direct determinant and certified absolute bound | `determinant`, `absolute_error_bound` |
| `ScalarWithErrorBound` | Opaque validated certificate | Paired scalar estimate, absolute bound, and outward endpoints | `estimate`, `absolute_error_bound`, `lower_bound`, `upper_bound` |
| `Lu<D>` | Inline factors + permutation | Factorization for solves/det | `solve`, `det` |
| `Ldlt<D>` | Inline factors | No-pivot SPD factorization for solves/det | `solve`, `det` |
| `Tolerance` | finite non-negative `f64` | Validated numerical threshold | `try_new`, `get` |
| `LaError` | typed variants and reasons | Structured, actionable failure reporting | See error semantics below |
| `DeterminantSign`¹ | enum | Exact determinant sign | `as_i8` |
| `ExactF64Conversion`¹ | trait | Strict or explicitly rounded conversion of exact results to `f64` | `try_to_f64`, `to_rounded_f64` |

`Matrix<D>` and `Vector<D>` use the intentional inline `f64` scalar model.
`IntervalMatrix<D>` retains inline fixed-size storage and uses a fixed 128-entry
stack workspace for its supported determinant dimensions. The exact-feature
rational types retain fixed-size outer arrays while their `BigRational` scalars
use arbitrary-precision integer storage.

For a runtime-selected interval dimension from 0 through
`MAX_INTERVAL_MATRIX_DIM` (7), `try_with_interval_matrix!` dispatches to a
concrete `IntervalMatrix<N>`. This is useful when stable Rust cannot express a
derived const dimension such as `D + 1`.

For a runtime dimension from 0 through `MAX_STACK_MATRIX_DISPATCH_DIM` (7),
`try_with_stack_matrix!` dispatches to a concrete `Matrix<N>` while preserving
inline stack storage. Larger dimensions produce
`LaError::UnsupportedDimension`, converted through `From<LaError>` into the
closure's declared `Result` error type; the macro does not introduce a
dynamically sized matrix representation.

With the exact feature, a runtime dimension from 0 through
`MAX_RATIONAL_MATRIX_DISPATCH_DIM` (8) can similarly be dispatched with
`try_with_rational_matrix!` to a concrete `RationalMatrix<N>`. Larger
dimensions produce `LaError::UnsupportedDimension`; the rational macro also
preserves the const-generic representation rather than introducing a
dynamically sized matrix type.

`Matrix<D>` key methods: `as_rows`, `into_rows`, `lu`, `ldlt`, `det`,
`det_direct`, `det_direct_with_errbound`, `det_errbound`,
`det_exact`¹, `det_exact_f64`¹, `det_exact_rounded_f64`¹, `det_sign_exact`¹,
`solve_exact`¹, `solve_exact_f64`¹, `solve_exact_rounded_f64`¹.
Matrix and vector constructors validate non-finite inputs at public API
boundaries. After construction, `Matrix<D>` and `Vector<D>` carry that
finite-storage invariant directly, so factorization kernels do not repeat an
O(D²) input scan. Computed factor matrices are still checked before they become
observable results.

`Matrix::as_rows` and `Vector::as_array` borrow their validated backing arrays;
`Matrix::into_rows` and `Vector::into_array` consume the value and return the
owned fixed-size arrays.

`Matrix::get(row, col)` returns `None` for out-of-bounds coordinates;
`Matrix::try_get` instead returns a structured `LaError` preserving those
coordinates. The single fallible `Matrix::set` validates both coordinates and
finiteness before mutating the matrix.

`LaError` and its reason/location enums are non-exhaustive. Numerical
singularity records the `FactorizationKind`,
observed pivot magnitude, and tolerance, while exact-arithmetic singularity is
identified separately. `LaError::NonFinite` retains the crate-wide non-finite
contract but uses `NonFiniteOrigin`, `NonFiniteLocation`, and
`ArithmeticOperation` to distinguish invalid inputs from computed overflow.
`LaError::InvertedInterval` preserves rejected finite bounds when the lower
endpoint exceeds the upper endpoint, and `LaError::IntervalRangeExhausted`
distinguishes finite-input interval range loss from a non-finite value.
`InvalidToleranceReason` distinguishes negative from non-finite tolerances, and
`PositiveSemidefiniteViolation` distinguishes negative LDLT pivots from a zero
pivot with nonzero coupling. Match these public enums with a wildcard and use
`..` for struct-style variants so future error context can be added without
breaking callers.

¹ Requires `features = ["exact"]`.

## 📊 Benchmarks (vs nalgebra/faer)

![LU solve (factor + solve): median time vs dimension][lu-solve-benchmark]

Raw data:
[docs/assets/bench/vs_linalg_lu_solve_median.csv](https://github.com/acgetchell/la-stack/blob/v0.4.5/docs/assets/bench/vs_linalg_lu_solve_median.csv)
Measurement provenance:
[docs/assets/bench/vs_linalg_lu_solve_median.provenance.json][benchmark-provenance]

Representative benchmark: `lu_solve` factors the matrix and solves one
right-hand side. Median time is lower-is-better, and the “la-stack vs
nalgebra/faer” columns show the % time reduction relative to each baseline
(positive means the recorded la-stack median is lower). These are descriptive
point-estimate ratios, not statistical significance claims or an aggregate score
across operations.

Timings count only when the implementation preserves the documented
correctness guarantees and invariants. Performance claims require comparable
before-and-after evidence using the same inputs, configuration, and environment.
This snapshot records the measured source state, available CPU model, operating system, Rust
toolchain, dependency lock and harness digests, Criterion command, and
correctness-gate result in the adjacent JSON sidecar. The publication workflow
requires complete canonical-dimension coverage and regenerates the CSV, SVG,
README table, and provenance together.

For the full per-kernel comparison methodology, algorithm citations, input
construction, and release-comparison workflow details, see
[docs/BENCHMARKING.md](https://github.com/acgetchell/la-stack/blob/v0.4.5/docs/BENCHMARKING.md).
For the current release-to-release performance snapshot, see
[docs/PERFORMANCE.md](https://github.com/acgetchell/la-stack/blob/v0.4.5/docs/PERFORMANCE.md).
The exact release suite includes the already-exact rational-input groups for
D=2 through D=8. Those rows report `RationalMatrix::det_sign`, `det`, and
`solve` alongside straightforward `BigRational` Gaussian determinant and solve
references. Releases produced with the rational-input harness include Criterion
point estimates and confidence intervals for these rows; comparisons against a
pre-API baseline retain them as explicit current-only measurements.

The focused `interval` Criterion suite covers conclusive and inconclusive
relative-coordinate lifted determinant signs at D=4 and the maximum supported
D=7 workload. Run it with `just bench-interval`; fixture validation stays
outside the timed closures.

The focused `linear_form` Criterion suite compares plain and certified dot
products and covers both well-separated and inconclusive dot/affine-difference
filters at D=4. Run it with `just bench-linear-form`; exact small-integer fixture
expectations are validated outside the timed closures.

<!-- BENCH_TABLE:lu_solve:median:new:BEGIN -->

| D | la-stack median (ns) | nalgebra median (ns) | faer median (ns) | reduction vs nalgebra (point est.) | reduction vs faer (point est.) |
|---:|--------------------:|--------------------:|----------------:|---------------------:|----------------:|
| 2 | 2.044 | 4.601 | 151.939 | +55.6% | +98.7% |
| 3 | 9.989 | 23.513 | 196.357 | +57.5% | +94.9% |
| 4 | 21.865 | 54.716 | 223.910 | +60.0% | +90.2% |
| 5 | 44.510 | 71.219 | 293.420 | +37.5% | +84.8% |
| 8 | 145.405 | 188.352 | 381.872 | +22.8% | +61.9% |
| 16 | 672.491 | 585.261 | 897.236 | -14.9% | +25.0% |
| 32 | 2,777.707 | 2,501.361 | 2,952.778 | -11.0% | +5.9% |
| 64 | 17,357.785 | 13,878.401 | 12,199.761 | -25.1% | -42.3% |

<!-- BENCH_TABLE:lu_solve:median:new:END -->

## 📋 Examples

The `examples/` directory contains small, runnable programs:

- **`solve_5x5`** — solve a 5×5 system via LU with partial pivoting
- **`det_5x5`** — determinant of a 5×5 matrix via LU
- **`ldlt_solve_3x3`** — solve a 3×3 symmetric positive definite system via LDLT
- **`const_det_4x4`** — compile-time 4×4 determinant via `det_direct()`
- **`exact_det_3x3`** — exact determinant value of a near-singular 3×3 matrix (requires `exact` feature)
- **`exact_sign_3x3`** — exact determinant sign of a near-singular 3×3 matrix (requires `exact` feature)
- **`exact_solve_3x3`** — exact solve of a near-singular 3×3 system vs f64 LU (requires `exact` feature)
- **`rational_input_5x5`** — exact rational solve of a 5×5 system that becomes singular as f64 (requires `exact` feature)

```bash
just examples
# or individually:
cargo run --example solve_5x5
cargo run --example det_5x5
cargo run --example ldlt_solve_3x3
cargo run --example const_det_4x4
cargo run --features exact --example exact_det_3x3
cargo run --features exact --example exact_sign_3x3
cargo run --features exact --example exact_solve_3x3
cargo run --features exact --example rational_input_5x5
```

## 🤝 Contributing

A short contributor workflow:

Install Rust 1.98.0 through [rustup](https://rustup.rs/), Git,
[GitHub CLI](https://cli.github.com/), Python 3.14,
[`uv` 0.12.5](https://docs.astral.sh/uv/), and `jq`. Then install the pinned
`just` release from its locked dependency graph:

```bash
cargo install --locked just --version 1.58.0
just setup        # install/verify dev tools + sync Python deps + build
just check        # lint/validate (non-mutating)
just fix          # apply auto-fixes (mutating)
just ci           # lint + tests + examples + bench compile
```

The repository uses `cargo-nextest` for runnable Rust tests, `cargo-machete`
for unused-dependency checks, `rumdl` for Markdown, `dprint` plus `yamllint`
for YAML/CFF, `taplo` for TOML, and `typos` for spelling. Python 3.14 support
tooling is locked with `uv` and checked by Ruff, Ty, and Semgrep. GitHub Actions
references are SHA-pinned, restricted to an explicit allowlist, and kept with
readable version comments for review.

CI runs `just ci` on Ubuntu, macOS, and Windows to keep platform coverage
aligned with the local comprehensive validation path.

For coverage commands and report locations, see
[`docs/COVERAGE.md`](https://github.com/acgetchell/la-stack/blob/v0.4.5/docs/COVERAGE.md).
For the full contributor workflow, see
[CONTRIBUTING.md](https://github.com/acgetchell/la-stack/blob/v0.4.5/CONTRIBUTING.md).

## 📝 Citation

If you use this library in academic work, please cite it using
[CITATION.cff](https://github.com/acgetchell/la-stack/blob/v0.4.5/CITATION.cff)
(or GitHub's "Cite this repository" feature). Tagged releases are archived on
Zenodo under the
[all-versions concept DOI](https://doi.org/10.5281/zenodo.18158926).

## 📚 References

For canonical references to the algorithms used by this crate, see
[REFERENCES.md](https://github.com/acgetchell/la-stack/blob/v0.4.5/REFERENCES.md).

## 🤖 AI Agents

AI coding assistants should read
[AGENTS.md](https://github.com/acgetchell/la-stack/blob/v0.4.5/AGENTS.md)
before proposing or applying changes. See
[CONTRIBUTING.md](https://github.com/acgetchell/la-stack/blob/v0.4.5/CONTRIBUTING.md)
for the repository's AI-assisted development note.

## 📄 License

BSD 3-Clause License. See [LICENSE](https://github.com/acgetchell/la-stack/blob/v0.4.5/LICENSE).

[audit-badge]: https://github.com/acgetchell/la-stack/actions/workflows/audit.yml/badge.svg
[audit-workflow]: https://github.com/acgetchell/la-stack/actions/workflows/audit.yml
[benchmark-provenance]: https://github.com/acgetchell/la-stack/blob/v0.4.5/docs/assets/bench/vs_linalg_lu_solve_median.provenance.json
[clippy-badge]: https://github.com/acgetchell/la-stack/actions/workflows/rust-clippy.yml/badge.svg
[clippy-workflow]: https://github.com/acgetchell/la-stack/actions/workflows/rust-clippy.yml
[lu-solve-benchmark]: https://raw.githubusercontent.com/acgetchell/la-stack/v0.4.5/docs/assets/bench/vs_linalg_lu_solve_median.svg
