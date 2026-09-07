# References and citations

## Contents

- [How to cite this library](#how-to-cite-this-library)
- [AI-Assisted Development Tools](#ai-assisted-development-tools)
- [Linear algebra algorithms](#linear-algebra-algorithms)
  - [Absolute error bound for closed-form determinants](#absolute-error-bound-for-closed-form-determinants)
  - [Certified fixed-vector reductions](#certified-fixed-vector-reductions)
  - [Exact determinant sign](#exact-determinant-sign-adaptive-precision-integer-arithmetic)
  - [Exact linear system solve](#exact-linear-system-solve-hybrid-bareiss--bigrational)
  - [Exact rational inputs](#exact-rational-inputs-and-row-denominator-clearing)
  - [Exact-to-binary64 conversion](#exact-to-binary64-conversion)
  - [f64 → integer decomposition](#f64--integer-decomposition-decompose_proven_finite_f64)
  - [Gram matrices and geometric measures](#gram-matrices-and-geometric-measures)
  - [LDLᵀ factorization](#ldlt-factorization)
  - [LU decomposition](#lu-decomposition-gaussian-elimination-with-partial-pivoting)
  - [Outward-rounded interval determinant sign](#outward-rounded-interval-determinant-sign)
  - [Scaled determinant products](#scaled-determinant-products)
  - [Scaled Euclidean vector norm](#scaled-euclidean-vector-norm)
- [References](#references)

## How to cite this library

If you use this library in your research or project, please cite it using the information in
[CITATION.cff](CITATION.cff). This file contains structured citation metadata that can be
processed by GitHub and other platforms.

Tagged releases are archived on Zenodo under the all-versions concept DOI
[10.5281/zenodo.18158926](https://doi.org/10.5281/zenodo.18158926).

## AI-Assisted Development Tools

- Anthropic. "Claude." <https://www.anthropic.com/claude>.
- CodeRabbit AI, Inc. "CodeRabbit." <https://coderabbit.ai/>.
- KiloCode. "KiloCode AI Engineering Assistant." <https://kilocode.ai/>.
- OpenAI. "ChatGPT." <https://openai.com/chatgpt>.
- OpenAI. "Codex." <https://openai.com/codex/>.
- Warp Dev, Inc. "WARP." <https://www.warp.dev/>.

All AI-generated output was reviewed and/or edited by the maintainer.
No generated content was used without human oversight.

## Linear algebra algorithms

### Absolute error bound for closed-form determinants

`Matrix::det_errbound()` returns a conservative Shewchuk-style absolute error bound \[[8]\]
for `Matrix::det_direct()` in dimensions 2–4 when every rounded intermediate is normal
or an exact structural zero. The returned bound is the rounded product
`fl(ERR_COEFF_D · p_hat)`, where `p_hat` is the computed approximation to
`p(|A|) = perm(|A|)`, the absolute Leibniz sum. The dimension-specific constants
`ERR_COEFF_2`, `ERR_COEFF_3`, and `ERR_COEFF_4` cover rounding in the determinant,
the permanent, and the final multiplication.
The [coefficient derivation](docs/mathematical_basis.md#derivation-of-the-returned-determinant-bound)
counts the longest rounding paths in both dense and sparse evaluation trees and
proves the returned binary64 bound using the arithmetic model in \[[9], [10], [11]\].
The method returns `None` when gradual underflow could violate the relative-error model.
The same bound is used internally by `det_sign_exact()`'s fast filter, but
`det_errbound()` itself is available without the `exact` feature, so downstream crates
can build custom adaptive-precision logic with pure f64 arithmetic.

### Certified fixed-vector reductions

`Vector::dot_with_errbound()` and `Vector::dot_difference_with_errbound()` use
deterministic left-to-right binary64 FMA reductions. When their rounded
intermediates stay normal or are exact zeros, the standard
`gamma_n = n·u / (1 - n·u)` model, with `u = 2^-53`, `n = D` for the dot
product, and `n = 2D` for the affine difference, bounds the absolute forward error by
`gamma_n Σ |a_i b_i|` \[[9], [10], [11]\]. The magnitude sum and final bound are
rounded upward, while `TwoSum` supplies outward endpoints. Gradual underflow or
proof-only range exhaustion makes the filter unavailable rather than turning an
inconclusive result into equality. The affine form evaluates alternating
`axis_i × left_i` and `-axis_i × right_i` FMAs, so its certificate covers the
original coordinates rather than an already-rounded difference vector.

### Exact determinant sign (adaptive-precision integer arithmetic)

`det_sign_exact()` uses a Shewchuk-style f64 error-bound filter \[[8]\] (the same bound exposed
by `det_errbound()` above) backed by exact `BigInt` arithmetic. Each f64 entry is decomposed
into `mantissa × 2^exponent` and scaled to a common integer base. Dimensions 0–4 use direct
integer determinant expansions; D ≥ 5 uses integer-only Bareiss elimination \[[7]\]. Neither
path constructs `BigRational` values or performs GCD normalization.
See `src/exact.rs` for the full architecture description.

### Exact linear system solve (hybrid Bareiss / BigRational)

`solve_exact()`, `solve_exact_f64()`, and `solve_exact_rounded_f64()` share the determinant
path's exact f64 decomposition and integer scaling. Matrix and RHS entries are decomposed via
IEEE 754 bit extraction \[[9]\]. Matrix and RHS scales start from their respective minimum
exponents. When `|e_rhs − e_matrix| ≤ 64`, both sides use `min(e_rhs, e_matrix)` as the shared
scale; when `|e_rhs − e_matrix| > 64`, they retain independent scales so one side is not
inflated excessively.
Forward elimination runs in `BigInt` using Bareiss fraction-free updates
\[[7]\]—no `BigRational` and no GCD normalisation in the `O(D³)` phase. The upper-triangular
result is then lifted into `BigRational` for back-substitution, where fractions are inherent
and the cost is only `O(D²)`. Row swaps from first-non-zero pivoting are applied to both the
matrix and RHS. After back-substitution, multiplying by the exact power-of-two scale ratio
`2^(e_rhs − e_matrix)` recovers the solution to the original `A x = b` system.

### Exact rational inputs and row denominator clearing

`RationalMatrix` and `RationalVector` accept exact rational coefficients and
canonicalize each quotient to lowest terms with a positive denominator. For row
`i`, a positive denominator LCM `s_i` makes `A_int[i, :] = s_i A[i, :]` integral.
The LCM is accumulated with `lcm(a, b) = (a / gcd(a, b)) b`, using Euclidean GCD.
Determinant multilinearity gives `det(A_int) = (Π s_i) det(A)` \[[12]\]: the
positive scales preserve sign, and determinant values divide by their product.
Solves include the RHS denominator in each row's LCM, preserving the solution
set of the augmented system. The resulting integer matrices reuse direct
expansions through D=4 and the Bareiss determinant/solve backend \[[7]\].
See the [rational-input construction](docs/mathematical_basis.md#exact-arithmetic-over-rational-inputs).

### Exact-to-binary64 conversion

Strict conversion accepts only dyadic rationals whose reduced significand and
exponent fit binary64 exactly. Rounded conversion uses round-to-nearest,
ties-to-even, including subnormal values and signed underflow to zero \[[9], [10]\].
The integer-and-exponent determinant path reads retained, guard, and sticky
bits directly from `BigInt`; rational-value rounding delegates to
`num-rational`'s `ToPrimitive::to_f64`. Both paths reject rounding overflow.
`RequiresRounding` means a strict conversion failed but finite rounded output
is available; `NotFinite` means rounding cannot produce finite output.
See the [representability and rounding criteria](docs/mathematical_basis.md#exact-to-binary64-conversion).

### f64 → integer decomposition (`decompose_proven_finite_f64`)

Both the determinant and solve paths convert their finite-by-construction entries via
`decompose_proven_finite_f64`, which extracts the IEEE 754 binary64 sign, unbiased exponent,
and significand \[[9]\]. For nonzero `x`, it strips trailing zeros from the
significand so `|x| = m · 2^e` with `m` odd; signed zeros use a separate zero
component. The integer matrix is then assembled by shifting each mantissa left by
`exp − e_min`, giving a GCD-free exact-integer starting point. Solves and D ≥ 5 determinants
then apply Bareiss elimination; D ≤ 4 determinants use direct expansions. The test-only
fallible wrapper `decompose_f64` verifies rejection of non-finite raw scalars, while the
test-only `f64_to_big_rational` helper packages the same decomposition into a single
`BigRational`. See Goldberg \[[10]\] for background on floating-point representation and
conversion.

### Gram matrices and geometric measures

A Gram matrix collects pairwise inner products: `G[i,j] = v_i · v_j`.
Writing the vectors as rows of `V` gives `G = V Vᵀ`. Its diagonal contains
squared lengths; off-diagonal entries describe angles through
`v_i · v_j = ||v_i|| ||v_j|| cos(θ)` for nonzero vectors.

In exact real arithmetic, `G` is positive semidefinite and is positive definite
exactly when the vectors are linearly independent. For `M ≤ N`, `det(G)` is
the squared M-dimensional volume spanned by the vectors. For simplex edge
vectors from a common vertex, the simplex volume is `sqrt(det(G)) / M!` \[[16]\].
This applies to triangles embedded in 3D and to higher-dimensional facets.

`gram_matrix` computes rounded binary64 entries with exact mirrored symmetry;
it does not certify rank, positive definiteness, or volume accuracy. See \[[9], [10], [11], [12]\]
for floating-point and conditioning background.

<a name="ldlt-factorization"></a>

### LDLᵀ factorization (exactly symmetric positive-definite inputs)

The no-pivot LDLT implementation targets `A = L D Lᵀ` for exactly symmetric
positive-definite inputs \[[4], [5], [11], [12]\]. Successful construction requires every
computed diagonal pivot to be positive and greater than the caller's tolerance.
Computed zero and tolerance-small positive pivots are therefore part of the
typed diagnostic domain, not returned in a usable factorization. Because the
pivots are computed in binary64, a successful factorization is not an exact
certificate that the represented matrix is positive definite.

For pivoted variants used for symmetric *indefinite* matrices, see \[[6], [11], [12]\].

### LU decomposition (Gaussian elimination with partial pivoting)

The LU implementation targets `P A = L U` and uses partial pivoting: each step
selects the remaining entry of largest magnitude in the active column. Partial
pivoting is a practical stability strategy, not an unconditional accuracy
guarantee; worst-case growth and average-case behavior are distinct concerns.

See \[[1], [2], [3], [11], [12]\] for stability analysis, finite-precision behavior, and standard
algorithmic background.

### Outward-rounded interval determinant sign

`Interval` uses IEEE-754 round-to-nearest binary64 operations plus adjacent
representable values to enclose exact-real addition, subtraction,
multiplication, and square results \[[9], [10], [11]\]. Addition and subtraction
use an error-free `TwoSum` residual \[[8]\]; multiplication independently compares
the exact integer-significand product with the rounded binary64 result,
including gradual underflow to zero. Results whose exact range cannot fit
between finite binary64 endpoints return a typed range failure rather than
storing infinity. For the broader standardized interval arithmetic model, see
\[[14]\]; this crate does not claim IEEE 1788 conformance.

`IntervalMatrix::det()` evaluates the Leibniz expansion with a division-free
column-subset dynamic program. It uses `2^D` inline interval states and
`D × 2^(D-1)` coefficient products through D=7. A determinant interval strictly
separated from zero certifies its sign; `[0, 0]` certifies zero; every other
overlap is explicitly inconclusive. The determinant identity is standard
linear algebra \[[12]\]; the interval evaluation and subset-DP
organization are implemented specifically for this crate's small
fixed-dimension scope.

### Scaled determinant products

`Lu::det` and `Ldlt::det` first multiply factor diagonals directly. A zero,
subnormal, or non-finite running product triggers a complete replay with
normalized mantissas and a separate power-of-two exponent. Each nonzero factor
has `1 ≤ mantissa < 2`; normalization therefore keeps intermediate mantissa
products away from overflow and underflow \[[9], [10]\].

The last factor is deferred so a subnormal final result needs only one final
rounding in the destination range. Earlier mantissa products still round:
scaling does not make the product exact or supply a certified error bound.
This crate-specific range-management implementation also retains LU permutation
parity. See the [scaled-product description](docs/mathematical_basis.md#scaled-determinant-products).

### Scaled Euclidean vector norm

`Vector::norm` maintains a scale and a sum of squares relative to that scale,
avoiding raw coordinate squares that would overflow or underflow. This follows
the scaled safe-norm approach described by Blue \[[15]\]. The implementation retains
a deterministic coordinate order and documents its binary64 rounding contract;
it does not generally claim a certified error bound or correct rounding. A
fixed-size integer sum of exact binary64 squares handles the upper range, using
the representation and nearest-even rounding model in \[[9], [10]\] to distinguish
finite results from overflow.

## References

Reference numbers are stable citation keys used throughout the code and documentation.
Each entry has a permanent anchor such as `#ref-8` for direct links.
The bibliography retains its thematic order; algorithm headings above are
alphabetized for navigation without renumbering citations.

1. <a name="ref-1"></a> Trefethen, Lloyd N., and Robert S. Schreiber. "Average-case stability of Gaussian elimination."
   *SIAM Journal on Matrix Analysis and Applications* 11.3 (1990): 335–360.
   [DOI](https://doi.org/10.1137/0611023) ·
   [PDF](https://people.maths.ox.ac.uk/trefethen/publication/PDF/1990_44.pdf)
2. <a name="ref-2"></a> Businger, P. A. "Monitoring the Numerical Stability of Gaussian Elimination."
   *Numerische Mathematik* 16.4 (1971): 360–361.
   [DOI](https://doi.org/10.1007/BF02165006) · [Full text](https://eudml.org/doc/132040)
3. <a name="ref-3"></a> Huang, Han, and K. Tikhomirov. "Average-case analysis of the Gaussian elimination with partial pivoting."
   *Probability Theory and Related Fields* 189 (2024): 501–567.
   [DOI](https://doi.org/10.1007/s00440-024-01276-2) ·
   [Open-access article](https://link.springer.com/article/10.1007/s00440-024-01276-2) ·
   [arXiv:2206.01726](https://arxiv.org/abs/2206.01726)
4. <a name="ref-4"></a> Cholesky, André-Louis. "Sur la résolution numérique des systèmes d'équations linéaires."
   *Bulletin de la Sabix* 39 (2005): 81–95. Manuscript dated 2 December 1910.
   [DOI](https://doi.org/10.4000/sabix.529)
5. <a name="ref-5"></a> Brezinski, Claude. "La méthode de Cholesky."
   *Revue d'histoire des mathématiques* 11.2 (2005): 205–238.
   [DOI](https://doi.org/10.24033/rhm.30) ·
   [Full text](https://www.numdam.org/articles/10.24033/rhm.30/)
6. <a name="ref-6"></a> Bunch, James R., Linda Kaufman, and Beresford N. Parlett. "Decomposition of a Symmetric Matrix."
   *Numerische Mathematik* 27 (1976): 95–109.
   [DOI](https://doi.org/10.1007/BF01399088) · [Full text](https://eudml.org/doc/132435)
7. <a name="ref-7"></a> Bareiss, Erwin H. "Sylvester's Identity and Multistep Integer-Preserving Gaussian
   Elimination." *Mathematics of Computation* 22.103 (1968): 565–578.
   [DOI](https://doi.org/10.1090/S0025-5718-1968-0226829-0) ·
   [PDF](https://www.ams.org/journals/mcom/1968-22-103/S0025-5718-1968-0226829-0/S0025-5718-1968-0226829-0.pdf)
8. <a name="ref-8"></a> Shewchuk, Jonathan Richard. "Adaptive Precision Floating-Point Arithmetic and Fast
   Robust Geometric Predicates." *Discrete & Computational Geometry* 18.3 (1997): 305–363.
   [DOI](https://doi.org/10.1007/PL00009321) ·
   [PDF](https://people.eecs.berkeley.edu/~jrs/papers/robustr.pdf)
   Also: Technical Report CMU-CS-96-140, Carnegie Mellon University, May 1996.
9. <a name="ref-9"></a> IEEE Computer Society. "IEEE Standard for Floating-Point Arithmetic." *IEEE Std 754-2019*
   (Revision of IEEE 754-2008), 2019.
   [DOI](https://doi.org/10.1109/IEEESTD.2019.8766229)
   Section 3.4 (binary64 format): 1 sign bit, 11 exponent bits (bias 1023), 52 trailing
   significand bits; subnormals have biased exponent 0 with implicit leading 0.
10. <a name="ref-10"></a> Goldberg, David. "What Every Computer Scientist Should Know About Floating-Point
    Arithmetic." *ACM Computing Surveys* 23.1 (1991): 5–48.
    [DOI](https://doi.org/10.1145/103162.103163) ·
    [Authorized HTML reprint](https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html)
    Comprehensive survey of floating-point representation, rounding, and conversion.
11. <a name="ref-11"></a> Higham, Nicholas J. *Accuracy and Stability of Numerical Algorithms*. 2nd ed.
    Society for Industrial and Applied Mathematics, 2002.
    [DOI](https://doi.org/10.1137/1.9780898718027)
12. <a name="ref-12"></a> Golub, Gene H., and Charles F. Van Loan. *Matrix Computations*. 4th ed.
    Johns Hopkins University Press, 2013.
    [DOI](https://doi.org/10.56021/9781421407944) ·
    [Publisher record](https://www.press.jhu.edu/books/title/10678/matrix-computations)
13. <a name="ref-13"></a> Kalibera, Tomas, and Richard Jones. "Rigorous Benchmarking in Reasonable Time."
    *Proceedings of the 2013 International Symposium on Memory Management* (ISMM '13),
    2013: 63–74. [DOI](https://doi.org/10.1145/2464157.2464160)
14. <a name="ref-14"></a> IEEE Computer Society. "IEEE Standard for Interval Arithmetic."
    *IEEE Std 1788-2015*, 2015: 1–97.
    [DOI](https://doi.org/10.1109/IEEESTD.2015.7140721) ·
    [IEEE record](https://standards.ieee.org/ieee/1788/4431/)
15. <a name="ref-15"></a> Blue, James L. "A Portable Fortran Program to Find the Euclidean Norm of a Vector."
    *ACM Transactions on Mathematical Software* 4.1 (1978): 15–23.
    [DOI](https://doi.org/10.1145/355769.355771)
16. <a name="ref-16"></a> Kock, Anders. "Square-densities, and volume forms." Notes, December 10, 2020.
    Introduction and §1.2 (Gram's formula).
    [Author's PDF](https://math.au.dk/~kock/heron4.pdf)

[1]: #ref-1
[2]: #ref-2
[3]: #ref-3
[4]: #ref-4
[5]: #ref-5
[6]: #ref-6
[7]: #ref-7
[8]: #ref-8
[9]: #ref-9
[10]: #ref-10
[11]: #ref-11
[12]: #ref-12
[14]: #ref-14
[15]: #ref-15
[16]: #ref-16
