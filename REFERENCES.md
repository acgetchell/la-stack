# References and citations

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

`Matrix::det_errbound()` returns a conservative Shewchuk-style absolute error bound [8]
for `Matrix::det_direct()` in dimensions 2–4 when every rounded intermediate is normal
or an exact structural zero. The bound has the form
`ERR_COEFF_D · p(|A|)`, where `p(|A|) = perm(|A|)` is the absolute Leibniz sum—the
combinatorial permanent of the entrywise-absolute matrix—and
`ERR_COEFF_D ∈ {ERR_COEFF_2, ERR_COEFF_3, ERR_COEFF_4}` is a dimension-specific constant
derived from the rounding-event count of `det_direct`.
The method returns `None` when gradual underflow could violate the relative-error model.
The same bound is used internally by `det_sign_exact()`'s fast filter, but
`det_errbound()` itself is available without the `exact` feature, so downstream crates
can build custom adaptive-precision logic with pure f64 arithmetic.

### Exact determinant sign (adaptive-precision integer arithmetic)

`det_sign_exact()` uses a Shewchuk-style f64 error-bound filter [8] (the same bound exposed
by `det_errbound()` above) backed by exact `BigInt` arithmetic. Each f64 entry is decomposed
into `mantissa × 2^exponent` and scaled to a common integer base. Dimensions 0–4 use direct
integer determinant expansions; D ≥ 5 uses integer-only Bareiss elimination [7]. Neither
path constructs `BigRational` values or performs GCD normalization.
See `src/exact.rs` for the full architecture description.

### Exact linear system solve (hybrid Bareiss / BigRational)

`solve_exact()`, `solve_exact_f64()`, and `solve_exact_rounded_f64()` share the determinant
path's exact f64 decomposition and integer scaling. Matrix and RHS entries are decomposed via
IEEE 754 bit extraction [9]. Each collection is scaled independently to its own minimum
exponent, producing a `BigInt` matrix and RHS without inflating one side to accommodate the
other's range. Forward elimination runs in `BigInt` using Bareiss fraction-free updates
[7]—no `BigRational` and no GCD normalisation in the `O(D³)` phase. The upper-triangular
result is then lifted into `BigRational` for back-substitution, where fractions are inherent
and the cost is only `O(D²)`. Row swaps from first-non-zero pivoting are applied to both the
matrix and RHS. After back-substitution, multiplying by the exact power-of-two scale ratio
`2^(e_rhs − e_matrix)` recovers the solution to the original `A x = b` system.

### f64 → integer decomposition (`decompose_proven_finite_f64`)

Both the determinant and solve paths convert their finite-by-construction entries via
`decompose_proven_finite_f64`, which extracts the IEEE 754 binary64 sign, unbiased exponent,
and significand [9]. For nonzero `x`, it strips trailing zeros from the
significand so `|x| = m · 2^e` with `m` odd; signed zeros use a separate zero
component. The integer matrix is then assembled by shifting each mantissa left by
`exp − e_min`, giving a GCD-free exact-integer starting point. Solves and D ≥ 5 determinants
then apply Bareiss elimination; D ≤ 4 determinants use direct expansions. The test-only
fallible wrapper `decompose_f64` verifies rejection of non-finite raw scalars, while the
test-only `f64_to_big_rational` helper packages the same decomposition into a single
`BigRational`. See Goldberg [10] for background on floating-point representation and
conversion.

### LDLᵀ factorization (exactly symmetric positive-definite inputs)

The no-pivot LDLT implementation targets `A = L D Lᵀ` for exactly symmetric
positive-definite inputs [4-5, 11-12]. Successful construction requires every
computed diagonal pivot to be positive and greater than the caller's tolerance.
Computed zero and tolerance-small positive pivots are therefore part of the
typed diagnostic domain, not returned in a usable factorization. Because the
pivots are computed in binary64, a successful factorization is not an exact
certificate that the represented matrix is positive definite.

For pivoted variants used for symmetric *indefinite* matrices, see [6, 11-12].

### LU decomposition (Gaussian elimination with partial pivoting)

The LU implementation targets `P A = L U` and uses partial pivoting: each step
selects the remaining entry of largest magnitude in the active column. Partial
pivoting is a practical stability strategy, not an unconditional accuracy
guarantee; worst-case growth and average-case behavior are distinct concerns.

See [1-3, 11-12] for stability analysis, finite-precision behavior, and standard
algorithmic background.

## References

1. Trefethen, Lloyd N., and Robert S. Schreiber. "Average-case stability of Gaussian elimination."
   *SIAM Journal on Matrix Analysis and Applications* 11.3 (1990): 335–360.
   [DOI](https://doi.org/10.1137/0611023) ·
   [PDF](https://people.maths.ox.ac.uk/trefethen/publication/PDF/1990_44.pdf)
2. Businger, P. A. "Monitoring the Numerical Stability of Gaussian Elimination."
   *Numerische Mathematik* 16.4 (1971): 360–361.
   [DOI](https://doi.org/10.1007/BF02165006) · [Full text](https://eudml.org/doc/132040)
3. Huang, Han, and K. Tikhomirov. "Average-case analysis of the Gaussian elimination with partial pivoting."
   *Probability Theory and Related Fields* 189 (2024): 501–567.
   [DOI](https://doi.org/10.1007/s00440-024-01276-2) ·
   [Open-access article](https://link.springer.com/article/10.1007/s00440-024-01276-2) ·
   [arXiv:2206.01726](https://arxiv.org/abs/2206.01726)
4. Cholesky, André-Louis. "Sur la résolution numérique des systèmes d'équations linéaires."
   *Bulletin de la Sabix* 39 (2005): 81–95. Manuscript dated 2 December 1910.
   [DOI](https://doi.org/10.4000/sabix.529)
5. Brezinski, Claude. "La méthode de Cholesky."
   *Revue d'histoire des mathématiques* 11.2 (2005): 205–238.
   [DOI](https://doi.org/10.24033/rhm.30) ·
   [Full text](https://www.numdam.org/articles/10.24033/rhm.30/)
6. Bunch, James R., Linda Kaufman, and Beresford N. Parlett. "Decomposition of a Symmetric Matrix."
   *Numerische Mathematik* 27 (1976): 95–109.
   [DOI](https://doi.org/10.1007/BF01399088) · [Full text](https://eudml.org/doc/132435)
7. Bareiss, Erwin H. "Sylvester's Identity and Multistep Integer-Preserving Gaussian
   Elimination." *Mathematics of Computation* 22.103 (1968): 565–578.
   [DOI](https://doi.org/10.1090/S0025-5718-1968-0226829-0) ·
   [PDF](https://www.ams.org/journals/mcom/1968-22-103/S0025-5718-1968-0226829-0/S0025-5718-1968-0226829-0.pdf)
8. Shewchuk, Jonathan Richard. "Adaptive Precision Floating-Point Arithmetic and Fast
   Robust Geometric Predicates." *Discrete & Computational Geometry* 18.3 (1997): 305–363.
   [DOI](https://doi.org/10.1007/PL00009321) ·
   [PDF](https://people.eecs.berkeley.edu/~jrs/papers/robustr.pdf)
   Also: Technical Report CMU-CS-96-140, Carnegie Mellon University, May 1996.
9. IEEE Computer Society. "IEEE Standard for Floating-Point Arithmetic." *IEEE Std 754-2019*
   (Revision of IEEE 754-2008), 2019.
   [DOI](https://doi.org/10.1109/IEEESTD.2019.8766229)
   Section 3.4 (binary64 format): 1 sign bit, 11 exponent bits (bias 1023), 52 trailing
   significand bits; subnormals have biased exponent 0 with implicit leading 0.
10. Goldberg, David. "What Every Computer Scientist Should Know About Floating-Point
    Arithmetic." *ACM Computing Surveys* 23.1 (1991): 5–48.
    [DOI](https://doi.org/10.1145/103162.103163) ·
    [Authorized HTML reprint](https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html)
    Comprehensive survey of floating-point representation, rounding, and conversion.
11. Higham, Nicholas J. *Accuracy and Stability of Numerical Algorithms*. 2nd ed.
    Society for Industrial and Applied Mathematics, 2002.
    [DOI](https://doi.org/10.1137/1.9780898718027)
12. Golub, Gene H., and Charles F. Van Loan. *Matrix Computations*. 4th ed.
    Johns Hopkins University Press, 2013.
    [DOI](https://doi.org/10.56021/9781421407944) ·
    [Publisher record](https://www.press.jhu.edu/books/title/10678/matrix-computations)
13. Kalibera, Tomas, and Richard Jones. "Rigorous Benchmarking in Reasonable Time."
    *Proceedings of the 2013 International Symposium on Memory Management* (ISMM '13),
    2013: 63–74. [DOI](https://doi.org/10.1145/2464157.2464160)
