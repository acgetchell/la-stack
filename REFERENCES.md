# References and citations

## How to cite this library

If you use this library in your research or project, please cite it using the information in
[CITATION.cff](CITATION.cff). This file contains structured citation metadata that can be
processed by GitHub and other platforms.

A Zenodo DOI will be added for tagged releases.

## AI-Assisted Development Tools

- Anthropic. "Claude." <https://www.anthropic.com/claude>.
- CodeRabbit AI, Inc. "CodeRabbit." <https://coderabbit.ai/>.
- KiloCode. "KiloCode AI Engineering Assistant." <https://kilocode.ai/>.
- OpenAI. "ChatGPT." <https://openai.com/chatgpt>.
- Warp Dev, Inc. "WARP." <https://www.warp.dev/>.

All AI-generated output was reviewed and/or edited by the maintainer.
No generated content was used without human oversight.

## Linear algebra algorithms

### Absolute error bound for closed-form determinants

`Matrix::det_errbound()` returns a conservative Shewchuk-style absolute error bound [8]
for `Matrix::det_direct()` in dimensions 2–4. The bound has the form
`ERR_COEFF_D · p(|A|)`, where `p(|A|)` is the absolute Leibniz sum (the cofactor-expansion
tree with `|·|` at every leaf) and `ERR_COEFF_D ∈ {ERR_COEFF_2, ERR_COEFF_3, ERR_COEFF_4}`
is a dimension-specific constant derived from the rounding-event count of `det_direct`.
The same bound is used internally by `det_sign_exact()`'s fast filter, but
`det_errbound()` itself is available without the `exact` feature, so downstream crates
can build custom adaptive-precision logic with pure f64 arithmetic.

### Exact determinant sign (adaptive-precision Bareiss)

`det_sign_exact()` uses a Shewchuk-style f64 error-bound filter [8] (the same bound exposed
by `det_errbound()` above) backed by integer-only Bareiss elimination [7] in `BigInt`. Each
f64 entry is decomposed into `mantissa × 2^exponent`, scaled to a common integer base, and
eliminated without any `BigRational` or GCD overhead.
See `src/exact.rs` for the full architecture description.

### Exact linear system solve (hybrid Bareiss / BigRational)

`solve_exact()` and `solve_exact_f64()` share the BigInt core used for determinants.  Matrix
and RHS entries are decomposed via IEEE 754 bit extraction [9] and scaled to a shared base
`2^e_min` so the augmented system `(A | b)` becomes a `BigInt` matrix.  Forward elimination
runs in `BigInt` using Bareiss fraction-free updates [7] — no `BigRational` and no GCD
normalisation in the `O(D³)` phase.  The upper-triangular result is then lifted into
`BigRational` for back-substitution, where fractions are inherent and the cost is only
`O(D²)`.  Row swaps from first-non-zero pivoting are applied to both the matrix and the
RHS; because power-of-two scaling is applied uniformly to both sides of `A x = b`, the
solution is unchanged by the scale factor.

### f64 → integer decomposition (`f64_decompose`)

Both the determinant and solve paths convert f64 entries via `f64_decompose`, which extracts
the IEEE 754 binary64 sign, unbiased exponent, and significand [9] and strips trailing zeros
from the significand so `|x| = m · 2^e` with `m` odd.  The integer matrix is then assembled
by shifting each mantissa left by `exp − e_min`, giving a GCD-free, Bareiss-ready starting
point.  A one-shot wrapper `f64_to_bigrational` (used only in tests) packages the same
decomposition into a single `BigRational`.  See Goldberg [10] for background on IEEE 754
representation and exact rational reconstruction.

### LDL^T factorization (symmetric SPD/PSD)

The LDL^T (often abbreviated "LDLT") implementation in `la-stack` is intended for symmetric positive
definite (SPD) and positive semi-definite (PSD) matrices (e.g. Gram matrices), and does not perform
pivoting.

For background on the SPD/PSD setting, see [4-5]. For pivoted variants used for symmetric *indefinite*
matrices, see [6].

### LU decomposition (Gaussian elimination with partial pivoting)

The LU implementation in `la-stack` follows the standard Gaussian elimination / LU factorization
approach with partial pivoting for numerical stability.

See references [1-3] below.

## References

1. Trefethen, Lloyd N., and Robert S. Schreiber. "Average-case stability of Gaussian elimination."
   *SIAM Journal on Matrix Analysis and Applications* 11.3 (1990): 335–360.
   [PDF](https://people.maths.ox.ac.uk/trefethen/publication/PDF/1990_44.pdf)
2. Businger, P. A. "Monitoring the Numerical Stability of Gaussian Elimination."
   *Numerische Mathematik* 16 (1970/71): 360–361.
   [Full text](https://eudml.org/doc/132040)
3. Huang, Han, and K. Tikhomirov. "Average-case analysis of the Gaussian elimination with partial pivoting."
   *Probability Theory and Related Fields* 189 (2024): 501–567.
   [Open-access PDF](https://link.springer.com/article/10.1007/s00440-024-01276-2) (also: [arXiv:2206.01726](https://arxiv.org/abs/2206.01726))
4. Cholesky, Andre-Louis. "On the numerical solution of systems of linear equations"
   (manuscript dated 2 Dec 1910; published 2005).
   Scan + English analysis: [BibNum](https://www.bibnum.education.fr/mathematiques/algebre/sur-la-resolution-numerique-des-systemes-d-equations-lineaires)
5. Brezinski, Claude. "La methode de Cholesky." (2005).
   [PDF](https://eudml.org/doc/252115)
6. Bunch, J. R., L. Kaufman, and B. N. Parlett. "Decomposition of a Symmetric Matrix."
   *Numerische Mathematik* 27 (1976/1977): 95–110.
   [Full text](https://eudml.org/doc/132435)
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
    [PDF](https://www.validlab.com/goldberg/paper.pdf)
    Comprehensive survey of IEEE 754 representation, rounding, and exact rational
    reconstruction of floating-point values.
