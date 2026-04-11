# References and citations

## How to cite this library

If you use this library in your research or project, please cite it using the information in
[CITATION.cff](CITATION.cff). This file contains structured citation metadata that can be
processed by GitHub and other platforms.

A Zenodo DOI will be added for tagged releases.

## Linear algebra algorithms

### LU decomposition (Gaussian elimination with partial pivoting)

The LU implementation in `la-stack` follows the standard Gaussian elimination / LU factorization
approach with partial pivoting for numerical stability.

See references [1-3] below.

### LDL^T factorization (symmetric SPD/PSD)

The LDL^T (often abbreviated "LDLT") implementation in `la-stack` is intended for symmetric positive
definite (SPD) and positive semi-definite (PSD) matrices (e.g. Gram matrices), and does not perform
pivoting.

For background on the SPD/PSD setting, see [4-5]. For pivoted variants used for symmetric *indefinite*
matrices, see [6].

### Exact determinant sign (adaptive-precision Bareiss)

`det_sign_exact()` uses a Shewchuk-style f64 error-bound filter [8] backed by integer-only
Bareiss elimination [7] in `BigInt`. Each f64 entry is decomposed into `mantissa × 2^exponent`,
scaled to a common integer base, and eliminated without any `BigRational` or GCD overhead.
See `src/exact.rs` for the full architecture description.

## References

### LU / Gaussian elimination

1. Trefethen, Lloyd N., and Robert S. Schreiber. "Average-case stability of Gaussian elimination."
   *SIAM Journal on Matrix Analysis and Applications* 11.3 (1990): 335–360.
   [PDF](https://people.maths.ox.ac.uk/trefethen/publication/PDF/1990_44.pdf)
2. Businger, P. A. "Monitoring the Numerical Stability of Gaussian Elimination."
   *Numerische Mathematik* 16 (1970/71): 360–361.
   [Full text](https://eudml.org/doc/132040)
3. Huang, Han, and K. Tikhomirov. "Average-case analysis of the Gaussian elimination with partial pivoting."
   *Probability Theory and Related Fields* 189 (2024): 501–567.
   [Open-access PDF](https://link.springer.com/article/10.1007/s00440-024-01276-2) (also: [arXiv:2206.01726](https://arxiv.org/abs/2206.01726))

### LDL^T / Cholesky (symmetric SPD/PSD)

4. Cholesky, Andre-Louis. "On the numerical solution of systems of linear equations"
   (manuscript dated 2 Dec 1910; published 2005).
   Scan + English analysis: [BibNum](https://www.bibnum.education.fr/mathematiques/algebre/sur-la-resolution-numerique-des-systemes-d-equations-lineaires)
5. Brezinski, Claude. "La methode de Cholesky." (2005).
   [PDF](https://eudml.org/doc/252115)

### Pivoted LDL^T (symmetric indefinite)

6. Bunch, J. R., L. Kaufman, and B. N. Parlett. "Decomposition of a Symmetric Matrix."
   *Numerische Mathematik* 27 (1976/1977): 95–110.
   [Full text](https://eudml.org/doc/132435)

### Exact determinant (`det_exact`, `det_exact_f64`, `det_sign_exact`)

7. Bareiss, Erwin H. "Sylvester's Identity and Multistep Integer-Preserving Gaussian
   Elimination." *Mathematics of Computation* 22.103 (1968): 565–578.
   [DOI](https://doi.org/10.1090/S0025-5718-1968-0226829-0) ·
   [PDF](https://www.ams.org/journals/mcom/1968-22-103/S0025-5718-1968-0226829-0/S0025-5718-1968-0226829-0.pdf)
8. Shewchuk, Jonathan Richard. "Adaptive Precision Floating-Point Arithmetic and Fast
   Robust Geometric Predicates." *Discrete & Computational Geometry* 18.3 (1997): 305–363.
   [DOI](https://doi.org/10.1007/PL00009321) ·
   [PDF](https://people.eecs.berkeley.edu/~jrs/papers/robustr.pdf)
   Also: Technical Report CMU-CS-96-140, Carnegie Mellon University, May 1996.

### f64 → BigRational conversion (`f64_to_bigrational`)

`f64_to_bigrational` converts an f64 to an exact `BigRational` by decomposing the IEEE 754
binary64 bit representation into its sign, exponent, and significand fields.  Because every
finite f64 is exactly `±m × 2^e` (where `m` is an integer), the rational can be constructed
directly via `BigRational::new_raw` without GCD normalization — trailing zeros in the
significand are stripped first so the fraction is already in lowest terms.

See references [9-10] below.

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
