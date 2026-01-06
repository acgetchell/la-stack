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
