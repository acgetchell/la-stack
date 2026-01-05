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

- Golub, Gene H., and Charles F. Van Loan. *Matrix Computations* (4th ed.). Johns Hopkins University Press, 2013.
- Higham, Nicholas J. *Accuracy and Stability of Numerical Algorithms* (2nd ed.). SIAM, 2002.
- Trefethen, Lloyd N., and David Bau III. *Numerical Linear Algebra*. SIAM, 1997.

### LDL^T factorization (symmetric SPD/PSD)

The LDL^T (often abbreviated "LDLT") implementation in `la-stack` is intended for symmetric positive
definite (SPD) and positive semi-definite (PSD) matrices (e.g. Gram matrices), and does not perform
pivoting.

- Golub, Gene H., and Charles F. Van Loan. *Matrix Computations* (4th ed.). Johns Hopkins University Press, 2013.
- Higham, Nicholas J. *Accuracy and Stability of Numerical Algorithms* (2nd ed.). SIAM, 2002.
- Trefethen, Lloyd N., and David Bau III. *Numerical Linear Algebra*. SIAM, 1997.
