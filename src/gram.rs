#![forbid(unsafe_code)]

//! Fixed-size Gram construction.

use crate::{LaError, Matrix, Vector};

/// Construct a stack-backed [`Matrix<M>`] of pairwise vector dot products.
///
/// A Gram matrix records pairwise inner products: diagonal entries are squared
/// vector lengths, and off-diagonal entries encode their relative angles.
/// The input contains `M` finite-by-construction [`Vector<N>`] values.
/// If `V` has these vectors as rows, the mathematical Gram matrix is `G = V Vᵀ`,
/// with `G[i,j] = vectors[i] · vectors[j]`. In exact real arithmetic and for
/// `M ≤ N`, its determinant is the squared volume of the spanned
/// parallelotope. For edges from one simplex vertex, the simplex volume is
/// `sqrt(det(G)) / M!`; this also handles facets embedded in higher dimensions.
/// See `REFERENCES.md` \[16\] for the Gram determinant and volume interpretation.
///
/// Each upper-triangle dot product is computed once using [`Vector::dot`]'s
/// left-to-right fused multiply-add reduction and copied to the other triangle,
/// giving bit-for-bit symmetry. No absolute rounding-error bound is provided.
/// Rounding and underflow can destroy positive semidefiniteness or rank; this
/// operation proves neither positive definiteness nor affine independence.
/// [`Matrix::ldlt`] retains its symmetry and positive-definiteness preconditions.
/// Forming a Gram matrix squares the spectral condition number of an exact input
/// with full row rank. See the floating-point discussion in `REFERENCES.md` \[9-11\].
///
/// `M` and `N` are independent, with no dimension cap (including dimensions
/// through 8); storage is `O(M²)` and work is `O(M²(N + 1))`, including output
/// initialization when `N = 0`. `M = 0` returns an empty
/// matrix; `N = 0` returns an all-zero matrix. No optional feature is required.
///
/// # Errors
/// Returns [`LaError::NonFinite`] if a dot-product accumulator overflows, even
/// when the exact result would be finite after cancellation. The error preserves
/// [`Vector::dot`]'s [`Computation`](crate::NonFiniteOrigin::Computation) origin,
/// [`VectorDotProduct`](crate::ArithmeticOperation::VectorDotProduct) operation,
/// and first failing reduction [`Step`](crate::NonFiniteLocation::Step).
/// The step indexes the vector coordinate, not the output matrix cell.
/// Pairs are visited in upper-triangle row order.
///
/// # Examples
/// ```
/// use la_stack::prelude::*;
/// # fn main() -> Result<(), LaError> {
/// let vectors = [Vector::try_new([1.0, 0.0, 0.0])?,
///                Vector::try_new([0.0, 2.0, 0.0])?];
/// let gram = gram_matrix(&vectors)?;
/// assert_eq!(gram.as_rows(), &[[1.0, 0.0], [0.0, 4.0]]);
/// assert_eq!(gram.ldlt(Tolerance::try_new(0.0)?)?.det()?, 4.0);
/// # Ok(())
/// # }
/// ```
pub const fn gram_matrix<const M: usize, const N: usize>(
    vectors: &[Vector<N>; M],
) -> Result<Matrix<M>, LaError> {
    let mut rows = [[0.0; M]; M];
    let mut i = 0;
    while i < M {
        let mut j = i;
        while j < M {
            let value = match vectors[i].dot(&vectors[j]) {
                Ok(value) => value,
                Err(error) => return Err(error),
            };
            rows[i][j] = value;
            rows[j][i] = value;
            j += 1;
        }
        i += 1;
    }
    Matrix::try_from_rows(rows)
}
