#![forbid(unsafe_code)]

//! Out-of-line public-operation probes for the Rust 1.99 compiler study.
//! Compile as a library with the commands in the adjacent study; these wrappers
//! expose assembly without adding an API or changing the timed Criterion kernels.

use la_stack::{
    BigRational, DEFAULT_SINGULAR_TOL, IntervalDeterminantSign, IntervalMatrix, LaError, Matrix,
    RationalMatrix, RationalVector, ScalarWithErrorBound, Vector,
};

/// Preserve the four-element ordered FMA reduction.
pub fn dot4(left: &Vector<4>, right: &Vector<4>) -> Result<f64, LaError> {
    left.dot(right)
}

/// Preserve the five-element squared-norm reduction.
pub fn norm_squared5(vector: &Vector<5>) -> Result<f64, LaError> {
    vector.norm_squared()
}

/// Exercise the scale-aware norm and exact overflow-boundary fallback.
pub fn norm4(vector: &Vector<4>) -> Result<f64, LaError> {
    vector.norm()
}

/// Exercise the four-dimensional direct determinant.
pub fn det4(matrix: Matrix<4>) -> Result<f64, LaError> {
    matrix.det()
}

/// Exercise pivoting and substitution with the benchmark tolerance.
pub fn lu_solve4(matrix: Matrix<4>, rhs: Vector<4>) -> Result<Vector<4>, LaError> {
    matrix.lu(DEFAULT_SINGULAR_TOL)?.solve(rhs)
}

/// Exercise inclusive-range LDLT updates and substitution.
pub fn ldlt_solve4(matrix: Matrix<4>, rhs: Vector<4>) -> Result<Vector<4>, LaError> {
    matrix.ldlt(DEFAULT_SINGULAR_TOL)?.solve(rhs)
}

/// Investigate the largest measured LDLT improvement at dimension five.
pub fn ldlt_solve5(matrix: Matrix<5>, rhs: Vector<5>) -> Result<Vector<5>, LaError> {
    matrix.ldlt(DEFAULT_SINGULAR_TOL)?.solve(rhs)
}

/// Investigate the dimension-five LU timing change.
pub fn lu_solve5(matrix: Matrix<5>, rhs: Vector<5>) -> Result<Vector<5>, LaError> {
    matrix.lu(DEFAULT_SINGULAR_TOL)?.solve(rhs)
}

/// Investigate the dimension-three squared-norm timing change.
pub fn norm_squared3(vector: &Vector<3>) -> Result<f64, LaError> {
    vector.norm_squared()
}

/// Exercise the certified dot-product error-bound path.
pub fn certified_dot4(
    left: &Vector<4>,
    right: &Vector<4>,
) -> Result<Option<ScalarWithErrorBound>, LaError> {
    left.dot_with_errbound(right)
}

/// Exercise the division-free interval sign filter.
pub fn interval_det4(matrix: &IntervalMatrix<4>) -> Result<IntervalDeterminantSign, LaError> {
    matrix.det_sign()
}

/// Exercise the arbitrary-precision Bareiss control.
pub fn exact_det5(matrix: &Matrix<5>) -> Result<BigRational, LaError> {
    matrix.det_exact()
}

/// Investigate row clearing and Bareiss solve cost independently of f64.
pub fn rational_solve5(
    matrix: &RationalMatrix<5>,
    rhs: &RationalVector<5>,
) -> Result<RationalVector<5>, LaError> {
    matrix.solve(rhs)
}
