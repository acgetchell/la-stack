#![forbid(unsafe_code)]

// ruleid: la-stack.rust.no-public-raw-linear-algebra-modules
pub mod matrix;

// ruleid: la-stack.rust.no-public-raw-linear-algebra-modules
pub mod vector;

// ok: la-stack.rust.no-public-raw-linear-algebra-modules
mod exact;

#[must_use]
pub struct Matrix<const D: usize> {
    // ruleid: la-stack.rust.no-public-matrix-vector-storage-fields
    pub rows: [[f64; D]; D],
}

#[must_use]
pub struct Vector<const D: usize> {
    // ruleid: la-stack.rust.no-public-matrix-vector-storage-fields
    pub(crate) data: [f64; D],
}

pub struct CleanMatrix<const D: usize> {
    // ok: la-stack.rust.no-public-matrix-vector-storage-fields
    rows: [[f64; D]; D],
}

impl<const D: usize> Matrix<D> {
    // ruleid: la-stack.rust.no-public-unchecked-finite-constructors
    pub const fn from_rows_unchecked(rows: [[f64; D]; D]) -> Self {
        Self { rows }
    }

    // ruleid: la-stack.rust.no-public-unchecked-finite-constructors
    pub(crate) const fn from_rows_unchecked_internal(rows: [[f64; D]; D]) -> Self {
        Self { rows }
    }

    // ruleid: la-stack.rust.no-public-unchecked-finite-constructors
    pub(crate) const fn rows_mut_unchecked(&mut self) -> &mut [[f64; D]; D] {
        &mut self.rows
    }

    // ruleid: la-stack.rust.no-public-unchecked-finite-constructors
    pub(super) const fn from_rows_unchecked_super(rows: [[f64; D]; D]) -> Self {
        Self { rows }
    }

    // ruleid: la-stack.rust.no-public-unchecked-finite-constructors
    pub(in crate::matrix) const fn rows_mut_unchecked_scoped(&mut self) -> &mut [[f64; D]; D] {
        &mut self.rows
    }

    // ok: la-stack.rust.no-public-unchecked-finite-constructors
    const fn new_unchecked_private(rows: [[f64; D]; D]) -> Self {
        Self { rows }
    }

    // ruleid: la-stack.rust.det-sign-exact-must-be-infallible
    pub fn det_sign_exact(&self) -> Result<DeterminantSign, LaError> {
        Ok(DeterminantSign::Zero)
    }
}

impl<const D: usize> Vector<D> {
    // ruleid: la-stack.rust.no-public-unchecked-finite-constructors
    pub const fn new_unchecked(data: [f64; D]) -> Self {
        Self { data }
    }

    // ruleid: la-stack.rust.no-public-unchecked-finite-constructors
    pub(crate) const fn new_unchecked_internal(data: [f64; D]) -> Self {
        Self { data }
    }
}

pub enum DeterminantSign {
    Zero,
}

pub enum LaError {}

pub struct InfallibleMatrix;

impl InfallibleMatrix {
    // ok: la-stack.rust.det-sign-exact-must-be-infallible
    pub fn det_sign_exact(&self) -> DeterminantSign {
        DeterminantSign::Zero
    }
}

// ruleid: la-stack.rust.exact-benchmark-validation-must-return-proof
pub fn validate_exact_fixture<const D: usize>(_input: &ExactInput<D>) {}

mod invalid_validation_return {
    use super::ExactInput;

    // ruleid: la-stack.rust.exact-benchmark-validation-must-return-proof
    pub fn validate_exact_fixture<const D: usize>(_input: ExactInput<D>) -> bool {
        true
    }
}

mod valid_validation_return {
    use super::{ExactInput, ValidatedExactInput};

    // ok: la-stack.rust.exact-benchmark-validation-must-return-proof
    pub fn validate_exact_fixture<const D: usize>(
        _input: ExactInput<D>,
    ) -> ValidatedExactInput<D> {
        todo!()
    }
}

pub struct ExactInput<const D: usize> {
    matrix: Matrix<D>,
    rhs: Vector<D>,
}

pub struct ValidatedExactInput<const D: usize> {
    // ruleid: la-stack.rust.validated-exact-input-fields-private
    pub matrix: Matrix<D>,
    rhs: Vector<D>,
}

pub struct CleanValidatedExactInput<const D: usize> {
    // ok: la-stack.rust.validated-exact-input-fields-private
    matrix: Matrix<D>,
    rhs: Vector<D>,
}

mod private_validated_fields {
    use super::{Matrix, Vector};

    pub struct ValidatedExactInput<const D: usize> {
        // ok: la-stack.rust.validated-exact-input-fields-private
        pub(self) matrix: Matrix<D>,
        rhs: Vector<D>,
    }
}

// ruleid: la-stack.rust.exact-benchmark-helpers-require-validated-input
fn bench_raw_exact_input<const D: usize>(_input: &ExactInput<D>) {}

// ok: la-stack.rust.exact-benchmark-helpers-require-validated-input
fn bench_validated_exact_input<const D: usize>(_input: &ValidatedExactInput<D>) {}

pub struct Lu<const D: usize>;

impl<const D: usize> Lu<D> {
    // ruleid: la-stack.rust.no-legacy-solve-vec-api
    pub fn solve_vec(&self, rhs: Vector<D>) -> Vector<D> {
        rhs
    }

    // ok: la-stack.rust.no-legacy-solve-vec-api
    pub fn solve(&self, rhs: Vector<D>) -> Vector<D> {
        rhs
    }
}

fn call_sites<const D: usize>(lu: &Lu<D>, rhs: Vector<D>) {
    // ruleid: la-stack.rust.no-legacy-solve-vec-api
    let _ = lu.solve_vec(rhs);
}
