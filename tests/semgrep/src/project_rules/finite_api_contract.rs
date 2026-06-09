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

    // ok: la-stack.rust.no-public-unchecked-finite-constructors
    pub(crate) const fn from_rows_unchecked_internal(rows: [[f64; D]; D]) -> Self {
        Self { rows }
    }
}

impl<const D: usize> Vector<D> {
    // ruleid: la-stack.rust.no-public-unchecked-finite-constructors
    pub const fn new_unchecked(data: [f64; D]) -> Self {
        Self { data }
    }

    // ok: la-stack.rust.no-public-unchecked-finite-constructors
    pub(crate) const fn new_unchecked_internal(data: [f64; D]) -> Self {
        Self { data }
    }
}

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
