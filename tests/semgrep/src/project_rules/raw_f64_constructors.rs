#![forbid(unsafe_code)]
#![allow(dead_code)]

#[must_use]
pub struct Matrix<const D: usize> {
    rows: [[f64; D]; D],
}

#[must_use]
pub struct Vector<const D: usize> {
    data: [f64; D],
}

#[non_exhaustive]
pub enum LaError {
    NonFinite,
}

impl<const D: usize> Matrix<D> {
    // ruleid: la-stack.rust.no-public-infallible-raw-f64-constructors
    pub const fn from_rows(rows: [[f64; D]; D]) -> Self {
        Self { rows }
    }

    // ok: la-stack.rust.no-public-infallible-raw-f64-constructors
    pub const fn try_from_rows(rows: [[f64; D]; D]) -> Result<Self, LaError> {
        Ok(Self { rows })
    }

    // ok: la-stack.rust.no-public-infallible-raw-f64-constructors
    pub(crate) const fn from_rows_literal(rows: [[f64; D]; D]) -> Self {
        Self { rows }
    }
}

impl<const D: usize> Vector<D> {
    // ruleid: la-stack.rust.no-public-infallible-raw-f64-constructors
    pub const fn new(data: [f64; D]) -> Self {
        Self { data }
    }

    // ok: la-stack.rust.no-public-infallible-raw-f64-constructors
    pub const fn try_new(data: [f64; D]) -> Result<Self, LaError> {
        Ok(Self { data })
    }

    // ok: la-stack.rust.no-public-infallible-raw-f64-constructors
    pub(crate) const fn new_literal(data: [f64; D]) -> Self {
        Self { data }
    }
}
