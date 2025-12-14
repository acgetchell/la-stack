//! Fixed-size, stack-allocated vectors.

/// Fixed-size vector of length `D`, stored inline.
#[must_use]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vector<const D: usize> {
    pub(crate) data: [f64; D],
}

impl<const D: usize> Vector<D> {
    /// Create a vector from a backing array.
    #[inline]
    pub const fn new(data: [f64; D]) -> Self {
        Self { data }
    }

    /// All-zeros vector.
    #[inline]
    pub const fn zero() -> Self {
        Self { data: [0.0; D] }
    }

    /// Borrow the backing array.
    #[inline]
    #[must_use]
    pub const fn as_array(&self) -> &[f64; D] {
        &self.data
    }

    /// Consume and return the backing array.
    #[inline]
    #[must_use]
    pub const fn into_array(self) -> [f64; D] {
        self.data
    }

    /// Dot product.
    #[inline]
    #[must_use]
    pub fn dot(self, other: Self) -> f64 {
        let mut acc = 0.0;
        for i in 0..D {
            acc = self.data[i].mul_add(other.data[i], acc);
        }
        acc
    }

    /// Squared Euclidean norm.
    #[inline]
    #[must_use]
    pub fn norm2_sq(self) -> f64 {
        self.dot(self)
    }
}

impl<const D: usize> Default for Vector<D> {
    #[inline]
    fn default() -> Self {
        Self::zero()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_approx(a: f64, b: f64, eps: f64) {
        assert!((a - b).abs() <= eps, "{a} !~= {b} (eps={eps})");
    }

    #[test]
    fn dot_and_norm2_sq() {
        let a = Vector::<3>::new([1.0, 2.0, 3.0]);
        let b = Vector::<3>::new([-2.0, 0.5, 4.0]);
        assert_approx(a.dot(b), 11.0, 0.0);
        assert_approx(a.norm2_sq(), 14.0, 0.0);
    }
}
