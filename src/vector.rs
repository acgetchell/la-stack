//! Fixed-size, stack-allocated vectors.

/// Fixed-size vector of length `D`, stored inline.
#[must_use]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vector<const D: usize> {
    pub(crate) data: [f64; D],
}

impl<const D: usize> Vector<D> {
    /// Create a vector from a backing array.
    ///
    /// # Examples
    /// ```
    /// #![allow(unused_imports)]
    /// use la_stack::prelude::*;
    ///
    /// let v = Vector::<3>::new([1.0, 2.0, 3.0]);
    /// assert_eq!(v.into_array(), [1.0, 2.0, 3.0]);
    /// ```
    #[inline]
    pub const fn new(data: [f64; D]) -> Self {
        Self { data }
    }

    /// All-zeros vector.
    ///
    /// # Examples
    /// ```
    /// #![allow(unused_imports)]
    /// use la_stack::prelude::*;
    ///
    /// let z = Vector::<2>::zero();
    /// assert_eq!(z.into_array(), [0.0, 0.0]);
    /// ```
    #[inline]
    pub const fn zero() -> Self {
        Self { data: [0.0; D] }
    }

    /// Borrow the backing array.
    ///
    /// # Examples
    /// ```
    /// #![allow(unused_imports)]
    /// use la_stack::prelude::*;
    ///
    /// let v = Vector::<2>::new([1.0, -2.0]);
    /// assert_eq!(v.as_array(), &[1.0, -2.0]);
    /// ```
    #[inline]
    #[must_use]
    pub const fn as_array(&self) -> &[f64; D] {
        &self.data
    }

    /// Consume and return the backing array.
    ///
    /// # Examples
    /// ```
    /// #![allow(unused_imports)]
    /// use la_stack::prelude::*;
    ///
    /// let v = Vector::<2>::new([1.0, 2.0]);
    /// let a = v.into_array();
    /// assert_eq!(a, [1.0, 2.0]);
    /// ```
    #[inline]
    #[must_use]
    pub const fn into_array(self) -> [f64; D] {
        self.data
    }

    /// Dot product.
    ///
    /// # Examples
    /// ```
    /// #![allow(unused_imports)]
    /// use la_stack::prelude::*;
    ///
    /// let a = Vector::<3>::new([1.0, 2.0, 3.0]);
    /// let b = Vector::<3>::new([-2.0, 0.5, 4.0]);
    /// assert!((a.dot(b) - 11.0).abs() <= 1e-12);
    /// ```
    #[inline]
    #[must_use]
    pub fn dot(self, other: Self) -> f64 {
        let mut acc = 0.0;
        let mut i = 0;
        while i < D {
            acc = self.data[i].mul_add(other.data[i], acc);
            i += 1;
        }
        acc
    }

    /// Squared Euclidean norm.
    ///
    /// # Examples
    /// ```
    /// #![allow(unused_imports)]
    /// use la_stack::prelude::*;
    ///
    /// let v = Vector::<3>::new([1.0, 2.0, 3.0]);
    /// assert!((v.norm2_sq() - 14.0).abs() <= 1e-12);
    /// ```
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
        // Use black_box to avoid constant-folding/inlining eliminating the actual dot loop,
        // which can make coverage tools report the mul_add line as uncovered.
        use core::hint::black_box;

        let a = Vector::<3>::new([black_box(1.0), black_box(2.0), black_box(3.0)]);
        let b = Vector::<3>::new([black_box(-2.0), black_box(0.5), black_box(4.0)]);

        // Call via (black_boxed) fn pointers to discourage inlining, improving line-level coverage
        // attribution for the loop body.
        let dot_fn: fn(Vector<3>, Vector<3>) -> f64 = black_box(Vector::<3>::dot);
        let norm2_sq_fn: fn(Vector<3>) -> f64 = black_box(Vector::<3>::norm2_sq);

        assert_approx(dot_fn(black_box(a), black_box(b)), 11.0, 0.0);
        assert_approx(norm2_sq_fn(black_box(a)), 14.0, 0.0);
    }

    #[test]
    fn zero_as_array_into_array_and_default() {
        let z = Vector::<3>::zero();
        for &x in z.as_array() {
            assert_approx(x, 0.0, 0.0);
        }
        for x in z.into_array() {
            assert_approx(x, 0.0, 0.0);
        }

        let d = Vector::<3>::default();
        for x in d.into_array() {
            assert_approx(x, 0.0, 0.0);
        }
    }
}
