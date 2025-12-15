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

    use core::hint::black_box;

    use approx::assert_abs_diff_eq;
    use pastey::paste;

    macro_rules! gen_public_api_vector_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<public_api_vector_new_as_array_into_array_ $d d>]() {
                    let arr = {
                        let mut arr = [0.0f64; $d];
                        let values = [1.0f64, 2.0, 3.0, 4.0, 5.0];
                        for (dst, src) in arr.iter_mut().zip(values.iter()) {
                            *dst = *src;
                        }
                        arr
                    };

                    let v = Vector::<$d>::new(arr);

                    for i in 0..$d {
                        assert_abs_diff_eq!(v.as_array()[i], arr[i], epsilon = 0.0);
                    }

                    let out = v.into_array();
                    for i in 0..$d {
                        assert_abs_diff_eq!(out[i], arr[i], epsilon = 0.0);
                    }
                }

                #[test]
                fn [<public_api_vector_zero_as_array_into_array_default_ $d d>]() {
                    let z = Vector::<$d>::zero();
                    for &x in z.as_array() {
                        assert_abs_diff_eq!(x, 0.0, epsilon = 0.0);
                    }
                    for x in z.into_array() {
                        assert_abs_diff_eq!(x, 0.0, epsilon = 0.0);
                    }

                    let d = Vector::<$d>::default();
                    for x in d.into_array() {
                        assert_abs_diff_eq!(x, 0.0, epsilon = 0.0);
                    }
                }

                #[test]
                fn [<public_api_vector_dot_and_norm2_sq_ $d d>]() {
                    // Use black_box to avoid constant-folding/inlining eliminating the actual dot loop,
                    // which can make coverage tools report the mul_add line as uncovered.

                    let a_arr = {
                        let mut arr = [0.0f64; $d];
                        let values = [1.0f64, 2.0, 3.0, 4.0, 5.0];
                        for (dst, src) in arr.iter_mut().zip(values.iter()) {
                            *dst = black_box(*src);
                        }
                        arr
                    };
                    let b_arr = {
                        let mut arr = [0.0f64; $d];
                        let values = [-2.0f64, 0.5, 4.0, -1.0, 2.0];
                        for (dst, src) in arr.iter_mut().zip(values.iter()) {
                            *dst = black_box(*src);
                        }
                        arr
                    };

                    let expected_dot = {
                        let mut acc = 0.0;
                        let mut i = 0;
                        while i < $d {
                            acc = a_arr[i].mul_add(b_arr[i], acc);
                            i += 1;
                        }
                        acc
                    };
                    let expected_norm2_sq = {
                        let mut acc = 0.0;
                        let mut i = 0;
                        while i < $d {
                            acc = a_arr[i].mul_add(a_arr[i], acc);
                            i += 1;
                        }
                        acc
                    };

                    let a = Vector::<$d>::new(black_box(a_arr));
                    let b = Vector::<$d>::new(black_box(b_arr));

                    // Call via (black_boxed) fn pointers to discourage inlining, improving line-level coverage
                    // attribution for the loop body.
                    let dot_fn: fn(Vector<$d>, Vector<$d>) -> f64 = black_box(Vector::<$d>::dot);
                    let norm2_sq_fn: fn(Vector<$d>) -> f64 = black_box(Vector::<$d>::norm2_sq);

                    assert_abs_diff_eq!(dot_fn(black_box(a), black_box(b)), expected_dot, epsilon = 1e-14);
                    assert_abs_diff_eq!(norm2_sq_fn(black_box(a)), expected_norm2_sq, epsilon = 1e-14);
                }
            }
        };
    }

    // Mirror delaunay-style multi-dimension tests.
    gen_public_api_vector_tests!(2);
    gen_public_api_vector_tests!(3);
    gen_public_api_vector_tests!(4);
    gen_public_api_vector_tests!(5);
}
