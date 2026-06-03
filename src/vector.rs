//! Fixed-size, stack-allocated vectors.

use crate::LaError;

/// Fixed-size vector of length `D`, stored inline.
#[must_use]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vector<const D: usize> {
    pub(crate) data: [f64; D],
}

/// Fixed-size vector whose stored entries are all finite.
///
/// This proof-carrying wrapper lets callers validate raw [`Vector`] values once
/// at an input boundary, then pass the finite invariant into numerical
/// algorithms without rediscovering that no stored entry is NaN or infinite.
#[must_use]
#[derive(Clone, Copy, Debug, PartialEq)]
#[allow(clippy::redundant_pub_crate)]
pub(crate) struct FiniteVector<const D: usize> {
    vector: Vector<D>,
}

impl<const D: usize> FiniteVector<D> {
    /// Construct a finite vector without checking the invariant.
    ///
    /// This is crate-internal so public callers must use [`try_new`](Self::try_new)
    /// or [`from_array`](Self::from_array), which preserve diagnostics for
    /// rejected raw entries.
    #[inline]
    pub(crate) const fn new_unchecked(vector: Vector<D>) -> Self {
        Self { vector }
    }

    /// Validate that every stored vector entry is finite.
    /// # Errors
    /// Returns [`LaError::NonFinite`] with `row: None` and the first offending
    /// entry index when `vector` contains NaN or infinity.
    #[inline]
    pub const fn try_new(vector: Vector<D>) -> Result<Self, LaError> {
        let mut i = 0;
        while i < D {
            if !vector.data[i].is_finite() {
                return Err(LaError::non_finite_at(i));
            }
            i += 1;
        }
        Ok(Self::new_unchecked(vector))
    }

    /// Validate raw vector storage and construct a finite vector.
    /// # Errors
    /// Returns [`LaError::NonFinite`] with `row: None` and the first offending
    /// entry index when `data` contains NaN or infinity.
    #[inline]
    pub const fn from_array(data: [f64; D]) -> Result<Self, LaError> {
        Self::try_new(Vector::new(data))
    }

    /// All-zeros finite vector.
    #[inline]
    pub const fn zero() -> Self {
        Self::new_unchecked(Vector::zero())
    }

    /// Consume the wrapper and return the underlying raw vector.
    #[inline]
    pub const fn into_vector(self) -> Vector<D> {
        self.vector
    }

    /// Borrow the backing array without revalidating stored entries.
    #[cfg(feature = "exact")]
    #[inline]
    pub(crate) const fn as_array(&self) -> &[f64; D] {
        self.vector.as_array()
    }

    /// Consume and return the backing array.
    #[inline]
    #[must_use]
    pub const fn into_array(self) -> [f64; D] {
        self.vector.into_array()
    }
}

impl<const D: usize> Default for FiniteVector<D> {
    #[inline]
    fn default() -> Self {
        Self::zero()
    }
}

impl<const D: usize> From<FiniteVector<D>> for Vector<D> {
    #[inline]
    fn from(value: FiniteVector<D>) -> Self {
        value.into_vector()
    }
}

impl<const D: usize> TryFrom<Vector<D>> for FiniteVector<D> {
    type Error = LaError;

    #[inline]
    fn try_from(value: Vector<D>) -> Result<Self, Self::Error> {
        Self::try_new(value)
    }
}

impl<const D: usize> TryFrom<[f64; D]> for FiniteVector<D> {
    type Error = LaError;

    #[inline]
    fn try_from(value: [f64; D]) -> Result<Self, Self::Error> {
        Self::from_array(value)
    }
}

impl<const D: usize> Vector<D> {
    /// Create a vector from a backing array.
    ///
    /// # Examples
    /// ```
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
    /// Terms are accumulated in `f64` using [`f64::mul_add`] at each index.
    /// Intermediate rounding occurs, and this method does not provide a
    /// certified absolute rounding bound for the returned dot product.  The
    /// returned [`Result`] is still checked for non-finite inputs and for
    /// non-finite accumulation.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let a = Vector::<3>::new([1.0, 2.0, 3.0]);
    /// let b = Vector::<3>::new([-2.0, 0.5, 4.0]);
    /// assert!((a.dot(b)? - 11.0).abs() <= 1e-12);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] when either input contains NaN or infinity,
    /// or when the accumulated dot product overflows to NaN or infinity.
    #[inline]
    pub const fn dot(self, other: Self) -> Result<f64, LaError> {
        let mut acc = 0.0;
        let mut i = 0;
        while i < D {
            if !self.data[i].is_finite() {
                return Err(LaError::non_finite_at(i));
            }
            if !other.data[i].is_finite() {
                return Err(LaError::non_finite_at(i));
            }
            acc = self.data[i].mul_add(other.data[i], acc);
            if !acc.is_finite() {
                return Err(LaError::non_finite_at(i));
            }
            i += 1;
        }
        Ok(acc)
    }

    /// Squared Euclidean norm.
    ///
    /// This is computed as `dot(self, self)`, so `norm2_sq` has the same
    /// `f64` [`mul_add`](f64::mul_add) accumulation behavior as [`dot`](Self::dot).
    /// Intermediate rounding occurs, and this method does not provide a
    /// certified absolute rounding bound for the returned squared norm.  The
    /// returned [`Result`] is still checked for non-finite inputs and for
    /// non-finite accumulation.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let v = Vector::<3>::new([1.0, 2.0, 3.0]);
    /// assert!((v.norm2_sq()? - 14.0).abs() <= 1e-12);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] when the input contains NaN or infinity,
    /// or when the accumulated norm overflows to NaN or infinity.
    #[inline]
    pub const fn norm2_sq(self) -> Result<f64, LaError> {
        let mut acc = 0.0;
        let mut i = 0;
        while i < D {
            if !self.data[i].is_finite() {
                return Err(LaError::non_finite_at(i));
            }
            acc = self.data[i].mul_add(self.data[i], acc);
            if !acc.is_finite() {
                return Err(LaError::non_finite_at(i));
            }
            i += 1;
        }
        Ok(acc)
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

    fn assert_array_abs_eq<const D: usize>(actual: &[f64; D], expected: &[f64; D]) {
        for i in 0..D {
            assert_abs_diff_eq!(actual[i], expected[i], epsilon = 0.0);
        }
    }

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
                    let dot_fn: fn(Vector<$d>, Vector<$d>) -> Result<f64, LaError> =
                        black_box(Vector::<$d>::dot);
                    let norm2_sq_fn: fn(Vector<$d>) -> Result<f64, LaError> =
                        black_box(Vector::<$d>::norm2_sq);

                    assert_abs_diff_eq!(
                        dot_fn(black_box(a), black_box(b)).unwrap(),
                        expected_dot,
                        epsilon = 1e-14
                    );
                    assert_abs_diff_eq!(
                        norm2_sq_fn(black_box(a)).unwrap(),
                        expected_norm2_sq,
                        epsilon = 1e-14
                    );
                }

                #[test]
                fn [<public_api_vector_dot_and_norm2_sq_reject_nonfinite_ $d d>]() {
                    let mut a_arr = [1.0f64; $d];
                    a_arr[$d - 1] = f64::NAN;
                    let a = Vector::<$d>::new(a_arr);
                    let b = Vector::<$d>::new([1.0; $d]);

                    assert_eq!(
                        a.dot(b),
                        Err(LaError::NonFinite {
                            row: None,
                            col: $d - 1,
                        })
                    );
                    assert_eq!(
                        a.norm2_sq(),
                        Err(LaError::NonFinite {
                            row: None,
                            col: $d - 1,
                        })
                    );
                }

                #[test]
                fn [<public_api_vector_dot_rejects_nonfinite_rhs_ $d d>]() {
                    let a = Vector::<$d>::new([1.0; $d]);
                    let mut b_arr = [1.0f64; $d];
                    b_arr[0] = f64::INFINITY;
                    let b = Vector::<$d>::new(b_arr);

                    assert_eq!(a.dot(b), Err(LaError::NonFinite { row: None, col: 0 }));
                }

                #[test]
                fn [<public_api_vector_dot_and_norm2_sq_reject_overflow_ $d d>]() {
                    let mut a_arr = [1.0f64; $d];
                    a_arr[0] = f64::MAX;
                    let a = Vector::<$d>::new(a_arr);

                    let mut b_arr = [1.0f64; $d];
                    b_arr[0] = 2.0;
                    let b = Vector::<$d>::new(b_arr);

                    assert_eq!(a.dot(b), Err(LaError::NonFinite { row: None, col: 0 }));
                    assert_eq!(a.norm2_sq(), Err(LaError::NonFinite { row: None, col: 0 }));
                }

                #[test]
                fn [<finite_vector_accepts_and_roundtrips_ $d d>]() {
                    let arr = {
                        let mut arr = [0.0f64; $d];
                        let values = [1.0f64, 2.0, 3.0, 4.0, 5.0];
                        for (dst, src) in arr.iter_mut().zip(values.iter()) {
                            *dst = *src;
                        }
                        arr
                    };

                    let finite = FiniteVector::<$d>::from_array(arr).unwrap();

                    assert_array_abs_eq(&finite.into_array(), &arr);
                    assert_array_abs_eq(&finite.into_vector().into_array(), &arr);
                    assert_array_abs_eq(&FiniteVector::<$d>::zero().into_array(), &[0.0; $d]);
                    assert_array_abs_eq(Vector::from(finite).as_array(), &arr);
                    assert_array_abs_eq(
                        &FiniteVector::<$d>::try_from(Vector::<$d>::new(arr)).unwrap().into_array(),
                        &arr
                    );
                    assert_array_abs_eq(&FiniteVector::<$d>::try_from(arr).unwrap().into_array(), &arr);
                }

                #[test]
                fn [<finite_vector_rejects_nonfinite_with_index_ $d d>]() {
                    let mut arr = [1.0f64; $d];
                    arr[$d - 1] = f64::INFINITY;

                    assert_eq!(
                        FiniteVector::<$d>::from_array(arr),
                        Err(LaError::NonFinite {
                            row: None,
                            col: $d - 1,
                        })
                    );
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
