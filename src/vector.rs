//! Fixed-size, stack-allocated vectors.

use core::hint::cold_path;

use crate::LaError;

/// Finite fixed-size vector of length `D`, stored inline.
///
/// Public construction rejects NaN and infinity through [`try_new`](Self::try_new),
/// and the storage field is private, so a `Vector` value carries the invariant
/// that every stored entry is finite. Algorithms therefore do not re-scan stored
/// entries before using a `Vector`; they only report non-finite values computed
/// during arithmetic, such as overflowed accumulators.
///
/// Direct field construction is intentionally unavailable to downstream callers:
///
/// ```compile_fail
/// use la_stack::Vector;
///
/// let _ = Vector::<2> {
///     data: [1.0, f64::NAN],
/// };
/// ```
#[must_use]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vector<const D: usize> {
    data: [f64; D],
}

/// Fixed-size vector whose stored entries are all finite.
///
/// This internal proof wrapper is used where carrying the invariant explicitly
/// makes algorithm boundaries clearer. Public callers use [`Vector`], which is
/// already finite by construction.
#[must_use]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteVector<const D: usize> {
    vector: Vector<D>,
}

impl<const D: usize> FiniteVector<D> {
    /// Construct a finite vector without checking the invariant.
    ///
    /// This crate-internal escape hatch is only for paths with a local proof
    /// that every stored entry is finite.
    #[inline]
    pub(crate) const fn new_unchecked(vector: Vector<D>) -> Self {
        Self { vector }
    }

    /// Consume the wrapper and return the underlying raw vector.
    #[inline]
    pub(crate) const fn into_vector(self) -> Vector<D> {
        self.vector
    }

    /// Borrow the backing array without revalidating stored entries.
    #[inline]
    #[must_use]
    pub(crate) const fn as_array(&self) -> &[f64; D] {
        self.vector.as_array()
    }

    /// Consume and return the backing array.
    #[inline]
    #[must_use]
    pub(crate) const fn into_array(self) -> [f64; D] {
        self.vector.into_array()
    }

    /// Dot product for already finite vectors.
    ///
    /// Stored entries are known finite, so this path only checks whether the
    /// accumulation overflows to NaN or infinity.
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] when the accumulated dot product overflows
    /// to NaN or infinity.
    #[inline]
    pub(crate) const fn dot(self, other: Self) -> Result<f64, LaError> {
        let lhs = self.as_array();
        let rhs = other.as_array();
        let mut acc = 0.0;
        let mut i = 0;
        while i < D {
            acc = lhs[i].mul_add(rhs[i], acc);
            if !acc.is_finite() {
                cold_path();
                return Err(LaError::non_finite_at(i));
            }
            i += 1;
        }
        Ok(acc)
    }

    /// Squared Euclidean norm for an already finite vector.
    ///
    /// Stored entries are known finite, so this path only checks whether the
    /// accumulation overflows to NaN or infinity.
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] when the accumulated norm overflows to NaN
    /// or infinity.
    #[inline]
    pub(crate) const fn norm2_sq(self) -> Result<f64, LaError> {
        let data = self.as_array();
        let mut acc = 0.0;
        let mut i = 0;
        while i < D {
            acc = data[i].mul_add(data[i], acc);
            if !acc.is_finite() {
                cold_path();
                return Err(LaError::non_finite_at(i));
            }
            i += 1;
        }
        Ok(acc)
    }
}

impl<const D: usize> Vector<D> {
    /// Test-only infallible constructor for finite literal fixtures.
    #[cfg(test)]
    #[inline]
    pub(crate) const fn new(data: [f64; D]) -> Self {
        match Self::try_new(data) {
            Ok(vector) => vector,
            Err(_) => panic!("Vector::new requires finite entries"),
        }
    }

    /// Try to create a finite vector from a backing array.
    ///
    /// This is the public raw-storage boundary for vectors. Successful
    /// construction makes the returned [`Vector`] a finite-storage proof.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let v = Vector::<3>::try_new([1.0, 2.0, 3.0])?;
    /// assert_eq!(v.into_array(), [1.0, 2.0, 3.0]);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] with the first offending entry index when
    /// `data` contains NaN or infinity.
    #[inline]
    pub const fn try_new(data: [f64; D]) -> Result<Self, LaError> {
        if let Some(index) = Self::first_non_finite_entry_in(&data) {
            Err(LaError::non_finite_at(index))
        } else {
            Ok(Self::new_unchecked(data))
        }
    }

    /// Construct a vector without checking that entries are finite.
    ///
    /// This crate-internal escape hatch is reserved for finite literals and
    /// algorithm outputs whose finite invariant is visible at the call site.
    /// Computed outputs must check their intermediates before using this
    /// constructor to create an observable [`Vector`].
    #[inline]
    pub(crate) const fn new_unchecked(data: [f64; D]) -> Self {
        Self { data }
    }

    /// Return the first non-finite stored entry in index order.
    ///
    /// Shared by the public raw-storage boundary and the internal finite-vector
    /// proof wrapper so both paths report the same first offending index with
    /// [`LaError::NonFinite`].
    const fn first_non_finite_entry_in(data: &[f64; D]) -> Option<usize> {
        let mut i = 0;
        while i < D {
            if !data[i].is_finite() {
                return Some(i);
            }
            i += 1;
        }
        None
    }

    /// All-zeros finite vector.
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
        Self::new_unchecked([0.0; D])
    }

    /// Borrow the finite backing array.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let v = Vector::<2>::try_new([1.0, -2.0])?;
    /// assert_eq!(v.as_array(), &[1.0, -2.0]);
    /// # Ok(())
    /// # }
    /// ```
    #[inline]
    #[must_use]
    pub const fn as_array(&self) -> &[f64; D] {
        &self.data
    }

    /// Consume and return the finite backing array.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let v = Vector::<2>::try_new([1.0, 2.0])?;
    /// let a = v.into_array();
    /// assert_eq!(a, [1.0, 2.0]);
    /// # Ok(())
    /// # }
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
    /// certified absolute rounding bound for the returned dot product. Raw
    /// `Vector` values are finite by construction, so this method only checks
    /// whether the accumulation overflows to NaN or infinity.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let a = Vector::<3>::try_new([1.0, 2.0, 3.0])?;
    /// let b = Vector::<3>::try_new([-2.0, 0.5, 4.0])?;
    /// assert!((a.dot(b)? - 11.0).abs() <= 1e-12);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] when the accumulated dot product overflows
    /// to NaN or infinity.
    #[inline]
    pub const fn dot(self, other: Self) -> Result<f64, LaError> {
        FiniteVector::new_unchecked(self).dot(FiniteVector::new_unchecked(other))
    }

    /// Squared Euclidean norm.
    ///
    /// This is computed as `dot(self, self)`, so `norm2_sq` has the same
    /// `f64` [`mul_add`](f64::mul_add) accumulation behavior as [`dot`](Self::dot).
    /// Intermediate rounding occurs, and this method does not provide a
    /// certified absolute rounding bound for the returned squared norm.
    /// `Vector` values are finite by construction, so this method only checks
    /// whether the accumulation overflows to NaN or infinity.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let v = Vector::<3>::try_new([1.0, 2.0, 3.0])?;
    /// assert!((v.norm2_sq()? - 14.0).abs() <= 1e-12);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] when the accumulated norm overflows to NaN
    /// or infinity.
    #[inline]
    pub const fn norm2_sq(self) -> Result<f64, LaError> {
        FiniteVector::new_unchecked(self).norm2_sq()
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
                fn [<public_api_vector_try_new_rejects_nonfinite_ $d d>]() {
                    let mut a_arr = [1.0f64; $d];
                    a_arr[$d - 1] = f64::NAN;

                    assert_eq!(
                        Vector::<$d>::try_new(a_arr),
                        Err(LaError::NonFinite {
                            row: None,
                            col: $d - 1,
                        })
                    );
                }

                #[test]
                fn [<public_api_vector_try_new_rejects_nonfinite_rhs_fixture_ $d d>]() {
                    let mut b_arr = [1.0f64; $d];
                    b_arr[0] = f64::INFINITY;

                    assert_eq!(
                        Vector::<$d>::try_new(b_arr),
                        Err(LaError::NonFinite { row: None, col: 0 })
                    );
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

                    let finite = FiniteVector::<$d>::new_unchecked(Vector::<$d>::new(arr));

                    assert_array_abs_eq(&finite.into_array(), &arr);
                    assert_array_abs_eq(&finite.into_vector().into_array(), &arr);
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
