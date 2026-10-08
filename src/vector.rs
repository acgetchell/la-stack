#![forbid(unsafe_code)]

//! Fixed-size, stack-allocated vectors.

use core::hint::cold_path;

use crate::norm::norm_near_overflow;
use crate::rounding::{compare_product_with_rounded, two_sum_error};
use crate::{ArithmeticOperation, LaError};

/// A scalar estimate paired with a certified absolute error bound.
///
/// Values of this type are produced by [`Vector::dot_with_errbound`] and
/// [`Vector::dot_difference_with_errbound`]. The exact-real value of the
/// corresponding expression over the stored binary64 inputs lies between
/// [`lower_bound`](Self::lower_bound) and [`upper_bound`](Self::upper_bound),
/// and differs from [`estimate`](Self::estimate) by at most
/// [`absolute_error_bound`](Self::absolute_error_bound).
///
/// The bound certifies floating-point roundoff in one specified arithmetic
/// tree. It is not a caller-selected numerical tolerance and does not classify
/// an interval containing zero as equality.
///
/// Callers cannot construct this type directly. Every value has a finite
/// estimate, a finite non-negative absolute error bound, and finite ordered
/// endpoints.
///
/// # Examples
/// ```
/// use core::assert_matches;
/// use la_stack::prelude::*;
///
/// # fn main() -> Result<(), LaError> {
/// let left = Vector::<2>::try_new([1.0, 2.0])?;
/// let right = Vector::<2>::try_new([3.0, 4.0])?;
/// let certificate = left.dot_with_errbound(&right)?;
/// assert_eq!(certificate.map(ScalarWithErrorBound::estimate), Some(11.0));
/// // Preserve None: without a certificate, an exact fallback is needed
/// // before making a claim about the exact-real result.
/// let enclosure = certificate.map(|bounded| (bounded.lower_bound(), bounded.upper_bound()));
/// assert_matches!(enclosure, Some((lower, upper)) if lower <= 11.0 && 11.0 <= upper);
/// # Ok(())
/// # }
/// ```
#[must_use]
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScalarWithErrorBound {
    estimate: f64,
    absolute_error_bound: f64,
}

impl ScalarWithErrorBound {
    /// Return the rounded scalar estimate.
    #[inline]
    #[must_use]
    pub const fn estimate(self) -> f64 {
        self.estimate
    }

    /// Return the certified absolute error bound.
    #[inline]
    #[must_use]
    pub const fn absolute_error_bound(self) -> f64 {
        self.absolute_error_bound
    }

    /// Return a finite outward-rounded lower bound on the exact-real value.
    #[expect(
        clippy::inline_always,
        reason = "avoid a measured out-of-line endpoint call for each projection"
    )]
    #[inline(always)]
    #[must_use]
    pub const fn lower_bound(self) -> f64 {
        if self.absolute_error_bound == 0.0 {
            return self.estimate;
        }
        let rounded = self.estimate - self.absolute_error_bound;
        if two_sum_error(self.estimate, -self.absolute_error_bound, rounded) < 0.0 {
            rounded.next_down()
        } else {
            rounded
        }
    }

    /// Return a finite outward-rounded upper bound on the exact-real value.
    #[expect(
        clippy::inline_always,
        reason = "avoid a measured out-of-line endpoint call for each projection"
    )]
    #[inline(always)]
    #[must_use]
    pub const fn upper_bound(self) -> f64 {
        if self.absolute_error_bound == 0.0 {
            return self.estimate;
        }
        let rounded = self.estimate + self.absolute_error_bound;
        if two_sum_error(self.estimate, self.absolute_error_bound, rounded) > 0.0 {
            rounded.next_up()
        } else {
            rounded
        }
    }

    /// Construct a validated public result with finite outward endpoints.
    ///
    /// The exact `|estimate| + bound` dominates both endpoint magnitudes.
    /// A rounded sum below `MAX` has a finite successor; at `MAX` its residual
    /// must be nonpositive. This proves both outward endpoints finite without
    /// computing an endpoint the caller may never request.
    #[inline]
    const fn try_new(estimate: f64, absolute_error_bound: f64) -> Option<Self> {
        if !estimate.is_finite() || !absolute_error_bound.is_finite() || absolute_error_bound < 0.0
        {
            return None;
        }

        if absolute_error_bound == 0.0 {
            return Some(Self {
                estimate,
                absolute_error_bound: 0.0,
            });
        }

        let magnitude = estimate.abs();
        let outer_rounded = magnitude + absolute_error_bound;
        if !outer_rounded.is_finite()
            || (outer_rounded >= f64::MAX
                && two_sum_error(magnitude, absolute_error_bound, outer_rounded) > 0.0)
        {
            return None;
        }

        Some(Self {
            estimate,
            absolute_error_bound,
        })
    }
}

/// State for one certified left-to-right FMA reduction.
///
/// The rounded estimate continues accumulating after `proof_available` becomes
/// false. This lets the public methods distinguish a non-finite estimate
/// (`Err`) from a finite estimate whose proof arithmetic is unavailable
/// (`Ok(None)`). `magnitude_upper` encloses the sum of exact product
/// magnitudes while that proof remains available.
/// While proof is available, `estimate` is normal or positive zero.
/// `single_exact_product` records one nonzero product with a unit multiplier;
/// adding zero products preserves that fact and the exact estimate.
#[derive(Clone, Copy, Debug, PartialEq)]
struct CertifiedReduction {
    estimate: f64,
    magnitude_upper: f64,
    proof_available: bool,
    single_exact_product: bool,
}

impl CertifiedReduction {
    const ZERO: Self = Self {
        estimate: 0.0,
        magnitude_upper: 0.0,
        proof_available: true,
        single_exact_product: false,
    };

    /// Add one exact-real product through one rounded FMA.
    ///
    /// A non-finite estimate becomes a typed public error. Underflow or range
    /// loss confined to proof construction instead clears `proof_available`,
    /// preserving the finite estimate for the eventual `Ok(None)` result.
    // Keep the reduction state in caller registers instead of returning a
    // Result-sized state aggregate across a call for every coordinate.
    #[expect(
        clippy::inline_always,
        reason = "measured per-coordinate reduction-state copy overhead"
    )]
    #[inline(always)]
    const fn add_product(
        mut self,
        left: f64,
        right: f64,
        operation: ArithmeticOperation,
        index: usize,
    ) -> Result<Self, LaError> {
        if left == 0.0 || right == 0.0 {
            // While proof is available, the estimate is normal or +0: the
            // initial estimate is +0 and exact cancellation rounds to +0.
            // Adding either signed zero therefore leaves every bit unchanged.
            // After proof loss an underflow may have produced -0, so retain
            // the FMA there to preserve the specified signed-zero result.
            if !self.proof_available {
                self.estimate = left.mul_add(right, self.estimate);
            }
            return Ok(self);
        }
        let prior = self.estimate;
        let estimate = left.mul_add(right, prior);
        if !estimate.is_normal() {
            if !estimate.is_finite() {
                cold_path();
                return Err(LaError::non_finite_computation_step(operation, index));
            }
            if self.proof_available {
                self.proof_available =
                    estimate == 0.0 && Self::fma_result_is_exact_zero(left, right, prior);
            }
        }
        if self.proof_available {
            match Self::add_product_magnitude_upper(self.magnitude_upper, left, right) {
                Some((magnitude_upper, single_exact_product)) => {
                    self.magnitude_upper = magnitude_upper;
                    self.single_exact_product = single_exact_product;
                }
                None => self.proof_available = false,
            }
        }
        self.estimate = estimate;
        Ok(self)
    }

    /// Return whether `left × right + addend` is exactly zero.
    ///
    /// Nonzero operands are guaranteed by `add_product`. This distinguishes
    /// exact cancellation from a nonzero value rounded to
    /// zero. Only exact zero is admissible under the relative-error model used
    /// by the public certified reductions.
    #[inline]
    const fn fma_result_is_exact_zero(left: f64, right: f64, addend: f64) -> bool {
        let rounded_product = left * right;
        let rounded_bits = rounded_product.to_bits();
        let negated_addend_bits = (-addend).to_bits();
        let same_rounded_value = rounded_bits == negated_addend_bits
            || (rounded_bits << 1 == 0 && negated_addend_bits << 1 == 0);
        rounded_product.is_finite()
            && same_rounded_value
            && compare_product_with_rounded(left, right, rounded_product) == 0
    }

    /// Add an upward-rounded bound on `|left × right|` to the magnitude sum.
    ///
    /// The caller supplies nonzero finite operands.
    /// The returned flag identifies a first unit product; later products clear it.
    ///
    /// `None` means a nonzero product was not normal or the product/sum could
    /// not be enclosed by a finite binary64 value, so the public result must be
    /// proof-unavailable.
    #[inline]
    const fn add_product_magnitude_upper(
        magnitude_upper: f64,
        left: f64,
        right: f64,
    ) -> Option<(f64, bool)> {
        let left_magnitude = left.abs();
        let right_magnitude = right.abs();
        let rounded_product = left_magnitude * right_magnitude;
        if !rounded_product.is_normal() {
            return None;
        }

        if magnitude_upper == 0.0 {
            let exact = left_magnitude.to_bits() == 1.0_f64.to_bits()
                || right_magnitude.to_bits() == 1.0_f64.to_bits();
            // Multiplication by ±1 is exact. Otherwise an unconditional upward
            // step bounds the product without computing its rounding residual.
            let product_upper = if exact {
                rounded_product
            } else {
                rounded_product.next_up()
            };
            return if product_upper.is_finite() {
                Some((product_upper, exact))
            } else {
                None
            };
        }
        // One upward-widened FMA encloses magnitude_upper + |left * right|.
        // It is no wider than first rounding the product upward, adding, and
        // widening that sum: round-to-nearest is monotone. Keep the separate
        // product-normality check above for the relative-error proof contract.
        let rounded_sum = left_magnitude.mul_add(right_magnitude, magnitude_upper);
        if !rounded_sum.is_finite() {
            return None;
        }
        let sum_upper = rounded_sum.next_up();
        if sum_upper.is_finite() {
            Some((sum_upper, false))
        } else {
            None
        }
    }

    /// Finish the reduction with an upward-rounded `gamma_n` error bound.
    ///
    /// Returns `None` when an earlier proof step failed, a needed nonzero bound
    /// has a term count outside the relative-error model, or bound/endpoints cannot
    /// remain finite. Otherwise the result satisfies every public
    /// [`ScalarWithErrorBound`] invariant.
    #[expect(
        clippy::cast_precision_loss,
        reason = "a usable gamma requires a term count below 2^53, where the cast is exact"
    )]
    #[inline]
    const fn finish(self, term_count: Option<usize>) -> Option<ScalarWithErrorBound> {
        if !self.proof_available {
            return None;
        }
        if self.magnitude_upper == 0.0 || self.single_exact_product {
            return ScalarWithErrorBound::try_new(self.estimate, 0.0);
        }

        let Some(term_count) = term_count else {
            return None;
        };
        let scaled_roundoff = (term_count as f64) * (f64::EPSILON / 2.0);
        if !scaled_roundoff.is_finite() || scaled_roundoff >= 1.0 {
            return None;
        }

        // The count conversion, multiplication by 2^-53, and subtraction from
        // one are exact throughout the usable range. Round the division and
        // final multiplication upward to retain a certified upper bound.
        let gamma = scaled_roundoff / (1.0 - scaled_roundoff);
        let gamma_upper = gamma.next_up();
        if !gamma_upper.is_finite() {
            return None;
        }
        let rounded_bound = gamma_upper * self.magnitude_upper;
        if !rounded_bound.is_finite() {
            return None;
        }
        let absolute_error_bound = if rounded_bound == 0.0 {
            0.0
        } else {
            rounded_bound.next_up()
        };
        ScalarWithErrorBound::try_new(self.estimate, absolute_error_bound)
    }
}

/// Finite fixed-size vector of length `D`, stored inline.
///
/// Public construction rejects NaN and infinity through [`try_new`](Self::try_new),
/// and the storage field is private, so a `Vector` value carries the invariant
/// that every stored entry is finite. Algorithms therefore do not re-scan stored
/// entries at every use; user-visible non-finite errors come from construction
/// boundaries or from values computed during arithmetic, such as overflowed
/// accumulators.
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
        if let Some(index) = Self::first_non_finite_entry(&data) {
            Err(LaError::non_finite_input_vector(index))
        } else {
            Ok(Self { data })
        }
    }

    /// Finalize vector storage produced by an arithmetic operation.
    ///
    /// Keeping this validation in the type that owns the finite-storage
    /// invariant prevents a new computation path from accidentally turning raw
    /// non-finite storage into a [`Vector`].
    #[inline]
    pub(crate) const fn from_computation(
        data: [f64; D],
        operation: ArithmeticOperation,
    ) -> Result<Self, LaError> {
        if let Some(index) = Self::first_non_finite_entry(&data) {
            Err(LaError::non_finite_computation_step(operation, index))
        } else {
            Ok(Self { data })
        }
    }

    /// Return the first non-finite stored entry in index order.
    ///
    /// Used by the public raw-storage boundary to report the first offending
    /// index with [`LaError::NonFinite`].
    const fn first_non_finite_entry(data: &[f64; D]) -> Option<usize> {
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
        Self { data: [0.0; D] }
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

    /// Unsigned angle to another nonzero vector, in radians in `[0, π]`.
    ///
    /// Shares the allocation-free numerical implementation and rounded semantics
    /// of [`crate::angle_between`], using this type's finite-storage invariant.
    /// Arbitrary finite magnitudes are supported even when a norm would overflow.
    /// No certified absolute error bound or exact parallelism test is provided.
    ///
    /// # Errors
    /// Returns [`LaError::EmptyVector`] for `D = 0`, [`LaError::ZeroVector`] for
    /// an all-zero operand (`self` first), or [`LaError::NonFinite`] for a
    /// non-finite computed angle, with vector-angle computation provenance.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let x = Vector::try_new([1.0, 0.0, 0.0])?;
    /// let y = Vector::try_new([1.0, 2.0_f64.powi(-30), 0.0])?;
    /// assert!((x.angle(&y)? - 2.0_f64.powi(-30)).abs() < 1e-24);
    /// # Ok(())
    /// # }
    /// ```
    #[inline]
    pub fn angle(&self, other: &Self) -> Result<f64, LaError> {
        crate::angle::angle_finite(&self.data, &other.data)
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
    /// assert!((a.dot(&b)? - 11.0).abs() <= 1e-12);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] when the accumulated dot product overflows
    /// to NaN or infinity.
    #[inline]
    pub const fn dot(&self, other: &Self) -> Result<f64, LaError> {
        self.dot_with_operation(other, ArithmeticOperation::VectorDotProduct)
    }

    /// Dot product with a certified absolute roundoff bound.
    ///
    /// The estimate uses the deterministic left-to-right recurrence
    /// `s[0] = 0` and `s[i + 1] = self[i].mul_add(other[i], s[i])`. When the
    /// relative-error model is valid, the returned certificate bounds the
    /// difference between `s[D]` and the exact-real expression
    /// `Σᵢ self[i] × other[i]` over the stored binary64 inputs.
    ///
    /// The general bound is `gamma_D × Σᵢ |self[i] × other[i]|`, where
    /// `gamma_D = D u / (1 - D u)` and `u = 2^-53`. The magnitude sum and the
    /// published bound are rounded upward. See `REFERENCES.md` \[9-11\].
    /// A single nonzero normal product with a unit (±1) multiplier receives a zero
    /// bound: every other FMA adds exact zero and preserves that product.
    ///
    /// `Ok(None)` means no certificate is available because a nonzero product
    /// or FMA result entered the subnormal range, the reduction dimension made
    /// `gamma_D` invalid, or a proof-only magnitude/bound calculation exhausted
    /// the finite binary64 range. It does not mean the exact dot product is
    /// zero. Unlike a user-selected tolerance, a returned error bound describes
    /// rounding in this specific arithmetic tree.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let left = Vector::<3>::try_new([1.0, 2.0, 3.0])?;
    /// let right = Vector::<3>::try_new([4.0, 5.0, 6.0])?;
    /// let positive = left.dot_with_errbound(&right)?.and_then(|bounded| {
    ///     if bounded.lower_bound() > 0.0 {
    ///         Some(true)
    ///     } else if bounded.upper_bound() <= 0.0 {
    ///         Some(false)
    ///     } else {
    ///         None // The enclosure cannot establish whether the result is positive.
    ///     }
    /// });
    /// assert_eq!(positive, Some(true));
    ///
    /// // Finite inputs can also lack a certificate; use an exact fallback
    /// // before deciding the sign in this case.
    /// let tiny = Vector::<2>::try_new([1e-200, 1e-200])?;
    /// assert_eq!(tiny.dot_with_errbound(&tiny)?, None);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] with the first failing reduction index and
    /// [`ArithmeticOperation::VectorDotProduct`] when an FMA estimate becomes
    /// non-finite.
    #[expect(
        clippy::inline_always,
        reason = "eliminate measured Result-return overhead in repeated projections"
    )]
    #[inline(always)]
    pub const fn dot_with_errbound(
        &self,
        other: &Self,
    ) -> Result<Option<ScalarWithErrorBound>, LaError> {
        let left = self.as_array();
        let right = other.as_array();
        let mut reduction = CertifiedReduction::ZERO;
        let mut i = 0;
        while i < D {
            reduction = match reduction.add_product(
                left[i],
                right[i],
                ArithmeticOperation::VectorDotProduct,
                i,
            ) {
                Ok(reduction) => reduction,
                Err(error) => return Err(error),
            };
            i += 1;
        }
        Ok(reduction.finish(Some(D)))
    }

    /// Certified dot product with an unrounded vector difference.
    ///
    /// This evaluates the exact-real expression
    /// `Σᵢ self[i] × (left[i] - right[i])` without first rounding
    /// `left - right` into a [`Vector`]. Its deterministic arithmetic tree is
    ///
    /// ```text
    /// s[0]       = 0
    /// s[2i + 1]  = self[i].mul_add(left[i], s[2i])
    /// s[2i + 2]  = (-self[i]).mul_add(right[i], s[2i + 1]).
    /// ```
    ///
    /// A returned certificate therefore includes all `2D` FMA rounding events
    /// in that tree and bounds the intended expression over the original
    /// binary64 coordinates. The general absolute bound is
    /// `gamma_2D × Σᵢ (|self[i] × left[i]| + |self[i] × right[i]|)`, where
    /// `gamma_2D = 2D u / (1 - 2D u)` and `u = 2^-53`; every magnitude and the
    /// final bound are rounded upward. Its
    /// [`lower_bound`](ScalarWithErrorBound::lower_bound) and
    /// [`upper_bound`](ScalarWithErrorBound::upper_bound) can certify a sign or
    /// separation from a caller's threshold. An overlapping endpoint range
    /// remains inconclusive and should trigger the caller's exact fallback.
    /// As with [`dot_with_errbound`](Self::dot_with_errbound), a single nonzero
    /// normal product with a unit (±1) multiplier receives a zero bound.
    ///
    /// `Ok(None)` has the same proof-unavailable meaning as in
    /// [`dot_with_errbound`](Self::dot_with_errbound), including gradual
    /// underflow and proof-only range exhaustion.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let axis = Vector::<2>::try_new([2.0, -1.0])?;
    /// let left = Vector::<2>::try_new([4.0, 1.0])?;
    /// let right = Vector::<2>::try_new([1.0, 3.0])?;
    /// let separated = axis.dot_difference_with_errbound(&left, &right)?.and_then(|bounded| {
    ///     if bounded.lower_bound() > 1.0 {
    ///         Some(true)
    ///     } else if bounded.upper_bound() <= 1.0 {
    ///         Some(false)
    ///     } else {
    ///         None // An exact fallback is needed to decide this threshold test.
    ///     }
    /// });
    /// // An unavailable certificate also remains None through and_then.
    /// assert_eq!(separated, Some(true));
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] with the first failing coordinate index
    /// and [`ArithmeticOperation::VectorDotDifference`] when either FMA for
    /// that coordinate produces a non-finite estimate.
    #[inline]
    pub const fn dot_difference_with_errbound(
        &self,
        left: &Self,
        right: &Self,
    ) -> Result<Option<ScalarWithErrorBound>, LaError> {
        let axis = self.as_array();
        let left = left.as_array();
        let right = right.as_array();
        let mut reduction = CertifiedReduction::ZERO;
        let mut i = 0;
        while i < D {
            reduction = match reduction.add_product(
                axis[i],
                left[i],
                ArithmeticOperation::VectorDotDifference,
                i,
            ) {
                Ok(reduction) => reduction,
                Err(error) => return Err(error),
            };
            reduction = match reduction.add_product(
                -axis[i],
                right[i],
                ArithmeticOperation::VectorDotDifference,
                i,
            ) {
                Ok(reduction) => reduction,
                Err(error) => return Err(error),
            };
            i += 1;
        }
        Ok(reduction.finish(D.checked_mul(2)))
    }

    /// Accumulate a dot product while retaining the public operation that owns it.
    const fn dot_with_operation(
        &self,
        other: &Self,
        operation: ArithmeticOperation,
    ) -> Result<f64, LaError> {
        let lhs = self.as_array();
        let rhs = other.as_array();
        let mut acc = 0.0;
        let mut i = 0;
        while i < D {
            acc = lhs[i].mul_add(rhs[i], acc);
            i += 1;
        }
        if acc.is_finite() {
            Ok(acc)
        } else {
            cold_path();
            Err(Self::dot_non_finite_error(lhs, rhs, operation))
        }
    }

    /// Replay a non-finite dot product to locate the first failing step.
    ///
    /// This runs only after the success-path traversal has produced a non-finite
    /// final accumulator. Stored entries are finite, so once a fused multiply-add
    /// produces a non-finite accumulator, later steps cannot make it finite again.
    /// Replaying the same left-to-right operations must therefore find the first
    /// failing index.
    #[cold]
    const fn dot_non_finite_error(
        lhs: &[f64; D],
        rhs: &[f64; D],
        operation: ArithmeticOperation,
    ) -> LaError {
        let mut acc = 0.0;
        let mut i = 0;
        let last = D.saturating_sub(1);
        while i < last {
            acc = lhs[i].mul_add(rhs[i], acc);
            if !acc.is_finite() {
                return LaError::non_finite_computation_step(operation, i);
            }
            i += 1;
        }

        LaError::non_finite_computation_step(operation, last)
    }

    /// Squared Euclidean norm.
    ///
    /// This is computed as `dot(self, self)`, so `norm_squared` has the same
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
    /// assert!((v.norm_squared()? - 14.0).abs() <= 1e-12);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] when the accumulated norm overflows to NaN
    /// or infinity.
    #[inline]
    pub const fn norm_squared(&self) -> Result<f64, LaError> {
        self.dot_with_operation(self, ArithmeticOperation::VectorSquaredNorm)
    }

    /// Overflow- and underflow-safe Euclidean norm.
    ///
    /// This computes `sqrt(Σᵢ self[i]²)` with a deterministic left-to-right
    /// scaled sum-of-squares recurrence. Each non-zero magnitude is divided by
    /// the largest magnitude seen so far before it is squared, so intermediate
    /// squares cannot overflow and an all-subnormal vector is scaled into the
    /// normal range. See `REFERENCES.md` \[15\].
    ///
    /// The divisions, fused multiply-adds, square root, and final rescaling are
    /// rounded in binary64. This method does not claim correct rounding or
    /// provide a certified absolute error bound. Because [`Vector`] entries are
    /// finite by construction, it returns a finite non-negative result unless
    /// the exact Euclidean norm rounds outside the finite binary64 range.
    /// Near that boundary, a fixed-size stack accumulator sums the coordinate
    /// squares exactly and compares squared rounding midpoints. This fallback
    /// prevents accumulated roundoff from causing or hiding overflow and does
    /// not require the `exact` feature.
    ///
    /// Unlike [`norm_squared`](Self::norm_squared), this method does not require the
    /// squared norm to be representable. For example, the norm of
    /// `[1.0e200, 1.0e200]` is finite even though its squared norm is not.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let ordinary = Vector::<2>::try_new([3.0, 4.0])?;
    /// assert_eq!(ordinary.norm()?, 5.0);
    /// assert_eq!(ordinary.norm_squared()?, 25.0);
    ///
    /// let large = Vector::<2>::try_new([1.0e200, 1.0e200])?;
    /// assert!(large.norm()?.is_finite());
    /// assert_eq!(
    ///     large.norm_squared(),
    ///     Err(LaError::non_finite_computation_step(ArithmeticOperation::VectorSquaredNorm, 0)),
    /// );
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] with
    /// [`ArithmeticOperation::VectorNorm`] when the exact Euclidean norm rounds
    /// to infinity under round-to-nearest, ties-to-even.
    #[inline]
    pub fn norm(&self) -> Result<f64, LaError> {
        let mut entries = self.as_array().iter();
        // The first coordinate establishes the scale without a division or FMA.
        // A zero (or absent) first coordinate preserves the empty-prefix state.
        let mut scale = entries.next().copied().unwrap_or(0.0).abs();
        let mut scaled_sum = 1.0;

        for &entry in entries {
            let magnitude = entry.abs();
            if magnitude == 0.0 {
                continue;
            }

            if scale < magnitude {
                let ratio = scale / magnitude;
                scaled_sum = (scaled_sum * ratio).mul_add(ratio, 1.0);
                scale = magnitude;
            } else {
                let ratio = magnitude / scale;
                scaled_sum = ratio.mul_add(ratio, scaled_sum);
            }
        }

        // With b = bit_length(D), D < 2^b. If scale <= 2^(1023-b),
        // even the L1 upper bound D*scale is below 2^1023. Both the exact
        // norm and the rounded recurrence therefore have ample range margin.
        // Checking scale, rather than only the computed norm, also catches
        // true overflow that rounding in the recurrence could hide.
        let dimension_bits = D.bit_width();
        let safe_scale = f64::from_bits(u64::from(2046 - dimension_bits) << 52);
        if scale > safe_scale {
            return norm_near_overflow(self.as_array(), scale);
        }
        Ok(scale * scaled_sum.sqrt())
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
    use core::hint::black_box;

    use approx::assert_abs_diff_eq;
    use pastey::paste;

    use super::*;

    fn assert_certified_proof_loss_survives_normal_terms<const D: usize>() {
        let mut left_data = [0.0; D];
        left_data[0] = f64::MIN_POSITIVE;
        left_data[D - 1] = 1.0;
        let mut right_data = [1.0; D];
        right_data[0] = 0.5;
        let left = Vector::new(left_data);
        let right = Vector::new(right_data);

        // The first product is subnormal; the last restores a normal estimate,
        // but cannot restore the lost certificate for the complete reduction.
        assert_abs_diff_eq!(left.dot(&right).unwrap(), 1.0, epsilon = 0.0);
        assert_eq!(left.dot_with_errbound(&right), Ok(None));
        assert_eq!(
            left.dot_difference_with_errbound(&right, &Vector::zero()),
            Ok(None)
        );
    }

    fn assert_certified_proof_loss_preserves_later_overflow<const D: usize>() {
        let mut axis_data = [0.0; D];
        axis_data[0] = f64::MIN_POSITIVE;
        axis_data[D - 1] = f64::MAX;
        let axis = Vector::new(axis_data);
        let mut left_data = [0.0; D];
        left_data[0] = 0.5;
        left_data[D - 1] = 2.0;

        assert_eq!(
            axis.dot_with_errbound(&Vector::new(left_data)),
            Err(LaError::non_finite_computation_step(
                ArithmeticOperation::VectorDotProduct,
                D - 1,
            ))
        );
        let difference_error = Err(LaError::non_finite_computation_step(
            ArithmeticOperation::VectorDotDifference,
            D - 1,
        ));
        assert_eq!(
            axis.dot_difference_with_errbound(&Vector::new(left_data), &Vector::zero()),
            difference_error
        );

        // Overflow in the second FMA must also remain observable after the
        // earlier coordinate has already made certification unavailable.
        left_data[D - 1] = 0.0;
        let mut right_data = [0.0; D];
        right_data[D - 1] = -2.0;
        assert_eq!(
            axis.dot_difference_with_errbound(&Vector::new(left_data), &Vector::new(right_data)),
            difference_error
        );
    }

    macro_rules! gen_certified_reduction_sequence_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<certified_proof_loss_survives_normal_terms_ $d d>]() {
                    assert_certified_proof_loss_survives_normal_terms::<$d>();
                }

                #[test]
                fn [<certified_proof_loss_preserves_later_overflow_ $d d>]() {
                    assert_certified_proof_loss_preserves_later_overflow::<$d>();
                }
            }
        };
    }

    gen_certified_reduction_sequence_tests!(2);
    gen_certified_reduction_sequence_tests!(3);
    gen_certified_reduction_sequence_tests!(4);
    gen_certified_reduction_sequence_tests!(5);

    macro_rules! gen_vector_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<vector_new_as_array_into_array_ $d d>]() {
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
                fn [<vector_zero_as_array_into_array_default_ $d d>]() {
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
                fn [<vector_dot_and_norm_squared_ $d d>]() {
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
                    let expected_norm_squared = {
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
                    let dot_fn: fn(&Vector<$d>, &Vector<$d>) -> Result<f64, LaError> =
                        black_box(Vector::<$d>::dot);
                    let norm_squared_fn: fn(&Vector<$d>) -> Result<f64, LaError> =
                        black_box(Vector::<$d>::norm_squared);

                    assert_abs_diff_eq!(
                        dot_fn(black_box(&a), black_box(&b)).unwrap(),
                        expected_dot,
                        epsilon = 1e-14
                    );
                    assert_abs_diff_eq!(
                        norm_squared_fn(black_box(&a)).unwrap(),
                        expected_norm_squared,
                        epsilon = 1e-14
                    );
                }

                #[test]
                fn [<vector_certified_dot_and_difference_ $d d>]() {
                    let mut left_data = [0.0; $d];
                    let mut right_data = [0.0; $d];
                    let left_values = [1.0, 2.0, 3.0, 4.0, 5.0];
                    let right_values = [2.0, 3.0, 4.0, 5.0, 6.0];
                    for (destination, source) in left_data.iter_mut().zip(left_values) {
                        *destination = source;
                    }
                    for (destination, source) in right_data.iter_mut().zip(right_values) {
                        *destination = source;
                    }
                    let left = Vector::<$d>::new(left_data);
                    let right = Vector::<$d>::new(right_data);
                    let zero = Vector::<$d>::zero();

                    let dot = left.dot(&right).unwrap();
                    let dot_bound = left.dot_with_errbound(&right).unwrap().unwrap();
                    assert_abs_diff_eq!(dot_bound.estimate(), dot, epsilon = 0.0);
                    assert!(dot_bound.absolute_error_bound() >= 0.0);
                    assert!(dot_bound.lower_bound() <= dot);
                    assert!(dot <= dot_bound.upper_bound());

                    let difference_bound = left
                        .dot_difference_with_errbound(&right, &zero)
                        .unwrap()
                        .unwrap();
                    assert_abs_diff_eq!(difference_bound.estimate(), dot, epsilon = 0.0);
                    assert!(difference_bound.lower_bound() <= dot);
                    assert!(dot <= difference_bound.upper_bound());
                }

                #[test]
                fn [<vector_try_new_rejects_non_finite_ $d d>]() {
                    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
                        let mut data = [1.0f64; $d];
                        data[$d - 1] = value;
                        assert_eq!(
                            Vector::<$d>::try_new(data),
                            Err(LaError::non_finite_input_vector($d - 1))
                        );
                    }

                    let mut data = [1.0f64; $d];
                    data[0] = f64::INFINITY;
                    data[$d - 1] = f64::NAN;
                    assert_eq!(
                        Vector::<$d>::try_new(data),
                        Err(LaError::non_finite_input_vector(0))
                    );
                }

                #[test]
                fn [<vector_from_computation_preserves_failure_provenance_ $d d>]() {
                    let mut data = [1.0f64; $d];
                    data[$d - 1] = f64::INFINITY;

                    assert_eq!(
                        Vector::<$d>::from_computation(
                            data,
                            ArithmeticOperation::LuSolve,
                        ),
                        Err(LaError::non_finite_computation_step(
                            ArithmeticOperation::LuSolve,
                            $d - 1,
                        ))
                    );
                }

                #[test]
                fn [<vector_dot_and_norm_squared_reject_overflow_ $d d>]() {
                    let mut a_arr = [1.0f64; $d];
                    a_arr[0] = f64::MAX;
                    let a = Vector::<$d>::new(a_arr);

                    let mut b_arr = [1.0f64; $d];
                    b_arr[0] = 2.0;
                    let b = Vector::<$d>::new(b_arr);

                    assert_eq!(
                        a.dot(&b),
                        Err(LaError::non_finite_computation_step(
                            ArithmeticOperation::VectorDotProduct,
                            0,
                        ))
                    );
                    assert_eq!(
                        a.dot_with_errbound(&b),
                        Err(LaError::non_finite_computation_step(
                            ArithmeticOperation::VectorDotProduct,
                            0,
                        ))
                    );
                    assert_eq!(
                        a.dot_difference_with_errbound(&b, &Vector::zero()),
                        Err(LaError::non_finite_computation_step(
                            ArithmeticOperation::VectorDotDifference,
                            0,
                        ))
                    );
                    assert_eq!(
                        a.norm_squared(),
                        Err(LaError::non_finite_computation_step(
                            ArithmeticOperation::VectorSquaredNorm,
                            0,
                        ))
                    );
                }

            }
        };
    }

    // Mirror delaunay-style multi-dimension tests.
    gen_vector_tests!(1);
    gen_vector_tests!(2);
    gen_vector_tests!(3);
    gen_vector_tests!(4);
    gen_vector_tests!(5);
    gen_vector_tests!(6);
    gen_vector_tests!(7);
    gen_vector_tests!(8);

    fn known_norm_input<const D: usize>() -> ([f64; D], f64) {
        let mut data = [0.0; D];
        if D == 1 {
            data[0] = -5.0;
        } else if D >= 2 {
            data[0] = -3.0;
            data[1] = 4.0;
        }
        (data, if D == 0 { 0.0 } else { 5.0 })
    }

    macro_rules! gen_vector_norm_known_answer_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<vector_norm_known_answer_ $d d>]() {
                    let (data, expected) = known_norm_input::<$d>();
                    let vector = Vector::<$d>::new(data);

                    assert_eq!(vector.norm(), Ok(expected));
                }
            }
        };
    }

    gen_vector_norm_known_answer_tests!(0);
    gen_vector_norm_known_answer_tests!(1);
    gen_vector_norm_known_answer_tests!(2);
    gen_vector_norm_known_answer_tests!(3);
    gen_vector_norm_known_answer_tests!(4);
    gen_vector_norm_known_answer_tests!(5);
    gen_vector_norm_known_answer_tests!(6);
    gen_vector_norm_known_answer_tests!(7);
    gen_vector_norm_known_answer_tests!(8);

    macro_rules! gen_vector_replay_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<vector_dot_and_norm_squared_report_last_overflowing_step_ $d d>]() {
                    let mut dot_lhs = [1.0f64; $d];
                    dot_lhs[$d - 1] = f64::MAX;
                    let mut dot_rhs = [1.0f64; $d];
                    dot_rhs[$d - 1] = 2.0;
                    let dot_lhs = Vector::<$d>::new(dot_lhs);
                    let dot_rhs = Vector::<$d>::new(dot_rhs);

                    assert_eq!(
                        dot_lhs.dot(&dot_rhs),
                        Err(LaError::non_finite_computation_step(
                            ArithmeticOperation::VectorDotProduct,
                            $d - 1,
                        ))
                    );

                    let mut norm_data = [1.0f64; $d];
                    norm_data[$d - 1] = f64::MAX;
                    let vector = Vector::<$d>::new(norm_data);

                    assert_eq!(
                        vector.norm_squared(),
                        Err(LaError::non_finite_computation_step(
                            ArithmeticOperation::VectorSquaredNorm,
                            $d - 1,
                        ))
                    );
                }
            }
        };
    }

    gen_vector_replay_tests!(2);
    gen_vector_replay_tests!(3);
    gen_vector_replay_tests!(4);
    gen_vector_replay_tests!(5);

    macro_rules! gen_vector_const_eval_tests {
        ($d:literal, $dot:literal, $norm_squared:literal) => {
            paste! {
                #[test]
                fn [<vector_dot_and_norm_squared_const_eval_ $d d>]() {
                    const DOT: Result<f64, LaError> = Vector::<$d>::new([1.0; $d])
                        .dot(&Vector::<$d>::new([2.0; $d]));
                    const NORM_SQUARED: Result<f64, LaError> =
                        Vector::<$d>::new([1.0; $d]).norm_squared();

                    assert_eq!(DOT, Ok($dot));
                    assert_eq!(NORM_SQUARED, Ok($norm_squared));
                }
            }
        };
    }

    gen_vector_const_eval_tests!(2, 4.0, 2.0);
    gen_vector_const_eval_tests!(3, 6.0, 3.0);
    gen_vector_const_eval_tests!(4, 8.0, 4.0);
    gen_vector_const_eval_tests!(5, 10.0, 5.0);

    #[test]
    fn vector_dot_and_norm_squared_overflow_const_eval() {
        const DOT: Result<f64, LaError> =
            Vector::<2>::new([f64::MAX; 2]).dot(&Vector::<2>::new([1.0; 2]));
        const NORM_SQUARED: Result<f64, LaError> = Vector::<2>::new([f64::MAX; 2]).norm_squared();

        assert_eq!(
            DOT,
            Err(LaError::non_finite_computation_step(
                ArithmeticOperation::VectorDotProduct,
                1,
            ))
        );
        assert_eq!(
            NORM_SQUARED,
            Err(LaError::non_finite_computation_step(
                ArithmeticOperation::VectorSquaredNorm,
                0,
            ))
        );
    }

    #[test]
    fn vector_dot_and_norm_squared_preserve_fma_and_left_to_right_order() {
        let dot_large = 9_007_199_254_740_992.0;
        let dot_lhs = Vector::<4>::new([dot_large, 1.0, 1.0, 1.0]);
        let dot_rhs = Vector::<4>::new([1.0; 4]);
        assert_eq!(dot_lhs.dot(&dot_rhs), Ok(dot_large));

        let fused_lhs = Vector::<2>::new([f64::MAX, f64::MAX]);
        let fused_rhs = Vector::<2>::new([-1.0, 2.0]);
        assert_eq!(fused_lhs.dot(&fused_rhs), Ok(f64::MAX));

        let norm_large = 134_217_728.0;
        let vector = Vector::<4>::new([norm_large, 1.0, 1.0, 1.0]);
        assert_eq!(vector.norm_squared(), Ok(norm_large * norm_large));
    }

    #[test]
    fn vector_norm_preserves_zero_sign_and_subnormal_magnitudes() {
        let signed_zero = Vector::<4>::new([-0.0, 0.0, -0.0, 0.0]);
        assert_eq!(signed_zero.norm().unwrap().to_bits(), 0.0f64.to_bits());

        let least_subnormal = f64::from_bits(1);
        let subnormal = Vector::<2>::new([3.0 * least_subnormal, -4.0 * least_subnormal]);
        assert_eq!(
            subnormal.norm().unwrap().to_bits(),
            (5.0 * least_subnormal).to_bits()
        );
    }

    #[test]
    fn vector_norm_handles_mixed_and_overflowing_magnitudes() {
        let large = Vector::<2>::new([1.0e200, -1.0e200]);
        let expected = 2.0f64.sqrt() * 1.0e200;
        assert_abs_diff_eq!(large.norm().unwrap(), expected, epsilon = 2.0e184);
        assert_eq!(
            large.norm_squared(),
            Err(LaError::non_finite_computation_step(
                ArithmeticOperation::VectorSquaredNorm,
                0,
            )),
        );

        let mixed = Vector::<4>::new([1.0e200, 1.0e-200, -f64::from_bits(1), 0.0]);
        assert_eq!(mixed.norm(), Ok(1.0e200));

        let unrepresentable = Vector::<2>::new([f64::MAX, f64::MAX]);
        assert_eq!(
            unrepresentable.norm(),
            Err(LaError::non_finite_computation_scalar(
                ArithmeticOperation::VectorNorm,
            ))
        );
    }

    #[test]
    fn vector_norm_accepts_largest_finite_norm() {
        let maximum = Vector::<2>::new([f64::MAX, 0.0]);

        assert_eq!(maximum.norm(), Ok(f64::MAX));
    }

    #[test]
    fn certified_dot_preserves_fma_estimate_and_withholds_range_exhausted_bound() {
        let left = Vector::<2>::new([f64::MAX, f64::MAX]);
        let right = Vector::<2>::new([-1.0, 2.0]);

        assert_eq!(left.dot(&right), Ok(f64::MAX));
        assert_eq!(left.dot_with_errbound(&right), Ok(None));
    }

    #[test]
    fn certified_dot_withholds_bound_when_finite_endpoints_cannot_be_published() {
        // The magnitude bound reaches MAX but stays finite, so rejection must
        // come from publishing the estimate plus its positive error bound.
        let maximum = Vector::<2>::new([f64::MAX.next_down(), 1.0]);
        let one = Vector::<2>::new([1.0, 1.0]);

        assert_eq!(maximum.dot(&one), Ok(f64::MAX.next_down()));
        assert_eq!(maximum.dot_with_errbound(&one), Ok(None));
    }

    fn assert_one_exact_normal_product<const D: usize>() {
        let zero = Vector::<D>::new([0.0; D]);
        for value in [f64::MIN_POSITIVE, -2.0, f64::MAX] {
            let left =
                Vector::<D>::new(core::array::from_fn(|i| if i == 1 { value } else { -0.0 }));
            let right = Vector::<D>::new([1.0; D]);
            for (result, expected) in [
                (left.dot_with_errbound(&right), value),
                (left.dot_difference_with_errbound(&right, &zero), value),
                (left.dot_difference_with_errbound(&zero, &right), -value),
            ] {
                let certificate = result.unwrap().unwrap();
                assert_eq!(certificate.estimate().to_bits(), expected.to_bits());
                assert_eq!(certificate.absolute_error_bound(), 0.0);
                assert_eq!(certificate.lower_bound().to_bits(), expected.to_bits());
                assert_eq!(certificate.upper_bound().to_bits(), expected.to_bits());
            }
        }
        let inexact = Vector::<D>::new(core::array::from_fn(|i| if i == 1 { 0.1 } else { 0.0 }))
            .dot_with_errbound(&Vector::new([0.3; D]))
            .unwrap()
            .unwrap();
        assert!(inexact.absolute_error_bound() > 0.0);
        assert_eq!(
            Vector::<D>::new(core::array::from_fn(|i| if i == 1 {
                f64::MIN_POSITIVE
            } else {
                0.0
            }))
            .dot_with_errbound(&Vector::new([0.5; D])),
            Ok(None),
        );
    }

    macro_rules! gen_single_exact_product_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<single_exact_product_has_zero_error_ $d d>]() {
                    assert_one_exact_normal_product::<$d>();
                }
            }
        };
    }

    gen_single_exact_product_tests!(2);
    gen_single_exact_product_tests!(3);
    gen_single_exact_product_tests!(4);
    gen_single_exact_product_tests!(5);
    gen_single_exact_product_tests!(6);

    #[test]
    fn scalar_certificate_checks_both_endpoint_ranges_before_construction() {
        const ENDPOINTS: Option<(f64, f64)> =
            match ScalarWithErrorBound::try_new(f64::MAX / 2.0, f64::MAX / 2.0) {
                Some(value) => Some((value.lower_bound(), value.upper_bound())),
                None => None,
            };
        let previous = f64::MAX.next_down();
        let gap = f64::MAX - previous;
        for sign in [-1.0, 1.0] {
            assert_eq!(
                ScalarWithErrorBound::try_new(sign * f64::MAX, f64::from_bits(1)),
                None,
            );
            assert_eq!(
                ScalarWithErrorBound::try_new(sign * previous, gap.next_up()),
                None,
            );
            for bound in [gap.next_down(), gap] {
                let certificate = ScalarWithErrorBound::try_new(sign * previous, bound)
                    .expect("both exact endpoints remain within the finite range");
                let (lower, upper) = if sign > 0.0 {
                    (previous.next_down(), f64::MAX)
                } else {
                    (-f64::MAX, -previous.next_down())
                };
                assert_eq!(certificate.lower_bound().to_bits(), lower.to_bits());
                assert_eq!(certificate.upper_bound().to_bits(), upper.to_bits());
            }
        }
        let symmetric = ScalarWithErrorBound::try_new(0.0, f64::MAX).unwrap();
        assert_eq!(symmetric.lower_bound(), -f64::MAX);
        assert_eq!(symmetric.upper_bound(), f64::MAX);
        let negative_zero = ScalarWithErrorBound::try_new(-0.0, -0.0).unwrap();
        assert_eq!(negative_zero.lower_bound().to_bits(), (-0.0_f64).to_bits());
        assert_eq!(negative_zero.upper_bound().to_bits(), (-0.0_f64).to_bits());

        assert_eq!(ENDPOINTS, Some((0.0, f64::MAX)));
    }

    #[test]
    fn certified_dot_withholds_bound_when_magnitude_sum_exhausts_range() {
        let maximum = Vector::<2>::new([f64::MAX, f64::MAX]);

        for factor in [0.5, 0.75] {
            let cancelling = Vector::<2>::new([factor, -factor]);
            let estimate = maximum
                .dot(&cancelling)
                .expect("the cancelling FMA estimate must remain finite");
            assert!(estimate.is_finite());
            assert_eq!(maximum.dot_with_errbound(&cancelling), Ok(None));
        }
    }

    #[test]
    fn certified_bounds_distinguish_conclusive_and_inconclusive_results() {
        let conclusive = Vector::<2>::new([1.0, 2.0])
            .dot_with_errbound(&Vector::new([3.0, 4.0]))
            .unwrap()
            .unwrap();
        assert_abs_diff_eq!(conclusive.estimate(), 11.0, epsilon = 0.0);
        assert!(conclusive.lower_bound() > 1.0);

        let inconclusive = Vector::<2>::new([1.0, 1.0])
            .dot_with_errbound(&Vector::new([1.0, -1.0]))
            .unwrap()
            .unwrap();
        assert_abs_diff_eq!(inconclusive.estimate(), 0.0, epsilon = 0.0);
        assert!(inconclusive.absolute_error_bound() > 0.0);
        assert!(inconclusive.lower_bound() < 0.0);
        assert!(inconclusive.upper_bound() > 0.0);
    }

    #[test]
    fn certified_dot_withholds_bound_when_finite_magnitude_sum_cannot_be_widened() {
        let large = Vector::<2>::new([f64::MAX / 2.0, f64::MAX / 2.0]);
        let cancelling = Vector::<2>::new([1.0, -1.0]);

        // The magnitude FMA reaches finite MAX; its upward successor is
        // infinite. Exact estimate cancellation cannot bypass proof loss.
        assert_eq!(large.dot(&cancelling), Ok(0.0));
        assert_eq!(large.dot_with_errbound(&cancelling), Ok(None));
    }

    #[test]
    fn certified_proof_loss_survives_repeated_underflow_and_a_normal_final_estimate() {
        let axis = Vector::<3>::new([f64::from_bits(1), f64::from_bits(1), 1.0]);
        let values = Vector::<3>::new([0.5, 0.5, 2.0]);
        let zero = Vector::<3>::zero();

        // Both tiny products round to zero. The final term produces a normal
        // estimate, but cannot restore a proof lost before that term.
        assert_eq!(axis.dot(&values), Ok(2.0));
        assert_eq!(axis.dot_with_errbound(&values), Ok(None));
        assert_eq!(axis.dot_difference_with_errbound(&values, &zero), Ok(None));
    }

    #[test]
    fn zero_products_preserve_proof_loss_and_signed_zero_estimates() {
        let operation = ArithmeticOperation::VectorDotProduct;
        let underflow = CertifiedReduction::ZERO
            .add_product(-f64::from_bits(1), 0.5, operation, 0)
            .unwrap();
        assert_eq!(underflow.estimate.to_bits(), (-0.0_f64).to_bits());
        let after_zero = underflow.add_product(0.0, 1.0, operation, 1).unwrap();
        assert_eq!(after_zero.estimate.to_bits(), 0.0_f64.to_bits());
        assert!(!after_zero.proof_available);
        assert_eq!(after_zero.finish(Some(2)), None);

        let normal = CertifiedReduction::ZERO
            .add_product(1.0, 2.0, operation, 0)
            .unwrap();
        let cancelled = normal.add_product(-1.0, 2.0, operation, 1).unwrap();
        assert_eq!(cancelled.estimate.to_bits(), 0.0_f64.to_bits());
        for state in [CertifiedReduction::ZERO, normal, cancelled] {
            for zero in [0.0, -0.0] {
                for other in [f64::MAX, -f64::MIN_POSITIVE] {
                    let next = state.add_product(zero, other, operation, 2).unwrap();
                    assert_eq!(next, state);
                    assert_eq!(
                        next.estimate.to_bits(),
                        zero.mul_add(other, state.estimate).to_bits()
                    );
                }
            }
        }
    }

    #[test]
    fn certified_zero_and_signed_zero_have_an_exact_zero_bound() {
        let left = Vector::<3>::new([-0.0, 0.0, -0.0]);
        let right = Vector::<3>::new([f64::MAX, -1.0, f64::MIN_POSITIVE]);
        let bounded = left.dot_with_errbound(&right).unwrap().unwrap();

        assert_eq!(bounded.estimate().to_bits(), 0.0_f64.to_bits());
        assert_eq!(bounded.absolute_error_bound().to_bits(), 0.0_f64.to_bits());
        assert_eq!(bounded.lower_bound().to_bits(), 0.0_f64.to_bits());
        assert_eq!(bounded.upper_bound().to_bits(), 0.0_f64.to_bits());
    }

    #[test]
    fn certified_reductions_withhold_bounds_for_subnormal_products() {
        let tiny = Vector::<1>::new([f64::MIN_POSITIVE]);
        let half = Vector::<1>::new([0.5]);
        assert_eq!(tiny.dot_with_errbound(&half), Ok(None));

        let min_subnormal = Vector::<1>::new([f64::from_bits(1)]);
        assert_eq!(
            min_subnormal.dot_with_errbound(&Vector::new([1.0])),
            Ok(None)
        );
        assert_eq!(
            min_subnormal.dot_with_errbound(&Vector::new([0.5])),
            Ok(None)
        );
        assert_eq!(
            tiny.dot_difference_with_errbound(&half, &Vector::zero()),
            Ok(None)
        );
    }

    #[test]
    fn certified_reduction_detects_fma_cancellation_below_subnormal_range() {
        let near_sqrt_min = f64::from_bits((512_u64 << 52) | 1);
        let rounded_product = near_sqrt_min * near_sqrt_min;
        assert!(rounded_product.is_normal());
        assert_abs_diff_eq!(
            near_sqrt_min.mul_add(near_sqrt_min, -rounded_product),
            0.0,
            epsilon = 0.0
        );

        let left = Vector::<2>::new([-rounded_product, near_sqrt_min]);
        let right = Vector::<2>::new([1.0, near_sqrt_min]);
        assert_eq!(left.dot(&right), Ok(0.0));
        assert_eq!(left.dot_with_errbound(&right), Ok(None));
    }

    #[test]
    fn certified_dot_handles_mixed_normal_magnitudes() {
        let left = Vector::<2>::new([1.0e100, 1.0e-100]);
        let right = Vector::<2>::new([1.0e-100, 1.0e100]);
        let bounded = left.dot_with_errbound(&right).unwrap().unwrap();

        assert_abs_diff_eq!(bounded.estimate(), 2.0, epsilon = 0.0);
        assert!(bounded.lower_bound() <= 2.0);
        assert!(bounded.upper_bound() >= 2.0);
    }

    #[test]
    fn certified_dot_difference_does_not_round_coordinates_first() {
        let scale = 18_014_398_509_481_984.0;
        let axis = Vector::<1>::new([scale]);
        let left = Vector::<1>::new([1.0]);
        let right = Vector::<1>::new([1.0 / scale]);
        assert_abs_diff_eq!(left.as_array()[0] - right.as_array()[0], 1.0, epsilon = 0.0);

        let bounded = axis
            .dot_difference_with_errbound(&left, &right)
            .unwrap()
            .unwrap();
        assert_abs_diff_eq!(bounded.estimate(), scale, epsilon = 0.0);
        assert!(bounded.lower_bound() < scale);
        assert!(bounded.upper_bound() >= scale);
    }

    #[test]
    fn certified_reductions_are_const_evaluable() {
        const DOT: Result<Option<ScalarWithErrorBound>, LaError> =
            Vector::<2>::new([1.0, 2.0]).dot_with_errbound(&Vector::<2>::new([3.0, 4.0]));
        const DIFFERENCE: Result<Option<ScalarWithErrorBound>, LaError> =
            Vector::<2>::new([2.0, -1.0]).dot_difference_with_errbound(
                &Vector::<2>::new([4.0, 1.0]),
                &Vector::<2>::new([1.0, 3.0]),
            );

        assert_abs_diff_eq!(DOT.unwrap().unwrap().estimate(), 11.0, epsilon = 0.0);
        assert_abs_diff_eq!(DIFFERENCE.unwrap().unwrap().estimate(), 8.0, epsilon = 0.0);
    }

    #[test]
    fn certified_dot_difference_reports_second_fma_overflow() {
        let axis = Vector::<2>::new([1.0, f64::MAX]);
        let left = Vector::<2>::new([0.0, 1.0]);
        let right = Vector::<2>::new([0.0, -1.0]);

        assert_eq!(
            axis.dot_difference_with_errbound(&left, &right),
            Err(LaError::non_finite_computation_step(
                ArithmeticOperation::VectorDotDifference,
                1,
            ))
        );
    }

    #[test]
    fn vector_dot_and_norm_squared_report_first_middle_overflowing_step() {
        let dot_lhs = Vector::<3>::new([f64::MAX, f64::MAX, 1.0]);
        let dot_rhs = Vector::<3>::new([1.0; 3]);
        assert_eq!(
            dot_lhs.dot(&dot_rhs),
            Err(LaError::non_finite_computation_step(
                ArithmeticOperation::VectorDotProduct,
                1,
            ))
        );

        let norm_large = 1.0e154;
        let vector = Vector::<3>::new([norm_large, norm_large, 1.0]);
        assert_eq!(
            vector.norm_squared(),
            Err(LaError::non_finite_computation_step(
                ArithmeticOperation::VectorSquaredNorm,
                1,
            ))
        );
    }

    #[test]
    fn zero_dimension_vector_has_zero_dot_and_norm() {
        let vector = Vector::<0>::try_new([]).unwrap();

        assert_eq!(vector.as_array().as_slice(), []);
        assert_eq!(vector.into_array().as_slice(), []);
        assert_eq!(vector.dot(&Vector::zero()), Ok(0.0));
        let dot_bound = vector.dot_with_errbound(&Vector::zero()).unwrap().unwrap();
        assert_abs_diff_eq!(dot_bound.absolute_error_bound(), 0.0, epsilon = 0.0);
        let difference_bound = vector
            .dot_difference_with_errbound(&Vector::zero(), &Vector::zero())
            .unwrap()
            .unwrap();
        assert_abs_diff_eq!(difference_bound.absolute_error_bound(), 0.0, epsilon = 0.0);
        assert_eq!(vector.norm_squared(), Ok(0.0));
    }

    #[test]
    fn certified_dot_withholds_bound_when_exact_product_exceeds_finite_range() {
        let left_value = f64::from_bits(0x7fe3_0319_b612_3729);
        let right_value = f64::from_bits(0x3ffa_ee21_bf46_bc00);
        let left = Vector::<1>::new([left_value]);
        let right = Vector::<1>::new([right_value]);

        assert!(left_value.mul_add(right_value, -f64::MAX) > 0.0);
        assert_eq!(left.dot(&right), Ok(f64::MAX));
        assert_eq!(left.dot_with_errbound(&right), Ok(None));
    }
}
