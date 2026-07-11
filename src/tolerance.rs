#![forbid(unsafe_code)]

//! Finite tolerance values used by numerical predicates and factorizations.

use crate::LaError;

/// Finite, non-negative tolerance used by numerical predicates and factorizations.
///
/// Construct with [`Tolerance::try_new`] when accepting raw caller input. Once
/// constructed, the stored value is guaranteed to be finite and `>= 0`, so
/// downstream algorithms do not need to revalidate the tolerance.
///
/// This is the crate-wide tolerance contract: raw negative, NaN, and infinite
/// values are rejected with [`LaError::InvalidTolerance`] at construction time.
#[must_use]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Tolerance {
    value: f64,
}

impl Tolerance {
    /// Construct a tolerance for a finite, non-negative module-local literal.
    const fn new_unchecked(value: f64) -> Self {
        Self { value }
    }

    /// Exact zero tolerance for crate-internal algorithms.
    pub(crate) const ZERO: Self = Self::new_unchecked(0.0);

    /// Construct a finite, non-negative tolerance.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let tol = Tolerance::try_new(1e-12)?;
    /// assert_eq!(tol.get(), 1e-12);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::InvalidTolerance`] with
    /// [`crate::InvalidToleranceReason::NotFinite`] for NaN/infinity or
    /// [`crate::InvalidToleranceReason::Negative`] for a finite negative
    /// value. Both signed-zero representations are accepted and preserved.
    #[inline]
    pub const fn try_new(value: f64) -> Result<Self, LaError> {
        if value >= 0.0 && value.is_finite() {
            Ok(Self::new_unchecked(value))
        } else {
            Err(LaError::invalid_tolerance(value))
        }
    }

    /// Return the raw finite, non-negative tolerance value.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let tol = Tolerance::try_new(0.0)?;
    /// assert_eq!(tol.get(), 0.0);
    /// # Ok(())
    /// # }
    /// ```
    #[inline]
    #[must_use]
    pub const fn get(self) -> f64 {
        self.value
    }
}

/// Default absolute threshold used for singularity/degeneracy detection.
///
/// This is intentionally conservative for geometric predicates and small systems.
///
/// Conceptually, this is an absolute bound for deciding when a scalar should be treated
/// as "numerically zero" (e.g. LU pivots, LDLT diagonal entries).
///
/// # Examples
/// ```
/// use la_stack::prelude::*;
///
/// # fn main() -> Result<(), LaError> {
/// let lu = Matrix::<2>::identity().lu(DEFAULT_SINGULAR_TOL)?;
/// assert_eq!(lu.det()?, 1.0);
/// # Ok(())
/// # }
/// ```
pub const DEFAULT_SINGULAR_TOL: Tolerance = Tolerance::new_unchecked(1e-12);

#[cfg(test)]
mod tests {
    use core::assert_matches;

    use approx::assert_abs_diff_eq;

    use super::*;

    #[test]
    fn default_singular_tol_is_expected() {
        assert_abs_diff_eq!(DEFAULT_SINGULAR_TOL.get(), 1e-12, epsilon = 0.0);
        assert_eq!(Tolerance::ZERO.get().to_bits(), 0.0f64.to_bits());
    }

    #[test]
    fn try_new_accepts_finite_non_negative_values() {
        assert_eq!(
            Tolerance::try_new(0.0).unwrap().get().to_bits(),
            0.0f64.to_bits()
        );
        assert_eq!(
            Tolerance::try_new(1e-12).unwrap().get().to_bits(),
            1e-12f64.to_bits()
        );
        assert_eq!(
            Tolerance::try_new(f64::MAX).unwrap().get().to_bits(),
            f64::MAX.to_bits()
        );
    }

    #[test]
    fn try_new_accepts_and_preserves_negative_zero() {
        let tolerance = Tolerance::try_new(-0.0).unwrap();

        assert_eq!(tolerance.get().to_bits(), (-0.0f64).to_bits());
    }

    #[test]
    fn try_new_rejects_negative_finite_values() {
        assert_eq!(
            Tolerance::try_new(-1.0),
            Err(LaError::InvalidTolerance {
                value: -1.0,
                reason: crate::InvalidToleranceReason::Negative,
            })
        );
    }

    #[test]
    fn try_new_rejects_non_finite_values_with_structured_reason() {
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_matches!(
                Tolerance::try_new(value),
                Err(LaError::InvalidTolerance {
                    value: observed,
                    reason: crate::InvalidToleranceReason::NotFinite,
                }) if observed.to_bits() == value.to_bits()
            );
        }
    }
}
