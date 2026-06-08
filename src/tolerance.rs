//! Finite tolerance values used by numerical predicates and factorizations.

use crate::LaError;

/// Finite, non-negative tolerance used by numerical predicates and factorizations.
///
/// Construct with [`Tolerance::new`] when accepting raw caller input. Once
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
    /// Construct a tolerance without checking the raw value.
    ///
    /// This crate-internal escape hatch is only for constants whose finite,
    /// non-negative value is visible at the call site. Public callers should
    /// use [`Tolerance::new`] so the returned value carries the validation
    /// proof.
    pub(crate) const fn new_unchecked(value: f64) -> Self {
        Self { value }
    }

    /// Construct a finite, non-negative tolerance.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let tol = Tolerance::new(1e-12)?;
    /// assert_eq!(tol.get(), 1e-12);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::InvalidTolerance`] when `value` is NaN, infinite, or
    /// negative.
    #[inline]
    pub const fn new(value: f64) -> Result<Self, LaError> {
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
    /// let tol = Tolerance::new(0.0)?;
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

/// Relative tolerance used to validate matrices for LDLT factorization.
///
/// This is crate-internal because LDLT callers provide the factorization
/// tolerance separately; the symmetry tolerance is a fixed domain check used to
/// parse a public [`Matrix`](crate::Matrix) into the internal symmetric proof
/// type before factorization.
pub const LDLT_SYMMETRY_REL_TOL: Tolerance = Tolerance::new_unchecked(1e-12);

#[cfg(test)]
mod tests {
    use core::assert_matches;

    use approx::assert_abs_diff_eq;

    use super::*;

    #[test]
    fn default_singular_tol_is_expected() {
        assert_abs_diff_eq!(DEFAULT_SINGULAR_TOL.get(), 1e-12, epsilon = 0.0);
    }

    #[test]
    fn tolerance_new_accepts_finite_non_negative_values() {
        assert_eq!(
            Tolerance::new(0.0).unwrap().get().to_bits(),
            0.0f64.to_bits()
        );
        assert_eq!(
            Tolerance::new(1e-12).unwrap().get().to_bits(),
            1e-12f64.to_bits()
        );
        assert_eq!(
            Tolerance::new(f64::MAX).unwrap().get().to_bits(),
            f64::MAX.to_bits()
        );
    }

    #[test]
    fn tolerance_new_rejects_negative_nan_and_infinity() {
        assert_eq!(
            Tolerance::new(-1.0),
            Err(LaError::InvalidTolerance { value: -1.0 })
        );
        assert_matches!(
            Tolerance::new(f64::NAN),
            Err(LaError::InvalidTolerance { value }) if value.is_nan()
        );
        assert_eq!(
            Tolerance::new(f64::INFINITY),
            Err(LaError::InvalidTolerance {
                value: f64::INFINITY,
            })
        );
        assert_eq!(
            Tolerance::new(f64::NEG_INFINITY),
            Err(LaError::InvalidTolerance {
                value: f64::NEG_INFINITY,
            })
        );
    }
}
