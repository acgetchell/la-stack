#![forbid(unsafe_code)]

//! Allocation-free scaled products for floating-point factor diagonals.

const SIGN_MASK: u64 = 1_u64 << 63;
const FRACTION_BITS: u32 = 52;
const FRACTION_MASK: u64 = (1_u64 << FRACTION_BITS) - 1;
const EXPONENT_MASK: u64 = 0x7ff;
const EXPONENT_BIAS: i128 = 1023;
const MIN_NORMAL_EXPONENT: i128 = -1022;
const MIN_SUBNORMAL_EXPONENT: i128 = -1074;

/// Result of multiplying one direct product step with its range proof attached.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum RangeCheckedProduct {
    /// The value is finite and normal, or is an exact zero caused by a zero
    /// operand, so direct accumulation may continue.
    Safe(f64),
    /// Direct multiplication overflowed or lost range through gradual
    /// underflow, so all factors must be recomputed with scaling.
    NeedsScaling,
}

/// Multiply one direct product step and retain only values proven safe for
/// sequential accumulation.
///
/// The common normal case needs only the result's binary exponent field. Zero
/// operands are checked only when that field is zero, preserving signed-zero
/// multiplication without treating underflow from two non-zero operands as an
/// exact zero.
#[inline]
pub(crate) const fn range_checked_product(accumulator: f64, factor: f64) -> RangeCheckedProduct {
    let product = accumulator * factor;
    let product_exponent = (product.to_bits() >> FRACTION_BITS) & EXPONENT_MASK;
    // Subtracting one maps the valid normal fields 1..=0x7fe to
    // 0..=0x7fd. Zero wraps high and 0x7ff maps to the exclusive upper
    // bound, so the common normal path needs one unsigned comparison.
    if product_exponent.wrapping_sub(1) < EXPONENT_MASK - 1 {
        return RangeCheckedProduct::Safe(product);
    }

    if product_exponent == 0 && (accumulator == 0.0 || factor == 0.0) {
        RangeCheckedProduct::Safe(product)
    } else {
        RangeCheckedProduct::NeedsScaling
    }
}

/// One non-zero finite factor normalized as `mantissa × 2^exponent`.
#[derive(Clone, Copy)]
struct NormalizedFactor {
    mantissa: f64,
    exponent: i128,
}

/// A finite product kept as `(-1)^negative × mantissa × 2^exponent`.
///
/// Non-zero finite factors are normalized to `1 ≤ mantissa < 2` before
/// multiplication. Consequently, no intermediate mantissa multiplication can
/// underflow or overflow; only [`Self::finish`] decides whether the final
/// rounded result is finite. The most recent factor stays deferred so a final
/// subnormal product can be formed directly in the binary64 destination range
/// without first rounding it as a normal mantissa.
pub(crate) struct ScaledProduct {
    mantissa: f64,
    exponent: i128,
    pending_factor: Option<NormalizedFactor>,
    negative: bool,
    zero: bool,
    non_finite: bool,
}

impl ScaledProduct {
    /// Start an empty product with the requested initial sign.
    #[inline]
    pub(crate) const fn new(negative: bool) -> Self {
        Self {
            mantissa: 1.0,
            exponent: 0,
            pending_factor: None,
            negative,
            zero: false,
            non_finite: false,
        }
    }

    /// Multiply by one factor while retaining a normalized mantissa.
    #[inline]
    pub(crate) const fn multiply(&mut self, factor: f64) {
        let bits = factor.to_bits();
        self.negative ^= bits & SIGN_MASK != 0;

        let magnitude = bits & !SIGN_MASK;
        let biased_exponent = (magnitude >> FRACTION_BITS) & EXPONENT_MASK;
        let fraction = magnitude & FRACTION_MASK;

        if biased_exponent == EXPONENT_MASK {
            self.non_finite = true;
            return;
        }
        if biased_exponent == 0 && fraction == 0 {
            self.zero = true;
            return;
        }
        if self.zero {
            return;
        }

        let (factor_mantissa, factor_exponent) = if biased_exponent == 0 {
            // A subnormal value is `fraction × 2^-1074`. Move its highest set
            // bit to the binary64 hidden-bit position to obtain a mantissa in
            // [1, 2), and compensate in the exponent.
            let highest_bit = fraction.ilog2();
            let shift = FRACTION_BITS - highest_bit;
            let significand = fraction << shift;
            (
                f64::from_bits((1023_u64 << FRACTION_BITS) | (significand & FRACTION_MASK)),
                (highest_bit as i128) + MIN_SUBNORMAL_EXPONENT,
            )
        } else {
            (
                f64::from_bits((1023_u64 << FRACTION_BITS) | fraction),
                (biased_exponent as i128) - EXPONENT_BIAS,
            )
        };

        if let Some(pending) = self.pending_factor {
            self.absorb_factor(pending);
        }
        self.pending_factor = Some(NormalizedFactor {
            mantissa: factor_mantissa,
            exponent: factor_exponent,
        });
    }

    /// Fold one normalized factor into the running product.
    #[inline]
    const fn absorb_factor(&mut self, factor: NormalizedFactor) {
        self.mantissa *= factor.mantissa;
        self.exponent += factor.exponent;

        // Both operands were in [1, 2), and even the two largest binary64
        // mantissas multiply to a value that rounds below 4.0. Therefore one
        // factor-of-two normalization is sufficient to preserve
        // `mantissa < 2`.
        if self.mantissa >= 2.0 {
            self.mantissa *= 0.5;
            self.exponent += 1;
        }
    }

    /// Round the accumulated product to binary64.
    ///
    /// Returns `None` only when the accumulated result rounds outside the
    /// finite binary64 range. Magnitudes below that range round to a signed
    /// zero or subnormal value with round-to-nearest, ties-to-even semantics.
    #[inline]
    #[expect(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "the preceding bounds prove the normal and deferred-final biased exponents fit u64"
    )]
    pub(crate) const fn finish(mut self) -> Option<f64> {
        if self.non_finite {
            return None;
        }

        let sign = if self.negative { SIGN_MASK } else { 0 };
        if self.zero {
            return Some(f64::from_bits(sign));
        }
        let Some(pending) = self.pending_factor else {
            return Some(f64::from_bits(sign | (1023_u64 << FRACTION_BITS)));
        };

        let final_exponent = self.exponent + pending.exponent;

        // The product of two mantissas in [1, 2) is strictly below 4. Values
        // below this exponent are therefore strictly below half the least
        // subnormal and round to signed zero.
        if final_exponent < MIN_SUBNORMAL_EXPONENT - 2 {
            return Some(f64::from_bits(sign));
        }

        if final_exponent < MIN_NORMAL_EXPONENT {
            // Scale both normal operands so their single multiplication lands
            // directly in the final subnormal range. Multiplying normalized
            // mantissas first would round once to 53 bits here and a second time
            // to the much coarser subnormal grid below.
            let left = f64::from_bits(
                (1_u64 << FRACTION_BITS) | (self.mantissa.to_bits() & FRACTION_MASK),
            );
            let right_biased_exponent =
                (final_exponent - MIN_NORMAL_EXPONENT + EXPONENT_BIAS) as u64;
            let right = f64::from_bits(
                (right_biased_exponent << FRACTION_BITS)
                    | (pending.mantissa.to_bits() & FRACTION_MASK),
            );
            let magnitude = left * right;
            return Some(f64::from_bits(sign | magnitude.to_bits()));
        }

        self.absorb_factor(pending);
        if self.exponent > EXPONENT_BIAS {
            return None;
        }

        let mantissa_bits = self.mantissa.to_bits();
        let fraction = mantissa_bits & FRACTION_MASK;
        let biased_exponent = (self.exponent + EXPONENT_BIAS) as u64;
        Some(f64::from_bits(
            sign | (biased_exponent << FRACTION_BITS) | fraction,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::{RangeCheckedProduct, SIGN_MASK, ScaledProduct, range_checked_product};

    const TWO_NEG_800: f64 = f64::from_bits(223_u64 << 52);
    const TWO_POS_800: f64 = f64::from_bits(1823_u64 << 52);
    const SUBNORMAL_ROUNDING_LEFT: f64 = f64::from_bits(0x3cb2_e219_27ac_435a);
    const SUBNORMAL_ROUNDING_RIGHT: f64 = f64::from_bits(0x0014_55e5_f80b_50eb);

    /// Return the exact bits produced by a two-factor scaled product.
    fn scaled_product_bits(left: f64, right: f64) -> Option<u64> {
        let mut product = ScaledProduct::new(false);
        product.multiply(left);
        product.multiply(right);
        product.finish().map(f64::to_bits)
    }

    #[test]
    fn mantissa_product_is_renormalized_once() {
        assert_eq!(scaled_product_bits(1.5, 1.5), Some(2.25_f64.to_bits()));
    }

    #[test]
    fn signed_zero_tracks_later_factor_signs() {
        let mut product = ScaledProduct::new(false);
        product.multiply(-0.0);
        product.multiply(-2.0);

        assert_eq!(product.finish().map(f64::to_bits), Some(0.0_f64.to_bits()));
    }

    #[test]
    fn empty_product_preserves_initial_sign() {
        assert_eq!(
            ScaledProduct::new(false).finish().map(f64::to_bits),
            Some(1.0_f64.to_bits())
        );
        assert_eq!(
            ScaledProduct::new(true).finish().map(f64::to_bits),
            Some((-1.0_f64).to_bits())
        );
    }

    #[test]
    fn non_finite_factors_make_the_product_unrepresentable() {
        for factor in [f64::INFINITY, f64::NEG_INFINITY, f64::NAN] {
            let mut product = ScaledProduct::new(false);
            product.multiply(factor);

            assert_eq!(product.finish(), None);
        }
    }

    #[test]
    fn balanced_extreme_factors_do_not_depend_on_storage_order() {
        let mut forward = ScaledProduct::new(false);
        for factor in [TWO_NEG_800, TWO_NEG_800, TWO_POS_800, TWO_POS_800] {
            forward.multiply(factor);
        }

        let mut reverse = ScaledProduct::new(false);
        for factor in [TWO_POS_800, TWO_POS_800, TWO_NEG_800, TWO_NEG_800] {
            reverse.multiply(factor);
        }

        assert_eq!(forward.finish(), Some(1.0));
        assert_eq!(reverse.finish(), Some(1.0));
    }

    #[test]
    fn final_range_decision_distinguishes_underflow_and_overflow() {
        let mut underflow = ScaledProduct::new(true);
        underflow.multiply(TWO_NEG_800);
        underflow.multiply(TWO_NEG_800);
        assert_eq!(underflow.finish().map(f64::to_bits), Some(1_u64 << 63));

        let mut overflow = ScaledProduct::new(false);
        overflow.multiply(TWO_POS_800);
        overflow.multiply(TWO_POS_800);
        assert_eq!(overflow.finish(), None);
    }

    #[test]
    fn final_subnormal_product_is_rounded_once_in_const_evaluation() {
        const POSITIVE: Option<f64> = {
            let mut product = ScaledProduct::new(false);
            product.multiply(SUBNORMAL_ROUNDING_LEFT);
            product.multiply(SUBNORMAL_ROUNDING_RIGHT);
            product.finish()
        };
        const NEGATIVE: Option<f64> = {
            let mut product = ScaledProduct::new(true);
            product.multiply(SUBNORMAL_ROUNDING_LEFT);
            product.multiply(SUBNORMAL_ROUNDING_RIGHT);
            product.finish()
        };

        assert_eq!(
            (SUBNORMAL_ROUNDING_LEFT * SUBNORMAL_ROUNDING_RIGHT).to_bits(),
            1
        );
        assert_eq!(POSITIVE.map(f64::to_bits), Some(1));
        assert_eq!(NEGATIVE.map(f64::to_bits), Some(SIGN_MASK | 1));
    }

    #[test]
    fn final_subnormal_product_is_rounded_once_after_earlier_range_loss() {
        let factors = [
            TWO_NEG_800,
            TWO_NEG_800,
            TWO_POS_800,
            TWO_POS_800,
            SUBNORMAL_ROUNDING_LEFT,
            SUBNORMAL_ROUNDING_RIGHT,
        ];
        let mut product = ScaledProduct::new(false);
        for factor in factors {
            product.multiply(factor);
        }

        assert_eq!(
            range_checked_product(TWO_NEG_800, TWO_NEG_800),
            RangeCheckedProduct::NeedsScaling
        );
        assert_eq!(product.finish().map(f64::to_bits), Some(1));
    }

    #[test]
    fn final_subnormal_product_preserves_ties_to_even_at_range_boundaries() {
        let least_subnormal = f64::from_bits(1);
        let three_subnormals = f64::from_bits(3);
        let largest_below_one = f64::from_bits(0x3fef_ffff_ffff_ffff);

        assert_eq!(scaled_product_bits(least_subnormal, 0.5), Some(0));
        assert_eq!(scaled_product_bits(three_subnormals, 0.5), Some(2));
        assert_eq!(
            scaled_product_bits(f64::MIN_POSITIVE, largest_below_one),
            Some(f64::MIN_POSITIVE.to_bits())
        );
    }

    #[test]
    fn direct_product_range_check_distinguishes_exact_zero_from_range_loss() {
        assert_eq!(
            range_checked_product(1.5, 2.0),
            RangeCheckedProduct::Safe(3.0)
        );
        assert_eq!(
            range_checked_product(-0.0, -2.0),
            RangeCheckedProduct::Safe(0.0)
        );
        assert_eq!(
            range_checked_product(TWO_NEG_800, TWO_NEG_800),
            RangeCheckedProduct::NeedsScaling
        );
        assert_eq!(
            range_checked_product(f64::MIN_POSITIVE, 0.5),
            RangeCheckedProduct::NeedsScaling
        );
        assert_eq!(
            range_checked_product(TWO_POS_800, TWO_POS_800),
            RangeCheckedProduct::NeedsScaling
        );
    }
}
