#![forbid(unsafe_code)]

//! Shared binary64 rounding primitives for certified arithmetic.
//!
//! Representation and rounding use the IEEE 754 model in `REFERENCES.md`
//! \[9-10\]. The interval and reduction certificates share these primitives.

/// Return the exact error in a rounded binary64 sum.
///
/// This is `FastTwoSum` with operands ordered by magnitude. With IEEE-754
/// round-to-nearest and gradual underflow, `rounded + error` equals the
/// exact-real sum whenever the rounded sum is finite. Ordering prevents an
/// intermediate overflow even at the finite range boundary; see
/// `REFERENCES.md` \[17\], Theorem 5.1. Callers supply finite operands and their
/// rounded sum.
#[inline]
pub(crate) const fn two_sum_error(left: f64, right: f64, rounded: f64) -> f64 {
    let (large, small) = if left.abs() >= right.abs() {
        (left, right)
    } else {
        (right, left)
    };
    let virtual_small = rounded - large;
    small - virtual_small
}

/// Decompose a nonzero finite binary64 magnitude as `significand × 2^exponent`.
#[inline]
pub(crate) const fn decompose_magnitude(value: f64) -> (u128, i64) {
    let magnitude_bits = value.to_bits() & 0x7fff_ffff_ffff_ffff;
    let biased_exponent = ((magnitude_bits >> 52) & 0x7ff).cast_signed();
    let fraction = magnitude_bits & 0x000f_ffff_ffff_ffff;

    if biased_exponent == 0 {
        (fraction as u128, -1074)
    } else {
        (
            (fraction | (1_u64 << 52)) as u128,
            biased_exponent - 1023 - 52,
        )
    }
}

/// Compare two positive values represented as `significand × 2^exponent`.
#[inline]
const fn compare_binary_magnitudes(
    left_significand: u128,
    left_exponent: i64,
    right_significand: u128,
    right_exponent: i64,
) -> i8 {
    let left_zeros = left_significand.trailing_zeros() as i64;
    let right_zeros = right_significand.trailing_zeros() as i64;
    let normalized_left = left_significand >> left_zeros.cast_unsigned();
    let normalized_right = right_significand >> right_zeros.cast_unsigned();
    let normalized_left_exponent = left_exponent + left_zeros;
    let normalized_right_exponent = right_exponent + right_zeros;

    let left_top = normalized_left_exponent + (normalized_left.bit_width() - 1) as i64;
    let right_top = normalized_right_exponent + (normalized_right.bit_width() - 1) as i64;
    if left_top < right_top {
        return -1;
    }
    if left_top > right_top {
        return 1;
    }

    let common_exponent = if normalized_left_exponent < normalized_right_exponent {
        normalized_left_exponent
    } else {
        normalized_right_exponent
    };
    let aligned_left =
        normalized_left << (normalized_left_exponent - common_exponent).cast_unsigned();
    let aligned_right =
        normalized_right << (normalized_right_exponent - common_exponent).cast_unsigned();
    if aligned_left < aligned_right {
        -1
    } else if aligned_left > aligned_right {
        1
    } else {
        0
    }
}

/// Compare the exact-real product `left × right` with its rounded result.
///
/// Nonzero finite operands contribute at most 53 significand bits each, so the
/// exact product fits in `u128`. Above the residual-underflow threshold,
/// `TwoProductFMA` gives its exact rounding error; see `REFERENCES.md` \[18\].
/// Integer comparison handles smaller products, including rounding to zero.
/// Callers supply nonzero finite operands and their finite rounded product.
#[inline]
pub(crate) const fn compare_product_with_rounded(left: f64, right: f64, rounded: f64) -> i8 {
    // If rounded has leading exponent E, the exact product has at most 106
    // bits and leading exponent at least E-1 (a binade-crossing round-up).
    // Its least significant bit is therefore at least 2^(E-106). E >= -968
    // puts that bit at or above 2^-1074, so the FMA residual cannot underflow.
    // A zero residual below this threshold would NOT establish exactness.
    const MIN_EXACT_RESIDUAL_PRODUCT: f64 = f64::from_bits((1023 - 968) << 52);
    if rounded.abs() >= MIN_EXACT_RESIDUAL_PRODUCT {
        // Scaling by a normal power of two only changes the exponent in this
        // range. This also avoids an FMA dependency for exact unit products.
        const FRACTION_MASK: u64 = 0x000f_ffff_ffff_ffff;
        if left.to_bits() & FRACTION_MASK == 0 || right.to_bits() & FRACTION_MASK == 0 {
            return 0;
        }
        let residual = left.mul_add(right, -rounded);
        return if residual < 0.0 {
            -1
        } else if residual > 0.0 {
            1
        } else {
            0
        };
    }

    compare_product_with_rounded_integer(left, right, rounded)
}

/// Keep the underflow-sensitive integer path out of each inlined normal product.
#[cold]
#[inline(never)]
const fn compare_product_with_rounded_integer(left: f64, right: f64, rounded: f64) -> i8 {
    let negative = left.is_sign_negative() != right.is_sign_negative();
    if rounded == 0.0 {
        return if negative { -1 } else { 1 };
    }

    let (left_significand, left_exponent) = decompose_magnitude(left);
    let (right_significand, right_exponent) = decompose_magnitude(right);
    let exact_significand = left_significand * right_significand;
    let exact_exponent = left_exponent + right_exponent;
    let (rounded_significand, rounded_exponent) = decompose_magnitude(rounded);
    let magnitude_relation = compare_binary_magnitudes(
        exact_significand,
        exact_exponent,
        rounded_significand,
        rounded_exponent,
    );

    if negative {
        -magnitude_relation
    } else {
        magnitude_relation
    }
}

#[cfg(test)]
mod tests {
    use super::{compare_binary_magnitudes, compare_product_with_rounded};

    #[test]
    fn product_comparison_preserves_tiny_residuals_and_const_evaluation() {
        const LEFT: f64 = 1.0_f64.next_up();
        const RELATION: i8 = compare_product_with_rounded(LEFT, LEFT, LEFT * LEFT);
        assert_eq!(RELATION, 1);
        assert_eq!(compare_product_with_rounded(-LEFT, LEFT, -LEFT * LEFT), -1);

        let tiny = f64::from_bits((1023 - 1000) << 52).next_up();
        let rounded = tiny * LEFT;
        assert!(rounded.is_normal());
        // The exact residual is positive but lies below the least subnormal.
        assert_eq!(tiny.mul_add(LEFT, -rounded).to_bits() << 1, 0);
        assert_eq!(compare_product_with_rounded(tiny, LEFT, rounded), 1);
        assert_eq!(compare_product_with_rounded(-tiny, LEFT, -rounded), -1);
    }

    #[test]
    fn binary_magnitude_comparison_orders_distinct_top_exponents() {
        assert_eq!(compare_binary_magnitudes(1, 1, 1, 0), 1);
        assert_eq!(compare_binary_magnitudes(1, 0, 1, 1), -1);
    }
}
