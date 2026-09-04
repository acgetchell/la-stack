#![forbid(unsafe_code)]

//! Shared binary64 rounding primitives for certified arithmetic.

/// Return the exact error in a rounded binary64 sum.
///
/// This is Knuth's `TwoSum` transform. With IEEE-754 round-to-nearest and
/// gradual underflow, `rounded + error` equals the exact-real sum whenever the
/// rounded sum is finite.
#[inline]
pub(crate) const fn two_sum_error(left: f64, right: f64, rounded: f64) -> f64 {
    let virtual_right = rounded - left;
    let virtual_left = rounded - virtual_right;
    let right_error = right - virtual_right;
    let left_error = left - virtual_left;
    left_error + right_error
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

    let left_top =
        normalized_left_exponent + (u128::BITS - normalized_left.leading_zeros() - 1) as i64;
    let right_top =
        normalized_right_exponent + (u128::BITS - normalized_right.leading_zeros() - 1) as i64;
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
#[inline]
pub(crate) const fn compare_product_with_rounded(left: f64, right: f64, rounded: f64) -> i8 {
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
    use super::compare_binary_magnitudes;

    #[test]
    fn binary_magnitude_comparison_orders_distinct_top_exponents() {
        assert_eq!(compare_binary_magnitudes(1, 1, 1, 0), 1);
        assert_eq!(compare_binary_magnitudes(1, 0, 1, 1), -1);
    }
}
