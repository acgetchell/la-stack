#![forbid(unsafe_code)]

//! Exact range-boundary fallback for the otherwise approximate Euclidean norm.

use core::cmp::Ordering;

use crate::rounding::decompose_magnitude;
use crate::{ArithmeticOperation, LaError};

/// Squared binary64 values are integer multiples of `2^-2148`.
const SQUARE_BASE_EXPONENT: i64 = -2148;
/// Each square is below `2^2048`, requiring 4196 bits at the chosen scale.
/// Summing at most `usize::MAX` squares needs at most `usize::BITS` extra bits.
const SUM_WORDS: usize = (4196 + usize::BITS as usize).div_ceil(u64::BITS as usize);

/// Non-negative exact sum of squares, in little-endian base-2^64 storage.
///
/// The fixed capacity covers every finite binary64 vector length representable
/// by `usize`. No term is discarded, including subnormal squares that can break
/// an otherwise exact rounding tie. Representation follows IEEE 754 \[9-10\].
struct SquareSum {
    words: [u64; SUM_WORDS],
}

impl SquareSum {
    const ZERO: Self = Self {
        words: [0; SUM_WORDS],
    };

    /// Add one word and propagate its carry within the proven sum capacity.
    fn add_word(&mut self, mut index: usize, mut word: u64) {
        while word != 0 {
            let (sum, carry) = self.words[index].overflowing_add(word);
            self.words[index] = sum;
            word = u64::from(carry);
            index += 1;
        }
    }

    /// Add `(significand × 2^exponent)²` without rounding.
    ///
    /// Inputs are binary64 magnitudes or upper-range rounding midpoints, so
    /// the significand has at most 54 bits and its square fits `u128`.
    /// Their exponents are at least -1074, so all shifts are non-negative;
    /// all shifted squares and carries fit the fixed storage above.
    #[expect(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "extracting low/high words intentionally truncates; finite binary64 square exponents give shifts in 0..=4090"
    )]
    fn add_square(&mut self, significand: u128, exponent: i64) {
        let square = significand * significand;
        let shift = (2 * exponent - SQUARE_BASE_EXPONENT) as usize;
        let index = shift / u64::BITS as usize;
        let offset = shift % u64::BITS as usize;
        let low = square as u64;
        let high = (square >> u64::BITS) as u64;

        self.add_word(index, low << offset);
        if offset == 0 {
            self.add_word(index + 1, high);
        } else {
            self.add_word(index + 1, (high << offset) | (low >> (64 - offset)));
            self.add_word(index + 2, high >> (64 - offset));
        }
    }

    /// Compare with one squared positive normal binary64 value.
    fn compare_square(&self, value: f64) -> Ordering {
        let (significand, exponent) = decompose_magnitude(value);
        self.compare_scaled_square(significand, exponent)
    }

    /// Compare with the squared midpoint above an upper-range normal value.
    ///
    /// If `lower = m × 2^e`, its next rounding boundary is
    /// `(2m + 1) × 2^(e - 1)`, including the boundary above `f64::MAX`.
    fn compare_midpoint(&self, lower: f64) -> Ordering {
        let (significand, exponent) = decompose_magnitude(lower);
        self.compare_scaled_square(2 * significand + 1, exponent - 1)
    }

    fn compare_scaled_square(&self, significand: u128, exponent: i64) -> Ordering {
        let mut other = Self::ZERO;
        other.add_square(significand, exponent);
        self.words.iter().rev().cmp(other.words.iter().rev())
    }
}

/// Round a potentially overflowing norm by comparing exact squared values.
///
/// `scale` is the largest input magnitude and is normal in this cold path.
/// The overflow midpoint is `f64::MAX + 2^970`. Its tie rounds to infinity;
/// smaller norms round to a finite value. Binary search and an exact squared
/// midpoint comparison then select the nearest finite value, ties to even.
/// This does not require a square-root approximation or optional dependencies.
#[cold]
pub(crate) fn norm_near_overflow<const D: usize>(
    values: &[f64; D],
    scale: f64,
) -> Result<f64, LaError> {
    let mut sum = SquareSum::ZERO;
    for &value in values {
        if value != 0.0 {
            let (significand, exponent) = decompose_magnitude(value);
            sum.add_square(significand, exponent);
        }
    }

    if sum.compare_midpoint(f64::MAX) != Ordering::Less {
        return Err(LaError::non_finite_computation_scalar(
            ArithmeticOperation::VectorNorm,
        ));
    }
    if sum.compare_square(f64::MAX) != Ordering::Less {
        return Ok(f64::MAX);
    }

    // The exact norm is at least its largest coordinate. Positive binary64 bit
    // patterns are monotonically ordered, and every candidate here is normal.
    let mut lower = scale.to_bits();
    let mut upper = f64::MAX.to_bits();
    while upper - lower > 1 {
        let middle = lower + (upper - lower) / 2;
        match sum.compare_square(f64::from_bits(middle)) {
            Ordering::Less => upper = middle,
            Ordering::Greater => lower = middle,
            Ordering::Equal => return Ok(f64::from_bits(middle)),
        }
    }

    let rounded = match sum.compare_midpoint(f64::from_bits(lower)) {
        Ordering::Less => lower,
        Ordering::Greater => upper,
        Ordering::Equal => lower + (lower & 1),
    };
    Ok(f64::from_bits(rounded))
}
