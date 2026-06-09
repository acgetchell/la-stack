#![forbid(unsafe_code)]

//! Shared helpers for exact-arithmetic benchmark input generation and tests.

use std::fmt::{self, Display};
use std::num::NonZeroU64;

/// Configuration errors for exact-arithmetic benchmark input generation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ExactBenchConfigError {
    /// The random input corpus length was zero.
    EmptyCorpus,
    /// An ordered inclusive range produced an invalid non-zero sampling width.
    InvalidRangeWidth {
        /// Inclusive lower bound.
        min: i16,
        /// Inclusive upper bound.
        max: i16,
        /// Computed inclusive width before conversion to the cached sampling width.
        width: i32,
    },
    /// The inclusive lower bound was greater than the inclusive upper bound.
    UnorderedRange {
        /// Inclusive lower bound.
        min: i16,
        /// Inclusive upper bound.
        max: i16,
    },
}

impl Display for ExactBenchConfigError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::EmptyCorpus => f.write_str("random input corpus must be nonempty"),
            Self::InvalidRangeWidth { min, max, width } => {
                write!(
                    f,
                    "random integer range {min}..={max} produced invalid sampling width {width}"
                )
            }
            Self::UnorderedRange { min, max } => {
                write!(f, "random integer range must be ordered: {min}..={max}")
            }
        }
    }
}

impl std::error::Error for ExactBenchConfigError {}

/// Inclusive integer range used by the fixed-seed exact benchmark generator.
#[derive(Clone, Copy)]
#[must_use]
pub struct I16Range {
    min: i16,
    width: NonZeroU64,
}

impl I16Range {
    /// Validate an inclusive `i16` range and cache its sampling width.
    ///
    /// # Errors
    ///
    /// Returns [`ExactBenchConfigError::UnorderedRange`] when `min > max`, or
    /// [`ExactBenchConfigError::InvalidRangeWidth`] if the inclusive range width
    /// cannot be represented as a non-zero sampling width.
    pub fn new(min: i16, max: i16) -> Result<Self, ExactBenchConfigError> {
        if min > max {
            return Err(ExactBenchConfigError::UnorderedRange { min, max });
        }

        let raw_width = i32::from(max) - i32::from(min) + 1;
        let width =
            u64::try_from(raw_width).map_err(|_| ExactBenchConfigError::InvalidRangeWidth {
                min,
                max,
                width: raw_width,
            })?;
        let width = NonZeroU64::new(width).ok_or(ExactBenchConfigError::InvalidRangeWidth {
            min,
            max,
            width: raw_width,
        })?;
        Ok(Self { min, width })
    }
}

/// Deterministic `SplitMix64` generator for reproducible benchmark corpora.
#[must_use]
pub struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    /// Initialize the generator with a fixed state.
    pub const fn new(state: u64) -> Self {
        Self { state }
    }

    /// Advance the generator and return the next 64 random bits.
    const fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    #[allow(clippy::cast_possible_truncation)]
    /// Draw a random `i16` inside a validated inclusive range.
    #[must_use]
    pub fn next_i16(&mut self, range: I16Range) -> i16 {
        let offset = (self.next_u64() % range.width.get()) as i32;
        let value = i32::from(range.min) + offset;
        value as i16
    }
}
