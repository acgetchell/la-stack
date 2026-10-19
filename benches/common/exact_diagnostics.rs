#![forbid(unsafe_code)]

//! Independently checked conversion and determinant diagnostics.

use core::array::from_fn;

use la_stack::{
    BigInt, BigRational, DeterminantSign, ExactF64Conversion, LaError, Matrix, RationalVector,
    UnrepresentableReason,
};
use num_bigint::Sign;

use crate::bench_utils::OrAbort;
use crate::rational_bench::rational_determinant_gaussian;

/// Canonical conversion workloads, including expected strict rejections.
#[derive(Clone, Copy, Debug)]
pub enum ConversionKind {
    /// Every component is exactly one half.
    Dyadic,
    /// A final one-third component needs rounding.
    NonDyadic,
    /// A final ratio has wide numerator and denominator storage.
    Wide256,
    /// A final ratio has still wider numerator and denominator storage.
    Wide1024,
    /// The final component is the smallest positive binary64 subnormal.
    MinSubnormal,
    /// A negative half-subnormal rounds to negative zero.
    NegativeUnderflow,
    /// An integer just below the overflow midpoint rounds to `f64::MAX`.
    BelowOverflow,
    /// The exact overflow midpoint has no finite rounded result.
    OverflowMidpoint,
}

impl ConversionKind {
    /// Every diagnostic family, in stable order.
    pub const ALL: [Self; 8] = [
        Self::Dyadic,
        Self::NonDyadic,
        Self::Wide256,
        Self::Wide1024,
        Self::MinSubnormal,
        Self::NegativeUnderflow,
        Self::BelowOverflow,
        Self::OverflowMidpoint,
    ];

    /// Stable group label.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Dyadic => "dyadic",
            Self::NonDyadic => "non_dyadic",
            Self::Wide256 => "wide256",
            Self::Wide1024 => "wide1024",
            Self::MinSubnormal => "min_subnormal",
            Self::NegativeUnderflow => "negative_underflow",
            Self::BelowOverflow => "below_overflow",
            Self::OverflowMidpoint => "overflow_midpoint",
        }
    }
}

type ConversionBits = Result<u64, UnrepresentableReason>;

/// A final component and independently known strict/rounded binary64 outcomes.
fn conversion_component(kind: ConversionKind) -> (BigRational, ConversionBits, ConversionBits) {
    let half = 0.5_f64.to_bits();
    let requires_rounding = Err(UnrepresentableReason::RequiresRounding);
    let not_finite = Err(UnrepresentableReason::NotFinite);
    match kind {
        ConversionKind::Dyadic => (BigRational::new(1.into(), 2.into()), Ok(half), Ok(half)),
        ConversionKind::NonDyadic => (
            BigRational::new(1.into(), 3.into()),
            requires_rounding,
            Ok(0x3fd5_5555_5555_5555),
        ),
        ConversionKind::Wide256 | ConversionKind::Wide1024 => {
            let bits = if matches!(kind, ConversionKind::Wide256) {
                256_u32
            } else {
                1024
            };
            (
                BigRational::new(
                    (BigInt::from(1_u8) << bits) + BigInt::from(1_u8),
                    (BigInt::from(1_u8) << (bits + 1)) - BigInt::from(1_u8),
                ),
                requires_rounding,
                Ok(half),
            )
        }
        ConversionKind::MinSubnormal => (
            BigRational::new(1.into(), BigInt::from(1_u8) << 1074_u32),
            Ok(1),
            Ok(1),
        ),
        ConversionKind::NegativeUnderflow => (
            BigRational::new((-1).into(), BigInt::from(1_u8) << 1075_u32),
            requires_rounding,
            Ok(1_u64 << 63),
        ),
        ConversionKind::BelowOverflow | ConversionKind::OverflowMidpoint => {
            // f64::MAX + half its final ULP = 2^1024 - 2^970.
            let midpoint = (BigInt::from(1_u8) << 1024_u32) - (BigInt::from(1_u8) << 970_u32);
            if matches!(kind, ConversionKind::BelowOverflow) {
                (
                    BigRational::from_integer(midpoint - 1_u8),
                    requires_rounding,
                    Ok(f64::MAX.to_bits()),
                )
            } else {
                (BigRational::from_integer(midpoint), not_finite, not_finite)
            }
        }
    }
}

/// Construct canonical values and verify strict/rounded bits and typed errors.
///
/// # Panics
/// Panics for an empty vector or a mismatch with the expected typed result.
pub fn canonical_conversion_input<const D: usize>(kind: ConversionKind) -> RationalVector<D> {
    assert!(D > 0);
    let mut data = from_fn(|_| BigRational::new(1.into(), 2.into()));
    let (component, strict, rounded) = conversion_component(kind);
    data[D - 1] = component;
    let expected = |outcome: ConversionBits| {
        outcome
            .map(|bits| {
                let mut values = [0.5_f64.to_bits(); D];
                values[D - 1] = bits;
                values
            })
            .map_err(|reason| LaError::unrepresentable(Some(D - 1), reason))
    };
    let vector = RationalVector::try_new(data).or_abort("canonical conversion input");
    assert_eq!(
        vector
            .try_to_f64()
            .map(|v| v.into_array().map(f64::to_bits)),
        expected(strict)
    );
    assert_eq!(
        vector
            .as_array()
            .try_to_f64()
            .map(|v| v.into_array().map(f64::to_bits)),
        expected(strict)
    );
    assert_eq!(
        vector
            .to_rounded_f64()
            .map(|v| v.into_array().map(f64::to_bits)),
        expected(rounded)
    );
    assert_eq!(
        vector
            .as_array()
            .to_rounded_f64()
            .map(|v| v.into_array().map(f64::to_bits)),
        expected(rounded)
    );
    vector
}

/// Small D=4 shapes that must remain distinct from dense determinant workloads.
#[derive(Clone, Copy, Debug)]
pub enum Det4Kind {
    /// All first-row coefficients are active.
    Dense,
    /// Only one first-row cofactor is needed.
    Sparse,
    /// Repeated rows force the exact sign fallback.
    Singular,
    /// A positive nonzero determinant that requires exact sign fallback.
    NearSingularPositive,
    /// A row swap reverses the near-singular determinant's exact sign.
    NearSingularNegative,
    /// Exactly scaled rows span subnormals through 2^900, preserving determinant.
    MixedExponents,
    /// Dense, extreme diagonal entries produce an exact result beyond binary64.
    LargeEntries,
}

impl Det4Kind {
    /// Every shape and adversarial family, in stable order.
    pub const ALL: [Self; 7] = [
        Self::Dense,
        Self::Sparse,
        Self::Singular,
        Self::NearSingularPositive,
        Self::NearSingularNegative,
        Self::MixedExponents,
        Self::LargeEntries,
    ];

    /// Stable group label.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Dense => "dense",
            Self::Sparse => "sparse",
            Self::Singular => "singular",
            Self::NearSingularPositive => "near_singular_positive",
            Self::NearSingularNegative => "near_singular_negative",
            Self::MixedExponents => "mixed_exponents",
            Self::LargeEntries => "large_entries",
        }
    }
}

const DENSE_ROWS: [[f64; 4]; 4] = [
    [11.0, 2.0, -3.0, 4.0],
    [2.0, 13.0, 5.0, -1.0],
    [3.0, -2.0, 17.0, 6.0],
    [-1.0, 4.0, 2.0, 19.0],
];

/// Binary exponent spreads from equal-scale rows to nearly the full f64 range.
pub const EXPONENT_SPREADS: [i32; 9] = [0, 8, 32, 64, 128, 256, 512, 1024, 2000];

/// Which axis carries the power-of-two scaling in determinant diagnostics.
#[derive(Clone, Copy, Debug)]
pub enum ExponentAxis {
    /// Each row has its own scale, removable by row-local integer lifting.
    Row,
    /// Each column has its own scale, retaining wide exponents within rows.
    Column,
}

impl ExponentAxis {
    /// Stable benchmark label.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Row => "row",
            Self::Column => "column",
        }
    }
}

/// Validate an exactly scaled dense matrix against rational Gaussian elimination.
///
/// Symmetric exponents sum to zero, so scaling preserves the determinant.
/// Column scaling is a control for exponent variation that row scaling cannot
/// remove. The base matrix is strictly diagonally dominant with integer entries.
///
/// # Panics
/// Panics outside D=2–5 and the declared spreads, or if any exact determinant
/// API disagrees with the independent oracle. Validation precedes all timing.
pub fn exponent_spread_input<const D: usize>(spread: i32, axis: ExponentAxis) -> Matrix<D> {
    assert!((2..=5).contains(&D));
    assert!(EXPONENT_SPREADS.contains(&spread));
    let dimension = i32::try_from(D).or_abort("small diagnostic dimension");
    let base: [[f64; D]; D] = from_fn(|row| {
        from_fn(|column| {
            if row == column {
                f64::from(2 * dimension + 3)
            } else {
                f64::from(i32::try_from((row + 2 * column) % 5).or_abort("small residue") - 2)
            }
        })
    });
    let scales: [f64; D] = from_fn(|index| {
        let index = i32::try_from(index).or_abort("small diagnostic index");
        let exponent = (spread / 2) * (2 * index - dimension + 1) / (dimension - 1);
        let biased = u64::try_from(exponent + 1023).or_abort("normal diagnostic scale");
        f64::from_bits(biased << 52)
    });
    let rows = from_fn(|row| {
        from_fn(|column| {
            let scale = match axis {
                ExponentAxis::Row => scales[row],
                ExponentAxis::Column => scales[column],
            };
            base[row][column] * scale
        })
    });
    let rational_rows = |entries: [[f64; D]; D]| {
        entries.map(|row| row.map(|value| BigRational::from_float(value).or_abort("finite entry")))
    };
    let expected = rational_determinant_gaussian(rational_rows(base));
    assert_eq!(rational_determinant_gaussian(rational_rows(rows)), expected);
    assert_eq!(expected.numer().sign(), Sign::Plus);
    let matrix = Matrix::try_from_rows(rows).or_abort("finite scaled matrix");
    let actual = matrix.det_exact().or_abort("scaled determinant");
    assert_eq!(actual.numer(), expected.numer());
    assert_eq!(actual.denom(), expected.denom());
    assert_eq!(matrix.det_sign_exact(), DeterminantSign::Positive);
    assert_eq!(matrix.det_exact_f64(), expected.try_to_f64());
    assert_eq!(matrix.det_exact_rounded_f64(), expected.to_rounded_f64());
    matrix
}

/// Build exactly representable determinant fixtures, including row scalings.
fn det4_rows(kind: Det4Kind) -> [[f64; 4]; 4] {
    let mut rows = DENSE_ROWS;
    match kind {
        Det4Kind::Dense => {}
        Det4Kind::Sparse => rows[0] = [0.0, 2.0, 0.0, 0.0],
        Det4Kind::Singular => rows[1] = rows[0],
        Det4Kind::NearSingularPositive | Det4Kind::NearSingularNegative => {
            let perturbation = f64::from_bits(0x3cd0_0000_0000_0000); // 2^-50
            rows = [
                [1.0, 1.0, 1.0, 1.0],
                [1.0, 1.0 + perturbation, 1.0, 1.0],
                [1.0, 1.0, 2.0, 1.0],
                [1.0, 1.0, 1.0, 2.0],
            ];
            if matches!(kind, Det4Kind::NearSingularNegative) {
                rows.swap(2, 3);
            }
        }
        Det4Kind::MixedExponents => {
            // Exact row scales 2^[900, -1074, 700, -526] multiply to one.
            let scales = [
                f64::from_bits((1023 + 900) << 52),
                f64::from_bits(1),
                f64::from_bits((1023 + 700) << 52),
                f64::from_bits((1023 - 526) << 52),
            ];
            for (row, scale) in rows.iter_mut().zip(scales) {
                for value in row {
                    *value *= scale;
                }
            }
        }
        Det4Kind::LargeEntries => {
            let big = f64::MAX / 2.0;
            rows = from_fn(|i| from_fn(|j| if i == j { big } else { 1.0 }));
        }
    }
    rows
}

/// Validate determinant-only fixtures against rational Gaussian elimination.
///
/// # Panics
/// Panics if either exact API disagrees with the independent determinant.
pub fn exact_det4_input(kind: Det4Kind) -> Matrix<4> {
    let rows = det4_rows(kind);
    let exact_rows =
        rows.map(|row| row.map(|value| BigRational::from_float(value).or_abort("exact input")));
    let expected = rational_determinant_gaussian(exact_rows);
    if matches!(kind, Det4Kind::MixedExponents) {
        let base = DENSE_ROWS
            .map(|row| row.map(|value| BigRational::from_float(value).or_abort("dense input")));
        assert_eq!(expected, rational_determinant_gaussian(base));
    }
    if matches!(kind, Det4Kind::LargeEntries) {
        assert_eq!(
            expected.try_to_f64(),
            Err(LaError::unrepresentable(
                None,
                UnrepresentableReason::NotFinite
            ))
        );
    }
    if matches!(
        kind,
        Det4Kind::NearSingularPositive | Det4Kind::NearSingularNegative
    ) {
        let numerator = if matches!(kind, Det4Kind::NearSingularPositive) {
            1
        } else {
            -1
        };
        assert_eq!(
            expected,
            BigRational::new(numerator.into(), BigInt::from(1_u8) << 50_u32)
        );
    }
    let expected_sign = match expected.numer().sign() {
        Sign::Minus => DeterminantSign::Negative,
        Sign::NoSign => DeterminantSign::Zero,
        Sign::Plus => DeterminantSign::Positive,
    };
    let matrix = Matrix::try_from_rows(rows).or_abort("determinant diagnostic input");
    assert_eq!(
        matrix.det_exact().or_abort("determinant diagnostic"),
        expected
    );
    assert_eq!(matrix.det_sign_exact(), expected_sign);
    if matches!(
        kind,
        Det4Kind::Singular | Det4Kind::NearSingularPositive | Det4Kind::NearSingularNegative
    ) {
        let estimate = matrix
            .det_direct_with_errbound()
            .or_abort("finite diagnostic filter");
        assert!(
            estimate.is_none_or(
                |estimate| estimate.determinant().abs() <= estimate.absolute_error_bound()
            )
        );
    }
    matrix
}
