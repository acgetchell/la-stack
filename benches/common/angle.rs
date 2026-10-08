#![forbid(unsafe_code)]

//! Independently checked fixtures shared by angle benchmarks and their tests.

use core::f64::consts::{FRAC_PI_4, PI};

use la_stack::{Vector, angle_between};

use crate::bench_utils::OrAbort;

/// The ordinary max-scaled Kahan control is valid for these vetted fixtures.
/// It is deliberately not a general-purpose API: pre-normalization can lose
/// subnormal differences. No input checks are included in this kernel control.
pub(crate) fn stable_control<const D: usize>(left: &[f64; D], right: &[f64; D]) -> f64 {
    let maximum = |values: &[f64; D]| values.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
    let left_scale = maximum(left);
    let right_scale = maximum(right);
    let mut left_norm = 0.0_f64;
    let mut right_norm = 0.0_f64;
    for (&x, &y) in left.iter().zip(right) {
        left_norm = left_norm.hypot(x / left_scale);
        right_norm = right_norm.hypot(y / right_scale);
    }
    let mut difference = 0.0_f64;
    let mut sum = 0.0_f64;
    for (&x, &y) in left.iter().zip(right) {
        let x = (x / left_scale) * right_norm;
        let y = (y / right_scale) * left_norm;
        difference = difference.hypot(x - y);
        sum = sum.hypot(x + y);
    }
    (2.0 * difference * sum).atan2((sum - difference) * (sum + difference))
}

/// Private construction makes numerical validation a prerequisite to timing.
pub(crate) struct AngleInput<const D: usize> {
    name: &'static str,
    left: Vector<D>,
    right: Vector<D>,
}

impl<const D: usize> AngleInput<D> {
    fn new(name: &'static str, left: [f64; D], right: [f64; D], expected: f64) -> Self {
        let left = Vector::try_new(left).or_abort("finite benchmark fixture");
        let right = Vector::try_new(right).or_abort("finite benchmark fixture");
        let tolerance = expected.abs() * (8.0 * f64::EPSILON);
        for actual in [
            angle_between(left.as_array(), right.as_array()).or_abort("valid angle"),
            left.angle(&right).or_abort("valid vector angle"),
            stable_control(left.as_array(), right.as_array()),
        ] {
            assert!(
                (actual - expected).abs() <= tolerance,
                "{name}: {actual:e} vs {expected:e}"
            );
        }
        Self { name, left, right }
    }

    pub(crate) const fn name(&self) -> &'static str {
        self.name
    }
    pub(crate) const fn left(&self) -> &Vector<D> {
        &self.left
    }
    pub(crate) const fn right(&self) -> &Vector<D> {
        &self.right
    }
}

#[expect(
    clippy::cast_precision_loss,
    reason = "benchmark dimensions are 3 through 6"
)]
pub(crate) fn fixtures<const D: usize>() -> [AngleInput<D>; 5] {
    assert!((3..=6).contains(&D));
    let planar = |x, y| {
        let mut values = [0.0; D];
        values[0] = x;
        values[D - 1] = y;
        values
    };
    let mut dense_right = [1.0; D];
    dense_right[0] = -1.0;
    // Integer dot is D-2; D-1 nonzero exterior minors each have magnitude 2.
    let dense_reference = (4.0 * (D - 1) as f64).sqrt().atan2((D - 2) as f64);
    let t = 2.0_f64.powi(-30);
    [
        AngleInput::new("dense", [1.0; D], dense_right, dense_reference),
        AngleInput::new("near_parallel", planar(1.0, 0.0), planar(1.0, t), t.atan()),
        AngleInput::new(
            "near_antipodal",
            planar(1.0, 0.0),
            planar(-1.0, t),
            PI - t.atan(),
        ),
        AngleInput::new(
            "mixed_scale",
            planar(f64::MAX, 0.0),
            planar(f64::MIN_POSITIVE, f64::MIN_POSITIVE),
            FRAC_PI_4,
        ),
        AngleInput::new(
            "subnormal",
            planar(1.0, 0.0),
            planar(1.0, f64::from_bits(1)),
            f64::from_bits(1),
        ),
    ]
}
