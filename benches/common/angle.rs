#![forbid(unsafe_code)]

//! Independently checked fixtures shared by angle benchmarks and their tests.

use core::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

use la_stack::{Vector, VectorAngle};

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
            left.as_array()
                .as_slice()
                .angle(right.as_array())
                .or_abort("valid angle"),
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
pub(crate) fn fixtures<const D: usize>() -> [AngleInput<D>; 14] {
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
    let mut adjacent_left = [0.0; D];
    let mut adjacent_right = [-0.0; D];
    adjacent_left[1] = 3.0;
    adjacent_left[2] = 4.0;
    adjacent_right[1] = -4.0;
    adjacent_right[2] = 3.0;
    let mut alternating_left = [0.0; D];
    let mut alternating_right = [0.0; D];
    for i in (0..D).step_by(2) {
        alternating_left[i] = 1.0;
        alternating_right[i] = if i == 0 { -1.0 } else { 1.0 };
    }
    let support = D.div_ceil(2) as f64;
    let alternating_reference = (4.0 * (support - 1.0)).sqrt().atan2(support - 2.0);
    // Deterministic varied dense coordinates, independently evaluated with
    // integer Gram products (no production minor or compensated reduction).
    let dense_left = core::array::from_fn(|i| i32::try_from(i + 1).or_abort("small dimension"));
    let dense_varied = core::array::from_fn(|i| {
        i32::try_from(D - i).or_abort("small dimension") * if i % 2 == 0 { 1 } else { -1 }
    });
    let dot: i32 = dense_left
        .iter()
        .zip(dense_varied)
        .map(|(x, y)| x * y)
        .sum();
    let left_squared: i32 = dense_left.iter().map(|x| x * x).sum();
    let right_squared: i32 = dense_varied.iter().map(|x| x * x).sum();
    let varied_reference = f64::from(left_squared * right_squared - dot.pow(2))
        .sqrt()
        .atan2(f64::from(dot));
    [
        AngleInput::new("dense", [1.0; D], dense_right, dense_reference),
        AngleInput::new(
            "dense_varied",
            dense_left.map(f64::from),
            dense_varied.map(f64::from),
            varied_reference,
        ),
        AngleInput::new("dense_identical", [1.0; D], [1.0; D], 0.0),
        AngleInput::new("dense_antipodal", [1.0; D], [-1.0; D], PI),
        AngleInput::new("sparse_adjacent", adjacent_left, adjacent_right, FRAC_PI_2),
        AngleInput::new(
            "sparse_alternating",
            alternating_left,
            alternating_right,
            alternating_reference,
        ),
        AngleInput::new("sparse_identical", planar(1.0, 0.0), planar(1.0, -0.0), 0.0),
        AngleInput::new("sparse_antipodal", planar(0.0, 1.0), planar(-0.0, -1.0), PI),
        AngleInput::new(
            "subnormal_axes",
            planar(f64::from_bits(1), 0.0),
            planar(0.0, f64::from_bits(1)),
            FRAC_PI_2,
        ),
        AngleInput::new(
            "maximum_small_angle",
            planar(f64::MAX, 0.0),
            planar(f64::MAX, f64::MAX * t),
            t.atan(),
        ),
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
