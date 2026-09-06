//! Allocation evidence for preserving canonical rational conversion proofs.

#![cfg(all(feature = "bench", feature = "exact"))]
#![forbid(unsafe_code)]

#[path = "../benches/common/bench_utils.rs"]
mod bench_utils;
#[path = "../benches/common/exact_diagnostics.rs"]
pub mod exact_diagnostics;
#[path = "../benches/common/rational.rs"]
pub mod rational_bench;

use std::hint::black_box;

use allocation_counter::measure;
use la_stack::ExactF64Conversion;
use pastey::paste;

use exact_diagnostics::{ConversionKind, canonical_conversion_input};

fn record<const D: usize>(kind: ConversionKind, operation: &str, run: impl Fn()) {
    let counts = measure(&run);
    assert_eq!(counts.count_current, 0);
    assert_eq!(counts.bytes_current, 0);
    for _ in 0..2 {
        assert_eq!(measure(&run), counts);
    }
    println!(
        "conversion_allocation,{},{D},{operation},{},{}",
        kind.name(),
        counts.count_total,
        counts.bytes_total
    );
}

fn conversion_allocations<const D: usize>() {
    for kind in ConversionKind::ALL {
        let input = canonical_conversion_input::<D>(kind);
        record::<D>(kind, "canonical", || {
            let _ = black_box(black_box(&input).try_to_f64());
        });
        record::<D>(kind, "raw", || {
            let _ = black_box(black_box(input.as_array()).try_to_f64());
        });
    }
}

macro_rules! gen_conversion_allocation_tests {
    ($d:literal) => {
        paste! {
            #[test]
            fn [<canonical_conversion_allocations_ $d d>]() {
                conversion_allocations::<$d>();
            }
        }
    };
}

gen_conversion_allocation_tests!(2);
gen_conversion_allocation_tests!(3);
gen_conversion_allocation_tests!(4);
gen_conversion_allocation_tests!(5);
