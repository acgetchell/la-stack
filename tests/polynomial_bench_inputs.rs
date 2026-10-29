#![forbid(unsafe_code)]

//! Exercise the full lifted consumer shape through geometric dimension six.

#![cfg(feature = "exact")]

#[path = "../benches/common/bench_utils.rs"]
mod bench_utils;
#[path = "../benches/common/polynomial.rs"]
mod polynomial;

use la_stack::polynomial_determinant;
use pastey::paste;
use polynomial::{BUDGET, coincident_lifted};

macro_rules! lifted_test {
    ($n:literal) => {
        paste! {
            #[test]
            fn [<coincident_lifted_order_ $n>]() {
                let mut input = coincident_lifted::<$n>();
                // Repeating a complete polynomial row must cancel every term.
                let (first, rest) = input.matrix.split_at_mut(1);
                rest[0].clone_from(&first[0]);
                assert!(polynomial_determinant(&input.matrix, BUDGET).unwrap().is_zero());
            }
        }
    };
}
lifted_test!(4);
lifted_test!(5);
lifted_test!(6);
lifted_test!(7);
lifted_test!(8);

#[test]
fn consumer_product_envelope_matches_documented_policy() {
    let maxima = [2_usize, 2, 2, 2, 2, 2, 13, 1];
    let mut total = 0;
    for mask in 0_usize..255 {
        let size = mask.count_ones() as usize;
        let partial = (1..=size).product::<usize>()
            * maxima
                .iter()
                .enumerate()
                .filter(|(i, _)| mask & (1 << i) != 0)
                .map(|(_, &count)| count)
                .product::<usize>();
        let choices: usize = maxima
            .iter()
            .enumerate()
            .filter(|(i, _)| mask & (1 << i) == 0)
            .map(|(_, &count)| count)
            .sum();
        total += partial * choices;
    }
    assert_eq!(total, BUDGET.max_term_products);
}
