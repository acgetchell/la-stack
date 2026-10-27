#![forbid(unsafe_code)]

//! Exact rational Gauss-Jordan oracle, independent of interval endpoints,
//! rounded quotient comparisons, and the production back-substitution path.

use la_stack::prelude::*;

pub(crate) fn exact(value: f64) -> BigRational {
    BigRational::from_f64(value).expect("finite oracle input")
}

pub(crate) fn assert_encloses<const D: usize>(bounds: &[Interval; D], values: &[BigRational; D]) {
    for (interval, value) in bounds.iter().zip(values) {
        assert!(
            exact(interval.lower()) <= *value && *value <= exact(interval.upper()),
            "{interval:?} excludes {value}"
        );
    }
}

pub(crate) fn exact_solve<const D: usize>(
    mut rows: [[BigRational; D]; D],
    mut rhs: [BigRational; D],
) -> Option<[BigRational; D]> {
    let original_rows = rows.clone();
    let original_rhs = rhs.clone();
    let zero = exact(0.0);
    for column in 0..D {
        let pivot_row = (column..D).find(|&i| rows[i][column] != zero)?;
        rows.swap(column, pivot_row);
        rhs.swap(column, pivot_row);
        let pivot = rows[column][column].clone();
        for x in &mut rows[column] {
            *x /= &pivot;
        }
        rhs[column] /= &pivot;
        let pivot_entries = rows[column].clone();
        for row in 0..D {
            if row == column {
                continue;
            }
            let factor = rows[row][column].clone();
            for (entry, pivot_entry) in rows[row].iter_mut().zip(&pivot_entries) {
                *entry -= &factor * pivot_entry;
            }
            let delta = &factor * &rhs[column];
            rhs[row] -= delta;
        }
    }
    // Independently substitute the oracle result in the original equations.
    for (row, b) in original_rows.iter().zip(original_rhs) {
        assert_eq!(
            row.iter()
                .zip(&rhs)
                .map(|(a, x)| a * x)
                .sum::<BigRational>(),
            b
        );
    }
    Some(rhs)
}
