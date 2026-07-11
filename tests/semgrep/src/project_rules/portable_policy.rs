#![allow(dead_code, unused_imports)]

use num_traits::NumCast;

// ruleid: la-stack.rust.no-module-scope-cfg-test-use
#[cfg(test)]
use crate::FixtureOnlyImport;

fn safe_f64(_value: u64) -> Option<f64> {
    Some(1.0)
}

pub fn silent_conversion_fallback(value: u64) -> f64 {
    // ruleid: la-stack.rust.no-silent-conversion-fallbacks, la-stack.rust.no-silent-conversion-fallbacks-in-public-samples
    NumCast::from(value).unwrap_or(0.0)
}

pub fn partial_cmp_ordering_default(left: f64, right: f64) -> std::cmp::Ordering {
    // ruleid: la-stack.rust.no-partial-cmp-ordering-defaults
    left.partial_cmp(&right).unwrap_or(std::cmp::Ordering::Equal)
}

pub fn function_local_use_fixture() {
    // ruleid: la-stack.rust.no-function-local-use-in-src
    use std::cmp::Ordering;

    let _ordering = Ordering::Equal;
}

// ruleid: la-stack.rust.no-public-api-cfg-test-shim
#[cfg(any(test, feature = "diagnostics"))]
pub fn public_api_test_cfg_shim_fixture() {}

// ruleid: la-stack.rust.borrowed-view-types-require-lifetime
pub struct OwnedView {
    value: usize,
}

// ok: la-stack.rust.borrowed-view-types-require-lifetime
pub struct BorrowedView<'a> {
    value: &'a usize,
}

pub struct UncheckedFixture;

impl UncheckedFixture {
    // ruleid: la-stack.rust.no-public-unchecked-apis
    pub fn from_unchecked_value() -> Self {
        Self
    }
}

pub fn ignored_question_mark_result(value: Result<u8, ()>) -> Result<(), ()> {
    // ruleid: la-stack.rust.no-ignored-fallible-results
    let _ = value?;
    Ok(())
}

// ruleid: la-stack.rust.no-box-dyn-error-in-src, la-stack.rust.no-box-dyn-error-in-examples-benches
type ErasedError = Box<dyn std::error::Error>;

// ruleid: la-stack.rust.no-clippy-allow-lints
#[allow(clippy::too_many_lines)]
fn clippy_allow_fixture() {}

// ok: la-stack.rust.no-clippy-allow-lints
#[expect(clippy::too_many_lines, reason = "fixture documents the suppression")]
fn clippy_expect_fixture() {}

// ruleid: la-stack.rust.no-ignored-tests
#[ignore = "use an explicit slow-test feature"]
fn ignored_test_fixture() {}

// ruleid: la-stack.rust.expect-requires-reason
#[expect(clippy::too_many_lines)]
fn undocumented_expectation_fixture() {}

// ok: la-stack.rust.expect-requires-reason
#[expect(clippy::too_many_lines, reason = "fixture documents the suppression")]
fn documented_expectation_fixture() {}

// ruleid: la-stack.rust.no-box-dyn-error-in-doctests
/// # fn main() -> Result<(), Box<dyn std::error::Error>> { Ok(()) }
fn doctest_erased_error_fixture() {}

// ruleid: la-stack.rust.prefer-assert-matches-in-doctests
/// assert!(matches!(value, Some(_)));
fn doctest_assert_matches_fixture() {}
