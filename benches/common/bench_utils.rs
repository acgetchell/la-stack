#![forbid(unsafe_code)]

//! Shared helpers for benchmark operations that cannot recover from failure.

use std::fmt::Display;

/// Convert a fallible benchmark operation into its successful value or abort.
pub(crate) trait OrAbort {
    /// Successful value produced by the operation.
    type Output;

    /// Return the successful value or panic with the named operation.
    ///
    /// # Panics
    ///
    /// Panics when the benchmark operation contains an error or no value.
    fn or_abort(self, operation: &str) -> Self::Output;
}

impl<T, E: Display> OrAbort for Result<T, E> {
    type Output = T;

    fn or_abort(self, operation: &str) -> Self::Output {
        match self {
            Ok(value) => value,
            Err(err) => panic!("{operation} failed: {err}"),
        }
    }
}

impl<T> OrAbort for Option<T> {
    type Output = T;

    fn or_abort(self, operation: &str) -> Self::Output {
        self.unwrap_or_else(|| panic!("{operation} returned no result"))
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn returns_result_and_option_values() {
        assert_eq!(
            <Result<_, &str> as super::OrAbort>::or_abort(Ok(7), "successful result"),
            7
        );
        assert_eq!(
            <Option<_> as super::OrAbort>::or_abort(Some(11), "present option"),
            11
        );
    }

    #[test]
    #[should_panic(expected = "fallible setup failed: fixture error")]
    fn result_error_preserves_context_and_error() {
        <Result<(), &str> as super::OrAbort>::or_abort(Err("fixture error"), "fallible setup");
    }

    #[test]
    #[should_panic(expected = "optional setup returned no result")]
    fn missing_option_preserves_context() {
        <Option<()> as super::OrAbort>::or_abort(None, "optional setup");
    }
}
