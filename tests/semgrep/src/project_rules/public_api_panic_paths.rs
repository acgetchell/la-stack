#![forbid(unsafe_code)]

pub enum Error {
    Invalid,
}

// ruleid: la-stack.rust.no-public-api-panic-paths
pub fn panics_on_input(value: usize) -> usize {
    if value == 0 {
        panic!("zero is invalid");
    }
    value
}

// ruleid: la-stack.rust.no-public-api-panic-paths
pub const fn asserts_on_input(value: usize) -> usize {
    assert!(value > 0);
    value
}

// ruleid: la-stack.rust.no-public-api-panic-paths
pub fn unwraps_on_input(value: Option<usize>) -> usize {
    value.unwrap()
}

// ruleid: la-stack.rust.no-public-api-panic-paths
pub async fn expects_on_input(value: Option<usize>) -> usize {
    value.expect("value is required")
}

// ok: la-stack.rust.no-public-api-panic-paths
pub fn fallible_result(value: usize) -> Result<usize, Error> {
    if value == 0 {
        Err(Error::Invalid)
    } else {
        Ok(value)
    }
}

// ok: la-stack.rust.no-public-api-panic-paths
pub fn fallible_option(value: usize) -> Option<usize> {
    if value == 0 { None } else { Some(value) }
}

// ok: la-stack.rust.no-public-api-panic-paths
pub fn total(value: usize) -> usize {
    value + 1
}

// ok: la-stack.rust.no-public-api-panic-paths
pub(crate) fn crate_private_literal_helper(value: usize) -> usize {
    assert!(value > 0);
    value
}

#[cfg(test)]
mod tests {
    // ok: la-stack.rust.no-public-api-panic-paths
    pub(crate) fn test_only_helper(value: usize) -> usize {
        value.expect("test fixture has a value")
    }
}
