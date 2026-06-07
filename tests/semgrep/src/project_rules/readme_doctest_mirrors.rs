#![forbid(unsafe_code)]
#![allow(dead_code)]

// ruleid: la-stack.rust.no-unwrap-expect-in-readme-doctest-mirrors
mod readme_doctests_unwrap {
    #[test]
    fn readme_mirror_uses_unwrap() {
        let _ = Some(1_u32).unwrap();
    }
}

// ruleid: la-stack.rust.no-unwrap-expect-in-readme-doctest-mirrors
mod readme_doctests_expect {
    #[test]
    fn readme_mirror_uses_expect() {
        let _ = Ok::<u32, &'static str>(1).expect("README mirrors should use ?");
    }
}

// ok: la-stack.rust.no-unwrap-expect-in-readme-doctest-mirrors
mod tests {
    #[test]
    fn ordinary_internal_tests_may_use_unwrap() {
        let _ = Some(1_u32).unwrap();
    }
}

// ok: la-stack.rust.no-unwrap-expect-in-readme-doctest-mirrors
mod readme_doctests_result {
    #[test]
    fn readme_mirror_uses_result() -> Result<(), &'static str> {
        let _ = Ok::<u32, &'static str>(1)?;
        Ok(())
    }
}
