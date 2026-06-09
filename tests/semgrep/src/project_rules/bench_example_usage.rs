#![forbid(unsafe_code)]

fn benchmark_and_example_fixture(result: Result<u32, &'static str>, value: Option<u32>) {
    // ruleid: la-stack.rust.no-unwrap-expect-in-benches-examples
    let _ = result.unwrap();

    // ruleid: la-stack.rust.no-unwrap-expect-in-benches-examples
    let _ = value.expect("benchmarks and examples should avoid expect");

    // ok: la-stack.rust.no-unwrap-expect-in-benches-examples
    let _ = result.unwrap_or(0);
}

fn criterion_closure_fixture(result: Result<u32, &'static str>) {
    let bench_body = || {
        // ruleid: la-stack.rust.no-unwrap-expect-in-benches-examples
        result.expect("criterion setup should surface typed errors")
    };
    let _ = bench_body;
}

fn result_returning_example_fixture(result: Result<u32, &'static str>) -> Result<(), &'static str> {
    // ruleid: la-stack.rust.no-unwrap-expect-in-benches-examples
    let _ = result.unwrap();
    Ok(())
}

fn explicit_error_handling_fixture(result: Result<u32, &'static str>) -> Result<u32, &'static str> {
    // ok: la-stack.rust.no-unwrap-expect-in-benches-examples
    let value = result?;
    Ok(value)
}
