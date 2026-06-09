#![forbid(unsafe_code)]

fn benchmark_and_example_fixture(result: Result<u32, &'static str>, value: Option<u32>) {
    // ruleid: la-stack.rust.no-unwrap-expect-in-benches-examples
    let _ = result.unwrap();

    // ruleid: la-stack.rust.no-unwrap-expect-in-benches-examples
    let _ = value.expect("benchmarks and examples should avoid expect");

    // ok: la-stack.rust.no-unwrap-expect-in-benches-examples
    let _ = result.unwrap_or(0);
}
