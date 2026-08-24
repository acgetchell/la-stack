pub fn forbidden_f64_operations(left: f64, right: f64) -> [f64; 5] {
    [
        // ruleid: la-stack.rust.no-algebraic-float-operations
        left.algebraic_add(right),
        // ruleid: la-stack.rust.no-algebraic-float-operations
        left.algebraic_sub(right),
        // ruleid: la-stack.rust.no-algebraic-float-operations
        left.algebraic_mul(right),
        // ruleid: la-stack.rust.no-algebraic-float-operations
        left.algebraic_div(right),
        // ruleid: la-stack.rust.no-algebraic-float-operations
        left.algebraic_rem(right),
    ]
}

pub fn forbidden_f64_associated_operation(left: f64, right: f64) -> f64 {
    // ruleid: la-stack.rust.no-algebraic-float-operations
    f64::algebraic_add(left, right)
}

pub fn forbidden_f64_function_item() -> fn(f64, f64) -> f64 {
    // ruleid: la-stack.rust.no-algebraic-float-operations
    f64::algebraic_sub
}

pub fn forbidden_f64_qualified_call(left: f64, right: f64) -> f64 {
    // ruleid: la-stack.rust.no-algebraic-float-operations
    <f64>::algebraic_add(left, right)
}

pub fn forbidden_f64_qualified_function_item() -> fn(f64, f64) -> f64 {
    // ruleid: la-stack.rust.no-algebraic-float-operations
    <f64>::algebraic_sub
}

pub fn forbidden_f64_reduction(values: &[f64]) -> Option<f64> {
    // ruleid: la-stack.rust.no-algebraic-float-operations
    values.iter().copied().reduce(f64::algebraic_mul)
}

pub fn permitted_f64_operations(left: f64, right: f64) -> [f64; 6] {
    [
        // ok: la-stack.rust.no-algebraic-float-operations
        left + right,
        // ok: la-stack.rust.no-algebraic-float-operations
        left - right,
        // ok: la-stack.rust.no-algebraic-float-operations
        left * right,
        // ok: la-stack.rust.no-algebraic-float-operations
        left / right,
        // ok: la-stack.rust.no-algebraic-float-operations
        left % right,
        // ok: la-stack.rust.no-algebraic-float-operations
        left.mul_add(right, 1.0),
    ]
}

pub fn permitted_f64_function_item() -> fn(f64, f64) -> f64 {
    fn add(left: f64, right: f64) -> f64 {
        left + right
    }

    // ok: la-stack.rust.no-algebraic-float-operations
    add
}

pub fn permitted_f64_qualified_fma(
    left: f64,
    right: f64,
) -> (f64, fn(f64, f64, f64) -> f64) {
    (
        // ok: la-stack.rust.no-algebraic-float-operations
        <f64>::mul_add(left, right, 1.0),
        // ok: la-stack.rust.no-algebraic-float-operations
        <f64>::mul_add,
    )
}
