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

pub fn forbidden_f32_associated_operation(left: f32, right: f32) -> f32 {
    // ruleid: la-stack.rust.no-algebraic-float-operations
    f32::algebraic_add(left, right)
}

pub fn forbidden_f64_associated_operation(left: f64, right: f64) -> f64 {
    // ruleid: la-stack.rust.no-algebraic-float-operations
    f64::algebraic_add(left, right)
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
