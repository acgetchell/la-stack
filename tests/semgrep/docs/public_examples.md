# Public example policy fixture

```rust
// ruleid: la-stack.rust.no-unwrap-expect-in-markdown-examples
let value = Some(1_u8).unwrap();
```

```rust
// ok: la-stack.rust.no-unwrap-expect-in-markdown-examples
let value = maybe_value?;
```

<!-- ruleid: la-stack.docs.check-before-fix-command-order -->

```bash
just fix
just check
```
