# Security Policy

## Supported Versions

Use the latest released crate version or the default branch. Security fixes
are not backported to older versions unless noted in a release.

This crate is pre-1.0 and under active development, so API compatibility and
security support are tied to the current release line.

## Reporting a Vulnerability

Please report vulnerabilities privately using GitHub private vulnerability
reporting:

<https://github.com/acgetchell/la-stack/security/advisories/new>

Do not open a public issue for suspected vulnerabilities.

Include:

- Affected crate version or commit.
- Enabled Cargo features, especially `exact` if exact arithmetic is involved.
- Steps to reproduce, ideally with a minimal Rust example or test.
- Expected and observed behavior.
- Security impact, such as a panic, denial of service, or incorrect result.
- Any relevant matrix, vector, or benchmark input shape, with sensitive data
  removed.
- A suggested fix or mitigation, if available.

For numerical correctness issues that are not security-sensitive, open a
normal GitHub issue with a minimal reproduction.

## Disclosure Process

- Reports are acknowledged as maintainer availability allows.
- The issue is triaged and its severity assessed on a best-effort basis.
- Accepted reports receive updates when there is meaningful progress or a
  material change in the assessment.
- For an accepted vulnerability, the project prepares a fix, publishes a
  GitHub Security Advisory, releases the fix, and requests a RustSec advisory
  when appropriate.
- If a report is declined, the reporter receives an explanation.

Please follow coordinated disclosure and avoid public disclosure until a fix
or mitigation is available.

## Scope

The crate uses `#![forbid(unsafe_code)]`, which reduces memory-safety risk.
Security-relevant correctness and availability issues can still exist. In
scope are:

- Panics or crashes triggered by malformed or adversarial matrices or vectors.
- CPU or memory denial of service caused by crafted inputs, including inputs
  to optional exact-arithmetic paths.
- Incorrect numerical results that affect security, data integrity, or
  availability when processing untrusted input.
- Violations of documented exact-arithmetic guarantees, such as silent
  precision loss, when they have a security impact.

Out of scope are:

- Documented floating-point limitations, conditioning behavior, or rounding
  bounds that do not create a security impact.
- Performance limitations that are not exploitable as denial of service.
- Issues caused by use outside the documented API contracts or supported
  problem scope.

## Patch and Advisory Policy

- Fixes are released on the latest supported release line. Older releases
  receive fixes only when explicitly noted.
- Releases are published to crates.io with corresponding GitHub releases.
- Accepted vulnerabilities are documented with GitHub Security Advisories and,
  when appropriate, RustSec advisories.
- Public technical detail may be limited until users have had a reasonable
  opportunity to update.

## RustSec

Applicable vulnerabilities may be disclosed through the
[RustSec Advisory Database](https://github.com/RustSec/advisory-db), enabling
detection with `cargo audit`.

## Safe Harbor

Good-faith security research is welcome. Avoid privacy violations, data
destruction, persistence, service disruption, and public disclosure before a
fix or mitigation is available. Reports that follow coordinated disclosure
and make a reasonable effort to avoid harm are treated as helpful
contributions.

## Acknowledgements

Responsible disclosure is appreciated. Reporters may be credited in
advisories or release notes unless anonymity is requested.

## Security Checks

This project uses GitHub CodeQL, Dependabot security updates, secret scanning
with push protection, `cargo audit`, zizmor, Clippy SARIF analysis, and
repository-owned Semgrep rules.
