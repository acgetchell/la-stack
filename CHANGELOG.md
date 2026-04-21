# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.4.1] - 2026-04-21

### Added

- Regression tests for solver and determinant overflow handling [`f763b11`](https://github.com/acgetchell/la-stack/commit/f763b119bcc57276b83f370b0bf7abce654c7eb8)
- Defensive-path test coverage for LU and LDLT solve_vec [`87d426f`](https://github.com/acgetchell/la-stack/commit/87d426fca1627445b804fd26b62fc7d9d4f0ae48)
- Const-ify Lu/Ldlt det + solve_vec and Matrix inf_norm + det_errbound [`81ecb35`](https://github.com/acgetchell/la-stack/commit/81ecb35bdaf159f1f44d1eb24274ecf82c6567d5)
- Fast-filter boundary proptests for exact determinant sign [`6357db3`](https://github.com/acgetchell/la-stack/commit/6357db35c70bca1b93e6bbf9a4fd231913631950)

### Changed

- Report infinite vs finite off-diagonal pairs as asymmetric [`1805779`](https://github.com/acgetchell/la-stack/commit/1805779dbca49183fbfa95c68ec00984966aa551)
- Finalize documentation, benchmarks, and error handling [`0b98d3f`](https://github.com/acgetchell/la-stack/commit/0b98d3f6dbdd74699c318c4744a2b2f9a1b78481)
- Consolidate and expand const-evaluability tests via macros [`f8d80a0`](https://github.com/acgetchell/la-stack/commit/f8d80a01e9913f87e1b19b2ad5ffbc0994e2bfdb)
- Refactor solve_exact to use hybrid Bareiss forward elimination [`ecbbe8a`](https://github.com/acgetchell/la-stack/commit/ecbbe8a571ccaeb9cfedbf0269b8db44d43a5773)
- Polish exact module (Component struct, errors, perf) [`53a5be6`](https://github.com/acgetchell/la-stack/commit/53a5be6abecc0af332398236ed6803ed75564b03)
- Add adversarial-input coverage for exact arithmetic [#80](https://github.com/acgetchell/la-stack/pull/80) [`5bf5815`](https://github.com/acgetchell/la-stack/commit/5bf5815cb165c3b6145c5592420a58085a66efaa)
- Expand exact-arithmetic re-exports and adversarial benchmarks [`b1a491d`](https://github.com/acgetchell/la-stack/commit/b1a491d902ccdaba6f9cd2e6f8e05514b6dfa3de)
- Update AGENTS.md [`1e0648d`](https://github.com/acgetchell/la-stack/commit/1e0648dad3147ef127c775d8969b7cd214a2a6ed)

### Documentation

- Strengthen LDLT symmetry precondition and add is_symmetric API [`1693307`](https://github.com/acgetchell/la-stack/commit/1693307c58e30a804e9b79abc783af99be8dc947)

### Fixed

- Propagate NaN in Matrix::inf_norm [#85](https://github.com/acgetchell/la-stack/pull/85) [`16ffa45`](https://github.com/acgetchell/la-stack/commit/16ffa45ade11b21a179cad8fcecc51d802086a1d)

### Maintenance

- Bump taiki-e/install-action from 2.73.0 to 2.75.7 [`a1d1c1e`](https://github.com/acgetchell/la-stack/commit/a1d1c1edba7d5bf6928de4e8cab86ba853685430)
- Bump actions/download-artifact from 4.3.0 to 8.0.1 [`8772ea2`](https://github.com/acgetchell/la-stack/commit/8772ea2f18cf1b94a21e269b4bc0c8e0a3b650eb)
- Bump rand from 0.9.2 to 0.9.4 [`b91670a`](https://github.com/acgetchell/la-stack/commit/b91670a46f82306dde5433ca7c406a8757e3fc59)
- Bump actions/upload-artifact from 7.0.0 to 7.0.1 [`7579274`](https://github.com/acgetchell/la-stack/commit/7579274d56ab8c46eae6906c07a3c0d0fd41c84f)
- Bump actions/github-script from 8.0.0 to 9.0.0 [`a522abd`](https://github.com/acgetchell/la-stack/commit/a522abdabf211d47f0f52cbfc2a349fa72ce87e2)
- Bump actions/cache from 5.0.4 to 5.0.5 [`b12d95c`](https://github.com/acgetchell/la-stack/commit/b12d95c8ce242ea6014c6a6aaa7845ff501d2e4d)
- Bump actions-rust-lang/setup-rust-toolchain [`c728774`](https://github.com/acgetchell/la-stack/commit/c72877460efba3d598f64ff253d2df8f3695acb1)
- Bump astral-sh/setup-uv from 8.0.0 to 8.1.0 [`7a18733`](https://github.com/acgetchell/la-stack/commit/7a18733ba0f39102678bbe30772fcde10f20e9d8)
- Bump taiki-e/install-action from 2.75.7 to 2.75.18 [`433cfc1`](https://github.com/acgetchell/la-stack/commit/433cfc1b01c9128e93c82cb553aa63d4091bace3)
- Bump MSRV to Rust 1.95 and adopt new stable features [`0ab3c33`](https://github.com/acgetchell/la-stack/commit/0ab3c336074f2b866256fbe5db8a8ec5306d580a)

### Performance

- Integer-only forward elimination for gauss_solve [#72](https://github.com/acgetchell/la-stack/pull/72) [`a1d3bdb`](https://github.com/acgetchell/la-stack/commit/a1d3bdbdb6fcb778a78a2a3d0cb66b79484e1472)

## [0.4.0] - 2026-04-11

### Added

- Add benchmark comparison tool and performance documentation [`d448fc3`](https://github.com/acgetchell/la-stack/commit/d448fc3e3c884d6145580e5a40eb265bf1f83cb8)

### Changed

- Update version to 0.3.0 and refine metadata keywords [`cc0e967`](https://github.com/acgetchell/la-stack/commit/cc0e967aedb4eabbcd98714b16ef67a8290c112e)
- Improve non-finite error reporting and coordinate accuracy [`d4b1452`](https://github.com/acgetchell/la-stack/commit/d4b14524a02b49a1a8e832818913ab875c4e5fe0)
- Refine benchmark comparison reporting and documentation [`e1b5955`](https://github.com/acgetchell/la-stack/commit/e1b5955fb5024232e34e9df9701bb24fb98efa15)
- Expand test coverage for benchmark comparison edge cases [`bced7d9`](https://github.com/acgetchell/la-stack/commit/bced7d988bd2f42cb6bb5af9c54dabdcf787a5fc)
- Update documentation and tests for integer-only Bareiss [`2ee3f05`](https://github.com/acgetchell/la-stack/commit/2ee3f05caecfdf1a23b61257a7465b3bb6d63614)
- Restrict benchmark baselines to main and improve reporting [`9a7caa2`](https://github.com/acgetchell/la-stack/commit/9a7caa241b5f476c2659772f3189468b10fcba2e)

### Maintenance

- Bump astral-sh/setup-uv from 7.3.1 to 7.5.0 [`0505061`](https://github.com/acgetchell/la-stack/commit/050506123b1da9bc79b04e7d8b44d8195d115199)
- Bump taiki-e/install-action from 2.68.25 to 2.68.31 [`dcff501`](https://github.com/acgetchell/la-stack/commit/dcff501354df2e6434fd6f817a7c15b86a40685f)
- Bump actions-rust-lang/setup-rust-toolchain [`4678d2e`](https://github.com/acgetchell/la-stack/commit/4678d2e8b121c3116d638b59b03eca95510888fe)
- Bump actions/cache from 5.0.3 to 5.0.4 [`633c6e2`](https://github.com/acgetchell/la-stack/commit/633c6e234dc4a94744fa433f41bceb84efe1c99c)
- Bump taiki-e/install-action from 2.68.34 to 2.69.6 [`89a306b`](https://github.com/acgetchell/la-stack/commit/89a306b64bf586ec201957953afac0a3d9c13eb7)
- Bump astral-sh/setup-uv from 7.5.0 to 7.6.0 [`676a6ff`](https://github.com/acgetchell/la-stack/commit/676a6ffe3df416ab316daa84119161893eeeadd0)
- Bump codecov/codecov-action from 5.5.2 to 5.5.3 [`2c15dc1`](https://github.com/acgetchell/la-stack/commit/2c15dc1a0c93f75dd23e0d540041bc01f2e3b71a)
- Bump the dependencies group with 2 updates [`23ade13`](https://github.com/acgetchell/la-stack/commit/23ade1343e575edead08ab49a4e361c9e37f5ca0)
- Bump taiki-e/install-action from 2.69.6 to 2.70.1 [`4073197`](https://github.com/acgetchell/la-stack/commit/40731976f465f69bacb0a6dc4b8520fb9bc9c8f4)
- Bump codecov/codecov-action from 5.5.3 to 6.0.0 [`b3e1380`](https://github.com/acgetchell/la-stack/commit/b3e1380e0b8df85478648036eeb1d9b2c79aaac5)
- Bump taiki-e/install-action from 2.70.1 to 2.73.0 [`7de720d`](https://github.com/acgetchell/la-stack/commit/7de720dfb8328d01843cfabd15d23086ee98832b)
- Bump astral-sh/setup-uv from 7.6.0 to 8.0.0 [`af40753`](https://github.com/acgetchell/la-stack/commit/af40753130fac56d48da7fce2f18b11dc391ebe6)
- Add performance regression detection for exact-arithmetic benchmarks [`44bce99`](https://github.com/acgetchell/la-stack/commit/44bce99bcae8f852dbf7500e71fa182730e08bca)
- Release v0.4.0 [`843bb3c`](https://github.com/acgetchell/la-stack/commit/843bb3cfcc88e114bdd5e00d242b605af9b2b434)

### Performance

- Stack-backed storage for exact arithmetic kernels [`b64ac60`](https://github.com/acgetchell/la-stack/commit/b64ac60e286bcdc1c21334ab271e6ffacfadc187)
- Custom f64 → BigRational via IEEE 754 bit decomposition [#63](https://github.com/acgetchell/la-stack/pull/63) [`0a8ce5b`](https://github.com/acgetchell/la-stack/commit/0a8ce5b45e21fa88e4e9832bc6ada38b3ba68eeb)
- Integer-only Bareiss determinant via BigInt [#64](https://github.com/acgetchell/la-stack/pull/64) [`d422b25`](https://github.com/acgetchell/la-stack/commit/d422b251ca86a914522f80285964d4513bca1817)

## [0.3.0] - 2026-03-12

### Added

- Expose exact determinant value (det_exact, det_exact_f64) [`b1a676c`](https://github.com/acgetchell/la-stack/commit/b1a676c8bd8fde98eb9915a9d2e528f1225f46bf)
- [**breaking**] Expose exact determinant value (det_exact, det_exact_f64) [`92ce476`](https://github.com/acgetchell/la-stack/commit/92ce476201d5759766d73221612cd27492bccbe5)
- Exact linear system solve (solve_exact, solve_exact_f64) [`d04fcd3`](https://github.com/acgetchell/la-stack/commit/d04fcd3b24843b7377a34849dbb2e6469bf9fc56)

### Changed

- Add LaError::Overflow, dimension-coverage tests, and doc fixes [`fed0082`](https://github.com/acgetchell/la-stack/commit/fed00823e552c60ec5a98bc74aeb61254f4e2999)
- Add #[non_exhaustive] to LaError [`d068c0e`](https://github.com/acgetchell/la-stack/commit/d068c0e62073a9746a36dab06761605635281481)
- Use det_direct() in exact_det_3x3 example [`5b69fa8`](https://github.com/acgetchell/la-stack/commit/5b69fa8011cd5471af10be5db8fcde64ec9a9dfa)
- Restructure gauss_solve to avoid break/continue [`125b259`](https://github.com/acgetchell/la-stack/commit/125b259608a2c1965aba486d2594ff9fa328609c)

### Documentation

- Clarify pivoting strategy and use idiomatic rev() in exact.rs [`d884581`](https://github.com/acgetchell/la-stack/commit/d8845816b4cc4cb282913f67d7b94b58fbf91bd2)
- Add LDLT example, solve_exact README snippet, and sync examples [`2bddccf`](https://github.com/acgetchell/la-stack/commit/2bddccf3f20156fadb4fee6a4093e76b1b631de1)

### Maintenance

- Release v0.3.0 [`a09bd6f`](https://github.com/acgetchell/la-stack/commit/a09bd6f9f5f5a9cfa97f56e37481601264a210d6)

## [0.2.2] - 2026-03-11

### Added

- Expose determinant error bounds (det_errbound) [`0643bcd`](https://github.com/acgetchell/la-stack/commit/0643bcdd5bd3634c99951acfe68152fd444c2d3f)
- Expose det_errbound without requiring exact feature [`d96f676`](https://github.com/acgetchell/la-stack/commit/d96f676bee8436bde622addcf3ac48c24111da75)

### Maintenance

- Bump actions/setup-node from 6.2.0 to 6.3.0 [`0f2fc34`](https://github.com/acgetchell/la-stack/commit/0f2fc348b06d72476610821c77945c265207ad05)
- Bump taiki-e/install-action from 2.68.20 to 2.68.22 [`c3e49bd`](https://github.com/acgetchell/la-stack/commit/c3e49bdd75e3ca26117942da243dac1ee59eb7a1)
- Release v0.2.2 [`432c7fe`](https://github.com/acgetchell/la-stack/commit/432c7fefea1039f85f4052ab4d12174ae7c2a69b)

## [0.2.1] - 2026-03-08

### Changed

- Add LDLT NonFinite coverage tests; update README example [`d3b7012`](https://github.com/acgetchell/la-stack/commit/d3b7012c80148e319cf34e3e6c3461177bdcd0f5)
- Update changelog generation and release metadata (internal) [`8f97a0f`](https://github.com/acgetchell/la-stack/commit/8f97a0f060dae57550ec73513983c43c79525696)

### Fixed

- [**breaking**] Return Result from det_sign_exact, rename NonFinite field [`717d5cf`](https://github.com/acgetchell/la-stack/commit/717d5cf13af77764e8941004cdf2c153c260690f)

### Maintenance

- Release v0.2.1 [`b4d0286`](https://github.com/acgetchell/la-stack/commit/b4d028677d33b1cb3fad1f6d79b9c8ba00e7a265)

## [0.2.0] - 2026-03-07

### Added

- Add det_sign_exact() behind "exact" feature flag [#33](https://github.com/acgetchell/la-stack/pull/33) [`d45572a`](https://github.com/acgetchell/la-stack/commit/d45572a05488dbcb73b66a4b36a0e0cfc4656ca9)

### Changed

- Rename det benchmark and add black_box to tests [`40cc2b7`](https://github.com/acgetchell/la-stack/commit/40cc2b7bde49cdd0a2c7a247fb045447447d19df)
- Clarify det_direct D=0 support and fix example description [`1308612`](https://github.com/acgetchell/la-stack/commit/130861273dcd9b459a5c1891656e686a6760756c)
- Polish exact feature — docs, CI, example, review fixes [`abfa297`](https://github.com/acgetchell/la-stack/commit/abfa297d1e8e771a60d1cc9a2a9227011876a965)
- Make Codacy analysis non-blocking on pull requests (internal) [`ecd0a30`](https://github.com/acgetchell/la-stack/commit/ecd0a30d1e580ae47f215b811a60110b8f94223a)
- Harden Codacy CI and expand exact feature tests (internal) [`d911fcc`](https://github.com/acgetchell/la-stack/commit/d911fcc51c02807e91372e1db94aaaa196d6ccb9)
- Expand det_sign_exact tests for large matrices (internal) [`2205f12`](https://github.com/acgetchell/la-stack/commit/2205f12a94743c49333c6f6a7228edde21a0d74f)
- Improve release safety and consolidate error handling [`0ef2b94`](https://github.com/acgetchell/la-stack/commit/0ef2b94798f5fa530a09d92d0d117d42a3bc2b8d)
- Update developer docs and optimize release tagging logic [`8aa129a`](https://github.com/acgetchell/la-stack/commit/8aa129ace22919510022a58def30b85c2639551d)
- Improve internal release scripts and expand test coverage [`5eb5064`](https://github.com/acgetchell/la-stack/commit/5eb5064b46b8cbe408ddc958b86bf01dfee1d62c)

### Dependencies

- Bump the dependencies group with 2 updates [`ad4b505`](https://github.com/acgetchell/la-stack/commit/ad4b505ab5a5ec4c2c960a07f4d43553877828cb)
- Bump taiki-e/install-action from 2.67.17 to 2.67.30 [`26f31c4`](https://github.com/acgetchell/la-stack/commit/26f31c42faf4067604e8f04ee17c6db03c63ccf6)
- Bump astral-sh/setup-uv from 7.2.0 to 7.3.0 [`f5a3fdd`](https://github.com/acgetchell/la-stack/commit/f5a3fddb8d17e4ff9b7c0dae130074f7edc12e18)
- Bump taiki-e/install-action from 2.67.30 to 2.68.7 [`4ec8732`](https://github.com/acgetchell/la-stack/commit/4ec87329f6e4a41e3dbe4e3c69c824282312bfeb)
- Bump astral-sh/setup-uv from 7.3.0 to 7.3.1 [`8d584ee`](https://github.com/acgetchell/la-stack/commit/8d584ee25c18496a181e46e60c4d2d9242aa7886)
- Bump taiki-e/install-action from 2.68.7 to 2.68.15 [`62a0f39`](https://github.com/acgetchell/la-stack/commit/62a0f39df68109c9849d98d8b5f199311295c392)
- Bump actions/upload-artifact from 6.0.0 to 7.0.0 [`d69f150`](https://github.com/acgetchell/la-stack/commit/d69f150b30374150786dd8dde493e24c48eec1ba)
- Bump actions-rust-lang/setup-rust-toolchain from 1.15.2 to 1.15.3 [`58c47f8`](https://github.com/acgetchell/la-stack/commit/58c47f804af43903867a39616fb40a4792e48b89)

### Documentation

- Add det_direct example, hoist 4×4 minors, update docs [`159d04e`](https://github.com/acgetchell/la-stack/commit/159d04edfb4c480d6551b970bea3bcae5cdbfff8)

### Fixed

- Harden tag_release and add release docs [`9ed8915`](https://github.com/acgetchell/la-stack/commit/9ed8915e992728ebc5bf53988a97a4de64f8dae3)

### Maintenance

- Replace cspell with typos-cli, bump MSRV to 1.94.0 [`6db12be`](https://github.com/acgetchell/la-stack/commit/6db12bebdc332f93f9663c0c0be8e785ab449fa6)
- Add changelog tooling (git-cliff + tag script) [`5475651`](https://github.com/acgetchell/la-stack/commit/547565142a374d6a7fed432dd6dbaa8b24cf98a1)
- Release v0.2.0 [`922006d`](https://github.com/acgetchell/la-stack/commit/922006d4a6d905bdb66b7eb729e0cb4eca1da70a)

### Performance

- Closed-form determinant for D=1–4 [#27](https://github.com/acgetchell/la-stack/pull/27) [`a8344e5`](https://github.com/acgetchell/la-stack/commit/a8344e5d14a63d3bafd1b4c285e23f71a296fa6a)

## [0.1.3] - 2026-01-31

### Added

- Adds Zenodo DOI badge and identifier [`927b500`](https://github.com/acgetchell/la-stack/commit/927b5007faa2ec51bbc6aab46c95b5dbb1942b81)

### Changed

- Update tooling, Rust version, and dependencies [`7cde40b`](https://github.com/acgetchell/la-stack/commit/7cde40bcfe26f35d237b6baea10ad777e892f0de)
- Enables benchmark feature for criterion [`d701b58`](https://github.com/acgetchell/la-stack/commit/d701b586e89dd3583a77cd7c96512a7e6d6343a9)

### Dependencies

- Bump astral-sh/setup-uv from 7.1.6 to 7.2.0 [`b96c5b8`](https://github.com/acgetchell/la-stack/commit/b96c5b849fb212878ad89832d14a1ca4b54ea5c2)
- Bump actions/setup-node from 6.1.0 to 6.2.0 [`7df54ec`](https://github.com/acgetchell/la-stack/commit/7df54ec3b21321f4f18ab436776ca225cb705c61)
- Bump taiki-e/install-action from 2.65.7 to 2.67.9 [`bcda58a`](https://github.com/acgetchell/la-stack/commit/bcda58a308cdd6122e6fc7d6a07e056af4b34848)
- Bump actions/cache from 5.0.1 to 5.0.2 [`9fbbdb5`](https://github.com/acgetchell/la-stack/commit/9fbbdb5fecfec9a587f1b9b6df72bf024047a022)
- Bump actions/checkout from 6.0.1 to 6.0.2 [`d490ec2`](https://github.com/acgetchell/la-stack/commit/d490ec2fc6cccf8bdeb9065c98e3e6360984a6c1)

### Maintenance

- Release v0.1.3 [`8f26d72`](https://github.com/acgetchell/la-stack/commit/8f26d72083adb1e76a72169b840455736f39390e)

## [0.1.2] - 2026-01-06

### Added

- Add LDLT factorization, citation metadata, and dev workflow updates [`defa9c9`](https://github.com/acgetchell/la-stack/commit/defa9c9a18cefc86dadc0925e40af656eab20acf)

### Changed

- Correctly filters checksum files in CI [`54533c8`](https://github.com/acgetchell/la-stack/commit/54533c83947458fac27f2dfd50009223ad1d053b)
- Adds `prettier` check to `yaml-fix` [`f5de9f2`](https://github.com/acgetchell/la-stack/commit/f5de9f2bf83ab3ba33f510cd6a62dc5436e3f1fe)

### Dependencies

- Bump pastey from 0.2.0 to 0.2.1 in the dependencies group [`ee47f50`](https://github.com/acgetchell/la-stack/commit/ee47f5027cefd27664bd472cfdca26ea0789c8a7)
- Bump astral-sh/setup-uv from 7.1.5 to 7.1.6 [`2ddd8fc`](https://github.com/acgetchell/la-stack/commit/2ddd8fc9c430b5fb245b8cce2567f896f77f5a7b)
- Bump codecov/codecov-action from 5.5.1 to 5.5.2 [`069c9da`](https://github.com/acgetchell/la-stack/commit/069c9dae2a625537cee609fe7f97109013563c43)
- Bump taiki-e/install-action from 2.62.63 to 2.65.0 [`71e7e0a`](https://github.com/acgetchell/la-stack/commit/71e7e0a721b0680e1597659814dd0b907577049d)
- Bump taiki-e/install-action from 2.65.0 to 2.65.12 [`93d02cf`](https://github.com/acgetchell/la-stack/commit/93d02cfb55d347ab88da8b7fd598f23dcec52163)

## [0.1.1] - 2025-12-16

### Added

- Benchmarks and plotting for performance analysis [`d7dd9b2`](https://github.com/acgetchell/la-stack/commit/d7dd9b2d0fbd9e4b9d4615fc011f86c04fcc4740)

### Changed

- Corrects nalgebra time display in benchmark table [`71908f2`](https://github.com/acgetchell/la-stack/commit/71908f25a790782d4a63a0579b6f4703ae59b5bd)
- Refactors benchmarks to compare against other crates [`8224d86`](https://github.com/acgetchell/la-stack/commit/8224d865a782460eef0332971f64041b20a65ba9)
- Correctly escapes backslashes in gnuplot strings [`732480a`](https://github.com/acgetchell/la-stack/commit/732480afe68a14b36ccb546d34ad969ef5b35f41)
- Refactors criterion plotting script for clarity [`be24d03`](https://github.com/acgetchell/la-stack/commit/be24d033ec052fcd25b9598a067c774036675bc2)
- Improves error handling for README table updates [`cd6873f`](https://github.com/acgetchell/la-stack/commit/cd6873f74aaaff4cd272223ad7f0b41f889f5870)
- Updates benchmark data and plot [`f157436`](https://github.com/acgetchell/la-stack/commit/f157436439b684e28e87e7438eec9a0736d94b5b)

### Maintenance

- Release v0.1.1 [`31cd334`](https://github.com/acgetchell/la-stack/commit/31cd3348d5816b2ee9268cd35f5ebd2dd68e8adc)

## [0.1.0] - 2025-12-15

### Added

- Implements CI/CD workflows and code quality tools [`1989c3d`](https://github.com/acgetchell/la-stack/commit/1989c3dc3603b3b7674c978826af2d4318da504a)
- Initial implementation of `la-stack` crate [`c3ab9c5`](https://github.com/acgetchell/la-stack/commit/c3ab9c59d48463a29992a6f3918d22b037059dfa)
- Configures CI workflow for cross-platform builds [`ff96949`](https://github.com/acgetchell/la-stack/commit/ff9694913a52d6c54a643e9186a8d7481d4d8eda)

### Changed

- Initial commit [`e9a3553`](https://github.com/acgetchell/la-stack/commit/e9a3553e823fab736b603c97712f484cfa990f7c)
- Improves linting and validation workflows [`f874ff7`](https://github.com/acgetchell/la-stack/commit/f874ff78552e1f61333e30ca07b20a914ddcc98f)
- Overhauls public API tests and examples [`5699ac5`](https://github.com/acgetchell/la-stack/commit/5699ac59cfe923cd4db7f5fca0420c450b5dc27a)
- Improves numerical stability in tests [`9665d94`](https://github.com/acgetchell/la-stack/commit/9665d946b774376f2e5c5d71310dcc4cdba74847)
- Improves code quality with Codacy and stricter lints [`6e5b004`](https://github.com/acgetchell/la-stack/commit/6e5b004e79c07ed3ca02e89f0af77f911e791b57)
- Improves spell check configuration [`0ad77f6`](https://github.com/acgetchell/la-stack/commit/0ad77f6d1383be1e6ddf885d66c1f8d2c5e8c252)

### Dependencies

- Bump actions/cache from 4.2.4 to 5.0.1 [`a0834f9`](https://github.com/acgetchell/la-stack/commit/a0834f9fe04d76c978ffe1cff5f291939eed65e4)
- Bump actions/upload-artifact from 5.0.0 to 6.0.0 [`9a61b4f`](https://github.com/acgetchell/la-stack/commit/9a61b4f91f5b50b9c34abce3852e8a8fbb9cdee4)

### Maintenance

- Add tarpaulin coverage upload [`7486dfd`](https://github.com/acgetchell/la-stack/commit/7486dfd54e16a6dbde41575c3f35a1acb65f57d2)

[0.4.1]: https://github.com/acgetchell/la-stack/compare/v0.4.0..v0.4.1
[0.4.0]: https://github.com/acgetchell/la-stack/compare/v0.3.0..v0.4.0
[0.3.0]: https://github.com/acgetchell/la-stack/compare/v0.2.2..v0.3.0
[0.2.2]: https://github.com/acgetchell/la-stack/compare/v0.2.1..v0.2.2
[0.2.1]: https://github.com/acgetchell/la-stack/compare/v0.2.0..v0.2.1
[0.2.0]: https://github.com/acgetchell/la-stack/compare/v0.1.3..v0.2.0
[0.1.3]: https://github.com/acgetchell/la-stack/compare/v0.1.2..v0.1.3
[0.1.2]: https://github.com/acgetchell/la-stack/compare/v0.1.1..v0.1.2
[0.1.1]: https://github.com/acgetchell/la-stack/compare/v0.1.0..v0.1.1
[0.1.0]: https://github.com/acgetchell/la-stack/tree/v0.1.0
