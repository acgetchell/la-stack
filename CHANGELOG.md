# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.2] - 2026-03-11

### Added

- Expose determinant error bounds (det_errbound) [`0643bcd`](https://github.com/acgetchell/la-stack/commit/0643bcdd5bd3634c99951acfe68152fd444c2d3f)
- Expose det_errbound without requiring exact feature [`d96f676`](https://github.com/acgetchell/la-stack/commit/d96f676bee8436bde622addcf3ac48c24111da75)

### Maintenance

- Bump actions/setup-node from 6.2.0 to 6.3.0 [`0f2fc34`](https://github.com/acgetchell/la-stack/commit/0f2fc348b06d72476610821c77945c265207ad05)
- Bump taiki-e/install-action from 2.68.20 to 2.68.22 [`c3e49bd`](https://github.com/acgetchell/la-stack/commit/c3e49bdd75e3ca26117942da243dac1ee59eb7a1)

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

[0.2.2]: https://github.com/acgetchell/la-stack/compare/v0.2.1..v0.2.2
[0.2.1]: https://github.com/acgetchell/la-stack/compare/v0.2.0..v0.2.1
[0.2.0]: https://github.com/acgetchell/la-stack/compare/v0.1.3..v0.2.0
[0.1.3]: https://github.com/acgetchell/la-stack/compare/v0.1.2..v0.1.3
[0.1.2]: https://github.com/acgetchell/la-stack/compare/v0.1.1..v0.1.2
[0.1.1]: https://github.com/acgetchell/la-stack/compare/v0.1.0..v0.1.1
[0.1.0]: https://github.com/acgetchell/la-stack/tree/v0.1.0
