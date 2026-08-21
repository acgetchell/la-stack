"""Stable content identity for release benchmark measurements."""

import hashlib
from typing import TYPE_CHECKING, Final

if TYPE_CHECKING:
    from pathlib import Path

BENCHMARK_CONTRACT_DIRS: Final[tuple[str, ...]] = ("benches",)
BENCHMARK_CONTRACT_FILES: Final[tuple[str, ...]] = (
    ".config/nextest.toml",
    "Cargo.toml",
    "Cargo.lock",
    "rust-toolchain.toml",
    "tests/exact_bench_config.rs",
    "tests/vs_linalg_inputs.rs",
)


def benchmark_contract_files(checkout: Path) -> list[Path]:
    """Return files that determine benchmark code, inputs, and dependencies."""
    files: list[Path] = []
    for relative in BENCHMARK_CONTRACT_FILES:
        path = checkout / relative
        if not path.is_file():
            msg = f"benchmark contract file is missing: {path}"
            raise FileNotFoundError(msg)
        files.append(path)
    for relative in BENCHMARK_CONTRACT_DIRS:
        directory = checkout / relative
        if not directory.is_dir():
            msg = f"benchmark contract directory is missing: {directory}"
            raise FileNotFoundError(msg)
        files.extend(path for path in directory.rglob("*") if path.is_file())
    return sorted(files, key=lambda path: path.relative_to(checkout).as_posix())


def benchmark_contract_digest(checkout: Path) -> str:
    """Hash benchmark code, inputs, dependency resolution, and toolchain."""
    digest = hashlib.sha256()
    for path in benchmark_contract_files(checkout):
        relative = path.relative_to(checkout).as_posix().encode()
        payload = path.read_bytes()
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        digest.update(len(payload).to_bytes(8, "big"))
        digest.update(payload)
    return digest.hexdigest()
