"""Real managed Cargo updates against a disposable local sparse registry."""

import hashlib
import io
import json
import shutil
import sys
import tarfile
import threading
import tomllib
from functools import partial
from http.server import HTTPServer, SimpleHTTPRequestHandler
from pathlib import Path
from typing import TYPE_CHECKING

import pytest
from research_repo_tools.process import run_command

if TYPE_CHECKING:
    from collections.abc import Iterator

REPO_ROOT = Path(__file__).resolve().parents[2]


def crate(registry: Path, name: str, version: str) -> dict[str, object]:
    """Create one deterministic registry candidate for native Cargo resolution."""
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz") as archive:
        for path, contents in {
            "Cargo.toml": f'[package]\nname="{name}"\nversion="{version}"\nedition="2024"\n',
            "src/lib.rs": "pub const VALUE: u8 = 1;\n",
        }.items():
            payload = contents.encode()
            member = tarfile.TarInfo(f"{name}-{version}/{path}")
            member.size = len(payload)
            archive.addfile(member, io.BytesIO(payload))
    payload = buffer.getvalue()
    destination = registry / "crates" / name / version / "download"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(payload)
    return {"name": name, "vers": version, "deps": [], "cksum": hashlib.sha256(payload).hexdigest(), "features": {}, "yanked": False}


@pytest.fixture
def consumer(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Iterator[Path]:
    """Keep actual tool declarations and recipes; substitute only registry inputs."""
    root = tmp_path / "cargo consumer"
    root.mkdir()
    for name in ("pyproject.toml", "uv.lock", ".python-version", "rust-toolchain.toml"):
        shutil.copyfile(REPO_ROOT / name, root / name)
    # Python upgrades have their own real-uv integration test. The aggregate still
    # executes the actual Cargo adapter here without contacting public registries.
    (root / "recipes.just").write_text(
        f"set allow-duplicate-recipes\nimport '{(REPO_ROOT / 'justfile').as_posix()}'\nupdate-python-dependencies:\n",
        encoding="utf-8",
    )
    dependencies = tomllib.loads((REPO_ROOT / "Cargo.toml").read_text(encoding="utf-8"))["dependencies"]
    versions = {name: dependencies[name]["version"] for name in ("num-bigint", "num-rational", "num-traits")}
    registry = tmp_path / "registry"
    registry.mkdir()
    server = HTTPServer(("127.0.0.1", 0), partial(SimpleHTTPRequestHandler, directory=str(registry)))
    host = f"http://127.0.0.1:{server.server_port}"
    index = registry / "index"
    index.mkdir()
    (index / "config.json").write_text(json.dumps({"dl": host + "/crates/{crate}/{version}/download"}), encoding="utf-8")
    for name, version in versions.items():
        candidates = [version, "99.0.0"]
        path = index / name[:2] / name[2:4] / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("".join(json.dumps(crate(registry, name, candidate)) + "\n" for candidate in candidates), encoding="utf-8")
    (root / ".cargo").mkdir()
    (root / ".cargo/config.toml").write_text(f'[registries.fixture]\nindex="sparse+{host}/index/"\n', encoding="utf-8")
    (root / "src").mkdir()
    (root / "src/lib.rs").write_text("pub const VALUE: u8 = 1;\n", encoding="utf-8")
    requirements = "".join(f'{name} = {{ version="{version}", registry="fixture" }}\n' for name, version in versions.items())
    (root / "Cargo.toml").write_text('[package]\nname="update-fixture"\nversion="0.1.0"\nedition="2024"\n[dependencies]\n' + requirements, encoding="utf-8")
    monkeypatch.setenv("UV_NO_SYNC", "1")
    monkeypatch.setenv("UV_PROJECT_ENVIRONMENT", sys.prefix)
    monkeypatch.delenv("CARGO_NET_OFFLINE", raising=False)
    monkeypatch.delenv("CARGO_REGISTRIES_FIXTURE_INDEX", raising=False)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield root
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


@pytest.mark.parametrize("recipe", ["update-cargo-dependencies", "update-dependencies"])
def test_native_updates_advance_eligible_dependencies_and_preserve_coupled_requirements(consumer: Path, recipe: str) -> None:
    """Cargo upgrade and cargo update run through the installed checked toolchain."""
    before = tomllib.loads((consumer / "Cargo.toml").read_text(encoding="utf-8"))["dependencies"]
    result = run_command(
        "just",
        ["--justfile", str(consumer / "recipes.just"), "--working-directory", str(consumer), recipe],
        cwd=consumer,
        timeout=120,
    )
    assert result.returncode == 0
    after = tomllib.loads((consumer / "Cargo.toml").read_text(encoding="utf-8"))["dependencies"]
    assert after["num-bigint"] == before["num-bigint"]
    assert after["num-rational"] == before["num-rational"]
    assert after["num-traits"]["version"] == "99.0.0"
    lock = tomllib.loads((consumer / "Cargo.lock").read_text(encoding="utf-8"))
    assert next(package["version"] for package in lock["package"] if package["name"] == "num-traits") == "99.0.0"
