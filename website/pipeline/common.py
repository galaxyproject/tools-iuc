"""Shared helpers for the IUC website data-extraction pipeline."""

from __future__ import annotations

import hashlib
import os
import subprocess
import time
from pathlib import Path
from typing import Any

import yaml

# Resolve paths relative to the repo root (website/ is two levels up from this file).
WEBSITE_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = WEBSITE_DIR.parent
DATA_DIR = WEBSITE_DIR / "data"

# Default IUC owner on the ToolShed.
IUC_OWNER = "iuc"


def ensure_data_dir(name: str) -> Path:
    d = DATA_DIR / name
    d.mkdir(parents=True, exist_ok=True)
    return d


def write_yaml(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        # allow_unicode so names like "Bérénice" survive; sort_keys off to keep human order
        yaml.safe_dump(data, fh, allow_unicode=True, sort_keys=False, default_flow_style=False, width=120)


def read_yaml(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def git_tracked_files(repo: Path, pathspec: str) -> list[Path]:
    """Return tracked files under `pathspec` via `git ls-files`."""
    result = subprocess.run(
        ["git", "ls-files", "--", pathspec],
        cwd=str(repo),
        capture_output=True,
        text=True,
        check=True,
    )
    return [repo / line for line in result.stdout.splitlines() if line]


def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def md5_files(paths: list[Path]) -> str:
    """Stable MD5 over file paths and contents for small dependency sets."""
    h = hashlib.md5()
    for path in sorted(paths):
        h.update(str(path).encode("utf-8"))
        h.update(b"\0")
        with path.open("rb") as fh:
            for chunk in iter(lambda: fh.read(65536), b""):
                h.update(chunk)
        h.update(b"\0")
    return h.hexdigest()


def git_head_commit(repo: Path = REPO_ROOT) -> str | None:
    """Return the current git commit SHA, or None outside a usable git checkout."""
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=str(repo),
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        return None
    return result.stdout.strip() or None


def env_bool(name: str, default: bool = False) -> bool:
    val = os.environ.get(name)
    if val is None:
        return default
    return val.lower() in ("1", "true", "yes", "on")


# --- status / progress reporting ---------------------------------------------

_start_times: dict[str, float] = {}


def step(name: str) -> None:
    """Announce the start of a pipeline step."""
    _start_times[name] = time.time()
    print(f"→ {name} …", flush=True)


def step_done(name: str, detail: str = "") -> None:
    """Announce completion of a pipeline step with elapsed time."""
    elapsed = time.time() - _start_times.pop(name, time.time())
    msg = f"  ✓ {name} ({elapsed:.1f}s)"
    if detail:
        msg += f" — {detail}"
    print(msg, flush=True)


def progress(current: int, total: int, every: int = 100, label: str = "tools") -> None:
    """Periodic progress counter for long loops. No-op when current % every != 0."""
    if total and (current % every == 0 or current == total):
        pct = current * 100 // total
        print(f"  … {current}/{total} {label} ({pct}%)", flush=True)
