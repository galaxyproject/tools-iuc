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
TOOLS_DIR = REPO_ROOT / "tools"
DATA_MANAGERS_DIR = REPO_ROOT / "data_managers"
SUITES_DIR = REPO_ROOT / "suites"

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


def git_last_commit_date(repo: Path, path: str) -> str | None:
    """ISO date of the last commit touching `path`."""
    result = subprocess.run(
        ["git", "log", "-1", "--format=%cI", "--", path],
        cwd=str(repo),
        capture_output=True,
        text=True,
    )
    out = result.stdout.strip()
    return out or None


def git_last_commit_dates_batch(repo: Path, pathspec: str = "tools/") -> dict[str, str]:
    """Batch-fetch last commit date for all files under `pathspec`.

    Returns a dict mapping relative file path → ISO date string.
    Uses a single `git log` call instead of one per file.
    """
    result = subprocess.run(
        ["git", "log", "--format=%cI", "--name-only", "--no-merges", "--", pathspec],
        cwd=str(repo),
        capture_output=True,
        text=True,
    )
    dates: dict[str, str] = {}
    current_date = None
    for line in result.stdout.splitlines():
        if not line:
            # Blank line separates commits — keep current_date for the
            # file paths that follow within the same commit entry.
            continue
        if line[0].isdigit() and "T" in line:
            # This is a date line (ISO format: 2026-08-03T15:23:52+02:00)
            current_date = line
        elif current_date and line not in dates:
            # File path — only set if not already set (git log goes newest-first)
            dates[line] = current_date
    return dates


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


_shed_cache: dict[str, tuple[str, str]] = {}


def owner_repo_from_shed_yml(tool_dir: Path) -> tuple[str, str]:
    """Best-effort (owner, repo) from a .shed.yml next to the tool XML.
    Results are cached per directory (many tools share one .shed.yml)."""
    key = str(tool_dir)
    if key in _shed_cache:
        return _shed_cache[key]
    shed = tool_dir / ".shed.yml"
    if shed.exists():
        try:
            data = read_yaml(shed)
            owner = data.get("owner") or IUC_OWNER
            repo = data.get("name") or tool_dir.name
            result = (owner, repo)
            _shed_cache[key] = result
            return result
        except Exception:
            pass
    result = (IUC_OWNER, tool_dir.name)
    _shed_cache[key] = result
    return result


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
