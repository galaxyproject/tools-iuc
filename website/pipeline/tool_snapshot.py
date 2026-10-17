"""Publish tool YAML, index, and metadata together through one directory pointer."""

from __future__ import annotations

import os
import shutil
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

from .common import read_yaml


def load_tool_index(data_dir: Path) -> list[dict]:
    return read_yaml(data_dir / "tools" / "index.yaml") or []


@contextmanager
def tool_snapshot(data_dir: Path) -> Iterator[Path]:
    snapshots = data_dir / ".tool-snapshots"
    snapshots.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(dir=snapshots))
    pointer = data_dir / "tools"
    previous = pointer.resolve() if pointer.exists() else None
    pending = data_dir / ".tools-next"
    try:
        yield stage
        pending.symlink_to(stage.relative_to(data_dir), target_is_directory=True)
        # Migrate the original flat directory once; subsequent switches are
        # a single atomic rename of a symlink to an immutable snapshot.
        if pointer.exists() and not pointer.is_symlink():
            previous = snapshots / f"legacy-{stage.name}"
            pointer.rename(previous)
        try:
            os.replace(pending, pointer)
        except BaseException:
            if previous and not pointer.exists():
                previous.rename(pointer)
            raise
    except BaseException:
        pending.unlink(missing_ok=True)
        shutil.rmtree(stage)
        raise
    # Keep the preceding generation so existing readers can finish using it.
    for old in snapshots.iterdir():
        if old not in {stage, previous}:
            shutil.rmtree(old)
