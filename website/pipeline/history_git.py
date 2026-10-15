"""Batched Git reads for the persistent website history index."""

from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

from .inventory import TOOL_ROOTS

HISTORY_ROOTS = (*TOOL_ROOTS, "macros")
EMAIL_IN_NAME = re.compile(r"[\w.+-]+@[\w.-]+\.[A-Za-z]{2,}")
NOREPLY = re.compile(r"(?:\d+\+)?([^@]+)@users\.noreply\.github\.com", re.I)


def git(repo: Path, *args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=repo, text=True).strip()


@dataclass
class Change:
    status: str
    path: str
    old_path: str | None = None
    tool_id: str | None = None


@dataclass
class Commit:
    sha: str
    parents: str
    name: str
    github: str | None
    authored: str
    date: str
    subject: str
    changes: list[Change] = field(default_factory=list)


def read_commits(repo: Path, head: str, exclude: str | None) -> list[Commit]:
    """Read one history range. Emails are discarded before records are returned."""
    command = [
        "git",
        "log",
        "--reverse",
        "--topo-order",
        "--root",
        "--find-renames",
        "--diff-merges=first-parent",
        "--name-status",
        "-z",
        "--format=%x00HISTORY-COMMIT%x00%H%x00%P%x00%an%x00%ae%x00%aI%x00%cI%x00%s%x00",
        head,
    ]
    if exclude:
        command.extend(["--not", exclude])
    command.extend(["--", *HISTORY_ROOTS])
    raw = subprocess.check_output(command, cwd=repo).decode("utf-8")
    tokens = iter(raw.split("\0"))
    commits: list[Commit] = []
    for token in tokens:
        token = token.lstrip("\n")
        if not token:
            continue
        if token == "HISTORY-COMMIT":
            sha, parents, name, email, authored, date, subject = (next(tokens) for _ in range(7))
            match = NOREPLY.fullmatch(email)
            name = EMAIL_IN_NAME.sub("", name).strip(" <>")
            commits.append(
                Commit(sha, parents, name, match.group(1).lower() if match else None, authored, date, subject)
            )
        else:
            path = next(tokens)
            if token.startswith(("R", "C")):
                commits[-1].changes.append(Change(token, next(tokens), path))
            else:
                commits[-1].changes.append(Change(token, path))
    return commits


class GitReader:
    """Read historical files through one cat-file process, rather than one per file."""

    def __init__(self, repo: Path):
        self.process = subprocess.Popen(
            ["git", "cat-file", "--batch"],
            cwd=repo,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )
        self.read = lru_cache(maxsize=8192)(self._read)

    def _read(self, sha: str, path: str) -> str | None:
        assert self.process.stdin is not None and self.process.stdout is not None
        self.process.stdin.write(f"{sha}:{path}\n".encode())
        self.process.stdin.flush()
        header = self.process.stdout.readline().split()
        if not header:
            raise RuntimeError("git cat-file exited while reading history")
        if header[-1] == b"missing":
            return None
        size = int(header[-1])
        data = self.process.stdout.read(size)
        if len(data) != size or self.process.stdout.read(1) != b"\n":
            raise RuntimeError("Incomplete git cat-file response")
        return data.decode("utf-8", errors="replace")

    def close(self) -> None:
        assert self.process.stdin is not None and self.process.stdout is not None
        self.process.stdin.close()
        self.process.stdout.close()
        self.process.wait()
        self.read.cache_clear()
