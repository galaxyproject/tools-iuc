"""Batched Git reads for the persistent website history index."""

from __future__ import annotations

import posixpath
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
        self._object = lru_cache(maxsize=8192)(self._read_object)
        self._tree = lru_cache(maxsize=8192)(self._read_tree)
        self.resolve = lru_cache(maxsize=8192)(self._resolve)
        self.read = lru_cache(maxsize=8192)(self._read)

    def _read(self, sha: str, path: str) -> str | None:
        oid, _paths = self.resolve(sha, path)
        obj = self._object(oid) if oid else None
        return obj[1].decode("utf-8", errors="replace") if obj and obj[0] == "blob" else None

    def dependencies(self, sha: str, path: str) -> tuple[str, ...]:
        """Include both the link paths and their targets at this historical commit."""
        return self.resolve(sha, path)[1]

    def _resolve(self, sha: str, path: str) -> tuple[str | None, tuple[str, ...]]:
        path = posixpath.normpath(path)
        dependencies = {path}
        visited: set[str] = set()
        # Match normal filesystem limits while keeping broken/cyclic links finite.
        for _ in range(41):
            if path.startswith(("../", "/")) or path in visited:
                break
            visited.add(path)
            tree = f"{sha}:"
            parts = path.split("/")
            for number, part in enumerate(parts):
                entry = self._tree(tree).get(part)
                if entry is None:
                    return None, tuple(sorted(dependencies))
                mode, oid = entry
                if mode == "120000":
                    link = "/".join(parts[: number + 1])
                    dependencies.add(link)
                    obj = self._object(oid)
                    assert obj is not None
                    target = obj[1].decode("utf-8", errors="replace")
                    path = posixpath.normpath(posixpath.join(posixpath.dirname(link), target, *parts[number + 1 :]))
                    dependencies.add(path)
                    break
                if number == len(parts) - 1:
                    return (oid if mode in {"100644", "100755"} else None), tuple(sorted(dependencies))
                tree = oid
        return None, tuple(sorted(dependencies))

    def _read_tree(self, spec: str) -> dict[str, tuple[str, str]]:
        obj = self._object(spec)
        if obj is None or obj[0] != "tree":
            return {}
        _kind, data, oid_size = obj
        entries = {}
        offset = 0
        # Git tree entries are mode + space + name + NUL + the binary object ID.
        while offset < len(data):
            space = data.index(b" ", offset)
            nul = data.index(b"\0", space)
            end = nul + 1 + oid_size
            name = data[space + 1 : nul].decode("utf-8", errors="replace")
            entries[name] = (data[offset:space].decode("ascii"), data[nul + 1 : end].hex())
            offset = end
        return entries

    def _read_object(self, spec: str) -> tuple[str, bytes, int] | None:
        assert self.process.stdin is not None and self.process.stdout is not None
        self.process.stdin.write(f"{spec}\n".encode())
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
        return header[1].decode("ascii"), data, len(header[0]) // 2

    def close(self) -> None:
        assert self.process.stdin is not None and self.process.stdout is not None
        self.process.stdin.close()
        self.process.stdout.close()
        self.process.wait()
        self.read.cache_clear()
        self.resolve.cache_clear()
        self._tree.cache_clear()
        self._object.cache_clear()
