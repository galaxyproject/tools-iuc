"""Version observations for only the commits newly introduced to the index."""

from __future__ import annotations

import html
import posixpath
import re
from dataclasses import asdict

from .history_git import Change, GitReader
from .stats.version_resolver import (
    IMPORT_RE,
    TOKEN_PATTERN,
    TOOL_ID_RE,
    TOOL_VERSION_RE,
    ResolvedVersion,
    _parse_macros,
)

NAME_RE = re.compile(r'<tool\s[^>]*?name=["\']([^"\']+)["\']', re.I)


def expand(value: str, tokens: dict[str, str]) -> str:
    for _ in range(10):
        updated = TOKEN_PATTERN.sub(lambda match: tokens.get(match[1], match[0]), value)
        if updated == value:
            break
        value = updated
    return value


def tool_state(reader: GitReader, sha: str, path: str) -> dict | None:
    source = reader.read(sha, path)
    identifier = TOOL_ID_RE.search(source or "")
    if not identifier:
        return None
    dependencies: set[str] = set()
    loaded: set[str] = set()
    tokens: dict[str, str] = {}

    def load(name: str) -> None:
        if name in loaded:
            return
        loaded.add(name)
        text = reader.read(sha, name)
        dependencies.update(reader.dependencies(sha, name))
        if text is None:
            return
        for imported in IMPORT_RE.findall(text):
            nested = posixpath.normpath(posixpath.join(posixpath.dirname(name), imported.strip()))
            if not nested.startswith(("../", "/")):
                load(nested)
        tokens.update(_parse_macros(text))

    load(path)
    tid = expand(identifier[1], tokens)
    if "@" in tid:
        tid = posixpath.splitext(posixpath.basename(path))[0]
    name = NAME_RE.search(source or "")
    version = TOOL_VERSION_RE.search(source or "")
    resolved = None
    if version:
        full = expand(version[1], tokens)
        if "@" not in full:
            upstream = next(
                (
                    expand(tokens[k], tokens)
                    for k in ("TOOL_VERSION", "VERSION")
                    if f"@{k}@" in version[1] and k in tokens
                ),
                None,
            )
            wrapper = next(
                (
                    expand(tokens[k], tokens)
                    for k in ("VERSION_SUFFIX", "WRAPPER_VERSION")
                    if f"@{k}@" in version[1] and k in tokens
                ),
                None,
            )
            resolved = asdict(ResolvedVersion(full, version[1], upstream_token_val=upstream, wrapper_token_val=wrapper))
            resolved.pop("tokens")
    return {
        "id": html.unescape(tid),
        "name": html.unescape(expand(name[1], tokens)) if name else tid,
        "dependencies": sorted(dependencies),
        "version": resolved,
    }


def observe_commit(
    reader: GitReader, sha: str, changes: list[Change], states: dict[str, dict]
) -> dict[str, dict | None]:
    changed = {path for change in changes for path in (change.path, change.old_path) if path}
    candidates = {change.path for change in changes if change.status != "D" and change.path.endswith(".xml")}
    candidates.update(path for path, state in states.items() if changed.intersection(state["dependencies"]))
    observations: dict[str, dict | None] = {}
    for change in changes:
        if change.status == "D" or change.status.startswith("R"):
            removed = change.old_path or change.path
            if removed in states:
                observations[removed] = None
                states.pop(removed)
    for path in sorted(candidates):
        state = tool_state(reader, sha, path)
        if state is not None or path in states:
            observations[path] = state
        if state is None:
            states.pop(path, None)
        else:
            states[path] = state
    return observations
