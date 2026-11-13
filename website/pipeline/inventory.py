"""The tracked tool files shared by website metadata and history extraction."""

from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

from .common import git_tracked_files
from .repositories import repository_for_tool

TOOL_ROOTS = ("tools", "tool_collections", "data_managers", "suites")


def tool_files(repo: Path, roots: tuple[str, ...] = TOOL_ROOTS) -> list[Path]:
    files = []
    for root in roots:
        for path in git_tracked_files(repo, root):
            relative = path.relative_to(repo)
            if path.suffix != ".xml" or any(p == "test-data" or p.startswith(".") for p in relative.parts):
                continue
            if is_tool_file(path):
                files.append(path)
    return sorted(files)


def is_tool_file(path: Path) -> bool:
    try:
        root = ET.parse(path).getroot()
        return root.tag == "tool" and bool(root.get("id"))
    except (ET.ParseError, OSError):
        return False


def tool_section(path: Path, repo: Path) -> tuple[str, str]:
    parts = path.relative_to(repo).parts
    section = parts[1] if parts[0] in {"tools", "tool_collections"} else parts[0]
    return f"{parts[0]}_{section}", section.replace("_", " ").title()


def published_tool_files(repo: Path, roots: tuple[str, ...] = TOOL_ROOTS) -> list[Path]:
    return [path for path in tool_files(repo, roots) if repository_for_tool(path) is not None]
