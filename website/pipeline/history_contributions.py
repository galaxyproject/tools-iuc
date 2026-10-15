"""Project cached public Git facts onto the current published tool inventory."""

from __future__ import annotations

from collections import defaultdict
from pathlib import PurePosixPath

from .history_index import HistoryIndex


def contributor_facts(index: HistoryIndex, tools: list[dict]) -> list[dict]:
    by_id: dict[str, set[str]] = defaultdict(set)
    by_path: dict[str, set[str]] = defaultdict(set)
    by_directory: dict[str, set[str]] = defaultdict(set)
    for tool in tools:
        path = tool["source_path"]
        by_id[tool["id"]].add(path)
        by_path[path].add(path)
        by_directory[tool.get("tool_dir", str(PurePosixPath(path).parent))].add(path)
        for dependency in tool.get("macro_paths", []):
            by_path[dependency].add(path)

    # Follow explicit renames for shared files. Direct wrapper edits carry their
    # historical tool ID, so moving a wrapper does not lose its contributors.
    previous_paths: dict[str, set[str]] = defaultdict(set)
    for row in index.db.execute(
        "SELECT ch.path, ch.old_path FROM changes ch JOIN commits c USING (sha) "
        "WHERE c.active = 1 AND ch.status LIKE 'R%'",
    ):
        previous_paths[row["path"]].add(row["old_path"])
    pending = list(by_path)
    while pending:
        path = pending.pop()
        for old in previous_paths[path]:
            added = by_path[path] - by_path[old]
            if added:
                by_path[old].update(added)
                pending.append(old)

    authors: dict[tuple[str, str | None], dict] = {}
    for row in index.db.execute(
        "SELECT c.sha, c.name, c.github, c.date, ch.path, ch.tool_id FROM commits c "
        "JOIN changes ch USING (sha) WHERE c.active = 1 AND instr(c.parents, ' ') = 0",
    ):
        # Publishing metadata alone does not count as a tool contribution or
        # contributor commit. Filter cached facts here so no reindex is needed.
        if PurePosixPath(row["path"]).name == ".shed.yml":
            continue
        if not row["name"] and not row["github"]:
            continue
        signal = (row["name"] or row["github"], row["github"])
        if signal not in authors:
            authors[signal] = {
                "name": row["name"] or row["github"],
                "github": row["github"],
                "commit_ids": set(),
                "first_commit": row["date"],
                "contributions": {},
            }
        author = authors[signal]
        author["commit_ids"].add(row["sha"])
        author["first_commit"] = min(author["first_commit"], row["date"])
        candidates = set()
        if row["tool_id"]:
            known = by_id.get(row["tool_id"], set())
            candidates.update(known if len(known) == 1 else known.intersection(by_path.get(row["path"], set())))
        else:
            candidates.update(by_path.get(row["path"], set()))
        # Scripts, test data, and documentation maintained beside a wrapper
        # contribute to the tools in that source directory. Other XML only
        # counts when it is an actual imported macro or wrapper.
        if not row["path"].endswith(".xml"):
            for parent in PurePosixPath(row["path"]).parents:
                candidates.update(by_directory.get(str(parent), set()))
        for path in candidates:
            evidence = author["contributions"].setdefault(path, {"commits": set(), "paths": set()})
            evidence["commits"].add(row["sha"])
            evidence["paths"].add(row["path"])
    return sorted(authors.values(), key=lambda author: (author["name"], author["github"] or ""))
