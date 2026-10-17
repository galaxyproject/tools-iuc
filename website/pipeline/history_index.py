"""A rebuildable SQLite checkpoint, persisted exclusively through the CI cache."""

from __future__ import annotations

import json
import os
import sqlite3
import subprocess
from contextlib import closing
from datetime import UTC, datetime
from pathlib import Path

from .common import REPO_ROOT, WEBSITE_DIR, step, step_done
from .history_git import Change, GitReader, git, read_commits
from .history_versions import observe_commit, tool_state
from .stats.version_resolver import ResolvedVersion, classify_change

SCHEMA = 1
HISTORY_PATH = WEBSITE_DIR / ".history-cache" / "history.sqlite"

DDL = """
CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
CREATE TABLE commits (
    sha TEXT PRIMARY KEY, parents TEXT NOT NULL, name TEXT NOT NULL, github TEXT,
    authored TEXT NOT NULL, date TEXT NOT NULL, subject TEXT NOT NULL, active INTEGER NOT NULL
);
CREATE TABLE changes (
    sha TEXT NOT NULL REFERENCES commits(sha), path TEXT NOT NULL, status TEXT NOT NULL,
    old_path TEXT, tool_id TEXT, PRIMARY KEY(sha, path)
);
CREATE TABLE mainline (sha TEXT PRIMARY KEY, position INTEGER UNIQUE NOT NULL);
CREATE TABLE observations (sha TEXT NOT NULL, path TEXT NOT NULL, state TEXT, PRIMARY KEY(sha, path));
CREATE TABLE processed_stats (sha TEXT PRIMARY KEY);
CREATE INDEX changes_path ON changes(path);
"""


class HistoryIndex:
    def __init__(self, repo: Path = REPO_ROOT, path: Path = HISTORY_PATH, *, rebuild: bool = False):
        self.repo = repo
        path.parent.mkdir(parents=True, exist_ok=True)
        self.db = sqlite3.connect(path)
        try:
            version = self.db.execute("PRAGMA user_version").fetchone()[0]
        except sqlite3.DatabaseError:
            version = -1
        if rebuild or version != SCHEMA:
            self.db.close()
            path.unlink()
            self.db = sqlite3.connect(path)
            self.db.executescript(DDL)
            self.db.execute(f"PRAGMA user_version = {SCHEMA}")
        self.db.execute("PRAGMA foreign_keys = ON")
        self.db.row_factory = sqlite3.Row

    def close(self) -> None:
        self.db.close()

    def checkpoint(self) -> str | None:
        row = self.db.execute("SELECT value FROM metadata WHERE key = 'head'").fetchone()
        return row[0] if row else None

    def update(self) -> int:
        head = git(self.repo, "rev-parse", "HEAD")
        previous = self.checkpoint()
        if previous == head:
            print("  History index: unchanged HEAD; no history extraction", flush=True)
            return 0
        if git(self.repo, "rev-parse", "--is-shallow-repository") == "true":
            raise RuntimeError("History indexing requires a full checkout (fetch-depth: 0)")
        step("Updating history index")
        reachable = None
        fast_forward = (
            previous
            and subprocess.run(
                ["git", "merge-base", "--is-ancestor", previous, head],
                cwd=self.repo,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            ).returncode
            == 0
        )
        exclude = previous if fast_forward else None
        if previous and not fast_forward:
            reachable = git(self.repo, "rev-list", "--topo-order", head).splitlines()
            known = {row[0] for row in self.db.execute("SELECT sha FROM commits")}
            exclude = next((sha for sha in reachable if sha in known), None)
        commits = read_commits(self.repo, head, exclude)
        mainline = git(self.repo, "rev-list", "--first-parent", "--reverse", head).splitlines()
        with self.db, closing(GitReader(self.repo)) as reader:
            for number, commit in enumerate(commits, 1):
                self.db.execute(
                    "INSERT INTO commits VALUES (?, ?, ?, ?, ?, ?, ?, 1) ON CONFLICT(sha) DO UPDATE SET active = 1",
                    (
                        commit.sha,
                        commit.parents,
                        commit.name,
                        commit.github,
                        commit.authored,
                        commit.date,
                        commit.subject,
                    ),
                )
                for change in commit.changes:
                    if change.status != "D" and change.path.endswith(".xml") and " " not in commit.parents:
                        state = tool_state(reader, commit.sha, change.path)
                        change.tool_id = state["id"] if state else None
                    self.db.execute(
                        "INSERT OR IGNORE INTO changes VALUES (?, ?, ?, ?, ?)",
                        (commit.sha, change.path, change.status, change.old_path, change.tool_id),
                    )
                if number % 1000 == 0:
                    print(f"  … {number}/{len(commits)} commits indexed", flush=True)
            if reachable is not None:
                self.db.execute("CREATE TEMP TABLE reachable (sha TEXT PRIMARY KEY)")
                self.db.executemany("INSERT INTO reachable VALUES (?)", ((sha,) for sha in reachable))
                self.db.execute(
                    "UPDATE commits SET active = EXISTS(SELECT 1 FROM reachable r WHERE r.sha = commits.sha)"
                )
                self.db.execute("DROP TABLE reachable")
            self.db.execute("DELETE FROM mainline")
            self.db.executemany("INSERT INTO mainline VALUES (?, ?)", ((sha, n) for n, sha in enumerate(mainline)))
            self._update_versions(reader, mainline)
            self.db.execute("INSERT OR REPLACE INTO metadata VALUES ('head', ?)", (head,))
        step_done("Updating history index", f"{len(commits)} commits read; checkpoint {head[:10]}")
        return len(commits)

    def _update_versions(self, reader: GitReader, mainline: list[str]) -> None:
        processed = {row[0] for row in self.db.execute("SELECT sha FROM processed_stats")}
        start = next((n for n, sha in enumerate(mainline) if sha not in processed), len(mainline))
        states: dict[str, dict] = {}
        for row in self.db.execute(
            "SELECT o.path, o.state FROM observations o JOIN mainline m USING (sha) "
            "WHERE m.position < ? ORDER BY m.position, o.path",
            (start,),
        ):
            if row["state"] is None:
                states.pop(row["path"], None)
            else:
                states[row["path"]] = json.loads(row["state"])
        for number, sha in enumerate(mainline[start:], 1):
            if sha in processed:
                observations = {
                    row["path"]: json.loads(row["state"]) if row["state"] else None
                    for row in self.db.execute("SELECT path, state FROM observations WHERE sha = ?", (sha,))
                }
                for path, state in observations.items():
                    if state is None:
                        states.pop(path, None)
                    else:
                        states[path] = state
            else:
                changes = [
                    Change(row["status"], row["path"], row["old_path"], row["tool_id"])
                    for row in self.db.execute("SELECT * FROM changes WHERE sha = ?", (sha,))
                ]
                observations = observe_commit(reader, sha, changes, states)
                self.db.executemany(
                    "INSERT INTO observations VALUES (?, ?, ?)",
                    ((sha, path, json.dumps(state) if state else None) for path, state in observations.items()),
                )
                self.db.execute("INSERT INTO processed_stats VALUES (?)", (sha,))
            if number % 1000 == 0:
                print(f"  … {number}/{len(mainline) - start} first-parent commits observed", flush=True)

    def last_commit_dates(self) -> dict[str, str]:
        dates: dict[str, str] = {}
        for row in self.db.execute(
            "SELECT ch.path, c.date FROM changes ch JOIN commits c USING (sha) "
            "WHERE c.active = 1 AND instr(c.parents, ' ') = 0 "
            "ORDER BY julianday(c.date) DESC, c.sha"
        ):
            dates.setdefault(row["path"], row["date"])
        return dates

    def events(self, tools: list[dict]) -> list[dict]:
        names = {tool["id"]: tool["name"] for tool in tools}
        versions: dict[str, ResolvedVersion] = {}
        events = []
        for row in self.db.execute(
            "SELECT o.path, o.state, c.sha, c.date, c.subject FROM observations o "
            "JOIN mainline m USING (sha) JOIN commits c USING (sha) "
            "WHERE o.state IS NOT NULL ORDER BY m.position, o.path",
        ):
            state = json.loads(row["state"])
            tid = state["id"]
            if tid not in names or state["version"] is None:
                continue
            version = ResolvedVersion(**state["version"])
            previous = versions.get(tid)
            change = classify_change(previous, version) if previous else "initial"
            if change == "no-change":
                continue
            events.append(
                {
                    "tool_id": tid,
                    "tool_name": names[tid],
                    "directory": str(Path(row["path"]).parent),
                    "event_type": "update" if previous else "new",
                    "version": version.full,
                    "old_version": previous.full if previous else None,
                    "commit_hash": row["sha"],
                    "commit_date": datetime.fromisoformat(row["date"]).astimezone(UTC).isoformat(),
                    "change_type": change,
                    "commit_message": row["subject"],
                }
            )
            versions[tid] = version
        return sorted(events, key=lambda event: (event["commit_date"], event["commit_hash"], event["tool_id"]))


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rebuild", action="store_true", help="Discard the cached index and bootstrap again")
    args = parser.parse_args()
    with closing(HistoryIndex(rebuild=args.rebuild)) as index:
        previous = index.checkpoint()
        index.update()
        changed = previous != index.checkpoint()
    if output := os.environ.get("GITHUB_OUTPUT"):
        with open(output, "a", encoding="utf-8") as stream:
            stream.write(f"changed={str(changed).lower()}\n")


if __name__ == "__main__":
    main()
