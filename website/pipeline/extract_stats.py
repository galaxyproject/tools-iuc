"""Compute statistics + recent additions/updates.

Vendors the git-history analyzer from PR #8152 (``pipeline/stats/``) — that PR
will not be merged; this website replaces it. The analyzer walks the full git
history of tools-iuc, resolves version macros at each commit, and classifies
every version change as new / upstream / wrapper / version, emitting
``data/stats/tool_events.json``. We transform that into ``data/stats.yaml`` for
the Vega-Lite charts and ``data/recent.yaml`` for the homepage.

Outputs:
  - data/stats/tool_events.json   (raw events, from the vendored analyzer)
  - data/stats.yaml               (aggregated monthly/cumulative/top-updated)
  - data/recent.yaml              (last N added + updated tool dirs)
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any

from .common import DATA_DIR, REPO_ROOT, env_bool, git_head_commit, read_yaml, step, step_done, write_yaml

RECENT_N = 5


def _stats_metadata_path() -> Path:
    return DATA_DIR / "stats" / "metadata.yaml"


def _cached_events_path() -> Path | None:
    """Return cached tool_events.json when it was generated for the current commit."""
    if env_bool("IUC_FORCE_EXTRACT_STATS"):
        return None
    events_path = DATA_DIR / "stats" / "tool_events.json"
    metadata_path = _stats_metadata_path()
    if not events_path.exists() or not metadata_path.exists():
        return None
    metadata = read_yaml(metadata_path) or {}
    current_commit = git_head_commit(REPO_ROOT)
    if current_commit and metadata.get("git_commit") == current_commit:
        return events_path
    return None


def _run_analyzer() -> Path:
    """Run the vendored #8152 analyzer over the repo history."""
    cached = _cached_events_path()
    if cached is not None:
        step("Analyzing git history (cached)")
        step_done("Analyzing git history (cached)", "reusing events for current HEAD")
        return cached

    step("Analyzing git history (vendored #8152 analyzer)")
    from .stats.main import run  # imported lazily (gitpython dep)

    stats_dir = DATA_DIR / "stats"
    stats_dir.mkdir(parents=True, exist_ok=True)
    import os
    workers = max(1, (os.cpu_count() or 2) - 1)
    events = run(
        repo=str(REPO_ROOT),
        tools_path="tools",
        output_dir=str(stats_dir),
        fmt="json",
        verbose=False,
        workers=workers,
        skip_plots=True,
    )
    # Clean up CSV if it exists from a previous run (we only need JSON)
    csv = stats_dir / "tool_events.csv"
    if csv.exists():
        csv.unlink()
    write_yaml(
        _stats_metadata_path(),
        {
            "git_commit": git_head_commit(REPO_ROOT),
            "event_count": len(events),
        },
    )
    step_done("Analyzing git history")
    return stats_dir / "tool_events.json"


def _build_stats(events_path: Path, index_by_id: dict[str, dict]) -> dict[str, Any]:
    events = json.loads(events_path.read_text(encoding="utf-8"))
    by_month: Counter[str] = Counter()
    new_per_month: Counter[str] = Counter()
    update_per_month: Counter[str] = Counter()
    change_types: Counter[str] = Counter()
    event_types: Counter[str] = Counter()
    first_seen: dict[str, str] = {}
    last_seen: dict[str, str] = {}
    for ev in events:
        date = (ev.get("commit_date") or "")[:7]  # YYYY-MM
        if date:
            by_month[date] += 1
        ct = ev.get("change_type") or "unknown"
        et = ev.get("event_type") or "unknown"
        # Fold "version" (custom/unknown tokens) into "upstream" — custom
        # version tokens are discouraged and should be replaced over time.
        if ct == "version":
            ct = "upstream"
        change_types[ct] += 1
        event_types[et] += 1
        if date and et == "new":
            new_per_month[date] += 1
        elif date:
            update_per_month[date] += 1
        tid = ev.get("tool_id", "")
        d = ev.get("commit_date", "")
        if tid and d:
            if tid not in first_seen or d < first_seen[tid]:
                first_seen[tid] = d
            if tid not in last_seen or d > last_seen[tid]:
                last_seen[tid] = d
    months = sorted(by_month)
    cumulative = []
    running = 0
    for m in months:
        running += by_month[m]
        cumulative.append({"month": m, "events": by_month[m], "cumulative": running})
    # new tools per month (first appearance)
    first_seen_per_month: Counter[str] = Counter()
    for tid, d in first_seen.items():
        first_seen_per_month[d[:7]] += 1
    cum_tools = []
    n = 0
    for m in months:
        n += first_seen_per_month.get(m, 0)
        cum_tools.append({"month": m, "new": first_seen_per_month.get(m, 0), "cumulative": n})
    # top updated by tool directory (suite), not individual tool.
    # Also count unique tools per directory for normalization.
    dir_update_counts: Counter[str] = Counter()
    dir_tool_counts: Counter[str] = Counter()
    for ev in events:
        d = ev.get("directory") or ev.get("tool_id", "").rsplit("_", 1)[0]
        tid = ev.get("tool_id", "")
        if tid:
            dir_tool_counts[d] += 1 if ev.get("event_type") == "new" else 0
        if ev.get("event_type") != "new":
            dir_update_counts[d] += 1
    # Normalize: updates per tool in the directory
    top_updated_suites = [
        {"directory": d, "updates": c, "tools": dir_tool_counts.get(d, 1),
         "updates_per_tool": round(c / max(dir_tool_counts.get(d, 1), 1), 1)}
        for d, c in dir_update_counts.most_common(15)
    ]
    return {
        "source": "pipeline/stats (vendored from #8152)",
        "total_events": len(events),
        "total_tools": len(first_seen),
        "by_month": [{"month": m, "count": by_month[m],
                       "new": new_per_month.get(m, 0), "updates": update_per_month.get(m, 0)} for m in months],
        "cumulative": cumulative,
        "cumulative_tools": cum_tools,
        "change_types": dict(change_types),
        "event_types": dict(event_types),
        "top_updated_suites": top_updated_suites,
    }


def _recent_from_events(index_by_path: dict[str, dict], index_by_id: dict[str, dict]) -> dict[str, Any]:
    """Last N added + last N updated, derived from the events + index."""
    events_path = DATA_DIR / "stats" / "tool_events.json"
    if not events_path.exists():
        return {"added": [], "updated": []}
    events = json.loads(events_path.read_text(encoding="utf-8"))
    # sort by commit_date desc
    events_sorted = sorted(events, key=lambda e: e.get("commit_date", ""), reverse=True)
    added: list[dict[str, Any]] = []
    updated: list[dict[str, Any]] = []
    seen_add: set[str] = set()
    seen_upd: set[str] = set()
    for ev in events_sorted:
        tid = ev.get("tool_id", "")
        entry = index_by_id.get(tid)
        if not entry:
            continue
        rec = {
            "date": (ev.get("commit_date") or "")[:10],
            "tool": entry,
            "version": ev.get("version"),
            "subject": ev.get("commit_message", "")[:120],
        }
        if ev.get("event_type") == "new" and tid not in seen_add:
            seen_add.add(tid)
            added.append(rec)
        elif ev.get("event_type") != "new" and tid not in seen_upd:
            seen_upd.add(tid)
            updated.append(rec)
        if len(added) >= RECENT_N and len(updated) >= RECENT_N:
            break
    return {"added": added[:RECENT_N], "updated": updated[:RECENT_N]}


def extract_stats() -> dict[str, Any]:
    step("Aggregating statistics")
    index_path = DATA_DIR / "tools_index.yaml"
    index = read_yaml(index_path) if index_path.exists() else []
    index_by_id = {r["id"]: r for r in (index or [])}
    index_by_path = {r["path"]: r for r in (index or [])}

    events_path = _run_analyzer()
    stats = _build_stats(events_path, index_by_id)
    stats["tool_count"] = len(index or [])
    write_yaml(DATA_DIR / "stats.yaml", stats)

    recent = _recent_from_events(index_by_path, index_by_id)
    write_yaml(DATA_DIR / "recent.yaml", recent)
    detail = (
        f"{stats['total_events']} events, {stats['total_tools']} tools, "
        f"{len(recent['added'])} added / {len(recent['updated'])} updated recently"
    )
    step_done("Aggregating statistics", detail)
    return stats


def main() -> None:
    extract_stats()
    print(f"  → {DATA_DIR / 'stats.yaml'}", flush=True)
    print(f"  → {DATA_DIR / 'recent.yaml'}", flush=True)


if __name__ == "__main__":
    main()
