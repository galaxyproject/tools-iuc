"""Fetch which tools are installed on each configured Galaxy server.

For each server in ``config/site.yaml``, calls ``GET /api/tools?in_panel=false``
(a flat list of all installed tools). Extracts the tool short ID (the last
segment of the GUID) and ``tool_shed_repository.{owner, name}`` to build:

  - ``tool_ids``: set of short tool IDs (e.g. ``hicexplorer_hicinfo``)
  - ``repos``: set of ``owner/repo`` strings (e.g. ``bgruening/hicexplorer_hicinfo``)

We match by tool ID first (most reliable — the short ID is the same regardless
of how the ToolShed repo is named), then fall back to ``owner/repo`` matching.

Output: ``data/tool_availability.yaml``

  servers:
    useGalaxy.eu:
      url: https://usegalaxy.eu
      tool_ids: [hicexplorer_hicinfo, abricate, ...]
      repos: [bgruening/hicexplorer_hicinfo, iuc/abricate, ...]
"""

from __future__ import annotations

import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from .common import DATA_DIR, WEBSITE_DIR, read_yaml, step, step_done, write_yaml


def _fetch_installed(server_url: str) -> tuple[set[str], set[str]] | None:
    """Return (tool_ids, repos) for all ToolShed tools on the server.
    tool_ids: short tool IDs (last segment of GUID).
    repos: 'owner/repo' strings.
    Returns None if the API is unreachable."""
    import json
    import urllib.request

    url = f"{server_url.rstrip('/')}/api/tools?in_panel=false"
    try:
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=60) as resp:
            data = json.loads(resp.read())
    except Exception as e:
        print(f"    WARNING: {server_url} — {e}", flush=True)
        return None

    tool_ids: set[str] = set()
    repos: set[str] = set()
    for t in data:
        tid = t.get("id", "")
        # Extract short tool ID from GUID:
        # toolshed.g2.bx.psu.edu/repos/<owner>/<repo>/<tool_id>/<version>
        # → tool_id is the second-to-last segment
        parts = tid.split("/")
        if len(parts) >= 6 and parts[1] == "repos":
            short_id = parts[-2]  # tool_id
            tool_ids.add(short_id)
        tsr = t.get("tool_shed_repository")
        if tsr and tsr.get("owner") and tsr.get("name"):
            repos.add(f"{tsr['owner']}/{tsr['name']}")
    return tool_ids, repos


CACHE_MAX_AGE_DAYS = 7


def _is_cache_fresh(path: Path) -> bool:
    """Check if the cached tool_availability.yaml is younger than CACHE_MAX_AGE_DAYS.
    Uses a ``fetched_at`` timestamp stored inside the YAML rather than file mtime,
    which is not reliable across CI caches, copies, or tar extraction."""
    if not path.exists():
        return False
    cached = read_yaml(path)
    fetched_at = cached.get("fetched_at") if cached else None
    if not fetched_at:
        return False
    try:
        from datetime import datetime
        age = time.time() - datetime.fromisoformat(fetched_at).timestamp()
        return age < CACHE_MAX_AGE_DAYS * 86400
    except Exception:
        return False


def extract_tool_availability() -> dict[str, Any]:
    out_path = DATA_DIR / "tool_availability.yaml"

    # Reuse cached data if fresh (server tool lists don't change often).
    if _is_cache_fresh(out_path):
        step("Fetching tool availability (cached)")
        cached = read_yaml(out_path)
        step_done("Fetching tool availability (cached)", "reusing data < 7 days old")
        return cached

    step("Fetching tool availability from Galaxy servers")
    site = read_yaml(WEBSITE_DIR / "config" / "site.yaml")
    servers = site.get("servers", [])

    import concurrent.futures

    result: dict[str, Any] = {"servers": {}, "fetched_at": datetime.now(UTC).isoformat()}

    def _fetch_one(srv):
        name = srv["name"]
        url = srv["url"]
        print(f"  {name} ({url}) …", flush=True)
        fetched = _fetch_installed(url)
        return name, url, fetched

    with concurrent.futures.ThreadPoolExecutor(max_workers=len(servers)) as pool:
        for name, url, fetched in pool.map(_fetch_one, servers):
            if fetched is not None:
                tool_ids, repos = fetched
                result["servers"][name] = {
                    "url": url,
                    "tool_ids": sorted(tool_ids),
                    "repos": sorted(repos),
                    "tool_count": len(tool_ids),
                    "repo_count": len(repos),
                }
                print(f"    {name}: {len(tool_ids)} tools, {len(repos)} repos", flush=True)
            else:
                print(f"    {name}: skipping (unreachable)", flush=True)

    write_yaml(DATA_DIR / "tool_availability.yaml", result)
    n_ok = len(result["servers"])
    step_done(
        "Fetching tool availability",
        f"{n_ok}/{len(servers)} servers reachable",
    )
    return result


def main() -> None:
    extract_tool_availability()
    print(f"  → {DATA_DIR / 'tool_availability.yaml'}", flush=True)


if __name__ == "__main__":
    main()
