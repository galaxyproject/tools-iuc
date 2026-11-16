"""Daily full pull of the ToolShed repository index -> data/toolshed_index.yaml.

Used by the /request/ form's client-side duplicate check (the TS API has no CORS
headers, so the browser can't query it directly). One batched paginated read per
day is negligible load on the TS.
"""

from __future__ import annotations

from typing import Any

from bioblend.toolshed import ToolShedInstance

from .common import DATA_DIR, step, step_done, write_yaml

TS_URL = "https://toolshed.g2.bx.psu.edu"


def _slim_repo(r: dict[str, Any]) -> dict[str, Any]:
    return {
        "name": r.get("name"),
        "owner": r.get("owner"),
        "description": (r.get("description") or "")[:300],
        "remote_repository_url": r.get("remote_repository_url"),
        "homepage_url": r.get("homepage_url"),
        "type": r.get("type"),
        "times_downloaded": r.get("times_downloaded"),
    }


def extract_toolshed() -> int:
    step("Pulling ToolShed repository index")
    ts = ToolShedInstance(url=TS_URL)
    # get_repositories is paginated server-side; bioblend fetches all.
    repos = ts.repositories.get_repositories()
    slim = [_slim_repo(r) for r in repos if r.get("name") and r.get("owner")]
    slim.sort(key=lambda r: (r["owner"] or "", r["name"] or ""))
    write_yaml(DATA_DIR / "toolshed_index.yaml", {"count": len(slim), "url": TS_URL, "repos": slim})
    step_done("Pulling ToolShed repository index", f"{len(slim)} repositories")
    return len(slim)


def main() -> None:
    extract_toolshed()
    print(f"  → {DATA_DIR / 'toolshed_index.yaml'}", flush=True)


if __name__ == "__main__":
    main()
