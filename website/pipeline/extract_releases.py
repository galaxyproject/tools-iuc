"""Ensure data/recent.yaml exists (last N added + updated tools).

extract_stats already writes recent.yaml from the vendored git-history analyzer;
this module is kept as a thin passthrough so the justfile `extract` target reads
naturally. Run extract_stats first.
"""

from __future__ import annotations

from .common import DATA_DIR, read_yaml, step, step_done


def extract_releases() -> dict:
    step("Extracting recent additions & updates")
    path = DATA_DIR / "recent.yaml"
    if not path.exists():
        from .extract_stats import extract_stats
        extract_stats()
    r = read_yaml(path) or {"added": [], "updated": []}
    step_done(
        "Extracting recent additions & updates",
        f"{len(r.get('added', []))} added, {len(r.get('updated', []))} updated",
    )
    return r


def main() -> None:
    extract_releases()
    print(f"  → {DATA_DIR / 'recent.yaml'}", flush=True)


if __name__ == "__main__":
    main()
