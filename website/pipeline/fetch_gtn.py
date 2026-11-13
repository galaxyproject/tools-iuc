"""Fetch the Galaxy Training Network (GTN) ``CONTRIBUTORS.yaml`` and
``ORGANISATIONS.yaml`` — the canonical source of truth for contributor
metadata across the Galaxy ecosystem (the Galaxy Hub mirrors these).

Both files are keyed by GitHub handle (userid) and contain curated name, ORCID,
bio, affiliations, location, and social-media links.  We use them as the
identity authority for de-duplicating IUC contributors.

Output: ``data/gtn_people.yaml``

  fetched_at: 2026-08-05T12:00:00+00:00
  contributors:
    <handle>: { name, orcid, bio, affiliations, ... }
  organisations:
    <handle>: { name, url, avatar, ror, ... }

A 7-day cache (via ``fetched_at`` stored inside the YAML) avoids refetching
on every pipeline run, mirroring the ``extract_tool_availability`` approach.
"""

from __future__ import annotations

import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from .common import DATA_DIR, read_yaml, step, step_done, write_yaml

# Raw GitHub URLs for the GTN's canonical contributor / organisation files.
GTN_RAW_BASE = "https://raw.githubusercontent.com/galaxyproject/training-material/main"
GTN_CONTRIBUTORS_URL = f"{GTN_RAW_BASE}/CONTRIBUTORS.yaml"
GTN_ORGANISATIONS_URL = f"{GTN_RAW_BASE}/ORGANISATIONS.yaml"

CACHE_MAX_AGE_DAYS = 7
OUT_PATH = DATA_DIR / "gtn_people.yaml"


def _fetch_yaml(url: str) -> Any:
    """Fetch and parse a YAML document from `url`. Returns None on failure."""
    import urllib.request

    try:
        req = urllib.request.Request(url, headers={"Accept": "text/yaml"})
        with urllib.request.urlopen(req, timeout=60) as resp:
            import yaml

            return yaml.safe_load(resp.read())
    except Exception as e:
        print(f"    WARNING: {url} — {e}", flush=True)
        return None


def _is_cache_fresh(path: Path) -> bool:
    """True when the cached gtn_people.yaml has a fetched_at < 7 days old."""
    if not path.exists():
        return False
    cached = read_yaml(path)
    fetched_at = cached.get("fetched_at") if cached else None
    if not fetched_at:
        return False
    try:
        age = time.time() - datetime.fromisoformat(fetched_at).timestamp()
        return age < CACHE_MAX_AGE_DAYS * 86400
    except Exception:
        return False


def _clean_contributors(raw: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Drop dummy / non-hall-of-fame entries and normalise the handle keys."""
    out: dict[str, dict[str, Any]] = {}
    for handle, v in (raw or {}).items():
        if not isinstance(v, dict):
            continue
        # Skip dummy / hidden entries (gtn-halloffame: 'no' or hub-halloffame: 'no').
        if v.get("gtn-halloffame") == "no" or v.get("hub-halloffame") == "no":
            continue
        # Skip entries without a name (can't display meaningfully).
        if not v.get("name"):
            continue
        cleaned = dict(v)
        cleaned.pop("email", None)
        out[handle] = cleaned
    return out


def _clean_organisations(raw: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Drop dummy / hidden organisation entries."""
    out: dict[str, dict[str, Any]] = {}
    for handle, v in (raw or {}).items():
        if not isinstance(v, dict):
            continue
        if v.get("gtn-halloffame") == "no" or v.get("hub-halloffame") == "no":
            continue
        if not v.get("name"):
            continue
        cleaned = dict(v)
        cleaned.pop("email", None)
        out[handle] = cleaned
    return out


def fetch_gtn_people() -> dict[str, Any]:
    """Return the GTN contributor + organisation data, fetching or using cache.

    The returned dict has the shape::

        { fetched_at: str, contributors: {...}, organisations: {...} }
    """
    if _is_cache_fresh(OUT_PATH):
        step("Fetching GTN people (cached)")
        cached = read_yaml(OUT_PATH)
        cached["contributors"] = _clean_contributors(cached.get("contributors", {}))
        cached["organisations"] = _clean_organisations(cached.get("organisations", {}))
        write_yaml(OUT_PATH, cached)
        n_c = len(cached.get("contributors", {}))
        n_o = len(cached.get("organisations", {}))
        step_done("Fetching GTN people (cached)", f"{n_c} contributors, {n_o} organisations")
        return cached

    step("Fetching GTN contributors & organisations")
    raw_contribs = _fetch_yaml(GTN_CONTRIBUTORS_URL)
    raw_orgs = _fetch_yaml(GTN_ORGANISATIONS_URL)

    if raw_contribs is None:
        # Fall back to existing cache if the fetch failed but we have one.
        if OUT_PATH.exists():
            print("    Using stale cache (fetch failed).", flush=True)
            cached = read_yaml(OUT_PATH)
            cached["contributors"] = _clean_contributors(cached.get("contributors", {}))
            cached["organisations"] = _clean_organisations(cached.get("organisations", {}))
            write_yaml(OUT_PATH, cached)
            return cached
        raise SystemExit(
            "Could not fetch GTN CONTRIBUTORS.yaml and no cache available. "
            f"Check connectivity to {GTN_CONTRIBUTORS_URL}"
        )

    contributors = _clean_contributors(raw_contribs)
    organisations = _clean_organisations(raw_orgs) if raw_orgs else {}

    result: dict[str, Any] = {
        "fetched_at": datetime.now(UTC).isoformat(),
        "contributors": contributors,
        "organisations": organisations,
    }
    write_yaml(OUT_PATH, result)
    step_done(
        "Fetching GTN contributors & organisations",
        f"{len(contributors)} contributors, {len(organisations)} organisations",
    )
    return result


def main() -> None:
    fetch_gtn_people()
    print(f"  → {OUT_PATH}", flush=True)


if __name__ == "__main__":
    main()
