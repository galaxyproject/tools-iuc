"""Clone the research-software-ecosystem/content repo for bio.tools EDAM data.

Galaxy's ``GitContentBiotoolsMetadataSource`` reads ``data/<biotoolsID>/<biotoolsID>.biotools.json``
from a local checkout of this repo. We do a blobless shallow clone (the JSON
files are tiny but the full repo has 172k files / 1.3 GB).

This module also provides a ``build_edam_label_index`` helper that extracts
``{edam_id: label}`` mappings from the same JSON files so the website can
display human-readable EDAM labels alongside URIs.
"""

from __future__ import annotations

import json
import subprocess

from .common import WEBSITE_DIR, step, step_done

REPO_URL = "https://github.com/research-software-ecosystem/content.git"
CLONE_DIR = WEBSITE_DIR / ".rse-content"


def ensure_content_directory() -> str:
    """Blobless shallow clone of rse-content. Returns the path for Galaxy's
    ``GitContentBiotoolsMetadataSource(content_directory=...)``.

    If a cached checkout exists (e.g. from CI cache), refresh it with
    ``git pull`` before use.
    """
    if CLONE_DIR.exists() and (CLONE_DIR / ".git").exists():
        step("Refreshing rse-content checkout")
        subprocess.run(
            ["git", "-C", str(CLONE_DIR), "pull", "--ff-only"],
            check=True, capture_output=True,
        )
        step_done("Refreshing rse-content checkout")
        return str(CLONE_DIR)
    step(f"Cloning rse-content (blobless) → {CLONE_DIR}")
    CLONE_DIR.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["git", "clone", "--filter=blob:none", "--depth=1",
         REPO_URL, str(CLONE_DIR)],
        check=True, capture_output=True,
    )
    step_done("Cloning rse-content")
    return str(CLONE_DIR)


def build_edam_label_index(biotools_ids: list[str]) -> dict[str, dict[str, str]]:
    """Return ``{biotoolsID: {edam_id: label}}`` for label lookup.

    Galaxy's ``expand_ontology_data`` returns EDAM IDs (e.g. ``operation_3198``)
    without human-readable labels. This helper reads the same biotools JSON
    files to build a ``{edam_id: label}`` dict per biotoolsID so the website
    can display labels.
    """
    ids = sorted({i for i in biotools_ids if i})
    if not ids:
        return {}
    result: dict[str, dict[str, str]] = {}
    for bid in ids:
        p = CLONE_DIR / "data" / bid / f"{bid}.biotools.json"
        if not p.exists():
            continue
        try:
            rec = json.loads(p.read_text())
        except Exception:
            continue
        labels: dict[str, str] = {}
        for t in rec.get("topic", []) or []:
            uri = t.get("uri", "")
            term = t.get("term", "")
            if uri and term:
                labels[uri.rsplit("/", 1)[-1]] = term
        for fn in rec.get("function", []) or []:
            for op in fn.get("operation", []) or []:
                uri = op.get("uri", "")
                term = op.get("term", "")
                if uri and term:
                    labels[uri.rsplit("/", 1)[-1]] = term
        result[bid] = labels
    return result
