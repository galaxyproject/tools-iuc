"""Resolve tool <creator> entries AND git commit authors into contributor records,
using the GTN ``CONTRIBUTORS.yaml`` as the canonical identity authority.

Sources (in order of authority):
  1. GTN CONTRIBUTORS.yaml (fetched by ``fetch_gtn``) — keyed by GitHub handle,
     contains curated name, ORCID, bio, affiliations, location, socials.
  2. ``<creator>`` elements in tool XML (parsed by extract_tools) — structured
     metadata: name, GitHub URL, ORCID/handle identifier.
  3. Git log over the tools-iuc history — all commit authors with name,
     transient email for GitHub noreply handle extraction, and commit count.

The GTN file is the source of truth for contributor metadata (bio, socials,
affiliations, location).  We do NOT duplicate that data — instead we store
just the ``gtn_handle`` (the key into gtn_people.yaml) on each contributor,
and the frontend joins at build time.

For IUC-specific data that GTN doesn't have (commit counts, tool lists,
member flag, creator-tag-only people), we store it directly.

Matching (no guessing — only exact matches):
  - GitHub handle from noreply email or creator URL → direct identity
  - ORCID (from creator identifier or URL) → GTN orcid map
  - Exact name → GTN name map
  - Git author name that IS a GTN handle (e.g. ``mvdbeek``) → direct identity
  Unmatched entries get their own record keyed by GitHub handle or name slug.

Outputs:
  - data/contributors.yaml   (IUC-specific: tools, commits, member, gtn_handle)
  - data/organisations.yaml  (from <creator class="Organization">)
  - data/members.yaml        (curated IUC members)
"""

from __future__ import annotations

import re
import subprocess
from typing import Any

from .common import DATA_DIR, REPO_ROOT, WEBSITE_DIR, read_yaml, step, step_done, write_yaml
from .fetch_gtn import fetch_gtn_people

# Curated aliases for contributors who changed public GitHub handles.
# Maps old public handle signals to a canonical GitHub handle.
_ALIASES_PATH = WEBSITE_DIR / "config" / "contributor_aliases.yaml"

# ORCID: 0000-0000-0000-0000 (last digit may be X), optionally with URL prefix.
_ORCID_RE = re.compile(
    r"^(?:https?://orcid\.org/)?(\d{4}-\d{4}-\d{4}-\d{3}[0-9X])$", re.I
)
# GitHub profile URL: https://github.com/<handle>
_GITHUB_URL_RE = re.compile(r"https?://github\.com/([^/]+)/?$", re.I)
# GitHub noreply email: <handle>@users.noreply.github.com or <id>+<handle>@...
_GITHUB_NOREPLY_RE = re.compile(r"(?:\d+\+)?([^@]+)@users\.noreply\.github\.com", re.I)
# Handle-like identifier: alphanumeric, dashes, underscores, short.
_HANDLE_RE = re.compile(r"^[a-zA-Z0-9_-]+$")


def _extract_github(url: str | None) -> str | None:
    if not url:
        return None
    m = _GITHUB_URL_RE.match(url.strip().rstrip("/"))
    return m.group(1) if m else None


def _extract_orcid(value: str | None) -> str | None:
    if not value:
        return None
    m = _ORCID_RE.match(str(value).strip())
    return m.group(1) if m else None


def _full_name(creator: dict[str, Any]) -> str:
    return (
        " ".join(filter(None, [creator.get("givenName"), creator.get("familyName")])).strip()
        or creator.get("name")
        or _extract_github(creator.get("url"))
        or "Unknown"
    )


def _name_slug(name: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", name.lower()).strip("-") or "unknown"
    return slug


def _build_gtn_maps(gtn_data: dict[str, Any]) -> dict[str, dict[str, str]]:
    """Build lookup maps from GTN contributor data.

    Each map returns the GTN handle (key) for a given lookup value.
    """
    contribs = gtn_data.get("contributors", {})
    by_handle: dict[str, str] = {}
    by_name: dict[str, str] = {}
    by_orcid: dict[str, str] = {}
    for handle, v in contribs.items():
        h = handle.lower()
        by_handle[h] = handle
        if v.get("name"):
            by_name[v["name"].lower()] = handle
        if v.get("orcid"):
            by_orcid[v["orcid"]] = handle
    return {
        "by_handle": by_handle,
        "by_name": by_name,
        "by_orcid": by_orcid,
    }


def _load_aliases() -> dict[str, dict[str, list[str]]]:
    """Load curated alias mappings from config/contributor_aliases.yaml.

    Returns ``{canonical_handle: {"github": [...], "name": [...]}}``
    with all values lowercased. Returns ``{}`` if the file is absent.
    """
    if not _ALIASES_PATH.exists():
        return {}
    raw = read_yaml(_ALIASES_PATH) or {}
    out: dict[str, dict[str, list[str]]] = {}
    for handle, v in raw.items():
        if not isinstance(v, dict):
            continue
        out[handle.lower()] = {
            "github": [g.lower() for g in (v.get("github") or [])],
            "name": [_name_slug(str(n)) for n in (v.get("name") or [])],
        }
    return out


def _build_alias_maps(
    aliases: dict[str, dict[str, list[str]]],
) -> tuple[dict[str, str], dict[str, str]]:
    """Build reverse-lookup maps from the alias config.

    Returns ``(github_alias→canonical, name_alias→canonical)``.
    """
    gh_map: dict[str, str] = {}
    name_map: dict[str, str] = {}
    for canonical, v in aliases.items():
        for gh in v.get("github", []):
            gh_map[gh] = canonical
        for name in v.get("name", []):
            name_map[name] = canonical
    return gh_map, name_map


def _resolve_key(
    gtn_maps: dict[str, dict[str, str]],
    alias_gh: dict[str, str],
    alias_name: dict[str, str],
    gh: str | None = None,
    orcid: str | None = None,
    name: str | None = None,
) -> tuple[str | None, str | None]:
    """Resolve available signals to a canonical contributor key.

    Tries each signal in order of reliability.  Returns ``(key, gtn_handle)``
    where ``key`` is the canonical identifier (GTN handle, GitHub handle, or
    name slug) and ``gtn_handle`` is the GTN handle if matched (for enrichment),
    or ``None``.
    """
    gh_lower = gh.lower() if gh else None

    # 0. Curated aliases (public handle → canonical handle)
    if gh_lower and gh_lower in alias_gh:
        canonical = alias_gh[gh_lower]
        if canonical in gtn_maps["by_handle"]:
            return gtn_maps["by_handle"][canonical], gtn_maps["by_handle"][canonical]
        return canonical, None
    if name:
        name_alias = _name_slug(name)
        if name_alias in alias_name:
            canonical = alias_name[name_alias]
            if canonical in gtn_maps["by_handle"]:
                return gtn_maps["by_handle"][canonical], gtn_maps["by_handle"][canonical]
            return canonical, None

    # 1. GitHub handle → GTN by_handle
    if gh_lower and gh_lower in gtn_maps["by_handle"]:
        h = gtn_maps["by_handle"][gh_lower]
        return h, h

    # 2. ORCID → GTN by_orcid
    if orcid and orcid in gtn_maps["by_orcid"]:
        h = gtn_maps["by_orcid"][orcid]
        return h, h

    # 3. Name → GTN by_name
    if name and name.lower() in gtn_maps["by_name"]:
        h = gtn_maps["by_name"][name.lower()]
        return h, h

    # 4. Name IS a GTN handle (git authors who commit using their handle)
    if name and name.lower() in gtn_maps["by_handle"]:
        h = gtn_maps["by_handle"][name.lower()]
        return h, h

    # 5. No match — key by GitHub handle (if we have one) or name slug
    if gh_lower:
        return gh_lower, None
    if name:
        return _name_slug(name), None

    return None, None


def _merge_local_key(
    key: str,
    gtn_handle: str | None,
    local_orcid_keys: dict[str, str],
    orcid: str | None = None,
) -> str:
    """Merge unmatched identities by exact stable local identifiers.

    GTN matches and curated aliases already have canonical keys. For everyone
    else, exact ORCID is a strong public identifier and prevents creator-tag
    records from splitting when spelling/accent variants are present.
    """
    if gtn_handle:
        return key

    if orcid and orcid in local_orcid_keys:
        key = local_orcid_keys[orcid]

    if orcid:
        local_orcid_keys.setdefault(orcid, key)
    return key


def _canonical_github(key: str, gh: str | None, name: str | None) -> str | None:
    """Return the public GitHub handle to expose for a resolved contributor."""
    if key != _name_slug(name):
        return key
    return gh



def _git_authors() -> list[dict[str, Any]]:
    """Return [{name, email, commits, first_commit}] from git log over tools/ history.

    Uses a single ``git log`` call with a custom format, then aggregates by
    (name, email). ``first_commit`` is the ISO date of the earliest commit.
    """
    result = subprocess.run(
        ["git", "log", "--format=%cI|%ae|%an", "--date-order", "HEAD", "--", "tools/"],
        cwd=str(REPO_ROOT), capture_output=True, text=True, check=True,
    )
    # date-order: reverse chronological, so the last entry for each author is
    # their earliest commit.
    agg: dict[tuple[str, str], dict[str, Any]] = {}
    for line in result.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split("|", 2)
        if len(parts) != 3:
            continue
        date_str, email, name = parts
        email = email.strip()
        name = name.strip()
        key = (name, email)
        if key not in agg:
            agg[key] = {"name": name, "email": email, "commits": 0, "first_commit": date_str}
        agg[key]["commits"] += 1
        # date-order outputs newest first, so keep updating first_commit —
        # the last value seen for this author is their earliest.
        agg[key]["first_commit"] = date_str
    return list(agg.values())


def extract_people() -> tuple[int, int]:
    step("Extracting contributors & organisations")
    tools_dir = DATA_DIR / "tools"
    if not tools_dir.exists():
        raise SystemExit("Run extract_tools first (data/tools/ missing).")

    # Load GTN as canonical source
    gtn_data = fetch_gtn_people()
    gtn_contribs = gtn_data.get("contributors", {})
    gtn_maps = _build_gtn_maps(gtn_data)

    # Load curated aliases (public handle changes → canonical handle)
    aliases = _load_aliases()
    alias_gh, alias_name = _build_alias_maps(aliases)

    # Load IUC members
    members_path = WEBSITE_DIR / "config" / "iuc_members.yaml"
    members = read_yaml(members_path) if members_path.exists() else []
    member_githubs = {m.get("github", "").lower() for m in (members or []) if m.get("github")}
    # Expand member set with alias canonical handles so that, e.g., erasche
    # being a member also marks hexylena as a member.
    member_githubs |= {alias_gh.get(g, g) for g in member_githubs}

    # Each contributor keyed by canonical key
    contribs: dict[str, dict[str, Any]] = {}
    local_orcid_keys: dict[str, str] = {}

    def _get_contrib(key: str) -> dict[str, Any]:
        if key not in contribs:
            contribs[key] = {
                "id": key,
                "name": "",
                "github": None,
                "orcid": None,
                "url": None,
                "tools": [],
                "commits": 0,
                "first_commit": None,
                # GTN handle — frontend joins gtn_people.yaml for full metadata
                "gtn_handle": None,
            }
        return contribs[key]

    # Source 1: <creator> elements (Persons only)
    step("  Parsing <creator> tags")
    for yml in sorted(tools_dir.glob("**/*.yaml")):
        tool = read_yaml(yml)
        if not tool:
            continue
        tid = tool.get("id")
        for c in tool.get("creators", []) or []:
            cls = c.get("class", "Person")
            if cls == "Organization":
                continue  # handled separately below

            name = _full_name(c)
            gh = _extract_github(c.get("url"))
            orcid = _extract_orcid(c.get("identifier")) or _extract_orcid(c.get("url"))
            url = c.get("url")
            # If URL is an ORCID URL, don't store it as a website
            if url and _extract_orcid(url):
                url = None

            # Handle-like identifier (e.g. "hexylena") — treat as GitHub handle
            ident = c.get("identifier")
            if ident and not orcid:
                s = str(ident).strip()
                if _HANDLE_RE.match(s) and len(s) < 40:
                    if not gh:
                        gh = s

            key, gtn_handle = _resolve_key(gtn_maps, alias_gh, alias_name, gh=gh, orcid=orcid, name=name)
            if not key:
                continue
            key = _merge_local_key(key, gtn_handle, local_orcid_keys, orcid=orcid)

            entry = _get_contrib(key)
            if gtn_handle:
                gtn_entry = gtn_contribs.get(gtn_handle, {})
                entry["name"] = gtn_entry.get("name") or name or entry["name"]
                entry["github"] = gtn_handle
                entry["gtn_handle"] = gtn_handle
            else:
                entry["name"] = name or entry["name"]
                entry["github"] = _canonical_github(key, gh, name) or entry["github"]
                entry["orcid"] = orcid or entry["orcid"]
                entry["url"] = url or entry["url"]

            if tid and tid not in entry["tools"]:
                entry["tools"].append(tid)

    # Source 2: git commit authors (with commit counts + first-commit dates)
    step("  Merging git commit authors")
    for author in _git_authors():
        name = author["name"]
        email = author["email"]
        commits = author["commits"]
        first_commit = author.get("first_commit")

        # Extract GitHub handle from noreply email
        gh = None
        m = _GITHUB_NOREPLY_RE.match(email)
        if m:
            gh = m.group(1)

        key, gtn_handle = _resolve_key(gtn_maps, alias_gh, alias_name, gh=gh, name=name)
        if not key:
            continue
        key = _merge_local_key(key, gtn_handle, local_orcid_keys)

        entry = _get_contrib(key)

        if gtn_handle and not entry["name"]:
            gtn_entry = gtn_contribs.get(gtn_handle, {})
            entry["name"] = gtn_entry.get("name") or name
            entry["github"] = gtn_handle
            entry["gtn_handle"] = gtn_handle
        elif gtn_handle and entry["name"]:
            pass  # already have data from creator tag
        else:
            if not entry["name"]:
                entry["name"] = name
            if not entry["github"]:
                entry["github"] = _canonical_github(key, gh, name)

        entry["commits"] += commits
        if first_commit:
            existing = entry.get("first_commit")
            if not existing or first_commit < existing:
                entry["first_commit"] = first_commit

    # Source 3: Organisations from <creator class="Organization">
    orgs: dict[str, dict[str, Any]] = {}
    for yml in sorted(tools_dir.glob("**/*.yaml")):
        tool = read_yaml(yml)
        if not tool:
            continue
        tid = tool.get("id")
        for c in tool.get("creators", []) or []:
            if c.get("class") != "Organization":
                continue
            name = c.get("name") or "Unknown Organization"
            key = _name_slug(name)
            if key not in orgs:
                orgs[key] = {
                    "id": key,
                    "name": name,
                    "url": c.get("url"),
                    "identifier": c.get("identifier"),
                    "tools": [],
                }
            if c.get("url"):
                orgs[key]["url"] = c.get("url") or orgs[key]["url"]
            if tid and tid not in orgs[key]["tools"]:
                orgs[key]["tools"].append(tid)

    # Build output
    contrib_list = []
    for key, c in sorted(contribs.items(), key=lambda kv: kv[1]["name"].lower()):
        gh = (c.get("github") or "").lower() if c.get("github") else ""
        c["id"] = key
        c["is_member"] = gh in member_githubs
        c["tool_count"] = len(c["tools"])
        contrib_list.append(c)

    org_list = []
    for key, o in sorted(orgs.items(), key=lambda kv: kv[1]["name"].lower()):
        o["tool_count"] = len(o["tools"])
        org_list.append(o)

    write_yaml(DATA_DIR / "contributors.yaml", contrib_list)
    write_yaml(DATA_DIR / "organisations.yaml", org_list)
    n_members = sum(1 for c in contrib_list if c["is_member"])
    n_gtn = sum(1 for c in contrib_list if c.get("gtn_handle"))
    step_done(
        "Extracting contributors & organisations",
        f"{len(contrib_list)} contributors ({n_members} members, {n_gtn} GTN-linked), "
        f"{len(org_list)} organisations",
    )
    return len(contrib_list), len(org_list)


def main() -> None:
    extract_people()
    print(f"  → {DATA_DIR / 'contributors.yaml'}", flush=True)
    print(f"  → {DATA_DIR / 'organisations.yaml'}", flush=True)


if __name__ == "__main__":
    main()
