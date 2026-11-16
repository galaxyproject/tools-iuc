# IUC Website

Static website for the Galaxy [Intergalactic Utilities Commission](https://github.com/galaxyproject/tools-iuc),
built with [Astro](https://astro.build) + Tailwind. A Python pipeline extracts
metadata from the IUC tool wrappers; Astro renders it into ~2,600 pages.

## Quick start

```bash
# from the repo root
just website-setup    # install Python deps (incl. Galaxy from git) + Node deps
just website-build    # extract all data + build the site + search index
just website-serve    # preview at http://localhost:4321/tools-iuc/
```

Requires [`just`](https://github.com/casey/just) and [`uv`](https://github.com/astral-sh/uv).
See `website/justfile` for all targets.

## Where things live

```
website/
  config/              Curated YAML config (human-edited)
    site.yaml          Galaxy servers, ToolShed URLs
    tool_repos.yaml    Other Galaxy tool repositories (front-page listing)
    iuc_members.yaml   Curated list of active IUC members (see below)
  pipeline/            Python data-extraction pipeline (code only)
    history_index.py   Incremental history checkpoint (SQLite in CI cache)
    stats/             Version classification and standalone history utilities
  astro/               Astro site (pages, components, layouts)
    src/lib/data.ts    YAML loaders + TypeScript types
  data/                Generated YAML (gitignored, rebuilt each run)
  pyproject.toml       Pipeline deps — Galaxy installed from upstream git
```

## Website theme

The site imports `@galaxyproject/brand-tokens/theme.css` after Tailwind in
`astro/src/styles/global.css`, sharing Galaxy's palette and font stack with the
Galaxy Hub and IWC. Use the package's `galaxy-primary` and `bay-of-many-*`
utilities for brand colors; links and buttons use the darker shades for text
contrast. Statistics charts read the same colors through CSS custom properties.
Update the package and npm lockfile to adopt future palette releases instead of
copying color values into the website.

## Data sources

| Source | What | How |
|--------|------|-----|
| **Tool XML** (all four tool roots) | Name, description, version, EDAM, requirements, inputs, outputs, help, creators | Parsed via Galaxy's `get_tool_source` + `source_store` populator |
| **bio.tools EDAM** | EDAM topic/operation labels (human-readable) | [`research-software-ecosystem/content`](https://github.com/research-software-ecosystem/content) repo, read by Galaxy's `GitContentBiotoolsMetadataSource` |
| **Publishing metadata** (`.shed.yml`) | Published owner/repository for each tool, including automatic repositories | Planemo configuration and file expansion |
| **Git history** | Version events, original authors, changed paths and explicit renames | Incremental SQLite index (`pipeline/history_index.py`); one batched scan on bootstrap |
| **Git contributors** | Tools contributed to, commit counts and commit/path evidence | Cached public author signals resolved through current GTN metadata and reviewed aliases |
| **ToolShed API** | All 7,700+ repos for duplicate-check in request form | Daily `bioblend` pull → `data/toolshed_index.yaml` |
| **Galaxy server APIs** | Which tools are installed on each public Galaxy server | `GET /api/tools?in_panel=false` per server → `data/tool_availability.yaml` |
| **`<creator>` XML tags** | Structured person/org metadata (name, GitHub, ORCID) | Parsed by Galaxy's `ToolSource.parse_creators()` |

### Event classification (from git history)

- **initial** — first version of a tool (new tool added)
- **upstream** — the underlying software version changed (e.g. samtools 1.17 → 1.18)
- **wrapper** — only the Galaxy wrapper revision changed (e.g. `galaxy0` → `galaxy1`)

The analyzer also detects a **version** type for changes via custom/unknown
tokens (e.g. `@STACKS_VERSION@` instead of the standard `TOOL_VERSION` macro).
These are folded into **upstream** in the stats because custom version tokens
are discouraged — tools should use the standard `TOOL_VERSION` / `VERSION`
macros so version changes are classified correctly. Over time, remaining custom
tokens should be replaced.

The website and history analyzer share the tracked, published tool inventory
under `tools/`, `tool_collections/`, `data_managers/`, and `suites/`. Test data,
untracked XML, and tools excluded by their publishing rules are omitted.
Current tool IDs are macro-expanded by Galaxy. The history walk analyzes `HEAD`,
including first-parent merge diffs, additions, updates, and moves between roots.

## IUC members

Active IUC **members** are a curated subset of all contributors, defined in:

```
website/config/iuc_members.yaml
```

Each entry: `name`, `github` (lowercase handle), optional `orcid` / `affiliation`.
Contributors are auto-discovered from `<creator>` tags + git commit history; this
file only flags which contributors are active IUC members (shown with a badge on
the contributors page). Keep alphabetical by GitHub handle.

## Contributor identity and duplicates

Contributors are merged only by public, stable identity signals:

- GitHub handle from a `<creator>` URL or a GitHub noreply commit address.
- ORCID from a `<creator>` identifier or ORCID URL.
- GTN `CONTRIBUTORS.yaml` handle, ORCID, or exact GTN name.
- Explicit public GitHub-handle and reviewed name-slug aliases in
  [`config/contributor_aliases.yaml`](config/contributor_aliases.yaml).

Email addresses from git history and `<creator>` tags are not written to the
contributor output and are not used as identity keys. GitHub noreply addresses
are parsed transiently only to recover the public GitHub handle. Do not add
normal/private email addresses to `contributor_aliases.yaml`.

To avoid duplicate contributor pages, prefer fixing the upstream public
metadata:

1. Add the contributor's GitHub URL to the relevant `<creator>` XML entry.
2. Add the contributor's ORCID to the relevant `<creator>` XML entry.
3. Add/update the contributor in GTN `CONTRIBUTORS.yaml` with GitHub/ORCID data.
4. If the person changed GitHub handles, add only the old public handle under
   `github:` for the canonical handle in `contributor_aliases.yaml`.
5. If no upstream public identifier is available, a maintainer may add a
   generated name slug under `name:` after personally verifying the identity.
   This is a fallback for legacy git authors or creator records; prefer
   GitHub/ORCID metadata when it can be added upstream.

After changing identity metadata, regenerate contributors and inspect duplicate
name clusters:

```bash
cd website
uv run python -m pipeline.extract_people
python - <<'PY'
import re, unicodedata, yaml
from collections import defaultdict

contributors = yaml.safe_load(open("data/contributors.yaml"))

def normalise(name):
    ascii_name = unicodedata.normalize("NFKD", name or "").encode("ascii", "ignore").decode()
    return re.sub(r"[^a-z0-9]+", "", ascii_name.lower())

groups = defaultdict(list)
for contributor in contributors:
    groups[normalise(contributor.get("name"))].append(contributor)

for key, group in sorted(groups.items()):
    if len(group) > 1:
        print(key, [(c["id"], c.get("github"), c.get("orcid")) for c in group])
PY
```

With the current reviewed alias file, the generated data has 0 duplicate
normalized-name clusters. If the check above prints new clusters, they need
human review before adding any `name:` alias. Do not merge by private email or
by fuzzy name matching.

## Other Galaxy tool repositories

The front page lists other Galaxy tool repositories beyond the IUC. To add one,
edit [`config/tool_repos.yaml`](config/tool_repos.yaml) and
append an entry with `name`, `url`, `description`, and optional `owner`.

## Pipeline stages

```
history_index         → .history-cache/history.sqlite (independent CI cache)
build_tool_conf       → tool_conf.xml (all tool XML paths)
extract_tools         → data/tools/ (snapshot: owner/repo/slug.yaml, index.yaml, metadata.yaml)
extract_toolshed      → data/toolshed_index.yaml
extract_tool_availability → data/tool_availability.yaml (which servers have each tool)
extract_stats         → data/stats.yaml, data/recent.yaml
extract_people        → data/contributors.yaml, data/organisations.yaml, data/creator_links.yaml, data/contributions.yaml
extract_releases      → validates data/recent.yaml (written by extract_stats)
```

## Galaxy packages

Installed from upstream git via `pyproject.toml` (no manual checkout or
`sys.path` hacking):

- `galaxy-tool-util[edam]` — parser, biotools, ontology_data
- `galaxy-app` — `source_store` (populator, index, SQLAlchemy store)

## Caching

- **rse-content** — cloned (blobless) to `website/.rse-content/`, refreshed via
  `git pull` on subsequent runs. CI caches this directory.
- **uv venv** — CI caches the venv keyed by `pyproject.toml` hash.
- **generated data** — `website/data/` and `tool_source_store.sqlite` are
  generated artifacts. CI restores them from cache for scheduled rebuilds and
  validates freshness before reusing expensive outputs.
- **Git history** — `.history-cache/history.sqlite` is restored independently
  with `actions/cache/restore`, updated, closed, and saved with
  `actions/cache/save` before data extraction or Astro runs. Keys use the history
  schema version, runner OS, commit SHA, and run ID/attempt; prefixes restore
  the latest compatible checkpoint within GitHub's branch cache scope. Only
  changed or rebuilt indexes are uploaded, so unchanged daily builds do not
  duplicate the database. Unique save keys also allow a corrupt cache entry
  to be replaced by a later valid snapshot without overwriting immutable keys. A later website
  failure therefore does not discard a completed index update.
- **Tool update dates** — read from the same index; metadata extraction does
  not run a second full-history Git log.
- **Stats** — `data/stats/tool_events.json` is exported from cached version
  observations. Only newly encountered first-parent commits require historical
  XML reads, using one `git cat-file` process. Re-exporting events, aggregating
  statistics, or changing the published tool inventory requires no history scan.
- **tool metadata** — `extract_tools` writes `data/tools/metadata.yaml` with a
  conservative signature over every tracked XML and `.shed.yml` file in the
  repository (including shared macros), `tool_conf.xml`, all pipeline Python
  files and curated YAML, `pyproject.toml`, and `uv.lock`. `HEAD` is also included
  because metadata contains git-derived commit dates. Additions, deletions, and
  working-tree edits invalidate the whole snapshot; no import graph or per-tool
  cache is maintained. Changes to unrelated XML/configuration can therefore
  cause an extra rebuild. Generated data and bytecode do not enter the signature.
  Set `IUC_FORCE_EXTRACT_TOOLS=1` to force a full Galaxy parse.

Tool extraction writes a fresh immutable generation in `data/.tool-snapshots/`
and switches the `data/tools` symlink only after the entire extraction succeeds.
The index and per-tool YAML therefore come from the same generation, and removed
tools cannot survive in contributor metadata. The previous generation is retained
for existing readers; older generations are removed on the next publication.
Legacy flat output migrates on the next successful extraction.

Contributor extraction reads only the published index and writes creator links
keyed by source path and creator position. Astro renders those resolved IDs;
identity matching and aliases live entirely in Python.

Run `just website-check` (or `bash website/check.sh`) for the same Python tests,
Ruff, JavaScript tests, and Astro type check that run in CI. After building the
site and Pagefind index, run `cd website/astro && npm run test:browser` for browser
search checks (`npx playwright install chromium` installs the browser).

## Incremental contributor and version history

The cache is a generated SQLite file, not a hosted database. Git remains the
source of truth. The initial run scans author/path records once and records
first-parent version observations. An unchanged `HEAD` skips history extraction;
normal updates read only the commits since the checkpoint. For rebases and force
pushes, the index reconciles reachable commit IDs, deactivates removed commits,
and reads changes since the closest indexed ancestor. This also works when the
runner no longer has the previous checkpoint's Git object. Updates and the
checkpoint are committed in one transaction.

GitHub may evict the cache. Missing, invalid, or incompatible database formats
bootstrap again. Keep `fetch-depth: 0` in CI. The history schema (`SCHEMA` in
`history_index.py`) and workflow cache-key version must be bumped together when
the stored interpretation changes. `just -f website/justfile history-rebuild`
explicitly rebuilds the index; ordinary `website-clean` retains it.

Profiles distinguish explicit **Authored tools** credits from **Tools contributed
to**. Contributions use historical wrapper IDs (so wrapper moves retain their
contributors), actual imported macro paths, and support files in a tool's source
directory, including scripts, tests, and documentation. Changes to `.shed.yml`
are excluded from contributor commit counts, first-contribution dates, tool
links, and evidence; mixed commits count only through their other changed files.
This filter applies to cached facts without reindexing. Explicit shared-file
renames are followed. Only currently published tools appear in the output. Merge
commits do not double-credit edits already attributed to their original authors.
Each contributor/tool association retains example commit hashes and changed paths
in `data/contributions.yaml`.

The cache stores public author names and GitHub handles; commit email fields are
used transiently to recover GitHub noreply handles and then discarded. Email
addresses embedded in author names are also removed. Identity aliases and GTN
metadata are applied when generating profiles, so correcting an identity does
not require another Git scan. Missing biographies remain missing; the site does
not invent them.
