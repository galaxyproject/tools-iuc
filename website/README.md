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
    stats/             Vendored git-history analyzer (from PR #8152)
  astro/               Astro site (pages, components, layouts)
    src/lib/data.ts    YAML loaders + TypeScript types
  data/                Generated YAML (gitignored, rebuilt each run)
  pyproject.toml       Pipeline deps — Galaxy installed from upstream git
```

## Data sources

| Source | What | How |
|--------|------|-----|
| **Tool XML** (`tools/**/*.xml`) | Name, description, version, EDAM, requirements, inputs, outputs, help, creators | Parsed via Galaxy's `get_tool_source` + `source_store` populator |
| **bio.tools EDAM** | EDAM topic/operation labels (human-readable) | [`research-software-ecosystem/content`](https://github.com/research-software-ecosystem/content) repo, read by Galaxy's `GitContentBiotoolsMetadataSource` |
| **Git history** | Tool version events (new/update), classified as initial / wrapper / version / upstream | Vendored PR #8152 analyzer (`pipeline/stats/`) |
| **Git authors** | Commit authors merged with `<creator>` entries | `git shortlog -sne` in `extract_people.py` |
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

> **Note on tool counts:** The git-history analyzer discovers tools by
> regex-scanning XML for `<tool id="...">` at each commit. Tools whose `id`
> attribute uses macro tokens (e.g. `bcftools_@EXECUTABLE@`) have their real
> ID derived from the XML filename. This brings the tracked count to ~1,600
> (up from ~1,400 with the raw regex). The remaining gap vs. the 2,073 current
> tools is due to merge-commit handling in the first-parent history walk —
> some tools added via squash-merges may not appear in the name-status diff.

## IUC members

Active IUC **members** are a curated subset of all contributors, defined in:

```
website/config/iuc_members.yaml
```

Each entry: `name`, `github` (lowercase handle), optional `orcid` / `affiliation`.
Contributors are auto-discovered from `<creator>` tags + git commit history; this
file only flags which contributors are active IUC members (shown with a badge on
the contributors page). Keep alphabetical by GitHub handle.

## Other Galaxy tool repositories

The front page lists other Galaxy tool repositories beyond the IUC. To add one,
edit [`config/tool_repos.yaml`](config/tool_repos.yaml) and
append an entry with `name`, `url`, `description`, and optional `owner`.

## Pipeline stages

```
build_tool_conf       → tool_conf.xml (all tool XML paths)
extract_tools         → data/tools/, data/tools_index.yaml
extract_toolshed      → data/toolshed_index.yaml
extract_tool_availability → data/tool_availability.yaml (which servers have each tool)
extract_stats         → data/stats.yaml, data/recent.yaml
extract_people        → data/contributors.yaml, data/organisations.yaml
extract_releases      → data/releases.yaml
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
- **stats history** — `data/stats/tool_events.json` is reused when its metadata
  says it was generated for the current `HEAD`. Set `IUC_FORCE_EXTRACT_STATS=1`
  to force a full history re-analysis.
- **tool metadata** — `extract_tools` writes `data/tools_metadata.yaml` with a
  signature over `tool_conf.xml`, tool XML files, and local `macros.xml` files.
  Set `IUC_FORCE_EXTRACT_TOOLS=1` to force a full Galaxy parse.
