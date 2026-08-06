"""Extract per-tool metadata into compact YAML, reusing Galaxy's ``source_store``
infrastructure plus ``galaxy.tool_util.parser.get_tool_source``.

We do NOT write a new XML parser. Flow:

  1. Galaxy packages (``galaxy-app``, ``galaxy-tool-util``) are installed from
     the upstream git repo's ``packages/`` subdirectories via pyproject.toml.
  2. A throwaway SQLite ``SqlAlchemyToolSourceStore`` is populated with one
     ``StoredToolSource`` per tool (macro-expanded source + hash).
  3. ``build_index_entry_from_source`` (the populator's extractor) builds a
     ``ToolIndexEntry`` per tool — the same object Galaxy's toolbox consumes,
     with curated EDAM/xref expansion via ``expand_ontology_data``.
  4. We re-parse a few detail-page-only fields (creators, inputs, outputs,
     rendered help) off the ``ToolSource`` and emit:
       - data/tools/<owner>/<repo>/<id>.yaml   (full per-tool metadata)
       - data/tools_index.yaml                 (slim aggregate)
"""

from __future__ import annotations

import datetime
import html
import re
import warnings
import xml.etree.ElementTree as ET
from html.parser import HTMLParser
from pathlib import Path
from typing import Any
from urllib.parse import quote, urlsplit

# whoosh (a transitive dep of Galaxy's source_store) emits SyntaxWarnings on
# Python 3.12+ for invalid regex escapes in its source. They're harmless and
# come from third-party code we can't fix, so suppress them before importing.
warnings.filterwarnings("ignore", category=SyntaxWarning, module="whoosh")

from galaxy.tool_util.biotools.source import GitContentBiotoolsMetadataSource  # noqa: E402
from galaxy.tool_util.parser import get_tool_source  # noqa: E402
from galaxy.tools.source_store.discover import DiscoveredTool  # noqa: E402
from galaxy.tools.source_store.index import ToolIndex, ToolIndexEntry  # noqa: E402
from galaxy.tools.source_store.interface import StoredToolSource  # noqa: E402
from galaxy.tools.source_store.populator import build_index_entry_from_source, compute_hash  # noqa: E402
from galaxy.tools.source_store.sqlalchemy import SqlAlchemyToolSourceStore  # noqa: E402
from galaxy.util.hash_util import md5_hash_file  # noqa: E402

from .common import (  # noqa: E402
    DATA_DIR,
    REPO_ROOT,
    WEBSITE_DIR,
    ensure_data_dir,
    env_bool,
    git_last_commit_dates_batch,
    md5_files,
    owner_repo_from_shed_yml,
    progress,
    read_yaml,
    step,
    step_done,
    write_yaml,
)
from .rse_content import build_edam_label_index, ensure_content_directory  # noqa: E402

STORE_URL = f"sqlite:///{(WEBSITE_DIR / 'tool_source_store.sqlite').as_posix()}"
TOOLS_METADATA_PATH = DATA_DIR / "tools_metadata.yaml"

# Populated once in extract_tools(); maps biotoolsID → {edam_id: label}.
# Used by _entry_to_full() to add human-readable labels to EDAM terms.
_EDAM_LABELS: dict[str, dict[str, str]] = {}


def _slug_id(tool_id: str) -> str:
    """Convert a tool ID into a filesystem/URL-safe slug.

    Replaces colons, spaces, and any char not in [a-zA-Z0-9_-] with ``_``.
    e.g. ``EMBOSS: antigenic1`` → ``EMBOSS__antigenic1``
    """
    return re.sub(r"[^a-zA-Z0-9_-]+", "_", tool_id)


def _check_slug_collision(
    slug_keys: dict[tuple[str, str, str], str],
    owner: str,
    repo: str,
    tool_id: str,
) -> str:
    slug = _slug_id(tool_id)
    key = (owner, repo, slug)
    existing = slug_keys.get(key)
    if existing and existing != tool_id:
        raise ValueError(
            f"Tool IDs {existing!r} and {tool_id!r} both map to slug {slug!r} in {owner}/{repo}"
        )
    slug_keys[key] = tool_id
    return slug

# Populated once in extract_tools(); maps file path → last commit ISO date.
_COMMIT_DATES: dict[str, str] = {}


def _is_static_help_link(value: str) -> bool:
    value = (value or "").strip()
    if not value:
        return False
    lowered = value.lower()
    if lowered.startswith("#"):
        return True
    scheme = urlsplit(value).scheme.lower()
    if scheme in {"http", "https", "mailto", "ftp"}:
        return True
    return False


def _sanitize_help_url(value: str) -> str | None:
    value = (value or "").strip()
    lowered = value.lower()
    if "$" in value or "@presto_base_url@" in lowered:
        return None
    if lowered.startswith(("root?", "/root?")):
        return None
    if _is_static_help_link(value):
        return value
    return None


class _HelpHtmlSanitizer(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.parts: list[str] = []
        self._skip_object_text_depth = 0

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        formatted, stripped_data = self._format_starttag(tag, attrs, close=False)
        self.parts.append(formatted)
        if tag == "object" and stripped_data:
            self._skip_object_text_depth = 1

    def handle_startendtag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        formatted, _stripped_data = self._format_starttag(tag, attrs, close=True)
        self.parts.append(formatted)

    def handle_endtag(self, tag: str) -> None:
        if self._skip_object_text_depth:
            if tag == "object":
                self._skip_object_text_depth = 0
            else:
                self._skip_object_text_depth += 1
        self.parts.append(f"</{tag}>")

    def handle_data(self, data: str) -> None:
        if self._skip_object_text_depth:
            return
        self.parts.append(html.escape(data, quote=False))

    def handle_comment(self, data: str) -> None:
        self.parts.append(f"<!--{data}-->")

    def handle_entityref(self, name: str) -> None:
        self.parts.append(f"&{name};")

    def handle_charref(self, name: str) -> None:
        self.parts.append(f"&#{name};")

    def _format_starttag(self, tag: str, attrs: list[tuple[str, str | None]], *, close: bool) -> tuple[str, bool]:
        kept_attrs = []
        stripped_data = False
        for name, value in attrs:
            if name in {"data", "href", "src"}:
                value = _sanitize_help_url(value or "")
                if value is None:
                    stripped_data = stripped_data or name == "data"
                    continue
            elif value and "$" in value:
                continue
            if value is None:
                kept_attrs.append(html.escape(name, quote=True))
            else:
                kept_attrs.append(f'{html.escape(name, quote=True)}="{html.escape(value, quote=True)}"')
        suffix = " /" if close else ""
        if kept_attrs:
            return f"<{tag} {' '.join(kept_attrs)}{suffix}>", stripped_data
        return f"<{tag}{suffix}>", stripped_data


def _sanitize_help_html(rendered: str) -> str:
    sanitizer = _HelpHtmlSanitizer()
    sanitizer.feed(rendered)
    sanitizer.close()
    return "".join(sanitizer.parts)


def _render_help(text: str, fmt: str | None) -> str:
    """Render tool help to HTML, matching Galaxy's rendering as closely as possible.

    For RST, Galaxy uses docutils html4css1 writer with doctitle_xform=False and
    a template that yields only the body fragment (see galaxy/util/rst_to_html.py).
    """
    if not text:
        return ""
    fmt = (fmt or "").lower()
    try:
        if fmt == "restructuredtext":
            import docutils.core
            import docutils.utils
            settings_overrides = {
                "embed_stylesheet": False,
                "doctitle_xform": False,
                "halt_level": docutils.utils.Reporter.SEVERE_LEVEL + 1,
                "report_level": 5,  # suppress all reports incl. unknown directives
                "output_encoding": "unicode",
            }
            parts = docutils.core.publish_parts(
                source=text, writer_name="html4css1",
                settings_overrides=settings_overrides,
            )
            return _sanitize_help_html(parts.get("fragment", "") or "")
        elif fmt == "markdown":
            import markdown
            return _sanitize_help_html(markdown.markdown(text, extensions=["tables", "fenced_code"]))
        else:
            return f"<pre>{html.escape(text)}</pre>"
    except Exception:
        return f"<pre>{html.escape(text)}</pre>"


def _compute_container_links(requirements: list[dict[str, Any]]) -> dict[str, Any]:
    """Compute conda, quay.io, and depot.galaxyproject.org links for a tool's
    package requirements. For single-package tools, link directly to the
    package repo. For multi-package tools, compute the mulled-v2 container
    name so we can link to the exact container on quay.io / depot.
    """
    pkg_reqs = [r for r in requirements if r.get("type") == "package" and r.get("name")]
    if not pkg_reqs:
        return {"conda": [], "docker": None, "singularity": None}

    # Conda links — one per package
    conda_links = []
    for r in pkg_reqs:
        name = r["name"]
        version = r.get("version") or ""
        conda_url = f"https://anaconda.org/bioconda/{name}"
        if version:
            conda_url += f"/files?version={version}"
        conda_links.append({"name": name, "version": version, "url": conda_url})

    # For containers:
    # Single-package: link to the package's container repo.
    # Multi-package: compute the mulled-v2 name for the exact combined container.
    if len(pkg_reqs) == 1:
        name = pkg_reqs[0]["name"]
        docker_url = _quay_biocontainer_url(name)
        singularity_url = f"https://depot.galaxyproject.org/singularity/{name}"
    else:
        mulled_name = _mulled_v2_name(pkg_reqs)
        if mulled_name:
            docker_url = _quay_biocontainer_url(mulled_name)
            singularity_url = f"https://depot.galaxyproject.org/singularity/{mulled_name}"
        else:
            docker_url = "https://quay.io/biocontainers/"
            singularity_url = "https://depot.galaxyproject.org/singularity/"

    return {
        "conda": conda_links,
        "docker": docker_url,
        "singularity": singularity_url,
    }


def _quay_biocontainer_url(container_name: str) -> str:
    repository, _, tag = container_name.partition(":")
    url = f"https://quay.io/repository/biocontainers/{repository}"
    if tag:
        url += f"?tab=tags&tag={quote(tag, safe='')}"
    return url


def _mulled_v2_name(pkg_reqs: list[dict[str, Any]]) -> str | None:
    """Compute the mulled-v2 container name for a set of package requirements.

    Uses Galaxy's ``v2_image_name`` which produces names like::
        mulled-v2-<packages_hash>:<versions_hash>

    Returns None if computation fails.
    """
    try:
        from galaxy.tool_util.deps.conda_util import CondaTarget
        from galaxy.tool_util.deps.mulled.util import v2_image_name

        targets = [
            CondaTarget(r["name"], r.get("version") or None)
            for r in sorted(pkg_reqs, key=lambda r: r["name"])
        ]
        return v2_image_name(targets)
    except Exception:
        return None


def _safe(func, *args, default=None, **kw):
    try:
        return func(*args, **kw)
    except Exception:
        return default


def _input_dict(inp) -> dict[str, Any]:
    """Build an input metadata dict from a parsed InputSource."""
    return {
        "name": _safe(inp.parse_name),
        "type": "data",
        "label": _safe(inp.parse_label),
        "help": _safe(inp.parse_help),
        "optional": _safe(inp.parse_optional, default=None),
        "extensions": _safe(inp.parse_extensions, default=[]),
    }


def _detail_fields(ts) -> dict[str, Any]:
    """Re-parse fields the index entry doesn't carry but the detail page needs."""
    creators = _safe(ts.parse_creator, default=[]) or []

    inputs: list[dict[str, Any]] = []

    def _collect_data_inputs(inp):
        """Recursively collect data-type inputs, recursing into conditionals/sections."""
        if _safe(inp.get, "type") == "data":
            inputs.append(_input_dict(inp))
            return
        # Try to recurse into conditional/section nested inputs
        nested_page = _safe(inp.parse_nested_inputs_source)
        if nested_page is not None:
            for n in _safe(nested_page.parse_input_sources, default=[]) or []:
                _collect_data_inputs(n)

    pages = _safe(ts.parse_input_pages)
    if pages is not None:
        for page in pages.page_sources:
            for inp in _safe(page.parse_input_sources, default=[]) or []:
                _collect_data_inputs(inp)

    outputs: list[dict[str, Any]] = []
    outs, _outcol = _safe(ts.parse_outputs, None, default=({}, {}))
    for name, out in (outs or {}).items():
        outputs.append(
            {
                "name": name,
                "format": getattr(out, "format", None),
                "label": getattr(out, "label", None),
                "hidden": getattr(out, "hidden", False),
            }
        )

    help_obj = _safe(ts.parse_help)
    help_text = ""
    help_format = None
    help_html = ""
    if help_obj is not None:
        help_text = getattr(help_obj, "content", "") or ""
        help_format = getattr(help_obj, "format", None)
        help_html = _render_help(help_text, help_format)

    return {
        "creators": creators,
        "inputs": inputs,
        "outputs": outputs,
        "help_html": help_html,
        "help_format": help_format,
        "license": _safe(ts.parse_license),
        "icon": _safe(ts.parse_icon),
    }


def _entry_to_full(entry: ToolIndexEntry, tool_dir: Path, ts) -> dict[str, Any]:
    """Combine a ToolIndexEntry + detail fields into the per-tool YAML dict."""
    detail = _detail_fields(ts)
    src_abs = entry.source_path or ""
    try:
        rel_path = str(Path(src_abs).resolve().relative_to(REPO_ROOT.resolve()))
    except ValueError:
        rel_path = src_abs
    owner, repo = owner_repo_from_shed_yml(tool_dir)
    # biotools / doi from xrefs
    biotools = None
    doi = None
    for x in entry.xrefs or []:
        if x.get("type") == "bio.tools":
            biotools = x.get("value")
        elif x.get("type") == "doi":
            doi = x.get("value")

    # Also check <citation type="doi"> tags (most tools use citations, not xrefs)
    if not doi:
        for c in _safe(ts.parse_citations, default=[]) or []:
            if getattr(c, "type", None) == "doi" and getattr(c, "content", None):
                doi = c.content.strip()
                break

    # EDAM terms from expand_ontology_data are plain IDs like "operation_3198".
    # Convert to {uri, label} objects, looking up labels from the rse-content
    # biotools JSON files.
    bt_labels = _EDAM_LABELS.get(biotools, {}) if biotools else {}

    def _edam_obj(edam_id: str) -> dict[str, str]:
        return {
            "uri": f"http://edamontology.org/{edam_id}",
            "label": bt_labels.get(edam_id, edam_id),
        }

    edam_topics = [_edam_obj(e) for e in (entry.edam_topics or [])]
    edam_operations = [_edam_obj(e) for e in (entry.edam_operations or [])]

    return {
        "id": entry.id,
        "name": entry.name,
        "version": entry.version,
        "description": entry.description,
        "profile": str(entry.profile),
        "tool_type": entry.tool_type,
        "license": detail["license"],
        "hidden": entry.hidden,
        "icon": detail["icon"],
        "edam_operations": edam_operations,
        "edam_topics": edam_topics,
        "xrefs": [dict(x) for x in (entry.xrefs or [])],
        "biotools": biotools,
        "doi": doi,
        "panel_section_id": entry.panel_section_id,
        "panel_section_name": entry.panel_section_name,
        "source_path": rel_path,
        "owner": owner,
        "repo": repo,
        "tool_dir": tool_dir.relative_to(REPO_ROOT).as_posix() if tool_dir.is_relative_to(REPO_ROOT) else str(tool_dir),
        "test_count": entry.test_count,
        "requirements": list(entry.requirements or []),
        "container_requirements": list(entry.container_requirements or []),
        "container_links": _compute_container_links(list(entry.requirements or [])),
        "creators": detail["creators"],
        "help_html": detail["help_html"],
        "help_format": detail["help_format"],
        "inputs": detail["inputs"],
        "outputs": detail["outputs"],
        "tags": list(entry.tags or []),
        "updated": _COMMIT_DATES.get(rel_path),
    }


def _slim(full: dict[str, Any]) -> dict[str, Any]:
    # Collect unique input/output datatypes for filtering
    input_types = sorted({e for inp in full.get("inputs", []) for e in (inp.get("extensions") or [])})
    output_types = sorted({o.get("format") for o in full.get("outputs", []) if o.get("format")})
    return {
        "id": full["id"],
        "slug": full["slug"],
        "name": full["name"],
        "version": full["version"],
        "description": full["description"],
        "edam_operations": full["edam_operations"],
        "edam_topics": full["edam_topics"],
        "biotools": full["biotools"],
        "doi": full["doi"],
        "panel": full["panel_section_name"],
        "panel_id": full["panel_section_id"],
        "owner": full["owner"],
        "repo": full["repo"],
        "path": full["source_path"],
        "tests": full["test_count"],
        "updated": full["updated"],
        "input_types": input_types,
        "output_types": output_types,
    }


def _iter_conf_tools(conf_path: Path):
    """Yield (abs_path, section_id, section_name) from the tool_conf.xml."""
    root = ET.parse(conf_path).getroot()
    for section in root.findall("section"):
        sid = section.get("id")
        sname = section.get("name")
        for tool in section.findall("tool"):
            p = tool.get("file")
            if p:
                yield Path(p), sid, sname
    for tool in root.findall("tool"):
        p = tool.get("file")
        if p:
            yield Path(p), None, None


def _tool_dependency_paths(tool_path: Path) -> list[Path]:
    """Small dependency set used for website extraction cache invalidation."""
    paths = [tool_path]
    macros = tool_path.parent / "macros.xml"
    if macros.exists():
        paths.append(macros)
    return paths


def _tool_dependency_hash(tool_path: Path) -> str:
    return md5_files(_tool_dependency_paths(tool_path))


def _tool_yaml_cache_fresh(path: Path, dependency_hash: str) -> bool:
    if not path.exists():
        return False
    try:
        cached = read_yaml(path) or {}
    except Exception:
        return False
    return cached.get("source_dependency_hash") == dependency_hash


def _tools_manifest_fresh(signature: str) -> bool:
    if env_bool("IUC_FORCE_EXTRACT_TOOLS"):
        return False
    if not TOOLS_METADATA_PATH.exists() or not (DATA_DIR / "tools_index.yaml").exists():
        return False
    try:
        metadata = read_yaml(TOOLS_METADATA_PATH) or {}
    except Exception:
        return False
    return metadata.get("tool_tree_signature") == signature


def _tool_tree_signature(conf_path: Path, conf_tools: list[tuple[Path, str | None, str | None]]) -> str:
    paths = [conf_path]
    for fpath, _sid, _sname in conf_tools:
        if not fpath.is_absolute():
            fpath = REPO_ROOT / fpath
        if fpath.exists():
            paths.extend(_tool_dependency_paths(fpath))
    return md5_files(paths)


def extract_tools(conf_path: Path | None = None) -> tuple[int, int]:
    if conf_path is None:
        conf_path = WEBSITE_DIR / "tool_conf.xml"
    if not conf_path.exists():
        raise SystemExit(f"tool_conf.xml not found at {conf_path}. Run build_tool_conf first.")

    conf_tools = list(_iter_conf_tools(conf_path))
    total = len(conf_tools)
    print(f"  {total} tool XML entries in tool_conf.xml", flush=True)

    signature = _tool_tree_signature(conf_path, conf_tools)
    if _tools_manifest_fresh(signature):
        cached = read_yaml(DATA_DIR / "tools_index.yaml") or []
        step("Building tool source store (cached)")
        step_done("Building tool source store (cached)", f"reusing {len(cached)} tools")
        return len(cached), 0

    step("Building tool source store")
    store = SqlAlchemyToolSourceStore(url=STORE_URL, read_only=False)
    index = ToolIndex()

    tools_out = ensure_data_dir("tools")

    # Pre-scan tool XMLs for bio.tools xrefs, then add legacy biotools
    # mappings from Galaxy's curated TSV (some tools have no <xref> tag but
    # are mapped via biotools_mappings.tsv). Clone the rse-content repo so
    # Galaxy's GitContentBiotoolsMetadataSource can expand EDAM terms. Also
    # build a label index for display.
    import re
    biotools_ids: set[str] = set()
    bt_re = re.compile(r'type="bio\.tools"[^>]*>([^<]+)<')
    for fpath, _sid, _sname in conf_tools:
        if not fpath.is_absolute():
            fpath = REPO_ROOT / fpath
        if not fpath.exists():
            continue
        try:
            text = fpath.read_text(errors="replace")
        except Exception:
            continue
        for m in bt_re.finditer(text):
            biotools_ids.add(m.group(1).strip())
    # Add legacy biotools mappings (tool_id → biotools_id)
    from galaxy.tool_util.ontologies.ontology_data import _biotools_mapping
    legacy = _biotools_mapping()
    for xrefs in legacy.values():
        biotools_ids.update(xrefs)

    biotools_source = None
    if biotools_ids:
        content_dir = ensure_content_directory()
        biotools_source = GitContentBiotoolsMetadataSource(content_dir)
        step(f"Building EDAM label index for {len(biotools_ids)} biotoolsIDs")
        global _EDAM_LABELS
        _EDAM_LABELS = build_edam_label_index(sorted(biotools_ids))
        step_done("Building EDAM label index")

    n = 0
    n_errors = 0
    seen: set[str] = set()
    slug_keys: dict[tuple[str, str, str], str] = {}
    index_rows: list[dict[str, Any]] = []

    # Batch-fetch last commit dates for all tool files (one git log call
    # instead of 2000+ individual subprocess calls).
    step("  Batch-fetching git commit dates")
    global _COMMIT_DATES
    _COMMIT_DATES = git_last_commit_dates_batch(REPO_ROOT, "tools/")
    step_done("  Batch-fetching git commit dates", f"{len(_COMMIT_DATES)} files")

    for i, (fpath, sid, sname) in enumerate(conf_tools, 1):
        if not fpath.is_absolute():
            fpath = REPO_ROOT / fpath
        if not fpath.exists():
            continue
        try:
            ts = get_tool_source(config_file=str(fpath))
            tool_id = ts.parse_id()
            if not tool_id or tool_id in seen:
                continue
            expanded = ts.to_string()
            file_hash = md5_hash_file(str(fpath))
            stored = StoredToolSource(
                hash=compute_hash(expanded),
                tool_source_class=type(ts).__name__,
                raw_source=expanded,
                tool_id=tool_id,
                tool_version=ts.parse_version(),
                tool_dir=str(fpath.parent),
                source_path=str(fpath),
                stored_at=datetime.datetime.now(datetime.UTC),
                metadata={"file_hash": file_hash},
            )
            store.store(stored)
            owner, repo = owner_repo_from_shed_yml(fpath.parent)
            disc = DiscoveredTool(
                path=str(fpath), tool_conf=None, tool_path=str(fpath), guid=None,
                is_shed_tool=False, tool_shed=None, repository_name=repo,
                repository_owner=owner, installed_changeset_revision=None,
                data_manager_id=None, hidden=False, labels=[], section_id=sid,
                section_name=sname, in_panel=True,
            )
            entry = build_index_entry_from_source(disc, stored, ts, biotools_source)
            if entry is None:
                n_errors += 1
                continue
            index.add_entry(entry)
            seen.add(tool_id)
            full = _entry_to_full(entry, fpath.parent, ts)
            slug = _check_slug_collision(slug_keys, owner, repo, tool_id)
            full["slug"] = slug
            full["source_dependency_hash"] = _tool_dependency_hash(fpath)
            write_yaml(tools_out / owner / repo / f"{slug}.yaml", full)
            index_rows.append(_slim(full))
            n += 1
        except Exception as e:
            n_errors += 1
            print(f"  ERROR parsing {fpath}: {e}", flush=True)
            continue
        progress(i, total, every=200, label="entries")

    store.store_index(index)
    store.close()
    index_rows.sort(key=lambda r: (r["owner"], r["repo"], r["id"]))
    write_yaml(DATA_DIR / "tools_index.yaml", index_rows)
    write_yaml(
        TOOLS_METADATA_PATH,
        {
            "tool_tree_signature": signature,
            "tool_count": n,
            "error_count": n_errors,
            "extracted_at": datetime.datetime.now(datetime.UTC).isoformat(),
        },
    )
    step_done("Building tool source store", f"{n} tools, {n_errors} errors")
    return n, n_errors


def main() -> None:
    step("Extracting tool metadata")
    n, errs = extract_tools()
    step_done("Extracting tool metadata", f"{n} tools, {errs} errors")
    print(f"  → {DATA_DIR / 'tools_index.yaml'}", flush=True)


if __name__ == "__main__":
    main()
