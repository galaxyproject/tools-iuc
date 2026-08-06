"""Synthesize a Galaxy toolbox-conf XML listing every IUC tool XML.

This is the input format Galaxy's tool-source populator expects. We group tools
into <section>s by their parent directory so the derived panel_section metadata
is meaningful, and we also surface data_managers and suites.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

from .common import (
    DATA_MANAGERS_DIR,
    REPO_ROOT,
    SUITES_DIR,
    TOOLS_DIR,
    WEBSITE_DIR,
    git_tracked_files,
    step,
    step_done,
)


def _section_name(dir_name: str) -> tuple[str, str]:
    """Map a tool directory name to a (section_id, section_name). Heuristic: use the dir name."""
    return dir_name, dir_name.replace("_", " ").title()


def _collect_xmls(base: Path) -> list[Path]:
    """Return tracked *.xml files under base, excluding test-data and macros."""
    if not base.exists():
        return []
    files = git_tracked_files(REPO_ROOT, str(base.relative_to(REPO_ROOT)))
    out: list[Path] = []
    for f in files:
        if not f.suffix == ".xml":
            continue
        if "test-data" in f.parts:
            continue
        name = f.name
        # skip macros/test-only files; keep the main tool xml + any .xml that has <tool id=...>
        if name in ("macros.xml", "tool_conf.xml", "repository_dependencies.xml", "tool_data_table_conf.xml"):
            continue
        out.append(f)
    return out


def _has_tool_id(path: Path) -> bool:
    try:
        tree = ET.parse(path)
        root = tree.getroot()
        return root.tag == "tool" and root.get("id") is not None
    except Exception:
        return False


def build_tool_conf(out_path: Path | None = None) -> Path:
    """Build a tool_conf.xml and return the path written."""
    step("Building tool_conf.xml")
    if out_path is None:
        out_path = WEBSITE_DIR / "tool_conf.xml"

    toolbox = ET.Element("toolbox")
    toolbox.set("tool_path", str(REPO_ROOT))

    # tools/ grouped by directory -> section
    for child in sorted(TOOLS_DIR.iterdir()) if TOOLS_DIR.exists() else []:
        if not child.is_dir():
            continue
        xmls = [f for f in _collect_xmls(child) if _has_tool_id(f)]
        if not xmls:
            continue
        sid, sname = _section_name(child.name)
        section = ET.SubElement(toolbox, "section", id=sid, name=sname)
        for xml in sorted(xmls):
            t = ET.SubElement(section, "tool", file=str(xml))
            t.set("dir", str(xml.parent.relative_to(REPO_ROOT)))

    # data_managers/ in a dedicated section
    dm_xmls = [f for f in _collect_xmls(DATA_MANAGERS_DIR) if _has_tool_id(f)]
    if dm_xmls:
        section = ET.SubElement(toolbox, "section", id="data_managers", name="Data Managers")
        for xml in sorted(dm_xmls):
            t = ET.SubElement(section, "tool", file=str(xml))
            t.set("dir", str(xml.parent.relative_to(REPO_ROOT)))

    # suites/
    suite_xmls = [f for f in _collect_xmls(SUITES_DIR) if _has_tool_id(f)]
    if suite_xmls:
        section = ET.SubElement(toolbox, "section", id="suites", name="Suites")
        for xml in sorted(suite_xmls):
            t = ET.SubElement(section, "tool", file=str(xml))
            t.set("dir", str(xml.parent.relative_to(REPO_ROOT)))

    ET.indent(toolbox, space="  ")
    xml_bytes = ET.tostring(toolbox, encoding="unicode", xml_declaration=False)
    out_path.write_text('<?xml version="1.0"?>\n' + xml_bytes + "\n", encoding="utf-8")
    n_sections = len(toolbox.findall("section"))
    n_tools = len(toolbox.findall(".//tool"))
    step_done("Building tool_conf.xml", f"{n_sections} sections, {n_tools} tools")
    return out_path


def main() -> None:
    path = build_tool_conf()
    print(f"Wrote tool_conf: {path}")


if __name__ == "__main__":
    main()
