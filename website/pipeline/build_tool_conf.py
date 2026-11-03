"""Build a Galaxy toolbox configuration from the shared tracked-tool inventory."""

from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

from .common import REPO_ROOT, WEBSITE_DIR, step, step_done
from .inventory import published_tool_files, tool_section


def build_tool_conf(out_path: Path | None = None) -> Path:
    step("Building tool_conf.xml")
    out_path = out_path or WEBSITE_DIR / "tool_conf.xml"
    toolbox = ET.Element("toolbox", tool_path=str(REPO_ROOT))
    sections: dict[str, ET.Element] = {}
    for path in published_tool_files(REPO_ROOT):
        sid, name = tool_section(path, REPO_ROOT)
        if sid not in sections:
            sections[sid] = ET.SubElement(toolbox, "section", id=sid, name=name)
        ET.SubElement(sections[sid], "tool", file=str(path), dir=str(path.parent.relative_to(REPO_ROOT)))
    ET.indent(toolbox, space="  ")
    ET.ElementTree(toolbox).write(out_path, encoding="utf-8", xml_declaration=True)
    step_done("Building tool_conf.xml", f"{len(sections)} sections, {len(toolbox.findall('.//tool'))} tools")
    return out_path


def main() -> None:
    print(f"Wrote tool_conf: {build_tool_conf()}")


if __name__ == "__main__":
    main()
