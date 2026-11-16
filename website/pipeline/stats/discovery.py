"""Discover history targets from the same tool inventory as metadata extraction."""

from pathlib import Path

from galaxy.tool_util.parser import get_tool_source

from ..common import read_yaml
from ..inventory import TOOL_ROOTS, published_tool_files
from .data_model import ToolInfo


def discover_tools(repo_path: str, tools_path: str | None = None) -> list[ToolInfo]:
    repo = Path(repo_path).resolve()
    tools = []
    for path in published_tool_files(repo, (tools_path,) if tools_path else TOOL_ROOTS):
        source = get_tool_source(config_file=str(path))
        shed = path.parent / ".shed.yml"
        config = read_yaml(shed) if shed.exists() else {}
        tools.append(ToolInfo(
            tool_id=source.parse_id(),
            tool_name=source.parse_name(),
            directory=path.parent.relative_to(repo).as_posix(),
            xml_files=[path.name],
            categories=(config or {}).get("categories", []),
        ))
    return tools
