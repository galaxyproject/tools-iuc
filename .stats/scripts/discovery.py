from pathlib import Path
import re
import subprocess
import yaml

from .data_model import ToolInfo

TOOL_ID_RE = re.compile(r'<tool\s+id=["\']([^"\']+)["\']', re.IGNORECASE)
TOOL_NAME_RE = re.compile(
    r'<tool\s[^>]*?name=["\']([^"\']+)["\']', re.IGNORECASE
)
SKIP_DIRS = {'test-data', '.vscode'}


def _tracked_xml_files(repo_path: str, tools_path: str) -> list[Path]:
    """Return tracked XML files under tools_path, excluding untracked artifacts."""
    result = subprocess.run(
        ['git', 'ls-files', '--', tools_path],
        cwd=repo_path,
        capture_output=True,
        text=True,
    )
    paths: list[Path] = []
    if result.returncode != 0:
        return paths
    repo = Path(repo_path)
    for line in result.stdout.splitlines():
        if not line.endswith('.xml'):
            continue
        paths.append(repo / line)
    return paths


def discover_tools(repo_path: str, tools_path: str = 'tools') -> list[ToolInfo]:
    tools_root = Path(repo_path) / tools_path
    xml_files = _tracked_xml_files(repo_path, tools_path)

    files_by_dir: dict[str, list[Path]] = {}
    for xml_path in xml_files:
        try:
            rel = xml_path.relative_to(tools_root)
        except ValueError:
            continue
        dir_name = str(rel.parent)
        if dir_name == '.' or any(
            part in SKIP_DIRS or part.startswith('.')
            for part in rel.parent.parts
        ):
            continue
        files_by_dir.setdefault(dir_name, []).append(xml_path)

    tools: list[ToolInfo] = []
    seen_ids: set[str] = set()

    for dir_name in sorted(files_by_dir):
        for xml_path in sorted(files_by_dir[dir_name]):
            text = xml_path.read_text(encoding='utf-8', errors='replace')
            m = TOOL_ID_RE.search(text)
            if not m:
                continue
            tid = m.group(1)
            if tid in seen_ids:
                continue
            seen_ids.add(tid)
            tname = _parse_tool_name(text)

            categories: list[str] = []
            shed_yml = tools_root / dir_name / '.shed.yml'
            if shed_yml.exists():
                try:
                    data = yaml.safe_load(
                        shed_yml.read_text(encoding='utf-8')
                    )
                    if data and 'categories' in data:
                        cats = data['categories']
                        if isinstance(cats, list):
                            categories = cats
                        else:
                            categories = [cats]
                except Exception:
                    pass

            tools.append(ToolInfo(
                tool_id=tid,
                tool_name=tname or tid,
                directory=dir_name,
                xml_files=[xml_path.name],
                categories=categories,
            ))

    return tools


def _parse_tool_name(xml_text: str) -> str | None:
    m = TOOL_NAME_RE.search(xml_text)
    return m.group(1) if m else None
