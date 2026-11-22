from functools import lru_cache
import re
from dataclasses import dataclass, field
from typing import Optional

from git import Repo
from git.exc import GitCommandError


_repo_cache: dict[str, Repo] = {}


def _repo(repo_path: str) -> Repo:
    if repo_path not in _repo_cache:
        _repo_cache[repo_path] = Repo(repo_path)
    return _repo_cache[repo_path]

TOKEN_PATTERN = re.compile(r'@(\w+)@')
TOOL_ID_RE = re.compile(r'<tool\s+id=["\']([^"\']+)["\']', re.IGNORECASE)
TOOL_VERSION_RE = re.compile(
    r'<tool\s[^>]*?version=["\']([^"\']+)["\']', re.IGNORECASE
)
IMPORT_RE = re.compile(r'<import>([^<]+)</import>')
TOKEN_DEF_RE = re.compile(
    r'<token\s+name=[\'"]@(\w+)@[\'"]>(.*?)</token>', re.DOTALL
)


@dataclass
class ResolvedVersion:
    full: str
    template: str
    tokens: dict[str, str] = field(default_factory=dict)
    upstream_token_val: str | None = None
    wrapper_token_val: str | None = None


_KNOWN_UPSTREAM_TOKENS = [
    'TOOL_VERSION', 'VERSION',
]
_KNOWN_WRAPPER_TOKENS = [
    'VERSION_SUFFIX', 'WRAPPER_VERSION',
]


def _parse_macros(content: str) -> dict[str, str]:
    tokens: dict[str, str] = {}
    for m in TOKEN_DEF_RE.finditer(content):
        name = m.group(1)
        value = m.group(2).strip()
        tokens[name] = value
    for _ in range(10):
        changed = False
        for name in list(tokens):
            val = tokens[name]
            unresolved = TOKEN_PATTERN.findall(val)
            if unresolved:
                for tok in unresolved:
                    if tok in tokens and tok != name:
                        val = val.replace(f'@{tok}@', tokens[tok])
                        changed = True
                tokens[name] = val
        if not changed:
            break
    return tokens


def _read_at_commit(repo_path: str, commit_hash: str, path: str) -> str | None:
    try:
        return _repo(repo_path).git.show(f'{commit_hash}:{path}')
    except GitCommandError:
        return None


@lru_cache(maxsize=4096)
def _get_dir_xml_files(
    repo_path: str,
    commit_hash: str,
    tool_dir: str,
) -> tuple[tuple[str, tuple[str, ...]], ...]:
    """Return a mapping of tool_id -> XML filenames in tools/{tool_dir}.

    Results are cached per (commit_hash, tool_dir) because many tools share
    the same directory. The cache is bounded to 4096 entries.
    """
    prefix = f'tools/{tool_dir}'
    mapping: dict[str, list[str]] = {}
    try:
        output = _repo(repo_path).git.ls_tree(
            '-r', '--name-only', commit_hash, prefix
        )
    except GitCommandError:
        return ()

    full_paths: list[str] = []
    for line in output.splitlines():
        if line.endswith('.xml'):
            rel = line[len(prefix) + 1:]
            if '/' not in rel:
                full_paths.append(line)

    for full_path in full_paths:
        content = _read_at_commit(repo_path, commit_hash, full_path)
        if content is None:
            continue
        m = TOOL_ID_RE.search(content)
        if m:
            tid = m.group(1)
            filename = full_path.split('/')[-1]
            mapping.setdefault(tid, []).append(filename)

    return tuple((tid, tuple(files)) for tid, files in mapping.items())


def find_xml_files_for_tool_id(
    repo_path: str,
    commit_hash: str,
    tool_dir: str,
    tool_id: str,
) -> list[str]:
    """Return XML filenames in tools/{tool_dir} at commit that define tool_id."""
    mapping = _get_dir_xml_files(repo_path, commit_hash, tool_dir)
    for tid, files in mapping:
        if tid == tool_id:
            return list(files)
    return []


def _collect_macros(
    repo_path: str, commit_hash: str, tool_dir: str, xml_files: list[str]
) -> dict[str, str]:
    macros: dict[str, str] = {}
    visited: set[str] = set()

    def load(path: str):
        if path in visited:
            return
        visited.add(path)
        content = _read_at_commit(repo_path, commit_hash, path)
        if content is not None:
            macros.update(_parse_macros(content))

    load(f'tools/{tool_dir}/macros.xml')

    for xml_name in sorted(xml_files):
        xml_content = _read_at_commit(
            repo_path, commit_hash, f'tools/{tool_dir}/{xml_name}'
        )
        if xml_content is None:
            continue

        macros.update(_parse_macros(xml_content))

        for imp in IMPORT_RE.finditer(xml_content):
            imported = imp.group(1)
            if imported == 'macros.xml':
                continue
            if imported.endswith('.xml'):
                load(f'tools/{tool_dir}/{imported}')

    shared = _read_at_commit(
        repo_path, commit_hash, 'macros/read_group_macros.xml'
    )
    if shared is not None:
        macros.update(_parse_macros(shared))

    return macros


def resolve_version(
    repo_path: str,
    commit_hash: str,
    tool_dir: str,
    xml_files: list[str],
) -> ResolvedVersion | None:
    macros = _collect_macros(repo_path, commit_hash, tool_dir, xml_files)

    for xml_name in sorted(xml_files):
        xml_content = _read_at_commit(
            repo_path, commit_hash, f'tools/{tool_dir}/{xml_name}'
        )
        if xml_content is None:
            continue

        m = TOOL_VERSION_RE.search(xml_content)
        if not m:
            continue

        template = m.group(1)

        resolved = template
        for _ in range(10):
            unresolved = TOKEN_PATTERN.findall(resolved)
            if not unresolved:
                break
            for tok in unresolved:
                if tok in macros:
                    resolved = resolved.replace(
                        f'@{tok}@', macros[tok]
                    )

        if '@' in resolved:
            continue

        upstream_val = None
        wrapper_val = None
        for tokname in TOKEN_PATTERN.findall(template):
            tokval = macros.get(tokname)
            if tokval is not None:
                if tokname in _KNOWN_UPSTREAM_TOKENS:
                    upstream_val = tokval
                if tokname in _KNOWN_WRAPPER_TOKENS:
                    wrapper_val = tokval

        return ResolvedVersion(
            full=resolved,
            template=template,
            tokens=dict(macros),
            upstream_token_val=upstream_val,
            wrapper_token_val=wrapper_val,
        )

    return None


def classify_change(
    old_v: ResolvedVersion,
    new_v: ResolvedVersion,
) -> str:
    if old_v.full == new_v.full:
        return 'no-change'

    old_up = old_v.upstream_token_val
    new_up = new_v.upstream_token_val
    old_wr = old_v.wrapper_token_val
    new_wr = new_v.wrapper_token_val

    if old_up is not None and new_up is not None:
        if old_up != new_up:
            return 'upstream'
        return 'wrapper'

    if old_wr is not None and new_wr is not None:
        if old_wr != new_wr:
            return 'wrapper'

    if old_v.template != new_v.template:
        return 'wrapper'

    return 'version'


def extract_tool_version_direct(xml_text: str) -> str | None:
    m = TOOL_VERSION_RE.search(xml_text)
    return m.group(1) if m else None
