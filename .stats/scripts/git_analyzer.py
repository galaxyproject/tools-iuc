from collections import defaultdict
from datetime import datetime, timezone
import subprocess
import sys
from typing import Optional

from git import Repo
from git.exc import GitCommandError

from .data_model import ToolInfo, ToolEvent
from .version_resolver import (
    find_xml_files_for_tool_id,
    resolve_version,
    classify_change,
)


StatusLine = tuple[str, str, Optional[str]]
CommitStatus = tuple[str, int, str, list[StatusLine]]


def _get_main_branch(repo_path: str) -> str:
    repo = Repo(repo_path)
    for candidate in ['origin/main', 'origin/master', 'main', 'master']:
        try:
            repo.git.rev_parse('--verify', candidate)
            return candidate
        except GitCommandError:
            continue
    raise RuntimeError('No main/master branch found')


def _dir_from_path(path: str, tools_path: str = 'tools') -> str | None:
    """Return the tool directory for a path like tools/dir/file.xml."""
    prefix = tools_path + '/'
    if not path.startswith(prefix):
        return None
    rel = path[len(prefix):]
    if '/' not in rel:
        return None
    return rel.rsplit('/', 1)[0]


def get_commit_data(
    repo_path: str, branch: str, tools_path: str = 'tools'
) -> list[CommitStatus]:
    """Return commit hash, timestamp, message, and name-status lines."""
    repo = Repo(repo_path)
    log_format = '%H||%ct||%s'
    args = [
        'log',
        '--first-parent',
        '--topo-order',
        branch,
        f'--format={log_format}',
        '--name-status',
        '--',
        tools_path,
    ]
    _, raw, _ = repo.git.execute(
        ['git'] + args,
        with_extended_output=True,
    )

    results: list[CommitStatus] = []
    current_hash: str | None = None
    current_ts: int = 0
    current_msg: str = ''
    current_statuses: list[StatusLine] = []

    for line in raw.splitlines():
        line = line.strip()
        if not line:
            continue
        if '||' in line:
            if current_hash is not None:
                results.append(
                    (current_hash, current_ts, current_msg, current_statuses)
                )
                current_statuses = []
            parts = line.split('||', 2)
            current_hash = parts[0]
            current_ts = int(parts[1])
            current_msg = parts[2] if len(parts) > 2 else ''
            continue

        # name-status line: status [old_path [new_path]]
        fields = line.split('\t')
        if len(fields) < 2:
            continue
        status = fields[0]
        old_path = fields[1]
        new_path = fields[2] if len(fields) > 2 else None
        current_statuses.append((status, old_path, new_path))

    if current_hash is not None:
        results.append(
            (current_hash, current_ts, current_msg, current_statuses)
        )

    return results


def _extract_tool_id_at_commit(
    repo_path: str, commit_hash: str, path: str
) -> str | None:
    """Read an XML file at a commit and extract its <tool id="...">."""
    try:
        content = subprocess.check_output(
            ['git', 'show', f'{commit_hash}:{path}'],
            cwd=repo_path,
            text=True,
            stderr=subprocess.DEVNULL,
            timeout=30,
        )
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return None
    from .version_resolver import TOOL_ID_RE
    m = TOOL_ID_RE.search(content)
    return m.group(1) if m else None


def _build_location_history(
    repo_path: str,
    commit_data: list[CommitStatus],
    tracked_tool_ids: set[str],
    tools_path: str = 'tools',
    verbose: bool = False,
) -> dict[str, dict[str, tuple[str, str]]]:
    """Map commit hash -> {tool_id -> (directory, xml_file)}.

    Processes commits chronologically, tracking file renames and directory
    moves. Only locations for tool ids discovered in the current tree are
    retained.
    """
    prefix = tools_path + '/'
    location_history: dict[str, dict[str, tuple[str, str]]] = {}
    current_map: dict[str, tuple[str, str]] = {}

    for seq, (chash, cts, cmsg, statuses) in enumerate(commit_data):
        for status, old_path, new_path in statuses:
            paths_to_add: list[str] = []
            paths_to_remove: list[str] = []

            if status == 'A':
                if new_path:
                    paths_to_add.append(new_path)
            elif status == 'D':
                if old_path:
                    paths_to_remove.append(old_path)
            elif status.startswith('R'):
                if old_path:
                    paths_to_remove.append(old_path)
                if new_path:
                    paths_to_add.append(new_path)
            elif status == 'M':
                if old_path:
                    paths_to_add.append(old_path)
            else:
                continue

            for path in paths_to_remove:
                if not path.startswith(prefix) or not path.endswith('.xml'):
                    continue
                rel = path[len(prefix):]
                for tid, (d, f) in list(current_map.items()):
                    if f'{d}/{f}' == rel:
                        del current_map[tid]
                        break

            for path in paths_to_add:
                if not path.startswith(prefix) or not path.endswith('.xml'):
                    continue
                tool_id = _extract_tool_id_at_commit(repo_path, chash, path)
                if tool_id is None or tool_id not in tracked_tool_ids:
                    continue
                rel = path[len(prefix):]
                if '/' not in rel:
                    continue
                dir_name, filename = rel.rsplit('/', 1)
                current_map[tool_id] = (dir_name, filename)

        location_history[chash] = current_map.copy()

    if verbose:
        total = sum(len(v) for v in location_history.values())
        avg = total / len(location_history) if location_history else 0
        print(
            f'Built location history for {len(location_history)} commits '
            f'({avg:.1f} tools/commit)',
            file=sys.stderr,
        )

    return location_history


def _process_tool(
    tool: ToolInfo,
    commits: list[tuple[str, int, str, int]],
    location_history: dict[str, dict[str, tuple[str, str]]],
    repo_path: str,
) -> list[ToolEvent]:
    """Reconstruct events for a single tool."""
    last_version = None
    tool_events: list[ToolEvent] = []

    for chash, cts, cmsg, _ in commits:
        location = location_history.get(chash, {}).get(tool.tool_id)

        if location is None:
            # Tool did not exist at this commit (e.g. before addition or
            # during a directory move window not captured for this id).
            continue

        directory, xml_name = location
        xml_files = [xml_name]
        resolved = resolve_version(repo_path, chash, directory, xml_files)

        # If the recorded XML filename did not exist at this commit (rare,
        # can happen with concurrent macro-only commits near a rename),
        # fall back to searching the directory for the tool id.
        if resolved is None:
            renamed = find_xml_files_for_tool_id(
                repo_path, chash, directory, tool.tool_id
            )
            if renamed:
                resolved = resolve_version(
                    repo_path, chash, directory, renamed
                )

        if resolved is None and last_version is None:
            continue

        if last_version is None and resolved is not None:
            tool_events.append(ToolEvent(
                tool_id=tool.tool_id,
                tool_name=tool.tool_name,
                directory=directory,
                event_type='new',
                version=resolved.full,
                commit_hash=chash,
                commit_date=datetime.fromtimestamp(cts, tz=timezone.utc),
                change_type='initial',
                commit_message=cmsg,
            ))
            last_version = resolved
            continue

        if resolved is not None:
            ctype = classify_change(last_version, resolved)
            if ctype != 'no-change':
                tool_events.append(ToolEvent(
                    tool_id=tool.tool_id,
                    tool_name=tool.tool_name,
                    directory=directory,
                    event_type='update',
                    version=resolved.full,
                    old_version=last_version.full,
                    commit_hash=chash,
                    commit_date=datetime.fromtimestamp(
                        cts, tz=timezone.utc
                    ),
                    change_type=ctype,
                    commit_message=cmsg,
                ))
                last_version = resolved

    return tool_events


def analyze_tools(
    repo_path: str,
    tools: list[ToolInfo],
    branch: str | None = None,
    tools_path: str = 'tools',
    verbose: bool = False,
    workers: int = 1,
) -> list[ToolEvent]:
    if branch is None:
        branch = _get_main_branch(repo_path)

    if verbose:
        print(f'Using branch: {branch}', file=sys.stderr)

    commit_data = get_commit_data(repo_path, branch, tools_path)
    commit_data.reverse()

    if verbose:
        print(
            f'Found {len(commit_data)} commits touching {tools_path}/ on '
            f'{branch}',
            file=sys.stderr,
        )

    tracked_tool_ids = {t.tool_id for t in tools}
    location_history = _build_location_history(
        repo_path, commit_data, tracked_tool_ids, tools_path, verbose
    )

    tool_commits: dict[str, list[tuple[str, int, str, int]]] = defaultdict(list)

    for seq, (chash, cts, cmsg, statuses) in enumerate(commit_data):
        touched_dirs: set[str] = set()
        for status, old_path, new_path in statuses:
            for path in (old_path, new_path):
                if path is None:
                    continue
                d = _dir_from_path(path, tools_path)
                if d is not None:
                    touched_dirs.add(d)
        for tool in tools:
            if tool.directory in touched_dirs:
                tool_commits[tool.tool_id].append((chash, cts, cmsg, seq))

    events: list[ToolEvent] = []
    tool_map = {t.tool_id: t for t in tools}

    tool_ids = sorted(tool_commits.keys())

    if workers > 1:
        import concurrent.futures
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
            futures = {
                ex.submit(
                    _process_tool,
                    tool_map[tid],
                    tool_commits[tid],
                    location_history,
                    repo_path,
                ): tid
                for tid in tool_ids
            }
            for future in concurrent.futures.as_completed(futures):
                tid = futures[future]
                tool_events = future.result()
                events.extend(tool_events)
                if verbose:
                    updates = len(
                        [e for e in tool_events if e.event_type == 'update']
                    )
                    print(
                        f'  {tid:30s}: '
                        f'{len(tool_commits[tid]):3d} commits, '
                        f'{updates} updates',
                        file=sys.stderr,
                    )
    else:
        for tid in tool_ids:
            tool = tool_map[tid]
            tool_events = _process_tool(
                tool, tool_commits[tid], location_history, repo_path
            )
            events.extend(tool_events)
            if verbose:
                updates = len(
                    [e for e in tool_events if e.event_type == 'update']
                )
                print(
                    f'  {tool.tool_id:30s}: '
                    f'{len(tool_commits[tid]):3d} commits, '
                    f'{updates} updates',
                    file=sys.stderr,
                )

    events.sort(key=lambda e: e.commit_date or datetime.min)
    return events
