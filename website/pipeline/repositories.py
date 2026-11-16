"""Adapt Planemo's publishing rules to a source-file → ToolShed repository map.

Planemo's expansion is used without realizing upload directories (which would
copy every suite's test data). Its internal helpers are confined to this
adapter and covered by publishing-rule integration tests.
"""

from __future__ import annotations

from functools import cache
from pathlib import Path

from planemo.context import PlanemoContext
from planemo.shed import RawRepositoryDirectory, _glob, _shed_config_excludes, shed_repo_config, shed_repo_type

from .common import IUC_OWNER


@cache
def _repository_files(directory: Path) -> dict[Path, set[tuple[str, str]]]:
    config = shed_repo_config(PlanemoContext(), str(directory))
    raw = RawRepositoryDirectory(str(directory), config, multiple=False)
    published: dict[Path, set[tuple[str, str]]] = {}
    for name in config["repositories"]:
        effective = raw._realize_config(name)
        if shed_repo_type(effective, name) != "unrestricted":
            continue
        excluded = {
            Path(path).absolute()
            for pattern in _shed_config_excludes(effective)
            for path in _glob(str(directory), pattern)
        }
        repository = (effective.get("owner") or IUC_OWNER, name)
        for file in raw._realized_files(name).files:
            source = Path(file.absolute_src)
            if source.suffix != ".xml" or source in excluded or raw._implicit_ignores(file.dest):
                continue
            published.setdefault(source, set()).add(repository)
    return published


def repository_for_tool(path: Path) -> tuple[str, str] | None:
    path = path.absolute()
    for directory in path.parents:
        if (directory / ".shed.yml").exists():
            repositories = _repository_files(directory).get(path, set())
            if len(repositories) > 1:
                raise ValueError(f"{path} is published in multiple repositories: {sorted(repositories)}")
            return next(iter(repositories), None)
        if (directory / ".git").exists():
            break
    raise ValueError(f"No .shed.yml found for {path}")
