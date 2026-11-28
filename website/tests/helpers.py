from __future__ import annotations

import os
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest import mock

# Every GIT_* variable is dropped: some redirect git to another repository, object store or
# config (hooks and `git rebase --exec` export them), others override the identity or dates
# the tests set up themselves (GIT_AUTHOR_*, GIT_COMMITTER_*).
_GIT_ENV_PREFIX = "GIT_"


def clean_git_env(env: dict[str, str] | None = None) -> dict[str, str]:
    """Return a copy of ``env`` (default ``os.environ``) without any GIT_* variable."""
    source = os.environ if env is None else env
    return {key: value for key, value in source.items() if not key.startswith(_GIT_ENV_PREFIX)}


class GitTestCase(unittest.TestCase):
    def setUp(self) -> None:
        # Also protects pipeline code under test that spawns git itself.
        patcher = mock.patch.dict(os.environ, clean_git_env(), clear=True)
        patcher.start()
        self.addCleanup(patcher.stop)
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.repo = Path(directory.name)
        self.git("init", "--quiet")
        # Local to the temporary repository; tests switch authors the same way.
        self.git("config", "user.name", "Website test")
        self.git("config", "user.email", "test@example.invalid")

    def git(self, *args: str) -> str:
        return subprocess.check_output(
            ["git", "-C", str(self.repo), *args], text=True, env=clean_git_env()
        ).strip()

    def write(self, name: str, content: str) -> Path:
        path = self.repo / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def commit(self, message: str) -> str:
        self.git("add", ".")
        self.git("commit", "--quiet", "-m", message)
        return self.git("rev-parse", "HEAD")
