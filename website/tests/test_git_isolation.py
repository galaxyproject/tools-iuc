from __future__ import annotations

import os
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from tests.helpers import GitTestCase, clean_git_env


class _Probe(GitTestCase):
    def runTest(self) -> None:  # pragma: no cover - driven by hand below
        pass


def _snapshot(git_dir: Path) -> tuple[str, str, str]:
    env = clean_git_env()

    def run(*args: str) -> str:
        return subprocess.check_output(["git", f"--git-dir={git_dir}", *args], text=True, env=env)

    return (
        (git_dir / "config").read_text(encoding="utf-8"),
        run("for-each-ref"),
        (git_dir / "HEAD").read_text(encoding="utf-8"),
    )


class GitIsolationTests(unittest.TestCase):
    def test_clean_git_env_drops_repository_redirecting_variables(self) -> None:
        env = {
            "PATH": "/usr/bin",
            "GIT_DIR": "/x",
            "GIT_WORK_TREE": "/x",
            "GIT_INDEX_FILE": "/x",
            "GIT_COMMON_DIR": "/x",
            "GIT_OBJECT_DIRECTORY": "/x",
            "GIT_ALTERNATE_OBJECT_DIRECTORIES": "/x",
            "GIT_CONFIG_COUNT": "1",
            "GIT_CONFIG_KEY_0": "core.bare",
            "GIT_CONFIG_GLOBAL": "/x",
        }

        self.assertEqual(clean_git_env(env), {"PATH": "/usr/bin"})

    def test_clean_git_env_drops_every_git_variable(self) -> None:
        env = {"HOME": "/h", "GIT_AUTHOR_NAME": "x", "GIT_COMMITTER_EMAIL": "x", "GIT_AUTHOR_DATE": "x"}

        self.assertEqual(clean_git_env(env), {"HOME": "/h"})

    def test_identity_variables_do_not_override_the_test_repository_identity(self) -> None:
        identity = {
            "GIT_AUTHOR_NAME": "Outer Author",
            "GIT_AUTHOR_EMAIL": "outer-author@example.invalid",
            "GIT_AUTHOR_DATE": "2001-01-01T00:00:00Z",
            "GIT_COMMITTER_NAME": "Outer Committer",
            "GIT_COMMITTER_EMAIL": "outer-committer@example.invalid",
            "GIT_COMMITTER_DATE": "2001-01-01T00:00:00Z",
        }
        with mock.patch.dict(os.environ, identity):
            probe = _Probe()
            probe.setUp()
            try:
                probe.write("a.txt", "a")
                probe.commit("Probe commit")
                author = probe.git("log", "-1", "--format=%an <%ae> %ad", "--date=short")
                committer = probe.git("log", "-1", "--format=%cn <%ce>")
            finally:
                probe.doCleanups()

        self.assertEqual(author.split(" <")[0], "Website test")
        self.assertNotIn("2001-01-01", author)
        self.assertEqual(committer, "Website test <test@example.invalid>")

    def test_helper_does_not_touch_a_repository_named_by_git_dir(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            scratch = Path(directory)
            subprocess.check_call(["git", "-C", str(scratch), "init", "--quiet"], env=clean_git_env())
            git_dir = scratch / ".git"
            before = _snapshot(git_dir)

            polluted = {"GIT_DIR": str(git_dir), "GIT_WORK_TREE": str(scratch), "GIT_INDEX_FILE": str(git_dir / "index")}
            with mock.patch.dict(os.environ, polluted):
                probe = _Probe()
                probe.setUp()
                try:
                    probe.write("tools/example/tool.xml", "<tool/>")
                    commit = probe.commit("Probe commit")
                    self.assertTrue(commit)
                    self.assertNotEqual(Path(probe.repo).resolve(), scratch.resolve())
                    self.assertTrue((probe.repo / ".git").is_dir())
                    self.assertNotIn("GIT_DIR", os.environ)
                finally:
                    probe.doCleanups()
                self.assertEqual(os.environ["GIT_DIR"], str(git_dir))

            self.assertEqual(_snapshot(git_dir), before)
            self.assertFalse((scratch / "tools").exists())


if __name__ == "__main__":
    unittest.main()
