from __future__ import annotations

import subprocess
import tempfile
import unittest
from pathlib import Path


class GitTestCase(unittest.TestCase):
    def setUp(self) -> None:
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.repo = Path(directory.name)
        self.git("init", "--quiet")
        self.git("config", "user.name", "Website test")
        self.git("config", "user.email", "test@example.invalid")

    def git(self, *args: str) -> str:
        return subprocess.check_output(["git", "-C", str(self.repo), *args], text=True).strip()

    def write(self, name: str, content: str) -> Path:
        path = self.repo / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def commit(self, message: str) -> str:
        self.git("add", ".")
        self.git("commit", "--quiet", "-m", message)
        return self.git("rev-parse", "HEAD")
