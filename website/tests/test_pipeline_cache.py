from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest import mock

from pipeline import extract_stats, extract_tools
from pipeline.common import md5_files, write_yaml


class CacheTests(unittest.TestCase):
    def test_md5_files_changes_when_dependency_changes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first = root / "tool.xml"
            second = root / "macros.xml"
            first.write_text("<tool id='x'/>", encoding="utf-8")
            second.write_text("<macros/>", encoding="utf-8")

            before = md5_files([first, second])
            second.write_text("<macros><token name='@VERSION@'>1</token></macros>", encoding="utf-8")

            self.assertNotEqual(before, md5_files([first, second]))

    def test_tool_yaml_cache_requires_matching_dependency_hash(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tool_yaml = Path(tmp) / "tool.yaml"
            write_yaml(tool_yaml, {"id": "example", "source_dependency_hash": "abc"})

            self.assertTrue(extract_tools._tool_yaml_cache_fresh(tool_yaml, "abc"))
            self.assertFalse(extract_tools._tool_yaml_cache_fresh(tool_yaml, "def"))

    def test_stats_cache_reuses_existing_events_for_same_commit(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            data_dir = Path(tmp)
            stats_dir = data_dir / "stats"
            stats_dir.mkdir()
            events = stats_dir / "tool_events.json"
            events.write_text("[]", encoding="utf-8")
            write_yaml(stats_dir / "metadata.yaml", {"git_commit": "abc123"})

            with mock.patch.object(extract_stats, "DATA_DIR", data_dir):
                with mock.patch.object(extract_stats, "git_head_commit", return_value="abc123"):
                    self.assertEqual(extract_stats._cached_events_path(), events)

    def test_stats_cache_misses_for_different_commit(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            data_dir = Path(tmp)
            stats_dir = data_dir / "stats"
            stats_dir.mkdir()
            events = stats_dir / "tool_events.json"
            events.write_text("[]", encoding="utf-8")
            write_yaml(stats_dir / "metadata.yaml", {"git_commit": "old"})

            with mock.patch.object(extract_stats, "DATA_DIR", data_dir):
                with mock.patch.object(extract_stats, "git_head_commit", return_value="new"):
                    self.assertIsNone(extract_stats._cached_events_path())


if __name__ == "__main__":
    unittest.main()
