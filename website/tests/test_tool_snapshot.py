from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from pipeline import extract_people
from pipeline.common import read_yaml, write_yaml
from pipeline.tool_snapshot import load_tool_index, tool_snapshot


class SnapshotTests(unittest.TestCase):
    def setUp(self) -> None:
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.data = Path(directory.name)

    def test_publication_replaces_the_entire_tool_generation(self) -> None:
        with tool_snapshot(self.data) as stage:
            write_yaml(stage / "index.yaml", [{"id": "old"}])
            write_yaml(stage / "old.yaml", {"id": "old"})
        old = (self.data / "tools").resolve()
        with tool_snapshot(self.data) as stage:
            write_yaml(stage / "index.yaml", [{"id": "new"}])
            self.assertEqual(load_tool_index(self.data), [{"id": "old"}])
        self.assertEqual(load_tool_index(self.data), [{"id": "new"}])
        self.assertFalse((self.data / "tools" / "old.yaml").exists())
        self.assertTrue(old.exists())

    def test_failed_generation_preserves_previous_snapshot_and_discards_staging(self) -> None:
        with tool_snapshot(self.data) as stage:
            write_yaml(stage / "index.yaml", [{"id": "old"}])
        with self.assertRaisesRegex(RuntimeError, "parse failed"):
            with tool_snapshot(self.data) as failed:
                write_yaml(failed / "index.yaml", [])
                raise RuntimeError("parse failed")
        self.assertEqual(load_tool_index(self.data), [{"id": "old"}])
        self.assertFalse(failed.exists())

    def test_legacy_directory_is_replaced_without_preserving_obsolete_tools(self) -> None:
        write_yaml(self.data / "tools" / "deleted.yaml", {"id": "deleted"})
        with tool_snapshot(self.data) as stage:
            write_yaml(stage / "index.yaml", [])
        self.assertTrue((self.data / "tools").is_symlink())
        self.assertFalse((self.data / "tools" / "deleted.yaml").exists())

    def test_failed_pointer_switch_restores_legacy_directory(self) -> None:
        write_yaml(self.data / "tools" / "old.yaml", {"id": "old"})
        with patch("pipeline.tool_snapshot.os.replace", side_effect=OSError("cannot publish")):
            with self.assertRaisesRegex(OSError, "cannot publish"):
                with tool_snapshot(self.data) as stage:
                    write_yaml(stage / "index.yaml", [])
        self.assertEqual(read_yaml(self.data / "tools" / "old.yaml"), {"id": "old"})
        self.assertFalse((self.data / ".tools-next").exists())

    def test_people_use_manifest_and_emit_canonical_creator_ids(self) -> None:
        row = {"id": "example", "owner": "iuc", "repo": "example", "slug": "example"}
        creators = [
            {"class": "Person", "name": "Different Display Name", "identifier": "0000-0001-6431-3442"},
            {"class": "Person", "name": "Legacy Name", "url": "https://github.com/old-handle"},
            {"class": "Organization", "name": "Example Organisation"},
        ]
        with tool_snapshot(self.data) as stage:
            write_yaml(stage / "index.yaml", [row])
            write_yaml(stage / "iuc/example/example.yaml", {
                "id": "example", "source_path": "tools/example/tool.xml", "creators": creators,
            })
            write_yaml(stage / "stale.yaml", {
                "id": "deleted", "creators": [{"name": "Deleted Creator"}],
            })
        gtn = {"contributors": {"canonical": {"name": "Canonical Name", "orcid": "0000-0001-6431-3442"}}}
        with (
            patch.object(extract_people, "DATA_DIR", self.data),
            patch.object(extract_people, "fetch_gtn_people", return_value=gtn),
            patch.object(extract_people, "_git_authors", return_value=[]),
            patch.object(extract_people, "_load_aliases", return_value={"new-handle": {"github": ["old-handle"]}}),
        ):
            extract_people.extract_people()
        contributors = read_yaml(self.data / "contributors.yaml")
        self.assertEqual({person["id"] for person in contributors}, {"canonical", "new-handle"})
        self.assertTrue(all(person["tools"] == ["example"] for person in contributors))
        self.assertEqual(read_yaml(self.data / "creator_links.yaml"), {
            "tools/example/tool.xml": [
                {"kind": "contributor", "id": "canonical"},
                {"kind": "contributor", "id": "new-handle"},
                {"kind": "organisation", "id": "example-organisation"},
            ],
        })
