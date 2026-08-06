from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest import mock

from pipeline import extract_stats, extract_tools
from pipeline.common import md5_files, write_yaml


class CacheTests(unittest.TestCase):
    def test_mulled_v2_container_link_targets_exact_quay_tag(self) -> None:
        requirements = [
            {"type": "package", "name": "bcftools", "version": "1.17"},
            {"type": "package", "name": "samtools", "version": "1.17"},
        ]

        links = extract_tools._compute_container_links(requirements)
        mulled_name = extract_tools._mulled_v2_name(requirements)
        self.assertIsNotNone(mulled_name)
        repository, tag = mulled_name.split(":", 1)

        self.assertEqual(
            links["docker"],
            f"https://quay.io/repository/biocontainers/{repository}?tab=tags&tag={tag}",
        )
        self.assertEqual(
            links["singularity"],
            f"https://depot.galaxyproject.org/singularity/{mulled_name}",
        )

    def test_edam_terms_are_enriched_with_human_readable_labels(self) -> None:
        previous_labels = extract_tools._EDAM_LABELS
        extract_tools._EDAM_LABELS = {
            "example-tool": {
                "topic_0622": "Genomics",
                "operation_3198": "Genome assembly",
            }
        }
        self.addCleanup(setattr, extract_tools, "_EDAM_LABELS", previous_labels)
        entry = mock.Mock(
            id="example",
            name="Example",
            version="1.0",
            description="Example tool",
            profile="24.0",
            tool_type="default",
            hidden=False,
            edam_topics=["topic_0622"],
            edam_operations=["operation_3198"],
            xrefs=[{"type": "bio.tools", "value": "example-tool"}],
            panel_section_id=None,
            panel_section_name=None,
            test_count=0,
            requirements=[],
            container_requirements=[],
            tags=[],
            source_path="tools/example/example.xml",
        )

        ts = mock.Mock()
        ts.parse_citations.return_value = []

        with (
            mock.patch.object(extract_tools, "_detail_fields", return_value={
                "license": None,
                "icon": None,
                "creators": [],
                "help_html": "",
                "help_format": None,
                "inputs": [],
                "outputs": [],
            }),
            mock.patch.object(extract_tools, "owner_repo_from_shed_yml", return_value=("owner", "repo")),
        ):
            full = extract_tools._entry_to_full(entry, Path.cwd(), ts)

        self.assertEqual(
            full["edam_topics"],
            [{"uri": "http://edamontology.org/topic_0622", "label": "Genomics"}],
        )
        self.assertEqual(
            full["edam_operations"],
            [{"uri": "http://edamontology.org/operation_3198", "label": "Genome assembly"}],
        )

    def test_tool_slug_collision_is_rejected_within_owner_repo(self) -> None:
        slug_keys: dict[tuple[str, str, str], str] = {}

        self.assertEqual(
            extract_tools._check_slug_collision(slug_keys, "owner", "repo", "tool:id"),
            "tool_id",
        )
        with self.assertRaises(ValueError):
            extract_tools._check_slug_collision(slug_keys, "owner", "repo", "tool/id")

    def test_tool_slug_collision_allows_same_slug_in_different_repos(self) -> None:
        slug_keys: dict[tuple[str, str, str], str] = {}

        self.assertEqual(
            extract_tools._check_slug_collision(slug_keys, "owner", "first", "tool:id"),
            "tool_id",
        )
        self.assertEqual(
            extract_tools._check_slug_collision(slug_keys, "owner", "second", "tool/id"),
            "tool_id",
        )

    def test_slim_tool_index_uses_full_tool_slug(self) -> None:
        row = extract_tools._slim(
            {
                "id": "tool:id",
                "slug": "tool_id",
                "name": "Tool",
                "version": "1.0",
                "description": "desc",
                "edam_operations": [],
                "edam_topics": [],
                "biotools": None,
                "doi": None,
                "panel_section_name": None,
                "panel_section_id": None,
                "owner": "owner",
                "repo": "repo",
                "source_path": "tools/tool.xml",
                "test_count": 0,
                "updated": None,
                "inputs": [],
                "outputs": [],
            }
        )

        self.assertEqual(row["slug"], "tool_id")

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
