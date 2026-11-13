from __future__ import annotations

import os
from unittest import mock

from pipeline import extract_tools
from pipeline.common import write_yaml
from tests.helpers import GitTestCase


class ToolMetadataCacheTests(GitTestCase):
    def setUp(self) -> None:
        super().setUp()
        self.tool = self.write("tools/example/tool.xml", '<tool id="example" name="Example" version="1.0"/>')
        for name in (
            "tools/example/macro.xml",
            "tools/example/macros_extra.xml",
            "macros/shared.xml",
            "tool_collections/example/macro.xml",
            "data_managers/example/macro.xml",
            "suites/example/macro.xml",
        ):
            self.write(name, "<macros/>")
        self.write("tools/example/.shed.yml", "owner: iuc\nname: example\n")
        self.write("website/pipeline/extract_tools.py", "# extractor\n")
        self.write("website/config/site.yaml", "servers: []\n")
        self.write("website/pyproject.toml", "[project]\nname = 'test'\n")
        self.write("website/uv.lock", "version = 1\n")
        self.conf = self.write("website/tool_conf.xml", f'<toolbox><tool file="{self.tool}"/></toolbox>')
        self.commit("Extraction inputs")
        self.data = self.repo / "website/data"
        overrides = mock.patch.multiple(
            extract_tools, REPO_ROOT=self.repo, WEBSITE_DIR=self.repo / "website", DATA_DIR=self.data,
        )
        overrides.start()
        self.addCleanup(overrides.stop)
        environment = mock.patch.dict(os.environ, {"IUC_FORCE_EXTRACT_TOOLS": ""})
        environment.start()
        self.addCleanup(environment.stop)
        self.signature = self.current_signature()
        write_yaml(self.data / "tools/index.yaml", [{"id": "example"}])
        write_yaml(self.data / "tools/metadata.yaml", {"tool_tree_signature": self.signature})

    def current_signature(self) -> str:
        return extract_tools._tool_tree_signature(self.conf, [(self.tool, None, None)])

    def assert_cache_invalidated(self) -> None:
        signature = self.current_signature()
        self.assertNotEqual(signature, self.signature)
        self.assertFalse(extract_tools._tools_manifest_fresh(signature))

    def test_unchanged_inputs_reuse_the_snapshot_without_galaxy_parsing(self) -> None:
        with mock.patch.object(extract_tools, "SqlAlchemyToolSourceStore") as store:
            self.assertEqual(extract_tools.extract_tools(self.conf), (1, 0))
        store.assert_not_called()

    def test_edits_to_any_source_category_invalidate_the_snapshot(self) -> None:
        for name in (
            "tools/example/tool.xml", "tools/example/macro.xml", "tools/example/macros_extra.xml",
            "macros/shared.xml", "tool_collections/example/macro.xml", "data_managers/example/macro.xml",
            "suites/example/macro.xml", "tools/example/.shed.yml", "website/pipeline/extract_tools.py",
            "website/config/site.yaml", "website/pyproject.toml", "website/uv.lock", "website/tool_conf.xml",
        ):
            with self.subTest(path=name):
                path = self.repo / name
                original = path.read_bytes()
                try:
                    path.write_bytes(original + b"\n")
                    self.assert_cache_invalidated()
                finally:
                    path.write_bytes(original)

    def test_new_tracked_xml_invalidates_the_snapshot(self) -> None:
        self.write("macros/new.xml", "<macros/>")
        self.git("add", "macros/new.xml")
        self.assert_cache_invalidated()

    def test_deleted_tracked_xml_invalidates_the_snapshot(self) -> None:
        (self.repo / "macros/shared.xml").unlink()
        self.assert_cache_invalidated()

    def test_renamed_tracked_xml_invalidates_the_snapshot(self) -> None:
        self.git("mv", "macros/shared.xml", "macros/renamed.xml")
        self.assert_cache_invalidated()

    def test_new_pipeline_and_config_files_invalidate_without_staging(self) -> None:
        for name in ("website/pipeline/new_module.py", "website/config/new_config.yml"):
            with self.subTest(path=name):
                path = self.write(name, "# new input\n")
                try:
                    self.assert_cache_invalidated()
                finally:
                    path.unlink()

    def test_generated_data_and_bytecode_do_not_invalidate_the_snapshot(self) -> None:
        self.write("website/data/generated.xml", "<generated/>")
        self.write("website/pipeline/__pycache__/extract_tools.pyc", "bytecode")
        self.assertEqual(self.current_signature(), self.signature)
        self.assertTrue(extract_tools._tools_manifest_fresh(self.signature))

    def test_head_changes_refresh_git_derived_metadata(self) -> None:
        self.write("README.md", "Documentation change\n")
        self.commit("New HEAD")
        self.assert_cache_invalidated()

    def test_changed_macro_enters_the_rebuild_path(self) -> None:
        self.write("tools/example/macro.xml", '<macros><token name="@VERSION@">2</token></macros>')
        with mock.patch.object(extract_tools, "SqlAlchemyToolSourceStore", side_effect=RuntimeError("rebuild")):
            with self.assertRaisesRegex(RuntimeError, "rebuild"):
                extract_tools.extract_tools(self.conf)

    def test_force_flag_bypasses_an_unchanged_snapshot(self) -> None:
        with mock.patch.dict(os.environ, {"IUC_FORCE_EXTRACT_TOOLS": "1"}):
            self.assertFalse(extract_tools._tools_manifest_fresh(self.signature))
