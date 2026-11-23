from __future__ import annotations

import xml.etree.ElementTree as ET
from unittest.mock import patch

from pipeline import build_tool_conf
from pipeline.inventory import TOOL_ROOTS, published_tool_files, tool_files
from pipeline.repositories import repository_for_tool
from pipeline.stats.discovery import discover_tools
from tests.helpers import GitTestCase

TOOL = '<tool id="{id}" name="Example" version="1.0"/>'


class InventoryTests(GitTestCase):
    def test_metadata_and_history_use_the_same_tracked_files_in_all_roots(self) -> None:
        expected = []
        for root in TOOL_ROOTS:
            expected.append(self.write(f"{root}/example/tool.xml", TOOL.format(id=root)))
            self.write(f"{root}/example/.shed.yml", f"name: {root}\nowner: iuc\n")
        self.write("tools/example/test-data/fixture.xml", TOOL.format(id="fixture"))
        self.write("tools/example/macros.xml", "<macros/>")
        self.commit("Add tools")
        self.write("tools/example/untracked.xml", TOOL.format(id="untracked"))

        self.assertEqual(tool_files(self.repo), sorted(expected))
        with patch.object(build_tool_conf, "REPO_ROOT", self.repo):
            conf = build_tool_conf.build_tool_conf(self.repo / "tool_conf.xml")
        configured = {element.get("file") for element in ET.parse(conf).findall(".//tool")}
        self.assertEqual(configured, {str(path) for path in expected})
        self.assertEqual({tool.tool_id for tool in discover_tools(str(self.repo))}, set(TOOL_ROOTS))


class RepositoryTests(GitTestCase):
    def test_single_repository_uses_declared_name_and_owner(self) -> None:
        tool = self.write("tools/example/tool.xml", TOOL.format(id="Example"))
        self.write("tools/example/.shed.yml", "name: published_name\nowner: devteam\n")
        self.assertEqual(repository_for_tool(tool), ("devteam", "published_name"))

    def test_auto_repositories_expand_templates_and_ignore_suite_and_shared_macros(self) -> None:
        tool = self.write("tools/example/one.xml", TOOL.format(id="UPPERCASE"))
        other = self.write("tools/example/two.xml", TOOL.format(id="other"))
        self.write("tools/example/macros.xml", "<macros/>")
        self.write("tools/example/.shed.yml", """owner: iuc
auto_tool_repositories:
  name_template: 'prefix_{{ tool_id }}'
suite:
  name: suite_example
""")
        self.assertEqual(repository_for_tool(tool), ("iuc", "prefix_uppercase"))
        self.assertEqual(repository_for_tool(other), ("iuc", "prefix_other"))

    def test_explicit_repositories_apply_includes_excludes_and_owner_overrides(self) -> None:
        one = self.write("tool_collections/example/one.xml", TOOL.format(id="one"))
        two = self.write("tool_collections/example/nested/two.xml", TOOL.format(id="two"))
        excluded = self.write("tool_collections/example/private.xml", TOOL.format(id="private"))
        self.write("tool_collections/example/.shed.yml", """owner: iuc
repositories:
  first:
    include: ['*.xml']
    exclude: ['private.xml']
  second:
    owner: devteam
    include:
      - source: nested/two.xml
        destination: tool.xml
""")
        self.assertEqual(repository_for_tool(one), ("iuc", "first"))
        self.assertEqual(repository_for_tool(two), ("devteam", "second"))
        self.assertIsNone(repository_for_tool(excluded))
        self.commit("Add publishing config")
        self.assertEqual(published_tool_files(self.repo), sorted([one, two]))
