from __future__ import annotations

from pipeline.stats.git_analyzer import _build_location_history, get_commit_data
from pipeline.stats.main import run
from tests.helpers import GitTestCase


class HistoryTests(GitTestCase):
    def test_never_modified_additions_are_discovered_in_every_tool_root(self) -> None:
        for root in ("tools", "tool_collections", "data_managers", "suites"):
            self.write(f"{root}/example/.shed.yml", f"name: {root}\nowner: iuc\n")
            self.write(f"{root}/example/tool.xml", f'<tool id="{root}" name="Example" version="1.0"/>')
        added = self.commit("Add tools")
        events = run(repo=str(self.repo), output_dir=str(self.repo / "output"))
        self.assertEqual({event.tool_id for event in events}, {"tools", "tool_collections", "data_managers", "suites"})
        self.assertTrue(all(event.event_type == "new" and event.commit_hash == added for event in events))

    def test_updates_follow_a_tool_renamed_between_roots(self) -> None:
        original = "tools/example/old.xml"
        self.write("tools/example/.shed.yml", "name: example\nowner: iuc\n")
        self.write(original, '<tool id="example" name="Example" version="1.0"/>')
        added = self.commit("Add tool")
        self.write(original, '<tool id="example" name="Example" version="1.1"/>')
        updated = self.commit("Update tool")
        moved = "tool_collections/example/new.xml"
        self.write("tool_collections/example/.shed.yml", "name: example\nowner: iuc\n")
        self.git("mv", original, moved)
        self.commit("Move tool")
        self.write(moved, '<tool id="example" name="Example" version="2.0"/>')
        updated_again = self.commit("Update moved tool")
        events = run(repo=str(self.repo), output_dir=str(self.repo / "output"))
        self.assertEqual([event.version for event in events], ["1.0", "1.1", "2.0"])
        self.assertEqual([event.commit_hash for event in events], [added, updated, updated_again])

    def test_deleted_files_are_removed_from_location_history(self) -> None:
        self.write("tools/example/tool.xml", '<tool id="example" name="Example" version="1.0"/>')
        added = self.commit("Add tool")
        self.git("rm", "tools/example/tool.xml")
        deleted = self.commit("Remove tool")
        commits = list(reversed(get_commit_data(str(self.repo), "HEAD")))
        locations = _build_location_history(str(self.repo), commits, {"example"})
        self.assertIn("example", locations[added])
        self.assertNotIn("example", locations[deleted])

    def test_merge_additions_are_counted_on_the_first_parent_history(self) -> None:
        self.write("tools/example/.shed.yml", "name: example\nowner: iuc\n")
        self.commit("Initialize publishing config")
        main = self.git("symbolic-ref", "--short", "HEAD")
        self.git("checkout", "--quiet", "-b", "feature")
        self.write("tools/example/tool.xml", '<tool id="example" name="Example" version="1.0"/>')
        self.commit("Add tool on feature branch")
        self.git("checkout", "--quiet", main)
        self.git("merge", "--quiet", "--no-ff", "feature", "-m", "Merge new tool")
        merged = self.git("rev-parse", "HEAD")
        events = run(repo=str(self.repo), output_dir=str(self.repo / "output"))
        self.assertEqual(len(events), 1)
        self.assertEqual((events[0].event_type, events[0].commit_hash), ("new", merged))
