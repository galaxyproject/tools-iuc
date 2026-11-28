from __future__ import annotations

from contextlib import closing
from unittest import mock

from pipeline.history_contributions import contributor_facts
from pipeline.history_index import SCHEMA, HistoryIndex
from tests.helpers import GitTestCase


class HistoryIndexTests(GitTestCase):
    def setUp(self) -> None:
        super().setUp()
        self.write(".gitignore", ".cache/\n")
        self.git("config", "user.name", "First Author")
        self.git("config", "user.email", "123+first@users.noreply.github.com")
        self.path = "tools/example/tool.xml"
        self.write("tools/example/macros.xml", '<macros><token name="@TOOL_VERSION@">1.0</token></macros>')
        self.write(
            self.path,
            """<tool id="example" name="Example" version="@TOOL_VERSION@">
<macros><import>macros.xml</import></macros></tool>""",
        )
        self.added = self.commit("Add tool")
        self.cache = self.repo / ".cache/history.sqlite"

    def tools(self, path: str | None = None) -> list[dict]:
        path = path or self.path
        return [
            {
                "id": "example",
                "name": "Example",
                "source_path": path,
                "tool_dir": path.rsplit("/", 1)[0],
                "macro_paths": ["tools/example/macros.xml"],
            }
        ]

    def index(self) -> HistoryIndex:
        return HistoryIndex(self.repo, self.cache)

    def assert_matches_bootstrap(self, index: HistoryIndex, tools: list[dict] | None = None) -> None:
        tools = tools or self.tools()
        with closing(HistoryIndex(self.repo, self.repo / ".cache/fresh.sqlite", rebuild=True)) as fresh:
            fresh.update()
            self.assertEqual(index.events(tools), fresh.events(tools))
            self.assertEqual(contributor_facts(index, tools), contributor_facts(fresh, tools))

    def test_initial_scan_and_same_head_make_no_more_history_reads(self) -> None:
        with closing(self.index()) as index:
            self.assertEqual(index.update(), 1)
            self.assertEqual(index.events(self.tools())[0]["version"], "1.0")
            with mock.patch("pipeline.history_index.read_commits", side_effect=AssertionError("history scanned")):
                with mock.patch("pipeline.history_index.GitReader", side_effect=AssertionError("files read")):
                    self.assertEqual(index.update(), 0)
        with closing(self.index()) as restored:
            with mock.patch(
                "pipeline.history_index.read_commits", side_effect=AssertionError("restored cache scanned")
            ):
                self.assertEqual(restored.update(), 0)

    def test_only_new_commits_are_extracted_and_macro_changes_credit_the_tool(self) -> None:
        with closing(self.index()) as index:
            index.update()
            self.git("config", "user.name", "Second Author")
            self.git("config", "user.email", "private@example.invalid")
            self.write("tools/example/macros.xml", '<macros><token name="@TOOL_VERSION@">2.0</token></macros>')
            updated = self.commit("Update macro")
            from pipeline.history_git import read_commits

            with mock.patch("pipeline.history_index.read_commits", wraps=read_commits) as scan:
                self.assertEqual(index.update(), 1)
                self.assertEqual(scan.call_args.args, (self.repo, updated, self.added))
            self.assertEqual({e["version"] for e in index.events(self.tools())}, {"1.0", "2.0"})
            facts = contributor_facts(index, self.tools())
            second = next(author for author in facts if author["name"] == "Second Author")
            self.assertEqual(second["contributions"][self.path]["commits"], {updated})
            self.assertIsNone(second["github"])
            self.assertNotIn("private@example.invalid", "\n".join(index.db.iterdump()))
            self.assert_matches_bootstrap(index)

    def test_shed_only_changes_do_not_count_even_from_existing_cache(self) -> None:
        with closing(self.index()) as index:
            index.update()
            before = index.events(self.tools())
            self.git("config", "user.name", "Metadata Author")
            self.git("config", "user.email", "metadata@users.noreply.github.com")
            self.write("tools/example/.shed.yml", "name: example\nowner: iuc\n")
            metadata = self.commit("Change publishing metadata only")
            self.git("rm", "tools/example/.shed.yml")
            deleted = self.commit("Remove publishing metadata")
            index.update()
            self.assertEqual(
                index.db.execute("SELECT count(*) FROM changes WHERE sha IN (?, ?)", (metadata, deleted)).fetchone()[0],
                2,
            )
            with mock.patch("pipeline.history_index.read_commits", side_effect=AssertionError("history scanned")):
                index.update()
                facts = contributor_facts(index, self.tools())
            self.assertEqual(len(facts), 1)
            self.assertEqual(facts[0]["commit_ids"], {self.added})
            self.assertEqual(index.events(self.tools()), before)

    def test_mixed_commit_credits_other_files_without_shed_tool_links(self) -> None:
        other = "tools/other/tool.xml"
        self.write(other, '<tool id="other" name="Other" version="1.0"/>')
        self.commit("Add other tool")
        with closing(self.index()) as index:
            index.update()
            self.git("config", "user.name", "Mixed Author")
            self.git("config", "user.email", "mixed@users.noreply.github.com")
            self.write("tools/example/.shed.yml", "categories: [Metagenomics]\n")
            self.write(other, '<tool id="other" name="Other" version="2.0"/>')
            mixed = self.commit("Edit wrapper and unrelated publishing metadata")
            index.update()
            tools = [*self.tools(), {"id": "other", "name": "Other", "source_path": other}]
            author = next(a for a in contributor_facts(index, tools) if a["name"] == "Mixed Author")
            self.assertEqual(author["commit_ids"], {mixed})
            self.assertEqual(author["contributions"], {other: {"commits": {mixed}, "paths": {other}}})
            self.assertEqual(next(e for e in index.events(tools) if e["version"] == "2.0")["commit_hash"], mixed)

    def test_rename_keeps_original_authors_and_deletion_hides_removed_tool(self) -> None:
        with closing(self.index()) as index:
            index.update()
            self.git("config", "user.name", "Mover")
            self.git("config", "user.email", "mover@users.noreply.github.com")
            moved = "tool_collections/example/renamed.xml"
            self.write("tool_collections/example/.shed.yml", "name: example\nowner: iuc\n")
            self.git("mv", self.path, moved)
            self.commit("Move wrapper")
            index.update()
            facts = contributor_facts(index, self.tools(moved))
            self.assertTrue(all(moved in author["contributions"] for author in facts))
            self.assert_matches_bootstrap(index, self.tools(moved))
            self.git("rm", moved)
            self.commit("Remove wrapper")
            index.update()
            self.assertEqual(index.events([]), [])
            self.assertTrue(all(not author["contributions"] for author in contributor_facts(index, [])))

    def test_rewrite_removes_unreachable_contributions_and_versions(self) -> None:
        with closing(self.index()) as index:
            index.update()
            self.git("config", "user.name", "Discarded Author")
            self.write("tools/example/macros.xml", '<macros><token name="@TOOL_VERSION@">9.0</token></macros>')
            removed = self.commit("Will be rewritten")
            index.update()
            self.git("reset", "--hard", self.added)
            self.git("config", "user.name", "Replacement Author")
            self.write("tools/example/macros.xml", '<macros><token name="@TOOL_VERSION@">3.0</token></macros>')
            self.commit("Replacement")
            index.update()
            self.assertEqual(index.db.execute("SELECT active FROM commits WHERE sha = ?", (removed,)).fetchone()[0], 0)
            self.assertEqual({e["version"] for e in index.events(self.tools())}, {"1.0", "3.0"})
            self.assertNotIn("Discarded Author", {a["name"] for a in contributor_facts(index, self.tools())})
            self.assert_matches_bootstrap(index)

    def test_merge_credits_original_author_once_and_updates_first_parent_stats(self) -> None:
        with closing(self.index()) as index:
            index.update()
            main = self.git("symbolic-ref", "--short", "HEAD")
            self.git("checkout", "--quiet", "-b", "feature")
            self.git("config", "user.name", "Feature Author")
            self.git("config", "user.email", "feature@users.noreply.github.com")
            self.write("tools/example/macros.xml", '<macros><token name="@TOOL_VERSION@">2.0</token></macros>')
            changed = self.commit("Feature version")
            self.git("checkout", "--quiet", main)
            self.git("config", "user.name", "Merger")
            self.git("merge", "--quiet", "--no-ff", "feature", "-m", "Merge feature")
            merged = self.git("rev-parse", "HEAD")
            index.update()
            author = next(a for a in contributor_facts(index, self.tools()) if a["name"] == "Feature Author")
            self.assertEqual(author["contributions"][self.path]["commits"], {changed})
            update = next(e for e in index.events(self.tools()) if e["version"] == "2.0")
            self.assertEqual(update["commit_hash"], merged)
            self.assert_matches_bootstrap(index)

    def test_new_tool_after_checkpoint_has_no_missing_past_events(self) -> None:
        with closing(self.index()) as index:
            index.update()
            self.write("tools/second/tool.xml", '<tool id="second" name="Second" version="1.0"/>')
            added = self.commit("Add second tool")
            index.update()
            tools = [*self.tools(), {"id": "second", "name": "Second", "source_path": "tools/second/tool.xml"}]
            self.assertEqual(next(e for e in index.events(tools) if e["tool_id"] == "second")["commit_hash"], added)
            self.assert_matches_bootstrap(index, tools)

    def test_failed_update_keeps_previous_checkpoint(self) -> None:
        with closing(self.index()) as index:
            index.update()
            self.write("tools/example/macros.xml", '<macros><token name="@TOOL_VERSION@">2.0</token></macros>')
            self.commit("Update")
            with mock.patch("pipeline.history_index.observe_commit", side_effect=RuntimeError("failed")):
                with self.assertRaisesRegex(RuntimeError, "failed"):
                    index.update()
            self.assertEqual(index.checkpoint(), self.added)
            self.assertEqual(index.db.execute("SELECT count(*) FROM commits").fetchone()[0], 1)
            index.update()
            self.assert_matches_bootstrap(index)

    def test_schema_change_and_corrupt_cache_bootstrap_again(self) -> None:
        with closing(self.index()) as index:
            index.update()
            index.db.execute(f"PRAGMA user_version = {SCHEMA + 1}")
        with closing(self.index()) as index:
            self.assertIsNone(index.checkpoint())
            index.update()
        self.cache.write_bytes(b"not a database")
        with closing(self.index()) as index:
            self.assertIsNone(index.checkpoint())
            index.update()

    def test_missing_old_git_object_after_force_push_still_uses_cached_ancestor(self) -> None:
        with closing(self.index()) as index:
            index.update()
            self.write("tools/example/macros.xml", '<macros><token name="@TOOL_VERSION@">9.0</token></macros>')
            removed = self.commit("Discard this commit")
            index.update()
            self.git("reset", "--hard", self.added)
            self.write("tools/example/macros.xml", '<macros><token name="@TOOL_VERSION@">2.0</token></macros>')
            current = self.commit("Replacement commit")
            self.git("reflog", "expire", "--expire=now", "--all")
            self.git("prune", "--expire=now")
            from pipeline.history_git import read_commits

            with mock.patch("pipeline.history_index.read_commits", wraps=read_commits) as scan:
                index.update()
                self.assertEqual(scan.call_args.args, (self.repo, current, self.added))
            self.assertEqual(index.db.execute("SELECT active FROM commits WHERE sha = ?", (removed,)).fetchone()[0], 0)
            self.assert_matches_bootstrap(index)

    def test_shared_macro_imports_and_support_files_deduplicate_commits(self) -> None:
        self.write("macros/shared.xml", '<macros><token name="@TOOL_VERSION@">1.0</token></macros>')
        self.write("tools/example/macros.xml", "<macros><import>../../macros/shared.xml</import></macros>")
        self.commit("Use shared macro")
        with closing(self.index()) as index:
            index.update()
            self.git("config", "user.name", "Shared Author")
            self.write("macros/shared.xml", '<macros><token name="@TOOL_VERSION@">2.0</token></macros>')
            self.write("tools/example/test-data/input.txt", "input fixture")
            self.write("tools/example/helper.py", "# script")
            changed = self.commit("Update macro and tests")
            index.update()
            tools = self.tools()
            tools[0]["macro_paths"].append("macros/shared.xml")
            author = next(a for a in contributor_facts(index, tools) if a["name"] == "Shared Author")
            self.assertEqual(author["commit_ids"], {changed})
            self.assertEqual(author["contributions"][self.path]["commits"], {changed})
            self.assertEqual(len(author["contributions"][self.path]["paths"]), 3)
            self.assertEqual({event["version"] for event in index.events(tools)}, {"1.0", "2.0"})
            self.assert_matches_bootstrap(index, tools)

    def test_linked_macro_target_updates_and_retargeting_are_incremental(self) -> None:
        macro = self.repo / "tools/example/macros.xml"
        macro.unlink()
        macro.symlink_to("../../macros/version.xml")
        self.write("macros/version.xml", '<macros><token name="@TOOL_VERSION@">1.0</token></macros>')
        self.commit("Link shared version macro")
        with closing(self.index()) as index:
            index.update()
            self.git("config", "user.name", "Macro Author")
            self.write("macros/version.xml", '<macros><token name="@TOOL_VERSION@">2.0</token></macros>')
            target_updated = self.commit("Only update symlink target")
            index.update()
            update = next(e for e in index.events(self.tools()) if e["version"] == "2.0")
            self.assertEqual(update["commit_hash"], target_updated)
            self.write("macros/other.xml", '<macros><token name="@TOOL_VERSION@">3.0</token></macros>')
            macro.unlink()
            macro.symlink_to("../../macros/other.xml")
            retargeted = self.commit("Retarget shared macro")
            index.update()
            update = next(e for e in index.events(self.tools()) if e["version"] == "3.0")
            self.assertEqual(update["commit_hash"], retargeted)
            self.assert_matches_bootstrap(index)
            with mock.patch("pipeline.history_index.GitReader", side_effect=AssertionError("historical files read")):
                self.assertEqual(index.update(), 0)

    def test_private_email_as_author_name_is_not_persisted(self) -> None:
        self.git("config", "user.name", "private-name@example.invalid")
        self.git("config", "user.email", "private-name@example.invalid")
        self.write("tools/example/helper.py", "# script")
        self.commit("Support script")
        with closing(self.index()) as index:
            index.update()
            self.assertNotIn("private-name@example.invalid", "\n".join(index.db.iterdump()))
            self.assertEqual(len(contributor_facts(index, self.tools())), 1)

    def test_cached_commit_reintroduced_by_merge_is_reactivated(self) -> None:
        with closing(self.index()) as index:
            index.update()
            self.git("checkout", "--quiet", "-b", "feature")
            self.write("tools/example/helper.py", "# feature")
            feature = self.commit("Feature contribution")
            index.update()
            self.git("checkout", "--quiet", "-b", "replacement", self.added)
            self.write("tools/example/README.md", "Documentation")
            self.commit("Independent contribution")
            index.update()
            self.assertEqual(index.db.execute("SELECT active FROM commits WHERE sha = ?", (feature,)).fetchone()[0], 0)
            self.git("merge", "--quiet", "--no-ff", "feature", "-m", "Merge cached feature")
            index.update()
            self.assertEqual(index.db.execute("SELECT active FROM commits WHERE sha = ?", (feature,)).fetchone()[0], 1)
            self.assert_matches_bootstrap(index)

    def test_ci_output_only_requests_save_when_checkpoint_changes(self) -> None:
        import os

        from pipeline.history_index import main

        output = self.repo / ".cache/github-output"
        output.parent.mkdir(parents=True, exist_ok=True)
        with (
            mock.patch("sys.argv", ["history_index"]),
            mock.patch.dict(os.environ, {"GITHUB_OUTPUT": str(output)}),
            mock.patch("pipeline.history_index.HistoryIndex", side_effect=lambda **kwargs: self.index()),
        ):
            main()
            main()
            self.write("website/style.css", "body {}")
            self.commit("Only website changed")
            main()
        self.assertEqual(output.read_text().splitlines(), ["changed=true", "changed=false", "changed=true"])

    def test_identity_alias_changes_reproject_without_git_extraction(self) -> None:
        from pipeline import extract_people
        from pipeline.common import read_yaml, write_yaml

        with closing(self.index()) as index:
            index.update()
        data = self.repo / ".cache/data"
        write_yaml(data / "tools/index.yaml", [{"id": "example", "owner": "iuc", "repo": "example", "slug": "example"}])
        write_yaml(data / "tools/iuc/example/example.yaml", self.tools()[0])
        with (
            mock.patch.object(extract_people, "DATA_DIR", data),
            mock.patch.object(extract_people, "fetch_gtn_people", return_value={"contributors": {}}),
            mock.patch.object(extract_people, "HistoryIndex", side_effect=lambda: self.index()),
            mock.patch("pipeline.history_index.read_commits", side_effect=AssertionError("history scanned")),
        ):
            with mock.patch.object(extract_people, "_load_aliases", return_value={}):
                extract_people.extract_people()
            first = read_yaml(data / "contributors.yaml")[0]
            with mock.patch.object(extract_people, "_load_aliases", return_value={"canonical": {"github": ["first"]}}):
                extract_people.extract_people()
            canonical = read_yaml(data / "contributors.yaml")[0]
        self.assertEqual((first["id"], canonical["id"]), ("first", "canonical"))
        self.assertEqual(first["contributed_tools"], canonical["contributed_tools"])
        self.assertEqual(canonical["contributed_tools"], [self.path])
        self.assertEqual(first["commits"], canonical["commits"])
        self.assertEqual(canonical["tools"], [])
        self.assertIn("canonical", read_yaml(data / "contributions.yaml"))

    def test_last_updated_dates_come_from_cached_changes(self) -> None:
        with closing(self.index()) as index:
            index.update()
            expected = self.git("log", "-1", "--format=%cI", "--", self.path)
            self.assertEqual(index.last_commit_dates()[self.path], expected)
            with mock.patch("pipeline.history_index.git", side_effect=AssertionError("git called")):
                self.assertEqual(index.last_commit_dates()[self.path], expected)
