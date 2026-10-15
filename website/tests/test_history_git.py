from __future__ import annotations

from contextlib import closing

from pipeline.history_git import GitReader
from tests.helpers import GitTestCase


class GitReaderTests(GitTestCase):
    def test_historical_file_and_directory_links_preserve_all_dependencies(self) -> None:
        self.write("macros/versions.xml", "<macros>historical version</macros>")
        self.write("tools/example/tool.xml", "<tool/>")
        (self.repo / "tools/example/lib").symlink_to("../../macros", target_is_directory=True)
        path = self.repo / "tools/example/macros.xml"
        path.symlink_to("lib/versions.xml")
        initial = self.commit("Add linked macros")
        self.write("macros/versions.xml", "<macros>current checkout</macros>")
        with closing(GitReader(self.repo)) as reader:
            self.assertEqual(reader.read(initial, "tools/example/macros.xml"), "<macros>historical version</macros>")
            self.assertEqual(
                set(reader.dependencies(initial, "tools/example/macros.xml")),
                {
                    "tools/example/macros.xml",
                    "tools/example/lib/versions.xml",
                    "tools/example/lib",
                    "macros/versions.xml",
                },
            )

    def test_broken_cyclic_and_external_links_do_not_resolve(self) -> None:
        self.write("tools/target.xml", "<macros/>")
        for name, target in {
            "broken": "missing.xml",
            "cycle_a": "cycle_b.xml",
            "cycle_b": "cycle_a.xml",
            "escape": "../../outside.xml",
            "absolute": str(self.repo / "tools/target.xml"),
        }.items():
            (self.repo / f"tools/{name}.xml").symlink_to(target)
        sha = self.commit("Add invalid links")
        with closing(GitReader(self.repo)) as reader:
            for name in ("broken", "cycle_a", "cycle_b", "escape", "absolute"):
                with self.subTest(name=name):
                    self.assertIsNone(reader.read(sha, f"tools/{name}.xml"))

    def test_plain_file_is_not_interpreted_as_a_link(self) -> None:
        self.write("tools/plain.xml", "target.xml")
        self.write("tools/target.xml", "<macros/>")
        sha = self.commit("Add ordinary files")
        with closing(GitReader(self.repo)) as reader:
            self.assertEqual(reader.read(sha, "tools/plain.xml"), "target.xml")
