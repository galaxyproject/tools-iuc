from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from pipeline.extract_tools import _render_help


class HelpLinkSanitizerTests(unittest.TestCase):
    def test_rst_disables_raw_html(self) -> None:
        html = _render_help(
            "Normal **help**.\n\n.. raw:: html\n\n    <script>alert('raw')</script>",
            "restructuredtext",
        )

        self.assertIn("<strong>help</strong>", html)
        self.assertNotIn("<script>", html)

    def test_rst_cannot_include_local_files(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "private.txt"
            source.write_text("private-file-content", encoding="utf-8")
            html = _render_help(
                f"Normal **help**.\n\n.. include:: {source}", "restructuredtext",
            )

        self.assertIn("<strong>help</strong>", html)
        self.assertNotIn("private-file-content", html)

    def test_removes_placeholder_image_sources_from_rendered_markdown(self) -> None:
        html = _render_help(
            "![example]($PATH_TO_IMAGES/example.png)\n"
            "![static](${static_path}/images/example.png)\n"
            "[presto](@PRESTO_BASE_URL@/en/stable)",
            "markdown",
        )

        self.assertNotIn("$PATH_TO_IMAGES", html)
        self.assertNotIn("${static_path}", html)
        self.assertNotIn("@PRESTO_BASE_URL@", html)
        self.assertNotIn("src=", html)
        self.assertNotIn("href=", html)

    def test_removes_placeholder_object_data_from_rendered_rst(self) -> None:
        html = _render_help(".. image:: $PATH_TO_IMAGES/example.svg", "restructuredtext")

        self.assertNotIn("$PATH_TO_IMAGES", html)
        self.assertNotIn("data=", html)

    def test_removes_local_wrapper_links_but_keeps_external_links(self) -> None:
        html = _render_help(
            "[paired tool](root?tool_id=other_tool)\n"
            "[bare domain](labs.primalscheme.com)\n"
            "[docs](https://example.org/docs)",
            "markdown",
        )

        self.assertNotIn("root?tool_id", html)
        self.assertNotIn("labs.primalscheme.com", html)
        self.assertIn('href="https://example.org/docs"', html)

    def test_strips_active_content_from_rendered_markdown(self) -> None:
        html = _render_help(
            "<script>alert(1)</script>\n\n"
            '<img src="https://example.org/a.png" onerror="alert(2)">\n\n'
            '<iframe src="https://example.org/frame"></iframe>\n\n'
            "<style>body { display: none }</style>\n\n"
            "[js](javascript:alert(3))\n\n"
            '<a href="https://example.org" onclick="alert(4)">ok</a>',
            "markdown",
        )

        lowered = html.lower()
        for needle in ("<script", "alert(", "onerror", "onclick", "<iframe", "<style", "javascript:"):
            self.assertNotIn(needle, lowered)
        self.assertIn('href="https://example.org"', html)

    def test_rst_raw_directive_is_dropped(self) -> None:
        html = _render_help(
            ".. raw:: html\n\n"
            "   <script>alert(1)</script>\n"
            "   <img src=x onerror=alert(2)>\n\n"
            "After.",
            "restructuredtext",
        )

        lowered = html.lower()
        for needle in ("<script", "alert(", "onerror", "<img"):
            self.assertNotIn(needle, lowered)
        self.assertIn("After.", html)

    def test_rst_file_insertion_directives_do_not_read_local_files(self) -> None:
        marker = "FILE-CONTENT-MARKER-7f3a9c"
        with tempfile.TemporaryDirectory() as directory:
            secret = Path(directory) / "secret.txt"
            secret.write_text(f"{marker}\n", encoding="utf-8")
            table = Path(directory) / "table.csv"
            table.write_text(f"{marker},b\n", encoding="utf-8")

            html = _render_help(
                f".. include:: {secret}\n\n"
                f".. literalinclude:: {secret}\n\n"
                f".. csv-table:: Data\n   :file: {table}\n\n"
                "After.",
                "restructuredtext",
            )

        self.assertNotIn(marker, html)
        self.assertIn("After.", html)

    def test_svg_objects_are_rewritten_to_images_and_other_objects_dropped(self) -> None:
        svg = _render_help(".. image:: https://example.org/a.svg", "restructuredtext")
        self.assertNotIn("<object", svg)
        self.assertIn('src="https://example.org/a.svg"', svg)
        self.assertIn('alt="https://example.org/a.svg"', svg)

        with_alt = _render_help(
            ".. image:: https://example.org/plot.svg\n   :alt: A plot", "restructuredtext"
        )
        self.assertNotIn("<object", with_alt)
        self.assertIn('alt="A plot"', with_alt)

        markdown = _render_help(
            '<object data="https://example.org/b.svg" type="image/svg+xml"></object>\n\n'
            '<object data="https://evil.example/x.html" type="text/html">fallback</object>\n\n'
            '<object data="https://evil.example/y.html">fallback</object>\n\n'
            '<object data="javascript:alert(1)" type="image/svg+xml">fallback</object>',
            "markdown",
        )
        self.assertNotIn("<object", markdown)
        self.assertNotIn("evil.example", markdown)
        self.assertNotIn("javascript:", markdown)
        self.assertIn('<img src="https://example.org/b.svg" alt="b.svg">', markdown)

    def test_strips_javascript_urls_from_rendered_rst(self) -> None:
        html = _render_help("`click <javascript:alert(1)>`_", "restructuredtext")

        self.assertNotIn("javascript:", html.lower())

    def test_keeps_tables_links_and_images(self) -> None:
        html = _render_help(
            "| a | b |\n| - | - |\n| 1 | 2 |\n\n"
            "[docs](https://example.org/docs 'Docs')\n\n"
            "![pic](https://example.org/pic.png)\n\n"
            "```\ncode\n```",
            "markdown",
        )

        for needle in ("<table>", "<th>a</th>", "<td>2</td>", "<pre>", "<code>"):
            self.assertIn(needle, html)
        self.assertIn('href="https://example.org/docs"', html)
        self.assertIn('src="https://example.org/pic.png"', html)

        rst = _render_help(
            "`docs <https://example.org/docs>`_\n\n"
            "=== ===\n"
            "a   b\n"
            "=== ===\n"
            "1   2\n"
            "=== ===\n\n"
            ".. image:: https://example.org/pic.png\n",
            "restructuredtext",
        )
        for needle in ("<table", "<td>2</td>", 'href="https://example.org/docs"', 'src="https://example.org/pic.png"'):
            self.assertIn(needle, rst)

    def test_rst_class_directive_keeps_only_docutils_classes(self) -> None:
        html = _render_help(
            ".. class:: infomark absolute inset-0 z-50\n\nMarked.\n",
            "restructuredtext",
        )

        self.assertIn('<p class="infomark">Marked.</p>', html)
        for needle in ("absolute", "inset-0", "z-50"):
            self.assertNotIn(needle, html)

    def test_markdown_html_cannot_set_page_classes_ids_or_names(self) -> None:
        html = _render_help(
            '<a class="absolute inset-0 z-50" id="theme-toggle" name="main" '
            'href="https://example.org">x</a>\n\n'
            '<div class="fixed" id="main">y</div>',
            "markdown",
        )

        for needle in ("class=", "absolute", "fixed", "name=", 'id="theme-toggle"', 'id="main"'):
            self.assertNotIn(needle, html)
        self.assertIn('id="help-theme-toggle"', html)
        self.assertIn('id="help-main"', html)
        self.assertIn('href="https://example.org"', html)

    def test_rst_section_ids_are_prefixed_and_internal_links_follow(self) -> None:
        html = _render_help(
            "See `Usage`_.\n\nUsage\n=====\n\nText.\n",
            "restructuredtext",
        )

        self.assertIn('id="help-usage"', html)
        self.assertIn('href="#help-usage"', html)
        self.assertNotIn('id="usage"', html)


if __name__ == "__main__":
    unittest.main()
