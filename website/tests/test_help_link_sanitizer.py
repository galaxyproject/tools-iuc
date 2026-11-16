from __future__ import annotations

import unittest

from pipeline.extract_tools import _render_help


class HelpLinkSanitizerTests(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
