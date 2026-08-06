from __future__ import annotations

import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class FrontendContractsTests(unittest.TestCase):
    def test_request_page_fetches_separate_toolshed_duplicate_endpoint(self) -> None:
        request_page = (ROOT / "astro" / "src" / "pages" / "request" / "index.astro").read_text(encoding="utf-8")
        endpoint = (ROOT / "astro" / "src" / "pages" / "toolshed-duplicates.json.ts").read_text(encoding="utf-8")

        self.assertIn("fetch(`${base}/toolshed-duplicates.json`)", request_page)
        self.assertNotIn("toolshed_index.yaml", request_page)
        self.assertIn("loadYaml", endpoint)
        self.assertIn("toolshed_index.yaml", endpoint)

    def test_link_check_uses_fixed_preview_port(self) -> None:
        justfile = (ROOT / "justfile").read_text(encoding="utf-8")
        workflow = (ROOT.parent / ".github" / "workflows" / "website.yml").read_text(encoding="utf-8")

        self.assertIn("--port 4325 --strictPort", justfile)
        self.assertIn("PORT=4325", justfile)
        self.assertIn("--port 4325 --strictPort", workflow)
        self.assertIn("PORT=4325", workflow)

    def test_link_check_builds_pagefind_assets_before_preview(self) -> None:
        justfile = (ROOT / "justfile").read_text(encoding="utf-8")

        self.assertIn("npx pagefind --site astro/dist", justfile)
        self.assertIn("cd astro && (npx astro preview", justfile)


if __name__ == "__main__":
    unittest.main()
