from __future__ import annotations

import unittest
from html.parser import HTMLParser
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = "/tools-iuc/"


class AssetReferenceParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.assets: set[str] = set()

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attr_map = dict(attrs)
        if tag in {"img", "script"} and attr_map.get("src"):
            self.assets.add(attr_map["src"] or "")
        if tag == "link" and attr_map.get("href"):
            rel = (attr_map.get("rel") or "").lower()
            if "stylesheet" in rel or "icon" in rel:
                self.assets.add(attr_map["href"] or "")


class FrontendContractsTests(unittest.TestCase):
    def test_request_page_fetches_separate_toolshed_duplicate_endpoint(self) -> None:
        request_page = (ROOT / "astro" / "src" / "pages" / "request" / "index.astro").read_text(encoding="utf-8")
        endpoint = (ROOT / "astro" / "src" / "pages" / "toolshed-duplicates.json.ts").read_text(encoding="utf-8")

        self.assertIn("fetch(`${base}/toolshed-duplicates.json`)", request_page)
        self.assertNotIn("toolshed_index.yaml", request_page)
        self.assertIn("loadYaml", endpoint)
        self.assertIn("toolshed_index.yaml", endpoint)
        self.assertNotIn("description", endpoint)

    def test_tools_listing_does_not_duplicate_visible_search_text_in_data_attributes(self) -> None:
        page = (ROOT / "astro" / "src" / "pages" / "tools" / "index.astro").read_text(encoding="utf-8")

        self.assertNotIn("data-name=", page)
        self.assertNotIn("data-desc=", page)
        self.assertIn("li.textContent", page)
        self.assertIn("data-id=", page)

    def test_generated_tool_yaml_keeps_rendered_help_only(self) -> None:
        extractor = (ROOT / "pipeline" / "extract_tools.py").read_text(encoding="utf-8")
        data_types = (ROOT / "astro" / "src" / "lib" / "data.ts").read_text(encoding="utf-8")

        self.assertIn('"help_html": help_html', extractor)
        self.assertNotIn('"help": help_text', extractor)
        self.assertNotIn("help: string;", data_types)

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

    def test_fork_pages_deploy_does_not_use_protected_environment(self) -> None:
        workflow = (ROOT.parent / ".github" / "workflows" / "website.yml").read_text(encoding="utf-8")
        canonical_job = workflow.split("  deploy-pages:", 1)[1].split("  deploy-pages-fork:", 1)[0]
        fork_job = workflow.split("  deploy-pages-fork:", 1)[1]

        self.assertIn("github.repository == 'galaxyproject/tools-iuc'", canonical_job)
        self.assertIn("environment:", canonical_job)
        self.assertIn("timeout: 600000", canonical_job)
        self.assertIn("github.repository == 'bgruening/tools-iuc'", fork_job)
        self.assertNotIn("environment:", fork_job)
        self.assertIn("timeout: 600000", fork_job)

    def test_workflow_does_not_publish_generated_data_to_gh_pages(self) -> None:
        workflow = (ROOT.parent / ".github" / "workflows" / "website.yml").read_text(encoding="utf-8")

        self.assertNotIn("commit-data:", workflow)
        self.assertNotIn("Download generated data", workflow)
        self.assertNotIn("git push origin HEAD:gh-pages", workflow)

    def test_contributor_tool_cards_render_edam_labels(self) -> None:
        page = (ROOT / "astro" / "src" / "pages" / "contributors" / "[id].astro").read_text(encoding="utf-8")

        self.assertIn("const edamLabel", page)
        self.assertIn("{edamLabel(e)}", page)
        self.assertNotIn('badge badge-slate">{e}</span>', page)

    def test_tool_creator_badges_link_to_hall_of_fame_when_known(self) -> None:
        page = (ROOT / "astro" / "src" / "pages" / "tools" / "[owner]" / "[repo]" / "[id].astro").read_text(
            encoding="utf-8"
        )

        self.assertIn("contributorIdForCreator", page)
        self.assertIn("organisationIdForCreator", page)
        self.assertIn("contributorLinks.byGithub", page)
        self.assertIn("contributorLinks.byOrcid", page)
        self.assertIn("contributorLinks.byName", page)
        self.assertIn('href={`${base}/contributors/${cid}/`}', page)
        self.assertIn('href={`${base}/organisations/${oid}/`}', page)

    def test_pagefind_indexes_tool_metadata_not_full_help(self) -> None:
        page = (ROOT / "astro" / "src" / "pages" / "tools" / "[owner]" / "[repo]" / "[id].astro").read_text(
            encoding="utf-8"
        )

        self.assertIn("data-pagefind-body", page)
        self.assertIn('data-pagefind-ignore="all"', page)

    def test_organisation_pages_use_generated_ids(self) -> None:
        detail = (ROOT / "astro" / "src" / "pages" / "organisations" / "[id].astro").read_text(encoding="utf-8")
        index = (ROOT / "astro" / "src" / "pages" / "organisations" / "index.astro").read_text(encoding="utf-8")

        self.assertIn("params: { id: o.id }", detail)
        self.assertIn('href={`${base}/organisations/${o.id}/`}', index)
        self.assertNotIn("o.name.toLowerCase().replace", detail + index)

    def test_built_pages_reference_existing_local_assets(self) -> None:
        dist = ROOT / "astro" / "dist"
        if not dist.exists():
            self.skipTest("Astro dist/ is not built")
        pages = list(dist.rglob("*.html"))
        if not pages:
            self.skipTest("Astro dist/ does not contain built HTML pages")

        missing: list[str] = []
        checked: set[str] = set()
        for page in pages:
            parser = AssetReferenceParser()
            parser.feed(page.read_text(encoding="utf-8"))
            for asset in parser.assets:
                if not asset.startswith(BASE_PATH):
                    continue
                asset_path = asset.removeprefix(BASE_PATH).split("?", 1)[0].split("#", 1)[0]
                if asset_path.startswith("pagefind/") and not (dist / "pagefind").exists():
                    continue
                checked.add(asset_path)
                if not (dist / asset_path).is_file():
                    missing.append(f"{page.relative_to(dist)} -> {asset}")

        self.assertGreater(len(checked), 0)
        self.assertEqual(missing, [])


if __name__ == "__main__":
    unittest.main()
