from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest import mock

from pipeline import extract_people as extract_people_module
from pipeline.common import http_url, http_url_or_handle, read_yaml, write_yaml
from pipeline.extract_people import _creator_website, _organisation_url, extract_people
from pipeline.extract_tools import _detail_fields
from pipeline.fetch_gtn import _clean_contributors, _clean_organisations

UNSAFE = (
    "javascript:alert(1)",
    " JavaScript:alert(1)",
    "java\tscript:alert(1)",
    "data:text/html,x",
    "vbscript:x",
    "https:example.org",
    "//example.org/path",
    "javascript://example.org/%0Aalert(1)",
    "data://x/y",
    "vbscript://example.org/x",
    "example.org",
    "",
    None,
)


class _CreatorSource:
    def __init__(self, creators: list[dict]) -> None:
        self._creators = creators

    def parse_creator(self) -> list[dict]:
        return self._creators

    def __getattr__(self, name: str):
        def unavailable(*_args, **_kwargs):
            raise NotImplementedError(name)
        return unavailable


class HttpUrlTests(unittest.TestCase):
    def test_keeps_http_and_https_urls(self) -> None:
        self.assertEqual(http_url("https://example.org/a?b=1"), "https://example.org/a?b=1")
        self.assertEqual(http_url(" HTTP://example.org "), "HTTP://example.org")

    def test_drops_other_schemes_and_non_urls(self) -> None:
        for value in UNSAFE:
            with self.subTest(value=value):
                self.assertIsNone(http_url(value))

    def test_handle_values_keep_their_current_form(self) -> None:
        self.assertEqual(http_url_or_handle("0000-0002-1825-0097"), "0000-0002-1825-0097")
        self.assertEqual(http_url_or_handle("@user@mastodon.social"), "@user@mastodon.social")
        self.assertEqual(http_url_or_handle("mastodon.social/@user"), "mastodon.social/@user")
        self.assertEqual(http_url_or_handle("https://orcid.org/0000-0002-1825-0097"), "https://orcid.org/0000-0002-1825-0097")
        self.assertIsNone(http_url_or_handle("javascript:alert(1)"))
        self.assertIsNone(http_url_or_handle(None))


class CreatorUrlTests(unittest.TestCase):
    def test_tool_creator_urls_keep_only_http(self) -> None:
        creators = _detail_fields(_CreatorSource([
            {"class": "Person", "name": "A", "url": "https://example.org/a"},
            {"class": "Person", "name": "B", "url": "javascript:alert(1)"},
            {"class": "Organization", "name": "C", "url": "data:text/html,x"},
        ]))["creators"]

        self.assertEqual(creators[0]["url"], "https://example.org/a")
        self.assertNotIn("url", creators[1])
        self.assertNotIn("url", creators[2])
        self.assertEqual(creators[1]["name"], "B")

    def test_contributor_website_keeps_only_http(self) -> None:
        self.assertEqual(_creator_website({"url": "https://example.org/me"}), "https://example.org/me")
        self.assertIsNone(_creator_website({"url": "https://orcid.org/0000-0002-1825-0097"}))
        self.assertIsNone(_creator_website({"url": "javascript:alert(1)"}))
        self.assertIsNone(_creator_website({}))

    def test_organisation_url_keeps_only_http(self) -> None:
        self.assertEqual(_organisation_url({"url": "https://example.org"}), "https://example.org")
        self.assertIsNone(_organisation_url({"url": "javascript:alert(1)"}))
        self.assertIsNone(_organisation_url({}))


class ExtractPeopleUrlTests(unittest.TestCase):
    def test_unsafe_creator_urls_are_dropped_from_people_and_organisations(self) -> None:
        creators = [
            {"class": "Person", "name": "Safe Person", "identifier": "0000-0002-1825-0097", "url": "https://example.org/me"},
            {"class": "Person", "name": "Unsafe Person", "identifier": "0000-0001-2345-6789", "url": "javascript://example.org/%0Aalert(1)"},
            {"class": "Organization", "name": "Safe Org", "url": "https://example.org"},
            {"class": "Organization", "name": "Unsafe Org", "url": "data://x/y"},
        ]
        with tempfile.TemporaryDirectory() as tmp:
            data_dir = Path(tmp)
            tool_dir = data_dir / "tools" / "owner" / "repo"
            tool_dir.mkdir(parents=True)
            write_yaml(data_dir / "tools" / "index.yaml", [{"owner": "owner", "repo": "repo", "slug": "t"}])
            write_yaml(tool_dir / "t.yaml", {"id": "t", "source_path": "tools/t/t.xml", "creators": creators})

            with (
                mock.patch.object(extract_people_module, "DATA_DIR", data_dir),
                mock.patch.object(extract_people_module, "fetch_gtn_people", return_value={}),
                mock.patch.object(extract_people_module, "_git_authors", return_value=[]),
            ):
                extract_people()

            people = {p["name"]: p for p in read_yaml(data_dir / "contributors.yaml")}
            orgs = {o["name"]: o for o in read_yaml(data_dir / "organisations.yaml")}

        self.assertEqual(people["Safe Person"]["url"], "https://example.org/me")
        self.assertIsNone(people["Unsafe Person"]["url"])
        self.assertEqual(orgs["Safe Org"]["url"], "https://example.org")
        self.assertIsNone(orgs["Unsafe Org"]["url"])


class GtnUrlTests(unittest.TestCase):
    def test_contributor_links_keep_only_http(self) -> None:
        cleaned = _clean_contributors({
            "safe": {
                "name": "Safe", "url": "https://example.org", "orcid": "0000-0002-1825-0097",
                "fediverse": "https://mastodon.social/@safe",
            },
            "handles": {"name": "Handles", "fediverse": "@user@mastodon.social"},
            "unsafe": {
                "name": "Unsafe", "url": "javascript:alert(1)", "website": "data:text/html,x",
                "orcid": "javascript:alert(2)", "fediverse": "vbscript:x",
            },
        })

        self.assertEqual(cleaned["safe"]["url"], "https://example.org")
        self.assertEqual(cleaned["safe"]["orcid"], "0000-0002-1825-0097")
        self.assertEqual(cleaned["safe"]["fediverse"], "https://mastodon.social/@safe")
        self.assertEqual(cleaned["handles"]["fediverse"], "@user@mastodon.social")
        for key in ("url", "website", "orcid", "fediverse"):
            self.assertNotIn(key, cleaned["unsafe"])
        self.assertEqual(cleaned["unsafe"]["name"], "Unsafe")

    def test_organisation_url_keeps_only_http(self) -> None:
        cleaned = _clean_organisations({
            "safe": {"name": "Safe", "url": "https://example.org"},
            "unsafe": {"name": "Unsafe", "url": "javascript:alert(1)"},
        })

        self.assertEqual(cleaned["safe"]["url"], "https://example.org")
        self.assertNotIn("url", cleaned["unsafe"])


if __name__ == "__main__":
    unittest.main()
