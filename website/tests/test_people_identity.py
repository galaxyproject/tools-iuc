from __future__ import annotations

import unittest

from pipeline.extract_people import _build_alias_maps, _canonical_github, _merge_local_key, _resolve_key


class PeopleIdentityTests(unittest.TestCase):
    def test_local_identity_does_not_merge_without_public_identifier(self) -> None:
        orcid_keys: dict[str, str] = {}

        self.assertEqual(
            _merge_local_key("github-handle", None, orcid_keys),
            "github-handle",
        )
        self.assertEqual(
            _merge_local_key("person-name", None, orcid_keys),
            "person-name",
        )

    def test_local_orcid_key_reuses_existing_contributor(self) -> None:
        orcid_keys: dict[str, str] = {}

        self.assertEqual(
            _merge_local_key("plain-name", None, orcid_keys, orcid="0000-0001-6431-3442"),
            "plain-name",
        )
        self.assertEqual(
            _merge_local_key("accented-name", None, orcid_keys, orcid="0000-0001-6431-3442"),
            "plain-name",
        )

    def test_gtn_identity_is_not_rewritten_by_local_maps(self) -> None:
        orcid_keys = {"0000-0001-6431-3442": "unmatched-person"}

        self.assertEqual(
            _merge_local_key("gtn-handle", "gtn-handle", orcid_keys, orcid="0000-0001-6431-3442"),
            "gtn-handle",
        )

    def test_reviewed_name_alias_resolves_to_canonical_key(self) -> None:
        alias_gh, alias_name = _build_alias_maps({"canonical-handle": {"github": [], "name": ["legacy-name"]}})

        self.assertEqual(
            _resolve_key(
                {"by_handle": {}, "by_name": {}, "by_orcid": {}},
                alias_gh,
                alias_name,
                name="Legacy Name",
            ),
            ("canonical-handle", None),
        )

    def test_reviewed_github_alias_exposes_canonical_handle(self) -> None:
        self.assertEqual(
            _canonical_github("canonical-handle", "legacy-handle", "Contributor Name"),
            "canonical-handle",
        )
