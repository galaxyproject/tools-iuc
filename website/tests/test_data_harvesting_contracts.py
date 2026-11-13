from __future__ import annotations

import unittest

from pipeline.common import WEBSITE_DIR, read_yaml


class DataHarvestingContractsTests(unittest.TestCase):
    def test_contributor_aliases_do_not_contain_email_values(self) -> None:
        aliases = read_yaml(WEBSITE_DIR / "config" / "contributor_aliases.yaml") or {}

        for canonical, values in aliases.items():
            self.assertNotIn("@", canonical)
            for field in ("github", "name"):
                for alias in values.get(field, []) or []:
                    self.assertNotIn("@", alias)


if __name__ == "__main__":
    unittest.main()
