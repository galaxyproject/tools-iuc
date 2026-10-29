"""Exercise the actual XML regex validators at their Galaxy matching boundary."""

import re
import unittest
from pathlib import Path
from xml.etree import ElementTree


MACROS = Path(__file__).resolve().parents[1] / "macros.xml"


class SingleLineValidatorTests(unittest.TestCase):
    def test_keys_and_csv_reject_line_breaks(self):
        root = ElementTree.parse(MACROS)
        for name in [
            "validator_required_key",
            "validator_optional_key",
            "validator_required_csv",
            "validator_optional_csv",
        ]:
            node = root.find(f".//xml[@name='{name}']/validator[@type='regex']")
            pattern = re.compile(node.text)
            for invalid in ["A\nB", "A\rB", "A\n", "A\r\n", "\n"]:
                with self.subTest(validator=name, value=repr(invalid)):
                    self.assertIsNone(pattern.match(invalid))
            with self.subTest(validator=name, value="internal spaces"):
                self.assertIsNotNone(pattern.match("Group with spaces"))
            if name.endswith("csv"):
                self.assertIsNotNone(pattern.match("Group with spaces,B"))
            with self.subTest(validator=name, value="empty"):
                if "optional" in name:
                    self.assertIsNotNone(pattern.match(""))
                else:
                    self.assertIsNone(pattern.match(""))


if __name__ == "__main__":
    unittest.main()
