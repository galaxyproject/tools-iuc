from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from galaxy.tool_util.parser.factory import get_tool_source
from galaxy.util.template import fill_template

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().parent


class EmptyTable:
    def get_fields(self):
        return []


class DataManagerWrapperTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        source = get_tool_source(
            str(ROOT / "data_manager/install_trackastra_models.xml")
        )
        cls.command = source.parse_command()
        cls.app = SimpleNamespace(
            tool_data_tables={"trackastra_models": EmptyTable()}
        )

    def render(self, model_source, output_file):
        return fill_template(
            self.command,
            context={
                "__tool_directory__": str(ROOT / "data_manager"),
                "__app__": self.app,
                "output_file": str(output_file),
                "model_source": model_source,
            },
        )

    def flattened_for_galaxy(self, command):
        return " ".join(
            line.strip() for line in command.splitlines() if line.strip()
        )

    def test_both_command_branches_have_valid_shell_syntax(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "output.json"
            branches = [
                {"source": "upstream", "models": "general_2d,ctc"},
                {
                    "source": "archive",
                    "archive": str(DATA / "tiny_test_model.zip"),
                    "value": "tiny_test",
                    "dimensionality": "2D",
                    "model_version": "test",
                },
            ]
            for branch in branches:
                command = self.render(branch, output)
                flattened = self.flattened_for_galaxy(command)
                self.assertNotIn(chr(92) + " --known-models", flattened)
                subprocess.run(["bash", "-n", "-c", command], check=True)
                subprocess.run(["bash", "-n", "-c", flattened], check=True)

    def test_rendered_custom_archive_command(self):
        with tempfile.TemporaryDirectory() as temp:
            temp_path = Path(temp)
            output = temp_path / "output.json"
            output.write_text(
                json.dumps(
                    {
                        "output_data": [
                            {"extra_files_path": str(temp_path / "extra")}
                        ]
                    }
                )
            )
            command = self.render(
                {
                    "source": "archive",
                    "archive": str(DATA / "tiny_test_model.zip"),
                    "value": "tiny_test",
                    "dimensionality": "2D",
                    "model_version": "test",
                },
                output,
            )
            subprocess.run(
                ["bash", "-c", self.flattened_for_galaxy(command)],
                check=True,
            )
            result = json.loads(output.read_text())
            row = result["data_tables"]["trackastra_models"][0]
            self.assertEqual(row["value"], "tiny_test")
            self.assertTrue(Path(row["path"]).is_dir())


if __name__ == "__main__":
    unittest.main()
