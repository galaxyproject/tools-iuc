from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
import zipfile
from pathlib import Path

SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "data_manager"
    / "fetch_models.py"
)
SPEC = importlib.util.spec_from_file_location("fetch_models", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class FetchModelsTest(unittest.TestCase):
    def test_safe_extract_and_install(self) -> None:
        archive = Path(__file__).with_name("tiny_test_model.zip")
        with tempfile.TemporaryDirectory() as temp:
            destination, archive_hash = MODULE.install_archive(
                archive, Path(temp), "tiny_test"
            )
            self.assertEqual(len(archive_hash), 64)
            installed_files = {item.name for item in destination.iterdir()}
            self.assertTrue(MODULE.REQUIRED_FILES.issubset(installed_files))

    def test_rejects_path_traversal(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            archive = Path(temp) / "unsafe.zip"
            with zipfile.ZipFile(archive, "w") as handle:
                handle.writestr("../outside.txt", "bad")
            with self.assertRaisesRegex(ValueError, "Unsafe archive member"):
                MODULE.safe_extract_zip(archive, Path(temp) / "output")

    def test_custom_archive_cli_output(self) -> None:
        archive = Path(__file__).with_name("tiny_test_model.zip")
        with tempfile.TemporaryDirectory() as temp:
            temp_path = Path(temp)
            data_manager_json = temp_path / "output.json"
            extra = temp_path / "extra"
            data_manager_json.write_text(
                json.dumps({"output_data": [{"extra_files_path": str(extra)}]})
            )
            import subprocess
            subprocess.run(
                [
                    "python",
                    str(SCRIPT),
                    str(data_manager_json),
                    "--archive",
                    str(archive),
                    "--value",
                    "tiny_test",
                    "--dimensionality",
                    "2D",
                    "--model-version",
                    "test",
                ],
                check=True,
            )
            result = json.loads(data_manager_json.read_text())
            row = result["data_tables"]["trackastra_models"][0]
            self.assertEqual(row["value"], "tiny_test")
            self.assertEqual(row["source"], "administrator_upload")
            self.assertTrue(Path(row["path"]).is_dir())


if __name__ == "__main__":
    unittest.main()
