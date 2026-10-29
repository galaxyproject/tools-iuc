"""Test LIANA resource identity and staged files independently of Galaxy internals."""

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class LianaResourceDownloaderTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.liana = load_module("liana_resource_test_stub", ROOT / "test-data/liana.py")
        with mock.patch.dict("sys.modules", {"liana": cls.liana}):
            cls.module = load_module(
                "liana_resources_downloader_tested",
                ROOT / "data_manager/liana_resources_downloader.py",
            )

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="ldm-")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.extra_files = self.directory / "x"
        self.output = self.directory / "output.json"

    def download(self):
        downloader = self.module.LianaResourcesDownloader(self.output, self.extra_files)
        downloader.run("ligand_receptor_consensus")
        return json.loads(self.output.read_text())["data_tables"]["liana_resources"][0]

    def test_same_content_retains_identity_and_changed_content_gets_new_identity(self):
        first = self.download()
        second = self.download()
        self.assertEqual(first, second)
        changed = self.liana.rs.select_resource("consensus")
        changed.loc[0, "receptor"] = "DIFFERENT_RECEPTOR"
        with mock.patch.object(self.liana.rs, "select_resource", return_value=changed):
            third = self.download()
        self.assertNotEqual(first["value"], third["value"])
        self.assertNotEqual(first["path"], third["path"])
        self.assertRegex(first["value"], r"^consensus_sha256_[a-f0-9]{64}$")
        self.assertEqual(list(self.extra_files.glob(".liana-*.tsv")), [])
        for entry in (first, third):
            file_path = self.extra_files / entry["path"] / f"{entry['value']}.tsv"
            self.assertTrue(file_path.is_file())
            self.assertEqual(set(entry), {"value", "name", "description", "path"})

    def test_existing_artifact_with_wrong_content_is_not_overwritten(self):
        entry = self.download()
        path = self.extra_files / entry["path"] / f"{entry['value']}.tsv"
        path.write_text("unexpected content\n")
        with self.assertRaisesRegex(ValueError, "does not match its content identity"):
            self.download()
        self.assertEqual(path.read_text(), "unexpected content\n")
        self.assertEqual(list(self.extra_files.glob(".liana-*.tsv")), [])


if __name__ == "__main__":
    unittest.main()
