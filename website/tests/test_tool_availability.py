from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from pipeline import extract_tool_availability as availability
from pipeline.common import read_yaml, write_yaml


class ToolAvailabilityTests(unittest.TestCase):
    def setUp(self) -> None:
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name)
        self.output = self.root / "tool_availability.yaml"
        write_yaml(self.root / "config/site.yaml", {
            "servers": [
                {"name": "EU", "url": "https://eu.example.org"},
                {"name": "US", "url": "https://us.example.org"},
            ],
        })
        for name in ("DATA_DIR", "WEBSITE_DIR"):
            patcher = patch.object(availability, name, self.root)
            patcher.start()
            self.addCleanup(patcher.stop)

    def test_complete_refresh_is_cached(self) -> None:
        with patch.object(availability, "_fetch_installed", return_value=({"tool"}, {"iuc/repo"})) as fetch:
            result = availability.extract_tool_availability()
            self.assertTrue(availability._is_cache_fresh(self.output))
            self.assertEqual(availability.extract_tool_availability(), result)

        self.assertEqual(fetch.call_count, 2)
        self.assertEqual(result["servers"]["EU"]["tool_ids"], ["tool"])

    def test_failed_refresh_retains_previous_server_and_retries(self) -> None:
        old = {"url": "https://us.example.org", "tool_ids": ["old"], "repos": ["iuc/old"]}
        write_yaml(self.output, {"fetched_at": "2000-01-01T00:00:00+00:00", "servers": {"US": old}})

        def fetch(url):
            return ({"new"}, {"iuc/new"}) if url == "https://eu.example.org" else None

        with patch.object(availability, "_fetch_installed", side_effect=fetch):
            result = availability.extract_tool_availability()

        self.assertEqual(result["servers"]["US"], old)
        self.assertEqual(result["servers"]["EU"]["tool_ids"], ["new"])
        self.assertEqual(read_yaml(self.output), result)
        self.assertFalse(availability._is_cache_fresh(self.output))

        with patch.object(availability, "_fetch_installed", return_value=({"recovered"}, set())) as fetch:
            recovered = availability.extract_tool_availability()

        self.assertEqual(fetch.call_count, 2)
        self.assertEqual(recovered["servers"]["US"]["tool_ids"], ["recovered"])
        self.assertTrue(availability._is_cache_fresh(self.output))

    def test_first_build_with_unreachable_servers_does_not_cache_empty_results(self) -> None:
        with patch.object(availability, "_fetch_installed", return_value=None):
            result = availability.extract_tool_availability()

        self.assertEqual(result["servers"], {})
        self.assertTrue(self.output.exists())
        self.assertFalse(availability._is_cache_fresh(self.output))
