from __future__ import annotations

import json
import tempfile
import unittest
from datetime import UTC, datetime
from pathlib import Path

from pipeline.extract_stats import _build_stats


class RecentActivityTests(unittest.TestCase):
    def build(self, events: list[dict]) -> dict:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "events.json"
            path.write_text(json.dumps(events))
            return _build_stats(path, {}, as_of=datetime(2026, 10, 15, 12, tzinfo=UTC))

    def test_last_30_days_counts_distinct_tools_with_exact_boundaries(self) -> None:
        events = [
            ("new", "added", "2026-09-16T00:00:00Z"),
            ("new", "added", "2026-09-17T00:00:00Z"),
            ("new", "boundary", "2026-09-15T14:00:00+02:00"),
            ("new", "too-old", "2026-09-15T11:59:59Z"),
            ("update", "existing", "2026-09-30T00:00:00Z"),
            ("update", "existing", "2026-10-01T00:00:00Z"),
            ("update", "existing", "2026-10-15T12:00:00Z"),
            ("update", "added", "2026-10-14T00:00:00Z"),
            ("update", "future", "2026-10-15T12:00:01Z"),
        ]
        stats = self.build([{"event_type": kind, "tool_id": tool, "commit_date": date} for kind, tool, date in events])
        self.assertEqual(
            stats["last_30_days"],
            {
                "start": "2026-09-15T12:00:00+00:00",
                "end": "2026-10-15T12:00:00+00:00",
                "added": 2,
                "updated": 2,
            },
        )

    def test_empty_history_reports_zero_recent_activity(self) -> None:
        period = self.build([])["last_30_days"]
        self.assertEqual((period["added"], period["updated"]), (0, 0))
