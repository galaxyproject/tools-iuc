import csv
import json
from datetime import datetime
from pathlib import Path
from typing import Any

from .data_model import ToolEvent


def _serialize(val: Any) -> str:
    if isinstance(val, datetime):
        return val.isoformat()
    if val is None:
        return ''
    return str(val)


def write_events_csv(events: list[ToolEvent], path: str):
    fields = [
        'tool_id', 'tool_name', 'directory', 'event_type',
        'version', 'old_version', 'change_type', 'commit_hash',
        'commit_date', 'commit_message',
    ]
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(fields)
        for e in events:
            w.writerow([
                e.tool_id, e.tool_name, e.directory, e.event_type,
                e.version, e.old_version or '', e.change_type or '',
                e.commit_hash,
                e.commit_date.isoformat() if e.commit_date else '',
                e.commit_message,
            ])


def write_events_json(events: list[ToolEvent], path: str):
    data = []
    for e in events:
        data.append({
            'tool_id': e.tool_id,
            'tool_name': e.tool_name,
            'directory': e.directory,
            'event_type': e.event_type,
            'version': e.version,
            'old_version': e.old_version,
            'change_type': e.change_type,
            'commit_hash': e.commit_hash,
            'commit_date': e.commit_date.isoformat()
            if e.commit_date else None,
            'commit_message': e.commit_message,
        })

    with open(path, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2)


def write_summary(events: list[ToolEvent], path: str):
    additions = [e for e in events if e.event_type == 'new']
    updates = [e for e in events if e.event_type == 'update']
    unique_tools = len({e.tool_id for e in events})

    categories: dict[str, int] = {}
    for e in additions:
        yr = e.commit_date.year if e.commit_date else 'unknown'
        categories[yr] = categories.get(yr, 0) + 1

    lines = [
        f'Total deployment events:  {len(events)}',
        f'New tool additions:        {len(additions)}',
        f'Tool version updates:      {len(updates)}',
        f'Unique tools tracked:      {unique_tools}',
        '',
        'New tools by year:',
    ]
    for yr in sorted(categories):
        lines.append(f'  {yr}: {categories[yr]}')

    update_counts: dict[str, int] = {}
    for e in updates:
        yr = e.commit_date.year if e.commit_date else 'unknown'
        update_counts[yr] = update_counts.get(yr, 0) + 1

    lines.append('')
    lines.append('Updates by year:')
    for yr in sorted(update_counts):
        lines.append(f'  {yr}: {update_counts[yr]}')

    with open(path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')
