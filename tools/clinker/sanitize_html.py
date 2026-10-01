#!/usr/bin/env python3
"""Escape clinker data before it is embedded in an HTML script element."""

import argparse
import json
from pathlib import Path


MARKER = "const data="


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--html", required=True, type=Path, help="clinker HTML output")
    return parser.parse_args()


def sanitize_html(path):
    """Escape HTML-sensitive characters in clinker\'s embedded data JSON."""
    if not path.is_file():
        raise SystemExit(f"HTML output does not exist: {path}")

    html = path.read_text(encoding="utf-8")
    if html.count(MARKER) != 1:
        raise SystemExit("Expected exactly one embedded clinker data object")

    start = html.index(MARKER) + len(MARKER)
    decoder = json.JSONDecoder()
    try:
        data, end = decoder.raw_decode(html[start:])
    except json.JSONDecodeError as error:
        raise SystemExit(f"Could not parse embedded clinker data: {error}") from error

    if start + end >= len(html) or html[start + end] != ";":
        raise SystemExit("Embedded clinker data is not followed by a semicolon")

    serialized = json.dumps(data, ensure_ascii=True, separators=(",", ":"))
    serialized = serialized.replace("&", r"\u0026")
    serialized = serialized.replace("<", r"\u003c")
    serialized = serialized.replace(">", r"\u003e")
    rewritten = html[:start] + serialized + html[start + end:]

    temporary = path.with_name(f".{path.name}.sanitized")
    try:
        temporary.write_text(rewritten, encoding="utf-8")
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def main():
    """Sanitize the requested clinker HTML output."""
    sanitize_html(parse_args().html)


if __name__ == "__main__":
    main()
