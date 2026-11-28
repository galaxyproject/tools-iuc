#!/usr/bin/env bash
# The same validation entry point is used locally and in CI.
set -euo pipefail
cd "$(dirname "$0")"
uv run python -m unittest discover -s tests
uv run ruff check pipeline
cd astro
npm test
npm run check
