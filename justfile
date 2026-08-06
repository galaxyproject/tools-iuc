# Repo-root justfile — proxies into website/justfile for convenience.
# Requires `just` (https://github.com/casey/just) and `uv` (https://github.com/astral-sh/uv).

# Show all available targets with descriptions
default:
    @just --list

# One-time setup: install Python deps (incl. Galaxy from git) + Node deps
website-setup:
    just -f website/justfile setup

# Start the Astro dev server with hot reload (extracts tool data first)
website-dev:
    just -f website/justfile dev

# Full build: extract all data + build static site + Pagefind search index
website-build:
    just -f website/justfile build

# Extract all data into website/data/ (tools, ToolShed, stats, people, recent)
website-extract:
    just -f website/justfile extract

# Extract only tool metadata (fastest, no network except ToolShed)
website-extract-tools:
    just -f website/justfile extract-tools

# Full build then preview the site locally
website-serve:
    just -f website/justfile serve

# Lint (ruff) + typecheck (astro check)
website-check:
    just -f website/justfile check

# Run website unit tests
website-test:
    just -f website/justfile test

# Remove all generated data and build artifacts
website-clean:
    just -f website/justfile clean
