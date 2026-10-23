# Repo-root justfile — proxies into website/justfile for convenience.
# Requires `just` (https://github.com/casey/just) and `uv` (https://github.com/astral-sh/uv).

# Show available commands
[private]
default:
    @just --list

# Install Python and Node dependencies
website-setup:
    just -f website/justfile setup

# Start the development server
website-dev:
    just -f website/justfile dev

# Extract data and build the searchable site
website-build:
    just -f website/justfile build

# Refresh all website data
website-extract:
    just -f website/justfile extract

# Refresh tool metadata
website-extract-tools:
    just -f website/justfile extract-tools

# Build and preview the site
website-serve:
    just -f website/justfile serve

# Run unit tests, lint, and type checks
website-check:
    just -f website/justfile check

# Run Python and JavaScript unit tests
website-test:
    just -f website/justfile test

# Remove generated data and build files
website-clean:
    just -f website/justfile clean

# Refresh the incremental Git history index
website-history:
    just -f website/justfile history
