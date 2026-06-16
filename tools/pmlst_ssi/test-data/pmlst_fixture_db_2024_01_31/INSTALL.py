#!/usr/bin/env python3

"""Fixture-only database installer placeholder.

The Galaxy fixture database is intentionally FASTA/BLAST-only and omits KMA
indexes to stay small. The maintained pMLST runtime requires this filename as a
root database marker, but Galaxy tests do not execute it.
"""


def main() -> None:
    print("This fixture database is intended for FASTA/BLAST Galaxy tests only.")


if __name__ == "__main__":
    main()
