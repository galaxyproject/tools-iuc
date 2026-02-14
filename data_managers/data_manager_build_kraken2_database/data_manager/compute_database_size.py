#!/usr/bin/env python3

import math
import pathlib
import sys

if __name__ == "__main__":
    database_path = sys.argv[1]
    all_files = [f for f in pathlib.Path(database_path).rglob("*") if f.is_file()]
    db_size = math.ceil(sum(f.stat().st_size for f in all_files) / (1024 * 1024))
    print(db_size, end="")
