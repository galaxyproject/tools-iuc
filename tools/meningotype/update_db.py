#!/usr/bin/env python3

import pathlib
import shutil
import subprocess

import meningotype

if __name__ == '__main__':
    db_path = 'db'
    original_db_path = pathlib.Path(meningotype.__file__).parent / 'db'
    shutil.copytree(original_db_path, db_path)
    cmd = ["meningotype", "--updatedb", "--db", db_path]
    subprocess.run(cmd)
