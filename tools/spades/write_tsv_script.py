#!/usr/bin/env python

import re
import sys

search_str = r"^>(NODE|\S+)_(\d+)(?:_|\s)length_(\d+)_cov_(\d+\.*\d*)(.*\$)?"

replace_str = r"\1_\2\t\3\t\4"

cmd = re.compile(search_str)

sys.stdout.write("#name\tlength\tcoverage\n")

for i, line in enumerate(sys.stdin):
    if cmd.match(line):
        sys.stdout.write(cmd.sub(replace_str, line))
