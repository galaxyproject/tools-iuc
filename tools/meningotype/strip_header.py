#!/usr/bin/env python3

import sys

if __name__ == '__main__':
    next(sys.stdin)  # skip first line
    for line in sys.stdin:
        sys.stdout.write(line)
