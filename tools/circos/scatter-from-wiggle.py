#!/usr/bin/env python
import wiggle
import logging
import sys

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

if __name__ == '__main__':
    walker = wiggle.Wiggle()
    for line in walker.walk(sys.stdin):
        sys.stdout.write('%s\t%s\t%s\t%s\n' % line)
