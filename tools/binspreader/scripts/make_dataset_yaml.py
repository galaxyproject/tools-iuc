#!/usr/bin/env python3
"""Write a SPAdes style dataset description for BinSPreader."""

import argparse
import sys


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, help='Path of the dataset description to write')
    parser.add_argument('--library-type', required=True, help='Type of the read library')
    parser.add_argument('--orientation', required=True, help='Orientation of the read library')
    parser.add_argument('--left', nargs='+', required=True, help='Forward reads')
    parser.add_argument('--right', nargs='+', required=True, help='Reverse reads')
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    left = ', '.join('"%s"' % path for path in args.left)
    right = ', '.join('"%s"' % path for path in args.right)
    with open(args.output, 'w') as handle:
        handle.write('[\n')
        handle.write('  {\n')
        handle.write('    orientation: "%s",\n' % args.orientation)
        handle.write('    type: "%s",\n' % args.library_type)
        handle.write('    left reads: [%s],\n' % left)
        handle.write('    right reads: [%s]\n' % right)
        handle.write('  }\n')
        handle.write(']\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
