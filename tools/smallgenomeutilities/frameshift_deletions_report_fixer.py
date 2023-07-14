"""Polish the output of the frameshift_deletions_check command.

- Drops the first index column, which is rather pointless to include
- Turns ref bases printed as literal bytes strings into plain output
- Removes [] around pos lists and spaces after comma separating list elements
- Turns None and empty list values into . as a cell placeholder
"""

import re
import sys


def matchrepl(matchobj):
    print(matchobj.groups())
    bytes_string_content = matchobj.group(1)
    if bytes_string_content is not None:
        return bytes_string_content
    list_content = matchobj.group(2)
    if list_content is not None:
        if list_content == '':
            return '.'
        return list_content.replace(', ', ',')
    none_cell = matchobj.group(3)
    if none_cell is not None:
        return '\t.\t'

    raise ValueError('Error in regex parsing code')


if __name__ == '__main__':
    regex = re.compile(r"b'(.+)'|\[([^\]]*)\]|\t(None)\t")
    with open(sys.argv[1]) as i:
        with open(sys.argv[2], 'w') as o:
            for line in i:
                line = line[line.index('\t') + 1:]
                o.write(regex.sub(matchrepl, line))
