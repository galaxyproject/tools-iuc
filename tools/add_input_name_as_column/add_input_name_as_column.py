import io
import argparse


def Parser():
    the_parser = argparse.ArgumentParser(description="add label to last column of file")
    the_parser.add_argument('--input', required=True, action="store", type=str, help="input tabular file")
    the_parser.add_argument('--output', required=True, action="store", type=str, help="output file path")
    the_parser.add_argument('--label', required=True, action="store", type=str, help="label to add in last column")
    the_parser.add_argument('--header', action="store", type=str, help="column label for last column")
    args = the_parser.parse_args()
    return args


args = Parser()


with io.open(args.input, encoding="utf-8") as input, io.open(args.output, 'w', encoding="utf-8") as output:
    for i, line in enumerate(input):
        line = line.strip('\n')
        if (i == 0) and args.header:
            line = "%s\t%s\n" % (line, args.header)
        else:
            line = "%s\t%s\n" % (line, args.label)
        output.write(line)
