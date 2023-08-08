import argparse

import pandas as pd


def parse_unequal_column_data(in_file, out_file, fillna="NA"):
    # Delimiter
    data_file_delimiter = '\t'

    # The max column count a line in the file could have
    largest_column_count = 0

    # Loop the data lines
    with open(in_file, 'r') as temp_f:
        # Read the lines
        lines = temp_f.readlines()

        for line in lines:
            # Count the column count for the current line
            column_count = len(line.split(data_file_delimiter)) + 1

            # Set the new most column count
            largest_column_count = column_count if \
                largest_column_count < column_count \
                else largest_column_count

    # Generate column names (will be 0, 1, 2, ..., largest_column_count - 1)
    column_names = [i for i in range(0, largest_column_count)]

    # Read csv
    df = pd.read_csv(in_file,
                     header=None,
                     delimiter=data_file_delimiter,
                     names=column_names)

    print(f"Generating a new table with shape: {df.shape}")

    df.dropna(axis=1, how='all', inplace=True)
    df.fillna(fillna, inplace=True)
    df.to_csv(out_file, sep="\t", index=False, header=None)

    print("Done.")


parser = argparse.ArgumentParser()
parser.add_argument('--input', dest='input', required=True)
parser.add_argument('--output', dest='output', required=True)
parser.add_argument('--fillna', dest='fillna', default="NA", required=False)
args = parser.parse_args()

parse_unequal_column_data(args.input, args.output, args.fillna)
