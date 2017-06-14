"""
genetrack.py

Input: Any combination of scidx and gff format of reads
Output: Called peaks in gff format
"""
import csv
import optparse
import os
import sys

sys.dont_write_bytecode = True

import genetrack_util

CHUNK_SIZE = 10000000


if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('-g', '--gff', dest='gff_inputs', type='string', action='append', nargs=1, help='Input GFF datasets')
    parser.add_option('-x', '--scidx', dest='scidx_inputs', type='string', action='append', nargs=1, help='Input ScIdx datasets')
    parser.add_option('-s', '--sigma', dest='sigma', type='int', default=5, help='Sigma.')
    parser.add_option('-e', '--exclusion', dest='exclusion', type='int', default=20, help='Exclusion zone.')
    parser.add_option('-u', '--up_width', dest='up_width', type='int', default=10, help='Upstream width of called peaks.')
    parser.add_option('-d', '--down_width', dest='down_width', type='int', default=10, help='Downstream width of called peaks.')
    parser.add_option('-f', '--filter', dest='filter', type='int', default=1, help='Absolute read filter.')
    parser.add_option('-o', '--output', dest='output', type='string', help='Output filename.')
    options, args = parser.parse_args()

    input_files = []
    # print(options)
    if options.gff_inputs is not None:
        input_files.extend([genetrack_util.sort_chromosome_reads_by_index(input_file) for input_file in options.gff_inputs])
    if options.scidx_inputs is not None:
        input_files.extend(options.scidx_inputs)

    with open(options.output, 'wt') as output_fh:
        for dataset_path in input_files:
            reader = csv.reader(open(dataset_path, 'rU'), delimiter='\t')
            writer = csv.writer(output_fh, delimiter='\t')
            width = options.sigma * 5
            manager = genetrack_util.ChromosomeManager(reader)
            while not manager.done:
                cname = manager.chromosome_name()
                # Should we process this chromosome?
                data = manager.load_chromosome()
                if not data:
                    continue
                keys = genetrack_util.make_keys(data)
                lo, hi = genetrack_util.get_range(data)
                for chunk in genetrack_util.get_chunks(lo, hi, size=CHUNK_SIZE, overlap=width):
                    (slice_start, slice_end), process_bounds = chunk
                    window = genetrack_util.get_window(data, slice_start, slice_end, keys)
                    genetrack_util.process(cname,
                                           window,
                                           writer,
                                           process_bounds,
                                           width,
                                           options.sigma,
                                           options.up_width,
                                           options.down_width,
                                           options.exclusion,
                                           options.filter)
