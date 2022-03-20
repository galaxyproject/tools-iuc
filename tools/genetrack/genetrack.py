"""
genetrack.py

Input: either scidx or gff format of reads
Output: Called peaks in gff format
"""
import csv
import optparse
import os

import genetrack_util

CHUNK_SIZE = 10000000


if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('-t', '--input_format', dest='input_format', type='string', help='Input format')
    parser.add_option('-i', '--input', dest='inputs', type='string', action='append', nargs=2, help='Input datasets')
    parser.add_option('-s', '--sigma', dest='sigma', type='int', default=5, help='Sigma.')
    parser.add_option('-e', '--exclusion', dest='exclusion', type='int', default=20, help='Exclusion zone.')
    parser.add_option('-u', '--up_width', dest='up_width', type='int', default=10, help='Upstream width of called peaks.')
    parser.add_option('-d', '--down_width', dest='down_width', type='int', default=10, help='Downstream width of called peaks.')
    parser.add_option('-f', '--filter', dest='filter', type='int', default=1, help='Absolute read filter.')
    options, args = parser.parse_args()

    os.mkdir('output')
    for (dataset_path, hid) in options.inputs:
        if options.input_format == 'gff':
            # Make sure the reads for each chromosome are sorted by index.
            input_path = genetrack_util.sort_chromosome_reads_by_index(dataset_path)
        else:
            # We're processing scidx data.
            input_path = dataset_path
        output_name = 's%se%su%sd%sF%s_on_data_%s' % (options.sigma,
                                                      options.exclusion,
                                                      options.up_width,
                                                      options.down_width,
                                                      options.filter,
                                                      hid)
        output_path = os.path.join('output', output_name)
        reader = csv.reader(open(input_path, 'rU'), delimiter='\t')
        writer = csv.writer(open(output_path, 'wt'), delimiter='\t')
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
                genetrack_util.process_chromosome(cname,
                                                  window,
                                                  writer,
                                                  process_bounds,
                                                  width,
                                                  options.sigma,
                                                  options.up_width,
                                                  options.down_width,
                                                  options.exclusion,
                                                  options.filter)
