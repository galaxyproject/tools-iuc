import bisect
import math
import re
import subprocess
import sys
import tempfile

import numpy
from six import Iterator

GFF_EXT = 'gff'
SCIDX_EXT = 'scidx'


def noop(data):
    return data


def zeropad_to_numeric(data):
    return re.sub(r'chr0(\d)', r'chr\1', data)


def numeric_to_zeropad(data):
    return re.sub(r'chr(\d([^\d]|$))', r'chr0\1', data)


FORMATS = ['zeropad', 'numeric']
IN_CONVERT = {'zeropad': zeropad_to_numeric, 'numeric': noop}
OUT_CONVERT = {'zeropad': numeric_to_zeropad, 'numeric': noop}


def conversion_functions(in_fmt, out_fmt):
    """
    Returns the proper list of functions to apply to perform a conversion
    """
    return [IN_CONVERT[in_fmt], OUT_CONVERT[out_fmt]]


def convert_data(data, in_fmt, out_fmt):
    for fn in conversion_functions(in_fmt, out_fmt):
        data = fn(data)
    return data


class ChromosomeManager(Iterator):
    """
    Manages a CSV reader of an index file to only load one chrom at a time
    """

    def __init__(self, reader):
        self.done = False
        self.reader = reader
        self.processed_chromosomes = []
        self.current_index = 0
        self.next_valid()

    def __next__(self):
        self.line = next(self.reader)

    def is_valid(self, line):
        if len(line) not in [4, 5, 9]:
            return False
        try:
            [int(i) for i in line[1:]]
            self.format = SCIDX_EXT
            return True
        except ValueError:
            try:
                if len(line) < 6:
                    return False
                [int(line[4]), int(line[5])]
                self.format = GFF_EXT
                return True
            except ValueError:
                return False

    def next_valid(self):
        """
        Advance to the next valid line in the reader
        """
        self.line = next(self.reader)
        s = 0
        while not self.is_valid(self.line):
            self.line = next(self.reader)
            s += 1
        if s > 0:
            # Skip initial line(s) of file
            pass

    def parse_line(self, line):
        if self.format == SCIDX_EXT:
            return [int(line[1]), int(line[2]), int(line[3])]
        else:
            return [int(line[3]), line[6], line[5]]

    def chromosome_name(self):
        """
        Return the name of the chromosome about to be loaded
        """
        return self.line[0]

    def load_chromosome(self, collect_data=True):
        """
        Load the current chromosome into an array and return it
        """
        cname = self.chromosome_name()
        if cname in self.processed_chromosomes:
            stop_err('File is not grouped by chromosome')
        self.data = []
        while self.line[0] == cname:
            if collect_data:
                read = self.parse_line(self.line)
                if read[0] < self.current_index:
                    msg = 'Reads in chromosome %s are not sorted by index. (At index %d)' % (cname, self.current_index)
                    stop_err(msg)
                self.current_index = read[0]
                self.add_read(read)
            try:
                next(self)
            except StopIteration:
                self.done = True
                break
        self.processed_chromosomes.append(cname)
        self.current_index = 0
        data = self.data
        # Don't retain reference anymore to save memory
        del self.data
        return data

    def add_read(self, read):
        if self.format == SCIDX_EXT:
            self.data.append(read)
        else:
            index, strand, value = read
            if value == '' or value == '.':
                value = 1
            else:
                value = int(value)
            if not self.data:
                self.data.append([index, 0, 0])
                current_read = self.data[-1]
            if self.data[-1][0] == index:
                current_read = self.data[-1]
            elif self.data[-1][0] < index:
                self.data.append([index, 0, 0])
                current_read = self.data[-1]
            else:
                msg = 'Reads in chromosome %s are not sorted by index. (At index %d)' % (self.chromosome_name(), index)
                stop_err(msg)
            if strand == '+':
                current_read[1] += value
            elif strand == '-':
                current_read[2] += value
            else:
                msg = 'Strand "%s" at chromosome "%s" index %d is not valid.' % (strand, self.chromosome_name(), index)
                stop_err(msg)

    def skip_chromosome(self):
        """
        Skip the current chromosome, discarding data
        """
        self.load_chromosome(collect_data=False)


class Peak(object):
    def __init__(self, index, pos_width, neg_width):
        self.index = index
        self.start = index - neg_width
        self.end = index + pos_width
        self.value = 0
        self.deleted = False
        self.safe = False

    def __repr__(self):
        return '[%d] %d' % (self.index, self.value)


def gff_row(cname, start, end, score, source, type='.', strand='.', phase='.', attrs={}):
    return (cname, source, type, start, end, score, strand, phase, gff_attrs(attrs))


def gff_attrs(d):
    if not d:
        return '.'
    return ';'.join('%s=%s' % item for item in d.items())


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit(1)


def is_int(i):
    try:
        int(i)
        return True
    except ValueError:
        return False


def make_keys(data):
    return [read[0] for read in data]


def make_peak_keys(peaks):
    return [peak.index for peak in peaks]


def get_window(data, start, end, keys):
    """
    Returns all reads from the data set with index between the two indexes
    """
    start_index = bisect.bisect_left(keys, start)
    end_index = bisect.bisect_right(keys, end)
    return data[start_index:end_index]


def get_index(value, keys):
    """
    Returns the index of the value in the keys using bisect
    """
    return bisect.bisect_left(keys, value)


def get_range(data):
    lo = min([item[0] for item in data])
    hi = max([item[0] for item in data])
    return lo, hi


def get_chunks(lo, hi, size, overlap=500):
    """
    Divides a range into chunks of maximum size size. Returns a list of
    2-tuples (slice_range, process_range), each a 2-tuple (start, end).
    process_range has zero overlap and should be given to process_chromosome
    as-is, and slice_range is overlapped and should be used to slice the
    data (using get_window) to be given to process_chromosome.
    """
    chunks = []
    for start_index in range(lo, hi, size):
        process_start = start_index
        # Don't go over upper bound
        process_end = min(start_index + size, hi)
        # Don't go under lower bound
        slice_start = max(process_start - overlap, lo)
        # Don't go over upper bound
        slice_end = min(process_end + overlap, hi)
        chunks.append(((slice_start, slice_end), (process_start, process_end)))
    return chunks


def allocate_array(data, width):
    """
    Allocates a new array with the dimensions required to fit all reads in
    the argument. The new array is totally empty. Returns the array and the
    shift (number to add to a read index to get the position in the array it
    should be at).
    """
    lo, hi = get_range(data)
    rng = hi - lo
    shift = width - lo
    return numpy.zeros(rng + width * 2, numpy.float), shift


def normal_array(width, sigma, normalize=True):
    """
    Returns an array of the normal distribution of the specified width
    """
    sigma2 = float(sigma)**2

    def normal_func(x):
        return math.exp(-x * x / (2 * sigma2))

    # width is the half of the distribution
    values = list(map(normal_func, range(-width, width)))
    values = numpy.array(values, numpy.float)
    # normalization
    if normalize:
        values = 1.0 / math.sqrt(2 * numpy.pi * sigma2) * values
    return values


def call_peaks(array, shift, data, keys, direction, down_width, up_width, exclusion):
    peaks = []

    def find_peaks():
        # Go through the array and call each peak
        results = (array > numpy.roll(array, 1)) & (array > numpy.roll(array, -1))
        indexes = numpy.where(results)
        for index in indexes[0]:
            pos = down_width or exclusion // 2
            neg = up_width or exclusion // 2
            # Reverse strand
            if direction == 2:
                # Swap positive and negative widths
                pos, neg = neg, pos
            peaks.append(Peak(int(index) - shift, pos, neg))
    find_peaks()

    def calculate_reads():
        # Calculate the number of reads in each peak
        for peak in peaks:
            reads = get_window(data, peak.start, peak.end, keys)
            peak.value = sum([read[direction] for read in reads])
            # Flat list of indexes with frequency
            indexes = [r for read in reads for r in [read[0]] * read[direction]]
            peak.stddev = numpy.std(indexes)
    calculate_reads()

    def perform_exclusion():
        # Process the exclusion zone
        peak_keys = make_peak_keys(peaks)
        peaks_by_value = peaks[:]
        peaks_by_value.sort(key=lambda peak: -peak.value)
        for peak in peaks_by_value:
            peak.safe = True
            window = get_window(peaks,
                                peak.index - exclusion // 2,
                                peak.index + exclusion // 2,
                                peak_keys)
            for excluded in window:
                if excluded.safe:
                    continue
                i = get_index(excluded.index, peak_keys)
                del peak_keys[i]
                del peaks[i]
    perform_exclusion()
    return peaks


def process_chromosome(cname, data, writer, process_bounds, width, sigma, down_width, up_width, exclusion, filter):
    """
    Process a chromosome. Takes the chromosome name, list of reads, a CSV
    writer to write processes results to, the bounds (2-tuple) to write
    results in, and options.
    """
    if not data:
        return
    keys = make_keys(data)
    # Create the arrays that hold the sum of the normals
    forward_array, forward_shift = allocate_array(data, width)
    reverse_array, reverse_shift = allocate_array(data, width)
    normal = normal_array(width, sigma)

    def populate_array():
        # Add each read's normal to the array
        for read in data:
            index, forward, reverse = read
            # Add the normals to the appropriate regions
            if forward:
                forward_array[index + forward_shift - width:index + forward_shift + width] += normal * forward
            if reverse:
                reverse_array[index + reverse_shift - width:index + reverse_shift + width] += normal * reverse
    populate_array()
    forward_peaks = call_peaks(forward_array, forward_shift, data, keys, 1, down_width, up_width, exclusion)
    reverse_peaks = call_peaks(reverse_array, reverse_shift, data, keys, 2, down_width, up_width, exclusion)
    # Convert chromosome name in preparation for writing output
    cname = convert_data(cname, 'zeropad', 'numeric')

    def write(cname, strand, peak):
        start = max(peak.start, 1)
        end = peak.end
        value = peak.value
        stddev = peak.stddev
        if value > filter:
            # This version of genetrack outputs only gff files.
            writer.writerow(gff_row(cname=cname,
                                    source='genetrack',
                                    start=start,
                                    end=end,
                                    score=value,
                                    strand=strand,
                                    attrs={'stddev': stddev}))

    for peak in forward_peaks:
        if process_bounds[0] < peak.index < process_bounds[1]:
            write(cname, '+', peak)
    for peak in reverse_peaks:
        if process_bounds[0] < peak.index < process_bounds[1]:
            write(cname, '-', peak)


def sort_chromosome_reads_by_index(input_path):
    """
    Return a gff file with chromosome reads sorted by index.
    """
    # Will this sort produce different results across platforms?
    output_path = tempfile.NamedTemporaryFile(delete=False).name
    command = 'sort -k 1,1 -k 4,4n "%s" > "%s"' % (input_path, output_path)
    p = subprocess.Popen(command, shell=True)
    p.wait()
    return output_path
