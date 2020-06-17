import bisect
import csv
import os
import sys
import traceback

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot  # noqa: I202,E402

# Data outputs
DETAILS = 'D'
MATCHED_PAIRS = 'MP'
ORPHANS = 'O'
# Data output formats
GFF_EXT = 'gff'
TABULAR_EXT = 'tabular'
# Statistics histograms output directory.
HISTOGRAM = 'H'
# Statistics outputs
FINAL_PLOTS = 'F'
PREVIEW_PLOTS = 'P'
STATS_GRAPH = 'C'

# Graph settings.
COLORS = 'krg'
Y_LABEL = 'Peak-pair counts'
X_LABEL = 'Peak-pair distance (bp)'
TICK_WIDTH = 3
ADJUST = [0.140, 0.9, 0.9, 0.1]
PLOT_FORMAT = 'pdf'
pyplot.rc('xtick.major', size=10.00)
pyplot.rc('ytick.major', size=10.00)
pyplot.rc('lines', linewidth=4.00)
pyplot.rc('axes', linewidth=3.00)
pyplot.rc('font', family='Bitstream Vera Sans', size=32.0)


class FrequencyDistribution(object):

    def __init__(self, start, end, binsize=10, d=None):
        self.start = start
        self.end = end
        self.dist = d or {}
        self.binsize = binsize

    def get_bin(self, x):
        """
        Returns the centre of the bin in which a data point falls
        """
        return self.start + (x - self.start) // self.binsize * self.binsize + self.binsize / 2.0

    def add(self, x):
        x = self.get_bin(x)
        self.dist[x] = self.dist.get(x, 0) + 1

    def graph_series(self):
        x = []
        y = []
        for i in range(self.start, self.end, self.binsize):
            center = self.get_bin(i)
            x.append(center)
            y.append(self.dist.get(center, 0))
        return x, y

    def mode(self):
        # There could be more than one mode for a frequency distribution,
        # return the median of the modes to be consistent
        max_frequency = max(self.dist.values())
        modes = sorted(_[0] for _ in self.dist.items() if _[1] == max_frequency)
        median_index = len(modes) // 2
        return modes[median_index]

    def size(self):
        return sum(self.dist.values())


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit(1)


def distance(peak1, peak2):
    return (peak2[1] + peak2[2]) / 2.0 - (peak1[1] + peak1[2]) / 2.0


def gff_row(cname, start, end, score, source, type='.', strand='.', phase='.', attrs={}):
    return (cname, source, type, start, end, score, strand, phase, gff_attrs(attrs))


def gff_attrs(d):
    if not d:
        return '.'
    return ';'.join('%s=%s' % item for item in d.items())


def parse_chromosomes(reader):
    # This version of cwpair2 accepts only gff format as input.
    chromosomes = {}
    for line in reader:
        line = line.rstrip("\r\n")
        if not line or line.startswith('#'):
            continue
        cname, _, _, start, end, value, strand, _, _ = line.split("\t")
        start = int(start)
        end = int(end)
        value = float(value)
        if cname not in chromosomes:
            chromosomes[cname] = []
        peaks = chromosomes[cname]
        peaks.append((strand, start, end, value))
    return chromosomes


def perc95(chromosomes):
    """
    Returns the 95th percentile value of the given chromosomes.
    """
    values = []
    for peaks in chromosomes.values():
        for peak in peaks:
            values.append(peak[3])
    values.sort()
    # Get 95% value
    return values[int(len(values) * 0.95)]


def peak_filter(chromosomes, threshold):
    """
    Filters the peaks to those above a threshold. Threshold < 1.0 is interpreted
    as a proportion of the maximum, >=1.0 as an absolute value.
    """
    if threshold < 1:
        p95 = perc95(chromosomes)
        threshold = p95 * threshold
        # Make the threshold a proportion of the
    for cname, peaks in chromosomes.items():
        chromosomes[cname] = [peak for peak in peaks if peak[3] > threshold]


def split_strands(chromosome):
    watson = [peak for peak in chromosome if peak[0] == '+']
    crick = [peak for peak in chromosome if peak[0] == '-']
    return watson, crick


def all_pair_distribution(chromosomes, up_distance, down_distance, binsize):
    dist = FrequencyDistribution(-up_distance, down_distance, binsize=binsize)
    for data in chromosomes.values():
        watson, crick = split_strands(data)
        crick.sort(key=lambda data: float(data[1]))
        keys = make_keys(crick)
        for peak in watson:
            for cpeak in get_window(crick, peak, up_distance, down_distance, keys):
                dist.add(distance(peak, cpeak))
    return dist


def make_keys(crick):
    return [(data[1] + data[2]) // 2 for data in crick]


def get_window(crick, peak, up_distance, down_distance, keys=None):
    """
    Returns a window of all crick peaks within a distance of a watson peak.
    crick strand MUST be sorted by distance
    """
    strand, start, end, value = peak
    midpoint = (start + end) // 2
    lower = midpoint - up_distance
    upper = midpoint + down_distance
    keys = keys or make_keys(crick)
    start_index = bisect.bisect_left(keys, lower)
    end_index = bisect.bisect_right(keys, upper)
    return [cpeak for cpeak in crick[start_index:end_index]]


def match_largest(window, peak):
    if not window:
        return None
    return max(window, key=lambda cpeak: cpeak[3])


def match_closest(window, peak):
    if not window:
        return None

    def key(cpeak):
        d = distance(peak, cpeak)
        # Search negative distances last
        if d < 0:
            # And then prefer less negative distances
            d = 10000 - d
        return d
    return min(window, key=key)


def match_mode(window, peak, mode):
    if not window:
        return None
    return min(window, key=lambda cpeak: abs(distance(peak, cpeak) - mode))


METHODS = {'mode': match_mode, 'closest': match_closest, 'largest': match_largest}


def frequency_plot(freqs, fname, labels=[], title=''):
    pyplot.clf()
    pyplot.figure(figsize=(10, 10))
    for i, freq in enumerate(freqs):
        x, y = freq.graph_series()
        pyplot.plot(x, y, '%s-' % COLORS[i])
    if len(freqs) > 1:
        pyplot.legend(labels)
    pyplot.xlim(freq.start, freq.end)
    pyplot.ylim(ymin=0)
    pyplot.ylabel(Y_LABEL)
    pyplot.xlabel(X_LABEL)
    pyplot.subplots_adjust(left=ADJUST[0], right=ADJUST[1], top=ADJUST[2], bottom=ADJUST[3])
    # Get the current axes
    ax = pyplot.gca()
    for l in ax.get_xticklines() + ax.get_yticklines():
        l.set_markeredgewidth(TICK_WIDTH)
    pyplot.savefig(fname)


def create_directories():
    # Output histograms in pdf.
    os.mkdir(HISTOGRAM)
    os.mkdir('data_%s' % DETAILS)
    os.mkdir('data_%s' % ORPHANS)
    os.mkdir('data_%s' % MATCHED_PAIRS)


def process_file(dataset_path, galaxy_hid, method, threshold, up_distance,
                 down_distance, binsize, output_files):
    if method == 'all':
        match_methods = METHODS.keys()
    else:
        match_methods = [method]
    statistics = []
    for match_method in match_methods:
        stats = perform_process(dataset_path,
                                galaxy_hid,
                                match_method,
                                threshold,
                                up_distance,
                                down_distance,
                                binsize,
                                output_files)
        statistics.append(stats)
    if output_files == 'all' and method == 'all':
        frequency_plot([s['dist'] for s in statistics],
                       statistics[0]['graph_path'],
                       labels=list(METHODS.keys()))
    return statistics


def perform_process(dataset_path, galaxy_hid, method, threshold, up_distance,
                    down_distance, binsize, output_files):
    output_details = output_files in ["all", "matched_pair_orphan_detail"]
    output_plots = output_files in ["all"]
    output_orphans = output_files in ["all", "matched_pair_orphan", "matched_pair_orphan_detail"]
    # Keep track of statistics for the output file
    statistics = {}
    fpath, fname = os.path.split(dataset_path)
    statistics['fname'] = '%s: data %s' % (method, str(galaxy_hid))
    statistics['dir'] = fpath
    if threshold >= 1:
        filter_string = 'fa%d' % threshold
    else:
        filter_string = 'f%d' % (threshold * 100)
    fname = '%s_%su%dd%d_on_data_%s' % (method, filter_string, up_distance, down_distance, galaxy_hid)

    def make_histogram_path(output_type, fname):
        return os.path.join(HISTOGRAM, 'histogram_%s_%s.%s' % (output_type, fname, PLOT_FORMAT))

    def make_path(output_type, extension, fname):
        # Returns the full path for an output.
        return os.path.join(output_type, '%s_%s.%s' % (output_type, fname, extension))

    def td_writer(output_type, extension, fname):
        # Returns a tab-delimited writer for a specified output.
        output_file_path = make_path(output_type, extension, fname)
        return csv.writer(open(output_file_path, 'wt'), delimiter='\t', lineterminator="\n")

    with open(dataset_path, 'rt') as input:
        try:
            chromosomes = parse_chromosomes(input)
        except Exception:
            stop_err('Unable to parse file "%s".\n%s' % (dataset_path, traceback.format_exc()))
    if output_details:
        # Details
        detailed_output = td_writer('data_%s' % DETAILS, TABULAR_EXT, fname)
        detailed_output.writerow(('chrom', 'start', 'end', 'value', 'strand') * 2 + ('midpoint', 'c-w reads sum', 'c-w distance (bp)'))
    if output_plots:
        # Final Plot
        final_plot_path = make_histogram_path(FINAL_PLOTS, fname)
    if output_orphans:
        # Orphans
        orphan_output = td_writer('data_%s' % ORPHANS, TABULAR_EXT, fname)
        orphan_output.writerow(('chrom', 'strand', 'start', 'end', 'value'))
    if output_plots:
        # Preview Plot
        preview_plot_path = make_histogram_path(PREVIEW_PLOTS, fname)
    # Matched Pairs.
    matched_pairs_output = td_writer('data_%s' % MATCHED_PAIRS, GFF_EXT, fname)
    statistics['stats_path'] = 'statistics.%s' % TABULAR_EXT
    if output_plots:
        statistics['graph_path'] = make_histogram_path(STATS_GRAPH, fname)
    statistics['perc95'] = perc95(chromosomes)
    if threshold > 0:
        # Apply peak_filter
        peak_filter(chromosomes, threshold)
    if method == 'mode':
        freq = all_pair_distribution(chromosomes, up_distance, down_distance, binsize)
        mode = freq.mode()
        statistics['preview_mode'] = mode
        if output_plots:
            frequency_plot([freq], preview_plot_path, title='Preview frequency plot')
    else:
        statistics['preview_mode'] = 'NA'
    dist = FrequencyDistribution(-up_distance, down_distance, binsize=binsize)
    orphans = 0
    # x will be used to archive the summary dataset
    x = []
    for cname, chromosome in chromosomes.items():
        # Each peak is (strand, start, end, value)
        watson, crick = split_strands(chromosome)
        # Sort by value of each peak
        watson.sort(key=lambda data: -float(data[3]))
        # Sort by position to facilitate binary search
        crick.sort(key=lambda data: float(data[1]))
        keys = make_keys(crick)
        for peak in watson:
            window = get_window(crick, peak, up_distance, down_distance, keys)
            if method == 'mode':
                match = match_mode(window, peak, mode)
            else:
                match = METHODS[method](window, peak)
            if match:
                midpoint = (match[1] + match[2] + peak[1] + peak[2]) // 4
                d = distance(peak, match)
                dist.add(d)
                # Simple output in gff format.
                x.append(gff_row(cname,
                                 source='cwpair',
                                 start=midpoint,
                                 end=midpoint + 1,
                                 score=peak[3] + match[3],
                                 attrs={'cw_distance': d}))
                if output_details:
                    detailed_output.writerow((cname,
                                              peak[1],
                                              peak[2],
                                              peak[3],
                                              '+',
                                              cname,
                                              match[1],
                                              match[2],
                                              match[3], '-',
                                              midpoint,
                                              peak[3] + match[3],
                                              d))
                i = bisect.bisect_left(keys, (match[1] + match[2]) / 2)
                del crick[i]
                del keys[i]
            else:
                if output_orphans:
                    orphan_output.writerow((cname, peak[0], peak[1], peak[2], peak[3]))
                # Keep track of orphans for statistics.
                orphans += 1
        # Remaining crick peaks are orphans
        if output_orphans:
            for cpeak in crick:
                orphan_output.writerow((cname, cpeak[0], cpeak[1], cpeak[2], cpeak[3]))
        # Keep track of orphans for statistics.
        orphans += len(crick)
    # Sort output descending by score.
    x.sort(key=lambda data: float(data[5]), reverse=True)
    # Writing a summary to gff format file
    for row in x:
        row_tmp = list(row)
        # Dataset in tuple cannot be modified in Python, so row will
        # be converted to list format to add 'chr'.
        if row_tmp[0] == "999":
            row_tmp[0] = 'chrM'
        elif row_tmp[0] == "998":
            row_tmp[0] = 'chrY'
        elif row_tmp[0] == "997":
            row_tmp[0] = 'chrX'
        else:
            row_tmp[0] = row_tmp[0]
        # Print row_tmp.
        matched_pairs_output.writerow(row_tmp)
    statistics['paired'] = dist.size() * 2
    statistics['orphans'] = orphans
    statistics['final_mode'] = dist.mode()
    if output_plots:
        frequency_plot([dist], final_plot_path, title='Frequency distribution')
    statistics['dist'] = dist
    return statistics
