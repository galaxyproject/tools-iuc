import bisect
import csv
import os
import shutil
import sys
import tempfile

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot  # noqa: I202,E402

# Graph settings
Y_LABEL = 'Counts'
X_LABEL = 'Number of matched replicates'
TICK_WIDTH = 3
# Amount to shift the graph to make labels fit, [left, right, top, bottom]
ADJUST = [0.180, 0.9, 0.9, 0.1]
# Length of tick marks, use TICK_WIDTH for width
pyplot.rc('xtick.major', size=10.00)
pyplot.rc('ytick.major', size=10.00)
pyplot.rc('lines', linewidth=4.00)
pyplot.rc('axes', linewidth=3.00)
pyplot.rc('font', family='Bitstream Vera Sans', size=32.0)

COLORS = 'krb'
ISPY2 = sys.version_info[0] == 2


class Replicate(object):

    def __init__(self, id, dataset_path):
        self.id = id
        self.dataset_path = dataset_path
        if ISPY2:
            fh = open(dataset_path, 'rb')
        else:
            fh = open(dataset_path, 'r', newline='')
        self.parse(csv.reader(fh, delimiter='\t'))

    def parse(self, reader):
        self.chromosomes = {}
        for line in reader:
            if line[0].startswith("#") or line[0].startswith('"'):
                continue
            cname, junk, junk, mid, midplus, value, strand, junk, attrs = line
            attrs = parse_gff_attrs(attrs)
            distance = int(attrs['cw_distance'])
            mid = int(mid)
            midplus = int(midplus)
            value = float(value)
            if cname not in self.chromosomes:
                self.chromosomes[cname] = Chromosome(cname)
            chrom = self.chromosomes[cname]
            chrom.add_peak(Peak(cname, mid, value, distance, self))
        for chrom in self.chromosomes.values():
            chrom.sort_by_index()

    def filter(self, up_limit, low_limit):
        for chrom in self.chromosomes.values():
            chrom.filter(up_limit, low_limit)

    def size(self):
        return sum([len(c.peaks) for c in self.chromosomes.values()])


class Chromosome(object):

    def __init__(self, name):
        self.name = name
        self.peaks = []

    def add_peak(self, peak):
        self.peaks.append(peak)

    def sort_by_index(self):
        self.peaks.sort(key=lambda peak: peak.midpoint)
        self.keys = make_keys(self.peaks)

    def remove_peak(self, peak):
        i = bisect.bisect_left(self.keys, peak.midpoint)
        # If the peak was actually found
        if i < len(self.peaks) and self.peaks[i].midpoint == peak.midpoint:
            del self.keys[i]
            del self.peaks[i]

    def filter(self, up_limit, low_limit):
        self.peaks = [p for p in self.peaks if low_limit <= p.distance <= up_limit]
        self.keys = make_keys(self.peaks)


class Peak(object):

    def __init__(self, chrom, midpoint, value, distance, replicate):
        self.chrom = chrom
        self.value = value
        self.midpoint = midpoint
        self.distance = distance
        self.replicate = replicate

    def normalized_value(self, med):
        return self.value * med / self.replicate.median


class PeakGroup(object):

    def __init__(self):
        self.peaks = {}

    def add_peak(self, repid, peak):
        self.peaks[repid] = peak

    @property
    def chrom(self):
        return list(self.peaks.values())[0].chrom

    @property
    def midpoint(self):
        return int(median([peak.midpoint for peak in self.peaks.values()]))

    @property
    def num_replicates(self):
        return len(self.peaks)

    @property
    def median_distance(self):
        return int(median([peak.distance for peak in self.peaks.values()]))

    @property
    def value_sum(self):
        return sum([peak.value for peak in self.peaks.values()])

    def normalized_value(self, med):
        values = []
        for peak in self.peaks.values():
            values.append(peak.normalized_value(med))
        return median(values)

    @property
    def peakpeak_distance(self):
        keys = list(self.peaks.keys())
        return abs(self.peaks[keys[0]].midpoint - self.peaks[keys[1]].midpoint)


class FrequencyDistribution(object):

    def __init__(self, d=None):
        self.dist = d or {}

    def add(self, x):
        self.dist[x] = self.dist.get(x, 0) + 1

    def graph_series(self):
        x = []
        y = []
        for key, val in self.dist.items():
            x.append(key)
            y.append(val)
        return x, y

    def mode(self):
        return max(self.dist.items(), key=lambda data: data[1])[0]

    def size(self):
        return sum(self.dist.values())


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit(1)


def median(data):
    """
    Find the integer median of the data set.
    """
    if not data:
        return 0
    sdata = sorted(data)
    if len(data) % 2 == 0:
        return (sdata[len(data) // 2] + sdata[len(data) // 2 - 1]) / 2
    else:
        return sdata[len(data) // 2]


def make_keys(peaks):
    return [data.midpoint for data in peaks]


def get_window(chromosome, target_peaks, distance):
    """
    Returns a window of all peaks from a replicate within a certain distance of
    a peak from another replicate.
    """
    lower = list(target_peaks)[0].midpoint
    upper = list(target_peaks)[0].midpoint
    for peak in target_peaks:
        lower = min(lower, peak.midpoint - distance)
        upper = max(upper, peak.midpoint + distance)
    start_index = bisect.bisect_left(chromosome.keys, lower)
    end_index = bisect.bisect_right(chromosome.keys, upper)
    return (chromosome.peaks[start_index: end_index], chromosome.name)


def match_largest(window, peak, chrum):
    if not window:
        return None
    if peak.chrom != chrum:
        return None
    return max(window, key=lambda cpeak: cpeak.value)


def match_closest(window, peak, chrum):
    if not window:
        return None
    if peak.chrom != chrum:
        return None
    return min(window, key=lambda match: abs(match.midpoint - peak.midpoint))


def frequency_histogram(freqs, dataset_path, labels=[], title=''):
    pyplot.clf()
    pyplot.figure(figsize=(10, 10))
    for i, freq in enumerate(freqs):
        xvals, yvals = freq.graph_series()
        # Go from high to low
        xvals.reverse()
        pyplot.bar([x - 0.4 + 0.8 / len(freqs) * i for x in xvals], yvals, width=0.8 / len(freqs), color=COLORS[i])
    pyplot.xticks(range(min(xvals), max(xvals) + 1), map(str, reversed(range(min(xvals), max(xvals) + 1))))
    pyplot.xlabel(X_LABEL)
    pyplot.ylabel(Y_LABEL)
    pyplot.subplots_adjust(left=ADJUST[0], right=ADJUST[1], top=ADJUST[2], bottom=ADJUST[3])
    ax = pyplot.gca()
    for l in ax.get_xticklines() + ax.get_yticklines():
        l.set_markeredgewidth(TICK_WIDTH)
    pyplot.savefig(dataset_path)


METHODS = {'closest': match_closest, 'largest': match_largest}


def gff_attrs(l):
    if len(l) == 0:
        return '.'
    return ';'.join('%s=%s' % (tup[0], tup[1]) for tup in l)


def parse_gff_attrs(s):
    d = {}
    if s == '.':
        return d
    for item in s.split(';'):
        key, val = item.split('=')
        d[key] = val
    return d


def gff_row(cname, start, end, score, source, stype='.', strand='.', phase='.', attrs=None):
    return (cname, source, stype, start, end, score, strand, phase, gff_attrs(attrs or []))


def get_temporary_plot_path():
    """
    Return the path to a temporary file with a valid image format
    file extension that can be used with bioformats.
    """
    tmp_dir = tempfile.mkdtemp(prefix='tmp-repmatch-')
    fd, name = tempfile.mkstemp(suffix='.pdf', dir=tmp_dir)
    os.close(fd)
    return name


def process_files(dataset_paths, galaxy_hids, method, distance, step, replicates, up_limit, low_limit, output_files,
                  output_matched_peaks, output_unmatched_peaks, output_detail, output_statistics_table, output_statistics_histogram):
    output_statistics_histogram_file = output_files in ["all"] and method in ["all"]
    if len(dataset_paths) < 2:
        return
    if method == 'all':
        match_methods = METHODS.keys()
    else:
        match_methods = [method]
    for match_method in match_methods:
        statistics = perform_process(dataset_paths,
                                     galaxy_hids,
                                     match_method,
                                     distance,
                                     step,
                                     replicates,
                                     up_limit,
                                     low_limit,
                                     output_files,
                                     output_matched_peaks,
                                     output_unmatched_peaks,
                                     output_detail,
                                     output_statistics_table,
                                     output_statistics_histogram)
    if output_statistics_histogram_file:
        tmp_statistics_histogram_path = get_temporary_plot_path()
        frequency_histogram([stat['distribution'] for stat in [statistics]],
                            tmp_statistics_histogram_path,
                            METHODS.keys())
        shutil.move(tmp_statistics_histogram_path, output_statistics_histogram)


def perform_process(dataset_paths, galaxy_hids, method, distance, step, num_required, up_limit, low_limit, output_files,
                    output_matched_peaks, output_unmatched_peaks, output_detail, output_statistics_table, output_statistics_histogram):
    output_detail_file = output_files in ["all"] and output_detail is not None
    output_statistics_table_file = output_files in ["all"] and output_statistics_table is not None
    output_unmatched_peaks_file = output_files in ["all", "matched_peaks_unmatched_peaks"] and output_unmatched_peaks is not None
    output_statistics_histogram_file = output_files in ["all"] and output_statistics_histogram is not None
    replicates = []
    for i, dataset_path in enumerate(dataset_paths):
        try:
            galaxy_hid = galaxy_hids[i]
            r = Replicate(galaxy_hid, dataset_path)
            replicates.append(r)
        except Exception as e:
            stop_err('Unable to parse file "%s", exception: %s' % (dataset_path, str(e)))
    attrs = 'd%sr%s' % (distance, num_required)
    if up_limit != 1000:
        attrs += 'u%d' % up_limit
    if low_limit != -1000:
        attrs += 'l%d' % low_limit
    if step != 0:
        attrs += 's%d' % step

    def td_writer(file_path):
        # Returns a tab-delimited writer for a certain output
        if ISPY2:
            fh = open(file_path, 'wb')
            return csv.writer(fh, delimiter='\t')
        else:
            fh = open(file_path, 'w', newline='')
            return csv.writer(fh, delimiter='\t', quoting=csv.QUOTE_NONE)

    labels = ('chrom',
              'median midpoint',
              'median midpoint+1',
              'median normalized reads',
              'replicates',
              'median c-w distance',
              'reads sum')
    for replicate in replicates:
        labels += ('chrom',
                   'median midpoint',
                   'median midpoint+1',
                   'c-w sum',
                   'c-w distance',
                   'replicate id')
    matched_peaks_output = td_writer(output_matched_peaks)
    if output_statistics_table_file:
        statistics_table_output = td_writer(output_statistics_table)
        statistics_table_output.writerow(('data', 'median read count'))
    if output_detail_file:
        detail_output = td_writer(output_detail)
        detail_output.writerow(labels)
    if output_unmatched_peaks_file:
        unmatched_peaks_output = td_writer(output_unmatched_peaks)
        unmatched_peaks_output.writerow(('chrom', 'midpoint', 'midpoint+1', 'c-w sum', 'c-w distance', 'replicate id'))
    # Perform filtering
    if up_limit < 1000 or low_limit > -1000:
        for replicate in replicates:
            replicate.filter(up_limit, low_limit)
    # Actually merge the peaks
    peak_groups = []
    unmatched_peaks = []
    freq = FrequencyDistribution()

    def do_match(reps, distance):
        # Copy list because we will mutate it, but keep replicate references.
        reps = reps[:]
        while len(reps) > 1:
            # Iterate over each replicate as "main"
            main = reps[0]
            reps.remove(main)
            for chromosome in list(main.chromosomes.values()):
                peaks_by_value = chromosome.peaks[:]
                # Sort main replicate by value
                peaks_by_value.sort(key=lambda peak: -peak.value)

                def search_for_matches(group):
                    # Here we use multiple passes, expanding the window to be
                    #  +- distance from any previously matched peak.
                    while True:
                        new_match = False
                        for replicate in reps:
                            if replicate.id in group.peaks:
                                # Stop if match already found for this replicate
                                continue
                            try:
                                # Lines changed to remove a major bug by Rohit Reja.
                                window, chrum = get_window(replicate.chromosomes[chromosome.name], list(group.peaks.values()), distance)
                                match = METHODS[method](window, peak, chrum)
                            except KeyError:
                                continue
                            if match:
                                group.add_peak(replicate.id, match)
                                new_match = True
                        if not new_match:
                            break
                # Attempt to enlarge existing peak groups
                for group in peak_groups:
                    old_peaks = list(group.peaks.values())
                    search_for_matches(group)
                    for peak in list(group.peaks.values()):
                        if peak not in old_peaks:
                            peak.replicate.chromosomes[chromosome.name].remove_peak(peak)
                # Attempt to find new peaks groups.  For each peak in the
                # main replicate, search for matches in the other replicates
                for peak in peaks_by_value:
                    matches = PeakGroup()
                    matches.add_peak(main.id, peak)
                    search_for_matches(matches)
                    # Were enough replicates matched?
                    if matches.num_replicates >= num_required:
                        for peak in list(matches.peaks.values()):
                            peak.replicate.chromosomes[chromosome.name].remove_peak(peak)
                        peak_groups.append(matches)
    # Zero or less = no stepping
    if step <= 0:
        do_match(replicates, distance)
    else:
        for d in range(0, distance, step):
            do_match(replicates, d)
    for group in peak_groups:
        freq.add(group.num_replicates)
    # Collect together the remaining unmatched_peaks
    for replicate in replicates:
        for chromosome in replicate.chromosomes.values():
            for peak in chromosome.peaks:
                freq.add(1)
                unmatched_peaks.append(peak)
    # Average the unmatched_peaks count in the graph by # replicates
    med = median([peak.value for group in peak_groups for peak in group.peaks.values()])
    for replicate in replicates:
        replicate.median = median([peak.value for group in peak_groups for peak in group.peaks.values() if peak.replicate == replicate])
        statistics_table_output.writerow((replicate.id, replicate.median))
    for group in peak_groups:
        # Output matched_peaks (matched pairs).
        matched_peaks_output.writerow(gff_row(cname=group.chrom,
                                              start=group.midpoint,
                                              end=group.midpoint + 1,
                                              score=group.normalized_value(med),
                                              source='repmatch',
                                              stype='.',
                                              strand='.',
                                              phase='.',
                                              attrs=[('median_distance', group.median_distance),
                                                     ('value_sum', group.value_sum),
                                                     ('replicates', group.num_replicates)]))
        if output_detail_file:
            matched_peaks = (group.chrom,
                             group.midpoint,
                             group.midpoint + 1,
                             group.normalized_value(med),
                             group.num_replicates,
                             group.median_distance,
                             group.value_sum)
            for peak in group.peaks.values():
                matched_peaks += (peak.chrom, peak.midpoint, peak.midpoint + 1, peak.value, peak.distance, peak.replicate.id)
            detail_output.writerow(matched_peaks)
    if output_unmatched_peaks_file:
        for unmatched_peak in unmatched_peaks:
            unmatched_peaks_output.writerow((unmatched_peak.chrom,
                                             unmatched_peak.midpoint,
                                             unmatched_peak.midpoint + 1,
                                             unmatched_peak.value,
                                             unmatched_peak.distance,
                                             unmatched_peak.replicate.id))
    if output_statistics_histogram_file:
        tmp_statistics_histogram_path = get_temporary_plot_path()
        frequency_histogram([freq], tmp_statistics_histogram_path)
        shutil.move(tmp_statistics_histogram_path, output_statistics_histogram)
    return {'distribution': freq}
