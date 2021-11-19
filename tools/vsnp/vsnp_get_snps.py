#!/usr/bin/env python

# Collect quality parsimonious SNPs from vcf files
# and output alignment files in fasta format.

import argparse
import multiprocessing
import os
import queue
import shutil
import sys
import time
from collections import OrderedDict
from datetime import datetime

import pandas
import vcf


def get_time_stamp():
    return datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H-%M-%S')


def set_num_cpus(num_files, processes):
    num_cpus = int(multiprocessing.cpu_count())
    if num_files < num_cpus and num_files < processes:
        return num_files
    if num_cpus < processes:
        half_cpus = int(num_cpus / 2)
        if num_files < half_cpus:
            return num_files
        return half_cpus
    return processes


def setup_all_vcfs(vcf_files, vcf_dirs):
    # Create the all_vcfs directory and link
    # all input vcf files into it for processing.
    all_vcfs_dir = 'all_vcf'
    os.makedirs(all_vcfs_dir)
    vcf_dirs.append(all_vcfs_dir)
    for vcf_file in vcf_files:
        file_name_base = os.path.basename(vcf_file)
        dst_file = os.path.join(all_vcfs_dir, file_name_base)
        os.symlink(vcf_file, dst_file)
    return vcf_dirs


class SnpFinder:

    def __init__(self, num_files, dbkey, input_excel, all_isolates, ac, min_mq, quality_score_n_threshold, min_quality_score, input_vcf_dir, output_json_avg_mq_dir, output_json_snps_dir, output_snps_dir, output_summary):
        # Allele count
        self.ac = ac
        # Create a group that will contain all isolates.
        self.all_isolates = all_isolates
        # Evolving positions dictionary.
        self.all_positions = None
        # Isolate groups.
        self.groups = []
        # Excel file for grouping.
        self.input_excel = input_excel
        # Directory of input zero coverage vcf files.
        self.input_vcf_dir = input_vcf_dir
        # Minimum map quality value.
        self.min_mq = min_mq
        # Minimum quality score value.
        self.min_quality_score = min_quality_score
        # Number of input zero coverage vcf files.
        self.num_files = num_files
        # Output directory for json average mq files.
        self.output_json_avg_mq_dir = output_json_avg_mq_dir
        # Output directory for json snps files.
        self.output_json_snps_dir = output_json_snps_dir
        # Output directory for snps files.
        self.output_snps_dir = output_snps_dir
        # Quality score N threshold value.
        self.quality_score_n_threshold = quality_score_n_threshold
        self.dbkey = dbkey
        self.start_time = get_time_stamp()
        self.summary_str = ""
        self.timer_start = datetime.now()
        self.initiate_summary(output_summary)

    def append_to_summary(self, html_str):
        # Append a string to the html summary output file.
        self.summary_str = "%s%s" % (self.summary_str, html_str)

    def bin_input_files(self, filename, samples_groups_dict, defining_snps, inverted_defining_snps, found_positions, found_positions_mix):
        # Categorize input files into closely related
        # isolate groups based on discovered SNPs, and
        # return a group dictionary.
        sample_groups_list = []
        table_name = self.get_sample_name(filename)
        defining_snp = False
        # Absolute positions in set union of two lists.
        for abs_position in list(defining_snps.keys() & (found_positions.keys() | found_positions_mix.keys())):
            group = defining_snps[abs_position]
            sample_groups_list.append(group)
            self.check_add_group(group)
            if len(list(defining_snps.keys() & found_positions_mix.keys())) > 0:
                table_name = self.get_sample_name(filename)
                table_name = '%s<font color="red">[[MIXED]]</font>' % table_name
            self.copy_file(filename, group)
            defining_snp = True
        if not set(inverted_defining_snps.keys()).intersection(found_positions.keys() | found_positions_mix.keys()):
            for abs_position in list(inverted_defining_snps.keys()):
                group = inverted_defining_snps[abs_position]
                sample_groups_list.append(group)
                self.check_add_group(group)
                self.copy_file(filename, group)
                defining_snp = True
        if defining_snp:
            samples_groups_dict[table_name] = sorted(sample_groups_list)
        else:
            samples_groups_dict[table_name] = ['<font color="red">No defining SNP</font>']
        return samples_groups_dict

    def check_add_group(self, group):
        # Add a group if it is npt already in the list.
        if group not in self.groups:
            self.groups.append(group)

    def copy_file(self, filename, dir):
        if not os.path.exists(dir):
            os.makedirs(dir)
        shutil.copy(filename, dir)

    def decide_snps(self, filename):
        # Find the SNPs in a vcf file to produce a pandas data
        # frame and a dictionary containing sample map qualities.
        positions_dict = self.all_positions
        sample_map_qualities = {}
        # Eliminate the path.
        file_name_base = self.get_sample_name(filename)
        vcf_reader = vcf.Reader(open(filename, 'r'))
        sample_dict = {}
        for record in vcf_reader:
            alt = str(record.ALT[0])
            record_position = "%s:%s" % (str(record.CHROM), str(record.POS))
            if record_position in positions_dict:
                if alt == "None":
                    sample_dict.update({record_position: "-"})
                else:
                    # On rare occassions MQM gets called "NaN", thus passing
                    # a string when a number is expected when calculating average.
                    mq_val = self.get_mq_val(record.INFO, filename)
                    if str(mq_val).lower() not in ["nan"]:
                        sample_map_qualities.update({record_position: mq_val})
                    if len(alt) == 1:
                        qual_val = self.val_as_int(record.QUAL)
                        ac = record.INFO['AC'][0]
                        ref = str(record.REF[0])
                        if ac == 2 and qual_val > self.quality_score_n_threshold:
                            # Add the SNP to a group.
                            sample_dict.update({record_position: alt})
                        elif ac == 1 and qual_val > self.quality_score_n_threshold:
                            # The position is ambiguous.
                            alt_ref = "%s%s" % (alt, ref)
                            if alt_ref == "AG":
                                sample_dict.update({record_position: "R"})
                            elif alt_ref == "CT":
                                sample_dict.update({record_position: "Y"})
                            elif alt_ref == "GC":
                                sample_dict.update({record_position: "S"})
                            elif alt_ref == "AT":
                                sample_dict.update({record_position: "W"})
                            elif alt_ref == "GT":
                                sample_dict.update({record_position: "K"})
                            elif alt_ref == "AC":
                                sample_dict.update({record_position: "M"})
                            elif alt_ref == "GA":
                                sample_dict.update({record_position: "R"})
                            elif alt_ref == "TC":
                                sample_dict.update({record_position: "Y"})
                            elif alt_ref == "CG":
                                sample_dict.update({record_position: "S"})
                            elif alt_ref == "TA":
                                sample_dict.update({record_position: "W"})
                            elif alt_ref == "TG":
                                sample_dict.update({record_position: "K"})
                            elif alt_ref == "CA":
                                sample_dict.update({record_position: "M"})
                            else:
                                sample_dict.update({record_position: "N"})
                            # Poor calls
                        elif qual_val <= 50:
                            # Call the reference allele.
                            # Do not coerce record.REF[0] to a string!
                            sample_dict.update({record_position: record.REF[0]})
                        elif qual_val <= self.quality_score_n_threshold:
                            sample_dict.update({record_position: "N"})
                        else:
                            # Insurance -- Will still report on a possible
                            # SNP even if missed with above statements.
                            # Do not coerce record.REF[0] to a string!
                            sample_dict.update({record_position: record.REF[0]})
        # Merge dictionaries and order
        merge_dict = {}
        merge_dict.update(positions_dict)
        merge_dict.update(sample_dict)
        sample_df = pandas.DataFrame(merge_dict, index=[file_name_base])
        return sample_df, file_name_base, sample_map_qualities

    def df_to_fasta(self, parsimonious_df, group):
        # Generate SNP alignment file from
        # the parsimonious_df data frame.
        snps_file = os.path.join(self.output_snps_dir, "%s.fasta" % group)
        test_duplicates = []
        has_sequence_data = False
        for index, row in parsimonious_df.iterrows():
            for pos in row:
                if len(pos) > 0:
                    has_sequence_data = True
                    break
        if has_sequence_data:
            with open(snps_file, 'w') as fh:
                for index, row in parsimonious_df.iterrows():
                    test_duplicates.append(row.name)
                    if test_duplicates.count(row.name) < 2:
                        print(f'>{row.name}', file=fh)
                        for pos in row:
                            print(pos, end='', file=fh)
                        print("", file=fh)
        return has_sequence_data

    def find_initial_positions(self, filename):
        # Find SNP positions in a vcf file.
        found_positions = {}
        found_positions_mix = {}
        vcf_reader = vcf.Reader(open(filename, 'r'))
        for record in vcf_reader:
            qual_val = self.val_as_int(record.QUAL)
            chrom = record.CHROM
            position = record.POS
            absolute_position = "%s:%s" % (str(chrom), str(position))
            alt = str(record.ALT[0])
            if alt != "None":
                mq_val = self.get_mq_val(record.INFO, filename)
                ac = record.INFO['AC'][0]
                if ac == self.ac and len(record.REF) == 1 and qual_val > self.min_quality_score and mq_val > self.min_mq:
                    found_positions.update({absolute_position: record.REF})
                if ac == 1 and len(record.REF) == 1 and qual_val > self.min_quality_score and mq_val > self.min_mq:
                    found_positions_mix.update({absolute_position: record.REF})
        return found_positions, found_positions_mix

    def gather_and_filter(self, prefilter_df, mq_averages, group_dir):
        # Group a data frame of SNPs.
        if self.input_excel is None:
            filtered_all_df = prefilter_df
            sheet_names = None
        else:
            # Filter positions to be removed from all.
            xl = pandas.ExcelFile(self.input_excel)
            sheet_names = xl.sheet_names
            # Use the first column to filter "all" postions.
            exclusion_list_all = self.get_position_list(sheet_names, 0)
            exclusion_list_group = self.get_position_list(sheet_names, group_dir)
            exclusion_list = exclusion_list_all + exclusion_list_group
            # Filters for all applied.
            filtered_all_df = prefilter_df.drop(columns=exclusion_list, errors='ignore')
        json_snps_file = os.path.join(self.output_json_snps_dir, "%s.json" % group_dir)
        parsimonious_df = self.get_parsimonious_df(filtered_all_df)
        samples_number, columns = parsimonious_df.shape
        if samples_number >= 4:
            # Sufficient samples have been found
            # to build a phylogenetic tree.
            has_sequence_data = self.df_to_fasta(parsimonious_df, group_dir)
            if has_sequence_data:
                json_avg_mq_file = os.path.join(self.output_json_avg_mq_dir, "%s.json" % group_dir)
                mq_averages.to_json(json_avg_mq_file, orient='split')
                parsimonious_df.to_json(json_snps_file, orient='split')
            else:
                msg = "<br/>No sequence data"
                if group_dir is not None:
                    msg = "%s for group: %s" % (msg, group_dir)
                self.append_to_summary("%s<br/>\n" % msg)
        else:
            msg = "<br/>Too few samples to build tree"
            if group_dir is not None:
                msg = "%s for group: %s" % (msg, group_dir)
            self.append_to_summary("%s<br/>\n" % msg)

    def get_sample_name(self, file_path):
        # Return the sample part of a file name.
        base_file_name = os.path.basename(file_path)
        if base_file_name.find(".") > 0:
            # Eliminate the extension.
            return os.path.splitext(base_file_name)[0]
        return base_file_name

    def get_mq_val(self, record_info, filename):
        # Get the MQ (gatk) or MQM (freebayes) value
        # from the record.INFO component of the vcf file.
        try:
            mq_val = record_info['MQM']
            return self.return_val(mq_val)
        except Exception:
            try:
                mq_val = record_info['MQ']
                return self.return_val(mq_val)
            except Exception:
                msg = "Invalid or unsupported vcf header %s in file: %s\n" % (str(record_info), filename)
                sys.exit(msg)

    def get_parsimonious_df(self, filtered_all_df):
        # Get the parsimonious SNPs data frame
        # from a data frame of filtered SNPs.
        try:
            ref_series = filtered_all_df.loc['root']
            # In all_vcf root needs to be removed.
            filtered_all_df = filtered_all_df.drop(['root'])
        except KeyError:
            pass
        parsimony = filtered_all_df.loc[:, (filtered_all_df != filtered_all_df.iloc[0]).any()]
        parsimony_positions = list(parsimony)
        parse_df = filtered_all_df[parsimony_positions]
        ref_df = ref_series.to_frame()
        ref_df = ref_df.T
        parsimonious_df = pandas.concat([parse_df, ref_df], join='inner')
        return parsimonious_df

    def get_position_list(self, sheet_names, group):
        # Get a list of positions defined by an excel file.
        exclusion_list = []
        try:
            filter_to_all = pandas.read_excel(self.input_excel, header=1, usecols=[group])
            for value in filter_to_all.values:
                value = str(value[0])
                if "-" not in value.split(":")[-1]:
                    exclusion_list.append(value)
                elif "-" in value:
                    try:
                        chrom, sequence_range = value.split(":")
                    except Exception as e:
                        sys.exit(str(e))
                    value = sequence_range.split("-")
                    for position in range(int(value[0].replace(',', '')), int(value[1].replace(',', '')) + 1):
                        exclusion_list.append(chrom + ":" + str(position))
            return exclusion_list
        except ValueError:
            return []

    def get_snps(self, task_queue, timeout):
        while True:
            try:
                group_dir = task_queue.get(block=True, timeout=timeout)
            except queue.Empty:
                break
            # Parse all vcf files to accumulate
            # the SNPs into a data frame.
            positions_dict = {}
            group_files = []
            for file_name in os.listdir(os.path.abspath(group_dir)):
                file_path = os.path.abspath(os.path.join(group_dir, file_name))
                group_files.append(file_path)
            for file_name in group_files:
                found_positions, found_positions_mix = self.find_initial_positions(file_name)
                positions_dict.update(found_positions)
            # Order before adding to file to match
            # with ordering of individual samples.
            # all_positions is abs_pos:REF
            self.all_positions = OrderedDict(sorted(positions_dict.items()))
            ref_positions_df = pandas.DataFrame(self.all_positions, index=['root'])
            all_map_qualities = {}
            df_list = []
            for file_name in group_files:
                sample_df, file_name_base, sample_map_qualities = self.decide_snps(file_name)
                df_list.append(sample_df)
                all_map_qualities.update({file_name_base: sample_map_qualities})
            all_sample_df = pandas.concat(df_list)
            # All positions have now been selected for each sample,
            # so select parisomony informative SNPs.  This removes
            # columns where all fields are the same.
            # Add reference to top row.
            prefilter_df = pandas.concat([ref_positions_df, all_sample_df], join='inner')
            all_mq_df = pandas.DataFrame.from_dict(all_map_qualities)
            mq_averages = all_mq_df.mean(axis=1).astype(int)
            self.gather_and_filter(prefilter_df, mq_averages, group_dir)
            task_queue.task_done()

    def group_vcfs(self, vcf_files):
        # Parse an excel file to produce a
        # grouping dictionary for SNPs.
        xl = pandas.ExcelFile(self.input_excel)
        sheet_names = xl.sheet_names
        ws = pandas.read_excel(self.input_excel, sheet_name=sheet_names[0])
        defining_snps = ws.iloc[0]
        defsnp_iterator = iter(defining_snps.iteritems())
        next(defsnp_iterator)
        defining_snps = {}
        inverted_defining_snps = {}
        for abs_pos, group in defsnp_iterator:
            if '!' in abs_pos:
                inverted_defining_snps[abs_pos.replace('!', '')] = group
            else:
                defining_snps[abs_pos] = group
        samples_groups_dict = {}
        for vcf_file in vcf_files:
            found_positions, found_positions_mix = self.find_initial_positions(vcf_file)
            samples_groups_dict = self.bin_input_files(vcf_file, samples_groups_dict, defining_snps, inverted_defining_snps, found_positions, found_positions_mix)
        # Output summary grouping table.
        self.append_to_summary('<br/>')
        self.append_to_summary('<b>Groupings with %d listed:</b><br/>\n' % len(samples_groups_dict))
        self.append_to_summary('<table  cellpadding="5" cellspaging="5" border="1">\n')
        for key, value in samples_groups_dict.items():
            self.append_to_summary('<tr align="left"><th>Sample Name</th>\n')
            self.append_to_summary('<td>%s</td>' % key)
            for group in value:
                self.append_to_summary('<td>%s</td>\n' % group)
            self.append_to_summary('</tr>\n')
        self.append_to_summary('</table><br/>\n')

    def initiate_summary(self, output_summary):
        # Output summary file handle.
        self.append_to_summary('<html>\n')
        self.append_to_summary('<head></head>\n')
        self.append_to_summary('<body style=\"font-size:12px;">')
        self.append_to_summary("<b>Time started:</b> %s<br/>" % get_time_stamp())
        self.append_to_summary("<b>Number of VCF inputs:</b> %d<br/>" % self.num_files)
        self.append_to_summary("<b>Reference:</b> %s<br/>" % self.dbkey)
        self.append_to_summary("<b>All isolates:</b> %s<br/>" % str(self.all_isolates))

    def return_val(self, val, index=0):
        # Handle element and single-element list values.
        if isinstance(val, list):
            return val[index]
        return val

    def val_as_int(self, val):
        # Handle integer value conversion.
        try:
            return int(val)
        except TypeError:
            # val is likely None here.
            return 0


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('--ac', action='store', dest='ac', type=int, help='Allele count value'),
    parser.add_argument('--all_isolates', action='store_true', dest='all_isolates', help='Create table with all isolates'),
    parser.add_argument('--input_excel', action='store', dest='input_excel', required=False, default=None, help='Optional Excel filter file'),
    parser.add_argument('--input_vcf_dir', action='store', dest='input_vcf_dir', help='Input vcf directory'),
    parser.add_argument('--min_mq', action='store', dest='min_mq', type=int, help='Minimum map quality value'),
    parser.add_argument('--min_quality_score', action='store', dest='min_quality_score', type=int, help='Minimum quality score value'),
    parser.add_argument('--output_json_avg_mq_dir', action='store', dest='output_json_avg_mq_dir', help='Output json average mq directory'),
    parser.add_argument('--output_json_snps_dir', action='store', dest='output_json_snps_dir', help='Output json snps directory'),
    parser.add_argument('--output_snps_dir', action='store', dest='output_snps_dir', help='Output snps directory'),
    parser.add_argument('--output_summary', action='store', dest='output_summary', help='Output summary html file'),
    parser.add_argument('--processes', action='store', dest='processes', type=int, help='Configured processes for job'),
    parser.add_argument('--quality_score_n_threshold', action='store', dest='quality_score_n_threshold', type=int, help='Minimum quality score N value for alleles'),
    parser.add_argument('--dbkey', action='store', dest='dbkey', help='Galaxy genome build dbkey'),

    args = parser.parse_args()

    # Build the list of all input zero coverage vcf
    # files, both the samples and the "database".
    vcf_files = []
    for file_name in os.listdir(args.input_vcf_dir):
        file_path = os.path.abspath(os.path.join(args.input_vcf_dir, file_name))
        vcf_files.append(file_path)

    multiprocessing.set_start_method('spawn')
    queue1 = multiprocessing.JoinableQueue()
    num_files = len(vcf_files)
    cpus = set_num_cpus(num_files, args.processes)
    # Set a timeout for get()s in the queue.
    timeout = 0.05

    # Initialize the snp_finder object.
    snp_finder = SnpFinder(num_files, args.dbkey, args.input_excel, args.all_isolates, args.ac, args.min_mq, args.quality_score_n_threshold, args.min_quality_score, args.input_vcf_dir, args.output_json_avg_mq_dir, args.output_json_snps_dir, args.output_snps_dir, args.output_summary)

    # Define and make the set of directories into which the input_zc_vcf
    # files will be placed.  Selected input values (e.g., the use of
    # an Excel file for grouping and filtering, creating a group with
    # all isolates) are used to define the directories.
    vcf_dirs = []
    if args.input_excel is None:
        vcf_dirs = setup_all_vcfs(vcf_files, vcf_dirs)
    else:
        if args.all_isolates:
            vcf_dirs = setup_all_vcfs(vcf_files, vcf_dirs)
        # Parse the Excel file to detemine groups for filtering.
        snp_finder.group_vcfs(vcf_files)
        # Append the list of group directories created by
        # the above call to the set of directories containing
        # vcf files for analysis.
        group_dirs = [d for d in os.listdir(os.getcwd()) if os.path.isdir(d) and d in snp_finder.groups]
        vcf_dirs.extend(group_dirs)

    # Populate the queue for job splitting.
    for vcf_dir in vcf_dirs:
        queue1.put(vcf_dir)

    # Complete the get_snps task.
    processes = [multiprocessing.Process(target=snp_finder.get_snps, args=(queue1, timeout, )) for _ in range(cpus)]
    for p in processes:
        p.start()
    for p in processes:
        p.join()
    queue1.join()

    # Finish summary log.
    snp_finder.append_to_summary("<br/><b>Time finished:</b> %s<br/>\n" % get_time_stamp())
    total_run_time = datetime.now() - snp_finder.timer_start
    snp_finder.append_to_summary("<br/><b>Total run time:</b> %s<br/>\n" % str(total_run_time))
    snp_finder.append_to_summary('</body>\n</html>\n')
    with open(args.output_summary, "w") as fh:
        fh.write("%s" % snp_finder.summary_str)
