#!/usr/bin/env python

import argparse
import multiprocessing
import os
import queue
import re

import pandas
import pandas.io.formats.excel
from Bio import SeqIO

# Maximum columns allowed in a LibreOffice
# spreadsheet is 1024.  Excel allows for
# 16,384 columns, but we'll set the lower
# number as the maximum.  Some browsers
# (e.g., Firefox on Linux) are configured
# to use LibreOffice for Excel spreadsheets.
MAXCOLS = 1024
OUTPUT_EXCEL_DIR = 'output_excel_dir'
INPUT_JSON_AVG_MQ_DIR = 'input_json_avg_mq_dir'
INPUT_JSON_DIR = 'input_json_dir'
INPUT_NEWICK_DIR = 'input_newick_dir'


def annotate_table(table_df, group, annotation_dict):
    for gbk_chrome, pro in list(annotation_dict.items()):
        ref_pos = list(table_df)
        ref_series = pandas.Series(ref_pos)
        ref_df = pandas.DataFrame(ref_series.str.split(':', expand=True).values, columns=['reference', 'position'])
        all_ref = ref_df[ref_df['reference'] == gbk_chrome]
        positions = all_ref.position.to_frame()
        # Create an annotation file.
        annotation_file = "%s_annotations.csv" % group
        with open(annotation_file, "a") as fh:
            for _, row in positions.iterrows():
                pos = row.position
                try:
                    aaa = pro.iloc[pro.index.get_loc(int(pos))][['chrom', 'locus', 'product', 'gene']]
                    try:
                        chrom, name, locus, tag = aaa.values[0]
                        print("{}:{}\t{}, {}, {}".format(chrom, pos, locus, tag, name), file=fh)
                    except ValueError:
                        # If only one annotation for the entire
                        # chromosome (e.g., flu) then having [0] fails
                        chrom, name, locus, tag = aaa.values
                        print("{}:{}\t{}, {}, {}".format(chrom, pos, locus, tag, name), file=fh)
                except KeyError:
                    print("{}:{}\tNo annotated product".format(gbk_chrome, pos), file=fh)
    # Read the annotation file into a data frame.
    annotations_df = pandas.read_csv(annotation_file, sep='\t', header=None, names=['index', 'annotations'], index_col='index')
    # Remove the annotation_file from disk since both
    # cascade and sort tables are built using the file,
    # and it is opened for writing in append mode.
    os.remove(annotation_file)
    # Process the data.
    table_df_transposed = table_df.T
    table_df_transposed.index = table_df_transposed.index.rename('index')
    table_df_transposed = table_df_transposed.merge(annotations_df, left_index=True, right_index=True)
    table_df = table_df_transposed.T
    return table_df


def excel_formatter(json_file_name, excel_file_name, group, annotation_dict):
    pandas.io.formats.excel.header_style = None
    table_df = pandas.read_json(json_file_name, orient='split')
    if annotation_dict is not None:
        table_df = annotate_table(table_df, group, annotation_dict)
    else:
        table_df = table_df.append(pandas.Series(name='no annotations'))
    writer = pandas.ExcelWriter(excel_file_name, engine='xlsxwriter')
    table_df.to_excel(writer, sheet_name='Sheet1')
    writer_book = writer.book
    ws = writer.sheets['Sheet1']
    format_a = writer_book.add_format({'bg_color': '#58FA82'})
    format_g = writer_book.add_format({'bg_color': '#F7FE2E'})
    format_c = writer_book.add_format({'bg_color': '#0000FF'})
    format_t = writer_book.add_format({'bg_color': '#FF0000'})
    format_normal = writer_book.add_format({'bg_color': '#FDFEFE'})
    formatlowqual = writer_book.add_format({'font_color': '#C70039', 'bg_color': '#E2CFDD'})
    format_ambigous = writer_book.add_format({'font_color': '#C70039', 'bg_color': '#E2CFDD'})
    format_n = writer_book.add_format({'bg_color': '#E2CFDD'})
    rows, cols = table_df.shape
    ws.set_column(0, 0, 30)
    ws.set_column(1, cols, 2.1)
    ws.freeze_panes(2, 1)
    format_annotation = writer_book.add_format({'font_color': '#0A028C', 'rotation': '-90', 'align': 'top'})
    # Set last row.
    ws.set_row(rows + 1, cols + 1, format_annotation)
    # Make sure that row/column locations don't overlap.
    ws.conditional_format(rows - 2, 1, rows - 1, cols, {'type': 'cell', 'criteria': '<', 'value': 55, 'format': formatlowqual})
    ws.conditional_format(2, 1, rows - 2, cols, {'type': 'cell', 'criteria': '==', 'value': 'B$2', 'format': format_normal})
    ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'A', 'format': format_a})
    ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'G', 'format': format_g})
    ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'C', 'format': format_c})
    ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'T', 'format': format_t})
    ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'S', 'format': format_ambigous})
    ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'Y', 'format': format_ambigous})
    ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'R', 'format': format_ambigous})
    ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'W', 'format': format_ambigous})
    ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'K', 'format': format_ambigous})
    ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'M', 'format': format_ambigous})
    ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'N', 'format': format_n})
    ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': '-', 'format': format_n})
    format_rotation = writer_book.add_format({})
    format_rotation.set_rotation(90)
    for column_num, column_name in enumerate(list(table_df.columns)):
        ws.write(0, column_num + 1, column_name, format_rotation)
    format_annotation = writer_book.add_format({'font_color': '#0A028C', 'rotation': '-90', 'align': 'top'})
    # Set last row.
    ws.set_row(rows, 400, format_annotation)
    writer.save()


def get_annotation_dict(gbk_file):
    gbk_dict = SeqIO.to_dict(SeqIO.parse(gbk_file, "genbank"))
    annotation_dict = {}
    tmp_file = "features.csv"
    # Create a file of chromosomes and features.
    for chromosome in list(gbk_dict.keys()):
        with open(tmp_file, 'w+') as fh:
            for feature in gbk_dict[chromosome].features:
                if "CDS" in feature.type or "rRNA" in feature.type:
                    try:
                        product = feature.qualifiers['product'][0]
                    except KeyError:
                        product = None
                    try:
                        locus = feature.qualifiers['locus_tag'][0]
                    except KeyError:
                        locus = None
                    try:
                        gene = feature.qualifiers['gene'][0]
                    except KeyError:
                        gene = None
                    fh.write("%s\t%d\t%d\t%s\t%s\t%s\n" % (chromosome, int(feature.location.start), int(feature.location.end), locus, product, gene))
        # Read the chromosomes and features file into a data frame.
        df = pandas.read_csv(tmp_file, sep='\t', names=["chrom", "start", "stop", "locus", "product", "gene"])
        # Process the data.
        df = df.sort_values(['start', 'gene'], ascending=[True, False])
        df = df.drop_duplicates('start')
        pro = df.reset_index(drop=True)
        pro.index = pandas.IntervalIndex.from_arrays(pro['start'], pro['stop'], closed='both')
        annotation_dict[chromosome] = pro
    return annotation_dict


def get_sample_name(file_path):
    base_file_name = os.path.basename(file_path)
    if base_file_name.find(".") > 0:
        # Eliminate the extension.
        return os.path.splitext(base_file_name)[0]
    return base_file_name


def output_cascade_table(cascade_order, mqdf, group, annotation_dict):
    cascade_order_mq = pandas.concat([cascade_order, mqdf], join='inner')
    output_table(cascade_order_mq, "cascade", group, annotation_dict)


def output_excel(df, type_str, group, annotation_dict, count=None):
    # Output the temporary json file that
    # is used by the excel_formatter.
    if count is None:
        if group is None:
            json_file_name = os.path.join(OUTPUT_EXCEL_DIR, "%s_order_mq.json" % type_str)
            excel_file_name = os.path.join(OUTPUT_EXCEL_DIR, "%s_table.xlsx" % type_str)
        else:
            json_file_name = os.path.join(OUTPUT_EXCEL_DIR, "%s_%s_order_mq.json" % (group, type_str))
            excel_file_name = os.path.join(OUTPUT_EXCEL_DIR, "%s_%s_table.xlsx" % (group, type_str))
    else:
        # The table has more columns than is allowed by the
        # MAXCOLS setting, so multiple files will be produced
        # as an output collection.
        if group is None:
            json_file_name = os.path.join(OUTPUT_EXCEL_DIR, "%s_order_mq_%d.json" % (type_str, count))
            excel_file_name = os.path.join(OUTPUT_EXCEL_DIR, "%s_table_%d.xlsx" % (type_str, count))
        else:
            json_file_name = os.path.join(OUTPUT_EXCEL_DIR, "%s_%s_order_mq_%d.json" % (group, type_str, count))
            excel_file_name = os.path.join(OUTPUT_EXCEL_DIR, "%s_%s_table_%d.xlsx" % (group, type_str, count))
    df.to_json(json_file_name, orient='split')
    # Output the Excel file.
    excel_formatter(json_file_name, excel_file_name, group, annotation_dict)


def output_sort_table(cascade_order, mqdf, group, annotation_dict):
    sort_df = cascade_order.T
    sort_df['abs_value'] = sort_df.index
    sort_df[['chrom', 'pos']] = sort_df['abs_value'].str.split(':', expand=True)
    sort_df = sort_df.drop(['abs_value', 'chrom'], axis=1)
    sort_df.pos = sort_df.pos.astype(int)
    sort_df = sort_df.sort_values(by=['pos'])
    sort_df = sort_df.drop(['pos'], axis=1)
    sort_df = sort_df.T
    sort_order_mq = pandas.concat([sort_df, mqdf], join='inner')
    output_table(sort_order_mq, "sort", group, annotation_dict)


def output_table(df, type_str, group, annotation_dict):
    if isinstance(group, str) and group.startswith("dataset"):
        # Inputs are single files, not collections,
        # so input file names are not useful for naming
        # output files.
        group_str = None
    else:
        group_str = group
    count = 0
    chunk_start = 0
    chunk_end = 0
    column_count = df.shape[1]
    if column_count >= MAXCOLS:
        # Here the number of columns is greater than
        # the maximum allowed by Excel, so multiple
        # outputs will be produced.
        while column_count >= MAXCOLS:
            count += 1
            chunk_end += MAXCOLS
            df_of_type = df.iloc[:, chunk_start:chunk_end]
            output_excel(df_of_type, type_str, group_str, annotation_dict, count=count)
            chunk_start += MAXCOLS
            column_count -= MAXCOLS
        count += 1
        df_of_type = df.iloc[:, chunk_start:]
        output_excel(df_of_type, type_str, group_str, annotation_dict, count=count)
    else:
        output_excel(df, type_str, group_str, annotation_dict)


def preprocess_tables(task_queue, annotation_dict, timeout):
    while True:
        try:
            tup = task_queue.get(block=True, timeout=timeout)
        except queue.Empty:
            break
        newick_file, json_file, json_avg_mq_file = tup
        avg_mq_series = pandas.read_json(json_avg_mq_file, typ='series', orient='split')
        # Map quality to dataframe.
        mqdf = avg_mq_series.to_frame(name='MQ')
        mqdf = mqdf.T
        # Get the group.
        group = get_sample_name(newick_file)
        snps_df = pandas.read_json(json_file, orient='split')
        with open(newick_file, 'r') as fh:
            for line in fh:
                line = re.sub('[:,]', '\n', line)
                line = re.sub('[)(]', '', line)
                line = re.sub(r'[0-9].*\.[0-9].*\n', '', line)
                line = re.sub('root\n', '', line)
        sample_order = line.split('\n')
        sample_order = list([_f for _f in sample_order if _f])
        sample_order.insert(0, 'root')
        tree_order = snps_df.loc[sample_order]
        # Count number of SNPs in each column.
        snp_per_column = []
        for column_header in tree_order:
            count = 0
            column = tree_order[column_header]
            for element in column:
                if element != column[0]:
                    count = count + 1
            snp_per_column.append(count)
        row1 = pandas.Series(snp_per_column, tree_order.columns, name="snp_per_column")
        # Count number of SNPS from the
        # top of each column in the table.
        snp_from_top = []
        for column_header in tree_order:
            count = 0
            column = tree_order[column_header]
            # for each element in the column
            # skip the first element
            for element in column[1:]:
                if element == column[0]:
                    count = count + 1
                else:
                    break
            snp_from_top.append(count)
        row2 = pandas.Series(snp_from_top, tree_order.columns, name="snp_from_top")
        tree_order = tree_order.append([row1])
        tree_order = tree_order.append([row2])
        # In pandas=0.18.1 even this does not work:
        # abc = row1.to_frame()
        # abc = abc.T --> tree_order.shape (5, 18), abc.shape (1, 18)
        # tree_order.append(abc)
        # Continue to get error: "*** ValueError: all the input arrays must have same number of dimensions"
        tree_order = tree_order.T
        tree_order = tree_order.sort_values(['snp_from_top', 'snp_per_column'], ascending=[True, False])
        tree_order = tree_order.T
        # Remove snp_per_column and snp_from_top rows.
        cascade_order = tree_order[:-2]
        # Output the cascade table.
        output_cascade_table(cascade_order, mqdf, group, annotation_dict)
        # Output the sorted table.
        output_sort_table(cascade_order, mqdf, group, annotation_dict)
        task_queue.task_done()


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


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--input_avg_mq_json', action='store', dest='input_avg_mq_json', required=False, default=None, help='Average MQ json file')
    parser.add_argument('--input_newick', action='store', dest='input_newick', required=False, default=None, help='Newick file')
    parser.add_argument('--input_snps_json', action='store', dest='input_snps_json', required=False, default=None, help='SNPs json file')
    parser.add_argument('--gbk_file', action='store', dest='gbk_file', required=False, default=None, help='Optional gbk file'),
    parser.add_argument('--processes', action='store', dest='processes', type=int, help='User-selected number of processes to use for job splitting')

    args = parser.parse_args()

    if args.gbk_file is not None:
        # Create the annotation_dict for annotating
        # the Excel tables.
        annotation_dict = get_annotation_dict(args.gbk_file)
    else:
        annotation_dict = None

    # The assumption here is that the list of files
    # in both INPUT_NEWICK_DIR and INPUT_JSON_DIR are
    # named such that they are properly matched if
    # the directories contain more than 1 file (i.e.,
    # hopefully the newick file names and json file names
    # will be something like Mbovis-01D6_* so they can be
    # sorted and properly associated with each other).
    if args.input_newick is not None:
        newick_files = [args.input_newick]
    else:
        newick_files = []
        for file_name in sorted(os.listdir(INPUT_NEWICK_DIR)):
            file_path = os.path.abspath(os.path.join(INPUT_NEWICK_DIR, file_name))
            newick_files.append(file_path)
    if args.input_snps_json is not None:
        json_files = [args.input_snps_json]
    else:
        json_files = []
        for file_name in sorted(os.listdir(INPUT_JSON_DIR)):
            file_path = os.path.abspath(os.path.join(INPUT_JSON_DIR, file_name))
            json_files.append(file_path)
    if args.input_avg_mq_json is not None:
        json_avg_mq_files = [args.input_avg_mq_json]
    else:
        json_avg_mq_files = []
        for file_name in sorted(os.listdir(INPUT_JSON_AVG_MQ_DIR)):
            file_path = os.path.abspath(os.path.join(INPUT_JSON_AVG_MQ_DIR, file_name))
            json_avg_mq_files.append(file_path)

    multiprocessing.set_start_method('spawn')
    queue1 = multiprocessing.JoinableQueue()
    queue2 = multiprocessing.JoinableQueue()
    num_files = len(newick_files)
    cpus = set_num_cpus(num_files, args.processes)
    # Set a timeout for get()s in the queue.
    timeout = 0.05

    for i, newick_file in enumerate(newick_files):
        json_file = json_files[i]
        json_avg_mq_file = json_avg_mq_files[i]
        queue1.put((newick_file, json_file, json_avg_mq_file))

    # Complete the preprocess_tables task.
    processes = [multiprocessing.Process(target=preprocess_tables, args=(queue1, annotation_dict, timeout, )) for _ in range(cpus)]
    for p in processes:
        p.start()
    for p in processes:
        p.join()
    queue1.join()

    if queue1.empty():
        queue1.close()
        queue1.join_thread()
