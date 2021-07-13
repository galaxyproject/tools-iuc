import argparse
import pathlib
import sys

import xlrd
import yaml

FILE_FORMAT = 'fastq'


def extract_data(xl_sheet, expected_columns):
    """
    1. Check that the columns I expect are present in the sheet
    (any order and mixed with others, it's just a verification that
    the user filled the correct template)
    2. Fill a dictionary with the rows data indexed by first column in list"""
    sheet_columns = {}
    for sh_col in range(xl_sheet.ncols):
        if xl_sheet.cell(0, sh_col).value in expected_columns:
            if xl_sheet.cell(0, sh_col).value in sheet_columns.keys():
                sys.exit("Duplicated columns")
            else:
                sheet_columns[xl_sheet.cell(0, sh_col).value] = sh_col
    for col in range(len(expected_columns)):
        assert expected_columns[col] in sheet_columns.keys(), \
            "Expected column %s not found" % expected_columns[col]

    # fetch rows in a dict
    data_dict = {}
    # the first of the expected columns will be the index
    index_col = sheet_columns[expected_columns[0]]
    # skip first 2 rows: column names + comments rows
    for row_id in range(2, xl_sheet.nrows):
        row_dict = {}
        for col in range(1, len(expected_columns)):
            sheet_col_index = sheet_columns[expected_columns[col]]
            row_dict[expected_columns[col]] = xl_sheet.cell(row_id, sheet_col_index).value
        # should check for duplicate alias/ids?
        if xl_sheet.cell(row_id, index_col).value in data_dict.keys():
            tmp = data_dict[xl_sheet.cell(row_id, index_col).value]
            data_dict[xl_sheet.cell(row_id, index_col).value] = [tmp]
            data_dict[xl_sheet.cell(row_id, index_col).value].append(row_dict)
        else:
            data_dict[xl_sheet.cell(row_id, index_col).value] = row_dict
    return data_dict


def paste_xls2yaml(xlsx_path):
    print('YAML -------------')
    xls = xlrd.open_workbook(xlsx_path)
    content_dict = {}
    for sheet_name in xls.sheet_names():
        if sheet_name == 'controlled_vocabulary':
            continue
        xls_sheet = xls.sheet_by_name(sheet_name)
        sheet_contents_dict = {}
        colnames = []
        for col in range(xls_sheet.ncols):
            colnames.append(xls_sheet.cell(0, col).value)
        # skip first 2 rows (column names and suggestions)
        for row_id in range(2, xls_sheet.nrows):
            row_dict = {}
            for col_id in range(0, xls_sheet.ncols):
                row_dict[colnames[col_id]] = xls_sheet.cell(row_id, col_id).value
            # should check for duplicate alias/ids?
            sheet_contents_dict[row_id] = row_dict
        content_dict[sheet_name] = sheet_contents_dict
    yaml.dump(content_dict, sys.stdout)
    print('YAML -------------')


parser = argparse.ArgumentParser()
parser.add_argument('--form', dest='xlsx_path', required=True)
parser.add_argument('--out_dir', dest='out_path', required=True)
parser.add_argument('--action', dest='action', required=True)
parser.add_argument('--vir', dest='viral_submission', required=False, action='store_true')
parser.add_argument('--verbose', dest='verbose', required=False, action='store_true')
args = parser.parse_args()

xl_workbook = xlrd.open_workbook(args.xlsx_path)

# PARSE STUDIES
#################
xl_sheet = xl_workbook.sheet_by_name('ENA_study')
if xl_sheet.nrows < 3:
    raise ValueError('No entries found in studies sheet')
studies_dict = {}
studies_col = ['alias', 'title', 'study_type', 'study_abstract']
studies_dict = extract_data(xl_sheet, studies_col)

# PARSE SAMPLES
#################
xl_sheet = xl_workbook.sheet_by_name('ENA_sample')
if xl_sheet.nrows < 3:
    raise ValueError('No entries found in samples')
if args.viral_submission:
    samples_cols = ['alias', 'title', 'scientific_name', 'sample_description',
                    'geographic location (country and/or sea)', 'host common name',
                    'host health state', 'host sex', 'host scientific name', 'collector name',
                    'collection date', 'collecting institution', 'isolate']
else:
    samples_cols = ['alias', 'title', 'scientific_name', 'sample_description']
samples_dict = extract_data(xl_sheet, samples_cols)

# PARSE EXPERIMENTS
#################
xl_sheet = xl_workbook.sheet_by_name('ENA_experiment')
if xl_sheet.nrows < 3:
    raise ValueError('No experiments found in experiments sheet')
exp_columns = ['alias', 'title', 'study_alias', 'sample_alias', 'design_description',
               'library_name', 'library_strategy', 'library_source', 'library_selection',
               'library_layout', 'insert_size', 'library_construction_protocol',
               'platform', 'instrument_model']

experiments_dict = extract_data(xl_sheet, exp_columns)

# PARSE RUNS SHEET
#################
xl_sheet = xl_workbook.sheet_by_name('ENA_run')
if xl_sheet.nrows < 3:
    raise ValueError('No entries found in runs sheet')
run_cols = ['alias', 'experiment_alias', 'file_name', 'file_format']
runs_dict = extract_data(xl_sheet, run_cols)

# WRITE HEADERS TO TABLES
studies_table = open(pathlib.Path(args.out_path) / 'studies.tsv', 'w')
studies_table.write('\t'.join(['alias', 'status', 'accession', 'title', 'study_type',
                               'study_abstract', 'pubmed_id', 'submission_date']) + '\n')
samples_table = open(pathlib.Path(args.out_path) / 'samples.tsv', 'w')
if args.viral_submission:
    samples_table.write('\t'.join(['alias', 'status', 'accession', 'title', 'scientific_name',
                                   'taxon_id', 'sample_description', 'collection_date',
                                   'geographic_location', 'host_common_name', 'host_subject_id',
                                   'host_health_state', 'host_sex', 'host_scientific_name',
                                   'collector_name', 'collecting_institution', 'isolate',
                                   'submission_date']) + '\n')
else:
    samples_table.write('\t'.join(['alias', 'status', 'accession', 'title', 'scientific_name',
                                   'taxon_id', 'sample_description', 'submission_date']) + '\n')

experiments_table = open(pathlib.Path(args.out_path) / 'experiments.tsv', 'w')
experiments_table.write('\t'.join(['alias', 'status', 'accession', 'title', 'study_alias',
                                   'sample_alias', 'design_description', 'library_name',
                                   'library_strategy', 'library_source', 'library_selection',
                                   'library_layout', 'insert_size', 'library_construction_protocol',
                                   'platform', 'instrument_model', 'submission_date']) + '\n')

runs_table = open(pathlib.Path(args.out_path) / 'runs.tsv', 'w')
runs_table.write('\t'.join(['alias', 'status', 'accession', 'experiment_alias', 'file_name',
                            'file_format', 'file_checksum', 'submission_date']) + '\n')
action = args.action

# WRITE  DICTIONARIES TO TABLE FILES

# ADD A TIMESTAMP TO THE ALIAS? SEEMS LIKE ENA REQUIRES ALL ENTRIES FOR A WEBIN TO HAVE UNIQUE IDS?
# dt_oobj = datetime.now(tz=None)
# timestamp = dt_oobj.strftime("%Y%m%d_%H:%M:%S")
runs_included = []
exp_included = []
for study_alias, study in studies_dict.items():
    # study_alias = study_alias + '_' + timestamp
    studies_table.write('\t'.join([study_alias, action, 'ENA_accession', study['title'],
                                   study['study_type'], study['study_abstract'], '',
                                   'ENA_submission_data']) + '\n')  # assuming no pubmed_id
for sample_alias, sample in samples_dict.items():
    # sample_alias = sample_alias + '_' + timestamp
    if args.viral_submission:
        if sample['collector name'] == '':
            sample['collector name'] = 'unknown'
        if isinstance(sample['collection date'], float):
            # excel saved it as date
            year, month, day, hour, minute, second = xlrd.xldate_as_tuple(sample['collection date'],
                                                                          xl_workbook.datemode)
            month = "{:02d}".format(month)
            day = "{:02d}".format(day)
            hour = "{:02d}".format(hour)
            minute = "{:02d}".format(minute)
            second = "{:02d}".format(second)
            # format it as 2008-01-23T19:23:10
            collection_date_str = str(year) + '-' + str(month) + '-' + str(day) + 'T' + str(hour) + \
                ':' + str(minute) + ':' + str(second)
        else:
            collection_date_str = str(sample['collection date'])
        samples_table.write('\t'.join([sample_alias, action, 'ena_accession', sample['title'],
                                       sample['scientific_name'], 'tax_id_updated_by_ENA',
                                       sample['sample_description'], collection_date_str,
                                       sample['geographic location (country and/or sea)'],
                                       sample['host common name'], 'host subject id',
                                       sample['host health state'], sample['host sex'],
                                       sample['host scientific name'], sample['collector name'],
                                       sample['collecting institution'], sample['isolate'],
                                       'ENA_submission_date']) + '\n')
    else:
        samples_table.write('\t'.join([sample_alias, action, 'ena_accession', sample['title'],
                                       sample['scientific_name'], 'tax_id_updated_by_ENA',
                                       sample['sample_description']]) + '\n')
    for exp_alias, exp in experiments_dict.items():
        # should I check here if any experiment has a study or sample alias that is incorrect?
        # (not listed in the samples or study dict)
        # process the experiments for this sample
        if exp['sample_alias'] == sample_alias:
            lib_alias = 'library_' + exp_alias + '_' + exp['sample_alias']
            experiments_table.write('\t'.join([exp_alias, action, 'accession_ena', exp['title'],
                                               exp['study_alias'], sample_alias,
                                               exp['design_description'], lib_alias,
                                               exp['library_strategy'], exp['library_source'],
                                               exp['library_selection'],
                                               exp['library_layout'].lower(),
                                               str(int(exp['insert_size'])),
                                               exp['library_construction_protocol'],
                                               exp['platform'], exp['instrument_model'],
                                               'submission_date_ENA']) + '\n')
            exp_included.append(exp_alias)
            for run_alias, run in runs_dict.items():
                # check that the experiments library_layout is set to paired
                # when multiple entries are associated with the same run alias
                if not isinstance(run, list):
                    runs_list = [run]
                else:
                    runs_list = run
                for run_entry in runs_list:
                    if run_entry['experiment_alias'] == exp_alias:
                        runs_table.write('\t'.join([run_alias, action, 'ena_run_accession',
                                                    exp_alias, run_entry['file_name'],
                                                    FILE_FORMAT, 'file_checksum',
                                                    'submission_date_ENA']) + '\n')
                runs_included.append(run_alias)

# check if any experiment or run was not associated with any sample
for run in runs_dict.keys():
    if run not in runs_included:
        print(f'The run {run} is listed in the runs section but not associated with any \
              used experiment')

for exp in experiments_dict.keys():
    if exp not in exp_included:
        print(f'The experiment {exp} is listed in the experiments section but not associated \
              with any used sample')

studies_table.close()
samples_table.close()
experiments_table.close()
runs_table.close()

if args.verbose:
    paste_xls2yaml(args.xlsx_path)
