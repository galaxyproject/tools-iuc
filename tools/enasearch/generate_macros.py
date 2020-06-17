#!/usr/bin/env python

import enasearch

spaces = '    '
operator_names = {
    "=": "equal",
    "!=": "different",
    "<": "lower",
    "<=": "equal or lower",
    ">": "higher",
    ">=": "equal or higher",
}


def format_name(name, alternative_name):
    """
    Format name to remove None name and & in name
    """
    if name is None:
        name = alternative_name
    name = name.replace("&", "and")
    return name


def sort_by_name(dict):
    """
    Sort a dictionary on the values
    """
    return sorted(dict, key=dict.get)


def write_analysis_fields():
    """
    Write the analysis fields
    """
    s = '%s<xml name="analysis_fields">\n' % (spaces)
    fields = enasearch.get_returnable_fields(result="analysis", verbose=False)
    for f in fields:
        s += '%s<option value="%s">%s</option>\n' % (2 * spaces, f, f)
    s += '%s</xml>\n' % (spaces)
    return s


def write_display_options():
    """
    Write the display options
    """
    s = '%s<xml name="display_options">\n' % (spaces)
    when_s = '%s<xml name="when_display_options">\n' % (spaces)
    options = enasearch.get_display_options(verbose=False)
    for opt in options:
        s += '%s<option value="%s">%s</option>\n' % (2 * spaces, opt, options[opt]['description'])
        when_s += '%s<when value="%s">\n' % (2 * spaces, opt)
        if opt == 'fasta' or opt == 'fastq':
            when_s += '%s<param name="range_start" argument="--subseq_range" type="integer" optional="true" label="Start integer for subsequences"/>\n' % (3 * spaces)
            when_s += '%s<param name="range_stop" argument="--subseq_range" type="integer" optional="true" label="Stop integer for subsequences"/>\n' % (3 * spaces)
        else:
            when_s += '%s<param argument="--offset" type="integer" optional="true" label="First record to get"/>\n' % (3 * spaces)
            when_s += '%s<param argument="--length" type="integer" optional="true" label="Number of records to retrieve"/>\n' % (3 * spaces)
        when_s += '%s</when>\n' % (2 * spaces)
    s += '%s</xml>\n' % (spaces)
    when_s += '%s</xml>\n' % (spaces)
    s += when_s
    return s


def write_run_fields():
    """
    Write the run fields
    """
    s = '%s<xml name="run_fields">\n' % (spaces)
    fields = enasearch.get_returnable_fields(result="read_run", verbose=False)
    for f in fields:
        s += '%s<option value="%s">%s</option>\n' % (2 * spaces, f, f)
    s += '%s</xml>\n' % (spaces)
    return s


def write_taxonomy_results():
    """
    Write the possible taxonomy results
    """
    s = '%s<xml name="taxonomy_results">\n' % (spaces)
    fields = enasearch.get_taxonomy_results(verbose=False)
    for f in fields:
        s += '%s<option value="%s">%s</option>\n' % (2 * spaces, f, fields[f]['description'])
    s += '%s</xml>\n' % (spaces)
    return s


def write_result_parameters(fts=False):
    """
    Write the parameters that are dependant of results
    """
    res = enasearch.get_results(verbose=False)
    options = enasearch.get_display_options(verbose=False)
    ft = enasearch.get_filter_types(verbose=False)
    # Format the filter type related parameters
    ft_parameters = {}
    for t in ft:
        s = ''
        if 'operators' in ft[t]:
            s = '%s<param name="operation" type="select" label="Operator">\n' % (7 * spaces)
            for o in ft[t]['operators']:
                on = o
                if o in operator_names:
                    on = operator_names[o]
                s += '%s<option value="%s">%s</option>\n' % (8 * spaces, on, on)
            s += '%s</param>\n' % (7 * spaces)
            if 'value' in ft[t]:
                value_format = 'float' if t == 'Number' else 'text'
                s += '%s<param name="value" type="%s" value="" label="%s"/>\n' % (7 * spaces, value_format, ft[t]['value'])
            elif 'values' in ft[t]:
                s += '%s<param name="value" type="select" label="Value">\n' % (7 * spaces)
                for v in ft[t]['values']:
                    s += '%s<option value="%s">%s</option>\n' % (8 * spaces, v, v)
                s += '%s</param>\n' % (7 * spaces)
        else:
            s += '%s<conditional name="op">\n' % (7 * spaces)
            s += '%s<param name="operation" type="select" label="Operation">\n' % (8 * spaces)
            for op in ft[t]:
                s += '%s<option value="%s">%s</option>\n' % (9 * spaces, op, ft[t][op]['description'])
            s += '%s</param>\n' % (8 * spaces)
            for op in ft[t]:
                s += '%s<when value="%s">\n' % (8 * spaces, op)
                s += '%s<param name="values" type="text" value="" label="%s" help="Values separated by simple comma"/>\n' % (9 * spaces, ",".join(ft[t][op]['parameters']))
                s += '%s</when>\n' % (8 * spaces)
            s += '%s</conditional>\n' % (7 * spaces)
        ft_parameters[t] = s
    # Start adding the conditional
    s = '%s<conditional name="res">\n' % (2 * spaces)
    # Add result parameter
    s += '%s<param argument="--result" type="select" label="Result to return">\n' % (3 * spaces)
    for r in res:
        s += '%s<option value="%s">%s</option>\n' % (4 * spaces, r, res[r]['description'])
    s += '%s</param>\n' % (3 * spaces)
    for r in res:
        sf = enasearch.get_sortable_fields(r)
        ff = res[r]['filter_fields']
        s += '%s<when value="%s">\n' % (3 * spaces, r)
        if not fts:
            s += '%s<repeat name="queries" title="Add a query">\n' % (4 * spaces)
            # Add combination operator
            s += '%s<param name="combination_operation" type="select" label="Combination operation">\n' % (5 * spaces)
            s += '%s<option value="AND">AND</option>\n' % (6 * spaces)
            s += '%s<option value="OR">OR</option>\n' % (6 * spaces)
            s += '%s<option value="NOT">NOT</option>\n' % (6 * spaces)
            s += '%s</param>\n' % (5 * spaces)
            s += '%s<conditional name="filter_field">\n' % (5 * spaces)
            s += '%s<param name="field" type="select" label="Field to query">\n' % (6 * spaces)
            for f in ff:
                s += '%s<option value="%s">%s</option>\n' % (7 * spaces, f, ff[f]['description'])
            s += '%s</param>\n' % (6 * spaces)
            for f in ff:
                # Add the correct parameter given the type of field
                typ = ff[f]['type'].capitalize()
                if typ not in ft_parameters:
                    if f == 'location':
                        typ = 'Geospatial'
                    else:
                        continue
                s += '%s<when value="%s">\n' % (6 * spaces, f)
                s += ft_parameters[typ]
                s += '%s</when>\n' % (6 * spaces)
            s += '%s</conditional>\n' % (5 * spaces)
            s += '%s</repeat>\n' % (4 * spaces)
        # Add display opt
        s += '%s<conditional name="display_opt">\n' % (4 * spaces)
        s += '%s<param argument="--display" type="select" label="Display option to specify the display format">\n' % (5 * spaces)
        s += '%s<expand macro="display_options"/>\n' % (6 * spaces)
        s += '%s</param>\n' % (5 * spaces)
        for opt in options:
            s += '%s<when value="%s"' % (5 * spaces, opt)
            if opt != 'fasta' and opt != 'fastq':
                s += '>\n'
                s += '%s<param argument="--offset" type="integer" optional="true" label="First record to get"/>\n' % (6 * spaces)
                s += '%s<param argument="--length" type="integer" optional="true" label="Number of records to retrieve"/>\n' % (6 * spaces)
                if opt == 'report':
                    s += '%s<param argument="--fields" type="select" multiple="true" label="Fields to return">\n' % (6 * spaces)
                    for f in res[r]['returnable_fields']:
                        s += '%s<option value="%s">%s</option>\n' % (7 * spaces, f, f)
                    s += '%s</param>\n' % (6 * spaces)
                    s += '%s<param argument="--sortfields" type="select" optional="true" multiple="true" label="Fields to sort the results">\n' % (6 * spaces)
                    for f in sf:
                        s += '%s<option value="%s">%s</option>\n' % (7 * spaces, f, sf[f]['description'])
                    s += '%s</param>\n' % (6 * spaces)
                s += '%s</when>\n' % (5 * spaces)
            else:
                s += '/>\n'
        s += '%s</conditional>\n' % (4 * spaces)
        s += '%s</when>\n' % (3 * spaces)
    s += '%s</conditional>\n' % (2 * spaces)
    return s


def write_search_data_parameters():
    """
    Write the parameters for search_data
    """
    fts = '%s<xml name="free_text_search">\n' % (spaces)
    fts += write_result_parameters(True)
    fts += '%s</xml>\n' % (spaces)
    cts = '%s<xml name="conditional_text_search">\n' % (spaces)
    cts += write_result_parameters(False)
    cts += '%s</xml>\n' % (spaces)
    return fts + cts


def generate_search_macros(filepath):
    """
    Generate the content of the macro file
    """
    s = '<?xml version="1.0" ?>\n'
    s += '<macros>\n'
    s += write_analysis_fields()
    s += write_display_options()
    s += write_run_fields()
    s += write_taxonomy_results()
    s += write_search_data_parameters()
    s += '</macros>\n'
    with open(filepath, "w") as file:
        file.write(s)


if __name__ == '__main__':
    generate_search_macros("search_macros.xml")
