#!/usr/bin/env python

import enasearch

spaces = '    '


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
        s += '%s<option value="%s">%s</option>\n' % (2*spaces, f, f)
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
        s += '%s<option value="%s">%s</option>\n' % (2*spaces, opt, options[opt]['description'])
        when_s += '%s<when value="%s">\n' % (2*spaces, opt)
        if opt == 'fasta' or opt == 'fastq':
            when_s += '%s<param name="range_start" argument="--subseq_range" type="integer" optional="true" label="Start integer for subsequences"/>\n'  % (3*spaces)
            when_s += '%s<param name="range_stop" argument="--subseq_range" type="integer" optional="true" label="Stop integer for subsequences"/>\n'  % (3*spaces)
        else:
            when_s += '%s<param argument="--offset" type="integer" optional="true" label="First record to get"/>\n'  % (3*spaces)
            when_s += '%s<param argument="--length" type="integer" optional="true" label="Number of records to retrieve"/>\n'  % (3*spaces)
        when_s += '%s</when>\n'  % (2*spaces)
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
        s += '%s<option value="%s">%s</option>\n' % (2*spaces, f, f)
    s += '%s</xml>\n' % (spaces)
    return s


def write_taxonomy_results():
    """
    Write the possible taxonomy results
    """
    s = '%s<xml name="taxonomy_results">\n' % (spaces)
    fields = enasearch.get_taxonomy_results(verbose=False)
    for f in fields:
        s += '%s<option value="%s">%s</option>\n' % (2*spaces, f, fields[f]['description'])
    s += '%s</xml>\n' % (spaces)
    return s


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
    s += '</macros>\n'
    with open(filepath, "w") as file:
        file.write(s)


if __name__ == '__main__':
    generate_search_macros("search_macros.xml")
