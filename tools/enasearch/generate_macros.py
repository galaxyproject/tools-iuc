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


def write_analysis_fields(filepath):
    """
    Write the analysis fields in a file
    """
    s = '%s<xml name="analysis_fields">\n' % (spaces)
    fields = enasearch.get_returnable_fields(result="analysis", verbose=False)
    for f in fields:
        s += '%s<option value="%s">%s</option>\n' % (2*spaces, f, f)
    s += '%s</xml>\n' % (spaces)


def generate_search_macros(filepath):
    """
    Generate the content of the macro file
    """
    s = '<?xml version="1.0" ?>\n'
    s += '<macros>\n'
    s += write_analysis_fields("tool-data/analysis_fields.loc")
    s += '</macros>\n'
    with open(filepath, "w") as file:
        file.write(s)


if __name__ == '__main__':
    generate_search_macros("search_macros.xml")
