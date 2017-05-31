#!/usr/bin/env python

import ebisearch

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


def write_domain_options(domains, domain_ids):
    """
    Write the domains as options
    """
    to_write = '%s<param argument="--domain" type="select" ' % (3 * spaces)
    to_write += 'label="Domain identifier in EBI">\n'
    for domain in domain_ids:
        to_write += '%s<option value="%s">%s</option>\n' % (
            4 * spaces,
            domain,
            format_name(domains[domain], domain))
    to_write += '%s</param>\n' % (3 * spaces)
    return to_write


def add_queries(searchable_fields):
    """
    Add (free text and constrained) queries on the searchable fields
    """
    to_write = ""
    if len(searchable_fields) == 0:
        return to_write
    to_write += '%s<conditional name="query">\n' % (4 * spaces)
    # Add the possibility to choose between a constrained or a free-text query
    to_write += '%s<param name="type" type="select" label="Type of query">\n' % (5 * spaces)
    to_write += '%s<option value="free">Free text query</option>\n' % (6 * spaces)
    to_write += '%s<option value="constrained">Constrained query</option>\n' % (6 * spaces)
    to_write += '%s</param>\n' % (5 * spaces)
    to_write += '%s<when value="free">\n' % (5 * spaces)
    to_write += '%s<param argument="--query" type="text" label="Query"/>\n' % (6 * spaces)
    to_write += '%s</when>\n' % (5 * spaces)
    to_write += '%s<when value="constrained">\n' % (5 * spaces)
    to_write += '%s<repeat name="queries" title="Add a query">\n' % (6 * spaces)
    # Parameter for the combination operator
    to_write += '%s<param name="combination_operation" type="select" label="Combination operation">\n' % (7 * spaces)
    to_write += '%s<option value="AND">AND</option>\n' % (8 * spaces)
    to_write += '%s<option value="OR">OR</option>\n' % (8 * spaces)
    to_write += '%s<option value="NOT">NOT</option>\n' % (8 * spaces)
    to_write += '%s</param>\n' % (7 * spaces)
    # Parameter with the searchable field
    to_write += '%s<param name="query_field" type="select" label="Fields">\n' % (7 * spaces)
    searchable_field_ids = sort_by_name(searchable_fields)
    for field in searchable_field_ids:
        to_write += '%s<option value="%s">%s</option>\n' % (
            8 * spaces,
            field,
            format_name(searchable_fields[field], field))
    to_write += '%s</param>\n' % (7 * spaces)
    # Parameter with the type of operations
    to_write += '%s<conditional name="comp_operation">\n' % (7 * spaces)
    to_write += '%s<param name="operation" type="select" label="Comparison operation">\n' % (8 * spaces)
    to_write += '%s<option value="equal">equal</option>\n' % (9 * spaces)
    to_write += '%s<option value="not">not</option>\n' % (9 * spaces)
    to_write += '%s<option value="range">range</option>\n' % (9 * spaces)
    to_write += '%s</param>\n' % (8 * spaces)
    to_write += '%s<when value="equal">\n' % (8 * spaces)
    to_write += '%s<param name="query_text" type="text" label="Searched term"/>\n' % (9 * spaces)
    to_write += '%s</when>\n' % (8 * spaces)
    to_write += '%s<when value="not">\n' % (8 * spaces)
    to_write += '%s<param name="query_text" type="text" label="Searched term"/>\n' % (9 * spaces)
    to_write += '%s<param name="not_query_text" type="text" label="Limiting term"/>\n' % (9 * spaces)
    to_write += '%s</when>\n' % (8 * spaces)
    to_write += '%s<when value="range">\n' % (8 * spaces)
    to_write += '%s<param name="min" type="text" label="From"/>\n' % (9 * spaces)
    to_write += '%s<param name="max" type="text" label="To"/>\n' % (9 * spaces)
    to_write += '%s</when>\n' % (8 * spaces)
    to_write += '%s</conditional>\n' % (7 * spaces)
    to_write += '%s</repeat>\n' % (6 * spaces)
    to_write += '%s</when>\n' % (5 * spaces)
    to_write += '%s</conditional>\n' % (4 * spaces)
    return to_write


def add_retrievable_field(retrievable_fields):
    """
    Add fields to extract
    """
    to_write = ''
    if len(retrievable_fields) == 0:
        return to_write
    retrievable_field_ids = sort_by_name(retrievable_fields)
    to_write += '%s<repeat name="fields" title="Add field to export">\n' % (4 * spaces)
    to_write += '%s<param argument="--field" type="select" label="Field to export">\n' % (5 * spaces)
    for field in retrievable_field_ids:
        to_write += '%s<option value="%s">%s</option>\n' % (
            6 * spaces,
            field,
            format_name(retrievable_fields[field], field))
    to_write += '%s</param>\n' % (5 * spaces)
    to_write += '%s</repeat>\n' % (4 * spaces)
    return to_write


def add_sorting(sortable_fields):
    """
    Add sorting possibility (simple or complex) if sortable fields
    """
    to_write = ''
    if len(sortable_fields) == 0:
        return to_write
    sortable_field_ids = sort_by_name(sortable_fields)
    to_write += '%s<conditional name="sort">\n' % (4 * spaces)
    to_write += '%s<param name="selection" type="select" label="Sorting?">\n' % (5 * spaces)
    to_write += '%s<option value="none">None</option>\n' % (6 * spaces)
    to_write += '%s<option value="simple">Simple</option>\n' % (6 * spaces)
    to_write += '%s<option value="complex">Complex</option>\n' % (6 * spaces)
    to_write += '%s</param>\n' % (5 * spaces)
    to_write += '%s<when value="none"/>\n' % (5 * spaces)
    # Simple sorting only on a field
    to_write += '%s<when value="simple">\n' % (5 * spaces)
    to_write += '%s<param argument="--sort_field" type="select" label="Field to sort">\n' % (6 * spaces)
    for field in sortable_field_ids:
        to_write += '%s<option value="%s">%s</option>\n' % (
            7 * spaces,
            field,
            format_name(sortable_fields[field], field))
    to_write += '%s</param>\n' % (6 * spaces)
    to_write += '%s<param argument="--order" type="select" label="Order to sort the results on the field">\n' % (6 * spaces)
    to_write += '%s<option value="ascending">Ascending</option>\n' % (7 * spaces)
    to_write += '%s<option value="descending">Descending</option>\n' % (7 * spaces)
    to_write += '%s</param>\n' % (6 * spaces)
    to_write += '%s</when>\n' % (5 * spaces)
    # Complex sorting only on multiple fields
    to_write += '%s<when value="complex">\n' % (5 * spaces)
    to_write += '%s<repeat name="complex_sort" title="Add a sorting criteria">\n' % (6 * spaces)
    to_write += '%s<param name="field" type="select" label="Field to sort">\n' % (7 * spaces)
    for field in sortable_field_ids:
        to_write += '%s<option value="%s">%s</option>\n' % (
            8 * spaces,
            field,
            format_name(sortable_fields[field], field))
    to_write += '%s</param>\n' % (7 * spaces)
    to_write += '%s<param name="order" type="select" label="Order to sort the results on the field">\n' % (7 * spaces)
    to_write += '%s<option value="ascending">Ascending</option>\n' % (8 * spaces)
    to_write += '%s<option value="descending">Descending</option>\n' % (8 * spaces)
    to_write += '%s</param>\n' % (7 * spaces)
    to_write += '%s</repeat>\n' % (6 * spaces)
    to_write += '%s</when>\n' % (5 * spaces)
    to_write += '%s</conditional>\n' % (4 * spaces)
    return to_write


def fill_when(domains, domain_ids):
    """
    Fill the conditional parameters for each domain
    """
    param_to_write = ""
    command_to_write = "#if "
    sep = ""
    for domain in domain_ids:
        fields = ebisearch.get_fields(domain, verbose=False)
        param_to_write += '%s<when value="%s">\n' % (3 * spaces, domain)
        param_to_write += add_queries(fields['searchable'])
        param_to_write += add_retrievable_field(fields['retrievable'])
        sorting_to_write = add_sorting(fields['sortable'])
        param_to_write += sorting_to_write
        param_to_write += '%s</when>\n' % (3 * spaces)
        if sorting_to_write != '':
            command_to_write += '%s$searched_domain.domain == "%s" ' % (
                sep,
                domain)
            sep = "or "

    command_to_write += '\n'
    command_to_write += '%s#if $searched_domain.sort.selection == "simple"\n' % (spaces)
    command_to_write += "%s--sort_field '$searched_domain.sort.sort_field'\n" % (2 * spaces)
    command_to_write += "%s--order '$searched_domain.sort.order'\n" % (2 * spaces)
    command_to_write += '%s#else if $searched_domain.sort.selection == "complex"\n' % (spaces)
    command_to_write += "%s#for $i, $s in enumerate( $searched_domain.sort.complex_sort )\n" % (2 * spaces)
    command_to_write += "%s--sort '$s.field:$s.order'\n" % (3 * spaces)
    command_to_write += "%s#end for\n" % (2 * spaces)
    command_to_write += "%s#end if\n" % (spaces)
    command_to_write += "#end if\n"
    return param_to_write, command_to_write


def generate_search_macros(filepath):
    """
    Generate the content of the macro file
    """
    domains = ebisearch.get_domains(verbose=False)
    domain_ids = sort_by_name(domains)
    to_write = '<?xml version="1.0" ?>\n'
    to_write += '<macros>\n'
    to_write += '%s<xml name="inputs">\n' % (spaces)
    to_write += '%s<conditional name="searched_domain">\n' % (2 * spaces)
    to_write += write_domain_options(domains, domain_ids)
    param_to_write, command_to_write = fill_when(domains, domain_ids)
    to_write += param_to_write
    to_write += '%s</conditional>\n' % (2 * spaces)
    to_write += '%s</xml>\n' % (spaces)
    to_write += '%s<token name="@SORTING@">\n' % (spaces)
    to_write += '<![CDATA[\n'
    to_write += '%s' % (command_to_write)
    to_write += ']]>\n'
    to_write += '%s</token>\n' % (spaces)
    to_write += '</macros>\n'
    with open(filepath, "w") as file:
        file.write(to_write)


if __name__ == '__main__':
    generate_search_macros("search_macros.xml")
