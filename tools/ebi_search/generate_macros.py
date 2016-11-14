#!/usr/bin/env python3

import ebeye_urllib3


def add_option(value, name, selected=False):
    to_write = "<option "
    to_write += "value=\"%s\"" % (value)
    if selected:
        to_write += " selected=\"true\""
    to_write += ">%s" % (name)
    to_write += "</option>\n"
    return to_write


def add_select_parameter(name, label, multiple=False):
    to_write = "<param "
    to_write += "name=\"%s\" " % (name)
    to_write += "type=\"select\" "
    to_write += "label=\"%s\" " % (label)
    if multiple:
        to_write += "multiple=\"true\""
    to_write += ">\n"
    return to_write


def write_macros_file(macros_filepath, domains_fields):
    spaces = "    "
    to_write = "<macros>\n"
    to_write += "%s<xml name=\"inputs\">\n" % (spaces)

    to_write += "%s<conditional name=\"searched_domain\">\n" % (2 * spaces)
    to_write += "%s%s" % (3 * spaces, add_select_parameter(
        "domain",
        "Domain to query"))

    sorted_domains = sorted(list(domains_fields.keys()))
    for domain in sorted_domains:
        to_write += "%s%s" % (4 * spaces, add_option(
            domain,
            domains_fields[domain]["name"]))

    to_write += "%s</param>\n\n" % (3 * spaces)

    for domain in sorted_domains:
        to_write += "%s<when value=\"%s\">\n" % (3 * spaces, domain)

        to_write += "%s%s" % (4 * spaces, add_select_parameter(
            "fields",
            "Fields to extract",
            multiple=True))
        for field in domains_fields[domain]["retrievable_fields"]:
            to_write += "%s%s" % (5 * spaces, add_option(
                field,
                field,
                selected=True))
        to_write += "%s<validator type=\"no_options\" message=\"Please select at least one field\" />\n" % (5 * spaces)
        to_write += "%s</param>\n" % (4 * spaces)

        to_write += "%s<repeat name=\"queries\" title=\"Add a query\">\n" % (
            4 * spaces)

        to_write += "%s%s" % (5 * spaces, add_select_parameter(
            "combination_operation",
            "Combination operation"))
        to_write += "%s%s" % (6 * spaces, add_option("AND", "AND"))
        to_write += "%s%s" % (6 * spaces, add_option("OR", "OR"))
        to_write += "%s%s" % (6 * spaces, add_option("NOT", "NOT"))
        to_write += "%s</param>\n" % (5 * spaces)

        to_write += "%s%s" % (5 * spaces, add_select_parameter(
            "query_field",
            "Fields"))
        for field in domains_fields[domain]["searchable_fields"]:
            to_write += "%s%s" % (6 * spaces, add_option(field, field))
        to_write += "%s</param>\n" % (5 * spaces)

        to_write += "%s<conditional name=\"comp_operation\">\n" % (5 * spaces)
        to_write += "%s%s" % (6 * spaces, add_select_parameter(
            "operation",
            "Comparison operation"))
        to_write += "%s%s" % (7 * spaces, add_option("equal", "equal"))
        to_write += "%s%s" % (7 * spaces, add_option("not", "not"))
        to_write += "%s%s" % (7 * spaces, add_option("range", "range"))
        to_write += "%s</param>\n" % (6 * spaces)

        to_write += "%s<when value=\"equal\">\n" % (6 * spaces)
        to_write += "%s<param name=\"query_text\" type=\"text\" label=\"Searched term\"/>\n" % (7 * spaces)
        to_write += "%s</when>\n" % (6 * spaces)

        to_write += "%s<when value=\"not\">\n" % (6 * spaces)
        to_write += "%s<param name=\"query_text\" type=\"text\" label=\"Searched term\"/>\n" % (7 * spaces)
        to_write += "%s<param name=\"not_query_text\" type=\"text\" label=\"Limiting term\"/>\n" % (7 * spaces)
        to_write += "%s</when>\n" % (6 * spaces)

        to_write += "%s<when value=\"range\">\n" % (6 * spaces)
        to_write += "%s<param name=\"min\" type=\"text\" label=\"From\"/>\n" % (7 * spaces)
        to_write += "%s<param name=\"max\" type=\"text\" label=\"To\"/>\n" % (
            7 * spaces)
        to_write += "%s</when>\n" % (6 * spaces)

        to_write += "%s</conditional>\n" % (5 * spaces)

        to_write += "%s</repeat>\n" % (4 * spaces)

        to_write += "%s</when>\n\n" % (3 * spaces)

    to_write += "%s</conditional>\n" % (2 * spaces)
    to_write += "%s</xml>\n" % (spaces)
    to_write += "</macros>\n"

    with open(macros_filepath, "w") as macros_file:
        macros_file.write(to_write)


def generate_macros():
    domains_fields = ebeye_urllib3.getDomainHierarchy()
    write_macros_file("macros.xml", domains_fields)


if __name__ == '__main__':
    generate_macros()
