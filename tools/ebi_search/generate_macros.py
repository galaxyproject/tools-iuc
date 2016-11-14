#!/usr/bin/env python3

import sys
import ebeye_urllib3
import writableObject


def add_option(value, name, selected="false"):
    to_write = "<option value=\"" + value + "\" selected=\"" + selected + "\">"
    to_write += name
    to_write += "</option>\n"
    return to_write


def add_select_parameter(name, label, multiple="false"):
    to_write = "<param name=\"" + name + "\" type=\"select\" "
    to_write += "label=\"" + label + "\" multiple=\"" + multiple + "\">\n"
    return to_write


def write_macros_file(macros_filepath, domains_fields):
    spaces = "    "
    to_write = "<macros>\n"
    to_write += spaces + "<xml name=\"inputs\">\n"

    to_write += 2*spaces + "<conditional name=\"searched_domain\">\n"
    to_write += 3*spaces + add_select_parameter("domain", "Domain to query")

    sorted_domains = sorted(list(domains_fields.keys()))
    for domain in sorted_domains:
        to_write += 4*spaces
        to_write += add_option(domain, domains_fields[domain]["name"])

    to_write += 3*spaces + "</param>\n\n"

    for domain in sorted_domains:
        to_write += 3*spaces + "<when value=\"" + domain + "\">\n"

        to_write += 4*spaces + add_select_parameter(
            "fields",
            "Fields to extract",
            multiple="true")
        for field in domains_fields[domain]["retrievable_fields"]:
            to_write += 5*spaces + add_option(field, field, selected="true")
        to_write += 5*spaces + "<validator type=\"no_options\" "
        to_write += "message=\"Please select at least one field\" />\n"
        to_write += 4*spaces + "</param>\n"

        to_write += 4*spaces + "<repeat name=\"queries\" "
        to_write += "title=\"Add a query\">\n"

        to_write += 5*spaces + add_select_parameter(
            "combination_operation",
            "Combination operation")
        to_write += 6*spaces + add_option("AND", "AND")
        to_write += 6*spaces + add_option("OR", "OR")
        to_write += 6*spaces + add_option("NOT", "NOT")
        to_write += 5*spaces + "</param>\n"

        to_write += 5*spaces + add_select_parameter(
            "query_field",
            "Fields")
        for field in domains_fields[domain]["searchable_fields"]:
            to_write += 6*spaces + add_option(field, field)
        to_write += 5*spaces + "</param>\n"

        to_write += 5*spaces + "<conditional name=\"comp_operation\">\n"
        to_write += 6*spaces + add_select_parameter(
            "operation",
            "Comparison operation")
        to_write += 7*spaces + add_option("equal", "equal")
        to_write += 7*spaces + add_option("not", "not")
        to_write += 7*spaces + add_option("range", "range")
        to_write += 6*spaces + "</param>\n"

        to_write += 6*spaces + "<when value=\"equal\">\n"
        to_write += 7*spaces + "<param name=\"query_text\" type=\"text\" "
        to_write += "label=\"Searched term\"/>\n"
        to_write += 6*spaces + "</when>\n"

        to_write += 6*spaces + "<when value=\"not\">\n"
        to_write += 7*spaces + "<param name=\"query_text\" type=\"text\" "
        to_write += "label=\"Searched term\"/>\n"
        to_write += 7*spaces + "<param name=\"not_query_text\" type=\"text\" "
        to_write += "label=\"Limiting term\"/>\n"
        to_write += 6*spaces + "</when>\n"

        to_write += 6*spaces + "<when value=\"range\">\n"
        to_write += 7*spaces + "<param name=\"min\" type=\"text\" "
        to_write += "label=\"From\"/>\n"
        to_write += 7*spaces + "<param name=\"max\" type=\"text\" "
        to_write += "label=\"To\"/>\n"
        to_write += 6*spaces + "</when>\n"

        to_write += 5*spaces + "</conditional>\n"

        to_write += 4*spaces + "</repeat>\n"

        to_write += 3*spaces + "</when>\n\n"

    to_write += 2*spaces + "</conditional>\n"
    to_write += spaces + "</xml>\n"
    to_write += "</macros>\n"

    with open(macros_filepath, "w") as macros_file:
        macros_file.write(to_write)


def generate_macros():
    #domains = extract_domains()
    #domains_fields = extract_fields(domains)
    domains_fields = ebeye_urllib3.getDomainHierarchy()
    write_macros_file("macros.xml", domains_fields)


if __name__ == '__main__':
    generate_macros()
