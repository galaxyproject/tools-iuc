#!/usr/bin/env python3

import sys
import ebeye_urllib3
import writableObject

def extract_domains():
    # redirection of sys.stdout
    domain_hierarchy = writableObject.WritableObject()
    sys.stdout = domain_hierarchy
    ebeye_urllib3.getDomainHierarchy()
    sys.stdout = sys.__stdout__

    domains = []
    last_level = 1
    tmp_domains = ""
    for line in domain_hierarchy.content:
        if line == "\n":
            continue
        domain = line.split("\t")[-1]
        level = line.count("\t")

        if level <= last_level and tmp_domains != "":
            domains.append(tmp_domains)
        tmp_domains = domain
        last_level = level
    domains.append(tmp_domains)
    return domains


def extract_fiel_info(field_info):
    searchable_fields = []
    retrievable_fields = []

    for line in field_info:
        if line == "\n" :
            continue

        if line.find("\t") == -1:
            continue

        if line.startswith("field_id") or line.startswith("$facets"):
            continue

        split_line = line.split("\t")
        field_id = split_line[0]
        searchable = split_line[1]
        retrievable = split_line[2]

        if searchable == "true":
            searchable_fields.append(field_id)

        if retrievable == "true":
            retrievable_fields.append(field_id)

    return searchable_fields, retrievable_fields


def extract_a_domain_fields(domain):
    domain_id = domain.split(":")[0]
    domain_name = domain.split(":")[1]

    searchable_fields = []
    retrievable_fields = []

    field_info = writableObject.WritableObject()
    sys.stdout = field_info
    ebeye_urllib3.getDomainDetails(domain_id)
    sys.stdout = sys.__stdout__

    searchable_fields, retrievable_fields = extract_fiel_info(field_info.content)

    domain_info = {}
    domain_info["name"] = domain_name
    domain_info["searchable_fields"] = sorted(searchable_fields)
    domain_info["retrievable_fields"] = sorted(retrievable_fields)

    return domain_id, domain_info


def extract_fields(domains):
    domains_fields = {}

    for domain in domains:
        domain_id, domain_info = extract_a_domain_fields(domain)
        domains_fields.setdefault(domain_id, domain_info)

    return domains_fields

def add_option(value, name):
    to_write = "<option value=\"" + value + "\">"
    to_write += name
    to_write += "</option>\n"
    return to_write

def add_select_parameter(name, label, multiple = "false"):
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
        to_write += 4*spaces + add_option(domain,
        domains_fields[domain]["name"])

    to_write += 3*spaces + "</param>\n\n"

    for domain in sorted_domains:
        to_write += 3*spaces + "<when value=\"" + domain + "\">\n"

        to_write += 4*spaces + add_select_parameter("fields",
        "Fields to extract", multiple = "true")
        for field in domains_fields[domain]["retrievable_fields"]:
            to_write += 5*spaces + add_option(field, field)

        to_write += 4*spaces + "</param>\n"

        to_write += 4*spaces + "<repeat name=\"queries\" "
        to_write += "title=\"Add a query\">\n"

        to_write += 5*spaces + add_select_parameter("combination_operation",
        "Combination operation")
        to_write += 6*spaces + add_option("AND", "AND")
        to_write += 6*spaces + add_option("OR", "OR")
        to_write += 6*spaces + add_option("NOT", "NOT")
        to_write += 5*spaces + "</param>\n"

        to_write += 5*spaces + add_select_parameter("query_field",
        "Fields")
        for field in domains_fields[domain]["searchable_fields"]:
            to_write += 6*spaces + add_option(field, field)
        to_write += 5*spaces + "</param>\n"

        to_write += 5*spaces + "<conditional name=\"comparison_operation\">\n"
        to_write += 6*spaces + add_select_parameter("comparison_operation",
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



def generate_macros ():
    domains = extract_domains()
    domains_fields = extract_fields(domains)
    write_macros_file("macros.xml", domains_fields)

if __name__ == '__main__':
    #parser = argparse.ArgumentParser()
    #parser.add_argument('', required=True)
    #args = parser.parse_args()

    generate_macros()
