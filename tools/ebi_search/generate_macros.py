#!/usr/bin/env python3

import sys
import ebeye_urllib3

class WritableObject:
    def __init__(self):
        self.content = []
    def write(self, string):
        self.content.append(string)

def extract_domains():
    # redirection of sys.stdout
    domain_hierarchy = WritableObject()
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

    field_info = WritableObject()
    sys.stdout = field_info
    ebeye_urllib3.getDomainDetails(domain_id)
    sys.stdout = sys.__stdout__

    searchable_fields, retrievable_fields = extract_fiel_info(field_info.content)

    domain_info = {}
    domain_info["name"] = domain_name
    domain_info["searchable_fields"] = searchable_fields
    domain_info["retrievable_fields"] = retrievable_fields

    return domain_id, domain_info


def extract_fields(domains):
    domains_fields = {}

    for domain in domains:
        domain_id, domain_info = extract_a_domain_fields(domain)
        domains_fields.setdefault(domain_id, domain_info)

    return domains_fields


def generate_macros ():
    domains = extract_domains()
    domains_fields = extract_fields(domains)

if __name__ == '__main__':
    #parser = argparse.ArgumentParser()
    #parser.add_argument('', required=True)
    #args = parser.parse_args()

    generate_macros()
