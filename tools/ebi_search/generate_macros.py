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
        print(line)
        domain = line[:-1].split("\t")[-1]
        level = line[:-1].count("\t")

        if level <= last_level and tmp_domains != "":
            domains.append(tmp_domains)
        tmp_domains = domain
        last_level = level
        domains.append(tmp_domains)
    return domains


def extract_fields(domains):
    domains_fields = {}

    for domain in domains:
        domain_id = domain.split(":")[0]
        domain_name = domain.split(":")[1]

        domains_fields.setdefault(domain_name, {})
        domains_fields[domain_name].setdefault("name", domain_name)
        domains_fields[domain_name].setdefault("searchable_fields", [])
        domains_fields[domain_name].setdefault("retrievable_fields", [])



    return domains_fields


def generate_macros ():
    domains = extract_domains()
    domains_fields = extract_fields(domains)

if __name__ == '__main__':
    #parser = argparse.ArgumentParser()
    #parser.add_argument('', required=True)
    #args = parser.parse_args()

    generate_macros()
