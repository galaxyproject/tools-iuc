"""Run ``bioheaders metadata`` from the parameters of the FAIR-bioHeaders metadata tool.

Usage: fair_bioheaders_metadata.py INPUTS_JSON OUTPUT

Galaxy writes the tool's parameters to INPUTS_JSON (``<configfiles><inputs
data_style="paths"/>``). The text values are read from there, unsanitized, and
passed to bioheaders as an argument list, without a shell, so every character
reaches bioheaders unchanged. Each value is joined to its option with ``=``, so
a value starting with ``-`` is never read as an option. The command is printed
to standard error, and the exit status is that of bioheaders.
"""

import json
import shlex
import subprocess
import sys
import unicodedata

# U+2028 and U+2029, which bioheaders refuses with the control characters
LINE_BREAKS = "\u2028\u2029"


class InputError(ValueError):
    pass


def text(value, label):
    """``value`` as a string, refusing what bioheaders refuses: control characters
    other than tab, and line breaks. Galaxy's validators check this too, but they
    may not run for a value connected in a workflow."""
    value = str(value)
    for char in value:
        category = unicodedata.category(char)
        if (
            (category == "Cc" and char != "\t")
            or char in LINE_BREAKS
            or category == "Cs"
        ):
            raise InputError(
                f"{label}: line breaks and other control characters are not allowed ({char!r})"
            )
    return value


def given(value):
    """Whether an optional value was given (Galaxy writes null or "" when it wasn't)."""
    return value is not None and value != ""


def author(item, label):
    name = text(item["name"], label)
    if given(item.get("orcid")):
        return f"{name} <{text(item['orcid'], label + ' ORCID')}>"
    return name


def software(item):
    name = text(item["name"], "Assembly software name")
    version = item.get("software_version")
    version = text(version, "Assembly software version") if given(version) else ""
    if "=" in name or "=" in version:
        raise InputError("Assembly software: the name and version cannot contain =")
    uri = item.get("uri")
    if given(uri):
        return f"{name}={version}={text(uri, 'Assembly software URI')}"
    if version:
        return f"{name}={version}"
    return name


# Optional single values in the "optional_fields" section, in command order
OPTIONAL = [
    ("assembly_protocol", "--assembly-protocol"),
    ("documentation", "--documentation"),
    ("accession", "--accession"),
    ("accession_url", "--accession-url"),
]
OPTIONAL_AFTER_IDENTIFIERS = [("reuse_conditions", "--reuse-conditions")]
OPTIONAL_AFTER_LINKS = [
    ("voucher_specimen", "--voucher-specimen"),
    ("scholarly_article", "--scholarly-article"),
    ("funding", "--funding"),
    ("seqcol_id", "--seqcol-id"),
]


def build_command(params, output):
    """The ``bioheaders metadata`` argument list for the Galaxy parameters ``params``."""
    command = ["bioheaders", "metadata"]

    def add(option, value, label):
        command.append(f"{option}={text(value, label)}")

    def add_optional(section, fields):
        for name, option in fields:
            if given(section.get(name)):
                add(option, section[name], option)

    optional = params.get("optional_fields") or {}
    add("--genome", params["genome"], "Genome name")
    for item in params.get("genome_synonyms") or []:
        add("--genome-synonym", item["value"], "Genome synonym")
    add("--taxon-name", params["taxon_name"], "Taxon name")
    add("--taxon-id", params["taxon_id"], "NCBI taxonomy ID")
    add("--version", params["assembly_version"], "Assembly version")
    for item in params.get("metadata_authors") or []:
        command.append("--metadata-author=" + author(item, "Metadata author"))
    for item in params.get("assembly_authors") or []:
        command.append("--assembly-author=" + author(item, "Assembly author"))
    add("--date-created", params.get("date_created") or "", "Date created")
    add("--masking", params.get("masking") or "", "Masking")
    for item in optional.get("instruments") or []:
        add("--instrument", item["value"], "Instrument")
    for item in optional.get("assembly_software") or []:
        command.append("--assembly-software=" + software(item))
    add_optional(optional, OPTIONAL)
    for item in optional.get("identifiers") or []:
        add("--identifier", item["value"], "Identifier")
    add_optional(optional, OPTIONAL_AFTER_IDENTIFIERS)
    for item in optional.get("related_links") or []:
        add("--related-link", item["value"], "Related link")
    add_optional(optional, OPTIONAL_AFTER_LINKS)
    statistics = params.get("statistics") or {}
    if statistics.get("compute") == "yes":
        add("--stats-from", statistics["fasta"], "Assembly FASTA")
    for item in optional.get("stats") or []:
        add(
            "--stat",
            f"{item['key']}={text(item['value'], 'Statistic value')}",
            "Statistic",
        )
    for item in optional.get("extras") or []:
        add(
            "--extra",
            f"{item['key']}={text(item['value'], 'Extension field value')}",
            "Extension field",
        )
    add("--format", params.get("output_format") or "yaml", "Output format")
    command.append("--output=" + output)
    return command


def main(argv):
    if len(argv) != 3:
        print(f"usage: {argv[0]} INPUTS_JSON OUTPUT", file=sys.stderr)
        return 2
    with open(argv[1], encoding="utf-8") as handle:
        params = json.load(handle)
    try:
        command = build_command(params, argv[2])
    except InputError as error:
        print(error, file=sys.stderr)
        return 1
    print(shlex.join(command), file=sys.stderr, flush=True)
    return subprocess.run(command).returncode


if __name__ == "__main__":
    sys.exit(main(sys.argv))
