#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Extract information from the Nextclade dataset list command that is
# needed to generate the Galaxy list of possible datasets and
# their columns names.

import argparse
import json
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile

from typing import TextIO

def get_nextclade_version() -> str:
    """
    Get the version of Nextclade installed on the system.
    Returns the version as a string.
    """

    nextclade_command = ["nextclade", "--version"]
    proc = subprocess.run(
        nextclade_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )

    if proc.returncode != 0:
        raise RuntimeError(f"Nextclade version check failed: {proc.stderr.decode()}")

    # Nextclade outputs the version in the format "nextclade X.Y.Z"
    version = proc.stdout.decode().strip().split()[-1]
    return version


def get_datasets() -> list[dict]:
    """
    Get the list of datasets from the Nextclade dataset list command.
    Returns a list of dictionaries, each representing a dataset.
    """

    nextclade_command = ["nextclade", "dataset", "list", "--json"]
    proc = subprocess.run(
        nextclade_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )

    if proc.returncode != 0:
        raise RuntimeError(f"Nextclade dataset list failed: {proc.stderr.decode()}")

    datasets_json = proc.stdout.decode()
    datasets = json.loads(datasets_json)

    return datasets


def fetch_reference(dataset_name: str) -> str:
    """
    Using nextclade dataset get command,
    fetch the reference file for the given dataset name.
    """

    tempdir = tempfile.TemporaryDirectory(delete=False)
    nextclade_command = [
        "nextclade",
        "dataset",
        "get",
        "--name",
        dataset_name,
        "--output-dir",
        tempdir.name,
    ]

    proc = subprocess.run(
        nextclade_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"Nextclade dataset get failed for dataset {dataset_name}: {proc.stderr.decode()}"
        )

    # The reference file is expected to be in the tempdir with the name "reference.fasta"
    tempdir_path = pathlib.Path(tempdir.name)
    dataset_info = json.load(open(tempdir_path / "pathogen.json"))
    assert (
        "files" in dataset_info
    ), f"Dataset {dataset_name} does not have 'files' key in pathogen.json"
    assert (
        "reference" in dataset_info["files"]
    ), f"Dataset {dataset_name} does not have 'reference' file in pathogen.json"
    reference_file = tempdir_path / dataset_info["files"]["reference"]
    assert (
        reference_file.exists()
    ), f"Reference file {reference_file} does not exist for dataset {dataset_name}"
    return str(reference_file)


def run_nextclade(dataset_name: str, reference_filename: str) -> list[str]:
    """
    Run the Nextclade command to process the dataset with the given name and reference file.
    Returns a list of columns in the Nextclade output.
    """

    output_filename = tempfile.NamedTemporaryFile(delete=False, suffix=".tsv")
    output_filename.close()
    dataset_dir = pathlib.Path(reference_filename).parent
    nextclade_command = [
        "nextclade",
        "run",
        "--input-dataset",
        dataset_dir,
        "--output-tsv",
        output_filename.name,
        reference_filename,
    ]
    proc = subprocess.run(
        nextclade_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )

    if proc.returncode != 0:
        raise RuntimeError(
            f"Nextclade run failed for dataset {dataset_name}: {proc.stderr.decode()}"
        )

    header = next(open(output_filename.name)).strip().split("\t")
    os.remove(output_filename.name)  # Clean up the temporary file
    tempdir = pathlib.Path(reference_filename).parent
    shutil.rmtree(
        tempdir, ignore_errors=True
    )  # Clean up the temporary dataset directory
    return header


def write_with_indent(file: TextIO, text: str, indent: int) -> None:
    """
    Write text to the file with the specified indentation.
    """
    file.write(" " * indent + text + "\n")


def generate_macros(
    macro_file: TextIO, nextclade_version: str, dataset_info: list[dict]
) -> None:
    """
    Generate the Galaxy metadata format for the given datasets and their columns.
    """

    indent_size = 4
    indent = 0
    indent += indent_size
    write_with_indent(macro_file, "<macros>", indent)
    indent += indent_size
    write_with_indent(
        macro_file, f'<token name="@TOOL_VERSION@">{nextclade_version}</token>', indent
    )
    write_with_indent(macro_file, '<xml name="dataset_selector">', indent)
    seen_first = False
    names_seen = set()
    descriptions_seen = set()
    column_info = {}
    for dataset in dataset_info:
        name = dataset["name"]
        if name == "sars-cov-2":
            column_info[name] = dataset["columns"]
        elif name == "MPXV":
            column_info[name] = dataset["columns"]
        description = "<![CDATA[" + dataset["description"] + "]]>"
        if name in names_seen:
            print("Warning: Duplicate name found:", name)
        names_seen.add(name)
        if description in descriptions_seen:
            print("Warning: Duplicate description found:", description)
        descriptions_seen.add(description)

        if not seen_first:
            selected = ' selected="true"'
            seen_first = True
        else:
            selected = ""
        dataset_selector = f'<option value="{name}"{selected}>{description}</option>'
        write_with_indent(macro_file, dataset_selector, indent)
    indent -= indent_size
    write_with_indent(macro_file, "</xml>", indent)
    write_with_indent(macro_file, '<xml name="output_columns">', indent)
    indent += indent_size
    write_with_indent(macro_file, '<conditional name="organism">', indent)
    indent += indent_size
    for dataset in dataset_info:
        write_with_indent(macro_file, f'<when value="{dataset["name"]}">', indent)
        indent += indent_size
        columns = dataset["columns"]
        columns_str = ",".join(columns)
        action = (
            f'<action name="column_names" type="metadata" default="{columns_str}"/>'
        )
        write_with_indent(macro_file, action, indent)
        indent -= indent_size
        write_with_indent(macro_file, "</when>", indent)
    indent -= indent_size
    write_with_indent(macro_file, "</conditional>", indent)
    indent -= indent_size
    write_with_indent(macro_file, "</xml>", indent)
    for name, columns in column_info.items():
        token_name = name.upper().replace("-", "_") + "_NUM_COLUMNS"
        write_with_indent(
            macro_file, f'<token name="@{token_name}@">{len(columns)}</token>', indent
        )
        token_name = (
            name.upper().replace("-", "_") + "_COLUMNS"
        )  # e.g. SARS_COV_2_COLUMNS
        columns_str = ",".join(columns)
        write_with_indent(
            macro_file, f'<token name="@{token_name}@">{columns_str}</token>', indent
        )
    indent -= indent_size
    write_with_indent(macro_file, "</macros>", indent)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert Nextclade dataset list JSON to Galaxy metadata format."
    )
    parser.add_argument(
        "output_file",
        type=argparse.FileType("w"),
        default=sys.stdout,
        nargs="?",
        help="Path to the output file for Galaxy metadata format.",
    )

    args = parser.parse_args()

    datasets = get_datasets()

    # While the datasets JSON has no defined schema, as of 3.15.3, it contains:
    # "path" - the long-form dataset name
    # "shortcuts" - a list of short-form dataset names
    # "attributes" - a dictionary of attributes, including the name and reference accession
    #  - "name" - a human-readable name for the dataset
    dataset_name_key = "path"
    dataset_aliases_key = "shortcuts"
    dataset_attributes_key = "attributes"
    dataset_attributes_name_key = "name"
    dataset_infos = []
    for dataset in datasets:
        name = dataset[dataset_name_key]
        if dataset_aliases_key in dataset:
            # Choose the first alias - typically the simplest "common name" - as the name
            name = dataset[dataset_aliases_key][0]
        assert (
            dataset_attributes_key in dataset
        ), f"Dataset {name} does not have attributes"
        assert (
            dataset_attributes_name_key in dataset[dataset_attributes_key]
        ), f"Dataset has attributes key but attributes does not have {dataset_attributes_name_key}"
        # since we're going to use this in a macro, we XML-escape the description
        # to avoid issues with special characters e.g. "&" or "<"
        description = dataset[dataset_attributes_key][dataset_attributes_name_key]
        # special case for the two influenza datasets where a description is duplicated
        if (
            description == "Influenza A H1N1pdm HA"
            and "CY121680" in dataset[dataset_name_key]
        ):
            description += " (broad)"
        elif (
            description == "Influenza A H3N2 HA"
            and "CY163680" in dataset[dataset_name_key]
        ):
            description += " (broad)"
        if name.startswith("community/"):
            description = description + " (community contributed)"
        reference_filename = fetch_reference(name)
        columns = run_nextclade(name, reference_filename)
        dataset_infos.append(
            {"name": name, "description": description, "columns": columns}
        )

    nextclade_version = get_nextclade_version()
    generate_macros(args.output_file, nextclade_version, dataset_infos)
