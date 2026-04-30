#!/usr/bin/env python3


import argparse
import sys
import tempfile
import xml.etree.ElementTree as ET
from csv import DictWriter
from pathlib import Path
from shutil import which
from subprocess import PIPE, Popen

from yaml import safe_load

_TIBERIUS_REPO = "https://github.com/Gaius-Augustus/Tiberius"
_YAML_KEYS = ["target_species", "tiberius_version", "date", "softmasking", "clamsa"]
_LOC_FILE_PATH = Path("tool-data/tiberius_models.loc.sample")


def get_model_dict(yaml_file):
    with open(yaml_file, "rt") as f:
        model_config = safe_load(f)
        model_dict = {k: model_config.get(k) for k in _YAML_KEYS}
        model_dict["model_cfg"] = yaml_file.name
    return model_dict


def get_version(macros_xml):
    tree = ET.parse(macros_xml)
    for token in tree.findall("token"):
        if token.attrib["name"] == "@VERSION@":
            return str(token.text)
    raise ValueError("No token element with name @VERSION@ found")


def pull_repo(tool_version):

    git = which("git")

    tempdir = tempfile.mkdtemp()

    # pull the tag
    pull_command = [
        Path(git),
        "clone",
        "--single-branch",
        "--branch",
        f"v{tool_version}",
        _TIBERIUS_REPO,
        tempdir,
    ]
    with Popen(pull_command, stdout=PIPE) as proc:
        print(proc.stdout.read(), sys.stdout)

    return tempdir


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("macros_xml", help="macros.xml file", type=Path)

    return parser.parse_args()


def make_display_name(model_dict):
    base_name = (
        model_dict["target_species"] + f' version {model_dict["tiberius_version"]}'
    )
    model_type = [
        "Soft-masked" if model_dict["softmasking"] else None,
        "ClaMSA" if model_dict["clamsa"] else None,
    ]
    model_type = [x for x in model_type if x is not None]

    if model_type:
        return base_name + f' ({", ".join(model_type)})'

    return base_name


def main():
    args = parse_arguments()
    tool_version = get_version(args.macros_xml)

    tool_repo = pull_repo(tool_version)

    yaml_files = Path(tool_repo, "model_cfg").glob("*.yaml")

    model_dicts = [get_model_dict(x) for x in yaml_files]

    for model_dict in model_dicts:
        model_dict["value"] = make_display_name(model_dict)
        model_dict["container_version"] = tool_version

    with open(_LOC_FILE_PATH, "wt") as f:
        loc_file_keys = sorted(set().union(*(d.keys() for d in model_dicts)))
        header = "\t".join(loc_file_keys)
        f.write(f"#{header}\n")
        fc = DictWriter(f, fieldnames=loc_file_keys, delimiter="\t")
        fc.writerows(model_dicts)


if __name__ == "__main__":

    main()
