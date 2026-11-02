#!/usr/bin/env python3
import argparse
import ast
import configparser
import json
import os
import re
import subprocess
from pathlib import Path


def load_json(path):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def clean_name(value, fallback):
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", (value or fallback).strip()).strip("._")
    return value or fallback


def output_designation(filename):
    for suffix in (".csv.gz", ".tsv.gz", ".txt.gz", ".text.gz", ".csv.bz2", ".tsv.bz2", ".txt.bz2", ".text.bz2"):
        if filename.endswith(suffix):
            return filename[: -len(suffix)]
    return Path(filename).stem


def unique_output_name(filename, used_designations, fallback):
    filename = clean_name(Path(str(filename).replace("\\", "/")).name, fallback)
    designation = output_designation(filename)
    if designation not in used_designations:
        used_designations.add(designation)
        return filename

    path = Path(filename)
    suffixes = "".join(path.suffixes)
    stem = filename[: -len(suffixes)] if suffixes else filename
    i = 2
    while f"{designation}_{i}" in used_designations:
        i += 1
    used_designations.add(f"{designation}_{i}")
    return f"{stem}_{i}{suffixes}"


def input_overrides(mappings, input_dir):
    input_dir.mkdir(exist_ok=True)
    overrides = []
    used = set()
    file_indexes = {}

    for i, mapping in enumerate(mappings, 1):
        node = mapping["node"].strip()
        parameter = mapping.get("parameter", "").strip()
        if parameter:
            match = re.fullmatch(r"files\.(\d+)\.filename", parameter)
            if match:
                file_indexes[node] = max(file_indexes.get(node, 0), int(match.group(1)) + 1)
        else:
            file_index = file_indexes.get(node, 0)
            parameter = f"files.{file_index}.filename"
            file_indexes[node] = file_index + 1

        ext = mapping.get("extension") or "data"
        identifier = mapping.get("identifier") or f"input_{i}.{ext}"
        staged_name = mapping.get("staged_name") or identifier
        if "." not in staged_name:
            staged_name = f"{staged_name}.{ext}"
        staged_name = clean_name(staged_name, f"input_{i}.{ext}")

        if staged_name in used:
            path = Path(staged_name)
            staged_name = f"{path.stem}_{i}{path.suffix}"
        used.add(staged_name)

        target = input_dir / staged_name
        if target.exists():
            target.unlink()
        os.symlink(mapping["path"], target)
        overrides.append(f"{node}.{parameter}={str(target)!r}")

    return overrides


def output_overrides(config_path, output_dir):
    output_dir.mkdir(exist_ok=True)
    overrides = []
    used_designations = set()
    config = configparser.ConfigParser(interpolation=None)
    config.optionxform = str
    config.read(config_path)

    for i, section in enumerate(config.sections(), 1):
        class_name = config[section].get("_class", "")
        if not class_name.startswith("Save") or "filename" not in config[section]:
            continue
        try:
            configured_filename = ast.literal_eval(config[section]["filename"])
        except (SyntaxError, ValueError):
            configured_filename = config[section]["filename"]
        filename = unique_output_name(configured_filename, used_designations, f"output_{i}.csv")
        overrides.append(f"{section}.filename={str(output_dir / filename)!r}")

    return overrides


def user_overrides(mappings):
    return [
        f"{m['node'].strip()}.{m['parameter'].strip()}={m['value'].strip()}"
        for m in mappings
        if m.get("node", "").strip() and m.get("parameter", "").strip() and m.get("value", "").strip()
    ]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--input-mappings", required=True)
    parser.add_argument("--extra-overrides", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()

    overrides = []
    overrides += input_overrides(load_json(args.input_mappings), Path("countess_inputs"))
    overrides += output_overrides(args.config, Path(args.output_dir))
    overrides += user_overrides(load_json(args.extra_overrides))

    command = ["countess_cmd", "--log", args.log_level]
    for override in overrides:
        command += ["--set", override]
    command.append(args.config)

    subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
