#!/usr/bin/env python3
import argparse
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


def input_overrides(mappings, input_dir):
    input_dir.mkdir(exist_ok=True)
    overrides = []
    used = set()

    for i, mapping in enumerate(mappings, 1):
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
        overrides.append(f"{mapping['node'].strip()}.{mapping['parameter'].strip()}={str(target)!r}")

    return overrides


def output_overrides(mappings, output_dir):
    output_dir.mkdir(exist_ok=True)
    overrides = []

    for i, mapping in enumerate(mappings, 1):
        filename = clean_name(mapping.get("filename"), f"output_{i}.csv")
        path = output_dir / filename
        overrides.append(f"{mapping['node'].strip()}.{mapping['parameter'].strip()}={str(path)!r}")

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
    parser.add_argument("--output-mappings", required=True)
    parser.add_argument("--extra-overrides", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()

    overrides = []
    overrides += input_overrides(load_json(args.input_mappings), Path("countess_inputs"))
    overrides += output_overrides(load_json(args.output_mappings), Path(args.output_dir))
    overrides += user_overrides(load_json(args.extra_overrides))

    command = ["countess_cmd", "--log", args.log_level]
    for override in overrides:
        command += ["--set", override]
    command.append(args.config)

    subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
