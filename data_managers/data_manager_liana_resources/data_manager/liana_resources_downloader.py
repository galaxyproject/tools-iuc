#!/usr/bin/env python
"""Download LIANA 1.8.1 resources and register them in Galaxy."""

import argparse
import json
import logging
import sys
from pathlib import Path

import liana as li

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger(__name__)

TARGET_ORGANISMS = (
    "anole_lizard", "c.elegans", "cat", "cattle", "chicken", "chimpanzee",
    "dog", "fruitfly", "horse", "macaque", "mouse", "opossum", "pig",
    "platypus", "rat", "s.cerevisiae", "s.pombe", "xenopus", "zebrafish",
)

METALINKS_PRESETS = {
    "metalinks_all": {
        "name": "Metalinks All",
        "description": "All Metalinks metabolite-protein edges",
        "biospecimen_location": None,
    },
    "metalinks_blood": {
        "name": "Metalinks Blood",
        "description": "Metalinks edges annotated to blood biospecimens",
        "biospecimen_location": ["Blood"],
    },
    "metalinks_plasma": {
        "name": "Metalinks Plasma",
        "description": "Metalinks edges annotated to plasma biospecimens",
        "biospecimen_location": ["Plasma"],
    },
    "metalinks_urine": {
        "name": "Metalinks Urine",
        "description": "Metalinks edges annotated to urine biospecimens",
        "biospecimen_location": ["Urine"],
    },
}

LR_NAMES = {
    "consensus": "Consensus",
    "cellchatdb": "CellChatDB",
    "cellphonedb": "CellPhoneDB",
    "connectomedb2020": "ConnectomeDB 2020",
    "mouseconsensus": "Mouse Consensus",
}


def hcop_id(organism):
    return f"hcop_human_{organism.replace('.', '_')}"


def hcop_metadata(organism):
    label = organism.replace("_", " ").replace("c.elegans", "C. elegans").replace("s.cerevisiae", "S. cerevisiae").replace("s.pombe", "S. pombe")
    return {
        "name": f"HCOP Human-{label.title()}",
        "description": f"HCOP ortholog mapping from human to {label}",
        "target_organism": organism,
    }


class LianaResourcesDownloader:
    """Download selected resources into one directory per data-table entry."""

    def __init__(self, output_json):
        self.output_json = Path(output_json)
        self.entries = []

    @staticmethod
    def available_lr_resources():
        resources = li.resource.show_resources()
        if hasattr(resources, "tolist"):
            resources = resources.tolist()
        return [str(resource) for resource in resources]

    def metadata(self, resource_id):
        if resource_id.startswith("hcop_human_"):
            for organism in TARGET_ORGANISMS:
                if resource_id == hcop_id(organism):
                    return hcop_metadata(organism)
        if resource_id in METALINKS_PRESETS:
            return METALINKS_PRESETS[resource_id]
        return {
            "name": LR_NAMES.get(resource_id, resource_id.replace("_", " ").title()),
            "description": f"LIANA ligand-receptor resource: {resource_id}",
        }

    @staticmethod
    def output_path(resource_id):
        directory = Path(resource_id)
        directory.mkdir(parents=True, exist_ok=True)
        return directory / f"{resource_id}.tsv"

    def download_lr(self, resource_id):
        frame = li.resource.select_resource(resource_name=resource_id)
        path = self.output_path(resource_id)
        frame.to_csv(path, sep="\t", index=False)
        return path

    def download_hcop(self, resource_id):
        metadata = self.metadata(resource_id)
        frame = li.resource.get_hcop_orthologs(
            target_organism=metadata["target_organism"],
            min_evidence=3,
        )
        path = self.output_path(resource_id)
        frame.to_csv(path, sep="\t", index=False)
        return path

    def download_metalinks(self, resource_id):
        metadata = METALINKS_PRESETS[resource_id]
        # LIANA 1.8.1 no longer accepts the deprecated hmdb/uniprot booleans.
        frame = li.resource.get_metalinks(
            biospecimen_location=metadata["biospecimen_location"],
        )
        path = self.output_path(resource_id)
        frame.to_csv(path, sep="\t", index=False)
        return path

    def resolve(self, selection):
        lr = self.available_lr_resources
        common_hcop = [hcop_id(org) for org in ("mouse", "rat", "chicken", "zebrafish")]
        all_hcop = [hcop_id(org) for org in TARGET_ORGANISMS]
        groups = {
            "ligand_receptor_consensus": ["consensus"],
            "ligand_receptor_main": ["consensus", "cellphonedb", "cellchatdb"],
            "ligand_receptor_single_cell": ["consensus", "cellphonedb", "cellchatdb", "connectomedb2020"],
            "ligand_receptor_all": lr,
            "orthologs_common": common_hcop,
            "orthologs_all_supported": all_hcop,
            "metalinks_bundle_all": list(METALINKS_PRESETS),
            "all": lambda: lr() + common_hcop + list(METALINKS_PRESETS),
        }
        if selection in groups:
            value = groups[selection]
            return value() if callable(value) else list(value)
        if selection in METALINKS_PRESETS or selection.startswith("hcop_human_"):
            return [selection]
        available = lr()
        if selection in available:
            return [selection]
        raise ValueError(f"Unknown resource selection: {selection}")

    def download(self, resource_id):
        if resource_id in METALINKS_PRESETS:
            return self.download_metalinks(resource_id)
        if resource_id.startswith("hcop_human_"):
            return self.download_hcop(resource_id)
        return self.download_lr(resource_id)

    def add_entry(self, resource_id, output_path):
        metadata = self.metadata(resource_id)
        self.entries.append({
            "value": resource_id,
            "name": metadata["name"],
            "description": metadata["description"],
            # The data-manager config moves this whole directory.
            "path": str(Path(output_path).parent),
        })

    def run(self, selection):
        resources = self.resolve(selection)
        log.info("Downloading %d LIANA resource(s)", len(resources))
        for resource_id in resources:
            log.info("Downloading %s", resource_id)
            output_path = self.download(resource_id)
            self.add_entry(resource_id, output_path)
        payload = {"data_tables": {"liana_resources": self.entries}}
        self.output_json.write_text(json.dumps(payload, indent=2) + "\n")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--resource", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    try:
        LianaResourcesDownloader(args.output).run(args.resource)
    except Exception as exc:
        log.exception("LIANA resource download failed: %s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
