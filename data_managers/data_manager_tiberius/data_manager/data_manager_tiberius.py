#!/usr/bin/env python

import sys
import argparse
import json
import os
import shutil
import tarfile
import tempfile
from pathlib import Path

try:
    # For Python 3.0 and later
    from urllib.request import Request, urlopen
except ImportError:
    # Fall back to Python 2 imports
    from urllib2 import Request, urlopen

# Base URL for downloads
BASE_URL = "https://bioinf.uni-greifswald.de/bioinf/tiberius/models/"

# Define filename patterns for each species and model version
FILENAME_MAP = {
    "Chlorophyta": {
        "unmasked": "chlorophyta.tar.gz"
    },
    "Bacillariophyta": {
        "softmasked": "diatoms_weights.tar.gz",
        "unmasked": "diatoms_unmasked.tar.gz"
    },
    "Eudicotyledons": {
        "softmasked": "eudicotyledons_weights.tar.gz"
    },
    "Fungi": {
        "unmasked": "fungi.tar.gz"
    },
    "Insecta": {
        "softmasked": "insecta_weights.tar.gz",
        "unmasked": "insecta_unmasked_v2.tar.gz"
    },
    "Mammalia": {
        "softmasked": "tiberius_weights_v2.tar.gz",
        "unmasked": "tiberius_nosm_weights_v2.tar.gz"
    },
    "Monocotyledonae": {
        "softmasked": "monocotyledonae_weights.tar.gz"
    },
    "Vertebrata": {
        "unmasked": "vertebrates_weights.tar.gz"
    }
}

# Display names with dates
NAMES = {
    "Chlorophyta_unmasked": "Chlorophyta (Green Algae) - Unmasked (2026-01-26)",
    "Bacillariophyta_softmasked": "Bacillariophyta (Diatoms) - Softmasked (2025-10-24)",
    "Bacillariophyta_unmasked": "Bacillariophyta (Diatoms) - Unmasked (2026-02-06)",
    "Eudicotyledons_softmasked": "Eudicotyledons (Eudicots) - Softmasked (2025-10-24)",
    "Fungi_unmasked": "Fungi - Unmasked (2026-02-06)",
    "Insecta_softmasked": "Insecta (Insects) - Softmasked (2025-12-03)",
    "Insecta_unmasked": "Insecta (Insects) - Unmasked (2026-02-15)",
    "Mammalia_softmasked": "Mammalia (Mammals) - Softmasked (2025-05-12)",
    "Mammalia_unmasked": "Mammalia (Mammals) - Unmasked (2025-05-12)",
    "Monocotyledonae_softmasked": "Monocotyledonae (Monocots) - Softmasked (2025-10-24)",
    "Vertebrata_unmasked": "Vertebrata (Vertebrates) - Unmasked (2025-11-10)"
}

# Define which species have which model version options
# Note: "both" is a valid option for species that have both softmasked and unmasked versions
SPECIES_MODEL_OPTIONS = {
    "all": ["softmasked", "unmasked", "both"],
    "Chlorophyta": ["unmasked"],  # Only unmasked available
    "Bacillariophyta": ["softmasked", "unmasked", "both"],  # Both versions + both
    "Eudicotyledons": ["softmasked"],  # Only softmasked available
    "Fungi": ["unmasked"],  # Only unmasked available
    "Insecta": ["softmasked", "unmasked", "both"],  # Both versions + both
    "Mammalia": ["softmasked", "unmasked", "both"],  # Both versions + both
    "Monocotyledonae": ["softmasked"],  # Only softmasked available
    "Vertebrata": ["unmasked"]  # Only unmasked available
}

# Define which species to download for bulk operations
BULK_SOFTMASKED_SPECIES = ["Bacillariophyta", "Eudicotyledons", "Insecta", "Mammalia", "Monocotyledonae"]
BULK_UNMASKED_SPECIES = ["Chlorophyta", "Bacillariophyta", "Fungi", "Insecta", "Mammalia", "Vertebrata"]


def get_filename(species, masking):
    """Get the filename for a given species and masking option"""
    if species not in FILENAME_MAP:
        raise ValueError(f"Unknown species: {species}")

    species_map = FILENAME_MAP[species]
    masking_str = str(masking)

    if masking_str not in species_map:
        raise ValueError(f"Masking option '{masking}' not available for species '{species}'")

    return species_map[masking_str]


def download_file(url, destination):
    """Download a file from URL to destination path"""
    print(f"Downloading from: {url}")

    req = Request(url)
    with urlopen(req) as src, open(destination, 'wb') as dst:
        while True:
            chunk = src.read(2**10)
            if chunk:
                dst.write(chunk)
            else:
                break


def extract_tar_gz(tarfname, extract_dir):
    """Extract tar.gz file safely"""
    with tarfile.open(tarfname, "r:gz") as tar:
        # Check for path traversal
        for member in tar.getmembers():
            member_path = os.path.join(extract_dir, member.name)
            if not os.path.realpath(member_path).startswith(os.path.realpath(extract_dir)):
                raise Exception("Attempted Path Traversal in Tar File")

        tar.extractall(extract_dir)
        return tar.getnames()[0]  # Return the top-level directory name


def download_model(species, masking, workdir):
    """
    Download Tiberius model weights to workdir and extract

    Returns the name of the resulting directory
    """
    filename = get_filename(species, masking)
    url = BASE_URL + filename
    tarfname = os.path.join(workdir, filename)

    download_file(url, tarfname)
    dirname = extract_tar_gz(tarfname, workdir)
    os.remove(tarfname)

    return dirname


def get_data_manager_entries(species, masking, target_directory):
    """Generate data manager entries for downloaded models"""
    entries = []

    if species == "all":
        # Handle bulk downloads
        if masking == "softmasked":
            # Download softmasked only for species that have them
            species_list = BULK_SOFTMASKED_SPECIES
        elif masking == "unmasked":
            # Download unmasked only for species that have them
            species_list = BULK_UNMASKED_SPECIES
        else:  # masking == "both"
            # Download both versions for all species that have them
            # First download softmasked
            for sp in BULK_SOFTMASKED_SPECIES:
                try:
                    entries.extend(download_and_create_entry(sp, "softmasked", target_directory))
                except Exception as e:
                    print(f"Warning: Failed to download {sp} (softmasked): {e}")

            # Then download unmasked
            for sp in BULK_UNMASKED_SPECIES:
                try:
                    entries.extend(download_and_create_entry(sp, "unmasked", target_directory))
                except Exception as e:
                    print(f"Warning: Failed to download {sp} (unmasked): {e}")

            return entries

        # For single masking option bulk download
        for sp in species_list:
            try:
                entries.extend(download_and_create_entry(sp, masking, target_directory))
            except Exception as e:
                print(f"Warning: Failed to download {sp}: {e}")

        return entries

    else:
        # Single species download
        # Handle the "both" masking option for a single species
        if masking == "both":
            # Download both versions for this species
            # Get both softmasked and unmasked versions
            available_versions = []
            if "softmasked" in FILENAME_MAP.get(species, {}):
                available_versions.append("softmasked")
            if "unmasked" in FILENAME_MAP.get(species, {}):
                available_versions.append("unmasked")
            
            for version in available_versions:
                try:
                    entries.extend(download_and_create_entry(species, version, target_directory))
                except Exception as e:
                    print(f"Warning: Failed to download {species} ({version}): {e}")
            
            return entries
        else:
            # Single version download
            return download_and_create_entry(species, masking, target_directory)


def download_and_create_entry(species, masking, target_directory):
    """Download a model and create a data manager entry"""
    entries = []
    workdir = tempfile.mkdtemp()

    try:
        # masking is already "softmasked" or "unmasked"
        masking_type = masking
        db_id = f"{species}_{masking_type}"

        # Download and extract
        path = download_model(species, masking, workdir)

        # Move files to target directory
        target_path = os.path.join(target_directory, species, masking_type)
        os.makedirs(target_path, exist_ok=True)

        # Move extracted contents
        extracted_path = os.path.join(workdir, path)
        for item in os.listdir(extracted_path):
            src = os.path.join(extracted_path, item)
            dst = os.path.join(target_path, item)
            if os.path.exists(dst):
                if os.path.isdir(dst):
                    shutil.rmtree(dst)
                else:
                    os.remove(dst)
            shutil.move(src, dst)

        # Create entry
        entry = {
            "value": db_id,
            "name": NAMES.get(db_id, f"{species} - {masking_type}"),
            "path": os.path.join(species, masking_type),
            "type": masking_type
        }

        entries.append(entry)
        print(f"Successfully installed {species} ({masking_type}) to {target_path}")

    finally:
        # Clean up temporary directory
        shutil.rmtree(workdir, ignore_errors=True)

    return entries


def main(species, masking, outjson, target_directory=None):
    """
    Main function to download and install Tiberius models

    Args:
        species: Species name or "all"
        masking: "softmasked", "unmasked", or "both"
        outjson: Path to output JSON file
        target_directory: Target directory for extracted files (optional)
    """

    # Validate inputs
    if species not in SPECIES_MODEL_OPTIONS:
        raise ValueError(f"Unknown species: {species}")

    if masking not in SPECIES_MODEL_OPTIONS[species]:
        raise ValueError(f"Masking option '{masking}' not available for species '{species}'")

    # Determine target directory
    if not target_directory:
        # Try to read from Galaxy data manager convention
        try:
            if os.path.exists(outjson):
                with open(outjson) as fh:
                    params = json.load(fh)
                target_directory = params['output_data'][0]['extra_files_path']
            else:
                # Fallback to current directory
                target_directory = os.getcwd()
        except (json.JSONDecodeError, KeyError, IndexError) as e:
            print(f"Warning: Could not parse input JSON: {e}")
            target_directory = os.getcwd()

    # Create target directory if it doesn't exist
    os.makedirs(target_directory, exist_ok=True)

    # Get all entries
    entries = get_data_manager_entries(species, masking, target_directory)

    if not entries:
        raise ValueError("No entries were created. Download may have failed.")

    # Create the data manager output
    data_manager_json = {
        "data_tables": {
            "tiberius": entries if len(entries) > 1 else entries[0]
        }
    }

    # Write output JSON
    with open(outjson, 'w') as fh:
        json.dump(data_manager_json, fh, sort_keys=False, indent=2)

    print(f"Successfully installed {len(entries)} model(s)")
    print(f"Output JSON written to: {outjson}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Data Manager for Tiberius models')
    parser.add_argument('--out', action='store', required=True,
                       help='JSON output filename')
    parser.add_argument('--species', action='store', required=True,
                       choices=list(SPECIES_MODEL_OPTIONS.keys()),
                       help='Target species or "all" for all species')
    parser.add_argument('--masking', action='store', required=True,
                       choices=['softmasked', 'unmasked', 'both'],
                       help='Masking option: softmasked, unmasked, or both')
    parser.add_argument('--target_dir', action='store',
                       help='Target directory for extracted files (optional, for testing)')

    args = parser.parse_args()

    try:
        main(args.species, args.masking, args.out, args.target_dir)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
