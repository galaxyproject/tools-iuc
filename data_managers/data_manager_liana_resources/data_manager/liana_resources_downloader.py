#!/usr/bin/env python
"""
Download LIANA resources (L-R, orthologs, metalinks) and register with Galaxy.

This script downloads selected resources from LIANA+ including:
- Ligand-receptor interaction databases
- Ortholog mappings (HCOP)
- Metabolite-protein interactions (Metalinks)

And generates Galaxy data manager JSON output for registration.
"""

import argparse
import json
import logging
import sys
from pathlib import Path

try:
    import liana as li
except ImportError as e:
    sys.exit(f"Error: Required packages not available: {e}\n"
             f"Please ensure 'liana' is installed.")

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
log = logging.getLogger(__name__)


class LianaResourcesDownloader:
    """Download all LIANA resource types and generate Galaxy data manager output."""

    # Ligand-Receptor resource metadata (17 databases)
    LR_RESOURCES_METADATA = {
        'consensus': {
            'name': 'Consensus',
            'description': 'Consensus L-R database (combined multiple sources - RECOMMENDED)',
            'type': 'ligand_receptor'
        },
        'cellchatdb': {
            'name': 'CellChat',
            'description': 'CellChat ligand-receptor database',
            'type': 'ligand_receptor'
        },
        'cellphonedb': {
            'name': 'CellPhoneDB',
            'description': 'CellPhoneDB ligand-receptor database',
            'type': 'ligand_receptor'
        },
        'connectomedb2020': {
            'name': 'Connectome',
            'description': 'Connectome ligand-receptor database',
            'type': 'ligand_receptor'
        },
        'baccin2019': {
            'name': 'Baccin et al. 2019',
            'description': 'Baccin et al. stromal cell interactions',
            'type': 'ligand_receptor'
        },
        'cellcall': {
            'name': 'CellCall',
            'description': 'CellCall ligand-receptor database',
            'type': 'ligand_receptor'
        },
        'cellinker': {
            'name': 'CellInker',
            'description': 'CellInker ligand-receptor database',
            'type': 'ligand_receptor'
        },
        'celltalkdb': {
            'name': 'CellTalkDB',
            'description': 'CellTalkDB ligand-receptor database',
            'type': 'ligand_receptor'
        },
        'embrace': {
            'name': 'EMBRACE',
            'description': 'EMBRACE extracellular matrix interactions',
            'type': 'ligand_receptor'
        },
        'guide2pharma': {
            'name': 'Guide to Pharmacology',
            'description': 'Guide to Pharmacology ligand-receptor pairs',
            'type': 'ligand_receptor'
        },
        'hpmr': {
            'name': 'HPMR',
            'description': 'Human protein-protein interactions',
            'type': 'ligand_receptor'
        },
        'icellnet': {
            'name': 'iCellNet',
            'description': 'iCellNet ligand-receptor database',
            'type': 'ligand_receptor'
        },
        'italk': {
            'name': 'iTALK',
            'description': 'iTALK ligand-receptor database',
            'type': 'ligand_receptor'
        },
        'kirouac2010': {
            'name': 'Kirouac et al. 2010',
            'description': 'Kirouac et al. signaling database',
            'type': 'ligand_receptor'
        },
        'lrdb': {
            'name': 'LRDB',
            'description': 'Ligand-Receptor Database',
            'type': 'ligand_receptor'
        },
        'mouseconsensus': {
            'name': 'Mouse Consensus',
            'description': 'Consensus L-R database for mouse (mus musculus)',
            'type': 'ligand_receptor'
        },
        'ramilowski2015': {
            'name': 'Ramilowski et al. 2015',
            'description': 'Ramilowski et al. ligand-receptor pairs',
            'type': 'ligand_receptor'
        },
    }

    # Ortholog mapping resources (HCOP)
    ORTHOLOG_RESOURCES_METADATA = {
        'hcop_human_mouse': {
            'name': 'HCOP Human-Mouse',
            'description': 'HCOP ortholog mappings: Human to Mouse (Mus musculus)',
            'type': 'ortholog_mapping',
            'url': 'https://ftp.ebi.ac.uk/pub/databases/genenames/hcop/human_mouse_hcop_fifteen_column.txt.gz'
        },
        'hcop_human_rat': {
            'name': 'HCOP Human-Rat',
            'description': 'HCOP ortholog mappings: Human to Rat (Rattus norvegicus)',
            'type': 'ortholog_mapping',
            'url': 'https://ftp.ebi.ac.uk/pub/databases/genenames/hcop/human_rat_hcop_fifteen_column.txt.gz'
        },
        'hcop_human_chicken': {
            'name': 'HCOP Human-Chicken',
            'description': 'HCOP ortholog mappings: Human to Chicken (Gallus gallus)',
            'type': 'ortholog_mapping',
            'url': 'https://ftp.ebi.ac.uk/pub/databases/genenames/hcop/human_chicken_hcop_fifteen_column.txt.gz'
        },
        'hcop_human_zebrafish': {
            'name': 'HCOP Human-Zebrafish',
            'description': 'HCOP ortholog mappings: Human to Zebrafish (Danio rerio)',
            'type': 'ortholog_mapping',
            'url': 'https://ftp.ebi.ac.uk/pub/databases/genenames/hcop/human_zebrafish_hcop_fifteen_column.txt.gz'
        },
    }

    # Metalinks resource
    METALINKS_RESOURCES_METADATA = {
        'metalinks_all': {
            'name': 'Metalinks All',
            'description': 'Metabolite-protein interactions from Metalinks database (all biospecimens)',
            'type': 'metalinks',
            'biospecimen_location': None,
            'source': None,
        },
        'metalinks_blood': {
            'name': 'Metalinks Blood',
            'description': 'Metabolite-protein interactions from blood biospecimens',
            'type': 'metalinks',
            'biospecimen_location': ['Blood'],
            'source': None,
        },
        'metalinks_plasma': {
            'name': 'Metalinks Plasma',
            'description': 'Metabolite-protein interactions from plasma biospecimens',
            'type': 'metalinks',
            'biospecimen_location': ['Plasma'],
            'source': None,
        },
        'metalinks_urine': {
            'name': 'Metalinks Urine',
            'description': 'Metabolite-protein interactions from urine biospecimens',
            'type': 'metalinks',
            'biospecimen_location': ['Urine'],
            'source': None,
        },
    }

    # Combined metadata for easy lookup
    ALL_RESOURCES_METADATA = {}

    # Resource groupings for selective download
    RESOURCE_GROUPS = {
        'ligand_receptor_all': None,  # all L-R databases discovered by show_resources()
        'ligand_receptor_consensus': ['consensus'],
        'ligand_receptor_main': ['consensus', 'cellphonedb', 'cellchatdb'],
        'ligand_receptor_single_cell': [
            'consensus', 'cellphonedb', 'cellchatdb', 'connectomedb2020'
        ],
        'orthologs_all': list(ORTHOLOG_RESOURCES_METADATA.keys()),
        'orthologs_common': ['hcop_human_mouse', 'hcop_human_rat'],
        'metalinks_bundle_all': list(METALINKS_RESOURCES_METADATA.keys()),
        'all': None,  # all resources of all types
    }

    def __init__(self, output_json):
        """Initialize downloader with output path."""
        self.output_json = output_json
        self.entries = []

        # Build combined metadata dictionary
        self.ALL_RESOURCES_METADATA = {}
        self.ALL_RESOURCES_METADATA.update(self.LR_RESOURCES_METADATA)
        self.ALL_RESOURCES_METADATA.update(self.ORTHOLOG_RESOURCES_METADATA)
        self.ALL_RESOURCES_METADATA.update(self.METALINKS_RESOURCES_METADATA)

    def get_available_lr_resources(self):
        """Get list of all available LIANA ligand-receptor resources."""
        try:
            resources = li.resource.show_resources()
            log.info(f"Found {len(resources)} available LIANA L-R resources")
            return resources
        except Exception as e:
            log.error(f"Failed to get available L-R resources: {e}")
            raise

    def download_lr_resource(self, resource_name):
        """Download a single ligand-receptor resource."""
        try:
            log.info(f"Downloading L-R resource: {resource_name}...")
            df = li.resource.select_resource(resource_name=resource_name)

            # Save as TSV
            output_path = Path(resource_name) / f"{resource_name}.tsv"
            output_path.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(output_path, sep='\t', index=False)

            log.info(f"  Saved {len(df)} interactions to {output_path}")
            return df, str(output_path)

        except Exception as e:
            log.error(f"Error downloading L-R resource {resource_name}: {e}")
            raise

    def download_ortholog_resource(self, resource_id, url):
        """Download ortholog mapping resource from HCOP."""
        try:
            log.info(f"Downloading ortholog mapping: {resource_id}...")
            df = li.resource.get_hcop_orthologs(url=url, min_evidence=3)

            # Save as TSV
            output_path = Path(resource_id) / f"{resource_id}.tsv"
            output_path.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(output_path, sep='\t', index=False)

            log.info(f"  Saved {len(df)} ortholog mappings to {output_path}")
            return df, str(output_path)

        except Exception as e:
            log.error(f"Error downloading ortholog resource {resource_id}: {e}")
            raise

    def download_metalinks_resource(self, resource_id, biospecimen_location=None, source=None):
        """Download metabolite-protein interaction resource."""
        try:
            log.info(f"Downloading metalinks: {resource_id}...")
            df = li.resource.get_metalinks(
                biospecimen_location=biospecimen_location,
                source=source
            )

            output_path = Path(resource_id) / f"{resource_id}.tsv"
            output_path.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(output_path, sep='\t', index=False)

            log.info(f"  Saved {len(df)} interactions to {output_path}")
            return df, str(output_path)

        except Exception as e:
            log.error(f"Error downloading metalinks {resource_id}: {e}")
            raise

    def create_data_table_entry(self, resource_id, output_path):
        """Create a Galaxy data table entry for a resource."""
        metadata = self.ALL_RESOURCES_METADATA.get(
            resource_id,
            {
                'name': resource_id.replace('_', ' ').title(),
                'description': f'LIANA resource: {resource_id}',
                'type': 'unknown'
            }
        )

        output_path = Path(output_path)

        return {
            'value': resource_id,
            'name': metadata['name'],
            'description': metadata['description'],
            'path': str(output_path.parent),
            'type': metadata.get('type', 'unknown')
        }

    def get_resources_to_download(self, resource_select):
        """Determine which resources to download based on selection."""
        if resource_select in self.RESOURCE_GROUPS:
            group = self.RESOURCE_GROUPS[resource_select]
            if group is None:
                if resource_select == 'all':
                    all_lr = self.get_available_lr_resources()
                    return all_lr + list(self.ORTHOLOG_RESOURCES_METADATA.keys()) + list(self.METALINKS_RESOURCES_METADATA.keys())
                if resource_select == 'ligand_receptor_all':
                    return self.get_available_lr_resources()
                return []
            return group

        if resource_select in self.ALL_RESOURCES_METADATA:
            return [resource_select]

        raise ValueError(
            f"Resource '{resource_select}' not found. "
            f"Available groups: {', '.join(self.RESOURCE_GROUPS.keys())}"
        )

    def download_resource(self, resource_id):
        """Download a single resource (detects type and calls appropriate method)."""
        metadata = self.ALL_RESOURCES_METADATA.get(resource_id)

        if not metadata:
            raise ValueError(f"Unknown resource: {resource_id}")

        resource_type = metadata.get('type')

        if resource_type == 'ligand_receptor':
            return self.download_lr_resource(resource_id)
        elif resource_type == 'ortholog_mapping':
            return self.download_ortholog_resource(resource_id, metadata['url'])
        elif resource_type == 'metalinks':
            return self.download_metalinks_resource(
                resource_id,
                metadata.get('biospecimen_location'),
                metadata.get('source')
            )
        else:
            raise ValueError(f"Unknown resource type: {resource_type}")

    def download_resources(self, resource_select):
        """Download selected resources and generate data table entries."""
        resources_to_download = self.get_resources_to_download(resource_select)

        log.info(f"Downloading {len(resources_to_download)} resource(s)...")

        for resource_id in resources_to_download:
            try:
                df, output_path = self.download_resource(resource_id)
                entry = self.create_data_table_entry(resource_id, output_path)
                self.entries.append(entry)

            except Exception as e:
                log.error(f"Failed to download {resource_id}: {e}")
                raise

        log.info(f"Successfully downloaded {len(self.entries)} resource(s)")
        return self.entries

    def write_output(self):
        """Write Galaxy data manager JSON output."""
        output_data = {
            'data_tables': {
                'liana_resources': self.entries
            }
        }

        with open(self.output_json, 'w') as f:
            json.dump(output_data, f, indent=2)

        log.info(f"Wrote data manager output to {self.output_json}")

    def run(self, resource_select):
        """Execute the download and output generation."""
        try:
            self.download_resources(resource_select)
            self.write_output()
            log.info("Data manager execution completed successfully")
            return 0
        except Exception as e:
            log.error(f"Data manager execution failed: {e}")
            return 1


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Download LIANA resources (L-R databases, orthologs, metalinks) for Galaxy',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
        Available --resource values:

        LIGAND-RECEPTOR RESOURCE GROUPS:
        ligand_receptor_consensus     - Consensus only
        ligand_receptor_main          - Consensus + CellPhoneDB + CellChat
        ligand_receptor_single_cell   - Consensus + CellPhoneDB + CellChat + ConnectomeDB2020
        ligand_receptor_all           - All available LIANA ligand-receptor resources

        ORTHOLOG RESOURCE GROUPS:
        orthologs_common              - Human-Mouse + Human-Rat
        orthologs_all                 - Human-Mouse + Human-Rat + Human-Chicken + Human-Zebrafish

        METALINKS RESOURCE GROUPS:
        metalinks_blood               - Blood biospecimens only
        metalinks_plasma              - Plasma biospecimens only
        metalinks_all                 - All biospecimens only
        metalinks_bundle_all          - All predefined Metalinks presets

        ALL RESOURCES:
        all                           - All ligand-receptor, ortholog, and Metalinks resources

        Examples:
        python liana_resources_downloader.py --resource ligand_receptor_consensus --output output.json
        python liana_resources_downloader.py --resource orthologs_common --output orthologs.json
        python liana_resources_downloader.py --resource metalinks_all --output metalinks.json
        '''
    )

    parser.add_argument(
        '--resource',
        default='ligand_receptor_consensus',
        help='Resources to download. See examples above. (default: ligand_receptor_consensus)'
    )

    parser.add_argument(
        '--output',
        required=True,
        help='Output JSON file for Galaxy data manager (required)'
    )

    return parser.parse_args()


def main():
    """Main entry point."""
    args = parse_arguments()

    try:
        downloader = LianaResourcesDownloader(args.output)
        return downloader.run(args.resource)
    except Exception as e:
        log.error(f"Fatal error: {e}")
        return 1


if __name__ == '__main__':
    sys.exit(main())
