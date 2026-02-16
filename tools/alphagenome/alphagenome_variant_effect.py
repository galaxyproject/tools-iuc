#!/usr/bin/env python
"""
AlphaGenome Variant Effect Predictor for Galaxy

POC tool that scores genetic variants using the real AlphaGenome API.
Uses predict_variant() to compute log-fold-change effect scores per output type.
"""

import argparse
import logging
import os
import sys

import cyvcf2
import numpy as np
from alphagenome.data import genome
from alphagenome.models import dna_client

__version__ = "0.5.1"

# Map CLI output type names to dna_client.OutputType enum values
OUTPUT_TYPE_MAP = {
    "RNA_SEQ": dna_client.OutputType.RNA_SEQ,
    "ATAC": dna_client.OutputType.ATAC,
    "CAGE": dna_client.OutputType.CAGE,
    "DNASE": dna_client.OutputType.DNASE,
    "CHIP_HISTONE": dna_client.OutputType.CHIP_HISTONE,
    "CHIP_TF": dna_client.OutputType.CHIP_TF,
    "SPLICE_SITES": dna_client.OutputType.SPLICE_SITES,
    "SPLICE_SITE_USAGE": dna_client.OutputType.SPLICE_SITE_USAGE,
    "SPLICE_JUNCTIONS": dna_client.OutputType.SPLICE_JUNCTIONS,
    "CONTACT_MAPS": dna_client.OutputType.CONTACT_MAPS,
    "PROCAP": dna_client.OutputType.PROCAP,
}

# Map output type names to VCF INFO field IDs
INFO_FIELD_MAP = {
    "RNA_SEQ": "AG_RNA_LFC",
    "ATAC": "AG_ATAC_LFC",
    "CAGE": "AG_CAGE_LFC",
    "DNASE": "AG_DNASE_LFC",
    "CHIP_HISTONE": "AG_HISTONE_LFC",
    "CHIP_TF": "AG_TF_LFC",
    "SPLICE_SITES": "AG_SPLICE_LFC",
    "SPLICE_SITE_USAGE": "AG_SPLICEUSE_LFC",
    "SPLICE_JUNCTIONS": "AG_SPLICEJNC_LFC",
    "CONTACT_MAPS": "AG_CONTACT_LFC",
    "PROCAP": "AG_PROCAP_LFC",
}

ORGANISM_MAP = {
    "human": dna_client.Organism.HOMO_SAPIENS,
    "mouse": dna_client.Organism.MUS_MUSCULUS,
}

SEQUENCE_LENGTH_MAP = {
    "16KB": 16_384,
    "128KB": 131_072,
    "512KB": 524_288,
    "1MB": 1_048_576,
}


def create_model(api_key, local_model=False):
    """Create an AlphaGenome model client.

    Swappable between API and local model for future TACC deployment.
    """
    if local_model:
        from alphagenome_research.model import dna_model
        return dna_model.create_from_huggingface("all_folds")
    return dna_client.create(api_key)


def compute_max_abs_lfc(ref_values, alt_values):
    """Compute max absolute log-fold-change between ref and alt predictions."""
    ref_vals = np.asarray(ref_values, dtype=np.float64)
    alt_vals = np.asarray(alt_values, dtype=np.float64)
    epsilon = 1e-6
    lfc = np.log2((alt_vals + epsilon) / (ref_vals + epsilon))
    return float(np.max(np.abs(lfc)))


def get_track_values(prediction_output, output_type_name):
    """Extract the .values array from a prediction output for a given type.

    The attribute name on the output object is the lowercase version of the
    OutputType enum name (e.g. OutputType.RNA_SEQ -> output.rna_seq).
    """
    attr_name = output_type_name.lower()
    track = getattr(prediction_output, attr_name, None)
    if track is None:
        return None
    return track.values


def run(args):
    """Main processing loop."""
    logging.info("AlphaGenome Variant Effect Predictor v%s", __version__)
    logging.info("Input: %s", args.input)
    logging.info("Output types: %s", ", ".join(args.output_types))
    logging.info("Organism: %s", args.organism)
    logging.info("Sequence length: %s", args.sequence_length)
    logging.info("Max variants: %d", args.max_variants)

    # Fixture mode for CI testing (bypasses API)
    fixture_lookup = None
    if args.test_fixture:
        import json
        with open(args.test_fixture) as f:
            fixture_data = json.load(f)
        fixture_lookup = {}
        for v in fixture_data["variants"]:
            key = (v["chrom"], v["pos"], v["ref"], v["alt"])
            fixture_lookup[key] = v["scores"]
        logging.info("Fixture mode: %d pre-computed variants", len(fixture_lookup))

    if fixture_lookup is None:
        # Resolve API key from CLI arg or environment variable
        api_key = args.api_key or os.environ.get("ALPHAGENOME_API_KEY")
        if not api_key and not args.local_model:
            logging.error("No API key provided. Set ALPHAGENOME_API_KEY or use --api-key")
            sys.exit(1)

        # Resolve parameters
        organism = ORGANISM_MAP[args.organism]
        seq_length = SEQUENCE_LENGTH_MAP[args.sequence_length]
        requested_outputs = [OUTPUT_TYPE_MAP[t] for t in args.output_types]
        ontology_terms = []
        if args.ontology_terms:
            ontology_terms = [t.strip() for t in args.ontology_terms.split(",") if t.strip()]

        # Create model
        logging.info("Connecting to AlphaGenome...")
        model = create_model(api_key, local_model=args.local_model)
        logging.info("Model ready.")

    # Open input VCF
    vcf_reader = cyvcf2.VCF(args.input)

    # Add INFO headers for each selected output type
    for otype in args.output_types:
        info_id = INFO_FIELD_MAP[otype]
        vcf_reader.add_info_to_header({
            "ID": info_id,
            "Number": "A",
            "Type": "Float",
            "Description": f"AlphaGenome {otype} max absolute log-fold-change",
        })
    vcf_reader.add_info_to_header({
        "ID": "AG_MAX_EFFECT",
        "Number": "A",
        "Type": "Float",
        "Description": "AlphaGenome max effect across all selected output types",
    })
    vcf_reader.add_to_header(f"##AlphaGenomeVariantEffectVersion={__version__}")

    # Open output VCF
    vcf_writer = cyvcf2.Writer(args.output, vcf_reader)

    stats = {"total": 0, "scored": 0, "errors": 0, "skipped": 0}

    try:
        for variant_num, record in enumerate(vcf_reader):
            stats["total"] += 1

            if variant_num >= args.max_variants:
                stats["skipped"] += 1
                vcf_writer.write_record(record)
                continue

            if variant_num > 0 and variant_num % 10 == 0:
                logging.info(
                    "Progress: %d/%d variants processed (%d scored, %d errors)",
                    variant_num, args.max_variants, stats["scored"], stats["errors"],
                )

            chrom = record.CHROM
            pos = record.POS  # 1-based, matches API expectation
            ref = record.REF

            # Process first ALT allele only for POC
            if not record.ALT or len(record.ALT) == 0:
                vcf_writer.write_record(record)
                continue

            alt = record.ALT[0]

            try:
                all_scores = []
                if fixture_lookup is not None:
                    fixture_scores = fixture_lookup.get((chrom, pos, ref, alt), {})
                    for otype in args.output_types:
                        if otype in fixture_scores:
                            score = fixture_scores[otype]
                            info_id = INFO_FIELD_MAP[otype]
                            record.INFO[info_id] = round(score, 6)
                            all_scores.append(score)
                else:
                    variant = genome.Variant(
                        chromosome=chrom,
                        position=pos,
                        reference_bases=ref,
                        alternate_bases=alt,
                    )
                    interval = variant.reference_interval.resize(seq_length)

                    outputs = model.predict_variant(
                        interval=interval,
                        variant=variant,
                        organism=organism,
                        ontology_terms=ontology_terms,
                        requested_outputs=requested_outputs,
                    )

                    for otype in args.output_types:
                        ref_vals = get_track_values(outputs.reference, otype)
                        alt_vals = get_track_values(outputs.alternate, otype)
                        if ref_vals is not None and alt_vals is not None:
                            score = compute_max_abs_lfc(ref_vals, alt_vals)
                            info_id = INFO_FIELD_MAP[otype]
                            record.INFO[info_id] = round(score, 6)
                            all_scores.append(score)
                        else:
                            logging.warning(
                                "No %s track in output for %s:%d %s>%s",
                                otype, chrom, pos, ref, alt,
                            )

                if all_scores:
                    record.INFO["AG_MAX_EFFECT"] = round(max(all_scores), 6)
                    stats["scored"] += 1
                else:
                    stats["errors"] += 1

            except Exception as e:
                logging.error("Error scoring %s:%d %s>%s: %s", chrom, pos, ref, alt, e)
                stats["errors"] += 1

            vcf_writer.write_record(record)

    finally:
        vcf_writer.close()
        vcf_reader.close()

    # Report
    logging.info("=" * 50)
    logging.info("DONE — %d total, %d scored, %d errors, %d skipped (over limit)",
                 stats["total"], stats["scored"], stats["errors"], stats["skipped"])
    logging.info("Output: %s", args.output)

    if stats["errors"] > 0 and stats["scored"] == 0:
        logging.error("All variants failed. Check API key and network connectivity.")
        sys.exit(1)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Score genetic variants using AlphaGenome predict_variant()",
    )
    parser.add_argument("--input", required=True, help="Input VCF file")
    parser.add_argument("--output", required=True, help="Output VCF file")
    parser.add_argument("--api-key", default=None, help="AlphaGenome API key (or set ALPHAGENOME_API_KEY)")
    parser.add_argument(
        "--organism", choices=["human", "mouse"], default="human",
        help="Organism (default: human)",
    )
    parser.add_argument(
        "--output-types", nargs="+", choices=list(OUTPUT_TYPE_MAP.keys()),
        default=["RNA_SEQ"],
        help="AlphaGenome output types to predict (default: RNA_SEQ)",
    )
    parser.add_argument(
        "--ontology-terms", default=None,
        help="Comma-separated ontology terms (e.g. UBERON:0001157,CL:0000746)",
    )
    parser.add_argument(
        "--sequence-length", choices=list(SEQUENCE_LENGTH_MAP.keys()), default="1MB",
        help="Prediction window size (default: 1MB)",
    )
    parser.add_argument(
        "--max-variants", type=int, default=100,
        help="Maximum variants to process (default: 100)",
    )
    parser.add_argument(
        "--local-model", action="store_true",
        help="Use local HuggingFace model instead of API",
    )
    parser.add_argument("--test-fixture", default=None,
                        help="Test fixture JSON for CI testing (bypasses API)")
    parser.add_argument("--verbose", action="store_true", help="Debug logging")
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}",
    )
    return parser.parse_args()


def main():
    args = parse_arguments()

    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stderr)],
    )

    try:
        run(args)
    except KeyboardInterrupt:
        logging.error("Interrupted")
        sys.exit(130)
    except Exception as e:
        logging.error("Fatal error: %s", e)
        if args.verbose:
            logging.exception("Details:")
        sys.exit(1)


if __name__ == "__main__":
    main()
