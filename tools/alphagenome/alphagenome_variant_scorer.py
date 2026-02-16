#!/usr/bin/env python
"""
AlphaGenome Variant Scorer for Galaxy

Uses score_variant() for server-side gene-level aggregation with spatial masking
and empirical quantile normalization. Outputs structured per-gene, per-track scores
via tidy_scores().
"""

import argparse
import logging
import os
import sys

import cyvcf2
from alphagenome.data import genome
from alphagenome.models import dna_client
from alphagenome.models.variant_scorers import RECOMMENDED_VARIANT_SCORERS, tidy_scores

__version__ = "0.5.1"

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
    if local_model:
        from alphagenome_research.model import dna_model
        return dna_model.create_from_huggingface("all_folds")
    return dna_client.create(api_key)


def run(args):
    logging.info("AlphaGenome Variant Scorer v%s", __version__)
    logging.info("Input: %s", args.input)
    logging.info("Scorers: %s", ", ".join(args.scorers))
    logging.info("Organism: %s", args.organism)
    logging.info("Sequence length: %s", args.sequence_length)
    logging.info("Max variants: %d", args.max_variants)

    if args.test_fixture:
        import json
        import pandas as pd
        with open(args.test_fixture) as f:
            fixture_data = json.load(f)
        df = pd.DataFrame(fixture_data["rows"], columns=fixture_data["columns"])
        df.to_csv(args.output, sep="\t", index=False)
        logging.info("Fixture mode: wrote %d rows to %s", len(df), args.output)
        return

    api_key = args.api_key or os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key and not args.local_model:
        logging.error("No API key provided. Set ALPHAGENOME_API_KEY or use --api-key")
        sys.exit(1)

    organism = ORGANISM_MAP[args.organism]
    seq_length = SEQUENCE_LENGTH_MAP[args.sequence_length]

    available_keys = list(RECOMMENDED_VARIANT_SCORERS.keys())
    selected_keys = args.scorers
    for key in selected_keys:
        if key not in RECOMMENDED_VARIANT_SCORERS:
            logging.error("Unknown scorer key: %s (available: %s)", key, ", ".join(available_keys))
            sys.exit(1)
    selected_scorers = [RECOMMENDED_VARIANT_SCORERS[k] for k in selected_keys]
    logging.info("Using %d scorers", len(selected_scorers))

    logging.info("Connecting to AlphaGenome...")
    model = create_model(api_key, local_model=args.local_model)
    logging.info("Model ready.")

    vcf_reader = cyvcf2.VCF(args.input)

    stats = {"total": 0, "scored": 0, "errors": 0, "skipped": 0}
    all_rows = []

    try:
        for variant_num, record in enumerate(vcf_reader):
            stats["total"] += 1

            if variant_num >= args.max_variants:
                stats["skipped"] += 1
                continue

            if variant_num > 0 and variant_num % 10 == 0:
                logging.info(
                    "Progress: %d/%d variants processed (%d scored, %d errors)",
                    variant_num, args.max_variants, stats["scored"], stats["errors"],
                )

            chrom = record.CHROM
            pos = record.POS
            ref = record.REF

            if not record.ALT:
                continue

            alt = record.ALT[0]
            variant_id = f"{chrom}:{pos}:{ref}>{alt}"

            try:
                variant = genome.Variant(
                    chromosome=chrom,
                    position=pos,
                    reference_bases=ref,
                    alternate_bases=alt,
                )
                interval = variant.reference_interval.resize(seq_length)

                scores = model.score_variant(
                    interval, variant, selected_scorers, organism=organism,
                )

                df = tidy_scores(scores)
                all_rows.append(df)
                stats["scored"] += 1

                logging.debug("Scored %s: %d rows", variant_id, len(df))

            except Exception as e:
                logging.error("Error scoring %s: %s", variant_id, e)
                stats["errors"] += 1

    finally:
        vcf_reader.close()

    if all_rows:
        import pandas as pd
        combined = pd.concat(all_rows, ignore_index=True)
        combined.to_csv(args.output, sep="\t", index=False)
        logging.info("Wrote %d rows to %s", len(combined), args.output)
    else:
        with open(args.output, "w") as f:
            f.write("variant_id\n")
        logging.warning("No variants scored successfully")

    logging.info("=" * 50)
    logging.info("DONE — %d total, %d scored, %d errors, %d skipped (over limit)",
                 stats["total"], stats["scored"], stats["errors"], stats["skipped"])

    if stats["errors"] > 0 and stats["scored"] == 0:
        logging.error("All variants failed. Check API key and network connectivity.")
        sys.exit(1)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Score variants using AlphaGenome score_variant() with gene-level aggregation",
    )
    parser.add_argument("--input", required=True, help="Input VCF file")
    parser.add_argument("--output", required=True, help="Output TSV file")
    parser.add_argument("--api-key", default=None, help="AlphaGenome API key (or set ALPHAGENOME_API_KEY)")
    parser.add_argument(
        "--organism", choices=["human", "mouse"], default="human",
    )
    parser.add_argument(
        "--scorers", nargs="+", default=["RNA_SEQ", "ATAC", "SPLICE_SITES"],
        help="Scorer keys from RECOMMENDED_VARIANT_SCORERS",
    )
    parser.add_argument(
        "--sequence-length", choices=list(SEQUENCE_LENGTH_MAP.keys()), default="1MB",
    )
    parser.add_argument(
        "--max-variants", type=int, default=100,
    )
    parser.add_argument("--local-model", action="store_true")
    parser.add_argument("--test-fixture", default=None,
                        help="Test fixture JSON for CI testing (bypasses API)")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
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
