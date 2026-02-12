#!/usr/bin/env python
"""
AlphaGenome ISM Scanner for Galaxy

In-silico saturation mutagenesis — systematically mutates every position in a
region to all 3 alt bases and scores each. Uses score_ism_variants() with
server-side chunking and parallelism.
"""

import argparse
import csv
import logging
import os
import sys

import numpy as np
from alphagenome.data import genome
from alphagenome.models import dna_client
from alphagenome.models.variant_scorers import RECOMMENDED_VARIANT_SCORERS

__version__ = "0.1.0"

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


def parse_bed(bed_path, max_regions, max_region_width):
    """Parse BED file and return list of (chrom, start, end, name) tuples."""
    regions = []
    with open(bed_path) as f:
        for line_num, line in enumerate(f):
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                logging.warning("Skipping malformed BED line %d: %s", line_num + 1, line)
                continue
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            name = fields[3] if len(fields) > 3 else f"{chrom}:{start}-{end}"
            width = end - start
            if width > max_region_width:
                logging.warning(
                    "Region %s is %dbp, exceeding max width %dbp — trimming to center %dbp",
                    name, width, max_region_width, max_region_width,
                )
                center = (start + end) // 2
                start = center - max_region_width // 2
                end = start + max_region_width
            if len(regions) >= max_regions:
                logging.warning("Reached max regions (%d), skipping remaining", max_regions)
                break
            regions.append((chrom, start, end, name))
    return regions


def run(args):
    logging.info("AlphaGenome ISM Scanner v%s", __version__)
    logging.info("Input: %s", args.input)
    logging.info("Scorers: %s", ", ".join(args.scorers))
    logging.info("Organism: %s", args.organism)
    logging.info("Sequence length: %s", args.sequence_length)
    logging.info("Max regions: %d, max region width: %dbp", args.max_regions, args.max_region_width)

    api_key = args.api_key or os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key and not args.local_model:
        logging.error("No API key provided. Set ALPHAGENOME_API_KEY or use --api-key")
        sys.exit(1)

    organism = ORGANISM_MAP[args.organism]
    seq_length = SEQUENCE_LENGTH_MAP[args.sequence_length]

    available_keys = list(RECOMMENDED_VARIANT_SCORERS.keys())
    for key in args.scorers:
        if key not in RECOMMENDED_VARIANT_SCORERS:
            logging.error("Unknown scorer key: %s (available: %s)", key, ", ".join(available_keys))
            sys.exit(1)
    selected_scorers = [RECOMMENDED_VARIANT_SCORERS[k] for k in args.scorers]

    regions = parse_bed(args.input, args.max_regions, args.max_region_width)
    if not regions:
        logging.error("No valid regions found in input BED file")
        sys.exit(1)
    logging.info("Loaded %d regions", len(regions))

    logging.info("Connecting to AlphaGenome...")
    model = create_model(api_key, local_model=args.local_model)
    logging.info("Model ready.")

    stats = {"regions": 0, "scored": 0, "errors": 0}
    row_count = 0

    with open(args.output, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow([
            "region", "position", "ref_base", "alt_base",
            "gene_id", "gene_name", "gene_type",
            "scorer", "track_name", "ontology_curie",
            "raw_score", "quantile_score",
        ])

        for region_num, (chrom, start, end, name) in enumerate(regions):
            stats["regions"] += 1
            width = end - start
            logging.info("Region %d/%d: %s (%s:%d-%d, %dbp, %d mutations)",
                         region_num + 1, len(regions), name, chrom, start, end, width, width * 3)

            try:
                interval = genome.Interval(chrom, start, end).resize(seq_length)
                ism_interval = genome.Interval(chrom, start, end, strand="+")

                results = model.score_ism_variants(
                    interval, ism_interval,
                    variant_scorers=selected_scorers,
                    organism=organism,
                    max_workers=args.max_workers,
                )

                # results is list[list[AnnData]] — outer=variants (3*width), inner=scorers
                # Each AnnData has: uns['variant'] with position/ref/alt,
                # X for raw scores, layers['quantiles'], obs for genes, var for tracks
                for var_results in results:
                    for scorer_idx, ad in enumerate(var_results):
                        variant_obj = ad.uns["variant"]
                        pos = variant_obj.position
                        ref_base = variant_obj.reference_bases
                        alt_base = variant_obj.alternate_bases
                        scorer_name = args.scorers[scorer_idx] if scorer_idx < len(args.scorers) else f"scorer_{scorer_idx}"

                        raw_scores = ad.X  # shape (n_genes, n_tracks)
                        quantile_scores = ad.layers.get("quantiles", None)

                        for gene_idx in range(ad.n_obs):
                            gene_row = ad.obs.iloc[gene_idx]
                            gene_id = str(gene_row.get("gene_id", ""))
                            gene_name = str(gene_row.get("gene_name", ""))
                            gene_type = str(gene_row.get("gene_type", ""))

                            for track_idx in range(ad.n_vars):
                                track_row = ad.var.iloc[track_idx]
                                track_name = str(track_row.get("name", ""))
                                ontology_curie = str(track_row.get("ontology_curie", ""))

                                raw = float(raw_scores[gene_idx, track_idx])
                                if np.isnan(raw):
                                    continue
                                quant = ""
                                if quantile_scores is not None:
                                    q = float(quantile_scores[gene_idx, track_idx])
                                    if not np.isnan(q):
                                        quant = f"{q:.6f}"

                                writer.writerow([
                                    name, pos, ref_base, alt_base,
                                    gene_id, gene_name, gene_type,
                                    scorer_name, track_name, ontology_curie,
                                    f"{raw:.6f}", quant,
                                ])
                                row_count += 1

                stats["scored"] += 1
                logging.info("Region %s: %d ISM variants scored", name, len(results))

            except Exception as e:
                logging.error("Error scanning region %s (%s:%d-%d): %s", name, chrom, start, end, e)
                stats["errors"] += 1

    logging.info("Wrote %d rows to %s", row_count, args.output)

    logging.info("=" * 50)
    logging.info("DONE — %d regions, %d scored, %d errors",
                 stats["regions"], stats["scored"], stats["errors"])

    if stats["errors"] > 0 and stats["scored"] == 0:
        logging.error("All regions failed. Check API key and network connectivity.")
        sys.exit(1)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="In-silico saturation mutagenesis using AlphaGenome score_ism_variants()",
    )
    parser.add_argument("--input", required=True, help="Input BED file")
    parser.add_argument("--output", required=True, help="Output TSV file")
    parser.add_argument("--api-key", default=None, help="AlphaGenome API key (or set ALPHAGENOME_API_KEY)")
    parser.add_argument(
        "--organism", choices=["human", "mouse"], default="human",
    )
    parser.add_argument(
        "--scorers", nargs="+", default=["RNA_SEQ", "ATAC"],
        help="Scorer keys from RECOMMENDED_VARIANT_SCORERS",
    )
    parser.add_argument(
        "--sequence-length", choices=list(SEQUENCE_LENGTH_MAP.keys()), default="1MB",
    )
    parser.add_argument("--max-regions", type=int, default=10)
    parser.add_argument("--max-region-width", type=int, default=200)
    parser.add_argument("--max-workers", type=int, default=5)
    parser.add_argument("--local-model", action="store_true")
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
