#!/usr/bin/env python
"""
AlphaGenome ISM Scanner for Galaxy

In-silico saturation mutagenesis — systematically mutates every position in a
region to all 3 alt bases and scores each. Uses score_ism_variants() with
server-side chunking and parallelism.
"""

import argparse
import csv
import json
import logging
import os
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
from alphagenome.data import genome
from alphagenome.models import dna_client
from alphagenome.models.variant_scorers import RECOMMENDED_VARIANT_SCORERS


def _col_as_list(df, col):
    """Extract a column as a list of strings; empty strings if the column is missing."""
    if col in df.columns:
        return df[col].astype(str).tolist()
    return [""] * len(df)


ISM_OUTPUT_COLUMNS = [
    "region", "position", "ref_base", "alt_base",
    "gene_id", "gene_name", "gene_type",
    "scorer", "track_name", "ontology_curie",
    "raw_score", "quantile_score",
]


def _write_chunked_tsv(outfile, columns):
    df = pd.DataFrame(columns)
    df.to_csv(outfile, sep="\t", header=False, index=False,
              float_format="%.6f", na_rep="",
              lineterminator="\r\n")
    return len(df)


def _write_region_results(outfile, name, results, scorers_arg, flush_rows=200_000):
    """Post-process ISM results for one region and stream TSV rows to outfile.

    Returns the number of non-NaN rows written. obs/var metadata is identical
    across variants for a given scorer, so extract once per (region, scorer)
    then flatten non-NaN cells into a DataFrame and let pandas' C path handle
    CSV formatting.
    """
    chunks = {col: [] for col in ISM_OUTPUT_COLUMNS}
    pending_rows = 0
    per_scorer_meta = {}
    row_count = 0

    def flush_pending():
        nonlocal pending_rows, row_count
        if pending_rows == 0:
            return
        row_count += _write_chunked_tsv(
            outfile,
            {col: np.concatenate(values) for col, values in chunks.items()},
        )
        for values in chunks.values():
            values.clear()
        pending_rows = 0

    for var_results in results:
        for scorer_idx, ad in enumerate(var_results):
            variant_obj = ad.uns["variant"]
            pos = variant_obj.position
            ref_base = variant_obj.reference_bases
            alt_base = variant_obj.alternate_bases
            scorer_name = (
                scorers_arg[scorer_idx]
                if scorer_idx < len(scorers_arg)
                else f"scorer_{scorer_idx}"
            )

            if scorer_idx not in per_scorer_meta:
                per_scorer_meta[scorer_idx] = (
                    np.asarray(_col_as_list(ad.obs, "gene_id")),
                    np.asarray(_col_as_list(ad.obs, "gene_name")),
                    np.asarray(_col_as_list(ad.obs, "gene_type")),
                    np.asarray(_col_as_list(ad.var, "name")),
                    np.asarray(_col_as_list(ad.var, "ontology_curie")),
                )
            gene_ids, gene_names, gene_types, track_names, track_curies = per_scorer_meta[scorer_idx]

            raw_scores = np.asarray(ad.X, dtype=float)
            quantile_scores = ad.layers.get("quantiles", None)
            if quantile_scores is not None:
                quantile_scores = np.asarray(quantile_scores, dtype=float)

            valid_mask = ~np.isnan(raw_scores)
            gi, ti = np.where(valid_mask)
            if gi.size == 0:
                continue

            if quantile_scores is not None:
                q_flat = quantile_scores[gi, ti]
            else:
                q_flat = np.full(gi.size, np.nan)

            size = gi.size
            chunks["region"].append(np.full(size, name, dtype=object))
            chunks["position"].append(np.full(size, pos))
            chunks["ref_base"].append(np.full(size, ref_base, dtype=object))
            chunks["alt_base"].append(np.full(size, alt_base, dtype=object))
            chunks["gene_id"].append(gene_ids[gi])
            chunks["gene_name"].append(gene_names[gi])
            chunks["gene_type"].append(gene_types[gi])
            chunks["scorer"].append(np.full(size, scorer_name, dtype=object))
            chunks["track_name"].append(track_names[ti])
            chunks["ontology_curie"].append(track_curies[ti])
            chunks["raw_score"].append(raw_scores[gi, ti])
            chunks["quantile_score"].append(q_flat)
            pending_rows += size
            if pending_rows >= flush_rows:
                flush_pending()
    flush_pending()
    return row_count


def _load_mock_ism_results(path):
    """Load JSON describing mock AnnData inputs for CI testing of the post-processing path.

    The real --test-fixture path bypasses post-processing entirely (it just dumps
    pre-computed TSV rows). This loader constructs the minimal AnnData interface
    consumed by _write_region_results so CI can exercise the vectorized code
    path against controlled multi-region / NaN / missing-quantile cases.
    """
    class _MockAd:
        def __init__(self, obs, var, X, quantiles, variant):
            self.obs = pd.DataFrame(obs)
            self.var = pd.DataFrame(var)
            self.X = np.asarray(X, dtype=float)
            self.layers = (
                {"quantiles": np.asarray(quantiles, dtype=float)}
                if quantiles is not None
                else {}
            )
            self.uns = {"variant": SimpleNamespace(**variant)}

    with open(path) as f:
        data = json.load(f)

    regions = []
    for region_data in data["regions"]:
        region_results = []
        for var_entry in region_data["variants"]:
            variant_meta = {
                "position": var_entry["position"],
                "reference_bases": var_entry["reference_bases"],
                "alternate_bases": var_entry["alternate_bases"],
            }
            region_results.append([
                _MockAd(
                    obs=scorer.get("obs", {}),
                    var=scorer.get("var", {}),
                    X=scorer["X"],
                    quantiles=scorer.get("quantiles"),
                    variant=variant_meta,
                )
                for scorer in var_entry["scorers"]
            ])
        regions.append((region_data["name"], region_results))
    return regions, data.get("scorers", [])


__version__ = "0.6.1"

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

    if args.test_fixture:
        import json
        with open(args.test_fixture) as f:
            fixture_data = json.load(f)
        with open(args.output, "w", newline="") as outfile:
            writer = csv.writer(outfile, delimiter="\t")
            writer.writerow(fixture_data["columns"])
            for row in fixture_data["rows"]:
                writer.writerow(row)
        logging.info("Fixture mode: wrote %d rows to %s", len(fixture_data["rows"]), args.output)
        return

    if args.mock_ism_results:
        mock_regions, mock_scorers = _load_mock_ism_results(args.mock_ism_results)
        with open(args.output, "w", newline="") as outfile:
            writer = csv.writer(outfile, delimiter="\t")
            writer.writerow([
                "region", "position", "ref_base", "alt_base",
                "gene_id", "gene_name", "gene_type",
                "scorer", "track_name", "ontology_curie",
                "raw_score", "quantile_score",
            ])
            row_count = 0
            for name, results in mock_regions:
                row_count += _write_region_results(outfile, name, results, mock_scorers)
        logging.info("Mock mode: wrote %d rows to %s", row_count, args.output)
        return

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

                # results is list[list[AnnData]] -- outer=variants (3*width), inner=scorers
                row_count += _write_region_results(outfile, name, results, args.scorers)

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
    parser.add_argument("--test-fixture", default=None,
                        help="Test fixture JSON for CI testing (bypasses API)")
    parser.add_argument("--mock-ism-results", default=None,
                        help="Mock AnnData JSON for CI testing the post-processing "
                             "path (bypasses API but exercises the vectorized loop)")
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
