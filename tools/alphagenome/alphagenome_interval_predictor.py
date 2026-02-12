#!/usr/bin/env python
"""
AlphaGenome Interval Predictor for Galaxy

Predicts regulatory tracks for genomic intervals — no variants, just baseline
characterization of the chromatin/expression landscape.
"""

import argparse
import csv
import logging
import os
import sys

import numpy as np
from alphagenome.data import genome
from alphagenome.models import dna_client

__version__ = "0.1.0"

OUTPUT_TYPE_MAP = {
    "RNA_SEQ": dna_client.OutputType.RNA_SEQ,
    "ATAC": dna_client.OutputType.ATAC,
    "CAGE": dna_client.OutputType.CAGE,
    "DNASE": dna_client.OutputType.DNASE,
    "CHIP_HISTONE": dna_client.OutputType.CHIP_HISTONE,
    "CHIP_TF": dna_client.OutputType.CHIP_TF,
    "SPLICE_SITES": dna_client.OutputType.SPLICE_SITES,
    "PROCAP": dna_client.OutputType.PROCAP,
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
    if local_model:
        from alphagenome_research.model import dna_model
        return dna_model.create_from_huggingface("all_folds")
    return dna_client.create(api_key)


def parse_bed(bed_path, max_intervals):
    """Parse BED file and return list of (chrom, start, end, name) tuples."""
    intervals = []
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
            if len(intervals) >= max_intervals:
                logging.warning("Reached max intervals (%d), skipping remaining", max_intervals)
                break
            intervals.append((chrom, start, end, name))
    return intervals


def extract_region_slice(values, interval_start, region_start, region_end):
    """Slice the prediction array to only cover the original BED region.

    The prediction covers the full resized interval; we want stats only for the
    user's original region within it.

    values: numpy array of shape (seq_length, num_tracks)
    interval_start: 0-based start of the resized prediction interval
    region_start: 0-based start of the user's BED region
    region_end: 0-based end of the user's BED region
    """
    offset_start = region_start - interval_start
    offset_end = region_end - interval_start
    # Clamp to valid bounds
    offset_start = max(0, offset_start)
    offset_end = min(values.shape[0], offset_end)
    return values[offset_start:offset_end]


def run(args):
    logging.info("AlphaGenome Interval Predictor v%s", __version__)
    logging.info("Input: %s", args.input)
    logging.info("Output types: %s", ", ".join(args.output_types))
    logging.info("Output mode: %s", args.output_mode)
    logging.info("Organism: %s", args.organism)
    logging.info("Sequence length: %s", args.sequence_length)
    logging.info("Max intervals: %d", args.max_intervals)

    api_key = args.api_key or os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key and not args.local_model:
        logging.error("No API key provided. Set ALPHAGENOME_API_KEY or use --api-key")
        sys.exit(1)

    organism = ORGANISM_MAP[args.organism]
    seq_length = SEQUENCE_LENGTH_MAP[args.sequence_length]
    requested_outputs = [OUTPUT_TYPE_MAP[t] for t in args.output_types]
    ontology_terms = []
    if args.ontology_terms:
        ontology_terms = [t.strip() for t in args.ontology_terms.split(",") if t.strip()]

    intervals = parse_bed(args.input, args.max_intervals)
    if not intervals:
        logging.error("No valid intervals found in input BED file")
        sys.exit(1)
    logging.info("Loaded %d intervals", len(intervals))

    logging.info("Connecting to AlphaGenome...")
    model = create_model(api_key, local_model=args.local_model)
    logging.info("Model ready.")

    stats = {"total": 0, "predicted": 0, "errors": 0}

    with open(args.output, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t")

        if args.output_mode == "summary":
            writer.writerow([
                "chrom", "start", "end", "name", "output_type",
                "track_name", "ontology_curie", "mean_signal", "max_signal",
            ])
        else:
            writer.writerow([
                "chrom", "bin_start", "bin_end", "region_name", "output_type",
                "track_name", "ontology_curie", "mean_signal",
            ])

        for interval_num, (chrom, start, end, name) in enumerate(intervals):
            stats["total"] += 1

            if interval_num > 0 and interval_num % 10 == 0:
                logging.info(
                    "Progress: %d/%d intervals (%d predicted, %d errors)",
                    interval_num, len(intervals), stats["predicted"], stats["errors"],
                )

            try:
                interval = genome.Interval(chrom, start, end).resize(seq_length)

                output = model.predict_interval(
                    interval, organism=organism,
                    requested_outputs=requested_outputs,
                    ontology_terms=ontology_terms,
                )

                for otype in args.output_types:
                    attr_name = otype.lower()
                    track_data = getattr(output, attr_name, None)
                    if track_data is None:
                        logging.warning("No %s data for %s", otype, name)
                        continue

                    values = track_data.values  # (seq_length, num_tracks)
                    metadata = track_data.metadata  # DataFrame with track info

                    # Slice to the original BED region
                    region_values = extract_region_slice(
                        values, interval.start, start, end,
                    )

                    num_tracks = region_values.shape[1] if region_values.ndim > 1 else 1
                    if region_values.ndim == 1:
                        region_values = region_values.reshape(-1, 1)

                    for track_idx in range(num_tracks):
                        track_vals = region_values[:, track_idx]
                        track_name = ""
                        ontology_curie = ""
                        if metadata is not None and len(metadata) > track_idx:
                            row = metadata.iloc[track_idx]
                            track_name = str(row.get("track_name", ""))
                            ontology_curie = str(row.get("ontology_curie", ""))

                        if args.output_mode == "summary":
                            mean_sig = float(np.mean(track_vals))
                            max_sig = float(np.max(track_vals))
                            writer.writerow([
                                chrom, start, end, name, otype,
                                track_name, ontology_curie,
                                f"{mean_sig:.6f}", f"{max_sig:.6f}",
                            ])
                        else:
                            # Binned mode
                            region_len = region_values.shape[0]
                            bin_size = args.bin_size
                            for bin_start_offset in range(0, region_len, bin_size):
                                bin_end_offset = min(bin_start_offset + bin_size, region_len)
                                bin_vals = track_vals[bin_start_offset:bin_end_offset]
                                mean_sig = float(np.mean(bin_vals))
                                writer.writerow([
                                    chrom, start + bin_start_offset,
                                    start + bin_end_offset, name, otype,
                                    track_name, ontology_curie,
                                    f"{mean_sig:.6f}",
                                ])

                stats["predicted"] += 1

            except Exception as e:
                logging.error("Error predicting %s (%s:%d-%d): %s", name, chrom, start, end, e)
                stats["errors"] += 1

    logging.info("=" * 50)
    logging.info("DONE — %d total, %d predicted, %d errors",
                 stats["total"], stats["predicted"], stats["errors"])
    logging.info("Output: %s", args.output)

    if stats["errors"] > 0 and stats["predicted"] == 0:
        logging.error("All intervals failed. Check API key and network connectivity.")
        sys.exit(1)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Predict regulatory tracks for genomic intervals using AlphaGenome",
    )
    parser.add_argument("--input", required=True, help="Input BED file")
    parser.add_argument("--output", required=True, help="Output TSV file")
    parser.add_argument("--api-key", default=None, help="AlphaGenome API key (or set ALPHAGENOME_API_KEY)")
    parser.add_argument(
        "--organism", choices=["human", "mouse"], default="human",
    )
    parser.add_argument(
        "--output-types", nargs="+", choices=list(OUTPUT_TYPE_MAP.keys()),
        default=["RNA_SEQ"],
    )
    parser.add_argument("--ontology-terms", default=None)
    parser.add_argument(
        "--sequence-length", choices=list(SEQUENCE_LENGTH_MAP.keys()), default="1MB",
    )
    parser.add_argument("--max-intervals", type=int, default=50)
    parser.add_argument(
        "--output-mode", choices=["summary", "binned"], default="summary",
    )
    parser.add_argument("--bin-size", type=int, default=128)
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
