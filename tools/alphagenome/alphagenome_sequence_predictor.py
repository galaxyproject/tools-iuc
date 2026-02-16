#!/usr/bin/env python
"""
AlphaGenome Sequence Predictor for Galaxy

Predicts regulatory tracks from raw DNA sequence — no genomic coordinates needed.
For synthetic biology (designed sequences) and non-reference assemblies.
"""

import argparse
import csv
import logging
import os
import sys

import numpy as np
from alphagenome.models import dna_client

__version__ = "0.5.1"

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


def parse_fasta(fasta_path, max_sequences):
    sequences = []
    current_id = None
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    sequences.append((current_id, "".join(current_seq)))
                    if len(sequences) >= max_sequences:
                        logging.warning("Reached max sequences (%d), skipping remaining", max_sequences)
                        current_id = None
                        break
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id is not None:
        sequences.append((current_id, "".join(current_seq)))

    return sequences


def prepare_sequence(seq, target_length):
    """Pad (N-centered) or center-trim to target_length. Returns (seq, content_start, content_end)."""
    seq = seq.upper()
    if len(seq) == target_length:
        return seq, 0, len(seq)
    elif len(seq) < target_length:
        pad_total = target_length - len(seq)
        pad_left = pad_total // 2
        prepared = "N" * pad_left + seq + "N" * (pad_total - pad_left)
        return prepared, pad_left, pad_left + len(seq)
    else:
        trim_start = (len(seq) - target_length) // 2
        prepared = seq[trim_start:trim_start + target_length]
        return prepared, 0, target_length


def run(args):
    logging.info("AlphaGenome Sequence Predictor v%s", __version__)
    logging.info("Input: %s", args.input)
    logging.info("Output types: %s", ", ".join(args.output_types))
    logging.info("Output mode: %s", args.output_mode)
    logging.info("Organism: %s", args.organism)
    logging.info("Sequence length: %s", args.sequence_length)
    logging.info("Max sequences: %d", args.max_sequences)

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

    api_key = args.api_key or os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key and not args.local_model:
        logging.error("No API key provided. Set ALPHAGENOME_API_KEY or use --api-key")
        sys.exit(1)

    organism = ORGANISM_MAP[args.organism]
    target_length = SEQUENCE_LENGTH_MAP[args.sequence_length]
    requested_outputs = [OUTPUT_TYPE_MAP[t] for t in args.output_types]
    ontology_terms = []
    if args.ontology_terms:
        ontology_terms = [t.strip() for t in args.ontology_terms.split(",") if t.strip()]

    sequences = parse_fasta(args.input, args.max_sequences)
    if not sequences:
        logging.error("No valid sequences found in input FASTA file")
        sys.exit(1)
    logging.info("Loaded %d sequences", len(sequences))

    logging.info("Connecting to AlphaGenome...")
    model = create_model(api_key, local_model=args.local_model)
    logging.info("Model ready.")

    stats = {"total": 0, "predicted": 0, "errors": 0}

    with open(args.output, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t")

        if args.output_mode == "summary":
            writer.writerow([
                "sequence_id", "sequence_length", "output_type",
                "track_name", "ontology_curie", "mean_signal", "max_signal",
            ])
        else:
            writer.writerow([
                "sequence_id", "bin_start", "bin_end", "output_type",
                "track_name", "ontology_curie", "mean_signal",
            ])

        for seq_num, (seq_id, raw_seq) in enumerate(sequences):
            stats["total"] += 1
            orig_length = len(raw_seq)

            if seq_num > 0 and seq_num % 5 == 0:
                logging.info(
                    "Progress: %d/%d sequences (%d predicted, %d errors)",
                    seq_num, len(sequences), stats["predicted"], stats["errors"],
                )

            logging.info("Sequence %d/%d: %s (%dbp)", seq_num + 1, len(sequences), seq_id, orig_length)

            try:
                prepared_seq, content_start, content_end = prepare_sequence(raw_seq, target_length)

                if orig_length != target_length:
                    if orig_length < target_length:
                        logging.info("  N-padded %dbp -> %dbp", orig_length, target_length)
                    else:
                        logging.info("  Center-trimmed %dbp -> %dbp", orig_length, target_length)

                output = model.predict_sequence(
                    prepared_seq, organism=organism,
                    requested_outputs=requested_outputs,
                    ontology_terms=ontology_terms,
                )

                for otype in args.output_types:
                    attr_name = otype.lower()
                    track_data = getattr(output, attr_name, None)
                    if track_data is None:
                        logging.warning("No %s data for %s", otype, seq_id)
                        continue

                    values = track_data.values
                    metadata = track_data.metadata

                    # Slice to the actual content region (non-N portion)
                    content_values = values[content_start:content_end]

                    num_tracks = content_values.shape[1] if content_values.ndim > 1 else 1
                    if content_values.ndim == 1:
                        content_values = content_values.reshape(-1, 1)

                    for track_idx in range(num_tracks):
                        track_vals = content_values[:, track_idx]
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
                                seq_id, orig_length, otype,
                                track_name, ontology_curie,
                                f"{mean_sig:.6f}", f"{max_sig:.6f}",
                            ])
                        else:
                            content_len = content_values.shape[0]
                            bin_size = args.bin_size
                            for bin_start in range(0, content_len, bin_size):
                                bin_end = min(bin_start + bin_size, content_len)
                                bin_vals = track_vals[bin_start:bin_end]
                                mean_sig = float(np.mean(bin_vals))
                                writer.writerow([
                                    seq_id, bin_start, bin_end, otype,
                                    track_name, ontology_curie,
                                    f"{mean_sig:.6f}",
                                ])

                stats["predicted"] += 1

            except Exception as e:
                logging.error("Error predicting %s: %s", seq_id, e)
                stats["errors"] += 1

    logging.info("=" * 50)
    logging.info("DONE — %d total, %d predicted, %d errors",
                 stats["total"], stats["predicted"], stats["errors"])
    logging.info("Output: %s", args.output)

    if stats["errors"] > 0 and stats["predicted"] == 0:
        logging.error("All sequences failed. Check API key and network connectivity.")
        sys.exit(1)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Predict regulatory tracks from DNA sequence using AlphaGenome",
    )
    parser.add_argument("--input", required=True, help="Input FASTA file")
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
        "--sequence-length", choices=list(SEQUENCE_LENGTH_MAP.keys()), default="16KB",
    )
    parser.add_argument("--max-sequences", type=int, default=20)
    parser.add_argument(
        "--output-mode", choices=["summary", "binned"], default="summary",
    )
    parser.add_argument("--bin-size", type=int, default=128)
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
