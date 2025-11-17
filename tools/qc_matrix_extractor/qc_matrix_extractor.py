#!/usr/bin/env python

import argparse
import json
import os
import sys

import pandas as pd


# Logging function #
def log(msg, log_file=None):
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


# CLI parsing #
def parse_args():
    parser = argparse.ArgumentParser(
        description="Galaxy-safe QC matrix extractor"
    )
    parser.add_argument("--input", "-i", required=True, help="MultiQC tabular stats file (TSV)")
    parser.add_argument("--attributes", "-a", default="", help="Comma-separated list of columns/metrics to extract")
    parser.add_argument("--thresholds", "-t", default="{}", help='JSON string with QC thresholds, e.g. {"Total_reads":1000000,"Coverage_gte_10x_pct":90,"Contam_pct":5}')
    parser.add_argument("--round", "-r", type=int, default=2, help="Number of decimals for numeric rounding")
    parser.add_argument("--derive_reads", action="store_true", help="Calculate derived read metrics from Kraken percentages")
    parser.add_argument("--outdir", "-o", default=".", help="Output directory")
    parser.add_argument("--list", action="store_true", help="List all available columns in the input file and exit")
    return parser.parse_args()


# QC decision function #
def qc_decision(row, thresholds):
    for metric, thresh in thresholds.items():
        if metric == "Contam_pct":
            if row.get(metric, 100) > thresh:
                return "Fail"
        else:
            if row.get(metric, 0) < thresh:
                return "Fail"
    return "Pass"


# Main #
def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    log_file = os.path.join(args.outdir, "qc_tool.log")

    try:
        df = pd.read_csv(args.input, sep="\t")
    except Exception as e:
        print(f"ERROR: Cannot read input file: {e}", file=sys.stderr)
        sys.exit(1)

    # List columns and exit #
    if args.list:
        print("\nAvailable columns in the input file:\n")
        for col in df.columns:
            print(f" - {col}")
        sys.exit(0)

    # Parse thresholds JSON #
    try:
        QC_THRESHOLDS = json.loads(args.thresholds)
    except json.JSONDecodeError as e:
        print(f"ERROR: Cannot parse thresholds JSON: {e}", file=sys.stderr)
        sys.exit(1)

    # Rename MultiQC columns to Galaxy-safe names #
    rename_map = {
        "qualimap_bamqc-total_reads": "Total_reads",
        "qualimap_bamqc-mapped_reads": "Mapped_reads",
        "qualimap_bamqc-percentage_aligned": "Mapping_pct",
        "qualimap_bamqc-median_coverage": "Median_depth",
        "qualimap_bamqc-10_x_pc": "Coverage_gte_10x_pct",
        "qualimap_bamqc-avg_gc": "GC_pct",
        "kraken-pct_top_one": "Kraken_top1_pct",
        "kraken-pct_unclassified": "Kraken_unclassified_pct"
    }
    df.rename(columns=rename_map, inplace=True)

    # Calculate contamination automatically #
    if "Kraken_top1_pct" in df.columns:
        df["Contam_pct"] = 100 - df["Kraken_top1_pct"]

    # Automatic read count scaling #
    for col in ["Total_reads", "Mapped_reads"]:
        if col in df.columns:
            max_val = df[col].max()
            if max_val < 1e5:
                log(f"{col} appears to be in millions — scaling by 1e6", log_file)
                df[col] = df[col] * 1e6
            else:
                log(f"{col} appears to be in raw counts — no scaling applied", log_file)

    # Optional derived reads #
    if args.derive_reads:
        if "Total_reads" in df.columns and "Kraken_top1_pct" in df.columns:
            df["MTB_reads"] = df["Total_reads"] * df["Kraken_top1_pct"] / 100
        if "Total_reads" in df.columns and "Kraken_unclassified_pct" in df.columns:
            df["Unclassified_reads"] = df["Total_reads"] * df["Kraken_unclassified_pct"] / 100

    # QC decisions #
    if QC_THRESHOLDS:
        df["QC_status"] = df.apply(lambda row: qc_decision(row, QC_THRESHOLDS), axis=1)
        for metric, thresh in QC_THRESHOLDS.items():
            col_name = f"{metric}_pass"
            if metric == "Contam_pct":
                df[col_name] = df.get(metric, 100) <= thresh
            else:
                df[col_name] = df.get(metric, 0) >= thresh

    # Columns selection #
    ALLOWED_COLUMNS = [
        "Sample", "Total_reads", "Mapped_reads", "Mapping_pct", "Median_depth",
        "Coverage_gte_10x_pct", "GC_pct", "Kraken_top1_pct", "Kraken_unclassified_pct",
        "Contam_pct", "QC_status", "Total_reads_pass", "Coverage_gte_10x_pct_pass",
        "Contam_pct_pass", "MTB_reads", "Unclassified_reads"
    ]

    if args.attributes and args.attributes.lower() != "none":
        ATTRIBUTES_OF_INTEREST = [c.strip() for c in args.attributes.split(",") if c.strip()]
        missing = [c for c in ATTRIBUTES_OF_INTEREST if c not in df.columns]
        if missing:
            log(f"WARNING: The following requested columns are missing: {missing}", log_file)
        cols_to_keep = [c for c in ATTRIBUTES_OF_INTEREST if c in df.columns]
    else:
        cols_to_keep = [c for c in ALLOWED_COLUMNS if c in df.columns]
        log("No specific columns selected — including all allowed columns.", log_file)

    df = df[cols_to_keep]

    # Round numeric columns #
    numeric_cols = df.select_dtypes(include="number").columns
    df[numeric_cols] = df[numeric_cols].round(args.round)

    # Save outputs #
    tsv_out = os.path.join(args.outdir, "QC_matrix_multiqc.tsv")
    csv_out = os.path.join(args.outdir, "QC_matrix_multiqc.csv")
    df.to_csv(tsv_out, sep="\t", index=False)
    df.to_csv(csv_out, sep=",", index=False, encoding="utf-8-sig")
    log("QC matrix saved successfully!", log_file)
    log(f"TSV: {tsv_out}", log_file)
    log(f"CSV: {csv_out}", log_file)


if __name__ == "__main__":
    main()
