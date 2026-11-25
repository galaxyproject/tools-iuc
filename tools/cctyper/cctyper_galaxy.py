#!/usr/bin/env python
"""Run CRISPRCasTyper (cctyper) and normalise its output directory for Galaxy.

``cctyper`` writes a directory of files and only creates a file when there is
data for it (for example ``CRISPR_Cas.tab`` is not written at all when no
CRISPR-Cas locus is found). Galaxy outputs on the other hand have to exist
after every successful run, otherwise the job fails. This wrapper therefore

* runs ``cctyper`` with the arguments handed over by the Galaxy command line,
* writes header-only versions of the three main tables (and an empty GFF3) when
  ``cctyper`` did not produce them,
* writes a placeholder SVG when plotting was requested but there was nothing
  to draw,
* moves every other file into ``results/extras/`` where Galaxy picks it up as
  a collection.

The Cas HMM profiles and the other reference data are shipped inside the cctyper
package, so no database download or configuration is needed.

It also works around two limits of ``cctyper`` 1.9.0 that matter on Galaxy:

* an input without a single Cas HMM hit aborts with a ``KeyError``. That is a
  valid input (no CRISPR-Cas system), so the wrapper writes the empty result
  tables, keeps the traceback in ``extras/cctyper_stderr.log`` and reports a
  warning instead of failing the job.
* pyrodigal cannot train on sequences shorter than 20 kb, which makes the
  default single mode fail on short contigs. With ``--prodigal-mode auto`` the
  metagenomic mode is used when no sequence reaches that length.
"""

import argparse
import gzip
import logging
import os
import shutil
import subprocess
import sys

__version__ = "1.9.0"

# Columns written by cctyper 1.9.0, used for header-only placeholder tables.
CRISPR_CAS_COLUMNS = (
    "Contig", "Operon", "Operon_Pos", "Prediction", "CRISPRs", "Distances",
    "Prediction_Cas", "Prediction_CRISPRs",
)
CAS_OPERON_COLUMNS = (
    "Contig", "Operon", "Start", "End", "Prediction", "Complete_Interference",
    "Complete_Adaptation", "Best_type", "Best_score", "Genes", "Positions",
    "E-values", "CoverageSeq", "CoverageHMM", "Strand_Interference",
    "Strand_Adaptation",
)
CRISPR_COLUMNS = (
    "Contig", "CRISPR", "Start", "End", "Consensus_repeat", "N_repeats",
    "Repeat_len", "Spacer_len_avg", "Repeat_identity", "Spacer_identity",
    "Spacer_len_sem", "Trusted", "Prediction", "Subtype",
    "Subtype_probability",
)

# Files copied to a fixed name in results/; everything else goes to extras/
MAIN_OUTPUTS = {
    "CRISPR_Cas.tab": "CRISPR_Cas.tab",
    "cas_operons.tab": "cas_operons.tab",
    "crisprs_all.tab": "crisprs_all.tab",
    "crisprs.gff": "crisprs.gff3",
    "plot.svg": "plot.svg",
}

# cctyper 1.9.0 aborts like this when the HMM search returns no hit at all
NO_CAS_HIT_ERROR = "['Hmm', 'ORF'] not in index"
SHORT_SEQUENCE_ERROR = "sequence must be at least 20000 characters"
PYRODIGAL_MIN_LENGTH = 20000

# cctyper calls sys.exit() without an exit code on several error paths (missing
# database files, unreadable input, ...), so a zero exit code alone does not mean
# the run completed. The last thing a successful run logs is this message.
COMPLETION_MARKER = "Removing temporary files"
MISSING_DB_ERROR = "Could not find HMM profile database"

PLACEHOLDER_SVG = """<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="640" height="60" version="1.1">
  <text x="12" y="36" font-family="sans-serif" font-size="14">\
No CRISPR-Cas locus or CRISPR array was detected, there is no map to draw.</text>
</svg>
"""


def setup_logging():
    logging.basicConfig(
        format="%(asctime)s %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.INFO,
        stream=sys.stderr,
    )


def is_gzip(path):
    with open(path, "rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def unpack_input(path, workdir):
    """Return a plain FASTA path, unpacking gzip input if needed."""
    if not is_gzip(path):
        return path
    target = os.path.join(workdir, "input.fasta")
    logging.info("Decompressing gzip input to %s", target)
    with gzip.open(path, "rt") as source, open(target, "w") as dest:
        shutil.copyfileobj(source, dest)
    return target


def longest_sequence(path):
    """Length of the longest sequence in a (optionally gzipped) FASTA file."""
    longest = current = 0
    opener = gzip.open if is_gzip(path) else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                longest = max(longest, current)
                current = 0
            else:
                current += len(line.strip())
    return max(longest, current)


def resolve_prodigal_mode(requested, fasta):
    if requested != "auto":
        return requested
    length = longest_sequence(fasta)
    if length < PYRODIGAL_MIN_LENGTH:
        logging.info(
            "Longest sequence is %d bp (< %d bp), pyrodigal cannot train on that: "
            "using the metagenomic gene prediction mode",
            length, PYRODIGAL_MIN_LENGTH)
        return "meta"
    logging.info("Longest sequence is %d bp: using the single genome gene prediction mode", length)
    return "single"


def run_cctyper(command):
    """Run cctyper, forward its log and return (returncode, captured stderr)."""
    captured = []
    logging.info("Running: %s", " ".join(command))
    process = subprocess.Popen(
        command, stdout=sys.stderr, stderr=subprocess.PIPE, text=True, bufsize=1,
    )
    for line in process.stderr:
        captured.append(line)
        sys.stderr.write(line)
        sys.stderr.flush()
    process.wait()
    return process.returncode, "".join(captured)


def write_table(path, columns):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        logging.warning("No data for %s, writing a header-only table", os.path.basename(path))
        with open(path, "w") as handle:
            handle.write("\t".join(columns) + "\n")


def write_gff(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        logging.warning("No CRISPR array features, writing an empty GFF3 file")
        with open(path, "w") as handle:
            handle.write("##gff-version 3\n")


def write_plot_placeholder(path):
    if not os.path.exists(path):
        logging.warning("Nothing to plot, writing a placeholder map")
        with open(path, "w") as handle:
            handle.write(PLACEHOLDER_SVG)


def count_rows(path):
    if not os.path.exists(path):
        return 0
    with open(path) as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def collect_outputs(cctyper_dir, results_dir, want_plot, stderr_text):
    """Move cctyper's output into the fixed Galaxy layout."""
    extras = os.path.join(results_dir, "extras")
    os.makedirs(extras, exist_ok=True)

    for source_name, target_name in MAIN_OUTPUTS.items():
        source = os.path.join(cctyper_dir, source_name)
        if os.path.exists(source):
            shutil.move(source, os.path.join(results_dir, target_name))

    # Everything left over (orphan/putative tables, genes, hmmer hits, spacers,
    # arguments, PNG plot, ...) is flattened into results/extras/
    for root, _dirs, files in os.walk(cctyper_dir):
        for name in files:
            source = os.path.join(root, name)
            target = os.path.join(extras, name)
            if os.path.exists(target):
                logging.warning("Skipping duplicate output file %s", name)
                continue
            shutil.move(source, target)

    if stderr_text:
        with open(os.path.join(extras, "cctyper_stderr.log"), "w") as handle:
            handle.write(stderr_text)

    write_table(os.path.join(results_dir, "CRISPR_Cas.tab"), CRISPR_CAS_COLUMNS)
    write_table(os.path.join(results_dir, "cas_operons.tab"), CAS_OPERON_COLUMNS)
    write_table(os.path.join(results_dir, "crisprs_all.tab"), CRISPR_COLUMNS)
    write_gff(os.path.join(results_dir, "crisprs.gff3"))
    if want_plot:
        write_plot_placeholder(os.path.join(results_dir, "plot.svg"))

    shutil.rmtree(cctyper_dir, ignore_errors=True)

    logging.info(
        "Result: %d CRISPR-Cas loci, %d Cas operons, %d CRISPR arrays",
        count_rows(os.path.join(results_dir, "CRISPR_Cas.tab")),
        count_rows(os.path.join(results_dir, "cas_operons.tab")),
        count_rows(os.path.join(results_dir, "crisprs_all.tab")),
    )


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(
        description="Galaxy wrapper for CRISPRCasTyper: run cctyper and normalise its output.",
    )
    parser.add_argument("--input", required=True, help="Nucleotide FASTA (optionally gzipped)")
    parser.add_argument("--prodigal-mode", default="auto", choices=["auto", "single", "meta"],
                        help="Gene prediction mode, 'auto' picks 'meta' for sequences shorter than 20 kb")
    parser.add_argument("--workdir", default="cctyper_out", help="Directory cctyper writes to")
    parser.add_argument("--results", default="results", help="Directory with the Galaxy outputs")
    parser.add_argument("--plot", action="store_true",
                        help="A map was requested, write a placeholder if nothing was drawn")
    parser.add_argument("cctyper_arguments", nargs=argparse.REMAINDER,
                        help="Arguments for cctyper, which have to be preceded by a '--' "
                             "separator because they use the same option names")
    args = parser.parse_args(argv)
    # cctyper uses the same option prefixes, so the command line separates the
    # two sets of arguments with '--'. argparse keeps that separator in the
    # remainder and cctyper would reject it.
    if args.cctyper_arguments[:1] == ["--"]:
        args.cctyper_arguments = args.cctyper_arguments[1:]
    return args


def main(argv=None):
    setup_logging()
    args = parse_arguments(argv)

    if not os.path.isfile(args.input):
        sys.exit("Input file %s does not exist" % args.input)

    shutil.rmtree(args.workdir, ignore_errors=True)
    os.makedirs(args.results, exist_ok=True)

    fasta = unpack_input(args.input, ".")
    mode = resolve_prodigal_mode(args.prodigal_mode, fasta)

    command = ["cctyper", fasta, args.workdir, "--prodigal", mode] + args.cctyper_arguments
    returncode, stderr_text = run_cctyper(command)

    failed = returncode != 0
    no_hits = False
    if failed and NO_CAS_HIT_ERROR in stderr_text:
        logging.warning(
            "CRISPRCasTyper found no Cas HMM hit in this input. cctyper 1.9.0 stops "
            "with an internal error in that case (it cannot report an empty result), so "
            "empty result tables are written instead. Note that CRISPR arrays are only "
            "searched for near Cas operons, i.e. none were searched for here. The "
            "original error is kept in 'cctyper_stderr.log'."
        )
        failed = False
        no_hits = True
    elif failed and SHORT_SEQUENCE_ERROR in stderr_text:
        sys.exit(
            "CRISPRCasTyper failed: pyrodigal cannot train a gene model on sequences "
            "shorter than 20 kb. Choose the metagenomic gene prediction mode instead."
        )
    elif failed:
        sys.exit("CRISPRCasTyper failed with exit code %d, see the log above" % returncode)

    if not failed and not no_hits and COMPLETION_MARKER not in stderr_text:
        if MISSING_DB_ERROR in stderr_text:
            sys.exit(
                "CRISPRCasTyper could not find its Cas HMM profile database. Neither "
                "$CCTYPER_DB nor the data directory shipped with the cctyper package "
                "contain cctyper_profiles.hmm[.zst]"
            )
        sys.exit(
            "CRISPRCasTyper stopped before it finished. cctyper 1.9.0 exits with "
            "status 0 on several errors, so this is a failure even though the exit "
            "code is 0, see the log above"
        )

    collect_outputs(args.workdir, args.results, args.plot,
                    stderr_text if returncode != 0 else "")


if __name__ == "__main__":
    main()
