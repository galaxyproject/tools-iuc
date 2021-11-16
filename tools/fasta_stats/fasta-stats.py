#!/usr/bin/env python


# python version of fasta-stats with some extra features
# written by anmol.kiran@gmail.com
# git: @codemeleon
# date: 10/11/2021


import click
import numpy as np
from Bio import SeqIO
import re
from os import path


@click.command()
@click.option(
    "--fasta",
    help="Scaffold File",
    type=str,
    default="./post_scaffold_newref.fa",
    show_default=True,
)
@click.option(
    "--stats",
    help="General stats csv output file. If not given, stdout",
    type=str,
    default=None,
    show_default=True,
)
@click.option(
    "--gap",
    help="Gap stats bed file. If not given, stdout",
    type=str,
    default=None,
    show_default=True,
)
def run(fasta, stats, gap):
    """Generates scaffold statistics."""
    if not fasta:
        exit("Input file not given.")
    if not path.isfile(fasta):
        exit(f"{fasta} path does not exist.")

    seq_len = {}
    bases_global = {"A": 0, "N": 0, "T": 0, "C": 0, "G": 0}
    bases_seq = {}
    seq_id_Ngaprange = {}
    nstart = 0
    contigs_len = []
    gap_count = 0
    for seq_record in SeqIO.parse(fasta, "fasta"):
        seq = str(seq_record.seq).upper()
        # print(len(seq))
        seq_len[seq_record.id] = len(seq)

        # NOTE: Nucleotide count
        bases_seq[seq_record.id] = {
            "A": seq.count("A"),
            "N": seq.count("N"),
            "T": seq.count("T"),
            "C": seq.count("C"),
            "G": seq.count("G"),
        }
        bases_global["A"] += bases_seq[seq_record.id]["A"]
        bases_global["N"] += bases_seq[seq_record.id]["N"]
        bases_global["T"] += bases_seq[seq_record.id]["T"]
        bases_global["C"] += bases_seq[seq_record.id]["C"]
        bases_global["G"] += bases_seq[seq_record.id]["G"]

        # NOTE: Gap count and their range
        range_gen = re.finditer("N+", seq)
        n_range = [match.span() for match in range_gen]
        for n_rng in n_range:
            if n_rng[0] == 0 or n_rng[1] == seq_len[seq_record.id]:
                continue
            else:
                gap_count += 1

        # NOTE: Contigs, their lenths from scaffold and their N gap range
        seq_id_Ngaprange[seq_record.id] = n_range
        n_range_len = len(n_range)
        if n_range_len > 0:
            n_range = (
                [(0, 0)] + n_range +
                [(seq_len[seq_record.id], seq_len[seq_record.id])]
            )
            for idx in range(n_range_len + 1):
                nstart = n_range[idx][1]
                nend = n_range[idx + 1][0]
                con_len = nend - nstart
                if con_len:
                    contigs_len.append(con_len)
        else:
            contigs_len.append(len(seq))

    seq_len_list = list(seq_len.values())
    scaffold_lens = list(seq_len_list)

    # NOTE: Scaffold statistics
    scaffold_lens.sort(reverse=True)
    scaffold_lens = np.array(scaffold_lens)
    scaffold_lens = np.cumsum(scaffold_lens)

    n50_len = scaffold_lens[-1] * 0.5
    n50_idx = np.where(scaffold_lens > n50_len)[0][0]
    n90_len = scaffold_lens[-1] * 0.9
    n90_idx = np.where(scaffold_lens > n90_len)[0][0]

    to_csv = []
    to_csv.append(["Fields", "Values"])
    to_csv.append(["Scaffold L50", f"{n50_idx + 1}"])
    to_csv.append(["Scaffold N50", f"{seq_len_list[n50_idx]}"])
    to_csv.append(["Scaffold L90", f"{n90_idx + 1}"])
    to_csv.append(["Scaffold N90", f"{seq_len_list[n90_idx]}"])
    to_csv.append(["Scaffold len_max", f"{np.max(seq_len_list)}"])
    to_csv.append(["Scaffold len_min", f"{np.min(seq_len_list)}"])
    to_csv.append(["Scaffold len_mean", f"{int(np.mean(seq_len_list))}"])
    to_csv.append(["Scaffold len_median", f"{int(np.median(seq_len_list))}"])
    to_csv.append(["Scaffold len_std", f"{int(np.std(seq_len_list))}"])
    to_csv.append(["Scaffold num_N", f'{bases_global["N"]}'])
    to_csv.append(["Scaffold num_A", f'{bases_global["A"]}'])
    to_csv.append(["Scaffold num_T", f'{bases_global["T"]}'])
    to_csv.append(["Scaffold num_C", f'{bases_global["C"]}'])
    to_csv.append(["Scaffold num_bp", f"{scaffold_lens[-1]}"])
    to_csv.append(["Scaffold num_bp_not_N",
                  f'{scaffold_lens[-1] - bases_global["N"]}'])
    to_csv.append(["Scaffold num_seq", f"{len(seq_len_list)}"])
    to_csv.append(
        [
            "Scaffold GC content overall",
            "%0.2f"
            % ((bases_global["G"] + bases_global["C"]) * 100.0 / sum(seq_len.values())),
        ]
    )
    to_csv.append(["Number of gaps", f"{gap_count}"])

    # NOTE: Contig statistics
    seq_len_list = list(contigs_len)
    contigs_len.sort(reverse=True)
    contigs_len = np.array(contigs_len)
    contigs_len = np.cumsum(contigs_len)
    n50_len = contigs_len[-1] * 0.5
    n50_idx = np.where(contigs_len > n50_len)[0][0]
    n90_len = contigs_len[-1] * 0.9
    n90_idx = np.where(contigs_len > n90_len)[0][0]
    to_csv.append(["Contig L50", f"{n50_idx + 1}"])
    to_csv.append(["Contig N50", f"{seq_len_list[n50_idx]}"])
    to_csv.append(["Contig L90", f"{n90_idx + 1}"])
    to_csv.append(["Contig N90", f"{seq_len_list[n90_idx]}"])
    to_csv.append(["Contig len_max", f"{np.max(seq_len_list)}"])
    to_csv.append(["Contig len_min", f"{np.min(seq_len_list)}"])
    to_csv.append(["Contig len_mean", f"{int(np.mean(seq_len_list))}"])
    to_csv.append(["Contig len_median", f"{int(np.median(seq_len_list))}"])
    to_csv.append(["Contig len_std", f"{int(np.std(seq_len_list))}"])
    to_csv.append(["Contig num_bp", f"{contigs_len[-1]}"])
    to_csv.append(["Contig num_seq", f"{len(contigs_len)}"])

    if stats:
        with open(stats, "w") as fout:
            for val in to_csv:
                print("\t".join(val), file=fout)
    else:
        print("General scaffolds and contigs stats.\n")
        for val in to_csv:
            print("\t".join(val))
    gaps_csv = [["Scaffold", "start", "end"]]
    for key in seq_id_Ngaprange:
        if not len(seq_id_Ngaprange[key]):
            continue
        for rng in seq_id_Ngaprange[key]:
            gaps_csv.append([f"{key}", f"{rng[0]}", f"{rng[1]}"])
    if gap:
        with open(gap, "w") as fout:
            for gp in gaps_csv:
                print("\t".join(gp), file=fout)

    else:
        print("\n\nGaps in scaffolds.\n")
        for gp in gaps_csv:
            print("\t".join(gp))


if __name__ == "__main__":
    run()
