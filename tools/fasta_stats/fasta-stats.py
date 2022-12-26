#!/usr/bin/env python


# python version of fasta-stats with some extra features
# written by anmol.kiran@gmail.com
# git: @codemeleon
# date: 10/11/2021

import argparse
import re
from os import path

import numpy as np
from Bio import SeqIO


def calculate_NG50(estimated_genome, total_length, sequence_lengths):
    temp = 0
    teoretical_NG50 = estimated_genome / 2.0
    NG50 = 0
    for seq in sequence_lengths:
        temp += seq
        if teoretical_NG50 < temp:
            NG50 = seq
            break
    return NG50


def run(fasta, stats_output, gaps_output, genome_size):
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
                [(0, 0)] + n_range + [(seq_len[seq_record.id], seq_len[seq_record.id])]
            )
            for idx in range(n_range_len + 1):
                nstart = n_range[idx][1]
                nend = n_range[idx + 1][0]
                con_len = nend - nstart
                if con_len:
                    contigs_len.append(con_len)
        else:
            contigs_len.append(len(seq))

    # NOTE: Scaffold statistics
    SEQ_LEN_LIST = sorted(seq_len.values(), reverse=True)
    scaffold_lens = np.array(SEQ_LEN_LIST)
    scaffold_lens_sum = np.cumsum(scaffold_lens)
    N50_len = scaffold_lens_sum[-1] * 0.5
    N50_idx = np.where(scaffold_lens_sum > N50_len)[0][0]
    N90_len = scaffold_lens_sum[-1] * 0.9
    N90_idx = np.where(scaffold_lens_sum > N90_len)[0][0]
    NG50 = calculate_NG50(genome_size, scaffold_lens_sum[-1], scaffold_lens)

    # NOTE: Contig statistics
    seq_len_list = sorted(contigs_len, reverse=True)
    contigs_len = np.array(seq_len_list)
    contigs_len_sum = np.cumsum(contigs_len)
    n50_len = contigs_len_sum[-1] * 0.5
    n50_idx = np.where(contigs_len_sum > n50_len)[0][0]
    n90_len = contigs_len_sum[-1] * 0.9
    n90_idx = np.where(contigs_len_sum > n90_len)[0][0]
    ng50 = calculate_NG50(genome_size, contigs_len_sum[-1], contigs_len)

    with open(stats_output, "w") as soutput:
        soutput.write("{}\t{}\n".format("Scaffold L50", N50_idx + 1))
        soutput.write("{}\t{}\n".format("Scaffold N50", SEQ_LEN_LIST[N50_idx]))
        soutput.write("{}\t{}\n".format("Scaffold L90", N90_idx + 1))
        soutput.write("{}\t{}\n".format("Scaffold N90", SEQ_LEN_LIST[N90_idx]))
        if genome_size != 0:
            soutput.write("{}\t{}\n".format("Scaffold NG50", NG50))
        soutput.write("{}\t{}\n".format("Scaffold len_max", SEQ_LEN_LIST[0]))
        soutput.write("{}\t{}\n".format("Scaffold len_min", SEQ_LEN_LIST[-1]))
        soutput.write(
            "{}\t{}\n".format("Scaffold len_mean", int(np.mean(SEQ_LEN_LIST)))
        )
        soutput.write(
            "{}\t{}\n".format("Scaffold len_median", int(np.median(SEQ_LEN_LIST)))
        )
        soutput.write("{}\t{}\n".format("Scaffold len_std", int(np.std(SEQ_LEN_LIST))))
        soutput.write("{}\t{}\n".format("Scaffold num_A", bases_global["A"]))
        soutput.write("{}\t{}\n".format("Scaffold num_T", bases_global["T"]))
        soutput.write("{}\t{}\n".format("Scaffold num_C", bases_global["C"]))
        soutput.write("{}\t{}\n".format("Scaffold num_G", bases_global["G"]))
        soutput.write("{}\t{}\n".format("Scaffold num_N", bases_global["N"]))
        soutput.write("{}\t{}\n".format("Scaffold num_bp", scaffold_lens_sum[-1]))
        soutput.write(
            "{}\t{}\n".format(
                "Scaffold num_bp_not_N", scaffold_lens_sum[-1] - bases_global["N"]
            )
        )
        soutput.write("{}\t{}\n".format("Scaffold num_seq", len(SEQ_LEN_LIST)))
        soutput.write(
            "{}\t{:.2f}\n".format(
                "Scaffold GC content overall",
                (
                    (bases_global["G"] + bases_global["C"])
                    * 100.0
                    / scaffold_lens_sum[-1]
                ),
            )
        )

        soutput.write("{}\t{}\n".format("Contig L50", n50_idx + 1))
        soutput.write("{}\t{}\n".format("Contig N50", seq_len_list[n50_idx]))
        soutput.write("{}\t{}\n".format("Contig L90", n90_idx + 1))
        soutput.write("{}\t{}\n".format("Contig N90", seq_len_list[n90_idx]))
        if genome_size != 0:
            soutput.write("{}\t{}\n".format("Contig NG50", ng50))
        soutput.write("{}\t{}\n".format("Contig len_max", seq_len_list[0]))
        soutput.write("{}\t{}\n".format("Contig len_min", seq_len_list[-1]))
        soutput.write("{}\t{}\n".format("Contig len_mean", int(np.mean(seq_len_list))))
        soutput.write(
            "{}\t{}\n".format("Contig len_median", int(np.median(seq_len_list)))
        )
        soutput.write("{}\t{}\n".format("Contig len_std", int(np.std(seq_len_list))))
        soutput.write("{}\t{}\n".format("Contig num_bp", contigs_len_sum[-1]))
        soutput.write("{}\t{}\n".format("Contig num_seq", len(contigs_len_sum)))
        soutput.write("{}\t{}\n".format("Number of gaps", gap_count))
        if gaps_output is not None:
            # NOTE: generate gaps statistics file
            with open(gaps_output, "w") as goutput:
                for key in seq_id_Ngaprange:
                    for rng in seq_id_Ngaprange[key]:
                        goutput.write("{}\t{}\t{}\n".format(key, rng[0], rng[1]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", required=True, help="FASTA file")
    parser.add_argument(
        "-z",
        "--genome_size",
        required=False,
        type=int,
        help="If provided, the NG50 statistic will be computed",
        default=0,
    )
    parser.add_argument(
        "-s",
        "--stats_output",
        required=True,
        help="File to store the general statistics",
    )
    parser.add_argument(
        "-r",
        "--gaps_output",
        required=False,
        help="File to store the gaps statistics",
        default=None,
    )
    args = parser.parse_args()

    run(
        args.fasta,
        args.stats_output,
        args.gaps_output,
        args.genome_size,
    )
