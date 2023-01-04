#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Grigorii Sukhorukov, Macha Nikolski

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import pathlib
import math
import random
from sklearn.utils import shuffle
import h5py


def reverse_complement(fragment):
    """
    provides reverse complement to sequences
    Input:
    sequences - list with SeqRecord sequences in fasta format
    Output:
    complementary_sequences -
    list with SeqRecord complementary sequences in fasta format
    """
    # complementary_sequences = []
    # for sequence in sequences:
    #     complementary_sequence = SeqRecord(
    #         seq=Seq(sequence.seq).reverse_complement(),
    #         id=sequence.id + "_reverse_complement",
    #     )
    #     complementary_sequences.append(complementary_sequence)
    fragment = fragment[::-1].translate(str.maketrans('ACGT', 'TGCA'))
    return fragment


def introduce_mutations(seqs, mut_rate, rs=None):
    """
    Function that mutates sequences in the entering fasta file
    A proportion of nucleotides are changed to other nucleotide
    Not yet taking account of mutation for gaps
    mut_rate - proportion from 0.0 to 1.0, float
    """
    random.seed(a=rs)
    assert 0.0 <= mut_rate <= 1.0
    mutated_seqs = []
    for seq in seqs:
        mut_seq = list(str(seq.seq))
        l_ = len(mut_seq)
        mutated_sites_i = random.sample(range(l_), int(mut_rate * l_))
        for mut_site_i in mutated_sites_i:
            mut_site = mut_seq[mut_site_i]
            mutations = ["A", "C", "T", "G"]
            if mut_site in mutations:
                mutations.remove(mut_site)
                mut_seq[mut_site_i] = random.sample(mutations, 1)[0]
        mutated_seq = SeqRecord(
            seq=Seq("".join(mut_seq)),
            id=seq.id + f"mut_{mut_rate}",
            name="",
            description="",
        )
        mutated_seqs.append(mutated_seq)
    return mutated_seqs


def separate_by_length(length_, seq_list, fold=None,):
    # TODO: add docs
    included = []
    to_process = []
    excluded = 0
    for seq_ in seq_list:
        l_ = len(seq_.seq)
        if l_ >= length_:
            if fold is None:
                included.append(seq_)
            elif l_ < length_ * fold:
                included.append(seq_)
            else:
                to_process.append(seq_)
        else:
            excluded += 1
    print(f"A total of {excluded} sequences was excluded due to being smaller than {length_}")
    return included, to_process


def chunks(lst, n):
    """Yield successive n-sized chunks from lst.
    https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks"""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def correct(frag):
    """
    leaves only unambiguous DNA code (ACTG-)
    Input:
    frag - string of nucleotides
    Output:
    pr_frag - corrected string of nucleotides
    """
    pr_frag = frag.upper()
    pr_frag_s = set(pr_frag)
    if pr_frag_s != {"A", "C", "G", "T", "-"}:
        for letter in pr_frag_s - {"A", "C", "G", "T", "-"}:
            pr_frag = pr_frag.replace(letter, "-")
    return pr_frag


def fragmenting(sequences, sl_wind_size, max_gap=0.05, sl_wind_step=None):
    """
    slices sequences in fragments by sliding window
    based on its size and step.
    last fragment is padded by '-'
    fragments have ambiguous bases replaced by '-'
    fragments with many '-' are discarded
    Input:
    sequences - list with SeqRecord sequences in fasta format
    max_gap - max allowed proportion of '-'
    sl_wind_size - sliding window step
    sl_wind_step - sliding window step, by default equals
    sliding window size (None is replaced by it)
    Output:
    fragments - list with sequence fragments
    """
    if sl_wind_step is None:
        sl_wind_step = sl_wind_size
    fragments = []
    fragments_rc = []
    out_sequences = []
    for sequence in sequences:
        seq = str(sequence.seq)
        n_fragments = 1 + max(0, math.ceil((len(seq) - sl_wind_size)/sl_wind_step))
        for n in range(n_fragments):
            if n + 1 != n_fragments:
                frag = seq[n * sl_wind_step: n * sl_wind_step + sl_wind_size]
            elif n_fragments == 1:
                # padding the shorter fragment to sl_wind_size
                frag_short = seq[n * sl_wind_step: n * sl_wind_step + sl_wind_size]
                frag = frag_short + (sl_wind_size - len(frag_short)) * "-"
            else:
                frag = seq[(len(seq) - sl_wind_size):]
            # replace ambiguous characters
            frag = correct(frag)
            assert len(frag) == sl_wind_size, f"{len(frag)} vs {sl_wind_size}"
            # skipping sequences with many gaps
            if frag.count("-") / sl_wind_size <= max_gap:
                fragments.append(frag)
                # generating reverse complement
                fragments_rc.append(reverse_complement(frag))
                fr_seq = SeqRecord(
                    seq=Seq(frag),
                    id=f"{sequence.id}_{n*sl_wind_step}_{sl_wind_size}",
                    name="",
                    description="",
                )
                out_sequences.append(fr_seq)
    return fragments, fragments_rc, out_sequences


def label_fasta_fragments(sequences, label):
    """
    Provides labels to generated fragments stored in fasta
    Input:
    sequences - list with SeqRecord sequences
    label - type of label (bacteria, virus, plant)
    Output:
    labeled_fragments - list with labeled SeqRecord sequences
    """
    # assert label in ["virus", "plant", "bacteria"]
    labeled_fragments = []
    for sequence in sequences:
        sequence.id = sequence.id + f"_{label}"
        labeled_fragments.append(sequence)
    return labeled_fragments


def one_hot_encode(fragments):
    """
    produces one-hot matrices from fragments and labels
    '-' is given all zeros
    Input:
    fragments - list with sequence fragments
    label - type of label (int <= depth)
    label_depth - number of possible labels
    Output:
    encoded_fragments - list with one-hot encoded fragments
    labels - list with one-hot encoded labels
    """
    import tensorflow as tf
    encoded_fragments = []
    map_dict = {"A": 0, "C": 1, "G": 2, "T": 3, "-": -1}
    for frag in fragments:
        frag_array = np.array(list(frag))
        integer_encoded = np.int8(np.vectorize(map_dict.get)(frag_array))
        one_hot_encoded = tf.one_hot(integer_encoded, depth=4, dtype=tf.int8).numpy()
        encoded_fragments.append(one_hot_encoded)
    encoded_fragments = np.stack(encoded_fragments)
    return encoded_fragments


def prepare_labels(fragments, label, label_depth):
    """
    produces one-hot labels
    '-' is given all zeros
    Input:
    fragments - list with sequence fragments
    label - type of label (int <= depth)
    label_depth - number of possible labels
    Output:
    labels - list with one-hot encoded labels
    """
    import tensorflow as tf
    n_fragments = len(fragments)
    labels = np.int8(np.full(n_fragments, label))
    labels = tf.one_hot(labels, depth=label_depth).numpy()
    return labels


# TODO: write docs for functions
def calculate_total_length(seq_path):
    """
    Calculate total length of the sequences in the fasta file.
    Needed for weighted sampling
    Input:
    seq_path - path to the file with sequences
    Output:
    seq_length - total length of all sequences in the file
    """
    seqs = list(SeqIO.parse(seq_path, "fasta"))
    seq_length = 0
    for seq in seqs:
        seq_length += len(seq.seq)
    return seq_length


def prepare_seq_lists(in_paths, n_fragments, weights=None,):
    """
    selects files with sequences based on extension
    and calculates number of fragments to be sampled
    Input:
    in_paths - list of paths to folder with sequence files. Can be a string also a string
    n_fragments - number of fragments to be sampled
    weights - upsampling of fragments. fractions should sum to one
    Output:
    seqs_list - list with path to files with sequences
    n_fragments_list - number of fragments to be sampled
    lists are zipped to work with ray iterators
    """
    # case when we recieve a single sequence file
    if type(in_paths) is str and in_paths.endswith(('.fna', '.fasta')):
        return [[in_paths, n_fragments]]
    else:
        # transform string to list
        if type(in_paths) is str or type(in_paths) is pathlib.PosixPath:
            in_paths = [in_paths]

        if weights:
            assert len(weights) == len(in_paths)
            assert 1.01 > round(sum(weights), 2) > 0.99
        else:
            l_ = len(in_paths)
            weights = [1/l_] * l_
        n_fragments_list_all = []
        seqs_list_all = []
        for in_paths, w_ in zip(in_paths, weights):
            seqs_list = []
            seq_length_list = []
            total_length = 0
            for file in os.listdir(in_paths):
                if file.endswith("fna") or file.endswith("fasta"):
                    seq_path = (os.path.join(in_paths, file))
                    seqs_length = calculate_total_length(seq_path)
                    seqs_list.append(seq_path)
                    seq_length_list.append(seqs_length)
                    total_length += seqs_length
            # + 1 may lead to a slightly bigger number than desired
            n_fragments_list = [((seq_length / total_length) * n_fragments * w_ + 1) for seq_length in seq_length_list]
            n_fragments_list_all.extend(n_fragments_list)
            seqs_list_all.extend(seqs_list)
        print("list calculation done")
        return list(zip(seqs_list_all, n_fragments_list_all))


def sample_fragments(seq_container, length, random_seed=1, limit=None, max_gap=0.05, sl_wind_step=None):
    """
    Randomly samples fragments from sequences in the list.
    Is a bit cumbersome written to work with ray.
    Input:
    seq_container - list with each entry containing path to sequence,
    and n samples from this sequence.
    length - desired length of sampled fragments
    Output:
    fragments - list with sequence fragments
    """
    random.seed(a=random_seed)
    total_fragments = []
    total_fragments_rc = []
    total_seqs = []
    for entry in seq_container:
        seq = list(SeqIO.parse(entry[0], "fasta"))
        n_fragments = entry[1]
        seqs = []
        fragments = []
        fragments_rc = []
        counter_1 = 0
        counter_2 = 0
        while counter_1 < n_fragments:
            # select chromosomes if there are any
            fragment_full = random.choice(seq)
            r_end = len(fragment_full.seq) - length
            try:
                r_start = random.randrange(r_end)
                fragment = SeqRecord(
                    seq=fragment_full.seq[r_start:(r_start + length)],
                    id=f"{fragment_full.id}_{length}_{r_start}",
                    name="",
                    description="",
                )
                temp_, temp_rc, _ = fragmenting([fragment], length, max_gap, sl_wind_step=sl_wind_step)
                if temp_ and temp_rc:
                    seqs.append(fragment)
                    fragments.extend(temp_)
                    fragments_rc.extend(temp_rc)
                    counter_1 += 1
            except ValueError:
                # print(f"{fragment_full.id} has length {len(fragment_full.seq)} and is too short to be sampled")
                pass
            counter_2 += 1
            if limit:
                assert counter_2 <= limit * n_fragments, f"While cycle iterated more than {limit}, data is ambiguous." \
                                                         f" Only {len(fragments)} fragments were sampled out of {n_fragments}"
        total_fragments.extend(fragments)
        total_fragments_rc.extend(fragments_rc)
        total_seqs.extend(seqs)
        # print("sequence sampling done")
    return total_fragments, total_fragments_rc, total_seqs


def prepare_ds_fragmenting(in_seq, label, label_int, fragment_length, sl_wind_step, max_gap=0.05, n_cpus=1):
    if sl_wind_step is None:
        sl_wind_step = int(fragment_length / 2)
    # generating viral fragments and labels
    seqs = list(SeqIO.parse(in_seq, "fasta"))
    frags, frags_rc, seqs_ = fragmenting(seqs, fragment_length, max_gap=max_gap, sl_wind_step=sl_wind_step)
    encoded = one_hot_encode(frags)
    encoded_rc = one_hot_encode(frags_rc)
    labs = prepare_labels(frags, label=label_int, label_depth=2)
    seqs_ = label_fasta_fragments(seqs_, label=label)
    # subsetting to unique fragments
    u_encoded, indices = np.unique(encoded, axis=0, return_index=True)
    u_encoded_rc = encoded_rc[indices]
    u_labs = labs[indices]
    u_seqs = [seqs_[i] for i in indices]
    assert (np.shape(u_encoded)[0] == np.shape(u_encoded_rc)[0])
    print(f"Encoding {label} sequences finished")
    # print(f"{np.shape(u_encoded)[0]} forward fragments generated")
    n_frags = np.shape(u_encoded)[0]
    return u_encoded, u_encoded_rc, u_labs, u_seqs, n_frags


def prepare_ds_sampling(in_seqs, fragment_length, n_frags, label, label_int, random_seed,  n_cpus=1, limit=100):
    # generating plant fragments and labels
    seqs_list = prepare_seq_lists(in_seqs, n_frags)
    frags, frags_rc, seqs_ = sample_fragments(seqs_list, fragment_length, random_seed, limit=limit, max_gap=0.05)
    frags, frags_rc, seqs_ = shuffle(frags, frags_rc, seqs_, random_state=random_seed, n_samples=int(n_frags))
    encoded = one_hot_encode(frags)
    encoded_rc = one_hot_encode(frags_rc)
    labs = prepare_labels(frags, label=label_int, label_depth=2)
    seqs_ = label_fasta_fragments(seqs_, label=label)
    assert (np.shape(encoded)[0] == np.shape(encoded_rc)[0])
    print(f"Encoding {label} sequences finished")
    # print(f"{np.shape(encoded)[0]} forward fragments generated")
    return encoded, encoded_rc, labs, seqs_, n_frags


def storing_encoded(encoded, encoded_rc, labs, out_path, ):
    f = h5py.File(out_path, "w")
    f.create_dataset("fragments", data=encoded)
    f.create_dataset("fragments_rc", data=encoded_rc)
    f.create_dataset("labels", data=labs)
    f.close()
