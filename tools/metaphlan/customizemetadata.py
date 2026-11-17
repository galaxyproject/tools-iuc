#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import bz2
import json
import pickle
import re
from pathlib import Path


def load_from_json(json_fp):
    '''
    Read JSON file with marker metadata

    :param json_fp: Path to JSON file
    '''
    with open(json_fp, 'r') as json_f:
        data = json.load(json_f)

    for m in data['markers']:
        data['markers'][m]['ext'] = set(data['markers'][m]['ext'])

    for t in data['taxonomy']:
        if isinstance(data['taxonomy'][t], list):
            data['taxonomy'][t] = tuple(data['taxonomy'][t])
    return data


def dump_to_json(data, json_fp):
    '''
    Dump marker metadata to JSON file

    :param json_fp: Path to JSON file
    '''
    for m in data['markers']:
        data['markers'][m]['ext'] = list(data['markers'][m]['ext'])

    with open(json_fp, 'w') as json_f:
        json.dump(data, json_f)


def transform_pkl_to_json(pkl_fp, json_fp):
    '''
    Read Pickle file and drop it to a JSON file

    :param pkl_fp: Path to input Pickle file
    :param json_fp: Path to output JSON file
    '''
    # load metadata from Pickle file
    with bz2.BZ2File(pkl_fp, 'r') as pkl_f:
        in_metadata = pickle.load(pkl_f)

    out_metadata = {
        'markers': in_metadata['markers'],
        'taxonomy': in_metadata['taxonomy'],
        'merged_taxon': {}
    }

    # transform merged_taxons tuple keys to string
    for k in in_metadata['merged_taxon']:
        n = ' , '.join(k)
        out_metadata[n] = in_metadata['merged_taxon'][k]

    # dump metadata to JSON file
    dump_to_json(out_metadata, json_fp)


def transform_json_to_pkl(json_fp, pkl_fp):
    '''
    Read JSON file and drop it to a Pickle file

    :param json_fp: Path to input JSON file
    :param pkl_fp: Path to output Pickle file
    '''
    # load metadata from JSON file
    in_metadata = load_from_json(json_fp)

    out_metadata = {
        'markers': in_metadata['markers'],
        'taxonomy': in_metadata['taxonomy'],
        'merged_taxon': {}
    }
    # transform merged_taxons keys to tuple
    for k in in_metadata['merged_taxon']:
        n = ' , '.split(k)
        out_metadata[n] = in_metadata['merged_taxon'][k]

    # Ensure that there are 8 taxonomy levels (for compatibility between Metaphlan v3 and v4)
    # v3 DB release encodes the taxids as: ('2|1224|1236|91347|543|547|354276', 4404432)
    # v4 DB release encodes the taxids as: ('2|1224|1236|91347|543|547|354276|', 4404432)
    for k in out_metadata['taxonomy']:
        if out_metadata['taxonomy'][k][0].count('|') == 6:
            out_metadata['taxonomy'][k] = (out_metadata['taxonomy'][k][0] + '|', out_metadata['taxonomy'][k][1])

    # dump metadata to Pickle file
    with bz2.BZ2File(pkl_fp, 'w') as pkl_f:
        pickle.dump(out_metadata, pkl_f)


def add_marker(in_json_fp, out_json_fp, name, m_length, g_length, gca, k_name, k_id, p_name, p_id, c_name, c_id, o_name, o_id, f_name, f_id, g_name, g_id, s_name, s_id, t_name):
    '''
    Add marker to JSON file

    :param in_json_fp: Path to input JSON file
    :param out_json_fp: Path to output JSON file
    :param name: Name of new marker
    :param m_length: Length of new marker
    :param g_length: List with lengths of genomes from which the new marker has been extracted
    :param gca: List with GCA of genomes from which the new marker has been extracted
    :param k_name: List with Name of Kingdom for genomes from which the new marker has been extracted
    :param k_id: List with NCBI id of Kingdom for genomes from which the new marker has been extracted
    :param p_name: List with Name of Phylum for genomes from which the new marker has been extracted
    :param p_id: List with NCBI id of Phylum for genomes from which the new marker has been extracted
    :param c_name: List with Name of Class for genomes from which the new marker has been extracted
    :param c_id: List with NCBI id of Class for genomes from which the new marker has been extracted
    :param o_name: List with Name of Order for genomes from which the new marker has been extracted
    :param o_id: List with NCBI id of Order for genomes from which the new marker has been extracted
    :param f_name: List with Name of Family for genomes from which the new marker has been extracted
    :param f_id: List with NCBI id of Family for genomes from which the new marker has been extracted
    :param g_name: List with Name of Genus for genomes from which the new marker has been extracted
    :param g_id: List with NCBI id of Genus for genomes from which the new marker has been extracted
    :param s_name: List with Name of Species for genomes from which the new marker has been extracted
    :param s_id: List with NCBI id of Species for genomes from which the new marker has been extracted
    :param t_name: List with Name of Strain for genomes from which the new marker has been extracted
    '''
    metadata = load_from_json(in_json_fp)

    # check that all lists have same size
    genome_n = len(g_length)
    if len(gca) != genome_n:
        raise ValueError("Missing/Extra values in GCA list")
    if len(k_name) != genome_n:
        raise ValueError("Missing/Extra values in Kingdom name list")
    if len(k_id) != genome_n:
        raise ValueError("Missing/Extra values in Kingdom ID list")
    if len(p_name) != genome_n:
        raise ValueError("Missing/Extra values in Phylum name list")
    if len(p_id) != genome_n:
        raise ValueError("Missing/Extra values in Phylum ID list")
    if len(c_name) != genome_n:
        raise ValueError("Missing/Extra values in Class name list")
    if len(c_id) != genome_n:
        raise ValueError("Missing/Extra values in Class ID list")
    if len(o_name) != genome_n:
        raise ValueError("Missing/Extra values in Order name list")
    if len(o_id) != genome_n:
        raise ValueError("Missing/Extra values in Order ID list")
    if len(f_name) != genome_n:
        raise ValueError("Missing/Extra values in Family name list")
    if len(f_id) != genome_n:
        raise ValueError("Missing/Extra values in Family ID list")
    if len(g_name) != genome_n:
        raise ValueError("Missing/Extra values in Genus name list")
    if len(g_id) != genome_n:
        raise ValueError("Missing/Extra values in Genus ID list")
    if len(s_name) != genome_n:
        raise ValueError("Missing/Extra values in Species name list")
    if len(s_id) != genome_n:
        raise ValueError("Missing/Extra values in Species ID list")
    if len(t_name) != genome_n:
        raise ValueError("Missing/Extra values in Strain name list")

    # create dictionary to aggregate genome taxonomies and identify marker taxonomy
    taxonomy = {
        'k': set(),
        'p': set(),
        'c': set(),
        'o': set(),
        'f': set(),
        'g': set(),
        's': set(),
        't': set(),
    }

    # parse genomes
    for i in range(genome_n):
        # add taxonomy of new genome
        g_taxo_names = "k__%s|p__%s|c__%s|o__%s|f__%s|g__%s|s__%s|t__%s" % (
            k_name[i],
            p_name[i],
            c_name[i],
            o_name[i],
            f_name[i],
            g_name[i],
            s_name[i],
            t_name[i]
        )
        g_taxo_ids = "%s|%s|%s|%s|%s|%s|%s" % (
            k_id[i],
            p_id[i],
            c_id[i],
            o_id[i],
            f_id[i],
            g_id[i],
            s_id[i]
        )
        metadata['taxonomy'][g_taxo_names] = (g_taxo_ids, g_length[i])
        # aggregate taxon levels using sets
        taxonomy['k'].add(k_name[i])
        taxonomy['p'].add(p_name[i])
        taxonomy['c'].add(c_name[i])
        taxonomy['o'].add(o_name[i])
        taxonomy['f'].add(f_name[i])
        taxonomy['g'].add(g_name[i])
        taxonomy['s'].add(s_name[i])
        taxonomy['t'].add(t_name[i])

    # extract clade and taxon of marker
    clade = ''  # last level before taxomy of genomes diverge
    taxon = ''  # combination of levels before divergence
    for level in ['k', 'p', 'c', 'o', 'f', 'g', 's', 't']:
        taxo = list(taxonomy[level])
        if len(taxo) == 1:
            clade = taxo[0]
            taxon = "%s|%s__%s" % (taxon, level, taxo)

    # add information about the new marker
    metadata['markers'][name] = {
        'clade': clade,
        'ext': set(gca),
        'len': m_length,
        'taxon': taxon
    }

    dump_to_json(metadata, out_json_fp)


def format_markers(marker_l):
    '''
    Format markers

    :param marker_l: list of markers
    '''
    markers = []
    for m in marker_l:
        m = m.rstrip()
        if ' ' in m:
            markers.append(m.split(' ')[0])
        else:
            markers.append(m)
    return markers


def get_markers(marker_fp):
    '''
    Get markers from a file

    :param marker_fp: Path to file with markers (1 per line)
    '''
    # load markers
    with open(marker_fp, 'r') as marker_f:
        markers = marker_f.readlines()

    # format markers
    markers = format_markers(markers)

    return markers


def check_not_found_markers(found_markers, original_markers):
    '''
    Check list of markers

    :param found_markers: list of found markers
    :param original_markers: list of original markers
    '''
    if len(found_markers) != len(original_markers):
        print('markers not found:')
        for m in original_markers:
            if m not in found_markers:
                print('- "%s"' % m)


def prune_taxonomy(in_taxonomy, taxon_s, gca_s):
    '''
    Prune taxonomy to keep only listed taxonomy

    :param in_taxonomy: dictionary with list of taxonomy
    :param taxon_s: set of taxons to keep
    :param gca_s: set of GCA ids to keep
    '''
    out_taxonomy = {}
    kept_taxonomy = set()
    kept_taxons = set()
    kept_gca = set()
    for t, v in in_taxonomy.items():
        # check if t match element in list of taxon_s
        kept_taxon = False
        for t_k in taxon_s:
            if t_k in t:
                kept_taxon = True
                out_taxonomy[t] = v
                kept_taxonomy.add(t)
                kept_taxons.add(t_k)
                break
        # check if GCA in the taxon id
        s = re.search(r'GCA_\d+$', t)
        if s:
            gca = s[0]
            # check if GCA in taxon id is in the list GCA to keep
            if gca in gca_s:
                kept_gca.add(gca)
                if not kept_taxon:
                    out_taxonomy[t] = v
                    kept_taxonomy.add(t)

    print('%s kept taxonomy' % len(kept_taxonomy))
    print('%s / %s taxons not found' % (len(taxon_s) - len(kept_taxons), len(taxon_s)))
    print('%s / %s GCA taxons not found' % (len(gca_s) - len(kept_gca), len(gca_s)))
    return out_taxonomy


def remove_markers(in_json_fp, marker_fp, out_json_fp, kept_marker_fp):
    '''
    Remove markers from JSON file

    :param in_json_fp: Path to input JSON file
    :param marker_fp: Path to file with markers to remove (1 per line)
    :param out_json_fp: Path to output JSON file
    :param kept_marker_fp: Path to file with kept markers
    '''
    in_metadata = load_from_json(in_json_fp)

    # load markers
    markers_to_remove = set(get_markers(marker_fp))
    print('%s markers to remove' % len(markers_to_remove))

    # keep merged_taxon
    out_metadata = {
        'markers': {},
        'taxonomy': {},
        'merged_taxon': in_metadata['merged_taxon']
    }

    # parse markers to keep
    removed_markers = []
    kept_markers = []
    taxons_to_keep = set()
    gca_to_keep = set()
    for m, v in in_metadata['markers'].items():
        if m not in markers_to_remove:
            out_metadata['markers'][m] = v
            kept_markers.append(m)
            taxons_to_keep.add(v['taxon'])
            gca_to_keep.update(v['ext'])
        else:
            removed_markers.append(m)
    print('%s removed markers' % len(removed_markers))

    # check markers that are not found
    check_not_found_markers(removed_markers, markers_to_remove)

    # keep only taxonomy in taxons_to_keep or with GCA in gca_to_keep
    out_metadata['taxonomy'] = prune_taxonomy(in_metadata['taxonomy'], taxons_to_keep, gca_to_keep)

    # save to JSON
    dump_to_json(out_metadata, out_json_fp)

    # write list of kept markers
    with open(kept_marker_fp, 'w') as kept_marker_f:
        for m in kept_markers:
            kept_marker_f.write("%s\n" % m)


def keep_markers(in_json_fp, marker_fp, out_json_fp):
    '''
    Keep markers from JSON file, others will be removed

    :param in_json_fp: Path to input JSON file
    :param marker_fp: Path to file with markers to keep (1 per line)
    :param out_json_fp: Path to output JSON file
    '''
    in_metadata = load_from_json(in_json_fp)

    # load markers
    markers_to_keep = set(get_markers(marker_fp))
    print('%s markers to keep' % len(markers_to_keep))

    # keep merged_taxon
    out_metadata = {
        'markers': {},
        'taxonomy': {},
        'merged_taxon': in_metadata['merged_taxon']
    }

    # parse markers to keep
    kept_markers = []
    taxons_to_keep = set()
    gca_to_keep = set()
    for m, v in in_metadata['markers'].items():
        if m in markers_to_keep:
            out_metadata['markers'][m] = v
            kept_markers.append(m)
            taxons_to_keep.add(v['taxon'])
            gca_to_keep.update(v['ext'])
    print('%s kept markers' % len(kept_markers))

    # check markers that are not found
    check_not_found_markers(kept_markers, markers_to_keep)

    # keep only taxonomy in taxons_to_keep or with GCA in gca_to_keep
    out_metadata['taxonomy'] = prune_taxonomy(in_metadata['taxonomy'], taxons_to_keep, gca_to_keep)

    # save to JSON
    dump_to_json(out_metadata, out_json_fp)


if __name__ == '__main__':
    # Read command line
    parser = argparse.ArgumentParser(description='Customize MetaPhlan database')
    subparsers = parser.add_subparsers(dest='function')
    # transform_pkl_to_json subcommand
    pkl_to_json_parser = subparsers.add_parser('transform_pkl_to_json', help='Transform Pickle to JSON to get marker metadata')
    pkl_to_json_parser.add_argument('--pkl', help="Path to input Pickle file")
    pkl_to_json_parser.add_argument('--json', help="Path to output JSON file")
    # transform_json_to_pkl subcommand
    json_to_pkl_parser = subparsers.add_parser('transform_json_to_pkl', help='Transform JSON to Pickle to push marker metadata')
    json_to_pkl_parser.add_argument('--json', help="Path to input JSON file")
    json_to_pkl_parser.add_argument('--pkl', help="Path to output Pickle file")
    # add_marker subcommand
    add_marker_parser = subparsers.add_parser('add_marker', help='Add new marker to JSON file')
    add_marker_parser.add_argument('--in_json', help="Path to input JSON file")
    add_marker_parser.add_argument('--out_json', help="Path to output JSON file")
    add_marker_parser.add_argument('--name', help="Name of new marker")
    add_marker_parser.add_argument('--m_length', help="Length of new marker")
    add_marker_parser.add_argument('--g_length', help="Length of genome from which the new marker has been extracted", action="append")
    add_marker_parser.add_argument('--gca', help="GCA of genome from which the new marker has been extracted", action="append")
    add_marker_parser.add_argument('--k_name', help="Name of Kingdom for genome from which the new marker has been extracted", action="append")
    add_marker_parser.add_argument('--k_id', help="NCBI id of Kingdom for genome from which the new marker has been extracted", action="append")
    add_marker_parser.add_argument('--p_name', help="Name of Phylum for genome from which the new marker has been extracted", action="append")
    add_marker_parser.add_argument('--p_id', help="NCBI id of Phylum for genome from which the new marker has been extracted", action="append")
    add_marker_parser.add_argument('--c_name', help="Name of Class for genome from which the new marker has been extracted", action="append")
    add_marker_parser.add_argument('--c_id', help="NCBI id of Class for genome from which the new marker has been extracted", action="append")
    add_marker_parser.add_argument('--o_name', help="Name of Order for genome from which the new marker has been extracted", action="append")
    add_marker_parser.add_argument('--o_id', help="NCBI id of Order for genome from which the new marker has been extracted", action="append")
    add_marker_parser.add_argument('--f_name', help="Name of Family for genome from which the new marker has been extracted", action="append")
    add_marker_parser.add_argument('--f_id', help="NCBI id of Family for genome from which the new marker has been extracted", action="append")
    add_marker_parser.add_argument('--g_name', help="Name of Genus for genome from which the new marker has been extracted", action="append")
    add_marker_parser.add_argument('--g_id', help="NCBI id of Genus for genome from which the new marker has been extracted", action="append")
    add_marker_parser.add_argument('--s_name', help="Name of Species for genome from which the new marker has been extracted", action="append")
    add_marker_parser.add_argument('--s_id', help="NCBI id of Species for genome from which the new marker has been extracted", action="append")
    add_marker_parser.add_argument('--t_name', help="Name of Strain for genome from which the new marker has been extracted", action="append")
    # remove_markers subcommand
    remove_markers_parser = subparsers.add_parser('remove_markers', help='Remove markers from JSON file')
    remove_markers_parser.add_argument('--in_json', help="Path to input JSON file")
    remove_markers_parser.add_argument('--markers', help="Path to file with markers to remove (1 per line)")
    remove_markers_parser.add_argument('--out_json', help="Path to output JSON file")
    remove_markers_parser.add_argument('--kept_markers', help="Path to file with kept markers")
    # keep_markers subcommand
    keep_markers_parser = subparsers.add_parser('keep_markers', help='Keep markers from JSON file, others will be removed')
    keep_markers_parser.add_argument('--in_json', help="Path to input JSON file")
    keep_markers_parser.add_argument('--markers', help="Path to file with markers to keep (1 per line)")
    keep_markers_parser.add_argument('--out_json', help="Path to output JSON file")

    args = parser.parse_args()

    if args.function == 'transform_pkl_to_json':
        transform_pkl_to_json(Path(args.pkl), Path(args.json))
    elif args.function == 'transform_json_to_pkl':
        transform_json_to_pkl(Path(args.json), Path(args.pkl))
    elif args.function == 'add_marker':
        add_marker(
            args.in_json,
            args.out_json,
            args.name,
            args.m_length,
            args.g_length,
            args.gca,
            args.k_name,
            args.k_id,
            args.p_name,
            args.p_id,
            args.c_name,
            args.c_id,
            args.o_name,
            args.o_id,
            args.f_name,
            args.f_id,
            args.g_name,
            args.g_id,
            args.s_name,
            args.s_id,
            args.t_name)
    elif args.function == 'remove_markers':
        remove_markers(args.in_json, args.markers, args.out_json, args.kept_markers)
    elif args.function == 'keep_markers':
        keep_markers(args.in_json, args.markers, args.out_json)
