#!/usr/bin/env python
import argparse
import collections
import csv
import sys

from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


def c(row, v, default=None):
    if v is None:
        return default

    try:
        _ = int(v)
        return row[int(v) - 1]
    except ValueError:
        return v


def tbl2gff3(
    table,
    rid,
    begin,
    end,
    source=None,
    type=None,
    score=None,
    frame=None,
    a=None,
    strand_column=None,
    strand_value=None,
    strand_infer=False,
):

    records = collections.OrderedDict()

    for row in csv.reader(table, delimiter="\t"):
        # print(', '.join(row))

        # if we haven't seen this record before, populate it.
        recid = c(row, rid)
        if recid not in records:
            records[recid] = SeqRecord(Seq("ACTG"), id=recid)

        r = records[recid]
        q = {}
        if c(row, score) is not None:
            q["score"] = float(c(row, score))

        q["source"] = c(row, source, "tbl2gff3")

        begin_i = int(c(row, begin))
        end_i = int(c(row, end))

        begin_f = min(begin_i, end_i)
        end_f = max(begin_i, end_i)

        _str = None
        if strand_column is not None:
            _str = int(c(row, strand_column))
        elif strand_value is not None:
            _str = int(strand_value)
        if strand_infer:
            if begin_i > begin_f:
                _str = -1
            else:
                _str = 1

        if a is not None:
            for x in a:
                k, v = x.split(":", 1)
                _v = c(row, v)
                if k in q:
                    q[k].append(_v)
                else:
                    q[k] = [_v]

        f = SeqFeature(
            FeatureLocation(begin_f, end_f),
            type=c(row, type),
            strand=_str,
            qualifiers=q,
        )
        r.features.append(f)

    return records


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert tables to gff3", epilog="")
    parser.add_argument("table", type=argparse.FileType("r"), help="Tabular Input")
    parser.add_argument("rid", help="id column")
    parser.add_argument("begin", help="begin column")
    parser.add_argument("end", help="end column")
    parser.add_argument("--type", help="feature type column")
    parser.add_argument("--score", help="score column")
    parser.add_argument("--source", help="source column")
    parser.add_argument("--strand_infer", action='store_true', help="infer strand")
    parser.add_argument("--strand_column", help="strand column")
    parser.add_argument("--strand_value", help="strand value")
    # parser.add_argument('--frame', help='frame column')
    parser.add_argument("-a", action="append", help="attribute column (-a k:v)")
    args = parser.parse_args()

    for rid, rec in tbl2gff3(**vars(args)).items():
        GFF.write([rec], sys.stdout)
