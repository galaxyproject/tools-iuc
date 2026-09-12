#!/usr/bin/env python3
"""PyEuk distance + SWEEP + REPORT driver for the Galaxy tool ``haplotype_pyeuk``.

Reads a haplotype sheet, computes the wIBS distance matrix with
``PyEukDistanceEngine.compute_revised_wibs_matrix`` and runs the OFFICIAL PyEuk
0.7.0 clustering output -- ``CyclosporaClusterFinder.cluster_sweep`` -- which
reports the count RANGE the cohort supports, its confidence, the stable cores and
the confidence tree, and writes a representative flat partition for downstream
tools. It then renders the ``pyeuk report`` dashboard (galaxy theme, no external
assets) from that sweep.

Everything numeric is the pyeuk 0.7.0 package; this file only marshals files in
and out, prints diagnostics, and keeps the sheet-validation + pairwise-completeness
checks the Galaxy arm has always carried. Single-k (the legacy one-partition
verdict) is available in the package CLI (``pyeuk cluster --single-k``) but is NOT
the tool's output: the sweep is the official answer.
"""

import argparse
import glob
import json
import os
import sys

import pandas as pd

from pyeuk.clustering import CyclosporaClusterFinder
from pyeuk.distance_engine import PyEukDistanceEngine, parse_locus_name
from pyeuk.report import render

EPSILON = 0.3072


def die(msg):
    sys.stderr.write("haplotype_pyeuk: %s\n" % msg)
    sys.exit(1)


def read_sheet(path):
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception as exc:  # noqa: BLE001
        die("could not parse the haplotype sheet as TSV: %s" % exc)
    df.columns = [str(c).strip() for c in df.columns]
    if df.shape[1] < 2:
        die("the haplotype sheet has %d column(s); expected 'Seq_ID' plus at "
            "least one haplotype column" % df.shape[1])
    if "Seq_ID" not in df.columns:
        die("the haplotype sheet has no 'Seq_ID' column (first column is %r)."
            % df.columns[0])
    df["Seq_ID"] = df["Seq_ID"].astype(str).str.strip()
    df = df[(df["Seq_ID"] != "") & (df["Seq_ID"].str.lower() != "nan")].copy()
    if df.empty:
        die("the haplotype sheet has no specimen rows")
    dup = df["Seq_ID"][df["Seq_ID"].duplicated()].unique().tolist()
    if dup:
        die("duplicate Seq_ID values in the haplotype sheet: %s" % ", ".join(dup[:10]))
    marker_cols = [c for c in df.columns if c != "Seq_ID"]
    values = set()
    for c in marker_cols:
        values.update(v for v in df[c].dropna().astype(str).unique())
    odd = sorted(v for v in values if v.strip() not in ("X", ""))
    if odd:
        sys.stderr.write("haplotype_pyeuk: warning: %d marker value(s) other than "
                         "'X' or empty present and treated as ABSENT: %s\n"
                         % (len(odd), ", ".join(repr(v) for v in odd[:10])))
    return df, marker_cols


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sheet", required=True)
    ap.add_argument("--matrix-out", required=True)
    ap.add_argument("--clusters-out", required=True, help="representative flat partition (TSV)")
    ap.add_argument("--sweep-out", required=True, help="SWEEP.json (official output)")
    ap.add_argument("--report-out", required=True, help="graphical HTML report")
    ap.add_argument("--min-completeness", type=float, default=0.10)
    ap.add_argument("--weight-mode", default="heterozygosity",
                    choices=["heterozygosity", "king", "none"])
    ap.add_argument("--project-psd", choices=["true", "false"], default="false")
    ap.add_argument("--k-min", type=int, default=2)
    ap.add_argument("--k-max", type=int, default=50)
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--linkage-method", default="ward",
                    choices=["ward", "single", "average", "complete"])
    ap.add_argument("--report-flavor", default="dashboard",
                    choices=["dashboard", "clinical", "narrative"])
    ap.add_argument("--report-theme", default="galaxy", choices=["studio", "galaxy"])
    ap.add_argument("--report-excluded", choices=["true", "false"], default="true")
    ap.add_argument("--float-format", default="%.10g")
    # Accepted for workflow compatibility; single-k tree-cut concepts, not used by the sweep.
    ap.add_argument("--cut", choices=["count", "distance"], default="count", help=argparse.SUPPRESS)
    ap.add_argument("--linkage-threshold", type=float, default=None, help=argparse.SUPPRESS)
    ap.add_argument("--gold", default=None, help=argparse.SUPPRESS)
    args = ap.parse_args()

    if args.k_min < 2:
        die("--k-min must be >= 2")
    if args.k_max < args.k_min:
        die("--k-max (%d) is below --k-min (%d)" % (args.k_max, args.k_min))
    if args.cut == "distance" or args.linkage_threshold is not None or args.gold:
        sys.stderr.write("haplotype_pyeuk: note: --cut/--linkage-threshold/--gold apply to the "
                         "legacy single-k path only; the sweep is unsupervised and reports a "
                         "count range, so they are ignored here.\n")

    df, marker_cols = read_sheet(args.sheet)
    all_ids = df["Seq_ID"].tolist()
    print("[haplotype_pyeuk] specimens in sheet: %d" % len(all_ids))
    print("[haplotype_pyeuk] haplotype columns : %d" % len(marker_cols))
    sys.stdout.flush()

    try:
        engine = PyEukDistanceEngine(epsilon=EPSILON,
                                     min_completeness=args.min_completeness,
                                     weight_mode=args.weight_mode,
                                     project_psd=(args.project_psd == "true"))
    except TypeError:
        die("this PyEuk build does not accept weight_mode=%r; expected PyEuk >= 0.4.0."
            % args.weight_mode)
    matrix = engine.compute_revised_wibs_matrix(df)
    if matrix.shape[0] == 0:
        die("no specimen passed the completeness filter (--min-completeness %g)."
            % args.min_completeness)

    # pairwise completeness diagnostic (unchanged) --------------------------------------
    import itertools as _it
    _called = {}
    for _sid in df["Seq_ID"]:
        _row = df[df["Seq_ID"] == _sid].iloc[0]
        _called[_sid] = {parse_locus_name(_c) for _c in marker_cols
                         if str(_row[_c]).strip().upper() == "X"}
    _ids = list(df["Seq_ID"]); _tot = _shared = 0; _sh = []
    for _a, _b in _it.combinations(_ids, 2):
        _tot += 1; _n = len(_called[_a] & _called[_b]); _sh.append(_n); _shared += 1 if _n else 0
    if _tot:
        _pct = 100.0 * _shared / _tot; _sh.sort()
        print("[haplotype_pyeuk] pairwise completeness : %d/%d pairs (%.1f%%) share >=1 called "
              "locus; median shared loci %d" % (_shared, _tot, _pct, _sh[len(_sh) // 2]))
        if _pct < 90.0:
            print("[haplotype_pyeuk] WARNING: %.1f%% of pairs share NO called locus; their "
                  "distance is the engine ceiling, not a measurement." % (100.0 - _pct))
    sys.stdout.flush()

    matrix.to_csv(args.matrix_out, sep="\t", index_label="Seq_ID", float_format=args.float_format)

    # ---- OFFICIAL OUTPUT: the sweep --------------------------------------------------
    finder = CyclosporaClusterFinder()
    workdir = os.path.join(os.getcwd(), "sweep_work")
    sweep = finder.cluster_sweep(matrix, k_min=args.k_min, k_max=args.k_max,
                                 n_boot=args.n_boot, linkage_method=args.linkage_method,
                                 output_dir=workdir)
    with open(args.sweep_out, "w") as fh:
        json.dump(sweep, fh, indent=1)

    # representative flat partition (from the sweep), + excluded specimens as -1
    cfiles = sorted(glob.glob(os.path.join(workdir, "*_RESULTING_CLUSTERS_*.txt")))
    if cfiles:
        rep = pd.read_csv(cfiles[-1], sep="\t")
    else:
        rep = pd.DataFrame({"Seq_ID": list(matrix.index), "Assigned_cluster": 1})
    rep["Seq_ID"] = rep["Seq_ID"].astype(str)
    if args.report_excluded == "true":
        missing = [s for s in all_ids if s not in set(rep["Seq_ID"])]
        if missing:
            rep = pd.concat([rep, pd.DataFrame({"Seq_ID": missing, "Assigned_cluster": -1})],
                            ignore_index=True)
    rep.to_csv(args.clusters_out, sep="\t", index=False)

    # graphical report (galaxy theme = system fonts, no external assets)
    html = render(sweep, dist_df=matrix, flavor=args.report_flavor, theme=args.report_theme)
    with open(args.report_out, "w") as fh:
        fh.write(html)

    cr = sweep.get("count_range"); pe = sweep.get("point_estimate")
    print("[haplotype_pyeuk] SWEEP count range   : %s" % (cr,))
    print("[haplotype_pyeuk] confident / point   : %s / %s" % (sweep.get("confident"), pe))
    print("[haplotype_pyeuk] representative k     : %s" % sweep.get("representative_k"))
    print("[haplotype_pyeuk] headline            : %s" % sweep.get("headline"))
    print("[haplotype_pyeuk] report              : %s (%s/%s)"
          % (args.report_out, args.report_flavor, args.report_theme))


if __name__ == "__main__":
    main()
