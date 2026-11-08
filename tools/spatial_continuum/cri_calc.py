#!/usr/bin/env python3
"""Biophysical Composite Resistance Index (CRI) and Reaction-Diffusion Solver.

Translational Modeling Tool for Perivascular Drug Penetration and Cohort Risk Stratification.
Research Use Only (RUO). This tool performs mathematical biophysical simulations and risk
scoring; it does not constitute clinical prescribing advice.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from typing import Dict, Iterable, Optional, Sequence, TextIO, Tuple


def bessel_k0_approx(x: float) -> float:
    """Accurate polynomial approximation of modified Bessel function of second kind K0(x).

    Uses Abramowitz and Stegun formulas 9.8.5 and 9.8.7 (error < 1e-7).
    """
    if x <= 0:
        return float("inf")
    if x <= 2.0:
        y = (x / 2.0) ** 2
        # I0 polynomial
        i0 = (
            1.0
            + 3.5156229 * y
            + 3.0899424 * (y**2)
            + 1.2067492 * (y**3)
            + 0.2659732 * (y**4)
            + 0.0360768 * (y**5)
            + 0.0045813 * (y**6)
        )
        k0 = (
            -math.log(x / 2.0) * i0
            - 0.57721566
            + 0.42278420 * y
            + 0.23069756 * (y**2)
            + 0.03488590 * (y**3)
            + 0.00262698 * (y**4)
            + 0.00010750 * (y**5)
            + 0.00000740 * (y**6)
        )
        return k0
    else:
        y = 2.0 / x
        k0_scaled = (
            1.25331414
            - 0.07832358 * y
            + 0.02189568 * (y**2)
            - 0.01062446 * (y**3)
            + 0.00587872 * (y**4)
            - 0.00251540 * (y**5)
            + 0.00053208 * (y**6)
        )
        return (math.exp(-x) / math.sqrt(x)) * k0_scaled


def compute_biophysical_sanctuary(
    pattern: str, lambda_pen: float = 34.1, r_vessel: float = 8.0
) -> Tuple[float, float, str]:
    """Calculate perivascular sanctuary depth and concentration ratio from reaction-diffusion.

    Returns:
        (sanctuary_depth_um, concentration_ratio, summary_text)
    """
    pattern_lower = str(pattern or "").strip().lower()

    # Estimate average tumor cord radius Rc based on IASLC histological architecture
    if "lepidic" in pattern_lower:
        r_cord = 10.0  # Alveolar monolayer, adjacent to capillary
    elif "solid" in pattern_lower or "micropapillary" in pattern_lower:
        r_cord = 100.0  # Expansile avascular nest
    elif "acinar" in pattern_lower or "papillary" in pattern_lower:
        r_cord = 35.0  # Intermediate cord radius
    elif "mixed" in pattern_lower:
        r_cord = 60.0
    else:
        r_cord = 45.0

    # Effective therapeutic penetration radius: r_eff = r_vessel + lambda_pen
    r_eff = r_vessel + lambda_pen
    sanctuary_depth = max(0.0, r_cord - r_eff)

    # Relative concentration at cord boundary: C(Rc) / C(rv) = K0(Rc/lambda) / K0(rv/lambda)
    k0_vessel = bessel_k0_approx(r_vessel / lambda_pen)
    k0_cord = bessel_k0_approx(r_cord / lambda_pen)
    conc_ratio = k0_cord / k0_vessel if k0_vessel > 0 else 0.0

    if sanctuary_depth <= 0.0:
        summary = "Sanctuary absent: full-thickness therapeutic perfusion (Delta r = 0 um)"
    elif sanctuary_depth < 30.0:
        summary = f"Borderline sanctuary shell: Delta r = {sanctuary_depth:.1f} um (C_min/C0 = {conc_ratio:.3f})"
    else:
        summary = f"Substantial sanctuary mantle: Delta r = {sanctuary_depth:.1f} um (C_min/C0 = {conc_ratio:.3f})"

    return round(sanctuary_depth, 1), round(conc_ratio, 4), summary


def score_spatial(pattern: object) -> Tuple[float, str]:
    val = str(pattern or "").strip().lower()
    if "lepidic" in val:
        return 0.00, "Lepidic architecture (perivascular capillary contact)"
    if "solid" in val or "micropapillary" in val:
        return 1.00, "Solid/Micropapillary architecture (severe vascular decoupling)"
    if "mixed" in val:
        return 0.80, "Mixed architectural pattern (intermediate cord expansion)"
    if "acinar" in val or "papillary" in val:
        return 0.40, "Acinar/Papillary architecture (moderate cord radius)"
    return 0.50, f"Unclassified pattern: {pattern}"


def score_genomic(tp53: object, rb1: object) -> Tuple[float, str]:
    tp53_str = str(tp53 or "").strip().lower()
    rb1_str = str(rb1 or "").strip().lower()

    tp53_is_wt = any(x in tp53_str for x in ["wt", "wild", "neg", "normal"]) and not any(
        x in tp53_str for x in ["mut", "alter", "loss", "del", "pos"]
    )
    rb1_is_wt = any(x in rb1_str for x in ["wt", "wild", "neg", "normal"]) and not any(
        x in rb1_str for x in ["mut", "alter", "loss", "del", "pos"]
    )

    tp53_mut = not tp53_is_wt and any(x in tp53_str for x in ["mut", "alter", "pos", "1", "true", "yes", "p."])
    rb1_mut = not rb1_is_wt and any(x in rb1_str for x in ["mut", "alter", "loss", "del", "pos", "1", "true", "yes"])

    if not tp53_mut and not rb1_mut:
        return 0.00, "Genomic: TP53 wild-type / RB1 wild-type"
    if tp53_mut and not rb1_mut:
        return 0.50, "Genomic: TP53 altered / RB1 wild-type"
    if tp53_mut and rb1_mut:
        return 1.00, "Genomic: TP53 altered / RB1 altered (dual disruption)"
    return 0.30, "Genomic: Isolated RB1 alteration or indeterminate"


def score_stromal(stroma: object) -> Tuple[float, str]:
    val = str(stroma or "").strip().lower()
    if any(x in val for x in ["low", "loose", "minimal", "normal"]):
        return 0.00, "Stroma: Loose / minimal desmoplasia"
    if any(x in val for x in ["dense", "severe", "high", "sclerot", "fibrot"]):
        return 1.00, "Stroma: Dense collagen desmoplasia"
    if any(x in val for x in ["moderate", "intermediate", "mod"]):
        return 0.50, "Stroma: Moderate interstitial desmoplasia"
    return 0.50, f"Stroma: Custom descriptor: {stroma}"


def calculate_cri(
    pattern: object,
    tp53: object = "wt",
    rb1: object = "wt",
    stroma: object = "moderate",
    w_spatial: float = 0.47,
    w_genomic: float = 0.33,
    w_stromal: float = 0.20,
    cutoff_low: float = 0.30,
    cutoff_high: float = 0.65,
    lambda_pen: float = 34.1,
    r_vessel: float = 8.0,
) -> Dict[str, object]:
    """Calculate normalized Composite Resistance Index and biophysical sanctuary parameters."""
    s_score, s_desc = score_spatial(pattern)
    g_score, g_desc = score_genomic(tp53, rb1)
    m_score, m_desc = score_stromal(stroma)

    total_w = w_spatial + w_genomic + w_stromal
    if total_w <= 0:
        total_w = 1.0

    cri = (w_spatial * s_score + w_genomic * g_score + w_stromal * m_score) / total_w

    if cri < cutoff_low:
        tier = "Tier 1 (Low Sanctuary Risk)"
    elif cri <= cutoff_high:
        tier = "Tier 2 (Intermediate Risk)"
    else:
        tier = "Tier 3 (High Sanctuary Risk)"

    sanctuary_depth, conc_ratio, bio_summary = compute_biophysical_sanctuary(
        str(pattern), lambda_pen=lambda_pen, r_vessel=r_vessel
    )

    return {
        "S": round(s_score, 2),
        "G": round(g_score, 2),
        "M": round(m_score, 2),
        "CRI": round(cri, 3),
        "tier": tier,
        "risk_tier": tier,
        "sanctuary_depth_um": sanctuary_depth,
        "concentration_ratio": conc_ratio,
        "bio_summary": bio_summary,
        "desc_S": s_desc,
        "desc_G": g_desc,
        "desc_M": m_desc,
    }


def find_column(fieldnames: Sequence[str], candidates: Iterable[str]) -> Optional[str]:
    norm = {f.strip().lower(): f for f in fieldnames}
    for c in candidates:
        if c.lower() in norm:
            return norm[c.lower()]
    return None


def process_batch(
    input_path: str,
    output_path: Optional[str] = None,
    delimiter: Optional[str] = None,
    w_spatial: float = 0.47,
    w_genomic: float = 0.33,
    w_stromal: float = 0.20,
    cutoff_low: float = 0.30,
    cutoff_high: float = 0.65,
    lambda_pen: float = 34.1,
    r_vessel: float = 8.0,
) -> None:
    """Process high-throughput cohort table and append quantitative biophysical metrics."""
    if not os.path.exists(input_path):
        raise FileNotFoundError(input_path)

    if delimiter is None:
        delimiter = "," if input_path.lower().endswith(".csv") else "\t"

    with open(input_path, encoding="utf-8-sig", newline="") as src:
        reader = csv.DictReader(src, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError("Input file has no header row")

        fields = list(reader.fieldnames)
        col_pattern = find_column(fields, ["histologic_subtype", "pattern", "histology", "subtype", "growth_pattern"])
        if col_pattern is None:
            raise ValueError("Input requires a pattern or histology column (e.g., 'histologic_subtype')")

        col_tp53 = find_column(fields, ["tp53_status", "tp53", "tp53_mut", "p53"])
        col_rb1 = find_column(fields, ["rb1_status", "rb1", "rb1_mut"])
        col_stroma = find_column(fields, ["stromal_desmoplasia", "stroma", "stromal_density", "desmoplasia"])

        new_cols = [
            "S_score",
            "G_score",
            "M_score",
            "Sanctuary_Depth_um",
            "C_min_Ratio",
            "CRI_Score",
            "Risk_Tier",
            "Biophysical_Summary",
        ]
        out_fields = fields + new_cols

        dest: TextIO = open(output_path, "w", encoding="utf-8", newline="") if output_path else sys.stdout
        try:
            writer = csv.DictWriter(dest, fieldnames=out_fields, delimiter="\t")
            writer.writeheader()
            for row in reader:
                pat = row.get(col_pattern, "")
                tp = row.get(col_tp53, "wt") if col_tp53 else "wt"
                rb = row.get(col_rb1, "wt") if col_rb1 else "wt"
                st = row.get(col_stroma, "moderate") if col_stroma else "moderate"

                res = calculate_cri(
                    pat,
                    tp,
                    rb,
                    st,
                    w_spatial=w_spatial,
                    w_genomic=w_genomic,
                    w_stromal=w_stromal,
                    cutoff_low=cutoff_low,
                    cutoff_high=cutoff_high,
                    lambda_pen=lambda_pen,
                    r_vessel=r_vessel,
                )

                out_row = dict(row)
                out_row["S_score"] = f"{res['S']:.2f}"
                out_row["G_score"] = f"{res['G']:.2f}"
                out_row["M_score"] = f"{res['M']:.2f}"
                out_row["Sanctuary_Depth_um"] = f"{res['sanctuary_depth_um']:.1f}"
                out_row["C_min_Ratio"] = f"{res['concentration_ratio']:.4f}"
                out_row["CRI_Score"] = f"{res['CRI']:.3f}"
                out_row["Risk_Tier"] = res["tier"]
                out_row["Biophysical_Summary"] = res["bio_summary"]
                writer.writerow(out_row)
        finally:
            if output_path and dest is not sys.stdout:
                dest.close()


def print_single_report(res: Dict[str, object], args: argparse.Namespace) -> None:
    print("========================================================================")
    print("  COMPOSITE RESISTANCE INDEX (CRI) BIOPHYSICAL ANALYSIS")
    print("  Translational Reaction-Diffusion Sanctuary Modeling for Solid Tumors")
    print("  Research Use Only (RUO) - Quantitative Model")
    print("========================================================================")
    print(f"  [1] Spatial Architecture (S):   {res['S']:.2f}  ({res['desc_S']})")
    print(f"  [2] Genomic Co-Mutation (G):    {res['G']:.2f}  ({res['desc_G']})")
    print(f"  [3] Stromal Desmoplasia (M):    {res['M']:.2f}  ({res['desc_M']})")
    print("------------------------------------------------------------------------")
    print(f"  Model Weights Applied:          w_S={args.w_spatial:.2f}, w_G={args.w_genomic:.2f}, w_M={args.w_stromal:.2f}")
    print(f"  Effective Penetration Lambda:   {args.lambda_pen:.1f} um (vessel radius: {args.r_vessel:.1f} um)")
    print(f"  Perivascular Sanctuary Depth:   Delta r = {res['sanctuary_depth_um']:.1f} um")
    print(f"  Boundary Concentration Ratio:   C_min / C0 = {res['concentration_ratio']:.4f}")
    print("------------------------------------------------------------------------")
    print(f"  COMPOSITE RESISTANCE INDEX:     {res['CRI']:.3f} / 1.000")
    print(f"  STRATIFICATION TIER:            {res['tier']}")
    print(f"  BIOPHYSICAL SUMMARY:            {res['bio_summary']}")
    print("========================================================================")
    print("  Disclaimer: For research and computational exploration only.")
    print("  Downstream survival modeling (e.g. Kaplan-Meier / Cox) should be")
    print("  performed on cohort datasets using Galaxy survival packages.")
    print("========================================================================")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Biophysical Composite Resistance Index and Reaction-Diffusion Solver (RUO)."
    )
    parser.add_argument("--pattern", type=str, default="solid", help="Histologic growth pattern")
    parser.add_argument("--tp53", type=str, default="wt", help="TP53 status (wt, mutated)")
    parser.add_argument("--rb1", type=str, default="wt", help="RB1 status (wt, mutated)")
    parser.add_argument("--stroma", type=str, default="moderate", help="Stromal density (low, moderate, dense)")

    # Customizable weights & cutoffs
    parser.add_argument("--w-spatial", type=float, default=0.47, help="Weight for spatial decoupling (default: 0.47)")
    parser.add_argument("--w-genomic", type=float, default=0.33, help="Weight for genomic co-mutations (default: 0.33)")
    parser.add_argument("--w-stromal", type=float, default=0.20, help="Weight for stromal desmoplasia (default: 0.20)")
    parser.add_argument("--cutoff-low", type=float, default=0.30, help="Tier 1 upper threshold (default: 0.30)")
    parser.add_argument("--cutoff-high", type=float, default=0.65, help="Tier 3 lower threshold (default: 0.65)")

    # Biophysical parameters
    parser.add_argument("--lambda-pen", type=float, default=34.1, help="Penetration distance lambda in um (default: 34.1)")
    parser.add_argument("--r-vessel", type=float, default=8.0, help="Microvessel radius in um (default: 8.0)")

    # Batch mode
    parser.add_argument("--batch", type=str, help="Input cohort TSV/CSV file path")
    parser.add_argument("--out", type=str, help="Output annotated cohort TSV file path")
    parser.add_argument("--delimiter", type=str, help="Field delimiter (comma, tab)")

    args = parser.parse_args()

    if args.batch:
        process_batch(
            args.batch,
            args.out,
            delimiter=args.delimiter,
            w_spatial=args.w_spatial,
            w_genomic=args.w_genomic,
            w_stromal=args.w_stromal,
            cutoff_low=args.cutoff_low,
            cutoff_high=args.cutoff_high,
            lambda_pen=args.lambda_pen,
            r_vessel=args.r_vessel,
        )
    else:
        res = calculate_cri(
            args.pattern,
            args.tp53,
            args.rb1,
            args.stroma,
            w_spatial=args.w_spatial,
            w_genomic=args.w_genomic,
            w_stromal=args.w_stromal,
            cutoff_low=args.cutoff_low,
            cutoff_high=args.cutoff_high,
            lambda_pen=args.lambda_pen,
            r_vessel=args.r_vessel,
        )
        print_single_report(res, args)


if __name__ == "__main__":
    main()
