"""PNG copies of the MIC-50-90 report charts, and a short answer table.

Reads results.json written by "mic-50-90 distribution". Each chart is drawn
by the program itself (mic_50_90.result_view.count_chart, the chart of the
HTML report) and only converted to PNG here, with the DejaVu Sans fonts that
the program ships. Every number in the table is copied from the results.
"""
import argparse
import json
import os
import re
import tempfile
from pathlib import Path

import mic_50_90
from fontTools.ttLib import TTFont
from mic_50_90.result_view import category_label, count_chart

FONTS = Path(mic_50_90.__file__).parent / "report_assets" / "fonts"
ANSWERS = {"supported": "yes", "contradicted": "no", "undetermined": "not enough information",
           "available": ""}
CRITERIA = {"<": "fewer than", "<=": "at most", ">": "more than", ">=": "at least"}


def use_program_fonts(directory):
    """Make fontconfig see only the DejaVu Sans fonts shipped with the program."""
    for woff in FONTS.glob("DejaVuSans*.woff"):
        font = TTFont(woff)
        font.flavor = None
        font.save(directory / (woff.stem + ".ttf"))
    conf = directory / "fonts.conf"
    conf.write_text(
        '<?xml version="1.0"?><!DOCTYPE fontconfig SYSTEM "fonts.dtd"><fontconfig>'
        f"<dir>{directory}</dir><cachedir>{directory / 'cache'}</cachedir></fontconfig>",
        encoding="utf-8")
    os.environ["FONTCONFIG_FILE"] = str(conf)


def count_text(components):
    return " or ".join(str(lo) if lo == hi else f"{lo}–{hi}" for lo, hi in components)


def percent_text(components, n):
    return " or ".join(f"{100 * lo / n:.2f}%" if lo == hi else f"{100 * lo / n:.2f}–{100 * hi / n:.2f}%"
                       for lo, hi in components)


def question_text(decision):
    return re.sub(r"^fraction with (recorded )?MIC ", "", decision.get("estimand", ""))


def criterion_text(decision):
    if not decision.get("decision_operator"):
        return ""
    words = CRITERIA.get(decision["decision_operator"], decision["decision_operator"])
    return f"{words} {100 * float(decision['decision_fraction']):g}%"


def through_text(cohort, label):
    """The suggested count covers every recorded group up to and including this one."""
    for category in cohort.get("panel", {}).get("categories", []):
        if category.get("label") == label and category.get("upper_bound") is not None:
            sign = "≤" if category.get("upper_closed") else "<"
            return f"isolates with MIC {sign} {category['upper_bound']:g} {cohort.get('unit', '')}".rstrip()
    return f"isolates in the MIC groups up to {category_label(label)}"


def cell(value):
    return re.sub(r"[\t\r\n]+", " ", str(value))


def chart_name(cohort_id, used):
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", cohort_id).strip("._") or "cohort"
    candidate, number = name, 1
    while candidate in used:
        number += 1
        candidate = f"{name}_{number}"
    used.add(candidate)
    return candidate


def summary_rows(cohort):
    name = cohort["cohort_id"]
    if cohort.get("status") != "ok":
        yield [name, "", "", "", "", "", cohort.get("status", ""), "", "", "", cohort.get("reason", "")]
        return
    n = cohort["n"]
    request = ""
    suggested = cohort.get("interpretation", {}).get("suggested_count")
    if suggested:
        request = (f"{through_text(cohort, suggested['through_label'])}: now "
                   f"{count_text(suggested['count_components'])}")
    note = ""
    rejected = [v.get("id", "") for v in cohort.get("rejected_variants", []) if isinstance(v, dict)]
    if rejected:
        kept = [v.get("id", "") if isinstance(v, dict) else str(v) for v in cohort.get("retained_variants", [])]
        note = f"kept variant {', '.join(kept)}; rejected {', '.join(rejected)}"
    decisions = cohort.get("decisions", [])
    if not decisions:
        yield [name, n, "", "", "", "", "", "", "", request, "; ".join(x for x in (note, "no question supplied") if x)]
    for decision in decisions:
        components = decision["sample"].get("count_components") or []
        population = decision.get("population", {})
        bounds = answer = ""
        if "lower" in population and "upper" in population:
            bounds = f"{100 * population['lower']:.2f}–{100 * population['upper']:.2f}%"
            answer = ANSWERS.get(population.get("status"), population.get("status", ""))
        yield [name, n, question_text(decision), count_text(components),
               percent_text(components, n) if components else "", criterion_text(decision),
               ANSWERS.get(decision["sample"]["status"], decision["sample"]["status"]),
               bounds, answer, request if decision["sample"]["status"] == "undetermined" else "", note]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", help="results.json written by mic-50-90 distribution")
    parser.add_argument("--charts", required=True, help="directory for the PNG charts")
    parser.add_argument("--summary", required=True, help="tab-separated answer table")
    parser.add_argument("--width", type=int, default=1600, help="width of each PNG in pixels")
    parser.add_argument("--no-header", dest="header", action="store_false",
                        help="omit the header line (Galaxy sets the column names as metadata)")
    args = parser.parse_args()
    cohorts = json.loads(Path(args.results).read_text(encoding="utf-8"))["cohorts"]
    lines = [["cohort", "isolates", "recorded MIC", "possible count", "possible share", "criterion",
              "answer for the sample", "population share", "answer for the population", "count worth requesting", "note"]]
    for cohort in cohorts:
        lines.extend(summary_rows(cohort))
    Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
    if not args.header:
        lines = lines[1:]
    Path(args.summary).write_text("".join("\t".join(map(cell, line)) + "\n" for line in lines), encoding="utf-8")
    charts = Path(args.charts)
    charts.mkdir(parents=True, exist_ok=True)
    used = set()
    with tempfile.TemporaryDirectory() as fonts:
        use_program_fonts(Path(fonts))
        import cairosvg
        for cohort in cohorts:
            if cohort.get("status") != "ok":
                continue
            svg = count_chart(cohort["sample"]["categories"], cohort["n"])
            name = chart_name(cohort["cohort_id"], used)
            cairosvg.svg2png(bytestring=svg.encode("utf-8"), write_to=str(charts / f"{name}.png"),
                             output_width=args.width)


if __name__ == "__main__":
    main()
