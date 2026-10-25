#!/usr/bin/env python3
"""
cost_of_off_110.py  (version 1)
===============================
Series 110 -- what does leaving channel loss OFF cost when the truth HAS channel loss?
READ-ONLY. Runs no tRIBS. Reads two saved tables and writes new files into its own folder.

WHAT THIS DOES (plain language)
--------------------------------
Part A already made 27 "truths" that include channel loss (cc ON), and already scored
two sets of 250 saved runs against each truth:
    "cc ON"  -- runs that had channel loss switched on (the matching structure)
    "cc OFF" -- runs with channel loss switched off (the structure you are considering)
This script asks, for each of the 27 truths: if I may only use the cc-OFF structure, how
well can it do, and how wrong are the answers?  For the BEST cc-OFF run of each truth it
reports

    - best KGE_2012 (and the same for the cc-ON runs, as a reference)
    - how far its Ks and f sit from the known truth values
    - volume error (PBIAS), peak error, peak-time error, rising-limb ratio,
      duration error, KGE beta and gamma
    - how many runs fit nearly as well (a rough measure of how wide the Ks-f "swoosh" is)
    - the share of runs with volume error under 2% ("feasible")

then counts how often the cc-OFF errors exceed the tolerances below, and breaks the
counts down by the true cc, Ks and f (the conditions the errors appear at).

The cc-ON numbers are printed beside the cc-OFF numbers on purpose. Both come from a
250-run sample, so some error is only sampling coarseness. If cc ON "exceeds" too, that
part is not caused by leaving loss off.

WHERE TO RUN IT
----------------
From workspaces/SMF/lab/ :
    python cost_of_off_110.py --check_only
    python cost_of_off_110.py

Defaults read the Series 110 truth-location tables (the original 12 Aug 2014 storm):
    ../calibration_work/03_comparisons/summary_tables/rescored_truth_location_110/
        rescored_location_summary_110.csv
        rescored_location_long_110.csv
and write to
    ../calibration_work/03_comparisons/summary_tables/cost_of_off_110/<label>/
The storm-size versions (x0.8 and x1.25 rain) can be read by pointing --summary, --long
and --label at those tables, if they have the same columns.

TOLERANCES (mine -- not from the literature, fixed before any result was looked at)
------------------------------------------------------------------------------------
    volume     |PBIAS|            above 10 %      (--vol_tol)
    peak       |peak error|       above 10 %      (--peak_tol)
    timing     |peak time error|  above 0.25 h    (--timing_tol_hr; 3 steps of 5 min)
    KGE gap    best cc ON minus best cc OFF above 0.03   (--gap_tol)
    near-best  within 0.01 KGE_2012 of the best (--near_tol)
    feasible   |PBIAS| under 2 % (--feas_pct), the Series 110 definition
Change them with the options; the values used are printed and saved in the PROVENANCE file.
Do not change a tolerance after looking at the result.

SIGN CONVENTIONS (from the saved tables)
-----------------------------------------
PBIAS, volume error and peak error are (model minus truth); positive = model too large.
Peak time error is model minus truth; positive = model peaks later.

GATES (the script stops and writes nothing if one fails)
---------------------------------------------------------
  1. Required columns exist in both tables, and the series names are "cc OFF" and "cc ON".
  2. Every (truth, series) in the summary table has the same number of rows in the long
     table as its n_scored.
  3. The best run recomputed from the long table has the same KGE_2012 as the summary
     table says (to 1e-9), and the same Ks/f offsets from the truth.

CAVEATS PRINTED WITH EVERY RUN
-------------------------------
  - Twin experiment: the only structural difference is the missing channel loss, so the
    costs are lower bounds on what real data would show.
  - Sampling: each best run is the best of 250 points in a box (Ks 4 to 10.5, f 0.004 to
    0.030). Truths near a box edge cannot show an honest shift; compare with the
    null medians in location_per_truth_110.csv.
  - Routing is pinned at the truth values; gauge error is not quantified; one storm size
    per run of this script.
  - Tolerances are mine.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

VERSION = "1"
SERIES_OFF = "cc OFF"
SERIES_ON = "cc ON"

REQ_SUMMARY = ["truth_id", "truth_Ks", "truth_f", "truth_cc", "series", "n_scored",
               "best_kge_2012", "best_run_id", "best_Ks", "best_f", "best_pbias_pct",
               "best_dKs_from_truth", "best_dlog10_f_from_truth"]
REQ_LONG = ["truth_id", "truth_Ks", "truth_f", "truth_cc", "series", "run_id", "Ks_mult",
            "f_RS_abs", "pbias_pct", "peak_error_pct", "peak_timing_error_hr", "kge_beta",
            "kge_gamma", "kge_2012", "rising_limb_steepness_ratio",
            "duration_above_thresh_error_min"]
EVENT = ["pbias_pct", "peak_error_pct", "peak_timing_error_hr", "rising_limb_steepness_ratio",
         "duration_above_thresh_error_min", "kge_beta", "kge_gamma"]

SUMMARY_REL = ("calibration_work", "03_comparisons", "summary_tables",
               "rescored_truth_location_110", "rescored_location_summary_110.csv")
LONG_REL = ("calibration_work", "03_comparisons", "summary_tables",
            "rescored_truth_location_110", "rescored_location_long_110.csv")
OUT_REL = ("calibration_work", "03_comparisons", "summary_tables", "cost_of_off_110")


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fail(msg):
    sys.exit("STOP: " + msg)


def best_table(long, series):
    """Best run (highest KGE_2012) for each truth in one series."""
    d = long[long["series"] == series].dropna(subset=["kge_2012"])
    idx = d.groupby("truth_id")["kge_2012"].idxmax()
    return d.loc[idx].set_index("truth_id")


def run_gates(summary, long):
    print("\nGATES")
    # gate 1
    miss_s = [c for c in REQ_SUMMARY if c not in summary.columns]
    miss_l = [c for c in REQ_LONG if c not in long.columns]
    if miss_s or miss_l:
        fail("missing columns. summary: %s   long: %s" % (miss_s, miss_l))
    for name, df in (("summary", summary), ("long", long)):
        got = set(df["series"].unique())
        if got != {SERIES_OFF, SERIES_ON}:
            fail("%s table has series %s, expected %s" % (name, sorted(got), [SERIES_OFF, SERIES_ON]))
    print("  gate 1 PASS: columns present; series are '%s' and '%s'" % (SERIES_OFF, SERIES_ON))

    # gate 2 and 3
    problems = []
    n_checked = 0
    worst_kge = 0.0
    worst_f = 0.0
    worst_ks = 0.0
    id_diff = 0
    for _, s in summary.iterrows():
        g = long[(long["truth_id"] == s["truth_id"]) & (long["series"] == s["series"])]
        if len(g) != int(s["n_scored"]):
            problems.append("%s / %s: %d rows in long table, n_scored says %d"
                            % (s["truth_id"], s["series"], len(g), int(s["n_scored"])))
            continue
        gv = g.dropna(subset=["kge_2012"])
        if gv.empty:
            problems.append("%s / %s: no valid KGE_2012" % (s["truth_id"], s["series"]))
            continue
        b = gv.loc[gv["kge_2012"].idxmax()]
        d_kge = abs(b["kge_2012"] - s["best_kge_2012"])
        d_ks = abs((b["Ks_mult"] - s["truth_Ks"]) - s["best_dKs_from_truth"])
        d_f = abs(np.log10(b["f_RS_abs"] / s["truth_f"]) - s["best_dlog10_f_from_truth"])
        worst_kge, worst_ks, worst_f = max(worst_kge, d_kge), max(worst_ks, d_ks), max(worst_f, d_f)
        if b["run_id"] != s["best_run_id"]:
            id_diff += 1
        if d_kge > 1e-9 or d_ks > 1e-9 or d_f > 1e-9:
            problems.append("%s / %s: recomputed best differs from summary (dKGE %.2e, dKs %.2e, dlogf %.2e)"
                            % (s["truth_id"], s["series"], d_kge, d_ks, d_f))
        n_checked += 1
    if problems:
        for p in problems[:15]:
            print("  FAIL:", p)
        fail("%d gate problem(s); nothing written." % len(problems))
    print("  gate 2 PASS: row counts match n_scored for all %d (truth, series) pairs" % n_checked)
    print("  gate 3 PASS: recomputed best KGE_2012 / Ks / f offsets match the summary "
          "(worst %.1e / %.1e / %.1e)" % (worst_kge, worst_ks, worst_f))
    if id_diff:
        print("  note: %d best run ids differ from the summary's (tie in KGE_2012 within 1e-9)" % id_diff)


def near_and_feasible(long, series, best_kge, near_tol, feas_pct):
    """Per truth: number of near-best runs, their spread in f (decades) and Ks, and % feasible."""
    out = {}
    d = long[long["series"] == series]
    for tid, g in d.groupby("truth_id"):
        gv = g.dropna(subset=["kge_2012"])
        n = gv[gv["kge_2012"] >= best_kge[tid] - near_tol]
        out[tid] = dict(
            n_near=len(n),
            f_decades=float(np.log10(n["f_RS_abs"].max() / n["f_RS_abs"].min())) if len(n) > 1 else 0.0,
            ks_span=float(n["Ks_mult"].max() - n["Ks_mult"].min()) if len(n) > 1 else 0.0,
            pct_feasible=float((gv["pbias_pct"].abs() < feas_pct).mean() * 100.0))
    return out


def build_per_truth(long, args):
    best_off = best_table(long, SERIES_OFF)
    best_on = best_table(long, SERIES_ON)
    near_off = near_and_feasible(long, SERIES_OFF, best_off["kge_2012"], args.near_tol, args.feas_pct)
    near_on = near_and_feasible(long, SERIES_ON, best_on["kge_2012"], args.near_tol, args.feas_pct)

    truths = (long[["truth_id", "truth_Ks", "truth_f", "truth_cc"]].drop_duplicates()
              .set_index("truth_id").sort_values(["truth_cc", "truth_Ks", "truth_f"]))
    rows = []
    for tid, t in truths.iterrows():
        r = dict(truth_id=tid, truth_Ks=t["truth_Ks"], truth_f=t["truth_f"], truth_cc=t["truth_cc"])
        for tag, b, nf in (("off", best_off, near_off), ("on", best_on, near_on)):
            x = b.loc[tid]
            r[tag + "_best_kge"] = x["kge_2012"]
            r[tag + "_best_Ks"] = x["Ks_mult"]
            r[tag + "_dKs"] = x["Ks_mult"] - t["truth_Ks"]
            r[tag + "_best_f"] = x["f_RS_abs"]
            r[tag + "_f_factor"] = x["f_RS_abs"] / t["truth_f"]
            for c in EVENT:
                r["%s_%s" % (tag, c)] = x[c]
            r[tag + "_n_near_best"] = nf[tid]["n_near"]
            r[tag + "_near_f_decades"] = nf[tid]["f_decades"]
            r[tag + "_near_Ks_span"] = nf[tid]["ks_span"]
            r[tag + "_pct_feasible"] = nf[tid]["pct_feasible"]
            r[tag + "_run_id"] = x["run_id"]
        r["kge_gap_on_minus_off"] = r["on_best_kge"] - r["off_best_kge"]
        for tag in ("off", "on"):
            r[tag + "_x_volume"] = abs(r[tag + "_pbias_pct"]) > args.vol_tol
            r[tag + "_x_peak"] = abs(r[tag + "_peak_error_pct"]) > args.peak_tol
            r[tag + "_x_timing"] = abs(r[tag + "_peak_timing_error_hr"]) > args.timing_tol_hr
            r[tag + "_x_any"] = bool(r[tag + "_x_volume"] or r[tag + "_x_peak"] or r[tag + "_x_timing"])
        r["off_x_gap"] = r["kge_gap_on_minus_off"] > args.gap_tol
        rows.append(r)
    return pd.DataFrame(rows)


def by_condition(pt):
    """Median gap and counts that exceed tolerance, by each true parameter level."""
    frames = []
    for col in ("truth_cc", "truth_Ks", "truth_f"):
        g = pt.groupby(col)
        t = pd.DataFrame({
            "n_truths": g.size(),
            "median_kge_gap": g["kge_gap_on_minus_off"].median(),
            "median_off_dKs": g["off_dKs"].median(),
            "median_off_pbias": g["off_pbias_pct"].median(),
            "median_off_peak_err": g["off_peak_error_pct"].median(),
            "n_off_x_volume": g["off_x_volume"].sum(),
            "n_off_x_peak": g["off_x_peak"].sum(),
            "n_off_x_timing": g["off_x_timing"].sum(),
            "n_off_x_gap": g["off_x_gap"].sum(),
            "n_off_x_any": g["off_x_any"].sum(),
            "n_on_x_any": g["on_x_any"].sum()})
        t.insert(0, "level", t.index)
        t.insert(0, "by", col)
        frames.append(t.reset_index(drop=True))
    out = pd.concat(frames, ignore_index=True)
    out["level"] = out["level"].map(lambda v: "%g" % v)
    return out


SHORT = {"n_truths": "n", "median_kge_gap": "gap_med", "median_off_dKs": "dKs_med",
         "median_off_pbias": "pbias_med", "median_off_peak_err": "peak_med",
         "n_off_x_volume": "OFF_vol", "n_off_x_peak": "OFF_peak", "n_off_x_timing": "OFF_time",
         "n_off_x_gap": "OFF_gap", "n_off_x_any": "OFF_any", "n_on_x_any": "ON_any"}


def fmt(df, nd=3):
    return df.rename(columns=SHORT).round(nd).to_string(index=False)


def report(pt, bc, args):
    n = len(pt)
    lines = []
    P = lines.append
    P("=" * 78)
    P("Series 110: cost of leaving channel loss OFF when the truth has it (version %s)" % VERSION)
    P("label %s   truths %d" % (args.label, n))
    P("tolerances (mine): volume %.0f%%, peak %.0f%%, timing %.2f h, KGE gap %.2f, near-best %.3f, feasible %.1f%%"
      % (args.vol_tol, args.peak_tol, args.timing_tol_hr, args.gap_tol, args.near_tol, args.feas_pct))
    P("=" * 78)

    P("\n1. BEST KGE_2012 FOR EACH STRUCTURE (cc OFF = no loss in the model, cc ON = loss in the model)")
    s = pd.DataFrame({"min": [pt.off_best_kge.min(), pt.on_best_kge.min(), pt.kge_gap_on_minus_off.min()],
                      "median": [pt.off_best_kge.median(), pt.on_best_kge.median(), pt.kge_gap_on_minus_off.median()],
                      "max": [pt.off_best_kge.max(), pt.on_best_kge.max(), pt.kge_gap_on_minus_off.max()]},
                     index=["cc OFF best", "cc ON best", "gap (ON - OFF)"])
    P(s.round(3).to_string())
    P("Truths where the gap is within %.2f: %d of %d" % (args.gap_tol, int((~pt.off_x_gap).sum()), n))

    P("\n2. SHIFT OF THE BEST RUN FROM THE KNOWN Ks AND f (median / min / max over truths)")
    sh = pd.DataFrame({
        "cc OFF": [pt.off_dKs.median(), pt.off_dKs.min(), pt.off_dKs.max(),
                   pt.off_f_factor.median(), pt.off_f_factor.min(), pt.off_f_factor.max()],
        "cc ON": [pt.on_dKs.median(), pt.on_dKs.min(), pt.on_dKs.max(),
                  pt.on_f_factor.median(), pt.on_f_factor.min(), pt.on_f_factor.max()]},
        index=["Ks shift, median", "Ks shift, min", "Ks shift, max",
               "f factor, median", "f factor, min", "f factor, max"])
    P(sh.round(3).to_string())
    P("(Ks shift = best Ks minus true Ks, in multiplier units. f factor = best f / true f; 1 = exact.)")
    P("Truths where cc OFF picked a HIGHER Ks than the truth: %d of %d; a higher f: %d of %d"
      % (int((pt.off_dKs > 0).sum()), n, int((pt.off_f_factor > 1).sum()), n))

    P("\n3. EVENT METRICS OF THE BEST cc-OFF RUN (median / min / max over truths; cc ON median for reference)")
    ev = []
    for c, lab in (("pbias_pct", "volume error, PBIAS (%)"), ("peak_error_pct", "peak error (%)"),
                   ("peak_timing_error_hr", "peak time error (h)"),
                   ("rising_limb_steepness_ratio", "rising-limb ratio (1 = same)"),
                   ("duration_above_thresh_error_min", "duration error (min)"),
                   ("kge_beta", "KGE beta (volume)"), ("kge_gamma", "KGE gamma (variability)")):
        ev.append([lab, pt["off_" + c].median(), pt["off_" + c].min(), pt["off_" + c].max(), pt["on_" + c].median()])
    P(pd.DataFrame(ev, columns=["metric", "OFF median", "OFF min", "OFF max", "ON median"]).round(3).to_string(index=False))

    P("\n4. HOW OFTEN THE BEST RUN EXCEEDS A TOLERANCE (of %d truths)" % n)
    ex = pd.DataFrame({
        "cc OFF": [int(pt.off_x_volume.sum()), int(pt.off_x_peak.sum()), int(pt.off_x_timing.sum()),
                   int(pt.off_x_gap.sum()), int(pt.off_x_any.sum())],
        "cc ON (sampling reference)": [int(pt.on_x_volume.sum()), int(pt.on_x_peak.sum()),
                                       int(pt.on_x_timing.sum()), "-", int(pt.on_x_any.sum())]},
        index=["volume", "peak", "timing", "KGE gap", "any of volume/peak/timing"])
    P(ex.to_string())

    P("\n5. DOES THE Ks-f 'SWOOSH' MOVE OR WIDEN? (median over truths)")
    sw = pd.DataFrame({
        "cc OFF": [pt.off_pct_feasible.median(), pt.off_n_near_best.median(),
                   pt.off_near_f_decades.median(), pt.off_near_Ks_span.median()],
        "cc ON": [pt.on_pct_feasible.median(), pt.on_n_near_best.median(),
                  pt.on_near_f_decades.median(), pt.on_near_Ks_span.median()]},
        index=["%% of runs with |PBIAS| < %.0f%%" % args.feas_pct, "runs within %.3f of best" % args.near_tol,
               "f spread of those runs (decades)", "Ks spread of those runs"])
    P(sw.round(2).to_string())
    P("(Moves: section 2. Widens: more near-best runs, a larger f or Ks spread, a larger feasible share.)")

    P("\n6. WHERE THE cc-OFF ERRORS APPEAR (by true cc, Ks and f)")
    P("   gap_med = median KGE gap (ON - OFF); dKs_med = median Ks shift of the best OFF run;")
    P("   pbias_med, peak_med = median volume and peak error of the best OFF run (%);")
    P("   OFF_vol, OFF_peak, OFF_time, OFF_gap, OFF_any = number of truths above each tolerance;")
    P("   ON_any = same 'any of volume/peak/timing' count for the cc-ON reference.")
    P(fmt(bc, 2))

    P("\nCAVEATS")
    P("  - Twin experiment: the only structural difference is the missing loss, so these costs are lower bounds.")
    P("  - Each best run is the best of 250 points; truths near a box edge cannot show an honest shift.")
    P("  - Routing pinned at the truth values; gauge error not quantified; one storm size per run of this script.")
    P("  - The tolerances are the author's, fixed before the result was looked at.")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description="Series 110: cost of leaving channel loss OFF (read-only).")
    ap.add_argument("--summary", type=Path, default=None, help="rescored summary CSV (default: Series 110 truth-location table)")
    ap.add_argument("--long", type=Path, default=None, help="rescored long CSV (default: Series 110 truth-location table)")
    ap.add_argument("--label", default="loc_x1p0", help="name for this run's output folder (default loc_x1p0)")
    ap.add_argument("--out_dir", type=Path, default=None)
    ap.add_argument("--vol_tol", type=float, default=10.0)
    ap.add_argument("--peak_tol", type=float, default=10.0)
    ap.add_argument("--timing_tol_hr", type=float, default=0.25)
    ap.add_argument("--gap_tol", type=float, default=0.03)
    ap.add_argument("--near_tol", type=float, default=0.01)
    ap.add_argument("--feas_pct", type=float, default=2.0)
    ap.add_argument("--check_only", action="store_true", help="run the gates and stop")
    ap.add_argument("--overwrite", action="store_true", help="allow writing into an existing output folder")
    args = ap.parse_args()

    root = Path.cwd().parent
    summary_path = args.summary or root.joinpath(*SUMMARY_REL)
    long_path = args.long or root.joinpath(*LONG_REL)
    out_dir = args.out_dir or root.joinpath(*OUT_REL) / args.label

    print("cost_of_off_110.py version %s   python %s   pandas %s   numpy %s"
          % (VERSION, sys.version.split()[0], pd.__version__, np.__version__))
    for p in (summary_path, long_path):
        if not Path(p).exists():
            fail("input not found: %s" % p)
    for p in (summary_path, long_path):
        if out_dir.resolve() == Path(p).resolve().parent or Path(p).resolve().parent in out_dir.resolve().parents:
            fail("output folder %s sits inside an input folder; refusing." % out_dir)
    print("summary:", summary_path)
    print("long   :", long_path)

    summary = pd.read_csv(summary_path)
    long = pd.read_csv(long_path)
    run_gates(summary, long)
    if args.check_only:
        print("\n--check_only: gates passed, nothing written.")
        return

    if out_dir.exists() and any(out_dir.iterdir()) and not args.overwrite:
        fail("%s already has files; use --overwrite or another --label." % out_dir)

    pt = build_per_truth(long, args)
    bc = by_condition(pt)
    text = report(pt, bc, args)
    print("\n" + text)

    out_dir.mkdir(parents=True, exist_ok=True)
    pt.to_csv(out_dir / "cost_of_off_per_truth_110.csv", index=False)
    bc.to_csv(out_dir / "cost_of_off_by_condition_110.csv", index=False)
    (out_dir / "cost_of_off_report_110.txt").write_text(text + "\n")
    prov = dict(script="cost_of_off_110.py", version=VERSION, label=args.label,
                script_md5=md5(__file__),
                inputs={"summary": str(summary_path), "summary_md5": md5(summary_path),
                        "long": str(long_path), "long_md5": md5(long_path)},
                tolerances=dict(vol_tol=args.vol_tol, peak_tol=args.peak_tol, timing_tol_hr=args.timing_tol_hr,
                                gap_tol=args.gap_tol, near_tol=args.near_tol, feas_pct=args.feas_pct),
                n_truths=int(len(pt)), python=sys.version.split()[0], pandas=pd.__version__, numpy=np.__version__)
    (out_dir / "PROVENANCE_cost_of_off_110.json").write_text(json.dumps(prov, indent=2) + "\n")
    print("\nWrote 4 files to", out_dir)


if __name__ == "__main__":
    main()
