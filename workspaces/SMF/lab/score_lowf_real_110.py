"""
score_lowf_real_110.py   (version 1)
====================================
Series 110 -- score the runs of a LOW-f sweep (run_lowf_sweep_110.py) against the REAL SMF gauge,
for ONE series at a time (cc OFF or cc ON), and answer one question:

    Is the best real-data fit now INSIDE the f range, or is it still at an edge of the box?

No tRIBS is started. This script only reads the stored compare CSVs, exactly as
rescore_real_gauge_110.py (version 2, the verified Stage 0 script) does, and it REUSES that script's
reading, filling, gate and metric code unchanged (it imports it; keep it in this folder together with
rescore_cc_truths_110.py, rescore_truth_location_110.py and run_sensitivity_single_interp_tribs6.py).
For the ON-versus-OFF comparison on a corrected box, run rescore_real_gauge_110.py itself later (see
"AFTER THIS" below); it needs both series.

THE TWO ARMS (same as Stage 0)
-------------------------------
  filled  all 241 five-minute bins; gaps in the event-triggered gauge are filled from the record (line up
          to 60 min, then hold the last reading up to 24 h).  PRIMARY (your decision of 2026-10-04).
  as-is   only the bins that hold a gauge record (23 of 241), what the scorer does in real-gauge mode.
          A sensitivity check.

THE READING RULE (fixed here, BEFORE any result of the new sweep has been seen)
--------------------------------------------------------------------------------
Position in the f box, u = (log f - log f_floor) / (log f_ceiling - log f_floor): 0 at the floor, 1 at the
ceiling. The verdict uses the PRIMARY arm; the other arm gets the same reading as a check.
  1. best run has u above %(ceil).2f                                   -> AT THE CEILING
  2. best run has u below %(floor).2f, or at least %(nn)d of the 5 best runs have u below %(nf).2f -> AT THE FLOOR
  3. else, if the best of the f bands (equal widths in log f) is NOT the lowest band and the lowest band's best
     KGE_2012 is at least %(tol).2f below it                         -> BRACKETED (fit rises, then falls toward the floor)
  4. else                                                             -> UNSETTLED (best run inside, but the curve is flat to
                                                                         within %(tol).2f toward the floor)
The numbers are the constants AT_CEIL_U, AT_FLOOR_U, NEAR_FLOOR_U, NEAR_FLOOR_COUNT and TURN_TOL below; they are
rules of thumb of mine. KGE_2012 gaps under %(tol).2f are treated as ties (the same tolerance as
rescore_cc_truths_110.py). A partial sweep (fewer runs than the design) can be scored; it is then labelled
PARTIAL, and the first runs of a design are an arbitrary subset, not a stratified one.

CHECKS ("gates", run before anything is scored; same code as Stage 0)
----------------------------------------------------------------------
  GATE 1   the real series (workbook read, duplicates, values, peak scale, cadence) and 1b the filled series.
  GATE 2A  every run's stored Observed column equals the synthetic truth in synth_truth/ as read here.
  GATE 2B  every run's stored metrics are reproduced by the metric code here.
  GATE 4   all runs share one regular 5-minute grid and every scored real bin is on it.
  Also: every run lies inside the box in the design file (a wrong box would give a wrong verdict).
  GATE 3 (the actual scorer in real-gauge mode on one point) is NOT repeated: it PASSED on 2026-10-04 for
  this pipeline, and these runs use the same scorer, builder and gauge reading, unchanged.

USAGE (from lab/)
------------------
    python score_lowf_real_110.py --label lf1 --series OFF --check_only   # gates only; writes nothing
    python score_lowf_real_110.py --label lf1 --series OFF                # the full scoring
    python score_lowf_real_110.py --label lf1 --series OFF --bands 8      # more f bands (default 6)
  It only READS CSVs, so it is safe to run while the sweep is still going (the sweep writes its CSV atomically);
  expect the result to be labelled PARTIAL.

WHAT YOU GET   (calibration_work/03_comparisons/summary_tables/lowf_real_gauge_110/<label>_<series>/)
  lowf_gates_110.csv   PROVENANCE_lowf_110.json
  filled/ and asis/ :  lowf_long_110.csv (every run, all metrics, position in the box)
                       lowf_top_runs_110.csv   lowf_f_bands_110.csv   lowf_verdict_110.csv
                       fig_lowf_kge_vs_f_110.png   fig_lowf_Ks_vs_f_110.png   fig_lowf_hydrograph_110.png
  The old sweep (Stage 0, f 0.004-0.030) is drawn alongside when real_gauge_110/<arm>/real_gauge_long_110.csv exists.

AFTER THIS
-----------
  BRACKETED   build the ON-versus-OFF box around the optimum: run_lowf_sweep_110.py with --mode off and --mode on,
              same --label, --n, --seed and box, then
                  python rescore_real_gauge_110.py --margin 0.05 --primary_arm filled --expect_n N \\
                      --stage2_csv <ON csv> --control_csv <OFF csv> --out_dir <new folder>
              (the ON and OFF points must be the same points; run_lowf_sweep_110.py guarantees that).
  AT THE FLOOR  lower the floor and sweep again under a new --label.
"""

import argparse
import json
import os
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

try:
    import rescore_real_gauge_110 as RG
    import rescore_cc_truths_110 as RS
    import rescore_truth_location_110 as TL
except ImportError as exc:
    sys.exit("This script needs rescore_real_gauge_110.py (version 2), rescore_cc_truths_110.py and "
             "rescore_truth_location_110.py in the same folder (it reuses their code): %s. Run from lab/." % exc)

for _name in ("read_gauge_workbook", "gauge_to_grid", "check_real_series", "build_filled", "check_filled",
              "check_grid", "Checks", "print_gates", "_style", "C_F", "C_KS", "NULLC", "INK", "INK2", "SURFACE", "GRID"):
    if not hasattr(RG, _name):
        sys.exit("rescore_real_gauge_110.py in this folder is not version 2 (missing %s). Replace it with the "
                 "current copy." % _name)

# ------------------------------------------------------------------
# Constants. The reading rule is mine and is fixed before any result is seen.
# ------------------------------------------------------------------
AT_CEIL_U = 0.90
AT_FLOOR_U = 0.10
NEAR_FLOOR_U = 0.15
NEAR_FLOOR_COUNT = 3        # of the 5 best runs
TURN_TOL = 0.01             # KGE_2012 gaps below this are ties (same value as rescore_cc_truths_110.DEFAULT_KGE_TOL)
OUT_DIRNAME = "lowf_real_gauge_110"
OLD_DIRNAME = "real_gauge_110"
SERIES_TEXT = {"OFF": "cc OFF", "ON": "cc ON"}
PERC = {"OFF": 0, "ON": 1}

VERDICTS = {
    "CEILING": "The best run sits at the TOP of this f box. That is unexpected after Stage 0 (which pointed to low f): "
               "check the design file, the Ks box and the figure before reading anything into it.",
    "FLOOR": "The best fit is at the LOWEST f of this box. The optimum is not bracketed: lower the floor and sweep "
             "again under a new --label. Do not read the channel-loss question from this box.",
    "BRACKETED": "The fit rises and then falls with f inside this box: the real-data optimum is bracketed. Build the "
                 "ON-versus-OFF box around it (see the top of this file).",
    "UNSETTLED": "The best run is inside the box, but the fit is flat to within the tie tolerance toward the floor. "
                 "The optimum is not clearly bracketed: read the figure, and consider a sweep that spans lower f "
                 "again.",
}


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
def write_csv_atomic(df, path):
    tmp = path.with_name(path.name + ".part")
    df.to_csv(tmp, index=False)
    os.replace(str(tmp), str(path))


def write_json_atomic(path, obj):
    tmp = path.with_name(path.name + ".part")
    tmp.write_text(json.dumps(obj, indent=2, default=str))
    os.replace(str(tmp), str(path))


def fnum(v, spec):
    try:
        if v is None or not np.isfinite(v):
            return "n/a"
    except TypeError:
        return "n/a"
    return format(v, spec)


def u_log(v, lo, hi):
    return (np.log10(np.asarray(v, dtype=float)) - np.log10(lo)) / (np.log10(hi) - np.log10(lo))


def u_lin(v, lo, hi):
    return (np.asarray(v, dtype=float) - lo) / (hi - lo)


def load_design(args, summary_dir):
    """The design file written by run_lowf_sweep_110.py, or the box from the command line."""
    ser = args.series
    p = summary_dir / ("lhs_design_lowf_%s_%s_110.json" % (args.label, ser))
    if p.exists():
        try:
            d = json.loads(p.read_text())
        except Exception as exc:
            sys.exit("Cannot read the design file %s: %s" % (p, exc))
        if d.get("series") != ser:
            sys.exit("The design file %s is for series %s, not %s." % (p.name, d.get("series"), ser))
        box = d["box"]
        return {"f_lo": float(box["f_RS_abs"]["lo"]), "f_hi": float(box["f_RS_abs"]["hi"]),
                "ks_lo": float(box["Ks_mult"]["lo"]), "ks_hi": float(box["Ks_mult"]["hi"]),
                "n": int(d["n"]), "seed": int(d["seed"]), "source": str(p), "design": d}
    need = [args.f_lo, args.f_hi, args.ks_lo, args.ks_hi]
    if any(v is None for v in need):
        sys.exit("No design file %s. Either run this on a sweep made by run_lowf_sweep_110.py, or give the box "
                 "with --f_lo --f_hi --ks_lo --ks_hi (and --expect_n)." % p.name)
    return {"f_lo": args.f_lo, "f_hi": args.f_hi, "ks_lo": args.ks_lo, "ks_hi": args.ks_hi,
            "n": args.expect_n, "seed": None, "source": "command line", "design": None}


def add_positions(m, box):
    m = m.copy()
    m["u_f"] = u_log(m["f_RS_abs"], box["f_lo"], box["f_hi"])
    m["u_Ks"] = u_lin(m["Ks_mult"], box["ks_lo"], box["ks_hi"])
    return m


def band_table(m, box, nb):
    rows = []
    ok = m[np.isfinite(m["kge_2012"])]
    for b in range(nb):
        lo_u, hi_u = b / float(nb), (b + 1) / float(nb)
        if b == nb - 1:
            sel = ok[(ok["u_f"] >= lo_u - 1e-12) & (ok["u_f"] <= hi_u + 1e-9)]
        else:
            sel = ok[(ok["u_f"] >= lo_u - 1e-12) & (ok["u_f"] < hi_u - 1e-12)]
        row = {"band": b + 1,
               "f_from": box["f_lo"] * (box["f_hi"] / box["f_lo"]) ** lo_u,
               "f_to": box["f_lo"] * (box["f_hi"] / box["f_lo"]) ** hi_u,
               "n": int(len(sel))}
        if len(sel):
            best = sel.loc[sel["kge_2012"].idxmax()]
            top3 = sel.nlargest(3, "kge_2012")
            row.update({"best_kge_2012": float(best["kge_2012"]), "median_kge_2012": float(sel["kge_2012"].median()),
                        "mean_top3_kge_2012": float(top3["kge_2012"].mean()), "f_of_best": float(best["f_RS_abs"]),
                        "Ks_of_best": float(best["Ks_mult"]), "pbias_of_best_pct": float(best["pbias_pct"]),
                        "median_Ks_of_top3": float(top3["Ks_mult"].median())})
        rows.append(row)
    return pd.DataFrame(rows)


def verdict(m, bands):
    """Apply the reading rule. Returns (code, facts dict)."""
    ok = m[np.isfinite(m["kge_2012"])].sort_values("kge_2012", ascending=False, kind="mergesort")
    best = ok.iloc[0]
    top5 = ok.head(5)
    n_near = int((top5["u_f"] < NEAR_FLOOR_U).sum())
    bb = bands.dropna(subset=["best_kge_2012"]).reset_index(drop=True)
    best_band = int(bb.loc[bb["best_kge_2012"].idxmax(), "band"])
    first_band = int(bb["band"].iloc[0])
    gap_low = float(bb["best_kge_2012"].max() - bb.loc[bb["band"] == first_band, "best_kge_2012"].iloc[0])
    facts = {"best_kge_2012": float(best["kge_2012"]), "best_f": float(best["f_RS_abs"]), "best_Ks": float(best["Ks_mult"]),
             "best_u_f": float(best["u_f"]), "best_u_Ks": float(best["u_Ks"]), "n_top5_near_floor": n_near,
             "best_band": best_band, "lowest_band": first_band, "lowest_band_gap_to_best_band": gap_low}
    if best["u_f"] > AT_CEIL_U:
        return "CEILING", facts
    if best["u_f"] < AT_FLOOR_U or n_near >= NEAR_FLOOR_COUNT:
        return "FLOOR", facts
    if best_band != first_band and gap_low >= TURN_TOL:
        return "BRACKETED", facts
    return "UNSETTLED", facts


def old_context(summary_dir, arm, series):
    """The Stage 0 runs of the same series and arm, or None."""
    p = summary_dir / OLD_DIRNAME / arm / "real_gauge_long_110.csv"
    if not p.exists():
        return None, "no %s (the old sweep is not drawn)" % p
    try:
        d = pd.read_csv(p)
    except Exception as exc:
        return None, "cannot read %s: %s" % (p, exc)
    if "series" not in d.columns:
        return None, "%s has no series column" % p
    d = d[d["series"] == SERIES_TEXT[series]].copy()
    d = d[np.isfinite(d["kge_2012"])]
    return (d if len(d) else None), str(p)


# ------------------------------------------------------------------
# Printing
# ------------------------------------------------------------------
def print_arm(arm, is_primary, m, bands, code, facts, box, n_scored, series, partial_note, old_df):
    print("\n" + "#" * 110)
    print("ARM: %s   %s   series %s%s" % (arm.upper(), "(PRIMARY: fixed before the run)" if is_primary else "(sensitivity arm)",
                                         series, partial_note))
    print("     %d real 5-minute points scored.  f box %g - %g (log), Ks box %g - %g." %
          (n_scored, box["f_lo"], box["f_hi"], box["ks_lo"], box["ks_hi"]))
    print("#" * 110)
    ok = m[np.isfinite(m["kge_2012"])]
    print("\nBest run:   KGE_2012 %.4f   Ks %.3f   f %.5f   PBIAS %+.1f%%   (position in the f box u = %.2f, in the Ks box %.2f)"
          % (facts["best_kge_2012"], facts["best_Ks"], facts["best_f"],
             float(ok.loc[ok["kge_2012"].idxmax(), "pbias_pct"]), facts["best_u_f"], facts["best_u_Ks"]))
    print("Median KGE_2012 over all %d runs: %.4f" % (len(ok), float(ok["kge_2012"].median())))
    top = ok.nlargest(10, "kge_2012")
    print("\nThe 10 best runs:")
    print("  %4s %9s %8s %10s %7s %8s %6s %6s %6s" % ("rank", "KGE_2012", "Ks", "f", "u_f", "PBIAS%", "r", "beta", "gamma"))
    for k, (_, r) in enumerate(top.iterrows()):
        print("  %4d %9.4f %8.3f %10.5f %7.2f %8.1f %6.3f %6.3f %6.3f" %
              (k + 1, r["kge_2012"], r["Ks_mult"], r["f_RS_abs"], r["u_f"], r["pbias_pct"], r["kge_r"], r["kge_beta"], r["kge_gamma"]))
    print("  (u_f: 0 = the f floor of this box, 1 = its ceiling.  r, beta, gamma are the three KGE components.)")

    print("\nBy band of f (equal widths in log f):")
    print("  %4s %21s %4s | %9s %9s | %10s %8s %8s" % ("band", "f from - to", "n", "best KGE", "median", "f of best", "Ks best", "PBIAS%"))
    for _, r in bands.iterrows():
        rng = "%.5f - %.5f" % (r["f_from"], r["f_to"])
        if r["n"] == 0:
            print("  %4d %21s %4d |   (no runs)" % (r["band"], rng, r["n"]))
        else:
            print("  %4d %21s %4d | %9.4f %9.4f | %10.5f %8.3f %8.1f" %
                  (r["band"], rng, r["n"], r["best_kge_2012"], r["median_kge_2012"], r["f_of_best"], r["Ks_of_best"], r["pbias_of_best_pct"]))
    top25 = ok.nlargest(25, "kge_2012")
    if len(top25) >= 8 and np.ptp(top25["f_RS_abs"]) > 0 and np.ptp(top25["Ks_mult"]) > 0:
        rho = stats.spearmanr(np.log10(top25["f_RS_abs"]), top25["Ks_mult"])
        sl = float(np.polyfit(np.log10(top25["f_RS_abs"]), top25["Ks_mult"], 1)[0])
        print("\nKs against f among the 25 best runs: rank correlation %+.2f (p %s), slope %+.2f Ks per tenfold change of f."
              % (rho[0], "<0.001" if rho[1] < 0.001 else "%.3f" % rho[1], sl))
        print("  A positive slope means the best Ks falls as f falls (the Ks-f ridge). A guide only: the best runs are not independent.")

    print("\n" + "=" * 110)
    print("READING (rule fixed before the run), %s ARM" % arm.upper())
    print("=" * 110)
    print("  best run u_f = %.2f;  of the 5 best runs, %d have u_f < %.2f;  best band = %d (lowest band = %d), the lowest band's best "
          "is %.4f below the best band's best (tie tolerance %.2f)"
          % (facts["best_u_f"], facts["n_top5_near_floor"], NEAR_FLOOR_U, facts["best_band"], facts["lowest_band"],
             facts["lowest_band_gap_to_best_band"], TURN_TOL))
    print("\n  VERDICT: %s" % code)
    print("  %s" % VERDICTS[code])
    if facts["best_u_Ks"] < 0.05 or facts["best_u_Ks"] > 0.95:
        print("\n  NOTE: the best run is also at the edge of the Ks box (position %.2f). The Ks range may be clipping the ridge too."
              % facts["best_u_Ks"])

    if old_df is not None:
        ob = old_df.loc[old_df["kge_2012"].idxmax()]
        nb = facts["best_kge_2012"]
        print("\nCompared with the OLD sweep (Stage 0, f 0.004-0.030, same event, same arm, same series):")
        print("  old best KGE_2012 %.4f (Ks %.3f, f %.5f);  new best %.4f;  new minus old %+.4f"
              % (ob["kge_2012"], ob["Ks_mult"], ob["f_RS_abs"], nb, nb - ob["kge_2012"]))
        zone = (box["f_lo"], box["f_hi"])
        lo_ov, hi_ov = max(0.004, zone[0]), min(0.030, zone[1])
        if lo_ov < hi_ov:
            o_in = old_df[(old_df["f_RS_abs"] >= lo_ov) & (old_df["f_RS_abs"] <= hi_ov)]
            n_in = ok[(ok["f_RS_abs"] >= lo_ov) & (ok["f_RS_abs"] <= hi_ov)]
            if len(o_in) and len(n_in):
                print("  Overlap of the two boxes, f %.4f - %.4f: old best %.4f (%d runs), new best %.4f (%d runs). Both are random "
                      "draws, so they should be close; a large gap points to the Ks ranges or the sampling, not to f."
                      % (lo_ov, hi_ov, o_in["kge_2012"].max(), len(o_in), n_in["kge_2012"].max(), len(n_in)))
    print("\n  Caveats: one event; routing pinned at the synthetic-truth values; the SMPHQ rainfall artifact is in every run; gauge "
          "error is not quantified; the filled arm assumes a silent logger means no change.")


# ------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------
def _f_axis(ax, box, old_df=None):
    """Log f axis covering the new box and, when it is drawn, the whole old sweep."""
    import matplotlib
    ax.set_xscale("log")
    lo = min(box["f_lo"], 0.004) * 0.85
    hi = max(box["f_hi"], 0.004) * 1.2
    if old_df is not None and len(old_df):
        lo = min(lo, float(old_df["f_RS_abs"].min()) * 0.9)
        hi = max(hi, float(old_df["f_RS_abs"].max()) * 1.12)
    ax.set_xlim(lo, hi)
    ticks = [t for t in (0.0005, 0.001, 0.002, 0.003, 0.004, 0.006, 0.01, 0.02, 0.03) if lo <= t <= hi]
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: "%g" % v))
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())


def fig_kge_vs_f(m, bands, old_df, box, label, series, arm, path, partial_note):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    ok = m[np.isfinite(m["kge_2012"])]
    ymax = float(ok["kge_2012"].max())
    if old_df is not None:
        ymax = max(ymax, float(old_df["kge_2012"].max()))
    ymin = ymax - 0.7
    pad = 0.04 * (ymax - ymin)
    fig, ax = plt.subplots(figsize=(9.2, 5.3), facecolor=RG.SURFACE)
    RG._style(ax)
    handles = []
    if old_df is not None:
        ax.scatter(old_df["f_RS_abs"], old_df["kge_2012"], s=20, marker="o", facecolors=RG.SURFACE, edgecolors=RG.NULLC,
                   linewidths=1.0, zorder=2)
        handles.append(Line2D([0], [0], marker="o", color=RG.NULLC, markerfacecolor=RG.SURFACE, markeredgewidth=1.0,
                              linestyle="none", markersize=5, label="old sweep, f 0.004-0.030 (%d runs)" % len(old_df)))
    ax.scatter(ok["f_RS_abs"], ok["kge_2012"], s=26, marker="o", color=RG.C_F, alpha=0.8, linewidths=0, zorder=3)
    handles.append(Line2D([0], [0], marker="o", color=RG.C_F, linestyle="none", markersize=5.5, alpha=0.8,
                          label="new sweep %s (%d runs)" % (label, len(ok))))
    top = ok.nlargest(10, "kge_2012")
    ax.scatter(top["f_RS_abs"], top["kge_2012"], s=75, marker="o", color=RG.C_F, edgecolors=RG.INK, linewidths=1.4, zorder=4)
    handles.append(Line2D([0], [0], marker="o", color=RG.C_F, markeredgecolor=RG.INK, markeredgewidth=1.4, linestyle="none",
                          markersize=8, label="the 10 best new runs"))
    bb = bands.dropna(subset=["best_kge_2012"])
    if len(bb) >= 2:
        xs = np.sqrt(bb["f_from"].to_numpy(float) * bb["f_to"].to_numpy(float))
        ax.plot(xs, bb["best_kge_2012"], color=RG.INK, linewidth=1.8, marker="D", markersize=5, zorder=5)
        handles.append(Line2D([0], [0], color=RG.INK, linewidth=1.8, marker="D", markersize=5, label="best of each f band (new sweep)"))
    ax.axvline(0.004, color=RG.INK2, linestyle=":", linewidth=1.2, zorder=1)
    handles.append(Line2D([0], [0], color=RG.INK2, linestyle=":", linewidth=1.2, label="f = 0.004, the floor of the old sweep"))
    _f_axis(ax, box, old_df)
    ax.set_ylim(ymin - pad, ymax + pad)
    hidden = int((ok["kge_2012"] < ymin).sum()) + (int((old_df["kge_2012"] < ymin).sum()) if old_df is not None else 0)
    ax.set_xlabel("f_RS_abs (1/mm), log scale", color=RG.INK2, fontsize=9)
    ax.set_ylabel("KGE_2012 against the real gauge (%s arm)" % arm, color=RG.INK2, fontsize=9)
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=8, labelcolor=RG.INK2)
    sub = "" if not hidden else "  (%d run(s) below the plotted range not shown)" % hidden
    ax.set_title("%s arm, series %s: fit to the real gauge against f%s%s" % (arm.upper(), series, partial_note, sub),
                 color=RG.INK, fontsize=10, loc="left")
    fig.tight_layout(rect=(0, 0.14, 1, 1))
    fig.savefig(path, dpi=150, facecolor=RG.SURFACE)
    plt.close(fig)


def fig_ks_vs_f(m, old_df, box, label, series, arm, path, partial_note):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    ok = m[np.isfinite(m["kge_2012"])]
    fig, ax = plt.subplots(figsize=(9.2, 5.3), facecolor=RG.SURFACE)
    RG._style(ax)
    handles = []
    ax.scatter(ok["f_RS_abs"], ok["Ks_mult"], s=14, color=RG.NULLC, alpha=0.35, linewidths=0, zorder=1)
    handles.append(Line2D([0], [0], marker="o", color=RG.NULLC, alpha=0.35, linestyle="none", markersize=5,
                          label="all new runs (where the sweep sampled)"))
    if old_df is not None:
        t_old = old_df.nlargest(25, "kge_2012")
        ax.scatter(t_old["f_RS_abs"], t_old["Ks_mult"], s=34, marker="o", facecolors=RG.SURFACE, edgecolors=RG.NULLC,
                   linewidths=1.2, zorder=2)
        handles.append(Line2D([0], [0], marker="o", color=RG.NULLC, markerfacecolor=RG.SURFACE, markeredgewidth=1.2,
                              linestyle="none", markersize=6, label="25 best runs of the old sweep"))
    t_new = ok.nlargest(25, "kge_2012")
    ax.scatter(t_new["f_RS_abs"], t_new["Ks_mult"], s=46, marker="o", color=RG.C_KS, alpha=0.9, linewidths=0, zorder=3)
    handles.append(Line2D([0], [0], marker="o", color=RG.C_KS, linestyle="none", markersize=6.5, label="25 best runs of the new sweep"))
    best = ok.loc[ok["kge_2012"].idxmax()]
    ax.scatter([best["f_RS_abs"]], [best["Ks_mult"]], s=110, marker="o", color=RG.C_KS, edgecolors=RG.INK, linewidths=1.6, zorder=4)
    handles.append(Line2D([0], [0], marker="o", color=RG.C_KS, markeredgecolor=RG.INK, markeredgewidth=1.6, linestyle="none",
                          markersize=9, label="the best new run"))
    ax.axvline(0.004, color=RG.INK2, linestyle=":", linewidth=1.2, zorder=1)
    handles.append(Line2D([0], [0], color=RG.INK2, linestyle=":", linewidth=1.2, label="f = 0.004, the floor of the old sweep"))
    _f_axis(ax, box, old_df)
    ax.set_ylim(min(box["ks_lo"], 4.0) - 0.3, max(box["ks_hi"], 10.5) + 0.3)
    ax.set_xlabel("f_RS_abs (1/mm), log scale", color=RG.INK2, fontsize=9)
    ax.set_ylabel("Ks_mult", color=RG.INK2, fontsize=9)
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=8, labelcolor=RG.INK2)
    ax.set_title("%s arm, series %s: where the best runs put Ks and f%s" % (arm.upper(), series, partial_note),
                 color=RG.INK, fontsize=10, loc="left")
    fig.tight_layout(rect=(0, 0.14, 1, 1))
    fig.savefig(path, dpi=150, facecolor=RG.SURFACE)
    plt.close(fig)


def fig_hydrograph(real_scored, record_pts, sim_best, best_row, sim_old, old_row, arm, series, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    fig, ax = plt.subplots(figsize=(9.4, 4.8), facecolor=RG.SURFACE)
    RG._style(ax)
    if sim_old is not None:
        ax.plot(sim_old.index, sim_old.values, color=RG.NULLC, linewidth=1.7, linestyle="--",
                label="best run of the old sweep (KGE_2012 %.3f, f %.4f)" % (old_row["kge_2012"], old_row["f_RS_abs"]))
    ax.plot(sim_best.index, sim_best.values, color=RG.C_F, linewidth=1.9,
            label="best new run (KGE_2012 %.3f, f %.5f, Ks %.2f)" % (best_row["kge_2012"], best_row["f_RS_abs"], best_row["Ks_mult"]))
    if arm == "asis":
        ax.plot(real_scored.index, real_scored.values, color=RG.INK, linewidth=0, marker="o", markersize=3.5,
                label="real gauge, as-is: the %d bins the scorer keeps" % len(real_scored))
    else:
        ax.plot(real_scored.index, real_scored.values, color=RG.INK, linewidth=1.5, label="real gauge, filled (%d points)" % len(real_scored))
        if record_pts is not None and len(record_pts):
            ax.plot(record_pts.index, record_pts.values, color=RG.INK, linewidth=0, marker="o", markersize=3.0,
                    label="bins that hold a gauge record (%d)" % len(record_pts))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %H:%M"))
    ax.set_xlabel("time (day of month and hour, August 2014)", color=RG.INK2, fontsize=9)
    ax.set_ylabel("discharge at the outlet (m$^3$/s)", color=RG.INK2, fontsize=9)
    ax.legend(loc="upper right", frameon=False, fontsize=8.5, labelcolor=RG.INK2)
    ax.set_title("%s arm, series %s: real gauge against the best new run (one event)" % (arm.upper(), series),
                 color=RG.INK, fontsize=10, loc="left")
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor=RG.SURFACE)
    plt.close(fig)


# ------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Series 110: score the runs of a low-f sweep against the REAL gauge, one series, and say whether the "
                    "best fit is inside the f range or at an edge.")
    ap.add_argument("--label", required=True, help="the --label of the sweep (for example lf1)")
    ap.add_argument("--series", required=True, choices=["OFF", "ON"], help="which series of that sweep")
    ap.add_argument("--primary_arm", choices=["asis", "filled"], default="filled",
                    help="the arm whose reading is the verdict (default filled, as decided 2026-10-04)")
    ap.add_argument("--gauge_mode", choices=["both", "asis", "filled"], default="both")
    ap.add_argument("--max_gap_min", type=float, default=RG.DEFAULT_MAX_GAP_MIN)
    ap.add_argument("--max_hold_hr", type=float, default=RG.DEFAULT_MAX_HOLD_HR)
    ap.add_argument("--bands", type=int, default=6, help="number of f bands, equal widths in log f (default 6)")
    ap.add_argument("--check_only", action="store_true", help="run the gates and stop; writes nothing")
    ap.add_argument("--min_points", type=int, default=30)
    ap.add_argument("--drop_bad_runs", action="store_true", help="exclude runs that fail gate 2A, 2B or 4 and continue")
    ap.add_argument("--no_plots", action="store_true")
    ap.add_argument("--results_csv", type=Path, default=None)
    ap.add_argument("--gauge_xlsx", type=Path, default=None)
    ap.add_argument("--out_dir", type=Path, default=None)
    ap.add_argument("--f_lo", type=float, default=None, help="only if there is no design file")
    ap.add_argument("--f_hi", type=float, default=None, help="only if there is no design file")
    ap.add_argument("--ks_lo", type=float, default=None, help="only if there is no design file")
    ap.add_argument("--ks_hi", type=float, default=None, help="only if there is no design file")
    ap.add_argument("--expect_n", type=int, default=None, help="only if there is no design file")
    args = ap.parse_args()
    if args.bands < 3 or args.bands > 20:
        ap.error("--bands must be between 3 and 20.")
    if args.max_gap_min <= 0 or args.max_hold_hr <= 0:
        ap.error("--max_gap_min and --max_hold_hr must be positive.")

    script_dir = Path.cwd()
    project_root = script_dir.parent
    calib_dir = project_root / "calibration_work"
    summary_dir = calib_dir / "03_comparisons" / "summary_tables"
    csv_dir = calib_dir / "03_comparisons" / "csv_exports"
    synth_dir = calib_dir / "synth_truth"
    series = args.series
    res_path = args.results_csv or (summary_dir / ("lhs_results_lowf_%s_%s_110.csv" % (args.label, series)))
    xlsx_path = args.gauge_xlsx or project_root.joinpath(*RG.GAUGE_REL)
    out_dir = args.out_dir or (summary_dir / OUT_DIRNAME / ("%s_%s" % (args.label, series)))

    arms_used = ["filled", "asis"] if args.gauge_mode == "both" else [args.gauge_mode]
    primary = args.primary_arm if args.gauge_mode == "both" else args.gauge_mode
    arms_used = [primary] + [a for a in arms_used if a != primary]

    print("\n" + "=" * 78)
    print("Series 110 -- score a LOW-f sweep against the REAL gauge (version 1)")
    print("label %s, series %s, python %s, pandas %s, numpy %s" % (args.label, series, sys.version.split()[0], pd.__version__, np.__version__))
    print("arms: %s   (primary: %s)" % (", ".join(arms_used), primary))
    print("=" * 78)
    if ".".join(pd.__version__.split(".")[:2]) != RG.EXPECT_PANDAS_MINOR:
        print("  *** NOTE: pandas here is %s; the runs were scored with pandas %s.x. Gate 2A will tell. ***"
              % (pd.__version__, RG.EXPECT_PANDAS_MINOR))

    # ---- inputs / safety -----------------------------------------------------------------------
    qouts = list(synth_dir.glob("*.qout")) if synth_dir.exists() else []
    if len(qouts) != 1:
        sys.exit("Expected exactly one *.qout in %s (the synthetic truth the runs were scored against), found %d: %s"
                 % (synth_dir, len(qouts), [q.name for q in qouts]))
    ccoff_path = qouts[0]
    print("synthetic truth (read-only, used only for gate 2): %s" % ccoff_path.name)
    if not xlsx_path.exists():
        sys.exit("The real gauge workbook was not found: %s (--gauge_xlsx points elsewhere)." % xlsx_path)
    out_res = out_dir.resolve()
    for guarded in (synth_dir, csv_dir):
        gd = guarded.resolve()
        if out_res == gd or gd in out_res.parents or out_res in gd.parents:
            sys.exit("Output folder %s overlaps an input folder (%s); refusing." % (out_dir, guarded))

    box = load_design(args, summary_dir)
    print("design: %s" % box["source"])
    print("box: f %g - %g (log), Ks %g - %g" % (box["f_lo"], box["f_hi"], box["ks_lo"], box["ks_hi"]))

    print("\nLoading the results table:")
    if not res_path.exists():
        sys.exit("Results file not found: %s. Run run_lowf_sweep_110.py --mode %s --label %s first." % (res_path, series.lower(), args.label))
    n_expect = box["n"] if box["n"] else len(pd.read_csv(res_path))
    res = RS.load_results(res_path, "%s %s" % (args.label, series), PERC[series], n_expect)
    res = res.reset_index(drop=True)
    n_design = box["n"]
    partial = bool(n_design) and len(res) < n_design
    partial_note = "   [PARTIAL: %d of %d runs]" % (len(res), n_design) if partial else ""
    if partial:
        print("  NOTE: %d of the %d designed runs are in the table (a partial sweep). The reading is labelled PARTIAL." % (len(res), n_design))

    # runs must be inside the box that the design file says
    tol = 1e-9
    out_f = res[(res["f_RS_abs"] < box["f_lo"] * (1 - tol)) | (res["f_RS_abs"] > box["f_hi"] * (1 + tol))]
    out_k = res[(res["Ks_mult"] < box["ks_lo"] * (1 - tol)) | (res["Ks_mult"] > box["ks_hi"] * (1 + tol))]
    ck = RG.Checks()
    if len(out_f) or len(out_k):
        ck.add("FAIL", "runs inside the box", "%d run(s) lie outside f %g-%g and %d outside Ks %g-%g. This is not the sweep "
               "the design file describes; the verdict would be wrong." % (len(out_f), box["f_lo"], box["f_hi"], len(out_k),
                                                                         box["ks_lo"], box["ks_hi"]))
    else:
        ck.add("PASS", "runs inside the box", "all %d runs lie inside f %g-%g and Ks %g-%g" % (len(res), box["f_lo"], box["f_hi"],
                                                                                              box["ks_lo"], box["ks_hi"]))
    if res.empty:
        sys.exit("The results table has no runs.")

    phase_fn = RS.import_phase_fn(script_dir)
    if phase_fn is None:
        sys.exit("The scorer's phase-metric function could not be imported (run_sensitivity_single_interp_tribs6.py must be in "
                 "this folder). Without it gate 2 cannot pass.")
    obs_synth = RS.read_truth_5min(ccoff_path)

    print("\nLoading compare CSVs (run membership = the run_id rows of the results table):")
    cache, bad = RS.load_all_compares(res, csv_dir)
    print("  %d/%d loaded" % (len(cache), len(res)))

    # ---- GATE 2 -------------------------------------------------------------------------------------
    allres = res.assign(_series=SERIES_TEXT[series])
    badA, worstA = RS.gate_a(allres, cache, obs_synth)
    badB, worstB = RS.gate_b(allres, cache, phase_fn)
    for dct in (badA, badB):
        for rid, why in dct.items():
            bad.setdefault(rid, why)
    n_checked = len(cache)
    ck.add("PASS" if not badA else "FAIL", "GATE 2A synthetic Observed",
           "%d/%d compare CSVs match the synthetic truth as read here (worst abs diff %.2e m3/s)" % (n_checked - len(badA), n_checked, worstA))
    if worstB:
        wk = max(worstB, key=worstB.get)
        ck.add("PASS" if not badB else "FAIL", "GATE 2B stored metrics",
               "%d/%d runs reproduce their stored metrics (worst abs deviation %.2e in %s)" % (n_checked - len(badB), n_checked, worstB[wk], wk))
    else:
        ck.add("PASS" if not badB else "FAIL", "GATE 2B stored metrics", "%d/%d" % (n_checked - len(badB), n_checked))
    if n_checked and len(badA) >= 0.9 * n_checked:
        print("\n  GATE 2A failed for nearly every run. The usual cause is a different pandas version from the one that scored the "
              "runs (here: pandas %s). Use the same python as for the runs." % pd.__version__)
    missing = {r: w for r, w in bad.items() if "missing" in w or "unreadable" in w}
    if missing:
        ck.add("FAIL", "compare CSVs", "%d compare CSV(s) missing or unreadable, e.g. %s" % (len(missing), ", ".join(list(missing)[:3])))

    # ---- GATE 1 ---------------------------------------------------------------------------------------
    print("\nReading the real gauge workbook (the whole 1993-2025 record is read, as the scorer does; this can take a minute) ...", flush=True)
    try:
        obs_df = RG.read_gauge_workbook(xlsx_path)
    except Exception as exc:
        sys.exit("Could not read the real gauge workbook the way the scorer does (%s: %s)." % (type(exc).__name__, str(exc)[:120]))
    real_5min = RG.gauge_to_grid(obs_df)
    peaks = [float(np.nanmax(d["Simulated"].to_numpy(dtype=float))) for d in cache.values()
             if len(d) and np.isfinite(d["Simulated"].to_numpy(dtype=float)).any()]
    sim_pk = float(np.median(peaks)) if peaks else None
    ck1, info1, real_valid, real_grid = RG.check_real_series(obs_df, real_5min, xlsx_path, args.min_points, sim_pk)
    if "asis" not in arms_used:
        for r in ck1.rows:
            if r["status"] == "FAIL" and r["check"] == "as-is window coverage":
                r["status"] = "WARN"
                r["detail"] += "  [as-is arm not in use]"
    filled, ftype, finfo = RG.build_filled(obs_df, real_5min, args.max_gap_min, args.max_hold_hr)
    t0w, t1w = pd.Timestamp(RG.EVENT_START), pd.Timestamp(RG.EVENT_END)
    raw_w = obs_df.loc[np.asarray((obs_df.index >= t0w) & (obs_df.index <= t1w))]
    filled_valid = filled.dropna()
    ck1f, infof = RG.check_filled(filled, ftype, finfo, real_valid, raw_w, args.min_points)
    if "filled" not in arms_used:
        for r in ck1f.rows:
            if r["status"] == "FAIL":
                r["status"] = "WARN"
                r["detail"] += "  [filled arm not in use]"
    ck.rows = ck1.rows + ck1f.rows + ck.rows
    gate3_state = ("NOT REPEATED here: PASSED on 2026-10-04 (Stage 0) for this pipeline; the scorer, builder and gauge reading are "
                   "unchanged")
    if ck1.n("FAIL") or ck1f.n("FAIL"):
        RG.print_gates(ck.rows, gate3_state)
        sys.exit("\nSTOPPED at gate 1: the real series is not usable (see FAIL above). Nothing was scored or written.")

    run_ids = [r for r in res["run_id"] if r in cache]
    if bad and not args.drop_bad_runs:
        RG.print_gates(ck.rows, gate3_state)
        print("\n  %d run(s) failed a gate or could not be loaded. First few:" % len(bad))
        for rid, why in list(bad.items())[:8]:
            print("    %s: %s" % (rid, why))
        sys.exit("\nSTOPPED: nothing was scored. Fix the cause (or re-run with --drop_bad_runs to exclude these runs).")
    if bad:
        for r in ck.rows:
            if r["status"] == "FAIL" and r["check"] in RG.DROPPABLE:
                r["status"] = "WARN"
                r["detail"] += "  [runs excluded by --drop_bad_runs]"
    run_ids = [r for r in run_ids if r not in bad]
    if not run_ids:
        sys.exit("No runs left after the gates; stopping.")

    # ---- GATE 4 ----------------------------------------------------------------------------------------
    valids = {"asis": real_valid, "filled": filled_valid}
    union_valid = pd.concat([valids[a] for a in arms_used]).groupby(level=0).first().sort_index()
    ck4, ref_idx, scored_union, bad4 = RG.check_grid(cache, run_ids, union_valid)
    if bad4:
        if not args.drop_bad_runs:
            ck.extend(ck4)
            RG.print_gates(ck.rows, gate3_state)
            sys.exit("\nSTOPPED at gate 4: stored grids differ between runs. Nothing was scored.")
        for rid, why in bad4.items():
            bad.setdefault(rid, why)
        run_ids = [r for r in run_ids if r not in bad4]
        ck.add("WARN", "GATE 4 runs excluded", "%d run(s) with a different stored grid were excluded" % len(bad4))
        ck4, ref_idx, scored_union, bad4 = RG.check_grid(cache, run_ids, union_valid)
    ck.extend(ck4)
    RG.print_gates(ck.rows, gate3_state)
    if ck.n("FAIL"):
        sys.exit("\nSTOPPED: a gate failed (see FAIL above). Nothing was written.")
    res = res[res["run_id"].isin(run_ids)].reset_index(drop=True)
    if res.empty:
        sys.exit("No runs left after the gates; stopping.")

    if args.check_only:
        print("\n--check_only: gates finished, nothing written, no fit shown. A full run scores %d runs: %s."
              % (len(res), "; ".join("%s arm on %d real points" % (a, len(valids[a].index.intersection(ref_idx))) for a in arms_used)))
        return

    # ---- score and analyse every arm in use -----------------------------------------------------------
    results = {}
    for arm in arms_used:
        valid = valids[arm]
        scored_idx = valid.index.intersection(ref_idx)
        real_scored = valid.loc[scored_idx]
        cache_g = {rid: pd.DataFrame({"Observed": real_scored, "Simulated": cache[rid]["Simulated"].reindex(scored_idx)})
                   for rid in res["run_id"]}
        m, failed = TL.score_truth(real_scored, res, cache_g, phase_fn)
        if failed:
            print("  WARNING (%s arm): %d run(s) could not be scored" % (arm, len(failed)))
        if m.empty or not np.isfinite(m["kge_2012"]).any():
            sys.exit("Nothing could be scored in the %s arm." % arm)
        m = add_positions(m, box)
        m["band"] = np.minimum((m["u_f"] * args.bands).astype(int), args.bands - 1) + 1
        bands = band_table(m, box, args.bands)
        code, facts = verdict(m, bands)
        old_df, old_note = old_context(summary_dir, arm, series)
        results[arm] = {"m": m, "bands": bands, "code": code, "facts": facts, "scored_idx": scored_idx,
                        "real_scored": real_scored, "old_df": old_df, "old_note": old_note}
        print_arm(arm, arm == primary, m, bands, code, facts, box, len(scored_idx), series, partial_note, old_df)
        if old_df is None:
            print("\n  (Old sweep not drawn: %s)" % old_note)

    # ---- agreement between the arms ------------------------------------------------------------------------
    if len(arms_used) == 2:
        a, b = arms_used
        ca, cb = results[a]["code"], results[b]["code"]
        print("\n" + "=" * 110)
        print("ARM AGREEMENT: %s (primary) reads %s; %s reads %s." % (a, ca, b, cb))
        print("=" * 110)
        if ca != cb:
            print("  The arms disagree. The primary arm's verdict is the one fixed beforehand; report the disagreement with it.")

    # ---- write ------------------------------------------------------------------------------------------------------
    out_dir.mkdir(parents=True, exist_ok=True)
    write_csv_atomic(pd.DataFrame(ck.rows), out_dir / "lowf_gates_110.csv")
    record_pts = valids["asis"].loc[valids["asis"].index.intersection(ref_idx)]
    for arm in arms_used:
        adir = out_dir / arm
        adir.mkdir(parents=True, exist_ok=True)
        r = results[arm]
        m = r["m"]
        write_csv_atomic(m, adir / "lowf_long_110.csv")
        keep = ["rank", "run_id", "Ks_mult", "f_RS_abs", "channelconductivity_mmhr", "u_f", "u_Ks", "kge_2012", "kge_r", "kge_beta",
                "kge_gamma", "pbias_pct", "nse", "rmse_m3s", "peak_error_pct", "peak_timing_error_hr"]
        top = m[np.isfinite(m["kge_2012"])].sort_values("kge_2012", ascending=False, kind="mergesort").head(25).reset_index(drop=True)
        top["rank"] = np.arange(1, len(top) + 1)
        write_csv_atomic(top[keep], adir / "lowf_top_runs_110.csv")
        write_csv_atomic(r["bands"], adir / "lowf_f_bands_110.csv")
        vrow = dict(r["facts"])
        vrow.update({"arm": arm, "primary": arm == primary, "verdict": r["code"], "label": args.label, "series": series,
                     "n_runs": int(len(m)), "n_design": n_design, "partial": partial, "f_lo": box["f_lo"], "f_hi": box["f_hi"],
                     "ks_lo": box["ks_lo"], "ks_hi": box["ks_hi"]})
        write_csv_atomic(pd.DataFrame([vrow]), adir / "lowf_verdict_110.csv")
        if args.no_plots:
            continue
        figs = [("fig_lowf_kge_vs_f_110.png", lambda p: fig_kge_vs_f(m, r["bands"], r["old_df"], box, args.label, series, arm, p, partial_note)),
                ("fig_lowf_Ks_vs_f_110.png", lambda p: fig_ks_vs_f(m, r["old_df"], box, args.label, series, arm, p, partial_note))]

        def _hydro(p, arm=arm, r=r, m=m):
            best = m.loc[m["kge_2012"].idxmax()]
            sim_best = cache[best["run_id"]]["Simulated"].reindex(ref_idx)
            sim_old, old_row = None, None
            if r["old_df"] is not None:
                old_row = r["old_df"].loc[r["old_df"]["kge_2012"].idxmax()]
                pth = csv_dir / ("%s_compare_obs_sim.csv" % old_row["run_id"])
                if pth.exists():
                    try:
                        sim_old = RS.load_compare(pth)["Simulated"].reindex(ref_idx)
                    except Exception:
                        sim_old = None
            fig_hydrograph(r["real_scored"], record_pts if arm == "filled" else None, sim_best, best, sim_old, old_row, arm, series, p)
        figs.append(("fig_lowf_hydrograph_110.png", _hydro))
        for name, fn in figs:
            try:
                fn(adir / name)
            except Exception as exc:
                print("\n  (figure %s/%s skipped: %s -- the CSVs are saved regardless)" % (arm, name, exc))

    scorer_py = script_dir / "run_sensitivity_single_interp_tribs6.py"
    prov = {
        "script": Path(__file__).name, "version": 1, "created_local": datetime.now().isoformat(timespec="seconds"),
        "python": sys.version.split()[0], "pandas": pd.__version__, "numpy": np.__version__,
        "label": args.label, "series": series, "arms": {"used": arms_used, "primary": primary,
                                                          "max_gap_min": args.max_gap_min, "max_hold_hr": args.max_hold_hr,
                                                          "filled_counts": finfo, "filled_checks": infof},
        "box": box, "n_runs_scored": int(len(res)), "partial": partial,
        "reading_rule": {"AT_CEIL_U": AT_CEIL_U, "AT_FLOOR_U": AT_FLOOR_U, "NEAR_FLOOR_U": NEAR_FLOOR_U,
                         "NEAR_FLOOR_COUNT": NEAR_FLOOR_COUNT, "TURN_TOL": TURN_TOL, "bands": args.bands},
        "verdicts": {a: results[a]["code"] for a in arms_used},
        "inputs": {"results_csv": {"path": str(res_path), "md5": RS.md5_of(res_path)},
                   "gauge_xlsx": {"path": str(xlsx_path), "md5": info1.get("xlsx_md5")},
                   "synthetic_truth": {"file": ccoff_path.name, "md5": RS.md5_of(ccoff_path)},
                   "scorer_py_md5": RS.md5_of(scorer_py) if scorer_py.exists() else None,
                   "rescore_real_gauge_md5": RS.md5_of(Path(RG.__file__)),
                   "old_sweep": {a: results[a]["old_note"] for a in arms_used}},
        "gates": {"2A_max_abs_diff_m3s": worstA, "2B_worst_abs_dev": worstB, "runs_checked": n_checked,
                  "runs_excluded": {k: v for k, v in bad.items()}, "3": gate3_state, "checks": ck.rows},
        "settings": {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()},
    }
    write_json_atomic(out_dir / "PROVENANCE_lowf_110.json", prov)
    print("\nSaved to: %s" % out_dir)
    print("  top folder: lowf_gates_110.csv   PROVENANCE_lowf_110.json")
    print("  per arm (%s): lowf_long_110.csv   lowf_top_runs_110.csv   lowf_f_bands_110.csv   lowf_verdict_110.csv   fig_lowf_*.png"
          % ", ".join("%s/" % a for a in arms_used))


# the docstring carries the rule constants; fill them in
if __doc__:
    __doc__ = __doc__ % {"ceil": AT_CEIL_U, "floor": AT_FLOOR_U, "nn": NEAR_FLOOR_COUNT, "nf": NEAR_FLOOR_U, "tol": TURN_TOL}

if __name__ == "__main__":
    main()
