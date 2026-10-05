"""
score_stageb_110.py   (version 1)
=================================
Series 110, STAGE B -- score the two sweeps of Stage B against the REAL SMF gauge and answer one question:

    With f held at one value (0.0007) and Ks free in BOTH series, does turning channel loss ON (cc)
    make the best fit to the real gauge better than leaving it OFF, by more than your margin (0.05 KGE_2012)?

The two sweeps come from run_lowf_sweep_110.py (version 1.2), same --label, f pinned with --f_fixed:
    OFF   --mode off   Ks 5.5 - 8.5    30 runs   (cc is written into the input file but is inert)
    ON    --mode on    Ks 3.5 - 8.5   120 runs   (cc 30 - 1000 mm/hr, log scale)
They are deliberately NOT the same points: OFF is a one-dimensional search in Ks, ON a two-dimensional search
in Ks and cc. Each series is allowed to find its own best Ks, so the comparison is best run against best run.
(Comparing the two at the same Ks would be unfair to ON: channel loss removes water, so ON needs a lower Ks
to reach the same volume. In the old 500-run sweep, going from cc 30 to cc 1000 moved the volume-matching Ks
from about 7.8 down to about 4.3.)

No tRIBS is started. This script only reads the stored compare CSVs, like score_lowf_real_110.py, and REUSES
the reading, filling, gate and metric code of rescore_real_gauge_110.py (version 2) unchanged: keep it in this
folder together with rescore_cc_truths_110.py, rescore_truth_location_110.py and
run_sensitivity_single_interp_tribs6.py.

THE TWO ARMS (same as Stage 0 and the low-f check)
---------------------------------------------------
  filled  all 241 five-minute bins; gaps in the event-triggered gauge are filled from the record (line up
          to 60 min, then hold the last reading up to 24 h).  PRIMARY (your decision of 2026-10-04).
  as-is   only the bins that hold a gauge record (23 of 241), what the scorer does in real-gauge mode.
          A sensitivity check.

THE READING RULE (fixed here, BEFORE any result of Stage B has been seen)
--------------------------------------------------------------------------
Your margin of %(margin).2f KGE_2012 (--margin) and the five-row table of the Stage 0 handoff, unchanged:
  gap       = best ON run minus best OFF run (each series free in Ks), PRIMARY arm.
  location  = where the best ON runs sit in the log cc range (the best %(topfrac)d%% of the ON runs; chance puts a
              third in each third). "Low" or "high" = that third holds more than chance at p < %(alpha).2f;
              "unstable" = the best runs span more than %(span).2f of the log cc range.
  gap below -margin                                          -> row 4  ON worse
  gap within the margin, best runs crowd the lowest third    -> row 1  ON about equal, cc at the floor
  gap within the margin, best runs do NOT crowd it           -> row 5  ON about equal, cc not at the floor
  FLAT override: if the best ON run of every cc band is within %(flat).2f of the others, row 1 is read as row 5
    (a fit that is the same at every cc cannot be evidence that the data prefer no loss; crowding is then chance)
  gap above the margin, best runs in an edge third or unstable -> row 3  ON better but not a clean signal
  gap above the margin, best runs interior and stable        -> row 2  ON better, a real signal to test further
CHECKS that must agree (stability): (i) the same rule applied to the MEAN of the best %(topfrac)d%% of runs in each
series (ON 12, OFF 3 when the sweeps are complete); (ii) ONLY when the gap is within the margin: the same rule with
the location read from the best %(topfrac2)d%% of the ON runs (24) instead of the best %(topfrac)d%%. The ON sweep puts
only a handful of runs close to the Ks-cc ridge, so which third of the cc range holds them can be chance. If a check
points to a different row than the best-run reading, the answer is reported as UNSETTLED. (The wider location check
is not used when the gap is above the margin: counting more runs then only dilutes a small region of good runs.)
NEAR-MARGIN flag (mine): the ON sweep covers a plane with 120 points while OFF covers a line with 30, so ON's best
point is a slight UNDER-estimate of what ON can reach. If the gap lands within %(nm).2f BELOW the margin, the
script says so, and a refinement of the ON sweep around its best point would be the way to settle it.
Everything is read from the primary (filled) arm; the as-is arm gets the same reading as a check, and any
disagreement between the arms is printed. A partial sweep can be scored but is labelled PARTIAL.

CHECKS ("gates", run before anything is scored; same code as Stage 0)
----------------------------------------------------------------------
  GATE 1   the real series (workbook read, duplicates, values, peak scale, cadence) and 1b the filled series.
  GATE 2A  every run's stored Observed column equals the synthetic truth in synth_truth/ as read here.
  GATE 2B  every run's stored metrics are reproduced by the metric code here.
  GATE 4   all runs share one regular 5-minute grid and every scored real bin is on it.
  Also: f was held at the same value in both design files and in every run; the two series share one label and
  routing; every run lies inside the box of its own design file; the ON cc box is the standard 30 - 1000 box;
  no run id appears in both series.
  GATE 3 (the actual scorer in real-gauge mode on one point) is NOT repeated: it PASSED on 2026-10-04 for this
  pipeline, and these runs use the same scorer, builder and gauge reading, unchanged.

USAGE (from lab/)
------------------
    python score_stageb_110.py --label sb1 --check_only      # gates only; writes nothing
    python score_stageb_110.py --label sb1                    # the full scoring
  It only READS CSVs, so it is safe to run while the ON sweep is still going; expect PARTIAL then.

WHAT YOU GET   (calibration_work/03_comparisons/summary_tables/stageb_real_gauge_110/<label>/)
  stageb_gates_110.csv   PROVENANCE_stageb_110.json
  filled/ and asis/ :  stageb_long_110.csv (every run of both series, all metrics)
                       stageb_summary_110.csv   stageb_top_runs_110.csv   stageb_cc_bands_110.csv
                       stageb_cc_location_110.csv   stageb_ridge_110.csv   stageb_verdict_110.csv
                       fig_stageb_kge_vs_cc_110.png   fig_stageb_kge_vs_Ks_110.png
                       fig_stageb_ridge_110.png       fig_stageb_hydrograph_110.png
                       fig_stageb_hydrograph_zoom_110.png
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
              "check_grid", "Checks", "print_gates", "_style", "_cc_axis", "cc_box_logs", "cc_unit", "cc_location",
              "decision_row", "series_summary", "top_rows", "score_arm", "fig_hydrographs", "print_cc_location",
              "print_arm_agreement", "fcc", "fmt_p", "ALPHA", "STABLE_SPAN", "FLOOR_TOL_MMHR", "C_CC", "C_KS", "NULLC",
              "INK", "INK2", "SURFACE", "GRID"):
    if not hasattr(RG, _name):
        sys.exit("rescore_real_gauge_110.py in this folder is not version 2 (missing %s). Replace it with the "
                 "current copy." % _name)

# ------------------------------------------------------------------
# Constants. The reading rule is fixed before any result is seen; the numbers are rules of thumb of mine.
# ------------------------------------------------------------------
DEFAULT_MARGIN = 0.05       # your margin (KGE_2012), the same as Stage 0
TOP_FRAC = 0.10             # the "best runs" of a series: its best 10 per cent (at least TOP_MIN runs)
TOP_MIN = 3
NEAR_MARGIN_TOL = 0.02      # a gap this far below the margin is flagged: ON's best is a slight under-estimate
LOW_CC_FRAC = 0.10          # "lowest cc" for the consistency check: the lowest tenth of the log cc range
LOW_CC_NOTE_TOL = 0.03      # ON at its lowest cc should reproduce OFF at the same Ks, to about this much KGE
LOW_CC_WARN_TOL = 0.10
FLAT_BANDS_TOL = 0.02      # the best ON run of every cc band within this much KGE of the others: cc is not pinned
RIDGE_FRAC = 0.20           # the Ks-cc ridge is read from the best 20 per cent of the ON runs (at least RIDGE_MIN)
RIDGE_MIN = 10
RIDGE_NEAR = 0.02           # figure: ON runs within this much KGE of the best ON run are drawn as "near-best"
VOL_R2_MIN = 0.80           # the straight-line volume fit is drawn / quoted only if R2 is at least this
OLD_PBIAS_PER_DECADE = -26.9    # Stage 0 (500 runs, f 0.004 - 0.03): ON minus OFF PBIAS per tenfold rise of cc
OUT_DIRNAME = "stageb_real_gauge_110"
SER_ON, SER_OFF = RS.SER_ON, RS.SER_OFF

READINGS = {
    1: "ON about equal to OFF, and the best ON runs crowd the low end of cc (where channel loss removes the least "
       "water): this event does not need channel loss to fit the real gauge. cc stays useful as a streambed-"
       "infiltration sensitivity tool, not as a calibration improvement.",
    2: "ON better than OFF by more than the margin, and the best ON runs agree on one interior cc: evidence that the "
       "data want channel loss. Next: check that it holds on another storm and with free routing before claiming it.",
    3: "ON better than OFF by more than the margin, but the best ON runs sit in an edge third of the cc range or are "
       "spread across it: extra flexibility, not a clear process signal. Do not claim it; check another storm first.",
    4: "ON worse than OFF by more than the margin. With Ks free in both series this should not happen unless the ON "
       "sweep missed its own best region (120 points in two dimensions against 30 in one). Read the cc-band table "
       "and the figures before concluding anything.",
    5: "ON about equal to OFF, but the best ON runs are NOT crowded at low cc: the real gauge cannot tell whether "
       "channel loss is present. Usually this is because a lower Ks makes up the water that cc removes (a Ks-cc ridge); "
       "it can also mean cc has no effect on the fit. Report cc as not identifiable from this event; read the ridge "
       "figure and the trade-off section.",
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


def u_lin(v, lo, hi):
    return (np.asarray(v, dtype=float) - lo) / (hi - lo)


def top_count(n, frac):
    return int(min(max(TOP_MIN, int(round(frac * n))), max(n, 1)))


def read_design(summary_dir, label, series):
    p = summary_dir / ("lhs_design_lowf_%s_%s_110.json" % (label, series))
    if not p.exists():
        sys.exit("No design file %s. Run run_lowf_sweep_110.py (version 1.2) with --mode %s --label %s --f_fixed 0.0007 "
                 "first (see the top of this file)." % (p.name, series.lower(), label))
    try:
        d = json.loads(p.read_text())
    except Exception as exc:
        sys.exit("Cannot read the design file %s: %s" % (p, exc))
    if d.get("series") != series:
        sys.exit("The design file %s is for series %s, not %s." % (p.name, d.get("series"), series))
    if d.get("label") != label:
        sys.exit("The design file %s has label %s, not %s." % (p.name, d.get("label"), label))
    if d.get("f_fixed") is None:
        sys.exit("The design file %s has no pinned f (f_fixed). This is a sweep that SAMPLES f, not a Stage B sweep. "
                 "Run run_lowf_sweep_110.py (version 1.2) with --f_fixed, or use score_lowf_real_110.py for a sampled-f "
                 "sweep." % p.name)
    box = d["box"]
    return {"path": p, "design": d, "n": int(d["n"]), "seed": int(d["seed"]), "f_fixed": float(d["f_fixed"]),
            "ks_lo": float(box["Ks_mult"]["lo"]), "ks_hi": float(box["Ks_mult"]["hi"]),
            "cc_lo": float(box["channelconductivity_mmhr"]["lo"]), "cc_hi": float(box["channelconductivity_mmhr"]["hi"]),
            "routing": d.get("routing_pinned")}


# ------------------------------------------------------------------
# Analysis (pure functions of the two metric tables)
# ------------------------------------------------------------------
def cc_band_table(m_on, off_best, nb):
    """The ON runs in nb bands of equal width in log cc: the best KGE of each, against the best OFF run."""
    ok = m_on[np.isfinite(m_on["kge_2012"])]
    u = RG.cc_unit(ok["channelconductivity_mmhr"].to_numpy(float))
    lo_log, hi_log = RG.cc_box_logs()
    rows = []
    for b in range(nb):
        lo_u, hi_u = b / float(nb), (b + 1) / float(nb)
        sel = ok[(u >= lo_u - 1e-12) & ((u <= hi_u + 1e-9) if b == nb - 1 else (u < hi_u - 1e-12))]
        row = {"band": b + 1, "cc_from": 10 ** (lo_log + lo_u * (hi_log - lo_log)),
               "cc_to": 10 ** (lo_log + hi_u * (hi_log - lo_log)), "n": int(len(sel))}
        if len(sel):
            best = sel.loc[sel["kge_2012"].idxmax()]
            row.update({"best_kge_2012": float(best["kge_2012"]), "median_kge_2012": float(sel["kge_2012"].median()),
                        "cc_of_best": float(best["channelconductivity_mmhr"]), "Ks_of_best": float(best["Ks_mult"]),
                        "pbias_of_best_pct": float(best["pbias_pct"]),
                        "best_minus_off_best": float(best["kge_2012"] - off_best)})
        rows.append(row)
    return pd.DataFrame(rows)


def ridge_fit(m_on):
    """Two readings of the Ks-cc trade-off, from the ON runs only.
    (1) Volume: PBIAS = a + b*Ks + c*log10(cc) by least squares over all ON runs; the Ks that gives PBIAS 0 at a given cc.
    (2) The best ON runs: how their Ks changes with log10(cc)."""
    ok = m_on[np.isfinite(m_on["kge_2012"]) & np.isfinite(m_on["pbias_pct"])]
    out = {"n_fit": int(len(ok))}
    nan = float("nan")
    out.update({"pbias_intercept": nan, "pbias_per_Ks": nan, "pbias_per_decade_cc": nan, "pbias_fit_r2": nan,
                "Ks_pbias0_at_cc_lo": nan, "Ks_pbias0_at_cc_hi": nan, "volume_line_drawn": False,
                "n_ridge": 0, "ridge_slope_Ks_per_decade": nan, "ridge_rho": nan, "ridge_p": nan})
    if len(ok) >= 8 and ok["Ks_mult"].nunique() > 2 and ok["channelconductivity_mmhr"].nunique() > 2:
        X = np.column_stack([np.ones(len(ok)), ok["Ks_mult"].to_numpy(float),
                             np.log10(ok["channelconductivity_mmhr"].to_numpy(float))])
        y = ok["pbias_pct"].to_numpy(float)
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        ss_res = float(np.sum((y - X @ coef) ** 2))
        ss_tot = float(np.sum((y - y.mean()) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else nan
        a, b, c = [float(v) for v in coef]
        out.update({"pbias_intercept": a, "pbias_per_Ks": b, "pbias_per_decade_cc": c, "pbias_fit_r2": r2})
        lo_cc, hi_cc = RS.BOX["channelconductivity_mmhr"][:2]
        if b < 0 and np.isfinite(r2) and r2 >= VOL_R2_MIN:
            out["Ks_pbias0_at_cc_lo"] = float((-a - c * np.log10(lo_cc)) / b)
            out["Ks_pbias0_at_cc_hi"] = float((-a - c * np.log10(hi_cc)) / b)
            out["volume_line_drawn"] = True
    n_r = int(min(max(RIDGE_MIN, int(round(RIDGE_FRAC * len(ok)))), len(ok)))
    top = ok.nlargest(n_r, "kge_2012")
    out["n_ridge"] = int(len(top))
    lcc = np.log10(top["channelconductivity_mmhr"].to_numpy(float))
    if len(top) >= 8 and np.ptp(lcc) > 0 and np.ptp(top["Ks_mult"].to_numpy(float)) > 0:
        rho = stats.spearmanr(lcc, top["Ks_mult"].to_numpy(float))
        out["ridge_slope_Ks_per_decade"] = float(np.polyfit(lcc, top["Ks_mult"].to_numpy(float), 1)[0])
        out["ridge_rho"] = float(rho[0])
        out["ridge_p"] = float(rho[1])
    return out


def low_cc_check(m_on, m_off):
    """ON runs at the very lowest cc lose almost no water, so they should fit about like OFF at the same Ks.
    Compares each such ON run with the OFF curve (KGE against Ks, straight lines between the OFF runs)."""
    ok_on = m_on[np.isfinite(m_on["kge_2012"])]
    ok_off = m_off[np.isfinite(m_off["kge_2012"])].sort_values("Ks_mult")
    out = {"n": 0, "median_diff": float("nan"), "min_diff": float("nan"), "max_diff": float("nan"),
           "cc_upper_mmhr": float("nan")}
    lo_log, hi_log = RG.cc_box_logs()
    out["cc_upper_mmhr"] = float(10 ** (lo_log + LOW_CC_FRAC * (hi_log - lo_log)))
    if len(ok_off) < 5:
        return out
    ks_off = ok_off["Ks_mult"].to_numpy(float)
    kge_off = ok_off["kge_2012"].to_numpy(float)
    sel = ok_on[(RG.cc_unit(ok_on["channelconductivity_mmhr"].to_numpy(float)) < LOW_CC_FRAC)
                & (ok_on["Ks_mult"] >= ks_off.min()) & (ok_on["Ks_mult"] <= ks_off.max())]
    if len(sel) < 3:
        return out
    d = sel["kge_2012"].to_numpy(float) - np.interp(sel["Ks_mult"].to_numpy(float), ks_off, kge_off)
    out.update({"n": int(len(sel)), "median_diff": float(np.median(d)), "min_diff": float(d.min()), "max_diff": float(d.max())})
    return out


def analyze_arm(ar, args, boxes):
    m_on, m_off = ar["m_on"].copy(), ar["m_off"].copy()
    m_on["u_cc"] = RG.cc_unit(m_on["channelconductivity_mmhr"].to_numpy(float))
    m_on["u_Ks"] = u_lin(m_on["Ks_mult"], boxes["on"]["ks_lo"], boxes["on"]["ks_hi"])
    m_off["u_Ks"] = u_lin(m_off["Ks_mult"], boxes["off"]["ks_lo"], boxes["off"]["ks_hi"])
    ok_on = m_on[np.isfinite(m_on["kge_2012"])]
    ok_off = m_off[np.isfinite(m_off["kge_2012"])]
    n_top_on = top_count(len(ok_on), args.top_frac)
    n_top_off = top_count(len(ok_off), args.top_frac)
    on_row = RG.series_summary(m_on, SER_ON, n_top_on)
    off_row = RG.series_summary(m_off, SER_OFF, n_top_off)
    gap = on_row["best_kge_2012"] - off_row["best_kge_2012"]
    gap_mean = on_row["top_n_mean_kge"] - off_row["top_n_mean_kge"]
    loc = RG.cc_location(m_on, n_top_on, args.floor_tol)
    n_wide = int(min(2 * n_top_on, len(ok_on)))
    loc_wide = RG.cc_location(m_on, n_wide, args.floor_tol)
    bands = cc_band_table(m_on, off_row["best_kge_2012"], args.cc_bands)
    bb = bands.dropna(subset=["best_kge_2012"])
    band_spread = float(bb["best_kge_2012"].max() - bb["best_kge_2012"].min()) if len(bb) >= 3 else float("nan")
    flat = bool(np.isfinite(band_spread) and band_spread <= FLAT_BANDS_TOL)

    def read_row(g, lc):
        r = RG.decision_row(g, args.margin, lc)
        return 5 if (flat and r == 1) else r
    row_best = read_row(gap, loc)
    row_mean = read_row(gap_mean, loc)
    within = bool(abs(gap) <= args.margin)
    row_wide = read_row(gap, loc_wide) if within else row_best
    flat_override = bool(flat and RG.decision_row(gap, args.margin, loc) == 1)
    near_margin = bool((args.margin - NEAR_MARGIN_TOL) <= gap <= args.margin)
    ridge = ridge_fit(m_on)
    lowcc = low_cc_check(m_on, m_off)
    near_best_off = ok_off[ok_off["kge_2012"] >= off_row["best_kge_2012"] - RIDGE_NEAR]
    diff_row = {"series": "ON minus OFF"}
    for key in ("best_kge_2012", "top_n_mean_kge", "top_n_median_kge", "median_kge_all_runs", "best_pbias_pct",
                "best_kge_r", "best_kge_beta", "best_kge_gamma"):
        diff_row[key] = on_row[key] - off_row[key]
    sdf = pd.concat([pd.DataFrame([on_row, off_row]), pd.DataFrame([diff_row])], ignore_index=True)
    # edges of the Ks boxes
    edge_notes = []
    for nm, row, bx in (("ON", on_row, boxes["on"]), ("OFF", off_row, boxes["off"])):
        u = float(u_lin(row["best_Ks"], bx["ks_lo"], bx["ks_hi"]))
        if u < 0.05 or u > 0.95:
            edge_notes.append("the best %s run is at the %s edge of its Ks box (Ks %.3f in %g - %g)%s"
                              % (nm, "lower" if u < 0.5 else "upper", row["best_Ks"], bx["ks_lo"], bx["ks_hi"],
                                 ": OFF may be UNDER-stated, which would favour ON" if nm == "OFF"
                                 else ": the Ks range may be clipping the ridge"))
    return {"m_on": m_on, "m_off": m_off, "on_row": on_row, "off_row": off_row, "sdf": sdf, "gap": gap,
            "gap_mean": gap_mean, "loc": loc, "loc_wide": loc_wide, "row_best": row_best, "row_mean": row_mean,
            "row_wide": row_wide, "within_margin": within, "flat": flat, "flat_override": flat_override,
            "near_margin": near_margin, "band_spread": band_spread, "bands": bands, "ridge": ridge, "lowcc": lowcc, "n_top_on": n_top_on, "n_top_off": n_top_off,
            "off_near_ks": (float(near_best_off["Ks_mult"].min()), float(near_best_off["Ks_mult"].max()))
            if len(near_best_off) else (float("nan"), float("nan")),
            "edge_notes": edge_notes}


# ------------------------------------------------------------------
# Printing
# ------------------------------------------------------------------
def print_arm(arm, is_primary, an, n_scored, boxes, partial_note, args):
    on_row, off_row = an["on_row"], an["off_row"]
    print("\n" + "#" * 110)
    print("ARM: %s   %s%s" % (arm.upper(), "(PRIMARY: fixed before the run)" if is_primary else "(sensitivity arm)", partial_note))
    print("     %d real 5-minute points scored.  f held at %g.  Ks box: OFF %g - %g, ON %g - %g.  cc box (ON): %g - %g."
          % (n_scored, boxes["f_fixed"], boxes["off"]["ks_lo"], boxes["off"]["ks_hi"],
             boxes["on"]["ks_lo"], boxes["on"]["ks_hi"], boxes["on"]["cc_lo"], boxes["on"]["cc_hi"]))
    print("     %s" % RG.ARM_TEXT[arm])
    print("#" * 110)

    print("\n" + "=" * 110)
    print("HEAD TO HEAD, %s ARM: the best run of each series (each series free in Ks)" % arm.upper())
    print("=" * 110)
    print("%-36s %14s %14s %12s" % ("", "cc ON", "cc OFF", "ON - OFF"))

    def line(label, key, spec):
        a, b = on_row[key], off_row[key]
        d = a - b if (np.isfinite(a) and np.isfinite(b)) else np.nan
        print("%-36s %14s %14s %12s" % (label, fnum(a, spec), fnum(b, spec), fnum(d, "+" + spec)))
    line("best KGE_2012", "best_kge_2012", ".4f")
    print("%-36s %14s %14s %12s" % ("runs in the series", on_row["n_valid_kge"], off_row["n_valid_kge"], ""))
    print("%-36s %14s %14s %12s" % ("mean KGE_2012 of the best 10%% (%d / %d)" % (an["n_top_on"], an["n_top_off"]),
                                    fnum(on_row["top_n_mean_kge"], ".4f"), fnum(off_row["top_n_mean_kge"], ".4f"),
                                    fnum(an["gap_mean"], "+.4f")))
    line("median KGE_2012, all runs", "median_kge_all_runs", ".4f")
    print("\nBest run in each series:")
    print("  %-6s %-26s %6s %9s | %7s %6s %6s %6s %9s %8s"
          % ("series", "run_id (SMF_20140812_110_...)", "Ks", "cc", "PBIAS%", "r", "beta", "gamma", "peak err%", "peak hr"))
    for nm, row in (("ON", on_row), ("OFF", off_row)):
        print("  %-6s %-26s %6.2f %9s | %7.1f %6.3f %6.3f %6.3f %9.1f %8.2f"
              % (nm, str(row["best_run_id"]).replace("SMF_20140812_110_", ""), row["best_Ks"], RG.fcc(row["best_cc"]) if nm == "ON" else "inert",
                 row["best_pbias_pct"], row["best_kge_r"], row["best_kge_beta"], row["best_kge_gamma"],
                 row["best_peak_error_pct"], row["best_peak_timing_error_hr"]))
    print("  (r, beta, gamma are the three parts of KGE; 1 is perfect. peak hr = simulated peak time minus real peak time; "
          "positive = the model peaks late.)")
    print("  OFF runs within %.2f KGE of the best OFF run cover Ks %.2f - %.2f." % (RIDGE_NEAR, an["off_near_ks"][0], an["off_near_ks"][1]))

    b = an["bands"]
    print("\n" + "=" * 110)
    print("BY BAND OF cc (ON runs; equal widths in log cc), %s ARM" % arm.upper())
    print("=" * 110)
    print("  %4s %19s %4s | %9s %9s | %9s %8s %8s | %12s"
          % ("band", "cc from - to", "n", "best KGE", "median", "cc of best", "Ks best", "PBIAS%", "best - OFF best"))
    for _, r in b.iterrows():
        rng = "%s - %s" % (RG.fcc(r["cc_from"]), RG.fcc(r["cc_to"]))
        if r["n"] == 0:
            print("  %4d %19s %4d |   (no runs)" % (r["band"], rng, r["n"]))
        else:
            print("  %4d %19s %4d | %9.4f %9.4f | %9s %8.3f %8.1f | %+12.4f"
                  % (r["band"], rng, r["n"], r["best_kge_2012"], r["median_kge_2012"], RG.fcc(r["cc_of_best"]),
                     r["Ks_of_best"], r["pbias_of_best_pct"], r["best_minus_off_best"]))
    print("  (a negative last column means that band's best ON run is worse than the best OFF run.)")
    if np.isfinite(an["band_spread"]):
        print("  Spread of the band bests (highest minus lowest): %.4f.%s" % (an["band_spread"],
              ("  That is within %.2f: the best fit is about the same at every cc, so this event does not pin cc." % FLAT_BANDS_TOL)
              if an["band_spread"] <= FLAT_BANDS_TOL else
              "  (Flat within %.2f would mean the event does not pin cc; a larger spread can also come from the sparse sampling of the ridge.)"
              % FLAT_BANDS_TOL))

    RG.print_cc_location([an["loc"], an["loc_wide"]] if an["within_margin"] else [an["loc"]], an["n_top_on"], arm)

    rd = an["ridge"]
    print("\n" + "=" * 110)
    print("THE Ks-cc TRADE-OFF (ON runs), %s ARM" % arm.upper())
    print("=" * 110)
    if np.isfinite(rd["pbias_per_decade_cc"]):
        print("  Volume: at a fixed Ks, every tenfold rise in cc changes PBIAS by %+.1f points, and every +1 in Ks by %+.1f points "
              "(straight-line fit to %d ON runs, R2 %.2f)." % (rd["pbias_per_decade_cc"], rd["pbias_per_Ks"], rd["n_fit"], rd["pbias_fit_r2"]))
        print("  The old 500-run sweep (f 0.004 - 0.03) gave %+.1f points per tenfold rise in cc." % OLD_PBIAS_PER_DECADE)
        if rd["volume_line_drawn"]:
            print("  The Ks that would give PBIAS 0 falls from %.2f at cc %g to %.2f at cc %g."
                  % (rd["Ks_pbias0_at_cc_lo"], RS.BOX["channelconductivity_mmhr"][0], rd["Ks_pbias0_at_cc_hi"],
                     RS.BOX["channelconductivity_mmhr"][1]))
        else:
            print("  (The straight-line volume fit is not good enough (R2 under %.2f, or Ks has no downward effect) to quote or "
                  "draw a volume-matching Ks.)" % VOL_R2_MIN)
        if rd["pbias_per_decade_cc"] > -5.0:
            print("  NOTE: cc hardly changes the volume in these ON runs. If channel loss were really acting, PBIAS would fall "
                  "with cc. Check that the ON runs used optpercolation 1 in their input files before reading anything else.")
    else:
        print("  Volume fit not possible (too few ON runs or no spread in Ks or cc).")
    if np.isfinite(rd["ridge_slope_Ks_per_decade"]):
        print("  Best %d ON runs: Ks changes by %+.2f per tenfold rise in cc (rank correlation %+.2f, p %s). A negative value "
              "means the best Ks falls as cc rises: the ridge. A guide only: the best runs are not independent."
              % (rd["n_ridge"], rd["ridge_slope_Ks_per_decade"], rd["ridge_rho"], RG.fmt_p(rd["ridge_p"])))
    lc = an["lowcc"]
    if lc["n"] >= 3:
        print("  Consistency: the %d ON runs with cc up to %s (the lowest tenth of the range) and Ks inside the OFF range differ "
              "from the OFF curve at the same Ks by a median of %+.3f KGE (range %+.3f to %+.3f)."
              % (lc["n"], RG.fcc(lc["cc_upper_mmhr"]), lc["median_diff"], lc["min_diff"], lc["max_diff"]))
        if abs(lc["median_diff"]) > LOW_CC_WARN_TOL:
            print("  *** WARNING: ON at its lowest cc should fit about like OFF at the same Ks, but the median differs by %.3f. "
                  "Check that the two series differ only in optpercolation (same f, same routing, same label). ***"
                  % abs(lc["median_diff"]))
        elif abs(lc["median_diff"]) > LOW_CC_NOTE_TOL:
            print("  NOTE: that is larger than the %.2f expected from the small loss at the lowest cc." % LOW_CC_NOTE_TOL)
    else:
        print("  Consistency check (ON at lowest cc against the OFF curve): too few ON runs there to say.")

    print("\n" + "=" * 110)
    print("READING (rule fixed before the run), %s ARM   (margin %g KGE_2012; ALPHA %g, STABLE_SPAN %g, floor within %g mm/hr)"
          % (arm.upper(), args.margin, RG.ALPHA, RG.STABLE_SPAN, args.floor_tol))
    print("=" * 110)
    loc = an["loc"]
    lw = an["loc_wide"]
    print("  Best ON minus best OFF                         : %+.4f   -> row %d" % (an["gap"], an["row_best"]))
    print("  Mean of best 10%% (ON %d, OFF %d) ON minus OFF   : %+.4f   -> row %d" % (an["n_top_on"], an["n_top_off"], an["gap_mean"], an["row_mean"]))
    print("  Best %d ON runs: lowest third p = %s, highest third p = %s, span %.2f of the log cc range"
          % (loc["top_n"], RG.fmt_p(loc["p_lowest_third"]), RG.fmt_p(loc["p_highest_third"]), loc["span_frac_of_box"]))
    if an["within_margin"]:
        print("  Best %d ON runs (location check): lowest third p = %s, highest third p = %s, span %.2f   -> row %d"
              % (lw["top_n"], RG.fmt_p(lw["p_lowest_third"]), RG.fmt_p(lw["p_highest_third"]), lw["span_frac_of_box"], an["row_wide"]))
    else:
        print("  (The wider location check is used only when the gap is within the margin.)")
    if an["flat_override"]:
        print("\n  FLAT OVERRIDE: by the counts the best ON runs crowd the lowest third (row 1), but the best ON run of every cc band is "
              "within %.2f of the others (spread %.4f), so that crowding is not evidence for no loss. Read as row 5." % (FLAT_BANDS_TOL, an["band_spread"]))
    print("\n  Row %d: %s" % (an["row_best"], READINGS[an["row_best"]]))
    rows_set = sorted(set([an["row_best"], an["row_mean"], an["row_wide"]]))
    if len(rows_set) == 1:
        print("\n  STABILITY: %s point to row %d."
              % ("the best run, the mean of the best 10%% and the wider location check all" if an["within_margin"]
                 else "the best run and the mean of the best 10% both", an["row_best"]))
    else:
        print("\n  STABILITY: the readings do NOT agree (best run row %d, mean of the best 10%% row %d%s). Treat the answer as UNSETTLED. "
              "The other row(s):" % (an["row_best"], an["row_mean"],
                                      (", location from the best 20%% row %d" % an["row_wide"]) if an["within_margin"] else ""))
        for r_other in rows_set:
            if r_other != an["row_best"]:
                print("        row %d: %s" % (r_other, READINGS[r_other]))
    if an["near_margin"]:
        print("\n  NOTE (near the margin): the gap is within %.2f below the margin of %g. The ON sweep samples a plane with %d "
              "points, OFF a line with %d, so ON's best point is a slight under-estimate of what ON can reach. A refinement of "
              "the ON sweep around its best point would settle it; do not call this either way yet."
              % (NEAR_MARGIN_TOL, args.margin, an["on_row"]["n_valid_kge"], an["off_row"]["n_valid_kge"]))
    for note in an["edge_notes"]:
        print("\n  NOTE: %s." % note)
    print("\n  Caveats: one event; routing pinned at the synthetic-truth values; the SMPHQ rainfall artifact is in every run; gauge "
          "error is not quantified; the filled arm assumes a silent logger means no change; f held at one value (%g)." % boxes["f_fixed"])


# ------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------
def _legend_below(fig, handles, ncol=2):
    fig.legend(handles=handles, loc="lower center", ncol=ncol, frameon=False, fontsize=8, labelcolor=RG.INK2)


def _title(ax, text, notes):
    notes = [n for n in notes if n]
    if notes:
        text = text + "\n" + "   ".join(notes)
    ax.set_title(text, color=RG.INK, fontsize=10, loc="left")


def fig_kge_vs_cc(an, margin, path, arm, partial_note, label, boxes):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    m_on = an["m_on"]
    ok = m_on[np.isfinite(m_on["kge_2012"])]
    off_best = an["off_row"]["best_kge_2012"]
    ymax = max(float(ok["kge_2012"].max()), off_best + margin)
    ymin = max(float(ok["kge_2012"].min()), ymax - 0.8)
    pad = 0.05 * (ymax - ymin)
    hidden = int((ok["kge_2012"] < ymin).sum())
    lo_log, hi_log = RG.cc_box_logs()
    edge = 10 ** (lo_log + (hi_log - lo_log) / 3.0)
    lo, hi = RS.BOX["channelconductivity_mmhr"][:2]
    fig, ax = plt.subplots(figsize=(9.0, 5.4), facecolor=RG.SURFACE)
    RG._style(ax)
    ax.axvspan(lo * 0.8, edge, color=RG.NULLC, alpha=0.13, linewidth=0, zorder=1)
    ax.axhline(off_best, color=RG.NULLC, linewidth=1.8, zorder=2)
    ax.axhline(off_best + margin, color=RG.INK, linewidth=1.3, linestyle="--", zorder=2)
    ax.scatter(ok["channelconductivity_mmhr"], ok["kge_2012"], s=22, marker="o", color=RG.C_CC, alpha=0.7, linewidths=0, zorder=3)
    top = ok.nlargest(an["n_top_on"], "kge_2012")
    ax.scatter(top["channelconductivity_mmhr"], top["kge_2012"], s=70, marker="o", color=RG.C_CC, edgecolors=RG.INK,
               linewidths=1.4, zorder=5)
    bb = an["bands"].dropna(subset=["best_kge_2012"])
    if len(bb) >= 2:
        xs = np.sqrt(bb["cc_from"].to_numpy(float) * bb["cc_to"].to_numpy(float))
        ax.plot(xs, bb["best_kge_2012"], color=RG.INK, linewidth=1.6, marker="D", markersize=5, zorder=4)
    RG._cc_axis(ax)
    ax.set_ylim(ymin - pad, ymax + pad)
    ax.set_xlabel("channel conductivity cc (mm/hr), log scale", color=RG.INK2, fontsize=9)
    ax.set_ylabel("KGE_2012 against the real gauge (%s arm)" % arm, color=RG.INK2, fontsize=9)
    handles = [
        Line2D([0], [0], marker="o", color=RG.C_CC, linestyle="none", markersize=5, alpha=0.7,
               label="ON runs, Ks free in %g - %g (%d)" % (boxes["on"]["ks_lo"], boxes["on"]["ks_hi"], len(ok))),
        Line2D([0], [0], marker="o", color=RG.C_CC, markeredgecolor=RG.INK, markeredgewidth=1.4, linestyle="none", markersize=8,
               label="the %d best ON runs" % len(top)),
        Line2D([0], [0], color=RG.INK, linewidth=1.6, marker="D", markersize=5, label="best ON run in each cc band"),
        Line2D([0], [0], color=RG.NULLC, linewidth=1.8, label="best OFF run (%.4f)" % off_best),
        Line2D([0], [0], color=RG.INK, linewidth=1.3, linestyle="--", label="best OFF + margin %g: ON must clear this line" % margin),
        Patch(facecolor=RG.NULLC, alpha=0.13, label="lowest third of the cc range"),
    ]
    _legend_below(fig, handles, ncol=2)
    _title(ax, "%s arm, sweep %s (f held at one value): does any ON run clear the bar?" % (arm.upper(), label),
           [partial_note.strip(), "%d run(s) below the y-range, not drawn" % hidden if hidden else ""])
    fig.tight_layout(rect=(0, 0.17, 1, 1))
    fig.savefig(path, dpi=150, facecolor=RG.SURFACE)
    plt.close(fig)


def fig_kge_vs_ks(an, margin, path, arm, partial_note, label):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    m_on, m_off = an["m_on"], an["m_off"]
    ok = m_on[np.isfinite(m_on["kge_2012"])]
    okf = m_off[np.isfinite(m_off["kge_2012"])].sort_values("Ks_mult")
    ymax = max(float(ok["kge_2012"].max()), float(okf["kge_2012"].max()))
    ymin = max(min(float(ok["kge_2012"].min()), float(okf["kge_2012"].min())), ymax - 0.8)
    pad = 0.05 * (ymax - ymin)
    hidden = int((ok["kge_2012"] < ymin).sum()) + int((okf["kge_2012"] < ymin).sum())
    lowcc = ok[ok["u_cc"] < 1.0 / 3.0]
    rest = ok[ok["u_cc"] >= 1.0 / 3.0]
    fig, ax = plt.subplots(figsize=(9.0, 5.4), facecolor=RG.SURFACE)
    RG._style(ax)
    ax.plot(okf["Ks_mult"], okf["kge_2012"], color=RG.NULLC, linewidth=1.4, zorder=2)
    ax.scatter(okf["Ks_mult"], okf["kge_2012"], s=34, marker="s", facecolors=RG.SURFACE, edgecolors=RG.INK2, linewidths=1.2, zorder=4)
    ax.scatter(rest["Ks_mult"], rest["kge_2012"], s=24, marker="o", color=RG.C_CC, alpha=0.75, linewidths=0, zorder=3)
    ax.scatter(lowcc["Ks_mult"], lowcc["kge_2012"], s=30, marker="o", facecolors=RG.SURFACE, edgecolors=RG.C_CC, linewidths=1.3, zorder=3)
    top = ok.nlargest(an["n_top_on"], "kge_2012")
    ax.scatter(top["Ks_mult"], top["kge_2012"], s=75, marker="o", facecolors="none", edgecolors=RG.INK, linewidths=1.4, zorder=5)
    ax.set_ylim(ymin - pad, ymax + pad)
    ax.set_xlabel("Ks_mult", color=RG.INK2, fontsize=9)
    ax.set_ylabel("KGE_2012 against the real gauge (%s arm)" % arm, color=RG.INK2, fontsize=9)
    handles = [
        Line2D([0], [0], marker="s", color=RG.NULLC, markerfacecolor=RG.SURFACE, markeredgecolor=RG.INK2, markeredgewidth=1.2,
               linewidth=1.4, markersize=6, label="channel loss OFF (%d runs; line joins them in order of Ks)" % len(okf)),
        Line2D([0], [0], marker="o", color=RG.C_CC, linestyle="none", markersize=5.5, alpha=0.75,
               label="ON, cc in the middle and upper thirds (%d)" % len(rest)),
        Line2D([0], [0], marker="o", color=RG.C_CC, markerfacecolor=RG.SURFACE, markeredgewidth=1.3, linestyle="none", markersize=6,
               label="ON, cc in the lowest third (%d): the least channel loss" % len(lowcc)),
        Line2D([0], [0], marker="o", color=RG.INK, markerfacecolor="none", markeredgewidth=1.4, linestyle="none", markersize=8,
               label="the %d best ON runs" % len(top)),
    ]
    _legend_below(fig, handles, ncol=1)
    _title(ax, "%s arm, sweep %s: fit against Ks, ON and OFF" % (arm.upper(), label),
           [partial_note.strip(), "%d run(s) below the y-range, not drawn" % hidden if hidden else ""])
    fig.tight_layout(rect=(0, 0.20, 1, 1))
    fig.savefig(path, dpi=150, facecolor=RG.SURFACE)
    plt.close(fig)


def fig_ridge(an, path, arm, partial_note, label, boxes):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    m_on = an["m_on"]
    ok = m_on[np.isfinite(m_on["kge_2012"])]
    best = float(ok["kge_2012"].max())
    near = ok[ok["kge_2012"] >= best - RIDGE_NEAR]
    top = ok.nlargest(an["n_top_on"], "kge_2012")
    rd = an["ridge"]
    fig, ax = plt.subplots(figsize=(9.0, 5.4), facecolor=RG.SURFACE)
    RG._style(ax)
    k0, k1 = an["off_near_ks"]
    if np.isfinite(k0) and np.isfinite(k1):
        ax.axhspan(k0, k1, color=RG.C_KS, alpha=0.14, linewidth=0, zorder=1)
    lo, hi = RS.BOX["channelconductivity_mmhr"][:2]
    if rd["volume_line_drawn"]:
        ax.plot([lo, hi], [rd["Ks_pbias0_at_cc_lo"], rd["Ks_pbias0_at_cc_hi"]], color=RG.INK2, linewidth=1.5, linestyle=":", zorder=2)
    rest = ok[~ok.index.isin(near.index)]
    ax.scatter(rest["channelconductivity_mmhr"], rest["Ks_mult"], s=16, color=RG.NULLC, alpha=0.4, linewidths=0, zorder=3)
    ax.scatter(near["channelconductivity_mmhr"], near["Ks_mult"], s=40, color=RG.C_CC, alpha=0.9, linewidths=0, zorder=4)
    ax.scatter(top["channelconductivity_mmhr"], top["Ks_mult"], s=80, marker="o", facecolors="none", edgecolors=RG.INK, linewidths=1.4, zorder=5)
    RG._cc_axis(ax)
    ax.set_ylim(boxes["on"]["ks_lo"] - 0.2, boxes["on"]["ks_hi"] + 0.2)
    ax.set_xlabel("channel conductivity cc (mm/hr), log scale", color=RG.INK2, fontsize=9)
    ax.set_ylabel("Ks_mult", color=RG.INK2, fontsize=9)
    handles = [
        Line2D([0], [0], marker="o", color=RG.NULLC, alpha=0.5, linestyle="none", markersize=5, label="other ON runs"),
        Line2D([0], [0], marker="o", color=RG.C_CC, linestyle="none", markersize=6.5,
               label="ON runs within %.2f of the best ON run (%d)" % (RIDGE_NEAR, len(near))),
        Line2D([0], [0], marker="o", color=RG.INK, markerfacecolor="none", markeredgewidth=1.4, linestyle="none", markersize=8,
               label="the %d best ON runs" % len(top)),
        Patch(facecolor=RG.C_KS, alpha=0.14, label="Ks range where OFF fits within %.2f of its best (%.2f - %.2f)" % (RIDGE_NEAR, k0, k1)),
    ]
    if rd["volume_line_drawn"]:
        handles.append(Line2D([0], [0], color=RG.INK2, linewidth=1.5, linestyle=":",
                              label="Ks that gives PBIAS 0 (straight-line fit to the ON runs)"))
    _legend_below(fig, handles, ncol=1)
    _title(ax, "%s arm, sweep %s: where the good ON runs put Ks and cc" % (arm.upper(), label), [partial_note.strip()])
    fig.tight_layout(rect=(0, 0.22, 1, 1))
    fig.savefig(path, dpi=150, facecolor=RG.SURFACE)
    plt.close(fig)


def fig_hydrograph_zoom(real_scored, record_pts, sim_on, sim_off, on_row, off_row, path, arm):
    """Same curves as the full-window figure, zoomed on the flood, with the three peak times in the legend."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    series = {"real": real_scored, "on": sim_on, "off": sim_off}
    pk = max(float(np.nanmax(v.values)) for v in series.values() if len(v) and np.isfinite(v.values).any())
    live = None
    for v in series.values():
        idx = v.index[np.asarray(v.values > 0.02 * pk)]
        if len(idx):
            live = (idx[0], idx[-1]) if live is None else (min(live[0], idx[0]), max(live[1], idx[-1]))
    if live is None:
        return
    z0, z1 = live[0] - pd.Timedelta(minutes=30), live[1] + pd.Timedelta(minutes=30)

    def pk_time(v):
        vv = v.dropna()
        return vv.idxmax().strftime("%H:%M") if len(vv) else "n/a"
    fig, ax = plt.subplots(figsize=(9.4, 4.8), facecolor=RG.SURFACE)
    RG._style(ax)
    ax.plot(sim_on.index, sim_on.values, color=RG.C_CC, linewidth=2.2,
            label="best cc-ON run (KGE_2012 %.3f, cc %s mm/hr; peak at %s)" % (on_row["best_kge_2012"], RG.fcc(on_row["best_cc"]), pk_time(sim_on)))
    ax.plot(sim_off.index, sim_off.values, color=RG.INK2, linewidth=1.5, linestyle=(0, (4, 3)),
            label="best cc-OFF run (KGE_2012 %.3f; peak at %s)" % (off_row["best_kge_2012"], pk_time(sim_off)))
    if arm == "asis":
        ax.plot(real_scored.index, real_scored.values, color=RG.INK, linewidth=0, marker="o", markersize=4.0,
                label="real gauge, as-is: the bins the scorer keeps (peak at %s)" % pk_time(real_scored))
    else:
        ax.plot(real_scored.index, real_scored.values, color=RG.INK, linewidth=1.5,
                label="real gauge, filled (peak at %s)" % pk_time(real_scored))
        if record_pts is not None and len(record_pts):
            ax.plot(record_pts.index, record_pts.values, color=RG.INK, linewidth=0, marker="o", markersize=3.5,
                    label="bins that hold a gauge record")
    ax.set_xlim(z0, z1)
    ax.set_ylim(-0.03 * pk, 1.08 * pk)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.set_xlabel("time of day (hour:minute), zoom on the flood, 12 August 2014", color=RG.INK2, fontsize=9)
    ax.set_ylabel("discharge at the outlet (m$^3$/s)", color=RG.INK2, fontsize=9)
    ax.legend(loc="upper right", frameon=False, fontsize=8.5, labelcolor=RG.INK2)
    ax.set_title("%s arm: the flood in detail, real gauge against the best run of each series" % arm.upper(),
                 color=RG.INK, fontsize=10, loc="left")
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor=RG.SURFACE)
    plt.close(fig)


# ------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Series 110 Stage B: score the pinned-f ON and OFF sweeps against the REAL gauge and compare the best "
                    "run of each (each free in Ks).")
    ap.add_argument("--label", required=True, help="the --label of the two sweeps (for example sb1)")
    ap.add_argument("--margin", type=float, default=DEFAULT_MARGIN,
                    help="KGE_2012 margin ON must beat OFF by (default %g, your Stage 0 margin)" % DEFAULT_MARGIN)
    ap.add_argument("--primary_arm", choices=["asis", "filled"], default="filled",
                    help="the arm whose reading is the verdict (default filled, as decided 2026-10-04)")
    ap.add_argument("--gauge_mode", choices=["both", "asis", "filled"], default="both")
    ap.add_argument("--max_gap_min", type=float, default=RG.DEFAULT_MAX_GAP_MIN)
    ap.add_argument("--max_hold_hr", type=float, default=RG.DEFAULT_MAX_HOLD_HR)
    ap.add_argument("--top_frac", type=float, default=TOP_FRAC, help="share of each series counted as its 'best runs' (default 0.10)")
    ap.add_argument("--cc_bands", type=int, default=5, help="number of cc bands, equal widths in log cc (default 5)")
    ap.add_argument("--floor_tol", type=float, default=RG.FLOOR_TOL_MMHR)
    ap.add_argument("--check_only", action="store_true", help="run the gates and stop; writes nothing")
    ap.add_argument("--min_points", type=int, default=30)
    ap.add_argument("--drop_bad_runs", action="store_true", help="exclude runs that fail gate 2A, 2B or 4 and continue")
    ap.add_argument("--no_plots", action="store_true")
    ap.add_argument("--on_csv", type=Path, default=None)
    ap.add_argument("--off_csv", type=Path, default=None)
    ap.add_argument("--gauge_xlsx", type=Path, default=None)
    ap.add_argument("--out_dir", type=Path, default=None)
    args = ap.parse_args()
    if not (0.0 < args.margin < 1.0):
        ap.error("--margin must be between 0 and 1.")
    if not (0.02 <= args.top_frac <= 0.5):
        ap.error("--top_frac must be between 0.02 and 0.5.")
    if args.cc_bands < 3 or args.cc_bands > 12:
        ap.error("--cc_bands must be between 3 and 12.")
    if args.max_gap_min <= 0 or args.max_hold_hr <= 0:
        ap.error("--max_gap_min and --max_hold_hr must be positive.")

    script_dir = Path.cwd()
    project_root = script_dir.parent
    calib_dir = project_root / "calibration_work"
    summary_dir = calib_dir / "03_comparisons" / "summary_tables"
    csv_dir = calib_dir / "03_comparisons" / "csv_exports"
    synth_dir = calib_dir / "synth_truth"
    on_path = args.on_csv or (summary_dir / ("lhs_results_lowf_%s_ON_110.csv" % args.label))
    off_path = args.off_csv or (summary_dir / ("lhs_results_lowf_%s_OFF_110.csv" % args.label))
    xlsx_path = args.gauge_xlsx or project_root.joinpath(*RG.GAUGE_REL)
    out_dir = args.out_dir or (summary_dir / OUT_DIRNAME / args.label)

    arms_used = ["filled", "asis"] if args.gauge_mode == "both" else [args.gauge_mode]
    primary = args.primary_arm if args.gauge_mode == "both" else args.gauge_mode
    arms_used = [primary] + [a for a in arms_used if a != primary]

    print("\n" + "=" * 78)
    print("Series 110, STAGE B -- channel loss ON against OFF, f held fixed, real gauge (version 1)")
    print("label %s, margin %g, python %s, pandas %s, numpy %s" % (args.label, args.margin, sys.version.split()[0], pd.__version__, np.__version__))
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

    d_off = read_design(summary_dir, args.label, "OFF")
    d_on = read_design(summary_dir, args.label, "ON")
    ck = RG.Checks()
    if not np.isclose(d_off["f_fixed"], d_on["f_fixed"], rtol=1e-12, atol=0.0):
        sys.exit("The two design files hold f at different values (OFF %g, ON %g). Stage B compares the series at ONE f."
                 % (d_off["f_fixed"], d_on["f_fixed"]))
    if d_off["routing"] != d_on["routing"]:
        sys.exit("The two design files pin the routing parameters at different values (%s against %s)." % (d_off["routing"], d_on["routing"]))
    lo_cc, hi_cc = RS.BOX["channelconductivity_mmhr"][:2]
    if not (np.isclose(d_on["cc_lo"], lo_cc) and np.isclose(d_on["cc_hi"], hi_cc)):
        sys.exit("The ON design samples cc over %g - %g, not the standard %g - %g box this script's cc reading is built on. "
                 "Refusing rather than reading cc against the wrong box." % (d_on["cc_lo"], d_on["cc_hi"], lo_cc, hi_cc))
    boxes = {"f_fixed": d_on["f_fixed"],
             "on": {"ks_lo": d_on["ks_lo"], "ks_hi": d_on["ks_hi"], "cc_lo": d_on["cc_lo"], "cc_hi": d_on["cc_hi"], "n": d_on["n"],
                    "seed": d_on["seed"]},
             "off": {"ks_lo": d_off["ks_lo"], "ks_hi": d_off["ks_hi"], "n": d_off["n"], "seed": d_off["seed"]}}
    print("designs: %s and %s" % (d_off["path"].name, d_on["path"].name))
    print("f held at %g.  OFF: Ks %g - %g, n %d.  ON: Ks %g - %g, cc %g - %g, n %d."
          % (boxes["f_fixed"], d_off["ks_lo"], d_off["ks_hi"], d_off["n"], d_on["ks_lo"], d_on["ks_hi"], d_on["cc_lo"], d_on["cc_hi"], d_on["n"]))
    ck.add("PASS", "design files", "same label %s, same pinned f %g, same routing in both; ON cc box is the standard %g - %g"
           % (args.label, boxes["f_fixed"], lo_cc, hi_cc))

    print("\nLoading the results tables:")
    res_off = RS.load_results(off_path, "%s OFF" % args.label, 0, d_off["n"]).reset_index(drop=True)
    res_on = RS.load_results(on_path, "%s ON" % args.label, 1, d_on["n"]).reset_index(drop=True)
    if res_off.empty or res_on.empty:
        sys.exit("A results table has no runs (OFF %d, ON %d)." % (len(res_off), len(res_on)))
    partial_bits = []
    if len(res_off) < d_off["n"]:
        partial_bits.append("OFF %d of %d" % (len(res_off), d_off["n"]))
    if len(res_on) < d_on["n"]:
        partial_bits.append("ON %d of %d" % (len(res_on), d_on["n"]))
    partial = bool(partial_bits)
    partial_note = "   [PARTIAL: %s]" % ", ".join(partial_bits) if partial else ""
    if partial:
        print("  NOTE: a partial sweep (%s). The reading is labelled PARTIAL; the first runs of a design are an arbitrary subset, "
              "not a stratified one." % ", ".join(partial_bits))

    shared = set(res_on["run_id"]) & set(res_off["run_id"])
    if shared:
        ck.add("FAIL", "run ids", "%d run id(s) appear in both series, e.g. %s. They would share one compare CSV."
               % (len(shared), ", ".join(sorted(shared)[:3])))
    else:
        ck.add("PASS", "run ids", "no run id is in both series")
    tolf = 1e-9
    for nm, r in (("OFF", res_off), ("ON", res_on)):
        bad_f = r[~np.isclose(r["f_RS_abs"].to_numpy(float), boxes["f_fixed"], rtol=1e-9, atol=0.0)]
        if len(bad_f):
            ck.add("FAIL", "f held fixed (%s)" % nm, "%d run(s) have f other than %g, e.g. %g. This is not a pinned-f sweep."
                   % (len(bad_f), boxes["f_fixed"], float(bad_f["f_RS_abs"].iloc[0])))
        else:
            ck.add("PASS", "f held fixed (%s)" % nm, "all %d runs have f = %g" % (len(r), boxes["f_fixed"]))
    for nm, r, bx in (("OFF", res_off, boxes["off"]), ("ON", res_on, boxes["on"])):
        out_k = r[(r["Ks_mult"] < bx["ks_lo"] * (1 - tolf)) | (r["Ks_mult"] > bx["ks_hi"] * (1 + tolf))]
        if nm == "ON":
            out_c = r[(r["channelconductivity_mmhr"] < bx["cc_lo"] * (1 - tolf)) | (r["channelconductivity_mmhr"] > bx["cc_hi"] * (1 + tolf))]
        else:
            out_c = r.iloc[0:0]
        if len(out_k) or len(out_c):
            ck.add("FAIL", "runs inside the box (%s)" % nm, "%d run(s) lie outside Ks %g-%g%s. This is not the sweep the design "
                   "file describes; the verdict would be wrong."
                   % (len(out_k), bx["ks_lo"], bx["ks_hi"], (" and %d outside cc %g-%g" % (len(out_c), bx["cc_lo"], bx["cc_hi"])) if nm == "ON" else ""))
        else:
            ck.add("PASS", "runs inside the box (%s)" % nm, "all %d runs lie inside Ks %g-%g%s"
                   % (len(r), bx["ks_lo"], bx["ks_hi"], (" and cc %g-%g" % (bx["cc_lo"], bx["cc_hi"])) if nm == "ON" else ""))

    phase_fn = RS.import_phase_fn(script_dir)
    if phase_fn is None:
        sys.exit("The scorer's phase-metric function could not be imported (run_sensitivity_single_interp_tribs6.py must be in "
                 "this folder). Without it gate 2 cannot pass.")
    obs_synth = RS.read_truth_5min(ccoff_path)

    print("\nLoading compare CSVs (run membership = the run_id rows of the results tables):")
    allres = pd.concat([res_on.assign(_series=SER_ON), res_off.assign(_series=SER_OFF)], ignore_index=True)
    cache, bad = RS.load_all_compares(allres, csv_dir)
    print("  %d/%d loaded" % (len(cache), len(allres)))

    # ---- GATE 2 -------------------------------------------------------------------------------------
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

    run_ids = [r for r in allres["run_id"] if r in cache]
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
    res_on = res_on[res_on["run_id"].isin(run_ids)].reset_index(drop=True)
    res_off = res_off[res_off["run_id"].isin(run_ids)].reset_index(drop=True)
    if res_on.empty or res_off.empty:
        sys.exit("A series has no runs left after the gates (OFF %d, ON %d); stopping." % (len(res_off), len(res_on)))

    if args.check_only:
        print("\n--check_only: gates finished, nothing written, no fit shown. A full run scores %d OFF and %d ON runs: %s."
              % (len(res_off), len(res_on),
                 "; ".join("%s arm on %d real points" % (a, len(valids[a].index.intersection(ref_idx))) for a in arms_used)))
        return

    # ---- score and analyse every arm in use -----------------------------------------------------------
    results = {}
    for arm in arms_used:
        ar = RG.score_arm(arm, valids[arm], ref_idx, res_on, res_off, cache, phase_fn)
        if ar["n_failed"]:
            print("  WARNING (%s arm): %d run(s) could not be scored" % (arm, ar["n_failed"]))
        if ar["m_on"].empty or ar["m_off"].empty or not np.isfinite(ar["m_on"]["kge_2012"]).any() \
                or not np.isfinite(ar["m_off"]["kge_2012"]).any():
            sys.exit("Nothing could be scored in one series of the %s arm." % arm)
        an = analyze_arm(ar, args, boxes)
        results[arm] = {"ar": ar, "an": an}
        print_arm(arm, arm == primary, an, len(ar["scored_idx"]), boxes, partial_note, args)

    # ---- agreement between the arms ------------------------------------------------------------------------
    ag = pd.DataFrame([{"arm": a, "primary": a == primary, "n_points": len(results[a]["ar"]["scored_idx"]),
                        "best_kge_on": results[a]["an"]["on_row"]["best_kge_2012"],
                        "best_kge_off": results[a]["an"]["off_row"]["best_kge_2012"],
                        "gap_best": results[a]["an"]["gap"], "gap_top_n_mean": results[a]["an"]["gap_mean"],
                        "row_best": results[a]["an"]["row_best"], "row_mean": results[a]["an"]["row_mean"],
                        "row_wide": results[a]["an"]["row_wide"]} for a in arms_used])
    if len(ag) == 2:
        RG.print_arm_agreement(ag, primary, args.margin)

    # ---- write ------------------------------------------------------------------------------------------------------
    out_dir.mkdir(parents=True, exist_ok=True)
    write_csv_atomic(pd.DataFrame(ck.rows), out_dir / "stageb_gates_110.csv")
    record_pts = valids["asis"].loc[valids["asis"].index.intersection(ref_idx)]
    for arm in arms_used:
        adir = out_dir / arm
        adir.mkdir(parents=True, exist_ok=True)
        ar, an = results[arm]["ar"], results[arm]["an"]
        long_df = pd.concat([an["m_on"].assign(series=SER_ON), an["m_off"].assign(series=SER_OFF)], ignore_index=True)
        write_csv_atomic(long_df, adir / "stageb_long_110.csv")
        write_csv_atomic(an["sdf"], adir / "stageb_summary_110.csv")
        keep = ["series", "rank", "run_id", "Ks_mult", "f_RS_abs", "channelconductivity_mmhr", "kge_2012", "kge_r", "kge_beta",
                "kge_gamma", "pbias_pct", "nse", "rmse_m3s", "peak_error_pct", "peak_timing_error_hr"]
        tabs = []
        for series, m in ((SER_ON, an["m_on"]), (SER_OFF, an["m_off"])):
            t = RG.top_rows(m, 25).reset_index(drop=True)
            t["rank"] = np.arange(1, len(t) + 1)
            t["series"] = series
            tabs.append(t[keep])
        write_csv_atomic(pd.concat(tabs, ignore_index=True), adir / "stageb_top_runs_110.csv")
        write_csv_atomic(an["bands"], adir / "stageb_cc_bands_110.csv")
        write_csv_atomic(pd.DataFrame([an["loc"], an["loc_wide"]]), adir / "stageb_cc_location_110.csv")
        rrow = dict(an["ridge"])
        rrow.update({"low_cc_n": an["lowcc"]["n"], "low_cc_median_diff": an["lowcc"]["median_diff"],
                     "low_cc_min_diff": an["lowcc"]["min_diff"], "low_cc_max_diff": an["lowcc"]["max_diff"]})
        write_csv_atomic(pd.DataFrame([rrow]), adir / "stageb_ridge_110.csv")
        vrow = {"arm": arm, "primary": arm == primary, "label": args.label, "f_fixed": boxes["f_fixed"], "margin": args.margin,
                "best_kge_on": an["on_row"]["best_kge_2012"], "best_kge_off": an["off_row"]["best_kge_2012"], "gap_best": an["gap"],
                "gap_top_mean": an["gap_mean"], "row_best": an["row_best"], "row_mean": an["row_mean"],
                "row_wide_location": an["row_wide"],
                "rows_agree": len(set([an["row_best"], an["row_mean"], an["row_wide"]])) == 1,
                "band_spread": an["band_spread"], "flat": an["flat"], "flat_override": an["flat_override"],
                "near_margin": an["near_margin"],
                "best_on_Ks": an["on_row"]["best_Ks"], "best_on_cc": an["on_row"]["best_cc"],
                "best_off_Ks": an["off_row"]["best_Ks"], "n_on": int(len(an["m_on"])), "n_off": int(len(an["m_off"])),
                "n_top_on": an["n_top_on"], "n_top_off": an["n_top_off"], "partial": partial,
                "p_lowest_third": an["loc"]["p_lowest_third"], "p_highest_third": an["loc"]["p_highest_third"],
                "span_frac_of_box": an["loc"]["span_frac_of_box"], "edge_notes": " | ".join(an["edge_notes"])}
        write_csv_atomic(pd.DataFrame([vrow]), adir / "stageb_verdict_110.csv")
        if args.no_plots:
            continue
        figs = [("fig_stageb_kge_vs_cc_110.png", lambda p, an=an, arm=arm: fig_kge_vs_cc(an, args.margin, p, arm, partial_note, args.label, boxes)),
                ("fig_stageb_kge_vs_Ks_110.png", lambda p, an=an, arm=arm: fig_kge_vs_ks(an, args.margin, p, arm, partial_note, args.label)),
                ("fig_stageb_ridge_110.png", lambda p, an=an, arm=arm: fig_ridge(an, p, arm, partial_note, args.label, boxes))]

        def _hydro(p, an=an, arm=arm, ar=ar):
            sim_on = cache[an["on_row"]["best_run_id"]]["Simulated"].reindex(ref_idx)
            sim_off = cache[an["off_row"]["best_run_id"]]["Simulated"].reindex(ref_idx)
            RG.fig_hydrographs(ar["real_scored"], sim_on, sim_off, an["on_row"], an["off_row"], p, arm,
                               record_pts if arm == "filled" else None)
        figs.append(("fig_stageb_hydrograph_110.png", _hydro))

        def _hydro_zoom(p, an=an, arm=arm, ar=ar):
            sim_on = cache[an["on_row"]["best_run_id"]]["Simulated"].reindex(ref_idx)
            sim_off = cache[an["off_row"]["best_run_id"]]["Simulated"].reindex(ref_idx)
            fig_hydrograph_zoom(ar["real_scored"], record_pts if arm == "filled" else None, sim_on, sim_off,
                                an["on_row"], an["off_row"], p, arm)
        figs.append(("fig_stageb_hydrograph_zoom_110.png", _hydro_zoom))
        for name, fn in figs:
            try:
                fn(adir / name)
            except Exception as exc:
                print("\n  (figure %s/%s skipped: %s -- the CSVs are saved regardless)" % (arm, name, exc))

    scorer_py = script_dir / "run_sensitivity_single_interp_tribs6.py"
    prov = {
        "script": Path(__file__).name, "version": 1, "created_local": datetime.now().isoformat(timespec="seconds"),
        "python": sys.version.split()[0], "pandas": pd.__version__, "numpy": np.__version__,
        "label": args.label, "f_fixed": boxes["f_fixed"], "margin": args.margin,
        "arms": {"used": arms_used, "primary": primary, "max_gap_min": args.max_gap_min, "max_hold_hr": args.max_hold_hr,
                 "filled_counts": finfo, "filled_checks": infof},
        "boxes": boxes, "n_runs_scored": {"ON": int(len(res_on)), "OFF": int(len(res_off))}, "partial": partial,
        "reading_rule": {"margin": args.margin, "TOP_FRAC": args.top_frac, "TOP_MIN": TOP_MIN, "NEAR_MARGIN_TOL": NEAR_MARGIN_TOL,
                         "ALPHA": RG.ALPHA, "STABLE_SPAN": RG.STABLE_SPAN, "FLOOR_TOL_MMHR": args.floor_tol,
                         "cc_bands": args.cc_bands},
        "verdicts": {a: {"row_best": int(results[a]["an"]["row_best"]), "row_mean": int(results[a]["an"]["row_mean"]),
                         "row_wide_location": int(results[a]["an"]["row_wide"]),
                         "gap_best": float(results[a]["an"]["gap"])} for a in arms_used},
        "inputs": {"on_csv": {"path": str(on_path), "md5": RS.md5_of(on_path)},
                   "off_csv": {"path": str(off_path), "md5": RS.md5_of(off_path)},
                   "on_design": {"path": str(d_on["path"]), "md5": RS.md5_of(d_on["path"])},
                   "off_design": {"path": str(d_off["path"]), "md5": RS.md5_of(d_off["path"])},
                   "gauge_xlsx": {"path": str(xlsx_path), "md5": info1.get("xlsx_md5")},
                   "synthetic_truth": {"file": ccoff_path.name, "md5": RS.md5_of(ccoff_path)},
                   "scorer_py_md5": RS.md5_of(scorer_py) if scorer_py.exists() else None,
                   "rescore_real_gauge_md5": RS.md5_of(Path(RG.__file__)),
                   "this_script_md5": RS.md5_of(Path(__file__))},
        "gates": {"2A_max_abs_diff_m3s": worstA, "2B_worst_abs_dev": worstB, "runs_checked": n_checked,
                  "runs_excluded": {k: v for k, v in bad.items()}, "3": gate3_state, "checks": ck.rows},
        "settings": {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()},
    }
    write_json_atomic(out_dir / "PROVENANCE_stageb_110.json", prov)
    print("\nSaved to: %s" % out_dir)
    print("  top folder: stageb_gates_110.csv   PROVENANCE_stageb_110.json")
    print("  per arm (%s): stageb_long_110.csv   stageb_summary_110.csv   stageb_top_runs_110.csv   stageb_cc_bands_110.csv   "
          "stageb_cc_location_110.csv   stageb_ridge_110.csv   stageb_verdict_110.csv   fig_stageb_*.png"
          % ", ".join("%s/" % a for a in arms_used))


# the docstring carries the rule constants; fill them in
if __doc__:
    __doc__ = __doc__ % {"margin": DEFAULT_MARGIN, "topfrac": int(round(TOP_FRAC * 100)), "topfrac2": int(round(2 * TOP_FRAC * 100)),
                         "alpha": RG.ALPHA, "span": RG.STABLE_SPAN, "nm": NEAR_MARGIN_TOL, "flat": FLAT_BANDS_TOL}

if __name__ == "__main__":
    main()
