"""
plot_series122_hydrographs.py
=============================
Series 122: LOOK at the hydrographs behind the numbers.  (Version 2: zooms in on the storm.)

The Series 122 analysis (analyze_series122_kchan.py) works with scores: KGE_2012 and its parts. This script draws the
outlet hydrographs themselves, to see WHEN in the storm the channel loss takes water and what a different Kchan does to
the shape. It runs NO tRIBS and changes none of your files. It only reads the kept outlet files, the two truths and the
results the analysis already wrote.

WHAT IS NEW IN VERSION 2
------------------------
The scoring window is 20 hours (12 Aug 16:00 to 13 Aug 12:00), but the action takes place in roughly 1.5 to 3 hours of
it, so version 1's figures were mostly flat line. This version draws every figure twice:
    * a ZOOM version (file names end in _zoom) on the storm itself, with a dot on every 5-minute point so you can see the
      real time resolution, and
    * the FULL 20-hour version, exactly as before (so nothing you have already looked at changes).
Nothing in the numbers changes: every volume, peak, score and ratio is still computed over the whole 20-hour scoring window,
and zooming only changes what is drawn. The one new thing is that the script now says how much of the lost volume falls
BEFORE, INSIDE and AFTER the zoom window, so a zoom can never hide where the loss happens.
Figure 3's title is also corrected: version 1 called it "what does Kchan alone do", but the three runs only share their
routing as far as the Series 121 dial can tell (cv, flowexp and n still differ a little), so it is "similar routing,
different Kchan". The cleaner test of Kchan alone is run_cc_clean_122_tribs6.py (plot_series122_clean_runs.py draws it).

THREE FIGURES (and the numbers behind them, printed and saved)
--------------------------------------------------------------
  1. Truth A (no channel loss) against Truth B (channel loss ON, Kchan 201).  The two truths use the same Ks, f, cv,
     flowexp and n, so the difference between them IS the channel loss.  Four panels:
        (a) the two hydrographs,  (b) Truth B minus Truth A through time,
        (c) the running total of the lost volume (% of Truth A's volume),
        (d) the ratio B / A at each flow level, rising limb and falling limb apart.  If the loss took the same SHARE
            of the flow everywhere, the ratio would be flat. If it takes more of the small flows, it slopes.
            (Panel d always uses every point of the whole scoring window above 5% of A's peak; it is not zoomed.)
  2. The best run against Truth B, and the best run against Truth A, each laid over its own truth, with the gap
     (run minus truth) underneath.
  3. THREE RUNS WITH SIMILAR ROUTING BUT DIFFERENT Kchan (low, middle, high), against Truth B. This is the closest the
     400 existing runs come to changing Kchan alone. It is NOT Kchan alone: cv, flowexp and n differ a little between the
     three (their values are printed under the figure).

HOW THE ZOOM WINDOW IS CHOSEN (a rule fixed in this script, computed from Truth A, not chosen by eye)
-----------------------------------------------------------------------------------------------------
  Add up Truth A's flow through the scoring window (that is its volume). The zoom starts 30 minutes before 1% of that
  volume has passed and ends 30 minutes after 97.5% has passed, rounded out to the 5-minute grid and kept inside the scoring
  window. One window is used for every figure so they can be compared side by side.
  To set it yourself:  --zoom_start "2014-08-12 17:00" --zoom_end "2014-08-12 21:30"   (either one alone also works; the
  other end then comes from the rule).

HOW THE RUNS ARE PICKED (rules fixed in this script, before looking at any hydrograph; not chosen by eye)
-------------------------------------------------------------------------------------------------------
  Best run against Truth B / against Truth A:   the run with the highest KGE_2012 against that truth.
  The three runs for figure 3:
     * "dial value" = the Series 121 dial, w1*cv + w2*flowexp + w3*n on standardized parameters (weights -0.36, +0.91,
       +0.21). Routing acts as ONE direction (Series 121 and 122), so runs with the same dial value have almost the
       same routing as far as the score is concerned.
     * candidates: runs whose dial value is within 0.5 dial units of the TRUTH's own dial value
       (cv 4.5, flowexp 0.24, n 0.026), so the hydrographs sit near the truth's and the Kchan effect is not buried
       under a routing misfit;
     * one run from each Kchan band:  low 3 to 15,  middle 120 to 330 (around the true 201),  high 600 to 1000 mm/hr;
     * of all such triples take the one whose three dial values are closest together (spreads within 0.03 count as a
       tie); ties go to the triple nearest the truth's dial value (within 0.10 counts as a tie); remaining ties go to
       the triple whose cv, flowexp and n are closest, each scaled by its box width;
     * if even the best triple differs by more than 0.15 dial units, or a band has no candidate, figure 3 is NOT drawn
       and the script says why.

WHAT IT READS  (only reads)
---------------------------
    calibration_work/03_comparisons/summary_tables/lhs_results_kchan_122.csv
    calibration_work/03_comparisons/sensitivity_plots/series122_kchan/series122_scored_both_truths.csv
                                                              (written by analyze_series122_kchan.py; run that first)
    calibration_work/kept_qout_122/<run_id>_Outlet.qout       the kept outlet files
    calibration_work/synth_truth/*.qout                       truth A
    calibration_work/synth_truth_location/SMF_20140812_110_cc201p0_loc_K7p0_f0p012_Outlet.qout    truth B

WHAT IT WRITES  (a NEW folder; it refuses to write into a folder that holds the results file)
---------------------------------------------------------------------------------------------
  calibration_work/03_comparisons/sensitivity_plots/series122_hydrographs/
    series122_h_report.txt              everything printed on screen
    series122_h_truthAB_numbers.csv     volumes, peaks, timing, share of the loss before the peak, and before / inside /
                                        after the zoom window
    series122_h_ratio_bands.csv         the B / A ratio by flow band and limb
    series122_h_truthAB_timeseries.csv  the two truth hydrographs, their difference and running totals (whole window)
    series122_h_picks.csv               the picked runs, their parameters, scores and gaps to the truth
    series122_h_picked_runs_timeseries.csv   the hydrograph of every picked run (to re-plot in a notebook)
    fig_series122_h1_truth_A_vs_B.png    fig_series122_h2_best_runs.png    fig_series122_h3_same_routing_trio.png
    fig_series122_h1_truth_A_vs_B_zoom.png   fig_series122_h2_best_runs_zoom.png   fig_series122_h3_same_routing_trio_zoom.png

HOW TO RUN  (from the lab/ folder, next to analyze_series122_kchan.py, which this script uses)
---------------------------------------------------------------------------------------------
    cd /workspaces/tRIBS-Fall2026/workspaces/SMF/lab
    python plot_series122_hydrographs.py --dry_run                    <- checks inputs, shows the picked runs and the zoom window, writes nothing
    python -u plot_series122_hydrographs.py 2>&1 | tee series122_hydrographs.log

It takes well under a minute.

BEFORE IT DRAWS ANYTHING IT CHECKS ITSELF, and stops if a check fails:
  * truth A and truth B are the exact files (checksums) the analysis used; every kept outlet file matches its recorded checksum;
  * for every picked run, its KGE_2012, r, gamma, beta and PBIAS recomputed from the outlet file match the values in the
    analysis table (to 1e-6), against both truths; and truth B scored against truth A reproduces what Series 110 recorded.

OPTIONS
-------
    --dry_run       check the inputs, show which runs would be picked and the zoom window, then stop
    --window W      zoom, full or both (default both): which versions of the figures to draw
    --zoom_start T  --zoom_end T    set the zoom window yourself (see above)
    --results PATH  the sweep's results file         --scored PATH   the analysis's series122_scored_both_truths.csv
    --out_dir PATH  where to write                    --no_plots      numbers and tables only        --dpi N   (default 150)

STATUS  (what has and has not been tested)
------------------------------------------
  NOT YET TESTED (this text is replaced once the tests have been run).

READ THIS BEFORE QUOTING ANYTHING FROM THE FIGURES
--------------------------------------------------
  * One synthetic truth pair, one storm (12 Aug 2014, 1.0x), Ks 7.0 and f 0.012 fixed, no noise. These are model
    outputs, not observations. A handful of hydrographs illustrates a pattern; it does not prove one.
  * The kept files are OUTLET only. They show what leaves the basin, not WHERE in the channel network the water is lost.
  * Figure 3's three runs share a dial value, not identical cv, flowexp and n. The script prints all three parameter sets.
  * A zoom is a choice of what to show. The zoom window follows a fixed rule and the script reports how much of the loss
    falls outside it, but a pattern seen only inside the zoom is not a claim about the rest of the 20 hours.
  * Why a pattern appears (for example, why a loss would take more of the small flows) is a hypothesis to check in the
    tRIBS code, not something these curves show.
"""

import argparse
import re
import sys
import time
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import analyze_series122_kchan as a122
except ImportError as _exc:                                                          # noqa: BLE001
    print("STOPPED: this script reuses the scoring code in analyze_series122_kchan.py and could not import it "
          f"({_exc}).\nPut both scripts in the same lab/ folder and run from there:\n"
          "    cd /workspaces/tRIBS-Fall2026/workspaces/SMF/lab")
    sys.exit(1)

DataProblem, GateFailure = a122.DataProblem, a122.GateFailure

# ======================================================================
# FILES
# ======================================================================
OUT_SUBFOLDER = "series122_hydrographs"
FILES = dict(
    report="series122_h_report.txt", truthab="series122_h_truthAB_numbers.csv", ratio="series122_h_ratio_bands.csv",
    truthab_ts="series122_h_truthAB_timeseries.csv", picks="series122_h_picks.csv",
    picks_ts="series122_h_picked_runs_timeseries.csv",
    fig1="fig_series122_h1_truth_A_vs_B.png", fig2="fig_series122_h2_best_runs.png",
    fig3="fig_series122_h3_same_routing_trio.png",
    fig1z="fig_series122_h1_truth_A_vs_B_zoom.png", fig2z="fig_series122_h2_best_runs_zoom.png",
    fig3z="fig_series122_h3_same_routing_trio_zoom.png",
)

# ======================================================================
# THE PICKING RULES: fixed before any hydrograph was looked at. Do not change these to get a different picture.
# ======================================================================
TRUE_ROUTING = {a122.P_CV: 4.5, a122.P_FE: 0.24, a122.P_N: 0.026}     # truth A and truth B
DIAL_POOL = 0.5                      # candidates: within this many dial units of the truth's dial value
KBANDS = (("low", 3.0, 15.0), ("middle", 120.0, 330.0), ("high", 600.0, 1000.0))     # mm/hr
DIAL_TIE = 0.03                      # triples whose dial spreads differ by less than this tie
DIAL_TIE_TRUTH = 0.10                # among those, triples whose mean dial distance from the truth's differs by less than this tie
MAX_TRIO_SPREAD = 0.15               # figure 3 is not drawn if the best triple's dial values differ by more than this
FLOOR_FRAC = 0.05                    # flows below 5% of the truth's peak are left out of the ratio (as in the phase metrics)
FLOW_BANDS = ((0.05, 0.20), (0.20, 0.60), (0.60, 1.00))     # as a fraction of Truth A's peak
MIN_BAND_POINTS = 5                  # a band needs this many 5-minute points to count
NEAR_PEAK_MIN = 15                   # "around the peak" means within this many minutes of the truth's peak
SCORE_COLS = ("kge_2012", "kge_r", "kge_gamma", "kge_beta", "pbias_pct")
CHECK_TOL = 1e-6

# ---- the zoom rule (also fixed before looking) ----
ZOOM_LO_FRAC = 0.01                  # the zoom starts ZOOM_PAD_MIN before this share of Truth A's volume has passed
ZOOM_HI_FRAC = 0.975                 # ... and ends ZOOM_PAD_MIN after this share has passed
ZOOM_PAD_MIN = 30
MIN_ZOOM_POINTS = 24                 # a zoom window set by hand must hold at least this many 5-minute points (2 hours)
MARK = dict(marker="o", ms=2.6, markeredgewidth=0)      # a dot on every 5-minute point (zoom figures only)

TEAL = "#12946a"     # a darker version of the Series 121 aqua; blue, teal and gold pass the palette validator together
BAND_COLOR = {"low": TEAL, "middle": a122.BLUE, "high": a122.GOLD}
BLUE, ORANGE, GOLD, INK, INK_SOFT = a122.BLUE, a122.ORANGE, a122.GOLD, a122.INK, a122.INK_SOFT


# ======================================================================
# Small helpers
# ======================================================================
def short_id(run_id):
    m = re.search(r"_d(\d+)_", str(run_id))
    return f"draw {int(m.group(1))}" if m else str(run_id)[-28:]


def fmt_t(ts):
    return "n/a" if pd.isna(ts) else ts.strftime("%d %b %H:%M")


def minutes_between(t1, t2):
    if pd.isna(t1) or pd.isna(t2):
        return float("nan")
    return (t2 - t1).total_seconds() / 60.0


def first_cross(q, thr):
    m = (q >= thr).to_numpy()
    return q.index[int(np.argmax(m))] if m.any() else pd.NaT


def time_at_fraction(cum, frac):
    tot = float(cum.iloc[-1])
    if not tot > 0:
        return pd.NaT
    ok = cum.to_numpy() >= frac * tot
    return cum.index[int(np.argmax(ok))] if ok.any() else pd.NaT


def window_frame(series, ev_start, ev_end):
    """Put several 5-minute series side by side and cut the event window, as the scorer does."""
    fr = pd.DataFrame(series).dropna().loc[ev_start:ev_end]
    if len(fr) < 10:
        raise DataProblem("The hydrographs do not overlap in the event window "
                          f"({ev_start} to {ev_end}); nothing to draw.")
    return fr


# ======================================================================
# The zoom window
# ======================================================================
def storm_zoom_window(qa, lo_frac=ZOOM_LO_FRAC, hi_frac=ZOOM_HI_FRAC, pad_min=ZOOM_PAD_MIN):
    """The rule: from pad_min before lo_frac of Truth A's volume has passed to pad_min after hi_frac has passed,
    rounded out to the 5-minute grid and kept inside the series. qa = Truth A on its 5-minute grid, scoring window only.
    Returns (start, end, info)."""
    tot = float(qa.sum())
    if not tot > 0:
        raise DataProblem("Truth A has no flow in the scoring window, so there is no storm to zoom in on.")
    cum = (qa.cumsum() / tot).to_numpy()
    t_lo = qa.index[int(np.argmax(cum >= lo_frac))]
    t_hi = qa.index[int(np.argmax(cum >= hi_frac))]
    pad = pd.Timedelta(minutes=pad_min)
    start = max((t_lo - pad).floor("5min"), qa.index[0])
    end = min((t_hi + pad).ceil("5min"), qa.index[-1])
    return start, end, dict(t_lo=t_lo, t_hi=t_hi)


def resolve_zoom(args, qa):
    """The zoom window to use: the rule, with --zoom_start / --zoom_end overriding either end. Returns ((start, end), info)."""
    start, end, info = storm_zoom_window(qa)
    info["source"] = "the rule (Truth A's volume)"
    info["rule_start"], info["rule_end"] = start, end
    manual = []
    for name, val in (("--zoom_start", args.zoom_start), ("--zoom_end", args.zoom_end)):
        if val is None:
            continue
        try:
            ts = pd.Timestamp(val)
        except (ValueError, TypeError):
            raise DataProblem(f"{name} {val!r} is not a date and time. Write it like \"2014-08-12 17:30\".")
        if pd.isna(ts):
            raise DataProblem(f"{name} {val!r} is not a date and time. Write it like \"2014-08-12 17:30\".")
        manual.append(name)
        if name == "--zoom_start":
            start = ts
        else:
            end = ts
    first, last = qa.index[0], qa.index[-1]
    if manual:
        info["source"] = " and ".join(manual) + " (set by you)" + ("" if len(manual) == 2 else "; the other end from the rule")
        if not (first <= start < end <= last):
            raise DataProblem(f"The zoom window {start} to {end} must lie inside the scoring window ({first} to {last}), "
                              "with the start before the end.")
        step = (qa.index[1] - qa.index[0]).total_seconds() / 60.0
        if (end - start).total_seconds() / 60.0 < MIN_ZOOM_POINTS * step:
            raise DataProblem(f"The zoom window {start} to {end} is shorter than {MIN_ZOOM_POINTS * step / 60:.0f} hours; "
                              "make it longer (the figures need enough points to show a shape).")
    return (start, end), info


def view_of(fr, zoom):
    """What a figure shows: the whole scoring window, or the zoom window."""
    if zoom is None:
        return dict(lo=fr.index[0], hi=fr.index[-1], zoomed=False)
    return dict(lo=zoom[0], hi=zoom[1], zoomed=True)


def nice_step_minutes(span_min, target=7):
    for step in (5, 10, 15, 30, 60, 120, 180, 360, 720):
        if span_min / step <= target:
            return step
    return 1440


def time_axis(ax, view):
    """Put sensible time ticks on the axis and cut it to the view (tick spacing follows the length of the window)."""
    import matplotlib.dates as mdates
    lo, hi = view["lo"], view["hi"]
    step = nice_step_minutes((hi - lo).total_seconds() / 60.0)
    if step < 60:
        ax.xaxis.set_major_locator(mdates.MinuteLocator(byminute=range(0, 60, step)))
    else:
        ax.xaxis.set_major_locator(mdates.HourLocator(interval=max(1, step // 60)))
    if lo.date() == hi.date():
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax.set_xlabel(f"time on {lo:%d %b %Y}", fontsize=8.5, color=INK_SOFT)
    else:
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b\n%H:%M"))
    ax.set_xlim(lo, hi)


def fit_y(ax, series, view, floor0=False, include=(), pad=0.07):
    """Set the y-axis from the data that is actually inside the view (matplotlib would use all of it).
    floor0: the bottom is zero (for discharge). include: extra values the axis must reach."""
    vals = []
    for s in series:
        v = s.loc[view["lo"]:view["hi"]].to_numpy(float)
        v = v[np.isfinite(v)]
        if len(v):
            vals.append(v)
    ext = [float(x) for x in include if x is not None and np.isfinite(x)]
    if not vals and not ext:
        return
    allv = np.concatenate(vals + [np.array(ext)]) if (vals or ext) else np.array([0.0])
    a, b = float(allv.min()), float(allv.max())
    if floor0:
        a = min(a, 0.0)
    span = (b - a) if b > a else max(abs(b), 1.0) * 0.1
    ax.set_ylim(0.0 if floor0 else a - pad * span, b + pad * span)


def title_with_view(base, view, zoom_note):
    if not view["zoomed"]:
        return base
    lo, hi = view["lo"], view["hi"]
    span = f"{lo:%H:%M} to {hi:%H:%M} on {lo:%d %b %Y}" if lo.date() == hi.date() else f"{lo:%d %b %H:%M} to {hi:%d %b %H:%M}"
    return f"{base}\nzoomed to the storm: {span}   ({zoom_note})"


# ======================================================================
# Locating and checking the inputs
# ======================================================================
def locate(args):
    root = Path.cwd().parent
    calib = root / "calibration_work"
    summ = calib / "03_comparisons" / "summary_tables"
    plots = calib / "03_comparisons" / "sensitivity_plots"
    return dict(
        calib=calib,
        results=args.results or (summ / a122.RESULTS_NAME),
        scored=args.scored or (plots / a122.OUT_SUBFOLDER / a122.FILES["scored"]),
        kept=calib / a122.KEPT_DIR_NAME,
        truth_a_dir=calib / "synth_truth",
        truth_b=calib / "synth_truth_location" / a122.TRUTH_B_NAME,
        out=args.out_dir or (plots / OUT_SUBFOLDER),
    )


def preflight(paths):
    """Read-only checks. Returns what the later steps need. Raises DataProblem on a real problem."""
    a122.sec("0. Inputs and checks (nothing is written yet)")
    df = a122.load_sweep(paths["results"])
    print(f"  results file    {Path(paths['results']).name}: {len(df)} rows")
    for c in a122.FOUR + ["kge_2012"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    bad = df[a122.FOUR + ["kge_2012"]].isna().any(axis=1)
    if bad.any():
        print(f"  NOTE  {int(bad.sum())} row(s) have a missing parameter or KGE_2012 and are left out.")
        df = df[~bad].copy()
    if df["run_id"].duplicated().any():
        raise DataProblem("A run_id appears more than once in the results file.")
    print("  ok    every run_id is unique")

    starts, ends = df["event_start"].astype(str).unique(), df["event_end"].astype(str).unique()
    if len(starts) != 1 or len(ends) != 1:
        raise DataProblem(f"The results file has more than one event window ({list(starts)} to {list(ends)}).")
    ev_start, ev_end = starts[0], ends[0]
    print(f"  ok    one event window: {ev_start} to {ev_end}")

    ta_dir = paths["truth_a_dir"]
    ta = sorted(ta_dir.glob("*.qout")) if ta_dir.exists() else []
    if len(ta) != 1:
        raise DataProblem(f"Expected exactly one *.qout file in {ta_dir} (truth A), found {len(ta)}.")
    ta_md5 = a122.md5_of(ta[0])
    if ta_md5 != a122.TRUTH_A_MD5:
        raise DataProblem(f"Truth A ({ta[0].name}) has checksum {ta_md5}, not the {a122.TRUTH_A_MD5} that Series 122 "
                          "was scored against. Do not continue with a different file.")
    print(f"  ok    truth A  {ta[0].name}  [md5 {ta_md5[:12]}]  is the file Series 122 used")
    tb = Path(paths["truth_b"])
    if not tb.exists():
        raise DataProblem(f"Truth B not found: {tb}")
    tb_md5 = a122.md5_of(tb)
    if tb_md5 != a122.TRUTH_B_MD5:
        raise DataProblem(f"Truth B has checksum {tb_md5}, not the verified {a122.TRUTH_B_MD5}. Do not continue.")
    print(f"  ok    truth B  {tb.name}  [md5 {tb_md5[:12]}]  is the file verify_truth_B_122_tribs6.py verified")

    kept = Path(paths["kept"])
    if not kept.is_dir():
        raise DataProblem(f"The kept outlet folder {kept} does not exist.")
    kept_paths, missing, md5_bad = {}, [], []
    for r in df.itertuples():
        f = kept / f"{r.run_id}_Outlet.qout"
        if not f.exists():
            missing.append(r.run_id)
            continue
        kept_paths[r.run_id] = f
        want = str(r.qout_md5)
        if want and want != "nan" and a122.md5_of(f) != want:
            md5_bad.append(r.run_id)
    if md5_bad:
        raise DataProblem(f"{len(md5_bad)} kept outlet file(s) do not match the checksum the sweep recorded (for "
                          f"example {md5_bad[0]}). The kept files must be the ones the sweep scored.")
    if missing:
        print(f"  WARN  {len(missing)} run(s) have no kept outlet file and are left out (for example {missing[0]})")
        df = df[~df["run_id"].isin(missing)].copy()
    print(f"  ok    {len(kept_paths)} kept outlet files found; each checksum equals the one the sweep recorded")

    sp = Path(paths["scored"])
    if not sp.exists():
        raise DataProblem(f"Could not find {sp}\nRun analyze_series122_kchan.py first (it writes that file), or point "
                          "to it with --scored.")
    try:
        sc = pd.read_csv(sp)
    except Exception as exc:                                                         # noqa: BLE001
        raise DataProblem(f"Could not read {sp} as a CSV: {exc}")
    need = ["run_id"] + [f"{m}_{s}" for m in SCORE_COLS for s in "AB"]
    miss = [c for c in need if c not in sc.columns]
    if miss:
        raise DataProblem("The analysis table is missing these column(s): " + ", ".join(miss)
                          + ". Is this series122_scored_both_truths.csv?")
    d = df.merge(sc[need], on="run_id", how="inner")
    if len(d) < len(df):
        print(f"  WARN  {len(df) - len(d)} run(s) are not in the analysis table and are left out")
    if len(d) < 100:
        raise DataProblem(f"Only {len(d)} runs are in both tables; need at least 100.")
    d = d.dropna(subset=["kge_2012_A", "kge_2012_B"]).reset_index(drop=True)
    d["Kchan"] = d[a122.P_K].astype(float)
    print(f"  ok    analysis table  {sp.name}: {len(d)} runs with scores against both truths")

    out_dir = Path(paths["out"])
    for guarded in (Path(paths["results"]).resolve().parent, sp.resolve().parent):
        if out_dir.resolve() == guarded:
            raise DataProblem(f"The output folder {out_dir} is a folder that holds an input table; refusing to write "
                              "there. Use the default or another --out_dir.")
    return dict(df=d, ev_start=ev_start, ev_end=ev_end, ta_path=ta[0], tb_path=tb, kept_paths=kept_paths,
                out_dir=out_dir)


# ======================================================================
# Picking the runs (rules above)
# ======================================================================
def dial_values(d):
    X = d[a122.ROUTING].to_numpy(float)
    Z, mu, sd = a122.standardize_cols(X)
    dial = Z @ a122.W121
    t = np.array([TRUE_ROUTING[p] for p in a122.ROUTING], float)
    return dial, float(((t - mu) / sd) @ a122.W121), float(np.std(dial, ddof=1))


def choose_trio(d, dial, d_true):
    K = d["Kchan"].to_numpy(float)
    pool = np.abs(dial - d_true) <= DIAL_POOL
    info = dict(pool_size=int(pool.sum()), ok=False, reason="", rows=[], bands=[])
    idx_by = []
    for name, lo, hi in KBANDS:
        idx = np.where(pool & (K >= lo) & (K <= hi))[0]
        info["bands"].append((name, lo, hi, int(len(idx))))
        idx_by.append(idx)
    empty = [b for b, i in zip(KBANDS, idx_by) if len(i) == 0]
    if empty:
        info["reason"] = ("no candidate run has Kchan in the " + " or the ".join(f"{b[0]} band ({b[1]:g} to {b[2]:g} mm/hr)"
                                                                                 for b in empty))
        return info
    R = np.column_stack([(d[p].to_numpy(float) - a122.BOX[p][0]) / (a122.BOX[p][1] - a122.BOX[p][0])
                         for p in a122.ROUTING])
    ia, ib, ic = idx_by
    da, db, dc = dial[ia], dial[ib], dial[ic]
    hi_ = np.maximum(np.maximum(da[:, None, None], db[None, :, None]), dc[None, None, :])
    lo_ = np.minimum(np.minimum(da[:, None, None], db[None, :, None]), dc[None, None, :])
    spread = hi_ - lo_

    def dist(x, y):
        return np.sqrt(((R[x][:, None, :] - R[y][None, :, :]) ** 2).sum(-1))

    dab, dac, dbc = dist(ia, ib), dist(ia, ic), dist(ib, ic)
    raw = np.maximum(np.maximum(dab[:, :, None], dac[:, None, :]), dbc[None, :, :])
    key1 = np.round(spread / DIAL_TIE).astype(int)
    mean_d = (da[:, None, None] + db[None, :, None] + dc[None, None, :]) / 3.0
    key2 = np.round(np.abs(mean_d - d_true) / DIAL_TIE_TRUTH).astype(int)
    order = np.lexsort((raw.ravel(), key2.ravel(), key1.ravel()))
    i, j, k = np.unravel_index(order[0], spread.shape)
    info["rows"] = [int(ia[i]), int(ib[j]), int(ic[k])]
    info["spread_dial"] = float(spread[i, j, k])
    info["spread_raw"] = float(raw[i, j, k])
    info["ok"] = info["spread_dial"] <= MAX_TRIO_SPREAD
    if not info["ok"]:
        info["reason"] = (f"even the best triple of runs differs by {info['spread_dial']:.3f} dial units "
                          f"(the limit is {MAX_TRIO_SPREAD}), so the three do not share their routing closely enough")
    return info


def choose_runs(d):
    a122.sec("1. Which runs the rules pick (fixed in the script; not chosen by eye)")
    dial, d_true, sd_dial = dial_values(d)
    d["dial"] = dial
    ib_, ia_ = int(d["kge_2012_B"].idxmax()), int(d["kge_2012_A"].idxmax())
    trio = choose_trio(d, dial, d_true)
    print(f"  dial value of the truth's routing (cv 4.5, flowexp 0.24, n 0.026): {d_true:+.3f}  "
          f"(the 400 runs span {dial.min():+.2f} to {dial.max():+.2f}; standard deviation {sd_dial:.2f})")
    print(f"  best run against truth B : {short_id(d.loc[ib_, 'run_id'])}  Kchan {a122.fk(d.loc[ib_, 'Kchan'])} mm/hr  "
          f"KGE_2012 {d.loc[ib_, 'kge_2012_B']:.4f}")
    print(f"  best run against truth A : {short_id(d.loc[ia_, 'run_id'])}  Kchan {a122.fk(d.loc[ia_, 'Kchan'])} mm/hr  "
          f"KGE_2012 {d.loc[ia_, 'kge_2012_A']:.4f}")
    print(f"  similar-routing trio: candidates are the {trio['pool_size']} runs within {DIAL_POOL} dial units of the truth; "
          "runs available per Kchan band: "
          + ", ".join(f"{n} ({lo:g} to {hi:g}): {c}" for n, lo, hi, c in trio["bands"]))
    if trio["rows"]:
        print(f"    best triple: dial values differ by {trio['spread_dial']:.3f} dial units "
              f"(limit {MAX_TRIO_SPREAD}); largest gap in cv, flowexp, n (each as a fraction of its box): "
              f"{trio['spread_raw']:.2f}")
        print(f"    {'band':<8}{'run':<14}{'Kchan':>9}{'cv':>8}{'flowexp':>9}{'n':>8}{'dial':>8}{'KGE vs B':>10}{'r vs B':>8}")
        for (name, _, _), row in zip(KBANDS, trio["rows"]):
            r = d.loc[row]
            print(f"    {name:<8}{short_id(r['run_id']):<14}{a122.fk(r['Kchan']):>9}{r[a122.P_CV]:>8.2f}"
                  f"{r[a122.P_FE]:>9.3f}{r[a122.P_N]:>8.3f}{r['dial']:>+8.3f}{r['kge_2012_B']:>10.3f}{r['kge_r_B']:>8.3f}")
    if not trio["ok"]:
        print(f"  NOTE  figure 3 will NOT be drawn: {trio['reason']}.")
    return dict(best_B=ib_, best_A=ia_, trio=trio, d_true=d_true)


def print_zoom(zoom, zinfo, qa, window_mode):
    a122.sec("2. The zoom window  (a rule on Truth A's volume; the numbers elsewhere still use the whole scoring window)")
    if zoom is None:
        print(f"  --window full: only the whole {len(qa) * 5 / 60:.0f}-hour scoring window is drawn; no zoom.")
        return
    lo, hi = zoom
    hours = (hi - lo).total_seconds() / 3600
    print(f"  Rule: from {ZOOM_PAD_MIN} min before {ZOOM_LO_FRAC * 100:g}% of Truth A's volume has passed "
          f"({fmt_t(zinfo['t_lo'])}) to {ZOOM_PAD_MIN} min after {ZOOM_HI_FRAC * 100:g}% has passed ({fmt_t(zinfo['t_hi'])}).")
    print(f"  Zoom window used: {fmt_t(lo)} to {fmt_t(hi)}  ({hours:.1f} hours, {int(round(hours * 12)) + 1} five-minute points), "
          f"from {zinfo['source']}.")
    if zinfo["source"] != "the rule (Truth A's volume)":
        print(f"  (the rule alone would have given {fmt_t(zinfo['rule_start'])} to {fmt_t(zinfo['rule_end'])})")
    full_h = (qa.index[-1] - qa.index[0]).total_seconds() / 3600
    print(f"  That is {100 * hours / full_h:.0f}% of the {full_h:.0f}-hour scoring window ({fmt_t(qa.index[0])} to {fmt_t(qa.index[-1])}).")
    vol_in = 100 * float(qa.loc[lo:hi].sum()) / float(qa.sum())
    print(f"  {vol_in:.1f}% of Truth A's volume passes inside the zoom window.")
    print(f"  Figures drawn: {'the zoom version and the full-window version' if window_mode == 'both' else 'the zoom version only'}.")


# ======================================================================
# Self-check (stops the script if its own scoring does not reproduce the analysis)
# ======================================================================
def self_check(pre, picks, gA, gB, runs):
    a122.sec("3. Self-check: do the outlet files give back the scores in the analysis table?")
    d = pre["df"].set_index("run_id")
    ids = list(runs)
    fails, shown = [], []
    print(f"  {'run':<30}{'largest difference over 10 checks':>44}   result   (5 scores x 2 truths)")
    for rid in ids:
        worst, bad = 0.0, False
        for suf, truth in (("A", gA), ("B", gB)):
            s = a122.score_pair(truth, runs[rid], pre["ev_start"], pre["ev_end"])
            if s is None:
                raise DataProblem(f"Run {rid} does not overlap the event window.")
            for m in SCORE_COLS:
                theirs = float(d.loc[rid, f"{m}_{suf}"])
                diff = abs(s[m] - theirs)
                worst = max(worst, diff)
                if not diff <= CHECK_TOL * (1 + abs(theirs)):
                    bad = True
                    fails.append(f"{short_id(rid)} {m} against truth {suf}: here {s[m]:.8f}, analysis table {theirs:.8f}")
        print(f"  {short_id(rid):<30}{worst:>44.3g}   {'FAIL' if bad else 'ok'}")
        shown.append(worst)
    bva = a122.score_pair(gA, gB, pre["ev_start"], pre["ev_end"])
    for m, rec in a122.REC_B_VS_A.items():
        ok = abs(bva[m] - rec) <= a122.REC_TOL * (1 + abs(rec))
        print(f"  truth B vs truth A  {m:<12} here {bva[m]:>14.8f}   Series 110 recorded {rec:>14.8f}   {'ok' if ok else 'FAIL'}")
        if not ok:
            fails.append(f"truth B vs A {m}: here {bva[m]:.8f}, Series 110 recorded {rec:.8f}")
    if fails:
        raise GateFailure("the outlet files do not reproduce the scores already on record, so nothing drawn from them "
                          "can be trusted yet:\n    - " + "\n    - ".join(fails)
                          + "\n  Nothing was drawn or written. Paste this section to Claude.")
    print("  PASS  the picked runs and the two truths reproduce the analysis table and Series 110's record.")


# ======================================================================
# Numbers
# ======================================================================
def truth_ab_analysis(fr, zoom):
    qA, qB = fr["A"], fr["B"]
    dt = (fr.index[1] - fr.index[0]).total_seconds()
    volA, volB = float(qA.sum() * dt), float(qB.sum() * dt)
    cumA, cumB = qA.cumsum() * dt, qB.cumsum() * dt
    lost = cumA - cumB
    total_loss = volA - volB
    peakA, peakB = float(qA.max()), float(qB.max())
    tA, tB = qA.idxmax(), qB.idxmax()
    num = dict(
        dt_seconds=dt, volume_A_m3=volA, volume_B_m3=volB, volume_lost_m3=total_loss,
        volume_lost_pct_of_A=100 * total_loss / volA, beta_B_over_A=volB / volA,
        peak_A_m3s=peakA, peak_B_m3s=peakB, peak_B_pct_of_A=100 * peakB / peakA,
        peak_time_A=tA, peak_time_B=tB, peak_shift_min=minutes_between(tA, tB),
        loss_share_before_A_peak_pct=(100 * float(lost.loc[tA]) / total_loss) if total_loss > 0 else float("nan"),
    )
    nan = float("nan")
    num.update(zoom_start=None, zoom_end=None, loss_share_before_zoom_pct=nan, loss_share_inside_zoom_pct=nan,
               loss_share_after_zoom_pct=nan)
    if zoom is not None:
        lo, hi = zoom

        def lost_at(t):
            s = lost.loc[:t]
            return float(s.iloc[-1]) if len(s) else 0.0

        before, at_end = lost_at(lo - pd.Timedelta(seconds=dt)), lost_at(hi)
        if total_loss > 0:
            num.update(loss_share_before_zoom_pct=100 * before / total_loss,
                       loss_share_inside_zoom_pct=100 * (at_end - before) / total_loss,
                       loss_share_after_zoom_pct=100 * (total_loss - at_end) / total_loss)
        num.update(zoom_start=lo, zoom_end=hi)
    for frac in (0.10, 0.50, 0.90):
        ta_, tb_ = time_at_fraction(cumA, frac), time_at_fraction(cumB, frac)
        num[f"t{int(frac * 100)}_A"], num[f"t{int(frac * 100)}_B"] = ta_, tb_
        num[f"t{int(frac * 100)}_shift_min"] = minutes_between(ta_, tb_)
    thr = FLOOR_FRAC * peakA
    fa, fb = first_cross(qA, thr), first_cross(qB, thr)
    num.update(first_flow_A=fa, first_flow_B=fb, first_flow_shift_min=minutes_between(fa, fb))

    mask = (qA >= thr).to_numpy()
    sub = fr[mask]
    ratio = sub["B"] / sub["A"]
    frac_of_peak = sub["A"] / peakA
    rising = (sub.index <= tA)
    pts = pd.DataFrame({"flow_A_m3s": sub["A"].to_numpy(), "ratio_B_over_A": ratio.to_numpy(),
                        "limb": np.where(rising, "rising", "falling"), "frac_of_peak": frac_of_peak.to_numpy()},
                       index=sub.index)
    rows = []
    for limb in ("rising", "falling"):
        for lo_f, hi_f in FLOW_BANDS:
            sel = (pts["limb"] == limb) & (pts["frac_of_peak"] >= lo_f) & ((pts["frac_of_peak"] < hi_f) | (hi_f >= 1.0))
            r = pts.loc[sel, "ratio_B_over_A"]
            rows.append(dict(limb=limb, band=f"{int(lo_f * 100)} to {int(hi_f * 100)}% of A's peak", n_points=int(len(r)),
                             median_ratio=float(r.median()) if len(r) else float("nan"),
                             min_ratio=float(r.min()) if len(r) else float("nan"),
                             max_ratio=float(r.max()) if len(r) else float("nan")))
    ratio_tab = pd.DataFrame(rows)
    ok = ratio_tab[ratio_tab["n_points"] >= MIN_BAND_POINTS]["median_ratio"]
    num["band_median_range"] = float(ok.max() - ok.min()) if len(ok) >= 2 else float("nan")
    num["n_bands_used"] = int(len(ok))

    ts = pd.DataFrame({"time": fr.index, "truth_A_m3s": qA.to_numpy(), "truth_B_m3s": qB.to_numpy(),
                       "B_minus_A_m3s": (qB - qA).to_numpy(), "cum_volume_A_m3": cumA.to_numpy(),
                       "cum_volume_B_m3": cumB.to_numpy(),
                       "cum_loss_pct_of_A_total": (100 * lost / volA).to_numpy()})
    return num, ratio_tab, pts, ts, lost


def compare_numbers(fr, truth_col, run_col, zoom=None):
    truth, run = fr[truth_col], fr[run_col]
    res = run - truth
    tp, rp = truth.idxmax(), run.idxmax()
    sq = res ** 2
    tot = float(sq.sum())
    iw = res.abs().idxmax()
    near = pd.Timedelta(minutes=NEAR_PEAK_MIN)
    phase = "rising limb" if iw < tp - near else ("falling limb" if iw > tp + near else "around the peak")
    thr = FLOOR_FRAC * float(truth.max())
    ft, fr_ = first_cross(truth, thr), first_cross(run, thr)
    inside = float("nan")
    if zoom is not None and tot > 0:
        inside = 100 * float(sq.loc[zoom[0]:zoom[1]].sum()) / tot
    return dict(
        volume_vs_truth_pct=100 * (float(run.sum()) / float(truth.sum()) - 1),
        peak_vs_truth_pct=100 * (float(run.max()) / float(truth.max()) - 1),
        peak_shift_min=minutes_between(tp, rp), first_flow_shift_min=minutes_between(ft, fr_),
        misfit_share_before_truth_peak_pct=(100 * float(sq.loc[:tp].sum()) / tot) if tot > 0 else float("nan"),
        misfit_share_inside_zoom_pct=inside,
        largest_gap_m3s=float(res.loc[iw]), largest_gap_time=iw, largest_gap_phase=phase,
        largest_gap_pct_of_truth_peak=100 * abs(float(res.loc[iw])) / float(truth.max()),
    )


# ======================================================================
# Printing
# ======================================================================
def print_truth_ab(num, ratio_tab):
    a122.sec("4. Truth A against Truth B: what does the channel loss do to the hydrograph?  (figure 1)")
    k = 1e3
    print("  Truth A has no channel loss. Truth B has it on (Kchan 201 mm/hr). Everything else is the same, so every "
          "difference below is the channel loss.")
    print(f"  Volume in the storm window:   A {num['volume_A_m3'] / k:>9,.1f}   B {num['volume_B_m3'] / k:>9,.1f}   "
          f"(thousand m3)    lost {num['volume_lost_m3'] / k:,.1f} = {num['volume_lost_pct_of_A']:.2f}% of A")
    print(f"  Peak discharge (m3/s):        A {num['peak_A_m3s']:>9.3f} at {fmt_t(num['peak_time_A'])}   "
          f"B {num['peak_B_m3s']:>9.3f} at {fmt_t(num['peak_time_B'])}   "
          f"(B's peak is {num['peak_B_pct_of_A']:.1f}% of A's; {abs(num['peak_shift_min']):.1f} min "
          f"{'later' if num['peak_shift_min'] > 0 else 'earlier' if num['peak_shift_min'] < 0 else 'at the same time'})")
    print(f"  First flow above {int(FLOOR_FRAC * 100)}% of A's peak:  A {fmt_t(num['first_flow_A'])}   "
          f"B {fmt_t(num['first_flow_B'])}   (B is {num['first_flow_shift_min']:+.1f} min from A)")
    for p in (10, 50, 90):
        print(f"  {p}% of the volume has passed by:  A {fmt_t(num[f't{p}_A'])}   B {fmt_t(num[f't{p}_B'])}   "
              f"(B is {num[f't{p}_shift_min']:+.1f} min from A)")
    print(f"  Share of the lost volume that is already gone by the time of A's peak: "
          f"{num['loss_share_before_A_peak_pct']:.1f}%  (the rest is lost after the peak)")
    if num["zoom_start"] is not None:
        print(f"  Where the lost volume falls relative to the zoom window ({fmt_t(num['zoom_start'])} to {fmt_t(num['zoom_end'])}): "
              f"{num['loss_share_before_zoom_pct']:.1f}% before it, {num['loss_share_inside_zoom_pct']:.1f}% inside it, "
              f"{num['loss_share_after_zoom_pct']:.1f}% after it.")
        if np.isfinite(num["loss_share_after_zoom_pct"]) and num["loss_share_after_zoom_pct"] >= 5:
            print(f"    NOTE  {num['loss_share_after_zoom_pct']:.0f}% of the loss falls after the zoom window ends. The zoom figures "
                  "show the running total still climbing at their right edge, with a dotted line at the full loss; "
                  "the full-window figures show where it levels off.")
    print("\n  Ratio B / A at each flow level (A's flow above 5% of its peak, whole scoring window). "
          "A ratio of 0.90 means B carries 90% of what A carries:")
    print(f"    {'limb':<9}{'flow band':<24}{'points':>8}{'median':>9}{'lowest':>9}{'highest':>9}")
    for r in ratio_tab.itertuples():
        print(f"    {r.limb:<9}{r.band:<24}{r.n_points:>8}{a122.fnum(r.median_ratio, 9, 3)}"
              f"{a122.fnum(r.min_ratio, 9, 3)}{a122.fnum(r.max_ratio, 9, 3)}")
    if np.isfinite(num["band_median_range"]):
        print(f"    whole-event ratio (beta) {num['beta_B_over_A']:.3f};  the band medians span {num['band_median_range']:.3f} "
              f"({num['n_bands_used']} bands with at least {MIN_BAND_POINTS} points)")
    print("  HOW TO READ THE RATIO (a key, not a verdict):")
    print("    * about the same in every band  -> the loss takes a similar SHARE of the flow at every flow level;")
    print("    * lower at small flows, nearer 1 at large flows  -> the loss behaves more like a fixed amount taken from the "
          "channel, so it removes a bigger share of the small flows;")
    print("    * lower at large flows  -> the loss grows with the flow;")
    print("    * rising limb different from falling limb  -> the loss depends on how long the water has been flowing.")
    print("    Why any of these happens is a hypothesis to check in the tRIBS code; these numbers only show that it does.")


def print_compare(title, role_rows):
    print(f"  {title}")
    for lab, run_id, K, kge, comps, c in role_rows:
        print(f"    {lab}:  {short_id(run_id)},  Kchan {a122.fk(K)} mm/hr,  KGE_2012 {kge:.4f}  "
              f"(r {comps[0]:.3f}, gamma {comps[1]:.3f}, beta {comps[2]:.3f})")
        print(f"        volume {c['volume_vs_truth_pct']:+.2f}% and peak {c['peak_vs_truth_pct']:+.2f}% of the truth's;  "
              f"peak {c['peak_shift_min']:+.1f} min and first flow {c['first_flow_shift_min']:+.1f} min from the truth's")
        zin = (f", {c['misfit_share_inside_zoom_pct']:.0f}% of it inside the zoom window"
               if np.isfinite(c.get("misfit_share_inside_zoom_pct", float("nan"))) else "")
        print(f"        largest gap {c['largest_gap_m3s']:+.3f} m3/s ({c['largest_gap_pct_of_truth_peak']:.1f}% of the truth's peak) "
              f"at {fmt_t(c['largest_gap_time'])}, {c['largest_gap_phase']};  "
              f"{c['misfit_share_before_truth_peak_pct']:.0f}% of the squared misfit comes before the truth's peak{zin}")


# ======================================================================
# Figures
# ======================================================================
def figure_1(fr, num, pts, lost, out_path, dpi, view, zoom_note):
    plt = a122._plt()
    from matplotlib.lines import Line2D
    mk = MARK if view["zoomed"] else {}
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 9.0 if view["zoomed"] else 8.8))
    fig.patch.set_facecolor("white")
    for ax in axes.ravel():
        a122.style_axes(ax)
    tA = num["peak_time_A"]
    volA = num["volume_A_m3"]
    cum_pct = 100 * lost / volA
    total_pct = 100 * num["volume_lost_m3"] / volA
    diff = fr["B"] - fr["A"]

    ax = axes[0, 0]
    ax.plot(fr.index, fr["A"], color=ORANGE, lw=2.0, label="Truth A (no channel loss)", **mk)
    ax.plot(fr.index, fr["B"], color=BLUE, lw=2.0, label="Truth B (channel loss on, Kchan 201)", **mk)
    ax.axvline(tA, color=INK_SOFT, linestyle=(0, (2, 3)), lw=1.0)
    ax.set_ylabel("discharge (m3/s)", fontsize=9.5, color=INK_SOFT)
    ax.set_title("(a) The two hydrographs", fontsize=10.5, color=INK, loc="left")
    time_axis(ax, view)
    fit_y(ax, [fr["A"], fr["B"]], view, floor0=True, pad=0.18)
    ax.legend(loc="upper right", frameon=False, fontsize=8.8)

    ax = axes[0, 1]
    ax.fill_between(fr.index, diff, 0, color=BLUE, alpha=0.22, linewidth=0)
    ax.plot(fr.index, diff, color=BLUE, lw=1.8, **mk)
    ax.axhline(0, color=INK_SOFT, lw=1.0)
    ax.axvline(tA, color=INK_SOFT, linestyle=(0, (2, 3)), lw=1.0)
    ax.set_ylabel("Truth B minus Truth A (m3/s)", fontsize=9.5, color=INK_SOFT)
    ax.set_title("(b) When is water missing?  (below zero = B has less)", fontsize=10.5, color=INK, loc="left")
    time_axis(ax, view)
    fit_y(ax, [diff], view, include=[0.0])

    ax = axes[1, 0]
    ax.plot(fr.index, cum_pct, color=BLUE, lw=2.0, **mk)
    ax.axvline(tA, color=INK_SOFT, linestyle=(0, (2, 3)), lw=1.0)
    time_axis(ax, view)
    if view["zoomed"]:
        ax.axhline(total_pct, color=INK_SOFT, linestyle=(0, (1, 2)), lw=1.1)
        fit_y(ax, [cum_pct], view, include=[0.0, total_pct])
        ax.text(view["lo"], total_pct, f" full loss over the 20 hours: {total_pct:.1f}%", fontsize=8.5, color=INK_SOFT,
                va="bottom", ha="left")
        if np.isfinite(num["loss_share_after_zoom_pct"]):
            ax.text(0.98, 0.20, f"{num['loss_share_after_zoom_pct']:.0f}% of the loss comes after this window ends",
                    transform=ax.transAxes, fontsize=9, color=INK, va="bottom", ha="right")
    else:
        fit_y(ax, [cum_pct], view, include=[0.0])
    ax.text(tA, ax.get_ylim()[1], "  A's peak", fontsize=8.5, color=INK_SOFT, va="top", ha="left")
    ax.text(0.98, 0.08, f"{num['loss_share_before_A_peak_pct']:.0f}% of the lost volume is gone by A's peak",
            transform=ax.transAxes, fontsize=9, color=INK, va="bottom", ha="right")
    ax.set_ylabel("running total lost (% of Truth A's volume)", fontsize=9.5, color=INK_SOFT)
    ax.set_title("(c) How the lost volume builds up", fontsize=10.5, color=INK, loc="left")

    ax = axes[1, 1]
    for limb, col in (("rising", BLUE), ("falling", GOLD)):
        s = pts[pts["limb"] == limb]
        ax.scatter(s["flow_A_m3s"], s["ratio_B_over_A"], s=16, color=col, alpha=0.85, linewidths=0, zorder=3)
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.axhline(num["beta_B_over_A"], color=INK_SOFT, linestyle=(0, (2, 3)), lw=1.0)
    if len(pts):
        lo = min(float(pts["ratio_B_over_A"].min()), num["beta_B_over_A"])
        hi = max(float(pts["ratio_B_over_A"].max()), num["beta_B_over_A"])
        mid, half = (lo + hi) / 2.0, max((hi - lo) / 2.0 * 1.08, 0.05)        # never zoom in closer than +-0.05
        ax.set_ylim(mid - half, mid + half)
        ax.set_xlim(float(pts["flow_A_m3s"].min()) * 0.9, float(pts["flow_A_m3s"].max()) * 1.1)
    ax.ticklabel_format(axis="y", useOffset=False)
    ax.set_xlabel("Truth A's discharge (m3/s, log scale)", fontsize=9.5, color=INK_SOFT)
    ax.set_ylabel("ratio  B / A", fontsize=9.5, color=INK_SOFT)
    ax.set_title("(d) Is the loss the same share of the flow at every flow level?\n(whole 20-hour window, not zoomed)",
                 fontsize=10.5, color=INK, loc="left")
    ax.legend(handles=[Line2D([0], [0], marker="o", linestyle="", color=BLUE, markersize=6, label="rising limb"),
                       Line2D([0], [0], marker="o", linestyle="", color=GOLD, markersize=6, label="falling limb"),
                       Line2D([0], [0], color=INK_SOFT, linestyle=(0, (2, 3)), lw=1.2,
                              label=f"whole event: {num['beta_B_over_A']:.3f}")],
              loc="lower right", frameon=False, fontsize=8.8)
    fig.suptitle(title_with_view("Series 122: what does a channel loss do to the hydrograph?  (Truth A against Truth B)",
                                 view, zoom_note), fontsize=12.5, color=INK, y=0.995)
    fig.text(0.5, 0.005, "Truth A and Truth B share Ks 7.0, f 0.012, cv 4.5, flowexp 0.24 and n 0.026; only the channel loss differs. "
             "One storm (12 Aug 2014). Model output, not observations. Dotted vertical line: Truth A's peak."
             + ("  A dot is one 5-minute point." if view["zoomed"] else ""),
             ha="center", va="bottom", fontsize=8.2, color=INK_SOFT)
    fig.tight_layout(rect=(0, 0.03, 1, 0.93 if view["zoomed"] else 0.96))
    a122.save_fig(fig, out_path, dpi)


def figure_2(frames, info, out_path, dpi, view, zoom_note):
    plt = a122._plt()
    mk = MARK if view["zoomed"] else {}
    fig, axes = plt.subplots(2, 2, figsize=(13.4, 8.8), gridspec_kw=dict(height_ratios=[3, 1.5]))
    fig.patch.set_facecolor("white")
    for j, (suf, tcol, key) in enumerate((("B", BLUE, "best_B"), ("A", ORANGE, "best_A"))):
        fr, row = frames[key], info[key]
        ax, axr = axes[0, j], axes[1, j]
        a122.style_axes(ax)
        a122.style_axes(axr)
        ax.plot(fr.index, fr["truth"], color=tcol, lw=5.5, alpha=0.40, solid_capstyle="round",
                label=f"Truth {suf}" + (" (Kchan 201)" if suf == "B" else " (no loss)"))
        ax.plot(fr.index, fr["run"], color=INK, lw=1.5,
                label=f"best run ({short_id(row['run_id'])}, Kchan {a122.fk(row['Kchan'])} mm/hr)", **mk)
        ax.set_ylabel("discharge (m3/s)", fontsize=9.5, color=INK_SOFT)
        ax.set_title(f"Best run against Truth {suf}:  KGE_2012 {row[f'kge_2012_{suf}']:.3f}\n"
                     f"r {row[f'kge_r_{suf}']:.3f}   gamma {row[f'kge_gamma_{suf}']:.3f}   beta {row[f'kge_beta_{suf}']:.3f}"
                     "   (scores use the whole 20 hours)",
                     fontsize=10, color=INK, loc="left")
        time_axis(ax, view)
        fit_y(ax, [fr["truth"], fr["run"]], view, floor0=True, pad=0.18)
        ax.legend(loc="upper right", frameon=False, fontsize=8.8)
        res = fr["run"] - fr["truth"]
        axr.fill_between(fr.index, res, 0, color=INK_SOFT, alpha=0.22, linewidth=0)
        axr.plot(fr.index, res, color=INK, lw=1.3, **mk)
        axr.axhline(0, color=INK_SOFT, lw=1.0)
        axr.set_ylabel("run minus truth (m3/s)", fontsize=9.5, color=INK_SOFT)
        time_axis(axr, view)
        fit_y(axr, [res], view, include=[0.0])
    fig.suptitle(title_with_view("Series 122: the best run against each truth", view, zoom_note),
                 fontsize=12.5, color=INK, y=0.995)
    fig.text(0.5, 0.005, "Thick pale line = the truth; thin dark line = the run. Same Ks 7.0 and f 0.012 in every run. "
             "One storm (12 Aug 2014). Model output, not observations." + ("  A dot is one 5-minute point." if view["zoomed"] else ""),
             ha="center", va="bottom", fontsize=8.2, color=INK_SOFT)
    fig.tight_layout(rect=(0, 0.03, 1, 0.93 if view["zoomed"] else 0.96))
    a122.save_fig(fig, out_path, dpi)


def figure_3(fr, trio_rows, num, spread, out_path, dpi, view, zoom_note):
    plt = a122._plt()
    mk = MARK if view["zoomed"] else {}
    fig = plt.figure(figsize=(13.4, 9.4))
    fig.patch.set_facecolor("white")
    gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1.0])
    ax1, ax2, ax3 = fig.add_subplot(gs[0, :]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
    for ax in (ax1, ax2, ax3):
        a122.style_axes(ax)
    ax1.plot(fr.index, fr["B"], color=INK_SOFT, lw=6.0, alpha=0.40, solid_capstyle="round", label="Truth B (Kchan 201)")
    ax1.plot(fr.index, fr["A"], color=INK, lw=1.5, linestyle=(0, (3, 2)), label="Truth A (no loss)")
    volB = float(fr["B"].sum())
    ax2.axhline(0, color=INK_SOFT, lw=1.0)
    ax3.axhline(0, color=INK_SOFT, lw=1.0)
    dt = (fr.index[1] - fr.index[0]).total_seconds()
    cumB = fr["B"].cumsum() * dt
    gapA = 100 * (fr["A"].cumsum() * dt - cumB) / (volB * dt)
    ax3.plot(fr.index, gapA, color=INK, lw=1.5, linestyle=(0, (3, 2)))
    shown1, shown2, shown3 = [fr["A"], fr["B"]], [], [gapA]
    for name, r in trio_rows:
        col = BAND_COLOR[name]
        lab = f"{name} Kchan {a122.fk(r['Kchan'])} mm/hr  (KGE_2012 vs B {r['kge_2012_B']:.3f})"
        ax1.plot(fr.index, fr[name], color=col, lw=1.9, label=lab, **mk)
        g2 = fr[name] - fr["B"]
        g3 = 100 * (fr[name].cumsum() * dt - cumB) / (volB * dt)
        ax2.plot(fr.index, g2, color=col, lw=1.8, **mk)
        ax3.plot(fr.index, g3, color=col, lw=1.8, **mk)
        shown1.append(fr[name])
        shown2.append(g2)
        shown3.append(g3)
    ax1.set_ylabel("discharge (m3/s)", fontsize=9.5, color=INK_SOFT)
    ax1.set_title("Three runs with similar routing (not identical) and different Kchan, against Truth B",
                  fontsize=10.5, color=INK, loc="left")
    ax2.set_ylabel("run minus Truth B (m3/s)", fontsize=9.5, color=INK_SOFT)
    ax2.set_title("Gap to Truth B", fontsize=10.5, color=INK, loc="left")
    ax3.set_ylabel("running volume gap to Truth B (% of B's volume)", fontsize=9.5, color=INK_SOFT)
    ax3.set_title("Running volume gap to Truth B (dashed: Truth A)", fontsize=10.5, color=INK, loc="left")
    for ax in (ax1, ax2, ax3):
        time_axis(ax, view)
    fit_y(ax1, shown1, view, floor0=True, pad=0.20)
    fit_y(ax2, shown2, view, include=[0.0])
    fit_y(ax3, shown3, view, include=[0.0])
    ax1.legend(loc="upper right", frameon=False, fontsize=8.8)
    fig.suptitle(title_with_view("Series 122: three runs with similar routing and different Kchan", view, zoom_note),
                 fontsize=12.5, color=INK, y=0.995)
    rt = "   |   ".join(f"{n}: cv {r[a122.P_CV]:.2f}, flowexp {r[a122.P_FE]:.3f}, n {r[a122.P_N]:.3f}" for n, r in trio_rows)
    fig.text(0.5, 0.005, f"This is NOT Kchan alone: the routing is matched on the Series 121 dial (the three dial values differ by {spread:.3f}), "
             f"so cv, flowexp and n still differ a little:\n{rt}.\n"
             "One storm (12 Aug 2014). Model output, not observations." + ("  A dot is one 5-minute point." if view["zoomed"] else ""),
             ha="center", va="bottom", fontsize=8.2, color=INK_SOFT, linespacing=1.5)
    fig.tight_layout(rect=(0, 0.065, 1, 0.93 if view["zoomed"] else 0.96))
    a122.save_fig(fig, out_path, dpi)


# ======================================================================
# Plumbing
# ======================================================================
def parse_args(argv=None):
    ap = argparse.ArgumentParser(description="Series 122: draw the hydrographs behind the Kchan results, zoomed on the storm "
                                             "(reads results; writes a new folder).")
    ap.add_argument("--dry_run", action="store_true")
    ap.add_argument("--window", choices=("zoom", "full", "both"), default="both",
                    help="which versions of the figures to draw (default both)")
    ap.add_argument("--zoom_start", default=None, help="start of the zoom window, e.g. \"2014-08-12 17:00\" (default: the rule)")
    ap.add_argument("--zoom_end", default=None, help="end of the zoom window, e.g. \"2014-08-12 21:30\" (default: the rule)")
    ap.add_argument("--results", type=Path, default=None)
    ap.add_argument("--scored", type=Path, default=None)
    ap.add_argument("--out_dir", type=Path, default=None)
    ap.add_argument("--no_plots", action="store_true")
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args(argv)
    if not (72 <= args.dpi <= 600):
        ap.error("--dpi must be between 72 and 600.")
    return args


def print_dry_run(pre, args):
    a122.sec("DRY RUN: nothing was drawn or written")
    print(f"  Event window {pre['ev_start']} to {pre['ev_end']}.")
    print(f"  Output folder (created only on a real run): {pre['out_dir']}")
    names = list(FILES.values())
    if args.window == "zoom":
        names = [n for n in names if not (n.startswith("fig_") and not n.endswith("_zoom.png"))]
    elif args.window == "full":
        names = [n for n in names if not n.endswith("_zoom.png")]
    print("  Files that would be written there: " + ", ".join(names))
    print("\nReady:  python -u plot_series122_hydrographs.py 2>&1 | tee series122_hydrographs.log")


def run(args, tee):
    paths = locate(args)
    pre = preflight(paths)
    d = pre["df"]
    picks = choose_runs(d)

    gA = a122.qout_to_grid(pre["ta_path"])
    qa_event = window_frame({"A": gA}, pre["ev_start"], pre["ev_end"])["A"]
    zoom, zinfo = None, None
    if args.window != "full":
        zoom, zinfo = resolve_zoom(args, qa_event)
    elif args.zoom_start or args.zoom_end:
        print("\n  NOTE  --zoom_start / --zoom_end are ignored with --window full.")
    print_zoom(zoom, zinfo, qa_event, args.window)
    if args.dry_run:
        print_dry_run(pre, args)
        return 0
    zoom_note = ("window set by you" if zinfo and zinfo["source"] != "the rule (Truth A's volume)"
                 else f"window from the rule: {ZOOM_PAD_MIN} min before {ZOOM_LO_FRAC * 100:g}% to {ZOOM_PAD_MIN} min after "
                      f"{ZOOM_HI_FRAC * 100:g}% of Truth A's volume")

    t0 = time.time()
    trio = picks["trio"]
    roles = {"best_B": d.loc[picks["best_B"]], "best_A": d.loc[picks["best_A"]]}
    trio_rows = []
    if trio["ok"]:
        trio_rows = [(n, d.loc[i]) for (n, _, _), i in zip(KBANDS, trio["rows"])]
    needed = {r["run_id"] for r in roles.values()} | {r["run_id"] for _, r in trio_rows}

    gB = a122.qout_to_grid(pre["tb_path"])
    runs = {rid: a122.qout_to_grid(pre["kept_paths"][rid]) for rid in sorted(needed)}
    self_check(pre, picks, gA, gB, runs)

    # ---- truth A against truth B ----
    fr_ab = window_frame({"A": gA, "B": gB}, pre["ev_start"], pre["ev_end"])
    num, ratio_tab, pts, ts_ab, lost = truth_ab_analysis(fr_ab, zoom)
    print_truth_ab(num, ratio_tab)

    # ---- best runs ----
    a122.sec("5. The best run against each truth  (figure 2)")
    frames, info, pick_rows = {}, {}, []
    role_print = []
    for key, suf, gT in (("best_B", "B", gB), ("best_A", "A", gA)):
        row = roles[key]
        frames[key] = window_frame({"truth": gT, "run": runs[row["run_id"]]}, pre["ev_start"], pre["ev_end"])
        info[key] = row
        c = compare_numbers(frames[key], "truth", "run", zoom)
        role_print.append((f"best against truth {suf}", row["run_id"], row["Kchan"], row[f"kge_2012_{suf}"],
                           (row[f"kge_r_{suf}"], row[f"kge_gamma_{suf}"], row[f"kge_beta_{suf}"]), c))
        pick_rows.append(dict(role=key, compared_with=f"truth {suf}", run_id=row["run_id"], draw=short_id(row["run_id"]),
                              Kchan=row["Kchan"], cv=row[a122.P_CV], flowexp=row[a122.P_FE], n=row[a122.P_N],
                              dial=row["dial"], **{f"{m}_{s}": row[f"{m}_{s}"] for m in SCORE_COLS for s in "AB"}, **c))
    print_compare("(compared with the truth each was picked for)", role_print)

    # ---- similar-routing trio ----
    a122.sec("6. Three runs with similar routing and different Kchan  (figure 3)")
    fr_trio = None
    if trio["ok"]:
        fr_trio = window_frame({"A": gA, "B": gB, **{n: runs[r["run_id"]] for n, r in trio_rows}},
                               pre["ev_start"], pre["ev_end"])
        tprint = []
        for n, r in trio_rows:
            c = compare_numbers(fr_trio, "B", n, zoom)
            tprint.append((f"{n} Kchan band", r["run_id"], r["Kchan"], r["kge_2012_B"],
                           (r["kge_r_B"], r["kge_gamma_B"], r["kge_beta_B"]), c))
            pick_rows.append(dict(role=f"trio_{n}", compared_with="truth B", run_id=r["run_id"], draw=short_id(r["run_id"]),
                                  Kchan=r["Kchan"], cv=r[a122.P_CV], flowexp=r[a122.P_FE], n=r[a122.P_N], dial=r["dial"],
                                  **{f"{m}_{s}": r[f"{m}_{s}"] for m in SCORE_COLS for s in "AB"}, **c))
        print("  All three are compared with truth B. Routing is matched on the Series 121 dial, not on cv, flowexp and n "
              "one by one, so those three can differ a little between the runs (printed in section 1). "
              "This is NOT Kchan alone.")
        print_compare("(compared with truth B)", tprint)
        print("  Reading it: the three share their routing as far as the dial can tell, so the differences between them are "
              "mostly the effect of Kchan, plus whatever small routing differences remain. Look at figure 3 for WHEN the gap opens. "
              "The clean test of Kchan alone is run_cc_clean_122_tribs6.py.")
    else:
        print(f"  Figure 3 is not drawn: {trio['reason']}.")

    a122.sec("7. CAUTIONS")
    for line in ("One synthetic truth pair, one storm (12 Aug 2014, 1.0x), Ks 7.0 and f 0.012 fixed, no noise. Model output, not observations.",
                 "A handful of hydrographs illustrates a pattern. It does not prove one.",
                 "The kept files are outlet-only: they show what leaves the basin, not where in the channel network the water is lost.",
                 "The zoom window follows a fixed rule on Truth A's volume. A pattern seen only inside it is not a claim about the other hours;"
                 " section 4 says how much of the loss falls outside it.",
                 "Why a pattern appears is a hypothesis to check in the tRIBS code, not something these curves show."):
        print("  * " + line)

    # ---- write ----
    a122.sec("Files written to " + str(pre["out_dir"]))
    out = pre["out_dir"]
    out.mkdir(parents=True, exist_ok=True)
    num_row = {k: (v.strftime("%Y-%m-%d %H:%M") if isinstance(v, pd.Timestamp) else v) for k, v in num.items()}
    pick_df = pd.DataFrame(pick_rows)
    for col in ("largest_gap_time",):
        if col in pick_df.columns:
            pick_df[col] = pick_df[col].map(lambda v: v.strftime("%Y-%m-%d %H:%M") if isinstance(v, pd.Timestamp) else v)
    ts_runs = pd.DataFrame({"time": fr_ab.index, "truth_A_m3s": fr_ab["A"].to_numpy(), "truth_B_m3s": fr_ab["B"].to_numpy()})
    for key in ("best_B", "best_A"):
        s = frames[key]["run"].reindex(fr_ab.index)
        ts_runs[f"{key}__{short_id(roles[key]['run_id']).replace(' ', '')}"] = s.to_numpy()
    if fr_trio is not None:
        for n, r in trio_rows:
            ts_runs[f"trio_{n}__{short_id(r['run_id']).replace(' ', '')}"] = fr_trio[n].reindex(fr_ab.index).to_numpy()
    for name, tab in (("truthab", pd.DataFrame([num_row])), ("ratio", ratio_tab), ("truthab_ts", ts_ab),
                      ("picks", pick_df), ("picks_ts", ts_runs)):
        tab.to_csv(out / FILES[name], index=False, float_format="%.8g")
        print(f"  saved {FILES[name]}  ({len(tab)} rows)")

    if not args.no_plots:
        print("\nDrawing the figures...")
        views = []
        if zoom is not None:
            views.append(("z", view_of(fr_ab, zoom)))
        if args.window in ("full", "both"):
            views.append(("", view_of(fr_ab, None)))
        figs = []
        for sfx, view in views:
            figs.append((f"fig1{sfx}", lambda v=view, s=sfx: figure_1(fr_ab, num, pts, lost, out / FILES[f"fig1{s}"], args.dpi, v, zoom_note)))
            figs.append((f"fig2{sfx}", lambda v=view, s=sfx: figure_2(frames, info, out / FILES[f"fig2{s}"], args.dpi, v, zoom_note)))
            if fr_trio is not None:
                figs.append((f"fig3{sfx}", lambda v=view, s=sfx: figure_3(fr_trio, trio_rows, num, trio["spread_dial"],
                                                                          out / FILES[f"fig3{s}"], args.dpi, v, zoom_note)))
        for name, fn in figs:
            try:
                fn()
            except Exception as exc:                                              # noqa: BLE001  a figure problem must not lose the tables
                print(f"  (figure {FILES[name]} skipped: {type(exc).__name__}: {exc})")
    try:
        (out / FILES["report"]).write_text(tee.buf.getvalue())
        print(f"\n  saved {FILES['report']}  (a copy of everything printed above)")
    except Exception as exc:                                                       # noqa: BLE001
        print(f"\n  (could not save {FILES['report']}: {exc})")
    print(f"\nDone in {time.time() - t0:.0f} s. This script only reads your results; it changed none of them.")
    return 0


def main(argv=None):
    args = parse_args(argv)
    tee = a122.Tee(sys.stdout)
    sys.stdout = tee
    print(f"Series 122 hydrographs (zoom version)   python {sys.version.split()[0]}, pandas {pd.__version__}, "
          f"numpy {np.__version__}   {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    code = 1
    try:
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        code = run(args, tee)
    except GateFailure as exc:
        print("\n" + "!" * 94)
        print(f"STOPPED (self-check failed): {exc}")
        print("!" * 94)
        code = 3
    except DataProblem as exc:
        print("\n" + "!" * 94)
        print(f"STOPPED: {exc}")
        print("!" * 94)
        print("Nothing was changed by this refusal. Fix the cause above, then run the command again.")
        code = 1
    except KeyboardInterrupt:
        print("\nStopped by you. Nothing in your results was changed.")
        code = 130
    finally:
        sys.stdout = tee.screen
    return code


if __name__ == "__main__":
    sys.exit(main())
