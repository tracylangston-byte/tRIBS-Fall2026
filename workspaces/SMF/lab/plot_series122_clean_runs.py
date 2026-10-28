"""
plot_series122_clean_runs.py
============================
Series 122 follow-up: READ and DRAW the twelve clean runs made by run_cc_clean_122_tribs6.py.

It runs NO tRIBS and changes none of your files. It only reads the results table and the kept outlet files that
run_cc_clean_122_tribs6.py wrote, and the two synthetic truths.

THE QUESTION
------------
The 400-run Series 122 sweep changed Kchan, cv, flowexp and n all at once, and those three routing parameters trade off
against each other (the "dial"), so the sweep could not show what Kchan alone does. The twelve clean runs change ONE thing
at a time, with Ks, f, cv, flowexp and n held at the truth's values:
    * 2 controls: ctlA (no channel loss, = truth A) and ctlB (loss ON, Kchan 201, = truth B);
    * 5 "rain error" runs: no loss, but the rain scaled down a little (x0.90 ... x0.99);
    * 5 "Kchan only" runs: loss ON, Kchan 10, 30, 100, 300, 1000 mm/hr.
So this script can ask four plain questions:
    1. DOSE-RESPONSE.  When only Kchan changes, how much volume and peak does the outlet lose?  (figure 1)
    2. LOSS OR RAIN ERROR?  Take the volume that truth B's channel loss removes (about 10%).  A rain error can remove the
       same volume.  Does it cut the PEAK by the same amount?  If a rain error cuts the peak much more than the loss does,
       the peak-and-volume pair tells the two apart.  (figure 2)
    3. TIMING AND SHAPE.  Draw the loss next to a rain error of the same volume, zoomed on the storm: WHEN is the water
       missing, and is it the same share of every flow (rising limb and recession kept apart)?  (figure 3)
    4. THE TAIL.  Cut truth A's storm into rise, recession and a long tail. How much of each run's missing water comes
       from each part, and does a run's outlet stop flowing for good after the flood?  (section 7, figure 4)

FOUR FIGURES (and the numbers behind them, printed and saved)
-------------------------------------------------------------
  1. fig_series122_c1_kchan_dose_response.png   volume and peak change against Kchan (log scale), and the KGE_2012 of each run
     against truth B and against truth A.  Kchan 201 is the ctlB run and is marked.
  2. fig_series122_c2_loss_vs_rain.png          percent of PEAK lost against percent of VOLUME lost, one line for the Kchan runs
     (loss) and one for the rain runs. A rain error that removes the same volume as truth B's loss is found by interpolating
     between the two rain runs on either side of it (a straight line between neighbours; the two runs used are printed).
  3. fig_series122_c3_same_volume_zoom.png      (and the full 20-hour version without _zoom) the loss (ctlB) and the rain run(s) with
     nearly the same volume, against truth A: hydrographs, what is missing through time, the running volume missing, and the
     ratio to truth A at each flow level (panel d: open dots = rising limb, filled dots = recession).
  4. fig_series122_c4_tail_log.png              the flow after the peak on a LOG scale: left, to the end of the scoring window;
     right, the first 3 hours from when truth A is down to 10% of its peak. A log scale cannot show zero, so a run's line stops
     where its flow is 0 and a thick bar along the bottom marks the time its outlet is dry.

HOW THE RAIN RUN(S) FOR FIGURE 3 ARE PICKED (a rule, fixed in the script before looking)
----------------------------------------------------------------------------------------
  Take the volume change of ctlB against truth A. The rain runs whose volume change lies just below and just above it are
  drawn (the nearer one solid, the other lighter). If ctlB's volume change is outside the range the rain runs cover, only the
  nearest rain run is drawn and the matched-volume interpolation is NOT done (the script says so; it never extrapolates).

THE ZOOM WINDOW
---------------
  Same rule as plot_series122_hydrographs.py (this script imports it): from 30 minutes before 1% of truth A's volume has passed to
  30 minutes after 97.5% has passed. Zooming changes only what is drawn; every number is computed over the whole 20-hour
  scoring window. --zoom_start / --zoom_end set it by hand. --window full|zoom|both chooses which version of figure 3 is drawn.

WHAT IT READS  (only reads)
---------------------------
    calibration_work/03_comparisons/summary_tables/series122_clean_runs.csv          (written by run_cc_clean_122_tribs6.py)
    calibration_work/03_comparisons/summary_tables/series122_clean_runs_design.json  (optional; used to list what was planned)
    calibration_work/kept_qout_122_clean/<run_id>_Outlet.qout                        the kept outlet files
    calibration_work/synth_truth/*.qout                                              truth A
    calibration_work/synth_truth_location/SMF_20140812_110_cc201p0_loc_K7p0_f0p012_Outlet.qout   truth B

WHAT IT WRITES  (a NEW folder; it refuses to write into the folder that holds the results table)
------------------------------------------------------------------------------------------------
  calibration_work/03_comparisons/sensitivity_plots/series122_clean_runs/
    series122_c_report.txt           everything printed on screen
    series122_c_table.csv            one row per run: volume, peak, timing, KGE_2012 and parts against both truths
    series122_c_matched_volume.csv   the matched-volume comparison of the loss and the rain error
    series122_c_ratio_bands.csv      run / truth A at each flow level and limb, for ctlB and the rain run(s) of figure 3
    series122_c_phases.csv           every run: its missing water by phase (rise / recession / long tail), and when its outlet goes dry
    series122_c_timeseries.csv       ctlA, ctlB and the rain run(s) of figure 3 through the 20-hour window
    fig_series122_c1_kchan_dose_response.png   fig_series122_c2_loss_vs_rain.png
    fig_series122_c3_same_volume_zoom.png      fig_series122_c3_same_volume.png      fig_series122_c4_tail_log.png

HOW TO RUN  (from the lab/ folder, next to analyze_series122_kchan.py and plot_series122_hydrographs.py, which this script uses)
----------------------------------------------------------------------------------------------------------------------------
    cd /workspaces/tRIBS-Fall2026/workspaces/SMF/lab
    python plot_series122_clean_runs.py --dry_run                          <- checks inputs and shows the zoom window; writes nothing
    python -u plot_series122_clean_runs.py 2>&1 | tee series122_clean_compare.log

It takes a few seconds. Run it AFTER run_cc_clean_122_tribs6.py has finished (python run_cc_clean_122_tribs6.py --audit says PASS).

BEFORE IT DRAWS ANYTHING IT CHECKS ITSELF, and stops if a check fails:
  * truth A and truth B are the exact files (checksums) Series 122 used; every kept outlet file matches the checksum recorded for it;
  * both controls passed in the run script (control_check PASS), and ctlB reproduces truth B and truth A's score of it matches Series 110;
  * for EVERY run, the scores recomputed here from the kept outlet file (against truth A and against truth B: KGE_2012, r, gamma,
    beta, PBIAS, volume ratio, peak ratio, peak-time shift) agree (to 1e-6) with what the tRIBS scorer and the run script
    recorded. Three separate pieces of code have to agree before a picture is drawn.

OPTIONS
-------
    --dry_run       check the inputs (including the re-scoring) and show the zoom window, then stop
    --window W      zoom, full or both (default both): which versions of figure 3 to draw
    --zoom_start T  --zoom_end T    set the zoom window yourself, e.g. "2014-08-12 17:00"
    --results PATH  the clean-runs results table      --out_dir PATH  where to write
    --no_plots      numbers and tables only           --dpi N   (default 150)

STATUS  (what has and has not been tested)
------------------------------------------
  VERSION 2.  Version 1 (without section 7 and figure 4) ran on your real files on 2026-10-27 and passed its own checks: every run
  re-scored to the recorded values and ctlB matched Series 110. This version 2 has NOT yet been run on your real files.
  What was checked for version 2:
    * the new phase and "outlet dry" calculations were run on your real timeseries file (ctlA, ctlB, rain x0.95, rain x0.97) and
      agree with separate plain-loop calculations: ctlB dry from 20:05 for 16.0 h, truth A carrying 2,799 m3 in that time;
    * on a stand-in lab (toy hydrographs) every number in the table, matched-volume, ratio-band (now by limb), phase and time-series
      files was recomputed separately from the raw kept files (plain loops, np.interp, the scorer's own reading and scoring code)
      and agreed to 1e-6, for the full 12-run set and for four odd sets (loss bigger than any rain run covers, 3 runs per curve, a
      failed run missing, 2 rain and 2 Kchan runs); edge cases of the phase and dry functions (no long tail, a run dry straight after
      the peak, a late flow ending dryness, never dry, no flow in truth A) were tested one by one;
    * it still stops (exit 1 or 3, nothing written, no input file changed) on the same 18 damaged inputs, and refuses bad options;
    * the figures were drawn from your real table and timeseries files and looked at, including the label placement in figure 2.
  NOT tested: the full version 2 on your real kept files. If it stops, it says why and changes nothing.
  Things the numbers can not tell you: why the outlet runs dry (a hypothesis to check in the tRIBS code), and whether a real gauge
  could see the difference (flows in the long tail are tiny; the gauge's resolution has not been checked).

READ THIS BEFORE QUOTING ANYTHING FROM THE FIGURES
--------------------------------------------------
  * One synthetic truth pair, one storm (12 Aug 2014), Ks 7.0 and f 0.012 fixed, no noise. Model outputs, not observations.
  * In these runs the ROUTING (cv, flowexp, n) is held at the truth. In a real calibration the routing parameters are free, and
    Series 121 and 122 found they can trade peak against shape (the "dial"). So a difference in peak-and-volume pairing found
    here may be weaker, or absorbed, once routing is free. These runs show what the loss and a rain error do to the hydrograph;
    they do not show that a calibration could tell them apart.
  * Rain is scaled by one factor over the whole storm and both gauges. A real rain error varies in time and space.
  * Interpolating between two rain runs assumes the relation is close to a straight line between them; the two runs are 2 to 5
    percent of rain apart, so the interpolated peak is an estimate, not a model run.
  * The kept files are OUTLET only: they show what leaves the basin, not WHERE in the channel network the water is lost.
  * Why a pattern appears is a hypothesis to check in the tRIBS code, not something these curves show.
"""

import argparse
import json
import sys
import time
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import analyze_series122_kchan as a122
    import plot_series122_hydrographs as h122
except ImportError as _exc:                                                          # noqa: BLE001
    print("STOPPED: this script reuses code from analyze_series122_kchan.py and plot_series122_hydrographs.py and could not import "
          f"one of them ({_exc}).\nPut all three scripts in the same lab/ folder and run from there:\n"
          "    cd /workspaces/tRIBS-Fall2026/workspaces/SMF/lab")
    sys.exit(1)

DataProblem, GateFailure = a122.DataProblem, a122.GateFailure
fmt_t, minutes_between, first_cross, window_frame = h122.fmt_t, h122.minutes_between, h122.first_cross, h122.window_frame

# ======================================================================
# FILES
# ======================================================================
RESULTS_NAME = "series122_clean_runs.csv"
DESIGN_NAME = "series122_clean_runs_design.json"
FAILED_NAME = "series122_clean_runs_FAILED.csv"
KEPT_DIR_NAME = "kept_qout_122_clean"
OUT_SUBFOLDER = "series122_clean_runs"
FILES = dict(
    report="series122_c_report.txt", table="series122_c_table.csv", matched="series122_c_matched_volume.csv",
    bands="series122_c_ratio_bands.csv", ts="series122_c_timeseries.csv", phases="series122_c_phases.csv",
    fig1="fig_series122_c1_kchan_dose_response.png", fig2="fig_series122_c2_loss_vs_rain.png",
    fig3="fig_series122_c3_same_volume.png", fig3z="fig_series122_c3_same_volume_zoom.png",
    fig4="fig_series122_c4_tail_log.png",
)

# ======================================================================
# RULES AND TOLERANCES: fixed before looking at any result. Do not change them to get a different picture.
# ======================================================================
CHECK_TOL = 1e-6                 # re-scoring here must agree with the recorded scores to this (relative to 1 + |value|)
CTL_PBIAS_TOL = 0.01             # a control must reproduce its truth: |PBIAS| below this (percent) ...
CTL_KGE_MIN = 0.9999             # ... and KGE_2012 above this
REC_PBIAS_TOL = 1e-3             # ctlB against truth A must match the Series 110 record to this (percent) ...
REC_KGE_TOL = 1e-5               # ... and in KGE_2012 to this
MIN_RUNS_PER_CURVE = 2           # fewer rain runs or Kchan runs than this and no curve can be drawn
PEAK_DISTINCT_PP = 2.0           # "the pairing differs": at matched volume the peak change differs by at least this many percentage points
MIN_VOL_FOR_RATIO_PCT = 0.5      # peak change / volume change is only quoted when the volume changed by at least this much (percent)
DRY_FLOW = 1e-6                  # m3/s: at or below this the outlet counts as dry (tRIBS writes exact zeros there)
TAIL_FRAC = 0.01                 # the "long tail" of truth A's storm begins once A has fallen below this fraction of its peak
TAIL_ZOOM_FROM = 0.10            # figure 4, right panel: starts when truth A has fallen to this fraction of its peak ...
TAIL_ZOOM_HOURS = 3.0            # ... and shows this many hours
SCORE_COLS = ("kge_2012", "kge_r", "kge_gamma", "kge_beta", "pbias_pct")
KIND_LABEL = {"control_noloss": "control A", "control_loss": "control B", "rain": "rain error", "kchan": "Kchan only"}
BLUE, GOLD, ORANGE, INK, INK_SOFT, TEAL = a122.BLUE, a122.GOLD, a122.ORANGE, a122.INK, a122.INK_SOFT, h122.TEAL
LOSS_COL, RAIN_COL = BLUE, GOLD          # figure 2 and 3: blue = channel loss, gold = rain error

# columns the run script's results table must hold
REQUIRED = (["run_id", "run_name", "run_kind", "kchan_asked", "rain_scale", "control_check", "truth_A_md5", "truth_B_md5",
             "qout_md5", "event_start", "event_end", "kge_2012", "pbias_pct", "kge_r", "kge_beta", "kge_gamma"]
            + [f"vs_{s}_{m}" for s in "AB" for m in ("pbias_pct", "kge_2012", "vol_ratio", "peak_ratio", "peak_shift_min")])


# ======================================================================
# Small helpers
# ======================================================================
def pct(x):
    return float("nan") if x is None or not np.isfinite(x) else 100.0 * (x - 1.0)


def label_of(row):
    k = row["run_kind"]
    if k == "kchan":
        return f"Kchan {float(row['kchan_asked']):g}"
    if k == "rain":
        return f"rain x{float(row['rain_scale']):g}"
    return "ctlA (no loss)" if k == "control_noloss" else f"ctlB (Kchan {float(row['kchan_asked']):g})"


def close(a, b):
    if (a is None or not np.isfinite(a)) and (b is None or not np.isfinite(b)):
        return True
    if a is None or b is None or not np.isfinite(a) or not np.isfinite(b):
        return False
    return abs(a - b) <= CHECK_TOL * (1.0 + abs(b))


# ======================================================================
# Locating and checking the inputs
# ======================================================================
def locate(args):
    root = Path.cwd().parent
    calib = root / "calibration_work"
    summ = calib / "03_comparisons" / "summary_tables"
    plots = calib / "03_comparisons" / "sensitivity_plots"
    return dict(
        calib=calib, summ=summ,
        results=args.results or (summ / RESULTS_NAME),
        design=summ / DESIGN_NAME, failed=summ / FAILED_NAME,
        kept=calib / KEPT_DIR_NAME,
        truth_a_dir=calib / "synth_truth",
        truth_b=calib / "synth_truth_location" / a122.TRUTH_B_NAME,
        out=args.out_dir or (plots / OUT_SUBFOLDER),
    )


def preflight(paths):
    """Read-only checks. Returns what the later steps need. Raises DataProblem on a real problem."""
    a122.sec("0. Inputs and checks (nothing is written yet)")
    rp = Path(paths["results"])
    if not rp.exists():
        raise DataProblem(f"Could not find the results table {rp}\nRun run_cc_clean_122_tribs6.py first (it writes that file), "
                          "or point to it with --results.")
    try:
        # checksums are read as text: a 32-character checksum made only of digits and one 'e' would otherwise be turned into a number
        df = pd.read_csv(rp, dtype={"truth_A_md5": str, "truth_B_md5": str, "qout_md5": str})
    except Exception as exc:                                                         # noqa: BLE001
        raise DataProblem(f"Could not read {rp} as a CSV: {exc}")
    miss = [c for c in REQUIRED if c not in df.columns]
    if miss:
        raise DataProblem("The results table is missing these column(s): " + ", ".join(miss)
                          + ". Is this series122_clean_runs.csv from run_cc_clean_122_tribs6.py?")
    print(f"  results table   {rp.name}: {len(df)} rows")
    if df["run_id"].duplicated().any():
        raise DataProblem("A run_id appears more than once in the results table.")
    print("  ok    every run_id is unique")
    bad_kind = sorted(set(df["run_kind"]) - set(KIND_LABEL))
    if bad_kind:
        raise DataProblem(f"Unknown run_kind value(s) {bad_kind}; expected {sorted(KIND_LABEL)}.")
    kinds = df["run_kind"].value_counts().to_dict()
    for k in ("control_noloss", "control_loss"):
        if kinds.get(k, 0) != 1:
            raise DataProblem(f"Expected exactly one {KIND_LABEL[k]} run ({k}) in the results table, found {kinds.get(k, 0)}.")
    for k, need in (("rain", "rain-error"), ("kchan", "Kchan-only")):
        if kinds.get(k, 0) < MIN_RUNS_PER_CURVE:
            raise DataProblem(f"Only {kinds.get(k, 0)} {need} run(s) in the results table; need at least {MIN_RUNS_PER_CURVE} "
                              "to draw a curve. Finish run_cc_clean_122_tribs6.py (python run_cc_clean_122_tribs6.py "
                              "--skip_existing) and run this again.")
        if kinds.get(k, 0) < 3:
            print(f"  NOTE  only {kinds[k]} {need} runs: the curve is a single line segment.")
    print(f"  ok    runs: 1 control A, 1 control B, {kinds['rain']} rain-error, {kinds['kchan']} Kchan-only")

    for c in ("kchan_asked", "rain_scale"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    if df.loc[df["run_kind"] == "kchan", "kchan_asked"].isna().any() or df.loc[df["run_kind"] == "rain", "rain_scale"].isna().any():
        raise DataProblem("A Kchan-only run has no Kchan value, or a rain run has no rain scale, in the results table.")
    loss_k = df.loc[df["run_kind"].isin(["kchan", "control_loss"]), "kchan_asked"]
    if loss_k.isna().any() or loss_k.duplicated().any() or df.loc[df["run_kind"] == "rain", "rain_scale"].duplicated().any():
        raise DataProblem("Two runs with the loss on share a Kchan (ctlB counts, at 201), or two rain runs share a rain scale, "
                          "or ctlB has no Kchan. Each point on a curve must be a different run.")

    starts, ends = df["event_start"].astype(str).unique(), df["event_end"].astype(str).unique()
    if len(starts) != 1 or len(ends) != 1:
        raise DataProblem(f"The results table has more than one event window ({list(starts)} to {list(ends)}).")
    ev_start, ev_end = starts[0], ends[0]
    print(f"  ok    one event window: {ev_start} to {ev_end}")

    ta_dir = paths["truth_a_dir"]
    ta = sorted(ta_dir.glob("*.qout")) if ta_dir.exists() else []
    if len(ta) != 1:
        raise DataProblem(f"Expected exactly one *.qout file in {ta_dir} (truth A), found {len(ta)}.")
    ta_md5 = a122.md5_of(ta[0])
    if ta_md5 != a122.TRUTH_A_MD5:
        raise DataProblem(f"Truth A ({ta[0].name}) has checksum {ta_md5}, not the {a122.TRUTH_A_MD5} that Series 122 was scored "
                          "against. Do not continue with a different file.")
    print(f"  ok    truth A  {ta[0].name}  [md5 {ta_md5[:12]}]  is the file Series 122 used")
    tb = Path(paths["truth_b"])
    if not tb.exists():
        raise DataProblem(f"Truth B not found: {tb}")
    tb_md5 = a122.md5_of(tb)
    if tb_md5 != a122.TRUTH_B_MD5:
        raise DataProblem(f"Truth B has checksum {tb_md5}, not the verified {a122.TRUTH_B_MD5}. Do not continue.")
    print(f"  ok    truth B  {tb.name}  [md5 {tb_md5[:12]}]  is the file verify_truth_B_122_tribs6.py verified")
    for col, want, nm in (("truth_A_md5", ta_md5, "A"), ("truth_B_md5", tb_md5, "B")):
        got = sorted(set(df[col].astype(str)))
        if got != [want]:
            raise DataProblem(f"The runs were scored against a different truth {nm} file than the one in the folder now "
                              f"(table says {got}, folder has {want}). Do not mix them.")
    print("  ok    every run in the table was scored against these two truth files")

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
        if a122.md5_of(f) != str(r.qout_md5):
            md5_bad.append(r.run_id)
    if missing:
        raise DataProblem(f"{len(missing)} run(s) have no kept outlet file (for example {missing[0]}). Run "
                          "python run_cc_clean_122_tribs6.py --audit to see which.")
    if md5_bad:
        raise DataProblem(f"{len(md5_bad)} kept outlet file(s) do not match the checksum the run script recorded (for example "
                          f"{md5_bad[0]}). The kept files must be the ones that were scored.")
    print(f"  ok    {len(kept_paths)} kept outlet files found; each checksum equals the one recorded")

    for k in ("control_noloss", "control_loss"):
        v = str(df.loc[df["run_kind"] == k, "control_check"].iloc[0])
        if v != "PASS":
            raise DataProblem(f"The {KIND_LABEL[k]} run's control_check is {v!r}, not PASS. The run script's gate did not pass, so "
                              "nothing here can be trusted. Run python run_cc_clean_122_tribs6.py --audit.")
    print("  ok    both controls have control_check PASS in the run script")

    notes = []
    dp = Path(paths["design"])
    if dp.exists():
        try:
            planned = json.loads(dp.read_text()).get("runs", [])
            have = set(df["run_id"])
            absent = [p.get("name", p.get("run_id")) for p in planned if p.get("run_id") not in have]
            if absent:
                notes.append(f"{len(absent)} planned run(s) have no result row yet: {', '.join(map(str, absent))}")
            print(f"  ok    design file lists {len(planned)} planned runs; {len(planned) - len(absent)} have results")
        except Exception as exc:                                                     # noqa: BLE001
            notes.append(f"could not read {dp.name}: {exc}")
    else:
        notes.append(f"{dp.name} not found (only used to list what was planned)")
    fp = Path(paths["failed"])
    if fp.exists():
        try:
            nf = len(pd.read_csv(fp))
            if nf:
                notes.append(f"{fp.name} lists {nf} failed or hung run(s); the curves below are missing them")
        except Exception:                                                            # noqa: BLE001
            pass
    for n in notes:
        print(f"  NOTE  {n}")

    out_dir = Path(paths["out"])
    if out_dir.resolve() == rp.resolve().parent:
        raise DataProblem(f"The output folder {out_dir} is the folder that holds the results table; refusing to write there. "
                          "Use the default or another --out_dir.")
    return dict(df=df.reset_index(drop=True), ev_start=ev_start, ev_end=ev_end, ta_path=ta[0], tb_path=tb,
                kept_paths=kept_paths, out_dir=out_dir)


# ======================================================================
# Scoring (recomputed here from the kept files)
# ======================================================================
def pair_metrics(gT, gR, ev_start, ev_end):
    """Score run gR against truth gT over the scoring window: the scorer's formulas (a122.score_pair) plus volume ratio,
    peak ratio and peak-time shift (run over truth) and the first-flow shift."""
    s = a122.score_pair(gT, gR, ev_start, ev_end)
    if s is None:
        raise DataProblem("A run does not overlap the scoring window; cannot score it.")
    fr = pd.DataFrame({"t": gT, "r": gR}).dropna().loc[ev_start:ev_end]
    t, r = fr["t"], fr["r"]
    thr = h122.FLOOR_FRAC * float(t.max())
    out = {k: s[k] for k in SCORE_COLS}
    out.update(vol_ratio=float(r.sum() / t.sum()), peak_ratio=float(r.max() / t.max()),
               peak_shift_min=minutes_between(t.idxmax(), r.idxmax()),
               first_flow_shift_min=minutes_between(first_cross(t, thr), first_cross(r, thr)))
    return out


def rescore_and_check(pre, gA, gB, grids):
    """Section 1. Recompute every run's scores from its kept file and compare with what the tRIBS scorer and the run script
    recorded. Returns {run_id: {"A": metrics, "B": metrics}}. Raises GateFailure on any disagreement."""
    a122.sec("1. Self-check: do the kept outlet files give back the scores on record?")
    df = pre["df"]
    scored, fails = {}, []
    print(f"  {'run':<22}{'largest difference':>20}   result   (scorer columns and run-script columns, against truth A and truth B)")
    for row in df.itertuples():
        rid = row.run_id
        scored[rid] = {s: pair_metrics(g, grids[rid], pre["ev_start"], pre["ev_end"]) for s, g in (("A", gA), ("B", gB))}
        worst, bad = 0.0, False
        pairs = [(f"scorer {m}", getattr(row, m), scored[rid]["A"][m]) for m in SCORE_COLS]      # the tRIBS scorer scores against truth A
        for s in "AB":
            for m in ("pbias_pct", "kge_2012", "vol_ratio", "peak_ratio", "peak_shift_min"):
                pairs.append((f"run script vs_{s}_{m}", getattr(row, f"vs_{s}_{m}"), scored[rid][s][m]))
        for name, theirs, mine in pairs:
            theirs = float(theirs)
            if (not np.isfinite(theirs)) and (not np.isfinite(mine)):
                continue
            diff = abs(mine - theirs) if np.isfinite(mine) and np.isfinite(theirs) else float("inf")
            worst = max(worst, diff)
            if not close(mine, theirs):
                bad = True
                fails.append(f"{row.run_name}: {name}: here {mine:.8g}, recorded {theirs:.8g}")
        print(f"  {row.run_name:<22}{worst:>20.3g}   {'FAIL' if bad else 'ok'}")
    # controls reproduce their truths; ctlB against truth A is Series 110's record
    for kind, suf, nm in (("control_noloss", "A", "ctlA against truth A"), ("control_loss", "B", "ctlB against truth B")):
        rid = df.loc[df["run_kind"] == kind, "run_id"].iloc[0]
        m = scored[rid][suf]
        ok = abs(m["pbias_pct"]) <= CTL_PBIAS_TOL and m["kge_2012"] >= CTL_KGE_MIN
        print(f"  {nm:<28} PBIAS {m['pbias_pct']:+.6f}%   KGE_2012 {m['kge_2012']:.8f}   {'ok' if ok else 'FAIL'}")
        if not ok:
            fails.append(f"{nm}: PBIAS {m['pbias_pct']:+.6f}% and KGE_2012 {m['kge_2012']:.8f} do not reproduce the truth")
    rid_b = df.loc[df["run_kind"] == "control_loss", "run_id"].iloc[0]
    mba = scored[rid_b]["A"]
    for m, tol in (("pbias_pct", REC_PBIAS_TOL), ("kge_2012", REC_KGE_TOL)):
        rec = a122.REC_B_VS_A[m]
        ok = abs(mba[m] - rec) <= tol
        print(f"  ctlB against truth A  {m:<10} here {mba[m]:>14.8f}   Series 110 recorded {rec:>14.8f}   {'ok' if ok else 'FAIL'}")
        if not ok:
            fails.append(f"ctlB against truth A {m}: here {mba[m]:.8f}, Series 110 recorded {rec:.8f}")
    if fails:
        raise GateFailure("the kept outlet files do not reproduce the scores on record, so nothing drawn from them can be trusted "
                          "yet:\n    - " + "\n    - ".join(fails) + "\n  Nothing was drawn or written. Paste this section to Claude.")
    print("  PASS  every run re-scores to the recorded values against both truths; both controls reproduce their truth; "
          "ctlB against truth A matches Series 110.")
    return scored


# ======================================================================
# The run table and the analysis numbers
# ======================================================================
def timing_metrics(gA, gR, ev_start, ev_end, zoom):
    """Where in time the run differs from truth A (whole scoring window)."""
    fr = window_frame({"A": gA, "R": gR}, ev_start, ev_end)
    dt = (fr.index[1] - fr.index[0]).total_seconds()
    gap = (fr["A"].cumsum() - fr["R"].cumsum()) * dt          # running volume missing compared with truth A (m3)
    tot = float(gap.iloc[-1])
    tA = fr["A"].idxmax()
    res = fr["R"] - fr["A"]
    iw = res.abs().idxmax()
    out = dict(volume_missing_m3=tot,
               gap_share_by_A_peak_pct=(100 * float(gap.loc[tA]) / tot) if tot > 0 else float("nan"),
               gap_share_after_zoom_pct=float("nan"),
               largest_gap_m3s=float(res.loc[iw]), largest_gap_time=iw)
    if zoom is not None and tot > 0:
        out["gap_share_after_zoom_pct"] = 100 * (tot - float(gap.loc[:zoom[1]].iloc[-1])) / tot
    return out


def build_table(pre, gA, grids, scored, zoom):
    rows = []
    for r in pre["df"].itertuples():
        a, b = scored[r.run_id]["A"], scored[r.run_id]["B"]
        tm = timing_metrics(gA, grids[r.run_id], pre["ev_start"], pre["ev_end"], zoom)
        curve = {"control_noloss": "anchor", "control_loss": "loss", "kchan": "loss", "rain": "rain"}[r.run_kind]
        row = dict(run_name=r.run_name, run_kind=r.run_kind, curve=curve, label=label_of(r._asdict()),
                   Kchan=r.kchan_asked if r.run_kind in ("kchan", "control_loss") else float("nan"),
                   rain_scale=r.rain_scale, run_id=r.run_id,
                   vol_change_vs_A_pct=pct(a["vol_ratio"]), peak_change_vs_A_pct=pct(a["peak_ratio"]),
                   peak_shift_vs_A_min=a["peak_shift_min"], first_flow_shift_vs_A_min=a["first_flow_shift_min"],
                   kge_2012_vs_A=a["kge_2012"], kge_r_vs_A=a["kge_r"], kge_gamma_vs_A=a["kge_gamma"], kge_beta_vs_A=a["kge_beta"],
                   pbias_vs_A_pct=a["pbias_pct"],
                   vol_change_vs_B_pct=pct(b["vol_ratio"]), peak_change_vs_B_pct=pct(b["peak_ratio"]),
                   peak_shift_vs_B_min=b["peak_shift_min"],
                   kge_2012_vs_B=b["kge_2012"], kge_r_vs_B=b["kge_r"], kge_gamma_vs_B=b["kge_gamma"], kge_beta_vs_B=b["kge_beta"],
                   pbias_vs_B_pct=b["pbias_pct"], **tm)
        v = row["vol_change_vs_A_pct"]
        row["peak_per_volume"] = (row["peak_change_vs_A_pct"] / v) if (np.isfinite(v) and abs(v) >= MIN_VOL_FOR_RATIO_PCT) else float("nan")
        rows.append(row)
    return pd.DataFrame(rows)


def match_volume(target, rain_tab):
    """The rain error that removes the same volume as the loss: linear interpolation of peak change (and rain scale) between
    the two neighbouring rain points, ctlA (no change) included as the end of the rain line. Never extrapolates."""
    pts = rain_tab.sort_values("vol_change_vs_A_pct").reset_index(drop=True)
    v = pts["vol_change_vs_A_pct"].to_numpy(float)
    out = dict(ok=False, target_vol_pct=float(target), reason="")
    if not (v.min() <= target <= v.max()):
        out["reason"] = (f"the loss changes the volume by {target:+.2f}%, outside the range the rain runs cover "
                         f"({v.min():+.2f}% to {v.max():+.2f}%), so nothing is interpolated")
        return out
    j = int(np.searchsorted(v, target, side="right") - 1)
    j = min(max(j, 0), len(v) - 2)
    lo, hi = pts.loc[j], pts.loc[j + 1]
    w = (target - lo["vol_change_vs_A_pct"]) / (hi["vol_change_vs_A_pct"] - lo["vol_change_vs_A_pct"])
    out.update(ok=True, w=float(w), lower_run=lo["label"], upper_run=hi["label"],
               lower_vol_pct=float(lo["vol_change_vs_A_pct"]), upper_vol_pct=float(hi["vol_change_vs_A_pct"]),
               peak_pct=float(lo["peak_change_vs_A_pct"] + w * (hi["peak_change_vs_A_pct"] - lo["peak_change_vs_A_pct"])),
               rain_scale=float(lo["rain_scale"] + w * (hi["rain_scale"] - lo["rain_scale"])))
    return out


def pick_rain_runs(target, rain_only):
    """The rain run(s) drawn in figure 3: the ones whose volume change lies just below and just above the loss's, the nearer first."""
    v = rain_only["vol_change_vs_A_pct"].to_numpy(float)
    below = np.where(v <= target)[0]
    above = np.where(v >= target)[0]
    chosen = []
    if len(below):
        chosen.append(int(below[np.argmax(v[below])]))
    if len(above):
        chosen.append(int(above[np.argmin(v[above])]))
    chosen = list(dict.fromkeys(chosen))
    chosen.sort(key=lambda i: abs(v[i] - target))
    return [rain_only.iloc[i] for i in chosen]


LIMBS = (("rising", "rising limb (up to A's peak)"), ("recession", "recession (after A's peak)"))


def ratio_bands(frs):
    """Run / truth A at each flow level (A's flow above 5% of its peak), whole window, rising limb and recession kept apart
    (the same flow level occurs twice in a storm, once going up and once coming down). frs = {name: (fr with 'A' and the run)}."""
    rows = []
    for name, (fr, col) in frs.items():
        peak = float(fr["A"].max())
        t_peak = fr["A"].idxmax()
        sub = fr[fr["A"] >= h122.FLOOR_FRAC * peak]
        ratio = sub[col] / sub["A"]
        frac = sub["A"] / peak
        for limb, _ in LIMBS:
            on_limb = (sub.index <= t_peak) if limb == "rising" else (sub.index > t_peak)
            for lo, hi in h122.FLOW_BANDS:
                sel = on_limb & ((frac >= lo) & ((frac < hi) | (hi >= 1.0))).to_numpy()
                r = ratio[sel]
                rows.append(dict(run=name, limb=limb, band=f"{int(lo * 100)} to {int(hi * 100)}% of A's peak", n_points=int(len(r)),
                                 median_ratio=float(r.median()) if len(r) else float("nan"),
                                 min_ratio=float(r.min()) if len(r) else float("nan"),
                                 max_ratio=float(r.max()) if len(r) else float("nan")))
    return pd.DataFrame(rows)


# ======================================================================
# The tail: which part of the storm the missing water comes from, and does the outlet stop flowing?
# ======================================================================
def phase_bounds(qa):
    """Cut truth A's storm into three phases, by a rule on truth A only, so every run is cut at the same times:
    rise = the start of the window up to A's peak (the peak included); recession = after the peak, while A is at or above TAIL_FRAC
    of its peak; long tail = the rest. Returns (t_peak, t_tail); t_tail is None if A never falls that low inside the window."""
    peak = float(qa.max())
    if not np.isfinite(peak) or peak <= 0:
        raise DataProblem("Truth A has no flow inside the scoring window, so its storm cannot be cut into phases.")
    t_peak = qa.idxmax()
    after = qa[qa.index > t_peak]
    low = after[after < TAIL_FRAC * peak]
    return t_peak, (low.index[0] if len(low) else None)


def run_phases(qa, qr, t_peak, t_tail):
    """Volume (m3) that run qr is missing compared with truth A in each phase (positive = the run has less water)."""
    dt = (qa.index[1] - qa.index[0]).total_seconds()
    gap = (qa - qr) * dt
    later = gap[gap.index > t_peak]
    if t_tail is None:
        return dict(rise_m3=float(gap[gap.index <= t_peak].sum()), recession_m3=float(later.sum()), tail_m3=0.0)
    return dict(rise_m3=float(gap[gap.index <= t_peak].sum()),
                recession_m3=float(later[later.index < t_tail].sum()),
                tail_m3=float(gap[gap.index >= t_tail].sum()))


def dry_info(qa, qr, t_peak):
    """Does the run's outlet stop flowing for good after truth A's peak?  dry_from = the first time from which the run's flow stays
    at or below DRY_FLOW to the end of the window (None if it never settles there). Also how long that lasts, and how much water
    truth A carries during it."""
    dt = (qa.index[1] - qa.index[0]).total_seconds()
    out = dict(dry_from=None, dry_hours=0.0, a_volume_dry_m3=0.0)
    after = qr[qr.index > t_peak]
    if len(after) == 0 or after.iloc[-1] > DRY_FLOW:
        return out
    wet = after[after > DRY_FLOW]
    start = after.index[0] if len(wet) == 0 else after.index[after.index > wet.index[-1]][0]
    out.update(dry_from=start, dry_hours=((qa.index[-1] - start).total_seconds() + dt) / 3600.0,
               a_volume_dry_m3=float(qa.loc[start:].sum() * dt))
    return out


def build_phases(pre, gA, grids, tab):
    """One row per run: where in truth A's storm the run's missing water is, and whether its outlet runs dry."""
    qa_all = window_frame({"A": gA}, pre["ev_start"], pre["ev_end"])["A"]
    t_peak, t_tail = phase_bounds(qa_all)
    dtA = (qa_all.index[1] - qa_all.index[0]).total_seconds()
    volA = float(qa_all.sum() * dtA)
    a_ph = run_phases(qa_all, qa_all * 0.0, t_peak, t_tail)
    rows = []
    for r in pre["df"].itertuples():
        fr = window_frame({"A": gA, "R": grids[r.run_id]}, pre["ev_start"], pre["ev_end"])
        ph = run_phases(fr["A"], fr["R"], t_peak, t_tail)
        dr = dry_info(fr["A"], fr["R"], t_peak)
        total = ph["rise_m3"] + ph["recession_m3"] + ph["tail_m3"]
        want = float(tab.loc[tab["run_id"] == r.run_id, "volume_missing_m3"].iloc[0])
        if abs(total - want) > 1e-6 * (1.0 + abs(want)):
            raise GateFailure(f"the three phases of {r.run_name} add up to {total:.4f} m3 but section 3 says {want:.4f} m3 missing; "
                              "the two calculations should agree.")
        share = (lambda v: 100.0 * v / total) if total > 0 else (lambda v: float("nan"))
        rows.append(dict(
            run_name=r.run_name, run_kind=r.run_kind, label=label_of(r._asdict()),
            rise_missing_m3=ph["rise_m3"], recession_missing_m3=ph["recession_m3"], tail_missing_m3=ph["tail_m3"],
            total_missing_m3=total, total_missing_pct_of_A=100.0 * total / volA,
            rise_share_of_gap_pct=share(ph["rise_m3"]), recession_share_of_gap_pct=share(ph["recession_m3"]),
            tail_share_of_gap_pct=share(ph["tail_m3"]),
            dry_from=dr["dry_from"].strftime("%Y-%m-%d %H:%M") if dr["dry_from"] is not None else "",
            dry_hours=dr["dry_hours"], A_volume_in_dry_period_pct=100.0 * dr["a_volume_dry_m3"] / volA,
            A_volume_in_dry_period_m3=dr["a_volume_dry_m3"],
            rise_end=t_peak.strftime("%Y-%m-%d %H:%M"),
            recession_end=(t_tail - pd.Timedelta(seconds=dtA)).strftime("%Y-%m-%d %H:%M") if t_tail is not None else "",
            A_rise_pct_of_A=100.0 * a_ph["rise_m3"] / volA, A_recession_pct_of_A=100.0 * a_ph["recession_m3"] / volA,
            A_tail_pct_of_A=100.0 * a_ph["tail_m3"] / volA))
    return pd.DataFrame(rows), dict(t_peak=t_peak, t_tail=t_tail, t_first=qa_all.index[0], t_last=qa_all.index[-1], dt=dtA,
                                    volA=volA, a_ph=a_ph)


# ======================================================================
# Printing
# ======================================================================
def print_table(tab):
    a122.sec("3. The runs  (every change is against truth A, the run with no loss and the true rain)")
    print(f"  {'run':<22}{'kind':<13}{'volume':>9}{'peak':>9}{'peak/vol':>10}{'peak':>8}{'KGE vs A':>10}{'KGE vs B':>10}"
          f"{'vol vs B':>10}{'peak vs B':>11}")
    print(f"  {'':<22}{'':<13}{'chg %':>9}{'chg %':>9}{'':>10}{'min':>8}{'':>10}{'':>10}{'chg %':>10}{'chg %':>11}")
    order = {"control_noloss": 0, "control_loss": 1, "rain": 2, "kchan": 3}
    t = tab.assign(_o=tab["run_kind"].map(order)).sort_values(["_o", "Kchan", "rain_scale"], ascending=[True, True, False])
    for r in t.itertuples():
        print(f"  {r.label:<22}{KIND_LABEL[r.run_kind]:<13}{r.vol_change_vs_A_pct:>+9.2f}{r.peak_change_vs_A_pct:>+9.2f}"
              f"{a122.fnum(r.peak_per_volume, 10, 2)}{r.peak_shift_vs_A_min:>+8.1f}{r.kge_2012_vs_A:>10.4f}{r.kge_2012_vs_B:>10.4f}"
              f"{r.vol_change_vs_B_pct:>+10.2f}{r.peak_change_vs_B_pct:>+11.2f}")
    print("  peak/vol = (percent change of the peak) / (percent change of the volume). 1.0 means the peak falls by the same percent as "
          "the volume; 0.5 would mean it falls half as much.")


def kchan_tables(tab):
    k = tab[tab["curve"] == "loss"].sort_values("Kchan").reset_index(drop=True)
    return k


def print_dose(k, ctlA):
    a122.sec("4. Kchan alone: how much does the outlet lose?  (figure 1)")
    print("  Ks, f, cv, flowexp, n and the rain are the truth's in every run below; only Kchan changes (ctlB is the Kchan-201 run).")
    print(f"  {'Kchan (mm/hr)':>14}{'volume change %':>17}{'peak change %':>15}{'extra volume lost per decade of Kchan':>40}")
    prev = None
    for r in k.itertuples():
        inc = ""
        if prev is not None:
            dec = np.log10(r.Kchan / prev.Kchan)
            inc = f"{-(r.vol_change_vs_A_pct - prev.vol_change_vs_A_pct) / dec:>+.2f} points"
        print(f"  {a122.fk(r.Kchan):>14}{r.vol_change_vs_A_pct:>+17.2f}{r.peak_change_vs_A_pct:>+15.2f}{inc:>40}")
        prev = r
    vol = k["vol_change_vs_A_pct"].to_numpy(float)
    pk = k["peak_change_vs_A_pct"].to_numpy(float)
    mono_v, mono_p = bool(np.all(np.diff(vol) <= 1e-9)), bool(np.all(np.diff(pk) <= 1e-9))
    print(f"\n  Volume falls steadily as Kchan rises:  {'yes' if mono_v else 'NO (it goes back up somewhere; look at the table)'}")
    print(f"  Peak falls steadily as Kchan rises:    {'yes' if mono_p else 'NO (it goes back up somewhere; look at the table)'}")
    span_v, span_p = float(vol.max() - vol.min()), float(pk.max() - pk.min())
    print(f"  Across Kchan {a122.fk(k['Kchan'].min())} to {a122.fk(k['Kchan'].max())} mm/hr the volume changes by {span_v:.2f} percentage points "
          f"and the peak by {span_p:.2f}.")
    if len(vol) >= 3 and abs(vol[-1] - vol[-2]) < 1.0:
        print(f"  The last step ({a122.fk(k['Kchan'].iloc[-2])} to {a122.fk(k['Kchan'].iloc[-1])} mm/hr) changes the volume by only "
              f"{abs(vol[-1] - vol[-2]):.2f} points: the effect is levelling off (the loss cannot take more than the water that arrives).")
    best_b = k.loc[k["kge_2012_vs_B"].idxmax()]
    print(f"  Highest KGE_2012 against truth B among these runs: Kchan {a122.fk(best_b['Kchan'])} mm/hr ({best_b['kge_2012_vs_B']:.4f}); "
          "truth B's own Kchan is 201 (that run scores 1 by construction).")
    off = k[k["run_kind"] != "control_loss"]
    if len(off):
        nxt = off.loc[off["kge_2012_vs_B"].idxmax()]
        print(f"  Best of the OTHER Kchan runs against truth B: Kchan {a122.fk(nxt['Kchan'])} mm/hr, KGE_2012 {nxt['kge_2012_vs_B']:.4f}, "
              f"volume {nxt['vol_change_vs_B_pct']:+.2f}% and peak {nxt['peak_change_vs_B_pct']:+.2f}% against truth B.")
    print("  HOW TO READ IT: a steep, steady curve means Kchan leaves a clear mark on the volume. A flat part means that range of Kchan "
          "is hard to tell apart from its neighbours by volume alone.")


def print_match(tab, match, loss_row, rain_rows):
    a122.sec("5. Loss or rain error?  Same volume, different peak?  (figure 2)")
    print(f"  Truth B's channel loss (ctlB, Kchan 201) removes {-loss_row['vol_change_vs_A_pct']:.2f}% of the volume and "
          f"{-loss_row['peak_change_vs_A_pct']:.2f}% of the peak (peak/volume {loss_row['peak_per_volume']:.2f}), "
          f"and moves the peak by {loss_row['peak_shift_vs_A_min']:+.1f} min.")
    print("  Rain runs (no loss; rain x0.90 ... x0.99 over the whole storm and both gauges):")
    for r in tab[tab["curve"] == "rain"].sort_values("rain_scale").itertuples():
        print(f"    {r.label:<14} volume {r.vol_change_vs_A_pct:+7.2f}%   peak {r.peak_change_vs_A_pct:+7.2f}%   "
              f"peak/volume {a122.fnum(r.peak_per_volume, 5, 2)}   peak shift {r.peak_shift_vs_A_min:+.1f} min")
    if not match["ok"]:
        print(f"\n  NOTE  {match['reason']}.")
        return
    gap = loss_row["peak_change_vs_A_pct"] - match["peak_pct"]
    print(f"\n  Matched volume: a rain error that removes the same {-match['target_vol_pct']:.2f}% of the volume is between "
          f"{match['lower_run']} ({match['lower_vol_pct']:+.2f}%) and {match['upper_run']} ({match['upper_vol_pct']:+.2f}%): "
          f"about rain x{match['rain_scale']:.3f}.")
    print(f"    Interpolated, it cuts the peak by {-match['peak_pct']:.2f}%.  The channel loss cuts it by {-loss_row['peak_change_vs_A_pct']:.2f}%.  "
          f"Difference: {abs(gap):.2f} percentage points ({'the rain error cuts the peak more' if gap > 0 else 'the loss cuts the peak more'}).")
    verdict = abs(gap) >= PEAK_DISTINCT_PP
    print(f"    Rule (a convention fixed in the script, not a statistical test): the peak-and-volume pairing 'differs' if the peak "
          f"change differs by at least {PEAK_DISTINCT_PP:g} percentage points at matched volume.  Here: {'it differs' if verdict else 'it does not differ by that much'}.")
    print(f"  HOW TO READ IT: with the routing at the truth, a hydrograph that has lost {-loss_row['vol_change_vs_A_pct']:.1f}% of its volume "
          "could be a channel loss or a rain error. If the peaks differ at that volume, the peak (and the shape, below) carry "
          "information that separates them. In a real calibration the routing parameters are free and can change the peak too "
          "(Series 121/122), so this separation can be weaker there.")


def print_pair(loss_row, rain_rows, bands, a_share_by_peak):
    a122.sec("6. Timing and shape of the loss next to a rain error of the same volume  (figure 3)")
    print("  Rain run(s) drawn next to the loss (rule: the runs whose volume lies just below and just above the loss's, the nearer first):")
    for i, r in enumerate(rain_rows):
        print(f"    {'nearer' if i == 0 else 'other '}: {r['label']}  volume {r['vol_change_vs_A_pct']:+.2f}% against the loss's "
              f"{loss_row['vol_change_vs_A_pct']:+.2f}%")
    print("  Everything below is against truth A (no loss, true rain), whole 20-hour scoring window.")
    names = [loss_row] + rain_rows
    print(f"  {'run':<20}{'volume':>9}{'peak':>9}{'volume missing':>16}{'gone by A peak':>16}{'after zoom':>12}{'largest gap':>15}")
    print(f"  {'':<20}{'chg %':>9}{'chg %':>9}{'(thousand m3)':>16}{'(% of gap)':>16}{'(% of gap)':>12}{'(m3/s)':>15}")
    for r in names:
        print(f"  {r['label']:<20}{r['vol_change_vs_A_pct']:>+9.2f}{r['peak_change_vs_A_pct']:>+9.2f}"
              f"{r['volume_missing_m3'] / 1e3:>16,.2f}{a122.fnum(r['gap_share_by_A_peak_pct'], 16, 1)}"
              f"{a122.fnum(r['gap_share_after_zoom_pct'], 12, 1)}{r['largest_gap_m3s']:>+15.3f}")
    print("  'gone by A peak' = the share of the missing volume that is already missing when truth A's peak passes. For reference, "
          f"{a_share_by_peak:.1f}% of truth A's OWN volume has passed by its peak.")
    print("\n  Run / truth A at each flow level (A's flow above 5% of its peak), rising limb and recession kept apart "
          "(the same flow level occurs twice in a storm). 0.90 means the run carries 90% of what A carries:")
    print(f"    {'run':<20}{'limb':<11}{'flow band':<24}{'points':>8}{'median':>9}{'lowest':>9}{'highest':>9}")
    for r in bands.itertuples():
        print(f"    {r.run:<20}{r.limb:<11}{r.band:<24}{r.n_points:>8}{a122.fnum(r.median_ratio, 9, 3)}{a122.fnum(r.min_ratio, 9, 3)}{a122.fnum(r.max_ratio, 9, 3)}")
    print("  HOW TO READ THE RATIO (a key, not a verdict):")
    print("    * the same in every band  -> that run takes about the same SHARE of the flow at every flow level;")
    print("    * lower at small flows, nearer 1 at large flows  -> it behaves like a fixed amount taken away, so it removes a bigger share of the small flows;")
    print("    * compare ctlB's pattern with the rain run's, limb by limb: if the patterns differ, shape separates the loss from a rain error.")
    print("    * each limb and band holds only a few 5-minute points (the 'points' column): read the dots in figure 3(d), not just the medians.")
    print("    Why a pattern appears is a hypothesis to check in the tRIBS code; these numbers only show that it does.")


def print_phases(ph, info, tab):
    a122.sec("7. The tail: where does the missing water come from, and does the outlet stop flowing?  (figure 4)")
    t_tail = info["t_tail"]
    dt = pd.Timedelta(seconds=info["dt"])
    print("  Truth A's storm is cut at fixed times (a rule on truth A only), so every run is cut the same way:")
    print(f"    rise        {fmt_t(info['t_first'])} to {fmt_t(info['t_peak'])}   (up to truth A's peak)")
    if t_tail is not None:
        print(f"    recession   {fmt_t(info['t_peak'] + dt)} to {fmt_t(t_tail - dt)}   (after the peak, while A is at or above {TAIL_FRAC * 100:g}% of its peak)")
        print(f"    long tail   {fmt_t(t_tail)} to {fmt_t(info['t_last'])}   (truth A below {TAIL_FRAC * 100:g}% of its peak)")
    else:
        print(f"    recession   {fmt_t(info['t_peak'] + dt)} to {fmt_t(info['t_last'])}   (truth A never falls below {TAIL_FRAC * 100:g}% of its peak here, so there is no long tail)")
    a = info["a_ph"]
    totA = a["rise_m3"] + a["recession_m3"] + a["tail_m3"]
    print(f"  For reference, truth A's own volume ({info['volA'] / 1e3:,.1f} thousand m3) is spread {100 * a['rise_m3'] / totA:.1f}% / "
          f"{100 * a['recession_m3'] / totA:.1f}% / {100 * a['tail_m3'] / totA:.1f}% over rise / recession / long tail.")
    print("\n  Share of each run's missing volume that falls in each phase, and whether its outlet stops flowing:")
    print(f"  {'run':<22}{'missing':>9}{'rise':>8}{'recession':>11}{'long tail':>11}   {'outlet dry from':<17}{'hours dry':>10}{'A carries then':>16}")
    print(f"  {'':<22}{'(1000 m3)':>9}{'(% gap)':>8}{'(% gap)':>11}{'(% gap)':>11}   {'':<17}{'':>10}{'(% of A vol)':>16}")
    order = {"control_noloss": 0, "control_loss": 1, "rain": 2, "kchan": 3}
    t = ph.assign(_o=ph["run_kind"].map(order), _k=tab.set_index("run_name").loc[ph["run_name"], "Kchan"].to_numpy(),
                  _r=tab.set_index("run_name").loc[ph["run_name"], "rain_scale"].to_numpy()).sort_values(["_o", "_k", "_r"], ascending=[True, True, False])
    for r in t.itertuples():
        dry = r.dry_from if r.dry_from else "never"
        if r.dry_from:
            dry = fmt_t(pd.Timestamp(r.dry_from))
        print(f"  {r.label:<22}{r.total_missing_m3 / 1e3:>9.2f}{a122.fnum(r.rise_share_of_gap_pct, 8, 1)}{a122.fnum(r.recession_share_of_gap_pct, 11, 1)}"
              f"{a122.fnum(r.tail_share_of_gap_pct, 11, 1)}   {dry:<17}{r.dry_hours:>10.1f}{r.A_volume_in_dry_period_pct:>16.2f}")
    print("  'outlet dry from' = the first time from which the run's outlet flow stays at zero to the end of the window ('never' = it keeps flowing). "
          "'A carries then' = the water truth A carries during that time, as a share of A's volume: a run that is dry carries none of it.")
    print("  Shares can be negative: a run that carries a little MORE water than truth A in a phase has a negative share there.")
    print("  HOW TO READ IT (a key, not a verdict): compare the long-tail share and the dry time of the loss with the rain runs'. "
          "Why the outlet runs dry (or not) is a hypothesis to check in the tRIBS code; these numbers only show that it does.")
    print("  CAUTION: the long tail is very small flow (see figure 4's scale). Whether a real gauge could resolve it has not been checked.")


# ======================================================================
# Figures
# ======================================================================
def figure_1(k, ctlA_row, out_path, dpi):
    plt = a122._plt()
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.8))
    fig.patch.set_facecolor("white")
    for ax in axes:
        a122.style_axes(ax)
    K = k["Kchan"].to_numpy(float)
    isb = (k["run_kind"] == "control_loss").to_numpy()
    ax = axes[0]
    ax.axhline(0, color=INK_SOFT, lw=1.0)
    ax.plot(K, k["vol_change_vs_A_pct"], color=BLUE, lw=2.0, marker="o", ms=6.5, label="volume")
    ax.plot(K, k["peak_change_vs_A_pct"], color=BLUE, lw=1.8, linestyle=(0, (4, 2)), marker="s", ms=6.0, mfc="white", mew=1.6, label="peak")
    if isb.any():
        ax.plot(K[isb], k["vol_change_vs_A_pct"][isb], marker="o", ms=13, mfc="none", mec=INK, mew=1.4, linestyle="")
        ax.annotate("Kchan 201 = truth B", (K[isb][0], k["vol_change_vs_A_pct"][isb].iloc[0]), xytext=(12, -22),
                    textcoords="offset points", fontsize=8.8, color=INK)
    ax.set_xscale("log")
    ax.set_xlabel("Kchan (mm/hr, log scale)   the only thing that changes between these runs", fontsize=9.5, color=INK_SOFT)
    ax.set_ylabel("change against truth A, no loss (%)", fontsize=9.5, color=INK_SOFT)
    ax.set_title("(a) How much volume and peak does the channel loss take?", fontsize=10.5, color=INK, loc="left")
    ax.legend(loc="lower left", frameon=False, fontsize=9)
    ax = axes[1]
    ax.plot(K, k["kge_2012_vs_B"], color=BLUE, lw=2.0, marker="o", ms=6.5, label="against truth B (loss, Kchan 201)")
    ax.plot(K, k["kge_2012_vs_A"], color=INK_SOFT, lw=1.8, linestyle=(0, (4, 2)), marker="s", ms=6.0, mfc="white", mew=1.6,
            label="against truth A (no loss)")
    ax.set_xscale("log")
    ax.set_xlabel("Kchan (mm/hr, log scale)", fontsize=9.5, color=INK_SOFT)
    ax.set_ylabel("KGE_2012 (1 = identical)", fontsize=9.5, color=INK_SOFT)
    ax.set_title("(b) How well does each run match each truth?", fontsize=10.5, color=INK, loc="left")
    ax.legend(loc="lower center", frameon=False, fontsize=9)
    fig.suptitle("Series 122 clean runs: what does Kchan alone do?  (Ks, f, cv, flowexp, n and rain at the truth's values)",
                 fontsize=12.5, color=INK, y=0.995)
    fig.text(0.5, 0.005, "One storm (12 Aug 2014); model output, not observations. Each dot is one tRIBS run. The 20-hour scoring window is used for every number.",
             ha="center", va="bottom", fontsize=8.2, color=INK_SOFT)
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))
    a122.save_fig(fig, out_path, dpi)


def place_labels(fig, ax, items, points_xy):
    """Put text labels beside points without letting them cover each other or the points. items: list of (text, (x, y), candidates);
    a candidate is (dx, dy, ha) in points from the point. Each label takes the first candidate whose box is clear of the labels
    already placed and of every point in points_xy (data coordinates); if none is clear, the first candidate is used.
    Call after the axes limits and the layout are final (it measures the boxes on screen)."""
    from matplotlib.transforms import Bbox
    try:
        renderer = fig.canvas.get_renderer()
    except Exception:                                                                # noqa: BLE001  fall back to the first candidate
        renderer = None
    box = dict(boxstyle="round,pad=0.15", fc=a122.SURFACE, ec="none", alpha=0.85)
    placed = []
    if renderer is not None:
        for x, y in points_xy:
            px, py = ax.transData.transform((x, y))
            placed.append(Bbox([[px - 6, py - 6], [px + 6, py + 6]]))
    for text, xy, cands in items:
        chosen = None
        for dx, dy, ha in cands:
            ann = ax.annotate(text, xy, xytext=(dx, dy), textcoords="offset points", ha=ha, fontsize=8.2, color=INK, zorder=5, bbox=box)
            if renderer is None:
                chosen = ann
                break
            bb = ann.get_window_extent(renderer)
            hit = any(bb.overlaps(p) for p in placed)
            inside = ax.get_window_extent(renderer).contains(bb.x0, bb.y0) and ax.get_window_extent(renderer).contains(bb.x1, bb.y1)
            if not hit and inside:
                placed.append(bb)
                chosen = ann
                break
            ann.remove()
        if chosen is None:
            dx, dy, ha = cands[0]
            chosen = ax.annotate(text, xy, xytext=(dx, dy), textcoords="offset points", ha=ha, fontsize=8.2, color=INK, zorder=5, bbox=box)
            if renderer is not None:
                placed.append(chosen.get_window_extent(renderer))


def figure_2(tab, match, loss_row, out_path, dpi):
    plt = a122._plt()
    fig, ax = plt.subplots(figsize=(9.6, 7.2))
    fig.patch.set_facecolor("white")
    a122.style_axes(ax)
    loss = tab[tab["curve"].isin(["loss", "anchor"])].copy()
    rain = tab[tab["curve"].isin(["rain", "anchor"])].copy()
    loss["x"], loss["y"] = -loss["vol_change_vs_A_pct"], -loss["peak_change_vs_A_pct"]
    rain["x"], rain["y"] = -rain["vol_change_vs_A_pct"], -rain["peak_change_vs_A_pct"]
    loss, rain = loss.sort_values("x"), rain.sort_values("x")
    top = float(max(loss["x"].max(), rain["x"].max(), loss["y"].max(), rain["y"].max())) * 1.12
    ax.plot([0, top], [0, top], color=INK_SOFT, lw=1.0, linestyle=(0, (2, 3)))
    ax.text(0.98, 0.04, "dotted diagonal: the peak falls by the same percent as the volume", transform=ax.transAxes,
            fontsize=8.5, color=INK_SOFT, ha="right", va="bottom")
    ax.plot(rain["x"], rain["y"], color=RAIN_COL, lw=2.0, marker="o", ms=6.5, label="rain error (no loss; rain scaled down)")
    ax.plot(loss["x"], loss["y"], color=LOSS_COL, lw=2.0, marker="s", ms=6.5, label="channel loss (Kchan only; rain at the truth)")
    ax.plot([-loss_row["vol_change_vs_A_pct"]], [-loss_row["peak_change_vs_A_pct"]], marker="s", ms=13, mfc="none", mec=INK, mew=1.4, linestyle="")
    if match["ok"]:
        x0 = -match["target_vol_pct"]
        ax.axvline(x0, color=INK_SOFT, lw=1.0, linestyle=(0, (2, 3)))
        ax.plot([x0], [-match["peak_pct"]], marker="o", ms=13, mfc="none", mec=INK, mew=1.4, linestyle="")
        ax.annotate("", xy=(x0, -loss_row["peak_change_vs_A_pct"]), xytext=(x0, -match["peak_pct"]),
                    arrowprops=dict(arrowstyle="<->", color=INK, lw=1.2))
        ax.text(x0 + top * 0.015, (-match["peak_pct"] - loss_row["peak_change_vs_A_pct"]) / 2.0,
                f"at the same {x0:.1f}% volume lost:\nchannel loss cuts the peak by {-loss_row['peak_change_vs_A_pct']:.1f}%\n"
                f"rain error cuts it by {-match['peak_pct']:.1f}% (interpolated)", fontsize=9, color=INK, va="center", ha="left")
    ax.set_xlim(-0.02 * top, top)
    ax.set_ylim(-0.05 * top, top)
    ax.set_xlabel("volume lost compared with truth A (%)", fontsize=9.8, color=INK_SOFT)
    ax.set_ylabel("peak discharge lost compared with truth A (%)", fontsize=9.8, color=INK_SOFT)
    ax.set_title("Same volume lost: does it cost the same peak?", fontsize=11, color=INK, loc="left")
    ax.legend(loc="upper left", frameon=False, fontsize=9)
    fig.suptitle("Series 122 clean runs: channel loss against rain error", fontsize=12.5, color=INK, y=0.995)
    cap2 = ("Circled square = truth B's channel loss (Kchan 201). Circle on the dotted line = the rain error at the same volume, "
            "a straight-line estimate between two rain runs." if match["ok"] else
            "Circled square = truth B's channel loss (Kchan 201). Its volume is outside the range the rain runs cover, "
            "so no rain point is interpolated.")
    fig.text(0.5, 0.005, "Routing (cv, flowexp, n), Ks and f at the truth in every run. One storm (12 Aug 2014); model output, not observations.\n" + cap2,
             ha="center", va="bottom", fontsize=8.2, color=INK_SOFT, linespacing=1.5)
    fig.tight_layout(rect=(0, 0.06, 1, 0.97))
    # labels last: they are placed by measuring the final layout, so none covers another label or a point
    pts = [(r["x"], r["y"]) for _, r in pd.concat([loss, rain]).iterrows()]
    items = []
    for _, r in loss[loss["curve"] == "loss"].iterrows():
        items.append((r["label"].replace("ctlB ", ""), (r["x"], r["y"]),
                      [(8, -12, "left"), (8, 8, "left"), (-8, -12, "right"), (-8, 8, "right"), (10, -26, "left"), (10, 22, "left")]))
    for _, r in rain[rain["curve"] == "rain"].iterrows():
        items.append((r["label"], (r["x"], r["y"]),
                      [(-8, 8, "right"), (-8, -14, "right"), (8, 8, "left"), (8, -14, "left"), (-8, 22, "right")]))
    place_labels(fig, ax, items, pts)
    a122.save_fig(fig, out_path, dpi)


def figure_3(fr, loss_row, rain_rows, bands_pts, out_path, dpi, view, zoom_note, loss_after_zoom, loss_dry=None):
    plt = a122._plt()
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
    mk = h122.MARK if view["zoomed"] else {}
    fig, axes = plt.subplots(2, 2, figsize=(13.4, 9.2))
    fig.patch.set_facecolor("white")
    for ax in axes.ravel():
        a122.style_axes(ax)
    dt = (fr.index[1] - fr.index[0]).total_seconds()
    volA = float(fr["A"].sum() * dt)
    series = [("ctlB", loss_row, LOSS_COL, 1.0, 2.2)]
    for i, r in enumerate(rain_rows):
        series.append((f"rain{i}", r, RAIN_COL, 1.0 if i == 0 else 0.55, 2.0 if i == 0 else 1.7))
    axA, axR, axC, axD = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]
    axA.plot(fr.index, fr["A"], color=INK, lw=1.5, linestyle=(0, (3, 2)), label="truth A (no loss, true rain)")
    shown_q, shown_res, shown_gap = [fr["A"]], [], []
    for col, r, color, alpha, lw in series:
        lab = f"{r['label']}:  volume {r['vol_change_vs_A_pct']:+.1f}%, peak {r['peak_change_vs_A_pct']:+.1f}%"
        axA.plot(fr.index, fr[col], color=color, lw=lw, alpha=alpha, label=lab, **mk)
        res = fr[col] - fr["A"]
        gap = 100 * (fr["A"].cumsum() - fr[col].cumsum()) * dt / volA
        axR.plot(fr.index, res, color=color, lw=lw, alpha=alpha, **mk)
        axC.plot(fr.index, gap, color=color, lw=lw, alpha=alpha, **mk)
        shown_q.append(fr[col]), shown_res.append(res), shown_gap.append(gap)
    axA.set_ylabel("discharge (m3/s)", fontsize=9.5, color=INK_SOFT)
    axA.set_title("(a) The hydrographs", fontsize=10.5, color=INK, loc="left")
    h122.time_axis(axA, view)
    h122.fit_y(axA, shown_q, view, floor0=True, pad=0.34)
    axA.legend(loc="upper right", frameon=False, fontsize=8.6)
    axR.axhline(0, color=INK_SOFT, lw=1.0)
    axR.set_ylabel("run minus truth A (m3/s)", fontsize=9.5, color=INK_SOFT)
    axR.set_title("(b) When is water missing?  (below zero = the run has less)", fontsize=10.5, color=INK, loc="left")
    h122.time_axis(axR, view)
    h122.fit_y(axR, shown_res, view, include=[0.0])
    axC.axhline(0, color=INK_SOFT, lw=1.0)
    axC.set_ylabel("running volume missing (% of truth A's volume)", fontsize=9.5, color=INK_SOFT)
    axC.set_title("(c) How the missing volume builds up", fontsize=10.5, color=INK, loc="left")
    h122.time_axis(axC, view)
    h122.fit_y(axC, shown_gap, view, include=[0.0])
    if view["zoomed"] and np.isfinite(loss_after_zoom):
        note = f"{loss_after_zoom:.0f}% of ctlB's missing volume comes after this window ends"
        if loss_dry is not None and loss_dry["dry_from"] is not None:
            note += f"\n(ctlB's outlet is dry from {loss_dry['dry_from']:%H:%M}; see figure 4)"
        axC.text(0.98, 0.12, note, transform=axC.transAxes, fontsize=8.8, color=INK, va="bottom", ha="right")
    for col, r, color, alpha, lw in series:
        s = bands_pts[col]
        up = s["rising"].to_numpy(bool)
        a_pt = 0.9 if alpha == 1.0 else 0.6
        axD.scatter(s["A"][~up], s["ratio"][~up], s=22, color=color, alpha=a_pt, linewidths=0, zorder=3)                # recession: filled
        axD.scatter(s["A"][up], s["ratio"][up], s=30, facecolors="none", edgecolors=color, alpha=a_pt, linewidths=1.4, zorder=3)   # rising: open
    axD.set_xscale("log")
    axD.xaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
    axD.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    axD.xaxis.set_minor_formatter(NullFormatter())
    allr = np.concatenate([bands_pts[c]["ratio"].to_numpy(float) for c, *_ in series])
    allq = np.concatenate([bands_pts[c]["A"].to_numpy(float) for c, *_ in series])
    if len(allr):
        lo, hi = float(allr.min()), float(allr.max())
        mid, half = (lo + hi) / 2.0, max((hi - lo) / 2.0 * 1.10, 0.05)
        axD.set_ylim(mid - half, mid + half)
        axD.set_xlim(float(allq.min()) * 0.9, float(allq.max()) * 1.1)
    axD.ticklabel_format(axis="y", useOffset=False)
    axD.set_xlabel("truth A's discharge at that time (m3/s, log scale)\nopen dots: rising limb (up to A's peak); filled dots: recession (after it)\n"
                   "every point above 5% of A's peak in the whole 20-hour window; not zoomed", fontsize=8.6, color=INK_SOFT)
    axD.set_ylabel("ratio  run / truth A", fontsize=9.5, color=INK_SOFT)
    axD.set_title("(d) Is the same share of the flow missing at every flow level?", fontsize=10.5, color=INK, loc="left")
    handles = [Line2D([0], [0], marker="o", linestyle="", color=c, alpha=a, markersize=6, label=r["label"]) for _, r, c, a, _ in series]
    axD.legend(handles=handles, loc="lower right", frameon=False, fontsize=8.6)
    fig.suptitle(h122.title_with_view("Series 122 clean runs: the channel loss next to a rain error of about the same volume",
                                      view, zoom_note), fontsize=12.5, color=INK, y=0.995)
    fig.text(0.5, 0.005, "Routing, Ks and f at the truth in every run; the only differences are the channel loss (ctlB) or the rain scale. "
             "One storm (12 Aug 2014); model output, not observations."
             + ("  A dot is one 5-minute point." if view["zoomed"] else ""), ha="center", va="bottom", fontsize=8.2, color=INK_SOFT)
    fig.tight_layout(rect=(0, 0.03, 1, 0.985))
    a122.save_fig(fig, out_path, dpi)


def figure_4(fr, loss_row, rain_rows, info, drys, out_path, dpi):
    """The tail after the flood on a log scale (zero cannot be drawn there, so a run's line stops where its flow is zero and a
    thick bar along the bottom marks the time its outlet is dry)."""
    plt = a122._plt()
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
    fig, axes = plt.subplots(1, 2, figsize=(14.0, 6.0), gridspec_kw=dict(width_ratios=[1.35, 1.0]))
    fig.patch.set_facecolor("white")
    for ax in axes:
        a122.style_axes(ax)
    t_peak, t_tail = info["t_peak"], info["t_tail"]
    peakA = float(fr["A"].max())
    t10 = fr["A"][(fr.index > t_peak) & (fr["A"] <= TAIL_ZOOM_FROM * peakA)]
    t10 = t10.index[0] if len(t10) else t_peak
    views = [(t_peak, fr.index[-1]), (t10, min(t10 + pd.Timedelta(hours=TAIL_ZOOM_HOURS), fr.index[-1]))]
    series = [("ctlB", loss_row, LOSS_COL, 1.0, 2.2)]
    for i, r in enumerate(rain_rows):
        series.append((f"rain{i}", r, RAIN_COL, 1.0 if i == 0 else 0.55, 2.0 if i == 0 else 1.7))
    pos = [float(v) for c in ["A"] + [s[0] for s in series] for v in fr[c][fr[c] > DRY_FLOW]]
    for (lo, hi), ax, panel in zip(views, axes, ("a", "b")):
        sub = fr.loc[lo:hi]
        vis = [float(v) for c in ["A"] + [s[0] for s in series] for v in sub[c][sub[c] > DRY_FLOW]] or pos
        ymin, ymax = min(vis) / 4.0, max(vis) * 1.6
        ax.set_yscale("log")
        ax.plot(fr.index, fr["A"].where(fr["A"] > DRY_FLOW), color=INK, lw=1.5, linestyle=(0, (3, 2)), label="truth A (no loss, true rain)",
                **(h122.MARK if panel == "b" else {}))
        for k, (col, r, color, alpha, lw) in enumerate(series):
            ax.plot(fr.index, fr[col].where(fr[col] > DRY_FLOW), color=color, lw=lw, alpha=alpha,
                    label=f"{r['label']}:  volume {r['vol_change_vs_A_pct']:+.1f}%, peak {r['peak_change_vs_A_pct']:+.1f}%",
                    **(h122.MARK if panel == "b" else {}))
        ax.set_ylim(ymin, ymax)
        h122.time_axis(ax, h122.view_of(fr, (lo, hi)))
        ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1.0,)))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.set_ylabel("discharge (m3/s, log scale: each gridline is 10 times the one below)", fontsize=8.8, color=INK_SOFT)
        # dry bars along the bottom, one level per run that goes dry
        level = ymin * 1.5
        for col, r, color, alpha, lw in series:
            d = drys[col]
            if d["dry_from"] is None:
                continue
            x0 = max(d["dry_from"], lo)
            if x0 <= hi:
                ax.plot([x0, hi], [level, level], color=color, lw=5, solid_capstyle="butt", alpha=alpha, zorder=4)
                ax.annotate(f"{r['label']}: outlet flow is 0 from {d['dry_from']:%H:%M} on {d['dry_from']:%d %b} to the end of the scoring window "
                            f"({d['dry_hours']:.1f} h)" if panel == "a" else f"{r['label']}: outlet dry from {d['dry_from']:%H:%M}",
                            (hi, level), xytext=(0, 9), textcoords="offset points", ha="right", fontsize=8.4, color=INK, zorder=5,
                            bbox=dict(boxstyle="round,pad=0.15", fc=a122.SURFACE, ec="none", alpha=0.85))
            level *= 1.7
        if t_tail is not None and lo < t_tail < hi:
            ax.axvline(t_tail, color=INK_SOFT, lw=1.0, linestyle=(0, (2, 3)))
    axes[0].set_title("(a) After the peak, to the end of the scoring window", fontsize=10.5, color=INK, loc="left")
    axes[1].set_title(f"(b) The first {TAIL_ZOOM_HOURS:g} hours from when truth A is down to {TAIL_ZOOM_FROM * 100:g}% of its peak",
                      fontsize=10.5, color=INK, loc="left")
    handles, labels = axes[0].get_legend_handles_labels()
    if t_tail is not None:
        handles.append(Line2D([0], [0], color=INK_SOFT, lw=1.0, linestyle=(0, (2, 3))))
        labels.append(f"dotted line: truth A is below {TAIL_FRAC * 100:g}% of its peak from here on")
    axes[0].legend(handles, labels, loc="upper right", frameon=False, fontsize=8.6)
    fig.suptitle("Series 122 clean runs: the tail after the flood (log scale)", fontsize=12.5, color=INK, y=0.995)
    fig.text(0.5, 0.005, "A log scale cannot show zero: a line stops where that run's flow is 0, and the thick bar along the bottom marks the time its outlet is dry. "
             "Routing, Ks and f at the truth in every run. One storm (12 Aug 2014); model output, not observations.\n"
             "Flows below about 0.05 m3/s are tiny next to the peak; whether a real gauge could resolve them has not been checked.",
             ha="center", va="bottom", fontsize=8.2, color=INK_SOFT, linespacing=1.5)
    fig.tight_layout(rect=(0, 0.06, 1, 0.97))
    a122.save_fig(fig, out_path, dpi)


# ======================================================================
# Plumbing
# ======================================================================
def parse_args(argv=None):
    ap = argparse.ArgumentParser(description="Series 122 clean runs: read the twelve one-at-a-time runs and draw what Kchan alone "
                                             "does, and how a channel loss differs from a rain error (reads results; writes a new folder).")
    ap.add_argument("--dry_run", action="store_true")
    ap.add_argument("--window", choices=("zoom", "full", "both"), default="both", help="which versions of figure 3 to draw (default both)")
    ap.add_argument("--zoom_start", default=None, help="start of the zoom window, e.g. \"2014-08-12 17:00\" (default: the rule)")
    ap.add_argument("--zoom_end", default=None, help="end of the zoom window, e.g. \"2014-08-12 21:30\" (default: the rule)")
    ap.add_argument("--results", type=Path, default=None)
    ap.add_argument("--out_dir", type=Path, default=None)
    ap.add_argument("--no_plots", action="store_true")
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args(argv)
    if not (72 <= args.dpi <= 600):
        ap.error("--dpi must be between 72 and 600.")
    return args


def print_dry_run(pre, args):
    a122.sec("DRY RUN: nothing was drawn or written")
    print(f"  Event window {pre['ev_start']} to {pre['ev_end']}.  {len(pre['df'])} runs read and re-scored.")
    print(f"  Output folder (created only on a real run): {pre['out_dir']}")
    names = list(FILES.values())
    if args.window == "zoom":
        names = [n for n in names if n != FILES["fig3"]]
    elif args.window == "full":
        names = [n for n in names if n != FILES["fig3z"]]
    print("  Files that would be written there: " + ", ".join(names))
    print("\nReady:  python -u plot_series122_clean_runs.py 2>&1 | tee series122_clean_compare.log")


def run(args, tee):
    paths = locate(args)
    pre = preflight(paths)
    gA, gB = a122.qout_to_grid(pre["ta_path"]), a122.qout_to_grid(pre["tb_path"])
    grids = {rid: a122.qout_to_grid(p) for rid, p in pre["kept_paths"].items()}
    scored = rescore_and_check(pre, gA, gB, grids)

    qa_event = window_frame({"A": gA}, pre["ev_start"], pre["ev_end"])["A"]
    zoom, zinfo = None, None
    if args.window != "full":
        zoom, zinfo = h122.resolve_zoom(args, qa_event)
    elif args.zoom_start or args.zoom_end:
        print("\n  NOTE  --zoom_start / --zoom_end are ignored with --window full.")
    h122.print_zoom(zoom, zinfo, qa_event, args.window)
    if args.dry_run:
        print_dry_run(pre, args)
        return 0
    manual = bool(zinfo) and zinfo["source"] != "the rule (Truth A's volume)"
    zoom_note = ("window set by you" if manual else
                 f"window from the rule: {h122.ZOOM_PAD_MIN} min before {h122.ZOOM_LO_FRAC * 100:g}% to {h122.ZOOM_PAD_MIN} min after "
                 f"{h122.ZOOM_HI_FRAC * 100:g}% of truth A's volume")

    t0 = time.time()
    tab = build_table(pre, gA, grids, scored, zoom)
    print_table(tab)

    kt = tab[tab["curve"] == "loss"].sort_values("Kchan").reset_index(drop=True)
    ctlA_row = tab[tab["run_kind"] == "control_noloss"].iloc[0]
    loss_row = tab[tab["run_kind"] == "control_loss"].iloc[0]
    print_dose(kt, ctlA_row)

    rain_only = tab[tab["run_kind"] == "rain"].reset_index(drop=True)
    rain_line = pd.concat([rain_only, tab[tab["run_kind"] == "control_noloss"]], ignore_index=True)
    match = match_volume(float(loss_row["vol_change_vs_A_pct"]), rain_line)
    print_match(tab, match, loss_row, rain_only)

    # ---- the pair for figure 3 (printed in section 6)
    rain_rows = pick_rain_runs(float(loss_row["vol_change_vs_A_pct"]), rain_only)
    fr_pair = window_frame({"A": gA, "ctlB": grids[loss_row["run_id"]],
                            **{f"rain{i}": grids[r["run_id"]] for i, r in enumerate(rain_rows)}}, pre["ev_start"], pre["ev_end"])
    cols = ["ctlB"] + [f"rain{i}" for i in range(len(rain_rows))]
    names = {"ctlB": loss_row["label"], **{f"rain{i}": r["label"] for i, r in enumerate(rain_rows)}}
    bands = ratio_bands({names[c]: (fr_pair, c) for c in cols})
    bands_pts = {}
    peakA = float(fr_pair["A"].max())
    t_peakA = fr_pair["A"].idxmax()
    sub = fr_pair[fr_pair["A"] >= h122.FLOOR_FRAC * peakA]
    for c in cols:
        bands_pts[c] = pd.DataFrame({"A": sub["A"].to_numpy(), "ratio": (sub[c] / sub["A"]).to_numpy(),
                                     "rising": np.asarray(sub.index <= t_peakA)}, index=sub.index)
    a_share_by_peak = 100.0 * float(fr_pair["A"].loc[:t_peakA].sum()) / float(fr_pair["A"].sum())
    print_pair(loss_row, rain_rows, bands, a_share_by_peak)

    # ---- the tail: phases for every run, and does the outlet run dry (section 7, figure 4)
    phases, pinfo = build_phases(pre, gA, grids, tab)
    print_phases(phases, pinfo, tab)
    drys = {c: dry_info(fr_pair["A"], fr_pair[c], pinfo["t_peak"]) for c in cols}

    a122.sec("8. CAUTIONS")
    for line in ("One synthetic truth pair, one storm (12 Aug 2014), Ks 7.0 and f 0.012 fixed, no noise. Model output, not observations.",
                 "Routing is held at the truth here. In a real calibration it is free and can change the peak, so a difference found at "
                 "fixed routing may be weaker, or absorbed, in practice.",
                 "Rain is scaled by one factor over the whole storm and both gauges; a real rain error varies in time and space.",
                 "The matched-volume rain peak is a straight-line estimate between two rain runs, not a model run.",
                 "The kept files are outlet-only: they show what leaves the basin, not where in the channel network the water is lost.",
                 "The zoom window follows a fixed rule on truth A's volume; the 'after zoom' column says how much of each run's missing volume falls outside it.",
                 "The phases (rise, recession, long tail) are cut by a rule on truth A only; 'dry' means a flow of 1e-6 m3/s or less. A different storm, "
                 "or a different Kchan, would give different times.",
                 "The long tail is very small flow. The loss and the rain error may differ most there, and a real gauge is usually least accurate there; "
                 "this experiment has no gauge error and the gauge's resolution has not been checked.",
                 "Why a pattern appears is a hypothesis to check in the tRIBS code, not something these curves show."):
        print("  * " + line)

    # ---- write ----
    a122.sec("Files written to " + str(pre["out_dir"]))
    out = pre["out_dir"]
    out.mkdir(parents=True, exist_ok=True)
    tab_out = tab.copy()
    tab_out["largest_gap_time"] = tab_out["largest_gap_time"].map(lambda v: v.strftime("%Y-%m-%d %H:%M") if isinstance(v, pd.Timestamp) else v)
    match_row = {k: v for k, v in match.items()}
    match_row.update(loss_run=loss_row["label"], loss_vol_pct=float(loss_row["vol_change_vs_A_pct"]),
                     loss_peak_pct=float(loss_row["peak_change_vs_A_pct"]),
                     peak_gap_points=(float(loss_row["peak_change_vs_A_pct"] - match["peak_pct"]) if match["ok"] else float("nan")),
                     pairing_differs_rule_points=PEAK_DISTINCT_PP, rain_runs_drawn=" | ".join(r["label"] for r in rain_rows))
    ts = pd.DataFrame({"time": fr_pair.index, "truth_A_m3s": fr_pair["A"].to_numpy(), "ctlB_m3s": fr_pair["ctlB"].to_numpy()})
    for i, r in enumerate(rain_rows):
        ts[f"{r['run_name']}_m3s"] = fr_pair[f"rain{i}"].to_numpy()
    for name, t in (("table", tab_out), ("matched", pd.DataFrame([match_row])), ("bands", bands), ("phases", phases), ("ts", ts)):
        t.to_csv(out / FILES[name], index=False, float_format="%.8g")
        print(f"  saved {FILES[name]}  ({len(t)} rows)")

    if not args.no_plots:
        print("\nDrawing the figures...")
        figs = [("fig1", lambda: figure_1(kt, ctlA_row, out / FILES["fig1"], args.dpi)),
                ("fig2", lambda: figure_2(tab, match, loss_row, out / FILES["fig2"], args.dpi))]
        after = float(loss_row["gap_share_after_zoom_pct"])
        if zoom is not None:
            figs.append(("fig3z", lambda: figure_3(fr_pair, loss_row, rain_rows, bands_pts, out / FILES["fig3z"], args.dpi,
                                                   h122.view_of(fr_pair, zoom), zoom_note, after, drys["ctlB"])))
        if args.window in ("full", "both"):
            figs.append(("fig3", lambda: figure_3(fr_pair, loss_row, rain_rows, bands_pts, out / FILES["fig3"], args.dpi,
                                                  h122.view_of(fr_pair, None), zoom_note, after, drys["ctlB"])))
        figs.append(("fig4", lambda: figure_4(fr_pair, loss_row, rain_rows, pinfo, drys, out / FILES["fig4"], args.dpi)))
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
    print(f"Series 122 clean runs: comparison   python {sys.version.split()[0]}, pandas {pd.__version__}, numpy {np.__version__}   "
          f"{datetime.now().strftime('%Y-%m-%d %H:%M')}")
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
