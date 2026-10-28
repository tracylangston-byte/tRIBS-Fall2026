"""
analyze_series122_kchan.py
==========================
Series 122 analysis ("Scenario 5"): can the channel-loss parameter Kchan be recovered, and is it falsely used?

The Series 122 sweep (run_lhs_kchan_122_tribs6.py) ran 400 tRIBS simulations with channel loss ON, Ks 7.0 and f 0.012
fixed, and four parameters varied together: cv (kinemvelcoef), flowexp, n (channelroughness) and Kchan
(channelconductivity_mmhr, 3 to 1000 mm/hr, log scale). The sweep scored every run against truth A (the no-loss
truth). This script scores the same 400 runs against truth B too (the Kchan = 201 mm/hr truth), from the outlet
files the sweep kept. It runs NO tRIBS, and it changes none of your files.

TWO QUESTIONS
-------------
  A. RECOVERY  (scored against truth B, which HAS channel loss, true Kchan 201):
       with Kchan on and routing free, do the best runs recover the true Kchan?
  B. FALSE USE (scored against truth A, which has NO channel loss):
       how much Kchan do good runs tolerate, and do the best runs pile up at the lowest Kchan?

THE READING RULES  (fixed before any result was seen; they are constants in this script, not options)
---------------------------------------------------------------------------------------------------
  RECOVERY (truth B). Recovered if BOTH hold:
      (1) the true Kchan (201) lies inside the middle 80% of Kchan among the best 20% of runs (by KGE_2012), and
      (2) that middle-80% span covers less than half of the log box. The box runs 3 to 1000 mm/hr, which is
          2.52 decades, so the span must be narrower than 1.26 decades (a factor of about 18 end to end).
  FALSE USE (truth A). Two numbers are reported: the largest Kchan among good runs (KGE_2012 of 0.95 or more),
      and the share of the best 20% with Kchan above 30 mm/hr. The pass/fail version: false use if more than 10%
      of the best 20% have Kchan above 30.
  KCHAN ADDS A DIRECTION. Yes if the held-out R squared of "the dial plus a curve in log Kchan" beats "the dial
      alone" by more than 0.02 (in Series 121 the largest model-form gap was 0.007). The handoff does not say
      which truth this is judged against, so it is reported against BOTH and each is labelled.
  NULL CONTROL. Random log-uniform Kchan values are attached to the 100 f0p012_true runs of
      lhs_results_anchor_ALL_121.csv (there channel loss was off, so Kchan did nothing). The same rules are applied
      to them, so you can see what "Kchan has no effect" looks like with these rules and this scoring.

Everything else in the output (quantile tables, volume-versus-shape split, partial correlations, the model ladder,
the shape-beyond-the-dial test) is INFORMATION to help read the rules. It is labelled as such and is not a rule.

WHAT IS THE "DIAL"?
-------------------
  Series 121 found that cv, flowexp and n act as ONE combined direction: score = g(w1*cv + w2*flowexp + w3*n),
  with g a smooth curve (here a cubic). Held-out R squared is the share of the variance in runs LEFT OUT of a fit
  that the model explains (repeated 5-fold cross-validation), so a model with more terms only wins if it truly
  captures more. "Dial plus a curve in log Kchan" is  g(dial) + h(log10 Kchan)  with h a cubic.

HOW THE SCORING IS DONE (so the truth-B scores can be trusted)
--------------------------------------------------------------
  The sweep's scorer resamples each outlet series to a 5-minute grid with time-interpolation, clips it to the event
  window, and computes KGE_2012 and the other metrics. This script repeats exactly that on the kept outlet files.
  Before it uses the result it checks itself in four ways, and STOPS if any fails:
     1. it re-scores every run against truth A and must reproduce the KGE_2012, PBIAS, r, beta, gamma and the
        phase metrics that are stored in lhs_results_kchan_122.csv (to 1e-6);
     2. it scores truth B against truth A and must reproduce what Series 110 recorded (PBIAS -9.98684%,
        KGE_2012 0.88414);
     3. every kept outlet file's checksum must equal the checksum the sweep recorded for it;
     4. truth A and truth B must be the exact files (checksums) that Series 122 and the truth-B check used.

WHAT IT READS  (only reads)
---------------------------
    calibration_work/03_comparisons/summary_tables/lhs_results_kchan_122.csv      the sweep's results (400 rows)
    calibration_work/kept_qout_122/<run_id>_Outlet.qout                            the 400 kept outlet files
    calibration_work/synth_truth/*.qout                                            truth A (the one file there)
    calibration_work/synth_truth_location/SMF_20140812_110_cc201p0_loc_K7p0_f0p012_Outlet.qout    truth B
    calibration_work/03_comparisons/summary_tables/lhs_results_anchor_ALL_121.csv  for the null control

WHAT IT WRITES  (a NEW folder; it refuses to write into the folder holding the results file)
-------------------------------------------------------------------------------------------
  calibration_work/03_comparisons/sensitivity_plots/series122_kchan/
    series122_report.txt               everything printed on screen
    series122_scored_both_truths.csv   the 400 runs with their metrics against truth A and against truth B
    series122_rules.csv                the reading-rule numbers and verdicts
    series122_best20_kchan.csv         Kchan quantiles among the best 20% (and among good runs)
    series122_scoring_gate.csv         the self-check comparison (per metric)
    series122_components.csv           volume versus shape: how Kchan relates to beta, r, gamma
    series122_partial_correlations.csv each parameter alone vs with the others held fixed
    series122_model_ladder.csv         held-out R squared of the models
    series122_dial_weights.csv         the fitted dial weights
    series122_shape_beyond_dial.csv    does Kchan move shape metrics beyond the dial?
    series122_null_control.csv         the null-control draws
    fig_series122_1_score_vs_kchan.png   fig_series122_2_best20_cumulative.png
    fig_series122_3_components.png       fig_series122_4_model_ladder.png
  An earlier run's files with the same names are overwritten (they are regenerable).

HOW TO RUN  (from the lab/ folder, like the other scripts)
----------------------------------------------------------
    cd /workspaces/tRIBS-Fall2026/workspaces/SMF/lab
    python analyze_series122_kchan.py --dry_run                      <- reads files and checks, writes nothing
    python -u analyze_series122_kchan.py 2>&1 | tee series122_analysis.log

It takes a few minutes (the cross-validation refits the dial many times).

STATUS  (what has and has not been tested)
------------------------------------------
  TESTED, on a stand-in lab (a toy hydrograph model, 400 fake runs, a verbatim copy of the real scorer's metric code):
    * the whole run, start to finish, in two scenarios (Kchan matters / Kchan inert), figures looked at;
    * the self-check reproduces the scorer to about 1e-16 on all 17 metrics for all 400 runs;
    * eleven bad-input cases (a changed outlet file, two different changed results values, a swapped truth A, two
      truth-A files, a missing column, duplicate run IDs, two event windows, a short results file, no results file,
      a wrong Series 110 record) are each refused with a plain STOPPED message, and the input files were
      byte-identical afterwards;
    * a missing kept outlet file, or a missing null-control file, gives a WARN line and the run carries on without
      that run (or without the null control), saying so in the report;
    * the rule arithmetic was recomputed independently, including the exact boundaries (10% share, 0.95 score,
      the 1.26-decade span, Kchan exactly 30);
    * .qout parsing in five layouts (spaces, commas, comma+space, Windows line ends, no header line);
    * the dial fit recovers a known direction; the same seed gives byte-identical tables.
    The null control was tested with the REAL Series 121 results file (its 100 f0p012_true rows).
  NOT TESTED: your real kept outlet files and your real lhs_results_kchan_122.csv. Two things could differ there:
    * reading the real .qout files with this script's own reader (the same reader the sweep and truth-B check used);
    * whether this script's time conversion matches pytRIBS's exactly.
  If either differs, the self-check (stage 2 of the report) fails and the script prints a STOPPED box and writes no
  results. It does not carry on with numbers it cannot vouch for.
  Every number from the stand-in tests is a test artefact, not a scientific result.

OPTIONS
-------
    --dry_run          check the inputs and say what would be written, then stop
    --results PATH     the sweep's results file        --null_csv PATH    the Series 121 file for the null control
    --out_dir PATH     where to write                   --no_null          skip the null control
    --cv_reps N        cross-validation repeats (default 2)       --cv_folds N   folds (default 5)
    --null_reps N      random Kchan draws for the cheap null statistics (default 500)
    --null_cv_draws N  random draws for the null's held-out gain (default 12)
    --seed N           random seed for splits and null draws (default 0; same seed, same numbers)
    --quick            1 repeat and 3 null draws (faster, noisier)
    --no_plots         tables only          --dpi N   figure resolution (default 150)

READ THIS BEFORE QUOTING ANY NUMBER
-----------------------------------
  * One synthetic truth pair, one storm (12 Aug 2014, 1.0x), Ks 7.0 and f 0.012 fixed at the truth values (Stage 1).
    Nothing here says what happens when Ks and f are also unknown.
  * A good Kchan estimate here does not mean Kchan is identifiable with real data, noise, or a different storm.
  * The 400 runs are one Latin hypercube draw (seed 301). A different draw would move every percentile a little.
  * The best run's Kchan is one number from one draw; the rules use the best 20% (80 runs) for that reason.
  * The null control uses 100 runs, the sweep uses 400, so the null is noisier than the real analysis.
  * The tRIBS mechanism behind any pattern is a hypothesis to check, not something these numbers show.
"""

import argparse
import hashlib
import io
import math
import re
import sys
import time
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize

# ======================================================================
# FILES AND CHECKSUMS
# ======================================================================
SERIES = 122
RESULTS_NAME = "lhs_results_kchan_122.csv"
KEPT_DIR_NAME = "kept_qout_122"
NULL_CSV_NAME = "lhs_results_anchor_ALL_121.csv"
OUT_SUBFOLDER = "series122_kchan"

TRUTH_A_MD5 = "4a433e441359313b6629734510f27813"      # SMF_20140812_60_Ks7p0x_truthcheck_Outlet.qout
TRUTH_B_NAME = "SMF_20140812_110_cc201p0_loc_K7p0_f0p012_Outlet.qout"
TRUTH_B_MD5 = "013097a414def7b59aca731f9665f632"      # verified by verify_truth_B_122_tribs6.py on 2026-10-26

FILES = dict(
    report="series122_report.txt", scored="series122_scored_both_truths.csv", rules="series122_rules.csv",
    best20="series122_best20_kchan.csv", gate="series122_scoring_gate.csv", components="series122_components.csv",
    partial="series122_partial_correlations.csv", ladder="series122_model_ladder.csv",
    dial="series122_dial_weights.csv", shape="series122_shape_beyond_dial.csv", null="series122_null_control.csv",
    fig1="fig_series122_1_score_vs_kchan.png", fig2="fig_series122_2_best20_cumulative.png",
    fig3="fig_series122_3_components.png", fig4="fig_series122_4_model_ladder.png",
)

# ======================================================================
# THE READING RULES: fixed before any result. Do not change these to get a different answer.
# ======================================================================
TRUE_KCHAN = 201.0                 # truth B
TOP_FRAC = 0.20                    # "the best 20%"
MIDDLE = (0.10, 0.90)              # "the middle 80%"
SPAN_MAX_FRACTION_OF_BOX = 0.5     # the middle-80% span must cover less than half of the log box
GOOD_KGE = 0.95                    # "good runs"
FALSE_USE_KCHAN = 30.0             # mm/hr
FALSE_USE_SHARE = 0.10             # false use if more than this share of the best 20% exceed FALSE_USE_KCHAN
DIRECTION_GAIN = 0.02              # held-out R squared gain that counts as "adds a direction"

# ======================================================================
# PARAMETERS AND THE SAMPLING BOX (same as the sweep script)
# ======================================================================
P_CV, P_FE, P_N, P_K = "kinemvelcoef", "flowexp", "channelroughness", "channelconductivity_mmhr"
ROUTING = [P_CV, P_FE, P_N]
FOUR = ROUTING + [P_K]
SHORT = {P_CV: "cv", P_FE: "flowexp", P_N: "n", P_K: "Kchan"}
BOX = {P_CV: (2.5, 6.5), P_FE: (0.18, 0.35), P_N: (0.02, 0.10), P_K: (3.0, 1000.0)}
BOX_TOL = 1e-6
I_CV, I_FE, I_N, I_K = 0, 1, 2, 3
BOX_DECADES = math.log10(BOX[P_K][1] / BOX[P_K][0])          # 2.5229
W121 = np.array([-0.36, 0.91, 0.21])                        # Series 121 dial weights (cv, flowexp, n), standardized
NULL_ANCHOR = "f0p012_true"

# ======================================================================
# SCORING (a copy of the sweep scorer's method: run_sensitivity_single_interp_tribs6.py)
# ======================================================================
ORIGIN = pd.Timestamp("2014-08-01")        # .qout Time column is hours from 2014-08-01 00:00
THRESHOLD_FRAC = 0.05                      # phase metrics use 5% of the observed peak
MIN_QOUT_ROWS = 1000

# metrics scored per truth; the first block is what the sweep CSV also holds (used for the self-check)
GATE_METRICS = ["kge_2012", "pbias_pct", "kge_r", "kge_beta", "kge_gamma", "kge_alpha", "kge", "nse", "rmse_m3s",
                "peak_error_pct", "peak_timing_error_hr", "volume_error_pct", "first_arrival_error_min",
                "rising_limb_steepness_ratio", "time_to_peak_from_exc_min", "duration_above_thresh_error_min",
                "recession_rate_ratio"]
GATE_TOL = 1e-6
# what Series 110 recorded for truth B scored against truth A (truth_location_110.csv, row cc201p0, Ks 7, f 0.012)
REC_B_VS_A = dict(pbias_pct=-9.98684101102667, kge_2012=0.8841407468110921, kge_beta=0.9001315898897334)
REC_TOL = 1e-6

# volume and shape metrics tested one by one for "does Kchan add a direction beyond the dial"
SHAPE_TEST = [
    ("pbias_pct", "PBIAS %", "volume"),
    ("kge_beta", "beta (volume ratio)", "volume"),
    ("kge_r", "r (correlation)", "shape"),
    ("kge_gamma", "gamma (variability ratio)", "shape"),
    ("peak_error_pct", "peak error %", "shape"),
    ("rising_limb_steepness_ratio", "rising-limb steepness", "shape"),
    ("first_arrival_error_min", "first arrival error (min)", "timing"),
    ("recession_rate_ratio", "recession rate ratio", "shape"),
]

# colours (the surface and ink of the other Series 120/121 figures)
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SOFT = "#52514e"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
GREY_DOT = "#b4b3ab"
GOLD = "#c08000"       # the Series 121 gold, darkened a little so it passes the palette contrast check on this surface
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"


class DataProblem(Exception):
    """A problem with the inputs that makes the analysis meaningless: printed as a plain message."""


class GateFailure(DataProblem):
    """The self-check failed: the script's scoring does not reproduce the sweep's own scores."""


# ======================================================================
# Small helpers
# ======================================================================
def md5_of(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def sec(title):
    print("\n" + "=" * 94)
    print(title)
    print("=" * 94)


def fnum(x, w=8, d=3, sign=False):
    if x is None or not np.isfinite(x):
        return "n/a".rjust(w)
    return (f"{x:+.{d}f}" if sign else f"{x:.{d}f}").rjust(w)


def fk(x):
    """Format a Kchan value in mm/hr."""
    if x is None or not np.isfinite(x):
        return "n/a"
    return f"{x:.1f}" if x < 100 else f"{x:.0f}"


class Tee(io.TextIOBase):
    """Send everything printed to the screen and also keep a copy for the report file."""

    def __init__(self, screen):
        self.screen, self.buf = screen, io.StringIO()

    def write(self, s):
        self.screen.write(s)
        self.buf.write(s)
        return len(s)

    def flush(self):
        self.screen.flush()


# ======================================================================
# Reading outlet files and scoring (copy of the sweep scorer's method)
# ======================================================================
def read_qout_arrays(path):
    """Parse a tRIBS outlet file by hand: one header line, then comma and/or whitespace separated rows.
    Returns (time_hr, q_m3s)."""
    t, q = [], []
    first = True
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            toks = [x for x in re.split(r"[,\s]+", line) if x]
            if first:
                first = False
                try:
                    float(toks[0])
                except ValueError:
                    continue                                   # header line
            try:
                tv, qv = float(toks[0]), float(toks[1])
            except (ValueError, IndexError):
                continue
            t.append(tv)
            q.append(qv)
    if len(t) < MIN_QOUT_ROWS:
        raise DataProblem(f"{Path(path).name} has only {len(t)} readable rows (expected thousands).")
    return np.array(t, dtype=float), np.array(q, dtype=float)


def qout_to_grid(path):
    """Read an outlet file and put it on the 5-minute grid exactly as the sweep scorer does."""
    t, q = read_qout_arrays(path)
    raw = pd.DataFrame({"Time_hr": t, "Q": q})
    raw["datetime"] = pd.to_datetime(raw["Time_hr"] * 3600, unit="s", origin=ORIGIN)
    raw = raw.set_index("datetime")
    return raw["Q"].resample("5min").interpolate(method="time")


def phase_metrics(obs, sim):
    """The five hydrograph-phase metrics, copied from the sweep scorer (threshold = 5% of the observed peak)."""
    threshold = THRESHOLD_FRAC * obs.max()
    dt_min = (obs.index[1] - obs.index[0]).total_seconds() / 60.0
    obs_above = obs >= threshold
    obs_cross = obs_above.idxmax() if obs_above.any() else None
    obs_tpeak = obs.idxmax()
    sim_above = sim >= threshold
    sim_cross = sim_above.idxmax() if sim_above.any() else None

    if obs_cross is not None and sim_cross is not None:
        first_arrival = (sim_cross - obs_cross).total_seconds() / 60.0
        obs_rise = max((obs_tpeak - obs_cross).total_seconds() / 60.0, dt_min)
        obs_slope = (obs.max() - threshold) / obs_rise
        sim_tpeak = sim.idxmax()
        sim_rise = max((sim_tpeak - sim_cross).total_seconds() / 60.0, dt_min)
        sim_slope = (sim.max() - threshold) / sim_rise
        steep = sim_slope / obs_slope if obs_slope > 0 else np.nan
        sim_ttp = (sim.idxmax() - sim_cross).total_seconds() / 60.0
        ttp = sim_ttp - (obs_tpeak - obs_cross).total_seconds() / 60.0
    else:
        first_arrival = steep = ttp = np.nan
    dur_err = sim_above.sum() * dt_min - obs_above.sum() * dt_min

    def rec_slope(series, tpeak):
        rec = series.loc[tpeak:]
        rec = rec[rec >= threshold]
        rec = rec[rec > 0]
        if len(rec) < 3:
            return np.nan
        return np.polyfit(np.arange(len(rec), dtype=float), np.log(rec.values), 1)[0]

    o, s = rec_slope(obs, obs_tpeak), rec_slope(sim, sim.idxmax())
    rec_ratio = s / o if (np.isfinite(o) and np.isfinite(s) and o != 0) else np.nan
    return dict(first_arrival_error_min=first_arrival, rising_limb_steepness_ratio=steep,
                time_to_peak_from_exc_min=ttp, duration_above_thresh_error_min=dur_err, recession_rate_ratio=rec_ratio)


def score_pair(obs5, sim5, event_start, event_end):
    """Score a simulated series against an observed one: the sweep scorer's formulas. Returns a dict, or None if the
    two do not overlap in the event window."""
    cmp_ = pd.DataFrame({"Observed": obs5, "Simulated": sim5}).dropna()
    ev = cmp_.loc[event_start:event_end]
    if len(ev) < 10:
        return None
    obs, sim = ev["Observed"], ev["Simulated"]
    dt = (ev.index[1] - ev.index[0]).total_seconds()
    obs_peak, sim_peak = obs.max(), sim.max()
    obs_vol, sim_vol = obs.sum() * dt, sim.sum() * dt
    with np.errstate(all="ignore"):
        r = np.corrcoef(sim, obs)[0, 1]
        alpha = np.std(sim) / np.std(obs)
        beta = np.mean(sim) / np.mean(obs)
        gamma = alpha / beta
        out = dict(
            kge_2012=1 - np.sqrt((r - 1) ** 2 + (gamma - 1) ** 2 + (beta - 1) ** 2),
            kge=1 - np.sqrt((r - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2),
            kge_r=r, kge_alpha=alpha, kge_beta=beta, kge_gamma=gamma,
            pbias_pct=100 * (np.sum(sim - obs) / np.sum(obs)),
            nse=1 - np.sum((sim - obs) ** 2) / np.sum((obs - obs.mean()) ** 2),
            rmse_m3s=np.sqrt(np.mean((sim - obs) ** 2)),
            peak_error_pct=(sim_peak - obs_peak) / obs_peak * 100,
            peak_timing_error_hr=(sim.idxmax() - obs.idxmax()).total_seconds() / 3600,
            volume_error_pct=(sim_vol - obs_vol) / obs_vol * 100,
        )
    out.update(phase_metrics(obs, sim))
    out["n_points"] = len(ev)
    return {k: (float(v) if np.isscalar(v) else v) for k, v in out.items()}


# ======================================================================
# Locating and checking the inputs
# ======================================================================
def locate(args):
    root = Path.cwd().parent
    calib = root / "calibration_work"
    summ = calib / "03_comparisons" / "summary_tables"
    p = dict(
        calib=calib,
        results=args.results or (summ / RESULTS_NAME),
        kept=calib / KEPT_DIR_NAME,
        truth_a_dir=calib / "synth_truth",
        truth_b=calib / "synth_truth_location" / TRUTH_B_NAME,
        null_csv=args.null_csv or (summ / NULL_CSV_NAME),
        out=args.out_dir or (calib / "03_comparisons" / "sensitivity_plots" / OUT_SUBFOLDER),
    )
    return p


def load_sweep(path):
    path = Path(path)
    if not path.exists():
        raise DataProblem(f"Could not find the results file {path}\n"
                          "Run this from the lab/ folder, or point to the file with --results.")
    try:
        df = pd.read_csv(path)
    except Exception as exc:                                            # noqa: BLE001
        raise DataProblem(f"Could not read {path} as a CSV: {exc}")
    need = ["run_id"] + FOUR + ["kge_2012", "pbias_pct", "event_start", "event_end", "qout_md5"]
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise DataProblem("The results file is missing these column(s): " + ", ".join(missing)
                          + ". Is this lhs_results_kchan_122.csv?")
    return df


def preflight(paths, args):
    """Read-only checks. Returns a dict of everything the later steps need. Raises DataProblem on a real problem."""
    sec("0. Inputs and checks (nothing is written yet)")
    notes = []
    df = load_sweep(paths["results"])
    print(f"  results file    {Path(paths['results']).name}: {len(df)} rows, {df.shape[1]} columns")

    for c in FOUR + ["kge_2012", "pbias_pct"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    bad = df[FOUR + ["kge_2012"]].isna().any(axis=1)
    if bad.any():
        print(f"  NOTE  {int(bad.sum())} row(s) have a missing parameter or KGE_2012 and are left out.")
        df = df[~bad].copy()
    if df["run_id"].duplicated().any():
        raise DataProblem("A run_id appears more than once in the results file; counting a run twice would distort "
                          "every number.")
    print("  ok    every run_id is unique")
    if "series" in df.columns and not (df["series"] == SERIES).all():
        print(f"  WARN  the 'series' column is not {SERIES} everywhere: {sorted(df['series'].unique().tolist())}")
    if "optpercolation" in df.columns and not (df["optpercolation"] == 1).all():
        print("  WARN  optpercolation is not 1 in every run (channel loss should be ON in every Series 122 run)")
    for p in FOUR:
        lo, hi = BOX[p]
        out = ((df[p] < lo * (1 - BOX_TOL) - BOX_TOL) | (df[p] > hi * (1 + BOX_TOL) + BOX_TOL)).sum()
        if out:
            print(f"  WARN  {int(out)} run(s) have {SHORT[p]} outside the sampling box {lo} to {hi}")
    print("  ok    parameters checked against the sampling box")

    starts, ends = df["event_start"].astype(str).unique(), df["event_end"].astype(str).unique()
    if len(starts) != 1 or len(ends) != 1:
        raise DataProblem(f"The results file has more than one event window ({list(starts)} to {list(ends)}). "
                          "This script assumes the sweep used a single window.")
    ev_start, ev_end = starts[0], ends[0]
    print(f"  ok    one event window: {ev_start} to {ev_end}")

    # truth A
    ta_dir = paths["truth_a_dir"]
    ta = sorted(ta_dir.glob("*.qout")) if ta_dir.exists() else []
    if len(ta) != 1:
        raise DataProblem(f"Expected exactly one *.qout file in {ta_dir} (truth A), found {len(ta)}.")
    ta_path = ta[0]
    ta_md5 = md5_of(ta_path)
    if ta_md5 != TRUTH_A_MD5:
        raise DataProblem(f"Truth A ({ta_path.name}) has checksum {ta_md5}, not the {TRUTH_A_MD5} that Series 122 was "
                          "scored against. Do not continue with a different file.")
    if "truth_md5" in df.columns and not (df["truth_md5"].astype(str) == TRUTH_A_MD5).all():
        raise DataProblem("The results file says some runs were scored against a truth with a different checksum.")
    print(f"  ok    truth A  {ta_path.name}  [md5 {ta_md5[:12]}]  matches the file Series 122 used")

    # truth B
    tb_path = Path(paths["truth_b"])
    if not tb_path.exists():
        raise DataProblem(f"Truth B not found: {tb_path}")
    tb_md5 = md5_of(tb_path)
    if tb_md5 != TRUTH_B_MD5:
        raise DataProblem(f"Truth B has checksum {tb_md5}, not the verified {TRUTH_B_MD5}. Do not continue.")
    print(f"  ok    truth B  {tb_path.name}  [md5 {tb_md5[:12]}]  is the file verify_truth_B_122_tribs6.py verified")

    # kept outlet files
    kept = Path(paths["kept"])
    if not kept.is_dir():
        raise DataProblem(f"The kept outlet folder {kept} does not exist.")
    paths_by_run, missing, md5_bad = {}, [], []
    for r in df.itertuples():
        f = kept / f"{r.run_id}_Outlet.qout"
        if not f.exists():
            missing.append(r.run_id)
            continue
        paths_by_run[r.run_id] = f
        want = str(r.qout_md5)
        if want and want != "nan" and md5_of(f) != want:
            md5_bad.append(r.run_id)
    if md5_bad:
        raise DataProblem(f"{len(md5_bad)} kept outlet file(s) do not match the checksum the sweep recorded (for "
                          f"example {md5_bad[0]}). The kept files must be the ones the sweep scored.")
    if missing:
        print(f"  WARN  {len(missing)} run(s) have no kept outlet file and are left out (for example {missing[0]})")
        df = df[~df["run_id"].isin(missing)].copy()
    if len(df) < 100:
        raise DataProblem(f"Only {len(df)} usable runs are left; the analysis needs at least 100.")
    print(f"  ok    {len(paths_by_run)} kept outlet files found; each checksum equals the one the sweep recorded")

    # null-control file
    null_ok = False
    if args.no_null:
        print("  NOTE  --no_null: the null control is skipped")
    elif Path(paths["null_csv"]).exists():
        null_ok = True
        print(f"  ok    null-control file  {Path(paths['null_csv']).name}")
    else:
        print(f"  WARN  null-control file {paths['null_csv']} not found; the null control will be skipped "
              "(point to it with --null_csv)")

    out_dir = Path(paths["out"])
    if out_dir.resolve() == Path(paths["results"]).resolve().parent:
        raise DataProblem(f"The output folder {out_dir} is the folder holding the results file; refusing to write there.")
    return dict(df=df.reset_index(drop=True), ev_start=ev_start, ev_end=ev_end, ta_path=ta_path, tb_path=tb_path,
                kept_paths=paths_by_run, null_ok=null_ok, out_dir=out_dir)


# ======================================================================
# Scoring all runs against both truths, and the self-check
# ======================================================================
def score_all(pre):
    sec("1. Scoring the 400 kept outlet files against truth A and truth B")
    df = pre["df"]
    obs_a = qout_to_grid(pre["ta_path"])
    obs_b = qout_to_grid(pre["tb_path"])
    rows = []
    t0 = time.time()
    for i, r in enumerate(df.itertuples(), 1):
        sim = qout_to_grid(pre["kept_paths"][r.run_id])
        sa = score_pair(obs_a, sim, pre["ev_start"], pre["ev_end"])
        sb = score_pair(obs_b, sim, pre["ev_start"], pre["ev_end"])
        if sa is None or sb is None:
            raise DataProblem(f"Run {r.run_id}: its outlet series does not overlap the event window "
                              f"({pre['ev_start']} to {pre['ev_end']}).")
        row = {"run_id": r.run_id}
        row.update({f"{k}_A": v for k, v in sa.items()})
        row.update({f"{k}_B": v for k, v in sb.items()})
        rows.append(row)
        if i % 100 == 0 or i == len(df):
            print(f"  scored {i} of {len(df)}   ({time.time() - t0:.0f} s)")
    scored = pd.DataFrame(rows)
    npts = scored[["n_points_A", "n_points_B"]].drop_duplicates()
    if len(npts) != 1:
        print(f"  WARN  the number of points in the event window differs between runs: {npts.to_dict('records')}")
    # truth B scored against truth A, for the second check
    b_vs_a = score_pair(obs_a, qout_to_grid(pre["tb_path"]), pre["ev_start"], pre["ev_end"])
    return scored, b_vs_a


def gate(pre, scored, b_vs_a):
    """Stop unless the script reproduces (1) the sweep's own scores against truth A and (2) the recorded numbers
    for truth B against truth A. Returns the per-metric comparison table."""
    sec("2. Self-check: does this script's scoring reproduce the numbers already on record?")
    df = pre["df"].merge(scored, on="run_id", how="left")
    rows, fails = [], []
    print("  (a) all runs vs truth A, against the sweep's own columns in the results file")
    print(f"      {'metric':<34}{'runs compared':>14}{'max |difference|':>20}   result")
    for m in GATE_METRICS:
        if m not in pre["df"].columns:
            continue
        mine, theirs = df[f"{m}_A"].to_numpy(float), pd.to_numeric(df[m], errors="coerce").to_numpy(float)
        both_nan = np.isnan(mine) & np.isnan(theirs)
        one_nan = np.isnan(mine) ^ np.isnan(theirs)
        ok_mask = ~np.isnan(mine) & ~np.isnan(theirs)
        diff = np.abs(mine - theirs)[ok_mask]
        tol = GATE_TOL + GATE_TOL * np.abs(theirs[ok_mask])
        n_bad = int((diff > tol).sum()) + int(one_nan.sum())
        mx = float(diff.max()) if len(diff) else 0.0
        status = "ok" if n_bad == 0 else f"FAIL ({n_bad} runs)"
        print(f"      {m:<34}{int(ok_mask.sum()):>14}{mx:>20.3g}   {status}")
        rows.append(dict(check="all runs vs truth A", metric=m, n_compared=int(ok_mask.sum()), max_abs_diff=mx,
                         n_fail=n_bad, passed=n_bad == 0))
        if n_bad:
            fails.append(f"{m}: {n_bad} run(s) differ by more than {GATE_TOL:g} (largest {mx:.3g})")
    print("  (b) truth B scored against truth A, against what Series 110 recorded")
    for m, rec in REC_B_VS_A.items():
        got = b_vs_a[m]
        ok = abs(got - rec) <= REC_TOL * (1 + abs(rec))
        print(f"      {m:<34} this script {got:>14.8f}   recorded {rec:>14.8f}   {'ok' if ok else 'FAIL'}")
        rows.append(dict(check="truth B vs truth A", metric=m, n_compared=1, max_abs_diff=abs(got - rec),
                         n_fail=0 if ok else 1, passed=ok))
        if not ok:
            fails.append(f"truth B vs A {m}: this script {got:.8f}, Series 110 recorded {rec:.8f}")
    tab = pd.DataFrame(rows)
    if fails:
        raise GateFailure("this script's scoring does NOT reproduce numbers that are already on record, so a truth-B "
                          "score from it cannot be trusted yet:\n    - " + "\n    - ".join(fails)
                          + "\n  Nothing else was computed. Paste this section to Claude; do not use any other output.")
    print("  PASS  the scoring reproduces the sweep's own truth-A scores and Series 110's truth-B record.")
    return tab


# ======================================================================
# The reading rules
# ======================================================================
def top_set(d, col, frac=TOP_FRAC):
    k = max(1, int(math.ceil(frac * len(d))))
    return d.nlargest(k, col)


def rule_recovery(d, col):
    """Rule 1 and 2 of the recovery rule, against one truth. d needs columns Kchan, logK and col."""
    top = top_set(d, col)
    lo_l, hi_l = np.quantile(top["logK"].to_numpy(), MIDDLE)
    span = float(hi_l - lo_l)
    inside = bool(lo_l <= math.log10(TRUE_KCHAN) <= hi_l)
    narrow = span < SPAN_MAX_FRACTION_OF_BOX * BOX_DECADES
    best = d.loc[d[col].idxmax()]
    return dict(n_best=len(top), kchan_p10=10 ** lo_l, kchan_p90=10 ** hi_l, span_decades=span,
                span_limit_decades=SPAN_MAX_FRACTION_OF_BOX * BOX_DECADES, true_inside=inside, span_narrow=narrow,
                recovered=bool(inside and narrow), best_run_kchan=float(best["Kchan"]), best_run_kge=float(best[col]),
                top_median_kchan=float(10 ** np.median(top["logK"])), top_min_kchan=float(top["Kchan"].min()),
                top_max_kchan=float(top["Kchan"].max()))


def rule_false_use(d, col):
    top = top_set(d, col)
    good = d[d[col] >= GOOD_KGE]
    share = float((top["Kchan"] > FALSE_USE_KCHAN).mean())
    lowest_fifth = BOX[P_K][0] * (BOX[P_K][1] / BOX[P_K][0]) ** 0.2
    lo_l, hi_l = np.quantile(top["logK"].to_numpy(), MIDDLE)
    best = d.loc[d[col].idxmax()]
    return dict(n_best=len(top), n_good=len(good), largest_kchan_good=float(good["Kchan"].max()) if len(good) else np.nan,
                share_best_above_30=share, false_use=bool(share > FALSE_USE_SHARE),
                share_best_in_lowest_fifth=float((top["Kchan"] < lowest_fifth).mean()), lowest_fifth_limit=lowest_fifth,
                best_run_kchan=float(best["Kchan"]), best_run_kge=float(best[col]),
                kchan_p10=10 ** lo_l, kchan_p90=10 ** hi_l, span_decades=float(hi_l - lo_l),
                spearman_logK_vs_kge=float(stats.spearmanr(d["logK"], d[col])[0]))


# ======================================================================
# The dial and the cross-validated models
# ======================================================================
def poly1(t, deg=3):
    return np.column_stack([t ** k for k in range(deg + 1)])


def poly_nc(u, deg=3):
    """Polynomial in u without the constant column."""
    return np.column_stack([u ** k for k in range(1, deg + 1)])


def polyfull(Z, deg=3):
    """Every polynomial term up to degree deg in all columns of Z, with every cross term (35 columns for 4 variables)."""
    import itertools
    cols = [np.ones(len(Z))]
    for dgr in range(1, deg + 1):
        for idx in itertools.combinations_with_replacement(range(Z.shape[1]), dgr):
            cols.append(np.prod(Z[:, list(idx)], axis=1))
    return np.column_stack(cols)


def _unit(v):
    n = float(np.linalg.norm(v))
    if n == 0 or not np.isfinite(n):
        out = np.zeros(len(v))
        out[0] = 1.0
        return out
    return v / n


def fit_index(Z, y, extra=None, n_start=None, polish=3, seed=2026):
    """Best direction w for  y = g(Z @ w) [+ extra terms], g a cubic. A coarse search over many directions (a fixed,
    seeded set, so the same data always gives the same answer), then a local polish from the three best.
    The sign is fixed by convention: the flowexp weight is positive. Returns (w, coefficients, SSE)."""
    d = Z.shape[1]
    rng = np.random.default_rng(seed)
    if n_start is None:
        n_start = 300 if d == 3 else 900
    cand = rng.normal(size=(n_start, d))
    cand /= np.linalg.norm(cand, axis=1, keepdims=True)
    cand = np.vstack([np.eye(d), -np.eye(d), cand])

    def design(w):
        X = poly1(Z @ w)
        return X if extra is None else np.column_stack([X, extra])

    def sse_w(w):
        X = design(w)
        b = np.linalg.lstsq(X, y, rcond=None)[0]
        r = y - X @ b
        return float(r @ r)

    vals = np.array([sse_w(c) for c in cand])
    order = np.argsort(vals)[:polish]
    best_w, best_v = cand[order[0]], float(vals[order[0]])
    for i in order:
        res = minimize(lambda v: sse_w(_unit(v)), cand[i], method="Nelder-Mead",
                       options=dict(xatol=1e-3, fatol=1e-10, maxiter=250 * d))
        if res.fun < best_v:
            best_w, best_v = _unit(res.x), float(res.fun)
    w = best_w if best_w[I_FE] >= 0 else -best_w
    X = design(w)
    b = np.linalg.lstsq(X, y, rcond=None)[0]
    r = y - X @ b
    return w, b, float(r @ r)


MODEL_LABELS = {
    "dial": "the dial alone (cv, flowexp, n)",
    "curve": "log Kchan alone (cubic)",
    "dial_curve": "the dial + a curve in log Kchan",
    "dial4": "ONE dial of all four (Kchan inside it)",
    "full4": "full cubic in all four (35 terms)",
}


def cv_models(Z, u, y, which, reps, folds, seed):
    """Repeated k-fold cross-validated R squared (held-out runs only), for the models in `which`, all on the SAME
    splits so differences between models are paired. Z = standardized cv, flowexp, n; u = standardized log10 Kchan."""
    n = len(y)
    rng = np.random.default_rng(seed)
    sst = float(((y - y.mean()) ** 2).sum())
    per_rep = {m: [] for m in which}
    for _ in range(reps):
        perm = rng.permutation(n)
        sse = {m: 0.0 for m in which}
        for f in np.array_split(perm, folds):
            tr = np.setdiff1d(perm, f)
            Ztr, utr, ytr, Zte, ute, yte = Z[tr], u[tr], y[tr], Z[f], u[f], y[f]
            if "dial" in which:
                w, b, _ = fit_index(Ztr, ytr)
                sse["dial"] += float(((yte - poly1(Zte @ w) @ b) ** 2).sum())
            if "dial_curve" in which:
                w, b, _ = fit_index(Ztr, ytr, extra=poly_nc(utr))
                pred = np.column_stack([poly1(Zte @ w), poly_nc(ute)]) @ b
                sse["dial_curve"] += float(((yte - pred) ** 2).sum())
            if "curve" in which:
                b = np.linalg.lstsq(poly1(utr), ytr, rcond=None)[0]
                sse["curve"] += float(((yte - poly1(ute) @ b) ** 2).sum())
            if "dial4" in which:
                Z4tr, Z4te = np.column_stack([Ztr, utr]), np.column_stack([Zte, ute])
                w, b, _ = fit_index(Z4tr, ytr)
                sse["dial4"] += float(((yte - poly1(Z4te @ w) @ b) ** 2).sum())
            if "full4" in which:
                Z4tr, Z4te = np.column_stack([Ztr, utr]), np.column_stack([Zte, ute])
                b = np.linalg.lstsq(polyfull(Z4tr), ytr, rcond=None)[0]
                sse["full4"] += float(((yte - polyfull(Z4te) @ b) ** 2).sum())
        for m in which:
            per_rep[m].append(1.0 - sse[m] / sst)
    return {m: float(np.mean(v)) for m, v in per_rep.items()}, per_rep


def standardize_cols(X):
    mu, sd = X.mean(axis=0), X.std(axis=0, ddof=1)
    sd = np.where(sd == 0, 1.0, sd)
    return (X - mu) / sd, mu, sd


def prepare_matrices(d):
    """Standardized routing parameters Z (n x 3) and standardized log10 Kchan u (n,)."""
    Z, _, _ = standardize_cols(d[ROUTING].to_numpy(float))
    u, _, _ = standardize_cols(d[["logK"]].to_numpy(float))
    return Z, u[:, 0]


def rule_direction(d, col, args, seed_offset=0):
    """KCHAN ADDS A DIRECTION: held-out R squared of (dial + curve in log Kchan) minus that of the dial alone,
    on the same cross-validation splits."""
    ok = d[col].notna().to_numpy()
    Z, u = prepare_matrices(d[ok])
    y = d.loc[ok, col].to_numpy(float)
    res, per = cv_models(Z, u, y, ["dial", "dial_curve"], args.cv_reps, args.cv_folds, args.seed + seed_offset)
    gains = np.array(per["dial_curve"]) - np.array(per["dial"])
    gain = float(res["dial_curve"] - res["dial"])
    return dict(r2_dial=res["dial"], r2_dial_curve=res["dial_curve"], gain=gain, gain_min=float(gains.min()),
                gain_max=float(gains.max()), adds_direction=bool(gain > DIRECTION_GAIN), n=int(ok.sum()))


# ======================================================================
# Information tables
# ======================================================================
def best20_table(d, col, label):
    top = top_set(d, col)
    good = d[d[col] >= GOOD_KGE]
    rows = []
    qs = [0, 0.10, 0.25, 0.50, 0.75, 0.90, 1.0]
    for set_name, sub in (("best 20%", top), (f"good runs (KGE_2012 >= {GOOD_KGE:g})", good)):
        if len(sub) == 0:
            continue
        row = dict(truth=label, set=set_name, n_runs=len(sub), kge_min=float(sub[col].min()), kge_max=float(sub[col].max()))
        for q in qs:
            row[f"kchan_q{int(round(q * 100)):03d}"] = float(10 ** np.quantile(sub["logK"].to_numpy(), q))
        for p in ROUTING:
            row[f"{SHORT[p]}_min"], row[f"{SHORT[p]}_max"] = float(sub[p].min()), float(sub[p].max())
        rows.append(row)
    return rows


def component_table(d, label, suffix):
    """Volume versus shape: how Kchan relates to beta (volume), r and gamma (shape), and which of the three limits the
    score. KGE_2012 = 1 - sqrt((r-1)^2 + (gamma-1)^2 + (beta-1)^2), so the three squared gaps add up."""
    rows = []
    Z, _, _ = standardize_cols(d[ROUTING].to_numpy(float))
    u = d["logK"].to_numpy(float)
    for key, name in (("kge_beta", "beta (volume ratio; 1 = right volume)"), ("pbias_pct", "PBIAS %"),
                      ("kge_r", "r (correlation; 1 = right timing/shape)"),
                      ("kge_gamma", "gamma (variability ratio; 1 = right shape)"), ("kge_2012", "KGE_2012")):
        y = d[f"{key}_{suffix}"].to_numpy(float)
        ok = np.isfinite(y)
        rho = float(stats.spearmanr(u[ok], y[ok])[0])
        # partial correlation of log Kchan with y, holding the three routing parameters (linear) fixed
        X = np.column_stack([np.ones(ok.sum()), Z[ok]])
        ru = u[ok] - X @ np.linalg.lstsq(X, u[ok], rcond=None)[0]
        ry = y[ok] - X @ np.linalg.lstsq(X, y[ok], rcond=None)[0]
        pr = float(np.corrcoef(ru, ry)[0, 1]) if ru.std() > 0 and ry.std() > 0 else np.nan
        rows.append(dict(truth=label, quantity=name, spearman_with_logK=rho, partial_r_with_logK_routing_fixed=pr,
                         value_median=float(np.median(y[ok])), value_min=float(y[ok].min()), value_max=float(y[ok].max())))
    # which component limits the score
    gaps = np.column_stack([(d[f"kge_beta_{suffix}"] - 1) ** 2, (d[f"kge_r_{suffix}"] - 1) ** 2,
                            (d[f"kge_gamma_{suffix}"] - 1) ** 2]).astype(float)
    tot = gaps.sum(axis=1)
    ok = tot > 0
    share = gaps[ok] / tot[ok][:, None]
    top_idx = np.argsort(-d[f"kge_2012_{suffix}"].to_numpy(float))[: int(math.ceil(TOP_FRAC * len(d)))]
    top_share = gaps[top_idx] / np.where(tot[top_idx] > 0, tot[top_idx], 1)[:, None]
    shares = dict(all_runs=np.median(share, axis=0), best_20=np.median(top_share, axis=0))
    return rows, shares


def partial_tables(d, suffix, label):
    """Each parameter alone (Pearson, Spearman) and with the other three held fixed. Linear model: signed partial r with
    a 95% interval. Quadratic model (each parameter and its square): partial R squared and unique R squared."""
    y = d[f"kge_2012_{suffix}"].to_numpy(float)
    ok = np.isfinite(y)
    y = y[ok]
    X = {p: (d["logK"] if p == P_K else d[p]).to_numpy(float)[ok] for p in FOUR}
    X = {p: (v - v.mean()) / v.std(ddof=1) for p, v in X.items()}
    n = len(y)

    def sse(cols):
        A = np.column_stack([np.ones(n)] + list(cols))
        b = np.linalg.lstsq(A, y, rcond=None)[0]
        r = y - A @ b
        return float(r @ r), b

    sst = float(((y - y.mean()) ** 2).sum())
    rows = []
    q_all = [c for p in FOUR for c in (X[p], X[p] ** 2)]
    sse_qfull = sse(q_all)[0]
    for p in FOUR:
        others = [X[o] for o in FOUR if o != p]
        s_r, _ = sse(others)
        s_f, b = sse(others + [X[p]])
        pr2 = max(0.0, (s_r - s_f) / s_r) if s_r > 0 else 0.0
        sign = 1.0 if b[-1] >= 0 else -1.0
        pr = sign * math.sqrt(pr2)
        dof = n - len(FOUR) - 1
        if abs(pr) < 1 and dof > 3:
            z, se = math.atanh(pr), 1.0 / math.sqrt(n - (len(FOUR) - 1) - 3)
            ci = (math.tanh(z - 1.96 * se), math.tanh(z + 1.96 * se))
        else:
            ci = (np.nan, np.nan)
        F = ((s_r - s_f) / 1) / (s_f / dof) if s_f > 0 and dof > 0 else np.nan
        pv = float(stats.f.sf(F, 1, dof)) if np.isfinite(F) else np.nan
        # quadratic block
        oq = [c for o in FOUR if o != p for c in (X[o], X[o] ** 2)]
        sq_r, _ = sse(oq)
        dofq = n - 2 * len(FOUR) - 1
        Fq = ((sq_r - sse_qfull) / 2) / (sse_qfull / dofq) if sse_qfull > 0 and dofq > 0 else np.nan
        pq = float(stats.f.sf(Fq, 2, dofq)) if np.isfinite(Fq) else np.nan
        rows.append(dict(truth=label, parameter=SHORT[p], n=n,
                         pearson_alone=float(np.corrcoef(X[p], y)[0, 1]),
                         spearman_alone=float(stats.spearmanr(X[p], y)[0]),
                         partial_r_linear=pr, partial_ci_lo=ci[0], partial_ci_hi=ci[1], partial_p_linear=pv,
                         r2_alone=float(np.corrcoef(X[p], y)[0, 1] ** 2),
                         partial_r2_with_square=max(0.0, (sq_r - sse_qfull) / sq_r) if sq_r > 0 else 0.0,
                         unique_r2_with_square=max(0.0, (sq_r - sse_qfull) / sst), p_block_square=pq))
    return rows


def ladder_and_dials(d, suffix, label, args):
    """Held-out R squared of the model ladder for KGE_2012, plus the fitted dial weights on the full sample."""
    col = f"kge_2012_{suffix}"
    ok = d[col].notna().to_numpy()
    sub = d[ok]
    Z, u = prepare_matrices(sub)
    y = sub[col].to_numpy(float)
    which = ["dial", "curve", "dial_curve", "dial4", "full4"]
    res, per = cv_models(Z, u, y, which, args.cv_reps, args.cv_folds, args.seed)
    ladder = [dict(truth=label, target="KGE_2012", model=m, model_label=MODEL_LABELS[m], held_out_r2=res[m],
                   held_out_r2_min_rep=float(min(per[m])), held_out_r2_max_rep=float(max(per[m])), n=len(y))
              for m in which]
    dials = []
    w3, _, sse3 = fit_index(Z, y)
    w3c, _, sse3c = fit_index(Z, y, extra=poly_nc(u))
    w4, _, sse4 = fit_index(np.column_stack([Z, u]), y)
    cosine = lambda a: float(abs(a @ W121) / (np.linalg.norm(a) * np.linalg.norm(W121)))
    dials.append(dict(truth=label, fit="dial alone", w_cv=w3[0], w_flowexp=w3[1], w_n=w3[2], w_logKchan=np.nan,
                      abs_cosine_with_series121_dial=cosine(w3)))
    dials.append(dict(truth=label, fit="dial inside (dial + log Kchan curve)", w_cv=w3c[0], w_flowexp=w3c[1], w_n=w3c[2],
                      w_logKchan=np.nan, abs_cosine_with_series121_dial=cosine(w3c)))
    dials.append(dict(truth=label, fit="one dial of all four", w_cv=w4[0], w_flowexp=w4[1], w_n=w4[2], w_logKchan=w4[3],
                      abs_cosine_with_series121_dial=cosine(w4[:3]) if np.linalg.norm(w4[:3]) > 0 else np.nan))
    return ladder, dials, res


def shape_beyond_dial(d, suffix, label, args):
    """For each metric: held-out R squared of the dial alone and of dial + curve in log Kchan, and the gain."""
    rows = []
    Z_all, u_all = prepare_matrices(d)
    for key, name, kind in SHAPE_TEST:
        col = f"{key}_{suffix}"
        y_all = d[col].to_numpy(float)
        ok = np.isfinite(y_all)
        if ok.sum() < 100 or np.ptp(y_all[ok]) == 0:
            print(f"  NOTE  {label}: {name}: too few finite values or no spread; skipped.")
            continue
        res, per = cv_models(Z_all[ok], u_all[ok], y_all[ok], ["dial", "dial_curve"], max(1, args.cv_reps - 1),
                             args.cv_folds, args.seed + 7)
        gain = res["dial_curve"] - res["dial"]
        rows.append(dict(truth=label, metric=key, metric_label=name, kind=kind, n=int(ok.sum()),
                         r2_dial=res["dial"], r2_dial_curve=res["dial_curve"], gain=gain,
                         gain_exceeds_yardstick=bool(gain > DIRECTION_GAIN)))
    return rows


# ======================================================================
# Null control
# ======================================================================
def null_control(path, args):
    """Attach random log-uniform Kchan to the 100 f0p012_true runs of the Series 121 file (channel loss was OFF there, so
    Kchan did nothing) and apply the Truth-A rules."""
    try:
        n121 = pd.read_csv(path)
    except Exception as exc:                                       # noqa: BLE001
        print(f"  NOTE  could not read {path} ({exc}); null control skipped.")
        return None
    need = ["anchor_label"] + ROUTING + ["kge_2012"]
    miss = [c for c in need if c not in n121.columns]
    if miss:
        print(f"  NOTE  {Path(path).name} lacks column(s) {miss}; null control skipped.")
        return None
    sub = n121[n121["anchor_label"] == NULL_ANCHOR].copy()
    for c in ROUTING + ["kge_2012"]:
        sub[c] = pd.to_numeric(sub[c], errors="coerce")
    sub = sub.dropna(subset=ROUTING + ["kge_2012"]).reset_index(drop=True)
    if len(sub) < 40:
        print(f"  NOTE  only {len(sub)} usable '{NULL_ANCHOR}' rows in {Path(path).name}; null control skipped.")
        return None
    if "optpercolation" in sub.columns and not (sub["optpercolation"] == 0).all():
        print("  WARN  optpercolation is not 0 in every null-control run; Kchan may not be inert there.")
    n = len(sub)
    rng = np.random.default_rng(args.seed)
    lo_l, hi_l = math.log10(BOX[P_K][0]), math.log10(BOX[P_K][1])
    rows = []
    for k in range(args.null_reps):
        logk = rng.uniform(lo_l, hi_l, n)
        dd = pd.DataFrame({"Kchan": 10 ** logk, "logK": logk, "score": sub["kge_2012"].to_numpy()})
        fu = rule_false_use(dd, "score")
        rows.append(dict(draw=k, largest_kchan_good=fu["largest_kchan_good"], share_best_above_30=fu["share_best_above_30"],
                         false_use=fu["false_use"], span_decades=fu["span_decades"],
                         share_best_in_lowest_fifth=fu["share_best_in_lowest_fifth"], gain=np.nan))
    nd = pd.DataFrame(rows)
    # held-out gain for the first few draws (the cross-validation is the slow part)
    Zs, _, _ = standardize_cols(sub[ROUTING].to_numpy(float))
    rng2 = np.random.default_rng(args.seed)
    for k in range(min(args.null_cv_draws, args.null_reps)):
        logk = rng2.uniform(lo_l, hi_l, n)
        u = (logk - logk.mean()) / logk.std(ddof=1)
        res, per = cv_models(Zs, u, sub["kge_2012"].to_numpy(float), ["dial", "dial_curve"], args.cv_reps,
                             args.cv_folds, args.seed + 100 + k)
        nd.loc[k, "gain"] = res["dial_curve"] - res["dial"]
        nd.loc[k, "r2_dial"], nd.loc[k, "r2_dial_curve"] = res["dial"], res["dial_curve"]
        print(f"    null control: held-out gain, draw {k + 1} of {min(args.null_cv_draws, args.null_reps)}: "
              f"{res['dial_curve'] - res['dial']:+.3f}")
    return dict(draws=nd, n_runs=n, kge=sub["kge_2012"].to_numpy(float), seed=args.seed, file=Path(path).name)


def summarize_null(nc):
    nd = nc["draws"]
    q = lambda s, p: float(np.nanquantile(nd[s].to_numpy(float), p))
    gains = nd["gain"].dropna()
    return dict(
        n_runs=nc["n_runs"], n_draws=len(nd),
        share_above30_median=q("share_best_above_30", 0.5), share_above30_lo=q("share_best_above_30", 0.025),
        share_above30_hi=q("share_best_above_30", 0.975), false_use_fraction_of_draws=float(nd["false_use"].mean()),
        span_median=q("span_decades", 0.5), span_lo=q("span_decades", 0.025), span_hi=q("span_decades", 0.975),
        largest_good_median=q("largest_kchan_good", 0.5), lowest_fifth_median=q("share_best_in_lowest_fifth", 0.5),
        n_gain=len(gains), gain_median=float(gains.median()) if len(gains) else np.nan,
        gain_max=float(gains.max()) if len(gains) else np.nan,
        gain_exceeds_n=int((gains > DIRECTION_GAIN).sum()) if len(gains) else 0)


# ======================================================================
# Printing the answer
# ======================================================================
def verdict_text(flag, yes, no):
    return yes if flag else no


def print_rules(rb, rfa, dir_a, dir_b, nsum, d):
    sec("3. THE READING RULES (fixed before any result was seen)")
    print(f"  Runs: {len(d)}.  'Best 20%' = the {rb['n_best']} runs with the highest KGE_2012 against the truth in question.")
    print(f"  Box: Kchan {BOX[P_K][0]:g} to {BOX[P_K][1]:g} mm/hr = {BOX_DECADES:.2f} decades on the log scale.")
    print()
    print("  QUESTION A - RECOVERY (scored against truth B: channel loss ON, true Kchan = 201 mm/hr)")
    print(f"    best 20% of runs:  Kchan from {fk(rb['top_min_kchan'])} to {fk(rb['top_max_kchan'])}, median {fk(rb['top_median_kchan'])} mm/hr")
    print(f"    middle 80% of the best 20%:  {fk(rb['kchan_p10'])} to {fk(rb['kchan_p90'])} mm/hr   "
          f"(span {rb['span_decades']:.2f} decades; the rule needs < {rb['span_limit_decades']:.2f})")
    print(f"    (1) is 201 inside that range?   {verdict_text(rb['true_inside'], 'YES', 'NO')}")
    print(f"    (2) is the range narrow enough? {verdict_text(rb['span_narrow'], 'YES', 'NO')}")
    print(f"    best single run: Kchan {fk(rb['best_run_kchan'])} mm/hr (KGE_2012 {rb['best_run_kge']:.4f})")
    print(f"    ==> RULE RESULT: {verdict_text(rb['recovered'], 'RECOVERED', 'NOT RECOVERED')}"
          + ("" if rb["recovered"] else
             "   (" + ("201 is outside the range" if not rb["true_inside"] else "the range is too wide") + ")"))
    print()
    print("  QUESTION B - FALSE USE (scored against truth A: NO channel loss, so the right Kchan is 'none')")
    print(f"    runs that score well (KGE_2012 >= {GOOD_KGE:g}): {rfa['n_good']}")
    print(f"    largest Kchan among them:                     {fk(rfa['largest_kchan_good'])} mm/hr")
    print(f"    share of the best 20% with Kchan above {FALSE_USE_KCHAN:g}:  {rfa['share_best_above_30']:.1%}   "
          f"(false use if more than {FALSE_USE_SHARE:.0%})")
    print(f"    ==> RULE RESULT: {verdict_text(rfa['false_use'], 'FALSE USE', 'NO FALSE USE')}")
    print("    pile-up at the lowest Kchan (information, not a rule):")
    print(f"      best 20% with Kchan below {fk(rfa['lowest_fifth_limit'])} mm/hr (the lowest fifth of the log box): "
          f"{rfa['share_best_in_lowest_fifth']:.1%}   (20% if Kchan made no difference)")
    print(f"      best single run: Kchan {fk(rfa['best_run_kchan'])} mm/hr (KGE_2012 {rfa['best_run_kge']:.4f});  "
          f"Spearman(log Kchan, KGE_2012) = {rfa['spearman_logK_vs_kge']:+.2f}")
    print()
    print(f"  KCHAN ADDS A DIRECTION (held-out R squared of dial + log-Kchan curve minus the dial alone; the rule: > {DIRECTION_GAIN:g})")
    print("    The handoff does not say which truth this is judged against, so both are shown.")
    for lab, dr in (("scored against truth B", dir_b), ("scored against truth A", dir_a)):
        print(f"    {lab:<24} dial {dr['r2_dial']:.3f}   dial + curve {dr['r2_dial_curve']:.3f}   "
              f"gain {dr['gain']:+.3f}  (repeats {dr['gain_min']:+.3f} to {dr['gain_max']:+.3f})   "
              f"==> {verdict_text(dr['adds_direction'], 'YES', 'NO')}")
    print("    Note: Kchan moves the VOLUME of the hydrograph directly, so a gain here is expected from volume alone.")
    print("    Whether it moves SHAPE beyond the dial is a different question; section 6 tests that.")
    if nsum:
        print()
        print("  NULL CONTROL for the same rules (details in section 7):")
        print(f"    with random Kchan (no effect), the share of the best 20% above {FALSE_USE_KCHAN:g}: median {nsum['share_above30_median']:.0%} "
              f"(95% of draws {nsum['share_above30_lo']:.0%} to {nsum['share_above30_hi']:.0%});  "
              f"your truth-A result: {rfa['share_best_above_30']:.1%}")
        print(f"    middle-80% span of the best 20% under the null: median {nsum['span_median']:.2f} decades "
              f"(95% {nsum['span_lo']:.2f} to {nsum['span_hi']:.2f});  your truth-B result: {rb['span_decades']:.2f}; "
              f"your truth-A result: {rfa['span_decades']:.2f}")
        if nsum["n_gain"]:
            print(f"    held-out gain under the null: median {nsum['gain_median']:+.3f}, largest {nsum['gain_max']:+.3f} "
                  f"({nsum['gain_exceeds_n']} of {nsum['n_gain']} draws above {DIRECTION_GAIN:g})")


def print_info(d, best20_rows, comp_rows, shares, partial_rows, ladder_rows, dial_rows, shape_rows):
    sec("4. Where the best runs sit (information)")
    print(f"  Kchan (mm/hr) among the best 20% and among good runs (KGE_2012 >= {GOOD_KGE:g}):")
    print(f"  {'truth':<10}{'set':<34}{'runs':>5}{'min':>9}{'p10':>9}{'p25':>9}{'median':>9}{'p75':>9}{'p90':>9}{'max':>9}")
    for r in best20_rows:
        print(f"  {r['truth']:<10}{r['set']:<34}{r['n_runs']:>5d}" + "".join(
            f"{fk(r[f'kchan_q{q:03d}']):>9}" for q in (0, 10, 25, 50, 75, 90, 100)))
    print("  Routing parameters among the same sets (min to max):")
    for r in best20_rows:
        print(f"    {r['truth']:<10}{r['set']:<34}" + "   ".join(
            f"{SHORT[p]} {r[f'{SHORT[p]}_min']:.3f}-{r[f'{SHORT[p]}_max']:.3f}" for p in ROUTING))

    sec("5. Volume versus shape (information)")
    print("  KGE_2012 = 1 - sqrt( (r-1)^2 + (gamma-1)^2 + (beta-1)^2 ). beta is the volume ratio; r and gamma describe shape and timing.")
    print("  Correlations with log Kchan (Spearman = across all runs; 'routing fixed' = partial r holding cv, flowexp, n fixed):")
    print(f"  {'truth':<10}{'quantity':<46}{'Spearman':>10}{'routing fixed':>15}{'median':>10}")
    for r in comp_rows:
        print(f"  {r['truth']:<10}{r['quantity']:<46}{fnum(r['spearman_with_logK'], 10, 2, True)}"
              f"{fnum(r['partial_r_with_logK_routing_fixed'], 15, 2, True)}{fnum(r['value_median'], 10, 3)}")
    print("\n  Which part of the score gap is largest (median share of the squared gap, volume / correlation / variability):")
    names = dict(all_runs="all runs", best_20="best 20%")
    for lab, sh in shares.items():
        for k, v in sh.items():
            print(f"    {lab:<10}{names.get(k, k):<10}  volume {v[0]:.0%}   correlation r {v[1]:.0%}   variability gamma {v[2]:.0%}")

    sec("5b. Each parameter alone and with the others held fixed (KGE_2012; information)")
    print("  'alone' = plain correlation with the score. 'held fixed' = partial correlation after the other three parameters are")
    print("  accounted for (straight-line model). Kchan's effect on a score that peaks at 201 is a curve, so its straight-line numbers")
    print("  understate it; the last column fits each parameter and its square.")
    print(f"  {'truth':<10}{'parameter':<10}{'alone r':>9}{'Spearman':>10}{'held fixed':>12}{'95% interval':>20}{'unique R2 (with square)':>26}")
    for r in partial_rows:
        print(f"  {r['truth']:<10}{r['parameter']:<10}{r['pearson_alone']:>+9.2f}{r['spearman_alone']:>+10.2f}{r['partial_r_linear']:>+12.2f}"
              f"{('[%+.2f, %+.2f]' % (r['partial_ci_lo'], r['partial_ci_hi'])):>20}{r['unique_r2_with_square']:>26.3f}")

    sec("5c. Model ladder: held-out R squared for KGE_2012 (information; the dial rows feed the direction rule)")
    for lab in dict.fromkeys(r["truth"] for r in ladder_rows):
        print(f"  Scored against {lab}:")
        for r in [x for x in ladder_rows if x["truth"] == lab]:
            print(f"    {r['model_label']:<44}{r['held_out_r2']:>8.3f}   (repeats {r['held_out_r2_min_rep']:.3f} to {r['held_out_r2_max_rep']:.3f})")
    print("  Fitted dial weights (standardized parameters; flowexp positive by convention). Series 121's dial was cv -0.36, flowexp +0.91, n +0.21:")
    print(f"    {'truth':<10}{'fit':<38}{'cv':>7}{'flowexp':>9}{'n':>7}{'log Kchan':>11}{'|cos| with Series 121':>24}")
    for r in dial_rows:
        print(f"    {r['truth']:<10}{r['fit']:<38}{r['w_cv']:>+7.2f}{r['w_flowexp']:>+9.2f}{r['w_n']:>+7.2f}"
              f"{fnum(r['w_logKchan'], 11, 2, True)}{fnum(r['abs_cosine_with_series121_dial'], 24, 2)}")

    sec("6. Does Kchan move SHAPE beyond the dial? (information)")
    print(f"  Same test as the direction rule, run metric by metric. Gains above {DIRECTION_GAIN:g} are flagged with *")
    print(f"  (the 0.02 yardstick is reused; it is not a pre-fixed rule for individual metrics).")
    print(f"  {'truth':<8}{'metric':<30}{'kind':<8}{'dial':>8}{'dial+curve':>12}{'gain':>9}")
    for r in shape_rows:
        print(f"  {r['truth']:<8}{r['metric_label']:<30}{r['kind']:<8}{r['r2_dial']:>8.3f}{r['r2_dial_curve']:>12.3f}"
              f"{r['gain']:>+9.3f} {'*' if r['gain_exceeds_yardstick'] else ''}")
    print("  Volume metrics (PBIAS, beta) are included for contrast. A flag on a shape or timing metric means Kchan carries")
    print("  information about the hydrograph's shape that the dial does not; it does not say why.")


def print_null(nc, nsum, rfa, rb):
    sec("7. NULL CONTROL: random Kchan attached to runs where Kchan did nothing")
    print(f"  {nc['n_runs']} runs from anchor {NULL_ANCHOR} in {nc['file']} (channel loss OFF, scored against truth A).")
    print(f"  Each draw gives every run a random Kchan, log-uniform from {BOX[P_K][0]:g} to {BOX[P_K][1]:g} mm/hr, "
          "unrelated to its score, and applies the Truth-A rules.")
    nd = nc["draws"]
    first = nd.iloc[0]
    print(f"  First draw (seed {nc['seed']}):  largest Kchan among good runs {fk(first['largest_kchan_good'])} mm/hr;  "
          f"share of best 20% above {FALSE_USE_KCHAN:g}: {first['share_best_above_30']:.0%};  "
          f"middle-80% span {first['span_decades']:.2f} decades")
    print(f"  Over {nsum['n_draws']} random draws:")
    print(f"    share of the best 20% above {FALSE_USE_KCHAN:g}: median {nsum['share_above30_median']:.0%}, "
          f"95% of draws {nsum['share_above30_lo']:.0%} to {nsum['share_above30_hi']:.0%};  "
          f"declared 'false use' in {nsum['false_use_fraction_of_draws']:.0%} of draws")
    print(f"    middle-80% span of the best 20%: median {nsum['span_median']:.2f} decades (95% {nsum['span_lo']:.2f} to {nsum['span_hi']:.2f}); "
          f"the recovery rule needs < {SPAN_MAX_FRACTION_OF_BOX * BOX_DECADES:.2f}")
    print(f"    largest Kchan among good runs: median {fk(nsum['largest_good_median'])} mm/hr;  "
          f"best-20% share in the lowest fifth of the box: median {nsum['lowest_fifth_median']:.0%}")
    if nsum["n_gain"]:
        print(f"    held-out gain from adding a random curve in log Kchan: median {nsum['gain_median']:+.3f}, largest "
              f"{nsum['gain_max']:+.3f} over {nsum['n_gain']} draws; {nsum['gain_exceeds_n']} above {DIRECTION_GAIN:g}")
    print("  HOW TO USE THIS: this is what the rules show when Kchan has NO effect. The more your real truth-A and truth-B")
    print("  numbers differ from these null ranges, the less they look like noise. The null has 100 runs, the sweep 400.")
    print(f"  Your truth-A share above {FALSE_USE_KCHAN:g}: {rfa['share_best_above_30']:.1%}  (null median {nsum['share_above30_median']:.0%}).  "
          f"Your truth-B span: {rb['span_decades']:.2f} decades (null median {nsum['span_median']:.2f}).")


# ======================================================================
# Figures (matplotlib; same look as the other Series 120/121 figures)
# ======================================================================
def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["font.family"] = "DejaVu Sans"
    return plt


def style_axes(ax, grid_axis="both"):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.grid(True, color=GRID, linewidth=0.8, axis=grid_axis)
    ax.set_axisbelow(True)
    ax.tick_params(colors=INK_SOFT, labelsize=8.5)


def save_fig(fig, path, dpi):
    import matplotlib.pyplot as plt
    fig.savefig(path, dpi=dpi, facecolor="white")
    plt.close(fig)
    print(f"  saved {Path(path).name}")


def figure_1(d, rb, rfa, out_path, dpi):
    plt = _plt()
    from matplotlib.lines import Line2D
    fig, axes = plt.subplots(1, 2, figsize=(14.0, 5.6), sharey=False)
    fig.patch.set_facecolor("white")
    specs = [("B", "kge_2012_B", rb, "Scored against truth B (channel loss ON, true Kchan = 201)"),
             ("A", "kge_2012_A", rfa, "Scored against truth A (NO channel loss)")]
    for ax, (suf, col, rr, title) in zip(axes, specs):
        style_axes(ax)
        top = top_set(d, col)
        ax.scatter(d["Kchan"], d[col], s=14, color=GREY_DOT, alpha=0.9, linewidths=0, zorder=3)
        ax.scatter(top["Kchan"], top[col], s=22, color=BLUE, linewidths=0, zorder=4)
        ax.axvspan(rr["kchan_p10"], rr["kchan_p90"], color=BLUE, alpha=0.10, zorder=1)
        ax.set_xscale("log")
        ymin = max(float(d[col].min()), -1.0)
        ax.set_ylim(ymin - 0.03, 1.03)
        if suf == "B":
            ax.axvline(TRUE_KCHAN, color=GOLD, linewidth=2.0, zorder=2)
            ax.text(TRUE_KCHAN, ymin - 0.02, " true Kchan 201", color=INK, fontsize=9, va="bottom")
        else:
            ax.axvline(FALSE_USE_KCHAN, color=INK_SOFT, linestyle=(0, (2, 3)), linewidth=1.0, zorder=2)
            ax.text(FALSE_USE_KCHAN, ymin - 0.02, " 30 mm/hr", color=INK_SOFT, fontsize=9, va="bottom")
            ax.axhline(GOOD_KGE, color=INK_SOFT, linestyle=(0, (2, 3)), linewidth=1.0, zorder=2)
            ax.text(BOX[P_K][1], GOOD_KGE + 0.008, "KGE_2012 = 0.95", color=INK_SOFT, fontsize=8.5, va="bottom", ha="right")
        ax.set_xlabel("Kchan (mm/hr, log scale)", fontsize=9.5, color=INK_SOFT)
        ax.set_ylabel(f"KGE_2012 against truth {suf}", fontsize=9.5, color=INK_SOFT)
        ax.set_title(title, fontsize=10.5, color=INK, loc="left")
        ax.legend(handles=[Line2D([0], [0], marker="o", linestyle="", color=GREY_DOT, markersize=6, label="all runs"),
                           Line2D([0], [0], marker="o", linestyle="", color=BLUE, markersize=6, label="best 20%"),
                           Line2D([0], [0], color=BLUE, alpha=0.3, linewidth=8,
                                  label=f"middle 80% of best 20%: {fk(rr['kchan_p10'])} to {fk(rr['kchan_p90'])}")],
                  loc="upper center", frameon=False, fontsize=8.5, ncol=3, bbox_to_anchor=(0.5, -0.14))
    fig.suptitle("Series 122: how well does each run score, by its Kchan?", fontsize=12.5, color=INK, y=0.995)
    fig.text(0.5, 0.005, "400 runs, Ks 7.0 and f 0.012 fixed, cv/flowexp/n/Kchan varied together. One storm (12 Aug 2014).",
             ha="center", va="bottom", fontsize=8.3, color=INK_SOFT)
    fig.tight_layout(rect=(0, 0.03, 1, 0.96))
    save_fig(fig, out_path, dpi)


def figure_2(d, nc, out_path, dpi):
    plt = _plt()
    from matplotlib.lines import Line2D
    fig, axes = plt.subplots(1, 2, figsize=(14.0, 5.4))
    fig.patch.set_facecolor("white")
    lo, hi = BOX[P_K]
    xs = np.logspace(math.log10(lo), math.log10(hi), 200)
    uniform = (np.log10(xs) - math.log10(lo)) / BOX_DECADES
    for ax, (suf, col, ttl) in zip(axes, [("B", "kge_2012_B", "Best 20% against truth B"), ("A", "kge_2012_A", "Best 20% against truth A")]):
        style_axes(ax)
        top = top_set(d, col)
        k = np.sort(top["Kchan"].to_numpy())
        cdf = np.arange(1, len(k) + 1) / len(k)
        if suf == "A" and nc is not None:
            rng = np.random.default_rng(nc["seed"])
            lo_l, hi_l = math.log10(lo), math.log10(hi)
            kg = nc["kge"]
            for _ in range(25):
                kk = np.sort(10 ** rng.uniform(lo_l, hi_l, len(kg))[np.argsort(-kg)[: int(math.ceil(TOP_FRAC * len(kg)))]])
                ax.step(kk, np.arange(1, len(kk) + 1) / len(kk), where="post", color=GREY_DOT, linewidth=0.8, alpha=0.5, zorder=2)
        ax.plot(xs, uniform, color=INK_SOFT, linestyle=(0, (2, 3)), linewidth=1.2, zorder=3)
        ax.step(k, cdf, where="post", color=BLUE if suf == "B" else ORANGE, linewidth=2.2, zorder=5)
        ax.set_xscale("log")
        ax.set_ylim(0, 1.02)
        if suf == "B":
            ax.axvline(TRUE_KCHAN, color=GOLD, linewidth=2.0, zorder=4)
            ax.text(TRUE_KCHAN, 0.02, " true Kchan 201", fontsize=9, color=INK)
        else:
            ax.axvline(FALSE_USE_KCHAN, color=INK_SOFT, linestyle=(0, (2, 3)), linewidth=1.0, zorder=4)
            ax.text(FALSE_USE_KCHAN, 0.02, " 30 mm/hr", fontsize=9, color=INK_SOFT)
        ax.set_xlabel("Kchan (mm/hr, log scale)", fontsize=9.5, color=INK_SOFT)
        ax.set_ylabel("share of the best 20% with Kchan at or below x", fontsize=9.5, color=INK_SOFT)
        ax.set_title(ttl, fontsize=10.5, color=INK, loc="left")
        handles = [Line2D([0], [0], color=BLUE if suf == "B" else ORANGE, linewidth=2.2, label="the 80 best runs"),
                   Line2D([0], [0], color=INK_SOFT, linestyle=(0, (2, 3)), linewidth=1.2, label="no information (Kchan unrelated to score)")]
        if suf == "A" and nc is not None:
            handles.append(Line2D([0], [0], color=GREY_DOT, linewidth=0.9, label="null control: random Kchan, 25 draws"))
        ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=8.5)
    fig.suptitle("Series 122: where do the best runs put Kchan?", fontsize=12.5, color=INK, y=0.995)
    fig.text(0.5, 0.005, "A curve that rises steeply around one value means the best runs agree on Kchan; a curve along the dotted line means they do not.",
             ha="center", va="bottom", fontsize=8.3, color=INK_SOFT)
    fig.tight_layout(rect=(0, 0.03, 1, 0.96))
    save_fig(fig, out_path, dpi)


def figure_3(d, out_path, dpi):
    plt = _plt()
    from matplotlib.lines import Line2D
    comps = [("kge_beta", "beta: volume ratio (1 = right volume)"), ("kge_r", "r: correlation (1 = right timing)"),
             ("kge_gamma", "gamma: variability ratio (1 = right shape)")]
    fig, axes = plt.subplots(2, 3, figsize=(15.0, 8.2))
    fig.patch.set_facecolor("white")
    for ri, (suf, ttl) in enumerate((("B", "against truth B"), ("A", "against truth A"))):
        top = top_set(d, f"kge_2012_{suf}")
        for ci, (key, lab) in enumerate(comps):
            ax = axes[ri][ci]
            style_axes(ax)
            col = f"{key}_{suf}"
            ax.scatter(d["Kchan"], d[col], s=12, color=GREY_DOT, alpha=0.9, linewidths=0, zorder=3)
            ax.scatter(top["Kchan"], top[col], s=18, color=BLUE, linewidths=0, zorder=4)
            ax.axhline(1.0, color=INK, linewidth=0.8, zorder=2)
            if suf == "B":
                ax.axvline(TRUE_KCHAN, color=GOLD, linewidth=1.8, zorder=2)
            ax.set_xscale("log")
            ax.set_xlabel("Kchan (mm/hr, log)", fontsize=9, color=INK_SOFT)
            ax.set_ylabel(lab, fontsize=8.8, color=INK_SOFT)
            ax.set_title(f"{lab.split(':')[0]} {ttl}", fontsize=10.5, color=INK, loc="left")
            if ri == 0 and ci == 0:
                ax.legend(handles=[Line2D([0], [0], marker="o", linestyle="", color=GREY_DOT, markersize=6, label="all runs"),
                                   Line2D([0], [0], marker="o", linestyle="", color=BLUE, markersize=6, label="best 20%"),
                                   Line2D([0], [0], color=GOLD, linewidth=1.8, label="true Kchan 201")],
                          loc="lower left", frameon=False, fontsize=8.3)
    fig.suptitle("Series 122: does Kchan move volume, or shape, or both?", fontsize=12.5, color=INK, y=0.995)
    fig.text(0.5, 0.005, "The score is 1 minus the combined gap of the three quantities from 1. The horizontal line marks 1.",
             ha="center", va="bottom", fontsize=8.3, color=INK_SOFT)
    fig.tight_layout(rect=(0, 0.02, 1, 0.965))
    save_fig(fig, out_path, dpi)


def figure_4(ladder_rows, dir_a, dir_b, nc, out_path, dpi):
    plt = _plt()
    from matplotlib.lines import Line2D
    fig, axes = plt.subplots(1, 3, figsize=(16.0, 5.2), gridspec_kw=dict(width_ratios=[1, 1, 0.9]))
    fig.patch.set_facecolor("white")
    keys = ["dial", "curve", "dial_curve", "dial4", "full4"]
    colors = {"dial": GREY_DOT, "curve": GREY_DOT, "dial_curve": BLUE, "dial4": GREY_DOT, "full4": "#8d8c85"}
    for ax, (lab, ttl) in zip(axes[:2], (("truth B", "Truth B: held-out R squared"), ("truth A", "Truth A: held-out R squared"))):
        style_axes(ax, grid_axis="x")
        rows = {r["model"]: r for r in ladder_rows if r["truth"] == lab}
        ypos = np.arange(len(keys))[::-1]
        vals = [rows[k]["held_out_r2"] for k in keys if k in rows]
        for yv, k in zip(ypos, [k for k in keys if k in rows]):
            ax.barh(yv, rows[k]["held_out_r2"], 0.62, color=colors[k], zorder=3)
            ax.text(max(rows[k]["held_out_r2"], 0) + 0.01, yv, f"{rows[k]['held_out_r2']:.3f}", va="center", fontsize=9, color=INK)
        ax.set_yticks(ypos[: len(vals)])
        ax.set_yticklabels([MODEL_LABELS[k] for k in keys if k in rows], fontsize=8.8)
        ax.set_xlim(min(0, min(vals) - 0.05), 1.12)
        ax.axvline(0, color=INK, linewidth=0.8)
        ax.set_xlabel("held-out R squared (KGE_2012)", fontsize=9, color=INK_SOFT)
        ax.set_title(ttl, fontsize=10.5, color=INK, loc="left")
    ax = axes[2]
    style_axes(ax, grid_axis="x")
    g = nc["draws"]["gain"].dropna().to_numpy() if (nc is not None and nc["draws"]["gain"].notna().any()) else np.array([])
    real = [("truth B", dir_b["gain"], BLUE), ("truth A", dir_a["gain"], ORANGE)]
    xlo = min([-0.02] + ([g.min() - 0.01] if len(g) else []))
    xhi = max([2.5 * DIRECTION_GAIN] + ([g.max() + 0.01] if len(g) else []))
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    if len(g):
        yj = np.random.default_rng(0).uniform(0.30, 0.55, len(g))
        ax.scatter(g, yj, s=36, color=GREY_DOT, linewidths=0, zorder=3)
        ax.text(0.02, 0.60, f"null control: {len(g)} random draws", transform=ax.transAxes, fontsize=8.8, color=INK_SOFT)
    ax.axvline(DIRECTION_GAIN, color=INK_SOFT, linestyle=(0, (2, 3)), linewidth=1.0, zorder=2)
    ax.text(DIRECTION_GAIN, 0.97, f" rule: {DIRECTION_GAIN:g}", fontsize=8.8, color=INK_SOFT, va="top")
    ypos = 0.85
    for name, val, col in real:
        if xlo <= val <= xhi:
            ax.axvline(val, color=col, linewidth=2.2, zorder=4)
            ax.text(val, ypos, f" {name}: {val:+.3f}", color=INK, fontsize=9, va="top")
        else:
            ax.annotate(f"{name}: {val:+.3f}  (off this scale)", xy=(xhi, ypos - 0.03), xytext=(xhi - 0.012, ypos - 0.03),
                        ha="right", va="center", fontsize=9, color=INK,
                        arrowprops=dict(arrowstyle="->", color=col, linewidth=2.0))
        ypos -= 0.08
    ax.set_xlabel("held-out gain from adding a curve in log Kchan", fontsize=9, color=INK_SOFT)
    ax.set_title("Gain vs the null control (zoomed in)", fontsize=10.5, color=INK, loc="left")
    fig.suptitle("Series 122: does Kchan add a direction beyond the dial?", fontsize=12.5, color=INK, y=0.995)
    fig.text(0.5, 0.005, "Held-out R squared from repeated cross-validation. Kchan moves volume directly, so a gain is expected; see section 6 for shape.",
             ha="center", va="bottom", fontsize=8.3, color=INK_SOFT)
    fig.tight_layout(rect=(0, 0.03, 1, 0.96))
    save_fig(fig, out_path, dpi)


# ======================================================================
# Main
# ======================================================================
def parse_args(argv=None):
    ap = argparse.ArgumentParser(description="Series 122 analysis: recovery and false use of Kchan (reads results; writes a new folder).")
    ap.add_argument("--dry_run", action="store_true")
    ap.add_argument("--results", type=Path, default=None)
    ap.add_argument("--null_csv", type=Path, default=None)
    ap.add_argument("--out_dir", type=Path, default=None)
    ap.add_argument("--no_null", action="store_true")
    ap.add_argument("--cv_reps", type=int, default=2)
    ap.add_argument("--cv_folds", type=int, default=5)
    ap.add_argument("--null_reps", type=int, default=500)
    ap.add_argument("--null_cv_draws", type=int, default=12)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--no_plots", action="store_true")
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args(argv)
    if args.quick:
        args.cv_reps, args.null_cv_draws = 1, 3
    if not (1 <= args.cv_reps <= 20):
        ap.error("--cv_reps must be between 1 and 20.")
    if not (3 <= args.cv_folds <= 20):
        ap.error("--cv_folds must be between 3 and 20.")
    if not (10 <= args.null_reps <= 20000):
        ap.error("--null_reps must be between 10 and 20000.")
    if not (0 <= args.null_cv_draws <= 200):
        ap.error("--null_cv_draws must be between 0 and 200.")
    if not (72 <= args.dpi <= 600):
        ap.error("--dpi must be between 72 and 600.")
    return args


def print_dry_run(pre, paths, args):
    sec("DRY RUN: nothing was scored or written")
    print(f"  {len(pre['df'])} runs would be scored against truth A and truth B from {len(pre['kept_paths'])} kept outlet files.")
    print(f"  Event window {pre['ev_start']} to {pre['ev_end']}.")
    print(f"  Output folder (created only on a real run): {pre['out_dir']}")
    print("  Files that would be written there: " + ", ".join(v for k, v in FILES.items() if not (args.no_plots and k.startswith("fig"))))
    print(f"  Cross-validation: {args.cv_reps} repeat(s) x {args.cv_folds} folds.  Null control: "
          + ("on" if pre["null_ok"] else "off") + f" ({args.null_reps} cheap draws, {args.null_cv_draws} with cross-validation).")
    print("\nReady:  python -u analyze_series122_kchan.py 2>&1 | tee series122_analysis.log")


def run(args, tee):
    paths = locate(args)
    pre = preflight(paths, args)
    if args.dry_run:
        print_dry_run(pre, paths, args)
        return 0

    t0 = time.time()
    scored, b_vs_a = score_all(pre)
    gate_tab = gate(pre, scored, b_vs_a)

    d = pre["df"].merge(scored, on="run_id", how="inner")
    d["Kchan"] = d[P_K].astype(float)
    d["logK"] = np.log10(d["Kchan"])
    d = d.dropna(subset=["kge_2012_A", "kge_2012_B"]).reset_index(drop=True)
    print(f"\n  {len(d)} runs have a finite KGE_2012 against both truths and go into the analysis.")

    # ---- the rules ----
    rb = rule_recovery(d, "kge_2012_B")
    rfa = rule_false_use(d, "kge_2012_A")
    print("\n  (running the cross-validation for the direction rule ...)")
    dir_b = rule_direction(d, "kge_2012_B", args)
    dir_a = rule_direction(d, "kge_2012_A", args)

    # the null control is printed as section 7, but section 3 compares with it, so it is computed first
    nc = nsum = None
    if pre["null_ok"]:
        print("\n  (computing the null control first, because section 3 compares with it; it is shown in section 7 ...)")
        nc = null_control(paths["null_csv"], args)
        if nc is not None:
            nsum = summarize_null(nc)

    print_rules(rb, rfa, dir_a, dir_b, nsum, d)

    # ---- information ----
    best20_rows = best20_table(d, "kge_2012_B", "truth B") + best20_table(d, "kge_2012_A", "truth A")
    comp_rows, shares = [], {}
    for lab, suf in (("truth B", "B"), ("truth A", "A")):
        rows, sh = component_table(d, lab, suf)
        comp_rows += rows
        shares[lab] = dict(all_runs=sh["all_runs"], best_20=sh["best_20"])
    partial_rows = partial_tables(d, "B", "truth B") + partial_tables(d, "A", "truth A")
    print("\n  (fitting the model ladder and the dials ...)")
    ladder_rows, dial_rows = [], []
    for lab, suf in (("truth B", "B"), ("truth A", "A")):
        lr, dr, _ = ladder_and_dials(d, suf, lab, args)
        ladder_rows += lr
        dial_rows += dr
    print("  (testing shape metrics one by one ...)")
    shape_rows = shape_beyond_dial(d, "B", "truth B", args) + shape_beyond_dial(d, "A", "truth A", args)
    print_info(d, best20_rows, comp_rows, shares, partial_rows, ladder_rows, dial_rows, shape_rows)

    if nc is not None:
        print_null(nc, nsum, rfa, rb)

    sec("8. CAUTIONS")
    for line in ("One synthetic truth pair, one storm (12 Aug 2014, 1.0x); Ks 7.0 and f 0.012 are fixed at their true values (Stage 1).",
                 "A good Kchan estimate here does not show Kchan is identifiable with real data, noise, other storms, or unknown Ks and f.",
                 "The 400 runs are one Latin hypercube draw (seed 301); a different draw would shift every percentile a little.",
                 "The reading rules were fixed before any result was seen. Everything marked 'information' is there to help read them, not to replace them.",
                 "The tRIBS mechanism behind any pattern is a hypothesis to check in the model code, not something these numbers show."):
        print("  * " + line)

    # ---- tables ----
    sec("Files written to " + str(pre["out_dir"]))
    out = pre["out_dir"]
    out.mkdir(parents=True, exist_ok=True)
    rules_rows = [
        dict(rule="recovery", truth="B", **{k: v for k, v in rb.items()}),
        dict(rule="false_use", truth="A", **{k: v for k, v in rfa.items()}),
        dict(rule="direction", truth="B", **dir_b),
        dict(rule="direction", truth="A", **dir_a),
    ]
    if nsum:
        rules_rows.append(dict(rule="null_control_summary", truth="A (Series 121 runs)", **nsum))
    keep = ["run_id"] + FOUR + [c for c in scored.columns if c != "run_id"]
    tables = [("scored", d[keep]), ("rules", pd.DataFrame(rules_rows)), ("best20", pd.DataFrame(best20_rows)),
              ("gate", gate_tab), ("components", pd.DataFrame(comp_rows)), ("partial", pd.DataFrame(partial_rows)),
              ("ladder", pd.DataFrame(ladder_rows)), ("dial", pd.DataFrame(dial_rows)), ("shape", pd.DataFrame(shape_rows))]
    if nc is not None:
        tables.append(("null", nc["draws"]))
    for key, tab in tables:
        if tab.empty:
            print(f"  (skipped {FILES[key]}: nothing to write)")
            continue
        tab.to_csv(out / FILES[key], index=False, float_format="%.8g")
        print(f"  saved {FILES[key]}  ({len(tab)} rows)")

    if not args.no_plots:
        print("\nDrawing the figures...")
        for name, fn in (("fig1", lambda: figure_1(d, rb, rfa, out / FILES["fig1"], args.dpi)),
                         ("fig2", lambda: figure_2(d, nc, out / FILES["fig2"], args.dpi)),
                         ("fig3", lambda: figure_3(d, out / FILES["fig3"], args.dpi)),
                         ("fig4", lambda: figure_4(ladder_rows, dir_a, dir_b, nc, out / FILES["fig4"], args.dpi))):
            try:
                fn()
            except Exception as exc:                                    # noqa: BLE001  a figure problem must not lose the tables
                print(f"  (figure {FILES[name]} skipped: {type(exc).__name__}: {exc})")
    try:
        (out / FILES["report"]).write_text(tee.buf.getvalue())
        print(f"\n  saved {FILES['report']}  (a copy of everything printed above)")
    except Exception as exc:                                             # noqa: BLE001
        print(f"\n  (could not save {FILES['report']}: {exc})")
    print(f"\nDone in {time.time() - t0:.0f} s. This script only reads your results; it changed none of them.")
    return 0


def main(argv=None):
    args = parse_args(argv)
    tee = Tee(sys.stdout)
    sys.stdout = tee
    print(f"Series 122 analysis   python {sys.version.split()[0]}, pandas {pd.__version__}, numpy {np.__version__}   "
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
