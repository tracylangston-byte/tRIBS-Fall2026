"""
rescore_real_gauge_110.py   (version 2)
=======================================
Series 110 -- STAGE 0 of the "does channel loss improve the fit to REAL data?" plan.
RE-SCORE the 500 stored Aug 12 2014 runs (250 cc ON, 250 cc OFF control) against the REAL
SMF gauge hydrograph. No tRIBS, no new simulations.

WHAT CHANGED SINCE VERSION 1
-----------------------------
Version 1 stopped at gate 1. The reason was the gauge, not the script. The Discharge sheet is not
a regular 5-minute record. It is written when the flow changes (every 1-3 minutes during the
flood, with gaps of 9 to 39 minutes on the recession) plus a "heartbeat" zero about every six
hours. Nothing is recorded between 21:11 and 07:59. Version 1 also expected a peak of about
122 cfs; the file says 1586 cfs (44.9 m3/s), the same order as the model's synthetic hydrographs,
so the 122 figure was wrong and that check is gone. It is replaced by a scale check against the
stored simulations. Version 2 scores the record two ways ("arms"), described next.

WHAT THIS DOES (plain language)
--------------------------------
Every one of the 500 runs has a stored hydrograph (the "Simulated" column of its
*_compare_obs_sim.csv). Until now each was scored against a SYNTHETIC truth. A simulation
does not depend on what it is scored against -- only the score does. So this script keeps
every stored "Simulated" series untouched, swaps in the REAL gauge record as "Observed",
and recomputes the same metrics with the same formulas the scorer uses.

TWO WAYS TO SCORE A GAUGE THAT HAS GAPS ("arms")
-------------------------------------------------
  as-is   Exactly what the scorer does in real-gauge mode: Date + Time joined, cfs x 0.0283168
          to m3/s, .resample("5min").mean(), empty bins dropped, clipped to the event window
          2014-08-12 16:00 to 2014-08-13 12:00. With this gauge that leaves 23 of the 241
          five-minute points (22 flood bins and one stray zero). The score is BLIND to any
          simulated flow in the other 218 bins: the false early flow from the SMPHQ rainfall
          (the double peak) and the dry hours are never compared with anything.
  filled  Keeps those same 23 values and fills the other 218 bins from the record itself.
          Between two records at most --max_gap_min apart (default 60) it draws a straight
          line. Across a longer silence it HOLDS the last reading (a logger that writes only
          on change writes nothing when nothing changes), for at most --max_hold_hr (default
          24 h); a bin with no reading that recent is left out and reported. All 241 points
          are scored, which is how every synthetic test in Series 100 and 110 was scored,
          so its KGE values and your margin are on the same footing as those. The assumption
          behind the hold is that a silent logger means zero change, not a dead logger; the
          file cannot tell those two apart.
The two arms use the same stored simulations and the same metric code. Only the "Observed"
column differs, and the as-is points are a strict subset of the filled ones, so any difference
between the arms comes from the fill alone.

WHICH ARM COUNTS (you must say, before the run)
------------------------------------------------
With --gauge_mode both (the default) both arms are scored and reported, but you must name the
PRIMARY arm with --primary_arm. That arm's reading is the pre-registered answer; the other is a
sensitivity check, and the report says whether the two agree. Like the margin, this has to be
fixed BEFORE any real-data result is shown (handoff section 4).

THREE QUESTIONS, ASKED OF EACH ARM, within this one event and this one (Ks, f, cc) box
------------------------------------------------------------------------------------------
  1. Best KGE_2012, and the mean of the 10 best, for cc ON versus cc OFF.
  2. Where do the best cc-ON runs put cc: lowest third of the range, at the 30 mm/hr floor,
     or interior? (Same binomial test as v4 section 3b: chance puts 1/3 in the lowest third.)
  3. The paired difference ON minus OFF at the same (Ks, f, cc) point, against log cc. Because
     the OFF run at a point has the same Ks and f, the pair controls for Ks and f exactly.
     Does fit get worse, flat or better as cc rises?

THE DECISION MARGIN (you must give one)
-----------------------------------------
Handoff section 4 says to fix the decision rule BEFORE looking at any real-data result. So the
full run refuses to start without --margin (in KGE_2012 units: how much better ON must be than
OFF before you call it a real gain). --check_only and --observed_only do not need it and show no
fit. For scale: in the synthetic work the ON-OFF gap was under 0.03 in 13 of 27 truths, and real
gauge error is not quantified.

SELF-CHECKS ("gates") -- run before anything is scored
--------------------------------------------------------
  GATE 1  The real series. Workbook read, duplicate timestamps, values, the peak (checked against
          the scale of the stored simulations, not against a number from the handoff), the
          logger's cadence and silences, how many 5-minute points the as-is arm keeps, and the
          workbook's checksum. 1b checks the filled series: it equals the as-is series at every
          bin that has a record, its volume agrees with the raw records' own trapezoid sum, and
          it reports how many bins were filled by line, by hold, or left uncovered.
  GATE 2  Scorer fidelity. Feeding the SYNTHETIC cc-OFF truth (the single .qout in
          synth_truth/) through this script's reading and metric code must reproduce every
          stored "Observed" column (2A) and every stored metric (2B) of all 500 runs. This
          is rescore_cc_truths_110.py's gate A and B, reused. It also confirms the metric code
          on the full 241-point grid the filled arm uses.
  GATE 3  Gauge-mode fidelity of the AS-IS path. One stored parameter point is re-run through
          the ACTUAL scorer in real-gauge mode (needs one tRIBS run, about 1 minute; this
          script never starts tRIBS). Point --gate3_compare at that run's
          *_compare_obs_sim.csv and --gate3_metrics at its *_metrics_summary.csv. Checked: the
          scorer's real Observed equals this script's as-is series (3a); the re-run's Simulated
          equals the stored one (3b); the scorer's metrics equal this script's (3c). Both arms
          read the workbook the same way, so this also vouches for the time stamps of the
          filled arm. Without --gate3_compare the results are labelled PROVISIONAL.
  GATE 4  Grid. All 500 stored Simulated series sit on one regular 5-minute grid, and every
          real 5-minute bin that is scored is on that grid.
If a gate fails the script stops and says why. (--drop_bad_runs excludes individual runs that
fail 2A, 2B or 4 and continues. Gate 1 and gate 3 failures always stop.)

BLIND SHARE (printed with the gates, no fit involved)
------------------------------------------------------
For each series: the median share of each run's simulated volume that falls in bins the as-is
arm can see, and in bins it cannot see, split into before the flood, inside the flood (holes in
the record) and after the flood. "Flood" = from the first to the last bin where the filled
series is above zero. A large share before the flood is the SMPHQ false-early-flow problem.

WHAT YOU GET   (calibration_work/03_comparisons/summary_tables/real_gauge_110/)
---------------------------------------------------------------------------------
  top folder
    real_gauge_observed_110.csv     the real series on the 5-minute grid (window): as-is value
                                    (blank where the scorer drops the bin), filled value, how each
                                    bin was filled (record / linear / hold / uncovered), m3/s and cfs
    real_gauge_blind_share_110.csv  the blind-share table
    real_gauge_arm_agreement_110.csv  per arm: best ON and OFF, gap, handoff row (full run only)
    real_gauge_gates_110.csv        every check with its status
    fig_real_gauge_observed_110.png raw records, as-is bins and the filled series, whole window and
                                    zoom on the flood
    PROVENANCE_real_gauge_110.json  checksums, versions, gate results, settings, margin, arm
  one folder per arm: filled/ and asis/   (full run only)
    real_gauge_long_110.csv         every run: all metrics vs the real gauge, plus its series
    real_gauge_summary_110.csv      one row per series (best run, top-N mean, feasible count) and
                                    an ON-minus-OFF row
    real_gauge_top_runs_110.csv     the best runs of each series with parameters and metrics
    real_gauge_cc_location_110.csv  where the best ON runs put cc (counts per third, floor,
                                    binomial p-values), for each --top_n
    real_gauge_paired_110.csv       the ON-OFF pairs: KGE_2012 of each, difference, log cc
    real_gauge_paired_stats_110.csv the difference against log cc: rank correlation, slope with a
                                    bootstrap interval, medians per cc third (all pairs, and the
                                    region where the OFF runs fit best)
    fig_real_gauge_hydrographs_110.png    real gauge with the best ON and best OFF run
    fig_real_gauge_kge_vs_cc_110.png      KGE_2012 against cc for all runs, best ON runs marked
    fig_real_gauge_paired_diff_110.png    ON-OFF difference against cc

HOW TO READ THE RESULT (rules of thumb, mine -- change the constants below before you run)
-------------------------------------------------------------------------------------------
  * The script prints which row of the handoff section 4 table the numbers fit. That uses your
    margin plus three rules of mine: a third of the cc range counts as "enriched" if the
    binomial p is below ALPHA; the best runs' cc counts as "stable" if their spread is at most
    STABLE_SPAN of the box (log scale); "the floor" is within FLOOR_TOL_MMHR of 30 mm/hr.
  * It also applies the same reading to the top-N MEAN gap and tells you if the two disagree.
    If they do, treat the answer as unsettled.
  * It then compares the two arms. If the arms point to different rows, the primary arm's row
    is still the pre-registered answer, but say that the other arm disagreed.
  * The "good region" in the paired analysis is the top --good_frac of the 250 points by how
    well the cc-OFF run fits. It is chosen from the control alone, where cc has no effect, so
    the choice cannot be influenced by cc.
  * A binomial test treats the best runs as independent draws. Neighbouring points in the
    sweep are not independent, so read p-values as a guide, not a verdict.

SAFETY
-------
  - Read-only on every existing file. Writes only inside its own output folder (refuses a folder
    that is, or overlaps, synth_truth/, csv_exports/ or any input file).
  - Never imports the builder, never touches current_run_config.json, never starts tRIBS.
  - Run membership comes only from the run_id rows of the two results CSVs; it never globs
    the compare folder (csv_exports/ holds extra _s2 files from earlier attempts).
  - Do not run it WHILE a sweep is still rewriting the results CSVs it reads.
  - Use the same environment as the runs: python 3.11.16, pandas 3.0.5, numpy 2.4.6.
    (.interpolate(method="time") differs between pandas 2 and 3; gate 2A will tell.)

USAGE (run from the lab/ directory; keep rescore_cc_truths_110.py, rescore_truth_location_110.py
       and run_sensitivity_single_interp_tribs6.py in the same folder)
-----------------------------------------------------------------------------------------------
    python rescore_real_gauge_110.py --check_only                  # gates, blind share; writes nothing
    python rescore_real_gauge_110.py --observed_only               # also writes the observed CSV + figure
    python rescore_real_gauge_110.py --margin 0.05 --primary_arm filled     # the full run (provisional)
    python rescore_real_gauge_110.py --margin 0.05 --gauge_mode filled      # one arm only
    python rescore_real_gauge_110.py --check_only \\
        --gate3_compare PATH/real_check_compare_obs_sim.csv \\
        --gate3_metrics PATH/real_check_metrics_summary.csv        # gate 3 only
    python rescore_real_gauge_110.py --margin 0.05 --primary_arm filled \\
        --gate3_compare ... --gate3_metrics ...                    # the full run with gate 3

CAVEATS (handoff section 5)
----------------------------
One event. The SMPHQ rainfall artifact (a double peak the real hydrograph lacks) is in every
run, ON and OFF alike, and parameters can compensate for it. Routing is pinned at the
synthetic-truth values. The (Ks, f, cc) box clipped the best-fit ridge in the cc-OFF-truth
analysis. Gauge error is not quantified, and the gauge clock and the model clock are assumed
to agree. In real-gauge scoring each 5-minute bin mean is labelled by its start while the
stored simulation is the value AT that time, so the real series sits about 2.5 minutes "early"
(the synthetic truths are read as instantaneous values, so they have no such offset). That is
small next to the 15-20 minute peak timing offset seen earlier. A small ON-OFF difference here
is not evidence either way.
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
    import rescore_cc_truths_110 as RS
    import rescore_truth_location_110 as TL
except ImportError:
    sys.exit("This script needs rescore_cc_truths_110.py and rescore_truth_location_110.py in the same folder "
             "(it reuses their reading, gate and scoring code). Put them in lab/ and run from there.")

_NEEDED_RS = ["load_results", "load_all_compares", "load_compare", "gate_a", "gate_b", "read_truth_5min",
              "compute_metrics", "import_phase_fn", "md5_of", "close", "RTOL", "ATOL", "STAGE2_NAME",
              "CONTROL_NAME", "SER_ON", "SER_OFF", "BOX", "FEASIBLE_PBIAS_PCT", "DEFAULT_EXPECT_N",
              "SCALAR_METRICS", "PHASE_METRICS", "OBS_RTOL", "OBS_ATOL"]
_missing = [n for n in _NEEDED_RS if not hasattr(RS, n)]
if _missing:
    sys.exit("rescore_cc_truths_110.py in this folder is an older version (missing: %s). Replace it with the "
             "current copy." % ", ".join(_missing))
if not hasattr(TL, "score_truth"):
    sys.exit("rescore_truth_location_110.py in this folder is an older version (missing: score_truth). Replace it "
             "with the current copy.")

SER_ON, SER_OFF = RS.SER_ON, RS.SER_OFF

# ------------------------------------------------------------------
# Constants. The first block must agree with the scorer; the second block is mine.
# ------------------------------------------------------------------
GAUGE_REL = ("init_data", "met", "SMF_Observations_1993-2025.xlsx")   # relative to the project root (lab/'s parent)
GAUGE_SHEET = "Discharge"
GAUGE_SKIPROWS = 6
CFS_TO_CMS = 0.0283168
EVENT_START = "2014-08-12 16:00"
EVENT_END = "2014-08-13 12:00"
GRID_FREQ = "5min"
EXPECT_PANDAS_MINOR = "3.0"
OUT_DIRNAME = "real_gauge_110"
PCOLS = ["Ks_mult", "f_RS_abs", "channelconductivity_mmhr"]

# --- rules of thumb of mine (edit before running if you disagree) ---------------------------
ALPHA = 0.05            # a third of the cc range is "enriched" among the best runs if the binomial p < ALPHA
STABLE_SPAN = 0.5       # best-runs cc spread, as a share of the log cc box, at or below which cc is "stable"
FLOOR_TOL_MMHR = 2.0    # "at the floor" = within this many mm/hr of the lowest cc in the box (as in v4 3b)
MIN_GOOD = 10           # fewest points in the "good region"
SCALE_BAND = (0.2, 5.0)     # real peak / median stored simulated peak outside this band => WARN (units? wrong event?)
DEFAULT_MAX_GAP_MIN = 60.0  # filled arm: straight line between records at most this many minutes apart
DEFAULT_MAX_HOLD_HR = 24.0  # filled arm: hold the last reading at most this long across a silence
MIN_ASIS_POINTS = 5         # fewer scorable as-is points than this is not a score at all (FAIL)
VOL_WARN_FRAC = 0.10        # filled volume vs the raw records' trapezoid sum: WARN beyond this
VOL_FAIL_FRAC = 0.30        # ... FAIL beyond this
BLIND_WARN_SHARE = 0.25     # NOTE when this much of the median run's volume is invisible to the as-is arm
GATE3_SIM_REL = 1e-6    # re-run Simulated must equal the stored one to this fraction of its peak
ARM_TEXT = {"asis": "as-is: the scorer's own reading, only the 5-minute bins that hold a record",
            "filled": "filled: all 5-minute bins, gaps filled from the record (line, then hold)"}
DROPPABLE = ("GATE 2A synthetic Observed", "GATE 2B stored metrics", "compare CSVs")   # failures --drop_bad_runs may exclude

READINGS = {
    1: "ON about equal to OFF, best cc at the floor: the outlet does not need channel loss to fit this event. "
       "cc stays valuable as a streambed-infiltration sensitivity tool, not as a calibration improvement.",
    2: "ON better than OFF by more than the margin, best cc interior and stable: evidence the data want channel "
       "loss. Go to Stage 2 to see whether it transfers.",
    3: "ON better, but best cc at a box edge or unstable between runs: extra flexibility, not process. Do not "
       "claim it.",
    4: "ON worse than OFF: the extra dimension costs search efficiency; check with Stage 1 before drawing any "
       "conclusion.",
    5: "ON about equal to OFF, but the best cc is NOT at the floor: no row of the handoff table fits exactly. "
       "Read the cc spread and the paired difference before saying anything.",
}

# palette (same as analyze_cc_pca_110.py)
C_KS, C_F, C_CC = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID, SURFACE, NULLC = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb", "#8a8984"


# ------------------------------------------------------------------
# Small helpers
# ------------------------------------------------------------------
class Checks(object):
    """A list of (status, check, detail) rows: PASS, WARN, FAIL, NOTE."""

    def __init__(self):
        self.rows = []

    def add(self, status, name, detail):
        self.rows.append({"status": status, "check": name, "detail": detail})

    def extend(self, other):
        self.rows.extend(other.rows)

    def n(self, status):
        return sum(1 for r in self.rows if r["status"] == status)

    def print_rows(self, rows=None):
        for r in (self.rows if rows is None else rows):
            print("  %-5s %-34s %s" % (r["status"], r["check"], r["detail"]))


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


def fmt_p(p):
    if p is None or not np.isfinite(p):
        return "n/a"
    return "<0.001" if p < 0.001 else "%.3f" % p


def fcc(v):
    return "n/a" if (v is None or not np.isfinite(v)) else ("%.0f" % v if v >= 100 else "%.1f" % v)


def cc_box_logs():
    lo, hi = RS.BOX["channelconductivity_mmhr"][:2]
    return float(np.log10(lo)), float(np.log10(hi))


def cc_unit(cc):
    """Position of a cc value in the log cc box: 0 = floor, 1 = ceiling."""
    lo, hi = cc_box_logs()
    return (np.log10(np.asarray(cc, dtype=float)) - lo) / (hi - lo)


# ------------------------------------------------------------------
# The real record, read exactly as the scorer reads it in real-gauge mode
# ------------------------------------------------------------------
def read_gauge_workbook(path):
    """Same calls as run_sensitivity_single_interp_tribs6.py, REAL GAUGE MODE."""
    obs_df = pd.read_excel(path, sheet_name=GAUGE_SHEET, skiprows=GAUGE_SKIPROWS)
    obs_df["datetime"] = pd.to_datetime(obs_df["Date"].astype(str) + " " + obs_df["Time"].astype(str))
    obs_df.set_index("datetime", inplace=True)
    obs_df["Observed_CMS"] = obs_df["cfs"] * CFS_TO_CMS
    return obs_df


def gauge_to_grid(obs_df):
    """The scorer's gauge-mode resampling: bin MEAN (not interpolation). Empty bins are NaN."""
    return obs_df["Observed_CMS"].resample(GRID_FREQ).mean()


def check_real_series(obs_df, real_5min, xlsx_path, min_points, sim_peak_med=None):
    """GATE 1 (the as-is series). Returns (Checks, info dict, real_valid Series on the window grid, real_grid Series incl. NaN bins).

    The record is an event-triggered logger, so thin coverage of the 5-minute grid is expected and is a WARN,
    not a FAIL. Only a record with no usable points at all stops the run."""
    ck = Checks()
    t0, t1 = pd.Timestamp(EVENT_START), pd.Timestamp(EVENT_END)
    grid = pd.date_range(t0, t1, freq=GRID_FREQ)
    info = {"xlsx": str(xlsx_path), "xlsx_md5": RS.md5_of(xlsx_path), "window_start": EVENT_START,
            "window_end": EVENT_END, "n_grid_points": int(len(grid)), "raw_rows_in_workbook": int(len(obs_df))}
    ck.add("NOTE", "workbook", "%s  (md5 %s, %d rows)" % (Path(xlsx_path).name, info["xlsx_md5"][:12], len(obs_df)))

    idx = obs_df.index
    inwin = np.asarray((idx >= t0) & (idx <= t1))
    raw_w = obs_df.loc[inwin]
    info["raw_rows_in_window"] = int(inwin.sum())
    real_grid = real_5min.reindex(grid)
    real_valid = real_grid.dropna()
    info["n_valid_points"] = int(len(real_valid))
    info["n_dropped_empty_bins"] = int(len(grid) - len(real_valid))

    if inwin.sum() == 0 or len(real_valid) == 0:
        ck.add("FAIL", "window coverage", "the workbook has no discharge records between %s and %s" % (EVENT_START, EVENT_END))
        return ck, info, real_valid, real_grid

    # duplicates and ordering
    dup_w = int(raw_w.index.duplicated().sum())
    dup_all = int(idx.duplicated().sum())
    info["duplicate_timestamps_in_window"] = dup_w
    info["duplicate_timestamps_whole_record"] = dup_all
    if dup_w:
        ck.add("WARN", "duplicate timestamps", "%d duplicated timestamp(s) inside the window; the scorer's .mean() "
               "silently averages them (so does the filled arm)" % dup_w)
    else:
        ck.add("PASS", "duplicate timestamps", "none inside the window (%d elsewhere in the record)" % dup_all)
    if not idx.is_monotonic_increasing:
        ck.add("WARN", "record order", "the workbook is not in time order; the scorer's resample copes with that, "
               "this script reproduces it")

    # values
    cfs = raw_w["cfs"]
    n_nan = int(cfs.isna().sum())
    n_neg = int((cfs < 0).sum())
    if n_nan or n_neg:
        ck.add("WARN", "discharge values", "%d missing and %d negative cfs value(s) inside the window" % (n_nan, n_neg))
    else:
        ck.add("PASS", "discharge values", "no missing or negative cfs inside the window")
    if cfs.notna().sum() == 0:
        ck.add("FAIL", "discharge values", "every cfs value inside the window is missing")
        return ck, info, real_valid, real_grid

    # the peak: a scale check against the stored simulations (catches a wrong unit or event), not a fit
    pk = float(cfs.max())
    pk_t = cfs.idxmax()
    pk_cms = pk * CFS_TO_CMS
    info["peak_cfs_raw"] = pk
    info["peak_m3s_raw"] = pk_cms
    info["peak_time_raw"] = str(pk_t)
    info["peak_m3s_on_grid"] = float(real_valid.max())
    if sim_peak_med is not None and np.isfinite(sim_peak_med) and sim_peak_med > 0:
        ratio = pk_cms / sim_peak_med
        info["stored_sim_median_peak_m3s"] = float(sim_peak_med)
        info["peak_scale_ratio"] = float(ratio)
        if SCALE_BAND[0] <= ratio <= SCALE_BAND[1]:
            ck.add("PASS", "peak (scale check)", "%.0f cfs = %.1f m3/s at %s; the stored runs' median simulated peak is "
                   "%.1f m3/s, so the real peak is the same order of magnitude (a units and event check, not a fit)"
                   % (pk, pk_cms, pk_t, sim_peak_med))
        else:
            ck.add("WARN", "peak (scale check)", "%.0f cfs = %.1f m3/s at %s is %.2g times the stored runs' median "
                   "simulated peak (%.1f m3/s); outside %.1f-%.1f times. Check the unit, the sheet and the event date"
                   % (pk, pk_cms, pk_t, ratio, sim_peak_med, SCALE_BAND[0], SCALE_BAND[1]))
    else:
        ck.add("NOTE", "peak (scale check)", "%.0f cfs = %.1f m3/s at %s (no stored simulations to compare the scale with)"
               % (pk, pk_cms, pk_t))

    # cadence, silences, coverage
    steps = pd.Series(raw_w.index.sort_values()).diff().dropna()
    steps = steps[steps > pd.Timedelta(0)]
    step_min = (steps.dt.total_seconds() / 60.0).to_numpy(dtype=float) if len(steps) else np.array([])
    cadence = float(np.median(step_min)) if len(step_min) else float("nan")
    info["raw_cadence_min"] = cadence
    info["raw_steps_over_15min"] = int((step_min > 15.0).sum()) if len(step_min) else 0
    info["raw_longest_step_min"] = float(step_min.max()) if len(step_min) else float("nan")
    info["raw_silences_over_3h"] = int((step_min > 180.0).sum()) if len(step_min) else 0
    ck.add("NOTE", "raw record", "%d records in the window; median step %.1f min, %d step(s) over 15 min, longest %.0f min "
           "(%d silence(s) over 3 h)" % (int(inwin.sum()), cadence, info["raw_steps_over_15min"],
                                         info["raw_longest_step_min"], info["raw_silences_over_3h"]))
    ck.add("NOTE", "as-is 5-minute points", "%d of the %d grid points in the window hold a record after the scorer's "
           "resampling; the scorer's .dropna() removes the other %d" % (len(real_valid), len(grid), len(grid) - len(real_valid)))
    sparse = len(real_valid) < 0.5 * len(grid)
    if sparse:
        ck.add("WARN", "event-triggered record", "the logger writes on change, so the as-is series keeps only %d of %d grid "
               "points and is blind to simulated flow in the other %d bins (false early flow, dry hours). The filled arm "
               "scores all %d points; see --gauge_mode and --primary_arm"
               % (len(real_valid), len(grid), len(grid) - len(real_valid), len(grid)))
    if len(real_valid) < MIN_ASIS_POINTS:
        ck.add("FAIL", "as-is window coverage", "only %d scorable as-is point(s) (minimum %d): not enough to score"
               % (len(real_valid), MIN_ASIS_POINTS))
    elif len(real_valid) < min_points:
        ck.add("WARN", "as-is window coverage", "only %d scorable as-is points (--min_points is %d): thin for KGE and NSE; "
               "read the as-is arm with that in mind" % (len(real_valid), min_points))
    else:
        ck.add("PASS", "as-is window coverage", "%d scorable as-is points, first %s, last %s"
               % (len(real_valid), real_valid.index[0], real_valid.index[-1]))
    return ck, info, real_valid, real_grid


# ------------------------------------------------------------------
# The FILLED series: the scorer's own bin means plus the empty bins filled from the record
# ------------------------------------------------------------------
def build_filled(obs_df, real_5min, max_gap_min, max_hold_hr):
    """Returns (filled Series on the window grid, fill_type Series, info dict).

    fill_type: 'record'    the scorer's own 5-minute bin mean (bin holds at least one reading)
               'linear'    empty bin between two readings at most max_gap_min apart: straight line, evaluated at the
                           bin's midpoint (exact for the bin mean, since a bin with no reading has no kink inside)
               'hold'      empty bin inside a longer silence: the last reading is held, if it is at most max_hold_hr old
               'uncovered' no earlier reading, or the last one is older than max_hold_hr: left out (NaN)
    """
    t0, t1 = pd.Timestamp(EVENT_START), pd.Timestamp(EVENT_END)
    grid = pd.date_range(t0, t1, freq=GRID_FREQ)
    raw = obs_df["Observed_CMS"].dropna()
    raw = raw.groupby(level=0).mean().sort_index()           # duplicate stamps averaged, as the scorer's .mean() does
    epoch = pd.Timestamp("1970-01-01")
    sec = np.asarray((raw.index - epoch) / pd.Timedelta(seconds=1), dtype=float)
    val = raw.to_numpy(dtype=float)
    asis = real_5min.reindex(grid)
    asis_v = asis.to_numpy(dtype=float)
    gsec = np.asarray((grid - epoch) / pd.Timedelta(seconds=1), dtype=float)
    half = 0.5 * pd.Timedelta(GRID_FREQ).total_seconds()
    max_gap = float(max_gap_min) * 60.0
    max_hold = float(max_hold_hr) * 3600.0
    fv = asis_v.copy()
    ft = np.empty(len(grid), dtype=object)
    for k in range(len(grid)):
        if np.isfinite(asis_v[k]):
            ft[k] = "record"
            continue
        c = gsec[k] + half
        i = int(np.searchsorted(sec, c))
        if i == 0:
            ft[k] = "uncovered"
            continue
        t_prev, v_prev = sec[i - 1], val[i - 1]
        if i < len(sec) and (sec[i] - t_prev) <= max_gap:
            w = (c - t_prev) / (sec[i] - t_prev)
            fv[k] = v_prev + w * (val[i] - v_prev)
            ft[k] = "linear"
        elif (c - t_prev) <= max_hold:
            fv[k] = v_prev
            ft[k] = "hold"
        else:
            ft[k] = "uncovered"
    filled = pd.Series(fv, index=grid, name="Observed_filled")
    ftype = pd.Series(ft, index=grid, name="fill_type")
    info = {"max_gap_min": float(max_gap_min), "max_hold_hr": float(max_hold_hr)}
    for kind in ("record", "linear", "hold", "uncovered"):
        info["n_" + kind] = int((ftype == kind).sum())
    pre = raw[raw.index < t0]
    if len(pre):
        info["last_record_before_window"] = {"time": str(pre.index[-1]), "cfs": float(pre.iloc[-1] / CFS_TO_CMS),
                                              "hours_before_window": float((t0 - pre.index[-1]) / pd.Timedelta(hours=1))}
    else:
        info["last_record_before_window"] = None
    return filled, ftype, info


def check_filled(filled, ftype, finfo, real_valid, raw_w, min_points):
    """GATE 1b. Returns (Checks, info)."""
    ck = Checks()
    info = {}
    n = len(filled)
    ck.add("NOTE", "filled arm bins", "%d bins: %d record (the scorer's own values), %d line (gap <= %g min), %d hold "
           "(silence > %g min, last reading <= %g h old), %d uncovered"
           % (n, finfo["n_record"], finfo["n_linear"], finfo["max_gap_min"], finfo["n_hold"], finfo["max_gap_min"],
              finfo["max_hold_hr"], finfo["n_uncovered"]))
    pre = finfo.get("last_record_before_window")
    if finfo["n_hold"] and pre is not None:
        sev = "WARN" if pre["cfs"] > 1.0 else "NOTE"
        ck.add(sev, "filled arm start", "the last reading before the window is %.2f cfs at %s (%.1f h before 16:00); "
               "it is held until the first reading inside the window%s"
               % (pre["cfs"], pre["time"], pre["hours_before_window"],
                  " -- that is not a dry start; check it" if sev == "WARN" else ""))
    if finfo["n_uncovered"]:
        first_unc = ftype[ftype == "uncovered"].index[0]
        hint = ""
        if pre is not None:
            hint = ("  The last reading before the window is at %s, %.1f h before 16:00; --max_hold_hr above %.1f would use it"
                    % (pre["time"], pre["hours_before_window"], pre["hours_before_window"]))
        else:
            hint = "  There is no reading at all before the window."
        ck.add("WARN", "filled arm uncovered", "%d bin(s) have no reading within --max_hold_hr %g h (first: %s); they are left "
               "out of the filled arm (it scores %d points, not %d).%s"
               % (finfo["n_uncovered"], finfo["max_hold_hr"], first_unc, int(filled.notna().sum()), n, hint))
    n_f = int(filled.notna().sum())
    if n_f < min_points:
        ck.add("FAIL", "filled arm coverage", "only %d scorable filled points (minimum %d, --min_points)" % (n_f, min_points))
    else:
        ck.add("PASS", "filled arm coverage", "%d scorable filled points of %d" % (n_f, n))
    # equal to the as-is series wherever the scorer has a value
    d = float((filled.reindex(real_valid.index) - real_valid).abs().max()) if len(real_valid) else 0.0
    info["max_abs_diff_vs_asis"] = d
    if d == 0.0:
        ck.add("PASS", "filled = as-is at records", "identical at all %d bins that hold a record" % len(real_valid))
    else:
        ck.add("FAIL", "filled = as-is at records", "the filled series differs from the as-is series at record bins "
               "(worst %.2e m3/s)" % d)
    # no negative fill
    fmin = float(np.nanmin(filled.to_numpy(dtype=float))) if n_f else float("nan")
    if np.isfinite(fmin) and fmin < 0:
        ck.add("WARN", "filled arm sign", "the filled series has negative discharge (minimum %.3f m3/s)" % fmin)
    # volume against the raw records' own trapezoid sum (independent of the binning)
    step_s = pd.Timedelta(GRID_FREQ).total_seconds()
    vol_f = float(filled.dropna().sum() * step_s)
    vol_a = float(real_valid.sum() * step_s)
    r = raw_w["Observed_CMS"].dropna().groupby(level=0).mean().sort_index()
    if len(r) >= 2:
        tt = np.asarray((r.index - r.index[0]) / pd.Timedelta(seconds=1), dtype=float)
        vv = r.to_numpy(dtype=float)
        vol_raw = float(np.sum(0.5 * (vv[1:] + vv[:-1]) * np.diff(tt)))
    else:
        vol_raw = float("nan")
    info.update({"volume_filled_m3": vol_f, "volume_asis_bins_m3": vol_a, "volume_raw_trapezoid_m3": vol_raw})
    if np.isfinite(vol_raw) and vol_raw > 0:
        rel = abs(vol_f - vol_raw) / vol_raw
        info["volume_rel_diff"] = float(rel)
        txt = ("filled volume %.0f m3 (%.1f acre-ft) vs the raw records' trapezoid sum %.0f m3: differ by %.1f%%; the %d "
               "as-is bins alone hold %.0f%% of the filled volume"
               % (vol_f, vol_f / 1233.48, vol_raw, 100 * rel, len(real_valid), 100 * vol_a / vol_f if vol_f > 0 else float("nan")))
        if rel > VOL_WARN_FRAC:
            txt += ("  [%d bin(s) were HELD across silences longer than --max_gap_min %g; a hold across a long gap with very "
                    "different readings at its ends changes the volume. Raise --max_gap_min or look at the observed figure]"
                    % (finfo["n_hold"], finfo["max_gap_min"]))
        ck.add("FAIL" if rel > VOL_FAIL_FRAC else ("WARN" if rel > VOL_WARN_FRAC else "PASS"), "filled arm volume", txt)
    else:
        ck.add("NOTE", "filled arm volume", "filled volume %.0f m3; the raw records give no positive volume to compare with" % vol_f)
    return ck, info


def blind_share(cache, run_ids, series_of, ref_idx, asis_idx, filled):
    """Median share of each run's simulated volume in bins the as-is arm cannot see, split by where they are.
    Returns (DataFrame, flood_start, flood_end). Pure function of the stored simulations and the observed series."""
    pos = filled[filled > 0]
    if len(pos) == 0:
        return None, None, None
    f0, f1 = pos.index[0], pos.index[-1]
    sees = np.asarray(ref_idx.isin(asis_idx))
    before = np.asarray(ref_idx < f0)
    after = np.asarray(ref_idx > f1)
    masks = {"seen_by_asis": sees,
             "blind_before_flood": (~sees) & before,
             "blind_inside_flood": (~sees) & (~before) & (~after),
             "blind_after_flood": (~sees) & after}
    rows = []
    for series in (SER_ON, SER_OFF):
        acc = {k: [] for k in masks}
        for rid in run_ids:
            if series_of.get(rid) != series:
                continue
            v = cache[rid]["Simulated"].reindex(ref_idx).to_numpy(dtype=float)
            v = np.where(np.isfinite(v), v, 0.0)
            tot = float(v.sum())
            if tot <= 0:
                continue
            for k, m in masks.items():
                acc[k].append(float(v[m].sum()) / tot)
        if not acc["seen_by_asis"]:
            continue
        row = {"series": series, "n_runs": len(acc["seen_by_asis"])}
        for k in masks:
            row["median_share_" + k] = float(np.median(acc[k]))
        blind = 1.0 - np.asarray(acc["seen_by_asis"])
        row["median_share_blind_total"] = float(np.median(blind))
        row["max_share_blind_total"] = float(blind.max())
        rows.append(row)
    return pd.DataFrame(rows), f0, f1


# ------------------------------------------------------------------
# GATE 4: the grid
# ------------------------------------------------------------------
def check_grid(cache, run_ids, real_valid):
    """Returns (Checks, ref_idx, scored_idx, bad dict). `bad` holds runs whose stored index differs from the majority."""
    ck = Checks()
    t0, t1 = pd.Timestamp(EVENT_START), pd.Timestamp(EVENT_END)
    sigs = {}
    for rid in run_ids:
        ix = cache[rid].index
        sigs.setdefault((len(ix), str(ix[0]), str(ix[-1])), []).append(rid)
    ref_key = max(sigs, key=lambda k: len(sigs[k]))
    ref_idx = cache[sigs[ref_key][0]].index
    bad = {}
    for rid in run_ids:
        if not cache[rid].index.equals(ref_idx):
            bad[rid] = "its stored 5-minute grid differs from the one shared by most runs"
    if bad:
        ck.add("FAIL", "GATE 4 same grid", "%d of %d run(s) have a different stored grid, e.g. %s"
               % (len(bad), len(run_ids), ", ".join(list(bad)[:3])))
    else:
        ck.add("PASS", "GATE 4 same grid", "all %d stored Simulated series share one grid of %d points (%s to %s)"
               % (len(run_ids), len(ref_idx), ref_idx[0], ref_idx[-1]))
    spacing = pd.Series(ref_idx).diff().dropna().unique()
    if len(spacing) == 1 and spacing[0] == pd.Timedelta(minutes=5):
        ck.add("PASS", "GATE 4 regular 5 minutes", "every step is exactly 5 minutes")
    else:
        ck.add("FAIL", "GATE 4 regular 5 minutes", "steps found: %s" % ", ".join(str(s) for s in spacing[:4]))
    if ref_idx[0] == t0 and ref_idx[-1] == t1:
        ck.add("PASS", "GATE 4 spans the window", "the stored grid runs from %s to %s" % (t0, t1))
    else:
        ck.add("WARN", "GATE 4 spans the window", "the stored grid runs %s to %s, not the full window %s to %s"
               % (ref_idx[0], ref_idx[-1], t0, t1))
    missing = real_valid.index.difference(ref_idx)
    if len(missing):
        ck.add("FAIL", "GATE 4 real bins on grid", "%d real 5-minute bin(s) are not on the stored grid, e.g. %s"
               % (len(missing), ", ".join(str(m) for m in missing[:3])))
    else:
        ck.add("PASS", "GATE 4 real bins on grid", "all %d real bins in the window are on the stored grid"
               % len(real_valid))
    scored_idx = real_valid.index.intersection(ref_idx)
    return ck, ref_idx, scored_idx, bad


# ------------------------------------------------------------------
# GATE 3: the actual scorer, in real-gauge mode, on one stored parameter point
# ------------------------------------------------------------------
def verify_gate3(compare_path, metrics_path, all_res, m_all, cache, real_scored, scored_idx):
    """Compare ONE real-gauge-mode scorer run with this script. Prints deviations only, never a fit value."""
    ck = Checks()
    info = {"compare_csv": str(compare_path), "metrics_csv": str(metrics_path) if metrics_path else None}
    try:
        df = RS.load_compare(compare_path)
    except Exception as e:
        ck.add("FAIL", "GATE 3 read", "cannot read %s: %s" % (compare_path, str(e)[:100]))
        return ck, info
    # 3a: observed
    if not df.index.equals(scored_idx):
        both = df.index.intersection(scored_idx)
        ck.add("FAIL", "GATE 3a real Observed", "the scorer's gauge-mode grid has %d points, this script's has %d "
               "(%d in common). The two read the real record differently" % (len(df), len(scored_idx), len(both)))
        return ck, info
    ref = real_scored.reindex(df.index).to_numpy(dtype=float)
    got = df["Observed"].to_numpy(dtype=float)
    d_obs = float(np.max(np.abs(ref - got)))
    info["obs_max_abs_diff_m3s"] = d_obs
    if np.allclose(ref, got, rtol=RS.OBS_RTOL, atol=RS.OBS_ATOL):
        ck.add("PASS", "GATE 3a real Observed", "same grid (%d points) and same values as the scorer's gauge mode "
               "(worst abs diff %.2e m3/s)" % (len(df), d_obs))
    else:
        ck.add("FAIL", "GATE 3a real Observed", "same grid but values differ (worst abs diff %.2e m3/s)" % d_obs)
        return ck, info

    # identify the stored run at the same parameter point
    pr = None
    if metrics_path is not None:
        try:
            mdf = pd.read_csv(metrics_path)
            pr = mdf.iloc[0]
        except Exception as e:
            ck.add("FAIL", "GATE 3 read", "cannot read %s: %s" % (metrics_path, str(e)[:100]))
            return ck, info
        mode = str(pr.get("obs_mode", ""))
        info["obs_mode"] = mode
        if mode != "gauge":
            ck.add("FAIL", "GATE 3 mode", "the scorer's metrics file says obs_mode = %r, not 'gauge': that run was scored "
                   "against a synthetic truth, so it cannot confirm the real-gauge path" % mode)
            return ck, info
        ck.add("PASS", "GATE 3 mode", "the scorer's metrics file says obs_mode = gauge")
    if pr is None:
        ck.add("WARN", "GATE 3c metrics", "no --gate3_metrics given: only the observed series and the simulation "
               "could be checked, not the metrics")
        return ck, info
    try:
        target = np.array([float(pr["Ks_mult"]), float(pr["f_RS_abs"]), float(pr["channelconductivity_mmhr"]),
                           float(pr["optpercolation"])])
    except Exception:
        ck.add("FAIL", "GATE 3 parameter point", "the metrics file lacks Ks_mult, f_RS_abs, channelconductivity_mmhr or "
               "optpercolation")
        return ck, info
    cols = PCOLS + ["optpercolation"]
    A = all_res[cols].to_numpy(float)
    hit = np.all(np.isclose(A, target[None, :], rtol=1e-4, atol=0.0), axis=1)
    if hit.sum() != 1:
        ck.add("FAIL", "GATE 3 parameter point", "the gauge-mode run's (Ks, f, cc, optpercolation) = (%g, %g, %g, %g) "
               "matches %d stored runs (need exactly 1). Re-run one of the 500 stored points"
               % (target[0], target[1], target[2], target[3], int(hit.sum())))
        return ck, info
    rid = all_res.loc[hit, "run_id"].iloc[0]
    info["matched_stored_run"] = rid
    ck.add("PASS", "GATE 3 parameter point", "matches stored run %s" % rid)

    # 3b: simulation
    stored = cache[rid]["Simulated"].reindex(df.index).to_numpy(dtype=float)
    sim = df["Simulated"].to_numpy(dtype=float)
    peak = float(np.nanmax(np.abs(stored))) or 1.0
    d_sim = float(np.nanmax(np.abs(stored - sim)))
    info["sim_max_abs_diff_m3s"] = d_sim
    info["sim_rel_to_peak"] = d_sim / peak
    if d_sim / peak <= GATE3_SIM_REL:
        ck.add("PASS", "GATE 3b re-run Simulated", "the re-run reproduces the stored simulation (largest difference "
               "%.2e m3/s, %.1e of its peak)" % (d_sim, d_sim / peak))
    else:
        ck.add("FAIL", "GATE 3b re-run Simulated", "the re-run differs from the stored simulation by %.2e m3/s "
               "(%.1e of its peak): today's build or inputs are not the ones that made the 500 runs"
               % (d_sim, d_sim / peak))
        return ck, info

    # 3c: metrics
    mine = m_all.set_index("run_id").loc[rid]
    names = RS.SCALAR_METRICS + RS.PHASE_METRICS
    off, worst, nchk = [], 0.0, 0
    for c in names:
        if c not in pr.index:
            continue
        try:
            stored_v = float(pr[c])
        except (TypeError, ValueError):
            continue
        new_v = float(mine[c])
        nchk += 1
        if not RS.close(new_v, stored_v, RS.RTOL, RS.ATOL):
            off.append("%s (scorer %.6g, here %.6g)" % (c, stored_v, new_v))
        elif np.isfinite(stored_v) and np.isfinite(new_v):
            worst = max(worst, abs(new_v - stored_v))
    info["metrics_checked"] = nchk
    info["metrics_worst_abs_dev"] = worst
    if off:
        ck.add("FAIL", "GATE 3c metrics", "%d of %d metric(s) differ: %s" % (len(off), nchk, "; ".join(off[:3])))
    else:
        ck.add("PASS", "GATE 3c metrics", "all %d metrics match the scorer's real-gauge run (worst abs deviation "
               "%.2e)" % (nchk, worst))
    return ck, info


# ------------------------------------------------------------------
# Analysis (pure functions of the two metric tables)
# ------------------------------------------------------------------
def top_rows(m, n, col="kge_2012"):
    ok = m[np.isfinite(m[col])]
    return ok.sort_values(col, ascending=False, kind="mergesort").head(n)


def series_summary(m, series, top_n):
    ok = m[np.isfinite(m["kge_2012"])]
    best = ok.loc[ok["kge_2012"].idxmax()]
    top = top_rows(m, top_n)
    feas = ok[np.abs(ok["pbias_pct"]) < RS.FEASIBLE_PBIAS_PCT]
    row = {"series": series, "n_scored": int(len(m)), "n_valid_kge": int(len(ok)),
           "best_kge_2012": float(best["kge_2012"]), "best_run_id": best["run_id"],
           "best_Ks": float(best["Ks_mult"]), "best_f": float(best["f_RS_abs"]),
           "best_cc": float(best["channelconductivity_mmhr"]) if series == SER_ON else np.nan,
           "best_pbias_pct": float(best["pbias_pct"]), "best_kge_r": float(best["kge_r"]),
           "best_kge_beta": float(best["kge_beta"]), "best_kge_gamma": float(best["kge_gamma"]),
           "best_nse": float(best["nse"]), "best_rmse_m3s": float(best["rmse_m3s"]),
           "best_peak_error_pct": float(best["peak_error_pct"]),
           "best_peak_timing_error_hr": float(best["peak_timing_error_hr"]),
           "top_n": int(len(top)), "top_n_mean_kge": float(top["kge_2012"].mean()),
           "top_n_median_kge": float(top["kge_2012"].median()),
           "median_kge_all_runs": float(ok["kge_2012"].median()),
           "n_feasible": int(len(feas)), "pct_feasible": 100.0 * len(feas) / len(ok)}
    return row


def cc_location(m_on, n, floor_tol):
    """Where do the n best cc-ON runs (by KGE_2012) sit in the cc box? Chance puts 1/3 in each third."""
    top = top_rows(m_on, n)
    cc = top["channelconductivity_mmhr"].to_numpy(float)
    u = cc_unit(cc)
    k = len(cc)
    lo_cc = RS.BOX["channelconductivity_mmhr"][0]
    hi_cc = RS.BOX["channelconductivity_mmhr"][1]
    ratio = (lo_cc + floor_tol) / lo_cc
    n_low = int((u < 1.0 / 3.0).sum())
    n_high = int((u >= 2.0 / 3.0).sum())
    ok_all = m_on[np.isfinite(m_on["kge_2012"])]
    u_all = cc_unit(ok_all["channelconductivity_mmhr"].to_numpy(float))
    lo_log, hi_log = cc_box_logs()
    third_edge = 10 ** (lo_log + (hi_log - lo_log) / 3.0)
    out = {"top_n": k, "n_lowest_third": n_low, "n_middle_third": k - n_low - n_high, "n_highest_third": n_high,
           "n_at_floor": int((cc <= lo_cc + floor_tol).sum()), "floor_tol_mmhr": floor_tol,
           "n_at_ceiling": int((cc >= hi_cc / ratio).sum()),
           "lowest_third_upper_edge_mmhr": float(third_edge),
           "cc_min": float(cc.min()), "cc_median": float(np.median(cc)), "cc_max": float(cc.max()),
           "span_frac_of_box": float((np.log10(cc.max()) - np.log10(cc.min())) / (hi_log - lo_log)),
           "p_lowest_third": float(stats.binom.sf(n_low - 1, k, 1.0 / 3.0)),
           "p_highest_third": float(stats.binom.sf(n_high - 1, k, 1.0 / 3.0)),
           "share_all_runs_lowest_third": float((u_all < 1.0 / 3.0).mean()),
           "share_all_runs_highest_third": float((u_all >= 2.0 / 3.0).mean()),
           "best_run_cc": float(cc[0])}
    return out


def pair_runs(m_on, m_off, rtol=1e-4):
    """Pair each ON run with the OFF run at the same (Ks, f, cc) point. Returns (DataFrame or None, note)."""
    A = m_on[PCOLS].to_numpy(float)
    B = m_off[PCOLS].to_numpy(float)
    close = np.all(np.isclose(A[:, None, :], B[None, :, :], rtol=rtol, atol=0.0), axis=2)
    n_match = close.sum(axis=1)
    taken = close.sum(axis=0)
    use = (n_match == 1)
    j = close.argmax(axis=1)
    use &= (taken[j] == 1)
    if use.sum() == 0:
        return None, "no ON run has exactly one OFF partner at the same (Ks, f, cc)"
    on = m_on.reset_index(drop=True)[use].reset_index(drop=True)
    off = m_off.reset_index(drop=True).iloc[j[use]].reset_index(drop=True)
    p = pd.DataFrame({"Ks_mult": on["Ks_mult"], "f_RS_abs": on["f_RS_abs"],
                      "channelconductivity_mmhr": on["channelconductivity_mmhr"],
                      "log10_cc": np.log10(on["channelconductivity_mmhr"].to_numpy(float)),
                      "run_id_on": on["run_id"], "run_id_off": off["run_id"],
                      "kge_on": on["kge_2012"], "kge_off": off["kge_2012"],
                      "pbias_on": on["pbias_pct"], "pbias_off": off["pbias_pct"]})
    p["diff"] = p["kge_on"] - p["kge_off"]
    p = p[np.isfinite(p["diff"])].reset_index(drop=True)
    note = "%d of %d ON runs paired with one OFF run" % (len(p), len(m_on))
    return p, note


def mark_good_region(paired, good_frac):
    """The top good_frac of pairs by the cc-OFF fit. Chosen from the control alone (cc is inert there)."""
    n_good = max(MIN_GOOD, int(round(good_frac * len(paired))))
    n_good = min(n_good, len(paired))
    order = paired["kge_off"].sort_values(ascending=False, kind="mergesort").index[:n_good]
    flag = np.zeros(len(paired), dtype=bool)
    flag[paired.index.get_indexer(order)] = True
    return flag


def _slope_boot(x, y, B, seed):
    n = len(x)
    rng = np.random.default_rng([seed, 9])
    ix = rng.integers(0, n, size=(B, n))
    xs, ys = x[ix], y[ix]
    xm = xs.mean(axis=1, keepdims=True)
    ym = ys.mean(axis=1, keepdims=True)
    with np.errstate(all="ignore"):
        sl = ((xs - xm) * (ys - ym)).sum(axis=1) / ((xs - xm) ** 2).sum(axis=1)
    return sl[np.isfinite(sl)]


def paired_stats(paired, good_flag, good_frac, boot, seed):
    rows = []
    regions = [("all pairs", np.ones(len(paired), dtype=bool)),
               ("good region (top %g%% by cc-OFF fit)" % (100.0 * good_frac), good_flag)]
    for name, sel in regions:
        sub = paired[sel]
        n = len(sub)
        row = {"region": name, "n_pairs": int(n), "median_diff": float(sub["diff"].median()),
               "mean_diff": float(sub["diff"].mean()), "share_on_better": float((sub["diff"] > 0).mean()),
               "median_kge_on": float(sub["kge_on"].median()), "median_kge_off": float(sub["kge_off"].median())}
        x = sub["log10_cc"].to_numpy(float)
        y = sub["diff"].to_numpy(float)
        if n >= 5 and np.ptp(x) > 0 and np.ptp(y) > 0:
            r = stats.spearmanr(x, y)
            row["spearman_rho_diff_vs_log_cc"] = float(r[0])
            row["spearman_p"] = float(r[1])
            slope = float(np.polyfit(x, y, 1)[0])
            sb = _slope_boot(x, y, boot, seed)
            row["slope_per_decade_cc"] = slope
            row["slope_ci_lo"] = float(np.percentile(sb, 2.5)) if len(sb) else np.nan
            row["slope_ci_hi"] = float(np.percentile(sb, 97.5)) if len(sb) else np.nan
        u = cc_unit(sub["channelconductivity_mmhr"].to_numpy(float))
        for lab, lo_u, hi_u in (("low", -1.0, 1.0 / 3.0), ("mid", 1.0 / 3.0, 2.0 / 3.0), ("high", 2.0 / 3.0, 2.0)):
            s2 = sub[(u >= lo_u) & (u < hi_u)]
            row["n_%s_third" % lab] = int(len(s2))
            row["median_diff_%s_third" % lab] = float(s2["diff"].median()) if len(s2) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def decision_row(gap, margin, loc):
    """Which row of handoff section 4 do the numbers fit? Returns an integer key of READINGS."""
    low = loc["p_lowest_third"] < ALPHA
    high = loc["p_highest_third"] < ALPHA
    if gap < -margin:
        return 4
    if abs(gap) <= margin:
        return 1 if low else 5
    if low or high or loc["span_frac_of_box"] > STABLE_SPAN:
        return 3
    return 2


# ------------------------------------------------------------------
# Printing
# ------------------------------------------------------------------
def print_gates(ck_rows, gate3_state):
    print("\n" + "-" * 110)
    print("GATES")
    print("-" * 110)
    for r in ck_rows:
        print("  %-5s %-34s %s" % (r["status"], r["check"], r["detail"]))
    print("  GATE 3 (gauge-mode fidelity): %s" % gate3_state)


def print_summary(sdf, margin, n_scored, arm):
    on = sdf[sdf.series == SER_ON].iloc[0]
    off = sdf[sdf.series == SER_OFF].iloc[0]
    print("\n" + "=" * 110)
    print("RE-SCORED AGAINST THE REAL GAUGE, %s ARM (%s to %s, %d points)" % (arm.upper(), EVENT_START, EVENT_END, n_scored))
    print("=" * 110)
    print("%-34s %14s %14s %12s" % ("", "cc ON", "cc OFF", "ON - OFF"))

    def line(label, key, spec):
        a, b = on[key], off[key]
        d = a - b if (np.isfinite(a) and np.isfinite(b)) else np.nan
        print("%-34s %14s %14s %12s" % (label, fnum(a, spec), fnum(b, spec), fnum(d, "+" + spec)))

    line("best KGE_2012", "best_kge_2012", ".4f")
    line("mean KGE_2012 of the best %d" % int(on["top_n"]), "top_n_mean_kge", ".4f")
    line("median KGE_2012, all runs", "median_kge_all_runs", ".4f")
    line("runs with |PBIAS| < %g%%" % RS.FEASIBLE_PBIAS_PCT, "n_feasible", ".0f")
    print("\nBest run in each series:")
    print("  %-8s %-26s %6s %8s %8s | %7s %6s %6s %6s %8s %8s"
          % ("series", "run_id", "Ks", "f", "cc", "PBIAS%", "r", "beta", "gamma", "peak err%", "peak hr"))
    for row in (on, off):
        print("  %-8s %-26s %6.2f %8.4f %8s | %7.1f %6.3f %6.3f %6.3f %8.1f %8.2f"
              % (row["series"].replace("cc ", ""), row["best_run_id"], row["best_Ks"], row["best_f"],
                 fcc(row["best_cc"]) if row["series"] == SER_ON else "inert", row["best_pbias_pct"],
                 row["best_kge_r"], row["best_kge_beta"], row["best_kge_gamma"], row["best_peak_error_pct"],
                 row["best_peak_timing_error_hr"]))
    print("  (peak hr = simulated peak time minus real peak time; positive = the model peaks late.)")


def print_cc_location(loc_rows, main_n, arm):
    print("\n" + "=" * 110)
    print("WHERE DO THE BEST cc-ON RUNS PUT cc, %s ARM?   (chance puts one third in each third of the log cc range)" % arm.upper())
    print("=" * 110)
    print("  %5s | %6s %6s %6s | %8s %8s | %-24s %5s | %9s %9s"
          % ("best", "lowest", "middle", "highest", "at floor", "at ceil.", "cc min / median / max", "span", "p lowest", "p highest"))
    for r in loc_rows:
        print("  %5d | %6d %6d %6d | %8d %8d | %-24s %5.2f | %9s %9s%s"
              % (r["top_n"], r["n_lowest_third"], r["n_middle_third"], r["n_highest_third"], r["n_at_floor"],
                 r["n_at_ceiling"], "%s / %s / %s" % (fcc(r["cc_min"]), fcc(r["cc_median"]), fcc(r["cc_max"])),
                 r["span_frac_of_box"], fmt_p(r["p_lowest_third"]), fmt_p(r["p_highest_third"]),
                 "   <- used for the reading below" if r["top_n"] == main_n else ""))
    r0 = loc_rows[0]
    print("  lowest third = cc up to %s mm/hr; the 'at floor' count is within %g mm/hr of %g mm/hr. Share of ALL runs in the "
          "lowest third: %.2f (chance 0.33); highest third: %.2f."
          % (fcc(r0["lowest_third_upper_edge_mmhr"]), r0["floor_tol_mmhr"], RS.BOX["channelconductivity_mmhr"][0],
             r0["share_all_runs_lowest_third"], r0["share_all_runs_highest_third"]))
    print("  span = how much of the log cc box the best runs cover (1 = all of it; small = they agree). p = chance of at least "
          "that many in that third if cc were irrelevant (binomial; the best runs are not independent, so a guide only).")


def print_paired(pst, note, arm):
    print("\n" + "=" * 110)
    print("PAIRED DIFFERENCE, %s ARM, ON minus OFF at the same (Ks, f, cc) point, against log cc   (%s)" % (arm.upper(), note))
    print("=" * 110)
    print("  %-40s %5s | %8s %8s %8s | %-26s | %-24s"
          % ("region", "n", "median", "mean", "ON>OFF", "rank corr with log cc [p]", "slope per decade of cc [95%]"))
    for _, r in pst.iterrows():
        rho = ("%+.2f [p %s]" % (r["spearman_rho_diff_vs_log_cc"], fmt_p(r["spearman_p"]))
               if "spearman_rho_diff_vs_log_cc" in r and np.isfinite(r.get("spearman_rho_diff_vs_log_cc", np.nan)) else "n/a")
        sl = ("%+.3f [%+.3f, %+.3f]" % (r["slope_per_decade_cc"], r["slope_ci_lo"], r["slope_ci_hi"])
              if "slope_per_decade_cc" in r and np.isfinite(r.get("slope_per_decade_cc", np.nan)) else "n/a")
        print("  %-40s %5d | %+8.4f %+8.4f %7.0f%% | %-26s | %-24s"
              % (r["region"], r["n_pairs"], r["median_diff"], r["mean_diff"], 100 * r["share_on_better"], rho, sl))
    print("  Median difference by third of the cc range (low / middle / high):")
    for _, r in pst.iterrows():
        print("    %-40s n = %d / %d / %d   median diff = %s / %s / %s"
              % (r["region"], r["n_low_third"], r["n_mid_third"], r["n_high_third"],
                 fnum(r["median_diff_low_third"], "+.4f"), fnum(r["median_diff_mid_third"], "+.4f"),
                 fnum(r["median_diff_high_third"], "+.4f")))
    print("  Positive = ON fits better. A slope near 0 = fit is flat as cc rises; negative = worse as cc rises. The OFF run does not "
          "depend on cc, so the trend comes from the ON runs.")


def print_reading(gap, gap_mean, margin, loc, row_best, row_mean, arm):
    print("\n" + "=" * 110)
    print("READING AGAINST HANDOFF SECTION 4, %s ARM   (your margin: %g KGE_2012;  mine: ALPHA %g, STABLE_SPAN %g, floor within %g mm/hr)"
          % (arm.upper(), margin, ALPHA, STABLE_SPAN, FLOOR_TOL_MMHR))
    print("=" * 110)
    print("  Best ON minus best OFF            : %+.4f   -> row %d" % (gap, row_best))
    print("  Mean of best %d ON minus OFF       : %+.4f   -> row %d" % (loc["top_n"], gap_mean, row_mean))
    print("  Best-run cc (n = %d): lowest third p = %s, highest third p = %s, span %.2f of the box"
          % (loc["top_n"], fmt_p(loc["p_lowest_third"]), fmt_p(loc["p_highest_third"]), loc["span_frac_of_box"]))
    print("\n  Row %d: %s" % (row_best, READINGS[row_best]))
    if row_mean != row_best:
        print("\n  NOTE: the top-%d mean points to row %d, not row %d. The two readings disagree, so treat the answer as "
              "unsettled:\n        row %d: %s" % (loc["top_n"], row_mean, row_best, row_mean, READINGS[row_mean]))


def print_blind_share(bs, f0, f1, n_asis, n_grid):
    print("\n" + "-" * 110)
    print("BLIND SHARE: how much of each run's SIMULATED volume sits in bins the as-is arm cannot see   (no fit involved)")
    print("-" * 110)
    if bs is None or not len(bs):
        print("  not available (no flow above zero in the filled series, or no simulated volume)")
        return
    print("  The as-is arm keeps %d of %d bins. 'Flood' = %s to %s (first to last bin above zero in the filled series)."
          % (n_asis, n_grid, f0, f1))
    print("  %-8s %5s | %10s | %12s %12s %12s | %12s %12s" % ("series", "runs", "seen as-is", "blind before", "blind inside",
                                                           "blind after", "blind total", "worst run"))
    for _, r in bs.iterrows():
        print("  %-8s %5d | %9.0f%% | %11.0f%% %11.0f%% %11.0f%% | %11.0f%% %11.0f%%"
              % (r["series"].replace("cc ", ""), r["n_runs"], 100 * r["median_share_seen_by_asis"],
                 100 * r["median_share_blind_before_flood"], 100 * r["median_share_blind_inside_flood"],
                 100 * r["median_share_blind_after_flood"], 100 * r["median_share_blind_total"],
                 100 * r["max_share_blind_total"]))
    worst = float(bs["median_share_blind_total"].max())
    if worst >= BLIND_WARN_SHARE:
        print("  NOTE: for the median run %.0f%% of the simulated volume is invisible to the as-is arm. Whatever the model does "
              "there (false early flow before the flood, flow where the record has holes) is never scored as-is." % (100 * worst))
    print("  (medians over runs; 'worst run' = the largest blind share of any single run.)")


def print_arm_agreement(ag, primary, margin):
    print("\n" + "=" * 110)
    print("ARM AGREEMENT   (primary arm, fixed before the run: %s;  margin %g KGE_2012)" % (primary.upper(), margin))
    print("=" * 110)
    print("  %-7s %7s | %9s %9s %9s | %9s | %8s %8s" % ("arm", "points", "best ON", "best OFF", "ON - OFF", "top-N gap",
                                                       "row best", "row mean"))
    for _, r in ag.iterrows():
        print("  %-7s %7d | %9.4f %9.4f %+9.4f | %+9.4f | %8d %8d%s"
              % (r["arm"], r["n_points"], r["best_kge_on"], r["best_kge_off"], r["gap_best"], r["gap_top_n_mean"],
                 r["row_best"], r["row_mean"], "   <- primary" if r["primary"] else ""))
    rows = set(int(x) for x in ag["row_best"])
    if len(ag) == 1:
        return
    if len(rows) == 1:
        print("\n  Both arms point to row %d of the handoff table." % rows.pop())
    else:
        prim = ag[ag["primary"]].iloc[0]
        oth = ag[~ag["primary"]].iloc[0]
        print("\n  THE ARMS DISAGREE: %s (primary) points to row %d, %s to row %d. The primary arm's row is the answer you fixed "
              "beforehand, but report the disagreement with it: the answer depends on how the gaps in the gauge are handled."
              % (prim["arm"], prim["row_best"], oth["arm"], oth["row_best"]))


# ------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------
def _style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def _cc_axis(ax):
    import matplotlib
    lo, hi = RS.BOX["channelconductivity_mmhr"][:2]
    ax.set_xscale("log")
    ax.set_xlim(lo * 0.8, hi * 1.25)
    ax.set_xticks([30, 100, 300, 1000])
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: "%g" % v))
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())


def fig_hydrographs(real_scored, sim_on, sim_off, row_on, row_off, path, arm, record_pts=None):
    """real_scored: the arm's real series. sim_on / sim_off: the best runs on the FULL stored grid, so a simulated
    flow in bins the arm cannot see is visible. record_pts: the as-is bins (drawn as markers on the filled arm)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    fig, ax = plt.subplots(figsize=(9.4, 4.8), facecolor=SURFACE)
    _style(ax)
    ax.plot(sim_off.index, sim_off.values, color=NULLC, linewidth=1.7, linestyle="--",
            label="best cc-OFF run (KGE_2012 %.3f)" % row_off["best_kge_2012"])
    ax.plot(sim_on.index, sim_on.values, color=C_CC, linewidth=1.9,
            label="best cc-ON run (KGE_2012 %.3f, cc %s mm/hr)" % (row_on["best_kge_2012"], fcc(row_on["best_cc"])))
    if arm == "asis":
        ax.plot(real_scored.index, real_scored.values, color=INK, linewidth=0, marker="o", markersize=3.5,
                label="real gauge, as-is: the %d bins the scorer keeps" % len(real_scored))
    else:
        ax.plot(real_scored.index, real_scored.values, color=INK, linewidth=1.5,
                label="real gauge, filled (%d points)" % len(real_scored))
        if record_pts is not None and len(record_pts):
            ax.plot(record_pts.index, record_pts.values, color=INK, linewidth=0, marker="o", markersize=3.0,
                    label="bins that hold a gauge record (%d)" % len(record_pts))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %H:%M"))
    ax.set_xlabel("time (day of month and hour, August 2014)", color=INK2, fontsize=9)
    ax.set_ylabel("discharge at the outlet (m$^3$/s)", color=INK2, fontsize=9)
    ax.legend(loc="upper right", frameon=False, fontsize=8.5, labelcolor=INK2)
    ax.set_title("%s arm: real gauge against the best run of each series (one event; runs drawn on the full grid)"
                 % arm.upper(), color=INK, fontsize=10, loc="left")
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def fig_observed(raw_w, real_valid, filled, ftype, path):
    """The data only, no fit: raw records, the as-is bins, the filled series. Whole window and a zoom on the flood."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from matplotlib.lines import Line2D
    raw = raw_w["Observed_CMS"].dropna()
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.6), facecolor=SURFACE, gridspec_kw={"width_ratios": [1.15, 1.0]})
    pos = filled[filled > 0]
    if len(pos):
        z0 = pos.index[0] - pd.Timedelta(minutes=30)
        z1 = pos.index[-1] + pd.Timedelta(minutes=30)
    else:
        z0, z1 = filled.index[0], filled.index[-1]
    for k, ax in enumerate(axes):
        _style(ax)
        ax.plot(filled.index, filled.values, color=INK, linewidth=1.6, zorder=2)
        hold = ftype[ftype == "hold"].index
        if len(hold) and k == 1:
            ax.plot(hold, filled.reindex(hold).values, color=INK, linewidth=0, marker="s", markersize=3.0, zorder=2)
        ax.plot(real_valid.index, real_valid.values, color=INK, linewidth=0, marker="o", markersize=5.0 if k else 3.5,
                markerfacecolor=SURFACE, markeredgewidth=1.2, zorder=3)
        ax.plot(raw.index, raw.values, color=INK2, linewidth=0, marker=".", markersize=4.0 if k else 2.5, zorder=4)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %H:%M" if k == 0 else "%H:%M"))
        ax.set_xlabel("time (day of month and hour, August 2014)" if k == 0 else "time, zoom on the flood (hour)",
                      color=INK2, fontsize=9)
        if k == 0:
            ax.set_ylabel("discharge at the outlet (m$^3$/s)", color=INK2, fontsize=9)
            ax.set_xlim(filled.index[0], filled.index[-1])
        else:
            ax.set_xlim(z0, z1)
            sub = raw[(raw.index >= z0) & (raw.index <= z1)]
            top = max(float(filled.reindex(pd.date_range(z0.ceil("5min"), z1, freq="5min")).max()),
                      float(sub.max()) if len(sub) else 0.0)
            ax.set_ylim(-0.04 * top, 1.08 * top)
    handles = [Line2D([0], [0], marker=".", color=INK2, linestyle="none", markersize=6, label="raw gauge records (%d in the window)" % len(raw)),
               Line2D([0], [0], marker="o", color=INK, markerfacecolor=SURFACE, markeredgewidth=1.2, linestyle="none", markersize=6,
                      label="as-is: the %d bins the scorer keeps" % len(real_valid)),
               Line2D([0], [0], color=INK, linewidth=1.6, label="filled series (%d points; squares in the zoom = held bins)" % int(filled.notna().sum()))]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=8.5, labelcolor=INK2)
    fig.suptitle("What gets scored: the real gauge as an event-triggered record, as-is and filled", color=INK, fontsize=10.5, x=0.01,
                 ha="left")
    fig.tight_layout(rect=(0, 0.07, 1, 0.95))
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def fig_kge_vs_cc(m_on, m_off, top_on, path, top_n, arm):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    both = pd.concat([m_on["kge_2012"], m_off["kge_2012"]])
    both = both[np.isfinite(both)]
    ymax = float(both.max())
    ymin = max(float(both.min()), ymax - 1.0)
    pad = 0.05 * (ymax - ymin)
    hidden = int((both < ymin).sum())
    lo_log, hi_log = cc_box_logs()
    edge = 10 ** (lo_log + (hi_log - lo_log) / 3.0)
    lo, hi = RS.BOX["channelconductivity_mmhr"][:2]
    fig, ax = plt.subplots(figsize=(9.0, 5.2), facecolor=SURFACE)
    _style(ax)
    ax.axvspan(lo * 0.8, edge, color=NULLC, alpha=0.13, linewidth=0, zorder=1)
    ax.scatter(m_off["channelconductivity_mmhr"], m_off["kge_2012"], s=22, marker="s", facecolors=SURFACE,
               edgecolors=NULLC, linewidths=1.1, zorder=2)
    ax.scatter(m_on["channelconductivity_mmhr"], m_on["kge_2012"], s=22, marker="o", color=C_CC, alpha=0.75,
               linewidths=0, zorder=3)
    ax.scatter(top_on["channelconductivity_mmhr"], top_on["kge_2012"], s=70, marker="o", color=C_CC,
               edgecolors=INK, linewidths=1.4, zorder=4)
    _cc_axis(ax)
    ax.set_ylim(ymin - pad, ymax + pad)
    ax.set_xlabel("channel conductivity cc (mm/hr)", color=INK2, fontsize=9)
    ax.set_ylabel("KGE_2012 against the real gauge (%s arm)" % arm, color=INK2, fontsize=9)
    handles = [
        Line2D([0], [0], marker="o", color=C_CC, linestyle="none", markersize=5, alpha=0.75, label="channel loss ON (one dot per run)"),
        Line2D([0], [0], marker="o", color=C_CC, markeredgecolor=INK, markeredgewidth=1.4, linestyle="none", markersize=8,
               label="the %d best ON runs" % len(top_on)),
        Line2D([0], [0], marker="s", color=NULLC, markerfacecolor=SURFACE, markeredgewidth=1.1, linestyle="none", markersize=5,
               label="control, loss OFF (same points; cc is inert)"),
        Patch(facecolor=NULLC, alpha=0.13, label="lowest third of the cc range"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=8, labelcolor=INK2)
    sub = "" if not hidden else "\n(%d run(s) below the plotted range are not shown)" % hidden
    ax.set_title("%s arm: where does the fit to the real gauge put cc?  (best %d ON runs ringed)%s"
                 % (arm.upper(), top_n, sub), color=INK, fontsize=10, loc="left")
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def fig_paired(paired, good_flag, good_frac, path, arm):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    d_good = paired.loc[good_flag, "diff"]
    lim = max(float(np.nanpercentile(np.abs(d_good), 98)) if len(d_good) else 0.05, 0.05) * 1.35
    inview = np.abs(paired["diff"]) <= lim
    hidden = int((~inview).sum())
    fig, ax = plt.subplots(figsize=(9.0, 5.2), facecolor=SURFACE)
    _style(ax)
    ax.axhline(0.0, color=INK2, linestyle=":", linewidth=1.1, zorder=1)
    rest = paired[~good_flag]
    ax.scatter(rest["channelconductivity_mmhr"], rest["diff"], s=20, color=NULLC, alpha=0.55, linewidths=0, zorder=2)
    good = paired[good_flag]
    ax.scatter(good["channelconductivity_mmhr"], good["diff"], s=34, color=C_CC, alpha=0.9, linewidths=0, zorder=3)
    u = cc_unit(good["channelconductivity_mmhr"].to_numpy(float))
    gd = good["diff"].to_numpy(float)
    lo_log, hi_log = cc_box_logs()
    xs, ys, nbins = [], [], 5
    for nb in (5, 3):
        edges = np.linspace(0.0, 1.0, nb + 1)
        edges[-1] = 1.0000001
        xs, ys = [], []
        for a, b in zip(edges[:-1], edges[1:]):
            sel = (u >= a) & (u < b)
            if sel.sum() >= 3:
                xs.append(10 ** (lo_log + (hi_log - lo_log) * (a + min(b, 1.0)) / 2.0))
                ys.append(float(np.median(gd[sel])))
        nbins = nb
        if len(xs) >= 3:
            break
    draw_median = len(xs) >= 2
    if draw_median:
        ax.plot(xs, ys, color=INK, linewidth=1.8, marker="D", markersize=5, zorder=4)
    _cc_axis(ax)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("channel conductivity cc of the ON run (mm/hr)", color=INK2, fontsize=9)
    ax.set_ylabel("KGE_2012: ON minus OFF at the same (Ks, f)", color=INK2, fontsize=9)
    handles = [
        Line2D([0], [0], marker="o", color=NULLC, linestyle="none", markersize=5, alpha=0.7, label="other pairs"),
        Line2D([0], [0], marker="o", color=C_CC, linestyle="none", markersize=6,
               label="good region: top %g%% of points by the cc-OFF fit (%d pairs)" % (100 * good_frac, len(good))),
        Line2D([0], [0], color=INK2, linestyle=":", linewidth=1.1, label="0 = no difference"),
    ]
    if draw_median:
        handles.insert(2, Line2D([0], [0], marker="D", color=INK, linewidth=1.8, markersize=5,
                                 label="median of the good region, per %s of the cc range (bins with 3 or more pairs)"
                                 % ("fifth" if nbins == 5 else "third")))
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=8, labelcolor=INK2)
    sub = "" if not hidden else "\n(%d pair(s) beyond the plotted range are not shown)" % hidden
    ax.set_title("%s arm: does channel loss help or hurt the fit as cc rises?  (positive = ON fits better)%s"
                 % (arm.upper(), sub), color=INK, fontsize=10, loc="left")
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


# ------------------------------------------------------------------
# Scoring and analysis of one arm
# ------------------------------------------------------------------
def score_arm(arm, valid, ref_idx, res_on, res_off, cache, phase_fn):
    """Score every stored run against one arm's real series (only the bins that arm has, and that the stored grid has)."""
    scored_idx = valid.index.intersection(ref_idx)
    real_scored = valid.loc[scored_idx]
    cache_g = {}
    for rid in list(res_on["run_id"]) + list(res_off["run_id"]):
        cache_g[rid] = pd.DataFrame({"Observed": real_scored, "Simulated": cache[rid]["Simulated"].reindex(scored_idx)})
    m_on, f_on = TL.score_truth(real_scored, res_on, cache_g, phase_fn)
    m_off, f_off = TL.score_truth(real_scored, res_off, cache_g, phase_fn)
    return {"arm": arm, "scored_idx": scored_idx, "real_scored": real_scored, "m_on": m_on, "m_off": m_off,
            "n_failed": int(len(f_on) + len(f_off))}


def analyze_arm(ar, args):
    """Everything the report needs for one arm. Pure function of the two metric tables and the settings."""
    m_on, m_off = ar["m_on"], ar["m_off"]
    main_n = args.top_n[0]
    srows = [series_summary(m_on, SER_ON, main_n), series_summary(m_off, SER_OFF, main_n)]
    on_row, off_row = srows[0], srows[1]
    diff_row = {"series": "ON minus OFF"}
    for key in ("best_kge_2012", "top_n_mean_kge", "top_n_median_kge", "median_kge_all_runs", "n_feasible",
                "best_pbias_pct", "best_kge_r", "best_kge_beta", "best_kge_gamma"):
        diff_row[key] = on_row[key] - off_row[key]
    sdf = pd.concat([pd.DataFrame(srows), pd.DataFrame([diff_row])], ignore_index=True)
    loc_rows = [cc_location(m_on, n, args.floor_tol) for n in args.top_n]
    paired, pnote = pair_runs(m_on, m_off)
    pst, good_flag = None, None
    if paired is not None:
        good_flag = mark_good_region(paired, args.good_frac)
        paired["in_good_region"] = good_flag
        pst = paired_stats(paired, good_flag, args.good_frac, args.boot, args.seed)
    gap = on_row["best_kge_2012"] - off_row["best_kge_2012"]
    gap_mean = on_row["top_n_mean_kge"] - off_row["top_n_mean_kge"]
    an = {"sdf": sdf, "on_row": on_row, "off_row": off_row, "loc_rows": loc_rows, "paired": paired, "pnote": pnote,
          "pst": pst, "good_flag": good_flag, "gap": gap, "gap_mean": gap_mean}
    if args.margin is not None:
        an["row_best"] = decision_row(gap, args.margin, loc_rows[0])
        an["row_mean"] = decision_row(gap_mean, args.margin, loc_rows[0])
    return an


def print_arm_report(arm, an, margin, n_scored, main_n, is_primary):
    print("\n" + "#" * 110)
    print("ARM: %s   %s" % (arm.upper(), "(PRIMARY: fixed before the run)" if is_primary else "(sensitivity arm)"))
    print("     %s" % ARM_TEXT[arm])
    print("#" * 110)
    print_summary(an["sdf"], margin, n_scored, arm)
    print_cc_location(an["loc_rows"], main_n, arm)
    if an["pst"] is not None:
        print_paired(an["pst"], an["pnote"], arm)
    else:
        print("\n  Paired analysis not possible: %s" % an["pnote"])
    print_reading(an["gap"], an["gap_mean"], margin, an["loc_rows"][0], an["row_best"], an["row_mean"], arm)


def write_arm(adir, arm, ar, an, args, cache, ref_idx, record_pts):
    adir.mkdir(parents=True, exist_ok=True)
    m_on, m_off = ar["m_on"], ar["m_off"]
    on_row, off_row = an["on_row"], an["off_row"]
    long_df = pd.concat([m_on.assign(series=SER_ON), m_off.assign(series=SER_OFF)], ignore_index=True)
    write_csv_atomic(long_df, adir / "real_gauge_long_110.csv")
    write_csv_atomic(an["sdf"], adir / "real_gauge_summary_110.csv")
    top_keep = ["series", "rank", "run_id", "Ks_mult", "f_RS_abs", "channelconductivity_mmhr", "kge_2012", "kge_r", "kge_beta",
                "kge_gamma", "pbias_pct", "nse", "rmse_m3s", "peak_error_pct", "peak_timing_error_hr"]
    top_tabs = []
    for series, m in ((SER_ON, m_on), (SER_OFF, m_off)):
        t = top_rows(m, max(args.top_n)).reset_index(drop=True)
        t["rank"] = np.arange(1, len(t) + 1)
        t["series"] = series
        top_tabs.append(t[top_keep])
    write_csv_atomic(pd.concat(top_tabs, ignore_index=True), adir / "real_gauge_top_runs_110.csv")
    write_csv_atomic(pd.DataFrame(an["loc_rows"]), adir / "real_gauge_cc_location_110.csv")
    if an["paired"] is not None:
        write_csv_atomic(an["paired"], adir / "real_gauge_paired_110.csv")
        write_csv_atomic(an["pst"], adir / "real_gauge_paired_stats_110.csv")
    if args.no_plots:
        return
    main_n = args.top_n[0]
    best_on_sim = cache[on_row["best_run_id"]]["Simulated"].reindex(ref_idx)
    best_off_sim = cache[off_row["best_run_id"]]["Simulated"].reindex(ref_idx)
    top_on = top_rows(m_on, main_n)
    figs = [("fig_real_gauge_hydrographs_110.png",
             lambda p: fig_hydrographs(ar["real_scored"], best_on_sim, best_off_sim, on_row, off_row, p, arm,
                                       record_pts if arm == "filled" else None)),
            ("fig_real_gauge_kge_vs_cc_110.png", lambda p: fig_kge_vs_cc(m_on, m_off, top_on, p, main_n, arm))]
    if an["paired"] is not None:
        figs.append(("fig_real_gauge_paired_diff_110.png",
                     lambda p: fig_paired(an["paired"], an["good_flag"], args.good_frac, p, arm)))
    for name, fn in figs:
        try:
            fn(adir / name)
        except Exception as e:
            print("\n  (figure %s/%s skipped: %s -- the CSVs are saved regardless)" % (arm, name, e))


def write_observed(out_dir, real_grid, filled, ftype, scored_by_arm, raw_w, valid_asis, no_plots):
    out_dir.mkdir(parents=True, exist_ok=True)
    obs_out = pd.DataFrame({"datetime": filled.index, "Observed_asis_m3s": real_grid.reindex(filled.index).values,
                            "Observed_filled_m3s": filled.values,
                            "Observed_filled_cfs": filled.values / CFS_TO_CMS,
                            "fill_type": ftype.values})
    for a in ("asis", "filled"):
        idx = scored_by_arm.get(a)
        obs_out["scored_" + a] = obs_out["datetime"].isin(idx) if idx is not None else False
    write_csv_atomic(obs_out, out_dir / "real_gauge_observed_110.csv")
    if not no_plots:
        try:
            fig_observed(raw_w, valid_asis, filled, ftype, out_dir / "fig_real_gauge_observed_110.png")
        except Exception as e:
            print("\n  (figure fig_real_gauge_observed_110.png skipped: %s -- the CSV is saved regardless)" % e)


# ------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Series 110, Stage 0 -- re-score the 500 stored Aug 12 2014 runs against the REAL gauge (no tRIBS).")
    ap.add_argument("--margin", type=float, default=None,
                    help="REQUIRED for the full run: how much better (KGE_2012) ON must be than OFF to count as a real "
                         "gain. Fix it before looking at any result (handoff section 4).")
    ap.add_argument("--gauge_mode", choices=["both", "asis", "filled"], default="both",
                    help="Which way to score the gauge: 'asis' (the scorer's own 5-minute bins that hold a record), 'filled' "
                         "(all bins, gaps filled from the record), or 'both' (default). See the top of this file.")
    ap.add_argument("--primary_arm", choices=["asis", "filled"], default=None,
                    help="REQUIRED for the full run when --gauge_mode is 'both': the arm whose reading is the answer. Fix it "
                         "before looking at any result, like the margin.")
    ap.add_argument("--max_gap_min", type=float, default=DEFAULT_MAX_GAP_MIN,
                    help="Filled arm: draw a straight line between two readings at most this many minutes apart (default %g)"
                         % DEFAULT_MAX_GAP_MIN)
    ap.add_argument("--max_hold_hr", type=float, default=DEFAULT_MAX_HOLD_HR,
                    help="Filled arm: across a longer silence, hold the last reading for at most this many hours (default %g)"
                         % DEFAULT_MAX_HOLD_HR)
    ap.add_argument("--check_only", action="store_true", help="Run the gates and stop; writes nothing, shows no fit.")
    ap.add_argument("--observed_only", action="store_true",
                    help="Run the gates, then write only the observed series (CSV and figure); shows no fit, needs no margin.")
    ap.add_argument("--gate3_compare", type=Path, default=None,
                    help="*_compare_obs_sim.csv of ONE stored parameter point re-run through the scorer in real-gauge mode")
    ap.add_argument("--gate3_metrics", type=Path, default=None,
                    help="*_metrics_summary.csv of that same gauge-mode run")
    ap.add_argument("--top_n", type=int, nargs="+", default=[10, 5, 25],
                    help="How many best ON runs to look at for cc; the first is the main one (default 10 5 25)")
    ap.add_argument("--floor_tol", type=float, default=FLOOR_TOL_MMHR,
                    help="'At the floor' means within this many mm/hr of the lowest cc in the box (default %g)" % FLOOR_TOL_MMHR)
    ap.add_argument("--good_frac", type=float, default=0.2,
                    help="Share of points (best by cc-OFF fit) forming the 'good region' (default 0.2)")
    ap.add_argument("--boot", type=int, default=2000, help="Bootstrap re-draws for the slope interval (default 2000)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--min_points", type=int, default=30,
                    help="Fewest scorable 5-minute points to accept for the filled arm; fewer as-is points only warns (default 30)")
    ap.add_argument("--expect_n", type=int, default=RS.DEFAULT_EXPECT_N,
                    help="Expected runs per series (default %d); a different count only warns" % RS.DEFAULT_EXPECT_N)
    ap.add_argument("--drop_bad_runs", action="store_true",
                    help="Exclude runs that fail gate 2A, 2B or 4 and continue, instead of stopping.")
    ap.add_argument("--no_plots", action="store_true")
    ap.add_argument("--stage2_csv", type=Path, default=None)
    ap.add_argument("--control_csv", type=Path, default=None)
    ap.add_argument("--gauge_xlsx", type=Path, default=None)
    ap.add_argument("--out_dir", type=Path, default=None)
    args = ap.parse_args()

    if args.check_only and args.observed_only:
        ap.error("use --check_only or --observed_only, not both.")
    full_run = not (args.check_only or args.observed_only)
    if full_run:
        if args.margin is None:
            ap.error("--margin is required for the full run. Handoff section 4: fix the decision rule BEFORE looking at "
                     "any real-data result. Give the KGE_2012 difference (ON minus OFF) you would need to see before "
                     "calling it a real gain, for example --margin 0.05. (--check_only needs no margin and shows no fit.)")
        if not np.isfinite(args.margin) or args.margin < 0:
            ap.error("--margin must be zero or a positive number (KGE_2012 units).")
        if args.gauge_mode == "both" and args.primary_arm is None:
            ap.error("--primary_arm is required when --gauge_mode is 'both' (asis or filled). Like the margin, fix it BEFORE "
                     "looking at any real-data result: it names the arm whose reading is the answer; the other is a "
                     "sensitivity check.")
    if args.gauge_mode != "both" and args.primary_arm not in (None, args.gauge_mode):
        ap.error("--primary_arm %s contradicts --gauge_mode %s." % (args.primary_arm, args.gauge_mode))
    if args.gate3_metrics is not None and args.gate3_compare is None:
        ap.error("--gate3_metrics needs --gate3_compare as well.")
    if not (0.0 < args.good_frac < 1.0):
        ap.error("--good_frac must be between 0 and 1.")
    if args.floor_tol <= 0:
        ap.error("--floor_tol must be positive.")
    if any(n < 3 for n in args.top_n):
        ap.error("--top_n values must be at least 3.")
    if args.max_gap_min <= 0 or args.max_hold_hr <= 0:
        ap.error("--max_gap_min and --max_hold_hr must be positive.")

    arms_used = ["asis", "filled"] if args.gauge_mode == "both" else [args.gauge_mode]
    primary = args.primary_arm or (args.gauge_mode if args.gauge_mode != "both" else None)
    if primary is not None:
        arms_used = [primary] + [a for a in arms_used if a != primary]

    script_dir = Path.cwd()
    project_root = script_dir.parent
    calib_dir = project_root / "calibration_work"
    summary_dir = calib_dir / "03_comparisons" / "summary_tables"
    csv_dir = calib_dir / "03_comparisons" / "csv_exports"
    synth_dir = calib_dir / "synth_truth"
    out_dir = args.out_dir or (summary_dir / OUT_DIRNAME)
    p_on = args.stage2_csv or (summary_dir / RS.STAGE2_NAME)
    p_off = args.control_csv or (summary_dir / RS.CONTROL_NAME)
    xlsx_path = args.gauge_xlsx or project_root.joinpath(*GAUGE_REL)

    print("\n" + "=" * 78)
    print("Series 110, Stage 0 (version 2) -- re-score the stored Aug 12 2014 runs against the REAL gauge")
    print("python %s, pandas %s, numpy %s" % (sys.version.split()[0], pd.__version__, np.__version__))
    print("arms: %s%s" % (", ".join(arms_used), "   (primary: %s)" % primary if primary else ""))
    print("=" * 78)
    if ".".join(pd.__version__.split(".")[:2]) != EXPECT_PANDAS_MINOR:
        print("  *** NOTE: pandas here is %s; the runs were scored with pandas %s.x. The 5-minute interpolation differs "
              "between pandas 2 and 3. Gate 2A will tell. ***" % (pd.__version__, EXPECT_PANDAS_MINOR))

    # ---- inputs / safety -------------------------------------------------------------
    qouts = list(synth_dir.glob("*.qout")) if synth_dir.exists() else []
    if len(qouts) != 1:
        sys.exit("Expected exactly one *.qout in %s (the cc-OFF truth the 500 runs were scored against), found %d: %s"
                 % (synth_dir, len(qouts), [q.name for q in qouts]))
    ccoff_path = qouts[0]
    print("cc-OFF synthetic truth (read-only, used only for gate 2): %s" % ccoff_path.name)
    if not xlsx_path.exists():
        sys.exit("The real gauge workbook was not found: %s\nIt is read from %s relative to the project root "
                 "(lab/'s parent); --gauge_xlsx points elsewhere." % (xlsx_path, "/".join(GAUGE_REL)))

    out_res = out_dir.resolve()
    for guarded in (synth_dir, csv_dir):
        gd = guarded.resolve()
        if out_res == gd or gd in out_res.parents or out_res in gd.parents:
            sys.exit("Output folder %s overlaps an input folder (%s); refusing." % (out_dir, guarded))
    for pth in (p_on, p_off, xlsx_path):
        rp = pth.resolve()
        if out_res in rp.parents or out_res == rp:
            sys.exit("Output folder %s contains an input file (%s); refusing." % (out_dir, pth))

    print("\nLoading results tables:")
    res_on = RS.load_results(p_on, SER_ON, 1, args.expect_n)
    res_off = RS.load_results(p_off, SER_OFF, 0, args.expect_n)
    phase_fn = RS.import_phase_fn(script_dir)
    if phase_fn is None:
        sys.exit("The scorer's phase-metric function could not be imported (run_sensitivity_single_interp_tribs6.py must "
                 "be in this folder). Without it the phase metrics would not be the scorer's, so gate 2 cannot pass.")
    obs_ccoff = RS.read_truth_5min(ccoff_path)

    print("\nLoading compare CSVs (run membership = run_id rows above):")
    cache_on, prob_on = RS.load_all_compares(res_on, csv_dir)
    cache_off, prob_off = RS.load_all_compares(res_off, csv_dir)
    print("  %s: %d/%d loaded;  %s: %d/%d loaded" % (SER_ON, len(cache_on), len(res_on), SER_OFF, len(cache_off), len(res_off)))
    cache = {}
    cache.update(cache_on)
    cache.update(cache_off)
    bad = {}
    bad.update(prob_on)
    bad.update(prob_off)

    # ---- GATE 2: scorer fidelity (rescore_cc_truths_110's gates A and B) ----------------------
    allres = pd.concat([res_on.assign(_series=SER_ON), res_off.assign(_series=SER_OFF)], ignore_index=True)
    badA, worstA = RS.gate_a(allres, cache, obs_ccoff)
    badB, worstB = RS.gate_b(allres, cache, phase_fn)
    for dct in (badA, badB):
        for rid, why in dct.items():
            bad.setdefault(rid, why)
    n_checked = len(cache)
    ck = Checks()
    ck.add("PASS" if not badA else "FAIL", "GATE 2A synthetic Observed",
           "%d/%d compare CSVs match the cc-OFF truth as read here (worst abs diff %.2e m3/s)"
           % (n_checked - len(badA), n_checked, worstA))
    if worstB:
        wk = max(worstB, key=worstB.get)
        ck.add("PASS" if not badB else "FAIL", "GATE 2B stored metrics",
               "%d/%d runs reproduce their stored metrics (worst abs deviation %.2e in %s; phase metrics included)"
               % (n_checked - len(badB), n_checked, worstB[wk], wk))
    else:
        ck.add("PASS" if not badB else "FAIL", "GATE 2B stored metrics", "%d/%d" % (n_checked - len(badB), n_checked))
    if n_checked and len(badA) >= 0.9 * n_checked:
        print("\n  GATE 2A failed for nearly every run. The usual cause is a different pandas version from the one that scored "
              "the runs (here: pandas %s). Use the same python as for the runs. Do NOT use --drop_bad_runs for this."
              % pd.__version__)
    missing_files = {r: w for r, w in bad.items() if "missing" in w or "unreadable" in w}
    if missing_files:
        ck.add("FAIL", "compare CSVs", "%d compare CSV(s) missing or unreadable, e.g. %s"
               % (len(missing_files), ", ".join(list(missing_files)[:3])))

    # ---- GATE 1: the real series (as-is) and 1b (filled) ----------------------------------------
    print("\nReading the real gauge workbook (the whole 1993-2025 record is read, as the scorer does; this can take a minute) ...",
          flush=True)
    try:
        obs_df = read_gauge_workbook(xlsx_path)
    except Exception as e:
        sys.exit("Could not read the real gauge workbook the way the scorer does (%s: %s). Expected sheet %r with %d header "
                 "rows skipped and columns Date, Time, cfs." % (type(e).__name__, str(e)[:120], GAUGE_SHEET, GAUGE_SKIPROWS))
    real_5min = gauge_to_grid(obs_df)
    peaks = [float(np.nanmax(df["Simulated"].to_numpy(dtype=float))) for df in cache.values()
             if len(df) and np.isfinite(df["Simulated"].to_numpy(dtype=float)).any()]
    sim_pk = float(np.median(peaks)) if peaks else None
    ck1, info1, real_valid, real_grid = check_real_series(obs_df, real_5min, xlsx_path, args.min_points, sim_pk)
    gate3_wanted = args.gate3_compare is not None
    if "asis" not in arms_used and not gate3_wanted:
        for r in ck1.rows:
            if r["status"] == "FAIL" and r["check"] == "as-is window coverage":
                r["status"] = "WARN"
                r["detail"] += "  [as-is arm not in use]"
    filled, ftype, finfo = build_filled(obs_df, real_5min, args.max_gap_min, args.max_hold_hr)
    t0w, t1w = pd.Timestamp(EVENT_START), pd.Timestamp(EVENT_END)
    raw_w = obs_df.loc[np.asarray((obs_df.index >= t0w) & (obs_df.index <= t1w))]
    filled_valid = filled.dropna()
    ck1f, infof = check_filled(filled, ftype, finfo, real_valid, raw_w, args.min_points)
    if "filled" not in arms_used:
        for r in ck1f.rows:
            if r["status"] == "FAIL":
                r["status"] = "WARN"
                r["detail"] += "  [filled arm not in use]"
    ck.rows = ck1.rows + ck1f.rows + ck.rows
    if ck1.n("FAIL") or ck1f.n("FAIL"):
        print_gates(ck.rows, "not reached")
        sys.exit("\nSTOPPED at gate 1: the real series is not usable (see FAIL above). Nothing was scored or written.")

    # ---- stop here if gate 2 failed, before spending time on gate 4 --------------------------------
    run_ids = [r for r in list(res_on["run_id"]) + list(res_off["run_id"]) if r in cache]
    if bad and not args.drop_bad_runs:
        print_gates(ck.rows, "not reached")
        print("\n  %d run(s) failed a gate or could not be loaded. First few:" % len(bad))
        for rid, why in list(bad.items())[:8]:
            print("    %s: %s" % (rid, why))
        sys.exit("\nSTOPPED: nothing was scored. Fix the cause (or re-run with --drop_bad_runs to exclude these runs and continue).")
    if bad:
        for r in ck.rows:
            if r["status"] == "FAIL" and r["check"] in DROPPABLE:
                r["status"] = "WARN"
                r["detail"] += "  [runs excluded by --drop_bad_runs]"
    run_ids = [r for r in run_ids if r not in bad]
    if not run_ids:
        sys.exit("No runs left after the gates; stopping.")

    # ---- GATE 4: the grid ------------------------------------------------------------------------------
    valids = {"asis": real_valid, "filled": filled_valid}
    grid_arms = [a for a in ("asis", "filled") if a in arms_used or (a == "asis" and gate3_wanted)]
    union_valid = pd.concat([valids[a] for a in grid_arms]).groupby(level=0).first().sort_index()
    ck4, ref_idx, scored_union, bad4 = check_grid(cache, run_ids, union_valid)
    if bad4:
        if not args.drop_bad_runs:
            ck.extend(ck4)
            print_gates(ck.rows, "not reached")
            sys.exit("\nSTOPPED at gate 4: stored grids differ between runs. Nothing was scored.")
        for rid, why in bad4.items():
            bad.setdefault(rid, why)
        run_ids = [r for r in run_ids if r not in bad4]
        ck.add("WARN", "GATE 4 runs excluded", "%d run(s) with a different stored grid were excluded (--drop_bad_runs)" % len(bad4))
        ck4, ref_idx, scored_union, bad4 = check_grid(cache, run_ids, union_valid)
    ck.extend(ck4)
    if ck4.n("FAIL"):
        print_gates(ck.rows, "not reached")
        sys.exit("\nSTOPPED at gate 4: the real 5-minute grid does not line up with the stored runs. Nothing was scored.")
    res_on = res_on[res_on["run_id"].isin(run_ids)].reset_index(drop=True)
    res_off = res_off[res_off["run_id"].isin(run_ids)].reset_index(drop=True)
    if res_on.empty or res_off.empty:
        sys.exit("No runs left in one of the series after the gates; stopping.")
    if bad:
        ck.add("WARN", "runs excluded", "%d run(s) excluded by --drop_bad_runs" % len(bad))

    # ---- score every arm in use ---------------------------------------------------------------------------
    arm_res = {}
    for a in grid_arms:
        if len(valids[a]) < (MIN_ASIS_POINTS if a == "asis" else 3):
            continue
        arm_res[a] = score_arm(a, valids[a], ref_idx, res_on, res_off, cache, phase_fn)
        if arm_res[a]["n_failed"]:
            print("  WARNING (%s arm): %d run(s) could not be scored (missing simulated values on the real grid)"
                  % (a, arm_res[a]["n_failed"]))
        if arm_res[a]["m_on"].empty or arm_res[a]["m_off"].empty:
            sys.exit("Nothing could be scored in one of the series (%s arm)." % a)
    for a in arms_used:
        if a not in arm_res:
            sys.exit("The %s arm could not be scored (too few real points)." % a)

    # ---- GATE 3 (the as-is path) -----------------------------------------------------------------------------
    gate3_state = ("NOT DONE -- results are PROVISIONAL until one run has been scored by the actual scorer in real-gauge mode"
                   " (it also vouches for how the workbook's time stamps are read, which both arms share)")
    gate3_info = {}
    if gate3_wanted:
        ar_a = arm_res["asis"]
        m_all = pd.concat([ar_a["m_on"], ar_a["m_off"]], ignore_index=True)
        ck3, gate3_info = verify_gate3(args.gate3_compare, args.gate3_metrics, allres[allres["run_id"].isin(run_ids)],
                                       m_all, cache, ar_a["real_scored"], ar_a["scored_idx"])
        ck.extend(ck3)
        gate3_state = "FAILED" if ck3.n("FAIL") else ("PASSED" if not ck3.n("WARN") else "PASSED WITH A WARNING")
    print_gates(ck.rows, gate3_state)
    if ck.n("FAIL"):
        sys.exit("\nSTOPPED: a gate failed (see FAIL above). Nothing was written.")

    # ---- blind share (stored simulations and the observed series only) -----------------------------------------
    asis_idx = valids["asis"].index.intersection(ref_idx)
    series_of = {}
    for rid in res_on["run_id"]:
        series_of[rid] = SER_ON
    for rid in res_off["run_id"]:
        series_of[rid] = SER_OFF
    bs, flood0, flood1 = blind_share(cache, run_ids, series_of, ref_idx, asis_idx, filled)
    print_blind_share(bs, flood0, flood1, len(asis_idx), len(ref_idx))

    scored_by_arm = {a: arm_res[a]["scored_idx"] for a in arm_res}
    if args.check_only:
        print("\n--check_only: gates finished, nothing written, no fit shown. A full run scores %d + %d runs: %s."
              % (len(res_on), len(res_off), "; ".join("%s arm on %d real points" % (a, len(arm_res[a]["scored_idx"]))
                                                      for a in arms_used)))
        return
    if args.observed_only:
        write_observed(out_dir, real_grid, filled, ftype, scored_by_arm, raw_w, valids["asis"], args.no_plots)
        if bs is not None:
            write_csv_atomic(bs, out_dir / "real_gauge_blind_share_110.csv")
        write_csv_atomic(pd.DataFrame(ck.rows), out_dir / "real_gauge_gates_110.csv")
        print("\n--observed_only: wrote the observed series (no fit, no scoring output) to %s" % out_dir)
        print("  real_gauge_observed_110.csv   fig_real_gauge_observed_110.png   real_gauge_blind_share_110.csv   real_gauge_gates_110.csv")
        return

    # ---- analysis ----------------------------------------------------------------------------------------------
    main_n = args.top_n[0]
    analyses = {a: analyze_arm(arm_res[a], args) for a in arms_used}

    # ---- report ----------------------------------------------------------------------------------------------------
    if gate3_state.startswith("NOT DONE"):
        print("\n" + "!" * 110)
        print("!! PROVISIONAL: gate 3 (the actual scorer in real-gauge mode) has not been run. Do not quote these numbers yet.")
        print("!" * 110)
    for a in arms_used:
        print_arm_report(a, analyses[a], args.margin, len(arm_res[a]["scored_idx"]), main_n, a == primary)
    ag_rows = []
    for a in arms_used:
        an = analyses[a]
        lc = an["loc_rows"][0]
        ag_rows.append({"arm": a, "primary": a == primary, "n_points": int(len(arm_res[a]["scored_idx"])),
                        "best_kge_on": float(an["on_row"]["best_kge_2012"]), "best_kge_off": float(an["off_row"]["best_kge_2012"]),
                        "gap_best": float(an["gap"]), "gap_top_n_mean": float(an["gap_mean"]),
                        "row_best": int(an["row_best"]), "row_mean": int(an["row_mean"]),
                        "best_cc": float(an["on_row"]["best_cc"]), "p_lowest_third": float(lc["p_lowest_third"]),
                        "span_frac_of_box": float(lc["span_frac_of_box"]), "margin": float(args.margin)})
    ag = pd.DataFrame(ag_rows)
    print_arm_agreement(ag, primary, args.margin)
    print("\n  Caveats (handoff section 5): one event; the SMPHQ rainfall artifact is in every run; routing is pinned at the "
          "synthetic-truth values; the box clipped the best-fit ridge earlier; gauge error is not quantified; the filled "
          "arm assumes a silent logger means no change.")

    # ---- write ----------------------------------------------------------------------------------------------------------
    write_observed(out_dir, real_grid, filled, ftype, scored_by_arm, raw_w, valids["asis"], args.no_plots)
    if bs is not None:
        write_csv_atomic(bs, out_dir / "real_gauge_blind_share_110.csv")
    write_csv_atomic(ag, out_dir / "real_gauge_arm_agreement_110.csv")
    write_csv_atomic(pd.DataFrame(ck.rows), out_dir / "real_gauge_gates_110.csv")
    record_pts = valids["asis"].loc[asis_idx]
    for a in arms_used:
        write_arm(out_dir / a, a, arm_res[a], analyses[a], args, cache, ref_idx, record_pts)

    inputs = {p_on.name: {"path": str(p_on), "md5": RS.md5_of(p_on), "runs": int(len(res_on))},
              p_off.name: {"path": str(p_off), "md5": RS.md5_of(p_off), "runs": int(len(res_off))}}
    prim_an = analyses[primary] if primary else None
    prov = {
        "script": Path(__file__).name,
        "version": 2,
        "created_local": datetime.now().isoformat(timespec="seconds"),
        "python": sys.version.split()[0], "pandas": pd.__version__, "numpy": np.__version__,
        "provisional": bool(gate3_state.startswith("NOT DONE")),
        "arms": {"used": arms_used, "primary": primary, "max_gap_min": args.max_gap_min, "max_hold_hr": args.max_hold_hr,
                 "filled_counts": finfo, "filled_checks": infof, "n_points": {a: int(len(arm_res[a]["scored_idx"])) for a in arms_used}},
        "decision": {"margin_kge_2012": args.margin, "alpha": ALPHA, "stable_span": STABLE_SPAN,
                     "floor_tol_mmhr": args.floor_tol, "main_top_n": main_n,
                     "per_arm": {a: {"gap_best": float(analyses[a]["gap"]), "gap_top_n_mean": float(analyses[a]["gap_mean"]),
                                     "row_best": int(analyses[a]["row_best"]), "row_top_n_mean": int(analyses[a]["row_mean"])}
                                 for a in arms_used},
                     "primary_row_best": int(prim_an["row_best"]) if prim_an else None,
                     "primary_reading": READINGS[prim_an["row_best"]] if prim_an else None},
        "blind_share": bs.to_dict(orient="records") if bs is not None else None,
        "real_series": info1,
        "cc_off_truth": {"file": ccoff_path.name, "md5": RS.md5_of(ccoff_path)},
        "inputs": inputs,
        "gates": {"2A_max_abs_diff_m3s": worstA, "2B_worst_abs_dev": worstB, "runs_checked": n_checked,
                  "runs_excluded": {k: v for k, v in bad.items()},
                  "3": {"state": gate3_state, "detail": gate3_info},
                  "grid_points_stored": int(len(ref_idx)),
                  "checks": ck.rows},
        "settings": {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()},
        "note": "Metrics recomputed from the stored Simulated column of each compare CSV against the real gauge series. "
                "as-is arm: read and resampled as the scorer does in real-gauge mode. filled arm: the same bin means plus "
                "the empty bins filled from the raw records (line, then hold). Scoring and metric code from "
                "rescore_truth_location_110.py and rescore_cc_truths_110.py.",
    }
    write_json_atomic(out_dir / "PROVENANCE_real_gauge_110.json", prov)

    print("\nSaved to: %s" % out_dir)
    print("  top folder: real_gauge_observed_110.csv   real_gauge_blind_share_110.csv   real_gauge_arm_agreement_110.csv   "
          "real_gauge_gates_110.csv")
    print("              fig_real_gauge_observed_110.png   PROVENANCE_real_gauge_110.json")
    print("  per arm (%s): real_gauge_summary_110.csv   real_gauge_top_runs_110.csv   real_gauge_cc_location_110.csv" %
          ", ".join("%s/" % a for a in arms_used))
    print("                real_gauge_paired_110.csv   real_gauge_paired_stats_110.csv   real_gauge_long_110.csv   fig_real_gauge_*.png")


if __name__ == "__main__":
    main()
