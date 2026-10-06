"""
run_lowf_sweep_110.py   (version 1.3)
===================================
Series 110 -- LOW-f sweep. One script for both series (channel loss OFF and ON).
Run it from the lab/ directory, like every other tribs6 script.

WHAT CHANGED IN VERSION 1.3 (for the multi-storm runs: one long window, storms cut out afterwards)
----------------------------------------------------------------------------------------------------
  * --runtime_hours H   runs every sample over a window of H hours from 1 Aug 2014 00:00 instead of the builder's
                  450 (to 19 Aug 18:00). 1464 hours reaches 1 Oct 2014 00:00 and so takes in 12 Aug, 19 Aug, 8 Sep
                  and 27 Sep. The Outlet .qout is written every 5 minutes for the whole window and is KEPT for
                  every run, so any storm in the window can be scored later without a new tRIBS run. The run is
                  still scored on the 12 Aug window (compare CSV and result row), exactly as before.
  * With --runtime_hours: the per-run timeout and the free-disk limit scale with the window (an explicit --timeout or
                  --min_free_gb still wins), the time estimate scales (an assumption: the first run prints the real
                  figure), and every run's Outlet .qout must reach the end of the window. A run whose .qout stops
                  early (tRIBS can end part-way and still exit cleanly) is reported FAILED and its raw folder is
                  kept; --skip_existing re-runs it.
  * The design file records the window, so a label cannot be reused with a different window by mistake.
  * Without --runtime_hours nothing changes (450 hours, timeout 300 s, no end-of-window check).
  * Not checked by this script: that the rainfall and weather forcing files run at least as long as the window.

WHAT CHANGED IN VERSION 1.2 (for Stage B, the ON-versus-OFF comparison)
------------------------------------------------------------------------
  * --f_fixed X   holds f_RS_abs at one value for every run instead of sampling it (use it instead of
                  --f_lo/--f_hi). Stage B fixes f at 0.0007: the first two sweeps showed the real-gauge fit is
                  flat in f from 0.0003 to about 0.0014, so f cannot be placed by this event.
  * Resume is stricter: with --skip_existing a run counts as DONE only if its compare CSV exists AND its Outlet
                  .qout was kept AND its row is in the results (or its metrics file exists). A run cut off between
                  those steps is re-run (it used to be skipped and silently lost its .qout). A finished run whose raw
                  folder was left behind is cleaned up.
  * The pairing check now only applies when the other series of the same label has the SAME design (same n, seed
                  and box). Stage B deliberately uses different designs for ON and OFF (OFF: Ks only; ON: Ks and
                  cc), so no point pairing is expected and none is required.
  * The time estimate is the measured 68 s per run (66 to 69 s on this project's runs), not 85 s.
  * Version 1.1 (kept): output is printed live even through a pipe or a log file.
Everything else is unchanged: the sampler, the builder and scorer calls, run names, output files.

WHY THIS SCRIPT EXISTS
-----------------------
Stage 0 scored the 500 stored Aug 12 2014 runs against the REAL gauge. The best fit sat at the
lowest f the sweep box allowed (f_RS_abs 0.004), in both series, and the fit was still
improving there. So the box was clipped: the real-data optimum probably lies at or below
f = 0.004 (synthetic truth: 0.012). Josh: f = 0.02 was only a starting estimate; lower or higher
is open to exploration. This script runs a NEW sweep with the f range (and the Ks and cc
ranges) set on the command line, so the box can be moved without editing code.

It is a copy of run_cc_stage2_control_110.py (OFF) and run_cc_stage2_joint_lhs_110.py (ON)
with these changes, and nothing else:
  1. The box comes from the command line. Defaults = the box suggested on 2026-10-04:
       Ks_mult 3.0 - 10.5 (linear)   f_RS_abs 0.001 - 0.006 (log)   cc 30 - 1000 mm/hr (log)
  2. --mode off|on picks the series (optpercolation 0 or 1). Same --label, --n, --seed and box
     give IDENTICAL sample points in both modes, so the two series pair point for point, as in
     Stage 2 and its control.
  3. Own run names and own output files, so nothing from Stage 2, the control or the
     "lf1" runs can be overwritten (run-id tag = label + mode, for example "lf1off").
  4. The OUTLET flow file (<run_id>_Outlet.qout, small) is COPIED to
     calibration_work/kept_qout_110/ before the bulky raw folder is deleted. The Stage 2 scripts
     deleted everything (about 270 MB per run), so no other storm in the 450-hour run could be
     scored afterwards. With the .qout kept, any window of 1 Aug 00:00 + 450 h (to 19 Aug 18:00)
     can be re-scored later without a new tRIBS run. If the copy fails, the raw folder is KEPT.
  5. Safer: it decides whether to skip a run BEFORE building it (the old scripts built first,
     which left empty folders and cost about 20 s per skipped run); it refuses to re-run runs that
     already exist unless you say --skip_existing (a re-run would cost hours); it checks free disk
     space before every run; it writes its CSVs atomically (a half-written file can never be read);
     Ctrl+C also stops the tRIBS run that is in progress.

WHAT IT DOES NOT CHANGE
------------------------
Routing stays pinned at the synthetic-truth values (kinemvelcoef 4.5, flowexp 0.24,
channelroughness 0.026). The builder (build_sensitivity_run_tribs6.py) and the scorer
(run_sensitivity_single_interp_tribs6.py) are used unchanged. optsnow stays 0. The run still needs
SYNTHETIC TRUTH MODE (exactly one *.qout in calibration_work/synth_truth/): every run is scored
against that truth by the scorer; the REAL-gauge scoring is done afterwards, on the stored
compare CSVs, by score_lowf_real_110.py (OFF or ON alone) or rescore_real_gauge_110.py v2 (both
series, paired). The KGE printed during the sweep is against the SYNTHETIC truth. It is not
the real-data fit.

USAGE (from lab/)
------------------
    python run_lowf_sweep_110.py --mode off --check_only         # checks only; builds and runs nothing
    python run_lowf_sweep_110.py --mode off --dry_run            # + prints the design and time estimate
    python run_lowf_sweep_110.py --mode off --limit 2            # first 2 samples only (a trial)
    python run_lowf_sweep_110.py --mode off --skip_existing      # the full 100; also resumes
    # other box (new label so nothing is overwritten):
    python run_lowf_sweep_110.py --mode off --label lf2 --f_lo 0.0005 --f_hi 0.004 --n 100
    # Stage B (f fixed; OFF samples Ks only in effect, ON samples Ks and cc), run each in the foreground:
    python run_lowf_sweep_110.py --mode off --label sb1 --f_fixed 0.0007 --ks_lo 5.5 --ks_hi 8.5 --n 30 --skip_existing
    python run_lowf_sweep_110.py --mode on  --label sb1 --f_fixed 0.0007 --ks_lo 3.5 --ks_hi 8.5 --n 120 --skip_existing
    # (a Codespace can stop a run if the tab is closed; the same command with --skip_existing picks up where it stopped)
    # Multi-storm window (version 1.3), a timing trial of one run, then the same command without --limit:
    python run_lowf_sweep_110.py --mode on --label sbl --f_fixed 0.0007 --ks_lo 3.5 --ks_hi 8.5 --n 120 --runtime_hours 1464 --limit 1

DO NOT run any other tRIBS build/run script in lab/ while this runs: all sweep scripts share
calibration_work/current_run_config.json with no locking. Scoring scripts that only read CSVs
(score_lowf_real_110.py, rescore_real_gauge_110.py) are safe to run during a sweep.

TIME: an estimate only, until the first run reports its own. Measured on this project's runs (build,
tRIBS and scoring together) for the default 450-hour window: 66 to 69 s per run. 30 runs is about 34 minutes;
100 runs about 1.9 hours; 120 runs about 2.3 hours; 250 runs about 4.7 hours. For a longer window the script
scales that figure in proportion to the window (an assumption, not a measurement).

OUTPUT (calibration_work/03_comparisons/summary_tables/)
    lhs_results_lowf_<label>_<OFF|ON>_110.csv          one row per finished run
    lhs_results_lowf_<label>_<OFF|ON>_FAILED_110.csv   runs that hung or failed (raw folder kept)
    lhs_design_lowf_<label>_<OFF|ON>_110.json          the settings (box, n, seed); read by the scorer
  and calibration_work/kept_qout_110/<run_id>_Outlet.qout   one per finished run
"""

import argparse
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import build_sensitivity_run_tribs6 as builder
except ImportError as exc:
    sys.exit("Could not import build_sensitivity_run_tribs6.py (%s). Run this script from the lab/ "
             "directory, where the builder lives." % exc)

# ------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------
LHS_SERIES = "110"
RUN_CATEGORY = "110_channel_conductivity"
SCORER_SCRIPT = "run_sensitivity_single_interp_tribs6.py"
KEEP_DIRNAME = "kept_qout_110"

# Routing pinned at the synthetic-truth values, exactly as in Stage 2 and the control.
ROUTING_TRUTH = {
    "kinemvelcoef":     4.5,
    "flowexp":          0.24,
    "channelroughness": 0.026,
}

DEFAULT_BOX = {
    "ks_lo": 3.0,   "ks_hi": 10.5,
    "f_lo":  0.001, "f_hi":  0.006,
    "cc_lo": 30.0,  "cc_hi": 1000.0,
}

EST_SEC_PER_RUN = 68.0        # measured on this project's runs: 66 to 69 s including the input-file build (450 h window)
MIN_FREE_GB_DEFAULT = 3.0     # one run needs about 0.3 GB while it runs; failed runs keep theirs
RAW_GB_PER_RUN = 0.3          # raw folder of one run in the 450-hour window (measured 272 MB)
BASE_TIMEOUT_SEC = 300        # per-run hard timeout for the 450-hour window
MATCH_RTOL = 1e-3             # tolerance when matching points to the other series
LABEL_RE = re.compile(r"^[A-Za-z0-9]{1,12}$")

# The builder's own window, read before any --runtime_hours override (450 in the real builder).
BASE_RUNTIME_HOURS = float(getattr(builder, "RUNTIME_HOURS", 450))
# A kept Outlet .qout must reach to within this many hours of the end of a --runtime_hours window.
COVER_TOL_H = 2.0
# Storm dates inside the 2014 forcing file (from Handoff_ChannelLossRealData_v7.md); only used to tell the user
# which of them a window takes in. Not used for any scoring.
KNOWN_STORM_DATES = [("12 Aug", "2014-08-12"), ("19 Aug", "2014-08-19"), ("8 Sep", "2014-09-08"), ("27 Sep", "2014-09-27")]


# ------------------------------------------------------------------
# Sampling (COPIED UNCHANGED from run_cc_stage2_control_110.py, so the same seed,
# n and box give the same points as that design)
# ------------------------------------------------------------------
def generate_lhs_samples(n, params, seed=None):
    """Log- or linear-stratified LHS. Each parameter is independently stratified into n bins and
    shuffled (stratified-independent sampling, not a true orthogonal LHS design, matching every
    earlier sweep script in this project). The parameters are drawn in dict order (Ks, f, cc),
    so Ks and f do not depend on whether cc is in the design."""
    rng = np.random.default_rng(seed)
    samples = {}
    for param, bounds in params.items():
        lo, hi = bounds["lo"], bounds["hi"]
        scale = bounds.get("scale", "linear")

        if scale in ("log", "fixed"):
            # a "fixed" parameter (lo == hi) goes through the same draws, so Ks and cc come out the same as
            # in a design with the same seed and n; its values are then replaced by the single fixed value
            log_lo, log_hi = np.log10(lo), np.log10(hi)
            intervals = np.linspace(log_lo, log_hi, n + 1)
            log_points = rng.uniform(intervals[:-1], intervals[1:])
            points = 10 ** log_points
        else:
            intervals = np.linspace(lo, hi, n + 1)
            points = rng.uniform(intervals[:-1], intervals[1:])

        rng.shuffle(points)
        if scale == "fixed":
            points = np.full(n, float(lo))
        samples[param] = points
    return pd.DataFrame(samples)


def make_params(args):
    f_fixed = getattr(args, "f_fixed", None)
    if f_fixed is not None:
        f_param = {"lo": f_fixed, "hi": f_fixed, "scale": "fixed"}
    else:
        f_param = {"lo": args.f_lo, "hi": args.f_hi, "scale": "log"}
    return {
        "Ks_mult":                  {"lo": args.ks_lo, "hi": args.ks_hi, "scale": "linear"},
        "f_RS_abs":                 f_param,
        "channelconductivity_mmhr": {"lo": args.cc_lo, "hi": args.cc_hi, "scale": "log"},
    }


PARAM_COLS = ["Ks_mult", "f_RS_abs", "channelconductivity_mmhr"]


# ------------------------------------------------------------------
# Names
# ------------------------------------------------------------------
def series_name(mode):
    return "ON" if mode == "on" else "OFF"


def run_tag(label, mode):
    return "%s%s" % (label, mode)


def planned_run_id(cc_val, label, mode):
    """The run_id the builder will give this sample. Computed BEFORE building, so a run that
    already exists can be skipped without paying for a build. Same recipe as
    build_sensitivity_run_tribs6.build_input_file: build_run_id(...) then '_<tag>'."""
    rid, _ = builder.build_run_id("channelconductivity_mmhr", cc_val)
    return "%s_%s" % (rid, run_tag(label, mode))


def result_paths(summary_dir, label, mode):
    ser = series_name(mode)
    return {
        "results": summary_dir / ("lhs_results_lowf_%s_%s_%s.csv" % (label, ser, LHS_SERIES)),
        "failed":  summary_dir / ("lhs_results_lowf_%s_%s_FAILED_%s.csv" % (label, ser, LHS_SERIES)),
        "design":  summary_dir / ("lhs_design_lowf_%s_%s_%s.json" % (label, ser, LHS_SERIES)),
    }


# ------------------------------------------------------------------
# Small helpers
# ------------------------------------------------------------------
def write_csv_atomic(df, path):
    tmp = path.with_name(path.name + ".part")
    df.to_csv(tmp, index=False)
    os.replace(str(tmp), str(path))


def write_json_atomic(path, obj):
    tmp = path.with_name(path.name + ".part")
    tmp.write_text(json.dumps(obj, indent=2, default=str))
    os.replace(str(tmp), str(path))


def free_gb(path):
    return shutil.disk_usage(str(path)).free / 1e9


def csv_already_exists(run_id, csv_dir):
    return (csv_dir / ("%s_compare_obs_sim.csv" % run_id)).exists()


def window_start():
    """Start of the model run (the builder's START_DATE, 1 Aug 2014 00:00)."""
    try:
        return pd.to_datetime(getattr(builder, "START_DATE", "08/01/2014/00/00"), format="%m/%d/%Y/%H/%M")
    except Exception:
        return pd.Timestamp("2014-08-01 00:00")


def event_end_hour():
    """Hours from the start of the run to the end of the 12 Aug scoring window (300 for 13 Aug 12:00)."""
    try:
        return (pd.Timestamp(getattr(builder, "EVENT_END", "2014-08-13 12:00")) - window_start()).total_seconds() / 3600.0
    except Exception:
        return 300.0


def qout_last_hour(path):
    """Time (hours from the start of the run) in the last data row of an Outlet .qout, or None if it cannot be
    read. Only the end of the file is read."""
    try:
        size = path.stat().st_size
        with open(str(path), "rb") as fh:
            fh.seek(max(0, size - 4096))
            tail = fh.read().decode("utf-8", errors="ignore")
    except OSError:
        return None
    for line in reversed(tail.splitlines()):
        parts = [p for p in re.split(r"[,\s]+", line.strip()) if p]
        if not parts:
            continue
        try:
            return float(parts[0])
        except ValueError:
            return None
    return None


def qout_coverage_problem(path, window_hours):
    """'' if the Outlet .qout reaches the end of the window (within COVER_TOL_H), else a sentence saying why not."""
    last = qout_last_hour(path)
    if last is None:
        return "could not read the last time in %s" % path.name
    if last < window_hours - COVER_TOL_H:
        return "the Outlet .qout stops at hour %.1f of the %g-hour window (tRIBS ended early)" % (last, window_hours)
    return ""


def kept_qout_ok(run_id, keep_dir, window_hours=None):
    """The kept Outlet .qout exists and is not empty; with window_hours, it must also reach the end of the window."""
    p = keep_dir / ("%s_Outlet.qout" % run_id)
    try:
        if not (p.exists() and p.stat().st_size > 0):
            return False
    except OSError:
        return False
    if window_hours is not None and qout_coverage_problem(p, window_hours):
        return False
    return True


def run_is_done(run_id, csv_dir, keep_dir, summary_dir, existing_ids, no_keep_qout, window_hours=None):
    """A run is DONE only if all three things are in place: its compare CSV, its kept Outlet .qout (unless
    --no_keep_qout), and its row in the results (or the metrics file the row is rebuilt from). Version 1.1 looked
    at the compare CSV alone, so a run cut off after that file was written (before the .qout copy) was skipped
    on resume and lost its .qout. With --runtime_hours (window_hours), the kept .qout must also reach the end of
    the window. Returns (done, reason_if_not)."""
    if not csv_already_exists(run_id, csv_dir):
        return False, "no compare CSV"
    if not no_keep_qout and not kept_qout_ok(run_id, keep_dir, window_hours):
        if window_hours is not None and kept_qout_ok(run_id, keep_dir):
            return False, "compare CSV exists but the kept Outlet .qout does not reach the end of the window"
        return False, "compare CSV exists but the Outlet .qout was not kept"
    if run_id not in existing_ids and not (summary_dir / ("%s_metrics_summary.csv" % run_id)).exists():
        return False, "compare CSV exists but there is no result row and no metrics file"
    return True, ""


def load_csv_if_any(path):
    if path.exists():
        try:
            return pd.read_csv(path)
        except Exception as exc:
            print("  Warning: could not read %s (%s)." % (path.name, exc))
    return pd.DataFrame()


def match_to_other(points_df, other_df, rtol=MATCH_RTOL):
    """Positional index of the row in other_df with the same (Ks, f, cc) within rtol, else -1."""
    out = np.full(len(points_df), -1, dtype=int)
    if other_df is None or other_df.empty or any(c not in other_df.columns for c in PARAM_COLS):
        return out
    other = other_df[PARAM_COLS].to_numpy(dtype=float)
    pts = points_df[PARAM_COLS].to_numpy(dtype=float)
    for i, row in enumerate(pts):
        ok = np.all(np.isclose(other, row, rtol=rtol, atol=0.0), axis=1)
        hits = np.flatnonzero(ok)
        if hits.size:
            out[i] = hits[0]
    return out


def design_dict(args, params):
    rt = getattr(args, "runtime_hours", None)
    return {
        "script": Path(__file__).name,
        "version": "1.3",
        "label": args.label,
        "mode": args.mode,
        "series": series_name(args.mode),
        "optpercolation": 1 if args.mode == "on" else 0,
        "n": int(args.n),
        "seed": int(args.seed),
        "f_fixed": None if getattr(args, "f_fixed", None) is None else float(args.f_fixed),
        "runtime_hours": None if rt is None else int(rt),
        "box": {k: {"lo": float(v["lo"]), "hi": float(v["hi"]), "scale": v["scale"]} for k, v in params.items()},
        "routing_pinned": ROUTING_TRUTH,
        "run_tag": run_tag(args.label, args.mode),
    }


def same_design(a, b):
    keys = ["label", "mode", "n", "seed", "box", "routing_pinned", "optpercolation", "runtime_hours"]
    return all(a.get(k) == b.get(k) for k in keys)


def same_design_but_mode(a, b):
    """True if two design files describe the same points (everything but the series): the pairing case."""
    keys = ["label", "n", "seed", "box", "routing_pinned"]
    return all(a.get(k) == b.get(k) for k in keys)


# ------------------------------------------------------------------
# Keeping the Outlet .qout, then removing the bulky raw folder
# ------------------------------------------------------------------
def find_raw_outlet_qout(run_id, calib_dir):
    """The run's Outlet .qout in its raw results folder. Returns (path or None, problem text)."""
    raw_dir = calib_dir / "02_results" / RUN_CATEGORY / run_id
    src = raw_dir / ("%s_Outlet.qout" % run_id)
    if not src.exists():
        cands = sorted(raw_dir.glob("*_Outlet.qout")) if raw_dir.exists() else []
        if len(cands) != 1:
            return None, "no single *_Outlet.qout in the raw folder (found %d)" % len(cands)
        src = cands[0]
    return src, ""


def keep_outlet_qout(run_id, calib_dir, keep_dir):
    """Copy <run_id>_Outlet.qout to keep_dir. Returns (destination or None, problem text)."""
    src, problem = find_raw_outlet_qout(run_id, calib_dir)
    if src is None:
        return None, problem
    size = src.stat().st_size
    if size == 0:
        return None, "the Outlet .qout is empty"
    try:
        keep_dir.mkdir(parents=True, exist_ok=True)
        dst = keep_dir / ("%s_Outlet.qout" % run_id)
        tmp = dst.with_name(dst.name + ".part")
        shutil.copy2(str(src), str(tmp))
        if tmp.stat().st_size != size:
            tmp.unlink()
            return None, "the copy has a different size from the original"
        os.replace(str(tmp), str(dst))
    except OSError as exc:
        return None, "copy failed: %s" % exc
    return dst, ""


def cleanup_raw_results(run_id, calib_dir):
    """Delete a finished run's raw results folder (pixel files, spatial snapshots). The .in file and
    the two CSVs live elsewhere and are untouched. Never raises."""
    raw_dir = calib_dir / "02_results" / RUN_CATEGORY / run_id
    try:
        if raw_dir.exists():
            shutil.rmtree(str(raw_dir))
            return True
    except Exception as exc:
        print("  (cleanup warning: could not remove %s: %s)" % (raw_dir, exc))
    return False


# ------------------------------------------------------------------
# Running one scored run (tRIBS + scorer) in its own process group
# ------------------------------------------------------------------
def kill_group(proc):
    try:
        os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
    except (ProcessLookupError, PermissionError):
        pass


def run_with_timeout(timeout_sec):
    """Same as the Stage 2 scripts, plus: Ctrl+C also kills the run in progress (it lives in its own
    process group, so the terminal's Ctrl+C would otherwise not reach it)."""
    proc = subprocess.Popen(
        [sys.executable, SCORER_SCRIPT],
        cwd=str(Path.cwd()),
        start_new_session=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    t0 = time.time()
    try:
        stdout, _ = proc.communicate(timeout=timeout_sec)
        timed_out = False
    except subprocess.TimeoutExpired:
        kill_group(proc)
        stdout, _ = proc.communicate()
        timed_out = True
    except KeyboardInterrupt:
        kill_group(proc)
        proc.communicate()
        raise
    elapsed = time.time() - t0

    if stdout:
        print(stdout, end="" if stdout.endswith("\n") else "\n")

    tribs_warning = bool(stdout) and "WARNING: tRIBS may have failed" in stdout
    returncode = None if timed_out else proc.returncode
    return returncode, elapsed, timed_out, tribs_warning


# ------------------------------------------------------------------
# Checks (no side effects)
# ------------------------------------------------------------------
class Checks(object):
    def __init__(self):
        self.rows = []

    def add(self, status, name, detail):
        self.rows.append((status, name, detail))

    def n(self, status):
        return sum(1 for r in self.rows if r[0] == status)

    def show(self):
        print("\nCHECKS")
        for st, name, detail in self.rows:
            print("  %-5s %-28s %s" % (st, name, detail))


def run_checks(args, calib_dir, summary_dir, csv_dir, keep_dir, samples, params, paths):
    ck = Checks()

    # synthetic truth mode
    synth_dir = calib_dir / "synth_truth"
    qouts = list(synth_dir.glob("*.qout")) if synth_dir.exists() else []
    if len(qouts) == 1:
        ck.add("PASS", "synthetic truth mode", "exactly one *.qout in synth_truth/: %s" % qouts[0].name)
    else:
        ck.add("FAIL", "synthetic truth mode",
               "expected exactly one *.qout in %s, found %d: %s. Do not run a sweep in gauge mode."
               % (synth_dir, len(qouts), [q.name for q in qouts]))

    # builder and scorer
    ok_b = hasattr(builder, "build_input_file") and hasattr(builder, "build_run_id")
    ck.add("PASS" if ok_b else "FAIL", "builder", "build_input_file and build_run_id found" if ok_b
           else "the builder in this folder lacks build_input_file or build_run_id")
    scorer = Path.cwd() / SCORER_SCRIPT
    ck.add("PASS" if scorer.exists() else "FAIL", "scorer", "%s found in %s" % (SCORER_SCRIPT, Path.cwd())
           if scorer.exists() else "%s not found in %s (run from lab/)" % (SCORER_SCRIPT, Path.cwd()))

    # box
    bad = []
    for name, p in params.items():
        if p["scale"] == "fixed":
            if not (np.isfinite(p["lo"]) and p["lo"] > 0 and p["lo"] == p["hi"]):
                bad.append("%s: a fixed value must be a positive number" % name)
        elif not (np.isfinite(p["lo"]) and np.isfinite(p["hi"]) and p["lo"] < p["hi"]):
            bad.append("%s: lo %g must be below hi %g" % (name, p["lo"], p["hi"]))
        elif p["scale"] == "log" and p["lo"] <= 0:
            bad.append("%s: a log range needs lo > 0" % name)
    f_txt = ("f fixed at %g" % args.f_fixed) if args.f_fixed is not None else ("f %g-%g (log)" % (args.f_lo, args.f_hi))
    ck.add("FAIL" if bad else "PASS", "box", "; ".join(bad) if bad else
           "Ks %g-%g, %s, cc %g-%g (log)" % (args.ks_lo, args.ks_hi, f_txt, args.cc_lo, args.cc_hi))
    f_top = args.f_fixed if args.f_fixed is not None else args.f_hi
    if f_top >= 0.0301:
        ck.add("WARN", "f range", "f above 0.030 is outside anything run so far in this series")

    # run window (version 1.3)
    win = getattr(args, "runtime_hours", None)
    scale = (win / BASE_RUNTIME_HOURS) if win else 1.0
    if win is not None:
        ev_end = event_end_hour()
        start = window_start()
        end = start + pd.Timedelta(hours=win)
        if win <= ev_end + 1:
            ck.add("FAIL", "run window", "%d h ends before the 12 Aug scoring window finishes (hour %g); use at least %d"
                   % (win, ev_end, int(ev_end) + 2))
        else:
            ck.add("PASS", "run window", "%d h: %s to %s" % (win, start.strftime("%d %b %Y %H:%M"), end.strftime("%d %b %Y %H:%M")))
            parts = []
            for lab, d in KNOWN_STORM_DATES:
                d0 = pd.Timestamp(d)
                d1 = d0 + pd.Timedelta(days=1)
                parts.append("%s %s" % (lab, "inside" if d1 <= end else ("not in" if d0 >= end else "partly in")))
            ck.add("NOTE", "storms in the window", "; ".join(parts))
            if win > 2400:
                ck.add("WARN", "run window", "%d h is longer than anything planned (about 100 days); check --runtime_hours" % win)
            ck.add("NOTE", "forcing coverage", "NOT checked here: the rainfall and weather files must run at least to %s"
                   % end.strftime("%d %b %Y %H:%M"))
            ck.add("NOTE", "time and disk",
                   "timeout %d s per run%s; free-disk limit %.1f GB; about %.0f s per run if run time grows in proportion "
                   "to the window (an assumption, the first run prints the real figure)"
                   % (args.timeout, " (scaled to the window)" if getattr(args, "timeout_scaled", False) else "",
                      args.min_free_gb, EST_SEC_PER_RUN * scale))

    # kept-qout folder must not sit inside synth_truth
    kd, sd = keep_dir.resolve(), synth_dir.resolve()
    if kd == sd or sd in kd.parents:
        ck.add("FAIL", "kept .qout folder", "%s is inside synth_truth/ and would break synthetic mode" % keep_dir)
    else:
        ck.add("PASS", "kept .qout folder", "%s%s" % (keep_dir, "" if not args.no_keep_qout else "   [OFF: --no_keep_qout]"))

    # disk
    gb = free_gb(calib_dir)
    if gb < args.min_free_gb:
        ck.add("FAIL", "free disk space", "%.1f GB free, below --min_free_gb %.1f" % (gb, args.min_free_gb))
    else:
        ck.add("PASS", "free disk space", "%.1f GB free (a run needs about %.1f GB while it runs)" % (gb, RAW_GB_PER_RUN * scale))

    # design collisions with an earlier design of the same label and mode
    design = design_dict(args, params)
    if paths["design"].exists():
        try:
            old = json.loads(paths["design"].read_text())
        except Exception:
            old = None
        if old is None:
            ck.add("FAIL", "earlier design file", "%s exists but cannot be read" % paths["design"].name)
        elif not same_design(old, design):
            ck.add("FAIL", "earlier design file",
                   "%s holds DIFFERENT settings (n, seed, box or mode). Using it would mix two designs under one "
                   "label. Pick a new --label." % paths["design"].name)
        else:
            ck.add("PASS", "earlier design file", "same settings as %s (resuming)" % paths["design"].name)
    else:
        ck.add("NOTE", "earlier design file", "none (a new design)")

    # runs that already exist
    ids = [planned_run_id(cc, args.label, args.mode) for cc in samples["channelconductivity_mmhr"]]
    if len(set(ids)) != len(ids):
        ck.add("FAIL", "unique run ids", "two samples would get the same run_id (same cc to 6 decimals)")
    else:
        ck.add("PASS", "unique run ids", "%d samples, %d different run ids, e.g. %s" % (len(ids), len(set(ids)), ids[0]))
    existing = sum(1 for rid in ids if csv_already_exists(rid, csv_dir))
    if existing and not args.skip_existing:
        ck.add("FAIL", "runs already exist",
               "%d of %d runs of this design already have a compare CSV. Add --skip_existing to resume, or pick a "
               "new --label. (Without it they would all be re-run and overwritten.)" % (existing, len(ids)))
    elif existing:
        old_df = load_csv_if_any(paths["results"])
        old_ids = set(old_df["run_id"].values) if not old_df.empty else set()
        done = sum(1 for rid in ids if run_is_done(rid, csv_dir, keep_dir, summary_dir, old_ids, args.no_keep_qout, win)[0])
        redo = existing - done
        txt = "%d of %d already done; --skip_existing will skip them" % (done, len(ids))
        if redo:
            txt += ("; %d more have a compare CSV but are incomplete (no kept Outlet .qout or no result row) and "
                    "will be re-run" % redo)
        ck.add("NOTE", "runs already exist", txt)
    else:
        ck.add("PASS", "runs already exist", "none of the %d runs exist yet" % len(ids))
    return ck, ids, design


def print_design(samples, ids, args, n_run):
    print("\nDESIGN   (label %s, series %s, n %d, seed %d)" % (args.label, series_name(args.mode), args.n, args.seed))
    print("  %5s  %9s  %10s  %9s   %s" % ("#", "Ks_mult", "f_RS_abs", "cc", "run_id"))
    show = min(len(samples), 12)
    for i in range(show):
        r = samples.iloc[i]
        print("  %5d  %9.3f  %10.5f  %9.1f   %s" % (i + 1, r["Ks_mult"], r["f_RS_abs"],
                                                   r["channelconductivity_mmhr"], ids[i]))
    if len(samples) > show:
        print("  ... (%d more)" % (len(samples) - show))
    for c in PARAM_COLS:
        print("  %-26s min %-10.5g max %-10.5g" % (c, samples[c].min(), samples[c].max()))
    if args.mode == "off":
        print("  (series OFF: cc is sampled and written into the .in file, but with OPTPERCOLATION=0 it has no effect)")
    if args.f_fixed is not None:
        print("  (f is FIXED at %g for every run)" % args.f_fixed)
    rt = getattr(args, "runtime_hours", None)
    est = EST_SEC_PER_RUN * ((rt / BASE_RUNTIME_HOURS) if rt else 1.0)
    hours = n_run * est / 3600.0
    print("\n  This invocation would run %d of %d samples: about %.1f hours at %.0f s per run (an estimate%s; the "
          "first run prints the real figure)." % (n_run, args.n, hours, est,
                                                  ", scaled to the window in proportion: unmeasured" if rt else ""))


# ------------------------------------------------------------------
def main():
    # Print as it happens even when the output goes through "| tee" or into a log file (nohup ... > log).
    # Python otherwise holds printed text back in an 8 KB buffer when stdout is not a terminal, so the
    # screen looks frozen and the last chunk is lost if the process is killed.
    for _stream in (sys.stdout, sys.stderr):
        try:
            _stream.reconfigure(line_buffering=True)
        except (AttributeError, ValueError):
            pass
    parser = argparse.ArgumentParser(
        description="Series 110 low-f sweep: LHS over Ks_mult / f_RS_abs / channelconductivity_mmhr with the box "
                    "on the command line; channel loss OFF or ON; keeps each run's Outlet .qout.")
    parser.add_argument("--mode", choices=["off", "on"], required=True,
                        help="off = optpercolation 0 (cc inert); on = optpercolation 1 (channel loss on)")
    parser.add_argument("--label", type=str, default="lf1",
                        help="Name of this design, letters and digits only (default lf1). A new box needs a new "
                             "label so nothing is overwritten.")
    parser.add_argument("--n", type=int, default=100, help="Number of LHS samples (default 100)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed (default 42)")
    parser.add_argument("--ks_lo", type=float, default=DEFAULT_BOX["ks_lo"])
    parser.add_argument("--ks_hi", type=float, default=DEFAULT_BOX["ks_hi"])
    parser.add_argument("--f_lo", type=float, default=None,
                        help="lowest f_RS_abs (default %g); not with --f_fixed" % DEFAULT_BOX["f_lo"])
    parser.add_argument("--f_hi", type=float, default=None,
                        help="highest f_RS_abs (default %g); not with --f_fixed" % DEFAULT_BOX["f_hi"])
    parser.add_argument("--f_fixed", type=float, default=None,
                        help="hold f_RS_abs at this one value for every run (Stage B uses 0.0007) instead of sampling it")
    parser.add_argument("--cc_lo", type=float, default=DEFAULT_BOX["cc_lo"])
    parser.add_argument("--cc_hi", type=float, default=DEFAULT_BOX["cc_hi"])
    parser.add_argument("--limit", type=int, default=None,
                        help="Run only the first K samples of the design (a trial); the design itself is unchanged")
    parser.add_argument("--skip_existing", action="store_true",
                        help="Skip samples whose compare CSV already exists (resume). Without it, the script "
                             "refuses to start if any run of the design already exists.")
    parser.add_argument("--runtime_hours", type=int, default=None,
                        help="Run every sample over this many hours from 1 Aug 2014 00:00 instead of the builder's "
                             "450 (1464 reaches 1 Oct 2014 and takes in 12 Aug, 19 Aug, 8 Sep and 27 Sep). The Outlet "
                             ".qout of the whole window is kept; scoring stays on the 12 Aug window.")
    parser.add_argument("--timeout", type=int, default=None,
                        help="Per-run hard timeout in seconds (default 300; with --runtime_hours, 300 scaled to the window)")
    parser.add_argument("--no_cleanup", action="store_true",
                        help="Keep every raw per-run folder (about 270 MB each) instead of deleting it after the run")
    parser.add_argument("--no_keep_qout", action="store_true",
                        help="Do not copy the Outlet .qout before cleanup (not recommended)")
    parser.add_argument("--min_free_gb", type=float, default=None,
                        help="Stop if free disk space falls below this many GB (default %g; with --runtime_hours, four "
                             "runs' worth of raw output if that is more)" % MIN_FREE_GB_DEFAULT)
    parser.add_argument("--check_only", action="store_true", help="Run the checks and stop; builds and runs nothing")
    parser.add_argument("--dry_run", action="store_true",
                        help="Run the checks, print the design and the time estimate, and stop; writes nothing")
    args = parser.parse_args()

    if args.f_fixed is not None and (args.f_lo is not None or args.f_hi is not None):
        parser.error("use either --f_fixed (one value) or --f_lo/--f_hi (a range), not both.")
    if args.f_fixed is not None and not (args.f_fixed > 0):
        parser.error("--f_fixed must be a positive number.")
    if args.f_lo is None:
        args.f_lo = DEFAULT_BOX["f_lo"]
    if args.f_hi is None:
        args.f_hi = DEFAULT_BOX["f_hi"]
    if not LABEL_RE.match(args.label):
        parser.error("--label must be 1 to 12 letters or digits (no underscores or spaces).")
    if args.n < 5:
        parser.error("--n must be at least 5.")
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be at least 1.")
    if args.runtime_hours is not None and args.runtime_hours < 1:
        parser.error("--runtime_hours must be a positive whole number of hours.")
    # window scaling (version 1.3): without --runtime_hours these resolve to the old defaults (300 s, 3.0 GB)
    scale = (args.runtime_hours / BASE_RUNTIME_HOURS) if args.runtime_hours else 1.0
    args.timeout_scaled = False
    if args.timeout is None:
        args.timeout = int(round(BASE_TIMEOUT_SEC * max(1.0, scale)))
        args.timeout_scaled = bool(args.runtime_hours) and scale > 1.0
    if args.min_free_gb is None:
        args.min_free_gb = max(MIN_FREE_GB_DEFAULT, round(4.0 * RAW_GB_PER_RUN * scale, 1)) if args.runtime_hours \
            else MIN_FREE_GB_DEFAULT
    if args.timeout < 30:
        parser.error("--timeout must be at least 30 seconds.")
    if args.runtime_hours is not None:
        builder.RUNTIME_HOURS = int(args.runtime_hours)    # the builder reads this when it writes each input file

    script_dir = Path.cwd()
    project_root = script_dir.parent
    calib_dir = project_root / "calibration_work"
    if not calib_dir.exists():
        sys.exit("Could not find %s. Run this script from the lab/ directory (one level below the project root)."
                 % calib_dir)
    summary_dir = calib_dir / "03_comparisons" / "summary_tables"
    csv_dir = calib_dir / "03_comparisons" / "csv_exports"
    keep_dir = calib_dir / KEEP_DIRNAME

    params = make_params(args)
    # an invalid box is reported by the checks below; sampling needs a valid one
    box_ok = all(np.isfinite(p["lo"]) and np.isfinite(p["hi"])
                 and (p["lo"] == p["hi"] and p["lo"] > 0 if p["scale"] == "fixed" else p["lo"] < p["hi"])
                 and (p["scale"] != "log" or p["lo"] > 0) for p in params.values())
    if not box_ok:
        print("Box is not valid:")
        for name, p in params.items():
            print("  %s: %g - %g (%s)" % (name, p["lo"], p["hi"], p["scale"]))
        sys.exit("Each range needs lo < hi, and log ranges need lo > 0.")
    # (a fixed value is a range with lo = hi)

    samples = generate_lhs_samples(args.n, params, seed=args.seed)
    paths = result_paths(summary_dir, args.label, args.mode)

    print("\n" + "=" * 78)
    print("Series %s -- LOW-f sweep (version 1.3) -- label %s, series %s"
          % (LHS_SERIES, args.label, series_name(args.mode)))
    print("python %s, pandas %s, numpy %s" % (sys.version.split()[0], pd.__version__, np.__version__))
    print("  %d samples, seed %d, timeout %d s, cleanup %s, keep Outlet .qout %s"
          % (args.n, args.seed, args.timeout, "off" if args.no_cleanup else "on", "off" if args.no_keep_qout else "on"))
    if args.runtime_hours is not None:
        print("  RUN WINDOW: %d hours from %s (builder default %g); every Outlet .qout must reach the end of it"
              % (args.runtime_hours, window_start().strftime("%d %b %Y %H:%M"), BASE_RUNTIME_HOURS))
    for pname, p in params.items():
        if p["scale"] == "fixed":
            print("  %-26s %g  (FIXED for every run)" % (pname, p["lo"]))
        else:
            print("  %-26s %g - %g  (%s)" % (pname, p["lo"], p["hi"], p["scale"]))
    print("  PINNED (truth routing): kinemvelcoef=%g  flowexp=%g  channelroughness=%g"
          % (ROUTING_TRUTH["kinemvelcoef"], ROUTING_TRUTH["flowexp"], ROUTING_TRUTH["channelroughness"]))
    print("  optpercolation=%d for every run%s" % (1 if args.mode == "on" else 0,
                                                    "  (cc is sampled but inert)" if args.mode == "off" else ""))
    print("  The KGE printed during the sweep is against the SYNTHETIC truth, not the real gauge.")
    print("=" * 78)

    ck, ids, design = run_checks(args, calib_dir, summary_dir, csv_dir, keep_dir, samples, params, paths)
    ck.show()

    n_run = args.n if args.limit is None else min(args.limit, args.n)
    if args.dry_run:
        print_design(samples, ids, args, n_run)
    if ck.n("FAIL"):
        sys.exit("\nSTOPPED: a check failed (see FAIL above). Nothing was built or run.")
    if args.check_only or args.dry_run:
        print("\n%s: checks passed; nothing was built, run or written." % ("--check_only" if args.check_only else "--dry_run"))
        return

    # ------------------------------------------------------------------
    # pairing check with the other series of the same label (read-only)
    # ------------------------------------------------------------------
    other_mode = "on" if args.mode == "off" else "off"
    other_paths = result_paths(summary_dir, args.label, other_mode)
    other_df = load_csv_if_any(other_paths["results"])
    other_design = None
    if other_paths["design"].exists():
        try:
            other_design = json.loads(other_paths["design"].read_text())
        except Exception:
            other_design = None
    if other_design is not None and not same_design_but_mode(other_design, design):
        print("\nPAIRING: the %s series of label %s uses a different design (n, seed or box), so no point-for-point "
              "pairing is expected. That is intended for Stage B." % (series_name(other_mode), args.label))
    elif not other_df.empty:
        tw = match_to_other(samples, other_df)
        n_m = int((tw >= 0).sum())
        print("\nPAIRING CHECK: %d of %d samples match a point of the %s series (%d rows there)."
              % (n_m, args.n, series_name(other_mode), len(other_df)))
        if n_m < 0.9 * min(args.n, len(other_df)):
            print("  *** WARNING: poor pairing. The two series should use the same --label, --n, --seed and box. "
                  "Press Ctrl+C now if that was not intended. ***")

    summary_dir.mkdir(parents=True, exist_ok=True)
    write_json_atomic(paths["design"], dict(design, created_local=datetime.now().isoformat(timespec="seconds"),
                                            python=sys.version.split()[0], pandas=pd.__version__,
                                            numpy=np.__version__))

    existing_df = load_csv_if_any(paths["results"])
    existing_ids = set(existing_df["run_id"].values) if not existing_df.empty else set()
    results = existing_df.to_dict("records") if not existing_df.empty else []
    failed_rows = []

    completed = skipped = hung = failed = cleaned = keep_failed = 0
    sweep_start = time.time()
    run_times = []
    interrupted = False
    stopped_for_disk = False

    try:
        for i in range(n_run):
            row = samples.iloc[i]
            ks_val = float(row["Ks_mult"])
            f_val = float(row["f_RS_abs"])
            cc_val = float(row["channelconductivity_mmhr"])
            run_id = ids[i]
            print("\n[%4d/%d]  Ks=%.3fx  f=%.5f  cc=%.1f mm/hr%s  -> %s"
                  % (i + 1, n_run, ks_val, f_val, cc_val, " (inert)" if args.mode == "off" else "", run_id))

            if args.skip_existing:
                done, why_not = run_is_done(run_id, csv_dir, keep_dir, summary_dir, existing_ids, args.no_keep_qout,
                                            args.runtime_hours)
                if done:
                    print("  SKIP (compare CSV exists)")
                    skipped += 1
                    metrics_file = summary_dir / ("%s_metrics_summary.csv" % run_id)
                    if run_id not in existing_ids:
                        try:
                            m = pd.read_csv(metrics_file).iloc[0].to_dict()
                            m.setdefault("Ks_mult", ks_val)
                            m.setdefault("f_RS_abs", f_val)
                            m.setdefault("channelconductivity_mmhr", cc_val)
                            m["optpercolation"] = 1 if args.mode == "on" else 0
                            m["design_label"] = args.label
                            m["qout_kept"] = bool(kept_qout_ok(run_id, keep_dir, args.runtime_hours))
                            results.append(m)
                            existing_ids.add(run_id)
                        except Exception as exc:
                            print("  (could not read %s: %s)" % (metrics_file.name, exc))
                    # a finished run whose raw folder was left behind (cut off after the copy): tidy it up
                    if (not args.no_cleanup and not args.no_keep_qout and kept_qout_ok(run_id, keep_dir, args.runtime_hours)
                            and cleanup_raw_results(run_id, calib_dir)):
                        cleaned += 1
                        print("  (removed the raw folder left behind by an earlier cut-off run)")
                    continue
                if why_not != "no compare CSV":
                    print("  RE-RUN: %s (an earlier run was cut off, or the copy failed)." % why_not)
                    cleanup_raw_results(run_id, calib_dir)     # start the re-run from a clean folder

            gb = free_gb(calib_dir)
            if gb < args.min_free_gb:
                print("  STOPPING: %.1f GB free, below --min_free_gb %.1f. Free some space (failed runs keep about "
                      "%.0f MB each in 02_results/%s/), then resume with --skip_existing."
                      % (gb, args.min_free_gb, 270 * ((args.runtime_hours / BASE_RUNTIME_HOURS) if args.runtime_hours else 1.0),
                         RUN_CATEGORY))
                stopped_for_disk = True
                break

            try:
                built_id, _, _ = builder.build_input_file(
                    "channelconductivity_mmhr", cc_val,
                    overrides={
                        "Ks_mult":          ks_val,
                        "f_RS_abs":         f_val,
                        "kinemvelcoef":     ROUTING_TRUTH["kinemvelcoef"],
                        "flowexp":          ROUTING_TRUTH["flowexp"],
                        "channelroughness": ROUTING_TRUTH["channelroughness"],
                        "optpercolation":   1 if args.mode == "on" else 0,
                    },
                    tag=run_tag(args.label, args.mode),
                )
            except Exception as exc:
                print("  BUILD FAILED: %s" % exc)
                failed += 1
                failed_rows.append({"run_id": run_id, "status": "BUILD_FAILED", "reason": str(exc)[:200],
                                    "elapsed_min": 0.0, "Ks_mult": ks_val, "f_RS_abs": f_val,
                                    "channelconductivity_mmhr": cc_val})
                write_csv_atomic(pd.DataFrame(failed_rows), paths["failed"])
                continue
            if built_id != run_id:
                print("  STOPPING: the builder named this run %s but %s was expected. The run-id recipe in this script "
                      "no longer matches the builder." % (built_id, run_id))
                sys.exit(2)

            returncode, elapsed, timed_out, tribs_warning = run_with_timeout(args.timeout)
            metrics_file = summary_dir / ("%s_metrics_summary.csv" % run_id)

            if timed_out:
                status, reason = "HANG", "wall-clock timeout"
            elif returncode != 0 or tribs_warning:
                status = "FAILED"
                reason = "tRIBS reported non-zero exit" if tribs_warning else "scorer exited %s" % returncode
            elif not metrics_file.exists():
                status, reason = "FAILED", "scorer finished but wrote no metrics file"
            else:
                status, reason = "SUCCESS", ""
                if args.runtime_hours is not None:
                    # a long window: the Outlet .qout must reach the end of it, or the later storms are silently lost
                    raw_q, raw_problem = find_raw_outlet_qout(run_id, calib_dir)
                    short = raw_problem if raw_q is None else qout_coverage_problem(raw_q, args.runtime_hours)
                    if short:
                        status, reason = "FAILED", short

            if status == "SUCCESS":
                metrics = pd.read_csv(metrics_file).iloc[0].to_dict()
                metrics.setdefault("Ks_mult", ks_val)
                metrics.setdefault("f_RS_abs", f_val)
                metrics.setdefault("channelconductivity_mmhr", cc_val)
                metrics["optpercolation"] = 1 if args.mode == "on" else 0
                metrics["design_label"] = args.label

                kept_ok = False
                if not args.no_keep_qout:
                    dst, why = keep_outlet_qout(run_id, calib_dir, keep_dir)
                    kept_ok = dst is not None
                    if kept_ok:
                        print("  kept: %s (%.0f kB)" % (dst.name, dst.stat().st_size / 1e3))
                    else:
                        keep_failed += 1
                        print("  WARNING: could not keep the Outlet .qout (%s). The raw folder is KEPT for this run." % why)
                metrics["qout_kept"] = bool(kept_ok)

                results = [r for r in results if r.get("run_id") != run_id]
                results.append(metrics)
                completed += 1
                run_times.append(elapsed)
                print("  vs SYNTHETIC truth: KGE_2012=%.3f  PBIAS=%+.1f%%"
                      % (float(metrics.get("kge_2012", float("nan"))), float(metrics.get("pbias_pct", float("nan")))))

                if not args.no_cleanup and (kept_ok or args.no_keep_qout):
                    if cleanup_raw_results(run_id, calib_dir):
                        cleaned += 1
            else:
                print("  %s: %s  (%s)" % (status, run_id, reason))
                if status == "HANG":
                    hung += 1
                else:
                    failed += 1
                failed_rows.append({"run_id": run_id, "status": status, "reason": reason,
                                    "elapsed_min": elapsed / 60.0, "Ks_mult": ks_val, "f_RS_abs": f_val,
                                    "channelconductivity_mmhr": cc_val})
                write_csv_atomic(pd.DataFrame(failed_rows), paths["failed"])

            todo_left = n_run - (i + 1)
            if run_times:
                per = (time.time() - sweep_start) / max(completed + hung + failed, 1)
                eta_h = per * todo_left / 3600.0
                print("  scorer time %.1f min  |  average %.0f s per run so far (incl. build)  |  about %.1f h left"
                      % (elapsed / 60.0, per, eta_h))

            if results:
                write_csv_atomic(pd.DataFrame(results), paths["results"])
    except KeyboardInterrupt:
        interrupted = True
        print("\nInterrupted (Ctrl+C). The run in progress was stopped; finished runs are saved.")

    if results:
        write_csv_atomic(pd.DataFrame(results), paths["results"])

    print("\nLow-f sweep %s: %d ran, %d skipped, %d hung, %d failed, %d raw folders cleaned, %d Outlet .qout not kept"
          % ("INTERRUPTED" if interrupted else ("STOPPED (disk)" if stopped_for_disk else "finished"),
             completed, skipped, hung, failed, cleaned, keep_failed))
    if results:
        final_df = pd.DataFrame(results)
        print("Saved: %s  (%d rows)" % (paths["results"].name, len(final_df)))
        for c in PARAM_COLS:
            if c in final_df.columns:
                print("  %-26s covered %.5g - %.5g" % (c, final_df[c].min(), final_df[c].max()))
        if not args.no_keep_qout and keep_dir.exists():
            files = list(keep_dir.glob("*_Outlet.qout"))
            tot = sum(f.stat().st_size for f in files) / 1e6
            print("Kept Outlet .qout files: %d in %s (%.1f MB)" % (len(files), keep_dir, tot))
        if args.runtime_hours is not None:
            print("Next: the kept Outlet .qout files cover the whole %d-hour window. The runs are scored on 12 Aug only "
                  "so far; scoring the other storms from the kept files needs a storm-by-storm scorer (not written yet)."
                  % args.runtime_hours)
        elif args.f_fixed is not None:
            print("Next: when BOTH series of label %s are finished (ON and OFF):  python score_stageb_110.py --label %s"
                  "   (scores them against the REAL gauge and compares ON with OFF)" % (args.label, args.label))
        else:
            print("Next: python score_lowf_real_110.py --label %s --series %s   (scores these runs against the REAL gauge)"
                  % (args.label, series_name(args.mode)))
    else:
        print("No results to save.")
    if hung or failed:
        print("Failed or hung runs are listed in %s; their raw folders were left for debugging." % paths["failed"].name)
    if interrupted:
        sys.exit(130)
    if stopped_for_disk:
        sys.exit(3)


if __name__ == "__main__":
    main()
