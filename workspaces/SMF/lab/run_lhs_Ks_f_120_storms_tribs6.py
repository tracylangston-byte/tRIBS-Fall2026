"""
run_lhs_Ks_f_120_storms_tribs6.py
=================================
Series 120, legs 2 and 3 -- the 0.8x and 1.25x STORM legs of the Ks_mult x f_RS_abs
synthetic-truth experiment, on tRIBS 6.0.0 / pytRIBS 1.0.0.

WHAT THIS IS
------------
The 1.0x leg (run_lhs_Ks_f_120_tribs6.py) is finished: 400 draws of Ks_mult and
f_RS_abs scored against the synthetic truth (Ks 7.0, f 0.012) for the original
12 Aug 2014 storm. This script repeats that experiment for two other storm
sizes. The storm is made 20% smaller (0.8x) or 25% larger (1.25x) by multiplying
the RAINFALL RATE in both rain-gauge files. Nothing else about the storm changes
(same dates, same timing, same offset between the two gauges).

For each storm the experiment needs its own "truth": the hydrograph tRIBS makes
at Ks 7.0 / f 0.012 under THAT storm's rain. The sweep is then scored against it.
So each leg has two jobs, done by two commands:

  1. --make_truth   scale the rain, build the storm's truth, check it, prove it
                    reproduces itself.                          (about 3 minutes)
  2. (default)      the 400-draw sweep, scored against that truth. (about 6 hours)

Only the rain changes. The truth parameters stay Ks 7.0 / f 0.012, routing stays
at truth (kinemvelcoef 4.5, flowexp r 0.24, channelroughness 0.026),
optpercolation stays 0 (channel loss OFF) and OPTSNOW stays 0 (set inside the
builder). That keeps the comparison between storms clean: the only thing that
differs from the 1.0x leg is the storm.

The sweep uses the SAME 400 draws as the 1.0x leg (seed 42, Ks 3-11 linear, f
0.003-0.05 log), so every draw has a twin in the 1.0x results. The old 5.3.1 storm
legs used n=200 and a narrower box (Ks 3.0-9.5 / 9.5, f 0.004-0.05), so they do
NOT pair point-for-point with these; old-vs-new comparison for the storm legs has
to be made on the fitted surface / ridge, not run by run.

HOW THE RUNS ARE NAMED AND WHERE THEY GO
----------------------------------------
(label = storm080 or storm125, short = st080 or st125)
    sweep run_id      SMF_20140812_120_Ks<Ks>x_f<f>_<short>
    truth run_id      SMF_20140812_120_Ks7p0x_f0p012_truth_<short>
    verify run_id     SMF_20140812_120_Ks7p0x_f0p012_verify_<short>
    raw results       calibration_work/02_results/120_lhs_synth_Ks_f_<short>/<run_id>/   (deleted after scoring)
    run inputs        calibration_work/01_run_inputs/120_lhs_synth_Ks_f_<short>/
    scaled rain       ../init_data/met/precip_<station>_<label>.mdf, Master_Precip_<label>.sdf
    the truth         calibration_work/synth_truth_120_<label>/<truth run_id>_Outlet.qout
    its provenance    calibration_work/synth_truth_120_<label>/PROVENANCE_truth_120_<label>.json
    forcing record    calibration_work/forcing_120_<label>_PROVENANCE.json
    results CSV       calibration_work/03_comparisons/summary_tables/lhs_results_synth_Ks_f_120_<label>.csv
    failed log        calibration_work/03_comparisons/summary_tables/lhs_results_synth_Ks_f_FAILED_120_<label>.csv

If scaled rain files with these names already exist (precip_*_storm080.mdf and
Master_Precip_storm080.sdf, the same for storm125), for example from an earlier script,
and they are identical to what this script would write, it simply reuses them ("already
present (identical)"). If one exists and DIFFERS, it stops and tells you; it never
overwrites a rain file unless you give --overwrite_forcing.

RAIN FILE FORMAT: the real files are comma-separated with the header
Year,Month,Day,Hour,R_mm (older files: Y M D H R, spaces). The R column is a rate in
mm/hr, 4 rows per hour. Only that column is multiplied; every other column and every
separator is copied through as text, so the scaled file has the same layout as the original.

Nothing in calibration_work/synth_truth/ is ever written or moved. The 1.0x truth
there is only READ (as the reference for the direction check below).

THE MODES
---------
--dry_run       Writes nothing, builds nothing, runs no tRIBS. Prints the rain
                scaling plan, what already exists, the truth status, the pairing
                check against the 1.0x results, and what the other commands would
                do. Takes seconds. Do this first.
--make_truth    1. writes (or reuses) the scaled rain files and checks them against
                   the originals;
                2. builds and runs ONE tRIBS run at the truth point (Ks 7.0 / f 0.012)
                   under the scaled rain, scored against the 1.0x truth;
                3. DIRECTION CHECK: a 0.8x storm must give a smaller peak AND volume
                   than the 1.0x truth, a 1.25x storm a bigger one. If not, tRIBS did
                   not really use the scaled rain and the script stops here;
                4. saves the outlet file as this storm's truth (copy checked by
                   checksum) and sanity-checks it (columns, time order, no negative or
                   missing flow, not dry);
                5. SELF-CONSISTENCY CHECK: runs the truth point a second time, scored
                   against the new truth file via the scorer's truth_file setting.
                   It must come back with PBIAS 0.0000% and KGE_2012 1.0000
                   (|PBIAS| < 0.01% and KGE_2012 > 0.9999). This proves the truth
                   file, the scorer and the truth_file setting all agree;
                6. writes the provenance file (checksums of the truth, the rain files
                   and the 1.0x truth; the check results);
                7. if Series 110 made a cc-OFF reference for this storm, compares the
                   new truth with it (information only; expected to be identical,
                   because it is the same build and the same parameters).
                Run it again on an existing truth: it keeps the truth, re-checks the
                checksums and repeats step 5. --remake_truth renames the old truth to
                *.superseded_<time> (never deletes it) and builds a new one.
--smoke_test    Only step 5 above, on an existing truth (about 1 minute).
(default)       The 400-draw sweep. Before the first run it re-checks that the truth
                file, the scaled rain files and the provenance record still agree.

SAFETY CHECKS IN THE SWEEP
--------------------------
  * Refuses to start unless the provenance file says the self-consistency check
    PASSED, the truth file's checksum matches the record, and the scaled rain files
    are exactly what this script would write from the current originals.
  * Every run is built with the scaled SDF (the builder's gauge_sdf argument). The
    builder's run configuration is read back and must name the scaled SDF, the Ks and
    f asked for, the pinned routing, optpercolation 0; the .in file must name the
    scaled SDF. Any mismatch STOPS the sweep.
  * Every run is scored against THIS storm's truth: the script writes truth_file into
    the run configuration after each build, and after each run checks (a) the scorer's
    own "Reading observed from:" line names this storm's truth, (b) the observed peak
    it used equals the truth's peak, (c) obs_mode is "synth". Any mismatch STOPS the
    sweep at once (every later run would be scored wrongly).
  * Each run is a subprocess with a hard timeout (default 300 s); a hang kills the
    whole process group. HANG / FAILED runs go to the FAILED csv and the sweep moves
    on. 8 failures in a row stop the script (that is a systematic problem).
  * Ctrl+C also kills the in-flight run, so no orphan tRIBS keeps writing the config.
  * Refuses to start if another tRIBS build/run script seems to be running (they share
    calibration_work/current_run_config.json with no locking).
  * The raw results folder is deleted only if it sits exactly at
    02_results/<category>/<run_id>.

DO NOT run any other tRIBS build/run script in lab/ while this is running.

USAGE (run from the lab/ directory, like every other script)
-------------------------------------------------------------
  0.8x storm:
    python run_lhs_Ks_f_120_storms_tribs6.py --rain_scale 0.8 --dry_run
    python -u run_lhs_Ks_f_120_storms_tribs6.py --rain_scale 0.8 --make_truth 2>&1 | tee series120_storm080_truth.log
    python -u run_lhs_Ks_f_120_storms_tribs6.py --rain_scale 0.8 --limit 3 2>&1 | tee series120_storm080_trial.log
    python -u run_lhs_Ks_f_120_storms_tribs6.py --rain_scale 0.8 --skip_existing 2>&1 | tee series120_storm080.log

  1.25x storm: the same four commands with --rain_scale 1.25 (and storm125 in the log names).

  After any interruption (Codespace stopped, tab closed, Ctrl+C) run the same
  sweep command again: --skip_existing resumes where it left off. Keep the terminal
  open and visible while it runs; typing and output count as activity for the
  Codespace idle timer, a background job does not. The -u matters: without it Python
  holds output back when it is piped to tee and the log looks stuck.

WHAT TO LOOK AT AFTER A SWEEP
-----------------------------
The end-of-run report prints the feasible count (|PBIAS| < 2%) against the extension
threshold of 20; the f-band table (runs on each side of zero bias); the ENCLOSURE
table (is the zero-bias ridge inside the sampled Ks range at every f, or does it leave
through the top or bottom); the Ks-upper-bound / lower-bound check; and a pairing with
the 1.0x results (same draws). The best run is printed for reference only: on an
equifinal ridge the best run wanders along the ridge, so it is not the estimate.
Extension rule (agreed before the run): extend only with a NEW-SEED batch, never a
larger n; add one if there are fewer than 20 feasible runs, or if any f band has no
runs on both sides of zero bias, or if the ridge is not enclosed. The last of those
is a finding as much as a gap, and acting on it is Tracy's decision.

STATUS (2026-10-25)
-------------------
Written without access to the Codespace. NOT YET RUN against the real tRIBS 6.0.0 build.
The first real checks are --dry_run (does the real SDF parse, what already exists) and
--make_truth (does the real builder honor gauge_sdf, does the real scorer honor
truth_file, does the direction check pass); both stop at once, with a plain message, if
something is not as assumed.

Tested against stand-ins for the builder and scorer (a toy hydrograph model; rain files in
the real layout Year,Month,Day,Hour,R_mm, comma-separated, named precip_SMF_2014-2014.mdf and
precip_SMPHQ_2014-2014.mdf; the same configuration keys the Series 110 script relies on).
The stand-ins are NOT tRIBS, and the real Master_Precip.sdf layout was not seen, only that
each station line contains the path of its .mdf file:
  * --dry_run writes nothing; it reports the rain plan, what exists, the truth state, and a
    competing job (without refusing).
  * --make_truth for 0.8x and 1.25x: rain written and checked, truth built, direction check,
    truth saved by verified copy, self-consistency run (PBIAS 0.0000%, KGE_2012 1.0000),
    provenance written. Run again on an existing truth: keeps it, re-checks it. --remake_truth:
    old files renamed *.superseded_<time>, new ones built.
  * Series 110 comparison: byte-identical / same numbers / different, all reported correctly.
  * Stops by itself, with nothing usable left behind, when: tRIBS ignores the scaled rain
    (direction check), the builder ignores gauge_sdf, the scorer ignores truth_file (in the
    self-check and in a sweep), obs_mode is not synth, a scored row does not carry the asked-for
    parameters, the truth file or a scaled rain file is changed after the truth was made, the
    provenance record is missing or not verified, the draws do not pair with the 1.0x results,
    another tRIBS script is running, 8 runs fail in a row.
  * Hung and failed runs are logged and skipped; a later resume re-runs them and clears the log.
  * A full 400-draw sweep killed abruptly (kill -9) part-way and resumed with --skip_existing:
    400 rows, 400 unique run ids; a planted half-finished draw (compare CSV but no metrics) was
    re-run, not skipped.
  * Every number in the end-of-run report (feasible count, f bands, enclosure table, best run,
    pairing counts, shift, rank correlation) was reproduced independently from the CSV, and the
    enclosure table matched an analytic ridge in every slice where the answer is clear.
  * The two legs coexist: separate run ids (..._st080 / ..._st125), categories, truth folders and
    CSVs; the 1.0x results file is only read.
"""

import argparse
import hashlib
import inspect
import json
import math
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

import build_sensitivity_run_tribs6 as builder

# ----------------------------------------------------------------------
# Series identity
# ----------------------------------------------------------------------
SERIES        = "120"
CATEGORY_BASE = "120_lhs_synth_Ks_f"          # the storm tag is added: ..._st080
SCORER_SCRIPT = "run_sensitivity_single_interp_tribs6.py"

# The two storm legs this script is for.
LEG_SCALES = (0.8, 1.25)

# Peak of the storm truth divided by the peak of the 1.0x truth on the OLD build
# (tRIBS 5.3.1: 56.877 -> 29.922 and 96.147 m3/s). Printed for comparison only; the
# direction check does not use them.
OLD_BUILD_PEAK_RATIO = {0.8: 0.526, 1.25: 1.690}

# ----------------------------------------------------------------------
# Truth and pinned values (same as the 1.0x leg)
# ----------------------------------------------------------------------
TRUTH_KS = 7.0
TRUTH_F  = 0.012
TRUTH_TOL = 1e-6

ROUTING_TRUTH = {
    "kinemvelcoef":     4.5,
    "flowexp":          0.24,
    "channelroughness": 0.026,
}

# Sampling box. MUST match the 1.0x leg exactly (ranges, scales, AND dict order --
# the order sets how the random generator is consumed) or the draws will not pair.
LHS_PARAMS = {
    "Ks_mult":  {"lo": 3.0,   "hi": 11.0, "scale": "linear"},
    "f_RS_abs": {"lo": 0.003, "hi": 0.05, "scale": "log"},
}
PARAM_COLS   = list(LHS_PARAMS.keys())
DEFAULT_N    = 400
DEFAULT_SEED = 42

# ----------------------------------------------------------------------
# Analysis thresholds (same as the 1.0x leg)
# ----------------------------------------------------------------------
FEASIBLE_PBIAS_PCT     = 2.0
MIN_FEASIBLE_PER_STORM = 20
N_F_BANDS              = 6        # the f-band bracket table (same as the 1.0x report)
N_F_SLICES_ENCL        = 10       # finer slices for the ridge-enclosure check
MIN_RUNS_PER_SLICE     = 8
MATCH_RTOL             = 1e-3
CAUTION_KS_LO = 8.85
CAUTION_F_LO  = 0.0437

SERIES120_CSV_NAME = "lhs_results_synth_Ks_f_120.csv"
FAILED_COLS = ["run_id", "status", "reason", "elapsed_min", "draw_idx", "Ks_mult", "f_RS_abs"]

# ----------------------------------------------------------------------
# Forcing / file-reading constants (same as run_cc_storm_110.py)
# ----------------------------------------------------------------------
SRC_SDF_DEFAULT = "../init_data/met/Master_Precip.sdf"
RAIN_DT_HR   = float(getattr(builder, "RAIN_INTERVAL", 0.25))   # one MDF row = this many hours
NO_DATA_MIN  = 9999.0                                           # tRIBS NO_DATA flag is 9999.99
OBS_PEAK_RTOL   = 0.02      # scorer's observed peak must match the expected file's within this
GATE_MIN_CHANGE = 0.005     # the storm truth must differ from the 1.0x truth by more than 0.5%
DRY_PROBLEM = "no flow at all"
EXPECTED_QOUT_COLS = 3

# Self-consistency check pass criteria (same precedent as verify_truth_provenance_110.py)
VERIFY_PBIAS_TOL = 0.01     # percent
VERIFY_KGE_MIN   = 0.9999

COMPETING_PATTERNS = (
    "run_sensitivity_single", "run_sensitivity_sweep", "build_sensitivity_run",
    "run_cc_", "run_lhs_", "run_crossversion", "verify_truth_provenance",
    "validate_truth_point", "rescore_series100", "run_truth_",
)


class Abort(RuntimeError):
    """A problem that would make every further run pointless: stop the whole script."""


class BuildFailed(RuntimeError):
    """The builder raised for one run: log it and carry on."""


# ======================================================================
# Names
# ======================================================================
def leg_names(scale):
    """(label, short, leg) for a rain scale: 0.8 -> ('storm080', 'st080', '0.80x')."""
    for s in LEG_SCALES:
        if abs(float(scale) - s) < 1e-9:
            pct = int(round(s * 100))
            return "storm%03d" % pct, "st%03d" % pct, "%.2fx" % s
    raise ValueError("--rain_scale must be 0.8 or 1.25 (the two Series 2 storm legs); got %r. The 1.0x "
                     "leg is run_lhs_Ks_f_120_tribs6.py." % (scale,))


def derived_name(path_str, suffix):
    """'../init_data/met/precip_SMF_2014-2014.mdf' + 'storm080' -> '../init_data/met/precip_SMF_2014-2014_storm080.mdf'."""
    head, sep, name = path_str.rpartition("/")
    stem, dot, ext = name.rpartition(".")
    if not dot:
        return "%s%s%s_%s" % (head, sep, name, suffix)
    return "%s%s%s_%s.%s" % (head, sep, stem, suffix, ext)


def label(v):
    """The builder's own number-to-text rule (7.0 -> 7p0, 0.012 -> 0p012)."""
    return builder.value_to_label(v)


# ======================================================================
# Sampler -- COPIED UNCHANGED from run_lhs_Ks_f_120_tribs6.py so the draws are identical
# ======================================================================
def generate_lhs_samples(n, params, seed=None):
    """Log- or linear-stratified sampling. Each parameter is split into n equal
    slices (log slices for log-scale parameters), one uniform draw per slice,
    then shuffled independently. Same seed + same n + same bounds + same dict
    order = identical draws."""
    rng     = np.random.default_rng(seed)
    samples = {}
    for param, bounds in params.items():
        lo, hi = bounds["lo"], bounds["hi"]
        scale  = bounds.get("scale", "linear")

        if scale == "log":
            log_lo, log_hi = np.log10(lo), np.log10(hi)
            intervals  = np.linspace(log_lo, log_hi, n + 1)
            log_points = rng.uniform(intervals[:-1], intervals[1:])
            points     = 10 ** log_points
        else:
            intervals = np.linspace(lo, hi, n + 1)
            points    = rng.uniform(intervals[:-1], intervals[1:])

        rng.shuffle(points)
        samples[param] = points
    return pd.DataFrame(samples)


# ======================================================================
# Teaching the (unedited) builder about series 120 and this leg's category
# ======================================================================
def install_series_120_hooks(category):
    """In-memory only: tell the builder that Ks_mult runs are series 120 and that
    series 120 lives in this leg's own category folder. The builder file on disk
    is not touched. Stops the script if the result is not what the rest of this
    script assumes."""
    cfg = dict(builder.PARAM_CONFIG["Ks_mult"])
    cfg["series"] = SERIES
    builder.PARAM_CONFIG["Ks_mult"] = cfg

    original = getattr(builder, "_orig_get_run_category_120", None) or builder.get_run_category
    builder._orig_get_run_category_120 = original

    def get_run_category_120(series_str):
        if str(series_str) == SERIES:
            return category
        return original(series_str)

    builder.get_run_category = get_run_category_120

    run_id, _ = builder.build_run_id("Ks_mult", 7.0)
    expected  = "%s_%s_%s_Ks7p0x" % (builder.LOCATION, builder.EVENT_DATE, SERIES)
    if run_id != expected:
        raise RuntimeError(
            "Series-120 hook self-check failed: build_run_id gave %r, expected %r. The builder's naming has "
            "changed; stop and look at build_sensitivity_run_tribs6.py before running anything." % (run_id, expected))
    cat = builder.get_run_category(builder.PARAM_CONFIG["Ks_mult"]["series"])
    if cat != category:
        raise RuntimeError("Series-120 hook self-check failed: category %r, expected %r." % (cat, category))


def predicted_run_id(ks, tag):
    run_id, _ = builder.build_run_id("Ks_mult", ks)
    return "%s_%s" % (run_id, tag) if tag else run_id


def sweep_tag(f, seed, short):
    """Suffix appended to the builder's run_id. Carries f (the builder's own run_id only
    encodes Ks), the seed if it is not the default, and the storm."""
    tag = "f%s" % label(f)
    if seed != DEFAULT_SEED:
        tag += "_sd%d" % seed
    return "%s_%s" % (tag, short)


def variant_suffix(n, seed):
    s = ""
    if n != DEFAULT_N:
        s += "_n%d" % n
    if seed != DEFAULT_SEED:
        s += "_seed%d" % seed
    return s


def is_truth_point(ks, f):
    return abs(ks - TRUTH_KS) < TRUTH_TOL and abs(f - TRUTH_F) < TRUTH_TOL


def in_caution_zone(ks, f):
    return ks >= CAUTION_KS_LO and f >= CAUTION_F_LO


def pinned_overrides(ks, f):
    """Everything that must be pinned for a Series 120 run."""
    return {
        "Ks_mult":          ks,
        "f_RS_abs":         f,
        "kinemvelcoef":     ROUTING_TRUTH["kinemvelcoef"],
        "flowexp":          ROUTING_TRUTH["flowexp"],
        "channelroughness": ROUTING_TRUTH["channelroughness"],
        "optpercolation":   0,
    }


# ======================================================================
# Small helpers
# ======================================================================
def md5_of(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def _num(x, default=np.nan):
    try:
        v = float(x)
        return v if np.isfinite(v) else default
    except (TypeError, ValueError):
        return default


def _truthy(v):
    if isinstance(v, str):
        return v.strip().lower() in ("true", "1")
    try:
        return bool(v) and not (isinstance(v, float) and np.isnan(v))
    except (TypeError, ValueError):
        return False


def _close(a, b, rtol=1e-2):
    """Equal within 1% (catches an ignored override; tolerates the builder rounding a value)."""
    try:
        return math.isclose(float(a), float(b), rel_tol=rtol, abs_tol=1e-12)
    except (TypeError, ValueError):
        return False


def _rank(a):
    return pd.Series(a).rank().to_numpy(dtype=float)


def write_json_atomic(path, obj):
    tmp = path.with_name(path.name + ".part")
    tmp.write_text(json.dumps(obj, indent=2, default=str))
    os.replace(tmp, path)


def _ancestor_pids():
    pids, pid = set(), os.getpid()
    for _ in range(32):
        if pid in pids or pid <= 1:
            break
        pids.add(pid)
        try:
            with open("/proc/%d/status" % pid) as fh:
                m = re.search(r"^PPid:\s*(\d+)", fh.read(), re.M)
            pid = int(m.group(1)) if m else 0
        except Exception:
            break
    return pids


def find_competing_jobs():
    """Other python processes running a script that shares current_run_config.json.
    Returns a list of (pid, command line), or None if the check could not be done."""
    try:
        out = subprocess.run(["ps", "-eo", "pid=,args="], capture_output=True,
                             text=True, timeout=10).stdout
    except Exception:
        return None
    mine = _ancestor_pids()
    found = []
    for line in out.splitlines():
        parts = line.strip().split(None, 1)
        if len(parts) != 2 or not parts[0].isdigit():
            continue
        pid, args = int(parts[0]), parts[1]
        if pid in mine:
            continue
        if not re.search(r"\bpython[\d.]*\b", args):
            continue
        if any(p in args for p in COMPETING_PATTERNS) and ".py" in args:
            found.append((pid, args.strip()))
    return found


def run_scorer(timeout_sec):
    """Run the scorer in its own process group so a hang can be killed as a unit.
    Returns (returncode, elapsed_sec, timed_out, tribs_warning, stdout)."""
    proc = subprocess.Popen(
        [sys.executable, SCORER_SCRIPT],
        cwd=Path.cwd(),
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
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        except (ProcessLookupError, PermissionError):
            pass
        stdout, _ = proc.communicate()
        timed_out = True
    except KeyboardInterrupt:
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        except (ProcessLookupError, PermissionError):
            pass
        proc.communicate()
        print("\nInterrupted: the running tRIBS job was stopped. Finished runs are kept; "
              "re-run the same command with --skip_existing to continue.")
        raise
    elapsed = time.time() - t0
    if stdout:
        print(stdout, end="" if stdout.endswith("\n") else "\n")
    tribs_warning = bool(stdout) and "WARNING: tRIBS may have failed" in stdout
    returncode = None if timed_out else proc.returncode
    return returncode, elapsed, timed_out, tribs_warning, stdout or ""


# ======================================================================
# Reading and checking outlet (.qout) files
# ======================================================================
def read_qout(path):
    """Parse a tRIBS outlet file by hand. One header line, then comma and/or whitespace
    separated rows. Returns (time_hr, q_m3s, column_counts_seen, n_unparseable_lines)."""
    t, q, ncols, bad = [], [], set(), 0
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
                    continue                      # header line
            try:
                tv, qv = float(toks[0]), float(toks[1])
            except (ValueError, IndexError):
                bad += 1
                continue
            t.append(tv)
            q.append(qv)
            ncols.add(len(toks))
    return np.array(t, dtype=float), np.array(q, dtype=float), ncols, bad


def describe_qout(path):
    """Sanity-check a .qout file. Never raises. qout_ok is True only if nothing looked wrong."""
    info = {"qout_ok": False, "qout_rows": 0, "qout_ncols": 0,
            "qout_dt_min": np.nan, "qout_peak_m3s": np.nan,
            "qout_tpeak_hr": np.nan, "qout_volume_m3": np.nan,
            "qout_problem": ""}
    try:
        t, flow, ncols, bad = read_qout(path)
    except Exception as e:
        info["qout_problem"] = "unreadable: %s" % e
        return info

    problems = []
    info["qout_rows"] = int(len(t))
    info["qout_ncols"] = int(max(ncols)) if ncols else 0
    if bad:
        problems.append("%d unparseable line(s)" % bad)
    if len(t) < 10:
        problems.append("only %d data rows" % len(t))
    elif ncols != {EXPECTED_QOUT_COLS}:
        problems.append("expected %d columns, saw %s" % (EXPECTED_QOUT_COLS, sorted(ncols)))
    if len(t) and not (np.isfinite(t).all() and np.isfinite(flow).all()):
        problems.append("non-finite values (nan/inf)")

    if len(t) >= 10 and np.isfinite(t).all() and np.isfinite(flow).all():
        if not np.all(np.diff(t) > 0):
            problems.append("time not strictly increasing")
        if (flow < 0).any():
            problems.append("negative discharge")
        info["qout_dt_min"] = float(np.median(np.diff(t)) * 60.0)
        info["qout_peak_m3s"] = float(np.max(flow))
        info["qout_tpeak_hr"] = float(t[int(np.argmax(flow))])
        info["qout_volume_m3"] = float(
            np.sum(0.5 * (flow[1:] + flow[:-1]) * np.diff(t)) * 3600.0)
        if not (info["qout_peak_m3s"] > 0):
            problems.append(DRY_PROBLEM)

    info["qout_problem"] = "; ".join(problems)
    info["qout_ok"] = not problems
    return info


def locate_raw_qout(raw_dir, run_id):
    expected = raw_dir / ("%s_Outlet.qout" % run_id)
    if expected.exists():
        return expected
    found = list(raw_dir.glob("*.qout")) if raw_dir.exists() else []
    return found[0] if len(found) == 1 else None


def copy_verified(src, dst):
    """Copy src -> dst via a temporary name, verify size and checksum, then move into place.
    Returns the md5. Raises on any mismatch (no partial file under the final name)."""
    tmp = dst.with_name(dst.name + ".part")
    shutil.copy2(src, tmp)
    if tmp.stat().st_size != Path(src).stat().st_size:
        tmp.unlink(missing_ok=True)
        raise IOError("size mismatch after copy")
    md5_src, md5_tmp = md5_of(src), md5_of(tmp)
    if md5_src != md5_tmp:
        tmp.unlink(missing_ok=True)
        raise IOError("checksum mismatch after copy")
    os.replace(tmp, dst)
    return md5_tmp


def cleanup_raw_results(run_id, calib_dir, category):
    """Delete a successful run's raw results folder (pixel files, spatial snapshots). Only
    deletes a folder located EXACTLY at 02_results/<category>/<run_id>. Never raises."""
    try:
        parent = (calib_dir / "02_results" / category).resolve()
        raw    = (parent / run_id).resolve()
        if raw.parent != parent or raw.name != run_id:
            print("  (cleanup skipped: %s is not directly inside %s)" % (raw, parent))
            return False
        if raw.exists():
            shutil.rmtree(raw)
            return True
    except Exception as e:
        print("  (cleanup warning for %s: %s)" % (run_id, e))
    return False


# ======================================================================
# FORCING: scaled copies of the two rain files and the SDF
# (split_keep / join_keep / scale_mdf_text / plan_forcing / write_forcing /
#  verify_forcing follow run_cc_storm_110.py, EXCEPT that the rain column is now found
#  by rain_column(), which accepts the real header  Year,Month,Day,Hour,R_mm  as well as
#  the old  Y M D H R. The 110 script's own header check does not accept the real header.)
# ======================================================================
MDF_TOKEN = re.compile(r"[^\s,]+\.mdf\b", re.IGNORECASE)     # a path: no spaces, no commas
_SEP = re.compile(r"(\s*,\s*|\s+)")                          # a column separator: comma or whitespace


def split_keep(line):
    """Split a line into (tokens, separators), KEEPING the exact separators so the line can be put
    back together unchanged. Works for space-, tab- and comma-delimited files."""
    parts = _SEP.split(line.strip())
    return parts[0::2], parts[1::2]


def join_keep(tokens, seps):
    out = [tokens[0]]
    for sep, tok in zip(seps, tokens[1:]):
        out.append(sep)
        out.append(tok)
    return "".join(out)


_R_NAMES = {"R", "RAIN", "RAINFALL", "PRECIP", "PRECIPITATION"}


def rain_column(header_tokens, where, first_lines):
    """Find the rainfall-rate column from a rain file's header, and confirm there is a year column.

    The header the tRIBS 6.0.0 rain files really have is  Year,Month,Day,Hour,R_mm  (comma-separated);
    older files had  Y M D H R  (spaces). Both are accepted, in any capitalisation, and so are the
    obvious variants (R_mm, Rain, Precip). Anything else stops with the first lines of the file
    printed, so nothing is guessed. The column is a RATE in mm/hr whatever it is called (the same
    reading check_storm_rain_110.py uses, which reproduced the storm totals)."""
    names = [h.strip("\"'").strip().upper() for h in header_tokens]
    r_idx = [i for i, n in enumerate(names)
             if n in _R_NAMES or n.startswith("R_") or n.startswith("RAIN") or n.startswith("PRECIP")]
    y_idx = [i for i, n in enumerate(names) if n == "Y" or n.startswith("YEAR")]
    if len(r_idx) != 1 or not y_idx:
        raise ValueError(
            "cannot tell which column is rainfall in %s. Its header reads %r, and I need exactly one "
            "rainfall column (R, R_mm, Rain, Precip) plus a year column (Y or Year). The first lines of "
            "the file are:\n    %s\nSend me those lines and I will adapt."
            % (where, [h.strip() for h in header_tokens], "\n    ".join(l.rstrip() for l in first_lines)))
    return r_idx[0]


def scale_mdf_text(src_path, factor):
    """Scale the R (rainfall rate) column of a rain gauge file by `factor`. Every other column
    is copied through as text. NO_DATA flags (>= 9999) and negative values are left alone.
    Returns (new_text, stats)."""
    raw = Path(src_path).read_text().splitlines()
    lines = [ln for ln in raw if ln.strip()]
    if not lines:
        raise ValueError("%s is empty" % src_path)
    header, _ = split_keep(lines[0])
    ri = rain_column(header, src_path, lines[:3])
    out = [lines[0]]
    st = {"rows": 0, "scaled": 0, "kept_flag": 0, "sum_before": 0.0, "sum_after": 0.0,
          "peak_before": 0.0, "peak_after": 0.0, "header": lines[0].strip(),
          "rain_col": header[ri].strip("\"'").strip()}
    for ln in lines[1:]:
        tk, sp = split_keep(ln)
        if len(tk) != len(header):
            raise ValueError("%s: a row has %d columns but the header has %d: %r"
                             % (src_path, len(tk), len(header), ln))
        try:
            r = float(tk[ri])
        except ValueError:
            raise ValueError("%s: cannot read %r as a rainfall value in the row %r" % (src_path, tk[ri], ln))
        st["rows"] += 1
        if r >= NO_DATA_MIN or r < 0:
            st["kept_flag"] += 1
            out.append(join_keep(tk, sp))
            continue
        r2 = r * factor
        tk[ri] = "%.6f" % r2
        out.append(join_keep(tk, sp))
        st["scaled"] += 1
        st["sum_before"] += r
        st["sum_after"] += r2
        st["peak_before"] = max(st["peak_before"], r)
        st["peak_after"] = max(st["peak_after"], r2)
    if st["rows"] == 0:
        raise ValueError("%s has a header but no data rows" % src_path)
    return "\n".join(out) + "\n", st


def plan_forcing(src_sdf_rel, suffix, scale):
    """Work out (in memory, writing nothing) the scaled rain files and SDF."""
    src_sdf = Path.cwd() / src_sdf_rel
    if not src_sdf.exists():
        raise FileNotFoundError(
            "rain gauge SDF not found: %s (looked relative to %s). Run this script from lab/, or give "
            "--src_sdf." % (src_sdf_rel, Path.cwd()))
    lines = src_sdf.read_text().splitlines()
    if not lines:
        raise ValueError("%s is empty" % src_sdf)
    out_lines, stations = [lines[0]], []
    for line in lines[1:]:
        m = MDF_TOKEN.search(line)
        if m is None:
            out_lines.append(line)
            continue
        tok = m.group(0)
        src_mdf = Path.cwd() / tok
        if not src_mdf.exists():
            raise FileNotFoundError(
                "%s line %r names the rain file %s, which does not exist relative to %s. tRIBS reads these "
                "paths from where it is run (lab/), so the scaled copy needs the same layout. If the file is "
                "somewhere else, send me the first lines of %s."
                % (src_sdf.name, line.strip()[:160], tok, Path.cwd(), src_sdf.name))
        new_tok = derived_name(tok, suffix)
        new_text, st = scale_mdf_text(src_mdf, scale)
        stations.append({"src": src_mdf, "dst": src_mdf.with_name(Path(new_tok).name),
                         "token": tok, "new_token": new_tok, "text": new_text, "stats": st})
        out_lines.append(line[:m.start()] + new_tok + line[m.end():])
    if not stations:
        raise ValueError("no .mdf file names found in %s; its first lines are %r" % (src_sdf, lines[:3]))
    try:
        declared = int(re.split(r"[\s,]+", lines[0].strip())[0])
        if declared != len(stations):
            raise ValueError("%s says %d station(s) but lists %d .mdf file(s)"
                             % (src_sdf.name, declared, len(stations)))
    except (ValueError, IndexError) as e:
        if "station" in str(e):
            raise
    return {
        "src_sdf": src_sdf, "dst_sdf_rel": derived_name(src_sdf_rel, suffix),
        "dst_sdf": Path.cwd() / derived_name(src_sdf_rel, suffix),
        "sdf_text": "\n".join(out_lines) + "\n", "stations": stations, "scale": scale,
    }


def print_forcing_plan(plan):
    print("  rain gauge files (only the rainfall rate R is multiplied; Y M D H are copied unchanged):")
    for s in plan["stations"]:
        st = s["stats"]
        print("    %-34s %5d rows   peak %7.3f -> %7.3f mm/hr   total %8.2f -> %8.2f mm  (x%g)"
              % (s["dst"].name, st["rows"], st["peak_before"], st["peak_after"],
                 st["sum_before"] * RAIN_DT_HR, st["sum_after"] * RAIN_DT_HR, plan["scale"]))
        if st["kept_flag"]:
            print("      (%d NO_DATA/negative row(s) left unscaled)" % st["kept_flag"])
    cols = sorted({s["stats"].get("rain_col", "?") for s in plan["stations"]})
    print("    rain column read as a rate in mm/hr: %s  (header seen: %s)"
          % (", ".join(cols), "; ".join(sorted({s["stats"]["header"] for s in plan["stations"]}))))
    print("    total assumes each row spans %g h" % RAIN_DT_HR)
    print("    SDF: %s  ->  %s" % (plan["src_sdf"].name, plan["dst_sdf"].name))


def planned_files(plan):
    """[(path, text)] for every file the forcing stage writes."""
    return [(s["dst"], s["text"]) for s in plan["stations"]] + [(plan["dst_sdf"], plan["sdf_text"])]


def forcing_status(plan):
    """[(file name, state)] without writing anything: absent / identical / DIFFERENT."""
    out = []
    for dst, text in planned_files(plan):
        if not dst.exists():
            out.append((dst.name, "will be written"))
        elif dst.read_bytes() == text.encode("utf-8"):
            out.append((dst.name, "already present, identical to what this script writes"))
        else:
            out.append((dst.name, "PRESENT BUT DIFFERENT (the script would stop; see --overwrite_forcing)"))
    return out


def write_forcing(plan, overwrite):
    """Write the planned files (atomically). A file that already exists with identical content
    is left alone; one that differs is an error unless overwrite. Returns {name: status}."""
    status = {}
    for dst, text in planned_files(plan):
        data = text.encode("utf-8")
        if dst.exists():
            if dst.read_bytes() == data:
                status[dst.name] = "already present (identical)"
                continue
            if not overwrite:
                raise Abort(
                    "%s already exists and is DIFFERENT from what this script would write. It may be an "
                    "older storm file with other contents. Move or rename it, or re-run with "
                    "--overwrite_forcing if you are sure." % dst)
        tmp = dst.with_name(dst.name + ".part")
        tmp.write_bytes(data)
        os.replace(tmp, dst)
        status[dst.name] = "written"
    return status


def verify_forcing(plan):
    """Re-read what is on disk and compare with the originals, independently of the code that
    wrote it. Returns a list of problems (empty = fine)."""
    problems = []
    factor = plan["scale"]
    for s in plan["stations"]:
        try:
            a = [split_keep(ln)[0] for ln in Path(s["src"]).read_text().splitlines() if ln.strip()]
            b = [split_keep(ln)[0] for ln in Path(s["dst"]).read_text().splitlines() if ln.strip()]
        except Exception as e:
            problems.append("%s: cannot re-read (%s)" % (s["dst"].name, e))
            continue
        if len(a) != len(b):
            problems.append("%s: %d lines vs %d in the original" % (s["dst"].name, len(b), len(a)))
            continue
        if a[0] != b[0]:
            problems.append("%s: header differs from the original" % s["dst"].name)
            continue
        try:
            ri = rain_column(a[0], s["src"], [",".join(t) for t in a[:3]])
        except ValueError as e:
            problems.append("%s: %s" % (s["dst"].name, e))
            continue
        n_bad = 0
        for ra, rb in zip(a[1:], b[1:]):
            if len(ra) != len(rb) or any(x != y for k, (x, y) in enumerate(zip(ra, rb)) if k != ri):
                n_bad += 1
                continue
            r0, r1 = float(ra[ri]), float(rb[ri])
            want = r0 if (r0 >= NO_DATA_MIN or r0 < 0) else r0 * factor
            if abs(r1 - want) > 1e-6:
                n_bad += 1
        if n_bad:
            problems.append("%s: %d row(s) differ from original x %g (or have changed Y/M/D/H)"
                            % (s["dst"].name, n_bad, factor))
    try:
        a = Path(plan["src_sdf"]).read_text().splitlines()
        b = Path(plan["dst_sdf"]).read_text().splitlines()
        if len(a) != len(b) or a[0] != b[0]:
            problems.append("%s: layout differs from the original SDF" % plan["dst_sdf"].name)
        for s in plan["stations"]:
            if s["new_token"] not in Path(plan["dst_sdf"]).read_text():
                problems.append("%s does not name %s" % (plan["dst_sdf"].name, s["new_token"]))
    except Exception as e:
        problems.append("%s: cannot re-read (%s)" % (plan["dst_sdf"].name, e))
    return problems


def forcing_md5s(plan):
    """{file name: md5} for the scaled files as they are on disk now."""
    return {p.name: md5_of(p) for p, _ in planned_files(plan) if p.exists()}


# ======================================================================
# Building one run, scoring it, and checking what was scored
# ======================================================================
def check_config(ctx, cfg, run_id, ks, f):
    """What the builder wrote into current_run_config.json versus what was asked for."""
    problems = []
    if cfg.get("run_id") != run_id:
        problems.append("run_id is %r, expected %r" % (cfg.get("run_id"), run_id))
    if "gauge_sdf" in cfg and cfg["gauge_sdf"] != ctx.gauge_sdf_rel:
        problems.append("gauge_sdf is %r, expected %r" % (cfg["gauge_sdf"], ctx.gauge_sdf_rel))
    wanted = [("Ks_mult", ks), ("f_RS_abs", f), ("optpercolation", 0),
              ("kinemvelcoef", ROUTING_TRUTH["kinemvelcoef"]), ("flowexp", ROUTING_TRUTH["flowexp"]),
              ("channelroughness", ROUTING_TRUTH["channelroughness"])]
    for key, want in wanted:
        if key not in cfg:
            problems.append("%s missing from the run configuration" % key)
        elif not _close(cfg[key], want):
            problems.append("%s is %r, expected %r" % (key, cfg[key], want))
    if cfg.get("run_category") not in (None, ctx.category):
        problems.append("run_category is %r, expected %r" % (cfg.get("run_category"), ctx.category))
    inp = cfg.get("input_file")
    if not inp:
        problems.append("no input_file in the run configuration")
    else:
        p = Path(inp) if Path(inp).is_absolute() else Path.cwd() / inp
        try:
            text = p.read_text(errors="replace")
            if Path(ctx.gauge_sdf_rel).name not in text:
                problems.append("the .in file does not name %s (the scaled rain would NOT be used)"
                                % Path(ctx.gauge_sdf_rel).name)
        except Exception as e:
            problems.append("cannot read the .in file %s (%s)" % (inp, e))
    return problems


def prepare_run(ctx, pred_id, ks, f, tag, truth_file_rel):
    """Build the input file with the scaled rain, check the configuration the builder wrote, and
    (for runs that must be scored against this storm's truth) add truth_file to it."""
    try:
        run_id, _, _ = builder.build_input_file("Ks_mult", ks, overrides=pinned_overrides(ks, f),
                                                tag=tag, gauge_sdf=ctx.gauge_sdf_rel)
    except Exception as e:
        raise BuildFailed(str(e)[:300])
    if run_id != pred_id:
        raise Abort("builder produced run_id %s, expected %s" % (run_id, pred_id))
    cfg_path = ctx.calib_dir / "current_run_config.json"
    try:
        cfg = json.loads(cfg_path.read_text())
    except Exception as e:
        raise Abort("cannot read %s after the build: %s" % (cfg_path, e))
    problems = check_config(ctx, cfg, run_id, ks, f)
    if problems:
        raise Abort("the builder's run configuration does not match what was asked for: " + "; ".join(problems))
    if truth_file_rel is None:
        if cfg.get("truth_file"):
            raise Abort("the run configuration already carries truth_file=%r although this run is meant to be "
                        "scored against the 1.0x truth in synth_truth/." % cfg.get("truth_file"))
    else:
        cfg["truth_file"] = truth_file_rel
        write_json_atomic(cfg_path, cfg)
        back = json.loads(cfg_path.read_text())
        if back.get("truth_file") != truth_file_rel:
            raise Abort("could not write truth_file into %s" % cfg_path)
    return run_id, cfg


_PEAK_CACHE = {}


def expected_obs_peak(path, event_start, event_end):
    """The observed peak the scorer WILL see if it reads this truth file: same reading, same
    5-minute time-interpolation, same event window."""
    key = (str(path), event_start, event_end)
    if key not in _PEAK_CACHE:
        t, q, _, _ = read_qout(path)
        idx = pd.DatetimeIndex(pd.to_datetime(t * 3600, unit="s", origin=pd.Timestamp("2014-08-01")))
        s = pd.Series(q, index=idx).resample("5min").interpolate(method="time")
        _PEAK_CACHE[key] = float(s.loc[event_start:event_end].max())
    return _PEAK_CACHE[key]


def check_observed(metrics, cfg, stdout, expected_obs_path):
    """Did the scorer really use the file we meant as 'observed'? Returns a list of problems."""
    problems = []
    if str(metrics.get("obs_mode")) != "synth":
        problems.append("obs_mode is %r, expected 'synth' (the scorer fell back to the real gauge)"
                        % metrics.get("obs_mode"))
    m = re.search(r"Reading observed from:\s*(.+)", stdout or "")
    if m and m.group(1).strip() != expected_obs_path.name:
        problems.append("the scorer read %r as observed, expected %r" % (m.group(1).strip(), expected_obs_path.name))
    elif not m:
        print("  (note: the scorer printed no 'Reading observed from:' line, so that check was skipped; "
              "the observed-peak check below still applies)")
    try:
        exp = expected_obs_peak(expected_obs_path, cfg["event_start"], cfg["event_end"])
        got = float(metrics["obs_peak_m3s"])
        if not (abs(got - exp) <= OBS_PEAK_RTOL * max(abs(exp), 1e-12)):
            problems.append("the observed peak the scorer used (%.4g m3/s) is not the peak of %s (%.4g m3/s)"
                            % (got, expected_obs_path.name, exp))
    except Exception as e:
        print("  (note: could not double-check the observed peak: %s)" % str(e)[:100])
    return problems


def readback_problems(m, ks, f):
    """Compare the scored metrics row against what this script asked for."""
    probs = []
    expected = {
        "Ks_mult":          ks,
        "f_RS_abs":         f,
        "kinemvelcoef":     ROUTING_TRUTH["kinemvelcoef"],
        "flowexp":          ROUTING_TRUTH["flowexp"],
        "channelroughness": ROUTING_TRUTH["channelroughness"],
        "optpercolation":   0,
    }
    for key, want in expected.items():
        got = _num(m.get(key))
        if not np.isfinite(got) or not np.isclose(got, want, rtol=1e-6, atol=1e-12):
            probs.append("%s: scored row says %r, asked for %r" % (key, m.get(key), want))
    return probs


def run_one(ctx, pred_id, ks, f, tag, truth_file_rel, expected_obs_path):
    """Build, score, and verify one run. Returns a dict with status in SUCCESS / BUILD_FAILED /
    HANG / FAILED, plus the metrics row. Raises Abort for a problem that would affect every run."""
    out = {"run_id": pred_id, "status": None, "reason": "", "elapsed": 0.0, "metrics": None}
    try:
        run_id, cfg = prepare_run(ctx, pred_id, ks, f, tag, truth_file_rel)
    except BuildFailed as e:
        out.update(status="BUILD_FAILED", reason=str(e))
        return out

    # Remove this run's old outputs so a stale file can never be mistaken for a new one.
    (ctx.csv_dir / ("%s_compare_obs_sim.csv" % run_id)).unlink(missing_ok=True)
    (ctx.summary_dir / ("%s_metrics_summary.csv" % run_id)).unlink(missing_ok=True)

    returncode, elapsed, timed_out, tribs_warning, stdout = run_scorer(ctx.args.timeout)
    out["elapsed"] = elapsed
    if timed_out:
        out.update(status="HANG", reason="wall-clock timeout")
        return out
    if returncode != 0 or tribs_warning:
        out.update(status="FAILED",
                   reason="tRIBS reported non-zero exit" if tribs_warning else "scorer exited %s" % returncode)
        return out
    metrics_file = ctx.summary_dir / ("%s_metrics_summary.csv" % run_id)
    if not metrics_file.exists():
        out.update(status="FAILED", reason="no metrics file %s" % metrics_file.name)
        return out
    try:
        metrics = pd.read_csv(metrics_file).iloc[0].to_dict()
    except Exception as e:
        out.update(status="FAILED", reason="metrics file unreadable: %s" % e)
        return out
    out["metrics"] = metrics
    problems = check_observed(metrics, cfg, stdout, expected_obs_path)
    if problems:
        out["abort"] = ("scoring did not use the intended 'observed' file for %s: %s. Every later run would be "
                        "scored the same wrong way, so the script is stopping." % (run_id, "; ".join(problems)))
        out["status"] = "OBSERVED_MISMATCH"
        out["reason"] = "; ".join(problems)
        return out
    probs = readback_problems(metrics, ks, f)
    if probs:
        out["abort"] = ("the scored row for %s does not match what this script asked for: %s. The override "
                        "plumbing is broken, so the remaining runs would be wasted." % (run_id, "; ".join(probs)))
        out["status"] = "READBACK_MISMATCH"
        out["reason"] = "; ".join(probs)
        return out
    out["status"] = "SUCCESS"
    return out


def note_outcome(ctx, ok):
    """Stop the whole script if too many runs in a row fail (a systematic problem)."""
    if ok:
        ctx.consec_fail = 0
        return
    ctx.consec_fail += 1
    if ctx.consec_fail >= ctx.args.max_consecutive_failures:
        raise Abort("%d runs in a row failed. That looks like a systematic problem, not bad luck. Look at "
                    "the last error above, fix it, and re-run the same command with --skip_existing."
                    % ctx.consec_fail)


def harvest_qout(ctx, run_id, dest_path):
    """Copy a finished run's outlet file to dest_path (verified). Returns
    (status, reason, md5, qinfo)."""
    raw_dir = ctx.calib_dir / "02_results" / ctx.category / run_id
    raw_qout = locate_raw_qout(raw_dir, run_id)
    if raw_qout is None:
        return "NO_QOUT", "no single *.qout found in %s" % raw_dir, None, None
    try:
        md5 = copy_verified(raw_qout, dest_path)
    except Exception as e:
        return "COPY_FAILED", str(e)[:300], None, None
    return "OK", "", md5, describe_qout(dest_path)


# ======================================================================
# Provenance of this storm's truth
# ======================================================================
def pytribs_version():
    try:
        from importlib import metadata
        for name in ("pytRIBS", "pytribs"):
            try:
                return metadata.version(name)
            except metadata.PackageNotFoundError:
                continue
    except Exception:
        pass
    return None


def load_truth_prov(ctx):
    if not ctx.prov_path.exists():
        return None
    try:
        return json.loads(ctx.prov_path.read_text())
    except Exception:
        return None


def save_truth_prov(ctx, prov):
    prov["updated_local"] = datetime.now().isoformat(timespec="seconds")
    ctx.truth_dir.mkdir(parents=True, exist_ok=True)
    write_json_atomic(ctx.prov_path, prov)


def truth_ids(ctx):
    """(truth run_id, truth tag, verify run_id, verify tag). Hooks must already be installed."""
    ttag = "f%s_truth_%s" % (label(TRUTH_F), ctx.short)
    vtag = "f%s_verify_%s" % (label(TRUTH_F), ctx.short)
    return predicted_run_id(TRUTH_KS, ttag), ttag, predicted_run_id(TRUTH_KS, vtag), vtag


# ======================================================================
# STAGE: forcing
# ======================================================================
def stage_forcing(ctx):
    print("\n" + "=" * 92)
    print("STEP 1 of 5 -- FORCING: scaled rain files for %s (x%g)" % (ctx.label, ctx.scale))
    print("=" * 92)
    status = write_forcing(ctx.forcing_plan, ctx.args.overwrite_forcing)
    for name, st in status.items():
        print("  %-36s %s" % (name, st))
    problems = verify_forcing(ctx.forcing_plan)
    if problems:
        raise Abort("the scaled rain files do not match the originals: " + "; ".join(problems))
    print("  CHECKED against the originals: same rows, same Y/M/D/H, rain rate = original x %g, "
          "SDF points at the scaled files." % ctx.scale)
    plan = ctx.forcing_plan
    prov = {
        "script": Path(__file__).name,
        "created_local": datetime.now().isoformat(timespec="seconds"),
        "rain_scale": ctx.scale, "label": ctx.label,
        "what_this_is": "Rain gauge forcing with the rainfall RATE column multiplied by a constant; "
                        "Y/M/D/H unchanged; NO_DATA flags unchanged; whole record scaled. Same method and "
                        "same file names as run_cc_storm_110.py.",
        "row_hours_assumed_for_totals": RAIN_DT_HR,
        "source_sdf": {"file": plan["src_sdf"].name, "md5": md5_of(plan["src_sdf"])},
        "scaled_sdf": {"file": plan["dst_sdf"].name, "md5": md5_of(plan["dst_sdf"])},
        "stations": [{
            "source": s["src"].name, "source_md5": md5_of(s["src"]),
            "scaled": s["dst"].name, "scaled_md5": md5_of(s["dst"]),
            "rows": s["stats"]["rows"], "peak_mm_hr": [s["stats"]["peak_before"], s["stats"]["peak_after"]],
            "total_mm": [s["stats"]["sum_before"] * RAIN_DT_HR, s["stats"]["sum_after"] * RAIN_DT_HR],
            "no_data_rows_unscaled": s["stats"]["kept_flag"],
        } for s in plan["stations"]],
        "python": sys.version.split()[0], "pandas": pd.__version__, "numpy": np.__version__,
    }
    name = "forcing_120_%s_PROVENANCE.json" % ctx.label
    write_json_atomic(ctx.calib_dir / name, prov)
    print("  Saved: %s" % name)


# ======================================================================
# STAGE: this storm's truth, with the direction check and self-consistency check
# ======================================================================
def gate_truth(ctx, metrics, qinfo_new=None):
    """The storm's truth must differ from the 1.0x truth in the right direction. Returns a dict
    of the numbers. Raises Abort if the direction is wrong."""
    obs_vol, sim_vol = _num(metrics.get("obs_volume_m3")), _num(metrics.get("sim_volume_m3"))
    obs_pk, sim_pk = _num(metrics.get("obs_peak_m3s")), _num(metrics.get("sim_peak_m3s"))
    vr = sim_vol / obs_vol if np.isfinite(obs_vol) and obs_vol else np.nan
    pr = sim_pk / obs_pk if np.isfinite(obs_pk) and obs_pk else np.nan
    print("  DIRECTION CHECK: this storm's truth has %.1f%% of the 1.0x truth's peak and %.1f%% of its volume "
          "(rain x%g)." % (100 * pr, 100 * vr, ctx.scale))
    print("    (for comparison, on the old build 5.3.1 the peak was %.1f%% of the 1.0x truth's)"
          % (100 * OLD_BUILD_PEAK_RATIO[ctx.scale]))
    if not (np.isfinite(vr) and np.isfinite(pr)):
        raise Abort("could not compute the truth's peak/volume relative to the 1.0x truth (the scorer's "
                    "obs_*/sim_* columns were missing or not finite).")
    direction = np.sign(ctx.scale - 1.0)
    if (np.sign(pr - 1.0) != direction or np.sign(vr - 1.0) != direction
            or abs(pr - 1.0) < GATE_MIN_CHANGE or abs(vr - 1.0) < GATE_MIN_CHANGE):
        raise Abort(
            "DIRECTION CHECK FAILED. With rain x%g the peak and volume should both move %s from the 1.0x "
            "truth's, but they are %.1f%% and %.1f%% of it. The scaled rainfall is apparently NOT being used "
            "by tRIBS (check that %s is what the .in file names, and that tRIBS reads it). Nothing was saved."
            % (ctx.scale, "DOWN" if direction < 0 else "UP", 100 * pr, 100 * vr, ctx.gauge_sdf_rel))
    print("  DIRECTION CHECK passed: the change is in the right direction.")
    return {"peak_ratio_vs_1p0": pr, "volume_ratio_vs_1p0": vr,
            "old_build_peak_ratio": OLD_BUILD_PEAK_RATIO[ctx.scale],
            "obs_peak_m3s": obs_pk, "sim_peak_m3s": sim_pk}


def crosscheck_series110(ctx):
    """Compare the new truth with Series 110's cc-OFF reference for this storm, if there is one
    (same build, same parameters, so it should be identical). Information only."""
    d = ctx.calib_dir / ("synth_truth_ccoff_%s" % ctx.label)
    files = sorted(d.glob("*_Outlet.qout")) if d.exists() else []
    if not files:
        print("  Series 110 cc-OFF reference for %s: not found (%s). Nothing to compare; that is fine."
              % (ctx.label, d.name))
        return {"status": "not found"}
    ref = next((p for p in files if "_ref_" in p.name), files[0])
    try:
        t0, q0, _, _ = read_qout(ctx.truth_path)
        t1, q1, _, _ = read_qout(ref)
    except Exception as e:
        print("  Series 110 reference %s could not be read (%s)." % (ref.name, e))
        return {"status": "unreadable", "file": ref.name}
    res = {"file": ref.name, "md5": md5_of(ref)}
    if md5_of(ctx.truth_path) == res["md5"]:
        print("  Series 110 cc-OFF reference (%s): BYTE-IDENTICAL to the new truth." % ref.name)
        res["status"] = "identical (bytes)"
        return res
    if len(t0) == len(t1) and np.allclose(t0, t1) and np.allclose(q0, q1, rtol=1e-6, atol=1e-9):
        print("  Series 110 cc-OFF reference (%s): same numbers as the new truth (not byte-identical, "
              "but equal to 1e-6)." % ref.name)
        res["status"] = "identical (numbers)"
        return res
    pk0, pk1 = (q0.max() if len(q0) else np.nan), (q1.max() if len(q1) else np.nan)
    v = lambda t, q: float(np.sum(0.5 * (q[1:] + q[:-1]) * np.diff(t)) * 3600.0) if len(t) > 1 else np.nan
    print("  Series 110 cc-OFF reference (%s) DIFFERS from the new truth: peak %.4g vs %.4g m3/s, volume %.4g vs "
          "%.4g m3, %d vs %d rows." % (ref.name, pk1, pk0, v(t1, q1), v(t0, q0), len(t1), len(t0)))
    print("  (Information only; both are valid truths for their own series. A difference would mean the build or "
          "the inputs changed between Series 110 and now. Worth telling me.)")
    res.update(status="DIFFERENT", peak_new=float(pk0), peak_110=float(pk1))
    return res


def run_selfcheck(ctx):
    """Run the truth point a second time and score it against this storm's truth via truth_file.
    Returns (passed, info dict). Raises Abort if the scorer used the wrong observed file."""
    truth_id, ttag, vid, vtag = truth_ids(ctx)
    print("\n" + "-" * 92)
    print("SELF-CONSISTENCY CHECK: run the truth point again, score it against %s" % ctx.truth_path.name)
    print("  Expect PBIAS +0.0000%% and KGE_2012 1.0000 (need |PBIAS| < %g%% and KGE_2012 > %g)."
          % (VERIFY_PBIAS_TOL, VERIFY_KGE_MIN))
    print("-" * 92, flush=True)
    res = run_one(ctx, vid, TRUTH_KS, TRUTH_F, vtag, ctx.truth_rel, ctx.truth_path)
    if res.get("abort"):
        raise Abort(res["abort"])
    info = {"run_id": vid, "when": datetime.now().isoformat(timespec="seconds"), "passed": False,
            "status": res["status"], "reason": res["reason"]}
    if res["status"] != "SUCCESS":
        print("  Self-check run %s: %s" % (res["status"], res["reason"]))
        return False, info
    m = res["metrics"]
    pb, kg = _num(m.get("pbias_pct")), _num(m.get("kge_2012"))
    info.update(pbias_pct=pb, kge_2012=kg, obs_mode=str(m.get("obs_mode")),
                elapsed_min=res["elapsed"] / 60.0)
    passed = bool(np.isfinite(pb) and np.isfinite(kg) and abs(pb) < VERIFY_PBIAS_TOL and kg > VERIFY_KGE_MIN)
    info["passed"] = passed
    print("\n  obs_mode=%s   PBIAS=%+.4f%%   KGE_2012=%.4f   (this run took %.2f min -> the 400-run sweep is "
          "roughly %.1f h)" % (m.get("obs_mode"), pb, kg, res["elapsed"] / 60.0, 400 * res["elapsed"] / 3600.0))
    if passed and not ctx.args.no_cleanup:
        cleanup_raw_results(vid, ctx.calib_dir, ctx.category)
    return passed, info


def stage_make_truth(ctx):
    truth_id, ttag, vid, vtag = truth_ids(ctx)
    stage_forcing(ctx)

    print("\n" + "=" * 92)
    print("STEP 2 of 5 -- THE TRUTH for %s: Ks %g, f %g, routing at truth, optpercolation 0, scaled rain"
          % (ctx.label, TRUTH_KS, TRUTH_F))
    print("=" * 92)
    ctx.truth_dir.mkdir(parents=True, exist_ok=True)
    prov = load_truth_prov(ctx)
    have = ctx.truth_path.exists()
    if have and ctx.args.remake_truth:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        for p in (ctx.truth_path, ctx.prov_path):
            if p.exists():
                new = p.with_name("%s.superseded_%s" % (p.name, stamp))
                os.replace(p, new)
                print("  --remake_truth: kept the old file as %s" % new.name)
        have, prov = False, None

    if have:
        print("  Truth already exists: %s" % ctx.truth_path.name)
        if prov is None or "truth" not in prov:
            raise Abort("the truth file exists but its provenance record (%s) is missing or unreadable. Do not "
                        "guess: re-run with --remake_truth to rebuild both (the old file is kept, renamed)."
                        % ctx.prov_path.name)
        md5 = md5_of(ctx.truth_path)
        if md5 != prov["truth"].get("md5"):
            raise Abort("the truth file's checksum (%s) does not match its provenance record (%s). The file has "
                        "changed since it was made. Use --remake_truth to rebuild (the old file is kept)."
                        % (md5[:12], str(prov["truth"].get("md5"))[:12]))
        qinfo = describe_qout(ctx.truth_path)
        if not qinfo["qout_ok"]:
            raise Abort("the existing truth file fails its sanity check: %s" % qinfo["qout_problem"])
        now = forcing_md5s(ctx.forcing_plan)
        changed = [n for n, m5 in prov.get("files_md5", {}).items() if now.get(n) != m5]
        if changed:
            raise Abort("the scaled rain file(s) %s have changed since this truth was made (checksum differs), so "
                        "the truth no longer matches the rain on disk. Use --remake_truth to rebuild the truth "
                        "(the old one is kept, renamed)." % ", ".join(changed))
        print("  Checksum matches the provenance record (%s) and the rain files are unchanged. Keeping it; "
              "repeating the self-consistency check only." % md5[:12])
    else:
        print("  -> %s" % truth_id, flush=True)
        # Scored against the 1.0x truth in synth_truth/ (no truth_file): that gives the direction check.
        res = run_one(ctx, truth_id, TRUTH_KS, TRUTH_F, ttag, None, ctx.storm1_ref)
        if res.get("abort"):
            raise Abort(res["abort"])
        if res["status"] != "SUCCESS":
            raise Abort("the truth run %s: %s. Cannot continue without it." % (res["status"], res["reason"]))
        gate = gate_truth(ctx, res["metrics"])
        st, why, md5, qinfo = harvest_qout(ctx, truth_id, ctx.truth_path)
        if st != "OK":
            raise Abort("the truth output could not be saved (%s: %s)." % (st, why))
        if not qinfo["qout_ok"]:
            ctx.truth_path.rename(ctx.truth_path.with_name(ctx.truth_path.name + ".REJECTED"))
            raise Abort("the new truth failed its sanity check (%s). It was saved as *.REJECTED, not used."
                        % qinfo["qout_problem"])
        print("  saved truth: %s  (%s bytes, md5 %s, %d rows, peak %.4g m3/s at hour %.2f, volume %.4g m3)"
              % (ctx.truth_path.name, format(ctx.truth_path.stat().st_size, ","), md5[:12], qinfo["qout_rows"],
                 qinfo["qout_peak_m3s"], qinfo["qout_tpeak_hr"], qinfo["qout_volume_m3"]))
        q1 = describe_qout(ctx.storm1_ref)
        if np.isfinite(q1["qout_tpeak_hr"]):
            print("  peak time vs the 1.0x truth: %+.1f min  (1.0x peak at hour %.2f)"
                  % ((qinfo["qout_tpeak_hr"] - q1["qout_tpeak_hr"]) * 60.0, q1["qout_tpeak_hr"]))
        prov = {
            "script": Path(__file__).name,
            "created_local": datetime.now().isoformat(timespec="seconds"),
            "what_this_is": "Synthetic truth for the %s storm leg of Series 120: tRIBS outlet hydrograph at "
                            "Ks_mult %g, f_RS_abs %g, routing at truth, optpercolation 0, under rain scaled by %g."
                            % (ctx.label, TRUTH_KS, TRUTH_F, ctx.scale),
            "storm_label": ctx.label, "rain_scale": ctx.scale,
            "truth_parameters": {"Ks_mult": TRUTH_KS, "f_RS_abs": TRUTH_F, "routing": ROUTING_TRUTH,
                                 "optpercolation": 0, "optsnow": "0 (hard-coded in the builder)"},
            "gauge_sdf": ctx.gauge_sdf_rel,
            "files_md5": forcing_md5s(ctx.forcing_plan),
            "original_1p0_truth": {"file": ctx.storm1_ref.name, "md5": md5_of(ctx.storm1_ref)},
            "truth": {"file": ctx.truth_path.name, "md5": md5, "bytes": int(ctx.truth_path.stat().st_size),
                      "run_id": truth_id, "rows": qinfo["qout_rows"], "peak_m3s": qinfo["qout_peak_m3s"],
                      "tpeak_hr": qinfo["qout_tpeak_hr"], "volume_m3": qinfo["qout_volume_m3"]},
            "direction_check": gate,
            "verify": {"passed": None},
            "python": sys.version.split()[0], "pandas": pd.__version__, "numpy": np.__version__,
            "pytRIBS": pytribs_version(),
        }
        save_truth_prov(ctx, prov)
        if not ctx.args.no_cleanup:
            cleanup_raw_results(truth_id, ctx.calib_dir, ctx.category)

    ctx.truth_md5 = prov["truth"]["md5"]

    print("\n" + "=" * 92)
    print("STEP 3 of 5 -- SELF-CONSISTENCY CHECK")
    print("=" * 92)
    passed, info = run_selfcheck(ctx)
    prov["verify"] = info
    save_truth_prov(ctx, prov)

    print("\n" + "=" * 92)
    print("STEP 4 of 5 -- COMPARE WITH THE SERIES 110 REFERENCE (information only)")
    print("=" * 92)
    prov["series110_crosscheck"] = crosscheck_series110(ctx)
    save_truth_prov(ctx, prov)

    print("\n" + "=" * 92)
    print("STEP 5 of 5 -- RESULT")
    print("=" * 92)
    if passed:
        print("TRUTH FOR %s: READY. Self-consistency PASSED (PBIAS %+.4f%%, KGE_2012 %.4f)."
              % (ctx.label.upper(), info["pbias_pct"], info["kge_2012"]))
        print("  truth file : %s" % ctx.truth_path)
        print("  provenance : %s" % ctx.prov_path)
        print("  Next: python -u %s --rain_scale %g --limit 3 2>&1 | tee series120_%s_trial.log"
              % (Path(__file__).name, ctx.scale, ctx.label))
        return 0
    print("TRUTH FOR %s: NOT READY. The self-consistency check FAILED (%s). Do not start the sweep."
          % (ctx.label.upper(), info.get("reason") or "PBIAS %+.4f%%, KGE_2012 %.4f"
             % (info.get("pbias_pct", float("nan")), info.get("kge_2012", float("nan")))))
    print("  verify_truth_provenance_110.py asks the same question directly. Send me this log.")
    return 1


def stage_smoke_test(ctx):
    problems, prov = preflight(ctx, require_verified=False)
    if problems:
        raise Abort("; ".join(problems))
    passed, info = run_selfcheck(ctx)
    print("\n" + "=" * 92)
    if passed:
        print("SMOKE TEST: PASS. Truth file, scorer and truth_file setting agree.")
        code = 0
    else:
        print("SMOKE TEST: FAIL. Do not start the sweep.")
        code = 1
    print("(This does not change the provenance record; --make_truth is what records a result.)")
    print("=" * 92)
    return code


def preflight(ctx, require_verified=True):
    """Everything that must still be true before a sweep. Returns (problems, prov)."""
    problems = []
    prov = load_truth_prov(ctx)
    if prov is None or "truth" not in prov:
        return ["no provenance record for %s (%s). Run --make_truth for this storm first."
                % (ctx.label, ctx.prov_path.name)], None
    if require_verified and prov.get("verify", {}).get("passed") is not True:
        problems.append("the provenance record does not show a PASSED self-consistency check. "
                        "Run --make_truth (it repeats the check on an existing truth).")
    if not ctx.truth_path.exists():
        problems.append("the truth file %s is missing" % ctx.truth_path.name)
    else:
        md5 = md5_of(ctx.truth_path)
        if md5 != prov["truth"].get("md5"):
            problems.append("the truth file's checksum (%s) does not match its record (%s): the file changed "
                            "since it was made" % (md5[:12], str(prov["truth"].get("md5"))[:12]))
        qinfo = describe_qout(ctx.truth_path)
        if not qinfo["qout_ok"]:
            problems.append("the truth file fails its sanity check: %s" % qinfo["qout_problem"])
        ctx.truth_md5 = md5
    # the scaled rain files must be exactly what this script would write from the current originals
    for name, state in forcing_status(ctx.forcing_plan):
        if not state.startswith("already present"):
            problems.append("rain file %s: %s" % (name, state))
    rec = prov.get("files_md5", {})
    now = forcing_md5s(ctx.forcing_plan)
    for name, md5 in rec.items():
        if name in now and now[name] != md5:
            problems.append("rain file %s has changed since the truth was made (checksum differs)" % name)
    # the 1.0x truth is only the reference for the direction check; a change there is worth a warning
    o = prov.get("original_1p0_truth", {})
    if o.get("md5") and ctx.storm1_ref.exists() and md5_of(ctx.storm1_ref) != o["md5"]:
        print("  WARNING: the 1.0x truth in synth_truth/ has changed since this storm's truth was made "
              "(checksum differs). The sweep is unaffected (it scores against this storm's own truth), but the "
              "recorded direction-check numbers refer to the old file.")
    return problems, prov


# ======================================================================
# Pairing with the 1.0x results (read-only)
# ======================================================================
def find_series120(summary_dir, explicit=None):
    p = Path(explicit) if explicit else summary_dir / SERIES120_CSV_NAME
    if not p.exists():
        return None, p
    try:
        df = pd.read_csv(p)
    except Exception as e:
        print("  Warning: could not read %s (%s)." % (p.name, e))
        return None, p
    if any(c not in df.columns for c in PARAM_COLS):
        print("  Warning: %s lacks %s; cannot pair." % (p.name, PARAM_COLS))
        return None, p
    return df, p


def match_to_ref(points_df, ref_df, rtol=MATCH_RTOL):
    """For each draw, positional index of the 1.0x row with the same (Ks_mult, f_RS_abs)
    within rtol, or -1."""
    out = np.full(len(points_df), -1, dtype=int)
    if ref_df is None:
        return out
    ref = ref_df[PARAM_COLS].to_numpy(dtype=float)
    pts = points_df[PARAM_COLS].to_numpy(dtype=float)
    for i, row in enumerate(pts):
        ok = np.all(np.isclose(ref, row, rtol=rtol, atol=0.0), axis=1)
        hits = np.flatnonzero(ok)
        if hits.size:
            out[i] = hits[0]
    return out


def pairing_check(ctx, samples, ref_df, ref_path):
    """Returns the array of twin indices. Raises Abort on poor pairing with default n and seed."""
    n, seed = ctx.args.n, ctx.args.seed
    if ref_df is None:
        print("  PAIRING CHECK vs the 1.0x results: %s not found in summary_tables/ (pass --series120_csv PATH). "
              "The sweep itself is unaffected; only the paired comparison in the report is skipped." % ref_path.name)
        return np.full(len(samples), -1, dtype=int)
    twin = match_to_ref(samples, ref_df)
    n_matched = int((twin >= 0).sum())
    print("  PAIRING CHECK vs the 1.0x results (%s, %d rows): %d/%d draws matched a twin."
          % (ref_path.name, len(ref_df), n_matched, len(samples)))
    if n_matched < 0.9 * min(len(samples), len(ref_df)):
        if n == DEFAULT_N and seed == DEFAULT_SEED and not ctx.args.allow_unpaired:
            raise Abort("poor pairing with the 1.0x results using the default n and seed: the box, the sampler "
                        "or the dict order here differs from the 1.0x script. Stopping before spending hours. "
                        "(--allow_unpaired overrides.)")
        print("  (Expected if n or seed differs from the 1.0x leg: the draws are new points, not twins.)")
    return twin


# ======================================================================
# End-of-run report
# ======================================================================
def enclosure_table(ks, f, pb, n_slices=N_F_SLICES_ENCL):
    """Is the zero-bias ridge inside the sampled Ks range at each f? Model-free. For each f slice
    take the k runs with the lowest Ks and the k with the highest Ks and look at their mean PBIAS.
    ENCLOSED = the two ends have opposite signs (PBIAS crosses zero somewhere inside the box).
    The ridge is NOT enclosed when both ends have the same sign: then the zero crossing, if any, lies
    outside the Ks range sampled. Which side is read from how PBIAS moves with Ks in that slice."""
    lo, hi = LHS_PARAMS["f_RS_abs"]["lo"], LHS_PARAMS["f_RS_abs"]["hi"]
    edges = np.logspace(np.log10(lo), np.log10(hi), n_slices + 1)
    rows = []
    for b in range(n_slices):
        last = (b == n_slices - 1)
        m = (f >= edges[b]) & ((f <= edges[b + 1]) if last else (f < edges[b + 1]))
        row = {"f_lo": edges[b], "f_hi": edges[b + 1], "runs": int(m.sum()),
               "low_end": np.nan, "high_end": np.nan, "verdict": "too few runs", "kind": "few"}
        if m.sum() >= MIN_RUNS_PER_SLICE:
            k_s, p_s = ks[m], pb[m]
            k = max(3, int(round(0.1 * m.sum())))
            order = np.argsort(k_s)
            low, high = float(p_s[order[:k]].mean()), float(p_s[order[-k:]].mean())
            row.update(low_end=low, high_end=high)
            rho = float(np.corrcoef(_rank(k_s), _rank(p_s))[0, 1])
            if low * high < 0:
                near = min(abs(low), abs(high)) < FEASIBLE_PBIAS_PCT
                row["verdict"] = ("enclosed (but one end is inside +/-%g%%: ridge near the box edge)"
                                  % FEASIBLE_PBIAS_PCT) if near else "enclosed"
                row["kind"] = "enclosed_edge" if near else "enclosed"
            else:
                sign = 1.0 if (low + high) > 0 else -1.0       # +: over-predicts at both ends
                # flow normally FALLS as Ks rises (rho < 0): over-prediction everywhere -> need MORE Ks
                above = (sign > 0) == (rho < 0)
                row["verdict"] = ("NOT enclosed: PBIAS %s at both ends, ridge ABOVE the Ks upper bound"
                                  if above else "NOT enclosed: PBIAS %s at both ends, ridge BELOW the Ks lower bound"
                                  ) % ("positive" if sign > 0 else "negative")
                row["kind"] = "above" if above else "below"
                if abs(rho) < 0.2:
                    row["verdict"] += " (weak Ks trend in this slice: side uncertain)"
        rows.append(row)
    return rows


def print_enclosure(ctx, rows):
    ks_lo, ks_hi = LHS_PARAMS["Ks_mult"]["lo"], LHS_PARAMS["Ks_mult"]["hi"]
    print("\n  RIDGE ENCLOSURE in Ks (box Ks %g-%g). Mean PBIAS (%%) of the lowest-Ks and highest-Ks runs in each f slice:"
          % (ks_lo, ks_hi))
    print("  %17s  %5s  %10s  %10s  verdict" % ("f slice", "runs", "low-Ks end", "high-Ks end"))
    for r in rows:
        le = "%+10.1f" % r["low_end"] if np.isfinite(r["low_end"]) else "%10s" % "-"
        he = "%+10.1f" % r["high_end"] if np.isfinite(r["high_end"]) else "%10s" % "-"
        print("  %7.4f-%7.4f  %5d  %s  %s  %s" % (r["f_lo"], r["f_hi"], r["runs"], le, he, r["verdict"]))
    above = [r for r in rows if r["kind"] == "above"]
    below = [r for r in rows if r["kind"] == "below"]
    few = [r for r in rows if r["kind"] == "few"]
    edge = [r for r in rows if r["kind"] == "enclosed_edge"]
    print("  Summary (%s):" % ctx.label)
    if not above and not below and not few:
        print("    Ridge is ENCLOSED in Ks at every f slice%s."
              % (" (near the box edge in %d slice(s))" % len(edge) if edge else ""))
    if above:
        print("    Ks UPPER BOUND (%g) CHECK: the ridge is above the box for f from %.4f to %.4f (%d slice(s))."
              % (ks_hi, min(r["f_lo"] for r in above), max(r["f_hi"] for r in above), len(above)))
    else:
        print("    Ks upper bound (%g) check: ridge not above the box at any f slice." % ks_hi)
    if below:
        print("    Ks LOWER BOUND (%g) CHECK: the ridge is below the box for f from %.4f to %.4f (%d slice(s))."
              % (ks_lo, min(r["f_lo"] for r in below), max(r["f_hi"] for r in below), len(below)))
    else:
        print("    Ks lower bound (%g) check: ridge not below the box at any f slice." % ks_lo)
    if few:
        print("    %d slice(s) had fewer than %d runs, so the check was not made there." % (len(few), MIN_RUNS_PER_SLICE))
    if above or below:
        print("    -> A finding about this storm as much as a sampling gap. Whether to extend the Ks range (a "
              "separate new-seed batch; it would have no twins in the 1.0x leg) is Tracy's decision.")


def report_final(ctx, final_df, ref_df, failed_rows, n_planned, n_full):
    print("\n" + "=" * 70)
    print("SERIES %s REPORT (storm leg %s = %s, channel loss OFF)" % (SERIES, ctx.leg, ctx.label))
    print("=" * 70)
    if "pbias_pct" not in final_df.columns:
        print("  No pbias_pct column in results; skipping report.")
        return
    pb = pd.to_numeric(final_df["pbias_pct"], errors="coerce").to_numpy(float)
    good = np.isfinite(pb)
    n_good = int(good.sum())
    print("  Scored runs: %d of %d planned; %d HANG/FAILED/mismatch record(s) in the FAILED log (see FAILED csv)."
          % (n_good, n_planned, len(failed_rows)))
    if n_good < 3:
        return
    for p in PARAM_COLS:
        v = pd.to_numeric(final_df[p], errors="coerce")
        print("  %s coverage: %.4f - %.4f" % (p, v.min(), v.max()))

    n_feas = int((np.abs(pb[good]) < FEASIBLE_PBIAS_PCT).sum())
    print("\n  Feasible (|PBIAS| < %g%%): %d/%d = %.1f%%" % (FEASIBLE_PBIAS_PCT, n_feas, n_good, 100.0 * n_feas / n_good))
    complete = n_good >= 0.9 * n_full
    if not complete:
        print("  -> partial sweep (%d of %d draws scored): the extension threshold of %d only applies to a "
              "complete sweep." % (n_good, n_full, MIN_FEASIBLE_PER_STORM))
    elif n_feas < MIN_FEASIBLE_PER_STORM:
        print("  -> BELOW the extension threshold of %d. Add a new-seed batch (never a larger n)." % MIN_FEASIBLE_PER_STORM)
    else:
        print("  -> meets the extension threshold of %d." % MIN_FEASIBLE_PER_STORM)

    ks_all = pd.to_numeric(final_df["Ks_mult"], errors="coerce").to_numpy(float)
    f_all = pd.to_numeric(final_df["f_RS_abs"], errors="coerce").to_numpy(float)
    lo, hi = LHS_PARAMS["f_RS_abs"]["lo"], LHS_PARAMS["f_RS_abs"]["hi"]
    edges = np.logspace(np.log10(lo), np.log10(hi), N_F_BANDS + 1)
    print("\n  f-band check (PBIAS = 100*sum(sim-obs)/sum(obs); <0 under-predicts, >0 over-predicts)")
    print("  %17s  %5s  %8s  %8s  %8s  note" % ("f band", "runs", "PBIAS<0", "PBIAS>0", "feasible"))
    band_gaps = 0
    for b in range(N_F_BANDS):
        last = (b == N_F_BANDS - 1)
        m = good & (f_all >= edges[b]) & ((f_all <= edges[b + 1]) if last else (f_all < edges[b + 1]))
        n_b, n_neg, n_pos = int(m.sum()), int((pb[m] < 0).sum()), int((pb[m] > 0).sum())
        n_fe = int((np.abs(pb[m]) < FEASIBLE_PBIAS_PCT).sum())
        note = ""
        if n_b == 0:
            note = "EMPTY"
            band_gaps += 1
        elif n_neg == 0 or n_pos == 0:
            note = "no bracket (all one side of zero bias)"
            band_gaps += 1
        print("  %7.4f-%7.4f  %5d  %8d  %8d  %8d  %s" % (edges[b], edges[b + 1], n_b, n_neg, n_pos, n_fe, note))
    print("  (A band with no bracket can mean the ridge leaves the Ks box there -- a finding to look at, "
          "not necessarily a sampling gap.)")

    if complete:
        print_enclosure(ctx, enclosure_table(ks_all[good], f_all[good], pb[good]))
    else:
        print("\n  (ridge-enclosure check skipped: it needs a complete sweep)")

    if "kge_2012" in final_df.columns:
        best = final_df.sort_values("kge_2012", ascending=False).iloc[0]
        print("\n  Best run (reference only -- the best run wanders along an equifinal ridge, it is not the estimate):")
        print("    Ks_mult=%.3f  f_RS_abs=%.4f  KGE_2012=%.3f  PBIAS=%+.2f%%"
              % (best["Ks_mult"], best["f_RS_abs"], best["kge_2012"], best["pbias_pct"]))

    if "caution_zone" in final_df.columns:
        cz = final_df[final_df["caution_zone"].astype(str) == "True"]
        print("\n  Caution-zone runs (Ks >= %g, f >= %g): %d" % (CAUTION_KS_LO, CAUTION_F_LO, len(cz)))
        for _, r in cz.head(15).iterrows():
            print("    Ks=%.3f  f=%.4f  KGE_2012=%.3f  PBIAS=%+.1f%%"
                  % (r["Ks_mult"], r["f_RS_abs"], _num(r.get("kge_2012")), _num(r.get("pbias_pct"))))

    # --- paired with the 1.0x leg (same build, same draws, a different storm and its own truth) ---
    if ref_df is None or "pbias_pct" not in ref_df.columns:
        print("\n  (1.0x results not available -- no paired comparison.)")
        return
    twin = match_to_ref(final_df, ref_df)
    ok = (twin >= 0) & good
    n_ok = int(ok.sum())
    print("\n  PAIRED WITH THE 1.0x LEG (same tRIBS 6.0.0 build, same draws): %d/%d points have a twin" % (n_ok, n_good))
    if n_ok < 3:
        return
    pb_s = pb[ok]
    pb_1 = pd.to_numeric(ref_df["pbias_pct"], errors="coerce").to_numpy(float)[twin[ok]]
    fin = np.isfinite(pb_1)
    pb_s, pb_1 = pb_s[fin], pb_1[fin]
    fe_s, fe_1 = np.abs(pb_s) < FEASIBLE_PBIAS_PCT, np.abs(pb_1) < FEASIBLE_PBIAS_PCT
    rho = np.corrcoef(_rank(pb_s), _rank(pb_1))[0, 1]
    print("    Feasible on the same %d points:  this storm %d   1.0x %d   BOTH %d   this storm only %d   1.0x only %d"
          % (len(pb_s), fe_s.sum(), fe_1.sum(), (fe_s & fe_1).sum(), (fe_s & ~fe_1).sum(), (~fe_s & fe_1).sum()))
    print("    PBIAS shift (this storm minus 1.0x, each against its own truth): median %+.2f, mean %+.2f pct-pts"
          % (np.median(pb_s - pb_1), np.mean(pb_s - pb_1)))
    print("    Rank correlation of PBIAS, this storm vs 1.0x: %.3f" % rho)
    print("    (The number feasible in BOTH is how much one storm's good runs overlap the other's. The full "
          "three-storm overlap needs both legs; that is a separate analysis.)")
    print("    (The old 5.3.1 storm legs used n=200 and a narrower box, so they do not pair with these runs. "
          "Compare with them on the fitted ridge, not run by run.)")


# ======================================================================
# Sweep
# ======================================================================
def load_existing_results(out_path):
    if out_path.exists():
        try:
            df = pd.read_csv(out_path)
            print("  Loaded existing results: %d rows from %s" % (len(df), out_path.name))
            return df
        except Exception as e:
            print("  Warning: could not load existing results (%s). Starting fresh." % e)
    return pd.DataFrame()


def upsert_failed(rows, new_row):
    """Add a FAILED-log record, replacing any earlier record for the same run_id."""
    return [r for r in rows if r.get("run_id") != new_row["run_id"]] + [new_row]


def run_sweep(ctx, samples, todo, twin_idx, ref_df):
    args = ctx.args
    variant = variant_suffix(args.n, args.seed)
    out_path = ctx.summary_dir / ("lhs_results_synth_Ks_f_%s_%s%s.csv" % (SERIES, ctx.label, variant))
    failed_path = ctx.summary_dir / ("lhs_results_synth_Ks_f_FAILED_%s_%s%s.csv" % (SERIES, ctx.label, variant))
    assert "_100" not in out_path.name, "output would collide with an old Series 100 file"
    assert out_path.name != SERIES120_CSV_NAME, "output would collide with the 1.0x results"

    existing_df = load_existing_results(out_path)
    existing_ids = set(existing_df["run_id"].values) if not existing_df.empty else set()
    results = existing_df.to_dict("records") if not existing_df.empty else []
    if not existing_df.empty and not args.skip_existing:
        print("  NOTE: existing results found and --skip_existing not given: every draw will be re-run "
              "(rows are replaced by run_id). Use --skip_existing to resume.")

    failed_log_rows = []
    if args.skip_existing and failed_path.exists():
        try:
            failed_log_rows = pd.read_csv(failed_path).to_dict("records")
            print("  Loaded %d earlier HANG/FAILED record(s) from %s" % (len(failed_log_rows), failed_path.name))
        except Exception as e:
            print("  Warning: could not read %s (%s); starting a new log." % (failed_path.name, e))
    n_total = len(todo)

    print("\n" + "=" * 70)
    print("Series %s -- Ks_mult x f_RS_abs vs SYNTHETIC TRUTH (%s), storm leg %s (%s, rain x%g), channel loss OFF"
          % (SERIES, ctx.truth_path.name, ctx.leg, ctx.label, ctx.scale))
    print("  %d of %d draws, seed=%d, timeout=%ds, cleanup=%s"
          % (n_total, len(samples), args.seed, args.timeout, "off" if args.no_cleanup else "on"))
    for p, b in LHS_PARAMS.items():
        print("  %s: %g - %g  (%s)" % (p, b["lo"], b["hi"], b["scale"]))
    print("  PINNED (truth): kinemvelcoef=%g  flowexp r=%g  channelroughness=%g  optpercolation=0"
          % (ROUTING_TRUTH["kinemvelcoef"], ROUTING_TRUTH["flowexp"], ROUTING_TRUTH["channelroughness"]))
    print("  Rain: %s   Truth: %s [md5 %s]" % (Path(ctx.gauge_sdf_rel).name, ctx.truth_path.name, ctx.truth_md5[:12]))
    print("  Output: %s" % out_path.name)
    print("=" * 70 + "\n", flush=True)

    completed = skipped = excluded = hung = failed = cleaned = caution_hits = 0
    aborted = None
    sweep_start = time.time()

    for i, row in todo.iterrows():
        ks, f_val, flagged = float(row["Ks_mult"]), float(row["f_RS_abs"]), bool(row["caution"])
        tag = sweep_tag(f_val, args.seed, ctx.short)
        run_id = predicted_run_id(ks, tag)

        if is_truth_point(ks, f_val):
            print("[%3d/%d]  EXCLUDED: draw equals the truth point." % (i + 1, len(samples)))
            excluded += 1
            continue

        print("\n[%3d/%d]  Ks=%.3fx  f=%.5f  -> %s%s"
              % (i + 1, len(samples), ks, f_val, run_id, "  [CAUTION ZONE]" if flagged else ""), flush=True)
        if flagged:
            caution_hits += 1

        extras = {"series": SERIES, "storm_leg": ctx.leg, "storm_label": ctx.label, "rain_scale": ctx.scale,
                  "truth_file": ctx.truth_path.name, "truth_md5": ctx.truth_md5, "seed": args.seed,
                  "draw_idx": int(i), "caution_zone": flagged,
                  "swept_param": "lhs_synth_Ks_f_%s_%s" % (SERIES, ctx.label)}

        csv_ok = (ctx.csv_dir / ("%s_compare_obs_sim.csv" % run_id)).exists()
        metrics_file = ctx.summary_dir / ("%s_metrics_summary.csv" % run_id)
        # Skip only a draw that is complete: its compare CSV exists AND its scored row is on record (in the
        # results file or as a metrics file). A run cut off between the two files is simply re-run.
        if args.skip_existing and csv_ok and (run_id in existing_ids or metrics_file.exists()):
            print("  SKIP (already complete): %s" % run_id)
            skipped += 1
            if metrics_file.exists() and run_id not in existing_ids:
                try:
                    m = pd.read_csv(metrics_file).iloc[0].to_dict()
                    m.update(extras)
                    results.append(m)
                    existing_ids.add(run_id)
                except Exception:
                    pass
            continue

        res = run_one(ctx, run_id, ks, f_val, tag, ctx.truth_rel, ctx.truth_path)
        status, reason = res["status"], res["reason"]

        if res.get("abort"):
            print("\n" + "!" * 70)
            print("STOPPING THE SWEEP (%s)" % status)
            print(res["abort"])
            print("Results so far are saved.")
            print("!" * 70)
            failed_log_rows = upsert_failed(failed_log_rows, {
                "run_id": run_id, "status": status, "reason": reason, "elapsed_min": res["elapsed"] / 60,
                "draw_idx": int(i), "Ks_mult": ks, "f_RS_abs": f_val})
            pd.DataFrame(failed_log_rows).to_csv(failed_path, index=False)
            failed += 1
            aborted = res["abort"]
            break

        if status == "SUCCESS":
            metrics = res["metrics"]
            metrics.update(extras)
            results = [r for r in results if r.get("run_id") != run_id]
            results.append(metrics)
            existing_ids.add(run_id)
            completed += 1
            note_outcome(ctx, True)
            if any(r.get("run_id") == run_id for r in failed_log_rows):
                failed_log_rows = [r for r in failed_log_rows if r.get("run_id") != run_id]
                pd.DataFrame(failed_log_rows, columns=FAILED_COLS).to_csv(failed_path, index=False)
            print("  %s (tRIBS 6.0):  KGE_2012=%.3f  PBIAS=%+.1f%%" % (ctx.label, _num(metrics.get("kge_2012")),
                                                                       _num(metrics.get("pbias_pct"))))
            tw = twin_idx[i] if i < len(twin_idx) else -1
            if ref_df is not None and tw >= 0:
                o = ref_df.iloc[tw]
                print("  1.0x twin:        KGE_2012=%.3f  PBIAS=%+.1f%%" % (_num(o.get("kge_2012")), _num(o.get("pbias_pct"))))
            if not args.no_cleanup and cleanup_raw_results(run_id, ctx.calib_dir, ctx.category):
                cleaned += 1
        else:
            print("  %s: %s  (%s)" % (status, run_id, reason))
            if status == "HANG":
                hung += 1
            else:
                failed += 1
            failed_log_rows = upsert_failed(failed_log_rows, {
                "run_id": run_id, "status": status, "reason": reason, "elapsed_min": res["elapsed"] / 60,
                "draw_idx": int(i), "Ks_mult": ks, "f_RS_abs": f_val})
            pd.DataFrame(failed_log_rows).to_csv(failed_path, index=False)
            note_outcome(ctx, False)

        attempted = completed + hung + failed
        remaining = n_total - attempted - skipped - excluded
        if attempted > 0:
            avg = (time.time() - sweep_start) / attempted
            eta_min = avg * max(remaining, 0) / 60
            print("  Run time: %.1f min  |  ETA: %.0f min (%.1f h) remaining"
                  % (res["elapsed"] / 60, eta_min, eta_min / 60), flush=True)

        if results:
            pd.DataFrame(results).to_csv(out_path, index=False)      # save after every run

    print("\nSweep %s: %d ran, %d skipped, %d excluded, %d hung, %d failed, %d raw folders cleaned, "
          "%d in the caution zone"
          % ("STOPPED EARLY" if aborted else "complete", completed, skipped, excluded, hung, failed, cleaned, caution_hits))

    if results:
        final_df = pd.DataFrame(results)
        if "kge_2012" in final_df.columns:
            final_df = final_df.sort_values("kge_2012", ascending=False)
        final_df.to_csv(out_path, index=False)
        print("Saved: %s  (%d rows)" % (out_path.name, len(final_df)))
        try:
            report_final(ctx, final_df, ref_df, failed_log_rows, n_total, len(samples))
        except Exception as e:
            print("\n(report skipped: %s -- results are saved regardless)" % e)
    else:
        print("No results to save.")
    if aborted:
        return 3
    return 0


# ======================================================================
# Main
# ======================================================================
class Ctx(object):
    pass


def print_dry_run(ctx, samples, ref_df, ref_path, competing):
    truth_id, ttag, vid, vtag = truth_ids(ctx)
    variant = variant_suffix(ctx.args.n, ctx.args.seed)
    print("\n" + "=" * 92)
    print("DRY RUN -- %s (rain x%g, storm leg %s). Nothing is written, built or run." % (ctx.label, ctx.scale, ctx.leg))
    print("=" * 92)
    print_forcing_plan(ctx.forcing_plan)
    print("\n  scaled rain files, state right now:")
    for name, state in forcing_status(ctx.forcing_plan):
        print("    %-36s %s" % (name, state))
    print("\n  1.0x truth (the reference for the direction check, read-only): %s  [md5 %s]"
          % (ctx.storm1_ref.name, md5_of(ctx.storm1_ref)[:12]))
    print("  This storm's truth will be: %s" % ctx.truth_path)
    print("    truth run_id %s" % truth_id)
    print("    self-check run_id %s" % vid)
    if ctx.truth_path.exists() or ctx.prov_path.exists():
        problems, prov = preflight(ctx, require_verified=True)
        if prov is None:
            print("    STATE: a truth file exists but there is no readable provenance record.")
        elif problems:
            print("    STATE: a truth exists but the sweep would NOT start:")
            for p in problems:
                print("      - %s" % p)
        else:
            v = prov.get("verify", {})
            print("    STATE: truth READY. Checksum OK, rain files OK, self-consistency PASSED "
                  "(PBIAS %+.4f%%, KGE_2012 %.4f)." % (_num(v.get("pbias_pct")), _num(v.get("kge_2012"))))
    else:
        print("    STATE: no truth yet. Next command: --make_truth")
    d110 = ctx.calib_dir / ("synth_truth_ccoff_%s" % ctx.label)
    n110 = len(list(d110.glob("*_Outlet.qout"))) if d110.exists() else 0
    print("  Series 110 cc-OFF reference for this storm: %s" % (("found (%d file)" % n110) if n110 else "not found"))
    print("\n  SWEEP PLAN: %d draws, seed %d, Ks %g-%g (linear), f %g-%g (log)"
          % (len(samples), ctx.args.seed, LHS_PARAMS["Ks_mult"]["lo"], LHS_PARAMS["Ks_mult"]["hi"],
             LHS_PARAMS["f_RS_abs"]["lo"], LHS_PARAMS["f_RS_abs"]["hi"]))
    for p in PARAM_COLS:
        print("    %s: draws span %.4f - %.4f" % (p, samples[p].min(), samples[p].max()))
    print("    caution-zone draws (Ks >= %g, f >= %g): %d" % (CAUTION_KS_LO, CAUTION_F_LO, int(samples["caution"].sum())))
    print("    first draw: Ks=%.4f  f=%.5f -> %s"
          % (samples.iloc[0]["Ks_mult"], samples.iloc[0]["f_RS_abs"],
             predicted_run_id(float(samples.iloc[0]["Ks_mult"]), sweep_tag(float(samples.iloc[0]["f_RS_abs"]),
                                                                            ctx.args.seed, ctx.short))))
    print("    output: lhs_results_synth_Ks_f_%s_%s%s.csv" % (SERIES, ctx.label, variant))
    print("    time: truth + self-check about 3 minutes; sweep about %.1f h (at ~0.9 min per run)"
          % (len(samples) * 0.9 / 60.0))
    if competing:
        print("\n  NOTE: another tRIBS build/run script appears to be running right now:")
        for pid, cmd in competing:
            print("    pid %d: %s" % (pid, cmd[:110]))
        print("  (A dry run is harmless, but --make_truth and the sweep would refuse to start.)")
    print("\n--dry_run: nothing written, built or run.")


def main():
    parser = argparse.ArgumentParser(
        description="Series 120 storm legs (0.8x and 1.25x): scaled rain, this storm's synthetic truth, "
                    "and the 400-draw Ks_mult x f_RS_abs sweep scored against it, on tRIBS 6.0.0.")
    parser.add_argument("--rain_scale", type=float, required=True,
                        help="Which storm leg: 0.8 or 1.25 (multiplier on the rainfall rate)")
    parser.add_argument("--dry_run", action="store_true",
                        help="Print the plan and the current state, then stop. Writes, builds and runs nothing")
    parser.add_argument("--make_truth", action="store_true",
                        help="Scale the rain, build this storm's truth, run the direction check and the "
                             "self-consistency check, write the provenance record")
    parser.add_argument("--remake_truth", action="store_true",
                        help="With --make_truth: rename the existing truth to *.superseded_<time> and build a new one")
    parser.add_argument("--smoke_test", action="store_true",
                        help="Only the self-consistency check, on an existing truth (about 1 minute)")
    parser.add_argument("--n", type=int, default=DEFAULT_N,
                        help="Number of draws (default %d; changing it changes EVERY draw, so they will not "
                             "pair with the 1.0x leg)" % DEFAULT_N)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED,
                        help="Random seed (default %d). A new seed gives a new batch of points -- use it for "
                             "extension batches" % DEFAULT_SEED)
    parser.add_argument("--limit", type=int, default=None,
                        help="Only process the first N draws (a short trial; resume later with --skip_existing)")
    parser.add_argument("--skip_existing", action="store_true",
                        help="Skip draws whose compare CSV already exists (resume)")
    parser.add_argument("--timeout", type=int, default=300, help="Per-run hard timeout in seconds (default 300)")
    parser.add_argument("--no_cleanup", action="store_true",
                        help="Keep raw per-run results folders instead of deleting them after a successful run")
    parser.add_argument("--series120_csv", type=str, default=None,
                        help="Path to the 1.0x results CSV (for pairing); default: look in summary_tables/")
    parser.add_argument("--allow_unpaired", action="store_true",
                        help="Continue even if the draws do not pair with the 1.0x results")
    parser.add_argument("--ignore_running_check", action="store_true",
                        help="Skip the check for other tRIBS scripts running (only if you are sure none are)")
    parser.add_argument("--overwrite_forcing", action="store_true",
                        help="Replace existing scaled rain files that differ from what this script writes")
    parser.add_argument("--max_consecutive_failures", type=int, default=8,
                        help="Stop if this many runs in a row fail (default 8)")
    parser.add_argument("--src_sdf", default=SRC_SDF_DEFAULT,
                        help="Original rain gauge SDF, relative to lab/ (default %s)" % SRC_SDF_DEFAULT)
    args = parser.parse_args()

    if sum([args.dry_run, args.make_truth, args.smoke_test]) > 1:
        parser.error("use only one of --dry_run, --make_truth, --smoke_test")
    if args.remake_truth and not args.make_truth:
        parser.error("--remake_truth only makes sense together with --make_truth")
    if args.n < 1:
        parser.error("--n must be at least 1")
    try:
        label_, short, leg = leg_names(args.rain_scale)
    except ValueError as e:
        parser.error(str(e))

    ctx = Ctx()
    ctx.args, ctx.label, ctx.short, ctx.leg, ctx.scale = args, label_, short, leg, float(args.rain_scale)
    ctx.consec_fail = 0
    ctx.truth_md5 = ""
    cwd = Path.cwd()
    ctx.calib_dir = cwd.parent / "calibration_work"
    ctx.summary_dir = ctx.calib_dir / "03_comparisons" / "summary_tables"
    ctx.csv_dir = ctx.calib_dir / "03_comparisons" / "csv_exports"
    synth_dir = ctx.calib_dir / "synth_truth"
    ctx.truth_dir = ctx.calib_dir / ("synth_truth_120_%s" % label_)
    ctx.prov_path = ctx.truth_dir / ("PROVENANCE_truth_120_%s.json" % label_)
    ctx.category = "%s_%s" % (CATEGORY_BASE, short)
    ctx.gauge_sdf_rel = derived_name(args.src_sdf, label_)

    print("Series 120 storm leg: %s  (rain x%g)   python %s, pandas %s"
          % (label_, ctx.scale, sys.version.split()[0], pd.__version__))

    # ------------------------------------------------------------------
    # HARD SAFETY CHECKS
    # ------------------------------------------------------------------
    if not (cwd / SCORER_SCRIPT).exists():
        raise RuntimeError("%s not found in %s. Run this script from the lab/ directory." % (SCORER_SCRIPT, cwd))
    if "gauge_sdf" not in inspect.signature(builder.build_input_file).parameters:
        raise RuntimeError("build_sensitivity_run_tribs6.build_input_file has no gauge_sdf argument, so it cannot "
                           "be pointed at scaled rain. This script needs the tRIBS-6 builder that has it.")
    qout_files = list(synth_dir.glob("*.qout")) if synth_dir.exists() else []
    if len(qout_files) != 1:
        raise RuntimeError(
            "Expected exactly one *.qout file in %s (the 1.0x truth), found %d: %s. Without it the scorer would "
            "silently fall back to real-gauge mode." % (synth_dir, len(qout_files), [f.name for f in qout_files]))
    ctx.storm1_ref = qout_files[0]
    print("1.0x truth, read-only: %s  [md5 %s]" % (ctx.storm1_ref.name, md5_of(ctx.storm1_ref)[:12]))

    synth_res = synth_dir.resolve()
    tres = ctx.truth_dir.resolve()
    if tres == synth_res or synth_res in tres.parents:
        raise RuntimeError("truth folder %s is, or is inside, %s; refusing." % (ctx.truth_dir, synth_dir))

    install_series_120_hooks(ctx.category)
    truth_id, _, _, _ = truth_ids(ctx)
    ctx.truth_path = ctx.truth_dir / ("%s_Outlet.qout" % truth_id)
    ctx.truth_rel = "%s/%s" % (ctx.truth_dir.name, ctx.truth_path.name)     # relative to calibration_work

    competing = None
    if args.ignore_running_check:
        print("Running-jobs check skipped (--ignore_running_check).")
    else:
        competing = find_competing_jobs()
        if competing is None:
            print("WARNING: could not check for other running tRIBS scripts (ps unavailable). Make sure none are "
                  "running -- they share current_run_config.json.")
        elif competing and not args.dry_run:
            print("\nANOTHER tRIBS BUILD/RUN SCRIPT APPEARS TO BE RUNNING:")
            for pid, cmd in competing:
                print("  pid %d: %s" % (pid, cmd[:110]))
            raise RuntimeError(
                "Refusing to start: these scripts share current_run_config.json with no locking, so running both "
                "can corrupt both. Wait for them to finish (or stop them with pkill -f <script>), or re-run with "
                "--ignore_running_check if you are sure they are stale.")
        elif not competing:
            print("No other tRIBS build/run script detected.")

    ctx.forcing_plan = plan_forcing(args.src_sdf, label_, ctx.scale)

    samples = generate_lhs_samples(args.n, LHS_PARAMS, seed=args.seed)
    samples["caution"] = [in_caution_zone(k, f) for k, f in zip(samples["Ks_mult"], samples["f_RS_abs"])]
    todo = samples.iloc[:args.limit] if args.limit else samples
    ref_df, ref_path = find_series120(ctx.summary_dir, args.series120_csv)

    if args.dry_run:
        pairing_check(ctx, samples, ref_df, ref_path)
        print_dry_run(ctx, samples, ref_df, ref_path, competing)
        return 0

    ctx.summary_dir.mkdir(parents=True, exist_ok=True)
    ctx.csv_dir.mkdir(parents=True, exist_ok=True)

    try:
        if args.make_truth:
            print("\nPLAN: --make_truth for %s: scaled rain, the truth run, direction check, self-consistency check, "
                  "provenance. About 3 minutes." % label_)
            print_forcing_plan(ctx.forcing_plan)
            return stage_make_truth(ctx)

        if args.smoke_test:
            return stage_smoke_test(ctx)

        # ---------------- SWEEP ----------------
        problems, prov = preflight(ctx, require_verified=True)
        if problems:
            raise Abort("the truth for %s is not ready for a sweep: %s" % (label_, "; ".join(problems)))
        print("\nPRE-FLIGHT OK: the truth file's checksum matches its record, the scaled rain files are exactly what "
              "this script writes, and the self-consistency check PASSED (PBIAS %+.4f%%, KGE_2012 %.4f)."
              % (_num(prov["verify"].get("pbias_pct")), _num(prov["verify"].get("kge_2012"))))
        twin_idx = pairing_check(ctx, samples, ref_df, ref_path)
        return run_sweep(ctx, samples, todo, twin_idx, ref_df)
    except Abort as e:
        print("\n" + "!" * 92)
        print("STOPPED: %s" % e)
        print("!" * 92)
        print("Everything finished before this point is kept. After fixing the cause, re-run the same command "
              "(with --skip_existing for a sweep).")
        return 3
    except KeyboardInterrupt:
        print("\nStopped by you. Finished runs are kept; re-run the same command with --skip_existing to continue.")
        return 130


if __name__ == "__main__":
    sys.exit(main())
