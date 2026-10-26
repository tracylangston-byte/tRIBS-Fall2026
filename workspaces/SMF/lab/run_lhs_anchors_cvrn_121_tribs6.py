"""
run_lhs_anchors_cvrn_121_tribs6.py
==================================
Series 121 -- the routing-parameter identifiability experiment ("Scenario 3") re-done
on tRIBS 6.0.0 / pytRIBS 1.0.0, on the 1.0x storm only.

WHAT THIS IS
------------
Scenario 1 (Series 120) found a ridge in the Ks_mult x f_RS_abs plane: many Ks/f pairs
give the same total runoff volume as the synthetic truth. This experiment asks what
happens to the three ROUTING parameters at points along that ridge:

    kinemvelcoef      (cv)      truth 4.5     sampled 2.5  - 6.5
    flowexp r                   truth 0.24    sampled 0.18 - 0.35
    channelroughness  (n)       truth 0.026   sampled 0.02 - 0.10

At each of five fixed "anchor" points (Ks and f held still, volume-matched), the three
routing parameters are varied together by Latin-hypercube-style sampling and every run is
scored against the 1.0x synthetic truth (Ks 7.0, f 0.012, cv 4.5, flowexp r 0.24,
n 0.026). The question: at each anchor, do cv, flowexp r and n affect the score, and
does that change from one anchor to the next?

THE FIVE ANCHORS (chosen 2026-10-25 from the fitted zero-bias ridges of all three
Series 120 storm legs, so each point also fits the 0.8x and 1.25x storm volumes)
-------------------------------------------------------------------------------
    label         Ks_mult    f_RS_abs   fitted PBIAS 0.8x / 1.0x / 1.25x   why
    f0p004         4.60       0.004      +1.7 / -0.7 / -0.8              flat low-f arm: volume says nothing about f
    f0p008         5.10       0.008      -1.5 / +0.0 / +1.0              below truth, high-KGE part of the ridge
    f0p012_true    7.00       0.012       0.0 /  0.0 /  0.0              the truth point itself
    f0p015         8.05       0.015      -0.8 / +0.9 / -0.3              above truth, high-KGE part of the ridge
    f0p02          8.80       0.020      -0.5 / +0.1 / +0.7              top of the all-storms region; hydrograph
                                                                          starts to be damped (KGE_2012 ~0.94)
The fitted values come from a thin-plate-spline surface through the 400 Series 120 runs
of each storm, not from runs at these exact points. That is why the --smoke_test below
runs every anchor once at truth routing first: it shows the volume error each anchor
really has.

THE DESIGN: 100 RUNS PER ANCHOR = 50 SHARED + 50 INDEPENDENT
------------------------------------------------------------
  * 50 "shared" draws (seed 42): the SAME 50 routing-parameter sets at every anchor. This
    makes anchors directly comparable, because any difference between anchors comes from
    the anchor, not from luck in the draws. (Same sampler, seed and ranges as old
    Series 101, so these 50 sets are the old Series 101 base draws.)
  * 50 "independent" draws: a different seed at every anchor. Old Series 101 found that
    pooling anchors that share draws is pseudoreplication (the same 50 numbers counted six
    times). The independent draws are what can be pooled honestly.
Both batches are kept in the results with a 'batch' column ('shared' or 'indep').
Seeds for the independent batches: f0p004 201, f0p008 202, f0p012_true 203, f0p015 204,
f0p02 205.

WHAT IS PINNED
--------------
Channel loss OFF (optpercolation = 0). OPTSNOW = 0 (set inside the builder, as in Series
120). Everything else as in Series 120. Storm: the original 12 Aug 2014 event (1.0x).
Scored against the 1.0x truth in calibration_work/synth_truth/ (exactly one *.qout whose
name contains 'ks7p0x'). The 0.8x / 1.25x truth files live in their own folders and are
not used here.

HOW THE RUNS ARE NAMED AND WHERE THEY GO
-----------------------------------------
Like the Series 120 scripts, this one teaches the (unedited) builder about series 121 in
memory only. Nothing on disk is changed.

    run_id        SMF_20140812_121_Ks<Ks>x_<anchor>_cv<cv>_r<r>_n<n>   e.g. ..._121_Ks4p6x_f0p004_cv5p257564_...
    raw results   calibration_work/02_results/121_lhs_anchors_cvrn/<run_id>/   (deleted after scoring)
    run inputs    calibration_work/01_run_inputs/121_lhs_anchors_cvrn/
    per anchor    calibration_work/03_comparisons/summary_tables/lhs_results_anchor_<label>_121.csv
    all anchors   calibration_work/03_comparisons/summary_tables/lhs_results_anchor_ALL_121.csv
    failed log    calibration_work/03_comparisons/summary_tables/lhs_results_anchor_FAILED_121.csv
    smoke test    calibration_work/03_comparisons/summary_tables/lhs_results_anchor_SMOKE_121.csv

The ALL file is rebuilt after every call from the per-anchor files on disk (all five
labels are named explicitly, and any that are missing are listed), so running one anchor
at a time never overwrites the others.

MODES
-----
1. --dry_run       Prints the anchors and the draws, checks that every run name is
                   unique, estimates the time. Builds nothing, runs no tRIBS.
2. --smoke_test    Builds, runs and scores each selected anchor ONCE at truth routing.
                   For the truth anchor expect PBIAS +0.0000% and KGE_2012 1.0000. For the
                   other anchors the volume error should be inside +/-2%. About 1 minute
                   per anchor. Do this before the sweep.
3. (default)       The sweep: 100 runs for each selected anchor.
4. --report_only   Runs nothing; rebuilds the ALL file and prints the report from the
                   files already on disk.

SAFETY CHECKS BUILT IN
-----------------------
  * Refuses to start if another tRIBS / sweep script is running (they all share
    calibration_work/current_run_config.json).
  * Refuses to run unless exactly one *.qout is in calibration_work/synth_truth/ and its
    name contains 'ks7p0x' (override with --allow_other_truth only if you are sure).
  * After every build, the run configuration the builder wrote is checked (right run name,
    right Ks/f/cv/flowexp r/n, channel loss off, baseline rain file, and NO truth_file
    setting, which would send the scorer to a different truth).
  * After every score, the scored row is read back: obs_mode must be 'synth', the
    scorer must report reading the 1.0x truth file, its observed peak must match that
    file's peak, and the row must carry the Ks, f, cv, flowexp r, n and optpercolation
    that were asked for. A mismatch STOPS the script (everything after it would be wrong
    in the same way).
  * Each run is a separate process with a hard time limit (default 300 s); a hung run is
    killed and logged and the script moves on. Ctrl+C also kills the run in flight.
  * 8 failures in a row stop the script (that is a systematic problem, not bad luck).
  * Raw result folders are deleted only if they sit exactly at 02_results/<category>/<run_id>.

DO NOT run any other tRIBS build/run script in lab/ while this is running.

USAGE (run from the lab/ directory, like every other script)
-------------------------------------------------------------
    python run_lhs_anchors_cvrn_121_tribs6.py --dry_run
    python run_lhs_anchors_cvrn_121_tribs6.py --smoke_test                    # all five anchors
    python run_lhs_anchors_cvrn_121_tribs6.py --anchors f0p012_true --limit 3 # short trial
    python run_lhs_anchors_cvrn_121_tribs6.py --anchors f0p012_true --skip_existing 2>&1 | tee series121_f0p012_true.log

    # Recommended: one anchor per command (about 1.6 h each), truth anchor first. Use -u
    # so output is not held back, and tee so you have a log:
    python -u run_lhs_anchors_cvrn_121_tribs6.py --anchors f0p012_true --skip_existing 2>&1 | tee -a series121_f0p012_true.log
    python -u run_lhs_anchors_cvrn_121_tribs6.py --anchors f0p008 --skip_existing 2>&1 | tee -a series121_f0p008.log
    ... and so on for f0p015, f0p004, f0p02. After an interruption, run the same command
    again; --skip_existing resumes.

    python run_lhs_anchors_cvrn_121_tribs6.py --report_only

STATUS (2026-10-25)
-------------------
Written without access to the Codespace. NOT yet run on the real builder and scorer.
Tested against stand-ins for the builder and scorer (fast fakes that mimic the files and
messages the real ones produce), all passing:
  * --dry_run: 500 run names, all unique, 7.8 h estimate.
  * --smoke_test: truth anchor PASS, other anchors OK; and the CHECK path when an anchor is
    not volume-matched.
  * A 6-run trial, then a resume with --skip_existing: 94 ran, 6 skipped, 100 rows, no repeats.
  * A hung run (killed at the time limit, no orphan process), a failed run and a build error
    all land in the FAILED log and the sweep carries on; a later success removes the record.
  * 8 failures in a row stop the script.
  * Each of these stops the script at once with a plain message: the scored row disagrees with
    what was asked, the scorer used the real gauge instead of the truth, the builder's
    configuration carries a truth_file setting or ignored an override, another tRIBS script
    is running, synth_truth/ holds no file / two files / a file with the wrong name.
  * Ctrl+C mid-run kills the run in flight, leaves no orphan and no half-written file, and
    --skip_existing resumes.
  * One anchor per command, then --report_only: the ALL file holds all five anchors.
  * The 50 shared draws are identical at every anchor and are the old Series 101 base draws
    (checked against the old sampler's code); the independent batches differ.
The first real use is therefore --dry_run, then --smoke_test.
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
from pathlib import Path

import numpy as np
import pandas as pd

import build_sensitivity_run_tribs6 as builder

# ----------------------------------------------------------------------
# Series identity
# ----------------------------------------------------------------------
SERIES        = "121"
CATEGORY      = "121_lhs_anchors_cvrn"
STORM_LEG     = "1.00x"
SCORER_SCRIPT = "run_sensitivity_single_interp_tribs6.py"
BASELINE_GAUGE_SDF = "../init_data/met/Master_Precip.sdf"

# ----------------------------------------------------------------------
# Truth (1.0x leg) and the pinned values
# ----------------------------------------------------------------------
TRUTH_KS  = 7.0
TRUTH_F   = 0.012
TRUTH_TOL = 1e-6
TRUTH_NAME_HINTS = ("ks7p0x", "ks7_f012")

ROUTING_TRUTH = {
    "kinemvelcoef":     4.5,
    "flowexp":          0.24,
    "channelroughness": 0.026,
}
ROUTING_COLS = ["kinemvelcoef", "flowexp", "channelroughness"]

# ----------------------------------------------------------------------
# Anchors (see the docstring). indep_seed = seed of the independent batch at this anchor.
# ----------------------------------------------------------------------
ANCHORS = [
    {"label": "f0p004",      "Ks_mult": 4.60, "f_RS_abs": 0.004, "indep_seed": 201,
     "note": "flat low-f arm"},
    {"label": "f0p008",      "Ks_mult": 5.10, "f_RS_abs": 0.008, "indep_seed": 202,
     "note": "below truth"},
    {"label": "f0p012_true", "Ks_mult": 7.00, "f_RS_abs": 0.012, "indep_seed": 203,
     "note": "the truth point"},
    {"label": "f0p015",      "Ks_mult": 8.05, "f_RS_abs": 0.015, "indep_seed": 204,
     "note": "above truth"},
    {"label": "f0p02",       "Ks_mult": 8.80, "f_RS_abs": 0.020, "indep_seed": 205,
     "note": "top of the all-storms region"},
]
ANCHOR_LABELS = [a["label"] for a in ANCHORS]
TRUE_ANCHOR_LABEL = "f0p012_true"

# ----------------------------------------------------------------------
# Routing-parameter sampling box. SAME as old Series 101 (ranges, order, sampler).
# ----------------------------------------------------------------------
LHS_PARAMS = {
    "kinemvelcoef":     {"lo": 2.5,  "hi": 6.5},    # truth 4.5
    "flowexp":          {"lo": 0.18, "hi": 0.35},   # truth 0.24
    "channelroughness": {"lo": 0.02, "hi": 0.10},   # truth 0.026
}
DEFAULT_N        = 100     # runs per anchor
DEFAULT_N_SHARED = 50      # of those, the shared (paired) draws
DEFAULT_SEED     = 42      # seed of the shared batch

FEASIBLE_PBIAS_PCT = 2.0
SMOKE_TRUTH_PBIAS_TOL = 0.01   # percent
SMOKE_TRUTH_KGE_MIN   = 0.9999
OBS_PEAK_RTOL  = 0.02
RTOL_ROUTING   = 1e-2          # read-back tolerance for cv / flowexp r / n (tolerates rounding)
RTOL_SOIL      = 1e-6          # read-back tolerance for Ks, f (anchor values are exact numbers)
MIN_PER_RUN_GUESS = 0.94       # minutes per run, from Series 120

FAILED_COLS = ["anchor_label", "run_id", "status", "reason", "elapsed_min", "draw_idx", "batch",
               "Ks_mult", "f_RS_abs", "kinemvelcoef", "flowexp", "channelroughness"]

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
# Sampler -- the same method as Series 91-101 (stratified-independent, linear)
# ======================================================================
def generate_lhs_samples(n, params, seed=None):
    """Each parameter range is cut into n equal slices, one uniform draw per slice, then
    shuffled independently. Same seed + same n + same ranges + same dict order =
    identical draws (the old Series 101 base draws for n=50, seed=42)."""
    if n <= 0:
        return pd.DataFrame({p: [] for p in params})
    rng = np.random.default_rng(seed)
    samples = {}
    for param, bounds in params.items():
        lo, hi = bounds["lo"], bounds["hi"]
        intervals = np.linspace(lo, hi, n + 1)
        points = rng.uniform(intervals[:-1], intervals[1:])
        rng.shuffle(points)
        samples[param] = points
    return pd.DataFrame(samples)


def build_draws(anchor, n_total, n_shared, seed):
    """The routing-parameter draws for one anchor: the shared batch, then the independent
    batch. Returns a DataFrame with the three parameters plus batch, seed and draw_idx."""
    n_shared = max(0, min(int(n_shared), int(n_total)))
    n_indep = int(n_total) - n_shared
    parts = []
    shared = generate_lhs_samples(n_shared, LHS_PARAMS, seed=seed)
    shared["batch"] = "shared"
    shared["seed"] = seed
    parts.append(shared)
    indep = generate_lhs_samples(n_indep, LHS_PARAMS, seed=anchor["indep_seed"])
    indep["batch"] = "indep"
    indep["seed"] = anchor["indep_seed"]
    parts.append(indep)
    df = pd.concat([p for p in parts if len(p)], ignore_index=True)
    df["draw_idx"] = np.arange(len(df), dtype=int)
    return df


# ======================================================================
# Teaching the (unedited) builder about series 121
# ======================================================================
_HOOKS_INSTALLED = False


def install_series_121_hooks():
    """In-memory only: tell the builder that Ks_mult runs are series 121 and that series
    121 lives in its own category folder. Idempotent. Stops the script if the result is not
    exactly what the rest of this script assumes."""
    global _HOOKS_INSTALLED
    if _HOOKS_INSTALLED:
        return
    cfg = dict(builder.PARAM_CONFIG["Ks_mult"])
    cfg["series"] = SERIES
    builder.PARAM_CONFIG["Ks_mult"] = cfg

    original = getattr(builder, "_orig_get_run_category_121", None) or builder.get_run_category
    builder._orig_get_run_category_121 = original

    def get_run_category_121(series_str):
        if str(series_str) == SERIES:
            return CATEGORY
        return original(series_str)

    builder.get_run_category = get_run_category_121

    run_id, _ = builder.build_run_id("Ks_mult", 7.0)
    expected = f"{builder.LOCATION}_{builder.EVENT_DATE}_{SERIES}_Ks7p0x"
    if run_id != expected:
        raise RuntimeError(f"Series-121 hook self-check failed: build_run_id gave {run_id!r}, expected "
                           f"{expected!r}. The builder's naming has changed; stop and look at "
                           f"build_sensitivity_run_tribs6.py before running anything.")
    cat = builder.get_run_category(builder.PARAM_CONFIG["Ks_mult"]["series"])
    if cat != CATEGORY:
        raise RuntimeError(f"Series-121 hook self-check failed: category {cat!r}, expected {CATEGORY!r}.")
    _HOOKS_INSTALLED = True


# ======================================================================
# Small helpers
# ======================================================================
def label(v):
    """The builder's own number-to-text rule (7.0 -> 7p0, 0.012 -> 0p012)."""
    return builder.value_to_label(v)


def make_tag(anchor_label, cv, r, n):
    """Run-name suffix. The anchor label already carries f (and the base run name carries Ks)."""
    return f"{anchor_label}_cv{label(cv)}_r{label(r)}_n{label(n)}"


def predicted_run_id(ks, tag):
    run_id, _ = builder.build_run_id("Ks_mult", ks)
    return f"{run_id}_{tag}" if tag else run_id


def is_truth_run(ks, f, cv, r, n):
    """True only if ALL five values equal the truth (a continuous draw will never do this)."""
    return (abs(ks - TRUTH_KS) < TRUTH_TOL and abs(f - TRUTH_F) < TRUTH_TOL and
            abs(cv - ROUTING_TRUTH["kinemvelcoef"]) < TRUTH_TOL and
            abs(r - ROUTING_TRUTH["flowexp"]) < TRUTH_TOL and
            abs(n - ROUTING_TRUTH["channelroughness"]) < TRUTH_TOL)


def pinned_overrides(ks, f, cv, r, n):
    return {
        "Ks_mult":          ks,
        "f_RS_abs":         f,
        "kinemvelcoef":     cv,
        "flowexp":          r,
        "channelroughness": n,
        "optpercolation":   0,
    }


def variant_suffix(n, n_shared, seed):
    s = ""
    if n != DEFAULT_N:
        s += f"_n{n}"
    if n_shared != DEFAULT_N_SHARED:
        s += f"_sh{n_shared}"
    if seed != DEFAULT_SEED:
        s += f"_seed{seed}"
    return s


def _num(x, default=np.nan):
    try:
        v = float(x)
        return v if np.isfinite(v) else default
    except (TypeError, ValueError):
        return default


def _close(a, b, rtol=1e-2):
    a, b = _num(a), _num(b)
    if not (np.isfinite(a) and np.isfinite(b)):
        return False
    return abs(a - b) <= rtol * abs(b) + 1e-12


def save_csv_atomic(df, path):
    """Write to a temporary file, then rename, so an interruption never leaves half a file."""
    tmp = path.with_name(path.name + ".part")
    df.to_csv(tmp, index=False)
    os.replace(tmp, path)


def load_csv(path, what="results"):
    if path.exists():
        try:
            df = pd.read_csv(path)
            print(f"  Loaded existing {what}: {len(df)} rows from {path.name}")
            return df
        except Exception as e:
            print(f"  Warning: could not load {path.name} ({e}). Starting fresh.")
    return pd.DataFrame()


def csv_already_exists(run_id, calib_dir):
    return (calib_dir / "03_comparisons" / "csv_exports" / f"{run_id}_compare_obs_sim.csv").exists()


def upsert_failed(rows, new_row):
    """Add a FAILED-log record, replacing any earlier record for the same run_id."""
    return [r for r in rows if r.get("run_id") != new_row["run_id"]] + [new_row]


def cleanup_raw_results(run_id, calib_dir):
    """Delete a successful run's raw results folder. Only deletes a folder located EXACTLY at
    02_results/<CATEGORY>/<run_id>. Never raises."""
    try:
        parent = (calib_dir / "02_results" / CATEGORY).resolve()
        raw = (parent / run_id).resolve()
        if raw.parent != parent or raw.name != run_id:
            print(f"  (cleanup skipped: {raw} is not directly inside {parent})")
            return False
        if raw.exists():
            shutil.rmtree(raw)
            return True
    except Exception as e:
        print(f"  (cleanup warning for {run_id}: {e})")
    return False


# ======================================================================
# Other jobs that share current_run_config.json
# ======================================================================
def _ancestor_pids():
    pids, pid = set(), os.getpid()
    for _ in range(32):
        if pid in pids or pid <= 1:
            break
        pids.add(pid)
        try:
            with open(f"/proc/{pid}/status") as fh:
                m = re.search(r"^PPid:\s*(\d+)", fh.read(), re.M)
            pid = int(m.group(1)) if m else 0
        except Exception:
            break
    return pids


def find_competing_jobs():
    """Other python processes running a script that shares current_run_config.json.
    Returns a list of (pid, command line), or None if the check could not be done."""
    try:
        out = subprocess.run(["ps", "-eo", "pid=,args="], capture_output=True, text=True, timeout=10).stdout
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


def check_no_competing_jobs():
    jobs = find_competing_jobs()
    if jobs is None:
        print("  (note: could not check for other running tRIBS scripts; make sure none is running)")
    elif jobs:
        lines = "\n".join(f"    pid {pid}: {args[:110]}" for pid, args in jobs)
        raise Abort("another tRIBS / sweep script appears to be running, and they all share "
                    "calibration_work/current_run_config.json, so running two at once would mix them up:\n"
                    f"{lines}\n  Wait for it to finish (or stop it), then run this again.")


# ======================================================================
# Reading the truth file (to check the scorer used it)
# ======================================================================
def read_qout(path):
    """Parse a tRIBS outlet file by hand. Returns (time_hr, q_m3s)."""
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
                    continue
            try:
                t.append(float(toks[0]))
                q.append(float(toks[1]))
            except (ValueError, IndexError):
                continue
    return np.array(t, dtype=float), np.array(q, dtype=float)


_PEAK_CACHE = {}


def expected_obs_peak(path, event_start, event_end):
    """The observed peak the scorer WILL see if it reads this truth file: same reading, same
    5-minute time-interpolation, same event window."""
    key = (str(path), event_start, event_end)
    if key not in _PEAK_CACHE:
        t, q = read_qout(path)
        idx = pd.DatetimeIndex(pd.to_datetime(t * 3600, unit="s", origin=pd.Timestamp("2014-08-01")))
        s = pd.Series(q, index=idx).resample("5min").interpolate(method="time")
        _PEAK_CACHE[key] = float(s.loc[event_start:event_end].max())
    return _PEAK_CACHE[key]


# ======================================================================
# Context: where things are
# ======================================================================
class Ctx:
    pass


def make_ctx(args):
    ctx = Ctx()
    ctx.args = args
    ctx.script_dir = Path.cwd()
    ctx.calib_dir = ctx.script_dir.parent / "calibration_work"
    ctx.summary_dir = ctx.calib_dir / "03_comparisons" / "summary_tables"
    ctx.csv_dir = ctx.calib_dir / "03_comparisons" / "csv_exports"
    ctx.summary_dir.mkdir(parents=True, exist_ok=True)
    ctx.suffix = variant_suffix(args.n, args.n_shared, args.seed)
    ctx.failed_path = ctx.summary_dir / f"lhs_results_anchor_FAILED_{SERIES}{ctx.suffix}.csv"
    ctx.all_path = ctx.summary_dir / f"lhs_results_anchor_ALL_{SERIES}{ctx.suffix}.csv"
    ctx.smoke_path = ctx.summary_dir / f"lhs_results_anchor_SMOKE_{SERIES}.csv"
    ctx.truth_path = None
    ctx.consec_fail = 0
    return ctx


def anchor_out_path(ctx, anchor_label):
    return ctx.summary_dir / f"lhs_results_anchor_{anchor_label}_{SERIES}{ctx.suffix}.csv"


def select_anchors(spec):
    if not spec:
        return list(ANCHORS)
    wanted = [s.strip() for s in spec.split(",") if s.strip()]
    bad = [w for w in wanted if w not in ANCHOR_LABELS]
    if bad:
        raise SystemExit(f"Unknown anchor label(s) {bad}. Valid labels: {ANCHOR_LABELS}")
    return [a for a in ANCHORS if a["label"] in wanted]


def check_truth_file(ctx, allow_other):
    """Exactly one *.qout in synth_truth/, with the 1.0x truth's name. Returns its path."""
    if not (ctx.script_dir / SCORER_SCRIPT).exists():
        raise Abort(f"{SCORER_SCRIPT} not found in {ctx.script_dir}. Run this script from the lab/ directory.")
    synth_dir = ctx.calib_dir / "synth_truth"
    qouts = sorted(synth_dir.glob("*.qout")) if synth_dir.exists() else []
    if len(qouts) != 1:
        raise Abort(f"Expected exactly one *.qout file in {synth_dir} to activate synthetic-truth mode, "
                    f"found {len(qouts)}: {[q.name for q in qouts]}. Without it the scorer silently "
                    f"falls back to the real gauge.")
    name = qouts[0].name
    print(f"\nSynthetic truth mode confirmed active: {name}")
    ok = any(h in name.lower() for h in TRUTH_NAME_HINTS)
    if ok:
        print("  The name looks like the Ks 7.0 / f 0.012 (1.0x) truth. (A name is only a hint; the "
              "smoke test is what proves the file's contents.)")
    elif allow_other:
        print(f"  NOTE: the name contains none of {', '.join(repr(h) for h in TRUTH_NAME_HINTS)}, but "
              f"--allow_other_truth was given, so continuing.")
    else:
        raise Abort(f"The truth file is {name!r}, whose name contains none of "
                    f"{', '.join(repr(h) for h in TRUTH_NAME_HINTS)}. This experiment uses the 1.0x "
                    f"truth (Ks 7.0, f 0.012). If you are sure this file is that truth, re-run with "
                    f"--allow_other_truth.")
    return qouts[0]


# ======================================================================
# Build, score, verify ONE run
# ======================================================================
def run_scorer(timeout_sec):
    """Run the scorer in its own process group so a hang can be killed as a unit.
    Returns (returncode, elapsed_sec, timed_out, tribs_warning, stdout)."""
    proc = subprocess.Popen([sys.executable, SCORER_SCRIPT], cwd=Path.cwd(), start_new_session=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    t0 = time.time()

    def kill():
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        except (ProcessLookupError, PermissionError):
            pass

    try:
        stdout, _ = proc.communicate(timeout=timeout_sec)
        timed_out = False
    except subprocess.TimeoutExpired:
        kill()
        stdout, _ = proc.communicate()
        timed_out = True
    except KeyboardInterrupt:
        kill()
        proc.communicate()
        print("\nInterrupted: the running tRIBS job was stopped. Finished runs are kept; re-run the "
              "same command with --skip_existing to continue.")
        raise
    elapsed = time.time() - t0
    if stdout:
        print(stdout, end="" if stdout.endswith("\n") else "\n")
    tribs_warning = bool(stdout) and "WARNING: tRIBS may have failed" in stdout
    returncode = None if timed_out else proc.returncode
    return returncode, elapsed, timed_out, tribs_warning, stdout or ""


def check_config(cfg, run_id, ks, f, cv, r, n):
    """What the builder wrote into current_run_config.json versus what was asked for."""
    problems = []
    if cfg.get("run_id") != run_id:
        problems.append(f"run_id is {cfg.get('run_id')!r}, expected {run_id!r}")
    wanted = [("Ks_mult", ks), ("f_RS_abs", f), ("kinemvelcoef", cv), ("flowexp", r),
              ("channelroughness", n)]
    for key, want in wanted:
        if key not in cfg:
            problems.append(f"{key} missing from the run configuration")
        elif not _close(cfg[key], want, 1e-3):
            problems.append(f"{key} is {cfg[key]!r}, expected {want!r}")
    if "optpercolation" not in cfg or _num(cfg["optpercolation"]) != 0:
        problems.append(f"optpercolation is {cfg.get('optpercolation')!r}, expected 0 (channel loss off)")
    if cfg.get("run_category") not in (None, CATEGORY):
        problems.append(f"run_category is {cfg.get('run_category')!r}, expected {CATEGORY!r}")
    if cfg.get("gauge_sdf") and Path(str(cfg["gauge_sdf"])).name != Path(BASELINE_GAUGE_SDF).name:
        problems.append(f"gauge_sdf is {cfg['gauge_sdf']!r}, expected the baseline {BASELINE_GAUGE_SDF!r} "
                        f"(the 1.0x storm; a scaled-storm rain file would be a different storm)")
    if cfg.get("truth_file"):
        problems.append(f"the configuration carries truth_file={cfg['truth_file']!r}, which would send the "
                        f"scorer to a different truth than the 1.0x truth in synth_truth/")
    return problems


def prepare_run(ctx, ks, f, cv, r, n, tag, pred_id):
    try:
        run_id, _, _ = builder.build_input_file("Ks_mult", ks, overrides=pinned_overrides(ks, f, cv, r, n),
                                                tag=tag)
    except Exception as e:
        raise BuildFailed(str(e)[:300])
    if run_id != pred_id:
        raise Abort(f"builder produced run_id {run_id}, expected {pred_id}")
    cfg_path = ctx.calib_dir / "current_run_config.json"
    try:
        cfg = json.loads(cfg_path.read_text())
    except Exception as e:
        raise Abort(f"cannot read {cfg_path} after the build: {e}")
    problems = check_config(cfg, run_id, ks, f, cv, r, n)
    if problems:
        raise Abort("the builder's run configuration does not match what was asked for: " + "; ".join(problems))
    return run_id, cfg


def check_observed(ctx, metrics, cfg, stdout):
    """Did the scorer really use the 1.0x truth file as 'observed'? Returns a list of problems."""
    problems = []
    if str(metrics.get("obs_mode")) != "synth":
        problems.append(f"obs_mode is {metrics.get('obs_mode')!r}, expected 'synth' "
                        f"(the scorer fell back to the real gauge)")
    m = re.search(r"Reading observed from:\s*(.+)", stdout or "")
    if m and m.group(1).strip() != ctx.truth_path.name:
        problems.append(f"the scorer read {m.group(1).strip()!r} as observed, expected {ctx.truth_path.name!r}")
    elif not m:
        print("  (note: the scorer printed no 'Reading observed from:' line, so that check was skipped; "
              "the observed-peak check below still applies)")
    try:
        exp = expected_obs_peak(ctx.truth_path, cfg["event_start"], cfg["event_end"])
        got = float(metrics["obs_peak_m3s"])
        if not (abs(got - exp) <= OBS_PEAK_RTOL * max(abs(exp), 1e-12)):
            problems.append(f"the observed peak the scorer used ({got:.4g} m3/s) is not the peak of "
                            f"{ctx.truth_path.name} ({exp:.4g} m3/s)")
    except Exception as e:
        print(f"  (note: could not double-check the observed peak: {str(e)[:100]})")
    return problems


def readback_problems(m, ks, f, cv, r, n):
    """Compare the scored metrics row against what this script asked for."""
    probs = []
    expected = [("Ks_mult", ks, RTOL_SOIL), ("f_RS_abs", f, RTOL_SOIL),
                ("kinemvelcoef", cv, RTOL_ROUTING), ("flowexp", r, RTOL_ROUTING),
                ("channelroughness", n, RTOL_ROUTING), ("optpercolation", 0, RTOL_SOIL)]
    for key, want, rtol in expected:
        got = _num(m.get(key))
        if not np.isfinite(got) or abs(got - want) > rtol * abs(want) + 1e-12:
            probs.append(f"{key}: scored row says {m.get(key)!r}, asked for {want!r}")
    return probs


def run_one(ctx, ks, f, cv, r, n, tag):
    """Build, score and verify one run. Returns a dict: status is SUCCESS, BUILD_FAILED, HANG or
    FAILED (log and carry on), or OBSERVED_MISMATCH / READBACK_MISMATCH (carries 'abort': stop).
    Raises Abort for a build-time problem that would affect every run."""
    pred_id = predicted_run_id(ks, tag)
    out = {"run_id": pred_id, "status": None, "reason": "", "elapsed": 0.0, "metrics": None}
    try:
        run_id, cfg = prepare_run(ctx, ks, f, cv, r, n, tag, pred_id)
    except BuildFailed as e:
        out.update(status="BUILD_FAILED", reason=str(e))
        return out

    # Remove this run's old outputs so a stale file can never be mistaken for a new one.
    (ctx.csv_dir / f"{run_id}_compare_obs_sim.csv").unlink(missing_ok=True)
    (ctx.summary_dir / f"{run_id}_metrics_summary.csv").unlink(missing_ok=True)

    returncode, elapsed, timed_out, tribs_warning, stdout = run_scorer(ctx.args.timeout)
    out["elapsed"] = elapsed
    if timed_out:
        out.update(status="HANG", reason="wall-clock timeout")
        return out
    if returncode != 0 or tribs_warning:
        out.update(status="FAILED",
                   reason="tRIBS reported non-zero exit" if tribs_warning else f"scorer exited {returncode}")
        return out
    metrics_file = ctx.summary_dir / f"{run_id}_metrics_summary.csv"
    if not metrics_file.exists():
        out.update(status="FAILED", reason=f"no metrics file {metrics_file.name}")
        return out
    try:
        metrics = pd.read_csv(metrics_file).iloc[0].to_dict()
    except Exception as e:
        out.update(status="FAILED", reason=f"metrics file unreadable: {e}")
        return out
    out["metrics"] = metrics

    problems = check_observed(ctx, metrics, cfg, stdout)
    if problems:
        out["abort"] = (f"scoring did not use the intended 'observed' file for {run_id}: "
                        f"{'; '.join(problems)}. Every later run would be scored the same wrong way, so "
                        f"the script is stopping.")
        out.update(status="OBSERVED_MISMATCH", reason="; ".join(problems))
        return out
    probs = readback_problems(metrics, ks, f, cv, r, n)
    if probs:
        out["abort"] = (f"the scored row for {run_id} does not match what this script asked for: "
                        f"{'; '.join(probs)}. The override plumbing is broken, so the remaining runs "
                        f"would be wasted.")
        out.update(status="READBACK_MISMATCH", reason="; ".join(probs))
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
        raise Abort(f"{ctx.consec_fail} runs in a row failed. That looks like a systematic problem, not "
                    f"bad luck. Look at the last error above, fix it, and re-run the same command with "
                    f"--skip_existing.")


# ======================================================================
# Reports
# ======================================================================
def corr_line(df, y="kge_2012"):
    parts = []
    for p in ROUTING_COLS:
        if p in df.columns and y in df.columns and len(df) >= 5:
            x = pd.to_numeric(df[p], errors="coerce")
            v = pd.to_numeric(df[y], errors="coerce")
            ok = x.notna() & v.notna()
            if ok.sum() >= 5 and x[ok].nunique() > 1 and v[ok].nunique() > 1:
                parts.append(f"{p}={np.corrcoef(x[ok], v[ok])[0, 1]:+.2f}")
                continue
        parts.append(f"{p}=n/a")
    return "  ".join(parts)


def report_anchor(anchor, df):
    print(f"\n  --- {anchor['label']}  (Ks {anchor['Ks_mult']}, f {anchor['f_RS_abs']}): "
          f"{len(df)} scored runs ---")
    if "kge_2012" not in df.columns or len(df) == 0:
        print("  (no kge_2012 column or no runs)")
        return
    k = pd.to_numeric(df["kge_2012"], errors="coerce")
    best = df.loc[k.idxmax()] if k.notna().any() else None
    print(f"  KGE_2012: min {k.min():.3f}   median {k.median():.3f}   max {k.max():.3f}")
    if best is not None:
        print(f"  best run: KGE_2012 {best['kge_2012']:.3f} at cv {best['kinemvelcoef']:.3f}, "
              f"flowexp r {best['flowexp']:.3f}, n {best['channelroughness']:.4f}  "
              f"(truth: cv 4.5, flowexp r 0.24, n 0.026)")
    print(f"  Pearson r with KGE_2012, all runs:  {corr_line(df)}")
    if "batch" in df.columns:
        for b in ("shared", "indep"):
            sub = df[df["batch"] == b]
            if len(sub) >= 5:
                print(f"    {b:>6} batch ({len(sub):3d} runs):  {corr_line(sub)}")


def rebuild_combined(ctx, announce=True):
    """Rebuild the ALL file from the per-anchor files on disk. All anchor labels are named
    explicitly; anchors with no file yet are listed as missing."""
    frames, rows = [], []
    for a in ANCHORS:
        p = anchor_out_path(ctx, a["label"])
        if p.exists():
            try:
                df = pd.read_csv(p)
            except Exception as e:
                rows.append((a["label"], f"UNREADABLE ({e})"))
                continue
            frames.append(df)
            rows.append((a["label"], f"{len(df)} runs"))
        else:
            rows.append((a["label"], "no file yet"))
    if announce:
        print(f"\n{'=' * 70}\nCombined file {ctx.all_path.name} is built from these anchor files:")
        for lab, what in rows:
            print(f"    {lab:<12} {what}")
    if not frames:
        print("  No anchor results on disk, so no combined file was written.")
        return None
    allf = pd.concat(frames, ignore_index=True)
    save_csv_atomic(allf, ctx.all_path)
    if announce:
        print(f"  Saved {ctx.all_path.name}  ({len(allf)} rows)")
    return allf


def report_all(ctx, allf, failed_rows):
    if allf is None or len(allf) == 0:
        return
    print(f"\n{'=' * 70}\nSERIES {SERIES} REPORT (1.0x storm, channel loss OFF)\n{'=' * 70}")
    print(f"  {'anchor':<12} {'runs':>5} {'KGE_2012 min':>13} {'median':>8} {'max':>7}   Pearson r with KGE_2012")
    for a in ANCHORS:
        sub = allf[allf["anchor_label"] == a["label"]] if "anchor_label" in allf.columns else allf.iloc[0:0]
        if len(sub) == 0:
            print(f"  {a['label']:<12} {'0':>5}   (not run yet)")
            continue
        k = pd.to_numeric(sub["kge_2012"], errors="coerce")
        print(f"  {a['label']:<12} {len(sub):>5} {k.min():>13.3f} {k.median():>8.3f} {k.max():>7.3f}   "
              f"{corr_line(sub)}")
    if failed_rows:
        st = pd.DataFrame(failed_rows)["status"].value_counts().to_dict()
        print(f"\n  {len(failed_rows)} HANG/FAILED/mismatch record(s) in {ctx.failed_path.name}: {st}")
    else:
        print("\n  No HANG/FAILED records.")
    print("  (Correlations here are a quick look. The full analysis script uses the ALL file.)")


# ======================================================================
# Dry run
# ======================================================================
def run_dry(ctx, anchors, args):
    install_series_121_hooks()
    print("\nDRY RUN -- nothing is built or run.")
    total, ids = 0, []
    for a in anchors:
        draws = build_draws(a, args.n, args.n_shared, args.seed)
        for _, d in draws.iterrows():
            tag = make_tag(a["label"], d["kinemvelcoef"], d["flowexp"], d["channelroughness"])
            ids.append(predicted_run_id(a["Ks_mult"], tag))
        total += len(draws)
    all_labels = [a["label"] for a in ANCHORS]
    dup = len(ids) - len(set(ids))
    print(f"\n  Anchors this call: {[a['label'] for a in anchors]}   (all five: {all_labels})")
    print(f"  Runs per anchor: {args.n}  ({args.n_shared} shared, seed {args.seed};  "
          f"{args.n - args.n_shared} independent, seeds per anchor)")
    print(f"  Total runs this call: {total}   run-name collisions: {dup}"
          f"{'   <-- PROBLEM' if dup else '   (all names unique)'}")
    est = total * MIN_PER_RUN_GUESS
    print(f"  Estimated time at {MIN_PER_RUN_GUESS} min/run: {est:.0f} min ({est / 60:.1f} h)")
    print(f"\n  {'anchor':<12} {'Ks':>6} {'f':>7} {'indep seed':>11} {'already on disk':>16}")
    for a in anchors:
        have = 0
        p = anchor_out_path(ctx, a["label"])
        if p.exists():
            try:
                have = len(pd.read_csv(p))
            except Exception:
                have = -1
        print(f"  {a['label']:<12} {a['Ks_mult']:>6.2f} {a['f_RS_abs']:>7.3f} {a['indep_seed']:>11} {have:>16}")
    a0 = anchors[0]
    d0 = build_draws(a0, args.n, args.n_shared, args.seed)
    print(f"\n  Draw ranges (all anchors share the first {args.n_shared}):")
    for p in ROUTING_COLS:
        print(f"    {p:<17} {d0[p].min():.4f} - {d0[p].max():.4f}   (box {LHS_PARAMS[p]['lo']} - "
              f"{LHS_PARAMS[p]['hi']}, truth {ROUTING_TRUTH[p]})")
    print("  First 3 shared draws:")
    for _, d in d0.head(3).iterrows():
        print(f"    [{int(d['draw_idx']) + 1:>3}] cv={d['kinemvelcoef']:.4f}  flowexp r={d['flowexp']:.4f}  "
              f"n={d['channelroughness']:.5f}")
    print(f"  Example run name: {ids[0] if ids else '(none)'}")
    qdir = ctx.calib_dir / "synth_truth"
    qs = [q.name for q in sorted(qdir.glob('*.qout'))] if qdir.exists() else []
    print(f"  Truth file(s) in synth_truth/ right now: {qs or 'none'}")
    return 1 if dup else 0


# ======================================================================
# Smoke test: every anchor once, at truth routing
# ======================================================================
def run_smoke(ctx, anchors, args):
    install_series_121_hooks()
    check_no_competing_jobs()
    ctx.truth_path = check_truth_file(ctx, args.allow_other_truth)
    cv, r, n = (ROUTING_TRUTH[k] for k in ROUTING_COLS)
    print(f"\n{'=' * 70}\nSMOKE TEST -- each anchor once, at truth routing "
          f"(cv {cv}, flowexp r {r}, n {n})\n{'=' * 70}")
    rows, all_ok = [], True
    for a in anchors:
        ks, f = a["Ks_mult"], a["f_RS_abs"]
        tag = f"{a['label']}_smoke"
        print(f"\n--- {a['label']}: Ks {ks}, f {f} ---")
        res = run_one(ctx, ks, f, cv, r, n, tag)
        rid = res["run_id"]
        if res.get("abort"):
            print(f"\nSMOKE TEST STOPPED: {res['abort']}")
            return 1
        if res["status"] != "SUCCESS":
            print(f"  {res['status']}: {res['reason']}")
            rows.append({"anchor_label": a["label"], "Ks_mult": ks, "f_RS_abs": f, "result": res["status"]})
            all_ok = False
            continue
        m = res["metrics"]
        pbias, kge = _num(m.get("pbias_pct")), _num(m.get("kge_2012"))
        is_true = a["label"] == TRUE_ANCHOR_LABEL
        if is_true:
            ok = abs(pbias) < SMOKE_TRUTH_PBIAS_TOL and kge > SMOKE_TRUTH_KGE_MIN
            verdict = "PASS" if ok else "FAIL (the truth point did not reproduce the truth file)"
        else:
            ok = abs(pbias) < FEASIBLE_PBIAS_PCT
            verdict = "OK" if ok else f"CHECK (|PBIAS| is not under {FEASIBLE_PBIAS_PCT:g}%: this anchor is not volume-matched)"
        all_ok = all_ok and ok
        print(f"  PBIAS {pbias:+.4f}%   KGE_2012 {kge:.4f}   tRIBS+scoring {res['elapsed'] / 60:.2f} min   -> {verdict}")
        rows.append({"anchor_label": a["label"], "Ks_mult": ks, "f_RS_abs": f, "pbias_pct": pbias,
                     "kge_2012": kge, "kge_r": _num(m.get("kge_r")), "kge_gamma": _num(m.get("kge_gamma")),
                     "peak_error_pct": _num(m.get("peak_error_pct")), "minutes": res["elapsed"] / 60,
                     "result": verdict.split(" ")[0]})
        if not args.no_cleanup:
            cleanup_raw_results(rid, ctx.calib_dir)
    if rows:
        save_csv_atomic(pd.DataFrame(rows), ctx.smoke_path)
        print(f"\nSaved {ctx.smoke_path.name}")
    mins = [x["minutes"] for x in rows if "minutes" in x]
    if mins:
        per = float(np.mean(mins))
        print(f"Average {per:.2f} min per run  ->  100 runs is about {100 * per / 60:.1f} h per anchor, "
              f"{500 * per / 60:.1f} h for all five.")
    print(f"\n{'=' * 70}")
    if all_ok and rows:
        print("SMOKE TEST: PASS. Next:  python -u run_lhs_anchors_cvrn_121_tribs6.py --anchors f0p012_true "
              "--limit 3     (optional short trial)")
        code = 0
    else:
        print("SMOKE TEST: not clean. Do not start the sweep. Send me the lines above.")
        code = 1
    print("=" * 70)
    return code


# ======================================================================
# The sweep
# ======================================================================
def run_sweep(ctx, anchors, args):
    install_series_121_hooks()
    check_no_competing_jobs()
    ctx.truth_path = check_truth_file(ctx, args.allow_other_truth)

    failed_rows = []
    if ctx.failed_path.exists():
        try:
            failed_rows = pd.read_csv(ctx.failed_path).to_dict("records")
            print(f"  Loaded {len(failed_rows)} earlier HANG/FAILED record(s) from {ctx.failed_path.name}")
        except Exception as e:
            print(f"  Warning: could not read {ctx.failed_path.name} ({e}); starting a new log.")

    plan = {}
    for a in anchors:
        d = build_draws(a, args.n, args.n_shared, args.seed)
        plan[a["label"]] = d.iloc[:args.limit] if args.limit else d
    n_planned = sum(len(d) for d in plan.values())

    print(f"\n{'=' * 70}")
    print(f"Series {SERIES} -- cv / flowexp r / n at fixed Ks-f anchors, scored against the 1.0x truth "
          f"({ctx.truth_path.name})")
    print(f"  Anchors this call: {[a['label'] for a in anchors]}   runs planned: {n_planned}   "
          f"timeout {args.timeout}s   cleanup {'off' if args.no_cleanup else 'on'}")
    print(f"  {args.n} runs per anchor = {args.n_shared} shared (seed {args.seed}) + "
          f"{args.n - args.n_shared} independent")
    for p in ROUTING_COLS:
        print(f"  {p}: {LHS_PARAMS[p]['lo']} - {LHS_PARAMS[p]['hi']}   (truth {ROUTING_TRUTH[p]})")
    print("  channel loss OFF (optpercolation 0), OPTSNOW 0 (set inside the builder)")
    print(f"{'=' * 70}\n")

    done_total = attempted_total = 0
    sweep_start = time.time()
    aborted = False
    summary = {}

    try:
        for a in anchors:
            label_, ks, f = a["label"], a["Ks_mult"], a["f_RS_abs"]
            draws = plan[label_]
            out_path = anchor_out_path(ctx, label_)
            existing = load_csv(out_path, what=f"{label_} results")
            results = existing.to_dict("records") if not existing.empty else []
            existing_ids = {r_.get("run_id") for r_ in results}
            if results and not args.skip_existing:
                print(f"  NOTE: existing {label_} results found and --skip_existing not given: every draw "
                      f"will be re-run (rows are replaced by run_id). Use --skip_existing to resume.")

            completed = skipped = excluded = hung = failed = cleaned = 0
            n_this = len(draws)
            print(f"\n--- Anchor '{label_}'  (Ks {ks}, f {f}, {n_this} runs) ---")

            for _, d in draws.iterrows():
                i = int(d["draw_idx"])
                cv, r_, n_ = float(d["kinemvelcoef"]), float(d["flowexp"]), float(d["channelroughness"])
                tag = make_tag(label_, cv, r_, n_)
                run_id = predicted_run_id(ks, tag)

                if is_truth_run(ks, f, cv, r_, n_):
                    print(f"[{label_} {i + 1:>3}/{n_this}]  EXCLUDED: all five values equal the truth.")
                    excluded += 1
                    continue

                print(f"\n[{label_} {i + 1:>3}/{n_this}]  ({d['batch']})  cv={cv:.3f}  flowexp r={r_:.3f}  "
                      f"n={n_:.4f}\n         -> {run_id}")

                extra = {"series": SERIES, "storm_leg": STORM_LEG, "anchor_label": label_,
                         "anchor_Ks": ks, "anchor_f": f, "batch": d["batch"], "seed": int(d["seed"]),
                         "draw_idx": i, "swept_param": f"lhs_anchor_{label_}"}

                # ---- resume ----
                if args.skip_existing and csv_already_exists(run_id, ctx.calib_dir):
                    print(f"  SKIP (compare CSV exists): {run_id}")
                    skipped += 1
                    mf = ctx.summary_dir / f"{run_id}_metrics_summary.csv"
                    if mf.exists() and run_id not in existing_ids:
                        try:
                            m = pd.read_csv(mf).iloc[0].to_dict()
                            m.update(extra)
                            results.append(m)
                            existing_ids.add(run_id)
                        except Exception:
                            pass
                    continue

                t_run = time.time()
                res = run_one(ctx, ks, f, cv, r_, n_, tag)
                attempted_total += 1

                row_log = {"anchor_label": label_, "run_id": run_id, "status": res["status"],
                           "reason": res["reason"], "elapsed_min": res["elapsed"] / 60 if res["elapsed"] else
                           (time.time() - t_run) / 60, "draw_idx": i, "batch": d["batch"],
                           "Ks_mult": ks, "f_RS_abs": f, "kinemvelcoef": cv, "flowexp": r_,
                           "channelroughness": n_}

                if res.get("abort"):
                    print("\n" + "!" * 70)
                    print("STOPPING THE SWEEP")
                    print(res["abort"])
                    print("Results so far are saved.")
                    print("!" * 70)
                    failed_rows = upsert_failed(failed_rows, row_log)
                    save_csv_atomic(pd.DataFrame(failed_rows, columns=FAILED_COLS), ctx.failed_path)
                    failed += 1
                    aborted = True
                    break

                if res["status"] == "SUCCESS":
                    metrics = res["metrics"]
                    metrics.update(extra)
                    results = [x for x in results if x.get("run_id") != run_id]
                    results.append(metrics)
                    existing_ids.add(run_id)
                    completed += 1
                    done_total += 1
                    if any(x.get("run_id") == run_id for x in failed_rows):
                        failed_rows = [x for x in failed_rows if x.get("run_id") != run_id]
                        save_csv_atomic(pd.DataFrame(failed_rows, columns=FAILED_COLS), ctx.failed_path)
                    print(f"  KGE_2012={_num(metrics.get('kge_2012')):.3f}  "
                          f"PBIAS={_num(metrics.get('pbias_pct')):+.1f}%")
                    if not args.no_cleanup and cleanup_raw_results(run_id, ctx.calib_dir):
                        cleaned += 1
                    note_outcome(ctx, True)
                else:
                    print(f"  {res['status']}: {run_id}  ({res['reason']})")
                    if res["status"] == "HANG":
                        hung += 1
                    else:
                        failed += 1
                    failed_rows = upsert_failed(failed_rows, row_log)
                    save_csv_atomic(pd.DataFrame(failed_rows, columns=FAILED_COLS), ctx.failed_path)
                    note_outcome(ctx, False)

                left = n_planned - attempted_total - skipped - excluded
                if attempted_total:
                    avg = (time.time() - sweep_start) / attempted_total
                    eta = avg * max(left, 0) / 60
                    print(f"  Run time: {(time.time() - t_run) / 60:.1f} min  |  ETA this call: "
                          f"{eta:.0f} min ({eta / 60:.1f} h)")

                if results:
                    df_save = pd.DataFrame(results)
                    if "draw_idx" in df_save.columns:
                        df_save = df_save.sort_values("draw_idx")
                    save_csv_atomic(df_save, out_path)

            print(f"\nAnchor '{label_}' {'STOPPED' if aborted else 'complete'}: {completed} ran, "
                  f"{skipped} skipped, {excluded} excluded, {hung} hung, {failed} failed, {cleaned} raw folders cleaned")
            if results:
                df_final = pd.DataFrame(results)
                if "draw_idx" in df_final.columns:
                    df_final = df_final.sort_values("draw_idx")
                save_csv_atomic(df_final, out_path)
                print(f"  Saved: {out_path.name}  ({len(df_final)} rows)")
                summary[label_] = df_final
                try:
                    report_anchor(a, df_final)
                except Exception as e:
                    print(f"  (anchor report skipped: {e} -- results are saved regardless)")
            if aborted:
                break
    except Abort as e:
        aborted = True
        print("\n" + "!" * 70 + f"\nSTOPPING: {e}\nResults so far are saved.\n" + "!" * 70)

    allf = rebuild_combined(ctx)
    try:
        report_all(ctx, allf, failed_rows)
    except Exception as e:
        print(f"\n(report skipped: {e} -- results are saved regardless)")
    return 1 if aborted else 0


# ======================================================================
# Main
# ======================================================================
def main():
    ap = argparse.ArgumentParser(
        description="Series 121 -- cv / flowexp r / n at five Ks-f anchors, 1.0x storm, tRIBS 6.0.0.")
    ap.add_argument("--anchors", type=str, default=None,
                    help=f"Comma-separated anchor labels to run this call (default: all). Valid: {', '.join(ANCHOR_LABELS)}")
    ap.add_argument("--n", type=int, default=DEFAULT_N, help=f"Runs per anchor (default {DEFAULT_N})")
    ap.add_argument("--n_shared", type=int, default=DEFAULT_N_SHARED,
                    help=f"Of those, how many are the shared (paired) draws (default {DEFAULT_N_SHARED})")
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED,
                    help=f"Seed of the shared batch (default {DEFAULT_SEED})")
    ap.add_argument("--dry_run", action="store_true", help="Print the plan and stop. Builds nothing.")
    ap.add_argument("--smoke_test", action="store_true",
                    help="Run each selected anchor ONCE at truth routing and stop")
    ap.add_argument("--report_only", action="store_true",
                    help="Run nothing: rebuild the ALL file and print the report from files on disk")
    ap.add_argument("--limit", type=int, default=None,
                    help="Only the first N draws at each anchor (short trial; resume with --skip_existing)")
    ap.add_argument("--skip_existing", action="store_true", help="Skip draws whose compare CSV exists (resume)")
    ap.add_argument("--timeout", type=int, default=300, help="Per-run hard timeout in seconds (default 300)")
    ap.add_argument("--no_cleanup", action="store_true",
                    help="Keep raw per-run results folders instead of deleting them after a successful run")
    ap.add_argument("--max_consecutive_failures", type=int, default=8,
                    help="Stop after this many failed runs in a row (default 8)")
    ap.add_argument("--allow_other_truth", action="store_true",
                    help="Allow a truth file whose name lacks 'ks7p0x' (only if you are sure it is the 1.0x truth)")
    args = ap.parse_args()

    if args.n < 1 or args.n_shared < 0 or args.n_shared > args.n:
        raise SystemExit("--n must be at least 1 and --n_shared must be between 0 and --n.")
    anchors = select_anchors(args.anchors)
    ctx = make_ctx(args)

    if args.dry_run:
        return run_dry(ctx, anchors, args)
    if args.report_only:
        install_series_121_hooks()
        failed = []
        if ctx.failed_path.exists():
            try:
                failed = pd.read_csv(ctx.failed_path).to_dict("records")
            except Exception:
                pass
        allf = rebuild_combined(ctx)
        report_all(ctx, allf, failed)
        if allf is not None:
            for a in ANCHORS:
                sub = allf[allf["anchor_label"] == a["label"]] if "anchor_label" in allf.columns else allf.iloc[0:0]
                if len(sub):
                    report_anchor(a, sub)
        return 0
    try:
        if args.smoke_test:
            return run_smoke(ctx, anchors, args)
        return run_sweep(ctx, anchors, args)
    except Abort as e:
        print("\n" + "!" * 70 + f"\nSTOPPED: {e}\n" + "!" * 70)
        return 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\nInterrupted (Ctrl+C). The run in flight was killed; every finished run is already saved. "
              "Resume with the same command plus --skip_existing.")
        sys.exit(130)
