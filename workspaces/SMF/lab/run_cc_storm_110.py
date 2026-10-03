"""
run_cc_storm_110.py
===================
Series 110 -- a SECOND (and, with --rain_scale 1.25, a THIRD) STORM, channel loss ON.

WHY THIS EXISTS (plain language)
----------------------------------
Everything in Series 110 so far used ONE storm: the August 12 2014 event as
recorded. A skeptic will ask whether the cc story (is channel conductivity
recoverable next to Ks and f?) is just a property of that one storm. A smaller
storm is the interesting case: when less water arrives, a larger share of it can
disappear into the channel bed, so cc may matter MORE.

This script makes the new storm by multiplying the rainfall (only the rainfall
rate, mm/hr) in both rain gauge files, SMF and SMPHQ, by a constant. Dates and
times are not touched, so the storm keeps its timing, its duration and the known
offset between the two gauges. Only its size changes. (The whole rainfall record
is scaled, including any rain before the event -- the same way the earlier
Series 100 storm siblings, storm080 and storm125, were made.)

WHAT IT DOES, IN ORDER
------------------------
  1. FORCING   writes scaled copies of the two rain files and of Master_Precip.sdf
               next to the originals (the originals are never changed):
                   ../init_data/met/precip_<station>_storm080.mdf   (one per gauge)
                   ../init_data/met/Master_Precip_storm080.sdf
               and checks them against the originals.
  2. REFERENCE builds ONE cc-OFF truth for the new storm at the standard Ks and f
               (7.0 / 0.012, routing at truth). It plays the role the file in
               synth_truth/ plays for the original storm. It is compared with the
               original storm's reference as a GATE: a smaller storm must give a
               smaller peak AND volume, a larger one bigger; if not, the scaled
               rainfall was not really used and the script stops here, before
               spending hours.
  3. TRUTHS    builds the 27 cc-ON truths (Ks 5/7/9 x f 0.007/0.012/0.020 x
               cc 95.1/201/425), the same ones as run_truth_location_110.py, but
               under the new storm.
  4. SWEEP     runs the same 250 cc-ON points as Stage 2 (LHS, seed 42; Ks 4-10.5,
               f 0.004-0.030 log, cc 30-1000 log; routing pinned at truth) under
               the new storm. Every simulation is stored (compare CSV), so any
               truth can be applied afterwards by re-scoring -- exactly as before.
  5. (OPTIONAL, run separately with --control) the same 250 points with channel
               loss OFF, under the new storm.

THE ONE THING THAT MATTERS FOR SAFETY (read this)
---------------------------------------------------
The scorer decides what "observed" is: by default the single file in
calibration_work/synth_truth/ -- which is the ORIGINAL storm's reference. For the
new storm that would be wrong. The tRIBS-6 builder has no switch for this, so this
script writes "truth_file" into current_run_config.json after each build (the
scorer already honors it) and then CHECKS, from the scorer's own output, that it
really scored against the new storm's reference. If it did not, the script stops
at once (every later run would be wrong too).
Other hard stops, all meant for when you are away: the builder's configuration
does not match what was asked for; the scaled rain file is not named in the .in
file; the reference fails its gate; 8 runs in a row fail.
Nothing in synth_truth/ is ever written. The new reference and truths live in their
own folders, and the script refuses a folder that is, or sits inside, synth_truth/.

WHAT YOU GET   (label = storm080 for 0.8x, storm125 for 1.25x)
--------------------------------------------------------------
  ../init_data/met/precip_*_<label>.mdf, Master_Precip_<label>.sdf   the scaled rain
  calibration_work/
      synth_truth_ccoff_<label>/<run_id>_Outlet.qout          the new storm's reference
      synth_truth_location_<label>/<run_id>_Outlet.qout        the 27 truths
      synth_truth_location_<label>/PROVENANCE_truth_location_110.json
      forcing_<label>_PROVENANCE.json                          checksums + rain totals
  calibration_work/03_comparisons/summary_tables/
      storm_ref_110_<label>.csv                    the reference (1 row)
      truth_location_110_<label>.csv               one row per truth (same columns as the original)
      lhs_results_joint_Ks_f_cc_<label>_110.csv    the 250-run cc-ON sweep
      lhs_results_joint_Ks_f_cc_CONTROL_<label>_110.csv   (only with --control)
      *_FAILED_*.csv                               only if something failed
  calibration_work/03_comparisons/csv_exports/<run_id>_compare_obs_sim.csv   every run's hydrograph

Run ids carry a storm tag (..._st080), so nothing from the original storm can be
overwritten:  reference ..._ref_st080,  truths ..._K5p0_f0p007_st080,
sweep ..._sw_st080,  control ..._ctl_st080.

TIME   about 1 minute per tRIBS run: reference + 27 truths ~ 30 minutes, the 250-run
sweep ~ 4-4.5 hours more. The control sweep is another ~4-4.5 hours.

USAGE (run from the lab/ directory; use the SAME python you used for Stage 2)
--------------------------------------------------------------------------------
    python run_cc_storm_110.py --dry_run                  # look, write and run nothing
    python -u run_cc_storm_110.py --skip_existing 2>&1 | tee -a storm080_110.log
                                                          # everything for the 0.8x storm
    # the same command is also how you RESUME after any interruption.

    python -u run_cc_storm_110.py --rain_scale 1.25 --skip_existing 2>&1 | tee -a storm125_110.log
                                                          # the third storm
    python -u run_cc_storm_110.py --control --skip_existing 2>&1 | tee -a storm080_ctl_110.log
                                                          # optional: channel loss OFF, same 250 points
    python run_cc_storm_110.py --only truths --skip_existing   # reference + truths, no sweep
    python run_cc_storm_110.py --n 12 --allow_unpaired ...     # a short test sweep (a different n does not
                                                               # pair with Stage 2, hence --allow_unpaired)

Never run two tRIBS build/run scripts at once: they share current_run_config.json
with no locking (this script refuses to start if it sees another one).
Overnight / unattended: keep the Codespaces tab open and always use
`python -u ... 2>&1 | tee file`.
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

# ------------------------------------------------------------------
# Truth values and design -- the SAME ones used in Stage 1, Stage 2, the control
# and the truth-location test. Passed explicitly on every build.
# ------------------------------------------------------------------
STANDARD_KS = 7.0
STANDARD_F = 0.012
ROUTING_TRUTH = {
    "kinemvelcoef":     4.5,
    "flowexp":          0.24,
    "channelroughness": 0.026,
}
CC_PARAM = "channelconductivity_mmhr"
SCORER_SCRIPT = "run_sensitivity_single_interp_tribs6.py"
REF_CC = 201.0          # written into the cc-OFF reference's .in file; INERT at optpercolation=0

# MUST match run_cc_stage2_joint_lhs_110.py exactly (ranges, scales, AND dict order:
# the order sets how the random generator is consumed) or the points will not pair.
LHS_PARAMS = {
    "Ks_mult":                  {"lo": 4.0,   "hi": 10.5,   "scale": "linear"},
    "f_RS_abs":                 {"lo": 0.004, "hi": 0.030,  "scale": "log"},
    "channelconductivity_mmhr": {"lo": 30.0,  "hi": 1000.0, "scale": "log"},
}
PARAM_COLS = list(LHS_PARAMS.keys())
LHS_SERIES = "110"
STAGE2_CSV_NAME = "lhs_results_joint_Ks_f_cc_%s.csv" % LHS_SERIES
CONTROL_CSV_NAME = "lhs_results_joint_Ks_f_cc_CONTROL_%s.csv" % LHS_SERIES
STORM1_TRUTH_CSV = "truth_location_110.csv"
FEASIBLE_PBIAS_PCT = 2.0
MATCH_RTOL = 1e-3

# Stage 2's sampled box. A truth must lie inside it for a recovery test to be fair.
BOX = {
    "Ks_mult":  (4.0,   10.5,   "lin"),
    "f_RS_abs": (0.004, 0.030,  "log"),
    "cc_mmhr":  (30.0,  1000.0, "log"),
}
DEFAULT_KS = [5.0, 7.0, 9.0]
DEFAULT_F = [0.007, 0.012, 0.020]
DEFAULT_CC = [95.1, 201.0, 425.0]

SRC_SDF_DEFAULT = "../init_data/met/Master_Precip.sdf"
RAIN_DT_HR = float(getattr(builder, "RAIN_INTERVAL", 0.25))    # one MDF row = this many hours
NO_DATA_MIN = 9999.0                                           # tRIBS NO_DATA flag is 9999.99
PROV_NAME = "PROVENANCE_truth_location_110.json"

EDGE_FRAC = 0.10            # truth within this fraction of a box end on any axis -> NEAR_BOX_EDGE
NEAR_DRY_VOL_RATIO = 0.10   # truth volume / reference volume below this -> NEAR_DRY
MONO_TOL = 1e-9
OBS_PEAK_RTOL = 0.02        # scorer's observed peak must match the expected file's within this
GATE_MIN_CHANGE = 0.005     # the reference must differ from the original storm's by more than 0.5%
DRY_PROBLEM = "no flow at all"
EXPECTED_QOUT_COLS = 3

# Other scripts that share current_run_config.json (checked before starting). This
# script's own name matches "run_cc_", so a second copy is caught too.
COMPETING_PATTERNS = (
    "run_sensitivity_single", "run_sensitivity_sweep", "build_sensitivity_run",
    "run_cc_", "run_lhs_", "run_crossversion", "verify_truth_provenance",
    "validate_truth_point", "rescore_series100", "run_truth_",
)


class Abort(RuntimeError):
    """A problem that would make every further run pointless: stop the whole script."""


class BuildFailed(RuntimeError):
    """The builder raised for one run: log it and carry on."""


# ------------------------------------------------------------------
# Names
# ------------------------------------------------------------------
def storm_names(scale):
    """(label, short) for a rain scale: 0.8 -> ('storm080', 'st080')."""
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("--rain_scale must be a positive number")
    if abs(scale - 1.0) < 1e-9:
        raise ValueError("--rain_scale 1.0 is the original storm (nothing to scale). Try 0.8 or 1.25.")
    if scale < 0.1 or scale > 3.0:
        raise ValueError("--rain_scale %g is outside 0.1-3.0; that is not a sensible storm change." % scale)
    pct = int(round(scale * 100))
    if abs(pct / 100.0 - scale) > 1e-9:
        raise ValueError("--rain_scale needs at most 2 decimal places (0.8, 0.85, 1.25), got %r" % scale)
    return "storm%03d" % pct, "st%03d" % pct


def label(v):
    """The builder's own number-to-text rule (7.0 -> 7p0, 0.012 -> 0p012)."""
    return builder.value_to_label(v)


def label_roundtrips(v):
    try:
        return math.isclose(float(label(v).replace("p", ".")), float(v), rel_tol=1e-12, abs_tol=0.0)
    except ValueError:
        return False


def derived_name(path_str, suffix):
    """'../init_data/met/precip_SMF_1.mdf' + 'storm080' -> '../init_data/met/precip_SMF_1_storm080.mdf'."""
    head, sep, name = path_str.rpartition("/")
    stem, dot, ext = name.rpartition(".")
    if not dot:
        return "%s%s%s_%s" % (head, sep, name, suffix)
    return "%s%s%s_%s.%s" % (head, sep, stem, suffix, ext)


def predicted_run_id(cc, tag):
    """The run_id the builder will produce (computed first so a finished run can be
    skipped without rebuilding it)."""
    base, _ = builder.build_run_id(CC_PARAM, cc)
    return "%s_%s" % (base, tag)


# ------------------------------------------------------------------
# Grid of truths
# ------------------------------------------------------------------
def box_pos(v, key):
    """0 at the low end of the sampled range, 1 at the high end (Ks linear, f and cc log)."""
    lo, hi, scale = BOX[key]
    if scale == "log":
        return float((np.log10(v) - np.log10(lo)) / (np.log10(hi) - np.log10(lo)))
    return float((v - lo) / (hi - lo))


def validate_axis(name, values, key):
    vals = [float(v) for v in values]
    if not vals:
        raise ValueError("no %s values given" % name)
    if any((not np.isfinite(v)) or v <= 0 for v in vals):
        raise ValueError("%s values must be finite and > 0: %s" % (name, vals))
    lo, hi, _ = BOX[key]
    outside = [v for v in vals if not (lo <= v <= hi)]
    if outside:
        raise ValueError(
            "%s value(s) %s lie outside the sampled box %g-%g. A truth outside the box can never be "
            "recovered by a sweep that did not look there. Pick values inside the box."
            % (name, ", ".join("%g" % v for v in outside), lo, hi))
    bad = [v for v in vals if not label_roundtrips(v)]
    if bad:
        raise ValueError("%s value(s) %s need more than 6 decimal places, which the run labels cannot "
                         "spell exactly. Round them." % (name, ", ".join("%r" % v for v in bad)))
    uniq = sorted(set(vals))
    if len(uniq) != len(vals):
        raise ValueError("%s has repeated values: %s" % (name, vals))
    return uniq


def make_truths(ks_vals, f_vals, cc_vals, short):
    """Every (Ks, f, cc) combination, centre-out. Each item is a dict."""
    c_k = box_pos(STANDARD_KS, "Ks_mult")
    c_f = box_pos(STANDARD_F, "f_RS_abs")
    items = []
    for ks in ks_vals:
        for f in f_vals:
            for cc in cc_vals:
                pk, pf, pc = box_pos(ks, "Ks_mult"), box_pos(f, "f_RS_abs"), box_pos(cc, "cc_mmhr")
                tag = "K%s_f%s_%s" % (label(ks), label(f), short)
                items.append({
                    "Ks_mult": ks, "f_RS_abs": f, "cc_mmhr": cc,
                    "tag": tag, "run_id": predicted_run_id(cc, tag),
                    "pos": (pk, pf, pc),
                    "cell_dist": float(np.hypot(pk - c_k, pf - c_f)),
                })
    items.sort(key=lambda d: (round(d["cell_dist"], 9), d["Ks_mult"], d["f_RS_abs"], d["cc_mmhr"]))
    return items


def generate_lhs_samples(n, params, seed=None):
    """Log- or linear-stratified LHS. COPIED UNCHANGED from run_cc_stage2_joint_lhs_110.py so
    the draws are identical (same seed and n give the same points)."""
    rng = np.random.default_rng(seed)
    samples = {}
    for param, bounds in params.items():
        lo, hi = bounds["lo"], bounds["hi"]
        scale = bounds.get("scale", "linear")
        if scale == "log":
            log_lo, log_hi = np.log10(lo), np.log10(hi)
            intervals = np.linspace(log_lo, log_hi, n + 1)
            log_points = rng.uniform(intervals[:-1], intervals[1:])
            points = 10 ** log_points
        else:
            intervals = np.linspace(lo, hi, n + 1)
            points = rng.uniform(intervals[:-1], intervals[1:])
        rng.shuffle(points)
        samples[param] = points
    return pd.DataFrame(samples)


# ------------------------------------------------------------------
# Small helpers (ported unchanged from run_truth_location_110.py)
# ------------------------------------------------------------------
def md5_of(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def _num(x):
    try:
        v = float(x)
        return v if np.isfinite(v) else np.nan
    except (TypeError, ValueError):
        return np.nan


def _truthy(v):
    """True for True / 'True' / 1 (values read back from a CSV may be bool, text or number)."""
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
    """Timeout-safe subprocess execution -- same pattern as every other Series 110 script.
    Own process group so a hang is killed as a unit. Returns the scorer's output too."""
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
        except ProcessLookupError:
            pass
        stdout, _ = proc.communicate()
        timed_out = True
    except KeyboardInterrupt:
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        except ProcessLookupError:
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


def is_dry(info_or_row):
    """True when the ONLY thing wrong with a truth file is that it carries no flow."""
    return str(info_or_row.get("qout_problem", "")).strip() == DRY_PROBLEM


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


def cleanup_raw_results(raw_dir):
    """Delete a run's raw results directory (pixel files, spatial snapshots). Never raises."""
    try:
        if raw_dir.exists():
            shutil.rmtree(raw_dir)
            return True
    except Exception as e:
        print("  (cleanup warning: could not remove %s: %s)" % (raw_dir, e))
    return False


def write_json_atomic(path, obj):
    tmp = path.with_name(path.name + ".part")
    tmp.write_text(json.dumps(obj, indent=2, default=str))
    os.replace(tmp, path)


# ------------------------------------------------------------------
# 1. FORCING: scaled copies of the two rain files and the SDF
# ------------------------------------------------------------------
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


def find_rain_column(header_tokens):
    """Index of the rainfall column in a rain gauge file header. Accepts the tRIBS name R and
    common variants (R_mm, Rain, Precip...). If no column is named like rainfall but the header
    is the classic 5 columns starting with a year column, the last one is rainfall."""
    names = [h.strip("\"'").upper() for h in header_tokens]
    cands = [i for i, n in enumerate(names)
             if n == "R" or n.startswith("R_") or n.startswith("RAIN") or n.startswith("PRECIP")]
    if len(cands) == 1:
        return cands[0]
    if not cands and len(names) == 5 and names[0].startswith("Y"):
        return 4
    raise ValueError("cannot tell which column is rainfall (columns: %s)" % ", ".join(header_tokens))


def scale_mdf_text(src_path, factor):
    """Scale the R (rainfall rate) column of a rain gauge file by `factor`. Every other column
    is copied through as text. NO_DATA flags (>= 9999) and negative values are left alone.
    Returns (new_text, stats)."""
    raw = Path(src_path).read_text().splitlines()
    lines = [ln for ln in raw if ln.strip()]
    if not lines:
        raise ValueError("%s is empty" % src_path)
    header, _ = split_keep(lines[0])
    try:
        ri = find_rain_column(header)
    except ValueError as e:
        raise ValueError("unexpected header in %s: %r (%s)" % (src_path, lines[0], e))
    out = [lines[0]]
    st = {"rows": 0, "scaled": 0, "kept_flag": 0, "sum_before": 0.0, "sum_after": 0.0,
          "peak_before": 0.0, "peak_after": 0.0, "header": lines[0].strip(), "rain_col": header[ri]}
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
        if (not math.isfinite(r)) or r >= NO_DATA_MIN or r < 0:
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
    print("  rain gauge files (only the rainfall column is multiplied; the date and time columns are copied unchanged):")
    for s in plan["stations"]:
        st = s["stats"]
        print("    %-34s %5d rows   column '%s'   peak %7.3f -> %7.3f   sum x %g h: %8.2f -> %8.2f   (x%g)"
              % (s["dst"].name, st["rows"], st["rain_col"], st["peak_before"], st["peak_after"], RAIN_DT_HR,
                 st["sum_before"] * RAIN_DT_HR, st["sum_after"] * RAIN_DT_HR, plan["scale"]))
        print("      header: %s" % st["header"])
        if st["kept_flag"]:
            print("      (%d NO_DATA/negative row(s) left unscaled)" % st["kept_flag"])
    print("    (values are in the file's own units; the 'sum x %g h' is only a total if they are mm/hr on %g-hour rows)"
          % (RAIN_DT_HR, RAIN_DT_HR))
    print("    SDF: %s  ->  %s" % (plan["src_sdf"].name, plan["dst_sdf"].name))


def write_forcing(plan, overwrite):
    """Write the planned files (atomically). A file that already exists with identical content
    is left alone; one that differs is an error unless overwrite. Returns {name: status}."""
    todo = [(s["dst"], s["text"]) for s in plan["stations"]] + [(plan["dst_sdf"], plan["sdf_text"])]
    status = {}
    for dst, text in todo:
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
        ri = find_rain_column(a[0])
        n_bad = 0
        for ra, rb in zip(a[1:], b[1:]):
            if len(ra) != len(rb) or any(x != y for k, (x, y) in enumerate(zip(ra, rb)) if k != ri):
                n_bad += 1
                continue
            r0, r1 = float(ra[ri]), float(rb[ri])
            want = r0 if ((not math.isfinite(r0)) or r0 >= NO_DATA_MIN or r0 < 0) else r0 * factor
            if (math.isfinite(want) or math.isfinite(r1)) and abs(r1 - want) > 1e-6:
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


# ------------------------------------------------------------------
# 2. RUN MACHINERY: build -> patch config -> score -> verify observed
# ------------------------------------------------------------------
def check_config(ctx, cfg, run_id, expect):
    """What the builder wrote into current_run_config.json versus what was asked for."""
    problems = []
    if cfg.get("run_id") != run_id:
        problems.append("run_id is %r, expected %r" % (cfg.get("run_id"), run_id))
    if "gauge_sdf" in cfg and cfg["gauge_sdf"] != ctx.gauge_sdf_rel:
        problems.append("gauge_sdf is %r, expected %r" % (cfg["gauge_sdf"], ctx.gauge_sdf_rel))
    wanted = [("Ks_mult", expect["Ks_mult"]), ("f_RS_abs", expect["f_RS_abs"]),
              (CC_PARAM, expect["cc"]), ("optpercolation", expect["optpercolation"]),
              ("kinemvelcoef", ROUTING_TRUTH["kinemvelcoef"]), ("flowexp", ROUTING_TRUTH["flowexp"]),
              ("channelroughness", ROUTING_TRUTH["channelroughness"])]
    for key, want in wanted:
        if key not in cfg:
            problems.append("%s missing from the run configuration" % key)
        elif not _close(cfg[key], want):
            problems.append("%s is %r, expected %r" % (key, cfg[key], want))
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


def prepare_run(ctx, pred_id, cc_val, overrides, tag, truth_file_rel, expect):
    """Build the input file with the scaled rain, check the configuration the builder wrote, and
    (for runs that must be scored against the new storm's reference) add truth_file to it."""
    try:
        run_id, _, _ = builder.build_input_file(CC_PARAM, cc_val, overrides=overrides, tag=tag,
                                                gauge_sdf=ctx.gauge_sdf_rel)
    except Exception as e:
        raise BuildFailed(str(e)[:300])
    if run_id != pred_id:
        raise Abort("builder produced run_id %s, expected %s" % (run_id, pred_id))
    cfg_path = ctx.calib_dir / "current_run_config.json"
    try:
        cfg = json.loads(cfg_path.read_text())
    except Exception as e:
        raise Abort("cannot read %s after the build: %s" % (cfg_path, e))
    problems = check_config(ctx, cfg, run_id, expect)
    if problems:
        raise Abort("the builder's run configuration does not match what was asked for: " + "; ".join(problems))
    if truth_file_rel is not None:
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
        problems.append("obs_mode is %r, expected 'synth'" % metrics.get("obs_mode"))
    m = re.search(r"Reading observed from:\s*(.+)", stdout or "")
    if m and m.group(1).strip() != expected_obs_path.name:
        problems.append("the scorer read %r as observed, expected %r" % (m.group(1).strip(), expected_obs_path.name))
    try:
        exp = expected_obs_peak(expected_obs_path, cfg["event_start"], cfg["event_end"])
        got = float(metrics["obs_peak_m3s"])
        if not (abs(got - exp) <= OBS_PEAK_RTOL * max(abs(exp), 1e-12)):
            problems.append("the observed peak the scorer used (%.4g m3/s) is not the peak of %s (%.4g m3/s)"
                            % (got, expected_obs_path.name, exp))
    except Exception as e:
        print("  (note: could not double-check the observed peak: %s)" % str(e)[:100])
    return problems


def run_one(ctx, pred_id, cc_val, overrides, tag, truth_file_rel, expect, expected_obs_path):
    """Build, score, and verify one run. Returns a dict with status in SUCCESS / BUILD_FAILED /
    HANG / FAILED. Raises Abort for a problem that would affect every run."""
    out = {"run_id": pred_id, "status": None, "reason": "", "elapsed": 0.0, "metrics": None}
    try:
        run_id, cfg = prepare_run(ctx, pred_id, cc_val, overrides, tag, truth_file_rel, expect)
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
    metrics = pd.read_csv(metrics_file).iloc[0].to_dict()
    problems = check_observed(metrics, cfg, stdout, expected_obs_path)
    if problems:
        raise Abort("scoring did not use the intended 'observed' file for %s: %s. Every later run would be "
                    "scored the same wrong way, so the script is stopping." % (run_id, "; ".join(problems)))
    out.update(status="SUCCESS", metrics=metrics)
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
    (status, reason, md5, qinfo, raw_dir)."""
    raw_dir = ctx.calib_dir / "02_results" / ctx.category / run_id
    raw_qout = locate_raw_qout(raw_dir, run_id)
    if raw_qout is None:
        return "NO_QOUT", "no single *.qout found in %s" % raw_dir, None, None, raw_dir
    try:
        md5 = copy_verified(raw_qout, dest_path)
    except Exception as e:
        return "COPY_FAILED", str(e)[:300], None, None, raw_dir
    return "OK", "", md5, describe_qout(dest_path), raw_dir


# ------------------------------------------------------------------
# Rows, flags, provenance
# ------------------------------------------------------------------
def make_flags(pos, vol_ratio, qout_ok):
    flags = []
    if any((p < EDGE_FRAC or p > 1.0 - EDGE_FRAC) for p in pos):
        flags.append("NEAR_BOX_EDGE")
    if np.isfinite(vol_ratio) and vol_ratio < NEAR_DRY_VOL_RATIO:
        flags.append("NEAR_DRY")
    if not qout_ok:
        flags.append("BAD_QOUT")
    return ";".join(flags)


def assemble_truth_row(ctx, t, run_id, metrics, cand_path, md5, qinfo):
    obs_vol = _num(metrics.get("obs_volume_m3"))
    sim_vol = _num(metrics.get("sim_volume_m3"))
    obs_peak = _num(metrics.get("obs_peak_m3s"))
    sim_peak = _num(metrics.get("sim_peak_m3s"))
    vol_ratio = sim_vol / obs_vol if obs_vol and np.isfinite(obs_vol) else np.nan
    peak_ratio = sim_peak / obs_peak if obs_peak and np.isfinite(obs_peak) else np.nan
    pk, pf, pc = t["pos"]
    row = {
        "run_id":                run_id,
        "storm_label":           ctx.label,
        "rain_scale":            ctx.scale,
        "Ks_mult":               t["Ks_mult"],
        "f_RS_abs":              t["f_RS_abs"],
        "cc_mmhr":               t["cc_mmhr"],
        "box_pos_Ks":            pk,
        "box_pos_f":             pf,
        "box_pos_cc":            pc,
        "flags":                 make_flags(t["pos"], vol_ratio, qinfo["qout_ok"] or is_dry(qinfo)),
        "reference_truth_ccoff": ctx.ref_path.name,
        "pbias_pct":             _num(metrics.get("pbias_pct")),
        "kge_2012":              _num(metrics.get("kge_2012")),
        "vol_ratio_vs_ccoff":    vol_ratio,
        "peak_ratio_vs_ccoff":   peak_ratio,
        "candidate_qout":        cand_path.name,
        "candidate_md5":         md5,
        "candidate_bytes":       int(cand_path.stat().st_size),
    }
    row.update(qinfo)
    for k, v in metrics.items():
        row.setdefault(k, v)
    row["optpercolation"] = 1
    return row


def assemble_ref_row(ctx, run_id, metrics, cand_path, md5, qinfo):
    obs_vol = _num(metrics.get("obs_volume_m3"))
    sim_vol = _num(metrics.get("sim_volume_m3"))
    obs_peak = _num(metrics.get("obs_peak_m3s"))
    sim_peak = _num(metrics.get("sim_peak_m3s"))
    row = {
        "run_id":                 run_id,
        "storm_label":            ctx.label,
        "rain_scale":             ctx.scale,
        "Ks_mult":                STANDARD_KS,
        "f_RS_abs":               STANDARD_F,
        "cc_mmhr":                REF_CC,
        "optpercolation":         0,
        "reference_storm1_ccoff": ctx.storm1_ref.name,
        "vol_ratio_vs_storm1":    sim_vol / obs_vol if obs_vol and np.isfinite(obs_vol) else np.nan,
        "peak_ratio_vs_storm1":   sim_peak / obs_peak if obs_peak and np.isfinite(obs_peak) else np.nan,
        "pbias_pct":              _num(metrics.get("pbias_pct")),
        "kge_2012":               _num(metrics.get("kge_2012")),
        "candidate_qout":         cand_path.name,
        "candidate_md5":          md5,
        "candidate_bytes":        int(cand_path.stat().st_size),
    }
    row.update(qinfo)
    for k, v in metrics.items():
        row.setdefault(k, v)
    row["optpercolation"] = 0
    return row


def write_truth_provenance(ctx, rows):
    prov = dict(ctx.truth_prov_base)
    prov["updated_local"] = datetime.now().isoformat(timespec="seconds")
    prov["files"] = {
        r["candidate_qout"]: {
            "Ks_mult":  float(r["Ks_mult"]),
            "f_RS_abs": float(r["f_RS_abs"]),
            "cc_mmhr":  float(r["cc_mmhr"]),
            "md5":      r["candidate_md5"],
            "bytes":    int(r["candidate_bytes"]),
            "run_id":   r["run_id"],
            "flags":    "" if pd.isna(r.get("flags", "")) else str(r.get("flags", "")),
        }
        for r in rows if r.get("candidate_qout")
    }
    write_json_atomic(ctx.loc_dir / PROV_NAME, prov)


# ------------------------------------------------------------------
# STAGE: forcing
# ------------------------------------------------------------------
def stage_forcing(ctx):
    print("\n" + "=" * 92)
    print("STAGE 1 of 4 -- FORCING: scaled rain files for %s (x%g)" % (ctx.label, ctx.scale))
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
                        "Y/M/D/H unchanged; NO_DATA flags unchanged; whole record scaled.",
        "row_hours_assumed_for_totals": RAIN_DT_HR,
        "source_sdf": {"file": plan["src_sdf"].name, "md5": md5_of(plan["src_sdf"])},
        "scaled_sdf": {"file": plan["dst_sdf"].name, "md5": md5_of(plan["dst_sdf"])},
        "stations": [{
            "source": s["src"].name, "scaled": s["dst"].name, "scaled_md5": md5_of(s["dst"]),
            "rows": s["stats"]["rows"], "peak_mm_hr": [s["stats"]["peak_before"], s["stats"]["peak_after"]],
            "total_mm": [s["stats"]["sum_before"] * RAIN_DT_HR, s["stats"]["sum_after"] * RAIN_DT_HR],
            "no_data_rows_unscaled": s["stats"]["kept_flag"],
        } for s in plan["stations"]],
        "python": sys.version.split()[0], "pandas": pd.__version__, "numpy": np.__version__,
    }
    write_json_atomic(ctx.calib_dir / ("forcing_%s_PROVENANCE.json" % ctx.label), prov)
    print("  Saved: %s" % ("forcing_%s_PROVENANCE.json" % ctx.label))


# ------------------------------------------------------------------
# STAGE: reference (cc OFF) for the new storm, with the direction gate
# ------------------------------------------------------------------
def load_rows(path):
    if path.exists():
        try:
            return pd.read_csv(path).to_dict("records")
        except Exception as e:
            print("  Warning: could not load %s (%s). Starting fresh." % (path.name, e))
    return []


def gate_ref(ctx, row):
    """The new storm's cc-OFF reference must differ from the original storm's in the right direction."""
    vr, pr = _num(row.get("vol_ratio_vs_storm1")), _num(row.get("peak_ratio_vs_storm1"))
    print("  REFERENCE GATE: this storm's cc-OFF reference is %.1f%% of the original storm's peak and "
          "%.1f%% of its volume (rain x%g)." % (100 * pr, 100 * vr, ctx.scale))
    if not _truthy(row.get("qout_ok", False)):
        raise Abort("the reference failed its sanity check (%s). Nothing sensible can be scored against it."
                    % row.get("qout_problem"))
    if not (np.isfinite(vr) and np.isfinite(pr)):
        raise Abort("could not compute the reference's peak/volume relative to the original storm.")
    direction = np.sign(ctx.scale - 1.0)
    if (np.sign(pr - 1.0) != direction or np.sign(vr - 1.0) != direction
            or abs(pr - 1.0) < GATE_MIN_CHANGE or abs(vr - 1.0) < GATE_MIN_CHANGE):
        raise Abort(
            "GATE FAILED. With rain x%g the peak and volume should both move %s from the original storm's, "
            "but they are %.1f%% and %.1f%% of it. The scaled rainfall is apparently NOT being used by tRIBS "
            "(check that %s is what the .in file names, and that tRIBS reads it). Nothing else was run."
            % (ctx.scale, "DOWN" if direction < 0 else "UP", 100 * pr, 100 * vr, ctx.gauge_sdf_rel))
    print("  REFERENCE GATE passed: the change is in the right direction.")


def stage_reference(ctx, allow_run):
    print("\n" + "=" * 92)
    print("STAGE 2 of 4 -- REFERENCE: cc-OFF truth for %s (Ks %g, f %g, routing at truth)"
          % (ctx.label, STANDARD_KS, STANDARD_F))
    print("=" * 92)
    tag = "ref_%s" % ctx.short
    pred_id = predicted_run_id(REF_CC, tag)
    ref_path = ctx.ref_dir / (pred_id + "_Outlet.qout")
    ctx.ref_path = ref_path
    rows = load_rows(ctx.ref_csv)
    row = next((r for r in rows if r.get("run_id") == pred_id), None)

    reusable = False
    if row is not None and ref_path.exists():
        rec = str(row.get("candidate_md5", ""))
        reusable = bool(rec) and md5_of(ref_path) == rec and _truthy(row.get("qout_ok", False))

    if reusable and (ctx.args.skip_existing or not allow_run):
        print("  Reference already complete (file present, checksum matches): %s" % ref_path.name)
    else:
        if not allow_run:
            raise Abort("the %s reference does not exist yet (or fails its checksum). Run this script "
                        "WITHOUT --only sweep / --control first." % ctx.label)
        ctx.ref_dir.mkdir(parents=True, exist_ok=True)
        print("  -> %s" % pred_id, flush=True)
        overrides = dict(ROUTING_TRUTH)
        overrides.update({"Ks_mult": STANDARD_KS, "f_RS_abs": STANDARD_F, "optpercolation": 0})
        expect = {"Ks_mult": STANDARD_KS, "f_RS_abs": STANDARD_F, "cc": REF_CC, "optpercolation": 0}
        # Scored against the ORIGINAL storm's reference (the file in synth_truth/): that is exactly
        # what gives the cross-storm gate below.
        res = run_one(ctx, pred_id, REF_CC, overrides, tag, None, expect, ctx.storm1_ref)
        if res["status"] != "SUCCESS":
            raise Abort("the reference run %s: %s. Cannot continue without it." % (res["status"], res["reason"]))
        st, reason, md5, qinfo, raw_dir = harvest_qout(ctx, pred_id, ref_path)
        if st != "OK":
            raise Abort("the reference output could not be saved (%s: %s)." % (st, reason))
        row = assemble_ref_row(ctx, pred_id, res["metrics"], ref_path, md5, qinfo)
        rows = [r for r in rows if r.get("run_id") != pred_id] + [row]
        pd.DataFrame(rows).to_csv(ctx.ref_csv, index=False)
        print("  saved reference: %s  (%s bytes, md5 %s, %d rows)"
              % (ref_path.name, format(row["candidate_bytes"], ","), md5[:12], qinfo["qout_rows"]))
        if not ctx.args.no_cleanup and (qinfo["qout_ok"] or is_dry(qinfo)):
            cleanup_raw_results(raw_dir)
    gate_ref(ctx, row)
    ctx.ref_row = row
    ctx.ref_md5 = str(row.get("candidate_md5"))


# ------------------------------------------------------------------
# STAGE: the 27 cc-ON truths under the new storm
# ------------------------------------------------------------------
def report_truths(ctx, df):
    print("\n" + "=" * 100)
    print("cc-ON TRUTHS UNDER %s (rain x%g).  Metrics are vs this storm's cc-OFF reference: %s"
          % (ctx.label, ctx.scale, ctx.ref_path.name))
    print("=" * 100)
    d = df.sort_values(["Ks_mult", "f_RS_abs", "cc_mmhr"]).reset_index(drop=True)
    d["flags"] = d["flags"].fillna("") if "flags" in d else ""

    s1 = None
    p1 = ctx.summary_dir / STORM1_TRUTH_CSV
    if p1.exists():
        try:
            s1 = pd.read_csv(p1)
            s1["_k"] = list(zip(s1["Ks_mult"].round(6), s1["f_RS_abs"].round(6), s1["cc_mmhr"].round(6)))
        except Exception:
            s1 = None

    if s1 is not None:
        print("%6s %8s %8s | %9s %9s | vol kept: original -> %-8s | peak kept: original -> %-8s  flags"
              % ("Ks", "f", "cc", "PBIAS %", "KGE_2012", ctx.label, ctx.label))
    else:
        print("%6s %8s %8s | %9s %9s | %10s %10s  flags" % ("Ks", "f", "cc", "PBIAS %", "KGE_2012", "vol kept", "peak kept"))
    last = None
    for _, r in d.iterrows():
        cell = (r["Ks_mult"], r["f_RS_abs"])
        if last is not None and cell != last:
            print("")
        last = cell
        vk, pk = _num(r.get("vol_ratio_vs_ccoff")), _num(r.get("peak_ratio_vs_ccoff"))
        base = "%6g %8g %8g | %+9.2f %9.3f |" % (r["Ks_mult"], r["f_RS_abs"], r["cc_mmhr"],
                                                 _num(r.get("pbias_pct")), _num(r.get("kge_2012")))
        if s1 is not None:
            key = (round(float(r["Ks_mult"]), 6), round(float(r["f_RS_abs"]), 6), round(float(r["cc_mmhr"]), 6))
            m = s1[s1["_k"] == key]
            ov = _num(m.iloc[0].get("vol_ratio_vs_ccoff")) if len(m) else np.nan
            op = _num(m.iloc[0].get("peak_ratio_vs_ccoff")) if len(m) else np.nan
            print("%s %7.1f%% -> %6.1f%%        | %7.1f%% -> %6.1f%%         %s"
                  % (base, 100 * ov, 100 * vk, 100 * op, 100 * pk, r["flags"]))
        else:
            print("%s %9.1f%% %9.1f%%  %s" % (base, 100 * vk, 100 * pk, r["flags"]))
    print("  (kept = the cc-ON truth's volume / peak as a percentage of the cc-OFF reference for the SAME storm, same "
          "event window.\n   They move with Ks and f as well as cc. A lower number under the new storm than under the "
          "original one means channel loss removes a larger share of the flow in that storm.)")

    if "sim_volume_m3" in d:
        rev = []
        for (ks, f), g in d.groupby(["Ks_mult", "f_RS_abs"]):
            g = g.sort_values("cc_mmhr")
            v, cc = g["sim_volume_m3"].astype(float).to_numpy(), g["cc_mmhr"].to_numpy(float)
            for i in range(len(v) - 1):
                if np.isfinite(v[i]) and np.isfinite(v[i + 1]) and v[i + 1] > v[i] * (1 + MONO_TOL):
                    rev.append((ks, f, cc[i], cc[i + 1]))
        if rev:
            print("\n  CHECK cc: outlet volume ROSE as cc rose at %d step(s) (it should fall):" % len(rev))
            for ks, f, a, b in rev[:6]:
                print("    Ks %g, f %g: cc %g -> %g" % (ks, f, a, b))
        else:
            print("\n  CHECK cc: at every Ks and f, outlet volume falls (or holds) as cc rises.")
        rev = []
        for (f, cc), g in d.groupby(["f_RS_abs", "cc_mmhr"]):
            g = g.sort_values("Ks_mult")
            v, ks = g["sim_volume_m3"].astype(float).to_numpy(), g["Ks_mult"].to_numpy(float)
            for i in range(len(v) - 1):
                if np.isfinite(v[i]) and np.isfinite(v[i + 1]) and v[i + 1] > v[i] * (1 + MONO_TOL):
                    rev.append((f, cc, ks[i], ks[i + 1]))
        if rev:
            print("  CHECK Ks: outlet volume ROSE as Ks rose at %d step(s) (it should fall):" % len(rev))
            for f, cc, a, b in rev[:6]:
                print("    f %g, cc %g: Ks %g -> %g" % (f, cc, a, b))
        else:
            print("  CHECK Ks: at every f and cc, outlet volume falls (or holds) as Ks rises.")

    n_dry = int(sum(1 for _, r in d.iterrows() if is_dry(r)))
    n_bad = int(sum(1 for _, r in d.iterrows() if not _truthy(r.get("qout_ok")) and not is_dry(r)))
    if n_dry:
        print("\n  %d truth(s) carry NO outlet flow at all under this storm (the channel takes everything). That is "
              "the physics, not a failure, but nothing can be scored against a flat-zero truth, so the re-scoring "
              "step will have to leave it out." % n_dry)
    if n_bad > 0:
        print("\n  WARNING: %d truth file(s) failed the sanity check (see qout_problem in the CSV)." % n_bad)
    if "obs_mode" in d and (d["obs_mode"] != "synth").any():
        print("  WARNING: some rows were NOT scored in synthetic-truth mode; their reference metrics are invalid.")
    flagged = d[d["flags"] != ""]
    if len(flagged):
        print("\n  %d truth(s) carry a flag: %s" % (len(flagged), "; ".join(
            "Ks %g f %g cc %g [%s]" % (r.Ks_mult, r.f_RS_abs, r.cc_mmhr, r.flags) for r in flagged.itertuples())))


def stage_truths(ctx):
    truths = ctx.truths
    print("\n" + "=" * 92)
    print("STAGE 3 of 4 -- TRUTHS: %d cc-ON truth(s) under %s" % (len(truths), ctx.label))
    print("=" * 92)
    ctx.loc_dir.mkdir(parents=True, exist_ok=True)
    rows = load_rows(ctx.truth_csv)
    if rows:
        print("  Loaded existing results: %d rows from %s" % (len(rows), ctx.truth_csv.name))
    done = {r.get("run_id"): r for r in rows}
    failed_rows = []
    completed = skipped = hung = failed = cleaned = 0
    t_start = time.time()
    truth_rel = "%s/%s" % (ctx.ref_dir.name, ctx.ref_path.name)

    for i, t in enumerate(truths, start=1):
        pred_id = t["run_id"]
        cand_path = ctx.loc_dir / ("%s_Outlet.qout" % pred_id)
        print("\n[%3d/%d]  Ks=%g  f=%g  cc=%g mm/hr  -> %s"
              % (i, len(truths), t["Ks_mult"], t["f_RS_abs"], t["cc_mmhr"], pred_id), flush=True)

        if ctx.args.skip_existing and pred_id in done and cand_path.exists():
            rec_md5 = str(done[pred_id].get("candidate_md5", ""))
            usable = _truthy(done[pred_id].get("qout_ok", False)) or is_dry(done[pred_id])
            if not usable:
                print("  earlier file failed its sanity check; re-running")
            elif rec_md5 and md5_of(cand_path) == rec_md5:
                print("  SKIP (already completed; file present and checksum matches)")
                skipped += 1
                continue
            else:
                print("  file present but its checksum does not match the record; re-running")

        overrides = dict(ROUTING_TRUTH)
        overrides.update({"Ks_mult": t["Ks_mult"], "f_RS_abs": t["f_RS_abs"], "optpercolation": 1})
        expect = {"Ks_mult": t["Ks_mult"], "f_RS_abs": t["f_RS_abs"], "cc": t["cc_mmhr"], "optpercolation": 1}
        res = run_one(ctx, pred_id, t["cc_mmhr"], overrides, t["tag"], truth_rel, expect, ctx.ref_path)
        status, reason = res["status"], res["reason"]

        if status == "SUCCESS":
            st, why, md5, qinfo, raw_dir = harvest_qout(ctx, pred_id, cand_path)
            if st != "OK":
                status, reason = st, why

        if status != "SUCCESS":
            print("  %s: %s  (%s)" % (status, pred_id, reason))
            if status == "HANG":
                hung += 1
            else:
                failed += 1
            failed_rows.append({"run_id": pred_id, "status": status, "reason": reason,
                                "elapsed_min": res["elapsed"] / 60, "Ks_mult": t["Ks_mult"],
                                "f_RS_abs": t["f_RS_abs"], "cc_mmhr": t["cc_mmhr"]})
            pd.DataFrame(failed_rows).to_csv(ctx.truth_failed_csv, index=False)
            note_outcome(ctx, False)
            continue
        note_outcome(ctx, True)

        row = assemble_truth_row(ctx, t, pred_id, res["metrics"], cand_path, md5, qinfo)
        rows = [r for r in rows if r.get("run_id") != pred_id] + [row]
        done[pred_id] = row
        completed += 1
        print("  vs this storm's cc-OFF reference: PBIAS=%+.2f%%  KGE_2012=%.3f  volume kept=%.1f%%  peak kept=%.1f%%%s"
              % (row["pbias_pct"], row["kge_2012"], 100 * row["vol_ratio_vs_ccoff"],
                 100 * row["peak_ratio_vs_ccoff"], ("  [%s]" % row["flags"]) if row["flags"] else ""))
        print("  saved truth: %s  (%s bytes, md5 %s, %d rows)"
              % (cand_path.name, format(row["candidate_bytes"], ","), md5[:12], qinfo["qout_rows"]))
        if not qinfo["qout_ok"]:
            print("  %s: %s" % ("NOTE" if is_dry(qinfo) else "WARNING: truth failed sanity check",
                                "no outlet flow at all under this storm" if is_dry(qinfo) else qinfo["qout_problem"]))
        if not ctx.args.no_cleanup and (qinfo["qout_ok"] or is_dry(qinfo)):
            if cleanup_raw_results(raw_dir):
                cleaned += 1

        remaining = len(truths) - i
        if completed > 0:
            eta = (time.time() - t_start) / completed * remaining / 60
            print("  Run time: %.1f min  |  ETA: %.0f min remaining" % (res["elapsed"] / 60, eta), flush=True)
        pd.DataFrame(rows).sort_values(["Ks_mult", "f_RS_abs", "cc_mmhr"]).to_csv(ctx.truth_csv, index=False)
        write_truth_provenance(ctx, rows)

    print("\nTruths done: %d ran, %d skipped, %d hung, %d failed, %d raw dirs cleaned"
          % (completed, skipped, hung, failed, cleaned))
    if rows:
        final_df = pd.DataFrame(rows).sort_values(["Ks_mult", "f_RS_abs", "cc_mmhr"]).reset_index(drop=True)
        final_df.to_csv(ctx.truth_csv, index=False)
        write_truth_provenance(ctx, rows)
        print("Saved: %s  (%d rows)" % (ctx.truth_csv, len(final_df)))
        try:
            report_truths(ctx, final_df)
        except Exception as e:
            print("\n(report skipped: %s -- results are saved regardless)" % e)
    missing = [t for t in truths if t["run_id"] not in done
               or not (_truthy(done[t["run_id"]].get("qout_ok", False)) or is_dry(done[t["run_id"]]))]
    if missing:
        print("\n%d of the %d truth(s) are NOT complete (failed, hung, or a file that failed its sanity check):"
              % (len(missing), len(truths)))
        for t in missing[:10]:
            print("    %s" % t["run_id"])
        print("See %s for the reasons. Re-run the same command with --skip_existing to retry just those."
              % ctx.truth_failed_csv.name)
    else:
        print("\nAll %d truth(s) in this grid are complete." % len(truths))


# ------------------------------------------------------------------
# STAGE: the 250-run sweep
# ------------------------------------------------------------------
def load_prior(path):
    """Read a storm-1 results CSV (never written). None if missing or lacking the parameter columns."""
    if not path.exists():
        return None
    try:
        df = pd.read_csv(path)
    except Exception as e:
        print("  Warning: could not read %s (%s)." % (path.name, e))
        return None
    if any(c not in df.columns for c in PARAM_COLS):
        return None
    return df


def match_rows(points_df, other_df, rtol=MATCH_RTOL):
    """For each row of points_df, the positional index of the row of other_df with the same
    (Ks_mult, f_RS_abs, cc) within rtol, or -1."""
    out = np.full(len(points_df), -1, dtype=int)
    if other_df is None or any(c not in points_df.columns for c in PARAM_COLS):
        return out
    s2 = other_df[PARAM_COLS].to_numpy(dtype=float)
    pts = points_df[PARAM_COLS].to_numpy(dtype=float)
    for i, row in enumerate(pts):
        ok = np.all(np.isclose(s2, row, rtol=rtol, atol=0.0), axis=1)
        hits = np.flatnonzero(ok)
        if hits.size:
            out[i] = hits[0]
    return out


def report_sweep(ctx, df, control, on_path, off_path):
    print("\n" + "=" * 70)
    print("SWEEP DIAGNOSTICS (%s, %s)" % ("channel loss OFF" if control else "channel loss ON", ctx.label))
    print("=" * 70)
    if "pbias_pct" not in df.columns:
        print("  No pbias_pct column; skipping.")
        return
    pb = pd.to_numeric(df["pbias_pct"], errors="coerce").to_numpy(float)
    good = np.isfinite(pb)
    n_good = int(good.sum())
    print("  Scored runs: %d of %d rows" % (n_good, len(df)))
    if n_good < 3:
        return
    kge = pd.to_numeric(df.get("kge_2012"), errors="coerce").to_numpy(float) if "kge_2012" in df else None
    if kge is not None and np.isfinite(kge).any():
        b = int(np.nanargmax(kge))
        print("  Best run by KGE_2012: Ks %.3f  f %.4f  cc %.1f   KGE_2012 %.3f  PBIAS %+.1f%%   (the reference "
              "was built at Ks %g, f %g with channel loss OFF)"
              % (df.iloc[b]["Ks_mult"], df.iloc[b]["f_RS_abs"], df.iloc[b][CC_PARAM], kge[b], pb[b],
                 STANDARD_KS, STANDARD_F))
    n_feas = int((np.abs(pb[good]) < FEASIBLE_PBIAS_PCT).sum())
    print("  Feasible (|PBIAS| < %g%% vs the cc-OFF reference of this storm): %d/%d = %.1f%%"
          % (FEASIBLE_PBIAS_PCT, n_feas, n_good, 100.0 * n_feas / n_good))
    if CC_PARAM in df.columns:
        logcc = np.log10(pd.to_numeric(df[CC_PARAM], errors="coerce").to_numpy(float))
        both = good & np.isfinite(logcc)
        r = np.corrcoef(pb[both], logcc[both])[0, 1]
        if control:
            print("  corr(PBIAS, log10 cc) = %+.3f   -> %s" % (
                r, "OK, cc looks inert" if abs(r) < 0.1 else "NOT NEAR ZERO: cc may not be inert at optpercolation=0"))
        else:
            print("  corr(PBIAS, log10 cc) = %+.3f   (negative = more cc, less outflow)" % r)
    on_df = pd.read_csv(on_path) if on_path.exists() else None
    off_df = pd.read_csv(off_path) if off_path.exists() else None
    if on_df is not None and off_df is not None and "pbias_pct" in on_df and "pbias_pct" in off_df:
        tw = match_rows(on_df, off_df)
        ok = tw >= 0
        if ok.any():
            pon = pd.to_numeric(on_df["pbias_pct"], errors="coerce").to_numpy(float)[ok]
            poff = pd.to_numeric(off_df["pbias_pct"], errors="coerce").to_numpy(float)[tw[ok]]
            fin = np.isfinite(pon) & np.isfinite(poff)
            if fin.any():
                sh = pon[fin] - poff[fin]
                f_on = int((np.abs(pon[fin]) < FEASIBLE_PBIAS_PCT).sum())
                f_off = int((np.abs(poff[fin]) < FEASIBLE_PBIAS_PCT).sum())
                print("  Paired cc ON vs OFF on the same %d points (this storm): mean PBIAS shift from turning cc ON "
                      "%+.2f, median %+.2f; feasible ON %d, OFF %d" % (int(fin.sum()), sh.mean(), np.median(sh), f_on, f_off))
                print("    (original storm, for comparison: mean -13.23, median -10.71; feasible ON 23, OFF 32)")


def stage_sweep(ctx, control):
    args = ctx.args
    optperc = 0 if control else 1
    tag = ("ctl_%s" if control else "sw_%s") % ctx.short
    word = "CONTROL (channel loss OFF)" if control else "cc ON"
    mid = "CONTROL_" if control else ""
    out_path = ctx.summary_dir / ("lhs_results_joint_Ks_f_cc_%s%s_%s.csv" % (mid, ctx.label, LHS_SERIES))
    failed_path = ctx.summary_dir / ("lhs_results_joint_Ks_f_cc_%sFAILED_%s_%s.csv" % (mid, ctx.label, LHS_SERIES))
    on_path = ctx.summary_dir / ("lhs_results_joint_Ks_f_cc_%s_%s.csv" % (ctx.label, LHS_SERIES))
    off_path = ctx.summary_dir / ("lhs_results_joint_Ks_f_cc_CONTROL_%s_%s.csv" % (ctx.label, LHS_SERIES))
    prior_path = ctx.summary_dir / (CONTROL_CSV_NAME if control else STAGE2_CSV_NAME)

    print("\n" + "=" * 92)
    print("STAGE 4 of 4 -- SWEEP: %d-run joint Ks/f/cc LHS, %s, under %s (rain x%g)"
          % (args.n, word, ctx.label, ctx.scale))
    print("=" * 92)
    samples = generate_lhs_samples(args.n, LHS_PARAMS, seed=args.seed)

    prior = load_prior(prior_path)
    twin = match_rows(samples, prior)
    if prior is None:
        print("  PAIRING CHECK: %s not found/readable -- running UNPAIRED with the original storm." % prior_path.name)
    else:
        n_matched = int((twin >= 0).sum())
        print("  PAIRING CHECK: %d/%d samples match a point of the original-storm %s (%d rows)."
              % (n_matched, args.n, prior_path.name, len(prior)))
        if n_matched < 0.9 * min(args.n, len(prior)) and not args.allow_unpaired:
            raise Abort("poor pairing with the original storm's sweep: this run is using a different --n or --seed "
                        "than Stage 2 (or the LHS code changed). Stopping before spending hours on points that do "
                        "not pair. (--allow_unpaired overrides.)")

    rows = load_rows(out_path)
    if rows:
        print("  Loaded existing results: %d rows from %s" % (len(rows), out_path.name))
    done_ids = {r.get("run_id") for r in rows}
    failed_rows = []
    truth_rel = "%s/%s" % (ctx.ref_dir.name, ctx.ref_path.name)
    completed = skipped = hung = failed = cleaned = 0
    t_start = time.time()
    print("  routing pinned at truth; optpercolation=%d; tag=%s; scored against %s"
          % (optperc, tag, ctx.ref_path.name), flush=True)

    for i, srow in samples.iterrows():
        ks, f, cc = float(srow["Ks_mult"]), float(srow["f_RS_abs"]), float(srow[CC_PARAM])
        pred_id = predicted_run_id(cc, tag)
        print("\n[%4d/%d]  Ks=%.3fx  f=%.4f  cc=%.1f mm/hr%s  -> %s"
              % (i + 1, args.n, ks, f, cc, " (inert)" if control else "", pred_id), flush=True)

        csv_ok = (ctx.csv_dir / ("%s_compare_obs_sim.csv" % pred_id)).exists()
        metrics_file = ctx.summary_dir / ("%s_metrics_summary.csv" % pred_id)
        if args.skip_existing and csv_ok and (pred_id in done_ids or metrics_file.exists()):
            print("  SKIP (compare CSV exists): %s" % pred_id)
            skipped += 1
            if pred_id not in done_ids:
                try:
                    m = pd.read_csv(metrics_file).iloc[0].to_dict()
                    m.setdefault("Ks_mult", ks)
                    m.setdefault("f_RS_abs", f)
                    m.setdefault(CC_PARAM, cc)
                    m["optpercolation"] = optperc
                    m["storm_label"], m["rain_scale"] = ctx.label, ctx.scale
                    m["reference_truth"] = ctx.ref_path.name
                    rows.append(m)
                    done_ids.add(pred_id)
                except Exception:
                    pass
            continue

        overrides = dict(ROUTING_TRUTH)
        overrides.update({"Ks_mult": ks, "f_RS_abs": f, "optpercolation": optperc})
        expect = {"Ks_mult": ks, "f_RS_abs": f, "cc": cc, "optpercolation": optperc}
        res = run_one(ctx, pred_id, cc, overrides, tag, truth_rel, expect, ctx.ref_path)

        if res["status"] == "SUCCESS":
            m = res["metrics"]
            m.setdefault("Ks_mult", ks)
            m.setdefault("f_RS_abs", f)
            m.setdefault(CC_PARAM, cc)
            m["optpercolation"] = optperc
            m["storm_label"], m["rain_scale"] = ctx.label, ctx.scale
            m["reference_truth"] = ctx.ref_path.name
            rows = [r for r in rows if r.get("run_id") != pred_id] + [m]
            done_ids.add(pred_id)
            completed += 1
            note_outcome(ctx, True)
            print("  %s: KGE_2012=%.3f  PBIAS=%+.1f%%"
                  % ("CONTROL (cc OFF)" if control else "cc ON", _num(m.get("kge_2012")), _num(m.get("pbias_pct"))))
            if not args.no_cleanup:
                if cleanup_raw_results(ctx.calib_dir / "02_results" / ctx.category / pred_id):
                    cleaned += 1
        else:
            print("  %s: %s  (%s)" % (res["status"], pred_id, res["reason"]))
            if res["status"] == "HANG":
                hung += 1
            else:
                failed += 1
            failed_rows.append({"run_id": pred_id, "status": res["status"], "reason": res["reason"],
                                "elapsed_min": res["elapsed"] / 60, "Ks_mult": ks, "f_RS_abs": f,
                                CC_PARAM: cc})
            pd.DataFrame(failed_rows).to_csv(failed_path, index=False)
            note_outcome(ctx, False)

        remaining = args.n - completed - skipped - hung - failed
        if completed > 0:
            eta_min = (time.time() - t_start) / completed * remaining / 60
            print("  Run time: %.1f min  |  ETA: %.0f min (%.1f h) remaining"
                  % (res["elapsed"] / 60, eta_min, eta_min / 60), flush=True)
        if rows:
            pd.DataFrame(rows).to_csv(out_path, index=False)      # save after every run

    print("\nSweep complete: %d ran, %d skipped, %d hung, %d failed, %d raw dirs cleaned"
          % (completed, skipped, hung, failed, cleaned))
    if rows:
        final_df = pd.DataFrame(rows)
        if "kge_2012" in final_df.columns:
            final_df = final_df.sort_values("kge_2012", ascending=False)
        final_df.to_csv(out_path, index=False)
        print("Saved: %s  (%d rows)" % (out_path, len(final_df)))
        try:
            report_sweep(ctx, final_df, control, on_path, off_path)
        except Exception as e:
            print("\n(diagnostics skipped: %s -- results are saved regardless)" % e)
    else:
        print("No results to save.")
    n_done = len(done_ids)
    if n_done < args.n:
        print("\n%d of %d sweep points are not complete. Re-run the same command with --skip_existing to retry "
              "them (see %s)." % (args.n - n_done, args.n, failed_path.name))


# ------------------------------------------------------------------
# main
# ------------------------------------------------------------------
class Ctx(object):
    pass


def check_pandas_note(summary_dir):
    """Soft check: is this python's pandas the one that scored Stage 2 (per the rescoring provenance)?"""
    p = summary_dir / "rescored_truth_location_110" / "PROVENANCE_rescored_truth_location_110.json"
    try:
        old = json.loads(p.read_text()).get("pandas")
    except Exception:
        return
    if not old:
        return
    if old.split(".")[:2] != pd.__version__.split(".")[:2]:
        print("  *** NOTE: pandas here is %s but the earlier rescoring recorded %s. The 5-minute interpolation "
              "differs slightly between pandas 2 and 3, so scores could differ from Stage 2's at the 1e-3 level. "
              "Use the same python you used before. ***" % (pd.__version__, old))
    else:
        print("  pandas %s matches the version recorded by the earlier rescoring." % pd.__version__)


def main():
    parser = argparse.ArgumentParser(
        description="Series 110 -- a second/third storm (rain scaled), channel loss ON: scaled forcing, "
                    "reference, 27 truths, 250-run sweep.")
    parser.add_argument("--rain_scale", type=float, default=0.8,
                        help="Multiplier on the rainfall rate (default 0.8 = the second storm; use 1.25 for the third)")
    parser.add_argument("--only", choices=["truths", "sweep"], default=None,
                        help="truths = forcing + reference + truths, no sweep; sweep = just the cc-ON sweep "
                             "(needs the reference from an earlier run)")
    parser.add_argument("--control", action="store_true",
                        help="Run only the CONTROL sweep (same 250 points, channel loss OFF); needs the reference")
    parser.add_argument("--n", type=int, default=250, help="Sweep size (default 250 -- must equal Stage 2's n)")
    parser.add_argument("--seed", type=int, default=42, help="LHS seed (default 42 -- must equal Stage 2's seed)")
    parser.add_argument("--ks_values", type=float, nargs="+", default=DEFAULT_KS)
    parser.add_argument("--f_values", type=float, nargs="+", default=DEFAULT_F)
    parser.add_argument("--cc_values", type=float, nargs="+", default=DEFAULT_CC)
    parser.add_argument("--skip_existing", action="store_true",
                        help="Skip work already completed (use this every time; it is also how you resume)")
    parser.add_argument("--timeout", type=int, default=300, help="Per-run hard timeout in seconds (default 300)")
    parser.add_argument("--no_cleanup", action="store_true", help="Keep each run's raw results directory")
    parser.add_argument("--dry_run", action="store_true", help="Run the checks and print the plan; write and run nothing")
    parser.add_argument("--ignore_running_check", action="store_true",
                        help="Skip the check for other tRIBS scripts running (only if you are sure none are)")
    parser.add_argument("--overwrite_forcing", action="store_true",
                        help="Replace existing scaled rain files that differ from what this script writes")
    parser.add_argument("--allow_unpaired", action="store_true",
                        help="Continue even if the sweep points do not pair with Stage 2's")
    parser.add_argument("--max_consecutive_failures", type=int, default=8,
                        help="Stop if this many runs in a row fail (default 8)")
    parser.add_argument("--src_sdf", default=SRC_SDF_DEFAULT,
                        help="Original rain gauge SDF, relative to lab/ (default %s)" % SRC_SDF_DEFAULT)
    args = parser.parse_args()

    if args.control and args.only == "truths":
        parser.error("--control and --only truths cannot be combined")
    try:
        label_, short = storm_names(args.rain_scale)
        ks_vals = validate_axis("Ks", args.ks_values, "Ks_mult")
        f_vals = validate_axis("f", args.f_values, "f_RS_abs")
        cc_vals = validate_axis("cc", args.cc_values, "cc_mmhr")
    except ValueError as e:
        parser.error(str(e))
    if args.n < 1:
        parser.error("--n must be at least 1")

    ctx = Ctx()
    ctx.args, ctx.label, ctx.short, ctx.scale = args, label_, short, float(args.rain_scale)
    ctx.consec_fail = 0
    cwd = Path.cwd()
    ctx.calib_dir = cwd.parent / "calibration_work"
    ctx.summary_dir = ctx.calib_dir / "03_comparisons" / "summary_tables"
    ctx.csv_dir = ctx.calib_dir / "03_comparisons" / "csv_exports"
    synth_dir = ctx.calib_dir / "synth_truth"
    ctx.ref_dir = ctx.calib_dir / ("synth_truth_ccoff_%s" % label_)
    ctx.loc_dir = ctx.calib_dir / ("synth_truth_location_%s" % label_)
    ctx.ref_csv = ctx.summary_dir / ("storm_ref_110_%s.csv" % label_)
    ctx.truth_csv = ctx.summary_dir / ("truth_location_110_%s.csv" % label_)
    ctx.truth_failed_csv = ctx.summary_dir / ("truth_location_FAILED_110_%s.csv" % label_)
    ctx.gauge_sdf_rel = derived_name(args.src_sdf, label_)
    ctx.category = builder.get_run_category(builder.PARAM_CONFIG[CC_PARAM]["series"])
    ctx.truths = make_truths(ks_vals, f_vals, cc_vals, short)
    ctx.ref_path = None

    ids = [t["run_id"] for t in ctx.truths]
    if len(set(ids)) != len(ids):
        raise RuntimeError("two truths map to the same run_id (labels collide); use more widely spaced values")

    print("Series 110 -- storm run: %s  (rain x%g)   python %s, pandas %s"
          % (label_, ctx.scale, sys.version.split()[0], pd.__version__))

    # ------------------------------------------------------------------
    # HARD SAFETY CHECKS
    # ------------------------------------------------------------------
    if "gauge_sdf" not in inspect.signature(builder.build_input_file).parameters:
        raise RuntimeError("build_sensitivity_run_tribs6.build_input_file has no gauge_sdf argument, so it cannot "
                           "be pointed at scaled rain. This script needs the tRIBS-6 builder that has it.")
    qout_files = list(synth_dir.glob("*.qout")) if synth_dir.exists() else []
    if len(qout_files) != 1:
        raise RuntimeError(
            "Expected exactly one *.qout file in %s (the ORIGINAL storm's cc-OFF truth), found %d: %s. Without it "
            "the scorer would silently fall back to real-gauge mode." % (synth_dir, len(qout_files), [f.name for f in qout_files]))
    ctx.storm1_ref = qout_files[0]
    print("Original storm's cc-OFF truth, read-only: %s  [md5 %s]" % (ctx.storm1_ref.name, md5_of(ctx.storm1_ref)[:12]))

    synth_res = synth_dir.resolve()
    for d in (ctx.ref_dir, ctx.loc_dir):
        dres = d.resolve()
        if dres == synth_res or synth_res in dres.parents:
            raise RuntimeError("truth folder %s is, or is inside, %s; refusing." % (d, synth_dir))
    if ctx.ref_dir.resolve() == ctx.loc_dir.resolve():
        raise RuntimeError("reference and truth folders must differ")

    if args.ignore_running_check:
        print("Running-jobs check skipped (--ignore_running_check).")
    else:
        competing = find_competing_jobs()
        if competing is None:
            print("WARNING: could not check for other running tRIBS scripts (ps unavailable). Make sure none are "
                  "running -- they share current_run_config.json.")
        elif competing:
            print("\nANOTHER tRIBS BUILD/RUN SCRIPT APPEARS TO BE RUNNING:")
            for pid, cmd in competing:
                print("  pid %d: %s" % (pid, cmd[:110]))
            raise RuntimeError(
                "Refusing to start: these scripts share current_run_config.json with no locking, so running both "
                "can corrupt both. Wait for them to finish (or stop them with pkill -f <script>), or re-run with "
                "--ignore_running_check if you are sure they are stale.")
        else:
            print("No other tRIBS build/run script detected.")
    check_pandas_note(ctx.summary_dir)

    # ------------------------------------------------------------------
    # Plan the forcing (in memory) and print the whole plan
    # ------------------------------------------------------------------
    ctx.forcing_plan = plan_forcing(args.src_sdf, label_, ctx.scale)
    if args.control:
        stages = ["forcing", "reference_check", "sweep_control"]
    elif args.only == "sweep":
        stages = ["forcing", "reference_check", "sweep"]
    elif args.only == "truths":
        stages = ["forcing", "reference", "truths"]
    else:
        stages = ["forcing", "reference", "truths", "sweep"]

    n_truth_runs = 1 + len(ctx.truths)
    est = (n_truth_runs if "truths" in stages else 0) + (args.n if any(s.startswith("sweep") for s in stages) else 0)
    print("\n" + "=" * 92)
    print("PLAN: %s -- rain x%g.   stages: %s" % (label_, ctx.scale, ", ".join(stages)))
    print("  scorer 'observed' for the reference: the original storm's truth (gives the cross-storm gate);")
    print("  scorer 'observed' for truths and sweep: this storm's reference (set via truth_file in the run config).")
    print("  routing at truth: cv=%g r=%g n=%g.  tags: ref_%s, K<Ks>_f<f>_%s, sw_%s, ctl_%s"
          % (ROUTING_TRUTH["kinemvelcoef"], ROUTING_TRUTH["flowexp"], ROUTING_TRUTH["channelroughness"],
             short, short, short, short))
    print("  about %d tRIBS runs, ~1 minute each = roughly %.1f hours" % (est, est / 60.0))
    print("  stops by itself if: the observed file is wrong, the builder config is wrong, the reference gate "
          "fails, or %d runs in a row fail." % args.max_consecutive_failures)
    print("=" * 92)
    print_forcing_plan(ctx.forcing_plan)
    if "truths" in stages:
        print("  truths (%d), centre first:" % len(ctx.truths))
        print("  %4s %6s %8s %8s   run_id" % ("#", "Ks", "f", "cc"))
        for i, t in enumerate(ctx.truths, start=1):
            print("  %4d %6g %8g %8g   %s" % (i, t["Ks_mult"], t["f_RS_abs"], t["cc_mmhr"], t["run_id"]))
    if any(s.startswith("sweep") for s in stages):
        print("  sweep: n=%d seed=%d  Ks %g-%g, f %g-%g (log), cc %g-%g (log)"
              % (args.n, args.seed, LHS_PARAMS["Ks_mult"]["lo"], LHS_PARAMS["Ks_mult"]["hi"],
                 LHS_PARAMS["f_RS_abs"]["lo"], LHS_PARAMS["f_RS_abs"]["hi"],
                 LHS_PARAMS[CC_PARAM]["lo"], LHS_PARAMS[CC_PARAM]["hi"]))
    for p in (ctx.ref_csv, ctx.truth_csv):
        if p.exists():
            print("  NOTE: %s exists; %s" % (p.name, "--skip_existing will keep what is complete."
                                            if args.skip_existing else "without --skip_existing those runs are RE-RUN."))

    if args.dry_run:
        print("\n--dry_run: nothing written, built or run.")
        return

    ctx.summary_dir.mkdir(parents=True, exist_ok=True)
    ctx.csv_dir.mkdir(parents=True, exist_ok=True)
    ctx.truth_prov_base = {
        "script": Path(__file__).name,
        "created_local": datetime.now().isoformat(timespec="seconds"),
        "what_this_is": "Synthetic cc-ON truths (outlet .qout) under rainfall scaled by %g (%s), one per (Ks, f, cc) "
                        "combination. Same grid as run_truth_location_110.py." % (ctx.scale, label_),
        "storm_label": label_, "rain_scale": ctx.scale,
        "routing_truth": ROUTING_TRUTH, "optpercolation": 1, "optsnow": 0,
        "cc_units": "mm/hr (assigned directly to CHANNELCONDUCTIVITY)",
        "grid": {"Ks_mult": ks_vals, "f_RS_abs": f_vals, "cc_mmhr": cc_vals},
        "sampled_box": {k: list(v[:2]) for k, v in BOX.items()},
        "gauge_sdf": ctx.gauge_sdf_rel,
        "anchor_status": "TEST TRUTHS ONLY -- not calibration anchors",
        "python": sys.version.split()[0], "pandas": pd.__version__, "numpy": np.__version__,
    }

    try:
        for s in stages:
            if s == "forcing":
                stage_forcing(ctx)
            elif s == "reference":
                stage_reference(ctx, allow_run=True)
            elif s == "reference_check":
                stage_reference(ctx, allow_run=False)
            elif s == "truths":
                ctx.truth_prov_base["reference_truth_ccoff"] = {"file": ctx.ref_path.name, "md5": ctx.ref_md5}
                stage_truths(ctx)
            elif s == "sweep":
                stage_sweep(ctx, control=False)
            elif s == "sweep_control":
                stage_sweep(ctx, control=True)
    except Abort as e:
        print("\n" + "!" * 92)
        print("STOPPED: %s" % e)
        print("!" * 92)
        print("Everything finished before this point is kept. After fixing the cause, re-run the same command "
              "with --skip_existing.")
        sys.exit(3)
    except KeyboardInterrupt:
        print("\nStopped by you. Finished runs are kept; re-run the same command with --skip_existing to continue.")
        sys.exit(130)

    print("\n" + "=" * 92)
    print("ALL REQUESTED STAGES FINISHED for %s." % label_)
    if "sweep" in stages:
        print("Optional next: the same 250 points with channel loss OFF:   python -u run_cc_storm_110.py "
              "--rain_scale %g --control --skip_existing 2>&1 | tee -a %s_ctl_110.log" % (ctx.scale, label_))
    print("Next (separate script, not written yet): re-score the stored runs against these truths.")
    print("=" * 92)


if __name__ == "__main__":
    main()
