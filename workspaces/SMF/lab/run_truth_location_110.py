"""
run_truth_location_110.py
=========================
Series 110 -- the TRUTH-LOCATION test, part 1: build synthetic cc-ON truths
whose TRUE Ks and f are NOT the centre of the sampled box.

WHY THIS EXISTS (plain language)
----------------------------------
Every cc-ON truth used so far was made at the same Ks and f (Ks_mult = 7.0,
f_RS_abs = 0.012), and the 250-run Stage 2 sweep was drawn over a box
(Ks 4-10.5, f 0.004-0.030, cc 30-1000) whose middle sits right on that point.
That is the easiest possible case: the answer is in the middle of the search
region. A skeptic will ask: "does cc still get recovered when the true Ks and
f are somewhere else?"

This script makes the new answer keys. It runs tRIBS once per truth, on a
grid of true Ks x true f x true cc, and keeps each run's outlet hydrograph.
A second script (rescore_truth_location_110.py) then re-scores the 250 + 250
runs you ALREADY HAVE against every one of these truths -- no old run is
repeated -- and a third (analyze_truth_location_110.py) asks whether the best
runs recover the cc, Ks and f that made each truth.

REMINDER OF THE IDEA: the 250 saved runs are candidate GUESSES (each with a
stored hydrograph). A truth is only the answer key they are scored against.
Each new answer key costs one tRIBS run; re-scoring against it takes seconds.

WHAT IS HELD FIXED
-------------------
  routing:   kinemvelcoef (cv) = 4.5, flowexp (r) = 0.24, channelroughness (n)
             = 0.026  -- the SAME routing values every Stage 2 run used. A truth
             with different routing would be a different test (not done here).
  channel loss: optpercolation = 1 (ON), optsnow = 0 (set in the builder)
  storm and forcing: the same August 12 2014 event and forcing as everything
             else in Series 110.
  Every truth lies INSIDE the sampled box (the script refuses otherwise): a
  truth outside the box could never be recovered by a sweep that never looked
  there, which would say nothing about cc.

DEFAULT GRID (27 truths, about 27-30 minutes at ~1 minute per tRIBS run)
-------------------------------------------------------------------------
    true Ks_mult:  5.0   7.0   9.0
    true f_RS_abs: 0.007 0.012 0.020
    true cc:       95.1  201   425      (mm/hr)
The Ks 7.0 / f 0.012 row is the ORIGINAL centre. Its three cc values (95.1, 201,
425) are exactly three of the 12 candidate truths made earlier, so those three
runs are also a reproducibility check: today's script must rebuild the same
hydrographs. The end-of-run report compares them with the earlier files.
Runs are done centre first, then outward, so an interrupted run still leaves
the most useful truths.

WHAT YOU GET
-------------
  calibration_work/synth_truth_location/
      SMF_20140812_110_cc<cc>_loc_K<Ks>_f<f>_Outlet.qout   (one per truth)
      PROVENANCE_truth_location_110.json                   (exact values + checksums)
  calibration_work/03_comparisons/summary_tables/
      truth_location_110.csv           one row per truth
      truth_location_FAILED_110.csv    only if something failed

Scorer metrics in the CSV (pbias_pct, kge_2012, ...) are measured against the
current cc-OFF truth in calibration_work/synth_truth/. They are NOT recovery
scores; they just show how far each new truth sits from the old one.

WHY THE FOLDER IS SEPARATE (safety)
-------------------------------------
Every sweep script refuses to run unless calibration_work/synth_truth/ holds
exactly ONE *.qout, and the scorer uses that file as "observed". New truths are
kept in their own folder and never go there. This script only READS synth_truth/
and refuses a truth folder that is, or sits inside, synth_truth/.

SAFETY (same rules as every other Series 110 script)
------------------------------------------------------
  - Needs exactly one *.qout in calibration_work/synth_truth/ (reference metrics
    only); refuses to run otherwise.
  - run_id tag is loc_K<Ks>_f<f>, so nothing from Stage 1, Stage 2, the control,
    the provenance check or the earlier cc candidates can be touched.
  - Refuses to start if another tRIBS build/run script appears to be running
    (they all share calibration_work/current_run_config.json with no locking).
    NEVER run two of these at once.
  - Each run is timeout-safe; failures are logged to the FAILED csv.
  - After each run the outlet file is copied, checksum-verified and
    sanity-checked BEFORE that run's raw output is deleted. --no_cleanup keeps
    everything.

USAGE (run from the lab/ directory, like every other script)
---------------------------------------------------------------
    python run_truth_location_110.py --dry_run          # look, don't run
    python -u run_truth_location_110.py 2>&1 | tee truth_location_110.log
    python -u run_truth_location_110.py --skip_existing 2>&1 | tee -a truth_location_110.log
                                                        # resume after an interruption

    python run_truth_location_110.py --ks_values 5 9 --f_values 0.007 0.02 \\
                                     --cc_values 95.1 425           # a smaller grid
    python run_truth_location_110.py --timeout 300

Overnight / unattended: keep the laptop awake and the Codespaces browser tab
open, and always use `python -u ... 2>&1 | tee file` so output keeps flowing
(Codespaces judges "idle" by terminal activity). If it is interrupted, rerun
the same command with --skip_existing; finished truths are kept.
"""

import argparse
import hashlib
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
# Truth values -- the SAME ones used in Stage 1, Stage 2, the control and the
# earlier cc candidates. Passed explicitly on every build (BASELINE is stale).
# Ks_mult and f_RS_abs are the ones this script VARIES; routing never changes.
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

TAG_PREFIX = "loc"
RESERVED_TAGS = {"s1", "s2", "ctl", "truthcheck", "volmatch", "cctruth"}

# Stage 2's sampled box. A truth must lie inside it for a recovery test to be fair.
BOX = {
    "Ks_mult":  (4.0,   10.5,   "lin"),
    "f_RS_abs": (0.004, 0.030,  "log"),
    "cc_mmhr":  (30.0,  1000.0, "log"),
}

DEFAULT_KS = [5.0, 7.0, 9.0]
DEFAULT_F = [0.007, 0.012, 0.020]
DEFAULT_CC = [95.1, 201.0, 425.0]

CAND_DIRNAME = "synth_truth_location"
OLD_CAND_DIRNAME = "synth_truth_cc_candidates"     # the earlier 12 cc candidates
OLD_CAND_TAG = "cctruth"
OUT_NAME = "truth_location_110.csv"
FAILED_NAME = "truth_location_FAILED_110.csv"
PROV_NAME = "PROVENANCE_truth_location_110.json"

# Rule-of-thumb flags (adjust freely).
EDGE_FRAC = 0.10            # within this fraction of a box end on any axis
NEAR_DRY_VOL_RATIO = 0.10   # truth volume / cc-OFF truth volume below this -> NEAR_DRY
MONO_TOL = 1e-9             # relative slack before a volume rise is called a reversal
REPRO_REL = 1e-6            # max |difference| / peak below this counts as "same hydrograph"

# Other scripts that share current_run_config.json (checked before starting).
COMPETING_PATTERNS = (
    "run_sensitivity_single", "run_sensitivity_sweep", "build_sensitivity_run",
    "run_cc_", "run_lhs_", "run_crossversion", "verify_truth_provenance",
    "validate_truth_point", "rescore_series100", "run_truth_",
)

EXPECTED_QOUT_COLS = 3      # Time_hr, Qstrm_m3_s, Hlev_m


# ------------------------------------------------------------------
# Grid
# ------------------------------------------------------------------
def label(v):
    """The builder's own number-to-text rule (7.0 -> 7p0, 0.012 -> 0p012)."""
    return builder.value_to_label(v)


def label_roundtrips(v):
    """True if the label spells the value exactly (the builder keeps 6 decimals)."""
    try:
        return math.isclose(float(label(v).replace("p", ".")), float(v), rel_tol=1e-12, abs_tol=0.0)
    except ValueError:
        return False


def box_pos(v, key):
    """0 at the low end of the sampled range, 1 at the high end (Ks linear, f and cc log)."""
    lo, hi, scale = BOX[key]
    if scale == "log":
        return float((np.log10(v) - np.log10(lo)) / (np.log10(hi) - np.log10(lo)))
    return float((v - lo) / (hi - lo))


def tag_for(ks, f):
    return "%s_K%s_f%s" % (TAG_PREFIX, label(ks), label(f))


def predicted_run_id(ks, f, cc):
    """The run_id the builder will produce (computed first so resume can skip a
    finished truth without rebuilding it)."""
    base, _ = builder.build_run_id(CC_PARAM, cc)
    return "%s_%s" % (base, tag_for(ks, f))


def old_candidate_run_id(cc):
    base, _ = builder.build_run_id(CC_PARAM, cc)
    return "%s_%s" % (base, OLD_CAND_TAG)


def validate_axis(name, values, key):
    """Sorted, de-duplicated values for one axis; raises ValueError with a plain message."""
    vals = [float(v) for v in values]
    if not vals:
        raise ValueError("no %s values given" % name)
    if any((not np.isfinite(v)) or v <= 0 for v in vals):
        raise ValueError("%s values must be finite and > 0: %s" % (name, vals))
    lo, hi, _ = BOX[key]
    outside = [v for v in vals if not (lo <= v <= hi)]
    if outside:
        raise ValueError(
            "%s value(s) %s lie outside the sampled box %g-%g. A truth outside the box can "
            "never be recovered by a sweep that did not look there, so it would say nothing "
            "about cc. Pick values inside the box." % (name, ", ".join("%g" % v for v in outside), lo, hi))
    bad = [v for v in vals if not label_roundtrips(v)]
    if bad:
        raise ValueError("%s value(s) %s need more than 6 decimal places, which the run labels "
                         "cannot spell exactly. Round them." % (name, ", ".join("%r" % v for v in bad)))
    uniq = sorted(set(vals))
    if len(uniq) != len(vals):
        raise ValueError("%s has repeated values: %s" % (name, vals))
    return uniq


def make_truths(ks_vals, f_vals, cc_vals):
    """Every (Ks, f, cc) combination, centre-out. Each item is a dict."""
    c_k = box_pos(STANDARD_KS, "Ks_mult")
    c_f = box_pos(STANDARD_F, "f_RS_abs")
    items = []
    for ks in ks_vals:
        for f in f_vals:
            for cc in cc_vals:
                pk, pf, pc = box_pos(ks, "Ks_mult"), box_pos(f, "f_RS_abs"), box_pos(cc, "cc_mmhr")
                items.append({
                    "Ks_mult": ks, "f_RS_abs": f, "cc_mmhr": cc,
                    "tag": tag_for(ks, f), "run_id": predicted_run_id(ks, f, cc),
                    "pos": (pk, pf, pc),
                    "cell_dist": float(np.hypot(pk - c_k, pf - c_f)),
                })
    items.sort(key=lambda d: (round(d["cell_dist"], 9), d["Ks_mult"], d["f_RS_abs"], d["cc_mmhr"]))
    return items


# ------------------------------------------------------------------
# Small helpers (ported unchanged from run_cc_truth_candidates_110.py)
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


def _ancestor_pids():
    """This process and its parents (so a wrapping shell is not mistaken for a
    competing job). Linux /proc only; returns what it can."""
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


def run_with_timeout(timeout_sec):
    """Timeout-safe subprocess execution -- same pattern as every other Series 110
    script. Own process group so a hang is killed as a unit."""
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
        # Ctrl-C: the scorer runs in its own process group, so the terminal would
        # not stop it. Kill it here so nothing is left running and a restart is
        # not refused by the "another tRIBS script is running" check.
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.communicate()
        print("\nInterrupted: the running tRIBS job was stopped. Finished truths are kept; "
              "re-run with --skip_existing to continue.")
        raise
    elapsed = time.time() - t0

    if stdout:
        print(stdout, end="" if stdout.endswith("\n") else "\n")

    tribs_warning = bool(stdout) and "WARNING: tRIBS may have failed" in stdout
    returncode = None if timed_out else proc.returncode
    return returncode, elapsed, timed_out, tribs_warning


def read_qout(path):
    """Parse a tRIBS outlet file by hand. Real layout (tRIBS 6.0.0): one header line
    'Time_hr,Qstrm_m3_s,Hlev_m', then COMMA-separated rows. Commas and/or whitespace
    are accepted, and a header line is skipped when its first token is not a number.
    Returns (time_hr, q_m3s, column_counts_seen, n_unparseable_lines)."""
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
    """Sanity-check a .qout file. Never raises. qout_ok is True only if nothing
    looked wrong."""
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
            problems.append("no flow at all")

    info["qout_problem"] = "; ".join(problems)
    info["qout_ok"] = not problems
    return info


def locate_raw_qout(raw_dir, run_id):
    """Where tRIBS wrote this run's outlet file. Expected name first; else a single
    *.qout in the run's folder; else None."""
    expected = raw_dir / ("%s_Outlet.qout" % run_id)
    if expected.exists():
        return expected
    found = list(raw_dir.glob("*.qout")) if raw_dir.exists() else []
    return found[0] if len(found) == 1 else None


def copy_verified(src, dst):
    """Copy src -> dst via a temporary name, verify size and checksum, then move into
    place. Returns the md5. Raises on any mismatch (no partial file under the final name)."""
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
    """Delete a run's raw results directory (pixel files, spatial snapshots). The .in
    file, logs and the two CSVs live elsewhere. Never raises."""
    try:
        if raw_dir.exists():
            shutil.rmtree(raw_dir)
            return True
    except Exception as e:
        print("  (cleanup warning: could not remove %s: %s)" % (raw_dir, e))
    return False


# ------------------------------------------------------------------
# Row assembly, flags, provenance
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


def assemble_row(t, run_id, metrics, cand_path, md5, qinfo, ref_truth_name):
    obs_vol = _num(metrics.get("obs_volume_m3"))
    sim_vol = _num(metrics.get("sim_volume_m3"))
    obs_peak = _num(metrics.get("obs_peak_m3s"))
    sim_peak = _num(metrics.get("sim_peak_m3s"))
    vol_ratio = sim_vol / obs_vol if obs_vol and np.isfinite(obs_vol) else np.nan
    peak_ratio = sim_peak / obs_peak if obs_peak and np.isfinite(obs_peak) else np.nan
    pk, pf, pc = t["pos"]
    row = {
        "run_id":                run_id,
        "Ks_mult":               t["Ks_mult"],
        "f_RS_abs":              t["f_RS_abs"],
        "cc_mmhr":               t["cc_mmhr"],
        "box_pos_Ks":            pk,
        "box_pos_f":             pf,
        "box_pos_cc":            pc,
        "flags":                 make_flags(t["pos"], vol_ratio, qinfo["qout_ok"]),
        "reference_truth_ccoff": ref_truth_name,
        "pbias_pct":             _num(metrics.get("pbias_pct")),
        "kge_2012":              _num(metrics.get("kge_2012")),
        "vol_ratio_vs_ccoff":    vol_ratio,
        "peak_ratio_vs_ccoff":   peak_ratio,
        "candidate_qout":        cand_path.name,
        "candidate_md5":         md5,
        "candidate_bytes":       int(cand_path.stat().st_size),
    }
    row.update(qinfo)
    for k, v in metrics.items():          # all remaining scorer columns, unchanged
        row.setdefault(k, v)
    row["optpercolation"] = 1
    return row


def write_provenance(cand_dir, base, rows):
    prov = dict(base)
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
    tmp = cand_dir / (PROV_NAME + ".part")
    tmp.write_text(json.dumps(prov, indent=2, default=str))
    os.replace(tmp, cand_dir / PROV_NAME)


# ------------------------------------------------------------------
# End-of-run report
# ------------------------------------------------------------------
def reproduction_check(df, calib_dir):
    """For truths at the ORIGINAL centre (Ks 7.0, f 0.012) whose cc matches one of the
    earlier cc candidates, compare the new file with the old one. Returns a list of
    (cc, verdict_text, ok_bool)."""
    old_dir = calib_dir / OLD_CAND_DIRNAME
    new_dir = calib_dir / CAND_DIRNAME
    out = []
    centre = df[np.isclose(df["Ks_mult"].astype(float), STANDARD_KS) &
                np.isclose(df["f_RS_abs"].astype(float), STANDARD_F)]
    for _, r in centre.sort_values("cc_mmhr").iterrows():
        cc = float(r["cc_mmhr"])
        old = old_dir / ("%s_Outlet.qout" % old_candidate_run_id(cc))
        new = new_dir / str(r["candidate_qout"])
        if not old.exists() or not new.exists():
            continue
        try:
            if md5_of(old) == md5_of(new):
                out.append((cc, "identical file (same bytes)", True))
                continue
            t0, q0, _, _ = read_qout(old)
            t1, q1, _, _ = read_qout(new)
            if len(q0) != len(q1) or not np.allclose(t0, t1):
                out.append((cc, "DIFFERENT: the two files do not have the same time grid", False))
                continue
            peak = float(np.max(np.abs(q0))) or 1.0
            d = float(np.max(np.abs(q0 - q1)))
            ok = d <= REPRO_REL * peak
            out.append((cc, ("same hydrograph; largest difference %.2e m3/s (%.1e of the peak)" % (d, d / peak))
                        if ok else
                        ("DIFFERENT: largest difference %.3g m3/s (%.1e of the peak)" % (d, d / peak)), ok))
        except Exception as e:
            out.append((cc, "could not compare (%s)" % str(e)[:80], False))
    return out


def report(df, ref_truth_name, calib_dir):
    print("\n" + "=" * 100)
    print("NEW cc-ON TRUTHS AT DIFFERENT TRUE Ks AND f  (metrics are vs the cc-OFF truth: %s)" % ref_truth_name)
    print("=" * 100)
    d = df.sort_values(["Ks_mult", "f_RS_abs", "cc_mmhr"]).reset_index(drop=True)
    d["flags"] = d["flags"].fillna("") if "flags" in d else ""
    print("%6s %8s %8s | %9s %9s %9s %10s  flags"
          % ("Ks", "f", "cc", "PBIAS %", "KGE_2012", "vol kept", "peak kept"))
    last = None
    for _, r in d.iterrows():
        cell = (r["Ks_mult"], r["f_RS_abs"])
        if last is not None and cell != last:
            print("")
        last = cell
        vk, pk = _num(r.get("vol_ratio_vs_ccoff")), _num(r.get("peak_ratio_vs_ccoff"))
        print("%6g %8g %8g | %+9.2f %9.3f %8.1f%% %9.1f%%  %s"
              % (r["Ks_mult"], r["f_RS_abs"], r["cc_mmhr"], _num(r.get("pbias_pct")),
                 _num(r.get("kge_2012")), 100 * vk if np.isfinite(vk) else np.nan,
                 100 * pk if np.isfinite(pk) else np.nan, r["flags"]))
    print("  (vol kept / peak kept = this truth relative to the cc-OFF truth, in the event window. "
          "They move with Ks and f as well as cc, so they are NOT a measure of the cc effect alone.)")

    # Sanity: more cc should remove more water at fixed Ks and f.
    rev_cc = []
    for (ks, f), g in d.groupby(["Ks_mult", "f_RS_abs"]):
        g = g.sort_values("cc_mmhr")
        v = g["sim_volume_m3"].astype(float).to_numpy() if "sim_volume_m3" in g else np.array([])
        cc = g["cc_mmhr"].to_numpy(float)
        for i in range(len(v) - 1):
            if np.isfinite(v[i]) and np.isfinite(v[i + 1]) and v[i + 1] > v[i] * (1 + MONO_TOL):
                rev_cc.append((ks, f, cc[i], cc[i + 1]))
    if "sim_volume_m3" in d:
        if rev_cc:
            print("\n  CHECK cc: outlet volume ROSE as cc rose at %d step(s) (it should fall):" % len(rev_cc))
            for ks, f, a, b in rev_cc[:6]:
                print("    Ks %g, f %g: cc %g -> %g" % (ks, f, a, b))
        else:
            print("\n  CHECK cc: at every Ks and f, outlet volume falls (or holds) as cc rises.")
        rev_ks = []
        for (f, cc), g in d.groupby(["f_RS_abs", "cc_mmhr"]):
            g = g.sort_values("Ks_mult")
            v = g["sim_volume_m3"].astype(float).to_numpy()
            ks = g["Ks_mult"].to_numpy(float)
            for i in range(len(v) - 1):
                if np.isfinite(v[i]) and np.isfinite(v[i + 1]) and v[i + 1] > v[i] * (1 + MONO_TOL):
                    rev_ks.append((f, cc, ks[i], ks[i + 1]))
        if rev_ks:
            print("  CHECK Ks: outlet volume ROSE as Ks rose at %d step(s) (it should fall):" % len(rev_ks))
            for f, cc, a, b in rev_ks[:6]:
                print("    f %g, cc %g: Ks %g -> %g" % (f, cc, a, b))
        else:
            print("  CHECK Ks: at every f and cc, outlet volume falls (or holds) as Ks rises.")

    rep = reproduction_check(d, calib_dir)
    if rep:
        print("\n  REPRODUCIBILITY: truths at the original centre (Ks %g, f %g) against the earlier cc candidates:"
              % (STANDARD_KS, STANDARD_F))
        for cc, text, ok in rep:
            print("    cc %-6g %s  %s" % (cc, "OK  " if ok else "FAIL", text))
        if not all(ok for _, _, ok in rep):
            print("    A failure here means today's build does not reproduce the earlier truth. Do not use the "
                  "new truths until that is explained (the re-scoring script will stop on it too).")
    else:
        print("\n  REPRODUCIBILITY: no truth at the original centre matches an earlier cc candidate "
              "(nothing to compare).")

    n_bad = int((~d["qout_ok"].astype(bool)).sum()) if "qout_ok" in d else 0
    if n_bad:
        print("\n  WARNING: %d truth file(s) failed the sanity check (see qout_problem in the CSV)." % n_bad)
    if "obs_mode" in d and (d["obs_mode"] != "synth").any():
        print("  WARNING: some rows were NOT scored in synthetic-truth mode (obs_mode != synth); their "
              "reference metrics are invalid.")
    flagged = d[d["flags"] != ""]
    if len(flagged):
        print("\n  %d truth(s) carry a flag: %s" % (len(flagged), "; ".join(
            "Ks %g f %g cc %g [%s]" % (r.Ks_mult, r.f_RS_abs, r.cc_mmhr, r.flags) for r in flagged.itertuples())))


# ------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Series 110 -- generate cc-ON synthetic truths at different TRUE Ks and f "
                    "(truth-location test, part 1).")
    parser.add_argument("--ks_values", type=float, nargs="+", default=DEFAULT_KS,
                        help="True Ks_mult values (default: %s)" % " ".join("%g" % v for v in DEFAULT_KS))
    parser.add_argument("--f_values", type=float, nargs="+", default=DEFAULT_F,
                        help="True f_RS_abs values (default: %s)" % " ".join("%g" % v for v in DEFAULT_F))
    parser.add_argument("--cc_values", type=float, nargs="+", default=DEFAULT_CC,
                        help="True cc values, mm/hr (default: %s)" % " ".join("%g" % v for v in DEFAULT_CC))
    parser.add_argument("--skip_existing", action="store_true",
                        help="Skip truths already completed (resume an interrupted run)")
    parser.add_argument("--timeout", type=int, default=300,
                        help="Per-run hard timeout in seconds (default: 300)")
    parser.add_argument("--no_cleanup", action="store_true",
                        help="Keep each run's raw results directory")
    parser.add_argument("--dry_run", action="store_true",
                        help="Run the safety checks and print the plan; build and run nothing")
    parser.add_argument("--ignore_running_check", action="store_true",
                        help="Skip the check for other tRIBS scripts running (only if you are sure none are)")
    args = parser.parse_args()

    script_dir = Path.cwd()
    project_root = script_dir.parent
    calib_dir = project_root / "calibration_work"
    summary_dir = calib_dir / "03_comparisons" / "summary_tables"
    synth_dir = calib_dir / "synth_truth"
    cand_dir = calib_dir / CAND_DIRNAME

    if TAG_PREFIX in RESERVED_TAGS:
        raise RuntimeError("TAG_PREFIX '%s' collides with an existing series tag" % TAG_PREFIX)

    # ------------------------------------------------------------------
    # Grid (checked first so a typo costs nothing)
    # ------------------------------------------------------------------
    try:
        ks_vals = validate_axis("Ks", args.ks_values, "Ks_mult")
        f_vals = validate_axis("f", args.f_values, "f_RS_abs")
        cc_vals = validate_axis("cc", args.cc_values, "cc_mmhr")
    except ValueError as e:
        parser.error(str(e))          # clean one-line message, exit code 2
    truths = make_truths(ks_vals, f_vals, cc_vals)
    run_ids = [t["run_id"] for t in truths]
    if len(set(run_ids)) != len(run_ids):
        raise RuntimeError("two truths map to the same run_id (labels collide); use more widely spaced values")

    # ------------------------------------------------------------------
    # HARD SAFETY CHECKS
    # ------------------------------------------------------------------
    qout_files = list(synth_dir.glob("*.qout")) if synth_dir.exists() else []
    if len(qout_files) != 1:
        raise RuntimeError(
            "Expected exactly one *.qout file in %s (the cc-OFF truth, used here only for reference "
            "metrics), found %d: %s. Without it the scorer would silently fall back to real-gauge mode."
            % (synth_dir, len(qout_files), [f.name for f in qout_files]))
    ref_truth = qout_files[0]
    ref_truth_md5 = md5_of(ref_truth)
    print("Reference (cc-OFF) truth, read-only: %s  [md5 %s]" % (ref_truth.name, ref_truth_md5[:12]))

    cand_res, synth_res = cand_dir.resolve(), synth_dir.resolve()
    if cand_res == synth_res or synth_res in cand_res.parents:
        raise RuntimeError("truth folder %s is, or is inside, %s; refusing." % (cand_dir, synth_dir))

    if args.ignore_running_check:
        print("Running-jobs check skipped (--ignore_running_check).")
    else:
        competing = find_competing_jobs()
        if competing is None:
            print("WARNING: could not check for other running tRIBS scripts (ps unavailable). Make sure "
                  "none are running -- they share current_run_config.json.")
        elif competing:
            print("\nANOTHER tRIBS BUILD/RUN SCRIPT APPEARS TO BE RUNNING:")
            for pid, cmd in competing:
                print("  pid %d: %s" % (pid, cmd[:110]))
            raise RuntimeError(
                "Refusing to start: these scripts share current_run_config.json with no locking, so "
                "running both can corrupt both. Wait for them to finish (or stop them with "
                "pkill -f <script>), or re-run with --ignore_running_check if you are sure they are stale.")
        else:
            print("No other tRIBS build/run script detected.")

    # ------------------------------------------------------------------
    # Plan
    # ------------------------------------------------------------------
    n_exist = sum(1 for t in truths if (cand_dir / ("%s_Outlet.qout" % t["run_id"])).exists())
    print("\n" + "=" * 92)
    print("Series 110 -- truth-location test: %d truth(s)   (Ks %s  x  f %s  x  cc %s)"
          % (len(truths), "/".join("%g" % v for v in ks_vals), "/".join("%g" % v for v in f_vals),
             "/".join("%g" % v for v in cc_vals)))
    print("  routing held at truth: cv=%g  r=%g  n=%g   optpercolation=1"
          % (ROUTING_TRUTH["kinemvelcoef"], ROUTING_TRUTH["flowexp"], ROUTING_TRUTH["channelroughness"]))
    print("  tag=%s_K<Ks>_f<f>  timeout=%ds  cleanup=%s" % (TAG_PREFIX, args.timeout,
                                                          "off" if args.no_cleanup else "on (after verified copy)"))
    print("  truth folder: %s" % cand_dir)
    print("  about %d minutes at ~1 minute per tRIBS run (centre first, then outward)" % len(truths))
    print("=" * 92)
    print("  %4s %6s %8s %8s   %-5s %-5s %-5s  run_id" % ("#", "Ks", "f", "cc", "posKs", "posf", "poscc"))
    for i, t in enumerate(truths, start=1):
        pk, pf, pc = t["pos"]
        centre = np.isclose(t["Ks_mult"], STANDARD_KS) and np.isclose(t["f_RS_abs"], STANDARD_F)
        print("  %4d %6g %8g %8g   %5.2f %5.2f %5.2f  %s%s"
              % (i, t["Ks_mult"], t["f_RS_abs"], t["cc_mmhr"], pk, pf, pc, t["run_id"],
                 "   (original centre: also a reproducibility check)" if centre else ""))
    print("  (pos = position in the sampled range: 0 = low end, 1 = high end)")
    if n_exist:
        if args.skip_existing:
            print("  %d truth file(s) already exist; --skip_existing will keep those that are complete." % n_exist)
        else:
            print("  NOTE: %d truth file(s) already exist and will be RE-RUN (use --skip_existing to resume "
                  "instead)." % n_exist)

    if args.dry_run:
        print("\n--dry_run: nothing built or run.")
        return

    summary_dir.mkdir(parents=True, exist_ok=True)
    cand_dir.mkdir(parents=True, exist_ok=True)

    out_path = summary_dir / OUT_NAME
    failed_log_path = summary_dir / FAILED_NAME
    category = builder.get_run_category(builder.PARAM_CONFIG[CC_PARAM]["series"])

    prov_base = {
        "script":           Path(__file__).name,
        "created_local":    datetime.now().isoformat(timespec="seconds"),
        "what_this_is":     "Synthetic cc-ON truths (outlet .qout) at different TRUE Ks_mult and f_RS_abs, "
                            "one per (Ks, f, cc) combination (truth-location test).",
        "routing_truth":    ROUTING_TRUTH,
        "standard_Ks_f":    {"Ks_mult": STANDARD_KS, "f_RS_abs": STANDARD_F},
        "optpercolation":   1,
        "optsnow":          0,
        "cc_units":         "mm/hr (assigned directly to CHANNELCONDUCTIVITY)",
        "run_tag_pattern":  "%s_K<Ks>_f<f>" % TAG_PREFIX,
        "grid":             {"Ks_mult": ks_vals, "f_RS_abs": f_vals, "cc_mmhr": cc_vals},
        "sampled_box":      {k: list(v[:2]) for k, v in BOX.items()},
        "reference_truth_ccoff": {"file": ref_truth.name, "md5": ref_truth_md5},
        "anchor_status":    "TEST TRUTHS ONLY -- not calibration anchors",
    }

    results = []
    if out_path.exists():
        try:
            results = pd.read_csv(out_path).to_dict("records")
            print("\n  Loaded existing results: %d rows from %s" % (len(results), out_path.name))
        except Exception as e:
            print("\n  Warning: could not load existing results (%s). Starting fresh." % e)
    done = {r.get("run_id"): r for r in results}

    failed_rows = []
    completed = skipped = hung = failed = cleaned = 0
    sweep_start = time.time()

    for i, t in enumerate(truths, start=1):
        pred_id = t["run_id"]
        cand_path = cand_dir / ("%s_Outlet.qout" % pred_id)
        print("\n[%3d/%d]  Ks=%g  f=%g  cc=%g mm/hr  -> %s"
              % (i, len(truths), t["Ks_mult"], t["f_RS_abs"], t["cc_mmhr"], pred_id), flush=True)

        if args.skip_existing and pred_id in done and cand_path.exists():
            rec_md5 = str(done[pred_id].get("candidate_md5", ""))
            if not _truthy(done[pred_id].get("qout_ok", False)):
                print("  earlier file failed its sanity check; re-running")
            elif rec_md5 and md5_of(cand_path) == rec_md5:
                print("  SKIP (already completed; file present and checksum matches)")
                skipped += 1
                continue
            else:
                print("  file present but its checksum does not match the record; re-running")

        try:
            run_id, _, log_file = builder.build_input_file(
                CC_PARAM, t["cc_mmhr"],
                overrides={**ROUTING_TRUTH,
                           "Ks_mult": t["Ks_mult"], "f_RS_abs": t["f_RS_abs"],
                           "optpercolation": 1},
                tag=t["tag"],
            )
        except Exception as e:
            print("  BUILD FAILED: %s" % e)
            failed += 1
            failed_rows.append({"run_id": pred_id, "status": "BUILD_FAILED", "reason": str(e)[:300],
                                "elapsed_min": 0.0, "Ks_mult": t["Ks_mult"], "f_RS_abs": t["f_RS_abs"],
                                "cc_mmhr": t["cc_mmhr"]})
            pd.DataFrame(failed_rows).to_csv(failed_log_path, index=False)
            continue
        if run_id != pred_id:
            raise RuntimeError("builder produced run_id %s, expected %s; refusing to continue." % (run_id, pred_id))

        returncode, elapsed, timed_out, tribs_warning = run_with_timeout(args.timeout)
        if timed_out:
            status, reason = "HANG", "wall-clock timeout"
        elif returncode != 0 or tribs_warning:
            status = "FAILED"
            reason = "tRIBS reported non-zero exit" if tribs_warning else "scorer exited %s" % returncode
        else:
            status, reason = "SUCCESS", ""

        raw_dir = calib_dir / "02_results" / category / run_id

        if status == "SUCCESS":
            metrics_file = summary_dir / ("%s_metrics_summary.csv" % run_id)
            raw_qout = locate_raw_qout(raw_dir, run_id)
            if not metrics_file.exists():
                status, reason = "FAILED", "no metrics file %s" % metrics_file.name
            elif raw_qout is None:
                status, reason = "NO_QOUT", "no single *.qout found in %s" % raw_dir

        if status == "SUCCESS":
            try:
                md5 = copy_verified(raw_qout, cand_path)
            except Exception as e:
                status, reason = "COPY_FAILED", str(e)[:300]

        if status != "SUCCESS":
            print("  %s: %s  (%s)" % (status, run_id, reason))
            if status == "HANG":
                hung += 1
            else:
                failed += 1
            failed_rows.append({"run_id": run_id, "status": status, "reason": reason,
                                "elapsed_min": elapsed / 60, "Ks_mult": t["Ks_mult"],
                                "f_RS_abs": t["f_RS_abs"], "cc_mmhr": t["cc_mmhr"]})
            pd.DataFrame(failed_rows).to_csv(failed_log_path, index=False)
            continue

        metrics = pd.read_csv(metrics_file).iloc[0].to_dict()
        qinfo = describe_qout(cand_path)
        row = assemble_row(t, run_id, metrics, cand_path, md5, qinfo, ref_truth.name)
        results = [r for r in results if r.get("run_id") != run_id]
        results.append(row)
        done[run_id] = row
        completed += 1

        print("  vs cc-OFF truth: PBIAS=%+.2f%%  KGE_2012=%.3f  volume kept=%.1f%%  peak kept=%.1f%%%s"
              % (row["pbias_pct"], row["kge_2012"], 100 * row["vol_ratio_vs_ccoff"],
                 100 * row["peak_ratio_vs_ccoff"], ("  [%s]" % row["flags"]) if row["flags"] else ""))
        print("  saved truth: %s  (%s bytes, md5 %s, %d rows)"
              % (cand_path.name, format(row["candidate_bytes"], ","), md5[:12], qinfo["qout_rows"]))
        if not qinfo["qout_ok"]:
            print("  WARNING: truth failed sanity check: %s" % qinfo["qout_problem"])

        # Delete raw output ONLY after the copy is verified AND looks sane.
        if not args.no_cleanup and qinfo["qout_ok"]:
            if cleanup_raw_results(raw_dir):
                cleaned += 1

        elapsed_total = time.time() - sweep_start
        remaining = len(truths) - i
        if completed > 0:
            eta = elapsed_total / completed * remaining / 60
            print("  Run time: %.1f min  |  ETA: %.0f min remaining" % (elapsed / 60, eta), flush=True)

        # Save after every run -- must survive an interruption at any point.
        pd.DataFrame(results).sort_values(["Ks_mult", "f_RS_abs", "cc_mmhr"]).to_csv(out_path, index=False)
        write_provenance(cand_dir, prov_base, results)

    print("\nDone: %d ran, %d skipped, %d hung, %d failed, %d raw dirs cleaned"
          % (completed, skipped, hung, failed, cleaned))

    if results:
        final_df = pd.DataFrame(results).sort_values(["Ks_mult", "f_RS_abs", "cc_mmhr"]).reset_index(drop=True)
        final_df.to_csv(out_path, index=False)
        write_provenance(cand_dir, prov_base, results)
        print("Saved: %s  (%d rows)" % (out_path, len(final_df)))
        print("Saved: %s" % (cand_dir / PROV_NAME))
        try:
            report(final_df, ref_truth.name, calib_dir)
        except Exception as e:
            print("\n(report skipped: %s -- results are saved regardless)" % e)
    else:
        print("No results to save.")
    missing = [t for t in truths
               if t["run_id"] not in done or not _truthy(done[t["run_id"]].get("qout_ok", False))]
    if missing:
        print("\n%d of the %d truth(s) in this grid are NOT complete (failed, hung, or a file that failed its "
              "sanity check):" % (len(missing), len(truths)))
        for t in missing[:10]:
            print("    %s" % t["run_id"])
        if failed or hung:
            print("See %s for the reasons." % failed_log_path)
        print("Re-run the same command with --skip_existing to retry just those.")
    else:
        print("\nAll %d truth(s) in this grid are complete." % len(truths))
        if failed_log_path.exists():
            print("(%s is a history of earlier failures that have since been re-run; it can be ignored.)"
                  % failed_log_path.name)
        print("Next: re-score the 250 + 250 saved runs against these truths (rescore_truth_location_110.py).")


if __name__ == "__main__":
    main()
