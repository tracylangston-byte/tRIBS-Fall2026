"""
run_cc_truth_candidates_110.py
==============================
Series 110 -- build a small set of CANDIDATE cc-ON synthetic truths.

WHAT THIS DOES (plain language)
--------------------------------
Every synthetic truth used so far was made with channel loss OFF
(optpercolation=0). To test recovery with channel loss ON, we need a truth
that was itself made with channel loss ON, at a known channelconductivity
(cc). We do not yet know which cc value makes a good truth, so this script
runs tRIBS once per cc value on a log-spaced grid, with every other
parameter held at the SAME truth values used throughout Series 100/101/110:

    Ks_mult = 7.0   f_RS_abs = 0.012   kinemvelcoef (cv) = 4.5
    flowexp (r) = 0.24   channelroughness (n) = 0.026
    optpercolation = 1   (channel loss ON)   optsnow = 0 (set in the builder)

and KEEPS each run's outlet hydrograph (the *_Outlet.qout file) in its own
folder. Each kept file is a candidate truth: "what the basin would produce if
the channel conductivity were exactly this value".

Nothing here is a calibration or a recovery test. It is the menu you pick a
truth from.

WHAT YOU GET
-------------
  calibration_work/synth_truth_cc_candidates/
      SMF_20140812_110_cc<value>_cctruth_Outlet.qout   (one per cc value)
      PROVENANCE_cc_truth_candidates_110.json          (what these files are)
  calibration_work/03_comparisons/summary_tables/
      cc_truth_candidates_110.csv          one row per cc value
      cc_truth_candidates_FAILED_110.csv   only if something failed

IMPORTANT: every scorer metric in the summary CSV (pbias_pct, kge_2012,
peak_error_pct, ...) is measured against the current cc-OFF truth sitting in
calibration_work/synth_truth/. That is deliberate. It tells you how much each
cc value changes the hydrograph relative to the no-loss truth -- the size of
the cc effect. They are NOT recovery scores.

WHY THE FOLDER IS SEPARATE
---------------------------
Every sweep script refuses to run unless calibration_work/synth_truth/
contains exactly ONE *.qout file, and the scorer uses that file as "observed".
Candidates must never be placed there until one is chosen on purpose. This
script only ever reads synth_truth/, and refuses to use a candidates folder
that is, or sits inside, synth_truth/.

HOW TO READ THE RESULT (rules of thumb, mine -- not from the literature)
-------------------------------------------------------------------------
  - A usable cc-ON truth sits in the INTERIOR of the Stage 2 cc box
    (30-1000 mm/hr, log scale). box_pos in the CSV is 0 at 30 and 1 at 1000;
    NEAR_BOX_EDGE flags values within 10% of either end.
  - If the cc effect is tiny (|PBIAS| < 5%, flag WEAK_EFFECT), the truth is
    hard to tell apart from the cc-OFF truth and will say little about cc.
  - If almost no water reaches the outlet (flag NEAR_DRY), the hydrograph is
    nearly empty and is a poor target.
  - cc should reduce outlet volume as it increases. The end-of-run check
    reports any reversal (PBIAS rising as cc rises).
The thresholds are constants near the top of this file; change them freely.
Whichever candidate is chosen still needs Josh's confirmation before it is
used as an anchor.

SAFETY (same rules as every other Series 110 script)
-----------------------------------------------------
  - Needs exactly one *.qout in calibration_work/synth_truth/ (for the
    reference metrics above); refuses to run otherwise.
  - run_id tag is "cctruth" (Stage 1 used s1, Stage 2 s2, control ctl,
    provenance check truthcheck), so nothing from those runs can be touched.
  - Refuses to start if another tRIBS build/run script appears to be running
    (all of them share calibration_work/current_run_config.json with no
    locking). Override with --ignore_running_check only if you are sure.
  - Each run is timeout-safe; failures are logged to the FAILED csv.
  - After each run the outlet file is copied, checksum-verified, and
    sanity-checked BEFORE the rest of that run's raw output is deleted.
    Raw output is kept if the copy or the check fails. --no_cleanup keeps all
    raw output.

USAGE (run from the lab/ directory, like every other script)
-------------------------------------------------------------
    python run_cc_truth_candidates_110.py --dry_run        # look, don't run
    python -u run_cc_truth_candidates_110.py 2>&1 | tee cc_truth_candidates_110.log
    python -u run_cc_truth_candidates_110.py --skip_existing 2>&1 | tee -a cc_truth_candidates_110.log

    python run_cc_truth_candidates_110.py --n 10 --cc_lo 50 --cc_hi 600
    python run_cc_truth_candidates_110.py --cc_values 60 100 200 300 500
    python run_cc_truth_candidates_110.py --timeout 300
    python run_cc_truth_candidates_110.py --recheck        # no tRIBS: re-evaluate
                                                           # files already saved

Default grid: 12 values, log-spaced from 45 to 700 mm/hr, rounded to 3
significant figures (45, 57.8, 74.1, 95.1, 122, 157, 201, 258, 331, 425, 545,
700). About 1 minute per run, so roughly 12-15 minutes in all.
"""

import argparse
import hashlib
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

import build_sensitivity_run_tribs6 as builder

# ------------------------------------------------------------------
# Truth values -- the SAME ones used in Stage 1, Stage 2, the control and the
# provenance check. Passed explicitly on every build (BASELINE is stale).
# ------------------------------------------------------------------
TRUTH_VALUES = {
    "Ks_mult":          7.0,
    "f_RS_abs":         0.012,
    "kinemvelcoef":     4.5,
    "flowexp":          0.24,
    "channelroughness": 0.026,
}
CC_PARAM = "channelconductivity_mmhr"

RUN_TAG       = "cctruth"
RESERVED_TAGS = {"s1", "s2", "ctl", "truthcheck", "volmatch"}

SCORER_SCRIPT = "run_sensitivity_single_interp_tribs6.py"

# Stage 2's sampled cc box (mm/hr, log scale). A truth must be interior to it
# for a Stage 2-style recovery test to be fair.
STAGE2_CC_BOX = (30.0, 1000.0)

# Default candidate grid (mm/hr). Interior to the box: the outer 10% of the
# log-width on each side is ~30-43 and ~708-1000.
DEFAULT_N     = 12
DEFAULT_CC_LO = 45.0
DEFAULT_CC_HI = 700.0

# Folder / file names.
CAND_DIRNAME = "synth_truth_cc_candidates"
OUT_NAME     = "cc_truth_candidates_110.csv"
FAILED_NAME  = "cc_truth_candidates_FAILED_110.csv"
PROV_NAME    = "PROVENANCE_cc_truth_candidates_110.json"

# Rule-of-thumb flags (see docstring). Adjust freely.
EDGE_FRAC          = 0.10   # within this fraction of the log-width of a box end
WEAK_EFFECT_PCT    = 5.0    # |PBIAS vs cc-OFF truth| below this -> WEAK_EFFECT
NEAR_DRY_VOL_RATIO = 0.10   # candidate volume / cc-OFF volume below this -> NEAR_DRY
MONO_TOL_PCT       = 0.05   # PBIAS rise (pct-pts) tolerated before calling a reversal

# Other scripts that share current_run_config.json (checked before starting).
COMPETING_PATTERNS = (
    "run_sensitivity_single", "run_sensitivity_sweep", "build_sensitivity_run",
    "run_cc_", "run_lhs_", "run_crossversion", "verify_truth_provenance",
    "validate_truth_point", "rescore_series100",
)


# ------------------------------------------------------------------
# Grid
# ------------------------------------------------------------------
def make_cc_grid(n, lo, hi):
    """Log-spaced cc values, rounded to 3 significant figures so run_ids and
    file names stay readable. Raises if rounding makes two values collide."""
    if n < 1:
        raise ValueError("--n must be at least 1")
    if not (0 < lo < hi) and not (n == 1 and 0 < lo <= hi):
        raise ValueError(f"need 0 < cc_lo < cc_hi (got {lo}, {hi})")
    raw = np.geomspace(lo, hi, n) if n > 1 else np.array([np.sqrt(lo * hi)])
    return validate_cc_values([float(f"{v:.3g}") for v in raw])


def validate_cc_values(values):
    vals = [float(v) for v in values]
    if not vals:
        raise ValueError("no cc values given")
    if any((not np.isfinite(v)) or v <= 0 for v in vals):
        raise ValueError(f"cc values must be finite and > 0: {vals}")
    if len(set(vals)) != len(vals):
        raise ValueError(
            f"{len(vals)} cc values collapse to only {len(set(vals))} distinct value(s) "
            f"(after rounding to 3 significant figures). Use fewer "
            f"values or a wider range.")
    return sorted(vals)


def box_position(cc, box=STAGE2_CC_BOX):
    """0 at the low end of the Stage 2 cc box, 1 at the high end (log scale)."""
    lo, hi = np.log10(box[0]), np.log10(box[1])
    return float((np.log10(cc) - lo) / (hi - lo))


def predicted_run_id(cc):
    """The run_id the builder will produce -- computed first so resume can
    skip a finished value without rebuilding it."""
    base, _ = builder.build_run_id(CC_PARAM, cc)
    return f"{base}_{RUN_TAG}"


# ------------------------------------------------------------------
# Small helpers
# ------------------------------------------------------------------
def md5_of(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def _ancestor_pids():
    """This process and its parents (so a wrapping shell is not mistaken for
    a competing job). Linux /proc only; returns what it can."""
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
    """Other python processes running a script that shares
    current_run_config.json. Returns a list of (pid, command line), or None
    if the check could not be done."""
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


def csv_already_exists(run_id, calib_dir):
    return (calib_dir / "03_comparisons" / "csv_exports"
            / f"{run_id}_compare_obs_sim.csv").exists()


def run_with_timeout(timeout_sec):
    """Timeout-safe subprocess execution -- same pattern as every other
    Series 110 script. Own process group so a hang is killed as a unit."""
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
    elapsed = time.time() - t0

    if stdout:
        print(stdout, end="" if stdout.endswith("\n") else "\n")

    tribs_warning = bool(stdout) and "WARNING: tRIBS may have failed" in stdout
    returncode    = None if timed_out else proc.returncode
    return returncode, elapsed, timed_out, tribs_warning


# ------------------------------------------------------------------
# The .qout file
# ------------------------------------------------------------------
EXPECTED_QOUT_COLS = 3     # Time_hr, Qstrm_m3_s, Hlev_m


def read_qout(path):
    """Parse a tRIBS outlet file by hand.

    Real file layout (seen in this Codespace, tRIBS 6.0.0): one header line
    'Time_hr,Qstrm_m3_s,Hlev_m', then COMMA-separated rows such as
    '0.0625,0.0000,0.0000'. Commas and/or whitespace are accepted, and a
    header line is skipped when its first token is not a number.

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
    """Sanity-check a .qout file. Never raises. qout_ok is True only if
    nothing looked wrong. (An earlier version of this check assumed
    whitespace-separated columns, copied from an out-of-date copy of the
    scorer, and wrongly flagged every real file.)"""
    info = {"qout_ok": False, "qout_rows": 0, "qout_ncols": 0,
            "qout_dt_min": np.nan, "qout_peak_m3s": np.nan,
            "qout_tpeak_hr": np.nan, "qout_volume_m3": np.nan,
            "qout_problem": ""}
    try:
        t, flow, ncols, bad = read_qout(path)
    except Exception as e:
        info["qout_problem"] = f"unreadable: {e}"
        return info

    problems = []
    info["qout_rows"]  = int(len(t))
    info["qout_ncols"] = int(max(ncols)) if ncols else 0
    if bad:
        problems.append(f"{bad} unparseable line(s)")
    if len(t) < 10:
        problems.append(f"only {len(t)} data rows")
    elif ncols != {EXPECTED_QOUT_COLS}:
        problems.append(f"expected {EXPECTED_QOUT_COLS} columns, saw "
                        f"{sorted(ncols)}")
    if len(t) and not (np.isfinite(t).all() and np.isfinite(flow).all()):
        problems.append("non-finite values (nan/inf)")

    if len(t) >= 10 and np.isfinite(t).all() and np.isfinite(flow).all():
        if not np.all(np.diff(t) > 0):
            problems.append("time not strictly increasing")
        if (flow < 0).any():
            problems.append("negative discharge")
        info["qout_dt_min"]    = float(np.median(np.diff(t)) * 60.0)
        info["qout_peak_m3s"]  = float(np.max(flow))
        info["qout_tpeak_hr"]  = float(t[int(np.argmax(flow))])
        info["qout_volume_m3"] = float(
            np.sum(0.5 * (flow[1:] + flow[:-1]) * np.diff(t)) * 3600.0)
        if not (info["qout_peak_m3s"] > 0):
            problems.append("no flow at all")

    info["qout_problem"] = "; ".join(problems)
    info["qout_ok"] = not problems
    return info


def locate_raw_qout(raw_dir, run_id):
    """Where tRIBS wrote this run's outlet file. Expected name first; else a
    single *.qout in the run's folder; else None."""
    expected = raw_dir / f"{run_id}_Outlet.qout"
    if expected.exists():
        return expected
    found = list(raw_dir.glob("*.qout")) if raw_dir.exists() else []
    return found[0] if len(found) == 1 else None


def copy_verified(src, dst):
    """Copy src -> dst via a temporary name, verify size and checksum, then
    move into place. Returns the md5. Raises on any mismatch (leaving no
    partial file under the final name)."""
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
    """Delete a run's raw results directory (pixel files, spatial snapshots).
    The .in file, logs and the two CSVs live elsewhere. Never raises."""
    try:
        if raw_dir.exists():
            shutil.rmtree(raw_dir)
            return True
    except Exception as e:
        print(f"  (cleanup warning: could not remove {raw_dir}: {e})")
    return False


# ------------------------------------------------------------------
# Row assembly, flags, provenance
# ------------------------------------------------------------------
def _num(x):
    try:
        v = float(x)
        return v if np.isfinite(v) else np.nan
    except (TypeError, ValueError):
        return np.nan


def make_flags(box_pos, pbias, vol_ratio, qout_ok):
    flags = []
    if not (0.0 <= box_pos <= 1.0):
        flags.append("OUTSIDE_STAGE2_BOX")
    elif box_pos < EDGE_FRAC or box_pos > 1.0 - EDGE_FRAC:
        flags.append("NEAR_BOX_EDGE")
    if np.isfinite(pbias) and abs(pbias) < WEAK_EFFECT_PCT:
        flags.append("WEAK_EFFECT")
    if np.isfinite(vol_ratio) and vol_ratio < NEAR_DRY_VOL_RATIO:
        flags.append("NEAR_DRY")
    if not qout_ok:
        flags.append("BAD_QOUT")
    return ";".join(flags)


def assemble_row(run_id, cc, metrics, cand_path, md5, qinfo, ref_truth_name):
    obs_vol  = _num(metrics.get("obs_volume_m3"))
    sim_vol  = _num(metrics.get("sim_volume_m3"))
    obs_peak = _num(metrics.get("obs_peak_m3s"))
    sim_peak = _num(metrics.get("sim_peak_m3s"))
    vol_ratio  = sim_vol / obs_vol   if obs_vol  and np.isfinite(obs_vol)  else np.nan
    peak_ratio = sim_peak / obs_peak if obs_peak and np.isfinite(obs_peak) else np.nan
    pbias = _num(metrics.get("pbias_pct"))
    bpos  = box_position(cc)

    row = {
        "run_id":                   run_id,
        "cc_mmhr":                  cc,
        "box_pos":                  bpos,
        "flags":                    make_flags(bpos, pbias, vol_ratio, qinfo["qout_ok"]),
        "reference_truth_ccoff":    ref_truth_name,
        "pbias_pct":                pbias,
        "kge_2012":                 _num(metrics.get("kge_2012")),
        "vol_ratio_vs_ccoff":       vol_ratio,
        "peak_ratio_vs_ccoff":      peak_ratio,
        "candidate_qout":           cand_path.name,
        "candidate_md5":            md5,
        "candidate_bytes":          int(cand_path.stat().st_size),
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


def recheck_existing(summary_dir, cand_dir):
    """--recheck: re-evaluate the candidate files ALREADY on disk. Runs no
    tRIBS, builds nothing, deletes nothing. Re-reads each file with the
    current sanity check, confirms its checksum still matches what was
    recorded when it was copied, recomputes qout_* columns and flags in the
    summary CSV, and refreshes the provenance JSON."""
    out_path = summary_dir / OUT_NAME
    if not out_path.exists():
        raise RuntimeError(f"nothing to recheck: {out_path} does not exist")
    df = pd.read_csv(out_path)
    if df.empty or "candidate_qout" not in df.columns:
        raise RuntimeError(f"{out_path.name} has no candidate rows to recheck")

    new = {k: [] for k in ("qout_ok", "qout_rows", "qout_ncols", "qout_dt_min",
                           "qout_peak_m3s", "qout_tpeak_hr", "qout_volume_m3",
                           "qout_problem", "flags")}
    n_changed = 0
    for _, r in df.iterrows():
        path = cand_dir / str(r["candidate_qout"])
        if path.exists():
            qinfo = describe_qout(path)
            if md5_of(path) != str(r["candidate_md5"]):
                qinfo["qout_ok"] = False
                qinfo["qout_problem"] = "; ".join(
                    x for x in (qinfo["qout_problem"],
                                "checksum differs from the one recorded at copy time")
                    if x)
        else:
            qinfo = {"qout_ok": False, "qout_rows": 0, "qout_ncols": 0,
                     "qout_dt_min": np.nan, "qout_peak_m3s": np.nan,
                     "qout_tpeak_hr": np.nan, "qout_volume_m3": np.nan,
                     "qout_problem": "candidate file missing"}
        flags = make_flags(_num(r.get("box_pos")), _num(r.get("pbias_pct")),
                           _num(r.get("vol_ratio_vs_ccoff")), qinfo["qout_ok"])
        old_flags = "" if pd.isna(r.get("flags", "")) else str(r.get("flags", ""))
        n_changed += int(flags != old_flags)
        for k in new:
            new[k].append(flags if k == "flags" else qinfo[k])
    for k, v in new.items():
        df[k] = v
    df = df.sort_values("cc_mmhr").reset_index(drop=True)
    df.to_csv(out_path, index=False)

    prov_path = cand_dir / PROV_NAME
    base = {}
    if prov_path.exists():
        try:
            base = json.loads(prov_path.read_text())
        except Exception:
            base = {}
    base.pop("files", None)
    base.pop("updated_local", None)
    base["rechecked_local"] = datetime.now().isoformat(timespec="seconds")
    write_provenance(cand_dir, base, df.to_dict("records"))

    print(f"Rechecked {len(df)} candidate file(s) in {cand_dir}")
    print(f"  flags changed on {n_changed} row(s); wrote {out_path.name} and {PROV_NAME}")
    ref = (str(df["reference_truth_ccoff"].iloc[0])
           if "reference_truth_ccoff" in df.columns else "(unknown)")
    report(df, ref)


# ------------------------------------------------------------------
# End-of-run report
# ------------------------------------------------------------------
def report(df, ref_truth_name):
    print(f"\n{'=' * 78}")
    print("CANDIDATE cc-ON TRUTHS  (metrics are vs the cc-OFF truth: "
          f"{ref_truth_name})")
    print(f"{'=' * 78}")
    d = df.sort_values("cc_mmhr").reset_index(drop=True)
    d["flags"] = d["flags"].fillna("") if "flags" in d else ""   # empty cells reload as NaN

    print(f"{'cc mm/hr':>9} {'box_pos':>8} {'PBIAS %':>9} {'KGE_2012':>9} "
          f"{'vol kept':>9} {'peak kept':>10}  flags")
    for _, r in d.iterrows():
        vk = r.get("vol_ratio_vs_ccoff", np.nan)
        pk = r.get("peak_ratio_vs_ccoff", np.nan)
        print(f"{r['cc_mmhr']:>9.4g} {r['box_pos']:>8.2f} "
              f"{_num(r.get('pbias_pct')):>+9.2f} {_num(r.get('kge_2012')):>9.3f} "
              f"{(100 * vk if np.isfinite(vk) else np.nan):>8.1f}% "
              f"{(100 * pk if np.isfinite(pk) else np.nan):>9.1f}%  "
              f"{r.get('flags', '')}")
    print("  (vol kept / peak kept = candidate relative to the cc-OFF truth, "
          "in the event window)")

    # Does more cc always remove more water?
    ok = d[np.isfinite(d["pbias_pct"].astype(float))]
    if len(ok) >= 3:
        cc = ok["cc_mmhr"].to_numpy(float)
        pb = ok["pbias_pct"].to_numpy(float)
        rho = np.corrcoef(pd.Series(cc).rank(), pd.Series(pb).rank())[0, 1]
        rev = [(cc[i], cc[i + 1], pb[i], pb[i + 1]) for i in range(len(cc) - 1)
               if pb[i + 1] > pb[i] + MONO_TOL_PCT]
        print(f"\n  Spearman rho (cc vs PBIAS) = {rho:+.3f}   "
              f"(expected strongly negative: more cc, less water)")
        if rev:
            print(f"  REVERSALS (PBIAS rose as cc rose) at {len(rev)} step(s):")
            for a, b, pa, pbb in rev:
                print(f"    cc {a:.4g} -> {b:.4g}:  PBIAS {pa:+.2f} -> {pbb:+.2f}")
        else:
            print("  No reversals: PBIAS falls (or holds) as cc rises, at every step.")

    n_bad = int((~d["qout_ok"].astype(bool)).sum()) if "qout_ok" in d else 0
    if n_bad:
        print(f"\n  WARNING: {n_bad} candidate file(s) failed the sanity check "
              f"(see qout_problem in the CSV).")
    if "obs_mode" in d and (d["obs_mode"] != "synth").any():
        print("  WARNING: some rows were NOT scored in synthetic-truth mode "
              "(obs_mode != synth); their reference metrics are invalid.")

    usable = d[(d["flags"] == "") & d["qout_ok"].astype(bool)]
    print(f"\n  {len(usable)} of {len(d)} candidates carry no flags: "
          f"{', '.join(f'{v:.4g}' for v in usable['cc_mmhr']) or 'none'}")
    print("  Final choice is a judgment call and needs Josh's confirmation "
          "before it is used as an anchor.")


# ------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Series 110 -- generate candidate cc-ON synthetic truths "
                    "at the standard truth values over a log-spaced cc grid.")
    parser.add_argument("--n", type=int, default=DEFAULT_N,
                        help=f"Number of log-spaced cc values (default: {DEFAULT_N})")
    parser.add_argument("--cc_lo", type=float, default=DEFAULT_CC_LO,
                        help=f"Low end of the cc grid, mm/hr (default: {DEFAULT_CC_LO:g})")
    parser.add_argument("--cc_hi", type=float, default=DEFAULT_CC_HI,
                        help=f"High end of the cc grid, mm/hr (default: {DEFAULT_CC_HI:g})")
    parser.add_argument("--cc_values", type=float, nargs="+", default=None,
                        help="Explicit cc values in mm/hr (overrides --n/--cc_lo/--cc_hi)")
    parser.add_argument("--skip_existing", action="store_true",
                        help="Skip cc values already completed (resume an "
                             "interrupted run)")
    parser.add_argument("--timeout", type=int, default=300,
                        help="Per-run hard timeout in seconds (default: 300)")
    parser.add_argument("--no_cleanup", action="store_true",
                        help="Keep each run's raw results directory")
    parser.add_argument("--dry_run", action="store_true",
                        help="Run the safety checks and print the plan; build "
                             "and run nothing")
    parser.add_argument("--ignore_running_check", action="store_true",
                        help="Skip the check for other tRIBS scripts running "
                             "(only if you are sure none are)")
    parser.add_argument("--recheck", action="store_true",
                        help="Do NOT run tRIBS. Re-evaluate the candidate files "
                             "already on disk (sanity check + checksum), refresh "
                             "the flags in the summary CSV and the provenance "
                             "JSON, and print the report. Deletes nothing.")
    args = parser.parse_args()

    script_dir   = Path.cwd()
    project_root = script_dir.parent
    calib_dir    = project_root / "calibration_work"
    summary_dir  = calib_dir / "03_comparisons" / "summary_tables"
    synth_dir    = calib_dir / "synth_truth"
    cand_dir     = calib_dir / CAND_DIRNAME

    if RUN_TAG in RESERVED_TAGS:
        raise RuntimeError(f"RUN_TAG '{RUN_TAG}' collides with an existing "
                           f"series tag {sorted(RESERVED_TAGS)}")

    if args.recheck:
        recheck_existing(summary_dir, cand_dir)
        return

    # ------------------------------------------------------------------
    # HARD SAFETY CHECKS
    # ------------------------------------------------------------------
    qout_files = list(synth_dir.glob("*.qout")) if synth_dir.exists() else []
    if len(qout_files) != 1:
        raise RuntimeError(
            f"Expected exactly one *.qout file in {synth_dir} (the cc-OFF "
            f"truth, used here only for reference metrics), found "
            f"{len(qout_files)}: {[f.name for f in qout_files]}. Without it "
            f"the scorer would silently fall back to real-gauge mode.")
    ref_truth      = qout_files[0]
    ref_truth_md5  = md5_of(ref_truth)
    print(f"Reference (cc-OFF) truth, read-only: {ref_truth.name}  "
          f"[md5 {ref_truth_md5[:12]}]")

    cand_res, synth_res = cand_dir.resolve(), synth_dir.resolve()
    if cand_res == synth_res or synth_res in cand_res.parents:
        raise RuntimeError(f"candidates folder {cand_dir} is, or is inside, "
                           f"{synth_dir}; refusing.")

    if args.ignore_running_check:
        print("Running-jobs check skipped (--ignore_running_check).")
    else:
        competing = find_competing_jobs()
        if competing is None:
            print("WARNING: could not check for other running tRIBS scripts "
                  "(ps unavailable). Make sure none are running -- they share "
                  "current_run_config.json.")
        elif competing:
            print("\nANOTHER tRIBS BUILD/RUN SCRIPT APPEARS TO BE RUNNING:")
            for pid, cmd in competing:
                print(f"  pid {pid}: {cmd[:110]}")
            raise RuntimeError(
                "Refusing to start: these scripts share current_run_config.json "
                "with no locking, so running both can corrupt both. Wait for "
                "them to finish (or stop them with pkill -f <script>), or "
                "re-run with --ignore_running_check if you are sure they are "
                "stale.")
        else:
            print("No other tRIBS build/run script detected.")

    # ------------------------------------------------------------------
    # Grid
    # ------------------------------------------------------------------
    try:
        if args.cc_values:
            grid = validate_cc_values(args.cc_values)
        else:
            grid = make_cc_grid(args.n, args.cc_lo, args.cc_hi)
    except ValueError as e:
        parser.error(str(e))      # clean one-line message, exit code 2
    run_ids = [predicted_run_id(cc) for cc in grid]
    if len(set(run_ids)) != len(run_ids):
        raise RuntimeError("two cc values map to the same run_id (labels "
                           "collide); use fewer or more widely spaced values")

    outside = [c for c in grid if not (STAGE2_CC_BOX[0] <= c <= STAGE2_CC_BOX[1])]
    if outside:
        print(f"\nNOTE: {len(outside)} value(s) lie outside the Stage 2 cc box "
              f"{STAGE2_CC_BOX[0]:g}-{STAGE2_CC_BOX[1]:g}: "
              f"{', '.join(f'{c:g}' for c in outside)}. A truth outside that "
              f"box could not be recovered by a Stage 2-style sweep.")

    print(f"\n{'=' * 70}")
    print(f"Series 110 -- candidate cc-ON truths  ({len(grid)} value(s))")
    print(f"  truth: Ks_mult={TRUTH_VALUES['Ks_mult']}  "
          f"f_RS_abs={TRUTH_VALUES['f_RS_abs']}  "
          f"cv={TRUTH_VALUES['kinemvelcoef']}  r={TRUTH_VALUES['flowexp']}  "
          f"n={TRUTH_VALUES['channelroughness']}  optpercolation=1")
    print(f"  tag={RUN_TAG}  timeout={args.timeout}s  "
          f"cleanup={'off' if args.no_cleanup else 'on (after verified copy)'}")
    print(f"  candidates folder: {cand_dir}")
    print(f"{'=' * 70}")
    print(f"  {'cc mm/hr':>9} {'box_pos':>8}  run_id")
    for cc, rid in zip(grid, run_ids):
        print(f"  {cc:>9.4g} {box_position(cc):>8.2f}  {rid}")

    if args.dry_run:
        print("\n--dry_run: nothing built or run.")
        return

    summary_dir.mkdir(parents=True, exist_ok=True)
    cand_dir.mkdir(parents=True, exist_ok=True)

    out_path        = summary_dir / OUT_NAME
    failed_log_path = summary_dir / FAILED_NAME
    category        = builder.get_run_category(builder.PARAM_CONFIG[CC_PARAM]["series"])

    prov_base = {
        "script":            Path(__file__).name,
        "created_local":     datetime.now().isoformat(timespec="seconds"),
        "what_this_is":      "Candidate cc-ON synthetic truths (outlet .qout) "
                             "at the standard truth values, one per cc value.",
        "truth_values":      TRUTH_VALUES,
        "optpercolation":    1,
        "optsnow":           0,
        "cc_units":          "mm/hr (assigned directly to CHANNELCONDUCTIVITY)",
        "run_tag":           RUN_TAG,
        "cc_grid_mmhr":      grid,
        "stage2_cc_box_mmhr": list(STAGE2_CC_BOX),
        "reference_truth_ccoff": {"file": ref_truth.name, "md5": ref_truth_md5},
        "anchor_status":     "CANDIDATES ONLY -- none chosen, none confirmed by Josh",
    }

    results = []
    if out_path.exists():
        try:
            results = pd.read_csv(out_path).to_dict("records")
            print(f"\n  Loaded existing results: {len(results)} rows from {out_path.name}")
        except Exception as e:
            print(f"\n  Warning: could not load existing results ({e}). Starting fresh.")
    done_ids = {r.get("run_id") for r in results}

    failed_rows = []
    completed = skipped = hung = failed = cleaned = 0
    sweep_start = time.time()

    for i, (cc, pred_id) in enumerate(zip(grid, run_ids), start=1):
        cand_path = cand_dir / f"{pred_id}_Outlet.qout"
        print(f"\n[{i:>3}/{len(grid)}]  cc={cc:g} mm/hr  (box_pos {box_position(cc):.2f})  "
              f"-> {pred_id}")

        if args.skip_existing and pred_id in done_ids and cand_path.exists():
            print("  SKIP (already completed; candidate file present)")
            skipped += 1
            continue

        try:
            run_id, _, log_file = builder.build_input_file(
                CC_PARAM, cc,
                overrides={**TRUTH_VALUES, "optpercolation": 1},
                tag=RUN_TAG,
            )
        except Exception as e:
            print(f"  BUILD FAILED: {e}")
            failed += 1
            failed_rows.append({"run_id": pred_id, "status": "BUILD_FAILED",
                                "reason": str(e)[:300], "elapsed_min": 0.0,
                                "cc_mmhr": cc})
            pd.DataFrame(failed_rows).to_csv(failed_log_path, index=False)
            continue
        if run_id != pred_id:
            raise RuntimeError(f"builder produced run_id {run_id}, expected "
                               f"{pred_id}; refusing to continue.")

        returncode, elapsed, timed_out, tribs_warning = run_with_timeout(args.timeout)
        if timed_out:
            status, reason = "HANG", "wall-clock timeout"
        elif returncode != 0 or tribs_warning:
            status = "FAILED"
            reason = ("tRIBS reported non-zero exit" if tribs_warning
                      else f"scorer exited {returncode}")
        else:
            status, reason = "SUCCESS", ""

        raw_dir = calib_dir / "02_results" / category / run_id

        if status == "SUCCESS":
            metrics_file = summary_dir / f"{run_id}_metrics_summary.csv"
            raw_qout = locate_raw_qout(raw_dir, run_id)
            if not metrics_file.exists():
                status, reason = "FAILED", f"no metrics file {metrics_file.name}"
            elif raw_qout is None:
                status, reason = "NO_QOUT", f"no single *.qout found in {raw_dir}"

        if status == "SUCCESS":
            try:
                md5 = copy_verified(raw_qout, cand_path)
            except Exception as e:
                status, reason = "COPY_FAILED", str(e)[:300]

        if status != "SUCCESS":
            print(f"  {status}: {run_id}  ({reason})")
            if status == "HANG":
                hung += 1
            else:
                failed += 1
            failed_rows.append({"run_id": run_id, "status": status,
                                "reason": reason, "elapsed_min": elapsed / 60,
                                "cc_mmhr": cc})
            pd.DataFrame(failed_rows).to_csv(failed_log_path, index=False)
            continue

        metrics = pd.read_csv(metrics_file).iloc[0].to_dict()
        qinfo   = describe_qout(cand_path)
        row     = assemble_row(run_id, cc, metrics, cand_path, md5, qinfo,
                               ref_truth.name)
        results = [r for r in results if r.get("run_id") != run_id]
        results.append(row)
        done_ids.add(run_id)
        completed += 1

        print(f"  vs cc-OFF truth: PBIAS={row['pbias_pct']:+.2f}%  "
              f"KGE_2012={row['kge_2012']:.3f}  "
              f"volume kept={100 * row['vol_ratio_vs_ccoff']:.1f}%  "
              f"peak kept={100 * row['peak_ratio_vs_ccoff']:.1f}%"
              + (f"  [{row['flags']}]" if row["flags"] else ""))
        print(f"  saved candidate: {cand_path.name}  "
              f"({row['candidate_bytes']:,} bytes, md5 {md5[:12]}, "
              f"{qinfo['qout_rows']} rows)")
        if not qinfo["qout_ok"]:
            print(f"  WARNING: candidate failed sanity check: {qinfo['qout_problem']}")

        # Delete raw output ONLY after the copy is verified AND looks sane.
        if not args.no_cleanup and qinfo["qout_ok"]:
            if cleanup_raw_results(raw_dir):
                cleaned += 1

        elapsed_total = time.time() - sweep_start
        remaining = len(grid) - completed - skipped - hung - failed
        if completed > 0:
            eta = elapsed_total / completed * remaining / 60
            print(f"  Run time: {elapsed / 60:.1f} min  |  ETA: {eta:.0f} min remaining")

        # Save after every run -- must survive an interruption at any point.
        pd.DataFrame(results).sort_values("cc_mmhr").to_csv(out_path, index=False)
        write_provenance(cand_dir, prov_base, results)

    print(f"\nDone: {completed} ran, {skipped} skipped, {hung} hung, "
          f"{failed} failed, {cleaned} raw dirs cleaned")

    if results:
        final_df = pd.DataFrame(results).sort_values("cc_mmhr").reset_index(drop=True)
        final_df.to_csv(out_path, index=False)
        write_provenance(cand_dir, prov_base, results)
        print(f"Saved: {out_path}  ({len(final_df)} rows)")
        print(f"Saved: {cand_dir / PROV_NAME}")
        try:
            report(final_df, ref_truth.name)
        except Exception as e:
            print(f"\n(report skipped: {e} -- results are saved regardless)")
    else:
        print("No results to save.")
    if failed or hung:
        print(f"\n{failed + hung} run(s) did not complete -- see {failed_log_path}")


if __name__ == "__main__":
    main()
