"""
run_gate3_gauge_110.py   (version 1)
=====================================
Series 110 -- the ONE tRIBS run that "gate 3" of rescore_real_gauge_110.py needs.

WHY THIS EXISTS
---------------
rescore_real_gauge_110.py (Stage 0) re-scores your 500 stored Aug 12 2014 runs against the REAL gauge
using its own copy of the scoring arithmetic. Gate 3 asks one question: does the project's ACTUAL
scorer (run_sensitivity_single_interp_tribs6.py), running in real-gauge mode, give the same numbers?
To answer it, one stored parameter point is built and run again through tRIBS and scored by the
actual scorer against the real gauge. Stage 0 then compares that run with its own numbers
(3a: the real "Observed" series, 3b: the re-run reproduces the stored simulation, 3c: the metrics).
It also vouches for how the workbook's time stamps are read, which both of Stage 0's arms share.
Without it every Stage 0 result stays labelled PROVISIONAL.

WHAT THIS SCRIPT DOES, IN ORDER
--------------------------------
  1. Picks the stored cc-ON run to repeat (rule below), and refuses the documented tRIBS hang zone.
  2. Reads your scorer, and makes a COPY in which the one line that chooses synthetic-truth mode
     is replaced by "synth_file = None" (real-gauge mode). It stops if that line is not found
     exactly once. Your scorer file is never changed.
  3. Builds the run with the project's own builder (build_sensitivity_run_tribs6.py), with the same
     settings Stage 2 used: Ks_mult, f_RS_abs and channelconductivity_mmhr from the stored point;
     kinemvelcoef 4.5, flowexp 0.24, channelroughness 0.026; optpercolation 1; tag "gate3".
  4. Checks that the builder's run configuration holds exactly those values, then points the two
     output folders of that configuration at its own folder (so nothing is added to csv_exports/
     or summary_tables/).
  5. Runs the scorer copy in its own process group with a hard time limit (a hang is killed).
  6. Checks the outputs (real-gauge mode, parameters echoed back, the re-run's simulation against
     the stored one) and prints the exact Stage 0 commands to run next.

WHY A COPY OF THE SCORER (and not simply "move the truth file")
----------------------------------------------------------------
The scorer switches to synthetic-truth mode whenever calibration_work/synth_truth/ holds exactly one
.qout file. That file (ks7_f012_Outlet.qout) must NEVER be moved out of synth_truth/. So the choice is
overridden in a copy instead. The script records every file in synth_truth/ (name and checksum)
before and after, and reports a failure if anything differs.

WHICH STORED POINT
-------------------
The cc-ON run whose (Ks_mult, log f_RS_abs, log channelconductivity_mmhr) is nearest the middle of
the sweep box (Ks 4-10.5 linear; f 0.004-0.030 log; cc 30-1000 log), among runs that are not in the
tRIBS hang zone (Ks 6.1-6.4 with f 0.0095-0.0115, a little wider than the documented 6.2-6.26 and
0.010-0.011) and not flagged as caution-zone. The rule uses only the parameters, never a score, so
the choice cannot be influenced by how well a run fits. Use --run_id to pick another stored run.

WHAT YOU GET   (calibration_work/03_comparisons/gate3_real_gauge_110/)
-----------------------------------------------------------------------
  <run_id>_compare_obs_sim.csv     the scorer's real-gauge comparison (Observed, Simulated)
  <run_id>_metrics_summary.csv     the scorer's metrics row (obs_mode must read "gauge")
  scorer_gauge_forced_110.py       the scorer copy that was run (diff of one line)
  gate3_scorer_stdout_110.txt      everything the scorer printed
  PROVENANCE_gate3_gauge_110.json  checksums, versions, the chosen point, checks, timings

SAFETY
-------
  - Never edits an existing file of yours, except calibration_work/current_run_config.json, which
    the builder itself rewrites for every run (this script changes only its two output-folder fields).
  - Writes only inside its own output folder, plus what the builder and tRIBS normally write for
    one run (the .in file, soil table, log, and the raw results under 02_results/110_channel_
    conductivity/<run_id>/; add --cleanup to delete that raw folder after a successful run).
  - Refuses an output folder that is, or lies inside, synth_truth/, csv_exports/ or summary_tables/.
  - Never starts anything without a time limit (default 600 s: tRIBS plus reading the 1993-2025
    workbook; the sweeps used 300 s for tRIBS alone).
  - --dry_run builds and runs nothing.
  - Do NOT run it while a sweep is running: the builder, this script and the scorer all share
    calibration_work/current_run_config.json.

USAGE (run from the lab/ directory, with the builder and the scorer in the same folder)
----------------------------------------------------------------------------------------
    python run_gate3_gauge_110.py --dry_run          # shows the chosen point and the patch; runs nothing
    python run_gate3_gauge_110.py                    # builds, runs, scores, checks
    python run_gate3_gauge_110.py --skip_existing    # outputs already there: re-check them, do not re-run
    python run_gate3_gauge_110.py --run_id <stored cc-ON run id> --timeout 900 --cleanup

THEN (the script prints these with the right paths)
-----------------------------------------------------
    python rescore_real_gauge_110.py --check_only --gate3_compare ... --gate3_metrics ...
    python rescore_real_gauge_110.py --margin 0.05 --primary_arm filled --gate3_compare ... --gate3_metrics ...

CAVEATS
--------
  * Stage 0's gate 3b is strict: the re-run must reproduce the stored simulation to 1e-6 of its
    peak. If it does not, today's build or inputs are not the ones that made the 500 runs. That is
    a finding about the runs, not necessarily a fault in either script.
  * One run, one point. It confirms the real-gauge path, not every run.
  * The scorer prints that one run's fit against the real gauge (KGE, PBIAS and so on). It is a single
    point, not an ON-versus-OFF comparison, and the decision margin (0.05) and the primary arm
    (filled) were fixed before this script was written. Nothing here uses that number.
  * Start with --dry_run: it checks your scorer for the one line to be replaced before anything is built.
  * Use the same environment as the sweeps (python 3.11.16, pandas 3.0.5, numpy 2.4.6).
"""

import argparse
import hashlib
import importlib
import json
import math
import os
import platform
import shutil
import signal
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

VERSION = 1
SERIES = "110"
TAG = "gate3"
PARAM_NAME = "channelconductivity_mmhr"
SCORER_DEFAULT = "run_sensitivity_single_interp_tribs6.py"
BUILDER_MODULE = "build_sensitivity_run_tribs6"
STAGE2_CSV = "lhs_results_joint_Ks_f_cc_110.csv"
OUT_FOLDER = "gate3_real_gauge_110"
PATCHED_NAME = "scorer_gauge_forced_110.py"
STDOUT_NAME = "gate3_scorer_stdout_110.txt"
PROV_NAME = "PROVENANCE_gate3_gauge_110.json"

# routing pinned at the synthetic truth, exactly as run_cc_stage2_joint_lhs_110.py did
ROUTING_TRUTH = {"kinemvelcoef": 4.5, "flowexp": 0.24, "channelroughness": 0.026}

# the Stage 2 sweep box (centre of this box is the default choice)
BOX_KS = (4.0, 10.5)
BOX_F = (0.004, 0.030)
BOX_CC = (30.0, 1000.0)

# documented tRIBS instability (methods-and-findings: Ks 6.2-6.26 with f 0.010-0.011), a little wider here
HANG_KS = (6.1, 6.4)
HANG_F = (0.0095, 0.0115)

# the scorer lines this script relies on
PATCH_OLD = "synth_file = _find_synth_truth_file(calib_dir, truth_file_override)"
PATCH_NEW = "synth_file = None  # GATE 3: real-gauge mode forced by run_gate3_gauge_110.py (original line removed)"
SCORER_MUST_CONTAIN = ["[GAUGE MODE]", "SMF_Observations_1993-2025.xlsx", "obs_mode"]

GATE3_SIM_REL = 1e-6    # same tolerance as Stage 0's gate 3b (fraction of the stored run's peak)


# ------------------------------------------------------------------
# small helpers
# ------------------------------------------------------------------
def say(msg=""):
    print(msg)
    sys.stdout.flush()


def md5_file(path):
    h = hashlib.md5()
    with open(str(path), "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def md5_text(text):
    return hashlib.md5(text.encode("utf-8")).hexdigest()


def listing(folder):
    """{relative name: md5} for every file under folder (empty if the folder does not exist)."""
    out = {}
    folder = Path(folder)
    if folder.is_dir():
        for p in sorted(folder.rglob("*")):
            if p.is_file():
                out[str(p.relative_to(folder))] = md5_file(p)
    return out


def inside(child, parent):
    c = Path(child).resolve()
    p = Path(parent).resolve()
    return c == p or p in c.parents


def jsonable(x):
    if isinstance(x, dict):
        return {str(k): jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jsonable(v) for v in x]
    if isinstance(x, (np.floating,)):
        return None if not np.isfinite(x) else float(x)
    if isinstance(x, float):
        return None if not math.isfinite(x) else x
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, Path):
        return str(x)
    return x


class Checks(object):
    def __init__(self):
        self.rows = []

    def add(self, status, name, detail):
        self.rows.append({"status": status, "name": name, "detail": detail})
        say("  %-5s %-34s %s" % (status, name, detail))

    def n(self, status):
        return sum(1 for r in self.rows if r["status"] == status)


def in_hang_zone(ks, f):
    return HANG_KS[0] <= ks <= HANG_KS[1] and HANG_F[0] <= f <= HANG_F[1]


def box_distance(ks, f, cc):
    k = (ks - BOX_KS[0]) / (BOX_KS[1] - BOX_KS[0])
    g = (math.log10(f) - math.log10(BOX_F[0])) / (math.log10(BOX_F[1]) - math.log10(BOX_F[0]))
    c = (math.log10(cc) - math.log10(BOX_CC[0])) / (math.log10(BOX_CC[1]) - math.log10(BOX_CC[0]))
    return math.sqrt((k - 0.5) ** 2 + (g - 0.5) ** 2 + (c - 0.5) ** 2)


# ------------------------------------------------------------------
# choosing the stored point
# ------------------------------------------------------------------
def load_stage2(path):
    df = pd.read_csv(str(path), float_precision="round_trip")
    need = ["run_id", "Ks_mult", "f_RS_abs", "channelconductivity_mmhr", "optpercolation"]
    miss = [c for c in need if c not in df.columns]
    if miss:
        raise SystemExit("ERROR: %s lacks column(s) %s." % (path.name, ", ".join(miss)))
    if df["run_id"].duplicated().any():
        raise SystemExit("ERROR: %s has repeated run_id values; this script will not guess which one to use." % path.name)
    return df


def candidate_table(df):
    """Every stored cc-ON run with finite parameters, with its distance from the box centre and why it may be excluded."""
    rows = []
    caution = None
    if "caution_zone" in df.columns:
        caution = df["caution_zone"].astype(str).str.strip().str.lower().isin(["true", "1"])
    for i, r in df.reset_index(drop=True).iterrows():
        ks, f, cc = float(r["Ks_mult"]), float(r["f_RS_abs"]), float(r["channelconductivity_mmhr"])
        why = []
        finite_ok = bool(np.isfinite(ks) and np.isfinite(f) and np.isfinite(cc) and f > 0 and cc > 0)
        if not finite_ok:
            why.append("non-finite or non-positive parameter")
        try:
            on = int(r["optpercolation"]) == 1
        except (TypeError, ValueError):
            on = False
        if not on:
            why.append("not a cc-ON run (optpercolation is not 1)")
        if finite_ok and in_hang_zone(ks, f):
            why.append("in the tRIBS hang zone")
        if caution is not None and bool(caution.iloc[i]):
            why.append("flagged caution-zone")
        d = box_distance(ks, f, cc) if finite_ok else float("nan")
        rows.append({"run_id": r["run_id"], "Ks_mult": ks, "f_RS_abs": f, "channelconductivity_mmhr": cc,
                     "distance": d, "excluded": "; ".join(why)})
    return pd.DataFrame(rows)


def choose_point(df, run_id):
    cand = candidate_table(df)
    if run_id is not None:
        hit = cand[cand["run_id"] == run_id]
        if hit.empty:
            raise SystemExit("ERROR: --run_id %s is not among the stored cc-ON runs in the Stage 2 results file." % run_id)
        row = hit.iloc[0]
        if row["excluded"]:
            raise SystemExit("ERROR: --run_id %s cannot be used: %s." % (run_id, row["excluded"]))
        return row, cand
    ok = cand[cand["excluded"] == ""].copy()
    if ok.empty:
        raise SystemExit("ERROR: no usable stored cc-ON run (every run is excluded; see the table above).")
    ok = ok.sort_values("distance", kind="mergesort")
    return ok.iloc[0], cand


def unique_match(df, ks, f, cc):
    """How many cc-ON stored rows match the point the way Stage 0's gate 3 will match it (rtol 1e-4)."""
    A = df[["Ks_mult", "f_RS_abs", "channelconductivity_mmhr"]].to_numpy(float)
    t = np.array([ks, f, cc])
    return int(np.sum(np.all(np.isclose(A, t[None, :], rtol=1e-4, atol=0.0), axis=1)))


# ------------------------------------------------------------------
# the scorer copy
# ------------------------------------------------------------------
def make_patched_scorer(scorer_path):
    """Return (original_text, patched_text). Stops unless the one line that picks the mode is found exactly once."""
    text = Path(scorer_path).read_text(encoding="utf-8")
    n = text.count(PATCH_OLD)
    if n != 1:
        raise SystemExit("ERROR: the line  %s  appears %d time(s) in %s (need exactly 1). The scorer has changed; this "
                         "script will not guess. Nothing was built or run." % (PATCH_OLD, n, Path(scorer_path).name))
    miss = [s for s in SCORER_MUST_CONTAIN if s not in text]
    if miss:
        raise SystemExit("ERROR: %s no longer contains %s, so the real-gauge path may have changed. Nothing was built or "
                         "run." % (Path(scorer_path).name, ", ".join(repr(m) for m in miss)))
    header = ("# GENERATED by run_gate3_gauge_110.py from %s (md5 %s). One line differs: %s  ->  %s\n"
              % (Path(scorer_path).name, md5_text(text)[:12], PATCH_OLD, "synth_file = None"))
    return text, header + text.replace(PATCH_OLD, PATCH_NEW)


# ------------------------------------------------------------------
# running the scorer
# ------------------------------------------------------------------
def kill_group(proc):
    try:
        os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
    except (ProcessLookupError, PermissionError, OSError):
        pass


def run_scorer(scorer_path, lab, timeout_sec):
    proc = subprocess.Popen([sys.executable, str(scorer_path)], cwd=str(lab), start_new_session=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    t0 = time.time()
    timed_out = False
    try:
        stdout, _ = proc.communicate(timeout=timeout_sec)
    except subprocess.TimeoutExpired:
        kill_group(proc)
        stdout, _ = proc.communicate()
        timed_out = True
    except KeyboardInterrupt:
        kill_group(proc)
        proc.communicate()
        raise
    elapsed = time.time() - t0
    stdout = stdout or ""
    tribs_warning = "WARNING: tRIBS may have failed" in stdout
    return (None if timed_out else proc.returncode), elapsed, timed_out, tribs_warning, stdout


# ------------------------------------------------------------------
# checks on what the builder and the scorer produced
# ------------------------------------------------------------------
def close(a, b, rel=1e-9):
    try:
        return math.isclose(float(a), float(b), rel_tol=rel, abs_tol=0.0)
    except (TypeError, ValueError):
        return False


def stored_window(csvd, stored_id):
    p = Path(csvd) / ("%s_compare_obs_sim.csv" % stored_id)
    if not p.exists():
        return None, None
    s = pd.read_csv(str(p), index_col=0, parse_dates=True)
    return s, p


def check_config(ck, cfg, run_id, point, stored_df):
    bad = []
    if cfg.get("run_id") != run_id:
        bad.append("run_id is %r, expected %r" % (cfg.get("run_id"), run_id))
    for key, want in (("Ks_mult", point["Ks_mult"]), ("f_RS_abs", point["f_RS_abs"]),
                      (PARAM_NAME, point["channelconductivity_mmhr"])):
        if not close(cfg.get(key), want):
            bad.append("%s is %r, expected %r" % (key, cfg.get(key), want))
    for key, want in ROUTING_TRUTH.items():
        if not close(cfg.get(key), want, 1e-12):
            bad.append("%s is %r, expected %r" % (key, cfg.get(key), want))
    try:
        if int(cfg.get("optpercolation", -1)) != 1:
            bad.append("optpercolation is %r, expected 1" % cfg.get("optpercolation"))
    except (TypeError, ValueError):
        bad.append("optpercolation is %r, expected 1" % cfg.get("optpercolation"))
    if bad:
        ck.add("FAIL", "builder configuration", "; ".join(bad))
        return False
    ck.add("PASS", "builder configuration", "run %s holds Ks %.6g, f %.6g, cc %.6g, optpercolation 1 and the truth routing"
           % (run_id, point["Ks_mult"], point["f_RS_abs"], point["channelconductivity_mmhr"]))
    if stored_df is not None and len(stored_df):
        a, b = str(stored_df.index[0]), str(stored_df.index[-1])
        es, ee = str(pd.Timestamp(cfg.get("event_start"))), str(pd.Timestamp(cfg.get("event_end")))
        if es != a or ee != b:
            ck.add("FAIL", "event window", "the run's window is %s to %s but the stored run's is %s to %s: the grids "
                   "will not match" % (es, ee, a, b))
            return False
        ck.add("PASS", "event window", "%s to %s, the same as the stored run" % (es, ee))
    return True


def check_outputs(ck, compare_csv, metrics_csv, run_id, point, stored_df, stdout):
    """Local checks only; the formal comparison is Stage 0's gate 3. Returns the preview numbers."""
    info = {}
    if "[GAUGE MODE]" in stdout and "[SYNTH MODE]" not in stdout:
        ck.add("PASS", "scorer mode (printout)", "the scorer printed [GAUGE MODE]")
    elif stdout:
        ck.add("FAIL", "scorer mode (printout)", "the scorer did not print [GAUGE MODE] (or printed [SYNTH MODE])")
    try:
        m = pd.read_csv(str(metrics_csv), float_precision="round_trip").iloc[0]
    except Exception as e:
        ck.add("FAIL", "metrics file", "cannot read %s: %s" % (Path(metrics_csv).name, str(e)[:100]))
        return info
    mode = str(m.get("obs_mode", ""))
    info["obs_mode"] = mode
    if mode == "gauge":
        ck.add("PASS", "metrics obs_mode", "obs_mode = gauge")
    else:
        ck.add("FAIL", "metrics obs_mode", "obs_mode = %r, not 'gauge': this run was scored against a synthetic truth and "
               "cannot serve as gate 3" % mode)
    bad = []
    if str(m.get("run_id")) != run_id:
        bad.append("run_id %r" % m.get("run_id"))
    for key, want in (("Ks_mult", point["Ks_mult"]), ("f_RS_abs", point["f_RS_abs"]),
                      (PARAM_NAME, point["channelconductivity_mmhr"])):
        if not close(m.get(key), want):
            bad.append("%s %r (want %r)" % (key, m.get(key), want))
    for key, want in ROUTING_TRUTH.items():
        if not close(m.get(key), want, 1e-12):
            bad.append("%s %r (want %r)" % (key, m.get(key), want))
    try:
        if int(m.get("optpercolation")) != 1:
            bad.append("optpercolation %r (want 1)" % m.get("optpercolation"))
    except (TypeError, ValueError):
        bad.append("optpercolation %r (want 1)" % m.get("optpercolation"))
    if bad:
        ck.add("FAIL", "metrics echo", "the metrics row does not echo the requested point: %s" % "; ".join(bad))
    else:
        ck.add("PASS", "metrics echo", "the metrics row echoes the requested point and routing")
    try:
        c = pd.read_csv(str(compare_csv), index_col=0, parse_dates=True)
    except Exception as e:
        ck.add("FAIL", "compare file", "cannot read %s: %s" % (Path(compare_csv).name, str(e)[:100]))
        return info
    if not {"Observed", "Simulated"}.issubset(c.columns) or c.empty or not np.all(np.isfinite(c[["Observed", "Simulated"]].to_numpy(float))):
        ck.add("FAIL", "compare file", "needs finite Observed and Simulated columns and at least one row")
        return info
    info["n_points"] = int(len(c))
    ck.add("PASS", "compare file", "%d real-gauge points; Observed peak %.2f m3/s, Simulated peak %.2f m3/s (Stage 0 checks "
           "the grid and values)" % (len(c), c["Observed"].max(), c["Simulated"].max()))
    if stored_df is None:
        ck.add("WARN", "re-run vs stored (preview)", "the stored run's compare CSV was not found, so the preview was skipped; "
               "Stage 0's gate 3b will do the check")
        return info
    sim_old = stored_df["Simulated"].reindex(c.index).to_numpy(float)
    if np.any(~np.isfinite(sim_old)):
        ck.add("FAIL", "re-run vs stored (preview)", "%d of the re-run's %d points are not on the stored run's grid"
               % (int(np.sum(~np.isfinite(sim_old))), len(c)))
        return info
    peak = float(np.max(np.abs(sim_old))) or 1.0
    d = float(np.max(np.abs(sim_old - c["Simulated"].to_numpy(float))))
    info["sim_max_abs_diff_m3s"] = d
    info["sim_rel_to_peak"] = d / peak
    if d / peak <= GATE3_SIM_REL:
        ck.add("PASS", "re-run vs stored (preview)", "the re-run reproduces the stored simulation: largest difference "
               "%.2e m3/s, %.1e of its peak (limit %.0e)" % (d, d / peak, GATE3_SIM_REL))
    else:
        ck.add("FAIL", "re-run vs stored (preview)", "the re-run differs from the stored simulation by %.2e m3/s (%.1e of "
               "its peak, limit %.0e): today's build or inputs are not the ones that made the stored runs"
               % (d, d / peak, GATE3_SIM_REL))
    return info


# ------------------------------------------------------------------
# main
# ------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Series 110: run one stored point through the ACTUAL scorer in real-gauge mode "
                                             "(the run that gate 3 of rescore_real_gauge_110.py needs).")
    ap.add_argument("--run_id", default=None, help="stored cc-ON run to repeat (default: the one nearest the middle of the "
                                                   "sweep box, outside the hang zone)")
    ap.add_argument("--timeout", type=int, default=600, help="hard time limit in seconds for tRIBS plus scoring (default 600)")
    ap.add_argument("--scorer", default=SCORER_DEFAULT, help="scorer script in this folder (default %s)" % SCORER_DEFAULT)
    ap.add_argument("--stage2_csv", default=None, help="Stage 2 results file (default summary_tables/%s)" % STAGE2_CSV)
    ap.add_argument("--out_dir", default=None, help="output folder (default calibration_work/03_comparisons/%s)" % OUT_FOLDER)
    ap.add_argument("--n_candidates", type=int, default=5, help="how many nearest candidates to list (default 5)")
    ap.add_argument("--dry_run", action="store_true", help="show the chosen point, the patch and the paths; build and run nothing")
    ap.add_argument("--skip_existing", action="store_true",
                    help="if this run's compare and metrics files already exist, re-check them and do not re-run")
    ap.add_argument("--cleanup", action="store_true", help="delete the run's raw tRIBS results folder after a successful run")
    args = ap.parse_args()

    started = datetime.now().isoformat(timespec="seconds")
    lab = Path.cwd()
    root = lab.parent
    calib = root / "calibration_work"
    summ = calib / "03_comparisons" / "summary_tables"
    csvd = calib / "03_comparisons" / "csv_exports"
    synth = calib / "synth_truth"
    workbook = root / "init_data" / "met" / "SMF_Observations_1993-2025.xlsx"
    scorer_path = lab / args.scorer
    stage2_path = Path(args.stage2_csv) if args.stage2_csv else summ / STAGE2_CSV
    out_dir = Path(args.out_dir) if args.out_dir else calib / "03_comparisons" / OUT_FOLDER

    say("=" * 100)
    say("Series 110 -- gate 3 run: one stored point through the ACTUAL scorer in real-gauge mode  (version %d)" % VERSION)
    say("=" * 100)
    say("python %s, pandas %s, numpy %s" % (platform.python_version(), pd.__version__, np.__version__))

    # ---- environment ----
    problems = []
    if not calib.is_dir():
        problems.append("no calibration_work/ folder next to the current folder (run this from lab/)")
    if not scorer_path.is_file():
        problems.append("%s not found in the current folder" % args.scorer)
    if not (lab / (BUILDER_MODULE + ".py")).is_file():
        problems.append("%s.py not found in the current folder" % BUILDER_MODULE)
    if not stage2_path.is_file():
        problems.append("Stage 2 results file not found: %s" % stage2_path)
    if not workbook.is_file():
        problems.append("real gauge workbook not found: %s" % workbook)
    for prot in (synth, csvd, summ):
        if inside(out_dir, prot):
            problems.append("the output folder %s is, or lies inside, %s" % (out_dir, prot.name))
    if problems:
        for p in problems:
            say("ERROR: %s" % p)
        sys.exit(1)

    # ---- the stored point ----
    df = load_stage2(stage2_path)
    point, cand = choose_point(df, args.run_id)
    ks, f, cc = float(point["Ks_mult"]), float(point["f_RS_abs"]), float(point["channelconductivity_mmhr"])
    n_match = unique_match(df, ks, f, cc)
    if n_match != 1:
        raise SystemExit("ERROR: the point (Ks %g, f %g, cc %g) matches %d stored cc-ON runs; Stage 0's gate 3 needs exactly 1."
                         % (ks, f, cc, n_match))
    say("synth_truth/ holds: %s   (the scorer would use it for synthetic mode; this script never touches it)"
        % (", ".join(sorted(listing(synth))) or "nothing"))
    say("")
    say("The stored cc-ON run to repeat: %s" % point["run_id"])
    say("  Ks_mult %.6g   f_RS_abs %.6g   channelconductivity_mmhr %.6g   (distance from the box centre %.3f)"
        % (ks, f, cc, float(point["distance"])))
    say("  chosen by: %s" % ("--run_id (you named it)" if args.run_id else "the nearest to the middle of the sweep box (parameters only)"))
    say("  hang zone Ks %g-%g with f %g-%g: this point is %s"
        % (HANG_KS[0], HANG_KS[1], HANG_F[0], HANG_F[1], "outside" if not in_hang_zone(ks, f) else "INSIDE"))
    ok = cand[cand["excluded"] == ""].sort_values("distance", kind="mergesort")
    n_excl = int((cand["excluded"] != "").sum())
    say("  %d of %d stored runs excluded (hang zone, caution zone or bad parameters)" % (n_excl, len(cand)))
    if args.n_candidates > 1:
        say("  nearest candidates:")
        for _, r in ok.head(args.n_candidates).iterrows():
            say("    %-34s Ks %.4g  f %.4g  cc %.4g  distance %.3f" % (r["run_id"], r["Ks_mult"], r["f_RS_abs"],
                                                                     r["channelconductivity_mmhr"], r["distance"]))
    say("  settings passed to the builder: Ks_mult, f_RS_abs, cc from the stored point; kinemvelcoef %g, flowexp %g, "
        "channelroughness %g; optpercolation 1; tag '%s'"
        % (ROUTING_TRUTH["kinemvelcoef"], ROUTING_TRUTH["flowexp"], ROUTING_TRUTH["channelroughness"], TAG))

    # ---- the scorer copy ----
    orig_text, patched_text = make_patched_scorer(scorer_path)
    say("")
    say("The scorer copy (%s is left untouched):" % args.scorer)
    say("  scorer md5 %s" % md5_text(orig_text)[:12])
    say("  one line replaced, found exactly once:")
    say("    was: %s" % PATCH_OLD)
    say("    now: %s" % PATCH_NEW)

    # ---- the builder (also gives the run id) ----
    sys.path.insert(0, str(lab))
    try:
        builder = importlib.import_module(BUILDER_MODULE)
    except Exception as e:
        raise SystemExit("ERROR: cannot import %s (%s: %s)" % (BUILDER_MODULE, type(e).__name__, str(e)[:150]))
    base_id, _ = builder.build_run_id(PARAM_NAME, cc)
    run_id = "%s_%s" % (base_id, TAG)
    compare_csv = out_dir / ("%s_compare_obs_sim.csv" % run_id)
    metrics_csv = out_dir / ("%s_metrics_summary.csv" % run_id)
    rel_cmp = os.path.relpath(str(compare_csv), str(lab))
    rel_met = os.path.relpath(str(metrics_csv), str(lab))
    say("")
    say("This run will be called %s" % run_id)
    say("  compare file : %s" % rel_cmp)
    say("  metrics file : %s" % rel_met)
    say("  time limit   : %d s" % args.timeout)

    stored_df, stored_path = stored_window(csvd, str(point["run_id"]))
    if stored_df is None:
        say("  note: the stored run's compare CSV was not found in csv_exports/; the re-run-vs-stored preview will be skipped.")

    if args.dry_run:
        say("")
        say("--dry_run: nothing was built or run, nothing was written.")
        return 0

    # ---- records before ----
    synth_before = listing(synth)
    out_dir.mkdir(parents=True, exist_ok=True)
    ck = Checks()
    status = "NOT RUN"
    build_s = run_s = None
    returncode = None
    stdout = ""
    info = {}
    cfg_path = calib / "current_run_config.json"
    skip = bool(args.skip_existing and compare_csv.exists() and metrics_csv.exists())

    say("")
    say("-" * 100)
    say("CHECKS")
    say("-" * 100)
    if skip:
        say("  --skip_existing: %s and %s exist; nothing is rebuilt or re-run." % (compare_csv.name, metrics_csv.name))
        status = "SKIPPED (existing files re-checked)"
        if (out_dir / STDOUT_NAME).exists():
            stdout = (out_dir / STDOUT_NAME).read_text(encoding="utf-8", errors="replace")
    else:
        # ---- build ----
        say("Building %s with the project's builder ..." % run_id)
        t0 = time.time()
        try:
            got_id, input_file, log_file = builder.build_input_file(
                PARAM_NAME, cc,
                overrides={"Ks_mult": ks, "f_RS_abs": f, "kinemvelcoef": ROUTING_TRUTH["kinemvelcoef"],
                           "flowexp": ROUTING_TRUTH["flowexp"], "channelroughness": ROUTING_TRUTH["channelroughness"],
                           "optpercolation": 1},
                tag=TAG)
        except Exception as e:
            say("")
            say("BUILD FAILED before tRIBS was launched: %s: %s" % (type(e).__name__, str(e)[:300]))
            status = "BUILD FAILED"
            return finish(args, started, lab, out_dir, point, run_id, orig_text, patched_text, scorer_path, synth, synth_before,
                          ck, status, None, None, None, False, info, rel_cmp, rel_met, compare_csv, metrics_csv, 1)
        build_s = time.time() - t0
        if got_id != run_id:
            ck.add("FAIL", "run id", "the builder returned run id %r but %r was expected" % (got_id, run_id))
            return finish(args, started, lab, out_dir, point, run_id, orig_text, patched_text, scorer_path, synth, synth_before,
                          ck, "BUILD CHECK FAILED", build_s, None, None, False, info, rel_cmp, rel_met, compare_csv, metrics_csv, 1)
        # ---- the builder's configuration ----
        try:
            cfg = json.loads(cfg_path.read_text())
        except Exception as e:
            ck.add("FAIL", "builder configuration", "cannot read %s: %s" % (cfg_path.name, str(e)[:100]))
            return finish(args, started, lab, out_dir, point, run_id, orig_text, patched_text, scorer_path, synth, synth_before,
                          ck, "BUILD CHECK FAILED", build_s, None, None, False, info, rel_cmp, rel_met, compare_csv, metrics_csv, 1)
        if not check_config(ck, cfg, run_id, point, stored_df):
            return finish(args, started, lab, out_dir, point, run_id, orig_text, patched_text, scorer_path, synth, synth_before,
                          ck, "BUILD CHECK FAILED", build_s, None, None, False, info, rel_cmp, rel_met, compare_csv, metrics_csv, 1)
        # point the two output folders at our own folder
        cfg["csv_export_dir"] = os.path.relpath(str(out_dir), str(lab))
        cfg["summary_export_dir"] = os.path.relpath(str(out_dir), str(lab))
        cfg_path.write_text(json.dumps(cfg, indent=2))
        ck.add("PASS", "output folders", "the scorer will write into %s" % os.path.relpath(str(out_dir), str(lab)))
        # ---- the scorer copy ----
        patched_path = out_dir / PATCHED_NAME
        patched_path.write_text(patched_text, encoding="utf-8")
        # ---- run ----
        say("")
        say("Built in %.1f s. Running tRIBS and the real-gauge scorer (limit %d s) ..." % (build_s, args.timeout))
        try:
            returncode, run_s, timed_out, tribs_warning, stdout = run_scorer(patched_path, lab, args.timeout)
        except KeyboardInterrupt:
            say("")
            say("Interrupted: the scorer's process group was killed. Nothing further was done.")
            return 130
        (out_dir / STDOUT_NAME).write_text(stdout, encoding="utf-8")
        if stdout:
            say("")
            say("----- what the scorer printed -----")
            say(stdout.rstrip("\n"))
            say("-----------------------------------")
        if timed_out:
            status = "HANG"
            ck.add("FAIL", "tRIBS and scorer", "no result within %d s: killed (a hang or a very slow run)" % args.timeout)
        elif returncode != 0 or tribs_warning:
            status = "FAILED"
            ck.add("FAIL", "tRIBS and scorer", "scorer exit code %s%s" % (returncode, "; tRIBS reported a failure" if tribs_warning else ""))
        else:
            status = "RAN"
            ck.add("PASS", "tRIBS and scorer", "finished cleanly in %.1f s (limit %d s)" % (run_s, args.timeout))
        if status != "RAN":
            say("  See the scorer printout above and the tRIBS log named in it.")
            return finish(args, started, lab, out_dir, point, run_id, orig_text, patched_text, scorer_path, synth, synth_before,
                          ck, status, build_s, run_s, returncode, timed_out, info, rel_cmp, rel_met, compare_csv, metrics_csv, 1)

    # ---- outputs ----
    if not compare_csv.exists() or not metrics_csv.exists():
        ck.add("FAIL", "outputs", "expected %s and %s in %s" % (compare_csv.name, metrics_csv.name, out_dir))
    else:
        info = check_outputs(ck, compare_csv, metrics_csv, run_id, point, stored_df, stdout)

    # ---- the truth folder ----
    synth_after = listing(synth)
    if synth_after == synth_before:
        ck.add("PASS", "synth_truth/ untouched", "%d file(s), same names and checksums before and after" % len(synth_after))
    else:
        ck.add("FAIL", "synth_truth/ untouched", "the contents of synth_truth/ changed: before %s, after %s. STOP and look at "
               "this folder before doing anything else" % (sorted(synth_before), sorted(synth_after)))

    code = 0 if ck.n("FAIL") == 0 else 1
    return finish(args, started, lab, out_dir, point, run_id, orig_text, patched_text, scorer_path, synth, synth_before,
                  ck, status if status != "RAN" else ("RAN, CHECKS PASSED" if code == 0 else "RAN, CHECKS FAILED"),
                  build_s, run_s, returncode, False, info, rel_cmp, rel_met, compare_csv, metrics_csv, code, cfg_path=cfg_path)


def finish(args, started, lab, out_dir, point, run_id, orig_text, patched_text, scorer_path, synth, synth_before, ck, status,
           build_s, run_s, returncode, timed_out, info, rel_cmp, rel_met, compare_csv, metrics_csv, code, cfg_path=None):
    """Print the verdict, write the provenance file, optionally clean up, and return the exit code."""
    synth_after = listing(synth)
    cleaned = None
    if args.cleanup and code == 0 and cfg_path is not None and Path(cfg_path).exists():
        try:
            cfg = json.loads(Path(cfg_path).read_text())
            raw = (lab / Path(cfg["output_prefix"]).parent).resolve()
            if raw.name == run_id and "02_results" in raw.parts and raw.is_dir():
                shutil.rmtree(str(raw))
                cleaned = str(raw)
        except Exception as e:
            say("  (cleanup warning: %s)" % str(e)[:120])
    next_check = ("python rescore_real_gauge_110.py --check_only --gate3_compare %s --gate3_metrics %s" % (rel_cmp, rel_met))
    next_full = ("python rescore_real_gauge_110.py --margin 0.05 --primary_arm filled --gate3_compare %s --gate3_metrics %s"
                 % (rel_cmp, rel_met))
    prov = {
        "script": "run_gate3_gauge_110.py", "version": VERSION, "started": started,
        "finished": datetime.now().isoformat(timespec="seconds"),
        "python": platform.python_version(), "pandas": pd.__version__, "numpy": np.__version__,
        "platform": platform.platform(), "status": status, "exit_code": code,
        "args": {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()},
        "stored_run": {"run_id": point["run_id"], "Ks_mult": point["Ks_mult"], "f_RS_abs": point["f_RS_abs"],
                       "channelconductivity_mmhr": point["channelconductivity_mmhr"],
                       "distance_from_box_centre": point["distance"],
                       "chosen_by": "--run_id" if args.run_id else "nearest the box centre (parameters only)"},
        "gate3_run_id": run_id,
        "settings": dict(ROUTING_TRUTH, optpercolation=1, tag=TAG),
        "hang_zone_excluded": {"Ks": list(HANG_KS), "f": list(HANG_F)},
        "scorer": {"file": str(scorer_path.name), "md5_original": md5_text(orig_text),
                   "md5_patched": md5_text(patched_text), "line_replaced": PATCH_OLD, "replacement": PATCH_NEW},
        "timing_s": {"build": build_s, "tribs_and_score": run_s}, "scorer_return_code": returncode, "timed_out": timed_out,
        "synth_truth": {"before": synth_before, "after": synth_after, "unchanged": synth_before == synth_after},
        "outputs": {"compare_csv": str(compare_csv), "metrics_csv": str(metrics_csv),
                    "compare_md5": md5_file(compare_csv) if Path(compare_csv).exists() else None,
                    "metrics_md5": md5_file(metrics_csv) if Path(metrics_csv).exists() else None},
        "preview": info, "checks": ck.rows, "raw_results_removed": cleaned,
        "next_commands": [next_check, next_full],
    }
    try:
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        (Path(out_dir) / PROV_NAME).write_text(json.dumps(jsonable(prov), indent=2))
    except Exception as e:
        say("  (could not write %s: %s)" % (PROV_NAME, str(e)[:100]))

    say("")
    say("=" * 100)
    say("RESULT: %s   (%d PASS, %d WARN, %d FAIL)" % (status, ck.n("PASS"), ck.n("WARN"), ck.n("FAIL")))
    say("=" * 100)
    if build_s is not None or run_s is not None:
        say("  build %s s, tRIBS and scoring %s s" % ("%.1f" % build_s if build_s is not None else "-",
                                                    "%.1f" % run_s if run_s is not None else "-"))
    if cleaned:
        say("  raw tRIBS results removed: %s" % cleaned)
    if code == 0:
        say("")
        say("Next (from this folder): the formal gate 3 check, which prints deviations only and shows no fit:")
        say("")
        say("    %s" % next_check)
        say("")
        say("If it says PASSED, the full Stage 0 run is:")
        say("")
        say("    %s" % next_full)
        say("")
        say("(Stage 0's gate 3 is the formal verdict; the checks above are this script's own, run first.)")
    else:
        say("")
        say("Do not run the full Stage 0 analysis yet. Send this printout (and %s) to Claude." % STDOUT_NAME)
    return code


if __name__ == "__main__":
    sys.exit(main())
