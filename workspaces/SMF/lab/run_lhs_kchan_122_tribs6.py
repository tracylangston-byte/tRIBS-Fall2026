"""
run_lhs_kchan_122_tribs6.py
===========================
Series 122 ("Scenario 5"), Stage 1 -- channel conductivity (Kchan) with routing free,
on tRIBS 6.0.0 / pytRIBS 1.0.0.

WHAT THIS DOES
--------------
Runs a 400-point Latin hypercube over FOUR parameters with channel loss ON
(optpercolation 1), and keeps the outlet hydrograph of every run:

    kinemvelcoef (cv)         2.5  - 6.5     linear     (the Series 121 box)
    flowexp (r)               0.18 - 0.35    linear     (the Series 121 box)
    channelroughness (n)      0.02 - 0.10    linear     (the Series 121 box)
    channelconductivity_mmhr  3    - 1000    log        (Kchan, mm/hr; NEW)

Ks and f are FIXED at the truth values (Ks 7.0, f 0.012). Seed 301.

Every run is scored by the project scorer against truth A (the Series 121 truth:
Ks 7.0, f 0.012, cv 4.5, flowexp 0.24, n 0.026, channel loss OFF), which is the single
*.qout in calibration_work/synth_truth/. Because every run's outlet file is kept, the same
runs can be scored later against truth B (the truth WITH channel loss, Kchan 201) and against
any later truth or noise test, with no new tRIBS runs. This script does not need truth B.

This script only RUNS the sweep and records it. The reading rules (recovery, false use,
"does Kchan add a direction") are applied afterwards by the analysis script, and were fixed
before any result existed (see the Series 121/122 handoff, section 5).

WHERE THINGS GO
---------------
    run id            SMF_20140812_122_Ks7p0x_d<draw>_cv<cv>_r<r>_n<n>_cc<Kchan>
    results CSV       calibration_work/03_comparisons/summary_tables/lhs_results_kchan_122.csv
    failed log        calibration_work/03_comparisons/summary_tables/lhs_results_kchan_FAILED_122.csv
    design record     calibration_work/03_comparisons/summary_tables/lhs_design_kchan_122.json
    kept outlet files calibration_work/kept_qout_122/<run id>_Outlet.qout   (one per run)
    raw tRIBS output  calibration_work/02_results/122_lhs_kchan/<run id>/    (deleted after the
                      outlet file has been copied and the copy checked by checksum)
    run inputs        calibration_work/01_run_inputs/122_lhs_kchan/
(If you use a different --n or --seed, "_n<N>_seed<S>" is added to the three summary_tables
names and to the kept folder name, so nothing from the default design is touched.)

Nothing in calibration_work/synth_truth/ is written or moved; it is only READ.

USAGE (run from the lab/ directory, like every other run script)
----------------------------------------------------------------
    cd /workspaces/tRIBS-Fall2026/workspaces/SMF/lab
    python -u run_lhs_kchan_122_tribs6.py --dry_run
    python -u run_lhs_kchan_122_tribs6.py --limit 3 2>&1 | tee series122_trial.log
    python -u run_lhs_kchan_122_tribs6.py --skip_existing 2>&1 | tee series122.log
    python -u run_lhs_kchan_122_tribs6.py --audit

  --dry_run   writes, builds and runs nothing; prints the design, the checks and the time estimate
  --limit N   only the first N draws (a trial; the full run later resumes with --skip_existing)
  --skip_existing   resume: a draw is skipped only if its compare CSV, its kept outlet file and
              its row in the results CSV all exist
  --audit     reads the results and the kept outlet files and checks them (no tRIBS)

After an interruption (Codespace stopped, tab closed, Ctrl+C) run the same sweep command again:
--skip_existing picks up where it stopped. Keep the terminal open and visible while it runs, and
do not run any other tRIBS build/run script in lab/ at the same time. The -u matters: without it
Python holds output back when it is piped to tee.

CHECKS BUILT IN (the sweep stops, with a plain message, when one fails)
-----------------------------------------------------------------------
  * exactly one *.qout in synth_truth/; it must be readable, not dry, and its last time must match
    the simulation window this script will use (so truth and sweep use the same window)
  * the builder must not divide Kchan by 3.6e6 (the Series 110 unit bug), and must have the
    keys this script overrides
  * after every build, the builder's run configuration must carry exactly the Ks, f, cv, r, n,
    Kchan and optpercolation asked for; the .in file is checked for optpercolation and Kchan where
    its layout allows
  * after every run, the scorer must have used truth A as "observed" (obs_mode synth, its own
    "Reading observed from:" line, and the observed peak equal to truth A's peak) and the scored
    row must carry the parameters asked for
  * the raw outlet file must reach the end of the window; the kept copy is verified by checksum
  * after 20 scored runs, Kchan must visibly move PBIAS (Spearman rank correlation of PBIAS with
    log10 Kchan below -0.5; Series 110 Stage 1 gave -0.994). If not, channel loss is not wired
    in and the sweep stops instead of wasting hours. (The -0.5 and the 20 are safety settings,
    not scientific reading rules.)
  * each run is a subprocess with a hard timeout (default 300 s); a hang kills its process group.
    HANG / FAILED runs go to the FAILED log and the sweep moves on; 8 in a row stop the script
  * refuses to start if another tRIBS build/run script seems to be running (they share
    calibration_work/current_run_config.json with no locking)
  * on resume, the draws must be identical to the saved design (same seed, n, box)
  * at least 2 GB of free disk before every run

STATUS
------
Written without access to the Codespace. NOT YET RUN against the real tRIBS 6.0.0 build.
Tested against stand-ins for the builder and the scorer (a toy hydrograph model, NOT tRIBS):
  passed  dry run; 3-run trial; full 12-run sweep + --audit; resume (--skip_existing); refusal to
          overwrite; changed draws; swapped truth file; two truth files; truth/window mismatch;
          .in file with channel loss off; inert Kchan (stop at run 20); scorer in real-gauge mode;
          hung run; run that reports failure; short outlet file; planted 3.6e6 unit bug; 8 failed
          builds in a row; competing job running; deleted / stray / altered kept outlet file.
  NOT tested: the real builder's .in layout, the real scorer's output, real timings, --runtime_hours
          against the real builder, and a real hang.
The first real checks are --dry_run (does the truth file window match, does the builder look right)
and a --limit 3 trial. In the trial output, look for a line starting "(note: could not read
OPTPERCOLATION" -- if it appears, the real .in file is laid out differently from what this script
expects and the .in-file check is being skipped (the inert-Kchan stop at run 20 still protects you).
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
SERIES        = "122"
CATEGORY      = "122_lhs_kchan"
SCORER_SCRIPT = "run_sensitivity_single_interp_tribs6.py"

# ----------------------------------------------------------------------
# Truth and fixed values
# ----------------------------------------------------------------------
TRUTH_KS = 7.0
TRUTH_F  = 0.012

# Truth routing values (for the record and the dry-run printout only; routing is FREE here).
TRUTH_ROUTING = {"kinemvelcoef": 4.5, "flowexp": 0.24, "channelroughness": 0.026}

# Sampling box. Order matters only for reproducing the draws (it sets how the random
# generator is consumed), so do not reorder. "dp" = decimal places the value is rounded to
# BEFORE it is run, so that the run id, the .in file and the results row all carry the same
# number.
LHS_PARAMS = {
    "kinemvelcoef":             {"lo": 2.5,  "hi": 6.5,    "scale": "linear", "dp": 4},
    "flowexp":                  {"lo": 0.18, "hi": 0.35,   "scale": "linear", "dp": 5},
    "channelroughness":         {"lo": 0.02, "hi": 0.10,   "scale": "linear", "dp": 5},
    "channelconductivity_mmhr": {"lo": 3.0,  "hi": 1000.0, "scale": "log",    "dp": 3},
}
PARAM_COLS   = list(LHS_PARAMS.keys())
DEFAULT_N    = 400
DEFAULT_SEED = 301            # Series 121 used 42 and 201 to 205

# ----------------------------------------------------------------------
# Safety settings (not scientific reading rules)
# ----------------------------------------------------------------------
DEFAULT_TIMEOUT      = 300
SEC_PER_MODEL_HOUR   = 0.152     # about 67 s per run at 450 h, 219 s at 1416 h (Series 110 timings)
QOUT_END_SLACK_HR    = 2.0       # an outlet file must reach to within this many hours of the window end
MIN_FREE_GB          = 2.0
INERT_CHECK_AT       = 20        # scored runs before the "does Kchan move PBIAS" check
INERT_RHO_MAX        = -0.5      # Spearman(PBIAS, log10 Kchan) must be below this
OBS_PEAK_RTOL        = 0.02
EXPECTED_QOUT_COLS   = 3
MAX_CONSEC_FAIL      = 8
FAILED_COLS = ["run_id", "status", "reason", "elapsed_min", "draw_idx"] + PARAM_COLS

COMPETING_PATTERNS = (
    "run_sensitivity_single", "run_sensitivity_sweep", "build_sensitivity_run",
    "run_cc_", "run_lhs_", "run_lowf_sweep", "run_crossversion", "verify_truth_provenance",
    "validate_truth_point", "rescore_series100", "run_truth_",
)


class Abort(RuntimeError):
    """A problem that would make every further run pointless: stop the whole script."""


class BuildFailed(RuntimeError):
    """The builder raised for one run: log it and carry on."""


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


def _close(a, b, rtol=1e-6):
    try:
        return math.isclose(float(a), float(b), rel_tol=rtol, abs_tol=1e-12)
    except (TypeError, ValueError):
        return False


def write_json_atomic(path, obj):
    tmp = path.with_name(path.name + ".part")
    tmp.write_text(json.dumps(obj, indent=2, default=str))
    os.replace(tmp, path)


def write_csv_atomic(df, path):
    tmp = path.with_name(path.name + ".part")
    df.to_csv(tmp, index=False)
    os.replace(tmp, path)


def label(v):
    """The builder's own number-to-text rule (7.0 -> 7p0, 0.012 -> 0p012)."""
    return builder.value_to_label(v)


def variant_suffix(n, seed):
    s = ""
    if n != DEFAULT_N:
        s += "_n%d" % n
    if seed != DEFAULT_SEED:
        s += "_seed%d" % seed
    return s


def spearman(x, y):
    x = pd.Series(np.asarray(x, dtype=float))
    y = pd.Series(np.asarray(y, dtype=float))
    ok = x.notna() & y.notna()
    if ok.sum() < 3:
        return np.nan
    return float(x[ok].rank().corr(y[ok].rank()))


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


# ======================================================================
# Sampler (Latin hypercube; log-stratified where the parameter is log-scale)
# ======================================================================
def generate_lhs_raw(n, params, seed):
    """Each parameter is split into n equal slices (log slices for log-scale parameters), one
    uniform draw per slice, then shuffled independently. Same seed + n + box + dict order
    = identical draws."""
    rng = np.random.default_rng(seed)
    raw = {}
    for name, b in params.items():
        lo, hi = b["lo"], b["hi"]
        if b["scale"] == "log":
            edges = np.linspace(np.log10(lo), np.log10(hi), n + 1)
            pts = 10 ** rng.uniform(edges[:-1], edges[1:])
        else:
            edges = np.linspace(lo, hi, n + 1)
            pts = rng.uniform(edges[:-1], edges[1:])
        rng.shuffle(pts)
        raw[name] = pts
    return pd.DataFrame(raw)


def lhs_check(raw, params, n):
    """Exactly one draw in each of the n slices, for every parameter. Returns {param: bool}."""
    out = {}
    for name, b in params.items():
        x = raw[name].to_numpy(float)
        if b["scale"] == "log":
            u = (np.log10(x) - np.log10(b["lo"])) / (np.log10(b["hi"]) - np.log10(b["lo"]))
        else:
            u = (x - b["lo"]) / (b["hi"] - b["lo"])
        idx = np.clip(np.floor(u * n).astype(int), 0, n - 1)
        out[name] = bool(len(np.unique(idx)) == n)
    return out


def round_samples(raw, params):
    s = raw.copy()
    for name, b in params.items():
        s[name] = np.clip(np.round(s[name].to_numpy(float), b["dp"]), b["lo"], b["hi"])
    return s


# ======================================================================
# Names, and teaching the (unedited) builder about series 122
# ======================================================================
def install_series_hooks():
    """In-memory only: Ks_mult runs are series 122 and series 122 lives in this category.
    The builder file on disk is not touched. Stops if the result is not what the rest of
    this script assumes."""
    cfg = dict(builder.PARAM_CONFIG["Ks_mult"])
    cfg["series"] = SERIES
    builder.PARAM_CONFIG["Ks_mult"] = cfg

    original = builder.get_run_category

    def get_run_category_122(series_str):
        if str(series_str) == SERIES:
            return CATEGORY
        return original(series_str)

    builder.get_run_category = get_run_category_122

    run_id, _ = builder.build_run_id("Ks_mult", TRUTH_KS)
    expected = "%s_%s_%s_Ks7p0x" % (builder.LOCATION, builder.EVENT_DATE, SERIES)
    if run_id != expected:
        raise RuntimeError("Series-122 hook self-check failed: build_run_id gave %r, expected %r. The builder's "
                           "naming has changed; look at build_sensitivity_run_tribs6.py before running anything."
                           % (run_id, expected))
    cat = builder.get_run_category(builder.PARAM_CONFIG["Ks_mult"]["series"])
    if cat != CATEGORY:
        raise RuntimeError("Series-122 hook self-check failed: category %r, expected %r." % (cat, CATEGORY))


def builder_unit_guard():
    """The Series 110 bug: the builder divided Kchan by 3.6e6 before writing it, so every Kchan
    collapsed to nearly zero. Stop if an executable line still does that."""
    src = inspect.getsource(builder.build_input_file)
    for ln in src.splitlines():
        code = ln.split("#", 1)[0]
        if "channelconductivity" in code and "3.6e6" in code:
            raise RuntimeError(
                "The builder still divides channel conductivity by 3.6e6 (line: %r). The .in field is already in "
                "mm/hr, so every Kchan would be written as nearly zero. Use the fixed builder." % ln.strip())
    needed = ("optpercolation", "channelconductivity_mmhr", "channelporosity",
              "kinemvelcoef", "flowexp", "channelroughness")
    missing = [k for k in needed if k not in builder.BASELINE]
    if missing:
        raise RuntimeError("The builder's BASELINE lacks %s, so this script cannot override them." % missing)
    params = inspect.signature(builder.build_input_file).parameters
    for p in ("overrides", "tag"):
        if p not in params:
            raise RuntimeError("builder.build_input_file has no %r argument. This script needs the tRIBS-6 builder." % p)


def run_tag(idx, row):
    return "d%03d_cv%s_r%s_n%s_cc%s" % (
        idx + 1, label(row["kinemvelcoef"]), label(row["flowexp"]),
        label(row["channelroughness"]), label(row["channelconductivity_mmhr"]))


def predicted_run_id(tag):
    run_id, _ = builder.build_run_id("Ks_mult", TRUTH_KS)
    return "%s_%s" % (run_id, tag)


def overrides_for(row):
    """Everything that must be pinned or set for a Series 122 run."""
    return {
        "Ks_mult":                  TRUTH_KS,
        "f_RS_abs":                 TRUTH_F,
        "kinemvelcoef":             float(row["kinemvelcoef"]),
        "flowexp":                  float(row["flowexp"]),
        "channelroughness":         float(row["channelroughness"]),
        "optpercolation":           1,
        "channelconductivity_mmhr": float(row["channelconductivity_mmhr"]),
    }


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
    info = {"qout_ok": False, "qout_rows": 0, "qout_ncols": 0, "qout_dt_min": np.nan,
            "qout_tlast_hr": np.nan, "qout_peak_m3s": np.nan, "qout_volume_m3": np.nan,
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
        info["qout_tlast_hr"] = float(t[-1])
        info["qout_peak_m3s"] = float(np.max(flow))
        info["qout_volume_m3"] = float(np.sum(0.5 * (flow[1:] + flow[:-1]) * np.diff(t)) * 3600.0)
        if not (info["qout_peak_m3s"] > 0):
            problems.append("no flow at all")
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


def cleanup_raw_results(run_id, calib_dir):
    """Delete a successful run's raw results folder. Only deletes a folder located EXACTLY at
    02_results/<category>/<run_id>. Never raises."""
    try:
        parent = (calib_dir / "02_results" / CATEGORY).resolve()
        raw = (parent / run_id).resolve()
        if raw.parent != parent or raw.name != run_id:
            print("  (cleanup skipped: %s is not directly inside %s)" % (raw, parent))
            return False
        if raw.exists():
            shutil.rmtree(raw)
            return True
    except Exception as e:
        print("  (cleanup warning for %s: %s)" % (run_id, e))
    return False


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


# ======================================================================
# Running the scorer
# ======================================================================
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
# One run: build, check the build, score, check the score, keep the outlet file
# ======================================================================
def in_file_value(text, keyword):
    """Value of a tRIBS .in keyword: the first number on the line after the line that starts with
    KEYWORD (the label line holds a description, the value sits on the next line). Returns None if
    the layout is not what this expects."""
    lines = text.splitlines()
    for i, ln in enumerate(lines):
        s = ln.strip()
        if s.upper().startswith(keyword.upper()):
            for cand in lines[i + 1:i + 3]:
                toks = cand.strip().split()
                if toks:
                    try:
                        return float(toks[0])
                    except ValueError:
                        return None
            return None
    return None


def check_config(ctx, cfg, run_id, row):
    """What the builder wrote into current_run_config.json versus what was asked for."""
    problems = []
    if cfg.get("run_id") != run_id:
        problems.append("run_id is %r, expected %r" % (cfg.get("run_id"), run_id))
    if cfg.get("run_category") not in (None, CATEGORY):
        problems.append("run_category is %r, expected %r" % (cfg.get("run_category"), CATEGORY))
    wanted = [("Ks_mult", TRUTH_KS), ("f_RS_abs", TRUTH_F), ("optpercolation", 1),
              ("kinemvelcoef", row["kinemvelcoef"]), ("flowexp", row["flowexp"]),
              ("channelroughness", row["channelroughness"]),
              ("channelconductivity_mmhr", row["channelconductivity_mmhr"])]
    for key, want in wanted:
        if key not in cfg:
            problems.append("%s missing from the run configuration" % key)
        elif not _close(cfg[key], want, 1e-9):
            problems.append("%s is %r, expected %r" % (key, cfg[key], want))
    if "runtime_hours" in cfg and not _close(cfg["runtime_hours"], ctx.runtime_hours, 1e-9):
        problems.append("runtime_hours is %r, expected %r" % (cfg["runtime_hours"], ctx.runtime_hours))
    if cfg.get("truth_file"):
        problems.append("the run configuration carries truth_file=%r, but this series is meant to be scored "
                        "against the single truth in synth_truth/" % cfg.get("truth_file"))
    inp = cfg.get("input_file")
    if not inp:
        problems.append("no input_file in the run configuration")
    else:
        p = Path(inp) if Path(inp).is_absolute() else Path.cwd() / inp
        try:
            text = p.read_text(errors="replace")
            got_opt = in_file_value(text, "OPTPERCOLATION")
            got_cc = in_file_value(text, "CHANNELCONDUCTIVITY")
            if got_opt is None or got_cc is None:
                if not ctx.in_note_printed:
                    print("  (note: could not read OPTPERCOLATION / CHANNELCONDUCTIVITY from %s in the layout I "
                          "expected, so the .in-file check was skipped; the configuration check above still applies)"
                          % p.name)
                    ctx.in_note_printed = True
            else:
                if int(round(got_opt)) != 1:
                    problems.append("the .in file has OPTPERCOLATION %r, expected 1" % got_opt)
                if not _close(got_cc, row["channelconductivity_mmhr"], 1e-3):
                    problems.append("the .in file has CHANNELCONDUCTIVITY %r, expected %r (mm/hr)"
                                    % (got_cc, row["channelconductivity_mmhr"]))
        except Exception as e:
            problems.append("cannot read the .in file %s (%s)" % (inp, e))
    return problems


def prepare_run(ctx, pred_id, tag, row):
    try:
        run_id, _, _ = builder.build_input_file("Ks_mult", TRUTH_KS, overrides=overrides_for(row), tag=tag)
    except Exception as e:
        raise BuildFailed(str(e)[:300])
    if run_id != pred_id:
        raise Abort("builder produced run_id %s, expected %s" % (run_id, pred_id))
    cfg_path = ctx.calib_dir / "current_run_config.json"
    try:
        cfg = json.loads(cfg_path.read_text())
    except Exception as e:
        raise Abort("cannot read %s after the build: %s" % (cfg_path, e))
    problems = check_config(ctx, cfg, run_id, row)
    if problems:
        raise Abort("the builder's run configuration does not match what was asked for: " + "; ".join(problems))
    return run_id, cfg


def check_observed(metrics, cfg, stdout, truth_path):
    """Did the scorer really use truth A as 'observed'? Returns a list of problems."""
    problems = []
    if str(metrics.get("obs_mode")) != "synth":
        problems.append("obs_mode is %r, expected 'synth' (the scorer fell back to the real gauge)"
                        % metrics.get("obs_mode"))
    m = re.search(r"Reading observed from:\s*(.+)", stdout or "")
    if m and m.group(1).strip() != truth_path.name:
        problems.append("the scorer read %r as observed, expected %r" % (m.group(1).strip(), truth_path.name))
    elif not m:
        print("  (note: the scorer printed no 'Reading observed from:' line, so that check was skipped; "
              "the observed-peak check below still applies)")
    try:
        exp = expected_obs_peak(truth_path, cfg["event_start"], cfg["event_end"])
        got = float(metrics["obs_peak_m3s"])
        if not (abs(got - exp) <= OBS_PEAK_RTOL * max(abs(exp), 1e-12)):
            problems.append("the observed peak the scorer used (%.4g m3/s) is not the peak of %s (%.4g m3/s)"
                            % (got, truth_path.name, exp))
    except Exception as e:
        print("  (note: could not double-check the observed peak: %s)" % str(e)[:100])
    return problems


def readback_problems(m, row):
    """Compare the scored metrics row against what this script asked for."""
    probs = []
    expected = {"Ks_mult": TRUTH_KS, "f_RS_abs": TRUTH_F, "optpercolation": 1,
                "kinemvelcoef": row["kinemvelcoef"], "flowexp": row["flowexp"],
                "channelroughness": row["channelroughness"],
                "channelconductivity_mmhr": row["channelconductivity_mmhr"]}
    for key, want in expected.items():
        got = _num(m.get(key))
        if not np.isfinite(got) or not np.isclose(got, want, rtol=1e-6, atol=1e-12):
            probs.append("%s: scored row says %r, asked for %r" % (key, m.get(key), want))
    return probs


def keep_outlet(ctx, run_id):
    """Copy a finished run's outlet file into kept_qout_122 (verified by checksum) and check it.
    Returns (status, reason, md5, qinfo). status OK or a failure code."""
    raw_dir = ctx.calib_dir / "02_results" / CATEGORY / run_id
    raw_qout = locate_raw_qout(raw_dir, run_id)
    if raw_qout is None:
        return "NO_QOUT", "no single *.qout found in %s" % raw_dir, None, None
    dest = ctx.kept_dir / ("%s_Outlet.qout" % run_id)
    try:
        md5 = copy_verified(raw_qout, dest)
    except Exception as e:
        return "COPY_FAILED", str(e)[:300], None, None
    qinfo = describe_qout(dest)
    if not qinfo["qout_ok"]:
        _set_aside_rejected(ctx, dest)
        return "QOUT_BAD", "outlet file failed its sanity check: %s" % qinfo["qout_problem"], md5, qinfo
    if qinfo["qout_tlast_hr"] < ctx.runtime_hours - QOUT_END_SLACK_HR:
        _set_aside_rejected(ctx, dest)
        return ("QOUT_SHORT", "outlet file ends at hour %.2f, but the window is %g h (a run that stopped early?)"
                % (qinfo["qout_tlast_hr"], ctx.runtime_hours), md5, qinfo)
    return "OK", "", md5, qinfo


def _set_aside_rejected(ctx, dest):
    """A copy that failed its checks must not sit among the good outlet files, where a later
    analysis (truth B scoring) could pick it up. Move it into kept_qout_122/_rejected/ (nothing is
    deleted; the original also stays in its raw results folder)."""
    try:
        rej = ctx.kept_dir / "_rejected"
        rej.mkdir(exist_ok=True)
        os.replace(dest, rej / dest.name)
    except Exception as e:
        print("  (could not move the rejected outlet file aside: %s -- delete %s by hand before analysing)" % (e, dest.name))


def check_disk(ctx):
    free_gb = shutil.disk_usage(ctx.calib_dir).free / 1e9
    if free_gb < MIN_FREE_GB:
        raise Abort("only %.1f GB of disk is free (need at least %.1f GB for one run). Free some space and re-run "
                    "the same command with --skip_existing." % (free_gb, MIN_FREE_GB))


def run_one(ctx, pred_id, tag, row):
    """Build, score, verify, and keep the outlet file of one run. Returns a dict with status in
    SUCCESS / BUILD_FAILED / HANG / FAILED / NO_QOUT / COPY_FAILED / QOUT_BAD / QOUT_SHORT /
    OBSERVED_MISMATCH / READBACK_MISMATCH, plus the metrics row. Raises Abort for a problem that
    would affect every run."""
    out = {"run_id": pred_id, "status": None, "reason": "", "elapsed": 0.0, "metrics": None, "qinfo": None,
           "qout_md5": None}
    check_disk(ctx)
    try:
        run_id, cfg = prepare_run(ctx, pred_id, tag, row)
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

    problems = check_observed(metrics, cfg, stdout, ctx.truth_path)
    if problems:
        out["abort"] = ("scoring did not use truth A as 'observed' for %s: %s. Every later run would be scored the "
                        "same wrong way, so the script is stopping." % (run_id, "; ".join(problems)))
        out.update(status="OBSERVED_MISMATCH", reason="; ".join(problems))
        return out
    probs = readback_problems(metrics, row)
    if probs:
        out["abort"] = ("the scored row for %s does not match what this script asked for: %s. The override plumbing "
                        "is broken, so the remaining runs would be wasted." % (run_id, "; ".join(probs)))
        out.update(status="READBACK_MISMATCH", reason="; ".join(probs))
        return out

    st, why, md5, qinfo = keep_outlet(ctx, run_id)
    out["qinfo"], out["qout_md5"] = qinfo, md5
    if st != "OK":
        out.update(status=st, reason=why)           # raw folder is left in place for a look
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
        raise Abort("%d runs in a row failed. That looks like a systematic problem, not bad luck. Look at the last "
                    "error above, fix it, and re-run the same command with --skip_existing." % ctx.consec_fail)


def inert_check(ctx, results):
    """After INERT_CHECK_AT scored runs, Kchan must visibly move PBIAS. If not, channel loss is not
    wired in. Raises Abort. Returns the rank correlation (or None if not enough runs yet)."""
    df = pd.DataFrame(results)
    if len(df) < INERT_CHECK_AT or "pbias_pct" not in df.columns:
        return None
    rho = spearman(np.log10(pd.to_numeric(df["channelconductivity_mmhr"], errors="coerce")),
                   pd.to_numeric(df["pbias_pct"], errors="coerce"))
    if not np.isfinite(rho) or rho > INERT_RHO_MAX:
        raise Abort(
            "KCHAN LOOKS INERT. Over %d scored runs the Spearman correlation of PBIAS with log10 Kchan is %.3f "
            "(it should be strongly negative, below %.1f; Series 110 Stage 1 gave -0.994). Channel loss is probably "
            "not reaching tRIBS (optpercolation, the Kchan unit, or the builder). Stopping before more hours are spent. "
            "Send me the first run's .in file and this log." % (len(df), rho, INERT_RHO_MAX))
    return rho


# ======================================================================
# Truth A and the window
# ======================================================================
def find_truth(synth_dir):
    qout_files = list(synth_dir.glob("*.qout")) if synth_dir.exists() else []
    if len(qout_files) != 1:
        raise RuntimeError(
            "Expected exactly one *.qout file in %s (truth A), found %d: %s. Without it the scorer would silently "
            "fall back to real-gauge mode." % (synth_dir, len(qout_files), [f.name for f in qout_files]))
    return qout_files[0]


def truth_window_problems(ctx):
    """Truth A and the sweep must use the same simulation window. Returns a list of problems."""
    problems = []
    q = ctx.truth_info
    if not q["qout_ok"]:
        problems.append("truth A fails its sanity check: %s" % q["qout_problem"])
        return problems
    if abs(q["qout_tlast_hr"] - ctx.runtime_hours) > QOUT_END_SLACK_HR:
        problems.append(
            "truth A ends at hour %.2f but this script will run %g hours. Truth and sweep must use the same window. "
            "If the truth is the 1416 h one, run with --runtime_hours 1416 (about 3.7 minutes per run, roughly 24 h for "
            "400 runs); if the truth is 450 h, the builder's default is right." % (q["qout_tlast_hr"], ctx.runtime_hours))
    if abs(q["qout_dt_min"] - 3.75) > 0.05:
        problems.append("truth A's median time step is %.3f min, not the usual 3.75 min" % q["qout_dt_min"])
    return problems


# ======================================================================
# Design record, audit, report
# ======================================================================
def design_dict(ctx, samples):
    return {
        "series": SERIES, "category": CATEGORY, "script": Path(__file__).name,
        "script_md5": md5_of(Path(__file__)),
        "created_local": datetime.now().isoformat(timespec="seconds"),
        "what_this_is": "Series 122 Stage 1: Latin hypercube over cv, flowexp, n and Kchan with Ks and f fixed at "
                        "truth, channel loss ON, scored against truth A; every outlet file kept.",
        "seed": ctx.args.seed, "n": ctx.args.n,
        "box": {k: {kk: vv for kk, vv in v.items()} for k, v in LHS_PARAMS.items()},
        "fixed": {"Ks_mult": TRUTH_KS, "f_RS_abs": TRUTH_F, "optpercolation": 1,
                  "channelporosity": builder.BASELINE.get("channelporosity"),
                  "channelwidthcoeff": builder.BASELINE.get("channelwidthcoeff"),
                  "As_value": builder.BASELINE.get("As_value"), "Au_value": builder.BASELINE.get("Au_value"),
                  "optsnow": "0 (set inside the builder)"},
        "runtime_hours": ctx.runtime_hours,
        "truth_A": {"file": ctx.truth_path.name, "md5": ctx.truth_md5,
                    "rows": ctx.truth_info["qout_rows"], "last_hr": ctx.truth_info["qout_tlast_hr"]},
        "builder_md5": md5_of(Path(inspect.getsourcefile(builder))),
        "scorer_md5": md5_of(Path.cwd() / SCORER_SCRIPT),
        "python": sys.version.split()[0], "pandas": pd.__version__, "numpy": np.__version__,
        "pytRIBS": pytribs_version(),
        "draws": samples[PARAM_COLS].round(8).to_dict("list"),
    }


def same_draws(saved, samples):
    try:
        for k in PARAM_COLS:
            a = np.asarray(saved["draws"][k], dtype=float)
            b = samples[k].to_numpy(float)
            if len(a) != len(b) or not np.allclose(a, b, rtol=0, atol=1e-7):
                return False
        return True
    except Exception:
        return False


def audit_results(ctx, df, n_planned):
    """Check the results CSV against the kept outlet files. Returns a list of problems."""
    problems = []
    if df.empty:
        return ["the results file has no rows"]
    if df["run_id"].duplicated().any():
        problems.append("%d duplicated run_id(s)" % int(df["run_id"].duplicated().sum()))
    if len(df) != n_planned:
        problems.append("%d rows, but %d draws were planned" % (len(df), n_planned))
    n_bad_q = n_bad_c = n_bad_box = 0
    for _, r in df.iterrows():
        run_id = r["run_id"]
        q = ctx.kept_dir / ("%s_Outlet.qout" % run_id)
        if not q.exists():
            n_bad_q += 1
        elif "qout_md5" in df.columns and md5_of(q) != str(r.get("qout_md5")):
            n_bad_q += 1
        if not (ctx.csv_dir / ("%s_compare_obs_sim.csv" % run_id)).exists():
            n_bad_c += 1
        for p, b in LHS_PARAMS.items():
            v = _num(r.get(p))
            if not np.isfinite(v) or v < b["lo"] - 1e-9 or v > b["hi"] + 1e-9:
                n_bad_box += 1
    if n_bad_q:
        problems.append("%d row(s) whose kept outlet file is missing or whose checksum differs" % n_bad_q)
    try:
        ids = set(df["run_id"].astype(str))
        stray = [f.name for f in ctx.kept_dir.glob("*_Outlet.qout")
                 if f.name[:-len("_Outlet.qout")] not in ids]
        if stray:
            problems.append("%d outlet file(s) in %s with no row in the results file (e.g. %s)"
                            % (len(stray), ctx.kept_dir.name, stray[0]))
    except Exception:
        pass
    if n_bad_c:
        problems.append("%d row(s) with no compare CSV" % n_bad_c)
    if n_bad_box:
        problems.append("%d parameter value(s) outside the sampled box" % n_bad_box)
    for col in ("kge_2012", "pbias_pct"):
        if col in df.columns and not np.isfinite(pd.to_numeric(df[col], errors="coerce")).all():
            problems.append("non-finite values in %s" % col)
    return problems


def report_final(ctx, df, n_planned, n_failed):
    print("\n" + "=" * 70)
    print("SERIES %s STAGE 1 REPORT (Ks %g and f %g fixed; cv, r, n, Kchan free; channel loss ON)"
          % (SERIES, TRUTH_KS, TRUTH_F))
    print("=" * 70)
    print("  Scored runs: %d of %d planned; %d HANG/FAILED record(s) in the FAILED log." % (len(df), n_planned, n_failed))
    if len(df) < 3 or "pbias_pct" not in df.columns:
        return
    for p in PARAM_COLS:
        v = pd.to_numeric(df[p], errors="coerce")
        print("  %-26s coverage %.4f - %.4f" % (p, v.min(), v.max()))
    pb = pd.to_numeric(df["pbias_pct"], errors="coerce")
    kc = pd.to_numeric(df["channelconductivity_mmhr"], errors="coerce")
    rho = spearman(np.log10(kc), pb)
    print("\n  PBIAS (%%) against truth A: min %.1f, median %.1f, max %.1f" % (pb.min(), pb.median(), pb.max()))
    print("  Spearman correlation of PBIAS with log10 Kchan: %.3f  (Series 110 Stage 1: -0.994)" % rho)
    edges = [3, 10, 30, 100, 300, 1000.0001]
    print("\n  PBIAS by Kchan band (mm/hr; truth A has NO channel loss, so more Kchan = more missing volume):")
    print("  %12s  %5s  %9s  %9s" % ("band", "runs", "mean", "range"))
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (kc >= lo) & (kc < hi)
        if m.sum():
            print("  %5g-%-6g  %5d  %+9.1f  %+.1f to %+.1f" % (lo, min(hi, 1000), int(m.sum()), pb[m].mean(),
                                                               pb[m].min(), pb[m].max()))
    if "kge_2012" in df.columns:
        kg = pd.to_numeric(df["kge_2012"], errors="coerce")
        best = df.loc[kg.idxmax()]
        print("\n  KGE_2012 against truth A: best %.3f, median %.3f, worst %.3f" % (kg.max(), kg.median(), kg.min()))
        print("  Best run (reference only; the analysis script applies the reading rules): cv %.3f  r %.4f  n %.4f  "
              "Kchan %.1f  KGE_2012 %.3f  PBIAS %+.1f%%"
              % (best["kinemvelcoef"], best["flowexp"], best["channelroughness"],
                 best["channelconductivity_mmhr"], best["kge_2012"], best["pbias_pct"]))
    print("\n  Next: the analysis script (extend analyze_series121_partials_and_dial.py to four parameters). It reads the "
          "results file and the kept outlet files; truth B is scored from the kept files with no new runs.")


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
    return [r for r in rows if r.get("run_id") != new_row["run_id"]] + [new_row]


def result_row(ctx, metrics, i, row, qinfo, qout_md5):
    m = dict(metrics)
    m.update({
        "series": int(SERIES), "storm_leg": "1.00x", "storm_label": "storm100", "rain_scale": 1.0,
        "anchor_label": "kchan_Ks7p0_f0p012", "anchor_Ks": TRUTH_KS, "anchor_f": TRUTH_F,
        "batch": "kchan_lhs_stage1", "seed": ctx.args.seed, "draw_idx": int(i),
        "swept_param": "lhs_kchan_%s" % SERIES,
        "log10_Kchan": float(np.log10(row["channelconductivity_mmhr"])),
        "truth_file": ctx.truth_path.name, "truth_md5": ctx.truth_md5,
        "runtime_hours": ctx.runtime_hours,
        "qout_kept": "%s/%s_Outlet.qout" % (ctx.kept_dir.name, m["run_id"]),
        "qout_md5": qout_md5,
        "qout_rows": qinfo["qout_rows"] if qinfo else np.nan,
        "qout_last_hr": qinfo["qout_tlast_hr"] if qinfo else np.nan,
    })
    return m


def run_sweep(ctx, samples, todo):
    args = ctx.args
    out_path, failed_path = ctx.out_path, ctx.failed_path
    existing_df = load_existing_results(out_path)
    existing_ids = set(existing_df["run_id"].values) if not existing_df.empty else set()
    results = existing_df.to_dict("records") if not existing_df.empty else []
    if results and not args.skip_existing:
        raise Abort("%s already holds %d result row(s). Re-run with --skip_existing to resume, or move that file "
                    "(and the kept outlet folder) aside to start again. This script will not overwrite a sweep."
                    % (out_path.name, len(results)))

    failed_log_rows = []
    if args.skip_existing and failed_path.exists():
        try:
            failed_log_rows = pd.read_csv(failed_path).to_dict("records")
            print("  Loaded %d earlier HANG/FAILED record(s) from %s" % (len(failed_log_rows), failed_path.name))
        except Exception as e:
            print("  Warning: could not read %s (%s); starting a new log." % (failed_path.name, e))

    n_total = len(todo)
    print("\n" + "=" * 70)
    print("Series %s -- Kchan with routing free, Ks %g / f %g fixed, scored against truth A (%s), channel loss ON"
          % (SERIES, TRUTH_KS, TRUTH_F, ctx.truth_path.name))
    print("  %d of %d draws, seed=%d, window %g h, timeout=%ds, cleanup=%s"
          % (n_total, len(samples), args.seed, ctx.runtime_hours, args.timeout, "off" if args.no_cleanup else "on"))
    for p, b in LHS_PARAMS.items():
        print("  %-26s %g - %g  (%s)" % (p, b["lo"], b["hi"], b["scale"]))
    print("  Truth A: %s [md5 %s]" % (ctx.truth_path.name, ctx.truth_md5[:12]))
    print("  Output: %s   Kept outlet files: %s/" % (out_path.name, ctx.kept_dir.name))
    print("=" * 70 + "\n", flush=True)

    completed = skipped = hung = failed = cleaned = 0
    aborted = None
    sweep_start = time.time()
    inert_done = False

    def save():
        if results:
            df = pd.DataFrame(results)
            if "draw_idx" in df.columns:
                df = df.sort_values("draw_idx")
            write_csv_atomic(df, out_path)

    for i, row in todo.iterrows():
        tag = run_tag(i, row)
        run_id = predicted_run_id(tag)
        print("\n[%3d/%d]  cv=%.4f  r=%.5f  n=%.5f  Kchan=%.3f  -> %s"
              % (i + 1, len(samples), row["kinemvelcoef"], row["flowexp"], row["channelroughness"],
                 row["channelconductivity_mmhr"], run_id), flush=True)

        kept = ctx.kept_dir / ("%s_Outlet.qout" % run_id)
        csv_ok = (ctx.csv_dir / ("%s_compare_obs_sim.csv" % run_id)).exists()
        metrics_file = ctx.summary_dir / ("%s_metrics_summary.csv" % run_id)
        if args.skip_existing and csv_ok and kept.exists() and (run_id in existing_ids or metrics_file.exists()):
            print("  SKIP (already complete): %s" % run_id)
            skipped += 1
            if run_id not in existing_ids and metrics_file.exists():
                try:
                    m = pd.read_csv(metrics_file).iloc[0].to_dict()
                    qinfo = describe_qout(kept)
                    results.append(result_row(ctx, m, i, row, qinfo, md5_of(kept)))
                    existing_ids.add(run_id)
                    print("  (recovered its row from the metrics file)")
                except Exception as e:
                    print("  (could not recover its row: %s; it will be re-run next time)" % e)
            continue
        if args.skip_existing and run_id in existing_ids and not kept.exists():
            print("  (a result row exists but the kept outlet file is missing: re-running)")

        res = run_one(ctx, run_id, tag, row)
        status, reason = res["status"], res["reason"]

        if res.get("abort"):
            print("\n" + "!" * 70)
            print("STOPPING THE SWEEP (%s)" % status)
            print(res["abort"])
            print("Results so far are saved.")
            print("!" * 70)
            failed_log_rows = upsert_failed(failed_log_rows, {
                "run_id": run_id, "status": status, "reason": reason, "elapsed_min": res["elapsed"] / 60,
                "draw_idx": int(i), **{p: row[p] for p in PARAM_COLS}})
            write_csv_atomic(pd.DataFrame(failed_log_rows), failed_path)
            failed += 1
            aborted = res["abort"]
            break

        if status == "SUCCESS":
            metrics = result_row(ctx, res["metrics"], i, row, res["qinfo"], res["qout_md5"])
            results = [r for r in results if r.get("run_id") != run_id]
            results.append(metrics)
            existing_ids.add(run_id)
            completed += 1
            note_outcome(ctx, True)
            if any(r.get("run_id") == run_id for r in failed_log_rows):
                failed_log_rows = [r for r in failed_log_rows if r.get("run_id") != run_id]
                write_csv_atomic(pd.DataFrame(failed_log_rows, columns=FAILED_COLS), failed_path)
            print("  KGE_2012=%.3f  PBIAS=%+.1f%%   (outlet file kept: %d rows, md5 %s)"
                  % (_num(metrics.get("kge_2012")), _num(metrics.get("pbias_pct")),
                     metrics["qout_rows"], str(res["qout_md5"])[:10]))
            if not args.no_cleanup and cleanup_raw_results(run_id, ctx.calib_dir):
                cleaned += 1
        else:
            print("  %s: %s  (%s)" % (status, run_id, reason))
            if status == "HANG":
                hung += 1
            else:
                failed += 1
            failed_log_rows = upsert_failed(failed_log_rows, {
                "run_id": run_id, "status": status, "reason": reason, "elapsed_min": res["elapsed"] / 60,
                "draw_idx": int(i), **{p: row[p] for p in PARAM_COLS}})
            write_csv_atomic(pd.DataFrame(failed_log_rows), failed_path)
            note_outcome(ctx, False)

        save()

        if status == "SUCCESS" and not inert_done and len(results) >= INERT_CHECK_AT:
            try:
                rho = inert_check(ctx, results)
            except Abort as e:
                aborted = str(e)
                print("\n" + "!" * 70)
                print("STOPPING THE SWEEP (KCHAN INERT)")
                print(aborted)
                print("Results so far are saved.")
                print("!" * 70)
                break
            inert_done = True
            print("  [check] Kchan moves PBIAS as expected: Spearman %.3f over %d runs (needs < %.1f)"
                  % (rho, len(results), INERT_RHO_MAX), flush=True)

        attempted = completed + hung + failed
        remaining = n_total - attempted - skipped
        if attempted > 0:
            avg = (time.time() - sweep_start) / attempted
            eta_min = avg * max(remaining, 0) / 60
            print("  Run time: %.1f min  |  ETA: %.0f min (%.1f h) remaining"
                  % (res["elapsed"] / 60, eta_min, eta_min / 60), flush=True)

    print("\nSweep %s: %d ran, %d skipped, %d hung, %d failed, %d raw folders cleaned"
          % ("STOPPED EARLY" if aborted else "complete", completed, skipped, hung, failed, cleaned))
    save()
    if results:
        final_df = pd.DataFrame(results).sort_values("draw_idx")
        print("Saved: %s  (%d rows)" % (out_path.name, len(final_df)))
        try:
            full = (not args.limit) and len(final_df) == len(samples)
            if full and not aborted:
                probs = audit_results(ctx, final_df, len(samples))
                print("AUDIT: %s" % ("PASS (every row has its kept outlet file with a matching checksum, a compare "
                                     "CSV, and parameters inside the box)" if not probs
                                     else "PROBLEMS: " + "; ".join(probs)))
            report_final(ctx, final_df, len(samples), len(failed_log_rows))
        except Exception as e:
            print("\n(report skipped: %s -- results are saved regardless)" % e)
    else:
        print("No results to save.")
    return 3 if aborted else 0


# ======================================================================
# Dry run
# ======================================================================
def print_dry_run(ctx, samples, raw, checks, competing, window_problems):
    n = len(samples)
    est_sec = SEC_PER_MODEL_HOUR * ctx.runtime_hours
    print("\n" + "=" * 92)
    print("DRY RUN -- Series %s Stage 1. Nothing is written, built or run." % SERIES)
    print("=" * 92)
    print("  Truth A (read-only): %s  [md5 %s]" % (ctx.truth_path.name, ctx.truth_md5[:12]))
    q = ctx.truth_info
    print("    %d rows, %d columns, step %.2f min, last time %.2f h, peak %.4g m3/s"
          % (q["qout_rows"], q["qout_ncols"], q["qout_dt_min"], q["qout_tlast_hr"], q["qout_peak_m3s"]))
    print("  Simulation window this script will use: %g h (builder default %g h%s)"
          % (ctx.runtime_hours, ctx.builder_runtime_default,
             "; overridden by --runtime_hours" if ctx.runtime_hours != ctx.builder_runtime_default else ""))
    if window_problems:
        print("  WINDOW CHECK: PROBLEM")
        for p in window_problems:
            print("    - %s" % p)
    else:
        print("  WINDOW CHECK: OK (truth A's last time %.2f h matches the %g h window)" % (q["qout_tlast_hr"], ctx.runtime_hours))
    print("  Builder: %s" % inspect.getsourcefile(builder))
    print("    unit guard OK (no division of Kchan by 3.6e6); overrides and tag accepted; "
          "run ids will look like %s" % predicted_run_id(run_tag(0, samples.iloc[0])))
    print("    held at the builder's values: channelporosity %s, channelwidthcoeff %s, As %s, Au %s; "
          "optsnow 0 (set inside the builder)"
          % (builder.BASELINE.get("channelporosity"), builder.BASELINE.get("channelwidthcoeff"),
             builder.BASELINE.get("As_value"), builder.BASELINE.get("Au_value")))
    print("\n  DESIGN: %d draws, seed %d. Ks %g and f %g fixed; optpercolation 1 (channel loss ON)."
          % (n, ctx.args.seed, TRUTH_KS, TRUTH_F))
    print("  LHS check (exactly one draw in each of the %d slices of each parameter): %s"
          % (n, ", ".join("%s %s" % (k, "PASS" if v else "FAIL") for k, v in checks.items())))
    for p, b in LHS_PARAMS.items():
        print("    %-26s box %g - %g (%s); draws span %.5g - %.5g; truth/true value %s"
              % (p, b["lo"], b["hi"], b["scale"], samples[p].min(), samples[p].max(),
                 TRUTH_ROUTING.get(p, "(no truth: truth A has channel loss OFF)")))
    kc = samples["channelconductivity_mmhr"]
    print("    Kchan draws: %d below 30, %d from 30 to 300, %d above 300 mm/hr; %d within a factor 2 of 201"
          % ((kc < 30).sum(), ((kc >= 30) & (kc <= 300)).sum(), (kc > 300).sum(), ((kc > 100.5) & (kc < 402)).sum()))
    ids = [predicted_run_id(run_tag(i, r)) for i, r in samples.iterrows()]
    print("    run ids: %d draws, %d different (longest %d characters)" % (n, len(set(ids)), max(len(x) for x in ids)))
    print("    first draw: %s" % ids[0])
    print("    state now: %s" % ctx.state_line)
    print("\n  TIME: about %.0f s per run at this window -> %d runs ~ %.1f h (a rough estimate from the Series 110 "
          "timings: 67 s at 450 h, 219 s at 1416 h)." % (est_sec, n, n * est_sec / 3600.0))
    print("  Disk: %.1f GB free (a run needs about 1 GB while it is running; each kept outlet file is under 1 MB)."
          % (shutil.disk_usage(ctx.calib_dir).free / 1e9))
    print("  Kept outlet files go to %s/ (not git-ignored; add by name or folder when you commit)." % ctx.kept_dir.name)
    if competing:
        print("\n  NOTE: another tRIBS build/run script appears to be running right now:")
        for pid, cmd in competing:
            print("    pid %d: %s" % (pid, cmd[:110]))
        print("  (A dry run is harmless, but the sweep would refuse to start.)")
    print("\n  Truth B (Kchan 201) is not used by this script; it is scored later from the kept outlet files.")
    print("\n--dry_run: nothing written, built or run.")
    if window_problems:
        print("\nNOT READY: fix the window problem above first.")
    else:
        print("\nReady for the trial:  python -u %s --limit 3 2>&1 | tee series122_trial.log" % Path(__file__).name)


# ======================================================================
# Main
# ======================================================================
class Ctx(object):
    pass


def main():
    parser = argparse.ArgumentParser(
        description="Series 122 Stage 1: Latin hypercube over cv, flowexp, n and Kchan (channel loss ON), Ks and f "
                    "fixed at truth, scored against truth A, every outlet file kept. tRIBS 6.0.0.")
    parser.add_argument("--dry_run", action="store_true",
                        help="Print the design, the checks and the time estimate, then stop. Writes, builds and runs nothing")
    parser.add_argument("--audit", action="store_true",
                        help="Check the results file against the kept outlet files (no tRIBS)")
    parser.add_argument("--n", type=int, default=DEFAULT_N,
                        help="Number of draws (default %d; changing it changes EVERY draw)" % DEFAULT_N)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED,
                        help="Random seed (default %d)" % DEFAULT_SEED)
    parser.add_argument("--limit", type=int, default=None,
                        help="Only process the first N draws (a short trial; resume later with --skip_existing)")
    parser.add_argument("--skip_existing", action="store_true",
                        help="Resume: skip draws whose compare CSV, kept outlet file and result row all exist")
    parser.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT, help="Per-run hard timeout in seconds (default %d)" % DEFAULT_TIMEOUT)
    parser.add_argument("--runtime_hours", type=float, default=None,
                        help="Simulation window in hours. Default: the builder's own value. Must match truth A's window")
    parser.add_argument("--no_cleanup", action="store_true",
                        help="Keep raw per-run results folders instead of deleting them after a successful run")
    parser.add_argument("--ignore_running_check", action="store_true",
                        help="Skip the check for other tRIBS scripts running (only if you are sure none are)")
    parser.add_argument("--max_consecutive_failures", type=int, default=MAX_CONSEC_FAIL,
                        help="Stop if this many runs in a row fail (default %d)" % MAX_CONSEC_FAIL)
    args = parser.parse_args()

    if args.n < 2:
        parser.error("--n must be at least 2")
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be at least 1")
    if sum([args.dry_run, args.audit]) > 1:
        parser.error("use only one of --dry_run, --audit")

    ctx = Ctx()
    ctx.args = args
    ctx.consec_fail = 0
    ctx.in_note_printed = False
    cwd = Path.cwd()
    ctx.calib_dir = cwd.parent / "calibration_work"
    ctx.summary_dir = ctx.calib_dir / "03_comparisons" / "summary_tables"
    ctx.csv_dir = ctx.calib_dir / "03_comparisons" / "csv_exports"
    variant = variant_suffix(args.n, args.seed)
    ctx.kept_dir = ctx.calib_dir / ("kept_qout_%s%s" % (SERIES, variant))
    ctx.out_path = ctx.summary_dir / ("%s%s.csv" % ("lhs_results_kchan_" + SERIES, variant))
    ctx.failed_path = ctx.summary_dir / ("%s%s.csv" % ("lhs_results_kchan_FAILED_" + SERIES, variant))
    ctx.design_path = ctx.summary_dir / ("%s%s.json" % ("lhs_design_kchan_" + SERIES, variant))

    print("Series %s Stage 1 (Kchan, routing free)   python %s, pandas %s, numpy %s"
          % (SERIES, sys.version.split()[0], pd.__version__, np.__version__))

    # ------------------------------------------------------------------
    # HARD SAFETY CHECKS
    # ------------------------------------------------------------------
    if not (cwd / SCORER_SCRIPT).exists():
        raise RuntimeError("%s not found in %s. Run this script from the lab/ directory." % (SCORER_SCRIPT, cwd))
    if not ctx.calib_dir.exists():
        raise RuntimeError("%s not found. Run this script from the lab/ directory (the folder that sits next to "
                           "calibration_work/)." % ctx.calib_dir)
    ctx.builder_runtime_default = float(builder.RUNTIME_HOURS)
    if args.runtime_hours is not None:
        builder.RUNTIME_HOURS = args.runtime_hours
    ctx.runtime_hours = float(builder.RUNTIME_HOURS)
    builder_unit_guard()
    install_series_hooks()

    ctx.truth_path = find_truth(ctx.calib_dir / "synth_truth")
    ctx.truth_md5 = md5_of(ctx.truth_path)
    ctx.truth_info = describe_qout(ctx.truth_path)
    print("Truth A, read-only: %s  [md5 %s]" % (ctx.truth_path.name, ctx.truth_md5[:12]))
    window_problems = truth_window_problems(ctx)

    competing = None
    if args.ignore_running_check:
        print("Running-jobs check skipped (--ignore_running_check).")
    else:
        competing = find_competing_jobs()
        if competing is None:
            print("WARNING: could not check for other running tRIBS scripts (ps unavailable). Make sure none are "
                  "running -- they share current_run_config.json.")
        elif competing and not (args.dry_run or args.audit):
            print("\nANOTHER tRIBS BUILD/RUN SCRIPT APPEARS TO BE RUNNING:")
            for pid, cmd in competing:
                print("  pid %d: %s" % (pid, cmd[:110]))
            raise RuntimeError(
                "Refusing to start: these scripts share current_run_config.json with no locking, so running both can "
                "corrupt both. Wait for them to finish (or stop them with pkill -f <script>), or re-run with "
                "--ignore_running_check if you are sure they are stale.")
        elif not competing:
            print("No other tRIBS build/run script detected.")

    raw = generate_lhs_raw(args.n, LHS_PARAMS, args.seed)
    checks = lhs_check(raw, LHS_PARAMS, args.n)
    samples = round_samples(raw, LHS_PARAMS)
    todo = samples.iloc[:args.limit] if args.limit else samples

    done_rows = 0
    if ctx.out_path.exists():
        try:
            done_rows = len(pd.read_csv(ctx.out_path))
        except Exception:
            done_rows = -1
    ctx.state_line = ("no results file yet" if not ctx.out_path.exists()
                      else "%s holds %s row(s)%s" % (ctx.out_path.name, "an unreadable number of" if done_rows < 0 else done_rows,
                                                    "" if args.skip_existing else " (a sweep needs --skip_existing to resume)"))

    if args.audit:
        df = load_existing_results(ctx.out_path)
        probs = audit_results(ctx, df, len(samples))
        print("\nAUDIT of %s: %s" % (ctx.out_path.name, "PASS" if not probs else "PROBLEMS"))
        for p in probs:
            print("  - %s" % p)
        return 0 if not probs else 1

    if args.dry_run:
        print_dry_run(ctx, samples, raw, checks, competing, window_problems)
        return 2 if window_problems else 0

    if window_problems:
        raise RuntimeError("the truth and the simulation window do not agree: " + "; ".join(window_problems))
    if not all(checks.values()):
        raise RuntimeError("the Latin hypercube self-check failed for %s" % [k for k, v in checks.items() if not v])

    ctx.summary_dir.mkdir(parents=True, exist_ok=True)
    ctx.csv_dir.mkdir(parents=True, exist_ok=True)
    ctx.kept_dir.mkdir(parents=True, exist_ok=True)

    try:
        # The design on disk must be the design in memory (resume safety).
        if ctx.design_path.exists():
            try:
                saved = json.loads(ctx.design_path.read_text())
            except Exception as e:
                raise Abort("cannot read the saved design %s (%s). Move it aside if you are sure you want a new sweep."
                            % (ctx.design_path.name, e))
            if not same_draws(saved, samples):
                raise Abort("the draws now (seed %d, n %d, box) differ from the saved design in %s. A resumed sweep must "
                            "use the same draws. Use the same --seed/--n, or move the old results aside."
                            % (args.seed, args.n, ctx.design_path.name))
            if saved.get("truth_A", {}).get("md5") != ctx.truth_md5:
                raise Abort("truth A (%s) has a different checksum from the one this sweep started with (%s). Do not "
                            "mix runs scored against two different truths." % (ctx.truth_md5[:12], str(saved.get("truth_A", {}).get("md5"))[:12]))
            if saved.get("runtime_hours") is not None and not _close(saved["runtime_hours"], ctx.runtime_hours, 1e-9):
                raise Abort("this sweep started with a %g h window; now it is %g h." % (saved["runtime_hours"], ctx.runtime_hours))
        else:
            write_json_atomic(ctx.design_path, design_dict(ctx, samples))
            print("Saved the design record: %s" % ctx.design_path.name)

        return run_sweep(ctx, samples, todo)
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


def cli():
    """Turn the script's own deliberate refusals into a readable STOPPED box instead of a traceback.
    Only RuntimeError is caught (that is what every refusal raises); a real bug still shows its traceback."""
    try:
        return main()
    except RuntimeError as e:
        print("\n" + "!" * 92)
        print("STOPPED: %s" % e)
        print("!" * 92)
        print("Nothing was changed by this refusal. Fix the cause above, then run the command again.")
        return 1


if __name__ == "__main__":
    sys.exit(cli())
