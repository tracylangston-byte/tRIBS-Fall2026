"""
run_cc_clean_122_tribs6.py
==========================
Series 122 follow-up -- twelve "clean" runs that change ONE thing at a time, on tRIBS 6.0.0 /
pytRIBS 1.0.0.

WHY
---
The 400-run Series 122 sweep changed four things at once (cv, flowexp, n and Kchan), so it could
not show what Kchan alone does to the hydrograph. And it cannot say whether the fingerprint of
channel loss (volume -10%, peak -4%, timing unchanged) is different from the fingerprint of an
error in the rain. These twelve runs answer both, with the routing parameters held at their truth
values every time:

    cv 4.5   flowexp 0.24   n 0.026   Ks 7.0   f 0.012      (the same truth as Series 121/122)

THE TWELVE RUNS (about 1 minute each, about 14 minutes in all)
--------------------------------------------------------------
  1. ctlA   control, channel loss OFF, rain as recorded. Must reproduce truth A.
  2. ctlB   control, channel loss ON, Kchan 201 mm/hr. Must reproduce truth B.
  3-7.      the rain comparison: channel loss OFF, the rain multiplied by 0.90, 0.93, 0.95, 0.97
            and 0.99 (--rain to change them). A rain that is a few percent too low also takes
            volume off the hydrograph. These runs show whether it takes it off in the same way
            as channel loss does (same effect on the peak and the timing) or in a different way.
  8-12.     the Kchan dose-response: channel loss ON, Kchan = 10, 30, 100, 300, 1000 mm/hr
            (--kchan to change them). Nothing else changes, so every difference from ctlA is
            caused by Kchan alone. (201 is ctlB; it is the sixth point of the same curve.)

Controls come first, then the rain runs (the part with new plumbing, so a problem shows up at run 3
and not at run 8), then the Kchan runs. A trial with --limit 3 runs both controls and the first
rain run (rain x0.90).

THE TWO CONTROLS ARE A SAFETY GATE
----------------------------------
The other ten runs only mean something if today's build still reproduces the two truths, so each
control must pass before the script goes on (it stops with a plain message if one fails):

  ctlA  scored against truth A: |PBIAS| < 0.01% and KGE_2012 > 0.9999 (the Series 110 rule).
  ctlB  (a) scored against truth B by this script's own formulas: |PBIAS| < 0.01% and
        KGE_2012 > 0.9999; and (b) scored against truth A by the project scorer, it must give the
        numbers Series 110 recorded for truth B (PBIAS -9.98684%, KGE_2012 0.88414).
  Both controls also print whether the new outlet file is byte-identical to the truth file
  (information only).

OTHER CHECKS (any one that fails stops the script, or fails that one run, as stated)
-------------------------------------------------------------------------------------
  * truth A (the single *.qout in synth_truth/) and truth B must have the md5 fingerprints they
    had when Series 122 was scored, and the same 450 h window as the runs. [stops]
  * the builder must not divide Kchan by 3.6e6 (the Series 110 unit bug) and must accept
    overrides, a tag and a gauge_sdf. [stops]
  * after every build the builder's run configuration is read back and must carry exactly the
    Ks, f, cv, flowexp, n, optpercolation (and Kchan for runs with loss) asked for, no truth_file,
    and (rain runs) the scaled rain file; the .in file is checked for optpercolation, Kchan and
    the rain file name where its layout allows. [stops]
  * the scorer must have used truth A as "observed" (obs_mode synth, its own "Reading observed
    from:" line, observed peak equal to truth A's) and the scored row must carry the parameters
    asked for. [stops]
  * the kept outlet file is checked by checksum and by re-scoring it here; this must agree with the
    scorer's own PBIAS and KGE_2012 (otherwise the kept file is not the run that was scored). [fails that run]
  * rain runs: the rain files are copies of the recorded rain with only the rain-rate column
    multiplied, checked against the originals row by row; and a rain run must give less volume AND
    a smaller peak than truth A (otherwise tRIBS did not use the scaled rain). [stops]
  * each run is a subprocess with a hard timeout (default 300 s); a hang kills its process
    group. HANG / FAILED runs go to the FAILED log and the script moves on; 4 in a row stop it.
  * refuses to start if another tRIBS build/run script seems to be running (they share
    calibration_work/current_run_config.json with no locking).
  * on resume, the design must be identical to the saved design (same runs, truths and rain files).
  * at least 2 GB of free disk before every run.

WHERE THINGS GO
---------------
    run id            SMF_20140812_122_Ks7p0x_clean_<what>
                         clean_ctlA_noloss | clean_ctlB_cc201p0 | clean_cc<Kchan> | clean_rain<scale>
                         (examples: ..._clean_cc30p0   ..._clean_rain0p95)
    results CSV       calibration_work/03_comparisons/summary_tables/series122_clean_runs.csv
    failed log        calibration_work/03_comparisons/summary_tables/series122_clean_runs_FAILED.csv
    design record     calibration_work/03_comparisons/summary_tables/series122_clean_runs_design.json
    kept outlet files calibration_work/kept_qout_122_clean/<run id>_Outlet.qout   (one per run)
    rain provenance   calibration_work/forcing_122_clean_PROVENANCE.json
    scaled rain files ../init_data/met/precip_<station>_rain<scale>.mdf and Master_Precip_rain<scale>.sdf
                      (one set per rain scale, e.g. Master_Precip_rain0p95.sdf; the originals are
                      never touched; an existing file with different contents is never overwritten)
    raw tRIBS output  calibration_work/02_results/122_clean_runs/<run id>/   (deleted after the
                      outlet file has been copied and the copy checked by checksum)
    run inputs        calibration_work/01_run_inputs/122_clean_runs/
The summary_tables/ folder is under the git-ignored 03_comparisons/ folder, so committing the CSV
needs `git add -f`. kept_qout_122_clean/ and the rain files are not ignored.

Every run is scored by the project scorer against truth A (the single *.qout in synth_truth/), the
same as the sweep. Nothing is scored against truth B by the scorer; the comparison script
(plot_series122_clean_runs.py) re-scores the kept files against both truths. Nothing in
synth_truth/ or synth_truth_location/ is written or moved; both are only READ.

The file is called run_cc_... on purpose: the other tRIBS scripts look for that prefix to detect
a running job, so they will refuse to start while this one is running (and vice versa).

USAGE (run from the lab/ directory, like every other run script)
----------------------------------------------------------------
    cd /workspaces/tRIBS-Fall2026/workspaces/SMF/lab
    python -u run_cc_clean_122_tribs6.py --dry_run
    python -u run_cc_clean_122_tribs6.py --limit 3 2>&1 | tee series122_clean_trial.log
    python -u run_cc_clean_122_tribs6.py --skip_existing 2>&1 | tee series122_clean.log
    python -u run_cc_clean_122_tribs6.py --audit

  --dry_run         writes, builds and runs nothing; prints the plan, the checks and the time
  --limit N         only the first N runs (the two controls come first, then the rain runs)
  --skip_existing   resume: a run is skipped only if its row, its kept outlet file (checksum
                    unchanged) and its compare CSV all exist; a control is skipped only if it passed
  --audit           reads the results and the kept outlet files and checks them (no tRIBS)
  --kchan 10,30,100,300,1000     the Kchan values (mm/hr) of the dose-response runs
  --rain 0.90,0.93,0.95,0.97,0.99   the rain multipliers of the rain runs
  (Changing --kchan or --rain after runs exist changes the design; the script then refuses to
   resume. Move the old results aside first.)

Do not run any other tRIBS build/run script at the same time. The -u matters: without it Python
holds output back when it is piped to tee.

EXIT CODES: 0 = done (or a clean dry run); 1 = refused to start; 2 = dry run found a problem;
            3 = stopped by a safety check; 130 = Ctrl+C.

STATUS
------
Written without access to the Codespace. NOT YET RUN against the real tRIBS 6.0.0 build.

What WAS tested (on stand-ins: a toy hydrograph model plus copies of the real scorer's
formulas, in a pretend lab folder -- not tRIBS):
  * dry run; a 3-run trial; a full 12-run pass; resume; --audit; every number cross-checked
    against an independent recomputation (agreement to ~1e-13).
  * stops and refusals: wrong or edited truth files, a missing truth B, changed rain files,
    a rain scale the model ignores (inert), a channel loss the model ignores, a builder that
    mis-handles units or lacks the gauge_sdf argument, a short outlet file, a scorer that
    disagrees with the recorded truth-B score, a hung run (--timeout), a failed build,
    refusing to overwrite old results, a changed --kchan/--rain list on resume, an edited
    kept file (that run is redone on resume), a competing run_lhs_* job, a low-disk guard,
    Ctrl+C (exit 130, scorer stopped, finished runs kept), and 6 ways of tampering with
    the results (all caught by --audit).
What was NOT tested (can only be learned on the real Codespace):
  * the real rain-file (.mdf/.sdf) layout, and whether the real builder honours gauge_sdf
    and the real scorer reads the scaled rain (the dry run and the two control runs will
    show this; the script stops if the rain scale turns out to be inert);
  * whether the real .in file layout allows the optpercolation/Kchan read-back check
    (if not, the script prints a note and skips only that check);
  * real run times (estimate ~68 s per run, ~14 min for 12 runs).
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
CATEGORY      = "122_clean_runs"
SCORER_SCRIPT = "run_sensitivity_single_interp_tribs6.py"
KEPT_DIR_NAME = "kept_qout_122_clean"
RESULTS_NAME  = "series122_clean_runs.csv"
FAILED_NAME   = "series122_clean_runs_FAILED.csv"
DESIGN_NAME   = "series122_clean_runs_design.json"
FORCING_PROV_NAME = "forcing_122_clean_PROVENANCE.json"

# ----------------------------------------------------------------------
# The truth every run is built around (routing at truth, Ks and f at truth)
# ----------------------------------------------------------------------
TRUTH_KS      = 7.0
TRUTH_F       = 0.012
TRUTH_ROUTING = {"kinemvelcoef": 4.5, "flowexp": 0.24, "channelroughness": 0.026}
TRUTH_B_KCHAN = 201.0

# Fingerprints of the two truth files as they were when Series 122 was scored.
TRUTH_A_MD5      = "4a433e441359313b6629734510f27813"
TRUTH_B_DIR_NAME = "synth_truth_location"
TRUTH_B_NAME     = "SMF_20140812_110_cc201p0_loc_K7p0_f0p012_Outlet.qout"
TRUTH_B_MD5      = "013097a414def7b59aca731f9665f632"
# What Series 110 recorded for truth B scored against truth A (the project scorer's numbers).
REC_B_VS_A = {"pbias_pct": -9.98684101102667, "kge_2012": 0.8841407468110921}

DEFAULT_KCHAN = (10.0, 30.0, 100.0, 300.0, 1000.0)
DEFAULT_RAIN  = (0.90, 0.93, 0.95, 0.97, 0.99)

# ----------------------------------------------------------------------
# Safety settings (not scientific reading rules)
# ----------------------------------------------------------------------
DEFAULT_TIMEOUT    = 300
SEC_PER_MODEL_HOUR = 0.152          # about 67 s per run at 450 h (Series 110 timings)
QOUT_END_SLACK_HR  = 2.0            # an outlet file must reach to within this many hours of the window end
MIN_FREE_GB        = 2.0
OBS_PEAK_RTOL      = 0.02
EXPECTED_QOUT_COLS = 3
MAX_CONSEC_FAIL    = 4
CTL_PBIAS_TOL      = 0.01           # percent   (the Series 110 provenance rule)
CTL_KGE_MIN        = 0.9999
REC_PBIAS_TOL      = 1e-3           # percentage points: ctlB scored against truth A vs the Series 110 record
REC_KGE_TOL        = 1e-5
LOCAL_PBIAS_TOL    = 0.01           # percentage points: this script's re-scoring vs the project scorer
LOCAL_KGE_TOL      = 1e-4
RAIN_MIN_CHANGE    = 0.001          # a rain run must move volume and peak by more than 0.1%
SERIES_ABS_TOL_M3S = 1e-6           # "numerically the same" over the whole run (information only)
ORIGIN             = "2014-08-01"   # the .qout Time column is hours from this date
KCHAN_RANGE        = (0.1, 1.0e5)
RAIN_RANGE         = (0.5, 1.5)

# ----------------------------------------------------------------------
# Rain-file constants (same as the Series 120 storm script)
# ----------------------------------------------------------------------
SRC_SDF_DEFAULT = "../init_data/met/Master_Precip.sdf"
RAIN_DT_HR  = float(getattr(builder, "RAIN_INTERVAL", 0.25))   # one MDF row = this many hours
NO_DATA_MIN = 9999.0                                           # tRIBS NO_DATA flag is 9999.99

FAILED_COLS = ["run_id", "status", "reason", "elapsed_min", "run_name", "run_kind", "kchan_asked", "rain_scale"]

COMPETING_PATTERNS = (
    "run_sensitivity_single", "run_sensitivity_sweep", "build_sensitivity_run",
    "run_cc_", "run_lhs_", "run_lowf_sweep", "run_crossversion", "verify_truth_provenance",
    "validate_truth_point", "rescore_series100", "run_truth_", "verify_truth_B_",
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


def _fmt(x, spec="%.4f"):
    v = _num(x)
    return spec % v if np.isfinite(v) else "-"


def write_json_atomic(path, obj):
    tmp = path.with_name(path.name + ".part")
    tmp.write_text(json.dumps(obj, indent=2, default=str))
    os.replace(tmp, path)


def write_csv_atomic(df, path):
    tmp = path.with_name(path.name + ".part")
    df.to_csv(tmp, index=False)
    os.replace(tmp, path)


def label(v):
    """The builder's own number-to-text rule (7.0 -> 7p0, 0.012 -> 0p012, 0.95 -> 0p95)."""
    return builder.value_to_label(v)


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


# ======================================================================
# Names, and teaching the (unedited) builder about series 122 and this category
# ======================================================================
def install_series_hooks():
    """In-memory only: Ks_mult runs are series 122, and series 122 lives in this category.
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
    collapsed to nearly zero. Stop if an executable line still does that. Also stop if the builder
    cannot do what this script needs."""
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
    for p in ("overrides", "tag", "gauge_sdf"):
        if p not in params:
            raise RuntimeError("builder.build_input_file has no %r argument. This script needs the tRIBS-6 builder "
                               "that has it (the rain runs need gauge_sdf)." % p)


def predicted_run_id(tag):
    run_id, _ = builder.build_run_id("Ks_mult", TRUTH_KS)
    return "%s_%s" % (run_id, tag)


def derived_name(path_str, suffix):
    """'../init_data/met/precip_SMF_2014-2014.mdf' + 'rain0p95' -> '../init_data/met/precip_SMF_2014-2014_rain0p95.mdf'."""
    head, sep, name = path_str.rpartition("/")
    stem, dot, ext = name.rpartition(".")
    if not dot:
        return "%s%s%s_%s" % (head, sep, name, suffix)
    return "%s%s%s_%s.%s" % (head, sep, stem, suffix, ext)


# ======================================================================
# The plan: twelve runs
# ======================================================================
def plan_runs(kchans, rains, src_sdf_rel):
    """The list of runs, controls first. Pure: builds names only, touches nothing."""
    runs = [
        {"name": "ctlA", "kind": "control_noloss", "optp": 0, "kchan": None, "rain": 1.0,
         "tag": "clean_ctlA_noloss"},
        {"name": "ctlB", "kind": "control_loss", "optp": 1, "kchan": TRUTH_B_KCHAN, "rain": 1.0,
         "tag": "clean_ctlB_cc%s" % label(TRUTH_B_KCHAN)},
    ]
    # The rain runs come before the Kchan runs on purpose: they are the part with new plumbing (scaled rain
    # files), so if something is wrong with it the script stops at run 3, not run 8.
    for s in rains:
        runs.append({"name": "rain%s" % label(s), "kind": "rain", "optp": 0, "kchan": None, "rain": float(s),
                     "tag": "clean_rain%s" % label(s)})
    for k in kchans:
        runs.append({"name": "K%s" % label(k), "kind": "kchan", "optp": 1, "kchan": float(k), "rain": 1.0,
                     "tag": "clean_cc%s" % label(k)})
    for i, r in enumerate(runs):
        r["idx"] = i
        r["run_id"] = predicted_run_id(r["tag"])
        r["rain_label"] = "rain%s" % label(r["rain"]) if r["kind"] == "rain" else None
        r["sdf_rel"] = derived_name(src_sdf_rel, r["rain_label"]) if r["kind"] == "rain" else None
    return runs


def parse_float_list(text, what, lo, hi):
    """'10,30,100' -> [10.0, 30.0, 100.0]. Sorted, no duplicates, inside (lo, hi)."""
    try:
        vals = [float(x) for x in re.split(r"[,\s]+", text.strip()) if x]
    except ValueError:
        raise ValueError("--%s must be numbers separated by commas, e.g. %s; got %r"
                         % (what, "10,30,100" if what == "kchan" else "0.9,0.95", text))
    if not vals:
        raise ValueError("--%s is empty" % what)
    for v in vals:
        if not np.isfinite(v) or v < lo or v > hi:
            raise ValueError("--%s value %g is outside the allowed range %g to %g" % (what, v, lo, hi))
    vals = sorted(vals)
    labels = [label(v) for v in vals]
    if len(set(labels)) != len(labels):
        raise ValueError("--%s has duplicate values (after rounding for the file names): %s" % (what, vals))
    return vals


def overrides_for(run):
    """Everything that must be pinned or set for a run (the rest comes from the builder's defaults)."""
    ov = {
        "Ks_mult":          TRUTH_KS,
        "f_RS_abs":         TRUTH_F,
        "kinemvelcoef":     TRUTH_ROUTING["kinemvelcoef"],
        "flowexp":          TRUTH_ROUTING["flowexp"],
        "channelroughness": TRUTH_ROUTING["channelroughness"],
        "optpercolation":   run["optp"],
    }
    if run["kchan"] is not None:
        ov["channelconductivity_mmhr"] = float(run["kchan"])
    return ov


# ======================================================================
# FORCING: scaled copies of the two rain files and the SDF
# (the same code as run_lhs_Ks_f_120_storms_tribs6.py, which handles the real header
#  Year,Month,Day,Hour,R_mm as well as the old  Y M D H R)
# ======================================================================
MDF_TOKEN = re.compile(r"[^\s,]+\.mdf\b", re.IGNORECASE)     # a path: no spaces, no commas
_SEP = re.compile(r"(\s*,\s*|\s+)")                          # a column separator: comma or whitespace
_R_NAMES = {"R", "RAIN", "RAINFALL", "PRECIP", "PRECIPITATION"}


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


def rain_column(header_tokens, where, first_lines):
    """Find the rainfall-rate column from a rain file's header, and confirm there is a year column.
    The tRIBS 6.0.0 header is  Year,Month,Day,Hour,R_mm  (comma-separated); older files had  Y M D H R.
    Anything else stops with the first lines of the file printed, so nothing is guessed."""
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
    print("    rain x%g: SDF %s -> %s" % (plan["scale"], plan["src_sdf"].name, plan["dst_sdf"].name))
    for s in plan["stations"]:
        st = s["stats"]
        print("      %-36s %5d rows   peak %7.3f -> %7.3f mm/hr   total %8.2f -> %8.2f mm"
              % (s["dst"].name, st["rows"], st["peak_before"], st["peak_after"],
                 st["sum_before"] * RAIN_DT_HR, st["sum_after"] * RAIN_DT_HR))
        if st["kept_flag"]:
            print("        (%d NO_DATA/negative row(s) left unscaled)" % st["kept_flag"])


def planned_files(plan):
    """[(path, text)] for every file the forcing stage writes."""
    return [(s["dst"], s["text"]) for s in plan["stations"]] + [(plan["dst_sdf"], plan["sdf_text"])]


def forcing_status(plan):
    """[(file name, state)] without writing anything: will be written / identical / DIFFERENT."""
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
                    "older rain file with other contents. Move or rename it, or re-run with "
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


def all_forcing_md5s(ctx):
    out = {}
    for plan in ctx.forcing_plans.values():
        out.update(forcing_md5s(plan))
    return out


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


# ----------------------------------------------------------------------
# Scoring a kept outlet file here (the project scorer's reading and formulas)
# ----------------------------------------------------------------------
_GRID_CACHE = {}


def grid5(path):
    """An outlet file on the 5-minute grid exactly as the project scorer builds it: time in hours
    from 2014-08-01, resampled to 5 minutes with time interpolation."""
    key = str(path)
    if key not in _GRID_CACHE:
        t, q, _, _ = read_qout(path)
        idx = pd.DatetimeIndex(pd.to_datetime(t * 3600, unit="s", origin=pd.Timestamp(ORIGIN)))
        _GRID_CACHE[key] = pd.Series(q, index=idx).resample("5min").interpolate(method="time")
    return _GRID_CACHE[key]


def window_compare(sim_path, ref_path, ev_start, ev_end):
    """Score sim against ref over the event window: PBIAS (percent), KGE_2012 and its parts, the
    volume and peak ratios (sim / ref) and the peak-time shift in minutes (sim minus ref)."""
    s = grid5(sim_path).loc[ev_start:ev_end]
    o = grid5(ref_path).loc[ev_start:ev_end]
    d = pd.concat([o.rename("o"), s.rename("s")], axis=1).dropna()
    if len(d) < 10:
        raise ValueError("the two outlet files share only %d points in the event window" % len(d))
    ov, sv = d["o"].to_numpy(float), d["s"].to_numpy(float)
    with np.errstate(all="ignore"):
        r = float(np.corrcoef(ov, sv)[0, 1])
        beta = float(sv.mean() / ov.mean())
        gamma = float((sv.std() / sv.mean()) / (ov.std() / ov.mean()))
        kge = float(1.0 - math.sqrt((r - 1.0) ** 2 + (beta - 1.0) ** 2 + (gamma - 1.0) ** 2))
        pbias = float(100.0 * (sv.sum() - ov.sum()) / ov.sum())
    t_o, t_s = d.index[int(np.argmax(ov))], d.index[int(np.argmax(sv))]
    return {"n": int(len(d)), "pbias": pbias, "kge": kge, "r": r, "beta": beta, "gamma": gamma,
            "vol_ratio": float(sv.sum() / ov.sum()), "peak_ratio": float(sv.max() / ov.max()),
            "peak_shift_min": float((t_s - t_o).total_seconds() / 60.0),
            "ref_peak": float(ov.max()), "sim_peak": float(sv.max())}


def expected_obs_peak(path, ev_start, ev_end):
    """The observed peak the scorer WILL see if it reads this truth file."""
    return float(grid5(path).loc[ev_start:ev_end].max())


def compare_series(fresh_path, truth_path):
    """Information only: is the new outlet file the same as the truth file, to the byte, and if not,
    how far apart are they over the WHOLE run (the scorer only looks at the storm)?"""
    out = {"fresh_md5": md5_of(fresh_path), "truth_md5": md5_of(truth_path)}
    out["byte_identical"] = out["fresh_md5"] == out["truth_md5"]
    try:
        tf, qf, _, _ = read_qout(fresh_path)
        tt, qt, _, _ = read_qout(truth_path)
        out["same_rows"] = bool(len(tf) == len(tt))
        if out["same_rows"] and len(tf):
            out["max_abs_time_diff_hr"] = float(np.max(np.abs(tf - tt)))
            out["max_abs_flow_diff_m3s"] = float(np.max(np.abs(qf - qt)))
            out["numerically_same"] = bool(out["max_abs_time_diff_hr"] <= 1e-9
                                           and out["max_abs_flow_diff_m3s"] <= SERIES_ABS_TOL_M3S)
        else:
            out["numerically_same"] = False
    except Exception as e:
        out["same_rows"] = None
        out["numerically_same"] = None
        out["compare_error"] = str(e)[:200]
    return out


def series_words(sc):
    if sc.get("byte_identical"):
        return "byte-identical to the truth file"
    if sc.get("numerically_same"):
        return "numerically the same as the truth file (largest flow difference %s m3/s)" \
               % _fmt(sc.get("max_abs_flow_diff_m3s"), "%.2g")
    if sc.get("same_rows") is False:
        return "NOT the same as the truth file (different number of rows)"
    return "NOT identical to the truth file (largest flow difference %s m3/s over the whole run)" \
           % _fmt(sc.get("max_abs_flow_diff_m3s"), "%.3g")


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
# Truths and the window
# ======================================================================
def find_truth_a(synth_dir):
    qout_files = list(synth_dir.glob("*.qout")) if synth_dir.exists() else []
    if len(qout_files) != 1:
        raise RuntimeError(
            "Expected exactly one *.qout file in %s (truth A), found %d: %s. Without it the scorer would silently "
            "fall back to real-gauge mode." % (synth_dir, len(qout_files), [f.name for f in qout_files]))
    return qout_files[0]


def check_truths(ctx):
    """Fingerprints and window for both truths. Sets ctx.truth_*_info. Returns a list of problems."""
    problems = []
    for tag, path, want_md5 in (("truth A", ctx.truth_a_path, TRUTH_A_MD5), ("truth B", ctx.truth_b_path, TRUTH_B_MD5)):
        if not path.exists():
            problems.append("%s not found: %s" % (tag, path))
            continue
        md5 = md5_of(path)
        if md5 != want_md5:
            problems.append("%s (%s) has md5 %s but Series 122 was scored against md5 %s. It is not the same file."
                            % (tag, path.name, md5, want_md5))
        info = describe_qout(path)
        if tag == "truth A":
            ctx.truth_a_md5, ctx.truth_a_info = md5, info
        else:
            ctx.truth_b_md5, ctx.truth_b_info = md5, info
        if not info["qout_ok"]:
            problems.append("%s fails its sanity check: %s" % (tag, info["qout_problem"]))
            continue
        if abs(info["qout_tlast_hr"] - ctx.runtime_hours) > QOUT_END_SLACK_HR:
            problems.append("%s ends at hour %.2f but the runs use a %g h window. Truth and runs must share one "
                            "window." % (tag, info["qout_tlast_hr"], ctx.runtime_hours))
        if abs(info["qout_dt_min"] - 3.75) > 0.05:
            problems.append("%s median time step is %.3f min, not the usual 3.75 min" % (tag, info["qout_dt_min"]))
    return problems


# ======================================================================
# One run: build, check the build, score, check the score, keep the outlet file, apply the gates
# ======================================================================
def in_file_value(text, keyword):
    """Value of a tRIBS .in keyword: the first number on the line after the line that starts with
    KEYWORD. Returns None if the layout is not what this expects."""
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


SCALED_NAME_RE = re.compile(r"_rain\d+p\d+")       # any scaled-rain file name this script writes


def check_config(ctx, cfg, run_id, run):
    """What the builder wrote into current_run_config.json (and the .in file) versus what was asked for."""
    problems = []
    if cfg.get("run_id") != run_id:
        problems.append("run_id is %r, expected %r" % (cfg.get("run_id"), run_id))
    if cfg.get("run_category") not in (None, CATEGORY):
        problems.append("run_category is %r, expected %r" % (cfg.get("run_category"), CATEGORY))
    wanted = [("Ks_mult", TRUTH_KS), ("f_RS_abs", TRUTH_F), ("optpercolation", run["optp"]),
              ("kinemvelcoef", TRUTH_ROUTING["kinemvelcoef"]), ("flowexp", TRUTH_ROUTING["flowexp"]),
              ("channelroughness", TRUTH_ROUTING["channelroughness"])]
    if run["kchan"] is not None:
        wanted.append(("channelconductivity_mmhr", run["kchan"]))
    for key, want in wanted:
        if key not in cfg:
            problems.append("%s missing from the run configuration" % key)
        elif not _close(cfg[key], want, 1e-9):
            problems.append("%s is %r, expected %r" % (key, cfg[key], want))
    if "runtime_hours" in cfg and not _close(cfg["runtime_hours"], ctx.runtime_hours, 1e-9):
        problems.append("runtime_hours is %r, expected %r" % (cfg["runtime_hours"], ctx.runtime_hours))
    if cfg.get("truth_file"):
        problems.append("the run configuration carries truth_file=%r, but this script scores every run against "
                        "the single truth in synth_truth/" % cfg.get("truth_file"))
    g = cfg.get("gauge_sdf")
    if run["kind"] == "rain":
        if g is not None and Path(str(g)).name != Path(run["sdf_rel"]).name:
            problems.append("gauge_sdf is %r, expected %r" % (g, run["sdf_rel"]))
    elif g and SCALED_NAME_RE.search(str(g)):
        problems.append("gauge_sdf is %r, a scaled-rain file, but this run is meant to use the recorded rain" % g)
    inp = cfg.get("input_file")
    if not inp:
        problems.append("no input_file in the run configuration")
        return problems
    p = Path(inp) if Path(inp).is_absolute() else Path.cwd() / inp
    try:
        text = p.read_text(errors="replace")
    except Exception as e:
        problems.append("cannot read the .in file %s (%s)" % (inp, e))
        return problems
    got_opt = in_file_value(text, "OPTPERCOLATION")
    got_cc = in_file_value(text, "CHANNELCONDUCTIVITY")
    if got_opt is None or (run["kchan"] is not None and got_cc is None):
        if not ctx.in_note_printed:
            print("  (note: could not read OPTPERCOLATION / CHANNELCONDUCTIVITY from %s in the layout I expected, "
                  "so that part of the .in-file check was skipped; the configuration check above still applies)"
                  % p.name)
            ctx.in_note_printed = True
    else:
        if int(round(got_opt)) != int(run["optp"]):
            problems.append("the .in file has OPTPERCOLATION %r, expected %d" % (got_opt, run["optp"]))
        if run["kchan"] is not None and not _close(got_cc, run["kchan"], 1e-3):
            problems.append("the .in file has CHANNELCONDUCTIVITY %r, expected %r (mm/hr)" % (got_cc, run["kchan"]))
    if run["kind"] == "rain":
        if Path(run["sdf_rel"]).name not in text:
            problems.append("the .in file does not name %s (the scaled rain would NOT be used)"
                            % Path(run["sdf_rel"]).name)
    else:
        if SCALED_NAME_RE.search(text):
            problems.append("the .in file names a scaled-rain file (%s), but this run is meant to use the "
                            "recorded rain" % SCALED_NAME_RE.search(text).group(0))
        elif Path(ctx.args.src_sdf).name not in text and not ctx.sdf_note_printed:
            print("  (note: the .in file does not mention %s; the recorded-rain check is skipped. The check that "
                  "no scaled-rain file is used still applies.)" % Path(ctx.args.src_sdf).name)
            ctx.sdf_note_printed = True
    return problems


def prepare_run(ctx, run):
    kwargs = {}
    if run["sdf_rel"]:
        kwargs["gauge_sdf"] = run["sdf_rel"]
    try:
        run_id, _, _ = builder.build_input_file("Ks_mult", TRUTH_KS, overrides=overrides_for(run),
                                                tag=run["tag"], **kwargs)
    except Exception as e:
        raise BuildFailed(str(e)[:300])
    if run_id != run["run_id"]:
        raise Abort("builder produced run_id %s, expected %s" % (run_id, run["run_id"]))
    cfg_path = ctx.calib_dir / "current_run_config.json"
    try:
        cfg = json.loads(cfg_path.read_text())
    except Exception as e:
        raise Abort("cannot read %s after the build: %s" % (cfg_path, e))
    problems = check_config(ctx, cfg, run_id, run)
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


def readback_problems(m, run):
    """Compare the scored metrics row against what this script asked for."""
    probs = []
    expected = {"Ks_mult": TRUTH_KS, "f_RS_abs": TRUTH_F, "optpercolation": run["optp"],
                "kinemvelcoef": TRUTH_ROUTING["kinemvelcoef"], "flowexp": TRUTH_ROUTING["flowexp"],
                "channelroughness": TRUTH_ROUTING["channelroughness"]}
    if run["kchan"] is not None:
        expected["channelconductivity_mmhr"] = run["kchan"]
    for key, want in expected.items():
        got = _num(m.get(key))
        if not np.isfinite(got) or not np.isclose(got, want, rtol=1e-6, atol=1e-12):
            probs.append("%s: scored row says %r, asked for %r" % (key, m.get(key), want))
    return probs


def keep_outlet(ctx, run_id):
    """Copy a finished run's outlet file into the kept folder (verified by checksum) and check it.
    Returns (status, reason, md5, qinfo, path). status OK or a failure code."""
    raw_dir = ctx.calib_dir / "02_results" / CATEGORY / run_id
    raw_qout = locate_raw_qout(raw_dir, run_id)
    if raw_qout is None:
        return "NO_QOUT", "no single *.qout found in %s" % raw_dir, None, None, None
    dest = ctx.kept_dir / ("%s_Outlet.qout" % run_id)
    try:
        md5 = copy_verified(raw_qout, dest)
    except Exception as e:
        return "COPY_FAILED", str(e)[:300], None, None, None
    qinfo = describe_qout(dest)
    if not qinfo["qout_ok"]:
        _set_aside_rejected(ctx, dest)
        return "QOUT_BAD", "outlet file failed its sanity check: %s" % qinfo["qout_problem"], md5, qinfo, None
    if qinfo["qout_tlast_hr"] < ctx.runtime_hours - QOUT_END_SLACK_HR:
        _set_aside_rejected(ctx, dest)
        return ("QOUT_SHORT", "outlet file ends at hour %.2f, but the window is %g h (a run that stopped early?)"
                % (qinfo["qout_tlast_hr"], ctx.runtime_hours), md5, qinfo, None)
    return "OK", "", md5, qinfo, dest


def _set_aside_rejected(ctx, dest):
    """A copy that failed its checks must not sit among the good outlet files, where a later analysis
    could pick it up. Move it into <kept folder>/_rejected/ (nothing is deleted; the original also
    stays in its raw results folder)."""
    try:
        rej = ctx.kept_dir / "_rejected"
        rej.mkdir(exist_ok=True)
        os.replace(dest, rej / dest.name)
    except Exception as e:
        print("  (could not move the rejected outlet file aside: %s -- delete %s by hand before analysing)"
              % (e, dest.name))


def check_disk(ctx):
    free_gb = shutil.disk_usage(ctx.calib_dir).free / 1e9
    if free_gb < MIN_FREE_GB:
        raise Abort("only %.1f GB of disk is free (need at least %.1f GB for one run). Free some space and re-run "
                    "the same command with --skip_existing." % (free_gb, MIN_FREE_GB))


def gate_run(ctx, run, cfg, metrics, kept):
    """Re-score the kept outlet file here, cross-check it against the project scorer, and apply the
    gate for this kind of run. Returns a dict: status, fail (a reason this run failed, others go on),
    abort (a reason to stop the whole script), extra (columns for the results row), line (text to print)."""
    res = {"status": "SUCCESS", "fail": None, "abort": None, "extra": {}, "line": "", "control_check": ""}
    ev0, ev1 = cfg["event_start"], cfg["event_end"]
    try:
        la = window_compare(kept, ctx.truth_a_path, ev0, ev1)
        lb = window_compare(kept, ctx.truth_b_path, ev0, ev1)
    except Exception as e:
        res.update(status="KEPT_MISMATCH", fail="could not re-score the kept outlet file here: %s" % str(e)[:200])
        return res
    res["extra"] = {
        "vs_A_pbias_pct": la["pbias"], "vs_A_kge_2012": la["kge"], "vs_A_vol_ratio": la["vol_ratio"],
        "vs_A_peak_ratio": la["peak_ratio"], "vs_A_peak_shift_min": la["peak_shift_min"],
        "vs_B_pbias_pct": lb["pbias"], "vs_B_kge_2012": lb["kge"], "vs_B_vol_ratio": lb["vol_ratio"],
        "vs_B_peak_ratio": lb["peak_ratio"], "vs_B_peak_shift_min": lb["peak_shift_min"],
    }
    res["line"] = ("vs A: KGE_2012 %.4f  PBIAS %+.2f%%  volume x%.4f  peak x%.4f  peak shift %+.0f min"
                   % (la["kge"], la["pbias"], la["vol_ratio"], la["peak_ratio"], la["peak_shift_min"]))

    # 1. the kept file must be the run that was scored
    sp, sk = _num(metrics.get("pbias_pct")), _num(metrics.get("kge_2012"))
    if (not (np.isfinite(sp) and np.isfinite(sk)) or abs(la["pbias"] - sp) > LOCAL_PBIAS_TOL
            or abs(la["kge"] - sk) > LOCAL_KGE_TOL):
        res.update(status="KEPT_MISMATCH",
                   fail="re-scoring the kept outlet file here gives PBIAS %+.4f%% / KGE_2012 %.5f but the project "
                        "scorer reported %s / %s. The kept file is not the run that was scored (or the two readings "
                        "of the file differ)." % (la["pbias"], la["kge"], _fmt(sp, "%+.4f%%"), _fmt(sk, "%.5f")))
        return res

    kind = run["kind"]
    if kind == "control_noloss":
        sc = compare_series(kept, ctx.truth_a_path)
        res["extra"]["series_vs_truth"] = series_words(sc)
        ok = abs(sp) < CTL_PBIAS_TOL and sk > CTL_KGE_MIN
        res["line"] += "\n  CONTROL A: %s -- the new outlet file is %s" % ("PASS" if ok else "FAIL", series_words(sc))
        if not ok:
            res.update(status="CONTROL_FAILED",
                       abort="control A (channel loss OFF, truth parameters) does not reproduce truth A: PBIAS %+.4f%% "
                             "and KGE_2012 %.4f against |PBIAS| < %g%% and KGE_2012 > %g. Today's build is not the "
                             "model that made truth A, so the other runs could not be compared with it. Send me "
                             "this log." % (sp, sk, CTL_PBIAS_TOL, CTL_KGE_MIN))
        else:
            res["control_check"] = "PASS"
    elif kind == "control_loss":
        sc = compare_series(kept, ctx.truth_b_path)
        res["extra"]["series_vs_truth"] = series_words(sc)
        ok_b = abs(lb["pbias"]) < CTL_PBIAS_TOL and lb["kge"] > CTL_KGE_MIN
        ok_rec = (abs(sp - REC_B_VS_A["pbias_pct"]) < REC_PBIAS_TOL and abs(sk - REC_B_VS_A["kge_2012"]) < REC_KGE_TOL)
        res["line"] += ("\n  CONTROL B: vs truth B PBIAS %+.4f%% KGE_2012 %.4f -> %s;  vs truth A the scorer gives "
                        "PBIAS %+.5f%% / KGE_2012 %.5f, Series 110 recorded %+.5f%% / %.5f -> %s;  the new outlet "
                        "file is %s"
                        % (lb["pbias"], lb["kge"], "PASS" if ok_b else "FAIL", sp, sk, REC_B_VS_A["pbias_pct"],
                           REC_B_VS_A["kge_2012"], "PASS" if ok_rec else "FAIL", series_words(sc)))
        if not (ok_b and ok_rec):
            res.update(status="CONTROL_FAILED",
                       abort="control B (channel loss ON, Kchan %g) does not reproduce truth B (%s). Today's build "
                             "or the Kchan setting is not the one that made truth B, so the dose-response could not "
                             "be tied to it. Send me this log."
                             % (TRUTH_B_KCHAN, "against truth B itself" if not ok_b else
                                "against the Series 110 record of truth B scored against truth A"))
        else:
            res["control_check"] = "PASS"
    elif kind == "rain":
        up = run["rain"] > 1.0
        moved = (abs(la["vol_ratio"] - 1.0) > RAIN_MIN_CHANGE and abs(la["peak_ratio"] - 1.0) > RAIN_MIN_CHANGE
                 and ((la["vol_ratio"] > 1.0) == up) and ((la["peak_ratio"] > 1.0) == up))
        if not moved:
            res.update(status="RAIN_NOT_USED",
                       abort="with the rain multiplied by %g the volume and the peak should both move %s from truth "
                             "A's, but they are x%.4f and x%.4f of it. tRIBS apparently did NOT use the scaled rain "
                             "(check that the .in file names %s and that tRIBS reads it)."
                             % (run["rain"], "UP" if up else "DOWN", la["vol_ratio"], la["peak_ratio"],
                                Path(run["sdf_rel"]).name))
    elif kind == "kchan":
        if not (la["pbias"] < 0):
            res["line"] += "\n  WARNING: channel loss is ON but the volume did not fall (PBIAS %+.3f%% against truth A)." \
                           % la["pbias"]
    return res


def run_one(ctx, run):
    """Build, score, verify, keep and gate one run. Returns a dict with status in SUCCESS /
    BUILD_FAILED / HANG / FAILED / NO_QOUT / COPY_FAILED / QOUT_BAD / QOUT_SHORT / KEPT_MISMATCH /
    OBSERVED_MISMATCH / READBACK_MISMATCH / CONTROL_FAILED / RAIN_NOT_USED, plus the metrics row.
    Raises Abort for a problem that would affect every run."""
    out = {"run_id": run["run_id"], "status": None, "reason": "", "elapsed": 0.0, "metrics": None, "qinfo": None,
           "qout_md5": None, "extra": {}, "line": "", "control_check": ""}
    check_disk(ctx)
    try:
        run_id, cfg = prepare_run(ctx, run)
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

    problems = check_observed(metrics, cfg, stdout, ctx.truth_a_path)
    if problems:
        out["abort"] = ("scoring did not use truth A as 'observed' for %s: %s. Every later run would be scored the "
                        "same wrong way, so the script is stopping." % (run_id, "; ".join(problems)))
        out.update(status="OBSERVED_MISMATCH", reason="; ".join(problems))
        return out
    probs = readback_problems(metrics, run)
    if probs:
        out["abort"] = ("the scored row for %s does not match what this script asked for: %s. The override plumbing "
                        "is broken, so the remaining runs would be wasted." % (run_id, "; ".join(probs)))
        out.update(status="READBACK_MISMATCH", reason="; ".join(probs))
        return out

    st, why, md5, qinfo, kept = keep_outlet(ctx, run_id)
    out["qinfo"], out["qout_md5"] = qinfo, md5
    if st != "OK":
        out.update(status=st, reason=why)           # raw folder is left in place for a look
        return out

    gate = gate_run(ctx, run, cfg, metrics, kept)
    out["extra"], out["line"], out["control_check"] = gate["extra"], gate["line"], gate["control_check"]
    if gate["abort"]:
        out["abort"] = gate["abort"]
        out.update(status=gate["status"], reason=gate["abort"])
        return out
    if gate["fail"]:
        # the kept file does not match what was scored: take it out of the good folder
        _set_aside_rejected(ctx, kept)
        out.update(status=gate["status"], reason=gate["fail"])
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


# ======================================================================
# Design record, audit, report
# ======================================================================
def run_record(r):
    return {"name": r["name"], "kind": r["kind"], "optpercolation": r["optp"], "kchan": r["kchan"],
            "rain_scale": r["rain"], "run_id": r["run_id"], "gauge_sdf": r["sdf_rel"]}


def design_dict(ctx, runs):
    return {
        "series": SERIES, "category": CATEGORY, "script": Path(__file__).name,
        "script_md5": md5_of(Path(__file__)),
        "created_local": datetime.now().isoformat(timespec="seconds"),
        "what_this_is": "Series 122 follow-up: twelve one-thing-at-a-time runs (two controls, a Kchan dose-response "
                        "with loss ON, a rain-scale series with loss OFF); Ks, f, cv, flowexp and n at truth; every "
                        "run scored against truth A and its outlet file kept.",
        "fixed": {"Ks_mult": TRUTH_KS, "f_RS_abs": TRUTH_F, "routing": TRUTH_ROUTING,
                  "channelporosity": builder.BASELINE.get("channelporosity"),
                  "channelwidthcoeff": builder.BASELINE.get("channelwidthcoeff"),
                  "As_value": builder.BASELINE.get("As_value"), "Au_value": builder.BASELINE.get("Au_value"),
                  "optsnow": "0 (set inside the builder)",
                  "kchan_of_no_loss_runs": "builder default %s (not used when optpercolation is 0)"
                                           % builder.BASELINE.get("channelconductivity_mmhr")},
        "runtime_hours": ctx.runtime_hours,
        "truth_A": {"file": ctx.truth_a_path.name, "md5": ctx.truth_a_md5},
        "truth_B": {"file": ctx.truth_b_path.name, "md5": ctx.truth_b_md5},
        "runs": [run_record(r) for r in runs],
        "forcing_md5s": all_forcing_md5s(ctx),
        "builder_md5": md5_of(Path(inspect.getsourcefile(builder))),
        "scorer_md5": md5_of(Path.cwd() / SCORER_SCRIPT),
        "python": sys.version.split()[0], "pandas": pd.__version__, "numpy": np.__version__,
        "pytRIBS": pytribs_version(),
    }


def design_problems(saved, ctx, runs):
    """Why a saved design does not match what this command would do now (empty list = it matches)."""
    probs = []
    try:
        a = [(x["name"], x["run_id"], int(x["optpercolation"]), x["kchan"], x["rain_scale"]) for x in saved["runs"]]
        b = [(r["name"], r["run_id"], int(r["optp"]), r["kchan"], r["rain"]) for r in runs]
        same = len(a) == len(b) and all(
            x[0] == y[0] and x[1] == y[1] and x[2] == y[2] and _close(x[4], y[4], 1e-12)
            and ((x[3] is None and y[3] is None) or (x[3] is not None and y[3] is not None and _close(x[3], y[3], 1e-12)))
            for x, y in zip(a, b))
        if not same:
            probs.append("the list of runs (--kchan / --rain) differs from the saved design")
    except Exception as e:
        probs.append("the saved design cannot be compared (%s)" % e)
    if saved.get("truth_A", {}).get("md5") != ctx.truth_a_md5:
        probs.append("truth A has a different checksum from the one this set started with")
    if saved.get("truth_B", {}).get("md5") != ctx.truth_b_md5:
        probs.append("truth B has a different checksum from the one this set started with")
    if saved.get("runtime_hours") is not None and not _close(saved["runtime_hours"], ctx.runtime_hours, 1e-9):
        probs.append("the window was %g h; now it is %g h" % (saved["runtime_hours"], ctx.runtime_hours))
    now = all_forcing_md5s(ctx)
    for name, md5 in (saved.get("forcing_md5s") or {}).items():
        if name in now and now[name] != md5:
            probs.append("rain file %s has changed since the design was saved" % name)
    return probs


def audit_results(ctx, df, runs):
    """Check the results CSV against the kept outlet files. Returns a list of problems."""
    problems = []
    if df.empty:
        return ["the results file has no rows"]
    if df["run_id"].duplicated().any():
        problems.append("%d duplicated run_id(s)" % int(df["run_id"].duplicated().sum()))
    want = {r["run_id"]: r for r in runs}
    have = set(df["run_id"].astype(str))
    missing = [want[k]["name"] for k in want if k not in have]
    extra = sorted(have - set(want))
    if missing:
        problems.append("%d planned run(s) have no row yet: %s" % (len(missing), ", ".join(missing)))
    if extra:
        problems.append("%d row(s) are not in the plan: %s" % (len(extra), ", ".join(extra[:3])))
    n_bad_q = n_bad_c = 0
    for _, r in df.iterrows():
        run_id = str(r["run_id"])
        q = ctx.kept_dir / ("%s_Outlet.qout" % run_id)
        if not q.exists():
            n_bad_q += 1
        elif "qout_md5" in df.columns and md5_of(q) != str(r.get("qout_md5")):
            n_bad_q += 1
        if not (ctx.csv_dir / ("%s_compare_obs_sim.csv" % run_id)).exists():
            n_bad_c += 1
        if run_id in want and want[run_id]["kind"].startswith("control") and str(r.get("control_check")) != "PASS":
            problems.append("control %s has no PASS on record" % want[run_id]["name"])
    if n_bad_q:
        problems.append("%d row(s) whose kept outlet file is missing or whose checksum differs" % n_bad_q)
    if n_bad_c:
        problems.append("%d row(s) with no compare CSV" % n_bad_c)
    try:
        stray = [f.name for f in ctx.kept_dir.glob("*_Outlet.qout")
                 if f.name[:-len("_Outlet.qout")] not in have]
        if stray:
            problems.append("%d outlet file(s) in %s with no row in the results file (e.g. %s)"
                            % (len(stray), ctx.kept_dir.name, stray[0]))
    except Exception:
        pass
    for col in ("kge_2012", "pbias_pct"):
        if col in df.columns and not np.isfinite(pd.to_numeric(df[col], errors="coerce")).all():
            problems.append("non-finite values in %s" % col)
    return problems


def _loglog_slope(x, y):
    x, y = np.log(np.asarray(x, float)), np.log(np.asarray(y, float))
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 2:
        return np.nan
    return float(np.polyfit(x[ok], y[ok], 1)[0])


def report_final(ctx, df, runs, n_failed):
    print("\n" + "=" * 100)
    print("SERIES %s CLEAN RUNS (Ks %g, f %g, cv %g, flowexp %g, n %g at truth; every run scored against truth A)"
          % (SERIES, TRUTH_KS, TRUTH_F, TRUTH_ROUTING["kinemvelcoef"], TRUTH_ROUTING["flowexp"],
             TRUTH_ROUTING["channelroughness"]))
    print("=" * 100)
    print("  Scored runs: %d of %d planned; %d record(s) in the FAILED log." % (len(df), len(runs), n_failed))
    if df.empty:
        return
    order = {r["run_id"]: r["idx"] for r in runs}
    d = df.copy()
    d["_o"] = d["run_id"].map(order)
    d = d.sort_values("_o")
    print("\n  %-8s %-14s %8s %7s | %9s %9s | %8s %8s %9s | %9s"
          % ("run", "kind", "Kchan", "rain x", "KGE vs A", "PBIAS A %", "vol x", "peak x", "shift min", "KGE vs B"))
    for _, r in d.iterrows():
        kc = _num(r.get("kchan_asked"))
        print("  %-8s %-14s %8s %7s | %9s %9s | %8s %8s %9s | %9s"
              % (str(r.get("run_name")), str(r.get("run_kind")), ("%g" % kc) if np.isfinite(kc) else "off",
                 "%g" % _num(r.get("rain_scale"), 1.0), _fmt(r.get("vs_A_kge_2012"), "%.4f"),
                 _fmt(r.get("vs_A_pbias_pct"), "%+.2f"), _fmt(r.get("vs_A_vol_ratio"), "%.4f"),
                 _fmt(r.get("vs_A_peak_ratio"), "%.4f"), _fmt(r.get("vs_A_peak_shift_min"), "%+.0f"),
                 _fmt(r.get("vs_B_kge_2012"), "%.4f")))
    print("  (vol x and peak x are the run's storm-window volume and peak divided by truth A's; shift = the run's "
          "peak time minus truth A's.)")

    ctl = d[d["run_kind"].astype(str).str.startswith("control")]
    print("\n  CONTROLS: %s" % (", ".join("%s %s" % (r["run_name"], "PASS" if str(r.get("control_check")) == "PASS"
                                                       else "not passed") for _, r in ctl.iterrows()) or "none yet"))

    loss = d[d["run_kind"].isin(["kchan", "control_loss"])].copy()
    loss["k"] = pd.to_numeric(loss["kchan_asked"], errors="coerce")
    loss = loss.sort_values("k")
    if len(loss) >= 3:
        v = pd.to_numeric(loss["vs_A_vol_ratio"], errors="coerce").to_numpy(float)
        ok = bool(np.all(np.diff(v) < 0))
        print("  DOSE-RESPONSE (Kchan, %d points): the volume falls steadily as Kchan rises: %s"
              % (len(loss), "yes" if ok else "NO -- look at the table (a bumpy dose-response needs an explanation)"))
        pk = pd.to_numeric(loss["vs_A_peak_ratio"], errors="coerce").to_numpy(float)
        print("    over Kchan %g to %g mm/hr the volume goes from x%.3f to x%.3f and the peak from x%.3f to x%.3f"
              % (loss["k"].iloc[0], loss["k"].iloc[-1], v[0], v[-1], pk[0], pk[-1]))
    rain = d[d["run_kind"] == "rain"].copy()
    rain["s"] = pd.to_numeric(rain["rain_scale"], errors="coerce")
    rain = rain.sort_values("s")
    if len(rain) >= 3:
        v = pd.to_numeric(rain["vs_A_vol_ratio"], errors="coerce").to_numpy(float)
        ok = bool(np.all(np.diff(v) > 0))
        print("  RAIN SERIES (%d points): the volume rises steadily with the rain multiplier: %s"
              % (len(rain), "yes" if ok else "NO -- look at the table"))
    if len(rain) >= 2:
        ev = _loglog_slope(rain["s"], pd.to_numeric(rain["vs_A_vol_ratio"], errors="coerce"))
        ep = _loglog_slope(rain["s"], pd.to_numeric(rain["vs_A_peak_ratio"], errors="coerce"))
        print("    rain elasticity (percent change in the result per percent change in the rain): volume %.2f, "
              "peak %.2f" % (ev, ep))
    b = d[d["run_name"] == "ctlB"]
    if len(b):
        print("  For reference, truth B (Kchan 201, loss ON): volume x%.4f, peak x%.4f, peak shift %+.0f min."
              % (_num(b.iloc[0].get("vs_A_vol_ratio")), _num(b.iloc[0].get("vs_A_peak_ratio")),
                 _num(b.iloc[0].get("vs_A_peak_shift_min"))))
    print("\n  Next: python plot_series122_clean_runs.py --dry_run   (the comparison and figures; no new tRIBS runs)")


# ======================================================================
# The run loop
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


def result_row(ctx, run, res):
    m = dict(res["metrics"])
    qinfo = res["qinfo"]
    m.update({
        "series": int(SERIES), "run_name": run["name"], "run_kind": run["kind"], "run_idx": int(run["idx"]),
        "optpercolation_asked": int(run["optp"]),
        "kchan_asked": float(run["kchan"]) if run["kchan"] is not None else np.nan,
        "rain_scale": float(run["rain"]), "gauge_sdf": run["sdf_rel"] or "recorded rain (builder default)",
        "control_check": res.get("control_check") or "",
        "truth_A_file": ctx.truth_a_path.name, "truth_A_md5": ctx.truth_a_md5,
        "truth_B_file": ctx.truth_b_path.name, "truth_B_md5": ctx.truth_b_md5,
        "runtime_hours": ctx.runtime_hours,
        "qout_kept": "%s/%s_Outlet.qout" % (ctx.kept_dir.name, m["run_id"]),
        "qout_md5": res["qout_md5"],
        "qout_rows": qinfo["qout_rows"] if qinfo else np.nan,
        "qout_last_hr": qinfo["qout_tlast_hr"] if qinfo else np.nan,
        "swept_param": "clean_runs_%s" % SERIES,
    })
    m.update(res.get("extra") or {})
    return m


def failed_row(run, res):
    return {"run_id": run["run_id"], "status": res["status"], "reason": res["reason"],
            "elapsed_min": res["elapsed"] / 60.0, "run_name": run["name"], "run_kind": run["kind"],
            "kchan_asked": run["kchan"], "rain_scale": run["rain"]}


def stage_forcing(ctx):
    """Write (or reuse) the scaled rain files, check them against the originals, record them."""
    print("\n" + "=" * 92)
    print("FORCING: scaled rain files for the %d rain runs" % len(ctx.forcing_plans))
    print("=" * 92)
    for scale, plan in ctx.forcing_plans.items():
        status = write_forcing(plan, ctx.args.overwrite_forcing)
        for name, st in status.items():
            print("  %-44s %s" % (name, st))
        problems = verify_forcing(plan)
        if problems:
            raise Abort("the scaled rain files for x%g do not match the originals: %s" % (scale, "; ".join(problems)))
    print("  CHECKED against the originals: same rows, same Y/M/D/H, rain rate = original x the multiplier, each SDF "
          "points at its own scaled files.")
    prov = {
        "script": Path(__file__).name,
        "created_local": datetime.now().isoformat(timespec="seconds"),
        "what_this_is": "Rain gauge forcing for the Series 122 clean runs: the rainfall RATE column multiplied by a "
                        "constant, Y/M/D/H unchanged, NO_DATA flags unchanged, whole record scaled. Same method as "
                        "run_lhs_Ks_f_120_storms_tribs6.py.",
        "row_hours_assumed_for_totals": RAIN_DT_HR,
        "scales": {},
        "python": sys.version.split()[0], "pandas": pd.__version__, "numpy": np.__version__,
    }
    for scale, plan in ctx.forcing_plans.items():
        prov["scales"]["%g" % scale] = {
            "source_sdf": {"file": plan["src_sdf"].name, "md5": md5_of(plan["src_sdf"])},
            "scaled_sdf": {"file": plan["dst_sdf"].name, "md5": md5_of(plan["dst_sdf"])},
            "stations": [{
                "source": s["src"].name, "source_md5": md5_of(s["src"]),
                "scaled": s["dst"].name, "scaled_md5": md5_of(s["dst"]),
                "rows": s["stats"]["rows"],
                "peak_mm_hr": [s["stats"]["peak_before"], s["stats"]["peak_after"]],
                "total_mm": [s["stats"]["sum_before"] * RAIN_DT_HR, s["stats"]["sum_after"] * RAIN_DT_HR],
                "no_data_rows_unscaled": s["stats"]["kept_flag"],
            } for s in plan["stations"]],
        }
    path = ctx.calib_dir / FORCING_PROV_NAME
    if path.exists():
        try:
            old = json.loads(path.read_text())
            same = {k: v for k, v in (old.get("scales") or {}).items()} == prov["scales"]
        except Exception:
            same = False
        if same:
            print("  Rain provenance record unchanged: %s" % path.name)
            return
        new = path.with_name("%s.superseded_%s" % (path.name, datetime.now().strftime("%Y%m%d_%H%M%S")))
        os.replace(path, new)
        print("  (an earlier rain provenance record was kept as %s)" % new.name)
    write_json_atomic(path, prov)
    print("  Saved: %s" % path.name)


def run_all(ctx, runs, todo):
    args = ctx.args
    out_path, failed_path = ctx.out_path, ctx.failed_path
    existing_df = load_existing_results(out_path)
    results = existing_df.to_dict("records") if not existing_df.empty else []
    if results and not args.skip_existing:
        raise Abort("%s already holds %d result row(s). Re-run with --skip_existing to resume, or move that file "
                    "(and %s/) aside to start again. This script will not overwrite earlier results."
                    % (out_path.name, len(results), ctx.kept_dir.name))
    by_id = {str(r.get("run_id")): r for r in results}

    failed_log_rows = []
    if args.skip_existing and failed_path.exists():
        try:
            failed_log_rows = pd.read_csv(failed_path).to_dict("records")
            print("  Loaded %d earlier FAILED record(s) from %s" % (len(failed_log_rows), failed_path.name))
        except Exception as e:
            print("  Warning: could not read %s (%s); starting a new log." % (failed_path.name, e))

    print("\n" + "=" * 92)
    print("Series %s clean runs -- %d of %d runs, window %g h, timeout %ds, cleanup %s"
          % (SERIES, len(todo), len(runs), ctx.runtime_hours, args.timeout, "off" if args.no_cleanup else "on"))
    print("  Truth A: %s [md5 %s]    Truth B: %s [md5 %s]"
          % (ctx.truth_a_path.name, ctx.truth_a_md5[:12], ctx.truth_b_path.name, ctx.truth_b_md5[:12]))
    print("  Output: %s   Kept outlet files: %s/" % (out_path.name, ctx.kept_dir.name))
    print("=" * 92 + "\n", flush=True)

    def save():
        if results:
            df = pd.DataFrame(results)
            if "run_idx" in df.columns:
                df = df.sort_values("run_idx")
            write_csv_atomic(df, out_path)

    def controls_passed():
        need = [r["run_id"] for r in runs if r["kind"].startswith("control")]
        return all(str(by_id.get(k, {}).get("control_check")) == "PASS" for k in need)

    completed = skipped = hung = failed = cleaned = 0
    aborted = None
    t_start = time.time()

    for run in todo:
        run_id = run["run_id"]
        print("\n[%2d/%d]  %-6s %s" % (run["idx"] + 1, len(runs), run["name"], run_id), flush=True)
        if run["kind"] == "kchan":
            print("        channel loss ON, Kchan %g mm/hr; routing, Ks and f at truth" % run["kchan"])
        elif run["kind"] == "rain":
            print("        channel loss OFF, rain x%g (%s); routing, Ks and f at truth"
                  % (run["rain"], Path(run["sdf_rel"]).name))
        elif run["kind"] == "control_loss":
            print("        CONTROL B: channel loss ON, Kchan %g mm/hr, rain as recorded" % run["kchan"])
        else:
            print("        CONTROL A: channel loss OFF, rain as recorded")

        kept = ctx.kept_dir / ("%s_Outlet.qout" % run_id)
        csv_ok = (ctx.csv_dir / ("%s_compare_obs_sim.csv" % run_id)).exists()
        row = by_id.get(run_id)
        if args.skip_existing and row is not None and kept.exists() and csv_ok:
            md5_ok = md5_of(kept) == str(row.get("qout_md5"))
            ctl_ok = (not run["kind"].startswith("control")) or str(row.get("control_check")) == "PASS"
            if md5_ok and ctl_ok:
                print("  SKIP (already complete): %s" % run_id)
                skipped += 1
                continue
            print("  (a result row exists but %s: re-running)"
                  % ("the kept outlet file's checksum has changed" if not md5_ok else "the control did not pass"))
        elif args.skip_existing and row is not None:
            print("  (a result row exists but the kept outlet file or compare CSV is missing: re-running)")

        if not run["kind"].startswith("control") and not controls_passed():
            aborted = ("the two controls (ctlA, ctlB) have not both passed, so %s will not be run. Run the controls "
                       "first (without --limit, or with --limit 2) and fix any control that fails." % run["name"])
            print("\n" + "!" * 92)
            print("STOPPING (CONTROLS NOT PASSED)\n" + aborted)
            print("!" * 92)
            break

        res = run_one(ctx, run)
        status, reason = res["status"], res["reason"]

        if res.get("abort"):
            print("\n" + "!" * 92)
            print("STOPPING (%s)" % status)
            if res.get("line"):
                print("  " + res["line"])
            print(res["abort"])
            print("Results so far are saved.")
            print("!" * 92)
            failed_log_rows = upsert_failed(failed_log_rows, failed_row(run, res))
            write_csv_atomic(pd.DataFrame(failed_log_rows), failed_path)
            failed += 1
            aborted = res["abort"]
            break

        if status == "SUCCESS":
            row = result_row(ctx, run, res)
            results = [r for r in results if r.get("run_id") != run_id]
            results.append(row)
            by_id[run_id] = row
            completed += 1
            note_outcome(ctx, True)
            if any(r.get("run_id") == run_id for r in failed_log_rows):
                failed_log_rows = [r for r in failed_log_rows if r.get("run_id") != run_id]
                write_csv_atomic(pd.DataFrame(failed_log_rows, columns=FAILED_COLS), failed_path)
            print("  scorer: KGE_2012=%.4f  PBIAS=%+.3f%%   (outlet file kept: %d rows, md5 %s)"
                  % (_num(row.get("kge_2012")), _num(row.get("pbias_pct")), row["qout_rows"], str(res["qout_md5"])[:10]))
            print("  " + res["line"])
            if not args.no_cleanup and cleanup_raw_results(run_id, ctx.calib_dir):
                cleaned += 1
        else:
            print("  %s: %s  (%s)" % (status, run_id, reason))
            if status == "HANG":
                hung += 1
            else:
                failed += 1
            failed_log_rows = upsert_failed(failed_log_rows, failed_row(run, res))
            write_csv_atomic(pd.DataFrame(failed_log_rows), failed_path)
            if run["kind"].startswith("control"):
                aborted = ("control %s did not finish (%s: %s). The other runs only mean something once both "
                           "controls have passed, so the script is stopping. Send me this log."
                           % (run["name"], status, reason))
                print("\n" + "!" * 92)
                print("STOPPING (CONTROL DID NOT RUN)\n" + aborted)
                print("!" * 92)
                save()
                break
            note_outcome(ctx, False)

        save()
        attempted = completed + hung + failed
        remaining = len(todo) - attempted - skipped
        if attempted > 0:
            avg = (time.time() - t_start) / attempted
            eta_min = avg * max(remaining, 0) / 60
            print("  Run time: %.1f min  |  ETA: %.0f min remaining" % (res["elapsed"] / 60, eta_min), flush=True)

    print("\nClean runs %s: %d ran, %d skipped, %d hung, %d failed, %d raw folders cleaned"
          % ("STOPPED EARLY" if aborted else "complete", completed, skipped, hung, failed, cleaned))
    save()
    if results:
        final_df = pd.DataFrame(results)
        if "run_idx" in final_df.columns:
            final_df = final_df.sort_values("run_idx")
        print("Saved: %s  (%d rows)" % (out_path.name, len(final_df)))
        try:
            full = (not args.limit) and len(final_df) == len(runs)
            if full and not aborted:
                probs = audit_results(ctx, final_df, runs)
                print("AUDIT: %s" % ("PASS (every row has its kept outlet file with a matching checksum and a compare "
                                     "CSV, and both controls passed)" if not probs else "PROBLEMS: " + "; ".join(probs)))
            report_final(ctx, final_df, runs, len(failed_log_rows))
        except Exception as e:
            print("\n(report skipped: %s -- results are saved regardless)" % e)
    else:
        print("No results to save.")
    return 3 if aborted else 0


# ======================================================================
# Dry run
# ======================================================================
def print_dry_run(ctx, runs, truth_problems, forcing_problems, competing):
    n = len(runs)
    est_sec = SEC_PER_MODEL_HOUR * ctx.runtime_hours
    print("\n" + "=" * 100)
    print("DRY RUN -- Series %s clean runs. Nothing is written, built or run." % SERIES)
    print("=" * 100)
    for tag, path, got, want in (("A", ctx.truth_a_path, getattr(ctx, "truth_a_md5", None), TRUTH_A_MD5),
                                 ("B", ctx.truth_b_path, getattr(ctx, "truth_b_md5", None), TRUTH_B_MD5)):
        print("  Truth %s (read-only): %s" % (tag, path))
        if got is None:
            print("    FILE NOT FOUND")
        else:
            print("    md5 %s -> %s" % (got, "matches the Series 122 fingerprint" if got == want
                                        else "DIFFERENT from the Series 122 fingerprint %s" % want))
    for tag, q in (("A", getattr(ctx, "truth_a_info", None)), ("B", getattr(ctx, "truth_b_info", None))):
        if q:
            print("    truth %s: %s rows, step %.2f min, last time %.2f h, peak %.4g m3/s"
                  % (tag, q["qout_rows"], q["qout_dt_min"], q["qout_tlast_hr"], q["qout_peak_m3s"]))
    print("  Simulation window used: %g h (the builder's own value)" % ctx.runtime_hours)
    print("  Builder: %s" % inspect.getsourcefile(builder))
    print("    unit guard OK (no division of Kchan by 3.6e6); overrides, tag and gauge_sdf accepted; "
          "optsnow 0 (set inside the builder)")
    print("    held at the builder's values: channelporosity %s, channelwidthcoeff %s, As %s, Au %s; Kchan of the "
          "no-loss runs %s (not used when channel loss is off)"
          % (builder.BASELINE.get("channelporosity"), builder.BASELINE.get("channelwidthcoeff"),
             builder.BASELINE.get("As_value"), builder.BASELINE.get("Au_value"),
             builder.BASELINE.get("channelconductivity_mmhr")))

    print("\n  RAIN FILES (originals are never touched; only the rain-rate column is multiplied):")
    if forcing_problems:
        for p in forcing_problems:
            print("    PROBLEM: %s" % p)
    else:
        for scale, plan in ctx.forcing_plans.items():
            print_forcing_plan(plan)
            for name, state in forcing_status(plan):
                print("        %-40s %s" % (name, state))

    print("\n  THE RUNS (%d; the controls come first):" % n)
    print("    %-3s %-8s %-15s %-12s %-8s  run id" % ("#", "name", "kind", "channel loss", "rain x"))
    for r in runs:
        loss = ("ON %g" % r["kchan"]) if r["optp"] == 1 else "off"
        print("    %-3d %-8s %-15s %-12s %-8g  %s" % (r["idx"] + 1, r["name"], r["kind"], loss, r["rain"], r["run_id"]))
    print("    Ks %g, f %g, cv %g, flowexp %g, n %g in every run"
          % (TRUTH_KS, TRUTH_F, TRUTH_ROUTING["kinemvelcoef"], TRUTH_ROUTING["flowexp"],
             TRUTH_ROUTING["channelroughness"]))
    ids = [r["run_id"] for r in runs]
    print("    run ids: %d runs, %d different (longest %d characters)"
          % (n, len(set(ids)), max(len(x) for x in ids)))
    print("    state now: %s" % ctx.state_line)
    print("\n  GATES: ctlA must reproduce truth A and ctlB must reproduce truth B (|PBIAS| < %g%%, KGE_2012 > %g; "
          "ctlB also the Series 110 numbers vs truth A); a rain run must lower the volume and the peak."
          % (CTL_PBIAS_TOL, CTL_KGE_MIN))
    print("\n  TIME: about %.0f s per run at this window -> %d runs ~ %.0f minutes (a rough estimate from the "
          "Series 110 timings: 67 s at 450 h)." % (est_sec, n, n * est_sec / 60.0))
    print("  Disk: %.1f GB free (a run needs about 1 GB while it is running; each kept outlet file is under "
          "1 MB)." % (shutil.disk_usage(ctx.calib_dir).free / 1e9))
    print("  Kept outlet files go to %s/ (not git-ignored); the scaled rain files go to %s."
          % (ctx.kept_dir.name, Path(ctx.args.src_sdf).parent))
    if competing:
        print("\n  NOTE: another tRIBS build/run script appears to be running right now:")
        for pid, cmd in competing:
            print("    pid %d: %s" % (pid, cmd[:110]))
        print("  (A dry run is harmless, but the real run would refuse to start.)")
    print("\n--dry_run: nothing written, built or run.")
    problems = list(truth_problems) + list(forcing_problems)
    if problems:
        print("\nNOT READY: fix these first:")
        for p in problems:
            print("  - %s" % p)
    else:
        print("\nReady for the trial (both controls and the first rain run):  python -u %s --limit 3 2>&1 | tee series122_clean_trial.log"
              % Path(__file__).name)
    return 2 if problems else 0


# ======================================================================
# Main
# ======================================================================
class Ctx(object):
    pass


def main():
    parser = argparse.ArgumentParser(
        description="Series 122 follow-up: twelve one-thing-at-a-time tRIBS runs (two controls, a Kchan "
                    "dose-response, a rain-scale series), each scored against truth A, every outlet file kept. "
                    "tRIBS 6.0.0.")
    parser.add_argument("--dry_run", action="store_true",
                        help="Print the plan, the checks and the time estimate, then stop. Writes, builds and runs nothing")
    parser.add_argument("--audit", action="store_true",
                        help="Check the results file against the kept outlet files (no tRIBS)")
    parser.add_argument("--kchan", default=",".join("%g" % k for k in DEFAULT_KCHAN),
                        help="Kchan values (mm/hr) of the dose-response runs, comma-separated (default %s)"
                             % ",".join("%g" % k for k in DEFAULT_KCHAN))
    parser.add_argument("--rain", default=",".join("%g" % s for s in DEFAULT_RAIN),
                        help="Rain multipliers of the rain runs, comma-separated (default %s)"
                             % ",".join("%g" % s for s in DEFAULT_RAIN))
    parser.add_argument("--limit", type=int, default=None,
                        help="Only the first N runs (the two controls come first, then the rain runs; resume later with --skip_existing)")
    parser.add_argument("--skip_existing", action="store_true",
                        help="Resume: skip runs whose result row, kept outlet file and compare CSV all exist")
    parser.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT,
                        help="Per-run hard timeout in seconds (default %d)" % DEFAULT_TIMEOUT)
    parser.add_argument("--no_cleanup", action="store_true",
                        help="Keep raw per-run results folders instead of deleting them after a successful run")
    parser.add_argument("--ignore_running_check", action="store_true",
                        help="Skip the check for other tRIBS scripts running (only if you are sure none are)")
    parser.add_argument("--overwrite_forcing", action="store_true",
                        help="Replace existing scaled rain files that differ from what this script writes")
    parser.add_argument("--max_consecutive_failures", type=int, default=MAX_CONSEC_FAIL,
                        help="Stop if this many runs in a row fail (default %d)" % MAX_CONSEC_FAIL)
    parser.add_argument("--src_sdf", default=SRC_SDF_DEFAULT,
                        help="Original rain gauge SDF, relative to lab/ (default %s)" % SRC_SDF_DEFAULT)
    args = parser.parse_args()

    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be at least 1")
    if sum([args.dry_run, args.audit]) > 1:
        parser.error("use only one of --dry_run, --audit")
    try:
        kchans = parse_float_list(args.kchan, "kchan", *KCHAN_RANGE)
        rains = parse_float_list(args.rain, "rain", *RAIN_RANGE)
    except ValueError as e:
        parser.error(str(e))
    if any(abs(k - TRUTH_B_KCHAN) < 1e-9 for k in kchans):
        parser.error("--kchan contains %g, which is control B already; leave it out" % TRUTH_B_KCHAN)
    if any(abs(s - 1.0) < 1e-9 for s in rains):
        parser.error("--rain contains 1, which is control A (rain as recorded); leave it out")

    ctx = Ctx()
    ctx.args = args
    ctx.consec_fail = 0
    ctx.in_note_printed = False
    ctx.sdf_note_printed = False
    cwd = Path.cwd()
    ctx.calib_dir = cwd.parent / "calibration_work"
    ctx.summary_dir = ctx.calib_dir / "03_comparisons" / "summary_tables"
    ctx.csv_dir = ctx.calib_dir / "03_comparisons" / "csv_exports"
    ctx.kept_dir = ctx.calib_dir / KEPT_DIR_NAME
    ctx.out_path = ctx.summary_dir / RESULTS_NAME
    ctx.failed_path = ctx.summary_dir / FAILED_NAME
    ctx.design_path = ctx.summary_dir / DESIGN_NAME

    print("Series %s clean runs   python %s, pandas %s, numpy %s"
          % (SERIES, sys.version.split()[0], pd.__version__, np.__version__))

    # ------------------------------------------------------------------
    # HARD SAFETY CHECKS
    # ------------------------------------------------------------------
    if not (cwd / SCORER_SCRIPT).exists():
        raise RuntimeError("%s not found in %s. Run this script from the lab/ directory." % (SCORER_SCRIPT, cwd))
    if not ctx.calib_dir.exists():
        raise RuntimeError("%s not found. Run this script from the lab/ directory (the folder that sits next to "
                           "calibration_work/)." % ctx.calib_dir)
    ctx.runtime_hours = float(builder.RUNTIME_HOURS)
    builder_unit_guard()
    install_series_hooks()

    ctx.truth_a_path = find_truth_a(ctx.calib_dir / "synth_truth")
    ctx.truth_b_path = (ctx.calib_dir / TRUTH_B_DIR_NAME / TRUTH_B_NAME).resolve()
    synth_res = (ctx.calib_dir / "synth_truth").resolve()
    if ctx.truth_b_path.parent == synth_res or synth_res in ctx.truth_b_path.parents:
        raise RuntimeError("truth B (%s) is inside synth_truth/. It must stay outside that folder (which must hold "
                           "exactly one file, truth A)." % ctx.truth_b_path)
    truth_problems = check_truths(ctx)
    print("Truth A, read-only: %s  [md5 %s]" % (ctx.truth_a_path.name, getattr(ctx, "truth_a_md5", "?")[:12]))
    print("Truth B, read-only: %s  [md5 %s]" % (ctx.truth_b_path.name, getattr(ctx, "truth_b_md5", "?")[:12]))

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

    runs = plan_runs(kchans, rains, args.src_sdf)
    todo = runs[:args.limit] if args.limit else runs

    ctx.forcing_plans, forcing_problems = {}, []
    for s in rains:
        try:
            ctx.forcing_plans[s] = plan_forcing(args.src_sdf, "rain%s" % label(s), s)
        except (FileNotFoundError, ValueError) as e:
            forcing_problems.append(str(e))
            break

    done_rows = 0
    if ctx.out_path.exists():
        try:
            done_rows = len(pd.read_csv(ctx.out_path))
        except Exception:
            done_rows = -1
    ctx.state_line = ("no results file yet" if not ctx.out_path.exists()
                      else "%s holds %s row(s)%s" % (ctx.out_path.name,
                                                    "an unreadable number of" if done_rows < 0 else done_rows,
                                                    "" if args.skip_existing else " (a run needs --skip_existing to resume)"))

    if args.audit:
        df = load_existing_results(ctx.out_path)
        probs = audit_results(ctx, df, runs)
        print("\nAUDIT of %s: %s" % (ctx.out_path.name, "PASS" if not probs else "PROBLEMS"))
        for p in probs:
            print("  - %s" % p)
        return 0 if not probs else 1

    if args.dry_run:
        return print_dry_run(ctx, runs, truth_problems, forcing_problems, competing)

    if truth_problems:
        raise RuntimeError("the truth files are not as expected: " + "; ".join(truth_problems))
    if forcing_problems:
        raise RuntimeError("the rain files cannot be prepared: " + "; ".join(forcing_problems))

    ctx.summary_dir.mkdir(parents=True, exist_ok=True)
    ctx.csv_dir.mkdir(parents=True, exist_ok=True)
    ctx.kept_dir.mkdir(parents=True, exist_ok=True)

    try:
        stage_forcing(ctx)
        # The design on disk must be the design in memory (resume safety).
        if ctx.design_path.exists():
            try:
                saved = json.loads(ctx.design_path.read_text())
            except Exception as e:
                raise Abort("cannot read the saved design %s (%s). Move it aside if you are sure you want a new set."
                            % (ctx.design_path.name, e))
            probs = design_problems(saved, ctx, runs)
            if probs:
                raise Abort("this command does not match the saved design in %s: %s. A resumed set must use the same "
                            "design; use the same --kchan/--rain, or move the old results aside."
                            % (ctx.design_path.name, "; ".join(probs)))
        else:
            write_json_atomic(ctx.design_path, design_dict(ctx, runs))
            print("Saved the design record: %s" % ctx.design_path.name)
        return run_all(ctx, runs, todo)
    except Abort as e:
        print("\n" + "!" * 92)
        print("STOPPED: %s" % e)
        print("!" * 92)
        print("Everything finished before this point is kept. After fixing the cause, re-run the same command "
              "(with --skip_existing).")
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
