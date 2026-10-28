"""
verify_truth_B_122_tribs6.py
============================
Series 122 -- is "truth B" really what its name says?

WHAT THIS IS
------------
Series 122 asks two questions of one 400-run sweep: (A) can the sweep recover a true Kchan
when the truth HAS channel loss, and (B) does it falsely use Kchan when the truth has none.
Question B was scored against truth A (the no-loss truth in synth_truth/). Question A needs
truth B: the same model with channel loss ON at Kchan 201 mm/hr.

Truth B is not a new run. It is a file Series 110 made on 2 Oct 2026:

    calibration_work/synth_truth_location/SMF_20140812_110_cc201p0_loc_K7p0_f0p012_Outlet.qout

Series 110 once had a truth file whose origin turned out to be wrong, so this script checks
the file directly instead of trusting its name or its records. It does two things.

  1. FINGERPRINTS (no tRIBS). The file's md5 must agree in four places: the file as it is
     now, the value this script was written with (013097a4...), Series 110's provenance
     record (PROVENANCE_truth_location_110.json) and Series 110's table
     (truth_location_110.csv). The file must also be 450 h long on the 3.75-minute grid,
     the same window as the Series 122 sweep.

  2. REBUILD AND SCORE (one tRIBS run, about 1 minute). It builds the truth-B point on the
     CURRENT build:  Ks 7.0, f 0.012, cv 4.5, flowexp 0.24, n 0.026, optpercolation 1,
     Kchan 201 mm/hr (optsnow 0 inside the builder), runs tRIBS, and scores the new
     hydrograph against the truth-B file. A file that really is that model comes back with
     PBIAS 0.0000% and KGE_2012 1.0000.

          PASS  = |PBIAS| < 0.01%  AND  KGE_2012 > 0.9999   (the Series 110 rule)

     It also checks that the scorer really used the truth-B file as "observed", that the
     parameters the scorer reports are the ones Series 110 recorded for this truth, and
     (information) whether the new outlet file is byte-identical to the truth-B file.

WHAT IT NEVER DOES
------------------
  * It never writes to, moves or renames anything in synth_truth/ or synth_truth_location/.
    Truth A and truth B are only READ.
  * It never touches the Series 122 results, the kept outlet files, or any Series 110 file.
  * It does not change the run configuration for anything that runs after it: the
    truth_file setting it adds to calibration_work/current_run_config.json is removed again
    at the end, even if the script stops part-way.

WHAT IT WRITES
--------------
  calibration_work/03_comparisons/summary_tables/truthB_provenance_check_122.json
      the result (PASS or FAIL) with all the numbers. Under the ignored 03_comparisons/
      folder, so committing it needs `git add -f`. An earlier record is renamed
      *.superseded_<time>, never overwritten.
  calibration_work/03_comparisons/csv_exports/<run_id>_compare_obs_sim.csv
  calibration_work/03_comparisons/summary_tables/<run_id>_metrics_summary.csv
      the scorer's usual outputs for the check run (regenerable).
  calibration_work/02_results/122_truthB_check/<run_id>/
      the raw tRIBS output of the check run. Deleted after a PASS; kept after a FAIL so
      it can be looked at.

  The check run is named  SMF_20140812_122_Ks7p0x_f0p012_truthBcheck_cc201p0

USAGE (run from the lab/ directory, like every other script)
------------------------------------------------------------
    cd /workspaces/tRIBS-Fall2026/workspaces/SMF/lab
    python verify_truth_B_122_tribs6.py --dry_run
    python -u verify_truth_B_122_tribs6.py 2>&1 | tee series122_truthB_check.log

Do not run any other tRIBS build/run script while it runs (they share
calibration_work/current_run_config.json with no locking).

EXIT CODES: 0 = PASS (or a clean dry run); 1 = FAIL / inconclusive / refused to start;
            2 = dry run found a problem; 3 = stopped by a safety check; 130 = Ctrl+C.

STATUS
------
Written without access to the Codespace. NOT YET RUN against the real tRIBS 6.0.0 build.

Tested against a stand-in lab (a toy hydrograph, a stub builder and a stub scorer; no real tRIBS):
  PASS path ........ byte-identical file verifies; PBIAS prints +0.0000%; old record renamed on re-run
  FAIL path ........ a file that is really Kchan 250 (records agree) -> NOT VERIFIED, raw output kept
  wrong file ....... tampered file, CSV disagreeing, window of 1416 h, builder default differing from
                     the record (porosity) -> all caught, before or after the run as designed
  scorer problems .. scorer ignoring truth_file -> STOPPED; run hanging (--timeout) and "tRIBS may
                     have failed" -> INCONCLUSIVE (says nothing about the truth file)
  safety ........... --truth_b inside synth_truth/, two files in synth_truth/, a planted 3.6e6 unit
                     bug, a real competing job -> all refused; Ctrl+C -> exit 130, tRIBS job killed
  clean-up ......... truth_file removed from current_run_config.json after every path above;
                     Truth A and truth B checksums unchanged after every path
  records .......... missing PROVENANCE json / truth_location csv -> warnings, check still runs

NOT tested (needs the real Codespace): the real builder's .in file, the real scorer's reported
parameter columns (any it does not report are skipped with a note, not failed), the real run time
(expected about 1 minute), and the real scorer honouring truth_file (the Series 120 storm script
uses the same mechanism; this script stops if the scorer did not use truth B as "observed").
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
# Identity
# ----------------------------------------------------------------------
SERIES        = "122"
CATEGORY      = "122_truthB_check"
SCORER_SCRIPT = "run_sensitivity_single_interp_tribs6.py"

# ----------------------------------------------------------------------
# Truth B: what it is claimed to be, and where Series 110 left it
# ----------------------------------------------------------------------
TRUTH_KS     = 7.0
TRUTH_F      = 0.012
TRUTH_KCHAN  = 201.0
TRUTH_ROUTING = {"kinemvelcoef": 4.5, "flowexp": 0.24, "channelroughness": 0.026}

TRUTH_B_DIR_NAME   = "synth_truth_location"
TRUTH_B_NAME       = "SMF_20140812_110_cc201p0_loc_K7p0_f0p012_Outlet.qout"
TRUTH_B_RUN_ID_110 = "SMF_20140812_110_cc201p0_loc_K7p0_f0p012"
# The md5 Tracy's own `md5sum` printed for this file on 26 Oct 2026 (it also appears in the
# Series 110 table and provenance record). Override with --expect_md5 only on purpose.
EXPECTED_MD5       = "013097a414def7b59aca731f9665f632"
PROV_JSON_NAME     = "PROVENANCE_truth_location_110.json"
TRUTH_CSV_NAME     = "truth_location_110.csv"
OUT_JSON_NAME      = "truthB_provenance_check_122.json"

# Parameters Series 110 recorded in truth_location_110.csv for every truth. The first seven
# are the ones this script sets on purpose; the rest come from the builder's defaults.
SET_BY_SCRIPT = ("Ks_mult", "f_RS_abs", "optpercolation", "channelconductivity_mmhr",
                 "kinemvelcoef", "flowexp", "channelroughness")
RECORDED_PARAM_COLS = SET_BY_SCRIPT + ("As_value", "Au_value", "thetaS_mult",
                                       "channelporosity", "channelwidthcoeff")

# ----------------------------------------------------------------------
# Pass rule (the Series 110 provenance rule) and other limits
# ----------------------------------------------------------------------
VERIFY_PBIAS_TOL    = 0.01      # percent
VERIFY_KGE_MIN      = 0.9999
DEFAULT_TIMEOUT     = 300
QOUT_END_SLACK_HR   = 2.0       # an outlet file must reach to within this many hours of the window end
OBS_PEAK_RTOL       = 0.02
EXPECTED_QOUT_COLS  = 3
SERIES_ABS_TOL_M3S  = 1e-6      # "numerically the same" for the full-series comparison (information)

COMPETING_PATTERNS = (
    "run_sensitivity_single", "run_sensitivity_sweep", "build_sensitivity_run",
    "run_cc_", "run_lhs_", "run_lowf_sweep", "run_crossversion", "verify_truth_provenance",
    "validate_truth_point", "rescore_series100", "run_truth_", "verify_truth_B_",
)


class Abort(RuntimeError):
    """A problem that makes the check meaningless: stop."""


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


def label(v):
    """The builder's own number-to-text rule (7.0 -> 7p0, 0.012 -> 0p012)."""
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
# Teaching the (unedited) builder about series 122 and this category
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


def check_run_tag():
    return "f%s_truthBcheck_cc%s" % (label(TRUTH_F), label(TRUTH_KCHAN))


def predicted_run_id():
    run_id, _ = builder.build_run_id("Ks_mult", TRUTH_KS)
    return "%s_%s" % (run_id, check_run_tag())


def overrides_now():
    """Everything this script sets for the check run (the rest comes from the builder's defaults)."""
    return {
        "Ks_mult":                  TRUTH_KS,
        "f_RS_abs":                 TRUTH_F,
        "kinemvelcoef":             TRUTH_ROUTING["kinemvelcoef"],
        "flowexp":                  TRUTH_ROUTING["flowexp"],
        "channelroughness":         TRUTH_ROUTING["channelroughness"],
        "optpercolation":           1,
        "channelconductivity_mmhr": TRUTH_KCHAN,
    }


def values_rebuild_will_use():
    """{column: value or None}: what the rebuild uses for each parameter Series 110 recorded."""
    ov = overrides_now()
    out = {}
    for col in RECORDED_PARAM_COLS:
        out[col] = ov[col] if col in ov else builder.BASELINE.get(col)
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
            "qout_tlast_hr": np.nan, "qout_peak_m3s": np.nan, "qout_tpeak_hr": np.nan,
            "qout_volume_m3": np.nan, "qout_problem": ""}
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
        info["qout_tpeak_hr"] = float(t[int(np.argmax(flow))])
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


def cleanup_raw_results(run_id, calib_dir):
    """Delete the check run's raw results folder. Only deletes a folder located EXACTLY at
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


def compare_series(fresh_path, truth_path):
    """Information only: is the new outlet file the same as the truth-B file, to the byte, and
    if not, how far apart are they over the WHOLE run (the scorer only looks at the storm)?"""
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
        print("\nInterrupted: the running tRIBS job was stopped.")
        raise
    elapsed = time.time() - t0
    if stdout:
        print(stdout, end="" if stdout.endswith("\n") else "\n")
    tribs_warning = bool(stdout) and "WARNING: tRIBS may have failed" in stdout
    returncode = None if timed_out else proc.returncode
    return returncode, elapsed, timed_out, tribs_warning, stdout or ""


# ======================================================================
# Facts about truth B (no tRIBS needed)
# ======================================================================
def gather_truth_facts(ctx):
    """Everything that can be known about truth B without running tRIBS. Returns a dict with
    'problems' (stop: the check cannot go on or would mean nothing), 'warnings' and the facts."""
    f = {"problems": [], "warnings": [], "notes": []}
    p = ctx.truth_path
    f["file"] = p.name
    if not p.exists():
        f["problems"].append("truth B file not found: %s" % p)
        return f
    f["bytes"] = int(p.stat().st_size)
    f["md5_live"] = md5_of(p)
    f["qinfo"] = describe_qout(p)
    q = f["qinfo"]
    if not q["qout_ok"]:
        f["problems"].append("truth B fails its sanity check: %s" % q["qout_problem"])

    # --- the four-way fingerprint ---
    sources = [("this script (EXPECTED_MD5 / --expect_md5)", ctx.expect_md5)]
    # Series 110 provenance record
    prov_path = p.parent / PROV_JSON_NAME
    f["prov_json"] = prov_path.name
    f["prov_record_md5"] = None
    f["truth_a_md5_in_record"] = None
    if prov_path.exists():
        try:
            prov = json.loads(prov_path.read_text())
            rec = (prov.get("files") or {}).get(p.name) or {}
            f["prov_record_md5"] = rec.get("md5")
            f["truth_a_md5_in_record"] = (prov.get("reference_truth_ccoff") or {}).get("md5")
            f["truth_a_name_in_record"] = (prov.get("reference_truth_ccoff") or {}).get("file")
            f["prov_routing_truth"] = prov.get("routing_truth")
            f["prov_optpercolation"] = prov.get("optpercolation")
            f["prov_created"] = prov.get("created_local")
            if f["prov_record_md5"] is None:
                f["problems"].append("%s has no record for %s" % (prov_path.name, p.name))
            else:
                sources.append((prov_path.name, f["prov_record_md5"]))
        except Exception as e:
            f["problems"].append("cannot read %s (%s)" % (prov_path.name, e))
    else:
        f["warnings"].append("%s not found next to the truth file; that cross-check was skipped" % prov_path.name)

    # Series 110 table
    f["csv_row"] = None
    csv_path = ctx.summary_dir / TRUTH_CSV_NAME
    f["csv_name"] = csv_path.name
    if csv_path.exists():
        try:
            df = pd.read_csv(csv_path)
            rows = df[df["run_id"].astype(str) == TRUTH_B_RUN_ID_110]
            if len(rows) == 0:
                f["problems"].append("%s has no row for %s" % (csv_path.name, TRUTH_B_RUN_ID_110))
            else:
                if len(rows) > 1:
                    f["warnings"].append("%s has %d rows for %s; the first is used"
                                         % (csv_path.name, len(rows), TRUTH_B_RUN_ID_110))
                row = rows.iloc[0].to_dict()
                f["csv_row"] = row
                if str(row.get("candidate_qout")) != p.name:
                    f["problems"].append("%s says this run's file is %r, not %r"
                                         % (csv_path.name, row.get("candidate_qout"), p.name))
                sources.append((csv_path.name + " (candidate_md5)", str(row.get("candidate_md5"))))
        except Exception as e:
            f["problems"].append("cannot read %s (%s)" % (csv_path.name, e))
    else:
        f["warnings"].append("%s not found in summary_tables/; that cross-check was skipped" % csv_path.name)

    f["fingerprints"] = [("the file as it is now", f["md5_live"])] + sources
    for name, value in sources:
        if value != f["md5_live"]:
            f["problems"].append("FINGERPRINT MISMATCH: the file's md5 is %s but %s says %s. The file is not "
                                 "the one that was recorded." % (f["md5_live"][:12], name, str(value)[:12]))

    # --- window ---
    if q["qout_ok"]:
        if abs(q["qout_tlast_hr"] - ctx.runtime_hours) > QOUT_END_SLACK_HR:
            f["problems"].append(
                "truth B ends at hour %.2f but the builder will run %g hours. The check, the truth and the "
                "Series 122 sweep must share one window. Bring this back to me before going further."
                % (q["qout_tlast_hr"], ctx.runtime_hours))
        if abs(q["qout_dt_min"] - 3.75) > 0.05:
            f["problems"].append("truth B's median time step is %.3f min, not the usual 3.75 min" % q["qout_dt_min"])

    # --- recorded parameters versus what the rebuild will use ---
    f["recorded_vs_rebuild"] = []
    row = f["csv_row"]
    now = values_rebuild_will_use()
    if row is not None:
        for col in RECORDED_PARAM_COLS:
            rec = row.get(col)
            use = now.get(col)
            if rec is None or (isinstance(rec, float) and np.isnan(rec)):
                f["recorded_vs_rebuild"].append((col, None, use, None))
                continue
            if use is None:
                f["recorded_vs_rebuild"].append((col, rec, None, None))
                f["warnings"].append("the builder has no default for %s, so the rebuild cannot be compared with "
                                     "the recorded value %r" % (col, rec))
                continue
            ok = _close(rec, use)
            f["recorded_vs_rebuild"].append((col, rec, use, ok))
            if not ok:
                if col in SET_BY_SCRIPT:
                    f["problems"].append(
                        "the Series 110 table records %s = %r for this truth, but this script rebuilds with %r. "
                        "The script's constants do not describe this file." % (col, rec, use))
                else:
                    f["warnings"].append(
                        "the Series 110 table records %s = %r for this truth, but the builder's default is now %r. "
                        "The rebuild will not be the same model, so expect the check to FAIL." % (col, rec, use))
        # Information: what truth B looked like against the no-loss truth
        for k in ("pbias_pct", "kge_2012", "vol_ratio_vs_ccoff", "peak_ratio_vs_ccoff"):
            f["vs_truth_a_" + k] = _num(row.get(k))
    return f


def truth_a_facts(ctx):
    """Truth A (the one file in synth_truth/): read only."""
    a = {"file": ctx.truth_a_path.name, "md5": md5_of(ctx.truth_a_path)}
    a["qinfo"] = describe_qout(ctx.truth_a_path)
    return a


# ======================================================================
# Printing
# ======================================================================
def print_truth_facts(ctx, f, a):
    q = f.get("qinfo", {})
    print("  Truth B (read-only): %s" % ctx.truth_path)
    if "md5_live" in f:
        print("    %s bytes, %s rows, step %.2f min, last time %.2f h, peak %.4g m3/s at %.2f h"
              % (format(f["bytes"], ","), q.get("qout_rows"), q.get("qout_dt_min", float("nan")),
                 q.get("qout_tlast_hr", float("nan")), q.get("qout_peak_m3s", float("nan")),
                 q.get("qout_tpeak_hr", float("nan"))))
        print("\n  FINGERPRINT (md5) in four places -- all must agree:")
        for name, value in f.get("fingerprints", []):
            mark = "ok" if value == f["md5_live"] else "DIFFERENT"
            print("    %-46s %s  %s" % (name, str(value)[:32], mark))
    print("\n  WINDOW: truth B ends at %.2f h; the builder will run %g h -> %s"
          % (q.get("qout_tlast_hr", float("nan")), ctx.runtime_hours,
             "OK" if abs(q.get("qout_tlast_hr", float("nan")) - ctx.runtime_hours) <= QOUT_END_SLACK_HR else "MISMATCH"))
    if f.get("recorded_vs_rebuild"):
        print("\n  PARAMETERS: what Series 110 recorded for this truth, against what the rebuild will use now:")
        print("    %-26s %-14s %-14s" % ("parameter", "recorded", "rebuild uses"))
        for col, rec, use, ok in f["recorded_vs_rebuild"]:
            tag = "" if ok is None else ("ok" if ok else "DIFFERENT")
            print("    %-26s %-14s %-14s %s"
                  % (col, "-" if rec is None else repr(rec), "-" if use is None else repr(use), tag))
    if f.get("prov_created"):
        print("    (Series 110 made this truth on %s; optpercolation %s, routing %s)"
              % (f["prov_created"], f.get("prov_optpercolation"), f.get("prov_routing_truth")))
    if f.get("vs_truth_a_pbias_pct") is not None and np.isfinite(f.get("vs_truth_a_pbias_pct", np.nan)):
        print("    Series 110 recorded truth B against the no-loss truth as: PBIAS %+.2f%%, volume ratio %.3f, "
              "peak ratio %.3f"
              % (f["vs_truth_a_pbias_pct"], f.get("vs_truth_a_vol_ratio_vs_ccoff", float("nan")),
                 f.get("vs_truth_a_peak_ratio_vs_ccoff", float("nan"))))
    print("\n  Truth A (read-only, never touched): %s  [md5 %s]" % (a["file"], a["md5"][:12]))
    rec_a = f.get("truth_a_md5_in_record")
    if rec_a:
        same = rec_a == a["md5"]
        print("    Series 110 made truth B against a no-loss reference with md5 %s -> %s"
              % (rec_a[:12], "the same file as truth A" if same else "NOT the same file as truth A"))
        if not same:
            f["warnings"].append("the no-loss reference recorded with truth B (md5 %s) is not the file now in "
                                 "synth_truth/ (md5 %s)" % (rec_a[:12], a["md5"][:12]))
    for w in f["warnings"]:
        print("  WARNING: %s" % w)


# ======================================================================
# The check run
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


def check_config(ctx, cfg, run_id):
    """What the builder wrote into current_run_config.json versus what was asked for."""
    problems = []
    if cfg.get("run_id") != run_id:
        problems.append("run_id is %r, expected %r" % (cfg.get("run_id"), run_id))
    if cfg.get("run_category") not in (None, CATEGORY):
        problems.append("run_category is %r, expected %r" % (cfg.get("run_category"), CATEGORY))
    for key, want in overrides_now().items():
        if key not in cfg:
            problems.append("%s missing from the run configuration" % key)
        elif not _close(cfg[key], want, 1e-9):
            problems.append("%s is %r, expected %r" % (key, cfg[key], want))
    if "runtime_hours" in cfg and not _close(cfg["runtime_hours"], ctx.runtime_hours, 1e-9):
        problems.append("runtime_hours is %r, expected %r" % (cfg["runtime_hours"], ctx.runtime_hours))
    if cfg.get("truth_file"):
        problems.append("the run configuration already carries truth_file=%r before this script set it"
                        % cfg.get("truth_file"))
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
                print("  (note: could not read OPTPERCOLATION / CHANNELCONDUCTIVITY from %s in the layout I "
                      "expected, so the .in-file check was skipped; the configuration check above still applies)"
                      % p.name)
            else:
                if int(round(got_opt)) != 1:
                    problems.append("the .in file has OPTPERCOLATION %r, expected 1" % got_opt)
                if not _close(got_cc, TRUTH_KCHAN, 1e-3):
                    problems.append("the .in file has CHANNELCONDUCTIVITY %r, expected %r (mm/hr)"
                                    % (got_cc, TRUTH_KCHAN))
        except Exception as e:
            problems.append("cannot read the .in file %s (%s)" % (inp, e))
    return problems


def check_observed(metrics, cfg, stdout, truth_path):
    """Did the scorer really use truth B as 'observed'? Returns a list of problems."""
    problems = []
    if str(metrics.get("obs_mode")) != "synth":
        problems.append("obs_mode is %r, expected 'synth' (the scorer fell back to the real gauge)"
                        % metrics.get("obs_mode"))
    m = re.search(r"Reading observed from:\s*(.+)", stdout or "")
    if m and m.group(1).strip() != truth_path.name:
        problems.append("the scorer read %r as observed, expected %r (it may have used truth A, the one file "
                        "in synth_truth/, instead of truth B)" % (m.group(1).strip(), truth_path.name))
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


def scored_params_vs_record(metrics, csv_row):
    """Compare the parameters the scorer reports for the check run with the ones Series 110 recorded
    for truth B. Returns (mismatches, compared, not_reported)."""
    mism, compared, missing = [], [], []
    if csv_row is None:
        return mism, compared, list(RECORDED_PARAM_COLS)
    for col in RECORDED_PARAM_COLS:
        rec = csv_row.get(col)
        got = metrics.get(col)
        if rec is None or (isinstance(rec, float) and np.isnan(rec)):
            continue
        if got is None or (isinstance(got, float) and np.isnan(got)):
            missing.append(col)
            continue
        compared.append(col)
        if not _close(rec, got):
            mism.append("%s: scorer reports %r, Series 110 recorded %r" % (col, got, rec))
    return mism, compared, missing


def set_truth_file_in_config(cfg_path, truth_rel):
    cfg = json.loads(cfg_path.read_text())
    cfg["truth_file"] = truth_rel
    write_json_atomic(cfg_path, cfg)
    back = json.loads(cfg_path.read_text())
    if back.get("truth_file") != truth_rel:
        raise Abort("could not write truth_file into %s" % cfg_path)


def clear_truth_file_in_config(cfg_path, truth_rel):
    """Leave the shared run configuration as we found it: remove the truth_file we added (only
    that one key, only if it is still ours). Never raises."""
    try:
        if not cfg_path.exists():
            return
        cfg = json.loads(cfg_path.read_text())
        if cfg.get("truth_file") == truth_rel:
            del cfg["truth_file"]
            write_json_atomic(cfg_path, cfg)
            print("  (removed the temporary truth_file setting from %s)" % cfg_path.name)
    except Exception as e:
        print("  WARNING: could not remove truth_file from %s (%s). Open it and delete the \"truth_file\" line "
              "by hand before running anything else." % (cfg_path.name, e))


def run_check(ctx, facts, truth_a):
    run_id = predicted_run_id()
    cfg_path = ctx.calib_dir / "current_run_config.json"
    print("\n" + "=" * 92)
    print("REBUILD AND SCORE: the truth-B point on today's build, scored against %s" % ctx.truth_path.name)
    print("  Expect PBIAS +0.0000%% and KGE_2012 1.0000 (need |PBIAS| < %g%% and KGE_2012 > %g)."
          % (VERIFY_PBIAS_TOL, VERIFY_KGE_MIN))
    print("=" * 92)
    print("  run id: %s" % run_id, flush=True)

    # Old outputs of an earlier attempt must never be mistaken for this one.
    cleanup_raw_results(run_id, ctx.calib_dir)
    (ctx.csv_dir / ("%s_compare_obs_sim.csv" % run_id)).unlink(missing_ok=True)
    (ctx.summary_dir / ("%s_metrics_summary.csv" % run_id)).unlink(missing_ok=True)

    try:
        built_id, _, _ = builder.build_input_file("Ks_mult", TRUTH_KS, overrides=overrides_now(), tag=check_run_tag())
    except Exception as e:
        raise Abort("the builder failed: %s" % str(e)[:300])
    if built_id != run_id:
        raise Abort("builder produced run_id %s, expected %s" % (built_id, run_id))
    try:
        cfg = json.loads(cfg_path.read_text())
    except Exception as e:
        raise Abort("cannot read %s after the build: %s" % (cfg_path, e))
    problems = check_config(ctx, cfg, run_id)
    if problems:
        raise Abort("the builder's run configuration does not match what was asked for: " + "; ".join(problems))

    set_truth_file_in_config(cfg_path, ctx.truth_rel)
    try:
        returncode, elapsed, timed_out, tribs_warning, stdout = run_scorer(ctx.args.timeout)
    finally:
        clear_truth_file_in_config(cfg_path, ctx.truth_rel)

    if timed_out:
        print("\nINCONCLUSIVE: the check run hung (no result in %d s) and was stopped. This says nothing about the "
              "truth file. Send me this log." % ctx.args.timeout)
        return 1, None
    if returncode != 0 or tribs_warning:
        print("\nINCONCLUSIVE: the check run did not finish (%s). This says nothing about the truth file. "
              "Send me this log." % ("tRIBS reported a failure" if tribs_warning else "scorer exited %s" % returncode))
        return 1, None
    metrics_file = ctx.summary_dir / ("%s_metrics_summary.csv" % run_id)
    if not metrics_file.exists():
        print("\nINCONCLUSIVE: no metrics file %s was written. Send me this log." % metrics_file.name)
        return 1, None
    try:
        metrics = pd.read_csv(metrics_file).iloc[0].to_dict()
    except Exception as e:
        print("\nINCONCLUSIVE: the metrics file could not be read (%s). Send me this log." % e)
        return 1, None

    obs_problems = check_observed(metrics, cfg, stdout, ctx.truth_path)
    if obs_problems:
        raise Abort("scoring did not use truth B as 'observed': %s. A PASS or FAIL from this run would mean "
                    "nothing." % "; ".join(obs_problems))

    pb, kg = _num(metrics.get("pbias_pct")), _num(metrics.get("kge_2012"))
    gate_ok = bool(np.isfinite(pb) and np.isfinite(kg) and abs(pb) < VERIFY_PBIAS_TOL and kg > VERIFY_KGE_MIN)
    mism, compared, not_reported = scored_params_vs_record(metrics, facts.get("csv_row"))

    raw_dir = ctx.calib_dir / "02_results" / CATEGORY / run_id
    raw_qout = locate_raw_qout(raw_dir, run_id)
    series = compare_series(raw_qout, ctx.truth_path) if raw_qout is not None else {"error": "no raw outlet file found"}

    reasons = []
    if not gate_ok:
        reasons.append("PBIAS %+.4f%% and KGE_2012 %.4f do not meet |PBIAS| < %g%% and KGE_2012 > %g"
                       % (_pb(pb), kg, VERIFY_PBIAS_TOL, VERIFY_KGE_MIN))
    if mism:
        reasons.append("the parameters of the rebuilt run differ from those Series 110 recorded: " + "; ".join(mism))
    passed = gate_ok and not mism

    print("\n" + "-" * 92)
    print("  obs_mode=%s   PBIAS=%+.4f%%   KGE_2012=%.4f   (kge_r %s, beta %s, gamma %s)"
          % (metrics.get("obs_mode"), _pb(pb), kg, _fmt(metrics.get("kge_r")), _fmt(metrics.get("kge_beta")),
             _fmt(metrics.get("kge_gamma"))))
    print("  the check run took %.2f min" % (elapsed / 60.0))
    print("  parameters: %d compared with the Series 110 record, %d differ%s"
          % (len(compared), len(mism), ("; not reported by the scorer: %s" % ", ".join(not_reported)) if not_reported else ""))
    if "byte_identical" in series:
        if series["byte_identical"]:
            print("  the new outlet file is BYTE-IDENTICAL to the truth-B file (md5 %s)" % series["fresh_md5"][:12])
        else:
            print("  the new outlet file is not byte-identical to the truth-B file: same rows %s, largest flow "
                  "difference %s m3/s over the whole run"
                  % (series.get("same_rows"), _fmt(series.get("max_abs_flow_diff_m3s"), "%.3g")))
    print("-" * 92)

    record = {
        "script": Path(__file__).name,
        "script_md5": md5_of(Path(__file__)),
        "when_local": datetime.now().isoformat(timespec="seconds"),
        "what_this_is": "Series 122 check that truth B (the Series 110 Ks 7 / f 0.012 / Kchan 201 centre-cell "
                        "truth) is what it is recorded to be: fingerprints agree, and a fresh rebuild on the "
                        "current build reproduces it.",
        "passed": passed,
        "pass_rule": {"abs_pbias_pct_below": VERIFY_PBIAS_TOL, "kge_2012_above": VERIFY_KGE_MIN},
        "truth_b": {"file": ctx.truth_path.name, "path": str(ctx.truth_path), "md5": facts["md5_live"],
                    "bytes": facts["bytes"], "rows": facts["qinfo"]["qout_rows"],
                    "last_hr": facts["qinfo"]["qout_tlast_hr"], "step_min": facts["qinfo"]["qout_dt_min"],
                    "peak_m3s": facts["qinfo"]["qout_peak_m3s"], "tpeak_hr": facts["qinfo"]["qout_tpeak_hr"],
                    "volume_m3": facts["qinfo"]["qout_volume_m3"]},
        "fingerprints": [{"source": n, "md5": v} for n, v in facts["fingerprints"]],
        "truth_a": {"file": truth_a["file"], "md5": truth_a["md5"],
                    "same_as_reference_in_series110_record": (facts.get("truth_a_md5_in_record") == truth_a["md5"])
                    if facts.get("truth_a_md5_in_record") else None},
        "rebuild": {"run_id": run_id, "set_by_script": overrides_now(),
                    "runtime_hours": ctx.runtime_hours, "optsnow": "0 (set inside the builder)",
                    "metrics": {k: metrics.get(k) for k in ("obs_mode", "pbias_pct", "kge_2012", "kge_r",
                                                            "kge_alpha", "kge_beta", "kge_gamma", "obs_peak_m3s",
                                                            "sim_peak_m3s", "peak_error_pct", "nse", "rmse_m3s")},
                    "elapsed_min": elapsed / 60.0},
        "parameters_compared_with_series110_record": compared,
        "parameters_differing": mism,
        "parameters_not_reported_by_scorer": not_reported,
        "full_series_comparison": series,
        "reasons_if_failed": reasons,
        "warnings": facts["warnings"],
        "builder_md5": md5_of(Path(inspect.getsourcefile(builder))),
        "scorer_md5": md5_of(Path.cwd() / SCORER_SCRIPT),
        "python": sys.version.split()[0], "pandas": pd.__version__, "numpy": np.__version__,
        "pytRIBS": pytribs_version(),
    }
    out = ctx.summary_dir / OUT_JSON_NAME
    if out.exists():
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        new = out.with_name("%s.superseded_%s" % (out.name, stamp))
        os.replace(out, new)
        print("  (an earlier record was kept as %s)" % new.name)
    write_json_atomic(out, record)
    print("  Saved: %s" % out.name)

    print("\n" + "=" * 92)
    if passed:
        print("TRUTH B: VERIFIED.  PBIAS %+.4f%%, KGE_2012 %.4f." % (_pb(pb), kg))
        print("  A fresh rebuild of Ks 7.0 / f 0.012 / cv 4.5 / flowexp 0.24 / n 0.026 / Kchan 201 (channel loss ON)")
        print("  on today's build reproduces %s." % ctx.truth_path.name)
        if series.get("byte_identical") is False and series.get("numerically_same") is False:
            print("  NOTE: the scorer only looks at the storm window, and the two files are not numerically the")
            print("  same over the whole run (largest flow difference %s m3/s). Tell me before relying on this."
                  % _fmt(series.get("max_abs_flow_diff_m3s"), "%.3g"))
        print("  It is safe to score the 400 kept Series 122 outlet files against this file (no new tRIBS runs).")
        if not ctx.args.no_cleanup:
            cleanup_raw_results(run_id, ctx.calib_dir)
        code = 0
    else:
        print("TRUTH B: NOT VERIFIED.")
        for r in reasons:
            print("  - %s" % r)
        print("  Do NOT use this file as truth B. The raw output of the check run was kept in")
        print("  %s for looking at. Send me this log." % raw_dir)
        code = 1
    print("=" * 92)
    return code, record


def _pb(x):
    """Percent bias for printing: a value that rounds to zero prints as +0.0000, not -0.0000."""
    return 0.0 if abs(x) < 5e-5 else x


def _fmt(x, spec="%.4f"):
    v = _num(x)
    return spec % v if np.isfinite(v) else "-"


# ======================================================================
# Dry run
# ======================================================================
def print_dry_run(ctx, facts, truth_a, competing):
    print("\n" + "=" * 92)
    print("DRY RUN -- Series 122 truth-B provenance check. Nothing is written, built or run.")
    print("=" * 92)
    print_truth_facts(ctx, facts, truth_a)
    print("\n  Builder: %s" % inspect.getsourcefile(builder))
    print("    unit guard OK (no division of Kchan by 3.6e6); overrides and tag accepted; optsnow 0 (set inside "
          "the builder)")
    print("    the check run will be named %s" % predicted_run_id())
    print("    it will be scored against %s via the scorer's truth_file setting (%s, relative to calibration_work)"
          % (ctx.truth_path.name, ctx.truth_rel))
    print("\n  WHAT THE CHECK WILL DO: build and run tRIBS once (about 1 minute) at Ks %g, f %g, cv %g, flowexp %g, "
          "n %g, optpercolation 1, Kchan %g mm/hr, window %g h,"
          % (TRUTH_KS, TRUTH_F, TRUTH_ROUTING["kinemvelcoef"], TRUTH_ROUTING["flowexp"],
             TRUTH_ROUTING["channelroughness"], TRUTH_KCHAN, ctx.runtime_hours))
    print("  and score it against truth B. PASS needs |PBIAS| < %g%% and KGE_2012 > %g, the scorer must report "
          "having used truth B as observed," % (VERIFY_PBIAS_TOL, VERIFY_KGE_MIN))
    print("  and the parameters the scorer reports must match the ones Series 110 recorded.")
    print("  It writes %s and the usual scorer outputs; nothing in synth_truth/ or synth_truth_location/ is touched."
          % OUT_JSON_NAME)
    state = ctx.summary_dir / OUT_JSON_NAME
    print("  state now: %s" % ("a record from an earlier check exists (%s); a new one will be saved beside it "
                               "(the old one renamed, not deleted)" % state.name if state.exists() else "no earlier check record"))
    if competing:
        print("\n  NOTE: another tRIBS build/run script appears to be running right now:")
        for pid, cmd in competing:
            print("    pid %d: %s" % (pid, cmd[:110]))
        print("  (A dry run is harmless, but the real check would refuse to start.)")
    print("\n--dry_run: nothing written, built or run.")
    if facts["problems"]:
        print("\nNOT READY: fix these first:")
        for p in facts["problems"]:
            print("  - %s" % p)
        return 2
    print("\nReady:  python -u %s 2>&1 | tee series122_truthB_check.log" % Path(__file__).name)
    return 0


# ======================================================================
# Main
# ======================================================================
class Ctx(object):
    pass


def main():
    parser = argparse.ArgumentParser(
        description="Series 122: check that truth B (Series 110, Ks 7 / f 0.012 / Kchan 201) is what it is "
                    "recorded to be, by fingerprints and by one fresh rebuild on tRIBS 6.0.0.")
    parser.add_argument("--dry_run", action="store_true",
                        help="Print the facts and the plan, then stop. Writes, builds and runs nothing")
    parser.add_argument("--truth_b", default=None,
                        help="Path to the truth-B .qout (default: calibration_work/%s/%s)"
                             % (TRUTH_B_DIR_NAME, TRUTH_B_NAME))
    parser.add_argument("--expect_md5", default=EXPECTED_MD5,
                        help="The md5 the truth-B file must have (default: the one recorded on 26 Oct 2026)")
    parser.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT,
                        help="Hard timeout for the check run in seconds (default %d)" % DEFAULT_TIMEOUT)
    parser.add_argument("--no_cleanup", action="store_true",
                        help="Keep the raw results folder of the check run even after a PASS")
    parser.add_argument("--ignore_running_check", action="store_true",
                        help="Skip the check for other tRIBS scripts running (only if you are sure none are)")
    args = parser.parse_args()

    ctx = Ctx()
    ctx.args = args
    ctx.expect_md5 = args.expect_md5.strip().lower()
    cwd = Path.cwd()
    ctx.calib_dir = cwd.parent / "calibration_work"
    ctx.summary_dir = ctx.calib_dir / "03_comparisons" / "summary_tables"
    ctx.csv_dir = ctx.calib_dir / "03_comparisons" / "csv_exports"
    ctx.synth_dir = ctx.calib_dir / "synth_truth"
    ctx.truth_path = Path(args.truth_b).resolve() if args.truth_b else \
        (ctx.calib_dir / TRUTH_B_DIR_NAME / TRUTH_B_NAME).resolve()

    print("Series 122 truth-B provenance check   python %s, pandas %s, numpy %s"
          % (sys.version.split()[0], pd.__version__, np.__version__))

    # ------------------------------------------------------------------
    # HARD SAFETY CHECKS
    # ------------------------------------------------------------------
    if not (cwd / SCORER_SCRIPT).exists():
        raise RuntimeError("%s not found in %s. Run this script from the lab/ directory (the folder that sits next "
                           "to calibration_work/)." % (SCORER_SCRIPT, cwd))
    if not ctx.calib_dir.exists():
        raise RuntimeError("%s not found. Run this script from the lab/ directory." % ctx.calib_dir)
    qouts = list(ctx.synth_dir.glob("*.qout")) if ctx.synth_dir.exists() else []
    if len(qouts) != 1:
        raise RuntimeError(
            "Expected exactly one *.qout file in %s (truth A), found %d: %s. Without it the scorer could silently "
            "fall back to real-gauge mode." % (ctx.synth_dir, len(qouts), [f.name for f in qouts]))
    ctx.truth_a_path = qouts[0]
    synth_res = ctx.synth_dir.resolve()
    if ctx.truth_path.parent == synth_res or synth_res in ctx.truth_path.parents:
        raise RuntimeError("the truth-B file %s is inside %s. Truth B must stay outside synth_truth/ (that folder "
                           "must hold exactly one file, truth A)." % (ctx.truth_path, ctx.synth_dir))
    try:
        ctx.truth_rel = os.path.relpath(ctx.truth_path, ctx.calib_dir.resolve()).replace(os.sep, "/")
    except ValueError:
        ctx.truth_rel = ""
    if not ctx.truth_rel or ctx.truth_rel.startswith(".."):
        raise RuntimeError("the truth-B file must sit inside calibration_work/ (the scorer's truth_file setting is "
                           "relative to that folder); it is at %s" % ctx.truth_path)

    install_series_hooks()
    builder_unit_guard()
    ctx.runtime_hours = float(getattr(builder, "RUNTIME_HOURS", 450))

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

    print("Truth A, read-only: %s  [md5 %s]" % (ctx.truth_a_path.name, md5_of(ctx.truth_a_path)[:12]))
    truth_a = truth_a_facts(ctx)
    facts = gather_truth_facts(ctx)

    if args.dry_run:
        return print_dry_run(ctx, facts, truth_a, competing)

    try:
        print()
        print_truth_facts(ctx, facts, truth_a)
        if facts["problems"]:
            raise Abort("the truth-B file cannot be checked yet: " + "; ".join(facts["problems"]))
        ctx.summary_dir.mkdir(parents=True, exist_ok=True)
        ctx.csv_dir.mkdir(parents=True, exist_ok=True)
        code, _ = run_check(ctx, facts, truth_a)
        return code
    except Abort as e:
        print("\n" + "!" * 92)
        print("STOPPED: %s" % e)
        print("!" * 92)
        print("Nothing in synth_truth/ or synth_truth_location/ was changed. After fixing the cause, run the same "
              "command again.")
        return 3
    except KeyboardInterrupt:
        clear_truth_file_in_config(ctx.calib_dir / "current_run_config.json", ctx.truth_rel)
        print("\nStopped by you. Nothing in synth_truth/ or synth_truth_location/ was changed.")
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
