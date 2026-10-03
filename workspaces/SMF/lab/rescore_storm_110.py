"""
rescore_storm_110.py
====================
Series 110 -- the SECOND-STORM test, part 2: RE-SCORE the stored storm runs
(250 cc-ON runs, and the 250 cc-OFF control runs if you made them) against every
cc-ON truth made by run_cc_storm_110.py under the same smaller (or larger) storm.
No tRIBS, no new simulations.

WHAT THIS DOES (plain language)
--------------------------------
run_cc_storm_110.py ran 250 joint (Ks, f, cc) simulations under the new storm and
saved every simulated hydrograph (the "Simulated" column of each run's
*_compare_obs_sim.csv). It also made 27 "answer keys" (truths): one tRIBS run for
each combination of true Ks, f and cc. While the sweep ran, every run was scored
against ONE reference (the cc-OFF truth for this storm). But a simulation does not
depend on which truth it is later compared to -- only its SCORE does. So each stored
hydrograph can be scored against any of the 27 truths in seconds.

This script keeps every run's stored "Simulated" series untouched, swaps in each
truth as "Observed" (read through the very same reading and resampling calls the
scorer uses), and recomputes the same metrics with the same formulas. The scoring and
summary functions are the ones from rescore_truth_location_110.py and
rescore_cc_truths_110.py (imported, not copied), so the original storm and the new
storm are scored by identical code and their tables have the same columns.

DRY TRUTHS
-----------
Under a smaller storm, a truth with a very high channel conductivity can lose all of
its water in the channel bed, so no flow reaches the outlet at all. A hydrograph of
zeros cannot be scored (KGE divides by the spread of the observed flow), so such a
truth is LEFT OUT and listed. That is a result in itself (it tells you where cc is so
large that the storm leaves no outlet flow), not an error. A truth that keeps less than
10% of the reference volume ("NEAR_DRY" in the table) is still scored, but treat it with
care: it is nearly an answer key of zeros.

SELF-CHECKS ("gates") -- run before anything is re-scored
-----------------------------------------------------------
  GATE A  This storm's cc-OFF reference (synth_truth_ccoff_<label>/), read by this
          script's code, must equal the "Observed" column stored in every compare CSV.
          (It proves every run was scored against this reference, and that this script
          reads a truth file the way the scorer did.)
  GATE B  Every metric recomputed from a compare CSV must equal the one stored in the
          sweep results CSV.
  GATE C  Each truth file must match the checksum recorded when it was made
          (PROVENANCE_truth_location_110.json in its folder), come from a run that passed
          its sanity check, lie inside the sampled box, and share its 5-minute time grid
          with the reference. Truths with no outlet flow are set aside as DRY (above).
  GATE D  Each truth file must equal the "Simulated" column of its OWN run's compare CSV.
          (Shows the saved truth file is the very simulation that was scored, and not a
          file from another run.) The difference must be below 1e-6 of the peak to count
          as exact; up to 1e-3 of the peak is noted and accepted (a small difference in how
          the two were read); anything larger stops the script.
If a gate fails the script stops and says why; nothing is re-scored.
(--drop_bad_runs excludes individual runs that fail A or B and continues. --ignore_selfcheck
continues past a gate D failure and records that in the provenance file -- only do that
if you understand why it failed.)
If GATE A fails for every run with a difference of a fraction of a m3/s, suspect a
different pandas version from the one that scored the runs (resampling a series with
.resample("5min").interpolate(method="time") gives slightly different values in pandas 2
and pandas 3). Use the same python you used for run_cc_storm_110.py; do not use
--drop_bad_runs or edit anything to get past it.

WHAT YOU GET   (calibration_work/03_comparisons/summary_tables/rescored_storm_110_<label>/)
------------------------------------------------------------------------------------------
  rescored_location_long_110_<label>.csv     every run x every truth: all metrics, plus the truth's
                                             own Ks, f and cc (truth_Ks, truth_f, truth_cc);
                                             same columns as the original storm's file, plus
                                             storm_label and rain_scale at the end
  rescored_location_summary_110_<label>.csv  one row per (truth, series): best KGE_2012, how many
                                             runs are feasible, where the near-best runs sit in cc,
                                             how close the nearest sample is to the truth point;
                                             plus storm_label, rain_scale, truth_flags
  PROVENANCE_rescored_storm_110_<label>.json checksums, versions, gate results, which truths were
                                             set aside and why

HOW TO READ THE TABLE (rules of thumb, mine)
----------------------------------------------
  - "ON best" is the best KGE_2012 any cc-ON run reached against that truth; "OFF best" is
    the best the cc-OFF runs reached (they can only move Ks and f). A small gap means Ks and f
    alone can imitate that truth. (No control sweep yet? Then the OFF columns show n/a.)
  - "near-best cc range" is the span of cc among cc-ON runs within --kge_tol of the best. A
    narrow range containing the true cc is what identifiable looks like.
  - "nearest" is the distance, in units of the sampled box (0 = same point, 1 = a full side),
    from the truth point to the closest of the 250 samples; "n<.15" counts samples within 0.15.
  - One storm, 250 points in 3 dimensions, routing pinned at truth, noise-free truth: this is
    a best case, not a realistic one. Feasible = |PBIAS| < 2%.

SAFETY
-------
  - Read-only on every existing file. Writes only inside its own output folder (refuses an
    output folder that is, or overlaps, any input location).
  - Never imports the builder, never touches current_run_config.json, never starts tRIBS.
  - Refuses to start while run_cc_storm_110.py is still running for the same storm (its results
    files are still changing). A runner for a DIFFERENT storm (for example the 1.25x one) does
    not touch this storm's files and is allowed.
  - Stops if the sweep or the truth set is incomplete (an unfinished run) -- re-run
    run_cc_storm_110.py with --skip_existing first (or use --allow_partial to look anyway).
  - Run membership comes only from the run_id rows in the results CSVs.

USAGE (run from the lab/ directory; keep rescore_cc_truths_110.py and
       rescore_truth_location_110.py in the same folder; use the SAME python as for the runs)
----------------------------------------------------------------------------------------------
    python rescore_storm_110.py --check_only        # gates only; writes nothing
    python rescore_storm_110.py                     # gates, then every usable truth (0.8x storm)
    python rescore_storm_110.py --rain_scale 1.25   # the third storm
    python rescore_storm_110.py --no_control        # ignore the control sweep even if it exists
    python rescore_storm_110.py --kge_tol 0.005
"""

import argparse
import json
import os
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import rescore_cc_truths_110 as RS
    import rescore_truth_location_110 as TL
except ImportError:
    sys.exit("This script needs rescore_cc_truths_110.py and rescore_truth_location_110.py in the same folder "
             "(it reuses their reading, gate and scoring code). Put them in lab/ and run from there.")

_NEEDED_TL = ["score_truth", "summarize", "truth_id", "inside_box", "report", "SER_ON", "SER_OFF"]
_missing = [n for n in _NEEDED_TL if not hasattr(TL, n)]
if _missing:
    sys.exit("rescore_truth_location_110.py in this folder is an older version (missing: %s). Replace it with "
             "the current copy." % ", ".join(_missing))
_NEEDED_RS = ["load_results", "load_all_compares", "load_compare", "gate_a", "gate_b", "read_truth_5min",
              "import_phase_fn", "md5_of", "DEFAULT_EXPECT_N", "DEFAULT_KGE_TOL", "FEASIBLE_PBIAS_PCT"]
_missing = [n for n in _NEEDED_RS if not hasattr(RS, n)]
if _missing:
    sys.exit("rescore_cc_truths_110.py in this folder is an older version (missing: %s). Replace it with the "
             "current copy." % ", ".join(_missing))

SER_ON, SER_OFF = RS.SER_ON, RS.SER_OFF
LOC_PROV = "PROVENANCE_truth_location_110.json"
LHS_SERIES = "110"
SELF_OK_REL = 1e-6       # gate D: a truth file equals its own run's Simulated column to this (relative to the peak)
SELF_FAIL_REL = 1e-3     # ... and is accepted with a note up to this; beyond it the script stops
DEFAULT_SCALE = 0.8      # run_cc_storm_110.py's default


# ------------------------------------------------------------------
# Small helpers
# ------------------------------------------------------------------
def storm_label(scale):
    """0.8 -> 'storm080' (the same rule run_cc_storm_110.py uses)."""
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("--rain_scale must be a positive number")
    if abs(scale - 1.0) < 1e-9:
        raise ValueError("--rain_scale 1.0 is the original storm; use rescore_truth_location_110.py for that.")
    pct = int(round(scale * 100))
    if abs(pct / 100.0 - scale) > 1e-9:
        raise ValueError("--rain_scale needs at most 2 decimal places (0.8, 0.85, 1.25), got %r" % scale)
    return "storm%03d" % pct


def truthy(v):
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    return str(v).strip().lower() in ("true", "1", "yes")


def write_csv_atomic(df, path):
    tmp = path.with_name(path.name + ".part")
    df.to_csv(tmp, index=False)
    os.replace(str(tmp), str(path))


def write_json_atomic(path, obj):
    tmp = path.with_name(path.name + ".part")
    tmp.write_text(json.dumps(obj, indent=2, default=str))
    os.replace(str(tmp), str(path))


def ancestor_pids():
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


def running_storm_runners():
    """Python processes running run_cc_storm_110.py: [(pid, command line, rain_scale)], or None if ps failed."""
    try:
        out = subprocess.run(["ps", "-eo", "pid=,args="], capture_output=True, text=True, timeout=10).stdout
    except Exception:
        return None
    mine = ancestor_pids()
    found = []
    for line in out.splitlines():
        parts = line.strip().split(None, 1)
        if len(parts) != 2 or not parts[0].isdigit():
            continue
        pid, args = int(parts[0]), parts[1]
        if pid in mine or not re.search(r"\bpython[\d.]*\b", args) or "run_cc_storm_110" not in args:
            continue
        m = re.search(r"--rain_scale[=\s]+([0-9]*\.?[0-9]+)", args)
        found.append((pid, args.strip(), float(m.group(1)) if m else DEFAULT_SCALE))
    return found


def check_storm_columns(df, name, label, scale):
    """If the file says which storm it belongs to, it must be this one (guards against a wrong path)."""
    if "storm_label" in df.columns:
        other = sorted(set(df["storm_label"].dropna().astype(str)) - {label})
        if other:
            sys.exit("%s says it belongs to %s, not %s. Wrong file for this storm?" % (name, ", ".join(other), label))
    if "rain_scale" in df.columns:
        vals = pd.to_numeric(df["rain_scale"], errors="coerce").dropna()
        if len(vals) and not np.allclose(vals.to_numpy(float), scale, rtol=0, atol=1e-9):
            sys.exit("%s has rain_scale values other than %g. Wrong file for this storm?" % (name, scale))


# ------------------------------------------------------------------
# Inputs
# ------------------------------------------------------------------
def load_reference(ref_dir, ref_csv, label, scale):
    """This storm's cc-OFF reference file (what every sweep run was scored against). Returns (path, info dict)."""
    if not ref_dir.exists():
        sys.exit("No folder %s. Run run_cc_storm_110.py first (and let it finish)." % ref_dir)
    info = {"recorded_md5": ""}
    name = None
    if ref_csv.exists():
        df = pd.read_csv(ref_csv)
        check_storm_columns(df, ref_csv.name, label, scale)
        if len(df) != 1 or "candidate_qout" not in df.columns:
            sys.exit("%s should hold exactly one reference row with a candidate_qout column; found %d row(s)."
                     % (ref_csv.name, len(df)))
        row = df.iloc[0]
        name = str(row["candidate_qout"])
        info["recorded_md5"] = "" if pd.isna(row.get("candidate_md5", "")) else str(row.get("candidate_md5", ""))
        problem = "" if pd.isna(row.get("qout_problem", "")) else str(row.get("qout_problem", "")).strip()
        if "qout_ok" in df.columns and not truthy(row["qout_ok"]):
            sys.exit("The %s reference failed its sanity check when it was made (%s). Re-run run_cc_storm_110.py "
                     "without --skip_existing to rebuild it." % (label, problem or "no reason recorded"))
    else:
        files = sorted(ref_dir.glob("*.qout"))
        if len(files) != 1:
            sys.exit("%s not found, and %s does not hold exactly one *.qout (found %d)." % (ref_csv.name, ref_dir.name, len(files)))
        name = files[0].name
        print("  NOTE: %s not found; using the only file in %s with no recorded checksum." % (ref_csv.name, ref_dir.name))
    path = ref_dir / name
    if not path.exists():
        sys.exit("The reference %s named in %s is missing from %s." % (name, ref_csv.name, ref_dir.name))
    info["md5"] = RS.md5_of(path)
    if info["recorded_md5"] and info["md5"] != info["recorded_md5"]:
        sys.exit("The reference %s does not match the checksum recorded when it was made. It was changed or "
                 "replaced after the runs were scored; refusing." % name)
    return path, info


def load_registry(loc_dir, label, scale):
    """The truths as recorded by run_cc_storm_110.py. Returns (registry list, unregistered file names, provenance dict)."""
    prov_path = loc_dir / LOC_PROV
    if not prov_path.exists():
        sys.exit("No %s in %s. Run run_cc_storm_110.py first (and let it finish)." % (LOC_PROV, loc_dir))
    try:
        prov = json.loads(prov_path.read_text())
    except Exception as e:
        sys.exit("Could not read %s: %s" % (prov_path, e))
    if prov.get("storm_label") not in (None, label):
        sys.exit("%s belongs to %s, not %s." % (prov_path.name, prov.get("storm_label"), label))
    files = prov.get("files", {})
    if not files:
        sys.exit("%s lists no truth files." % prov_path.name)
    reg = []
    for name, info in files.items():
        try:
            ks, f, cc = float(info["Ks_mult"]), float(info["f_RS_abs"]), float(info["cc_mmhr"])
        except (KeyError, TypeError, ValueError):
            sys.exit("%s: the entry for %s has no usable Ks_mult / f_RS_abs / cc_mmhr." % (prov_path.name, name))
        reg.append({"name": name, "Ks": ks, "f": f, "cc": cc, "md5": str(info.get("md5", "") or ""),
                    "flags": str(info.get("flags", "") or ""), "run_id": str(info.get("run_id", "") or ""),
                    "id": TL.truth_id(ks, f, cc)})
    ids = [r["id"] for r in reg]
    if len(set(ids)) != len(ids):
        sys.exit("%s lists two truths with the same (Ks, f, cc); refusing." % prov_path.name)
    on_disk = {p.name for p in loc_dir.glob("*.qout")}
    unregistered = sorted(on_disk - set(files))
    reg.sort(key=lambda r: (r["Ks"], r["f"], r["cc"]))
    return reg, unregistered, prov


def expected_truth_count(prov):
    grid = prov.get("grid", {})
    try:
        return int(len(grid["Ks_mult"]) * len(grid["f_RS_abs"]) * len(grid["cc_mmhr"]))
    except (KeyError, TypeError):
        return None


def selfcheck_truth(tr, s, csv_dir):
    """Gate D for one truth: its file (as read here) vs the Simulated column of its own run's compare CSV."""
    entry = {"truth_id": tr["id"], "run_id": tr["run_id"]}
    p = (csv_dir / ("%s_compare_obs_sim.csv" % tr["run_id"])) if tr["run_id"] else None
    if p is None or not p.exists():
        entry.update({"status": "not checked", "detail": "no compare CSV for this truth's own run"})
        return entry
    try:
        df = RS.load_compare(p)
    except Exception as e:
        entry.update({"status": "not checked", "detail": "its compare CSV is unreadable: %s" % str(e)[:80]})
        return entry
    ref = s.reindex(df.index).to_numpy(float)
    sim = df["Simulated"].to_numpy(float)
    if np.isnan(ref).any():
        entry.update({"status": "FAIL", "detail": "the truth file has no value at some timestamps of its own compare CSV"})
        return entry
    peak = float(np.nanmax(np.abs(sim))) or 1.0
    d = float(np.nanmax(np.abs(ref - sim)))
    rel = d / peak
    status = "ok" if rel <= SELF_OK_REL else ("minor" if rel <= SELF_FAIL_REL else "FAIL")
    entry.update({"status": status, "max_abs_diff_m3s": d, "rel_to_peak": rel,
                  "detail": "largest difference %.2e m3/s (%.1e of the peak)" % (d, rel)})
    return entry


# ------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Series 110 -- re-score the second/third storm's stored runs against every cc-ON truth of that storm "
                    "(no tRIBS).")
    ap.add_argument("--rain_scale", type=float, default=DEFAULT_SCALE,
                    help="Which storm: the rain multiplier given to run_cc_storm_110.py (default 0.8; 1.25 for the third storm)")
    ap.add_argument("--check_only", action="store_true", help="Run the gates and stop; writes nothing.")
    ap.add_argument("--kge_tol", type=float, default=RS.DEFAULT_KGE_TOL,
                    help="Near-best tolerance in KGE_2012 units (default %g)" % RS.DEFAULT_KGE_TOL)
    ap.add_argument("--expect_n", type=int, default=RS.DEFAULT_EXPECT_N,
                    help="Expected runs per series (default %d)" % RS.DEFAULT_EXPECT_N)
    ap.add_argument("--allow_partial", action="store_true",
                    help="Continue even if the sweep has fewer than --expect_n runs or fewer truths than the grid.")
    ap.add_argument("--no_control", action="store_true",
                    help="Ignore the control (cc OFF) sweep even if its file exists.")
    ap.add_argument("--drop_bad_runs", action="store_true",
                    help="Exclude runs that fail gate A or B and continue, instead of stopping.")
    ap.add_argument("--ignore_selfcheck", action="store_true",
                    help="Continue even if gate D fails. Recorded in the provenance file.")
    ap.add_argument("--ignore_running_check", action="store_true",
                    help="Skip the check for a run_cc_storm_110.py still running for this storm.")
    ap.add_argument("--sweep_csv", type=Path, default=None, help="cc-ON sweep results CSV (default: the storm's file)")
    ap.add_argument("--control_csv", type=Path, default=None, help="cc-OFF control results CSV (default: the storm's file)")
    ap.add_argument("--out_dir", type=Path, default=None)
    args = ap.parse_args()

    try:
        label = storm_label(args.rain_scale)
    except ValueError as e:
        ap.error(str(e))
    scale = float(args.rain_scale)

    script_dir = Path.cwd()
    calib_dir = script_dir.parent / "calibration_work"
    summary_dir = calib_dir / "03_comparisons" / "summary_tables"
    csv_dir = calib_dir / "03_comparisons" / "csv_exports"
    synth_dir = calib_dir / "synth_truth"
    ref_dir = calib_dir / ("synth_truth_ccoff_%s" % label)
    loc_dir = calib_dir / ("synth_truth_location_%s" % label)
    ref_csv = summary_dir / ("storm_ref_110_%s.csv" % label)
    truth_csv = summary_dir / ("truth_location_110_%s.csv" % label)
    forcing_prov = calib_dir / ("forcing_%s_PROVENANCE.json" % label)
    out_dir = args.out_dir or (summary_dir / ("rescored_storm_110_%s" % label))
    long_name = "rescored_location_long_110_%s.csv" % label
    summ_name = "rescored_location_summary_110_%s.csv" % label
    prov_name = "PROVENANCE_rescored_storm_110_%s.json" % label
    p_on = args.sweep_csv or (summary_dir / ("lhs_results_joint_Ks_f_cc_%s_%s.csv" % (label, LHS_SERIES)))
    p_off = args.control_csv or (summary_dir / ("lhs_results_joint_Ks_f_cc_CONTROL_%s_%s.csv" % (label, LHS_SERIES)))

    print("\n" + "=" * 78)
    print("Series 110 -- re-score %s (rain x%g) runs against the %s truths" % (label, scale, label))
    print("python %s, pandas %s, numpy %s" % (sys.version.split()[0], pd.__version__, np.__version__))
    print("=" * 78)

    # ---- safety: is the runner still writing these files? --------------------
    if args.ignore_running_check:
        print("Running-jobs check skipped (--ignore_running_check).")
    else:
        runners = running_storm_runners()
        if runners is None:
            print("WARNING: could not check whether run_cc_storm_110.py is still running (ps unavailable). Make sure it "
                  "has finished for this storm.")
        else:
            same = [r for r in runners if abs(r[2] - scale) < 1e-9]
            other = [r for r in runners if abs(r[2] - scale) >= 1e-9]
            if same:
                print("\nrun_cc_storm_110.py IS STILL RUNNING for this storm (x%g):" % scale)
                for pid, cmd, _ in same:
                    print("  pid %d: %s" % (pid, cmd[:110]))
                sys.exit("Refusing to start: its results files are still being written, and the set of runs and "
                         "truths is not final. Wait for it to finish (or --ignore_running_check if you are sure it is "
                         "a stale entry).")
            if other:
                print("NOTE: run_cc_storm_110.py is running for a different storm (x%g); it does not write this "
                      "storm's files, so this is allowed." % other[0][2])
            else:
                print("No run_cc_storm_110.py running.")

    # ---- inputs / output-folder safety ------------------------------------------
    ref_path, ref_info = load_reference(ref_dir, ref_csv, label, scale)
    print("%s cc-OFF reference (read-only): %s  [md5 %s%s]"
          % (label, ref_path.name, ref_info["md5"][:12], "" if ref_info["recorded_md5"] else ", no recorded checksum"))

    out_res = out_dir.resolve()
    for guarded in (synth_dir, ref_dir, loc_dir, csv_dir):
        gd = guarded.resolve()
        if out_res == gd or gd in out_res.parents or out_res in gd.parents:
            sys.exit("Output folder %s overlaps an input folder (%s); refusing." % (out_dir, guarded))
    for pth in (p_on, p_off, truth_csv, ref_csv):
        rp = pth.resolve()
        if out_res in rp.parents or out_res == rp:
            sys.exit("Output folder %s contains an input file (%s); refusing." % (out_dir, pth))

    registry, unregistered, loc_prov = load_registry(loc_dir, label, scale)
    print("Truth registry: %d truth(s) recorded in %s/%s" % (len(registry), loc_dir.name, LOC_PROV))
    if unregistered:
        print("  NOTE: %d .qout file(s) in %s are not in the registry and are ignored: %s"
              % (len(unregistered), loc_dir.name, ", ".join(unregistered[:4]) + (" ..." if len(unregistered) > 4 else "")))
    n_grid = expected_truth_count(loc_prov)
    if n_grid is not None and len(registry) < n_grid:
        msg = "only %d of the %d truths in the grid are recorded (the runner has not finished them, or some failed)" \
              % (len(registry), n_grid)
        if not args.allow_partial:
            sys.exit("STOPPED: %s. Re-run run_cc_storm_110.py with --skip_existing first (or --allow_partial to look anyway)."
                     % msg)
        print("  WARNING (--allow_partial): %s." % msg)

    if truth_csv.exists():
        tdf = pd.read_csv(truth_csv)
        check_storm_columns(tdf, truth_csv.name, label, scale)
    runner_pandas = loc_prov.get("pandas")
    if runner_pandas:
        if str(runner_pandas).split(".")[:2] != pd.__version__.split(".")[:2]:
            print("  *** NOTE: pandas here is %s but the truths were made with pandas %s. The 5-minute interpolation "
                  "differs slightly between pandas 2 and 3; use the same python as for the runs. Gate A will tell. ***"
                  % (pd.__version__, runner_pandas))
        else:
            print("  pandas %s matches the version that made the truths." % pd.__version__)

    print("\nLoading results tables:")
    if not p_on.exists():
        sys.exit("The cc-ON sweep results %s do not exist yet. Run run_cc_storm_110.py (the sweep stage) first." % p_on.name)
    res_on = RS.load_results(p_on, SER_ON, 1, args.expect_n)
    check_storm_columns(res_on, p_on.name, label, scale)
    use_control = (not args.no_control) and p_off.exists()
    res_off = None
    if use_control:
        res_off = RS.load_results(p_off, SER_OFF, 0, args.expect_n)
        check_storm_columns(res_off, p_off.name, label, scale)
    elif args.no_control:
        print("  Control sweep ignored (--no_control).")
    else:
        print("  No control (cc OFF) sweep for %s yet (%s): scoring the cc-ON runs only." % (label, p_off.name))
    for nm, rdf in ((SER_ON, res_on), (SER_OFF, res_off)):
        if rdf is not None and len(rdf) != args.expect_n:
            msg = "%s has %d runs, expected %d (unfinished or failed runs)" % (nm, len(rdf), args.expect_n)
            if not args.allow_partial:
                sys.exit("STOPPED: %s. Re-run run_cc_storm_110.py%s with --skip_existing to finish them (or "
                         "--allow_partial to look anyway; for an unfinished control, --no_control)."
                         % (msg, " --control" if nm == SER_OFF else ""))
            print("  WARNING (--allow_partial): %s." % msg)

    phase_fn = RS.import_phase_fn(script_dir)
    obs_ref = RS.read_truth_5min(ref_path)

    print("\nLoading compare CSVs (run membership = run_id rows above):")
    cache_on, prob_on = RS.load_all_compares(res_on, csv_dir)
    parts = [res_on.assign(_series=SER_ON)]
    cache, probs = dict(cache_on), dict(prob_on)
    msg = "  %s: %d/%d loaded" % (SER_ON, len(cache_on), len(res_on))
    if res_off is not None:
        cache_off, prob_off = RS.load_all_compares(res_off, csv_dir)
        cache.update(cache_off)
        probs.update(prob_off)
        parts.append(res_off.assign(_series=SER_OFF))
        msg += ";  %s: %d/%d loaded" % (SER_OFF, len(cache_off), len(res_off))
    print(msg)

    # ---- gates A and B -----------------------------------------------------------
    bad = dict(probs)
    allres = pd.concat(parts, ignore_index=True)
    badA, worstA = RS.gate_a(allres, cache, obs_ref)
    badB, worstB = RS.gate_b(allres, cache, phase_fn)
    for d in (badA, badB):
        for rid, why in d.items():
            bad.setdefault(rid, why)

    n_checked = len(cache)
    print("\nGATE A  observed series : %d/%d compare CSVs match this storm's cc-OFF reference as read here "
          "(worst abs diff %.2e m3/s)" % (n_checked - len(badA), n_checked, worstA))
    if worstB:
        wk = max(worstB, key=worstB.get)
        print("GATE B  stored metrics  : %d/%d runs reproduce their stored metrics (worst abs deviation %.2e in %s; %s)"
              % (n_checked - len(badB), n_checked, worstB[wk], wk,
                 "phase metrics included" if phase_fn else "phase metrics skipped"))
    else:
        print("GATE B  stored metrics  : %d/%d" % (n_checked - len(badB), n_checked))

    if n_checked and len(badA) >= 0.9 * n_checked:
        print("\n  GATE A failed for nearly every run. The usual cause is a different pandas version from the one that "
              "made the truths and scored the runs (here: pandas %s; the truths were made with pandas %s). Use the "
              "same python as for run_cc_storm_110.py. Do NOT use --drop_bad_runs for this."
              % (pd.__version__, runner_pandas or "unknown"))
    if bad:
        print("\n  %d run(s) failed a gate or could not be loaded. First few:" % len(bad))
        for rid, why in list(bad.items())[:8]:
            print("    %s: %s" % (rid, why))
        if not args.drop_bad_runs:
            sys.exit("\nSTOPPED: nothing was re-scored. Fix the cause (or re-run with --drop_bad_runs to exclude "
                     "these runs and continue).")
        print("  --drop_bad_runs set: these runs are excluded.")
    res_on = res_on[~res_on["run_id"].isin(bad)].reset_index(drop=True)
    if res_off is not None:
        res_off = res_off[~res_off["run_id"].isin(bad)].reset_index(drop=True)
    if res_on.empty or (res_off is not None and res_off.empty):
        sys.exit("No runs left in one of the series after the gates; stopping.")
    win_idx = cache[res_on["run_id"].iloc[0]].index

    # ---- gate C: the truths --------------------------------------------------------
    usable, dry, skipped, truth_obs = [], [], [], {}
    for tr in registry:
        p = loc_dir / tr["name"]
        lab = "Ks %g f %g cc %g" % (tr["Ks"], tr["f"], tr["cc"])
        if not p.exists():
            skipped.append((lab, "file missing")); continue
        if "BAD_QOUT" in tr["flags"]:
            skipped.append((lab, "the file failed its sanity check when it was made (BAD_QOUT); re-run "
                                 "run_cc_storm_110.py with --skip_existing to rebuild it")); continue
        if not TL.inside_box(tr["Ks"], tr["f"], tr["cc"]):
            skipped.append((lab, "outside the sampled box; a sweep that never looked there cannot recover it")); continue
        if not tr["md5"]:
            skipped.append((lab, "no recorded checksum")); continue
        if RS.md5_of(p) != tr["md5"]:
            skipped.append((lab, "checksum differs from the one recorded when it was made (GATE C)")); continue
        try:
            s = RS.read_truth_5min(p)
        except Exception as e:
            skipped.append((lab, "unreadable: %s" % str(e)[:80])); continue
        if not (s.index.equals(obs_ref.index) and (s.isna() == obs_ref.isna()).all()):
            skipped.append((lab, "5-minute grid differs from the reference's (GATE C)")); continue
        vals = s.reindex(win_idx).to_numpy(float)
        fin = vals[np.isfinite(vals)]
        if fin.size == 0 or float(fin.max()) <= 0.0:
            dry.append((lab, "DRY: no outlet flow inside the scoring window under this storm (left out; a hydrograph "
                             "of zeros cannot be scored)"))
            tr["_dry"] = True
            continue
        usable.append(tr)
        truth_obs[tr["id"]] = s
    n_other = len(skipped)
    print("GATE C  truths          : %d/%d scored (checksum, sanity flag, box, 5-minute grid);  %d dry (no outlet flow);  "
          "%d other problem(s)" % (len(usable), len(registry), len(dry), n_other))
    for name, why in dry:
        print("    set aside %s: %s" % (name, why))
    for name, why in skipped:
        print("    SKIPPED (needs attention) %s: %s" % (name, why))
    near_dry = [tr for tr in usable if "NEAR_DRY" in tr["flags"]]
    if near_dry:
        print("    note: %d scored truth(s) are NEAR_DRY (keep <10%% of the reference volume): %s"
              % (len(near_dry), ", ".join("Ks %g f %g cc %g" % (t["Ks"], t["f"], t["cc"]) for t in near_dry[:6])
                 + (" ..." if len(near_dry) > 6 else "")))
    if not usable:
        sys.exit("No usable truths.")

    # ---- gate D: each truth file is the simulation its own run produced ---------------
    selfc = [selfcheck_truth(tr, truth_obs[tr["id"]], csv_dir) for tr in usable]
    cnt = {k: sum(1 for e in selfc if e["status"] == k) for k in ("ok", "minor", "FAIL", "not checked")}
    rels = [e["rel_to_peak"] for e in selfc if "rel_to_peak" in e]
    print("GATE D  truth vs its own run: %d exact, %d minor (<=%.0e of peak), %d FAIL, %d not checked%s"
          % (cnt["ok"], cnt["minor"], SELF_FAIL_REL, cnt["FAIL"], cnt["not checked"],
             ("  (worst %.1e of the peak)" % max(rels)) if rels else ""))
    for e in selfc:
        if e["status"] in ("FAIL", "minor", "not checked"):
            print("    %-11s %s  (%s)" % (e["status"], e["truth_id"], e["detail"]))
    d_overridden = False
    if cnt["FAIL"]:
        if not args.ignore_selfcheck:
            sys.exit("\nSTOPPED: a saved truth file does not match the simulation its own run produced, so it may be a "
                     "file from a different run. Find out why first. --ignore_selfcheck overrides this, and says so in "
                     "the provenance file.")
        d_overridden = True
        print("  --ignore_selfcheck set: continuing despite the failure.")

    n_on = len(res_on)
    n_off = 0 if res_off is None else len(res_off)
    if args.check_only:
        print("\n--check_only: gates finished, nothing written. A full run would score %d + %d runs x %d truth(s)."
              % (n_on, n_off, len(usable)))
        return

    # ---- re-score -----------------------------------------------------------------------
    out_dir.mkdir(parents=True, exist_ok=True)
    print("\nRe-scoring %d + %d runs x %d truth(s)..." % (n_on, n_off, len(usable)), flush=True)
    series_list = [(SER_ON, res_on)] + ([(SER_OFF, res_off)] if res_off is not None else [])
    long_rows, summ_rows, truth_failed = [], [], []
    for k, tr in enumerate(usable, start=1):
        obs_new = truth_obs[tr["id"]]
        for series, res in series_list:
            mdf, failed = TL.score_truth(obs_new, res, cache, phase_fn)
            for rid in failed:
                truth_failed.append((tr["id"], rid))
            if mdf.empty:
                continue
            mdf.insert(0, "series", series)
            mdf.insert(0, "truth_cc", tr["cc"])
            mdf.insert(0, "truth_f", tr["f"])
            mdf.insert(0, "truth_Ks", tr["Ks"])
            mdf.insert(0, "truth_id", tr["id"])
            mdf["storm_label"] = label
            mdf["rain_scale"] = scale
            long_rows.append(mdf)
            srow = TL.summarize(mdf, series, tr, args.kge_tol)
            srow["storm_label"] = label
            srow["rain_scale"] = scale
            srow["truth_flags"] = tr["flags"]
            summ_rows.append(srow)
        if k % 5 == 0 or k == len(usable):
            print("  %d/%d truths scored" % (k, len(usable)), flush=True)
    if not long_rows:
        sys.exit("Nothing could be scored (every (truth, run) pair failed).")

    long_df = pd.concat(long_rows, ignore_index=True)
    summ_df = pd.DataFrame(summ_rows)
    write_csv_atomic(long_df, out_dir / long_name)
    write_csv_atomic(summ_df, out_dir / summ_name)

    inputs = {p_on.name: {"path": str(p_on), "md5": RS.md5_of(p_on), "runs": int(len(res_on))}}
    if res_off is not None:
        inputs[p_off.name] = {"path": str(p_off), "md5": RS.md5_of(p_off), "runs": int(len(res_off))}
    for extra in (truth_csv, ref_csv, loc_dir / LOC_PROV, forcing_prov):
        if extra.exists():
            inputs[extra.name] = {"path": str(extra), "md5": RS.md5_of(extra)}
    prov = {
        "script": Path(__file__).name,
        "created_local": datetime.now().isoformat(timespec="seconds"),
        "storm_label": label, "rain_scale": scale,
        "python": sys.version.split()[0], "pandas": pd.__version__, "numpy": np.__version__,
        "runner_python_pandas_numpy": [loc_prov.get("python"), loc_prov.get("pandas"), loc_prov.get("numpy")],
        "reference_cc_off": {"file": ref_path.name, "md5": ref_info["md5"], "recorded_md5": ref_info["recorded_md5"]},
        "control_included": res_off is not None,
        "inputs": inputs,
        "gates": {"A_max_abs_diff_m3s": worstA, "B_worst_abs_dev": worstB, "runs_checked": n_checked,
                  "runs_excluded": {k: v for k, v in bad.items()},
                  "phase_metrics_from_scorer": phase_fn is not None,
                  "D_selfcheck": {"counts": cnt, "worst_rel_to_peak": (max(rels) if rels else None), "per_truth": selfc,
                                  "tolerances_rel_to_peak": {"exact": SELF_OK_REL, "accepted_with_note": SELF_FAIL_REL}},
                  "D_overridden": d_overridden},
        "truths": {tr["id"]: {"file": tr["name"], "run_id": tr["run_id"], "Ks_mult": tr["Ks"], "f_RS_abs": tr["f"],
                              "cc_mmhr": tr["cc"], "md5": RS.md5_of(loc_dir / tr["name"]), "flags": tr["flags"]}
                   for tr in usable},
        "dry_truths": [{"truth": n, "reason": w} for n, w in dry],
        "skipped_truths": [{"truth": n, "reason": w} for n, w in skipped],
        "unregistered_files": unregistered,
        "kge_tol": args.kge_tol, "feasible_pbias_pct": RS.FEASIBLE_PBIAS_PCT,
        "note": "Metrics recomputed from the stored Simulated column of each compare CSV against each truth, read "
                "exactly as the scorer reads a truth (rescore_cc_truths_110.py's code); scoring and summary functions "
                "from rescore_truth_location_110.py. Dry truths (no outlet flow) are left out.",
    }
    write_json_atomic(out_dir / prov_name, prov)

    print("\n" + "-" * 128)
    print("STORM %s (rain x%g): %d truth(s) scored, %d dry set aside, %d other problem(s).  Series: %s."
          % (label, scale, len(usable), len(dry), n_other, "cc ON + cc OFF control" if res_off is not None else "cc ON only"))
    TL.report(summ_df, dry + skipped)
    if truth_failed:
        print("\n  WARNING: %d (truth, run) pair(s) could not be scored; first: %s" % (len(truth_failed), truth_failed[0]))
    print("\nSaved to: %s" % out_dir)
    print("  %s   %s   %s" % (long_name, summ_name, prov_name))
    print("Next: look at the table above; the storm-aware analysis (which reads these files) is the next step.")


if __name__ == "__main__":
    main()
