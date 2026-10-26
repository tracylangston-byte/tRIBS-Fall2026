"""
plot_series120_all.py
=====================
Runs all five Series 120 plot scripts one after another, so you only have to type one command.

    plot_series120_storms.py              the three-storm maps, feasible bands and pairwise comparison
    plot_series120_best_kge.py            the best-KGE runs for the 1.0x storm
    plot_series120_ridge_old_vs_new.py    the zero-bias ridge, old build vs tRIBS 6.0.0
    plot_series120_hydrograph_shape.py    the hydrographs of the volume-matched runs and how they miss the truth's shape
    plot_series120_late_rise.py           the small late rise of some volume-matched runs, drawn next to the rain forcing

How to run it (from the lab/ folder, the same place you run your other scripts):

    python plot_series120_all.py

All five scripts must sit in the same folder as this one. Nothing else to set up.

What it does, step by step
    1. Checks that it is being run from lab/ (it looks for ../calibration_work).
    2. Runs each script in turn. You see each script's own messages as it goes.
    3. After each script it checks that the files that script should have made are really there and are new
       (a file left over from an earlier run does not count). This matters because one script
       (the hydrograph one) can finish without an error yet skip a figure when it cannot find the hydrograph files.
    4. If one script fails, it carries on with the rest, and the summary at the end tells you which one failed and why.
    5. Prints a short summary table and exits with an error code only if something failed or was incomplete.

What it does NOT do
    It never changes your results CSVs, your run folders or your model inputs. The scripts only READ the Series 120 results
    files (and, for the hydrograph script, the per-run compare files in csv_exports) and WRITE new figures and small tables
    into  calibration_work/03_comparisons/sensitivity_plots/ . So it is safe to run while the Series 121 sweep is going
    (it uses a little processor time, nothing more). Running it again overwrites the previous figures with the same names.

Options (all optional)
    --only a,b           run only some steps. The names are: storms, best_kge, ridge, shape, late
                         (example:  python plot_series120_all.py --only shape)
    --skip_hydrographs   make only the part of the hydrograph script that needs no hydrograph files (figure 2);
                         the late-rise step needs them, so it is left out
    --dpi N              picture sharpness passed on to every script (default 150; 72 to 600)
    --list               print what would run and which files it should make, then stop

Things to keep in mind
    * plot_series120_ridge_old_vs_new.py still draws the ridge with the older degree-4 polynomial. The thin-plate spline
      we moved to later is not in that script yet, so treat that one figure as the "old method" picture.
    * The hydrograph script's figure 1 needs the per-run files   calibration_work/03_comparisons/csv_exports/<run_id>_compare_obs_sim.csv
      for the volume-matched runs. If they are not there you get figure 2 only, and this script will report the step as INCOMPLETE.
    * The late-rise script needs those same files and also the two rain files in init_data/met (precip_SMF_2014-2014.mdf and
      precip_SMPHQ_2014-2014.mdf). If either is missing it stops with a plain message and this script reports the step as FAILED.
"""

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

# ----------------------------------------------------------------------
# What to run. Each step: key, script, the folder (under sensitivity_plots/) it writes to, the files it should make.
# `needs_hyd` marks the files that need the per-run hydrograph files (so --skip_hydrographs does not expect them).
# ----------------------------------------------------------------------
STEPS = [
    {"key": "storms", "script": "plot_series120_storms.py", "folder": "series120_storms",
     "what": "three-storm maps, feasible bands, pairwise comparison",
     "files": [("fig_series120_storms_1_pbias_maps.png", False),
               ("fig_series120_storms_2_feasible_bands.png", False),
               ("fig_series120_storms_3_pbias_pairs.png", False),
               ("series120_storms_summary.csv", False),
               ("series120_storms_pairs.csv", False)]},
    {"key": "best_kge", "script": "plot_series120_best_kge.py", "folder": "series120",
     "what": "best-KGE runs, 1.0x storm",
     "files": [("fig_series120_best_kge.png", False),
               ("series120_top_kge_runs.csv", False)]},
    {"key": "ridge", "script": "plot_series120_ridge_old_vs_new.py", "folder": "series120",
     "what": "zero-bias ridge, old build vs tRIBS 6.0.0 (older polynomial ridge)",
     "files": [("fig_series120_ridge_old_vs_new.png", False)]},
    {"key": "shape", "script": "plot_series120_hydrograph_shape.py", "folder": "series120_hydrograph_shape",
     "what": "hydrographs of volume-matched runs vs the truth's shape",
     "files": [("fig_series120_hydrograph_shape_1_hydrographs.png", True),
               ("fig_series120_hydrograph_shape_2_shape_vs_f.png", False),
               ("series120_hydrograph_shape_runs.csv", False),
               ("series120_hydrograph_shape_bands.csv", False)]},
    {"key": "late", "script": "plot_series120_late_rise.py", "folder": "series120_late_rise",
     "what": "the small late rise of some volume-matched runs, next to the rain",
     "files": [("fig_series120_late_rise_1_night.png", True),
               ("fig_series120_late_rise_2_size_vs_f.png", True),
               ("series120_late_rise_runs.csv", True),
               ("series120_late_rise_summary.csv", True)]},
]
KEYS = [s["key"] for s in STEPS]


def line(char="=", n=78):
    print(char * n)


def expected_files(step, plots_dir, skip_hyd):
    """Full paths this step should have written."""
    return [plots_dir / step["folder"] / name for name, needs_hyd in step["files"] if not (skip_hyd and needs_hyd)]


def run_step(step, here, plots_dir, args):
    """Run one script. Returns (status, note). status is OK, INCOMPLETE or FAILED."""
    script = here / step["script"]
    if not script.is_file():
        return "FAILED", f"{step['script']} was not found next to this script ({here})."
    cmd = [sys.executable, "-u", str(script), "--dpi", str(args.dpi)]
    if step["key"] == "shape" and args.skip_hydrographs:
        cmd.append("--skip_hydrographs")
    want = expected_files(step, plots_dir, args.skip_hydrographs and step["key"] == "shape")
    started = time.time() - 2.0           # small allowance for file-time rounding
    env = dict(os.environ, MPLBACKEND="Agg", PYTHONUNBUFFERED="1")
    try:
        result = subprocess.run(cmd, cwd=Path.cwd(), env=env)
    except OSError as err:
        return "FAILED", f"could not start the script: {err}"
    if result.returncode != 0:
        return "FAILED", f"the script stopped with exit code {result.returncode}. The last lines of its output, just above, give the reason."
    missing = [p.name for p in want if not p.is_file()]
    stale = [p.name for p in want if p.is_file() and p.stat().st_mtime < started]
    if missing or stale:
        parts = []
        if missing:
            parts.append("not made: " + ", ".join(missing))
        if stale:
            parts.append("not refreshed (old file left over): " + ", ".join(stale))
        return "INCOMPLETE", "; ".join(parts) + ". The script's own messages above say why."
    return "OK", f"{len(want)} file(s) written to {(plots_dir / step['folder']).relative_to(plots_dir.parent.parent.parent)}"


def main():
    ap = argparse.ArgumentParser(description="Run all five Series 120 plot scripts in one go.")
    ap.add_argument("--only", default=None, help="comma-separated steps to run: " + ", ".join(KEYS))
    ap.add_argument("--skip_hydrographs", action="store_true",
                    help="hydrograph script: draw only figure 2 (needs no per-run hydrograph files)")
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--list", action="store_true", help="print what would run, then stop")
    args = ap.parse_args()

    if not (72 <= args.dpi <= 600):
        ap.error("--dpi must be between 72 and 600.")
    chosen = KEYS
    if args.only:
        chosen = [k.strip() for k in args.only.split(",") if k.strip()]
        bad = [k for k in chosen if k not in KEYS]
        if bad or not chosen or len(set(chosen)) != len(chosen):
            ap.error(f"--only must be a comma-separated list drawn from {KEYS} with no repeats (got '{args.only}').")
        chosen = [k for k in KEYS if k in chosen]     # keep the usual order

    here = Path(__file__).resolve().parent
    calib_dir = Path.cwd().parent / "calibration_work"
    plots_dir = calib_dir / "03_comparisons" / "sensitivity_plots"
    steps = [s for s in STEPS if s["key"] in chosen]
    if args.skip_hydrographs and any(s["key"] == "late" for s in steps):
        steps = [s for s in steps if s["key"] != "late"]
        print("NOTE: --skip_hydrographs given, so the late-rise step (which needs the hydrograph files) is left out.\n")
        if not steps:
            print("Nothing is left to run. Drop --skip_hydrographs to run the late-rise step.")
            return 1

    if args.list:
        print(f"Would run from {Path.cwd()}  (scripts in {here}):\n")
        for s in steps:
            print(f"  {s['key']:9s} {s['script']}\n            {s['what']}")
            for p in expected_files(s, plots_dir, args.skip_hydrographs and s["key"] == "shape"):
                print(f"            -> {p}")
        return 0

    if not calib_dir.is_dir():
        print(f"STOPPED: I expected a folder called calibration_work next to the folder you are in, at\n  {calib_dir}\n"
              f"but it is not there. Run this from the lab/ folder, for example:\n"
              f"  cd /workspaces/tRIBS-Fall2026/workspaces/SMF/lab\n  python plot_series120_all.py")
        return 1

    line()
    print(f"Series 120 plots: running {len(steps)} script(s) from {Path.cwd()}")
    line()
    results = []
    t_all = time.time()
    try:
        for i, step in enumerate(steps, 1):
            print(f"\n[{i}/{len(steps)}] {step['script']}  -  {step['what']}")
            line("-")
            t0 = time.time()
            status, note = run_step(step, here, plots_dir, args)
            results.append((step, status, note, time.time() - t0))
            line("-")
            print(f"[{i}/{len(steps)}] {status}")
    except KeyboardInterrupt:
        print("\n\nStopped by you (Ctrl+C). Figures from the steps that had already finished are saved.")
        return 130

    print()
    line()
    print(f"SUMMARY  ({time.time() - t_all:.0f} s in total)")
    line()
    for step, status, note, secs in results:
        print(f"  {status:10s} {step['script']}  ({secs:.0f} s)")
        print(f"             {note}")
    bad = [r for r in results if r[1] != "OK"]
    print()
    if bad:
        print(f"{len(bad)} of {len(results)} step(s) did not finish cleanly. Scroll up to that script's messages for the reason, "
              "fix it, and rerun just that step, for example:\n"
              f"  python plot_series120_all.py --only {bad[0][0]['key']}")
        return 1
    folders = sorted({r[0]["folder"] for r in results})
    print("All steps finished. Figures are in calibration_work/03_comparisons/sensitivity_plots/ "
          f"(folder{'s' if len(folders) > 1 else ''} {', '.join(folders)}).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
