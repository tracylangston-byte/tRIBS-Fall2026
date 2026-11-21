#!/usr/bin/env python3
"""Runoff made on the flood week (Jan 18-25, 2017), split by mechanism, from the .rft files that already exist.

READ-ONLY. It only reads .rft (and .mrf) files in results/. It does not run tRIBS, change an input file, or change a result.
Without --csv it writes nothing at all; it only prints.

WHAT IT ANSWERS
  For each run, how much of the runoff the model MADE ON THE LAND during the flood week came from each mechanism:
      Hsrf    infiltration excess      (tRIBS docs, Table 8.8: "Infiltration-excess Runoff")
      Sbsrf   saturation excess        (tRIBS docs: "Saturation-excess Runoff")
      Psrf    perched return flow      (tRIBS docs: "Perched Return Flow")
      Satsrf  groundwater exfiltration (tRIBS docs: "Groundwater Exfiltration")
  NOTE ON NAMES: some older Project scripts label Sbsrf "subsurface stormflow" and Satsrf "saturation excess". That does not match the
  tRIBS documentation. This script uses the documentation's names.

  The numbers are runoff GENERATED on the land, all cells added, BEFORE routing (Srf in the .mrf is the sum of these four columns).
  They are not flow at the outlet.

THE WINDOW
  The flood week is Jan 18 2017 00:00 to Jan 25 2017 00:00 (168 hours). Time_hr in a .rft file counts hours from the START OF THAT RUN
  (a January run starts Jan 18 00:00; the Oct-to-Mar runs start Oct 1 2016 00:00, so the flood week is Time_hr 2616 to 2784).
  A row at Time_hr = t is counted as the hour ending at t, so the rows used are lo < Time_hr <= hi.
  The script prints the window, the number of rows used (168 expected), and a boundary check (how much any share would change if
  the window were shifted by one hour), so you can see the answer does not depend on that choice.

  Run start: a run name starting Sycamore_oct is taken to start 2016-10-01 00:00; one starting Sycamore_jan, 2017-01-18 00:00.
  Any other name needs --run-start. The script prints the start it used.

HOW FILES ARE FOUND
  The .rft name is the run name, then the run hours, then _00 (no dot before the hours): run Sycamore_jan -> Sycamore_jan0168_00.rft;
  run Sycamore_octmar_ks0p2_cv6_psib0p5_nocl -> Sycamore_octmar_ks0p2_cv6_psib0p5_nocl4368_00.rft. The script matches the run name followed
  directly by digits, so Sycamore_jan does NOT pick up Sycamore_jan_ks0p3 files. If zero or more than one file fits, it says so and skips that
  run (use --rft to name a file yourself).

STOPS INSTEAD OF GUESSING
  A run is skipped, with the reason, if: the file is missing; a column is missing; the rows are not hourly; Time_hr repeats; the file ends
  before the window ends. Nothing is estimated or filled in.

COMMANDS (from workspaces/Sycamore/lab/)
  python rft_split_flood_sycamore.py --list                 # just list the .rft files that exist in results/ (reads no data)
  python rft_split_flood_sycamore.py --dry-run              # finds each default run's file, reads its header, prints the window; no volumes
  python rft_split_flood_sycamore.py                        # the seven default runs (see DEFAULT_RUNS below)
  python rft_split_flood_sycamore.py --runs Sycamore_jan,Sycamore_jan_ks0p3
  python rft_split_flood_sycamore.py --runs Sycamore_jan --rft results/Sycamore_jan0168_00.rft     # name the file yourself (one run only)
  python rft_split_flood_sycamore.py --csv score_out/rft_split_flood.csv                          # also write the table to a new CSV file

WORDS
  share        a type's volume divided by the sum of the four types' volumes, in percent.
  M m3         million cubic meters.
  mm over basin  volume divided by basin area (424.6 km2), so it can be compared with rain in mm.
"""
import argparse
import datetime as dt
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

TYPES = ["Hsrf", "Sbsrf", "Psrf", "Satsrf"]
COLS = ["{}_m3_s".format(t) for t in TYPES]
MEANING = {
    "Hsrf": "infiltration excess",
    "Sbsrf": "saturation excess",
    "Psrf": "perched return flow",
    "Satsrf": "groundwater exfiltration",
}
BASIN_KM2 = 424.6
FLOOD_START = dt.datetime(2017, 1, 18, 0, 0)
FLOOD_END = dt.datetime(2017, 1, 25, 0, 0)
OCT_START = dt.datetime(2016, 10, 1, 0, 0)

# The runs asked about: run 1 (first January run, Ks x1), its snow-on repeat, run 2 (Ks x0.3, loss ON and the loss-OFF repeat),
# run 3 (Ks x0.2 with cv 6, January-only and Oct-to-Mar control) and the PsiB x0.5 run.
DEFAULT_RUNS = [
    "Sycamore_jan",
    "Sycamore_jan_snow",
    "Sycamore_jan_ks0p3",
    "Sycamore_jan_ks0p3_nocl",
    "Sycamore_jan_ks0p2_cv6_nocl",
    "Sycamore_octmar_ks0p2_cv6_nocl",
    "Sycamore_octmar_ks0p2_cv6_psib0p5_nocl",
]


def fmt(x, d=3):
    return "n/a" if x is None or not np.isfinite(x) else "{:.{}f}".format(x, d)


def parse_dt(s):
    for f in ("%Y-%m-%d %H:%M", "%Y-%m-%d"):
        try:
            return dt.datetime.strptime(s, f)
        except ValueError:
            pass
    raise SystemExit("STOPPED: could not read the date '{}'. Use YYYY-MM-DD or 'YYYY-MM-DD HH:MM'.".format(s))


def run_start_for(run, override):
    if override is not None:
        return override, "from --run-start"
    if run.startswith("Sycamore_oct"):
        return OCT_START, "assumed from the name (Sycamore_oct...): 2016-10-01 00:00"
    if run.startswith("Sycamore_jan"):
        return FLOOD_START, "assumed from the name (Sycamore_jan...): 2017-01-18 00:00"
    return None, "unknown"


def find_rft(results_dir, run):
    pat = re.compile(r"^{}(\d+)_(\d+)\.rft$".format(re.escape(run)))
    return sorted(p for p in results_dir.glob("*.rft") if pat.match(p.name))


def read_rft(path):
    """Returns (DataFrame, None) or (None, reason)."""
    try:
        df = pd.read_csv(path)
    except Exception as e:  # noqa: BLE001 - report any read problem to the user
        return None, "could not read the file ({}: {})".format(type(e).__name__, e)
    df.columns = [str(c).strip() for c in df.columns]
    missing = [c for c in ["Time_hr"] + COLS if c not in df.columns]
    if missing:
        return None, "missing column(s) {}; the columns found are: {}".format(missing, list(df.columns))
    if df["Time_hr"].isna().any():
        return None, "Time_hr has blank or unreadable values"
    if df["Time_hr"].duplicated().any():
        return None, "Time_hr repeats (the same hour appears more than once)"
    if len(df) < 2:
        return None, "fewer than two rows"
    steps = np.unique(np.round(np.diff(df["Time_hr"].values), 6))
    if not (len(steps) == 1 and abs(steps[0] - 1.0) < 1e-6):
        return None, "rows are not one hour apart (spacings found: {}); volumes would need a different time step, so none are computed".format(list(steps[:6]))
    return df, None


def window_hours(start):
    lo = (FLOOD_START - start).total_seconds() / 3600.0
    hi = (FLOOD_END - start).total_seconds() / 3600.0
    return lo, hi


def type_volumes(df, lo, hi, inclusive_start=False):
    """Million m3 by type for rows in the window; also the number of rows and how many type values were blank."""
    t = df["Time_hr"].values
    mask = ((t >= lo) & (t < hi)) if inclusive_start else ((t > lo) & (t <= hi))
    sub = df.loc[mask]
    vols = {}
    blanks = 0
    for ty, c in zip(TYPES, COLS):
        col = sub[c]
        blanks += int(col.isna().sum())
        vols[ty] = float(col.sum(skipna=True)) * 3600.0 / 1e6
    return vols, len(sub), blanks, sub


def shares(vols):
    tot = sum(vols.values())
    if tot <= 0:
        return None, tot
    return {ty: 100.0 * vols[ty] / tot for ty in TYPES}, tot


def mrf_check(rft_path, lo, hi, tot_Mm3):
    mp = rft_path.with_suffix(".mrf")
    if not mp.exists():
        return None, "skipped (no matching .mrf file next to the .rft)"
    try:
        m = pd.read_csv(mp, usecols=lambda c: c.strip() in ("Time_hr", "Srf_m3_s"))
    except Exception as e:  # noqa: BLE001
        return None, "skipped (could not read the .mrf: {}: {})".format(type(e).__name__, e)
    m.columns = [c.strip() for c in m.columns]
    if "Srf_m3_s" not in m.columns or "Time_hr" not in m.columns:
        return None, "skipped (the .mrf has no Time_hr / Srf_m3_s columns)"
    sub = m[(m["Time_hr"] > lo) & (m["Time_hr"] <= hi)]
    srf = float(sub["Srf_m3_s"].sum(skipna=True)) * 3600.0 / 1e6
    if srf <= 0:
        return None, "skipped (the .mrf Srf volume in the window is zero)"
    return tot_Mm3 / srf, "sum of the four .rft columns / .mrf Srf over the same rows = {:.4f} (1.0000 expected; .mrf Srf volume {:.3f} M m3)".format(tot_Mm3 / srf, srf)


def process(run, results_dir, args, explicit_rft=None):
    out = {"run": run, "status": None}
    start, start_note = run_start_for(run, args.run_start)
    if start is None:
        out["status"] = "SKIPPED: cannot tell when this run starts; give --run-start"
        return out
    if explicit_rft is not None:
        files = [Path(explicit_rft)]
        if not files[0].exists():
            out["status"] = "NOT FOUND: {} does not exist".format(explicit_rft)
            return out
    else:
        files = find_rft(results_dir, run)
        if not files:
            out["status"] = "NOT FOUND: no .rft file named {}<digits>_<digits>.rft in {}".format(run, results_dir)
            return out
        if len(files) > 1:
            out["status"] = "SKIPPED: more than one file fits ({}); use --rft to name one".format(", ".join(f.name for f in files))
            return out
    path = files[0]
    out["file"] = str(path)
    lo, hi = window_hours(start)
    out.update(start=start, start_note=start_note, lo=lo, hi=hi)
    if lo < 0:
        out["status"] = "SKIPPED: the run starts after the flood week begins (start {})".format(start)
        return out
    df, why = read_rft(path)
    if df is None:
        out["status"] = "SKIPPED: " + why
        return out
    out["last_hr"] = float(df["Time_hr"].max())
    out["first_hr"] = float(df["Time_hr"].min())
    if out["last_hr"] < hi - 1e-6:
        out["status"] = "SKIPPED: the file ends at Time_hr {:g} but the flood week needs hour {:g}; the window is not fully covered, so no share is given".format(out["last_hr"], hi)
        return out
    if args.dry_run:
        out["status"] = "DRY-RUN OK"
        out["rows_total"] = len(df)
        return out
    vols, n, blanks, _ = type_volumes(df, lo, hi)
    sh, tot = shares(vols)
    # M m3 -> mm over the basin: (tot * 1e6 m3) / (BASIN_KM2 * 1e6 m2) * 1000 mm/m = tot * 1000 / BASIN_KM2
    out.update(vols=vols, rows=n, blanks=blanks, shares=sh, total=tot, mm=tot * 1000.0 / BASIN_KM2)
    expected = int(round(hi - lo))
    out["expected_rows"] = expected
    # boundary check: window shifted by one hour (start inclusive, end exclusive)
    vols2, n2, _, _ = type_volumes(df, lo, hi, inclusive_start=True)
    sh2, _ = shares(vols2)
    out["boundary_pp"] = None if (sh is None or sh2 is None) else max(abs(sh[t] - sh2[t]) for t in TYPES)
    out["mrf_ratio"], out["mrf_note"] = mrf_check(path, lo, hi, tot)
    out["status"] = "OK"
    return out


def print_result(r):
    print("\n" + "=" * 100)
    print("RUN {}".format(r["run"]))
    if r.get("file"):
        print("  file: {}".format(r["file"]))
    if r["status"] not in ("OK", "DRY-RUN OK"):
        print("  {}".format(r["status"]))
        return
    print("  run start: {}  ({})".format(r["start"].strftime("%Y-%m-%d %H:%M"), r["start_note"]))
    print("  window: {} to {}  = Time_hr > {:g} and <= {:g}".format(FLOOD_START.strftime("%Y-%m-%d %H:%M"), FLOOD_END.strftime("%Y-%m-%d %H:%M"), r["lo"], r["hi"]))
    print("  file covers Time_hr {:g} to {:g}".format(r["first_hr"], r["last_hr"]))
    if r["status"] == "DRY-RUN OK":
        print("  (dry run: header and window checked; no volumes computed)")
        return
    flag = "" if r["rows"] == r["expected_rows"] else "   <-- NOT the expected {} rows; read the numbers with care".format(r["expected_rows"])
    print("  rows used: {}{}".format(r["rows"], flag))
    if r["blanks"]:
        print("  WARNING: {} blank/NaN values among the four columns in the window; they were left out of the sums".format(r["blanks"]))
    if r["shares"] is None:
        print("  The four columns add to zero in the window: no runoff was made, so there are no shares.")
        return
    print("  {:8s} {:26s} {:>12s} {:>9s}".format("type", "meaning (tRIBS docs)", "volume M m3", "share %"))
    for t in TYPES:
        print("  {:8s} {:26s} {:>12s} {:>9s}".format(t, MEANING[t], fmt(r["vols"][t], 4), fmt(r["shares"][t], 2)))
    print("  {:8s} {:26s} {:>12s} {:>9s}   = {} mm over the basin ({} km2)".format("total", "", fmt(r["total"], 4), "100.00", fmt(r["mm"], 2), BASIN_KM2))
    print("  check against the .mrf: {}".format(r["mrf_note"]))
    if r["boundary_pp"] is not None:
        print("  boundary check: shifting the window by one hour changes the largest share by {:.3f} percentage points".format(r["boundary_pp"]))


def summary_table(results):
    ok = [r for r in results if r["status"] == "OK" and r["shares"] is not None]
    print("\n" + "=" * 100)
    print("SUMMARY (share of runoff made on the land, flood week Jan 18-25 2017, percent; before routing)")
    print("  {:42s} {:>8s} {:>8s} {:>8s} {:>8s} {:>11s}".format("run", "Hsrf", "Sbsrf", "Psrf", "Satsrf", "total M m3"))
    for r in results:
        if r in ok:
            s = r["shares"]
            print("  {:42s} {:>8s} {:>8s} {:>8s} {:>8s} {:>11s}".format(r["run"], fmt(s["Hsrf"], 2), fmt(s["Sbsrf"], 2), fmt(s["Psrf"], 2), fmt(s["Satsrf"], 2), fmt(r["total"], 3)))
        else:
            print("  {:42s} {}".format(r["run"], r["status"].split(":")[0] if r["status"] else "?"))
    print("  Hsrf = infiltration excess, Sbsrf = saturation excess, Psrf = perched return flow, Satsrf = groundwater exfiltration (tRIBS docs, Table 8.8).")


def write_csv(results, path):
    rows = []
    for r in results:
        row = {"run": r["run"], "status": r["status"], "file": r.get("file", "")}
        if r["status"] == "OK":
            row.update(window_start=FLOOD_START.strftime("%Y-%m-%d %H:%M"), window_end=FLOOD_END.strftime("%Y-%m-%d %H:%M"),
                       hour_lo=r["lo"], hour_hi=r["hi"], rows=r["rows"], total_Mm3=r["total"], total_mm=r["mm"], mrf_ratio=r["mrf_ratio"])
            for t in TYPES:
                row["{}_Mm3".format(t)] = r["vols"][t]
                row["{}_pct".format(t)] = None if r["shares"] is None else r["shares"][t]
        rows.append(row)
    p = Path(path)
    if p.exists():
        raise SystemExit("STOPPED: {} already exists; this script does not overwrite files. Pick a new name.".format(p))
    p.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(p, index=False, float_format="%.6f")
    print("\nwrote {}".format(p))


def main(argv=None):
    ap = argparse.ArgumentParser(description="Flood-week runoff by mechanism from existing .rft files (read-only).")
    ap.add_argument("--root", default=".", help="folder that holds results/ (default: the current folder; run from workspaces/Sycamore/lab/)")
    ap.add_argument("--runs", default=None, help="comma list of run names (default: the seven in DEFAULT_RUNS)")
    ap.add_argument("--rft", default=None, help="explicit .rft file; use with exactly one run in --runs")
    ap.add_argument("--run-start", default=None, help="run start 'YYYY-MM-DD HH:MM' (needed only for run names not starting Sycamore_oct / Sycamore_jan)")
    ap.add_argument("--list", action="store_true", help="only list the .rft files in results/")
    ap.add_argument("--dry-run", action="store_true", help="find files, check headers and the window; print no volumes")
    ap.add_argument("--csv", default=None, help="also write the table to this NEW csv file (never overwrites)")
    args = ap.parse_args(argv)
    if args.run_start:
        args.run_start = parse_dt(args.run_start)

    results_dir = Path(args.root) / "results"
    if not results_dir.is_dir():
        print("STOPPED: cannot find the folder {}. Run this from workspaces/Sycamore/lab/ or give --root.".format(results_dir))
        return 2

    if args.list:
        files = sorted(results_dir.glob("*.rft"))
        print("{} .rft file(s) in {}".format(len(files), results_dir))
        for f in files:
            print("  {:80s} {:>10.1f} KB".format(f.name, f.stat().st_size / 1024.0))
        return 0

    runs = [s.strip() for s in args.runs.split(",") if s.strip()] if args.runs else list(DEFAULT_RUNS)
    if args.rft and len(runs) != 1:
        print("STOPPED: --rft names one file, so give exactly one run in --runs.")
        return 2

    results = []
    for run in runs:
        r = process(run, results_dir, args, explicit_rft=args.rft)
        results.append(r)
        print_result(r)
    if not args.dry_run:
        summary_table(results)
        if args.csv:
            write_csv(results, args.csv)
    return 0


if __name__ == "__main__":
    sys.exit(main())
