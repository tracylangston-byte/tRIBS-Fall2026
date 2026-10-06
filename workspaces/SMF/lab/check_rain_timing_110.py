"""
check_rain_timing_110.py   (version 1)

PURPOSE
  Put four things side by side, in 15-minute blocks, for the four storms (12 Aug, 19 Aug, 8 Sep, 27 Sep 2014):
    1. the rain the MODEL WAS FED at each rain gauge (SMF and SMPHQ forcing files, mm/hr),
    2. the raw rain record from the gauge workbook (tipping-bucket records, summed into the same 15-minute blocks),
    3. the real streamflow at the SMF gauge (cfs),
    4. the simulated outlet flow of two ms1 runs (cfs), read from their kept .qout files.
  It answers three plain questions:
    - Is what the model was fed the same as the raw gauge record, block by block?   (forcing vs workbook)
    - Is there rain in the forcing at about the time the real flow starts?           (timing)
    - On 8 Sep, which rain period lines up with the model's third pulse near 07:00?  (compare with the model column)

WHAT IT READS (nothing is changed)
  ../init_data/met/precip_SMF_2014-2014.mdf and precip_SMPHQ_2014-2014.mdf     the forcing (header Year,Month,Day,Hour,R_mm)
  ../init_data/met/SMF_Observations_1993-2025.xlsx    sheets Discharge, SMF Rain, SMPHQ Rain (as the other scripts read them)
  ../calibration_work/kept_qout_110/<run id>_Outlet.qout   two ms1 runs (see --run_id)

WHAT IT WRITES (two NEW files, in a new folder; it refuses to overwrite unless you add --overwrite)
  ../calibration_work/03_comparisons/summary_tables/rain_timing_110/rain_timing_110.csv      every 15-minute block of all four windows
  ../calibration_work/03_comparisons/summary_tables/rain_timing_110/fig_rain_timing_110.png  the picture

RUN (from the lab/ folder)
  python check_rain_timing_110.py
  python check_rain_timing_110.py --run_id SMF_20140812_110_cc722p343129_ms1on    (choose other runs; any number)
  python check_rain_timing_110.py --no_plots                                      (table and CSV only)

HOW THE FORCING IS READ (same as check_storm_rain_110.py, which reproduced the storm list)
  - R is a RATE in mm/hr. Four rows share each hour (15-minute data), evenly spaced: the first row of an hour is
    at :00, then :15, :30, :45. A block's rain in mm is rate x 0.25 h. A block is labelled by its START time.
  - Workbook rain: the column 'Incremental inches' x 25.4 = mm, each tip placed in the 15-minute block that
    contains its time stamp.
  - Real flow in a block: the HIGHEST reading in the block (a dot means there was a reading; '.' means none).
  - Model flow in a block: the highest value in the block (the .qout has a value every 3.75 minutes).

THE ONLY PASS/FAIL IN THIS SCRIPT (a plain tolerance, fixed before it was run)
  Forcing and workbook AGREE in a storm window if no 15-minute block differs by more than 0.01 mm of rain.
  Everything else is printed as numbers and left for you to read: this is a look at the data, not a test.

WHAT IT CANNOT TELL YOU
  Two gauge points are not the whole basin: rain over the basin can differ from both gauges, and the model spreads
  these two gauges over the whole basin. Real flow comes from wherever the rain really fell.
"""
import argparse
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

VERSION = "1"
MET = Path("../init_data/met")
FORCING = {"SMF": MET / "precip_SMF_2014-2014.mdf", "SMPHQ": MET / "precip_SMPHQ_2014-2014.mdf"}
XLSX = MET / "SMF_Observations_1993-2025.xlsx"
SHEET_Q = "Discharge"
SHEET_RAIN = {"SMF": "SMF Rain", "SMPHQ": "SMPHQ Rain"}
SKIPROWS = 6
KEEP = Path("../calibration_work/kept_qout_110")
OUT = Path("../calibration_work/03_comparisons/summary_tables/rain_timing_110")
CSV_NAME, FIG_NAME = "rain_timing_110.csv", "fig_rain_timing_110.png"

CFS_TO_CMS = 0.0283168
MM_PER_IN = 25.4
ORIGIN = pd.Timestamp("2014-08-01")
QOUT_NAMES = ["Time_hr", "Qstrm_m3s", "Hlev_m"]
BLOCK = "15min"
BLOCK_HR = 0.25
TOL_BLOCK_MM = 0.01              # forcing vs workbook: largest allowed difference in any 15-minute block (mm)
FLOW_CFS = 1.0                   # "flow starts" = first value at or above this
QUIET_CFS = 0.5                  # a block is quiet if no rain, no real reading above 0 and the model is under this
GATE_DAY_TOL_MM = 0.15           # day totals against the numbers recorded in handoff v15

# display windows (start, end exclusive); the storm windows used for scoring are longer
STORMS = [("12 Aug", "2014-08-12 16:00", "2014-08-12 22:00", "real peak 1586 cfs"),
          ("19 Aug", "2014-08-19 04:00", "2014-08-19 12:00", "real peak 13 cfs"),
          ("8 Sep", "2014-09-08 00:00", "2014-09-08 10:00", "three real pulses, highest 442 cfs"),
          ("27 Sep", "2014-09-27 11:00", "2014-09-27 18:00", "real peak 252 cfs")]
# day totals in mm (SMF, SMPHQ) found by check_storm_rain_110.py (handoff v15)
KNOWN_DAY_MM = {"2014-08-12": (61.0, 65.8), "2014-08-19": (15.0, 23.1), "2014-09-08": (92.5, 106.4), "2014-09-27": (30.5, 21.6)}
DEFAULT_RUNS = ["SMF_20140812_110_cc722p343129_ms1on",   # best on the combined 12 Aug + 8 Sep score
                "SMF_20140812_110_cc107p236726_ms1on"]   # best on 12 Aug alone (low cc)

# colours (checked for colour-blind separation): rain gauges teal and orange; real flow black; model runs blue
C_SMF, C_SMPHQ = "#00897b", "#c2570c"
C_REAL, C_MODEL = "#1a1a1a", ["#2a6fd6", "#7aa7e6"]
BG, GRID, TXT, TXT2 = "#fcfcfb", "#e6e6e3", "#1a1a1a", "#5c5c58"


# ------------------------------------------------------------------
# Reading
# ------------------------------------------------------------------
def read_forcing(path):
    """Rain rate (mm/hr) indexed by block start time, plus the number of rows that sit in hours without exactly 4 rows."""
    if not path.exists():
        sys.exit("Cannot find %s. Run this from the lab/ folder." % path)
    df = pd.read_csv(path)
    df.columns = [str(c).strip() for c in df.columns]
    if "Year" not in df.columns:                      # an older, whitespace-delimited file
        df = pd.read_csv(path, sep=r"\s+")
        df.columns = [str(c).strip() for c in df.columns]
    rcols = [c for c in df.columns if c.upper().startswith("R")]
    if not rcols or not {"Year", "Month", "Day", "Hour"} <= set(df.columns):
        sys.exit("%s: expected the columns Year, Month, Day, Hour and R_mm; found %s." % (path, list(df.columns)))
    rate = pd.to_numeric(df[rcols[0]], errors="coerce")
    rate = rate.where(rate < 9999, np.nan)            # 9999.99 is the NO_DATA flag
    keys = ["Year", "Month", "Day", "Hour"]
    grp = df.groupby(keys, sort=False)
    n_in_hour = grp["Year"].transform("size")
    pos = grp.cumcount()
    minute = pos * 60.0 / n_in_hour
    t = pd.to_datetime(dict(year=df["Year"], month=df["Month"], day=df["Day"], hour=df["Hour"])) + pd.to_timedelta(minute, unit="m")
    s = pd.Series(rate.to_numpy(float), index=pd.DatetimeIndex(t)).sort_index()
    return s, int((n_in_hour != 4).sum())


def read_sheet(sheet):
    if not XLSX.exists():
        sys.exit("Cannot find %s. Run this from the lab/ folder." % XLSX)
    d = pd.read_excel(XLSX, sheet_name=sheet, skiprows=SKIPROWS)
    d.columns = [str(c).strip() for c in d.columns]
    d["datetime"] = pd.to_datetime(d["Date"].astype(str) + " " + d["Time"].astype(str))
    return d.dropna(subset=["datetime"]).set_index("datetime").sort_index()


def read_workbook_rain_mm(sheet):
    d = read_sheet(sheet)
    col = [c for c in d.columns if c.lower().startswith("incremental")]
    if not col:
        sys.exit("Sheet %s has no 'Incremental inches' column; found %s." % (sheet, list(d.columns)))
    mm = (pd.to_numeric(d[col[0]], errors="coerce") * MM_PER_IN).dropna()
    return mm.groupby(mm.index.floor(BLOCK)).sum()


def read_real_flow():
    d = read_sheet(SHEET_Q)
    return pd.to_numeric(d["cfs"], errors="coerce").dropna()


def read_qout_cfs(path):
    """Simulated outlet flow in cfs at the .qout's own time step (hours from 1 Aug 2014 00:00; tRIBS 6.0.0 is comma-delimited)."""
    with open(path) as fh:
        fh.readline()
        first = fh.readline()
    if "," in first:
        raw = pd.read_csv(path, sep=",", header=0, names=QOUT_NAMES)
    else:
        raw = pd.read_csv(path, sep=r"\s+", skiprows=1, names=QOUT_NAMES)
    t = pd.to_datetime(raw["Time_hr"] * 3600, unit="s", origin=ORIGIN)
    return pd.Series(raw["Qstrm_m3s"].to_numpy(float) / CFS_TO_CMS, index=pd.DatetimeIndex(t))


def short_name(run_id):
    m = re.search(r"cc(\d+)p", run_id)
    return "cc %s" % m.group(1) if m else run_id


# ------------------------------------------------------------------
# One storm window
# ------------------------------------------------------------------
def storm_frame(t0, t1, forcing, rawmm, flow, models):
    idx = pd.date_range(t0, t1, freq=BLOCK, inclusive="left")
    f = pd.DataFrame(index=idx)
    for g in ("SMF", "SMPHQ"):
        f[g + "_forcing_mm_hr"] = forcing[g].reindex(idx)
        f[g + "_workbook_mm_hr"] = rawmm[g].reindex(idx).fillna(0.0) / BLOCK_HR
    fl = flow[(flow.index >= t0) & (flow.index < t1)]
    f["real_cfs_max"] = fl.groupby(fl.index.floor(BLOCK)).max().reindex(idx)
    f["real_n_readings"] = fl.groupby(fl.index.floor(BLOCK)).count().reindex(idx).fillna(0).astype(int)
    for name, s in models.items():
        ss = s[(s.index >= t0) & (s.index < t1)]
        f["model_%s_cfs_max" % name] = ss.groupby(ss.index.floor(BLOCK)).max().reindex(idx)
    return f


def fnum(v, nd=1):
    return "n/a" if pd.isna(v) else ("%.*f" % (nd, v))


def fflow(v):
    if pd.isna(v):
        return "."
    return "%.0f" % v if v >= 10 else "%.1f" % v


def hhmm(ts):
    return pd.Timestamp(ts).strftime("%H:%M")


def first_at_or_above(s, level):
    s = s.dropna()
    hit = s[s >= level]
    return (hit.index[0], float(hit.iloc[0])) if len(hit) else None


def print_storm(name, note, t0, t1, f, flow, models):
    print("\n" + "=" * 100)
    print("%s 2014   window %s to %s     [%s]" % (name, hhmm(t0), hhmm(t1), note))
    print("=" * 100)
    runs = list(models)
    head = "  %-6s %10s %11s %10s" % ("block", "SMF rain", "SMPHQ rain", "real flow") + "".join(" %12s" % ("model " + short_name(r)) for r in runs)
    unit = "  %-6s %10s %11s %10s" % ("start", "mm/hr", "mm/hr", "cfs (max)") + "".join(" %12s" % "cfs (max)" for _ in runs)
    print(head)
    print(unit)
    quiet_run = False
    for ts, r in f.iterrows():
        model_hi = max([r["model_%s_cfs_max" % m] for m in runs if not pd.isna(r["model_%s_cfs_max" % m])] or [0.0])
        rain_hi = max(r["SMF_forcing_mm_hr"] if not pd.isna(r["SMF_forcing_mm_hr"]) else 0.0,
                      r["SMPHQ_forcing_mm_hr"] if not pd.isna(r["SMPHQ_forcing_mm_hr"]) else 0.0)
        active = rain_hi > 0 or (r["real_n_readings"] > 0 and (r["real_cfs_max"] or 0) > 0) or model_hi >= QUIET_CFS
        if not active:
            if not quiet_run:
                print("  ...    (quiet blocks skipped: no rain, no real flow, model under %.1f cfs)" % QUIET_CFS)
                quiet_run = True
            continue
        quiet_run = False
        print("  %-6s %10s %11s %10s" % (hhmm(ts), fnum(r["SMF_forcing_mm_hr"]), fnum(r["SMPHQ_forcing_mm_hr"]), fflow(r["real_cfs_max"]))
              + "".join(" %12s" % fflow(r["model_%s_cfs_max" % m]) for m in runs))
    print("  (real flow '.' = no gauge reading in that block; model '.' = run not available)")

    # ---- summary numbers -----------------------------------------------------------------------
    print("\n  SUMMARY")
    for g in ("SMF", "SMPHQ"):
        fc, wb = f[g + "_forcing_mm_hr"], f[g + "_workbook_mm_hr"]
        tot_f, tot_w = np.nansum(fc.to_numpy()) * BLOCK_HR, wb.sum() * BLOCK_HR
        diff = ((fc - wb) * BLOCK_HR).abs()
        worst = diff.max()
        verdict = "AGREE" if (not pd.isna(worst) and worst <= TOL_BLOCK_MM and fc.notna().all()) else "DIFFER"
        pk = fc.idxmax() if fc.notna().any() else None
        n_missing = int(fc.isna().sum())
        how = "largest difference %s mm" % fnum(worst, 2)
        if n_missing:
            how = "largest difference among blocks that are there %s mm; %d block%s missing in the forcing" % (
                fnum(worst, 2), n_missing, "" if n_missing == 1 else "s")
        print("  %-6s rain in window: forcing %6.1f mm, workbook %6.1f mm | block by block: %s (%s) | "
              "biggest block %s at %s" % (g, tot_f, tot_w, verdict, how, (fnum(fc.max()) + " mm/hr") if pk is not None else "n/a",
                                         hhmm(pk) if pk is not None else "n/a"))
        if verdict == "DIFFER":
            bad = diff[diff > TOL_BLOCK_MM].index.tolist() + fc[fc.isna()].index.tolist()
            for b in sorted(set(bad))[:8]:
                print("           block %s: forcing %s mm/hr, workbook %s mm/hr" % (hhmm(b), fnum(fc[b]), fnum(wb[b])))
        first_rain = fc[fc > 0]
        if len(first_rain):
            print("         first block with rain in the forcing: %s;  last: %s" % (hhmm(first_rain.index[0]), hhmm(first_rain.index[-1])))
    fl = flow[(flow.index >= t0) & (flow.index < t1)]
    st = first_at_or_above(fl, FLOW_CFS)
    if st is None:
        print("  real flow: no reading of %g cfs or more in the window" % FLOW_CFS)
    else:
        print("  real flow: first reading of %g cfs or more at %s; highest reading %.0f cfs at %s" % (FLOW_CFS, hhmm(st[0]), fl.max(), hhmm(fl.idxmax())))
    for m in runs:
        s = models[m]
        ss = s[(s.index >= t0) & (s.index < t1)]
        if ss.empty:
            print("  model %-8s: no values in the window" % short_name(m))
            continue
        st = first_at_or_above(ss, FLOW_CFS)
        if st is None:
            print("  model %-8s: never reaches %g cfs; highest %.2f cfs" % (short_name(m), FLOW_CFS, ss.max()))
        else:
            print("  model %-8s: first reaches %g cfs at %s; highest %.0f cfs at %s" % (short_name(m), FLOW_CFS, hhmm(st[0]), ss.max(), hhmm(ss.idxmax())))


# ------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------
def make_figure(path, frames, flow, models):
    runs = list(models)
    fig, axes = plt.subplots(2, len(frames), figsize=(18, 8.2), gridspec_kw={"height_ratios": [1, 2.2], "hspace": 0.12, "wspace": 0.22})
    fig.patch.set_facecolor(BG)
    rain_max = max(np.nanmax(f[["SMF_forcing_mm_hr", "SMPHQ_forcing_mm_hr"]].to_numpy(float)) for _, _, f in frames) * 1.1
    w = 7.0 / 1440.0
    for j, (name, note, f) in enumerate(frames):
        t0, t1 = f.index[0], f.index[-1] + pd.Timedelta(BLOCK)
        x = mdates.date2num(f.index.to_pydatetime())
        a, b = axes[0, j], axes[1, j]
        for ax in (a, b):
            ax.set_facecolor(BG)
            ax.grid(True, color=GRID, lw=0.8)
            ax.set_axisbelow(True)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            for sp in ("left", "bottom"):
                ax.spines[sp].set_color(GRID)
            ax.tick_params(colors=TXT2, labelsize=9)
        a.bar(x, f["SMF_forcing_mm_hr"].fillna(0.0), width=w, align="edge", color=C_SMF)
        a.bar(x + 7.5 / 1440.0, f["SMPHQ_forcing_mm_hr"].fillna(0.0), width=w, align="edge", color=C_SMPHQ)
        a.set_ylim(0, rain_max)
        a.set_title("%s\n%s" % (name, note), loc="left", fontsize=11, color=TXT, pad=8)
        if j == 0:
            a.set_ylabel("rain in the forcing\n(mm/hr, 15-min blocks)", color=TXT2, fontsize=10)
        else:
            a.tick_params(labelleft=False)
        fl = flow[(flow.index >= t0) & (flow.index < t1)]
        for k, m in enumerate(runs):
            s = models[m]
            ss = s[(s.index >= t0) & (s.index < t1)]
            b.plot(mdates.date2num(ss.index.to_pydatetime()), ss.to_numpy(), color=C_MODEL[k % len(C_MODEL)], lw=2.4 if k == 0 else 1.8,
                   ls="-" if k == 0 else "--", zorder=3)
        b.plot(mdates.date2num(fl.index.to_pydatetime()), fl.to_numpy(), color=C_REAL, lw=1.4, marker="o", ms=3.5, zorder=4)
        b.set_ylim(bottom=0)
        if j == 0:
            b.set_ylabel("flow (cfs)", color=TXT2, fontsize=10)
        lo, hi = mdates.date2num(t0), mdates.date2num(t1)
        for ax in (a, b):
            ax.set_xlim(lo, hi)
            ax.xaxis_date()
            ax.xaxis.set_major_locator(mdates.HourLocator(interval=1 if (t1 - t0) <= pd.Timedelta(hours=7) else 2))
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        a.tick_params(labelbottom=False)
    handles = [Patch(color=C_SMF, label="SMF rain (forcing)"), Patch(color=C_SMPHQ, label="SMPHQ rain (forcing)"),
               Line2D([], [], color=C_REAL, lw=1.4, marker="o", ms=3.5, label="real gauge (each reading)")]
    for k, m in enumerate(runs):
        handles.append(Line2D([], [], color=C_MODEL[k % len(C_MODEL)], lw=2.4 if k == 0 else 1.8, ls="-" if k == 0 else "--",
                              label="model run %s (%s)" % (short_name(m), "best on 12 Aug + 8 Sep" if k == 0 else "run 2")))
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), frameon=False, fontsize=10, labelcolor=TXT, bbox_to_anchor=(0.5, 0.005))
    fig.suptitle("Rain the model was fed, the real flow, and the model's flow, storm by storm", x=0.012, ha="left", fontsize=14, color=TXT, y=0.985)
    fig.text(0.012, 0.055, "Each rain bar is one 15-minute block that starts at its tick position; SMF bar on the left half, SMPHQ on the right half. "
             "The rain scale is the same in all four storms.", fontsize=9, color=TXT2)
    fig.subplots_adjust(left=0.07, right=0.99, top=0.88, bottom=0.14)
    fig.savefig(path, dpi=130, facecolor=BG)
    plt.close(fig)


# ------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Rain in the forcing, raw rain record, real flow and model flow, 15-minute blocks, four storms.")
    ap.add_argument("--run_id", nargs="*", default=DEFAULT_RUNS, help="ms1 run ids whose kept .qout to overlay (default: two runs)")
    ap.add_argument("--overwrite", action="store_true", help="replace the two output files if they already exist")
    ap.add_argument("--no_plots", action="store_true")
    args = ap.parse_args()

    print("=" * 100)
    print("Series 110 -- rain timing check, 15-minute blocks (version %s)    pandas %s  numpy %s" % (VERSION, pd.__version__, np.__version__))
    print("=" * 100)
    csv_path, fig_path = OUT / CSV_NAME, OUT / FIG_NAME
    if not args.overwrite:
        exists = [p.name for p in (csv_path, fig_path) if (p.exists() and not (p == fig_path and args.no_plots))]
        if exists:
            sys.exit("%s already %s in %s. Nothing was changed. Add --overwrite to replace %s." % (
                " and ".join(exists), "exists" if len(exists) == 1 else "exist", OUT, "it" if len(exists) == 1 else "them"))

    # ---- inputs ----------------------------------------------------------------------------------
    forcing, bad_rows = {}, {}
    for g, p in FORCING.items():
        forcing[g], bad_rows[g] = read_forcing(p)
    rawmm = {g: read_workbook_rain_mm(SHEET_RAIN[g]) for g in SHEET_RAIN}
    flow = read_real_flow()

    print("\nINPUT CHECKS")
    for g in forcing:
        s = forcing[g]
        print("  forcing %-5s: %d rows, %s to %s; rows in hours that do not have exactly 4 rows: %d; rows with the NO_DATA flag: %d"
              % (g, len(s), s.index[0].strftime("%d %b %H:%M"), s.index[-1].strftime("%d %b %H:%M"), bad_rows[g], int(s.isna().sum())))
    worst_gate = 0.0
    for day, (a, b) in KNOWN_DAY_MM.items():
        d0 = pd.Timestamp(day)
        got = [forcing[g][d0:d0 + pd.Timedelta(hours=23, minutes=59)].sum() * BLOCK_HR for g in ("SMF", "SMPHQ")]
        worst_gate = max(worst_gate, abs(got[0] - a), abs(got[1] - b))
    print("  [%s] day totals of the forcing on the four storm days against handoff v15 (SMF / SMPHQ mm): largest difference %.2f mm (limit %.2f)"
          % ("PASS" if worst_gate <= GATE_DAY_TOL_MM else "WARN", worst_gate, GATE_DAY_TOL_MM))
    if worst_gate > GATE_DAY_TOL_MM:
        print("         the forcing files are not the ones that gave 61.0/65.8, 15.0/23.1, 92.5/106.4 and 30.5/21.6 mm; read the blocks below with care")
    models = {}
    for rid in args.run_id:
        p = KEEP / ("%s_Outlet.qout" % rid)
        if not p.exists():
            print("  [NOTE] run %s: %s not found, left out" % (rid, p.name))
            continue
        models[rid] = read_qout_cfs(p)
        s = models[rid]
        print("  model %-38s: %d values, hours %.2f to %.2f" % (rid, len(s), (s.index[0] - ORIGIN) / pd.Timedelta(hours=1), (s.index[-1] - ORIGIN) / pd.Timedelta(hours=1)))

    # ---- the four storms -----------------------------------------------------------------------------
    frames = []
    for name, a, b, note in STORMS:
        t0, t1 = pd.Timestamp(a), pd.Timestamp(b)
        f = storm_frame(t0, t1, forcing, rawmm, flow, models)
        frames.append((name, note, f))
        print_storm(name, note, t0, t1, f, flow, models)

    # ---- files ----------------------------------------------------------------------------------------
    OUT.mkdir(parents=True, exist_ok=True)
    allf = pd.concat([f.assign(storm=name).rename_axis("block_start").reset_index() for name, _, f in frames], ignore_index=True)
    allf = allf[["storm", "block_start"] + [c for c in allf.columns if c not in ("storm", "block_start")]]
    tmp = csv_path.with_suffix(".csv.tmp")
    allf.to_csv(tmp, index=False, float_format="%.4f")
    os.replace(tmp, csv_path)
    made = [CSV_NAME]
    if not args.no_plots:
        try:
            make_figure(fig_path, frames, flow, models)
            made.append(FIG_NAME)
        except Exception as exc:
            print("\n  (figure skipped: %s -- the table and CSV are saved regardless)" % exc)
    print("\nSaved to: %s" % OUT.resolve())
    print("  " + "   ".join(made))


if __name__ == "__main__":
    main()
