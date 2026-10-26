"""
check_late_flow_120.py
======================
Looks at the small extra flow that some high-f runs show hours after the main flood
(about 07:10 to 07:45 on 13 Aug 2014, roughly 1% of the peak, never in the synthetic truth)
and answers ONE plain question:

    Is there RAIN in the forcing just before that late flow starts?

If there is, the late flow is probably a response to a small late rain that only some runs
turn into runoff. If there is no rain at all, the late flow is not rain-driven: it is water
the run stored during the storm and is letting out later.

It READS three things and CHANGES NOTHING (it prints to the screen and writes no files):
    ../init_data/met/precip_SMF_2014-2014.mdf        the rain the model was fed, SMF gauge
    ../init_data/met/precip_SMPHQ_2014-2014.mdf      the same for the SMPHQ gauge
    calibration_work/03_comparisons/sensitivity_plots/series120_hydrograph_shape/series120_hydrograph_shape_runs.csv
        (written by plot_series120_hydrograph_shape.py; it says which run has the largest late flow)
    calibration_work/03_comparisons/csv_exports/<run_id>_compare_obs_sim.csv    those runs' hydrographs

HOW TO RUN (from the lab/ folder, like your other scripts)
    python Series110_plots/check_late_flow_120.py
(or just  python check_late_flow_120.py  if the script is in lab/ itself)
Run  plot_series120_hydrograph_shape.py  (or plot_series120_all.py)  first so the runs table exists.

WHAT YOU GET
    1. For each rain gauge: how much rain the forcing has between 21:00 on 12 Aug and 09:00 on 13 Aug,
       and the first and last 15-minute block with rain.
    2. For each storm (0.8x, 1.0x, 1.25x): a table in 15-minute blocks of the rain at both gauges, the
       truth's flow, and the flow of the ONE volume-matched run with the largest late flow. Only blocks
       with rain or flow are printed; quiet stretches are skipped.
    3. For that run: where its late rise starts and peaks (found from the run's own 5-minute flow, not from a flow
       level), the strongest rain block after midnight (both gauges added), how many minutes after that block begins
       the flow gets 0.10 m3/s above where the rise started, the rain in the hour before the rise, and how much extra
       water the run puts out after midnight compared with the truth.

THINGS TO KNOW
    * The 0.8x and 1.25x storms use the rain multiplied by 0.8 and 1.25. A zero stays zero, so WHEN it rains
      is the same in all three storms; only the amounts change. This script reads the unscaled record.
    * The rain file holds a RATE (mm/hr), four rows per hour (:00, :15, :30, :45). A block is labelled by
      its START time. Flow in a block is the highest 5-minute value in it.
    * Two gauges are not the whole basin, and this check cannot say what the water did inside the model.
      It only says whether rain was there to cause the flow.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

WINDOW = ("2014-08-12 21:00", "2014-08-13 09:00")     # start, end (end not included)
BLOCK = "15min"
BLOCK_HR = 0.25
LATE_FROM = "2014-08-13 05:00"  # the late rise is looked for from this time on
FOOT_LOOKBACK_MIN = 120         # the start of the rise = the lowest (run minus truth) in this many minutes before its peak
STEP_UP_M3S = 0.10              # "the flow has really started to rise" = this far above where the rise started
STEP_MIN = 5                    # the hydrograph files are on a 5-minute grid
STORMS = ["0.8x storm", "1.0x storm", "1.25x storm"]
RUN_FILE_SUFFIX = "_compare_obs_sim.csv"
NO_DATA = 9999


def read_forcing(path):
    """Rain rate (mm/hr) indexed by block start time."""
    if not path.exists():
        sys.exit(f"Cannot find {path}\nRun this from the lab/ folder (the folder next to init_data and calibration_work).")
    df = pd.read_csv(path)
    df.columns = [str(c).strip() for c in df.columns]
    if "Year" not in df.columns:                       # an older, whitespace-delimited file
        df = pd.read_csv(path, sep=r"\s+")
        df.columns = [str(c).strip() for c in df.columns]
    rcols = [c for c in df.columns if c.upper().startswith("R")]
    if not rcols or not {"Year", "Month", "Day", "Hour"} <= set(df.columns):
        sys.exit(f"{path.name}: expected the columns Year, Month, Day, Hour and R_mm; found {list(df.columns)}.")
    rate = pd.to_numeric(df[rcols[0]], errors="coerce")
    rate = rate.where(rate < NO_DATA, np.nan)          # 9999.99 is the NO_DATA flag
    keys = ["Year", "Month", "Day", "Hour"]
    grp = df.groupby(keys, sort=False)
    n_in_hour = grp["Year"].transform("size")
    minute = grp.cumcount() * 60.0 / n_in_hour
    t = pd.to_datetime(dict(year=df["Year"], month=df["Month"], day=df["Day"], hour=df["Hour"])) + pd.to_timedelta(minute, unit="m")
    return pd.Series(rate.to_numpy(float), index=pd.DatetimeIndex(t)).sort_index()


def read_run(csv_dir, run_id):
    p = csv_dir / (run_id + RUN_FILE_SUFFIX)
    if not p.exists():
        return None
    d = pd.read_csv(p, index_col=0)
    d.index = pd.to_datetime(d.index)
    return d


def hhmm(ts):
    return pd.Timestamp(ts).strftime("%H:%M")


def main():
    ap = argparse.ArgumentParser(description="Is there rain just before the late flow of the high-f runs?")
    ap.add_argument("--runs_csv", type=Path, default=None, help="the runs table written by plot_series120_hydrograph_shape.py")
    ap.add_argument("--csv_dir", type=Path, default=None, help="folder with the <run_id>_compare_obs_sim.csv files")
    ap.add_argument("--met_dir", type=Path, default=None, help="folder with the two precip_*.mdf rain files")
    args = ap.parse_args()

    root = Path.cwd().parent
    calib = root / "calibration_work"
    runs_csv = args.runs_csv or (calib / "03_comparisons" / "sensitivity_plots" / "series120_hydrograph_shape" / "series120_hydrograph_shape_runs.csv")
    csv_dir = args.csv_dir or (calib / "03_comparisons" / "csv_exports")
    met = args.met_dir or (root / "init_data" / "met")

    if not runs_csv.exists():
        sys.exit(f"Cannot find {runs_csv}\nRun plot_series120_hydrograph_shape.py (or plot_series120_all.py) first, from the lab/ folder.")
    runs = pd.read_csv(runs_csv)
    need = {"storm", "run_id", "f_group", "flow_outside_window_max_m3s", "flow_outside_window_time"}
    if not need <= set(runs.columns):
        sys.exit(f"{runs_csv.name} does not have the late-flow columns {sorted(need - set(runs.columns))}. "
                 "Re-run the newest plot_series120_hydrograph_shape.py to rebuild it.")

    t0, t1 = pd.Timestamp(WINDOW[0]), pd.Timestamp(WINDOW[1])
    idx = pd.date_range(t0, t1, freq=BLOCK, inclusive="left")
    forcing = {"SMF": read_forcing(met / "precip_SMF_2014-2014.mdf"), "SMPHQ": read_forcing(met / "precip_SMPHQ_2014-2014.mdf")}
    rain = {g: s.reindex(idx) for g, s in forcing.items()}

    print("=" * 96)
    print(f"Late flow of the high-f runs against the rain in the forcing, {t0:%d %b %H:%M} to {t1:%d %b %H:%M} 2014")
    print("=" * 96)
    print("\n1. RAIN IN THE FORCING (unscaled record)")
    any_rain_after_storm = False
    for g, s in rain.items():
        n_missing = int(s.isna().sum())
        wet = s[s > 0]
        total = float(np.nansum(s.to_numpy())) * BLOCK_HR
        if len(wet):
            print(f"  {g:<6} {total:6.1f} mm in the window; first rain block {hhmm(wet.index[0])} ({wet.index[0]:%d %b}), "
                  f"last {hhmm(wet.index[-1])} ({wet.index[-1]:%d %b}); biggest block {s.max():.1f} mm/hr at {hhmm(s.idxmax())}")
            if (wet.index >= pd.Timestamp("2014-08-13 00:00")).any():
                any_rain_after_storm = True
        else:
            print(f"  {g:<6} no rain at all in the window")
        if n_missing:
            print(f"         ({n_missing} block(s) have no data in the forcing)")

    for storm in STORMS:
        sub = runs[(runs["storm"] == storm) & (runs["flow_outside_window_max_m3s"] > 0)]
        print("\n" + "-" * 96)
        print(f"2. {storm}")
        print("-" * 96)
        if sub.empty:
            print("  none of this storm's volume-matched runs has flow outside the plotted window; nothing to look at.")
            continue
        row = sub.sort_values("flow_outside_window_max_m3s", ascending=False).iloc[0]
        run_id = row["run_id"]
        d = read_run(csv_dir, run_id)
        print(f"  {len(sub)} run(s) have late flow ({', '.join(sorted(sub['f_group'].unique()))} f). Looking at the one with the most:")
        print(f"  {run_id}   f {row['f_RS_abs']:.4f}, Ks {row['Ks_mult']:.2f}; late flow up to "
              f"{row['flow_outside_window_max_m3s']:.2f} m3/s at {row['flow_outside_window_time']}")
        if d is None:
            print(f"  (its hydrograph file {run_id}{RUN_FILE_SUFFIX} was not found in {csv_dir}; skipped)")
            continue
        sim = d["Simulated"].groupby(d.index.floor(BLOCK)).max().reindex(idx)
        tru = d["Observed"].groupby(d.index.floor(BLOCK)).max().reindex(idx)
        print(f"\n  {'block':<6} {'SMF rain':>9} {'SMPHQ rain':>11} {'truth flow':>11} {'run flow':>9}")
        print(f"  {'start':<6} {'mm/hr':>9} {'mm/hr':>11} {'m3/s':>11} {'m3/s':>9}")
        skipped = False
        for ts in idx:
            r1, r2, q0, q1 = rain["SMF"][ts], rain["SMPHQ"][ts], tru[ts], sim[ts]
            active = (r1 > 0) or (r2 > 0) or (q0 >= 0.01) or (q1 >= 0.01)
            if not active:
                if not skipped:
                    print("  ...    (quiet blocks skipped: no rain, truth and run both under 0.01 m3/s)")
                    skipped = True
                continue
            skipped = False
            f = lambda v, nd=1: "n/a" if pd.isna(v) else f"{v:.{nd}f}"
            print(f"  {hhmm(ts):<6} {f(r1):>9} {f(r2):>11} {f(q0, 2):>11} {f(q1, 2):>9}")
        # ---- where the late rise starts and peaks, from the run's own 5-minute flow ----
        mid = pd.Timestamp("2014-08-13 00:00")
        diff = (d["Simulated"] - d["Observed"]).to_numpy(float)       # run minus truth, every 5 minutes
        tix = d.index
        after = np.flatnonzero(np.asarray(tix >= pd.Timestamp(LATE_FROM)))
        if len(after) == 0:
            print(f"\n  The hydrograph file has no data after {LATE_FROM}; cannot look for the late rise.")
            continue
        ip = int(after[int(np.argmax(diff[after]))])                    # where the run is furthest above the truth
        lo = max(0, ip - int(round(FOOT_LOOKBACK_MIN / STEP_MIN)))
        seg = diff[lo:ip + 1]
        i0 = lo + int(np.flatnonzero(seg == seg.min())[-1])             # the lowest point just before the peak
        sim_v, tru_v = d["Simulated"].to_numpy(float), d["Observed"].to_numpy(float)
        print(f"\n  Main late rise: it starts at {hhmm(tix[i0])} (run {sim_v[i0]:.2f} m3/s, truth {tru_v[i0]:.2f}) and peaks at "
              f"{hhmm(tix[ip])} (run {sim_v[ip]:.2f} m3/s, truth {tru_v[ip]:.2f}; the run is {diff[ip]:.2f} above the truth).")
        step_idx = np.flatnonzero(diff[i0:ip + 1] - diff[i0] >= STEP_UP_M3S)
        t_step = tix[i0 + int(step_idx[0])] if len(step_idx) else None
        for g, s_g in rain.items():
            before = s_g[(s_g.index < tix[i0]) & (s_g.index >= tix[i0] - pd.Timedelta(minutes=60))]
            print(f"  {g:<6} rain in the 60 min before the rise starts: {float(np.nansum(before.to_numpy())) * BLOCK_HR:.1f} mm")
        tot = (rain["SMF"].fillna(0.0) + rain["SMPHQ"].fillna(0.0))
        tot = tot[tot.index >= mid]
        if len(tot) and float(tot.max()) > 0:
            blk = tot.idxmax()
            print(f"  Strongest rain block after midnight (both gauges added): the block starting {hhmm(blk)} "
                  f"(SMF {rain['SMF'][blk]:.1f}, SMPHQ {rain['SMPHQ'][blk]:.1f} mm/hr; unscaled record).")
            if t_step is not None:
                lead = (t_step - blk).total_seconds() / 60.0
                when = f"{lead:.0f} min after" if lead > 0 else ("when" if lead == 0 else f"{-lead:.0f} min before")
                print(f"  The run gets {STEP_UP_M3S:.2f} m3/s above its starting level at {hhmm(t_step)}, {when} that block begins; it peaks "
                      f"{(tix[ip] - blk).total_seconds() / 60.0:.0f} min after the block begins.")
            else:
                print(f"  The rise never gets {STEP_UP_M3S:.2f} m3/s above its starting level, so no timing is given.")
        else:
            print("  There is no rain after midnight in the window, so the rise cannot be compared with a rain block.")
        day = np.asarray(tix >= mid)
        extra = float(diff[day].sum()) * STEP_MIN * 60.0
        vol = float(tru_v.sum()) * STEP_MIN * 60.0
        print(f"  Extra water the run puts out after midnight compared with the truth: {extra:.0f} m3 = {100.0 * extra / vol:.1f}% of the truth's volume.")

    print("\n" + "=" * 96)
    if any_rain_after_storm:
        print("There is some rain after midnight in the forcing (see section 1 and the tables). Compare its times with the start of the late flow.")
    else:
        print("The forcing has no rain after midnight in this window, so the late flow cannot be a response to rain in the record.")
    print("This looks at two gauge points only; it says whether rain was there and when, not what the water did inside the model.")
    print("Only the one run per storm with the biggest late rise is shown. plot_series120_late_rise.py draws all the volume-matched runs.")


if __name__ == "__main__":
    main()
