#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
wy2020_slide_data.py   (v1.0)
=============================
READ-ONLY. Pulls the small tables needed for the "Original vs Re-tuned, water year 2020" slide out of result files
that already exist. It runs no model, changes no input file and overwrites nothing.

THE SPAN. Oct 1 2019 00:00 to May 31 2020 23:00 (244 days), the span that was scored. The WY2020 runs start on Sep 1 2019 as spin-up; September is left out.

WHAT IT READS (from workspaces/Sycamore/lab/, or change the folders with --results-dir / --obs-file)
  results/Sycamore_WY2020_Original_Outlet.qout          Original, outlet flow         (Time_hr 0 = Sep 1 2019 00:00)
  results/Sycamore_WY2020_Retuned_Outlet.qout           Re-tuned, outlet flow
  results/Sycamore_WY2020_<run>6575_00.rft / .mrf       runoff made on the land, by type, and the basin rain (both runs). Found by name
                                                        pattern. If a file is missing the script still writes the flow table and says what it skipped.
  ../init_data/obs/sycamore_09510200_obs_flow_WY2020.csv   USGS gauge 09510200 (datetime, q_cms)

WHAT IT WRITES (a new folder, wy2020_slide_data/ by default)
  wy2020_daily_flow.csv      one row per day, Oct 1 2019 - May 31 2020: daily mean and daily highest hourly flow (m3/s) for the gauge,
                             Original and Re-tuned, how many gauge hours each day has, and the basin rain (mm per day)
  wy2020_runoff_by_type.csv  one row per run and month (Oct-May), plus a "span" row: Hsrf, Sbsrf, Psrf, Satsrf and their total (million m3)
  wy2020_checks.txt          the checks printed below, saved so you can attach them

RUN IT from workspaces/Sycamore/lab/:
    python wy2020_slide_data.py

WORDS
  Hsrf   infiltration excess     Sbsrf  saturation excess     Psrf  perched return flow     Satsrf  groundwater exfiltration
  runoff made   runoff GENERATED on the land, all cells added, BEFORE routing (not flow at the outlet)
  hourly mean   the model and gauge files are averaged into hours first, as in the project's other scripts
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

T0 = pd.Timestamp("2019-09-01 00:00")        # Time_hr 0 of every WY2020 run
RUN_HOURS = 6575
START = pd.Timestamp("2019-10-01 00:00")
END = pd.Timestamp("2020-06-01 00:00")
RUNS = [("Original", "Sycamore_WY2020_Original"),
        ("Re-tuned", "Sycamore_WY2020_Retuned")]
TYPES = ["Hsrf", "Sbsrf", "Psrf", "Satsrf"]
AFTER_DAYS = 3
STORMS = [("Nov 19-22", "2019-11-19", "2019-11-22"), ("Nov 28-30", "2019-11-28", "2019-11-30"), ("Dec 23-28", "2019-12-23", "2019-12-28"),
          ("Feb 10-11", "2020-02-10", "2020-02-11"), ("Feb 22-23", "2020-02-22", "2020-02-23"),
          ("Mar 10-13", "2020-03-10", "2020-03-13"), ("Mar 18-19", "2020-03-18", "2020-03-19")]
# Scored storm-window volumes (M m3) on record in the WY2020 write-up: gauge, Re-tuned, Original
EXPECTED = {"Nov 19-22": (0.015, 2.800, 3.165), "Nov 28-30": (1.473, 3.235, 3.538), "Dec 23-28": (4.069, 0.218, 0.154),
            "Feb 10-11": (0.386, 1.767, 1.765), "Feb 22-23": (4.131, 4.839, 4.981), "Mar 10-13": (5.626, 0.875, 0.493),
            "Mar 18-19": (3.317, 1.031, 0.841)}
# Re-tuned runoff made over the seven scored storms, share by type, as printed in the deck's mechanism slide (percent)
DECK_SHARES = {"Hsrf": 88.9, "Sbsrf": 7.2, "Satsrf": 3.9}

lines = []


def log(msg=""):
    print(msg, flush=True)
    lines.append(msg)


def find_obs(path):
    p = Path(path)
    if p.is_file():
        return p
    hits = sorted(Path(".").rglob(p.name))
    if hits:
        log(f"  (did not find {p}; using {hits[0]})")
        return hits[0]
    sys.exit(f"STOPPED: cannot find the gauge file {p.name}\nRun this from workspaces/Sycamore/lab/, or use --obs-file.")


def find_exact(folder, name):
    p = Path(folder) / name
    if p.is_file():
        return p
    hits = sorted(Path(".").rglob(name))
    if hits:
        log(f"  (did not find {p}; using {hits[0]})")
        return hits[0]
    sys.exit(f"STOPPED: cannot find {name}\nRun this from workspaces/Sycamore/lab/, or use --results-dir to point at the folder that holds it.")


def find_pattern(folder, stem, ext):
    """results/<stem><run hours...>.<ext>: the stem must be followed by a digit, so the channel-loss variants are not picked up."""
    pat = re.compile("^" + re.escape(stem) + r"\d.*\." + ext + "$")
    hits = sorted(p for p in Path(folder).glob(stem + "*." + ext) if pat.match(p.name))
    if not hits:
        hits = sorted(p for p in Path(".").rglob(stem + "*." + ext) if pat.match(p.name))
    if len(hits) > 1:
        pref = [p for p in hits if str(RUN_HOURS) in p.name]
        if len(pref) == 1:
            hits = pref
    return hits


def read_qout(path):
    raw = pd.read_csv(path).iloc[:, :2]
    raw.columns = ["t", "Q"]
    raw = raw.apply(pd.to_numeric, errors="coerce").dropna()
    s = pd.Series(raw.Q.to_numpy(float), index=T0 + pd.to_timedelta(raw.t.to_numpy(), unit="h"))
    return s[~s.index.duplicated(keep="last")].sort_index()


def read_hourly_table(path, need):
    d = pd.read_csv(path)
    missing = [c for c in need if c not in d.columns]
    if missing:
        sys.exit(f"STOPPED: {path} has no column(s) {missing}. Columns found: {list(d.columns)}")
    d = d[need].apply(pd.to_numeric, errors="coerce").dropna()
    big = (d[need[1:]] >= 9000).any(axis=1) | (d["Time_hr"] > RUN_HOURS + 0.5)    # an unfilled last row is flagged 9999 in these files
    d = d[~big].copy()
    d.index = T0 + pd.to_timedelta(d["Time_hr"].to_numpy(float), unit="h")      # label = END of the hour
    return d[~d.index.duplicated(keep="last")].sort_index()


def hourly(q, a, b):
    s = q[(q.index >= pd.Timestamp(a)) & (q.index < pd.Timestamp(b))]
    return s.resample("60min").mean().dropna()


def vol(h):
    return float(h.sum() * 3600.0 / 1e6) if len(h) else float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--obs-file", default="../init_data/obs/sycamore_09510200_obs_flow_WY2020.csv")
    ap.add_argument("--out-dir", default="wy2020_slide_data")
    a = ap.parse_args()
    out = Path(a.out_dir)
    if out.exists() and any(out.iterdir()):
        sys.exit(f"STOPPED: {out}/ already has files in it; this script does not overwrite. Use --out-dir with a new name.")
    out.mkdir(parents=True, exist_ok=True)

    log("wy2020_slide_data.py v1.0 (read-only)")
    gpath = find_obs(a.obs_file)
    g = pd.read_csv(gpath, parse_dates=["datetime"]).set_index("datetime")["q_cms"]
    log(f"  gauge file: {gpath}")
    q, mrf, rft, skipped = {}, {}, {}, []
    for lab, stem in RUNS:
        q[lab] = read_qout(find_exact(a.results_dir, f"{stem}_Outlet.qout"))
        last_h = (q[lab].index.max() - T0) / pd.Timedelta(hours=1)
        log(f"  {lab}: outlet file reaches model hour {last_h:.1f} of {RUN_HOURS}")
        for ext, store, need in (("mrf", mrf, ["Time_hr", "MAP_mm_hr", "Srf_m3_s"]), ("rft", rft, ["Time_hr"] + [f"{t}_m3_s" for t in TYPES])):
            hits = find_pattern(a.results_dir, stem, ext)
            if len(hits) == 1:
                store[lab] = read_hourly_table(hits[0], need)
                log(f"  {lab} .{ext}: {hits[0].name} ({len(store[lab])} hourly rows)")
            elif not hits:
                skipped.append(f"{lab} .{ext}: NOT FOUND (looked for {stem}<digits>*.{ext} in {a.results_dir} and below)")
            else:
                skipped.append(f"{lab} .{ext}: MORE THAN ONE file fits: {', '.join(p.name for p in hits)}")
    for s_ in skipped:
        log("  SKIPPED " + s_)
    have_rft = all(l in rft for l, _ in RUNS)
    have_mrf = all(l in mrf for l, _ in RUNS)

    # ---- 1. daily flow ---------------------------------------------------------------------------------
    days = pd.date_range(START, END - pd.Timedelta(days=1), freq="D")
    hg, h1, h2 = hourly(g, START, END), hourly(q["Original"], START, END), hourly(q["Re-tuned"], START, END)
    rain_h = mrf["Original"]["MAP_mm_hr"] if have_mrf else None
    rows = []
    for d0 in days:
        d1 = d0 + pd.Timedelta(days=1)
        sel = lambda h: h[(h.index >= d0) & (h.index < d1)]
        gd, od, rd = sel(hg), sel(h1), sel(h2)
        rn = rain_h[(rain_h.index > d0) & (rain_h.index <= d1)].sum() if rain_h is not None else np.nan   # hour-END labels
        rows.append(dict(date=d0.strftime("%Y-%m-%d"),
                         gauge_hours=len(gd), gauge_mean_m3s=gd.mean() if len(gd) else np.nan, gauge_max_m3s=gd.max() if len(gd) else np.nan,
                         original_mean_m3s=od.mean(), original_max_m3s=od.max(),
                         retuned_mean_m3s=rd.mean(), retuned_max_m3s=rd.max(), rain_mm=rn))
    D = pd.DataFrame(rows)
    D.to_csv(out / "wy2020_daily_flow.csv", index=False, float_format="%.5f")

    # ---- 2. runoff made by type, by month --------------------------------------------------------------
    months = pd.period_range("2019-10", "2020-05", freq="M")
    T = None
    if have_rft:
        trows = []
        for lab, _ in RUNS:
            d = rft[lab]
            hs = d.copy()
            hs.index = hs.index - pd.Timedelta(hours=1)                     # label by the START of each hour (calendar months)
            hs = hs[(hs.index >= START) & (hs.index < END)]
            per = hs.index.to_period("M")
            sums = {t: hs[f"{t}_m3_s"].groupby(per).sum() * 3600 / 1e6 for t in TYPES}
            for p in months:
                v = {t: float(sums[t].get(p, 0.0)) for t in TYPES}
                trows.append(dict(run=lab, month=str(p), **{t: v[t] for t in TYPES}, total=sum(v.values())))
            yv = {t: float(hs[f"{t}_m3_s"].sum() * 3600 / 1e6) for t in TYPES}
            trows.append(dict(run=lab, month="span", **yv, total=sum(yv.values())))
        T = pd.DataFrame(trows)
        T.to_csv(out / "wy2020_runoff_by_type.csv", index=False, float_format="%.5f")

    # ---- 3. checks -------------------------------------------------------------------------------------
    log("\nCHECKS (the 'on record' numbers come from the WY2020 write-up and the deck; small differences in the last digit are fine)")
    log("  A. Volumes in the seven scored storm windows, M m3 (window = first rain day 00:00 to 4 days after the last rain day)")
    log(f"     {'storm':10s} {'gauge':>8s} {'record':>7s}   {'Re-tuned':>8s} {'record':>7s}   {'Original':>8s} {'record':>7s}")
    hg_all, h1_all, h2_all = hourly(g, T0, END), hourly(q["Original"], T0, END), hourly(q["Re-tuned"], T0, END)
    bad = 0
    sums7 = {"gauge": 0.0, "Re-tuned": 0.0, "Original": 0.0}
    for name, f, l in STORMS:
        lo, hi = pd.Timestamp(f), pd.Timestamp(l) + pd.Timedelta(days=1 + AFTER_DAYS)
        vg, vr, vo = vol(hourly(g, lo, hi)), vol(hourly(q["Re-tuned"], lo, hi)), vol(hourly(q["Original"], lo, hi))
        eg, er, eo = EXPECTED[name]
        ok = all(abs(x - y) <= 0.005 for x, y in ((vg, eg), (vr, er), (vo, eo)))
        bad += 0 if ok else 1
        sums7["gauge"] += vg; sums7["Re-tuned"] += vr; sums7["Original"] += vo
        log(f"     {name:10s} {vg:8.3f} {eg:7.3f}   {vr:8.3f} {er:7.3f}   {vo:8.3f} {eo:7.3f}   {'ok' if ok else 'DIFFERS'}")
    log(f"     seven-storm totals here: gauge {sums7['gauge']:.2f}, Re-tuned {sums7['Re-tuned']:.2f}, Original {sums7['Original']:.2f}   (on record: gauge 19.02, Re-tuned 14.76)")
    if have_rft:
        log("  B. Re-tuned runoff made over the same seven storm windows, share by type (the deck's mechanism slide: infiltration excess 88.9, saturation excess 7.2, groundwater exfiltration 3.9)")
        for lab, _ in RUNS:
            tot = {t: 0.0 for t in TYPES}
            for name, f, l in STORMS:
                lo, hi = pd.Timestamp(f), pd.Timestamp(l) + pd.Timedelta(days=1 + AFTER_DAYS)
                for t in TYPES:
                    tot[t] += vol(hourly(rft[lab][f"{t}_m3_s"], lo, hi))
            s_ = sum(tot.values())
            sh = {t: 100 * tot[t] / s_ for t in TYPES}
            msg = f"     {lab:9s} Hsrf {sh['Hsrf']:5.1f}%  Sbsrf {sh['Sbsrf']:5.1f}%  Psrf {sh['Psrf']:5.1f}%  Satsrf {sh['Satsrf']:5.1f}%   (total {s_:.3f} M m3)"
            if lab == "Re-tuned":
                ok_ = all(abs(sh[t] - DECK_SHARES[t]) <= 0.15 for t in DECK_SHARES)
                msg += "   " + ("ok" if ok_ else "DIFFERS")
            log(msg)
    if have_rft and have_mrf:
        log("  C. The four types add up to the .mrf total runoff (Srf), whole run, M m3")
        for lab, _ in RUNS:
            tot = float(sum(rft[lab][f"{t}_m3_s"].sum() for t in TYPES) * 3600 / 1e6)
            srf = float(mrf[lab]["Srf_m3_s"].sum() * 3600 / 1e6)
            log(f"     {lab:9s} types {tot:8.3f}   Srf {srf:8.3f}   {'ok' if abs(tot - srf) <= 0.02 else 'DIFFERS'}")
    log(f"\n  {bad} of 7 storm windows differ from the record.")

    log("\nMONTHLY OUTLET VOLUME, M m3 (sums of hourly mean flows). The gauge counts only hours that have data.")
    log(f"  {'month':9s} {'gauge':>8s} {'Original':>9s} {'Re-tuned':>9s}   gauge hours / hours in month")
    for p in months:
        a0, b0 = p.start_time, p.end_time.ceil("D")
        hgm = hourly(g, a0, b0)
        log(f"  {str(p):9s} {vol(hgm):8.3f} {vol(hourly(q['Original'], a0, b0)):9.3f} {vol(hourly(q['Re-tuned'], a0, b0)):9.3f}   {len(hgm)} / {int((b0 - a0) / pd.Timedelta(hours=1))}")
    log(f"  {'Oct-May':9s} {vol(hg):8.3f} {vol(h1):9.3f} {vol(h2):9.3f}")
    if T is not None:
        log("\nRUNOFF MADE ON THE LAND, Oct 1 2019 - May 31 2020, M m3 (before routing)")
        log(f"  {'run':9s} {'Hsrf':>8s} {'Sbsrf':>8s} {'Psrf':>8s} {'Satsrf':>8s} {'total':>8s}")
        for lab, _ in RUNS:
            r = T[(T.run == lab) & (T.month == "span")].iloc[0]
            log(f"  {lab:9s} {r.Hsrf:8.3f} {r.Sbsrf:8.3f} {r.Psrf:8.3f} {r.Satsrf:8.3f} {r.total:8.3f}")
    else:
        log("\nRUNOFF MADE BY TYPE: skipped, because a .rft file was not found (see the SKIPPED lines at the top). The flow table was still written.")
    (out / "wy2020_checks.txt").write_text("\n".join(lines) + "\n")
    log(f"\nWrote {out}/wy2020_daily_flow.csv, " + (f"{out}/wy2020_runoff_by_type.csv, " if T is not None else "") + f"{out}/wy2020_checks.txt")
    log("Next step: attach those files to the chat.")


if __name__ == "__main__":
    main()
