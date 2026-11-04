#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
model_state_analysis_wy2017.py   (v1.0, 2026-11-03)
====================================================
READ-ONLY look at the four finished WY2017 whole-year runs (W5), storm by storm. It runs no model and changes no input file.
It reads, from ONE folder (--in-dir):
    <run>8759_00.mrf                                  the basin-average results tRIBS writes each hour (rain, runoff, water table, soil moisture, channel loss)
    <run>_Outlet.qout                                 routed outlet flow
    hourly_basin_wy2017.csv                           radar basin-mean rain per hour (from storm_features_sycamore.py)
    sycamore_09510200_obs_flow_WY2017.csv             the USGS 09510200 gauge file
and writes, into --out-dir (default model_state_out/):
    state_at_storm_start.csv      per storm and run: water table depth, soil moisture, channel loss to date, outlet flow at the start of the storm
    storm_windows.csv             storm volumes (gauge and each run) for windows that end 0, 1, 3, 7 days after the last rain day
                                  (the 7-day windows can run into the next storm, so read them as a sensitivity check, not as storm volumes)
    events_5mm_winter.csv         EVERY winter rain event of 5 mm or more (the 20 mm rule keeps only 5 of them): rain, intensity, season-to-date rain,
                                  gauge volume, model volume, runoff ratios
    fig_daily_gauge_vs_model_wy2017.png   daily basin rain and daily outlet volume (log scale), gauge vs model, Nov 1 - Mar 31
and prints the checks (does the model's own basin rain equal the radar basin mean? does the channel loss add up to the old report?).

WORDS
  runoff ratio   volume at the outlet divided by the rain that fell on the basin (both as mm over the basin; 1 mm over 424.6 km2 = 0.4246 M m3)
  window         the stretch of time a storm's volume is added up over; here it starts 00:00 on the first rain day and ends k days after the last rain day
  tail           the slow part of the recession (the flow that keeps coming for days after the peak)
  season-to-date rain   basin rain since Oct 1, a plain index of how wet the basin already was

USAGE (from the folder that holds the ten files, or give --in-dir)
  python model_state_analysis_wy2017.py --in-dir . --out-dir model_state_out
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

VERSION = 1
START = pd.Timestamp("2016-10-01 00:00")
AREA_KM2 = 424.6
MM_TO_MM3 = AREA_KM2 / 1000.0                  # 1 mm over the basin = 0.4246 M m3
RUNS = [("re-tuned", "W5_d012_k10_f0p15_as300_ksx1p03_nocl"),
        ("as drawn, loss off", "W5_d012_k10_base_nocl"),
        ("as drawn, Kchan 0.5", "W5_d012_k10_base_kchan0p5"),
        ("as drawn, Kchan 0.2", "W5_d012_k10_base_kchan0p2")]
STORMS = [("Dec 21-25", "2016-12-21", "2016-12-25"), ("Dec 30-Jan 2", "2016-12-30", "2017-01-02"), ("Jan 19-21", "2017-01-19", "2017-01-21"),
          ("Feb 18-19", "2017-02-18", "2017-02-19"), ("Feb 27-28", "2017-02-27", "2017-02-28")]
NEAR = [("Dec 15-17", "2016-12-15", "2016-12-17"), ("Jan 23-24", "2017-01-23", "2017-01-24")]
GAUGE_GRAY, INK, INK_2, GRID = "#8a8984", "#0b0b0b", "#52514e", "#d9d8d3"
RAIN_BLUE, RETUNED_COLOR, ASDRAWN_COLOR = "#8db4e2", "#6b3fa0", "#2a78d6"


def log(msg: str = "") -> None:
    print(msg, flush=True)


def read_qout(path: Path) -> pd.Series:
    raw = pd.read_csv(path).iloc[:, :2]
    raw.columns = ["t", "Q"]
    s = pd.Series(raw.Q.to_numpy(float), index=START + pd.to_timedelta(raw.t.to_numpy(), unit="h"))
    return s[~s.index.duplicated(keep="last")].sort_index()


def read_mrf(path: Path) -> pd.DataFrame:
    d = pd.read_csv(path)
    d = d[d["Time_hr"] <= 8759.5].copy()                 # the last row (hour 8760) is a filler row (RainMin 9999.99)
    d.index = START + pd.to_timedelta(d["Time_hr"].to_numpy(), unit="h")   # label = END of the hour
    return d


def hourly(q: pd.Series, a, b) -> pd.Series:
    s = q[(q.index >= pd.Timestamp(a)) & (q.index < pd.Timestamp(b))]
    return s.resample("60min").mean().dropna()


def vol(h: pd.Series) -> float:
    return float(h.sum() * 3600.0 / 1e6) if len(h) else np.nan


def find_events(daily: pd.Series, min_day=0.5, min_total=5.0):
    ev, run = [], []
    for d, v in daily.items():
        if v >= min_day:
            run.append((d, v))
        else:
            if run and sum(x for _, x in run) >= min_total:
                ev.append((run[0][0], run[-1][0], sum(x for _, x in run)))
            run = []
    if run and sum(x for _, x in run) >= min_total:
        ev.append((run[0][0], run[-1][0], sum(x for _, x in run)))
    return ev


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--in-dir", default=".")
    p.add_argument("--out-dir", default="model_state_out")
    p.add_argument("--no-figure", action="store_true")
    a = p.parse_args(argv)
    src, out = Path(a.in_dir).expanduser(), Path(a.out_dir).expanduser()
    out.mkdir(parents=True, exist_ok=True)
    mrf, qout = {}, {}
    for lab, rid in RUNS:
        mrf[lab] = read_mrf(src / f"{rid}8759_00.mrf")
        qout[lab] = read_qout(src / f"{rid}_Outlet.qout")
    gfile = pd.read_csv(src / "sycamore_09510200_obs_flow_WY2017.csv", parse_dates=["datetime"]).set_index("datetime")
    gq = gfile["q_cms"]
    hb = pd.read_csv(src / "hourly_basin_wy2017.csv", parse_dates=["label_end_mst", "interval_start_mst"]).set_index("interval_start_mst")
    first = RUNS[0][0]

    # ---- 1. checks ----------------------------------------------------------------------------------
    log("CHECK 1: the model's own basin rain (.mrf MAP_mm_hr) vs the radar basin mean from storm_features_sycamore.py")
    m = mrf[first]["MAP_mm_hr"]
    r = hb["basin_mean_mm"].copy()
    r.index = r.index + pd.Timedelta(hours=1)           # label by hour END, like the .mrf
    r = r.reindex(m.index)
    d = (m - r).abs()
    log(f"  year total: model {m.sum():.1f} mm, radar {r.sum():.1f} mm; hourly correlation {m.corr(r):.4f}; largest hourly difference {d.max():.3f} mm")
    same = all(np.allclose(mrf[k]["MAP_mm_hr"], m) for k in mrf)
    log(f"  the same rain went into all four runs: {same}")
    log("CHECK 2: channel loss (.mrf ChannelPerc_m3 is m3 PER HOUR, not a running total) and runoff generated vs delivered")
    for lab, _ in RUNS:
        gen = mrf[lab]["Srf_m3_s"].sum() * 3600 / 1e6
        dlv = vol(hourly(qout[lab], START, "2017-10-01"))
        log(f"  {lab:22s} channel loss {mrf[lab]['ChannelPerc_m3'].sum() / 1e6:6.3f} M m3;  runoff generated {gen:6.2f};  delivered at the outlet {dlv:6.2f} M m3")

    # ---- 2. volumes at the standard window and the sensitivity to the window ------------------------
    S = sorted(STORMS + NEAR, key=lambda x: x[1])
    rows = []
    for k in (0, 1, 3, 7):
        for n, f, l in S:
            a0, b0 = pd.Timestamp(f), pd.Timestamp(l) + pd.Timedelta(days=1 + k)
            gh = hourly(gq, a0, b0)
            row = dict(storm=n, kind="listed" if (n, f, l) in STORMS else "near-miss", k_days_after=k, window_start=a0, window_end=b0,
                       gauge_hours=len(gh), hours_expected=int((b0 - a0) / pd.Timedelta(hours=1)), gauge_Mm3=vol(gh))
            for lab, _ in RUNS:
                row[lab + " Mm3"] = vol(hourly(qout[lab], a0, b0))
            rows.append(row)
    W = pd.DataFrame(rows)
    W.to_csv(out / "storm_windows.csv", index=False, float_format="%.4f")
    log("\nVOLUMES by storm and window (M m3; re-tuned run / gauge). k = days after the last rain day; k = 3 is the window used so far and in the pre-registration")
    piv = W.assign(txt=lambda x: x[first + " Mm3"].map("{:.2f}".format) + " / " + x["gauge_Mm3"].map("{:.2f}".format)).pivot(index="storm", columns="k_days_after", values="txt")
    log(piv.reindex([n for n, _, _ in S]).to_string())
    a0, b0 = pd.Timestamp("2016-12-30"), pd.Timestamp("2017-01-06")
    gh = hourly(gq, a0, b0)
    log("\nLIKE FOR LIKE, Dec 30-Jan 2 (the gauge has no data Jan 1 00:00 - Jan 2 00:15): model counted only over the gauge's hours")
    for lab, _ in RUNS[:2]:
        mh = hourly(qout[lab], a0, b0)
        log(f"  {lab:22s} gauge {vol(gh):.3f} M m3 over {len(gh)} h;  model over the same hours {vol(mh.reindex(gh.index)):.3f};  ratio {vol(mh.reindex(gh.index)) / vol(gh):.2f}  (model over all 168 h: {vol(mh):.3f})")

    # ---- 3. state at storm start ---------------------------------------------------------------------
    srows = []
    for n, f, l in S:
        t0 = pd.Timestamp(f)
        for lab, _ in RUNS:
            mm = mrf[lab]
            s0 = mm.loc[t0]
            srows.append(dict(storm=n, run=lab, water_table_depth_mm=s0["MDGW_mm"], soil_moisture_top100_frac=s0["MSM100_[]"], soil_moisture_root_frac=s0["MSMRt_[]"],
                              saturated_area_pct=s0["SatPercent_[]"], channel_loss_to_date_Mm3=mm["ChannelPerc_m3"][mm.index <= t0].sum() / 1e6,
                              runoff_generated_prev24h_Mm3=mm["Srf_m3_s"][(mm.index > t0 - pd.Timedelta(hours=24)) & (mm.index <= t0)].sum() * 3600 / 1e6,
                              outlet_flow_at_start_m3s=float(qout[lab][qout[lab].index <= t0].iloc[-1])))
    ST = pd.DataFrame(srows)
    ST.to_csv(out / "state_at_storm_start.csv", index=False, float_format="%.4f")
    log("\nMODEL STATE AT THE START OF EACH STORM (00:00 on the first rain day)")
    for col, title in (("water_table_depth_mm", "water table depth, mm below the surface (smaller = shallower)"), ("soil_moisture_root_frac", "root-zone soil moisture (fraction)"),
                       ("channel_loss_to_date_Mm3", "channel loss so far, M m3"), ("outlet_flow_at_start_m3s", "outlet flow at the start, m3/s")):
        log(f"  {title}")
        log(ST.pivot(index="storm", columns="run", values=col).reindex([n for n, _, _ in S])[[x for x, _ in RUNS]].round(3).to_string())

    # ---- 4. every winter rain event of 5 mm or more ---------------------------------------------------
    dr = hb["basin_mean_mm"].resample("D").sum()
    dr = dr[(dr.index >= "2016-10-01") & (dr.index < "2017-06-01")]
    ev = find_events(dr)
    listed_starts = {pd.Timestamp(f) for _, f, _ in STORMS}
    erows = []
    for i, (f, l, tot) in enumerate(ev):
        b0 = l + pd.Timedelta(days=4)
        if i + 1 < len(ev):
            b0 = min(b0, ev[i + 1][0])                    # a window never runs into the next event
        sub = hb[(hb.index >= f) & (hb.index < l + pd.Timedelta(days=1))]
        season = hb["basin_mean_mm"][(hb.index >= START) & (hb.index < f)].sum()
        gh = hourly(gq, f, b0)
        row = dict(event=f"{f:%b} {f.day}" + ("" if f == l else f"-{l.day}" if f.month == l.month else f"-{l:%b} {l.day}"), start=f, on_20mm_list=f in listed_starts,
                   basin_mm=sub["basin_mean_mm"].sum(), peak_1h_mm=sub["basin_mean_mm"].max(), excess_4_mm=sub["basin_excess_4_mm"].sum(),
                   excess_2p2_mm=sub["basin_excess_2p2_mm"].sum(), season_to_date_mm=season, window_days=(b0 - f).days,
                   gauge_hours=f"{len(gh)}/{int((b0 - f) / pd.Timedelta(hours=1))}", gauge_Mm3=vol(gh),
                   water_table_depth_mm=mrf[first]["MDGW_mm"].loc[f])
        for lab, _ in RUNS:
            row[lab + " Mm3"] = vol(hourly(qout[lab], f, b0))
        row["gauge_runoff_ratio_pct"] = 100 * row["gauge_Mm3"] / (MM_TO_MM3 * row["basin_mm"])
        row["model_runoff_ratio_pct (re-tuned)"] = 100 * row[first + " Mm3"] / (MM_TO_MM3 * row["basin_mm"])
        erows.append(row)
    EV = pd.DataFrame(erows)
    EV.to_csv(out / "events_5mm_winter.csv", index=False, float_format="%.4f")
    log("\nEVERY WINTER RAIN EVENT OF 5 MM OR MORE (Oct 1 - May 31; windows end 3 days after the last rain day or at the next event, whichever comes first)")
    show = EV[["event", "on_20mm_list", "basin_mm", "excess_4_mm", "season_to_date_mm", "gauge_Mm3", first + " Mm3", "gauge_runoff_ratio_pct", "model_runoff_ratio_pct (re-tuned)"]]
    log(show.round(2).to_string(index=False))
    nl = EV[~EV.on_20mm_list]
    log(f"  not on the 20 mm list: {len(nl)} events, {nl.basin_mm.sum():.0f} mm of rain; gauge {nl.gauge_Mm3.sum():.2f} M m3, re-tuned model {nl[first + ' Mm3'].sum():.2f} M m3")

    # ---- 5. figure -------------------------------------------------------------------------------------
    if not a.no_figure:
        make_figure(out / "fig_daily_gauge_vs_model_wy2017.png", hb, gq, qout, mrf)
        log(f"\nFigure: {out / 'fig_daily_gauge_vs_model_wy2017.png'}")
    log(f"Wrote to {out}/ : state_at_storm_start.csv, storm_windows.csv, events_5mm_winter.csv")
    return 0


def make_figure(path: Path, hb, gq, qout, mrf) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    days = pd.date_range("2016-11-01", "2017-03-31", freq="D")
    rain = hb["basin_mean_mm"].resample("D").sum().reindex(days)
    floor = 0.002

    def daily(q):
        return pd.Series([vol(hourly(q, d0, d0 + pd.Timedelta(days=1))) for d0 in days], index=days)

    g, r1, r2 = daily(gq), daily(qout["re-tuned"]), daily(qout["as drawn, loss off"])
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(12.5, 7.4), sharex=True, gridspec_kw={"height_ratios": [1, 2.6], "hspace": 0.06})
    ax0.bar(days, rain, width=0.85, color=RAIN_BLUE, lw=0)
    ax0.set_ylabel("Basin rain\n(mm per day)")
    for n, f, l in STORMS:
        a0, b0 = pd.Timestamp(f), pd.Timestamp(l) + pd.Timedelta(days=4)
        for ax in (ax0, ax1):
            ax.axvspan(a0, b0, color="#efeeea", zorder=0, lw=0)
        row_y = 66 if n in ("Dec 30-Jan 2", "Feb 27-28") else 56         # two rows so neighboring labels do not collide
        ax0.text(a0 + (b0 - a0) / 2, row_y, n, ha="center", va="top", fontsize=8.5, color=INK_2)
    ax0.set_ylim(0, 68)
    for ax in (ax0, ax1):
        ax.grid(True, axis="y", color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        for s_ in ("top", "right"):
            ax.spines[s_].set_visible(False)
        ax.tick_params(labelsize=9)
    ax1.plot(days, g.clip(lower=floor), color=GAUGE_GRAY, lw=2.6, solid_capstyle="round", zorder=3)
    ax1.plot(days, r2.clip(lower=floor), color=ASDRAWN_COLOR, lw=1.2, zorder=4)
    ax1.plot(days, r1.clip(lower=floor), color=RETUNED_COLOR, lw=1.2, zorder=5)
    ax1.set_yscale("log")
    ax1.set_ylim(floor, 12)
    ax1.set_ylabel("Outlet volume per day (M m$^3$, log scale)")
    ax1.text(pd.Timestamp("2016-12-09"), floor * 1.6, "axis bottom = no flow", fontsize=8, color=INK_2, va="bottom")
    ax1.annotate("no gauge\ndata Jan 1", (pd.Timestamp("2016-12-27"), 1.3), fontsize=8, color=INK_2, ha="left", va="bottom")
    ax1.xaxis.set_major_locator(mdates.MonthLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    ax1.xaxis.set_minor_locator(mdates.WeekdayLocator(byweekday=0))
    leg = [Line2D([], [], color=GAUGE_GRAY, lw=2.6, label="USGS gauge"), Line2D([], [], color=RETUNED_COLOR, lw=1.4, label="model, re-tuned (As 300, f x0.15)"),
           Line2D([], [], color=ASDRAWN_COLOR, lw=1.4, label="model, as drawn, loss off"), Patch(color="#efeeea", label="storm window (last rain day + 3 days)")]
    ax1.legend(handles=leg, loc="upper left", frameon=False, fontsize=9, ncol=2, bbox_to_anchor=(0.0, 1.0))
    for d_, t_, y_ in (("2016-11-03", "Nov 3", 0.35), ("2016-11-20", "Nov 20-21", 0.1)):
        ax1.annotate(t_ + ":\nmodel flows,\ngauge zero", (pd.Timestamp(d_) + pd.Timedelta(days=4), y_), fontsize=8, color=INK_2, ha="left", va="bottom")
    ax1.annotate("Jan 23-24:\ngauge flows more\nthan the model", (pd.Timestamp("2017-01-26"), 1.0), fontsize=8, color=INK_2, ha="left", va="bottom")
    fig.suptitle("WY2017, Nov-Mar: daily outlet volume, gauge vs model. The model's recessions drop away faster than the gauge's", x=0.01, ha="left", fontsize=12)
    fig.subplots_adjust(left=0.075, right=0.99, top=0.93, bottom=0.07)
    fig.savefig(path, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
