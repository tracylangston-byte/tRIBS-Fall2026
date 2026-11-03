"""
radar_burst_check_sycamore.py  --  Read-only look at the radar versus the rain gauges for the Feb 27-28 and Feb 18-19 storms, hour by hour.
It reads two files you already have and changes nothing:
    radar_gauge_hourly_pairs.csv   (written by radar_gauge_check_sycamore.py; hours labeled by the hour START and END, Arizona time)
    sycamore_09510200_obs_flow_WY2017.csv   (the USGS gauge, 15-minute flow, optional)
QUESTIONS (from the mini re-write, test 2):
    1. Feb 27-28: was the radar rain for the first (evening) burst too heavy compared with the gauges? If it was, part of the early shoulder
       at the model's outlet is a forcing problem. If it was not, the shoulder comes from how the model turns that rain into flow.
    2. Feb 18-19: does the radar line up with the gauges, and is the 7.6 h centre-of-rain offset at Sycamore Upper (Stage 0) a clock/timing
       offset or a missing burst?
"Good gauges" = Sycamore Cr Upper and Iron Dike (Mt. Ord is shown but left out of the averages, as in Stage 0: its bucket misses snow).
Run (from lab/):   python radar_burst_check_sycamore.py --pairs <path to radar_gauge_hourly_pairs.csv> --flow <path to the USGS flow csv>
Everything is in hours ENDING (the radar file convention), Arizona time. The burst windows below can be changed with --burst1/--burst2.
"""
import argparse
import numpy as np
import pandas as pd

SHORT = {"Mt. Ord (2000)": "MtOrd", "Sycamore Cr Upper (82200)": "SycUp", "Iron Dike (83000)": "IronD"}
GOOD = ["SycUp", "IronD"]


def load(pairs):
    p = pd.read_csv(pairs, parse_dates=["hour_start", "hour_ending"])
    p["g"] = p["gauge"].map(SHORT)
    return p


def wide(p, a, b):
    d = p[(p.hour_ending > pd.Timestamp(a)) & (p.hour_ending <= pd.Timestamp(b))]
    out = pd.DataFrame(index=sorted(d.hour_ending.unique()))
    for g, s in d.groupby("g"):
        s = s.set_index("hour_ending")
        out[f"{g} gauge"], out[f"{g} radar3x3"] = s.gauge_mm, s.radar_3x3_mean_mm
    return out


def totals(p, a, b, label):
    d = p[(p.hour_ending > pd.Timestamp(a)) & (p.hour_ending <= pd.Timestamp(b))]
    parts = [f"{g} gauge {s.gauge_mm.sum():.1f} / radar {s.radar_3x3_mean_mm.sum():.1f} mm" for g, s in d.groupby("g")]
    gd = d[d.g.isin(GOOD)].groupby("hour_ending")[["gauge_mm", "radar_3x3_mean_mm"]].mean().sum()
    print(f"  {label:<34} hours ending {pd.Timestamp(a):%b %d %H:%M} to {pd.Timestamp(b):%b %d %H:%M}\n      " + "; ".join(parts)
          + f"\n      good-gauge mean: gauge {gd.gauge_mm:.1f} mm, radar {gd.radar_3x3_mean_mm:.1f} mm, radar/gauge {gd.radar_3x3_mean_mm / gd.gauge_mm:.2f}")


def centre(p, g, a, b, col):
    d = p[(p.g == g) & (p.hour_ending > pd.Timestamp(a)) & (p.hour_ending <= pd.Timestamp(b))]
    w = d[col].clip(lower=0).to_numpy()
    t = ((d.hour_ending - pd.Timestamp(a)) / pd.Timedelta(hours=1)).to_numpy() - 0.5
    return float((t * w).sum() / w.sum()) if w.sum() > 0 else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--flow", default=None)
    ap.add_argument("--burst1", nargs=2, default=["2017-02-27 14:00", "2017-02-28 02:00"], help="hours-ending window of the first Feb 27 burst")
    ap.add_argument("--burst2", nargs=2, default=["2017-02-28 03:00", "2017-02-28 11:00"], help="hours-ending window of the second (morning) burst")
    a = ap.parse_args()
    p = load(a.pairs)
    q = None
    if a.flow:
        q = pd.read_csv(a.flow, parse_dates=["datetime"]).set_index("datetime")["q_cms"]
    pd.set_option("display.width", 220)
    print("1. FEB 27-28, hour by hour (mm per hour, hours ending; radar = mean of the 3x3 cells)")
    w = wide(p, "2017-02-27 12:00", "2017-02-28 14:00")
    if q is not None:
        w["USGS Q m3/s"] = q.resample("1h", closed="right", label="right").max().reindex(w.index)
    print(w.round(1).to_string())
    print("\n  BURST TOTALS")
    totals(p, *a.burst1, "first burst (Feb 27 evening)")
    totals(p, *a.burst2, "second burst (Feb 28 morning)")
    if q is not None:
        seg = q[a.burst1[0]:a.burst2[1]]
        print(f"  USGS outlet flow from {a.burst1[0]} to the end of burst 2: min {seg.min():.2f}, max before 04:00 Feb 28 {q[a.burst1[0]:'2017-02-28 04:00'].max():.2f} m3/s")
    print("\n2. FEB 18-19, hour by hour (hours with at least 0.3 mm at a good gauge)")
    w = wide(p, "2017-02-18 00:00", "2017-02-20 12:00")
    if q is not None:
        w["USGS Q m3/s"] = q.resample("1h", closed="right", label="right").max().reindex(w.index)
    keep = (w[[f"{g} {k}" for g in GOOD for k in ("gauge", "radar3x3")]].fillna(0) > 0.3).any(axis=1)
    print(w[keep].round(1).to_string())
    print("\n  TOTALS")
    totals(p, "2017-02-18 00:00", "2017-02-19 12:00", "Feb 18 to midday Feb 19")
    totals(p, "2017-02-19 12:00", "2017-02-20 12:00", "midday Feb 19 to midday Feb 20")
    totals(p, "2017-02-18 00:00", "2017-02-20 12:00", "whole Feb 18-19 storm")
    print("\n  CENTRE OF RAIN at Sycamore Upper, hours after Feb 18 00:00 (the Stage 0 offset was radar 7.6 h EARLIER than the gauge)")
    for lab, (x, y) in {"whole storm to Feb 20 12:00": ("2017-02-18 00:00", "2017-02-20 12:00"), "without Feb 19 after 12:00": ("2017-02-18 00:00", "2017-02-19 12:00")}.items():
        cg, cr = centre(p, "SycUp", x, y, "gauge_mm"), centre(p, "SycUp", x, y, "radar_3x3_mean_mm")
        print(f"      {lab:<30} gauge {cg:5.1f} h, radar {cr:5.1f} h, radar minus gauge {cr - cg:+.1f} h")


if __name__ == "__main__":
    main()
