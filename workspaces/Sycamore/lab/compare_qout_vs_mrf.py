"""
compare_qout_vs_mrf.py   (version 1)   READ-ONLY: writes nothing.

Question it answers: in the Sycamore runs, does the ROUTED outlet flow (the .qout file,
column Qstrm) respond to channel loss, and how does it differ from the .mrf column Srf
that the compare cell has been plotting against the gauge?

Run it from lab/ (same folder the compare cell runs from):

    python compare_qout_vs_mrf.py Sycamore_jan_ks0p3 Sycamore_jan_ks0p3_nocl

Optional: --start '2017-01-18 00:00'  --gauge ../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv
          --qout_a PATH --qout_b PATH   (use if the script says it found more than one .qout file)
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

VERSION = "1"


def find_qout(name, others):
    """All .qout files whose name starts with this run name, minus files that belong to a longer run name."""
    hits = [p for p in glob.glob("**/*.qout*", recursive=True) if os.path.basename(p).startswith(name)]
    longer = [o for o in others if o != name and o.startswith(name)]
    return sorted(p for p in hits if not any(os.path.basename(p).startswith(o) for o in longer))


def read_qout(path):
    """tRIBS 6.0.0 .qout: comma-delimited with a header line. Older style: whitespace, one header line."""
    with open(path) as fh:
        fh.readline()
        second = fh.readline()
    if "," in second:
        raw = pd.read_csv(path, sep=",", header=0)
    else:
        raw = pd.read_csv(path, sep=r"\s+", skiprows=1, header=None)
    raw = raw.iloc[:, :2].copy()
    raw.columns = ["Time_hr", "Q"]
    return raw.apply(pd.to_numeric, errors="coerce").dropna()


def read_mrf(name):
    files = sorted(glob.glob("**/%s[0-9]*.mrf" % name, recursive=True))
    if not files:
        sys.exit("STOP: no .mrf file found for %r" % name)
    m = pd.read_csv(files[-1])
    m = m[m["RainMin_mm_hr"] < 9000]
    return m


def trapz_volume(t_hr, q):
    t = np.asarray(t_hr, float) * 3600.0
    q = np.asarray(q, float)
    ok = ~np.isnan(q[1:]) & ~np.isnan(q[:-1])
    return float(np.sum(np.where(ok, 0.5 * (q[1:] + q[:-1]) * np.diff(t), 0.0)))


def summarize(t_hr, q):
    q = np.asarray(q, float)
    i = int(np.nanargmax(q))
    return {"peak_m3s": q[i], "t_peak_h": float(np.asarray(t_hr)[i]), "volume_Mm3": trapz_volume(t_hr, q) / 1e6}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_a")
    ap.add_argument("run_b")
    ap.add_argument("--start", default="2017-01-18 00:00")
    ap.add_argument("--gauge", default="../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv")
    ap.add_argument("--qout_a")
    ap.add_argument("--qout_b")
    a = ap.parse_args()
    print("compare_qout_vs_mrf.py version %s (read-only)" % VERSION)

    paths = {}
    for run, given in ((a.run_a, a.qout_a), (a.run_b, a.qout_b)):
        found = [given] if given else find_qout(run, [a.run_a, a.run_b])
        if len(found) != 1:
            allq = sorted(glob.glob("**/*.qout*", recursive=True))
            sys.exit("STOP: expected exactly one .qout file for %r, found %d: %s\n"
                     "All .qout files under this folder: %s\n"
                     "Pass the right one with --qout_a / --qout_b." % (run, len(found), found, allq[:20]))
        paths[run] = found[0]
        print("qout for %-28s %s" % (run, paths[run]))

    rows = {}
    series = {}
    for run in (a.run_a, a.run_b):
        q = read_qout(paths[run])
        m = read_mrf(run)
        series[run] = (q, m)
        rows["%s | qout Qstrm (routed)" % run] = summarize(q["Time_hr"], q["Q"])
        rows["%s | mrf Srf" % run] = summarize(m["Time_hr"], m["Srf_m3_s"])

    t0 = pd.Timestamp(a.start)
    if os.path.exists(a.gauge):
        g = pd.read_csv(a.gauge, parse_dates=["datetime"]).set_index("datetime")["q_cms"]
        end =t0 + pd.Timedelta(hours=max(series[a.run_a][0]["Time_hr"].max(), 1))
        gw = g.loc[t0:end]
        gt = (gw.index - t0).total_seconds() / 3600.0
        rows["gauge"] = summarize(gt, gw.values)
    else:
        print("(gauge file not found at %s; skipping the gauge row)" % a.gauge)

    tab = pd.DataFrame(rows).T
    if "gauge" in rows:
        tab["peak vs gauge (h)"] = tab["t_peak_h"] - rows["gauge"]["t_peak_h"]
    print("\nPEAK AND VOLUME OVER THE WHOLE FILE (time = hours since the run start, %s)" % a.start)
    print(tab.round(3).to_string())

    qa, ma = series[a.run_a]
    qb, mb = series[a.run_b]
    qb_on_a = np.interp(qa["Time_hr"], qb["Time_hr"], qb["Q"])
    d = np.abs(qa["Q"].to_numpy() - qb_on_a)
    print("\nROUTED OUTLET FLOW, %s minus %s" % (a.run_a, a.run_b))
    print("  largest difference in Qstrm (m3/s): %.6g" % d.max())
    print("  volume difference (m3):              %.6g" % (trapz_volume(qa["Time_hr"], qa["Q"]) - trapz_volume(qb["Time_hr"], qb["Q"])))
    print("  ChannelPerc total (m3): %s = %.6g   %s = %.6g" % (a.run_a, ma["ChannelPerc_m3"].sum(), a.run_b, mb["ChannelPerc_m3"].sum()))
    for run in (a.run_a, a.run_b):
        q, m = series[run]
        lag = rows["%s | qout Qstrm (routed)" % run]["t_peak_h"] - rows["%s | mrf Srf" % run]["t_peak_h"]
        print("  %-28s Qstrm peak minus Srf peak: %+.2f h" % (run, lag))
    print("\nHow to read it: if the Qstrm difference is above zero and the channel-loss totals differ, the routed")
    print("flow feels the loss and Srf is the wrong series to compare with the gauge. If Qstrm is also identical,")
    print("the loss is not reaching the outlet flow at all.")


if __name__ == "__main__":
    main()
