#!/usr/bin/env python3
"""
plot_overlays_stage2_sycamore.py -- Stage 2 hydrograph check: overlays for the biggest-moving knobs.

For each knob it makes one figure (overlay_<knob>.png) with every screened run of that knob drawn against the gauge.

Soil knobs (ks, f, psib, m, ts, an, br) come from the W2 soil screen (Jan 18 - Mar 12):
  1. flood week Jan 18-25, linear flow axis        2. flood week, log flow axis
  3. zoom Feb 18-19 (Feb 18-21 so the recession shows)   4. zoom Feb 27 - Mar 2
  5. cumulative volume over the whole W2 window, model and gauge
Routing knobs (cv, flowexp, n, wc) come from the W1 routing screen (Jan 18-25, one week):
  1. flood week, linear     2. flood week, log flow axis
  3. zoom on the flood peak (Jan 20 18:00 - Jan 21 12:00), where timing and sharpness show
  4. Jan 23-25 (the second January storm and the gauge's slow tail)
  5. cumulative volume over the week, model and gauge
The baseline run is drawn thick; the others are colored by knob value (cool to warm = low to high); the gauge is black.

Run it from workspaces/Sycamore/lab/ (it reads results/<run>_Outlet.qout and the gauge file; nothing is run or changed):
    python plot_overlays_stage2_sycamore.py                     # ks, f, psib, cv, flowexp (the five biggest movers)
    python plot_overlays_stage2_sycamore.py --knobs cv,flowexp   # just some knobs
Needs score_run_sycamore.py in the same folder (it uses that script's own file readers), screen_soil_summary.csv and
screen_routing_summary.csv.
Writes: stage2_overlays/overlay_<knob>.png (the folder is created if it is not there).
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import score_run_sycamore as sc

ROUTING_KNOBS = {"cv", "rfe", "n", "wc"}          # the W1 routing screen; every other knob is a W2 soil knob
ALIASES = {"flowexp": "rfe"}                      # the routing summary calls flowexp "rfe"
DISPLAY = {"rfe": "flowexp"}                      # what the title and the file name call it
BASELINES = {"W1": "W1_ks0p2_cv6_psib0p5_nocl", "W2": "W2_ks0p2_cv6_psib0p5_nocl"}
PANELS = {
    "W2": [("Flood week Jan 18-25 (linear)", "2017-01-18", "2017-01-25", False),
           ("Flood week Jan 18-25 (log flow)", "2017-01-18", "2017-01-25", True),
           ("Feb 18-21 (zoom)", "2017-02-18", "2017-02-21", False),
           ("Feb 27 - Mar 2 (zoom)", "2017-02-27", "2017-03-02", False)],
    "W1": [("Flood week Jan 18-25 (linear)", "2017-01-18", "2017-01-25", False),
           ("Flood week Jan 18-25 (log flow)", "2017-01-18", "2017-01-25", True),
           ("Flood peak zoom (Jan 20 18:00 - Jan 21 12:00)", "2017-01-20 18:00", "2017-01-21 12:00", False),
           ("Jan 23-25 (second January storm and tail)", "2017-01-23", "2017-01-25", False)],
}
LOG_FLOOR = 1e-3


def cum_volume(q: pd.Series, a, b) -> pd.Series:
    """Cumulative volume (million m3) from a to b, trapezoid rule on the series' own time steps."""
    s = q[(q.index >= a) & (q.index <= b)].dropna()
    t = (s.index - s.index[0]).total_seconds().to_numpy()
    c = sc.cumulative(t, s.to_numpy(float)) / 1e6
    return pd.Series(c, index=s.index)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--knobs", default="ks,f,psib,cv,flowexp")
    ap.add_argument("--summary", default="screen_soil_summary.csv", help="W2 soil screen summary")
    ap.add_argument("--routing-summary", default="screen_routing_summary.csv", help="W1 routing screen summary")
    ap.add_argument("--results", default="results")
    ap.add_argument("--obs", default="../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv")
    a = ap.parse_args()
    gauge_all = sc.read_obs(a.obs)
    for knob in [ALIASES.get(k.strip(), k.strip()) for k in a.knobs.split(",") if k.strip()]:
        win = "W1" if knob in ROUTING_KNOBS else "W2"
        panels, baseline = PANELS[win], BASELINES[win]
        start, end = pd.Timestamp(sc.RUN_WINDOWS[win][0]), pd.Timestamp(sc.RUN_WINDOWS[win][1])
        gauge = gauge_all[(gauge_all.index >= start) & (gauge_all.index <= end)]
        summ_path = a.routing_summary if win == "W1" else a.summary
        if not Path(summ_path).is_file():
            print(f"missing {summ_path}, so {knob} was skipped")
            continue
        summ = pd.read_csv(summ_path)
        name = DISPLAY.get(knob, knob)
        rows = summ[summ["knob"] == knob].sort_values("value")
        if rows.empty:
            print(f"no rows for knob {knob!r} in {summ_path}")
            continue
        cmap = plt.get_cmap("coolwarm")
        vals = rows["value"].to_numpy(float)
        fig, axes = plt.subplots(2, 3, figsize=(17, 9))
        axf = axes.ravel()
        for pos, (_, r) in enumerate(rows.iterrows()):
            qf = Path(a.results) / f"{r['run_id']}_Outlet.qout"
            if not qf.is_file():
                print(f"missing {qf}")
                continue
            q = sc.read_qout(qf, start)
            is_base = r["run_id"] == baseline
            frac = 0.5 if len(vals) == 1 else np.argsort(np.argsort(vals))[pos] / (len(vals) - 1)
            col = "tab:green" if is_base else cmap(frac)
            lw = 2.6 if is_base else 1.3
            label = f"{name} = {r['value']:g}" + (" (baseline)" if is_base else "")
            for ax, (title, t0, t1, logy) in zip(axf[:4], panels):
                s = q[(q.index >= t0) & (q.index <= pd.Timestamp(t1))]
                y = np.maximum(s.to_numpy(float), LOG_FLOOR) if logy else s.to_numpy(float)
                ax.plot(s.index, y, color=col, lw=lw, label=label)
            c = cum_volume(q, start, end)
            axf[4].plot(c.index, c.to_numpy(), color=col, lw=lw, label=label)
        for ax, (title, t0, t1, logy) in zip(axf[:4], panels):
            g = gauge[(gauge.index >= t0) & (gauge.index <= pd.Timestamp(t1))]
            gy = np.maximum(g.to_numpy(float), LOG_FLOOR) if logy else g.to_numpy(float)
            ax.plot(g.index, gy, color="black", lw=1.8, label="gauge")
            ax.set_title(title)
            ax.set_ylabel("flow (m3/s)")
            if logy:
                ax.set_yscale("log")
                ax.set_ylim(LOG_FLOOR, None)
            ax.tick_params(axis="x", rotation=30)
            ax.grid(alpha=0.3)
        gc = cum_volume(gauge, start, end)
        axf[4].plot(gc.index, gc.to_numpy(), color="black", lw=1.8, label="gauge")
        axf[4].set_title(f"Cumulative volume, {start:%b %d} - {end:%b %d}")
        axf[4].set_ylabel("million m3")
        axf[4].tick_params(axis="x", rotation=30)
        axf[4].grid(alpha=0.3)
        axf[5].axis("off")
        h, l = axf[0].get_legend_handles_labels()
        axf[5].legend(h, l, loc="center", fontsize=10, title=f"{name}: green = baseline, black = gauge")
        fig.suptitle(f"Stage 2 overlay ({win}): {name} (all other settings at the baseline)", fontsize=14)
        fig.tight_layout()
        out_dir = Path("stage2_overlays")
        out_dir.mkdir(exist_ok=True)
        out = out_dir / f"overlay_{name}.png"
        fig.savefig(out, dpi=110)
        plt.close(fig)
        print("wrote", out)


if __name__ == "__main__":
    main()
