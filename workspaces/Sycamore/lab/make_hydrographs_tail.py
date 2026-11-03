"""Hydrographs of the most telling W4 tail-screen runs vs the USGS 09510200 gauge. Reads the outlet files; writes 5 PNGs.

Run it from workspaces/Sycamore/lab/ (it only reads files and writes pictures; no model is run, nothing is changed):
    python make_hydrographs_tail.py
    python make_hydrographs_tail.py --out-dir my_figs           # put the pictures somewhere else
Needs: results/W4_d012_k10_*_Outlet.qout for the 7 runs listed in RUNS below, and the gauge file
../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv. Use --results-dir / --obs-file if yours are elsewhere.
Writes (default folder tail_figs/): hydro_fig1_overview_log.png, hydro_fig2_tail_panels.png, hydro_fig3_tail_overlay.png,
hydro_fig4_flood_and_second_event.png, hydro_fig5_recession_log.png. Also prints the second-event peak times and daily means.
"""
import argparse, os
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

_ap = argparse.ArgumentParser(description="Hydrographs of the W4 tail-screen runs against the gauge.")
_ap.add_argument("--results-dir", default="results", help="folder holding the *_Outlet.qout files (default: results)")
_ap.add_argument("--obs-file", default="../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv", help="gauge flow file")
_ap.add_argument("--out-dir", default="tail_figs", help="where the PNG files go (created if missing)")
_args = _ap.parse_args()
RES, OBS, OUT = _args.results_dir, _args.obs_file, _args.out_dir
os.makedirs(OUT, exist_ok=True)
START = pd.Timestamp("2016-10-01 00:00")
PRE = "W4_d012_k10_"
# entity -> (file stem, label, color); colors = reference categorical slots 1-5, fixed per run in every figure
RUNS = {
    "base":   (PRE + "base_nocl",            "baseline (f x0.5874, As 1)",   "#2a78d6"),
    "f002":   (PRE + "f0p02_nocl",           "f x0.02, As 1",                "#6e6d69"),
    "hi300":  (PRE + "f0p5874_as300_nocl",   "f x0.5874, As 300",            "#52514e"),
    "f15a300":(PRE + "f0p15_as300_nocl",     "f x0.15, As 300",              "#eb6834"),
    "f02a300":(PRE + "f0p02_as300_nocl",     "f x0.02, As 300",              "#1baf7a"),
    "f02a500":(PRE + "f0p02_as500_nocl",     "f x0.02, As 500",              "#eda100"),
    "f02a900":(PRE + "f0p02_as900_nocl",     "f x0.02, As 900",              "#e87ba4"),
}
INK, INK2, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1", "#fcfcfb"
plt.rcParams.update({"figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF, "text.color": INK,
                     "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "axes.edgecolor": GRID,
                     "font.size": 10, "axes.spines.top": False, "axes.spines.right": False})

def read_q(stem):
    d = pd.read_csv(f"{RES}/{stem}_Outlet.qout")
    s = pd.Series(d["Qstrm_m3_s"].to_numpy(float), index=START + pd.to_timedelta(d["Time_hr"].to_numpy(), unit="h"))
    return s[~s.index.duplicated(keep="last")].sort_index().resample("60min").mean()

g = pd.read_csv(OBS, parse_dates=["datetime"]).set_index("datetime")["q_cms"]
g = pd.to_numeric(g, errors="coerce"); g = g[(g.index >= START) & (g.index < "2017-04-01")]
Q = {k: read_q(v[0]) for k, v in RUNS.items()}
PEAK = pd.Timestamp("2017-01-21 01:35")
BAND = (0.3, 1.0)

def style(ax):
    ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)

# ---------- Figure 1: whole record, log scale
fig, ax = plt.subplots(figsize=(11, 4.6))
for k in ("base", "f15a300", "f02a300"):
    ax.plot(Q[k].index, Q[k].clip(lower=1e-3), color=RUNS[k][2], lw=1.2, label=RUNS[k][1])
ax.plot(g.index, g.clip(lower=1e-3), color=INK, lw=1.4, label="gauge (USGS 09510200)")
ax.axhspan(*BAND, color="#2a78d6", alpha=0.06, lw=0)
ax.text(pd.Timestamp("2016-10-04"), 0.55, "0.3-1 m³/s band", color=INK2, fontsize=8.5, va="center")
ax.set_yscale("log"); ax.set_ylim(1e-2, 400)
ax.set_ylabel("Outlet flow, m³/s (log scale)")
ax.xaxis.set_major_locator(mdates.MonthLocator()); ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
ax.set_xlim(START, pd.Timestamp("2017-04-01")); style(ax)
ax.axvspan(START, pd.Timestamp("2017-01-01"), color="#e87ba4", alpha=0.07, lw=0)
ax.text(pd.Timestamp("2016-11-15"), 250, "Oct-Dec: gauge essentially dry, model flows (false flow)", ha="center", fontsize=8.5, color=INK2)
ax.set_title("Oct 2016 - Mar 2017: gauge vs three runs (flows under 0.001 m³/s are not drawn)", loc="left", fontsize=11, color=INK)
ax.legend(loc="upper right", frameon=False, fontsize=8.5, ncol=2)
fig.tight_layout(); fig.savefig(f"{OUT}/hydro_fig1_overview_log.png", dpi=170); plt.close(fig)

# ---------- Figure 2: tail close-up, one panel per run (small multiples, shared axes)
order = ["base", "f002", "hi300", "f15a300", "f02a300", "f02a500", "f02a900"]
x0, x1 = pd.Timestamp("2017-01-24"), pd.Timestamp("2017-02-20")
fig, axes = plt.subplots(2, 4, figsize=(13, 6.2), sharex=True, sharey=True)
for ax, k in zip(axes.flat, order):
    ax.axhspan(*BAND, color="#2a78d6", alpha=0.08, lw=0)
    ax.plot(g.index, g, color=INK, lw=1.3, label="gauge")
    ax.plot(Q[k].index, Q[k], color=RUNS[k][2] if k not in ("f002","hi300") else RUNS[k][2], lw=1.8, label=RUNS[k][1])
    for d in (7, 14, 21):
        ax.axvline(PEAK.normalize() + pd.Timedelta(days=d) + pd.Timedelta(hours=3), color=GRID, lw=1, ls=(0, (3, 3)), zorder=0)
    ax.set_xlim(x0, x1); ax.set_ylim(0, 2.4); style(ax)
    ax.set_title(RUNS[k][1], loc="left", fontsize=10, color=INK)
    ax.xaxis.set_major_locator(mdates.DayLocator(bymonthday=[25, 1, 8, 15])); ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax.tick_params(axis="x", labelsize=8)
ax = axes.flat[7]; ax.axis("off")
ax.plot([], [], color=INK, lw=1.3, label="gauge (black)"); ax.plot([], [], color=INK2, lw=1.8, label="model run (colored)")
ax.fill_between([], [], [], color="#2a78d6", alpha=0.15, label="0.3-1 m³/s band")
ax.plot([], [], color=MUTED, lw=1, ls=(0, (3, 3)), label="7, 14, 21 days after\nthe gauge flood peak")
ax.legend(loc="center left", frameon=False, fontsize=9)
for ax in axes[:, 0]: ax.set_ylabel("m³/s")
fig.suptitle("The tail, Jan 24 - Feb 20 2017: gauge (black) against each run. Flow above 2.4 m³/s is cut off.", x=0.01, ha="left", fontsize=11.5, color=INK)
fig.tight_layout(rect=(0, 0, 1, 0.96)); fig.savefig(f"{OUT}/hydro_fig2_tail_panels.png", dpi=170); plt.close(fig)

# ---------- Figure 3: tail overlay of the candidates on one axis
fig, ax = plt.subplots(figsize=(11, 4.8))
ax.axhspan(*BAND, color="#2a78d6", alpha=0.08, lw=0)
for k in ("base", "f15a300", "f02a300", "f02a500"):
    ax.plot(Q[k].index, Q[k], color=RUNS[k][2], lw=1.8, label=RUNS[k][1])
ax.plot(g.index, g, color=INK, lw=2.0, label="gauge")
ax.set_xlim(x0, x1); ax.set_ylim(0, 2.4); style(ax)
ax.xaxis.set_major_locator(mdates.DayLocator(interval=3)); ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
ax.set_ylabel("Outlet flow, m³/s")
ax.set_title("Tail: baseline and the three most useful runs against the gauge (shaded: 0.3-1 m³/s band)", loc="left", fontsize=11, color=INK)
ax.legend(frameon=False, fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=5)
fig.tight_layout(); fig.savefig(f"{OUT}/hydro_fig3_tail_overlay.png", dpi=170); plt.close(fig)

# ---------- Figure 4: the flood peak (A) and the second event on Jan 24 (B)
fig, (axa, axb) = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw={"width_ratios": [1, 1]})
for ax, (w0, w1, ymax, ttl) in zip((axa, axb), [("2017-01-20 12:00", "2017-01-23 00:00", 230, "A. The flood peak, Jan 20-22"),
                                                  ("2017-01-23 00:00", "2017-01-26 12:00", 22, "B. The second event, Jan 23-26 (note the smaller scale)")]):
    for k in ("base", "f15a300", "f02a300", "f02a900"):
        ax.plot(Q[k].index, Q[k], color=RUNS[k][2], lw=1.5, label=RUNS[k][1])
    ax.plot(g.index, g, color=INK, lw=1.8, label="gauge")
    ax.set_xlim(pd.Timestamp(w0), pd.Timestamp(w1)); ax.set_ylim(0, ymax); style(ax)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d %Hh")); ax.set_title(ttl, loc="left", fontsize=10.5, color=INK)
    ax.set_ylabel("Outlet flow, m³/s")
axa.xaxis.set_major_locator(mdates.HourLocator(byhour=[0, 12])); axb.xaxis.set_major_locator(mdates.HourLocator(byhour=[0, 12]))
axb.legend(frameon=False, fontsize=8.5, loc="upper right")
fig.tight_layout(); fig.savefig(f"{OUT}/hydro_fig4_flood_and_second_event.png", dpi=170); plt.close(fig)

# ---------- Figure 5: the whole recession on a log axis
fig, ax = plt.subplots(figsize=(11, 4.8))
for k in ("base", "f15a300", "f02a300", "f02a500"):
    ax.plot(Q[k].index, Q[k].clip(lower=1e-3), color=RUNS[k][2], lw=1.6, label=RUNS[k][1])
ax.plot(g.index, g.clip(lower=1e-3), color=INK, lw=1.9, label="gauge")
ax.axhspan(*BAND, color="#2a78d6", alpha=0.08, lw=0)
ax.set_yscale("log"); ax.set_ylim(0.03, 300); ax.set_xlim(pd.Timestamp("2017-01-20"), x1); style(ax)
ax.xaxis.set_major_locator(mdates.DayLocator(interval=3)); ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
ax.set_ylabel("Outlet flow, m³/s (log scale)")
ax.set_title("Recession after the January flood on a log axis: how far the models fall below the gauge", loc="left", fontsize=11, color=INK)
ax.legend(frameon=False, fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=5)
fig.tight_layout(); fig.savefig(f"{OUT}/hydro_fig5_recession_log.png", dpi=170); plt.close(fig)

# numbers for the write-up
print("gauge peak m3/s", g.loc["2017-01-20":"2017-01-23"].max(), g.loc["2017-01-20":"2017-01-23"].idxmax())
for k in order:
    q = Q[k]
    print(f"{k:8s} peak {q.loc['2017-01-20':'2017-01-23'].max():7.1f}  at {q.loc['2017-01-20':'2017-01-23'].idxmax()}  d3-5 mean {q.loc['2017-01-24':'2017-01-26'].mean():6.2f}  Jan26-Feb2 mean {q.loc['2017-01-26':'2017-02-02'].mean():5.2f}  Feb5-12 mean {q.loc['2017-02-05':'2017-02-12'].mean():5.2f}")
print("gauge    ", f"d3-5 mean {g.loc['2017-01-24':'2017-01-26'].mean():6.2f}  Jan26-Feb2 mean {g.loc['2017-01-26':'2017-02-02'].mean():5.2f}  Feb5-12 mean {g.loc['2017-02-05':'2017-02-12'].mean():5.2f}")

# second-event timing and size
for name, q in [("gauge", g)] + [(k, Q[k]) for k in order]:
    w = q.loc["2017-01-23 12:00":"2017-01-26 00:00"]
    print(f"{name:8s} second-event peak {w.max():6.2f} m3/s at {w.idxmax()}")
for name, q in [("gauge", g)] + [(k, Q[k]) for k in ("base","f15a300","f02a300","f02a900")]:
    w = q.loc["2017-01-20 18:00":"2017-01-21 06:00"]
    first = w[w > 5].index[0] if (w > 5).any() else None
    print(f"{name:8s} first hour above 5 m3/s: {first}")
for name, q in [("gauge", g)] + [(k, Q[k]) for k in ("base","f15a300","f02a300","f02a500")]:
    for d in ("2017-01-28","2017-02-04","2017-02-11","2017-02-17"):
        pass
    print(f"{name:8s} mean flow Jan28 {q.loc['2017-01-28'].mean():6.2f}  Feb4 {q.loc['2017-02-04'].mean():5.2f}  Feb11 {q.loc['2017-02-11'].mean():5.2f}  Feb17 {q.loc['2017-02-17'].mean():5.2f}")
