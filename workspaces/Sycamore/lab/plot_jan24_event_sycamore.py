"""Jan 23-24, 2017 event check for the W4 tail-screen runs: rain, runoff the model MAKES, and flow that REACHES the outlet.

WHAT IT MAKES (2 pictures + a printed table), all read-only:
  fig_jan24_rain_vs_flow.png         gauge rain and radar rain (hourly bars) over observed and modeled outlet flow, Jan 22-25
  fig_jan24_generated_vs_outlet.png  rain the model used (basin mean) / runoff generated on the land (all cells added, before
                                     routing) / flow at the outlet versus the gauge
  printed table                      generated volume versus outlet volume (M m3) for the flood (Jan 20 12:00 - Jan 23 06:00)
                                     and for the Jan 23-24 event (Jan 23 06:00 - Jan 25 12:00), plus the Jan 18-25 week ratio

RUN IT from workspaces/Sycamore/lab/ (nothing is run or changed; it only reads files and writes pictures):
    python plot_jan24_event_sycamore.py --pairs radar_gauge_hourly_pairs.csv
    python plot_jan24_event_sycamore.py --pairs path/to/radar_gauge_hourly_pairs.csv --out-dir my_figs

FILES IT READS (all must exist; the run names are the W4 tail-screen runs, change --runs to use others):
  results/W4_d012_k10_<run>_nocl4368_00.mrf   basin-average hourly series (Srf = runoff generated; MAP = basin mean rain)
  results/W4_d012_k10_<run>_nocl_Outlet.qout  outlet flow
  ../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv   the gauge
  radar_gauge_hourly_pairs.csv                from the radar-versus-gauge check (hourly gauge and radar rain; --pairs)
Run names in --runs are the part after W4_d012_k10_ and before _nocl, for example base, f0p15_as300, f0p02_as300.

WORDS
  generated runoff   the sum over all cells of the runoff made that hour (the Srf column of the .mrf file; it equals the four
                     runoff types in the .rft file added together). It is NOT the flow at the outlet.
  hour labels        model files label an hour by its END (Time_hr 2736 = Jan 23 00:00 means the hour before it). The rain
                     pairs file labels hours by their START. The pictures use the labels as the files give them; the gap is
                     one hour and does not change what you see at this scale.
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

PRE, SUF = "W4_d012_k10_", "_nocl"
T0 = pd.Timestamp("2016-10-01 00:00")
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
LABELS = {"base": "Base", "f0p15_as300": "f x0.15, As 300", "f0p02_as300": "f x0.02, As 300",
          "f0p02_as500": "f x0.02, As 500", "f0p02_as900": "f x0.02, As 900"}

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--pairs", default="radar_gauge_hourly_pairs.csv", help="the hourly gauge/radar rain file")
ap.add_argument("--results-dir", default="results")
ap.add_argument("--obs-file", default="../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv")
ap.add_argument("--runs", default="base,f0p15_as300,f0p02_as300", help="comma list; up to 5 runs")
ap.add_argument("--out-dir", default="tail_figs")
args = ap.parse_args()
runs = [r.strip() for r in args.runs.split(",") if r.strip()]
if len(runs) > 5:
    sys.exit("Use at most 5 runs (the figure has 5 colors).")
os.makedirs(args.out_dir, exist_ok=True)


def need(path):
    if not os.path.isfile(path):
        sys.exit(f"Cannot find {path}\nRun this from workspaces/Sycamore/lab/ or use the --results-dir / --obs-file / --pairs options.")
    return path


def read_outlet(run):
    d = pd.read_csv(need(f"{args.results_dir}/{PRE}{run}{SUF}_Outlet.qout"))
    return pd.Series(d["Qstrm_m3_s"].to_numpy(float), index=T0 + pd.to_timedelta(d["Time_hr"].to_numpy(float), unit="h"))


def read_mrf(run):
    d = pd.read_csv(need(f"{args.results_dir}/{PRE}{run}{SUF}4368_00.mrf"), usecols=["Time_hr", "Srf_m3_s", "MAP_mm_hr"])
    d.index = T0 + pd.to_timedelta(d["Time_hr"].to_numpy(float), unit="h")
    return d


obs = pd.read_csv(need(args.obs_file), parse_dates=["datetime"]).set_index("datetime")["q_cms"]
pairs = pd.read_csv(need(args.pairs), parse_dates=["hour_start"])
pairs["g"] = pairs["gauge"].str.replace(r"\s*\(.*", "", regex=True)
good = pairs[pairs["g"] != "Mt. Ord"]   # Mt. Ord undercounts snow; left out (same choice as the radar check)
gauge_rain = good.groupby("hour_start")["gauge_mm"].mean()
radar_rain = good.groupby("hour_start")["radar_3x3_mean_mm"].mean()

outlet = {r: read_outlet(r) for r in runs}
mrf = {r: read_mrf(r) for r in runs}
name = lambda r: LABELS.get(r, r)

# ---------- figure 1: rain over flow ----------
a, b = pd.Timestamp("2017-01-22 12:00"), pd.Timestamp("2017-01-25 12:00")
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 7), sharex=True, gridspec_kw={"height_ratios": [1, 2.2], "hspace": 0.06})
w = pd.Timedelta(minutes=24)
rg, rr = gauge_rain[a:b], radar_rain[a:b]
ax1.bar(rg.index + pd.Timedelta(minutes=18), rg.values, width=w, color="#444444", label="Gauges (Sycamore Upper + Iron Dike mean)")
ax1.bar(rr.index + pd.Timedelta(minutes=42), rr.values, width=w, color="#7fb3e8", label="Radar (3x3 mean at those gauges)")
ax1.invert_yaxis(); ax1.set_ylabel("Rain (mm/h)"); ax1.legend(loc="lower right", fontsize=8, frameon=False)
ax2.plot(obs[a:b].index, obs[a:b].values, color="black", lw=2.2, label="Gauge 09510200 (observed)")
for r, c in zip(runs, COLORS):
    ax2.plot(outlet[r][a:b].index, outlet[r][a:b].values, color=c, lw=1.6, label=name(r))
ax2.set_ylabel("Outlet flow (m$^3$/s)"); ax2.legend(loc="upper left", fontsize=8, frameon=False)
ax2.set_xlabel("Phoenix time (MST); rain bars sit at the START of each hour")
for ax in (ax1, ax2):
    ax.grid(alpha=0.25)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
ax1.set_title("Jan 23-24, 2017: rain input versus observed and modeled outlet flow", loc="left", fontsize=11)
f1 = os.path.join(args.out_dir, "fig_jan24_rain_vs_flow.png")
fig.savefig(f1, dpi=150, bbox_inches="tight"); plt.close(fig)

# ---------- figure 2: generated versus outlet ----------
a, b = pd.Timestamp("2017-01-22 18:00"), pd.Timestamp("2017-01-25 12:00")
fig, ax = plt.subplots(3, 1, figsize=(11, 9), sharex=True, gridspec_kw={"height_ratios": [1, 1.6, 2], "hspace": 0.07})
x = mrf[runs[0]]["MAP_mm_hr"][a:b]
ax[0].bar(x.index - pd.Timedelta(minutes=30), x.values, width=pd.Timedelta(minutes=50), color="#7fb3e8")
ax[0].invert_yaxis(); ax[0].set_ylabel("Rain the model got\n(basin mean, mm/h)")
for r, c in zip(runs, COLORS):
    y = mrf[r]["Srf_m3_s"][a:b]
    ax[1].plot(y.index, y.values, color=c, lw=1.8, label=name(r))
    ax[2].plot(outlet[r][a:b].index, outlet[r][a:b].values, color=c, lw=1.6, label=name(r))
ax[1].set_ylabel("Runoff generated on the land\n(all cells added, m$^3$/s)"); ax[1].legend(frameon=False, fontsize=8, loc="upper right")
ax[2].plot(obs[a:b].index, obs[a:b].values, color="black", lw=2.2, label="Gauge 09510200 (observed)")
ax[2].set_ylabel("Flow at the outlet (m$^3$/s)"); ax[2].legend(frameon=False, fontsize=8, loc="upper left")
ax[2].set_xlabel("Phoenix time (MST)")
for k in ax:
    k.grid(alpha=0.25)
    for s in ("top", "right"):
        k.spines[s].set_visible(False)
ax[0].set_title("Jan 23-24, 2017: runoff the model makes versus runoff that reaches the outlet", loc="left", fontsize=11)
f2 = os.path.join(args.out_dir, "fig_jan24_generated_vs_outlet.png")
fig.savefig(f2, dpi=150, bbox_inches="tight"); plt.close(fig)

# ---------- the numbers ----------
def vol(series, lo, hi):
    """Volume in million m3 from a series of m3/s values (hourly means), between lo and hi."""
    return series[lo:hi].resample("1h").mean().sum() * 3600 / 1e6


print()
print("Volumes in million m3 (hourly sums; 'generated' = runoff made on the land, 'outlet' = flow that reached the outlet)")
spans = [("Flood   Jan 20 12:00 - Jan 23 06:00", "2017-01-20 12:00", "2017-01-23 06:00"),
         ("Jan 23-24 event Jan 23 06:00 - Jan 25 12:00", "2017-01-23 06:00", "2017-01-25 12:00")]
print(f"  gauge: " + "   ".join(f"{lab.split('  ')[0].strip()} {vol(obs, lo, hi):.2f}" for lab, lo, hi in spans)
      + f"   week Jan 18-25 {vol(obs, '2017-01-18', '2017-01-25'):.2f}")
for r in runs:
    g = mrf[r]["Srf_m3_s"]
    parts = []
    for lab, lo, hi in spans:
        ge, ou = vol(g, lo, hi), vol(outlet[r], lo, hi)
        parts.append(f"{lab.split('  ')[0].strip()}: generated {ge:.2f}, outlet {ou:.2f}")
    wk = vol(outlet[r], "2017-01-18", "2017-01-25") / vol(obs, "2017-01-18", "2017-01-25")
    print(f"  {name(r):<16} " + " | ".join(parts) + f" | week model/gauge {wk:.2f}")
print("  (the outlet volumes include flow that was already there before the rain; for the base run the Jan 23-24 outlet volume")
print("   is about what its own recession gives, so it delivers essentially none of what it generates)")
rain_model = mrf[runs[0]]["MAP_mm_hr"]["2017-01-23 06:00":"2017-01-25 12:00"].sum()
print(f"Rain the model used, Jan 23 06:00 - Jan 25 12:00: {rain_model:.1f} mm (basin mean)")
print(f"Rain at the two good gauges, Jan 23 09:00 - Jan 24 04:00 (start-of-hour labels): gauges {gauge_rain['2017-01-23 09:00':'2017-01-24 03:00'].sum():.1f} mm, radar {radar_rain['2017-01-23 09:00':'2017-01-24 03:00'].sum():.1f} mm")
print(f"\nWrote {f1}\n      {f2}")
