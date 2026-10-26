"""
plot_series120_hydrograph_shape.py
==================================
Answers one question with pictures: when a run gets the TOTAL VOLUME of a storm right
(within 2% of the synthetic truth), does it also get the SHAPE of the hydrograph right?

It uses the three storms of Series 120 (tRIBS 6.0.0, South Mountain Fan, synthetic truth
Ks 7.0 / f 0.012):

    0.8x storm   (rain forcing multiplied by 0.8)
    1.0x storm   (rain as recorded)
    1.25x storm  (rain forcing multiplied by 1.25)

For each storm it picks the runs whose volume is within 2% of that storm's truth (the
runs sitting on the zero-bias "ridge" in the Ks-f plane; the circled dots in
fig_series120_storms_1_pbias_maps.png) and draws their hydrographs against the true one.

THE TWO FIGURES (both PNG)
    1  fig_series120_hydrograph_shape_1_hydrographs.png
       One column per storm, three rows:
         top     the main flood wave. Black = the synthetic truth. Each coloured line
                 is one volume-matched run. (The window covers the main flood only. If
                 a run carries a small extra rise hours away from it, the window does
                 not stretch to include it: the script prints it and notes it under the
                 figure instead. See "flow outside the window" below.)
         middle  a zoom on the top of the peak (about 35 minutes either side), where a
                 miss in peak height and peak timing is easiest to see.
         bottom  the miss itself: run minus truth, as a percent of the truth's peak flow,
                 at every time step. Zero = the run is exactly on the truth. Above zero =
                 the run is carrying more flow than the truth at that moment.
       Colour = the run's f:  blue = f below 0.007,  grey = f from 0.007 to 0.021 (the
       truth's f is 0.012),  red = f at or above 0.021. The cut-offs are in SETTINGS and can
       be changed with --f_low and --f_high.

    2  fig_series120_hydrograph_shape_2_shape_vs_f.png
       The same volume-matched runs as numbers instead of curves: for each storm, how big
       the miss is against f, as four scores (peak error, rising-limb steepness ratio,
       recession-rate ratio, KGE_2012). This figure needs ONLY the results CSVs, so it is
       drawn even if the hydrograph files for figure 1 are not available.

It also writes two small tables you can open in Excel:
    series120_hydrograph_shape_runs.csv    one row per volume-matched run per storm
    series120_hydrograph_shape_bands.csv   one row per storm and f group (counts and medians)

FLOW OUTSIDE THE WINDOW
    The top and bottom rows zoom to the main flood: every time any curve is above 1% of the
    truth's peak, joined together when the gaps between them are under 60 minutes, plus 30
    minutes of margin. Any flow above that 1% line outside this window (a late small rise,
    say) is not drawn. The script instead prints, for each storm, how many runs have it, how
    big it is (m3/s and as a percent of the truth's peak) and when, and writes the same two
    numbers per run into series120_hydrograph_shape_runs.csv (flow_outside_window_max_m3s and
    flow_outside_window_time). It does not say WHY the flow is there; that needs a look at the
    rain forcing and the run itself.

WHERE THE HYDROGRAPHS COME FROM
    Every run you scored left a file  calibration_work/03_comparisons/csv_exports/
    <run_id>_compare_obs_sim.csv  with two columns, Observed (here the synthetic truth)
    and Simulated, every 5 minutes for the 20-hour event window. Figure 1 reads those
    files. It does NOT run tRIBS and does not need the big per-run results folders (the
    sweeps delete those after each run; these small CSV files are kept).

WHAT THIS SCRIPT CHECKS BEFORE IT DRAWS (so the picture cannot quietly show the wrong runs)
    For every hydrograph file it recomputes, from the curves themselves, the run's
    PBIAS, peak flow and the truth's peak flow and volume, and compares them with the
    numbers in your results CSV. It also checks that every file for a storm carries the
    SAME truth curve and the same time steps. A file that fails any check is left out of
    figure 1, and the script says which file and why. Nothing is ever changed or deleted.

WHAT THIS SCRIPT DOES NOT DO
    It does not run tRIBS and it does not change any file you already have. It only
    READS the results CSVs and the hydrograph CSVs and WRITES PNGs and two small CSVs
    into its own folder. It never touches current_run_config.json or any run folder, so
    it is safe to run while a sweep is still going.

HOW TO RUN (from the lab/ folder, like your other scripts)
    python plot_series120_hydrograph_shape.py

It looks in calibration_work/03_comparisons/summary_tables/ for
    lhs_results_synth_Ks_f_120_storm080.csv     (0.8x)
    lhs_results_synth_Ks_f_120.csv              (1.0x)
    lhs_results_synth_Ks_f_120_storm125.csv     (1.25x)
and in calibration_work/03_comparisons/csv_exports/ for the hydrograph files, and saves into
    calibration_work/03_comparisons/sensitivity_plots/series120_hydrograph_shape/

OPTIONS (all optional)
    --legs 080,100,125   which storms to use (default all three; one is fine)
    --csv_080 PATH  --csv_100 PATH  --csv_125 PATH   use a different results file
    --summary_dir PATH   look for the results CSVs in a different folder
    --csv_dir PATH       look for the hydrograph files in a different folder
    --out_dir PATH       save somewhere else
    --feasible_pct 2     "volume-matched" means within this many percent of the truth volume
    --f_low 0.007        runs with f below this are drawn blue
    --f_high 0.021       runs with f at or above this are drawn red
    --peak_window_min 35 half-width of the middle (peak zoom) row, in minutes
    --skip_hydrographs   draw only figure 2 (needs no hydrograph files)
    --watermark TEXT     print a faint diagonal word across both figures (e.g. PRELIMINARY)
    --dpi 150            picture sharpness (72 to 600)

HOW TO READ THE SCORES IN FIGURE 2 (they come from run_sensitivity_single_interp_tribs6.py)
    Peak error                  (run peak - truth peak) / truth peak, in percent.
                                Negative = the run's peak is too low.
    Rising-limb steepness ratio the run's rate of rise divided by the truth's, measured from
                                the first time flow passes 5% of the truth's peak up to the
                                peak. 1 = same. Below 1 = the run rises more gently.
    Recession-rate ratio        the run's rate of fall divided by the truth's, on the falling
                                limb above 5% of the truth's peak (a straight-line fit to the
                                log of the flow). 1 = same. Below 1 = the run falls more slowly.
    KGE_2012                    Kling-Gupta efficiency (2012 version), 1 = perfect. Because the
                                runs are volume-matched, what pulls it below 1 is mostly shape.

THINGS TO KEEP IN MIND WHEN READING IT
    * The data are synthetic and noise-free: the truth is a model run, not a gauge. Real
      gauge error and rain-gauge error are NOT in these pictures and could be larger than
      the differences shown here.
    * Routing parameters were held at their truth values in every Series 120 run, so every
      difference is caused by Ks and f. (The script reads the routing values from the CSVs
      and prints them under figure 1; it says so if they were not all the same.)
    * Along the zero-bias ridge, Ks rises as f rises, so f is also a stand-in for "where on
      the ridge". The colour shows f, but each run has its own Ks as well.
    * The f cut-offs (0.007 and 0.021) are choices, not measurements: they sit where the
      KGE_2012 of the volume-matched runs changes in the Series 120 tables. If you change
      them, the colours and the counts change, nothing else.
    * Between f about 0.015 and 0.021 there are few or no volume-matched runs in these
      sweeps, so where the shape starts to degrade at high f is not pinned down.
"""

import argparse
import textwrap
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")          # draw to a file; no window needed (Codespaces has no screen)
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams["font.family"] = "DejaVu Sans"   # ships with matplotlib: same look everywhere

# ======================================================================
# SETTINGS you might want to change
# ======================================================================
TRUTH_KS = 7.0
TRUTH_F = 0.012
FEASIBLE_PCT = 2.0             # "volume-matched" = within 2% of the truth volume
F_LOW = 0.007                  # f below this: the "low f" group (blue)
F_HIGH = 0.021                 # f at or above this: the "high f" group (red)
PEAK_WINDOW_MIN = 35           # the middle row shows the truth's peak time +/- this many minutes
ZOOM_FRAC = 0.01               # rows 1 and 3 zoom to where any curve is above this share of the truth peak ...
PAD_MIN = 30                   # ... plus this many minutes of margin on each side
GAP_MIN = 60                   # curves above ZOOM_FRAC separated by more than this many minutes are separate events;
                               # only the one holding the truth's peak is drawn (the others are reported in words)
EXPECTED_DRAWS = 400
STEP_MIN = 5                   # the hydrograph files are on a 5-minute grid
CHECK_PBIAS_TOL = 1e-3         # file-vs-CSV check: PBIAS may differ by this many percentage points
CHECK_REL_TOL = 1e-5           # file-vs-CSV check: peaks and volume may differ by this share
TRUTH_SAME_TOL = 1e-6          # all files of a storm must carry the same truth to within this share of the peak

LEGS = {
    "080": dict(scale=0.8, name="0.8x storm", note="rain x 0.8", csv="lhs_results_synth_Ks_f_120_storm080.csv"),
    "100": dict(scale=1.0, name="1.0x storm", note="rain as recorded", csv="lhs_results_synth_Ks_f_120.csv"),
    "125": dict(scale=1.25, name="1.25x storm", note="rain x 1.25", csv="lhs_results_synth_Ks_f_120_storm125.csv"),
}
OUT_SUBFOLDER = "series120_hydrograph_shape"
FIG1 = "fig_series120_hydrograph_shape_1_hydrographs.png"
FIG2 = "fig_series120_hydrograph_shape_2_shape_vs_f.png"
RUNS_CSV = "series120_hydrograph_shape_runs.csv"
BANDS_CSV = "series120_hydrograph_shape_bands.csv"

# Colours: the light-mode chart surface and ink of the other Series 120 figures. The three f groups use the
# reference palette's blue and red as the two ends and a mid-grey for "near the truth f" (a diverging idea:
# below / near / above the truth value). I ran the palette checker on blue, grey and red (all pairs, light
# mode): every check passes except its chroma floor, which is expected for a deliberately neutral middle group.
# Marker shape in figure 2 is a second channel (circle / diamond / triangle) so colour is never the only cue.
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SOFT = "#52514e"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
GOLD = "#f0a30a"
BAND_ORDER = ["low", "near", "high"]
BAND_STYLE = {
    "low":  dict(color="#2a78d6", marker="o", lw=1.2, z=3.2),
    "near": dict(color="#898781", marker="D", lw=2.4, z=3.0),
    "high": dict(color="#e34948", marker="^", lw=1.2, z=3.4),
}
F_TICKS = [0.003, 0.005, 0.01, 0.02, 0.03, 0.05]
RUN_FILE_SUFFIX = "_compare_obs_sim.csv"


class DataProblem(Exception):
    """Something is wrong with an input file. Printed as a plain message, without a traceback."""


class FileProblem(Exception):
    """One hydrograph file cannot be used. The message says why."""


# ======================================================================
# Results CSVs: one per storm
# ======================================================================
REQUIRED = ["run_id", "Ks_mult", "f_RS_abs", "pbias_pct", "peak_error_pct", "kge_2012"]
OPTIONAL = ["obs_peak_m3s", "sim_peak_m3s", "obs_volume_m3", "peak_timing_error_hr", "first_arrival_error_min",
            "rising_limb_steepness_ratio", "recession_rate_ratio", "kge_gamma", "kinemvelcoef", "flowexp",
            "channelroughness", "optpercolation", "draw_idx", "obs_peak_time", "sim_peak_time"]


def load_leg(key, path):
    """Read one results CSV and check it really is this storm."""
    spec = LEGS[key]
    path = Path(path)
    if not path.exists():
        raise DataProblem(f"Could not find {path}")
    df = pd.read_csv(path)
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        raise DataProblem(f"{path.name} is missing column(s) {missing}.")
    n_raw = len(df)
    df = df.dropna(subset=["run_id", "Ks_mult", "f_RS_abs", "pbias_pct"]).copy()
    if len(df) < n_raw:
        print(f"  NOTE: {path.name}: dropped {n_raw - len(df)} row(s) with missing values.")
    if df["run_id"].duplicated().any():
        raise DataProblem(f"{path.name} has repeated run_id values: the same run appears twice.")
    if "obs_mode" in df.columns and not (df["obs_mode"] == "synth").all():
        raise DataProblem(f"{path.name}: not every row is obs_mode = synth, so this is not a synthetic-truth file.")
    if "rain_scale" in df.columns:
        vals = df["rain_scale"].astype(float).unique()
        if len(vals) != 1 or abs(vals[0] - spec["scale"]) > 1e-6:
            raise DataProblem(f"{path.name} has rain_scale {[float(v) for v in vals]} but this slot is the {spec['name']} "
                              f"(rain_scale {spec['scale']}). Wrong file for this slot?")
        how = f"rain_scale = {vals[0]:g}"
    elif "storm_leg" in df.columns:
        vals = df["storm_leg"].astype(str).unique()
        if len(vals) != 1 or vals[0] != f"{spec['scale']:.2f}x":
            raise DataProblem(f"{path.name} has storm_leg {list(vals)} but this slot is the {spec['name']}.")
        how = f"storm_leg = {vals[0]}"
    else:
        how = "no rain_scale or storm_leg column, so the storm size could not be confirmed"
    print(f"  {spec['name']:<12} {len(df):3d} runs  ({path.name};  {how})")
    if len(df) < EXPECTED_DRAWS:
        print(f"  WARNING: {spec['name']} has {len(df)} of {EXPECTED_DRAWS} runs. This looks like a sweep that "
              f"is not finished, so everything below is PARTIAL.")
    for c in OPTIONAL:
        if c not in df.columns:
            df[c] = np.nan
    return df.reset_index(drop=True)


def band_of(f, f_low, f_high):
    return "low" if f < f_low else ("high" if f >= f_high else "near")


def pick_matched(df, pct, f_low, f_high):
    """The volume-matched runs of one storm, lowest f first, each tagged with its f group."""
    m = df[df["pbias_pct"].abs() < pct].copy()
    m["f_band"] = [band_of(f, f_low, f_high) for f in m["f_RS_abs"]]
    return m.sort_values("f_RS_abs").reset_index(drop=True)


def partial_text(keys, data):
    short = [f"{LEGS[k]['name']} has {len(data[k])} of {EXPECTED_DRAWS} runs" for k in keys
             if len(data[k]) < EXPECTED_DRAWS]
    return ("PARTIAL RESULTS: " + "; ".join(short) + ". Treat as preliminary.") if short else ""


# ======================================================================
# Hydrograph files: reading and checking
# ======================================================================
def read_run_file(path):
    """Read one <run_id>_compare_obs_sim.csv. Returns (DatetimeIndex, observed array, simulated array)."""
    if not path.exists():
        raise FileProblem("file not found")
    try:
        df = pd.read_csv(path, index_col=0)
    except Exception as exc:                                   # a broken or half-written file
        raise FileProblem(f"could not be read ({str(exc)[:60]})")
    for c in ("Observed", "Simulated"):
        if c not in df.columns:
            raise FileProblem(f"has no '{c}' column (columns found: {list(df.columns)[:4]})")
    try:
        idx = pd.DatetimeIndex(pd.to_datetime(df.index))
    except Exception:
        raise FileProblem("its first column is not a date and time")
    obs = pd.to_numeric(df["Observed"], errors="coerce").to_numpy(float)
    sim = pd.to_numeric(df["Simulated"], errors="coerce").to_numpy(float)
    if len(df) < 10:
        raise FileProblem(f"has only {len(df)} rows")
    if not (np.isfinite(obs).all() and np.isfinite(sim).all()):
        raise FileProblem("has blank or non-numeric flow values")
    steps = np.diff(idx.to_numpy()).astype("timedelta64[s]").astype(float) / 60.0     # minutes between rows
    if not np.allclose(steps, STEP_MIN):
        raise FileProblem(f"time step is not a steady {STEP_MIN} minutes")
    return idx, obs, sim


def check_against_row(row, idx, obs, sim):
    """Recompute the run's numbers from its curves and compare with the results CSV. Returns '' if all agree."""
    problems = []
    pbias = 100.0 * float(np.sum(sim - obs)) / float(np.sum(obs))
    if abs(pbias - float(row["pbias_pct"])) > CHECK_PBIAS_TOL:
        problems.append(f"PBIAS from the file is {pbias:+.4f}% but the results CSV says {float(row['pbias_pct']):+.4f}%")
    if np.isfinite(row["obs_peak_m3s"]) and abs(obs.max() - row["obs_peak_m3s"]) > CHECK_REL_TOL * row["obs_peak_m3s"]:
        problems.append(f"truth peak {obs.max():.4f} differs from the CSV's {row['obs_peak_m3s']:.4f}")
    if np.isfinite(row["sim_peak_m3s"]) and abs(sim.max() - row["sim_peak_m3s"]) > CHECK_REL_TOL * max(1.0, row["sim_peak_m3s"]):
        problems.append(f"run peak {sim.max():.4f} differs from the CSV's {row['sim_peak_m3s']:.4f}")
    if np.isfinite(row["obs_volume_m3"]):
        vol = float(obs.sum()) * STEP_MIN * 60.0
        if abs(vol - row["obs_volume_m3"]) > CHECK_REL_TOL * row["obs_volume_m3"]:
            problems.append(f"truth volume {vol:.1f} m3 differs from the CSV's {row['obs_volume_m3']:.1f} m3")
    return "; ".join(problems)


def load_hydrographs(matched, csv_dir):
    """Read and check the hydrograph file of every volume-matched run of one storm.

    Returns a dict:  status  (list, one text per matched run: 'ok' or the reason it is left out)
                     idx, truth  (the shared time axis and truth curve; None if no file was usable)
                     runs   {row position in `matched`: simulated array}
                     truth_spread  largest difference between any used file's truth and the shared truth
    """
    status = ["" for _ in range(len(matched))]
    good = {}
    for i, row in matched.iterrows():
        path = csv_dir / (str(row["run_id"]) + RUN_FILE_SUFFIX)
        try:
            idx, obs, sim = read_run_file(path)
        except FileProblem as exc:
            status[i] = "file not found" if "not found" in str(exc) else str(exc)
            continue
        why = check_against_row(row, idx, obs, sim)
        if why:
            status[i] = "left out: " + why
            continue
        good[i] = (idx, obs, sim)
    out = dict(status=status, idx=None, truth=None, runs={}, truth_spread=np.nan)
    if not good:
        return out
    # all files of a storm must share one grid and one truth; the most common one wins
    sig = {i: (len(v[0]), str(v[0][0]), str(v[0][-1]), round(float(v[1].sum()), 6), round(float(v[1].max()), 6))
           for i, v in good.items()}
    common = Counter(sig.values()).most_common(1)[0][0]
    ref_i = [i for i in good if sig[i] == common][0]
    ref_idx, ref_obs = good[ref_i][0], good[ref_i][1]
    spread = 0.0
    for i, (idx, obs, sim) in good.items():
        same_grid = len(idx) == len(ref_idx) and bool((idx == ref_idx).all())
        diff = float(np.max(np.abs(obs - ref_obs))) / float(ref_obs.max()) if same_grid else np.inf
        if not same_grid:
            status[i] = "left out: its time steps differ from the other files of this storm"
        elif diff > TRUTH_SAME_TOL:
            status[i] = f"left out: its truth curve differs from the other files of this storm (by {100 * diff:.2f}% of the peak)"
        else:
            out["runs"][i] = sim
            status[i] = "ok"
            spread = max(spread, diff)
    out["idx"], out["truth"], out["truth_spread"] = ref_idx, ref_obs, spread
    return out


# ======================================================================
# Small drawing helpers
# ======================================================================
def style_axes(ax):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.grid(True, color=GRID, linewidth=0.9)
    ax.set_axisbelow(True)
    ax.tick_params(colors=INK_SOFT, labelsize=10.5)


def set_time_axis(ax, lo, hi):
    """x limits and clock-time ticks (hh:mm) that suit the width of the window."""
    ax.set_xlim(lo, hi)
    span_h = (pd.Timestamp(hi) - pd.Timestamp(lo)).total_seconds() / 3600.0
    if span_h <= 1.6:
        loc = mdates.MinuteLocator(byminute=range(0, 60, 15))
    elif span_h <= 4.5:
        loc = mdates.MinuteLocator(byminute=[0, 30])
    elif span_h <= 9:
        loc = mdates.HourLocator()
    else:
        loc = mdates.HourLocator(byhour=range(0, 24, 2))
    ax.xaxis.set_major_locator(loc)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))


def time_label(lo, hi):
    """x-axis label that names the date, or both dates if the window runs past midnight."""
    lo = pd.Timestamp(mdates.num2date(lo).replace(tzinfo=None)) if not isinstance(lo, pd.Timestamp) else lo
    hi = pd.Timestamp(mdates.num2date(hi).replace(tzinfo=None)) if not isinstance(hi, pd.Timestamp) else hi
    if lo.date() == hi.date():
        return f"Time on {lo:%d %b %Y} (hh:mm)"
    return f"Time (hh:mm), {lo:%d %b} to {hi:%d %b %Y}"


def panel_heading(ax, title, sub, sub_lines=2):
    ax.annotate(title, xy=(0, 1), xycoords="axes fraction", xytext=(0, 12 + 14 * sub_lines + 4),
                textcoords="offset points", fontsize=13.5, fontweight="bold", color=INK, va="bottom")
    ax.annotate(sub, xy=(0, 1), xycoords="axes fraction", xytext=(0, 8), textcoords="offset points",
                fontsize=10.5, color=INK_SOFT, va="bottom", linespacing=1.3)


MIN_FIG_W = 13.0               # inches: never narrower than this, so titles and footnotes fit even for one storm


def wrap_text(text, fig_w_in, left_in, fontsize):
    """Wrap every line of `text` so it fits between the left margin and the right edge of the figure."""
    chars = max(40, int((fig_w_in - left_in - 0.15) * 72.0 / (fontsize * 0.56)))
    return "\n".join(textwrap.fill(p, width=chars) if p.strip() else "" for p in text.split("\n"))


def n_lines(text):
    return text.count("\n") + 1


def top_block_in(sub_text, warn, heading_lines):
    """Inches from the top of the figure to the top of the panels: title, subtitle, optional warning, panel headings."""
    title_block = 0.58 + n_lines(sub_text) * 0.21 + (0.30 if warn else 0.0)
    return title_block + (12 + 14 * heading_lines + 4 + 16) / 72.0 + 0.10


def bottom_block_in(notes_text, fontsize):
    """Inches from the bottom of the figure to the bottom of the panels (axis labels, legend, footnotes).
    Returns (that height, the height of the footnote block alone)."""
    notes_h = n_lines(notes_text) * fontsize * 1.45 / 72.0
    return 0.60 + 0.35 + 0.55 + 0.20 + notes_h + 0.15, notes_h


def fig_title(fig, title, sub, left_frac, warn=""):
    """Title, subtitle (already wrapped) and an optional warning line, placed in inches from the top edge."""
    H = fig.get_figheight()
    fig.text(left_frac, 1 - 0.37 / H, title, fontsize=17, fontweight="bold", color=INK, ha="left", va="center")
    fig.text(left_frac, 1 - 0.58 / H, sub, fontsize=11, color=INK_SOFT, ha="left", va="top", linespacing=1.35)
    if warn:
        fig.text(left_frac, 1 - (0.58 + n_lines(sub) * 0.21 + 0.08) / H, warn, fontsize=11.5, fontweight="bold",
                 color=INK, ha="left", va="top")


def add_watermark(fig, text):
    if text:
        fig.text(0.5, 0.5, text, fontsize=64, color=INK, alpha=0.09, rotation=25, ha="center", va="center",
                 fontweight="bold", zorder=100)


def band_handles(f_low, f_high, kind):
    """Legend entries for the three f groups ('line' for figure 1, 'marker' for figure 2)."""
    names = {"low": f"f below {f_low:g}",
             "near": f"f from {f_low:g} to {f_high:g}  (truth f = {TRUTH_F:g})",
             "high": f"f at or above {f_high:g}"}
    hs = []
    for b in BAND_ORDER:
        st = BAND_STYLE[b]
        if kind == "line":
            hs.append(Line2D([], [], color=st["color"], linewidth=max(st["lw"], 2.2), label=names[b]))
        else:
            hs.append(Line2D([], [], marker=st["marker"], color="none", markerfacecolor=st["color"],
                             markeredgecolor=SURFACE, markersize=9, label=names[b]))
    return hs


def routing_text(data, keys):
    """One sentence on the routing parameters, read from the results CSVs (not typed in)."""
    parts, same = [], True
    for col, label in (("kinemvelcoef", "cv"), ("flowexp", "flowexp r"), ("channelroughness", "n")):
        vals = pd.concat([data[k][col] for k in keys]).dropna().round(9).unique()
        if len(vals) == 0:
            return ""
        if len(vals) == 1:
            parts.append(f"{label} {vals[0]:g}")
        else:
            same = False
    if not same:
        return "Routing parameters (cv, flowexp r, n) were NOT the same in every run, so differences here are not only Ks and f."
    loss = pd.concat([data[k]["optpercolation"] for k in keys]).dropna().unique()
    loss_txt = ", channel loss off (optpercolation 0)" if (len(loss) == 1 and loss[0] == 0) else ""
    return "Routing parameters were the same in every run (" + ", ".join(parts) + ")" + loss_txt + \
           ", so these differences come from Ks and f only."


def nice_ceil(x, step):
    return step * np.ceil(x / step)


def save(fig, out_png, dpi):
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=dpi, facecolor=fig.get_facecolor())
    plt.close(fig)


# ======================================================================
# Figure 1: the hydrographs
# ======================================================================
def window_of(idx, truth, sims, frac, pad_min, gap_min=GAP_MIN):
    """Time limits of the MAIN flood: every time any curve is above `frac` of the truth peak, joined into one event
    when the gaps are shorter than `gap_min`; the event that holds the truth's peak is kept, plus a margin.
    A separate small rise hours later does not stretch the window."""
    thr = frac * truth.max()
    on = truth >= thr
    for s in sims:
        on = on | (s >= thr)
    if not on.any():
        return idx[0], idx[-1]
    where = np.flatnonzero(on)
    breaks = np.flatnonzero(np.diff(where) * STEP_MIN > gap_min)
    starts = np.r_[0, breaks + 1]
    ends = np.r_[breaks, len(where) - 1]
    pk = int(np.argmax(truth))
    first, last = where[0], where[-1]
    for a, b in zip(starts, ends):
        if where[a] <= pk <= where[b]:
            first, last = where[a], where[b]
            break
    pad = pd.Timedelta(minutes=pad_min)
    return max(idx[0], idx[first] - pad), min(idx[-1], idx[last] + pad)


def analyse_windows(h):
    """Store the plotting window of one storm and the flow that lies outside it (above ZOOM_FRAC of the truth peak).
    Adds to h:  window (lo, hi);  outside {run position: (largest flow outside, its time)};  truth_outside (same, for the truth)."""
    idx, truth, sims = h["idx"], h["truth"], h["runs"]
    lo, hi = window_of(idx, truth, sims.values(), ZOOM_FRAC, PAD_MIN)
    thr = ZOOM_FRAC * float(truth.max())
    away = (idx < lo) | (idx > hi)

    def biggest(y):
        if not away.any() or float(y[away].max()) < thr:
            return (0.0, None)
        j = int(np.argmax(np.where(away, y, -np.inf)))
        return (float(y[j]), idx[j])

    h["window"] = (lo, hi)
    h["outside"] = {i: biggest(s) for i, s in sims.items()}
    h["truth_outside"] = biggest(truth)


def when_text(t):
    return f"{pd.Timestamp(t):%H:%M} on {pd.Timestamp(t):%d %b}" if t is not None else "-"


def outside_lines(keys, hyd, matched):
    """Facts about flow outside the plotted window, one dict per storm that has any (empty list if none anywhere)."""
    out = []
    for k in keys:
        h = hyd[k]
        if h["truth"] is None or "outside" not in h:
            continue
        P = float(h["truth"].max())
        hit = {i: v for i, v in h["outside"].items() if v[1] is not None}
        t_val, t_time = h["truth_outside"]
        if not hit and t_time is None:
            continue
        groups = Counter(matched[k].loc[i, "f_band"] for i in hit)
        d = dict(name=LEGS[k]["name"], n=len(hit), of=len(h["runs"]), truth_val=t_val, truth_time=t_time, peak=P,
                 groups=", ".join(f"{b} f {groups[b]}" for b in BAND_ORDER if groups.get(b)))
        if hit:
            i_big = max(hit, key=lambda i: hit[i][0])
            d["vmax"], d["tmax"] = hit[i_big]
            times = sorted(pd.Timestamp(x[1]) for x in hit.values())
            d["t_first"], d["t_last"] = times[0], times[-1]
        out.append(d)
    return out


def outside_console(d):
    P = d["peak"]
    if d["n"]:
        span = (f"{d['t_first']:%H:%M} to {d['t_last']:%H:%M}" if d["t_first"] != d["t_last"] else f"{d['t_first']:%H:%M}")
        head = (f"{d['n']} of {d['of']} runs ({d['groups']}); largest {d['vmax']:.2f} m3/s = {100 * d['vmax'] / P:.1f}% of the truth peak, "
                f"at {when_text(d['tmax'])}; the times of the largest values span {span}")
    else:
        head = "none of the runs"
    tr = (f"truth itself: {d['truth_val']:.2f} m3/s at {when_text(d['truth_time'])}" if d["truth_time"] is not None
          else "the truth has none")
    return head, tr


def outside_note(found):
    """One footnote sentence for figure 1 (empty text if there is nothing to say)."""
    if not found:
        return ""
    parts = []
    for d in found:
        if d["n"]:
            parts.append(f"{d['name']}: {d['n']} of {d['of']} runs ({d['groups']}), up to {100 * d['vmax'] / d['peak']:.1f}% of the truth "
                         f"peak near {d['tmax']:%H:%M} on {d['tmax']:%d %b}")
        else:
            parts.append(f"{d['name']}: the truth only")
    tru = [d for d in found if d["truth_time"] is not None]
    tail = ("; the truth itself has some there too (see the console)" if tru else "; the truth has none there")
    return ("Not drawn: a small extra rise away from the main flood (above 1% of the truth peak). " + "; ".join(parts) + tail + ".")


def figure_hydrographs(keys, matched, hyd, data, csv_dir, f_low, f_high, args, out_png):
    n = len(keys)
    W = max(5.4 * n + 1.3, MIN_FIG_W)
    left_in = 1.3
    left = left_in / W
    n_all = sum(len(matched[k]) for k in keys)
    n_ok = sum(len(hyd[k]["runs"]) for k in keys)
    warn = partial_text(keys, data)
    sub_text = wrap_text(
        f"Black = the synthetic truth (Ks {TRUTH_KS:g}, f {TRUTH_F:g}). Each coloured line is one run whose total volume is within "
        f"{args.feasible_pct:g}% of the truth's.\nThe colour is the run's f. Bottom row: the run's flow minus the truth's, "
        f"as a percent of the truth's peak (zero = exactly on the truth).", W, left_in, 11)
    notes = []
    rt = routing_text(data, keys)
    if rt:
        notes.append(rt)
    spreads = [hyd[k]["truth_spread"] for k in keys if hyd[k]["truth"] is not None]
    notes.append(f"{n_ok} of {n_all} hydrograph files read from {_short(csv_dir)} and checked against the results CSVs "
                 f"(each storm's files carry one identical truth curve; largest difference "
                 f"{(max(spreads) if spreads else 0.0):.1g} of the peak).")
    on = outside_note(outside_lines(keys, hyd, matched))
    if on:
        notes.append(on)
    notes.append("Along the zero-bias ridge Ks rises with f, so each run also has its own Ks. Noise-free synthetic data: real gauge and rain-gauge error "
                 "are not included.")
    notes_text = wrap_text("\n".join(notes), W, left_in, 9.6)
    top_in = top_block_in(sub_text, warn, 3)
    bottom_in, notes_h = bottom_block_in(notes_text, 9.6)
    H = top_in + 7.2 + bottom_in
    fig = plt.figure(figsize=(W, H))
    fig.patch.set_facecolor(SURFACE)
    gs = fig.add_gridspec(3, n, height_ratios=[1.15, 1.0, 1.0], left=left, right=0.985, top=1 - top_in / H,
                          bottom=bottom_in / H, hspace=0.40, wspace=0.20)

    # one shared y range for the bottom row (percent of each storm's own truth peak)
    max_abs = 3.0
    for k in keys:
        h = hyd[k]
        if h["truth"] is None:
            continue
        P = h["truth"].max()
        for s in h["runs"].values():
            max_abs = max(max_abs, float(np.max(np.abs(s - h["truth"]))) / P * 100.0)
    res_lim = nice_ceil(max_abs * 1.12, 1.0)

    for col, k in enumerate(keys):
        spec, h, m = LEGS[k], hyd[k], matched[k]
        ax1 = fig.add_subplot(gs[0, col])
        ax2 = fig.add_subplot(gs[1, col])
        ax3 = fig.add_subplot(gs[2, col], sharex=ax1)
        for ax in (ax1, ax2, ax3):
            style_axes(ax)
        used = {b: sum(1 for i in h["runs"] if m.loc[i, "f_band"] == b) for b in BAND_ORDER}
        n_used = len(h["runs"])
        line_runs = f"{len(m)} runs within {args.feasible_pct:g}% of the truth volume"
        line_used = (f"in the plot: {used['low']} low f, {used['near']} near, {used['high']} high f"
                     + ("" if n_used == len(m) else f"  ({len(m) - n_used} left out, see the table)"))
        if h["truth"] is None:
            panel_heading(ax1, f"{spec['name']}  ({spec['note']})", line_runs + "\n" + line_used)
            for ax in (ax1, ax2, ax3):
                ax.set_xticks([])
                ax.set_yticks([])
                ax.text(0.5, 0.5, "no usable hydrograph files\nfor these runs in\n" + _short(csv_dir), ha="center",
                        va="center", transform=ax.transAxes, color=INK_SOFT, fontsize=11)
            continue
        idx, truth = h["idx"], h["truth"]
        P = float(truth.max())
        t_pk = idx[int(np.argmax(truth))]
        sims = {i: s for i, s in h["runs"].items()}
        panel_heading(ax1, f"{spec['name']}  ({spec['note']})",
                      f"truth peak {P:.1f} m³/s at {t_pk:%H:%M}\n{line_runs}\n{line_used}", sub_lines=3)
        x = idx.to_numpy()
        lo, hi = h["window"]
        # draw order: near (thick grey) first, then low and high, so coloured lines stay visible; truth last
        order = sorted(sims, key=lambda i: (BAND_STYLE[m.loc[i, "f_band"]]["z"]))
        for i in order:
            st = BAND_STYLE[m.loc[i, "f_band"]]
            s = sims[i]
            for ax, y in ((ax1, s), (ax2, s), (ax3, (s - truth) / P * 100.0)):
                ax.plot(x, y, color=st["color"], linewidth=st["lw"], alpha=0.85, zorder=st["z"], solid_capstyle="round")
        for ax, y in ((ax1, truth), (ax2, truth)):
            ax.plot(x, y, color=INK, linewidth=1.5, zorder=6)
        ax3.axhline(0, color=INK, linewidth=1.5, zorder=6)

        # row 1: the whole flood wave
        top = max([float(truth.max())] + [float(s.max()) for s in sims.values()])
        ax1.set_ylim(0, top * 1.07)
        set_time_axis(ax1, lo, hi)
        ax1.set_ylabel("Flow (m³/s)", color=INK_SOFT, fontsize=11.5)

        # row 2: zoom on the top of the peak
        w = pd.Timedelta(minutes=args.peak_window_min)
        set_time_axis(ax2, max(idx[0], t_pk - w), min(idx[-1], t_pk + w))
        pk_sims = [float(s.max()) for s in sims.values()]
        ymin = min(0.80 * P, 0.92 * min(pk_sims))             # the zoom starts at about 80% of the peak
        ymax = 1.05 * max([P] + pk_sims)
        ax2.set_ylim(max(ymin, 0.0), ymax)
        ax2.set_ylabel("Flow near the peak (m³/s)", color=INK_SOFT, fontsize=11.5)
        ax2.set_xlabel(time_label(ax2.get_xlim()[0], ax2.get_xlim()[1]), color=INK_SOFT, fontsize=10.5)

        # row 3: the miss
        ax3.set_ylim(-res_lim, res_lim)
        set_time_axis(ax3, lo, hi)
        ax3.axvline(t_pk, color=AXIS, linewidth=1.0, linestyle=(0, (2, 3)), zorder=1)
        ax3.set_ylabel("Run minus truth\n(% of truth peak)", color=INK_SOFT, fontsize=11.5)
        ax3.set_xlabel(time_label(lo, hi), color=INK_SOFT, fontsize=10.5)

    fig_title(fig, "Do runs that match the truth volume also match its shape?  tRIBS 6.0.0, Series 120", sub_text, left, warn=warn)
    handles = [Line2D([], [], color=INK, linewidth=2.0, label=f"synthetic truth (Ks {TRUTH_KS:g}, f {TRUTH_F:g})")]
    handles += band_handles(f_low, f_high, "line")
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=10.5, labelcolor=INK_SOFT,
               bbox_to_anchor=(0.5, (0.15 + notes_h + 0.20) / H))
    fig.text(left, (0.15 + notes_h) / H, notes_text, fontsize=9.6, color=INK_SOFT, ha="left", va="top", linespacing=1.45)
    add_watermark(fig, args.watermark)
    save(fig, out_png, args.dpi)


def _short(p):
    """A folder path shortened to the part you recognise (from calibration_work on)."""
    parts = Path(p).parts
    return "/".join(parts[parts.index("calibration_work"):]) if "calibration_work" in parts else str(p)


# ======================================================================
# Figure 2: the same runs as scores against f
# ======================================================================
SCORE_ROWS = [      # (column in the results CSV, y-axis label, the value of a perfect match)
    ("peak_error_pct", "Peak error (%)\n0 = same peak as the truth", 0.0),
    ("rising_limb_steepness_ratio", "Rising-limb steepness ratio\n1 = same rise as the truth", 1.0),
    ("recession_rate_ratio", "Recession-rate ratio\n1 = same fall as the truth", 1.0),
    ("kge_2012", "KGE_2012\n1 = perfect", 1.0),
]


def figure_scores(keys, matched, data, f_low, f_high, args, out_png):
    n = len(keys)
    nr = len(SCORE_ROWS)
    W = max(5.4 * n + 1.3, MIN_FIG_W)
    left_in = 2.0
    left = left_in / W
    warn = partial_text(keys, data)
    sub_text = wrap_text(
        f"Each marker is one run within {args.feasible_pct:g}% of its storm's truth volume (the same runs as figure 1). "
        f"Black line = a perfect match. Gold dashed line = the truth f ({TRUTH_F:g}).", W, left_in, 11)
    notes = [
        "Peak error = (run peak - truth peak) / truth peak; negative = the run's peak is too low.",
        "Rising-limb ratio = run's rate of rise / truth's, from first flow above 5% of the truth peak up to the peak; below 1 = gentler rise.",
        "Recession-rate ratio = run's rate of fall / truth's on the falling limb above 5% of the truth peak; below 1 = slower fall.",
        "KGE_2012 = Kling-Gupta efficiency, 1 = perfect (these runs have the volume right, so what lowers it is mostly shape).",
        "Noise-free synthetic data: real gauge and rain-gauge error are not included. The dotted grey lines are the f cut-offs between the three groups; "
        "they are choices, not measurements.",
    ]
    notes_text = wrap_text("\n".join(notes), W, left_in, 9.4)
    top_in = top_block_in(sub_text, warn, 2)
    bottom_in, notes_h = bottom_block_in(notes_text, 9.4)
    H = top_in + 9.0 + bottom_in
    fig, axes = plt.subplots(nr, n, figsize=(W, H), sharex=True, sharey="row", squeeze=False)
    fig.patch.set_facecolor(SURFACE)
    fig.subplots_adjust(left=left, right=0.985, top=1 - top_in / H, bottom=bottom_in / H, hspace=0.12, wspace=0.07)
    for col, k in enumerate(keys):
        spec, m = LEGS[k], matched[k]
        counts = {b: int((m["f_band"] == b).sum()) for b in BAND_ORDER}
        for r, (colname, ylabel, ref) in enumerate(SCORE_ROWS):
            ax = axes[r, col]
            style_axes(ax)
            ax.set_xscale("log")
            ax.set_xlim(0.003 / 1.07, 0.05 * 1.10)
            ax.set_xticks(F_TICKS)
            ax.set_xticklabels([f"{t:g}" for t in F_TICKS])
            ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
            for edge in (f_low, f_high):
                ax.axvline(edge, color=AXIS, linewidth=1.0, linestyle=(0, (2, 3)), zorder=1)
            ax.axvline(TRUTH_F, color=GOLD, linewidth=1.6, linestyle=(0, (5, 3)), zorder=2)
            ax.axhline(ref, color=INK, linewidth=1.1, zorder=2)
            for b in BAND_ORDER:
                sub = m[m["f_band"] == b]
                st = BAND_STYLE[b]
                ax.scatter(sub["f_RS_abs"], sub[colname], s=58, marker=st["marker"], color=st["color"],
                           edgecolors=SURFACE, linewidths=0.7, zorder=4)
            if col == 0:
                ax.set_ylabel(ylabel, color=INK_SOFT, fontsize=10.5)
            else:
                ax.tick_params(labelleft=False)
            if r < nr - 1:
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xlabel("f (1/mm, log scale)", color=INK_SOFT, fontsize=11.5)
            if np.isfinite(m[colname]).sum() < len(m):
                ax.text(0.99, 0.03, f"{int((~np.isfinite(m[colname])).sum())} run(s) have no value",
                        transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5, color=INK_SOFT)
        panel_heading(axes[0, col], f"{spec['name']}  ({spec['note']})",
                      f"{len(m)} runs within {args.feasible_pct:g}% of the truth volume\n"
                      f"{counts['low']} low f, {counts['near']} near, {counts['high']} high f", sub_lines=2)
    fig_title(fig, "How far the volume-matched runs miss the shape, by f: tRIBS 6.0.0, Series 120", sub_text, left, warn=warn)
    handles = band_handles(f_low, f_high, "marker")
    handles.append(Line2D([], [], color=GOLD, linewidth=1.8, linestyle=(0, (5, 3)), label=f"truth f ({TRUTH_F:g})"))
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=10.5, labelcolor=INK_SOFT,
               bbox_to_anchor=(0.5, (0.15 + notes_h + 0.20) / H))
    fig.text(left, (0.15 + notes_h) / H, notes_text, fontsize=9.4, color=INK_SOFT, ha="left", va="top", linespacing=1.45)
    add_watermark(fig, args.watermark)
    save(fig, out_png, args.dpi)


# ======================================================================
# Tables
# ======================================================================
def write_tables(keys, matched, hyd, out_dir, hydro_asked):
    rows = []
    for k in keys:
        m = matched[k]
        for i, r in m.iterrows():
            outside = hyd[k].get("outside", {}).get(i) if hydro_asked else None
            if not hydro_asked:
                st = "not read (--skip_hydrographs)"
            else:
                st = hyd[k]["status"][i] or "not read"
            rows.append({"storm": LEGS[k]["name"], "run_id": r["run_id"], "Ks_mult": r["Ks_mult"], "f_RS_abs": r["f_RS_abs"],
                         "f_group": r["f_band"], "pbias_pct": r["pbias_pct"], "peak_error_pct": r["peak_error_pct"],
                         "peak_timing_error_hr": r["peak_timing_error_hr"],
                         "first_arrival_error_min": r["first_arrival_error_min"],
                         "rising_limb_steepness_ratio": r["rising_limb_steepness_ratio"],
                         "recession_rate_ratio": r["recession_rate_ratio"], "kge_2012": r["kge_2012"],
                         "kge_gamma": r["kge_gamma"], "hydrograph_file": st,
                         "flow_outside_window_max_m3s": outside[0] if outside else np.nan,
                         "flow_outside_window_time": (str(outside[1]) if outside and outside[1] is not None else "")})
    runs = pd.DataFrame(rows)
    runs.to_csv(out_dir / RUNS_CSV, index=False)
    brows = []
    for (storm, grp), g in runs.groupby(["storm", "f_group"], sort=False):
        brows.append({"storm": storm, "f_group": grp, "n_runs": len(g), "f_min": g["f_RS_abs"].min(), "f_max": g["f_RS_abs"].max(),
                      "peak_error_pct_median": g["peak_error_pct"].median(), "peak_error_pct_min": g["peak_error_pct"].min(),
                      "peak_error_pct_max": g["peak_error_pct"].max(),
                      "rising_limb_ratio_median": g["rising_limb_steepness_ratio"].median(),
                      "recession_rate_ratio_median": g["recession_rate_ratio"].median(),
                      "kge_2012_median": g["kge_2012"].median(), "kge_2012_min": g["kge_2012"].min()})
    bands = pd.DataFrame(brows)
    order = {b: j for j, b in enumerate(BAND_ORDER)}
    if len(bands):
        bands = bands.sort_values(["storm", "f_group"], key=lambda s: s.map(order) if s.name == "f_group" else s,
                                  kind="stable").reset_index(drop=True)
    bands.to_csv(out_dir / BANDS_CSV, index=False)
    return runs, bands


# ======================================================================
def main():
    global FEASIBLE_PCT
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--legs", default="080,100,125", help="comma-separated storms to use (default 080,100,125)")
    for k in LEGS:
        ap.add_argument(f"--csv_{k}", type=Path, default=None)
    ap.add_argument("--summary_dir", type=Path, default=None)
    ap.add_argument("--csv_dir", type=Path, default=None, help="folder holding the <run_id>_compare_obs_sim.csv files")
    ap.add_argument("--out_dir", type=Path, default=None)
    ap.add_argument("--feasible_pct", type=float, default=FEASIBLE_PCT)
    ap.add_argument("--f_low", type=float, default=F_LOW)
    ap.add_argument("--f_high", type=float, default=F_HIGH)
    ap.add_argument("--peak_window_min", type=float, default=PEAK_WINDOW_MIN)
    ap.add_argument("--skip_hydrographs", action="store_true", help="draw only figure 2 (needs no hydrograph files)")
    ap.add_argument("--watermark", default="", help="a faint word across both figures, e.g. PRELIMINARY")
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args()
    if not (72 <= args.dpi <= 600):
        ap.error("--dpi must be between 72 and 600.")
    if not (0.1 <= args.feasible_pct <= 20.0):
        ap.error("--feasible_pct must be between 0.1 and 20.")
    if not (0.003 < args.f_low < args.f_high < 0.05):
        ap.error("--f_low and --f_high must satisfy 0.003 < f_low < f_high < 0.05.")
    if not (10 <= args.peak_window_min <= 120):
        ap.error("--peak_window_min must be between 10 and 120.")
    FEASIBLE_PCT = args.feasible_pct

    keys = [k.strip() for k in args.legs.split(",") if k.strip()]
    bad = [k for k in keys if k not in LEGS]
    if bad or len(set(keys)) != len(keys) or not keys:
        ap.error(f"--legs must be a comma-separated list drawn from {list(LEGS)} with no repeats (got {args.legs!r}).")
    keys = [k for k in LEGS if k in keys]                       # smallest rain first, whatever order was typed

    calib_dir = Path.cwd().parent / "calibration_work"
    summary_dir = args.summary_dir or (calib_dir / "03_comparisons" / "summary_tables")
    csv_dir = args.csv_dir or (calib_dir / "03_comparisons" / "csv_exports")
    out_dir = args.out_dir or (calib_dir / "03_comparisons" / "sensitivity_plots" / OUT_SUBFOLDER)
    # never write into a folder that holds inputs
    out_res = out_dir.resolve()
    for guarded in (csv_dir, summary_dir):
        g = guarded.resolve()
        if out_res == g:
            print(f"The output folder {out_dir} is the same as an input folder ({guarded}); refusing to write there.")
            raise SystemExit(1)

    paths = {k: (getattr(args, f"csv_{k}") or summary_dir / LEGS[k]["csv"]) for k in keys}
    missing = [p for p in paths.values() if not Path(p).exists()]
    if missing:
        print("Could not find these results file(s):")
        for p in missing:
            print(f"    {p}")
        print("Run this from the lab/ folder, or point to the file(s) with --csv_080 / --csv_100 / --csv_125 or --summary_dir.")
        raise SystemExit(1)

    print(f"Loading {len(keys)} storm(s) from {summary_dir}:")
    try:
        data = {k: load_leg(k, paths[k]) for k in keys}
    except DataProblem as exc:
        print(f"\nSTOPPED: {exc}")
        raise SystemExit(1)

    matched = {k: pick_matched(data[k], args.feasible_pct, args.f_low, args.f_high) for k in keys}
    print(f"\nVolume-matched runs (|PBIAS| < {args.feasible_pct:g}%):")
    for k in keys:
        m = matched[k]
        c = m["f_band"].value_counts()
        print(f"  {LEGS[k]['name']:<12} {len(m):3d} runs   low f {int(c.get('low', 0))}, near {int(c.get('near', 0))}, high f {int(c.get('high', 0))}")
    empty = [LEGS[k]["name"] for k in keys if len(matched[k]) == 0]
    if empty:
        print(f"\nSTOPPED: no volume-matched runs for {', '.join(empty)}. Try a larger --feasible_pct.")
        raise SystemExit(1)

    # ---- the hydrograph files ----
    hydro_asked = not args.skip_hydrographs
    hyd = {k: dict(status=["" for _ in range(len(matched[k]))], idx=None, truth=None, runs={}, truth_spread=np.nan) for k in keys}
    if hydro_asked:
        print(f"\nReading the hydrograph files from {csv_dir}:")
        if not csv_dir.exists():
            print("  NOTE: that folder does not exist. Figure 1 will be left out. "
                  "(Point to the right folder with --csv_dir, or use --skip_hydrographs to silence this.)")
        else:
            for k in keys:
                hyd[k] = load_hydrographs(matched[k], csv_dir)
                st = pd.Series(hyd[k]["status"])
                n_ok = int((st == "ok").sum())
                n_nf = int((st == "file not found").sum())
                print(f"  {LEGS[k]['name']:<12} {n_ok:3d} of {len(st)} usable   ({n_nf} not found, {len(st) - n_ok - n_nf} left out)")
                shown = 0
                for i, s in enumerate(hyd[k]["status"]):
                    if s not in ("ok", "file not found") and shown < 8:
                        print(f"      {matched[k].loc[i, 'run_id']}: {s}")
                        shown += 1
                if int((st != "ok").sum()) - n_nf > shown:
                    print(f"      ... and {int((st != 'ok').sum()) - n_nf - shown} more (all listed in {RUNS_CSV})")
                if n_nf and n_nf == len(st):
                    print(f"      None of this storm's files were found. Looking for names like {matched[k].loc[0, 'run_id']}{RUN_FILE_SUFFIX}")
                if hyd[k]["truth"] is not None:
                    analyse_windows(hyd[k])
        found = outside_lines(keys, hyd, matched)
        if found:
            print(f"\nFlow outside the plotted window (above {100 * ZOOM_FRAC:g}% of the truth peak, more than {GAP_MIN} min away from the main flood).")
            print("  It is not drawn in figure 1; it is listed in " + RUNS_CSV + ".")
            for d in found:
                head, tr = outside_console(d)
                print(f"  {d['name']:<12} {head}")
                print(f"  {'':<12} {tr}")
        elif any(hyd[k]["truth"] is not None for k in keys):
            print(f"\nFlow outside the plotted window: none (no curve rises above {100 * ZOOM_FRAC:g}% of the truth peak away from the main flood).")

    # ---- figures and tables ----
    out_dir.mkdir(parents=True, exist_ok=True)
    wrote = []
    n_any = sum(len(hyd[k]["runs"]) for k in keys)
    if hydro_asked and n_any > 0:
        figure_hydrographs(keys, matched, hyd, data, csv_dir, args.f_low, args.f_high, args, out_dir / FIG1)
        wrote.append(FIG1)
    elif hydro_asked:
        print("\nFigure 1 was NOT drawn: no usable hydrograph files were found, so there is nothing to plot.")
        print("  Check where your <run_id>_compare_obs_sim.csv files are and pass that folder with --csv_dir.")
    figure_scores(keys, matched, data, args.f_low, args.f_high, args, out_dir / FIG2)
    wrote.append(FIG2)
    runs, bands = write_tables(keys, matched, hyd, out_dir, hydro_asked)
    wrote += [RUNS_CSV, BANDS_CSV]

    print("\nMedian scores of the volume-matched runs, by storm and f group (full table in " + BANDS_CSV + "):")
    show = bands[["storm", "f_group", "n_runs", "peak_error_pct_median", "rising_limb_ratio_median",
                  "recession_rate_ratio_median", "kge_2012_median"]].copy()
    show.columns = ["storm", "f group", "runs", "peak err %", "rise ratio", "recession ratio", "KGE_2012"]
    print(show.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print(f"\nSaved to {out_dir}:")
    for w in wrote:
        print(f"    {w}")
    if hydro_asked and n_any > 0 and any(s not in ("ok",) for k in keys for s in hyd[k]["status"]):
        print("\nNOTE: some hydrograph files were not used in figure 1 (reasons above; per-run list in " + RUNS_CSV + ").")
    print("\nREMINDER: synthetic, noise-free data. The scores describe what these runs do against a model-made truth, not against a real gauge.")


if __name__ == "__main__":
    main()