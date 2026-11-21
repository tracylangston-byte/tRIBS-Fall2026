#!/usr/bin/env python
"""
plot_storm_overlay_110.py   (version 1)
=======================================
Draws the real gauge hydrograph for the 2014 storms with model runs laid on top, so you can
SEE what "within a factor of 2" looks like. READ-ONLY: it changes nothing you already have.
It writes two NEW files into a NEW folder, and refuses to overwrite them unless you ask.

WHAT IS DRAWN (one picture, three panels: 12 Aug, 8 Sep, 27 Sep)
  * black dots, joined by a thin line : the real gauge readings, as the logger recorded them
  * thin gray lines                   : every run that is "in the ballpark" on BOTH 12 Aug and
                                        8 Sep (simulated peak / real peak AND simulated volume /
                                        real volume both between 0.5 and 2.0, and a peak of at
                                        least 1 cfs). The 0.5-2.0 band is the author's own choice.
  * coloured lines                    : the best combined run, the best run on 12 Aug, the best
                                        run on 8 Sep and the best run on 27 Sep (if any run flows)
                                        by KGE_2012 in the scores file.

HOW TO RUN (from workspaces/SMF/lab/):
    python plot_storm_overlay_110.py --label ms2 --dry_run      # checks only, writes nothing
    python plot_storm_overlay_110.py --label ms2                # draws the figure
    python plot_storm_overlay_110.py --label ms1 --dry_run
    python plot_storm_overlay_110.py --label ms1

OPTIONS
    --label LABEL     ms1 or ms2 (required): which sweep's scores file to use
    --storms LIST     storms to draw, comma-separated (default 12Aug,8Sep,27Sep; 19Aug also works)
    --units cfs|cms   units on the y axis (default cfs, the gauge's own unit)
    --gauge PATH      the gauge workbook (default ../init_data/met/SMF_Observations_1993-2025.xlsx)
    --overwrite       replace the two output files if they already exist
    --dry_run         run every check, print what would be drawn, write nothing

READS (all existing files, never changed)
    calibration_work/03_comparisons/summary_tables/storm_scores_<label>_110/storm_scores_filled_110.csv
    calibration_work/kept_qout_110/<run_id>...Outlet.qout      (the model's outlet flow, one per run)
    ../init_data/met/SMF_Observations_1993-2025.xlsx           (sheet "Discharge")
WRITES (new files only)
    calibration_work/03_comparisons/sensitivity_plots/storm_overlay_110/fig_storm_overlay_<label>_110.png
    calibration_work/03_comparisons/sensitivity_plots/storm_overlay_110/storm_overlay_runs_<label>_110.csv

CHECKS (the script STOPS with a plain message if one fails)
    1. the scores file exists, has the columns it needs, and every run name belongs to --label
    2. the gauge workbook has readings in every storm window, and its highest reading in each window
       equals the real peak the Series 110 record gives (1586, 13, 442, 252 cfs; 2% tolerance)
    3. every run to be drawn has an outlet file that can be read and covers the storm window
    A softer CHECK (printed, does not stop): the simulated peak recomputed from the outlet file is
    compared with the peak in the scores file. They can differ by a few percent, because the scores
    use a 5-minute grid and the outlet file steps every 3.75 minutes.

STATUS: TESTED ON STAND-IN DATA ONLY (made-up outlet files and a made-up gauge workbook, plus your
real ms1 and ms2 scores files). NOT yet run on your real outlet files or your real gauge workbook.
Caveats that go with the picture: same two rain gauges, one season, routing pinned, Kchan on.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.lines import Line2D

# ----------------------------------------------------------------------------------------------
# Settings (all in one place)
# ----------------------------------------------------------------------------------------------
CFS_TO_CMS = 0.0283168
MODEL_ORIGIN = pd.Timestamp("2014-08-01 00:00")   # outlet-file times are hours from here
BAND_LO, BAND_HI = 0.5, 2.0                         # the "ballpark" band (the author's choice)
FLOW_MIN_CFS = 1.0                                  # "made flow" = peak of at least 1 cfs
REAL_PEAK_TOL = 0.02                                # gauge check tolerance (2%)
SIM_PEAK_TOL = 0.05                                 # soft check tolerance (5%)

# storm tag -> (window start, window end, label, real peak in cfs from the Series 110 record)
STORMS = {
    "12Aug": ("2014-08-12 16:00", "2014-08-13 12:00", "12 Aug 2014", 1586.0),
    "19Aug": ("2014-08-19 05:00", "2014-08-19 17:00", "19 Aug 2014", 13.0),
    "8Sep":  ("2014-09-08 00:00", "2014-09-08 16:00", "8 Sep 2014", 442.0),
    "27Sep": ("2014-09-27 12:00", "2014-09-27 23:00", "27 Sep 2014", 252.0),
}

NAVY, BLUE, RED, ORANGE, GRAY, INK = "#0C234B", "#378DBD", "#AB0520", "#E08A00", "#8A8A8A", "#222222"


def die(msg):
    print("\nSTOPPED: " + msg)
    sys.exit(1)


# ----------------------------------------------------------------------------------------------
# Reading the three kinds of input
# ----------------------------------------------------------------------------------------------
def read_scores(path, label, storms):
    if not path.exists():
        die(f"The scores file was not found:\n  {path}\nRun this from workspaces/SMF/lab/ and check --label.")
    df = pd.read_csv(path)
    need = ["run_id", "Ks_mult", "f_RS_abs", "channelconductivity_mmhr", "combined_ALL"]
    for tg in ["12Aug", "8Sep"] + list(storms):
        for c in ("sim_peak_cfs", "made_flow", "peak_ratio", "vol_ratio", "peak_dt_min", "kge_2012"):
            need.append(f"{tg}_{c}")
    missing = sorted(set(c for c in need if c not in df.columns))
    if missing:
        die("The scores file is missing columns this script needs: " + ", ".join(missing[:8]))
    bad = [r for r in df["run_id"] if f"_{label}on" not in r]
    if bad:
        die(f"{len(bad)} run names in the scores file do not contain '_{label}on' (first: {bad[0]}). "
            f"Is --label {label} the right sweep?")
    return df


def read_gauge(path):
    if not path.exists():
        die(f"The gauge workbook was not found:\n  {path}\nUse --gauge to point at it.")
    g = pd.read_excel(path, sheet_name="Discharge", skiprows=6)
    g["datetime"] = pd.to_datetime(g["Date"].astype(str) + " " + g["Time"].astype(str))
    s = pd.Series(pd.to_numeric(g["cfs"], errors="coerce").to_numpy(), index=pd.DatetimeIndex(g["datetime"]))
    s = s.dropna().sort_index()
    return s


def find_qout(kept_dir, run_id):
    hits = sorted(kept_dir.glob(run_id + "*.qout"))
    if len(hits) > 1:
        narrowed = [h for h in hits if "Outlet" in h.name]
        hits = narrowed or hits
    if not hits:
        examples = [p.name for p in sorted(kept_dir.glob("*.qout"))[:3]]
        die(f"No outlet file found for run {run_id} in\n  {kept_dir}\n"
            f"(files that ARE there start like: {examples})")
    return hits[0]


def read_qout(path):
    """The outlet file is comma-delimited with a header and 3 columns (time in hours from
    2014-08-01 00:00, discharge in m3/s, stage). Returns discharge in m3/s indexed by time."""
    df = pd.read_csv(path)
    if df.shape[1] < 3:
        df = pd.read_csv(path, sep=r"\s+")
    if df.shape[1] < 3:
        die(f"Outlet file {path.name} does not have 3 columns, so it is not the format this script expects.")
    t = pd.to_numeric(df.iloc[:, 0], errors="coerce")
    q = pd.to_numeric(df.iloc[:, 1], errors="coerce")
    ok = t.notna() & q.notna()
    idx = pd.DatetimeIndex(MODEL_ORIGIN + pd.to_timedelta(t[ok].to_numpy() * 3600.0, unit="s")).round("s")
    return pd.Series(q[ok].to_numpy(), index=idx).sort_index()


# ----------------------------------------------------------------------------------------------
# Choosing which runs to draw
# ----------------------------------------------------------------------------------------------
def in_band(df, tg):
    return (df[f"{tg}_made_flow"].astype(bool)
            & df[f"{tg}_peak_ratio"].between(BAND_LO, BAND_HI)
            & df[f"{tg}_vol_ratio"].between(BAND_LO, BAND_HI))


def pick_runs(df, storms):
    """Returns (ensemble_ids, roles) where roles maps run_id -> list of role names."""
    ens_mask = in_band(df, "12Aug") & in_band(df, "8Sep")
    ensemble = list(df.loc[ens_mask, "run_id"])
    roles = {}

    def add(run_id, role):
        roles.setdefault(run_id, []).append(role)

    add(df.loc[df["combined_ALL"].idxmax(), "run_id"], "best combined")
    add(df.loc[df["12Aug_kge_2012"].idxmax(), "run_id"], "best on 12 Aug")
    add(df.loc[df["8Sep_kge_2012"].idxmax(), "run_id"], "best on 8 Sep")
    if "27Sep" in storms and df["27Sep_kge_2012"].notna().any():
        add(df.loc[df["27Sep_kge_2012"].idxmax(), "run_id"], "best on 27 Sep")
    return ensemble, roles


STYLE = {   # role -> (colour, linestyle, linewidth)
    "best combined": (NAVY, "-", 2.4),
    "best on 12 Aug": (BLUE, "--", 1.9),
    "best on 8 Sep": (RED, "-.", 1.9),
    "best on 27 Sep": (ORANGE, ":", 2.2),
}


def style_for(role_list):
    for r in ("best combined", "best on 12 Aug", "best on 8 Sep", "best on 27 Sep"):
        if r in role_list:
            return STYLE[r]
    return (GRAY, "-", 0.8)


# ----------------------------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Gauge-vs-model overlay for the Series 110 real storms.")
    ap.add_argument("--label", required=True, help="ms1 or ms2")
    ap.add_argument("--storms", default="12Aug,8Sep,27Sep")
    ap.add_argument("--units", default="cfs", choices=["cfs", "cms"])
    ap.add_argument("--gauge", default=None)
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--dry_run", action="store_true")
    args = ap.parse_args()

    storms = [s.strip() for s in args.storms.split(",") if s.strip()]
    for s in storms:
        if s not in STORMS:
            die(f"Unknown storm '{s}'. Choose from {', '.join(STORMS)}.")
    for s in ("12Aug", "8Sep"):
        if s not in storms:
            die("--storms must include 12Aug and 8Sep (they define the factor-of-2 set).")

    lab_dir = Path.cwd()
    ws_dir = lab_dir.parent
    calib = ws_dir / "calibration_work"
    scores_path = calib / "03_comparisons" / "summary_tables" / f"storm_scores_{args.label}_110" / "storm_scores_filled_110.csv"
    kept_dir = calib / "kept_qout_110"
    gauge_path = Path(args.gauge) if args.gauge else ws_dir / "init_data" / "met" / "SMF_Observations_1993-2025.xlsx"
    out_dir = calib / "03_comparisons" / "sensitivity_plots" / "storm_overlay_110"
    fig_path = out_dir / f"fig_storm_overlay_{args.label}_110.png"
    csv_path = out_dir / f"storm_overlay_runs_{args.label}_110.csv"

    print("=" * 78)
    print(f"Storm overlay, sweep {args.label}   ({'DRY RUN: nothing will be written' if args.dry_run else 'drawing'})")
    print("=" * 78)

    # ---- check 1: scores file -------------------------------------------------------------
    df = read_scores(scores_path, args.label, storms)
    print(f"CHECK 1 ok: scores file has {len(df)} runs, all named for {args.label}.")

    # ---- check 2: gauge ---------------------------------------------------------------------
    gauge = read_gauge(gauge_path)
    windows, real_win = {}, {}
    for tg in storms:
        w0, w1 = pd.Timestamp(STORMS[tg][0]), pd.Timestamp(STORMS[tg][1])
        gw = gauge.loc[w0:w1]
        if len(gw) < 3:
            die(f"The gauge workbook has only {len(gw)} readings in the {STORMS[tg][2]} window.")
        expect = STORMS[tg][3]
        got = float(gw.max())
        if abs(got - expect) / expect > REAL_PEAK_TOL:
            die(f"Gauge peak for {STORMS[tg][2]} is {got:.0f} cfs but the record says {expect:.0f} cfs. "
                f"Wrong workbook, or the window is wrong.")
        windows[tg], real_win[tg] = (w0, w1), gw
        print(f"CHECK 2 ok: {STORMS[tg][2]}: {len(gw)} gauge readings, highest {got:.0f} cfs (record {expect:.0f}).")

    # ---- choose runs ------------------------------------------------------------------------
    ensemble, roles = pick_runs(df, storms)
    to_draw = list(dict.fromkeys(ensemble + list(roles)))
    print(f"\nRuns in the ballpark on BOTH 12 Aug and 8 Sep: {len(ensemble)} of {len(df)}")
    for rid, rl in roles.items():
        row = df[df["run_id"] == rid].iloc[0]
        print(f"  {' + '.join(rl):<34s} Ks {row.Ks_mult:.3f}  f {row.f_RS_abs:.5f}  Kchan {row.channelconductivity_mmhr:.1f}   ({rid})")
    if not ensemble:
        print("  (none: the gray set will be empty)")

    # ---- check 3: outlet files ---------------------------------------------------------------
    if not kept_dir.exists():
        die(f"The folder of kept outlet files was not found:\n  {kept_dir}")
    sims, soft_bad, soft_n = {}, [], 0
    for rid in to_draw:
        s = read_qout(find_qout(kept_dir, rid))
        for tg in storms:
            w0, w1 = windows[tg]
            if s.index.min() > w0 or s.index.max() < w1:
                die(f"Outlet file for {rid} does not cover the {STORMS[tg][2]} window "
                    f"({s.index.min()} to {s.index.max()}).")
        sims[rid] = s
        row = df[df["run_id"] == rid].iloc[0]
        for tg in storms:
            w0, w1 = windows[tg]
            recomputed = float(s.loc[w0:w1].max()) / CFS_TO_CMS
            listed = float(row[f"{tg}_sim_peak_cfs"])
            soft_n += 1
            if listed >= FLOW_MIN_CFS:
                if abs(recomputed - listed) / listed > SIM_PEAK_TOL:
                    soft_bad.append((rid, tg, listed, recomputed))
            elif recomputed >= FLOW_MIN_CFS * 1.5:
                soft_bad.append((rid, tg, listed, recomputed))
    print(f"\nCHECK 3 ok: {len(to_draw)} outlet files found, readable, and cover every storm window.")
    print(f"CHECK (soft): simulated peaks recomputed from the outlet files agree with the scores file "
          f"(within {SIM_PEAK_TOL:.0%}) for {soft_n - len(soft_bad)} of {soft_n} run-storm pairs.")
    for rid, tg, listed, rec in soft_bad[:6]:
        print(f"   differs: {rid} {tg}: scores file {listed:.1f} cfs, outlet file {rec:.1f} cfs")

    if args.dry_run:
        print("\nDry run finished. Nothing was written. Run again without --dry_run to draw the figure.")
        return

    # ---- refuse to overwrite ------------------------------------------------------------------
    out_dir.mkdir(parents=True, exist_ok=True)
    if not args.overwrite:
        for p in (fig_path, csv_path):
            if p.exists():
                die(f"{p.name} already exists. Use --overwrite to replace it, or move it first.")

    # ---- draw ------------------------------------------------------------------------------------
    to_unit_sim = (1.0 / CFS_TO_CMS) if args.units == "cfs" else 1.0
    to_unit_real = 1.0 if args.units == "cfs" else CFS_TO_CMS
    ylab = "discharge (cfs)" if args.units == "cfs" else "discharge (m3/s)"

    fig, axes = plt.subplots(1, len(storms), figsize=(5.6 * len(storms), 5.6), squeeze=False)
    axes = axes[0]
    for ax, tg in zip(axes, storms):
        w0, w1 = windows[tg]
        gw = real_win[tg]
        active = gw[gw >= FLOW_MIN_CFS]
        if len(active):
            x0 = max(w0, active.index.min() - pd.Timedelta(hours=1))
            x1 = min(w1, active.index.max() + pd.Timedelta(hours=2))
        else:
            x0, x1 = w0, w1
        pad = pd.Timedelta(minutes=10)
        ymax = float(gw.loc[x0:x1].max()) * to_unit_real
        for rid in to_draw:
            if rid in roles or rid not in ensemble:
                continue
            seg = sims[rid].loc[x0 - pad:x1 + pad] * to_unit_sim
            ax.plot(seg.index, seg.values, color=GRAY, lw=0.8, alpha=0.30, zorder=2)
            ymax = max(ymax, float(seg.max()))
        for rid, rl in roles.items():
            colr, ls, lw = style_for(rl)
            seg = sims[rid].loc[x0 - pad:x1 + pad] * to_unit_sim
            ax.plot(seg.index, seg.values, color=colr, ls=ls, lw=lw, zorder=4)
            ymax = max(ymax, float(seg.max()))
        gseg = gw.loc[x0 - pad:x1 + pad] * to_unit_real
        ax.plot(gseg.index, gseg.values, color=INK, lw=0.9, alpha=0.7, zorder=5)
        ax.plot(gseg.index, gseg.values, "o", color=INK, ms=3.2, zorder=6)
        ax.set_xlim(x0, x1)
        ax.set_ylim(0, ymax * 1.10 if ymax > 0 else 1)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax.xaxis.set_major_locator(mdates.AutoDateLocator(minticks=4, maxticks=7))
        ax.set_title(STORMS[tg][2], fontsize=12, color=INK)
        ax.set_ylabel(ylab, fontsize=10)
        ax.set_xlabel("time of day (local, as recorded)", fontsize=9)
        ax.grid(alpha=0.25)
        ax.text(0.98, 0.97, f"gauge peak {float(gw.max()) * to_unit_real:,.0f} {args.units}",
                transform=ax.transAxes, ha="right", va="top", fontsize=9, color=INK)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)

    handles = [Line2D([0], [0], color=INK, marker="o", ms=4, lw=0.9, label="gauge readings (as recorded)")]
    if ensemble:
        handles.append(Line2D([0], [0], color=GRAY, lw=1.2, alpha=0.6,
                              label=f"runs within a factor of 2 on 12 Aug and 8 Sep (n = {len(ensemble)})"))
    for rid, rl in roles.items():
        colr, ls, lw = style_for(rl)
        row = df[df["run_id"] == rid].iloc[0]
        handles.append(Line2D([0], [0], color=colr, ls=ls, lw=lw,
                              label=f"{' + '.join(rl)} (Ks {row.Ks_mult:.2f}, f {row.f_RS_abs:.5f}, Kchan {row.channelconductivity_mmhr:.0f})"))
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.045), ncol=2, fontsize=8.5, frameon=False)
    fig.suptitle(f"Series 110, sweep {args.label}: model runs against the real gauge", fontsize=13, color=INK)
    fig.text(0.5, 0.005,
             "Model: tRIBS 6.0.0, Kchan on, routing pinned, 2014. Gauge: event-triggered readings joined by straight lines for the eye "
             "(the scores use a filled 5-minute series). The 0.5-2.0 band is the author's choice. Model output, not a calibration.",
             ha="center", va="bottom", fontsize=7.5, color="#555555")
    fig.tight_layout(rect=(0, 0.20, 1, 0.95))
    fig.savefig(fig_path, dpi=200)
    plt.close(fig)

    # ---- the table of runs that were drawn ----------------------------------------------------
    rows = []
    for rid in to_draw:
        r = df[df["run_id"] == rid].iloc[0]
        rec = {"run_id": rid, "role": " + ".join(roles[rid]) if rid in roles else "factor-of-2 set",
               "Ks_mult": r.Ks_mult, "f_RS_abs": r.f_RS_abs, "Kchan_mmhr": r.channelconductivity_mmhr,
               "combined_ALL": r.combined_ALL}
        for tg in storms:
            for c in ("kge_2012", "peak_ratio", "vol_ratio", "peak_dt_min"):
                rec[f"{tg}_{c}"] = r[f"{tg}_{c}"]
        rows.append(rec)
    pd.DataFrame(rows).to_csv(csv_path, index=False)

    print(f"\nWrote:\n  {fig_path}\n  {csv_path}")
    print("Nothing else was changed.")


if __name__ == "__main__":
    main()
