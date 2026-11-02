#!/usr/bin/env python3
"""
score_mainstorm_lhs_volume_sycamore.py  --  Stage 3: re-score the 150 W1 runs on the MAIN STORM ONLY.

WHY. The flood-week volume ratio adds up the whole week, Jan 18-25. Your hydrograph figure showed the model is too wet during the
main storm (Jan 20-23) and then almost misses the second storm (Jan 23 12:00 onward), so the two errors cancel in the week total.
A run can sit at "volume 1.0" for the week while the main storm is 10-35% too wet.

WHAT THIS DOES. It scores each run again over a shorter window, Jan 18 00:00 to Jan 23 12:00 (the main storm and its recession,
with the second storm left out), using the same scorer functions as score_run_sycamore.py:

    main-storm volume ratio    model volume / gauge volume in the window (1.0 = same amount of water)
    peak error (%)             model peak / gauge peak - 1
    first-arrival error (h)    when the flow first rises through 20% of the gauge peak, model minus gauge (negative = model early)
    KGE_2012, PBIAS            read only AFTER the volume ratio, never alone

It then prints how many runs pass the volume band in the main-storm window versus the whole week, which runs pass one but not the
other, what Ks, f and PsiB look like for the runs that pass, and the Ks that gives main-storm volume 1.0 at the baseline f and PsiB.

OUTPUT (into lhs_volume_plots/; it only READS results/, the gauge file and the CSV, and never changes any model file):
    lhs_volume_mainstorm.csv            one row per run with the main-storm numbers
    mainstorm_vs_week_ratio.png         every run: week ratio across, main-storm ratio up; shaded = in band
    mainstorm_band_hydrographs.png      the runs that pass the main-storm band against the gauge (whole week and close-up of the peak)

Run it from workspaces/Sycamore/lab/ (the folder with score_run_sycamore.py and results/):
    python -u score_mainstorm_lhs_volume_sycamore.py --dry-run     # checks the files and scores ONE run; writes nothing
    python -u score_mainstorm_lhs_volume_sycamore.py               # scores all runs, writes the CSV and two figures

Options: --csv FILE (lhs_volume_summary.csv)  --results-dir DIR (results)  --obs FILE  --out-dir DIR (lhs_volume_plots)
         --start "2017-01-18 00:00"  --end "2017-01-23 12:00"  --band 0.10  --dpi 150
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")                      # no screen needed (Codespace); files only
import matplotlib.dates as mdates          # noqa: E402
import matplotlib.pyplot as plt            # noqa: E402
import numpy as np                         # noqa: E402
import pandas as pd                        # noqa: E402
from matplotlib import cm                  # noqa: E402
from matplotlib.colors import Normalize    # noqa: E402

DEFAULT_OBS = "../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv"
WINDOW = "W1"
WEEK = ("2017-01-18 00:00", "2017-01-25 00:00")
PEAK_ZOOM = ("2017-01-20 18:00", "2017-01-21 12:00")
WEEK_VOL = "flood_week__vol_ratio"
BASE = {"KS_MULT": 0.2, "F_MULT": 1.0, "PSIB_MULT": 0.5}


def log(msg: str = "") -> None:
    print(msg, flush=True)


def die(msg: str) -> None:
    print(f"\nERROR: {msg}\n", file=sys.stderr, flush=True)
    sys.exit(2)


def import_scorer():
    for p in (Path.cwd(), Path(__file__).resolve().parent):
        if str(p) not in sys.path:
            sys.path.insert(0, str(p))
    try:
        import score_run_sycamore as sc
    except ImportError:
        die("cannot import score_run_sycamore.py. Run this from the lab/ folder.")
    for name in ("read_obs", "read_qout", "RUN_WINDOWS", "Flow", "arrival", "kge_parts"):
        if not hasattr(sc, name):
            die(f"score_run_sycamore.py has no '{name}'. Is it the current copy from lab/?")
    return sc


def load_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        die(f"summary CSV not found: {path}. Run this from lab/, or give --csv.")
    df = pd.read_csv(path)
    need = ["run_id", "KS_MULT", "F_MULT", "PSIB_MULT", WEEK_VOL]
    miss = [c for c in need if c not in df.columns]
    if miss:
        die(f"{path} is missing columns {miss}")
    if "status" in df.columns:
        df = df[df["status"] == "done"]
    df = df[np.isfinite(df[WEEK_VOL])].copy()
    df["draw"] = df["run_id"].str.extract(r"_d(\d+)_", expand=False).astype(float).astype("Int64")
    return df.reset_index(drop=True)


def read_run(sc, results_dir: Path, rid: str):
    f = results_dir / f"{rid}_Outlet.qout"
    if not f.exists():
        hits = sorted(results_dir.glob(f"**/{rid}_Outlet.qout"))
        f = hits[0] if len(hits) == 1 else None
    if f is None:
        return None
    q = sc.read_qout(f, pd.Timestamp(sc.RUN_WINDOWS[WINDOW][0]))
    return None if q.empty else q


def score_window(sc, model_q, G, a, b):
    """Main-storm numbers for one run. G is the scorer's Flow of the whole gauge record."""
    M = sc.Flow(model_q)
    mw, gw = M.window(a, b), G.window(a, b)
    vm = M.volume(a, b, match_to=G) / 1e6          # only where the gauge has data, as the scorer does
    vg = G.volume(a, b) / 1e6
    pm, pg = float(mw.max()), float(gw.max())
    tg, fg = sc.arrival(gw, 0.2 * pg)
    tm, fm = sc.arrival(mw, 0.2 * pg)
    arr = (tm - tg) / pd.Timedelta(hours=1) if (fg == "ok" and fm == "ok") else np.nan
    kp = sc.kge_parts(mw, gw)
    return dict(main_vol_model_Mm3=vm, main_vol_gauge_Mm3=vg, main_vol_ratio=vm / vg if vg > 0 else np.nan,
                main_peak_err_pct=(pm / pg - 1) * 100, main_arrival_err_h=arr, main_kge2012=kp["kge"], main_pbias_pct=kp["pbias"])


def clip(s, rng):
    return s[(s.index >= pd.Timestamp(rng[0])) & (s.index <= pd.Timestamp(rng[1]))]


# ------------------------------------------------------------------------------------------------------------------
def fit_ks_for_one(t: pd.DataFrame):
    """ln(main ratio) = a + b ln Ks + c ln f + d ln PsiB (least squares); returns (R2, Ks that gives 1.0 at baseline f and PsiB, coefficients)."""
    X = np.column_stack([np.ones(len(t)), np.log(t["KS_MULT"]), np.log(t["F_MULT"]), np.log(t["PSIB_MULT"])])
    y = np.log(t["main_vol_ratio"].to_numpy(float))
    ok = np.isfinite(y)
    beta, *_ = np.linalg.lstsq(X[ok], y[ok], rcond=None)
    pred = X[ok] @ beta
    r2 = 1 - ((y[ok] - pred) ** 2).sum() / ((y[ok] - y[ok].mean()) ** 2).sum()
    a, b, c, d = beta
    ks1 = float(np.exp(-(a + c * np.log(BASE["F_MULT"]) + d * np.log(BASE["PSIB_MULT"])) / b)) if b != 0 else np.nan
    return float(r2), ks1, beta


def fig_scatter(t, lo, hi, out, dpi):
    wk, mn = t[WEEK_VOL], t["main_vol_ratio"]
    fig, ax = plt.subplots(figsize=(7.4, 6.6))
    ax.axvspan(lo, hi, color="#4477AA", alpha=0.12, lw=0)
    ax.axhspan(lo, hi, color="#CC6677", alpha=0.12, lw=0)
    both = (wk.between(lo, hi)) & (mn.between(lo, hi))
    only_main = mn.between(lo, hi) & ~wk.between(lo, hi)
    only_week = wk.between(lo, hi) & ~mn.between(lo, hi)
    neither = ~(both | only_main | only_week)
    for mask, col, lab in ((neither, "#BBBBBB", "neither band"), (only_week, "#4477AA", "week band only"),
                           (only_main, "#CC6677", "main-storm band only"), (both, "#222222", "both bands")):
        ax.scatter(wk[mask], mn[mask], s=26, color=col, label=f"{lab} ({int(mask.sum())})", edgecolor="white", lw=0.4, zorder=3)
    ax.axhline(1, color="gray", lw=0.8)
    ax.axvline(1, color="gray", lw=0.8)
    ax.set_xlabel("Whole-week volume ratio (model / gauge)   [blue shading = in the week band]", fontsize=9)
    ax.set_ylabel("Main-storm volume ratio (model / gauge)   [red shading = in the main-storm band]", fontsize=9)
    ax.set_title("Same runs, two volume tests", fontsize=11)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, loc="upper left", frameon=True)
    fig.tight_layout()
    fig.savefig(out, dpi=dpi)
    plt.close(fig)


def fig_hydro(sc, t, flows, gauge, lo, hi, a, b, out, dpi):
    sub = t[t["main_vol_ratio"].between(lo, hi)].sort_values("main_vol_ratio")
    norm, cmap = Normalize(vmin=lo, vmax=hi), plt.get_cmap("viridis")
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2))
    for ax, rng in zip(axes, (WEEK, PEAK_ZOOM)):
        for rid, v in zip(sub["run_id"], sub["main_vol_ratio"]):
            q = clip(flows[rid], rng)
            ax.plot(q.index, q.values, color=cmap(norm(v)), lw=0.9, alpha=0.7)
        g = clip(gauge.dropna(), rng)
        ax.plot(g.index, g.values, color="black", lw=2.0, zorder=10)
        span_h = (pd.Timestamp(rng[1]) - pd.Timestamp(rng[0])) / pd.Timedelta(hours=1)
        if span_h <= 48:
            ax.xaxis.set_major_locator(mdates.HourLocator(interval=4))
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d\n%H:%M"))
        else:
            ax.xaxis.set_major_locator(mdates.DayLocator())
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
            ax.axvspan(pd.Timestamp(b), pd.Timestamp(rng[1]), color="gray", alpha=0.12, lw=0)
            ax.text(pd.Timestamp(b) + pd.Timedelta(hours=2), ax.get_ylim()[1] * 0.45, "not scored here\n(second storm)", fontsize=8, color="dimgray")
        ax.set_xlim(pd.Timestamp(rng[0]), pd.Timestamp(rng[1]))
        ax.set_ylabel("Outlet streamflow (m³/s)", fontsize=9)
        ax.grid(alpha=0.25)
        ax.tick_params(labelsize=8)
    axes[0].set_title("Whole week", fontsize=10, loc="left")
    axes[1].set_title("Close-up of the main peak (Jan 20 18:00 to Jan 21 12:00)", fontsize=10, loc="left")
    axes[0].legend(handles=[plt.Line2D([], [], color="black", lw=2, label="USGS gauge")], loc="upper left", fontsize=9, frameon=False)
    fig.suptitle(f"{len(sub)} runs inside the MAIN-STORM volume band ({lo:.1f} to {hi:.1f}, scored {pd.Timestamp(a):%b %d} to {pd.Timestamp(b):%b %d %H:%M}); color = main-storm volume ratio",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 0.90, 0.95))
    cax = fig.add_axes([0.925, 0.2, 0.012, 0.6])
    cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cb.set_label("main-storm volume ratio", fontsize=9)
    cb.ax.tick_params(labelsize=8)
    fig.savefig(out, dpi=dpi)
    plt.close(fig)


# ------------------------------------------------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description="Re-score the Stage 3 W1 runs on the main storm only (see the top of this file).")
    ap.add_argument("--csv", default="lhs_volume_summary.csv")
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--obs", default=DEFAULT_OBS)
    ap.add_argument("--out-dir", default="lhs_volume_plots")
    ap.add_argument("--start", default="2017-01-18 00:00")
    ap.add_argument("--end", default="2017-01-23 12:00")
    ap.add_argument("--band", type=float, default=0.10)
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    lo, hi = 1 - a.band, 1 + a.band

    sc = import_scorer()
    df = load_table(Path(a.csv))
    gauge = sc.read_obs(a.obs)
    G = sc.Flow(gauge)
    results_dir = Path(a.results_dir)
    log(f"score_mainstorm_lhs_volume_sycamore.py   {len(df)} finished runs in {a.csv}")
    log(f"  window scored: {a.start} to {a.end}   (volume band {lo:.2f} to {hi:.2f})")
    gw = G.window(a.start, a.end)
    if gw.dropna().empty:
        die("the gauge has no data in that window")
    log(f"  gauge in the window: peak {gw.max():.1f} m3/s, volume {G.volume(a.start, a.end) / 1e6:.3f} M m3")

    if a.dry_run:
        r = df.iloc[0]
        q = read_run(sc, results_dir, r["run_id"])
        if q is None:
            die(f"cannot read results/{r['run_id']}_Outlet.qout (is --results-dir right?)")
        m = score_window(sc, q, G, a.start, a.end)
        log(f"  test run {r['run_id']}:")
        log(f"    week volume ratio (from the CSV) {r[WEEK_VOL]:.3f};  main-storm volume ratio {m['main_vol_ratio']:.3f};  "
            f"peak error {m['main_peak_err_pct']:+.1f}%;  arrival {m['main_arrival_err_h']:+.2f} h;  KGE {m['main_kge2012']:.3f}")
        log("\nDRY RUN: files found and one run scored. Nothing was written.")
        return 0

    rows, flows, missing = [], {}, []
    for i, r in df.iterrows():
        q = read_run(sc, results_dir, r["run_id"])
        if q is None:
            missing.append(r["run_id"])
            continue
        flows[r["run_id"]] = q
        rows.append({**r[["run_id", "draw", "KS_MULT", "F_MULT", "PSIB_MULT", WEEK_VOL]].to_dict(), **score_window(sc, q, G, a.start, a.end)})
        if (i + 1) % 50 == 0:
            log(f"  scored {i + 1} of {len(df)}")
    if missing:
        log(f"  WARNING: {len(missing)} run(s) had no readable outlet file and are left out, for example {missing[:3]}")
    if not rows:
        die(f"no outlet files could be read from {results_dir}")
    t = pd.DataFrame(rows)
    t["in_week_band"] = t[WEEK_VOL].between(lo, hi)
    t["in_main_band"] = t["main_vol_ratio"].between(lo, hi)

    wk, mn = t["in_week_band"], t["in_main_band"]
    log(f"\nHOW MANY RUNS PASS THE VOLUME BAND ({lo:.1f} to {hi:.1f})")
    log(f"  whole week:               {int(wk.sum())} of {len(t)}")
    log(f"  main storm only:          {int(mn.sum())} of {len(t)}")
    log(f"  both:                     {int((wk & mn).sum())}")
    log(f"  week band only:           {int((wk & ~mn).sum())}   (right total, but the main storm is off)")
    log(f"  main-storm band only:     {int((mn & ~wk).sum())}   (main storm right, week total off because the second storm is missing)")
    log(f"  main-storm ratio of the {int(wk.sum())} week-band runs: {t.loc[wk, 'main_vol_ratio'].min():.2f} to {t.loc[wk, 'main_vol_ratio'].max():.2f}")
    if mn.any():
        s = t[mn]
        log(f"\nTHE {len(s)} RUNS INSIDE THE MAIN-STORM BAND")
        log(f"  {'':<22}{'min':>8}{'median':>8}{'max':>8}")
        for lab, col, nd in (("Ks multiplier", "KS_MULT", 3), ("f multiplier", "F_MULT", 2), ("PsiB multiplier", "PSIB_MULT", 2),
                             ("week volume ratio", WEEK_VOL, 2), ("peak error (%)", "main_peak_err_pct", 1),
                             ("arrival error (h)", "main_arrival_err_h", 2), ("KGE_2012", "main_kge2012", 3)):
            log(f"  {lab:<22}{s[col].min():>8.{nd}f}{s[col].median():>8.{nd}f}{s[col].max():>8.{nd}f}")
        top = s.sort_values("main_kge2012", ascending=False).head(5)
        log("  five best KGE inside the band (a tie-break only, read after the volume test):")
        for _, r in top.iterrows():
            log(f"    draw {int(r['draw']):>3}  Ks {r['KS_MULT']:.3f}  f {r['F_MULT']:.2f}  PsiB {r['PSIB_MULT']:.2f}   volume {r['main_vol_ratio']:.2f}  "
                f"peak {r['main_peak_err_pct']:+.0f}%  arrival {r['main_arrival_err_h']:+.1f} h  KGE {r['main_kge2012']:.3f}")
    r2, ks1, beta = fit_ks_for_one(t)
    log(f"\nFIT  ln(main-storm ratio) = a + b ln Ks + c ln f + d ln PsiB:   R2 {r2:.3f};   b {beta[1]:+.2f} (Ks), c {beta[2]:+.2f} (f), d {beta[3]:+.2f} (PsiB)")
    log(f"  Ks that gives main-storm volume 1.0 at the baseline f = 1 and PsiB = 0.5:  {ks1:.3f}   (baseline Ks is 0.2)")
    log("  A Ks outside the sampled range 0.15 to 0.30 is an extrapolation: treat it as a hint, not a result.")

    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    t.to_csv(out / "lhs_volume_mainstorm.csv", index=False)
    fig_scatter(t, lo, hi, out / "mainstorm_vs_week_ratio.png", a.dpi)
    if mn.any():
        fig_hydro(sc, t, flows, gauge, lo, hi, a.start, a.end, out / "mainstorm_band_hydrographs.png", a.dpi)
    log(f"\nWROTE {out / 'lhs_volume_mainstorm.csv'}, {out / 'mainstorm_vs_week_ratio.png'}" + (f", {out / 'mainstorm_band_hydrographs.png'}" if mn.any() else "")
        + ".\nNothing in results/ was touched.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
