"""
run_w5_sycamore.py  --  ONE whole-water-year run (W5: Oct 1 2016 - Sep 30 2017 23:00, 8,759 h) of the WORKING setting, channel loss off.
The first run of Apr - Sep (the monsoon) in this calibration, and the first real timing of a W5 run (the 68 minutes in the plan is scaled, never timed).

THE RUN (the same setting as the finished W4 working run W4_d012_k10_f0p15_as300_ksx1p03_nocl, with a longer window):
    Ks x1.03 (KS_MULT 0.19261), f x0.15, As 300, Au 1, PsiB x0.2939, routing k10 (cv 6.39, flowexp 0.2973, n 0.067), width coefficient 2.33,
    bedrock +1 m, water table 0.95 x bedrock at Oct 1, snow ON, channel loss OFF.     Run name: W5_d012_k10_f0p15_as300_ksx1p03_nocl
This script is built on the tested W4 scripts (run_tail_screen_sycamore.py and run_retune_w4_sycamore.py, both unchanged and both must be in lab/): it only
changes the window to the whole year, so the notebook build, the checks, the tRIBS run and the scoring are the ones already used on W4.

WHY: nobody has looked at Apr - Sep. The plan notes that July radar rain averaged 147 mm while the gauge never went above 0.45 m3/s, so false flow in the
monsoon could be larger than the Oct - Dec false flow (about 3.3 M m3). This run answers that, and gives the whole-year score and the first real W5 time.

ONE THING TO KNOW. The weather files end at 09/30/2017 23:00, so a full 8,760 h run would end one hour past them (the notebook warns about this). To avoid losing a
run to that, this script runs 8,759 h BY DEFAULT: Oct 1 2016 00:00 to Sep 30 2017 23:00. The scorer accepts a run that ends within 1.5 h of the window end, so every
whole-year metric is still computed (the last hour of September is simply not in the model). --hours 8760 is there only for a test of the full length.

WHAT IT CHECKS BEFORE IT STARTS (the W4 script's checks, plus one): the W4 working run is finished and its .in file is there (the new .in must equal it in every
keyword except the window and the output name), the radar folder has all the hourly files the whole year needs, free disk space, no other tRIBS running.

WHAT IT PRINTS AFTER THE RUN (also with --summary-only):
    SCORE        the scorer's headline numbers for the W5 run next to the W4 working run (the first 4,368 h should be the same run).
    CONSISTENCY  the W5 outlet flow over Oct 1 - Apr 1 against the W4 working run's (hourly means): the largest difference and the volume ratio.
    MONTHS       for each month: model and gauge volume (M m3), model and gauge maximum flow, days above 0.13 m3/s, hours at zero flow.
    SEASONS      Oct-Dec, Jan-Mar, Apr-Jun, Jul-Sep: model and gauge volume. The plan's reading rule is printed as a comparison of Jul-Sep with Oct-Dec; it decides nothing.
    BIGGEST DAYS the ten days of Apr - Sep with the most model flow, with the gauge's flow on the same days.
    PROCESS      build and tRIBS minutes (the first real W5 time), output size, disk space left.
It writes w5_months.csv, and with --plot w5_figs/fig_w5.png. It changes no setting, makes no recommendation, and does not touch the design doc.

HOW TO RUN. Put it in workspaces/Sycamore/lab/ next to run_tail_screen_sycamore.py AND run_retune_w4_sycamore.py. Foreground, one command at a time:
    python -u run_w5_sycamore.py --dry-run                    # checks only, builds nothing
    python -u run_w5_sycamore.py --build-only                 # builds the inputs for the whole year and checks the .in file, does NOT run tRIBS (a few minutes)
    python -u run_w5_sycamore.py                              # the run (expect roughly 60 - 80 minutes of tRIBS); safe to run again if the connection drops
    python -u run_w5_sycamore.py --summary-only --plot        # tables and figure from the finished run
Do NOT run the notebook (or anything that rebuilds data/model/) while it runs. The radar-vs-gauge check (read-only) may run in a second terminal.
"""
from __future__ import annotations

import argparse
import math
import os
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_retune_w4_sycamore as R          # noqa: E402  (importing it builds nothing)
T = R.T                                      # the tested W4 script that builds, runs and scores

# --------------------------------------------------------------------------------------------
# Names and the window
# --------------------------------------------------------------------------------------------
W5 = dict(name="W5", start="10/01/2016/00/00", hours=8759, label="full_year", gw=0.95)
WORK_KS = 1.03
SUMMARY_NAME = "w5_months.csv"
FIG_DIR = "w5_figs"
THR_DAYS = 0.13                      # m3/s: the gauge's own Oct-Dec maximum; a day counts as "above" when its highest hourly flow is over this
ZERO_THR = 0.001                     # m3/s: an hour whose mean outlet flow is under this counts as zero-flow
RADAR_DIR_DEFAULT = "~/mrms_wk_wy2017/rain"
Y0, Y_W4_END, Y1 = "2016-10-01 00:00", "2017-04-01 00:00", "2017-10-01 00:00"
SEASONS = [("Oct-Dec", "2016-10-01 00:00", "2017-01-01 00:00"), ("Jan-Mar", "2017-01-01 00:00", "2017-04-01 00:00"),
           ("Apr-Jun", "2017-04-01 00:00", "2017-07-01 00:00"), ("Jul-Sep", "2017-07-01 00:00", "2017-10-01 00:00")]
NAN = float("nan")


def w4_working_job() -> dict:
    """The W4 working job, exactly as the re-tune script made it (built BEFORE the window is changed, so its name is the W4 name)."""
    base = T.an_job(0.15, 300)
    return R._with_tag(base, f"ksx{T.vtag(WORK_KS)}", {"KS_MULT": round(T.BASE["KS_MULT"] * WORK_KS, 6)}, "ks")


W4_JOB = w4_working_job()
W4_ID = W4_JOB["run_id"]
W5_ID = W4_ID.replace("W4_", "W5_", 1)


def w5_job() -> dict:
    j = dict(W4_JOB)
    j["run_id"] = W5_ID
    j["step"] = "ks"
    return j


def describe(job: dict) -> str:
    return "working setting (Ks x1.03, f x0.15, As 300, Au 1), channel loss off, whole water year"


def configure(hours: int = 8759) -> None:
    T.WIN.clear()
    T.WIN.update(W5)
    T.WIN["hours"] = int(hours)
    T.BASELINE_ID = W4_ID                 # the W5 .in is checked against the W4 working run's .in, and that run must be finished
    T.SUMMARY_CSV = "w5_summary.csv"
    T.describe = describe


# --------------------------------------------------------------------------------------------
# Numbers from the outlet file and the gauge
# --------------------------------------------------------------------------------------------
def hourly(q: pd.Series, a, b) -> pd.Series:
    s = q[(q.index >= pd.Timestamp(a)) & (q.index < pd.Timestamp(b))]
    return s.resample("60min").mean().dropna()


def vol_mm3(h: pd.Series) -> float:
    return float(h.sum() * 3600.0 / 1e6) if len(h) else NAN


def days_above(h: pd.Series, thr: float) -> int:
    return int((h.resample("1D").max() > thr).sum()) if len(h) else 0


def zero_hours_pct(h: pd.Series) -> float:
    return float(100.0 * (h < ZERO_THR).mean()) if len(h) else NAN


def month_rows(model_q: pd.Series, gauge_q: pd.Series) -> list:
    rows = []
    for p in pd.period_range("2016-10", "2017-09", freq="M"):
        a, b = p.start_time, (p + 1).start_time
        mh, gh = hourly(model_q, a, b), hourly(gauge_q, a, b)
        rows.append(dict(month=p.strftime("%b %Y"), start=a, end=b, model_Mm3=vol_mm3(mh), gauge_Mm3=vol_mm3(gh),
                         model_max=float(mh.max()) if len(mh) else NAN, gauge_max=float(gh.max()) if len(gh) else NAN,
                         model_days_above=days_above(mh, THR_DAYS), gauge_days_above=days_above(gh, THR_DAYS),
                         model_zero_pct=zero_hours_pct(mh), gauge_zero_pct=zero_hours_pct(gh), hours_model=len(mh), hours_gauge=len(gh)))
    return rows


def season_rows(model_q: pd.Series, gauge_q: pd.Series) -> list:
    rows = []
    for name, a, b in SEASONS:
        mh, gh = hourly(model_q, a, b), hourly(gauge_q, a, b)
        rows.append(dict(season=name, model_Mm3=vol_mm3(mh), gauge_Mm3=vol_mm3(gh), model_days_above=days_above(mh, THR_DAYS), gauge_days_above=days_above(gh, THR_DAYS),
                         model_max=float(mh.max()) if len(mh) else NAN, gauge_max=float(gh.max()) if len(gh) else NAN))
    return rows


def biggest_days(model_q: pd.Series, gauge_q: pd.Series, a, b, n: int = 10) -> list:
    mh, gh = hourly(model_q, a, b), hourly(gauge_q, a, b)
    md = mh.resample("1D").sum() * 3600.0 / 1e6
    out = []
    for day in md.sort_values(ascending=False).head(n).index:
        mg = mh[(mh.index >= day) & (mh.index < day + pd.Timedelta(days=1))]
        gg = gh[(gh.index >= day) & (gh.index < day + pd.Timedelta(days=1))]
        out.append(dict(day=day, model_Mm3=float(md[day]), model_max=float(mg.max()), gauge_Mm3=vol_mm3(gg), gauge_max=float(gg.max()) if len(gg) else NAN))
    return sorted(out, key=lambda r: r["day"])


def consistency(w5_q: pd.Series, w4_q: pd.Series) -> dict:
    a, b = pd.Timestamp(Y0), pd.Timestamp(Y_W4_END)
    h5, h4 = hourly(w5_q, a, b), hourly(w4_q, a, b)
    both = pd.concat([h5.rename("w5"), h4.rename("w4")], axis=1).dropna()
    if both.empty:
        return {}
    d = (both["w5"] - both["w4"]).abs()
    return dict(hours=len(both), max_abs_diff=float(d.max()), t_max=d.idxmax(), vol_w5=vol_mm3(both["w5"]), vol_w4=vol_mm3(both["w4"]),
                ratio=float(both["w5"].sum() / both["w4"].sum()) if both["w4"].sum() > 0 else NAN)


# --------------------------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------------------------
def fm(x, nd=2, w=8) -> str:
    return f"{'-':>{w}}" if (x is None or not math.isfinite(x)) else f"{x:{w}.{nd}f}"


HEADLINES = [("flood_week__vol_ratio", "flood-week volume ratio", 3), ("flood_week__peak_err_pct", "flood-week peak error %", 1),
             ("feb27_28__vol_ratio", "Feb 27-28 volume ratio", 3), ("feb27_28__peak_err_pct", "Feb 27-28 peak error %", 1),
             ("feb18_19__vol_ratio", "Feb 18-19 volume ratio", 2), ("tail_jan26_feb2__tail_vol_ratio", "tail volume Jan 26 - Feb 2", 2),
             ("octdec__vol_model_Mm3", "Oct-Dec model volume, M m3", 2), ("julsep__vol_model_Mm3", "Jul-Sep model volume, M m3", 2),
             ("oct_mar__vol_ratio", "Oct-Mar volume ratio", 3), ("span__vol_ratio", "whole-span volume ratio", 3),
             ("span__kge2012", "KGE, whole span (gate first)", 3), ("proc__channel_loss_Mm3", "channel loss, M m3", 3)]


def print_report(man: dict, lab: Path, args) -> None:
    rec = man["runs"].get(W5_ID) or {}
    if rec.get("status") not in ("done", "unscored"):
        T.log(f"{W5_ID}: not finished yet (status: {rec.get('status', 'not started')}), so there is nothing to report.")
        return
    qf = lab / "results" / f"{W5_ID}_Outlet.qout"
    if not qf.is_file():
        T.log(f"{W5_ID}: the outlet file {qf.name} is missing.")
        return
    try:
        gauge_q = T.read_obs(lab / T.INIT_DATA / T.OBS_FILE)
    except Exception as exc:
        T.log(f"  (gauge file could not be read, so the tables are skipped: {exc})")
        return
    w5_q = T.read_qout(qf, T.RUN_START_TS)
    csv_scores = T.load_scores(lab)
    s5, s4 = R.scores_row(man, lab, W5_ID, csv_scores), R.scores_row(man, lab, W4_ID, csv_scores)

    T.log()
    T.log(f"W5 WHOLE-YEAR RUN  {W5_ID}   (Oct 1 2016 - Sep 30 2017 23:00, {T.WIN['hours']} h, working setting, channel loss off)")
    T.log("  Volumes in MONTHS, SEASONS and BIGGEST DAYS are sums of hourly mean flows, M m3 (the scorer integrates the raw steps, so its numbers differ slightly).")
    T.log("  'days above' = days whose highest hourly flow is over 0.13 m3/s (the gauge's Oct-Dec maximum). The gauge has no data on Jan 1-2 (outage): its January volume lacks those two days, the model's does not.")
    T.log()
    T.log("  SCORE   (the first 4,368 h are the W4 working run's window, so the W4 metrics should agree)")
    T.log(f"    {'metric':<34}{'W4 working':>12}{'W5':>10}")
    for key, name, nd in HEADLINES:
        a, b = T.fnum(s4, key), T.fnum(s5, key)
        T.log(f"    {name:<34}{fm(a, nd, 12)}{fm(b, nd, 10)}")
    T.log()

    w4f = lab / "results" / f"{W4_ID}_Outlet.qout"
    T.log("  CONSISTENCY   (W5 outlet over Oct 1 - Apr 1 against the W4 working run, hourly means)")
    if w4f.is_file():
        c = consistency(w5_q, T.read_qout(w4f, T.RUN_START_TS))
        if c:
            T.log(f"    {c['hours']} hours compared; largest difference {c['max_abs_diff']:.4f} m3/s (at {c['t_max']:%b %d %H:%M}); volume W5 {c['vol_w5']:.3f} vs W4 {c['vol_w4']:.3f} M m3 (ratio {c['ratio']:.4f})")
        else:
            T.log("    no hours in common")
    else:
        T.log(f"    the W4 working run's outlet file {w4f.name} is not in results/, so no comparison.")
    T.log()

    months = month_rows(w5_q, gauge_q)
    T.log("  MONTHS")
    T.log(f"    {'month':<9}{'model Mm3':>10}{'gauge Mm3':>10}{'model max':>10}{'gauge max':>10}{'days>0.13 mod':>14}{'gauge':>6}{'zero-h% mod':>12}{'gauge':>7}")
    for r in months:
        T.log(f"    {r['month']:<9}{fm(r['model_Mm3'], 3, 10)}{fm(r['gauge_Mm3'], 3, 10)}{fm(r['model_max'], 2, 10)}{fm(r['gauge_max'], 2, 10)}"
              f"{r['model_days_above']:>14d}{r['gauge_days_above']:>6d}{fm(r['model_zero_pct'], 1, 12)}{fm(r['gauge_zero_pct'], 1, 7)}")
    tm, tg = sum(r["model_Mm3"] for r in months), sum(r["gauge_Mm3"] for r in months)
    T.log(f"    whole year: model {tm:.2f} M m3, gauge {tg:.2f} M m3 (ratio {tm / tg:.3f})" if tg > 0 else f"    whole year: model {tm:.2f} M m3")
    T.log()

    seasons = season_rows(w5_q, gauge_q)
    T.log("  SEASONS")
    T.log(f"    {'season':<9}{'model Mm3':>10}{'gauge Mm3':>10}{'model max':>10}{'gauge max':>10}{'days>0.13 mod':>14}{'gauge':>6}")
    for r in seasons:
        T.log(f"    {r['season']:<9}{fm(r['model_Mm3'], 3, 10)}{fm(r['gauge_Mm3'], 3, 10)}{fm(r['model_max'], 2, 10)}{fm(r['gauge_max'], 2, 10)}{r['model_days_above']:>14d}{r['gauge_days_above']:>6d}")
    od = next(r for r in seasons if r["season"] == "Oct-Dec")
    js = next(r for r in seasons if r["season"] == "Jul-Sep")
    word = "LARGER than" if js["model_Mm3"] > od["model_Mm3"] else "smaller than (or equal to)"
    T.log(f"    Plan's reading rule: Jul-Sep model volume {js['model_Mm3']:.2f} M m3 is {word} Oct-Dec {od['model_Mm3']:.2f} M m3 "
          f"(gauge: {js['gauge_Mm3']:.2f} and {od['gauge_Mm3']:.2f}). If larger, ask what drives it before Stage 7. Nothing is decided here.")
    T.log()

    T.log("  BIGGEST DAYS  (Apr 1 - Oct 1: the ten days with the most model volume; the gauge's flow on the same day)")
    T.log(f"    {'day':<8}{'model Mm3':>10}{'model max':>10}{'gauge Mm3':>10}{'gauge max':>10}")
    for r in biggest_days(w5_q, gauge_q, "2017-04-01 00:00", Y1):
        T.log(f"    {r['day']:%b %d}  {fm(r['model_Mm3'], 3, 10)}{fm(r['model_max'], 2, 10)}{fm(r['gauge_Mm3'], 3, 10)}{fm(r['gauge_max'], 2, 10)}")
    T.log()

    hrs = T.WIN["hours"]
    tm_ = rec.get("tribs_min")
    T.log("  PROCESS")
    T.log(f"    build {rec.get('build_min', '?')} min; tRIBS {tm_ if tm_ is not None else '?'} min"
          + (f" = {tm_ / hrs * 1000:.2f} min per 1,000 model hours (the plan used 7.3 - 9.5; scaled W5 estimate about 68 min)" if tm_ else "")
          + f"; output {rec.get('output_mb', '?')} MB; {T.disk_free_gb(lab):.1f} GB free now")
    w4rec = man["runs"].get(W4_ID) or {}
    if w4rec.get("tribs_min"):
        T.log(f"    (the W4 working run took {w4rec['tribs_min']} min for 4,368 h = {w4rec['tribs_min'] / 4368 * 1000:.2f} min per 1,000 h; output {w4rec.get('output_mb', '?')} MB)")
    pd.DataFrame(months).drop(columns=["start", "end"]).to_csv(lab / SUMMARY_NAME, index=False)
    T.log(f"  Wrote {SUMMARY_NAME}. These tables decide nothing.")


def make_plot(man: dict, lab: Path, out: Path) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    w5_q = T.read_qout(lab / "results" / f"{W5_ID}_Outlet.qout", T.RUN_START_TS)
    gauge_q = T.read_obs(lab / T.INIT_DATA / T.OBS_FILE)
    mh, gh = hourly(w5_q, Y0, Y1), hourly(gauge_q, Y0, Y1)
    months = month_rows(w5_q, gauge_q)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 7), gridspec_kw={"height_ratios": [3, 2]})
    ax1.plot(gh.index, gh.clip(lower=1e-3), color="#8a8984", lw=1.0, label="USGS gauge")
    ax1.plot(mh.index, mh.clip(lower=1e-3), color="#2a78d6", lw=0.9, label="model (W5, working setting)")
    ax1.set_yscale("log")
    ax1.set_ylim(1e-3, 600)
    ax1.set_ylabel("Outlet flow, hourly mean (m$^3$/s), log scale")
    ax1.set_title("Whole water year 2017 (Oct 1 - Oct 1)", loc="left")
    ax1.legend(frameon=False, loc="upper right")
    x = np.arange(len(months))
    ax2.bar(x - 0.2, [r["gauge_Mm3"] for r in months], 0.4, color="#8a8984", label="gauge")
    ax2.bar(x + 0.2, [r["model_Mm3"] for r in months], 0.4, color="#2a78d6", label="model")
    ax2.set_xticks(x)
    ax2.set_xticklabels([r["month"].split()[0] for r in months])
    ax2.set_ylabel("Monthly volume (M m$^3$)")
    ax2.legend(frameon=False)
    for ax in (ax1, ax2):
        ax.grid(alpha=0.25)
        for s_ in ("top", "right"):
            ax.spines[s_].set_visible(False)
    out.parent.mkdir(exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return f"wrote {out.name}"


# --------------------------------------------------------------------------------------------
# Checks, arguments, main
# --------------------------------------------------------------------------------------------
def radar_check(radar_dir: str, hours: int) -> tuple:
    d = Path(radar_dir).expanduser()
    if not d.is_dir():
        return ("FAIL", f"radar folder not found: {d}")
    t0 = pd.Timestamp(Y0)
    need = [f"radar_{(t0 + pd.Timedelta(hours=h)):%m%d%Y%H}.asc" for h in range(1, hours + 1)]
    have = set(os.listdir(d))
    miss = [n for n in need if n not in have]
    if miss:
        return ("FAIL", f"radar folder {d} lacks {len(miss)} of the {len(need)} hourly files the whole year needs (first: {miss[0]})")
    return ("OK", f"radar folder {d} has all {len(need)} hourly files the whole year needs ({need[0]} to {need[-1]})")


def suggested_min_free_gb(man: dict) -> tuple:
    mb = (man["runs"].get(W4_ID) or {}).get("output_mb")
    if mb:
        return max(4.0, round(2.2 * float(mb) / 1000.0 + 1.0, 1)), f"the W4 working run wrote {float(mb):.0f} MB; a whole-year run writes about twice that, plus a margin"
    return 8.0, "the size of the W4 working run's output is not in the manifest, so a cautious default"


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description="One whole-water-year (W5) run of the working setting (see the top of this file).")
    ap.add_argument("--dry-run", action="store_true", help="check everything, list the plan, build nothing")
    ap.add_argument("--summary-only", action="store_true", help="print the tables from the finished run, then stop")
    ap.add_argument("--build-only", action="store_true", help="build and check the inputs but do not run tRIBS")
    ap.add_argument("--plot", action="store_true", help=f"make {FIG_DIR}/fig_w5.png from the finished run")
    ap.add_argument("--hours", type=int, default=8759, help="run length in hours (default 8759: the weather files end 1 h before a full 8760 h year would; 8760 is allowed)")
    ap.add_argument("--radar-dir", default=RADAR_DIR_DEFAULT, help="folder of the hourly radar files, for the check that all of them are there")
    ap.add_argument("--redo", default=None, help="comma-separated run ids to run again even if done")
    ap.add_argument("--timeout-min", type=float, default=None, help="stop tRIBS after this many minutes (default: 3.5 times the expected time, about 3.7 hours)")
    ap.add_argument("--cell-timeout", type=int, default=1800, help="seconds allowed for one notebook cell")
    ap.add_argument("--min-free-gb", type=float, default=None, help="stop if less than this much disk is free (default: worked out from the W4 run's output size)")
    ap.add_argument("--heartbeat-min", type=float, default=T.HEARTBEAT_MIN, help="print a progress line every this many minutes while tRIBS runs (default 5)")
    ap.add_argument("--band-lo", type=float, default=T.BAND_LO)
    ap.add_argument("--band-hi", type=float, default=T.BAND_HI)
    ap.add_argument("--notebook", default=None, help="notebook to run (default: Make_Sycamore_Model.ipynb, else Make_Sycamore_Model_clean.ipynb)")
    ap.add_argument("--scorer", default=T.SCORER)
    ap.add_argument("--tribs", default=os.environ.get("TRIBS_BIN") or (T.TRIBS_DEFAULT if os.path.isfile(T.TRIBS_DEFAULT) else shutil.which("tRIBS") or T.TRIBS_DEFAULT))
    ap.add_argument("--skip-kernel-check", action="store_true", help="skip the pytRIBS import test in --dry-run")
    a = ap.parse_args(argv)
    a.steps = "w5"
    a.limit = None
    return a


def main(argv=None) -> int:
    lost = R.missing_from_w4_script()
    if lost:
        print("run_tail_screen_sycamore.py in this folder is an older copy: it lacks " + ", ".join(lost) + ". Put the newest run_tail_screen_sycamore.py in lab/ and run again.")
        return 2
    args = parse_args(argv)
    if not (24 <= args.hours <= 8760):
        print(f"--hours {args.hours}: must be between 24 and 8760")
        return 2
    configure(args.hours)
    job = w5_job()
    lab = Path(__file__).resolve().parent
    os.chdir(lab)
    (lab / T.LOG_DIR).mkdir(exist_ok=True)
    T._LOGFILE = str(lab / T.LOG_DIR / "batch_log.txt")
    man = T.load_manifest()

    def finish_outputs() -> None:
        if args.plot:
            try:
                T.log("  Figure: " + make_plot(man, lab, lab / FIG_DIR / "fig_w5.png"))
            except Exception as exc:
                T.log(f"  Figure could not be made: {type(exc).__name__}: {exc}")

    if args.summary_only:
        print_report(man, lab, args)
        finish_outputs()
        return 0

    args.notebook, nb_note = T.pick_notebook(lab, args.notebook)
    args.redo_ids = set(x.strip() for x in (args.redo or "").split(",") if x.strip())
    if args.min_free_gb is None:
        args.min_free_gb, basis = suggested_min_free_gb(man)
        free_note = f"free disk needed: {args.min_free_gb:g} GB ({basis})"
    else:
        free_note = f"free disk needed: {args.min_free_gb:g} GB (given)"
    preview = [job]
    done = T.is_finished(man, job, args)
    lo, hi = T.WIN["hours"] * T.MIN_PER_1000H / 1000.0, T.WIN["hours"] * T.MIN_PER_1000H_SLOW / 1000.0
    T.log(f"Working folder: {lab}")
    T.log(f"Notebook: {args.notebook}")
    T.log(f"Plan: one whole-year run, window W5 ({T.WIN['hours']} h from {T.WIN['start']}, water table {T.WIN['gw']:g} x bedrock): "
          + ("already finished, skipped." if done else "1 run to make."))
    T.log(f"  {job['run_id']:56s} {describe(job):70s} {'already done, skipped' if done else 'NEW RUN'}")
    T.log(f"  (it is checked against the W4 working run {W4_ID}, which must be finished in the manifest)")
    T.log(f"  tRIBS time alone: about {lo:.0f}-{hi:.0f} min (the plan's scaled estimate is about 68 min; this is the first real W5 timing); the notebook build takes a few minutes longer than on W4.")
    T.log(f"  {free_note}")

    results = T.preflight(args, lab, preview, nb_note, man)
    results.append(radar_check(args.radar_dir, T.WIN["hours"]))
    w4rec = man["runs"].get(W4_ID) or {}
    results.append(("OK" if w4rec.get("status") == "done" else "FAIL", f"the W4 working run {W4_ID} is finished in the manifest (status: {w4rec.get('status', 'not found')})"))
    results.append(("OK" if (lab / f"{W4_ID}.in").is_file() else "FAIL", f"the W4 working run's .in file is in lab/: {W4_ID}.in"))
    if args.hours == 8760:
        results.append(("WARN", "--hours 8760 runs 1 h past the end of the weather files (the notebook warns about this); tRIBS may fail on the last hour. The default 8759 avoids it."))
    else:
        results.append(("OK", f"run length {args.hours} h ends Sep 30 2017 23:00, the last hour of the weather files (the scorer still counts it as the whole year)"))
    T.log()
    for level, msg in results:
        T.log(f"  {level:5s} {msg}")
    fails = [m for lv, m in results if lv == "FAIL"]
    T.log()
    if args.dry_run:
        T.log("DRY RUN: nothing was built or run. " + ("All checks passed." if not fails else f"{len(fails)} check(s) FAILED."))
        return 1 if fails else 0
    if fails:
        T.log(f"{len(fails)} check(s) FAILED, so nothing was started. Fix them and run again.")
        return 1

    try:
        T.take_lock()
    except T.Abort as exc:
        T.log(str(exc))
        return 1
    man["baseline_id"] = T.BASELINE_ID
    man["w5_last_started"] = T.now_iso()
    T.save_manifest(man)
    T.install_term_handler()

    t_batch, code = time.time(), 0
    try:
        try:
            T.do_run(job, args, lab, man)
        except T.RunError as exc:
            man["runs"].setdefault(W5_ID, {}).update(status="failed", error=str(exc), finished=T.now_iso())
            T.save_manifest(man)
            T.log(f"{W5_ID}: FAILED. {exc}")
            code = 1
    except T.Abort as exc:
        T.log(f"BATCH STOPPED: {exc}")
        code = 1
    except KeyboardInterrupt:
        T.log("INTERRUPTED. Run the same command again to resume.")
        for rid, rec in man["runs"].items():
            if rec.get("status") in ("building", "running"):
                rec["status"] = "interrupted"
        T.save_manifest(man)
        code = 130
    finally:
        T.save_manifest(man)

    T.log()
    T.log(f"Batch time: {T.fmt_min(time.time() - t_batch)}")
    rec = man["runs"].get(W5_ID) or {}
    T.log(f"{W5_ID:56s} {rec.get('status', '?'):16s} build {rec.get('build_min', '')} min  tRIBS {rec.get('tribs_min', '')} min")
    if not args.build_only:
        print_report(man, lab, args)
        finish_outputs()
    return code


if __name__ == "__main__":
    sys.exit(main())
