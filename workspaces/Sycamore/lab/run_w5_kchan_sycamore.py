#!/usr/bin/env python3
"""
run_w5_kchan_sycamore.py  --  THREE whole-water-year (W5: Oct 1 2016 - Sep 30 2017 23:00, 8,759 h) runs of draw 12 / k10 AS DRAWN, one after another in one
command: channel loss OFF (the control), then channel loss ON at Kchan 0.5 mm/hr, then at Kchan 0.2 mm/hr.   About 3.2 - 3.5 hours in all.

THE SETTING ("draw 12 / k10 as drawn" = the baseline of the W4 scripts, with NO re-tune):
    Ks x0.187, f x0.5874, PsiB x0.2939, As 1, Au 1, routing k10 (cv 6.39, flowexp 0.2973, n 0.067), width coefficient 2.33, bedrock +1 m,
    water table 0.95 x bedrock at Oct 1, snow ON.      (NOT the re-tuned working setting: no f x0.15, no As 300, no Ks x1.03.)
THE THREE RUNS (they differ ONLY in the two channel-loss keywords, so any difference between them is Kchan's):
    W5_d012_k10_base_nocl          OPTPERCOLATION 0                          (the control; runs first)
    W5_d012_k10_base_kchan0p5      OPTPERCOLATION 1, CHANNELCONDUCTIVITY 0.5
    W5_d012_k10_base_kchan0p2      OPTPERCOLATION 1, CHANNELCONDUCTIVITY 0.2
Your earlier W5 run (W5_d012_k10_f0p15_as300_ksx1p03_nocl, the re-tuned setting) is not touched.

WHY ONE AFTER ANOTHER, NOT AT THE SAME TIME: every run rebuilds the same data/model/ folder, and the shared lock refuses a second tRIBS while one is
running. So the three runs are chained in one command (about 62-65 min each, the last W5 run took 62.1 min).

WHAT QUESTION IT ANSWERS: the W5 run at the re-tuned setting turned July rain into outflow at about the January rate while the gauge shows almost
nothing in Jul-Sep (model 19.0 M m3, gauge 0.001). Channel loss is the one process that can remove flow in the channel. Does a constant Kchan remove
the summer false flow MORE than it removes the January - March flow the gauge confirms?

READING RULES (written BEFORE the runs; you accepted them as proposed on 2026-11-03 03:48). Each Kchan run is compared with the loss-off run above:
    RULE A  Jul-Sep volume falls at least twice as many percent as the Jan-Mar volume (and it does fall).
    RULE B  Jan-Mar volume stays at least 75% of the loss-off run's Jan-Mar volume.
    The script prints PASS / MISS for A and B. It decides nothing. If you change a rule, change it here and in the log BEFORE you look at results.
    Volumes are sums of hourly mean flows from the outlet file, M m3 (the scorer integrates the raw steps, so its numbers differ slightly).
    The report also PRINTS, without judging: Oct-Dec and Jul-Sep volumes next to the gauge, zero-flow hours, days above 0.13 m3/s, the control's ten
    biggest Apr-Sep days in every run, month-by-month volumes, channel loss, and the scorer's whole-year ratio and KGE.

WHAT IT CHECKS BEFORE IT STARTS (the W5 script's checks, plus): the W4 baseline run (W4_d012_k10_base_nocl) is finished in the manifest and its .in file is
in lab/ (every new .in must equal it in every keyword except the window, the output name and the two channel-loss keywords); the radar folder has all the
hourly files the whole year needs; free disk space; no other tRIBS running. After a Kchan run the scorer's channel loss must be above zero, or a warning prints.

HOW TO RUN. Put it in workspaces/Sycamore/lab/ next to run_w5_sycamore.py, run_retune_w4_sycamore.py and run_tail_screen_sycamore.py (all three must be there).
Foreground, one command at a time:
    python -u run_w5_kchan_sycamore.py --dry-run                      # checks only, builds nothing
    python -u run_w5_kchan_sycamore.py --build-only --limit 1         # builds the control's inputs and checks its .in file; does NOT run tRIBS (a few minutes)
    python -u run_w5_kchan_sycamore.py                                # the three runs (about 3.2 - 3.5 hours); safe to run again if the connection drops
    python -u run_w5_kchan_sycamore.py --summary-only --plot          # tables and figure from the finished runs
Do NOT run the notebook (or anything that rebuilds data/model/) while it runs. It changes no setting, makes no recommendation, and does not touch the design doc.
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
import run_w5_sycamore as W                      # noqa: E402  (importing it builds nothing; it brings the tested W4 scripts with it)
R = W.R
T = W.T                                           # the tested W4 script that builds, runs and scores

# --------------------------------------------------------------------------------------------
# Names
# --------------------------------------------------------------------------------------------
BASE_W4_ID = f"{T.RUN_PREFIX}_base_nocl"          # W4_d012_k10_base_nocl: draw 12 / k10 as drawn, W4, loss off (finished; the .in reference)
CONTROL_ID = BASE_W4_ID.replace("W4_", "W5_", 1)  # W5_d012_k10_base_nocl
DEFAULT_KCHAN = [0.5, 0.2]                        # mm/hr, in the order they run after the control
RULE_A_FACTOR = 2.0
RULE_B_MIN_FRACTION = 0.75
SUMMARY_NAME = "w5_kchan_summary.csv"
MONTHS_NAME = "w5_kchan_months.csv"
FIG_DIR = "w5_kchan_figs"
SEASON_RANGE = {n: (a, b) for n, a, b in W.SEASONS}
NAN = float("nan")
PLAN: list = []


def vtag(x) -> str:
    return T.vtag(x)


def describe(job: dict) -> str:
    """The text the W4 script prints for a run (it looks this function up by name, so it is replaced here)."""
    base = "draw 12 / k10 as drawn (Ks x0.187, f x0.5874, As 1)"
    if job["step"] == "kchan":
        return f"{base}, channel loss ON, Kchan {job['kchan']:g} mm/hr, whole water year"
    return f"{base}, channel loss off (control), whole water year"


def control_job() -> dict:
    j = dict(T.make_base_job())
    j["run_id"] = CONTROL_ID
    j["step"] = "w5"                  # not 'base': a 'base' job is not compared with the reference .in
    return j


def kchan_job(k: float) -> dict:
    j = dict(control_job())
    j["changes"] = {"OPTPERCOLATION": 1, "CHANNELCONDUCTIVITY": float(k)}
    j["run_id"] = CONTROL_ID[:-len("_nocl")] + f"_kchan{vtag(k)}"
    j["step"] = "kchan"
    j["kchan"] = float(k)
    return j


def set_plan(kchans: list) -> None:
    PLAN.clear()
    PLAN.append(control_job())
    PLAN.extend(kchan_job(k) for k in kchans)


def planned() -> list:
    return list(PLAN)


def label_of(job: dict) -> str:
    return "loss off" if job["step"] == "w5" else f"Kchan {job['kchan']:g}"


def configure(hours: int = 8759) -> None:
    T.WIN.clear()
    T.WIN.update(W.W5)
    T.WIN["hours"] = int(hours)
    T.BASELINE_ID = BASE_W4_ID        # every new .in is checked against the W4 baseline's .in, and that run must be finished
    T.SUMMARY_CSV = SUMMARY_NAME
    T.describe = describe


def parse_kchan(text: str) -> list:
    out = []
    for part in str(text).split(","):
        part = part.strip()
        if not part:
            continue
        try:
            v = float(part)
        except ValueError:
            raise ValueError(f"'{part}' is not a number")
        if not (math.isfinite(v) and 0 < v <= 1000):
            raise ValueError(f"{part} is not allowed: Kchan must be above 0 and at most 1000 mm/hr")
        if v not in out:
            out.append(v)
    if not out or len(out) > 6:
        raise ValueError("give between 1 and 6 Kchan values, for example 0.5,0.2")
    return out


# --------------------------------------------------------------------------------------------
# Numbers from the outlet files
# --------------------------------------------------------------------------------------------
def pct_change(new, old) -> float:
    return 100.0 * (new / old - 1.0) if (new is not None and old is not None and math.isfinite(new) and math.isfinite(old) and old != 0) else NAN


def season_vols(q: pd.Series) -> dict:
    return {n: W.vol_mm3(W.hourly(q, a, b)) for n, (a, b) in SEASON_RANGE.items()}


def day_stats(q: pd.Series, day) -> tuple:
    h = W.hourly(q, day, day + pd.Timedelta(days=1))
    return (W.vol_mm3(h), float(h.max())) if len(h) else (NAN, NAN)


def rules(ctrl_v: dict, k_v: dict) -> dict:
    """Rules A and B for one Kchan run against the control. Each value: (True / False / None, text)."""
    res = {}
    js_c, js_k = ctrl_v.get("Jul-Sep", NAN), k_v.get("Jul-Sep", NAN)
    jm_c, jm_k = ctrl_v.get("Jan-Mar", NAN), k_v.get("Jan-Mar", NAN)
    js_fall, jm_fall = -pct_change(js_k, js_c), -pct_change(jm_k, jm_c)
    if math.isfinite(js_fall) and math.isfinite(jm_fall):
        res["A"] = (bool(js_fall > 0 and js_fall >= RULE_A_FACTOR * jm_fall),
                    f"Jul-Sep volume falls {js_fall:.1f}% ({js_c:.2f} -> {js_k:.2f} M m3); Jan-Mar volume falls {jm_fall:.1f}% ({jm_c:.2f} -> {jm_k:.2f}); "
                    f"needed: Jul-Sep fall at least {RULE_A_FACTOR:g} x the Jan-Mar fall")
    else:
        res["A"] = (None, "not computed (a volume is missing or the control's Jul-Sep volume is zero)")
    if math.isfinite(jm_c) and math.isfinite(jm_k) and jm_c > 0:
        res["B"] = (bool(jm_k >= RULE_B_MIN_FRACTION * jm_c),
                    f"Jan-Mar volume {jm_k:.2f} M m3 against {jm_c:.2f} loss-off = {100.0 * jm_k / jm_c:.1f}% (needed: at least {100 * RULE_B_MIN_FRACTION:.0f}%)")
    else:
        res["B"] = (None, "not computed")
    return res


def yn3(v) -> str:
    return "n.a." if v is None else ("PASS" if v else "MISS")


def pct_text(x, w=12) -> str:
    return f"{'-':>{w}}" if not math.isfinite(x) else f"{x:>+{w - 1}.1f}%"


def f2(x, nd=2, w=9) -> str:
    return f"{'-':>{w}}" if (x is None or not math.isfinite(x)) else f"{x:{w}.{nd}f}"


# --------------------------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------------------------
def print_report(man: dict, lab: Path, args) -> None:
    csv_scores = T.load_scores(lab)
    done = []
    for j in planned():
        rec = man["runs"].get(j["run_id"]) or {}
        qf = lab / "results" / f"{j['run_id']}_Outlet.qout"
        if rec.get("status") in ("done", "unscored") and qf.is_file():
            done.append((j, rec, T.read_qout(qf, T.RUN_START_TS)))
    T.log()
    T.log("W5 KCHAN COMPARISON  draw 12 / k10 as drawn (Ks x0.187, f x0.5874, PsiB x0.2939, As 1, Au 1), whole water year 2017 (Oct 1 2016 - Sep 30 2017 23:00)")
    T.log("  'loss off' = the control (OPTPERCOLATION 0). 'Kchan x' = channel loss ON at x mm/hr. Volumes: sums of hourly mean flows, M m3.")
    if not done:
        T.log("  (no finished runs yet, so there is nothing to report)")
        return
    try:
        gauge_q = T.read_obs(lab / T.INIT_DATA / T.OBS_FILE)
    except Exception as exc:
        T.log(f"  (gauge file could not be read, so the tables are skipped: {exc})")
        return
    ctrl = next((d for d in done if d[0]["step"] == "w5"), None)
    vols = {d[0]["run_id"]: season_vols(d[2]) for d in done}
    gv = season_vols(gauge_q)
    months = {d[0]["run_id"]: W.month_rows(d[2], gauge_q) for d in done}
    mon0 = next(iter(months.values()))
    labels = [label_of(d[0]) for d in done]
    scores = {d[0]["run_id"]: R.scores_row(man, lab, d[0]["run_id"], csv_scores) for d in done}

    T.log()
    T.log("  RUNS")
    T.log(f"    {'run':<10}{'status':>9}{'build min':>10}{'tRIBS min':>10}{'output MB':>10}{'chan loss':>11}{'year vol ratio':>16}{'KGE (gate first)':>18}")
    for j, rec, q in done:
        sc = scores[j["run_id"]]
        T.log(f"    {label_of(j):<10}{rec.get('status', '?'):>9}{str(rec.get('build_min', '-')):>10}{str(rec.get('tribs_min', '-')):>10}{str(rec.get('output_mb', '-')):>10}"
              f"{f2(T.fnum(sc, 'proc__channel_loss_Mm3'), 3, 11)}{f2(T.fnum(sc, 'span__vol_ratio'), 3, 16)}{f2(T.fnum(sc, 'span__kge2012'), 3, 18)}")
    for j, rec, q in done:
        if j["step"] == "kchan":
            cl = T.fnum(scores[j["run_id"]], "proc__channel_loss_Mm3")
            if math.isfinite(cl) and cl <= 0:
                T.log(f"    WARNING: {label_of(j)} shows channel loss of {cl:g} M m3, so the loss was not applied; do not read this run as a loss run.")
    T.log("    (chan loss = M m3 removed by the channel, from the scorer / .mrf)")

    T.log()
    T.log("  CONSISTENCY  (the control over Oct 1 - Apr 1 against the finished W4 baseline run, hourly means)")
    w4f = lab / "results" / f"{BASE_W4_ID}_Outlet.qout"
    if ctrl is None:
        T.log("    the control run is not finished, so no comparison.")
    elif not w4f.is_file():
        T.log(f"    the W4 baseline's outlet file {w4f.name} is not in results/, so no comparison.")
    else:
        c = W.consistency(ctrl[2], T.read_qout(w4f, T.RUN_START_TS))
        if c:
            T.log(f"    {c['hours']} hours compared; largest difference {c['max_abs_diff']:.4f} m3/s (at {c['t_max']:%b %d %H:%M}); volume W5 {c['vol_w5']:.3f} vs W4 {c['vol_w4']:.3f} M m3 (ratio {c['ratio']:.4f})")
        else:
            T.log("    no hours in common")

    T.log()
    T.log("  SEASONS  (M m3)")
    T.log(f"    {'season':<9}{'gauge':>9}" + "".join(f"{lb:>12}" for lb in labels))
    for name in SEASON_RANGE:
        T.log(f"    {name:<9}{f2(gv[name], 3)}" + "".join(f"{f2(vols[d[0]['run_id']][name], 3, 12)}" for d in done))
    tg = sum(r["gauge_Mm3"] for r in mon0)
    T.log(f"    {'year':<9}{f2(tg, 3)}" + "".join(f"{f2(sum(r['model_Mm3'] for r in months[d[0]['run_id']]), 3, 12)}" for d in done))

    if ctrl is not None and len(done) > 1:
        T.log()
        T.log("  CHANGE FROM THE LOSS-OFF RUN, % of its volume")
        T.log(f"    {'season':<9}" + "".join(f"{label_of(d[0]):>12}" for d in done if d[0]["step"] == "kchan"))
        cv = vols[ctrl[0]["run_id"]]
        for name in SEASON_RANGE:
            T.log(f"    {name:<9}" + "".join(pct_text(pct_change(vols[d[0]['run_id']][name], cv[name])) for d in done if d[0]["step"] == "kchan"))

    T.log()
    if ctrl is None:
        T.log("  RULES: the loss-off control is not finished, so the rules cannot be computed.")
    else:
        T.log("  RULES  (written before the runs; PASS / MISS only prints; you decide)")
        for j, rec, q in done:
            if j["step"] != "kchan":
                continue
            res = rules(vols[ctrl[0]["run_id"]], vols[j["run_id"]])
            T.log(f"    {label_of(j)}")
            for key in ("A", "B"):
                st, txt = res[key]
                T.log(f"      RULE {key}  {yn3(st):>4}  {txt}")
            both = [res[k][0] for k in ("A", "B")]
            T.log("      A and B together: " + ("both PASS" if all(v is True for v in both) else ("at least one MISS" if any(v is False for v in both) else "not all computed")))
        if not any(d[0]["step"] == "kchan" for d in done):
            T.log("    (no Kchan run is finished yet)")

    T.log()
    T.log("  FALSE FLOW and ZERO-FLOW  (printed, not judged. days above = days whose highest hourly flow is over 0.13 m3/s; zero-flow = hourly mean under 0.001 m3/s)")
    T.log(f"    {'':<22}{'gauge':>9}" + "".join(f"{lb:>12}" for lb in labels))
    for name in ("Oct-Dec", "Jul-Sep"):
        a, b = SEASON_RANGE[name]
        gh = W.hourly(gauge_q, a, b)
        T.log(f"    {name + ' volume, M m3':<22}{f2(gv[name], 3)}" + "".join(f"{f2(vols[d[0]['run_id']][name], 3, 12)}" for d in done))
        T.log(f"    {name + ' days above':<22}{W.days_above(gh, W.THR_DAYS):>9d}" + "".join(f"{W.days_above(W.hourly(d[2], a, b), W.THR_DAYS):>12d}" for d in done))
        T.log(f"    {name + ' zero-flow %':<22}{f2(W.zero_hours_pct(gh), 1)}" + "".join(f"{f2(W.zero_hours_pct(W.hourly(d[2], a, b)), 1, 12)}" for d in done))
    T.log(f"    {'Jul-Sep max, m3/s':<22}{f2(float(W.hourly(gauge_q, *SEASON_RANGE['Jul-Sep']).max()), 2)}"
          + "".join(f"{f2(float(W.hourly(d[2], *SEASON_RANGE['Jul-Sep']).max()), 2, 12)}" for d in done))

    T.log()
    T.log("  MONTHS  (model volume, M m3; the gauge's January lacks Jan 1-2, an outage)")
    T.log(f"    {'month':<9}{'gauge':>9}" + "".join(f"{lb:>12}" for lb in labels))
    for i, r in enumerate(mon0):
        T.log(f"    {r['month']:<9}{f2(r['gauge_Mm3'], 3)}" + "".join(f"{f2(months[d[0]['run_id']][i]['model_Mm3'], 3, 12)}" for d in done))

    if ctrl is not None:
        T.log()
        T.log("  BIGGEST DAYS  (the control's ten biggest Apr 1 - Oct 1 days, and the same days in each run: volume M m3 / highest hourly flow m3/s)")
        T.log(f"    {'day':<8}{'gauge':>14}" + "".join(f"{lb:>16}" for lb in labels))
        for r in W.biggest_days(ctrl[2], gauge_q, "2017-04-01 00:00", W.Y1):
            day = r["day"]
            g = day_stats(gauge_q, day)
            T.log(f"    {day:%b %d}  {g[0]:>6.3f} /{g[1]:>5.2f}" + "".join(f"{day_stats(d[2], day)[0]:>9.3f} /{day_stats(d[2], day)[1]:>5.1f}" for d in done))

    T.log()
    tot = sum(float(rec.get("tribs_min") or 0) + float(rec.get("build_min") or 0) for _, rec, _ in done)
    T.log(f"  PROCESS   {len(done)} run(s) finished; build plus tRIBS time {tot:.1f} min in all; {T.disk_free_gb(lab):.1f} GB free now")
    rows = []
    for j, rec, q in done:
        sc = scores[j["run_id"]]
        v = vols[j["run_id"]]
        rows.append(dict(run_id=j["run_id"], label=label_of(j), kchan=j.get("kchan", ""), **{f"{n}_Mm3": v[n] for n in SEASON_RANGE},
                         year_Mm3=sum(r["model_Mm3"] for r in months[j["run_id"]]), span_vol_ratio=T.fnum(sc, "span__vol_ratio"),
                         span_kge2012=T.fnum(sc, "span__kge2012"), channel_loss_Mm3=T.fnum(sc, "proc__channel_loss_Mm3"),
                         tribs_min=rec.get("tribs_min", ""), build_min=rec.get("build_min", "")))
    pd.DataFrame(rows).to_csv(lab / SUMMARY_NAME, index=False)
    mrows = []
    for i, r in enumerate(mon0):
        row = dict(month=r["month"], gauge_Mm3=r["gauge_Mm3"])
        for d in done:
            row[f"model_Mm3_{label_of(d[0]).replace(' ', '_')}"] = months[d[0]["run_id"]][i]["model_Mm3"]
        mrows.append(row)
    pd.DataFrame(mrows).to_csv(lab / MONTHS_NAME, index=False)
    T.log(f"  Wrote {SUMMARY_NAME} and {MONTHS_NAME}. These tables decide nothing.")


def make_plot(man: dict, lab: Path, out: Path) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    gauge_q = T.read_obs(lab / T.INIT_DATA / T.OBS_FILE)
    runs = []
    for j in planned():
        qf = lab / "results" / f"{j['run_id']}_Outlet.qout"
        if (man["runs"].get(j["run_id"]) or {}).get("status") in ("done", "unscored") and qf.is_file():
            runs.append((j, T.read_qout(qf, T.RUN_START_TS)))
    if not runs:
        return "no finished runs, so no figure"
    colors = ["#2a78d6", "#eb6834", "#4a9d6b", "#8b5fbf"]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 7.5), gridspec_kw={"height_ratios": [2, 3]})
    mon = {j["run_id"]: W.month_rows(q, gauge_q) for j, q in runs}
    first = next(iter(mon.values()))
    x = np.arange(len(first))
    n = len(runs) + 1
    wd = 0.8 / n
    ax1.bar(x - 0.4 + wd / 2, [r["gauge_Mm3"] for r in first], wd, color="#8a8984", label="USGS gauge")
    for i, (j, q) in enumerate(runs):
        ax1.bar(x - 0.4 + wd * (i + 1.5), [r["model_Mm3"] for r in mon[j["run_id"]]], wd, color=colors[i % len(colors)], label=label_of(j))
    ax1.set_xticks(x)
    ax1.set_xticklabels([r["month"].split()[0] for r in first])
    ax1.set_ylabel("Monthly volume (M m$^3$)")
    ax1.set_title("Draw 12 / k10 as drawn, whole water year 2017: channel loss off and on", loc="left")
    ax1.legend(frameon=False, ncol=n)
    a, b = SEASON_RANGE["Jul-Sep"]
    gh = W.hourly(gauge_q, a, b)
    ax2.plot(gh.index, gh.clip(lower=1e-3), color="#8a8984", lw=1.0, label="USGS gauge")
    for i, (j, q) in enumerate(runs):
        mh = W.hourly(q, a, b)
        ax2.plot(mh.index, mh.clip(lower=1e-3), color=colors[i % len(colors)], lw=0.9, label=label_of(j))
    ax2.set_yscale("log")
    ax2.set_ylim(1e-3, 600)
    ax2.set_ylabel("Outlet flow, hourly mean (m$^3$/s), log scale")
    ax2.set_title("July - September", loc="left")
    ax2.legend(frameon=False, loc="upper right", ncol=n)
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
# Arguments, checks, main
# --------------------------------------------------------------------------------------------
def suggested_min_free_gb(man: dict, n_runs: int) -> tuple:
    mb = (man["runs"].get(BASE_W4_ID) or {}).get("output_mb")
    if mb:
        return max(4.0, round(n_runs * 2.2 * float(mb) / 1000.0 + 1.0, 1)), f"the W4 baseline wrote {float(mb):.0f} MB; each whole-year run writes about twice that, for {n_runs} runs, plus a margin"
    return 8.0, "the size of the W4 baseline's output is not in the manifest, so a cautious default"


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description="Three whole-water-year (W5) runs of draw 12 / k10 as drawn: loss off, Kchan 0.5, Kchan 0.2 (see the top of this file).")
    ap.add_argument("--kchan", default=",".join(f"{k:g}" for k in DEFAULT_KCHAN), help="Kchan value(s) in mm/hr, comma-separated, run in this order after the control (default %(default)s). Finished runs are skipped")
    ap.add_argument("--dry-run", action="store_true", help="check everything, list the plan, build nothing")
    ap.add_argument("--summary-only", action="store_true", help="print the tables from the finished runs, then stop")
    ap.add_argument("--build-only", action="store_true", help="build and check the inputs but do not run tRIBS (use with --limit 1)")
    ap.add_argument("--limit", type=int, default=None, help="make only the first N runs that still need work (finished runs do not count)")
    ap.add_argument("--plot", action="store_true", help=f"make {FIG_DIR}/fig_w5_kchan.png from the finished runs")
    ap.add_argument("--hours", type=int, default=8759, help="run length in hours (default 8759: the weather files end 1 h before a full 8760 h year would; 8760 is allowed)")
    ap.add_argument("--radar-dir", default=W.RADAR_DIR_DEFAULT, help="folder of the hourly radar files, for the check that all of them are there")
    ap.add_argument("--redo", default=None, help="comma-separated run ids to run again even if done")
    ap.add_argument("--timeout-min", type=float, default=None, help="stop a tRIBS run after this many minutes (default: 3.5 times the expected time, about 3.7 hours)")
    ap.add_argument("--cell-timeout", type=int, default=1800, help="seconds allowed for one notebook cell")
    ap.add_argument("--min-free-gb", type=float, default=None, help="stop if less than this much disk is free (default: worked out from the W4 baseline's output size)")
    ap.add_argument("--heartbeat-min", type=float, default=T.HEARTBEAT_MIN, help="print a progress line every this many minutes while tRIBS runs (default 5)")
    ap.add_argument("--band-lo", type=float, default=T.BAND_LO)
    ap.add_argument("--band-hi", type=float, default=T.BAND_HI)
    ap.add_argument("--notebook", default=None, help="notebook to run (default: Make_Sycamore_Model.ipynb, else Make_Sycamore_Model_clean.ipynb)")
    ap.add_argument("--scorer", default=T.SCORER)
    ap.add_argument("--tribs", default=os.environ.get("TRIBS_BIN") or (T.TRIBS_DEFAULT if os.path.isfile(T.TRIBS_DEFAULT) else shutil.which("tRIBS") or T.TRIBS_DEFAULT))
    ap.add_argument("--skip-kernel-check", action="store_true", help="skip the pytRIBS import test in --dry-run")
    a = ap.parse_args(argv)
    a.steps = "w5kchan"
    return a


def main(argv=None) -> int:
    lost = R.missing_from_w4_script()
    if lost:
        print("run_tail_screen_sycamore.py in this folder is an older copy: it lacks " + ", ".join(lost) + ". Put the newest run_tail_screen_sycamore.py in lab/ and run again.")
        return 2
    args = parse_args(argv)
    try:
        kchans = parse_kchan(args.kchan)
    except ValueError as exc:
        print(f"--kchan: {exc}")
        return 2
    if not (24 <= args.hours <= 8760):
        print(f"--hours {args.hours}: must be between 24 and 8760")
        return 2
    configure(args.hours)
    set_plan(kchans)
    lab = Path(__file__).resolve().parent
    os.chdir(lab)
    (lab / T.LOG_DIR).mkdir(exist_ok=True)
    T._LOGFILE = str(lab / T.LOG_DIR / "batch_log.txt")
    man = T.load_manifest()

    def finish_outputs() -> None:
        if args.plot:
            try:
                T.log("  Figure: " + make_plot(man, lab, lab / FIG_DIR / "fig_w5_kchan.png"))
            except Exception as exc:
                T.log(f"  Figure could not be made: {type(exc).__name__}: {exc}")

    if args.summary_only:
        print_report(man, lab, args)
        finish_outputs()
        return 0

    args.notebook, nb_note = T.pick_notebook(lab, args.notebook)
    args.redo_ids = set(x.strip() for x in (args.redo or "").split(",") if x.strip())
    preview = planned()
    n_done = sum(1 for j in preview if T.is_finished(man, j, args))
    n_todo = len(preview) - n_done
    if args.min_free_gb is None:
        args.min_free_gb, basis = suggested_min_free_gb(man, max(n_todo, 1))
        free_note = f"free disk needed: {args.min_free_gb:g} GB ({basis})"
    else:
        free_note = f"free disk needed: {args.min_free_gb:g} GB (given)"
    per_lo, per_hi = T.WIN["hours"] * T.MIN_PER_1000H / 1000.0, T.WIN["hours"] * T.MIN_PER_1000H_SLOW / 1000.0
    T.log(f"Working folder: {lab}")
    T.log(f"Notebook: {args.notebook}")
    T.log(f"Plan: draw 12 / k10 as drawn, window W5 ({T.WIN['hours']} h from {T.WIN['start']}, water table {T.WIN['gw']:g} x bedrock): loss-off control, then Kchan "
          f"{', '.join(f'{k:g}' for k in kchans)} mm/hr; {n_todo} run(s) to make" + (f", {n_done} already finished" if n_done else "") + ".")
    for j in preview:
        what = "already done, skipped" if T.is_finished(man, j, args) else "NEW RUN"
        T.log(f"  {j['run_id']:42s} {describe(j):100s} {what}")
    T.log(f"  (the new .in files are checked against the finished W4 baseline run {BASE_W4_ID}, keyword by keyword)")
    T.log(f"  tRIBS time alone: about {n_todo * per_lo / 60:.1f}-{n_todo * per_hi / 60:.1f} hours for {n_todo} run(s) (the script's estimate is {per_lo:.0f}-{per_hi:.0f} min per run; "
          f"the last W5 run took 62 min); each run also needs a notebook build of about 1 min.")
    T.log(f"  {free_note}")

    results = T.preflight(args, lab, preview, nb_note, man)
    results.append(W.radar_check(args.radar_dir, T.WIN["hours"]))
    w4rec = man["runs"].get(BASE_W4_ID) or {}
    results.append(("OK" if w4rec.get("status") == "done" else "FAIL", f"the W4 baseline run {BASE_W4_ID} is finished in the manifest (status: {w4rec.get('status', 'not found')})"))
    results.append(("OK" if (lab / f"{BASE_W4_ID}.in").is_file() else "FAIL", f"the W4 baseline run's .in file is in lab/: {BASE_W4_ID}.in"))
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
    man["w5_kchan_last_started"] = T.now_iso()
    T.save_manifest(man)
    T.install_term_handler()

    t_batch, code = time.time(), 0
    budget = {"left": args.limit}

    def run_batch(todo: list) -> int:
        consecutive_fail, rc = 0, 0
        for j in todo:
            rid = j["run_id"]
            if not T.is_finished(man, j, args):
                if budget["left"] is not None:
                    if budget["left"] <= 0:
                        T.log(f"--limit reached; {rid} and the runs after it were not started.")
                        break
                    budget["left"] -= 1
            try:
                T.do_run(j, args, lab, man)
                consecutive_fail = 0
            except T.RunError as exc:
                man["runs"].setdefault(rid, {}).update(status="failed", error=str(exc), finished=T.now_iso())
                T.save_manifest(man)
                T.log(f"{rid}: FAILED. {exc}")
                consecutive_fail += 1
                rc = 1
                if consecutive_fail >= 2:
                    T.log("Two runs failed in a row, so the batch stops (the cause is probably the same).")
                    break
        return rc

    try:
        code = run_batch([PLAN[0]])                   # the control always comes first: the Kchan runs are only readable against it
        if code != 0:
            T.log("The Kchan runs were not started because the control failed or was interrupted. Fix that, then run the same command again.")
        else:
            code = run_batch(PLAN[1:]) or code
    except T.Abort as exc:
        T.log(f"BATCH STOPPED: {exc}")
        code = 1
    except KeyboardInterrupt:
        T.log("BATCH INTERRUPTED. Run the same command again to resume; finished runs are skipped.")
        for rid, rec in man["runs"].items():
            if rec.get("status") in ("building", "running"):
                rec["status"] = "interrupted"
        T.save_manifest(man)
        code = 130
    finally:
        T.save_manifest(man)

    T.log()
    T.log(f"Batch time: {T.fmt_min(time.time() - t_batch)}")
    T.log(f"{'run':42s} {'status':16s} {'build':>7s} {'tRIBS':>7s}")
    for j in planned():
        rec = man["runs"].get(j["run_id"]) or {}
        T.log(f"{j['run_id']:42s} {rec.get('status', 'not started'):16s} {str(rec.get('build_min', '')):>7s} {str(rec.get('tribs_min', '')):>7s}")
    if not args.build_only:
        print_report(man, lab, args)
        finish_outputs()
    return code


if __name__ == "__main__":
    sys.exit(main())
