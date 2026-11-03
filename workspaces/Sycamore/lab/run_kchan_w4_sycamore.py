"""
run_kchan_w4_sycamore.py  --  The channel-conductivity (Kchan) probe on the WORKING setting, on the full tail-screen window W4
(Oct 1 2016 - Apr 1 2017, 4,368 h, water table 0.95 x bedrock). About 36 minutes per run.

THE WORKING SETTING (the loss-off run this probe is compared with, already finished on your machine):
    Ks x1.03 (KS_MULT 0.19261), f x0.15, As 300, Au 1, PsiB x0.2939, routing k10 (cv 6.39, flowexp 0.2973, n 0.067), width coefficient 2.33,
    bedrock +1 m, snow ON, channel loss OFF.   Run name: W4_d012_k10_f0p15_as300_ksx1p03_nocl
THE PROBE RUN: the same setting with channel loss ON (OPTPERCOLATION 1) and the streambed conductivity set to Kchan (mm/hr).
    Default: Kchan 0.5 mm/hr.   Run name: W4_d012_k10_f0p15_as300_ksx1p03_kchan0p5
    (This replaces run_w4_kchan_sycamore.py, which builds the OLD setting: f x0.5874, As 1. Do not run that one.)

WHY: the Oct-Dec false flow (the model makes flow in dry months when the gauge has almost none) is the largest error left. Loss-off, it is
about 3.3 times the gauge's volume. Channel loss is the physical process that removes small flows on their way down the channel. The question
the probe answers is whether Kchan can cut the Oct-Dec false flow WITHOUT also damaging the Feb storms, the tail and the zero-flow record.
The Stage 3 test (a different setting) said Kchan 0.5 cut a fall storm by 25% but a February storm by only 8-12%, and left the outlet at zero
flow 82.5% of hours (4.5% with loss off). That is the pattern to look for, or not, on the working setting.

READING RULES (written BEFORE the run, from the mini re-write; you accepted these numbers on 2026-11-02 21:25):
    RULE A  Oct-Dec volume falls at least twice as much, in percent, as the Feb 27-28 volume.     (against the loss-off run above)
    RULE B  Zero-flow hours stay under 50% of the 4,368 hours of W4.                            (the loss-off run's own figure is printed beside it)
    RULE C  Tail volume (Jan 26 - Feb 2) stays at least half of the loss-off run's (0.37 -> 0.185 or more).
    SHOULDER (no number, you judge): the early shoulder on Feb 28 (about 21 m3/s in the loss-off run; the gauge stays at or under 6.1 until 08:05)
            should fall clearly toward the gauge's. Printed for both runs.
    The script prints PASS / MISS for A, B, C. It decides nothing: whether a run is worth keeping, and what Kchan to try next, is your call.
    A second Kchan value is chosen AFTER you read the first (for example 2, if 0.5 is too harsh, or 0.25, if it is too mild): give it with --kchan.
    The script uses the numbers above only to print PASS / MISS; if you change a rule, change it here and in the log BEFORE you look at results.

WHAT IT PRINTS (after the runs, or with --summary-only):
    STORMS         flood week, Feb 27-28, Feb 18-19, second January storm: volume, peak, scorer arrival; and the change from the loss-off run.
    RULES          A, B, C above, with the numbers behind them.
    FALSE FLOW     Oct-Dec model volume against the gauge, days above 0.13 m3/s, zero-flow hours (whole W4 and Oct-Dec), channel loss volume.
    TIMING         hours, model minus gauge, for t20 (the scorer's rule), t50 (the proposed arrival rule) and t80, plus the highest flow BEFORE the
                   main rise, for model and gauge. Here "before the main rise" means before the time the GAUGE first reaches 10% of its own peak
                   in that event (the same cut time for both). A cut at 50% would only measure the rising limb itself, about half the peak, and
                   could never show an early shoulder. The cut time is printed in the CSV (kx_f27_gauge_rise_start).
    TAIL           24-hour mean flows 7, 14 and 21 days after the gauge flood peak; days in the 0.3-1 m3/s band (gauge 16 of 23); tail volumes.
    PROCESS        water table, channel loss, build and tRIBS minutes.
It also writes kchan_w4_summary.csv and (with --plot) kchan_w4_figs/fig_kchan_w4.png.

WHAT IT DOES NOT DO. It does not change Ks, f, As, routing or the width coefficient. It makes no recommendation. It does not touch the design doc.
Zero-flow hours here means hours whose mean outlet flow is under 0.001 m3/s, counted by this script for BOTH runs the same way. The Stage 3 numbers
(82.5% and 4.5%) were counted at a different setting and possibly a slightly different threshold, so compare run with run, not with those.

HOW TO RUN. Put it in workspaces/Sycamore/lab/ next to run_tail_screen_sycamore.py AND run_retune_w4_sycamore.py (both are imported), and run
in the foreground, one command at a time:
    python -u run_kchan_w4_sycamore.py --dry-run                     # checks only, builds nothing; lists the runs
    python -u run_kchan_w4_sycamore.py --build-only --limit 1        # builds the probe's inputs, checks them, does NOT run tRIBS (about 1-2 minutes)
    python -u run_kchan_w4_sycamore.py                               # the probe run, about 36 minutes (finished runs are skipped)
    python -u run_kchan_w4_sycamore.py --kchan 0.5,2                 # probe 0.5 (skipped if done) and then 2
    python -u run_kchan_w4_sycamore.py --summary-only --plot         # tables and figure from what is finished
Do NOT run the notebook (or anything that rebuilds data/model/) while it runs.
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_retune_w4_sycamore as R          # noqa: E402  (the re-tune script: its helpers read the timing; importing it builds nothing)
T = R.T                                      # the tested W4 script that builds, runs and scores

# --------------------------------------------------------------------------------------------
# Names
# --------------------------------------------------------------------------------------------
FIG_DIR = "kchan_w4_figs"
SUMMARY_NAME = "kchan_w4_summary.csv"
WORK_KS = 1.03                       # the working Ks factor on the baseline Ks 0.187
DEFAULT_KCHAN = [0.5]                # mm/hr
ZERO_THR = 0.001                     # m3/s: an hour whose mean outlet flow is under this counts as zero-flow
RULE_A_FACTOR = 2.0                  # Oct-Dec volume must fall at least this many times as much (percent) as the Feb 27-28 volume
RULE_B_MAX_ZERO_PCT = 50.0
RULE_C_MIN_FRACTION = 0.5
W4_START, W4_END = "2016-10-01 00:00", "2017-04-01 00:00"
OCTDEC = ("2016-10-01 00:00", "2017-01-01 00:00")
NAN = float("nan")


def describe(job: dict) -> str:
    """The text the W4 script prints for a run (it looks this function up by name, so it is replaced here)."""
    if job["step"] == "base":
        return "baseline (no change)"
    ch = job["changes"]
    parts = [f"f x{job['f']:g}", f"As {job['an']:g}"]
    if "KS_MULT" in ch:
        parts.append(f"Ks x{ch['KS_MULT'] / T.BASE['KS_MULT']:.2f} (KS_MULT {ch['KS_MULT']:g})")
    if ch.get("OPTPERCOLATION") == 1:
        parts.append(f"channel loss ON, Kchan {ch['CHANNELCONDUCTIVITY']:g} mm/hr")
    else:
        parts.append("channel loss off")
    return ", ".join(parts)


def configure() -> None:
    T.SUMMARY_CSV = SUMMARY_NAME
    T.describe = describe


configure()


# --------------------------------------------------------------------------------------------
# The runs
# --------------------------------------------------------------------------------------------
def work_job() -> dict:
    """The loss-off working run: the same job (and run id) the re-tune script made, so a finished one is recognized and skipped."""
    base = T.an_job(0.15, 300)
    return R._with_tag(base, f"ksx{T.vtag(WORK_KS)}", {"KS_MULT": round(T.BASE["KS_MULT"] * WORK_KS, 6)}, "ks")


def kchan_job(k: float) -> dict:
    w = work_job()
    j = dict(w)
    j["changes"] = {**w["changes"], "OPTPERCOLATION": 1, "CHANNELCONDUCTIVITY": float(k)}
    j["run_id"] = w["run_id"][:-len("_nocl")] + f"_kchan{T.vtag(k)}"
    j["step"] = "kchan"
    j["kchan"] = float(k)
    return j


PLAN: dict = {}                     # filled by set_plan: base, work, kchan jobs


def set_plan(kchans: list) -> None:
    PLAN.clear()
    PLAN["base"] = T.make_base_job()
    PLAN["work"] = work_job()
    PLAN["kchan"] = [kchan_job(k) for k in sorted(set(kchans))]


def planned() -> list:
    return [PLAN["base"], PLAN["work"]] + PLAN["kchan"]


def show_order() -> dict:
    return {j["run_id"]: i for i, j in enumerate(planned())}


# --------------------------------------------------------------------------------------------
# Extra numbers computed from each finished run's outlet file (kept in the manifest once)
# --------------------------------------------------------------------------------------------
def zero_pct(q: pd.Series, a, b) -> float:
    """Percent of the hours in [a, b) whose mean outlet flow is under ZERO_THR."""
    h = q[(q.index >= pd.Timestamp(a)) & (q.index < pd.Timestamp(b))].resample("60min").mean().dropna()
    return float(100.0 * (h < ZERO_THR).mean()) if len(h) else NAN


PRE_RISE_FRACTION = 0.10             # the gauge's main rise is taken to begin when the gauge first reaches this fraction of its own peak


def pre_rise_max(series: pd.Series, a, b, t_cut):
    """Highest flow of `series` in [a, t_cut): everything BEFORE the gauge's main rise has begun. NaN if t_cut is missing or nothing lies before it."""
    if pd.isna(t_cut):
        return NAN
    s = R._cut(series, a, b).dropna()
    before = s[s.index < t_cut]
    return float(before.max()) if len(before) else NAN


def extra_metrics(q: pd.Series, gauge_q: pd.Series) -> dict:
    out = {"zero_pct_w4": zero_pct(q, W4_START, W4_END), "zero_pct_octdec": zero_pct(q, *OCTDEC)}
    for tag, (a, b) in R.RISE_WINDOWS.items():
        gw = R._cut(gauge_q, a, b).dropna()
        if gw.empty:
            continue
        pg = float(gw.max())
        t_cut = R.cross(gw, PRE_RISE_FRACTION * pg)               # when the gauge's main rise begins (the same cut for model and gauge)
        out[f"{tag}_gauge_peak_m3s"] = pg
        out[f"{tag}_gauge_rise_start"] = t_cut.strftime("%d %H:%M") if pd.notna(t_cut) else ""
        out[f"{tag}_pre_rise_max_model"] = pre_rise_max(q, a, b, t_cut)
        out[f"{tag}_pre_rise_max_gauge"] = pre_rise_max(gauge_q, a, b, t_cut)
    return out


def ensure_extra(man: dict, lab: Path, rid: str, gauge_q) -> bool:
    rec = man["runs"].get(rid) or {}
    if rec.get("status") != "done" or rec.get("kx") or gauge_q is None:
        return False
    qf = lab / "results" / f"{rid}_Outlet.qout"
    if not qf.is_file():
        return False
    try:
        rec["kx"] = extra_metrics(T.read_qout(qf, T.RUN_START_TS), gauge_q)
        return True
    except Exception as exc:                                        # the extra numbers never fail a run
        rec["kx_error"] = f"{type(exc).__name__}: {exc}"
        return False


# --------------------------------------------------------------------------------------------
# The tables
# --------------------------------------------------------------------------------------------
def label_of(rec: dict) -> str:
    if rec.get("step") == "base":
        return "baseline (As 1, loss off)"
    s = rec.get("settings") or {}
    ks = s.get("KS_MULT")
    txt = f"Ks x{ks / T.BASE['KS_MULT']:.2f}" if ks is not None else "Ks ?"
    if s.get("OPTPERCOLATION") == 1:
        return f"{txt} Kchan {s.get('CHANNELCONDUCTIVITY'):g}"
    return f"{txt} loss off"


def pct_change(new, old):
    return 100.0 * (new / old - 1.0) if (new is not None and old is not None and math.isfinite(new) and math.isfinite(old) and old != 0) else NAN


def sorted_runs(man: dict) -> list:
    order = show_order()
    return sorted([(rid, rec) for rid, rec in man["runs"].items() if rid in order], key=lambda t: order[t[0]])


def scores_of(man, lab, rid, csv_scores):
    return R.scores_row(man, lab, rid, csv_scores)


def rules(work_row: dict, work_tail: dict, k_row: dict, k_kx: dict) -> dict:
    """Rules A, B, C for one probe run against the loss-off working run. Each value: (state True/False/None, text)."""
    o_w, o_k = T.fnum(work_row, "octdec__vol_model_Mm3"), T.fnum(k_row, "octdec__vol_model_Mm3")
    f_w, f_k = T.fnum(work_row, "feb27_28__vol_ratio"), T.fnum(k_row, "feb27_28__vol_ratio")
    o_fall, f_fall = -pct_change(o_k, o_w), -pct_change(f_k, f_w)
    res = {}
    if math.isfinite(o_fall) and math.isfinite(f_fall):
        res["A"] = (bool(o_fall >= RULE_A_FACTOR * f_fall and o_fall > 0),
                    f"Oct-Dec volume falls {o_fall:.1f}% ({o_w:.2f} -> {o_k:.2f} M m3); Feb 27-28 volume falls {f_fall:.1f}% ({f_w:.3f} -> {f_k:.3f}); needed: Oct-Dec fall at least {RULE_A_FACTOR:g} x the Feb fall")
    else:
        res["A"] = (None, "not computed (a volume is missing)")
    z = (k_kx or {}).get("zero_pct_w4", NAN)
    res["B"] = ((z < RULE_B_MAX_ZERO_PCT) if math.isfinite(z) else None, f"zero-flow hours {z:.1f}% of W4 (limit: under {RULE_B_MAX_ZERO_PCT:g}%)" if math.isfinite(z) else "not computed")
    t_w, t_k = T.fnum(work_row, "tail_jan26_feb2__tail_vol_ratio"), T.fnum(k_row, "tail_jan26_feb2__tail_vol_ratio")
    if math.isfinite(t_w) and math.isfinite(t_k):
        res["C"] = (bool(t_k >= RULE_C_MIN_FRACTION * t_w), f"tail volume Jan 26 - Feb 2: {t_k:.2f} against {t_w:.2f} loss-off (needed: at least {RULE_C_MIN_FRACTION * t_w:.3f})")
    else:
        res["C"] = (None, "not computed")
    return res


def yn3(v):
    return "n.a." if v is None else ("PASS" if v else "MISS")


def print_report(man: dict, lab: Path, args) -> None:
    csv_scores = T.load_scores(lab)
    recs = sorted_runs(man)
    if not recs:
        T.log("No runs yet.")
        return
    gauge_q = None
    try:
        gauge_q = T.read_obs(lab / T.INIT_DATA / T.OBS_FILE)
    except Exception as exc:
        T.log(f"  (gauge file could not be read, so timing numbers are skipped: {exc})")
    new = False
    for rid, rec in recs:
        if rec.get("status") == "done":
            if not rec.get("tail") and (lab / "results" / f"{rid}_Outlet.qout").is_file():
                try:
                    rec["tail"] = T.compute_tail(lab, rid, args)
                    new = True
                except Exception:
                    pass
            new = R.ensure_rise(man, lab, rid, gauge_q) or new
            new = ensure_extra(man, lab, rid, gauge_q) or new
    if new:
        T.save_manifest(man)
    done = [(rid, rec, scores_of(man, lab, rid, csv_scores)) for rid, rec in recs if rec.get("status") == "done"]
    wid = PLAN["work"]["run_id"]
    work = next(((rec, row) for rid, rec, row in done if rid == wid), None)
    T.log()
    T.log("KCHAN PROBE on W4 (Oct 1 2016 - Apr 1 2017), working setting Ks x1.03 / f x0.15 / As 300 / k10. 'loss off' = the working run; 'Kchan x' = channel loss ON at x mm/hr.")
    T.log("  vol = model volume / gauge volume. pk = model peak vs gauge peak, %. arr = the scorer's arrival error (20% rule), hours; negative = early.")
    if not done:
        T.log("  (no scored runs yet)")
        return

    T.log()
    T.log("  STORMS")
    T.log("    " + f"{'run':<26}| {'fw vol':>6} {'fw pk%':>7} {'fw arr':>6} | {'f27 vol':>7} {'f27 pk%':>8} {'f27 arr':>7} | {'f18 vol':>7} {'2nd Jan vol':>11} | {'storms_ok':>9} {'peaks_ok':>8}")
    for rid, rec, row in done:
        T.log("    " + f"{label_of(rec):<26}| {T.fm(T.fnum(row, 'flood_week__vol_ratio'), 3):>6} {T.fm(T.fnum(row, 'flood_week__peak_err_pct'), 1):>7} "
                       f"{T.fm(T.fnum(row, 'flood_week__arrival_err_h')):>6} | {T.fm(T.fnum(row, 'feb27_28__vol_ratio'), 3):>7} "
                       f"{T.fm(T.fnum(row, 'feb27_28__peak_err_pct'), 1):>8} {T.fm(T.fnum(row, 'feb27_28__arrival_err_h')):>7} | "
                       f"{T.fm(T.fnum(row, 'feb18_19__vol_ratio'), 2):>7} {T.fm(T.fnum(row, 'jan_second_storm__vol_ratio'), 2):>11} | "
                       f"{T.yn(T.storms_ok(row)):>9} {T.yn(R.peaks_ok(row)):>8}")
    if work:
        T.log("    change from the loss-off working run, % of volume:")
        for rid, rec, row in done:
            if rid == wid or rec.get("step") == "base":
                continue
            ch = [("flood week", "flood_week__vol_ratio"), ("Feb 27-28", "feb27_28__vol_ratio"), ("Feb 18-19", "feb18_19__vol_ratio"),
                  ("2nd Jan", "jan_second_storm__vol_ratio"), ("Oct-Dec", "octdec__vol_model_Mm3"), ("tail 1", "tail_jan26_feb2__tail_vol_ratio")]
            T.log("      " + f"{label_of(rec):<26}" + "  ".join(f"{n} {T.fm(pct_change(T.fnum(row, c), T.fnum(work[1], c)), 1, True)}%" for n, c in ch))

    if work:
        probes = [(rid, rec, row) for rid, rec, row in done if rec.get("step") == "kchan"]
        if probes:
            T.log()
            T.log("  RULES (written before the run; PASS / MISS only prints; you decide)")
            for rid, rec, row in probes:
                res = rules(work[1], work[0].get("tail") or {}, row, rec.get("kx") or {})
                T.log(f"    {label_of(rec)}")
                for key in ("A", "B", "C"):
                    st, txt = res[key]
                    T.log(f"      RULE {key}  {yn3(st):>4}  {txt}")
                kx_w, kx_k = work[0].get("kx") or {}, rec.get("kx") or {}
                T.log(f"      SHOULDER (highest Feb 27-28 flow before the gauge's main rise, cut at {kx_k.get('f27_gauge_rise_start', '?')}): loss off {T.fm(kx_w.get('f27_pre_rise_max_model'), 1)} m3/s, "
                      f"this run {T.fm(kx_k.get('f27_pre_rise_max_model'), 1)} m3/s, gauge {T.fm(kx_k.get('f27_pre_rise_max_gauge'), 1)} m3/s  (you judge 'clearly toward the gauge')")
                all3 = [res[k][0] for k in ("A", "B", "C")]
                T.log(f"      A, B, C together: " + ("all PASS" if all(v is True for v in all3) else ("at least one MISS" if any(v is False for v in all3) else "not all computed")))
    else:
        T.log("  (the loss-off working run is not finished, so the rules and changes cannot be computed)")

    T.log()
    T.log(f"  FALSE FLOW and ZERO-FLOW  (Oct 1 - Dec 31; the gauge is almost dry. threshold for 'days above' 0.13 m3/s; zero-flow = hourly mean under {ZERO_THR:g} m3/s)")
    T.log("    " + f"{'run':<26}| {'Oct-Dec M m3':>12} {'gauge':>7} {'x gauge':>8} {'goal*':>6} | {'days >0.13':>10} {'max m3/s':>9} | {'zero-flow % W4':>14} {'% Oct-Dec':>10} | {'chan loss M m3':>14}")
    work_octdec = T.fnum(work[1], "octdec__vol_model_Mm3") if work else NAN
    for rid, rec, row in done:
        kx = rec.get("kx") or {}
        v, gv = T.fnum(row, "octdec__vol_model_Mm3"), T.fnum(row, "octdec__vol_gauge_Mm3")
        T.log("    " + f"{label_of(rec):<26}| {T.fm(v, 2):>12} {T.fm(gv, 3):>7} {T.fm(v / gv if gv else NAN, 1):>8} {T.fm(0.5 * work_octdec, 2):>6} | "
                       f"{T.fm(T.fnum(row, 'octdec__days_above_model'), 0):>10} {T.fm(T.fnum(row, 'octdec__max_model_m3s'), 2):>9} | "
                       f"{T.fm(kx.get('zero_pct_w4'), 1):>14} {T.fm(kx.get('zero_pct_octdec'), 1):>10} | {T.fm(T.fnum(row, 'proc__channel_loss_Mm3'), 3):>14}")
    T.log("    (* goal = half of the loss-off working run's Oct-Dec volume; the gauge's own Oct-Dec maximum is " + T.fm(T.fnum(done[0][2], 'octdec__max_gauge_m3s'), 3) + " m3/s)")

    rises = [(rid, rec) for rid, rec, _ in done if rec.get("rise")]
    if rises:
        T.log()
        T.log("  TIMING  (hours, model minus gauge; negative = early. t20 = the scorer's rule, t50 = the proposed arrival rule, t80; 'before' = highest flow, m3/s, before the gauge first reaches 10% of its own peak)")
        T.log("    " + f"{'run':<26}| {'FLOOD WEEK t20':>14} {'t50':>6} {'t80':>6} {'before: model':>13} {'gauge':>6} | {'FEB 27-28 t20':>13} {'t50':>6} {'t80':>6} {'before: model':>13} {'gauge':>6}")
        for rid, rec in rises:
            r, kx = rec["rise"], rec.get("kx") or {}
            T.log("    " + f"{label_of(rec):<26}| {T.fm(r.get('fw_t20_h'), 2, True):>14} {T.fm(r.get('fw_t50_h'), 2, True):>6} {T.fm(r.get('fw_t80_h'), 2, True):>6} "
                           f"{T.fm(kx.get('fw_pre_rise_max_model'), 1):>13} {T.fm(kx.get('fw_pre_rise_max_gauge'), 1):>6} | {T.fm(r.get('f27_t20_h'), 2, True):>13} "
                           f"{T.fm(r.get('f27_t50_h'), 2, True):>6} {T.fm(r.get('f27_t80_h'), 2, True):>6} {T.fm(kx.get('f27_pre_rise_max_model'), 1):>13} {T.fm(kx.get('f27_pre_rise_max_gauge'), 1):>6}")

    first_tail = next((rec["tail"] for _, rec, _ in done if rec.get("tail")), None)
    T.log()
    if first_tail:
        T.log(f"  TAIL  (mean flow in 24 h blocks, m3/s, {T.TAIL_OFFSETS_D[0]}, {T.TAIL_OFFSETS_D[1]}, {T.TAIL_OFFSETS_D[2]} days after the gauge flood peak ({first_tail.get('gauge_peak_time', '')}); "
              f"band {args.band_lo:g}-{args.band_hi:g} m3/s; days in band counts the {first_tail.get('days_total', '?')} blocks from {T.TAIL_START[:10]} to {T.TAIL_END[:10]})")
        T.log("    " + f"{'run':<26}| " + " ".join(f"{'d' + str(d):>7}" for d in T.TAIL_OFFSETS_D) + f" | {'days in band':>12} | {'tail1 vol':>9} {'tail2 vol':>9}")
        T.log("    " + f"{'gauge':<26}| " + " ".join(f"{T.fm(first_tail.get(f'gauge_d{d}_m3s'), 3):>7}" for d in T.TAIL_OFFSETS_D)
              + f" | {str(first_tail.get('gauge_days_in_band', '?')) + '/' + str(first_tail.get('days_total', '?')):>12} |")
        for rid, rec, row in done:
            tail = rec.get("tail")
            if tail:
                T.log("    " + f"{label_of(rec):<26}| " + " ".join(f"{T.fm(tail.get(f'model_d{d}_m3s'), 3):>7}" for d in T.TAIL_OFFSETS_D)
                      + f" | {str(tail.get('model_days_in_band', '?')) + '/' + str(tail.get('days_total', '?')):>12} | {T.fm(T.fnum(row, 'tail_jan26_feb2__tail_vol_ratio')):>9} {T.fm(T.fnum(row, 'tail_mar3_mar11__tail_vol_ratio')):>9}")
    else:
        T.log("  TAIL  (no tail numbers yet)")

    T.log()
    T.log("  PROCESS  (water table = depth below the surface, mm)")
    T.log("    " + f"{'run':<26}| {'water table mm: first':>21} {'shallowest':>10} {'deepest':>8} | {'snow max cm':>11} | {'build min':>9} {'tRIBS min':>9}")
    for rid, rec, row in done:
        T.log("    " + f"{label_of(rec):<26}| {T.fm(T.fnum(row, 'proc__mdgw_first_mm'), 0):>21} {T.fm(T.fnum(row, 'proc__mdgw_min_mm'), 0):>10} {T.fm(T.fnum(row, 'proc__mdgw_max_mm'), 0):>8} | "
                       f"{T.fm(T.fnum(row, 'proc__snow_max_cm'), 1):>11} | {str(rec.get('build_min', '')):>9} {str(rec.get('tribs_min', '')):>9}")
    failed = [rid for rid, rec in man["runs"].items() if rec.get("status") == "failed" and rid in show_order()]
    if failed:
        T.log()
        T.log("  FAILED runs (see the log): " + ", ".join(failed))
    write_summary_csv(man, lab, csv_scores)
    T.log()
    T.log(f"  Wrote {SUMMARY_NAME}. These tables decide nothing: whether Kchan is worth keeping, and the next value, is your call.")


def write_summary_csv(man: dict, lab: Path, csv_scores: dict) -> None:
    rows = []
    for rid, rec in sorted_runs(man):
        row_sc = scores_of(man, lab, rid, csv_scores)
        s = rec.get("settings") or {}
        r = dict(run_id=rid, step=rec.get("step"), status=rec.get("status"), F_MULT=s.get("F_MULT"), ANISOTROPY=s.get("ANISOTROPY"),
                 KS_MULT=s.get("KS_MULT"), OPTPERCOLATION=s.get("OPTPERCOLATION"), CHANNELCONDUCTIVITY=s.get("CHANNELCONDUCTIVITY"),
                 storms_ok=(T.storms_ok(row_sc) if rec.get("status") == "done" else ""), peaks_ok=(R.peaks_ok(row_sc) if rec.get("status") == "done" else ""))
        for c in T.SCORE_COLS:
            r[c] = T.fnum(row_sc, c) if row_sc else NAN
        r.update({k: v for k, v in (rec.get("tail") or {}).items() if k not in ("band_lo", "band_hi")})
        r.update({f"rise_{k}": v for k, v in (rec.get("rise") or {}).items()})
        r.update({f"kx_{k}": v for k, v in (rec.get("kx") or {}).items()})
        r.update(build_min=rec.get("build_min", ""), tribs_min=rec.get("tribs_min", ""), output_mb=rec.get("output_mb", ""))
        rows.append(r)
    cols = []
    for r in rows:
        for c in r:
            if c not in cols:
                cols.append(c)
    with open(lab / SUMMARY_NAME, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: ("" if (isinstance(r.get(c), float) and not math.isfinite(r[c])) else r.get(c, "")) for c in cols})


# --------------------------------------------------------------------------------------------
# The figure
# --------------------------------------------------------------------------------------------
def make_plot(man: dict, lab: Path, out: Path) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    style = [("#eb6834", "-"), ("#2a78d6", "-"), ("#1baf7a", "--"), ("#7b4fc4", ":")]
    runs = []
    for j in [PLAN["work"]] + PLAN["kchan"]:
        rec = man["runs"].get(j["run_id"]) or {}
        qf = lab / "results" / f"{j['run_id']}_Outlet.qout"
        if rec.get("status") == "done" and qf.is_file():
            runs.append((label_of(rec), T.read_qout(qf, T.RUN_START_TS)))
    if not runs:
        return "no finished run with an outlet file; no figure made"
    gauge = T.read_obs(lab / T.INIT_DATA / T.OBS_FILE)
    panels = [("A. Flood week", "2017-01-18", "2017-01-25", False, None), ("B. Feb 27-28 storm", "2017-02-27 12:00", "2017-03-01 12:00", False, None),
              ("C. Oct-Dec (log scale; the false flow)", "2016-10-01", "2017-01-01", True, (0.001, 30)),
              ("D. The tail (log scale; shaded = the gauge's 0.3-1 m3/s band)", "2017-01-26", "2017-02-18", True, (0.001, 30))]
    fig, axes = plt.subplots(2, 2, figsize=(14, 8.5), gridspec_kw={"wspace": 0.2, "hspace": 0.4})
    for ax, (title, a, b, logy, ylim) in zip(axes.ravel(), panels):
        g = R._cut(gauge, a, b)
        ax.plot(g.index, g.values, color="black", lw=2.0, label="gauge")
        for k, (lbl, q) in enumerate(runs):
            color, ls = style[k % len(style)]
            s = R._cut(q, a, b)
            ax.plot(s.index, s.values, color=color, ls=ls, lw=1.5, label=lbl)
        if logy:
            ax.set_yscale("log")
            ax.set_ylim(*ylim)
            if title.startswith("D"):
                ax.axhspan(0.3, 1.0, color="#8a8984", alpha=0.15, lw=0)
        ax.set_title(title, loc="left", fontsize=11)
        ax.set_ylabel("Outlet flow (m$^3$/s)")
        ax.grid(alpha=0.25)
        ax.tick_params(axis="x", labelrotation=30)
        for s_ in ("top", "right"):
            ax.spines[s_].set_visible(False)
    axes[0, 0].legend(frameon=False, fontsize=8, loc="upper right")
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return f"wrote {out.name} with {len(runs)} run(s) and the gauge"


# --------------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------------
def parse_kchan(text: str) -> list:
    vals = []
    for x in text.split(","):
        x = x.strip()
        if not x:
            continue
        v = float(x)
        if not (0.001 <= v <= 1000):
            raise ValueError(f"Kchan {v:g} mm/hr is outside 0.001-1000")
        vals.append(v)
    if not vals:
        raise ValueError("no Kchan value given")
    return vals


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description="Kchan probe on the working setting on W4 (see the top of this file).")
    ap.add_argument("--kchan", default=",".join(f"{k:g}" for k in DEFAULT_KCHAN), help="Kchan value(s) in mm/hr, comma-separated (default %(default)s). Finished runs are skipped")
    ap.add_argument("--dry-run", action="store_true", help="check everything, list the plan, build nothing")
    ap.add_argument("--summary-only", action="store_true", help="print the tables from the runs already scored, then stop")
    ap.add_argument("--build-only", action="store_true", help="build and check the inputs but do not run tRIBS (use with --limit 1)")
    ap.add_argument("--limit", type=int, default=None, help="make only the first N runs that still need work (finished runs do not count)")
    ap.add_argument("--plot", action="store_true", help=f"make {FIG_DIR}/fig_kchan_w4.png from the finished runs")
    ap.add_argument("--band-lo", type=float, default=T.BAND_LO, help="slow-flow band, low end, m3/s")
    ap.add_argument("--band-hi", type=float, default=T.BAND_HI, help="slow-flow band, high end, m3/s")
    ap.add_argument("--redo", default=None, help="comma-separated run ids to run again even if done (you almost never want this)")
    ap.add_argument("--timeout-min", type=float, default=None, help="stop a tRIBS run after this many minutes (default: 3.5 times the expected time, about 2 hours for W4)")
    ap.add_argument("--cell-timeout", type=int, default=1800, help="seconds allowed for one notebook cell")
    ap.add_argument("--min-free-gb", type=float, default=3.0, help="stop if less than this much disk is free (default 3)")
    ap.add_argument("--heartbeat-min", type=float, default=T.HEARTBEAT_MIN, help="print a progress line every this many minutes while tRIBS runs (default 5)")
    ap.add_argument("--notebook", default=None, help="notebook to run (default: Make_Sycamore_Model.ipynb, else Make_Sycamore_Model_clean.ipynb)")
    ap.add_argument("--scorer", default=T.SCORER)
    ap.add_argument("--tribs", default=os.environ.get("TRIBS_BIN") or (T.TRIBS_DEFAULT if os.path.isfile(T.TRIBS_DEFAULT) else shutil.which("tRIBS") or T.TRIBS_DEFAULT))
    ap.add_argument("--skip-kernel-check", action="store_true", help="skip the pytRIBS import test in --dry-run")
    a = ap.parse_args(argv)
    a.steps = "kchan"
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
    set_plan(kchans)
    lab = Path(__file__).resolve().parent
    os.chdir(lab)
    (lab / T.LOG_DIR).mkdir(exist_ok=True)
    T._LOGFILE = str(lab / T.LOG_DIR / "batch_log.txt")
    man = T.load_manifest()

    def finish_outputs() -> None:
        if args.plot:
            try:
                T.log("  Figure: " + make_plot(man, lab, lab / FIG_DIR / "fig_kchan_w4.png"))
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
    per_lo, per_hi = T.WIN["hours"] * T.MIN_PER_1000H / 1000.0, T.WIN["hours"] * T.MIN_PER_1000H_SLOW / 1000.0
    T.log(f"Working folder: {lab}")
    T.log(f"Notebook: {args.notebook}")
    T.log(f"Plan: Kchan probe(s) {', '.join(f'{k:g}' for k in kchans)} mm/hr on the working setting, window W4 ({T.WIN['hours']} h from {T.WIN['start']}, water table {T.WIN['gw']:g} x bedrock): "
          f"{n_todo} run(s) to make" + (f", {n_done} already finished" if n_done else "") + ".")
    for j in preview:
        what = "already done, skipped" if T.is_finished(man, j, args) else "NEW RUN"
        T.log(f"  {j['run_id']:56s} {describe(j):70s} {what}")
    T.log("  (the baseline and the loss-off working run (Ks x1.03) are the earlier scripts' runs; they are skipped when the shared manifest says they are finished)")
    T.log(f"  tRIBS time alone: about {n_todo * per_lo / 60:.1f}-{n_todo * per_hi / 60:.1f} hours ({per_lo:.0f}-{per_hi:.0f} min per run, a loss-on run may be a little slower); each run also needs a notebook build of about 1 min.")

    results = T.preflight(args, lab, preview, nb_note, man)
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
    man["kchan_w4_last_started"] = T.now_iso()
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
        code = run_batch([PLAN["base"]])                         # the baseline always comes first: every other run's .in file is checked against it
        if code != 0:
            T.log("The remaining runs were not started because the baseline failed or was interrupted. Fix that, then run the same command again.")
        else:
            code = run_batch([PLAN["work"]])                     # the loss-off working run is normally finished and skipped; if not, it is made first
            if code != 0:
                T.log("The probe was not started because the loss-off working run failed or was interrupted. Fix that, then run the same command again.")
            else:
                new_runs = [describe(j) for j in PLAN["kchan"] if not T.is_finished(man, j, args)]
                if new_runs:
                    T.log()
                    T.log("probe: " + "; ".join(new_runs))
                code = run_batch(PLAN["kchan"]) or code
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
    T.log(f"{'run':56s} {'status':16s} {'build':>7s} {'tRIBS':>7s}")
    for rid, rec in man["runs"].items():
        if rid in show_order():
            T.log(f"{rid:56s} {rec.get('status', '?'):16s} {str(rec.get('build_min', '')):>7s} {str(rec.get('tribs_min', '')):>7s}")
    print_report(man, lab, args)
    finish_outputs()
    return code


if __name__ == "__main__":
    sys.exit(main())
