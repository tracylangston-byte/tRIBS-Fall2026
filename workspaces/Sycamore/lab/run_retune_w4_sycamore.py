#!/usr/bin/env python3
"""
run_retune_w4_sycamore.py  --  The small Ks re-tune at f x0.15 / As 300, on the FULL tail-screen window W4 (Oct 1 2016 - Apr 1 2017).

WHY W4. The W2 trial (run_retune_w2_sycamore.py, 2026-11-02) showed that the short window does not carry the slow flow the way W4 does:
f x0.15 / As 300 gave 16 of 23 days in the slow-flow band on W4 but 5 of 23 on W2, and its Feb 27-28 numbers moved by a different amount
than the baseline's. So the Ks re-tune is done on W4, where the earlier results were made. A W4 run takes about 34 minutes.

WHAT IT DOES. It makes three new runs that differ from the W4 run "f x0.15 / As 300" only in Ks:
    Ks x1.06  (KS_MULT 0.19822)      Ks x1.03  (KS_MULT 0.19261)      Ks x1.10  (KS_MULT 0.2057)
(factors on the baseline Ks 0.187; Au = 1; everything else as in the W4 run). The baseline run and the "f x0.15 / As 300" run
already exist from the tail screen and are REUSED, not repeated. The script shares the tail screen's manifest, log folder, score
file and lock, so it knows which runs are finished (the same command run again resumes, and a finished run is never repeated).
Its own report goes to retune_w4_summary.csv (the tail screen's tail_screen_w4_summary.csv is not touched).

WHY THESE SIZES. A regression on earlier runs at As 1 (NOT tested at As 300) says the Feb 27-28 volume goes as Ks^-1.9 and the flood
week as Ks^-0.9. From the W4 run f x0.15 / As 300 (flood-week volume 1.078, Feb 27-28 volume 1.065, Feb 27-28 peak +15.5%):
    Ks x1.03 -> Feb 27-28 volume about 1.01, flood week about 1.05
    Ks x1.06 -> Feb 27-28 volume about 0.95, flood week about 1.02
    Ks x1.10 -> Feb 27-28 volume about 0.89, flood week about 0.99
The run results replace this prior: the script prints the exponents it actually finds beside the prior. Ks x1.10 is there to show
where it goes too far (your anchor rule keeps Feb 27-28 volume between 0.90 and 1.10).

WHAT IT PRINTS (all read from files; it decides nothing):
    STORMS   flood week, 2nd January storm, Feb 18-19, Feb 27-28: volume ratio, peak error %, scorer arrival error.
             storms_ok = your anchor rule (flood-week volume 0.90-1.12 and Feb 27-28 volume 0.90-1.10); peaks_ok = both peak errors within +-15%.
    TAIL     mean flow in 24 h blocks 7, 14, 21 days after the gauge flood peak; days in the 0.3-1 m3/s band; tail volume ratios.
             (The question for the tail: does the Ks increase cost the slow flow? On W4, f x0.15 / As 300 gave 16 of 23 days.)
    TIMING   rising limb of the flood week and of Feb 27-28 (20, 50, 80% of the gauge peak, and the peak time) and the early Feb 28 shoulder.
    Ks RESPONSE   the runs side by side with change from the Ks x1.00 run and the exponents found.
    PROCESS  infiltration-excess share, water table depth, depth-decay number.
    --plot makes retune_w4_figs/fig_w4_overlays.png (flood week, Feb 27-28, tail).

WHAT IT DOES NOT DO. It does not re-tune PsiB or routing (if Ks cannot fix the Feb 27-28 peak, those are the next, smaller steps), and
it does not test the Oct - Dec false flow, the Jan 23-24 storm, rain scaling or snow.

Put it in workspaces/Sycamore/lab/ next to run_tail_screen_sycamore.py (both are needed) and run it in the foreground:
    python -u run_retune_w4_sycamore.py --dry-run                    # checks only, builds nothing; lists the runs that are new
    python -u run_retune_w4_sycamore.py --build-only --limit 1       # build the first new run and check it, no tRIBS
    python -u run_retune_w4_sycamore.py --limit 1                    # first new run (Ks x1.06), about 35 minutes
    python -u run_retune_w4_sycamore.py                              # all three new runs (about 1.7 hours); finished runs are skipped
    python -u run_retune_w4_sycamore.py --summary-only --plot        # tables and figure from what is finished
Do NOT run the notebook (or anything that rebuilds data/model/) while it runs.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_tail_screen_sycamore as T          # noqa: E402  (the tested W4 script; its functions do the building, running and scoring)

# --------------------------------------------------------------------------------------------
# Names. The window, run names, manifest, logs and score file are the tail screen's own (W4 is the default there).
# --------------------------------------------------------------------------------------------
PREFIX = T.RUN_PREFIX                           # W4_d012_k10 (the scorer reads the window W4 from the name)
FIG_DIR = "retune_w4_figs"
SUMMARY_NAME = "retune_w4_summary.csv"
SCRIPT_NAME = "run_retune_w4_sycamore.py"
STEP_ORDER = ["ks"]
DEFAULT_STEPS = "ks"


def describe(job: dict) -> str:
    """The text the W4 script prints for a run (it looks this function up by name, so it is replaced here)."""
    if job["step"] == "base":
        return "baseline (no change)"
    ch = job["changes"]
    parts = [f"f x{job['f']:g}", f"As {job['an']:g}"]
    if "KS_MULT" in ch:
        parts.append(f"Ks x{ch['KS_MULT'] / T.BASE['KS_MULT']:.2f} (KS_MULT {ch['KS_MULT']:g})")
    return ", ".join(parts)


def configure() -> None:
    T.SUMMARY_CSV = SUMMARY_NAME
    T.describe = describe


configure()

NEEDED_FROM_W4_SCRIPT = ["do_run", "preflight", "make_base_job", "an_job", "au_job", "tail_metrics", "read_qout", "read_obs", "compute_tail", "load_scores",
                         "expected_files", "settings_record", "settings_for", "is_finished", "take_lock", "save_manifest", "load_manifest", "pick_notebook",
                         "install_term_handler", "fnum", "fm", "yn", "sci", "storms_ok", "log", "now_iso", "fmt_min", "vtag", "au_of", "SCORE_COLS",
                         "ANCHOR_FW", "ANCHOR_F27", "TAIL_START", "TAIL_END", "TAIL_OFFSETS_D", "BAND_LO", "BAND_HI", "HEARTBEAT_MIN", "TRIBS_DEFAULT",
                         "SCORER", "INIT_DATA", "OBS_FILE", "MIN_PER_1000H", "MIN_PER_1000H_SLOW", "RunError", "Abort", "MANIFEST", "LOG_DIR", "WIN"]


def missing_from_w4_script() -> list:
    return [n for n in NEEDED_FROM_W4_SCRIPT if not hasattr(T, n)]


# --------------------------------------------------------------------------------------------
# The runs
# --------------------------------------------------------------------------------------------
KS_PLAN = [(0.15, 300, 1.06), (0.15, 300, 1.03), (0.15, 300, 1.10)]   # (f, As, Ks factor on the baseline Ks); order = order of running
KS_PRIOR = {"fw_vol": -0.9, "f27_vol": -1.9}      # exponents of volume against Ks from the regression at As 1 (150 W1 and 19 W2 runs)


def job_tag(job: dict) -> str:
    """'W4_d012_k10_f0p15_as300_ksx1p06_nocl' -> 'f0p15_as300_ksx1p06'"""
    return job["run_id"][len(PREFIX) + 1:-len("_nocl")]


def _with_tag(job: dict, tag: str, changes: dict, step: str) -> dict:
    j = dict(job)
    j["changes"] = {**job["changes"], **changes}
    j["run_id"] = job["run_id"].replace("_nocl", f"_{tag}_nocl")
    j["step"] = step
    return j


def build_jobs() -> tuple:
    """Returns (jobs by tag, jobs of each step in the order they run, every job in the order the tables show them)."""
    J = {"base": T.make_base_job()}
    ref = T.an_job(0.15, 300)
    J["f0p15_as300"] = ref
    ks_jobs = []
    for fv, a, factor in KS_PLAN:
        j = _with_tag(T.an_job(fv, a), f"ksx{T.vtag(factor)}", {"KS_MULT": round(T.BASE["KS_MULT"] * factor, 6)}, "ks")
        j["ks_factor"] = factor
        ks_jobs.append(j)
    steps = {"ks": [J["f0p15_as300"]] + ks_jobs}
    show = [J["base"], J["f0p15_as300"]] + sorted(ks_jobs, key=lambda j: j["ks_factor"])
    for j in ks_jobs:
        J[job_tag(j)] = j
    return J, steps, show


JOBS, STEP_JOBS, SHOW_JOBS = build_jobs()
SHOW_ORDER = {j["run_id"]: i for i, j in enumerate(SHOW_JOBS)}


def planned(steps: list) -> list:
    """The baseline first, then the runs of each chosen step, no run twice (the reference runs are shared)."""
    out, seen = [JOBS["base"]], {JOBS["base"]["run_id"]}
    for st in STEP_ORDER:
        if st in steps:
            for j in STEP_JOBS[st]:
                if j["run_id"] not in seen:
                    seen.add(j["run_id"])
                    out.append(j)
    return out

# --------------------------------------------------------------------------------------------
# The numbers that are read from each finished run
# --------------------------------------------------------------------------------------------
NAN = float("nan")


def scores_row(man: dict, lab: Path, rid: str, csv_scores: dict | None = None) -> dict:
    """Score columns of one run: the manifest's copy, overlaid by the scorer's CSV row when it exists."""
    row = dict((man["runs"].get(rid) or {}).get("scores") or {})
    sc = (csv_scores if csv_scores is not None else T.load_scores(lab)).get(rid)
    if sc:
        row.update(sc)
    return row


def key_metrics(man: dict, lab: Path, rid: str, csv_scores: dict | None = None):
    """The six headline numbers of a finished run, or None if the run is not done or one number is missing."""
    rec = man["runs"].get(rid) or {}
    if rec.get("status") != "done":
        return None
    row = scores_row(man, lab, rid, csv_scores)
    tail = rec.get("tail") or {}
    days = tail.get("model_days_in_band")
    vals = dict(fw_vol=T.fnum(row, "flood_week__vol_ratio"), fw_pk=T.fnum(row, "flood_week__peak_err_pct"),
                f27_vol=T.fnum(row, "feb27_28__vol_ratio"), f27_pk=T.fnum(row, "feb27_28__peak_err_pct"),
                days=float(days) if days is not None else NAN, tail1=T.fnum(row, "tail_jan26_feb2__tail_vol_ratio"))
    return vals if all(math.isfinite(v) for v in vals.values()) else None


def peaks_ok(row: dict):
    a, b = T.fnum(row, "flood_week__peak_err_pct"), T.fnum(row, "feb27_28__peak_err_pct")
    return bool(abs(a) <= 15 and abs(b) <= 15) if (math.isfinite(a) and math.isfinite(b)) else None


def cross(s: pd.Series, level: float):
    """First time the series rises through `level` (linear between samples); NaT if it never does."""
    s = s.dropna()
    if s.empty:
        return pd.NaT
    above = s.to_numpy() >= level
    r = np.nonzero(above[1:] & ~above[:-1])[0] + 1
    if len(r) == 0:
        return pd.NaT
    i = int(r[0])
    tp, tn = s.index[i - 1], s.index[i]
    qp, qn = float(s.iloc[i - 1]), float(s.iloc[i])
    return tn if qn == qp else tp + (tn - tp) * ((level - qp) / (qn - qp))


def _cut(q: pd.Series, a, b) -> pd.Series:
    return q[(q.index >= pd.Timestamp(a)) & (q.index < pd.Timestamp(b))]


RISE_WINDOWS = {"fw": ("2017-01-18 00:00", "2017-01-25 00:00"), "f27": ("2017-02-27 00:00", "2017-03-03 00:00")}
SHOULDER = ("2017-02-28 00:00", "2017-02-28 09:30")


def hours(td) -> float:
    return float(td / pd.Timedelta(hours=1)) if pd.notna(td) else NAN


def rise_metrics(model_q: pd.Series, gauge_q: pd.Series) -> dict:
    """Timing of the rising limb vs the gauge, same windows as the scorer's events. Times are hours, model minus gauge (negative = early)."""
    out = {}
    for tag, (a, b) in RISE_WINDOWS.items():
        gw, mw = _cut(gauge_q, a, b).dropna(), _cut(model_q, a, b).dropna()
        if gw.empty or mw.empty:
            continue
        pg = float(gw.max())
        for p in (0.2, 0.5, 0.8):
            out[f"{tag}_t{int(p * 100)}_h"] = hours(cross(mw, p * pg) - cross(gw, p * pg))
        out[f"{tag}_pk_dt_h"] = hours(mw.idxmax() - gw.idxmax())
    sm, sg = _cut(model_q, *SHOULDER).dropna(), _cut(gauge_q, *SHOULDER).dropna()
    if not sm.empty and not sg.empty:
        out["shoulder_max_m3s"], out["gauge_shoulder_max_m3s"] = float(sm.max()), float(sg.max())
        for key, s in (("shoulder_first5", sm), ("gauge_first5", sg)):
            t = cross(s, 5.0)
            if float(s.iloc[0]) >= 5.0:
                out[key] = "before 00:00"                           # already above 5 m3/s when the window opens
            else:
                out[key] = t.strftime("%d %H:%M") if pd.notna(t) else "not by 09:30"
    return out


def ensure_rise(man: dict, lab: Path, rid: str, gauge_q) -> bool:
    """Compute the timing numbers of a finished run once and keep them in the manifest. Returns True if something new was stored."""
    rec = man["runs"].get(rid) or {}
    if rec.get("status") != "done" or rec.get("rise") or gauge_q is None:
        return False
    qf = lab / "results" / f"{rid}_Outlet.qout"
    if not qf.is_file():
        return False
    try:
        rec["rise"] = rise_metrics(T.read_qout(qf, T.RUN_START_TS), gauge_q)
        return True
    except Exception as exc:                                        # the timing table is extra; never fail over it
        rec["rise_error"] = f"{type(exc).__name__}: {exc}"
        return False




# --------------------------------------------------------------------------------------------
# The tables
# --------------------------------------------------------------------------------------------
def label_of(rec: dict) -> str:
    if rec.get("step") == "base":
        return "baseline (As 1)"
    s = rec.get("settings") or {}
    f, a = s.get("F_MULT", rec.get("f")), s.get("ANISOTROPY", rec.get("an"))
    parts = [f"f x{f:g}", f"As {a:g}"]
    ks = s.get("KS_MULT")
    if ks is not None and abs(ks - T.BASE["KS_MULT"]) > 1e-9:
        parts.append(f"Ks x{ks / T.BASE['KS_MULT']:.2f}")
    return " ".join(parts)


def sorted_runs(man: dict) -> list:
    recs = [(rid, rec) for rid, rec in man["runs"].items() if rid in SHOW_ORDER]
    return sorted(recs, key=lambda t: SHOW_ORDER.get(t[0], 10 ** 6))


def ks_fit(points: list, key: str):
    """Slope of ln(value) against ln(Ks factor) (an exponent) through (factor, value) points; None with fewer than two points."""
    pts = [(f, v) for f, v in points if v is not None and math.isfinite(v) and v > 0]
    if len(set(f for f, _ in pts)) < 2:
        return None
    return float(np.polyfit([math.log(f) for f, _ in pts], [math.log(v) for _, v in pts], 1)[0])


def ks_peak_slope(points: list):
    """Change in Feb 27-28 peak error (points) per +1% in Ks, from a straight line through (factor, peak error) points."""
    pts = [(f, v) for f, v in points if v is not None and math.isfinite(v)]
    if len(set(f for f, _ in pts)) < 2:
        return None
    return float(np.polyfit([100 * math.log(f) for f, _ in pts], [v for _, v in pts], 1)[0])


def print_ks_response(man: dict, lab: Path, csv_scores: dict) -> None:
    printed = False
    for fv, a in dict.fromkeys((f, a) for f, a, _ in KS_PLAN):
        ref_tag = f"f{T.vtag(fv)}_as{T.vtag(a)}"
        members = [(1.0, JOBS[ref_tag])] + [(j["ks_factor"], j) for j in STEP_JOBS["ks"][1:] if j["f"] == fv and j["an"] == a]
        members.sort(key=lambda t: t[0])
        got = [(fac, j, key_metrics(man, lab, j["run_id"], csv_scores)) for fac, j in members]
        got = [t for t in got if t[2] is not None]
        if len(got) < 2 or not any(fac != 1.0 for fac, _, _ in got):
            continue
        if not printed:
            T.log()
            T.log("  Ks RESPONSE  (the same f and As, only Ks changes; Ks x1.00 is the baseline Ks 0.187, i.e. the W4 run f x0.15 / As 300. d = change from the Ks x1.00 run)")
            printed = True
        T.log(f"    f x{fv:g}, As {a:g}")
        T.log("      " + f"{'Ks factor':>9} {'KS_MULT':>8} | {'fw vol':>6} {'d':>7} | {'f27 vol':>7} {'d':>7} | {'f27 pk%':>7} {'d':>6} | {'fw pk%':>6} | {'days':>4} {'tail1':>6} | {'storms_ok':>9} {'peaks_ok':>8}")
        r0 = got[0][2] if got[0][0] == 1.0 else None
        for fac, j, m in got:
            row = scores_row(man, lab, j["run_id"], csv_scores)
            d = (lambda k: (m[k] - r0[k]) if r0 is not None else NAN)
            T.log("      " + f"{fac:>9.2f} {j['changes'].get('KS_MULT', T.BASE['KS_MULT']):>8.5g} | {T.fm(m['fw_vol'], 3):>6} {T.fm(d('fw_vol'), 3, True):>7} | "
                             f"{T.fm(m['f27_vol'], 3):>7} {T.fm(d('f27_vol'), 3, True):>7} | {T.fm(m['f27_pk'], 1):>7} {T.fm(d('f27_pk'), 1, True):>6} | "
                             f"{T.fm(m['fw_pk'], 1):>6} | {T.fm(m['days'], 0):>4} {T.fm(m['tail1'], 2):>6} | {T.yn(T.storms_ok(row)):>9} {T.yn(peaks_ok(row)):>8}")
        pts = [(fac, m) for fac, _, m in got]
        e_fw, e_f27 = ks_fit([(f, m["fw_vol"]) for f, m in pts], "fw_vol"), ks_fit([(f, m["f27_vol"]) for f, m in pts], "f27_vol")
        slope = ks_peak_slope([(f, m["f27_pk"]) for f, m in pts])
        T.log("      found: " + ", ".join(x for x in (
            f"flood-week volume ~ Ks^{e_fw:+.2f} (prior {KS_PRIOR['fw_vol']:+.1f})" if e_fw is not None else "",
            f"Feb 27-28 volume ~ Ks^{e_f27:+.2f} (prior {KS_PRIOR['f27_vol']:+.1f})" if e_f27 is not None else "",
            f"Feb 27-28 peak error {slope:+.2f} points per +1% Ks" if slope is not None else "") if x))
    if printed:
        T.log("      (the exponents come from the runs listed, 2 to 4 points; the prior came from 169 runs at As 1. A large difference means the As 1 prior does not carry over. Three points are a rough fit.)")


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
            new = ensure_rise(man, lab, rid, gauge_q) or new
    if new:
        T.save_manifest(man)
    done = [(rid, rec, scores_row(man, lab, rid, csv_scores)) for rid, rec in recs if rec.get("status") == "done"]
    T.log()
    T.log(f"W4 RE-TUNE (Oct 1 2016 - Apr 1 2017, {T.WIN['hours']} h from a {T.WIN['gw']:g} x bedrock water table; channel loss off). Baseline = draw 12 / k10 (Ks x{T.BASE['KS_MULT']:g}, f x{T.BASE['F_MULT']:g}, As 1).")
    T.log("  vol = model volume / gauge volume. pk = model peak vs gauge peak, %. arr = scorer arrival error, hours (negative = early; see the TIMING table for what it hides).")
    T.log(f"  storms_ok = your anchor rule (flood-week vol {T.ANCHOR_FW[0]:g}-{T.ANCHOR_FW[1]:g}, Feb 27-28 vol {T.ANCHOR_F27[0]:g}-{T.ANCHOR_F27[1]:g}). peaks_ok = flood-week and Feb 27-28 peak errors both within 15%.")
    if not done:
        T.log()
        T.log("  (no scored runs yet)")
        return

    T.log()
    T.log("  STORMS")
    T.log("    " + f"{'run':<26}| {'fw vol':>6} {'fw pk%':>7} {'fw arr':>6} | {'f27 vol':>7} {'f27 pk%':>8} {'f27 arr':>7} | {'f18 vol':>7} {'2nd Jan vol':>11} | {'storms_ok':>9} {'peaks_ok':>8} | note")
    for rid, rec, row in done:
        note = "IDENTICAL to baseline" if rec.get("identical_to_baseline") else ""
        T.log("    " + f"{label_of(rec):<26}| {T.fm(T.fnum(row, 'flood_week__vol_ratio'), 3):>6} {T.fm(T.fnum(row, 'flood_week__peak_err_pct'), 1):>7} "
                       f"{T.fm(T.fnum(row, 'flood_week__arrival_err_h')):>6} | {T.fm(T.fnum(row, 'feb27_28__vol_ratio'), 3):>7} "
                       f"{T.fm(T.fnum(row, 'feb27_28__peak_err_pct'), 1):>8} {T.fm(T.fnum(row, 'feb27_28__arrival_err_h')):>7} | "
                       f"{T.fm(T.fnum(row, 'feb18_19__vol_ratio'), 2):>7} {T.fm(T.fnum(row, 'jan_second_storm__vol_ratio'), 2):>11} | "
                       f"{T.yn(T.storms_ok(row)):>9} {T.yn(peaks_ok(row)):>8} | {note}")

    first_tail = next((rec["tail"] for _, rec, _ in done if rec.get("tail")), None)
    T.log()
    if first_tail:
        T.log(f"  TAIL  (mean flow in 24 h blocks, m3/s. Blocks start {T.TAIL_OFFSETS_D[0]}, {T.TAIL_OFFSETS_D[1]}, {T.TAIL_OFFSETS_D[2]} days after the gauge flood peak "
              f"({first_tail.get('gauge_peak_time', '')}); band {args.band_lo:g}-{args.band_hi:g} m3/s; days in band counts the {first_tail.get('days_total', '?')} blocks of 24 h from {T.TAIL_START[:10]} to {T.TAIL_END[:10]})")
        T.log("    " + f"{'run':<26}| " + " ".join(f"{'d' + str(d):>7}" for d in T.TAIL_OFFSETS_D) + f" | {'days in band':>12} | {'tail1 vol':>9} {'tail2 vol':>9} | {'tail_ok':>7} | {'Hsrf share':>10} {'decay x base':>12}")
        T.log("    " + f"{'gauge':<26}| " + " ".join(f"{T.fm(first_tail.get(f'gauge_d{d}_m3s'), 3):>7}" for d in T.TAIL_OFFSETS_D)
              + f" | {str(first_tail.get('gauge_days_in_band', '?')) + '/' + str(first_tail.get('days_total', '?')):>12} |")
        base_decay = ((man["runs"].get(T.BASELINE_ID) or {}).get("decay") or {}).get("decay_median")
        for rid, rec, row in done:
            tail = rec.get("tail")
            if not tail:
                continue
            dm = (rec.get("decay") or {}).get("decay_median")
            rel = (dm / base_decay) if (dm is not None and base_decay) else NAN
            T.log("    " + f"{label_of(rec):<26}| " + " ".join(f"{T.fm(tail.get(f'model_d{d}_m3s'), 3):>7}" for d in T.TAIL_OFFSETS_D)
                  + f" | {str(tail.get('model_days_in_band', '?')) + '/' + str(tail.get('days_total', '?')):>12} | {T.fm(T.fnum(row, 'tail_jan26_feb2__tail_vol_ratio')):>9} "
                  + f"{T.fm(T.fnum(row, 'tail_mar3_mar11__tail_vol_ratio')):>9} | {T.yn(tail.get('tail_ok')):>7} | {T.fm(T.fnum(row, 'proc__hsrf_share'), 3):>10} {T.sci(rel):>12}")
    else:
        T.log("  TAIL  (no tail numbers yet)")

    rises = [(rid, rec) for rid, rec, _ in done if rec.get("rise")]
    if rises:
        g = next(r["rise"] for _, r in rises)
        T.log()
        T.log("  TIMING  (hours, model minus gauge; negative = early. t20/t50/t80 = when the flow first reaches 20/50/80% of the gauge peak; pk = time of the peak)")
        T.log("    " + f"{'run':<26}| {'FLOOD WEEK  pk':>14} {'t20':>6} {'t50':>6} {'t80':>6} | {'FEB 27-28  pk':>13} {'t20':>6} {'t50':>6} {'t80':>6} | {'Feb 28 shoulder: max m3/s':>25} {'first 5 m3/s':>13}")
        for rid, rec in rises:
            r = rec["rise"]
            T.log("    " + f"{label_of(rec):<26}| {T.fm(r.get('fw_pk_dt_h'), 2, True):>14} {T.fm(r.get('fw_t20_h'), 2, True):>6} {T.fm(r.get('fw_t50_h'), 2, True):>6} {T.fm(r.get('fw_t80_h'), 2, True):>6} | "
                           f"{T.fm(r.get('f27_pk_dt_h'), 2, True):>13} {T.fm(r.get('f27_t20_h'), 2, True):>6} {T.fm(r.get('f27_t50_h'), 2, True):>6} {T.fm(r.get('f27_t80_h'), 2, True):>6} | "
                           f"{T.fm(r.get('shoulder_max_m3s'), 1):>25} {r.get('shoulder_first5', ''):>13}")
        T.log("    " + f"{'gauge':<26}| {'':>14} {'':>6} {'':>6} {'':>6} | {'':>13} {'':>6} {'':>6} {'':>6} | {T.fm(g.get('gauge_shoulder_max_m3s'), 1):>25} {g.get('gauge_first5', ''):>13}")
        T.log("    (a shoulder far above the gauge's, reached hours before the gauge's first 5 m3/s, is what pulls the scorer's f27 arrival early; compare t50/t80)")

    print_ks_response(man, lab, csv_scores)

    T.log()
    T.log("  PROCESS  (water table = depth below the surface in mm; 'decay' is exp(-f x depth) at the water table, median over the basin, f per mm assumed)")
    T.log("    " + f"{'run':<26}| {'water table mm: first':>21} {'shallowest':>10} {'deepest':>8} | {'f median':>8} {'decay':>9} | {'build min':>9} {'tRIBS min':>9}")
    for rid, rec, row in done:
        d = rec.get("decay") or {}
        T.log("    " + f"{label_of(rec):<26}| {T.fm(T.fnum(row, 'proc__mdgw_first_mm'), 0):>21} {T.fm(T.fnum(row, 'proc__mdgw_min_mm'), 0):>10} {T.fm(T.fnum(row, 'proc__mdgw_max_mm'), 0):>8} | "
                       f"{T.fm(d.get('f_median'), 4):>8} {T.sci(d.get('decay_median', NAN) if d.get('decay_median') is not None else NAN):>9} | {str(rec.get('build_min', '')):>9} {str(rec.get('tribs_min', '')):>9}")
    failed = [rid for rid, rec in man["runs"].items() if rec.get("status") == "failed"]
    if failed:
        T.log()
        T.log("  FAILED runs (see the log): " + ", ".join(failed))
    write_summary_csv(man, lab, csv_scores)
    T.log()
    T.log(f"  Wrote {T.SUMMARY_CSV}. These tables decide nothing: whether a Ks step is worth its cost to the tail and the Feb 27-28 peak is your call.")


def write_summary_csv(man: dict, lab: Path, csv_scores: dict) -> None:
    rows = []
    for rid, rec in sorted_runs(man):
        row_sc = scores_row(man, lab, rid, csv_scores)
        s = rec.get("settings") or {}
        r = dict(run_id=rid, tag=rid[len(PREFIX) + 1:-len("_nocl")] if rid.endswith("_nocl") else rid, step=rec.get("step"), status=rec.get("status"),
                 F_MULT=s.get("F_MULT"), ANISOTROPY=s.get("ANISOTROPY"), ANISOTROPY_U=s.get("ANISOTROPY_U"), KS_MULT=s.get("KS_MULT"),
                 ks_factor=(s.get("KS_MULT") / T.BASE["KS_MULT"]) if s.get("KS_MULT") else "", 
                 storms_ok=(T.storms_ok(row_sc) if rec.get("status") == "done" else ""), peaks_ok=(peaks_ok(row_sc) if rec.get("status") == "done" else ""),
                 tail_ok=((rec.get("tail") or {}).get("tail_ok", "")))
        for c in T.SCORE_COLS:
            r[c] = T.fnum(row_sc, c) if row_sc else NAN
        r.update({k: v for k, v in (rec.get("tail") or {}).items() if k not in ("band_lo", "band_hi")})
        r.update({f"rise_{k}": v for k, v in (rec.get("rise") or {}).items()})
        d = rec.get("decay") or {}
        r.update(f_median=d.get("f_median"), gw_depth_median_mm=d.get("gw_depth_median_mm"), decay_median=d.get("decay_median"),
                 identical_to_baseline=bool(rec.get("identical_to_baseline")), build_min=rec.get("build_min", ""), tribs_min=rec.get("tribs_min", ""),
                 output_mb=rec.get("output_mb", ""))
        rows.append(r)
    cols = []
    for r in rows:
        for c in r:
            if c not in cols:
                cols.append(c)
    with open(lab / T.SUMMARY_CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: ("" if (isinstance(r.get(c), float) and not math.isfinite(r[c])) else r.get(c, "")) for c in cols})


# --------------------------------------------------------------------------------------------
# The figure
# --------------------------------------------------------------------------------------------
DEFAULT_PLOT = ["base", "f0p15_as300", "f0p15_as300_ksx1p03", "f0p15_as300_ksx1p06"]
PLOT_STYLE = [("#2a78d6", "-"), ("#eb6834", "-"), ("#7b4fc4", "--"), ("#1baf7a", ":"), ("#b8860b", "-."), ("#c43e7e", "--")]


def make_plot(man: dict, lab: Path, tags: list, out: Path) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    runs = []
    for t in tags:
        job = JOBS.get(t)
        if job is None:
            raise ValueError(f"--plot-runs: '{t}' is not a run of this script (choices: {', '.join(JOBS)})")
        rec = man["runs"].get(job["run_id"]) or {}
        qf = lab / "results" / f"{job['run_id']}_Outlet.qout"
        if rec.get("status") == "done" and qf.is_file():
            runs.append((label_of(rec), T.read_qout(qf, T.RUN_START_TS)))
    if not runs:
        return "no finished run with an outlet file among the plotted ones; no figure made"
    gauge = T.read_obs(lab / T.INIT_DATA / T.OBS_FILE)
    panels = [("A. Flood week", "2017-01-18", "2017-01-25", False), ("B. Feb 27-28 storm", "2017-02-27 12:00", "2017-03-01 12:00", False),
              ("C. The tail (log scale; shaded = the gauge's 0.3-1 m3/s band)", "2017-01-26", "2017-02-18", True)]
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8), gridspec_kw={"wspace": 0.22})
    for ax, (title, a, b, logy) in zip(axes, panels):
        g = _cut(gauge, a, b)
        ax.plot(g.index, g.values, color="black", lw=2.2, label="gauge")
        for k, (lbl, q) in enumerate(runs):
            color, ls = PLOT_STYLE[k % len(PLOT_STYLE)]
            s = _cut(q, a, b)
            ax.plot(s.index, s.values, color=color, ls=ls, lw=1.7, label=lbl)
        if logy:
            ax.set_yscale("log")
            ax.set_ylim(0.01, 30)
            ax.axhspan(0.3, 1.0, color="#8a8984", alpha=0.15, lw=0)
        ax.set_title(title, loc="left", fontsize=11)
        ax.set_ylabel("Outlet flow (m$^3$/s)")
        ax.grid(alpha=0.25)
        ax.tick_params(axis="x", labelrotation=30)
        for s_ in ("top", "right"):
            ax.spines[s_].set_visible(False)
    axes[0].legend(frameon=False, fontsize=8, loc="upper right")
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return f"wrote {out.name} with {len(runs)} run(s) and the gauge"


# --------------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------------
def parse_args(argv=None):
    ap = argparse.ArgumentParser(description="Small Ks re-tune at f x0.15 / As 300 on the full window W4 (see the top of this file).")
    ap.add_argument("--dry-run", action="store_true", help="check everything, list the plan, build nothing")
    ap.add_argument("--summary-only", action="store_true", help="print the tables from the runs already scored, then stop")
    ap.add_argument("--build-only", action="store_true", help="build and check the inputs but do not run tRIBS (use with --limit 1)")
    ap.add_argument("--limit", type=int, default=None, help="make only the first N runs that still need work (finished runs do not count)")
    ap.add_argument("--steps", default=DEFAULT_STEPS, help=f"comma-separated steps (only ks exists; default {DEFAULT_STEPS}; the baseline is checked first)")
    ap.add_argument("--plot", action="store_true", help="make retune_w4_figs/fig_w4_overlays.png from the finished runs")
    ap.add_argument("--plot-runs", default=",".join(DEFAULT_PLOT), help="runs to draw, as the tag part of the run name (default %(default)s)")
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
    return ap.parse_args(argv)


def main(argv=None) -> int:
    lost = missing_from_w4_script()
    if lost:
        print("run_tail_screen_sycamore.py in this folder is an older copy: it lacks " + ", ".join(lost) + ". Put the newest run_tail_screen_sycamore.py in lab/ and run again.")
        return 2
    args = parse_args(argv)
    lab = Path(__file__).resolve().parent
    os.chdir(lab)
    (lab / T.LOG_DIR).mkdir(exist_ok=True)
    T._LOGFILE = str(lab / T.LOG_DIR / "batch_log.txt")
    steps = [x.strip() for x in args.steps.split(",") if x.strip()]
    bad = [x for x in steps if x not in STEP_ORDER]
    if bad or not steps:
        print(f"--steps must be a comma-separated list of {', '.join(STEP_ORDER)} (got {bad or 'nothing'}).")
        return 2
    plot_tags = [x.strip() for x in args.plot_runs.split(",") if x.strip()]
    badtags = [t for t in plot_tags if t not in JOBS]
    if badtags:
        print(f"--plot-runs: {', '.join(badtags)} is not a run of this script (choices: {', '.join(JOBS)}).")
        return 2
    man = T.load_manifest()

    def finish_outputs() -> None:
        if args.plot:
            try:
                T.log("  Figure: " + make_plot(man, lab, plot_tags, lab / FIG_DIR / "fig_w4_overlays.png"))
            except Exception as exc:
                T.log(f"  Figure could not be made: {type(exc).__name__}: {exc}")

    if args.summary_only:
        print_report(man, lab, args)
        finish_outputs()
        return 0

    args.notebook, nb_note = T.pick_notebook(lab, args.notebook)
    args.redo_ids = set(x.strip() for x in (args.redo or "").split(",") if x.strip())
    preview = planned(steps)
    n_done = sum(1 for j in preview if T.is_finished(man, j, args))
    n_todo = len(preview) - n_done
    per_lo, per_hi = T.WIN["hours"] * T.MIN_PER_1000H / 1000.0, T.WIN["hours"] * T.MIN_PER_1000H_SLOW / 1000.0
    T.log(f"Working folder: {lab}")
    T.log(f"Notebook: {args.notebook}")
    T.log(f"Plan: steps {','.join(steps)}, window W4 ({T.WIN['hours']} h from {T.WIN['start']}, water table {T.WIN['gw']:g} x bedrock): {n_todo} run(s) to make"
          + (f", {n_done} already finished" if n_done else "") + ".")
    for j in preview:
        what = "already done, skipped" if T.is_finished(man, j, args) else "NEW RUN"
        T.log(f"  {j['run_id']:48s} {describe(j):52s} {what}")
    T.log("  (the baseline and f x0.15 / As 300 are the tail screen's own runs; they are skipped when the shared manifest says they are finished)")
    T.log(f"  tRIBS time alone: about {n_todo * per_lo / 60:.1f}-{n_todo * per_hi / 60:.1f} hours ({per_lo:.0f}-{per_hi:.0f} min per run); each run also needs a notebook build of about 1 min.")

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
    man["retune_w4_last_started"] = T.now_iso()
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
        code = run_batch([JOBS["base"]])                       # the baseline always comes first: every other run's .in file is checked against it
        for st in STEP_ORDER:
            if st not in steps:
                continue
            if code != 0:
                T.log(f"The {st} step was not started because the baseline failed or was interrupted. Fix that, then run the same command again.")
                continue
            todo = list(STEP_JOBS[st])
            new_runs = [describe(j) for j in todo if not T.is_finished(man, j, args)]
            if new_runs:
                T.log()
                T.log(f"{st} step: " + "; ".join(new_runs))
            code = run_batch(todo) or code
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
    T.log(f"{'run':52s} {'status':16s} {'build':>7s} {'tRIBS':>7s}")
    for rid, rec in man["runs"].items():
        T.log(f"{rid:52s} {rec.get('status', '?'):16s} {str(rec.get('build_min', '')):>7s} {str(rec.get('tribs_min', '')):>7s}")
    print_report(man, lab, args)
    finish_outputs()
    return code


if __name__ == "__main__":
    sys.exit(main())
