#!/usr/bin/env python3
"""
run_retune_w2_sycamore.py  --  The small Ks re-tune at f x0.15 / As 300, on the SHORT window W2 (Jan 18 - Mar 12 2017).

WHERE THIS COMES FROM. The tail screen (run_tail_screen_sycamore.py, on W4) showed that Josh's two levers work together: a lower f
plus a higher SATURATED anisotropy ratio As lifts the slow flow after the January flood (6 of 23 days in the gauge's 0.3-1 m3/s band
with the draw 12 / k10 setting, 16 of 23 with f x0.15 and As 300). The price is a higher Feb 27-28 peak and volume. You decided to
keep f and As and re-tune Ks (and, if needed, PsiB and routing) by a much smaller amount than before. A W4 run takes about 34 minutes;
a W2 run takes about 10 to 13. This script does that re-tune on W2.

WHAT W2 CANNOT SEE. W2 starts cold on Jan 18 (water table 0.82 x bedrock) and has no autumn, so it cannot show the Oct - Dec false
flow. That stays a W4 question. W2 is also not guaranteed to rank settings the way W4 does, because its water-table start state is
different. So the script FIRST checks, with the five settings that already have W4 results, whether W2 can stand in for W4. That
check has a rule that is written here, BEFORE any run, so it cannot be bent after the numbers are in (the W1B lesson).

IT REUSES run_tail_screen_sycamore.py. Both files must be in lab/. Every build, check, tRIBS run and scoring step is the tested
W4-script code; this file only changes the window (to W2), the run names and the file names, defines the run lists, and prints
its own tables. Nothing here changes your notebook or any existing result.

WINDOW AND START STATE (set in this file)
    W2: START 01/18/2017/00/00, 1,272 hours (to Mar 12), water table 0.82 x bedrock, bedrock +1 m, channel loss OFF, snow ON.
    Starting values are the draw 12 / k10 setting: Ks x0.187, PsiB x0.2939, f x0.5874, m x1, thetaS x1, As = Au = 1,
    cv 6.39, flowexp 0.2973, n 0.067, width coefficient 2.33. The baseline run is the first run of every batch.
    Run names start with W2_rt_ (so the scorer reads the window W2 from the name; "rt" = re-tune, to keep them apart from
    any W2 runs you already have).

STEPS (choose with --steps; default proxy,ks)
    proxy  5 runs, about 1 hour. The baseline, then f x0.15 As 300, f x0.02 As 300, f x0.5874 As 300, f x0.02 As 500 (Au = 1).
           These five settings already have W4 results. The script sets W2 beside W4 and applies THE RULE below.
    ks     5 runs, about 1 hour. f x0.15 / As 300 at Ks x1.03, x1.06, x1.10 of the baseline Ks (KS_MULT 0.19261, 0.19822, 0.2057), and
           f x0.02 / As 300 at Ks x1.06 and x1.10. The f x0.15 and f x0.02 As 300 runs of the proxy step are the "Ks x1.00" points.
           Why these sizes: the regression priors (Feb 27-28 volume goes as Ks^-1.9, flood week as Ks^-0.9, found at As 1 and NOT
           tested at As 300) say Ks x1.03 to x1.06 should bring the Feb 27-28 volume from 1.07 to about 1.0 to 0.95. The run
           results replace the prior: the script prints the Ks exponents it actually finds, next to the prior.
    grid   (opt-in) 3 runs. f x0.15 with As 200 and As 500, and f x0.25 with As 300. The neighborhood of the working setting:
           how fragile is it? (At f x0.02 the response to As is a narrow hump: 100 -> 4 days in band, 300 -> 21, 900 -> 2.)
    phys   (opt-in) 4 runs, all at f x0.15 / As 300. Start state: water table 0.70 and 0.95 x bedrock. Bedrock: +0 m and +2 m
           (the water table is a fraction of bedrock depth, so the bedrock change also moves the water table; Josh's untested
           third point, and relevant to the WY2020 start state). Is the slow flow a property of the soil setting or of the start?

THE RULE (written before the runs; the same text is printed with the verdict)
    The five proxy settings have W4 results (the Oct 1 - Apr 1 tail screen). Six numbers are compared, W2 against W4:
        flood-week volume ratio, flood-week peak error, Feb 27-28 volume ratio, Feb 27-28 peak error (points),
        days in the 0.3-1 m3/s band (Jan 26 - Feb 18), first-tail volume ratio (Jan 26 - Feb 2).
    (Span KGE is NOT compared: the two windows cover different spans.)
    R1 ORDER. Take every pair of the five settings that W4 separates by more than a resolution (0.025 for volumes, 2 points for
       peak errors, 3 days, 0.03 for the tail). Each such pair must keep its order in W2. One flipped pair on any metric = FAIL.
    R2 OFFSET SPREAD. offset = W2 minus W4, per setting. The spread (largest minus smallest offset over the five settings) of
       the four storm numbers must stay within: volumes 0.04, peak errors 4 points. (For days in band and tail volume the spread
       is printed but does not decide, because the Ks step reads them as changes inside W2.)
    R3 CONSTANT OFFSET. If R1 and R2 hold but the mean offset of a storm number is larger than 0.02 (volumes) or 2 points (peaks),
       the verdict is still PASS, and that metric is read as "W4-equivalent = W2 value minus the mean offset" in a separately
       labeled table (marked W4eq). No other correction is ever applied.
    VERDICTS.  PASS: R1 and R2 hold; the ks and grid steps may run. ORDER ONLY: R1 holds, R2 fails; W2 can rank settings but its
               values cannot be read; the ks and grid steps are NOT started, and the Ks re-tune is done on W4. FAIL: R1 fails;
               same consequence. INCOMPLETE: fewer than five proxy runs scored; nothing is decided.
    Override with --ignore-verdict (it is printed in the log). Whatever W2 says, the best one or two settings are confirmed on W4
    (about 34 minutes each) before they are adopted.

WHAT IT PRINTS (all read from files; nothing is judged for you except the rule above)
    proxy check   the five settings, W4 beside W2, the order test, the offset spreads, the verdict.
    storms        flood week, 2nd January storm, Feb 18-19, Feb 27-28: volume ratio, peak error %, scorer arrival error.
                  storms_ok = your anchor rule (flood-week volume 0.90-1.12 and Feb 27-28 volume 0.90-1.10); peaks_ok = both peak
                  errors within +-15%.
    tail          mean flow in the 24 h blocks 7, 14, 21 days after the gauge flood peak; days in the band; tail volume ratios.
    timing        the rising limb of the flood week and of Feb 27-28: when the flow reaches 20, 50, 80% of the gauge peak and when
                  the peak itself comes, vs the gauge (hours; negative = early). Also the early SHOULDER on Feb 28: the highest flow
                  between 00:00 and 09:30 and the time the flow first reaches 5 m3/s (gauge: 6.3 m3/s at most, first 5 at 08:05).
                  The scorer's "arrival" (first 20% crossing) is thrown off by that shoulder; the 50% and 80% times are not.
    Ks response   for each (f, As) group: the runs side by side, change from the Ks x1.00 run, and the exponents found.
    process       infiltration-excess share, water table depth, the depth-decay number.
    W4eq table    only when R3 flags a constant offset.
    --plot makes retune_w2_figs/fig_w2_overlays.png (flood week, Feb 27-28 and the tail, runs against the gauge).

WHAT IT DOES NOT DO. It does not re-tune PsiB or routing (the Ks step comes first; if Ks cannot fix the Feb 27-28 peak, PsiB and
routing are the next, smaller steps). It does not test the Oct - Dec false flow, the Jan 23-24 storm, rain scaling, or snow.

Run it from workspaces/Sycamore/lab/, in the foreground (keep the terminal open: the progress line printed every 5 minutes also
keeps the Codespace awake):
    python -u run_retune_w2_sycamore.py --dry-run                      # checks only, builds nothing
    python -u run_retune_w2_sycamore.py --build-only --limit 1         # build the W2 baseline and check it, no tRIBS
    python -u run_retune_w2_sycamore.py --limit 2                      # baseline + f x0.15 As 300, then look at the lines
    python -u run_retune_w2_sycamore.py                                # proxy then ks (the ks step starts only if the verdict is PASS)
    python -u run_retune_w2_sycamore.py --steps grid,phys              # the opt-in steps
    python -u run_retune_w2_sycamore.py --summary-only --plot          # tables and figure from what is finished
Run the same command again to resume after any stop: finished runs are skipped.

Files it makes (in lab/): retune_w2_manifest.json, retune_w2.lock (only while running), retune_w2_logs/ (batch_log.txt and one build
output and score file per run), score_out_retune_w2/scores_retune_w2.csv, retune_w2_summary.csv, retune_w2_figs/, and for each run
<run>.in, simulation_<run>.log and results/<run>*. Do NOT run the notebook (or anything that rebuilds data/model/) while it runs.
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
# The W2 window and the file names (these replace the W4 values inside the imported module)
# --------------------------------------------------------------------------------------------
W2 = dict(name="W2", start="01/18/2017/00/00", hours=1272, label="jan_to_mar", gw=0.82)
W2_START_TS = "2017-01-18 00:00"
PREFIX = "W2_rt_d012_k10"                       # must start with W2_ so the scorer reads the window from the name
W4_PREFIX = "W4_d012_k10"
FILES = dict(LOG_DIR="retune_w2_logs", SCORE_DIR="score_out_retune_w2", SCORE_CSV="scores_retune_w2.csv",
             SUMMARY_CSV="retune_w2_summary.csv", MANIFEST="retune_w2_manifest.json", LOCK="retune_w2.lock")
FIG_DIR = "retune_w2_figs"
W4_SUMMARY = "tail_screen_w4_summary.csv"
SCRIPT_NAME = "run_retune_w2_sycamore.py"
STEP_ORDER = ["proxy", "ks", "grid", "phys"]
DEFAULT_STEPS = "proxy,ks"
GATED_STEPS = ("ks", "grid")                    # these need a PASS from the proxy check (phys does not)


def describe(job: dict) -> str:
    """The text the W4 script prints for a run (it looks this function up by name, so it is replaced here)."""
    if job["step"] == "base":
        return "baseline (no change)"
    ch = job["changes"]
    parts = [f"f x{job['f']:g}", f"As {job['an']:g}"]
    if "KS_MULT" in ch:
        parts.append(f"Ks x{ch['KS_MULT'] / T.BASE['KS_MULT']:.2f} (KS_MULT {ch['KS_MULT']:g})")
    if "GW_FRACTION" in ch:
        parts.append(f"start water table {ch['GW_FRACTION']:g} x bedrock")
    if "BEDROCK_ADD_M" in ch:
        parts.append(f"bedrock +{ch['BEDROCK_ADD_M']:g} m")
    return ", ".join(parts)


def configure() -> None:
    T.WIN = dict(W2)
    T.RUN_START_TS = W2_START_TS
    T.BASE["GW_FRACTION"] = W2["gw"]
    T.RUN_PREFIX = PREFIX
    T.BASELINE_ID = f"{PREFIX}_base_nocl"
    for name, value in FILES.items():
        setattr(T, name, value)
    T.describe = describe


configure()

# The names this file takes from run_tail_screen_sycamore.py. If one is missing, lab/ holds an older copy of that script.
NEEDED_FROM_W4_SCRIPT = ["do_run", "preflight", "make_base_job", "an_job", "au_job", "tail_metrics", "read_qout", "read_obs", "compute_tail", "load_scores",
                         "expected_files", "settings_record", "settings_for", "is_finished", "take_lock", "save_manifest", "load_manifest", "pick_notebook",
                         "install_term_handler", "fnum", "fm", "yn", "sci", "storms_ok", "log", "now_iso", "fmt_min", "vtag", "au_of", "SCORE_COLS",
                         "ANCHOR_FW", "ANCHOR_F27", "TAIL_START", "TAIL_END", "TAIL_OFFSETS_D", "BAND_LO", "BAND_HI", "HEARTBEAT_MIN", "TRIBS_DEFAULT",
                         "SCORER", "INIT_DATA", "OBS_FILE", "MIN_PER_1000H", "MIN_PER_1000H_SLOW", "RunError", "Abort"]


def missing_from_w4_script() -> list:
    return [n for n in NEEDED_FROM_W4_SCRIPT if not hasattr(T, n)]


# --------------------------------------------------------------------------------------------
# The runs
# --------------------------------------------------------------------------------------------
PROXY_RUN_ORDER = ["base", "f0p15_as300", "f0p02_as300", "f0p5874_as300", "f0p02_as500"]   # order they are RUN in
PROXY_SHOW_ORDER = ["base", "f0p5874_as300", "f0p15_as300", "f0p02_as300", "f0p02_as500"]  # order they are SHOWN in (W4 Feb 27-28 peak, low to high)
PROXY_SETTINGS = {"f0p5874_as300": (0.5874, 300), "f0p15_as300": (0.15, 300), "f0p02_as300": (0.02, 300), "f0p02_as500": (0.02, 500)}
KS_PLAN = [(0.15, 300, 1.06), (0.15, 300, 1.03), (0.15, 300, 1.10), (0.02, 300, 1.06), (0.02, 300, 1.10)]   # (f, As, Ks factor on the baseline Ks)
GRID_PLAN = [(0.15, 200), (0.15, 500), (0.25, 300)]
PHYS_PLAN = [("gw0p7", {"GW_FRACTION": 0.70}), ("gw0p95", {"GW_FRACTION": 0.95}),
             ("bed0", {"BEDROCK_ADD_M": 0.0}), ("bed2", {"BEDROCK_ADD_M": 2.0})]                      # all at f x0.15, As 300
KS_PRIOR = {"fw_vol": -0.9, "f27_vol": -1.9}      # exponents of volume against Ks from the regression at As 1 (150 W1 and 19 W2 runs)


def job_tag(job: dict) -> str:
    """'W2_rt_d012_k10_f0p15_as300_ksx1p06_nocl' -> 'f0p15_as300_ksx1p06'"""
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
    for tag, (fv, a) in PROXY_SETTINGS.items():
        j = T.an_job(fv, a)
        j["step"] = "proxy"
        J[tag] = j
    ks_jobs = []
    for fv, a, factor in KS_PLAN:
        j = _with_tag(T.an_job(fv, a), f"ksx{T.vtag(factor)}", {"KS_MULT": round(T.BASE["KS_MULT"] * factor, 6)}, "ks")
        j["ks_factor"] = factor
        ks_jobs.append(j)
    grid_jobs = []
    for fv, a in GRID_PLAN:
        j = T.an_job(fv, a)
        j["step"] = "grid"
        grid_jobs.append(j)
    phys_jobs = [_with_tag(T.an_job(0.15, 300), tag, ch, "phys") for tag, ch in PHYS_PLAN]
    refs = [J["f0p15_as300"], J["f0p02_as300"]]
    steps = {"proxy": [J[t] for t in PROXY_RUN_ORDER],
             "ks": refs + ks_jobs,
             "grid": [J["f0p15_as300"]] + grid_jobs,
             "phys": [J["f0p15_as300"]] + phys_jobs}
    show = [J[t] for t in PROXY_SHOW_ORDER] + sorted(ks_jobs, key=lambda j: (-j["f"], j["ks_factor"])) + grid_jobs + phys_jobs
    for j in ks_jobs + grid_jobs + phys_jobs:
        J[job_tag(j)] = j
    return J, steps, show


JOBS, STEP_JOBS, SHOW_JOBS = build_jobs()
SHOW_ORDER = {j["run_id"]: i for i, j in enumerate(SHOW_JOBS)}


def planned(steps: list) -> list:
    """The baseline first, then the runs of each chosen step, no run twice (the reference runs are shared between steps)."""
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
    """The six numbers of the proxy rule for a finished run, or None if the run is not done or one number is missing."""
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
# The proxy check: does W2 rank the five settings like W4?
# --------------------------------------------------------------------------------------------
# W4 values of the five proxy settings, copied from Tail_Screen_W4_Result_v1.md (scorer numbers; days in band and tail from the script's own tables).
# If tail_screen_w4_summary.csv is in lab/ and holds a setting, the CSV values are used instead for that setting.
EMBEDDED_W4 = {
    "base":          dict(fw_vol=1.101, fw_pk=-5.1, f27_vol=0.934, f27_pk=5.4,  days=6,  tail1=0.27),
    "f0p5874_as300": dict(fw_vol=1.110, fw_pk=-4.9, f27_vol=1.001, f27_pk=11.6, days=5,  tail1=0.254),
    "f0p15_as300":   dict(fw_vol=1.078, fw_pk=-9.4, f27_vol=1.065, f27_pk=15.5, days=16, tail1=0.35),
    "f0p02_as300":   dict(fw_vol=1.081, fw_pk=-9.7, f27_vol=1.113, f27_pk=20.3, days=21, tail1=0.43),
    "f0p02_as500":   dict(fw_vol=1.161, fw_pk=-4.8, f27_vol=1.198, f27_pk=30.2, days=19, tail1=0.588),
}
# key, label, W4 pair resolution (R1), offset-spread tolerance (R2), constant-offset tolerance (R3), does R2 decide?
KEY_METRICS = [
    ("fw_vol",  "flood-week volume",             0.025, 0.04, 0.02, True),
    ("fw_pk",   "flood-week peak error, points", 2.0,   4.0,  2.0,  True),
    ("f27_vol", "Feb 27-28 volume",              0.025, 0.04, 0.02, True),
    ("f27_pk",  "Feb 27-28 peak error, points",  2.0,   4.0,  2.0,  True),
    ("days",    "days in the slow-flow band",    3,     3,    2,    False),
    ("tail1",   "first-tail volume",             0.03,  0.05, 0.03, False),
]
CSV_COLS = dict(fw_vol="flood_week__vol_ratio", fw_pk="flood_week__peak_err_pct", f27_vol="feb27_28__vol_ratio",
                f27_pk="feb27_28__peak_err_pct", days="model_days_in_band", tail1="tail_jan26_feb2__tail_vol_ratio")


def w4_reference(lab: Path, csv_path: str | None = None) -> tuple:
    """(reference dict by proxy tag, list of notes). A setting found in the W4 summary CSV takes its values from there."""
    ref = {t: dict(v, source="Result doc") for t, v in EMBEDDED_W4.items()}
    notes = []
    p = Path(csv_path) if csv_path else lab / W4_SUMMARY
    if p.is_file():
        try:
            df = pd.read_csv(p)
            used = []
            for t in ref:
                row = df[df["run_id"] == f"{W4_PREFIX}_{t}_nocl"] if "run_id" in df.columns else df.iloc[0:0]
                if len(row) != 1:
                    continue
                r = row.iloc[0]
                vals = {k: pd.to_numeric(r.get(c), errors="coerce") for k, c in CSV_COLS.items()}
                if all(pd.notna(v) for v in vals.values()):
                    ref[t].update({k: float(v) for k, v in vals.items()}, source=p.name)
                    used.append(t)
            notes.append(f"W4 values from {p.name} for {len(used)} of 5 settings ({', '.join(used) or 'none'}); the rest from the W4 result document")
        except Exception as exc:
            notes.append(f"{p.name} could not be read ({type(exc).__name__}: {exc}); all W4 values are from the W4 result document")
    else:
        notes.append(f"{p.name} not found in lab/; all W4 values are from the W4 result document")
    return ref, notes


def order_test(w4: list, w2: list, res: float) -> tuple:
    """(pairs tested, pairs flipped, list of flipped (i, j)). A pair is tested when W4 separates it by at least `res`."""
    pairs, bad = 0, []
    for i in range(len(w4)):
        for j in range(i + 1, len(w4)):
            d4 = w4[j] - w4[i]
            if abs(d4) < res - 1e-9:
                continue
            pairs += 1
            if not ((w2[j] - w2[i]) * d4 > 0):
                bad.append((i, j))
    return pairs, len(bad), bad


RULE_TEXT = ("R1 order: every pair of settings that W4 separates by 0.025 (volumes), 2 points (peaks), 3 days or 0.03 (tail) keeps its order in W2.  "
             "R2 spread: the W2-minus-W4 offset varies over the five settings by at most 0.04 (volumes) and 4 points (peaks).  "
             "R3: a mean offset above 0.02 / 2 points is read as a flagged W4-equivalent correction.  "
             "PASS = R1 and R2; ORDER ONLY = R1 only; FAIL = R1 broken.")


def proxy_verdict(man: dict, lab: Path, csv_path: str | None = None) -> dict:
    ref, notes = w4_reference(lab, csv_path)
    csv_scores = T.load_scores(lab)
    w2, missing = {}, []
    for t in PROXY_SHOW_ORDER:
        m = key_metrics(man, lab, JOBS[t]["run_id"], csv_scores)
        if m is None:
            missing.append(t)
        else:
            w2[t] = m
    out = dict(status="INCOMPLETE", missing=missing, notes=notes, rows=[], ref=ref, w2=w2, corrections={}, reasons=[])
    if missing:
        out["reasons"].append(f"{len(PROXY_SHOW_ORDER) - len(missing)} of 5 proxy runs are scored; not finished: {', '.join(missing)}")
        return out
    tags = PROXY_SHOW_ORDER
    for key, label, res, tol, flat_tol, gated in KEY_METRICS:
        a4, a2 = [ref[t][key] for t in tags], [w2[t][key] for t in tags]
        pairs, flips, bad = order_test(a4, a2, res)
        offs = [y - x for x, y in zip(a4, a2)]
        spread, mean = max(offs) - min(offs), sum(offs) / len(offs)
        out["rows"].append(dict(key=key, label=label, pairs=pairs, flips=flips, bad=[(tags[i], tags[j]) for i, j in bad], offsets=offs,
                                spread=spread, tol=tol, spread_ok=spread <= tol + 1e-9, gated=gated, mean=mean, flat_tol=flat_tol,
                                flat=abs(mean) <= flat_tol + 1e-9))
    broke = [r for r in out["rows"] if r["flips"]]
    wide = [r for r in out["rows"] if r["gated"] and not r["spread_ok"]]
    if broke:
        out["status"] = "FAIL"
        out["reasons"] += [f"{r['label']}: {r['flips']} of {r['pairs']} pairs flipped ({'; '.join(a + ' vs ' + b for a, b in r['bad'])})" for r in broke]
    elif wide:
        out["status"] = "ORDER ONLY"
        out["reasons"] += [f"{r['label']}: offset spread {r['spread']:.3g} is above {r['tol']:g}" for r in wide]
    else:
        out["status"] = "PASS"
        out["corrections"] = {r["key"]: r["mean"] for r in out["rows"] if r["gated"] and not r["flat"]}
    return out


CONSEQUENCE = {
    "PASS": "W2 ranks the five settings like W4 and its offset from W4 is steady enough to read values. The ks and grid steps may run. "
            "The best one or two settings are still confirmed on W4 before they are adopted.",
    "ORDER ONLY": "W2 ranks the settings like W4, but its offset from W4 changes from setting to setting, so its values cannot be read against the targets. "
                  "The ks and grid steps are NOT started. Use W2 only to rank settings; do the Ks re-tune on W4.",
    "FAIL": "W2 does not rank the settings the way W4 does. The ks and grid steps are NOT started. Do the Ks re-tune on W4.",
    "INCOMPLETE": "The rule cannot be applied yet. The ks and grid steps are NOT started until all five proxy runs are scored.",
}


def print_verdict(v: dict) -> None:
    T.log("  PROXY CHECK  (does W2 rank the five settings like W4? Rule written before the runs; see the top of the script)")
    for n in v["notes"]:
        T.log(f"    {n}")
    if v["rows"]:
        tags = PROXY_SHOW_ORDER
        T.log("    " + f"{'setting':<18}| " + " | ".join(f"{lab[:14]:^15}" for lab in ("flood-wk vol", "flood-wk pk", "Feb27-28 vol", "Feb27-28 pk", "days in band", "tail1 vol")))
        T.log("    " + f"{'':<18}| " + " | ".join(f"{'W4':>7}{'W2':>8}" for _ in KEY_METRICS))
        for t in tags:
            cells = []
            for key, _, _, _, _, _ in KEY_METRICS:
                nd = 0 if key == "days" else (1 if key.endswith("pk") else 3)
                cells.append(f"{T.fm(v['ref'][t][key], nd):>7}{T.fm(v['w2'][t][key], nd):>8}")
            T.log("    " + f"{PROXY_LABELS[t]:<18}| " + " | ".join(cells))
        T.log()
        T.log("    " + f"{'metric':<32}| {'pairs tested':>12} {'flipped':>8} | {'offset spread':>13} {'tolerance':>9} {'decides':>8} | {'mean offset':>11} {'constant?':>10}")
        for r in v["rows"]:
            nd = 0 if r["key"] == "days" else (1 if r["key"].endswith("pk") else 3)
            T.log("    " + f"{r['label']:<32}| {r['pairs']:>12} {r['flips']:>8} | {T.fm(r['spread'], max(nd, 2)):>13} {r['tol']:>9g} {('yes' if r['gated'] else 'no'):>8} | "
                           f"{T.fm(r['mean'], max(nd, 2), sign=True):>11} {('within ' + format(r['flat_tol'], 'g')) if r['flat'] else ('ABOVE ' + format(r['flat_tol'], 'g')):>10}")
    T.log()
    T.log(f"    VERDICT: {v['status']}")
    for rs in v["reasons"]:
        T.log(f"      - {rs}")
    T.log(f"    {CONSEQUENCE[v['status']]}")
    if v["status"] == "PASS" and v["corrections"]:
        T.log("    Constant offsets above the R3 tolerance (read these W2 numbers as W4-equivalent = W2 minus the offset; see the W4eq table): "
              + ", ".join(f"{k} {c:+.3g}" for k, c in v["corrections"].items()))


PROXY_LABELS = {"base": "baseline (As 1)", "f0p5874_as300": "f x0.5874 As 300", "f0p15_as300": "f x0.15 As 300",
                "f0p02_as300": "f x0.02 As 300", "f0p02_as500": "f x0.02 As 500"}


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
    gw = s.get("GW_FRACTION")
    if gw is not None and abs(gw - W2["gw"]) > 1e-9:
        parts.append(f"GW {gw:g}")
    bd = s.get("BEDROCK_ADD_M")
    if bd is not None and abs(bd - T.BASE["BEDROCK_ADD_M"]) > 1e-9:
        parts.append(f"bed +{bd:g}")
    return " ".join(parts)


def sorted_runs(man: dict) -> list:
    recs = [(rid, rec) for rid, rec in man["runs"].items() if rid.startswith(PREFIX)]
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
        members = [(1.0, JOBS[ref_tag])] + [(j["ks_factor"], j) for j in STEP_JOBS["ks"][2:] if j["f"] == fv and j["an"] == a]
        members.sort(key=lambda t: t[0])
        got = [(fac, j, key_metrics(man, lab, j["run_id"], csv_scores)) for fac, j in members]
        got = [t for t in got if t[2] is not None]
        if len(got) < 2 or not any(fac != 1.0 for fac, _, _ in got):
            continue
        if not printed:
            T.log()
            T.log("  Ks RESPONSE  (each group: the same f and As, only Ks changes; Ks x1.00 is the baseline Ks 0.187. delta = change from the Ks x1.00 run)")
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
        T.log("      (the exponents come from the runs listed, 2 to 4 points; the prior came from 169 runs at As 1. A large difference means the As 1 prior does not carry over.)")


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
    T.log(f"W2 RE-TUNE (Jan 18 - Mar 12 2017, {W2['hours']} h from a {W2['gw']:g} x bedrock water table; channel loss off). Baseline = draw 12 / k10 (Ks x{T.BASE['KS_MULT']:g}, f x{T.BASE['F_MULT']:g}, As 1).")
    T.log("  vol = model volume / gauge volume. pk = model peak vs gauge peak, %. arr = scorer arrival error, hours (negative = early; see the TIMING table for what it hides).")
    T.log(f"  storms_ok = your anchor rule (flood-week vol {T.ANCHOR_FW[0]:g}-{T.ANCHOR_FW[1]:g}, Feb 27-28 vol {T.ANCHOR_F27[0]:g}-{T.ANCHOR_F27[1]:g}). peaks_ok = flood-week and Feb 27-28 peak errors both within 15%.")
    v = proxy_verdict(man, lab, getattr(args, "w4_summary", None))
    T.log()
    print_verdict(v)
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

    if v["status"] == "PASS" and v["corrections"]:
        T.log()
        T.log("  W4eq  (W4-EQUIVALENT ESTIMATE = W2 value minus the mean W2-minus-W4 offset of the five proxy settings; only the metrics R3 flagged are changed, marked *; the uncertainty is about half the offset spread)")
        T.log("    " + f"{'run':<26}| " + " ".join(f"{k + ('*' if k in v['corrections'] else ''):>9}" for k, *_ in KEY_METRICS))
        for rid, rec, row in done:
            m = key_metrics(man, lab, rid, csv_scores)
            if m is None:
                continue
            T.log("    " + f"{label_of(rec):<26}| " + " ".join(f"{T.fm(m[k] - v['corrections'].get(k, 0.0), 0 if k == 'days' else (1 if k.endswith('pk') else 3)):>9}" for k, *_ in KEY_METRICS))

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
    write_summary_csv(man, lab, csv_scores, v)
    T.log()
    T.log(f"  Wrote {T.SUMMARY_CSV}. These tables do not decide anything except the proxy rule: whether a Ks step is worth its cost to the tail and the Feb 27-28 peak is your call.")


def write_summary_csv(man: dict, lab: Path, csv_scores: dict, verdict: dict) -> None:
    rows = []
    for rid, rec in sorted_runs(man):
        row_sc = scores_row(man, lab, rid, csv_scores)
        s = rec.get("settings") or {}
        r = dict(run_id=rid, tag=rid[len(PREFIX) + 1:-len("_nocl")] if rid.endswith("_nocl") else rid, step=rec.get("step"), status=rec.get("status"),
                 F_MULT=s.get("F_MULT"), ANISOTROPY=s.get("ANISOTROPY"), ANISOTROPY_U=s.get("ANISOTROPY_U"), KS_MULT=s.get("KS_MULT"),
                 ks_factor=(s.get("KS_MULT") / T.BASE["KS_MULT"]) if s.get("KS_MULT") else "", GW_FRACTION=s.get("GW_FRACTION"), BEDROCK_ADD_M=s.get("BEDROCK_ADD_M"),
                 storms_ok=(T.storms_ok(row_sc) if rec.get("status") == "done" else ""), peaks_ok=(peaks_ok(row_sc) if rec.get("status") == "done" else ""),
                 tail_ok=((rec.get("tail") or {}).get("tail_ok", "")), proxy_verdict=verdict["status"])
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
DEFAULT_PLOT = ["base", "f0p15_as300", "f0p15_as300_ksx1p06", "f0p02_as300"]
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
    ap = argparse.ArgumentParser(description="Small Ks re-tune at f x0.15 / As 300 on the short window W2, with a W2-vs-W4 proxy check first (see the top of this file).")
    ap.add_argument("--dry-run", action="store_true", help="check everything, list the plan, build nothing")
    ap.add_argument("--summary-only", action="store_true", help="print the tables from the runs already scored, then stop")
    ap.add_argument("--build-only", action="store_true", help="build and check the inputs but do not run tRIBS (use with --limit 1)")
    ap.add_argument("--limit", type=int, default=None, help="make only the first N runs that still need work (finished runs do not count)")
    ap.add_argument("--steps", default=DEFAULT_STEPS, help=f"comma-separated steps: proxy, ks, grid, phys (default {DEFAULT_STEPS}; the baseline always runs first)")
    ap.add_argument("--ignore-verdict", action="store_true", help="start the ks and grid steps even if the proxy check did not PASS (the override is logged)")
    ap.add_argument("--plot", action="store_true", help="make retune_w2_figs/fig_w2_overlays.png from the finished runs")
    ap.add_argument("--plot-runs", default=",".join(DEFAULT_PLOT), help="runs to draw, as the tag part of the run name (default %(default)s)")
    ap.add_argument("--w4-summary", default=None, help=f"path of the W4 summary CSV (default lab/{W4_SUMMARY}, if it exists)")
    ap.add_argument("--band-lo", type=float, default=T.BAND_LO, help="slow-flow band, low end, m3/s")
    ap.add_argument("--band-hi", type=float, default=T.BAND_HI, help="slow-flow band, high end, m3/s")
    ap.add_argument("--redo", default=None, help="comma-separated run ids to run again even if done (you almost never want this)")
    ap.add_argument("--timeout-min", type=float, default=None, help="stop a tRIBS run after this many minutes (default: 3.5 times the expected time, about 33 minutes for W2)")
    ap.add_argument("--cell-timeout", type=int, default=1800, help="seconds allowed for one notebook cell")
    ap.add_argument("--min-free-gb", type=float, default=3.0, help="stop if less than this much disk is free (default 3)")
    ap.add_argument("--heartbeat-min", type=float, default=T.HEARTBEAT_MIN, help="print a progress line every this many minutes while tRIBS runs (default 5)")
    ap.add_argument("--notebook", default=None, help="notebook to run (default: Make_Sycamore_Model.ipynb, else Make_Sycamore_Model_clean.ipynb)")
    ap.add_argument("--scorer", default=T.SCORER)
    ap.add_argument("--tribs", default=os.environ.get("TRIBS_BIN") or (T.TRIBS_DEFAULT if os.path.isfile(T.TRIBS_DEFAULT) else shutil.which("tRIBS") or T.TRIBS_DEFAULT))
    ap.add_argument("--skip-kernel-check", action="store_true", help="skip the pytRIBS import test in --dry-run")
    return ap.parse_args(argv)


def gate(step: str, man: dict, lab: Path, args) -> tuple:
    """(may the step start?, message). ks and grid need a PASS from the proxy check, unless --ignore-verdict."""
    if step not in GATED_STEPS:
        return True, ""
    v = proxy_verdict(man, lab, getattr(args, "w4_summary", None))
    man["proxy_verdict"] = dict(status=v["status"], reasons=v["reasons"], checked=T.now_iso(), ignored=bool(args.ignore_verdict and v["status"] != "PASS"))
    if v["status"] == "PASS":
        return True, "proxy check PASS"
    if args.ignore_verdict:
        return True, f"proxy check is {v['status']}; --ignore-verdict was given, so the {step} step starts anyway (read its numbers with that in mind)"
    return False, f"the {step} step was NOT started: proxy check is {v['status']}. " + CONSEQUENCE[v["status"]]


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
                T.log("  Figure: " + make_plot(man, lab, plot_tags, lab / FIG_DIR / "fig_w2_overlays.png"))
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
    per_lo, per_hi = W2["hours"] * T.MIN_PER_1000H / 1000.0, W2["hours"] * T.MIN_PER_1000H_SLOW / 1000.0
    T.log(f"Working folder: {lab}")
    T.log(f"Notebook: {args.notebook}")
    T.log(f"Plan: steps {','.join(steps)}, window W2 ({W2['hours']} h from {W2['start']}, water table {W2['gw']:g} x bedrock): {n_todo} run(s) to make"
          + (f", {n_done} already finished" if n_done else "") + ".")
    for j in preview:
        what = "already done, skipped" if T.is_finished(man, j, args) else "NEW RUN"
        T.log(f"  {j['run_id']:48s} {describe(j):52s} {what}")
    gated = [s for s in steps if s in GATED_STEPS]
    if gated:
        T.log(f"  (the {' and '.join(gated)} step{'s' if len(gated) > 1 else ''} start{'' if len(gated) > 1 else 's'} only if the proxy check PASSES, unless --ignore-verdict; a step that is not started is simply listed above and skipped)")
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
    man["script"] = SCRIPT_NAME
    man["baseline_id"] = T.BASELINE_ID
    man["last_started"] = T.now_iso()
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
            if st != "proxy":
                if args.build_only:
                    T.log(f"The {st} step is skipped with --build-only.")
                    continue
                if code != 0:
                    T.log(f"The {st} step was not started because an earlier run failed or was interrupted. Fix that, then run the same command again.")
                    continue
                if budget["left"] is not None and budget["left"] <= 0:
                    T.log(f"The {st} step was not started: --limit was used up. Run again without --limit (or with a larger one) to continue.")
                    continue
                if st in GATED_STEPS:
                    may, msg = gate(st, man, lab, args)
                    T.save_manifest(man)
                    T.log()
                    T.log((f"{st} step: {msg}") if may else msg)
                    if not may:
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
