#!/usr/bin/env python3
"""
run_w4_kchan_sycamore.py  --  Stage 5, Step 8: a small CHANNEL-LOSS PROBE. Draw 12 with the W1B-chosen routing (sweep point k10) on the
long window W4 (Oct 1 2016 - Apr 1 2017), with channel loss turned ON at two strengths: Kchan 1 and Kchan 2 mm/hr.

Why. The W4 run with the k10 routing (Stage5_Routing_Log.md, Section 8.8) fixed arrival time on the flood week and on Feb 27-28, but it left
four things untouched: the Oct-Dec false flow (3.71 million m3 against a gauge of 0), the Feb 18-19 overshoot (+312%), the second January
storm (volume ratio 0.05) and the recession tails (volume ratios 0.27 and 0.08). Routing cannot fix those. The one knob left that removes
water on its way to the outlet is the channel loss. The Stage 3 test on draw 99 (baseline routing) showed that a CONSTANT loss takes a
larger share of small fall storms than of the big winter storms, but not by enough to be a fix, and that it also erases small real events
and tails. This probe asks the same question for draw 12 at the k10 routing, with only two values, so it costs about 75 minutes:
    Kchan 1 mm/hr   (Stage 3: a 14-19% cut to the February storms on draw 99)
    Kchan 2 mm/hr   (Stage 3: a 24-32% cut to the February storms on draw 99)
It decides nothing. It shows, run by run and storm by storm, how much of each part of the record the loss removes, so we can see whether the
fall false flow shrinks faster than the real storms and tails do.

Which runs. TWO runs, in this order:  Kchan 1, then Kchan 2.  Both are anchor draw 12 (Ks 0.187, f 0.5874, PsiB 0.2939 x baseline), routing
cv (KINEMVELCOEF) 6.39, flowexp (FLOWEXP) 0.2973, n (CHANNELROUGHNESS) 0.067 (row k=10 of rt1_list.csv), channel width coefficient 2.33,
snow on (OPTSNOW 1), thetaS, m, anisotropy and bedrock at baseline (bedrock = SOLUS + 1 m), warm start (water table 0.95 x bedrock),
Oct 1 2016, 4,368 h. The ONLY differences from the loss-off run W4_ks0p187_cv6p39_psib0p2939_f0p5874_d012_fe0p2973_n0p067_nocl are
OPTPERCOLATION (0 -> 1) and CHANNELCONDUCTIVITY (5, unused while the loss is off -> 1 or 2 mm/hr). The new runs are named
W4_ks0p187_cv6p39_psib0p2939_f0p5874_d012_fe0p2973_n0p067_kc1  and  ..._kc2.

This script IMPORTS run_w4_draw12_sycamore.py (keep it in the same folder) and reuses its tested build, check, tRIBS and scoring functions.
It does NOT import run_w4_routing_sycamore.py. It changes only: (1) the settings of the planned runs (routing as above, OPTPERCOLATION 1,
CHANNELCONDUCTIVITY 1 or 2); (2) every output file has its own name (w4kc_*), so the w4d12_* and w4rt_* files are never touched; (3) the
run names; (4) the checks of the .in file (below); (5) the summary, which puts the two new runs NEXT TO the loss-off k10 run (read only) and
shows how much of each storm the loss removed.

Checks added here:
  * The routing values are read from row k=10 of rt1_list.csv (the frozen sweep list, md5 c1e799f82ff16d30d91ff377b63356b7) and must equal
    the values above. If the list was changed or row 10 differs, nothing starts.
  * The new .in file must show OPTPERCOLATION 1 and CHANNELCONDUCTIVITY equal to the planned Kchan, and match the W4 baseline's .in in
    every other keyword except the run-specific ones.
  * It is also compared with the .in of the loss-off k10 run (if it is in lab/): every keyword must be equal except OPTPERCOLATION,
    CHANNELCONDUCTIVITY and the output name. (If that file is missing, the comparison is made with the earlier baseline-routing draw-12
    .in, which also allows the three routing keywords to differ.)
  * After each run, the scorer's channel-loss volume (proc__channel_loss_Mm3, read from the model's .mrf file) must be above zero, and the
    outlet file must not be byte-identical to the loss-off k10 run's. If either fails, the batch STOPS (loss was not really on, so the next
    run would be wasted). The finished run is kept and scored.

NOTHING IS REPEATED OR OVERWRITTEN. Each run has a unique name, and a planned run whose name already has files in results/ stops the
start-up check. The loss-off k10 run, the earlier draw-12 run at baseline routing and the W4 baseline are only read.

Run it from workspaces/Sycamore/lab/, in this order:
    python -u run_w4_kchan_sycamore.py --dry-run               # checks only, builds nothing
    python -u run_w4_kchan_sycamore.py --build-only --limit 1  # builds Kchan 1 and checks it, no tRIBS (about 1 min)
    python -u run_w4_kchan_sycamore.py                         # both runs (about 36 min each, about 75 min in all), then paste the tables
    python -u run_w4_kchan_sycamore.py --summary-only          # just reprint the tables from what is finished

Options (as in run_w4_draw12_sycamore.py): --limit N  --redo ID  --timeout-min N  --heartbeat-min N  --notebook FILE
(--anchors exists but only 12 is accepted.)

Files it makes (all in lab/): w4kc_list.csv (the frozen list), w4kc_manifest.json, w2.lock (only while running), w4kc_logs/,
score_out_w4kc/scores_w4kc.csv, w4kc_summary.csv, and <run>.in, simulation_<run>.log and results/<run>* for each run.

Do NOT run the notebook (or any cell that rebuilds data/model/) while this runs: every build rewrites lab/data/model/. Keep the terminal
open: the progress lines also stop the Codespace from timing out for being idle. If the connection drops, run the same command again: a
finished run is skipped.
"""
from __future__ import annotations

import csv
import hashlib
import math
import sys
from pathlib import Path

import numpy as np

import run_w4_draw12_sycamore as base      # the tested W4 runner: its functions are reused here

# --------------------------------------------------------------------------------------------
# What this batch does
# --------------------------------------------------------------------------------------------
K_POINT = 10                                                       # sweep point 10 of rt1_list.csv (the W1B winner at both anchors)
ROUTING = {"KINEMVELCOEF": 6.39, "FLOWEXP": 0.2973, "CHANNELROUGHNESS": 0.067}
ROUTING_KEYS = tuple(ROUTING)
RT1_LIST = "rt1_list.csv"
RT1_LIST_MD5 = "c1e799f82ff16d30d91ff377b63356b7"                  # the frozen sweep list (Output, Nov 1)
KCHAN_VALUES = (1.0, 2.0)                                          # mm/hr, in run order
DRAW = 12
ONLY_ANCHORS = "12"
LOSS_KEYS = ("OPTPERCOLATION", "CHANNELCONDUCTIVITY")
GROUP_PREFIX = "k10_kchan_"
SCRIPT_NAME = "run_w4_kchan_sycamore.py"

OFF_SCORES = "score_out_w4rt/scores_w4rt.csv"                      # the loss-off k10 run (routing script's batch): read only
D12_SCORES = "score_out_w4d12/scores_w4d12.csv"                    # the earlier draw-12 run at baseline routing: read only
CRITERIA_ARRIVAL_H, CRITERIA_PEAK_PCT = 1.0, 15.0                  # the closing criteria (reported only)
OCTDEC_GOAL_MM3 = 1.25                                             # the design's goal for the Oct-Dec false flow (reported only)
C_LOSS = "proc__channel_loss_Mm3"

# Output names of THIS batch (the earlier batches' w4d12_* and w4rt_* names are never used here).
base.LOG_DIR = "w4kc_logs"
base.SCORE_DIR = "score_out_w4kc"
base.SCORE_CSV = "scores_w4kc.csv"
base.SUMMARY_CSV = "w4kc_summary.csv"
base.DRAWS_CSV = "w4kc_list.csv"
base.MANIFEST = "w4kc_manifest.json"

_orig_read_anchors = base.read_anchors
_orig_settings_for = base.settings_for
_orig_describe = base.describe
_orig_preflight = base.preflight
_orig_save_manifest = base.save_manifest
_orig_do_run = base.do_run


def log(msg: str = "") -> None:
    base.log(msg)


def file_md5(path: Path) -> str:
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


# --------------------------------------------------------------------------------------------
# Names and the routing source
# --------------------------------------------------------------------------------------------
def kc_tag(kchan: float) -> str:
    return base.vtag(kchan)


def group_for(kchan: float) -> str:
    return f"{GROUP_PREFIX}{kc_tag(kchan)}"


GROUP_KCHAN = {group_for(k): k for k in KCHAN_VALUES}


def soil_part(ks, f, psib, idx) -> str:
    v = base.vtag
    return (f"W4_ks{v(ks)}_cv{v(ROUTING['KINEMVELCOEF'])}_psib{v(psib)}_f{v(f)}_d{idx:03d}"
            f"_fe{v(ROUTING['FLOWEXP'])}_n{v(ROUTING['CHANNELROUGHNESS'])}")


def new_run_id(ks, f, psib, idx, kchan) -> str:
    return f"{soil_part(ks, f, psib, idx)}_kc{kc_tag(kchan)}"


def loss_off_id(job: dict) -> str:
    """The k10 run of the SAME soil with the loss off (the routing script's run): the main comparison, read only."""
    return f"{soil_part(job['ks'], job['f'], job['psib'], job['idx'])}_nocl"


def baseline_routing_id(job: dict) -> str:
    """The earlier W4 run of the SAME soil at baseline routing and loss off, e.g. W4_ks0p187_cv6_psib0p2939_f0p5874_d012_nocl."""
    return base.run_id_for(job["ks"], job["f"], job["psib"], job["idx"])


def kchan_of(job: dict) -> float:
    k = GROUP_KCHAN.get(job.get("group"))
    if k is None:
        raise base.RunError(f"{job.get('run_id')}: its group {job.get('group')!r} does not name a planned Kchan {sorted(GROUP_KCHAN)}.")
    return k


def check_routing_source(lab: Path) -> str:
    """Row K_POINT of the frozen sweep list must hold exactly the routing values above. Raises base.Abort."""
    p = lab / RT1_LIST
    if not p.is_file():
        raise base.Abort(f"{RT1_LIST} not found in lab/. It is the frozen routing list the W1B check was made from.")
    md5 = file_md5(p)
    if md5 != RT1_LIST_MD5:
        raise base.Abort(f"{RT1_LIST} has md5 {md5}, not the frozen {RT1_LIST_MD5}. The list must not change; restore it.")
    with open(p, newline="") as fh:
        rows = {int(r["k"]): r for r in csv.DictReader(fh)}
    r = rows.get(K_POINT)
    if r is None:
        raise base.Abort(f"{RT1_LIST} has no row k = {K_POINT}.")
    got = {"KINEMVELCOEF": float(r["cv"]), "FLOWEXP": float(r["flowexp"]), "CHANNELROUGHNESS": float(r["n"])}
    bad = [f"{k}: list has {got[k]!r}, this script has {ROUTING[k]!r}" for k in ROUTING if abs(got[k] - ROUTING[k]) > 1e-9]
    if bad:
        raise base.Abort(f"row k = {K_POINT} of {RT1_LIST} does not match the routing values in this script:\n  " + "\n  ".join(bad))
    return (f"routing = row k={K_POINT} of {RT1_LIST} (md5 {md5[:8]}...): cv {ROUTING['KINEMVELCOEF']:g}, flowexp {ROUTING['FLOWEXP']:g}, "
            f"n {ROUTING['CHANNELROUGHNESS']:g}")


# --------------------------------------------------------------------------------------------
# Overrides of the base functions (the base module looks these up by name at call time)
# --------------------------------------------------------------------------------------------
def read_anchors(lab: Path, anchors_arg: str) -> list:
    """Draw 12 is read and checked by the base (Stage 3 design table and handoff values), then used once per Kchan."""
    if anchors_arg.replace(" ", "") != ONLY_ANCHORS:
        raise base.Abort(f"--anchors {anchors_arg}: this probe was designed for draw {DRAW} only (--anchors {ONLY_ANCHORS}).")
    check_routing_source(lab)
    draws = _orig_read_anchors(lab, ONLY_ANCHORS)
    if len(draws) != 1 or draws[0]["idx"] != DRAW:
        raise base.Abort(f"expected exactly draw {DRAW} from the design table, got {[d['idx'] for d in draws]}.")
    d0 = draws[0]
    out = []
    for k in KCHAN_VALUES:
        d = dict(d0)
        d["run_id"] = new_run_id(d["ks"], d["f"], d["psib"], d["idx"], k)
        d["group"] = group_for(k)
        out.append(d)
    return out


def settings_for(job: dict) -> dict:
    """The base settings, with the k10 routing and the channel loss ON for the planned runs (the baseline reference job is unchanged)."""
    s = _orig_settings_for(job)
    if job.get("kind") == "run":
        s.update(ROUTING)
        s.update(OPTPERCOLATION=1, CHANNELCONDUCTIVITY=float(kchan_of(job)))
    return s


def describe(job: dict) -> str:
    text = _orig_describe(job)
    if job.get("kind") == "run":
        text += (f"; routing cv {ROUTING['KINEMVELCOEF']:g}, flowexp {ROUTING['FLOWEXP']:g}, n {ROUTING['CHANNELROUGHNESS']:g}; "
                 f"channel loss ON (OPTPERCOLATION 1), Kchan {kchan_of(job):g} mm/hr")
    return text


def save_manifest(man: dict) -> None:
    man["script"] = SCRIPT_NAME
    man["routing"] = dict(ROUTING, source=f"{RT1_LIST} row k={K_POINT}", rt1_list_md5=RT1_LIST_MD5)
    man["kchan_mm_hr"] = list(KCHAN_VALUES)
    _orig_save_manifest(man)


def check_in_file(new_in: Path, ref_in: Path, job: dict) -> None:
    """(1) Against the W4 baseline .in: equal in every keyword except the run-specific ones, with the routing keywords, OPTPERCOLATION and
    CHANNELCONDUCTIVITY at exactly the planned values. (2) Against the loss-off k10 run's .in, if it is in lab/: equal in every keyword
    except OPTPERCOLATION, CHANNELCONDUCTIVITY and the output name."""
    s = base.settings_for(job)
    gw_name, bed_name = base.expected_files(s)
    new, ref = base.read_in_file(new_in), base.read_in_file(ref_in)
    expect = {"RUNTIME": base.WIN["hours"], "STARTDATE": base.WIN["start"], "OUTFILENAME": f"results/{job['run_id']}",
              "GWATERFILE": f"{base.MODEL_SOIL}/{gw_name}", "BEDROCKFILE": f"{base.MODEL_SOIL}/{bed_name}"}
    for kw in base.ROUTING_KWS + list(LOSS_KEYS):
        expect[kw] = s[kw]
    problems = []
    for key in sorted(set(ref) | set(new)):
        b, n = ref.get(key), new.get(key)
        if key in expect:
            if n is None or not base.same_value(n, expect[key]):
                problems.append(f"{key}: written {n!r}, expected {expect[key]!r}")
        elif b is None or n is None or not base.same_value(b, n):
            problems.append(f"{key}: reference {b!r} but this run {n!r}  (this keyword may not change)")
    if problems:
        hint = ""
        if any(p.startswith("CHANNELCONDUCTIVITY") for p in problems):
            hint = ("\n  (If CHANNELCONDUCTIVITY is written in a different unit from the Settings value, send me the .in line and the build "
                    "output file; do not run tRIBS.)")
        raise base.RunError("The .in file differs from the reference in more than it should:\n  " + "\n  ".join(problems) + hint)

    off = Path(f"{loss_off_id(job)}.in")
    earlier = Path(f"{baseline_routing_id(job)}.in")
    if off.is_file():
        other, allowed, what = base.read_in_file(off), set(LOSS_KEYS) | {"OUTFILENAME"}, f"{off.name} (same soil and routing, loss off)"
    elif earlier.is_file():
        other, allowed = base.read_in_file(earlier), set(LOSS_KEYS) | set(ROUTING_KEYS) | {"OUTFILENAME"}
        what = f"{earlier.name} (same soil, baseline routing, loss off; {off.name} is not in lab/)"
    else:
        log(f"{job['run_id']}: note: neither {off.name} nor {earlier.name} is in lab/, so the same-soil .in comparison was skipped")
        return
    probs2 = []
    for key in sorted(set(new) | set(other)):
        if key in allowed:
            continue
        a, b = other.get(key), new.get(key)
        if a is None or b is None or not base.same_value(a, b):
            probs2.append(f"{key}: {what.split(' ')[0]} has {a!r} but this run {b!r}")
    if probs2:
        raise base.RunError(f"The .in file differs from {what} in more than it should:\n  " + "\n  ".join(probs2))
    log(f"{job['run_id']}: same-soil comparison OK: the .in matches {what.split(' (')[0]} in every keyword except OPTPERCOLATION, "
        f"CHANNELCONDUCTIVITY" + (", routing" if allowed & set(ROUTING_KEYS) else "") + " and the output name; "
        f"OPTPERCOLATION {new.get('OPTPERCOLATION')}, CHANNELCONDUCTIVITY {new.get('CHANNELCONDUCTIVITY')}")


def preflight(args, lab: Path, jobs: list, base_problems: list, nb_note: str, man: dict) -> list:
    out = _orig_preflight(args, lab, jobs, base_problems, nb_note, man)
    try:
        out.append(("OK", check_routing_source(lab)))
    except base.Abort as exc:
        out.append(("FAIL", str(exc)))
    offs = read_scores(lab / OFF_SCORES)
    seen = set()
    for j in jobs:
        rid = loss_off_id(j)
        if rid in seen:
            continue
        seen.add(rid)
        has_in = (lab / f"{rid}.in").is_file()
        has_q = (lab / "results" / f"{rid}_Outlet.qout").is_file()
        has_sc = rid in offs
        if has_in and has_q and has_sc:
            out.append(("OK", f"loss-off comparison run {rid}: .in, outlet file and scores are there (read only)"))
        else:
            missing = [n for n, ok in ((".in file", has_in), ("outlet file", has_q), (f"score row in {OFF_SCORES}", has_sc)) if not ok]
            out.append(("WARN", f"loss-off comparison run {rid} is missing: {', '.join(missing)}. The runs themselves are not affected; "
                                f"the tables will show '-' for it and the same-soil .in comparison may be skipped or weaker."))
    out.append(("OK", "channel loss will be ON in the planned runs: OPTPERCOLATION 1, Kchan "
                      + ", ".join(f"{kchan_of(j):g}" for j in jobs) + " mm/hr (one planned run each)"))
    return out


def do_run(job: dict, args, lab: Path, man: dict, ref_in: Path) -> None:
    """The base run, then two checks that the loss was really on. A failure stops the batch (the next run would be wasted)."""
    _orig_do_run(job, args, lab, man, ref_in)
    rid = job["run_id"]
    rec = man["runs"].get(rid) or {}
    if job.get("kind") != "run" or rec.get("status") != "done" or args.build_only:
        return
    scores = rec.get("scores") or {}
    loss = scores.get(C_LOSS)
    off_q = lab / "results" / f"{loss_off_id(job)}_Outlet.qout"
    this_md5 = (rec.get("outlet") or {}).get("md5")
    same_as_off = bool(this_md5) and off_q.is_file() and file_md5(off_q) == this_md5
    rec["channel_loss_Mm3"] = loss
    rec["identical_to_loss_off"] = same_as_off
    base.save_manifest(man)
    if loss is None:
        log(f"{rid}: WARNING: the scorer did not report {C_LOSS} (no .mrf file found?), so the channel-loss volume cannot be checked. "
            f"Check the .in file ({base.read_in_file(lab / (rid + '.in')).get('OPTPERCOLATION')!r} for OPTPERCOLATION) and the table.")
    else:
        log(f"{rid}: channel loss in the model: {loss:.3f} million m3 over the window (from the .mrf file)")
    if same_as_off or (loss is not None and not loss > 0):
        why = []
        if same_as_off:
            why.append("the outlet file is byte-identical to the loss-off k10 run's")
        if loss is not None and not loss > 0:
            why.append(f"the model reports {loss:g} million m3 of channel loss")
        raise base.Abort(f"{rid}: channel loss does not look ON ({'; '.join(why)}). The run is kept and scored, but the batch stops here so no "
                         f"more hours are spent. Send me the .in file and this log.")


# --------------------------------------------------------------------------------------------
# The tables
# --------------------------------------------------------------------------------------------
def read_scores(path: Path) -> dict:
    out = {}
    if Path(path).is_file():
        with open(path, newline="") as fh:
            for r in csv.DictReader(fh):
                out[r["run_id"]] = r
    return out


def f_(x, nd=2, sign=False) -> str:
    return base.fm(x, nd, sign)


def n_(sc, col) -> float:
    return base.fnum(sc, col)


def zero_share(lab: Path, rid: str) -> float:
    """Share (percent) of the outlet file's rows where the model flow is exactly zero (read only; nan if unreadable)."""
    p = lab / "results" / f"{rid}_Outlet.qout"
    try:
        q = []
        with open(p) as fh:
            rd = csv.reader(fh)
            next(rd, None)
            for r in rd:
                if len(r) >= 2:
                    try:
                        q.append(float(r[1]))
                    except ValueError:
                        continue
        q = np.array(q)
        return float(100.0 * np.mean(q <= 1e-9)) if len(q) else float("nan")
    except Exception:
        return float("nan")


def pct_txt(x) -> str:
    """'30%' for a cut, '-5%' when the loss made the value larger, '-' when it cannot be computed."""
    return "-" if not math.isfinite(x) else f"{x:.0f}%"


def removed(off, kc) -> float:
    """Percent of the loss-off value that the loss removed: 100 * (1 - kc / off). nan if off is not above zero."""
    if off is None or kc is None or not (math.isfinite(off) and math.isfinite(kc)) or off <= 0:
        return float("nan")
    return 100.0 * (1.0 - kc / off)


STORMS = [("flood_week", "flood week Jan 18-25"), ("jan_second_storm", "second Jan storm"), ("feb18_19", "Feb 18-19"),
          ("feb27_28", "Feb 27-28")]
TAILS = [("tail_jan26_feb2", "tail Jan 26-Feb 2"), ("tail_mar3_mar11", "tail Mar 3-11")]


def print_summary(man: dict, lab: Path, draws: list) -> None:
    mine = base.load_scores(lab)
    rt = read_scores(lab / OFF_SCORES)
    d12 = read_scores(lab / D12_SCORES)
    get = lambda rid: mine.get(rid) or rt.get(rid) or d12.get(rid)
    log()
    log(f"W4 CHANNEL-LOSS PROBE (Oct 1 2016 - Apr 1 2017). Draw {DRAW}, k10 routing (cv {ROUTING['KINEMVELCOEF']:g}, flowexp {ROUTING['FLOWEXP']:g}, "
        f"n {ROUTING['CHANNELROUGHNESS']:g}), loss ON at Kchan " + " and ".join(f"{k:g}" for k in KCHAN_VALUES) + " mm/hr. Warm start; snow on.")
    log("  Oct-Dec = model volume Oct 1 - Dec 31 in million m3 (gauge: 0). days = days the model flow is above 0.13 m3/s; max = largest Oct-Dec flow (m3/s).")
    log("  Oct-Mar = whole-window volume ratio (model / gauge); KGE = KGE_2012 over the whole span (not a pass/fail here).")
    log("  loss = channel loss in the model (million m3, from the .mrf file). zero% = share of the 4,368 outlet rows where the model flow is exactly 0.")
    log("  Nothing is promoted or decided here: the tables show what the loss does to each part of the record.")

    ref = get(base.BASELINE_ID)
    if ref is not None:
        log(f"  Gauge, Oct-Dec: volume {f_(n_(ref, base.C_FFG), 3)} M m3, largest flow {f_(n_(ref, base.C_FFGMAX), 2)} m3/s, "
            f"{f_(n_(ref, base.C_FFGDAYS), 0)} days above {base.FF_THRESHOLD:g}. Design goal for the model's Oct-Dec flow: at most {OCTDEC_GOAL_MM3:g} M m3 (reported only).")
    log()
    log("  TABLE 1  false flow and whole window")
    log(f"  {'run':<46} {'OctDec':>7} {'days':>5} {'max':>6} {'Oct-Mar':>8} {'KGE':>6} {'loss':>6} {'zero%':>6}")

    def line1(label, sc, rid):
        if sc is None:
            return f"  {label:<46} (no scores)"
        z = zero_share(lab, rid)
        return (f"  {label:<46} {f_(n_(sc, base.C_FFV), 2):>7} {f_(n_(sc, base.C_FFDAYS), 0):>5} {f_(n_(sc, base.C_FFMAX), 1):>6} "
                f"{f_(n_(sc, base.C_OM), 3):>8} {f_(n_(sc, 'span__kge2012'), 3):>6} {f_(n_(sc, C_LOSS), 2):>6} {f_(z, 1):>6}")

    if ref is not None:
        log(line1("W4 baseline (Ks 0.2, f 1, PsiB 0.5; cv 6)", ref, base.BASELINE_ID))

    rows, finished = [], []
    for d in draws:
        rid, kc = d["run_id"], GROUP_KCHAN.get(d["group"], float("nan"))
        job = dict(ks=d["ks"], f=d["f"], psib=d["psib"], idx=d["idx"])
        off_id, bl_id = loss_off_id(job), baseline_routing_id(job)
        rec = man["runs"].get(rid) or {}
        status = rec.get("status", "not run")
        if not finished and not rows:
            log(line1(f"draw {d['idx']}, baseline routing, loss off", get(bl_id), bl_id))
            log(line1(f"draw {d['idx']}, k10 routing, loss off", get(off_id), off_id))
        sc = mine.get(rid) if status in ("done", "unscored") else None
        if sc is None:
            log(f"  {'draw ' + str(d['idx']) + ', k10 routing, Kchan ' + format(kc, 'g'):<46} ({status})")
        else:
            log(line1(f"draw {d['idx']}, k10 routing, Kchan {kc:g} mm/hr", sc, rid))
            finished.append((d, kc, sc, off_id))
        row = dict(idx=d["idx"], run_id=rid, group=d["group"], KS_MULT=d["ks"], F_MULT=d["f"], PSIB_MULT=d["psib"],
                   cv=ROUTING["KINEMVELCOEF"], flowexp=ROUTING["FLOWEXP"], n=ROUTING["CHANNELROUGHNESS"], kchan_mmhr=kc,
                   status=status, loss_off_run_id=off_id)
        if sc is not None:
            for k, v in sc.items():
                if k not in ("run_id", "window"):
                    row[k] = v
            row["zero_share_pct"] = zero_share(lab, rid)
            row["build_min"], row["tribs_min"], row["output_mb"] = rec.get("build_min", ""), rec.get("tribs_min", ""), rec.get("output_mb", "")
        rows.append(row)

    if finished:
        off_sc = get(finished[0][3])
        cols_runs = [("loss off", off_sc)] + [(f"Kchan {kc:g}", sc) for _, kc, sc, _ in finished]

        log()
        log("  TABLE 2  by storm   (vol = volume ratio, model / gauge; pk = peak error %; arr = first-arrival error h, negative = early)")
        head1 = f"  {'':<22}" + "".join(f" | {lab_:^21}" for lab_, _ in cols_runs)
        head2 = f"  {'storm':<22}" + "".join(f" | {'vol':>6} {'pk %':>7} {'arr h':>6}" for _ in cols_runs)
        log(head1)
        log(head2)

        def cells(s, key):
            if s is None:
                return f"{'-':>6} {'-':>7} {'-':>6}"
            return (f"{f_(n_(s, key + '__vol_ratio'), 3):>6} {f_(n_(s, key + '__peak_err_pct'), 1, True):>7} "
                    f"{f_(n_(s, key + '__arrival_err_h'), 2, True):>6}")
        for key, label in STORMS:
            log(f"  {label:<22}" + "".join(f" | {cells(s, key)}" for _, s in cols_runs))
        for key, label in TAILS:
            parts = [f"{lab_} {f_(n_(s, key + '__tail_vol_ratio'), 3)} ({f_(n_(s, key + '__recession_rate_ratio'), 2)})" for lab_, s in cols_runs]
            log(f"  {label:<22} tail volume ratio (recession-rate ratio): " + ";  ".join(parts))
        mo = []
        for tag, label in (("2017_01", "Jan"), ("2017_02", "Feb"), ("2017_03", "Mar")):
            mo.append(f"{label} " + " -> ".join(f_(n_(s, f"month_{tag}__vol_ratio"), 3) for _, s in cols_runs))
        log("  monthly volume ratio (" + " -> ".join(lab_ for lab_, _ in cols_runs) + "): " + "; ".join(mo))

        log()
        log("  TABLE 3  how much of the loss-off volume the loss removed, in percent (100 x (1 - Kchan run / loss-off run))")
        log("           For the loss to help, the Oct-Dec (false flow) row has to be clearly larger than the rows for real flow below it.")
        log(f"  {'part of the record':<36}" + "".join(f" {'Kchan ' + format(kc, 'g'):>10}" for _, kc, _, _ in finished) + f"   {'loss-off value':>14}")

        def rem_row(label, col, nd_off=3):
            off_v = n_(off_sc, col) if off_sc else float("nan")
            vals = "".join(f" {pct_txt(removed(off_v, n_(sc, col))):>10}" for _, _, sc, _ in finished)
            return f"  {label:<36}{vals}   {f_(off_v, nd_off):>14}"
        log(rem_row("Oct-Dec false flow (M m3)", base.C_FFV, 2))
        for key, label in STORMS:
            log(rem_row(label + " (vol ratio)", f"{key}__vol_ratio"))
        for key, label in TAILS:
            log(rem_row(label + " (tail ratio)", f"{key}__tail_vol_ratio"))
        log(rem_row("Oct-Mar whole window (ratio)", base.C_OM))
        log("           Rows with a small loss-off value (the second Jan storm, the tails) swing a lot in percent from a small amount of water.")
        log("           For reference, Stage 3 (draw 99, BASELINE routing, different soil), loss off -> Kchan 1 / Kchan 2: the fall storms lost 38% / 57%,")
        log("           Feb 18-19 lost 19% / 32% and Feb 27-28 lost 14% / 24%; small real events and the tails were erased at Kchan 0.5 and above.")

        log()
        log(f"  CLOSING CRITERIA (arrival within {CRITERIA_ARRIVAL_H:g} h, peak within {CRITERIA_PEAK_PCT:g}%): reported only")
        for lab_, s in cols_runs:
            parts = []
            for key, label in (("flood_week", "flood week"), ("feb27_28", "Feb 27-28")):
                arr, pk = n_(s, f"{key}__arrival_err_h"), n_(s, f"{key}__peak_err_pct")
                ok_a = math.isfinite(arr) and abs(arr) <= CRITERIA_ARRIVAL_H
                ok_p = math.isfinite(pk) and abs(pk) <= CRITERIA_PEAK_PCT
                parts.append(f"{label}: arrival {f_(arr, 2, True)} h {'met' if ok_a else 'NOT met'}, peak {f_(pk, 1, True)}% {'met' if ok_p else 'NOT met'}")
            log(f"    {lab_:<10} " + ";  ".join(parts))
        log()
        for lab_, s in cols_runs:
            v = n_(s, base.C_FFV)
            log(f"  Oct-Dec design goal (at most {OCTDEC_GOAL_MM3:g} M m3), {lab_}: {f_(v, 2)} M m3 -> {'met' if math.isfinite(v) and v <= OCTDEC_GOAL_MM3 else 'NOT met'}")

    log()
    log(f"  {len(finished)} of {len(draws)} run(s) scored; {sum(1 for r in rows if r['status'] == 'failed')} failed; "
        f"{sum(1 for r in rows if r['status'] == 'not run')} not run yet.")
    rejected = [rid for rid, rec in man["runs"].items() if rec.get("status") == "reuse_rejected"]
    if rejected:
        log("  NOT used (failed their checks) and NOT rerun: " + ", ".join(rejected))
    if rows:
        cols = []
        for r in rows:
            for c in r:
                if c not in cols:
                    cols.append(c)
        with open(lab / base.SUMMARY_CSV, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            for r in rows:
                w.writerow({c: r.get(c, "") for c in cols})
        log()
        log(f"  Wrote {base.SUMMARY_CSV} (all score columns for the new run(s)).")


base.read_anchors = read_anchors
base.settings_for = settings_for
base.describe = describe
base.save_manifest = save_manifest
base.check_in_file = check_in_file
base.preflight = preflight
base.do_run = do_run
base.print_summary = print_summary


def main(argv=None) -> int:
    return base.main(argv)


if __name__ == "__main__":
    sys.exit(main())
