#!/usr/bin/env python3
"""
run_kchan_w4_sycamore.py  --  Channel-loss check: one soil set (a draw already run on W4 with loss OFF), channel loss ON at a few
channel conductivity (Kchan) values, on W4 (Oct 1 2016 - Apr 1 2017).

Why. Every Sycamore run so far has had channel loss off (OPTPERCOLATION 0, "nocl"). On W4 the soil dial (Ks, f, PsiB) cannot remove
the Oct-Dec false flow (the model makes flow in the fall; the gauge records none) without also cutting the January flood. Channel
loss is water sinking into the streambed as it travels downstream. Small flows lose the largest share of their water and big floods
lose a smaller share, so loss could remove the false flow while leaving the flood mostly alone. This batch asks whether it does.

What changes. ONLY two settings: OPTPERCOLATION goes from 0 to 1, and CHANNELCONDUCTIVITY (Kchan, mm/hr) takes each value asked for
(default 0.5, 1 and 2). Ks, f and PsiB stay at the chosen draw's values (default draw 99). Everything else stays at the baseline:
thetaS, m, anisotropy, bedrock +1 m, width coefficient, the routing knobs (cv 6, flowexp 0.3, n 0.04), snow on, and the warm start
(water table 0.95 x bedrock). Soil is NOT re-tuned here: loss only removes water, so the whole-window volume will fall below the
loss-off run. That is expected. The question is how the false flow and the flood change relative to each other.

The reference. The loss-off W4 run of the same draw (for draw 99: W4_ks0p2524_cv6_psib0p5743_f1p165_d099_nocl). It is checked
read-only (it must be marked done in w4_manifest.json, its outlet file must match the md5 and checksum recorded then, and its .in
file must show loss off) and only scored into this batch's own score file. It is never rerun or changed.

Run names: W4_ks0p2524_cv6_psib0p5743_f1p165_d099_cc0p5 (cc0p5 = channel loss on, Kchan 0.5 mm/hr).

The list is frozen in kchan_w4_list.csv on the first start. A later start may ADD Kchan values (for example --kchan 0.5,1,2,5,10
after the first three are done), but it stops if the draw or an earlier value would change, so the batch never silently changes.
Values above 20 mm/hr need --allow-high (a guard against typing 50 for 0.5; at Sycamore a much larger Kchan dried every flow).

NOTHING IS REPEATED OR OVERWRITTEN. Every run has a unique name, and a planned run whose name already has files in results/ stops
the start-up check.

For each run the script:
    1. Copies the Settings cell of the Make_Sycamore_Model notebook, changes ONLY OPTPERCOLATION and CHANNELCONDUCTIVITY (plus the
       run names, the window, and the soil values of the draw, which are the ones the notebook already needs), and runs the whole
       notebook headlessly (no browser) with nbclient. The notebook file is never changed.
    2. Checks the build: the new .in file must match the loss-off reference .in in every keyword except the run-specific ones and
       the two loss keywords, and every soil grid in data/model/soil must equal its SOLUS grid times the multiplier.
    3. Runs tRIBS (one at a time, in the foreground) and checks that it finished.
    4. Scores the run with score_run_sycamore.py (window W4) and flags a run whose outlet file is byte-identical to the reference.
    5. Writes everything to kchan_w4_manifest.json, so a stopped batch resumes by running the same command.

READING RULE (set before any result is seen). A Kchan value helps only if it cuts the false flow by a clearly larger percent than it
cuts the flood-week volume. The tables print that "cut ratio" for each run, and, for scale, the same ratio for the soil dial alone
(each earlier W4 run against the W4 baseline). Kchan 0.5 should look almost the same as the loss-off run; if it does not, turning
loss on changes something besides the loss itself, and that comes first. Nothing is selected or dropped by the script.

Run it from workspaces/Sycamore/lab/. Typical use, in this order:
    python -u run_kchan_w4_sycamore.py --dry-run               # checks only, builds nothing
    python -u run_kchan_w4_sycamore.py --limit 1 --build-only  # builds the first run (Kchan 0.5) and checks it, no tRIBS
    python -u run_kchan_w4_sycamore.py --limit 1               # one full run (about 35 min), then look at the output
    python -u run_kchan_w4_sycamore.py --prune-soil            # everything left; run it again to resume after any stop
    python -u run_kchan_w4_sycamore.py --kchan 0.5,1,2,5,10 --prune-soil   # later: add 5 and 10 (only if the first three look promising)
    python -u run_kchan_w4_sycamore.py --summary-only          # just print the tables from what is finished

Options: --draw 99  --kchan 0.5,1,2  --allow-high  --limit N  --redo ID  --timeout-min N  --heartbeat-min N  --notebook FILE

Files it makes (all in lab/): kchan_w4_list.csv (the frozen list), kchan_w4_manifest.json, kchan_w4.lock (only while running),
kchan_w4_logs/, score_out_kchan_w4/scores_kchan_w4.csv, kchan_w4_summary.csv, and for each run <run>.in, simulation_<run>.log and
results/<run>*.

Do NOT run the notebook (or any cell that rebuilds data/model/) while this runs: every build rewrites lab/data/model/ and
all the runs share it. Keep the terminal open: the progress lines also stop the Codespace from timing out for being idle.
If the connection drops, run the same command again: finished runs are skipped.
"""
from __future__ import annotations

import argparse
import ast
import atexit
import csv
import datetime as dt
import filecmp
import hashlib
import json
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------------------------
# What the batch does
# --------------------------------------------------------------------------------------------
WARM_GW = 0.95     # W4 starts Oct 1: water table at 95% of bedrock depth (the start-state rule from Stage 1)
WIN = dict(name="W4", start="10/01/2016/00/00", hours=4368, label="oct_to_mar", gw=WARM_GW)   # Oct 1 2016 to Apr 1 2017; label = the RUN_CHOICE text printed (as in run_proxy_sycamore.py)

# Everything that is not drawn (the baseline, with the warm-start water table).
BASE = {
    "GW_FRACTION": WIN["gw"], "BEDROCK_ADD_M": 1.0, "KS_MULT": 0.2, "PSIB_MULT": 0.5, "F_MULT": 1,
    "THETAS_MULT": 1, "ANISOTROPY": 1, "M_MULT": 1,
    "KINEMVELCOEF": 6, "FLOWEXP": 0.3, "CHANNELROUGHNESS": 0.04, "CHANNELWIDTHCOEFF": 2.33,
    "OPTSNOW": 1, "OPTPERCOLATION": 0, "CHANNELCONDUCTIVITY": 5,
    "REBUILD_SOIL": False,
}
SIG_FIGS = 4
ROUTING_KWS = ["KINEMVELCOEF", "FLOWEXP", "CHANNELROUGHNESS", "CHANNELWIDTHCOEFF"]   # never changed here; must stay at baseline

W4_BASELINE_ID = "Sycamore_octmar_ks0p2_cv6_psib0p5_nocl"   # the W4 baseline (Ks 0.2, f 1, PsiB 0.5): used only as a context row and for the soil-dial scale
BASELINE_ID = ""      # the REFERENCE run: the loss-off W4 run of the chosen draw. Set in main() once the draw has been read.
ANCHOR = {}           # the chosen draw (idx, ks, f, psib, ref_id). Set in main().
PROXY_MANIFEST = "proxy_manifest.json"
NOTEBOOK_CANDIDATES = ["Make_Sycamore_Model.ipynb", "Make_Sycamore_Model_clean.ipynb"]
SCORER = "score_run_sycamore.py"
TRIBS_DEFAULT = "/opt/tRIBS/build/tRIBS"
INIT_DATA = "../init_data"           # SOLUS_DIR in the notebook
MODEL_SOIL = "data/model/soil"
LOG_DIR = "kchan_w4_logs"
SCORE_DIR = "score_out_kchan_w4"
SCORE_CSV = "scores_kchan_w4.csv"
SUMMARY_CSV = "kchan_w4_summary.csv"
DRAWS_CSV = "kchan_w4_list.csv"                        # the frozen list of runs
MANIFEST = "kchan_w4_manifest.json"
LOCK = "kchan_w4.lock"
DRAW_DEFAULT = 99                                      # the soil set to test (found in w4_list.csv); draw 99 is the wettest of the four anchors
KCHAN_DEFAULT = "0.5,1,2"                              # mm/hr, in run order
KCHAN_SOFT_MAX = 20.0                                  # mm/hr; anything higher needs --allow-high
W4_LIST_CSV = "w4_list.csv"                            # the anchors that finished on W4 (frozen by run_w4_anchors_sycamore.py)
W4_MANIFEST = "w4_manifest.json"                       # where those W4 runs are recorded as done
W4_SCORES_CSV = "score_out_w4/scores_w4.csv"           # read only: the W4 baseline and anchors, for the "soil dial alone" scale
W1_DESIGN_CSV = "lhs_volume_list.csv"                  # the Stage 3 design (the values the draw must match)
MIN_PER_1000H = 7.8          # measured tRIBS speed on this machine: minutes per 1,000 model hours
HEARTBEAT_MIN = 5.0
ROWS_PER_MODEL_HOUR = 16     # the outlet file gets 16 rows per model hour
F_FLOOR = 1e-4               # the notebook keeps f at or above this
THETAS_CAP = 0.99            # the notebook keeps thetaS and porosity below this

# Score columns shown in the tables (all are columns of the scorer's CSV).
C_FFV, C_FFG = "octdec__vol_model_Mm3", "octdec__vol_gauge_Mm3"            # Oct-Dec flow volume, M m3 (model, gauge)
C_FFMAX, C_FFGMAX = "octdec__max_model_m3s", "octdec__max_gauge_m3s"      # largest Oct-Dec flow, m3/s
C_FFDAYS, C_FFGDAYS = "octdec__days_above_model", "octdec__days_above_gauge"   # days above the 0.13 m3/s threshold
C_VOL, C_PK, C_ARR = "flood_week__vol_ratio", "flood_week__peak_err_pct", "flood_week__arrival_err_h"
C_OM, C_SPAN = "oct_mar__vol_ratio", "span__vol_ratio"
C_F18, C_F27 = "feb18_19__vol_ratio", "feb27_28__vol_ratio"
C_T1, C_T2 = "tail_jan26_feb2__tail_vol_ratio", "tail_mar3_mar11__tail_vol_ratio"
FF_THRESHOLD = 0.13          # m3/s; the design's false-flow threshold (the gauge's Oct-Dec maximum)

_LOGFILE = None


def log(msg: str = "") -> None:
    stamp = dt.datetime.now().strftime("%H:%M:%S")
    line = f"[{stamp}] {msg}" if msg else ""
    print(line, flush=True)
    if _LOGFILE is not None:
        try:
            with open(_LOGFILE, "a") as fh:
                fh.write(line + "\n")
        except OSError:
            pass


def now_iso() -> str:
    return dt.datetime.now().isoformat(timespec="seconds")


def fmt_min(seconds: float) -> str:
    return f"{seconds / 60:.1f} min"


def fmt_hm(seconds: float) -> str:
    m = int(round(seconds / 60))
    return f"{m // 60} h {m % 60:02d} min" if m >= 60 else f"{m} min"


class RunError(Exception):
    """A problem with one run. The batch records it and moves on to the next run."""


class Abort(Exception):
    """A problem with the whole batch (competing job, no disk, changed design). The batch stops."""


# --------------------------------------------------------------------------------------------
# The list of runs (the anchors, taken from w2_list.csv and frozen in w4_list.csv)
# --------------------------------------------------------------------------------------------
def vtag(x) -> str:
    """0.18 -> '0p18', 1.5 -> '1p5', 15 -> '15', 0.2374 -> '0p2374'."""
    return f"{x:g}".replace(".", "p")


def run_id_for(ks, f, psib, idx, kchan=None) -> str:
    """W4_ks0p2524_cv6_psib0p5743_f1p165_d099_cc0p5: window, Ks, cv, PsiB, then the other changes (f, draw number), channel setting last.
    nocl = channel loss off; cc<N> = channel loss on with Kchan N mm/hr."""
    tail = "nocl" if kchan is None else f"cc{vtag(kchan)}"
    return f"{WIN['name']}_ks{vtag(ks)}_cv6_psib{vtag(psib)}_f{vtag(f)}_d{idx:03d}_{tail}"


def parse_kchan(arg: str, allow_high: bool) -> list:
    """--kchan 0.5,1,2 -> [0.5, 1.0, 2.0]. Raises Abort on anything that is not a list of distinct positive numbers."""
    try:
        vals = [float(t) for t in arg.replace(" ", "").split(",") if t]
    except ValueError:
        raise Abort(f"--kchan must be numbers (mm/hr) separated by commas (for example 0.5,1,2); got {arg!r}.")
    if not vals:
        raise Abort("--kchan is empty.")
    if any((not math.isfinite(v)) or v <= 0 for v in vals):
        raise Abort(f"--kchan values must be positive numbers (mm/hr); got {arg!r}. (Use the loss-off run for 'no loss'.)")
    if len({vtag(v) for v in vals}) != len(vals):
        raise Abort(f"--kchan lists a value twice: {arg}")
    high = [v for v in vals if v > KCHAN_SOFT_MAX]
    if high and not allow_high:
        raise Abort(f"--kchan has value(s) above {KCHAN_SOFT_MAX:g} mm/hr ({', '.join(f'{v:g}' for v in high)}). That is much higher than the plan, and a large "
                    f"Kchan dried every flow at Sycamore. If you really mean it, add --allow-high.")
    return vals


def read_anchor(lab: Path, draw: int) -> dict:
    """Look up the chosen draw in w4_list.csv (it must already have run on W4 with loss off) and check it against the Stage 3 design."""
    p = lab / W4_LIST_CSV
    if not p.is_file():
        raise Abort(f"{p} not found. It is the list of anchor runs that finished on W4 (made by run_w4_anchors_sycamore.py).")
    with open(p, newline="") as fh:
        rows = {int(r["idx"]): r for r in csv.DictReader(fh)}
    r = rows.get(draw)
    if r is None:
        raise Abort(f"draw {draw} is not in {W4_LIST_CSV} (it holds draws {sorted(rows)}). The reference has to be a draw that already ran on W4 with loss off.")
    ks, f, psib = float(r["KS_MULT"]), float(r["F_MULT"]), float(r["PSIB_MULT"])
    ref_id = run_id_for(ks, f, psib, draw)
    problems = []
    if r["run_id"] != ref_id:
        problems.append(f"draw {draw}: the run id in {W4_LIST_CSV} ({r['run_id']}) does not match the name made from its values ({ref_id})")
    dp = lab / W1_DESIGN_CSV
    if dp.is_file():
        design = {}
        for dr in csv.DictReader(open(dp, newline="")):
            design[int(dr["idx"])] = (float(dr["KS_MULT"]), float(dr["F_MULT"]), float(dr["PSIB_MULT"]), dr["run_id"])
        want = design.get(draw)
        if want is None:
            problems.append(f"draw {draw} is not in {dp.name}")
        elif (want[0], want[1], want[2]) != (ks, f, psib) or want[3].replace("W1_", "W4_", 1) != ref_id:
            problems.append(f"draw {draw}: Ks/f/PsiB or run id differ from {dp.name} (draw {ks:g}, {f:g}, {psib:g}; design {want[0]:g}, {want[1]:g}, {want[2]:g})")
    if problems:
        raise Abort("the draw has problems; nothing was started:\n  " + "\n  ".join(problems))
    return dict(idx=draw, ks=ks, f=f, psib=psib, ref_id=ref_id)


def kchan_draws(anchor: dict, kchans: list) -> list:
    return [dict(idx=anchor["idx"], ks=anchor["ks"], f=anchor["f"], psib=anchor["psib"], kchan=float(k),
                 run_id=run_id_for(anchor["ks"], anchor["f"], anchor["psib"], anchor["idx"], k), group="kchan",
                 ref_run_id=anchor["ref_id"]) for k in kchans]


def write_draws(path: Path, draws: list) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["idx", "KS_MULT", "F_MULT", "PSIB_MULT", "KCHAN_MMHR", "run_id", "group", "ref_run_id"])
        for d in draws:
            w.writerow([d["idx"], d["ks"], d["f"], d["psib"], d["kchan"], d["run_id"], d["group"], d["ref_run_id"]])


def read_draws(path: Path) -> list:
    with open(path, newline="") as fh:
        return [dict(idx=int(r["idx"]), ks=float(r["KS_MULT"]), f=float(r["F_MULT"]), psib=float(r["PSIB_MULT"]), kchan=float(r["KCHAN_MMHR"]),
                     run_id=r["run_id"], group=r["group"], ref_run_id=r["ref_run_id"])
                for r in csv.DictReader(fh)]


def load_plan(lab: Path, draw: int | None, kchan_arg: str | None, allow_high: bool, write: bool = True) -> tuple:
    """Read the draw and the Kchan values, then read or make the frozen list. Returns (anchor, draws, note). Raises Abort if the frozen
    list would change. Adding Kchan values to the end of the list is allowed; changing the draw or dropping a value is not.
    With write=False (--dry-run, --summary-only) the list file is never created or changed.
    If --draw or --kchan is not given, the frozen list decides (so a restart, or --summary-only, needs no options); with no list yet, the defaults apply."""
    path = lab / DRAWS_CSV
    first = read_draws(path) if path.is_file() else None
    if draw is None:
        draw = first[0]["idx"] if first else DRAW_DEFAULT
    if kchan_arg is None:
        kchan_arg = ",".join("%g" % d["kchan"] for d in first) if first else KCHAN_DEFAULT
    anchor = read_anchor(lab, draw)
    wanted = kchan_draws(anchor, parse_kchan(kchan_arg, allow_high))
    if not path.is_file():
        if write:
            write_draws(path, wanted)
        return anchor, wanted, (f"draw {draw}, Kchan {kchan_arg}: checked against {W4_LIST_CSV} and {W1_DESIGN_CSV}, and "
                                + (f"froze the list in {DRAWS_CSV}" if write else f"would freeze the list in {DRAWS_CSV} (a real start does that)"))
    saved = read_draws(path)
    other = [d["run_id"] for d in saved if d["idx"] != draw or d["ref_run_id"] != anchor["ref_id"]
             or (d["ks"], d["f"], d["psib"]) != (anchor["ks"], anchor["f"], anchor["psib"])]
    if other:
        raise Abort(f"{DRAWS_CSV} is already there but is for a different draw or soil set than --draw {draw} now. Nothing was started. If you want a "
                    f"different set on purpose, rename the old {DRAWS_CSV} and {MANIFEST} first; do not mix two lists.")
    want_ids = {d["run_id"] for d in wanted}
    dropped = [d["run_id"] for d in saved if d["run_id"] not in want_ids]
    if dropped:
        raise Abort(f"{DRAWS_CSV} already holds {len(saved)} run(s), and --kchan {kchan_arg} no longer lists {len(dropped)} of them ({', '.join(dropped)}). "
                    f"Nothing was started. To add values, give all the earlier values plus the new ones (for example --kchan 0.5,1,2,5,10).")
    saved_ids = {d["run_id"] for d in saved}
    new = [d for d in wanted if d["run_id"] not in saved_ids]
    if new:
        merged = saved + new
        if write:
            write_draws(path, merged)
        new_txt = ", ".join("%g" % d["kchan"] for d in new)
        verb = "added" if write else "would add"
        return anchor, merged, (f"{verb} {len(new)} new Kchan value(s) ({new_txt}) to the frozen list in {DRAWS_CSV}; "
                                f"the earlier {len(saved)} are unchanged")
    return anchor, saved, f"read the frozen list from {DRAWS_CSV} ({len(saved)} runs; it matches --kchan {kchan_arg})"


# --------------------------------------------------------------------------------------------
# Jobs and settings
# --------------------------------------------------------------------------------------------
def make_job(d: dict | None) -> dict:
    """d = one Kchan run of the list. None gives the REFERENCE job: the loss-off run of the same draw (never rebuilt)."""
    if d is None:
        a = ANCHOR
        return dict(idx=a["idx"], ks=a["ks"], f=a["f"], psib=a["psib"], kchan=None, run_id=a["ref_id"], kind="baseline", group="reference")
    return dict(idx=d["idx"], ks=d["ks"], f=d["f"], psib=d["psib"], kchan=d["kchan"], run_id=d["run_id"], kind="run", group=d["group"])


def settings_for(job: dict) -> dict:
    """Every NAME = value line of the Settings cell that this run sets. Channel loss is on only for a job with a Kchan."""
    s = dict(BASE)
    s.update(KS_MULT=job["ks"], F_MULT=job["f"], PSIB_MULT=job["psib"],
             RUN_CHOICE=WIN["label"], INPUT_FILE=f"{job['run_id']}.in", RESULTS_NAME=job["run_id"])
    if job.get("kchan") is not None:
        s.update(OPTPERCOLATION=1, CHANNELCONDUCTIVITY=float(job["kchan"]))
    return s


RECORD_KEYS = ("GW_FRACTION", "BEDROCK_ADD_M", "KS_MULT", "PSIB_MULT", "F_MULT", "M_MULT", "THETAS_MULT", "ANISOTROPY", "KINEMVELCOEF",
               "FLOWEXP", "CHANNELROUGHNESS", "CHANNELWIDTHCOEFF", "OPTSNOW", "OPTPERCOLATION", "CHANNELCONDUCTIVITY")


def settings_record(job: dict) -> dict:
    """What the manifest remembers about the settings, to notice a changed setting on resume."""
    s = settings_for(job)
    rec = {k: s[k] for k in RECORD_KEYS}
    rec.update(START_DATE=WIN["start"], RUN_HOURS=WIN["hours"])
    return rec


def same_record(old, cur) -> bool:
    if not old:
        return False
    filled = {k: BASE[k] for k in RECORD_KEYS if k in BASE and k not in old}
    return {**filled, **old} == cur


def expected_files(s):
    """Names the notebook gives the bedrock and water-table grids for these settings."""
    bed = s["BEDROCK_ADD_M"]
    gw_name = f"InitGW_{int(round(s['GW_FRACTION'] * 100))}pct_plus{bed:g}m_mm.asc"
    bed_name = "bedrock_depth_m.asc" if bed == 0 else f"bedrock_depth_plus{bed:g}m.asc"
    return gw_name, bed_name


def timeout_for(args) -> float:
    if args.timeout_min:
        return args.timeout_min
    return max(15.0, 3.5 * WIN["hours"] * MIN_PER_1000H / 1000.0)


def describe(job: dict) -> str:
    base = f"draw {job['idx']} (Ks {job['ks']:g}, f {job['f']:g}, PsiB {job['psib']:g})"
    if job["kind"] == "baseline":
        return f"reference: {base}, channel loss OFF (already run)"
    return f"{base}, channel loss ON, Kchan {job['kchan']:g} mm/hr"


# --------------------------------------------------------------------------------------------
# Patching the Settings cell
# --------------------------------------------------------------------------------------------
def find_settings_cell(nb):
    for i, cell in enumerate(nb.cells):
        if cell.cell_type != "code":
            continue
        src = cell.source
        if "# SETTINGS: everything specific" in src or ("GW_FRACTION = " in src and "INPUT_FILE = " in src):
            return i
    raise RunError("Could not find the Settings cell in the notebook.")


def patch_settings(src: str, values: dict) -> str:
    """Change only the listed NAME = value lines (keeps the trailing # comment)."""
    for name, value in values.items():
        pat = re.compile(rf"^({re.escape(name)}[ \t]*=[ \t]*)(.*?)([ \t]*(?:#.*)?)$", re.M)
        hits = pat.findall(src)
        if len(hits) != 1:
            raise RunError(f"Settings cell: expected exactly one line '{name} = ...' but found {len(hits)}.")
        src = pat.sub(lambda m, v=value: m.group(1) + repr(v) + m.group(3), src, count=1)
    return src


def patch_window(src: str, start: str, hours: int) -> str:
    """Replace the line 'START_DATE, RUN_HOURS = RUNS[RUN_CHOICE]' by the window's own start and length."""
    pat = re.compile(r"^(START_DATE[ \t]*,[ \t]*RUN_HOURS[ \t]*=[ \t]*)(.*?)([ \t]*(?:#.*)?)$", re.M)
    hits = pat.findall(src)
    if len(hits) != 1:
        raise RunError(f"Settings cell: expected exactly one line 'START_DATE, RUN_HOURS = ...' but found {len(hits)}.")
    return pat.sub(lambda m: m.group(1) + repr((start, hours)) + m.group(3), src, count=1)


def verify_patched(src: str, values: dict, window: tuple | None = None) -> None:
    """Parse the patched Settings cell and confirm each patched name is assigned once, to the right value."""
    tree = ast.parse(src)
    found, win = {}, []
    for node in tree.body:
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        tgt = node.targets[0]
        if isinstance(tgt, ast.Name) and tgt.id in values:
            found.setdefault(tgt.id, []).append(ast.literal_eval(node.value))
        elif isinstance(tgt, ast.Tuple) and [getattr(e, "id", None) for e in tgt.elts] == ["START_DATE", "RUN_HOURS"]:
            try:
                win.append(tuple(ast.literal_eval(node.value)))
            except ValueError:
                win.append(None)
    for name, want in values.items():
        got = found.get(name, [])
        if len(got) != 1 or got[0] != want or type(got[0]) is not type(want):
            raise RunError(f"Settings cell after patching: {name} = {got}, wanted exactly one value {want!r}.")
    if window is not None and win != [tuple(window)]:
        raise RunError(f"Settings cell after patching: START_DATE, RUN_HOURS = {win}, wanted exactly {tuple(window)!r}.")


def patched_source(src: str, job: dict) -> str:
    s = settings_for(job)
    out = patch_window(patch_settings(src, s), WIN["start"], WIN["hours"])
    verify_patched(out, s, (WIN["start"], WIN["hours"]))
    return out


# --------------------------------------------------------------------------------------------
# Running the notebook headlessly
# --------------------------------------------------------------------------------------------
def make_kernel_manager():
    """A kernel that uses THIS python (so it sees the same packages as this script)."""
    from jupyter_client.kernelspec import KernelSpec
    from jupyter_client.manager import KernelManager
    km = KernelManager(kernel_name="python3")
    km._kernel_spec = KernelSpec(argv=[sys.executable, "-m", "ipykernel_launcher", "-f", "{connection_file}"],
                                 display_name="lhs", language="python")
    return km


def execute_nb(nb, lab: Path, cell_timeout: int):
    """Execute all cells of an nbformat notebook in a fresh kernel started in lab/. Returns the notebook."""
    import nbclient
    km = make_kernel_manager()
    client = nbclient.NotebookClient(nb, km=km, timeout=cell_timeout, startup_timeout=180,
                                     allow_errors=False, resources={"metadata": {"path": str(lab)}})
    try:
        client.execute()
    finally:
        try:
            km.shutdown_kernel(now=True)
        except Exception:
            pass
        try:
            km.cleanup_resources()
        except Exception:
            pass
    return nb


def outputs_text(nb) -> str:
    """All text a notebook printed (images dropped)."""
    chunks = []
    for cell in nb.cells:
        if cell.cell_type != "code":
            continue
        for out in cell.get("outputs", []):
            t = out.get("output_type")
            if t == "stream":
                chunks.append(out.get("text", ""))
            elif t in ("execute_result", "display_data"):
                txt = out.get("data", {}).get("text/plain")
                if txt:
                    chunks.append(txt + "\n")
            elif t == "error":
                chunks.append("\n".join(out.get("traceback", [])) + "\n")
    return "".join(chunks)


def build_run(job: dict, args, lab: Path, log_dir: Path) -> tuple:
    """Patch + execute the notebook. Returns (output_text, warnings). Raises RunError on any problem."""
    import nbformat
    rid = job["run_id"]
    nb = nbformat.read(str(lab / args.notebook), as_version=4)
    for cell in nb.cells:                               # start from empty outputs
        if cell.cell_type == "code":
            cell.outputs = []
            cell.execution_count = None
    idx = find_settings_cell(nb)
    nb.cells[idx].source = patched_source(nb.cells[idx].source, job)
    try:
        execute_nb(nb, lab, args.cell_timeout)
    except Exception as exc:                            # a cell failed, or the kernel did not start
        failed = log_dir / f"{rid}_build_FAILED_executed.ipynb"
        try:
            nbformat.write(nb, str(failed))
        except Exception:
            pass
        partial = outputs_text(nb)
        (log_dir / f"{rid}_build_output.txt").write_text(partial)
        tail = re.sub(r"\x1b\[[0-9;]*m", "", "\n".join(partial.splitlines()[-12:]))
        raise RunError(f"The notebook stopped with an error: {type(exc).__name__}: {str(exc)[:600]}\n"
                       f"  (partial notebook saved as {failed.name}; last output lines follow)\n{tail}")
    text = outputs_text(nb)
    (log_dir / f"{rid}_build_output.txt").write_text(text)
    return text, check_build_text(text, job)


def check_build_text(text: str, job: dict) -> list:
    """Lines the notebook should have printed. Returns a list of warnings; raises RunError on a problem."""
    s = settings_for(job)
    gw_name, bed_name = expected_files(s)
    must = [
        f"Run choice: {s['RUN_CHOICE']} | start {WIN['start']} | {WIN['hours']} hours",
        "Every expected setting is in the input file.",
        "OK: every hourly file the run needs is there.",
        f"(KS_MULT = {s['KS_MULT']:g},", f"(PSIB_MULT = {s['PSIB_MULT']:g},", f"(F_MULT = {s['F_MULT']:g},",
        f"(M_MULT = {s['M_MULT']:g},", f"(THETAS_MULT = {s['THETAS_MULT']:g},", f"(ANISOTROPY = {s['ANISOTROPY']:g})",
        "equal to the thetaS grid in use: True",
        f"Bedrock file in use: {INIT_DATA}/{bed_name}",
        f"GWATERFILE: {MODEL_SOIL}/{gw_name}",
        f"BEDROCKFILE: {MODEL_SOIL}/{bed_name}",
    ]
    missing = [m for m in must if m not in text]
    if not re.search(r"All \d+ checks passed\.", text):
        missing.append("All <n> checks passed.")
    if ("OK: thetaS is above theta*_s everywhere in the basin." not in text
            and "NOTE: the unscaled SOLUS thetaS is already at or below theta*_s" not in text):
        missing.append("the thetaS against theta*_s check (OK: or NOTE: line)")
    if missing:
        raise RunError("The notebook output is missing line(s) it should have printed:\n  " + "\n  ".join(missing)
                       + "\n  (If these are the thetaS or anisotropy lines, the notebook in lab/ is the old one: replace it with the new one.)")
    if re.search(r"^PROBLEM:", text, re.M):
        raise RunError("The notebook printed a PROBLEM line (see the build output file).")
    warns = [ln.strip() for ln in text.splitlines() if re.match(r"^CHECK\s", ln)]
    m = re.search(r"cells held at the [\d.]+ cap: thetaS (\d+), porosity (\d+)", text)
    if m and (int(m.group(1)) or int(m.group(2))):
        warns.append(f"thetaS cells held at the {THETAS_CAP:g} cap: {m.group(1)}, porosity {m.group(2)}")
    return warns


# --------------------------------------------------------------------------------------------
# Checks on what the build produced
# --------------------------------------------------------------------------------------------
def read_in_file(path: Path) -> dict:
    lines = Path(path).read_text().splitlines()
    found = {}
    for i, line in enumerate(lines):
        m = re.match(r"^([A-Z0-9]+):\s", line + " ")
        if m and i + 1 < len(lines):
            found[m.group(1)] = lines[i + 1].strip()
    return found


def same_value(a, b) -> bool:
    try:
        return float(a) == float(b)
    except (TypeError, ValueError):
        return str(a).strip() == str(b).strip()


def check_in_file(new_in: Path, ref_in: Path, job: dict) -> None:
    """The new .in must equal the reference .in except for the run-specific lines, and the routing keywords must have exactly
    their baseline values. Ks, f and PsiB do not appear in the .in file; check_grids covers them."""
    s = settings_for(job)
    gw_name, bed_name = expected_files(s)
    new, ref = read_in_file(new_in), read_in_file(ref_in)
    expect = {"RUNTIME": WIN["hours"], "STARTDATE": WIN["start"], "OUTFILENAME": f"results/{job['run_id']}",
              "GWATERFILE": f"{MODEL_SOIL}/{gw_name}", "BEDROCKFILE": f"{MODEL_SOIL}/{bed_name}"}
    for kw in ROUTING_KWS:
        expect[kw] = s[kw]
    expect["OPTPERCOLATION"] = s["OPTPERCOLATION"]              # the two loss keywords: 0 and 5 for the reference, 1 and Kchan for a loss run
    expect["CHANNELCONDUCTIVITY"] = s["CHANNELCONDUCTIVITY"]
    problems = []
    for key in sorted(set(ref) | set(new)):
        b, n = ref.get(key), new.get(key)
        if key in expect:
            if n is None or not same_value(n, expect[key]):
                problems.append(f"{key}: written {n!r}, expected {expect[key]!r}")
        elif b is None or n is None or not same_value(b, n):
            problems.append(f"{key}: reference {b!r} but this run {n!r}  (this keyword may not change)")
    if problems:
        raise RunError("The .in file differs from the reference in more than it should:\n  " + "\n  ".join(problems))


def read_asc(path: Path):
    """Read an ESRI ASCII grid. Returns (header dict, 2-D array)."""
    header, nhead = {}, 0
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) == 2 and parts[0].lower() in ("ncols", "nrows", "xllcorner", "yllcorner", "xllcenter",
                                                        "yllcenter", "cellsize", "nodata_value"):
                header[parts[0].lower()] = float(parts[1])
                nhead += 1
            else:
                break
    data = np.loadtxt(path, skiprows=nhead, dtype=float)
    return header, data.reshape(int(header["nrows"]), int(header["ncols"]))


def read_gdf(path: Path) -> dict:
    with open(path) as fh:
        return {r["Variable"].strip(): r["BasePath"].strip() for r in csv.DictReader(fh)}


def check_one_grid(lab: Path, var: str, base_model: str, base_src: str, mult, t0: float) -> str:
    """One soil grid in data/model/soil against its SOLUS grid: a plain copy, or the SOLUS values times the multiplier."""
    model_file, src_file = lab / f"{base_model}.asc", lab / f"{base_src}.asc"
    name = os.path.basename(base_model)
    if not model_file.is_file():
        raise RunError(f"Soil grid {var}: {model_file} not found.")
    if not src_file.is_file():
        raise RunError(f"Soil grid {var}: the SOLUS grid {src_file} not found.")
    if mult is None or mult == 1:
        if "_x" in os.path.basename(base_model):
            raise RunError(f"scgrid.gdf row {var} is '{base_model}', which has a multiplier in its name but none is set.")
        if not filecmp.cmp(model_file, src_file, shallow=False):
            raise RunError(f"Soil grid {var}: {model_file.name} differs from the SOLUS grid, but no multiplier is set.")
        return f"{var} x1"
    if not name.endswith(f"_x{mult:g}"):
        raise RunError(f"scgrid.gdf row {var} is '{base_model}', which does not match multiplier {mult:g}.")
    if model_file.stat().st_mtime < t0 - 2:
        raise RunError(f"Soil grid {var}: {model_file.name} is older than this build (stale file).")
    hm, gm = read_asc(model_file)
    hs, gs = read_asc(src_file)
    if gm.shape != gs.shape:
        raise RunError(f"Soil grid {var}: the scaled grid and the SOLUS grid have different shapes.")
    nd_s, nd_m = hs.get("nodata_value", -9999), hm.get("nodata_value", -9999)
    ok = (gs != nd_s) & (gm != nd_m) & (gs != 0)
    if var == "FD":
        ok &= gs * mult > F_FLOOR * 1.02                  # cells at the f floor do not scale
    if var in ("TS", "PO"):
        ok &= gs * mult < THETAS_CAP * 0.999              # cells at the cap do not scale
    if int(ok.sum()) < 10:
        raise RunError(f"Soil grid {var}: fewer than 10 cells left to check the multiplier on.")
    med = float(np.median(gm[ok] / gs[ok]))
    tot = float(gm[ok].sum() / gs[ok].sum())
    if abs(med / mult - 1) > 0.01 or abs(tot / mult - 1) > 0.01:
        raise RunError(f"Soil grid {var}: the scaled grid is {med:.4g} x the SOLUS grid (median) and {tot:.4g} x (total); "
                       f"the multiplier is {mult:g}.")
    return f"{var} x{mult:g}"


def check_grids(lab: Path, job: dict, t_build_start: float) -> str:
    """Bedrock, water-table and every soil grid / the soil table in data/model/soil must be fresh and correct."""
    s = settings_for(job)
    gw_name, bed_name = expected_files(s)
    soil = lab / MODEL_SOIL
    gw_path, bed_path = soil / gw_name, soil / bed_name
    for p in (gw_path, bed_path):
        if not p.is_file():
            raise RunError(f"Expected file not found: {p}")
    if gw_path.stat().st_mtime < t_build_start - 2:
        raise RunError(f"{gw_name} in data/model/soil is older than this build (stale file).")
    if not filecmp.cmp(bed_path, lab / INIT_DATA / bed_name, shallow=False):
        raise RunError(f"{bed_name} in data/model/soil differs from the copy in {INIT_DATA}.")
    hb, bed = read_asc(bed_path)
    hg, gw = read_asc(gw_path)
    nb_, ng_ = hb.get("nodata_value", -9999), hg.get("nodata_value", -9999)
    if bed.shape != gw.shape:
        raise RunError("Bedrock and water-table grids have different shapes.")
    ok = (bed != nb_) & (gw != ng_) & (bed > 0)
    if not ok.any():
        raise RunError("No valid cells to compare in the bedrock and water-table grids.")
    orig_h, orig = read_asc(lab / INIT_DATA / "bedrock_depth_m.asc")
    ov = (orig != orig_h.get("nodata_value", -9999)) & (bed != nb_)
    diff = (bed - orig)[ov]
    if diff.size == 0 or np.max(np.abs(diff - s["BEDROCK_ADD_M"])) > 0.011:
        raise RunError(f"Bedrock grid is not the SOLUS grid + {s['BEDROCK_ADD_M']:g} m "
                       f"(difference ranges {diff.min():.3f} to {diff.max():.3f} m).")
    med = float(np.median(gw[ok] / bed[ok]))
    if abs(med - s["GW_FRACTION"] * 1000) <= 5:
        unit = "mm"
    elif abs(med - s["GW_FRACTION"]) <= 0.005:
        unit = "m"
    else:
        raise RunError(f"Water-table grid / bedrock grid has median {med:.4f}; expected "
                       f"{s['GW_FRACTION'] * 1000:g} (mm per m) for GW_FRACTION = {s['GW_FRACTION']}.")
    model_rows, src_rows = read_gdf(soil / "scgrid.gdf"), read_gdf(lab / INIT_DATA / "scgrid.gdf")
    for var in ("KS", "TR", "TS", "PB", "PI", "FD", "PO"):
        if var not in model_rows or var not in src_rows:
            raise RunError(f"scgrid.gdf has no row {var} (model folder: {sorted(model_rows)}, SOLUS: {sorted(src_rows)}).")
    if model_rows["TS"] != model_rows["PO"]:
        raise RunError(f"thetaS and porosity point at different grids ({model_rows['TS']} and {model_rows['PO']}); "
                       f"they must be the same file so that the cap never cancels the thetaS multiplier.")
    if model_rows["PI"] in (model_rows["TS"], model_rows["PO"]):
        raise RunError(f"the pore-size index m (row PI: {model_rows['PI']}) shares a file with thetaS or porosity; it is expected to have its own grid.")
    mults = {"KS": s["KS_MULT"], "PB": s["PSIB_MULT"], "FD": s["F_MULT"], "PI": s["M_MULT"], "TS": s["THETAS_MULT"],
             "PO": s["THETAS_MULT"], "TR": None}
    notes = [check_one_grid(lab, var, model_rows[var], src_rows[var], mults[var], t_build_start) for var in mults]
    sdt = soil / "soils.sdt"
    if not sdt.is_file() or sdt.stat().st_mtime < t_build_start - 2:
        raise RunError("data/model/soil/soils.sdt is missing or older than this build.")
    table = [ln.split(",") for ln in sdt.read_text().splitlines()[1:] if ln.strip()]
    if not table or any(len(r) != 12 for r in table):
        raise RunError("soils.sdt does not have 12 comma-separated columns on every row.")
    for r in table:
        if float(r[7]) != float(s["ANISOTROPY"]) or float(r[8]) != float(s["ANISOTROPY"]):
            raise RunError(f"soils.sdt row {r[0]} has As = {r[7]}, Au = {r[8]}; ANISOTROPY is {s['ANISOTROPY']:g}.")
    return (f"Ks x{s['KS_MULT']:g}, PsiB x{s['PSIB_MULT']:g}, f x{s['F_MULT']:g}, m x{s['M_MULT']:g}, thetaS x{s['THETAS_MULT']:g}, "
            f"As = Au = {s['ANISOTROPY']:g}; bedrock +{s['BEDROCK_ADD_M']:g} m; water table = {s['GW_FRACTION']:g} x bedrock ({unit})")


def dir_mb(path: Path) -> float:
    total = 0
    for p in Path(path).glob("*"):
        try:
            if p.is_file():
                total += p.stat().st_size
        except OSError:
            pass
    return total / 1e6


def prune_scaled_grids(lab: Path) -> int:
    """Opt-in (--prune-soil). Delete scaled soil grids (names ending _x<number>.asc) that the CURRENT scgrid.gdf does not use.
    They can always be rebuilt by the notebook. Returns how many files were removed."""
    soil = lab / MODEL_SOIL
    try:
        used = {os.path.basename(v) + ".asc" for v in read_gdf(soil / "scgrid.gdf").values()}
    except Exception:
        return 0
    n = 0
    for p in soil.glob("*_x*.asc"):
        if re.search(r"_x[0-9.eE+-]+\.asc$", p.name) and p.name not in used:
            try:
                p.unlink()
                n += 1
            except OSError:
                pass
    return n


# --------------------------------------------------------------------------------------------
# Running tRIBS and scoring
# --------------------------------------------------------------------------------------------
def competing_tribs() -> list:
    out = subprocess.run(["ps", "-eo", "pid=,args="], capture_output=True, text=True).stdout
    hits = []
    for line in out.splitlines():
        line = line.strip()
        if not line:
            continue
        pid, _, rest = line.partition(" ")
        words = rest.split()
        if pid.isdigit() and int(pid) != os.getpid() and words and os.path.basename(words[0]) == "tRIBS":
            hits.append(line)
    return hits


def disk_free_gb(path=".") -> float:
    return shutil.disk_usage(path).free / 1e9


def count_rows(path: Path) -> int:
    try:
        with open(path, "rb") as fh:
            return sum(chunk.count(b"\n") for chunk in iter(lambda: fh.read(1 << 20), b""))
    except OSError:
        return 0


def kill_group(proc) -> None:
    try:
        os.killpg(proc.pid, signal.SIGTERM)
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait(timeout=10)
    except (ProcessLookupError, PermissionError):
        pass


def _on_term(signum, frame):
    raise KeyboardInterrupt


def install_term_handler() -> None:
    """Make 'kill' (SIGTERM) behave like Ctrl-C: stop tRIBS, mark the run interrupted, remove the lock. Running a notebook
    resets this handler (nbclient does it), so it is installed again before every tRIBS run."""
    signal.signal(signal.SIGTERM, _on_term)


def run_tribs(tribs: str, infile: str, logname: str, timeout_min: float,
              progress_file: Path | None = None, hours: int | None = None,
              heartbeat_min: float = HEARTBEAT_MIN, label: str = "") -> tuple:
    """Run tRIBS to completion. Returns (return code, seconds). Kills the whole process group on timeout.
    tRIBS prints nothing while it runs. Every `heartbeat_min` minutes this prints a progress line, which also counts as
    terminal activity, so a Codespace (30 minute idle timeout) does not shut down in the middle of a run."""
    t0 = time.time()
    deadline = t0 + timeout_min * 60
    beat = max(heartbeat_min, 0.01) * 60
    expected = hours * ROWS_PER_MODEL_HOUR if hours else None
    with open(logname, "w") as fh:
        proc = subprocess.Popen([tribs, infile], stdout=fh, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            while True:
                try:
                    rc = proc.wait(timeout=max(0.1, min(beat, deadline - time.time())))
                    break
                except subprocess.TimeoutExpired:
                    if time.time() >= deadline:
                        kill_group(proc)
                        raise RunError(f"tRIBS did not finish within {timeout_min:g} minutes and was stopped.")
                    msg = f"{label}: tRIBS still running, {(time.time() - t0) / 60:.1f} min so far"
                    if progress_file is not None:
                        n = count_rows(progress_file)
                        msg += f"; outlet file has {n} rows" + (f" of about {expected} ({100 * n / expected:.0f}%)" if expected else "")
                    log(msg)
        except BaseException:
            kill_group(proc)
            raise
    return rc, time.time() - t0


def check_outlet(path: Path, hours: int) -> dict:
    if not path.is_file():
        raise RunError(f"tRIBS did not write {path}.")
    rows = []
    with open(path) as fh:
        reader = csv.reader(fh)
        next(reader, None)
        for r in reader:
            if len(r) >= 2:
                try:
                    rows.append((float(r[0]), float(r[1])))
                except ValueError:
                    continue
    if len(rows) < 0.99 * hours:
        raise RunError(f"{path.name} has only {len(rows)} rows; expected about {hours}.")
    arr = np.array(rows)
    if not np.isfinite(arr[:, 1]).all():
        raise RunError(f"{path.name} contains NaN or infinite flows.")
    if arr[-1, 0] < hours - 2:
        raise RunError(f"{path.name} stops at hour {arr[-1, 0]:g}; the run should reach hour {hours}.")
    return {"rows": len(rows), "last_hour": float(arr[-1, 0]), "peak_m3s": float(arr[:, 1].max()),
            "md5": hashlib.md5(path.read_bytes()).hexdigest(), "dry_frac": float(np.mean(arr[:, 1] <= 1e-9))}


def run_scorer(rid: str, scorer: str, score_dir: Path, log_dir: Path, write_checksum: bool = True) -> dict:
    """Score one run on its window; returns its whole row of the scorer's CSV (all values as text)."""
    cmd = [sys.executable, scorer, "score", "--run", rid, "--window", WIN["name"], "--results-dir", "results",
           "--out-dir", str(score_dir), "--csv-name", SCORE_CSV]
    if write_checksum:
        cmd.append("--write-checksum")
    res = subprocess.run(cmd, capture_output=True, text=True)
    (log_dir / f"{rid}_score.txt").write_text((res.stdout or "") + (res.stderr or ""))
    if res.returncode != 0:
        raise RunError(f"The scorer exited with code {res.returncode}: {(res.stderr or res.stdout)[-400:]}")
    csv_path = score_dir / SCORE_CSV
    if csv_path.is_file():
        with open(csv_path, newline="") as fh:
            for row in csv.DictReader(fh):
                if row.get("run_id") == rid:
                    return dict(row)
    raise RunError(f"The scorer ran but {rid} is not in {csv_path}.")


def fnum(row, col) -> float:
    try:
        v = float((row or {}).get(col, ""))
        return v if math.isfinite(v) else float("nan")
    except (TypeError, ValueError):
        return float("nan")


def show_scores(row: dict) -> dict:
    """Keep every numeric score column the scorer wrote (rounded) and its arrival flags, in the manifest."""
    out = {}
    for c, v in (row or {}).items():
        if c in ("run_id", "window"):
            continue
        x = fnum(row, c)
        if math.isfinite(x):
            out[c] = round(x, 5)
        elif c.endswith("_flag") and v:
            out[c] = v
    return out


# --------------------------------------------------------------------------------------------
# Manifest, lock
# --------------------------------------------------------------------------------------------
def load_manifest() -> dict:
    if Path(MANIFEST).is_file():
        try:
            return json.loads(Path(MANIFEST).read_text())
        except json.JSONDecodeError:
            shutil.copy(MANIFEST, MANIFEST + ".corrupt")
    return {"runs": {}}


def save_manifest(man: dict) -> None:
    tmp = MANIFEST + ".tmp"
    Path(tmp).write_text(json.dumps(man, indent=2))
    os.replace(tmp, MANIFEST)


def take_lock() -> None:
    p = Path(LOCK)
    if p.is_file():
        try:
            pid = int(p.read_text().strip())
            os.kill(pid, 0)
            raise Abort(f"Another copy of this script is running (process {pid}). Lock file: {LOCK}")
        except (ValueError, ProcessLookupError, PermissionError):
            pass                                        # stale lock from a stopped run
    p.write_text(str(os.getpid()))
    atexit.register(lambda: p.unlink() if p.is_file() and p.read_text().strip() == str(os.getpid()) else None)


# --------------------------------------------------------------------------------------------
# The reference run (the W4 baseline: checked read-only, never rerun)
# --------------------------------------------------------------------------------------------
def verify_baseline_in(path: Path) -> None:
    """The baseline's .in file must carry the baseline routing values, loss off, snow on, and the W4 window."""
    got = read_in_file(path)
    gw_name, _ = expected_files(dict(BASE))
    want = {"STARTDATE": WIN["start"], "RUNTIME": WIN["hours"], "OPTPERCOLATION": BASE["OPTPERCOLATION"],
            "OPTSNOW": BASE["OPTSNOW"], "CHANNELCONDUCTIVITY": BASE["CHANNELCONDUCTIVITY"],
            "KINEMVELCOEF": BASE["KINEMVELCOEF"], "FLOWEXP": BASE["FLOWEXP"],
            "CHANNELROUGHNESS": BASE["CHANNELROUGHNESS"], "CHANNELWIDTHCOEFF": BASE["CHANNELWIDTHCOEFF"]}
    bad = [f"{k} is {got.get(k)!r}, wanted {v!r}" for k, v in want.items() if k not in got or not same_value(got[k], v)]
    gwf = got.get("GWATERFILE", "")
    if not gwf.endswith(gw_name):
        bad.append(f"GWATERFILE is {gwf!r}, wanted a file ending in {gw_name}")
    if bad:
        raise RunError(f"{path.name} is not a baseline-routing {WIN['name']} run: " + "; ".join(bad))


def reuse_check(lab: Path, job: dict) -> tuple:
    """Read-only check of the loss-off reference run (the same draw on W4). Returns (problems, notes). No problems = it can be used."""
    rid = job["run_id"]
    problems, notes = [], []
    q = lab / "results" / f"{rid}_Outlet.qout"
    if not q.is_file():
        return [f"results/{q.name} not found"], notes
    try:
        outlet = check_outlet(q, WIN["hours"])
    except RunError as exc:
        return [str(exc)], notes
    side = q.with_suffix(q.suffix + ".sha256")
    if side.is_file():
        saved = side.read_text().split()[0]
        if saved != hashlib.sha256(q.read_bytes()).hexdigest():
            problems.append(f"{q.name} no longer matches its .sha256 checksum file (the kept file has changed)")
        else:
            notes.append("checksum file matches")
    else:
        notes.append("no .sha256 file yet (it is written when the run is scored)")
    wm = Path(W4_MANIFEST)
    if not wm.is_file():
        problems.append(f"{W4_MANIFEST} not found, so there is no record that this run finished and was checked")
    else:
        try:
            rec = json.loads(wm.read_text())["runs"].get(rid)
        except Exception:
            rec = None
            problems.append(f"{W4_MANIFEST} could not be read")
        else:
            if rec is None:
                problems.append(f"{W4_MANIFEST} has no record of this run")
            else:
                if rec.get("status") != "done":
                    problems.append(f"{W4_MANIFEST} says this run's status is {rec.get('status')!r}, not 'done'")
                md5 = (rec.get("outlet") or {}).get("md5")
                if not md5:
                    problems.append(f"{W4_MANIFEST} has no outlet md5 for this run")
                elif md5 != outlet["md5"]:
                    problems.append(f"the outlet file's md5 ({outlet['md5']}) differs from the one recorded in {W4_MANIFEST} ({md5})")
                else:
                    notes.append(f"md5 matches {W4_MANIFEST}")
                if not same_record(rec.get("settings"), settings_record(job)):
                    problems.append(f"{W4_MANIFEST} records settings {rec.get('settings')}, expected {settings_record(job)}")
    inp = lab / f"{rid}.in"
    if not inp.is_file():
        problems.append(f"{inp.name} not found in lab/")
    else:
        try:
            verify_baseline_in(inp)
            notes.append(f"{inp.name} shows channel loss off and the baseline routing values")
        except RunError as exc:
            problems.append(str(exc))
    return problems, notes


def do_reuse(job: dict, args, lab: Path, man: dict) -> bool:
    """Check the loss-off reference read-only and score it into this batch's own score CSV. It is never rerun. Returns True if usable."""
    rid = job["run_id"]
    log_dir, score_dir = lab / LOG_DIR, lab / SCORE_DIR
    log_dir.mkdir(exist_ok=True)
    score_dir.mkdir(exist_ok=True)
    rec = man["runs"].setdefault(rid, {})
    cur = settings_record(job)
    if rec.get("status") == "reused" and same_record(rec.get("settings"), cur) and rec.get("scores"):
        log(f"{rid}: reference already checked and scored, skipping")
        return True
    rec.clear()
    rec.update(kind="baseline", status="checking", started=now_iso())
    problems, notes = reuse_check(lab, job)
    if not problems:
        try:
            q = lab / "results" / f"{rid}_Outlet.qout"
            has_sidecar = q.with_suffix(q.suffix + ".sha256").is_file()
            row = run_scorer(rid, str(lab / args.scorer), score_dir, log_dir, write_checksum=not has_sidecar)
            rec["scores"] = show_scores(row)
            rec["outlet"] = check_outlet(q, WIN["hours"])
        except RunError as exc:
            problems.append(str(exc))
    if problems:
        rec.update(status="reuse_rejected", problems=problems, finished=now_iso())
        save_manifest(man)
        log(f"{rid}: can NOT be used as the reference, and will NOT be rerun:")
        for p in problems:
            log(f"    {p}")
        return False
    rec.update(status="reused", finished=now_iso(), settings=cur)
    save_manifest(man)
    for n in notes:
        log(f"{rid}: {n}")
    log(f"{rid}: reference reused, not rerun. {one_line(rec['scores'])}")
    return True


# --------------------------------------------------------------------------------------------
# Preflight
# --------------------------------------------------------------------------------------------
def pick_notebook(lab: Path, given: str | None) -> tuple:
    if given:
        return given, ""
    present = [n for n in NOTEBOOK_CANDIDATES if (lab / n).is_file()]
    if not present:
        return NOTEBOOK_CANDIDATES[0], ""
    note = ""
    if len(present) > 1:
        note = (f"both {present[0]} and {present[1]} are in lab/; using {present[0]} (the one the earlier runners use). "
                f"They should have the same Settings cell. Use --notebook to choose the other.")
    return present[0], note


def preflight(args, lab: Path, jobs: list, base_problems: list, nb_note: str, man: dict) -> list:
    """Checks before anything is built. Returns a list of (level, message)."""
    out = []

    def add(level, msg):
        out.append((level, msg))

    for f in (args.notebook, args.scorer):
        add("OK" if (lab / f).is_file() else "FAIL", f"file in lab/: {f}")
    if nb_note:
        add("WARN", nb_note)
    if base_problems:
        add("FAIL", f"the loss-off reference run {BASELINE_ID} cannot be used: " + "; ".join(base_problems)
            + " (nothing was made, and nothing is rerun automatically)")
    else:
        add("OK", f"reference run {BASELINE_ID} passes its checks (reused, not rerun)")
    obs = lab / INIT_DATA / "obs" / "sycamore_09510200_obs_flow_WY2017.csv"
    add("OK" if obs.is_file() else "FAIL", f"gauge file: {INIT_DATA}/obs/{obs.name}")
    for f in ("../init_data/bedrock_depth_m.asc", "../init_data/scgrid.gdf", "../init_data/soils.sdt",
              "../init_data/met/Master_Met.sdf", "../init_data/LandUse.asc"):
        add("OK" if (lab / f).is_file() else "FAIL", f"input: {f}")
    tribs = args.tribs
    add("OK" if (os.path.isfile(tribs) and os.access(tribs, os.X_OK)) else "FAIL", f"tRIBS program: {tribs}")
    free = disk_free_gb(lab)
    add("OK" if free >= args.min_free_gb else "FAIL", f"free disk space: {free:.1f} GB (need {args.min_free_gb:g})")
    comp = competing_tribs()
    add("FAIL" if comp else "OK", "no other tRIBS running" + (": " + comp[0] if comp else ""))
    try:
        import nbformat
        nb = nbformat.read(str(lab / args.notebook), as_version=4)
        idx = find_settings_cell(nb)
        src = nb.cells[idx].source
        old = []
        for name in ("M_MULT", "THETAS_MULT", "ANISOTROPY"):
            if not re.search(rf"^{name}[ \t]*=", src, re.M):
                old.append(f"no line '{name} = ...' in the Settings cell")
        if not re.search(r"^SOIL_CONSTANTS[ \t]*=.*'As':\s*ANISOTROPY.*'Au':\s*ANISOTROPY", src, re.M):
            old.append("SOIL_CONSTANTS does not use ANISOTROPY")
        all_code = "\n".join(c.source for c in nb.cells if c.cell_type == "code")
        if "THETAS_CAP" not in all_code or "theta*_s" not in all_code:
            old.append("the thetaS cap or the theta*_s check is missing")
        if old:
            add("FAIL", f"{args.notebook} is the OLD notebook ({'; '.join(old)}). Replace it with the new Make_Sycamore_Model.ipynb, then run again.")
        else:
            for j in jobs:
                patched_source(src, j)
            add("OK", f"notebook Settings cell (cell {idx}) is the new one and patches correctly for all {len(jobs)} planned runs")
        cur = {n: re.search(rf"^{n}[ \t]*=[ \t]*([^#\n]*)", src, re.M) for n in ("F_MULT", "KS_MULT", "PSIB_MULT")}
        add("OK", "notebook Settings cell now reads " + ", ".join(f"{n} = {m.group(1).strip()}" for n, m in cur.items() if m)
            + " (the script overrides all three for every run, so these only matter when you run the notebook by hand)")
    except ImportError:
        add("FAIL", "nbformat is not installed: pip install nbformat nbclient ipykernel")
    except Exception as exc:
        add("FAIL", f"Settings cell check: {exc}")
    ref = lab / f"{BASELINE_ID}.in"
    try:
        got = read_in_file(ref)
        need = ROUTING_KWS + ["OPTPERCOLATION", "OPTSNOW", "CHANNELCONDUCTIVITY", "RUNTIME", "STARTDATE", "GWATERFILE",
                              "BEDROCKFILE", "OUTFILENAME"]
        miss = [k for k in need if k not in got]
        add("FAIL" if miss else "OK", f"reference .in ({ref.name}) has the keywords the comparison needs" + (f" (missing: {miss})" if miss else ""))
    except Exception as exc:
        add("FAIL", f"reference .in: {exc}")
    try:
        res = subprocess.run([sys.executable, str(lab / args.scorer), "score", "--help"], capture_output=True, text=True)
        flags = ["--run", "--window", "--results-dir", "--out-dir", "--csv-name", "--write-checksum"]
        miss = [f for f in flags if f not in res.stdout]
        add("FAIL" if (res.returncode != 0 or miss) else "OK",
            "scorer accepts " + " ".join(flags) + (f" (missing: {miss})" if miss else ""))
    except Exception as exc:
        add("FAIL", f"scorer check: {exc}")
    if not args.skip_kernel_check:
        try:
            from nbformat.v4 import new_code_cell, new_notebook
            nb = new_notebook(cells=[new_code_cell(
                "import sys, os\nfrom pytRIBS.classes import *\nprint('KERNEL', sys.executable, os.getcwd())")])
            execute_nb(nb, lab, 300)
            txt = outputs_text(nb).strip()
            add("OK", f"notebook kernel imports pytRIBS ({txt.replace('KERNEL ', '')})")
        except Exception as exc:
            add("FAIL", f"notebook kernel could not import pytRIBS: {type(exc).__name__}: {str(exc)[:300]}")
    results_dir = lab / "results"
    for j in jobs:
        stale = list(results_dir.glob(j["run_id"] + "*")) if results_dir.is_dir() else []
        if stale and j["run_id"] not in man["runs"]:
            add("FAIL", f"results/ already has {len(stale)} file(s) named {j['run_id']}*; this script never overwrites existing results")
    base_bytes = sum(p.stat().st_size for p in results_dir.glob(BASELINE_ID + "*") if p.is_file()) if results_dir.is_dir() else 0
    n_new = sum(1 for j in jobs if (man["runs"].get(j["run_id"]) or {}).get("status") != "done")
    need = base_bytes * n_new / 1e9
    add("FAIL" if (base_bytes and need > free - args.min_free_gb) else "OK",
        f"disk: one {WIN['name']} run's results are about {base_bytes / 1e6:.0f} MB (the reference's); {n_new} runs need about {need:.1f} GB; {free:.1f} GB free")
    soil = lab / MODEL_SOIL
    if soil.is_dir():
        scaled = [p for p in soil.glob("*_x*.asc") if re.search(r"_x[0-9.eE+-]+\.asc$", p.name)]
        one = (sum(p.stat().st_size for p in scaled) / len(scaled) / 1e6) if scaled else 0.0
        add("OK", f"data/model/soil holds {len(scaled)} scaled grid file(s) now"
                  + (f" (about {one:.1f} MB each). Each run adds up to 3 (Ks, f, PsiB), so {n_new} runs could add up to "
                     f"{3 * n_new * one / 1e3:.1f} GB unless the notebook replaces them; --prune-soil removes the unused ones." if one else
                     ". The script logs its size after each build."))
    kernels = subprocess.run(["ps", "-eo", "args="], capture_output=True, text=True).stdout.count("ipykernel_launcher")
    if kernels:
        add("WARN", f"{kernels} notebook kernel(s) are running (normal if VS Code is open). Do not run cells in "
                    f"{args.notebook} or anything that rebuilds data/model/ while the batch runs.")
    return out


# --------------------------------------------------------------------------------------------
# One run
# --------------------------------------------------------------------------------------------
def remove_stale(lab: Path, rid: str) -> None:
    for p in list((lab / "results").glob(rid + "*")):
        try:
            p.unlink() if p.is_file() else shutil.rmtree(p)
        except OSError:
            pass
    for p in (lab / f"{rid}.in", lab / f"simulation_{rid}.log"):
        if p.is_file():
            p.unlink()


def output_mb(lab: Path, rid: str) -> float:
    return round(sum(p.stat().st_size for p in (lab / "results").glob(rid + "*") if p.is_file()) / 1e6, 1)


def do_run(job: dict, args, lab: Path, man: dict, ref_in: Path) -> None:
    rid, s = job["run_id"], settings_for(job)
    rec = man["runs"].setdefault(rid, {})
    log_dir, score_dir = lab / LOG_DIR, lab / SCORE_DIR
    log_dir.mkdir(exist_ok=True)
    score_dir.mkdir(exist_ok=True)
    cur = settings_record(job)

    if rec.get("status") in ("done", "unscored") and not same_record(rec.get("settings"), cur) and rid not in args.redo_ids:
        raise RunError(f"{rid} is finished, but its recorded settings differ from the ones now planned "
                       f"({rec.get('settings')} -> {cur}). It is NOT rerun automatically; use --redo {rid} if you really want that.")
    if rec.get("status") == "done" and rid not in args.redo_ids:
        log(f"{rid}: already done, skipping")
        return
    if rec.get("status") == "unscored" and rid not in args.redo_ids:
        log(f"{rid}: ran earlier but was not scored; scoring again (tRIBS is not rerun)")
        rec["scores"] = show_scores(run_scorer(rid, str(lab / args.scorer), score_dir, log_dir))
        rec["status"] = "done"
        rec["finished"] = now_iso()
        save_manifest(man)
        return

    free = disk_free_gb(lab)
    if free < args.min_free_gb:
        raise Abort(f"Only {free:.1f} GB free; need {args.min_free_gb:g} GB. Free some space and run again to resume.")
    comp = competing_tribs()
    if comp:
        raise Abort("Another tRIBS is running, so the batch stops (two at once is not allowed): " + comp[0])

    rec.clear()
    rec.update(status="building", started=now_iso(), kind=job["kind"], group=job.get("group", ""), idx=job["idx"], ks=job["ks"], f=job["f"], psib=job["psib"], kchan=job.get("kchan"),
               settings=cur)
    save_manifest(man)
    remove_stale(lab, rid)
    log(f"=== draw {job['idx']:03d}  {rid}: {describe(job)} on {WIN['name']} ({WIN['hours']} h from {WIN['start']}) ===")

    t_build = time.time()
    text, warns = build_run(job, args, lab, log_dir)
    check_in_file(lab / s["INPUT_FILE"], ref_in, job)
    grid_note = check_grids(lab, job, t_build)
    soil_mb = dir_mb(lab / MODEL_SOIL)
    rec.update(build_min=round((time.time() - t_build) / 60, 2), warnings=warns, grid_check=grid_note, soil_dir_mb=round(soil_mb, 1),
               status="built")
    save_manifest(man)
    log(f"{rid}: build OK in {fmt_min(time.time() - t_build)}; .in file matches the reference except the allowed lines; "
        f"data/model/soil is {soil_mb:.0f} MB")
    log(f"{rid}: {grid_note}")
    for wmsg in warns:
        log(f"{rid}: notebook said: {wmsg}")
    if args.build_only:
        log(f"{rid}: --build-only, not running tRIBS")
        return

    install_term_handler()
    log(f"{rid}: starting tRIBS ({WIN['hours']} model hours, expect about {WIN['hours'] * MIN_PER_1000H / 1000:.0f} to {WIN['hours'] * MIN_PER_1000H / 500:.0f} min)")
    rec["status"] = "running"
    save_manifest(man)
    rc, secs = run_tribs(args.tribs, s["INPUT_FILE"], f"simulation_{rid}.log", timeout_for(args),
                         progress_file=lab / "results" / f"{rid}_Outlet.qout", hours=WIN["hours"],
                         heartbeat_min=args.heartbeat_min, label=rid)
    rec["tribs_min"] = round(secs / 60, 2)
    rec["tribs_returncode"] = rc
    if rc != 0:
        raise RunError(f"tRIBS exited with code {rc}; see simulation_{rid}.log")
    outlet = check_outlet(lab / "results" / f"{rid}_Outlet.qout", WIN["hours"])
    rec["outlet"] = outlet
    rec["status"] = "unscored"
    rec["output_mb"] = output_mb(lab, rid)
    save_manifest(man)
    log(f"{rid}: tRIBS finished in {fmt_min(secs)}; outlet peak {outlet['peak_m3s']:.1f} m3/s over {outlet['rows']} rows; "
        f"files {rec['output_mb']:.0f} MB; {disk_free_gb(lab):.1f} GB free")
    bmd5 = ((man["runs"].get(BASELINE_ID) or {}).get("outlet") or {}).get("md5")
    if bmd5 and outlet["md5"] == bmd5:
        rec["identical_to_baseline"] = True
        log(f"{rid}: WARNING: the outlet file is byte-identical to the loss-off reference's. tRIBS ran, but these settings gave exactly "
            f"the loss-off flow. Either channel loss has no effect at this Kchan, or tRIBS did not read OPTPERCOLATION / CHANNELCONDUCTIVITY.")
    try:
        rec["scores"] = show_scores(run_scorer(rid, str(lab / args.scorer), score_dir, log_dir))
        rec["status"] = "done"
    except RunError as exc:
        rec["score_error"] = str(exc)
        log(f"{rid}: scoring failed (the run itself is fine): {exc}")
    rec["finished"] = now_iso()
    save_manifest(man)
    if rec.get("scores"):
        log(f"{rid}: {one_line(rec['scores'])}")
    if args.prune_soil:
        n = prune_scaled_grids(lab)
        if n:
            log(f"{rid}: --prune-soil removed {n} unused scaled soil grid file(s); data/model/soil is now {dir_mb(lab / MODEL_SOIL):.0f} MB")


# --------------------------------------------------------------------------------------------
# The tables and the summary CSV
# --------------------------------------------------------------------------------------------
def fm(x, nd=2, sign=False) -> str:
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "-"
    if abs(x) < 0.5 * 10 ** -nd:
        x = 0.0
    return f"{x:+.{nd}f}" if sign else f"{x:.{nd}f}"


def load_scores(lab: Path) -> dict:
    p = lab / SCORE_DIR / SCORE_CSV
    out = {}
    if p.is_file():
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh):
                out[r["run_id"]] = r
    return out


def load_w4_scores(lab: Path) -> dict:
    """The W4 batch's own score file (read only): the W4 baseline and the anchors, loss off. {} if it is missing."""
    p = lab / W4_SCORES_CSV
    out = {}
    if p.is_file():
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh):
                out[r["run_id"]] = r
    return out


def one_line(scores: dict) -> str:
    names = [(C_FFV, "Oct-Dec flow Mm3"), (C_FFDAYS, "days>0.13"), (C_FFMAX, "max m3/s"), (C_VOL, "flood vol"), (C_PK, "peak %"),
             (C_ARR, "arrival h"), (C_OM, "Oct-Mar vol")]
    return ", ".join(f"{label} {scores[c]}" for c, label in names if c in scores)


def pctchg(x: float, ref: float) -> float:
    """Percent change of x relative to ref (nan if either is missing or ref is 0)."""
    if not (math.isfinite(x) and math.isfinite(ref)) or ref == 0:
        return float("nan")
    return 100.0 * (x - ref) / ref


def cut_ratio(false_chg: float, flood_chg: float) -> float:
    """False-flow cut (percent) divided by flood-week-volume cut (percent). Only when both fall by at least 1%."""
    if math.isfinite(false_chg) and math.isfinite(flood_chg) and false_chg <= -1.0 and flood_chg <= -1.0:
        return false_chg / flood_chg
    return float("nan")


def dry_pct(rec: dict) -> float:
    v = ((rec or {}).get("outlet") or {}).get("dry_frac")
    return float("nan") if v is None else 100.0 * v


def print_summary(man: dict, lab: Path, draws: list) -> None:
    scores = load_scores(lab)
    ref = scores.get(BASELINE_ID)
    ref_rec = man["runs"].get(BASELINE_ID) or {}
    a = ANCHOR
    log()
    log(f"CHANNEL-LOSS CHECK on W4 (Oct 1 2016 - Apr 1 2017). Draw {a['idx']} (Ks {a['ks']:g}, f {a['f']:g}, PsiB {a['psib']:g}); everything else at the baseline; "
        f"warm start (water table 0.95 x bedrock).")
    log("  Only OPTPERCOLATION (0 -> 1) and Kchan (mm/hr) change. 'off' is the same draw with channel loss off (already run, only read here).")
    log("  Oct-Dec flow = model volume Oct 1 - Dec 31 in million m3. days = days the model flow is above 0.13 m3/s. max = largest Oct-Dec flow (m3/s).")
    log("  dry h% = share of the 4,368 hours with exactly zero flow at the outlet. flood vol, pk, arr = flood-week volume ratio, peak error % and arrival error h.")
    log("  Feb18, Feb27 = volume ratio of the Feb 18-19 and Feb 27-28 storms. Oct-Mar = whole-window volume ratio (1.00 = the gauge's volume).")
    if ref is not None:
        log(f"  Gauge, Oct-Dec: volume {fm(fnum(ref, C_FFG), 3)} M m3, largest flow {fm(fnum(ref, C_FFGMAX), 2)} m3/s, "
            f"{fm(fnum(ref, C_FFGDAYS), 0)} days above {FF_THRESHOLD:g}.")
    else:
        log(f"  (the reference {BASELINE_ID} has no score yet)")
    log()
    head = (f"  {'Kchan':>7} | {'OctDec Mm3':>10} {'days':>5} {'max':>6} {'dry h%':>6} | {'flood vol':>9} {'pk %':>6} {'arr h':>6} | "
            f"{'Feb18':>6} {'Feb27':>6} | {'Oct-Mar':>7}")
    log(head)

    def line(label, sc, rec, tail=""):
        return (f"  {label:>7} | {fm(fnum(sc, C_FFV), 2):>10} {fm(fnum(sc, C_FFDAYS), 0):>5} {fm(fnum(sc, C_FFMAX), 1):>6} {fm(dry_pct(rec), 1):>6} | "
                f"{fm(fnum(sc, C_VOL), 2):>9} {fm(fnum(sc, C_PK), 1):>6} {fm(fnum(sc, C_ARR), 2):>6} | "
                f"{fm(fnum(sc, C_F18), 2):>6} {fm(fnum(sc, C_F27), 2):>6} | {fm(fnum(sc, C_OM), 2):>7}{tail}")

    rows = []
    if ref is not None:
        log(line("off", ref, ref_rec, "   (same draw, loss off)"))
        row = dict(idx=a["idx"], run_id=BASELINE_ID, group="reference", KS_MULT=a["ks"], F_MULT=a["f"], PSIB_MULT=a["psib"], KCHAN_MMHR="off", status="reference")
        row.update({k: v for k, v in ref.items() if k not in ("run_id", "window")})
        row["dry_pct"] = dry_pct(ref_rec)
        rows.append(row)
    w4s = load_w4_scores(lab)
    wb = w4s.get(W4_BASELINE_ID)
    if wb is not None:
        log(line("W4 base", wb, {}, "   (context: the W4 baseline, Ks 0.2, f 1, PsiB 0.5, loss off)"))
    n_done, done = 0, []
    for d in draws:
        rec = man["runs"].get(d["run_id"]) or {}
        status = rec.get("status", "not run")
        sc = scores.get(d["run_id"]) if status in ("done", "unscored") else None
        row = dict(idx=d["idx"], run_id=d["run_id"], group=d["group"], KS_MULT=d["ks"], F_MULT=d["f"], PSIB_MULT=d["psib"], KCHAN_MMHR=d["kchan"], status=status)
        if sc is None:
            log(f"  {d['kchan']:>7g} | ({status})")
            rows.append(row)
            continue
        n_done += 1
        extra = "   IDENTICAL to the loss-off run" if rec.get("identical_to_baseline") else ""
        log(line(f"{d['kchan']:g}", sc, rec, extra))
        row.update({k: v for k, v in sc.items() if k not in ("run_id", "window")})
        row["dry_pct"] = dry_pct(rec)
        row["identical_to_reference"] = bool(rec.get("identical_to_baseline"))
        row["build_min"], row["tribs_min"], row["output_mb"] = rec.get("build_min", ""), rec.get("tribs_min", ""), rec.get("output_mb", "")
        rows.append(row)
        done.append((d["kchan"], sc, row))

    log()
    log(f"  {n_done} of {len(draws)} runs scored; {sum(1 for r in rows if r['status'] == 'failed')} failed; "
        f"{sum(1 for r in rows if r['status'] == 'not run')} not run yet.")

    if ref is not None and done:
        r_ff, r_fv, r_pk = fnum(ref, C_FFV), fnum(ref, C_VOL), fnum(ref, C_PK)
        log()
        log("  CHANGE FROM THE LOSS-OFF RUN (percent of the loss-off value; flood peak in points of peak error):")
        log(f"  {'Kchan':>7} | {'OctDec flow':>11} {'flood vol':>9} {'flood pk':>8} {'Feb18':>7} {'Feb27':>7} {'Oct-Mar':>8} | {'cut ratio':>9}")
        changes = []
        for k, sc, row in done:
            ff = pctchg(fnum(sc, C_FFV), r_ff)
            fv = pctchg(fnum(sc, C_VOL), r_fv)
            pk = fnum(sc, C_PK) - r_pk
            f18, f27, om = pctchg(fnum(sc, C_F18), fnum(ref, C_F18)), pctchg(fnum(sc, C_F27), fnum(ref, C_F27)), pctchg(fnum(sc, C_OM), fnum(ref, C_OM))
            cr = cut_ratio(ff, fv)
            log(f"  {k:>7g} | {fm(ff, 0, True):>10}% {fm(fv, 0, True):>8}% {fm(pk, 0, True):>7}  {fm(f18, 0, True):>6}% {fm(f27, 0, True):>6}% {fm(om, 0, True):>7}% | {fm(cr, 1):>9}")
            row.update(chg_octdec_pct=ff, chg_flood_vol_pct=fv, chg_flood_peak_pts=pk, chg_feb18_pct=f18, chg_feb27_pct=f27, chg_octmar_pct=om, cut_ratio=cr)
            changes.append((k, ff, fv, pk))
        log("  cut ratio = (percent cut in Oct-Dec false flow) / (percent cut in flood-week volume), shown when both fall by 1% or more.")
        log("  A ratio well above the soil dial's (below) means loss removes the false flow more selectively than the soil knobs do.")
        if wb is not None:
            b_ff, b_fv, b_pk = fnum(wb, C_FFV), fnum(wb, C_VOL), fnum(wb, C_PK)
            lines = []
            for rid, r in w4s.items():
                if rid == W4_BASELINE_ID or not rid.startswith("W4_"):
                    continue
                ff, fv = pctchg(fnum(r, C_FFV), b_ff), pctchg(fnum(r, C_VOL), b_fv)
                m = re.search(r"_d(\d+)_", rid)
                lines.append(f"    draw {int(m.group(1)) if m else '?'}: false flow {fm(ff, 0, True)}%, flood volume {fm(fv, 0, True)}%, "
                             f"flood peak {fm(fnum(r, C_PK) - b_pk, 0, True)} points, cut ratio {fm(cut_ratio(ff, fv), 1)}")
            if lines:
                log("  For scale, the soil dial alone (each W4 run with loss off, against the W4 baseline):")
                for ln in lines:
                    log(ln)
        if len(changes) >= 2:
            lo, hi = changes[0], changes[-1]
            dec = math.log10(hi[0] / lo[0]) if lo[0] > 0 else 0.0
            if dec >= 0.3:
                log(f"  Rough sensitivity, Kchan {lo[0]:g} to {hi[0]:g} ({dec:.1f} tenfold steps): per tenfold increase in Kchan the false flow changes by "
                    f"{fm((hi[1] - lo[1]) / dec, 0, True)} points of its loss-off value, the flood volume by {fm((hi[2] - lo[2]) / dec, 0, True)} points, "
                    f"and the flood peak error by {fm((hi[3] - lo[3]) / dec, 0, True)} points. (A straight-line guess from the two end runs: the fewer tenfold steps apart they are, the rougher it is.)")
        log("  Check: Kchan 0.5 should look almost the same as the loss-off run. If it does not, turning loss on changes something besides the loss itself.")
        log("  Nothing is selected or dropped here. Reading it is the next step.")
    rejected = [rid for rid, rec in man["runs"].items() if rec.get("status") == "reuse_rejected"]
    if rejected:
        log()
        log("  NOT used (failed their checks) and NOT rerun: " + ", ".join(rejected))
    if rows:
        cols = []
        for r in rows:
            for c in r:
                if c not in cols:
                    cols.append(c)
        with open(lab / SUMMARY_CSV, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            for r in rows:
                w.writerow({c: r.get(c, "") for c in cols})
        log()
        log(f"  Wrote {SUMMARY_CSV} (all score columns for the reference and every run, plus the changes from the loss-off run).")


# --------------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------------
def parse_args(argv=None):
    ap = argparse.ArgumentParser(description="Channel-loss (Kchan) check on W4 for one soil set (see the top of this file).")
    ap.add_argument("--dry-run", action="store_true", help="check everything, list the plan, build nothing")
    ap.add_argument("--summary-only", action="store_true", help="print the tables from the runs already scored, then stop")
    ap.add_argument("--build-only", action="store_true", help="build and check the inputs but do not run tRIBS")
    ap.add_argument("--limit", type=int, default=None, help="make only the first N runs that still need work (finished runs do not count)")
    ap.add_argument("--draw", type=int, default=None, help=f"draw number of the soil set to test (default: the draw in {DRAWS_CSV} if it exists, else {DRAW_DEFAULT}); it must be one of the runs in {W4_LIST_CSV}")
    ap.add_argument("--kchan", default=None, help=f"channel conductivity values in mm/hr, in run order (default: the values in {DRAWS_CSV} if it exists, else {KCHAN_DEFAULT}); to add values, list all earlier ones plus the new ones")
    ap.add_argument("--allow-high", action="store_true", help=f"allow Kchan values above {KCHAN_SOFT_MAX:g} mm/hr")
    ap.add_argument("--redo", default=None, help="comma-separated run ids to run again even if done (you almost never want this)")
    ap.add_argument("--timeout-min", type=float, default=None, help="stop a tRIBS run after this many minutes (default: 3.5 x the expected time, at least 15)")
    ap.add_argument("--cell-timeout", type=int, default=1800, help="seconds allowed for one notebook cell")
    ap.add_argument("--min-free-gb", type=float, default=2.0)
    ap.add_argument("--heartbeat-min", type=float, default=HEARTBEAT_MIN,
                    help="print a progress line every this many minutes while tRIBS runs (default 5; keeps a Codespace awake)")
    ap.add_argument("--prune-soil", action="store_true",
                    help="after each run, delete scaled soil grids in data/model/soil that the current scgrid.gdf does not use (opt-in; the notebook rebuilds them)")
    ap.add_argument("--notebook", default=None, help="notebook to run (default: Make_Sycamore_Model.ipynb, else Make_Sycamore_Model_clean.ipynb)")
    ap.add_argument("--scorer", default=SCORER)
    ap.add_argument("--tribs", default=os.environ.get("TRIBS_BIN") or (TRIBS_DEFAULT if os.path.isfile(TRIBS_DEFAULT)
                                                                       else shutil.which("tRIBS") or TRIBS_DEFAULT))
    ap.add_argument("--skip-kernel-check", action="store_true", help="skip the pytRIBS import test in --dry-run")
    return ap.parse_args(argv)


def main(argv=None) -> int:
    global _LOGFILE, ANCHOR, BASELINE_ID
    args = parse_args(argv)
    lab = Path(__file__).resolve().parent
    os.chdir(lab)
    (lab / LOG_DIR).mkdir(exist_ok=True)
    _LOGFILE = str(lab / LOG_DIR / "batch_log.txt")
    try:
        ANCHOR, draws, draw_note = load_plan(lab, args.draw, args.kchan, args.allow_high, write=not (args.dry_run or args.summary_only))
    except Abort as exc:
        print(str(exc))
        return 1
    BASELINE_ID = ANCHOR["ref_id"]
    man = load_manifest()
    meta = dict(window=WIN["name"], draw=ANCHOR["idx"], reference=BASELINE_ID)
    if man.get("meta") and man["meta"] != meta:
        print(f"{MANIFEST} was made with a different draw ({man['meta']}) than this script now uses ({meta}). "
              f"Nothing was started. Rename the old manifest and {DRAWS_CSV} first; do not mix two lists.")
        return 1

    if args.summary_only:
        print_summary(man, lab, draws)
        return 0

    args.notebook, nb_note = pick_notebook(lab, args.notebook)
    args.redo_ids = set(x.strip() for x in (args.redo or "").split(",") if x.strip())
    base_problems, base_notes = reuse_check(lab, make_job(None))
    jobs = [make_job(d) for d in draws]

    def done_already(j):
        rec = man["runs"].get(j["run_id"], {})
        return rec.get("status") == "done" and same_record(rec.get("settings"), settings_record(j)) and j["run_id"] not in args.redo_ids

    n_done = sum(1 for j in jobs if done_already(j))
    n_todo = len(jobs) - n_done
    log(f"Working folder: {lab}")
    log(f"Notebook: {args.notebook}")
    log(f"Design: {draw_note}")
    kchan_txt = ", ".join("%g" % j["kchan"] for j in jobs)
    log(f"Plan: draw {ANCHOR['idx']} with channel loss ON at Kchan {kchan_txt} mm/hr on {WIN['name']} ({WIN['hours']} h): "
        f"{n_todo} to make, {n_done} already finished and skipped.")
    log(f"Reference run (loss off): existing {BASELINE_ID} ({'; '.join(base_notes) if not base_problems else 'CANNOT be used: ' + '; '.join(base_problems)})")
    t_run = WIN["hours"] * MIN_PER_1000H / 1000.0
    log(f"  Time: each run is a notebook build (about 1 min) + tRIBS (about {t_run:.0f} min, up to {2 * t_run:.0f} if the machine is slower), "
        f"so {n_todo} runs take roughly {fmt_hm(n_todo * (t_run + 1.0) * 60)} to {fmt_hm(n_todo * (2 * t_run + 1.0) * 60)}.")
    for j in jobs:
        log(f"    run: {j['run_id']:52s} {describe(j)}")

    results = preflight(args, lab, jobs, base_problems, nb_note, man)
    log()
    for level, msg in results:
        log(f"  {level:5s} {msg}")
    fails = [m for lv, m in results if lv == "FAIL"]
    log()
    if args.dry_run:
        log("DRY RUN: nothing was built or run. " + ("All checks passed." if not fails else f"{len(fails)} check(s) FAILED."))
        return 1 if fails else 0
    if fails:
        log(f"{len(fails)} check(s) FAILED, so nothing was started. Fix them and run again.")
        return 1

    try:
        take_lock()
    except Abort as exc:
        log(str(exc))
        return 1

    man["script"] = "run_kchan_w4_sycamore.py"
    man["meta"] = meta
    man["baseline_id"] = BASELINE_ID
    man["last_started"] = now_iso()
    save_manifest(man)
    install_term_handler()

    t_batch, code = time.time(), 0
    args.draws_total = len(jobs)
    try:
        if not do_reuse(make_job(None), args, lab, man):
            log("The loss-off reference run cannot be used, so nothing was run. Nothing is rerun automatically: "
                "tell me what the messages above say and we decide together.")
            return 1
        ref_in = lab / f"{BASELINE_ID}.in"
        todo = [j for j in jobs if not done_already(j)]
        if args.limit is not None:
            todo = todo[: args.limit]
        if not todo:
            log("Nothing left to run: every Kchan value is already done.")
        consecutive_fail, finished_secs = 0, []
        for k, j in enumerate(todo, 1):
            rid = j["run_id"]
            t_run = time.time()
            try:
                do_run(j, args, lab, man, ref_in)
                consecutive_fail = 0
                if (man["runs"].get(rid) or {}).get("status") == "done":
                    finished_secs.append(time.time() - t_run)
            except RunError as exc:
                man["runs"].setdefault(rid, {}).update(status="failed", error=str(exc), finished=now_iso())
                save_manifest(man)
                log(f"{rid}: FAILED. {exc}")
                consecutive_fail += 1
                code = 1
                if consecutive_fail >= 2:
                    log("Two runs failed in a row, so the batch stops (the cause is probably the same).")
                    break
            if finished_secs:
                left = len(todo) - k
                avg = sum(finished_secs) / len(finished_secs)
                log(f"--- {k} of {len(todo)} in this session done; average {fmt_min(avg)} per run; about {fmt_hm(left * avg)} left ---")
    except Abort as exc:
        log(f"BATCH STOPPED: {exc}")
        code = 1
    except KeyboardInterrupt:
        log("BATCH INTERRUPTED. Run the same command again to resume; finished runs are skipped.")
        for rid, rec in man["runs"].items():
            if rec.get("status") in ("building", "running", "checking"):
                rec["status"] = "interrupted"
        save_manifest(man)
        code = 130
    finally:
        save_manifest(man)

    log()
    log(f"Batch time: {fmt_hm(time.time() - t_batch)}")
    counts = {}
    for rec in man["runs"].values():
        counts[rec.get("status", "?")] = counts.get(rec.get("status", "?"), 0) + 1
    log("Run status counts: " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))
    bad = [(rid, rec) for rid, rec in man["runs"].items() if rec.get("status") in ("failed", "interrupted", "unscored", "reuse_rejected")]
    for rid, rec in bad:
        log(f"  {rid}: {rec.get('status')}" + (f" - {str(rec.get('error', ''))[:160]}" if rec.get("error") else ""))
    print_summary(man, lab, draws)
    return code


if __name__ == "__main__":
    sys.exit(main())