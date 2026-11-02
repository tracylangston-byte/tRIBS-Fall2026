#!/usr/bin/env python3
"""
run_w1b_validate_sycamore.py  --  W1-B, Step A (validation): does a short run that starts Feb 21 give the same Feb 27-28 volume as W2?

THE IDEA. W2 (Jan 18 to Mar 12, 1,272 h) takes about 12 minutes per run. If a short run that starts on Feb 21 and stops on Mar 3
(240 h, about 2.5 minutes) gives the same Feb 27-28 storm volume, then all 150 Latin hypercube settings can be run on Feb 27-28
in about 6 to 7 hours instead of about 31. A run that starts in the middle of the winter starts "cold" (it forgets the January flood
and the Feb 18-19 storm), so this script FIRST checks that the short run agrees with W2.

THE TEST (the pass rule was written before any run):
  Run the settings we already have on W2 (19 Latin hypercube draws + the W2 baseline) on the short window, then compare each run's
  Feb 27-28 volume ratio (model / gauge) with its W2 value. PASS = every run within 3 percent of its W2 value.
  If it fails, nothing is decided here: we look at the pattern (a constant offset is different from scatter) and the next step is
  to try an earlier start (Feb 17, which includes the Feb 18-19 storm), not to loosen the rule.

THE WINDOW. Start Feb 21 2017 00:00, 240 h, end Mar 3 00:00 (the end of the scorer's Feb 27-28 event window). Water table at
0.82 x bedrock depth (the same cold-start value W1 and W2 use). Everything else is the Oct 21 baseline: cv 6, flowexp 0.3, n 0.04,
width 2.33, channel loss off, snow on, bedrock +1 m. Only Ks, f and PsiB change, exactly as in W2.

WHAT IT READS. w1b_candidates.csv (the 19 draws; the baseline is added automatically as the first run), lhs_volume_list.csv (the
Stage 3 design, to check the values), the W2 summaries w2_summary.csv and w2b_summary.csv (the W2 Feb 27-28 values to compare with),
and the existing W2 baseline's .in file (W2_ks0p2_cv6_psib0p5_nocl.in, read-only, as the reference for the .in comparison).
It does NOT change the scorer: the scorer is told the start time (--start), and the run names start with W1B_ so it does not
mistake them for W1.

WHAT IT MAKES (all in lab/): w1b_list.csv (the frozen list), w1b_manifest.json, w1b_logs/, score_out_w1b/scores_w1b.csv,
w1b_summary.csv, and for each run <run>.in, simulation_<run>.log and results/<run>*. It shares w2.lock with the W2 runners so two
batches can never run at once. Nothing existing is overwritten (a planned run whose name already has files in results/ stops the check).

Run it from workspaces/Sycamore/lab/, in this order:
    python -u run_w1b_validate_sycamore.py --dry-run               # checks only, builds nothing
    python -u run_w1b_validate_sycamore.py --limit 1 --build-only  # builds the W1-B baseline and checks its .in file against W2's
    python -u run_w1b_validate_sycamore.py --limit 1               # the baseline run (about 4 min), prints its Feb 27-28 volume next to W2's
    python -u run_w1b_validate_sycamore.py --prune-soil            # everything left (about 1 h 10 min); run it again to resume after any stop
    python -u run_w1b_validate_sycamore.py --summary-only          # print the comparison table from what is finished

Do NOT run the notebook (or any cell that rebuilds data/model/) while this runs: every build rewrites lab/data/model/ and all the
runs share it. Keep the terminal open: the progress lines also stop the Codespace from timing out for being idle.
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
WIN = dict(name="W1B", start="02/21/2017/00/00", hours=240, label="feb21_to_mar3", gw=0.82, scorer_start="2017-02-21 00:00")   # Feb 21 to Mar 3, 10 days; label = the RUN_CHOICE text printed; scorer_start = the start time given to the scorer

# Everything that is not drawn (the baseline, with the cold-start water table).
BASE = {
    "GW_FRACTION": WIN["gw"], "BEDROCK_ADD_M": 1.0, "KS_MULT": 0.2, "PSIB_MULT": 0.5, "F_MULT": 1,
    "THETAS_MULT": 1, "ANISOTROPY": 1, "M_MULT": 1,
    "KINEMVELCOEF": 6, "FLOWEXP": 0.3, "CHANNELROUGHNESS": 0.04, "CHANNELWIDTHCOEFF": 2.33,
    "OPTSNOW": 1, "OPTPERCOLATION": 0, "CHANNELCONDUCTIVITY": 5,
    "REBUILD_SOIL": False,
}
SIG_FIGS = 4
ROUTING_KWS = ["KINEMVELCOEF", "FLOWEXP", "CHANNELROUGHNESS", "CHANNELWIDTHCOEFF"]   # never changed here; must stay at baseline

REF_ID = "W2_ks0p2_cv6_psib0p5_nocl"               # the existing W2 baseline: its .in file is the READ-ONLY reference
BASELINE_ID = "W1B_ks0p2_cv6_psib0p5_f1_d000_nocl"     # the W1-B baseline run made by this script (draw 0)
W2_SUMMARIES = ["w2_summary.csv", "w2b_summary.csv"]
PASS_TOL = 0.03              # Step A rule: every run within 3 percent of its W2 Feb 27-28 volume ratio
W2_REF = {}                  # draw number -> the W2 values to compare with (filled in main)
NOTEBOOK_CANDIDATES = ["Make_Sycamore_Model.ipynb", "Make_Sycamore_Model_clean.ipynb"]
SCORER = "score_run_sycamore.py"
TRIBS_DEFAULT = "/opt/tRIBS/build/tRIBS"
INIT_DATA = "../init_data"           # SOLUS_DIR in the notebook
MODEL_SOIL = "data/model/soil"
LOG_DIR = "w1b_logs"
SCORE_DIR = "score_out_w1b"
SCORE_CSV = "scores_w1b.csv"
SUMMARY_CSV = "w1b_summary.csv"
DRAWS_CSV = "w1b_list.csv"                              # the frozen list of runs
MANIFEST = "w1b_manifest.json"
LOCK = "w2.lock"
CANDIDATES_DEFAULT = "w1b_candidates.csv"   # the 19 draws to validate (written by hand, see the top of this file)
W1_DESIGN_CSV = "lhs_volume_list.csv"                       # the Stage 3 design (the values every candidate must match)
CONTROL_EVERY = 4            # every 4th run in the order is a control (when there are controls), so a stopped batch has both kinds
MIN_PER_1000H = 10.5         # measured on W1 (1.74 min for 168 h); short runs carry a fixed start-up cost, so this is higher than for W2
HEARTBEAT_MIN = 5.0
ROWS_PER_MODEL_HOUR = 16     # the outlet file gets 16 rows per model hour
F_FLOOR = 1e-4               # the notebook keeps f at or above this
THETAS_CAP = 0.99            # the notebook keeps thetaS and porosity below this

# Proposed tolerances (from the design doc). Used only to LABEL rows in the tables; nothing is selected or removed by them.
VOL_LO, VOL_HI = 0.90, 1.10  # model / gauge volume
PEAK_TOL = 15.0              # percent
ARR_TOL = 1.0                # hours (model minus gauge)
TAIL_MIN = 0.7               # tail volume ratio, at least
FLOOD_BAND = (0.90, 1.12)       # NEW anchor rule: flood-week volume ratio
FEB27_BAND = (0.90, 1.10)       # NEW anchor rule: Feb 27-28 volume ratio (Feb 18-19 is not scored)
DEAD_VOL = 0.05              # a flood-week volume ratio below this is called "near-dead" in the tables

C_VOL, C_PK, C_ARR = "flood_week__vol_ratio", "flood_week__peak_err_pct", "flood_week__arrival_err_h"
C_KGE, C_SPAN = "flood_week__kge2012", "span__vol_ratio"
C_F18, C_F27 = "feb18_19__vol_ratio", "feb27_28__vol_ratio"
C_BAL, C_BALG = "balance__feb27_28_over_feb18_19__model", "balance__feb27_28_over_feb18_19__gauge"
C_T1, C_T2 = "tail_jan26_feb2__tail_vol_ratio", "tail_mar3_mar11__tail_vol_ratio"
C_J2 = "jan_second_storm__vol_ratio"
# (label, score column, kind). kind: band = volume ratio in VOL_LO..VOL_HI; abs_pk = within PEAK_TOL; abs_arr = within ARR_TOL; tail = at least TAIL_MIN
TARGETS = [
    ("flood vol", C_VOL, "band"), ("flood peak", C_PK, "abs_pk"), ("flood arr", C_ARR, "abs_arr"),
    ("Feb18 vol", C_F18, "band"), ("Feb27 vol", C_F27, "band"),
    ("tail 1", C_T1, "tail"), ("tail 2", C_T2, "tail"), ("window vol", C_SPAN, "band"),
]

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
# The list of runs (read from w1b_candidates.csv, frozen in w1b_list.csv)
# --------------------------------------------------------------------------------------------
def vtag(x) -> str:
    """0.18 -> '0p18', 1.5 -> '1p5', 15 -> '15', 0.2374 -> '0p2374'."""
    return f"{x:g}".replace(".", "p")


def run_id_for(ks, f, psib, idx) -> str:
    """W2_ks0p2374_cv6_psib0p8123_f1p234_d007_nocl: window, Ks, cv, PsiB, then the other changes (f, draw number), loss setting last."""
    return f"{WIN['name']}_ks{vtag(ks)}_cv6_psib{vtag(psib)}_f{vtag(f)}_d{idx:03d}_nocl"


def interleave(main: list, ctrl: list) -> list:
    """Main-group runs in their file order, with a control after every (CONTROL_EVERY - 1) of them; leftovers at the end."""
    out, mi, ci = [], 0, 0
    while mi < len(main) or ci < len(ctrl):
        for _ in range(CONTROL_EVERY - 1):
            if mi < len(main):
                out.append(main[mi])
                mi += 1
        if ci < len(ctrl):
            out.append(ctrl[ci])
            ci += 1
        elif mi >= len(main):
            break
    return out


def read_candidates(path: Path, w1_design: Path) -> list:
    """Read w1b_candidates.csv and check each run against the Stage 3 design. Returns the draws in run order. Raises Abort on a problem."""
    if not path.is_file():
        raise Abort(f"{path} not found. Put w1b_candidates.csv in lab/ (the same folder as this script), or give --candidates.")
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    need = ["group", "draw", "run_id", "KS_MULT", "F_MULT", "PSIB_MULT"]
    if not rows or any(c not in rows[0] for c in need):
        raise Abort(f"{path} needs the columns {need}; found {list(rows[0]) if rows else 'an empty file'}.")
    design = {}
    if w1_design.is_file():
        for r in csv.DictReader(open(w1_design, newline="")):
            design[int(r["idx"])] = (float(r["KS_MULT"]), float(r["F_MULT"]), float(r["PSIB_MULT"]), r["run_id"])
    main, ctrl, problems = [], [], []
    for r in rows:
        try:
            idx = int(float(r["draw"]))
            ks, f, psib = float(r["KS_MULT"]), float(r["F_MULT"]), float(r["PSIB_MULT"])
        except ValueError:
            problems.append(f"unreadable row {r}")
            continue
        group = r["group"].strip()
        if group not in ("main", "control"):
            problems.append(f"draw {idx}: group is {group!r}, expected 'main' or 'control'")
            continue
        w1id = r["run_id"].strip()
        d = dict(idx=idx, ks=ks, f=f, psib=psib, group=group, w1_run_id=w1id, run_id=run_id_for(ks, f, psib, idx))
        if d["run_id"].replace(WIN["name"] + "_", "W1_", 1) != w1id:
            problems.append(f"draw {idx}: its W1 run id {w1id} does not match the name made from its values ({d['run_id']})")
        if design:
            want = design.get(idx)
            if want is None:
                problems.append(f"draw {idx} is not in {w1_design.name}")
            elif (want[0], want[1], want[2]) != (ks, f, psib) or want[3] != w1id:
                problems.append(f"draw {idx}: Ks/f/PsiB or run id differ from {w1_design.name} (candidates {ks:g}, {f:g}, {psib:g}; design {want[0]:g}, {want[1]:g}, {want[2]:g})")
        (main if group == "main" else ctrl).append(d)
    draws = interleave(main, ctrl)
    ids = [d["run_id"] for d in draws]
    if len(set(ids)) != len(ids):
        problems.append("a run appears twice in the candidates file")
    if not draws:
        problems.append("the candidates file has no runs")
    if problems:
        raise Abort(f"{path} has problems; nothing was started:\n  " + "\n  ".join(problems))
    return draws


def write_draws(path: Path, draws: list) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["idx", "KS_MULT", "F_MULT", "PSIB_MULT", "run_id", "group"])
        for d in draws:
            w.writerow([d["idx"], d["ks"], d["f"], d["psib"], d["run_id"], d["group"]])


def read_draws(path: Path) -> list:
    with open(path, newline="") as fh:
        return [dict(idx=int(r["idx"]), ks=float(r["KS_MULT"]), f=float(r["F_MULT"]), psib=float(r["PSIB_MULT"]), run_id=r["run_id"], group=r["group"])
                for r in csv.DictReader(fh)]


def load_candidates(lab: Path, cand_arg: str) -> tuple:
    """Read the candidates, then read or make the frozen list. Returns (draws, note). Raises Abort if the frozen list no longer matches."""
    cpath = Path(cand_arg)
    cpath = cpath if cpath.is_absolute() else lab / cpath
    base = dict(idx=0, ks=BASE["KS_MULT"], f=float(BASE["F_MULT"]), psib=BASE["PSIB_MULT"], group="baseline", w1_run_id="", run_id=BASELINE_ID)
    fresh = [base] + read_candidates(cpath, lab / W1_DESIGN_CSV)
    path = lab / DRAWS_CSV
    if not path.is_file():
        write_draws(path, fresh)
        return fresh, f"read {len(fresh)} runs from {cpath.name}, checked them against {W1_DESIGN_CSV}, and froze the list in {DRAWS_CSV}"
    saved = read_draws(path)
    same = len(saved) == len(fresh) and all(a["run_id"] == b["run_id"] and a["ks"] == b["ks"] and a["f"] == b["f"]
                                           and a["psib"] == b["psib"] and a["group"] == b["group"] for a, b in zip(saved, fresh))
    if not same:
        raise Abort(f"{DRAWS_CSV} is already there but does not match {cpath.name} now. Nothing was started. If you changed the candidates on "
                    f"purpose, rename the old {DRAWS_CSV} and {MANIFEST} first; do not mix two lists.")
    return saved, f"read the frozen list from {DRAWS_CSV} ({len(saved)} runs; it matches {cpath.name})"


# --------------------------------------------------------------------------------------------
# Jobs and settings
# --------------------------------------------------------------------------------------------
def make_job(d: dict | None) -> dict:
    if d is None:
        return dict(idx=0, ks=BASE["KS_MULT"], f=BASE["F_MULT"], psib=BASE["PSIB_MULT"], run_id=BASELINE_ID, kind="baseline", group="baseline")
    return dict(idx=d["idx"], ks=d["ks"], f=d["f"], psib=d["psib"], run_id=d["run_id"], kind="run", group=d["group"])


def settings_for(job: dict) -> dict:
    """Every NAME = value line of the Settings cell that this run sets."""
    s = dict(BASE)
    s.update(KS_MULT=job["ks"], F_MULT=job["f"], PSIB_MULT=job["psib"],
             RUN_CHOICE=WIN["label"], INPUT_FILE=f"{job['run_id']}.in", RESULTS_NAME=job["run_id"])
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
    if job["kind"] == "baseline":
        return "baseline (no change)"
    return f"{job['group']}: Ks {job['ks']:g}, f {job['f']:g}, PsiB {job['psib']:g}"


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
    problems = []
    for key in sorted(set(ref) | set(new)):
        b, n = ref.get(key), new.get(key)
        if key in expect:
            if n is None or not same_value(n, expect[key]):
                problems.append(f"{key}: written {n!r}, expected {expect[key]!r}")
        elif b is None or n is None or not same_value(b, n):
            problems.append(f"{key}: reference {b!r} but this run {n!r}  (this keyword may not change)")
    if problems:
        raise RunError("The .in file differs from the reference in more than it should:\n  " + "\n  ".join(problems)
                       + "\n  (The reference is the W2 baseline's .in. A start-date or run-length line is expected to differ and is already allowed. "
                         "If a different line differs, it may depend on the dates: do not edit anything; copy this message to Claude.)")


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
            "md5": hashlib.md5(path.read_bytes()).hexdigest()}


def run_scorer(rid: str, scorer: str, score_dir: Path, log_dir: Path, write_checksum: bool = True) -> dict:
    """Score one run on its window; returns its whole row of the scorer's CSV (all values as text)."""
    cmd = [sys.executable, scorer, "score", "--run", rid, "--start", WIN["scorer_start"], "--results-dir", "results",
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
# The reference files (the W2 baseline's .in file and the W2 Feb 27-28 values: read only)
# --------------------------------------------------------------------------------------------
def verify_reference_in(path: Path) -> None:
    """The reference .in (the W2 baseline's) must exist and carry the baseline routing values, loss off, snow on, 0.82 water table.
    Its start date and run length are NOT checked: they belong to W2, and the new runs are compared with it on every other line."""
    if not path.is_file():
        raise RunError(f"{path.name} not found in lab/ (the existing W2 baseline's .in file is the reference for the .in comparison)")
    got = read_in_file(path)
    gw_name, _ = expected_files(dict(BASE))
    want = {"OPTPERCOLATION": BASE["OPTPERCOLATION"], "OPTSNOW": BASE["OPTSNOW"], "CHANNELCONDUCTIVITY": BASE["CHANNELCONDUCTIVITY"],
            "KINEMVELCOEF": BASE["KINEMVELCOEF"], "FLOWEXP": BASE["FLOWEXP"], "CHANNELROUGHNESS": BASE["CHANNELROUGHNESS"],
            "CHANNELWIDTHCOEFF": BASE["CHANNELWIDTHCOEFF"]}
    bad = [f"{k} is {got.get(k)!r}, wanted {v!r}" for k, v in want.items() if k not in got or not same_value(got[k], v)]
    if not got.get("GWATERFILE", "").endswith(gw_name):
        bad.append(f"GWATERFILE is {got.get('GWATERFILE')!r}, wanted a file ending in {gw_name}")
    if bad:
        raise RunError(f"{path.name} is not a baseline-routing run with the 0.82 water table: " + "; ".join(bad))


def load_w2_reference(lab: Path) -> dict:
    """draw number -> the W2 Feb 27-28 values (volume ratio, peak error %, arrival error h), from the W2 summaries; draw 0 = the W2 baseline."""
    ref = {}
    for fn in W2_SUMMARIES:
        p = lab / fn
        if p.is_file():
            with open(p, newline="") as fh:
                for r in csv.DictReader(fh):
                    try:
                        idx = int(float(r["idx"]))
                    except (KeyError, ValueError):
                        continue
                    v = fnum(r, C_F27)
                    if math.isfinite(v):
                        ref[idx] = dict(vol=v, pk=fnum(r, "feb27_28__peak_err_pct"), arr=fnum(r, "feb27_28__arrival_err_h"))
    for d, fn in (("score_out_w2b", "scores_w2b.csv"), ("score_out_w2", "scores_w2.csv")):
        p = lab / d / fn
        if p.is_file() and 0 not in ref:
            with open(p, newline="") as fh:
                for r in csv.DictReader(fh):
                    if r.get("run_id") == REF_ID and math.isfinite(fnum(r, C_F27)):
                        ref[0] = dict(vol=fnum(r, C_F27), pk=fnum(r, "feb27_28__peak_err_pct"), arr=fnum(r, "feb27_28__arrival_err_h"))
                        break
    return ref


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
        add("FAIL", f"the reference .in ({REF_ID}.in) cannot be used: " + "; ".join(base_problems)
            + " (nothing was made)")
    else:
        add("OK", f"reference .in ({REF_ID}.in, the W2 baseline's) has the baseline routing values; it is only read")
    obs = lab / INIT_DATA / "obs" / "sycamore_09510200_obs_flow_WY2017.csv"
    add("OK" if obs.is_file() else "FAIL", f"gauge file: {INIT_DATA}/obs/{obs.name}")
    for f in ("../init_data/bedrock_depth_m.asc", "../init_data/scgrid.gdf", "../init_data/soils.sdt",
              "../init_data/met/Master_Met.sdf", "../init_data/LandUse.asc"):
        add("OK" if (lab / f).is_file() else "FAIL", f"input: {f}")
    tribs = args.tribs
    add("OK" if (os.path.isfile(tribs) and os.access(tribs, os.X_OK)) else "FAIL", f"tRIBS program: {tribs}")
    free = disk_free_gb(lab)
    add("OK" if free >= args.min_free_gb else "FAIL", f"free disk space: {free:.1f} GB (need {args.min_free_gb:g})")
    lacking = [j["idx"] for j in jobs if j["idx"] not in W2_REF]
    add("FAIL" if lacking else "OK", f"W2 Feb 27-28 values to compare with: found for {len(jobs) - len(lacking)} of {len(jobs)} planned runs"
        + (f" (missing for draws {lacking}; need {' and '.join(W2_SUMMARIES)} and score_out_w2b/scores_w2b.csv in lab/)" if lacking else ""))
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
    ref = lab / f"{REF_ID}.in"
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
        flags = ["--run", "--start", "--results-dir", "--out-dir", "--csv-name", "--write-checksum"]
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
    base_bytes = (sum(p.stat().st_size for p in results_dir.glob(REF_ID + "*") if p.is_file()) * WIN["hours"] / 1272) if results_dir.is_dir() else 0
    n_new = sum(1 for j in jobs if (man["runs"].get(j["run_id"]) or {}).get("status") != "done")
    need = base_bytes * n_new / 1e9
    add("FAIL" if (base_bytes and need > free - args.min_free_gb) else "OK",
        f"disk: one {WIN['name']} run's results are about {base_bytes / 1e6:.0f} MB (the W2 baseline's, scaled to {WIN['hours']} h); {n_new} runs need about {need:.1f} GB; {free:.1f} GB free")
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
    rec.update(status="building", started=now_iso(), kind=job["kind"], group=job.get("group", ""), idx=job["idx"], ks=job["ks"], f=job["f"], psib=job["psib"],
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
    if bmd5 and outlet["md5"] == bmd5 and rid != BASELINE_ID:
        rec["identical_to_baseline"] = True
        log(f"{rid}: WARNING: the outlet file is byte-identical to the baseline's. tRIBS ran, but these settings gave exactly "
            f"the baseline flow. Either the three knobs have no pathway to the outlet here, or tRIBS does not read the setting.")
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
        w2 = W2_REF.get(job["idx"])
        v = fnum(rec["scores"], C_F27)
        if w2 and math.isfinite(v) and w2["vol"] > 0:
            r = v / w2["vol"]
            log(f"{rid}: Feb 27-28 volume {v:.4f} on W1-B vs {w2['vol']:.4f} on W2: ratio {r:.4f} "
                f"({'inside' if abs(r - 1) <= PASS_TOL else 'OUTSIDE'} the {PASS_TOL:.0%} rule)")
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


def one_line(scores: dict) -> str:
    names = [(C_F27, "Feb27 vol"), ("feb27_28__peak_err_pct", "Feb27 peak %"), ("feb27_28__arrival_err_h", "Feb27 arrival h"), (C_F27 + "_x", "")]
    return ", ".join(f"{label} {scores[c]}" for c, label in names if c in scores and label)


def spearman(a: list, b: list) -> float:
    if len(a) < 3:
        return float("nan")
    ra, rb = np.argsort(np.argsort(a)), np.argsort(np.argsort(b))
    return float(np.corrcoef(ra, rb)[0, 1])


def print_summary(man: dict, lab: Path, draws: list) -> None:
    scores, w2 = load_scores(lab), load_w2_reference(lab)
    log()
    log(f"STEP A: DOES THE SHORT RUN (Feb 21 - Mar 3, {WIN['hours']} h) GIVE THE SAME FEB 27-28 VOLUME AS W2?")
    log("  vol = Feb 27-28 volume ratio (model / gauge). W2 = from the long run; W1B = from the short run. ratio = W1B / W2.")
    log(f"  RULE (fixed before any run): every run's ratio within {PASS_TOL:.0%} of 1 (between {1 - PASS_TOL:.2f} and {1 + PASS_TOL:.2f}).")
    log()
    log(f"  {'draw':>4} {'Ks':>7} {'f':>6} {'PsiB':>6} | {'W2 vol':>7} {'W1B vol':>8} {'ratio':>6} {'rule':>7} | {'pk% W2':>7} {'pk% W1B':>8} | {'arr W2':>7} {'arr W1B':>8}")
    rows, pairs = [], []
    for d in draws:
        rec = man["runs"].get(d["run_id"]) or {}
        status = rec.get("status", "not run")
        sc = scores.get(d["run_id"]) if status in ("done", "unscored") else None
        ref = w2.get(d["idx"])
        row = dict(idx=d["idx"], run_id=d["run_id"], group=d["group"], KS_MULT=d["ks"], F_MULT=d["f"], PSIB_MULT=d["psib"], status=status)
        pre = f"  {d['idx']:>4} {d['ks']:>7g} {d['f']:>6g} {d['psib']:>6g} | "
        if sc is None:
            if status != "not run":
                log(pre + f"({status})")
            rows.append(row)
            continue
        v, pk, ar = fnum(sc, C_F27), fnum(sc, "feb27_28__peak_err_pct"), fnum(sc, "feb27_28__arrival_err_h")
        wv = ref["vol"] if ref else float("nan")
        ratio = v / wv if (math.isfinite(v) and math.isfinite(wv) and wv > 0) else float("nan")
        ok = math.isfinite(ratio) and abs(ratio - 1) <= PASS_TOL
        log(pre + f"{fm(wv, 4):>7} {fm(v, 4):>8} {fm(ratio, 3):>6} {('inside' if ok else 'OUTSIDE'):>7} | "
                  f"{fm(ref['pk'] if ref else float('nan'), 1):>7} {fm(pk, 1):>8} | {fm(ref['arr'] if ref else float('nan'), 2):>7} {fm(ar, 2):>8}")
        row.update(w2_feb27_vol=wv, w1b_feb27_vol=v, ratio_w1b_over_w2=ratio, inside_rule=ok, w2_feb27_peak_err_pct=(ref or {}).get("pk"),
                   w1b_feb27_peak_err_pct=pk, w2_feb27_arrival_err_h=(ref or {}).get("arr"), w1b_feb27_arrival_err_h=ar,
                   tribs_min=rec.get("tribs_min", ""), build_min=rec.get("build_min", ""))
        rows.append(row)
        if math.isfinite(ratio):
            pairs.append((d["idx"], wv, v, ratio, ok))
    n_done = sum(1 for r in rows if "ratio_w1b_over_w2" in r)
    waiting = [str(r["idx"]) for r in rows if r["status"] == "not run"]
    if waiting:
        log(f"  not run yet: draws {', '.join(waiting)}")
    log()
    log(f"  {n_done} of {len(draws)} runs scored; {sum(1 for r in rows if r['status'] == 'failed')} failed; {sum(1 for r in rows if r['status'] == 'not run')} not run yet.")
    if pairs:
        ratios = [p[3] for p in pairs]
        worst = max(pairs, key=lambda p: abs(p[3] - 1))
        n_ok = sum(1 for p in pairs if p[4])
        rho = spearman([p[1] for p in pairs], [p[2] for p in pairs])
        log(f"  Inside the {PASS_TOL:.0%} rule: {n_ok} of {len(pairs)}.  Mean ratio {np.mean(ratios):.3f} (1 = no systematic shift).  "
            f"Largest miss: draw {worst[0]} at {worst[3]:.3f}.  Rank agreement W1B vs W2 (Spearman): {('%.3f' % rho) if math.isfinite(rho) else 'n/a (needs 3 or more runs)'}.")
        if n_done < len(draws):
            log("  RESULT: not final (some runs are not scored yet).")
        elif n_ok == len(pairs):
            log("  RESULT: PASS. The short window reproduces W2's Feb 27-28 volume within the rule for every run.")
        else:
            spread = float(np.std(ratios))
            log(f"  RESULT: FAIL ({len(pairs) - n_ok} run(s) outside the rule). Spread of the ratios (sd) {spread:.3f}; mean {np.mean(ratios):.3f}. "
                "A mean far from 1 with a small spread means a roughly constant shift; a large spread means the short window changes how the settings act. "
                "Nothing is loosened or decided here: send this table and we choose the next step (the planned one is a start on Feb 17).")
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
        log(f"  Wrote {SUMMARY_CSV}.")


# --------------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------------
def parse_args(argv=None):
    ap = argparse.ArgumentParser(description="W1-B Step A: run the W2 settings on the short Feb 21 - Mar 3 window and compare (see the top of this file).")
    ap.add_argument("--dry-run", action="store_true", help="check everything, list the plan, build nothing")
    ap.add_argument("--summary-only", action="store_true", help="print the tables from the runs already scored, then stop")
    ap.add_argument("--build-only", action="store_true", help="build and check the inputs but do not run tRIBS")
    ap.add_argument("--limit", type=int, default=None, help="make only the first N runs that still need work (finished runs do not count)")
    ap.add_argument("--candidates", default=CANDIDATES_DEFAULT, help=f"the chosen runs (default {CANDIDATES_DEFAULT})")
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
    global _LOGFILE, W2_REF
    args = parse_args(argv)
    lab = Path(__file__).resolve().parent
    os.chdir(lab)
    (lab / LOG_DIR).mkdir(exist_ok=True)
    _LOGFILE = str(lab / LOG_DIR / "batch_log.txt")
    try:
        draws, draw_note = load_candidates(lab, args.candidates)
    except Abort as exc:
        print(str(exc))
        return 1
    man = load_manifest()
    meta = dict(window=WIN["name"], n=len(draws), run_ids_md5=hashlib.md5(",".join(d["run_id"] for d in draws).encode()).hexdigest())
    if man.get("meta") and man["meta"] != meta:
        print(f"{MANIFEST} was made with a different list of runs ({man['meta']}) than this script now uses ({meta}). "
              f"Nothing was started. Rename the old manifest and {DRAWS_CSV} first; do not mix two lists.")
        return 1

    if args.summary_only:
        print_summary(man, lab, draws)
        return 0

    args.notebook, nb_note = pick_notebook(lab, args.notebook)
    args.redo_ids = set(x.strip() for x in (args.redo or "").split(",") if x.strip())
    W2_REF = load_w2_reference(lab)
    base_problems, base_notes = [], ["the existing W2 baseline's .in file, read-only, is the reference for the .in comparison"]
    try:
        verify_reference_in(lab / f"{REF_ID}.in")
    except RunError as exc:
        base_problems.append(str(exc))
    jobs = [make_job(d) for d in draws]

    def done_already(j):
        rec = man["runs"].get(j["run_id"], {})
        return rec.get("status") == "done" and same_record(rec.get("settings"), settings_record(j)) and j["run_id"] not in args.redo_ids

    n_done = sum(1 for j in jobs if done_already(j))
    n_todo = len(jobs) - n_done
    log(f"Working folder: {lab}")
    log(f"Notebook: {args.notebook}")
    log(f"Design: {draw_note}")
    log(f"Plan: {len(jobs)} runs on {WIN['name']} (start {WIN['start']}, {WIN['hours']} h): the W1-B baseline + {len(jobs) - 1} Latin hypercube draws already run on W2; "
        f"{n_todo} to make, {n_done} already finished and skipped.")
    log(f"Reference .in: existing {REF_ID}.in ({'; '.join(base_notes) if not base_problems else 'CANNOT be used: ' + '; '.join(base_problems)})")
    for grp in ("baseline", "main"):
        g = [d for d in draws if d["group"] == grp]
        if g:
            log(f"  {grp:<8} group: {len(g)} runs; Ks {min(d['ks'] for d in g):g} to {max(d['ks'] for d in g):g}, f {min(d['f'] for d in g):g} to "
                f"{max(d['f'] for d in g):g}, PsiB {min(d['psib'] for d in g):g} to {max(d['psib'] for d in g):g}")
    t_run = WIN["hours"] * MIN_PER_1000H / 1000.0
    log(f"  Time: each run is a notebook build (about 1 min) + tRIBS (about {t_run:.0f} min, up to {2 * t_run:.0f} if the machine is slower), "
        f"so {n_todo} runs take roughly {fmt_hm(n_todo * (t_run + 1.0) * 60)} to {fmt_hm(n_todo * (2 * t_run + 1.0) * 60)}.")
    for j in jobs[: min(len(jobs), 5)]:
        log(f"    first runs: {j['run_id']:48s} {describe(j)}")

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

    man["script"] = "run_w1b_validate_sycamore.py"
    man["meta"] = meta
    man["baseline_id"] = BASELINE_ID
    man["last_started"] = now_iso()
    save_manifest(man)
    install_term_handler()

    t_batch, code = time.time(), 0
    args.draws_total = len(jobs)
    try:
        ref_in = lab / f"{REF_ID}.in"
        todo = [j for j in jobs if not done_already(j)]
        if args.limit is not None:
            todo = todo[: args.limit]
        if not todo:
            log("Nothing left to run: every draw is already done.")
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
    bad = [(rid, rec) for rid, rec in man["runs"].items() if rec.get("status") in ("failed", "interrupted", "unscored")]
    for rid, rec in bad:
        log(f"  {rid}: {rec.get('status')}" + (f" - {str(rec.get('error', ''))[:160]}" if rec.get("error") else ""))
    print_summary(man, lab, draws)
    return code


if __name__ == "__main__":
    sys.exit(main())
