#!/usr/bin/env python3
"""
run_tail_screen_sycamore.py  --  Does a lower f (and then a higher SATURATED anisotropy ratio) lift the slow flow after the big storm?

Question it answers: after the January flood, the gauge keeps flowing at about 0.3 to 1 m3/s for weeks, and the model drops
to near zero within about 2 days. Josh's suggestion is to decrease f and increase the saturated anisotropy ratio (As). This
script tests that on W4 (Oct 1 2016 - Apr 1 2017, 4,368 h, water table at 0.95 x bedrock at the start), starting from the
draw 12 / k10 setting:

    Ks x0.187, f x0.5874, PsiB x0.2939, m x1, thetaS x1, As = Au = 1, bedrock +1 m,
    cv 6.39, flowexp 0.2973, n 0.067, width coefficient 2.33, channel loss OFF (OPTPERCOLATION 0), snow ON (OPTSNOW 1).

W4 is used (not W2) because it includes the dry months Oct - Dec, when the gauge shows essentially no flow. A lower f and a
higher As could make the model's false flow in those months worse, and only a W4 run shows it. The water table has the whole
autumn to settle before the January flood, so the flood-week numbers here will not equal the W2 numbers.

THREE STEPS, one script (each step is skipped if its runs are already finished, so you can stop and resume any time):

    step 1  base   the draw 12 / k10 baseline on W4 (1 run). This is the "before" picture, and every later run is checked
                   against its .in file.
    step 2  f      f screen, one knob at a time: F_MULT = 0.40, 0.25, 0.15, 0.08, 0.04, 0.02 (6 runs; baseline is 0.5874).
                   Everything else stays at the baseline. The lowest value runs first, so a problem at the extreme shows up early.
    step 3  an     saturated-anisotropy screen at ONE lower f: As = 10, 100, 300, 900 with Au kept at 1 (4 runs; baseline is 1;
                   900 is the value Josh used in his Beaver Creek model), PLUS one diagnostic run with As = Au = 900 at the same f.
                   If the diagnostic run gives the same outlet flow as the As = 900, Au = 1 run, Au does not matter here; if it
                   differs, Au matters and the reason needs to be found. Which f? The rule, printed and saved in the manifest:
                       take the f-screen runs whose flood-week volume ratio stays between 0.6 and 1.6 (a wide sanity band;
                       change it with --guard-lo / --guard-hi) and use the LOWEST F_MULT among them. If none qualifies, step 3
                       stops and tells you; then choose yourself with --an-f.
                   --an-f 0.15 picks the f yourself. The run at that f with As = 1 is the f-screen run itself, so no run
                   is repeated (if you give an f that is not in the screen, one extra f run is added for it).

Total: 12 runs (1 + 6 + 5), each about 32 to 42 minutes of tRIBS plus about 1 minute of notebook build, so about 6.5 to 8.5 hours.

What each run is scored on (all read from files; nothing is judged for you):
    storms  flood week, 2nd January storm, Feb 18-19, Feb 27-28: volume ratio, peak error %, arrival error h (score_run_sycamore.py)
    tails   the two scorer tail windows (Jan 26 - Feb 2, Mar 3 - 11): volume ratio and recession-rate ratio, plus THIS script's
            own tail numbers, which answer the question directly:
              model and gauge mean flow (m3/s) in the 24 h blocks that start 7, 14 and 21 days after the gauge's flood peak
              the number of 24 h blocks between Jan 26 and Feb 18 whose mean flow is inside the 0.3 to 1 m3/s band
              (--band-lo / --band-hi), for the model and for the gauge
    dry months  Oct 1 - Dec 31 false flow (model volume in M m3, days above 0.13 m3/s, highest flow), and the Oct - Mar total
            volume ratio. The gauge is essentially dry in Oct - Dec, so any model flow there is false flow.
    process runoff type (the Hsrf share = infiltration-excess), the water table (MDGW, depth below the surface in mm: first, shallowest, deepest), and an
            indicative number for how much the depth decay of conductivity still allows lateral flow at the water table:
            exp(-f x depth), median over the basin, with f taken from the scaled f grid and the depth from the water-table
            grid. It ASSUMES f is per mm and the depth is in mm. It is a guide to how far the f change moves the answer, not a result.
    tail_ok / storms_ok in the table are two yes/no flags: tail_ok = the model's mean flow is inside the band on all of the
            7, 14, 21 day blocks; storms_ok = your anchor rule (flood-week volume ratio 0.90-1.12 and Feb 27-28 ratio 0.90-1.10,
            written for W2, so read it on W4 as a guide only).

WHAT THIS DOES NOT DO. Ks is NOT re-tuned: a lower f changes the storms, so some runs will leave the anchor band, and that is
expected here. Re-balancing Ks along the Ks-f ridge is the next step (the f x As grid and Ks re-compensation).
As and Au are set separately: the notebook's ANISOTROPY is As (saturated) and its ANISOTROPY_U is Au (unsaturated; None means
"same as ANISOTROPY"). The notebook needs the ANISOTROPY_U line; the start-up check says so if it is missing.
What As and Au each do inside tRIBS is not written in the tRIBS manual (it only lists them); As is taken to be the saturated
(groundwater) lateral path and Au the unsaturated one, from their names and the project's parameter reference.

How a run is made (copied from run_screen_soil_sycamore.py, which has been working):
    1. Copies the Settings cell of the Make_Sycamore_Model notebook, changes only what this run changes (plus the run names and
       the window), and runs the whole notebook headlessly with nbclient. The notebook file is never changed.
    2. Checks the build: the new .in file must match the baseline's .in in every keyword except the run-specific ones, and every
       grid in data/model/soil must equal its SOLUS grid times the multiplier (and soils.sdt must carry As and Au).
    3. Runs tRIBS (one at a time, in the foreground) and checks that it finished.
    4. Scores the run with score_run_sycamore.py (window W4) and computes the tail numbers above.
    5. Writes everything to tail_screen_w4_manifest.json, so a stopped batch can be resumed by running the same command.
Nothing in results/ is ever overwritten: a planned run whose name already has files in results/ stops the start-up check.

Run it from workspaces/Sycamore/lab/. Typical use (foreground only, keep the terminal open: the progress line printed every
5 minutes also keeps the Codespace awake):
    python -u run_tail_screen_sycamore.py --dry-run                 # checks only, builds nothing
    python -u run_tail_screen_sycamore.py --limit 1 --build-only    # build the baseline and check it, no tRIBS
    python -u run_tail_screen_sycamore.py --limit 2                 # baseline + the extreme f run, then look at the lines
    python -u run_tail_screen_sycamore.py                           # everything; run it again to resume after any stop
    python -u run_tail_screen_sycamore.py --steps base,f            # stop before the anisotropy screen
    python -u run_tail_screen_sycamore.py --steps an --an-f 0.15    # only the anisotropy screen, at f x0.15
    python -u run_tail_screen_sycamore.py --summary-only            # just print the tables from what is finished
    python -u run_tail_screen_sycamore.py --redo-decay              # fill in the f / decay columns of finished runs from the grids on disk (no run repeated)
    python -u run_tail_screen_sycamore.py --extra-runs 0.5874:300,0.15:300,0.02:200,0.02:500
                                                                    # extra As runs (Au = 1) at chosen f values, written f:As; only these are made

Files it makes (all in lab/): tail_screen_w4_manifest.json, tail_screen_w4.lock (only while running), tail_screen_w4_logs/,
score_out_tail_screen_w4/scores_tail_screen_w4.csv, tail_screen_w4_summary.csv, and for each run <run>.in,
simulation_<run>.log and results/<run>*. A W4 run writes more output than a W2 run: the start-up check needs 4 GB free
(--min-free-gb), and each finished run prints its output size.

Do NOT run the notebook (or any cell that rebuilds data/model/) while this runs: every build rewrites lab/data/model/, and all the
runs share it. If the connection drops, run the same command again: finished runs are skipped.
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
import pandas as pd

# --------------------------------------------------------------------------------------------
# What the screen does
# --------------------------------------------------------------------------------------------
WIN = dict(name="W4", start="10/01/2016/00/00", hours=4368, label="oct_to_mar", gw=0.95)   # label = the RUN_CHOICE text
RUN_START_TS = "2016-10-01 00:00"                     # the same start as WIN["start"], in pandas form

# The draw 12 / k10 setting (everything that is not being screened).
BASE = {
    "GW_FRACTION": WIN["gw"], "BEDROCK_ADD_M": 1.0, "KS_MULT": 0.187, "PSIB_MULT": 0.2939, "F_MULT": 0.5874,
    "THETAS_MULT": 1, "ANISOTROPY": 1, "ANISOTROPY_U": 1, "M_MULT": 1,
    "KINEMVELCOEF": 6.39, "FLOWEXP": 0.2973, "CHANNELROUGHNESS": 0.067, "CHANNELWIDTHCOEFF": 2.33,
    "OPTSNOW": 1, "OPTPERCOLATION": 0, "CHANNELCONDUCTIVITY": 5,
    "REBUILD_SOIL": False,
}
F_VALUES = [0.40, 0.25, 0.15, 0.08, 0.04, 0.02]       # F_MULT values of the f screen (baseline 0.5874)
AN_VALUES = [10, 100, 300, 900]                       # ANISOTROPY (= As) values of the screen, Au kept at 1 (baseline 1); the largest is also run with Au = As
ROUTING_KWS = ["KINEMVELCOEF", "FLOWEXP", "CHANNELROUGHNESS", "CHANNELWIDTHCOEFF"]
RUN_PREFIX = "W4_d012_k10"
BASELINE_ID = f"{RUN_PREFIX}_base_nocl"

GUARD_LO, GUARD_HI = 0.6, 1.6        # wide sanity band on the flood-week volume ratio for choosing the f of step 3
BAND_LO, BAND_HI = 0.3, 1.0          # the gauge's slow-flow band, m3/s
TAIL_OFFSETS_D = (7, 14, 21)         # days after the gauge's flood peak
TAIL_START, TAIL_END = "2017-01-26 00:00", "2017-02-18 00:00"   # the 24 h blocks used for the days-in-band count
FLOOD_WEEK = ("2017-01-18 00:00", "2017-01-25 00:00")          # where the gauge's flood peak is looked for
ANCHOR_FW = (0.90, 1.12)             # your anchor rule, flood-week volume ratio (written for W1/W2; a guide on W4)
ANCHOR_F27 = (0.90, 1.10)            # your anchor rule, Feb 27-28 volume ratio (written for W2; a guide on W4)

NOTEBOOK_CANDIDATES = ["Make_Sycamore_Model.ipynb", "Make_Sycamore_Model_clean.ipynb"]
SCORER = "score_run_sycamore.py"
TRIBS_DEFAULT = "/opt/tRIBS/build/tRIBS"
INIT_DATA = "../init_data"           # SOLUS_DIR in the notebook
OBS_FILE = "obs/sycamore_09510200_obs_flow_WY2017.csv"
MODEL_SOIL = "data/model/soil"
LOG_DIR = "tail_screen_w4_logs"
SCORE_DIR = "score_out_tail_screen_w4"
SCORE_CSV = "scores_tail_screen_w4.csv"
SUMMARY_CSV = "tail_screen_w4_summary.csv"
MANIFEST = "tail_screen_w4_manifest.json"
LOCK = "tail_screen_w4.lock"
MIN_PER_1000H = 7.3          # tRIBS speed on this machine, minutes per 1,000 model hours, fast end (a W4 run took about 32 min = 7.3)
MIN_PER_1000H_SLOW = 9.5     # slow end (the slowest W2 runs were about 10 per 1,000 h; allow for some slowing on the lowest-f runs)
HEARTBEAT_MIN = 5.0
ROWS_PER_MODEL_HOUR = 16     # the outlet file gets 16 rows per model hour
F_FLOOR = 1e-4               # the notebook keeps f at or above this
THETAS_CAP = 0.99            # the notebook keeps thetaS and porosity below this

# Score columns taken from the scorer's CSV. (short name, scorer column)
EVENTS = [("fw", "flood_week", "FLOOD WEEK"), ("s2", "jan_second_storm", "2ND JAN STORM"),
          ("f18", "feb18_19", "FEB 18-19"), ("f27", "feb27_28", "FEB 27-28")]
SCORE_COLS = []
for _s, _c, _t in EVENTS:
    SCORE_COLS += [f"{_c}__vol_ratio", f"{_c}__peak_err_pct", f"{_c}__arrival_err_h"]
SCORE_COLS += ["tail_jan26_feb2__tail_vol_ratio", "tail_mar3_mar11__tail_vol_ratio",
               "tail_jan26_feb2__recession_rate_ratio", "tail_mar3_mar11__recession_rate_ratio",
               "span__vol_ratio", "span__kge2012", "flood_week__kge2012", "oct_mar__vol_ratio",
               "octdec__vol_model_Mm3", "octdec__vol_gauge_Mm3", "octdec__days_above_model", "octdec__days_above_gauge",
               "octdec__max_model_m3s", "octdec__max_gauge_m3s",
               "proc__hsrf_share", "proc__mdgw_first_mm", "proc__mdgw_min_mm",
               "proc__mdgw_max_mm", "proc__snow_max_cm", "proc__channel_loss_Mm3"]

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


class RunError(Exception):
    """A problem with one run. The batch records it and moves on to the next run."""


class Abort(Exception):
    """A problem with the whole batch (competing job, no disk). The batch stops."""


# --------------------------------------------------------------------------------------------
# Names and settings for each job
# --------------------------------------------------------------------------------------------
def vtag(x) -> str:
    """0.15 -> '0p15', 0.5874 -> '0p5874', 100 -> '100'."""
    return f"{x:g}".replace(".", "p")


def au_of(s: dict):
    """The unsaturated anisotropy ratio a settings dict gives: ANISOTROPY_U, or ANISOTROPY when ANISOTROPY_U is None."""
    return s["ANISOTROPY"] if s.get("ANISOTROPY_U") is None else s["ANISOTROPY_U"]


def make_base_job() -> dict:
    return dict(step="base", changes={}, run_id=BASELINE_ID, f=BASE["F_MULT"], an=BASE["ANISOTROPY"], au=au_of(BASE))


def f_job(v: float) -> dict:
    return dict(step="f", changes={"F_MULT": v}, run_id=f"{RUN_PREFIX}_f{vtag(v)}_nocl", f=v, an=BASE["ANISOTROPY"], au=au_of(BASE))


def an_job(fv: float, a) -> dict:
    """As = a, Au stays at the baseline (1)."""
    return dict(step="an", changes={"F_MULT": fv, "ANISOTROPY": a}, run_id=f"{RUN_PREFIX}_f{vtag(fv)}_as{vtag(a)}_nocl",
                f=fv, an=a, au=au_of(BASE))


def au_job(fv: float, a) -> dict:
    """The diagnostic run: As = Au = a at this f. 'partner' is the run that differs from it only in Au (As = a, Au = 1)."""
    return dict(step="an", changes={"F_MULT": fv, "ANISOTROPY": a, "ANISOTROPY_U": a},
                run_id=f"{RUN_PREFIX}_f{vtag(fv)}_as{vtag(a)}au{vtag(a)}_nocl", f=fv, an=a, au=a, partner=an_job(fv, a)["run_id"])


def an_base_id(fv: float) -> str:
    """The run that is 'anisotropy 1' at this f: the baseline itself, or the f-screen run."""
    return BASELINE_ID if abs(fv - BASE["F_MULT"]) < 1e-12 else f_job(fv)["run_id"]


def f_jobs() -> list:
    """Lowest f first (the extreme runs first, as a smoke test), then the rest from high to low."""
    vals = sorted(F_VALUES)
    return [f_job(v) for v in [vals[0]] + sorted(vals[1:], reverse=True)]


def an_jobs(fv: float) -> list:
    """Largest As first, then the As = Au diagnostic run right after it (so the Au question is answered early), then the rest ascending."""
    vals = sorted(AN_VALUES)
    return [an_job(fv, vals[-1]), au_job(fv, vals[-1])] + [an_job(fv, a) for a in vals[:-1]]


def settings_for(job: dict) -> dict:
    """Every NAME = value line of the Settings cell that this run sets."""
    s = dict(BASE)
    s.update(job["changes"])
    s.update(RUN_CHOICE=WIN["label"], INPUT_FILE=f"{job['run_id']}.in", RESULTS_NAME=job["run_id"])
    return s


RECORD_KEYS = ("GW_FRACTION", "BEDROCK_ADD_M", "KS_MULT", "PSIB_MULT", "F_MULT", "M_MULT", "THETAS_MULT", "ANISOTROPY", "ANISOTROPY_U", "KINEMVELCOEF",
               "FLOWEXP", "CHANNELROUGHNESS", "CHANNELWIDTHCOEFF", "OPTSNOW", "OPTPERCOLATION", "CHANNELCONDUCTIVITY")


def settings_record(job: dict) -> dict:
    """What the manifest remembers about the settings, to notice a changed setting on resume."""
    s = settings_for(job)
    rec = {k: s[k] for k in RECORD_KEYS}
    rec.update(START_DATE=WIN["start"], RUN_HOURS=WIN["hours"])
    return rec


def expected_files(s):
    """Names the notebook gives the bedrock and water-table grids for these settings."""
    bed = s["BEDROCK_ADD_M"]
    gw_name = f"InitGW_{int(round(s['GW_FRACTION'] * 100))}pct_plus{bed:g}m_mm.asc"
    bed_name = "bedrock_depth_m.asc" if bed == 0 else f"bedrock_depth_plus{bed:g}m.asc"
    return gw_name, bed_name


def describe(job: dict) -> str:
    if job["step"] == "base":
        return "baseline (no change)"
    if job["step"] == "f":
        return f"f x{job['f']:g}"
    if job["au"] == job["an"]:
        return f"f x{job['f']:g}, As = Au = {job['an']:g}"
    return f"f x{job['f']:g}, As {job['an']:g}, Au {job['au']:g}"


def timeout_for(args) -> float:
    if args.timeout_min:
        return args.timeout_min
    return max(30.0, 3.5 * WIN["hours"] * MIN_PER_1000H / 1000.0)          # W4: about 119 min, 3.5 times the expected time


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
                                 display_name="tail_screen", language="python")
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
        f"Au setting: ANISOTROPY_U = {s['ANISOTROPY_U']!r} (None = same as ANISOTROPY)",
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
                       + "\n  (If these are the thetaS or anisotropy lines, or the 'Au setting' line, the notebook in lab/ is not the one with ANISOTROPY_U: replace it with the new one.)")
    if re.search(r"^PROBLEM:", text, re.M):
        raise RunError("The notebook printed a PROBLEM line (see the build output file).")
    warns = [ln.strip() for ln in text.splitlines() if re.match(r"^CHECK\s", ln)]
    m = re.search(r"cells held at the [\d.]+ cap: thetaS (\d+), porosity (\d+)", text)
    if m and (int(m.group(1)) or int(m.group(2))):
        warns.append(f"thetaS cells held at the {THETAS_CAP:g} cap: {m.group(1)}, porosity {m.group(2)} "
                     f"(the scaling was limited in those cells)")
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


def check_in_file(new_in: Path, ref_in: Path | None, job: dict) -> None:
    """The new .in must carry the values this run's settings call for (window, output name, water-table and bedrock files,
    the routing keywords, snow and channel-loss settings). If a reference .in is given (the baseline's), every other keyword
    must equal the reference. The soil knobs do not appear in the .in file; check_grids covers them."""
    s = settings_for(job)
    gw_name, bed_name = expected_files(s)
    new = read_in_file(new_in)
    expect = {"RUNTIME": WIN["hours"], "STARTDATE": WIN["start"], "OUTFILENAME": f"results/{job['run_id']}",
              "GWATERFILE": f"{MODEL_SOIL}/{gw_name}", "BEDROCKFILE": f"{MODEL_SOIL}/{bed_name}",
              "OPTPERCOLATION": s["OPTPERCOLATION"], "OPTSNOW": s["OPTSNOW"], "CHANNELCONDUCTIVITY": s["CHANNELCONDUCTIVITY"]}
    for kw in ROUTING_KWS:
        expect[kw] = s[kw]
    problems = []
    for key, want in expect.items():
        n = new.get(key)
        if n is None or not same_value(n, want):
            problems.append(f"{key}: written {n!r}, expected {want!r}")
    if ref_in is not None:
        ref = read_in_file(ref_in)
        for key in sorted(set(ref) | set(new)):
            if key in expect:
                continue
            b, n = ref.get(key), new.get(key)
            if b is None or n is None or not same_value(b, n):
                problems.append(f"{key}: reference {b!r} but this run {n!r}  (this keyword may not change)")
    if problems:
        raise RunError("The .in file differs from what it should be:\n  " + "\n  ".join(problems))


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


def grid_axes(h: dict, shape: tuple):
    """Cell-centre x (west to east) and y (north to south: row 0 is the top row) of an ESRI ASCII grid, plus the cell size."""
    cs = h.get("cellsize")
    if cs is None:
        raise ValueError("the grid header has no cellsize")
    x0 = h["xllcenter"] if "xllcenter" in h else h.get("xllcorner", 0.0) + cs / 2.0
    y0 = h["yllcenter"] if "yllcenter" in h else h.get("yllcorner", 0.0) + cs / 2.0
    nrows, ncols = shape
    return x0 + cs * np.arange(ncols), y0 + cs * (nrows - 1 - np.arange(nrows)), cs


def decay_diagnostic(fd_file: Path, gw: np.ndarray, hg: dict, unit: str) -> dict:
    """How much the depth decay exp(-f x depth) still allows at the water table, median over the basin. ASSUMES f is per mm and
    that the water-table grid is in mm (or in m, converted). An indicative number, not a result.
    The f grid is read at the centre of every water-table cell: directly when the two grids are the same size and place, otherwise
    the nearest f cell. Returns a dict; if it cannot be computed the dict is {"error": reason}. It never raises."""
    try:
        if not Path(fd_file).is_file():
            return {"error": f"f grid not found: {fd_file}"}
        hf, f = read_asc(Path(fd_file))
        f_nd, g_nd = hf.get("nodata_value", -9999), hg.get("nodata_value", -9999)
        same = (f.shape == gw.shape and all(abs(hf.get(k, 0) - hg.get(k, 0)) < 1e-6
                                            for k in ("xllcorner", "yllcorner", "xllcenter", "yllcenter", "cellsize")))
        if same:
            ok = (f != f_nd) & (gw != g_nd) & (gw > 0) & (f > 0)
            fv, depth = f[ok], gw[ok]
            how = "same grid"
        else:
            gx, gy, _ = grid_axes(hg, gw.shape)
            fx, fy, fcs = grid_axes(hf, f.shape)
            rr, cc = np.nonzero((gw != g_nd) & (gw > 0))
            col = np.rint((gx[cc] - fx[0]) / fcs).astype(int)
            row = np.rint((fy[0] - gy[rr]) / fcs).astype(int)
            inside = (col >= 0) & (col < f.shape[1]) & (row >= 0) & (row < f.shape[0])
            rr, cc, col, row = rr[inside], cc[inside], col[inside], row[inside]
            fsamp = f[row, col]
            keep = (fsamp != f_nd) & (fsamp > 0)
            fv, depth = fsamp[keep], gw[rr[keep], cc[keep]]
            how = f"nearest f cell (f grid {f.shape[0]}x{f.shape[1]}, water-table grid {gw.shape[0]}x{gw.shape[1]})"
        if fv.size < 10:
            return {"error": f"only {fv.size} cells where both the f grid and the water-table grid have values ({how})"}
        depth_mm = depth * (1.0 if unit == "mm" else 1000.0)
        d = np.exp(-fv * depth_mm)
        return {"f_median": float(np.median(fv)), "gw_depth_median_mm": float(np.median(depth_mm)),
                "decay_median": float(np.median(d)), "n_cells": int(fv.size), "sampled": how}
    except Exception as exc:                                           # an extra diagnostic must never fail a run
        return {"error": f"{type(exc).__name__}: {exc}"}


def decay_ok(diag) -> bool:
    return bool(diag) and "error" not in diag


def redo_decay(man: dict, lab: Path) -> int:
    """Recompute the f / decay diagnostic of every finished run from the grids already in data/model/soil (no run is repeated).
    Returns the number of runs that got a number."""
    src_fd = read_gdf(lab / INIT_DATA / "scgrid.gdf")["FD"]
    stem = os.path.basename(src_fd)
    got = 0
    for rid, rec in man["runs"].items():
        s = rec.get("settings")
        if rec.get("status") not in ("done", "unscored") or not s:
            continue
        mult = s["F_MULT"]
        gw_name, _ = expected_files(s)
        if mult == 1:
            fd_file = lab / f"{src_fd}.asc"
        else:
            cands = sorted((lab / MODEL_SOIL).glob(f"{stem}*_x{mult:g}.asc"))
            if len(cands) != 1:
                rec["decay"] = {"error": f"expected one scaled f grid like {stem}*_x{mult:g}.asc in {MODEL_SOIL}, found {len(cands)}"}
                log(f"{rid}: decay not computed: {rec['decay']['error']}")
                continue
            fd_file = cands[0]
        try:
            hg, gw = read_asc(lab / MODEL_SOIL / gw_name)
        except Exception as exc:
            rec["decay"] = {"error": f"{type(exc).__name__}: {exc}"}
            log(f"{rid}: decay not computed: {rec['decay']['error']}")
            continue
        valid = gw[(gw != hg.get("nodata_value", -9999)) & (gw > 0)]
        unit = "mm" if valid.size and float(np.median(valid)) > 50 else "m"     # depths here are about 1,500 mm; in metres they would be about 1.5
        diag = decay_diagnostic(fd_file, gw, hg, unit)
        rec["decay"] = diag
        if decay_ok(diag):
            got += 1
            log(f"{rid}: f median {diag['f_median']:.4g}, water-table depth median {diag['gw_depth_median_mm']:.0f} mm, decay {diag['decay_median']:.2e} ({diag['sampled']})")
        else:
            log(f"{rid}: decay not computed: {diag['error']}")
    save_manifest(man)
    return got


def check_grids(lab: Path, job: dict, t_build_start: float) -> tuple:
    """Bedrock, water-table and every soil grid / the soil table in data/model/soil must be fresh and correct.
    Returns (a text note, the decay diagnostic dict; it holds an 'error' entry if the diagnostic could not be computed)."""
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
    # the soil grids: every row of scgrid.gdf against its SOLUS grid
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
    for var in mults:
        check_one_grid(lab, var, model_rows[var], src_rows[var], mults[var], t_build_start)
    # the soil table: As and Au
    sdt = soil / "soils.sdt"
    if not sdt.is_file() or sdt.stat().st_mtime < t_build_start - 2:
        raise RunError("data/model/soil/soils.sdt is missing or older than this build.")
    table = [ln.split(",") for ln in sdt.read_text().splitlines()[1:] if ln.strip()]
    if not table or any(len(r) != 12 for r in table):
        raise RunError("soils.sdt does not have 12 comma-separated columns on every row.")
    for r in table:
        if float(r[7]) != float(s["ANISOTROPY"]) or float(r[8]) != float(au_of(s)):
            raise RunError(f"soils.sdt row {r[0]} has As = {r[7]}, Au = {r[8]}; this run calls for As = {s['ANISOTROPY']:g}, Au = {au_of(s):g}.")
    diag = decay_diagnostic(lab / f"{model_rows['FD']}.asc", gw, hg, unit)
    note = (f"Ks x{s['KS_MULT']:g}, PsiB x{s['PSIB_MULT']:g}, f x{s['F_MULT']:g}, m x{s['M_MULT']:g}, thetaS x{s['THETAS_MULT']:g}, "
            f"As = {s['ANISOTROPY']:g}, Au = {au_of(s):g}; bedrock +{s['BEDROCK_ADD_M']:g} m; water table = {s['GW_FRACTION']:g} x bedrock "
            f"({unit}), median depth {np.median(gw[ok]):.0f} {unit}")
    if decay_ok(diag):
        note += (f"; f grid median {diag['f_median']:.4g} (per mm assumed), decay exp(-f x depth) at the water table "
                 f"{diag['decay_median']:.2e}")
    else:
        note += f"; f decay diagnostic skipped ({diag['error']})"
    return note, diag


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
    terminal activity, so a Codespace (30 minute idle timeout) does not shut down in the middle of a run.
    """
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


def run_scorer(rid: str, scorer: str, score_dir: Path, log_dir: Path) -> dict:
    """Score one run; returns its whole row of the scorer's CSV (all values as text)."""
    cmd = [sys.executable, scorer, "score", "--run", rid, "--window", WIN["name"], "--results-dir", "results",
           "--out-dir", str(score_dir), "--csv-name", SCORE_CSV, "--write-checksum"]
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
    out = {}
    for c in SCORE_COLS:
        v = fnum(row, c)
        if math.isfinite(v):
            out[c] = round(v, 4)
    return out


def one_line(scores: dict, tail: dict | None = None) -> str:
    names = [("flood_week__vol_ratio", "flood vol"), ("flood_week__peak_err_pct", "peak %"),
             ("flood_week__arrival_err_h", "arrival h"), ("feb27_28__vol_ratio", "Feb27-28 vol"),
             ("tail_jan26_feb2__tail_vol_ratio", "tail1"), ("tail_mar3_mar11__tail_vol_ratio", "tail2"),
             ("proc__hsrf_share", "Hsrf share"), ("octdec__vol_model_Mm3", "Oct-Dec M m3"), ("oct_mar__vol_ratio", "Oct-Mar vol")]
    text = ", ".join(f"{label} {scores[c]}" for c, label in names if c in scores)
    if tail and tail.get("model_days_in_band") is not None:
        text += (f"; slow flow in band {tail['model_days_in_band']} of {tail['days_total']} days "
                 f"(gauge {tail['gauge_days_in_band']})")
        d14 = tail.get("model_d14_m3s")
        if d14 is not None and math.isfinite(d14):
            text += f", day-14 flow {d14:.3f} m3/s (gauge {tail['gauge_d14_m3s']:.3f})"
    return text


# --------------------------------------------------------------------------------------------
# The tail numbers (read from the outlet file and the gauge file)
# --------------------------------------------------------------------------------------------
def read_qout(path: Path, start) -> pd.Series:
    """Routed outlet flow (m3/s) as a time series (comma-delimited, header line, columns Time_hr,Q), like the scorer reads it."""
    with open(path) as fh:
        fh.readline()
        second = fh.readline()
    if "," in second:
        raw = pd.read_csv(path, sep=",", header=0)
    else:
        raw = pd.read_csv(path, sep=r"\s+", skiprows=1, header=None)
    raw = raw.iloc[:, :2].copy()
    raw.columns = ["Time_hr", "Q"]
    raw = raw.apply(pd.to_numeric, errors="coerce").dropna()
    s = pd.Series(raw["Q"].to_numpy(float), index=pd.Timestamp(start) + pd.to_timedelta(raw["Time_hr"].to_numpy(), unit="h"))
    return s[~s.index.duplicated(keep="last")].sort_index()


def read_obs(path: Path) -> pd.Series:
    obs = pd.read_csv(path, parse_dates=["datetime"]).set_index("datetime")
    s = pd.to_numeric(obs["q_cms"], errors="coerce")
    return s[~s.index.duplicated(keep="last")].sort_index()


def block_mean(q: pd.Series, a, b, min_hours: int = 18) -> float:
    """Mean flow over [a, b), from hourly means. NaN if fewer than min_hours of the hours have data."""
    s = q[(q.index >= pd.Timestamp(a)) & (q.index < pd.Timestamp(b))].dropna()
    if s.empty:
        return float("nan")
    hourly = s.resample("60min").mean().dropna()
    return float(hourly.mean()) if len(hourly) >= min_hours else float("nan")


def tail_metrics(model_q: pd.Series, gauge_q: pd.Series, lo: float = BAND_LO, hi: float = BAND_HI) -> dict:
    """The tail numbers. Blocks are 24 h long. d7, d14, d21 start that many days after the gauge's flood peak (floored to the
    hour); the days-in-band count uses every 24 h block from TAIL_START to TAIL_END."""
    out = {"band_lo": lo, "band_hi": hi}
    fw = gauge_q[(gauge_q.index >= pd.Timestamp(FLOOD_WEEK[0])) & (gauge_q.index < pd.Timestamp(FLOOD_WEEK[1]))].dropna()
    peak_t = fw.idxmax() if len(fw) else None
    out["gauge_peak_time"] = str(peak_t) if peak_t is not None else ""
    t0, t1 = pd.Timestamp(TAIL_START), pd.Timestamp(TAIL_END)
    day = pd.Timedelta(hours=24)
    for d in TAIL_OFFSETS_D:
        m = g = float("nan")
        out[f"d{d}_start"] = ""
        if peak_t is not None:
            a = (peak_t + pd.Timedelta(days=d)).floor("h")
            out[f"d{d}_start"] = str(a)
            if a >= t0 and a + day <= t1:
                m, g = block_mean(model_q, a, a + day), block_mean(gauge_q, a, a + day)
        out[f"model_d{d}_m3s"], out[f"gauge_d{d}_m3s"] = m, g
    blocks = list(pd.date_range(t0, t1, freq="24h", inclusive="left"))
    mm = [block_mean(model_q, a, a + day) for a in blocks]
    gg = [block_mean(gauge_q, a, a + day) for a in blocks]
    out["days_total"] = len(blocks)
    out["model_days_in_band"] = int(sum(1 for v in mm if math.isfinite(v) and lo <= v <= hi))
    out["gauge_days_in_band"] = int(sum(1 for v in gg if math.isfinite(v) and lo <= v <= hi))
    out["model_days_valid"] = int(sum(1 for v in mm if math.isfinite(v)))
    out["gauge_days_valid"] = int(sum(1 for v in gg if math.isfinite(v)))
    vals = [out[f"model_d{d}_m3s"] for d in TAIL_OFFSETS_D]
    ok_vals = [v for v in vals if math.isfinite(v)]
    out["tail_ok"] = bool(len(ok_vals) >= 2 and all(lo <= v <= hi for v in ok_vals))
    return out


def storms_ok(row: dict) -> bool | None:
    fw, f27 = fnum(row, "flood_week__vol_ratio"), fnum(row, "feb27_28__vol_ratio")
    if not (math.isfinite(fw) and math.isfinite(f27)):
        return None
    return bool(ANCHOR_FW[0] <= fw <= ANCHOR_FW[1] and ANCHOR_F27[0] <= f27 <= ANCHOR_F27[1])


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


def load_scores(lab: Path) -> dict:
    p = lab / SCORE_DIR / SCORE_CSV
    out = {}
    if p.is_file():
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh):
                out[r["run_id"]] = r
    return out


def is_finished(man: dict, job: dict, args=None) -> bool:
    rec = man["runs"].get(job["run_id"], {})
    redo = getattr(args, "redo_ids", set()) if args is not None else set()
    return rec.get("status") == "done" and rec.get("settings") == settings_record(job) and job["run_id"] not in redo


# --------------------------------------------------------------------------------------------
# Preflight
# --------------------------------------------------------------------------------------------
def pick_notebook(lab: Path, given: str | None) -> tuple:
    """(notebook name, note). With no --notebook, Make_Sycamore_Model.ipynb if it exists, else Make_Sycamore_Model_clean.ipynb."""
    if given:
        return given, ""
    present = [n for n in NOTEBOOK_CANDIDATES if (lab / n).is_file()]
    if not present:
        return NOTEBOOK_CANDIDATES[0], ""
    note = ""
    if len(present) > 1:
        note = (f"both {present[0]} and {present[1]} are in lab/; using {present[0]}. "
                f"They should have the same Settings cell. Use --notebook to choose the other.")
    return present[0], note


def preflight(args, lab: Path, jobs: list, nb_note: str, man: dict) -> list:
    """Checks before anything is built. Returns a list of (level, message)."""
    out = []

    def add(level, msg):
        out.append((level, msg))

    for f in (args.notebook, args.scorer):
        add("OK" if (lab / f).is_file() else "FAIL", f"file in lab/: {f}")
    if nb_note:
        add("WARN", nb_note)
    obs = lab / INIT_DATA / OBS_FILE
    add("OK" if obs.is_file() else "FAIL", f"gauge file: {INIT_DATA}/{OBS_FILE}")
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
        for name in ("M_MULT", "THETAS_MULT", "ANISOTROPY", "ANISOTROPY_U"):
            if not re.search(rf"^{name}[ \t]*=", src, re.M):
                old.append(f"no line '{name} = ...' in the Settings cell")
        if not re.search(r"^SOIL_CONSTANTS[ \t]*=.*'As':\s*ANISOTROPY\b[^_].*'Au':[^}]*ANISOTROPY_U", src, re.M):
            old.append("SOIL_CONSTANTS does not take Au from ANISOTROPY_U")
        all_code = "\n".join(c.source for c in nb.cells if c.cell_type == "code")
        if "THETAS_CAP" not in all_code or "theta*_s" not in all_code:
            old.append("the thetaS cap or the theta*_s check is missing")
        if old:
            add("FAIL", f"{args.notebook} is not the notebook with ANISOTROPY_U ({'; '.join(old)}). Replace it with the new Make_Sycamore_Model.ipynb, then run again.")
        else:
            test_jobs = list(jobs)
            if not any(j["step"] == "an" for j in test_jobs):
                test_jobs += [an_job(min(F_VALUES), max(AN_VALUES)), au_job(min(F_VALUES), max(AN_VALUES))]   # samples, so both anisotropy patches are tested
            for j in test_jobs:
                patched_source(src, j)
            add("OK", f"Settings cell (cell {idx}) has M_MULT, THETAS_MULT, ANISOTROPY and ANISOTROPY_U, the cap and theta*_s cells are there, "
                      f"and it patches for all {len(test_jobs)} runs checked")
    except ImportError:
        add("FAIL", "nbformat is not installed: pip install nbformat nbclient ipykernel")
    except Exception as exc:
        add("FAIL", f"Settings cell check: {exc}")
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
    # nothing may be overwritten: a new run's name must not already have files in results/
    results_dir = lab / "results"
    for j in jobs:
        stale = list(results_dir.glob(j["run_id"] + "*")) if results_dir.is_dir() else []
        if stale and j["run_id"] not in man["runs"]:
            add("FAIL", f"results/ already has {len(stale)} file(s) named {j['run_id']}*; this script never overwrites existing results")
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


def compute_tail(lab: Path, rid: str, args) -> dict:
    model_q = read_qout(lab / "results" / f"{rid}_Outlet.qout", RUN_START_TS)
    gauge_q = read_obs(lab / INIT_DATA / OBS_FILE)
    return tail_metrics(model_q, gauge_q, args.band_lo, args.band_hi)


def do_run(job: dict, args, lab: Path, man: dict) -> None:
    rid, s = job["run_id"], settings_for(job)
    rec = man["runs"].setdefault(rid, {})
    log_dir, score_dir = lab / LOG_DIR, lab / SCORE_DIR
    log_dir.mkdir(exist_ok=True)
    score_dir.mkdir(exist_ok=True)
    cur = settings_record(job)
    ref_in = None if job["step"] == "base" else lab / f"{BASELINE_ID}.in"

    if rec.get("status") in ("done", "unscored") and rec.get("settings") != cur and rid not in args.redo_ids:
        raise RunError(f"{rid} is finished, but its recorded settings differ from the ones now planned "
                       f"({rec.get('settings')} -> {cur}). It is NOT rerun automatically; use --redo {rid} if you really want that.")
    if rec.get("status") == "done" and rid not in args.redo_ids:
        log(f"{rid}: already done, skipping")
        return
    if rec.get("status") == "unscored" and rid not in args.redo_ids:
        log(f"{rid}: ran earlier but was not scored; scoring again (tRIBS is not rerun)")
        rec["scores"] = show_scores(run_scorer(rid, str(lab / args.scorer), score_dir, log_dir))
        rec["tail"] = compute_tail(lab, rid, args)
        rec["status"] = "done"
        rec["finished"] = now_iso()
        save_manifest(man)
        return

    if job["step"] != "base" and (man["runs"].get(BASELINE_ID) or {}).get("status") != "done":
        raise Abort(f"The baseline run {BASELINE_ID} is not finished, and every other run is checked against its .in file. "
                    f"Look at the baseline's status above (it always runs first) and run the script again.")
    if ref_in is not None and not ref_in.is_file():
        raise Abort(f"The baseline .in file {ref_in.name} is missing from lab/.")

    free = disk_free_gb(lab)
    if free < args.min_free_gb:
        raise Abort(f"Only {free:.1f} GB free; need {args.min_free_gb:g} GB. Free some space and run again to resume.")
    comp = competing_tribs()
    if comp:
        raise Abort("Another tRIBS is running, so the batch stops (two at once is not allowed): " + comp[0])

    rec.clear()
    rec.update(status="building", started=now_iso(), step=job["step"], f=job["f"], an=job["an"], au=job["au"],
               partner=job.get("partner"), changes=job["changes"], settings=cur)
    save_manifest(man)
    remove_stale(lab, rid)
    log(f"=== {rid}: {describe(job)} on {WIN['name']} ({WIN['hours']} h from {WIN['start']}, water table {WIN['gw']:g} x bedrock) ===")

    t_build = time.time()
    text, warns = build_run(job, args, lab, log_dir)
    check_in_file(lab / s["INPUT_FILE"], ref_in, job)
    grid_note, diag = check_grids(lab, job, t_build)
    rec.update(build_min=round((time.time() - t_build) / 60, 2), warnings=warns, grid_check=grid_note, decay=diag, status="built")
    save_manifest(man)
    log(f"{rid}: build OK in {fmt_min(time.time() - t_build)}; the .in file checks out")
    log(f"{rid}: {grid_note}")
    for wmsg in warns:
        log(f"{rid}: notebook said: {wmsg}")
    if args.build_only:
        log(f"{rid}: --build-only, not running tRIBS")
        return

    install_term_handler()
    log(f"{rid}: starting tRIBS ({WIN['hours']} model hours, expect about {WIN['hours'] * MIN_PER_1000H / 1000:.0f}-{WIN['hours'] * MIN_PER_1000H_SLOW / 1000:.0f} min)")
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
    if job["step"] != "base" and bmd5 and outlet["md5"] == bmd5:
        rec["identical_to_baseline"] = True
        log(f"{rid}: WARNING: the outlet file is byte-identical to the baseline's. tRIBS ran, but this change had no effect on "
            f"the outlet flow. Either the knob has no pathway to the outlet here, or tRIBS does not read the setting.")
    partner = job.get("partner")
    if partner:
        pmd5 = ((man["runs"].get(partner) or {}).get("outlet") or {}).get("md5")
        if pmd5 is None:
            log(f"{rid}: its comparison run {partner} (the same As, Au = {BASE['ANISOTROPY_U']:g}) is not finished yet; the Au comparison is printed in the summary once both are done.")
        elif outlet["md5"] == pmd5:
            rec["identical_to_partner"] = True
            log(f"{rid}: the outlet file is byte-identical to {partner} (same As, Au = {BASE['ANISOTROPY_U']:g}). Au has NO effect on the outlet flow in this setup.")
        else:
            log(f"{rid}: the outlet file DIFFERS from {partner} (same As, Au = {BASE['ANISOTROPY_U']:g}), so Au does change the outlet flow here; see the Au diagnostic in the summary.")

    try:
        rec["scores"] = show_scores(run_scorer(rid, str(lab / args.scorer), score_dir, log_dir))
        rec["status"] = "done"
    except RunError as exc:
        rec["score_error"] = str(exc)
        log(f"{rid}: scoring failed (the run itself is fine): {exc}")
    try:
        rec["tail"] = compute_tail(lab, rid, args)
    except Exception as exc:                              # the tail numbers are extra; never fail a run over them
        rec["tail_error"] = f"{type(exc).__name__}: {exc}"
        log(f"{rid}: tail numbers could not be computed: {exc}")
    rec["finished"] = now_iso()
    save_manifest(man)
    if rec.get("scores"):
        log(f"{rid}: {one_line(rec['scores'], rec.get('tail'))}")


# --------------------------------------------------------------------------------------------
# Choosing the f for the anisotropy screen
# --------------------------------------------------------------------------------------------
def pick_an_f(man: dict, scores: dict, args) -> tuple:
    """Returns (info dict or None, reason). --an-f wins; then an earlier choice kept in the manifest (unless --repick);
    otherwise the rule: the lowest F_MULT among the f-screen runs whose flood-week volume ratio is inside the guard band."""
    if args.an_f is not None:
        return dict(value=args.an_f, source="--an-f"), ""
    if man.get("an_f") and not args.repick:
        return man["an_f"], ""
    cands, seen = [], []
    for v in F_VALUES:
        rid = f_job(v)["run_id"]
        sc = scores.get(rid)
        fw = fnum(sc, "flood_week__vol_ratio") if sc else float("nan")
        seen.append([v, round(fw, 4) if math.isfinite(fw) else None])
        if math.isfinite(fw) and args.guard_lo <= fw <= args.guard_hi:
            cands.append((v, fw))
    if not cands:
        return None, (f"no f-screen run has a flood-week volume ratio between {args.guard_lo:g} and {args.guard_hi:g} "
                      f"(seen: {seen}); choose with --an-f")
    v = min(c[0] for c in cands)
    return dict(value=v, source="rule", rule="lowest F_MULT with flood-week volume ratio inside the guard band",
                guard=[args.guard_lo, args.guard_hi], f_screen_flood_week_vol=seen), ""


# --------------------------------------------------------------------------------------------
# The tables and the summary CSV
# --------------------------------------------------------------------------------------------
def fm(x, nd=2, sign=False) -> str:
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "-"
    if abs(x) < 0.5 * 10 ** -nd:                       # avoid printing -0.00
        x = 0.0
    return f"{x:+.{nd}f}" if sign else f"{x:.{nd}f}"


def yn(v) -> str:
    return "-" if v is None else ("yes" if v else "no")


def sci(x) -> str:
    return "-" if x is None or not math.isfinite(x) else f"{x:.1e}"


def run_label(rec: dict) -> str:
    if rec.get("step") == "base":
        return "baseline"
    if rec.get("step") == "f":
        return f"f x{rec['f']:g}"
    an, au = rec["an"], rec.get("au", 1)
    if au == an:
        return f"f x{rec['f']:g} As=Au {an:g}"
    return f"f x{rec['f']:g} As {an:g}"


def au_diagnostic(man: dict, rows: list) -> None:
    """For each As = Au diagnostic run: set it beside the run that differs from it only in Au (As = same, Au = 1)."""
    by_id = {rid: (rec, sc, tail) for rid, rec, sc, tail in rows}
    pairs = [(rid, rec["partner"]) for rid, rec, _, _ in rows if rec.get("partner") and rec["partner"] in by_id]
    if not pairs:
        return
    log()
    log("  AU DIAGNOSTIC  (the two runs differ only in Au: does the unsaturated anisotropy ratio change anything?)")
    log("    " + f"{'run':<22}| {'fw vol':>6} {'f27 vol':>7} | " + " ".join(f"{'d' + str(d):>7}" for d in TAIL_OFFSETS_D) + f" | {'days in band':>12} | {'Oct-Dec M m3':>12} | outlet file")
    for diag_id, part_id in pairs:
        drec, dsc, dtail = by_id[diag_id]
        prec, psc, ptail = by_id[part_id]
        for rec, sc, tail, tag in ((prec, psc, ptail, f"As {prec['an']:g}, Au {prec.get('au', 1):g}"), (drec, dsc, dtail, f"As {drec['an']:g}, Au {drec['au']:g}")):
            if not sc:
                log("    " + f"{tag + ' (' + str(rec.get('status', 'not run')) + ')':<22}|")
                continue
            d = " ".join(f"{fm(tail[f'model_d{x}_m3s'], 3) if tail else '-':>7}" for x in TAIL_OFFSETS_D)
            band = (str(tail['model_days_in_band']) + '/' + str(tail['days_total'])) if tail else '-'
            log("    " + f"{tag:<22}| {fm(fnum(sc, 'flood_week__vol_ratio'), 3):>6} {fm(fnum(sc, 'feb27_28__vol_ratio'), 3):>7} | {d} | {band:>12} | "
                         f"{fm(fnum(sc, 'octdec__vol_model_Mm3'), 2):>12} | {'identical to the As-only run' if rec.get('identical_to_partner') else ''}")
        pm, dm = (prec.get("outlet") or {}).get("md5"), (drec.get("outlet") or {}).get("md5")
        if pm and dm and pm == dm:
            log(f"    -> byte-identical outlet files: Au has no effect on the outlet flow here (As {drec['an']:g}, f x{drec['f']:g}).")
        elif pm and dm:
            small = False
            if psc and dsc and ptail and dtail:
                dv = abs(fnum(dsc, "flood_week__vol_ratio") - fnum(psc, "flood_week__vol_ratio"))
                dd = abs(dtail["model_d14_m3s"] - ptail["model_d14_m3s"]) if (math.isfinite(dtail["model_d14_m3s"]) and math.isfinite(ptail["model_d14_m3s"])) else float("nan")
                small = (dv < 0.005) and (dd < 0.005 if math.isfinite(dd) else False) and dtail["model_days_in_band"] == ptail["model_days_in_band"]
            log("    -> the outlet files differ, so Au does change the flow. " + ("The numbers above differ by under 0.005 in flood-week volume ratio and day-14 flow, with the same days in band: a very small effect."
                                                                             if small else "Compare the two lines above; if the difference matters, the next question is why Au acts here."))
        else:
            log("    -> one of the two runs has no outlet file yet; run the script again to finish it.")


def print_summary(man: dict, lab: Path, args) -> None:
    scores = load_scores(lab)
    order = {"base": 0, "f": 1, "an": 2}
    recs = [(rid, rec) for rid, rec in man["runs"].items() if rec.get("step") in order]
    recs.sort(key=lambda t: (order[t[1]["step"]], -float(t[1].get("f", 0)), float(t[1].get("an", 1)), float(t[1].get("au", 1))))
    if not recs:
        log("No runs yet.")
        return
    base_rec = man["runs"].get(BASELINE_ID) or {}
    base_decay = (base_rec.get("decay") or {}).get("decay_median")
    rows, tails = [], {}
    for rid, rec in recs:
        sc = scores.get(rid)
        tail = rec.get("tail")
        if tail is None and sc and (lab / "results" / f"{rid}_Outlet.qout").is_file():
            try:
                tail = compute_tail(lab, rid, args)
            except Exception:
                tail = None
        tails[rid] = tail
        rows.append((rid, rec, sc, tail))

    log()
    log("TAIL SCREEN ON W4 (Oct 1 2016 - Apr 1 2017). Draw 12 / k10 baseline; each row changes only f (step 2) or f and the saturated anisotropy As (step 3; Au stays 1 except on the 'As=Au' diagnostic row).")
    log("  vol = model volume / gauge volume. pk = model peak vs gauge peak, %. arr = first-arrival error, hours (negative = early).")
    log(f"  storms_ok = your anchor rule: flood-week vol {ANCHOR_FW[0]:g}-{ANCHOR_FW[1]:g} and Feb 27-28 vol {ANCHOR_F27[0]:g}-{ANCHOR_F27[1]:g}. "
        "It is expected to fail for some low-f runs here, because Ks is not re-tuned in this script. (The rule was written for W2; on W4 read it as a guide.)")
    log()
    log("  STORMS")
    log("    " + f"{'run':<22}| {'fw vol':>6} {'fw pk%':>7} {'fw arr':>6} | {'f27 vol':>7} {'f27 pk%':>8} | {'span vol':>8} | {'storms_ok':>9} | note")
    for rid, rec, sc, tail in rows:
        if not sc:
            log("    " + f"{run_label(rec):<22}| ({rec.get('status', 'not run')})")
            continue
        note = "IDENTICAL to baseline" if rec.get("identical_to_baseline") else ""
        log("    " + f"{run_label(rec):<22}| {fm(fnum(sc, 'flood_week__vol_ratio'), 3):>6} {fm(fnum(sc, 'flood_week__peak_err_pct'), 1):>7} "
                     f"{fm(fnum(sc, 'flood_week__arrival_err_h')):>6} | {fm(fnum(sc, 'feb27_28__vol_ratio'), 3):>7} "
                     f"{fm(fnum(sc, 'feb27_28__peak_err_pct'), 1):>8} | {fm(fnum(sc, 'span__vol_ratio'), 3):>8} | {yn(storms_ok(sc)):>9} | {note}")
    first_tail = next((t for _, _, _, t in rows if t), None)
    log()
    if first_tail:
        log(f"  TAIL  (mean flow in 24 h blocks, m3/s. Blocks start {TAIL_OFFSETS_D[0]}, {TAIL_OFFSETS_D[1]}, {TAIL_OFFSETS_D[2]} days after the gauge flood peak "
            f"({first_tail.get('gauge_peak_time', '')}); band {args.band_lo:g}-{args.band_hi:g} m3/s; days-in-band counts the 24 h blocks "
            f"{TAIL_START[:10]} to {TAIL_END[:10]}, {first_tail['days_total']} blocks)")
        g = first_tail
        log("    " + f"{'run':<22}| " + " ".join(f"{'d' + str(d):>7}" for d in TAIL_OFFSETS_D) + f" | {'days in band':>12} | {'tail1 vol':>9} {'tail2 vol':>9} {'k ratio1':>8} | {'tail_ok':>7}")
        log("    " + f"{'gauge':<22}| " + " ".join(f"{fm(g[f'gauge_d{d}_m3s'], 3):>7}" for d in TAIL_OFFSETS_D)
            + f" | {str(g['gauge_days_in_band']) + '/' + str(g['days_total']):>12} |")
    else:
        log("  TAIL  (no scored run yet)")
    for rid, rec, sc, tail in rows:
        if not sc or not tail:
            continue
        log("    " + f"{run_label(rec):<22}| " + " ".join(f"{fm(tail[f'model_d{d}_m3s'], 3):>7}" for d in TAIL_OFFSETS_D)
            + f" | {str(tail['model_days_in_band']) + '/' + str(tail['days_total']):>12} | "
            + f"{fm(fnum(sc, 'tail_jan26_feb2__tail_vol_ratio')):>9} {fm(fnum(sc, 'tail_mar3_mar11__tail_vol_ratio')):>9} "
            + f"{fm(fnum(sc, 'tail_jan26_feb2__recession_rate_ratio')):>8} | {yn(tail['tail_ok']):>7}")
    log("    (tail vol = model / gauge volume on the Jan 26-Feb 2 and Mar 3-11 tails; k ratio1 = recession-rate ratio on the first tail, 1 = same recession, above 1 = model dries faster)")
    log()
    log("  PROCESS  (is any water arriving by a slow route? Water table = depth below the surface in mm, so a smaller number is a higher water table. 'decay' is exp(-f x depth) at the water table, median over the basin; assumes f per mm)")
    log("    " + f"{'run':<22}| {'Hsrf share':>10} | {'water table depth mm: first':>28} {'shallowest':>10} {'deepest':>8} | {'f median':>8} {'decay':>9} {'x baseline':>11} | snow max cm")
    for rid, rec, sc, tail in rows:
        if not sc:
            continue
        d = rec.get("decay") or {}
        dm = d.get("decay_median")
        rel = (dm / base_decay) if (dm is not None and base_decay) else float("nan")
        log("    " + f"{run_label(rec):<22}| {fm(fnum(sc, 'proc__hsrf_share'), 3):>10} | {fm(fnum(sc, 'proc__mdgw_first_mm'), 0):>28} "
                     f"{fm(fnum(sc, 'proc__mdgw_min_mm'), 0):>10} {fm(fnum(sc, 'proc__mdgw_max_mm'), 0):>8} | "
                     f"{fm(d.get('f_median'), 4):>8} {sci(dm if dm is not None else float('nan')):>9} {sci(rel):>11} | {fm(fnum(sc, 'proc__snow_max_cm'), 1)}")
    bad_decay = {(rec.get("decay") or {}).get("error") for rid, rec, sc, tail in rows if sc and (rec.get("decay") or {}).get("error")}
    if bad_decay:
        log("    (f / decay columns are blank for some runs: " + "; ".join(sorted(bad_decay)) + ". Run with --redo-decay to try again from the grids on disk.)")
    log()
    log(f"  DRY MONTHS AND WHOLE PERIOD  (the gauge is essentially dry Oct 1 - Dec 31: model flow there is false flow. Design target: at most "
        f"{0.5 * 2.5:.2f} M m3, half of the old baseline's 2.50)")
    log("    " + f"{'run':<22}| {'Oct-Dec model M m3':>18} {'gauge':>6} | {'days > 0.13 m3/s (model)':>25} {'(gauge)':>8} | {'max model m3/s':>14} | {'Oct-Mar vol':>11} {'span KGE':>9}")
    for rid, rec, sc, tail in rows:
        if not sc:
            continue
        log("    " + f"{run_label(rec):<22}| {fm(fnum(sc, 'octdec__vol_model_Mm3'), 2):>18} {fm(fnum(sc, 'octdec__vol_gauge_Mm3'), 3):>6} | "
                     f"{fm(fnum(sc, 'octdec__days_above_model'), 0):>25} {fm(fnum(sc, 'octdec__days_above_gauge'), 0):>8} | "
                     f"{fm(fnum(sc, 'octdec__max_model_m3s'), 3):>14} | {fm(fnum(sc, 'oct_mar__vol_ratio'), 3):>11} {fm(fnum(sc, 'span__kge2012'), 3):>9}")
    au_diagnostic(man, rows)
    if man.get("an_f"):
        a = man["an_f"]
        log()
        log(f"  Step 3 f: x{a['value']:g} ({a.get('source')}" + (f"; {a['rule']}; guard {a['guard']}" if a.get("rule") else "") + ")")
    failed = [rid for rid, rec in man["runs"].items() if rec.get("status") == "failed"]
    if failed:
        log()
        log("  FAILED runs (see the log): " + ", ".join(failed))

    out_rows = []
    for rid, rec, sc, tail in rows:
        d = rec.get("decay") or {}
        row = dict(step=rec.get("step"), run_id=rid, F_MULT=rec.get("f"), ANISOTROPY=rec.get("an"), ANISOTROPY_U=rec.get("au"), status=rec.get("status"),
                   storms_ok=(storms_ok(sc) if sc else ""), tail_ok=(tail["tail_ok"] if tail else ""))
        for c in SCORE_COLS:
            row[c] = fnum(sc, c) if sc else float("nan")
        if tail:
            row.update({k: v for k, v in tail.items() if k not in ("band_lo", "band_hi")})
        row.update(f_median=d.get("f_median"), gw_depth_median_mm=d.get("gw_depth_median_mm"), decay_median=d.get("decay_median"),
                   identical_to_baseline=bool(rec.get("identical_to_baseline")), identical_to_partner=bool(rec.get("identical_to_partner")), build_min=rec.get("build_min", ""),
                   tribs_min=rec.get("tribs_min", ""), output_mb=rec.get("output_mb", ""))
        out_rows.append(row)
    cols = []
    for r in out_rows:
        for c in r:
            if c not in cols:
                cols.append(c)
    with open(lab / SUMMARY_CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in out_rows:
            w.writerow({c: ("" if (isinstance(r.get(c), float) and not math.isfinite(r[c])) else r.get(c, "")) for c in cols})
    log()
    log(f"  Wrote {SUMMARY_CSV}. These tables do not decide anything: whether a lower f is worth its cost to the storms is the next step.")


# --------------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------------
def parse_args(argv=None):
    ap = argparse.ArgumentParser(description="Tail screen on W4: baseline, f screen, saturated-anisotropy screen with an Au diagnostic (see the top of this file).")
    ap.add_argument("--dry-run", action="store_true", help="check everything, list the plan, build nothing")
    ap.add_argument("--summary-only", action="store_true", help="print the tables from the runs already scored, then stop")
    ap.add_argument("--build-only", action="store_true", help="build and check the inputs but do not run tRIBS")
    ap.add_argument("--limit", type=int, default=None, help="make only the first N runs that still need work (finished runs do not count)")
    ap.add_argument("--steps", default=None, help="comma-separated steps to run: base,f,an (default all three, or none when --extra-runs is given; the baseline always runs first if it is not done)")
    ap.add_argument("--extra-runs", default=None, help="extra As runs as f:As pairs, Au kept at 1, for example 0.5874:300,0.15:300,0.02:200,0.02:500. "
                    "Each is the baseline with that F_MULT and that ANISOTROPY (As). Runs already finished are skipped. With this option and no --steps, only these runs are made.")
    ap.add_argument("--redo-decay", action="store_true", help="recompute the f / decay columns of the finished runs from the grids already on disk (no run is repeated), then print the tables")
    ap.add_argument("--an-f", type=float, default=None, help="F_MULT for the As screen and the Au diagnostic (default: chosen by the rule after the f screen)")
    ap.add_argument("--repick", action="store_true", help="choose the step 3 f again by the rule instead of keeping the earlier choice")
    ap.add_argument("--guard-lo", type=float, default=GUARD_LO, help="rule: lowest flood-week volume ratio allowed for the step 3 f")
    ap.add_argument("--guard-hi", type=float, default=GUARD_HI, help="rule: highest flood-week volume ratio allowed for the step 3 f")
    ap.add_argument("--band-lo", type=float, default=BAND_LO, help="slow-flow band, low end, m3/s")
    ap.add_argument("--band-hi", type=float, default=BAND_HI, help="slow-flow band, high end, m3/s")
    ap.add_argument("--redo", default=None, help="comma-separated run ids to run again even if done (you almost never want this)")
    ap.add_argument("--timeout-min", type=float, default=None, help="stop a tRIBS run after this many minutes (default: 3.5 times the expected time, about 2 hours for W4)")
    ap.add_argument("--cell-timeout", type=int, default=1800, help="seconds allowed for one notebook cell")
    ap.add_argument("--min-free-gb", type=float, default=4.0, help="stop if less than this much disk is free (W4 runs write more than W2 runs; default 4)")
    ap.add_argument("--heartbeat-min", type=float, default=HEARTBEAT_MIN,
                    help="print a progress line every this many minutes while tRIBS runs (default 5; keeps a Codespace awake)")
    ap.add_argument("--notebook", default=None, help="notebook to run (default: Make_Sycamore_Model.ipynb, else Make_Sycamore_Model_clean.ipynb)")
    ap.add_argument("--scorer", default=SCORER)
    ap.add_argument("--tribs", default=os.environ.get("TRIBS_BIN") or (TRIBS_DEFAULT if os.path.isfile(TRIBS_DEFAULT)
                                                                       else shutil.which("tRIBS") or TRIBS_DEFAULT))
    ap.add_argument("--skip-kernel-check", action="store_true", help="skip the pytRIBS import test in --dry-run")
    return ap.parse_args(argv)


def an_phase_jobs(fv: float) -> list:
    """The runs of step 3 at f = fv: the anisotropy runs, preceded by the f run at this f (the 'anisotropy 1' run). That f run is
    skipped when it is already finished, so nothing is repeated; it is only made if the f screen did not include it."""
    jobs = []
    if abs(fv - BASE["F_MULT"]) > 1e-12:
        jobs.append(f_job(fv))
    return jobs + an_jobs(fv)


def parse_extra_runs(text) -> list:
    """'0.15:300,0.02:500' -> As-only jobs (Au stays 1) at those f values. Raises ValueError with a plain message."""
    jobs, seen = [], set()
    for part in [p.strip() for p in (text or "").split(",") if p.strip()]:
        try:
            f_text, a_text = part.split(":")
            fv, a = float(f_text), float(a_text)
        except ValueError:
            raise ValueError(f"--extra-runs: '{part}' is not an f:As pair (for example 0.15:300)") from None
        if not (math.isfinite(fv) and math.isfinite(a)) or fv <= 0 or a <= 0:
            raise ValueError(f"--extra-runs: '{part}' needs f above 0 and As above 0")
        a = int(a) if a == int(a) else a
        job = an_job(fv, a)
        if job["run_id"] not in seen:
            seen.add(job["run_id"])
            jobs.append(job)
    return jobs


def planned_now(args, man: dict) -> list:
    """The jobs that can be listed right now: the baseline, the f screen, step 3 if its f is already known, and the extra runs."""
    steps = args.steps
    jobs = [make_base_job()]
    if "f" in steps:
        jobs += f_jobs()
    if "an" in steps:
        fv = args.an_f if args.an_f is not None else (man.get("an_f") or {}).get("value")
        if fv is not None:
            jobs += [j for j in an_phase_jobs(fv) if j["run_id"] not in {x["run_id"] for x in jobs}]
    have = {x["run_id"] for x in jobs}
    jobs += [j for j in getattr(args, "extra_jobs", []) if j["run_id"] not in have]
    return jobs


def main(argv=None) -> int:
    global _LOGFILE
    args = parse_args(argv)
    lab = Path(__file__).resolve().parent
    os.chdir(lab)
    (lab / LOG_DIR).mkdir(exist_ok=True)
    _LOGFILE = str(lab / LOG_DIR / "batch_log.txt")
    try:
        args.extra_jobs = parse_extra_runs(args.extra_runs)
    except ValueError as exc:
        print(str(exc))
        return 2
    if args.steps is None:
        args.steps = "" if args.extra_jobs else "base,f,an"
    args.steps = [x.strip() for x in args.steps.split(",") if x.strip()]
    bad = [x for x in args.steps if x not in ("base", "f", "an")]
    if bad or (not args.steps and not args.extra_jobs):
        print(f"--steps must be a comma-separated list of base, f, an (got {bad or 'nothing'}).")
        return 2
    if args.an_f is not None and args.an_f <= 0:
        print("--an-f must be above 0.")
        return 2
    man = load_manifest()

    if args.summary_only:
        print_summary(man, lab, args)
        return 0
    if args.redo_decay:
        n = redo_decay(man, lab)
        log(f"The f / decay numbers were recomputed for {n} finished run(s); nothing was rerun.")
        print_summary(man, lab, args)
        return 0

    args.notebook, nb_note = pick_notebook(lab, args.notebook)
    args.redo_ids = set(x.strip() for x in (args.redo or "").split(",") if x.strip())
    preview = planned_now(args, man)
    n_done = sum(1 for j in preview if is_finished(man, j, args))
    n_todo = len(preview) - n_done
    an_pending = "an" in args.steps and not any(j["step"] == "an" for j in preview)

    log(f"Working folder: {lab}")
    log(f"Notebook: {args.notebook}")
    log(f"Plan: steps {args.steps or 'none'}{' + ' + str(len(args.extra_jobs)) + ' extra run(s)' if args.extra_jobs else ''}, window {WIN['name']} ({WIN['hours']} h): {n_todo} run(s) to make"
        + (f", {n_done} already finished" if n_done else "")
        + (f" + {len(AN_VALUES) + 1} step 3 runs whose f is chosen after the f screen" if an_pending else "") + ".")
    for j in preview:
        what = "already done, skipped" if is_finished(man, j, args) else "NEW RUN"
        log(f"  {j['run_id']:44s} {describe(j):26s} {what}")
    if an_pending:
        log(f"  (step 3: {len(AN_VALUES)} runs at As {AN_VALUES} with Au = 1, plus one run with As = Au = {max(AN_VALUES)}; their f is the lowest "
            f"f-screen value whose flood-week volume ratio is between {args.guard_lo:g} and {args.guard_hi:g}, or the one you give with --an-f)")
    n_est = n_todo + (len(AN_VALUES) + 1 if an_pending else 0)
    log(f"  tRIBS time alone: about {n_est * WIN['hours'] * MIN_PER_1000H / 60000:.1f}-{n_est * WIN['hours'] * MIN_PER_1000H_SLOW / 60000:.1f} hours "
        f"({WIN['hours']} h runs take about {WIN['hours'] * MIN_PER_1000H / 1000:.0f}-{WIN['hours'] * MIN_PER_1000H_SLOW / 1000:.0f} min each); each run also needs a notebook build of about 1 min.")

    results = preflight(args, lab, preview, nb_note, man)
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

    man["script"] = "run_tail_screen_sycamore.py"
    man["baseline_id"] = BASELINE_ID
    man["last_started"] = now_iso()
    save_manifest(man)
    install_term_handler()

    t_batch, code = time.time(), 0
    budget = {"left": args.limit}

    def run_batch(todo: list) -> int:
        consecutive_fail, rc = 0, 0
        for j in todo:
            rid = j["run_id"]
            if not is_finished(man, j, args):
                if budget["left"] is not None:
                    if budget["left"] <= 0:
                        log(f"--limit reached; {rid} and the runs after it were not started.")
                        break
                    budget["left"] -= 1
            try:
                do_run(j, args, lab, man)
                consecutive_fail = 0
            except RunError as exc:
                man["runs"].setdefault(rid, {}).update(status="failed", error=str(exc), finished=now_iso())
                save_manifest(man)
                log(f"{rid}: FAILED. {exc}")
                consecutive_fail += 1
                rc = 1
                if consecutive_fail >= 2:
                    log("Two runs failed in a row, so the batch stops (the cause is probably the same).")
                    break
        return rc

    try:
        phase1 = [make_base_job()] + (f_jobs() if "f" in args.steps else [])
        code = run_batch(phase1)
        if "an" in args.steps and args.build_only:
            log("Step 3 is skipped with --build-only.")
        elif "an" in args.steps and code != 0:
            log("Step 3 was not started because a run in the earlier steps failed or was interrupted. Fix that, then run the same command again.")
        elif "an" in args.steps and budget["left"] is not None and budget["left"] <= 0:
            log("Step 3 was not started: --limit was used up. Run again without --limit (or with a larger one) to continue.")
        elif "an" in args.steps:
            info, why = pick_an_f(man, load_scores(lab), args)
            if info is None:
                log()
                log(f"Step 3 not started: {why}.")
            else:
                if man.get("an_f") != info:
                    man["an_f"] = info
                    save_manifest(man)
                log()
                log(f"Step 3: As screen (Au = 1) plus the As = Au diagnostic, at f x{info['value']:g} ({info.get('source')}"
                    + (f"; rule: {info['rule']}; guard band {info['guard']}; f-screen flood-week volume ratios {info['f_screen_flood_week_vol']}" if info.get("rule") else "")
                    + ")")
                code = run_batch(an_phase_jobs(info["value"])) or code
        if args.extra_jobs:
            if args.build_only:
                log("The extra runs are skipped with --build-only.")
            elif code != 0:
                log("The extra runs were not started because an earlier run failed or was interrupted. Fix that, then run the same command again.")
            elif budget["left"] is not None and budget["left"] <= 0:
                log("The extra runs were not started: --limit was used up. Run again without --limit (or with a larger one) to continue.")
            else:
                log()
                log("Extra As runs (Au = 1): " + ", ".join(f"f x{j['f']:g} As {j['an']:g}" for j in args.extra_jobs))
                code = run_batch(args.extra_jobs) or code
    except Abort as exc:
        log(f"BATCH STOPPED: {exc}")
        code = 1
    except KeyboardInterrupt:
        log("BATCH INTERRUPTED. Run the same command again to resume; finished runs are skipped.")
        for rid, rec in man["runs"].items():
            if rec.get("status") in ("building", "running"):
                rec["status"] = "interrupted"
        save_manifest(man)
        code = 130
    finally:
        save_manifest(man)

    log()
    log(f"Batch time: {fmt_min(time.time() - t_batch)}")
    log(f"{'run':46s} {'status':16s} {'build':>7s} {'tRIBS':>7s}")
    for rid, rec in man["runs"].items():
        log(f"{rid:46s} {rec.get('status', '?'):16s} {str(rec.get('build_min', '')):>7s} {str(rec.get('tribs_min', '')):>7s}")
    print_summary(man, lab, args)
    return code


if __name__ == "__main__":
    sys.exit(main())
