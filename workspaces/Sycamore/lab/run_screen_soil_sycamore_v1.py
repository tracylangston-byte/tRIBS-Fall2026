#!/usr/bin/env python3
"""
run_screen_soil_sycamore.py  --  Stage 2, piece 2: the soil-knob screen on W2.

Question it answers: when ONE soil knob is moved and everything else stays at the baseline, how much do the hydrograph
targets change on W2 (Jan 18 - Mar 12 2017, 1,272 h, water table at 0.82 x bedrock): the flood week, the second January
storm, the two February storms and the two recession tails? The seven soil knobs and the values screened (the design
doc's ranges table; the baseline value is never run again):

    ks    saturated conductivity multiplier (KS_MULT)               0.1  0.3  1  3  10      (baseline 0.2), then a second pass
    f     conductivity decay multiplier (F_MULT)                    0.25  0.5  2  4         (baseline 1)
    psib  air-entry pressure multiplier (PSIB_MULT)                 0.25  1.0  1.5          (baseline 0.5)
    m     pore-size index multiplier (M_MULT)                       0.8  1.2               (baseline 1)
    ts    saturated water content multiplier (THETAS_MULT)          1.05  1.10  1.15        (baseline 1)
    an    anisotropy ratio, As and Au together (ANISOTROPY)         50  100  300            (baseline 1)
    br    meters added to the SOLUS bedrock depth (BEDROCK_ADD_M)   0.5  2                (baseline 1)

The starting water table is not screened (Stage 1 found it inert): it stays at 0.82 x bedrock on W2.

NO RUN IS REPEATED. The W2 baseline and the three runs the Stage 1 proxy check already made (Ks x0.1, Ks x0.3 and
PsiB x1) are checked read-only and reused: their outlet files, checksums, md5 values and .in files are verified, and then
they are only scored. Nothing in results/ is ever overwritten: a planned run whose name already has files in results/
stops the start-up check. If a reused run fails its checks, it is left out and reported; it is not rerun.
That leaves 19 new runs: ts 3, an 3, m 2, ks 3, f 4, psib 2, br 2. After them the Ks second pass adds 3 more (see below).
The pore-size index m has its own grid (m_0-5cm; scgrid.gdf row PI), separate from porosity, so it is scaled by itself.

thetaS and porosity. pytRIBS writes the porosity grid as the SAME file as the thetaS grid, so capping thetaS at the
porosity grid (the design doc's wording) would cancel every upward step. The notebook therefore scales both together,
and this script checks both grids.

Ks second pass. The design doc says "a second pass refines where the response bends" and leaves the rule to judgment.
The rule used here: take the first-pass Ks values (the baseline included), find the pair of neighbors across which the
whole-window volume ratio (span__vol_ratio) changes most, and put 3 more values between them, evenly spaced on a log
scale. The interval, the values and the first-pass numbers are printed and saved in the manifest. Use --no-refine to skip it.

For each run to be made, the script:
    1. Copies the Settings cell of the Make_Sycamore_Model notebook, changes ONLY the one knob (plus the run names and
       the window), and runs the whole notebook headlessly (no browser) with nbclient. The notebook file is never changed.
    2. Checks the build: the new .in file must match the reference .in in every keyword except the run-specific ones
       (RUNTIME, STARTDATE, OUTFILENAME, GWATERFILE, BEDROCKFILE). Because the soil knobs do not appear in the .in file, it
       also checks the soil files themselves: every grid in data/model/soil against its SOLUS grid (the ratio must equal
       the multiplier), thetaS and porosity being the same file, and As and Au in soils.sdt.
    3. Runs tRIBS (one at a time, never two at once, in the foreground) and checks that it finished.
    4. Scores the run with score_run_sycamore.py (window W2) and flags a run whose outlet file is byte-identical to the
       baseline's (that would mean the knob changed nothing tRIBS reads).
    5. Writes everything to screen_soil_manifest.json, so a stopped batch can be resumed by running the same command.
The two new levers (ts and an) run first, so a problem with them shows up in the first half hour.

Run it from workspaces/Sycamore/lab/. Typical use:
    python -u run_screen_soil_sycamore.py --dry-run          # checks only, builds nothing
    python -u run_screen_soil_sycamore.py                    # everything; run it again to resume after any stop
    python -u run_screen_soil_sycamore.py --summary-only     # just print the tables from what is finished
    python -u run_screen_soil_sycamore.py --knobs ts,an      # only some knobs (ks,f,psib,m,ts,an,br)
    python -u run_screen_soil_sycamore.py --limit 1 --build-only   # build the next run and check it, no tRIBS

Files it makes (all in lab/): screen_soil_manifest.json, screen_soil.lock (only while running), screen_soil_logs/,
score_out_screen_soil/scores_screen_soil.csv, screen_soil_summary.csv, and for each new run <run>.in,
simulation_<run>.log and results/<run>*.

Do NOT run the notebook (or any cell that rebuilds data/model/) while this runs: every build rewrites lab/data/model/,
and all the runs share it. Keep the terminal open: the progress line printed every 5 minutes also stops the Codespace
from timing out for being idle. If the connection drops, run the same command again: finished runs are skipped.
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
# What the screen does
# --------------------------------------------------------------------------------------------
WIN = dict(name="W2", start="01/18/2017/00/00", hours=1272, label="jan_to_mar", gw=0.82)   # label = the RUN_CHOICE text printed

# Everything that is not being screened (the Oct 21 baseline, with the W2 cold-start water table).
BASE = {
    "GW_FRACTION": WIN["gw"], "BEDROCK_ADD_M": 1.0, "KS_MULT": 0.2, "PSIB_MULT": 0.5, "F_MULT": 1,
    "THETAS_MULT": 1, "ANISOTROPY": 1, "M_MULT": 1,
    "KINEMVELCOEF": 6, "FLOWEXP": 0.3, "CHANNELROUGHNESS": 0.04, "CHANNELWIDTHCOEFF": 2.33,
    "OPTSNOW": 1, "OPTPERCOLATION": 0, "CHANNELCONDUCTIVITY": 5,
    "REBUILD_SOIL": False,
}
# key = what you type after --knobs; tag = what goes in the run name; kw = the Settings name. Listed in the order the tables print.
KNOBS = [
    dict(key="ks",   kw="KS_MULT",       tag="ks",   base=0.2, values=[0.1, 0.2, 0.3, 1, 3, 10],
         label="saturated conductivity multiplier"),
    dict(key="f",    kw="F_MULT",        tag="f",    base=1,   values=[0.25, 0.5, 1, 2, 4],
         label="conductivity decay (f) multiplier"),
    dict(key="psib", kw="PSIB_MULT",     tag="psib", base=0.5, values=[0.25, 0.5, 1.0, 1.5],
         label="air-entry pressure (PsiB) multiplier"),
    dict(key="m",    kw="M_MULT",        tag="m",    base=1,   values=[0.8, 1.2],
         label="pore-size index (m) multiplier"),
    dict(key="ts",   kw="THETAS_MULT",   tag="ts",   base=1,   values=[1, 1.05, 1.10, 1.15],
    dict(key="an",   kw="ANISOTROPY",    tag="an",   base=1,   values=[1, 50, 100, 300],
         label="anisotropy ratio (As and Au together)"),
    dict(key="br",   kw="BEDROCK_ADD_M", tag="br",   base=1.0, values=[0.5, 1.0, 2.0],
         label="meters added to the SOLUS bedrock depth"),
]
RUN_ORDER = ["ts", "an", "m", "ks", "f", "psib", "br"]  # the order runs are made in (the new levers first)
SMOKE_FIRST = ["ts", "an", "m"]                              # the largest step of each new lever runs before anything else
ROUTING_KWS = ["KINEMVELCOEF", "FLOWEXP", "CHANNELROUGHNESS", "CHANNELWIDTHCOEFF"]   # never screened here; must stay at baseline

BASELINE_ID = "W2_ks0p2_cv6_psib0p5_nocl"              # the existing W2 baseline run (from the proxy check)
# Runs the Stage 1 proxy check already made with exactly one screened knob changed. They are reused, never rerun.
PROXY_TWINS = {
    ("ks", 0.1): "W2_ks0p1_cv6_psib0p5_nocl",
    ("ks", 0.3): "W2_ks0p3_cv6_psib0p5_nocl",
    ("psib", 1.0): "W2_ks0p2_cv6_psib1p0_nocl",
}
PROXY_MANIFEST = "proxy_manifest.json"
NOTEBOOK_CANDIDATES = ["Make_Sycamore_Model.ipynb", "Make_Sycamore_Model_clean.ipynb"]
SCORER = "score_run_sycamore.py"
TRIBS_DEFAULT = "/opt/tRIBS/build/tRIBS"
INIT_DATA = "../init_data"           # SOLUS_DIR in the notebook
MODEL_SOIL = "data/model/soil"
LOG_DIR = "screen_soil_logs"
SCORE_DIR = "score_out_screen_soil"
SCORE_CSV = "scores_screen_soil.csv"
SUMMARY_CSV = "screen_soil_summary.csv"
MANIFEST = "screen_soil_manifest.json"
LOCK = "screen_soil.lock"
MIN_PER_1000H = 7.8          # measured tRIBS speed on this machine: minutes per 1,000 model hours (W2 took 10.0-12.9 min)
HEARTBEAT_MIN = 5.0
ROWS_PER_MODEL_HOUR = 16     # the outlet file gets 16 rows per model hour
F_FLOOR = 1e-4               # the notebook keeps f at or above this
THETAS_CAP = 0.99            # the notebook keeps thetaS and porosity below this
REFINE_METRIC = "span__vol_ratio"   # the Ks second pass looks at this score column
REFINE_N = 3                 # runs in the Ks second pass

# Score columns: (short name, scorer column). The first group is shown in the tables.
EVENTS = [("fw", "flood_week", "FLOOD WEEK Jan 18-25"), ("s2", "
         label="saturated water content (thetaS) multiplier; porosity follows"),jan_second_storm", "2ND JAN STORM"),
          ("f18", "feb18_19", "FEB 18-19"), ("f27", "feb27_28", "FEB 27-28")]
METRICS = []
for _short, _col, _title in EVENTS:
    METRICS += [(f"{_short}_vol", f"{_col}__vol_ratio"), (f"{_short}_pk", f"{_col}__peak_err_pct"),
                (f"{_short}_arr", f"{_col}__arrival_err_h")]
METRICS += [("tail1", "tail_jan26_feb2__tail_vol_ratio"), ("tail2", "tail_mar3_mar11__tail_vol_ratio"),
            ("bal", "balance__feb27_28_over_feb18_19__model"), ("span", "span__vol_ratio"),
            ("kge_fw", "flood_week__kge2012")]
# Saved in the CSV and the manifest but not in the printed tables.
EXTRA_COLS = ["flood_week__peak_time_diff_h", "jan_second_storm__kge2012", "feb18_19__kge2012", "feb27_28__kge2012",
              "jan_flood_tail__vol_ratio", "feb18_19_tail__vol_ratio", "feb27_28_tail__vol_ratio",
              "tail_jan26_feb2__recession_rate_ratio", "tail_mar3_mar11__recession_rate_ratio",
              "balance__feb27_28_tail_over_feb18_19_tail__model", "balance__feb27_28_over_feb18_19__gauge",
              "proc__mdgw_first_mm", "proc__mdgw_min_mm", "proc__mdgw_max_mm", "proc__channel_loss_Mm3",
              "proc__hsrf_share", "proc__snow_max_cm", "span__vol_model_Mm3_matched", "span__vol_gauge_Mm3"]
FLAG_COLS = ["flood_week__arrival_flag", "jan_second_storm__arrival_flag", "feb18_19__arrival_flag", "feb27_28__arrival_flag"]
SHOW_COLS = [c for _, c in METRICS] + EXTRA_COLS

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
    """0.18 -> '0p18', 1.5 -> '1p5', 15 -> '15', 0.10 -> '0p1', 0.0 -> '0'."""
    return f"{x:g}".replace(".", "p")


def run_id_for(knob: dict | None, value=None) -> str:
    """W2_ks0p2_cv6_psib0p5_nocl is the baseline. Ks and PsiB go in their usual slots; the other knobs go after psib,
    before the loss setting. The runs the proxy check already made keep their old names (they are reused, not rerun)."""
    if knob is None:
        return BASELINE_ID
    twin = PROXY_TWINS.get((knob["key"], value))
    if twin:
        return twin
    ks, psib, extra = vtag(BASE["KS_MULT"]), vtag(BASE["PSIB_MULT"]), ""
    if knob["key"] == "ks":
        ks = vtag(value)
    elif knob["key"] == "psib":
        psib = vtag(value)
    else:
        extra = f"_{knob['tag']}{vtag(value)}"
    return f"{WIN['name']}_ks{ks}_cv6_psib{psib}{extra}_nocl"


def make_job(knob: dict | None, value=None) -> dict:
    if knob is None:
        return dict(knob=None, kw=None, value=None, run_id=BASELINE_ID, kind="baseline")
    kind = "reuse" if (knob["key"], value) in PROXY_TWINS else "run"
    return dict(knob=knob["key"], kw=knob["kw"], value=value, run_id=run_id_for(knob, value), kind=kind)


def knob_by_key(key: str) -> dict:
    return next(k for k in KNOBS if k["key"] == key)


def screened_values(knob: dict, man: dict | None = None) -> list:
    """The values that need a run: every screened value except the baseline one (plus the Ks second pass once chosen)."""
    vals = list(knob["values"])
    if knob["key"] == "ks" and man and man.get("ks_refine"):
        vals += man["ks_refine"]["values"]
    return sorted(v for v in set(vals) if v != knob["base"])


def settings_for(job: dict) -> dict:
    """Every NAME = value line of the Settings cell that this run sets."""
    s = dict(BASE)
    s.update(RUN_CHOICE=WIN["label"], INPUT_FILE=f"{job['run_id']}.in", RESULTS_NAME=job["run_id"])
    if job["kw"]:
        s[job["kw"]] = job["value"]
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
    """Compare a settings record from the manifest with the current one. A record written before a knob existed
    (for example M_MULT) is read as having that knob at its baseline value."""
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
    return max(30.0, 3.5 * WIN["hours"] * MIN_PER_1000H / 1000.0)


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
                                 display_name="screen", language="python")
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


def check_in_file(new_in: Path, ref_in: Path, job: dict) -> None:
    """The new .in must equal the reference .in except for the run-specific lines (and the bedrock file name for a bedrock
    run), and the routing keywords must have exactly their baseline values. The soil knobs do not appear in the .in file;
    check_grids covers them."""
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
    notes = [check_one_grid(lab, var, model_rows[var], src_rows[var], mults[var], t_build_start) for var in mults]
    # the soil table: As and Au
    sdt = soil / "soils.sdt"
    if not sdt.is_file() or sdt.stat().st_mtime < t_build_start - 2:
        raise RunError("data/model/soil/soils.sdt is missing or older than this build.")
    table = [ln.split(",") for ln in sdt.read_text().splitlines()[1:] if ln.strip()]
    if not table or any(len(r) != 12 for r in table):
        raise RunError("soils.sdt does not have 12 comma-separated columns on every row.")
    for r in table:
        if float(r[7]) != float(s["ANISOTROPY"]) or float(r[8]) != float(s["ANISOTROPY"]):
            raise RunError(f"soils.sdt row {r[0]} has As = {r[7]}, Au = {r[8]}; ANISOTROPY is {s['ANISOTROPY']:g}.")
    return (f"Ks x{s['KS_MULT']:g}, PsiB x{s['PSIB_MULT']:g}, f x{s['F_MULT']:g}, m x{s['M_MULT']:g}, thetaS x{s['THETAS_MULT']:g} (porosity the same file), "
            f"As = Au = {s['ANISOTROPY']:g}; bedrock +{s['BEDROCK_ADD_M']:g} m (median {np.median(bed[ok]):.2f} m); "
            f"water table = {s['GW_FRACTION']:g} x bedrock ({unit}), median depth {np.median(gw[ok]):.0f} {unit}")


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


def run_scorer(rid: str, scorer: str, score_dir: Path, log_dir: Path, write_checksum: bool = True) -> dict:
    """Score one W2 run; returns its whole row of the scorer's CSV (all values as text)."""
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
    out = {}
    for c in SHOW_COLS:
        v = fnum(row, c)
        if math.isfinite(v):
            out[c] = round(v, 4)
    for c in FLAG_COLS:
        if (row or {}).get(c):
            out[c] = row[c]
    return out


def one_line(scores: dict) -> str:
    names = [("flood_week__vol_ratio", "flood vol"), ("flood_week__peak_err_pct", "peak %"),
             ("flood_week__arrival_err_h", "arrival h"), ("tail_jan26_feb2__tail_vol_ratio", "tail1"),
             ("tail_mar3_mar11__tail_vol_ratio", "tail2"), ("span__vol_ratio", "W2 vol")]
    return ", ".join(f"{label} {scores[c]}" for c, label in names if c in scores)


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
# The reference run and the other existing W2 runs (checked read-only, never rerun)
# --------------------------------------------------------------------------------------------
def verify_baseline_in(path: Path) -> None:
    """A reused run's .in file must carry the baseline routing values, loss off, snow on, and the W2 window."""
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
        raise RunError(f"{path.name} is not a baseline-routing W2 run: " + "; ".join(bad))


def reuse_check(lab: Path, rid: str, ks_mult: float, psib_mult: float) -> tuple:
    """Read-only check of a W2 run the proxy check already made. Returns (problems, notes). No problems = it can be reused."""
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
    pm = Path(PROXY_MANIFEST)
    if pm.is_file():
        try:
            prec = json.loads(pm.read_text())["runs"].get(rid)
        except Exception:                                                  # an unreadable manifest is not a reason to reject
            prec = None
            notes.append(f"{PROXY_MANIFEST} could not be read")
        if prec:
            if prec.get("status") != "done":
                problems.append(f"{PROXY_MANIFEST} says this run's status is {prec.get('status')!r}, not 'done'")
            md5 = (prec.get("outlet") or {}).get("md5")
            if md5 and md5 != outlet["md5"]:
                problems.append(f"the outlet file's md5 ({outlet['md5']}) differs from the one recorded in {PROXY_MANIFEST} ({md5})")
            elif md5:
                notes.append(f"md5 matches {PROXY_MANIFEST}")
            want = {"GW_FRACTION": WIN["gw"], "BEDROCK_ADD_M": BASE["BEDROCK_ADD_M"], "KS_MULT": ks_mult,
                    "PSIB_MULT": psib_mult, "START_DATE": WIN["start"], "RUN_HOURS": WIN["hours"]}
            if prec.get("settings") and prec["settings"] != want:
                problems.append(f"{PROXY_MANIFEST} records settings {prec['settings']}, expected {want}")
        elif "could not be read" not in " ".join(notes):
            notes.append(f"{PROXY_MANIFEST} has no record of this run")
    else:
        notes.append(f"{PROXY_MANIFEST} not found, so the run's origin is checked only by its file, checksum and .in file")
    inp = lab / f"{rid}.in"
    if not inp.is_file():
        problems.append(f"{inp.name} not found in lab/")
    else:
        try:
            verify_baseline_in(inp)
            notes.append(f"{inp.name} has the baseline routing values")
        except RunError as exc:
            problems.append(str(exc))
    return problems, notes


def do_reuse(job: dict, args, lab: Path, man: dict, ref_in: Path | None) -> bool:
    """Check an existing run read-only and score it into the screen's own score CSV. It is never rerun.
    Returns True if it can be used."""
    rid, s = job["run_id"], settings_for(job)
    log_dir, score_dir = lab / LOG_DIR, lab / SCORE_DIR
    log_dir.mkdir(exist_ok=True)
    score_dir.mkdir(exist_ok=True)
    rec = man["runs"].setdefault(rid, {})
    cur = settings_record(job)
    if rec.get("status") in ("reused",) and same_record(rec.get("settings"), cur) and rid not in args.redo_ids and rec.get("scores"):
        log(f"{rid}: already checked and scored as a reused run, skipping")
        return True
    rec.clear()
    rec.update(kind="reuse", knob=job["knob"], kw=job["kw"], value=job["value"], status="checking", started=now_iso())
    problems, notes = reuse_check(lab, rid, s["KS_MULT"], s["PSIB_MULT"])
    if not problems and ref_in is not None:
        try:
            check_in_file(lab / f"{rid}.in", ref_in, job)
            notes.append("its .in file matches the baseline's except the run-specific lines")
        except RunError as exc:
            problems.append(str(exc))
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
        log(f"{rid}: can NOT be reused, and will NOT be rerun:")
        for p in problems:
            log(f"    {p}")
        return False
    rec.update(status="reused", finished=now_iso(), settings=cur)
    save_manifest(man)
    for n in notes:
        log(f"{rid}: {n}")
    log(f"{rid}: reused, not rerun. {one_line(rec['scores'])}")
    return True


# --------------------------------------------------------------------------------------------
# Preflight
# --------------------------------------------------------------------------------------------
def pick_notebook(lab: Path, given: str | None) -> tuple:
    """(notebook name, note). With no --notebook, Make_Sycamore_Model.ipynb if it exists (that is what built the baseline),
    otherwise Make_Sycamore_Model_clean.ipynb."""
    if given:
        return given, ""
    present = [n for n in NOTEBOOK_CANDIDATES if (lab / n).is_file()]
    if not present:
        return NOTEBOOK_CANDIDATES[0], ""
    note = ""
    if len(present) > 1:
        note = (f"both {present[0]} and {present[1]} are in lab/; using {present[0]} (the one the proxy runner uses). "
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
        add("FAIL", f"the W2 baseline {BASELINE_ID} cannot be the reference run: " + "; ".join(base_problems)
            + " (nothing was made, and nothing is rerun automatically)")
    else:
        add("OK", f"reference run {BASELINE_ID} passes its checks (reused, not rerun)")
    for j in jobs:
        if j["kind"] == "reuse":
            probs, _ = reuse_check(lab, j["run_id"], settings_for(j)["KS_MULT"], settings_for(j)["PSIB_MULT"])
            add("WARN" if probs else "OK",
                f"existing run {j['run_id']} ({describe(j)}) " + ("can NOT be reused (it will be left out, not rerun): " + "; ".join(probs)
                                                               if probs else "passes its checks (reused, not rerun)"))
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
            every = [make_job(k, v) for k in KNOBS for v in screened_values(k, man) if (k["key"], v) not in PROXY_TWINS]
            for j in every:
                patched_source(src, j)
            add("OK", f"new notebook: Settings cell (cell {idx}) has M_MULT, THETAS_MULT and ANISOTROPY, the cap and theta*_s cells are there, "
                      f"and it patches for all {len(every)} planned runs")
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
    # nothing may be overwritten: a new run's name must not already have files in results/
    results_dir = lab / "results"
    for j in jobs:
        if j["kind"] != "run":
            continue
        stale = list(results_dir.glob(j["run_id"] + "*")) if results_dir.is_dir() else []
        if stale and j["run_id"] not in man["runs"]:
            add("FAIL", f"results/ already has {len(stale)} file(s) named {j['run_id']}*; this script never overwrites existing results")
    # disk: the baseline's files give the size of one run
    base_bytes = sum(p.stat().st_size for p in results_dir.glob(BASELINE_ID + "*") if p.is_file()) if results_dir.is_dir() else 0
    n_new = sum(1 for j in jobs if j["kind"] == "run" and (man["runs"].get(j["run_id"]) or {}).get("status") != "done")
    if "ks" in args.knob_keys and not args.no_refine and not man.get("ks_refine"):
        n_new += REFINE_N
    need = base_bytes * n_new / 1e9
    add("FAIL" if (base_bytes and need > free - args.min_free_gb) else "OK",
        f"disk: one run's files are about {base_bytes / 1e6:.0f} MB (the baseline's); {n_new} runs need about {need:.1f} GB; "
        f"{free:.1f} GB free")
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


def describe(job: dict) -> str:
    return "baseline (no change)" if job["kw"] is None else f"{job['knob']} = {job['value']:g}"


def output_mb(lab: Path, rid: str) -> float:
    return round(sum(p.stat().st_size for p in (lab / "results").glob(rid + "*") if p.is_file()) / 1e6, 1)


def do_run(job: dict, args, lab: Path, man: dict, ref_in: Path) -> None:
    if job["kind"] == "reuse":
        do_reuse(job, args, lab, man, ref_in)
        return
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
    rec.update(status="building", started=now_iso(), kind=job["kind"], knob=job["knob"], kw=job["kw"], value=job["value"],
               settings=cur)
    save_manifest(man)
    remove_stale(lab, rid)
    log(f"=== {rid}: {describe(job)} on {WIN['name']} ({WIN['hours']} h from {WIN['start']}, water table {WIN['gw']:g} x bedrock) ===")

    t_build = time.time()
    text, warns = build_run(job, args, lab, log_dir)
    check_in_file(lab / s["INPUT_FILE"], ref_in, job)
    grid_note = check_grids(lab, job, t_build)
    rec.update(build_min=round((time.time() - t_build) / 60, 2), warnings=warns, grid_check=grid_note, status="built")
    save_manifest(man)
    log(f"{rid}: build OK in {fmt_min(time.time() - t_build)}; the .in file matches the reference except the allowed lines")
    log(f"{rid}: {grid_note}")
    for wmsg in warns:
        log(f"{rid}: notebook said: {wmsg}")
    if args.build_only:
        log(f"{rid}: --build-only, not running tRIBS")
        return

    install_term_handler()
    log(f"{rid}: starting tRIBS ({WIN['hours']} model hours, expect about {WIN['hours'] * MIN_PER_1000H / 1000:.0f}-13 min)")
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
        log(f"{rid}: WARNING: the outlet file is byte-identical to the baseline's. tRIBS ran, but this change had no effect on "
            f"the outlet flow. Either the knob has no pathway to the outlet here, or tRIBS does not read the setting.")

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


# --------------------------------------------------------------------------------------------
# The table and the summary CSV
# --------------------------------------------------------------------------------------------


def pick_ks_refinement(scores: dict, man: dict) -> tuple:
    """The Ks second pass. Find the neighbors (first-pass Ks values, the baseline included) across which REFINE_METRIC changes
    most, and return REFINE_N values evenly spaced between them on a log scale. Returns (info or None, reason)."""
    knob = knob_by_key("ks")
    pts = []
    for v in sorted(set(knob["values"])):
        rid = BASELINE_ID if v == knob["base"] else run_id_for(knob, v)
        sc = scores.get(rid)
        val = fnum(sc, REFINE_METRIC) if sc else float("nan")
        if not math.isfinite(val):
            return None, f"no {REFINE_METRIC} for Ks x{v:g} ({rid})"
        pts.append((v, val))
    steps = [(abs(pts[i + 1][1] - pts[i][1]), i) for i in range(len(pts) - 1)]
    d, i = max(steps, key=lambda t: (t[0], -t[1]))
    a, b = pts[i][0], pts[i + 1][0]
    vals = []
    for sig in (2, 3, 4):
        vals = [float(f"{a * (b / a) ** (j / (REFINE_N + 1)):.{sig}g}") for j in range(1, REFINE_N + 1)]
        if len(set(vals)) == REFINE_N and all(a < x < b for x in vals):
            break
    return dict(metric=REFINE_METRIC, interval=[a, b], change=round(d, 4), pass1=[[x, round(y, 4)] for x, y in pts],
                values=vals), ""


def load_scores(lab: Path) -> dict:
    p = lab / SCORE_DIR / SCORE_CSV
    out = {}
    if p.is_file():
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh):
                out[r["run_id"]] = r
    return out


def fm(x, nd=2, sign=False) -> str:
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "-"
    if abs(x) < 0.5 * 10 ** -nd:                       # avoid printing -0.00
        x = 0.0
    return f"{x:+.{nd}f}" if sign else f"{x:.{nd}f}"


def print_summary(man: dict, lab: Path, knob_keys: list) -> None:
    scores = load_scores(lab)
    base = scores.get(BASELINE_ID)
    out_rows = []
    log()
    log("SOIL SCREEN ON W2 (Jan 18 - Mar 12 2017). Everything not screened is at the baseline.")
    log("  vol = model volume / gauge volume.  pk = model peak vs gauge peak, %.  arr = first-arrival error, hours "
        "(model minus gauge; negative = early).")
    log("  tail1, tail2 = model / gauge volume on the Jan 26-Feb 2 and Mar 3-11 recession tails (proposed tolerance: 0.7 or more).")
    log("  bal = runoff ratio of Feb 27-28 over Feb 18-19 (the gauge's is 3.4).  span = whole-window volume ratio.  "
        "kge = flood-week KGE_2012.")
    log("  chg = the run minus the baseline row (volume ratio points, peak percentage points, hours).")
    log("  Proposed tolerances: volume within 10%, peak within 15%, arrival within 1 h.")
    if base is None:
        log(f"  (the baseline run {BASELINE_ID} has no score yet, so the chg columns are blank)")
    for k in KNOBS:
        if k["key"] not in knob_keys:
            continue
        values = sorted(set(k["values"]) | set(screened_values(k, man)))
        log()
        log(f"  {k['key']:<5}{k['label']} ({k['kw']}), baseline {k['base']:g}")
        head = "    " + f"{'value':>7} | " + " | ".join(f"{t:<21}" for _, _, t in EVENTS)
        log(head)
        log("    " + f"{'':>7} | " + " | ".join(f"{'vol':>6} {'pk %':>7} {'arr h':>6}" for _ in EVENTS))
        rows_info = []
        for v in values:
            is_base = (v == k["base"])
            rid = BASELINE_ID if is_base else run_id_for(k, v)
            sc = scores.get(rid)
            rec = man["runs"].get(rid) or {}
            status = rec.get("status", "not run")
            vals = {n: fnum(sc, c) for n, c in METRICS} if sc else None
            if vals is None:
                log("    " + f"{v:>7g} | {'(' + ('baseline: ' if is_base else '') + status + ')':^70}")
            else:
                groups = [f"{fm(vals[s + '_vol'], 3):>6} {fm(vals[s + '_pk'], 1):>7} {fm(vals[s + '_arr']):>6}" for s, _, _ in EVENTS]
                log("    " + f"{v:>7g} | " + " | ".join(groups))
            rows_info.append((v, is_base, rid, sc, vals, rec, status))
        log("    " + f"{'value':>7} | {'tail1':>6} {'tail2':>6} {'bal':>6} {'span':>6} {'kge':>6} | {'chg fw vol':>10} {'chg fw pk':>9} "
                     f"{'chg fw arr':>10} {'chg span':>8} | note")
        for v, is_base, rid, sc, vals, rec, status in rows_info:
            row = dict(knob=k["key"], kw=k["kw"], value=v, run_id=rid, is_baseline=is_base, status=status, kind=rec.get("kind", ""))
            if vals is None:
                log("    " + f"{v:>7g} | {'(' + ('baseline: ' if is_base else '') + status + ')':^60}")
                out_rows.append(row)
                continue
            d = {n: (vals[n] - fnum(base, c) if (base is not None and not is_base) else float("nan")) for n, c in METRICS}
            notes = []
            if is_base:
                notes.append("baseline (reused)")
            elif rec.get("kind") == "reuse":
                notes.append("reused from the proxy check")
            if rec.get("identical_to_baseline"):
                notes.append("IDENTICAL to baseline")
            bad = [lbl for lbl, c in zip(("fw", "s2", "f18", "f27"), FLAG_COLS) if sc.get(c) not in ("ok", "", None)]
            if bad:
                notes.append("arrival n/a: " + ",".join(bad))
            log("    " + f"{v:>7g} | {fm(vals['tail1']):>6} {fm(vals['tail2']):>6} {fm(vals['bal']):>6} {fm(vals['span'], 3):>6} "
                         f"{fm(vals['kge_fw'], 3):>6} | {fm(d['fw_vol'], 3, True):>10} {fm(d['fw_pk'], 1, True):>9} "
                         f"{fm(d['fw_arr'], 2, True):>10} {fm(d['span'], 3, True):>8} | {'; '.join(notes)}")
            for n, c in METRICS:
                row[c] = vals[n]
                row["change_" + c] = d[n]
            for c in EXTRA_COLS:
                row[c] = fnum(sc, c)
            for c in FLAG_COLS:
                row[c] = sc.get(c, "")
            row["identical_to_baseline"] = bool(rec.get("identical_to_baseline"))
            row["build_min"], row["tribs_min"], row["output_mb"] = rec.get("build_min", ""), rec.get("tribs_min", ""), rec.get("output_mb", "")
            out_rows.append(row)
    ref = man.get("ks_refine")
    if ref and "ks" in knob_keys:
        log()
        log(f"  Ks second pass: across Ks x{ref['interval'][0]:g} to x{ref['interval'][1]:g} the {ref['metric']} changed most "
            f"({ref['change']:.3f}); extra values {ref['values']}.")
    rejected = [rid for rid, rec in man["runs"].items() if rec.get("status") == "reuse_rejected"]
    if rejected:
        log()
        log("  NOT reused (failed their checks) and NOT rerun: " + ", ".join(rejected))
    if out_rows:
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
        log(f"  Wrote {SUMMARY_CSV}. These tables do not label any knob active, inert or trade-off: that is the analysis step.")


# --------------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------------
def parse_args(argv=None):
    ap = argparse.ArgumentParser(description="Stage 2 piece 2: the soil-knob screen on W2 (see the top of this file).")
    ap.add_argument("--dry-run", action="store_true", help="check everything, list the plan, build nothing")
    ap.add_argument("--summary-only", action="store_true", help="print the tables from the runs already scored, then stop")
    ap.add_argument("--build-only", action="store_true", help="build and check the inputs but do not run tRIBS")
    ap.add_argument("--limit", type=int, default=None, help="make only the first N runs that still need work (finished and reused runs do not count)")
    ap.add_argument("--knobs", default=None, help="comma-separated knobs: ks,f,psib,m,ts,an,br (default: all seven)")
    ap.add_argument("--no-refine", action="store_true", help="skip the Ks second pass")
    ap.add_argument("--redo", default=None, help="comma-separated run ids to run again even if done (you almost never want this)")
    ap.add_argument("--timeout-min", type=float, default=None, help="stop a tRIBS run after this many minutes (default: about 35)")
    ap.add_argument("--cell-timeout", type=int, default=1800, help="seconds allowed for one notebook cell")
    ap.add_argument("--min-free-gb", type=float, default=2.0)
    ap.add_argument("--heartbeat-min", type=float, default=HEARTBEAT_MIN,
                    help="print a progress line every this many minutes while tRIBS runs (default 5; keeps a Codespace awake)")
    ap.add_argument("--notebook", default=None, help="notebook to run (default: Make_Sycamore_Model.ipynb, else Make_Sycamore_Model_clean.ipynb)")
    ap.add_argument("--scorer", default=SCORER)
    ap.add_argument("--tribs", default=os.environ.get("TRIBS_BIN") or (TRIBS_DEFAULT if os.path.isfile(TRIBS_DEFAULT)
                                                                       else shutil.which("tRIBS") or TRIBS_DEFAULT))
    ap.add_argument("--skip-kernel-check", action="store_true", help="skip the pytRIBS import test in --dry-run")
    return ap.parse_args(argv)


def pick(arg, allowed, label):
    if not arg:
        return list(allowed)
    want = [x.strip() for x in arg.split(",") if x.strip()]
    unknown = [x for x in want if x not in allowed]
    if unknown:
        raise SystemExit(f"Unknown {label} {unknown}. Choose from {list(allowed)}")
    return [a for a in allowed if a in want]


def planned(knob_keys: list, man: dict) -> list:
    """Reused runs first (they cost nothing), then the largest step of each new lever, then everything else in RUN_ORDER."""
    reuse, runs = [], []
    for key in RUN_ORDER:
        if key not in knob_keys:
            continue
        knob = knob_by_key(key)
        jobs = [make_job(knob, v) for v in screened_values(knob, man)]
        reuse += [j for j in jobs if j["kind"] == "reuse"]
        runs += [j for j in jobs if j["kind"] == "run"]
    first = []
    for key in SMOKE_FIRST:
        cand = [j for j in runs if j["knob"] == key]
        if cand:
            first.append(max(cand, key=lambda j: j["value"]))
    return reuse + first + [j for j in runs if j not in first]


def main(argv=None) -> int:
    global _LOGFILE
    args = parse_args(argv)
    lab = Path(__file__).resolve().parent
    os.chdir(lab)
    (lab / LOG_DIR).mkdir(exist_ok=True)
    _LOGFILE = str(lab / LOG_DIR / "batch_log.txt")
    args.knob_keys = knob_keys = pick(args.knobs, [k["key"] for k in KNOBS], "knob")
    man = load_manifest()

    if args.summary_only:
        print_summary(man, lab, knob_keys)
        return 0

    args.notebook, nb_note = pick_notebook(lab, args.notebook)
    args.redo_ids = set(x.strip() for x in (args.redo or "").split(",") if x.strip())

    base_problems, base_notes = reuse_check(lab, BASELINE_ID, BASE["KS_MULT"], BASE["PSIB_MULT"])
    preview = planned(knob_keys, man)
    n_run = sum(1 for j in preview if j["kind"] == "run")
    n_reuse = sum(1 for j in preview if j["kind"] == "reuse")
    refine_note = ("" if ("ks" not in knob_keys or args.no_refine) else
                   (f" + the Ks second pass already chosen ({len(man['ks_refine']['values'])} runs, in the list)" if man.get("ks_refine")
                    else f" + {REFINE_N} Ks second-pass runs chosen after the first pass"))
    log(f"Working folder: {lab}")
    log(f"Notebook: {args.notebook}")
    log(f"Plan: knobs {knob_keys}, window {WIN['name']} ({WIN['hours']} h): {n_run} run(s) to make{refine_note}; "
        f"{n_reuse} existing run(s) reused without rerunning.")
    log(f"Reference run: existing {BASELINE_ID} ({'; '.join(base_notes) if not base_problems else 'CANNOT be used: ' + '; '.join(base_problems)})")
    for j in preview:
        gw_name, _ = expected_files(settings_for(j))
        what = "REUSE (existing run)" if j["kind"] == "reuse" else "new run"
        log(f"  {j['run_id']:44s} {describe(j):18s} {what}")
    extra_runs = n_run + (0 if ("ks" not in knob_keys or args.no_refine or man.get("ks_refine")) else REFINE_N)
    log(f"  tRIBS time alone: about {extra_runs * WIN['hours'] * MIN_PER_1000H / 1000:.0f}-{extra_runs * 13:.0f} min "
        f"(W2 runs took 10-13 min); each run also needs a notebook build of about 1 min.")

    results = preflight(args, lab, preview, base_problems, nb_note, man)
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

    man["script"] = "run_screen_soil_sycamore.py"
    man["baseline_id"] = BASELINE_ID
    man["last_started"] = now_iso()
    save_manifest(man)

    install_term_handler()

    t_batch, code = time.time(), 0

    def run_batch(todo: list, ref_in: Path) -> int:
        consecutive_fail, rc = 0, 0
        for j in todo:
            rid = j["run_id"]
            try:
                do_run(j, args, lab, man, ref_in)
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
        base_job = make_job(None)
        if not do_reuse(base_job, args, lab, man, None):
            log("The W2 baseline cannot be the reference run, so nothing was run. Nothing is rerun automatically: "
                "tell me what the messages above say and we decide together.")
            return 1
        ref_in = lab / f"{BASELINE_ID}.in"
        todo = preview
        if args.limit is not None:
            def finished(j):
                rec = man["runs"].get(j["run_id"], {})
                return rec.get("status") == "done" and same_record(rec.get("settings"), settings_record(j)) and j["run_id"] not in args.redo_ids
            # the reused runs cost nothing, so they are always checked and scored; --limit counts only runs that need tRIBS
            fresh = [j for j in todo if j["kind"] == "run" and not finished(j)][:args.limit]
            todo = [j for j in todo if j["kind"] == "reuse"] + fresh
            if not fresh:
                log("Nothing left to run: every run in this selection is already done.")
        code = run_batch(todo, ref_in)
        if (code == 0 and "ks" in knob_keys and not args.no_refine and args.limit is None and not args.build_only):
            if not man.get("ks_refine"):
                info, why = pick_ks_refinement(load_scores(lab), man)
                if info:
                    man["ks_refine"] = info
                    save_manifest(man)
                    log()
                    log(f"Ks second pass: the {info['metric']} changes most between Ks x{info['interval'][0]:g} and "
                        f"x{info['interval'][1]:g} ({info['change']:.3f}), so the 3 extra values are {info['values']}.")
                else:
                    log(f"Ks second pass skipped: {why}")
            if man.get("ks_refine"):
                extra = [make_job(knob_by_key("ks"), v) for v in man["ks_refine"]["values"]]
                code = run_batch(extra, ref_in) or code
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
    log(f"Batch time: {fmt_min(time.time() - t_batch)}")
    log(f"{'run':46s} {'status':16s} {'build':>7s} {'tRIBS':>7s}")
    for rid, rec in man["runs"].items():
        log(f"{rid:46s} {rec.get('status', '?'):16s} {str(rec.get('build_min', '')):>7s} {str(rec.get('tribs_min', '')):>7s}")
    print_summary(man, lab, knob_keys)
    return code


if __name__ == "__main__":
    sys.exit(main())
