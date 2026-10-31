#!/usr/bin/env python3
"""
run_screen_forcing_sycamore.py  --  Stage 2, piece 3: the forcing probes on W2.

Question it answers: when ONE input to the model is changed (the rain, the air temperature, the snow module, the weather clock
or the number of weather stations) and everything else stays at the baseline, how much do the hydrograph targets change on W2
(Jan 18 - Mar 12 2017, 1,272 h, water table at 0.82 x bedrock)? These are PROBES, not calibration knobs: the point is to learn
how much each input could matter, and whether any of them can explain what the soil knobs could not (the Feb 27-28 / Feb 18-19
balance, the early Feb 18-19 arrival, the recession tails).

The eight probes (the baseline is never run again; it is the existing W2_ks0p2_cv6_psib0p5_nocl run):

    hour0    weather hour label 0 instead of -1      weather notebook rerun with LABEL_SHIFT_HOURS = 0 (every weather row moves 1 h later)
    rain1p2  rain x1.2                               a scaled copy of the 1,272 W2 radar hours
    snowoff  snow module off                         OPTSNOW = 0 in the build notebook
    tempp2   +2 C at the highest station             weather notebook rerun; air temperature +2 C at the highest NLDAS station only
    tempm2   -2 C at the highest station             same, -2 C
    onestn   one NLDAS station instead of three      weather notebook rerun with one station (default: the middle-east cell, 1,247 m)
    rain0p8  rain x0.8
    rain1p4  rain x1.4

Rain. The radar grids are ESRI ASCII files, one per hour, labeled by the END of the hour. For each rain probe the script makes a
scaled copy of exactly the 1,272 hours W2 needs, in ~/forcing_probe_rain/rain_x<factor>/. The header of every file is copied
unchanged and every rain cell is multiplied by the factor. It then re-reads all 1,272 files against the originals and checks each
header, each cell (to within the 4-decimal rounding) and the total. The ORIGINAL radar folder is only read, never written.

Weather. A weather probe reruns the weather notebook (Generate_Met_Forcing_Sycamore.ipynb) headlessly, with a few lines changed
IN MEMORY (the notebook file is never changed): LABEL_SHIFT_HOURS, STATIONS (one-station probe), and the output folder, so the
files go to ../init_data/met_probe_<probe>/ and the live ../init_data/met/ files are never touched. For the temperature probes one cell
is added in memory that changes the air temperature at the highest station; by default relative humidity is kept as it was, so
only the air temperature changes in the files (--temp-keeps q keeps specific humidity instead). The NLDAS downloads are cached in
~/nldas_cache; the script blocks every download (it never asks for the Earthdata login, and no credentials are involved): if a piece
is missing from the cache it stops and says so.
BEFORE the first weather probe the script proves the route works: it reruns the notebook with the baseline settings and requires
the new weather files to be BYTE-IDENTICAL to the live ones. If they are not, no weather probe is run and it says how they differ.
After each weather rerun it compares the new files with the live ones and requires that ONLY the intended change is there (for
example: +2.00 C on every row of one station and nothing else; or every row equal to the baseline's row one hour earlier).

For each probe the script:
    1. Makes and checks the changed forcing files (above).
    2. Copies the Settings cell of the Make_Sycamore_Model notebook, changes ONLY that probe's lines (plus the run names and the
       window), and runs the whole notebook headlessly (no browser) with nbclient. The notebook file is never changed.
    3. Checks the build: the new .in file must match the baseline's .in in every keyword except the run-specific ones (RUNTIME,
       STARTDATE, OUTFILENAME, GWATERFILE, BEDROCKFILE) and the one line the probe changes (RAINFILE, HYDROMETSTATIONS, or OPTSNOW
       with SNOWFILENAME). It also checks that the soil grids are the baseline's.
    4. Runs tRIBS (one at a time, never two at once, in the foreground) and checks that it finished.
    5. Scores the run with score_run_sycamore.py (window W2) and notes a run whose outlet file is byte-identical to the baseline's.
    6. Writes everything to screen_forcing_manifest.json, so a stopped batch can be resumed by running the same command.

NO RUN IS REPEATED and nothing in results/ is ever overwritten: a planned run whose name already has files in results/ stops the
start-up check. The W2 baseline is checked read-only (outlet file, checksum, md5, .in file) and then only scored.

Run it from workspaces/Sycamore/lab/. Typical use:
    python -u run_screen_forcing_sycamore.py --dry-run          # checks only, builds nothing (includes a one-minute cache check)
    python -u run_screen_forcing_sycamore.py                    # everything; run it again to resume after any stop
    python -u run_screen_forcing_sycamore.py --summary-only     # just print the table from what is finished
    python -u run_screen_forcing_sycamore.py --only hour0,rain1p2 --build-only    # make the files and build 2 runs, no tRIBS
    python -u run_screen_forcing_sycamore.py --limit 1          # make the next probe only
    python -u run_screen_forcing_sycamore.py --one-station north-east   # choose a different cell for the one-station probe

Files it makes (all in lab/ unless noted): screen_forcing_manifest.json, screen_forcing.lock (only while running),
screen_forcing_logs/, score_out_screen_forcing/scores_screen_forcing.csv, screen_forcing_summary.csv, and for each probe
<run>.in, simulation_<run>.log and results/<run>*. Outside lab/: ../init_data/met_probe_<probe>/ (weather probes) and
~/forcing_probe_rain/ and ~/forcing_probe_weather/ (scaled radar copies; unused NLDAS rain files).

Do NOT run the notebook (or any cell that rebuilds data/model/) while this runs: every build rewrites lab/data/model/, and all the
runs share it. Keep the terminal open: the progress line printed every 5 minutes also stops the Codespace from timing out for being
idle. If the connection drops, run the same command again: finished runs are skipped.
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
# What the probes do
# --------------------------------------------------------------------------------------------
WIN = dict(name="W2", start="01/18/2017/00/00", hours=1272, label="jan_to_mar", gw=0.82)   # label = the RUN_CHOICE text printed

# Everything that is not being probed (the Oct 21 baseline, with the W2 cold-start water table).
BASE = {
    "GW_FRACTION": WIN["gw"], "BEDROCK_ADD_M": 1.0, "KS_MULT": 0.2, "PSIB_MULT": 0.5, "F_MULT": 1,
    "THETAS_MULT": 1, "ANISOTROPY": 1, "M_MULT": 1,
    "KINEMVELCOEF": 6, "FLOWEXP": 0.3, "CHANNELROUGHNESS": 0.04, "CHANNELWIDTHCOEFF": 2.33,
    "OPTSNOW": 1, "OPTPERCOLATION": 0, "CHANNELCONDUCTIVITY": 5,
    "REBUILD_SOIL": False,
}
ROUTING_KWS = ["KINEMVELCOEF", "FLOWEXP", "CHANNELROUGHNESS", "CHANNELWIDTHCOEFF"]   # must stay at baseline in every probe

BASELINE_ID = "W2_ks0p2_cv6_psib0p5_nocl"              # the existing W2 baseline run (from the proxy check); reused, never rerun
PROXY_MANIFEST = "proxy_manifest.json"
NOTEBOOK_CANDIDATES = ["Make_Sycamore_Model.ipynb", "Make_Sycamore_Model_clean.ipynb"]
WEATHER_NOTEBOOK = "Generate_Met_Forcing_Sycamore.ipynb"
SCORER = "score_run_sycamore.py"
TRIBS_DEFAULT = "/opt/tRIBS/build/tRIBS"
INIT_DATA = "../init_data"           # SOLUS_DIR in the notebook
MODEL_SOIL = "data/model/soil"
LOG_DIR = "screen_forcing_logs"
SCORE_DIR = "score_out_screen_forcing"
SCORE_CSV = "scores_screen_forcing.csv"
SUMMARY_CSV = "screen_forcing_summary.csv"
MANIFEST = "screen_forcing_manifest.json"
LOCK = "screen_forcing.lock"
MIN_PER_1000H = 7.8          # measured tRIBS speed on this machine: minutes per 1,000 model hours (W2 took 10.0-12.9 min)
HEARTBEAT_MIN = 5.0
ROWS_PER_MODEL_HOUR = 16     # the outlet file gets 16 rows per model hour
F_FLOOR = 1e-4               # the notebook keeps f at or above this
THETAS_CAP = 0.99            # the notebook keeps thetaS and porosity below this

# Where the probe inputs go. Nothing here is inside the original radar folder or the original weather folder.
RAIN_OUT = "~/forcing_probe_rain"                    # rain_x1p2/, rain_x0p8/, rain_x1p4/ (copies of the W2 radar hours, scaled)
WEATHER_SCRATCH = "~/forcing_probe_weather"          # the weather notebook also writes NLDAS rain files nobody uses; they go here
LIVE_MET_SDF = "../init_data/met/Master_Met.sdf"     # the weather files the baseline run uses (never touched)
MET_PROBE_DIR = "../init_data/met_probe_{tag}"       # one folder per weather probe; Master_Met.sdf and the .mdf files go inside
BASE_SHIFT = -1                                      # the baseline's weather hour-label shift (LABEL_SHIFT_HOURS), from the Oct 30 handoff
DEFAULT_ONE_STATION = "middle-east"                  # the NLDAS cell kept for the one-station probe (see --one-station)
CHECK_TAG = "check"                                  # the baseline-regeneration test writes to met_probe_check
DT_TOL = 0.011                                       # C: the weather files are rounded to 0.01, so +2.00 may show as 1.99 to 2.01

# The probes, in the order they are made. The first two exercise the two new mechanisms (a weather rebuild, then a rain
# copy) so a problem shows up in the first half hour. tag = what goes in the run name.
PROBES = [
    dict(key="hour0",   kind="weather", tag="hour0",   shift=0,   label="weather hour label 0 instead of -1"),
    dict(key="rain1p2", kind="rain",    tag="rain1p2", factor=1.2, label="rain x1.2"),
    dict(key="snowoff", kind="snow",    tag="snowoff",            label="snow module off (OPTSNOW 0)"),
    dict(key="tempp2",  kind="weather", tag="tempp2",  dT=+2.0,   label="+2 C at the highest station"),
    dict(key="tempm2",  kind="weather", tag="tempm2",  dT=-2.0,   label="-2 C at the highest station"),
    dict(key="onestn",  kind="weather", tag="onestn",  one_station=True, label="one NLDAS station instead of three"),
    dict(key="rain0p8", kind="rain",    tag="rain0p8", factor=0.8, label="rain x0.8"),
    dict(key="rain1p4", kind="rain",    tag="rain1p4", factor=1.4, label="rain x1.4"),
]

# Score columns: (short name, scorer column). The first group is shown in the tables.
EVENTS = [("fw", "flood_week", "FLOOD WEEK Jan 18-25"), ("s2", "jan_second_storm", "2ND JAN STORM"),
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


def probe_by_key(key: str) -> dict:
    return next(p for p in PROBES if p["key"] == key)


def run_id_for(probe: dict | None) -> str:
    """W2_ks0p2_cv6_psib0p5_nocl is the baseline. A probe goes after psib, before the loss setting, as in the soil screen."""
    if probe is None:
        return BASELINE_ID
    return f"{WIN['name']}_ks{vtag(BASE['KS_MULT'])}_cv6_psib{vtag(BASE['PSIB_MULT'])}_{probe['tag']}_nocl"


def make_job(probe: dict | None) -> dict:
    if probe is None:
        return dict(probe=None, key=None, knob=None, kw=None, value=None, run_id=BASELINE_ID, kind="baseline")
    return dict(probe=probe, key=probe["key"], knob=probe["key"], kw=None, value=None, run_id=run_id_for(probe), kind="run")


def rain_dir_for(probe: dict) -> str:
    """Absolute folder of the scaled radar copy for a rain probe (it does not have to exist yet)."""
    return str(Path(os.path.expanduser(RAIN_OUT)) / f"rain_x{vtag(probe['factor'])}")


def met_sdf_for(tag: str) -> str:
    """Station list of a weather probe, written the same relative way as the baseline's ('../init_data/...')."""
    return MET_PROBE_DIR.format(tag=tag) + "/Master_Met.sdf"


def radar_name(label: dt.datetime) -> str:
    return f"radar_{label:%m%d%Y%H}.asc"


def window_labels() -> list:
    """The end-labeled radar hours a W2 run needs: start+1 h ... start+hours."""
    start = dt.datetime.strptime(WIN["start"], "%m/%d/%Y/%H/%M")
    return [start + dt.timedelta(hours=h) for h in range(1, WIN["hours"] + 1)]


def settings_for(job: dict) -> dict:
    """Every NAME = value line of the Settings cell that this run sets."""
    s = dict(BASE)
    s.update(RUN_CHOICE=WIN["label"], INPUT_FILE=f"{job['run_id']}.in", RESULTS_NAME=job["run_id"])
    p = job.get("probe")
    if p is None:
        return s
    if p["kind"] == "rain":
        labels = window_labels()
        s.update(RADAR_DIR=rain_dir_for(p), RADAR_FIRST_FILE=radar_name(labels[0]), RADAR_LAST_FILE=radar_name(labels[-1]))
    elif p["kind"] == "weather":
        s.update(MET_STATIONS=met_sdf_for(p["tag"]))
    elif p["kind"] == "snow":
        s.update(OPTSNOW=0)
    return s


RECORD_KEYS = ("GW_FRACTION", "BEDROCK_ADD_M", "KS_MULT", "PSIB_MULT", "F_MULT", "M_MULT", "THETAS_MULT", "ANISOTROPY", "KINEMVELCOEF",
               "FLOWEXP", "CHANNELROUGHNESS", "CHANNELWIDTHCOEFF", "OPTSNOW", "OPTPERCOLATION", "CHANNELCONDUCTIVITY",
               "RADAR_DIR", "MET_STATIONS")


def settings_record(job: dict, args=None) -> dict:
    """What the manifest remembers about the settings, to notice a changed setting on resume."""
    s = settings_for(job)
    rec = {k: s[k] for k in RECORD_KEYS if k in s}
    rec.update(START_DATE=WIN["start"], RUN_HOURS=WIN["hours"])
    p = job.get("probe")
    if p is not None:
        for k in ("factor", "shift", "dT"):
            if k in p:
                rec["probe_" + k] = p[k]
        if p.get("one_station") and args is not None:
            rec["probe_station"] = args.one_station
        if "dT" in p and args is not None:
            rec["probe_temp_keeps"] = args.temp_keeps
    return rec


def same_record(old, cur) -> bool:
    return bool(old) and old == cur


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
    p = job.get("probe")
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
    if p and p["kind"] == "rain":
        must += [f"Radar rain      : {s['RADAR_DIR']}/radar_*.asc",
                 f"This run needs {WIN['hours']} files: first {s['RADAR_FIRST_FILE']}, last {s['RADAR_LAST_FILE']}"]
    if p and p["kind"] == "weather":
        must.append(f"Weather stations: {s['MET_STATIONS']}")
    missing = [m for m in must if m not in text]
    if not re.search(r"All \d+ checks passed\.", text):
        missing.append("All <n> checks passed.")
    if ("OK: thetaS is above theta*_s everywhere in the basin." not in text
            and "NOTE: the unscaled SOLUS thetaS is already at or below theta*_s" not in text):
        missing.append("the thetaS against theta*_s check (OK: or NOTE: line)")
    n_expect = 1 if (p and p.get("one_station")) else len(read_sdf_rows(LIVE_MET_SDF))
    n_ok = len(re.findall(r"^OK\s+station \d+:", text, re.M))
    if n_ok != n_expect:
        missing.append(f"{n_expect} line(s) 'OK  station <n>: ...' from the weather-station check (found {n_ok})")
    if missing:
        raise RunError("The notebook output is missing line(s) it should have printed:\n  " + "\n  ".join(missing)
                       + "\n  (If these are the thetaS or anisotropy lines, the notebook in lab/ is the old one: replace it with the new one.)")
    if re.search(r"^PROBLEM:", text, re.M):
        raise RunError("The notebook printed a PROBLEM line (see the build output file).")
    bad = re.findall(r"^(?:CHECK|MISSING)\s+(?:station|Radar|Weather)[^\n]*$", text, re.M)
    if bad:
        raise RunError("The notebook flagged the forcing files:\n  " + "\n  ".join(b.strip() for b in bad[:6]))
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


def same_path(a, b) -> bool:
    """Two paths are the same if they point to the same place (a relative path counts from the current folder, lab/)."""
    try:
        return os.path.normpath(os.path.abspath(os.path.expanduser(str(a).strip()))) == \
               os.path.normpath(os.path.abspath(os.path.expanduser(str(b).strip())))
    except (TypeError, ValueError):
        return False


def check_in_file(new_in: Path, ref_in: Path, job: dict) -> list:
    """The new .in must equal the reference .in except for the run-specific lines and the ONE forcing line this probe changes
    (RAINFILE for a rain probe, HYDROMETSTATIONS for a weather probe, OPTSNOW and SNOWFILENAME for snow off). Routing keywords
    must have exactly their baseline values. Returns notes; raises RunError on any other difference."""
    s = settings_for(job)
    p = job.get("probe")
    gw_name, bed_name = expected_files(s)
    new, ref = read_in_file(new_in), read_in_file(ref_in)
    expect = {"RUNTIME": WIN["hours"], "STARTDATE": WIN["start"], "OUTFILENAME": f"results/{job['run_id']}",
              "GWATERFILE": f"{MODEL_SOIL}/{gw_name}", "BEDROCKFILE": f"{MODEL_SOIL}/{bed_name}"}
    for kw in ROUTING_KWS:
        expect[kw] = s[kw]
    expect_path, allow, notes = {}, set(), []
    if p and p["kind"] == "rain":
        expect_path["RAINFILE"] = f"{s['RADAR_DIR']}/radar_"
    elif p and p["kind"] == "weather":
        expect_path["HYDROMETSTATIONS"] = s["MET_STATIONS"]
    elif p and p["kind"] == "snow":
        expect["OPTSNOW"] = 0
        allow.add("SNOWFILENAME")            # the snow-parameter file is only written when the snow module is on
    problems = []
    for key in sorted(set(ref) | set(new)):
        b, n = ref.get(key), new.get(key)
        if key in expect:
            if n is None or not same_value(n, expect[key]):
                problems.append(f"{key}: written {n!r}, expected {expect[key]!r}")
        elif key in expect_path:
            if n is None or not same_path(n, expect_path[key]):
                problems.append(f"{key}: written {n!r}, expected {expect_path[key]!r}")
        elif key in allow:
            if b != n:
                notes.append(f"{key} differs from the baseline's ({b!r} -> {n!r}); allowed for this probe")
        elif b is None or n is None or not same_value(b, n):
            problems.append(f"{key}: reference {b!r} but this run {n!r}  (this keyword may not change)")
    if problems:
        raise RunError("The .in file differs from the reference in more than it should:\n  " + "\n  ".join(problems))
    return notes


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




# --------------------------------------------------------------------------------------------
# Rain probes: scaled copies of the W2 radar hours (the original radar folder is only ever READ)
# --------------------------------------------------------------------------------------------
def rain_source(ref_in: Path) -> tuple:
    """(folder, base name, extension) of the radar files the baseline run reads, from its .in file."""
    kw = read_in_file(ref_in)
    rf, ext = kw.get("RAINFILE"), kw.get("RAINEXTENSION")
    if not rf or not ext:
        raise RunError(f"{ref_in.name} has no RAINFILE / RAINEXTENSION line, so the radar folder is unknown.")
    folder, base = os.path.split(os.path.expanduser(rf))
    if base != "radar_" or ext.strip() != "asc":
        raise RunError(f"The baseline reads radar files named {base}<MMDDYYYYHH>.{ext}; this script expects radar_<MMDDYYYYHH>.asc "
                       f"(the build notebook's RADAR_BASENAME and RADAR_EXT).")
    return folder, base, ext.strip()


def read_grid_text(path: Path) -> tuple:
    """ESRI ASCII grid -> (the 6 header lines as text, header dict, 2-D array). Raises RunError if it is malformed."""
    lines = Path(path).read_text().splitlines()
    header = lines[:6]
    try:
        hdr = {ln.split()[0].lower(): float(ln.split()[1]) for ln in header}
        body = np.array(" ".join(lines[6:]).split(), dtype=float)
        ncols, nrows = int(hdr["ncols"]), int(hdr["nrows"])
    except (IndexError, ValueError, KeyError) as exc:
        raise RunError(f"{Path(path).name}: not a readable ESRI ASCII grid ({exc}).")
    if len(header) < 6 or body.size != ncols * nrows:
        raise RunError(f"{Path(path).name}: expected {ncols * nrows} values, found {body.size}.")
    return header, hdr, body.reshape(nrows, ncols)


def write_scaled_grid(src: Path, dst: Path, factor: float) -> None:
    header, hdr, arr = read_grid_text(src)
    nd = hdr.get("nodata_value", -9999)
    if (arr == nd).any() or (arr < 0).any() or not np.isfinite(arr).all():
        raise RunError(f"{src.name}: has NODATA, negative or non-finite cells in the data; the radar files are expected to have none.")
    with open(dst, "w") as fh:
        fh.write("\n".join(header) + "\n")
        np.savetxt(fh, arr * factor, fmt="%.4f")


def verify_rain_dir(dest: Path, folder: str, names: list, factor: float) -> dict:
    """Every scaled file must exist, share its source's header, and equal (source x factor) to within the 4-decimal rounding.
    No other file may be in the folder."""
    if not dest.is_dir():
        raise RunError(f"{dest} does not exist.")
    extra = sorted(set(os.listdir(dest)) - set(names))
    if extra:
        raise RunError(f"{dest} has {len(extra)} file(s) that are not W2 radar hours, e.g. {extra[:3]}.")
    max_err, src_sum, dst_sum, cells = 0.0, 0.0, 0.0, 0
    for n in names:
        f = dest / n
        if not f.is_file():
            raise RunError(f"{dest} is missing {n}.")
        h0, _, a0 = read_grid_text(Path(folder) / n)
        h1, _, a1 = read_grid_text(f)
        if h0 != h1:
            raise RunError(f"{n}: the header differs from the original's.")
        err = float(np.max(np.abs(a1 - a0 * factor)))
        max_err = max(max_err, err)
        if err > 1.1e-4:
            raise RunError(f"{n}: a cell differs from (original x {factor:g}) by {err:.2e}.")
        src_sum += float(a0.sum()); dst_sum += float(a1.sum()); cells = a0.size
    ratio = dst_sum / src_sum if src_sum > 0 else float("nan")
    if not (abs(ratio / factor - 1) < 2e-4):
        raise RunError(f"Total rain in {dest.name} is {ratio:.5f} x the original's; expected {factor:g}.")
    return dict(files=len(names), max_err=max_err, ratio=ratio, src_mean_mm=src_sum / cells, dst_mean_mm=dst_sum / cells)


def prep_rain(probe: dict, ref_in: Path) -> str:
    folder, base, ext = rain_source(ref_in)
    names = [f"{base}{lab:%m%d%Y%H}.{ext}" for lab in window_labels()]
    missing = [n for n in names if not (Path(folder) / n).is_file()]
    if missing:
        raise RunError(f"{len(missing)} W2 radar file(s) are missing from {folder}, e.g. {missing[:3]}.")
    f, dest = float(probe["factor"]), Path(rain_dir_for(probe))
    root = Path(os.path.expanduser(RAIN_OUT))
    try:
        stats, how = verify_rain_dir(dest, folder, names, f), "found and checked again, not rewritten"
    except RunError:
        root.mkdir(parents=True, exist_ok=True)
        tmp = dest.with_name(dest.name + ".tmp")
        for old in (tmp, dest):
            if old.exists():
                if old.parent.resolve() != root.resolve() or not old.name.startswith("rain_x"):
                    raise RunError(f"Refusing to delete {old}: it is not one of this script's own scaled folders.")
                shutil.rmtree(old)
        tmp.mkdir()
        for n in names:
            write_scaled_grid(Path(folder) / n, tmp / n, f)
        stats = verify_rain_dir(tmp, folder, names, f)
        tmp.rename(dest)
        how = "written and checked"
    return (f"rain x{f:g}: {stats['files']} W2 radar files, {how} in {dest}. Grid-mean rain over the window "
            f"{stats['src_mean_mm']:.1f} mm -> {stats['dst_mean_mm']:.1f} mm (ratio {stats['ratio']:.4f}); "
            f"largest cell difference from exactly x{f:g}: {stats['max_err']:.1e}. Original radar folder only read: {folder}")


# --------------------------------------------------------------------------------------------
# Weather probes: the weather notebook is rerun headlessly with a few lines changed in memory (the file is never changed)
# --------------------------------------------------------------------------------------------
MDF_TIME = ["Year", "Month", "Day", "Hour"]
MDF_DATA = ["PA_mb", "RH_pct", "XC_tenths", "US_m/s", "TA_C", "IS_W/m2", "TS_C"]
WEATHER_BLOCKED = None            # set to a message when the baseline-regeneration test fails; weather probes then wait


def read_sdf_rows(path) -> list:
    """Rows of a tRIBS station list: ID,DataFile,Northing,Easting,Elevation (one header line first)."""
    lines = Path(path).read_text().splitlines()
    rows = []
    for ln in lines[1:]:
        if not ln.strip():
            continue
        p = [x.strip() for x in ln.split(",")]
        if len(p) < 5:
            raise RunError(f"{path}: cannot read the station line {ln!r}.")
        rows.append(dict(id=p[0], file=p[1], north=p[2], east=p[3], elev=p[4]))
    if not rows:
        raise RunError(f"{path} lists no stations.")
    return rows


def read_mdf(path):
    import pandas as pd
    df = pd.read_csv(path)
    missing = [c for c in MDF_TIME + MDF_DATA if c not in df.columns]
    if missing:
        raise RunError(f"{path}: missing column(s) {missing}.")
    return df


def find_weather_settings_cell(nb) -> int:
    for i, cell in enumerate(nb.cells):
        if cell.cell_type != "code":
            continue
        s = cell.source
        if all(re.search(rf"^{n}[ \t]*=", s, re.M) for n in ("LABEL_SHIFT_HOURS", "STATIONS", "MET_SDF", "PRECIP_SDF")):
            return i
    raise RunError("Could not find the Settings cell of the weather notebook (the cell with LABEL_SHIFT_HOURS, STATIONS, MET_SDF and PRECIP_SDF).")


def weather_live_settings(lab: Path, wnb: str) -> dict:
    """Read the weather notebook's Settings as plain values (no code is run)."""
    import nbformat
    nb = nbformat.read(str(lab / wnb), as_version=4)
    out = {}
    for node in ast.parse(nb.cells[find_weather_settings_cell(nb)].source).body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            try:
                out[node.targets[0].id] = ast.literal_eval(node.value)
            except ValueError:
                pass
    return out


def weather_spec(probe: dict | None, args, live: dict) -> dict:
    """What to change in the weather notebook for a probe (probe None = the baseline-regeneration test)."""
    tag = CHECK_TAG if probe is None else probe["tag"]
    spec = dict(tag=tag, shift=BASE_SHIFT, stations=None, dT=None, keep=args.temp_keeps,
                met_sdf=met_sdf_for(tag), precip_sdf=f"{WEATHER_SCRATCH}/{tag}/Precip_unused.sdf")
    if probe is not None:
        spec["shift"] = probe.get("shift", BASE_SHIFT)
        spec["dT"] = probe.get("dT")
        if probe.get("one_station"):
            if args.one_station not in live.get("STATIONS", []):
                raise RunError(f"--one-station {args.one_station!r} is not one of the notebook's STATIONS {live.get('STATIONS')}.")
            spec["stations"] = [args.one_station]
    return spec


def build_weather_nb(lab: Path, wnb: str, spec: dict, mode: str = "full"):
    """The weather notebook with its Settings patched in memory, the download blocked, and (for a temperature probe) one cell added.
    mode 'check' stops right after the download cell: it only counts how many NLDAS pieces are missing from the cache."""
    import nbformat
    from nbformat.v4 import new_code_cell
    nb = nbformat.read(str(lab / wnb), as_version=4)
    for cell in nb.cells:
        if cell.cell_type == "code":
            cell.outputs = []
            cell.execution_count = None
    idx = find_weather_settings_cell(nb)
    values = {"LABEL_SHIFT_HOURS": spec["shift"], "MET_SDF": spec["met_sdf"], "PRECIP_SDF": spec["precip_sdf"]}
    if spec.get("stations"):
        values["STATIONS"] = list(spec["stations"])
    src = patch_settings(nb.cells[idx].source, values)
    verify_patched(src, values)
    nb.cells[idx].source = src
    # No downloads, ever: if anything is missing from the cache the notebook stops before it would ask for the Earthdata login.
    gi = [i for i, c in enumerate(nb.cells) if c.cell_type == "code" and "earthaccess.login" in c.source]
    pat = re.compile(r"^if todo:[ \t]*$", re.M)
    if len(gi) != 1 or len(pat.findall(nb.cells[gi[0]].source)) != 1:
        raise RunError("The weather notebook's download cell is not in the expected form (one cell with 'earthaccess.login' and one line "
                       "'if todo:'), so the no-download guard cannot be installed. Nothing was run.")
    guard = ("if todo:\n    raise RuntimeError(f'DOWNLOAD GUARD: {len(todo)} NLDAS piece(s) are not in the cache. This script never "
             "downloads. Run the weather notebook yourself to fill the cache, then run this script again.')\nif todo:")
    nb.cells[gi[0]].source = pat.sub(lambda m: guard, nb.cells[gi[0]].source, count=1)
    if mode == "check":
        nb.cells = nb.cells[:gi[0] + 1]
        return nb
    if spec.get("dT") is not None:
        ci = [i for i, c in enumerate(nb.cells) if c.cell_type == "code" and "convert_and_write_nldas_timeseries(" in c.source]
        fi = [i for i, c in enumerate(nb.cells) if c.cell_type == "code" and "s['data'] = df" in c.source]
        if len(ci) != 1 or len(fi) != 1 or fi[0] >= ci[0]:
            raise RunError("The weather notebook's cells are not in the expected form (one cell that sets s['data'] = df, followed by one "
                           "cell that calls convert_and_write_nldas_timeseries), so the temperature change cannot be placed.")
        code = (
            "# --- forcing probe (added in memory by run_screen_forcing_sycamore.py) ---\n"
            f"_dT = {float(spec['dT'])!r}\n"
            f"_keep = {spec['keep']!r}\n"
            "_hi = max(stations, key=lambda s: s['elev'])\n"
            "_d = _hi['data']\n"
            "if _keep == 'rh':   # keep relative humidity as it is: scale specific humidity by the ratio of saturation vapor pressures\n"
            "    _es = lambda T: 6.11 * np.exp((2.453e6 / 461) * ((1 / 273.15) - (1 / T)))\n"
            "    _d['humidity'] = _d['humidity'] * _es(_d['temp'] + _dT) / _es(_d['temp'])\n"
            "_d['temp'] = _d['temp'] + _dT\n"
            "print(f\"PROBE: air temperature {_dT:+g} C at {_hi['label']} (NLDAS_{_hi['n']}, {_hi['elev']} m); \"\n"
            "      + ('relative humidity kept as it was' if _keep == 'rh' else 'specific humidity kept, so relative humidity changes'))\n")
        nb.cells.insert(ci[0], new_code_cell(code))
    return nb


def run_weather_nb(lab: Path, args, spec: dict, mode: str, log_dir: Path) -> str:
    """Run the patched weather notebook. Returns everything it printed; raises RunError on a failure or a problem line."""
    import nbformat
    nb = build_weather_nb(lab, args.weather_notebook, spec, mode)
    name = f"weather_{spec['tag']}_{mode}"
    try:
        execute_nb(nb, lab, args.cell_timeout)
    except Exception as exc:
        partial = outputs_text(nb)
        (log_dir / f"{name}_output.txt").write_text(partial)
        try:
            nbformat.write(nb, str(log_dir / f"{name}_FAILED_executed.ipynb"))
        except Exception:
            pass
        tail = re.sub(r"\x1b\[[0-9;]*m", "", "\n".join(partial.splitlines()[-10:]))
        raise RunError(f"The weather notebook stopped with an error: {type(exc).__name__}: {str(exc)[:500]}\n  (last output lines follow)\n{tail}")
    text = outputs_text(nb)
    (log_dir / f"{name}_output.txt").write_text(text)
    m = re.search(r"(\d+) still to download", text)
    if not m:
        raise RunError("The weather notebook did not print its '... still to download' line, so the cache could not be checked.")
    if int(m.group(1)) != 0:
        raise RunError(f"{m.group(1)} NLDAS piece(s) are not in the cache (~/nldas_cache); this script never downloads.")
    if mode == "full":
        if "No problems found." not in text or "PROBLEMS FOUND" in text:
            raise RunError("The weather notebook's own checks (its section 7) reported problems; see " + str(log_dir / f"{name}_output.txt"))
        if spec.get("dT") is not None and "PROBE: air temperature" not in text:
            raise RunError("The temperature change cell did not run (no 'PROBE: air temperature' line in the output).")
    return text


def describe_mdf_difference(pa: str, pb: str) -> str:
    """Say how two weather files differ (used when the regenerated baseline does not match the live files)."""
    a, b = read_mdf(pa), read_mdf(pb)
    if len(a) != len(b):
        return f"{len(a)} rows in the live file but {len(b)} in the regenerated one"
    out = []
    for c in MDF_DATA:
        d = float(np.max(np.abs(a[c].to_numpy(float) - b[c].to_numpy(float))))
        if d > 0:
            out.append(f"{c} differs by up to {d:g}")
    if not out:
        out.append("the numbers are the same, so only the file text differs (for example the date columns or line endings)")
    x, y = a[MDF_DATA].to_numpy(float), b[MDF_DATA].to_numpy(float)
    if np.array_equal(y[:-1], x[1:]):
        out.append("HINT: the regenerated file equals the live file moved one hour earlier, which is what you get when the live files were "
                   f"made with a label shift of {BASE_SHIFT + 1}, not {BASE_SHIFT}")
    elif np.array_equal(y[1:], x[:-1]):
        out.append("HINT: the regenerated file equals the live file moved one hour later, which is what you get when the live files were "
                   f"made with a label shift of {BASE_SHIFT - 1}, not {BASE_SHIFT}")
    return "; ".join(out)


def compare_baseline_copy(new_sdf: str) -> str:
    """The regenerated baseline weather files must be BYTE-identical to the live ones, and list the same stations."""
    live, new = read_sdf_rows(LIVE_MET_SDF), read_sdf_rows(new_sdf)
    if len(live) != len(new):
        raise RunError(f"The regenerated station list has {len(new)} station(s); the live Master_Met.sdf has {len(live)}.")
    for i, (a, b) in enumerate(zip(live, new), start=1):
        for k in ("north", "east", "elev"):
            if not same_value(a[k], b[k]):
                raise RunError(f"Station {i}: {k} is {b[k]} in the regenerated list but {a[k]} in the live one.")
        if not Path(a["file"]).is_file() or not Path(b["file"]).is_file():
            raise RunError(f"Station {i}: weather file not found ({a['file']} or {b['file']}).")
        if not filecmp.cmp(a["file"], b["file"], shallow=False):
            raise RunError(f"Station {i} ({a['elev']} m): the regenerated baseline weather file is NOT identical to the live one "
                           f"({Path(a['file']).name}): {describe_mdf_difference(a['file'], b['file'])}.")
    return (f"the notebook, rerun with the baseline settings (label shift {BASE_SHIFT}, no change), rewrote all {len(live)} weather files "
            f"byte-for-byte identical to the live ones in {os.path.dirname(LIVE_MET_SDF)}")


def weather_signature(lab: Path, wnb: str) -> str:
    h = hashlib.sha256()
    h.update((lab / wnb).read_bytes())
    for r in read_sdf_rows(LIVE_MET_SDF):
        h.update(r["file"].encode())
        h.update(Path(r["file"]).read_bytes())
    h.update(Path(LIVE_MET_SDF).read_bytes())
    return h.hexdigest()[:16]


def validate_weather_route(args, lab: Path, man: dict, log_dir: Path) -> str:
    """Before the first weather probe: prove that rerunning the weather notebook reproduces the baseline's weather files exactly.
    If it does not, no weather probe can be trusted, so none is run."""
    global WEATHER_BLOCKED
    if WEATHER_BLOCKED:
        raise RunError("Weather probes are blocked: " + WEATHER_BLOCKED)
    rec = man.setdefault("weather_validation", {})
    sig = weather_signature(lab, args.weather_notebook)
    chk = met_sdf_for(CHECK_TAG)
    if rec.get("ok") and rec.get("sig") == sig:
        try:
            return compare_baseline_copy(chk) + " (checked again from the files; the notebook was not rerun)"
        except (RunError, OSError):
            pass
    live = weather_live_settings(lab, args.weather_notebook)
    log(f"Weather route test: rerunning the weather notebook with the baseline settings (label shift {BASE_SHIFT}); its live value is "
        f"{live.get('LABEL_SHIFT_HOURS')!r}. This takes a few minutes and writes only to {os.path.dirname(chk)}.")
    t0 = time.time()
    try:
        run_weather_nb(lab, args, weather_spec(None, args, live), "full", log_dir)
        note = compare_baseline_copy(chk)
    except RunError as exc:
        WEATHER_BLOCKED = str(exc)
        rec.update(ok=False, sig=sig, problem=str(exc), when=now_iso())
        save_manifest(man)
        raise RunError("The weather route test FAILED, so no weather probe will be run. " + str(exc))
    rec.update(ok=True, sig=sig, note=note, minutes=round((time.time() - t0) / 60, 2), when=now_iso())
    save_manifest(man)
    return note + f" ({fmt_min(time.time() - t0)})"


def verify_probe_weather(probe: dict, sdf: str, args, live_stations: list) -> str:
    """Compare a probe's weather files with the live (baseline) ones and say what changed. Raises RunError if anything else changed."""
    live, new = read_sdf_rows(LIVE_MET_SDF), read_sdf_rows(sdf)
    if probe.get("one_station"):
        k = live_stations.index(args.one_station)
        if len(new) != 1 or k >= len(live):
            raise RunError(f"Expected one station in {sdf}; found {len(new)}.")
        a, b = read_mdf(live[k]["file"]), read_mdf(new[0]["file"])
        for c in ("north", "east", "elev"):
            if not same_value(live[k][c], new[0][c]):
                raise RunError(f"The one station's {c} is {new[0][c]} but the baseline's station {k + 1} has {live[k][c]}.")
        if a.shape != b.shape or not np.array_equal(a.to_numpy(float), b.to_numpy(float)):
            raise RunError(f"The one-station file is not identical to the baseline's station {k + 1} ({args.one_station}) file.")
        return (f"one station only: the {args.one_station} cell (NLDAS_{k + 1} in the baseline, {live[k]['elev']} m). Its file equals the baseline's "
                f"NLDAS_{k + 1} file row for row, so the other {len(live) - 1} station(s) are simply left out (nearest-station rule: "
                f"every node now takes this station's weather).")
    if len(new) != len(live):
        raise RunError(f"{sdf} lists {len(new)} stations; the baseline lists {len(live)}.")
    hi = max(range(len(live)), key=lambda i: float(live[i]["elev"]))
    parts = []
    for i, (la, ne) in enumerate(zip(live, new)):
        for c in ("north", "east", "elev"):
            if not same_value(la[c], ne[c]):
                raise RunError(f"Station {i + 1}: {c} is {ne[c]} but the baseline has {la[c]}.")
        a, b = read_mdf(la["file"]), read_mdf(ne["file"])
        if len(a) != len(b) or not np.array_equal(a[MDF_TIME].to_numpy(), b[MDF_TIME].to_numpy()):
            raise RunError(f"Station {i + 1}: the date and hour columns differ from the baseline's.")
        x, y = a[MDF_DATA].to_numpy(float), b[MDF_DATA].to_numpy(float)
        if "shift" in probe:
            step = int(probe["shift"] - BASE_SHIFT)
            if step < 1 or not np.array_equal(y[step:], x[:-step]):
                raise RunError(f"Station {i + 1}: expected every row to equal the baseline's row {step} h earlier; it does not.")
            parts.append(f"NLDAS_{i + 1}: every row equals the baseline's row {step} h earlier")
        else:
            dT = float(probe["dT"])
            if i != hi:
                if not filecmp.cmp(la["file"], ne["file"], shallow=False):
                    raise RunError(f"Station {i + 1} should be unchanged but its file differs from the baseline's.")
                parts.append(f"NLDAS_{i + 1} ({la['elev']} m): identical to the baseline")
                continue
            dta = b["TA_C"].to_numpy(float) - a["TA_C"].to_numpy(float)
            if np.max(np.abs(dta - dT)) > DT_TOL:
                raise RunError(f"Station {i + 1}: air temperature changed by {dta.min():+.3f} to {dta.max():+.3f}, expected {dT:+g} on every row.")
            same_cols = [c for c in MDF_DATA if c not in ("TA_C", "RH_pct")]
            if not np.array_equal(a[same_cols].to_numpy(float), b[same_cols].to_numpy(float)):
                raise RunError(f"Station {i + 1}: a column other than TA and RH changed.")
            drh = b["RH_pct"].to_numpy(float) - a["RH_pct"].to_numpy(float)
            if args.temp_keeps == "rh" and np.max(np.abs(drh)) > DT_TOL:
                raise RunError(f"Station {i + 1}: relative humidity changed by up to {np.max(np.abs(drh)):.3f} points but should have been kept.")
            parts.append(f"NLDAS_{i + 1} ({la['elev']} m, the highest): air temperature {dT:+g} C on all {len(a)} rows "
                         f"(measured {dta.min():+.2f} to {dta.max():+.2f}); relative humidity change {drh.min():+.2f} to {drh.max():+.2f} points "
                         f"(hours above 100%: {int((a['RH_pct'] > 100).sum())} before, {int((b['RH_pct'] > 100).sum())} after); "
                         f"pressure, wind and sunshine identical")
    return "; ".join(parts)


def prep_weather(probe: dict, args, lab: Path, man: dict, log_dir: Path) -> str:
    route = validate_weather_route(args, lab, man, log_dir)
    live = weather_live_settings(lab, args.weather_notebook)
    spec = weather_spec(probe, args, live)
    sdf = spec["met_sdf"]
    try:
        desc, how = verify_probe_weather(probe, sdf, args, list(live.get("STATIONS", []))), "found and checked again, notebook not rerun"
    except (RunError, OSError, ValueError, KeyError, IndexError):
        t0 = time.time()
        log(f"{probe['tag']}: rerunning the weather notebook ({probe['label']}) ...")
        run_weather_nb(lab, args, spec, "full", log_dir)
        desc, how = verify_probe_weather(probe, sdf, args, list(live.get("STATIONS", []))), f"notebook rerun in {fmt_min(time.time() - t0)}"
    return f"{probe['label']}: weather files in {os.path.dirname(sdf)} ({how}). {desc.rstrip('. ')}. Route test: {route}"


# --------------------------------------------------------------------------------------------
# Preflight
# --------------------------------------------------------------------------------------------
def describe(job: dict) -> str:
    return "baseline (no change)" if job.get("probe") is None else job["probe"]["label"]


def preflight(args, lab: Path, jobs: list, base_problems: list, nb_note: str, man: dict) -> list:
    """Checks before anything is built. Returns a list of (level, message)."""
    out = []

    def add(level, msg):
        out.append((level, msg))

    for f in (args.notebook, args.scorer, args.weather_notebook):
        add("OK" if (lab / f).is_file() else "FAIL", f"file in lab/: {f}")
    if nb_note:
        add("WARN", nb_note)
    if base_problems:
        add("FAIL", f"the W2 baseline {BASELINE_ID} cannot be the reference run: " + "; ".join(base_problems)
            + " (nothing was made, and nothing is rerun automatically)")
    else:
        add("OK", f"reference run {BASELINE_ID} passes its checks (reused, not rerun)")
    obs = lab / INIT_DATA / "obs" / "sycamore_09510200_obs_flow_WY2017.csv"
    add("OK" if obs.is_file() else "FAIL", f"gauge file: {INIT_DATA}/obs/{obs.name}")
    for f in ("../init_data/bedrock_depth_m.asc", "../init_data/scgrid.gdf", "../init_data/soils.sdt", "../init_data/LandUse.asc"):
        add("OK" if (lab / f).is_file() else "FAIL", f"input: {f}")
    live_ok = False
    try:
        live_rows = read_sdf_rows(lab / LIVE_MET_SDF)
        gone = [r["file"] for r in live_rows if not (lab / r["file"]).is_file()]
        live_ok = not gone
        add("FAIL" if gone else "OK", f"live weather stations: {len(live_rows)} listed in {LIVE_MET_SDF}"
            + (f"; missing files: {gone}" if gone else "; every station file is there"))
    except Exception as exc:
        add("FAIL", f"live weather station list {LIVE_MET_SDF}: {exc}")
    tribs = args.tribs
    add("OK" if (os.path.isfile(tribs) and os.access(tribs, os.X_OK)) else "FAIL", f"tRIBS program: {tribs}")
    free = disk_free_gb(lab)
    add("OK" if free >= args.min_free_gb else "FAIL", f"free disk space: {free:.1f} GB (need {args.min_free_gb:g})")
    comp = competing_tribs()
    add("FAIL" if comp else "OK", "no other tRIBS running" + (": " + comp[0] if comp else ""))
    # the build notebook: new enough, and every planned run patches cleanly
    try:
        import nbformat
        nb = nbformat.read(str(lab / args.notebook), as_version=4)
        idx = find_settings_cell(nb)
        src = nb.cells[idx].source
        old = []
        for name in ("M_MULT", "THETAS_MULT", "ANISOTROPY", "RADAR_DIR", "RADAR_FIRST_FILE", "RADAR_LAST_FILE", "MET_STATIONS", "OPTSNOW"):
            if not re.search(rf"^{name}[ \t]*=", src, re.M):
                old.append(f"no line '{name} = ...' in the Settings cell")
        if not re.search(r"^SOIL_CONSTANTS[ \t]*=.*'As':\s*ANISOTROPY.*'Au':\s*ANISOTROPY", src, re.M):
            old.append("SOIL_CONSTANTS does not use ANISOTROPY")
        all_code = "\n".join(c.source for c in nb.cells if c.cell_type == "code")
        if "THETAS_CAP" not in all_code or "theta*_s" not in all_code:
            old.append("the thetaS cap or the theta*_s check is missing")
        if old:
            add("FAIL", f"{args.notebook} is not the notebook this script expects ({'; '.join(old)}).")
        else:
            for j in jobs:
                patched_source(src, j)
            add("OK", f"build notebook: Settings cell (cell {idx}) has every line this script patches, and it patches for all {len(jobs)} planned runs")
    except ImportError:
        add("FAIL", "nbformat is not installed: pip install nbformat nbclient ipykernel")
    except Exception as exc:
        add("FAIL", f"Settings cell check: {exc}")
    # the reference .in: the keywords the comparison needs, and the radar folder it points to
    ref = lab / f"{BASELINE_ID}.in"
    try:
        got = read_in_file(ref)
        need = ROUTING_KWS + ["OPTPERCOLATION", "OPTSNOW", "CHANNELCONDUCTIVITY", "RUNTIME", "STARTDATE", "GWATERFILE",
                              "BEDROCKFILE", "OUTFILENAME", "RAINFILE", "RAINEXTENSION", "HYDROMETSTATIONS"]
        miss = [k for k in need if k not in got]
        add("FAIL" if miss else "OK", f"reference .in ({ref.name}) has the keywords the comparison needs" + (f" (missing: {miss})" if miss else ""))
        if not miss:
            if not same_path(got["HYDROMETSTATIONS"], LIVE_MET_SDF):
                add("WARN", f"the baseline .in points at {got['HYDROMETSTATIONS']!r}, not {LIVE_MET_SDF!r}; the weather probes are compared with {LIVE_MET_SDF}")
            folder, base, ext = rain_source(ref)
            names = [f"{base}{lab_:%m%d%Y%H}.{ext}" for lab_ in window_labels()]
            have = set(os.listdir(folder)) if os.path.isdir(folder) else set()
            lacking = [n for n in names if n not in have]
            add("FAIL" if lacking else "OK", f"radar folder the baseline reads: {folder}; " +
                (f"{len(lacking)} of the {len(names)} W2 hours are missing, e.g. {lacking[:2]}" if lacking else f"all {len(names)} W2 hours are there (only read, never changed)"))
    except Exception as exc:
        add("FAIL", f"reference .in: {exc}")
    # the weather notebook: the anchors the script needs, and the cache
    if live_ok and (lab / args.weather_notebook).is_file():
        try:
            live = weather_live_settings(lab, args.weather_notebook)
            shift = live.get("LABEL_SHIFT_HOURS")
            add("OK" if shift == BASE_SHIFT else "WARN",
                f"weather notebook Settings: LABEL_SHIFT_HOURS = {shift!r}, STATIONS = {live.get('STATIONS')}, window {live.get('LOCAL_START')} to "
                f"{live.get('LOCAL_END')}" + ("" if shift == BASE_SHIFT else
                f". The baseline is meant to be {BASE_SHIFT}; this script sets {BASE_SHIFT} itself for the baseline test, and that test (it compares "
                f"with the live files byte for byte) decides whether the live weather files really are a {BASE_SHIFT} build"))
            n_w = 0
            for p in [q for q in PROBES if q["kind"] == "weather"]:
                build_weather_nb(lab, args.weather_notebook, weather_spec(p, args, live), "full")
                n_w += 1
            build_weather_nb(lab, args.weather_notebook, weather_spec(None, args, live), "full")
            add("OK", f"weather notebook patches in memory for the baseline test and all {n_w} weather probes (Settings lines, download guard, "
                      f"temperature cell); the notebook file itself is never changed")
            hi = max(range(len(live_rows)), key=lambda i: float(live_rows[i]["elev"]))
            add("OK", f"highest station for the temperature probes: NLDAS_{hi + 1} at {live_rows[hi]['elev']} m (of {len(live_rows)}); "
                      f"one-station probe keeps the {args.one_station!r} cell")
            if not args.skip_weather_check:
                try:
                    t0 = time.time()
                    run_weather_nb(lab, args, weather_spec(None, args, live), "check", lab / LOG_DIR)
                    add("OK", f"NLDAS cache is complete for the weather notebook's stations (0 pieces to download; check took {fmt_min(time.time() - t0)}); "
                              f"no downloads and no Earthdata login are ever needed")
                except RunError as exc:
                    add("FAIL", f"NLDAS cache check: {exc}")
            else:
                add("WARN", "NLDAS cache check skipped (--skip-weather-check); the download guard still stops any weather run that would need a download")
        except Exception as exc:
            add("FAIL", f"weather notebook check: {type(exc).__name__}: {exc}")
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
    base_bytes = sum(p.stat().st_size for p in results_dir.glob(BASELINE_ID + "*") if p.is_file()) if results_dir.is_dir() else 0
    n_new = sum(1 for j in jobs if j["kind"] == "run" and (man["runs"].get(j["run_id"]) or {}).get("status") != "done")
    need = base_bytes * n_new / 1e9 + 0.2        # plus the scaled radar copies (about 14 MB each) and the weather folders (a few MB each)
    add("FAIL" if (base_bytes and need > free - args.min_free_gb) else "OK",
        f"disk: one run's files are about {base_bytes / 1e6:.0f} MB (the baseline's); {n_new} runs need about {need:.1f} GB; {free:.1f} GB free")
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
    rid, s, p = job["run_id"], settings_for(job), job["probe"]
    rec = man["runs"].setdefault(rid, {})
    log_dir, score_dir = lab / LOG_DIR, lab / SCORE_DIR
    log_dir.mkdir(exist_ok=True)
    score_dir.mkdir(exist_ok=True)
    cur = settings_record(job, args)

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
    if p["kind"] == "weather" and WEATHER_BLOCKED:
        rec.clear()
        rec.update(status="blocked", error="weather route test failed earlier in this batch: " + WEATHER_BLOCKED, finished=now_iso(), kind="probe",
                   probe=p["key"], settings=cur)
        save_manifest(man)
        log(f"{rid}: not run, because the weather route test failed earlier in this batch (see above).")
        return

    free = disk_free_gb(lab)
    if free < args.min_free_gb:
        raise Abort(f"Only {free:.1f} GB free; need {args.min_free_gb:g} GB. Free some space and run again to resume.")
    comp = competing_tribs()
    if comp:
        raise Abort("Another tRIBS is running, so the batch stops (two at once is not allowed): " + comp[0])

    rec.clear()
    rec.update(status="preparing", started=now_iso(), kind="probe", probe=p["key"], probe_kind=p["kind"], settings=cur)
    save_manifest(man)
    remove_stale(lab, rid)
    log(f"=== {rid}: {p['label']} on {WIN['name']} ({WIN['hours']} h from {WIN['start']}, water table {WIN['gw']:g} x bedrock) ===")

    # 1. make the changed forcing files and check them against the baseline's
    t_prep = time.time()
    if p["kind"] == "rain":
        prep = prep_rain(p, ref_in)
    elif p["kind"] == "weather":
        prep = prep_weather(p, args, lab, man, log_dir)
    else:
        prep = ("snow module off: OPTSNOW = 0 in the Settings cell. The weather and radar files are the baseline's; the snow-parameter "
                "file is simply not written.")
    rec.update(input_change=prep, prep_min=round((time.time() - t_prep) / 60, 2))
    save_manifest(man)
    log(f"{rid}: {prep}")

    # 2. build the model with the notebook and check what it wrote
    rec["status"] = "building"
    save_manifest(man)
    t_build = time.time()
    text, warns = build_run(job, args, lab, log_dir)
    notes = check_in_file(lab / s["INPUT_FILE"], ref_in, job)
    grid_note = check_grids(lab, job, t_build)
    rec.update(build_min=round((time.time() - t_build) / 60, 2), warnings=warns + notes, grid_check=grid_note, status="built")
    save_manifest(man)
    log(f"{rid}: build OK in {fmt_min(time.time() - t_build)}; the .in file matches the reference except the allowed lines")
    for n in notes:
        log(f"{rid}: {n}")
    log(f"{rid}: soil unchanged: {grid_note}")
    for wmsg in warns:
        log(f"{rid}: notebook said: {wmsg}")
    if args.build_only:
        log(f"{rid}: --build-only, not running tRIBS")
        return

    # 3. run tRIBS and score
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
        log(f"{rid}: NOTE: the outlet file is byte-identical to the baseline's. tRIBS ran, but this change did not reach the outlet flow "
            f"(either it has no pathway to the outlet in this window, or tRIBS does not read the changed file).")

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



# --------------------------------------------------------------------------------------------
# The table and the summary CSV
# --------------------------------------------------------------------------------------------
def print_summary(man: dict, lab: Path) -> None:
    scores = load_scores(lab)
    base = scores.get(BASELINE_ID)
    out_rows = []
    log()
    log("FORCING PROBES ON W2 (Jan 18 - Mar 12 2017). Everything not probed is at the baseline.")
    log("  vol = model volume / gauge volume.  pk = model peak vs gauge peak, %.  arr = first-arrival error, hours "
        "(model minus gauge; negative = early).")
    log("  tail1, tail2 = model / gauge volume on the Jan 26-Feb 2 and Mar 3-11 recession tails (proposed tolerance: 0.7 or more).")
    log("  bal = runoff ratio of Feb 27-28 over Feb 18-19 (the gauge's is 4.57; the design doc's 3.4 is the tail-window version).  span = whole-window volume ratio.  "
        "kge = flood-week KGE_2012.")
    log("  chg = the probe minus the baseline row (volume ratio points, peak percentage points, hours).")
    log("  Proposed tolerances: volume within 10%, peak within 15%, arrival within 1 h.")
    if base is None:
        log(f"  (the baseline run {BASELINE_ID} has no score yet, so the chg columns are blank)")
    rows_info = [dict(key="baseline", label="baseline", rid=BASELINE_ID, is_base=True)]
    rows_info += [dict(key=p["key"], label=p["label"], rid=run_id_for(p), is_base=False) for p in PROBES]
    log()
    log("    " + f"{'probe':<10}| " + " | ".join(f"{t:<21}" for _, _, t in EVENTS))
    log("    " + f"{'':<10}| " + " | ".join(f"{'vol':>6} {'pk %':>7} {'arr h':>6}" for _ in EVENTS))
    parsed = []
    for r in rows_info:
        sc = scores.get(r["rid"])
        rec = man["runs"].get(r["rid"]) or {}
        status = rec.get("status", "not run")
        vals = {n: fnum(sc, c) for n, c in METRICS} if sc else None
        if vals is None:
            log("    " + f"{r['key']:<10}| {'(' + ('baseline: ' if r['is_base'] else '') + status + ')':^70}")
        else:
            groups = [f"{fm(vals[s_ + '_vol'], 3):>6} {fm(vals[s_ + '_pk'], 1):>7} {fm(vals[s_ + '_arr']):>6}" for s_, _, _ in EVENTS]
            log("    " + f"{r['key']:<10}| " + " | ".join(groups))
        parsed.append((r, sc, rec, status, vals))
    log()
    log("    " + f"{'probe':<10}| {'tail1':>6} {'tail2':>6} {'bal':>6} {'span':>6} {'kge':>6} | {'chg fw vol':>10} {'chg fw pk':>9} "
                 f"{'chg fw arr':>10} {'chg span':>8} | note")
    for r, sc, rec, status, vals in parsed:
        row = dict(probe=r["key"], description=r["label"], run_id=r["rid"], is_baseline=r["is_base"], status=status, kind=rec.get("probe_kind", ""))
        if vals is None:
            log("    " + f"{r['key']:<10}| {'(' + ('baseline: ' if r['is_base'] else '') + status + ')':^60}")
            out_rows.append(row)
            continue
        d = {n: (vals[n] - fnum(base, c) if (base is not None and not r["is_base"]) else float("nan")) for n, c in METRICS}
        notes = []
        if r["is_base"]:
            notes.append("baseline (reused)")
        if rec.get("identical_to_baseline"):
            notes.append("IDENTICAL to baseline")
        bad = [lbl for lbl, c in zip(("fw", "s2", "f18", "f27"), FLAG_COLS) if sc.get(c) not in ("ok", "", None)]
        if bad:
            notes.append("arrival n/a: " + ",".join(bad))
        log("    " + f"{r['key']:<10}| {fm(vals['tail1']):>6} {fm(vals['tail2']):>6} {fm(vals['bal']):>6} {fm(vals['span'], 3):>6} "
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
        row["input_change"] = rec.get("input_change", "")
        row["build_min"], row["tribs_min"], row["output_mb"] = rec.get("build_min", ""), rec.get("tribs_min", ""), rec.get("output_mb", "")
        out_rows.append(row)
    log()
    log("  What each probe changed (measured by this script before the run):")
    for r, sc, rec, status, vals in parsed:
        if rec.get("input_change"):
            log(f"    {r['key']:<8} {rec['input_change']}")
    blocked = [r["key"] for r, _, rec, status, _ in parsed if status == "blocked"]
    if blocked:
        log()
        log("  NOT run because the weather route test failed: " + ", ".join(blocked))
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
        log(f"  Wrote {SUMMARY_CSV}. These tables do not say which probes matter: that is the reading step, done with the written rule.")


# --------------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------------
def parse_args(argv=None):
    ap = argparse.ArgumentParser(description="Stage 2 piece 3: the forcing probes on W2 (see the top of this file).")
    ap.add_argument("--dry-run", action="store_true", help="check everything, list the plan, build nothing")
    ap.add_argument("--summary-only", action="store_true", help="print the table from the runs already scored, then stop")
    ap.add_argument("--build-only", action="store_true", help="make the forcing files and build the inputs, but do not run tRIBS")
    ap.add_argument("--limit", type=int, default=None, help="make only the first N probes that still need work (finished ones do not count)")
    ap.add_argument("--only", default=None, help="comma-separated probes: " + ",".join(p["key"] for p in PROBES) + " (default: all eight)")
    ap.add_argument("--one-station", default=DEFAULT_ONE_STATION,
                    help=f"NLDAS cell kept for the one-station probe (default {DEFAULT_ONE_STATION}); must be one of the weather notebook's STATIONS")
    ap.add_argument("--temp-keeps", choices=["rh", "q"], default="rh",
                    help="temperature probes: 'rh' (default) keeps relative humidity as it was, so only air temperature changes in the files; "
                         "'q' keeps specific humidity, so relative humidity changes with temperature (and can pass 100%%)")
    ap.add_argument("--redo", default=None, help="comma-separated run ids to run again even if done (you almost never want this)")
    ap.add_argument("--timeout-min", type=float, default=None, help="stop a tRIBS run after this many minutes (default: about 35)")
    ap.add_argument("--cell-timeout", type=int, default=1800, help="seconds allowed for one notebook cell")
    ap.add_argument("--min-free-gb", type=float, default=2.0)
    ap.add_argument("--heartbeat-min", type=float, default=HEARTBEAT_MIN,
                    help="print a progress line every this many minutes while tRIBS runs (default 5; keeps a Codespace awake)")
    ap.add_argument("--notebook", default=None, help="build notebook (default: Make_Sycamore_Model.ipynb, else Make_Sycamore_Model_clean.ipynb)")
    ap.add_argument("--weather-notebook", default=WEATHER_NOTEBOOK)
    ap.add_argument("--scorer", default=SCORER)
    ap.add_argument("--tribs", default=os.environ.get("TRIBS_BIN") or (TRIBS_DEFAULT if os.path.isfile(TRIBS_DEFAULT)
                                                                       else shutil.which("tRIBS") or TRIBS_DEFAULT))
    ap.add_argument("--skip-kernel-check", action="store_true", help="skip the pytRIBS import test in --dry-run")
    ap.add_argument("--skip-weather-check", action="store_true", help="skip the NLDAS cache check (about a minute) in --dry-run and at the start")
    return ap.parse_args(argv)


def main(argv=None) -> int:
    global _LOGFILE
    args = parse_args(argv)
    lab = Path(__file__).resolve().parent
    os.chdir(lab)
    (lab / LOG_DIR).mkdir(exist_ok=True)
    _LOGFILE = str(lab / LOG_DIR / "batch_log.txt")
    keys = [p["key"] for p in PROBES]
    want = [x.strip() for x in (args.only or "").split(",") if x.strip()]
    unknown = [x for x in want if x not in keys]
    if unknown:
        print(f"Unknown probe(s) {unknown}. Choose from {keys}")
        return 2
    sel = [k for k in keys if (not want or k in want)]
    man = load_manifest()

    if args.summary_only:
        print_summary(man, lab)
        return 0

    args.notebook, nb_note = pick_notebook(lab, args.notebook)
    args.redo_ids = set(x.strip() for x in (args.redo or "").split(",") if x.strip())

    base_problems, base_notes = reuse_check(lab, BASELINE_ID, BASE["KS_MULT"], BASE["PSIB_MULT"])
    preview = [make_job(probe_by_key(k)) for k in sel]

    def done_already(j):
        rec = man["runs"].get(j["run_id"], {})
        return rec.get("status") == "done" and same_record(rec.get("settings"), settings_record(j, args))

    n_done = sum(1 for j in preview if done_already(j))
    n_run = sum(1 for j in preview if not done_already(j))
    log(f"Working folder: {lab}")
    log(f"Notebook: {args.notebook}; weather notebook: {args.weather_notebook}")
    log(f"Plan: {len(preview)} probe(s) on window {WIN['name']} ({WIN['hours']} h): {n_run} run(s) to make"
        + (f"; {n_done} already finished, skipped." if n_done else "."))
    log(f"Reference run: existing {BASELINE_ID} ({'; '.join(base_notes) if not base_problems else 'CANNOT be used: ' + '; '.join(base_problems)})")
    for j in preview:
        what = "already done, skipped" if done_already(j) else "NEW RUN"
        log(f"  {j['run_id']:44s} {describe(j):42s} {what}")
    log(f"  tRIBS time alone: about {n_run * WIN['hours'] * MIN_PER_1000H / 1000:.0f}-{n_run * 13:.0f} min "
        f"(W2 runs took 10-13 min); each run also needs about 1 min to build and, for a weather probe, a few minutes to rerun the weather notebook.")

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

    man["script"] = "run_screen_forcing_sycamore.py"
    man["baseline_id"] = BASELINE_ID
    man["last_started"] = now_iso()
    save_manifest(man)
    install_term_handler()

    t_batch, code = time.time(), 0
    try:
        base_job = make_job(None)
        if not do_reuse(base_job, args, lab, man, None):
            log("The W2 baseline cannot be the reference run, so nothing was run. Nothing is rerun automatically: "
                "tell me what the messages above say and we decide together.")
            return 1
        ref_in = lab / f"{BASELINE_ID}.in"
        todo = preview
        if args.limit is not None:
            fresh = [j for j in todo if not (done_already(j) and j["run_id"] not in args.redo_ids)][:args.limit]
            todo = fresh
            if not fresh:
                log("Nothing left to run: every probe in this selection is already done.")
        consecutive_fail = 0
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
                code = 1
                if consecutive_fail >= 2:
                    log("Two runs failed in a row, so the batch stops (the cause is probably the same).")
                    break
    except Abort as exc:
        log(f"BATCH STOPPED: {exc}")
        code = 1
    except KeyboardInterrupt:
        log("BATCH INTERRUPTED. Run the same command again to resume; finished runs are skipped.")
        for rid, rec in man["runs"].items():
            if rec.get("status") in ("preparing", "building", "running", "checking"):
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
    print_summary(man, lab)
    return code


if __name__ == "__main__":
    sys.exit(main())
