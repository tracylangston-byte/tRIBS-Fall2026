#!/usr/bin/env bash
# run_wy2020_screens_s11.sh   (v1.0, 2026-11-09)
# =============================================================================
# Builds the WY2020 model inputs and runs tRIBS once for the screen you name. The default is S11:
#     S11  S1 with the air-entry pressure (PsiB) multiplier 1.0 (SOLUS values unchanged; every earlier run: 0.2939)
#          setting S11-Ks1.0-PsiB1.0   results Sycamore_WY2020_S11_Ks1p0_PsiB1p0
# S1 itself = Ks x1.0, f x0.15, As 300, Au 1, bedrock SOLUS + 1 m, starting water table 0.95 x bedrock depth, PsiB x0.2939, channel loss off, snow on.
# S11 differs from S1 in exactly one value. It starts Sep 1 2019 and runs 6,575 hours (to May 31 2020 23:00).
#
# This is a NEW script, based on run_wy2020_screens_s9_s10.sh v1.0. It does not change that script, or any other.
# It needs the notebook made by apply_wy2020_s11_patch.py: Make_Sycamore_Model_WY2020_screens_S11.ipynb.
#
# For the screen it:
#   a. checks that the notebook's Settings cell gives the values listed above (Ks, f, PsiB, As, Au, bedrock, water table,
#      channel loss) and stops if anything differs, before spending any time
#   b. runs the notebook from top to bottom (SYCAMORE_SETTING tells it which screen), which writes lab/Sycamore_WY2020_<tag>.in
#   c. checks the keywords in that input file (start date, run length, results name, channel loss off)
#   d. NEW in this script: compares the new input file with S1's (Sycamore_WY2020_S1_Ks1p0.in). They must be identical apart from the
#      results name (the PsiB grid is not in the input file, it is in the soil folder, which step e checks)
#   e. NEW in this script: checks the soil grids the notebook just wrote in data/model/soil: Ks and PsiB must be plain, unscaled copies of
#      the SOLUS grids (byte for byte), and f must be the SOLUS grid x0.15. This is the proof that PsiB x1.0 reached the model.
#   f. runs tRIBS on that input file, in this terminal, and prints a progress line every 5 minutes
#   g. checks that the outlet file reached the end of the run (row count and last hour only)
#
# It does NOT open, score, plot or compare any flow. It prints no flow values.
#
# HOW TO RUN (from workspaces/Sycamore/lab/):
#     chmod +x run_wy2020_screens_s11.sh                     # once
#     CHECK_ONLY=1 ./run_wy2020_screens_s11.sh               # builds the input file and runs the checks, does NOT run tRIBS (a dry run, a few minutes)
#     ./run_wy2020_screens_s11.sh                            # the real thing (about an hour; an estimate)
#
# Keep the terminal open and the Codespace tab open while it runs. Do not run any other notebook that writes to data/model/.
# Do not run it at the same time as a hand-started tRIBS run (it checks for one and stops).
#
# OPTIONAL SETTINGS (put them in front of the command, for example  TRIBS_BIN=/path/tRIBS ./run_wy2020_screens_s11.sh):
#     NOTEBOOK       the notebook to run          (default Make_Sycamore_Model_WY2020_screens_S11.ipynb)
#     TRIBS_BIN      the tRIBS program            (default /opt/tRIBS/build/tRIBS, else whatever "which tRIBS" finds)
#     RADAR_DIR      the folder with the radar_*.asc files (default: searched for)
#     PYTHON         python to use                (default python)
#     SETTINGS       which screen                 (default "S11")
#     CHECK_ONLY     1 = build inputs and check only
#     TIMEOUT_MIN    stop a tRIBS run after this many minutes (default 240)
#
# VERSION HISTORY
#   v1.0 (2026-11-09)  first version. Tested with stand-in programs (a fake tRIBS, a small notebook, small soil grids), not in your Codespace.
# =============================================================================
set -u -o pipefail

NOTEBOOK="${NOTEBOOK:-Make_Sycamore_Model_WY2020_screens_S11.ipynb}"
PYTHON="${PYTHON:-python}"
SETTINGS="${SETTINGS:-S11}"
CHECK_ONLY="${CHECK_ONLY:-0}"
TIMEOUT_MIN="${TIMEOUT_MIN:-240}"
HEARTBEAT_SEC="${HEARTBEAT_SEC:-300}"
RUN_HOURS="${RUN_HOURS:-6575}"
RADAR_LAST="${RADAR_LAST:-radar_0531202023.asc}"
RADAR_FIRST="${RADAR_FIRST:-radar_0901201901.asc}"
SDF="${SDF:-../init_data/met_wy2020/Master_Met.sdf}"
LOGDIR="logs"
ROWS_PER_HOUR=16          # the outlet file gets 16 rows per model hour (same as the earlier runs)

say()  { printf '\n[%s] %s\n' "$(date '+%H:%M:%S')" "$*"; }
die()  { printf '\nSTOPPED: %s\n' "$*" >&2; exit 1; }

# ---------------------------------------------------------------- 1. where am I, what is there
cd "$(dirname "${BASH_SOURCE[0]}")" || die "cannot change to the script's folder."
[ "$(basename "$PWD")" = "lab" ] || die "run this from workspaces/Sycamore/lab/ (now in $PWD)."
[ -f "$NOTEBOOK" ] || die "$NOTEBOOK is not in $PWD. Run  python apply_wy2020_s11_patch.py  first (it makes it from your S10 notebook)."
[ -f "$SDF" ] || die "weather station list $SDF not found. Run Generate_Met_Forcing_Sycamore_WY2020.ipynb first."
mkdir -p "$LOGDIR"

"$PYTHON" -c "import nbclient, nbformat, jupyter_client, ipykernel, pytRIBS" 2>/dev/null \
  || die "this python ($PYTHON) cannot import nbclient, nbformat, jupyter_client, ipykernel and pytRIBS. Try:  pip install nbclient nbformat ipykernel   (pytRIBS should already be installed)."

TRIBS="${TRIBS_BIN:-}"
if [ -z "$TRIBS" ]; then
  if [ -x /opt/tRIBS/build/tRIBS ]; then TRIBS=/opt/tRIBS/build/tRIBS; else TRIBS="$(command -v tRIBS || true)"; fi
fi
if [ "$CHECK_ONLY" != "1" ]; then
  [ -n "$TRIBS" ] && [ -x "$TRIBS" ] || die "cannot find the tRIBS program. Give it with  TRIBS_BIN=/path/to/tRIBS ./run_wy2020_screens_s11.sh"
  if pgrep -x tRIBS >/dev/null 2>&1; then
    die "a tRIBS run is already going ($(pgrep -x tRIBS | tr '\n' ' ')). Wait for it to finish; two runs must not overlap."
  fi
fi

FREE_GB=$(df -P . | awk 'NR==2 {printf "%.0f", $4/1048576}')
[ "${FREE_GB:-0}" -ge 5 ] || die "only ${FREE_GB} GB of disk free here; need at least 5."

# ---------------------------------------------------------------- 2. the radar folder
if [ -z "${RADAR_DIR:-}" ]; then
  say "Looking for the WY2020 radar folder (the one that has $RADAR_LAST). This can take a minute."
  mapfile -t HITS < <(find "$HOME" /workspaces /mnt -name "$RADAR_LAST" -type f 2>/dev/null | xargs -r -n1 dirname | sort -u)
  if [ "${#HITS[@]}" -eq 0 ]; then
    die "no folder holding $RADAR_LAST was found under $HOME, /workspaces or /mnt. The WY2020 radar rain files are not on this machine.
       If they are somewhere else, run:  RADAR_DIR=/folder/with/the/files ./run_wy2020_screens_s11.sh"
  elif [ "${#HITS[@]}" -gt 1 ]; then
    printf '  found in more than one place:\n'; printf '    %s\n' "${HITS[@]}"
    die "pick one:  RADAR_DIR=<one of the folders above> ./run_wy2020_screens_s11.sh"
  fi
  RADAR_DIR="${HITS[0]}"
fi
[ -f "$RADAR_DIR/$RADAR_FIRST" ] || die "$RADAR_DIR has no $RADAR_FIRST (the first hour of the run)."
[ -f "$RADAR_DIR/$RADAR_LAST" ]  || die "$RADAR_DIR has no $RADAR_LAST (the last hour of the run)."
NRAD=$(find "$RADAR_DIR" -maxdepth 1 -name 'radar_*.asc' | wc -l)
say "Radar folder: $RADAR_DIR  ($NRAD radar files)"
export SYCAMORE_RADAR_DIR="$RADAR_DIR"

# ---------------------------------------------------------------- helpers
check_settings_cell() {   # $1 = full setting name, $2 = tag, $3... = KEY=value pairs the setting must give
  # Runs ONLY the Settings cell of the notebook (nothing else, no files) and compares what it gives with what this screen is supposed to be.
  local setting="$1" tag="$2"; shift 2
  SYCAMORE_SETTING="$setting" "$PYTHON" - "$NOTEBOOK" "$tag" "$@" <<'PY'
import contextlib, io, json, os, sys
nb_path, tag, *pairs = sys.argv[1:]
nb = json.load(open(nb_path, encoding='utf-8'))
cells = [c for c in nb['cells'] if c.get('id') == 'b5dfbf73']
if len(cells) != 1:
    sys.exit(f'the notebook has {len(cells)} cells with the Settings cell id (needed 1), so it is not the notebook I expect')
ns = {'os': os}
try:
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(''.join(cells[0]['source']), '<settings cell>', 'exec'), ns)
except ValueError as e:
    sys.exit(f'the notebook does not have this setting ({e}). Did apply_wy2020_s11_patch.py run, and is NOTEBOOK the file it wrote?')
bad = []
for p in pairs:
    k, v = p.split('=')
    if abs(float(ns[k]) - float(v)) > 1e-12:
        bad.append(f'{k} is {ns[k]} (this screen needs {v})')
want_in, want_res = f'Sycamore_WY2020_{tag}.in', f'Sycamore_WY2020_{tag}'
if ns['INPUT_FILE'] != want_in or ns['RESULTS_NAME'] != want_res:
    bad.append(f"names are {ns['INPUT_FILE']} / {ns['RESULTS_NAME']} (needed {want_in} / {want_res})")
if bad:
    sys.exit('the Settings cell does not give the values this screen needs:\n   - ' + '\n   - '.join(bad))
print('   ' + ', '.join(f'{k}={ns[k]}' for k in ('KS_MULT', 'F_MULT', 'PSIB_MULT', 'ANISOTROPY', 'ANISOTROPY_U', 'BEDROCK_ADD_M', 'GW_FRACTION', 'OPTPERCOLATION')))
PY
}

build_inputs() {   # $1 = setting name, $2 = where to save the executed copy; builds the input file by running the notebook headlessly
  local setting="$1" out="$2"
  SYCAMORE_SETTING="$setting" "$PYTHON" - "$NOTEBOOK" "$out" "$PWD" <<'PY'
import sys, nbformat, nbclient
from jupyter_client.kernelspec import KernelSpec
from jupyter_client.manager import KernelManager
nb_in, nb_out, lab = sys.argv[1:4]
nb = nbformat.read(nb_in, as_version=4)
km = KernelManager(kernel_name="python3")          # a kernel that uses THIS python, so it sees the same packages
km._kernel_spec = KernelSpec(argv=[sys.executable, "-m", "ipykernel_launcher", "-f", "{connection_file}"],
                             display_name="wy2020", language="python")
client = nbclient.NotebookClient(nb, km=km, timeout=3600, startup_timeout=180, allow_errors=False,
                                 resources={"metadata": {"path": lab}})
try:
    client.execute()
finally:
    nbformat.write(nb, nb_out)                     # written even if a cell failed, so you can see where
    try:
        km.shutdown_kernel(now=True)
    except Exception:
        pass
PY
}

keyword_value() {  # $1 = input file, $2 = KEYWORD: the line after the "KEYWORD:" line
  awk -v k="$2:" '$1==k {getline v; print v; exit}' "$1"
}

compare_with_s1() {   # $1 = new input file, $2 = its tag. The only difference from S1's input file may be the results name.
  local infile="$1" tag="$2" s1in="Sycamore_WY2020_S1_Ks1p0.in" tmp d
  if [ ! -f "$s1in" ]; then
    say "NOT CHECKED: S1's input file $s1in is not in $PWD, so $infile was not compared with it."
    return 0
  fi
  tmp="$LOGDIR/S1_input_with_${tag}_name.tmp"
  sed "s/WY2020_S1_Ks1p0/WY2020_${tag}/g" "$s1in" > "$tmp"
  if d=$(diff "$tmp" "$infile"); then
    say "$infile is identical to S1's input file ($s1in) apart from the results name. Good: nothing else in the input file changed."
    rm -f "$tmp"
    return 0
  fi
  printf '\n%s differs from S1'"'"'s input file in more than the results name (first lines; "<" is S1 with the name swapped, ">" is the new file):\n' "$infile" >&2
  printf '%s\n' "$d" | head -n 20 >&2
  return 1
}

check_soil_grids() {   # $1 = when the build started (seconds), $2 = the executed notebook, $3... = VAR=multiplier (1 = plain copy of the SOLUS grid)
  # Looks at the soil grids the notebook just wrote in data/model/soil. Reads files only; prints no flow values.
  local t0="$1" executed="$2"; shift 2
  "$PYTHON" - "$t0" "$executed" "$@" <<'PY'
import csv, filecmp, json, os, sys
from pathlib import Path
import numpy as np

lab = Path.cwd()
MODEL_SOIL, INIT_DATA = 'data/model/soil', '../init_data'
t0, executed = float(sys.argv[1]), sys.argv[2]
want = {}
for p in sys.argv[3:]:
    k, v = p.split('=')
    want[k] = float(v)


def fail(msg):
    sys.exit('   ' + msg)


def read_gdf(path):
    if not path.is_file():
        fail(f'{path} not found.')
    with open(path) as fh:
        rd = csv.DictReader(fh)
        rows = [{(k or '').strip(): (v or '').strip() for k, v in r.items()} for r in rd]
        cols = [c.strip() for c in (rd.fieldnames or [])]
    if not rows or 'Variable' not in cols or 'BasePath' not in cols:
        fail(f'{path} does not have the columns Variable and BasePath (it has {cols}), so this check cannot read it. Nothing was run.')
    return {r['Variable']: (r['BasePath'], (r.get('FileExtension') or 'asc').lstrip('.')) for r in rows}


def read_asc(path):
    header, nhead = {}, 0
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) == 2 and parts[0].lower() in ('ncols', 'nrows', 'xllcorner', 'yllcorner', 'xllcenter', 'yllcenter', 'cellsize', 'nodata_value'):
                header[parts[0].lower()] = float(parts[1])
                nhead += 1
            else:
                break
    data = np.loadtxt(path, skiprows=nhead, dtype=float)
    return header, data.reshape(int(header['nrows']), int(header['ncols']))


model_rows = read_gdf(lab / MODEL_SOIL / 'scgrid.gdf')
src_rows = read_gdf(lab / INIT_DATA / 'scgrid.gdf')
for var, mult in want.items():
    if var not in model_rows or var not in src_rows:
        fail(f'scgrid.gdf has no row {var} (model folder: {sorted(model_rows)}, SOLUS: {sorted(src_rows)}).')
    mb, mext = model_rows[var]
    sb, sext = src_rows[var]
    mf, sf = lab / f'{mb}.{mext}', lab / f'{sb}.{sext}'
    name = os.path.basename(mb)
    if not mf.is_file():
        fail(f'soil grid {var}: {mf} not found.')
    if not sf.is_file():
        fail(f'soil grid {var}: the SOLUS grid {sf} not found.')
    if mult == 1:
        if '_x' in name:
            fail(f"scgrid.gdf row {var} is '{mb}', which has a multiplier in its name, but this screen needs x1.")
        if not filecmp.cmp(mf, sf, shallow=False):
            fail(f'soil grid {var}: {mf.name} differs from the SOLUS grid {sf.name}, but this screen needs x1 (an unchanged copy).')
        print(f'   {var} x1: {mf.name} in data/model/soil is a byte-for-byte copy of the SOLUS grid {sf.name}')
    else:
        if not name.endswith(f'_x{mult:g}'):
            fail(f"scgrid.gdf row {var} is '{mb}', which does not match the multiplier {mult:g} this screen needs.")
        if mf.stat().st_mtime < t0 - 2:
            fail(f'soil grid {var}: {mf.name} is older than this build (a stale file).')
        hm, gm = read_asc(mf)
        hs, gs = read_asc(sf)
        if gm.shape != gs.shape:
            fail(f'soil grid {var}: the scaled grid and the SOLUS grid have different shapes.')
        nd_s, nd_m = hs.get('nodata_value', -9999), hm.get('nodata_value', -9999)
        ok = (gs != nd_s) & (gm != nd_m) & (gs != 0)
        if var == 'FD':
            ok &= gs * mult > 1e-4 * 1.02            # cells at the f floor do not scale
        if int(ok.sum()) < 10:
            fail(f'soil grid {var}: fewer than 10 cells left to check the multiplier on.')
        ratio = gm[ok] / gs[ok]
        if abs(float(np.median(ratio)) - mult) > 0.005 * mult or float(np.max(np.abs(ratio - mult))) > 0.02 * mult:
            fail(f'soil grid {var}: scaled / SOLUS ratio has median {np.median(ratio):.5f}, range {ratio.min():.5f} to {ratio.max():.5f} (needed {mult:g}).')
        print(f'   {var} x{mult:g}: {mf.name} is the SOLUS grid x{mult:g} (median ratio {np.median(ratio):.4f} over {int(ok.sum())} cells)')

# The notebook's own printout of the soil values (basin cells only), copied from the executed notebook. Information, not a check.
try:
    nb = json.load(open(executed, encoding='utf-8'))
    seen = []
    for c in nb.get('cells', []):
        for o in c.get('outputs', []):
            txt = ''.join(o.get('text', [])) if o.get('output_type') == 'stream' else ''
            for line in txt.splitlines():
                if 'inside the basin' in line or 'psib_0-5cm' in line or 'm_0-5cm' in line:
                    seen.append(line.strip())
    if seen:
        print("   the notebook's own printout of soil values (basin cells):")
        for line in seen[:12]:
            print('     ' + line)
except Exception as e:                                   # information only
    print(f'   (could not read the executed notebook for the soil printout: {e})')
PY
}

TRIBS_PID=""
cleanup() { if [ -n "$TRIBS_PID" ] && kill -0 "$TRIBS_PID" 2>/dev/null; then say "Stopping tRIBS (pid $TRIBS_PID)."; kill -TERM -- "-$TRIBS_PID" 2>/dev/null || kill -TERM "$TRIBS_PID" 2>/dev/null; fi; }
trap 'cleanup; die "interrupted."' INT TERM

run_tribs() {      # $1 = input file, $2 = log file, $3 = outlet file
  local infile="$1" logf="$2" outlet="$3" t0 now elapsed rows expected
  expected=$(( RUN_HOURS * ROWS_PER_HOUR ))
  t0=$(date +%s)
  setsid "$TRIBS" "$infile" > "$logf" 2>&1 &
  TRIBS_PID=$!
  while kill -0 "$TRIBS_PID" 2>/dev/null; do
    sleep "$HEARTBEAT_SEC" &
    wait $! 2>/dev/null
    kill -0 "$TRIBS_PID" 2>/dev/null || break
    now=$(date +%s); elapsed=$(( (now - t0) / 60 ))
    rows=0; [ -f "$outlet" ] && rows=$(wc -l < "$outlet")
    say "tRIBS still running: ${elapsed} min so far; outlet file has ${rows} rows of about ${expected} ($(( 100 * rows / expected ))%)"
    if [ "$elapsed" -ge "$TIMEOUT_MIN" ]; then
      cleanup; TRIBS_PID=""; die "tRIBS did not finish within ${TIMEOUT_MIN} minutes and was stopped."
    fi
  done
  wait "$TRIBS_PID"; local rc=$?
  TRIBS_PID=""
  return $rc
}

check_outlet() {   # counts only: rows, last hour, non-numbers. Prints no flow values.
  local outlet="$1" rows last bad expected
  [ -f "$outlet" ] || die "tRIBS did not write $outlet."
  expected=$(( RUN_HOURS * ROWS_PER_HOUR ))
  rows=$(awk -F, 'NR>1 && NF>=2 {n++} END {print n+0}' "$outlet")
  last=$(awk -F, 'NR>1 && NF>=2 {t=$1} END {print t+0}' "$outlet")
  bad=$(awk -F, 'NR>1 && NF>=2 && ($2 ~ /[nN][aA][nN]|[iI][nN][fF]/) {n++} END {print n+0}' "$outlet")
  say "Outlet file check: $rows rows (about $expected expected), last model hour $last (run is $RUN_HOURS), rows with NaN/Inf: $bad"
  [ "$rows" -ge $(( expected * 99 / 100 )) ] || die "$outlet has only $rows rows; expected about $expected."
  [ "$bad" -eq 0 ] || die "$outlet has $bad rows with NaN or Inf."
  awk -v l="$last" -v r="$RUN_HOURS" 'BEGIN {exit !(l >= r - 2)}' || die "$outlet stops at hour $last; the run should reach hour $RUN_HOURS."
}

# ---------------------------------------------------------------- 3. the screens, one after the other
# First make sure every name in SETTINGS is known and that none of the results files exist yet, so nothing is wasted late.
declare -A FULL TAGS PAIRS GRIDS
FULL[S11]="S11-Ks1.0-PsiB1.0";    TAGS[S11]="S11_Ks1p0_PsiB1p0"
# The screen is S1 with one value changed: Ks x1.0, f x0.15, PsiB x1.0 (S1: x0.2939), As 300, Au 1, bedrock + 1 m, water table 0.95, channel loss off.
PAIRS[S11]="KS_MULT=1 F_MULT=0.15 PSIB_MULT=1 ANISOTROPY=300 ANISOTROPY_U=1 BEDROCK_ADD_M=1 GW_FRACTION=0.95 OPTPERCOLATION=0"
# The soil grids the notebook must have written (variable=multiplier): 1 means a plain copy of the SOLUS grid.
GRIDS[S11]="KS=1 PB=1 FD=0.15"

for S in $SETTINGS; do
  [ -n "${FULL[$S]:-}" ] || die "unknown screen '$S' (use S11)."
  if [ "$CHECK_ONLY" != "1" ] && [ -e "results/Sycamore_WY2020_${TAGS[$S]}_Outlet.qout" ]; then
    die "results/Sycamore_WY2020_${TAGS[$S]}_Outlet.qout already exists. Move or delete the old $S results first (this script will not overwrite them). Nothing has been rebuilt yet."
  fi
done

START_ALL=$(date +%s)
SUMMARY=()
for S in $SETTINGS; do
  SETTING="${FULL[$S]}"; TAG="${TAGS[$S]}"
  INFILE="Sycamore_WY2020_${TAG}.in"

  say "=== $S ($SETTING): checking the Settings cell ==="
  # shellcheck disable=SC2086
  check_settings_cell "$SETTING" "$TAG" ${PAIRS[$S]} || die "the notebook $NOTEBOOK is not set up for $S (see the lines above). Nothing was run."

  say "=== $S: building the input file ($NOTEBOOK, headless; several minutes) ==="
  rm -f "$INFILE"
  T_BUILD=$(date +%s)
  build_inputs "$SETTING" "$LOGDIR/executed_Make_Sycamore_Model_WY2020_${TAG}.ipynb" > "$LOGDIR/notebook_${TAG}.log" 2>&1 \
    || die "the notebook failed for $S. Open $LOGDIR/executed_Make_Sycamore_Model_WY2020_${TAG}.ipynb to see the cell that stopped (details in $LOGDIR/notebook_${TAG}.log)."
  [ -f "$INFILE" ] || die "the notebook finished but $INFILE was not written."

  RT=$(keyword_value "$INFILE" RUNTIME); SD=$(keyword_value "$INFILE" STARTDATE); OUTP=$(keyword_value "$INFILE" OUTFILENAME)
  say "$INFILE says: STARTDATE $SD, RUNTIME $RT, OUTFILENAME $OUTP"
  [ "$RT" = "$RUN_HOURS" ] || die "$INFILE has RUNTIME $RT; expected $RUN_HOURS."
  [ "$SD" = "09/01/2019/00/00" ] || die "$INFILE has STARTDATE $SD; expected 09/01/2019/00/00."
  case "$OUTP" in *"WY2020_${TAG}") ;; *) die "$INFILE writes results to '$OUTP'; expected a name ending in WY2020_${TAG}." ;; esac
  PERC=$(keyword_value "$INFILE" OPTPERCOLATION)
  say "$INFILE says: OPTPERCOLATION $PERC (channel loss $([ "$PERC" = 1 ] && echo ON || echo off))"
  [ "$PERC" = "0" ] || die "$INFILE has OPTPERCOLATION $PERC; $S needs 0 (channel loss off)."
  OUTLET="${OUTP}_Outlet.qout"

  say "=== $S: comparing $INFILE with S1's input file ==="
  compare_with_s1 "$INFILE" "$TAG" || die "$INFILE is not S1's input file plus a new results name. tRIBS was not run."
  say "=== $S: checking the soil grids the notebook wrote in data/model/soil ==="
  # shellcheck disable=SC2086
  check_soil_grids "$T_BUILD" "$LOGDIR/executed_Make_Sycamore_Model_WY2020_${TAG}.ipynb" ${GRIDS[$S]} \
    || die "the soil grids in data/model/soil are not what $S needs (see the lines above). tRIBS was not run."

  if [ "$CHECK_ONLY" = "1" ]; then
    say "CHECK_ONLY: $S input file is built and checked. Not running tRIBS."
    SUMMARY+=("$S: input file built ($INFILE), tRIBS not run")
    continue
  fi
  [ ! -e "$OUTLET" ] || die "$OUTLET already exists. Move or delete the old $S results first (this script will not overwrite them)."

  say "=== $S: running tRIBS (log: $LOGDIR/simulation_WY2020_${TAG}.log) ==="
  T1=$(date +%s)
  run_tribs "$INFILE" "$LOGDIR/simulation_WY2020_${TAG}.log" "$OUTLET" || die "tRIBS exited with an error for $S; see $LOGDIR/simulation_WY2020_${TAG}.log (last lines: $(tail -n 3 "$LOGDIR/simulation_WY2020_${TAG}.log" | tr '\n' ' '))."
  check_outlet "$OUTLET"
  SUMMARY+=("$S: finished in $(( ($(date +%s) - T1) / 60 )) min; outlet file $OUTLET")
done

say "ALL DONE in $(( ($(date +%s) - START_ALL) / 60 )) minutes."
printf '  %s\n' "${SUMMARY[@]}"
cat <<'MSG'

This run has not been opened, scored or plotted by this script.
Score it the same way as S5 (the same seven storms, same scoring and edge-check scripts), side by side with Re-tuned and S1 (its base).
MSG
