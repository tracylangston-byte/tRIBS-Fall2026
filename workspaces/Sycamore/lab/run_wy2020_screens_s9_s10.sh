#!/usr/bin/env bash
# run_wy2020_screens_s9_s10.sh   (v1.0, 2026-11-09)
# =============================================================================
# Builds the WY2020 model inputs and runs tRIBS once for each screen you list, one after the other. The default is S9, then S10:
#     S9   S1 with the bedrock SOLUS + 0.5 m   (S1: + 1 m)        setting S9-Ks1.0-Bedrock0.5   results Sycamore_WY2020_S9_Ks1p0_Bed0p5
#     S10  S1 with the starting water table 0.5 x bedrock depth   (S1: 0.95)    setting S10-Ks1.0-GW0.5   results Sycamore_WY2020_S10_Ks1p0_GW0p5
# S1 itself = Ks x1.0, f x0.15, As 300, Au 1, bedrock SOLUS + 1 m, starting water table 0.95 x bedrock depth, channel loss off, snow on.
# Each run differs from S1 in exactly one value. Both start Sep 1 2019 and run 6,575 hours (to May 31 2020 23:00).
# (Caution for reading S9: the starting water table is 0.95 x bedrock depth, so it moves up with the bedrock. S9 changes soil depth AND
#  starting water table together; S10 changes only the starting water table.)
#
# This is a NEW script, based on run_wy2020_screens.sh v1.0 (the S6 to S8 one). It does not change that script, or run_wy2020_both.sh.
# It needs the notebook made by apply_wy2020_s9_s10_patch.py: Make_Sycamore_Model_WY2020_screens_S10.ipynb.
#
# For each screen it:
#   a. checks that the notebook's Settings cell gives the values listed above for that screen (Ks, f, As, Au, bedrock, water table,
#      channel loss) and stops if anything differs, before spending any time
#   b. runs the notebook from top to bottom (SYCAMORE_SETTING tells it which screen), which writes lab/Sycamore_WY2020_<tag>.in
#   c. checks the keywords in that input file (start date, run length, results name, channel loss off)
#   d. runs tRIBS on that input file, in this terminal, and prints a progress line every 5 minutes
#   e. checks that the outlet file reached the end of the run (row count and last hour only)
#
# It runs them one at a time, never together, because every screen writes its soil files into data/model/.
# It does NOT open, score, plot or compare any flow. It prints no flow values.
# If one screen fails, the script stops there and does not start the next one.
#
# HOW TO RUN (from workspaces/Sycamore/lab/):
#     chmod +x run_wy2020_screens_s9_s10.sh                  # once
#     CHECK_ONLY=1 ./run_wy2020_screens_s9_s10.sh            # builds the two input files, does NOT run tRIBS (a dry run, a few minutes each)
#     ./run_wy2020_screens_s9_s10.sh                         # the real thing: S9, then S10 (about an hour each; an estimate)
#     SETTINGS="S10" ./run_wy2020_screens_s9_s10.sh          # just one (also S9)
#
# Keep the terminal open and the Codespace tab open while it runs. Do not run any other notebook that writes to data/model/.
# Do not run it at the same time as a hand-started tRIBS run (it checks for one and stops).
#
# OPTIONAL SETTINGS (put them in front of the command, for example  TRIBS_BIN=/path/tRIBS ./run_wy2020_screens_s9_s10.sh):
#     NOTEBOOK       the notebook to run          (default Make_Sycamore_Model_WY2020_screens_S10.ipynb)
#     TRIBS_BIN      the tRIBS program            (default /opt/tRIBS/build/tRIBS, else whatever "which tRIBS" finds)
#     RADAR_DIR      the folder with the radar_*.asc files (default: searched for)
#     PYTHON         python to use                (default python)
#     SETTINGS       which screens, in order      (default "S9 S10")
#     CHECK_ONLY     1 = build inputs only
#     TIMEOUT_MIN    stop a tRIBS run after this many minutes (default 240)
#
# VERSION HISTORY
#   v1.0 (2026-11-09)  first version. Tested with stand-in programs (a fake tRIBS, a small notebook), not in your Codespace.
# =============================================================================
set -u -o pipefail

NOTEBOOK="${NOTEBOOK:-Make_Sycamore_Model_WY2020_screens_S10.ipynb}"
PYTHON="${PYTHON:-python}"
SETTINGS="${SETTINGS:-S9 S10}"
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
[ -f "$NOTEBOOK" ] || die "$NOTEBOOK is not in $PWD. Run  python apply_wy2020_s9_s10_patch.py  first (it makes it from your S8 notebook)."
[ -f "$SDF" ] || die "weather station list $SDF not found. Run Generate_Met_Forcing_Sycamore_WY2020.ipynb first."
mkdir -p "$LOGDIR"

"$PYTHON" -c "import nbclient, nbformat, jupyter_client, ipykernel, pytRIBS" 2>/dev/null \
  || die "this python ($PYTHON) cannot import nbclient, nbformat, jupyter_client, ipykernel and pytRIBS. Try:  pip install nbclient nbformat ipykernel   (pytRIBS should already be installed)."

TRIBS="${TRIBS_BIN:-}"
if [ -z "$TRIBS" ]; then
  if [ -x /opt/tRIBS/build/tRIBS ]; then TRIBS=/opt/tRIBS/build/tRIBS; else TRIBS="$(command -v tRIBS || true)"; fi
fi
if [ "$CHECK_ONLY" != "1" ]; then
  [ -n "$TRIBS" ] && [ -x "$TRIBS" ] || die "cannot find the tRIBS program. Give it with  TRIBS_BIN=/path/to/tRIBS ./run_wy2020_screens_s9_s10.sh"
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
       If they are somewhere else, run:  RADAR_DIR=/folder/with/the/files ./run_wy2020_screens_s9_s10.sh"
  elif [ "${#HITS[@]}" -gt 1 ]; then
    printf '  found in more than one place:\n'; printf '    %s\n' "${HITS[@]}"
    die "pick one:  RADAR_DIR=<one of the folders above> ./run_wy2020_screens_s9_s10.sh"
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
    sys.exit(f'the notebook does not have this setting ({e}). Did apply_wy2020_s9_s10_patch.py run, and is NOTEBOOK the file it wrote?')
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
print('   ' + ', '.join(f'{k}={ns[k]}' for k in ('KS_MULT', 'F_MULT', 'ANISOTROPY', 'ANISOTROPY_U', 'BEDROCK_ADD_M', 'GW_FRACTION', 'OPTPERCOLATION')))
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
declare -A FULL TAGS PAIRS
FULL[S9]="S9-Ks1.0-Bedrock0.5";   TAGS[S9]="S9_Ks1p0_Bed0p5"
FULL[S10]="S10-Ks1.0-GW0.5";      TAGS[S10]="S10_Ks1p0_GW0p5"
# Each screen is S1 with one value changed: Ks x1.0, f x0.15, As 300, Au 1, bedrock + 1 m, water table 0.95, channel loss off.
PAIRS[S9]="KS_MULT=1 F_MULT=0.15 ANISOTROPY=300 ANISOTROPY_U=1 BEDROCK_ADD_M=0.5 GW_FRACTION=0.95 OPTPERCOLATION=0"
PAIRS[S10]="KS_MULT=1 F_MULT=0.15 ANISOTROPY=300 ANISOTROPY_U=1 BEDROCK_ADD_M=1 GW_FRACTION=0.5 OPTPERCOLATION=0"

for S in $SETTINGS; do
  [ -n "${FULL[$S]:-}" ] || die "unknown screen '$S' (use S9 or S10)."
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

These runs have not been opened, scored or plotted by this script.
Score them the same way as S5 (the same seven storms, same scoring and edge-check scripts), side by side with Re-tuned and S1 (their base).
MSG
