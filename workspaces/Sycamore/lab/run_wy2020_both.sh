#!/usr/bin/env bash
# run_wy2020_both.sh   (v1.1, 2026-11-04)
# =============================================================================
# Builds the WY2020 model inputs and runs tRIBS once for each setting you list, one after the other. The settings:
#     Re-tuned           the frozen setting, channel loss off                       (the default run 1; done 2026-11-04)
#     Original           what older notes call "as drawn", channel loss off          (the default run 2; done 2026-11-04)
#     Original-Kchan0.5  the Original with channel loss ON, Kchan 0.5 mm/hr          (added in v1.1)
#     Original-Kchan0.2  the Original with channel loss ON, Kchan 0.2 mm/hr          (added in v1.1)
# All start Sep 1 2019 and run 6,575 hours (to May 31 2020 23:00).
#
# For each setting it:
#   a. runs Make_Sycamore_Model_WY2020.ipynb from top to bottom (SYCAMORE_SETTING tells it which setting),
#      which writes lab/Sycamore_WY2020_<Retuned|Original>.in
#   b. runs tRIBS on that input file, in this terminal, and prints a progress line every 5 minutes
#   c. checks that the outlet file reached the end of the run (row count and last hour only)
#
# It runs them one at a time, never together, because both settings write their soil files into data/model/.
# It does NOT open, score, plot or compare any flow. WY2020 was unsealed on 2026-11-04 under the logged scoring plan; this script still opens nothing.
#
# HOW TO RUN (from workspaces/Sycamore/lab/):
#     chmod +x run_wy2020_both.sh          # once
#     CHECK_ONLY=1 ./run_wy2020_both.sh    # builds both input files, does NOT run tRIBS (a dry run, a few minutes each)
#     ./run_wy2020_both.sh                 # the real thing for Re-tuned then Original (roughly 2 hours in all; an estimate)
#     SETTINGS="Original-Kchan0.5 Original-Kchan0.2" ./run_wy2020_both.sh     # the two channel-loss runs (v1.1)
#     CHECK_ONLY=1 SETTINGS="Original-Kchan0.5 Original-Kchan0.2" ./run_wy2020_both.sh   # their dry run first
#
# Keep the terminal open and the Codespace tab open while it runs. Do not run any other notebook that writes to data/model/.
#
# OPTIONAL SETTINGS (put them in front of the command, for example  TRIBS_BIN=/path/tRIBS ./run_wy2020_both.sh):
#     TRIBS_BIN      the tRIBS program            (default /opt/tRIBS/build/tRIBS, else whatever "which tRIBS" finds)
#     RADAR_DIR      the folder with the radar_*.asc files (default: searched for)
#     PYTHON         python to use                (default python)
#     SETTINGS       which settings, in order     (default "Re-tuned Original"; also Original-Kchan0.5, Original-Kchan0.2)
#     CHECK_ONLY     1 = build inputs only
#     TIMEOUT_MIN    stop a tRIBS run after this many minutes (default 240)
# =============================================================================
set -u -o pipefail

NOTEBOOK="${NOTEBOOK:-Make_Sycamore_Model_WY2020.ipynb}"
PYTHON="${PYTHON:-python}"
SETTINGS="${SETTINGS:-Re-tuned Original}"
CHECK_ONLY="${CHECK_ONLY:-0}"
TIMEOUT_MIN="${TIMEOUT_MIN:-240}"
HEARTBEAT_SEC="${HEARTBEAT_SEC:-300}"
RUN_HOURS="${RUN_HOURS:-6575}"
RADAR_LAST="${RADAR_LAST:-radar_0531202023.asc}"
RADAR_FIRST="${RADAR_FIRST:-radar_0901201901.asc}"
SDF="${SDF:-../init_data/met_wy2020/Master_Met.sdf}"
LOGDIR="logs"
ROWS_PER_HOUR=16          # the outlet file gets 16 rows per model hour (same as the WY2017 runs)

say()  { printf '\n[%s] %s\n' "$(date '+%H:%M:%S')" "$*"; }
die()  { printf '\nSTOPPED: %s\n' "$*" >&2; exit 1; }

# ---------------------------------------------------------------- 1. where am I, what is there
cd "$(dirname "${BASH_SOURCE[0]}")" || die "cannot change to the script's folder."
[ "$(basename "$PWD")" = "lab" ] || die "run this from workspaces/Sycamore/lab/ (now in $PWD)."
[ -f "$NOTEBOOK" ] || die "$NOTEBOOK is not in $PWD."
[ -f "$SDF" ] || die "weather station list $SDF not found. Run Generate_Met_Forcing_Sycamore_WY2020.ipynb first."
mkdir -p "$LOGDIR"

"$PYTHON" -c "import nbclient, nbformat, jupyter_client, ipykernel, pytRIBS" 2>/dev/null \
  || die "this python ($PYTHON) cannot import nbclient, nbformat, jupyter_client, ipykernel and pytRIBS. Try:  pip install nbclient nbformat ipykernel   (pytRIBS should already be installed)."

TRIBS="${TRIBS_BIN:-}"
if [ -z "$TRIBS" ]; then
  if [ -x /opt/tRIBS/build/tRIBS ]; then TRIBS=/opt/tRIBS/build/tRIBS; else TRIBS="$(command -v tRIBS || true)"; fi
fi
if [ "$CHECK_ONLY" != "1" ]; then
  [ -n "$TRIBS" ] && [ -x "$TRIBS" ] || die "cannot find the tRIBS program. Give it with  TRIBS_BIN=/path/to/tRIBS ./run_wy2020_both.sh"
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
       If they are somewhere else, run:  RADAR_DIR=/folder/with/the/files ./run_wy2020_both.sh
       If they were never made here, they have to be built first (build_mrms_radar_rainfall_tribs6.py makes them); tell Claude and we will do that step."
  elif [ "${#HITS[@]}" -gt 1 ]; then
    printf '  found in more than one place:\n'; printf '    %s\n' "${HITS[@]}"
    die "pick one:  RADAR_DIR=<one of the folders above> ./run_wy2020_both.sh"
  fi
  RADAR_DIR="${HITS[0]}"
fi
[ -f "$RADAR_DIR/$RADAR_FIRST" ] || die "$RADAR_DIR has no $RADAR_FIRST (the first hour of the run)."
[ -f "$RADAR_DIR/$RADAR_LAST" ]  || die "$RADAR_DIR has no $RADAR_LAST (the last hour of the run)."
NRAD=$(find "$RADAR_DIR" -maxdepth 1 -name 'radar_*.asc' | wc -l)
say "Radar folder: $RADAR_DIR  ($NRAD radar files)"
export SYCAMORE_RADAR_DIR="$RADAR_DIR"

# ---------------------------------------------------------------- helpers
build_inputs() {   # $1 = setting name; builds the input file by running the notebook headlessly
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

check_outlet() {   # counts only: rows, last hour, non-numbers. Prints no flow values (WY2020 is sealed).
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

# ---------------------------------------------------------------- 3. the two runs, one after the other
START_ALL=$(date +%s)
SUMMARY=()
for SETTING in $SETTINGS; do
  case "$SETTING" in
    Re-tuned)          TAG=Retuned;           WANT_PERC=0; WANT_KCHAN="" ;;
    Original)          TAG=Original;          WANT_PERC=0; WANT_KCHAN="" ;;
    Original-Kchan0.5) TAG=Original_kchan0p5; WANT_PERC=1; WANT_KCHAN=0.5 ;;
    Original-Kchan0.2) TAG=Original_kchan0p2; WANT_PERC=1; WANT_KCHAN=0.2 ;;
    *) die "unknown setting '$SETTING' (use Re-tuned, Original, Original-Kchan0.5 or Original-Kchan0.2)." ;;
  esac
  INFILE="Sycamore_WY2020_${TAG}.in"
  if [ "$CHECK_ONLY" != "1" ] && [ -e "results/Sycamore_WY2020_${TAG}_Outlet.qout" ]; then
    die "results/Sycamore_WY2020_${TAG}_Outlet.qout already exists. Move or delete the old $SETTING results first (this script will not overwrite them). Nothing has been rebuilt yet."
  fi
  say "=== $SETTING: building the input file ($NOTEBOOK, headless; several minutes) ==="
  rm -f "$INFILE"
  build_inputs "$SETTING" "$LOGDIR/executed_Make_Sycamore_Model_WY2020_${TAG}.ipynb" > "$LOGDIR/notebook_${TAG}.log" 2>&1 \
    || die "the notebook failed for $SETTING. Open $LOGDIR/executed_Make_Sycamore_Model_WY2020_${TAG}.ipynb to see the cell that stopped (details in $LOGDIR/notebook_${TAG}.log)."
  [ -f "$INFILE" ] || die "the notebook finished but $INFILE was not written."

  RT=$(keyword_value "$INFILE" RUNTIME); SD=$(keyword_value "$INFILE" STARTDATE); OUTP=$(keyword_value "$INFILE" OUTFILENAME)
  say "$INFILE says: STARTDATE $SD, RUNTIME $RT, OUTFILENAME $OUTP"
  [ "$RT" = "$RUN_HOURS" ] || die "$INFILE has RUNTIME $RT; expected $RUN_HOURS."
  [ "$SD" = "09/01/2019/00/00" ] || die "$INFILE has STARTDATE $SD; expected 09/01/2019/00/00."
  case "$OUTP" in *"WY2020_${TAG}") ;; *) die "$INFILE writes results to '$OUTP'; expected a name ending in WY2020_${TAG}." ;; esac
  PERC=$(keyword_value "$INFILE" OPTPERCOLATION); KC=$(keyword_value "$INFILE" CHANNELCONDUCTIVITY)
  say "$INFILE says: OPTPERCOLATION $PERC (channel loss $([ "$PERC" = 1 ] && echo ON || echo off)), CHANNELCONDUCTIVITY $KC"
  [ "$PERC" = "$WANT_PERC" ] || die "$INFILE has OPTPERCOLATION $PERC; $SETTING needs $WANT_PERC."
  if [ "$WANT_PERC" = "1" ]; then
    awk -v a="$KC" -v b="$WANT_KCHAN" 'BEGIN {exit !(a+0 == b+0)}' || die "$INFILE has CHANNELCONDUCTIVITY $KC; $SETTING needs $WANT_KCHAN."
  fi
  OUTLET="${OUTP}_Outlet.qout"

  if [ "$CHECK_ONLY" = "1" ]; then
    say "CHECK_ONLY: $SETTING input file is built and checked. Not running tRIBS."
    SUMMARY+=("$SETTING: input file built ($INFILE), tRIBS not run")
    continue
  fi
  [ ! -e "$OUTLET" ] || die "$OUTLET already exists. Move or delete the old $SETTING results first (this script will not overwrite them)."

  say "=== $SETTING: running tRIBS (log: $LOGDIR/simulation_WY2020_${TAG}.log) ==="
  T1=$(date +%s)
  run_tribs "$INFILE" "$LOGDIR/simulation_WY2020_${TAG}.log" "$OUTLET" || die "tRIBS exited with an error for $SETTING; see $LOGDIR/simulation_WY2020_${TAG}.log (last lines: $(tail -n 3 "$LOGDIR/simulation_WY2020_${TAG}.log" | tr '\n' ' '))."
  check_outlet "$OUTLET"
  SUMMARY+=("$SETTING: finished in $(( ($(date +%s) - T1) / 60 )) min; outlet file $OUTLET")
done

say "ALL DONE in $(( ($(date +%s) - START_ALL) / 60 )) minutes."
printf '  %s\n' "${SUMMARY[@]}"
cat <<'MSG'

These runs have not been opened, scored or plotted by this script.
Score them with the logged scoring plan (the same seven storms), side by side with Re-tuned and Original.
They are comparison runs on the validation year: do not use them to pick or tune a WY2020 setting.
MSG
