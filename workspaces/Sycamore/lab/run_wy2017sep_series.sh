#!/usr/bin/env bash
# run_wy2017sep_series.sh   (v1.0, 2026-11-21)
# =============================================================================
# Builds the model inputs and runs tRIBS for TEN settings, one after the other, for water year 2017 with a September 2016 spin-up:
#     start Sep 1 2016 00:00, 6,551 hours (to May 31 2017 23:00).  September (the first 720 hours) is spin-up and is not scored.
#
#     Retuned   Re-tuned                  Ks x0.19261, f x0.15, As 300, Au 1                      results Sycamore_WY2017sep_Retuned
#     Original  Original                  Ks x0.187,   f x0.5874, As 1, Au 1                     results Sycamore_WY2017sep_Original
#     S1        S1-Ks1.0                  Ks x1.0, f x0.15, As 300                                results Sycamore_WY2017sep_S1_Ks1p0
#     S5        S5-Ks1.0-Bedrock2.0       S1 with bedrock SOLUS + 2 m                             results Sycamore_WY2017sep_S5_Ks1p0_Bed2p0
#     S6        S6-Ks0.5-Bedrock2.0       S5 with Ks x0.5                                         results Sycamore_WY2017sep_S6_Ks0p5_Bed2p0
#     S7        S7-Ks1.0-Bedrock2.0-As500 S5 with As 500                                          results Sycamore_WY2017sep_S7_Ks1p0_Bed2p0_As500
#     S8        S8-Ks1.0-Bedrock2.0-f0.02 S5 with f x0.02                                         results Sycamore_WY2017sep_S8_Ks1p0_Bed2p0_f0p02
#     S9        S9-Ks1.0-Bedrock0.5       S1 with bedrock SOLUS + 0.5 m                           results Sycamore_WY2017sep_S9_Ks1p0_Bed0p5
#     S10       S10-Ks1.0-GW0.5           S1 with starting water table 0.5 x bedrock depth        results Sycamore_WY2017sep_S10_Ks1p0_GW0p5
#     S11       S11-Ks1.0-PsiB1.0         S1 with PsiB x1.0 (every other setting: x0.2939)        results Sycamore_WY2017sep_S11_Ks1p0_PsiB1p0
# Every setting also has: bedrock SOLUS + 1 m (unless it says otherwise), starting water table 0.95 x bedrock depth, channel loss off, snow on.
# The new result names all start Sycamore_WY2017sep_, so they cannot collide with the older WY2017 results (W5_d012_k10_...) or with WY2020.
#
# This is a NEW script, based on run_wy2020_screens_s11.sh v1.0 and run_wy2020_both.sh. It does not change any other script.
# It needs the notebook made by apply_wy2017sep_patch.py:  Make_Sycamore_Model_WY2017sep.ipynb
#
# For each setting it:
#   a. checks that the notebook's Settings cell gives that setting's values (Ks, f, PsiB, As, Au, bedrock, water table, channel loss, snow,
#      run window, radar file names, weather list, run names). If anything differs it stops before spending any time.
#   b. runs the notebook from top to bottom, headless (SYCAMORE_SETTING says which setting), which writes lab/Sycamore_WY2017sep_<tag>.in
#   c. reads the printout of the executed notebook: stops on any line that starts with CHECK, MISSING or PROBLEM, and checks the As and Au
#      the notebook read back from the soil table
#   d. checks the keywords in the input file (start date, run length, results name, channel loss off, snow on, rain and weather inputs)
#   e. compares the input file with Re-tuned's: they must be identical apart from the results name, and (only when the setting changes the
#      bedrock or the starting water table) the BEDROCKFILE and GWATERFILE lines
#   f. checks the soil grids the notebook just wrote in data/model/soil: each one is either a byte-for-byte copy of the SOLUS grid (x1) or
#      the SOLUS grid times the setting's multiplier, cell by cell (f has a floor of 0.0001)
#   g. checks the bedrock and starting-water-table grids the input file points to: bedrock = SOLUS + the setting's meters, water table =
#      the setting's fraction x bedrock depth (in mm), and that both files were copied in during this build
#   h. runs tRIBS on that input file, in this terminal, with a progress line every 5 minutes
#   i. checks that the outlet file reached the end of the run (row count and last hour only)
# It runs them one at a time, never together, because every setting writes its soil files into data/model/.
# It does NOT open, score, plot or compare any flow. It prints no flow values.
# If one setting fails, the script stops there and does not start the next one.
#
# HOW TO RUN (from workspaces/Sycamore/lab/):
#     chmod +x run_wy2017sep_series.sh                       # once
#     CHECK_ONLY=1 ./run_wy2017sep_series.sh                 # a dry run: builds and checks all ten input files, does NOT run tRIBS (about 15 to 20 min; an estimate)
#     ./run_wy2017sep_series.sh                              # the real thing: all ten, one after the other (very roughly 9 to 10 hours; an estimate)
#     SETTINGS="S5 S6" ./run_wy2017sep_series.sh             # only some of them, in the order you list
#
# IF IT STOPS PART WAY (the terminal closed, the Codespace slept, a check failed): fix the cause, then run the same command again.
# Settings whose outlet file is finished are skipped. A setting that was interrupted has its leftover files moved (not deleted) to
# results/_partial_<date>/ and is run again from the start.
# Keep the terminal and the Codespace browser tab open while it runs. Do not run any other notebook that writes to data/model/.
# Do not run it at the same time as a hand-started tRIBS run (it checks for one and stops, in a dry run too).
# Everything it prints is also saved in logs/series_<date>.log.
#
# OPTIONAL SETTINGS (put them in front of the command, for example  TRIBS_BIN=/path/tRIBS ./run_wy2017sep_series.sh):
#     NOTEBOOK       the notebook to run          (default Make_Sycamore_Model_WY2017sep.ipynb)
#     TRIBS_BIN      the tRIBS program            (default /opt/tRIBS/build/tRIBS, else whatever "which tRIBS" finds)
#     RADAR_DIR      the folder with the radar_*.asc files (default ~/mrms_wk_wy2017_sep/rain)
#     PYTHON         python to use                (default python)
#     SETTINGS       which settings, in order     (default "Retuned Original S1 S5 S6 S7 S8 S9 S10 S11")
#     CHECK_ONLY     1 = build inputs and check only
#     TIMEOUT_MIN    stop a tRIBS run after this many minutes (default 240)
#     HEARTBEAT_SEC  seconds between progress lines (default 300)
#
# WHAT THE CHECKS CANNOT PROVE: they read your files, but they do not know what tRIBS does with them. They do not look at any flow, so a run
# that finishes can still be a bad run. And they trust the SOLUS grids in ../init_data (they compare against them, not against the world).
#
# RUN LENGTH: 6,551 hours = Sep 1 2016 00:00 to May 31 2017 23:00 (273 days = 6,552 hourly weather rows, less one hour, as for WY2020's 6,575 of 6,576).
# SCORING LATER: Time_hr 0 in these results is Sep 1 2016 00:00, so Oct 1 2016 00:00 is hour 720 and the scored window is hours 720 to 6551
# (5,832 hours). The older scorers assume an Oct 1 start, so they will need that 720-hour offset. This script does not score anything.
#
# VERSION HISTORY
#   v1.0 (2026-11-21)  first version. Tested with stand-in programs (a fake tRIBS, a stand-in notebook, small fake soil and bedrock grids, a fake
#                      radar folder), not in your Codespace.
# =============================================================================

# ---------------------------------------------------------------- 0. keep a log of everything (run this script again, once, through tee)
if [ -z "${WY17_INNER:-}" ]; then
  cd "$(dirname "${BASH_SOURCE[0]}")" || { echo "STOPPED: cannot change to the script's folder." >&2; exit 1; }
  mkdir -p logs
  WY17_LOG="logs/series_$(date +%Y%m%d_%H%M%S).log"
  export WY17_INNER=1 WY17_LOG
  bash "$(basename "${BASH_SOURCE[0]}")" "$@" 2>&1 | tee -a "$WY17_LOG"
  exit "${PIPESTATUS[0]}"
fi

set -u -o pipefail

NOTEBOOK="${NOTEBOOK:-Make_Sycamore_Model_WY2017sep.ipynb}"
PYTHON="${PYTHON:-python}"
SETTINGS="${SETTINGS:-Retuned Original S1 S5 S6 S7 S8 S9 S10 S11}"
CHECK_ONLY="${CHECK_ONLY:-0}"
TIMEOUT_MIN="${TIMEOUT_MIN:-240}"
HEARTBEAT_SEC="${HEARTBEAT_SEC:-300}"
MIN_FREE_GB="${MIN_FREE_GB:-5}"
RUN_HOURS="${RUN_HOURS:-6551}"
START_EXPECT="09/01/2016/00/00"
RADAR_FIRST="${RADAR_FIRST:-radar_0901201601.asc}"
RADAR_LAST="${RADAR_LAST:-radar_0531201723.asc}"
SDF="${SDF:-../init_data/met_wy2017_sep/Master_Met.sdf}"
RADAR_HOME="${RADAR_HOME:-$HOME/mrms_wk_wy2017_sep}"
LOGDIR="logs"
ROWS_PER_HOUR=16          # the outlet file gets 16 rows per model hour (same as the earlier runs)
PREFIX="Sycamore_WY2017sep_"
BASE_SET="Retuned"        # the setting the others' input files are compared with
export WANT_START="$START_EXPECT" WANT_HOURS="$RUN_HOURS" WANT_FIRST="$RADAR_FIRST" WANT_LAST="$RADAR_LAST" WANT_SDF="$SDF"

say()  { printf '\n[%s] %s\n' "$(date '+%H:%M:%S')" "$*"; }
die()  { printf '\nSTOPPED: %s\n' "$*" >&2; exit 1; }
fmt_min() { local m=$1; if [ "$m" -ge 60 ]; then printf '%d h %02d min' $((m / 60)) $((m % 60)); else printf '%d min' "$m"; fi; }

# ---------------------------------------------------------------- 1. where am I, what is there
cd "$(dirname "${BASH_SOURCE[0]}")" || die "cannot change to the script's folder."
[ "$(basename "$PWD")" = "lab" ] || die "run this from workspaces/Sycamore/lab/ (now in $PWD)."
[ -f "$NOTEBOOK" ] || die "$NOTEBOOK is not in $PWD. Run  python apply_wy2017sep_patch.py  first (it makes it from your S11 notebook)."
[ -f "$SDF" ] || die "weather station list $SDF not found. The September-2016-start weather files (met_wy2017_sep) are needed."
mkdir -p "$LOGDIR" results

"$PYTHON" -c "import nbclient, nbformat, jupyter_client, ipykernel, pytRIBS, numpy" 2>/dev/null \
  || die "this python ($PYTHON) cannot import nbclient, nbformat, jupyter_client, ipykernel, numpy and pytRIBS. Try:  pip install nbclient nbformat ipykernel   (pytRIBS should already be installed)."

TRIBS="${TRIBS_BIN:-}"
if [ -z "$TRIBS" ]; then
  if [ -x /opt/tRIBS/build/tRIBS ]; then TRIBS=/opt/tRIBS/build/tRIBS; else TRIBS="$(command -v tRIBS || true)"; fi
fi
if [ "$CHECK_ONLY" != "1" ]; then
  [ -n "$TRIBS" ] && [ -x "$TRIBS" ] || die "cannot find the tRIBS program. Give it with  TRIBS_BIN=/path/to/tRIBS ./run_wy2017sep_series.sh"
fi
if pgrep -x tRIBS >/dev/null 2>&1; then
  die "a tRIBS run is already going ($(pgrep -x tRIBS | tr '\n' ' ')). Wait for it to finish (or stop it); this script must not overlap with it."
fi

free_gb() { df -P . | awk 'NR==2 {printf "%.0f", $4/1048576}'; }
FREE_GB=$(free_gb)
[ "${FREE_GB:-0}" -ge "$MIN_FREE_GB" ] || die "only ${FREE_GB} GB of disk free here; need at least $MIN_FREE_GB."

# ---------------------------------------------------------------- 2. the radar folder
RADAR_DIR="${RADAR_DIR:-$RADAR_HOME/rain}"
if [ ! -f "$RADAR_DIR/$RADAR_LAST" ]; then
  say "$RADAR_DIR does not have $RADAR_LAST. Looking under $RADAR_HOME."
  mapfile -t HITS < <(find "$RADAR_HOME" -name "$RADAR_LAST" -type f 2>/dev/null | xargs -r -n1 dirname | sort -u)
  if [ "${#HITS[@]}" -eq 0 ]; then
    die "no folder under $RADAR_HOME holds $RADAR_LAST. The combined September-to-May radar rain files are not where I expected.
       If they are somewhere else, run:  RADAR_DIR=/folder/with/the/files ./run_wy2017sep_series.sh"
  elif [ "${#HITS[@]}" -gt 1 ]; then
    printf '  found in more than one place:\n'; printf '    %s\n' "${HITS[@]}"
    die "pick one:  RADAR_DIR=<one of the folders above> ./run_wy2017sep_series.sh"
  fi
  RADAR_DIR="${HITS[0]}"
fi
RADAR_DIR="$(cd "$RADAR_DIR" && pwd)" || die "cannot open the radar folder $RADAR_DIR."
export SYCAMORE_RADAR_DIR="$RADAR_DIR"

NRAD=$(find "$RADAR_DIR" -maxdepth 1 -name 'radar_*.asc' | wc -l)
say "Radar folder: $RADAR_DIR  ($NRAD radar files)"
"$PYTHON" - "$RADAR_DIR" <<'PY' || die "the radar folder does not hold every hourly file the run needs (see above). Nothing was run."
import os, sys
from datetime import datetime, timedelta
folder = sys.argv[1]
start = datetime.strptime(os.environ['WANT_START'], '%m/%d/%Y/%H/%M')
hours = int(os.environ['WANT_HOURS'])
names = [f"radar_{(start + timedelta(hours=h)):%m%d%Y%H}.asc" for h in range(1, hours + 1)]
have = set(os.listdir(folder))
missing = [n for n in names if n not in have]
empty = [n for n in names if n in have and os.path.getsize(os.path.join(folder, n)) == 0]
print(f'   the run needs {len(names)} radar files ({names[0]} to {names[-1]}): {len(names) - len(missing)} present, {len(missing)} missing, {len(empty)} empty')
if names[0] != os.environ['WANT_FIRST'] or names[-1] != os.environ['WANT_LAST']:
    sys.exit(f"   the start date and run length need {names[0]} to {names[-1]}, but this script expects {os.environ['WANT_FIRST']} to {os.environ['WANT_LAST']}.")
if missing:
    sys.exit(f'   MISSING, for example: {missing[:5]}')
if empty:
    sys.exit(f'   EMPTY (0 bytes), for example: {empty[:5]}')
PY

# ---------------------------------------------------------------- 3. the settings
declare -A FULL TAGS PAIRS GRIDS BEDADD GWFRAC ASV AUV
define() {   # short-name full-name tag Ks f PsiB As Au bedrock-add-m GW-fraction
  local s="$1"
  FULL[$s]="$2"; TAGS[$s]="$3"; ASV[$s]="$7"; AUV[$s]="$8"; BEDADD[$s]="$9"; GWFRAC[$s]="${10}"
  PAIRS[$s]="KS_MULT=$4 F_MULT=$5 PSIB_MULT=$6 ANISOTROPY=$7 ANISOTROPY_U=$8 BEDROCK_ADD_M=$9 GW_FRACTION=${10} OPTPERCOLATION=0 OPTSNOW=1 M_MULT=1 THETAS_MULT=1"
  GRIDS[$s]="KS=$4 PB=$6 FD=$5 PI=1 TS=1 PO=1 TR=1"        # variable=multiplier; 1 means a plain copy of the SOLUS grid
}
#       short     full name                  tag                      Ks      f     PsiB    As   Au  bedrock GW
define  Retuned   Re-tuned                   Retuned                  0.19261 0.15  0.2939  300  1   1       0.95
define  Original  Original                   Original                 0.187   0.5874 0.2939 1    1   1       0.95
define  S1        S1-Ks1.0                   S1_Ks1p0                 1       0.15  0.2939  300  1   1       0.95
define  S5        S5-Ks1.0-Bedrock2.0        S5_Ks1p0_Bed2p0          1       0.15  0.2939  300  1   2       0.95
define  S6        S6-Ks0.5-Bedrock2.0        S6_Ks0p5_Bed2p0          0.5     0.15  0.2939  300  1   2       0.95
define  S7        S7-Ks1.0-Bedrock2.0-As500  S7_Ks1p0_Bed2p0_As500    1       0.15  0.2939  500  1   2       0.95
define  S8        S8-Ks1.0-Bedrock2.0-f0.02  S8_Ks1p0_Bed2p0_f0p02    1       0.02  0.2939  300  1   2       0.95
define  S9        S9-Ks1.0-Bedrock0.5        S9_Ks1p0_Bed0p5          1       0.15  0.2939  300  1   0.5     0.95
define  S10       S10-Ks1.0-GW0.5            S10_Ks1p0_GW0p5          1       0.15  0.2939  300  1   1       0.5
define  S11       S11-Ks1.0-PsiB1.0          S11_Ks1p0_PsiB1p0        1       0.15  1       300  1   1       0.95

norm_name() {   # Re-tuned / retuned -> Retuned, s5 -> S5, and so on; prints nothing for a name it does not know
  local s="${1,,}"; s="${s//-/}"
  case "$s" in
    retuned) echo Retuned ;;
    original) echo Original ;;
    s1|s5|s6|s7|s8|s9|s10|s11) echo "${s^^}" ;;
  esac
}
PLAN=()
for w in $SETTINGS; do
  n="$(norm_name "$w")"
  [ -n "$n" ] || die "unknown setting '$w'. Use any of: Retuned Original S1 S5 S6 S7 S8 S9 S10 S11."
  for p in ${PLAN[@]+"${PLAN[@]}"}; do [ "$p" != "$n" ] || die "setting $n is listed twice."; done
  PLAN+=("$n")
done
[ "${#PLAN[@]}" -gt 0 ] || die "no settings to run."

# ---------------------------------------------------------------- helpers
check_settings_cell() {   # $1 = full setting name, $2 = tag, $3... = KEY=value pairs the setting must give
  # Runs ONLY the Settings cell of the notebook (nothing else, no files) and compares what it gives with what this setting is supposed to be.
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
    sys.exit(f'the notebook does not have this setting ({e}). Did apply_wy2017sep_patch.py run, and is NOTEBOOK the file it wrote?')
bad = []
for p in pairs:
    k, v = p.split('=')
    try:
        if abs(float(ns[k]) - float(v)) > 1e-12:
            bad.append(f'{k} is {ns[k]} (this setting needs {v})')
    except (TypeError, ValueError, KeyError):
        bad.append(f'{k} is {ns.get(k)!r} (this setting needs {v})')
for key, want in [('START_DATE', os.environ['WANT_START']), ('RUN_HOURS', int(os.environ['WANT_HOURS'])), ('RUN_CHOICE', 'wy2017sep'),
                  ('RADAR_FIRST_FILE', os.environ['WANT_FIRST']), ('RADAR_LAST_FILE', os.environ['WANT_LAST']), ('MET_STATIONS', os.environ['WANT_SDF'])]:
    if ns.get(key) != want:
        bad.append(f'{key} is {ns.get(key)!r} (needed {want!r})')
want_in, want_res = f'Sycamore_WY2017sep_{tag}.in', f'Sycamore_WY2017sep_{tag}'
if ns['INPUT_FILE'] != want_in or ns['RESULTS_NAME'] != want_res:
    bad.append(f"names are {ns['INPUT_FILE']} / {ns['RESULTS_NAME']} (needed {want_in} / {want_res})")
if bad:
    sys.exit('the Settings cell does not give the values this setting needs:\n   - ' + '\n   - '.join(bad))
print('   ' + ', '.join(f'{k}={ns[k]}' for k in ('KS_MULT', 'F_MULT', 'PSIB_MULT', 'ANISOTROPY', 'ANISOTROPY_U', 'BEDROCK_ADD_M', 'GW_FRACTION', 'OPTPERCOLATION')))
print(f"   run window: {ns['START_DATE']}, {ns['RUN_HOURS']} hours; radar {ns['RADAR_FIRST_FILE']} to {ns['RADAR_LAST_FILE']}; weather {ns['MET_STATIONS']}")
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
                             display_name="wy2017sep", language="python")
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

scan_executed() {   # $1 = executed notebook, $2 = As, $3 = Au. Reads the notebook's own printout; prints no flow values.
  "$PYTHON" - "$1" "$2" "$3" <<'PY'
import json, re, sys
path, want_as, want_au = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
nb = json.load(open(path, encoding='utf-8'))
lines = []
for c in nb.get('cells', []):
    for o in c.get('outputs', []):
        if o.get('output_type') == 'stream':
            t = o.get('text', '')
            lines += (''.join(t) if isinstance(t, list) else t).splitlines()
        elif o.get('output_type') == 'error':
            lines.append(f"CELL ERROR: {o.get('ename')}: {o.get('evalue')}")
bad = [l for l in lines if l.startswith(('CHECK', 'MISSING', 'PROBLEM', 'CELL ERROR')) or l.startswith('WARNING: the working folder')]
if bad:
    print('   the executed notebook printed lines that need a look:')
    for l in bad[:15]:
        print('     ' + l[:200])
    sys.exit(1)
info = re.compile(r'^(OK\s+(Weather|station|Radar|Radar first|Radar last)|OK: every hourly|Run window:|\d+ radar files|This run needs|SETTING:|Run choice:|Ks 0-5|PsiB 0-5|f inside|m \(pore|thetaS 0-5|BEDROCKFILE:|GWATERFILE:|  depth to bedrock, m|  starting water table|Every expected setting)')
shown = [l for l in lines if info.match(l)]
print('   from the notebook\'s own printout:')
for l in shown[:30]:
    print('     ' + l.strip()[:190])
cand = [l for l in lines if 'As and Au in the soil table in use:' in l]
if not cand:
    print('   the executed notebook did not print the line "As and Au in the soil table in use", so As and Au could not be read back.')
    sys.exit(1)
m = re.search(r'in use:\s*\[(.*)\]\s*\(ANISOTROPY\s*=\s*([-+0-9.eE]+)\)', cand[0])
if not m:
    print('   could not read this line: ' + cand[0][:200])
    sys.exit(1)
table_txt = re.sub(r'(?:np|numpy)\.\w+\(', '(', m.group(1))      # numpy 2 prints np.float64(300.0): keep just the number
nums = [float(x) for x in re.findall(r'[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?', table_txt)]
if len(nums) != 2:
    print(f'   the soil table in use has more than one As/Au pair or none: {m.group(1)}')
    sys.exit(1)
if abs(nums[0] - want_as) > 1e-9 or abs(nums[1] - want_au) > 1e-9 or abs(float(m.group(2)) - want_as) > 1e-9:
    print(f'   the soil table in use has As {nums[0]:g}, Au {nums[1]:g} (ANISOTROPY {m.group(2)}); this setting needs As {want_as:g}, Au {want_au:g}.')
    sys.exit(1)
print(f'   As and Au read back from the soil table in the model folder: As {nums[0]:g}, Au {nums[1]:g}. OK')
PY
}

check_input_keywords() {   # $1 = input file, $2 = tag. Prints the lines it checks.
  "$PYTHON" - "$1" "$2" "$RADAR_DIR" <<'PY'
import os, re, sys
infile, tag, radar_dir = sys.argv[1:4]
want_hours, want_start = int(os.environ['WANT_HOURS']), os.environ['WANT_START']


def read_keywords(path):
    with open(path) as fh:
        lines = fh.read().splitlines()
    found = {}
    for i, line in enumerate(lines):
        m = re.match(r'^([A-Z0-9]+):\s', line + ' ')
        if m and i + 1 < len(lines):
            found[m.group(1)] = lines[i + 1].strip()
    return found


def same(a, b):
    try:
        return float(a) == float(b)
    except (TypeError, ValueError):
        return str(a) == str(b)


kw = read_keywords(infile)
bad = []
for key, want in [('RAINSOURCE', 1), ('RAININTRVL', 1), ('METDATAOPTION', 1), ('OPTSOILTYPE', 1), ('OPTBEDROCK', 1), ('OPTGWFILE', 0),
                  ('OPTLANDUSE', 0), ('OPTMESHINPUT', 1), ('OPTSNOW', 1), ('OPTPERCOLATION', 0), ('STARTDATE', want_start), ('RUNTIME', want_hours)]:
    got = kw.get(key)
    if got is None or not same(got, want):
        bad.append(f'{key} is {got!r} (needed {want!r})')
rf = kw.get('RAINFILE')
if rf is None or os.path.abspath(rf) != os.path.abspath(radar_dir + '/radar_'):
    bad.append(f'RAINFILE is {rf!r} (needed {radar_dir}/radar_)')
if not kw.get('HYDROMETSTATIONS', '').endswith('met_wy2017_sep/Master_Met.sdf'):
    bad.append(f"HYDROMETSTATIONS is {kw.get('HYDROMETSTATIONS')!r} (needed a path ending met_wy2017_sep/Master_Met.sdf)")
if not kw.get('OUTFILENAME', '').endswith(f'Sycamore_WY2017sep_{tag}'):
    bad.append(f"OUTFILENAME is {kw.get('OUTFILENAME')!r} (needed a name ending Sycamore_WY2017sep_{tag})")
print(f"   STARTDATE {kw.get('STARTDATE')}, RUNTIME {kw.get('RUNTIME')}, OPTPERCOLATION {kw.get('OPTPERCOLATION')} (channel loss "
      f"{'off' if str(kw.get('OPTPERCOLATION')) == '0' else 'ON'}), OPTSNOW {kw.get('OPTSNOW')}")
print(f"   OUTFILENAME {kw.get('OUTFILENAME')}")
print(f"   RAINFILE {kw.get('RAINFILE')} | HYDROMETSTATIONS {kw.get('HYDROMETSTATIONS')}")
if bad:
    print('   the input file is not what this run needs:')
    for b in bad:
        print('     - ' + b)
    sys.exit(1)
PY
}

keyword_value() {  # $1 = input file, $2 = KEYWORD: the line after the "KEYWORD:" line
  awk -v k="$2:" '$1==k {getline v; print v; exit}' "$1"
}

compare_inputs() {   # $1 new input file, $2 base input file, $3 new tag, $4 base tag, $5 new bedrock m, $6 new GW fraction, $7 base bedrock m, $8 base GW fraction
  "$PYTHON" - "$@" <<'PY'
import re, sys
new, base, tag_new, tag_base, bed_new, gw_new, bed_base, gw_base = sys.argv[1:9]


def read_keywords(path):
    with open(path) as fh:
        lines = fh.read().splitlines()
    found = {}
    for i, line in enumerate(lines):
        m = re.match(r'^([A-Z0-9]+):\s', line + ' ')
        if m and i + 1 < len(lines):
            found[m.group(1)] = lines[i + 1].strip()
    return found


def norm(d, t):
    return {k: v.replace(f'Sycamore_WY2017sep_{t}', 'Sycamore_WY2017sep_<RUN>') for k, v in d.items()}


a, b = norm(read_keywords(base), tag_base), norm(read_keywords(new), tag_new)
bed_moved = float(bed_new) != float(bed_base)
gw_moved = bed_moved or float(gw_new) != float(gw_base)
allowed = set()
if bed_moved:
    allowed.add('BEDROCKFILE')
if gw_moved:
    allowed.add('GWATERFILE')
problems, changed = [], []
if set(a) != set(b):
    problems.append(f'the keyword lists differ: only in the base {sorted(set(a) - set(b))}, only in the new file {sorted(set(b) - set(a))}')
for k in sorted(set(a) & set(b)):
    if a[k] != b[k]:
        if k in allowed:
            changed.append(k)
        else:
            problems.append(f'{k}: base "{a[k]}", new "{b[k]}"')
for k in sorted(allowed):
    if k in a and k in b and a[k] == b[k]:
        problems.append(f'{k} is the same as in the base input file, but this setting moves it')
if problems:
    print('   the input file differs from the base (Re-tuned) input file in more than it should:')
    for p in problems[:12]:
        print('     - ' + p)
    sys.exit(1)
if changed:
    print(f"   differs from the base input file only in the results name and {' and '.join(changed)} (as expected for this setting): "
          + '; '.join(f'{k} {b[k].split("/")[-1]}' for k in changed))
else:
    print('   identical to the base input file apart from the results name. Good: nothing else in the input file changed.')
PY
}

check_soil_grids() {   # $1 = when the build started (seconds), $2... = VAR=multiplier (1 = plain copy of the SOLUS grid)
  # Looks at the soil grids the notebook just wrote in data/model/soil. Reads files only; prints no flow values.
  local t0="$1"; shift
  "$PYTHON" - "$t0" "$@" <<'PY'
import csv, filecmp, os, sys
from pathlib import Path
import numpy as np

lab = Path.cwd()
MODEL_SOIL, INIT_DATA, F_FLOOR = 'data/model/soil', '../init_data', 1e-4
t0 = float(sys.argv[1])
want = {}
for p in sys.argv[2:]:
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
            fail(f"scgrid.gdf row {var} is '{mb}', which has a multiplier in its name, but this setting needs x1.")
        if not filecmp.cmp(mf, sf, shallow=False):
            fail(f'soil grid {var}: {mf.name} differs from the SOLUS grid {sf.name}, but this setting needs x1 (an unchanged copy).')
        print(f'   {var} x1: {mf.name} in data/model/soil is a byte-for-byte copy of the SOLUS grid {sf.name}')
        continue
    if not name.endswith(f'_x{mult:g}'):
        fail(f"scgrid.gdf row {var} is '{mb}', which does not end in _x{mult:g}, the multiplier this setting needs.")
    if mf.stat().st_mtime < t0 - 2:
        fail(f'soil grid {var}: {mf.name} is older than this build (a stale file).')
    hm, gm = read_asc(mf)
    hs, gs = read_asc(sf)
    if gm.shape != gs.shape:
        fail(f'soil grid {var}: the scaled grid and the SOLUS grid have different shapes.')
    nd_s, nd_m = hs.get('nodata_value', -9999), hm.get('nodata_value', -9999)
    valid = gs != nd_s
    if int(valid.sum()) < 10:
        fail(f'soil grid {var}: fewer than 10 cells with data to check.')
    if np.any(gm[valid] == nd_m):
        fail(f'soil grid {var}: some cells with data in the SOLUS grid are no-data in the scaled grid.')
    exp = gs * mult
    if var == 'FD':
        exp = np.maximum(exp, F_FLOOR)                     # f has a floor of 0.0001: cells below it are raised to it
    tol = np.maximum(0.02 * np.abs(exp), 1e-7 if var == 'FD' else 2e-4)      # f is written with 8 decimals, the others with 4
    dev = np.abs(gm - exp)
    nbad = int((valid & (dev > tol)).sum())
    if nbad:
        worst = np.where(valid, dev - tol, -np.inf).argmax()
        fail(f'soil grid {var}: {nbad} of {int(valid.sum())} cells are not the SOLUS grid x{mult:g} (worst: SOLUS {gs.flat[worst]:.6g}, '
             f'scaled {gm.flat[worst]:.6g}, expected {exp.flat[worst]:.6g}).')
    note = ''
    if var == 'FD':
        note = f'; {int((valid & (exp <= F_FLOOR * 1.0000001)).sum())} of {int(valid.sum())} cells sit at the 0.0001 floor'
    print(f'   {var} x{mult:g}: {mf.name} is the SOLUS grid x{mult:g} in all {int(valid.sum())} cells (largest difference {float(dev[valid].max()):.3g}){note}')
PY
}

check_depth_grids() {   # $1 = input file, $2 = bedrock add (m), $3 = GW fraction, $4 = when the build started (seconds)
  "$PYTHON" - "$@" <<'PY'
import os, re, sys
from pathlib import Path
import numpy as np

infile, add, frac, t0 = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
lab = Path.cwd()


def fail(msg):
    sys.exit('   ' + msg)


def read_keywords(path):
    with open(path) as fh:
        lines = fh.read().splitlines()
    found = {}
    for i, line in enumerate(lines):
        m = re.match(r'^([A-Z0-9]+):\s', line + ' ')
        if m and i + 1 < len(lines):
            found[m.group(1)] = lines[i + 1].strip()
    return found


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


kw = read_keywords(infile)
for key in ('BEDROCKFILE', 'GWATERFILE'):
    if not kw.get(key):
        fail(f'{infile} has no {key} line, so the bedrock and starting water table could not be checked.')
bed_path, gw_path = lab / kw['BEDROCKFILE'], lab / kw['GWATERFILE']
orig_path = lab / '../init_data/bedrock_depth_m.asc'
for p in (bed_path, gw_path, orig_path):
    if not p.is_file():
        fail(f'{p} not found.')
for p in (bed_path, gw_path):
    if p.stat().st_mtime < t0 - 2:
        fail(f'{p.name} in the model folder is older than this build (a stale file).')

fmt = f'{add:g}'.replace('.', 'p')
want_bed, want_gw = f'bedrock_depth_plus{fmt}m.asc', f'InitGW_{int(round(frac * 100))}pct_plus{fmt}m_mm.asc'
for got, want in ((bed_path.name, want_bed), (gw_path.name, want_gw)):
    if got != want:
        print(f'   NOTE: the file is called {got}, not {want} (the values are checked below either way).')

hb, bed = read_asc(bed_path)
ho, orig = read_asc(orig_path)
hg, gw = read_asc(gw_path)
nd_b, nd_o, nd_g = hb.get('nodata_value', -9999), ho.get('nodata_value', -9999), hg.get('nodata_value', -9999)
if not (bed.shape == orig.shape == gw.shape):
    fail(f'the bedrock, SOLUS bedrock and water table grids have different shapes: {bed.shape}, {orig.shape}, {gw.shape}.')
valid = orig != nd_o
if np.any(bed[valid] == nd_b):
    fail('some cells with SOLUS bedrock data are no-data in the model bedrock grid.')
diff = bed[valid] - orig[valid]
off = np.abs(diff - add)
if float(off.max()) > 0.011:
    fail(f'the bedrock grid is SOLUS + {float(np.median(diff)):.3f} m (range {diff.min():.3f} to {diff.max():.3f} m), but this setting needs + {add:g} m.')
use = valid & (bed > 0)
if int(use.sum()) < 10:
    fail('fewer than 10 cells with bedrock depth above 0 to check the water table on.')
if np.any(gw[use] == nd_g):
    fail('some cells with bedrock data are no-data in the starting water table grid.')
exp = frac * 1000.0 * bed[use]
dev = np.abs(gw[use] - exp)
tol = 0.005 * np.abs(exp) + 1.0
if np.any(dev > tol):
    i = (dev - tol).argmax()
    fail(f'the starting water table is not {frac:g} x the bedrock depth (in mm) in {int((dev > tol).sum())} cells (worst: bedrock {bed[use][i]:.2f} m, '
         f'water table {gw[use][i]:.1f} mm, expected {exp[i]:.1f} mm).')
print(f'   bedrock {bed_path.name}: SOLUS + {add:g} m in all {int(valid.sum())} cells with data (median depth {np.median(bed[valid]):.2f} m)')
print(f'   starting water table {gw_path.name}: {frac:g} x bedrock depth (mm) in all {int(use.sum())} cells above 0 m (median {np.median(gw[use]):.0f} mm); '
      'both files were written into the model folder in this build')
PY
}

run_files() {   # the result files of one run (not the Original-Kchan runs, whose names start NAME_kchan)
  local name="$1" f
  for f in results/"${name}"_* results/"${name}".*; do
    [ -e "$f" ] || continue
    case "$f" in results/"${name}"_kchan*) continue ;; esac
    printf '%s\n' "$f"
  done
}

OUTLET_MSG=""
outlet_ok() {   # $1 = outlet file. Counts only: rows, last hour, non-numbers. Sets OUTLET_MSG. Returns 0 if the run reached its end.
  local outlet="$1" rows last bad expected
  [ -f "$outlet" ] || { OUTLET_MSG="$outlet does not exist"; return 1; }
  expected=$(( RUN_HOURS * ROWS_PER_HOUR ))
  rows=$(awk -F, 'NR>1 && NF>=2 {n++} END {print n+0}' "$outlet")
  last=$(awk -F, 'NR>1 && NF>=2 {t=$1} END {print t+0}' "$outlet")
  bad=$(awk -F, 'NR>1 && NF>=2 && ($2 ~ /[nN][aA][nN]|[iI][nN][fF]/) {n++} END {print n+0}' "$outlet")
  OUTLET_MSG="$rows rows (about $expected expected), last model hour $last (run is $RUN_HOURS), rows with NaN/Inf: $bad"
  [ "$rows" -ge $(( expected * 99 / 100 )) ] || { OUTLET_MSG="$OUTLET_MSG; too few rows"; return 1; }
  [ "$bad" -eq 0 ] || { OUTLET_MSG="$OUTLET_MSG; NaN or Inf in the outlet file"; return 1; }
  awk -v l="$last" -v r="$RUN_HOURS" 'BEGIN {exit !(l >= r - 2)}' || { OUTLET_MSG="$OUTLET_MSG; stops before the end of the run"; return 1; }
  return 0
}

TRIBS_PID=""
cleanup() {
  if [ -n "$TRIBS_PID" ] && kill -0 "$TRIBS_PID" 2>/dev/null; then
    kill -TERM -- "-$TRIBS_PID" 2>/dev/null || kill -TERM "$TRIBS_PID" 2>/dev/null
    printf '\n[%s] Stopped tRIBS (pid %s).\n' "$(date '+%H:%M:%S')" "$TRIBS_PID"
  fi
}
trap cleanup EXIT
trap 'printf "\nSTOPPED: interrupted.\n" >&2; exit 130' INT TERM

CUR_LABEL=""
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
    say "$CUR_LABEL: tRIBS still running: ${elapsed} min so far; outlet file has ${rows} rows of about ${expected} ($(( 100 * rows / expected ))%)"
    if [ "$elapsed" -ge "$TIMEOUT_MIN" ]; then
      cleanup; TRIBS_PID=""
      die "tRIBS did not finish within ${TIMEOUT_MIN} minutes and was stopped."
    fi
  done
  wait "$TRIBS_PID"; local rc=$?
  TRIBS_PID=""
  return $rc
}

# ---------------------------------------------------------------- 4. the plan: what is finished, what is left
declare -A STATUS
N_TODO=0
say "Plan: ${#PLAN[@]} setting(s), in this order. Run window $START_EXPECT, $RUN_HOURS hours (September 2016 is spin-up, not scored)."
printf '  %-3s %-9s %-26s %-44s %s\n' "#" "short" "setting" "results name" "now"
i=0
for S in "${PLAN[@]}"; do
  i=$((i + 1))
  NAME="${PREFIX}${TAGS[$S]}"
  OUTLET="results/${NAME}_Outlet.qout"
  if [ -f "$OUTLET" ] && outlet_ok "$OUTLET"; then STATUS[$S]="COMPLETE"
  elif [ -n "$(run_files "$NAME")" ]; then STATUS[$S]="PARTIAL"
  else STATUS[$S]="TODO"; fi
  [ "${STATUS[$S]}" = "COMPLETE" ] || N_TODO=$((N_TODO + 1))
  printf '  %-3s %-9s %-26s %-44s %s\n' "$i" "$S" "${FULL[$S]}" "$NAME" "${STATUS[$S]}"
done
printf '\n  COMPLETE = its outlet file reached the end of the run, so it is skipped. PARTIAL = leftover files from a run that did not finish: they are moved to\n'
printf '  results/_partial_<date>/ and the run starts again. TODO = nothing there yet.\n'
if [ "$CHECK_ONLY" = "1" ]; then
  printf '  CHECK_ONLY=1: inputs are built and checked, tRIBS is not run, and nothing in results/ is moved.\n'
else
  say "About 50 to 60 minutes per run plus 1 to 2 minutes to build the inputs (an estimate from the WY2020 timing, not a measurement), so $N_TODO run(s) is very roughly $(( N_TODO * 55 / 60 )) to $(( N_TODO * 62 / 60 + 1 )) hours."
fi

# ---------------------------------------------------------------- 5. the settings, one after the other
START_ALL=$(date +%s)
SUMMARY=()
N_DONE=0
DONE_SECS=0
IDX=0
for S in "${PLAN[@]}"; do
  IDX=$((IDX + 1))
  SETTING="${FULL[$S]}"; TAG="${TAGS[$S]}"
  NAME="${PREFIX}${TAG}"; INFILE="${NAME}.in"; OUTLET="results/${NAME}_Outlet.qout"
  CUR_LABEL="[$IDX/${#PLAN[@]}] $S"
  if [ "${STATUS[$S]}" = "COMPLETE" ]; then
    say "=== $CUR_LABEL ($SETTING): already finished ($OUTLET); skipping ==="
    SUMMARY+=("$S: already finished before this run (skipped)")
    continue
  fi
  T_SETTING=$(date +%s)
  [ "$(free_gb)" -ge "$MIN_FREE_GB" ] || die "only $(free_gb) GB of disk free; need at least $MIN_FREE_GB. Finished runs are kept: free some space and run the same command again."

  say "=== $CUR_LABEL ($SETTING): checking the Settings cell ==="
  # shellcheck disable=SC2086
  check_settings_cell "$SETTING" "$TAG" ${PAIRS[$S]} || die "the notebook $NOTEBOOK is not set up for $S (see the lines above). Nothing was run for it."

  say "=== $CUR_LABEL: building the input file ($NOTEBOOK, headless; a minute or two) ==="
  rm -f "$INFILE"
  T_BUILD=$(date +%s)
  EXECUTED="$LOGDIR/executed_Make_Sycamore_Model_WY2017sep_${TAG}.ipynb"
  build_inputs "$SETTING" "$EXECUTED" > "$LOGDIR/notebook_${TAG}.log" 2>&1 \
    || die "the notebook failed for $S. Open $EXECUTED to see the cell that stopped (details in $LOGDIR/notebook_${TAG}.log)."
  [ -f "$INFILE" ] || die "the notebook finished but $INFILE was not written."

  say "=== $CUR_LABEL: reading the executed notebook's printout ==="
  scan_executed "$EXECUTED" "${ASV[$S]}" "${AUV[$S]}" || die "the executed notebook ($EXECUTED) has a problem for $S (see the lines above). tRIBS was not run."

  say "=== $CUR_LABEL: checking the keywords in $INFILE ==="
  check_input_keywords "$INFILE" "$TAG" || die "$INFILE is not what $S needs (see the lines above). tRIBS was not run."
  OUTP=$(keyword_value "$INFILE" OUTFILENAME)
  OUTLET="${OUTP}_Outlet.qout"

  say "=== $CUR_LABEL: comparing $INFILE with the base ($BASE_SET) input file ==="
  BASE_IN="${PREFIX}${TAGS[$BASE_SET]}.in"
  if [ "$S" = "$BASE_SET" ]; then
    echo "   this is the base input file; nothing to compare it with. The other settings are compared with it."
  elif [ ! -f "$BASE_IN" ]; then
    echo "   NOT CHECKED: the base input file $BASE_IN is not in $PWD (it is built when Retuned is in the list), so $INFILE was not compared with it."
    SUMMARY+=("$S: NOT compared with the Re-tuned input file (it was not there)")
  else
    compare_inputs "$INFILE" "$BASE_IN" "$TAG" "${TAGS[$BASE_SET]}" "${BEDADD[$S]}" "${GWFRAC[$S]}" "${BEDADD[$BASE_SET]}" "${GWFRAC[$BASE_SET]}" \
      || die "$INFILE is not the Re-tuned input file plus what $S changes. tRIBS was not run."
  fi

  say "=== $CUR_LABEL: checking the soil grids the notebook wrote in data/model/soil ==="
  # shellcheck disable=SC2086
  check_soil_grids "$T_BUILD" ${GRIDS[$S]} || die "the soil grids in data/model/soil are not what $S needs (see the lines above). tRIBS was not run."

  say "=== $CUR_LABEL: checking the bedrock and starting water table grids ==="
  check_depth_grids "$INFILE" "${BEDADD[$S]}" "${GWFRAC[$S]}" "$T_BUILD" || die "the bedrock or water table grid is not what $S needs (see the lines above). tRIBS was not run."

  if [ "$CHECK_ONLY" = "1" ]; then
    say "CHECK_ONLY: $S input file is built and checked ($(( $(date +%s) - T_SETTING )) s). Not running tRIBS."
    SUMMARY+=("$S: input file built and checked ($INFILE), tRIBS not run")
    continue
  fi

  if [ "${STATUS[$S]}" = "PARTIAL" ] || [ -n "$(run_files "$NAME")" ]; then
    MOVED="results/_partial_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$MOVED" || die "cannot make $MOVED."
    while IFS= read -r f; do mv -- "$f" "$MOVED/" || die "cannot move $f to $MOVED."; done < <(run_files "$NAME")
    say "Leftover files of an unfinished $S run were moved to $MOVED/ (not deleted)."
  fi
  [ ! -e "$OUTLET" ] || die "$OUTLET already exists. Move or delete the old $S results first (this script will not overwrite them)."

  say "=== $CUR_LABEL: running tRIBS (log: $LOGDIR/simulation_WY2017sep_${TAG}.log) ==="
  T1=$(date +%s)
  run_tribs "$INFILE" "$LOGDIR/simulation_WY2017sep_${TAG}.log" "$OUTLET" \
    || die "tRIBS exited with an error for $S; see $LOGDIR/simulation_WY2017sep_${TAG}.log (last lines: $(tail -n 3 "$LOGDIR/simulation_WY2017sep_${TAG}.log" | tr '\n' ' '))."
  outlet_ok "$OUTLET" || die "outlet file check for $S failed: $OUTLET_MSG"
  say "Outlet file check: $OUTLET_MSG"
  RUN_MIN=$(( ($(date +%s) - T1) / 60 ))
  SIZE=$(run_files "$NAME" | xargs -r du -ch 2>/dev/null | tail -n 1 | awk '{print $1}')
  SUMMARY+=("$S: finished, tRIBS ${RUN_MIN} min; outlet file $OUTLET; results about ${SIZE:-?}")
  N_DONE=$((N_DONE + 1)); DONE_SECS=$(( DONE_SECS + $(date +%s) - T_SETTING ))
  LEFT=0
  for S2 in "${PLAN[@]:$IDX}"; do [ "${STATUS[$S2]}" = "COMPLETE" ] || LEFT=$((LEFT + 1)); done
  MSG="$CUR_LABEL finished in $(fmt_min $(( ($(date +%s) - T_SETTING) / 60 ))). Elapsed $(fmt_min $(( ($(date +%s) - START_ALL) / 60 ))); $LEFT run(s) left"
  [ "$LEFT" -gt 0 ] && MSG="$MSG, about $(fmt_min $(( DONE_SECS / N_DONE * LEFT / 60 ))) more (an estimate from the runs so far)"
  say "$MSG. Free disk: $(free_gb) GB."
done

if [ "$CHECK_ONLY" = "1" ]; then
  say "DRY RUN DONE in $(fmt_min $(( ($(date +%s) - START_ALL) / 60 ))): every input file was built and checked, and tRIBS was not run."
else
  say "ALL DONE in $(fmt_min $(( ($(date +%s) - START_ALL) / 60 )))."
fi
printf '  %s\n' ${SUMMARY[@]+"${SUMMARY[@]}"}
cat <<'MSG'

These runs have not been opened, scored or plotted by this script.
Time_hr 0 in these results is Sep 1 2016 00:00 (September is spin-up). Oct 1 2016 00:00 is hour 720; the scored window is hours 720 to 6551
(5,832 hours). Scoring scripts that assume an Oct 1 start need that 720-hour offset.
MSG
