"""
apply_wy2020_screens_patch.py  --  one-time patch that adds the four one-change screens (S1 to S4) to
Make_Sycamore_Model_WY2020.ipynb

WHAT IT DOES
  Reads your notebook and writes a NEW file, Make_Sycamore_Model_WY2020_screens.ipynb. Your original notebook is never
  changed, and the script will not overwrite an existing output file unless you add --force.

  The new notebook has four more named settings. Each is the Re-tuned setting with ONE value changed:
      S1-Ks1.0        Ks multiplier 1.0 (SOLUS as is)                   Re-tuned: 0.19261
      S2-Ks0.5        Ks multiplier 0.5                                 Re-tuned: 0.19261
      S3-Bedrock0.5   bedrock depth = SOLUS + 0.5 m                     Re-tuned: SOLUS + 1 m
      S4-GW0.5        starting water table = 0.5 x bedrock depth        Re-tuned: 0.95 x bedrock depth
  The four existing settings (Re-tuned, Original, Original-Kchan0.5, Original-Kchan0.2) keep exactly the values they have now.

  To make S3 and S4 possible, two numbers that used to be fixed in the Settings cell (BEDROCK_ADD_M = 1.0 and
  GW_FRACTION = 0.95) now come from the chosen setting, with those same values as the default. Two file names also change
  for non-whole numbers only (0.5 becomes 0p5, so a file name has no extra dot). Whole numbers such as 1 give the same
  names as before.

WHY A PATCH AND NOT A WHOLE NEW NOTEBOOK
  I could read the text of your notebook but did not have the file itself, and retyping a notebook this size by hand risks
  a typo slipping in. This script edits your own copy instead and checks the result.

HOW IT PROTECTS YOU
  - Every edit looks for one exact piece of text in one named cell and must find it exactly once. If your notebook differs
    from the copy I read (for example you edited that cell since Nov 4), the script says which edit failed and writes nothing.
  - After editing, it checks that each changed code cell still parses, then runs ONLY the Settings cell (no files, no tRIBS) once
    per setting, and checks that:
        * the four existing settings give the same values as before the patch, and
        * each screen differs from Re-tuned in exactly one value.
    It prints the table so you can read it. If any check fails, nothing is written.

RUN IT from workspaces/Sycamore/lab/ (nothing else is touched):
    python apply_wy2020_screens_patch.py
    python apply_wy2020_screens_patch.py --in Make_Sycamore_Model_WY2020.ipynb --out Make_Sycamore_Model_WY2020_screens.ipynb
    python apply_wy2020_screens_patch.py --force        # allow overwriting the output file
"""
import argparse
import ast
import contextlib
import io
import json
import os
import sys

# The cells that get edited, by their notebook cell ids (as in the copy read on 2026-11-07).
CELL_INTRO = '69ebac8e'      # markdown: how to run it
CELL_SETTINGS_MD = '2fb78b48'  # markdown: notes on the Settings cell
CELL_SETTINGS = 'b5dfbf73'   # code: the Settings cell
CELL_BEDROCK = '7e90bc2a'    # code: deepen the bedrock grid
CELL_GW = '2ee2f5dd'         # code: initial groundwater

SCREEN_BLOCK = """# One-change screens from Re-tuned (added 2026-11-07; the plan is in Soil_Storage_Running_Log.md). Each one is built from the
# Re-tuned entry with ONE value changed, so nothing else can drift. They are probes, not candidate settings: nothing is re-fit on WY2020.
_retuned = SETTINGS['Re-tuned']
SETTINGS['S1-Ks1.0'] = dict(_retuned, KS_MULT=1.0, tag='S1_Ks1p0')                  # Ks = SOLUS as is (Re-tuned: 0.19261)
SETTINGS['S2-Ks0.5'] = dict(_retuned, KS_MULT=0.5, tag='S2_Ks0p5')                  # Ks = half of SOLUS
SETTINGS['S3-Bedrock0.5'] = dict(_retuned, BEDROCK_ADD_M=0.5, tag='S3_Bed0p5')      # bedrock = SOLUS + 0.5 m (Re-tuned: + 1 m)
SETTINGS['S4-GW0.5'] = dict(_retuned, GW_FRACTION=0.5, tag='S4_GW0p5')              # starting water table = 0.5 x bedrock depth (Re-tuned: 0.95)
"""

# (cell id, text to find, text to put in its place, required, what it does)
EDITS = [
    (CELL_SETTINGS, "if SETTING not in SETTINGS:",
     SCREEN_BLOCK + "if SETTING not in SETTINGS:", True,
     "add the four screens to the SETTINGS table"),
    (CELL_SETTINGS, "GW_FRACTION = 0.95",
     "GW_FRACTION = SETTINGS[SETTING].get('GW_FRACTION', 0.95)", True,
     "starting water table fraction now comes from the setting (default 0.95)"),
    (CELL_SETTINGS, "BEDROCK_ADD_M = 1.0",
     "BEDROCK_ADD_M = SETTINGS[SETTING].get('BEDROCK_ADD_M', 1.0)", True,
     "meters added to SOLUS bedrock depth now come from the setting (default 1.0)"),
    (CELL_SETTINGS, "As {ANISOTROPY}, Au {ANISOTROPY_U}')",
     "As {ANISOTROPY}, Au {ANISOTROPY_U}')\n"
     "print(f'  bedrock SOLUS + {BEDROCK_ADD_M:g} m | starting water table {GW_FRACTION:g} x bedrock depth | channel loss OPTPERCOLATION = {OPTPERCOLATION}')", True,
     "print the bedrock and water table values with the setting, so you can see which setting is in use"),
    (CELL_SETTINGS, "# Two more settings, Original-Kchan0.5 and Original-Kchan0.2, are the Original with channel loss ON (added 2026-11-04).",
     "# Two more settings, Original-Kchan0.5 and Original-Kchan0.2, are the Original with channel loss ON (added 2026-11-04).\n"
     "# Four more, S1-Ks1.0, S2-Ks0.5, S3-Bedrock0.5 and S4-GW0.5 (added 2026-11-07), are one-change screens from Re-tuned (see below).", False,
     "comment above SETTING"),
    (CELL_BEDROCK, "bedrock_path = f'{SOLUS_DIR}/bedrock_depth_plus{BEDROCK_ADD_M:g}m.asc'",
     "bed_tag = f'{BEDROCK_ADD_M:g}'.replace('.', 'p')    # 0.5 -> 0p5, so the file name has no extra dot\n"
     "    bedrock_path = f'{SOLUS_DIR}/bedrock_depth_plus{bed_tag}m.asc'", True,
     "bedrock file name without an extra dot for 0.5 m"),
    (CELL_GW, "gw_name = f'InitGW_{int(round(GW_FRACTION * 100))}pct_plus{BEDROCK_ADD_M:g}m_mm.asc'",
     "bed_tag = f'{BEDROCK_ADD_M:g}'.replace('.', 'p')    # 0.5 -> 0p5, so the file name has no extra dot\n"
     "gw_name = f'InitGW_{int(round(GW_FRACTION * 100))}pct_plus{bed_tag}m_mm.asc'", True,
     "water table file name without an extra dot for 0.5 m"),
    (CELL_INTRO, "and two settings, `Re-tuned` and `Original`, chosen by",
     "and named settings (`Re-tuned`, `Original`, the two Kchan runs, and the one-change screens `S1-Ks1.0`, `S2-Ks0.5`, `S3-Bedrock0.5` and `S4-GW0.5`), chosen by", False,
     "intro text"),
    (CELL_SETTINGS_MD, "`SETTING` picks Re-tuned or Original; the multipliers, anisotropy, input file name and results name all follow from it.",
     "`SETTING` picks the named setting (Re-tuned, Original, the two Kchan runs, or a one-change screen S1 to S4); the multipliers, anisotropy, "
     "bedrock depth, starting water table, channel loss, input file name and results name all follow from it. The screens are probes, not "
     "candidate settings: agree how you will read their results before you open them.", False,
     "Settings notes"),
    (CELL_SETTINGS_MD, "off for Re-tuned and Original, on for Original-Kchan0.5",
     "off for Re-tuned, Original and the S1 to S4 screens, on for Original-Kchan0.5", False,
     "Settings notes (channel loss line)"),
]

# The values compared between settings.
KEYS = ['KS_MULT', 'F_MULT', 'ANISOTROPY', 'ANISOTROPY_U', 'BEDROCK_ADD_M', 'GW_FRACTION', 'OPTPERCOLATION', 'CHANNELCONDUCTIVITY']
EXISTING = ['Re-tuned', 'Original', 'Original-Kchan0.5', 'Original-Kchan0.2']
SCREENS = ['S1-Ks1.0', 'S2-Ks0.5', 'S3-Bedrock0.5', 'S4-GW0.5']


def cell_by_id(nb, cid):
    hits = [c for c in nb['cells'] if c.get('id') == cid]
    if len(hits) != 1:
        raise SystemExit(f'STOP: expected exactly one cell with id {cid}, found {len(hits)}. Nothing was written.')
    return hits[0]


def get_text(cell):
    return ''.join(cell['source'])


def set_text(cell, text):
    """Store text the way Jupyter does: a list of lines, each ending in a newline except the last."""
    parts = text.split('\n')
    lines = [p + '\n' for p in parts[:-1]]
    if parts[-1] != '':
        lines.append(parts[-1])
    cell['source'] = lines


def settings_table(code):
    """Run ONLY the Settings cell (it needs nothing but os) once per setting; return {name: values}. No files are read."""
    out = {}
    names = EXISTING + [n for n in SCREENS]
    for name in names:
        os.environ['SYCAMORE_SETTING'] = name
        ns = {'os': os}
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf):
                exec(compile(code, '<settings cell>', 'exec'), ns)
        except ValueError:
            out[name] = None            # the setting does not exist in this version of the cell
            continue
        out[name] = {k: ns[k] for k in KEYS}
        out[name]['INPUT_FILE'] = ns['INPUT_FILE']
        out[name]['RESULTS_NAME'] = ns['RESULTS_NAME']
    os.environ.pop('SYCAMORE_SETTING', None)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--in', dest='src', default='Make_Sycamore_Model_WY2020.ipynb', help='the notebook to read (never changed)')
    ap.add_argument('--out', dest='dst', default='Make_Sycamore_Model_WY2020_screens.ipynb', help='the new notebook to write')
    ap.add_argument('--force', action='store_true', help='allow overwriting the output file')
    args = ap.parse_args()

    if not os.path.isfile(args.src):
        sys.exit(f'STOP: {args.src} not found. Run this from workspaces/Sycamore/lab/, or give the path with --in.')
    if os.path.exists(args.dst) and not args.force:
        sys.exit(f'STOP: {args.dst} already exists. Use --out for another name, or --force to overwrite it.')
    if os.path.abspath(args.src) == os.path.abspath(args.dst):
        sys.exit('STOP: the input and output are the same file. This script never overwrites the original.')

    with open(args.src, encoding='utf-8') as fh:
        nb = json.load(fh)

    # ---- before: what the existing settings give today -------------------------------------------------------
    before = settings_table(get_text(cell_by_id(nb, CELL_SETTINGS)))
    for name in EXISTING:
        if before[name] is None:
            sys.exit(f"STOP: the setting '{name}' is not in this notebook, so it is not the copy I read. Nothing was written.")
    if any(before[n] is not None for n in SCREENS):
        sys.exit('STOP: this notebook already has one of the S1 to S4 settings. Nothing was written.')

    # ---- edits ----------------------------------------------------------------------------------------------
    print('Edits:')
    problems = []
    for cid, old, new, required, what in EDITS:
        cell = cell_by_id(nb, cid)
        text = get_text(cell)
        n = text.count(old)
        if n != 1:
            msg = f'cell {cid}: looked for the text for "{what}" and found it {n} times (needed exactly 1)'
            if required:
                problems.append(msg)
                print('  FAILED  ' + msg)
            else:
                print('  skipped ' + msg + ' (wording only; the patch still works)')
            continue
        set_text(cell, text.replace(old, new))
        print(f'  done    cell {cid}: {what}')
    if problems:
        sys.exit('\nSTOP: the notebook is not the copy I read, so nothing was written. Send me the lines marked FAILED.')

    # ---- checks ---------------------------------------------------------------------------------------------
    for cid in (CELL_SETTINGS, CELL_BEDROCK, CELL_GW):
        ast.parse(get_text(cell_by_id(nb, cid)))          # raises SyntaxError if an edit broke a cell
    print('\nEvery edited code cell still parses.')

    after = settings_table(get_text(cell_by_id(nb, CELL_SETTINGS)))
    failures = []
    for name in EXISTING:
        if after[name] != before[name]:
            failures.append(f'{name} changed: before {before[name]} | after {after[name]}')
    base = after['Re-tuned']
    for name in SCREENS:
        if after[name] is None:
            failures.append(f'{name} is missing after the patch')
            continue
        diff = [k for k in KEYS if after[name][k] != base[k]]
        if len(diff) != 1:
            failures.append(f'{name} differs from Re-tuned in {len(diff)} values {diff} (needed exactly 1)')
    names_out = [after[n]['INPUT_FILE'] for n in EXISTING + SCREENS if after[n]]
    if len(set(names_out)) != len(names_out):
        failures.append('two settings would write the same input file name')

    print('\nSettings after the patch (a * marks a value that differs from Re-tuned):')
    width = max(len(n) for n in EXISTING + SCREENS)
    print(f"{'setting':{width}s}  " + '  '.join(f'{k:>13s}' for k in KEYS))
    for name in EXISTING + SCREENS:
        row = after[name]
        if row is None:
            continue
        cells = []
        for k in KEYS:
            mark = '*' if row[k] != base[k] else ' '
            cells.append(f'{str(row[k]):>12s}{mark}')
        print(f'{name:{width}s}  ' + '  '.join(cells))
    print('\nInput files the settings write:')
    for name in EXISTING + SCREENS:
        print(f'  {name:{width}s} -> {after[name]["INPUT_FILE"]}')

    if failures:
        print('\nCHECKS FAILED, nothing was written:')
        for f in failures:
            print('  - ' + f)
        sys.exit(1)
    print('\nChecks passed: the four existing settings are unchanged, and each screen differs from Re-tuned in exactly one value.')

    # ---- write ------------------------------------------------------------------------------------------------
    try:
        import nbformat
        nbformat.validate(nbformat.from_dict(nb))
        print('nbformat says the notebook is valid.')
    except ImportError:
        print('(nbformat is not installed, so the notebook structure was not validated that way.)')
    with open(args.dst, 'w', encoding='utf-8') as fh:
        json.dump(nb, fh, indent=1, ensure_ascii=False, sort_keys=True)
        fh.write('\n')
    print(f'\nWrote {args.dst}. Your original {args.src} was not changed.')


if __name__ == '__main__':
    main()
