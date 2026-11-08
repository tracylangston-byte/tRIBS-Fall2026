"""
add_sep2019_gauge_sycamore.py  (v1.0)  --  put the September 2019 USGS readings in front of your WY2020 gauge file

WHAT IT DOES
  Your gauge file (sycamore_09510200_obs_flow_WY2020.csv) starts Oct 1 2019, but every WY2020 model run starts Sep 1 2019. This script reads a
  USGS download of the same gauge (09510200, discharge, "instantaneous values") that covers September 2019, converts it to the format of your
  file, and writes a NEW file that is September 2019 followed by your file, row for row.

  Same columns as your file:   datetime,q_cms,q_cfs,tier,note
    datetime  Arizona time (MST = UTC - 7 hours, no daylight saving), written as YYYY-MM-DD HH:MM
    q_cfs     the USGS value, cubic feet per second
    q_cms     q_cfs x 0.028316846592 (exactly 0.3048 cubed), rounded to 6 decimals (cubic meters per second)
    tier      'measured' (every row of the USGS file is approved data with no qualifier; the script stops if it is not)
    note      empty

  The USGS rows are kept at their own times (every 15 minutes, every 5 minutes during the storm on Sep 23), the same way your file keeps its
  own times. Nothing is filled in or smoothed. If the USGS file has no reading for a time, that time is simply not in the new file (your file
  has such holes too). The script lists the holes it finds.

THE TEST OF THE TIME CONVERSION
  The USGS download usually ends a little past Oct 1, so it overlaps the first day of your file. The script compares every overlapping USGS
  reading with the row at the same time in your file (flow in CFS and in m3/s). If even one differs, it stops and writes nothing. If
  the clock conversion were off by an hour, they would not match, so a match means the September times are right.
  If the two files do not overlap at all, it also stops, because then nothing proves the conversion.

YOUR FILES ARE NEVER CHANGED
  Your gauge file and the USGS file are only read. The new file is written under the name you give; the script will not overwrite an existing
  file unless you add --force, and will not write over either input. Your rows are copied into the new file as the exact same bytes.

RUN IT (from anywhere; give the paths):
    python add_sep2019_gauge_sycamore.py --usgs <USGS download>.csv --existing <your gauge file>.csv --out sycamore_09510200_obs_flow_Sep2019_to_May2020.csv
  Then give the new file to the scripts that take --obs:
    python score_wy2020_run_sycamore.py --obs sycamore_09510200_obs_flow_Sep2019_to_May2020.csv ...
    python edge_check_wy2020_sycamore.py --obs sycamore_09510200_obs_flow_Sep2019_to_May2020.csv ...

THE USGS FILE MUST HAVE these columns (the "primary-time-series" download has them): monitoring_location_id, parameter_code, statistic_id, time,
  value, unit_of_measure, approval_status, qualifier. Times must end in +00:00 (UTC). It may be in any time order.

WORDS
  UTC                the world clock (Arizona time + 7 hours).         instantaneous value   one reading at one moment, not a daily average.
  parameter 00060    discharge (flow).                                   statistic 00011       instantaneous values.
  CFS                cubic feet per second (m3/s x 35.3147).             overlap               the times that are in both files.

VERSION HISTORY
  v1.0 (2026-11-07)  first version.
"""
import argparse
import math
import os
import sys

import numpy as np
import pandas as pd

CFS_TO_CMS = 0.028316846592                           # exactly 0.3048 cubed (1 foot = 0.3048 m); your file's q_cms follows this factor to the last digit
AZ_OFFSET = pd.Timedelta(hours=7)                     # Arizona = UTC - 7, all year
HEADER = 'datetime,q_cms,q_cfs,tier,note'
NEED_USGS = ['monitoring_location_id', 'parameter_code', 'statistic_id', 'time', 'value', 'unit_of_measure', 'approval_status', 'qualifier']
GAP_MIN = pd.Timedelta(minutes=16)                    # a step longer than this is reported as a hole (readings are 15 or 5 minutes apart)


def stop(msg):
    sys.exit('STOP: ' + msg + '\nNothing was written.')


def fmt_cfs(v):
    s = repr(float(v))
    if 'e' in s or 'n' in s:
        stop('a flow value ({}) would be written in scientific notation or as nan; this script cannot format it.'.format(s))
    return s


def fmt_cms(cfs):
    s = repr(round(float(cfs) * CFS_TO_CMS, 6))
    if 'e' in s or 'n' in s:
        stop('a flow value ({} cfs) would be written in scientific notation or as nan; this script cannot format it.'.format(cfs))
    return s


def read_usgs(path):
    try:
        u = pd.read_csv(path, dtype=str, keep_default_na=False)
    except Exception as e:                                   # noqa: BLE001  (any read problem is reported the same way)
        stop('could not read {} as a CSV file ({}).'.format(path, e))
    miss = [c for c in NEED_USGS if c not in u.columns]
    if miss:
        stop('{} is missing column(s) {}. Columns found: {}'.format(path, miss, list(u.columns)[:15]))
    if len(u) == 0:
        stop('{} has no rows.'.format(path))
    one = {'monitoring_location_id': 'USGS-09510200', 'parameter_code': '00060', 'statistic_id': '00011', 'unit_of_measure': 'ft^3/s'}
    for col, want in one.items():
        vals = sorted(set(u[col]))
        if vals != [want]:
            stop('column {} must be {} on every row (this script is for Sycamore Wash discharge, instantaneous, in ft^3/s); found {}.'.format(col, want, vals[:6]))
    appr = u['approval_status'].value_counts().to_dict()
    if set(appr) != {'Approved'}:
        stop('approval_status is not "Approved" on every row ({}). Your file labels rows measured/estimated by its own rules, '
             'so tell me how you want unapproved rows labeled before using them.'.format(appr))
    if (u['qualifier'].str.strip() != '').any():
        q = u['qualifier'].value_counts().to_dict()
        stop('some rows have a qualifier ({}). Tell me how you want those rows labeled before using them.'.format(q))
    bad_t = ~u['time'].str.endswith('+00:00')
    if bad_t.any():
        stop('{} time(s) do not end in +00:00 (UTC), for example {!r}. This script converts from UTC only.'.format(int(bad_t.sum()), u.loc[bad_t, 'time'].iloc[0]))
    try:
        t_utc = pd.to_datetime(u['time'], utc=True)
    except Exception as e:                                   # noqa: BLE001
        stop('could not read the times in {} ({}).'.format(path, e))
    v = pd.to_numeric(u['value'], errors='coerce')
    if v.isna().any():
        stop('{} value(s) are empty or not numbers, for example {!r} at {}.'.format(int(v.isna().sum()), u.loc[v.isna(), 'value'].iloc[0], u.loc[v.isna(), 'time'].iloc[0]))
    if (~np.isfinite(v)).any() or (v < 0).any():
        stop('some flow values are negative or not finite (first: {} at {}).'.format(v[(~np.isfinite(v)) | (v < 0)].iloc[0],
                                                                                     u.loc[(~np.isfinite(v)) | (v < 0), 'time'].iloc[0]))
    az_exact = (t_utc.dt.tz_convert(None) - AZ_OFFSET)
    az = az_exact.dt.floor('min')
    n_floored = int((az != az_exact).sum())
    floored = [(str(a), str(b)) for a, b in zip(az_exact[az != az_exact], az[az != az_exact])]
    out = pd.DataFrame({'az': az.to_numpy(), 'cfs': v.to_numpy(float)}).sort_values('az').reset_index(drop=True)
    dup = out['az'].duplicated()
    if dup.any():
        stop('{} reading(s) share a minute with another reading (first: {}). I will not guess which to keep.'.format(int(dup.sum()), out.loc[dup, 'az'].iloc[0]))
    return out, n_floored, floored


def read_existing(path):
    raw = open(path, 'rb').read()
    if not raw:
        stop('{} is empty.'.format(path))
    eol = b'\r\n' if b'\r\n' in raw else b'\n'
    lines = raw.split(eol)
    if lines and lines[-1] == b'':
        lines = lines[:-1]
    head = lines[0].decode('utf-8-sig').strip()
    if head != HEADER:
        stop('the first line of {} is {!r}, not {!r}.'.format(path, head, HEADER))
    body = lines[1:]
    if not body:
        stop('{} has a header but no rows.'.format(path))
    rows = []
    for i, ln in enumerate(body, start=2):
        parts = ln.decode('utf-8').split(',')
        if len(parts) != 5:
            stop('line {} of {} does not have 5 columns: {!r}'.format(i, path, ln[:80]))
        rows.append(parts)
    d = pd.DataFrame(rows, columns=['datetime', 'q_cms', 'q_cfs', 'tier', 'note'])
    d['t'] = pd.to_datetime(d['datetime'], format='%Y-%m-%d %H:%M', errors='coerce')
    if d['t'].isna().any():
        stop('line {} of {} has a datetime I cannot read: {!r}'.format(int(np.flatnonzero(d['t'].isna())[0]) + 2, path, d['datetime'][d['t'].isna()].iloc[0]))
    if not d['t'].is_monotonic_increasing or d['t'].duplicated().any():
        stop('the times in {} are not in strictly increasing order.'.format(path))
    return eol, body, d


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--usgs', required=True, help='the USGS download (CSV, UTC, ft^3/s)')
    ap.add_argument('--existing', required=True, help='your gauge file (never changed)')
    ap.add_argument('--out', required=True, help='the new file to write')
    ap.add_argument('--force', action='store_true', help='allow overwriting the output file')
    a = ap.parse_args()

    for label, p in (('--usgs', a.usgs), ('--existing', a.existing)):
        if not os.path.isfile(p):
            stop('{} file not found: {}'.format(label, p))
    ap_out = os.path.abspath(a.out)
    if ap_out in (os.path.abspath(a.usgs), os.path.abspath(a.existing)):
        stop('--out is the same file as an input. This script never writes over an input.')
    if os.path.exists(a.out) and not a.force:
        stop('{} already exists. Use another --out, or add --force to overwrite it.'.format(a.out))

    print('add_sep2019_gauge_sycamore.py v1.0 (reads your two files; writes only the new file)')
    usgs, n_floor, floored = read_usgs(a.usgs)
    eol, body, ex = read_existing(a.existing)
    first_ex = ex['t'].iloc[0]
    print('  USGS file : {} readings, {} to {} (Arizona time)'.format(len(usgs), usgs['az'].iloc[0], usgs['az'].iloc[-1]))
    print('  your file : {} rows, {} to {}'.format(len(ex), first_ex, ex['t'].iloc[-1]))

    # ---- the overlap test --------------------------------------------------------------------------------
    ov = usgs[usgs['az'] >= first_ex].copy()
    if len(ov) == 0:
        stop('the USGS file does not reach your file (it ends {} and your file starts {}), so there is nothing to compare and the clock '
             'conversion cannot be checked. Download a few hours past {} too.'.format(usgs['az'].iloc[-1], first_ex, first_ex))
    ex_idx = ex.set_index('t')
    bad, same_text, n_cmp = [], 0, 0
    for t, cfs in zip(ov['az'], ov['cfs']):
        if t not in ex_idx.index:
            continue                                         # your file has no row at that time (counted below)
        n_cmp += 1
        row = ex_idx.loc[t]
        try:
            e_cfs, e_cms = float(row['q_cfs']), float(row['q_cms'])
        except ValueError:
            bad.append((t, cfs, row['q_cfs'], row['q_cms'], 'your row has no number'))
            continue
        if not (math.isclose(e_cfs, cfs, rel_tol=0, abs_tol=1e-6) and math.isclose(e_cms, round(cfs * CFS_TO_CMS, 6), rel_tol=0, abs_tol=5e-7)):
            bad.append((t, cfs, row['q_cfs'], row['q_cms'], 'differs'))
        elif row['q_cfs'] == fmt_cfs(cfs) and row['q_cms'] == fmt_cms(cfs):
            same_text += 1
    if n_cmp == 0:
        stop('none of the {} USGS readings from {} on has a row at the same minute in your file, so the clock conversion cannot be checked.'.format(len(ov), first_ex))
    if bad:
        lines = ['  {}  USGS {} cfs | your file {} cfs, {} cms ({})'.format(t, c, ec, em, why) for t, c, ec, em, why in bad[:10]]
        stop('{} of the {} overlapping readings differ from your file, so the September times or units may be wrong. First differences:\n{}'.format(
            len(bad), n_cmp, '\n'.join(lines)))
    only_usgs = int(sum(1 for t in ov['az'] if t not in ex_idx.index))
    span_end = ov['az'].iloc[-1]
    only_ex = int(((ex['t'] >= first_ex) & (ex['t'] <= span_end) & ~ex['t'].isin(ov['az'])).sum())
    print('  overlap   : {} readings from {} to {}; all {} that exist in both files match your file in CFS and m3/s ({} also match the text exactly)'.format(
        len(ov), ov['az'].iloc[0], span_end, n_cmp, same_text))
    if only_usgs or only_ex:
        print('              ({} USGS reading(s) have no row in your file and {} row(s) of your file have no USGS reading in that span; your rows are kept as they are)'.format(
            only_usgs, only_ex))

    # ---- the new rows ------------------------------------------------------------------------------------
    new = usgs[usgs['az'] < first_ex].reset_index(drop=True)
    if len(new) == 0:
        stop('the USGS file has no readings before {}, the start of your file, so there is nothing to add.'.format(first_ex))
    junction = first_ex - new['az'].iloc[-1]
    new_lines = ['{:%Y-%m-%d %H:%M},{},{},measured,'.format(t, fmt_cms(c), fmt_cfs(c)).encode('utf-8') for t, c in zip(new['az'], new['cfs'])]
    print('  new rows  : {} rows, {} to {}'.format(len(new), new['az'].iloc[0], new['az'].iloc[-1]))
    if junction > pd.Timedelta(minutes=60):
        print('  WARNING   : the last new row is {} before the first row of your file (a hole at the join).'.format(junction))
    steps = new['az'].diff()
    holes = [(new['az'].iloc[i - 1], new['az'].iloc[i]) for i in np.flatnonzero((steps > GAP_MIN).to_numpy())]
    if holes:
        print('  holes     : {} place(s) where the USGS file has no reading for more than 15 minutes (nothing was filled in):'.format(len(holes)))
        for lo, hi in holes[:20]:
            print('              no readings between {} and {} ({:g} minutes apart)'.format(lo, hi, (hi - lo) / pd.Timedelta(minutes=1)))
    else:
        print('  holes     : none')
    if n_floor:
        print('  seconds   : {} USGS time(s) had seconds; your file has none, so the seconds were dropped:'.format(n_floor))
        for ex_t, fl in floored[:10]:
            print('              {} -> {}'.format(ex_t, fl))
    pk = new.loc[new['cfs'].idxmax()]
    h = new.set_index('az')['cfs'].mul(CFS_TO_CMS).resample('60min').mean().dropna()
    print('  check     : peak {:g} cfs ({:.1f} m3/s) at {};  hourly-mean volume {:.3f} million m3;  {} readings above 1 cfs'.format(
        pk['cfs'], pk['cfs'] * CFS_TO_CMS, pk['az'], h.sum() * 3600 / 1e6, int((new['cfs'] > 1).sum())))

    # ---- write -------------------------------------------------------------------------------------------
    nl = eol
    payload = nl.join([HEADER.encode()] + new_lines + body) + nl
    with open(a.out, 'wb') as fh:
        fh.write(payload)

    # ---- read it back --------------------------------------------------------------------------------------
    try:
        eol2, body2, out = read_existing(a.out)
    except SystemExit as e:                                  # the new file cannot even be read back
        try:
            os.remove(a.out)
        except OSError:
            pass
        sys.exit('{}\nI removed the new file.'.format(e))
    problems = []
    if len(body2) != len(new) + len(body):
        problems.append('row count {} is not {} + {}'.format(len(body2), len(new), len(body)))
    if body2[len(new):] != body:
        problems.append('your rows are not copied as the same bytes')
    if not out['t'].is_monotonic_increasing or out['t'].duplicated().any():
        problems.append('times are not strictly increasing')
    back = out.iloc[:len(new)]
    if not (back['t'].to_numpy() == new['az'].to_numpy()).all() or not np.allclose(back['q_cfs'].astype(float), new['cfs']):
        problems.append('new rows do not read back the same')
    if problems:
        try:
            os.remove(a.out)
        except OSError:
            pass
        stop('the file I wrote did not pass the read-back check ({}); I removed it.'.format('; '.join(problems)))
    print('\nWrote {}: {} rows = {} new + {} of yours (copied byte for byte). Times strictly increase; read back OK.'.format(a.out, len(body2), len(new), len(body)))
    print('Your {} and {} were not changed.'.format(os.path.basename(a.existing), os.path.basename(a.usgs)))


if __name__ == '__main__':
    main()
