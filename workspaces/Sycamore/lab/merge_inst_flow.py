#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
merge_inst_flow.py
==================

Join the pieces of a USGS instantaneous-value (15-minute) discharge download
into ONE tidy CSV file, and report anything that looks wrong.

WHY THIS EXISTS
---------------
The USGS web service sometimes refuses a large request (HTTP 503), so a year of
15-minute data may have been downloaded in pieces: one file per quarter, month
or day. This script stitches the pieces back together in time order.

WHAT IT DOES
------------
1. READS every file matching --glob in --indir (default iv_*.rdb). Files that
   are really web error pages are skipped and named.
2. JOINS them, sorts by time and drops exact duplicate readings (for example
   where two pieces overlap). If two copies of the same time disagree, it keeps
   the first and warns.
3. CHECKS
     * time zone code (tz_cd) is what you expect (Arizona: MST, no daylight saving)
     * readings that are not numbers (USGS writes things like Eqp, Ice, Dis, Mnt
       in the value column when the gauge could not measure) are kept as blanks
       and counted
     * quality codes (A approved, e estimated, P provisional)
     * GAPS: any stretch longer than --gap-minutes between valid readings
     * days with fewer readings than expected (--min-per-day)
4. CONVERTS cubic feet per second (cfs) to cubic metres per second (m3/s).
5. WRITES
     <out>              one CSV, one row per reading
     <out>_gaps.csv     one row per gap
6. OPTIONAL cross-check (--daily-rdb): recomputes daily means from the joined
   readings and compares them with the USGS daily-mean file. Small differences
   are normal (USGS averages the readings its own way); big ones are worth a look.

It only READS your downloaded files. It never changes or deletes them.

QUICK START
-----------
    cd workspaces/Sycamore/lab
    python merge_inst_flow.py
    python merge_inst_flow.py --daily-rdb ../init_data/obs/sycamore_09510200_daily_flow_WY2016_2026.rdb

OUTPUT COLUMNS
--------------
    datetime_local   the clock time exactly as USGS wrote it (Arizona: MST)
    tz_cd            time zone code from the file
    q_cfs            discharge, cubic feet per second (blank if not a number)
    q_cms            discharge, cubic metres per second
    qual_code        A, A:e, P, P:e ... as USGS wrote it
    flag             the text USGS wrote instead of a number (Eqp, Ice ...), else blank
    source_file      which downloaded piece the row came from

REQUIRES: numpy, pandas.

VERSION HISTORY
---------------
v1.0 (2026-10-09)  first version. Tested on stand-in files shaped like the USGS
                   RDB format; NOT yet run on the real WY2017 pieces.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_VERSION = "1.0"
CFS_TO_CMS = 0.028316846592      # 1 cubic foot = 0.028316846592 cubic metres (exact)

USER_SETTINGS = {
    "indir": "../init_data/obs",
    "glob": "iv_*.rdb",
    "out": "../init_data/obs/sycamore_09510200_inst_flow_WY2017.csv",
    "start": "2016-10-01 00:00",      # first moment you expect data (local clock)
    "end": "2017-09-30 23:59",        # last moment you expect data
    "expect_tz": "MST",
    "gap_minutes": 30.0,              # a stretch longer than this counts as a gap
    "min_per_day": 90,                # normal day = 96 readings (every 15 min)
    "daily_rdb": None,
    "daily_tol": 5.0,                 # percent; bigger daily differences get listed
}


# =============================================================================
# reading
# =============================================================================
def read_piece(path: Path):
    """Read one USGS RDB file. Returns (DataFrame, None) or (None, reason).

    An RDB file is tab-separated text: comment lines starting with '#', one
    header line starting with 'agency_cd', one line describing column widths
    (starts with '5s'), then the data rows.
    """
    text = path.read_text(encoding="utf-8", errors="replace")
    if text.lstrip()[:20].lower().startswith(("<!doctype", "<html")):
        return None, "looks like a web error page, not data"
    lines = text.splitlines()
    header_i = next((i for i, ln in enumerate(lines) if ln.startswith("agency_cd")), None)
    if header_i is None:
        return None, "no header line starting with 'agency_cd'"
    cols = lines[header_i].split("\t")
    if "datetime" not in cols:
        return None, "no 'datetime' column"
    val_idx = [i for i, c in enumerate(cols) if re.fullmatch(r"\d+_00060", c)]
    cd_idx = [i for i, c in enumerate(cols) if re.fullmatch(r"\d+_00060_cd", c)]
    if not val_idx:
        return None, "no discharge (parameter 00060) column"
    i_dt, i_val = cols.index("datetime"), val_idx[0]
    i_cd = cd_idx[0] if cd_idx else None
    i_tz = cols.index("tz_cd") if "tz_cd" in cols else None
    need = max(i for i in (i_dt, i_val, i_cd, i_tz) if i is not None)

    rec = []
    for ln in lines[header_i + 1:]:
        if not ln.strip() or ln.startswith("#"):
            continue
        p = ln.split("\t")
        if not p[0].isalpha() or len(p) <= need:      # skips the '5s 15s ...' line
            continue
        rec.append((p[i_dt].strip(),
                    p[i_tz].strip() if i_tz is not None else "",
                    p[i_val].strip(),
                    p[i_cd].strip() if i_cd is not None else ""))
    if not rec:
        return None, "no data rows"
    df = pd.DataFrame(rec, columns=["dt_text", "tz_cd", "value_text", "qual_code"])
    df["datetime_local"] = pd.to_datetime(df["dt_text"], errors="coerce")
    df["source_file"] = path.name
    df.attrs["n_value_cols"] = len(val_idx)
    return df, None


def read_daily_rdb(path: Path) -> pd.Series:
    """USGS daily-mean RDB file -> Series of cfs indexed by date."""
    out = {}
    for ln in Path(path).read_text(errors="replace").splitlines():
        if ln.startswith("#") or ln.startswith("agency_cd"):
            continue
        p = ln.split("\t")
        if len(p) >= 4 and p[0].isalpha():
            try:
                out[pd.Timestamp(p[2])] = float(p[3])
            except ValueError:
                pass
    return pd.Series(out, dtype="float64").sort_index()


# =============================================================================
# checks
# =============================================================================
class Report:
    """Collects PASS / WARN / INFO lines so they print together."""
    def __init__(self):
        self.lines, self.n_warn = [], 0

    def add(self, level: str, msg: str):
        self.lines.append(f"[{level:<4}] {msg}")
        if level == "WARN":
            self.n_warn += 1

    def blank(self):
        self.lines.append("")


def find_gaps(times: pd.Series, start: pd.Timestamp, end: pd.Timestamp, gap_min: float):
    """Gaps = stretches longer than gap_min between consecutive valid readings.

    The window edges count too: if the first reading comes long after `start`,
    that is a gap.
    """
    t = pd.Series([start, *times.tolist(), end])
    d = t.diff().dt.total_seconds() / 60.0
    rows = []
    for i in np.flatnonzero(d.values > gap_min):
        a, b = t.iloc[i - 1], t.iloc[i]
        rows.append({"last_reading_before": a, "next_reading_after": b,
                     "gap_hours": round(d.iloc[i] / 60.0, 2),
                     "approx_missing_15min_readings": int(round(d.iloc[i] / 15.0)) - 1})
    return pd.DataFrame(rows, columns=["last_reading_before", "next_reading_after",
                                       "gap_hours", "approx_missing_15min_readings"])


def daily_means_from_readings(times: pd.Series, values: pd.Series, gap_min: float,
                              min_valid_minutes: int = 1380) -> pd.Series:
    """Daily mean (cfs) rebuilt from the readings.

    Readings are joined by straight lines to a 1-minute grid. A minute only
    counts if the readings on both sides of it are at most gap_min apart. Days
    with fewer than min_valid_minutes counted minutes (of 1440) are left out.
    """
    t = times.values.astype("datetime64[m]").astype("int64")
    v = values.values.astype("float64")
    m0, m1 = (t[0] // 1440) * 1440, (t[-1] // 1440 + 1) * 1440
    g = np.arange(m0, m1)
    q = np.interp(g, t, v)
    p = np.searchsorted(t, g, side="right")
    prev_i, next_i = p - 1, p
    has_prev, has_next = prev_i >= 0, next_i < len(t)
    span = np.full(len(g), np.inf)
    both = has_prev & has_next
    span[both] = t[next_i[both]] - t[prev_i[both]]
    exact = has_prev & (t[np.clip(prev_i, 0, len(t) - 1)] == g)
    valid = exact | (both & (span <= gap_min))
    day = pd.to_datetime(g * 60, unit="s").normalize()
    df = pd.DataFrame({"day": day, "q": np.where(valid, q, np.nan)})
    grp = df.groupby("day")["q"]
    out = grp.mean()
    out[grp.count() < min_valid_minutes] = np.nan
    return out.dropna()


# =============================================================================
# main
# =============================================================================
def build_parser() -> argparse.ArgumentParser:
    U = USER_SETTINGS
    p = argparse.ArgumentParser(
        description="Join USGS 15-minute discharge pieces into one CSV (see the top of this file).",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--indir", default=U["indir"], help="folder holding the downloaded pieces")
    p.add_argument("--glob", default=U["glob"], help="which files in --indir to join")
    p.add_argument("--out", default=U["out"], help="output CSV path")
    p.add_argument("--start", default=U["start"], help="first expected moment, local clock")
    p.add_argument("--end", default=U["end"], help="last expected moment, local clock")
    p.add_argument("--expect-tz", default=U["expect_tz"], help="expected tz_cd (Arizona: MST)")
    p.add_argument("--gap-minutes", type=float, default=U["gap_minutes"])
    p.add_argument("--min-per-day", type=int, default=U["min_per_day"])
    p.add_argument("--daily-rdb", default=U["daily_rdb"],
                   help="USGS daily-mean RDB file, for the cross-check")
    p.add_argument("--daily-tol", type=float, default=U["daily_tol"],
                   help="percent difference above which a day is listed")
    return p


def main(argv=None) -> int:
    a = build_parser().parse_args(argv)
    rep = Report()
    indir = Path(a.indir)
    files = sorted(indir.glob(a.glob))
    print(f"merge_inst_flow.py v{SCRIPT_VERSION}")
    print(f"  looking in {indir.resolve()} for {a.glob}: {len(files)} files")
    if not files:
        print("ERROR: no files found. Check --indir and --glob.", file=sys.stderr)
        return 1

    # ---- 1. read ---------------------------------------------------------------
    parts, skipped = [], []
    for f in files:
        df, why = read_piece(f)
        if df is None:
            skipped.append((f.name, why))
        else:
            parts.append(df)
            if df.attrs.get("n_value_cols", 1) > 1:
                rep.add("WARN", f"{f.name} has more than one discharge column; used the first")
    for name, why in skipped:
        rep.add("WARN", f"skipped {name}: {why}")
    if not parts:
        print("\n".join(rep.lines))
        print("ERROR: no usable data in any file.", file=sys.stderr)
        return 1
    raw = pd.concat(parts, ignore_index=True)
    rep.add("INFO", f"read {len(parts)} files, {len(raw):,} rows"
                    + (f" ({len(skipped)} skipped)" if skipped else ""))

    bad_t = raw["datetime_local"].isna()
    if bad_t.any():
        rep.add("WARN", f"{int(bad_t.sum())} rows had a time that could not be read; dropped")
        raw = raw[~bad_t]

    # ---- 2. time zone ------------------------------------------------------------
    tzs = raw["tz_cd"].value_counts().to_dict()
    if set(tzs) <= {a.expect_tz}:
        rep.add("PASS", f"every row is in time zone {a.expect_tz}")
    else:
        odd = raw[raw["tz_cd"] != a.expect_tz]["source_file"].value_counts().to_dict()
        rep.add("WARN", f"time zone codes found: {tzs}; expected only {a.expect_tz}. "
                        f"Files with other codes: {odd}. Times were NOT converted.")

    # ---- 3. numbers, codes -----------------------------------------------------
    raw["q_cfs"] = pd.to_numeric(raw["value_text"], errors="coerce")
    nonnum = raw["q_cfs"].isna()
    raw["flag"] = np.where(nonnum, raw["value_text"], "")
    if nonnum.any():
        counts = raw.loc[nonnum, "flag"].replace("", "(blank)").value_counts().to_dict()
        rep.add("WARN", f"{int(nonnum.sum())} readings are not numbers and were left blank: {counts}")
    else:
        rep.add("PASS", "every reading is a number")
    neg = (raw["q_cfs"] < 0).sum()
    if neg:
        rep.add("WARN", f"{int(neg)} readings are negative")

    # ---- 4. sort, remove duplicates ------------------------------------------
    raw = raw.sort_values(["datetime_local"], kind="stable").reset_index(drop=True)
    dup = raw.duplicated("datetime_local", keep=False)
    if dup.any():
        first = raw[dup].groupby("datetime_local")
        disagree = sum(1 for _, g in first
                       if g["q_cfs"].nunique(dropna=False) > 1 or g["qual_code"].nunique() > 1)
        n_extra = int(raw.duplicated("datetime_local", keep="first").sum())
        rep.add("INFO", f"{n_extra} duplicate readings removed (same time in more than one piece)")
        if disagree:
            rep.add("WARN", f"{disagree} duplicated times had DIFFERENT values or codes; kept the first")
        raw = raw.drop_duplicates("datetime_local", keep="first").reset_index(drop=True)

    # ---- 5. window ------------------------------------------------------------
    start, end = pd.Timestamp(a.start), pd.Timestamp(a.end)
    outside = (raw["datetime_local"] < start) | (raw["datetime_local"] > end)
    if outside.any():
        rep.add("WARN", f"{int(outside.sum())} readings fall outside {a.start} to {a.end}; dropped")
        raw = raw[~outside].reset_index(drop=True)
    if raw.empty:
        print("\n".join(rep.lines))
        print("ERROR: no readings inside the requested window.", file=sys.stderr)
        return 1
    raw["q_cms"] = raw["q_cfs"] * CFS_TO_CMS

    # ---- 6. quality codes ---------------------------------------------------------
    codes = raw["qual_code"].replace("", "(blank)").value_counts().to_dict()
    rep.add("INFO", f"quality codes: {codes}")
    n_prov = int(raw["qual_code"].str.startswith("P").sum())
    if n_prov:
        rep.add("WARN", f"{n_prov} readings are provisional (P) and can still be revised")
    else:
        rep.add("PASS", "no provisional readings")
    n_est = int(raw["qual_code"].str.contains("e").sum())
    rep.add("INFO", f"{n_est} readings ({100 * n_est / len(raw):.1f}%) are marked estimated (e)")

    # ---- 7. gaps and thin days ------------------------------------------------
    valid = raw[raw["q_cfs"].notna()]
    gaps = find_gaps(valid["datetime_local"], start, end, a.gap_minutes)
    if gaps.empty:
        rep.add("PASS", f"no gaps longer than {a.gap_minutes:g} minutes")
    else:
        tot_h = gaps["gap_hours"].sum()
        rep.add("WARN", f"{len(gaps)} gaps longer than {a.gap_minutes:g} minutes, "
                        f"{tot_h:.1f} hours in total (listed below and in the gaps file)")
        for r in gaps.itertuples():
            rep.add("INFO", f"   gap: {r.last_reading_before:%Y-%m-%d %H:%M} -> "
                            f"{r.next_reading_after:%Y-%m-%d %H:%M}   ({r.gap_hours:g} h)")
    per_day = valid.groupby(valid["datetime_local"].dt.normalize()).size()
    all_days = pd.date_range(start.normalize(), end.normalize(), freq="D")
    per_day = per_day.reindex(all_days, fill_value=0)
    thin = per_day[per_day < a.min_per_day]
    if thin.empty:
        rep.add("PASS", f"every day has at least {a.min_per_day} readings")
    else:
        # Show consecutive thin days as ranges to keep the report short.
        runs, run = [], [thin.index[0]]
        for d in thin.index[1:]:
            if (d - run[-1]).days == 1:
                run.append(d)
            else:
                runs.append(run)
                run = [d]
        runs.append(run)
        txt = "; ".join(f"{r[0]:%Y-%m-%d}" + (f" to {r[-1]:%Y-%m-%d}" if len(r) > 1 else "")
                        + f" (fewest {int(thin[r].min())})" for r in runs)
        rep.add("WARN", f"{len(thin)} days have fewer than {a.min_per_day} readings: {txt}")

    # ---- 8. summary numbers ------------------------------------------------------
    rep.blank()
    pk = valid.loc[valid["q_cfs"].idxmax()]
    rep.add("INFO", f"largest reading: {pk['q_cfs']:,.0f} cfs ({pk['q_cms']:,.1f} m3/s) at "
                    f"{pk['datetime_local']:%Y-%m-%d %H:%M} {pk['tz_cd']}")
    rep.add("INFO", f"first reading {valid['datetime_local'].iloc[0]:%Y-%m-%d %H:%M}, "
                    f"last {valid['datetime_local'].iloc[-1]:%Y-%m-%d %H:%M}; "
                    f"{len(valid):,} valid readings")
    mon = valid.groupby(valid["datetime_local"].dt.to_period("M"))["q_cfs"].agg(
        ["count", "max", "mean"])
    rep.add("INFO", "by month (readings, largest cfs, mean cfs):")
    for m, r in mon.iterrows():
        rep.add("INFO", f"   {m}: {int(r['count']):>6,}  {r['max']:>9,.1f}  {r['mean']:>8,.2f}")

    # ---- 9. optional cross-check against the daily file ------------------------------
    if a.daily_rdb:
        rep.blank()
        dpath = Path(a.daily_rdb)
        if not dpath.exists():
            rep.add("WARN", f"daily file not found: {dpath}")
        else:
            daily = read_daily_rdb(dpath)
            mine = daily_means_from_readings(valid["datetime_local"], valid["q_cfs"], a.gap_minutes)
            j = pd.concat([mine.rename("mine"), daily.rename("usgs")], axis=1, join="inner").dropna()
            if j.empty:
                rep.add("WARN", "no days to compare (check dates in the daily file)")
            else:
                j["pct_diff"] = 100 * (j["mine"] - j["usgs"]) / j["usgs"].where(j["usgs"] > 0.5)
                ok = (j["pct_diff"].abs() <= a.daily_tol) | j["pct_diff"].isna()
                rep.add("INFO", f"daily cross-check: {len(j)} days compared (days with big gaps are left out); "
                                f"{int(ok.sum())} agree within {a.daily_tol:g}% (days with USGS mean below 0.5 cfs not scored)")
                worst = j[~ok].reindex(j[~ok]["pct_diff"].abs().sort_values(ascending=False).index).head(10)
                if worst.empty:
                    rep.add("PASS", "all compared days agree with the USGS daily means")
                else:
                    rep.add("WARN", f"{int((~ok).sum())} days differ by more than {a.daily_tol:g}%; the largest:")
                    for d, r in worst.iterrows():
                        rep.add("INFO", f"   {d:%Y-%m-%d}: mine {r['mine']:.2f} cfs, USGS {r['usgs']:.2f} cfs "
                                        f"({r['pct_diff']:+.1f}%)")

    # ---- write ------------------------------------------------------------------
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    table = raw[["datetime_local", "tz_cd", "q_cfs", "q_cms", "qual_code", "flag", "source_file"]].copy()
    table["datetime_local"] = table["datetime_local"].dt.strftime("%Y-%m-%d %H:%M")
    table["q_cms"] = table["q_cms"].round(4)
    table.to_csv(out, index=False)
    gap_path = out.with_name(out.stem + "_gaps.csv")
    gaps_out = gaps.copy()
    for c in ("last_reading_before", "next_reading_after"):
        gaps_out[c] = gaps_out[c].dt.strftime("%Y-%m-%d %H:%M")
    gaps_out.to_csv(gap_path, index=False)

    rep.blank()
    verdict = "PASS" if rep.n_warn == 0 else f"{rep.n_warn} warning(s) - read them above"
    rep.lines.append(f"OVERALL: {verdict}")
    print()
    print("\n".join(rep.lines))
    print(f"\n  wrote {len(table):,} rows to {out}")
    print(f"  wrote {len(gaps)} gaps to {gap_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
