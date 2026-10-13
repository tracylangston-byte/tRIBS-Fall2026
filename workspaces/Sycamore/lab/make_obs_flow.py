#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_obs_flow.py
================

Turn the joined USGS 15-minute record into the ONE observed-flow file the model
will be scored against: gaps handled by written rules, every row tagged, units
in m3/s, and a printed check of the result.

INPUTS (only read, never changed)
    the joined file from merge_inst_flow.py   (--merged)
    the USGS daily-mean RDB file              (--daily-rdb)

WHAT IT DOES
------------
1. FINDS every gap: a stretch longer than --gap-minutes between valid readings.
2. DECIDES what to do with each gap by looking up the USGS daily mean for the
   days the gap covers. The first rule that fits wins:

   a. DRY DAY      daily mean is 0.00 on every gap day, and the flow was
                   (nearly) zero at both ends of the gap
                   -> fill with 0 cfs                      tier: filled_zero
   b. SHORT GAP    1 hour or less, and the flows at the two ends are within
                   25% of each other (or both under 1 cfs)
                   -> straight line                        tier: filled_interp
   c. RECESSION    USGS daily values for the gap days are approved (no 'e' or
                   'P'), the gap is 36 hours or less, and the daily means
                   change by 25% or less from one day to the next across the
                   gap (a smooth recession)
                   -> straight line                        tier: filled_interp
   d. EVERYTHING ELSE (storm days, estimated days, long gaps, no daily value)
                   -> leave blank                          tier: missing

   Why blank rather than guessed: anything invented for a storm day would
   become the thing the model is tuned to match.
3. TAGS every real reading: measured (quality code A), estimated (code with
   'e'), provisional (code with 'P').
4. CONVERTS cfs to m3/s (x 0.028316846592).
5. CHECKS the result and prints it: counts per tier, how each gap was handled,
   whether each straight-line day reproduces the USGS daily mean, the file's
   total volume against the USGS daily total (and how much of any difference
   the blank stretches explain), the peak, and the largest flow in each month.

OUTPUTS (written next to --out)
    <out>              one row per time slot:
                         datetime  q_cms  q_cfs  tier  note
                       Arizona time, no time-zone label (matches the model
                       clock); q_cms/q_cfs are blank where tier is 'missing'.
    <out stem>_fills.csv   one row per gap: what was decided and why

HOW THE MODEL CODE SHOULD USE IT
    obs = pd.read_csv(path, index_col="datetime", parse_dates=True)
    observed = obs["q_cms"]
    aligned = results.align_streamflow(observed, freq="15min", start=..., end=...)
    (leave interpolate OFF - it would draw a line across the blank stretches)

QUICK START
    cd workspaces/Sycamore/lab
    python make_obs_flow.py

REQUIRES: numpy, pandas.

VERSION HISTORY
v1.0 (2026-10-09)  first version. Tested on stand-in files with every case
                   planted; NOT yet run on the real WY2017 files.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_VERSION = "1.0"
CFS_TO_CMS = 0.028316846592      # 1 cubic foot = 0.028316846592 cubic metres (exact)
ZERO_CFS = 0.005                 # a USGS daily mean at or below this is "dry all day"

USER_SETTINGS = {
    "merged": "../init_data/obs/sycamore_09510200_inst_flow_WY2017.csv",
    "daily_rdb": "../init_data/obs/sycamore_09510200_daily_flow_WY2016_2026.rdb",
    "out": "../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv",
    "start": "2016-10-01 00:00",
    "end": "2017-09-30 23:59",
    "gap_minutes": 30.0,         # a stretch longer than this is a gap
    "slot_minutes": 15,          # blank/filled rows are written on this grid
    "smooth_pct": 25.0,          # rules b and c: how smooth counts as smooth (percent)
    "short_gap_hours": 1.0,      # rule b
    "max_interp_hours": 36.0,    # rule c
    "zero_edge_cfs": 0.5,        # rule a: flow at the gap ends must be at most this
    "day_tol": 5.0,              # percent: a straight-line day this far from USGS gets a warning
    "vol_tol": 2.0,              # percent: unexplained volume difference that gets a warning
}


# =============================================================================
# reading
# =============================================================================
def read_daily(path: Path) -> pd.DataFrame:
    """USGS daily-mean RDB -> DataFrame indexed by date with columns cfs, code."""
    rows = []
    for ln in Path(path).read_text(errors="replace").splitlines():
        if ln.startswith("#") or ln.startswith("agency_cd"):
            continue
        p = ln.split("\t")
        if len(p) >= 4 and p[0].isalpha():
            try:
                rows.append((pd.Timestamp(p[2]), float(p[3]), p[4].strip() if len(p) > 4 else ""))
            except ValueError:
                pass            # blank or non-numeric daily value: skipped
    if not rows:
        raise SystemExit(f"ERROR: no daily values found in {path}")
    return (pd.DataFrame(rows, columns=["date", "cfs", "code"])
            .drop_duplicates("date").set_index("date").sort_index())


def read_merged(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    need = {"datetime_local", "q_cfs", "qual_code"}
    if not need <= set(df.columns):
        raise SystemExit(f"ERROR: {path} is missing columns {sorted(need - set(df.columns))}. "
                         "Is it the output of merge_inst_flow.py?")
    df["datetime"] = pd.to_datetime(df["datetime_local"])
    df["q_cfs"] = pd.to_numeric(df["q_cfs"], errors="coerce")
    if "flag" not in df.columns:
        df["flag"] = ""
    return df.sort_values("datetime").drop_duplicates("datetime").reset_index(drop=True)


# =============================================================================
# gaps
# =============================================================================
def find_gaps(times, start, end, gap_min):
    """List of (a, b): last reading before the gap and first reading after it.
    a or b is None when the gap runs into the edge of the window."""
    t = list(times)
    if not t:
        return [(None, None)]
    out = []
    if (t[0] - start).total_seconds() / 60 > gap_min:
        out.append((None, t[0]))
    for x, y in zip(t[:-1], t[1:]):
        if (y - x).total_seconds() / 60 > gap_min:
            out.append((x, y))
    if (end - t[-1]).total_seconds() / 60 > gap_min:
        out.append((t[-1], None))
    return out


def gap_slots(a, b, start, end, slot_min):
    """Grid times (every slot_min minutes) that fall inside the gap."""
    lo = start if a is None else a
    hi = end if b is None else b
    freq = f"{slot_min}min"
    slots = pd.date_range(lo.ceil(freq), hi, freq=freq)
    return [s for s in slots
            if (a is None or s > a) and (b is None or s < b) and start <= s <= end]


def classify(a, b, qa, qb, days, daily, S):
    """Return (tier, reason) for one gap. See the top of the file for the rules."""
    dur_h = None if (a is None or b is None) else (b - a).total_seconds() / 3600
    present = [d for d in days if d in daily.index]
    means = [float(daily.loc[d, "cfs"]) for d in present]
    edges = [x for x in (qa, qb) if x is not None]

    # a. dry day
    zero_reason = None
    if len(present) == len(days) and all(m <= ZERO_CFS for m in means):
        if all(x <= S["zero_edge_cfs"] for x in edges):
            return "filled_zero", "USGS daily mean is 0.00 on every gap day and flow was ~0 at both ends"
        zero_reason = (f"USGS daily mean is 0.00 but flow was {max(edges):.2f} cfs at the edge "
                       "of the gap, so a zero fill is not safe")

    # b. short gap
    if dur_h is not None and dur_h <= S["short_gap_hours"]:
        hi, lo = max(qa, qb), min(qa, qb)
        if hi < 1.0 or (hi - lo) <= S["smooth_pct"] / 100 * hi:
            return "filled_interp", f"short gap ({dur_h:.2f} h) with similar flow at both ends"

    # c. smooth recession
    if len(present) < len(days):
        missing_days = ", ".join(f"{d:%Y-%m-%d}" for d in days if d not in daily.index)
        return "missing", f"no USGS daily value for {missing_days}"
    codes = [str(daily.loc[d, "code"]) for d in days]
    if any(("e" in c) or ("P" in c) for c in codes):
        return "missing", ("USGS daily values are estimated/provisional "
                           f"({', '.join(sorted(set(codes)))}) with mean flow "
                           f"{', '.join(f'{m:.1f}' for m in means)} cfs; not safe to draw a line")
    if dur_h is None:
        return "missing", "gap runs into the edge of the window"
    if dur_h > S["max_interp_hours"]:
        return "missing", f"gap of {dur_h:.1f} h is longer than the {S['max_interp_hours']:g} h limit for a straight line"
    prev_d, next_d = days[0] - pd.Timedelta(days=1), days[-1] + pd.Timedelta(days=1)
    if prev_d not in daily.index or next_d not in daily.index:
        return "missing", "no USGS daily value for the day before or after the gap"
    seq_days = [prev_d] + list(days) + [next_d]
    seq = [float(daily.loc[d, "cfs"]) for d in seq_days]
    smooth = True
    for i in range(len(seq) - 1):
        denom = seq[i + 1] if seq_days[i] == prev_d else seq[i]     # the gap day's mean
        if denom <= ZERO_CFS or abs(seq[i + 1] - seq[i]) > S["smooth_pct"] / 100 * denom:
            smooth = False
    if smooth:
        return "filled_interp", ("approved, smooth recession in the USGS daily means "
                                 f"({' -> '.join(f'{v:.1f}' for v in seq)} cfs)")
    why = zero_reason or ("daily means jump too much to be a smooth recession "
                          f"({' -> '.join(f'{v:.1f}' for v in seq)} cfs)")
    return "missing", why


# =============================================================================
# checks
# =============================================================================
def daily_means_from_series(times, values, gap_min, min_valid_minutes=1380):
    """Daily mean (cfs) rebuilt from a series: straight lines between readings,
    minutes inside a stretch longer than gap_min are not counted, days with fewer
    than min_valid_minutes counted minutes are left out."""
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


class Report:
    def __init__(self):
        self.lines, self.n_warn = [], 0

    def add(self, level, msg):
        self.lines.append(f"[{level:<4}] {msg}")
        if level == "WARN":
            self.n_warn += 1

    def blank(self):
        self.lines.append("")


# =============================================================================
# main
# =============================================================================
def build_parser():
    U = USER_SETTINGS
    p = argparse.ArgumentParser(
        description="Build the model-ready observed flow file (see the top of this file).",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--merged", default=U["merged"], help="joined 15-minute CSV from merge_inst_flow.py")
    p.add_argument("--daily-rdb", default=U["daily_rdb"], help="USGS daily-mean RDB file")
    p.add_argument("--out", default=U["out"], help="output CSV path")
    p.add_argument("--start", default=U["start"])
    p.add_argument("--end", default=U["end"])
    p.add_argument("--gap-minutes", type=float, default=U["gap_minutes"])
    p.add_argument("--slot-minutes", type=int, default=U["slot_minutes"])
    p.add_argument("--smooth-pct", type=float, default=U["smooth_pct"])
    p.add_argument("--short-gap-hours", type=float, default=U["short_gap_hours"])
    p.add_argument("--max-interp-hours", type=float, default=U["max_interp_hours"])
    p.add_argument("--zero-edge-cfs", type=float, default=U["zero_edge_cfs"])
    p.add_argument("--day-tol", type=float, default=U["day_tol"])
    p.add_argument("--vol-tol", type=float, default=U["vol_tol"])
    return p


def main(argv=None) -> int:
    a = build_parser().parse_args(argv)
    S = {"smooth_pct": a.smooth_pct, "short_gap_hours": a.short_gap_hours,
         "max_interp_hours": a.max_interp_hours, "zero_edge_cfs": a.zero_edge_cfs}
    rep = Report()
    start, end = pd.Timestamp(a.start), pd.Timestamp(a.end)

    print(f"make_obs_flow.py v{SCRIPT_VERSION}")
    for label, path in (("joined 15-minute file", a.merged), ("USGS daily file", a.daily_rdb)):
        if not Path(path).exists():
            print(f"ERROR: {label} not found: {path}", file=sys.stderr)
            return 1
    merged = read_merged(Path(a.merged))
    daily = read_daily(Path(a.daily_rdb))
    inwin = (merged["datetime"] >= start) & (merged["datetime"] <= end)
    merged = merged[inwin].reset_index(drop=True)
    valid = merged[merged["q_cfs"].notna()].reset_index(drop=True)
    if valid.empty:
        print("ERROR: no valid readings in the window.", file=sys.stderr)
        return 1
    rep.add("INFO", f"read {len(valid):,} valid readings from {Path(a.merged).name} "
                    f"and {len(daily):,} daily values from {Path(a.daily_rdb).name}")
    if "tz_cd" in merged.columns:
        tzs = sorted(set(merged["tz_cd"]))
        rep.add("PASS" if len(tzs) == 1 else "WARN", f"time zone codes in the joined file: {tzs}")

    # ---- real readings -> rows -----------------------------------------------------
    def tier_of(code):
        code = str(code)
        if "P" in code:
            return "provisional"
        return "estimated" if "e" in code else "measured"

    rows = {}
    for r in valid.itertuples():
        rows[r.datetime] = (r.q_cfs, tier_of(r.qual_code), "")
    for r in merged[merged["q_cfs"].isna()].itertuples():           # e.g. 'Eqp' readings
        rows[r.datetime] = (np.nan, "missing", f"USGS wrote '{r.flag}' instead of a number" if r.flag else "blank reading")

    # ---- gaps -----------------------------------------------------------------------
    gaps = find_gaps(valid["datetime"], start, end, a.gap_minutes)
    qmap = dict(zip(valid["datetime"], valid["q_cfs"]))
    fills = []
    for gi, (ga, gb) in enumerate(gaps, 1):
        slots = gap_slots(ga, gb, start, end, a.slot_minutes)
        if not slots:
            continue
        days = sorted({s.normalize() for s in slots})
        qa = qmap.get(ga) if ga is not None else None
        qb = qmap.get(gb) if gb is not None else None
        tier, reason = classify(ga, gb, qa, qb, days, daily, S)
        if tier == "filled_zero":
            vals = np.zeros(len(slots))
        elif tier == "filled_interp":
            xs = [ga.value, gb.value]
            vals = np.interp([s.value for s in slots], xs, [qa, qb])
        else:
            vals = np.full(len(slots), np.nan)
        for s, v in zip(slots, vals):
            rows[s] = (v, tier, f"gap {gi} (reason in the _fills.csv file)")
        fills.append({"gap": gi,
                      "last_reading_before": "" if ga is None else f"{ga:%Y-%m-%d %H:%M}",
                      "next_reading_after": "" if gb is None else f"{gb:%Y-%m-%d %H:%M}",
                      "hours": "" if (ga is None or gb is None) else round((gb - ga).total_seconds() / 3600, 2),
                      "slots": len(slots), "tier": tier, "reason": reason,
                      "_days": days, "_a": ga, "_b": gb})

    out = pd.DataFrame([(t, *v) for t, v in sorted(rows.items())],
                       columns=["datetime", "q_cfs", "tier", "note"])
    out["q_cms"] = out["q_cfs"] * CFS_TO_CMS

    # ---- report: each gap -------------------------------------------------------------
    rep.blank()
    if not fills:
        rep.add("PASS", "no gaps to handle")
    else:
        rep.add("INFO", f"{len(fills)} gaps found; what was decided:")
        for f in fills:
            span = (f"{f['last_reading_before'] or '(window start)'} -> {f['next_reading_after'] or '(window end)'}")
            rep.add("INFO", f"  gap {f['gap']}: {span}  ({f['hours']} h, {f['slots']} slots)")
            rep.add("INFO", f"      -> {f['tier']}: {f['reason']}")

    # ---- report: counts ----------------------------------------------------------------
    rep.blank()
    counts = out["tier"].value_counts()
    order = ["measured", "estimated", "provisional", "filled_zero", "filled_interp", "missing"]
    rep.add("INFO", f"{len(out):,} rows written; by tier:")
    for k in order:
        if k in counts:
            rep.add("INFO", f"   {k:<14} {int(counts[k]):>7,}")
    rep.add("PASS" if int(counts.sum()) == len(out) else "WARN",
            f"tier counts add up to the row total ({int(counts.sum()):,} of {len(out):,})")
    n_prov = int(counts.get("provisional", 0))
    if n_prov:
        rep.add("WARN", f"{n_prov} provisional readings (can still be revised by USGS)")
    n_miss = int(counts.get("missing", 0))
    if n_miss:
        miss = out[out["tier"] == "missing"]["datetime"]
        rep.add("INFO", f"{n_miss} rows are blank (tier 'missing'): first {miss.iloc[0]:%Y-%m-%d %H:%M}, "
                        f"last {miss.iloc[-1]:%Y-%m-%d %H:%M}. These times are not scored.")

    # ---- report: straight-line days against USGS ------------------------------------------
    filled_days = sorted({d for f in fills if f["tier"] == "filled_interp" for d in f["_days"]})
    pts = out.dropna(subset=["q_cfs"])
    if filled_days:
        rep.blank()
        dm = daily_means_from_series(pts["datetime"], pts["q_cfs"], a.gap_minutes)
        for d in filled_days:
            u = float(daily.loc[d, "cfs"]) if d in daily.index else np.nan
            if d in dm.index and not np.isnan(u) and u > 0:
                pct = 100 * (dm.loc[d] - u) / u
                lvl = "PASS" if abs(pct) <= a.day_tol else "WARN"
                rep.add(lvl, f"straight-line day {d:%Y-%m-%d}: file mean {dm.loc[d]:.2f} cfs vs USGS {u:.2f} cfs ({pct:+.1f}%)")
            else:
                rep.add("INFO", f"straight-line day {d:%Y-%m-%d}: could not be compared")

    # ---- report: volume --------------------------------------------------------------------
    rep.blank()
    t = pts["datetime"].values.astype("datetime64[s]").astype("int64")
    q_cms = pts["q_cms"].values
    dt = np.diff(t)
    ok = dt <= a.gap_minutes * 60
    file_vol = float(np.sum(0.5 * (q_cms[:-1] + q_cms[1:]) * dt * ok))
    days_all = pd.date_range(start.normalize(), end.normalize(), freq="D")
    dd = daily.reindex(days_all)
    usgs_vol = float(np.nansum(dd["cfs"].values) * 86400 * CFS_TO_CMS)
    n_nodaily = int(dd["cfs"].isna().sum())
    est_missing = 0.0
    for i in np.flatnonzero(~ok):                                 # blank stretches, valued at the USGS daily mean
        t0, t1 = pd.Timestamp(t[i], unit="s"), pd.Timestamp(t[i + 1], unit="s")
        cur = t0
        while cur < t1:
            nxt = min(t1, cur.normalize() + pd.Timedelta(days=1))
            if cur.normalize() in daily.index:
                est_missing += float(daily.loc[cur.normalize(), "cfs"]) * CFS_TO_CMS * (nxt - cur).total_seconds()
            cur = nxt
    ACRE_FT = 1233.48184
    rep.add("INFO", f"volume in the file:       {file_vol:>14,.0f} m3 ({file_vol / ACRE_FT:,.0f} acre-ft)")
    rep.add("INFO", f"USGS daily-mean total:    {usgs_vol:>14,.0f} m3 ({usgs_vol / ACRE_FT:,.0f} acre-ft)"
                    + (f"   [{n_nodaily} days have no daily value]" if n_nodaily else ""))
    short = usgs_vol - file_vol
    pct_short = 100 * short / usgs_vol
    rep.add("INFO", f"file is {abs(pct_short):.2f}% {'below' if pct_short >= 0 else 'above'} USGS; the blank stretches, "
                    f"valued at the USGS daily means, account for {100 * est_missing / usgs_vol:.2f}%")
    unexplained = 100 * (short - est_missing) / usgs_vol
    rep.add("PASS" if abs(unexplained) <= a.vol_tol else "WARN",
            f"unexplained volume difference: {unexplained:+.2f}% (limit {a.vol_tol:g}%)")

    # ---- report: peak and months -----------------------------------------------------------
    rep.blank()
    pk = pts.loc[pts["q_cfs"].idxmax()]
    rep.add("INFO", f"largest value: {pk['q_cfs']:,.0f} cfs ({pk['q_cms']:,.1f} m3/s) at "
                    f"{pk['datetime']:%Y-%m-%d %H:%M}, tier {pk['tier']}")
    mon = pts.groupby(pts["datetime"].dt.to_period("M"))["q_cfs"].max()
    rep.add("INFO", "largest flow in each month (cfs): "
                    + ", ".join(f"{m}: {v:,.1f}" for m, v in mon.items()))

    # ---- write ------------------------------------------------------------------------------
    outp = Path(a.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    table = pd.DataFrame({
        "datetime": out["datetime"].dt.strftime("%Y-%m-%d %H:%M"),
        "q_cms": out["q_cms"].round(6),
        "q_cfs": out["q_cfs"].round(4),
        "tier": out["tier"],
        "note": out["note"],
    })
    table.to_csv(outp, index=False)
    fills_path = outp.with_name(outp.stem + "_fills.csv")
    pd.DataFrame([{k: v for k, v in f.items() if not k.startswith("_")} for f in fills],
                 columns=["gap", "last_reading_before", "next_reading_after", "hours", "slots", "tier", "reason"]
                 ).to_csv(fills_path, index=False)

    rep.blank()
    rep.lines.append("OVERALL: " + ("PASS" if rep.n_warn == 0 else f"{rep.n_warn} warning(s) - read them above"))
    print()
    print("\n".join(rep.lines))
    print(f"\n  wrote {len(table):,} rows to {outp}")
    print(f"  wrote {len(fills)} gap decisions to {fills_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
