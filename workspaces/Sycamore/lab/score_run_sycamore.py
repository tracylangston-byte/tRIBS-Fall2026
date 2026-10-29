#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
score_run_sycamore.py   (v1.0, 2026-10-27)
==========================================

The standalone run scorer for the Lower Sycamore calibration design (v1). It scores ROUTED OUTLET FLOW
(results/<run>_Outlet.qout) against the USGS 09510200 gauge file, with every metric in the design's metrics table,
by window and by event, and writes ONE ROW PER RUN to a CSV. It runs no model and changes none of your files.

COMMANDS
--------
  score        score one or more runs and update the CSV           (this is also the RESCORE mode: it only reads kept files)
  fingerprint  Kchan fingerprints (Stage 7) for a pair: loss-off run vs loss-on run, or loss-off run vs the gauge
  dump-profile write the built-in WY2017 window list to a JSON file, to edit and reuse with --profile-file
  self-test    build synthetic runs with known answers and check every metric

  python score_run_sycamore.py score --run W5_ks0p2_cv6_nocl --window W5 --results-dir results
  python score_run_sycamore.py score --run A,B,C --window W3 --results-dir results
  python score_run_sycamore.py score --scan results --out-dir score_out      # every *_Outlet.qout under results
  python score_run_sycamore.py score --run W5_ks0p2_cv6_nocl --window W5 --results-dir results --dry-run
  python score_run_sycamore.py fingerprint --ref W3_base_nocl --test gauge --results-dir results --event flood_week
  python score_run_sycamore.py self-test

WHAT IT NEEDS FOR A RUN
-----------------------
  results/<run>_Outlet.qout   REQUIRED. comma-delimited, header line, columns Time_hr,Q (m3/s), Time_hr counted from the run start
  <run><digits>.mrf           optional. Time_hr, RainMin_mm_hr, MAP_mm_hr (labeled by hour END), Srf_m3_s, ChannelPerc_m3,
                              MDGW_mm, AvSWE_cm. Gives the rain for runoff ratios and the water-table / channel-loss checks.
  <run><digits>.rft           optional. Looked for columns whose names contain Hsrf and Sbsrf (see PROCESS CHECKS below).
  ../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv   (default --obs): datetime, q_cms, ... (Arizona time)
  the run's start date and time: from --window (W1..W5, V) or --start. A run id that begins with W1_ ... W5_ or V_
  (the design's run-id pattern) is read for its window automatically, so --scan needs no --window.

RUN WINDOWS (from the design, Stages 0 to 8)
--------------------------------------------
  W1 Jan 18 - Jan 25 2017 (168 h)     W2 Jan 18 - Mar 12 2017 (about 1,270 h)    W3 Oct 1 2016 - Jan 25 2017 (2,784 h)
  W4 Oct 1 2016 - Apr 1 2017 (4,368 h)    W5 Oct 1 2016 - Oct 1 2017 (8,760 h)   V  Sep 1 2019 - Jun 1 2020 (6,576 h)
A metric is computed only if the run covers the whole metric window; otherwise it is left blank (NaN). So a W1 run
gets the flood-week metrics and nothing else; a W5 run gets everything.

THE SCORING PROFILE
-------------------
  wy2017 (default): the windows below, from the design doc and the earlier cells. V runs use --profile wy2020, which has
  only the scored period, the months and the whole-period metrics; add WY2020 storm windows later with --profile-file,
  after looking at the WY2020 GAUGE record (not model results).
  Every window can be changed in the JSON written by `dump-profile`.

METRICS (design: "Which metric, when, and where")
-------------------------------------------------
  volume ratio        model / gauge volume; on every event, every month, the whole scored span. The model volume is counted
                      only over times when the gauge has data (so a gauge blank can never make the model look too high), and
                      Jan 1-2 2017 is NOT excluded from either. The raw model volume is also saved.
  peak error (%)      model peak / gauge peak - 1 inside the window (not an arrival time)
  first-arrival (h)   time the flow first RISES THROUGH 20% of that window's gauge peak: model minus gauge. Blank if the flow never
                      gets there, or starts above that level and never dips below it (the flag column says which; a second storm
                      that arrives while the first storm's recession is still above the level can have no arrival time)
  tail volume ratio   model / gauge volume on each recession tail; and recession rate ratio = the model's decay rate divided by
                      the gauge's. A decay rate is the slope of ln(flow) against time over the hourly flows above 0.01 m3/s
                      (my definition: the design names the metric, not the formula). 1 = same recession; above 1 = model dries faster
  false flow          Oct-Dec (and Jul-Sep) model volume; the number of days whose highest flow is above the threshold
                      (0.13 m3/s, the gauge's Oct-Dec maximum; the gauge's real maximum is saved next to it)
  storm balance       runoff ratio (outflow volume / basin rain volume) per storm from the .mrf, and the later storm over the
                      earlier one (design: gauge 3.4 for Feb 27-28 over Feb 18-19; no pass or fail yet)
  KGE_2012 + r, beta, gamma, NSE, PBIAS   per event and over the span. The model is interpolated onto the gauge's 15-minute
                      times first. KGE is read only after the volume gate and never alone: the report prints it, with the gate state
  process checks      water table (.mrf MDGW), channel loss and snow (.mrf), Hsrf and Sbsrf shares (.rft). They are not scores.

TOLERANCES in the report are the design's PROPOSALS (volume within 10%, peak within 15%, arrival within 1 h, tail ratio at least 0.7,
Oct-Dec false flow at least halved from 2.50 M m3). Change them in TOLERANCES at the top of this file once they are agreed.

KNOWN LIMITS
------------
  * Never run on the real outlet files or the real gauge file (written without access to them). `self-test` checks every metric on
    synthetic runs with known answers. Run `score --dry-run` and then the baseline first (see the REFERENCE CHECKS block of the report).
  * The tail windows (Jan 26 - Feb 2 and Mar 3 - Mar 11, each INCLUDING the last named day) are my reading of the design. The report
    prints the gauge volumes for them next to the numbers in the design doc (1.17 and 1.47 M m3); if they do not match, the windows differ.
  * The .rft column names (Hsrf, Sbsrf) are assumed. If they are not found the report prints the column names it saw.
  * Fingerprint rules (Series 122) come from synthetic SMF truths with no gauge or rain error; they may not hold at Sycamore.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
import shutil
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

VERSION = "1.0"
HOUR = pd.Timedelta(hours=1)
BASIN_KM2 = 424.6                       # watershed area from the build (424.6151 km2)
COVER_TOL_H = 1.5                       # a run "covers" a window if it reaches within this many hours of both ends
MAX_STEP_H = 3.0                        # a gap between two flow samples longer than this counts as a blank (no volume)
RECESSION_FLOOR_M3S = 0.01              # recession rate uses hourly flows at or above this
DEFAULT_OBS = "../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv"

# ----------------------------------------------------------------------------------------------------------------
# RUN WINDOWS (where a run starts and stops). From the design doc: W1 168 h, W2 about 1,270 h, W3 2,784 h, W4 4,368 h,
# W5 8,760 h, V 6,576 h. Start/stop are the local clock times those hours work out to.
# ----------------------------------------------------------------------------------------------------------------
RUN_WINDOWS = {
    "W1": ("2017-01-18 00:00", "2017-01-25 00:00"),
    "W2": ("2017-01-18 00:00", "2017-03-12 00:00"),
    "W3": ("2016-10-01 00:00", "2017-01-25 00:00"),
    "W4": ("2016-10-01 00:00", "2017-04-01 00:00"),
    "W5": ("2016-10-01 00:00", "2017-10-01 00:00"),
    "V":  ("2019-09-01 00:00", "2020-06-01 00:00"),
}

# ----------------------------------------------------------------------------------------------------------------
# SCORING PROFILES: every window the metrics use. All times are local Arizona time, [start, end).
#   events        volume, peak, first arrival, KGE (the storm windows WITHOUT recession tails)
#   volume_windows  volume and runoff ratio (the storm windows WITH their tails, as in storm_rain_cell.py)
#   tails         tail volume ratio and recession rate ratio
#   false_flow    dry-period false flow
#   balance_pairs [later, earlier] pairs for the storm-to-storm runoff-ratio contrast
#   named_volume  other windows whose volume is wanted (the Oct-Mar total used in the baseline reproduction check)
# ----------------------------------------------------------------------------------------------------------------
PROFILES = {
    "wy2017": {
        "period": ["2016-10-01 00:00", "2017-10-01 00:00"],
        "exclude": [],            # none: the Jan 1 gap is blank in the gauge file, so matching already skips it; 
                                  # Jan 2-3 estimated values count
        "events": {
            "flood_week":       ["2017-01-18 00:00", "2017-01-25 00:00"],
            "jan_second_storm": ["2017-01-23 12:00", "2017-01-25 00:00"],
            "feb18_19":         ["2017-02-18 00:00", "2017-02-21 00:00"],
            "feb27_28":         ["2017-02-27 00:00", "2017-03-03 00:00"],
        },
        "volume_windows": {
            "jan_flood_tail":   ["2017-01-18 00:00", "2017-02-03 00:00"],
            "feb18_19_tail":    ["2017-02-18 00:00", "2017-02-27 00:00"],
            "feb27_28_tail":    ["2017-02-27 00:00", "2017-03-12 00:00"],
        },
        "tails": {
            "tail_jan26_feb2":  ["2017-01-26 00:00", "2017-02-03 00:00"],
            "tail_mar3_mar11":  ["2017-03-03 00:00", "2017-03-12 00:00"],
        },
        "false_flow": {
            "octdec": ["2016-10-01 00:00", "2017-01-01 00:00"],
            "julsep": ["2017-07-01 00:00", "2017-10-01 00:00"],
        },
        "false_flow_threshold_m3s": 0.13,
        "balance_pairs": [["feb27_28_tail", "feb18_19_tail"], ["feb27_28", "feb18_19"]],
        "named_volume": {"oct_mar": ["2016-10-01 00:00", "2017-04-01 00:00"]},
        "months": True,
    },
    "wy2020": {
        "period": ["2019-10-01 00:00", "2020-06-01 00:00"],               # scored Oct 1 2019 - May 31 2020 (the Bush Fire was June 2020)
        "exclude": [],
        "events": {}, "volume_windows": {}, "tails": {}, "false_flow": {}, "balance_pairs": [], "named_volume": {},
        "false_flow_threshold_m3s": None,
        "months": True,
    },
}

# Gauge numbers and baseline numbers quoted in the design doc, used ONLY for the reference-check block of the report.
DESIGN_REFERENCE = {
    "oct_mar_gauge_total_Mm3": 17.87,
    "oct_mar_baseline_model_total_Mm3": 15.30,
    "tail_jan26_feb2_gauge_Mm3": 1.17, "tail_jan26_feb2_baseline_model_Mm3": 0.26,
    "tail_mar3_mar11_gauge_Mm3": 1.47, "tail_mar3_mar11_baseline_model_Mm3": 0.09,
    "octdec_baseline_model_false_flow_Mm3": 2.50,
    "feb27_28_over_feb18_19_runoff_ratio_gauge": 3.4,
}

# ----------------------------------------------------------------------------------------------------------------
# TOLERANCES: the design's PROPOSALS, to be fixed before the first run. (metric group, limit)
# ----------------------------------------------------------------------------------------------------------------
TOLERANCES = {
    "volume_ratio_within": 0.10,            # |model/gauge - 1| <= 0.10 per big event and for the span
    "peak_error_within_pct": 15.0,
    "arrival_error_within_h": 1.0,
    "tail_ratio_at_least": 0.70,
    "false_flow_must_fall_to_fraction": 0.5,   # at least halve ...
    "false_flow_baseline_Mm3": 2.50,           # ... this baseline volume
}


# ================================================================================================================
# small helpers
# ================================================================================================================
def die(msg: str, code: int = 2) -> None:
    print(f"\nERROR: {msg}\n", file=sys.stderr)
    sys.exit(code)


def ts(x) -> pd.Timestamp:
    return pd.Timestamp(x)


def nz(x):
    return np.nan if x is None else x


def run_window_from_id(run_id: str):
    """W1_..., W3-..., V_... at the start of a run id -> 'W1'.. 'V', else None."""
    m = re.match(r"^(W[1-5]|V)(?=[_\-\s.]|$)", run_id, flags=re.IGNORECASE)
    return m.group(1).upper() if m else None


def load_profile(name: str, profile_file: str | None):
    if profile_file:
        prof = json.loads(Path(profile_file).read_text())
        base = copy.deepcopy(PROFILES["wy2017"])
        base.update(prof)
        return base
    if name not in PROFILES:
        die(f"unknown profile {name!r}; choose from {sorted(PROFILES)} or give --profile-file")
    return copy.deepcopy(PROFILES[name])


# ================================================================================================================
# reading the files
# ================================================================================================================
def read_qout(path: Path, start) -> pd.Series:
    """Routed outlet flow (m3/s) as a time series. tRIBS 6.0.0 writes a comma-delimited file with a header line;
    an older whitespace-delimited layout is also accepted."""
    with open(path) as fh:
        fh.readline()
        second = fh.readline()
    if "," in second:
        raw = pd.read_csv(path, sep=",", header=0)
    else:
        raw = pd.read_csv(path, sep=r"\s+", skiprows=1, header=None)
    raw = raw.iloc[:, :2].copy()
    raw.columns = ["Time_hr", "Q"]
    raw = raw.apply(pd.to_numeric, errors="coerce").dropna()
    s = pd.Series(raw["Q"].to_numpy(float), index=ts(start) + pd.to_timedelta(raw["Time_hr"].to_numpy(), unit="h"))
    s = s[~s.index.duplicated(keep="last")].sort_index()
    return s


def read_obs(path: str) -> pd.Series:
    if not Path(path).exists():
        die(f"gauge file not found: {path} (use --obs)")
    obs = pd.read_csv(path, parse_dates=["datetime"]).set_index("datetime")
    s = pd.to_numeric(obs["q_cms"], errors="coerce")
    s = s[~s.index.duplicated(keep="last")].sort_index()
    return s


def read_mrf(path: Path, start) -> pd.DataFrame:
    sim = pd.read_csv(path)
    sim = sim[sim["RainMin_mm_hr"] < 9000]                       # drop the unfilled last row
    sim.index = ts(start) + pd.to_timedelta(sim["Time_hr"], unit="h")
    return sim


def read_rft(path: Path, start):
    """Hsrf / Sbsrf columns from a .rft file, or (None, column names seen). Column names are ASSUMED (not seen in the project)."""
    try:
        df = pd.read_csv(path, sep=",")
        if df.shape[1] < 2:
            df = pd.read_csv(path, sep=r"\s+")
    except Exception as exc:  # noqa: BLE001
        return None, [f"unreadable: {exc}"]
    cols = list(df.columns)
    low = {c.lower().strip(): c for c in cols}
    hs = next((low[k] for k in low if "hsrf" in k), None)
    sb = next((low[k] for k in low if "sbsrf" in k), None)
    tcol = next((low[k] for k in low if k.startswith("time")), cols[0])
    if hs is None or sb is None:
        return None, cols
    t = pd.to_numeric(df[tcol], errors="coerce")
    out = pd.DataFrame({"Hsrf": pd.to_numeric(df[hs], errors="coerce"), "Sbsrf": pd.to_numeric(df[sb], errors="coerce")})
    out.index = ts(start) + pd.to_timedelta(t, unit="h")
    return out[out.index.notna()], cols


def find_run_files(run_id: str, results_dir: Path):
    qout = sorted(results_dir.glob(f"**/{run_id}_Outlet.qout"))
    if len(qout) != 1:
        return None, None, None, f"expected exactly one {run_id}_Outlet.qout under {results_dir}, found {len(qout)}"
    mrf = sorted(results_dir.glob(f"**/{run_id}[0-9]*.mrf"))
    rft = sorted(results_dir.glob(f"**/{run_id}[0-9]*.rft"))
    if not rft:
        alt = sorted(results_dir.glob(f"**/{run_id}*.rft"))
        rft = alt if len(alt) == 1 else []
    return qout[0], (mrf[-1] if mrf else None), (rft[-1] if rft else None), ""


def checksum_state(qout: Path, write: bool):
    side = qout.with_suffix(qout.suffix + ".sha256")
    h = hashlib.sha256(qout.read_bytes()).hexdigest()
    if write:
        side.write_text(f"{h}  {qout.name}\n")
        return "written", h
    if side.exists():
        saved = side.read_text().split()[0]
        return ("ok" if saved == h else "MISMATCH"), h
    return "no sidecar", h


# ================================================================================================================
# flow series with volumes (the trapezoid rule, as in compare_outlet_cell.py, plus blank handling)
# ================================================================================================================
class Flow:
    """A flow series with its trapezoid steps. A step is blank (no volume) if either end is NaN or the gap is over MAX_STEP_H."""

    def __init__(self, q: pd.Series):
        q = q[~q.index.duplicated(keep="last")].sort_index()
        self.q = q
        t = q.index.values.astype("datetime64[s]").astype("int64")
        v = q.to_numpy(float)
        self.t0, self.t1 = t[:-1], t[1:]
        dt = np.diff(t)
        self.valid = ~np.isnan(v[1:]) & ~np.isnan(v[:-1]) & (dt <= MAX_STEP_H * 3600)
        self.vol = np.where(self.valid, 0.5 * (v[1:] + v[:-1]) * dt, 0.0)
        self.n_blank = int((~self.valid).sum())

    def _keep(self, a, b, exclude):
        a_s, b_s = ts(a).value // 10**9, ts(b).value // 10**9
        keep = (self.t0 >= a_s) & (self.t0 < b_s)
        for ea, eb in exclude:
            keep &= ~((self.t0 >= ts(ea).value // 10**9) & (self.t0 < ts(eb).value // 10**9))
        return keep

    def volume(self, a, b, exclude=(), match_to: "Flow | None" = None) -> float:
        """Volume (m3) of the steps that START in [a, b), minus the excluded times. With match_to, only steps whose midpoint
        falls in a valid (non-blank) step of that other series count."""
        keep = self._keep(a, b, exclude) & self.valid
        if match_to is not None and len(match_to.t0):
            mid = 0.5 * (self.t0 + self.t1)
            i = np.searchsorted(match_to.t0, mid, side="right") - 1
            ok = (i >= 0)
            ii = np.clip(i, 0, len(match_to.t0) - 1)
            ok &= (mid < match_to.t1[ii]) & match_to.valid[ii]
            keep &= ok
        return float(self.vol[keep].sum())

    def has_data(self, a, b, exclude=()) -> bool:
        return bool((self._keep(a, b, exclude) & self.valid).any())

    def window(self, a, b, exclude=()) -> pd.Series:
        s = self.q[(self.q.index >= ts(a)) & (self.q.index < ts(b))]
        for ea, eb in exclude:
            s = s[~((s.index >= ts(ea)) & (s.index < ts(eb)))]
        return s


def covers(flow: Flow, a, b) -> bool:
    return (flow.q.index[0] <= ts(a) + pd.Timedelta(hours=COVER_TOL_H)) and (flow.q.index[-1] >= ts(b) - pd.Timedelta(hours=COVER_TOL_H))


# ================================================================================================================
# metric functions
# ================================================================================================================
def kge_parts(model_q: pd.Series, obs_q: pd.Series):
    """KGE (2012 version) with r, beta, gamma, plus NSE and PBIAS (%). The model is interpolated onto the gauge's times first
    (the model writes its own, finer, time steps). Same method as scores_cell.py."""
    out = dict(kge=np.nan, r=np.nan, beta=np.nan, gamma=np.nan, nse=np.nan, pbias=np.nan, n=0)
    o = obs_q.dropna()
    if len(o) < 3 or model_q.dropna().empty:
        return out
    both = model_q.index.union(o.index)
    m = model_q.reindex(both).interpolate(method="time", limit_area="inside").reindex(o.index)
    ok = m.notna() & o.notna()
    m, o = m[ok].to_numpy(float), o[ok].to_numpy(float)
    out["n"] = int(len(o))
    if len(o) < 3 or o.mean() == 0 or o.std() == 0:
        return out
    denom = ((o - o.mean()) ** 2).sum()
    out["nse"] = float(1 - ((m - o) ** 2).sum() / denom)
    out["pbias"] = float(100 * (m.sum() - o.sum()) / o.sum())
    out["beta"] = float(m.mean() / o.mean())
    if m.std() == 0 or m.mean() == 0:                      # no flow, or constant flow: r and gamma are not defined
        return out
    out["r"] = float(np.corrcoef(m, o)[0, 1])
    out["gamma"] = float((m.std() / m.mean()) / (o.std() / o.mean()))
    out["kge"] = float(1 - np.sqrt((out["r"] - 1) ** 2 + (out["beta"] - 1) ** 2 + (out["gamma"] - 1) ** 2))
    return out


def arrival(q: pd.Series, level: float):
    """First time the flow RISES THROUGH `level` (linear interpolation between samples). Returns (time, flag)."""
    s = q.dropna()
    if s.empty:
        return pd.NaT, "no data"
    above = s.to_numpy() >= level
    rising = np.nonzero(above[1:] & ~above[:-1])[0] + 1          # first sample at or above the level, after one below it
    if len(rising) == 0:
        return pd.NaT, ("already above at start" if above[0] else "never reaches")
    i = int(rising[0])
    t_prev, t_next = s.index[i - 1], s.index[i]
    q_prev, q_next = float(s.iloc[i - 1]), float(s.iloc[i])
    if (t_next - t_prev) > pd.Timedelta(hours=MAX_STEP_H) or q_next == q_prev:
        return t_next, "ok"
    return t_prev + (t_next - t_prev) * ((level - q_prev) / (q_next - q_prev)), "ok"


def recession_rate(q: pd.Series, a, b, floor=RECESSION_FLOOR_M3S):
    """Decay rate k (per hour) = minus the slope of ln(flow) against time over the hourly mean flows >= floor. Returns (k, n)."""
    s = q[(q.index >= ts(a)) & (q.index < ts(b))].dropna()
    if len(s) < 3:
        return np.nan, 0
    hourly = s.resample("60min").mean().dropna()
    hourly = hourly[hourly >= floor]
    if len(hourly) < 12:
        return np.nan, int(len(hourly))
    t = np.asarray((hourly.index - hourly.index[0]) / HOUR, dtype=float)
    slope = np.polyfit(t, np.log(hourly.to_numpy(float)), 1)[0]
    return float(-slope), int(len(hourly))


def days_above(q: pd.Series, a, b, thr, exclude=()):
    s = q[(q.index >= ts(a)) & (q.index < ts(b))]
    for ea, eb in exclude:
        s = s[~((s.index >= ts(ea)) & (s.index < ts(eb)))]
    s = s.dropna()
    if s.empty:
        return np.nan
    return int((s.resample("D").max().dropna() > thr).sum())


# ================================================================================================================
# scoring one run
# ================================================================================================================
def month_windows(period):
    start, end = ts(period[0]), ts(period[1])
    cur = start.replace(day=1)
    while cur < end:
        nxt = (cur + pd.offsets.MonthBegin(1))
        if cur >= start and nxt <= end + pd.Timedelta(hours=COVER_TOL_H):
            yield f"month_{cur:%Y_%m}", cur, nxt
        cur = nxt


def score_run(run_id, model_q, gauge_q, mrf, rft_df, rft_cols, profile, basin_km2):
    """Returns (row dict, list of (group, label, value, ok/None, note)) for the tolerance table, and the report lines."""
    M, G = Flow(model_q), Flow(gauge_q)
    excl = [tuple(e) for e in profile["exclude"]]
    row, checks, rep = {}, [], []
    rain_hr = None
    if mrf is not None:
        rain_hr = mrf["MAP_mm_hr"].copy()
        rain_hr.index = rain_hr.index - HOUR                 # .mrf labels each hour by its END; hours are labeled by their START here

    def ratio(a, b):
        return a / b if (b is not None and np.isfinite(b) and b > 0) else np.nan

    def vol_pair(a, b):
        if not covers(M, a, b):
            return np.nan, np.nan, np.nan
        g = G.volume(a, b, excl) / 1e6 if G.has_data(a, b, excl) else np.nan
        mm = M.volume(a, b, excl, match_to=G) / 1e6       # matched: only where the gauge has data, Jan 1-2 excluded
        mr = M.volume(a, b) / 1e6                         # raw: every hour of the model, nothing excluded
        return mm, mr, g

    # ---- run description ----
    row["model_first"], row["model_last"] = model_q.index[0], model_q.index[-1]
    row["model_blank_steps"], row["gauge_blank_steps"] = M.n_blank, G.n_blank

    # ---- the scored span ----
    p0, p1 = ts(profile["period"][0]), ts(profile["period"][1])
    s0 = max(p0, model_q.index[0])
    s1 = p1 if model_q.index[-1] >= p1 - pd.Timedelta(hours=COVER_TOL_H) else model_q.index[-1]
    row["span__start"], row["span__end"] = s0, s1
    row["span__is_full_period"] = bool(s0 == p0 and s1 == p1)
    mm = M.volume(s0, s1, excl, match_to=G) / 1e6
    mr = M.volume(s0, s1) / 1e6
    g = G.volume(s0, s1, excl) / 1e6 if G.has_data(s0, s1, excl) else np.nan
    row.update({"span__vol_model_Mm3_matched": mm, "span__vol_model_Mm3_raw": mr, "span__vol_gauge_Mm3": g, "span__vol_ratio": ratio(mm, g)})
    kp = kge_parts(M.window(s0, s1, excl), G.window(s0, s1, excl))
    for k, v in kp.items():
        row[f"span__{ 'kge2012' if k == 'kge' else ('kge_' + k if k in ('r', 'beta', 'gamma') else k) }"] = v
    span_full = row["span__is_full_period"]
    checks.append(("volume", f"span ({s0:%b %d %Y} to {s1:%b %d %Y})" + ("" if span_full else ", partial"), row["span__vol_ratio"],
                   None if not np.isfinite(row["span__vol_ratio"]) else abs(row["span__vol_ratio"] - 1) <= TOLERANCES["volume_ratio_within"],
                   f"within {TOLERANCES['volume_ratio_within']:.0%} of 1"))

    # ---- months ----
    if profile.get("months"):
        for name, a, b in month_windows(profile["period"]):
            mm, mr, g = vol_pair(a, b)
            row[f"{name}__vol_model_Mm3_matched"], row[f"{name}__vol_model_Mm3_raw"] = mm, mr
            row[f"{name}__vol_gauge_Mm3"], row[f"{name}__vol_ratio"] = g, ratio(mm, g)
            if rain_hr is not None and covers(M, a, b):
                row[f"{name}__rain_mm"] = float(rain_hr[(rain_hr.index >= a) & (rain_hr.index < b)].sum())

    # ---- named volume windows (Oct-Mar total etc.) ----
    for name, (a, b) in profile.get("named_volume", {}).items():
        mm, mr, g = vol_pair(a, b)
        row[f"{name}__vol_model_Mm3_matched"], row[f"{name}__vol_model_Mm3_raw"] = mm, mr
        row[f"{name}__vol_gauge_Mm3"], row[f"{name}__vol_ratio"] = g, ratio(mm, g)
        if np.isfinite(row[f"{name}__vol_ratio"]):
            checks.append(("volume", f"{name} ({ts(a):%b %d} to {ts(b):%b %d})", row[f"{name}__vol_ratio"],
                           abs(row[f"{name}__vol_ratio"] - 1) <= TOLERANCES["volume_ratio_within"], f"within {TOLERANCES['volume_ratio_within']:.0%} of 1"))

    # ---- volume windows with tails + runoff ratio ----
    vol_wins = dict(profile.get("volume_windows", {}))
    for name, (a, b) in vol_wins.items():
        mm, mr, g = vol_pair(a, b)
        row[f"{name}__vol_model_Mm3_matched"], row[f"{name}__vol_model_Mm3_raw"] = mm, mr
        row[f"{name}__vol_gauge_Mm3"], row[f"{name}__vol_ratio"] = g, ratio(mm, g)
        if name.startswith("feb") and np.isfinite(row[f"{name}__vol_ratio"]):
            checks.append(("volume", f"{name} (each February storm, with tail)", row[f"{name}__vol_ratio"],
                           abs(row[f"{name}__vol_ratio"] - 1) <= TOLERANCES["volume_ratio_within"], f"within {TOLERANCES['volume_ratio_within']:.0%} of 1"))

    # ---- events: volume, peak, arrival, KGE ----
    gate_flood = None
    for name, (a, b) in profile.get("events", {}).items():
        mm, mr, g = vol_pair(a, b)
        row[f"{name}__vol_model_Mm3_matched"], row[f"{name}__vol_model_Mm3_raw"] = mm, mr
        row[f"{name}__vol_gauge_Mm3"], row[f"{name}__vol_ratio"] = g, ratio(mm, g)
        ev_m, ev_g = M.window(a, b, excl), G.window(a, b, excl)
        cov = covers(M, a, b)
        pm = float(ev_m.max()) if (cov and ev_m.notna().any()) else np.nan
        pg = float(ev_g.max()) if ev_g.notna().any() else np.nan
        row[f"{name}__peak_model_m3s"], row[f"{name}__peak_gauge_m3s"] = pm, pg
        row[f"{name}__peak_err_pct"] = (pm / pg - 1) * 100 if (np.isfinite(pm) and np.isfinite(pg) and pg > 0) else np.nan
        row[f"{name}__peak_time_diff_h"] = ((ev_m.idxmax() - ev_g.idxmax()) / HOUR) if (np.isfinite(pm) and np.isfinite(pg)) else np.nan
        flag, err = "not covered", np.nan
        if cov and np.isfinite(pg) and pg > 0:
            level = 0.2 * pg
            tg, fg = arrival(ev_g, level)
            tm, fm = arrival(ev_m, level)
            row[f"{name}__arrival_level_m3s"] = level
            if fg != "ok":
                flag = f"gauge: {fg}"
            elif fm != "ok":
                flag = f"model: {fm}"
            else:
                err, flag = (tm - tg) / HOUR, "ok"
        row[f"{name}__arrival_err_h"], row[f"{name}__arrival_flag"] = err, flag
        kp = kge_parts(ev_m, ev_g) if cov else dict(kge=np.nan, r=np.nan, beta=np.nan, gamma=np.nan, nse=np.nan, pbias=np.nan, n=0)
        row[f"{name}__kge2012"], row[f"{name}__kge_r"], row[f"{name}__kge_beta"], row[f"{name}__kge_gamma"] = kp["kge"], kp["r"], kp["beta"], kp["gamma"]
        row[f"{name}__nse"], row[f"{name}__pbias_pct"] = kp["nse"], kp["pbias"]
        v = row[f"{name}__vol_ratio"]
        if np.isfinite(v):
            checks.append(("volume", f"{name}", v, abs(v - 1) <= TOLERANCES["volume_ratio_within"], f"within {TOLERANCES['volume_ratio_within']:.0%} of 1"))
        if name == "flood_week":
            gate_flood = (abs(v - 1) <= TOLERANCES["volume_ratio_within"]) if np.isfinite(v) else None
        pe = row[f"{name}__peak_err_pct"]
        if np.isfinite(pe):
            checks.append(("peak", f"{name}", pe, abs(pe) <= TOLERANCES["peak_error_within_pct"], f"within +/-{TOLERANCES['peak_error_within_pct']:.0f}%"))
        if np.isfinite(err):
            checks.append(("arrival", f"{name}", err, abs(err) <= TOLERANCES["arrival_error_within_h"], f"within +/-{TOLERANCES['arrival_error_within_h']:.0f} h"))
        elif flag not in ("not covered",):
            checks.append(("arrival", f"{name}", np.nan, None, flag))

    # ---- tails ----
    for name, (a, b) in profile.get("tails", {}).items():
        mm, mr, g = vol_pair(a, b)
        row[f"{name}__vol_model_Mm3"], row[f"{name}__vol_gauge_Mm3"] = mm, g
        tr = ratio(mm, g)
        row[f"{name}__tail_vol_ratio"] = tr
        if covers(M, a, b):
            km, nm = recession_rate(model_q, a, b)
            kg, ng = recession_rate(gauge_q, a, b)
        else:
            km, nm, kg, ng = np.nan, 0, np.nan, 0
        row[f"{name}__k_model_per_h"], row[f"{name}__k_gauge_per_h"] = km, kg
        row[f"{name}__n_hours_model"], row[f"{name}__n_hours_gauge"] = nm, ng
        row[f"{name}__recession_rate_ratio"] = km / kg if (np.isfinite(km) and np.isfinite(kg) and kg > 0) else np.nan
        if np.isfinite(tr):
            checks.append(("tail", f"{name}", tr, tr >= TOLERANCES["tail_ratio_at_least"], f"at least {TOLERANCES['tail_ratio_at_least']:.1f}"))

    # ---- false flow ----
    thr = profile.get("false_flow_threshold_m3s")
    for name, (a, b) in profile.get("false_flow", {}).items():
        if not covers(M, a, b):
            continue
        mr = M.volume(a, b, excl) / 1e6
        g = G.volume(a, b, excl) / 1e6 if G.has_data(a, b, excl) else np.nan
        row[f"{name}__vol_model_Mm3"], row[f"{name}__vol_gauge_Mm3"] = mr, g
        wm, wg = M.window(a, b, excl), G.window(a, b, excl)
        row[f"{name}__max_model_m3s"] = float(wm.max()) if wm.notna().any() else np.nan
        row[f"{name}__max_gauge_m3s"] = float(wg.max()) if wg.notna().any() else np.nan
        row[f"{name}__threshold_m3s"] = nz(thr)
        row[f"{name}__days_above_model"] = days_above(model_q, a, b, thr, excl) if thr is not None else np.nan
        row[f"{name}__days_above_gauge"] = days_above(gauge_q, a, b, thr, excl) if thr is not None else np.nan
        if name == "octdec":
            lim = TOLERANCES["false_flow_must_fall_to_fraction"] * TOLERANCES["false_flow_baseline_Mm3"]
            checks.append(("false flow", f"{name} model volume (M m3)", mr, mr <= lim, f"at most {lim:.2f} (half of {TOLERANCES['false_flow_baseline_Mm3']:.2f})"))

    # ---- storm balance (runoff ratios) ----
    rr_model, rr_gauge = {}, {}
    if rain_hr is not None:
        windows_for_rr = {**profile.get("events", {}), **vol_wins}
        for name, (a, b) in windows_for_rr.items():
            if not covers(M, a, b):
                continue
            rain_mm = float(rain_hr[(rain_hr.index >= ts(a)) & (rain_hr.index < ts(b))].sum())
            rain_vol = rain_mm * basin_km2 / 1000.0                    # mm over km2 -> M m3
            row[f"{name}__rain_mm"] = rain_mm
            vm, vg = row.get(f"{name}__vol_model_Mm3_matched", np.nan), row.get(f"{name}__vol_gauge_Mm3", np.nan)
            row[f"{name}__runoff_ratio_model"] = vm / rain_vol if rain_vol > 0 else np.nan
            row[f"{name}__runoff_ratio_gauge"] = vg / rain_vol if rain_vol > 0 else np.nan
            rr_model[name], rr_gauge[name] = row[f"{name}__runoff_ratio_model"], row[f"{name}__runoff_ratio_gauge"]
    for late, early in profile.get("balance_pairs", []):
        key = f"balance__{late}_over_{early}"
        for lab, rr in (("model", rr_model), ("gauge", rr_gauge)):
            a_, b_ = rr.get(late, np.nan), rr.get(early, np.nan)
            row[f"{key}__{lab}"] = a_ / b_ if (np.isfinite(a_) and np.isfinite(b_) and b_ > 0) else np.nan

    # ---- process checks ----
    if mrf is not None:
        wt = mrf["MDGW_mm"]
        row["proc__mdgw_first_mm"] = float(wt.iloc[0])
        row["proc__mdgw_min_mm"], row["proc__mdgw_max_mm"] = float(wt.min()), float(wt.max())
        fw = profile.get("events", {}).get("flood_week")
        if fw and wt.index[0] <= ts(fw[0]) <= wt.index[-1]:
            row["proc__mdgw_at_flood_week_start_mm"] = float(wt.asof(ts(fw[0])))
        row["proc__channel_loss_Mm3"] = float(mrf["ChannelPerc_m3"].sum() / 1e6)
        row["proc__snow_max_cm"] = float(mrf["AvSWE_cm"].max())
    if rft_df is not None:
        a, b = s0, s1
        part = rft_df[(rft_df.index >= a) & (rft_df.index < b)]
        if len(part):
            hs, sb = float(part["Hsrf"].sum()), float(part["Sbsrf"].sum())
            row["proc__hsrf_mean"], row["proc__sbsrf_mean"] = float(part["Hsrf"].mean()), float(part["Sbsrf"].mean())
            row["proc__hsrf_share"] = hs / (hs + sb) if (hs + sb) > 0 else np.nan
        fw = profile.get("events", {}).get("flood_week")
        if fw:
            pf = rft_df[(rft_df.index >= ts(fw[0])) & (rft_df.index < ts(fw[1]))]
            if len(pf):
                hs, sb = float(pf["Hsrf"].sum()), float(pf["Sbsrf"].sum())
                row["proc__hsrf_share_flood_week"] = hs / (hs + sb) if (hs + sb) > 0 else np.nan

    # ---- KGE gate state ----
    span_ok = row["span__vol_ratio"]
    gate = []
    if gate_flood is not None:
        gate.append(gate_flood)
    if np.isfinite(row.get("oct_mar__vol_ratio", np.nan)):
        gate.append(abs(row["oct_mar__vol_ratio"] - 1) <= TOLERANCES["volume_ratio_within"])
    row["kge_volume_gate"] = ("open" if all(gate) else "CLOSED") if gate else "not evaluated"

    # ---- reference checks vs the design doc ----   (label, my value, design value, kind: "gauge" must match, "baseline" only for the Oct 21 baseline)
    ref = []
    def add_ref(label, value, design, kind, extra=None):
        if np.isfinite(value) or (extra is not None and np.isfinite(extra)):
            ref.append((label, value, design, kind, extra))
    add_ref("gauge Oct-Mar total (M m3)", row.get("oct_mar__vol_gauge_Mm3", np.nan), DESIGN_REFERENCE["oct_mar_gauge_total_Mm3"], "gauge")
    add_ref("model Oct-Mar total, every hour (M m3)", row.get("oct_mar__vol_model_Mm3_raw", np.nan), DESIGN_REFERENCE["oct_mar_baseline_model_total_Mm3"], "baseline")
    for tname in ("tail_jan26_feb2", "tail_mar3_mar11"):
        add_ref(f"gauge volume, {tname} (M m3)", row.get(f"{tname}__vol_gauge_Mm3", np.nan), DESIGN_REFERENCE[f"{tname}_gauge_Mm3"], "gauge")
        add_ref(f"model volume, {tname} (M m3)", row.get(f"{tname}__vol_model_Mm3", np.nan), DESIGN_REFERENCE[f"{tname}_baseline_model_Mm3"], "baseline")
    add_ref("model Oct-Dec false-flow volume (M m3)", row.get("octdec__vol_model_Mm3", np.nan), DESIGN_REFERENCE["octdec_baseline_model_false_flow_Mm3"], "baseline")
    add_ref("gauge runoff ratio, Feb 27-28 over Feb 18-19: long windows (first-days windows in the next number)",
            row.get("balance__feb27_28_tail_over_feb18_19_tail__gauge", np.nan), DESIGN_REFERENCE["feb27_28_over_feb18_19_runoff_ratio_gauge"], "gauge",
            extra=row.get("balance__feb27_28_over_feb18_19__gauge", np.nan))
    row["_reference"] = ref
    return row, checks


# ================================================================================================================
# the text report for one run
# ================================================================================================================
def fmt(x, nd=2):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "-"
    if isinstance(x, (pd.Timestamp,)):
        return f"{x:%b %d %H:%M}"
    return f"{x:.{nd}f}" if isinstance(x, (float, np.floating)) else str(x)


def report_lines(run_id, window, row, checks, profile, checksum_note, notes):
    L = []
    L.append(f"RUN {run_id}   (window {window or '?'})   scorer v{VERSION}")
    L.append(f"  model flow {row['model_first']:%b %d %Y %H:%M} to {row['model_last']:%b %d %Y %H:%M}; scored span {row['span__start']:%b %d %Y} to "
             f"{row['span__end']:%b %d %Y}" + ("" if row["span__is_full_period"] else "  (PARTIAL: not the whole scoring period)"))
    L.append(f"  blank steps skipped: model {row['model_blank_steps']}, gauge {row['gauge_blank_steps']}; excluded: "
             + (", ".join(f"{a} to {b}" for a, b in profile["exclude"]) or "none") + f"; outlet file checksum: {checksum_note}")
    for n in notes:
        L.append("  note: " + n)
    L.append("")
    # volume table
    names = (list(profile.get("events", {})) + list(profile.get("volume_windows", {})) + list(profile.get("named_volume", {})) + ["span"])
    L.append("VOLUME  (M m3; model counted only where the gauge has data; ratio = model / gauge)")
    L.append(f"  {'window':<22}{'model':>9}{'gauge':>9}{'ratio':>8}{'model raw':>11}")
    for n in names:
        if f"{n}__vol_ratio" not in row:
            continue
        L.append(f"  {n:<22}{fmt(row[f'{n}__vol_model_Mm3_matched']):>9}{fmt(row[f'{n}__vol_gauge_Mm3']):>9}{fmt(row[f'{n}__vol_ratio'], 3):>8}"
                 f"{fmt(row[f'{n}__vol_model_Mm3_raw']):>11}")
    months = [k[:-len("__vol_ratio")] for k in row if k.startswith("month_") and k.endswith("__vol_ratio")]
    shown = [m for m in months if max(np.nan_to_num(row[f"{m}__vol_gauge_Mm3"]), np.nan_to_num(row[f"{m}__vol_model_Mm3_matched"])) >= 0.005]
    if shown:
        L.append("  by month (months with under 0.005 M m3 on both sides are left out; the CSV has them all):")
        for m in shown:
            extra = f"   rain {fmt(row.get(f'{m}__rain_mm'), 1)} mm" if f"{m}__rain_mm" in row else ""
            L.append(f"    {m[6:].replace('_', '-'):<10}{fmt(row[f'{m}__vol_model_Mm3_matched']):>9}{fmt(row[f'{m}__vol_gauge_Mm3']):>9}{fmt(row[f'{m}__vol_ratio'], 3):>8}{extra}")
    L.append("")
    evs = list(profile.get("events", {}))
    if evs:
        L.append("EVENTS  (peak error %, first-arrival error h, KGE_2012 with r, beta, gamma, NSE; arrival level = 20% of the gauge peak)")
        L.append(f"  {'event':<18}{'peak model':>11}{'peak gauge':>11}{'peak err %':>11}{'arrival h':>10}{'KGE':>7}{'r':>6}{'beta':>7}{'gamma':>7}{'NSE':>7}   arrival note")
        for n in evs:
            L.append(f"  {n:<18}{fmt(row[f'{n}__peak_model_m3s'], 1):>11}{fmt(row[f'{n}__peak_gauge_m3s'], 1):>11}{fmt(row[f'{n}__peak_err_pct'], 1):>11}"
                     f"{fmt(row[f'{n}__arrival_err_h'], 2):>10}{fmt(row[f'{n}__kge2012']):>7}{fmt(row[f'{n}__kge_r']):>6}{fmt(row[f'{n}__kge_beta']):>7}"
                     f"{fmt(row[f'{n}__kge_gamma']):>7}{fmt(row[f'{n}__nse']):>7}   {row[f'{n}__arrival_flag']}")
        L.append("")
    tl = list(profile.get("tails", {}))
    if tl:
        L.append("TAILS  (volume M m3; recession rate k per hour = slope of ln(flow); ratio 1 = same recession, above 1 = model dries faster)")
        L.append(f"  {'tail':<18}{'model':>8}{'gauge':>8}{'ratio':>7}{'k model':>10}{'k gauge':>10}{'k ratio':>9}{'hours used (model/gauge)':>28}")
        for n in tl:
            L.append(f"  {n:<18}{fmt(row[f'{n}__vol_model_Mm3']):>8}{fmt(row[f'{n}__vol_gauge_Mm3']):>8}{fmt(row[f'{n}__tail_vol_ratio']):>7}"
                     f"{fmt(row[f'{n}__k_model_per_h'], 4):>10}{fmt(row[f'{n}__k_gauge_per_h'], 4):>10}{fmt(row[f'{n}__recession_rate_ratio']):>9}"
                     f"{str(row[f'{n}__n_hours_model']) + ' / ' + str(row[f'{n}__n_hours_gauge']):>28}")
        L.append("")
    ff = [n for n in profile.get("false_flow", {}) if f"{n}__vol_model_Mm3" in row]
    if ff:
        thr = profile.get("false_flow_threshold_m3s")
        L.append(f"FALSE FLOW  (the gauge is dry here; threshold {thr} m3/s = the gauge's Oct-Dec maximum per the design)")
        L.append(f"  {'window':<10}{'model M m3':>12}{'gauge M m3':>12}{'model max':>11}{'gauge max':>11}{'days above (model)':>20}{'(gauge)':>9}")
        for n in ff:
            L.append(f"  {n:<10}{fmt(row[f'{n}__vol_model_Mm3']):>12}{fmt(row[f'{n}__vol_gauge_Mm3'], 3):>12}{fmt(row[f'{n}__max_model_m3s'], 3):>11}"
                     f"{fmt(row[f'{n}__max_gauge_m3s'], 3):>11}{fmt(row[f'{n}__days_above_model'], 0):>20}{fmt(row[f'{n}__days_above_gauge'], 0):>9}")
        L.append("")
    rr_names = [n for n in {**profile.get("events", {}), **profile.get("volume_windows", {})} if f"{n}__runoff_ratio_model" in row]
    if rr_names:
        L.append(f"STORM BALANCE  (runoff ratio = outflow volume / rain volume over {BASIN_KM2:g} km2, rain = basin mean from the .mrf; report only, no pass/fail yet)")
        L.append(f"  {'window':<18}{'rain mm':>9}{'model':>8}{'gauge':>8}")
        for n in rr_names:
            L.append(f"  {n:<18}{fmt(row[f'{n}__rain_mm'], 1):>9}{fmt(row[f'{n}__runoff_ratio_model'], 3):>8}{fmt(row[f'{n}__runoff_ratio_gauge'], 3):>8}")
        for late, early in profile.get("balance_pairs", []):
            k = f"balance__{late}_over_{early}"
            L.append(f"  {late} / {early}:  model {fmt(row.get(k + '__model'))}   gauge {fmt(row.get(k + '__gauge'))}")
        L.append("")
    proc = {k: v for k, v in row.items() if k.startswith("proc__")}
    if proc:
        L.append("PROCESS CHECKS  (not scores: is the flood made by a sensible mechanism, and does it change with the start state?)")
        for k, v in proc.items():
            L.append(f"  {k[6:]:<32}{fmt(v, 3)}")
        L.append("")
    L.append(f"KGE_2012 over the span: {fmt(row['span__kge2012'], 3)}  (r {fmt(row['span__kge_r'])}, beta {fmt(row['span__kge_beta'])}, gamma {fmt(row['span__kge_gamma'])});  "
             f"NSE {fmt(row['span__nse'])} (sanity only);  PBIAS {fmt(row['span__pbias'], 1)} %.   VOLUME GATE for reading KGE: {row['kge_volume_gate']}")
    L.append("")
    L.append("TOLERANCE CHECK  (the design's PROPOSED tolerances; OK / MISS / n.a.)")
    for grp, label, val, ok, note in checks:
        status = "n.a." if ok is None else ("OK  " if ok else "MISS")
        L.append(f"  {status}  {grp:<12}{label:<46}{fmt(val, 3):>9}   ({note})")
    ref = row.get("_reference", [])
    if ref:
        L.append("")
        L.append("REFERENCE CHECKS  (numbers quoted in the design doc)")
        L.append("  'gauge' lines depend only on the gauge file and my window choices: if one does not match, my window differs from the one the design used.")
        L.append("  'baseline' lines match only when this run IS the Oct 21 baseline (Ks x0.2, PsiB x0.5, f x1, cv 6, flowexp 0.3, n 0.04, loss off, bedrock +1 m).")
        for label, value, design, kind, extra in ref:
            both = f"{fmt(value, 2)} | {fmt(extra, 2)}" if extra is not None else fmt(value, 2)
            near_ = np.isfinite(value) and abs(value - design) <= 0.01 * max(1.0, abs(design))
            if kind == "gauge":
                verdict = "matches the design to within 1%" if near_ else "DIFFERS from the design by more than 1%: check the window"
            else:
                verdict = "matches the design baseline to within 1%" if near_ else "differs from the design baseline (expected unless this run is the baseline)"
            L.append(f"  [{kind}] {label}: {both}   design {design:.2f}   {verdict}")
    return L


# ================================================================================================================
# CSV upsert
# ================================================================================================================
def upsert_csv(path: Path, rows: list):
    path.parent.mkdir(parents=True, exist_ok=True)
    new = pd.DataFrame(rows)
    if path.exists():
        old = pd.read_csv(path)
        old = old[~old["run_id"].isin(new["run_id"])]
        new = pd.concat([old, new], ignore_index=True, sort=False)
    first = ["run_id", "window", "scored_at", "scorer_version", "profile", "kge_volume_gate"]
    ordered = [c for c in first if c in new.columns] + [c for c in new.columns if c not in first]
    new[ordered].to_csv(path, index=False)


def flat(row):
    return {k: (v.isoformat(sep=" ") if isinstance(v, pd.Timestamp) else v) for k, v in row.items() if not k.startswith("_")}


# ================================================================================================================
# the `score` command
# ================================================================================================================
def collect_runs(args):
    runs = []
    if args.runs_file:
        df = pd.read_csv(args.runs_file, comment="#", skipinitialspace=True, dtype=str, keep_default_na=False)
        df.columns = [c.strip().lower() for c in df.columns]
        if "run" not in df.columns:
            die(f"{args.runs_file} needs a column named run (optional: window, start); found {list(df.columns)}")
        for r in df.itertuples():
            runs.append((r.run, (getattr(r, "window", "") or "").upper() or None, getattr(r, "start", "") or None))
    if args.run:
        runs += [(n.strip(), args.window.upper() if args.window else None, args.start) for n in args.run.split(",") if n.strip()]
    if args.scan:
        d = Path(args.scan)
        if not d.is_dir():
            die(f"--scan folder not found: {d}")
        for p in sorted(d.glob("**/*_Outlet.qout")):
            rid = p.name[:-len("_Outlet.qout")]
            runs.append((rid, args.window.upper() if args.window else None, args.start))
    if not runs:
        die("nothing to score: give --run, --runs-file or --scan")
    seen, out = set(), []
    for rid, w, st in runs:
        if rid not in seen:
            seen.add(rid)
            out.append((rid, w, st))
    return out


def score_one(run_id, window, start, args, profile, gauge_q, results_dir, write_files=True):
    window = window or run_window_from_id(run_id)
    start = start or (RUN_WINDOWS[window][0] if window in RUN_WINDOWS else None)
    if window and window not in RUN_WINDOWS:
        die(f"unknown window {window!r}; choose from {sorted(RUN_WINDOWS)}")
    if start is None:
        die(f"run {run_id!r}: cannot tell when it starts. Give --window (W1 to W5 or V), start the run id with W1_ ... V_, or give --start.")
    qpath, mpath, rpath, err = find_run_files(run_id, results_dir)
    if qpath is None:
        return None, err
    notes = []
    cs_state, _ = checksum_state(qpath, args.write_checksum)
    if cs_state == "MISMATCH":
        notes.append(f"CHECKSUM MISMATCH: {qpath.name} differs from its .sha256 sidecar; the kept file may have changed")
    model_q = read_qout(qpath, start)
    if model_q.empty:
        return None, f"{qpath} has no readable rows"
    exp_end = ts(RUN_WINDOWS[window][1]) if window in RUN_WINDOWS else None
    if exp_end is not None and model_q.index[-1] < exp_end - pd.Timedelta(hours=COVER_TOL_H):
        notes.append(f"run ends {model_q.index[-1]:%b %d %Y %H:%M}, before the {window} window end {exp_end:%b %d %Y %H:%M}: metrics beyond that are blank")
    if (model_q < 0).any():
        notes.append(f"{int((model_q < 0).sum())} negative flows in the outlet file")
    mrf = read_mrf(mpath, start) if mpath else None
    if mrf is None:
        notes.append("no .mrf file found: runoff ratios and water-table checks are blank")
    rft_df, rft_cols = (read_rft(rpath, start) if rpath else (None, None))
    if rpath and rft_df is None:
        notes.append(f".rft found but no Hsrf/Sbsrf columns; columns seen: {rft_cols}")
    elif rpath is None:
        notes.append("no .rft file found: Hsrf / Sbsrf shares are blank")
    row, checks = score_run(run_id, model_q, gauge_q, mrf, rft_df, rft_cols, profile, args.basin_km2)
    row.update({"run_id": run_id, "window": window or "", "scored_at": pd.Timestamp.now().isoformat(timespec="seconds"),
                "scorer_version": VERSION, "profile": args.profile if not args.profile_file else args.profile_file,
                "run_start": start, "qout_file": str(qpath), "mrf_file": str(mpath) if mpath else "", "rft_file": str(rpath) if rpath else ""})
    lines = report_lines(run_id, window, row, checks, profile, cs_state, notes)
    return (row, checks, lines), ""


def cmd_score(args) -> int:
    profile = load_profile(args.profile, args.profile_file)
    gauge_q = read_obs(args.obs)
    results_dir = Path(args.results_dir or args.scan or ".")
    runs = collect_runs(args)
    out_dir = Path(args.out_dir)
    print(f"score_run_sycamore.py v{VERSION}   gauge {args.obs}: {gauge_q.index[0]:%b %d %Y} to {gauge_q.index[-1]:%b %d %Y} ({int(gauge_q.notna().sum())} values)")
    print(f"  profile: {args.profile if not args.profile_file else args.profile_file};  runs: {len(runs)};  results folder: {results_dir}")
    if args.dry_run:
        for rid, w, st in runs:
            w2 = w or run_window_from_id(rid)
            q, m, r, err = find_run_files(rid, results_dir)
            start = st or (RUN_WINDOWS[w2][0] if w2 in RUN_WINDOWS else None)
            print(f"\n  run {rid}: window {w2 or '?'}, start {start or '?'}")
            if q is None:
                print(f"    PROBLEM: {err}")
                continue
            mq = read_qout(q, start) if start else None
            print(f"    outlet file {q}  ->  {'' if mq is None else f'{len(mq)} rows, {mq.index[0]:%b %d %Y %H:%M} to {mq.index[-1]:%b %d %Y %H:%M}'}")
            print(f"    .mrf {m if m else 'NOT FOUND'};  .rft {r if r else 'NOT FOUND'}")
            if mq is not None:
                M = Flow(mq)
                print("    windows this run covers fully:")
                allw = {"span": profile["period"], **profile["events"], **profile["volume_windows"], **profile["tails"], **profile["false_flow"], **profile["named_volume"]}
                for n, (a, b) in allw.items():
                    print(f"      {n:<18} {ts(a):%b %d %Y %H:%M} to {ts(b):%b %d %Y %H:%M}: {'covered' if covers(M, a, b) else 'not covered (metrics blank)'}")
        print("\nDRY RUN: nothing was scored or written.")
        return 0
    rows, all_lines, table = [], [], []
    for rid, w, st in runs:
        res, err = score_one(rid, w, st, args, profile, gauge_q, results_dir)
        if res is None:
            print(f"\n  SKIPPED {rid}: {err}")
            continue
        row, checks, lines = res
        rows.append(flat(row))
        all_lines += lines + ["", "-" * 100, ""]
        if args.quiet:
            print(f"  scored {rid}")
        else:
            print("\n" + "\n".join(lines))
        table.append((rid, row))
    if not rows:
        die("no run could be scored")
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / args.csv_name
    upsert_csv(csv_path, rows)
    (out_dir / "score_report.txt").write_text("\n".join(all_lines))
    if len(table) > 1:
        print("\nSUMMARY OF THIS BATCH (volume gate state; KGE_2012 only meaningful when the gate is open)")
        print(f"  {'run':<40}{'flood vol':>10}{'flood peak %':>13}{'arrival h':>10}{'oct-mar vol':>12}{'tail 1':>8}{'tail 2':>8}{'octdec M m3':>12}{'KGE span':>9}  gate")
        for rid, row in table:
            print(f"  {rid:<40}{fmt(row.get('flood_week__vol_ratio'), 3):>10}{fmt(row.get('flood_week__peak_err_pct'), 1):>13}{fmt(row.get('flood_week__arrival_err_h')):>10}"
                  f"{fmt(row.get('oct_mar__vol_ratio'), 3):>12}{fmt(row.get('tail_jan26_feb2__tail_vol_ratio')):>8}{fmt(row.get('tail_mar3_mar11__tail_vol_ratio')):>8}"
                  f"{fmt(row.get('octdec__vol_model_Mm3')):>12}{fmt(row.get('span__kge2012'), 3):>9}  {row['kge_volume_gate']}")
    print(f"\nWROTE {csv_path} ({len(rows)} run(s) added or replaced) and {out_dir / 'score_report.txt'}")
    return 0


# ================================================================================================================
# the `fingerprint` command (Stage 7)
# ================================================================================================================
def cumulative(t_s, q):
    return np.concatenate([[0.0], np.cumsum(0.5 * (q[1:] + q[:-1]) * np.diff(t_s))])


def fingerprints(ref_q: pd.Series, test_q: pd.Series, a, b, step_min=15, dry_thr=1e-6):
    """Series 122 fingerprints of `test` relative to `ref` over [a, b]: ref = the loss-off run, test = the loss-on run or the gauge.
    Both are put on one regular grid (linear interpolation)."""
    grid = pd.date_range(ts(a), ts(b), freq=f"{step_min}min")
    def on_grid(s):
        s = s.dropna()
        both = s.index.union(grid)
        return s.reindex(both).interpolate(method="time", limit_area="inside").reindex(grid)
    r, t = on_grid(ref_q), on_grid(test_q)
    ok = r.notna() & t.notna()
    r, t = r[ok], t[ok]
    if len(r) < 8:
        return None
    ts_s = r.index.values.astype("datetime64[s]").astype("int64")
    rv, tv = r.to_numpy(float), t.to_numpy(float)
    cr, ct = cumulative(ts_s, rv), cumulative(ts_s, tv)
    V_ref, V_test = cr[-1], ct[-1]
    out = {"n_points": int(len(r)), "grid_start": r.index[0], "grid_end": r.index[-1]}
    if V_ref <= 0:
        return out
    ipk = int(np.argmax(rv))
    out["ref_peak_m3s"], out["test_peak_m3s"] = float(rv[ipk]), float(tv.max())
    out["ref_peak_time"] = r.index[ipk]
    out["volume_change_pct"] = (V_test / V_ref - 1) * 100
    out["peak_change_pct"] = (tv.max() / rv[ipk] - 1) * 100
    out["peak_per_volume"] = out["peak_change_pct"] / out["volume_change_pct"] if abs(out["volume_change_pct"]) >= 0.5 else np.nan
    gap_total = V_ref - V_test
    out["gap_total_Mm3"] = gap_total / 1e6
    gap_cum = cr - ct
    out["share_of_gap_by_ref_peak"] = float(gap_cum[ipk] / gap_total) if abs(gap_total) > 0 else np.nan
    out["ref_own_volume_share_by_peak"] = float(cr[ipk] / V_ref)
    # phases cut by the reference: rise (to the peak), recession (until the reference is under 1% of its peak), tail
    below = np.nonzero((rv < 0.01 * rv[ipk]) & (np.arange(len(rv)) > ipk))[0]
    i_end = int(below[0]) if len(below) else len(rv) - 1
    if abs(gap_total) > 0:
        out["gap_share_rise"] = float(gap_cum[ipk] / gap_total)
        out["gap_share_recession"] = float((gap_cum[i_end] - gap_cum[ipk]) / gap_total)
        out["gap_share_tail"] = float((gap_cum[-1] - gap_cum[i_end]) / gap_total)
    # flow-band ratios, split by limb
    bands = {"5-20%": (0.05, 0.20), "20-60%": (0.20, 0.60), "60-100%": (0.60, 1.0001)}
    for limb, idx in (("rising", np.arange(len(rv)) <= ipk), ("recession", np.arange(len(rv)) > ipk)):
        for bname, (lo, hi) in bands.items():
            sel = idx & (rv >= lo * rv[ipk]) & (rv < hi * rv[ipk])
            out[f"ratio_{limb}_{bname}"] = float(np.median(tv[sel] / rv[sel])) if sel.sum() else np.nan
            out[f"n_{limb}_{bname}"] = int(sel.sum())
    # dry spell at the end
    for lab, v in (("ref", rv), ("test", tv)):
        wet = np.nonzero(v > dry_thr)[0]
        last_wet = int(wet[-1]) if len(wet) else -1
        out[f"{lab}_dry_hours_at_end"] = float((len(v) - 1 - last_wet) * step_min / 60.0)
        out[f"{lab}_dry_from"] = r.index[min(last_wet + 1, len(v) - 1)] if last_wet < len(v) - 1 else pd.NaT
    return out


def fingerprint_reading(f):
    """The Series 122 checklist (SMF synthetic truths, no gauge or rain error): numbers are priors, not rules."""
    L = []
    def line(text):
        L.append("  " + text)
    ppv = f.get("peak_per_volume", np.nan)
    if np.isfinite(ppv):
        tag = "loss-like" if ppv < 0.7 else ("rain-error-like" if ppv > 0.9 else "in between")
        line(f"peak change / volume change = {ppv:.2f}   [SMF: 0.16-0.53 for a channel loss, about 1.0 for a rain error]  -> {tag}")
    else:
        line("peak change / volume change: not defined (volume change under 0.5%)")
    sh, own = f.get("share_of_gap_by_ref_peak", np.nan), f.get("ref_own_volume_share_by_peak", np.nan)
    if np.isfinite(sh):
        tag = "back-loaded (loss-like)" if sh < own else "front-loaded (rain-error-like)"
        line(f"share of the missing water gone by the reference peak = {sh:.0%} (the reference itself has passed {own:.0%} of its volume by then)   "
             f"[SMF: about 35% for a loss, 70-74% for a rain error]  -> {tag}")
    sm, pk = f.get("ratio_recession_5-20%", np.nan), f.get("ratio_rising_60-100%", np.nan)
    pk2 = f.get("ratio_recession_60-100%", np.nan)
    peak_band = np.nanmean([pk, pk2]) if (np.isfinite(pk) or np.isfinite(pk2)) else np.nan
    if np.isfinite(sm) and np.isfinite(peak_band):
        tag = "small flows lose a bigger share than the peak (loss-like)" if sm < peak_band - 0.02 else "small flows lose no more than the peak (rain-error-like)"
        line(f"test / reference on the low recession (5-20% of the peak) = {sm:.2f}; near the peak (60-100%) = {peak_band:.2f}  -> {tag}")
    else:
        line("flow-band comparison: not enough points in the bands")
    dt, dr = f.get("test_dry_hours_at_end", np.nan), f.get("ref_dry_hours_at_end", np.nan)
    line(f"flow stays at or under the dry threshold for the last {dt:.1f} h in the test and {dr:.1f} h in the reference   "
         f"[SMF: a channel loss leaves the outlet at exactly zero for hours after the recession; a rain error does not]")
    line("These thresholds come from synthetic SMF truths with no gauge or rain error. At Sycamore, treat them as a prompt to look at the hydrograph, not a verdict.")
    return L


def cmd_fingerprint(args) -> int:
    profile = load_profile(args.profile, args.profile_file)
    results_dir = Path(args.results_dir or ".")
    def load(name):
        if name.lower() == "gauge":
            return read_obs(args.obs), "gauge"
        w = (args.ref_window if name == args.ref else args.test_window) or run_window_from_id(name)
        st = (args.ref_start if name == args.ref else args.test_start) or (RUN_WINDOWS[w][0] if w in RUN_WINDOWS else None)
        if st is None:
            die(f"run {name!r}: give its window or start (--ref-window/--ref-start or --test-window/--test-start, or start the run id with W1_ ... V_)")
        q, m, r, err = find_run_files(name, results_dir)
        if q is None:
            die(err)
        return read_qout(q, st), name
    ref_q, _ = load(args.ref)
    test_q, _ = load(args.test)
    if args.event:
        if args.event not in profile["events"]:
            die(f"unknown event {args.event!r}; events in the profile: {sorted(profile['events'])}")
        a, b = profile["events"][args.event]
    elif args.start_time and args.end_time:
        a, b = args.start_time, args.end_time
    else:
        die("give --event (a window name from the profile) or both --start-time and --end-time")
    f = fingerprints(ref_q, test_q, a, b, args.step_min, args.dry_threshold)
    if f is None or "volume_change_pct" not in f:
        die("not enough overlapping flow in that window to compute fingerprints")
    print(f"KCHAN FINGERPRINTS   reference (loss OFF): {args.ref}   test: {args.test}   window {ts(a):%b %d %Y %H:%M} to {ts(b):%b %d %Y %H:%M}")
    print(f"  grid {args.step_min} min, {f['n_points']} points; reference peak {f['ref_peak_m3s']:.2f} m3/s at {f['ref_peak_time']:%b %d %H:%M}, test peak {f['test_peak_m3s']:.2f} m3/s")
    print(f"  volume change {f['volume_change_pct']:+.2f}%   peak change {f['peak_change_pct']:+.2f}%   missing volume {f['gap_total_Mm3']:.3f} M m3")
    if "gap_share_rise" in f:
        print(f"  where the missing water is: rise {f['gap_share_rise']:.0%}, recession {f['gap_share_recession']:.0%}, long tail {f['gap_share_tail']:.0%}")
    print("  test / reference by flow level (median; n points):")
    for limb in ("rising", "recession"):
        print("    " + limb.ljust(10) + "  ".join(f"{b}: {f[f'ratio_{limb}_{b}']:.2f} (n {f[f'n_{limb}_{b}']})" if np.isfinite(f[f'ratio_{limb}_{b}']) else f"{b}: - (n {f[f'n_{limb}_{b}']})" for b in ("5-20%", "20-60%", "60-100%")))
    print("\nREADING (Series 122 checklist)")
    print("\n".join(fingerprint_reading(f)))
    if args.out_csv:
        p = Path(args.out_csv)
        p.parent.mkdir(parents=True, exist_ok=True)
        rec = {"ref": args.ref, "test": args.test, "window_start": ts(a), "window_end": ts(b), **{k: (v.isoformat(sep=" ") if isinstance(v, pd.Timestamp) else v) for k, v in f.items()}}
        pd.DataFrame([rec]).to_csv(p, mode="a", header=not p.exists(), index=False)
        print(f"\nappended to {p}")
    return 0


# ================================================================================================================
# self-test: synthetic gauge and synthetic runs with known answers
# ================================================================================================================
def _pulse(t_h, t0_h, peak, rise_h, k_per_h):
    """A flood: linear rise over rise_h hours, then exponential recession with rate k_per_h."""
    x = t_h - t0_h
    return np.where(x < 0, 0.0, np.where(x <= rise_h, peak * x / rise_h, peak * np.exp(-k_per_h * (x - rise_h))))


def _synthetic_world(tmp: Path):
    """Gauge: 15-min flows for WY2017 with 4 storms. Returns paths and the truth hydrograph function."""
    t0 = pd.Timestamp("2016-10-01 00:00")
    idx = pd.date_range(t0, "2017-10-01 00:00", freq="15min", inclusive="left")
    h = np.asarray((idx - t0) / HOUR, dtype=float)

    def hours_from(s):
        return float((pd.Timestamp(s) - t0) / HOUR)

    q = np.zeros(len(idx))
    storms = [("2017-01-01 06:00", 5.0, 2, 0.10),      # a small storm inside the gauge outage (Jan 1-2)
              ("2017-01-21 03:00", 220.0, 6, 0.045),   # flood: peak 220 m3/s, slow recession (tail carries water)
              ("2017-01-23 18:00", 18.0, 3, 0.08),
              ("2017-02-18 14:00", 40.0, 3, 0.05),
              ("2017-02-27 20:00", 80.0, 3, 0.04)]
    for s, pk, rise, k in storms:
        q += _pulse(h, hours_from(s), pk, rise, k)
    q[(idx >= "2016-11-10") & (idx < "2016-11-12")] += 0.10      # small Oct-Dec flow: the gauge's Oct-Dec maximum becomes 0.10 m3/s
    q = np.where(q < 0.003, 0.0, q)
    obs = pd.DataFrame({"datetime": idx, "q_cms": q, "q_cfs": q * 35.3147, "tier": "A", "note": ""})
    obs.loc[(idx >= "2017-01-01") & (idx < "2017-01-03"), "q_cms"] = np.nan          # the outage
    obs_path = tmp / "obs.csv"
    obs.to_csv(obs_path, index=False)
    return {"idx": idx, "q": q, "t0": t0, "obs_path": obs_path, "storms": storms, "hours_from": hours_from}


def _write_run(tmp: Path, run_id, start, q_fn, rows_hours, rain_mm_fn=None, with_rft=True, md_fn=None):
    """Write <run>_Outlet.qout (+ .mrf, .rft) from a function of hours-since-start."""
    d = tmp / "results"
    d.mkdir(exist_ok=True)
    th = np.arange(0.0, rows_hours + 1e-9, 0.25)
    th = np.sort(np.concatenate([th, th[:-1] + 0.1]))          # irregular steps, like the model's own
    pd.DataFrame({"Time_hr": th, "Q": q_fn(th)}).to_csv(d / f"{run_id}_Outlet.qout", index=False)
    hr = np.arange(1.0, rows_hours + 1e-9, 1.0)                # .mrf: hourly, labeled by hour END
    rain = rain_mm_fn(hr) if rain_mm_fn else np.zeros(len(hr))
    mrf = pd.DataFrame({"Time_hr": hr, "RainMin_mm_hr": 0.0, "MAP_mm_hr": rain, "Srf_m3_s": 0.0, "ChannelPerc_m3": 0.0,
                        "MDGW_mm": md_fn(hr) if md_fn else 1300.0, "AvSWE_cm": 0.0})
    mrf = pd.concat([mrf, pd.DataFrame({"Time_hr": [hr[-1] + 1], "RainMin_mm_hr": 9999.0, "MAP_mm_hr": 0.0, "Srf_m3_s": 0.0,
                                         "ChannelPerc_m3": 0.0, "MDGW_mm": 0.0, "AvSWE_cm": 0.0}) ], ignore_index=True)
    mrf.to_csv(d / f"{run_id}001.mrf", index=False)
    if with_rft:
        pd.DataFrame({"Time_hr": hr, "Hsrf": 3.0, "Sbsrf": 1.0}).to_csv(d / f"{run_id}001.rft", index=False)
    return d


def do_self_test() -> int:
    ok = True

    def check(cond, text):
        nonlocal ok
        print(("  PASS  " if cond else "  FAIL  ") + text)
        ok = ok and bool(cond)

    def near(a, b, tol):
        return bool(np.isfinite(a) and abs(a - b) <= tol)

    print("SELF-TEST: synthetic gauge and runs with known answers (nothing here is real data)")
    tmp = Path(tempfile.mkdtemp(prefix="score_selftest_"))
    try:
        w = _synthetic_world(tmp)
        t0, hours_from = w["t0"], w["hours_from"]
        storms = w["storms"]
        total_h = float((pd.Timestamp("2017-10-01") - t0) / HOUR)

        def truth(th):
            q = np.zeros_like(th)
            for s, pk, rise, k in storms:
                q += _pulse(th, hours_from(s), pk, rise, k)
            q[(th >= hours_from("2016-11-10")) & (th < hours_from("2016-11-12"))] += 0.10
            return np.where(q < 0.003, 0.0, q)

        # --- run P: the perfect run, rain chosen so the runoff ratios are known
        def rain_fn(hr):
            r = np.zeros(len(hr))
            for s_, mm in (("2017-02-18 12:00", 20.0), ("2017-02-27 18:00", 20.0)):        # 20 mm in the hour that STARTS then
                r[np.isclose(hr, hours_from(s_) + 1.0)] = mm                                # the .mrf labels that hour by its END
            return r
        _write_run(tmp, "W5_perfect", t0, truth, total_h, rain_mm_fn=rain_fn)
        # --- run B: a "bad" run: peaks x0.85, flood volume up, one hour late, tails cut (faster recession), false flow in Oct-Dec
        def bad(th):
            q = np.zeros_like(th)
            for s, pk, rise, k in storms:
                kk = k * (1.6 if pk > 100 else 1.0)
                q += _pulse(th - 1.0, hours_from(s), pk * 0.85, rise, kk)                  # one hour late
            q = np.where(q < 0.003, 0.0, q)
            q[(th >= 0) & (th < hours_from("2017-01-01"))] += 0.5                           # false flow all Oct-Dec
            return q
        _write_run(tmp, "W5_bad", t0, bad, total_h, rain_mm_fn=rain_fn)
        # --- run S: a W1 run (only the flood week), the perfect hydrograph
        t1 = pd.Timestamp("2017-01-18 00:00")
        def truth_w1(th):
            return truth(th + float((t1 - t0) / HOUR))
        _write_run(tmp, "W1_perfect", t1, truth_w1, 168.0)

        class A:  # argument stand-in
            obs = str(w["obs_path"]); profile = "wy2017"; profile_file = None; basin_km2 = BASIN_KM2
            write_checksum = False
        gauge_q = read_obs(A.obs)
        prof = load_profile("wy2017", None)

        # 1. perfect W5 run
        res, err = score_one("W5_perfect", None, None, A, prof, gauge_q, tmp / "results")
        check(res is not None, f"W5_perfect found and scored ({err})")
        row, checks, lines = res
        check(near(row["flood_week__vol_ratio"], 1.0, 0.002), f"perfect run: flood-week volume ratio is 1 ({row['flood_week__vol_ratio']:.4f})")
        check(near(row["flood_week__peak_err_pct"], 0.0, 0.2), f"perfect run: peak error is 0 ({row['flood_week__peak_err_pct']:.3f} %)")
        check(near(row["flood_week__arrival_err_h"], 0.0, 0.1), f"perfect run: arrival error is 0 ({row['flood_week__arrival_err_h']:.3f} h)")
        check(near(row["span__kge2012"], 1.0, 0.002), f"perfect run: KGE over the span is 1 ({row['span__kge2012']:.4f})")
        check(near(row["tail_jan26_feb2__tail_vol_ratio"], 1.0, 0.01) and near(row["tail_mar3_mar11__tail_vol_ratio"], 1.0, 0.01),
              "perfect run: both tail volume ratios are 1")
        check(near(row["tail_jan26_feb2__recession_rate_ratio"], 1.0, 0.02), f"perfect run: recession rate ratio is 1 ({row['tail_jan26_feb2__recession_rate_ratio']:.3f})")
        check(near(row["octdec__vol_model_Mm3"], 0.10 * 2 * 86400 / 1e6, 0.002), f"perfect run: Oct-Dec false-flow volume equals the gauge's own small flow ({row['octdec__vol_model_Mm3']:.4f} M m3)")
        check(row["octdec__days_above_model"] == 0, "perfect run: no day above 0.13 m3/s in Oct-Dec")
        check(near(row["octdec__max_gauge_m3s"], 0.10, 1e-6), "gauge Oct-Dec maximum reported (0.10 m3/s in this synthetic gauge)")
        # Jan 1-2: gauge blank + excluded in both: the matched model volume must be below the raw model volume for the span
        check(row["span__vol_model_Mm3_raw"] > row["span__vol_model_Mm3_matched"], "model volume over the Jan 1-2 outage counts in 'raw' but not in the matched volume")
        check(near(row["span__vol_ratio"], 1.0, 0.002), f"perfect run: span volume ratio is 1 ({row['span__vol_ratio']:.4f}) even with the gauge outage")
        # runoff ratios: 20 mm over 424.6 km2 = 8.492 M m3 of rain each
        rv = 20.0 * BASIN_KM2 / 1000.0
        exp_rr_gauge = row["feb18_19_tail__vol_gauge_Mm3"] / rv
        check(near(row["feb18_19_tail__runoff_ratio_gauge"], exp_rr_gauge, 1e-6) and near(row["feb18_19_tail__rain_mm"], 20.0, 1e-6),
              f"runoff ratio uses the .mrf rain shifted to hour-start labels (rain {row['feb18_19_tail__rain_mm']:.1f} mm; ratio {row['feb18_19_tail__runoff_ratio_gauge']:.3f})")
        check(near(row["balance__feb27_28_tail_over_feb18_19_tail__model"], row["balance__feb27_28_tail_over_feb18_19_tail__gauge"], 0.01),
              "perfect run: storm balance model = gauge")
        check(near(row["proc__hsrf_share"], 0.75, 1e-9), "Hsrf share read from the .rft (3 of 4 = 0.75)")
        check(row["kge_volume_gate"] == "open", f"perfect run: KGE volume gate is open ({row['kge_volume_gate']})")
        check(all(c[3] in (True, None) for c in checks), "perfect run: no tolerance MISS")

        # 2. bad W5 run: everything must show up
        res, err = score_one("W5_bad", None, None, A, prof, gauge_q, tmp / "results")
        rowb, checksb, _ = res
        check(near(rowb["flood_week__peak_err_pct"], -15.0, 2.0), f"bad run: flood peak about -15% ({rowb['flood_week__peak_err_pct']:.1f} %)")
        exp_arr = 1.0 + 6 * 0.2 * (1 / 0.85 - 1)           # one hour late, plus a slower climb to a lower peak
        check(near(rowb["flood_week__arrival_err_h"], exp_arr, 0.1), f"bad run: arrival error is +{exp_arr:.2f} h ({rowb['flood_week__arrival_err_h']:.2f} h)")
        check(rowb["tail_jan26_feb2__tail_vol_ratio"] < 0.7, f"bad run: Jan 26-Feb 2 tail ratio under 0.7 ({rowb['tail_jan26_feb2__tail_vol_ratio']:.2f})")
        check(rowb["tail_jan26_feb2__recession_rate_ratio"] > 1.4, f"bad run: model recedes faster, rate ratio about 1.6 ({rowb['tail_jan26_feb2__recession_rate_ratio']:.2f})")
        check(near(rowb["octdec__vol_model_Mm3"], 0.5 * 92 * 86400 / 1e6 + 0.0, 0.25), f"bad run: Oct-Dec false flow about 0.5 m3/s over 92 days ({rowb['octdec__vol_model_Mm3']:.2f} M m3)")
        check(rowb["octdec__days_above_model"] >= 85, f"bad run: nearly every Oct-Dec day above 0.13 m3/s ({rowb['octdec__days_above_model']})")
        misses = [c for c in checksb if c[3] is False]
        check(len(misses) >= 3, f"bad run: tolerance table shows MISS on several metrics ({len(misses)})")
        check(rowb["kge_volume_gate"] in ("open", "CLOSED"), f"bad run: KGE volume gate evaluated ({rowb['kge_volume_gate']})")
        check(rowb["span__kge2012"] < row["span__kge2012"], "bad run scores a lower KGE than the perfect run")

        # 3. W1 run: only the flood week is covered
        res, err = score_one("W1_perfect", None, None, A, prof, gauge_q, tmp / "results")
        r1, c1, l1 = res
        check(near(r1["flood_week__vol_ratio"], 1.0, 0.002), "W1 run: flood-week volume ratio is 1")
        check(near(r1["flood_week__arrival_err_h"], 0.0, 0.1), "W1 run: arrival 0")
        check(not np.isfinite(r1["feb18_19__vol_ratio"]) and not np.isfinite(r1["tail_jan26_feb2__tail_vol_ratio"]), "W1 run: windows it does not cover are left blank")
        check(not np.isfinite(r1["span__kge2012"]) or r1["span__is_full_period"] is False, "W1 run: span is marked partial")
        check(not r1["span__is_full_period"], "W1 run: span flagged as not the whole scoring period")

        # 4. CSV: two runs, rescoring replaces rather than duplicates
        out = tmp / "out"
        out.mkdir()
        csv = out / "scores_wide.csv"
        upsert_csv(csv, [flat(row), flat({**rowb})])
        upsert_csv(csv, [flat(row)])
        df = pd.read_csv(csv)
        check(len(df) == 2 and sorted(df["run_id"]) == ["W5_bad", "W5_perfect"], f"CSV holds one row per run after a rescore ({len(df)} rows)")
        check(df.columns[0] == "run_id", "run_id is the first column")

        # 5. scan + command-line path (window from the run id), dry run, checksum sidecar
        ns = argparse.Namespace(obs=A.obs, profile="wy2017", profile_file=None, basin_km2=BASIN_KM2, write_checksum=True, runs_file=None, run=None,
                                window=None, start=None, scan=str(tmp / "results"), results_dir=str(tmp / "results"), out_dir=str(out), dry_run=False,
                                quiet=True, csv_name="scan_scores.csv")
        rc = cmd_score(ns)
        d2 = pd.read_csv(out / "scan_scores.csv")
        check(rc == 0 and sorted(d2["run_id"]) == ["W1_perfect", "W5_bad", "W5_perfect"], f"--scan scores every outlet file and reads the window from the run id ({sorted(d2['run_id'])})")
        side = list((tmp / "results").glob("*.sha256"))
        check(len(side) == 3, f"--write-checksum wrote a sidecar for each outlet file ({len(side)})")
        ns.write_checksum = False
        # damage one kept file: the checksum must complain
        qfile = tmp / "results" / "W1_perfect_Outlet.qout"
        qfile.write_text(qfile.read_text() + "999.0,1.0\n")
        st, _ = checksum_state(qfile, False)
        check(st == "MISMATCH", "an altered outlet file is caught by its checksum")
        ns.dry_run = True
        ns.quiet = False
        rc = cmd_score(ns)
        check(rc == 0, "--dry-run runs and writes nothing")
        # runs file with an explicit window and start, and a run id that does not begin with a window name
        shutil.copy(tmp / "results" / "W5_perfect_Outlet.qout", tmp / "results" / "mystery_run_Outlet.qout")
        shutil.copy(tmp / "results" / "W5_perfect001.mrf", tmp / "results" / "mystery_run001.mrf")
        (tmp / "runs.csv").write_text("run,window,start\nmystery_run,W5,\n")
        ns2 = argparse.Namespace(**{**vars(ns), "dry_run": False, "quiet": True, "scan": None, "run": None, "runs_file": str(tmp / "runs.csv"), "csv_name": "runs_scores.csv"})
        cmd_score(ns2)
        d3 = pd.read_csv(out / "runs_scores.csv")
        check(list(d3["run_id"]) == ["mystery_run"] and d3["window"].iloc[0] == "W5", "--runs-file with an explicit window works for a run id with no window prefix")

        # 6. fingerprints: a loss-like run and a rain-error-like run against the perfect run
        def loss_like(th):
            return np.maximum(truth(th) - 8.0, 0.0)             # a flat loss: takes up to 8 m3/s from every flow, so small flows lose the most
        def rain_like(th):
            return truth(th) * 0.9
        _write_run(tmp, "W5_loss", t0, loss_like, total_h, with_rft=False)
        _write_run(tmp, "W5_rain", t0, rain_like, total_h, with_rft=False)
        a_, b_ = prof["events"]["flood_week"]
        qp = read_qout(tmp / "results" / "W5_perfect_Outlet.qout", t0)
        fl = fingerprints(qp, read_qout(tmp / "results" / "W5_loss_Outlet.qout", t0), a_, b_)
        fr = fingerprints(qp, read_qout(tmp / "results" / "W5_rain_Outlet.qout", t0), a_, b_)
        check(near(fr["peak_per_volume"], 1.0, 0.05), f"fingerprint: a 10% rain error gives peak change / volume change = 1 ({fr['peak_per_volume']:.2f})")
        check(fl["peak_per_volume"] < 0.5, f"fingerprint: a flat loss gives peak change / volume change well under 1 ({fl['peak_per_volume']:.2f})")
        check(fr["share_of_gap_by_ref_peak"] > fr["ref_own_volume_share_by_peak"] - 0.02, f"fingerprint: rain error is not back-loaded ({fr['share_of_gap_by_ref_peak']:.0%} by the peak)")
        check(fl["share_of_gap_by_ref_peak"] < fr["share_of_gap_by_ref_peak"], f"fingerprint: the loss is more back-loaded than the rain error ({fl['share_of_gap_by_ref_peak']:.0%} vs {fr['share_of_gap_by_ref_peak']:.0%})")
        check(fl["ratio_recession_5-20%"] < fl["ratio_recession_60-100%"] if np.isfinite(fl["ratio_recession_60-100%"]) else True,
              "fingerprint: the loss takes a bigger share of small recession flows than of large ones")
        check(fl["test_dry_hours_at_end"] > fr["test_dry_hours_at_end"] or fl["test_dry_hours_at_end"] >= 1.0, f"fingerprint: loss run goes dry for hours ({fl['test_dry_hours_at_end']:.1f} h)")
        rc = cmd_fingerprint(argparse.Namespace(ref="W5_perfect", test="W5_loss", results_dir=str(tmp / "results"), event="flood_week", profile="wy2017",
                                                profile_file=None, obs=A.obs, ref_window=None, test_window=None, ref_start=None, test_start=None,
                                                start_time=None, end_time=None, step_min=15, dry_threshold=1e-6, out_csv=str(out / "fp.csv")))
        check(rc == 0 and (out / "fp.csv").exists(), "fingerprint command runs and appends a CSV row")
        rc = cmd_fingerprint(argparse.Namespace(ref="W5_perfect", test="gauge", results_dir=str(tmp / "results"), event="flood_week", profile="wy2017",
                                                profile_file=None, obs=A.obs, ref_window=None, test_window=None, ref_start=None, test_start=None,
                                                start_time=None, end_time=None, step_min=15, dry_threshold=1e-6, out_csv=None))
        check(rc == 0, "fingerprint against the gauge runs")
        # 7. dump-profile round trip
        pj = out / "profile.json"
        pj.write_text(json.dumps(PROFILES["wy2017"], indent=2))
        check(load_profile("wy2017", str(pj))["events"] == PROFILES["wy2017"]["events"], "profile JSON round trip")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("\nSELF-TEST " + ("PASSED" if ok else "FAILED"))
    return 0 if ok else 1


# ================================================================================================================
# command line
# ================================================================================================================
def build_parser():
    p = argparse.ArgumentParser(description="Standalone run scorer for the Lower Sycamore calibration design (see the text at the top of this file).",
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--obs", default=DEFAULT_OBS, help=f"gauge file (default {DEFAULT_OBS})")
    common.add_argument("--profile", default="wy2017", choices=sorted(PROFILES), help="window list: wy2017 (default) or wy2020")
    common.add_argument("--profile-file", help="JSON window list (from dump-profile) instead of a built-in profile")
    common.add_argument("--results-dir", help="folder that holds the kept run files (searched recursively)")
    s = sub.add_parser("score", parents=[common], help="score runs (also the rescore mode)")
    s.add_argument("--run", help="run id, or several separated by commas")
    s.add_argument("--runs-file", help="CSV with columns run[,window][,start]")
    s.add_argument("--scan", help="score every *_Outlet.qout under this folder")
    s.add_argument("--window", help="W1 W2 W3 W4 W5 or V (default: read from the run id)")
    s.add_argument("--start", help="run start, e.g. '2016-10-01 00:00' (default: the window's start)")
    s.add_argument("--out-dir", default="score_out")
    s.add_argument("--csv-name", default="scores_wide.csv")
    s.add_argument("--basin-km2", type=float, default=BASIN_KM2)
    s.add_argument("--write-checksum", action="store_true", help="write a .sha256 sidecar next to each outlet file (and verify nothing)")
    s.add_argument("--quiet", action="store_true", help="do not print each run's report (still written to score_report.txt)")
    s.add_argument("--dry-run", action="store_true", help="show which files and windows would be used; score nothing")
    f = sub.add_parser("fingerprint", parents=[common], help="Stage 7 Kchan fingerprints")
    f.add_argument("--ref", required=True, help="the loss-OFF run id")
    f.add_argument("--test", required=True, help="the loss-ON run id, or the word gauge")
    f.add_argument("--event", help="a window name from the profile (default window set: flood_week)")
    f.add_argument("--start-time", help="window start if not using --event")
    f.add_argument("--end-time", help="window end if not using --event")
    f.add_argument("--ref-window", help="W1..W5 or V for the reference run (default: from its run id)")
    f.add_argument("--test-window", help="W1..W5 or V for the test run (default: from its run id)")
    f.add_argument("--ref-start", help="start time of the reference run, if no window")
    f.add_argument("--test-start", help="start time of the test run, if no window")
    f.add_argument("--step-min", type=int, default=15)
    f.add_argument("--dry-threshold", type=float, default=1e-6, help="flow at or under this counts as dry (m3/s); use something like 0.01 for a gauge")
    f.add_argument("--out-csv", help="append the fingerprint numbers to this CSV")
    d = sub.add_parser("dump-profile", help="write a built-in profile to JSON for editing")
    d.add_argument("--profile", default="wy2017", choices=sorted(PROFILES))
    d.add_argument("--out", default="profile_wy2017.json")
    sub.add_parser("self-test", help="check every metric on synthetic runs")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    warnings.filterwarnings("ignore", category=FutureWarning)
    if args.cmd == "self-test":
        return do_self_test()
    if args.cmd == "dump-profile":
        Path(args.out).write_text(json.dumps(PROFILES[args.profile], indent=2) + "\n")
        print(f"wrote {args.out}")
        return 0
    if args.cmd == "fingerprint":
        return cmd_fingerprint(args)
    return cmd_score(args)


if __name__ == "__main__":
    sys.exit(main())
