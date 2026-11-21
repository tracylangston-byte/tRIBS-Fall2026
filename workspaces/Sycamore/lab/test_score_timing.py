#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_score_timing.py   (v1.0, 2026-11-16)
=========================================
Stand-in test of score_timing_sycamore.py: no tRIBS, no real data. Builds a fake lab folder in a temp directory with
synthetic gauge and outlet files whose volumes equal the volumes on record, then checks the formulas against answers
worked out independently, the stop rules, the blank-with-reason rule, and that no existing file is touched.

Run (from the folder that holds both scripts):   python test_score_timing.py
Mutation check (breaks a copy of the script on purpose, expects the test to notice):   python test_score_timing.py --mutations
"""
import contextlib
import hashlib
import importlib.util
import io
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SCRIPT = Path(os.environ.get("ST_SCRIPT", HERE / "score_timing_sycamore.py"))
FAILS = []
HOUR = pd.Timedelta(hours=1)


def check(cond, text):
    print(("  PASS  " if cond else "  FAIL  ") + text)
    if not cond:
        FAILS.append(text)


def load(path):
    spec = importlib.util.spec_from_file_location("st_under_test", path)
    m = importlib.util.module_from_spec(spec)
    sys.modules["st_under_test"] = m
    spec.loader.exec_module(m)
    return m


st = load(SCRIPT)


def hs(values, start="2020-01-01 00:00"):
    idx = pd.date_range(start, periods=len(values), freq="60min")
    return pd.Series(np.asarray(values, float), index=idx)


# ============================================================================================ unit tests
def unit_tests():
    print("\nUNIT TESTS (formulas against answers worked out by hand or by a different method)")
    lo = pd.Timestamp("2020-01-01 00:00")
    # 1. rectangular pulse: 1 m3/s for hours 10..19 of a 48 h window -> t10 = 11, t50 = 15, t90 = 19
    v = np.zeros(48)
    v[10:20] = 1.0
    tot, t, why = st.volume_times(hs(v), lo)
    check(abs(tot - 10 * 3600 / 1e6) < 1e-12, "volume of a 10 h rectangle = 0.036 M m3")
    check(t is not None and abs(t[0.1] - 11) < 1e-9 and abs(t[0.5] - 15) < 1e-9 and abs(t[0.9] - 19) < 1e-9,
          f"rectangle: t10, t50, t90 = 11, 15, 19 h (got {None if t is None else [round(t[k], 6) for k in (0.1, 0.5, 0.9)]})")
    # 2. a different method on a random hydrograph: inverse interpolation of the cumulative volume at hour boundaries
    rng = np.random.default_rng(7)
    v = np.abs(rng.normal(0, 1, 96)) * np.exp(-((np.arange(96) - 30) / 15.0) ** 2)
    tot, t, why = st.volume_times(hs(v), lo)
    cum = np.concatenate([[0.0], np.cumsum(v)])
    bounds = np.arange(97, dtype=float)
    alt = [float(np.interp(f * cum[-1], cum, bounds)) for f in (0.1, 0.5, 0.9)]
    check(t is not None and np.allclose([t[0.1], t[0.5], t[0.9]], alt, atol=1e-9),
          f"random hydrograph: matches inverse interpolation of the cumulative curve ({[round(x, 4) for x in alt]})")
    # 3. a gap: hours 20..29 missing add no volume; times still come out and the gap is found
    v = np.ones(48)
    s = hs(v).drop(hs(v).index[20:30])
    tot, t, why = st.volume_times(s, lo)
    check(abs(tot - 38 * 3600 / 1e6) < 1e-12, "gap of 10 h: volume counts only the 38 hours with data")
    runs = st.missing_runs(s, lo, lo + 48 * HOUR)
    check(len(runs) == 1 and (runs[0][1] - runs[0][0]) == 10 * HOUR, "missing_runs finds one 10 h gap")
    # 4. zero volume and empty
    tot, t, why = st.volume_times(hs(np.zeros(10)), lo)
    check(t is None and "zero" in why, "zero volume: times undefined, reason given")
    tot, t, why = st.volume_times(hs([]), lo)
    check(t is None and "no hourly" in why, "no data: times undefined, reason given")
    # 5. tail identity: the gauge tail after its own t90 is 10% of its volume
    v = np.abs(rng.normal(0, 1, 120)) * np.exp(-((np.arange(120) - 25) / 20.0) ** 2) + 0.01
    h = hs(v)
    tot, t, why = st.volume_times(h, lo)
    tail = st.volume_after(h, lo, lo + 120 * HOUR, lo + pd.Timedelta(hours=t[0.9]))
    check(abs(tail - 0.1 * tot) < 1e-9, f"tail after t90 is 10% of the volume ({tail:.6f} vs {0.1 * tot:.6f})")
    # 6. lag sign and size: model = gauge shifted LATER by 3 h -> lag +3; shifted earlier by 2 -> -2
    base = np.exp(-((np.arange(120) - 50) / 6.0) ** 2) + 0.4 * np.exp(-((np.arange(120) - 80) / 4.0) ** 2)
    g = hs(base)
    for shift in (3, -2, 0, 12):
        m = hs(np.roll(base, shift))
        lag, r, atl, why = st.best_lag(g, m)
        check(lag == shift and r > 0.999, f"best lag of a series shifted by {shift:+d} h is {shift:+d} (got {lag})")
    # 7. Pearson
    r, n, why = st.pearson(g, g * 3.0)
    check(abs(r - 1) < 1e-12 and n == 120, "Pearson r of a series with 3 x itself is 1")
    r, n, why = st.pearson(g, hs(np.ones(120)))
    check(r is None and "model" in why and "constant" in why, "Pearson r with a constant MODEL series is blank with a reason (not zero)")
    r, n, why = st.pearson(hs(np.ones(120)), g)
    check(r is None and "gauge" in why and "constant" in why, "Pearson r with a constant GAUGE series is blank with a reason (not zero)")
    # 8. second bump
    def bump(a, dipv, b, n=40):
        v = np.zeros(n)
        v[10] = a
        v[11:20] = np.linspace(a, dipv, 9)
        v[20:29] = np.linspace(dipv, b, 9)
        v[29] = b
        v[30:] = np.linspace(b, 0, n - 30)
        v[:10] = np.linspace(0, a, 10)
        return hs(v)
    check(st.second_bump(bump(10, 1, 9.5)), "two bumps (10 and 9.5, dip to 1): flagged")
    check(not st.second_bump(bump(10, 1, 7.0)), "second bump only 70% as high: not flagged")
    check(not st.second_bump(bump(10, 9.5, 9.0)), "second bump 90% as high but the dip is only 5%: not flagged")
    check(not st.second_bump(hs(np.exp(-((np.arange(60) - 30) / 8.0) ** 2))), "single smooth peak: not flagged")
    # 9. peak flags: edge and flat top
    idx = pd.date_range(lo, periods=48 * 4, freq="15min")
    q = pd.Series(np.exp(-(((idx - lo) / HOUR).to_numpy(float) - 3) ** 2 / 4.0), index=idx)
    pk, why = st.peak_info(q, st.hourly(q, lo, lo + 48 * HOUR), lo, lo + 48 * HOUR)
    check(pk is not None and "edge" in pk["flags"], "peak at hour 3 of the window is flagged as near the edge")
    q = pd.Series(np.exp(-(((idx - lo) / HOUR).to_numpy(float) - 24) ** 2 / 4.0), index=idx)
    pk, why = st.peak_info(q, st.hourly(q, lo, lo + 48 * HOUR), lo, lo + 48 * HOUR)
    check(pk is not None and pk["flags"] == "" and pk["time"] == lo + 24 * HOUR, "peak at hour 24: no flag, right time")
    q = pd.Series(np.where(((idx - lo) / HOUR).to_numpy(float).round(6).clip(0) % 48 < 30, 0.0, 0.0), index=idx)
    q[(idx >= lo + 20 * HOUR) & (idx < lo + 23 * HOUR)] = 5.0
    pk, why = st.peak_info(q, st.hourly(q, lo, lo + 48 * HOUR), lo, lo + 48 * HOUR)
    check(pk is not None and "flat top" in pk["flags"], "a 3 h plateau is flagged as flat-topped")
    # 10. arrival
    idx = pd.date_range(lo, periods=40, freq="15min")
    ramp = pd.Series(np.maximum(0, (np.arange(40) - 8) * 0.5), index=idx)             # zero until sample 8, then +0.5 per 15 min
    t, why = st.arrival(ramp, 2.0, lo, lo + 10 * HOUR)
    check(t is not None and abs(t - 3.0) < 1e-9, f"arrival through 2.0 on a ramp: 3.0 h (got {t})")
    t, why = st.arrival(ramp, 1e9, lo, lo + 10 * HOUR)
    check(t is None and "never" in why, "arrival never reached: blank with reason")
    t, why = st.arrival(ramp + 5, 2.0, lo, lo + 10 * HOUR)
    check(t is None and "already" in why, "arrival already above at the start: blank with reason")
    gap = ramp.drop(ramp.index[10:30])
    t, why = st.arrival(gap, 3.0, lo, lo + 10 * HOUR)
    check(t is None and "gap" in why, "first rise inside a reading gap: blank with reason")

    # 11. windows
    print("  -- windows")
    wa = {w.key: w for w in st.make_windows("as-scored")}
    wc = {w.key: w for w in st.make_windows("cut")}
    wu = {w.key: w for w in st.make_windows("uncut")}
    exp_hours_as = {("WY2017", "Dec 21-25"): 192, ("WY2017", "Dec 30-Jan 2"): 168, ("WY2017", "Jan 19-21"): 96,
                    ("WY2017", "Feb 18-19"): 120, ("WY2017", "Feb 27-28"): 120,
                    ("WY2020", "Nov 19-22"): 168, ("WY2020", "Nov 28-30"): 144, ("WY2020", "Dec 23-28"): 216,
                    ("WY2020", "Feb 10-11"): 120, ("WY2020", "Feb 22-23"): 120, ("WY2020", "Mar 10-13"): 168, ("WY2020", "Mar 18-19"): 120}
    check(all(wa[k].hours == h for k, h in exp_hours_as.items()), "as-scored window lengths (h) equal the storm_table / scorer lengths")
    check(wa[("WY2020", "Nov 28-30")].hi == pd.Timestamp("2019-12-04"), "as-scored: WY2020 Nov 28-30 is NOT cut (ends Dec 4)")
    check(wa[("WY2017", "Jan 19-21")].hi == pd.Timestamp("2017-01-23"), "as-scored: WY2017 Jan 19-21 is cut at Jan 23")
    check(wc[("WY2020", "Nov 28-30")].hi == pd.Timestamp("2019-12-03") and wu[("WY2017", "Jan 19-21")].hi == pd.Timestamp("2017-01-25"),
          "cut mode ends Nov 28-30 at Dec 3; uncut mode ends Jan 19-21 at Jan 25")
    diff_cut = [k for k in wa if wa[k].hi != wc[k].hi]
    diff_unc = [k for k in wa if wa[k].hi != wu[k].hi]
    check(diff_cut == [("WY2020", "Nov 28-30")] and diff_unc == [("WY2017", "Jan 19-21")],
          "the three window rules differ in exactly one storm each, as expected")
    print("  -- fmt")
    runs = [(lo, lo + 24 * HOUR)]
    check("24 h" in st.fmt_runs(runs), "fmt_runs prints the length of the gap")


# ============================================================================================ synthetic lab
RAIN = {  # (year, first day): basin_mm, peak_1h, peak_3h, peak_6h, top6h_share, t50_frac, dry_days, prev7, prev30
    ("WY2017", "2016-12-21"): ("47.803", "4.429", "11.1693", "16.69516", "0.46", "0.27", "3", "17.694", "41.998"),
    ("WY2017", "2016-12-30"): ("48.178", "4.909", "8.65767", "12.13796", "0.4", "0.38", "5", "19.662", "65.503"),
    ("WY2017", "2017-01-19"): ("63.4", "6.906", "20.38584", "32.03875", "0.51", "0.62", "3", "16.872", "115.89"),
    ("WY2017", "2017-02-18"): ("46.41", "5.23", "9.94524", "17.18367", "0.43", "0.47", "24", "1.35", "83.018"),
    ("WY2017", "2017-02-27"): ("41.5", "6.596", "13.63578", "17.95637", "0.6", "0.59", "7", "0.253", "48.025"),
    ("WY2020", "2019-11-19"): ("71.545", "5.658", "13.88666", "20.96068", "0.36", "0.35", "53", "0.0", "0.765"),
    ("WY2020", "2019-11-28"): ("44.039", "6.489", "16.88779", "22.87426", "0.55", "0.41", "6", "17.798", "72.476"),
    ("WY2020", "2019-12-23"): ("75.196", "3.595", "7.62042", "10.81779", "0.22", "0.41", "14", "0.0", "70.22"),
    ("WY2020", "2020-02-10"): ("30.319", "4.346", "9.62573", "14.16476", "0.64", "0.51", "19", "0.0", "9.103"),
    ("WY2020", "2020-02-22"): ("52.637", "4.83", "12.90156", "24.43584", "0.46", "0.24", "10", "0.0", "30.643"),
    ("WY2020", "2020-03-10"): ("64.518", "4.831", "9.98611", "13.89467", "0.3", "0.39", "7", "0.0", "86.525"),
    ("WY2020", "2020-03-18"): ("24.61", "4.363", "11.09771", "17.56747", "0.71", "0.36", "4", "51.502", "120.724"),
}
SHIFT = {"retuned": 2, "original": 1, "kchan05": 1, "kchan02": 1}      # model later than gauge by this many hours (others 3)
OUTAGES = {"WY2017": [("2017-01-01 00:00", "2017-01-02 00:00")],
           "WY2020": [("2019-12-02 00:00", "2019-12-17 00:00"), ("2019-12-27 00:00", "2019-12-30 00:00")]}
ENDS = {"WY2017": "2017-04-01 00:00", "WY2020": "2020-06-01 00:00"}
PULSE_OFFSET, RISE, FALL = 18, 6, 24


def pulse(idx, t0, rise, fall):
    x = ((idx - t0) / HOUR).to_numpy(float)
    y = np.where((x >= 0) & (x < rise), x / rise, np.where((x >= rise) & (x < rise + fall), 1 - (x - rise) / fall, 0.0))
    return pd.Series(y, index=idx)


def scaled(unit, lo, hi, target, mask=None):
    s = unit.copy()
    if mask is not None:
        s[mask] = np.nan
    v = st.vol_mm3(st.hourly(s, lo, hi))
    return (unit * (target / v) if v > 0 else unit * 0.0), (target / v if v > 0 else 0.0)


def build_lab(tmp, with_tier=True, jitter=None):
    """Synthetic lab. Returns (lab path, info dict). jitter: {(year, key): factor} multiplies one outlet file."""
    base = Path(tmp) / "Sycamore"
    lab, res, obs = base / "lab", base / "lab" / "results", base / "init_data" / "obs"
    for d in (lab, res, obs):
        d.mkdir(parents=True, exist_ok=True)
    info = {"gauge_factor": {}, "rain": RAIN}
    wins = {w.key: w for w in st.make_windows("as-scored")}
    refs = st.REFS
    for yr in ("WY2017", "WY2020"):
        start = st.YEARS[yr]["start"]
        end = pd.Timestamp(ENDS[yr])
        idx = pd.date_range(start, end, freq="15min", inclusive="left")
        masks = pd.Series(False, index=idx)
        for a, b in OUTAGES[yr]:
            masks[(idx >= pd.Timestamp(a)) & (idx < pd.Timestamp(b))] = True
        gq = pd.Series(0.0, index=idx)
        for w in [x for x in wins.values() if x.year == yr]:
            unit = pulse(idx, w.lo + PULSE_OFFSET * HOUR, RISE, FALL)
            ref = [r for r in refs[w.key] if r[0] == w.hi][0][1]["gauge"][0]
            sc, f = scaled(unit, w.lo, w.hi, ref, masks)
            info["gauge_factor"][w.key] = f
            gq += sc
        if yr == "WY2017":        # the Jan 23-24 storm that only the uncut Jan 19-21 window contains
            unit = pulse(idx, pd.Timestamp("2017-01-23 12:00"), 4, 14)
            w = wins[("WY2017", "Jan 19-21")]
            sc, _ = scaled(unit, pd.Timestamp("2017-01-23"), pd.Timestamp("2017-01-25"), st.WY2017_JAN_UNCUT[0] - 5.5667)
            gq += sc
        tier = np.array(["measured"] * len(idx), dtype=object)
        gq_out = gq.copy()
        gq_out[masks] = np.nan
        tier[masks.to_numpy()] = "missing"
        if yr == "WY2017":        # known quality: in the Feb 27-28 window mark 20 readings estimated and 10 filled_interp
            w = wins[("WY2017", "Feb 27-28")]
            pos = np.nonzero(np.asarray((idx >= w.lo) & (idx < w.hi)))[0]
            tier[pos[100:120]] = "estimated"
            tier[pos[200:210]] = "filled_interp"
        df = pd.DataFrame({"datetime": idx.strftime("%Y-%m-%d %H:%M"), "q_cms": gq_out.to_numpy(), "q_cfs": gq_out.to_numpy() * st.CFS_PER_CMS})
        if with_tier:
            df["tier"] = tier
        name = "sycamore_09510200_obs_flow_WY2017.csv" if yr == "WY2017" else "sycamore_09510200_obs_flow_Sep2019_to_May2020.csv"
        df.to_csv(obs / name, index=False)
        # model files
        for key, label, stem in st.RUNS[yr]:
            step = 1 / 16 if (yr, key) in (("WY2020", "retuned"), ("WY2017", "original")) else 1 / 4
            mi = start + pd.to_timedelta(np.arange(0, (end - start) / HOUR, step), unit="h")
            shift = SHIFT.get(key, 3)
            mq = pd.Series(0.0, index=mi)
            for w in [x for x in wins.values() if x.year == yr]:
                unit = pulse(mi, w.lo + (PULSE_OFFSET + shift) * HOUR, RISE, FALL)
                rr = [r for r in refs[w.key] if r[0] == w.hi][0][1]
                if key in rr:
                    target = rr[key][0]
                elif key in st.WY2020_TOTALS:
                    target = st.WY2020_TOTALS[key] / 7.0
                else:
                    target = 1.0
                if (w.key == ("WY2017", "Dec 30-Jan 2")) and key == "retuned":
                    # part A on gauge hours (ratio 0.58 against the gauge), part B inside the outage
                    ua = pulse(mi, w.lo + 19 * HOUR, 3, 12)
                    sa, _ = scaled(ua, w.lo, w.hi, 0.58 * rr["gauge"][0])
                    ub = pulse(mi, pd.Timestamp("2017-01-01 06:00"), 3, 6)
                    sb, _ = scaled(ub, w.lo, w.hi, target - 0.58 * rr["gauge"][0])
                    mq += sa + sb
                    continue
                sc, _ = scaled(unit, w.lo, w.hi, target)
                mq += sc
            if yr == "WY2017":
                unit = pulse(mi, pd.Timestamp("2017-01-23 12:00"), 4, 14)
                extra = (st.WY2017_JAN_UNCUT[1] - 7.3129) if key == "retuned" else 0.3
                sc, _ = scaled(unit, pd.Timestamp("2017-01-23"), pd.Timestamp("2017-01-25"), extra)
                mq += sc
            if jitter and (yr, key) in jitter:
                mq = mq * jitter[(yr, key)]
            hrs = (mi - start) / HOUR
            out = pd.DataFrame({"Time_hr": hrs.to_numpy(float).round(4), "Q": mq.to_numpy(), "Qs": mq.to_numpy() * 0.5})
            out.to_csv(res / f"{stem}_Outlet.qout", index=False, float_format="%.6f")
    # storm table with the short headers
    rows = []
    cols = ["year", "event", "first_day", "last_day", "basin_mm", "peak_1h_mm", "peak_3h_mm", "peak_6h_mm", "top6h_share", "t50_frac", "dry_days_before", "prev7_mm", "prev30_mm"]
    for (y, n, f, l, nx) in st.STORMS:
        rows.append([y, n, f, l, *RAIN[(y, f)]])
    rows.append(["WY2020", "Dec 3-4", "2019-12-03", "2019-12-04", "5.067", "0.683", "1.68977", "2.68984", "", "", "3", "44.203", "116.515"])
    pd.DataFrame(rows, columns=cols).to_csv(lab / "storm_table.csv", index=False)
    readable = {"year": "Water year", "event": "Event", "first_day": "First rain day", "last_day": "Last rain day",
                "basin_mm": "Basin rain (mm)", "peak_1h_mm": "Peak 1 h (mm)", "peak_3h_mm": "Peak 3 h (mm)", "peak_6h_mm": "Peak 6 h (mm)",
                "top6h_share": "Top-6 h share (fraction)", "t50_frac": "Half-rain point (fraction)", "dry_days_before": "Dry days before (d)",
                "prev7_mm": "Prior 7 days (mm)", "prev30_mm": "Prior 30 days (mm)"}
    pd.DataFrame(rows, columns=cols).rename(columns=readable).to_csv(lab / "storm_table_readable.csv", index=False)
    (lab / "score_run_sycamore.py").write_text(
        "def arrival(q, level):  # first time flow RISES THROUGH level; returns (time, flag)\n    return None, 'never'\n\n"
        "def main():\n    pg = 1.0\n    level = 0.2 * pg   # same level for gauge and model\n    return arrival(None, level)\n")
    shutil.copy(SCRIPT, lab / "score_timing_sycamore.py")
    return lab, info


def run(lab, *args):
    p = subprocess.run([sys.executable, str(lab / "score_timing_sycamore.py"), *args], cwd=lab, capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def csv(path):
    return pd.read_csv(path, comment="#", dtype=str, keep_default_na=False)


def tree_hash(root, skip=("timing_out",)):
    h = {}
    for p in sorted(Path(root).rglob("*")):
        if p.is_file() and not any(s in p.parts for s in skip):
            h[str(p.relative_to(root))] = hashlib.md5(p.read_bytes()).hexdigest()
    return h


def fl(x):
    return float(x) if x not in ("", None) else float("nan")


# ============================================================================================ integration tests
def integration_tests():
    print("\nINTEGRATION TESTS (synthetic lab; volumes equal the volumes on record)")
    tmp = Path(tempfile.mkdtemp(prefix="stt_"))
    try:
        lab, info = build_lab(tmp)
        before = tree_hash(lab.parent)
        # --- dry run
        rc, out = run(lab, "--pass", "1", "--dry-run")
        check(rc == 0 and "DRY RUN" in out and out.count("last rain day + 4 d") + out.count("cut at the next storm") == 12, "dry run: 12 windows listed, exit 0")
        check(not (lab / "timing_out").exists(), "dry run writes nothing")
        # --- pass 1
        rc, out = run(lab, "--pass", "1")
        check(rc == 0, "pass 1 exits 0" + ("" if rc == 0 else "\n" + out[-1500:]))
        check("VOLUME CHECK: MATCH" in out, "pass 1: the volume check reproduces the record (MATCH)")
        p1 = lab / "timing_out" / "pass1_gauge_table.csv"
        check(p1.is_file() and p1.read_text().startswith("# score_timing_sycamore.py"), "pass 1 CSV exists and starts with the definitions header")
        d = csv(p1)
        check(len(d) == 12, "pass 1: 12 rows (one per storm, both years)")
        check(not any(c.startswith("model") for c in d.columns), "pass 1 has no model columns")
        row = {(r.water_year, r.storm): r for r in d.itertuples()}
        r = row[("WY2017", "Dec 21-25")]
        check(r.gauge_t50_h == "" and "volume is zero" in r.blank_reasons, "Dec 21-25 (gauge volume 0.00): times blank with the reason, not zero")
        r = row[("WY2017", "Dec 30-Jan 2")]
        check(r.jan1_2_outage_flag != "" and r.gauge_hours_with_data == "144" and "Jan 01 00:00-Jan 02 00:00 (24 h)" in r.gauge_missing_intervals,
              "Dec 30-Jan 2: outage flag set, 144 gauge hours, the 24 h gap listed")
        r = row[("WY2017", "Dec 21-25")]
        check(r.jan1_2_outage_flag == "", "Dec 21-25 does not touch the outage")
        r = row[("WY2020", "Dec 23-28")]
        check("Dec 27 00:00-Dec 30 00:00 (72 h)" in r.gauge_missing_intervals and r.gap_before_t90 == "no", "Dec 23-28: the Dec 27 dropout is listed, after t90")
        r = row[("WY2020", "Nov 28-30")]
        check(r.window_end == "2019-12-04 00:00" and "Dec 02 00:00-Dec 04 00:00 (48 h)" in r.gauge_missing_intervals, "Nov 28-30: uncut window to Dec 4, with the Dec 2 dropout listed")
        r = row[("WY2017", "Jan 19-21")]
        check(r.window_end == "2017-01-23 00:00" and r.window_hours == "96", "Jan 19-21: cut at Jan 23 (96 h)")
        r = row[("WY2017", "Feb 27-28")]
        check(abs(fl(r.pct_interp_or_est) - 100 * 30 / 480) < 1e-3, f"Feb 27-28: 30 of 480 readings interpolated/estimated = 6.25% (got {r.pct_interp_or_est})")
        # gauge peak and time against the construction
        okp = True
        for (y, n, f_, l_, nx) in st.STORMS:
            rr = row[(y, n)]
            fac = info["gauge_factor"][(y, n)]
            if fac == 0:
                continue
            lo = pd.Timestamp(f_)
            exp_pk = fac * st.CFS_PER_CMS
            okp &= abs(fl(rr.gauge_peak_cfs) - exp_pk) <= 1e-3 * max(1, exp_pk) and rr.gauge_peak_time == f"{lo + pd.Timedelta(hours=PULSE_OFFSET + RISE):%Y-%m-%d %H:%M}"
        check(okp, "gauge peak (cfs) and its time equal the construction for every storm")
        check(row[("WY2020", "Mar 18-19")].rain_basin_mm == "24.61" and row[("WY2017", "Dec 21-25")].rain_prev30_mm == "41.998", "rain columns appended from storm_table.csv")
        check(fl(row[("WY2020", "Feb 10-11")].gauge_left_edge_cfs) == 0.0, "left edge flow reported (0 before the pulse)")
        # --- the lab's readable headers: same rain values, no blank-columns note
        rc, out = run(lab, "--pass", "1", "--table", "storm_table_readable.csv", "--out-dir", "timing_out/readable")
        check(rc == 0 and "those rain columns are blank" not in out, "readable-header storm table: pass 1 exits 0 and gives no blank-rain-columns note" + ("" if rc == 0 else "\n" + out[-1500:]))
        dr = csv(lab / "timing_out" / "readable" / "pass1_gauge_table.csv")
        rowr = {(r.water_year, r.storm): r for r in dr.itertuples()}
        check(rowr[("WY2020", "Mar 18-19")].rain_basin_mm == "24.61" and rowr[("WY2017", "Dec 21-25")].rain_prev30_mm == "41.998"
              and rowr[("WY2017", "Dec 21-25")].rain_dry_days_before == row[("WY2017", "Dec 21-25")].rain_dry_days_before,
              "readable-header storm table: the nine rain columns are filled with the same values as the short-header table")
        # a table with a rain column missing: that one column blank, the rest filled, and the note names it
        tb = pd.read_csv(lab / "storm_table_readable.csv", dtype=str, keep_default_na=False).drop(columns=["Prior 30 days (mm)"])
        tb.to_csv(lab / "timing_out" / "storm_table_missing30.csv", index=False)
        rc, out = run(lab, "--pass", "1", "--table", "timing_out/storm_table_missing30.csv", "--out-dir", "timing_out/missing30")
        dm = csv(lab / "timing_out" / "missing30" / "pass1_gauge_table.csv")
        rowm = {(r.water_year, r.storm): r for r in dm.itertuples()}
        check(rc == 0 and "prev30_mm" in out and rowm[("WY2020", "Mar 18-19")].rain_prev30_mm == "" and rowm[("WY2020", "Mar 18-19")].rain_basin_mm == "24.61",
              "a table missing one rain column: that column blank, the others filled, the note names it")
        # --- pass 2 needs a floor
        rc, out = run(lab, "--pass", "2")
        check(rc != 0 and "--floor-mm3" in out, "pass 2 without --floor-mm3 stops and says why")
        # --- pass 2
        rc, out = run(lab, "--pass", "2", "--floor-mm3", "0.05")
        check(rc == 0, "pass 2 exits 0" + ("" if rc == 0 else "\n" + out[-2000:]))
        p2 = lab / "timing_out" / "pass2_timing_long.csv"
        d2 = csv(p2)
        check(len(d2) == 4 * 5 + 12 * 7, f"pass 2: {4 * 5 + 12 * 7} rows (4 runs x 5 storms + 12 runs x 7 storms), got {len(d2)}")
        r2 = {(r.water_year, r.run_key, r.storm): r for r in d2.itertuples()}
        def get(y, k, s):
            return r2[(y, k, s)]
        good = [("WY2017", "retuned", "Dec 21-25"), ]
        ok_exact, ok_grid, n_exact, n_grid = True, True, 0, 0
        fine = {("WY2020", "retuned"), ("WY2017", "original")}          # these two files are written every 3.75 min, the gauge every 15 min
        for (y, n, f_, l_, nx) in st.STORMS:
            if (y, n) in (("WY2017", "Dec 21-25"), ("WY2017", "Dec 30-Jan 2")):
                continue
            for k, sh in (("retuned", 2), ("original", 1), ("kchan05", 1), ("kchan02", 1)) + ((("S1", 3), ("S11", 3)) if y == "WY2020" else ()):
                if (y, k) not in {(yy, kk) for (yy, kk, _) in r2}:
                    continue
                r = get(y, k, n)
                if r.t50_model_minus_gauge_h == "":      # blanked by the floor
                    continue
                rmin = 0.999 if (y, k) in fine else 0.99999
                good = (abs(fl(r.peak_time_model_minus_gauge_h) - sh) < 1e-6 and int(r.xcorr_best_lag_h) == sh and fl(r.xcorr_best_r) > rmin)
                tt = [fl(r.t10_model_minus_gauge_h), fl(r.t50_model_minus_gauge_h), fl(r.t90_model_minus_gauge_h)]
                if (y, k) in fine:
                    n_grid += 1
                    ok_grid &= good and all(sh - 0.12 < x < sh for x in tt)       # sampling offset: model looks 0.09 h early
                else:
                    n_exact += 1
                    ok_exact &= good and all(abs(x - sh) < 1e-6 for x in tt)
        check(ok_exact and n_exact > 30, f"model = gauge shape shifted by 1, 2 or 3 h on the same 15-min grid: t10, t50, t90, peak time and lag equal the shift exactly ({n_exact} rows)")
        check(ok_grid and n_grid > 8, f"same, but the model written every 3.75 min: peak time and lag exact, volume times 0.09 h early from the sampling convention ({n_grid} rows)")
        r = get("WY2020", "retuned", "Feb 22-23")
        check(0.5 < fl(r.pearson_r_hourly) < 0.995, f"Pearson r at lag 0 is below the best-lag correlation for a shifted model (r {r.pearson_r_hourly}, best {r.xcorr_best_r})")
        check(abs(fl(r.tail_gauge_Mm3) - 0.1 * fl(r.gauge_vol_Mm3)) < 1e-3, "gauge tail volume is 10% of the gauge volume")
        check(fl(r.model_tail_fraction) > 0 and fl(r.tail_volume_ratio) > 0, "tail volume ratio and model tail fraction are filled for a normal storm")
        # floor behavior
        r = get("WY2020", "retuned", "Nov 19-22")
        check(r.below_floor == "gauge" and r.t50_model_minus_gauge_h == "" and "below the floor" in r.blank_reasons and r.peak_ratio_model_over_gauge != "",
              "Nov 19-22 (gauge 0.015 < floor 0.05): timing blank with the floor reason, peak ratio still shown")
        r = get("WY2020", "kchan05", "Dec 23-28")
        check(r.below_floor == "" and r.t50_model_minus_gauge_h != "" and fl(r.model_vol_Mm3) < 0.05 < fl(r.gauge_vol_Mm3),
              "Kchan 0.5 Dec 23-28 (model 0.044, under the floor 0.05; gauge above it): the floor is on the gauge only, so timing is still given")
        check("% of gauge" in out, "the printed table notes a run whose own volume is a small fraction of the gauge's")
        check(all(not (r.below_floor and "model" in r.below_floor) for r in d2.itertuples()), "with --floor-mm3 alone, no row is blanked because of the model's volume")
        # --model-floor-mm3: the model side as well (option A)
        rc, out2 = run(lab, "--pass", "2", "--floor-mm3", "0.05", "--model-floor-mm3", "0.05", "--out-dir", "timing_out/modelfloor")
        check(rc == 0, "pass 2 with --model-floor-mm3 exits 0" + ("" if rc == 0 else "\n" + out2[-1500:]))
        dm2 = csv(lab / "timing_out" / "modelfloor" / "pass2_timing_long.csv")
        rm2 = {(r.water_year, r.run_key, r.storm): r for r in dm2.itertuples()}
        r = rm2[("WY2020", "kchan05", "Dec 23-28")]
        check(r.below_floor == "model" and r.t50_model_minus_gauge_h == "" and "model floor" in r.blank_reasons, "with --model-floor-mm3 0.05: Kchan 0.5 Dec 23-28 timing blank, below_floor says model")
        r = rm2[("WY2020", "retuned", "Nov 19-22")]
        check(r.below_floor in ("gauge", "gauge and model") and r.t50_model_minus_gauge_h == "", "with --model-floor-mm3: a gauge-small storm is still blank")
        check("MODEL FLOOR 0.05" in (lab / "timing_out" / "modelfloor" / "pass2_timing_long.csv").read_text().split("\n\n")[0] or "MODEL FLOOR 0.05" in (lab / "timing_out" / "modelfloor" / "pass2_timing_long.csv").read_text()[:20000],
              "the CSV definitions header states the model floor")
        check("No model floor was set" in (lab / "timing_out" / "pass2_timing_long.csv").read_text()[:20000], "the default CSV header says no model floor was set")
        rc, out3 = run(lab, "--pass", "2", "--model-floor-mm3", "0.05", "--out-dir", "timing_out/x")
        check(rc != 0 and "--floor-mm3" in out3, "--model-floor-mm3 without --floor-mm3 stops and says why")
        r = get("WY2017", "retuned", "Dec 21-25")
        check(r.t50_model_minus_gauge_h == "" and r.blank_reasons != "" and fl(r.vol_ratio_model_over_gauge) != fl(r.vol_ratio_model_over_gauge) or r.vol_ratio_model_over_gauge == "",
              "Dec 21-25 (gauge volume 0): ratio and timing blank with reasons")
        r = get("WY2017", "retuned", "Dec 30-Jan 2")
        check(r.jan1_2_outage_flag != "" and r.volume_check == "MATCH", "Dec 30-Jan 2 row carries the outage flag and its volume check")
        check(set(d2.role) == {"calibration", "validation"} and set(d2[d2.water_year == "WY2017"].role) == {"calibration"}, "role column: WY2017 calibration, WY2020 validation")
        check("CHECK 5" in out and "-0.58" in out and "volume t50" in out, "CHECK 5 is printed with the old -0.58 h figure")
        check(all(d2.volume_check.isin(["MATCH", "NOT CHECKED"])), "no row says MISMATCH when the record is reproduced")
        check(p2.read_text().count("\n# ") >= 15 and "TIME DIFFERENCES are always model minus gauge" in p2.read_text(), "pass 2 CSV carries the definitions header")
        # --- read-only
        check(tree_hash(lab.parent) == before, "no file outside timing_out/ changed (md5 of every input and script)")
        # --- re-run overwrites its own file; refuses a foreign file
        rc, out = run(lab, "--pass", "1")
        check(rc == 0, "re-running pass 1 overwrites its own output")
        (lab / "timing_out" / "pass1_gauge_table.csv").write_text("a,b\n1,2\n")
        rc, out = run(lab, "--pass", "1")
        check(rc != 0 and "not written by this script" in out, "refuses to overwrite a file it did not write")
        shutil.rmtree(lab / "timing_out")
        # --- windows modes
        rc, out = run(lab, "--pass", "1", "--windows", "uncut")
        check(rc == 0 and "VOLUME CHECK: MATCH" in out, "--windows uncut: Jan 19-21 to Jan 25 reproduces the 7.1497 / 7.6066 storm-window volumes")
        rc, out = run(lab, "--pass", "1", "--windows", "cut")
        d = csv(lab / "timing_out" / "pass1_gauge_table.csv")
        check(rc == 0 and d[(d.storm == "Nov 28-30")].window_end.iloc[0] == "2019-12-03 00:00" and "not checked (the WY2020 windows are not the scorer's)" in out,
              "--windows cut: Nov 28-30 ends Dec 3, and the seven-storm totals are skipped")
        shutil.rmtree(lab / "timing_out")
        # --- show arrival
        rc, out = run(lab, "--show-arrival")
        check(rc == 0 and "level = 0.2 * pg" in out and "md5" in out, "--show-arrival prints the arrival function and the level line from score_run_sycamore.py")
        # --- missing files
        lab2 = Path(tmp) / "lab2"
        shutil.copytree(lab.parent, lab2)
        res2 = lab2 / "lab" / "results"
        (res2 / "Sycamore_WY2020_Original_kchan0p2_Outlet.qout").unlink()
        (res2 / "W5_d012_k10_base_kchan0p5_Outlet.qout").unlink()
        rc, out = run(lab2 / "lab", "--pass", "2", "--floor-mm3", "0.05")
        d2b = csv(lab2 / "lab" / "timing_out" / "pass2_timing_long.csv")
        check(rc == 0 and "NOT FOUND" in out and "Original + Kchan 0.2" in out and len(d2b) == 3 * 5 + 11 * 7, "two outlet files missing: they are named, the rest are scored")
        (res2 / "W5_d012_k10_base_nocl_Outlet.qout").unlink()
        (res2 / "W5_d012_k10_f0p15_as300_ksx1p03_nocl_Outlet.qout").unlink()
        rc, out = run(lab2 / "lab", "--pass", "1")
        check(rc != 0 and "priority runs" in out, "neither Re-tuned nor Original found for a year: stops")
        shutil.copy(lab / "results" / "W5_d012_k10_base_nocl_Outlet.qout", res2 / "my_original_Outlet.qout")
        shutil.copy(lab / "results" / "W5_d012_k10_f0p15_as300_ksx1p03_nocl_Outlet.qout", res2 / "my_retuned_Outlet.qout")
        rc, out = run(lab2 / "lab", "--pass", "1", "--stem", "2017:original=my_original", "--stem", "2017:retuned=my_retuned_Outlet.qout")
        check(rc == 0 and "my_original_Outlet.qout" in out and "my_retuned_Outlet.qout" in out, "--stem names another outlet file for a run (with or without _Outlet.qout)")
        # --- volume mismatch stops
        lab3, _ = build_lab(Path(tmp) / "m", jitter={("WY2020", "original"): 1.05})
        b3 = tree_hash(lab3.parent)
        rc, out = run(lab3, "--pass", "1")
        check(rc != 0 and "MISMATCH" in out and "STOPPED" in out and not (lab3 / "timing_out").exists(), "a 5% change in one outlet file: the volume check MISMATCHes, the script stops, nothing is written")
        rc, out = run(lab3, "--pass", "2", "--floor-mm3", "0.05", "--ignore-check-failure")
        d3 = csv(lab3 / "timing_out" / "pass2_timing_long.csv")
        check(rc == 0 and "MISMATCH" in set(d3.volume_check) and "MISMATCH, carried on" in (lab3 / "timing_out" / "pass2_timing_long.csv").read_text(),
              "--ignore-check-failure carries on, and the rows say MISMATCH")
        # --- no tier column
        lab4, _ = build_lab(Path(tmp) / "n", with_tier=False)
        rc, out = run(lab4, "--pass", "1")
        d4 = csv(lab4 / "timing_out" / "pass1_gauge_table.csv")
        check(rc == 0 and set(d4.pct_interp_or_est) == {""} and "no 'tier' column" in d4.blank_reasons.iloc[0], "gauge file without a tier column: percentages blank with a reason, no crash")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def gap_before_test():
    print("\nGAP BEFORE T90 (a gap inside the active part of the storm)")
    lo = pd.Timestamp("2020-01-01")
    idx = pd.date_range(lo, lo + 96 * HOUR, freq="15min", inclusive="left")
    q = pd.Series(pulse(idx, lo + 10 * HOUR, 6, 24).to_numpy(), index=idx)
    q[(idx >= lo + 14 * HOUR) & (idx < lo + 20 * HOUR)] = np.nan
    w = st.Win("WY2020", "test", "2020-01-01", "2020-01-01", None, "uncut")
    w.hi = lo + 96 * HOUR
    w.hours = 96
    G = st.gauge_bundle(w, pd.DataFrame({"q_cms": q.to_numpy()}, index=idx))
    check(G["gap_before_t90"] == "yes" and len(G["runs"]) == 1, "a 6 h gap before the 90% time sets gap_before_t90 = yes")
    check(any("tier" in x for x in G["why"]), "no tier column: reason recorded")


def short_model_test():
    print("\nSHORT MODEL FILE (the run ends inside the window)")
    lo = pd.Timestamp("2020-01-01")
    idx = pd.date_range(lo, lo + 96 * HOUR, freq="15min", inclusive="left")
    q = pulse(idx, lo + 10 * HOUR, 6, 24)
    w = st.Win("WY2020", "test", "2020-01-01", "2020-01-01", None, "uncut")
    w.hi, w.hours = lo + 96 * HOUR, 96
    G = st.gauge_bundle(w, pd.DataFrame({"q_cms": q.to_numpy()}, index=idx))
    row = st.pass2_row(w, G, "x", "X", "stem", q[q.index < lo + 60 * HOUR], 0.0, {}, {})
    check(row["t50_model_minus_gauge_h"] is None and "only 60 of the window's 96" in row["blank_reasons"], "a model that ends 36 h early: every model value blank, reason names the hours")
    row = st.pass2_row(w, G, "x", "X", "stem", q * 0.5, 0.0, {}, {})
    check(row["t50_model_minus_gauge_h"] == 0.0 and abs(row["vol_ratio_model_over_gauge"] - 0.5) < 1e-9 and row["xcorr_best_lag_h"] == 0,
          "half-size copy of the gauge: timing differences 0, volume ratio 0.5, lag 0")
    gt = G["tot"]
    row = st.pass2_row(w, G, "x", "X", "stem", q * 0.05, 0.5 * gt, {}, {})
    check(not row["below_floor"] and row["t50_model_minus_gauge_h"] == 0.0 and abs(row["vol_ratio_model_over_gauge"] - 0.05) < 1e-9,
          "gauge floor only: a model with 5% of the gauge volume still gets timing (and the volume ratio shows 0.05)")
    row = st.pass2_row(w, G, "x", "X", "stem", q * 0.05, 0.5 * gt, {}, {}, model_floor=0.5 * gt)
    check(row["below_floor"] == "model" and row["t50_model_minus_gauge_h"] is None and "model floor" in row["blank_reasons"],
          "model floor given: the same model's timing is blank and below_floor says model")
    row = st.pass2_row(w, G, "x", "X", "stem", q, 2.0 * gt, {}, {})
    check(row["below_floor"] == "gauge" and row["t50_model_minus_gauge_h"] is None and row["peak_ratio_model_over_gauge"] == 1.0,
          "gauge under the floor: timing blank, peak ratio still given")


def main():
    mut = "--mutations" in sys.argv
    if mut:
        return mutations()
    unit_tests()
    gap_before_test()
    short_model_test()
    if os.environ.get("ST_SKIP_INTEGRATION") != "1":
        integration_tests()
    print("\n" + ("ALL TESTS PASSED" if not FAILS else f"{len(FAILS)} TEST(S) FAILED:\n  - " + "\n  - ".join(FAILS)))
    return 0 if not FAILS else 1


MUTATIONS = [
    ("lag sign flipped", "m2.index = m2.index - pd.Timedelta(hours=L)", "m2.index = m2.index + pd.Timedelta(hours=L)", True),
    ("difference sign flipped", 'row[f"{c}_model_minus_gauge_h"] = r4(mt[f] - gt[f])', 'row[f"{c}_model_minus_gauge_h"] = r4(gt[f] - mt[f])', False),
    ("window ends 3 d after the last rain day", "self.hi_uncut = self.last + pd.Timedelta(days=4)", "self.hi_uncut = self.last + pd.Timedelta(days=3)", True),
    ("volume uses 3000 s per hour", "return float(h.sum() * 3600.0 / 1e6) if len(h) else float(\"nan\")", "return float(h.sum() * 3000.0 / 1e6) if len(h) else float(\"nan\")", False),
    ("tail starts at t50", "t0 = w.lo + pd.Timedelta(hours=gt[0.9])", "t0 = w.lo + pd.Timedelta(hours=gt[0.5])", False),
    ("no interpolation inside the hour", "out[f] = float(t0[k] + (target - prev) / v[k])", "out[f] = float(t0[k])", True),
    ("constant series gives r = 0 instead of blank", 'return None, len(idx), "gauge hourly flow is constant in the window"', "return 0.0, len(idx), \"\"", True),
    ("second bump ignores the dip", "if float(v[a:b + 1].min()) <= DIP_FRAC * v[j]:", "if True:", True),
    ("missing hours filled with zero", "return x.resample(\"60min\").mean().dropna()", "return x.resample(\"60min\").mean().fillna(0.0)", True),
    ("model volume gates timing even without --model-floor-mm3", "if model_floor is not None and (not isnum(mtot) or mtot < model_floor):", "if not isnum(mtot) or mtot < (model_floor if model_floor is not None else floor):", True),
    ("percent estimated counts only 'estimated'", 'quality["pct_interp_or_est"] = 100.0 * float(t.isin(["estimated", "filled_interp"]).sum()) / n', 'quality["pct_interp_or_est"] = 100.0 * float(t.isin(["estimated"]).sum()) / n', False),
]


def mutations():
    print("MUTATION CHECK: each broken copy of the script must make this test fail")
    src = SCRIPT.read_text()
    bad = 0
    tmp = Path(tempfile.mkdtemp(prefix="stm_"))
    try:
        for name, old, new, unit_only in MUTATIONS:
            if src.count(old) != 1:
                print(f"  ??    {name}: the text to change was found {src.count(old)} times (expected 1)")
                bad += 1
                continue
            p = tmp / "score_timing_sycamore.py"
            p.write_text(src.replace(old, new))
            env = dict(os.environ, ST_SCRIPT=str(p), ST_SKIP_INTEGRATION="1" if unit_only else "0")
            r = subprocess.run([sys.executable, str(Path(__file__).resolve())], env=env, capture_output=True, text=True)
            caught = r.returncode != 0
            print(("  CAUGHT  " if caught else "  MISSED  ") + name + f"   (exit {r.returncode})")
            bad += 0 if caught else 1
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("\n" + ("ALL MUTATIONS CAUGHT" if not bad else f"{bad} MUTATION(S) NOT CAUGHT"))
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
