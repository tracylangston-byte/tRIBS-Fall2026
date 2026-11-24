#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_score_timing.py   (v1.2, 2026-11-22)
=========================================
Stand-in test of score_timing_sycamore.py v1.4: no tRIBS, no real data. Builds a fake lab folder in a temp directory with
synthetic gauge and outlet files, then checks the formulas against answers worked out independently, the Step 0 windows
and their hour check, the per-run start and length guard, the volume-ratio check, the medians file, the old window rules
(--windows as-scored), and that no existing file is touched.

Run (from the folder that holds the scripts):   python test_score_timing.py
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
WINDOWS_CSV = Path(os.environ.get("ST_WINDOWS", HERE / "step0_windows.csv"))
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


def stops(fn, *args, **kw):
    """Run fn; return (True, message) if it stopped with SystemExit, else (False, '')."""
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            fn(*args, **kw)
    except SystemExit as exc:
        return True, str(exc)
    return False, ""


# ============================================================================================ unit tests: formulas (unchanged from v1.1)
def unit_formula_tests():
    print("\nUNIT TESTS: FORMULAS (definitions are unchanged in v1.4)")
    lo = pd.Timestamp("2020-01-01 00:00")
    v = np.zeros(48)
    v[10:20] = 1.0
    tot, t, why = st.volume_times(hs(v), lo)
    check(abs(tot - 10 * 3600 / 1e6) < 1e-12, "volume of a 10 h rectangle = 0.036 M m3")
    check(t is not None and abs(t[0.1] - 11) < 1e-9 and abs(t[0.5] - 15) < 1e-9 and abs(t[0.9] - 19) < 1e-9,
          f"rectangle: t10, t50, t90 = 11, 15, 19 h (got {None if t is None else [round(t[k], 6) for k in (0.1, 0.5, 0.9)]})")
    rng = np.random.default_rng(7)
    v = np.abs(rng.normal(0, 1, 96)) * np.exp(-((np.arange(96) - 30) / 15.0) ** 2)
    tot, t, why = st.volume_times(hs(v), lo)
    cum = np.concatenate([[0.0], np.cumsum(v)])
    alt = [float(np.interp(f * cum[-1], cum, np.arange(97, dtype=float))) for f in (0.1, 0.5, 0.9)]
    check(t is not None and np.allclose([t[0.1], t[0.5], t[0.9]], alt, atol=1e-9), "random hydrograph: matches inverse interpolation of the cumulative curve")
    s = hs(np.ones(48)).drop(hs(np.ones(48)).index[20:30])
    tot, t, why = st.volume_times(s, lo)
    check(abs(tot - 38 * 3600 / 1e6) < 1e-12, "gap of 10 h: volume counts only the 38 hours with data")
    runs = st.missing_runs(s, lo, lo + 48 * HOUR)
    check(len(runs) == 1 and (runs[0][1] - runs[0][0]) == 10 * HOUR, "missing_runs finds one 10 h gap")
    tot, t, why = st.volume_times(hs(np.zeros(10)), lo)
    check(t is None and "zero" in why, "zero volume: times undefined, reason given")
    tot, t, why = st.volume_times(hs([]), lo)
    check(t is None and "no hourly" in why, "no data: times undefined, reason given")
    v = np.abs(rng.normal(0, 1, 120)) * np.exp(-((np.arange(120) - 25) / 20.0) ** 2) + 0.01
    h = hs(v)
    tot, t, why = st.volume_times(h, lo)
    tail = st.volume_after(h, lo, lo + 120 * HOUR, lo + pd.Timedelta(hours=t[0.9]))
    check(abs(tail - 0.1 * tot) < 1e-9, "tail after t90 is 10% of the volume")
    base = np.exp(-((np.arange(120) - 50) / 6.0) ** 2) + 0.4 * np.exp(-((np.arange(120) - 80) / 4.0) ** 2)
    g = hs(base)
    for shift in (3, -2, 0, 12):
        lag, r, atl, why = st.best_lag(g, hs(np.roll(base, shift)))
        check(lag == shift and r > 0.999, f"best lag of a series shifted by {shift:+d} h is {shift:+d} (got {lag})")
    r, n, why = st.pearson(g, g * 3.0)
    check(abs(r - 1) < 1e-12 and n == 120, "Pearson r of a series with 3 x itself is 1")
    r, n, why = st.pearson(g, hs(np.ones(120)))
    check(r is None and "model" in why and "constant" in why, "Pearson r with a constant MODEL series is blank with a reason (not zero)")
    r, n, why = st.pearson(hs(np.ones(120)), g)
    check(r is None and "gauge" in why and "constant" in why, "Pearson r with a constant GAUGE series is blank with a reason (not zero)")

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
    idx = pd.date_range(lo, periods=48 * 4, freq="15min")
    q = pd.Series(np.exp(-(((idx - lo) / HOUR).to_numpy(float) - 3) ** 2 / 4.0), index=idx)
    pk, why = st.peak_info(q, st.hourly(q, lo, lo + 48 * HOUR), lo, lo + 48 * HOUR)
    check(pk is not None and "edge" in pk["flags"], "peak at hour 3 of the window is flagged as near the edge")
    q = pd.Series(np.exp(-(((idx - lo) / HOUR).to_numpy(float) - 24) ** 2 / 4.0), index=idx)
    pk, why = st.peak_info(q, st.hourly(q, lo, lo + 48 * HOUR), lo, lo + 48 * HOUR)
    check(pk is not None and pk["flags"] == "" and pk["time"] == lo + 24 * HOUR, "peak at hour 24: no flag, right time")
    q = pd.Series(np.zeros(len(idx)), index=idx)
    q[(idx >= lo + 20 * HOUR) & (idx < lo + 23 * HOUR)] = 5.0
    pk, why = st.peak_info(q, st.hourly(q, lo, lo + 48 * HOUR), lo, lo + 48 * HOUR)
    check(pk is not None and "flat top" in pk["flags"], "a 3 h plateau is flagged as flat-topped")
    idx = pd.date_range(lo, periods=40, freq="15min")
    ramp = pd.Series(np.maximum(0, (np.arange(40) - 8) * 0.5), index=idx)
    t, why = st.arrival(ramp, 2.0, lo, lo + 10 * HOUR)
    check(t is not None and abs(t - 3.0) < 1e-9, f"arrival through 2.0 on a ramp: 3.0 h (got {t})")
    t, why = st.arrival(ramp, 1e9, lo, lo + 10 * HOUR)
    check(t is None and "never" in why, "arrival never reached: blank with reason")
    t, why = st.arrival(ramp + 5, 2.0, lo, lo + 10 * HOUR)
    check(t is None and "already" in why, "arrival already above at the start: blank with reason")
    t, why = st.arrival(ramp.drop(ramp.index[10:30]), 3.0, lo, lo + 10 * HOUR)
    check(t is None and "gap" in why, "first rise inside a reading gap: blank with reason")
    check("24 h" in st.fmt_runs([(lo, lo + 24 * HOUR)]), "fmt_runs prints the length of the gap")


# ============================================================================================ unit tests: what is new in v1.4
def unit_v14_tests():
    print("\nUNIT TESTS: v1.4 (windows, runs, ratio source, medians, length guard)")
    tmp = Path(tempfile.mkdtemp(prefix="stu_"))
    try:
        print("  -- Step 0 windows")
        w0 = st.load_step0_windows(WINDOWS_CSV)
        req = [188, 153, 85, 135, 121, 160, 141, 210, 132, 121, 187, 120]
        check([w.hours for w in w0] == req, "step0_windows.csv: the 12 window lengths equal 188 153 85 135 121 160 141 210 132 121 187 120")
        check([w.gauge_response for w in w0] == ["none", "yes", "yes", "yes", "yes", "none", "yes", "yes", "yes", "yes", "yes", "yes"],
              "gauge_response = none for Dec 21-25 2016 and Nov 19-22 2019, yes for the other ten")
        check([w.no_gauge_response for w in w0].count("yes") == 2 and w0[0].no_gauge_response == "yes" and w0[5].no_gauge_response == "yes",
              "no_gauge_response is yes for exactly those two storms")
        b = {w.key: w for w in w0}
        check(b[("WY2017", "Dec 21-25")].lo == pd.Timestamp("2016-12-21 14:00") and b[("WY2020", "Dec 23-28")].hi == pd.Timestamp("2020-01-01 15:00"),
              "window edges are read as written (Dec 21 14:00; Dec 23-28 ends Jan 1 2020 15:00)")
        base = pd.read_csv(WINDOWS_CSV, comment="#", dtype=str, keep_default_na=False)

        def variant(name, mutate):
            d = base.copy()
            mutate(d)
            p = tmp / f"{name}.csv"
            d.to_csv(p, index=False)
            return p
        i30 = base.index[base.storm == "Dec 30-Jan 2"][0]

        def one_hour_longer(d):
            d.loc[i30, "end"] = "2017-01-06 04:00"
        s, msg = stops(st.load_step0_windows, variant("long", one_hour_longer))
        check(s and "154 h" in msg and "153 h is required" in msg, "a window one hour too long: the script stops and names the storm and both lengths")
        s, msg = stops(st.load_step0_windows, variant("drop", lambda d: d.drop(index=[3], inplace=True)))
        check(s and "missing windows" in msg and "Feb 18-19" in msg, "a missing storm: stops and names it")

        def dup(d):
            d.loc[len(d)] = d.loc[0]
        s, msg = stops(st.load_step0_windows, variant("dup", dup))
        check(s and "listed twice" in msg, "a storm listed twice: stops")

        def odd_minute(d):
            d.loc[0, "start"] = "2016-12-21 14:30"
        s, msg = stops(st.load_step0_windows, variant("min", odd_minute))
        check(s and "whole hours" in msg, "a start that is not on the hour: stops")

        def bad_resp(d):
            d.loc[2, "gauge_response"] = "maybe"
        s, msg = stops(st.load_step0_windows, variant("resp", bad_resp))
        check(s and "yes or none" in msg, "gauge_response other than yes or none: stops")

        def odd_name(d):
            d.loc[2, "storm"] = "Jan 20-22"
        s, msg = stops(st.load_step0_windows, variant("name", odd_name))
        check(s and "not one of the 12" in msg, "an unknown storm name: stops")
        s, msg = stops(st.load_step0_windows, tmp / "nothing.csv")
        check(s and "cannot find the windows table" in msg, "no windows file: stops with a plain message")
        print("  -- old window rules (--windows as-scored|cut|uncut) are as in v1.3")
        wa = {w.key: w for w in st.make_windows("as-scored")}
        wc = {w.key: w for w in st.make_windows("cut")}
        wu = {w.key: w for w in st.make_windows("uncut")}
        exp = {("WY2017", "Dec 21-25"): 192, ("WY2017", "Dec 30-Jan 2"): 168, ("WY2017", "Jan 19-21"): 96, ("WY2017", "Feb 18-19"): 120,
               ("WY2017", "Feb 27-28"): 120, ("WY2020", "Nov 19-22"): 168, ("WY2020", "Nov 28-30"): 144, ("WY2020", "Dec 23-28"): 216,
               ("WY2020", "Feb 10-11"): 120, ("WY2020", "Feb 22-23"): 120, ("WY2020", "Mar 10-13"): 168, ("WY2020", "Mar 18-19"): 120}
        check(all(wa[k].hours == h for k, h in exp.items()), "as-scored window lengths (h) equal the v1.3 lengths")
        check(wa[("WY2020", "Nov 28-30")].hi == pd.Timestamp("2019-12-04") and wa[("WY2017", "Jan 19-21")].hi == pd.Timestamp("2017-01-23"),
              "as-scored: WY2020 Nov 28-30 is not cut (Dec 4); WY2017 Jan 19-21 is cut at Jan 23")
        check([k for k in wa if wa[k].hi != wc[k].hi] == [("WY2020", "Nov 28-30")] and [k for k in wa if wa[k].hi != wu[k].hi] == [("WY2017", "Jan 19-21")],
              "the three old rules differ in exactly one storm each")
        check(all(wa[k].no_gauge_response == ("yes" if k in st.NO_RESPONSE else "no") for k in wa), "old rules: no_gauge_response set for the same two storms")
        print("  -- runs: per-run start and expected hours")
        runs = st.build_runs()
        check(len(runs) == 20 and sum(r.year == "WY2017" for r in runs) == 10 and sum(r.year == "WY2020" for r in runs) == 10,
              "default run list: Re-tuned, Original and 8 screens in both years (20 runs, no Kchan)")
        check(all(r.start == pd.Timestamp("2016-09-01") and r.hours == 6551 for r in runs if r.year == "WY2017")
              and all(r.start == pd.Timestamp("2019-09-01") and r.hours == 6575 for r in runs if r.year == "WY2020"),
              "all default runs start 1 September (WY2017 2016-09-01, 6,551 h; WY2020 2019-09-01, 6,575 h)")
        check([r.label for r in runs if r.year == "WY2020"][2:] == [f"Screen {n} ({k})" for n, k in [(1, "S1"), (2, "S5"), (3, "S6"), (4, "S7"), (5, "S8"), (6, "S9"), (7, "S10"), (8, "S11")]],
              "screen labels: Screen 1 (S1), 2 (S5), 3 (S6), 4 (S7), 5 (S8), 6 (S9), 7 (S10), 8 (S11)")
        check(not any("kchan" in r.key for r in runs), "the four Kchan runs are left out of the default list")
        rk = st.build_runs(include_kchan=True)
        k17 = [r for r in rk if r.year == "WY2017" and "kchan" in r.key]
        k20 = [r for r in rk if r.year == "WY2020" and "kchan" in r.key]
        check(len(rk) == 24 and len(k17) == 2 and len(k20) == 2, "--include-kchan adds the four Kchan runs (24 runs)")
        check(all(r.start == pd.Timestamp("2016-10-01") and r.hours == 8759 for r in k17) and all(r.start == pd.Timestamp("2019-09-01") and r.hours == 6575 for r in k20),
              "the two WY2017 Kchan runs keep their own 1 October start (8,759 h), mixed with 1 September runs: start is set per run, not per year")
        ro = st.build_runs("oct1", True)
        check([r.key for r in ro if r.year == "WY2017"] == ["retuned", "original", "kchan05", "kchan02"] and all(r.start == pd.Timestamp("2016-10-01") for r in ro if r.year == "WY2017"),
              "--wy2017-runs oct1 --include-kchan: the four old W5 runs, all from 2016-10-01")
        ro2 = st.build_runs("oct1", False)
        check([r.key for r in ro2 if r.year == "WY2017"] == ["retuned", "original"], "--wy2017-runs oct1 without Kchan: Re-tuned and Original only")
        ov = st.build_runs(overrides=["2017:original=my_orig_Outlet.qout"])
        o = [r for r in ov if r.year == "WY2017" and r.key == "original"][0]
        check(o.stem == "my_orig" and o.start == pd.Timestamp("2016-09-01") and o.hours == 6551, "--stem renames the file, keeps the run's start and length")
        s, msg = stops(st.build_runs, overrides=["2017:nosuch=x"])
        check(s and "no run 'nosuch'" in msg, "--stem for a run that is not in the list: stops")
        s, msg = stops(st.build_runs, overrides=["2017 original x"])
        check(s and "not in the form" in msg, "--stem in the wrong form: stops")
        print("  -- run design")
        every = {(y, k) for y in ("WY2017", "WY2020") for k in st.ALL_KEYS}
        check(set(st.RUN_DESIGN) == every and len(every) == 24, "RUN_DESIGN covers exactly the 24 run keys (12 per year)")
        check(st.run_design("WY2020", "retuned") == ("validation", "yes") and st.run_design("WY2020", "kchan05") == ("validation", "no")
              and st.run_design("WY2020", "S9") == ("post_validation_screen", "no") and st.run_design("WY2017", "S9") == ("post_validation_screen", "no")
              and st.run_design("WY2017", "retuned") == ("calibration", "n/a"), "roles: WY2020 Re-tuned validation/yes; screens (both years) post_validation_screen/no; WY2017 Re-tuned calibration/n/a")
        print("  -- run file length guard")
        r1 = st.Run("WY2017", "retuned", "Re-tuned", "x", st.SEP1_2016, 6551)
        good = pd.DataFrame({"Time_hr": np.arange(0, 6551 + 0.0625, 0.0625), "Q": 0.0})
        (tmp / "good.qout").write_text("x\n")
        inf = st.run_file_info(r1, tmp / "good.qout", good)
        check(inf["ok"] and inf["last"] == 6551.0, "a file that ends at 6,551 h matches a 1 September WY2017 run")
        old = pd.DataFrame({"Time_hr": np.arange(0, 8759 + 0.25, 0.25), "Q": 0.0})
        inf = st.run_file_info(r1, tmp / "good.qout", old)
        check(not inf["ok"] and "8,759 h is the old 1 October whole-year run" in inf["why"], "an 8,759 h file given to a 1 September run: length does not match, and the message says why")
        r2 = st.Run("WY2017", "kchan05", "K", "x", st.OCT1_2016, 8759)
        check(st.run_file_info(r2, tmp / "good.qout", old)["ok"], "the same 8,759 h file matches the old 1 October Kchan run")
        s20 = pd.DataFrame({"Time_hr": np.arange(0, 6575 + 0.0625, 0.0625), "Q": 0.0})
        check(st.run_file_info(st.Run("WY2020", "retuned", "R", "x", st.SEP1_2019, 6575), tmp / "good.qout", s20)["ok"], "a 6,575 h file matches a WY2020 run")
        check(not st.run_file_info(st.Run("WY2020", "retuned", "R", "x", st.SEP1_2019, 6575), tmp / "good.qout", good)["ok"], "a 6,551 h file given to a WY2020 run: length does not match")
        print("  -- ratio source (window_scores CSV)")
        hdr = "setting,window,kind,name,years,runs,vol_gauge_Mm3,vol_model_Mm3,pbias_pct"
        rows_ = ["Retuned,storm02_2017_Dec30-Jan2,storm,Dec 30-Jan 2,2017,x,0.819931,0.476662,-41.8",
                 "Retuned,storm01_2017_Dec21-25,storm,Dec 21-25,2017,x,0,1.18309,",
                 "Retuned,storm08_2020_Dec23-28,storm,Dec 23-28,2020,x,4.20157,0.224831,-94.6",
                 "Retuned,AllStorms,all_storms,AllStorms,2017+2020,x,31.0983,29.8344,-4.06",
                 "Retuned_old,storm02_2017_Dec30-Jan2,storm,Dec 30-Jan 2,2017,x,0.819931,0.473105,-42.3",
                 "S1,storm02_2017_Dec30-Jan2,storm,Dec 30-Jan 2,2017,x,0.819931,0.0101296,-98.7",
                 "Original,storm08_2020_Dec23-28,storm,Dec 23-28,2020,x,4.20157,0.155322,-96.3"]
        pth = tmp / "ws.csv"
        pth.write_text("# comment line\n" + hdr + "\n" + "\n".join(rows_) + "\n")
        rs = st.read_ratio_source(pth)
        check(abs(rs[("retuned", "WY2017", "Dec 30-Jan 2")] - 0.476662 / 0.819931) < 1e-9 and abs(rs[("original", "WY2020", "Dec 23-28")] - 0.155322 / 4.20157) < 1e-9,
              "ratio source: vol_model_Mm3 / vol_gauge_Mm3 per storm for Retuned and Original")
        check(("retuned", "WY2017", "Dec 21-25") not in rs and len(rs) == 3, "ratio source: zero-gauge storms, AllStorms rows, Retuned_old and screens are not taken")
        check(all(len(v) == 10 for v in st.RATIO_TYPED.values()) and len(st.RATIO_STORMS) == 10 and st.RATIO_TYPED["retuned"][0] == 0.58 and st.RATIO_TYPED["original"][5] == 0.04,
              "typed ratios: 10 per run, in your storm order (Re-tuned 0.58 ... 0.31; Original 0.78 ... 0.25)")
        s, msg = stops(st.read_ratio_source, tmp / "nofile.csv")
        check(s and "cannot read the ratio source" in msg, "no ratio source file: stops")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ============================================================================================ pulses and single-storm tests
def pulse(idx, t0, rise, fall):
    x = ((idx - t0) / HOUR).to_numpy(float)
    y = np.where((x >= 0) & (x < rise), x / rise, np.where((x >= rise) & (x < rise + fall), 1 - (x - rise) / fall, 0.0))
    return pd.Series(y, index=idx)


def one_storm_window(rule="test", resp="yes"):
    lo = pd.Timestamp("2020-01-01")
    return st.Win("WY2020", "test", lo, lo + 96 * HOUR, rule, resp), lo


def row_tests():
    print("\nSINGLE-STORM TESTS (pass2_row: volume ratios, no floor, blanks)")
    w, lo = one_storm_window()
    idx = pd.date_range(lo, lo + 96 * HOUR, freq="15min", inclusive="left")
    q = pulse(idx, lo + 10 * HOUR, 6, 24)
    G = st.gauge_bundle(w, pd.DataFrame({"q_cms": q.to_numpy()}, index=idx))
    check(G["gap_before_t90"] == "no" and any("tier" in x for x in G["why"]), "no tier column: reason recorded; no gap")
    row = st.pass2_row(w, G, "x", "X", "stem", q, None, {})
    check(row["t50_model_minus_gauge_h"] == 0.0 and row["vol_ratio_gauge_hours"] == 1.0 and row["vol_ratio_model_over_gauge"] == 1.0 and row["no_gauge_response"] == "no",
          "model = gauge: timing differences 0, both volume ratios 1.0, no_gauge_response no")
    row = st.pass2_row(w, G, "x", "X", "stem", q * 0.5, None, {})
    check(abs(row["vol_ratio_model_over_gauge"] - 0.5) < 1e-9 and row["xcorr_best_lag_h"] == 0 and row["t50_model_minus_gauge_h"] == 0.0,
          "half-size copy of the gauge: volume ratio 0.5, lag 0, timing differences 0")
    # a gauge with 10 missing hours inside the active part: the two ratios differ, the gauge-hours ratio is 1
    qg = q.copy()
    qg[(idx >= lo + 14 * HOUR) & (idx < lo + 24 * HOUR)] = np.nan
    Gg = st.gauge_bundle(w, pd.DataFrame({"q_cms": qg.to_numpy()}, index=idx))
    row = st.pass2_row(w, Gg, "x", "X", "stem", q, None, {})
    check(abs(row["vol_ratio_gauge_hours"] - 1.0) < 1e-4 and row["vol_ratio_model_over_gauge"] > 1.05 and Gg["gap_before_t90"] == "yes",
          f"gauge with a 10 h gap: vol_ratio_gauge_hours = 1.0 (matched hours), vol_ratio_model_over_gauge > 1 ({row['vol_ratio_model_over_gauge']}); gap before t90 noted")
    check(abs(row["model_vol_gauge_hours_Mm3"] - row["gauge_vol_Mm3"]) < 1e-4, "model_vol_gauge_hours_Mm3 equals the gauge volume when the model equals the gauge on the gauge's hours")
    # no floor: a tiny gauge volume is still scored
    tiny = pd.DataFrame({"q_cms": (q * 1e-4).to_numpy()}, index=idx)
    Gt = st.gauge_bundle(w, tiny)
    row = st.pass2_row(w, Gt, "x", "X", "stem", q * 1e-4 * 3, None, {})
    check(Gt["tot"] < 0.001 and not row["below_floor"] and row["t50_model_minus_gauge_h"] == 0.0 and row["tail_volume_ratio"] is not None,
          "no floor (default): a gauge window under 0.001 M m3 still gets every timing number")
    check(row["blank_reasons"] == "" or "floor" not in row["blank_reasons"], "no floor: no blank reason mentions a floor")
    row = st.pass2_row(w, Gt, "x", "X", "stem", q * 1e-4, 0.5, {})
    check(row["below_floor"] == "gauge" and row["t50_model_minus_gauge_h"] is None, "--floor-mm3 still works when asked for: gauge under it, timing blank")
    wn, _ = one_storm_window(resp="none")
    row = st.pass2_row(wn, G, "x", "X", "stem", q, None, {})
    check(row["no_gauge_response"] == "yes", "a window marked gauge_response none: no_gauge_response = yes on the row")
    # model ends inside the window
    row = st.pass2_row(w, G, "x", "X", "stem", q[q.index < lo + 60 * HOUR], None, {})
    check(row["t50_model_minus_gauge_h"] is None and row["vol_ratio_gauge_hours"] is None and "only 60 of the window's 96" in row["blank_reasons"],
          "a model that ends 36 h early: every model value (including the ratios) blank, reason names the hours")
    # zero gauge volume
    z = st.gauge_bundle(w, pd.DataFrame({"q_cms": np.zeros(len(idx))}, index=idx))
    row = st.pass2_row(w, z, "x", "X", "stem", q, None, {})
    check(row["vol_ratio_gauge_hours"] is None and "gauge volume is zero" in row["blank_reasons"] and row["peak_ratio_model_over_gauge"] is None,
          "zero gauge volume: ratios and peak ratio blank with reasons, not zero")
    check(not any(c.startswith("rain_") for c in st.P2_COLS) and not any(c.startswith("rain_") for c in st.P1_COLS), "no rain_* columns in either output")
    check(all(c in st.P2_COLS for c in ["no_gauge_response", "vol_ratio_gauge_hours", "vol_ratio_model_over_gauge", "model_vol_gauge_hours_Mm3"]), "new columns are in the pass 2 column list")


def medians_tests():
    print("\nMEDIANS FILE TESTS (hand-worked answers)")
    order = list(st.STEP0_HOURS)           # 12 storms, in order
    rows = []
    for i, (yr, nm) in enumerate(order):
        r = {c: None for c in st.P2_COLS}
        r.update({"water_year": yr, "storm": nm, "run_key": "retuned", "run_label": "Re-tuned",
                  "no_gauge_response": "yes" if (yr, nm) in st.NO_RESPONSE else "no"})
        r["peak_ratio_model_over_gauge"] = float(i + 1)          # 1..12 in storm order
        r["peak_time_model_minus_gauge_h"] = -float(i + 1)
        r["t10_model_minus_gauge_h"] = 0.5
        r["t50_model_minus_gauge_h"] = 1.0
        r["t90_model_minus_gauge_h"] = 2.0
        r["xcorr_best_lag_h"] = 0
        r["tail_volume_ratio"] = None if (yr, nm) == st.OUTAGE_STORM else 0.25      # undefined for the outage storm
        r["pearson_r_hourly"] = 0.5
        rows.append(r)
    # a second run: only WY2020 rows
    for r in [dict(x) for x in rows if x["water_year"] == "WY2020"]:
        r.update({"run_key": "S1", "run_label": "Screen 1 (S1)"})
        rows.append(r)
    med = st.build_medians(rows)
    got = {(m["run_key"], m["years"], m["storm_set"], m["measure"]): m for m in med}

    def g(years, sid, measure, key="retuned"):
        return got[(key, years, sid, measure)]
    check(g("both years", "a", "peak ratio")["median"] == 7.5 and g("both years", "a", "peak ratio")["n"] == 10,
          "set (a) both years, peak ratio: median of 2,3,4,5,7,8,9,10,11,12 = 7.5, n = 10")
    check(g("both years", "b", "peak ratio")["median"] == 8.0 and g("both years", "b", "peak ratio")["n"] == 9,
          "set (b) both years: drops the outage storm too, median of 3,4,5,7,8,9,10,11,12 = 8, n = 9")
    check(g("both years", "c", "peak ratio")["median"] == 6.5 and g("both years", "c", "peak ratio")["n"] == 12, "set (c) both years: all 12, median 6.5, n = 12")
    check(g("WY2017", "a", "peak ratio")["median"] == 3.5 and g("WY2017", "b", "peak ratio")["median"] == 4.0 and g("WY2017", "c", "peak ratio")["median"] == 3.0,
          "WY2017 alone: (a) 3.5, (b) 4, (c) 3")
    check(g("WY2020", "a", "peak ratio")["median"] == 9.5 and g("WY2020", "b", "peak ratio")["median"] == 9.5 and g("WY2020", "c", "peak ratio")["median"] == 9.0,
          "WY2020 alone: (a) 9.5, (b) 9.5 (the outage storm is a 2017 storm), (c) 9")
    m = g("both years", "a", "peak ratio")
    check(m["min"] == 2.0 and m["max"] == 12.0 and m["n_storms_in_set"] == 10, "min and max over the same storms (2 and 12); the set holds 10 storms")
    m = g("both years", "a", "tail ratio")
    check(m["n"] == 9 and m["n_storms_in_set"] == 10 and m["median"] == 0.25, "a value undefined for one storm: n = 9 of the 10 in the set, median over the 9 only")
    check(g("both years", "b", "tail ratio")["n"] == 9 and g("both years", "b", "tail ratio")["n_storms_in_set"] == 9, "set (b) takes the outage storm out, so the undefined tail ratio no longer lowers n")
    check(g("both years", "a", "peak time (h)")["median"] == -7.5 and g("both years", "a", "lag (h)")["median"] == 0.0, "peak time (model minus gauge) and lag medians come out with their signs")
    check(("S1", "WY2017", "a", "peak ratio") not in got and g("both years", "a", "peak ratio", "S1")["n_storms_in_set"] == 6 and g("both years", "a", "peak ratio", "S1")["years_scored"] == "WY2020",
          "a run scored for one year only: no WY2017 rows, 'both years' holds its WY2020 storms (6 with a response), years_scored says WY2020")
    keys = {(m["run_key"], m["years"]) for m in med}
    check(keys == {("retuned", "both years"), ("retuned", "WY2017"), ("retuned", "WY2020"), ("S1", "both years"), ("S1", "WY2020")}, "one block per run and years: Re-tuned 3, screen run 2")
    check(len([m for m in med if m["run_key"] == "retuned"]) == 3 * 3 * 8, "Re-tuned: 3 year scopes x 3 storm sets x 8 measures = 72 rows")
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        st.print_medians(med)
    txt = buf.getvalue()
    check("Re-tuned" in txt and "storm sets:" in txt and "(a) storms with a gauge response" in txt, "print_medians prints the Re-tuned rows and the storm-set key")
    print("\nCOMPARE-OLD TESTS")
    tmp = Path(tempfile.mkdtemp(prefix="stc_"))
    try:
        rr = []
        for k, (y, n) in enumerate(order[:4]):
            r = {c: None for c in st.P2_COLS}
            r.update({"water_year": y, "run_key": "retuned", "run_label": "Re-tuned", "storm": n, "peak_ratio_model_over_gauge": 1.5 + k, "blank_reasons": ""})
            rr.append(r)
        old = pd.DataFrame(rr, columns=st.P2_COLS, dtype=object)
        p = tmp / "old.csv"
        with open(p, "w") as fh:
            fh.write("# old file header\n")
            old.to_csv(fh, index=False)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            lines, ok = st.compare_old(rr, p)
        check(ok and any("IDENTICAL" in x for x in lines), "compare_old: identical rows are reported IDENTICAL")
        rr2 = [dict(r) for r in rr]
        rr2[2]["peak_ratio_model_over_gauge"] = 9.9
        lines, ok = st.compare_old(rr2, p)
        check((not ok) and any("cells that differ: 1" in x for x in lines) and any("peak_ratio_model_over_gauge" in x and "9.9" in x for x in lines), "compare_old: one changed cell is found and shown")
        lines, ok = st.compare_old(rr[:3], p)
        check((not ok) and any("runs only in the old file" in x or "cells that differ: 0" in x for x in lines), "compare_old: a missing row means NOT IDENTICAL")
        lines, ok = st.compare_old(rr, tmp / "none.csv")
        check(not ok and "not found" in lines[-1], "compare_old: a missing old file is said plainly")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ============================================================================================ synthetic lab
GAUGE_VOL = {("WY2017", "Dec 21-25"): 0.0, ("WY2017", "Dec 30-Jan 2"): 0.819931, ("WY2017", "Jan 19-21"): 5.52942,
             ("WY2017", "Feb 18-19"): 1.30391, ("WY2017", "Feb 27-28"): 3.915, ("WY2020", "Nov 19-22"): 0.0148829,
             ("WY2020", "Nov 28-30"): 1.45586, ("WY2020", "Dec 23-28"): 4.20157, ("WY2020", "Feb 10-11"): 0.416057,
             ("WY2020", "Feb 22-23"): 4.15998, ("WY2020", "Mar 10-13"): 5.95779, ("WY2020", "Mar 18-19"): 3.32388}
SHIFT = {"retuned": 2, "original": 1, "kchan05": 1, "kchan02": 1}          # model later than gauge by this many hours (others 3)
FINE_STEMS = {"Sycamore_WY2020_Retuned", "Sycamore_WY2017sep_Original"}    # written every 3.75 min; the others every 15 min
OUTAGES = {"WY2017": [("2017-01-01 00:00", "2017-01-02 00:00")], "WY2020": [("2019-12-29 00:00", "2019-12-31 00:00")]}
GAUGE_SPAN = {"WY2017": ("2016-12-01", "2017-04-01"), "WY2020": ("2019-11-01", "2020-04-30")}
PULSE_OFFSET, RISE, FALL = 18, 6, 24
RATIO_IDX = {k: i for i, k in enumerate(st.RATIO_STORMS)}
OUTAGE_KEY = ("WY2017", "Dec 30-Jan 2")


def scale_to(unit, lo, hi, target, hours=None):
    """Scale a model pulse so its volume (over the given hours, if any) equals target (M m3)."""
    h = st.hourly(unit, lo, hi)
    if hours is not None:
        h = h.reindex(hours).dropna()
    v = st.vol_mm3(h)
    return unit * (target / v) if v > 0 else unit * 0.0


def build_lab(tmp, with_tier=True, jitter=None, write_scores=True):
    """Synthetic lab. Returns (lab path, info). jitter: {(year, run key): factor} multiplies one outlet file (all stems of that key)."""
    base = Path(tmp) / "Sycamore"
    lab, res, obs = base / "lab", base / "lab" / "results", base / "init_data" / "obs"
    for d in (lab, res, obs):
        d.mkdir(parents=True, exist_ok=True)
    wins = {w.key: w for w in st.load_step0_windows(WINDOWS_CSV)}
    info = {"gfac": {}, "ghours": {}, "wins": wins}
    for yr in ("WY2017", "WY2020"):
        idx = pd.date_range(*GAUGE_SPAN[yr], freq="15min", inclusive="left")
        mask = pd.Series(False, index=idx)
        for a, b in OUTAGES[yr]:
            mask[(idx >= pd.Timestamp(a)) & (idx < pd.Timestamp(b))] = True
        gq = pd.Series(0.0, index=idx)
        for w in [x for x in wins.values() if x.year == yr]:
            unit = pulse(idx, w.lo + PULSE_OFFSET * HOUR, RISE, FALL)
            u = unit.copy()
            u[mask] = np.nan
            v = st.vol_mm3(st.hourly(u, w.lo, w.hi))
            share = 0.5 if w.key == OUTAGE_KEY else 1.0          # the outage storm also has a second pulse after the outage, so the gap comes before the gauge t90
            f = share * GAUGE_VOL[w.key] / v if v > 0 else 0.0
            gq += unit * f
            info["gfac"][w.key] = f
            if w.key == OUTAGE_KEY:
                unit2 = pulse(idx, pd.Timestamp("2017-01-02 12:00"), RISE, FALL)
                u2 = unit2.copy()
                u2[mask] = np.nan
                gq += unit2 * ((1 - share) * GAUGE_VOL[w.key] / st.vol_mm3(st.hourly(u2, w.lo, w.hi)))
        out = gq.copy()
        out[mask] = np.nan
        for w in [x for x in wins.values() if x.year == yr]:
            info["ghours"][w.key] = st.hourly(out, w.lo, w.hi).index
        tier = np.array(["measured"] * len(idx), dtype=object)
        tier[mask.to_numpy()] = "missing"
        if yr == "WY2017":        # known quality: in the Feb 27-28 window mark 20 readings estimated and 10 filled_interp
            w = wins[("WY2017", "Feb 27-28")]
            pos = np.nonzero(np.asarray((idx >= w.lo) & (idx < w.hi)))[0]
            tier[pos[100:120]] = "estimated"
            tier[pos[200:210]] = "filled_interp"
        df = pd.DataFrame({"datetime": idx.strftime("%Y-%m-%d %H:%M"), "q_cms": out.to_numpy(), "q_cfs": out.to_numpy() * st.CFS_PER_CMS})
        if with_tier:
            df["tier"] = tier
        name = "sycamore_09510200_obs_flow_WY2017.csv" if yr == "WY2017" else "sycamore_09510200_obs_flow_Sep2019_to_May2020.csv"
        df.to_csv(obs / name, index=False)
    allruns = {r.fname: r for r in st.build_runs("sep1", True) + st.build_runs("oct1", True)}
    for r in allruns.values():
        step = 1 / 16 if r.stem in FINE_STEMS else 1 / 4
        ti = np.arange(0, r.hours + step / 2, step)
        mi = r.start + pd.to_timedelta(ti, unit="h")
        mq = pd.Series(0.0, index=mi)
        shift = SHIFT.get(r.key, 3)
        for w in [x for x in wins.values() if x.year == r.year]:
            unit = pulse(mi, w.lo + (PULSE_OFFSET + shift) * HOUR, RISE, FALL)
            if r.key in ("retuned", "original") and w.key in RATIO_IDX:
                mq += scale_to(unit, w.lo, w.hi, st.RATIO_TYPED[r.key][RATIO_IDX[w.key]] * GAUGE_VOL[w.key], info["ghours"][w.key])
            else:
                mq += scale_to(unit, w.lo, w.hi, 1.0)
        if jitter and (r.year, r.key) in jitter:
            mq = mq * jitter[(r.year, r.key)]
        out = pd.DataFrame({"Time_hr": ti.round(4), "Q": mq.to_numpy(), "Qs": mq.to_numpy() * 0.5})
        out.to_csv(res / r.fname, index=False, float_format="%.6f")
    if write_scores:
        hdr = "setting,window,kind,name,years,runs,start,end,vol_gauge_Mm3,vol_model_Mm3,pbias_pct"
        rows = []
        for key, setting in (("retuned", "Retuned"), ("original", "Original")):
            for (yr, nm), i in RATIO_IDX.items():
                g = GAUGE_VOL[(yr, nm)]
                rows.append(f"{setting},storm_x,storm,{nm},{yr[2:]},x,,,{g},{st.RATIO_TYPED[key][i] * g},")
            rows.append(f"{setting},storm01_2017_Dec21-25,storm,Dec 21-25,2017,x,,,0,1.2,")
            rows.append(f"{setting},AllStorms,all_storms,AllStorms,2017+2020,x,,,31.0,30.0,")
        rows.append("Retuned_old,storm_x,storm,Dec 30-Jan 2,2017,x,,,0.819931,0.473105,")
        (lab / "window_scores_R4.csv").write_text("# synthetic stand-in in the format of score_windows_sycamore.py\n" + hdr + "\n" + "\n".join(rows) + "\n")
    (lab / "score_run_sycamore.py").write_text("def arrival(q, level):\n    return None, 'never'\n\ndef main():\n    pg = 1.0\n    level = 0.2 * pg\n    return arrival(None, level)\n")
    shutil.copy(SCRIPT, lab / "score_timing_sycamore.py")
    shutil.copy(WINDOWS_CSV, lab / "step0_windows.csv")
    return lab, info


def run(lab, *args):
    p = subprocess.run([sys.executable, str(lab / "score_timing_sycamore.py"), *args], cwd=lab, capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def csv(path):
    return pd.read_csv(path, comment="#", dtype=str, keep_default_na=False)


def tree_hash(root, skip=()):
    h = {}
    for p in sorted(Path(root).rglob("*")):
        if p.is_file() and not any(s in p.parts for s in skip):
            h[str(p.relative_to(root))] = hashlib.md5(p.read_bytes()).hexdigest()
    return h


def fl(x):
    return float(x) if x not in ("", None) else float("nan")


# ============================================================================================ integration tests
def integration_main(tmp):
    print("\nINTEGRATION TESTS 1: dry run, pass 1, pass 2 (Step 0 windows, 1 September runs, no floor)")
    lab, info = build_lab(tmp)
    wins = info["wins"]
    # two "old" output folders that must not be touched (stand-ins for the folders she will rename)
    for d in ("timing_out", "timing_out_all_storms"):
        (lab / d).mkdir()
        (lab / d / "pass2_timing_long.csv").write_text("# score_timing_sycamore.py v1.3 (stand-in)\nwater_year,run_key,storm\nWY2017,retuned,Dec 21-25\n")
    before = tree_hash(lab.parent, skip=("timing_out_step0",))
    # --- dry run
    rc, out = run(lab, "--dry-run")
    check(rc == 0 and "[DRY RUN]" in out and "DRY RUN: 12 windows, 20 of 20 run files found, every length matches" in out, "dry run: exit 0, 12 windows, 20 of 20 run files found, every length matches")
    check("WINDOW TABLE" in out and "RUN INVENTORY" in out and out.index("WINDOW TABLE") < out.index("VOLUME CHECK"), "dry run prints the window table and the run inventory before the volume check")
    check(all(f"{w.lo:%Y-%m-%d %H:%M}" in out and f"{w.hi:%Y-%m-%d %H:%M}" in out for w in wins.values()), "dry run: every window's start and end is in the table")
    check(out.count("Time_hr 0 =") == 1 and "2016-09-01 00:00" in out and "2019-09-01 00:00" in out and "6551.0000" in out and "6575.0000" in out,
          "dry run: the inventory shows each run's start clock, first and last Time_hr (6551 and 6575)")
    check("Screen 1 (S1)" in out and "Screen 2 (S5)" in out and "Screen 8 (S11)" in out and "Kchan" not in out.split("RUN INVENTORY")[1].split("VOLUME CHECK")[0],
          "dry run: screens are labeled Screen 1 (S1) ... Screen 8 (S11); the Kchan runs are not in the list")
    check("candidates" in out and "<- used as WY2017 S1" in out and "(not in the run list)" in out, "dry run: lists the WY2017 outlet files in results/ for confirmation, marking which are used")
    check("VOLUME CHECK: MATCH" in out and "20 run-and-storm values reproduce the record" in out and "window_scores_R4.csv" in out, "dry run: the ratio check MATCHes against the window_scores source (20 values)")
    check(not (lab / "timing_out_step0").exists(), "dry run writes nothing")
    # --- pass 1
    rc, out = run(lab, "--pass", "1")
    check(rc == 0, "pass 1 exits 0" + ("" if rc == 0 else "\n" + out[-1500:]))
    p1 = lab / "timing_out_step0" / "pass1_gauge_table.csv"
    txt1 = p1.read_text()
    d = csv(p1)
    check(len(d) == 12 and txt1.startswith("# score_timing_sycamore.py v1.4"), "pass 1: 12 rows; the file starts with the v1.4 header")
    check("WINDOW SET USED IN THIS FILE: Step 0 windows" in txt1 and "rain_" not in "".join(d.columns) and "RAIN COLUMNS: the rain_* columns of v1.3 are DROPPED" in txt1,
          "pass 1 header names the window set and says the rain_* columns are dropped; there are none")
    row = {(r.water_year, r.storm): r for r in d.itertuples()}
    check(all(int(row[k].window_hours) == w.hours and row[k].window_start == f"{w.lo:%Y-%m-%d %H:%M}" and row[k].window_end == f"{w.hi:%Y-%m-%d %H:%M}" for k, w in wins.items()),
          "pass 1: window start, end and hours of all 12 storms equal the table")
    check(row[("WY2017", "Dec 21-25")].no_gauge_response == "yes" and row[("WY2020", "Nov 19-22")].no_gauge_response == "yes" and sum(r.no_gauge_response == "yes" for r in d.itertuples()) == 2,
          "pass 1: no_gauge_response = yes for exactly Dec 21-25 2016 and Nov 19-22 2019")
    r = row[("WY2017", "Dec 21-25")]
    check(r.gauge_t50_h == "" and "volume is zero" in r.blank_reasons, "Dec 21-25: gauge volume 0, times blank with the reason (not zero)")
    r = row[("WY2017", "Dec 30-Jan 2")]
    check(r.jan1_2_outage_flag != "" and r.gauge_hours_with_data == "129" and "Jan 01 00:00-Jan 02 00:00 (24 h)" in r.gauge_missing_intervals and r.gap_before_t90 == "yes",
          "outage storm is scored: flag set, 129 of 153 gauge hours, the 24 h gap listed, gap before t90")
    r = row[("WY2020", "Dec 23-28")]
    check("Dec 29 00:00-Dec 31 00:00 (48 h)" in r.gauge_missing_intervals and r.gap_before_t90 == "no" and r.window_hours == "210", "Dec 23-28: the Dec 29-31 dropout is listed, after t90; 210 h")
    r = row[("WY2017", "Feb 27-28")]
    check(abs(fl(r.pct_interp_or_est) - 100 * 30 / 484) < 1e-3, f"Feb 27-28: 30 of 484 readings interpolated/estimated = {100 * 30 / 484:.3f}% (got {r.pct_interp_or_est})")
    okp = True
    for k, w in wins.items():
        fac = info["gfac"][k]
        if fac == 0:
            continue
        rr = row[k]
        okp &= abs(fl(rr.gauge_peak_cfs) - fac * st.CFS_PER_CMS) <= 1e-3 * max(1, fac * st.CFS_PER_CMS) and rr.gauge_peak_time == f"{w.lo + pd.Timedelta(hours=PULSE_OFFSET + RISE):%Y-%m-%d %H:%M}"
    check(okp, "pass 1: gauge peak (cfs) and its time equal the construction for every storm")
    # --- pass 2
    rc, out = run(lab, "--pass", "2")
    check(rc == 0, "pass 2 exits 0 with no --floor-mm3" + ("" if rc == 0 else "\n" + out[-2000:]))
    p2 = lab / "timing_out_step0" / "pass2_timing_long.csv"
    pm = lab / "timing_out_step0" / "pass2_medians.csv"
    d2, dm = csv(p2), csv(pm)
    txt2 = p2.read_text()
    check(len(d2) == 10 * 5 + 10 * 7, f"pass 2: 120 rows (10 runs x 5 storms + 10 runs x 7 storms), got {len(d2)}")
    check(not any(c.startswith("rain_") for c in d2.columns) and {"no_gauge_response", "vol_ratio_gauge_hours", "vol_ratio_model_over_gauge"} <= set(d2.columns), "pass 2: no rain_* columns; the new columns are there")
    check(txt2.startswith("# score_timing_sycamore.py v1.4") and "WINDOW SET USED IN THIS FILE: Step 0 windows" in txt2 and "FLOOR: none." in txt2 and "WINDOW TABLE:" in txt2 and "RUN FILE:" in txt2,
          "pass 2 header: version, window set, 'FLOOR: none', the window table and one RUN FILE line per run")
    check("THIS FILE (pass2_medians.csv)" in pm.read_text() and "WINDOW SET USED IN THIS FILE: Step 0 windows" in pm.read_text(), "the medians file's header names the window set too")
    check(len(dm) == 10 * 3 * 3 * 8, f"medians file: 720 rows (10 runs x 3 year scopes x 3 storm sets x 8 measures), got {len(dm)}")
    sw = csv(lab / "timing_out_step0" / "step0_windows.csv")
    check(len(sw) == 12 and list(sw.columns) == ["water_year", "storm", "start", "end", "hours", "gauge_response"] and list(sw.hours.astype(int)) == [188, 153, 85, 135, 121, 160, 141, 210, 132, 121, 187, 120],
          "step0_windows.csv is saved in the output folder with the hours")
    check(set(d2.window_rule) == {"Step 0 window (step0_windows.csv)"}, "every row says which window set it used")
    r2 = {(r.water_year, r.run_key, r.storm): r for r in d2.itertuples()}
    # exactness: model = gauge shape shifted by 1, 2 or 3 h on the same grid
    n_exact = n_fine = 0
    ok_exact = ok_fine = True
    for (y, k, n), r in r2.items():
        if n in ("Dec 21-25", "Dec 30-Jan 2"):
            continue
        sh = SHIFT.get(k, 3)
        stem = r.run_stem
        tt = [fl(r.t10_model_minus_gauge_h), fl(r.t50_model_minus_gauge_h), fl(r.t90_model_minus_gauge_h)]
        good = abs(fl(r.peak_time_model_minus_gauge_h) - sh) < 1e-6 and int(r.xcorr_best_lag_h) == sh
        if stem in FINE_STEMS:
            n_fine += 1
            ok_fine &= good and all(sh - 0.12 < x < sh for x in tt) and fl(r.xcorr_best_r) > 0.999
        else:
            n_exact += 1
            ok_exact &= good and all(abs(x - sh) < 1e-6 for x in tt) and fl(r.xcorr_best_r) > 0.99999
    check(ok_exact and n_exact > 80, f"model = gauge shape shifted 1, 2 or 3 h: t10, t50, t90, peak time and lag equal the shift exactly ({n_exact} rows)")
    check(ok_fine and n_fine >= 8, f"same, model written every 3.75 min: peak time and lag exact; volume times 0.09 h early (sampling convention) ({n_fine} rows)")
    tl = [(r.tail_gauge_Mm3, r.gauge_vol_Mm3, r.model_tail_fraction, r.tail_volume_ratio) for r in r2.values() if r.tail_gauge_Mm3 != ""]
    check(len(tl) == 110 and all(abs(fl(t) - 0.1 * fl(g)) < 2e-4 for t, g, _, _ in tl) and all(fl(f) >= 0 and fl(q) >= 0 for _, _, f, q in tl)
          and all(fl(r.model_tail_fraction) > 0 and fl(r.tail_volume_ratio) > 0 for r in r2.values() if r.tail_gauge_Mm3 != "" and r.storm != "Dec 30-Jan 2"),
          "tail: the gauge tail after its own t90 is 10% of the gauge volume in all 110 scored rows (every row except the zero-volume storm); tail ratio and model tail fraction are filled (0 for the outage storm, whose model has no flow after the outage)")
    r = r2[("WY2017", "retuned", "Dec 21-25")]
    check(r.no_gauge_response == "yes" and r.t50_model_minus_gauge_h == "" and "volume is zero" in r.blank_reasons and r.below_floor == "",
          "Dec 21-25 (gauge volume 0): no_gauge_response yes, blanks carry reasons, no floor is mentioned")
    r = r2[("WY2020", "retuned", "Nov 19-22")]
    check(r.no_gauge_response == "yes" and r.below_floor == "" and r.t50_model_minus_gauge_h != "" and r.tail_volume_ratio != "",
          "Nov 19-22 (gauge 0.015 M m3, no gauge response): NOT blanked for small volume, every timing number is given")
    r = r2[("WY2020", "retuned", "Dec 23-28")]
    check(r.t50_model_minus_gauge_h != "" and fl(r.model_vol_Mm3) < 0.3 and r.below_floor == "", "Dec 23-28 (model 5% of the gauge): timing still given, no floor")
    r = r2[("WY2017", "retuned", "Dec 30-Jan 2")]
    check(r.jan1_2_outage_flag != "" and r.gap_before_t90 == "yes" and abs(fl(r.vol_ratio_gauge_hours) - 0.58) < 5e-4 and fl(r.vol_ratio_model_over_gauge) > fl(r.vol_ratio_gauge_hours) + 0.05,
          f"outage storm: scored, vol_ratio_gauge_hours = 0.58 (matched hours), vol_ratio_model_over_gauge larger ({r.vol_ratio_model_over_gauge}), flags kept")
    # the thesis volume ratio and the check on all ten storms, Re-tuned and Original
    okr = all(abs(fl(r2[(y, k, n)].vol_ratio_gauge_hours) - st.RATIO_TYPED[k][i]) < 6e-4 for (y, n), i in RATIO_IDX.items() for k in ("retuned", "original"))
    check(okr, "vol_ratio_gauge_hours of Re-tuned and Original equals the per-storm ratios on record for all 10 storms")
    chk = {(r.water_year, r.run_key, r.storm): r.volume_check for r in d2.itertuples()}
    check(all(chk[(y, k, n)] == "MATCH" for (y, n) in RATIO_IDX for k in ("retuned", "original")), "volume_check = MATCH on the 20 checked run-and-storm values")
    check(all(v == "NOT CHECKED" for (y, k, n), v in chk.items() if k not in ("retuned", "original") or (y, n) not in RATIO_IDX), "volume_check = NOT CHECKED for the screens and the two no-response storms")
    check(all(fl(r.model_vol_gauge_hours_Mm3) <= fl(r.model_vol_Mm3) + 1e-9 for r in d2.itertuples() if r.model_vol_gauge_hours_Mm3 != ""), "model volume over the gauge's hours never exceeds the model's whole-window volume")
    check(set(d2[d2.water_year == "WY2020"].role) == {"validation", "post_validation_screen"} and set(d2[d2.water_year == "WY2017"].role) == {"calibration", "post_validation_screen"}
          and set(d2[d2.run_key.isin(["retuned", "original"]) & (d2.water_year == "WY2020")].fixed_before_wy2020_results) == {"yes"}, "role and fixed_before_wy2020_results columns are as in v1.3 (screens in both years are post_validation_screen)")
    # medians file: independent recomputation
    ok_med = True
    for key in ("retuned", "original", "S5"):
        for ml, col in st.MEASURES:
            for sid, sel in (("a", lambda x: x.no_gauge_response == "no"), ("b", lambda x: (x.no_gauge_response == "no") & ~((x.water_year == "WY2017") & (x.storm == "Dec 30-Jan 2"))), ("c", lambda x: x.water_year != "none")):
                sub = d2[(d2.run_key == key)]
                sub = sub[sel(sub)]
                v = pd.to_numeric(sub[col], errors="coerce").dropna()
                m = dm[(dm.run_key == key) & (dm.years == "both years") & (dm.storm_set == sid) & (dm.measure == ml)]
                ok_med &= len(m) == 1 and int(m.n.iloc[0]) == len(v) and (len(v) == 0 or abs(fl(m["median"].iloc[0]) - v.median()) < 1e-3)
    check(ok_med, "medians file: median and n for 3 runs x 8 measures x 3 storm sets (both years) equal an independent pandas recomputation from the long file")
    check("MEDIANS (pass2_medians.csv rows for Re-tuned and Original" in out and "  Re-tuned" in out and "  Original" in out and "(a) storms with a gauge response" in out, "console prints the Re-tuned and Original medians and the storm-set key")
    check("CHECK 5" in out and "Wrote 720 rows" in out, "CHECK 5 and the medians file are reported")
    # old folders untouched
    check(tree_hash(lab.parent, skip=("timing_out_step0",)) == before, "nothing outside timing_out_step0/ changed (md5 of every input, script and the two old output folders)")
    rc, out = run(lab, "--pass", "2")
    check(rc == 0, "re-running pass 2 overwrites its own v1.4 files")
    (lab / "timing_out_step0" / "pass2_medians.csv").write_text("a,b\n1,2\n")
    rc, out = run(lab, "--pass", "2")
    check(rc != 0 and "was not written by score_timing_sycamore.py v1.4" in out, "refuses to overwrite a file that is not its own v1.4 output")
    check(len(csv(p2)) == 120 and txt2 == p2.read_text(), "...and it checked before writing anything: the long file it had already written was left as it was")
    shutil.rmtree(lab / "timing_out_step0")
    return lab, info


def clone(lab, tmp, name):
    dst = Path(tmp) / name
    shutil.copytree(lab.parent, dst)
    return dst / "lab"


def integration_variants(tmp, lab, info):
    print("\nINTEGRATION TESTS 2: length guard, ratio check, floor, Kchan, missing files, other options")
    # --- A. a wrong-length file for a run
    la = clone(lab, tmp, "a")
    res = la / "results"
    shutil.copy(res / "W5_d012_k10_f0p15_as300_ksx1p03_nocl_Outlet.qout", res / "Sycamore_WY2017sep_Retuned_Outlet.qout")
    rc, out = run(la, "--dry-run")
    check(rc != 0 and "LENGTH CHECK FAILED" in out and "8,759 h is the old 1 October whole-year run" in out and "LENGTH DOES NOT MATCH" in out,
          "a 1 October 8,759 h file in the place of a 1 September WY2017 run: the dry run says so (and exits non-zero)")
    rc, out = run(la, "--pass", "2")
    check(rc != 0 and "STOPPED" in out and "shift every time" in out and not (la / "timing_out_step0").exists(), "the same wrong file in a real run: stops with a clear message, nothing written")
    df = pd.read_csv(res / "Sycamore_WY2020_S1_Ks1p0_Outlet.qout")
    df[df.Time_hr <= 6000].to_csv(res / "Sycamore_WY2020_S1_Ks1p0_Outlet.qout", index=False)
    shutil.copy(res / "Sycamore_WY2017sep_S1_Ks1p0_Outlet.qout", res / "Sycamore_WY2017sep_Retuned_Outlet.qout")      # fix the first problem
    rc, out = run(la, "--pass", "2")
    check(rc != 0 and "WY2020 Screen 1 (S1)" in out and "last Time_hr 6000" in out, "a WY2020 file cut short at 6,000 h: stops and names the run and the length found")
    rc, out = run(la, "--pass", "2", "--ignore-check-failure", "--out-dir", "timing_out_ign")
    check("--ignore-check-failure was given" in out, "--ignore-check-failure can carry on past a length problem (diagnosis only, it says so)")
    # --- B. Kchan: start is per run (1 October WY2017 Kchan files next to 1 September runs)
    rc, out = run(lab, "--pass", "2", "--include-kchan", "--out-dir", "timing_out_k")
    dk = csv(lab / "timing_out_k" / "pass2_timing_long.csv")
    check(rc == 0 and len(dk) == 12 * 5 + 12 * 7, "--include-kchan: 24 runs, 144 rows" + ("" if rc == 0 else "\n" + out[-1500:]))
    k = dk[(dk.run_key == "kchan05") & (dk.water_year == "WY2017") & (dk.storm == "Jan 19-21")].iloc[0]
    check(k.run_stem == "W5_d012_k10_base_kchan0p5" and abs(fl(k.peak_time_model_minus_gauge_h) - 1.0) < 1e-6 and abs(fl(k.t50_model_minus_gauge_h) - 1.0) < 1e-6 and k.volume_check == "NOT CHECKED",
          "a WY2017 Kchan run (1 October start, 8,759 h) scored beside 1 September runs: times exact, so its own start was used")
    check("RUN FILE: WY2017 kchan05" in (lab / "timing_out_k" / "pass2_timing_long.csv").read_text() and "Time_hr 0 = 2016-10-01 00:00" in (lab / "timing_out_k" / "pass2_timing_long.csv").read_text(),
          "the CSV header records each run's start clock and file")
    # --- C. ratio check: mismatch stops; tolerance edge
    lj, _ = build_lab(Path(tmp) / "j", jitter={("WY2020", "original"): 1.01})
    rc, out = run(lj, "--dry-run")
    check(rc != 0 and "VOLUME CHECK: MISMATCH" in out and "WY2020" in out and "original" in out, "Original WY2020 volumes 1% high: the dry run reports MISMATCH")
    rc, out = run(lj, "--pass", "2")
    check(rc != 0 and "STOPPED" in out and "the volumes do not reproduce" in out and not (lj / "timing_out_step0").exists(), "ratio mismatch in a real run: stops, writes nothing")
    rc, out = run(lj, "--pass", "2", "--ignore-check-failure")
    dj = csv(lj / "timing_out_step0" / "pass2_timing_long.csv")
    check(rc == 0 and "MISMATCH" in set(dj.volume_check) and "MISMATCH, carried on" in (lj / "timing_out_step0" / "pass2_timing_long.csv").read_text(),
          "--ignore-check-failure carries on; the rows say MISMATCH and the header says so")
    lj2, _ = build_lab(Path(tmp) / "j2", jitter={("WY2020", "original"): 1.001, ("WY2017", "retuned"): 0.999})
    rc, out = run(lj2, "--dry-run")
    check(rc == 0 and "VOLUME CHECK: MATCH" in out, "a 0.1% change in volume stays inside the 0.01 tolerance: MATCH")
    # --- D. the record: source file, typed values, bad sources
    lt = clone(lab, tmp, "t")
    (lt / "window_scores_R4.csv").unlink()
    rc, out = run(lt, "--dry-run")
    check(rc == 0 and "typed per-storm ratios" in out and "VOLUME CHECK: MATCH" in out, "no window_scores file: the typed ratios are used, the printout says so, and they MATCH")
    ls = clone(lab, tmp, "s")
    sc = []
    for x in (ls / "window_scores_R4.csv").read_text().splitlines():
        if x.startswith("Original,storm_x,storm,Feb 22-23,2020"):
            x = "Original,storm_x,storm,Feb 22-23,2020,x,,,4.15998,6.23997,"        # ratio 1.5 instead of 1.20
        sc.append(x)
    (ls / "window_scores_R4.csv").write_text("\n".join(sc) + "\n")
    rc, out = run(ls, "--dry-run")
    check(rc != 0 and "WY2020 Feb 22-23 original" in out and "MISMATCH" in out, "a changed ratio in the source file: that storm and run are named in the MISMATCH list")
    rc, out = run(lab, "--dry-run", "--ratio-source", "nofile.csv")
    check(rc != 0 and "cannot find the ratio source" in out, "--ratio-source pointing at nothing: stops")
    lines = [x for x in (lab / "window_scores_R4.csv").read_text().splitlines() if "Mar 18-19" not in x]
    (ls / "window_scores_R4.csv").write_text("\n".join(lines) + "\n")
    rc, out = run(ls, "--dry-run")
    check(rc != 0 and "has no storm row for" in out and "Mar 18-19" in out, "a source file missing a storm: stops and names it")
    # --- E. the floor is optional
    rc, out = run(lab, "--pass", "2", "--floor-mm3", "0.05", "--out-dir", "timing_out_floor")
    df_ = csv(lab / "timing_out_floor" / "pass2_timing_long.csv")
    check(rc == 0 and (df_.below_floor == "gauge").sum() == 20 and set(df_[df_.below_floor == "gauge"].storm) == {"Dec 21-25", "Nov 19-22"}
          and "FLOOR: 0.05 M m3 (your choice)" in (lab / "timing_out_floor" / "pass2_timing_long.csv").read_text(),
          "--floor-mm3 0.05 still works when asked for: 20 rows (the two small storms x 10 runs) blank their timing; the header states the floor")
    rc, out = run(lab, "--pass", "2", "--model-floor-mm3", "0.05", "--out-dir", "timing_out_x")
    check(rc != 0 and "--floor-mm3" in out, "--model-floor-mm3 without --floor-mm3 stops and says why")
    # --- F. missing files
    lm = clone(lab, tmp, "m")
    (lm / "results" / "Sycamore_WY2017sep_S10_Ks1p0_GW0p5_Outlet.qout").unlink()
    (lm / "results" / "Sycamore_WY2020_S1_Ks1p0_Outlet.qout").unlink()
    rc, out = run(lm, "--pass", "2")
    dmz = csv(lm / "timing_out_step0" / "pass2_timing_long.csv")
    check(rc == 0 and "NOT FOUND" in out and "Screen 7 (S10)" in out and "Screen 1 (S1)" in out and len(dmz) == 9 * 5 + 9 * 7, "two outlet files missing: named as NOT FOUND, the other runs are scored (108 rows)")
    (lm / "results" / "Sycamore_WY2020_Retuned_Outlet.qout").unlink()
    rc, out = run(lm, "--pass", "2", "--out-dir", "timing_out_y")
    check(rc != 0 and "both Re-tuned and Original are needed for WY2020" in out, "WY2020 Re-tuned missing: stops (the ratio check covers both runs)")
    rc, out = run(lm, "--dry-run")
    check(rc != 0 and "WOULD STOP" in out and "NOT FOUND" in out, "the dry run reports the same thing as WOULD STOP instead of stopping")
    # --- G. other options
    ln = clone(lab, tmp, "n")
    shutil.copy(ln / "results" / "Sycamore_WY2017sep_Retuned_Outlet.qout", ln / "results" / "my_ret_Outlet.qout")
    (ln / "results" / "Sycamore_WY2017sep_Retuned_Outlet.qout").unlink()
    rc, out = run(ln, "--dry-run", "--stem", "2017:retuned=my_ret_Outlet.qout")
    check(rc == 0 and "my_ret" in out and "<- used as WY2017 retuned" in out, "--stem names another outlet file for a run (the dry run shows it, and the candidate list marks it as used)")
    rc, out = run(lab, "--definitions")
    check(rc == 0 and "WINDOW SET USED IN THIS FILE: Step 0 windows" in out and "RAIN COLUMNS" in out and "vol_ratio_gauge_hours" in out, "--definitions prints the definitions without any files")
    os.rename(lab / "step0_windows.csv", lab / "step0_windows.csv.off")
    rc, out = run(lab, "--dry-run")
    os.rename(lab / "step0_windows.csv.off", lab / "step0_windows.csv")
    check(rc != 0 and "cannot find the windows table" in out, "no step0_windows.csv: stops with a plain message")
    lnt, _ = build_lab(Path(tmp) / "nt", with_tier=False)
    rc, out = run(lnt, "--pass", "1")
    d4 = csv(lnt / "timing_out_step0" / "pass1_gauge_table.csv")
    check(rc == 0 and set(d4.pct_interp_or_est) == {""} and "no 'tier' column" in d4.blank_reasons.iloc[0], "gauge file without a tier column: percentages blank with a reason, no crash")
    # whitespace-separated outlet file
    p = Path(tmp) / "ws.qout"
    p.write_text("Time_hr Q\n0.0 0.0\n0.0625 1.5\n0.125 2.5\n")
    fr = st.read_qout_frame(p)
    check(list(fr.Q) == [0.0, 1.5, 2.5] and list(fr.Time_hr) == [0.0, 0.0625, 0.125], "an outlet file separated by spaces is read the same way as a comma-separated one")


def integration_old_rules(tmp, lab, info):
    print("\nINTEGRATION TESTS 3: --windows as-scored (v1.3 behavior) and --compare-old")
    # The synthetic lab holds volumes for the Step 0 windows, not the old record, so the old-rule volume check must MISMATCH here;
    # that is itself the test that the old check still runs and still stops.
    rc, out = run(lab, "--pass", "2", "--windows", "as-scored", "--wy2017-runs", "oct1", "--include-kchan", "--floor-mm3", "0", "--out-dir", "timing_out_stop")
    check(rc != 0 and "VOLUME CHECK: MISMATCH" in out and "STOPPED" in out and not (lab / "timing_out_stop").exists(),
          "as-scored on volumes that are not the old record: the v1.3 volume check MISMATCHes, stops, writes nothing")
    args = ["--pass", "2", "--windows", "as-scored", "--wy2017-runs", "oct1", "--include-kchan", "--floor-mm3", "0", "--ignore-check-failure"]
    rc, out = run(lab, *args, "--out-dir", "timing_out_old_a")
    pa = lab / "timing_out_old_a" / "pass2_timing_long.csv"
    da = csv(pa)
    check(rc == 0 and len(da) == 4 * 5 + 12 * 7, "as-scored, oct1 runs, Kchan, floor 0: 104 rows (4 WY2017 runs x 5 storms + 12 WY2020 runs x 7), as in the 20 Nov file" + ("" if rc == 0 else "\n" + out[-1500:]))
    w = {(r.water_year, r.storm): r for r in da.itertuples() if r.run_key == "retuned"}
    check(w[("WY2017", "Jan 19-21")].window_end == "2017-01-23 00:00" and w[("WY2020", "Nov 28-30")].window_end == "2019-12-04 00:00" and w[("WY2017", "Dec 21-25")].window_start == "2016-12-21 00:00",
          "as-scored windows are the v1.3 ones (start 00:00 on the first rain day; Jan 19-21 cut at Jan 23; Nov 28-30 uncut to Dec 4)")
    txt = pa.read_text()
    check("WINDOW SET USED IN THIS FILE: old rule, as-scored" in txt and set(da.window_rule) <= {"last rain day + 4 d", "cut at the next storm's first day"},
          "the old-rule CSV header says which window set was used; window_rule is the v1.3 wording")
    check(set(da[da.water_year == "WY2017"].run_stem) == {"W5_d012_k10_f0p15_as300_ksx1p03_nocl", "W5_d012_k10_base_nocl", "W5_d012_k10_base_kchan0p5", "W5_d012_k10_base_kchan0p2"},
          "--wy2017-runs oct1 reads the four old whole-year WY2017 files (Time_hr 0 = 1 October)")
    k = da[(da.run_key == "retuned") & (da.water_year == "WY2017") & (da.storm == "Feb 27-28")].iloc[0]
    check(abs(fl(k.peak_time_model_minus_gauge_h) - 2.0) < 1e-6 and int(k.xcorr_best_lag_h) == 2, "old whole-year file (1 October start): peak time and lag equal the 2 h construction, so the start used was right")
    check("rain_" not in "".join(da.columns), "the old-rule run has no rain_* columns either (dropped in v1.4, for every window set)")
    # default 1 September runs under the old rules: their model volumes are not on record
    rc, out = run(lab, "--pass", "2", "--windows", "as-scored", "--floor-mm3", "0", "--ignore-check-failure", "--out-dir", "timing_out_old_b")
    check(rc == 0 and "Model volumes NOT CHECKED for" in out and "WY2017 retuned" in out, "as-scored with the 1 September WY2017 runs: their model volumes are marked NOT CHECKED (the record belongs to the old runs)")
    # compare-old: identical on a second run, different between window sets
    rc, out = run(lab, *args, "--out-dir", "timing_out_old_c", "--compare-old", str(pa))
    check(rc == 0 and "COMPARE WITH OLD FILE" in out and "RESULT: IDENTICAL on every compared cell" in out and "cells that differ: 0" in out,
          "--compare-old: a second as-scored run reproduces the first line by line (IDENTICAL, 104 rows)")
    rc, out = run(lab, "--pass", "2", "--include-kchan", "--wy2017-runs", "oct1", "--out-dir", "timing_out_new_cmp", "--compare-old", str(pa))
    check(rc == 0 and "NOT IDENTICAL" in out and "cells that differ" in out and "window_start" in out, "--compare-old against a different window set: NOT IDENTICAL, with the differing cells listed" + ("" if rc == 0 else "\n" + out[-1200:]))
    # an old v1.3-style file: it has the rain_* columns and lacks the new ones
    old13 = pd.DataFrame([{"water_year": "WY2017", "run_key": "retuned", "run_label": "Re-tuned", "storm": "Jan 19-21", "rain_basin_mm": "63.4", "t50_model_minus_gauge_h": "2.0"}])
    p13 = lab / "old_v13_stub.csv"
    with open(p13, "w") as fh:
        fh.write("# score_timing_sycamore.py v1.3 (stand-in)\n")
        old13.to_csv(fh, index=False)
    rc, out = run(lab, *args, "--out-dir", "timing_out_old_d", "--compare-old", str(p13))
    check(rc == 0 and "columns only in the new file" in out and "columns only in the old file (not compared): rain_basin_mm" in out and "runs only in the new file" in out,
          "--compare-old with a v1.3-style file: the columns that exist on one side only are listed, and the rows missing from the old file are named")
    rc, out = run(lab, "--pass", "1", "--compare-old", str(pa))
    check(rc != 0 and "--compare-old compares pass 2 rows" in out, "--compare-old with pass 1: stops and says why")


# ============================================================================================ main and mutation check
def main():
    if "--mutations" in sys.argv:
        return mutations()
    unit_formula_tests()
    unit_v14_tests()
    row_tests()
    medians_tests()
    if os.environ.get("ST_SKIP_INTEGRATION") != "1":
        tmp = tempfile.mkdtemp(prefix="stt_")
        try:
            lab, info = integration_main(tmp)
            if os.environ.get("ST_FAST") != "1":
                integration_variants(tmp, lab, info)
                integration_old_rules(tmp, lab, info)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
    print("\n" + ("ALL TESTS PASSED" if not FAILS else f"{len(FAILS)} TEST(S) FAILED:\n  - " + "\n  - ".join(FAILS)))
    return 0 if not FAILS else 1


# (name, text to change, replacement, how much of the test is needed to notice: unit / fast / full)
MUTATIONS = [
    ("lag sign flipped", "m2.index = m2.index - pd.Timedelta(hours=L)", "m2.index = m2.index + pd.Timedelta(hours=L)", "unit"),
    ("difference sign flipped", 'row[f"{c}_model_minus_gauge_h"] = r4(mt[f] - gt[f])', 'row[f"{c}_model_minus_gauge_h"] = r4(gt[f] - mt[f])', "fast"),
    ("volume uses 3000 s per hour", 'return float(h.sum() * 3600.0 / 1e6) if len(h) else float("nan")', 'return float(h.sum() * 3000.0 / 1e6) if len(h) else float("nan")', "unit"),
    ("tail starts at t50", "t0 = w.lo + pd.Timedelta(hours=gt[0.9])", "t0 = w.lo + pd.Timedelta(hours=gt[0.5])", "fast"),
    ("no interpolation inside the hour", "out[f] = float(t0[k] + (target - prev) / v[k])", "out[f] = float(t0[k])", "unit"),
    ("constant series gives r = 0 instead of blank", 'return None, len(idx), "gauge hourly flow is constant in the window"', 'return 0.0, len(idx), ""', "unit"),
    ("second bump ignores the dip", "if float(v[a:b + 1].min()) <= DIP_FRAC * v[j]:", "if True:", "unit"),
    ("missing hours filled with zero", 'return x.resample("60min").mean().dropna()', 'return x.resample("60min").mean().fillna(0.0)', "unit"),
    ("window includes the end hour", "x = s[(s.index >= lo) & (s.index < hi)]", "x = s[(s.index >= lo) & (s.index <= hi)]", "fast"),
    ("vol_ratio_gauge_hours uses all model hours", 'row["vol_ratio_gauge_hours"] = r4(mgh / gtot)', 'row["vol_ratio_gauge_hours"] = r4(mtot / gtot)', "fast"),
    ("no_gauge_response flipped", 'return "yes" if self.gauge_response == "none" else "no"', 'return "no" if self.gauge_response == "none" else "yes"', "unit"),
    ("a floor is on by default", 'ap.add_argument("--floor-mm3", type=float, default=None,', 'ap.add_argument("--floor-mm3", type=float, default=0.05,', "fast"),
    ("Step 0 hour check switched off", "if hours != STEP0_HOURS[key]:", "if False:", "unit"),
    ("run length guard switched off", "ok = (first <= tol) and abs(last - run.hours) <= tol", "ok = True", "unit"),
    ("start set per year: Kchan WY2017 read from 1 September",
     'elif key in ("kchan05", "kchan02"):\n                runs.append(Run(yr, key, label, STEMS_2017_OCT1[key], OCT1_2016, HOURS_OCT1_2017))',
     'elif key in ("kchan05", "kchan02"):\n                runs.append(Run(yr, key, label, STEMS_2017_OCT1[key], SEP1_2016, HOURS_OCT1_2017))', "full"),
    ("medians set (b) keeps the outage storm", 'return (df["no_gauge_response"] == "no") & ~((df["water_year"] == OUTAGE_STORM[0]) & (df["storm"] == OUTAGE_STORM[1]))', 'return df["no_gauge_response"] == "no"', "unit"),
    ("medians n counts undefined values", '"n": int(len(v))', '"n": int(len(s3))', "unit"),
    ("ratio check tolerance 0.5", "RATIO_TOL = 0.01 ", "RATIO_TOL = 0.5 ", "full"),
    ("Re-tuned not marked fixed before results", '("WY2020", "retuned"): ("validation", "yes"),', '("WY2020", "retuned"): ("validation", "no"),', "unit"),
    ("compare-old sees no differences", "if str(o.loc[k, c]) != str(n.loc[k, c]):", "if False:", "unit"),
    ("percent estimated counts only 'estimated'", 'quality["pct_interp_or_est"] = 100.0 * float(t.isin(["estimated", "filled_interp"]).sum()) / n', 'quality["pct_interp_or_est"] = 100.0 * float(t.isin(["estimated"]).sum()) / n', "fast"),
]


def mutations():
    print("MUTATION CHECK: each broken copy of the script must make this test fail")
    src = SCRIPT.read_text()
    bad = 0
    tmp = Path(tempfile.mkdtemp(prefix="stm_"))
    try:
        wanted = os.environ.get("ST_MUT_LEVELS", "unit,fast,full").split(",")
        for name, old, new, level in MUTATIONS:
            if level not in wanted or os.environ.get("ST_MUT_ONLY", "") not in name:
                continue
            if src.count(old) != 1:
                print(f"  ??    {name}: the text to change was found {src.count(old)} times (expected 1)")
                bad += 1
                continue
            p = tmp / "score_timing_sycamore.py"
            p.write_text(src.replace(old, new))
            env = dict(os.environ, ST_SCRIPT=str(p), ST_SKIP_INTEGRATION="1" if level == "unit" else "0", ST_FAST="1" if level == "fast" else "0")
            r = subprocess.run([sys.executable, str(Path(__file__).resolve())], env=env, capture_output=True, text=True)
            caught = r.returncode != 0
            print(("  CAUGHT  " if caught else "  MISSED  ") + name + f"   (exit {r.returncode}, level {level})")
            bad += 0 if caught else 1
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("\n" + ("ALL MUTATIONS CAUGHT" if not bad else f"{bad} MUTATION(S) NOT CAUGHT"))
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
