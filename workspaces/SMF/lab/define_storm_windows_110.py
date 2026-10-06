"""
define_storm_windows_110.py   (version 1)

PURPOSE
  Fix the four storm windows for the storm-by-storm scoring of the long run (ms1), using the REAL gauge only.
  It never looks at any model result. Windows are fixed before results so they cannot be bent to fit them.

THE RULE (approved by Tracy, 4 Oct 2026, 21:32)
  - 12 Aug keeps the Stage B window: 12 Aug 16:00 to 13 Aug 12:00 (hours 280 to 300 from 1 Aug 00:00),
    so the Stage B numbers carry over.
  - 19 Aug, 8 Sep, 27 Sep: the window starts 2 h before the first real gauge reading at or above 1 cfs and
    ends 6 h after the last real reading at or above 1 cfs. Both ends are then rounded OUTWARD to the whole hour
    (start down, end up).
  - For 12 Aug the same recipe is also printed, as a check only. (It gives a shorter window inside the Stage B one.)
  - A "storm" is a run of real readings at or above 1 cfs. Two such runs belong to the same storm unless they are
    more than 12 h apart.

SAFETY STOPS (nothing is written if any of these fires)
  - The gauge file is missing, or has no readings in the model window.
  - The readings at or above 1 cfs do not form exactly four storms, one starting on each of 12 Aug, 19 Aug, 8 Sep
    and 27 Sep (calendar days of the gauge clock).
  - A window runs outside the model run (1 Aug 00:00 for 1416 h, which ends 29 Sep 00:00).
  - Two windows overlap.

WHAT IT WRITES
  One small new file, storm_windows_110.json, in the folder you run it from. It refuses to overwrite an
  existing one unless you add --overwrite. Nothing else is created or changed.

RUN (from the lab/ folder; the gauge workbook is large, so reading it takes up to a minute)
  python define_storm_windows_110.py

NOTES
  - The gauge logger writes a reading only when the flow changes, so a reading is held until the next one.
    Volumes below are that step-function sum ("hold-last-value", the same idea as the filled arm), in acre-feet.
  - Gauge clock and model clock are taken to be the same (as in the Stage B event window).
"""
import argparse
import json
import sys
from pathlib import Path

import pandas as pd

VERSION = "1"
GAUGE = Path("../init_data/met/SMF_Observations_1993-2025.xlsx")
GAUGE_SHEET = "Discharge"
GAUGE_SKIPROWS = 6
OUT_NAME = "storm_windows_110.json"

MODEL_START = pd.Timestamp("2014-08-01 00:00")
RUNTIME_HOURS = 1416
MODEL_END = MODEL_START + pd.Timedelta(hours=RUNTIME_HOURS)

THRESHOLD_CFS = 1.0
PAD_BEFORE_H = 2
PAD_AFTER_H = 6
EPISODE_GAP_H = 12
STAGE_B_WINDOW = (pd.Timestamp("2014-08-12 16:00"), pd.Timestamp("2014-08-13 12:00"))
STORMS = [("2014-08-12", "12 Aug"), ("2014-08-19", "19 Aug"), ("2014-09-08", "8 Sep"), ("2014-09-27", "27 Sep")]
CFS_TO_ACFT_PER_SEC = 1.0 / 43560.0   # 1 cfs for 1 second is 1/43560 acre-feet


def read_gauge(path):
    obs = pd.read_excel(path, sheet_name=GAUGE_SHEET, skiprows=GAUGE_SKIPROWS)
    obs["datetime"] = pd.to_datetime(obs["Date"].astype(str) + " " + obs["Time"].astype(str))
    obs["cfs"] = pd.to_numeric(obs["cfs"], errors="coerce")
    obs = obs.dropna(subset=["datetime", "cfs"]).sort_values("datetime")
    return obs[["datetime", "cfs"]].reset_index(drop=True)


def find_storms(g):
    """Group readings at or above the threshold into storms (gap of more than EPISODE_GAP_H starts a new one)."""
    hi = g[g["cfs"] >= THRESHOLD_CFS]
    storms, cur = [], []
    for t in hi["datetime"]:
        if cur and (t - cur[-1]) > pd.Timedelta(hours=EPISODE_GAP_H):
            storms.append(cur)
            cur = []
        cur.append(t)
    if cur:
        storms.append(cur)
    return [(s[0], s[-1]) for s in storms]


def recipe_window(first, last):
    ws = (first - pd.Timedelta(hours=PAD_BEFORE_H)).floor("h")
    we = (last + pd.Timedelta(hours=PAD_AFTER_H)).ceil("h")
    return ws, we


def hold_volume_acft(g, ws, we):
    """Step-function volume over [ws, we]: each reading holds until the next one (the one before ws holds from ws)."""
    before = g[g["datetime"] <= ws]
    q0 = float(before["cfs"].iloc[-1]) if len(before) else 0.0
    inside = g[(g["datetime"] > ws) & (g["datetime"] < we)]
    times = [ws] + list(inside["datetime"]) + [we]
    vals = [q0] + list(inside["cfs"])
    total = 0.0
    for i, q in enumerate(vals):
        total += q * (times[i + 1] - times[i]).total_seconds()
    return total * CFS_TO_ACFT_PER_SEC


def hours_from_start(t):
    return (t - MODEL_START).total_seconds() / 3600.0


def fmt(t):
    return t.strftime("%d %b %H:%M")


def main():
    ap = argparse.ArgumentParser(description="Fix the four storm windows from the real gauge (read-only apart from one small JSON).")
    ap.add_argument("--gauge", default=str(GAUGE), help="gauge workbook (default %s)" % GAUGE)
    ap.add_argument("--out", default=OUT_NAME, help="output JSON name (default %s)" % OUT_NAME)
    ap.add_argument("--overwrite", action="store_true", help="replace an existing output file")
    a = ap.parse_args()

    print("define_storm_windows_110.py  version %s   (reads the real gauge only; looks at no model result)" % VERSION)
    gp, outp = Path(a.gauge), Path(a.out)
    if not gp.exists():
        sys.exit("STOP: cannot find %s. Run this from the lab/ folder." % gp)
    if outp.exists() and not a.overwrite:
        sys.exit("STOP: %s already exists. Nothing was changed. Add --overwrite only if you mean to replace it." % outp)

    print("Reading %s (up to a minute) ..." % gp.name)
    g_all = read_gauge(gp)
    g = g_all[(g_all["datetime"] >= MODEL_START) & (g_all["datetime"] < MODEL_END)].reset_index(drop=True)
    if g.empty:
        sys.exit("STOP: no gauge readings between %s and %s." % (MODEL_START, MODEL_END))
    print("  %d readings in the model window (%s to %s); negative values: %d"
          % (len(g), fmt(MODEL_START), fmt(MODEL_END), int((g["cfs"] < 0).sum())))

    episodes = find_storms(g)
    print("\nStorms found (runs of readings at or above %g cfs, joined if less than %d h apart): %d"
          % (THRESHOLD_CFS, EPISODE_GAP_H, len(episodes)))
    for i, (f, l) in enumerate(episodes, 1):
        print("  %d. first %s   last %s" % (i, f.strftime("%Y-%m-%d %H:%M:%S"), l.strftime("%Y-%m-%d %H:%M:%S")))

    expect = [d for d, _ in STORMS]
    got = [f.strftime("%Y-%m-%d") for f, _ in episodes]
    if got != expect:
        print("\nSTOP: expected exactly four storms starting on %s;\n      found %d starting on %s. Nothing was written."
              % (", ".join(expect), len(got), ", ".join(got) if got else "(none)"))
        sys.exit(2)

    rows = []
    for (date, label), (first, last) in zip(STORMS, episodes):
        rws, rwe = recipe_window(first, last)
        if date == STORMS[0][0]:
            ws, we = STAGE_B_WINDOW
            how = "Stage B window kept"
        else:
            ws, we = rws, rwe
            how = "rule"
        sub = g[(g["datetime"] >= ws) & (g["datetime"] <= we)]
        pk = sub.loc[sub["cfs"].idxmax()]
        rows.append(dict(date=date, label=label, how=how,
                         window_start=ws, window_end=we,
                         start_hour=hours_from_start(ws), end_hour=hours_from_start(we),
                         first_ge1=first, last_ge1=last,
                         recipe_start=rws, recipe_end=rwe,
                         peak_cfs=float(pk["cfs"]), peak_time=pk["datetime"],
                         n_readings=int(len(sub)), n_ge1=int((sub["cfs"] >= THRESHOLD_CFS).sum()),
                         volume_acft=hold_volume_acft(g, ws, we)))

    problems = []
    for r in rows:
        if r["window_start"] < MODEL_START or r["window_end"] > MODEL_END:
            problems.append("%s window runs outside the model run" % r["label"])
    for a_, b_ in zip(rows, rows[1:]):
        if a_["window_end"] > b_["window_start"]:
            problems.append("%s and %s windows overlap" % (a_["label"], b_["label"]))
    if problems:
        print("\nSTOP: " + "; ".join(problems) + ". Nothing was written.")
        sys.exit(2)

    print("\nWINDOWS (hours counted from 1 Aug 2014 00:00; this is how the scorer will cut the long run)")
    print("  %-7s %-17s %-17s %6s %6s %5s %9s %-12s %5s %9s"
          % ("storm", "starts", "ends", "from h", "to h", "len h", "peak cfs", "peak time", "reads", "acre-ft"))
    for r in rows:
        print("  %-7s %-17s %-17s %6.0f %6.0f %5.0f %9.0f %-12s %5d %9.2f"
              % (r["label"], fmt(r["window_start"]), fmt(r["window_end"]), r["start_hour"], r["end_hour"],
                 r["end_hour"] - r["start_hour"], r["peak_cfs"], r["peak_time"].strftime("%d %H:%M"),
                 r["n_readings"], r["volume_acft"]))
    print("  (acre-ft: real flow summed over the window, each reading held until the next)")

    r0 = rows[0]
    print("\n12 Aug check: the rule's own window would be %s to %s (hours %.0f to %.0f), %s the Stage B window "
          "(hours %.0f to %.0f), which is kept."
          % (fmt(r0["recipe_start"]), fmt(r0["recipe_end"]), hours_from_start(r0["recipe_start"]), hours_from_start(r0["recipe_end"]),
             "inside" if (r0["recipe_start"] >= r0["window_start"] and r0["recipe_end"] <= r0["window_end"]) else "NOT inside",
             r0["start_hour"], r0["end_hour"]))

    covered = pd.Series(False, index=g.index)
    for r in rows:
        covered |= (g["datetime"] >= r["window_start"]) & (g["datetime"] <= r["window_end"])
    outside = g[(g["cfs"] > 0) & ~covered]
    print("Readings above 0 cfs that fall outside all four windows: %d%s"
          % (len(outside), "" if len(outside) == 0 else "  <-- look at these"))
    if len(outside):
        print(outside.to_string(index=False))

    doc = dict(version=VERSION, rule=dict(threshold_cfs=THRESHOLD_CFS, pad_before_h=PAD_BEFORE_H, pad_after_h=PAD_AFTER_H,
                                          episode_gap_h=EPISODE_GAP_H, rounding="outward to whole hour",
                                          storm_1="Stage B window kept"),
               model_start=str(MODEL_START), runtime_hours=RUNTIME_HOURS, gauge_file=gp.name,
               created=str(pd.Timestamp.now().floor("s")),
               storms=[dict(date=r["date"], label=r["label"], how=r["how"],
                            window_start=str(r["window_start"]), window_end=str(r["window_end"]),
                            start_hour=r["start_hour"], end_hour=r["end_hour"],
                            first_reading_ge_threshold=str(r["first_ge1"]), last_reading_ge_threshold=str(r["last_ge1"]),
                            real_peak_cfs=r["peak_cfs"], real_peak_time=str(r["peak_time"]),
                            n_readings=r["n_readings"], real_volume_acft_hold_last=round(r["volume_acft"], 3))
                       for r in rows])
    outp.write_text(json.dumps(doc, indent=2))
    print("\nWrote %s (the only file created). Look at the table above; if any window looks wrong, delete that file "
          "and tell me before anything is scored." % outp)


if __name__ == "__main__":
    main()
