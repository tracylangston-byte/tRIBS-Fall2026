#!/usr/bin/env python3
"""
analyze_stage2_screens_sycamore.py -- Stage 2, step 4: what each screened knob does to each target, against tolerance.

Reads screen_soil_summary.csv (W2) and screen_routing_summary.csv (W1), both made by the screen runners. Runs nothing.

For every knob and every target it reports the baseline value, the value at each screened setting, the change from the
baseline, and whether that change is larger than the tolerance ("moves the target"). Then each knob gets a label:

  inert      no target moves by more than its tolerance at any screened value   -> drop it and fix it at baseline
  trade-off  at some screened value one target moves TOWARD its goal and another moves AWAY from its goal
  active     at least one target moves by more than its tolerance, and nothing moves in opposite directions

Targets and tolerances (from the design doc's metrics table; all are proposals):
  volume ratio, model / gauge      goal 1.0    moves if it changes by more than 0.10
  peak error, %                    goal 0      moves if it changes by more than 15 points
  first-arrival error, hours       goal 0      moves if it changes by more than 1 h
  tail volume ratio                goal >= 0.7 moves if it changes by more than 0.10
Toward / away is read from the distance to the goal (for a tail ratio, the shortfall below 0.7). A change in distance
smaller than a quarter of the tolerance counts as neither. KGE, the storm-to-storm balance and the process checks are
reported but are not used for the label (the doc says KGE is a tie-break only, and gives the balance no pass or fail).

usage:  python analyze_stage2_screens_sycamore.py --soil screen_soil_summary.csv --routing screen_routing_summary.csv
Writes: stage2_effects_long.csv (every knob x value x target), stage2_knob_labels.csv (one row per knob).
"""
import argparse
import numpy as np
import pandas as pd

TOL = {"vol": 0.10, "pk": 15.0, "arr": 1.0, "tail": 0.10}
# (short name, column, kind, window)
SOIL_TARGETS = [
    ("flood wk vol", "flood_week__vol_ratio", "vol"), ("flood wk peak", "flood_week__peak_err_pct", "pk"),
    ("flood wk arrival", "flood_week__arrival_err_h", "arr"),
    ("2nd Jan vol", "jan_second_storm__vol_ratio", "vol"), ("2nd Jan peak", "jan_second_storm__peak_err_pct", "pk"),
    ("Feb18-19 vol", "feb18_19__vol_ratio", "vol"), ("Feb18-19 peak", "feb18_19__peak_err_pct", "pk"),
    ("Feb18-19 arrival", "feb18_19__arrival_err_h", "arr"),
    ("Feb27-28 vol", "feb27_28__vol_ratio", "vol"), ("Feb27-28 peak", "feb27_28__peak_err_pct", "pk"),
    ("Feb27-28 arrival", "feb27_28__arrival_err_h", "arr"),
    ("tail Jan26-Feb2", "tail_jan26_feb2__tail_vol_ratio", "tail"), ("tail Mar3-11", "tail_mar3_mar11__tail_vol_ratio", "tail"),
    ("whole window vol", "span__vol_ratio", "vol"),
]
ROUTING_TARGETS = [
    ("flood wk vol", "flood_week__vol_ratio", "vol"), ("flood wk peak", "flood_week__peak_err_pct", "pk"),
    ("flood wk arrival", "flood_week__arrival_err_h", "arr"),
    ("2nd Jan vol", "jan_second_storm__vol_ratio", "vol"), ("2nd Jan peak", "jan_second_storm__peak_err_pct", "pk"),
]
# columns shown but not used for the label
REPORT_ONLY = {"kge flood wk": "flood_week__kge2012", "balance Feb27-28/Feb18-19": "balance__feb27_28_over_feb18_19__model"}


def gap(x, kind):
    if not np.isfinite(x):
        return np.nan
    if kind == "vol":
        return abs(x - 1.0)
    if kind == "tail":
        return max(0.0, 0.7 - x)
    return abs(x)


def analyze(df, targets, window):
    rows, labels = [], []
    for knob, g in df.groupby("knob", sort=False):
        g = g.sort_values("value")
        base = g[g["is_baseline"] == True]
        if base.empty:
            # a knob whose baseline value is not a row (m): take the baseline row from any other knob (same run)
            base = df[df["is_baseline"] == True].iloc[[0]]
        b = base.iloc[0]
        moved_any, toward, away, per_value = set(), {}, {}, {}
        for _, r in g.iterrows():
            if r["is_baseline"] == True:
                continue
            v = r["value"]
            for name, col, kind in targets:
                x, x0 = r[col], b[col]
                d = x - x0 if np.isfinite(x) and np.isfinite(x0) else np.nan
                dgap = gap(x, kind) - gap(x0, kind) if np.isfinite(d) else np.nan
                moves = bool(np.isfinite(d) and abs(d) > TOL[kind])
                small = TOL[kind] * 0.25
                direction = ""
                if moves:
                    direction = "toward" if dgap < -small else ("away" if dgap > small else "neither")
                rows.append(dict(window=window, knob=knob, value=v, target=name, baseline=x0, at_value=x, change=d,
                                 tolerance=TOL[kind], moves=moves, direction=direction, run_id=r["run_id"]))
                if moves:
                    moved_any.add(name)
                    per_value.setdefault(v, {"toward": [], "away": []})
                    if direction in ("toward", "away"):
                        per_value[v][direction].append(name)
        both = {v: d for v, d in per_value.items() if d["toward"] and d["away"]}
        if not moved_any:
            label = "inert"
        elif both:
            label = "trade-off"
        else:
            label = "active"
        labels.append(dict(window=window, knob=knob, label=label, n_targets_moved=len(moved_any),
                           targets_moved="; ".join(sorted(moved_any)),
                           trade_off_at=" | ".join(f"{v:g}: toward [{', '.join(d['toward'])}] / away [{', '.join(d['away'])}]"
                                                   for v, d in sorted(both.items()))))
    return pd.DataFrame(rows), pd.DataFrame(labels)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--soil", required=True)
    ap.add_argument("--routing", required=True)
    a = ap.parse_args()
    soil, rout = pd.read_csv(a.soil), pd.read_csv(a.routing)
    e1, l1 = analyze(soil, SOIL_TARGETS, "W2")
    e2, l2 = analyze(rout, ROUTING_TARGETS, "W1")
    eff, lab = pd.concat([e2, e1]), pd.concat([l2, l1])
    eff.to_csv("stage2_effects_long.csv", index=False)
    lab.to_csv("stage2_knob_labels.csv", index=False)
    pd.set_option("display.width", 250, "display.max_colwidth", 200)
    print(lab.to_string(index=False))


if __name__ == "__main__":
    main()
