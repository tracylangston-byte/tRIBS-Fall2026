"""
rescore_cc_truths_110.py
========================
Series 110 -- RE-SCORE the existing Stage 2 (cc ON) and control (cc OFF) runs
against each CANDIDATE cc-ON synthetic truth. No tRIBS, no new simulations.

WHAT THIS DOES (plain language)
--------------------------------
Stage 2 ran 250 joint (Ks, f, cc) simulations with channel loss ON, and the
control ran the same 250 points with channel loss OFF. Every one of those 500
runs was scored against the cc-OFF truth. But a simulation does not depend on
which truth it is later compared to -- only the SCORE does. So the 500
simulated hydrographs we already have can be compared to a different truth
(one of the cc-ON candidates) in seconds, with no new model runs.

Why it is exact and not an approximation: each run's *_compare_obs_sim.csv
already holds that run's simulated series, resampled to 5 minutes and clipped
to the event window by the scorer. This script keeps that "Simulated" column
untouched, replaces "Observed" with a candidate truth passed through the very
same reading and resampling calls the scorer uses, and recomputes the same
metrics with the same formulas.

SELF-CHECKS ("gates") -- run before anything is re-scored
-----------------------------------------------------------
The script does not trust itself; it checks against your existing outputs.
  GATE A  Feed the current cc-OFF truth (the one in synth_truth/) through this
          script's reading/resampling code. The result must equal the
          "Observed" column stored in every one of the 500 compare CSVs.
          Catches: wrong truth file in synth_truth/, a different resampling.
  GATE B  Recompute every metric from each compare CSV and compare it with
          the metric stored in the Stage 2 / control results CSV (KGE_2012,
          PBIAS, beta, gamma, r, NSE, RMSE, peak and volume errors, and the
          scorer's phase metrics). Catches: any formula or file mismatch, a
          compare CSV that was overwritten by a later run.
  GATE C  Each candidate truth file's checksum must match the one recorded
          when it was copied (PROVENANCE_cc_truth_candidates_110.json), and
          its 5-minute grid must line up with the cc-OFF truth's grid.
If a gate fails the script stops and says why; nothing is re-scored.
(--drop_bad_runs excludes individual runs that fail A/B and continues.)

WHAT YOU GET   (calibration_work/03_comparisons/summary_tables/rescored_cc_truths_110/)
----------------------------------------------------------------------------------------
  rescored_truth_summary_110.csv   one row per (truth, series): best KGE_2012,
                                   how many runs are feasible, where the near-best
                                   runs sit in cc, how close the nearest sample is
                                   to the truth point
  rescored_long_110.csv            every run x every truth: all metrics
  stage2_ccON_vs_truth_cc<V>.csv   per truth: the 250 cc-ON runs, original
  control_ccOFF_vs_truth_cc<V>.csv   columns kept, metrics replaced. These two
                                   files have exactly the layout analyze_cc_pareto_110.py
                                   reads, so:
                                   python analyze_cc_pareto_110.py \\
                                     --stage2_csv <stage2_ccON_vs_truth_cc122.csv> \\
                                     --control_csv <control_ccOFF_vs_truth_cc122.csv> \\
                                     --out_dir <somewhere new>
  PROVENANCE_rescored_cc_truths_110.json   checksums, versions, gate results

SAFETY
-------
  - Read-only on every existing file. Writes only inside its own output folder
    (refuses an output folder that is, or is inside, any input location).
  - Never imports the builder, never touches current_run_config.json, never
    starts tRIBS, so it cannot collide with a sweep's shared state. (Do not run
    it WHILE a sweep is still rewriting the results CSVs it reads.)
  - Run membership comes only from the run_id rows in the two results CSVs; it
    never globs the compare folder (extra _s2 files from earlier attempts exist).

HOW TO READ THE RESULT (rules of thumb, mine)
-----------------------------------------------
  - "ON best KGE" is the best score any cc-ON run reached against that truth.
    "OFF best KGE" is the best a model WITHOUT channel loss reached against the
    same truth by moving Ks and f alone. A small gap means Ks/f can imitate this
    truth well enough that the data cannot easily say "channel loss is present".
  - "near-best cc range" is the span of cc among cc-ON runs within --kge_tol of
    the best. A narrow range containing the true cc is what identifiable looks
    like; a range covering most of the box is what not-identifiable looks like.
  - "nearest sample" is the distance, in units of the sampled box (0 = same
    point, 1 = opposite corners of one axis), from the truth point to the
    closest of the 250 samples. Large distance = the sweep never looked near the
    truth, so a poor recovery there says little. (n within 0.15 is the count of
    samples inside that radius; a uniform sweep would put about 3-4 there.)
  - One storm, 250 points in 3 dimensions, routing pinned at truth, noise-free
    truth: this is the best case, not a realistic one. Feasible = |PBIAS| < 2%.

USAGE (run from the lab/ directory, like every other script)
-------------------------------------------------------------
    python rescore_cc_truths_110.py --check_only      # gates only; writes nothing
    python rescore_cc_truths_110.py                   # gates, then all candidates
    python rescore_cc_truths_110.py --cc_values 122 157 425
    python rescore_cc_truths_110.py --kge_tol 0.005
    python rescore_cc_truths_110.py --drop_bad_runs   # skip runs that fail a gate
"""

import argparse
import hashlib
import json
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

# ------------------------------------------------------------------
# Constants (must agree with the scorer and the Stage 2 / control scripts)
# ------------------------------------------------------------------
ORIGIN       = pd.Timestamp("2014-08-01")      # .qout time is hours from here (scorer)
QOUT_NAMES   = ["Time_hr", "Qstrm_m3s", "Hlev_m"]

STAGE2_NAME  = "lhs_results_joint_Ks_f_cc_110.csv"
CONTROL_NAME = "lhs_results_joint_Ks_f_cc_CONTROL_110.csv"
CAND_DIRNAME = "synth_truth_cc_candidates"
CAND_SUMMARY = "cc_truth_candidates_110.csv"
CAND_PROV    = "PROVENANCE_cc_truth_candidates_110.json"
OUT_DIRNAME  = "rescored_cc_truths_110"

SER_ON, SER_OFF = "cc ON", "cc OFF"

# Stage 2 sampled box (Ks linear; f and cc log) -- from LHS_PARAMS.
BOX = {
    "Ks_mult":                  (4.0,   10.5,   "lin"),
    "f_RS_abs":                 (0.004, 0.030,  "log"),
    "channelconductivity_mmhr": (30.0,  1000.0, "log"),
}
TRUTH_KS, TRUTH_F = 7.0, 0.012

FEASIBLE_PBIAS_PCT = 2.0     # same definition as the Stage 2 analysis
NEAR_RADIUS        = 0.15    # box-fraction radius for "samples near the truth"
DEFAULT_EXPECT_N   = 250
DEFAULT_KGE_TOL    = 0.01

RTOL, ATOL = 1e-7, 1e-9      # gate B: stored vs recomputed metrics
OBS_RTOL, OBS_ATOL = 1e-9, 1e-12   # gate A: observed series

SCALAR_METRICS = [
    "obs_peak_m3s", "sim_peak_m3s", "peak_error_m3s", "peak_error_pct",
    "peak_timing_error_hr", "obs_volume_m3", "sim_volume_m3", "volume_error_pct",
    "rmse_m3s", "nse", "pbias_pct", "kge", "kge_r", "kge_alpha", "kge_beta",
    "kge_gamma", "kge_2012",
]
TIME_METRICS  = ["obs_peak_time", "sim_peak_time"]       # replaced, not gated
PHASE_METRICS = [
    "threshold_m3s", "first_arrival_error_min", "rising_limb_steepness_ratio",
    "time_to_peak_from_exc_min", "duration_above_thresh_error_min",
    "recession_rate_ratio",
]
REQUIRED_COLS = ["run_id", "Ks_mult", "f_RS_abs", "channelconductivity_mmhr",
                 "optpercolation"]


# ------------------------------------------------------------------
# Small helpers
# ------------------------------------------------------------------
def md5_of(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def label_for_cc(cc):
    return f"{cc:g}".replace(".", "p")


def _num(x):
    try:
        v = float(x)
        return v if np.isfinite(v) else np.nan
    except (TypeError, ValueError):
        return np.nan


def import_phase_fn(script_dir):
    """The scorer's own phase-metric function (so those numbers are the
    scorer's, not a copy). Its module-level imports are only numpy/pandas;
    pytRIBS is imported inside run_and_score(), which is never called here."""
    sys.path.insert(0, str(script_dir))
    try:
        import run_sensitivity_single_interp_tribs6 as scorer
        return scorer._compute_phase_metrics
    except Exception as e:
        print(f"  NOTE: could not import the scorer's phase-metric function ({e}).\n"
              f"        Phase metrics will be left blank in the outputs and skipped in gate B.")
        return None


# ------------------------------------------------------------------
# Reading a truth the way the scorer does
# ------------------------------------------------------------------
def read_truth_5min(path):
    """Read a .qout and resample to the 5-minute grid, using the same calls as
    run_sensitivity_single_interp_tribs6.py (format sniffing, time origin,
    resample('5min').interpolate(method='time')). Returns a Series."""
    with open(path) as fh:
        fh.readline()                       # header
        first_data = fh.readline()
    if "," in first_data:                   # tRIBS 6.0.0: comma-delimited
        raw = pd.read_csv(path, sep=",", header=0, names=QOUT_NAMES)
    else:                                   # old tRIBS: whitespace-delimited
        raw = pd.read_csv(path, sep=r"\s+", skiprows=1, names=QOUT_NAMES)
    raw["datetime"] = pd.to_datetime(raw["Time_hr"] * 3600, unit="s", origin=ORIGIN)
    raw = raw.set_index("datetime")
    return raw["Qstrm_m3s"].resample("5min").interpolate(method="time")


def load_compare(path):
    """A scorer-written compare CSV: first column is the time index (its header
    is blank), then Observed and Simulated."""
    df = pd.read_csv(path, index_col=0, parse_dates=True)
    missing = {"Observed", "Simulated"} - set(df.columns)
    if missing:
        raise ValueError(f"missing column(s) {sorted(missing)}")
    df.index = pd.DatetimeIndex(df.index)
    if len(df) < 3:
        raise ValueError(f"only {len(df)} rows")
    return df[["Observed", "Simulated"]]


# ------------------------------------------------------------------
# Metrics -- the scorer's formulas, copied (gate B proves they still agree)
# ------------------------------------------------------------------
def compute_metrics(obs, sim, phase_fn):
    """obs, sim: pandas Series on the same 5-minute DatetimeIndex (the event
    window). Same formulas as run_and_score() in the scorer."""
    with np.errstate(all="ignore"):
        obs_peak, sim_peak   = obs.max(), sim.max()
        obs_tpeak, sim_tpeak = obs.idxmax(), sim.idxmax()

        dt_seconds = (obs.index[1] - obs.index[0]).total_seconds()
        obs_vol    = obs.sum() * dt_seconds
        sim_vol    = sim.sum() * dt_seconds

        peak_err_m3s = sim_peak - obs_peak
        rmse  = np.sqrt(np.mean((sim - obs) ** 2))
        nse   = 1 - (np.sum((sim - obs) ** 2) / np.sum((obs - obs.mean()) ** 2))
        pbias = 100 * (np.sum(sim - obs) / np.sum(obs))

        r     = np.corrcoef(sim, obs)[0, 1]
        alpha = np.std(sim) / np.std(obs)
        beta  = np.mean(sim) / np.mean(obs)
        kge   = 1 - np.sqrt((r - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2)
        gamma = alpha / beta
        kge_2012 = 1 - np.sqrt((r - 1) ** 2 + (gamma - 1) ** 2 + (beta - 1) ** 2)

        out = {
            "obs_peak_m3s":         obs_peak,
            "sim_peak_m3s":         sim_peak,
            "peak_error_m3s":       peak_err_m3s,
            "peak_error_pct":       (peak_err_m3s / obs_peak) * 100,
            "obs_peak_time":        str(obs_tpeak),
            "sim_peak_time":        str(sim_tpeak),
            "peak_timing_error_hr": (sim_tpeak - obs_tpeak).total_seconds() / 3600,
            "obs_volume_m3":        obs_vol,
            "sim_volume_m3":        sim_vol,
            "volume_error_pct":     ((sim_vol - obs_vol) / obs_vol) * 100,
            "rmse_m3s":             rmse,
            "nse":                  nse,
            "pbias_pct":            pbias,
            "kge":                  kge,
            "kge_r":                r,
            "kge_alpha":            alpha,
            "kge_beta":             beta,
            "kge_gamma":            gamma,
            "kge_2012":             kge_2012,
        }
        if phase_fn is not None:
            try:
                out.update(phase_fn(obs, sim))
            except Exception:
                out.update({k: np.nan for k in PHASE_METRICS})
        else:
            out.update({k: np.nan for k in PHASE_METRICS})
    return {k: (float(v) if k not in TIME_METRICS else v) for k, v in out.items()}


# ------------------------------------------------------------------
# Loading inputs
# ------------------------------------------------------------------
def load_results(path, label, expect_perc, expect_n):
    if not path.exists():
        sys.exit(f"[{label}] results file not found: {path}")
    df = pd.read_csv(path)
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        sys.exit(f"[{label}] {path.name} is missing required columns: {missing}")
    if df["run_id"].duplicated().any():
        dup = df.loc[df["run_id"].duplicated(), "run_id"].tolist()[:3]
        sys.exit(f"[{label}] duplicate run_id values in {path.name} (e.g. {dup}); refusing.")
    for c in ["Ks_mult", "f_RS_abs", "channelconductivity_mmhr", "optpercolation"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    vals = sorted(df["optpercolation"].dropna().unique().tolist())
    if vals != [expect_perc]:
        sys.exit(f"[{label}] {path.name}: optpercolation values are {vals}, expected "
                 f"[{expect_perc}]. Wrong file for this series?")
    note = "" if len(df) == expect_n else f"   <-- expected {expect_n}"
    print(f"  [{label}] {path.name}: {len(df)} runs{note}")
    return df.reset_index(drop=True)


def load_all_compares(results, csv_dir):
    cache, problems = {}, {}
    for rid in results["run_id"]:
        p = csv_dir / f"{rid}_compare_obs_sim.csv"
        if not p.exists():
            problems[rid] = "compare CSV missing"
            continue
        try:
            cache[rid] = load_compare(p)
        except Exception as e:
            problems[rid] = f"compare CSV unreadable: {str(e)[:100]}"
    return cache, problems


def close(a, b, rtol, atol):
    return bool(np.isclose(a, b, rtol=rtol, atol=atol, equal_nan=True))


# ------------------------------------------------------------------
# Gates
# ------------------------------------------------------------------
def gate_a(results, cache, obs_ccoff):
    """Observed column in each compare CSV == the cc-OFF truth, as read here."""
    bad, worst = {}, 0.0
    for rid in results["run_id"]:
        if rid not in cache:
            continue
        df = cache[rid]
        ref = obs_ccoff.reindex(df.index).to_numpy(dtype=float)
        got = df["Observed"].to_numpy(dtype=float)
        if np.isnan(ref).any():
            bad[rid] = "cc-OFF truth has no value at some compare-CSV timestamps"
            continue
        d = float(np.max(np.abs(ref - got)))
        worst = max(worst, d)
        if not np.allclose(ref, got, rtol=OBS_RTOL, atol=OBS_ATOL):
            bad[rid] = f"Observed differs from the cc-OFF truth (max abs diff {d:.3g} m3/s)"
    return bad, worst


def gate_b(results, cache, phase_fn):
    """Metrics recomputed from the compare CSV == metrics stored in the results CSV."""
    bad, worst = {}, {}
    check_cols = SCALAR_METRICS + (PHASE_METRICS if phase_fn is not None else [])
    check_cols = [c for c in check_cols if c in results.columns]
    for _, row in results.iterrows():
        rid = row["run_id"]
        if rid not in cache:
            continue
        df = cache[rid]
        m = compute_metrics(df["Observed"], df["Simulated"], phase_fn)
        off = []
        for c in check_cols:
            stored = _num(row[c])
            new    = _num(m[c])
            if not close(new, stored, RTOL, ATOL):
                off.append(f"{c} (stored {stored:.6g}, recomputed {new:.6g})")
            elif np.isfinite(stored) and np.isfinite(new):
                worst[c] = max(worst.get(c, 0.0), abs(new - stored))
        if off:
            bad[rid] = "stored metric differs: " + "; ".join(off[:3])
    return bad, worst


# ------------------------------------------------------------------
# Candidate truths
# ------------------------------------------------------------------
def find_candidates(cand_dir, summary_dir):
    recorded, prov_files = {}, {}
    prov_path = cand_dir / CAND_PROV
    if prov_path.exists():
        try:
            prov_files = json.loads(prov_path.read_text()).get("files", {})
        except Exception:
            prov_files = {}
    rows = {}
    csv_path = summary_dir / CAND_SUMMARY
    if csv_path.exists():
        df = pd.read_csv(csv_path)
        for _, r in df.iterrows():
            name = str(r.get("candidate_qout", ""))
            if name:
                fl = "" if pd.isna(r.get("flags", "")) else str(r.get("flags", ""))
                rows[name] = {"name": name, "cc": float(r["cc_mmhr"]),
                              "md5": str(r.get("candidate_md5", "")), "flags": fl}
    for name, info in prov_files.items():
        rows.setdefault(name, {"name": name, "cc": float(info["cc_mmhr"]),
                               "md5": info.get("md5", ""), "flags": info.get("flags", "")})
        if not rows[name]["md5"]:
            rows[name]["md5"] = info.get("md5", "")
    for p in sorted(cand_dir.glob("*_cctruth_Outlet.qout")):
        if p.name not in rows:
            m = re.search(r"_cc([0-9p]+)_cctruth", p.name)
            if m:
                rows[p.name] = {"name": p.name, "cc": float(m.group(1).replace("p", ".")),
                                "md5": "", "flags": ""}
    return sorted(rows.values(), key=lambda d: d["cc"])


# ------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------
def _unit(value, lo, hi, scale):
    if scale == "log":
        return (np.log10(value) - np.log10(lo)) / (np.log10(hi) - np.log10(lo))
    return (value - lo) / (hi - lo)


def summarize(m, series, truth_cc, kge_tol):
    """One summary row for one (truth, series). m has one row per run."""
    out = {"truth_cc_mmhr": truth_cc, "series": series, "n_scored": int(len(m))}
    ok = m[np.isfinite(m["kge_2012"])]
    out["n_valid_kge"] = int(len(ok))
    if ok.empty:
        return out
    best = ok.loc[ok["kge_2012"].idxmax()]
    out.update({
        "best_kge_2012": float(best["kge_2012"]),
        "best_run_id":   best["run_id"],
        "best_Ks":       float(best["Ks_mult"]),
        "best_f":        float(best["f_RS_abs"]),
        "best_pbias_pct": float(best["pbias_pct"]),
        "best_dKs_from_truth":    float(best["Ks_mult"] - TRUTH_KS),
        "best_dlog10_f_from_truth": float(np.log10(best["f_RS_abs"] / TRUTH_F)),
    })
    feas = ok[np.abs(ok["pbias_pct"]) < FEASIBLE_PBIAS_PCT]
    out["n_feasible"]   = int(len(feas))
    out["pct_feasible"] = 100.0 * len(feas) / len(ok)
    near = ok[ok["kge_2012"] >= best["kge_2012"] - kge_tol]
    out["n_near_best"] = int(len(near))
    out["kge_tol"]     = kge_tol

    if series == SER_ON:
        cc = near["channelconductivity_mmhr"].astype(float)
        lo, hi = float(cc.min()), float(cc.max())
        span = np.log10(hi / lo) / np.log10(BOX["channelconductivity_mmhr"][1]
                                            / BOX["channelconductivity_mmhr"][0])
        out.update({
            "near_cc_min": lo, "near_cc_median": float(cc.median()), "near_cc_max": hi,
            "near_cc_logspan_frac_of_box": float(span),
            "truth_cc_inside_near_range": bool(lo <= truth_cc <= hi),
            "best_cc": float(best["channelconductivity_mmhr"]),
            "best_dlog10_cc_from_truth": float(np.log10(best["channelconductivity_mmhr"] / truth_cc)),
        })
        pts = np.column_stack([
            (ok["Ks_mult"].to_numpy(float) - TRUTH_KS) / (BOX["Ks_mult"][1] - BOX["Ks_mult"][0]),
            _unit(ok["f_RS_abs"].to_numpy(float), *BOX["f_RS_abs"][:2], "log")
            - _unit(TRUTH_F, *BOX["f_RS_abs"][:2], "log"),
            _unit(ok["channelconductivity_mmhr"].to_numpy(float),
                  *BOX["channelconductivity_mmhr"][:2], "log")
            - _unit(truth_cc, *BOX["channelconductivity_mmhr"][:2], "log"),
        ])
        dist = np.sqrt((pts ** 2).sum(axis=1))
        i = int(np.argmin(dist))
        out.update({
            "nearest_sample_dist":   float(dist[i]),
            "nearest_sample_run_id": ok.iloc[i]["run_id"],
            "nearest_sample_kge":    float(ok.iloc[i]["kge_2012"]),
            "n_samples_within_0p15": int((dist <= NEAR_RADIUS).sum()),
        })
    return out


# ------------------------------------------------------------------
# Output helpers
# ------------------------------------------------------------------
def write_pair(res_df, mdf, cc, path, phase_available):
    """Original results columns kept; every scorer metric column replaced by
    the re-scored value. Same layout analyze_cc_pareto_110.py reads."""
    out = res_df[res_df["run_id"].isin(mdf["run_id"])].copy()
    mi = mdf.set_index("run_id")
    replace = SCALAR_METRICS + TIME_METRICS + PHASE_METRICS
    for col in replace:
        if col in out.columns and col in mi.columns:
            if col in PHASE_METRICS and not phase_available:
                out[col] = np.nan
            else:
                out[col] = out["run_id"].map(mi[col])
    out["rescored_vs_truth_cc_mmhr"] = cc
    out.to_csv(path, index=False)


def fmt(v, spec):
    return "   n/a" if (v is None or (isinstance(v, float) and not np.isfinite(v))) else format(v, spec)


def report(summ, skipped):
    on  = summ[summ["series"] == SER_ON].set_index("truth_cc_mmhr")
    off = summ[summ["series"] == SER_OFF].set_index("truth_cc_mmhr")
    print(f"\n{'=' * 118}")
    print("RE-SCORED AGAINST EACH CANDIDATE cc-ON TRUTH   (feasible = |PBIAS| < "
          f"{FEASIBLE_PBIAS_PCT:g}%;  near-best = within kge_tol of the best KGE_2012)")
    print(f"{'=' * 118}")
    print(f"{'truth cc':>8} | {'ON best':>7} {'OFF best':>8} {'gap':>6} | {'feasible':>9} {'ON/OFF':>7} | "
          f"{'near-best cc range (ON)':>26} {'span':>5} {'in?':>4} | "
          f"{'nearest':>7} {'n<.15':>5}  flags")
    for cc in sorted(on.index):
        a = on.loc[cc]
        b = off.loc[cc] if cc in off.index else None
        gap = (a["best_kge_2012"] - b["best_kge_2012"]) if b is not None else np.nan
        rng = (f"{a['near_cc_min']:.0f}-{a['near_cc_max']:.0f} (med {a['near_cc_median']:.0f})"
               if "near_cc_min" in a and np.isfinite(a.get("near_cc_min", np.nan)) else "n/a")
        feas_on  = int(a["n_feasible"]) if np.isfinite(a.get("n_feasible", np.nan)) else -1
        feas_off = int(b["n_feasible"]) if b is not None and np.isfinite(b.get("n_feasible", np.nan)) else -1
        print(f"{cc:>8.4g} | {fmt(a['best_kge_2012'], '7.4f')} "
              f"{fmt(b['best_kge_2012'], '8.4f') if b is not None else '     n/a'} "
              f"{fmt(gap, '6.3f')} | {feas_on:>4d}/{feas_off:<4d} "
              f"{'' if feas_off <= 0 else format(feas_on / feas_off, '7.2f'):>7} | "
              f"{rng:>26} {fmt(a.get('near_cc_logspan_frac_of_box', np.nan), '5.2f')} "
              f"{'yes' if a.get('truth_cc_inside_near_range') else 'no':>4} | "
              f"{fmt(a.get('nearest_sample_dist', np.nan), '7.3f')} "
              f"{int(a.get('n_samples_within_0p15', 0)):>5d}  {a.get('candidate_flags', '')}")
    print("\n  ON best / OFF best : best KGE_2012 of the cc-ON runs / of the cc-OFF control runs, both vs this truth")
    print("  gap                : ON best minus OFF best (small = Ks and f alone can imitate this truth)")
    print("  feasible ON/OFF    : number of runs with |PBIAS| < 2% in each series (same 250 points)")
    print("  span / in?         : share of the cc box (log scale) covered by the near-best runs / does it contain the true cc")
    print("  nearest / n<.15    : distance from the truth point to the closest sample, in box units / samples within 0.15")
    print("  Caveat: one storm, 250 points in 3-D, routing pinned at truth, noise-free truth -- a best case.")
    if skipped:
        print("\n  SKIPPED TRUTHS:")
        for name, why in skipped:
            print(f"    {name}: {why}")


# ------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Series 110 -- re-score the Stage 2 and control runs against "
                    "each candidate cc-ON truth (no tRIBS).")
    ap.add_argument("--check_only", action="store_true",
                    help="Run the gates and stop; writes nothing.")
    ap.add_argument("--cc_values", type=float, nargs="+", default=None,
                    help="Only these candidate cc values (mm/hr). Default: all candidates found.")
    ap.add_argument("--kge_tol", type=float, default=DEFAULT_KGE_TOL,
                    help=f"Near-best tolerance in KGE_2012 units (default {DEFAULT_KGE_TOL})")
    ap.add_argument("--expect_n", type=int, default=DEFAULT_EXPECT_N,
                    help=f"Expected runs per series (default {DEFAULT_EXPECT_N}); a different count only warns")
    ap.add_argument("--drop_bad_runs", action="store_true",
                    help="Exclude runs that fail a gate and continue, instead of stopping.")
    ap.add_argument("--no_pairs", action="store_true",
                    help="Skip the per-truth stage2/control CSV pairs.")
    ap.add_argument("--stage2_csv", type=Path, default=None)
    ap.add_argument("--control_csv", type=Path, default=None)
    ap.add_argument("--out_dir", type=Path, default=None)
    args = ap.parse_args()

    script_dir   = Path.cwd()
    calib_dir    = script_dir.parent / "calibration_work"
    summary_dir  = calib_dir / "03_comparisons" / "summary_tables"
    csv_dir      = calib_dir / "03_comparisons" / "csv_exports"
    synth_dir    = calib_dir / "synth_truth"
    cand_dir     = calib_dir / CAND_DIRNAME
    out_dir      = args.out_dir or (summary_dir / OUT_DIRNAME)
    p_on         = args.stage2_csv or (summary_dir / STAGE2_NAME)
    p_off        = args.control_csv or (summary_dir / CONTROL_NAME)

    print(f"\n{'=' * 78}\nSeries 110 -- re-score Stage 2 / control against candidate cc-ON truths\n{'=' * 78}")

    # ---- inputs / safety -------------------------------------------------
    qouts = list(synth_dir.glob("*.qout")) if synth_dir.exists() else []
    if len(qouts) != 1:
        sys.exit(f"Expected exactly one *.qout in {synth_dir} (the cc-OFF truth the 500 runs were "
                 f"scored against), found {len(qouts)}: {[q.name for q in qouts]}")
    ccoff_path = qouts[0]
    print(f"cc-OFF truth (read-only): {ccoff_path.name}")

    out_res = out_dir.resolve()
    for guarded in (synth_dir, cand_dir, csv_dir):
        g = guarded.resolve()
        if out_res == g or g in out_res.parents or out_res in g.parents:
            sys.exit(f"Output folder {out_dir} overlaps an input folder ({guarded}); refusing.")
    for pth in (p_on, p_off):
        if out_res in pth.resolve().parents or out_res == pth.resolve():
            sys.exit(f"Output folder {out_dir} contains an input file ({pth}); refusing.")

    print("\nLoading results tables:")
    res_on  = load_results(p_on,  SER_ON,  1, args.expect_n)
    res_off = load_results(p_off, SER_OFF, 0, args.expect_n)

    phase_fn = import_phase_fn(script_dir)
    obs_ccoff = read_truth_5min(ccoff_path)

    print("\nLoading compare CSVs (run membership = run_id rows above):")
    cache_on,  prob_on  = load_all_compares(res_on,  csv_dir)
    cache_off, prob_off = load_all_compares(res_off, csv_dir)
    print(f"  {SER_ON}: {len(cache_on)}/{len(res_on)} loaded;  "
          f"{SER_OFF}: {len(cache_off)}/{len(res_off)} loaded")

    # ---- gates A and B ---------------------------------------------------
    bad = {}
    for rid, why in {**prob_on, **prob_off}.items():
        bad[rid] = why
    cache = {**cache_on, **cache_off}
    allres = pd.concat([res_on.assign(_series=SER_ON), res_off.assign(_series=SER_OFF)],
                       ignore_index=True)
    badA, worstA = gate_a(allres, cache, obs_ccoff)
    badB, worstB = gate_b(allres, cache, phase_fn)
    for d in (badA, badB):
        for rid, why in d.items():
            bad.setdefault(rid, why)

    n_checked = len(cache)
    print(f"\nGATE A  observed series : {n_checked - len(badA)}/{n_checked} compare CSVs match the cc-OFF "
          f"truth as read here (worst abs diff {worstA:.2e} m3/s)")
    if worstB:
        wk = max(worstB, key=worstB.get)
        print(f"GATE B  stored metrics  : {n_checked - len(badB)}/{n_checked} runs reproduce their stored "
              f"metrics (worst abs deviation {worstB[wk]:.2e} in {wk}; "
              f"{'phase metrics included' if phase_fn else 'phase metrics skipped'})")
    else:
        print(f"GATE B  stored metrics  : {n_checked - len(badB)}/{n_checked}")

    if bad:
        print(f"\n  {len(bad)} run(s) failed a gate or could not be loaded. First few:")
        for rid, why in list(bad.items())[:8]:
            print(f"    {rid}: {why}")
        if not args.drop_bad_runs:
            sys.exit("\nSTOPPED: nothing was re-scored. Fix the cause (or re-run with "
                     "--drop_bad_runs to exclude these runs and continue).")
        print("  --drop_bad_runs set: these runs are excluded.")
    res_on  = res_on[~res_on["run_id"].isin(bad)].reset_index(drop=True)
    res_off = res_off[~res_off["run_id"].isin(bad)].reset_index(drop=True)
    if res_on.empty or res_off.empty:
        sys.exit("No runs left in one of the series after the gates; stopping.")

    # ---- candidates (gate C) ---------------------------------------------
    cands = find_candidates(cand_dir, summary_dir)
    if args.cc_values:
        keep = []
        available = ", ".join("%g" % c["cc"] for c in cands) or "none"
        for v in args.cc_values:
            hit = [c for c in cands if np.isclose(c["cc"], v, rtol=1e-3)]
            if not hit:
                sys.exit(f"--cc_values: no candidate truth near cc={v:g}. "
                         f"Available: {available}")
            keep.append(hit[0])
        cands = sorted({c["name"]: c for c in keep}.values(), key=lambda d: d["cc"])
    if not cands:
        sys.exit(f"No candidate truths found in {cand_dir}")

    usable, skipped, cand_obs = [], [], {}
    for c in cands:
        p = cand_dir / c["name"]
        if not p.exists():
            skipped.append((c["name"], "file missing")); continue
        if c["md5"] and md5_of(p) != c["md5"]:
            skipped.append((c["name"], "checksum differs from the one recorded at copy time (GATE C)")); continue
        if not c["md5"]:
            print(f"  NOTE: no recorded checksum for {c['name']}; provenance not verified.")
        try:
            s = read_truth_5min(p)
        except Exception as e:
            skipped.append((c["name"], f"unreadable: {str(e)[:80]}")); continue
        if not (s.index.equals(obs_ccoff.index) and (s.isna() == obs_ccoff.isna()).all()):
            skipped.append((c["name"], "5-minute grid differs from the cc-OFF truth's (GATE C)")); continue
        usable.append(c); cand_obs[c["name"]] = s
    print(f"GATE C  candidate truths: {len(usable)}/{len(cands)} verified "
          f"(checksum and 5-minute grid)")
    for name, why in skipped:
        print(f"    skipped {name}: {why}")
    if not usable:
        sys.exit("No usable candidate truths.")

    if args.check_only:
        print("\n--check_only: gates finished, nothing written.")
        return

    # ---- re-score ----------------------------------------------------------
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"\nRe-scoring {len(res_on)} + {len(res_off)} runs x {len(usable)} truth(s)...")
    long_rows, summ_rows, truth_failed = [], [], []
    for c in usable:
        cc, obs_new = c["cc"], cand_obs[c["name"]]
        per_series = {}
        for series, res in ((SER_ON, res_on), (SER_OFF, res_off)):
            rows = []
            for _, r in res.iterrows():
                df = cache[r["run_id"]]
                o = obs_new.reindex(df.index)
                if o.isna().any():
                    truth_failed.append((c["name"], r["run_id"], "no truth value at some timestamps"))
                    continue
                m = compute_metrics(o, df["Simulated"], phase_fn)
                rows.append({"truth_cc_mmhr": cc, "series": series, "run_id": r["run_id"],
                             "Ks_mult": r["Ks_mult"], "f_RS_abs": r["f_RS_abs"],
                             "channelconductivity_mmhr": r["channelconductivity_mmhr"],
                             "optpercolation": r["optpercolation"], **m})
            mdf = pd.DataFrame(rows)
            per_series[series] = mdf
            long_rows.append(mdf)
            s = summarize(mdf, series, cc, args.kge_tol)
            s["candidate_flags"] = c["flags"] if series == SER_ON else ""
            summ_rows.append(s)
        if not args.no_pairs:
            lab = label_for_cc(cc)
            write_pair(res_on,  per_series[SER_ON],  cc, out_dir / f"stage2_ccON_vs_truth_cc{lab}.csv",
                       phase_fn is not None)
            write_pair(res_off, per_series[SER_OFF], cc, out_dir / f"control_ccOFF_vs_truth_cc{lab}.csv",
                       phase_fn is not None)

    long_df = pd.concat(long_rows, ignore_index=True)
    summ_df = pd.DataFrame(summ_rows)
    long_df.to_csv(out_dir / "rescored_long_110.csv", index=False)
    summ_df.to_csv(out_dir / "rescored_truth_summary_110.csv", index=False)

    prov = {
        "script": Path(__file__).name,
        "created_local": datetime.now().isoformat(timespec="seconds"),
        "pandas": pd.__version__, "numpy": np.__version__,
        "cc_off_truth": {"file": ccoff_path.name, "md5": md5_of(ccoff_path)},
        "inputs": {STAGE2_NAME: {"path": str(p_on), "md5": md5_of(p_on), "runs": int(len(res_on))},
                   CONTROL_NAME: {"path": str(p_off), "md5": md5_of(p_off), "runs": int(len(res_off))}},
        "gates": {"A_max_abs_diff_m3s": worstA, "B_worst_abs_dev": worstB,
                  "runs_checked": n_checked,
                  "runs_excluded": {k: v for k, v in bad.items()},
                  "phase_metrics_from_scorer": phase_fn is not None},
        "candidate_truths": {c["name"]: {"cc_mmhr": c["cc"], "md5": md5_of(cand_dir / c["name"])}
                             for c in usable},
        "skipped_truths": [{"file": n, "reason": w} for n, w in skipped],
        "kge_tol": args.kge_tol, "feasible_pbias_pct": FEASIBLE_PBIAS_PCT,
        "note": "Metrics recomputed from the stored Simulated column of each compare CSV against each "
                "candidate truth read exactly as the scorer reads a truth.",
    }
    (out_dir / "PROVENANCE_rescored_cc_truths_110.json").write_text(json.dumps(prov, indent=2, default=str))

    report(summ_df, skipped)
    if truth_failed:
        print(f"\n  WARNING: {len(truth_failed)} (truth, run) pair(s) could not be scored; first: {truth_failed[0]}")
    print(f"\nSaved to: {out_dir}")
    print("  rescored_truth_summary_110.csv   rescored_long_110.csv   PROVENANCE_rescored_cc_truths_110.json")
    if not args.no_pairs:
        print("  stage2_ccON_vs_truth_cc<V>.csv / control_ccOFF_vs_truth_cc<V>.csv  (one pair per truth; "
              "feed to analyze_cc_pareto_110.py)")


if __name__ == "__main__":
    main()
