"""
rescore_truth_location_110.py
=============================
Series 110 -- the TRUTH-LOCATION test, part 2: RE-SCORE the existing Stage 2
(cc ON) and control (cc OFF) runs against every new truth made by
run_truth_location_110.py. No tRIBS, no new simulations.

WHAT THIS DOES (plain language)
--------------------------------
Stage 2 ran 250 joint (Ks, f, cc) simulations with channel loss ON, and the
control ran the same 250 points with channel loss OFF. Each of those 500 runs
is a candidate GUESS with a stored hydrograph. A simulation does not depend on
which truth it is later compared to -- only its SCORE does. So the 500 stored
hydrographs can be compared to any new truth in seconds.

This script takes each new truth (a different true Ks, f and cc), keeps every
run's stored "Simulated" series untouched, swaps in the new truth as
"Observed" (read through the very same reading and resampling calls the scorer
uses), and recomputes the same metrics with the same formulas. It reuses the
reading, gate and scoring code of rescore_cc_truths_110.py, so the numbers are
computed exactly the way your earlier re-scoring computed them.

SELF-CHECKS ("gates") -- run before anything is re-scored
-----------------------------------------------------------
  GATE A  The current cc-OFF truth (synth_truth/), read by this script's code,
          must equal the "Observed" column stored in all 500 compare CSVs.
  GATE B  Every metric recomputed from a compare CSV must equal the one stored
          in the Stage 2 / control results CSV.
  GATE C  Each new truth file must still match the checksum recorded when it was
          made (PROVENANCE_truth_location_110.json), come from a run that passed
          its sanity check, lie inside the sampled box, and share its 5-minute
          time grid with the cc-OFF truth.
  GATE D  REPRODUCIBILITY. New truths at the original centre (Ks 7.0, f 0.012)
          whose cc matches one of the earlier cc candidates must reproduce those
          earlier files (largest difference below 1e-6 of the peak). This shows
          today's build is the build that made the earlier truths, so the new
          truths and the 500 stored runs are comparable.
If a gate fails the script stops and says why; nothing is re-scored.
(--drop_bad_runs excludes individual runs that fail A/B and continues;
--ignore_reproduction continues past a gate D failure, and records that in
the provenance file -- only do that if you understand why it failed.)
If GATE A fails for every run with a difference of a fraction of a m3/s, suspect a
change of pandas version since the Stage 2 runs were scored (resampling a series with
.resample("5min").interpolate(method="time") gives slightly different values in pandas 2
and pandas 3). Re-scoring needs the same pandas that scored the runs; do not use
--drop_bad_runs or edit anything to get past it.

WHAT YOU GET   (calibration_work/03_comparisons/summary_tables/rescored_truth_location_110/)
----------------------------------------------------------------------------------------------
  rescored_location_long_110.csv       every run x every truth: all metrics, plus the
                                       truth's own Ks, f and cc (truth_Ks, truth_f, truth_cc)
  rescored_location_summary_110.csv    one row per (truth, series): best KGE_2012, how many
                                       runs are feasible, where the near-best runs sit in cc,
                                       how close the nearest sample is to the truth point
  PROVENANCE_rescored_truth_location_110.json   checksums, versions, gate results

The next script, analyze_truth_location_110.py, reads these files.

HOW TO READ THE TABLE (rules of thumb, mine)
----------------------------------------------
  - "ON best" is the best KGE_2012 any cc-ON run reached against that truth; "OFF
    best" is the best the cc-OFF runs reached (they can only move Ks and f). A small
    gap means Ks and f alone can imitate that truth.
  - "near-best cc range" is the span of cc among cc-ON runs within --kge_tol of the
    best. A narrow range containing the true cc is what identifiable looks like.
  - "nearest" is the distance, in units of the sampled box (0 = same point, 1 = a
    full side of the box), from the truth point to the closest of the 250 samples;
    "n<.15" counts samples within 0.15 of it (a uniform sweep puts about 3-4 there
    for a truth in the middle, fewer near an edge). A truth with a large "nearest"
    and a poor best KGE says more about the sweep's coverage than about cc.
  - One storm, 250 points in 3 dimensions, routing pinned at truth, noise-free
    truth: this is a best case, not a realistic one. Feasible = |PBIAS| < 2%.

SAFETY
-------
  - Read-only on every existing file. Writes only inside its own output folder
    (refuses an output folder that is, or overlaps, any input location).
  - Never imports the builder, never touches current_run_config.json, never starts
    tRIBS. Do not run it WHILE run_truth_location_110.py is still running: it needs
    the finished set of truths.
  - Run membership comes only from the run_id rows in the two results CSVs.

USAGE (run from the lab/ directory; keep rescore_cc_truths_110.py in the same folder)
----------------------------------------------------------------------------------------
    python rescore_truth_location_110.py --check_only     # gates only; writes nothing
    python rescore_truth_location_110.py                  # gates, then every truth
    python rescore_truth_location_110.py --kge_tol 0.005
"""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import rescore_cc_truths_110 as RS
except ImportError:
    sys.exit("This script needs rescore_cc_truths_110.py in the same folder (it reuses its reading, gate "
             "and scoring code). Put both files in lab/ and run from there.")

_NEEDED = ["load_results", "load_all_compares", "gate_a", "gate_b", "read_truth_5min", "compute_metrics",
           "import_phase_fn", "md5_of", "find_candidates", "STAGE2_NAME", "CONTROL_NAME", "CAND_DIRNAME",
           "SER_ON", "SER_OFF", "BOX", "FEASIBLE_PBIAS_PCT", "NEAR_RADIUS", "DEFAULT_EXPECT_N",
           "DEFAULT_KGE_TOL", "_unit"]
_missing = [n for n in _NEEDED if not hasattr(RS, n)]
if _missing:
    sys.exit("rescore_cc_truths_110.py in this folder is an older version (missing: %s). Replace it with the "
             "current copy." % ", ".join(_missing))

SER_ON, SER_OFF = RS.SER_ON, RS.SER_OFF
LOC_DIRNAME = "synth_truth_location"
LOC_PROV = "PROVENANCE_truth_location_110.json"
OUT_DIRNAME = "rescored_truth_location_110"
LONG_NAME = "rescored_location_long_110.csv"
SUMMARY_NAME = "rescored_location_summary_110.csv"
PROV_NAME = "PROVENANCE_rescored_truth_location_110.json"
STANDARD_KS, STANDARD_F = 7.0, 0.012
REPRO_REL = 1e-6


# ------------------------------------------------------------------
# Small helpers
# ------------------------------------------------------------------
def _g(v):
    return ("%g" % v).replace(".", "p")


def truth_id(ks, f, cc):
    return "K%s_f%s_cc%s" % (_g(ks), _g(f), _g(cc))


def inside_box(ks, f, cc):
    vals = {"Ks_mult": ks, "f_RS_abs": f, "channelconductivity_mmhr": cc}
    return all(RS.BOX[k][0] * (1 - 1e-9) <= vals[k] <= RS.BOX[k][1] * (1 + 1e-9) for k in vals)


def box_unit(value, key):
    lo, hi, scale = RS.BOX[key]
    return RS._unit(value, lo, hi, "log" if scale == "log" else "lin")


def load_registry(loc_dir):
    """The new truths as recorded when they were made. Returns (list of dicts, unregistered file names)."""
    prov_path = loc_dir / LOC_PROV
    if not prov_path.exists():
        sys.exit("No %s in %s. Run run_truth_location_110.py first (and let it finish)." % (LOC_PROV, loc_dir))
    try:
        files = json.loads(prov_path.read_text()).get("files", {})
    except Exception as e:
        sys.exit("Could not read %s: %s" % (prov_path, e))
    if not files:
        sys.exit("%s lists no truth files." % prov_path.name)
    reg = []
    for name, info in files.items():
        try:
            ks, f, cc = float(info["Ks_mult"]), float(info["f_RS_abs"]), float(info["cc_mmhr"])
        except (KeyError, TypeError, ValueError):
            sys.exit("%s: the entry for %s has no usable Ks_mult / f_RS_abs / cc_mmhr." % (prov_path.name, name))
        reg.append({"name": name, "Ks": ks, "f": f, "cc": cc, "md5": str(info.get("md5", "")),
                    "flags": str(info.get("flags", "") or ""), "id": truth_id(ks, f, cc)})
    ids = [r["id"] for r in reg]
    if len(set(ids)) != len(ids):
        sys.exit("%s lists two truths with the same (Ks, f, cc); refusing." % prov_path.name)
    on_disk = {p.name for p in loc_dir.glob("*.qout")}
    unregistered = sorted(on_disk - set(files))
    reg.sort(key=lambda r: (r["Ks"], r["f"], r["cc"]))
    return reg, unregistered


def score_truth(obs_new, res, cache, phase_fn):
    """Metrics of every run in `res` against one truth. Same loop as rescore_cc_truths_110.py."""
    rows, failed = [], []
    for _, r in res.iterrows():
        df = cache[r["run_id"]]
        o = obs_new.reindex(df.index)
        if o.isna().any():
            failed.append(r["run_id"])
            continue
        m = RS.compute_metrics(o, df["Simulated"], phase_fn)
        rows.append({"run_id": r["run_id"], "Ks_mult": r["Ks_mult"], "f_RS_abs": r["f_RS_abs"],
                     "channelconductivity_mmhr": r["channelconductivity_mmhr"],
                     "optpercolation": r["optpercolation"], **m})
    return pd.DataFrame(rows), failed


def summarize(m, series, tr, kge_tol):
    """One summary row for one (truth, series). m has one row per run."""
    out = {"truth_id": tr["id"], "truth_Ks": tr["Ks"], "truth_f": tr["f"], "truth_cc": tr["cc"],
           "series": series, "n_scored": int(len(m))}
    ok = m[np.isfinite(m["kge_2012"])]
    out["n_valid_kge"] = int(len(ok))
    if ok.empty:
        return out
    best = ok.loc[ok["kge_2012"].idxmax()]
    out.update({
        "best_kge_2012": float(best["kge_2012"]), "best_run_id": best["run_id"],
        "best_Ks": float(best["Ks_mult"]), "best_f": float(best["f_RS_abs"]),
        "best_pbias_pct": float(best["pbias_pct"]),
        "best_dKs_from_truth": float(best["Ks_mult"] - tr["Ks"]),
        "best_dlog10_f_from_truth": float(np.log10(best["f_RS_abs"] / tr["f"])),
    })
    feas = ok[np.abs(ok["pbias_pct"]) < RS.FEASIBLE_PBIAS_PCT]
    out["n_feasible"] = int(len(feas))
    out["pct_feasible"] = 100.0 * len(feas) / len(ok)
    near = ok[ok["kge_2012"] >= best["kge_2012"] - kge_tol]
    out["n_near_best"] = int(len(near))
    out["kge_tol"] = kge_tol
    if series == SER_ON:
        cc = near["channelconductivity_mmhr"].astype(float)
        lo, hi = float(cc.min()), float(cc.max())
        cc_box = RS.BOX["channelconductivity_mmhr"]
        span = np.log10(hi / lo) / np.log10(cc_box[1] / cc_box[0])
        out.update({
            "near_cc_min": lo, "near_cc_median": float(cc.median()), "near_cc_max": hi,
            "near_cc_logspan_frac_of_box": float(span),
            "truth_cc_inside_near_range": bool(lo <= tr["cc"] <= hi),
            "best_cc": float(best["channelconductivity_mmhr"]),
            "best_dlog10_cc_from_truth": float(np.log10(best["channelconductivity_mmhr"] / tr["cc"])),
        })
        pts = np.column_stack([
            box_unit(ok["Ks_mult"].to_numpy(float), "Ks_mult") - box_unit(tr["Ks"], "Ks_mult"),
            box_unit(ok["f_RS_abs"].to_numpy(float), "f_RS_abs") - box_unit(tr["f"], "f_RS_abs"),
            box_unit(ok["channelconductivity_mmhr"].to_numpy(float), "channelconductivity_mmhr")
            - box_unit(tr["cc"], "channelconductivity_mmhr"),
        ])
        dist = np.sqrt((pts ** 2).sum(axis=1))
        i = int(np.argmin(dist))
        out.update({
            "nearest_sample_dist": float(dist[i]),
            "nearest_sample_run_id": ok.iloc[i]["run_id"],
            "nearest_sample_kge": float(ok.iloc[i]["kge_2012"]),
            "n_samples_within_0p15": int((dist <= RS.NEAR_RADIUS).sum()),
        })
    return out


def fmt(v, spec):
    return "   n/a" if (v is None or (isinstance(v, float) and not np.isfinite(v))) else format(v, spec)


def report(summ, skipped):
    on = summ[summ["series"] == SER_ON].set_index("truth_id")
    off = summ[summ["series"] == SER_OFF].set_index("truth_id")
    print("\n" + "=" * 128)
    print("RE-SCORED AGAINST EACH NEW TRUTH   (feasible = |PBIAS| < %g%%;  near-best = within kge_tol of the best KGE_2012)"
          % RS.FEASIBLE_PBIAS_PCT)
    print("=" * 128)
    print("%5s %7s %6s | %7s %8s %6s | %9s | %-24s %5s %4s | %7s %5s"
          % ("Ks", "f", "cc", "ON best", "OFF best", "gap", "feas ON/OFF", "near-best cc range (ON)", "span", "in?",
             "nearest", "n<.15"))
    order = on.sort_values(["truth_Ks", "truth_f", "truth_cc"]).index
    last = None
    for tid in order:
        a = on.loc[tid]
        b = off.loc[tid] if tid in off.index else None
        cell = (a["truth_Ks"], a["truth_f"])
        if last is not None and cell != last:
            print("")
        last = cell
        gap = (a["best_kge_2012"] - b["best_kge_2012"]) if b is not None else np.nan
        rng = ("%.0f-%.0f (med %.0f)" % (a["near_cc_min"], a["near_cc_max"], a["near_cc_median"])
               if np.isfinite(a.get("near_cc_min", np.nan)) else "n/a")
        feas_on = int(a["n_feasible"]) if np.isfinite(a.get("n_feasible", np.nan)) else -1
        feas_off = int(b["n_feasible"]) if b is not None and np.isfinite(b.get("n_feasible", np.nan)) else -1
        print("%5g %7g %6g | %s %s %s | %4d/%-4d | %-24s %s %4s | %s %5d"
              % (a["truth_Ks"], a["truth_f"], a["truth_cc"], fmt(a["best_kge_2012"], "7.4f"),
                 fmt(b["best_kge_2012"], "8.4f") if b is not None else "     n/a", fmt(gap, "6.3f"),
                 feas_on, feas_off, rng, fmt(a.get("near_cc_logspan_frac_of_box", np.nan), "5.2f"),
                 "yes" if a.get("truth_cc_inside_near_range") else "no",
                 fmt(a.get("nearest_sample_dist", np.nan), "7.3f"), int(a.get("n_samples_within_0p15", 0))))
    print("\n  ON best / OFF best : best KGE_2012 of the cc-ON runs / of the cc-OFF control runs, both vs this truth")
    print("  gap                : ON best minus OFF best (small = Ks and f alone can imitate this truth)")
    print("  feas ON/OFF        : number of runs with |PBIAS| < %g%% in each series (same 250 points)" % RS.FEASIBLE_PBIAS_PCT)
    print("  span / in?         : share of the cc box (log scale) covered by the near-best runs / does it contain the true cc")
    print("  nearest / n<.15    : distance from the truth point to the closest sample, in box units / samples within 0.15")
    print("  Caveat: one storm, 250 points in 3-D, routing pinned at truth, noise-free truth -- a best case.")
    if skipped:
        print("\n  SKIPPED TRUTHS:")
        for name, why in skipped:
            print("    %s: %s" % (name, why))


# ------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Series 110 -- re-score the Stage 2 and control runs against every truth-location truth (no tRIBS).")
    ap.add_argument("--check_only", action="store_true", help="Run the gates and stop; writes nothing.")
    ap.add_argument("--kge_tol", type=float, default=RS.DEFAULT_KGE_TOL,
                    help="Near-best tolerance in KGE_2012 units (default %g)" % RS.DEFAULT_KGE_TOL)
    ap.add_argument("--expect_n", type=int, default=RS.DEFAULT_EXPECT_N,
                    help="Expected runs per series (default %d); a different count only warns" % RS.DEFAULT_EXPECT_N)
    ap.add_argument("--drop_bad_runs", action="store_true",
                    help="Exclude runs that fail gate A or B and continue, instead of stopping.")
    ap.add_argument("--ignore_reproduction", action="store_true",
                    help="Continue even if gate D (reproduction of the earlier cc candidates) fails. Recorded in the "
                         "provenance file.")
    ap.add_argument("--stage2_csv", type=Path, default=None)
    ap.add_argument("--control_csv", type=Path, default=None)
    ap.add_argument("--out_dir", type=Path, default=None)
    args = ap.parse_args()

    script_dir = Path.cwd()
    calib_dir = script_dir.parent / "calibration_work"
    summary_dir = calib_dir / "03_comparisons" / "summary_tables"
    csv_dir = calib_dir / "03_comparisons" / "csv_exports"
    synth_dir = calib_dir / "synth_truth"
    old_cand_dir = calib_dir / RS.CAND_DIRNAME
    loc_dir = calib_dir / LOC_DIRNAME
    out_dir = args.out_dir or (summary_dir / OUT_DIRNAME)
    p_on = args.stage2_csv or (summary_dir / RS.STAGE2_NAME)
    p_off = args.control_csv or (summary_dir / RS.CONTROL_NAME)

    print("\n" + "=" * 78 + "\nSeries 110 -- re-score Stage 2 / control against the truth-location truths\n" + "=" * 78)

    # ---- inputs / safety -------------------------------------------------
    qouts = list(synth_dir.glob("*.qout")) if synth_dir.exists() else []
    if len(qouts) != 1:
        sys.exit("Expected exactly one *.qout in %s (the cc-OFF truth the 500 runs were scored against), "
                 "found %d: %s" % (synth_dir, len(qouts), [q.name for q in qouts]))
    ccoff_path = qouts[0]
    print("cc-OFF truth (read-only): %s" % ccoff_path.name)

    out_res = out_dir.resolve()
    for guarded in (synth_dir, old_cand_dir, loc_dir, csv_dir):
        gd = guarded.resolve()
        if out_res == gd or gd in out_res.parents or out_res in gd.parents:
            sys.exit("Output folder %s overlaps an input folder (%s); refusing." % (out_dir, guarded))
    for pth in (p_on, p_off):
        if out_res in pth.resolve().parents or out_res == pth.resolve():
            sys.exit("Output folder %s contains an input file (%s); refusing." % (out_dir, pth))

    registry, unregistered = load_registry(loc_dir)
    print("Truth-location registry: %d truth(s) recorded in %s" % (len(registry), LOC_PROV))
    if unregistered:
        print("  NOTE: %d .qout file(s) in %s are not in the registry and are ignored: %s"
              % (len(unregistered), loc_dir.name, ", ".join(unregistered[:4]) + (" ..." if len(unregistered) > 4 else "")))

    print("\nLoading results tables:")
    res_on = RS.load_results(p_on, SER_ON, 1, args.expect_n)
    res_off = RS.load_results(p_off, SER_OFF, 0, args.expect_n)

    phase_fn = RS.import_phase_fn(script_dir)
    obs_ccoff = RS.read_truth_5min(ccoff_path)

    print("\nLoading compare CSVs (run membership = run_id rows above):")
    cache_on, prob_on = RS.load_all_compares(res_on, csv_dir)
    cache_off, prob_off = RS.load_all_compares(res_off, csv_dir)
    print("  %s: %d/%d loaded;  %s: %d/%d loaded" % (SER_ON, len(cache_on), len(res_on),
                                                     SER_OFF, len(cache_off), len(res_off)))

    # ---- gates A and B ---------------------------------------------------
    bad = {}
    for rid, why in {**prob_on, **prob_off}.items():
        bad[rid] = why
    cache = {**cache_on, **cache_off}
    allres = pd.concat([res_on.assign(_series=SER_ON), res_off.assign(_series=SER_OFF)], ignore_index=True)
    badA, worstA = RS.gate_a(allres, cache, obs_ccoff)
    badB, worstB = RS.gate_b(allres, cache, phase_fn)
    for d in (badA, badB):
        for rid, why in d.items():
            bad.setdefault(rid, why)

    n_checked = len(cache)
    print("\nGATE A  observed series : %d/%d compare CSVs match the cc-OFF truth as read here (worst abs diff %.2e m3/s)"
          % (n_checked - len(badA), n_checked, worstA))
    if worstB:
        wk = max(worstB, key=worstB.get)
        print("GATE B  stored metrics  : %d/%d runs reproduce their stored metrics (worst abs deviation %.2e in %s; %s)"
              % (n_checked - len(badB), n_checked, worstB[wk], wk,
                 "phase metrics included" if phase_fn else "phase metrics skipped"))
    else:
        print("GATE B  stored metrics  : %d/%d" % (n_checked - len(badB), n_checked))

    if bad:
        print("\n  %d run(s) failed a gate or could not be loaded. First few:" % len(bad))
        for rid, why in list(bad.items())[:8]:
            print("    %s: %s" % (rid, why))
        if not args.drop_bad_runs:
            sys.exit("\nSTOPPED: nothing was re-scored. Fix the cause (or re-run with --drop_bad_runs to exclude "
                     "these runs and continue).")
        print("  --drop_bad_runs set: these runs are excluded.")
    res_on = res_on[~res_on["run_id"].isin(bad)].reset_index(drop=True)
    res_off = res_off[~res_off["run_id"].isin(bad)].reset_index(drop=True)
    if res_on.empty or res_off.empty:
        sys.exit("No runs left in one of the series after the gates; stopping.")

    # ---- gate C: the new truths ------------------------------------------
    usable, skipped, truth_obs = [], [], {}
    for tr in registry:
        p = loc_dir / tr["name"]
        label = "Ks %g f %g cc %g" % (tr["Ks"], tr["f"], tr["cc"])
        if not p.exists():
            skipped.append((label, "file missing")); continue
        if "BAD_QOUT" in tr["flags"]:
            skipped.append((label, "the file failed its sanity check when it was made (BAD_QOUT); re-run "
                                   "run_truth_location_110.py with --skip_existing to rebuild it")); continue
        if not inside_box(tr["Ks"], tr["f"], tr["cc"]):
            skipped.append((label, "outside the sampled box; a sweep that never looked there cannot recover it")); continue
        if not tr["md5"]:
            skipped.append((label, "no recorded checksum")); continue
        if RS.md5_of(p) != tr["md5"]:
            skipped.append((label, "checksum differs from the one recorded when it was made (GATE C)")); continue
        try:
            s = RS.read_truth_5min(p)
        except Exception as e:
            skipped.append((label, "unreadable: %s" % str(e)[:80])); continue
        if not (s.index.equals(obs_ccoff.index) and (s.isna() == obs_ccoff.isna()).all()):
            skipped.append((label, "5-minute grid differs from the cc-OFF truth's (GATE C)")); continue
        usable.append(tr)
        truth_obs[tr["id"]] = s
    print("GATE C  new truths      : %d/%d verified (checksum, sanity flag, box, 5-minute grid)"
          % (len(usable), len(registry)))
    for name, why in skipped:
        print("    skipped %s: %s" % (name, why))
    if not usable:
        sys.exit("No usable truths.")

    # ---- gate D: reproduction of the earlier candidates --------------------
    repro = []
    old_cands = RS.find_candidates(old_cand_dir, summary_dir) if old_cand_dir.exists() else []
    for tr in usable:
        if not (np.isclose(tr["Ks"], STANDARD_KS) and np.isclose(tr["f"], STANDARD_F)):
            continue
        hit = [c for c in old_cands if np.isclose(c["cc"], tr["cc"], rtol=1e-3)]
        if not hit:
            continue
        old_path = old_cand_dir / hit[0]["name"]
        if not old_path.exists():
            continue
        entry = {"truth_id": tr["id"], "cc": tr["cc"], "earlier_file": hit[0]["name"]}
        try:
            old_s = RS.read_truth_5min(old_path)
            new_s = truth_obs[tr["id"]]
            if not old_s.index.equals(new_s.index):
                entry.update({"ok": False, "detail": "the two files have different time grids"})
            else:
                peak = float(np.nanmax(np.abs(old_s.to_numpy(float)))) or 1.0
                d = float(np.nanmax(np.abs(old_s.to_numpy(float) - new_s.to_numpy(float))))
                entry.update({"ok": bool(d <= REPRO_REL * peak), "max_abs_diff_m3s": d, "rel_to_peak": d / peak,
                              "detail": "largest difference %.2e m3/s (%.1e of the peak)" % (d, d / peak)})
        except Exception as e:
            entry.update({"ok": False, "detail": "could not compare: %s" % str(e)[:80]})
        repro.append(entry)
    if not repro:
        n_centre = sum(1 for tr in usable if np.isclose(tr["Ks"], STANDARD_KS) and np.isclose(tr["f"], STANDARD_F))
        if not n_centre:
            why = "no new truth at the original centre (Ks %g, f %g)" % (STANDARD_KS, STANDARD_F)
        elif not old_cands:
            why = "%d new truth(s) at the original centre, but no earlier cc candidates were found in %s" \
                  % (n_centre, old_cand_dir.name)
        else:
            why = "%d new truth(s) at the original centre, but none has a cc that matches an earlier candidate" % n_centre
        print("GATE D  reproduction    : NOT CHECKED (%s)" % why)
    else:
        n_ok = sum(1 for e in repro if e["ok"])
        print("GATE D  reproduction    : %d/%d centre truths reproduce the earlier cc candidates" % (n_ok, len(repro)))
        for e in repro:
            print("    cc %-6g %s  (%s)" % (e["cc"], "ok  " if e["ok"] else "FAIL", e["detail"]))
        if n_ok < len(repro):
            if not args.ignore_reproduction:
                sys.exit("\nSTOPPED: today's build does not reproduce the earlier truths, so the new truths may not be "
                         "comparable with the 500 stored runs. Find out why first (different builder, parameters, or "
                         "tRIBS build?). --ignore_reproduction overrides this, and says so in the provenance file.")
            print("  --ignore_reproduction set: continuing despite the failure.")

    if args.check_only:
        print("\n--check_only: gates finished, nothing written.")
        return

    # ---- re-score ----------------------------------------------------------
    out_dir.mkdir(parents=True, exist_ok=True)
    print("\nRe-scoring %d + %d runs x %d truth(s)..." % (len(res_on), len(res_off), len(usable)), flush=True)
    long_rows, summ_rows, truth_failed = [], [], []
    for k, tr in enumerate(usable, start=1):
        obs_new = truth_obs[tr["id"]]
        for series, res in ((SER_ON, res_on), (SER_OFF, res_off)):
            mdf, failed = score_truth(obs_new, res, cache, phase_fn)
            for rid in failed:
                truth_failed.append((tr["id"], rid))
            if mdf.empty:
                continue
            mdf.insert(0, "series", series)
            mdf.insert(0, "truth_cc", tr["cc"])
            mdf.insert(0, "truth_f", tr["f"])
            mdf.insert(0, "truth_Ks", tr["Ks"])
            mdf.insert(0, "truth_id", tr["id"])
            long_rows.append(mdf)
            summ_rows.append(summarize(mdf, series, tr, args.kge_tol))
        if k % 5 == 0 or k == len(usable):
            print("  %d/%d truths scored" % (k, len(usable)), flush=True)

    long_df = pd.concat(long_rows, ignore_index=True)
    summ_df = pd.DataFrame(summ_rows)
    long_df.to_csv(out_dir / LONG_NAME, index=False)
    summ_df.to_csv(out_dir / SUMMARY_NAME, index=False)

    prov = {
        "script": Path(__file__).name,
        "created_local": datetime.now().isoformat(timespec="seconds"),
        "pandas": pd.__version__, "numpy": np.__version__,
        "cc_off_truth": {"file": ccoff_path.name, "md5": RS.md5_of(ccoff_path)},
        "inputs": {RS.STAGE2_NAME: {"path": str(p_on), "md5": RS.md5_of(p_on), "runs": int(len(res_on))},
                   RS.CONTROL_NAME: {"path": str(p_off), "md5": RS.md5_of(p_off), "runs": int(len(res_off))}},
        "gates": {"A_max_abs_diff_m3s": worstA, "B_worst_abs_dev": worstB, "runs_checked": n_checked,
                  "runs_excluded": {k: v for k, v in bad.items()},
                  "phase_metrics_from_scorer": phase_fn is not None,
                  "D_reproduction": repro, "D_overridden": bool(args.ignore_reproduction and repro and
                                                                not all(e["ok"] for e in repro))},
        "truths": {tr["id"]: {"file": tr["name"], "Ks_mult": tr["Ks"], "f_RS_abs": tr["f"], "cc_mmhr": tr["cc"],
                              "md5": RS.md5_of(loc_dir / tr["name"])} for tr in usable},
        "skipped_truths": [{"truth": n, "reason": w} for n, w in skipped],
        "unregistered_files": unregistered,
        "kge_tol": args.kge_tol, "feasible_pbias_pct": RS.FEASIBLE_PBIAS_PCT,
        "note": "Metrics recomputed from the stored Simulated column of each compare CSV against each truth, read "
                "exactly as the scorer reads a truth (rescore_cc_truths_110.py's code).",
    }
    (out_dir / PROV_NAME).write_text(json.dumps(prov, indent=2, default=str))

    report(summ_df, skipped)
    if truth_failed:
        print("\n  WARNING: %d (truth, run) pair(s) could not be scored; first: %s" % (len(truth_failed), truth_failed[0]))
    print("\nSaved to: %s" % out_dir)
    print("  %s   %s   %s" % (LONG_NAME, SUMMARY_NAME, PROV_NAME))
    print("Next: python analyze_truth_location_110.py")


if __name__ == "__main__":
    main()
