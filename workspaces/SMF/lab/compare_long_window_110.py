"""
compare_long_window_110.py   (version 1, read-only: it writes nothing)

PURPOSE
  Check that running tRIBS over the long window (label sbl, 1416 h) gives the same 12 Aug result as the
  450-hour Stage B run (label sb1) for the SAME parameter set (Ks 6.304, f 0.0007, cc 62.351 mm/hr).
  Also prints the highest model flow on each of the four 2014 storm dates as a sanity look (not a score).

RUN (from the lab/ folder)
  python compare_long_window_110.py

READING RULE (fixed before looking)
  The two kept Outlet .qout files should agree for every hour they share, or differ only in the last digits.
  If they differ by more than that (more than 0.001 m3/s anywhere in hours 0-300), STOP: the long window changes
  the 12 Aug result and that has to be understood before the big run.
"""
import sys
from pathlib import Path

import pandas as pd

CC_TAG = "cc62p351454"
CALIB = Path("../calibration_work")
KEEP = CALIB / "kept_qout_110"
SUMM = CALIB / "03_comparisons" / "summary_tables"
TOL = 1e-3          # m3/s
START = pd.Timestamp("2014-08-01")


def read_qout(path):
    """Time_hr and Qstrm_m3s of a kept Outlet .qout (comma- or whitespace-delimited; first line is a header)."""
    if not path.exists():
        sys.exit("Cannot find %s. Run this from the lab/ folder." % path)
    with open(path) as fh:
        fh.readline()
        first = fh.readline().rstrip()
    print("  %s: first data line is  %s" % (path.name, first))
    df = pd.read_csv(path, sep=r"[,\s]+", engine="python", skiprows=1, header=None,
                     usecols=[0, 1], names=["t", "q"])
    df["t"] = df["t"].round(4)
    return df


a = read_qout(KEEP / ("SMF_20140812_110_%s_sb1on_Outlet.qout" % CC_TAG))
b = read_qout(KEEP / ("SMF_20140812_110_%s_sblon_Outlet.qout" % CC_TAG))
print()
print("sb1on (450 h run): %d rows, last hour %.2f" % (len(a), a.t.max()))
print("sblon (long run):  %d rows, last hour %.2f" % (len(b), b.t.max()))

j = a.merge(b, on="t", suffixes=("_sb1", "_sbl"))
dq = (j.q_sb1 - j.q_sbl).abs()
print("hours in common: %d   biggest flow difference over all of them: %.3g m3/s" % (len(j), dq.max()))
ev = j[(j.t >= 280) & (j.t <= 300)]
evd = (ev.q_sb1 - ev.q_sbl).abs().max()
print("12 Aug scoring window (hours 280-300): %d rows, biggest difference %.3g m3/s, peaks %.4f vs %.4f m3/s"
      % (len(ev), evd, ev.q_sb1.max(), ev.q_sbl.max()))
early = j[j.t <= 300]
print("hours 0-300: biggest difference %.3g m3/s" % (early.q_sb1 - early.q_sbl).abs().max())
if len(j) == 0 or len(ev) == 0:
    print("READING: the two files share no hours; something is wrong with the time column. Send me this output.")
elif (early.q_sb1 - early.q_sbl).abs().max() <= TOL:
    print("READING: AGREE (within %.0e m3/s). The long window does not change the 12 Aug result." % TOL)
else:
    print("READING: DIFFER by more than %.0e m3/s in hours 0-300. Do not start the big run; send me this output." % TOL)

print()
print("Stored scores of the two runs (from the results tables):")
for name in ("lhs_results_lowf_sb1_ON_110.csv", "lhs_results_lowf_sbl_ON_110.csv"):
    r = pd.read_csv(SUMM / name)
    r = r[r.run_id.str.contains(CC_TAG)]
    cols = [c for c in ("run_id", "kge_2012", "pbias_pct", "sim_volume_m3") if c in r.columns]
    print(r[cols].to_string(index=False))

print()
print("Sanity look at the long run (not a score): highest model flow on each storm date")
for lab, day in (("12 Aug", "2014-08-12"), ("19 Aug", "2014-08-19"), ("8 Sep", "2014-09-08"), ("27 Sep", "2014-09-27")):
    h0 = (pd.Timestamp(day) - START).total_seconds() / 3600
    w = b[(b.t >= h0) & (b.t < h0 + 24)]
    if w.empty:
        print("  %-7s no output in the file for this date" % lab)
        continue
    i = w.q.idxmax()
    print("  %-7s highest flow %.3f m3/s (%.1f cfs) at %s"
          % (lab, w.q.max(), w.q.max() / 0.0283168,
             (START + pd.Timedelta(hours=float(b.t[i]))).strftime("%d %b %H:%M")))
