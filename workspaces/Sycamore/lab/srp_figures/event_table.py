#!/usr/bin/env python
"""Event numbers for the five WY2017 storms: rain, peak flow, lag, direct runoff and runoff ratio.

Inputs: init_data/obs/sycamore_09510200_obs_flow_WY2017.csv (15-minute flow, MST)
        init_data/obs/sycamore_09510200_daily_flow_WY2016_2026.rdb (USGS daily means, for the whole-year volume)
        out/basin_mean_hourly_rain_wy2017.csv (from load_radar.py)
Writes out/event_table_wy2017.csv. Imported by plot5_event_table.py; run alone it prints the table and checks.

Direct runoff = flow above a straight baseline drawn from the start of the rise to N days after the peak,
N = 0.827 * A**0.2 (A in square miles; the Linsley-Kohler-Paulhus rule of thumb for how long direct runoff lasts).
"""
import numpy as np
import pandas as pd
import srp_common as c

AREA_KM2 = 424.6153                      # delineated basin (data_sources.md / Make_Sycamore_Model.ipynb)
AREA_MI2 = AREA_KM2 / 2.589988
N_DAYS = 0.827 * AREA_MI2 ** 0.2
CFS_TO_CMS = 0.0283168
ACFT = 1233.48                           # m3 per acre-foot


def _trapz(y, t):
    return float(np.trapezoid(y, t)) if hasattr(np, "trapezoid") else float(np.trapz(y, t))


def load():
    q = pd.read_csv(c.INIT / "obs" / "sycamore_09510200_obs_flow_WY2017.csv", parse_dates=["datetime"]).set_index("datetime")
    r = pd.read_csv(c.OUT / "basin_mean_hourly_rain_wy2017.csv", parse_dates=["interval_start_mst"]).set_index("interval_start_mst").basin_mean_mm
    import io
    lines = [ln for ln in open(c.INIT / "obs" / "sycamore_09510200_daily_flow_WY2016_2026.rdb") if not ln.startswith("#")]
    del lines[1]                              # the RDB format row ("5s 15s 20d ...")
    dv = pd.read_csv(io.StringIO("".join(lines)), sep="\t")
    dv["datetime"] = pd.to_datetime(dv["datetime"], format="%Y-%m-%d")
    dv = dv.set_index("datetime").loc["2016-10-01":"2017-09-30"]
    cfs = dv[[k for k in dv.columns if k.endswith("00003")][0]].astype(float)
    assert len(cfs) == 365 and cfs.notna().all()
    return q, r, cfs


def direct_runoff(q, peak_t, n_days=N_DAYS, frac=0.005, look_back_h=36):
    """Start of rise, end of direct runoff, and the runoff volume (m3) above a straight baseline."""
    pre = q.loc[peak_t - pd.Timedelta(hours=look_back_h):peak_t]
    qmin, qpk = pre.min(), q.loc[peak_t]
    before = q.loc[:peak_t]
    start = before[before <= qmin + frac * (qpk - qmin)].index[-1]
    end = peak_t + pd.Timedelta(days=n_days)
    w = q.loc[start:end]
    t = (w.index - start).total_seconds().to_numpy()
    base = np.interp(t, [0, (end - start).total_seconds()], [w.iloc[0], w.iloc[-1]])
    vol = _trapz(w.values - base, t)
    return start, end, float(vol), float(w.iloc[0]), float(w.iloc[-1])


def compute(n_days=N_DAYS):
    q, r, cfs = load()
    rows = []
    for lab, name, t0, t1 in c.EVENTS:
        rw = r.loc[t0:t1]
        wh = rw.idxmax()
        qw = q.loc[t0:t1]
        pk = qw.q_cms.idxmax()
        row = dict(storm=lab, name=name, window_start=t0, window_end=t1, rain_mm=rw.sum(), wettest_hr_start=wh, wettest_hr_mm=rw.max(),
                   peak_time=pk, peak_cms=qw.q_cms.max(), peak_cfs=qw.q_cfs.max())
        if lab == "1":      # USGS outage covers the rise and the peak: nothing honest to report beyond the first reading
            row.update(peak_recorded=False, lag_h=np.nan, runoff_mm=np.nan, runoff_m3=np.nan, runoff_ratio=np.nan)
        else:
            s, e, vol, qs, qe = direct_runoff(q.q_cms, pk, n_days)
            row.update(peak_recorded=True, lag_h=(pk - wh) / pd.Timedelta(hours=1) - 0.5,     # to the middle of the wettest hour
                       rise_start=s, direct_end=e, base_start_cms=qs, base_end_cms=qe,
                       runoff_m3=vol, runoff_mm=vol / (AREA_KM2 * 1e6) * 1000)
            row["runoff_ratio"] = row["runoff_mm"] / row["rain_mm"]
        rows.append(row)
    ev = pd.DataFrame(rows)
    year = dict(rain_mm=float(r.sum()), wettest_hr_mm=float(r.max()), wettest_hr_start=r.idxmax(),
                flow_m3=float((cfs * CFS_TO_CMS * 86400).sum()), peak_cms=float(q.q_cms.max()), peak_cfs=float(q.q_cfs.max()),
                peak_time=q.q_cms.idxmax())
    year["flow_mm"] = year["flow_m3"] / (AREA_KM2 * 1e6) * 1000
    year["flow_acft"] = year["flow_m3"] / ACFT
    year["ratio"] = year["flow_mm"] / year["rain_mm"]
    return ev, year


if __name__ == "__main__":
    pd.set_option("display.width", 220, "display.max_columns", 30)
    ev, yr = compute()
    print(f"A = {AREA_MI2:.1f} mi2 -> N = {N_DAYS:.2f} days = {N_DAYS*24:.0f} h")
    print(ev[["storm", "name", "rain_mm", "wettest_hr_start", "wettest_hr_mm", "peak_time", "peak_cms", "peak_cfs", "lag_h"]].round(2).to_string())
    print(ev[["storm", "rise_start", "direct_end", "base_start_cms", "base_end_cms", "runoff_m3", "runoff_mm", "runoff_ratio"]].round(3).to_string())
    print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in yr.items()})
    ro = ev.runoff_mm.sum()
    print(f"sum of direct runoff, storms 2-5: {ro:.2f} mm = {ro/yr['flow_mm']*100:.0f}% of the year's flow ({yr['flow_mm']:.1f} mm)")
    print("storm-2 share of year's flow: %.0f%%" % (ev.runoff_mm[1] / yr["flow_mm"] * 100))
    print("\nsensitivity of runoff ratio to the end-of-direct-runoff rule (days after peak):")
    for nd in (1.5, N_DAYS, 3.5, 5.0):
        e2, _ = compute(nd)
        print(f"  N={nd:.2f}", e2.runoff_ratio.round(3).tolist(), "storm-2 share %.0f%%" % (e2.runoff_mm[1] / yr["flow_mm"] * 100))
    out = ev.copy()
    out.to_csv(c.OUT / "event_table_wy2017.csv", index=False)
    print("wrote", c.OUT / "event_table_wy2017.csv")
