#!/usr/bin/env python
"""Load the MRMS hourly rain grids (radar_MMDDYYYYHH.asc) and make basin-mean hourly rain.

Use --label start for the tar built 2026-10-09 and --label end for the 2026-10-13 relabel. Label convention in the 2026-10-09 tar: HH is the START of the hour, Arizona time (MST), so
radar_0121201701.asc is the rain from 01:00 to 02:00 on 2017-01-21. (The 2026-10-13 end-of-hour relabel
renames the same grids to the hour that ends: this file would be radar_0121201702.asc.)
Amounts are mm per hour. Outputs: out/radar_wy2017_cache.npz and out/basin_mean_hourly_rain_wy2017.csv
"""
import re, sys
from pathlib import Path
import numpy as np
import pandas as pd
import srp_common as c

RAIN_DIR = Path("/home/claude/srp_plots/radar_in/mrms_wk_wy2017/rain")   # change, or pass --rain-dir
EXT = (440000.0, 472000.0, 3720000.0, 3769000.0)     # L, R, B, T from the file header
SHAPE = (49, 32)


def load(boundary=None, rain_dir=None, label="start"):
    files = sorted(Path(rain_dir or RAIN_DIR).glob("radar_*.asc"))
    pat = re.compile(r"radar_(\d{2})(\d{2})(\d{4})(\d{2})\.asc")
    times, stack = [], []
    for f in files:
        m = pat.fullmatch(f.name)
        mo, dy, yr, hr = map(int, m.groups())
        ts = pd.Timestamp(year=yr, month=mo, day=dy, hour=hr)
        times.append(ts - pd.Timedelta(hours=1) if label == "end" else ts)   # always store the START of the hour
        txt = f.read_text().split()
        hdr = dict(zip([t.lower() for t in txt[0:12:2]], map(float, txt[1:12:2])))
        assert (hdr["ncols"], hdr["nrows"], hdr["xllcorner"], hdr["yllcorner"], hdr["cellsize"]) == (32, 49, 440000, 3720000, 1000), (f.name, hdr)
        stack.append(np.array(txt[12:], dtype=np.float32).reshape(SHAPE))
    t = pd.DatetimeIndex(times)
    order = np.argsort(t)
    t, rain = t[order], np.stack(stack)[order]
    assert len(t) == 8760 and t.is_unique and (np.diff(t.values) == np.timedelta64(1, "h")).all()
    rings = c.load_boundary(boundary or c.BOUNDARY)
    cover = c.grid_cover(rings, EXT, SHAPE, ss=50)          # fraction of each 1 km cell inside the basin
    return t, rain, cover


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--rain-dir", default=str(RAIN_DIR), help="folder holding radar_MMDDYYYYHH.asc (start-of-hour labels)")
    ap.add_argument("--boundary", default=str(c.BOUNDARY))
    ap.add_argument("--label", choices=["start", "end"], default="start",
                    help="hour in the file names: start of the rain hour (2026-10-09 build) or end of it (2026-10-13 relabel)")
    a = ap.parse_args()
    t, rain, cover = load(a.boundary, a.rain_dir, a.label)
    w = cover / cover.sum()
    basin = (rain * w).sum(axis=(1, 2))                      # mm/hr, area-weighted over the basin
    grid = rain.mean(axis=(1, 2))
    print("basin area from 1 km coverage: %.1f km2 (outline 424.62)" % cover.sum())
    print("start-label range:", t[0], "->", t[-1])
    print("whole-grid year total: %.1f mm (audit 550.6) | basin-mean year total: %.1f mm" % (grid.sum(), basin.sum()))
    s = pd.Series(basin, index=t); g = pd.Series(grid, index=t)
    mo = pd.DataFrame({"grid_mm": g.resample("MS").sum().round(1), "basin_mm": s.resample("MS").sum().round(1)})
    print(mo.T.to_string())
    print("wettest basin hours:\n", s.nlargest(5).round(2).to_string())
    pd.DataFrame({"interval_start_mst": t, "interval_end_mst": t + pd.Timedelta(hours=1),
                  "basin_mean_mm": basin, "whole_grid_mean_mm": grid}).to_csv(c.OUT / "basin_mean_hourly_rain_wy2017.csv", index=False)
    np.savez_compressed(c.OUT / "radar_wy2017_cache.npz", t=t.values.astype("datetime64[m]"), rain=rain, cover=cover)
    # manifest cross-check: filled hours and hourly grid means
    mpath = Path(a.rain_dir).parent / "manifest.csv"
    if mpath.exists():
        man = pd.read_csv(mpath, parse_dates=["label_local"]).set_index("label_local").reindex(t)
        print("manifest status counts:", man.status.value_counts().to_dict())
        print("max |manifest domain_mean_mm - my whole-grid mean| = %.4f mm" % np.abs(man.domain_mean_mm.values - grid).max())
