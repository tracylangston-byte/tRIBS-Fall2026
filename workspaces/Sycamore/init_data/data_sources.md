# Data sources: Lower Sycamore Wash tRIBS model

This log records every dataset used to build the Lower Sycamore Wash model (USGS gauge 09510200): where it came from, when it was pulled, its coordinate system and resolution, and how it was processed. It also serves as the data-requirements record for Task 4.

## Digital elevation model (DEM)

- **Dataset:** USGS 3DEP 1/3 arc-second DEM, tiles n34w112 and n35w112 (vintage 2026-06-01)
- **Source:** The National Map (USGS). Tile links from the National Map API listing:
  - https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/13/TIFF/historical/n34w112/USGS_13_n34w112_20260601.tif
  - https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/13/TIFF/historical/n35w112/USGS_13_n35w112_20260601.tif
- **Date pulled:** 2026-09-30
- **CRS:** source geographic NAD83; output NAD 1983 UTM Zone 12N (EPSG:26912)
- **Resolution:** 1/3 arc-second source (about 10 m); output 10 m cells
- **Processing (ArcGIS Pro):** Mosaic To New Raster (32-bit float, 1 band, MEAN); Project Raster to EPSG:26912, bilinear, 10 m x 10 m; processing extent = basin plus 5 km buffer. Output: `Sycamore_DEM.tif`
- **Extent (EPSG:26912):** left 441789.5, right 470614.7, bottom 3720273.3, top 3767962.6
- **Checks:** DEM at the gauge reads 537.55 m (expected about 537 m); maximum elevation 2329.72 m

## Gauge location (outlet point)

- **Dataset:** USGS 09510200, Sycamore Creek near Fort McDowell, AZ
- **Source:** USGS site page, 33°41′39″ N, 111°32′28″ W (NAD27)
- **Date pulled:** 2026-09-30
- **CRS:** NAD27 geographic converted to EPSG:26912 (NAD_1927_To_NAD_1983_NADCON)
- **Result:** x = 449721, y = 3728390
- **Checks:** sits on the stream; 65 m from the USGS gauge symbol and 134 m from the StreamStats outlet pin; snap_dis = 100 m

## Watershed boundary

- **StreamStats basin polygon:** pulled 2026-09-30; 164 mi² (about 425 km²); WGS 84
- **HUC-10 / HUC-12 (Watershed Boundary Dataset, WBD):**
  - **Dataset:** USGS Watershed Boundary Dataset (WBD), 2-digit hydrologic unit 15 (Lower Colorado), file geodatabase, published 2025-01-08
  - **Source:** The National Map downloader (https://apps.nationalmap.gov/downloader/#/)
    - Download link: https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/WBD/HU2/GDB/WBD_15_HU2_GDB.zip
    - How to retrieve: select Hydrography, then Watershed Boundary Dataset (WBD). Set Data Extent to HU-2 Region and File Formats to FileGDB. Zoom the map to the basin and click Search Products. Sycamore Creek is in region 15.
  - **Date pulled:** 2026-10-08
  - **Local copy:** `Sycamore_GIS\wbd_raw\[geodatabase name]` (outside the repo)
  - **Feature load date:** 2024-08-15 (LoadDate field on these features)
  - **CRS:** as delivered, geographic GCS_North_American_1983 (datum D_North_American_1983, EPSG:4269); projected to EPSG:26912 for the overlay in ArcGIS Pro
  - **Resolution:** vector polygons; the WBD is delineated at 1:24,000 scale or better
  - **Processing:** from layers WBDHU10 and WBDHU12, selected HUC10 = 1506020306 and HUC12 starting with 1506020306, then exported the selected records to `Sycamore_Map.gdb` as `huc10_sycamore_wbd` and `huc12_sycamore_wbd`
  - **Kept:** HUC-10 1506020306 Sycamore Creek, 498.72 km² (`huc10_sycamore_wbd`; the value used, chosen 2026-10-15). Its five HUC-12s (`huc12_sycamore_wbd`): 150602030601 Upper Sycamore Creek 159.99 km², 150602030602 Mud Springs-Rock Creek 39.90, 150602030603 Mesquite Wash 51.30, 150602030604 Middle Sycamore Creek 133.18, 150602030605 Lower Sycamore Creek 114.72
  - **Gauge location:** the gauge point falls in HUC-12 150602030605 (Lower Sycamore Creek), inside HUC-10 1506020306
  - **Checks:** the five HUC-12 areas sum to 499.09 km², which is 0.37 km² more than the HUC-10 area above (498.72 km²). Not explained; probably two different area methods. The HUC-10 is about 74 km² larger than the gauge's drainage area (424.8 km², 164 mi²), so it extends past the gauge. In the ArcGIS Pro overlay, the HUC-10 and the StreamStats basin differ by 0.162 km² along their shared divide.
  - **Citation:** U.S. Geological Survey (USGS), 20250108, USGS Watershed Boundary Dataset (WBD) for 2-digit Hydrologic Unit - 15 (published 20250108) FileGDB: U.S. Geological Survey (USGS).

## Observed streamflow
### Streamflow: USGS 09510200 Sycamore Creek near Fort McDowell, AZ (daily mean discharge)

- **What:** Daily mean discharge. Parameter 00060 (cubic feet per second), statistic 00003 (mean). USGS time series 5732.
- **Source:** USGS NWIS daily-values web service, tab-delimited (RDB) format.
- **URL:** https://waterservices.usgs.gov/nwis/dv/?format=rdb&sites=09510200&parameterCd=00060&statCd=00003&startDT=2015-10-01&endDT=2026-10-07
- **Pulled:** 2026-10-08. The file header says retrieved 2026-10-09 02:52:24 -04:00, which is 23:52 Arizona time on 10-08.
- **Saved as:** workspaces/Sycamore/init_data/obs/sycamore_09510200_daily_flow_WY2016_2026.rdb  (confirm this path)
- **Period / size:** 2015-10-01 to 2026-10-07 (WY2016 to WY2026, plus 7 days of WY2027). 4,006 rows.
- **Missing days:** 19, from 2026-07-17 to 2026-08-04 (WY2026). No other gaps.
- **Quality codes in file:** A 3,144 | A:e 453 | P 227 | P:e 182. A = approved, e = estimated, P = provisional. Approved through 2025-08-05. Everything from 2025-08-06 on is provisional and can be revised.
- **Gauge:** Drainage area 164 mi² (about 425 km²). USGS site page coordinates are 33°41'39" N, 111°32'28" W (NAD27). Model outlet is 449721 E, 3728390 N (EPSG 26912).
- **Units / time:** cfs, not yet converted (1 cfs = 0.0283168 m³/s). Daily mean for the Arizona calendar day (MST, no daylight saving).
- **Used for:** Water-year screening (2026-10-08). WY2017 chosen as the calibration year (2026-10-09).
- **WY2017 (2016-10-01 to 2017-09-30):** 365 of 365 days, all approved (252 A, 113 A:e). Largest daily mean 1,920 cfs on 2017-01-21. Volume about 15,000 acre-feet. Essentially no flow Oct–Dec 2016 or Jun–Sep 2017. Estimated days that matter: 2017-01-01 to 01-04 and 2017-01-25 (both event days). The rest are zero-flow days in Dec and May–Sep.
- **Limits:** Daily means understate flood peaks. The 15-minute (instantaneous) series is pulled for WY2017 (see the next entry). Upstream diversions and regulation are not checked.

### Streamflow, 15-minute: USGS 09510200 Sycamore Creek near Fort McDowell, AZ (WY2017)

- **What:** Instantaneous (15-minute or finer) discharge. Parameter 00060 (cubic feet per second). USGS time series 6802, the one the daily means in the entry above are computed from.
- **Source:** USGS NWIS instantaneous-values web service, tab-delimited (RDB) format.
- **URL pattern:** https://waterservices.usgs.gov/nwis/iv/?format=rdb&sites=09510200&parameterCd=00060&startDT=YYYY-MM-DD&endDT=YYYY-MM-DD
- **Pulled:** 2026-10-09 (Arizona time), with curl, in pieces. Requests for the whole year, for Jan–Mar and for Mar 1–15 returned HTTP 503 "service temporarily unavailable", so the year was pulled in smaller pieces.
- **Saved in:** workspaces/Sycamore/init_data/obs/
- **Period covered:** 2016-10-01 to 2017-09-30 (WY2017), except the gaps below.
- **Files and row counts (data rows, `grep -c '^USGS'`):**
  - iv_2016-10-01_to_2016-12-31.rdb: 8,340
  - iv_2017-01-01_to_2017-01-31.rdb: 3,028
  - iv_2017-02-01_to_2017-02-28.rdb: 2,715
  - iv_2017-03-01.rdb to iv_2017-03-15.rdb: one file per day, no file for 03-12 (counts not recorded yet)
  - iv_2017-03-16_to_2017-03-31.rdb: 1,536
  - iv_2017-04-01_to_2017-06-30.rdb: 8,736
  - iv_2017-07-01_to_2017-09-30.rdb: 8,832
- **Columns:** agency_cd, site_no, datetime, tz_cd, 6802_00060, 6802_00060_cd. Time zone code was MST on all 34,568 rows (checked 2026-10-15). Arizona has no daylight saving.
- **Spacing:** Readings are every 15 minutes in calm weather and as often as every 5 minutes during rises (2017-01-21 had 148 readings; a normal day has 96).
- **Gaps:** (confirmed 2026-10-09): Dec 9–13 2016 (USGS estimated zero flow, dry); Jan 1 2017 all day + Jan 2 00:00–00:30 (real outage during a storm, daily means 168/232 cfs estimated; exclude from scoring); Mar 12 2017 (failed download, daily mean 31.3 cfs; retry); three 1.25 h dry-day dropouts Nov 27, Dec 2, Dec 6.
- **Checks:** 34,568 readings, all MST, none provisional; 355 of 355 compared days agree with USGS daily means within 5%; peak 7,800 cfs (220.9 m³/s) at 2017-01-21 01:35 MST.
- **Quality codes:** A 24,620 and A:e 9,948 across all 15-minute files (checked 2026-10-15); no provisional codes.
- **Example of why this file is needed:** on 2017-01-21 the daily mean is 1,920 cfs, while the 15-minute record peaks at 7,800 cfs at 01:35 MST. Confirm this from the merged file.
**Units:** cfs in the raw files (1 cfs = 0.0283168 m3/s). The model-ready file made by make_obs_flow.py is in m3/s, Arizona time.- **Used for:** WY2017 calibration target for storm timing and peaks. The daily-mean file is used for volumes and for screening.
- **Done:** merge_inst_flow.py joined the pieces and make_obs_flow.py wrote the model-ready file (Mar 12 interpolated from the daily mean; Jan 1-2 outage excluded from scoring).
- **Open items:** retry 2017-03-12 (try format=json if rdb keeps failing). Status as of 2026-10-15: No 2017-03-12 file exists in the folder (checked 2026-10-15).

## Rainfall
### Radar rainfall, pre-Oct 2020 (candidate)

- Source: Iowa Environmental Mesonet MTArchive (NCEP MRMS files)
- URL pattern: https://mtarchive.geol.iastate.edu/YYYY/MM/DD/mrms/ncep/GaugeCorr_QPE_01H/GaugeCorr_QPE_01H_00.00_YYYYMMDD-HH0000.grib2.gz
- Product: Local Gauge Bias Corrected radar QPE, 1-hour, mm (older MRMS algorithm,
  before v12 on 2020-10-14; not the same as the AWS Pass 2 product)
- Accessed: 2026-10-08 (local)
- Tested: one file only (2019-07-20 21:00 UTC), saved outside the repo in ~/mrms_test.
  Grid is 7000 x 3500, 0.01 deg, lat/lon. Valid time matches the file name.
  -3 is the only negative value and means no coverage. Cells at the gauge read 0.0 (dry hour).
- Downloaded for WY2017 only (see the WY2017 entry below).
- Open items:
  - Time label (start vs end of the hour): the archive stamp is assumed to mark the END of the accumulation hour (the GRIB metadata does not say). tRIBS reads the hour in a rain file name as the past hour (Josh Cederstrom, 2026-10-13), so the WY2017 files are labeled by the end of the hour (see the WY2017 entry below). Still to confirm against a gauge on a real storm.
  - GaugeCorr folder absent on 2015-01-15, present on 2015-07-20; start date unknown, so use WY2016 onward.
  - 2014-10-15 returned a 404.
  - Spot-checked 2015-07-20, 2017-08-15 and 2019-07-20: all had 24 of 24 hours. Not a completeness guarantee.
  - IEM has said it will cut retention of its hourly zip files; the per-date folders above were still up when checked.

### Radar rainfall, Oct 2020 onward (candidate, untested)

- Source: NOAA AWS bucket noaa-mrms-pds (MultiSensor_QPE_01H_Pass2)
- Archive starts 2020-10-14 (MRMS v12.0). The first usable date should be confirmed from the bucket listing.
- build_mrms_radar_rainfall_tribs6.py has not been run against the real bucket.

### Radar rainfall, hourly gridded (WY2017)

- Product: NCEP MRMS hourly rain accumulation `GaugeCorr_QPE_01H_00.00` (mm), from the Iowa State Mesonet MTArchive copy. This is the older, locally gauge-corrected algorithm. SRP uses the current MRMS product (available from 2020-10-14), and the two have not been spliced.
- URL pattern: https://mtarchive.geol.iastate.edu/YYYY/MM/DD/mrms/ncep/GaugeCorr_QPE_01H/GaugeCorr_QPE_01H_00.00_YYYYMMDD-HH0000.grib2.gz
- Pulled 2026-10-09 with build_mrms_radar_rainfall_tribs6.py v1.1: `download`, then `all` (process + audit). Re-labeled 2026-10-13 by re-running `process` and `audit` only (no re-download) with `--stamp-convention end`. Settings of the current files: --source iem --bbox -111.64 33.62 -111.31 34.06 --buffer-km 0 --start "2016-10-01 01:00" --end "2017-10-01 00:00" --stamp-convention end. Working folder: ~/mrms_wk_wy2017
- Period: 8,760 hourly files, labeled 2016-10-01 01:00 to 2017-10-01 00:00 Arizona time (MST, UTC-7). The label is the END of the hour, so radar_1001201601.asc is the rain from 00:00 to 01:00 on 10/01/2016 and the last file, radar_1001201700.asc, is 23:00 to 24:00 on 09/30/2017. Archive file stamps run 2016-10-01 08:00 to 2017-10-01 07:00 UTC.
- Grid: 32 x 49 cells of 1,000 m, EPSG:26912, x 440,000-472,000, y 3,720,000-3,769,000. Area-average resampling. Units mm/hr (an hourly accumulation in mm is numerically mm/hr), no conversion.
- Files: ~/mrms_wk_wy2017/rain/radar_MMDDYYYYHH.asc (HH = the hour that ENDS, 00-23. The hour ending at midnight is written as 00 of the next day; per the tRIBS 6.0.0 source (read 2026-10-14), tRIBS rolls hour 24 over to 00, so this is right, though not yet confirmed by a run) [copy to the model folder: path]. Also kept in the same folder: manifest.csv, download_log.csv, audit/audit_report.txt, run_all.log. The earlier start-labeled grids are in ~/mrms_wk_wy2017/rain_start_label_OLD/ (to delete after a successful tRIBS run).
- Hour-label convention: END of the hour (changed 2026-10-13). Rule from Josh Cederstrom, relayed by Tracy (secondhand; the tRIBS documentation pages checked do not state it): in tRIBS the rainfall hour in a file name means the past hour. The MRMS archive stamp is assumed to mark the end of its accumulation hour (the GRIB metadata does not say). The first build (2026-10-09) used start-of-hour labels and was re-labeled; the rain amounts did not change, only the file names. Check of 2026-10-09: FCDMC reports ending 17:00, 18:00 and 19:00 MST on 07/23/17 showed all three Sycamore gauges peaking 17:00-18:00 MST, the hour now in the file labeled 18:00. One storm, hourly resolution; rules out a one-hour shift only.
- Missing from the archive (HTTP 404): 17 of 8,760 hours (0.19%), written as all-zero grids and flagged filled_zero_missing_hour in manifest.csv. Archive stamps (UTC): 2016-10-11 16:00-17:00; 2016-11-06 01:00; 2016-12-06 17:00-23:00; 2017-04-27 12:00-13:00; 2017-05-09 18:00; 2017-05-11 21:00-22:00; 2017-05-14 21:00; 2017-07-31 18:00. None fall in January-March.
- Audit: PASS with 1 warning (the zero-filled hours), also after the 2026-10-13 relabel. All 8,760 files present and parse, one shared header, no NaN or negative values.
- Basin mean for the year: 567 mm (the value on the data inventory slide, adopted 2026-10-15; computed from the radar grids clipped to the basin outline). Whole-grid total, not basin-only: 550.6 mm (21.7 in). Monthly whole-grid means (mm): Oct 17.3, Nov 39.5, Dec 88.7, Jan 111.9, Feb 81.4, Mar 10.2, Apr 1.2, May 7.1, Jun 0.6, Jul 147.3, Aug 31.4, Sep 14.1. Wettest winter hour: 2017-01-20 21:00 (6.17 mm/hr). Wettest hour: 2017-07-23 18:00 (8.17 mm/hr). December and January differ slightly from the first build because the relabel moves the last hour of each month into the next month.
- Plausibility check against FCDMC ALERT gauges (Mt. Ord 2000, Sycamore Cr Upper 82200, Iron Dike 83000), using five by-hand historic reports (12/24/16, 12/31/16, 01/20/17, 07/14/17, 07/24/17 at 23:59 MST, from alert.fcd.maricopa.gov/php/prpt_hist.php): radar grid-mean 72-hr totals were 69-87% of the mean of the three gauges in four widespread storms; 07/14/17 was a localized storm (only 82200 recorded notable rain). Mt. Ord read ~0 on 01/20/17 (outage or unmelted snow; left out of that comparison). Comparison used the whole-grid mean, not radar at the gauge cells.
- Known limits: tRIBS 6.0.0 reads gridded rain hourly only, which smooths convective peaks. The monthly values and the 550.6 mm total use the whole grid, not a basin outline.

## Other weather
### Other weather: NLDAS-2 point stations (written 2026-10-13)

- Source: NLDAS-2 hourly primary forcing (NLDAS_FORA0125_H v2.0, 0.125 degree), NASA GES DISC, retrieved through the Giovanni time-series service with pytRIBS Met.get_nldas_point (pytRIBS 1.0.0).
- Variables: air temperature, specific humidity, surface pressure, 10 m wind (u and v), downward shortwave radiation. NLDAS rain was also downloaded but is NOT used (rain comes from the MRMS/IEM radar grids).
- Window: 2016-10-01 00:00 to 2017-09-30 23:00 Arizona time (UTC-7); requested 2016-09-30 to 2017-10-01 UTC. Timestamp shift: gmt = -8 (LABEL_SHIFT_HOURS = -1; start-versus-end-of-hour meaning NOT yet confirmed).
- Conversion: pytRIBS convert_and_write_nldas_timeseries (Pa to hPa, K to C, specific humidity to RH, 10 m wind to 2 m, UTC to local). Gaps up to 3 hours filled by linear interpolation.
- Stations (nearest-station rule, no extra lapse-rate correction; station elevation = NLDAS cell mean elevation):
  - NLDAS_1 (south-west cell): easting 450970.2, northing 3730800.7, elevation 657 m, lat 33.7161, lon -111.5292
  - NLDAS_2 (middle-east cell): easting 458992.2, northing 3741416.2, elevation 1247 m, lat 33.8122, lon -111.4431
  - NLDAS_3 (north-east cell): easting 457187.2, northing 3754435.2, elevation 1435 m, lat 33.9295, lon -111.4632
- Files: ../init_data/met/Master_Met.sdf and the met_NLDAS_*.mdf files beside it. NLDAS rain files written outside the repo (~/nldas_unused_precip/Precip_NLDAS_unused.sdf) and not used.
- tRIBS check (source read 2026-10-14): TEMPLAPSE (default -0.0065 K/m) adjusts air temperature to each node's elevation, but not pressure. Whether 3 stations are enough across 540-2,330 m is open (in the plan doc’s Josh table). Whether the weather rows follow the "past hour" rule is open (in the plan doc’s Josh table).

## Soils
### SOLUS100 soil properties and bedrock depth

- Source: SOLUS100 (Soil Landscapes of the United States, 100 m), read from the public `solus100pub` Google Cloud bucket through pytRIBS (`Soil.get_solus_grids`, `Soil.get_solus_bedrock`).
- Variables: sand, silt, clay, bulk density at 0-5, 5-15 and 15-30 cm; depth to bedrock.
- Accessed: 2026-10-13. Cell size 100 m; downloaded in a larger upright box with zero padding outside the basin.
- Processing: `Soil.run_soil_workflow(source='SOLUS')`. Ks, thetaS, thetaR, m and PsiB come from the Rosetta v3 pedotransfer function, f from a fit of Ks over 0-30 cm. tRIBS uses only the 0-5 cm Ks grid; the 5-15 and 15-30 cm Ks grids only feed the fit for f. scgrid.gdf has 7 rows (KS, TR, TS, PB, PI, FD, PO). Basin medians: Ks 11.53 mm/hr; f 0.00302 per mm, with 9.3% of cells at the 0.0001 floor.
- Used with optsoiltype = 1 (gridded). Josh says SRP uses SOLUS.
- Citation: add once confirmed (ask Josh).

## Bedrock depth and initial groundwater

- **Dataset:** SOLUS100 depth to bedrock, the same product as the soil properties (see the SOLUS100 entry above). The initial groundwater grid is derived from it.
- **Source:** public `solus100pub` Google Cloud bucket, read through pytRIBS (`Soil.get_solus_bedrock`).
- **Date pulled:** 2026-10-13
- **Resolution:** 100 m. Downloaded in a larger upright box with zero padding outside the basin, as for the soil grids.
- **Processing:**
  1. `Soil.get_solus_bedrock` converts depth from cm to m and writes `init_data/bedrock_depth_m.asc` (and `.prj`).
  2. `Soil.generate_initial_groundwater('../init_data/bedrock_depth_m.asc', 0.95, filename='../init_data/InitGW_95pct_mm.asc')` writes the initial water table depth in mm: 0.95 x bedrock depth (m) x 1000. This is the same 95% rule as the SMF model.
- **Files:** `init_data/bedrock_depth_m.asc`, `init_data/InitGW_95pct_mm.asc` (each with a `.prj`).
- **Checks (inside the basin outline, 42,449 cells):** no missing cells; bedrock depth 0.12 to 1.73 m, median 0.60 m, none at 0 m. Initial water table 114 to 1,644 mm below the surface, median 570 mm.
- **Limits:** SOLUS depth to bedrock is a modeled estimate, not a measurement. The 95% rule is an assumption that has not been tested; where bedrock is shallow (0.12 m) it gives a very shallow starting water table. The spin-up length is still open (question for Josh).
- **Used for:** starting water table (`gwaterfile`) for the WY2017 run. Not yet wired into the model; that happens when the model is assembled.
- **Citation:** same as the SOLUS100 entry; add once confirmed (ask Josh).

## Land cover
### Land cover: LANDFIRE 2016 Remap, Existing Vegetation Type (EVT)

- **Product:** LANDFIRE 2016 Remap, CONUS, 30 m. Existing Vegetation Type (EVT) builds the land cover map. Existing Vegetation Cover (EVC) and Existing Vegetation Height (EVH) set vegetation fraction (V) and height (h) in the land use table. Base map circa 2016 (imagery about 2013-2017).
- **Source:** LANDFIRE data viewer (landfire.gov); each download arrived as an emailed zip.
- **Pulled:** EVT 2026-10-13. EVC and EVH 2026-10-14.
- **Why this version:** its date fits WY2017 and predates the June 2020 Bush Fire.
- **Download box (EPSG 3857):** x -12,440,073.9 to -12,373,957.6, y 3,965,047.9 to 4,038,942.6. It is larger than the DEM. EVC and EVH reused the GeoJSON from the EVT download (`lf_download_extent.geojson`), so they cover the same area. They sit on the same 30 m grid but are one cell smaller on every side (2128 x 2321 cells versus 2130 x 2323), which does not matter because the basin is well inside.
- **Raw files:** `Sycamore_GIS\landcover_raw\` on the local computer (`.tif`, `.tif.vat.dbf`, `.tif.ovr`, `.tfw`, `.dbf`, GeoJSON of the box). The scripts read copies of the EVT, EVC and EVH `.tif` and `.tif.vat.dbf` files in `init_data/landcover_raw/`. If `.tif` files are git-ignored, those copies are not in the repo; re-download from the viewer.
- **CRS of the raw files:** NAD83 Contiguous USA Albers (EPSG 5070). Nodata 32767.
- **Processing:**
  1. ArcGIS Pro, to count classes only: Project Raster to EPSG 26912 (nearest neighbour, 10 m, snapped to the DEM), Clip Raster to the DEM extent, Extract by Mask with the StreamStats basin, Export Table. Result: 45 EVT types, 4,243,489 cells = 424.35 km2 in the basin.
  2. `lab/make_sycamore_landuse_map.py` (adapted from Josh Cederstrom's `make_landuse_map.py`): reads the raw EVT `.tif` and `.tif.vat.dbf`, resamples by nearest neighbour onto the 10 m DEM grid (`Sycamore_DEM.tif`), groups the EVT types into 8 classes, and writes `init_data/LandUse.asc` (NODATA -9999) and `init_data/LandUse_crosswalk.csv`. This is the step the model uses.
  3. `lab/summarize_landuse_cover_height.py` (2026-10-14): moves EVC and EVH onto the `LandUse.asc` grid by nearest neighbour, keeps cells inside the basin, decodes the codes (tree/shrub/herb cover % and height in m), and prints cover and height per class. It writes nothing.
  4. `lab/make_sycamore_landuse_table.py` (2026-10-14): holds the 8 class rows, checks the class IDs, value ranges and stress thresholds against the soil grids, and writes `init_data/land_use_params.ldt` (14 columns, 8 rows).
- **Classes (share of basin, ArcGIS counts inside the StreamStats polygon):** 1 woodland 46.7%, 2 desert scrub 35.4%, 3 chaparral and oak shrub 10.3%, 4 bare rock 2.9%, 5 recently burned or logged 2.6%, 6 developed 0.8%, 7 riparian and water 0.7%, 8 grassland 0.7%. Areas in km2: 198.4, 150.1, 43.6, 12.1, 11.0, 3.4, 3.0, 2.8 (classes 1 to 8).
- **Vegetation values (basin means from EVC and EVH):** V and h in m. 1 woodland 0.28, 6.3. 2 desert scrub 0.23, 0.7. 3 chaparral 0.24, 0.5. 4 bare rock 0.05 (all cells), no heights. 5 burned or logged 0.29, 4.7. 6 developed none. 7 riparian and water 0.31, 4.3. 8 grassland 0.27, 0.3. The table does not use the raw LANDFIRE means for classes 4 and 6: both get V = 0.01 and h = 0.01 (Josh's barren and developed rows).
- **Land use table:** the other columns come from Josh Cederstrom's South Mountain table (shrub, grassland, developed and barren rows). Woodland and burned LAI is 1.8 (Grier et al. 1992). Stress thresholds are 0.33 (evaporation) and 0.30 (transpiration) for every class. They must sit above the largest thetaR (0.09) and below the smallest thetaS (0.36) near the basin, so South Mountain's 0.37 would stop tRIBS. RZD_m is left at 9999.99, which tRIBS reads as its 1 m default (never deeper than bedrock). The v6.0.0 table has no a or b1 columns. Final values: see "Building the Land Use Table: Lower Sycamore Wash Workflow".
- **Limits:** the 30 m source is resampled to 10 m, so cells are blocky and add no detail. The grouping is by vegetation type and has not been tested for hydrologic effect. The 2.6% mapped as recently burned or logged in the 2016 base map has not been checked against fire history before 2017. LANDFIRE cover counts only the dominant lifeform of each pixel, so V under woodland is probably higher than 0.28. Chaparral is mapped at 24% cover and 0.5 m, but the LANDFIRE description of Mogollon chaparral says 40-80% cover and 1-2.5 m shrubs; its V is the first value to test in calibration.
- **Open questions for Josh (optional, none blocks a run):** keep chaparral V at 0.24 or raise it? Would he vary the stress thresholds by class? Is there a source for his grassland and barren rows? Answered without him: vegetation values come from LANDFIRE (above), and land cover does not have to match the DEM cell size because tRIBS resamples each grid onto the mesh.

## Fire history
- **Status:** Pulled 2026-10-09 from Forest Service reports and news articles (project doc "Choosing the Water Year: Data Pulled and How We Screened It", Step 5). Acres are as reported and sources disagree. Nothing has been clipped to the model basin boundary.
- **Used for:** Choosing the water year. The June 2020 Bush Fire burned part of this basin, so WY2021 and later are post-fire. WY2016-2019, plus WY2020 from October to early June, are pre-fire. WY2017 was chosen because it is pre-fire, so one set of land-cover and soil conditions describes the whole year.
- **Why it matters:** Burned areas can shed water faster and may produce more runoff, so parameters calibrated on a post-fire year might not carry over to a pre-fire basin.
- **Bush Fire (June 13 to about July 6, 2020):** 195,107 acres in the USFS burned-area memo (Jul 21, 2020); 193,455 acres on Wikipedia. Tonto National Forest. Burn severity for the whole fire: unburned 9.4%, low 41.3%, moderate 44.7%, high 4.6%.
- **Bush Fire in the Sycamore Creek sub-watersheds (HUC-12, acres burned of total):** Upper Sycamore 15,007 of 39,506; Middle Sycamore 11,388 of 32,885; Lower Sycamore 4,997 of 28,327. Together about 31,400 acres, roughly 30% of the 424 km² basin. Treat this as an upper bound, because part of Lower Sycamore lies below the gauge. The USFS report disagrees with itself: 15,007 ÷ 39,506 is 38%, but its table says 41.5%. A flood-model table in the same report gives 51,364 of 119,157 acres (43%) for a different set of areas. The percentages are ballpark only.
- **Later fires (as reported):**
  - Bullet Fire, Jun 4-7, 2023: about 3,355 acres (the same article also gives 3,335 and 4,178), Lower Sycamore area near SR 87.
  - Diamond Fire, started Jul 22, 2023: 5,060 acres (an earlier update says 5,101), near Bartlett Lake, southwest of Sunflower.
  - Sand Stone Fire, started Jul 25, 2024: over 27,000 acres in the article text, but the web address says nearly 13,000. Between Hwy 87 and the Verde River, near Sunflower.
  - Sycamore Fire, June 2026: nearly 12,000 acres, about 15 miles north of Globe. A different watershed.
- **Not done yet:**
  - Overlap with the model boundary is not measured. Plan: intersect the NIFC (National Interagency Fire Center) perimeter history with the basin polygon in ArcGIS Pro, and report acres and percent burned by year.
  - The 2012 Sunflower Fire is not checked. An older note put it in the upper Sycamore area at about 17,500 acres.
  - Fires in WY2025 are not checked.
  - The Diamond and Sand Stone articles place both fires near Sunflower but do not show either inside the basin. That is not the same as confirming they are outside.
  - Land cover must match 2016-17. The LandUse map (2016 base) has 2.6% mapped as recently burned or logged, still unchecked against pre-2017 fire history.
- **Sources:** Bush Fire Watershed Specialist Report, USFS (https://www.fs.usda.gov/sites/default/files/2022-07/230%20Watershed%20Specialist%20Report_BushFire.pdf); Bush Fire BAER funding request, USFS (https://www.fs.usda.gov/sites/default/files/2022-07/209%20Bush%20Fire%20BAER%20Funding%20Request_interim_2%20signed.pdf); Bush Fire, Wikipedia (https://en.wikipedia.org/wiki/Bush_Fire_%28Arizona%29); Bullet Fire, Fox10 (https://fox10phoenix.com/news/fire-burning-near-sr-87-in-tonto-national-forest-shutters-portion-of-highway); Diamond Fire, KJZZ (https://kjzz.org/content/1853142/diamond-fire-latest-wildfire-near-sunflower-arizona); Sand Stone Fire, ABC15 (https://abc15.com/news/wildfires/sand-stone-fire-shuts-down-south-sr-87-after-growing-to-nearly-13-000-acres); Sycamore Fire, KJZZ (https://www.kjzz.org/2026-06-18/sycamore-fire-the-latest-on-the-wildfire-in-the-tonto-national-forest).