SRP figure scripts: Lower Sycamore Wash, water year 2017

Setup
  Copy these .py files into workspaces/Sycamore/lab/ (they find init_data/ above them).
  Put Sycamore_boundary.shp/.shx/.dbf/.prj in lab/data/preprocessing/.
  pip install pyshp   (the other packages are numpy, pandas, matplotlib, scipy, pillow)

Run, in this order, from the lab folder
  python load_radar.py --rain-dir <folder of radar_MMDDYYYYHH.asc files> --label start
        --label start = the 2026-10-09 build (hour in the name = START of the rain hour)
        --label end   = the 2026-10-13 end-of-hour relabel
        writes out/basin_mean_hourly_rain_wy2017.csv and out/radar_wy2017_cache.npz
  python plot1_rain_and_flow.py     -> out/fig1_rain_and_flow.png
  python plot2_landcover.py         -> out/fig2_landcover.png
  python plot3_soils.py             -> out/fig3_soils.png
  python plot4_rainfall.py          -> out/fig4_rainfall.png
  python plot5_event_table.py       -> out/fig5_event_table.png and out/event_table_wy2017.csv
  python plot6_data_inventory.py    -> out/fig6_data_inventory.png and out/data_inventory.csv
  python plot7_fire_history.py      -> out/fig7_fire_history.png
  python plot8_year_choices.py      -> out/fig8_year_choices.png and out/year_screen_wy2016_2026.csv

Notes
  Everything reads from the repo's init_data/ and never writes there.
  Figures 2, 3, 6, 7 and 8 need no radar. Figures 1, 4 and 5 need load_radar.py to have run first.
  event_table.py holds the event numbers (rain, peak flow, lag, direct runoff, runoff ratio); run it alone to
  print them and the sensitivity check. The five storm windows live in srp_common.py (EVENTS).
  plot6_data_inventory.py has its content typed in from data_sources.md; edit the ROWS list there to update it.
  plot7_fire_history.py has its content typed in from the project doc "Choosing the Water Year" (fire section);
  edit the YEARS list and the numbers there. Fire acres are as reported and are not clipped to the basin.
  plot8_year_choices.py has its content typed in from the same doc (Step 3 flow screen table and Step 6);
  edit the ROWS list there. Peaks are counted on USGS daily means, so a peak is not the same as a storm.
