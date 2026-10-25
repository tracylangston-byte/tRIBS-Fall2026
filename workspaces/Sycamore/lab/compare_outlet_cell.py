# ---- Compare Sycamore runs with the gauge, using the ROUTED OUTLET FLOW (version 2 of the compare cell) ----------
# What changed from compare_spinup_cell.py:
#   * The model flow is now Qstrm from results/{RESULTS_NAME}_Outlet.qout (the flow that has been routed down the
#     channels, with channel loss applied). The old cell used Srf_m3_s from the .mrf, which is runoff generated
#     BEFORE routing: channel loss cannot change it, and it peaks hours before the gauge.
#   * Srf is still shown, as a thin dashed line (and in the table, marked "unrouted"), for reference only.
# Read-only: this cell writes nothing.
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# Each run: its RESULTS_NAME, and the date and time it started from (Arizona time).
RUNS = {
    'Sycamore_jan':                '2017-01-18 00:00',   # snow off, original bedrock, water table 0.95 of bedrock
    'Sycamore_jan_snow':           '2017-01-18 00:00',   # snow on
    'Sycamore_jan_snow_bed1':      '2017-01-18 00:00',   # bedrock +1 m (water table file was stale)
    'Sycamore_jan_snow_bed1gw':    '2017-01-18 00:00',   # bedrock +1 m, water table recomputed
    'Sycamore_octjan_snow_bed1gw': '2016-10-01 00:00',   # started Oct 1 (spin-up)
    'Sycamore_jan_ks0p3':          '2017-01-18 00:00',   # Ks x 0.3, water table 0.82 of bedrock, channel loss ON (70 mm/hr)
    'Sycamore_jan_ks0p3_nocl':     '2017-01-18 00:00',   # same, channel loss OFF
}
WINDOW_START, WINDOW_END = '2017-01-18 00:00', '2017-01-25 00:00'    # the January storm week
SECOND_STORM_FROM = '2017-01-23 12:00'     # the table also reports the biggest flow after this time
OBS_FILE = '../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv'


def read_qout(path, start):
    """Routed outlet flow (m3/s) as a time series. tRIBS 6.0.0 writes a comma-delimited file with a header line."""
    with open(path) as fh:
        fh.readline()
        second = fh.readline()
    if ',' in second:
        raw = pd.read_csv(path, sep=',', header=0)
    else:
        raw = pd.read_csv(path, sep=r'\s+', skiprows=1, header=None)
    raw = raw.iloc[:, :2].copy()
    raw.columns = ['Time_hr', 'Q']
    raw = raw.apply(pd.to_numeric, errors='coerce').dropna()
    return pd.Series(raw['Q'].to_numpy(), index=pd.Timestamp(start) + pd.to_timedelta(raw['Time_hr'].to_numpy(), unit='h'))


def load_run(name, start):
    mrf_files = sorted(glob.glob(f'**/{name}[0-9]*.mrf', recursive=True))
    if not mrf_files:
        raise FileNotFoundError(f'No .mrf file found for {name!r}. Files found: {glob.glob("**/*.mrf", recursive=True)}')
    sim = pd.read_csv(mrf_files[-1])
    sim = sim[sim['RainMin_mm_hr'] < 9000]                     # drop the unfilled last row
    sim.index = pd.Timestamp(start) + pd.to_timedelta(sim['Time_hr'], unit='h')
    qout_files = glob.glob(f'**/{name}_Outlet.qout', recursive=True)
    if len(qout_files) != 1:
        raise FileNotFoundError(f'Expected exactly one {name}_Outlet.qout, found {qout_files}. '
                                f'All .qout files: {sorted(glob.glob("**/*.qout", recursive=True))}')
    return sim, read_qout(qout_files[0], start)


def volume_pieces(q):
    """Water volume (m3) of each time step: trapezoid rule on a flow series (m3/s) indexed by time.
    Each value is labeled with the time its step STARTS. Steps that touch a blank are skipped."""
    t = q.index.values.astype('datetime64[s]').astype('int64')
    v = q.values.astype(float)
    good = ~np.isnan(v[1:]) & ~np.isnan(v[:-1])
    return pd.Series(np.where(good, 0.5 * (v[1:] + v[:-1]) * np.diff(t), 0.0), index=q.index[:-1])


def volume(q, a=None, b=None):
    """Total volume (m3) from time a up to (not including) time b. Leave a and b out for the whole series."""
    p = volume_pieces(q)
    keep = np.ones(len(p), dtype=bool)
    if a is not None:
        keep &= p.index >= pd.Timestamp(a)
    if b is not None:
        keep &= p.index < pd.Timestamp(b)
    return float(p[keep].sum())


loaded = {name: load_run(name, start) for name, start in RUNS.items()}
runs = {name: sim for name, (sim, q) in loaded.items()}          # the .mrf tables (basin averages, hourly)
qouts = {name: q for name, (sim, q) in loaded.items()}           # the routed outlet flow (m3/s)
for name in runs:
    q = qouts[name]
    print(f'{name}: outlet flow {q.index[0]:%b %d %Y %H:%M} to {q.index[-1]:%b %d %Y %H:%M} ({len(q)} rows); '
          f'.mrf {len(runs[name])} hourly rows')

obs = pd.read_csv(OBS_FILE, parse_dates=['datetime']).set_index('datetime')
q_obs = obs['q_cms'].loc[min(r.index[0] for r in runs.values()):max(r.index[-1] for r in runs.values())]
q_obs_week = q_obs.loc[WINDOW_START:WINDOW_END]

# ---- Table 1: the January storm week, every run vs the gauge ---------------------------------------------
t_gauge = q_obs_week.idxmax()
gauge_peak = q_obs_week.max()
gauge_vol = volume(q_obs, WINDOW_START, WINDOW_END) / 1e6
rows = {'gauge': {'peak (m3/s)': gauge_peak, 'peak vs gauge (h)': 0.0,
                  'volume (M m3)': gauge_vol,
                  'biggest flow after 2nd-storm time (m3/s)': q_obs.loc[SECOND_STORM_FROM:WINDOW_END].max()}}
for name in runs:
    q = qouts[name].loc[WINDOW_START:WINDOW_END]
    week = runs[name].loc[WINDOW_START:WINDOW_END]
    if q.empty or week.empty:
        print(f'(skipping {name} in the January table: the run does not cover {WINDOW_START} to {WINDOW_END})')
        continue
    vol = volume(q, WINDOW_START, WINDOW_END) / 1e6
    rows[name] = {'peak (m3/s)': q.max(),
                  'peak vs gauge (h)': (q.idxmax() - t_gauge).total_seconds() / 3600,
                  'peak error vs gauge (%)': (q.max() / gauge_peak - 1) * 100,
                  'volume (M m3)': vol,
                  'volume error vs gauge (%)': (vol / gauge_vol - 1) * 100,
                  'biggest flow after 2nd-storm time (m3/s)': q.loc[SECOND_STORM_FROM:].max(),
                  'channel loss (M m3)': week['ChannelPerc_m3'].sum() / 1e6,
                  'UNROUTED Srf peak (m3/s)': week['Srf_m3_s'].max(),
                  'UNROUTED Srf volume (M m3)': volume(week['Srf_m3_s'], WINDOW_START, WINDOW_END) / 1e6,
                  'water table on Jan 18 (mm)': week['MDGW_mm'].iloc[0],
                  'water table, deepest in week (mm)': week['MDGW_mm'].max(),
                  'water table, shallowest in week (mm)': week['MDGW_mm'].min(),
                  'snow max (cm)': week['AvSWE_cm'].max()}
print('\nJANUARY STORM WEEK (Jan 18 to Jan 25), model flow = routed outlet discharge from the .qout file')
print(pd.DataFrame(rows).to_string(float_format=lambda x: f'{x:.2f}'))      # one column per run
print('Rows marked UNROUTED are the .mrf Srf column: runoff generated on the hillslopes before routing. Reference only.')
print('water table = mean depth to the water table below the surface (mm); a bigger number is deeper.')

# ---- Table 2: the spin-up months, model vs gauge (only runs that start before Jan 18) --------------------------------
for name, sim in runs.items():
    if sim.index[0] >= pd.Timestamp(WINDOW_START):
        continue
    month_model = volume_pieces(qouts[name]).resample('MS').sum() / 1e6
    month_srf = volume_pieces(sim['Srf_m3_s']).resample('MS').sum() / 1e6
    month_gauge = volume_pieces(q_obs).resample('MS').sum() / 1e6
    month_rain = sim['MAP_mm_hr'].resample('MS').sum()                # hourly mm/hr x 1 hour = mm
    month_table = pd.DataFrame({'rain (mm)': month_rain, 'gauge flow (M m3)': month_gauge,
                                'model outlet flow (M m3)': month_model, 'unrouted Srf (M m3)': month_srf})
    month_table.index = month_table.index.strftime('%b %Y')
    month_table.loc['total'] = month_table.sum()
    month_table['model / gauge'] = month_table['model outlet flow (M m3)'] / month_table['gauge flow (M m3)']
    print(f'\nMONTH BY MONTH, {name} (last hour of the run: {sim.index[-1]:%b %d %H:%M})')
    print(month_table.to_string(float_format=lambda x: f'{x:.2f}'))
    wt = sim['MDGW_mm']
    print(f'water table: starts {wt.iloc[0]:.0f} mm, is {wt.asof(pd.Timestamp(WINDOW_START)):.0f} mm on Jan 18, '
          f'deepest {wt.max():.0f} mm, shallowest {wt.min():.0f} mm.')
print('\n(The gauge total skips any blank stretch, such as the Jan 1 outage, so it can only be too small.)')

# ---- Pictures: rain, flow, water-table depth ----------------------------------------------------------------------
palette = ['tab:orange', 'tab:red', 'tab:green', 'tab:purple', 'tab:brown', 'tab:cyan', 'tab:pink', 'tab:olive', 'tab:gray']
colors = dict(zip(runs, palette))


def draw(start, end, title, date_format, shade=False):
    fig, (ax_r, ax_q, ax_w) = plt.subplots(3, 1, figsize=(11, 9), sharex=True, gridspec_kw={'height_ratios': [1, 3, 2]})
    longest = max(runs.values(), key=len)
    rain = longest.loc[start:end]
    ax_r.bar(rain.index, rain['MAP_mm_hr'], width=-1/24, color='tab:blue', align='edge')   # each bar covers the hour ENDING at its time
    ax_r.invert_yaxis()
    ax_r.set_ylabel('Rain\n(mm/hr)')
    ax_r.set_title(title)
    qo = q_obs.loc[start:end]
    ax_q.plot(qo.index, qo.values, color='black', lw=1.5, label='gauge (15 min)')
    for name, sim in runs.items():
        part = sim.loc[start:end]
        q = qouts[name].loc[start:end]
        if part.empty:
            continue
        ax_q.plot(q.index, q.values, lw=1.6, color=colors[name], label=name)                       # routed outlet flow
        ax_q.plot(part.index, part['Srf_m3_s'], lw=0.8, ls='--', color=colors[name], alpha=0.8)   # unrouted Srf, same color
        ax_w.plot(part.index, part['MDGW_mm'], lw=1.6, color=colors[name], label=name)
        ax_w.plot(part.index[0], part['MDGW_mm'].iloc[0], 'o', color=colors[name])           # dot = where the run started
    ax_q.plot([], [], lw=0.8, ls='--', color='gray', label='thin dashed = Srf, unrouted runoff (same colors)')
    if shade:                                                                                # mark the January week
        for ax in (ax_r, ax_q, ax_w):
            ax.axvspan(pd.Timestamp(WINDOW_START), pd.Timestamp(WINDOW_END), color='gold', alpha=0.25)
    ax_q.set_ylabel('Flow (m3/s)')
    ax_q.legend(fontsize=8)
    ax_q.grid(alpha=0.3)
    ax_w.set_ylabel('Mean depth to\nwater table (mm)')
    ax_w.invert_yaxis()
    ax_w.grid(alpha=0.3)
    ax_w.xaxis.set_major_formatter(mdates.DateFormatter(date_format))
    plt.tight_layout()
    plt.show()


first_start = min(r.index[0] for r in runs.values())
draw(first_start, WINDOW_END, 'Whole spin-up run vs USGS 09510200 gauge (shaded = the January storm week)',
     '%b %d', shade=True)
draw(WINDOW_START, WINDOW_END, 'January storm week: routed outlet flow (solid) and unrouted Srf (dashed) vs USGS 09510200 gauge',
     '%b %d\n%H:%M')
