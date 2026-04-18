"""
Regenerate all figures for CTH_ET manuscript using updated 40yr CLM5 data.
Data: h0 files averaged over last 10 years (30yr spinup + 10yr analysis).
All 8 cases, stored in G:/Hangkai/CTH_ET project/Final_data/
"""
import netCDF4
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
from matplotlib.colors import ListedColormap
import sys
sys.stdout.reconfigure(encoding='utf-8')

OUT_DIR = 'C:/Users/hyou34/Documents/CTH_ET/figures_new_version/'
FINAL_DATA = 'G:/Hangkai/CTH_ET project/Final_data/'

# ============================================================
# 1. Load all data
# ============================================================
scenario_mapping = {
    'CLM Default': 'CLM Default.nc',
    'Default 0.8': 'Default 0.8.nc',
    'Default 0.9': 'Default 0.9.nc',
    'Default 1.1': 'Default 1.1.nc',
    'Default 1.2': 'Default 1.2.nc',
    'GEDI Mean': 'GEDI Mean.nc',
    'GEDI Max': 'GEDI Max.nc',
    'GEDI Median': 'GEDI Median.nc',
}

variable_names = ['FCEV', 'FCTR', 'FGEV']

data_all = {}
pft_meta = {}

for scenario, fname in scenario_mapping.items():
    data_all[scenario] = {}
    nc_file = os.path.join(FINAL_DATA, fname)
    with netCDF4.Dataset(nc_file, 'r') as ds:
        for var in variable_names:
            raw = ds.variables[var][:]
            # Handle both (pft,) and (time, pft) shapes
            if raw.ndim == 2:
                data_all[scenario][var] = np.mean(raw, axis=0)
            else:
                data_all[scenario][var] = raw
        # Compute ET
        data_all[scenario]['ET'] = data_all[scenario]['FCEV'] + data_all[scenario]['FCTR'] + data_all[scenario]['FGEV']

        pft_meta[scenario] = {
            'pfts1d_itype_veg': ds.variables['pfts1d_itype_veg'][:],
            'pfts1d_ixy': ds.variables['pfts1d_ixy'][:],
            'pfts1d_jxy': ds.variables['pfts1d_jxy'][:],
            'pfts1d_wtgcell': ds.variables['pfts1d_wtgcell'][:],
            'area': ds.variables['area'][:],
        }

# Use area from CLM Default
area = pft_meta['CLM Default']['area']
grid_shape = area.shape

print(f'Data loaded. Grid shape: {grid_shape}')

# ============================================================
# 2. PFT definitions
# ============================================================
pft_names = {
    "NET Temp": 1, "NET Bor": 2, "NDT Bor": 3,
    "BET Trop": 4, "BET Temp": 5,
    "BDT Trop": 6, "BDT Temp": 7, "BDT Bor": 8,
}

# ============================================================
# 3. Grid PFT-level data to 2D
# ============================================================
annual_data_2d = {var: {} for var in variable_names + ['ET']}

for scenario in scenario_mapping:
    meta = pft_meta[scenario]
    for variable in variable_names + ['ET']:
        grid_2d = np.zeros(grid_shape)
        values = data_all[scenario][variable] * meta['pfts1d_wtgcell']

        for pft_name, pft_id in pft_names.items():
            pft_idx = np.where(meta['pfts1d_itype_veg'] == pft_id)[0]
            for i, idx in enumerate(pft_idx):
                row = int(meta['pfts1d_jxy'][idx]) - 1
                col = int(meta['pfts1d_ixy'][idx]) - 1
                grid_2d[row, col] += values[idx]

        annual_data_2d[variable][scenario] = grid_2d

# Get lat/lon from any file
with netCDF4.Dataset(os.path.join(FINAL_DATA, 'CLM Default.nc'), 'r') as ds:
    if 'LATIXY' in ds.variables:
        lat2d = ds.variables['LATIXY'][:]
        lon2d = ds.variables['LONGXY'][:] - 180
    elif 'lat' in ds.variables:
        lat1d = ds.variables['lat'][:]
        lon1d = ds.variables['lon'][:] - 180
        lon2d, lat2d = np.meshgrid(lon1d, lat1d)

print('2D grids computed.')

# ============================================================
# 4. Percentage comparison analysis
# ============================================================
gedi_scenarios = ['GEDI Mean', 'GEDI Max', 'GEDI Median']
combinations = [
    ('GEDI Mean', 'CLM Default'),
    ('GEDI Max', 'CLM Default'),
    ('GEDI Max', 'GEDI Mean'),
    ('GEDI Median', 'GEDI Mean')
]
variables_compare = ['FCEV', 'FCTR', 'FGEV', 'ET']

# Area-weighted percentage differences
weighted_means = {}
for scenario in scenario_mapping:
    weighted_means[scenario] = {}
    for var in variables_compare:
        data = annual_data_2d[var][scenario]
        weighted_means[scenario][var] = np.nansum(data * area) / np.nansum(area)

print('\n===== PERCENTAGE COMPARISON (NEW 40yr data, last 10yr avg) =====')
results_pct = {}
for s1, s2 in combinations:
    results_pct[(s1, s2)] = {}
    print(f'  {s1} vs {s2}:')
    for var in variables_compare:
        pct = 100 * (weighted_means[s1][var] - weighted_means[s2][var]) / weighted_means[s2][var]
        results_pct[(s1, s2)][var] = pct
        print(f'    {var}: {pct:+.5f}%')

# Also show sensitivity scenarios
sens_combos = [
    ('Default 0.8', 'CLM Default'),
    ('Default 0.9', 'CLM Default'),
    ('Default 1.1', 'CLM Default'),
    ('Default 1.2', 'CLM Default'),
]
print('\n===== SENSITIVITY SCENARIOS =====')
for s1, s2 in sens_combos:
    print(f'  {s1} vs {s2}:')
    for var in variables_compare:
        pct = 100 * (weighted_means[s1][var] - weighted_means[s2][var]) / weighted_means[s2][var]
        print(f'    {var}: {pct:+.5f}%')

# ============================================================
# 5. Comparison with old results
# ============================================================
old_results = {
    ('GEDI Mean', 'CLM Default'): {'FCEV': -1.03792, 'FCTR': 2.78438, 'FGEV': -1.89215, 'ET': 1.17332},
    ('GEDI Max', 'CLM Default'):  {'FCEV': 2.51836, 'FCTR': -2.10603, 'FGEV': 1.46937, 'ET': -0.448814},
    ('GEDI Max', 'GEDI Mean'):    {'FCEV': 3.59358, 'FCTR': -4.75794, 'FGEV': 3.42635, 'ET': -1.60332},
    ('GEDI Median', 'GEDI Mean'): {'FCEV': 0.235559, 'FCTR': -0.389009, 'FGEV': 0.345889, 'ET': -0.136397},
}

print('\n===== OLD vs NEW COMPARISON =====')
for combo in combinations:
    print(f'  {combo[0]} vs {combo[1]}:')
    for var in variables_compare:
        new_v = results_pct[combo][var]
        old_v = old_results[combo][var]
        print(f'    {var}: old={old_v:+.5f}%  new={new_v:+.5f}%  delta={new_v-old_v:+.5f}%')

# ============================================================
# 6. Figure 4: ET comparison maps (GEDI scenarios)
# ============================================================
print('\nGenerating Figure 4...')
projection = ccrs.PlateCarree(central_longitude=180)

turbo_cmap = plt.get_cmap('coolwarm')
colors_cmap = turbo_cmap(np.linspace(0, 1, 200))
colors_cmap[95:105] = (1, 1, 1, 1)
white_center_cmap = ListedColormap(colors_cmap)

variables_plot = ['ET', 'FCEV', 'FCTR', 'FGEV']
ylabel_map = {'ET': 'ET', 'FCEV': 'Canopy EV', 'FCTR': 'Canopy TR', 'FGEV': 'Ground EV'}

fig, axes = plt.subplots(len(variables_plot), len(combinations), figsize=(22, 18),
                         subplot_kw=dict(projection=projection), dpi=300)

for j, variable in enumerate(variables_plot):
    for i, (s1, s2) in enumerate(combinations):
        ax = axes[j, i]

        data1 = annual_data_2d[variable][s1] if variable != 'ET' else \
                annual_data_2d['FCEV'][s1] + annual_data_2d['FCTR'][s1] + annual_data_2d['FGEV'][s1]
        data2 = annual_data_2d[variable][s2] if variable != 'ET' else \
                annual_data_2d['FCEV'][s2] + annual_data_2d['FCTR'][s2] + annual_data_2d['FGEV'][s2]

        base = data2.copy()
        base[base == 0] = np.nan
        plot_data = np.ma.masked_invalid((data1 - base) / base * 100)
        levels = np.linspace(-10, 10, num=11)

        cs = ax.contourf(lon2d, lat2d, plot_data, levels=levels, cmap=white_center_cmap,
                         extend='both', transform=projection)

        gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,
                          linewidth=0.5, color='gray', alpha=0.5,
                          xlocs=[-120, 0, 120], ylocs=[-90, -45, 0, 45, 90])
        gl.xlabel_style['size'] = 0
        gl.ylabel_style['size'] = 0
        gl.bottom_labels = False
        gl.right_labels = False
        ax.tick_params(labelsize=20)
        ax.coastlines()

        if j == 0:
            ax.set_title(f'{s1} vs {s2}', fontsize=18)
            gl.xlabel_style['size'] = 14
        if i == 0:
            ax.annotate(ylabel_map[variable], xy=(-0.2, 0.5), xycoords='axes fraction',
                        fontsize=18, rotation=90, va='center')
            gl.ylabel_style['size'] = 14

cbar_ax = fig.add_axes([0.39, 0.22, 0.25, 0.015])
cbar = fig.colorbar(cs, cax=cbar_ax, orientation='horizontal', ticks=levels)
cbar.set_label("Relative Difference (%)", fontsize=18)
cbar.ax.tick_params(labelsize=14)
cbar.outline.set_visible(True)
plt.subplots_adjust(hspace=-0.67, wspace=0.08)
plt.savefig(os.path.join(OUT_DIR, 'Figure_4_ET_comparison_maps.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(OUT_DIR, 'Figure_4_ET_comparison_maps.tiff'), dpi=300, bbox_inches='tight')
plt.close()
print('Figure 4 saved.')

# ============================================================
# 7. Figure 5: Bar chart by PFT
# ============================================================
print('Generating Figure 5...')

components = ['FCEV', 'FCTR', 'FGEV', 'ET']
baseline = 'CLM Default'
pft_stats = {comp: {pft: {} for pft in pft_names} for comp in components}

for comp in components:
    for pft_name, pft_id in pft_names.items():
        base_meta = pft_meta[baseline]
        base_idx = np.where(base_meta['pfts1d_itype_veg'] == pft_id)[0]
        if len(base_idx) == 0:
            continue

        base_vals = data_all[baseline][comp] * base_meta['pfts1d_wtgcell']
        base_pft = base_vals[base_idx]

        area_flat = area.flatten()
        aw = np.zeros(len(base_idx))
        for k, idx in enumerate(base_idx):
            row = int(base_meta['pfts1d_jxy'][idx]) - 1
            col = int(base_meta['pfts1d_ixy'][idx]) - 1
            fi = row * grid_shape[1] + col
            aw[k] = area_flat[fi] if fi < len(area_flat) else 1.0

        for scenario in gedi_scenarios:
            scen_meta = pft_meta[scenario]
            scen_idx = np.where(scen_meta['pfts1d_itype_veg'] == pft_id)[0]
            if len(scen_idx) == 0:
                pft_stats[comp][pft_name][scenario] = (0, 0, 0)
                continue

            scen_vals = data_all[scenario][comp] * scen_meta['pfts1d_wtgcell']
            scen_pft = scen_vals[scen_idx]

            n = min(len(base_pft), len(scen_pft))
            b, s, w = base_pft[:n], scen_pft[:n], aw[:n]
            valid = np.abs(b) > 1e-10
            if np.sum(valid) == 0:
                pft_stats[comp][pft_name][scenario] = (0, 0, 0)
                continue

            pct = (s[valid] - b[valid]) / b[valid] * 100
            wv = w[valid]
            wmean = np.average(pct, weights=wv)
            q25 = np.percentile(pct, 25)
            q75 = np.percentile(pct, 75)
            pft_stats[comp][pft_name][scenario] = (wmean, q25, q75)

fig, axes = plt.subplots(2, 2, figsize=(16, 10), dpi=300)
panel_order = ['FCEV', 'FGEV', 'FCTR', 'ET']
panel_titles = ['(a) Canopy Evaporation', '(b) Ground Evaporation',
                '(c) Canopy Transpiration', '(d) Total Evapotranspiration']
colors = {'GEDI Max': '#E8A838', 'GEDI Mean': '#E07070', 'GEDI Median': '#7CAE7A'}
pft_list = list(pft_names.keys())
x = np.arange(len(pft_list))
bar_width = 0.25

for idx, (comp, title) in enumerate(zip(panel_order, panel_titles)):
    ax = axes[idx // 2, idx % 2]
    for j, scenario in enumerate(gedi_scenarios):
        means, lo_err, hi_err = [], [], []
        for pft in pft_list:
            if scenario in pft_stats[comp][pft]:
                m, q25, q75 = pft_stats[comp][pft][scenario]
                means.append(m)
                lo_err.append(abs(m - q25))
                hi_err.append(abs(q75 - m))
            else:
                means.append(0); lo_err.append(0); hi_err.append(0)
        offset = (j - 1) * bar_width
        ax.bar(x + offset, means, bar_width, label=scenario,
               color=colors[scenario], edgecolor='black', linewidth=0.5,
               yerr=[lo_err, hi_err], capsize=2, error_kw={'linewidth': 0.8})

    ax.axhline(y=0, color='black', linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(pft_list, rotation=45, ha='right', fontsize=10)
    ax.set_ylabel('Relative Difference (%)', fontsize=12)
    ax.set_title(title, fontsize=13, fontweight='bold', loc='left')
    ax.tick_params(axis='y', labelsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if idx == 0:
        ax.legend(fontsize=10, loc='upper right', frameon=True, edgecolor='gray')

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, 'Figure_5_bar_chart.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(OUT_DIR, 'Figure_5_bar_chart.tiff'), dpi=300, bbox_inches='tight')
plt.close()
print('Figure 5 saved.')

# ============================================================
# 8. Sensitivity scenario maps (Default 0.8/0.9/1.1/1.2 vs CLM Default)
# ============================================================
print('Generating sensitivity scenario maps...')

sens_combos_map = [
    ('Default 0.8', 'CLM Default'),
    ('Default 0.9', 'CLM Default'),
    ('Default 1.1', 'CLM Default'),
    ('Default 1.2', 'CLM Default'),
]

fig3, axes3 = plt.subplots(len(variables_plot), len(sens_combos_map), figsize=(22, 18),
                           subplot_kw=dict(projection=projection), dpi=300)

for j, variable in enumerate(variables_plot):
    for i, (s1, s2) in enumerate(sens_combos_map):
        ax = axes3[j, i]

        data1 = annual_data_2d[variable][s1] if variable != 'ET' else \
                annual_data_2d['FCEV'][s1] + annual_data_2d['FCTR'][s1] + annual_data_2d['FGEV'][s1]
        data2 = annual_data_2d[variable][s2] if variable != 'ET' else \
                annual_data_2d['FCEV'][s2] + annual_data_2d['FCTR'][s2] + annual_data_2d['FGEV'][s2]

        base = data2.copy()
        base[base == 0] = np.nan
        plot_data = np.ma.masked_invalid((data1 - base) / base * 100)
        levels = np.linspace(-10, 10, num=11)

        cs = ax.contourf(lon2d, lat2d, plot_data, levels=levels, cmap=white_center_cmap,
                         extend='both', transform=projection)
        ax.coastlines()

        gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,
                          linewidth=0.5, color='gray', alpha=0.5,
                          xlocs=[-120, 0, 120], ylocs=[-90, -45, 0, 45, 90])
        gl.xlabel_style['size'] = 0
        gl.ylabel_style['size'] = 0
        gl.bottom_labels = False
        gl.right_labels = False

        if j == 0:
            ax.set_title(f'{s1} vs {s2}', fontsize=18)
            gl.xlabel_style['size'] = 14
        if i == 0:
            ax.annotate(ylabel_map[variable], xy=(-0.2, 0.5), xycoords='axes fraction',
                        fontsize=18, rotation=90, va='center')
            gl.ylabel_style['size'] = 14

cbar_ax3 = fig3.add_axes([0.39, 0.22, 0.25, 0.015])
cbar3 = fig3.colorbar(cs, cax=cbar_ax3, orientation='horizontal', ticks=levels)
cbar3.set_label("Relative Difference (%)", fontsize=18)
cbar3.ax.tick_params(labelsize=14)
plt.subplots_adjust(hspace=-0.67, wspace=0.08)
plt.savefig(os.path.join(OUT_DIR, 'Variations in Evapotranspiration in experimental scenarios.tif'), dpi=300, bbox_inches='tight')
plt.close()
print('Sensitivity map saved.')

# ============================================================
# 9. GEDI scenario maps (same as percentage notebook cell 6)
# ============================================================
print('Generating GEDI scenario maps...')

fig4, axes4 = plt.subplots(len(variables_plot), len(combinations), figsize=(22, 18),
                           subplot_kw=dict(projection=projection), dpi=300)

for j, variable in enumerate(variables_plot):
    for i, (s1, s2) in enumerate(combinations):
        ax = axes4[j, i]

        data1 = annual_data_2d[variable][s1] if variable != 'ET' else \
                annual_data_2d['FCEV'][s1] + annual_data_2d['FCTR'][s1] + annual_data_2d['FGEV'][s1]
        data2 = annual_data_2d[variable][s2] if variable != 'ET' else \
                annual_data_2d['FCEV'][s2] + annual_data_2d['FCTR'][s2] + annual_data_2d['FGEV'][s2]

        base = data2.copy()
        base[base == 0] = np.nan
        plot_data = np.ma.masked_invalid((data1 - base) / base * 100)
        levels = np.linspace(-10, 10, num=11)

        cs = ax.contourf(lon2d, lat2d, plot_data, levels=levels, cmap=white_center_cmap,
                         extend='both', transform=projection)
        ax.coastlines()

        gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,
                          linewidth=0.5, color='gray', alpha=0.5,
                          xlocs=[-120, 0, 120], ylocs=[-90, -45, 0, 45, 90])
        gl.xlabel_style['size'] = 0
        gl.ylabel_style['size'] = 0
        gl.bottom_labels = False
        gl.right_labels = False

        if j == 0:
            ax.set_title(f'{s1} vs {s2}', fontsize=18)
            gl.xlabel_style['size'] = 14
        if i == 0:
            ax.annotate(ylabel_map[variable], xy=(-0.2, 0.5), xycoords='axes fraction',
                        fontsize=18, rotation=90, va='center')
            gl.ylabel_style['size'] = 14

cbar_ax4 = fig4.add_axes([0.39, 0.22, 0.25, 0.015])
cbar4 = fig4.colorbar(cs, cax=cbar_ax4, orientation='horizontal', ticks=levels)
cbar4.set_label("Relative Difference (%)", fontsize=18)
cbar4.ax.tick_params(labelsize=14)
plt.subplots_adjust(hspace=-0.67, wspace=0.08)
plt.savefig(os.path.join(OUT_DIR, 'Variations of Evapotranspiration under Multisource Scenarios.tif'), dpi=300, bbox_inches='tight')
plt.close()
print('GEDI scenario map saved.')

print('\n===== ALL DONE =====')
print(f'All figures saved to: {OUT_DIR}')
