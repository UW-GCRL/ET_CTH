"""
Regenerate spatial mapping figures with fixed title spacing.
- Figure "experimental scenarios": 4 rows x 4 cols (CLM Default + 3 GEDI scenarios)
- Figure "multisource scenarios": 4 rows x 4 cols (all % difference)
- Figure absolute difference: 4 rows x 4 cols (CLM Default + 3 abs diff)
- Figure absolute values: 4 rows x 4 cols (all absolute mm/month)
- Figure 4 (ET comparison): 4 rows x 4 cols with subplot labels
"""
import netCDF4
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
from matplotlib.colors import ListedColormap
from matplotlib import ticker
import sys
sys.stdout.reconfigure(encoding='utf-8')

OUT_DIR = 'C:/Users/hyou34/Documents/CTH_ET/figures_new_version/'
FINAL_DATA = 'G:/Hangkai/CTH_ET project/Final_data/'

# ============================================================
# Load data
# ============================================================
scenario_mapping_EX = {
    'CLM Default': 'CLM Default.nc',
    'GEDI Mean': 'GEDI Mean.nc',
    'GEDI Max': 'GEDI Max.nc',
    'GEDI Median': 'GEDI Median.nc',
}

variable_names = ['FCEV', 'FCTR', 'FGEV']

pft_names = {
    "NET Temperate": 1, "NET Boreal": 2, "NDT Boreal": 3,
    "BET Tropical": 4, "BET Temperate": 5,
    "BDT Tropical": 6, "BDT Temperate": 7, "BDT Boreal": 8,
}

data_EX = {}
pft_meta = {}
for scenario, fname in scenario_mapping_EX.items():
    data_EX[scenario] = {}
    nc_file = os.path.join(FINAL_DATA, fname)
    with netCDF4.Dataset(nc_file, 'r') as ds:
        for var in variable_names:
            raw = ds.variables[var][:]
            data_EX[scenario][var] = np.mean(raw, axis=0) if raw.ndim == 2 else raw
        data_EX[scenario]['ET'] = data_EX[scenario]['FCEV'] + data_EX[scenario]['FCTR'] + data_EX[scenario]['FGEV']
        pft_meta[scenario] = {
            'pfts1d_itype_veg': ds.variables['pfts1d_itype_veg'][:],
            'pfts1d_ixy': ds.variables['pfts1d_ixy'][:],
            'pfts1d_jxy': ds.variables['pfts1d_jxy'][:],
            'pfts1d_wtgcell': ds.variables['pfts1d_wtgcell'][:],
            'area': ds.variables['area'][:],
        }

area = pft_meta['CLM Default']['area']
grid_shape = area.shape

# Get lat/lon
with netCDF4.Dataset(os.path.join(FINAL_DATA, 'CLM Default.nc'), 'r') as ds:
    if 'LATIXY' in ds.variables:
        lat2d = ds.variables['LATIXY'][:]
        lon2d = ds.variables['LONGXY'][:] - 180
    elif 'lat' in ds.variables:
        lat1d = ds.variables['lat'][:]
        lon1d = ds.variables['lon'][:] - 180
        lon2d, lat2d = np.meshgrid(lon1d, lat1d)

# Grid to 2D
annual_data_EX = {var: {} for var in variable_names + ['ET']}
for scenario in scenario_mapping_EX:
    meta = pft_meta[scenario]
    for variable in variable_names + ['ET']:
        grid_2d = np.zeros(grid_shape)
        values = data_EX[scenario][variable] * meta['pfts1d_wtgcell']
        for pft_name, pft_id in pft_names.items():
            pft_idx = np.where(meta['pfts1d_itype_veg'] == pft_id)[0]
            for i, idx in enumerate(pft_idx):
                row = int(meta['pfts1d_jxy'][idx]) - 1
                col = int(meta['pfts1d_ixy'][idx]) - 1
                grid_2d[row, col] += values[idx]
        annual_data_EX[variable][scenario] = grid_2d

print('Data loaded and gridded.')

# Colormaps
turbo_cmap = plt.get_cmap('coolwarm')
colors_cmap = turbo_cmap(np.linspace(0, 1, 200))
colors_cmap[95:105] = (1, 1, 1, 1)
white_center_cmap = ListedColormap(colors_cmap)

projection = ccrs.PlateCarree(central_longitude=180)
variables = ['ET', 'FCEV', 'FCTR', 'FGEV']
ylabel_map = {'ET': 'ET', 'FCEV': 'Canopy EV', 'FCTR': 'Canopy TR', 'FGEV': 'Ground EV'}
row_labels = ['a', 'b', 'c', 'd']

scenario_order = ['CLM Default', 'GEDI Mean', 'GEDI Max', 'GEDI Median']


def add_subplot_label(ax, row_idx, col_idx):
    """Add (a.1) style label to top-left of subplot."""
    label = f'({row_labels[row_idx]}.{col_idx+1})'
    ax.text(0.02, 0.95, label, transform=ax.transAxes, fontsize=13,
            fontweight='bold', va='top', ha='left',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='gray'))


def add_gridlines(ax, show_xlabel=False, show_ylabel=False):
    """Add gridlines with consistent styling."""
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,
                      linewidth=0.5, color='gray', alpha=0.5,
                      xlocs=[-120, 0, 120], ylocs=[-90, -45, 0, 45, 90])
    gl.xlabel_style['size'] = 13 if show_xlabel else 0
    gl.ylabel_style['size'] = 13 if show_ylabel else 0
    gl.bottom_labels = False
    gl.right_labels = False
    return gl


def get_et_data(scenario):
    """Get ET as sum of components."""
    return (annual_data_EX['FCEV'][scenario] +
            annual_data_EX['FCTR'][scenario] +
            annual_data_EX['FGEV'][scenario])


# ============================================================
# Figure 1: Experimental Scenarios (4 rows x 4 cols)
# Col 0: CLM Default absolute (mm/month)
# Cols 1-3: GEDI Mean, Max, Median vs Default (%)
# ============================================================
print('Generating experimental scenarios figure...')

fig, axes = plt.subplots(4, 4, figsize=(24, 20),
                         subplot_kw=dict(projection=projection), dpi=200)

for j, variable in enumerate(variables):
    for i, scenario in enumerate(scenario_order):
        ax = axes[j, i]

        data_et = get_et_data(scenario) if variable == 'ET' else annual_data_EX[variable][scenario]

        if scenario == 'CLM Default':
            base_data = data_et.copy()
            base_data[data_et == 0] = np.nan
            plot_data = np.ma.masked_invalid(base_data * 1.0562)
            levels = np.linspace(np.nanmin(plot_data), np.nanmax(plot_data), num=21)
            mn, mx = np.nanmin(plot_data), np.nanmax(plot_data)
            mid = 0.5 * (mn + mx)
            cbar_ticks = [mn, 0.5*(mid+mn), mid, 0.5*(mid+mx), mx]
            colormap = plt.get_cmap('afmhot_r')
            title = 'CLM Default\n(mm/month)'
        else:
            base_data = get_et_data('CLM Default') if variable == 'ET' else annual_data_EX[variable]['CLM Default'].copy()
            base_data = base_data.copy()
            base_data[base_data == 0] = np.nan
            plot_data = (data_et - base_data) / base_data * 100
            plot_data = np.ma.masked_invalid(plot_data)
            levels = np.linspace(-10, 10, num=21)
            colormap = white_center_cmap
            cbar_ticks = [-10, -5, 0, 5, 10]
            title = f'{scenario}\nvs Default (%)'

        cs = ax.contourf(lon2d, lat2d, plot_data, levels=levels, cmap=colormap,
                         extend='both', transform=projection)
        ax.coastlines()
        add_gridlines(ax, show_xlabel=(j == 0), show_ylabel=(i == 0))

        cbar = plt.colorbar(cs, ax=ax, orientation='horizontal', pad=0.005, shrink=0.6, ticks=cbar_ticks)
        cbar.formatter = ticker.FormatStrFormatter('%.d')
        cbar.ax.tick_params(labelsize=14)

        add_subplot_label(ax, j, i)

        if j == 0:
            ax.set_title(title, fontsize=16, pad=10)

        if i == 0:
            ax.annotate(ylabel_map[variable], xy=(-0.23, 0.5), xycoords='axes fraction',
                        fontsize=22, rotation=90, va='center')

plt.subplots_adjust(hspace=-0.55, wspace=0.08)
plt.savefig(os.path.join(OUT_DIR, 'Variations in Evapotranspiration in experimental scenarios.tif'),
            dpi=300, bbox_inches='tight')
plt.close()
print('  Saved.')

# ============================================================
# Figure 2: Multisource Scenarios (4 rows x 4 cols, all % diff)
# ============================================================
print('Generating multisource scenarios figure...')

combinations = [
    ('GEDI Mean', 'CLM Default'),
    ('GEDI Max', 'CLM Default'),
    ('GEDI Max', 'GEDI Mean'),
    ('GEDI Median', 'GEDI Mean')
]

fig2, axes2 = plt.subplots(4, 4, figsize=(24, 20),
                           subplot_kw=dict(projection=projection), dpi=200)

for j, variable in enumerate(variables):
    for i, (s1, s2) in enumerate(combinations):
        ax = axes2[j, i]

        data1 = get_et_data(s1) if variable == 'ET' else annual_data_EX[variable][s1]
        base = get_et_data(s2) if variable == 'ET' else annual_data_EX[variable][s2].copy()
        base = base.copy()
        base[base == 0] = np.nan
        plot_data = np.ma.masked_invalid((data1 - base) / base * 100)
        levels = np.linspace(-10, 10, num=21)

        cs = ax.contourf(lon2d, lat2d, plot_data, levels=levels, cmap=white_center_cmap,
                         extend='both', transform=projection)
        ax.coastlines()
        add_gridlines(ax, show_xlabel=(j == 0), show_ylabel=(i == 0))
        add_subplot_label(ax, j, i)

        if j == 0:
            ax.set_title(f'{s1}\nvs {s2} (%)', fontsize=16, pad=10)

        if i == 0:
            ax.annotate(ylabel_map[variable], xy=(-0.2, 0.5), xycoords='axes fraction',
                        fontsize=22, rotation=90, va='center')

cbar_ax = fig2.add_axes([0.30, 0.24, 0.42, 0.015])
cbar = fig2.colorbar(cs, cax=cbar_ax, orientation='horizontal', ticks=levels)
cbar.ax.tick_params(labelsize=14)
plt.subplots_adjust(hspace=-0.70, wspace=0.08)
plt.savefig(os.path.join(OUT_DIR, 'Variations of Evapotranspiration under Multisource Scenarios.tif'),
            dpi=300, bbox_inches='tight')
plt.close()
print('  Saved.')

# ============================================================
# Figure 3: Absolute difference (W/m2), 4 rows x 4 cols
# Col 0: CLM Default (mm/month), Cols 1-3: abs diff (W/m2)
# ============================================================
print('Generating absolute difference figure...')

fig3, axes3 = plt.subplots(4, 4, figsize=(24, 20),
                           subplot_kw=dict(projection=projection), dpi=200)

for j, variable in enumerate(variables):
    for i, scenario in enumerate(scenario_order):
        ax = axes3[j, i]

        data_et = get_et_data(scenario) if variable == 'ET' else annual_data_EX[variable][scenario]

        if scenario == 'CLM Default':
            base_data = data_et.copy()
            base_data[data_et == 0] = np.nan
            plot_data = np.ma.masked_invalid(base_data * 1.0562)
            levels = np.linspace(np.nanmin(plot_data), np.nanmax(plot_data), num=21)
            mn, mx = np.nanmin(plot_data), np.nanmax(plot_data)
            mid = 0.5 * (mn + mx)
            cbar_ticks = [mn, 0.5*(mid+mn), mid, 0.5*(mid+mx), mx]
            colormap = plt.get_cmap('afmhot_r')
            title = 'CLM Default\n(mm/month)'
        else:
            base_data = get_et_data('CLM Default') if variable == 'ET' else annual_data_EX[variable]['CLM Default'].copy()
            base_data = base_data.copy()
            diff = data_et - base_data
            diff[base_data == 0] = np.nan
            plot_data = np.ma.masked_invalid(diff)
            vmax = 5 if variable in ['ET', 'FCTR'] else 2
            levels = np.linspace(-vmax, vmax, num=21)
            colormap = white_center_cmap
            cbar_ticks = [-vmax, -vmax/2, 0, vmax/2, vmax]
            title = f'{scenario}\nvs Default (W/m$^2$)'

        cs = ax.contourf(lon2d, lat2d, plot_data, levels=levels, cmap=colormap,
                         extend='both', transform=projection)
        ax.coastlines()
        add_gridlines(ax, show_xlabel=(j == 0), show_ylabel=(i == 0))

        cbar = plt.colorbar(cs, ax=ax, orientation='horizontal', pad=0.005, shrink=0.6, ticks=cbar_ticks)
        cbar.formatter = ticker.FormatStrFormatter('%.d')
        cbar.ax.tick_params(labelsize=14)

        add_subplot_label(ax, j, i)

        if j == 0:
            ax.set_title(title, fontsize=16, pad=10)

        if i == 0:
            ax.annotate(ylabel_map[variable], xy=(-0.23, 0.5), xycoords='axes fraction',
                        fontsize=22, rotation=90, va='center')

plt.subplots_adjust(hspace=-0.55, wspace=0.08)
plt.savefig(os.path.join(OUT_DIR, 'Figure_absolute_difference.tiff'),
            dpi=300, bbox_inches='tight')
plt.close()
print('  Saved.')

# ============================================================
# Figure 4: Absolute values for all scenarios (mm/month)
# ============================================================
print('Generating absolute values figure...')

fig4, axes4 = plt.subplots(4, 4, figsize=(24, 20),
                           subplot_kw=dict(projection=projection), dpi=200)

for j, variable in enumerate(variables):
    for i, scenario in enumerate(scenario_order):
        ax = axes4[j, i]

        data_et = get_et_data(scenario) if variable == 'ET' else annual_data_EX[variable][scenario]
        base_data = data_et.copy()
        base_data[data_et == 0] = np.nan
        plot_data = np.ma.masked_invalid(base_data * 1.0562)
        levels = np.linspace(np.nanmin(plot_data), np.nanmax(plot_data), num=21)
        mn, mx = np.nanmin(plot_data), np.nanmax(plot_data)
        mid = 0.5 * (mn + mx)
        cbar_ticks = [mn, 0.5*(mid+mn), mid, 0.5*(mid+mx), mx]
        colormap = plt.get_cmap('afmhot_r')

        if scenario == 'CLM Default':
            title = 'CLM Default\n(mm/month)'
        else:
            title = f'{scenario}\n(mm/month)'

        cs = ax.contourf(lon2d, lat2d, plot_data, levels=levels, cmap=colormap,
                         extend='both', transform=projection)
        ax.coastlines()
        add_gridlines(ax, show_xlabel=(j == 0), show_ylabel=(i == 0))

        cbar = plt.colorbar(cs, ax=ax, orientation='horizontal', pad=0.005, shrink=0.6, ticks=cbar_ticks)
        cbar.formatter = ticker.FormatStrFormatter('%.d')
        cbar.ax.tick_params(labelsize=14)

        add_subplot_label(ax, j, i)

        if j == 0:
            ax.set_title(title, fontsize=16, pad=10)

        if i == 0:
            ax.annotate(ylabel_map[variable], xy=(-0.23, 0.5), xycoords='axes fraction',
                        fontsize=22, rotation=90, va='center')

plt.subplots_adjust(hspace=-0.55, wspace=0.08)
plt.savefig(os.path.join(OUT_DIR, 'Figure_absolute_values.tiff'),
            dpi=300, bbox_inches='tight')
plt.close()
print('  Saved.')

# ============================================================
# Figure 4 (ET comparison maps): 4 rows x 4 cols with labels
# Same as multisource but with the old Figure 4 styling
# ============================================================
print('Generating Figure 4 ET comparison maps...')

fig5, axes5 = plt.subplots(4, 4, figsize=(24, 20),
                           subplot_kw=dict(projection=projection), dpi=300)

for j, variable in enumerate(variables):
    for i, (s1, s2) in enumerate(combinations):
        ax = axes5[j, i]

        data1 = get_et_data(s1) if variable == 'ET' else annual_data_EX[variable][s1]
        base = get_et_data(s2) if variable == 'ET' else annual_data_EX[variable][s2].copy()
        base = base.copy()
        base[base == 0] = np.nan
        plot_data = np.ma.masked_invalid((data1 - base) / base * 100)
        levels = np.linspace(-10, 10, num=11)

        cs = ax.contourf(lon2d, lat2d, plot_data, levels=levels, cmap=white_center_cmap,
                         extend='both', transform=projection)
        ax.coastlines()
        add_gridlines(ax, show_xlabel=(j == 0), show_ylabel=(i == 0))
        add_subplot_label(ax, j, i)

        if j == 0:
            ax.set_title(f'{s1}\nvs {s2}', fontsize=16, pad=10)

        if i == 0:
            ax.annotate(ylabel_map[variable], xy=(-0.2, 0.5), xycoords='axes fraction',
                        fontsize=22, rotation=90, va='center')

cbar_ax = fig5.add_axes([0.30, 0.24, 0.42, 0.015])
cbar = fig5.colorbar(cs, cax=cbar_ax, orientation='horizontal', ticks=levels)
cbar.set_label("Relative Difference (%)", fontsize=16)
cbar.ax.tick_params(labelsize=14)
plt.subplots_adjust(hspace=-0.70, wspace=0.08)
plt.savefig(os.path.join(OUT_DIR, 'Figure_4_ET_comparison_maps.tiff'),
            dpi=300, bbox_inches='tight')
plt.close()
print('  Saved.')

print('\nAll figures regenerated!')
