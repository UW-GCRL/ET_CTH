"""
Generate PFT heatmaps as 2x2 PNG, no suptitle, with significance test results.
Uses individual h1 yearly files (2031-2040) for paired t-tests.
"""
import netCDF4
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import sys
sys.stdout.reconfigure(encoding='utf-8')

OUT_DIR = 'C:/Users/hyou34/Documents/CTH_ET/figures_new_version/'
FINAL_DATA = 'G:/Hangkai/CTH_ET project/Final_data/'
H1_DIR = 'C:/Users/hyou34/Downloads/temp_h1/'

pft_names = {
    "NET Temp": 1, "NET Bor": 2, "NDT Bor": 3,
    "BET Trop": 4, "BET Temp": 5,
    "BDT Trop": 6, "BDT Temp": 7, "BDT Bor": 8,
}
pft_list = list(pft_names.keys())
pft_vars = ['FCEV', 'FCTR', 'FGEV', 'ET']
var_labels = {'FCEV': 'Canopy Evaporation', 'FCTR': 'Canopy Transpiration',
              'FGEV': 'Ground Evaporation', 'ET': 'Total ET'}

# Case name mapping
case_names = {
    'CLM Default': 'Default_case',
    'Default 0.8': 'Default_0.8_case',
    'Default 0.9': 'Default_0.9_case',
    'Default 1.1': 'Default_1.1_case',
    'Default 1.2': 'Default_1.2_case',
    'GEDI Mean': 'Mean_case',
    'GEDI Max': 'Max_case',
    'GEDI Median': 'Median_case',
}

scenarios_order = ['CLM Default', 'Default 0.8', 'Default 0.9', 'Default 1.1', 'Default 1.2',
                   'GEDI Mean', 'GEDI Max', 'GEDI Median']

# ============================================================
# Load yearly h1 data for all cases (2031-2040)
# ============================================================
print('Loading yearly h1 data...')
yearly_data = {}  # {case_dir: {year: {var: pft_array}}}

for scenario, case_dir in case_names.items():
    yearly_data[scenario] = {}
    for yr in range(2031, 2041):
        fname = f'{case_dir}.clm2.h1.{yr}-01-01-00000.nc'
        fpath = os.path.join(H1_DIR, case_dir, fname)
        if os.path.exists(fpath):
            with netCDF4.Dataset(fpath, 'r') as ds:
                yearly_data[scenario][yr] = {}
                for var in ['FCEV', 'FCTR', 'FGEV']:
                    raw = ds.variables[var][:]
                    yearly_data[scenario][yr][var] = raw.squeeze()
                yearly_data[scenario][yr]['ET'] = (yearly_data[scenario][yr]['FCEV'] +
                                                   yearly_data[scenario][yr]['FCTR'] +
                                                   yearly_data[scenario][yr]['FGEV'])
                if 'pfts1d_itype_veg' not in yearly_data[scenario]:
                    yearly_data[scenario]['pfts1d_itype_veg'] = ds.variables['pfts1d_itype_veg'][:]
                    yearly_data[scenario]['pfts1d_wtgcell'] = ds.variables['pfts1d_wtgcell'][:]
                    yearly_data[scenario]['area'] = ds.variables['area'][:]
                    yearly_data[scenario]['pfts1d_jxy'] = ds.variables['pfts1d_jxy'][:]
                    yearly_data[scenario]['pfts1d_ixy'] = ds.variables['pfts1d_ixy'][:]

print('Yearly data loaded.')

# ============================================================
# Compute PFT-scale yearly area-weighted means
# ============================================================
area = yearly_data['CLM Default']['area']
grid_shape = area.shape
area_flat = area.flatten()

# {scenario: {pft: {var: [10 yearly values]}}}
pft_yearly = {}
for scenario in scenarios_order:
    meta = yearly_data[scenario]
    pft_yearly[scenario] = {}
    for pft_name, pft_id in pft_names.items():
        pft_idx = np.where(meta['pfts1d_itype_veg'] == pft_id)[0]
        if len(pft_idx) == 0:
            pft_yearly[scenario][pft_name] = {v: np.zeros(10) for v in pft_vars}
            continue

        aw = np.zeros(len(pft_idx))
        for k, idx in enumerate(pft_idx):
            row = int(meta['pfts1d_jxy'][idx]) - 1
            col = int(meta['pfts1d_ixy'][idx]) - 1
            fi = row * grid_shape[1] + col
            aw[k] = area_flat[fi] if fi < len(area_flat) else 1.0

        pft_yearly[scenario][pft_name] = {}
        for var in pft_vars:
            yr_vals = []
            for yr in range(2031, 2041):
                vals = meta[yr][var][pft_idx]
                weighted = vals * meta['pfts1d_wtgcell'][pft_idx]
                wmean = np.average(weighted, weights=aw)
                yr_vals.append(wmean)
            pft_yearly[scenario][pft_name][var] = np.array(yr_vals)

print('PFT yearly means computed.')

# ============================================================
# Compute 10yr means and paired t-tests
# ============================================================
pft_combos = [
    ('GEDI Mean', 'CLM Default'),
    ('GEDI Max', 'CLM Default'),
    ('GEDI Median', 'CLM Default'),
    ('Default 0.8', 'CLM Default'),
    ('Default 0.9', 'CLM Default'),
    ('Default 1.1', 'CLM Default'),
    ('Default 1.2', 'CLM Default'),
]
combo_labels = ['Mean\nvs Def', 'Max\nvs Def', 'Med\nvs Def',
                '0.8x\nvs Def', '0.9x\nvs Def', '1.1x\nvs Def', '1.2x\nvs Def']

# Compute absolute means
pft_abs = {}
for scenario in scenarios_order:
    pft_abs[scenario] = {}
    for pft_name in pft_list:
        pft_abs[scenario][pft_name] = {}
        for var in pft_vars:
            pft_abs[scenario][pft_name][var] = np.mean(pft_yearly[scenario][pft_name][var])

# Significance test results
# {(s1,s2): {pft: {var: p_value}}}
sig_results = {}
for s1, s2 in pft_combos:
    sig_results[(s1, s2)] = {}
    for pft_name in pft_list:
        sig_results[(s1, s2)][pft_name] = {}
        for var in pft_vars:
            arr1 = pft_yearly[s1][pft_name][var]
            arr2 = pft_yearly[s2][pft_name][var]
            if np.std(arr1 - arr2) < 1e-15:
                p_val = 1.0
            else:
                t_stat, p_val = stats.ttest_rel(arr1, arr2)
            sig_results[(s1, s2)][pft_name][var] = p_val

print('Significance tests computed.')

# Print summary
for combo in pft_combos[:3]:
    s1, s2 = combo
    print(f'\n{s1} vs {s2}:')
    for pft_name in pft_list:
        for var in pft_vars:
            p = sig_results[combo][pft_name][var]
            sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
            pct = 100 * (pft_abs[s1][pft_name][var] - pft_abs[s2][pft_name][var]) / pft_abs[s2][pft_name][var] if abs(pft_abs[s2][pft_name][var]) > 1e-10 else 0
            if sig != 'ns':
                print(f'  {pft_name:>12} {var:>5}: {pct:+.2f}% p={p:.4f} {sig}')

# ============================================================
# Heatmap helper
# ============================================================
def make_heatmap_2x2(matrix_dict, sig_dict, pft_list, combo_labels, pft_combos, var_labels,
                     value_fmt, cmap, vmin_vmax_func, cbar_label, filename):
    """Generate a 2x2 heatmap figure."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), dpi=300)

    for idx, var in enumerate(pft_vars):
        ax = axes[idx // 2, idx % 2]
        matrix = np.zeros((len(pft_list), len(pft_combos)))
        sigs = np.empty((len(pft_list), len(pft_combos)), dtype=object)

        for j, (s1, s2) in enumerate(pft_combos):
            for i, pft in enumerate(pft_list):
                matrix[i, j] = matrix_dict[var][(s1, s2)][pft]
                p = sig_dict[(s1, s2)][pft][var]
                if p < 0.001:
                    sigs[i, j] = '***'
                elif p < 0.01:
                    sigs[i, j] = '**'
                elif p < 0.05:
                    sigs[i, j] = '*'
                else:
                    sigs[i, j] = ''

        vmin, vmax = vmin_vmax_func(matrix)
        im = ax.imshow(matrix, cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
        cbar = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
        cbar.set_label(cbar_label, fontsize=10)
        cbar.ax.tick_params(labelsize=9)

        # Annotate values + significance
        for i in range(len(pft_list)):
            for j in range(len(pft_combos)):
                val = matrix[i, j]
                sig_marker = sigs[i, j]
                color = 'white' if abs(val - (vmin+vmax)/2) > (vmax-vmin)*0.3 else 'black'
                text = value_fmt(val)
                if sig_marker:
                    text += f'\n{sig_marker}'
                ax.text(j, i, text, ha='center', va='center', fontsize=7, color=color, fontweight='bold' if sig_marker else 'normal')

        ax.set_xticks(range(len(pft_combos)))
        ax.set_xticklabels(combo_labels, fontsize=8, ha='center')
        ax.set_yticks(range(len(pft_list)))
        ax.set_yticklabels(pft_list, fontsize=10)
        panel_label = chr(ord('a') + idx)
        ax.set_title(f'({panel_label}) {var_labels[var]}', fontsize=13, fontweight='bold', loc='left')
        ax.axvline(x=2.5, color='black', linewidth=1.5, linestyle='--')

    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, filename), dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  Saved: {filename}')

# ============================================================
# Prepare data matrices
# ============================================================
# Percentage differences
pct_matrices = {var: {} for var in pft_vars}
abs_diff_matrices = {var: {} for var in pft_vars}

for s1, s2 in pft_combos:
    for var in pft_vars:
        pct_matrices[var][(s1, s2)] = {}
        abs_diff_matrices[var][(s1, s2)] = {}
        for pft in pft_list:
            base_v = pft_abs[s2][pft][var]
            new_v = pft_abs[s1][pft][var]
            pct = 100 * (new_v - base_v) / base_v if abs(base_v) > 1e-10 else 0
            pct_matrices[var][(s1, s2)][pft] = pct
            abs_diff_matrices[var][(s1, s2)][pft] = new_v - base_v

# ============================================================
# Heatmap 1: PFT percentage differences (2x2 PNG)
# ============================================================
print('Generating PFT percentage difference heatmap...')
make_heatmap_2x2(
    pct_matrices, sig_results, pft_list, combo_labels, pft_combos, var_labels,
    value_fmt=lambda v: f'{v:+.1f}',
    cmap='RdBu_r',
    vmin_vmax_func=lambda m: (-max(abs(m.min()), abs(m.max()), 1), max(abs(m.min()), abs(m.max()), 1)),
    cbar_label='(%)',
    filename='Heatmap_PFT_percentage_differences.png'
)

# ============================================================
# Heatmap 2: PFT absolute differences (2x2 PNG)
# ============================================================
print('Generating PFT absolute difference heatmap...')
make_heatmap_2x2(
    abs_diff_matrices, sig_results, pft_list, combo_labels, pft_combos, var_labels,
    value_fmt=lambda v: f'{v:+.2f}' if abs(v) >= 0.01 else f'{v:+.3f}',
    cmap='RdBu_r',
    vmin_vmax_func=lambda m: (-max(abs(m.min()), abs(m.max()), 0.01), max(abs(m.min()), abs(m.max()), 0.01)),
    cbar_label='(W/m$^2$)',
    filename='Heatmap_PFT_absolute_differences.png'
)

# ============================================================
# Heatmap 3: PFT absolute values (2x2 PNG)
# ============================================================
print('Generating PFT absolute value heatmap...')

scen_labels = ['Default', '0.8x', '0.9x', '1.1x', '1.2x', 'GEDI\nMean', 'GEDI\nMax', 'GEDI\nMedian']

fig, axes = plt.subplots(2, 2, figsize=(14, 10), dpi=300)
for idx, var in enumerate(pft_vars):
    ax = axes[idx // 2, idx % 2]
    matrix = np.zeros((len(pft_list), len(scenarios_order)))
    for j, scenario in enumerate(scenarios_order):
        for i, pft in enumerate(pft_list):
            matrix[i, j] = pft_abs[scenario][pft][var]

    im = ax.imshow(matrix, cmap='YlOrRd', aspect='auto')
    cbar = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
    cbar.set_label('(W/m$^2$)', fontsize=10)
    cbar.ax.tick_params(labelsize=9)

    for i in range(len(pft_list)):
        for j in range(len(scenarios_order)):
            val = matrix[i, j]
            color = 'white' if val > (matrix.max() + matrix.min()) / 2 else 'black'
            ax.text(j, i, f'{val:.2f}', ha='center', va='center', fontsize=7, color=color)

    ax.set_xticks(range(len(scenarios_order)))
    ax.set_xticklabels(scen_labels, fontsize=9, rotation=0, ha='center')
    ax.set_yticks(range(len(pft_list)))
    ax.set_yticklabels(pft_list, fontsize=10)
    panel_label = chr(ord('a') + idx)
    ax.set_title(f'({panel_label}) {var_labels[var]}', fontsize=13, fontweight='bold', loc='left')

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, 'Heatmap_PFT_absolute_values.png'), dpi=300, bbox_inches='tight')
plt.close()
print('  Saved: Heatmap_PFT_absolute_values.png')

# Clean up old tiff versions
for old_tiff in ['Heatmap_PFT_percentage_differences.tiff',
                  'Heatmap_PFT_absolute_differences.tiff',
                  'Heatmap_PFT_absolute_values.tiff']:
    old_path = os.path.join(OUT_DIR, old_tiff)
    if os.path.exists(old_path):
        os.remove(old_path)
        print(f'  Removed old: {old_tiff}')

print('\nAll done!')
