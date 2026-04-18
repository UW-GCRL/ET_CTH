"""
Generate global-scale and PFT-scale summary tables for CTH_ET manuscript.
Outputs CSV files for easy inclusion in manuscript/supplementary.
"""
import netCDF4
import os
import numpy as np
import csv
import sys
sys.stdout.reconfigure(encoding='utf-8')

OUT_DIR = 'C:/Users/hyou34/Documents/CTH_ET/figures_new_version/'
FINAL_DATA = 'G:/Hangkai/CTH_ET project/Final_data/'

# ============================================================
# Load data
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

pft_names = {
    "NET Temperate": 1, "NET Boreal": 2, "NDT Boreal": 3,
    "BET Tropical": 4, "BET Temperate": 5,
    "BDT Tropical": 6, "BDT Temperate": 7, "BDT Boreal": 8,
}

data_all = {}
pft_meta = {}

for scenario, fname in scenario_mapping.items():
    data_all[scenario] = {}
    nc_file = os.path.join(FINAL_DATA, fname)
    with netCDF4.Dataset(nc_file, 'r') as ds:
        for var in variable_names:
            raw = ds.variables[var][:]
            data_all[scenario][var] = np.mean(raw, axis=0) if raw.ndim == 2 else raw
        data_all[scenario]['ET'] = data_all[scenario]['FCEV'] + data_all[scenario]['FCTR'] + data_all[scenario]['FGEV']
        # Also load FSH
        raw_fsh = ds.variables['FSH'][:]
        data_all[scenario]['FSH'] = np.mean(raw_fsh, axis=0) if raw_fsh.ndim == 2 else raw_fsh

        pft_meta[scenario] = {
            'pfts1d_itype_veg': ds.variables['pfts1d_itype_veg'][:],
            'pfts1d_ixy': ds.variables['pfts1d_ixy'][:],
            'pfts1d_jxy': ds.variables['pfts1d_jxy'][:],
            'pfts1d_wtgcell': ds.variables['pfts1d_wtgcell'][:],
            'area': ds.variables['area'][:],
        }

area = pft_meta['CLM Default']['area']
grid_shape = area.shape

# Grid to 2D
all_vars = ['FCEV', 'FCTR', 'FGEV', 'ET', 'FSH']
annual_data_2d = {var: {} for var in all_vars}

for scenario in scenario_mapping:
    meta = pft_meta[scenario]
    for variable in all_vars:
        grid_2d = np.zeros(grid_shape)
        values = data_all[scenario][variable] * meta['pfts1d_wtgcell']
        for pft_name, pft_id in pft_names.items():
            pft_idx = np.where(meta['pfts1d_itype_veg'] == pft_id)[0]
            for i, idx in enumerate(pft_idx):
                row = int(meta['pfts1d_jxy'][idx]) - 1
                col = int(meta['pfts1d_ixy'][idx]) - 1
                grid_2d[row, col] += values[idx]
        annual_data_2d[variable][scenario] = grid_2d

print('Data loaded and gridded.')

# ============================================================
# Table 1: Global-scale absolute values (W/m2)
# ============================================================
scenarios_order = ['CLM Default', 'Default 0.8', 'Default 0.9', 'Default 1.1', 'Default 1.2',
                   'GEDI Mean', 'GEDI Max', 'GEDI Median']
table_vars = ['FCEV', 'FCTR', 'FGEV', 'ET', 'FSH']

print('\n' + '=' * 80)
print('TABLE 1: Global-Scale Absolute Values (W/m2)')
print('=' * 80)
header = f'{"Scenario":<18}' + ''.join(f'{v:>10}' for v in table_vars)
print(header)
print('-' * 68)

rows_t1 = []
global_abs = {}
for scenario in scenarios_order:
    vals = {}
    for var in table_vars:
        d = annual_data_2d[var][scenario]
        vals[var] = np.nansum(d * area) / np.nansum(area)
    global_abs[scenario] = vals
    row_str = f'{scenario:<18}' + ''.join(f'{vals[v]:>10.4f}' for v in table_vars)
    print(row_str)
    rows_t1.append([scenario] + [f'{vals[v]:.4f}' for v in table_vars])

# Save CSV
with open(os.path.join(OUT_DIR, 'Table_1_global_absolute_values.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['Scenario'] + [f'{v} (W/m2)' for v in table_vars])
    w.writerows(rows_t1)

# ============================================================
# Table 2: Global-scale percentage differences
# ============================================================
combinations = [
    ('GEDI Mean', 'CLM Default'),
    ('GEDI Max', 'CLM Default'),
    ('GEDI Median', 'CLM Default'),
    ('GEDI Max', 'GEDI Mean'),
    ('GEDI Median', 'GEDI Mean'),
    ('Default 0.8', 'CLM Default'),
    ('Default 0.9', 'CLM Default'),
    ('Default 1.1', 'CLM Default'),
    ('Default 1.2', 'CLM Default'),
]

print('\n' + '=' * 80)
print('TABLE 2: Global-Scale Percentage Differences (%)')
print('=' * 80)
header2 = f'{"Comparison":<35}' + ''.join(f'{v:>10}' for v in table_vars)
print(header2)
print('-' * 85)

rows_t2 = []
for s1, s2 in combinations:
    vals = {}
    for var in table_vars:
        pct = 100 * (global_abs[s1][var] - global_abs[s2][var]) / global_abs[s2][var]
        vals[var] = pct
    row_str = f'{s1 + " vs " + s2:<35}' + ''.join(f'{vals[v]:>+10.4f}' for v in table_vars)
    print(row_str)
    rows_t2.append([f'{s1} vs {s2}'] + [f'{vals[v]:+.4f}' for v in table_vars])

with open(os.path.join(OUT_DIR, 'Table_2_global_percentage_differences.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['Comparison'] + [f'{v} (%)' for v in table_vars])
    w.writerows(rows_t2)

# ============================================================
# Table 3: Global-scale absolute differences (W/m2)
# ============================================================
print('\n' + '=' * 80)
print('TABLE 3: Global-Scale Absolute Differences (W/m2)')
print('=' * 80)
print(header2)
print('-' * 85)

rows_t3 = []
for s1, s2 in combinations:
    vals = {}
    for var in table_vars:
        diff = global_abs[s1][var] - global_abs[s2][var]
        vals[var] = diff
    row_str = f'{s1 + " vs " + s2:<35}' + ''.join(f'{vals[v]:>+10.4f}' for v in table_vars)
    print(row_str)
    rows_t3.append([f'{s1} vs {s2}'] + [f'{vals[v]:+.4f}' for v in table_vars])

with open(os.path.join(OUT_DIR, 'Table_3_global_absolute_differences.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['Comparison'] + [f'{v} (W/m2)' for v in table_vars])
    w.writerows(rows_t3)

# ============================================================
# Table 4: PFT-scale absolute values (W/m2)
# Area-weighted mean per PFT for each scenario
# ============================================================
print('\n' + '=' * 80)
print('TABLE 4: PFT-Scale Absolute Values (W/m2) - Area-Weighted Mean')
print('=' * 80)

gedi_scenarios = ['CLM Default', 'GEDI Mean', 'GEDI Max', 'GEDI Median']
pft_vars = ['FCEV', 'FCTR', 'FGEV', 'ET']
area_flat = area.flatten()

# Compute PFT-scale weighted means
pft_abs = {}  # {scenario: {pft: {var: value}}}
for scenario in scenarios_order:
    meta = pft_meta[scenario]
    pft_abs[scenario] = {}
    for pft_name, pft_id in pft_names.items():
        pft_idx = np.where(meta['pfts1d_itype_veg'] == pft_id)[0]
        if len(pft_idx) == 0:
            pft_abs[scenario][pft_name] = {v: 0 for v in pft_vars}
            continue

        # Area weights for this PFT
        aw = np.zeros(len(pft_idx))
        for k, idx in enumerate(pft_idx):
            row = int(meta['pfts1d_jxy'][idx]) - 1
            col = int(meta['pfts1d_ixy'][idx]) - 1
            fi = row * grid_shape[1] + col
            aw[k] = area_flat[fi] if fi < len(area_flat) else 1.0

        pft_abs[scenario][pft_name] = {}
        for var in pft_vars:
            vals = data_all[scenario][var][pft_idx]
            # Weight by pft area fraction and grid area
            weighted = vals * meta['pfts1d_wtgcell'][pft_idx]
            wmean = np.average(weighted, weights=aw)
            pft_abs[scenario][pft_name][var] = wmean

# Print PFT absolute values for GEDI scenarios
rows_t4 = []
for var in pft_vars:
    print(f'\n--- {var} (W/m2) ---')
    header4 = f'{"PFT":<18}' + ''.join(f'{s:>16}' for s in scenarios_order)
    print(header4)
    print('-' * (18 + 16 * len(scenarios_order)))
    for pft_name in pft_names:
        row_str = f'{pft_name:<18}'
        row_data = [var, pft_name]
        for scenario in scenarios_order:
            v = pft_abs[scenario][pft_name][var]
            row_str += f'{v:>16.4f}'
            row_data.append(f'{v:.4f}')
        print(row_str)
        rows_t4.append(row_data)

with open(os.path.join(OUT_DIR, 'Table_4_pft_absolute_values.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['Variable', 'PFT'] + [f'{s} (W/m2)' for s in scenarios_order])
    w.writerows(rows_t4)

# ============================================================
# Table 5: PFT-scale percentage differences
# ============================================================
print('\n' + '=' * 80)
print('TABLE 5: PFT-Scale Percentage Differences (%)')
print('=' * 80)

pft_combos = [
    ('GEDI Mean', 'CLM Default'),
    ('GEDI Max', 'CLM Default'),
    ('GEDI Median', 'CLM Default'),
    ('GEDI Max', 'GEDI Mean'),
    ('GEDI Median', 'GEDI Mean'),
]

rows_t5 = []
for s1, s2 in pft_combos:
    print(f'\n--- {s1} vs {s2} ---')
    header5 = f'{"PFT":<18}' + ''.join(f'{v:>10}' for v in pft_vars)
    print(header5)
    print('-' * (18 + 10 * len(pft_vars)))
    for pft_name in pft_names:
        row_str = f'{pft_name:<18}'
        row_data = [f'{s1} vs {s2}', pft_name]
        for var in pft_vars:
            base_v = pft_abs[s2][pft_name][var]
            new_v = pft_abs[s1][pft_name][var]
            if abs(base_v) > 1e-10:
                pct = 100 * (new_v - base_v) / base_v
            else:
                pct = 0
            row_str += f'{pct:>+10.2f}'
            row_data.append(f'{pct:+.2f}')
        print(row_str)
        rows_t5.append(row_data)

with open(os.path.join(OUT_DIR, 'Table_5_pft_percentage_differences.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['Comparison', 'PFT'] + [f'{v} (%)' for v in pft_vars])
    w.writerows(rows_t5)

# ============================================================
# Table 6: PFT-scale absolute differences (W/m2)
# ============================================================
print('\n' + '=' * 80)
print('TABLE 6: PFT-Scale Absolute Differences (W/m2)')
print('=' * 80)

rows_t6 = []
for s1, s2 in pft_combos:
    print(f'\n--- {s1} vs {s2} ---')
    header6 = f'{"PFT":<18}' + ''.join(f'{v:>10}' for v in pft_vars)
    print(header6)
    print('-' * (18 + 10 * len(pft_vars)))
    for pft_name in pft_names:
        row_str = f'{pft_name:<18}'
        row_data = [f'{s1} vs {s2}', pft_name]
        for var in pft_vars:
            diff = pft_abs[s1][pft_name][var] - pft_abs[s2][pft_name][var]
            row_str += f'{diff:>+10.4f}'
            row_data.append(f'{diff:+.4f}')
        print(row_str)
        rows_t6.append(row_data)

with open(os.path.join(OUT_DIR, 'Table_6_pft_absolute_differences.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['Comparison', 'PFT'] + [f'{v} (W/m2)' for v in pft_vars])
    w.writerows(rows_t6)

print('\n===== ALL TABLES SAVED =====')
print(f'CSV files saved to: {OUT_DIR}')
