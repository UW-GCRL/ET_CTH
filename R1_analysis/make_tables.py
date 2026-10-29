"""Table 1 and Tables S1-S4 of the revised manuscript as CSV files (tables/).
Run after an3_spatial.py (reads an3_results.json). Means are forest-area-weighted (patch land area,
wtgcell x cell area x land fraction); relative differences are ratios of these means, as in Table 1."""
import os, csv, json
import numpy as np
from common import *

OUTT = 'tables/'
os.makedirs(OUTT, exist_ok=True)
m = load_meta(); f = m['forest']
w = m['wt'] * m['area'][m['jx'], m['ix']]
D = {s: load_scen(s) for s in SCEN}
a3 = json.load(open('an3_results.json'))
F4 = ['FCEV', 'FCTR', 'FGEV', 'ET']


def wm(x, sel):
    return float(np.nansum((x * w)[sel]) / np.nansum(w[sel]))


def z(s):                                   # '-0.00' -> '0.00'
    return s[1:] if s.startswith('-') and float(s) == 0 else s


def f1(x): return z(f'{x:.1f}')
def f2(x): return z(f'{x:.2f}')
def f3(x): return z(f'{x:.3f}')


def write(name, rows):
    with open(OUTT + name, 'w', encoding='utf-8', newline='') as fh:
        csv.writer(fh).writerows(rows)
    print(name, len(rows) - 1, 'rows')


# Table 1: global relative differences with 95% block-bootstrap intervals; dagger = not significant in either test
rows = [['Comparison', 'ET', 'FCEV', 'FCTR', 'FGEV']]
for c in ['Default 0.8 vs CLM Default', 'Default 0.9 vs CLM Default', 'Default 1.1 vs CLM Default', 'Default 1.2 vs CLM Default',
          'GEDI Mean vs CLM Default', 'GEDI Max vs CLM Default', 'GEDI Median vs CLM Default', 'GEDI Max vs GEDI Mean', 'GEDI Median vs GEDI Mean']:
    r = [c]
    for v in ['ET', 'FCEV', 'FCTR', 'FGEV']:
        x = a3['global'][f'{c}|{v}']
        ns = ' †' if max(x['q_neff_t'], x['q_block_wilcoxon']) >= 0.05 else ''
        r.append(f"{x['rel_pct']:+.2f}{ns} [{x['ci_lo']:+.2f}, {x['ci_hi']:+.2f}]".replace('-', '−'))
    rows.append(r)
write('Table1_global_differences.csv', rows)

# Table S1: global forest fluxes per unit forest area (W m-2)
rows = [['Scenario', 'FCEV (W m-2)', 'FCTR (W m-2)', 'FGEV (W m-2)', 'ET (W m-2)', 'FSH (W m-2)']]
for s in SCEN:
    rows.append([s] + [f2(wm(D[s][v], f)) for v in ['FCEV', 'FCTR', 'FGEV', 'ET', 'FSH']])
write('TableS1_global_forest_fluxes.csv', rows)

# Table S2: PFT-level base fluxes, relative and absolute differences
rows = [['PFT', 'Comparison', 'CLM Default FCEV (W m-2)', 'CLM Default FCTR (W m-2)', 'CLM Default FGEV (W m-2)', 'CLM Default ET (W m-2)',
         'FCEV relative difference (%)', 'FCTR relative difference (%)', 'FGEV relative difference (%)', 'ET relative difference (%)',
         'FCEV absolute difference (W m-2)', 'FCTR absolute difference (W m-2)', 'FGEV absolute difference (W m-2)', 'ET absolute difference (W m-2)']]
for p, k in PFTS.items():
    sel = f & (m['iv'] == k)
    base = [wm(D['CLM Default'][v], sel) for v in F4]
    for a, b in [('Default 0.8', 'CLM Default'), ('Default 1.2', 'CLM Default'), ('GEDI Mean', 'CLM Default'), ('GEDI Max', 'CLM Default'),
                 ('GEDI Median', 'CLM Default'), ('GEDI Max', 'GEDI Mean')]:
        ma = [wm(D[a][v], sel) for v in F4]; mb = [wm(D[b][v], sel) for v in F4]
        rows.append([p, f'{a} vs {b}'] + [f2(x) for x in base] + [f2(100 * (x / y - 1)) for x, y in zip(ma, mb)] +
                    [f3(x - y) for x, y in zip(ma, mb)])
write('TableS2_PFT_fluxes_and_differences.csv', rows)

# Table S3: PFT-level significance (from an3_spatial.py)
rows = [['PFT', 'Comparison', 'Flux', 'Relative difference (%)', '95% interval low (%)', '95% interval high (%)', 'Effective sample size',
         'q, t-test with effective sample size', 'q, block-level Wilcoxon test', 'q, cell-level Wilcoxon test']]
for p in PFTS:
    for key, v in a3[p].items():
        c, fl = key.split('|')
        rows.append([p, c, fl, f2(v['rel_pct']), f2(v['ci_lo']), f2(v['ci_hi']), str(int(round(v['n_eff']))),
                     f"{v['q_neff_t']:.2g}", f"{v['q_block_wilcoxon']:.2g}", f"{v['q_cell_wilcoxon']:.2g}"])
write('TableS3_PFT_significance.csv', rows)

# Table S4: canopy height inputs where GEDI-based heights were assigned, and share of PFT area that kept the default
H = {s: patch_height(s, m) for s in ['CLM Default', 'GEDI Max', 'GEDI Mean', 'GEDI Median']}
ret = f & (np.abs(H['GEDI Mean'] - H['CLM Default']) < 1e-4)
rows = [['PFT', 'CLM Default (m)', 'GEDI Max (m)', 'GEDI Mean (m)', 'GEDI Median (m)', 'Area with default height retained (%)']]
for p, sel in [(p, f & (m['iv'] == k)) for p, k in PFTS.items()] + [('All forest', f)]:
    rows.append([p] + [f1(wm(H[s], sel & ~ret)) for s in H] + [f1(100 * w[sel & ret].sum() / w[sel].sum())])
write('TableS4_canopy_height_inputs.csv', rows)
