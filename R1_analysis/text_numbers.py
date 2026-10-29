"""Every model-derived number quoted in the text of the revised manuscript (Sections 3 and 4, Abstract, Conclusions,
Text S3), recomputed from the model output. Run after an1-an4 (reads their *_results.json and emulator_coef.npz).
Output: text_numbers.json (one entry per quantity; keys name the comparison, PFT and flux)."""
import json
import numpy as np
from scipy import stats
from common import *

J = {k: json.load(open(k + '_results.json')) for k in ['an1', 'an2', 'an3', 'an4']}
m = load_meta(); f = m['forest']; iv = m['iv']
A = np.nan_to_num(m['area'])                      # land area of each cell (common.load_meta)
w = m['wt'] * A[m['jx'], m['ix']]
plat = m['pfts1d_lat']
S = {s: load_scen(s) for s in SCEN}
for s in S:
    S[s]['GPP'] = S[s]['FPSN']
Q = {}


def wm(x, sel):
    return float(np.nansum((x * w)[sel]) / np.nansum(w[sel]))


def rel(a, b, v, sel=f):
    return 100 * (wm(S[a][v], sel) / wm(S[b][v], sel) - 1)


def dabs(a, b, v, sel=f):
    return wm(S[a][v], sel) - wm(S[b][v], sel)


PF = {p: f & (iv == k) for p, k in PFTS.items()}
COMPS = [('Default 0.8', 'CLM Default'), ('Default 0.9', 'CLM Default'), ('Default 1.1', 'CLM Default'), ('Default 1.2', 'CLM Default'),
         ('GEDI Mean', 'CLM Default'), ('GEDI Max', 'CLM Default'), ('GEDI Median', 'CLM Default'), ('GEDI Max', 'GEDI Mean'), ('GEDI Median', 'GEDI Mean')]
FL = ['FCEV', 'FCTR', 'FGEV', 'ET', 'FSH', 'Rnet', 'GPP']
for a, b in COMPS:
    for v in FL:
        Q[f'rel|{a}|{b}|{v}'] = rel(a, b, v)
        Q[f'abs|{a}|{b}|{v}'] = dabs(a, b, v)
        for p, sel in PF.items():
            Q[f'relP|{p}|{a}|{b}|{v}'] = rel(a, b, v, sel)
            Q[f'absP|{p}|{a}|{b}|{v}'] = dabs(a, b, v, sel)
for s in SCEN:
    for v in FL:
        Q[f'base|{s}|{v}'] = wm(S[s][v], f)
# water-flux units (mm/yr) for selected differences
for a, b in [('GEDI Max', 'GEDI Mean'), ('Default 1.2', 'CLM Default')]:
    for v in ['T', 'ET_w', 'E_can', 'E_gnd']:
        Q[f'mm|{a}|{b}|{v}'] = dabs(a, b, v)
# compensation / offsets (absolute, W m-2)
def offs(a, b, sel=f):
    return -(dabs(a, b, 'FCEV', sel) + dabs(a, b, 'FGEV', sel)) / dabs(a, b, 'FCTR', sel)
Q['comp|Default 1.2|global'] = offs('Default 1.2', 'CLM Default')
for p in PF:
    Q[f'comp|GEDI Mean|{p}'] = offs('GEDI Mean', 'CLM Default', PF[p])
sel3 = PF['BET Trop'] | PF['BET Temp'] | PF['BDT Temp']
Q['comp|GEDI Mean|3 broadleaf'] = offs('GEDI Mean', 'CLM Default', sel3)
# regional canopy evaporation, Default 1.2
for name, sel in [('NH>=50', f & (plat >= 50)), ('|lat|>=50', f & (np.abs(plat) >= 50)), ('tropics', f & (np.abs(plat) < 23.5)),
                  ('temperate', f & (np.abs(plat) >= 23.5) & (np.abs(plat) < 50))]:
    Q[f'regFCEV|{name}'] = rel('Default 1.2', 'CLM Default', 'FCEV', sel)
# seasonal shares, NH >= 50 (an2 'Boreal')
for k, v in J['an2']['seasonal']['Default 1.2 vs CLM Default'].items():
    mo = np.array(v['Ecan_monthly_mm']); tot = mo.sum()
    Q[f'season|{k}|JunNov'] = 100 * mo[5:11].sum() / tot
    Q[f'season|{k}|DJF'] = 100 * (mo[11] + mo[0] + mo[1]) / tot
    Q[f'season|{k}|Ecan_rel'] = v['Ecan_rel_pct']
Q['runoff|Default 1.2 mm'] = J['an2']['implied_runoff']['Default 1.2 vs CLM Default']['dET_mm_per_forest_area']
# snow share of the reduction in water reaching the ground (column level, forest-weighted), global and NH >= 50
z = {s: np.load(R1 + SCEN[s][2] + '_annual.npz') for s in ['CLM Default', 'Default 1.2', 'Default 0.8', 'GEDI Max', 'GEDI Mean', 'GEDI Median']}
ci = m['ci']
for a, b in [('Default 1.2', 'CLM Default'), ('Default 0.8', 'CLM Default'), ('GEDI Max', 'CLM Default'), ('GEDI Mean', 'CLM Default'), ('GEDI Max', 'GEDI Mean')]:
    dL = np.nanmean(z[a]['QFLX_LIQ_GRND'][30:40] - z[b]['QFLX_LIQ_GRND'][30:40], 0)
    dSn = np.nanmean(z[a]['QFLX_SNOW_GRND'][30:40] - z[b]['QFLX_SNOW_GRND'][30:40], 0)
    for name, sel in [('global', f), ('NH>=50', f & (plat >= 50))]:
        l = np.nansum((dL[ci] * w)[sel]); sn = np.nansum((dSn[ci] * w)[sel])
        Q[f'snowshare|{a}|{b}|{name}'] = 100 * sn / (l + sn)
Q['snowshare_wb_range'] = [100 * min(v['frac_from_snow'] for v in J['an1']['canopy_water_balance_mm_yr'].values()),
                           100 * max(v['frac_from_snow'] for v in J['an1']['canopy_water_balance_mm_yr'].values())]
# energy: per unit forest area
Q['Rnet_rel|Default 1.2'] = rel('Default 1.2', 'CLM Default', 'Rnet')
# heights
H = {s: patch_height(s, m) for s in SCEN}
ret = f & (np.abs(H['GEDI Mean'] - H['CLM Default']) < 1e-4)
asg = f & ~ret
Q['retained|all'] = 100 * w[ret].sum() / w[f].sum()
for p, sel in PF.items():
    Q[f'retained|{p}'] = 100 * w[sel & ret].sum() / w[sel].sum()
for s in ['CLM Default', 'GEDI Max', 'GEDI Mean', 'GEDI Median']:
    Q[f'hS4|all|{s}'] = wm(H[s], asg)
    Q[f'hall|all|{s}'] = wm(H[s], f)
    for p, sel in PF.items():
        Q[f'hS4|{p}|{s}'] = wm(H[s], sel & asg)
for p in list(PF) + ['all']:
    for s in ['GEDI Max', 'GEDI Mean', 'GEDI Median']:
        Q[f'hS4pct|{p}|{s}'] = 100 * (Q[f'hS4|{p}|{s}'] / Q[f'hS4|{p}|CLM Default'] - 1)
        Q[f'hS4m|{p}|{s}'] = Q[f'hS4|{p}|{s}'] - Q[f'hS4|{p}|CLM Default']


def wq(x, ww, p):
    o = np.argsort(x); cw = np.cumsum(ww[o]) / ww.sum()
    return float(x[o][np.searchsorted(cw, p)])


for p, sel in PF.items():                             # Fig. 3e medians (area-weighted, incl. retained patches)
    for s in ['CLM Default', 'GEDI Max', 'GEDI Mean', 'GEDI Median']:
        Q[f'med3e|{p}|{s}'] = wq(H[s][sel], w[sel], 0.5)
        Q[f'mean3e|{p}|{s}'] = wm(H[s], sel)
# Fig. 3b-d text: cell-level forest-weighted heights; shares and regional means by forest area
ff = forest_frac(m); msk = ff > 0
Hg = {s: grid_sum(np.nan_to_num(H[s]), m) / np.where(msk, ff, np.nan) for s in H}
fa = ff * A                                            # forest area per cell
LAT = np.repeat(m['lat'][:, None], A.shape[1], 1)
rg = {s: 100 * (Hg[s] / Hg['CLM Default'] - 1) for s in ['GEDI Max', 'GEDI Mean', 'GEDI Median']}
Q['share_Max_gt5'] = 100 * np.nansum(fa[msk & (rg['GEDI Max'] > 5)]) / np.nansum(fa[msk])
regs = {'eq10': np.abs(LAT) <= 10, 'NH23.5-60': (LAT >= 23.5) & (LAT < 60), 'N60': LAT >= 60,
        'trop10-23.5': (np.abs(LAT) > 10) & (np.abs(LAT) < 23.5), 'SH<-23.5': LAT <= -23.5}
for s in rg:
    for rn, rs in regs.items():
        sel = msk & rs
        Q[f'reg3|{s}|{rn}|cellmean'] = float(np.nansum((rg[s] * fa)[sel]) / np.nansum(fa[sel]))
        Q[f'reg3|{s}|{rn}|ratio'] = float(100 * (np.nansum((Hg[s] * fa)[sel]) / np.nansum((Hg['CLM Default'] * fa)[sel]) - 1))
sel = msk & regs['N60']
Q['N60_Max_up_share'] = 100 * np.nansum(fa[sel & (rg['GEDI Max'] > 5)]) / np.nansum(fa[sel])
Q['N60_Max_down_share'] = 100 * np.nansum(fa[sel & (rg['GEDI Max'] < -5)]) / np.nansum(fa[sel])
# emulator: share of forest area beyond +-20% and clamped emulation
for s in ['GEDI Max', 'GEDI Mean', 'GEDI Median']:
    rr = H[s] / H['CLM Default']
    Q[f'beyond20|{s}'] = 100 * w[f & ((rr < 0.8) | (rr > 1.2))].sum() / w[f].sum()
rr = np.stack([H[s] / H['CLM Default'] for s in ['GEDI Max', 'GEDI Mean', 'GEDI Median']])
Q['beyond20|any'] = 100 * w[f & np.any((rr < 0.8) | (rr > 1.2), 0)].sum() / w[f].sum()
emu = np.load('emulator_coef.npz')
for clamp in ['const', 'linear']:
    for a, b in [('GEDI Mean', 'CLM Default'), ('GEDI Max', 'CLM Default'), ('GEDI Median', 'CLM Default'), ('GEDI Max', 'GEDI Mean')]:
        y = {}
        for s in (a, b):
            lr = np.log(np.where(f, H[s] / H['CLM Default'], 1.0)); lc = np.clip(lr, np.log(0.6), np.log(1.6))
            tot = 0
            for v in ['FCEV', 'FCTR', 'FGEV']:
                c = emu[v]
                yy = c[0] + c[1] * lc + c[2] * lc ** 2 + (c[1] * (lr - lc) if clamp == 'linear' else 0)
                tot = tot + np.where(f, yy, 0)
            y[s] = tot
        em = np.nansum(((y[a] - y[b]) * w)[f]) / np.nansum((S[b]['ET'] * w)[f]) * 100
        Q[f'emu|{clamp}|{a}|{b}'] = float(em)
        Q[f'emudev|{clamp}|{a}|{b}'] = float(100 * (em / Q[f'rel|{a}|{b}|ET'] - 1))
Q['emu_r_cells'] = [min(v2['r_cells'] for v1 in J['an4']['emulator_validation'].values() for v2 in v1.values()),
                    max(v2['r_cells'] for v1 in J['an4']['emulator_validation'].values() for v2 in v1.values())]
Q['emu_json_ET'] = {k: (v['ET']['sim'], v['ET']['emulated']) for k, v in J['an4']['emulator_validation'].items()}
# significance and robustness
g = J['an3']['global']
Q['max_q_neff_global'] = max(v['q_neff_t'] for v in g.values())
Q['q_block|GEDI Max vs CLM Default|ET'] = g['GEDI Max vs CLM Default|ET']['q_block_wilcoxon']
Q['halfwidth_ET_GEDIvsDef'] = [(g[f'{s} vs CLM Default|ET']['ci_hi'] - g[f'{s} vs CLM Default|ET']['ci_lo']) / 2 for s in ['GEDI Mean', 'GEDI Max', 'GEDI Median']]
Q['halfwidth_all'] = [min((v['ci_hi'] - v['ci_lo']) / 2 for v in g.values()), max((v['ci_hi'] - v['ci_lo']) / 2 for v in g.values())]
Q['pft_sig_counts'] = {t: sum(1 for p in J['an3'] if p != 'global' for k in J['an3'][p] if J['an3'][p][k][t] < 0.05) for t in ['q_neff_t', 'q_block_wilcoxon', 'q_cell_wilcoxon']}
fy = J['an1']['forcing_year_robustness']
Q['fy_max_sd'] = max(fy[k][v]['sd'] for k in fy for v in fy[k])
Q['fy_all_same_sign'] = all(fy[k][v]['same_sign'] == 10 for k in fy for v in fy[k])
Q['fy_max_c3_minus_mean'] = max(abs(fy[k][v]['cycle3'] - fy[k][v]['mean']) for k in fy for v in fy[k])
dr = J['an1']['drift_cycle4_vs_cycle3_pct']
Q['drift_ET_range'] = [min(dr[s]['ET'] for s in dr), max(dr[s]['ET'] for s in dr)]
Q['spread_ET_pp'] = g['GEDI Mean vs CLM Default|ET']['rel_pct'] - g['GEDI Max vs CLM Default|ET']['rel_pct']
# block-size sensitivity of the block Wilcoxon test and bootstrap (global family, BH within family)
dG = {}
for a, b in COMPS:
    for v in ['FCEV', 'FCTR', 'FGEV', 'ET']:
        dG[(a, b, v)] = (grid_sum(S[a][v], m), grid_sum(S[b][v], m))
nlat, nlon = A.shape
rng = np.random.default_rng(0)
for bl, (BLAT, BLON) in {'5deg': (3, 2), '10deg': (5, 4), '20deg': (11, 8)}.items():
    bid = (np.arange(nlat)[:, None] // BLAT) * 1000 + (np.arange(nlon)[None, :] // BLON)
    blocks = np.unique(bid[msk])
    P = {}; CI = {}
    for (a, b, v), (ga, gb) in dG.items():
        d = ga - gb
        bd = np.array([np.sum((d * A)[(bid == k) & msk]) for k in blocks]); bb = np.array([np.sum((gb * A)[(bid == k) & msk]) for k in blocks])
        ba = np.array([np.sum(A[(bid == k) & msk]) for k in blocks])
        P[(a, b, v)] = stats.wilcoxon(bd / np.where(ba > 0, ba, np.nan)).pvalue
        idx = rng.integers(0, len(blocks), size=(2000, len(blocks)))
        bs = 100 * bd[idx].sum(1) / bb[idx].sum(1)
        CI[(a, b, v)] = np.percentile(bs, [2.5, 97.5])
    keys = list(P); p = np.array([P[k] for k in keys]); o = np.argsort(p); n = len(p)
    q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    for k, qq in zip(keys, q):
        Q[f'blk|{bl}|{k[0]} vs {k[1]}|{k[2]}'] = [float(P[k]), float(min(qq, 1)), float(CI[k][0]), float(CI[k][1])]
    Q[f'blk|{bl}|nblocks'] = int(len(blocks))
# potential-temperature pathway (Text S3 / A20)
R = np.where((iv == 4) | (iv == 5), 0.075, 0.055) + 0.67
for a in ['Default 1.2', 'Default 0.8']:
    Q[f'dtheta|{a}'] = wm(0.0098 * R * (H[a] - H['CLM Default']), f)
json.dump(Q, open('text_numbers.json', 'w'), indent=0, default=float)
print('wrote', len(Q))
