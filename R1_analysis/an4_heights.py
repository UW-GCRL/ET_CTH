"""R1 analysis 4: canopy-height inputs on the model grid.
(a) Default vs GEDI Max/Mean/Median per PFT (area-weighted mean, bias, RMSD, correlation, where
    the Default lies between Mean and Max);
(b) 1-m fill values (PFT present in CLM surface data but absent from the 500-m PFT map);
(c) patch-level height-flux emulator from the uniform-scaling runs, validated against the
    GEDI simulations.
"""
import numpy as np, json
from common import *

m = load_meta(); f = m['forest']
w = m['wt'] * m['area'][m['jx'], m['ix']]
H = {s: patch_height(s, m) for s in SCEN}
res = {}

# (a) height statistics per PFT
stats_h = {}
for pn, pid in PFTS.items():
    sel = f & (m['iv'] == pid)
    ww = w[sel] / w[sel].sum()
    d = {}
    for s in ['CLM Default', 'GEDI Max', 'GEDI Mean', 'GEDI Median']:
        d[s] = dict(wmean=float(np.sum(H[s][sel] * ww)), median=float(np.median(H[s][sel])))
    for s in ['GEDI Max', 'GEDI Mean', 'GEDI Median']:
        diff = H[s][sel] - H['CLM Default'][sel]
        d[s].update(bias=float(np.sum(diff * ww)), rmsd=float(np.sqrt(np.sum(diff ** 2 * ww))),
                    r=float(np.corrcoef(H[s][sel], H['CLM Default'][sel])[0, 1]))
    # position of Default between Mean (0) and Max (1)
    mx, mn, df = d['GEDI Max']['wmean'], d['GEDI Mean']['wmean'], d['CLM Default']['wmean']
    d['position_default_between_mean0_max1'] = float((df - mn) / (mx - mn))
    # fraction of patches where Default is closer to Max than to Mean
    closer_max = np.abs(H['CLM Default'][sel] - H['GEDI Max'][sel]) < np.abs(H['CLM Default'][sel] - H['GEDI Mean'][sel])
    d['frac_patch_area_closer_to_max'] = float(np.sum(ww[closer_max]))
    stats_h[pn] = d
# all forest
ww = w[f] / w[f].sum()
allf = {s: float(np.sum(H[s][f] * ww)) for s in ['CLM Default', 'GEDI Max', 'GEDI Mean', 'GEDI Median']}
allf['position_default_between_mean0_max1'] = (allf['CLM Default'] - allf['GEDI Mean']) / (allf['GEDI Max'] - allf['GEDI Mean'])
closer_max = np.abs(H['CLM Default'][f] - H['GEDI Max'][f]) < np.abs(H['CLM Default'][f] - H['GEDI Mean'][f])
allf['frac_area_closer_to_max'] = float(np.sum(ww[closer_max]))
stats_h['All forest'] = allf
res['height_stats'] = stats_h
for k, v in stats_h.items():
    print(k, json.dumps(v)[:600])

# (b) fill values
fill = {}
for s in ['GEDI Max', 'GEDI Mean', 'GEDI Median']:
    isfill = f & (np.abs(H[s] - 1.0) < 1e-6)
    fill[s] = dict(n_patches=int(isfill.sum()), area_frac=float(w[isfill].sum() / w[f].sum()),
                   by_pft={pn: int((isfill & (m['iv'] == pid)).sum()) for pn, pid in PFTS.items()})
res['fill_1m'] = fill
print('fill:', fill)

# (c) emulator: per patch, fit flux(k) = a + b ln k + c (ln k)^2 using scaling runs k=0.8..1.2
D = {s: load_scen(s) for s in SCEN}
ks = {'Default 0.8': 0.8, 'Default 0.9': 0.9, 'CLM Default': 1.0, 'Default 1.1': 1.1, 'Default 1.2': 1.2}
X = np.log(np.array(list(ks.values())))
Amat = np.vstack([np.ones_like(X), X, X ** 2]).T
pinv = np.linalg.pinv(Amat)
emu = {}
for v in ['FCEV', 'FCTR', 'FGEV']:
    Y = np.vstack([D[s][v] for s in ks])          # (5, pft)
    coef = pinv @ np.nan_to_num(Y)                # (3, pft)
    emu[v] = coef
# linearity check: fraction of variance explained by linear term across the 5 runs (global)
lin = {}
for v in ['FCEV', 'FCTR', 'FGEV']:
    g = np.array([gmean(grid_sum(D[s][v], m), m) for s in ks])
    p1 = np.polyfit(X, g, 1); p2 = np.polyfit(X, g, 2)
    lin[v] = dict(elasticity_pct_per_lnh=float(100 * p1[0] / g[2]), curvature_rel=float(p2[0] / p2[1]),
                  max_abs_resid_linear_pct=float(np.max(np.abs(np.polyval(p1, X) - g)) / g[2] * 100))
res['global_linearity'] = lin
print('global linearity:', lin)

pred = {}
for s in ['GEDI Max', 'GEDI Mean', 'GEDI Median']:
    lr = np.log(np.where(f, H[s] / H['CLM Default'], 1.0))
    lr_c = np.clip(lr, np.log(0.6), np.log(1.6))  # limit quadratic extrapolation
    pred[s] = {}
    for v in ['FCEV', 'FCTR', 'FGEV']:
        c = emu[v]
        yhat = c[0] + c[1] * lr_c + c[2] * lr_c ** 2
        # linear-only extrapolation beyond the clip range
        yhat = yhat + c[1] * (lr - lr_c)
        pred[s][v] = np.where(f, yhat, D['CLM Default'][v])
    pred[s]['ET'] = pred[s]['FCEV'] + pred[s]['FCTR'] + pred[s]['FGEV']
val = {}
for a, b in [('GEDI Mean', 'CLM Default'), ('GEDI Max', 'CLM Default'), ('GEDI Median', 'CLM Default'), ('GEDI Max', 'GEDI Mean')]:
    val[f'{a} vs {b}'] = {}
    for v in ['FCEV', 'FCTR', 'FGEV', 'ET']:
        sim = 100 * (gmean(grid_sum(D[a][v], m), m) - gmean(grid_sum(D[b][v], m), m)) / gmean(grid_sum(D[b][v], m), m)
        pa = pred[a][v]; pb = pred[b][v] if b != 'CLM Default' else D['CLM Default'][v]
        emu_ = 100 * (gmean(grid_sum(pa, m), m) - gmean(grid_sum(pb, m), m)) / gmean(grid_sum(pb, m), m)
        # cell-level skill
        ga = grid_sum(D[a][v], m) - grid_sum(D[b][v], m); ge = grid_sum(pa, m) - grid_sum(pb, m)
        msk = forest_frac(m) > 0
        r = np.corrcoef(ga[msk], ge[msk])[0, 1]
        val[f'{a} vs {b}'][v] = dict(sim=float(sim), emulated=float(emu_), r_cells=float(r))
res['emulator_validation'] = val
for k, d in val.items():
    print('EMU', k, {v: f"sim={x['sim']:+.2f} emu={x['emulated']:+.2f} r={x['r_cells']:.2f}" for v, x in d.items()})
# per-PFT elasticities (% per unit ln h) for later use
el = {}
for pn, pid in PFTS.items():
    sel = f & (m['iv'] == pid)
    el[pn] = {}
    for v in ['FCEV', 'FCTR', 'FGEV', 'ET']:
        g = np.array([np.sum((D[s][v] * w)[sel]) for s in ks])
        el[pn][v] = float(100 * np.polyfit(X, g, 1)[0] / g[2])
g = {v: np.array([gmean(grid_sum(D[s][v], m), m) for s in ks]) for v in ['FCEV', 'FCTR', 'FGEV', 'ET']}
el['All forest'] = {v: float(100 * np.polyfit(X, g[v], 1)[0] / g[v][2]) for v in g}
res['elasticities_pct_per_lnh'] = el
print('elasticities:', el)
json.dump(res, open('an4_results.json', 'w'), indent=1)
np.savez('emulator_coef.npz', **{v: emu[v] for v in emu})
