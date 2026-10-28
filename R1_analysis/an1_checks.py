"""R1 analysis 1: snow-burial check, spin-up drift, forcing-year robustness,
canopy water balance closure, energy diagnostics."""
import numpy as np, json
from common import *

m = load_meta()
mm = np.load(R1 + 'meta.npz')
assert np.array_equal(mm['pfts1d_itype_veg'], m['pfts1d_itype_veg'])
f = m['forest']
names = list(SCEN)
cases = {s: SCEN[s][2] for s in names}
V = ['FCEV', 'FCTR', 'FGEV', 'QVEGE', 'QVEGT', 'QFLX_EVAP_TOT', 'QFLX_LIQ_GRND', 'QFLX_SNOW_GRND',
     'ELAI', 'ESAI', 'FSH', 'Rnet', 'TV', 'FSH_V', 'FSH_G', 'SABV', 'FIRE']
A = {}
for s in names:
    z = np.load(R1 + cases[s] + '_annual.npz')
    A[s] = {v: z[v].astype(np.float64) for v in V}
years = np.arange(2001, 2041)
res = {}

# 1. ELAI/ESAI identical across scenarios? (snow burial depends on htop/hbot)
for v in ['ELAI', 'ESAI']:
    dmax = {s: float(np.nanmax(np.abs(A[s][v][30:, f] - A['CLM Default'][v][30:, f]))) for s in names}
    res[v + '_maxabsdiff_vs_default'] = dmax
print('ELAI/ESAI max |diff| vs Default (analysis decade, forest patches):', res['ELAI_maxabsdiff_vs_default'], res['ESAI_maxabsdiff_vs_default'])

# helper: global forest mean (manuscript convention) for an annual (pft,) vector
def G(vec):
    return gmean(grid_sum(vec, m), m)

# 2. spin-up drift: cycle 3 (2021-2030) vs cycle 4 (2031-2040) climatological global means
drift = {}
for s in names:
    drift[s] = {}
    for v in ['FCEV', 'FCTR', 'FGEV']:
        c3 = G(np.nanmean(A[s][v][20:30], axis=0)); c4 = G(np.nanmean(A[s][v][30:40], axis=0))
        drift[s][v] = 100 * (c4 - c3) / c3
    et3 = G(np.nanmean(A[s]['FCEV'][20:30] + A[s]['FCTR'][20:30] + A[s]['FGEV'][20:30], axis=0))
    et4 = G(np.nanmean(A[s]['FCEV'][30:40] + A[s]['FCTR'][30:40] + A[s]['FGEV'][30:40], axis=0))
    drift[s]['ET'] = 100 * (et4 - et3) / et3
res['drift_cycle4_vs_cycle3_pct'] = drift
print('Drift cycle4 vs cycle3 (%):', {s: {k: round(v, 3) for k, v in d.items()} for s, d in drift.items()})

# also: scenario difference in cycle 3 vs cycle 4 (stability of the signal)
combos = [('GEDI Mean', 'CLM Default'), ('GEDI Max', 'CLM Default'), ('GEDI Max', 'GEDI Mean'),
          ('GEDI Median', 'GEDI Mean'), ('Default 0.8', 'CLM Default'), ('Default 0.9', 'CLM Default'),
          ('Default 1.1', 'CLM Default'), ('Default 1.2', 'CLM Default')]

def et(s, sl):
    return A[s]['FCEV'][sl] + A[s]['FCTR'][sl] + A[s]['FGEV'][sl]

# 3. forcing-year robustness: relative difference for each of the 10 analysis years
yr_stats = {}
for a, b in combos:
    yr_stats[f'{a} vs {b}'] = {}
    for v in ['FCEV', 'FCTR', 'FGEV', 'ET']:
        vals = []
        for iy in range(30, 40):
            xa = A[a][v][iy] if v != 'ET' else et(a, iy)
            xb = A[b][v][iy] if v != 'ET' else et(b, iy)
            ga, gb = G(xa), G(xb)
            vals.append(100 * (ga - gb) / gb)
        vals = np.array(vals)
        c3a = G(np.nanmean(A[a][v][20:30], 0) if v != 'ET' else np.nanmean(et(a, slice(20, 30)), 0))
        c3b = G(np.nanmean(A[b][v][20:30], 0) if v != 'ET' else np.nanmean(et(b, slice(20, 30)), 0))
        yr_stats[f'{a} vs {b}'][v] = dict(min=float(vals.min()), max=float(vals.max()),
                                          mean=float(vals.mean()), sd=float(vals.std(ddof=1)),
                                          same_sign=int(np.sum(np.sign(vals) == np.sign(vals.mean()))),
                                          cycle3=float(100 * (c3a - c3b) / c3b))
res['forcing_year_robustness'] = yr_stats
for k, d in yr_stats.items():
    print(k, {v: f"[{x['min']:+.2f},{x['max']:+.2f}] sd={x['sd']:.2f} sign={x['same_sign']}/10 c3={x['cycle3']:+.2f}" for v, x in d.items()})

# 4. canopy water balance (column scale). Forcing precipitation is identical across runs, and
#    canopy storage change ~0 over a decade, so Delta(water reaching ground) = -Delta(sum wtcol*E_can)
ci = m['ci']
ncol = len(m['cols1d_ixy'])
def col_sum(vec_patch):
    out = np.zeros(ncol); np.add.at(out, ci, vec_patch * m['wtcol']); return out
# natural vegetation columns that contain forest patches
colw = m['cols1d_wtgcell']
has_forest = np.zeros(ncol, bool); has_forest[ci[f]] = True
wb = {}
for a, b in combos:
    dE = col_sum(np.nanmean(A[a]['QVEGE'][30:40] - A[b]['QVEGE'][30:40], 0)) * SEC_YR
    dL = np.nanmean(A[a]['QFLX_LIQ_GRND'][30:40] - A[b]['QFLX_LIQ_GRND'][30:40], 0) * SEC_YR
    dS = np.nanmean(A[a]['QFLX_SNOW_GRND'][30:40] - A[b]['QFLX_SNOW_GRND'][30:40], 0) * SEC_YR
    sel = has_forest & np.isfinite(dL) & np.isfinite(dS)
    # area-weighted (column area = colw * cell area)
    jx = m['cols1d_jxy'].astype(int) - 1; ix = m['cols1d_ixy'].astype(int) - 1
    w = colw * m['area'][jx, ix]; w = np.where(np.isfinite(w), w, 0) * sel
    dEm, dLm, dSm = [np.nansum(x * w) / np.nansum(w) for x in (dE, dL, dS)]
    # closure across columns
    resid = dE + dL + dS
    r = np.corrcoef(dE[sel], -(dL + dS)[sel])[0, 1]
    wb[f'{a} vs {b}'] = dict(dEcan=dEm, dLiqGround=dLm, dSnowGround=dSm, residual=dEm + dLm + dSm,
                             r_columns=float(r), frac_from_snow=float(-dSm / dEm) if dEm != 0 else np.nan,
                             resid_p99=float(np.nanpercentile(np.abs(resid[sel]), 99)))
res['canopy_water_balance_mm_yr'] = wb
for k, d in wb.items():
    print('WB', k, {x: round(y, 4) for x, y in d.items()})

# 5. energy diagnostics for the scaling runs (global forest means, analysis decade)
en = {}
for s in names:
    en[s] = {v: G(np.nanmean(A[s][v][30:40], 0)) for v in ['FSH', 'Rnet', 'TV', 'FSH_V', 'FSH_G', 'SABV', 'FIRE']}
    en[s]['LH'] = G(np.nanmean(et(s, slice(30, 40)), 0))
res['energy_global'] = en
for s in names:
    print('EN', s, {k: round(v, 3) for k, v in en[s].items()})
json.dump(res, open('an1_results.json', 'w'), indent=1)
