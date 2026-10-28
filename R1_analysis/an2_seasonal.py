"""R1 analysis 2: snow-burial extent, seasonal/regional decomposition of flux changes,
canopy temperature diagnostics, implied runoff change."""
import numpy as np, json
from common import *

m = load_meta(); f = m['forest']
names = list(SCEN); cases = {s: SCEN[s][2] for s in names}
V = ['FCEV', 'FCTR', 'FGEV', 'QVEGE', 'QVEGT', 'QFLX_EVAP_TOT', 'QFLX_LIQ_GRND', 'QFLX_SNOW_GRND',
     'ELAI', 'ESAI', 'TV', 'TG', 'FSH', 'Rnet']
M = {}
for s in names:
    z = np.load(R1 + cases[s] + '_monclim.npz')
    M[s] = {v: z[v].astype(np.float64) for v in V}
res = {}
plat = m['pfts1d_lat']
ffrac = forest_frac(m)

# --- 1. snow burial: ELAI change vs height ------------------------------------------
hdef = patch_height('CLM Default', m)
for s in ['Default 0.8', 'Default 1.2', 'GEDI Mean', 'GEDI Max', 'GEDI Median']:
    h = patch_height(s, m)
    dl = np.nanmax(np.abs(M[s]['ELAI'] - M['CLM Default']['ELAI']), axis=0)  # max over months
    aff = f & (dl > 0.01)
    w = m['wt'] * m['area'][m['jx'], m['ix']]
    res[f'burial_{s}'] = dict(n_patches=int(aff.sum()), n_forest=int(f.sum()),
                              area_frac=float(np.nansum(w[aff]) / np.nansum(w[f])),
                              htop_scen_median=float(np.nanmedian(h[aff])) if aff.any() else None,
                              htop_def_median=float(np.nanmedian(hdef[aff])) if aff.any() else None,
                              lat_median=float(np.nanmedian(plat[aff])) if aff.any() else None)
    ann_dl = np.nanmean(M[s]['ELAI'] - M['CLM Default']['ELAI'], 0)
    res[f'burial_{s}']['global_mean_dELAI'] = float(gmean(grid_sum(ann_dl, m), m) / gmean(ffrac, m))
print({k: v for k, v in res.items() if k.startswith('burial')})

# --- 2. seasonal/regional decomposition -------------------------------------------
regions = {'Tropics (|lat|<23.5)': np.abs(plat) < 23.5,
           'Temperate (23.5-50)': (np.abs(plat) >= 23.5) & (np.abs(plat) < 50),
           'Boreal (>=50)': np.abs(plat) >= 50}
seasons = {'DJF': [11, 0, 1], 'MAM': [2, 3, 4], 'JJA': [5, 6, 7], 'SON': [8, 9, 10]}
def region_mean_monthly(s, v, sel):
    x = M[s][v]  # (12, pft)
    out = []
    for k in range(12):
        out.append(gmean(grid_sum(x[k], m, sel), m))
    return np.array(out)
seas = {}
for a, b in [('Default 1.2', 'CLM Default'), ('GEDI Max', 'GEDI Mean'), ('GEDI Mean', 'CLM Default')]:
    seas[f'{a} vs {b}'] = {}
    for rn, rsel in regions.items():
        sel = f & rsel
        d = {}
        for v, nm in [('QVEGE', 'Ecan'), ('QVEGT', 'T')]:
            da = region_mean_monthly(a, v, sel) - region_mean_monthly(b, v, sel)
            d[nm + '_monthly_mm'] = (da * 86400 * 365 / 12).tolist()   # mm per (average) month
            d[nm + '_annual_mm'] = float(da.mean() * SEC_YR)
            # share of annual change by season (NH seasons; for SH tropics this is calendar season)
            d[nm + '_season_share'] = {sn: float(da[idx].sum() / da.sum()) for sn, idx in seasons.items()}
        ge = (region_mean_monthly(a, 'QFLX_EVAP_TOT', sel) - region_mean_monthly(a, 'QVEGE', sel) - region_mean_monthly(a, 'QVEGT', sel)) - \
             (region_mean_monthly(b, 'QFLX_EVAP_TOT', sel) - region_mean_monthly(b, 'QVEGE', sel) - region_mean_monthly(b, 'QVEGT', sel))
        d['Egnd_annual_mm'] = float(ge.mean() * SEC_YR)
        base_e = region_mean_monthly(b, 'QVEGE', sel).mean() * SEC_YR
        d['Ecan_base_annual_mm'] = float(base_e)
        d['Ecan_rel_pct'] = 100 * d['Ecan_annual_mm'] / base_e
        seas[f'{a} vs {b}'][rn] = d
res['seasonal'] = seas
for k, d in seas.items():
    for rn, x in d.items():
        print(k, rn, f"dEcan={x['Ecan_annual_mm']:+.3f} mm/yr ({x['Ecan_rel_pct']:+.2f}%), dT={x['T_annual_mm']:+.3f}, dEgnd={x['Egnd_annual_mm']:+.3f}",
              'Ecan season share', {s: round(v, 2) for s, v in x['Ecan_season_share'].items()})

# BDT Bor specifically: monthly change in canopy evaporation, transpiration (W/m2, forest-PFT mean)
def pft_monthly(s, v, pid):
    sel = f & (m['iv'] == pid)
    w = m['wt'] * m['area'][m['jx'], m['ix']]
    x = M[s][v][:, sel]
    return (x * w[sel]).sum(1) / w[sel].sum()
bdt = {}
for a, b in [('GEDI Mean', 'CLM Default'), ('Default 1.2', 'CLM Default'), ('GEDI Max', 'GEDI Mean')]:
    bdt[f'{a} vs {b}'] = {v: (pft_monthly(a, v, 8) - pft_monthly(b, v, 8)).round(3).tolist() for v in ['FCEV', 'FCTR', 'FGEV']}
    bdt[f'{a} vs {b}']['base_FCEV'] = pft_monthly(b, 'FCEV', 8).round(2).tolist()
    bdt[f'{a} vs {b}']['base_ELAI'] = pft_monthly(b, 'ELAI', 8).round(2).tolist()
    bdt[f'{a} vs {b}']['base_ESAI'] = pft_monthly(b, 'ESAI', 8).round(2).tolist()
res['BDT_Bor_monthly'] = bdt
for k, d in bdt.items():
    print('BDT Bor', k)
    for v, x in d.items():
        print('   ', v, x)

# --- 3. canopy temperature diagnostics (forest-weighted means) ----------------------
w = m['wt'] * m['area'][m['jx'], m['ix']]
def fw(vec, sel=f):
    return float(np.nansum(vec[sel] * w[sel]) / np.nansum(w[sel]))
tv = {s: fw(np.nanmean(M[s]['TV'], 0)) for s in names}
res['TV_forest_mean_K'] = tv
print('TV forest mean (K):', {s: round(v - tv['CLM Default'], 4) for s, v in tv.items()})

# --- 4. implied runoff+drainage change = -dET (forcing precipitation identical) ------
ro = {}
for a, b in [('Default 0.8', 'CLM Default'), ('Default 1.2', 'CLM Default'), ('GEDI Mean', 'CLM Default'),
             ('GEDI Max', 'CLM Default'), ('GEDI Max', 'GEDI Mean')]:
    det = np.nanmean(M[a]['QFLX_EVAP_TOT'] - M[b]['QFLX_EVAP_TOT'], 0) * SEC_YR
    # per unit forest area
    ro[f'{a} vs {b}'] = dict(dET_mm_per_forest_area=fw(det), implied_dR=-fw(det),
                             base_ET_mm=fw(np.nanmean(M[b]['QFLX_EVAP_TOT'], 0) * SEC_YR))
res['implied_runoff'] = ro
print('implied runoff:', ro)
json.dump(res, open('an2_results.json', 'w'), indent=1)
