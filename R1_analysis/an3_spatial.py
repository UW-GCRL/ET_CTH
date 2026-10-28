"""R1 analysis 3: spatial-autocorrelation-aware significance testing.
For each scenario pair and flux, grid-cell differences (forest flux per unit grid-cell area, as in
the global means) are tested with:
  (1) the original cell-level Wilcoxon signed-rank test (for reference);
  (2) Moran's I of the differences and an effective sample size n_eff (Dutilleul-type
      approximation n_eff = n (1-rho)/(1+rho) with rho the lag-1 neighbour correlation);
  (3) a spatial block bootstrap (blocks of ~10x10 degrees, i.e. 5x4 grid cells) of the
      area-weighted global relative difference -> 95% CI;
  (4) a Wilcoxon test on block-aggregated differences (one value per block).
Global and per-PFT results are written to an3_results.json / CSV.
"""
import numpy as np, json, csv
from scipy import stats
from common import *

rng = np.random.default_rng(20261019)
m = load_meta(); f = m['forest']
D = {s: load_scen(s) for s in SCEN}
area = np.where(np.isfinite(m['area']), m['area'], 0.0)
nlat, nlon = area.shape
lat = m['lat']; lon = m['lon']
combos = [('GEDI Mean', 'CLM Default'), ('GEDI Max', 'CLM Default'), ('GEDI Median', 'CLM Default'),
          ('GEDI Max', 'GEDI Mean'), ('GEDI Median', 'GEDI Mean'), ('Default 0.8', 'CLM Default'),
          ('Default 0.9', 'CLM Default'), ('Default 1.1', 'CLM Default'), ('Default 1.2', 'CLM Default')]
VARS4 = ['FCEV', 'FCTR', 'FGEV', 'ET']
BLAT, BLON = 5, 4   # 5 x 1.89 = 9.5 deg lat, 4 x 2.5 = 10 deg lon
block_id = (np.arange(nlat)[:, None] // BLAT) * 1000 + (np.arange(nlon)[None, :] // BLON)


def morans_I(x, mask):
    """Moran's I with rook+queen (8-neighbour) binary weights on the lat-lon grid, longitude periodic."""
    z = np.where(mask, x - x[mask].mean(), 0.0)
    num = 0.0; wsum = 0.0
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di == 0 and dj == 0:
                continue
            zs = np.roll(np.roll(z, di, 0), dj, 1)
            ms = np.roll(np.roll(mask, di, 0), dj, 1)
            if di == 1: ms[0, :] = False
            if di == -1: ms[-1, :] = False
            pair = mask & ms
            num += np.sum(z[pair] * zs[pair]); wsum += pair.sum()
    n = mask.sum()
    I = (n / wsum) * num / np.sum(z[mask] ** 2)
    return float(I), int(n)


def neighbour_corr(x, mask):
    """lag-1 correlation between each cell and the mean of its valid 8 neighbours."""
    s = np.zeros_like(x); c = np.zeros_like(x)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di == 0 and dj == 0: continue
            xs = np.roll(np.roll(np.where(mask, x, 0.0), di, 0), dj, 1)
            ms = np.roll(np.roll(mask, di, 0), dj, 1).astype(float)
            if di == 1: ms[0, :] = 0
            if di == -1: ms[-1, :] = 0
            s += xs * ms; c += ms
    ok = mask & (c > 0)
    return float(np.corrcoef(x[ok], (s[ok] / c[ok]))[0, 1])


def analyse(sel_patch, label):
    out = {}
    fr = forest_frac(m, sel_patch)
    mask = fr > 0
    blocks = np.unique(block_id[mask])
    for a, b in combos:
        for v in VARS4:
            ga = grid_sum(D[a][v], m, sel_patch); gb = grid_sum(D[b][v], m, sel_patch)
            d = ga - gb
            rel = 100 * np.sum(d[mask] * area[mask]) / np.sum(gb[mask] * area[mask])
            # (1) original cell-level Wilcoxon
            try:
                p_cell = stats.wilcoxon(d[mask]).pvalue
            except ValueError:
                p_cell = np.nan
            # (2) autocorrelation
            I, n = morans_I(d, mask)
            rho = neighbour_corr(d, mask)
            rho1 = max(min(rho, 0.999), 0.0)
            n_eff = n * (1 - rho1) / (1 + rho1)
            # t-test on the area-weighted mean difference with n_eff
            wts = area[mask] / area[mask].sum()
            mu = np.sum(wts * d[mask]); sd = np.sqrt(np.sum(wts * (d[mask] - mu) ** 2))
            t = mu / (sd / np.sqrt(max(n_eff, 2)))
            p_neff = 2 * stats.t.sf(abs(t), df=max(n_eff - 1, 1))
            # (3) block bootstrap of the global relative difference
            bsum_d = np.array([np.sum((d * area)[(block_id == k) & mask]) for k in blocks])
            bsum_b = np.array([np.sum((gb * area)[(block_id == k) & mask]) for k in blocks])
            idx = rng.integers(0, len(blocks), size=(2000, len(blocks)))
            boot = 100 * bsum_d[idx].sum(1) / bsum_b[idx].sum(1)
            lo, hi = np.percentile(boot, [2.5, 97.5])
            # (4) Wilcoxon across block means (area-weighted mean difference per block)
            barea = np.array([np.sum(area[(block_id == k) & mask]) for k in blocks])
            bmean = bsum_d / barea
            try:
                p_block = stats.wilcoxon(bmean).pvalue
            except ValueError:
                p_block = np.nan
            frac_pos = float(np.mean(d[mask] > 0))
            out[f'{a} vs {b}|{v}'] = dict(rel_pct=float(rel), ci_lo=float(lo), ci_hi=float(hi),
                                          p_cell_wilcoxon=float(p_cell), morans_I=I, rho_nb=rho,
                                          n=n, n_eff=float(n_eff), p_neff_t=float(p_neff),
                                          n_blocks=int(len(blocks)), p_block_wilcoxon=float(p_block),
                                          frac_cells_positive=frac_pos)
    return out


res = {'global': analyse(f, 'global')}
for pn, pid in PFTS.items():
    res[pn] = analyse(f & (m['iv'] == pid), pn)
json.dump(res, open('an3_results.json', 'w'), indent=1)

# raw p-values -> BH-FDR within each family (global; each PFT) for the cell, n_eff and block tests
from statsmodels.stats.multitest import multipletests
rows = []
for fam, d in res.items():
    keys = list(d)
    for test in ['p_cell_wilcoxon', 'p_neff_t', 'p_block_wilcoxon']:
        p = np.array([d[k][test] for k in keys]); ok = np.isfinite(p)
        q = np.full(p.shape, np.nan); q[ok] = multipletests(p[ok], method='fdr_bh')[1]
        for k, qq in zip(keys, q):
            d[k][test.replace('p_', 'q_')] = float(qq)
    for k in keys:
        comp, v = k.split('|')
        rows.append([fam, comp, v] + [d[k][x] for x in ['rel_pct', 'ci_lo', 'ci_hi', 'frac_cells_positive', 'morans_I', 'n', 'n_eff',
                                                        'q_cell_wilcoxon', 'q_neff_t', 'n_blocks', 'q_block_wilcoxon']])
with open('an3_results.csv', 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['domain', 'comparison', 'variable', 'rel_pct', 'ci95_lo', 'ci95_hi', 'frac_cells_positive', 'morans_I',
                'n_cells', 'n_eff', 'q_cell_wilcoxon', 'q_neff_t', 'n_blocks', 'q_block_wilcoxon'])
    w.writerows(rows)
json.dump(res, open('an3_results.json', 'w'), indent=1)
for k, x in res['global'].items():
    print(f"{k:38s} rel={x['rel_pct']:+.2f} CI=[{x['ci_lo']:+.2f},{x['ci_hi']:+.2f}] I={x['morans_I']:.2f} n={x['n']} neff={x['n_eff']:.0f} "
          f"q_cell={x['q_cell_wilcoxon']:.1e} q_neff={x['q_neff_t']:.1e} q_block={x['q_block_wilcoxon']:.1e} (nb={x['n_blocks']}) pos={x['frac_cells_positive']:.2f}")
