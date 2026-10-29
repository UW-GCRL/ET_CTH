"""Figures 3 and 5 of the revised manuscript (AFM R1).
Conventions:
 * Figure 3: map values are the forest-area-weighted mean canopy height of the forest PFTs in each grid cell,
   drawn cell by cell (no interpolation) on the CLM5 1.9x2.5 grid;
 * Figure 5: relative differences of the PFT-area-weighted decade means, from CLM5's native water-flux outputs
   (QVEGE, QVEGT, and QFLX_EVAP_TOT - QVEGE - QVEGT for ground evaporation); whiskers are the area-weighted
   interquartile range of the grid-cell relative differences.
The maps of Figures 2 and 4 and Figs. S2 and S3 are drawn by scripts in the repository root (see the README).
"""
import os, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import cartopy.crs as ccrs
from common import *

OUT = FIG
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.linewidth': 0.6})

m = load_meta(); f = m['forest']
D = {s: load_scen(s) for s in SCEN}
ff = forest_frac(m)
mask = ff > 0
lat, lon = m['lat'], m['lon']
# cell edges (CLM 1.9x2.5: pole-centred latitudes)
late = np.concatenate([[-90], (lat[:-1] + lat[1:]) / 2, [90]])
lone = np.concatenate([lon - 1.25, [lon[-1] + 1.25]])
lone = np.where(lone > 180, lone - 360, lone)
# reorder longitudes to -180..180
order = np.argsort(np.where(lon > 180, lon - 360, lon))
lon_s = np.where(lon > 180, lon - 360, lon)[order]
lone_s = np.concatenate([lon_s - 1.25, [lon_s[-1] + 1.25]])

FLX = {'ET': 'ET_w', 'FCEV': 'E_can', 'FCTR': 'T', 'FGEV': 'E_gnd'}
LBL = {'ET': 'ET', 'FCEV': 'Canopy evaporation (FCEV)', 'FCTR': 'Transpiration (FCTR)', 'FGEV': 'Ground evaporation (FGEV)'}


def fgrid(vec):
    """forest-area-weighted mean per cell (NaN where no forest)."""
    g = grid_sum(vec, m) / np.where(mask, ff, np.nan)
    return g[:, order]


def setup(ax, gedi=False):
    ax.set_extent([-180, 180, -58, 84], crs=ccrs.PlateCarree())
    ax.coastlines(linewidth=0.3, color='0.35')
    if gedi:    # GEDI latitudinal limit, only on panels that involve the GEDI-based inputs
        ax.plot([-180, 180], [51.6, 51.6], color='0.3', lw=0.4, ls=(0, (3, 2)), transform=ccrs.PlateCarree())
        ax.plot([-180, 180], [-51.6, -51.6], color='0.3', lw=0.4, ls=(0, (3, 2)), transform=ccrs.PlateCarree())
    for s in ax.spines.values():
        s.set_linewidth(0.4)


def pmesh(ax, data, cmap, norm):
    return ax.pcolormesh(lone_s, late, np.ma.masked_invalid(data), cmap=cmap, norm=norm,
                         transform=ccrs.PlateCarree(), rasterized=True, shading='flat')


def save(fig, name):
    fig.savefig(OUT + name + '.png', dpi=300, bbox_inches='tight')
    fig.savefig(OUT + name + '.tif', dpi=300, bbox_inches='tight', pil_kwargs={'compression': 'tiff_lzw'})
    plt.close(fig)


# ---------------------------------------------------------------- Fig. 3 (heights on the model grid)
def fig3():
    H = {s: patch_height(s, m) for s in ['CLM Default', 'GEDI Max', 'GEDI Mean', 'GEDI Median']}
    Hg = {s: fgrid(np.nan_to_num(H[s])) for s in H}
    W, HT = 7.2, 6.22
    fig = plt.figure(figsize=(W, HT))

    def ax_in(x, y, w, h, **kw):    # axes placed in inches from the lower-left corner
        return fig.add_axes([x / W, y / HT, w / W, h / HT], **kw)

    mw = (W - 0.10 - 0.10 - 0.15) / 2; mh = mw * 142 / 360     # map extent is 360 x 142 degrees
    xl, xr = 0.10, 0.10 + mw + 0.15
    y1 = HT - 0.22 - mh; y2 = y1 - 0.26 - mh
    pos = [(xl, y1), (xr, y1), (xl, y2), (xr, y2)]
    hl = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40])
    hc = ListedColormap(plt.get_cmap('Greens')(np.linspace(0.15, 1, len(hl) - 1))); hn = BoundaryNorm(hl, hc.N)
    panels = [('CLM Default', None), ('GEDI Max', 'CLM Default'), ('GEDI Mean', 'CLM Default'), ('GEDI Median', 'CLM Default')]
    levels = np.array([-40, -30, -20, -10, -5, 5, 10, 20, 30, 40])
    b2 = plt.get_cmap('RdBu_r')(np.linspace(0, 1, len(levels) + 1)); b2[len(levels) // 2] = (0.93, 0.93, 0.93, 1)
    dc = ListedColormap(b2[1:-1]); dc.set_under(b2[0]); dc.set_over(b2[-1]); dn = BoundaryNorm(levels, dc.N)
    letters = 'abcd'
    for k, (s, ref) in enumerate(panels):
        ax = ax_in(*pos[k], mw, mh, projection=ccrs.PlateCarree()); setup(ax, gedi=ref is not None)
        if ref is None:
            im0 = pmesh(ax, Hg[s], hc, hn); ax.set_title(f'({letters[k]}) CLM Default forest canopy height (m)', fontsize=7.5, loc='left', pad=3)
        else:
            with np.errstate(all='ignore'):
                im1 = pmesh(ax, 100 * (Hg[s] - Hg[ref]) / Hg[ref], dc, dn)
            ax.set_title(f'({letters[k]}) {s} vs CLM Default (%)', fontsize=7.5, loc='left', pad=3)
    # (a) colour bar inset over the empty South Pacific (lon ~ -169 to -97, lat ~ -34 to -28)
    c0 = ax_in(xl + 0.03 * mw, y1 + 0.17 * mh, 0.20 * mw, 0.045 * mh)
    cb0 = fig.colorbar(im0, cax=c0, orientation='horizontal', extend='max', ticks=[0, 10, 20, 30, 40])
    cb0.ax.tick_params(labelsize=5.5, length=1.5, width=0.4, pad=1); cb0.outline.set_linewidth(0.4); cb0.ax.minorticks_off()
    # shared colour bar for (b-d), centred below the maps
    c1 = ax_in((W - 3.4) / 2, y2 - 0.20, 3.4, 0.08)
    cb1 = fig.colorbar(im1, cax=c1, orientation='horizontal', extend='both', ticks=levels)
    cb1.set_label('(b–d) Relative difference in forest-area-weighted canopy height (%)', fontsize=7, labelpad=2)
    cb1.ax.tick_params(labelsize=6, length=2, width=0.4, pad=1); cb1.outline.set_linewidth(0.4)
    # (e) PFT distributions (area-weighted quantiles)
    ax = ax_in(0.55, 0.30, W - 0.55 - 0.10, 2.0)
    w = m['wt'] * m['area'][m['jx'], m['ix']]
    cols = {'CLM Default': '#8c8c8c', 'GEDI Max': '#2a78d6', 'GEDI Mean': '#eb6834', 'GEDI Median': '#1baf7a'}
    for i, (pn, pid) in enumerate(PFTS.items()):
        sel = f & (m['iv'] == pid)
        for j, s in enumerate(cols):
            x = H[s][sel]; ww = w[sel]
            o = np.argsort(x); cw = np.cumsum(ww[o]) / ww.sum()
            q = [x[o][np.searchsorted(cw, p)] for p in (0.1, 0.25, 0.5, 0.75, 0.9)]
            xc = i + (j - 1.5) * 0.19
            ax.plot([xc, xc], [q[0], q[4]], color=cols[s], lw=0.8)
            ax.add_patch(plt.Rectangle((xc - 0.07, q[1]), 0.14, q[3] - q[1], facecolor=cols[s], edgecolor='white', lw=0.5))
            ax.plot([xc - 0.07, xc + 0.07], [q[2], q[2]], color='white', lw=1.0)
            ax.plot(xc, np.sum(x * ww) / ww.sum(), 'o', ms=2.6, mfc='white', mec='k', mew=0.5)
    ax.set_xticks(range(8)); ax.set_xticklabels(list(PFTS), fontsize=7); ax.set_ylabel('Canopy height (m)', fontsize=7.5)
    ax.set_xlim(-0.6, 7.6); ax.set_ylim(0, 42); ax.set_yticks(range(0, 41, 10)); ax.tick_params(labelsize=7, length=2.5, width=0.6)
    ax.spines[['top', 'right']].set_visible(False); ax.grid(axis='y', lw=0.3, color='0.85'); ax.set_axisbelow(True)
    ax.set_title('(e) Patch-level canopy height inputs by PFT', fontsize=7.5, loc='left', pad=4)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in cols.values()]
    ax.legend(handles, list(cols), ncol=4, frameon=False, fontsize=7, loc='lower right', bbox_to_anchor=(1, 1.0),
              handlelength=1.2, handleheight=0.8, columnspacing=1.2, borderaxespad=0.1)
    save(fig, 'Figure3_canopy_height_inputs')


# ---------------------------------------------------------------- Fig. 5 (PFT bars, ratio of means)
def fig5():
    w = m['wt'] * m['area'][m['jx'], m['ix']]
    cols = {'GEDI Max': '#2a78d6', 'GEDI Mean': '#eb6834', 'GEDI Median': '#1baf7a'}
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 4.8))
    out = {}
    for k, v in enumerate(['FCEV', 'FGEV', 'FCTR', 'ET']):
        ax = axes[k // 2, k % 2]
        for i, (pn, pid) in enumerate(PFTS.items()):
            sel = f & (m['iv'] == pid)
            for j, s in enumerate(cols):
                A = D[s][FLX[v]][sel]; B = D['CLM Default'][FLX[v]][sel]
                r = 100 * (np.sum(A * w[sel]) - np.sum(B * w[sel])) / np.sum(B * w[sel])
                ok = np.abs(B) > 1.0
                cell = 100 * (A[ok] - B[ok]) / np.abs(B[ok])
                o = np.argsort(cell); cw = np.cumsum(w[sel][ok][o]) / w[sel][ok].sum()
                q25, q75 = cell[o][np.searchsorted(cw, 0.25)], cell[o][np.searchsorted(cw, 0.75)]
                xc = i + (j - 1) * 0.26
                ax.bar(xc, r, 0.24, color=cols[s], edgecolor='white', lw=0.4, label=s if i == 0 else None)
                ax.plot([xc, xc], [q25, q75], color='0.15', lw=0.6)
                out[f'{pn}|{s}|{v}'] = dict(rel=float(r), q25=float(q25), q75=float(q75))
        ax.axhline(0, color='0.2', lw=0.5)
        ax.set_xticks(range(8)); ax.set_xticklabels(list(PFTS), rotation=40, ha='right', fontsize=7)
        ax.set_ylabel('Relative difference vs CLM Default (%)', fontsize=7)
        ax.set_title(f"({'abcd'[k]}) {LBL[v]}", fontsize=8, loc='left')
        ax.spines[['top', 'right']].set_visible(False); ax.grid(axis='y', lw=0.3, color='0.88')
        if k == 0: ax.legend(frameon=False, fontsize=7, ncol=3, loc='upper left')
    fig.tight_layout()
    save(fig, 'Figure5_PFT_bars')
    json.dump(out, open(OUT + 'Figure5_values.json', 'w'), indent=1)


if __name__ == '__main__':
    import sys
    which = sys.argv[1:] or ['3', '5']
    for w_ in which:
        {'3': fig3, '5': fig5}[w_]()
        print('done', w_, flush=True)
