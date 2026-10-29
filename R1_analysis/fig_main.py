"""Revised main-text figures 2-5 and related SI figures for AFM R1.
Conventions (stated in captions):
 * fluxes are water fluxes from CLM5's native water-flux outputs (QVEGE, QVEGT, and
   QFLX_EVAP_TOT - QVEGE - QVEGT for ground evaporation), mm yr-1;
 * map values are per unit forest area within each grid cell (forest-PFT-weighted mean);
 * maps are drawn cell by cell (no interpolation) on the CLM5 1.9x2.5 grid.
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
SHORT = {'ET': 'ET', 'FCEV': 'FCEV', 'FCTR': 'FCTR', 'FGEV': 'FGEV'}


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


div_levels = np.array([-10, -8, -6, -4, -2, -0.5, 0.5, 2, 4, 6, 8, 10])
base = plt.get_cmap('RdBu_r')(np.linspace(0, 1, len(div_levels) + 1))
base[len(div_levels) // 2] = (0.93, 0.93, 0.93, 1)   # neutral grey midpoint
div_cmap = ListedColormap(base[1:-1]); div_cmap.set_under(base[0]); div_cmap.set_over(base[-1])
div_norm = BoundaryNorm(div_levels, div_cmap.N)

abs_levels = np.array([0, 50, 100, 200, 300, 400, 500, 600, 800, 1000, 1200, 1400])
ac = plt.get_cmap('Purples')(np.linspace(0.12, 1, len(abs_levels) - 1))
abs_cmap = ListedColormap(ac); abs_cmap.set_under('#b35806'); abs_cmap.set_over(ac[-1])
abs_norm = BoundaryNorm(abs_levels, abs_cmap.N)


def rel(a, b, v):
    A = fgrid(D[a][FLX[v]]); B = fgrid(D[b][FLX[v]])
    with np.errstate(all='ignore'):
        r = 100 * (A - B) / np.abs(B)
    r[np.abs(B) < 1.0] = np.nan     # undefined where the baseline flux is < 1 mm/yr in magnitude
    return r


def save(fig, name):
    fig.savefig(OUT + name + '.png', dpi=300, bbox_inches='tight')
    fig.savefig(OUT + name + '.tif', dpi=300, bbox_inches='tight', pil_kwargs={'compression': 'tiff_lzw'})
    plt.close(fig)


# ---------------------------------------------------------------- Fig. 2
def fig2():
    cols = ['Default 0.8', 'Default 0.9', 'CLM Default', 'Default 1.1', 'Default 1.2']
    fig, axes = plt.subplots(4, 5, figsize=(7.2, 3.35), subplot_kw=dict(projection=ccrs.PlateCarree()))
    for i, v in enumerate(['ET', 'FCEV', 'FCTR', 'FGEV']):
        for j, s in enumerate(cols):
            ax = axes[i, j]; setup(ax)
            if s == 'CLM Default':
                im_a = pmesh(ax, fgrid(D[s][FLX[v]]), abs_cmap, abs_norm)
            else:
                im_d = pmesh(ax, rel(s, 'CLM Default', v), div_cmap, div_norm)
            if i == 0:
                ax.set_title('CLM Default (mm yr$^{-1}$)' if s == 'CLM Default' else f'{s} vs CLM Default', fontsize=7)
            if j == 0:
                ax.text(-0.06, 0.5, SHORT[v], transform=ax.transAxes, rotation=90, va='center', ha='right', fontsize=8)
    cax1 = fig.add_axes([0.1, 0.075, 0.36, 0.022]); cax2 = fig.add_axes([0.58, 0.075, 0.36, 0.022])
    cb = fig.colorbar(im_d, cax=cax1, orientation='horizontal', extend='both', ticks=[-10, -6, -2, 2, 6, 10])
    cb.set_label('Relative difference from CLM Default (%)', fontsize=7); cb.ax.tick_params(labelsize=6)
    cb2 = fig.colorbar(im_a, cax=cax2, orientation='horizontal', extend='both', ticks=[0, 200, 400, 600, 1000, 1400])
    cb2.set_label('CLM Default forest flux (mm yr$^{-1}$; brown: net condensation < 0)', fontsize=7); cb2.ax.tick_params(labelsize=6)
    fig.subplots_adjust(left=0.04, right=0.995, top=0.94, bottom=0.17, wspace=0.03, hspace=0.06)
    save(fig, 'Figure2_sensitivity_maps')


# ---------------------------------------------------------------- Fig. 4
def fig4():
    combos = [('GEDI Mean', 'CLM Default'), ('GEDI Max', 'CLM Default'), ('GEDI Max', 'GEDI Mean'), ('GEDI Median', 'GEDI Mean')]
    fig, axes = plt.subplots(4, 4, figsize=(7.2, 4.0), subplot_kw=dict(projection=ccrs.PlateCarree()))
    for i, v in enumerate(['ET', 'FCEV', 'FCTR', 'FGEV']):
        for j, (a, b) in enumerate(combos):
            ax = axes[i, j]; setup(ax, gedi=True)
            im = pmesh(ax, rel(a, b, v), div_cmap, div_norm)
            if i == 0: ax.set_title(f'{a} vs {b}', fontsize=7)
            if j == 0: ax.text(-0.06, 0.5, SHORT[v], transform=ax.transAxes, rotation=90, va='center', ha='right', fontsize=8)
    cax = fig.add_axes([0.3, 0.07, 0.42, 0.02])
    cb = fig.colorbar(im, cax=cax, orientation='horizontal', extend='both', ticks=[-10, -6, -2, 2, 6, 10])
    cb.set_label('Relative difference (%)', fontsize=7); cb.ax.tick_params(labelsize=6)
    fig.subplots_adjust(left=0.05, right=0.99, top=0.95, bottom=0.14, wspace=0.03, hspace=0.06)
    save(fig, 'Figure4_GEDI_maps')


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


# ---------------------------------------------------------------- SI Fig. S4 (PFT heatmap, ratio of means)
def figS4():
    w = m['wt'] * m['area'][m['jx'], m['ix']]
    comps = [('Default 0.8', 'CLM Default'), ('Default 0.9', 'CLM Default'), ('Default 1.1', 'CLM Default'), ('Default 1.2', 'CLM Default'),
             ('GEDI Mean', 'CLM Default'), ('GEDI Max', 'CLM Default'), ('GEDI Median', 'CLM Default'),
             ('GEDI Max', 'GEDI Mean'), ('GEDI Median', 'GEDI Mean')]
    rlab = [f'{a} vs {b}' for a, b in comps]
    lv = np.array([-8, -6, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 6, 8])
    b = plt.get_cmap('RdBu_r')(np.linspace(0, 1, len(lv) + 1)); b[len(lv) // 2] = (0.93, 0.93, 0.93, 1)
    cm = ListedColormap(b[1:-1]); cm.set_under(b[0]); cm.set_over(b[-1]); nm = BoundaryNorm(lv, cm.N)
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.4))
    out = {}
    for k, v in enumerate(['FCEV', 'FGEV', 'FCTR', 'ET']):
        ax = axes[k // 2, k % 2]
        M = np.zeros((len(comps), len(PFTS)))
        for i, (a, bb) in enumerate(comps):
            for j, (pn, pid) in enumerate(PFTS.items()):
                sel = f & (m['iv'] == pid)
                A = np.sum((D[a][FLX[v]] * w)[sel]); B = np.sum((D[bb][FLX[v]] * w)[sel])
                M[i, j] = 100 * (A - B) / B
                out[f'{pn}|{a} vs {bb}|{v}'] = float(M[i, j])
        im = ax.imshow(M, cmap=cm, norm=nm, aspect='auto')
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                ax.text(j, i, '0.0' if abs(M[i, j]) < 0.05 else f'{M[i, j]:+.1f}'.replace('-', '−'), ha='center', va='center', fontsize=5.5,
                        color='white' if abs(M[i, j]) >= 4 else '0.15')
        ax.axhline(3.5, color='white', lw=2); ax.axhline(6.5, color='white', lw=2)
        ax.set_xticks(range(len(PFTS))); ax.set_xticklabels(list(PFTS), rotation=40, ha='right', fontsize=6.5)
        ax.set_yticks(range(len(comps))); ax.set_yticklabels(rlab if k % 2 == 0 else [], fontsize=6.5)
        ax.tick_params(length=0)
        for s_ in ax.spines.values():
            s_.set_visible(False)
        ax.set_title(f"({'abcd'[k]}) {LBL[v]}", fontsize=8, loc='left', pad=4)
    cax = fig.add_axes([0.3, 0.035, 0.42, 0.015])
    cb = fig.colorbar(im, cax=cax, orientation='horizontal', extend='both', ticks=lv)
    cb.set_label('Relative difference (%)', fontsize=7); cb.ax.tick_params(labelsize=6, length=2); cb.outline.set_linewidth(0.4)
    fig.subplots_adjust(left=0.2, right=0.99, top=0.96, bottom=0.15, wspace=0.04, hspace=0.34)
    save(fig, 'FigureS4_PFT_relative_heatmap')
    json.dump(out, open(OUT + 'FigureS4_values.json', 'w'), indent=1)


# ---------------------------------------------------------------- SI: absolute maps for all scenarios
def figS_abs():
    sc = ['CLM Default', 'GEDI Mean', 'GEDI Max', 'GEDI Median']
    fig, axes = plt.subplots(4, 4, figsize=(7.2, 4.0), subplot_kw=dict(projection=ccrs.PlateCarree()))
    for i, v in enumerate(['ET', 'FCEV', 'FCTR', 'FGEV']):
        for j, s in enumerate(sc):
            ax = axes[i, j]; setup(ax); im = pmesh(ax, fgrid(D[s][FLX[v]]), abs_cmap, abs_norm)
            if i == 0: ax.set_title(s, fontsize=7)
            if j == 0: ax.text(-0.06, 0.5, SHORT[v], transform=ax.transAxes, rotation=90, va='center', ha='right', fontsize=8)
    cax = fig.add_axes([0.3, 0.07, 0.42, 0.02])
    cb = fig.colorbar(im, cax=cax, orientation='horizontal', extend='both', ticks=[0, 200, 400, 600, 1000, 1400])
    cb.set_label('Forest flux (mm yr$^{-1}$; brown: net condensation < 0)', fontsize=7); cb.ax.tick_params(labelsize=6)
    fig.subplots_adjust(left=0.05, right=0.99, top=0.95, bottom=0.14, wspace=0.03, hspace=0.06)
    save(fig, 'FigureS_absolute_maps')


def figS_absdiff():
    combos = [('GEDI Mean', 'CLM Default'), ('GEDI Max', 'CLM Default'), ('GEDI Median', 'CLM Default'), ('GEDI Max', 'GEDI Mean')]
    lv = np.array([-60, -40, -20, -10, -5, -1, 1, 5, 10, 20, 40, 60])
    b = plt.get_cmap('BrBG')(np.linspace(0, 1, len(lv) + 1)); b[len(lv) // 2] = (0.93, 0.93, 0.93, 1)
    cm = ListedColormap(b[1:-1]); cm.set_under(b[0]); cm.set_over(b[-1]); nm = BoundaryNorm(lv, cm.N)
    fig, axes = plt.subplots(4, 4, figsize=(7.2, 4.0), subplot_kw=dict(projection=ccrs.PlateCarree()))
    for i, v in enumerate(['ET', 'FCEV', 'FCTR', 'FGEV']):
        for j, (a, bb) in enumerate(combos):
            ax = axes[i, j]; setup(ax)
            im = pmesh(ax, fgrid(D[a][FLX[v]]) - fgrid(D[bb][FLX[v]]), cm, nm)
            if i == 0: ax.set_title(f'{a} vs {bb}', fontsize=7)
            if j == 0: ax.text(-0.06, 0.5, SHORT[v], transform=ax.transAxes, rotation=90, va='center', ha='right', fontsize=8)
    cax = fig.add_axes([0.3, 0.07, 0.42, 0.02])
    cb = fig.colorbar(im, cax=cax, orientation='horizontal', extend='both', ticks=lv)
    cb.set_label('Absolute difference in forest flux (mm yr$^{-1}$)', fontsize=7); cb.ax.tick_params(labelsize=6)
    fig.subplots_adjust(left=0.05, right=0.99, top=0.95, bottom=0.14, wspace=0.03, hspace=0.06)
    save(fig, 'FigureS_absolute_difference_maps')


if __name__ == '__main__':
    import sys
    which = sys.argv[1:] or ['2', '3', '4', '5', 'abs', 'absdiff']
    for w_ in which:
        {'2': fig2, '3': fig3, '4': fig4, '5': fig5, 'S4': figS4, 'abs': figS_abs, 'absdiff': figS_absdiff}[w_]()
        print('done', w_, flush=True)
