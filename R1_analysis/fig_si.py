"""New SI figures for AFM R1: canopy water balance & seasonality, forcing-year robustness,
emulator validation, stomatal (GPP vs transpiration) diagnostic."""
import os, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import *

OUT = FIG
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.linewidth': 0.6})
m = load_meta(); f = m['forest']
w = m['wt'] * m['area'][m['jx'], m['ix']]
CAT = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948']


def save(fig, name):
    fig.savefig(OUT + name + '.png', dpi=300, bbox_inches='tight')
    fig.savefig(OUT + name + '.tif', dpi=300, bbox_inches='tight', pil_kwargs={'compression': 'tiff_lzw'})
    plt.close(fig)


def style(ax):
    ax.spines[['top', 'right']].set_visible(False); ax.grid(lw=0.3, color='0.9')


# ------------------------------------------------ canopy water balance + seasonality
def fig_wb():
    A = {}
    for s in ['CLM Default', 'Default 1.2', 'GEDI Max', 'GEDI Mean']:
        z = np.load(R1 + SCEN[s][2] + '_annual.npz')
        A[s] = {v: np.nanmean(z[v][30:40].astype(float), 0) for v in ['QVEGE', 'QFLX_LIQ_GRND', 'QFLX_SNOW_GRND']}
    Mo = {}
    for s in ['CLM Default', 'Default 1.2']:
        z = np.load(R1 + SCEN[s][2] + '_monclim.npz')
        Mo[s] = {v: z[v].astype(float) for v in ['QVEGE', 'QVEGT', 'QFLX_EVAP_TOT']}
    ci = m['ci']; ncol = len(m['cols1d_ixy'])
    has_f = np.zeros(ncol, bool); has_f[ci[f]] = True
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.5))
    ax = axes[0]
    for k, (a, b, c) in enumerate([('Default 1.2', 'CLM Default', CAT[0]), ('GEDI Max', 'GEDI Mean', CAT[1])]):
        dE = np.zeros(ncol); np.add.at(dE, ci, (A[a]['QVEGE'] - A[b]['QVEGE']) * m['wtcol'])
        dG = (A[a]['QFLX_LIQ_GRND'] + A[a]['QFLX_SNOW_GRND'] - A[b]['QFLX_LIQ_GRND'] - A[b]['QFLX_SNOW_GRND'])
        sel = has_f & np.isfinite(dG)
        ax.scatter(dE[sel] * SEC_YR, -dG[sel] * SEC_YR, s=2, color=c, alpha=0.5, lw=0, label=f'{a} vs {b}')
    lim = 25; ax.plot([-lim, lim], [-lim, lim], color='0.3', lw=0.6, ls='--')
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel('Change in canopy evaporation\n(column mean, mm yr$^{-1}$)'); ax.set_ylabel('Reduction in water reaching ground\n(throughfall + drip + unloading, mm yr$^{-1}$)')
    ax.set_title('(a) Canopy water balance closure', loc='left', fontsize=8); ax.legend(frameon=False, fontsize=6, markerscale=3, loc='upper left'); style(ax)
    # seasonal cycle by region (per unit forest area)
    plat = m['pfts1d_lat']
    regs = {'Tropics': np.abs(plat) < 23.5, 'Temperate': (np.abs(plat) >= 23.5) & (np.abs(plat) < 50), 'Boreal (NH ≥ 50°N)': plat >= 50}
    for j, (v, ttl) in enumerate([('QVEGE', '(b) Δ canopy evaporation, +20% height'), ('QVEGT', '(c) Δ transpiration, +20% height')]):
        ax = axes[j + 1]
        for k, (rn, rs) in enumerate(regs.items()):
            sel = f & rs
            d = [(np.sum((Mo['Default 1.2'][v][mo] - Mo['CLM Default'][v][mo])[sel] * w[sel]) / np.sum(w[sel])) * 86400 * 365 / 12 for mo in range(12)]
            ax.plot(range(1, 13), d, '-o', ms=2.5, lw=1.2, color=CAT[k], label=rn)
        ax.axhline(0, color='0.3', lw=0.5); ax.set_xticks(range(1, 13)); ax.set_xticklabels('JFMAMJJASOND')
        ax.set_ylabel('mm month$^{-1}$ (per unit forest area)'); ax.set_title(ttl, loc='left', fontsize=8); style(ax)
        if j == 0: ax.legend(frameon=False, fontsize=6)
    fig.tight_layout(); save(fig, 'FigureS_canopy_water_balance_seasonal')


# ------------------------------------------------ forcing-year robustness
def fig_years():
    r = json.load(open('an1_results.json'))['forcing_year_robustness']
    keys = list(r)
    fig, axes = plt.subplots(1, 4, figsize=(7.2, 2.6), sharey=True)
    for j, v in enumerate(['FCEV', 'FCTR', 'FGEV', 'ET']):
        ax = axes[j]
        for i, k in enumerate(keys):
            x = r[k][v]
            ax.plot([x['min'], x['max']], [i, i], color='0.3', lw=1.6, solid_capstyle='butt')
            ax.plot(x['mean'], i, 'o', ms=3, color=CAT[0])
            ax.plot(x['cycle3'], i, 'x', ms=3, color=CAT[1])
        ax.axvline(0, color='0.3', lw=0.5); ax.set_title(v, fontsize=8); style(ax)
        ax.set_xlabel('Relative difference (%)')
    axes[0].set_yticks(range(len(keys))); axes[0].set_yticklabels(keys, fontsize=6.5); axes[0].invert_yaxis()
    axes[-1].plot([], [], 'o', color=CAT[0], ms=3, label='10-yr mean (cycle 4)'); axes[-1].plot([], [], 'x', color=CAT[1], ms=3, label='10-yr mean (cycle 3)')
    axes[-1].plot([], [], color='0.3', lw=1.6, label='range of single years')
    axes[-1].legend(frameon=False, fontsize=6, loc='lower right')
    fig.tight_layout(); save(fig, 'FigureS_forcing_year_robustness')


# ------------------------------------------------ emulator validation + stomatal diagnostic
def fig_emu():
    D = {s: load_scen(s) for s in SCEN}
    emu = np.load('emulator_coef.npz')
    H = {s: patch_height(s, m) for s in ['CLM Default', 'GEDI Max', 'GEDI Mean', 'GEDI Median']}
    ff = forest_frac(m); msk = ff > 0
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.5))
    for j, v in enumerate(['FCTR', 'FCEV']):
        ax = axes[j]
        for k, s in enumerate(['GEDI Max', 'GEDI Mean', 'GEDI Median']):
            lr = np.log(np.where(f, H[s] / H['CLM Default'], 1.0)); lrc = np.clip(lr, np.log(0.6), np.log(1.6))
            c = emu[v]; yh = c[0] + c[1] * lrc + c[2] * lrc ** 2 + c[1] * (lr - lrc)
            yh = np.where(f, yh, D['CLM Default'][v])
            sim = (grid_sum(D[s][v], m) - grid_sum(D['CLM Default'][v], m))[msk] / ff[msk]
            em = (grid_sum(yh, m) - grid_sum(D['CLM Default'][v], m))[msk] / ff[msk]
            ax.scatter(em, sim, s=2, lw=0, alpha=0.5, color=CAT[k], label=s)
        lim = np.nanpercentile(np.abs(sim), 99.5) * 1.1
        ax.plot([-lim, lim], [-lim, lim], 'k--', lw=0.5); ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        ax.set_xlabel(f'Emulated Δ{v} (W m$^{{-2}}$)'); ax.set_ylabel(f'Simulated Δ{v} (W m$^{{-2}}$)')
        ax.set_title(f"({'ab'[j]}) {v}: GEDI scenario − CLM Default", loc='left', fontsize=8); style(ax)
        if j == 0: ax.legend(frameon=False, fontsize=6, markerscale=3)
    ax = axes[2]
    r = json.load(open('an6_phs_diag.json'))
    names = list(PFTS)
    x = np.arange(len(names))
    ax.bar(x - 0.18, [r[n]['dT_p20'] for n in names], 0.34, color=CAT[0], label='Transpiration')
    ax.bar(x + 0.18, [r[n]['dGPP_p20'] for n in names], 0.34, color=CAT[2], label='GPP')
    ax.axhline(0, color='0.3', lw=0.5); ax.set_xticks(x); ax.set_xticklabels(names, rotation=45, ha='right', fontsize=6.5)
    ax.set_ylabel('Default 1.2 vs CLM Default (%)'); ax.set_title('(c) Transpiration and GPP response', loc='left', fontsize=8)
    ax.legend(frameon=False, fontsize=6, loc='lower right', bbox_to_anchor=(1.0, 0.0)); style(ax)
    fig.tight_layout(); save(fig, 'FigureS_emulator_and_stomatal')


if __name__ == '__main__':
    fig_wb(); print('wb'); fig_years(); print('years'); fig_emu(); print('emu')
