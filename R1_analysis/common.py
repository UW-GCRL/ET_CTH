"""Shared loaders for the R1 analyses (10-yr mean PFT-level CLM5 output)."""
import os
import numpy as np
import netCDF4 as nc

# Data locations (edit here or set the environment variables)
FINAL = os.environ.get('ETCTH_FINAL', 'G:/Hangkai/CTH_ET project/Final_data/')   # 10-yr mean PFT-level output, one file per scenario
SURF = os.environ.get('ETCTH_SURF', 'G:/Hangkai/CTH_ET/surfacedata/')            # surface datasets of the eight scenarios
R1 = os.environ.get('ETCTH_WORK', 'G:/Hangkai/CTH_ET/R1_analysis/')              # npz files written by extract_raw.py
FIG = os.environ.get('ETCTH_FIGURES', R1 + 'figures/')                           # figure output
SCEN = {  # display name: (Final_data file, surfdata file, raw case)
    'CLM Default': ('CLM Default.nc', 'Default.nc', 'Default_case'),
    'Default 0.8': ('Default 0.8.nc', 'Default_0.8.nc', 'Default_0.8_case'),
    'Default 0.9': ('Default 0.9.nc', 'Default_0.9.nc', 'Default_0.9_case'),
    'Default 1.1': ('Default 1.1.nc', 'Default_1.1.nc', 'Default_1.1_case'),
    'Default 1.2': ('Default 1.2.nc', 'Default_1.2.nc', 'Default_1.2_case'),
    'GEDI Mean': ('GEDI Mean.nc', 'Mean.nc', 'Mean_case'),
    'GEDI Max': ('GEDI Max.nc', 'Max.nc', 'Max_case'),
    'GEDI Median': ('GEDI Median.nc', 'Median.nc', 'Median_case'),
}
PFTS = {'NET Temp': 1, 'NET Bor': 2, 'NDT Bor': 3, 'BET Trop': 4, 'BET Temp': 5,
        'BDT Trop': 6, 'BDT Temp': 7, 'BDT Bor': 8}
HVAP = 2.501e6   # J/kg, CLM SHR_CONST_LATVAP
HSUB = 2.8345e6  # J/kg, CLM latvap+latice
SEC_YR = 365 * 86400.0  # CLM noleap calendar

VARS = ['FCEV', 'FCTR', 'FGEV', 'QVEGE', 'QVEGT', 'QFLX_EVAP_TOT', 'QFLX_LIQ_GRND',
        'QFLX_SNOW_GRND', 'FSH', 'Rnet', 'TV', 'ELAI', 'ESAI', 'FSH_V', 'FSH_G', 'TG',
        'QH2OSFC', 'FPSN', 'SABV', 'SABG', 'FIRE', 'EFLX_LH_TOT', 'TSA', 'RH', 'Q2M']


def load_meta():
    with nc.Dataset(FINAL + 'CLM Default.nc') as d:
        m = {k: np.asarray(np.ma.filled(d[k][:], np.nan)) for k in
             ['pfts1d_itype_veg', 'pfts1d_ixy', 'pfts1d_jxy', 'pfts1d_wtgcell', 'area',
              'landfrac', 'lat', 'lon', 'pfts1d_lat', 'pfts1d_lon', 'pfts1d_wtcol',
              'pfts1d_ci', 'cols1d_ixy', 'cols1d_jxy', 'cols1d_wtgcell', 'cols1d_itype_lunit']}
    m['iv'] = m['pfts1d_itype_veg'].astype(int)
    m['ix'] = m['pfts1d_ixy'].astype(int) - 1
    m['jx'] = m['pfts1d_jxy'].astype(int) - 1
    m['wt'] = m['pfts1d_wtgcell']
    # CLM 'area' is the whole grid cell (land + ocean) and pfts1d_wtgcell is relative to the land part of the
    # cell, so the land area of a patch is wtgcell * area * landfrac. From here on m['area'] is the land area of
    # the cell; the whole-cell area is kept as m['area_cell'].
    m['area_cell'] = m['area'].copy()
    m['area'] = np.where(np.isfinite(m['area']), m['area'], 0.0) * np.nan_to_num(m['landfrac'])
    m['forest'] = (m['iv'] >= 1) & (m['iv'] <= 8) & (m['wt'] > 0)
    m['ci'] = m['pfts1d_ci'].astype(int) - 1          # patch -> column index (0-based)
    m['wtcol'] = m['pfts1d_wtcol']
    return m


def load_scen(name, vars_=VARS):
    with nc.Dataset(FINAL + SCEN[name][0]) as d:
        out = {}
        for v in vars_:
            if v in d.variables:
                x = np.asarray(np.ma.filled(d[v][:].astype(np.float64), np.nan))
                out[v] = x.mean(axis=0) if x.ndim == 2 else x
    # water fluxes in mm/yr; ground evaporation as residual of total minus vegetation
    out['E_can'] = out['QVEGE'] * SEC_YR
    out['T'] = out['QVEGT'] * SEC_YR
    out['E_gnd'] = (out['QFLX_EVAP_TOT'] - out['QVEGE'] - out['QVEGT']) * SEC_YR
    out['ET_w'] = out['QFLX_EVAP_TOT'] * SEC_YR
    out['ET'] = out['FCEV'] + out['FCTR'] + out['FGEV']
    # column-level (natural-vegetation soil column) water reaching the ground after interception
    out['LIQ_GRND_c'] = out.pop('QFLX_LIQ_GRND') * SEC_YR
    out['SNOW_GRND_c'] = out.pop('QFLX_SNOW_GRND') * SEC_YR
    out['GRND_IN_c'] = out['LIQ_GRND_c'] + out['SNOW_GRND_c']
    for v in ['QH2OSFC', 'TG']:
        if v in out:
            out[v + '_c'] = out.pop(v)
    return out


def grid_sum(vec, m, sel=None):
    """Sum of flux*wtgcell over selected patches per grid cell (per unit grid-cell area)."""
    sel = m['forest'] if sel is None else sel
    g = np.zeros(m['area'].shape)
    np.add.at(g, (m['jx'][sel], m['ix'][sel]), (vec * m['wt'])[sel])
    return g


def forest_frac(m, sel=None):
    sel = m['forest'] if sel is None else sel
    g = np.zeros(m['area'].shape)
    np.add.at(g, (m['jx'][sel], m['ix'][sel]), m['wt'][sel])
    return g


def gmean(grid, m, mask=None):
    """Area-weighted mean across grid cells (weights = land area of the cell, m['area'])."""
    a = np.where(np.isfinite(m['area']), m['area'], 0.0)
    if mask is not None:
        a = a * mask
    return np.nansum(grid * a) / np.nansum(a)


def load_heights(name, month=0):
    with nc.Dataset(SURF + SCEN[name][1]) as d:
        ht = np.asarray(d['MONTHLY_HEIGHT_TOP'][month])  # (lsmpft, lat, lon)
        hb = np.asarray(d['MONTHLY_HEIGHT_BOT'][month])
        pct = np.asarray(d['PCT_NAT_PFT'][:])
    return ht, hb, pct


def patch_height(name, m):
    ht, hb, _ = load_heights(name)
    h = np.full(m['iv'].shape, np.nan)
    sel = m['forest']
    h[sel] = ht[m['iv'][sel], m['jx'][sel], m['ix'][sel]]
    return h
