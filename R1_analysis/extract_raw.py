"""Extract PFT-level annual means (all 40 yr) and the monthly climatology of the
analysis decade (2031-2040) from the zipped raw CLM5 h0 output of each case.
Output: <ETCTH_WORK>/<case>_annual.npz, <case>_monclim.npz, meta.npz
"""
import sys, os, zipfile
import numpy as np
import netCDF4 as nc
from multiprocessing import Pool

SRC = os.environ.get('ETCTH_RAW', 'G:/Hangkai/CTH_ET/')          # zipped raw h0/h1 output of each case
OUT = os.environ.get('ETCTH_WORK', 'G:/Hangkai/CTH_ET/R1_analysis/')
CASES = ['Default_case', 'Default_0.8_case', 'Default_0.9_case', 'Default_1.1_case',
         'Default_1.2_case', 'Max_case', 'Mean_case', 'Median_case']
VARS = ['FCEV', 'FCTR', 'FGEV', 'QVEGE', 'QVEGT', 'QFLX_EVAP_TOT', 'QFLX_LIQ_GRND',
        'QFLX_SNOW_GRND', 'QH2OSFC', 'FSH', 'FSH_V', 'FSH_G', 'Rnet', 'TV', 'TG', 'TSA',
        'ELAI', 'ESAI', 'FPSN', 'SABV', 'SABG', 'FIRE', 'FLDS', 'FSDS', 'FSR',
        'EFLX_LH_TOT', 'EFLX_SOIL_GRND', 'Q2M', 'RH']
YEARS = list(range(2001, 2041))
ANALYSIS = list(range(2031, 2041))


def run(case):
    z = zipfile.ZipFile(SRC + case + '_h0_h1.zip')
    ann, mon = {}, {}
    for iy, yr in enumerate(YEARS):
        name = f'{case}.clm2.h0.{yr}-02-01-00000.nc'
        d = nc.Dataset('m.nc', memory=z.read(name))
        tb = d['time_bounds'][:]
        ndays = (tb[:, 1] - tb[:, 0]).astype(np.float64)
        w = ndays / ndays.sum()
        for v in VARS:
            x = np.ma.filled(d[v][:].astype(np.float64), np.nan)   # (12, n) n = pft/column/gridcell
            if v not in ann:
                ann[v] = np.full((len(YEARS), x.shape[1]), np.nan, dtype=np.float32)
                mon[v] = np.zeros((12, x.shape[1]), dtype=np.float64)
            ann[v][iy] = np.nansum(x * w[:, None], axis=0) if np.isfinite(x).any() else np.nan
            # keep NaN where all months NaN
            allnan = np.all(~np.isfinite(x), axis=0)
            ann[v][iy][allnan] = np.nan
            if yr in ANALYSIS:
                mon[v] += x / len(ANALYSIS)
        if iy == 0 and case == 'Default_case':
            meta = {k: d[k][:] for k in ['pfts1d_itype_veg', 'pfts1d_ixy', 'pfts1d_jxy',
                                         'pfts1d_wtgcell', 'pfts1d_wtcol', 'pfts1d_lat', 'pfts1d_lon',
                                         'pfts1d_ci', 'cols1d_ixy', 'cols1d_jxy', 'cols1d_wtgcell',
                                         'cols1d_itype_lunit', 'area', 'landfrac', 'lat', 'lon']}
            np.savez(OUT + 'meta.npz', **{k: np.ma.filled(v, np.nan) for k, v in meta.items()})
        d.close()
        print(case, yr, flush=True)
    np.savez_compressed(OUT + case + '_annual.npz', years=np.array(YEARS), **ann)
    np.savez_compressed(OUT + case + '_monclim.npz', **{v: mon[v].astype(np.float32) for v in VARS})
    return case


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    with Pool(8) as p:
        for c in p.imap_unordered(run, CASES):
            print('DONE', c, flush=True)
