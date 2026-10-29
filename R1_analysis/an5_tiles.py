"""R1 analysis 5: sensitivity of aggregated canopy height to the intermediate scale, the treatment
of zero-height pixels and retrieval outliers, using native ETH 10-m tiles (3x3 deg, 1/12000 deg px).

Reference chain (as used for the GEDI scenarios): 10 m -> 15-arcsec (~500 m) blocks using
max / mean / median of valid pixels (0 < h < 255), then the mean over forest blocks
(the 3x3 deg tile stands in for a model grid cell).
Forest blocks: >= COVER of the valid pixels have h >= 5 m (proxy for the MODIS forest-PFT cells).
Variants:
  scale   : intermediate block of 10, 25, 50 (reference), 100 px  (~100 m, 250 m, 500 m, 1 km)
  direct  : no intermediate step (mean / median of all valid pixels in forest 500-m blocks)
  zeros   : h == 0 pixels counted as height 0 (mean / median)
  outliers: pixels > 60 m excluded; robust maxima P99 and P95 instead of max
Usage: python -I an5_tiles.py <tile.tif> <out.json>
"""
import sys, json
import numpy as np
import rasterio
from rasterio.windows import Window

path, outp = sys.argv[1], sys.argv[2]
COVERS = [0.3, 0.6]
STRIP = 100


def block_view(a, b):
    r, c = a.shape
    return a.reshape(r // b, b, c // b, b).transpose(0, 2, 1, 3).reshape(r // b, c // b, b * b)


def sorted_stats(blocks, keep, qs=()):
    """blocks: uint8 (..., n); keep: bool mask of pixels to use. Returns n, max, mean, median, quantiles
    (numpy 'linear' definition), NaN where n == 0."""
    x = np.where(keep, blocks, 255).astype(np.uint8)
    n = keep.sum(-1)
    s = np.sort(x, axis=-1)
    out = {'n': n}
    with np.errstate(all='ignore'):
        tot = np.where(keep, blocks.astype(np.float64), 0).sum(-1)
        out['mean'] = np.where(n > 0, tot / np.maximum(n, 1), np.nan)
        idx_max = np.clip(n - 1, 0, None)
        out['max'] = np.where(n > 0, np.take_along_axis(s, idx_max[..., None], -1)[..., 0].astype(float), np.nan)
        for name, q in [('median', 0.5)] + list(qs):
            pos = q * (n - 1)
            lo = np.floor(pos).astype(int); hi = np.ceil(pos).astype(int)
            lo = np.clip(lo, 0, None); hi = np.clip(hi, 0, None)
            vlo = np.take_along_axis(s, lo[..., None], -1)[..., 0].astype(float)
            vhi = np.take_along_axis(s, hi[..., None], -1)[..., 0].astype(float)
            out[name] = np.where(n > 0, vlo + (pos - lo) * (vhi - vlo), np.nan)
    return out


acc = {}
hist = {cv: {'valid': np.zeros(256, np.int64), 'incl0': np.zeros(256, np.int64)} for cv in COVERS}
counts = {'valid': 0, 'zero': 0, 'nodata': 0, 'gt60': 0, 'gt90': 0}


def add(cv, key, vals):
    vals = np.asarray(vals, float); vals = vals[np.isfinite(vals)]
    s = acc.setdefault((cv, key), [0.0, 0])
    s[0] += float(vals.sum()); s[1] += int(vals.size)


with rasterio.open(path) as src:
    H, W = src.height, src.width
    for r0 in range(0, H, STRIP):
        a = src.read(1, window=Window(0, r0, W, STRIP))
        counts['nodata'] += int((a == 255).sum()); counts['zero'] += int((a == 0).sum())
        valid = (a > 0) & (a < 255)
        counts['valid'] += int(valid.sum()); counts['gt60'] += int((valid & (a > 60)).sum())
        counts['gt90'] += int((valid & (a > 90)).sum())
        # ---- 500 m blocks
        b50 = block_view(a, 50); v50 = block_view(valid, 50)
        nv = v50.sum(-1)
        tall = (block_view((a >= 5) & valid, 50)).sum(-1)
        cov50 = np.where(nv > 0, tall / np.maximum(nv, 1), 0)
        st = sorted_stats(b50, v50, qs=[('p99', 0.99), ('p95', 0.95)])
        stz = sorted_stats(b50, b50 < 255)
        st6 = sorted_stats(b50, v50 & (b50 <= 60))
        nz = (b50 == 0).sum(-1)
        res50 = {'max': st['max'], 'mean': st['mean'], 'median': st['median'], 'p99': st['p99'], 'p95': st['p95'],
                 'mean_incl0': stz['mean'], 'median_incl0': stz['median'],
                 'max_le60': st6['max'], 'mean_le60': st6['mean'], 'median_le60': st6['median'],
                 'zero_frac_in_forest': nz / np.maximum(nz + nv, 1)}
        for cv in COVERS:
            fm = cov50 >= cv
            for k, v in res50.items():
                add(cv, f'500m_{k}', v[fm])
            add(cv, 'forest_blocks_500m', np.ones(int(fm.sum())))
            fpx = np.repeat(np.repeat(fm, 50, 0), 50, 1)
            hist[cv]['valid'] += np.bincount(a[fpx & valid], minlength=256)
            hist[cv]['incl0'] += np.bincount(a[fpx & (a < 255)], minlength=256)
        # ---- other intermediate scales (forest mask defined at that scale)
        for b in [10, 25, 100]:
            bb = block_view(a, b); vb = block_view(valid, b)
            nvb = vb.sum(-1)
            covb = np.where(nvb > 0, block_view((a >= 5) & valid, b).sum(-1) / np.maximum(nvb, 1), 0)
            sb = sorted_stats(bb, vb)
            for cv in COVERS:
                fm = covb >= cv
                for k in ['max', 'mean', 'median']:
                    add(cv, f'{b * 10}m_{k}', sb[k][fm])


def hist_stats(hh):
    vals = np.arange(256); n = hh.sum()
    if n == 0:
        return np.nan, np.nan
    cdf = np.cumsum(hh) / n
    return float((hh * vals).sum() / n), float(vals[np.searchsorted(cdf, 0.5)])


out = {'tile': path.replace('\\', '/').split('/')[-1], 'counts': counts}
for cv in COVERS:
    d = {k: acc[(c, k)][0] / max(acc[(c, k)][1], 1) for (c, k) in acc if c == cv}
    d['n_forest_blocks_500m'] = acc[(cv, 'forest_blocks_500m')][1]
    d['direct_mean'], d['direct_median'] = hist_stats(hist[cv]['valid'])
    d['direct_mean_incl0'], d['direct_median_incl0'] = hist_stats(hist[cv]['incl0'])
    out[f'cover{cv}'] = d
json.dump(out, open(outp, 'w'), indent=1)
print(json.dumps(out)[:2000])
