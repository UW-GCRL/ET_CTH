# Canopy Height Preprocessing (MATLAB)

MATLAB scripts for aggregating 10 m ETH Global Canopy Height data to CLM5 grid resolution (1.9x2.5 degrees) per PFT.

## Height Thresholds

In all aggregation methods (max, mean, median), only pixels with:
- **Canopy height > 0 m** (excludes non-vegetated areas and no-data pixels)
- **Canopy height < 90 m** (excludes unrealistic outliers)

are included. See `cth_with_pft.m`, line 27: `if((CTH(i,j)>0&&(CTH(i,j)<90)))`

## Scripts

| Script | Description |
|---|---|
| `cth_with_pft.m` | Aggregates 10 m canopy height to 2.5x1.9 degree grid per PFT (8 forest PFTs). Applies height thresholds (>0, <90 m). Computes weighted mean height per grid cell. |
| `upscale_cth.m` | Upscales canopy height with PFT masking at 500 m resolution. |
| `CTH_for_surface_data_replace.m` | Replaces `MONTHLY_HEIGHT_TOP` in CLM5 surface data file with GEDI-derived canopy heights per PFT. |
| `CTHcompare.m` | Compares canopy height across scenarios. |
| `canopyheightcompare.m` | Additional canopy height comparison utilities. |

## Data Sources

- **ETH Global Canopy Height**: 10 m resolution, 2020 version (Lang et al., 2023)
- **PFT map**: 500 m resolution, derived from MODIS land cover (2020)
- **CLM5 surface data**: `surfdata_1.9x2.5_2005.nc` with `MONTHLY_HEIGHT_TOP` variable
