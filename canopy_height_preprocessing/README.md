# Canopy Height Preprocessing

Reference scripts used to aggregate the 10 m ETH Global Canopy Height
(Lang et al., 2023) to the CLM5 1.9 x 2.5 degree grid per PFT.

## Workflow

1. **Stage A (10 m -> 500 m, external):** The 10 m canopy height is
   aggregated to 500 m using three statistics: `max`, `mean`, and `median`,
   keeping only valid pixels (height > 0 m). This produces three 500 m
   GeoTIFFs used as the three GEDI scenarios (Max, Mean, Median).

2. **Stage B (500 m -> CLM5 grid):** The 500 m canopy height is further
   aggregated to the CLM5 1.9 x 2.5 degree grid using simple area-weighted
   mean within each grid cell, stratified by PFT (1-8 forest types). The
   500 m PFT map is derived from MODIS via
   [build_surface_dataset_for_ELM](https://github.com/daleihao/build_surface_dataset_for_ELM)
   (`ELM_PFT_{year}-WGS84-merged.tif`).

3. **Stage C (surface data update):** The aggregated heights replace
   `MONTHLY_HEIGHT_TOP` in the CLM5 surface dataset. `MONTHLY_HEIGHT_BOT`
   is derived from top height using fixed PFT-specific ratios:

   | PFT | Name | HBOT/HTOP |
   |---|---|---|
   | 1 | NET Temperate | 0.5000 |
   | 2 | NET Boreal | 0.5000 |
   | 3 | NDT Boreal | 0.5000 |
   | 4 | BET Tropical | 0.0286 |
   | 5 | BET Temperate | 0.0286 |
   | 6 | BDT Tropical | 0.5556 |
   | 7 | BDT Temperate | 0.5750 |
   | 8 | BDT Boreal | 0.5750 |

## Files

| File | Description |
|---|---|
| `input_data_generation_hangkai.py` | Python script: reads 500 m canopy height + 500 m PFT map, aggregates to CLM5 1.9 x 2.5 grid per PFT, writes modified surfdata NC. Also generates standalone PFT and canopy-height NC products. |
| `cht_with_pft.m` | MATLAB equivalent of the same aggregation logic. |

Both scripts implement the same Stage B + Stage C workflow; they are kept
together for transparency and reproducibility.

See also: [https://github.com/youhangkai/clm_upscale_with_pft](https://github.com/youhangkai/clm_upscale_with_pft)
