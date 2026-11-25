# Canopy Height Preprocessing

Scripts used to aggregate the 10 m ETH Global Canopy Height (Lang et al.,
2023) to the CLM5 1.9 x 2.5 degree grid per PFT and to write the heights into
the CLM5 surface datasets.

## Workflow

1. **Stage A (10 m -> 500 m, external):** The 10 m canopy height is
   aggregated to 500 m using three statistics: `max`, `mean`, and `median`,
   keeping only valid pixels (height > 0 m). This produces three 500 m
   GeoTIFFs used as the three GEDI scenarios (Max, Mean, Median). The 500 m
   height maps and the 500 m PFT map are on a regular latitude-longitude grid
   of 80,150 x 40,076 pixels (about 0.0045 degrees, or 500 m at the equator).
   The PFT map is derived from MODIS via
   [build_surface_dataset_for_ELM](https://github.com/daleihao/build_surface_dataset_for_ELM)
   (`ELM_PFT_{year}-WGS84-merged.tif`).

2. **Stage B (500 m -> model grid, `cth_with_pft.m`):** For each forest PFT
   (1-8), the 500 m heights of all pixels of that PFT are averaged with equal
   weight within 144 x 96 equal blocks of 2.5 x 1.875 degrees (about
   557 x 417 pixels; one block for each model grid cell), counted from
   180 degrees W and from the north. The block means are stored as whole
   meters (`uint8`, rounded); 1 m is the lowest stored value and occurs in a
   few patches covering less than 0.04% of the forest area.

3. **Stage C (surface data update, `CTH_for_surface_data_replace.m`):** The
   block heights are reordered from 180 W-180 E to the 0-360 degree order of
   the CLM5 surface dataset and replace `MONTHLY_HEIGHT_TOP` of the forest
   PFTs in all 12 months wherever both the block height and the CLM5 default
   height are positive. Where a forest PFT of the CLM5 surface dataset is
   absent from the 500 m PFT map within the block, the CLM5 default height is
   kept. `MONTHLY_HEIGHT_BOT` is derived from top height using fixed
   PFT-specific ratios:

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

## Position of the blocks relative to the model grid

The blocks start at multiples of 2.5 degrees longitude, whereas the model
cells are centered on multiples of 2.5 degrees (0, 2.5, ... 357.5), so no
whole-column shift aligns the two grids. `CTH_for_surface_data_replace.m`
(`line = j-71`) pairs each block with the model cell centered on the eastern
edge of the block: the center of the block is 1.25 degrees (half a cell) west
of the center of the cell. In latitude, block rows (1.875 degrees) and model
rows (1.895 degrees apart) are paired by index, so the block center is closer
to the Equator than the cell center by up to 0.94 degrees at the poles (0.62
degrees at 60 degrees latitude, 0 at the Equator). The offset is the same in
the three GEDI-based inputs. The uniform-scaling simulations use the CLM5
default heights and are not affected.

## Files

| File | Description |
|---|---|
| `cth_with_pft.m` | MATLAB, used for the simulations (Stage B): 500 m canopy height + 500 m PFT map -> block means per PFT. Run once per scenario with the paths changed; the file is the last saved version (Median). |
| `CTH_for_surface_data_replace.m` | MATLAB, used for the simulations (Stage C): writes the block heights into the CLM5 surface dataset. Run once per scenario with the paths changed; the file is the last saved version (Max). |
| `input_data_generation_hangkai.py` | Earlier Python version of Stages B and C, not used for the simulations. It assumes a 1/240 degree grid (blocks of 450 x 600 pixels), writes 1 m where a PFT is absent, and reorders the longitudes by 72 columns. |
| `cht_with_pft.m` | Earlier version of `cth_with_pft.m`, kept for reference. |
