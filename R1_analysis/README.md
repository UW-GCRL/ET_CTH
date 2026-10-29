# Analysis for the revised manuscript

These scripts produce all results, tables and figures of the revised manuscript from the CLM5 output of the eight
simulations (CLM Default, Default 0.8, 0.9, 1.1 and 1.2, GEDI Max, GEDI Mean and GEDI Median).

## Area weighting

Global and PFT-level means are weighted by the land area of each patch:

    patch land area = pfts1d_wtgcell x area x landfrac

In CLM5 history files, `area` is the area of the whole grid cell (land and ocean; km2), `landfrac` is the land
fraction of the cell, and `pfts1d_wtgcell` is the weight of the patch relative to the land part of the cell (the
weights of all patches of a cell sum to 1). The history files state that none of the variables is weighted by land
fraction. The weighting is set once, in `common.load_meta`: `m['area']` is the land area of each cell and the
whole-cell area is kept as `m['area_cell']`. Earlier versions of these scripts weighted by the whole-cell area,
which counted cells that are partly ocean as if they were all land; all numbers of the revised manuscript use the
land-area weights.

## Inputs

Set the paths at the top of `common.py`, or with environment variables:

| Variable | Contents |
|---|---|
| `ETCTH_FINAL` | One netCDF file per scenario with the 10-year mean (2031-2040) PFT-level output |
| `ETCTH_SURF` | The surface datasets of the eight scenarios (`MONTHLY_HEIGHT_TOP`, `MONTHLY_HEIGHT_BOT`) |
| `ETCTH_RAW` | Zipped raw h0/h1 output of each case (read by `extract_raw.py` only) |
| `ETCTH_WORK` | Annual means (2001-2040), monthly climatology (2031-2040) and metadata written by `extract_raw.py` |
| `ETCTH_FIGURES` | Output folder for the figures (default: `<ETCTH_WORK>/figures/`) |

`an5_tiles.py` reads native 10 m tiles of the ETH Global Canopy Height (Lang et al., 2023).

## Run order

Run from this folder (the analysis scripts write their results next to themselves):

    python extract_raw.py
    python an1_checks.py
    python an2_seasonal.py
    python an4_heights.py
    python an3_spatial.py
    python an6_stomatal.py
    python fig_main.py
    python fig_si.py
    python make_tables.py
    python text_numbers.py
    python an5_tiles.py <tile.tif> <out.json>      # once per 10 m tile

## Scripts

| Script | Content | Used for |
|---|---|---|
| `common.py` | Scenario names, loaders, patch weights, grid helpers | all scripts |
| `extract_raw.py` | PFT-level annual means and monthly climatology from the raw output | inputs of `an1`, `an2`, `fig_si` |
| `an1_checks.py` | Snow burial of the canopy, spin-up drift (third vs fourth forcing cycle), forcing-year robustness, canopy water balance, energy diagnostics | Section 3.3, Fig. S5, Text S3 |
| `an2_seasonal.py` | Seasonal and regional decomposition of the flux changes, implied change in runoff | Section 3.1, Fig. S4 |
| `an3_spatial.py` | Moran's I, t-test with an effective sample size, spatial block bootstrap, block-level Wilcoxon test, Benjamini-Hochberg FDR | Table 1, Table S3, Section 3.3 |
| `an4_heights.py` | Canopy height inputs by PFT; patch-level height-flux emulator fitted to the uniform-scaling runs | Section 3.2, Table S4, Fig. S6a,b |
| `an5_tiles.py` | Sensitivity of the aggregated heights to the intermediate scale, zero-height pixels and outliers | Section 4.3 |
| `an6_stomatal.py` | Transpiration and GPP responses by PFT | Section 3.1, Fig. S6c |
| `fig_main.py` | Figures 2 to 5 and the maps of Figs. S2 and S3 | Figures 2-5, Figs. S2-S3 |
| `fig_si.py` | Figs. S4 to S6 | Figs. S4-S6 |
| `make_tables.py` | Table 1 and Tables S1 to S4 | `tables/` |
| `text_numbers.py` | Every model-derived number quoted in the text, written to `text_numbers.json` | Abstract, Sections 3 and 4, Conclusions |

## Tables (`tables/`)

| File | Content |
|---|---|
| `Table1_global_differences.csv` | Global relative differences (%) with 95% block-bootstrap intervals (Table 1) |
| `TableS1_global_forest_fluxes.csv` | Global forest fluxes per unit forest area (W m-2) |
| `TableS2_PFT_fluxes_and_differences.csv` | PFT-level fluxes, relative and absolute differences |
| `TableS3_PFT_significance.csv` | PFT-level differences, intervals, effective sample sizes and FDR-adjusted p-values |
| `TableS4_canopy_height_inputs.csv` | Canopy height inputs by PFT and the area with the default height retained |
| `TableS5_CLM5_height_parameters.csv` | Height-related parameters of the forest PFTs (`clm50_params.c240208.nc`) |
