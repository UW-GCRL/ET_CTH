# Structural Controls on Forest Evapotranspiration from Canopy Height

Code and analysis for the manuscript:

**"Structural Controls on Forest Evapotranspiration from Canopy Height"**

Hangkai You, Fa Li, Dalei Hao, Benjamin Dechant, Zhenpeng Zuo, Taejin Park, Meng Luo, Min Chen*

*Corresponding author: min.chen@wisc.edu*

The results, tables and figures of the revised manuscript are produced by the scripts in [`R1_analysis/`](R1_analysis/). The scripts in the repository root and the notebooks are those of the first submission (see [First-submission scripts](#first-submission-scripts)).

## Abstract

Forest canopy height is a fundamental control on aerodynamic conductance and surface energy partitioning, but its global influence on evapotranspiration (ET) remains under-investigated. We combined a Penman–Monteith-based analysis with controlled sensitivity experiments in the Community Land Model version 5 (CLM5) to quantify how prescribed canopy height affects the magnitude and partitioning of global forest ET. Uniformly scaling canopy height by up to ±20% in a fixed CLM5 configuration showed that a 20% canopy height increase reduces global forest transpiration by 4.35% and increases canopy and ground evaporation by 2.87% and 2.14%, reducing total forest ET by 1.94%. We then prescribed canopy height fields derived from GEDI (circa 2020), aggregated with the maximum, mean and median, and compared them with the default CLM5 heights derived from ICESat-1 GLAS (2005). Differences among canopy height inputs changed global forest ET and transpiration by up to 1.99% and 4.95%, and the choice of aggregation statistic alone produced differences comparable to a uniform 20% change in canopy height. Every global difference had the same sign in each of ten forcing years, and all but one remained significant after accounting for spatial autocorrelation. Responses differed among plant functional types: for the same relative change in canopy height, transpiration and total ET responded most strongly in tropical broadleaf evergreen forests (−5.7% and −2.7% for a 20% increase), whereas in needleleaf forests the changes in evaporation and transpiration largely offset, leaving total ET nearly unchanged. These results quantify the sensitivity of a land surface model to canopy height inputs and identify the data source and aggregation of canopy height as a source of structural uncertainty in simulated ET that should be tested across models.

## Model Configuration

- **Model**: CLM5 (Community Land Model version 5) in Satellite Phenology (SP) mode: CTSM tag `ctsm5.2.005` with the CLM5.0 physics option (compset `2000_DATM%GSWP3v1_CLM50%SP_SICE_SOCN_MOSART_SGLC_SWAV_SESP`, parameter file `clm50_params.c240208.nc`). The model code was used unchanged except for one compatibility change in `src/main/surfrdMod.F90` ([`model_mods/surfrdMod.F90.diff`](model_mods/surfrdMod.F90.diff)), which sets `PCT_OCEAN` to zero when it is missing from the surface dataset, as it is in the CLM5.0-format datasets used here.
- **Resolution**: 1.9 x 2.5 degrees
- **Forcing**: GSWP3 (3-hourly, 2001-2010, looped every 10 years)
- **Run length**: 40 years (30-year spinup + 10-year analysis)
- **Analysis period**: Last 10 years (climatological mean over one full forcing cycle)
- **Height modification**: Both canopy top height (`MONTHLY_HEIGHT_TOP`) and canopy bottom height (`MONTHLY_HEIGHT_BOT`) are modified in the CLM5 surface dataset. Bottom height is derived from top height using a fixed PFT-specific ratio, preserving the default relationship for each PFT:

  | PFT | Name | HBOT/HTOP ratio |
  |---|---|---|
  | 1 | NET Temperate | 0.5000 |
  | 2 | NET Boreal | 0.5000 |
  | 3 | NDT Boreal | 0.5000 |
  | 4 | BET Tropical | 0.0286 |
  | 5 | BET Temperate | 0.0286 |
  | 6 | BDT Tropical | 0.5556 |
  | 7 | BDT Temperate | 0.5750 |
  | 8 | BDT Boreal | 0.5750 |
- **Scenarios** (8 total):
  - **Height scaling**: Default, 0.8x, 0.9x, 1.1x, 1.2x (uniform scaling of the heights of all natural vegetation PFTs: trees, shrubs and grasses)
  - **Height aggregation**: GEDI Max, GEDI Mean, GEDI Median (ETH canopy height aggregated at 500 m; only the heights of the eight forest PFTs are replaced)

## Canopy Height Data

- **Default (ICESat)**: Simard et al. (2011), based on ICESat-1 GLAS (2003-2009)
- **GEDI-derived**: Lang et al. (2023), ETH Global Canopy Height at 10 m resolution (circa 2020)
- **PFT map**: 500 m WGS84, derived from MODIS land cover using [build_surface_dataset_for_ELM](https://github.com/daleihao/build_surface_dataset_for_ELM) (`ELM_PFT_{year}-WGS84-merged.tif`; PFTs 1-8 are trees)
- Two-step aggregation: (1) 10 m canopy height aggregated to the 500 m grid using max, mean, or median; (2) for each forest PFT, the 500 m heights of all pixels of that PFT averaged with equal weight within fixed blocks of 450 x 600 pixels (1.875 x 2.5 degrees at the 1/240 degree pixel size; 96 x 144 blocks, one for each model grid cell). Where a forest PFT of the CLM5 surface dataset was absent from the 500 m PFT map within the block, the CLM5 default height was retained.
- Height threshold applied during Step 1: > 0 m (excludes non-vegetated pixels and no-data values)

## Repository Structure

### Revised manuscript (`R1_analysis/`)

Scripts for all results, tables and figures of the revised manuscript: diagnostics and robustness checks, spatial-autocorrelation-aware significance tests, canopy height statistics and the height-flux emulator, the figure scripts, and `make_tables.py`, which writes Table 1 and Tables S1 to S4 to [`R1_analysis/tables/`](R1_analysis/tables/). See the folder's README for the inputs, the run order and which script produces which result.

Global and PFT-level means are weighted by the land area of each patch (`pfts1d_wtgcell` x `area` x `landfrac`). In CLM5 output, `area` is the area of the whole grid cell and `pfts1d_wtgcell` is a fraction of the land part of the cell.

### Canopy Height Preprocessing (`canopy_height_preprocessing/`)

Reference scripts used to aggregate ETH 10 m canopy height to the CLM5
grid per PFT and to modify the CLM5 surface data. See the folder's README
for full workflow details.

| File | Description |
|---|---|
| `input_data_generation_hangkai.py` | Python reference: 500 m canopy height + 500 m PFT map -> CLM5 1.9 x 2.5 grid per PFT; writes modified surfdata. |
| `cht_with_pft.m` | MATLAB equivalent of the same aggregation logic. |

### First-submission scripts

The scripts and notebooks below produced the results of the first submission. Those that compute global or PFT-level means (`generate_tables.py`, `generate_heatmaps_v2.py`, `generate_all_figures.py`, `Percentage_Calculation.ipynb` and `Figure_5_Bar_Chart.ipynb`) weight grid cells by the whole grid-cell `area`, without the land fraction, so their means differ slightly from those of the revised manuscript; use `R1_analysis/` for the revised results.

| Script | Description |
|---|---|
| `generate_maps_v3.py` | Spatial maps: experimental scenarios (CLM Default absolute + % diff), multisource scenarios, absolute difference (W/m2), absolute values (mm/month) |
| `generate_heatmaps_v2.py` | PFT-level heatmaps (2x2) with paired t-test significance markers |
| `generate_tables.py` | Global and PFT-scale statistical summary tables (CSV) |
| `generate_all_figures.py` | Combined figure generation and old-vs-new comparison analysis |

| Notebook | Description |
|---|---|
| `NC_File_Averaging.ipynb` | Averaging CLM5 output NC files |
| `Percentage_Calculation.ipynb` | Global area-weighted percentage differences |
| `Figure_1_Canopy_Height.ipynb` | Forest canopy height comparison across datasets |
| `Figure_3_Sensitivity_Maps.ipynb` | Spatial maps of ET sensitivity |
| `Figure_4_ET_Spatial_Maps.ipynb` | ET differences across GEDI vs. CLM Default |
| `Figure_5_Bar_Chart.ipynb` | Bar chart of ET component differences across PFTs |

Output of the first submission (`figures_new_version/`):

- Spatial maps (TIFF): experimental scenarios, multisource scenarios, absolute values, absolute differences
- PFT heatmaps (PNG): percentage differences, absolute differences, absolute values with significance markers
- Statistical tables (CSV): global and PFT-level absolute values, percentage and absolute differences

## Statistical Methods (revised manuscript)

- **Global and PFT-level means**: weighted by patch land area (`pfts1d_wtgcell` x `area` x `landfrac`); relative differences are ratios of these means
- **Spatial autocorrelation**: Moran's I of the grid-cell differences and an effective sample size from the correlation of each cell with the mean of its eight neighbors
- **Significance testing**: t-test with the effective sample size, and a Wilcoxon signed-rank test on block-mean differences (blocks of about 10 x 10 degrees), with Benjamini-Hochberg FDR adjustment; the cell-level Wilcoxon test is reported for reference
- **Intervals**: 95% intervals of the global relative differences from a spatial block bootstrap (2,000 resamples of 289 blocks)
- **Robustness**: sign and spread of the differences in each of the ten forcing years, and agreement between the third and fourth forcing cycles
- **Spatial variability**: Interquartile range (IQR) of grid-cell-level relative differences

## Variables Analyzed

| Variable | Description | Units |
|---|---|---|
| FCEV | Canopy evaporation | W/m2 |
| FCTR | Canopy transpiration | W/m2 |
| FGEV | Ground evaporation | W/m2 |
| ET | Total evapotranspiration (FCEV + FCTR + FGEV) | W/m2 |
| FSH | Sensible heat flux | W/m2 |
