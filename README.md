# Structural Controls on Forest Evapotranspiration from Canopy Height

Code and analysis for the manuscript:

**"Structural Controls on Forest Evapotranspiration from Canopy Height"**

Hangkai You, Fa Li, Dalei Hao, Benjamin Dechant, Zhenpeng Zuo, Taejin Park, Meng Luo, Min Chen*

*Corresponding author: min.chen@wisc.edu*

The tables, figures and model-derived numbers of the revised manuscript are produced by the scripts listed under [Figures and tables of the revised manuscript](#figures-and-tables-of-the-revised-manuscript). Most of them are in [`R1_analysis/`](R1_analysis/); the maps of Figures 2 and 4 and Figs. S2 and S3 are drawn by scripts in the repository root. The other scripts and notebooks in the repository root are those of the first submission and of earlier drafts (see [First-submission scripts](#first-submission-scripts)).

## Abstract

Forest canopy height is a fundamental control on aerodynamic conductance and surface energy partitioning, but its global influence on evapotranspiration (ET) remains under-investigated. We combined a Penman–Monteith-based analysis with controlled sensitivity experiments in the Community Land Model version 5 (CLM5) to quantify how prescribed canopy height affects the magnitude and partitioning of global forest ET. Uniformly scaling canopy height by up to ±20% in a fixed CLM5 configuration showed that a 20% canopy height increase reduces global forest transpiration by 4.35% and increases canopy and ground evaporation by 2.87% and 2.14%, reducing total forest ET by 1.94%. We then prescribed canopy height fields derived from GEDI (circa 2020), aggregated with the maximum, mean and median, and compared them with the default CLM5 heights derived from ICESat-1 GLAS (2005). Differences among canopy height inputs changed global forest ET and transpiration by up to 1.99% and 4.95%, and the choice of aggregation statistic alone produced differences comparable to a uniform 20% change in canopy height. Every global difference had the same sign in each of ten forcing years, and all but one remained significant after accounting for spatial autocorrelation. Responses differed among plant functional types: for the same relative change in canopy height, transpiration and total ET responded most strongly in tropical broadleaf evergreen forests (-5.7% and -2.7% for a 20% increase), whereas in needleleaf forests the changes in evaporation and transpiration largely offset, leaving total ET nearly unchanged. These results quantify the sensitivity of a land surface model to canopy height inputs and identify the data source and aggregation of canopy height as a source of structural uncertainty in simulated ET that should be tested across models.

## Figures and tables of the revised manuscript

| Item | Produced by | Output |
|---|---|---|
| Figure 1 | Schematic (not drawn from model output) | |
| Figure 2 | `Figure_2_Sensitivity_Maps.ipynb` | Notebook output |
| Figure 3 | `R1_analysis/fig_main.py` | `Figure3_canopy_height_inputs` |
| Figure 4 | `generate_all_figures.py` | [`figures_new_version/Figure_4_ET_comparison_maps.png`](figures_new_version/Figure_4_ET_comparison_maps.png) |
| Figure 5 | `R1_analysis/fig_main.py` | `Figure5_PFT_bars` |
| Table 1 | `R1_analysis/an3_spatial.py`, `R1_analysis/make_tables.py` | [`R1_analysis/tables/Table1_global_differences.csv`](R1_analysis/tables/Table1_global_differences.csv) |
| Fig. S1 | Schematic (not drawn from model output) | |
| Fig. S2 | `generate_maps_v3.py` | [`figures_new_version/Figure_absolute_values.tiff`](figures_new_version/Figure_absolute_values.tiff) |
| Fig. S3 | `generate_maps_v3.py` | [`figures_new_version/Figure_absolute_difference.tiff`](figures_new_version/Figure_absolute_difference.tiff) |
| Fig. S4 | `R1_analysis/fig_si.py` | `FigureS4_canopy_water_balance_seasonal` |
| Fig. S5 | `R1_analysis/fig_si.py` | `FigureS5_forcing_year_robustness` |
| Fig. S6 | `R1_analysis/fig_si.py` | `FigureS6_emulator_and_stomatal` |
| Tables S1 to S4 | `R1_analysis/make_tables.py` | [`R1_analysis/tables/`](R1_analysis/tables/) |
| Table S5 | Parameter values of `clm50_params.c240208.nc` (not computed by the scripts) | [`R1_analysis/tables/TableS5_CLM5_height_parameters.csv`](R1_analysis/tables/TableS5_CLM5_height_parameters.csv) |
| Numbers in the text | `R1_analysis/text_numbers.py` and the analysis scripts (see [`R1_analysis/README.md`](R1_analysis/README.md)) | `text_numbers.json` |

The maps of Figures 2 and 4 and Figs. S2 and S3 show the values of individual grid cells, which do not depend on the area weights of the global and PFT-level means (see [Statistical Methods](#statistical-methods-revised-manuscript)). `generate_maps_v3.py` also writes a version of `Figure_4_ET_comparison_maps` with panel letters and longer row labels; Figure 4 is the version written by `generate_all_figures.py`.

## Model Configuration

- **Model**: CLM5 (Community Land Model version 5) in Satellite Phenology (SP) mode: CTSM tag `ctsm5.2.005` with the CLM5.0 physics option (compset `2000_DATM%GSWP3v1_CLM50%SP_SICE_SOCN_MOSART_SGLC_SWAV_SESP`, parameter file `clm50_params.c240208.nc`). The model code was used unchanged except for one compatibility change in `src/main/surfrdMod.F90` ([`model_mods/surfrdMod.F90.diff`](model_mods/surfrdMod.F90.diff)), which sets `PCT_OCEAN` to zero when it is missing from the surface dataset, as it is in the CLM5.0-format datasets used here.
- **Resolution**: 1.9 x 2.5 degrees
- **Forcing**: GSWP3 v1 (3-hourly), with the 2001-2010 block cycled four times
- **Run length**: 40 years from a cold start (30-year spinup + 10-year analysis)
- **Analysis period**: Last 10 years (the fourth 2001-2010 forcing cycle)
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

- **Default (ICESat)**: Simard et al. (2011), a global 1 km map produced from 2005 observations of ICESat-1 GLAS (RH100 waveform metric)
- **GEDI-derived**: Lang et al. (2023), ETH Global Canopy Height at 10 m resolution, from GEDI RH98 observations and Sentinel-2 imagery (circa 2020)
- **PFT map**: 500 m, derived from the MODIS land cover type product and WorldClim V1 climate data on Google Earth Engine using [build_surface_dataset_for_ELM](https://github.com/daleihao/build_surface_dataset_for_ELM) (`ELM_PFT_{year}-WGS84-merged.tif`; PFTs 1-8 are trees)
- Two-step aggregation: (1) 10 m canopy height aggregated to the 500 m grid using max, mean, or median; (2) for each forest PFT, the 500 m heights of all pixels of that PFT averaged with equal weight within fixed blocks of 450 x 600 pixels (1.875 x 2.5 degrees at the 1/240 degree pixel size; 96 x 144 blocks, one for each model grid cell). Where a forest PFT of the CLM5 surface dataset was absent from the 500 m PFT map within the block, the CLM5 default height was retained.
- Height threshold applied during Step 1: > 0 m (excludes zero-height pixels and no-data values, so GEDI Mean and GEDI Median describe vegetated pixels)

## Repository Structure

### Revised manuscript (`R1_analysis/`)

Scripts for the results, tables and figures of the revised manuscript other than the maps of Figures 2 and 4 and Figs. S2 and S3: diagnostics and robustness checks, spatial-autocorrelation-aware significance tests, canopy height statistics and the height-flux emulator, the figure scripts, and `make_tables.py`, which writes Table 1 and Tables S1 to S4 to [`R1_analysis/tables/`](R1_analysis/tables/) (the folder also holds Table S5). See the folder's README for the inputs, the run order and which script produces which result.

Global and PFT-level means are weighted by the land area of each patch (`pfts1d_wtgcell` x `area` x `landfrac`). In CLM5 output, `area` is the area of the whole grid cell and `pfts1d_wtgcell` is a fraction of the land part of the cell.

### Canopy Height Preprocessing (`canopy_height_preprocessing/`)

Reference scripts used to aggregate ETH 10 m canopy height to the CLM5
grid per PFT and to modify the CLM5 surface data. See the folder's README
for full workflow details.

| File | Description |
|---|---|
| `input_data_generation_hangkai.py` | Python reference: 500 m canopy height + 500 m PFT map -> CLM5 1.9 x 2.5 grid per PFT; writes modified surfdata. |
| `cht_with_pft.m` | MATLAB equivalent of the same aggregation logic. |

### CLM5 output (`clm_output/`)

Decade means (model years 31-40) of the monthly h0 output of the eight simulations, one file per scenario. See the folder's README.

### First-submission scripts

The scripts and notebooks below produced the results of the first submission and of earlier drafts. Three of them also draw maps of the revised manuscript (see the table above): `Figure_2_Sensitivity_Maps.ipynb` (Figure 2), `generate_all_figures.py` (Figure 4) and `generate_maps_v3.py` (Figs. S2 and S3). Those that compute global or PFT-level means (`generate_tables.py`, `generate_heatmaps_v2.py`, `generate_all_figures.py`, `Percentage_Calculation.ipynb` and `Figure_5_Bar_Chart.ipynb`) weight grid cells by the whole grid-cell `area`, without the land fraction, so their means differ slightly from those of the revised manuscript; use `R1_analysis/` for the revised results. The notebooks marked "earlier draft" are from drafts that used 5-year means (2010-2014) of shorter simulations.

| Script | Description |
|---|---|
| `generate_maps_v3.py` | Spatial maps: absolute values (mm/month; Fig. S2), absolute differences (W/m2; Fig. S3), experimental scenarios (CLM Default absolute + % diff) and multisource scenarios |
| `generate_heatmaps_v2.py` | PFT-level heatmaps (2x2) with paired t-test significance markers |
| `generate_tables.py` | Global and PFT-scale statistical summary tables (CSV) |
| `generate_all_figures.py` | Combined figure generation, including Figure 4 (`Figure_4_ET_comparison_maps`), and old-vs-new comparison analysis |

| Notebook | Description |
|---|---|
| `Figure_2_Sensitivity_Maps.ipynb` | Figure 2: maps of the uniform-scaling experiments per unit forest area |
| `NC_File_Averaging.ipynb` | Averaging CLM5 output NC files |
| `Percentage_Calculation.ipynb` | Global area-weighted percentage differences (earlier draft) |
| `Percentage_Calculation_Original.ipynb` | Original percentage calculation and scenario maps (earlier draft) |
| `Figure_1_Canopy_Height.ipynb` | Forest canopy height comparison across datasets (earlier draft) |
| `Figure_3_Sensitivity_Maps.ipynb` | Spatial maps of ET sensitivity (earlier draft) |
| `Figure_4_ET_Spatial_Maps.ipynb` | ET differences across GEDI vs. CLM Default (earlier draft) |
| `Figure_5_Bar_Chart.ipynb` | Bar chart of ET component differences across PFTs (earlier draft) |
| `Figure_5_PFT_Boxplots.ipynb` | PFT-level box plots of ET components (earlier draft) |

Stored outputs (`figures_new_version/`):

- Maps of the revised manuscript: Figure 4 (`Figure_4_ET_comparison_maps.png`), Fig. S2 (`Figure_absolute_values.tiff`) and Fig. S3 (`Figure_absolute_difference.tiff`)
- Scenario maps of the first submission (TIFF): experimental scenarios and multisource scenarios
- Bar chart (TIFF) and PFT heatmaps (PNG) of the first submission: percentage differences, absolute differences, absolute values with significance markers
- Statistical tables of the first submission (CSV): global and PFT-level absolute values, percentage and absolute differences

## Statistical Methods (revised manuscript)

- **Forested grid cells**: the 3,625 grid cells that contain one or more of the eight forest PFTs
- **Global and PFT-level means**: weighted by patch land area (`pfts1d_wtgcell` x `area` x `landfrac`); relative differences are ratios of these means. Energy fluxes (W/m2) are the CLM5 latent heat outputs; water fluxes (mm/yr) are the CLM5 water-flux outputs (365-day calendar)
- **Spatial autocorrelation**: Moran's I of the grid-cell differences and an effective sample size n_eff = n(1 - ρ)/(1 + ρ), where ρ is the correlation of each cell with the mean of its eight neighbors
- **Significance testing**: t-test with the effective sample size, and a Wilcoxon signed-rank test on block-mean differences, with Benjamini-Hochberg FDR adjustment (α = 0.05); the cell-level Wilcoxon test is reported for reference (Table S3), together with the fraction of forested cells whose difference has the same sign as the global mean
- **Intervals**: 95% intervals of the global relative differences from a spatial block bootstrap (2,000 resamples of the 289 blocks of 5 x 4 grid cells, about 10 x 10 degrees, that contain forest); the block tests were repeated with blocks of about 5 and 20 degrees
- **Robustness**: sign and spread of the differences in each of the ten forcing years, and agreement between the third and fourth forcing cycles
- **Spatial variability**: Interquartile range (IQR) of grid-cell-level relative differences
- **Diagnostics**: canopy water balance (change in canopy evaporation against the change in liquid and snow water reaching the ground, per soil column); relative changes in transpiration and GPP; and a patch-level emulator, a quadratic function of the log height ratio fitted to CLM Default and the four uniform-scaling runs, held constant or extended linearly in log height beyond height ratios of 0.6 and 1.6

## Variables Analyzed

| Variable | Description | Units |
|---|---|---|
| FCEV | Canopy evaporation | W/m2 |
| FCTR | Canopy transpiration | W/m2 |
| FGEV | Ground evaporation | W/m2 |
| ET | Total evapotranspiration (FCEV + FCTR + FGEV) | W/m2 |
| FSH | Sensible heat flux | W/m2 |
| QVEGE, QVEGT, QFLX_EVAP_TOT | Water fluxes: canopy evaporation, transpiration and total ET (ground evaporation = QFLX_EVAP_TOT - QVEGE - QVEGT) | mm/s |
| QFLX_LIQ_GRND, QFLX_SNOW_GRND | Liquid and snow water reaching the ground (canopy water balance) | mm/s |
| FPSN | Photosynthesis (GPP) | µmol m-2 s-1 |
