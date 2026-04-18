# Structural Controls on Forest Evapotranspiration from Canopy Height

Code and analysis for the manuscript:

**"Structural Controls on Forest Evapotranspiration from Canopy Height"**

Hangkai You, Fa Li, Dalei Hao, Benjamin Dechant, Zhenpeng Zuo, Taejin Park, Meng Luo, Min Chen*

*Corresponding author: min.chen@wisc.edu*

Target journal: Global Change Biology

## Abstract

Forest canopy height is a fundamental control on aerodynamic conductance and surface energy partitioning, but its influence on evapotranspiration (ET) at a global scale remains under-investigated. We combined biophysical theory with Community Land Model version 5 (CLM5) simulations to quantify how canopy height affects the magnitude and partitioning of global forest ET. By uniformly scaling canopy height by up to +/-20%, we first establish the intrinsic sensitivity of ET to canopy structure: a 20% increase in canopy height reduces global transpiration by 4.15% while increasing canopy evaporation by 2.94%, as enhanced aerodynamic conductance shifts energy partitioning between latent and sensible heat fluxes. We then prescribe canopy height fields derived from GEDI (circa 2020) and ICESat-1 GLAS (2003-2009) to assess how differences in data source, observation epoch, and spatial aggregation statistics (maximum, mean, median) propagate into simulated ET. The ET response to canopy height is strongly biome dependent: tropical forests are comparatively insensitive despite large height contrasts, whereas semi-arid and temperate forests show much stronger responses because aerodynamic effects interact with moisture limitation and stomatal regulation. Differences among canopy height datasets and aggregation statistics produced up to 1.76% divergence in global ET and up to 4.93% in individual flux components such as transpiration, with substantially larger differences at regional and biome scales, revealing a previously overlooked source of structural uncertainty. Our results highlight the need for updated canopy height constraints and standardized aggregation protocols to reduce structural uncertainty in the representation of vegetation-atmosphere coupling in Earth system models.

## Model Configuration

- **Model**: CLM5 (Community Land Model version 5) in Satellite Phenology (SP) mode
- **Resolution**: 1.9 x 2.5 degrees
- **Forcing**: GSWP3 (3-hourly, 2001-2010, looped every 10 years)
- **Run length**: 40 years (30-year spinup + 10-year analysis)
- **Analysis period**: Last 10 years (climatological mean over one full forcing cycle)
- **Height modification**: Both canopy top height (`MONTHLY_HEIGHT_TOP`) and canopy bottom height (`MONTHLY_HEIGHT_BOT`) are modified in the CLM5 surface dataset. Bottom height is derived from top height using a fixed PFT-specific ratio, preserving the default relationship for each PFT.
- **Scenarios** (8 total):
  - **Height scaling**: Default, 0.8x, 0.9x, 1.1x, 1.2x (uniform scaling of all forest PFT heights)
  - **Height aggregation**: GEDI Max, GEDI Mean, GEDI Median (ETH canopy height aggregated at 500 m)

## Canopy Height Data

- **Default (ICESat)**: Simard et al. (2011), based on ICESat-1 GLAS (2003-2009)
- **GEDI-derived**: Lang et al. (2023), ETH Global Canopy Height at 10 m resolution (circa 2020)
- Upscaled to CLM5 grid using three aggregation methods (max, mean, median) per PFT
- Height thresholds applied during aggregation: > 0 m and < 90 m (excludes non-vegetated pixels and outliers)

## Repository Structure

### Python Analysis Scripts

| Script | Description |
|---|---|
| `generate_maps_v3.py` | Spatial maps: experimental scenarios (CLM Default absolute + % diff), multisource scenarios, absolute difference (W/m2), absolute values (mm/month) |
| `generate_heatmaps_v2.py` | PFT-level heatmaps (2x2) with paired t-test significance markers |
| `generate_tables.py` | Global and PFT-scale statistical summary tables (CSV) |
| `generate_all_figures.py` | Combined figure generation and old-vs-new comparison analysis |

### Canopy Height Preprocessing

| Script | Description |
|---|---|
| `canopy_height_aggregation.py` | Aggregates 10 m ETH canopy height to CLM5 grid per PFT using max, mean, and median. Applies height thresholds (>0 m, <90 m). Modifies both `MONTHLY_HEIGHT_TOP` and `MONTHLY_HEIGHT_BOT` in CLM5 surface data. |

### Legacy Notebooks

| Notebook | Description |
|---|---|
| `NC_File_Averaging.ipynb` | Averaging CLM5 output NC files |
| `Percentage_Calculation.ipynb` | Global area-weighted percentage differences |
| `Figure_1_Canopy_Height.ipynb` | Forest canopy height comparison across datasets |
| `Figure_3_Sensitivity_Maps.ipynb` | Spatial maps of ET sensitivity |
| `Figure_4_ET_Spatial_Maps.ipynb` | ET differences across GEDI vs. CLM Default |
| `Figure_5_Bar_Chart.ipynb` | Bar chart of ET component differences across PFTs |

### Output (`figures_new_version/`)

- Spatial maps (TIFF): experimental scenarios, multisource scenarios, absolute values, absolute differences
- PFT heatmaps (PNG): percentage differences, absolute differences, absolute values with significance markers
- Statistical tables (CSV): global and PFT-level absolute values, percentage and absolute differences

## Statistical Methods

- **Global comparisons**: Area-weighted means across forested grid cells
- **PFT-level comparisons**: Weighted by PFT fractional area (`pfts1d_wtgcell`) and grid cell area
- **Significance testing**: Paired Wilcoxon signed-rank tests (global) and paired t-tests (PFT-level) using 10 annual means from the analysis period
- **Spatial variability**: Interquartile range (IQR) of grid-cell-level relative differences

## Variables Analyzed

| Variable | Description | Units |
|---|---|---|
| FCEV | Canopy evaporation | W/m2 |
| FCTR | Canopy transpiration | W/m2 |
| FGEV | Ground evaporation | W/m2 |
| ET | Total evapotranspiration (FCEV + FCTR + FGEV) | W/m2 |
| FSH | Sensible heat flux | W/m2 |
