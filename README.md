# Canopy Height Effects on Evapotranspiration (ET_CTH)

Code repository for the manuscript: **"Structural Controls on Forest Evapotranspiration from Canopy Height in an Earth System Model"**

## Contents

- `NC_File_Averaging.ipynb` — Averaging CLM5 output NC files over the last 5 simulation years (2010–2014)
- `Percentage_Calculation.ipynb` — Global area-weighted percentage differences in ET components across canopy height scenarios
- `Percentage_Calculation_Original.ipynb` — Original percentage calculation notebook with full outputs
- `Figure_1_Canopy_Height.ipynb` — Figure 1: Forest canopy height comparison across datasets
- `Figure_3_Sensitivity_Maps.ipynb` — Figure 3: Spatial maps of ET sensitivity to canopy height changes
- `Figure_4_ET_Spatial_Maps.ipynb` — Figure 4: ET differences across GEDI vs. CLM Default scenarios
- `Figure_5_Bar_Chart.ipynb` — Figure 5: Bar chart of ET component differences across PFTs
- `Figure_5_PFT_Boxplots.ipynb` — Figure 5: PFT-level boxplots of ET components

## Model Configuration

- **Model**: CLM5 (Community Land Model version 5) in Satellite Phenology (SP) mode
- **Resolution**: 1.9° × 2.5°
- **Forcing**: GSWP3 (3-hourly, 2001–2014)
- **Spin-up**: 9 years; Analysis: 5-year average (2010–2014)
- **Scenarios**: CLM Default (ICESat), GEDI Max, GEDI Mean, GEDI Median, ±10%/±20% sensitivity
