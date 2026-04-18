# CLM5 Output Data

h0 monthly output files averaged over the last 10 years (years 31-40) of
the 40-year CLM5 simulations. Each file contains PFT-level and grid-level
fluxes (FCEV, FCTR, FGEV, FSH, etc.) at 1.9 x 2.5 degree resolution.

## Files

| File | Scenario |
|---|---|
| `Default_case_h0_avg_last10yr.nc` | CLM Default (ICESat-1 GLAS canopy heights) |
| `Default_0.8_case_h0_avg_last10yr.nc` | Default x 0.8 (uniform height scaling) |
| `Default_0.9_case_h0_avg_last10yr.nc` | Default x 0.9 |
| `Default_1.1_case_h0_avg_last10yr.nc` | Default x 1.1 |
| `Default_1.2_case_h0_avg_last10yr.nc` | Default x 1.2 |
| `Mean_case_h0_avg_last10yr.nc` | GEDI Mean (ETH 10m aggregated by mean) |
| `Max_case_h0_avg_last10yr.nc` | GEDI Max (ETH 10m aggregated by max) |
| `Median_case_h0_avg_last10yr.nc` | GEDI Median (ETH 10m aggregated by median) |

## Format

- Dimensions: `pft=48932` (1D PFT vector), `lat=96`, `lon=144`
- Time-averaged (no time dimension)
- All 92 CLM5 output variables preserved from the raw h0 files

## Surface Data Input Files

The CLM5 surface data input files (`surfdata_1.9x2.5_*_modified.nc`, 140 MB each,
1.1 GB total) containing the modified `MONTHLY_HEIGHT_TOP` and `MONTHLY_HEIGHT_BOT`
fields are not included in this repository due to GitHub file size limits.
They are available from the corresponding author on request.
