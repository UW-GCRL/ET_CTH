"""
Canopy Height Aggregation for CLM5 Surface Data

Two-step aggregation of 10 m ETH Global Canopy Height (Lang et al., 2023)
to CLM5 grid resolution (1.9 x 2.5 degrees) per PFT:

  Step 1: 10 m canopy height -> 500 m PFT grid
          Within each 500 m PFT cell, aggregate all valid 10 m pixels
          (height > 0 m and < 90 m) using max, mean, or median.

  Step 2: 500 m PFT grid -> CLM5 1.9 x 2.5 degree grid
          Average the 500 m aggregated heights within each CLM5 grid cell
          for each PFT.

Both MONTHLY_HEIGHT_TOP and MONTHLY_HEIGHT_BOT are modified. Bottom height
is derived from top height using a fixed PFT-specific ratio, preserving the
default relationship for each PFT.

Usage:
    python canopy_height_aggregation.py \
        --canopy_height <path_to_10m_canopy_height.tif> \
        --pft_map <path_to_500m_pft_map.tif> \
        --surfdata <path_to_clm5_surfdata.nc> \
        --output_dir <output_directory>

Outputs:
    - Max_surfdata.nc   (aggregation using maximum height at 500 m)
    - Mean_surfdata.nc  (aggregation using mean height at 500 m)
    - Median_surfdata.nc (aggregation using median height at 500 m)
"""
import numpy as np
import argparse
import os
import sys

try:
    import rasterio
except ImportError:
    rasterio = None

try:
    import netCDF4
except ImportError:
    netCDF4 = None


# CLM5 grid dimensions (1.9 x 2.5 degrees)
NLAT = 96
NLON = 144

# Forest PFTs (1-indexed as in CLM5)
FOREST_PFTS = {
    1: "NET Temperate",
    2: "NET Boreal",
    3: "NDT Boreal",
    4: "BET Tropical",
    5: "BET Temperate",
    6: "BDT Tropical",
    7: "BDT Temperate",
    8: "BDT Boreal",
}

# Height thresholds (meters)
HEIGHT_MIN = 0.0   # exclusive: height > 0
HEIGHT_MAX = 90.0  # exclusive: height < 90


def load_raster(filepath):
    """Load a GeoTIFF raster and return data + transform info."""
    if rasterio is None:
        raise ImportError("rasterio is required. Install with: pip install rasterio")
    with rasterio.open(filepath) as src:
        data = src.read(1).astype(np.float32)
        transform = src.transform
        height, width = data.shape
    return data, height, width, transform


def step1_aggregate_10m_to_500m(canopy_height_10m, pft_map_500m, pft_id, method='mean'):
    """
    Step 1: Aggregate 10 m canopy height to 500 m PFT grid.

    For each 500 m PFT pixel that matches pft_id, collect all overlapping
    10 m canopy height pixels with valid height (>0 m and <90 m), then
    compute the max, mean, or median.

    Parameters
    ----------
    canopy_height_10m : 2D array, canopy height in meters (10 m resolution)
    pft_map_500m : 2D array, PFT classification (500 m resolution)
    pft_id : int, PFT index to aggregate
    method : str, one of 'max', 'mean', 'median'

    Returns
    -------
    height_500m : 2D array (same shape as pft_map_500m), aggregated height
    """
    pft_h, pft_w = pft_map_500m.shape
    ch_h, ch_w = canopy_height_10m.shape

    # Ratio of 10 m pixels per 500 m pixel (500/10 = 50 pixels per side)
    ratio_y = ch_h / pft_h
    ratio_x = ch_w / pft_w

    height_500m = np.zeros((pft_h, pft_w), dtype=np.float32)

    for j in range(pft_h):
        if j % (pft_h // 10) == 0:
            print(f"    {int(j / pft_h * 100)}%...")

        for i in range(pft_w):
            # Only process pixels matching the target PFT
            if pft_map_500m[j, i] != pft_id:
                continue

            # Determine the 10 m pixel range for this 500 m cell
            y_start = int(j * ratio_y)
            y_end = int((j + 1) * ratio_y)
            x_start = int(i * ratio_x)
            x_end = int((i + 1) * ratio_x)

            # Clip to bounds
            y_end = min(y_end, ch_h)
            x_end = min(x_end, ch_w)

            # Extract 10 m pixels within this 500 m cell
            block = canopy_height_10m[y_start:y_end, x_start:x_end]

            # Apply height thresholds
            valid = block[(block > HEIGHT_MIN) & (block < HEIGHT_MAX)]

            if len(valid) == 0:
                continue

            if method == 'max':
                height_500m[j, i] = np.max(valid)
            elif method == 'mean':
                height_500m[j, i] = np.mean(valid)
            elif method == 'median':
                height_500m[j, i] = np.median(valid)

    return height_500m


def step2_aggregate_500m_to_clm(height_500m, pft_map_500m, pft_id):
    """
    Step 2: Aggregate 500 m PFT height grid to CLM5 1.9 x 2.5 degree grid.

    For each CLM5 grid cell, compute the mean of all valid 500 m pixels
    that match the target PFT.

    Parameters
    ----------
    height_500m : 2D array, aggregated canopy height at 500 m
    pft_map_500m : 2D array, PFT classification at 500 m
    pft_id : int, PFT index

    Returns
    -------
    grid : 2D array (NLON, NLAT), mean height per CLM5 grid cell
    count : 2D array (NLON, NLAT), number of valid 500 m pixels per cell
    """
    pft_h, pft_w = pft_map_500m.shape

    grid_sum = np.zeros((NLON, NLAT))
    grid_count = np.zeros((NLON, NLAT), dtype=np.int64)

    # Mask: valid height and matching PFT
    valid = (height_500m > 0) & (pft_map_500m == pft_id)
    rows, cols = np.where(valid)

    for idx in range(len(rows)):
        j, i = rows[idx], cols[idx]
        # Map 500 m pixel to CLM5 grid cell
        gx = int((i / pft_w) * NLON)
        gy = int((j / pft_h) * NLAT)
        gx = min(gx, NLON - 1)
        gy = min(gy, NLAT - 1)

        grid_sum[gx, gy] += height_500m[j, i]
        grid_count[gx, gy] += 1

    # Compute mean
    grid = np.zeros((NLON, NLAT))
    mask = grid_count > 0
    grid[mask] = grid_sum[mask] / grid_count[mask]

    # Flip latitude (raster convention -> CLM convention)
    grid = np.fliplr(grid)
    grid_count = np.fliplr(grid_count)

    return grid, grid_count


def update_surfdata(surfdata_path, output_path, height_grids):
    """
    Update CLM5 surface data with new canopy heights.

    Both MONTHLY_HEIGHT_TOP and MONTHLY_HEIGHT_BOT are updated.
    Bottom height is derived from top height using the original
    PFT-specific ratio (bottom/top) from the default surface data.

    Parameters
    ----------
    surfdata_path : str, path to original surfdata NC file
    output_path : str, path to output NC file
    height_grids : dict, {pft_id: 2D array (NLON, NLAT) of heights}
    """
    if netCDF4 is None:
        raise ImportError("netCDF4 is required. Install with: pip install netCDF4")

    import shutil
    shutil.copy2(surfdata_path, output_path)

    with netCDF4.Dataset(output_path, 'r+') as ds:
        htop = ds.variables['MONTHLY_HEIGHT_TOP'][:]
        hbot = ds.variables['MONTHLY_HEIGHT_BOT'][:]

        for pft_id, new_height in height_grids.items():
            # PFT dimension index (bare soil=0, so forest PFT 1 = index 1)
            pft_idx = pft_id

            for a in range(NLON):
                for b in range(NLAT):
                    if htop[a, b, pft_idx, 0] > 0 and new_height[a, b] > 0:
                        # Preserve bottom/top ratio from original data
                        old_top = htop[a, b, pft_idx, 0]
                        old_bot = hbot[a, b, pft_idx, 0]
                        ratio = old_bot / old_top if old_top > 0 else 0

                        new_top = new_height[a, b]
                        new_bot = new_top * ratio

                        # Apply to all 12 months
                        for m in range(12):
                            htop[a, b, pft_idx, m] = new_top
                            hbot[a, b, pft_idx, m] = new_bot

        ds.variables['MONTHLY_HEIGHT_TOP'][:] = htop
        ds.variables['MONTHLY_HEIGHT_BOT'][:] = hbot

    print(f"  Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Two-step aggregation of 10m canopy height to CLM5 grid per PFT"
    )
    parser.add_argument('--canopy_height', required=True,
                        help='Path to 10m canopy height GeoTIFF')
    parser.add_argument('--pft_map', required=True,
                        help='Path to 500m PFT map GeoTIFF (e.g. ELM_PFT_{year}-WGS84-merged.tif from build_surface_dataset_for_ELM)')
    parser.add_argument('--surfdata', required=True,
                        help='Path to CLM5 surface data NC file')
    parser.add_argument('--output_dir', required=True,
                        help='Output directory for modified surfdata files')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("Loading 10m canopy height raster...")
    canopy_height_10m, ch_h, ch_w, _ = load_raster(args.canopy_height)
    print(f"  Shape: {canopy_height_10m.shape}")

    print("Loading 500m PFT map...")
    pft_map_500m, pft_h, pft_w, _ = load_raster(args.pft_map)
    print(f"  Shape: {pft_map_500m.shape}")

    methods = ['max', 'mean', 'median']

    for method in methods:
        print(f"\n{'='*60}")
        print(f"Aggregation method: {method.upper()}")
        print(f"{'='*60}")
        height_grids = {}

        for pft_id, pft_name in FOREST_PFTS.items():
            print(f"\n  PFT {pft_id} ({pft_name}):")

            # Step 1: 10m -> 500m
            print(f"  Step 1: Aggregating 10m -> 500m ({method})...")
            height_500m = step1_aggregate_10m_to_500m(
                canopy_height_10m, pft_map_500m, pft_id, method=method
            )
            n_valid_500m = np.sum(height_500m > 0)
            if n_valid_500m > 0:
                print(f"    Valid 500m pixels: {n_valid_500m}, "
                      f"Mean: {np.mean(height_500m[height_500m > 0]):.1f} m")
            else:
                print(f"    No valid pixels found.")

            # Step 2: 500m -> CLM5 grid
            print(f"  Step 2: Aggregating 500m -> CLM5 grid (mean)...")
            grid, count = step2_aggregate_500m_to_clm(
                height_500m, pft_map_500m, pft_id
            )
            n_valid_clm = np.sum(count > 0)
            if n_valid_clm > 0:
                print(f"    Valid CLM5 cells: {n_valid_clm}, "
                      f"Mean: {np.mean(grid[grid > 0]):.1f} m")
            else:
                print(f"    No valid CLM5 cells.")

            height_grids[pft_id] = grid

        # Step 3: Update surface data
        output_path = os.path.join(args.output_dir, f"{method.capitalize()}_surfdata.nc")
        print(f"\n  Step 3: Updating CLM5 surface data...")
        update_surfdata(args.surfdata, output_path, height_grids)

    print(f"\n{'='*60}")
    print("Done! Output files:")
    for method in methods:
        print(f"  {os.path.join(args.output_dir, f'{method.capitalize()}_surfdata.nc')}")


if __name__ == '__main__':
    main()
