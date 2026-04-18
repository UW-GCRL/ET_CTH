"""
Canopy Height Aggregation for CLM5 Surface Data

Aggregates 10 m ETH Global Canopy Height (Lang et al., 2023) to CLM5 grid
resolution (1.9 x 2.5 degrees) per PFT using max, mean, and median statistics.

Height thresholds:
  - Lower bound: > 0 m (excludes non-vegetated / no-data pixels)
  - Upper bound: < 90 m (excludes unrealistic outliers)

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
    - Max_surfdata.nc   (aggregation using maximum height)
    - Mean_surfdata.nc  (aggregation using mean height)
    - Median_surfdata.nc (aggregation using median height)
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


def aggregate_to_clm_grid(canopy_height, pft_map, pft_id, method='mean'):
    """
    Aggregate 10m canopy height to CLM5 grid for a specific PFT.

    Parameters
    ----------
    canopy_height : 2D array, canopy height in meters (10m resolution)
    pft_map : 2D array, PFT classification (same resolution as canopy_height)
    pft_id : int, PFT index to aggregate
    method : str, one of 'max', 'mean', 'median'

    Returns
    -------
    grid : 2D array (NLAT, NLON), aggregated height values
    count : 2D array (NLAT, NLON), number of valid pixels per grid cell
    """
    nrows, ncols = canopy_height.shape
    grid = np.zeros((NLON, NLAT))
    count = np.zeros((NLON, NLAT), dtype=np.int64)

    # For max/median, we need to collect all values per grid cell
    if method in ('max', 'median'):
        cell_values = [[[] for _ in range(NLAT)] for _ in range(NLON)]

    # Mask: valid height and matching PFT
    valid = (canopy_height > HEIGHT_MIN) & (canopy_height < HEIGHT_MAX) & (pft_map == pft_id)

    # Map pixels to CLM grid cells
    rows, cols = np.where(valid)
    for idx in range(len(rows)):
        i, j = rows[idx], cols[idx]
        # Map to CLM grid (assuming global coverage)
        gx = int((j / ncols) * NLON)
        gy = int((i / nrows) * NLAT)
        gx = min(gx, NLON - 1)
        gy = min(gy, NLAT - 1)

        h = canopy_height[i, j]

        if method == 'mean':
            grid[gx, gy] += h
            count[gx, gy] += 1
        elif method in ('max', 'median'):
            cell_values[gx][gy].append(h)
            count[gx, gy] += 1

    # Finalize
    if method == 'mean':
        mask = count > 0
        grid[mask] = grid[mask] / count[mask]
    elif method == 'max':
        for gx in range(NLON):
            for gy in range(NLAT):
                if cell_values[gx][gy]:
                    grid[gx, gy] = np.max(cell_values[gx][gy])
    elif method == 'median':
        for gx in range(NLON):
            for gy in range(NLAT):
                if cell_values[gx][gy]:
                    grid[gx, gy] = np.median(cell_values[gx][gy])

    # Flip latitude (raster convention -> CLM convention)
    grid = np.fliplr(grid)
    count = np.fliplr(count)

    return grid, count


def update_surfdata(surfdata_path, output_path, height_grids):
    """
    Update CLM5 surface data with new canopy heights.

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
            # CLM PFT index is 0-based in the array, but pft_id is 1-based
            pft_idx = pft_id  # PFT dimension index (bare soil=0, so PFT1=index 1)

            for a in range(NLON):
                for b in range(NLAT):
                    if htop[a, b, pft_idx, 0] > 0 and new_height[a, b] > 0:
                        # Compute ratio of bottom/top from original
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
        description="Aggregate 10m canopy height to CLM5 grid per PFT"
    )
    parser.add_argument('--canopy_height', required=True,
                        help='Path to 10m canopy height GeoTIFF')
    parser.add_argument('--pft_map', required=True,
                        help='Path to PFT classification GeoTIFF (500m)')
    parser.add_argument('--surfdata', required=True,
                        help='Path to CLM5 surface data NC file')
    parser.add_argument('--output_dir', required=True,
                        help='Output directory for modified surfdata files')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("Loading canopy height raster...")
    canopy_height, ch_h, ch_w, _ = load_raster(args.canopy_height)
    print(f"  Shape: {canopy_height.shape}")

    print("Loading PFT map...")
    pft_map, pft_h, pft_w, _ = load_raster(args.pft_map)
    print(f"  Shape: {pft_map.shape}")

    # Resize PFT map to match canopy height if needed
    if pft_map.shape != canopy_height.shape:
        from scipy.ndimage import zoom
        zoom_factors = (ch_h / pft_h, ch_w / pft_w)
        print(f"  Resampling PFT map with zoom factors: {zoom_factors}")
        pft_map = zoom(pft_map, zoom_factors, order=0)  # nearest neighbor

    methods = ['max', 'mean', 'median']

    for method in methods:
        print(f"\nAggregating with method: {method}")
        height_grids = {}

        for pft_id, pft_name in FOREST_PFTS.items():
            print(f"  PFT {pft_id} ({pft_name})...")
            grid, count = aggregate_to_clm_grid(canopy_height, pft_map, pft_id, method=method)
            n_valid = np.sum(count > 0)
            print(f"    Valid grid cells: {n_valid}, Mean height: {np.mean(grid[grid > 0]):.1f} m")
            height_grids[pft_id] = grid

        output_path = os.path.join(args.output_dir, f"{method.capitalize()}_surfdata.nc")
        print(f"  Updating surface data...")
        update_surfdata(args.surfdata, output_path, height_grids)

    print("\nDone! Output files:")
    for method in methods:
        print(f"  {os.path.join(args.output_dir, f'{method.capitalize()}_surfdata.nc')}")


if __name__ == '__main__':
    main()
