from netCDF4 import Dataset
import numpy as np
from rasterio.warp import calculate_default_transform, reproject, Resampling
from netCDF4 import num2date, date2num
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
from skimage.transform import resize
## scikit-image
from rasterio.transform import Affine
from skimage.measure import block_reduce
import matplotlib.pyplot as plt
import pandas as pd
from rasterio import features
import tarfile
import random
import rasterio
from rasterio.plot import show
from rasterio.merge import merge
from geopandas.tools import sjoin
from geopandas import GeoDataFrame
import geopandas as gpd
from shapely.geometry import Point
from scipy import interpolate
import elevation
import richdem as rd
from scipy import interpolate
import os
# from shapely import speedups
# speedups.disable()
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
# import cartopy.io.shapereader as shpreader
# from osgeo import gdal
# from pyhdf.SD import SD, SDC
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
# import cartopy.io.shapereader as shpreader
import matplotlib
import matplotlib.ticker as mticker
# from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
def CLM_Surface_data_generation():
    source_dir=r'H:\CLM input\source/'
    target_dir=r'H:\CLM input\target/'
    source_file_name = source_dir+'surfdata_1.9x2.5_hist_16pfts_Irrig_CMIP6_simyr2000_c190304.nc'
    target_file_name = target_dir + 'surfdata_1.9x2.5_hist_16pfts_Irrig_CMIP6_simyr2000_c190304.nc'
    src = Dataset(source_file_name, 'r')
    dst=Dataset(target_file_name, "w")
    # variable_list=nc.variables.keys()
    # nc.close()
    # print(variable_list)

    # with Dataset(source_file_name) as src, Dataset(target_file_name, "w") as dst:
    # copy global attributes all at once via dictionary
    dst.setncatts(src.__dict__)
    # copy dimensions
    for name, dimension in src.dimensions.items():
        dst.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))
    # copy all file data except for the excluded
    for name, variable in src.variables.items():
        x = dst.createVariable(name, variable.datatype, variable.dimensions)
        dst[name][:] = src[name][:]
        # copy variable attributes all at once via dictionary
        dst[name].setncatts(src[name].__dict__)


    CTH = src['MONTHLY_LAI'][:]/100
    src.close()
    CTH_var = dst.createVariable('Newtop', "f8", ("time","lsmpft","lsmlat", "lsmlon",))
    CTH_var.units = 'meter'
    CTH_var.long_name = 'GEDI canopy top height'
    CTH_var[:] = CTH[:]
    dst.close()
    print('CLM_Surface_data_generation done!')

def CTH_nc_data_generation():
    pft_file = 'H:/CLM input/PFT_500m_2020.tif'
    src_pft = rasterio.open(pft_file, "r")
    PFT_data = src_pft.read(1)  # .astype(np.float32)

    CTH_file = 'H:/ETH_GlobalCanopyHeight_10m_2020_version1/ETH_GlobalCanopyHeight_10m_2020_version1/Global Canopy Height/cth2020_GEE1.tif'
    # ci_file = f'F:/lifa/uwm/clumpingIndex/ShanshanWei/2006/Full_Extend_Global_MODIS_Clumping_Index_A2006001_500m_V001_8days_Comp_Postprocess_Complete.tif'
    # ci_file = f'F:/lifa/uwm/clumpingIndex/LimingHe/Global_Clumping_Index_1531/data/global_clumping_index.tif'
    src_CTH = rasterio.open(CTH_file, "r")
    CTH_data = src_CTH.read(1)
    # CI_data[CI_data==255]=0
    # show(CI_data)

    nodataval=np.nan
    pft_num=16
    degree = 2  # 2 for 1.9x2.5
    ###for 1 degree
    if degree == 1:
        delta_lat = 240
        delta_lon = 240
        nlats = 180
        nlons = 360
    elif degree == 0.5:
        delta_lat = 120
        delta_lon = 120
        nlats = 360
        nlons = 720
    elif degree == 0.05:
        delta_lat = 12
        delta_lon = 12
        nlats = 3600
        nlons = 7200
    ###for 1.9x2.5 degree
    elif degree == 2:
        delta_lat = 450
        delta_lon = 600
        nlats = 96
        nlons = 144


    CTH_PFT = np.full((pft_num, nlats, nlons), nodataval, dtype='float32')
    for lat_idx in range(0, src_pft.height, delta_lat):
        for lon_idx in range(0, src_pft.width, delta_lon):
            temp_CTH = CTH_data[lat_idx:(lat_idx + delta_lat), lon_idx:(lon_idx + delta_lon)]
            temp_PFT = PFT_data[lat_idx:(lat_idx + delta_lat), lon_idx:(lon_idx + delta_lon)]
            for pft_idx in range(0, pft_num):
                temp_CTH_Pft = temp_CTH[((temp_PFT==pft_idx)&(temp_CTH<255)&(temp_CTH>0))].reshape(-1)
                if len(temp_CTH_Pft) > 0:
                    temp_CTH_Pft = np.mean(temp_CTH_Pft)
                else:
                    temp_CTH_Pft = 0
                CTH_PFT[pft_idx, int(lat_idx / delta_lat), int(lon_idx / delta_lon)] = temp_CTH_Pft
    print('stats completed!')
    lats=np.linspace(90,-90,CTH_PFT.shape[1])
    lons = np.linspace(-180, 180, CTH_PFT.shape[2])
    dir = r'H:\new_canopy_top_height_nc\nc/'
    save_file = dir + 'CTH_2020-221209.nc'
    nc_fid2 = Dataset(save_file, 'w', format="NETCDF4")
    nc_fid2.createDimension('lat', len(lats))
    nc_fid2.createDimension('lon', len(lons))
    nc_fid2.createDimension('pft', CTH_PFT.shape[0])

    latitudes = nc_fid2.createVariable('lat', 'f4', ('lat',))
    longitudes = nc_fid2.createVariable('lon', 'f4', ('lon',))
    type_var = nc_fid2.createVariable("pft", "f4", ("pft",))
    CTH_var = nc_fid2.createVariable('CTH', "f4", ("pft", "lat", "lon",), zlib=True)
    type_var[:] = range(pft_num)
    # time_day.units = "days since 1993-08-01 00:00:00.0"
    # BA_var.units = 'acres'
    # time_day.calendar = "gregorian"
    # dates = [datetime(2006, 1, 1) + relativedelta(years=+n) for n in range(Fire_data.shape[0])]
    # time_day[:] = date2num(dates, units=time_day.units, calendar=time_day.calendar)
    latitudes[:] = lats[:]
    longitudes[:] = lons[:]
    CTH_var[:] = CTH_PFT[:]
    nc_fid2.close()
    print('CTH_nc_data_generation done!')

def upscale_CTH_nc_data_generation():
    pft_file = 'H:/CLM input/PFT_500m_2020.tif'
    src_pft = rasterio.open(pft_file, "r")

    cth_file = 'H:/ETH_GlobalCanopyHeight_10m_2020_version1/ETH_GlobalCanopyHeight_10m_2020_version1/Global Canopy Height/cth2020_GEE1.tif'
    src_cth = rasterio.open(cth_file, "r")
    CTH_data = src_cth.read(1)
    # CI_data[CI_data==255]=0
    # show(CI_data)

    nodataval=np.nan
    degree = 2 # 2 for 1.9x2.5
    ###for 1 degree
    if degree == 1:
        delta_lat = 240
        delta_lon = 240
        nlats = 180
        nlons = 360
    elif degree==0.5:
        delta_lat = 120
        delta_lon = 120
        nlats = 360
        nlons = 720
    elif degree==0.05:
        delta_lat = 12
        delta_lon = 12
        nlats = 3600
        nlons = 7200
    ###for 1.9x2.5 degree
    elif degree==2:
        delta_lat = 450
        delta_lon = 600
        nlats = 96
        nlons = 144


    CTH_PFT = np.full((nlats, nlons), nodataval, dtype='float32')
    for lat_idx in range(0, src_pft.height, delta_lat):
        for lon_idx in range(0, src_pft.width, delta_lon):
            temp_CTH = CTH_data[lat_idx:(lat_idx + delta_lat), lon_idx:(lon_idx + delta_lon)]
            temp_CTH = temp_CTH[(temp_CTH>0)&(temp_CTH<255)].reshape(-1)
            if len(temp_CTH) > 0:
                temp_CTH = np.mean(temp_CTH)
            else:
                temp_CTH = np.nan
            CTH_PFT[int(lat_idx / delta_lat), int(lon_idx / delta_lon)] = temp_CTH
    print('stats completed!')
    lats=np.linspace(90,-90,CTH_PFT.shape[0])
    lons = np.linspace(-180, 180, CTH_PFT.shape[1])
    dir = r'H:\new_canopy_top_height_nc\upscale_nc/'
    save_file = dir + 'CTH_2020-221209.nc'
    nc_fid2 = Dataset(save_file, 'w', format="NETCDF4")
    nc_fid2.createDimension('lat', len(lats))
    nc_fid2.createDimension('lon', len(lons))
    latitudes = nc_fid2.createVariable('lat', 'f4', ('lat',))
    longitudes = nc_fid2.createVariable('lon', 'f4', ('lon',))
    CTH_var = nc_fid2.createVariable('CTH', "f4", ("lat", "lon",), zlib=True)
    # time_day.units = "days since 1993-08-01 00:00:00.0"
    # BA_var.units = 'acres'
    # time_day.calendar = "gregorian"
    # dates = [datetime(2006, 1, 1) + relativedelta(years=+n) for n in range(Fire_data.shape[0])]
    # time_day[:] = date2num(dates, units=time_day.units, calendar=time_day.calendar)
    latitudes[:] = lats[:]
    longitudes[:] = lons[:]
    CTH_var[:] = CTH_PFT[:]
    nc_fid2.close()
    print('upscale_CTH_nc_data_generation done!')






def upscale_PFT_nc_data_generation():
    pft_file = 'H:/CLM input/PFT_500m_2020.tif'
    src_pft = rasterio.open(pft_file, "r")
    pft_data=src_pft.read(1)

    nodataval=np.nan
    degree = 2 # 2 for 1.9x2.5
    ###for 1 degree
    if degree == 1:
        delta_lat = 240
        delta_lon = 240
        nlats = 180
        nlons = 360
    elif degree==0.5:
        delta_lat = 120
        delta_lon = 120
        nlats = 360
        nlons = 720
    elif degree==0.05:
        delta_lat = 12
        delta_lon = 12
        nlats = 3600
        nlons = 7200
    ###for 1.9x2.5 degree
    elif degree==2:
        delta_lat = 450
        delta_lon = 600
        nlats = 96
        nlons = 144

    pft_num=16
    Upscaled_PFT = np.full((pft_num,nlats, nlons), nodataval, dtype='float32')
    for lat_idx in range(0, src_pft.height, delta_lat):
        for lon_idx in range(0, src_pft.width, delta_lon):
            temp_CTH = pft_data[lat_idx:(lat_idx + delta_lat), lon_idx:(lon_idx + delta_lon)]
            for pft_idx in range(pft_num):
                temp_pft = temp_CTH[(temp_CTH==pft_idx)].reshape(-1)
                temp_pft = len(temp_pft)/(delta_lat*delta_lon)
                Upscaled_PFT[pft_idx,int(lat_idx / delta_lat), int(lon_idx / delta_lon)] = temp_pft
    print('stats completed!')
    lats=np.linspace(90,-90,Upscaled_PFT.shape[1])
    lons = np.linspace(-180, 180, Upscaled_PFT.shape[2])
    dir = r'H:\new_canopy_top_height_nc\PFT/'
    save_file = dir + 'PFT_2020-{degree}degree221209.nc'
    nc_fid2 = Dataset(save_file, 'w', format="NETCDF4")
    nc_fid2.createDimension('lat', len(lats))
    nc_fid2.createDimension('lon', len(lons))
    nc_fid2.createDimension('pft', pft_num)
    latitudes = nc_fid2.createVariable('lat', 'f4', ('lat',))
    longitudes = nc_fid2.createVariable('lon', 'f4', ('lon',))
    pfts = nc_fid2.createVariable('pft', 'f4', ('pft',))
    pft_var = nc_fid2.createVariable('PFT_Percentage', "f4", ("pft","lat", "lon",), zlib=True)
    # time_day.units = "days since 1993-08-01 00:00:00.0"
    # BA_var.units = 'acres'
    # time_day.calendar = "gregorian"
    # dates = [datetime(2006, 1, 1) + relativedelta(years=+n) for n in range(Fire_data.shape[0])]
    # time_day[:] = date2num(dates, units=time_day.units, calendar=time_day.calendar)
    latitudes[:] = lats[:]
    longitudes[:] = lons[:]
    pfts[:]=range(pft_num)
    pft_var[:] = Upscaled_PFT[:]
    nc_fid2.close()
    print('upscale_PFT_nc_data_generation done!')

# def CI_Ziti_data_reproject():
#     dir = 'G:/UWM/ClumpingIndex/ZitiJiao/global_tif/sin/'
#     files_all=os.listdir(dir)
#     dst_crs="EPSG:4326"
#     pft_file = 'F:/lifa/uwm/PFT_MODIS/tif/wgs84/2006.tif'
#     # pft_file = f'F:/lifa/uwm/PFT_MODIS/tif/sin_projection/{year}.tif'
#     src_pft = rasterio.open(pft_file, "r")
#     for tempfile in files_all[:]:
#         tempfile_dir = dir + tempfile
#         src = rasterio.open(tempfile_dir, "r")
#
#         out_meta=src.meta.copy()
#         out_meta.update({"driver": "GTiff",
#                          "height": src_pft.height,
#                          "width": src_pft.width,
#                          "transform": src_pft.transform,  #
#                          "count": 1,
#                          'compress': 'lzw',
#                          "crs": dst_crs,  # "+proj=utm +zone=35 +ellps=GRS80 +units=m +no_defs"
#                          # 'nodata': nodataval,
#                          })  #
#         out_fp = 'G:/UWM/ClumpingIndex/ZitiJiao/global_tif/sin/wgs84/{tempfile}'
#         # Write the mosaic raster to disk
#         with rasterio.open(out_fp, 'w', **out_meta) as dst:
#             for i in range(1, src.count + 1):
#                 reproject(
#                     source=rasterio.band(src, i),
#                     destination=rasterio.band(dst, i),
#                     src_transform=src.transform,
#                     src_crs=src.crs,
#                     dst_transform=src_pft.transform,
#                     dst_crs=dst_crs,
#                     resampling=Resampling.nearest)
#         print(tempfile, 'reprojection done!')

def CLM_surface_data_with_CTH_generation():

    raw_CTH_file=r'H:\new_canopy_top_height_nc\nc\CTH_2020-221209.nc'
    src = Dataset(raw_CTH_file, 'r')
    CTH_data=src.variables['CTH'][:]
    src.close()
    CTH_data=CTH_data*1
    CTH_data[CTH_data==0]=1###########assign pixels without specific pft to 1


    CTH_data=CTH_data[np.newaxis,:]
    CTH_data=np.repeat(CTH_data,12,axis=0)
    last_pft_cth=np.ones((12,1,CTH_data.shape[2],CTH_data.shape[3]))
    CTH_data=np.concatenate((CTH_data,last_pft_cth),axis=1)
    CTH_data_new=np.full(CTH_data.shape,np.nan)
    half_width=int(CTH_data.shape[3]/2)
    # CI_data_new[:]=CI_data[:]

    CTH_data_new[:,:,:,:half_width]=CTH_data[:,:,:,half_width:]
    CTH_data_new[:, :, :, half_width:] = CTH_data[:, :, :, :half_width]
    CTH_data_new[:] = CTH_data_new[:, :, ::-1]

    #####################################
    #####generate surface data with CTH
    source_dir = r'H:\CLM input\source/'
    target_dir = r'H:\CLM input\target/surface_data/'
    source_file_name = source_dir + 'surfdata_1.9x2.5_hist_16pfts_Irrig_CMIP6_simyr2000_c190304.nc'
    target_file_name = target_dir + 'surfdata_1.9x2.5_hist_16pfts_Irrig_CMIP6_simyr2000_c190304_new_canopy_top_at_2020_20221209.nc'
    src = Dataset(source_file_name, 'r')
    dst = Dataset(target_file_name, "w")
    # variable_list=nc.variables.keys()
    # nc.close()
    # print(variable_list)

    # with Dataset(source_file_name) as src, Dataset(target_file_name, "w") as dst:
    # copy global attributes all at once via dictionary
    dst.setncatts(src.__dict__)
    # copy dimensions
    for name, dimension in src.dimensions.items():
        dst.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))
    # copy all file data except for the excluded
    for name, variable in src.variables.items():
        x = dst.createVariable(name, variable.datatype, variable.dimensions)
        dst[name][:] = src[name][:]
        # copy variable attributes all at once via dictionary
        dst[name].setncatts(src[name].__dict__)
    src.close()
    CTH_var = dst.createVariable('NEW_Canopy_height', "f4", ("time", "lsmpft", "lsmlat", "lsmlon",))
    CTH_var.units = 'meters'
    CTH_var.long_name = '2020 monthly height top'
    CTH_var[:] = CTH_data_new[:]
    dst.close()
    print('CLM_surface_data_with_CTH_generation done!')

# def interpolate_upscale_CI_nc_data_ziti():
#     dir = r'G:\UWM\ClumpingIndex\ZitiJiao\nc/'
#     degree = 0.5
#     nc_file = dir + 'CI_2006-{degree}degree.nc'
#     src = Dataset(nc_file, 'r')
#     CI_data = src.variables['CI'][:]
#     src.close()
#     nan_num_max = 6
#     for lat_idx in range(CI_data.shape[1]):
#         for lon_idx in range(CI_data.shape[2]):
#             nan_num = np.sum(np.isnan(CI_data[:, lat_idx, lon_idx]))
#             if nan_num > 0 and nan_num < nan_num_max:
#                 cell_data = CI_data[:, lat_idx, lon_idx]
#                 nonan_index = np.where(~np.isnan(cell_data))[0]
#                 # print(mask[i,j])
#                 x = nonan_index
#                 y = cell_data[nonan_index]
#                 f = interpolate.interp1d(x, y, kind='nearest', fill_value="extrapolate")  # linear,nearest
#                 x_new = np.where(np.isnan(cell_data))[0]
#                 # print(x_new)
#                 y_new = f(x_new)
#                 CI_data[x_new, lat_idx, lon_idx] = y_new
#             elif nan_num != 0:
#                 CI_data[:, lat_idx, lon_idx] = np.nan
#     lats = np.linspace(90, -90, CI_data.shape[1])
#     lons = np.linspace(-180, 180, CI_data.shape[2])
#     dir = r'G:\UWM\ClumpingIndex\ZitiJiao\nc/'
#     save_file = dir + 'CI_2006-{degree}degree_interpolated.nc'
#     nc_fid2 = Dataset(save_file, 'w', format="NETCDF4")
#     nc_fid2.createDimension('lat', len(lats))
#     nc_fid2.createDimension('lon', len(lons))
#     nc_fid2.createDimension('time', CI_data.shape[0])
#     latitudes = nc_fid2.createVariable('lat', 'f4', ('lat',))
#     longitudes = nc_fid2.createVariable('lon', 'f4', ('lon',))
#     time_day = nc_fid2.createVariable("time", "f4", ("time",))
#     CI_var = nc_fid2.createVariable('CI', "f4", ("time", "lat", "lon",), zlib=True)
#     time_day.units = "days since 1993-08-01 00:00:00.0"
#     time_day.calendar = "gregorian"
#     # dates = [datetime(2006, 1, 1) + relativedelta(days=+n*8) for n in range(CI_PFT.shape[0])]
#     dates = [datetime(2006, 1, 1) + relativedelta(months=+n) for n in range(CI_data.shape[0])]
#     time_day[:] = date2num(dates, units=time_day.units, calendar=time_day.calendar)
#     latitudes[:] = lats[:]
#     longitudes[:] = lons[:]
#     CI_var[:] = CI_data[:]
#     nc_fid2.close()
#     print('done!')

# def GPP_nc_data_generation_GOSIF():
#     gpp_file = 'G:/UWM/ClumpingIndex/GPP_data/GOSIF/tiff/annual/GOSIF_GPP_2006_Mean.tif'
#     src_pft = rasterio.open(gpp_file, "r")
#     gpp_data = src_pft.read(1)  # .astype(np.float32)
#
#     scale_factor=0.1
#     nodataval=np.nan
#     degree = 2 # 2 for 1.9x2.5
#     ###for 1 degree
#     if degree == 1:
#         delta_lat = 240
#         delta_lon = 240
#         nlats = 180
#         nlons = 360
#     ###for 1.9x2.5 degree
#     else:
#         delta_lat = 36
#         delta_lon = 50
#         nlats = 100
#         nlons = 144
#
#
#     CI_PFT = np.full((nlats, nlons), nodataval, dtype='float32')
#     for lat_idx in range(0, src_pft.height, delta_lat):
#         for lon_idx in range(0, src_pft.width, delta_lon):
#             temp_gpp = gpp_data[lat_idx:(lat_idx + delta_lat), lon_idx:(lon_idx + delta_lon)]
#             # temp_CI_Pft = temp_gpp[((temp_gpp>=1)&(temp_gpp<255))].reshape(-1)
#             temp_CI_Pft = temp_gpp[((temp_gpp < 50000))].reshape(-1)
#             if len(temp_CI_Pft) > 0:
#                 temp_CI_Pft = np.mean(temp_CI_Pft)
#             else:
#                 temp_CI_Pft = 0
#             CI_PFT[int(lat_idx / delta_lat), int(lon_idx / delta_lon)] = temp_CI_Pft
#     CI_PFT_resized=resize(CI_PFT, (96 , 144),anti_aliasing=False, order=1, mode='edge')*scale_factor
#     CI_PFT_new=np.full(CI_PFT_resized.shape,np.nan)
#     CI_PFT_new[:,:(int(CI_PFT_resized.shape[1]/2))]=CI_PFT_resized[:,(int(CI_PFT_resized.shape[1]/2)):]
#     CI_PFT_new[:, (int(CI_PFT_resized.shape[1] / 2)):] = CI_PFT_resized[:, :(int(CI_PFT_resized.shape[1] / 2))]
#     CI_PFT_new = CI_PFT_new[::-1]
#     CI_PFT_resized=CI_PFT_new
#     print('stats completed!')
#     lats=np.linspace(-90,90,CI_PFT_resized.shape[0])
#     lons = np.linspace(0, 360, CI_PFT_resized.shape[1])
#     dir = r'G:\UWM\ClumpingIndex\GPP_data\GOSIF\tiff\annual/'
#     save_file = dir + 'GPP_2006-{degree}degree_T62.nc'
#     nc_fid2 = Dataset(save_file, 'w', format="NETCDF4")
#     nc_fid2.createDimension('lat', len(lats))
#     nc_fid2.createDimension('lon', len(lons))
#     # nc_fid2.createDimension('pft', CI_PFT.shape[0])
#     latitudes = nc_fid2.createVariable('lat', 'f4', ('lat',))
#     longitudes = nc_fid2.createVariable('lon', 'f4', ('lon',))
#     # type_var = nc_fid2.createVariable("pft", "f4", ("pft",))
#     # CI_var = nc_fid2.createVariable('CI', "f4", ("pft", "lat", "lon",), zlib=True)
#     CI_var = nc_fid2.createVariable('GPP', "f4", ( "lat", "lon",), zlib=True)
#     # type_var[:] = range(pft_num)
#     # time_day.units = "days since 1993-08-01 00:00:00.0"
#     CI_var.units = 'gC m-2 year-1'
#     # time_day.calendar = "gregorian"
#     # dates = [datetime(2006, 1, 1) + relativedelta(years=+n) for n in range(Fire_data.shape[0])]
#     # time_day[:] = date2num(dates, units=time_day.units, calendar=time_day.calendar)
#     latitudes[:] = lats[:]
#     longitudes[:] = lons[:]
#     CI_var[:] = CI_PFT_resized[:]
#     nc_fid2.close()
#     print('done!')







if __name__=='__main__':

    CLM_Surface_data_generation()

    ####################generate pft-based canopy top height from tif
    CTH_nc_data_generation()

    upscale_CTH_nc_data_generation()

    ####################upscale canopy top hight from tif for data uncertainty analysis


    upscale_PFT_nc_data_generation()

    CLM_surface_data_with_CTH_generation()

    print('ALL steps finished!')



