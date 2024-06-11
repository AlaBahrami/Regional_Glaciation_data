"""
Created on Thu May 30 14:52:22 2024

@author: AlaBahrami

Purpose of this script is to combine glacier region produced from Regional Glaciation 
Model for Western Canada, which are saved in .mat format and save data to netcdf 
and GeoTIFF formats by considering its original CRS projection information provided in readme file. 
As a test case example, I used CanESM-320km/RCP85 data; however any datasets produced by  
RGM, inclduing CSIRO-320km, GFDL-320km, HadGEM-320km, MIROC-320km, MPI-320km can be used 
and converted to a requested format. 
  

The time-mean ice thickness can be calculated for any time period of interest, 
for instance 2026 to 2055, and then save the averaged time-mean. Later, GeoTIFF
raster images can be used for calculating zonal statistics for any subbasin of 
interest covering western Canada. 
   
"""
#%% import modules 
import numpy as np
import os 
import h5py
import matplotlib.pyplot as plt
from scipy.io import loadmat
import netCDF4 as nc
import geopandas as gpd

#%% set directory of data 
directory        = 'Input/Glacier/CanESM-320km/RCP85/'
directory2       = 'Input/Glacier/Grids/All'
input_shape      = 'Input/Fraser/Basin_Boundary/Entire_Fraser_Basin.shp'

#%% list all mat files in the directory 
mat_files = [file for file in os.listdir(directory) if file.endswith('.mat')]
latlon_files = [file for file in os.listdir(directory2) if file.endswith('.mat')]

#%% functions 
def region_boundary(input_files):
    n = len(input_files)
    bound = np.zeros((n,4))
    for i in range(n):
        data = h5py.File(os.path.join(directory, input_files[i]))
        bound[i,0] = int(np.squeeze(data['map_western_edge']))
        bound[i,1] = int(np.squeeze(data['map_eastern_edge']))
        bound[i,2] = int(np.squeeze(data['map_northern_edge']))
        bound[i,3] = int(np.squeeze(data['map_southern_edge']))
    Xw = np.min(bound[:,0])
    Xe = np.max(bound[:,1])
    Yn = np.max(bound[:,2])
    Ys = np.min(bound[:,3])
    
    dx = int(np.squeeze(data['dx']))
    dy = int(np.squeeze(data['dy']))
    row = abs(int((Yn - Ys)/dy))
    col = abs(int((Xe - Xw)/dx))
    
    # boundary for each region
    index = np.zeros((n,4))  
    for i in range(n):
        index[i,0] = abs(int((Yn - bound[i,2])/dy))
        index[i,1] = abs(int((Yn - bound[i,3])/dy))
        index[i,2] = abs(int((Xw - bound[i,0])/dx))
        index[i,3] = abs(int((Xw - bound[i,1])/dx))
    index = index.astype(int)
    return bound, row, col, index

def region_mask(input_files, row, col, index):
    n = len(input_files)
    mask = np.zeros((row,col))
    for i in range(n):
        data = h5py.File(os.path.join(directory, input_files[i]))
        R_mask      = np.transpose(data['R_mask'])
        mask[index[i,0]:index[i,1],index[i,2]:index[i,3]] += R_mask
        R_mask[R_mask>1] = 1
    return mask

def ice_mask(input_files, row, col, index, t_st, t_fin):
    n = len(input_files)
    data = h5py.File(os.path.join(directory, input_files[0]))
    year = np.squeeze(data['year']).astype(int)
    rs = np.where(year == t_st)[0]
    rf = np.where(year == t_fin)[0]
    time_len = (rf - rs)[0] + 1  
    # surface array 
    S = np.zeros((time_len, row,col))
    # bed topography
    B = np.zeros((row,col))
    for i in range(n):
        data = h5py.File(os.path.join(directory, input_files[i]))
        B_region      = np.transpose(data['B'])
        B[index[i,0]:index[i,1],index[i,2]:index[i,3]] = B_region
        # surface topography
        S_region      = np.transpose(data['S'])
        S[:,index[i,0]:index[i,1],index[i,2]:index[i,3]] = S_region[rs[0]:(rf[0]+1),:,:]
    # average over time period
    S_mean    =  np.mean(S, axis = 0)
    H_ice = S_mean - B
    # set ice thickness less than 0.25 m to 0 and more than it to 1
    H_ice[H_ice < 0.25] = 0
    H_ice[H_ice > 0.25] = 1
    return H_ice

def lat_lon_extract(input_files, bound):
        # read data     
        data = loadmat(os.path.join(directory2, input_files[0]))
        lat_all      = data['latitude']
        lon_all      = data['longitude']
        # boundary merged regions 
        Xw = np.min(bound[:,0])
        Xe = np.max(bound[:,1])
        Yn = np.max(bound[:,2])
        Ys = np.min(bound[:,3])
        # pixel resolution 
        dx = int(np.squeeze(data['dx']))
        dy = int(np.squeeze(data['dy']))
        # boundary a large region containing the merged regions
        Xw_all = int(np.squeeze(data['map_western_edge']))
        Yn_all = int(np.squeeze(data['map_northern_edge']))
        # find out boundary for merged regions 
        row_s = abs(int((Yn_all - Yn)/dy))
        row_f = abs(int((Yn_all - Ys)/dy))
        col_s = abs(int((Xw_all - Xw)/dx))
        col_f = abs(int((Xw_all - Xe)/dx))
        # extract lat-lon for domain of interest
        lat = lat_all[row_s:row_f, col_s:col_f]
        lon = lon_all[row_s:row_f, col_s:col_f]
        return lat, lon, dx, dy

def netcdf_write(bound, dx, dy, longitude_array, latitude_array, input_array, output_nc_file):
    # domain boundary 
    map_western_edge = np.min(bound[:,0])
    map_eastern_edge = np.max(bound[:,1])
    map_northern_edge = np.max(bound[:,2])
    map_southern_edge = np.min(bound[:,3])
    
    # Generate coordinates for center of each array 
    x = np.arange(map_western_edge + 0.5 * dx, map_eastern_edge + 0.5 * dx, dx)
    y = np.arange(map_northern_edge - 0.5 * dy, map_southern_edge - 0.5 * dy, -dy)
    
    # Create NetCDF file
    ncfile = nc.Dataset(output_nc_file, mode='w', format='NETCDF4')
    
    # Create dimensions
    ncfile.createDimension('x', len(x))
    ncfile.createDimension('y', len(y))
    
    # Add global attribute
    ncfile.Conventions = "CF-1.5"
    
    # Create coordinate variables
    x_var = ncfile.createVariable('x', 'f8', ('x',))
    y_var = ncfile.createVariable('y', 'f8', ('y',))
    
    # Assign attributes to coordinate variables
    x_var.standard_name = "projection_x_coordinate"
    x_var.long_name = "x coordinate of projection"
    x_var.units = "m"
    x_var.axis = "X"
    
    y_var.standard_name = "projection_y_coordinate"
    y_var.long_name = "y coordinate of projection"
    y_var.units = "m"
    y_var.axis = "Y"
    
    # Write data to coordinate variables
    x_var[:] = x
    y_var[:] = y
    
    # Create the grid mapping variable (obtained data readme and data info)
    proj = ncfile.createVariable('lambert_conformal_conic', 'i4')
    proj.grid_mapping_name = "lambert_conformal_conic"
    proj.longitude_of_central_meridian = -107.0
    proj.false_easting = 0.0
    proj.false_northing = 0.0
    proj.latitude_of_projection_origin = 50.0
    proj.standard_parallel = [50.0, 50.0]
    proj.longitude_of_prime_meridian = 0.0
    proj.earth_radius = 6371200.0
    proj.longitudeOfFirstGridPointInDegrees = longitude_array[0,0]
    proj.latitudeOfFirstGridPointInDegrees = latitude_array[0,0]
    
    # Create the data variable
    band1 = ncfile.createVariable('Band1', 'f4', ('y', 'x'), fill_value=0.0)
    band1.long_name = "GDAL Band Number 1"
    band1.grid_mapping = "lambert_conformal_conic"
    band1.missing_value = 0.0
    
    # Write data to the data variable
    band1[:] = input_array
    
    # Close the NetCDF file
    ncfile.close()
    return    
#%% call function 
bound, row, col, index = region_boundary(mat_files)
R_mask = region_mask(mat_files, row, col, index)
# ice thinkess mask for baseline (2009)
H_ice_mask_2009 = ice_mask(mat_files, row, col, index, 2010, 2010)
# ice thinkess mask for 2026-2055 (2040)
H_ice_mask_2040 = ice_mask(mat_files, row, col, index, 2027, 2056)
# ice thinkess mask for 2026-2055 (2085)
H_ice_mask_2085 = ice_mask(mat_files, row, col, index, 2072, 2101)

lat,lon, dx, dy = lat_lon_extract(latlon_files, bound)

#%% show results with lat/lon with Fraser Basin boundary 
Fraser = gpd.read_file(input_shape)

lon_min = np.round(np.min(lon[:]))
lon_max = np.round(np.max(lon[:]))
lat_min = np.round(np.min(lat[:]))
lat_max = np.round(np.max(lat[:]))

fig, axs = plt.subplots(1 , 1, figsize=(20,20))
 # set axis 
axs.set_title('Ice mask time_mean 2071 - 2100 (2085)', fontsize=14, fontname='Times New Roman', fontweight='bold')
plt.setp(axs.get_xticklabels(), fontsize=14, fontname='Times New Roman', fontweight='bold')
plt.setp(axs.get_yticklabels(), fontsize=14, fontname='Times New Roman', fontweight='bold')
plt.setp(axs, xlim = [lon_min , lon_max])
plt.setp(axs, ylim = [lat_min , lat_max])
axs.set_xlabel('Longitude [\N{DEGREE SIGN}]', fontname = 'Times New Roman', 
                fontweight = 'bold', fontsize=14)

axs.set_ylabel('Latitude [\N{DEGREE SIGN}]', fontname = 'Times New Roman', 
                fontweight = 'bold', fontsize=14)

pcm  = axs.pcolormesh(lon, lat, H_ice_mask_2085, cmap='Greys', vmin = 0, vmax = 1) 
#cbar = fig.colorbar(pcm, ax = axs, location='right')

Fraser.plot(color='none', edgecolor='black',ax = axs)
plt.show() 
fs = "output/Ice_mask_2085_time_mean.png"
plt.savefig(fs, dpi=300)
plt.close()

#%% write data to netcdf
netcdf_write(bound, dx, dy, lon, lat, H_ice_mask_2009, 'output/Ice_mask_2009.nc')
netcdf_write(bound, dx, dy, lon, lat, H_ice_mask_2040, 'output/Ice_mask_2040.nc')
netcdf_write(bound, dx, dy, lon, lat, H_ice_mask_2085, 'output/Ice_mask_2085.nc')