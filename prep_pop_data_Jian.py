# -*- coding: utf-8 -*-
"""
prep_pop_data_Jian.py
    python script to regrid population data to the WRF Chem grid
    written by Katelyn O'Dell
    adapted from https://github.com/kaodell/smoke_HIA/blob/main/regrid_population.py
    02.03.22 - modified to regrid to Shobha's 4 km grid
    03.03.22 - modfied back to the 12 km grid
    03.08.22 - modified to make 2020 pop data for Molly Robertson on the kriging grid
    03.17.23 - finalized version for Jian GeoXO project, 12km WRF-chem grid
"""

#%% user inputs
# folder with population data and parent folder for outputs
prj_folder = '/Users/kodell/Library/CloudStorage/GoogleDrive-kodell@email.gwu.edu/My Drive/Ongoing Projects/GeoXO/'

# file with grid to regrid population data to, available from Dr. Jian He
to_grid_fn = '/Users/kodell/Desktop/GeoXO_modelling/Jian_data_20221011/outdir_12kmCONUS_COV/wrfchem_data_2020_07_ts.nc'

# NASA SECDAC population density file
# accessible here : https://sedac.ciesin.columbia.edu/data/collection/gpw-v4
pop_fn = 'population_data/NASA_SEDAC_pop/gpw-v4-population-density-rev11_totpop_2pt5_min_nc/gpw_v4_population_density_rev11_2pt5_min.nc'

# US national shapefile to sum US pop to check we are close and to make masks for plotting
# accessible here : https://www.census.gov/geographies/mapping-files/time-series/geo/carto-boundary-file.html
shp_fn = 'population_data/cb_2018_us_nation_5m/cb_2018_us_nation_5m'

# where to place regridded file
out_fn = 'population_data/NASA_SEDAC_pop/regridded/regrided_2020pop_12km_test.nc'

#%% load modules
import numpy as np
import netCDF4
import datetime
from mpl_toolkits.basemap import interp
import shapefile
from ODell_udf_Jian import plt_map, haversine
import matplotlib as mplt

#%% load data
# load  'to' grid, the 12 km WRF Chem grid from Jian
nc_fid = netCDF4.Dataset(to_grid_fn)
to_lat = nc_fid.variables['latitude'][:]
to_lon = nc_fid.variables['longitude'][:]
nc_fid.close()

# replace fill values with nans
to_glat = np.where(to_lat.data==to_lat.fill_value,np.nan,to_lat)
to_glon = np.where(to_lon.data==to_lon.fill_value,np.nan,to_lon)

# set grid area (have previously tested calculating grid area for other works
# as it is not exactly 12kmx12km, but it does not make a significant difference.)
grid_area = 12.0*12.0

# load population density file
# details of this dataset available here : https://sedac.ciesin.columbia.edu/data/collection/gpw-v4/documentation
nc_fid = netCDF4.Dataset(prj_folder+pop_fn) 
pop_density_data = nc_fid.variables['Population Density, v4.11 (2000, 2005, 2010, 2015, 2020): 2.5 arc-minutes'][:] # people per km2 
og_lon = nc_fid.variables['longitude'][:].data
og_lat = nc_fid.variables['latitude'][:].data
nc_fid.close()

# pull the data we need
pop_density_wm = pop_density_data[4] # 4 = 2020
land_area_wm = pop_density_data[8] + pop_density_data[9] # add water for pop counts to match the regridded calcs
national_ids = pop_density_data[10] # country IDs to pull US data for comparing totals

# replace fill values with zero
pop_density = np.ma.getdata(pop_density_wm)
land_area = np.ma.getdata(land_area_wm)
mask = np.ma.getmask(pop_density_wm)
la_mask = np.ma.getmask(land_area_wm)
pop_density[mask] = 0.0
land_area[la_mask] = np.nan

# cut to the US - use this instead of national ids to maintain grid structure
US_ind_latmax = np.max(np.where(og_lat > 20))
US_ind_latmin = np.min(np.where(og_lat < 60))
US_ind_lonmax = np.max(np.where(og_lon < -60))
US_ind_lonmin = np.min(np.where(og_lon > -130))
US_og_lon = np.copy(og_lon[US_ind_lonmin:US_ind_lonmax])
US_og_lat = np.copy(og_lat[US_ind_latmin:US_ind_latmax])
US_pop_density = np.copy(pop_density[ US_ind_latmin:US_ind_latmax,US_ind_lonmin:US_ind_lonmax])
US_land_area = np.copy(land_area[ US_ind_latmin:US_ind_latmax,US_ind_lonmin:US_ind_lonmax])
US_ids = np.copy(national_ids[ US_ind_latmin:US_ind_latmax,US_ind_lonmin:US_ind_lonmax])

#%% run grid interpolation
# lat has to be increasing for the function, so flip using pn.flipud
# function documentation: https://numpy.org/doc/stable/reference/generated/numpy.flipud.html
US_og_lat_flip = np.flipud(US_og_lat)
US_pop_density_flip = np.flipud(US_pop_density)

# use interp to regrid
# documentation: https://matplotlib.org/basemap/api/basemap_api.html
pop_density_regrid = interp(US_pop_density_flip, US_og_lon, US_og_lat_flip, to_glon, to_glat, order=1) 
# oredr = 1 is a bi-linear interpolation

#%% calculate grid area and population
# use grid area to get counts
pop_regrid = pop_density_regrid * grid_area

#%% make figures to compare
# first create meshgrid of lons and lats
US_og_lon_m, US_og_lat_m = np.meshgrid(US_og_lon,US_og_lat)
# plot original population density
# function plt_map is from ODell_udf_Jian and uses pcolormesh.
# refer to ODell_udf_jian.py for details
plt_map(US_og_lon_m, US_og_lat_m, US_pop_density,1, 'magma', 'population density', 
        'NASA SEDAC pop density 2020',clim=[0,100])
# regridded population density
plt_map(to_glon, to_glat, pop_density_regrid,1, 'magma', 'population density', 
        'regridded pop density 2020',clim=[0,100])

#%% print US pop totals with both grids to check
# load US shapefile
shps_file = shapefile.Reader(prj_folder+shp_fn)
shps_shp = shps_file.shape(0)
shps_records = shps_file.records()
shps_shapes = shps_file.shapes()

# create a US mask for the to-grid
si = 0
# designate array to fill
area_mask = np.zeros(to_glon.shape)
# create shp points of lat lon array
points = np.array((to_glon.flatten(), to_glat.flatten())).T
# loop through shape parts
for j in range(len(shps_records)):
        area_shp = shps_shapes[j]
        # identify inds for points for shape part
        for i in range(len(area_shp.parts)):
            i0 = area_shp.parts[i]
            if i < len(area_shp.parts)-1:
            		i1 = area_shp.parts[i+1] - 1
            else:
            		i1 = len(area_shp.points)
            seg = area_shp.points[i0:i1+1]
            # create path of shape part
            mpath = mplt.path.Path(seg)
            # create mask of points within shape part
            mask = mpath.contains_points(points).reshape(to_glon.shape)
            # fill in area mask for inds within the shape
            area_inds = np.where(mask==True)
            area_mask[area_inds] = 1
        
# check mask
plt_map(to_glon, to_glat, area_mask,1, 'jet', 'mask', 'US border mask',clim=[0,1])

# get US inds for NASA SECDAC population density grid
us_inds = np.where(US_ids == 840)
us_pop_og = US_land_area[us_inds[0],us_inds[1]]*US_pop_density[us_inds[0],us_inds[1]]

# total pop check
print('estimated contig. us pop regrid',np.sum(area_mask*pop_regrid))           
print('estimated contig. us pop SEDAC',np.sum(us_pop_og))           
print(100.0*(np.sum(area_mask*pop_regrid)-np.sum(us_pop_og))/np.sum(us_pop_og),'% difference')

#%% write to file
nc_w_fid = netCDF4.Dataset(prj_folder+out_fn, 'w', format='NETCDF4')
nc_w_fid.description = 'Regridded population from NASA SEDAC in 2020'
nc_w_fid.history = 'Created' + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# define file dimensions
nc_w_fid.createDimension('gridx', to_glon.shape[0])
nc_w_fid.createDimension('gridy', to_glon.shape[1])

lat_w = nc_w_fid.createVariable('lat', np.float32, ('gridx','gridy',))
lon_w = nc_w_fid.createVariable('lon', np.float32, ('gridx','gridy'))
pop_w = nc_w_fid.createVariable('population', np.float32, ('gridx','gridy',))
pop_density_w = nc_w_fid.createVariable('population_density', np.float32, ('gridx','gridy',))
grid_area_w = nc_w_fid.createVariable('grid_area', np.float32, ('gridx','gridy',))
us_mask_w = nc_w_fid.createVariable('us_mask', np.float32, ('gridx','gridy',))

lon_w[:] = to_lon
lat_w[:] = to_lat
pop_w[:] = pop_regrid
pop_density_w[:] = pop_density_regrid
grid_area_w[:] = grid_area
us_mask_w[:] = area_mask

nc_w_fid.close()

print('regridded population saved to netCDF')

