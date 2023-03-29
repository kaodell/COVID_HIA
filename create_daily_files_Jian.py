#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
create_daily_files_Jian.py
    python script to transfer Jian's monthly files to daily files so they can work with my previous code
Created on Fri Oct 28 11:22:20 2022
written by Katelyn O'Dell
"""
#%% import modules
from netCDF4 import Dataset
import netCDF4
import numpy as np
import sys

#%% user inputs
# filenames to load for this analysis, note filenames for COV and BAU files are
# the same but they are stored in different folders.
# files span the months 04/2020 to 03/2021
# These data can be obtained through correspondance with Dr. Jian He
file_names_2load = ['wrfchem_data_2020_04_ts.nc',
                    'wrfchem_data_2020_05_ts.nc',
                    'wrfchem_data_2020_06_ts.nc',
                    'wrfchem_data_2020_07_ts.nc',
                    'wrfchem_data_2020_08_ts.nc',
                    'wrfchem_data_2020_09_ts.nc',
                    'wrfchem_data_2020_10_ts.nc',
                    'wrfchem_data_2020_11_ts.nc',
                    'wrfchem_data_2020_12_ts.nc',
                    'wrfchem_data_2021_01_ts.nc',
                    'wrfchem_data_2021_02_ts.nc',
                    'wrfchem_data_2021_03_ts.nc']
# indicate if you wish to run code for the COV of BAU simulations (have to run code sprately for both)
version = 'COV'
# path to data
data_path = '/Users/kodell/Desktop/GeoXO_modelling/Jian_data_20221011/outdir_12kmCONUS_'+version+'/'

#%% loop through files and create daily files
diff = [] # check to make sure we are assigning the days we expect for each daily file
for fn in file_names_2load:
    file_load = data_path + fn
    print('loading ',fn)
    # load data
    nc_fid = Dataset(file_load)
    lat = nc_fid['latitude'][:]
    lon = nc_fid['longitude'][:]
    hour = nc_fid['Time'][:]
    #print(nc_fid.summary)

    day_num = 1
    for i in np.arange(0,len(hour),24):
        # create netCDF file for the day 
        outfn = file_load.split('.')[0][:-2] + str(day_num).zfill(2) + '_ts.nc'
        nc_w_fid = netCDF4.Dataset(outfn, 'w', format='NETCDF4',clobber=True)
        nc_w_fid.description = 'Hourly WRFChem data from Jian He, '+version
        
        # Need to define dimensions that will be used in the file
        nc_w_fid.createDimension('hour',24)
        nc_w_fid.createDimension('lonx', lon.shape[0])
        nc_w_fid.createDimension('lony', lon.shape[1])
    
        lon_nc = nc_w_fid.createVariable('longitude', 'f8', ('lonx','lony'))
        lat_nc = nc_w_fid.createVariable('latitude', 'f8', ('lonx','lony'))
        hour_nc = nc_w_fid.createVariable('hour', 'f8', ('hour'))
    
        lat_nc.setncatts({'units':'degrees','long_name':'Latitude (centers) [degrees]',\
                       'var_desc':'degrees latitude for data grid centers'})
        lon_nc.setncatts({'units':'degrees','long_name':'Longitude (centers) [degrees]',\
                       'var_desc':'degrees longitude for data grid centers'})
        hour_nc.setncatts({'units':'hours','long_name':'hour of day UTC',\
                       'var_desc':'hour of day'})
        lat_nc[:,:] = lat
        lon_nc[:,:] = lon
        hour_nc[:] = hour[i:i+24]
        
        # now pull pollution data for the day
        for var_str in nc_fid.variables.keys():
            if var_str in ['longitude', 'latitude', 'Time','Times']:
                continue
            # create variable for new file
            var_nc = nc_w_fid.createVariable(var_str, 'f8', ('hour', 'lonx','lony'))
            # pull descriptors from monthly file
            var_nc.setncatts({'units':nc_fid[var_str].units,
                           'var_desc':nc_fid[var_str].description})    
            # add variable information for that day to the file
            var_nc[:,:,:] = nc_fid[var_str][:][i:i+24,:,:]
        nc_w_fid.close()

        print(day_num,'saved')
        # check that code is indexing properly
        #print(i)        
        #print(hour[i:i+24][0],hour[i:i+24][-1])
        
        # check that the file we just saved is the day we expected
        day_file = int(np.unique(nc_fid['Times'][:][i:i+24,8])[0]+np.unique(nc_fid['Times'][:][i:i+24,9])[0])
        if (day_file - day_num) != 0.0: # these should always be the same
            sys.exit('day mistatch')
        day_num += 1


