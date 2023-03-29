#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
create_loc_daily_avg.py
    python script to read in hourly WRFChem data from Jian He daily files in UTC time,
    convert to local time, and create local daily average
    
    revised to work with new monthly output from Jian and to calcualte wet pm2.5 mass using
    dust emissons from the BAU case
Created on Tue Sep  7 16:37:50 2021
@author: kodell
"""
#%% user inputs
var_avg = 'pm25w' # needs to match variable name in the WRF Chem file, but pm needs to be pm25w
avg_type = '24hravg' # use 24hravg for pm calculation and mda8 for ozone
sim_type = 'BAU' # which WRF-Chem simulation to use
desc = 'full_annual_ts'

# indicate file locations
in_fp = '/Users/kodell/Desktop/GeoXO_modelling/Jian_data_20221011/outdir_12kmCONUS_'+sim_type+'/'
out_fp = '/Users/kodell/Library/CloudStorage/GoogleDrive-kodell@email.gwu.edu/My Drive/Ongoing Projects/GeoXO/concentration_data/Jian_WRF_processed_lda/'

#%% load modules
import netCDF4
from netCDF4 import Dataset
import numpy as np
from timezonefinder import TimezoneFinder
import pytz
import datetime as dt
import sys as sys
# functions I've written
from ODell_udf_Jian import calc_pmw

#%% pull lat, lon and averaging variable units
print('calculating local daily average for ',var_avg)

# load a lat lon and get units of variable to average
nc_fid = Dataset(in_fp+'wrfchem_data_2020_06_01_ts.nc')
lon = nc_fid['longitude'][:]
lat = nc_fid['latitude'][:]
if var_avg == 'pm25w': # calculate pm by hand
    var_units = 'ugm-3_STP'
else:
    var_units = nc_fid[var_avg].units
nc_fid.close()

#%% create grid of time zones to use to claculate utc offset for each day below (have to do for each day because of DST)
tz_grid = np.empty(lat.shape,dtype='<U40')
utc_diff = np.empty(lat.shape,dtype='float')
tf = TimezoneFinder()
for li in range(lon.shape[0]):
    for lj in range(lat.shape[1]):
        tz_str = tf.timezone_at(lng=lon[li,lj],lat=lat[li,lj])
        tz = pytz.timezone(tz_str)    
        tz_grid[li,lj] = tz_str

#%% loop through and load files, do calculations
# set dates to do calculation over, for the paper it is 04/01/2020 - 03/31/2021
start_date = dt.datetime(year=2020, month=4,day=1)
end_date = dt.datetime(year=2021, month=3,day=31)

# define arrays to fill
da_sfc_var = np.empty([(end_date-start_date).days,lon.shape[0],lon.shape[1]])
da_sfc_var[:] = np.nan
dates = ['9999-99-99']
offset_check = np.empty([(end_date-start_date).days,lon.shape[0],lon.shape[1]])
offset_check[:] = np.nan

# loop through days
ti = 0
today_dt = start_date
while today_dt < end_date:
    tmr_dt = today_dt + dt.timedelta(days=1)
    date_load_today = str(today_dt.year)+'_'+str(today_dt.month).zfill(2)+'_'+str(today_dt.day).zfill(2)
    date_load_tmr = str(tmr_dt.year)+'_'+str(tmr_dt.month).zfill(2)+'_'+str(tmr_dt.day).zfill(2)
        
    nc_fid_today = Dataset(in_fp+'wrfchem_data_'+ date_load_today+'_ts.nc')
    nc_fid_tmr = Dataset(in_fp+'wrfchem_data_'+ date_load_tmr+'_ts.nc')

    if var_avg == 'pm25w':
        if sim_type == 'BAU':
            today_var = calc_pmw(nc_fid_today,'p25') # calcualte wet pm2.5, but always use p25_bau per Jian instructions
            tmr_var = calc_pmw(nc_fid_tmr,'p25')
        else:
            today_var = calc_pmw(nc_fid_today,'p25_bau') # calcualte wet pm2.5, but always use p25_bau per Jian instructions
            tmr_var = calc_pmw(nc_fid_tmr,'p25_bau')
    else:
        today_var = nc_fid_today[var_avg][:]
        tmr_var  = nc_fid_tmr[var_avg][:]
    nc_fid_today.close()
    nc_fid_tmr.close()
    
    var_loc_time = np.empty(tmr_var.shape)
    var_loc_time[:] = np.nan
    # calculate var local time with extra hours for the mda8 ozone calculation
    var_loc_time_o3 = np.empty([31,tmr_var.shape[1],tmr_var.shape[2]])
    var_loc_time_o3[:] = np.nan
    hours = np.empty(tmr_var.shape)
    
    # create local time arrays
    for tz_str in np.unique(tz_grid):
        tz = pytz.timezone(tz_str)   
        offset = (tz.localize(today_dt)-pytz.timezone('utc').localize(today_dt)).seconds/3600.0
        inds = np.where(tz_grid==tz_str)
        var_loc_time[:,inds[0],inds[1]] = np.vstack([today_var[int(offset):,inds[0],inds[1]],
                           tmr_var[:int(offset),inds[0],inds[1]]])
        var_loc_time_o3[:,inds[0],inds[1]] = np.vstack([today_var[int(offset):,inds[0],inds[1]],
                           tmr_var[:int(offset)+7,inds[0],inds[1]]])
        # save offest hours to check our offset assignment to the grid
        hours[:,inds[0],inds[1]] = np.transpose(len(inds[0])*[np.arange(int(offset),int(offset)+24)])

    # calc average and append
    if avg_type == '24hravg':            
        da_sfc_var[ti,:,:] = np.mean(var_loc_time,axis=0)
    elif avg_type =='mda8':
        # calc MDA8 - maxiumum daily 8 hour average
        da8_all = np.empty(var_loc_time.shape)
        da8_all[:] = -999
        for hi in range(0,24):
            DA8 = np.mean(var_loc_time_o3[hi:hi+8],axis=0)
            da8_all[hi,:,:] = DA8
        mda8 = np.max(da8_all,axis=0)
        da_sfc_var[ti,:,:] = mda8
    else:
        print('averaging type not available')
        sys.exit()
    # offset check - assign the offset to the day so we can check DST works
    offset_check[ti,:,:] = hours[0,:,:]
    dates.append(date_load_today)
    today_dt = tmr_dt
    ti+=1
    print(date_load_today+' completed')
    
# create numpy arrays of the variables we created
dates = np.array(dates[1:]) # remove the initial 999 to start the array of dates
da_sfc_var = np.array(da_sfc_var)

#%% write to file
outfn = out_fp + 'WRFChem_' + sim_type + '_' + var_avg +'_' + avg_type + '_' +desc+'.nc'
nc_w_fid = netCDF4.Dataset(outfn, 'w', format='NETCDF4',clobber=True)
nc_w_fid.description = 'Daily average '+var_avg+'on local time created from WRFChem data from Jian He'

# Need to define dimensions that will be used in the file
nc_w_fid.createDimension('date',len(dates))
nc_w_fid.createDimension('lonx', lon.shape[0])
nc_w_fid.createDimension('lony', lon.shape[1])

# define variables for netcdf file
lon_nc = nc_w_fid.createVariable('lon', 'f8', ('lonx','lony'))
lat_nc = nc_w_fid.createVariable('lat', 'f8', ('lonx','lony'))
date_nc = nc_w_fid.createVariable('date','S10',('date'))
data_nc = nc_w_fid.createVariable(var_avg, 'f8', ('date', 'lonx','lony'))
data2_nc = nc_w_fid.createVariable('annual', 'f8', ('lonx','lony'))

lat_nc.setncatts({'units':'degrees','long_name':'Latitude (centers) [degrees]',\
               'var_desc':'degrees latitude for data grid centers?'})
lon_nc.setncatts({'units':'degrees','long_name':'Longitude (centers) [degrees]',\
               'var_desc':'degrees longitude for data grid centers?'})
date_nc.setncatts({'units':'days','long_name':'local date',\
               'var_desc':'date for local daily average'})
data_nc.setncatts({'units':var_units,'long_name':var_avg + ' daily average',\
               'var_desc':'daily average '+var_avg+' on local time'})
data2_nc.setncatts({'units':var_units,'long_name':var_avg + ' annual average',\
               'var_desc':'annual average '+var_avg+' on local time'})

# data
lat_nc[:,:] = lat
lon_nc[:,:] = lon
date_nc[:] = dates
data_nc[:,:,:] = da_sfc_var

nc_w_fid.close()
print('data saved to netCDF')











