#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mk_figures_4Jian.py
    python script to make figures for Jian He's 2023 manuscript
    code uses output run_HIA_gridded.py, which runs an HIA on the 
    WRF Chem output from Jian
Created on Tue Jul 19 10:25:40 2022
09.22.22
    updated to add mortality calculations
written by: Katelyn O'Dell
"""

#%% user inputs
prj_folder = '/Users/kodell/Library/CloudStorage/GoogleDrive-kodell@email.gwu.edu/My Drive/Ongoing Projects/GeoXO/'
data_folder = prj_folder + 'HIA_results/Jian/'
figure_folder = prj_folder + 'figures/Jian/'
pop_fn = 'population_data/NASA_SEDAC_pop/regridded/regrided_2020pop_12km_final.nc'
shp_fn = prj_folder + 'population_data/cb_2018_us_nation_5m/'

#%% import modules
import numpy as np
from netCDF4 import Dataset
from ODell_udf_Jian import plt_map
import shapefile
import matplotlib as mplt

#%% load data
# mortality pm
fid = np.load(data_folder + 'pm25w_mort_hia_full_cntybr_final.npz')
pm_bau_mort_tot = fid['bau_mort_tot'][:]
pm_cov_mort_tot = fid['cov_mort_tot'][:]
pm_poll_bau_avg = fid['poll_bau_avg'][:]
pm_poll_cov_avg = fid['poll_cov_avg'][:]
fid.close()

# mortality 03
fid = np.load(data_folder + 'o3_mort_hia_full_cntybr_final.npz')
o3_bau_mort_tot = fid['bau_mort_tot'][:]
o3_cov_mort_tot = fid['cov_mort_tot'][:]
o3_poll_bau_avg = fid['poll_bau_avg'][:]
o3_poll_cov_avg = fid['poll_cov_avg'][:]
fid.close()

# population
nc_fid = Dataset(prj_folder+pop_fn)
lat = nc_fid['lat'][:]
lon = nc_fid['lon'][:]
pop = nc_fid['population'][:].data
nc_fid.close()

# shapefile for cutting to the US
shps_file = shapefile.Reader(shp_fn+'cb_2018_us_nation_5m')
shps_shp = shps_file.shape(0)
shps_records = shps_file.records()
shps_shapes = shps_file.shapes()

#%% make US mask for figures
si = 0
area_mask = np.zeros(lon.shape)
area_mask[:] = np.nan
for j in range(len(shps_records)):
        area_shp = shps_shapes[j]
        for i in range(len(area_shp.parts)):
            i0 = area_shp.parts[i]
            if i < len(area_shp.parts)-1:
            		i1 = area_shp.parts[i+1] - 1
            else:
            		i1 = len(area_shp.points)
            seg = area_shp.points[i0:i1+1]
            mpath = mplt.path.Path(seg)
            points = np.array((lon.flatten(), lat.flatten())).T
            mask = mpath.contains_points(points).reshape(lon.shape)
            area_inds = np.where(mask==True)
            area_mask[area_inds] = 1
pop = pop*area_mask
#plt_map(lon,lat,area_mask,1,'magma','area flag','area mask check',clim=[0,1])

#%% Mortality map for manuscript
# calculate mort per 10^5 people
pm_bau_mort_tot_pp = area_mask*(10**5.0)*pm_bau_mort_tot[0,:,:]/pop
pm_cov_mort_tot_pp = area_mask*(10**5.0)*pm_cov_mort_tot[0,:,:]/pop
o3_bau_mort_tot_pp = area_mask*(10**5.0)*o3_bau_mort_tot[0,:,:]/pop
o3_cov_mort_tot_pp = area_mask*(10**5.0)*o3_cov_mort_tot[0,:,:]/pop

# make places where pop is 0 also 0 for the totals, rather than nans for plotting
pm_bau_mort_tot_pp_4plot = np.where(pop==0,0,pm_bau_mort_tot_pp)
pm_cov_mort_tot_pp_4plot = np.where(pop==0,0,pm_cov_mort_tot_pp)
o3_bau_mort_tot_pp_4plot = np.where(pop==0,0,o3_bau_mort_tot_pp)
o3_cov_mort_tot_pp_4plot = np.where(pop==0,0,o3_cov_mort_tot_pp)

# calculate totals to add to plot titles
pm_sum_bau_mort = np.nansum(area_mask*pm_bau_mort_tot[0,:,:])
pm_sum_cov_mort = np.nansum(area_mask*pm_cov_mort_tot[0,:,:])
o3_sum_bau_mort = np.nansum(area_mask*o3_bau_mort_tot[0,:,:])
o3_sum_cov_mort = np.nansum(area_mask*o3_cov_mort_tot[0,:,:])

plt_map(lon,lat,
        [o3_bau_mort_tot_pp_4plot,o3_cov_mort_tot_pp_4plot,
         o3_cov_mort_tot_pp_4plot-o3_bau_mort_tot_pp_4plot,
         pm_bau_mort_tot_pp_4plot,pm_cov_mort_tot_pp_4plot,
         pm_cov_mort_tot_pp_4plot-pm_bau_mort_tot_pp_4plot],
        1,
        ['magma','magma','seismic','magma','magma','seismic'],
        ['Deaths y$^{-1}$ per 10$^5$ People']*6,
        ['BAU\nO$_3$ '+str(int(round(o3_sum_bau_mort,-2))),
         'COV\nO$_3$ '+str(int(round(o3_sum_cov_mort,-2))),
         'Difference\nO$_3$ '+str(int(round(o3_sum_bau_mort-o3_sum_cov_mort,0))),
         'PM$_{2.5}$ '+str(int(round(pm_sum_bau_mort,-2))),
          'PM$_{2.5}$ '+str(int(round(pm_sum_cov_mort,-2))),
          'PM$_{2.5}$ '+str(int(round(pm_sum_bau_mort-pm_sum_cov_mort,0)))],
        clim = [[0,50],[0,50],[-5,5],
                [0,50],[0,50],[-5,5]],
        multi=[2,3],
        outname = figure_folder + 'total_mortality_mapsv2_final.png')
