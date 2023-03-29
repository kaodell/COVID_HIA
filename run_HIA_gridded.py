#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_HIA_gridded.py
    python script to run the HIA at the gridcell level
    09.22.22 - added annual average for mortality HIA to code
written by: Katelyn O'Dell
Created on Mon Nov  1 15:16:11 2021
"""
#%% user inputs
# local path to project folder
prj_folder = '/Users/kodell/Library/CloudStorage/GoogleDrive-kodell@email.gwu.edu/My Drive/Ongoing Projects/GeoXO/'

# WRF simulation to use
sim_base = 'BAU'
sim_cov = 'COV'
base_name = 'Baseline'
cov_name = 'COVID'

# pollutant
poll = 'o3'# o3, or pm25w
# file description from create_loc_daily_avg.py output files
desc = 'full_annual_ts'

# gridded population
pop_fn = 'population_data/NASA_SEDAC_pop/regridded/regrided_2020pop_12km_final.nc'

# gridded mortality baseline rates 
br_mort_fn = prj_folder + 'health_data/gridded_cnty_mort_2_final.csv'

# designate paths to save outputs
out_path_mort = prj_folder + 'HIA_results/Jian/'+poll + '_' + 'mort_hia_full_cntybr_final'

#%% import modules
import pandas as pd
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from ODell_udf_Jian import plt_map, chronic_HIA
import sys

#%% set additional parameters based on pollutant choice
if poll == 'pm25w':
    tmrl =  2.8 # theretical minimum risk level for mortality, minimum conc estiamted in Turner at al., 2016
    mort_RRs = [1.06,1.04,1.08] # relative risks [central, l95%ci, u95%ci] from Turner et al., 2016; 
                                # Table E10; single pollutant model, same as Maria paper
    avg_time = '24hravg' # use mda8 for ozone and 24hravg for PM
elif poll =='o3':
    tmrl = 26.7 # same as above but for ozone
    mort_RRs = [1.02,1.01,1.04] # same as above but for ozone (per 10 ppb, annual mean of mda8)
    avg_time = 'mda8'
else:
    sys.exit('pollutant name not recognized')
    
#%% load and prep data
# gridded pollutant: base case
nc_fid = Dataset(prj_folder+'concentration_data/Jian_WRF_processed_lda/WRFChem_'+sim_base+'_'+poll+'_'+avg_time+'_'+desc+'.nc')
lat = nc_fid['lat'][:]
lon = nc_fid['lon'][:]
dates = nc_fid['date'][:]
poll_bau = nc_fid[poll][:]   
nc_fid.close()

# baseline mortality rate
br_mort_df = pd.read_csv(br_mort_fn)
br1 = br_mort_df['cnty_mort_rate'].values/100000.0
br2 = br1.reshape([lat.shape[0],lat.shape[1]])
br = np.where(br2==-7777.7,np.nan,br2)

# gridded pollutant: covid case
nc_fid = Dataset(prj_folder+'concentration_data/Jian_WRF_processed_lda/WRFChem_'+sim_cov+'_'+poll+'_'+avg_time+'_'+desc+'.nc')
poll_cov = nc_fid[poll][:]  
 
nc_fid.close()

# gridded population counts
nc_fid = Dataset(prj_folder+pop_fn)
pop_all = nc_fid['population'][:].data
grid_area = nc_fid['grid_area'][:].data
us_mask = nc_fid['us_mask'][:].data
nc_fid.close()
# cut population to the US
pop = pop_all*us_mask

#%% do HIA for long-term exposure - pm and ozone mortality
# calculate annual average
poll_bau_avg = np.nanmean(poll_bau,axis=0)
poll_cov_avg = np.nanmean(poll_cov,axis=0)

# calculate betas from RRs
mort_betas = np.log(np.array(mort_RRs))/10

# calculate attributable events for the bau case
[paf_avg_bau, bau_mort_tot_out, events_tot_pp_out, events_tot_pk_out] = chronic_HIA(poll_bau_avg, tmrl, pop, 
                                                                                 br,mort_betas, grid_area)
# and for the cov case
[cov_avg_bau, cov_mort_tot_out, events_tot_pp_out, events_tot_pk_out] = chronic_HIA(poll_cov_avg, tmrl, pop, 
                                                                                 br,mort_betas, grid_area)

# print output
bau_mort_sum = np.nansum(bau_mort_tot_out[0,:,:]) # 0 = central estimate 
cov_mort_sum = np.nansum(cov_mort_tot_out[0,:,:])
print('base - adj deaths:',cov_mort_sum-bau_mort_sum,'(',
      100.0*(cov_mort_sum-bau_mort_sum)/bau_mort_sum,'% )')

#%% save output to npz file
np.savez(out_path_mort,
         bau_mort_tot = bau_mort_tot_out,
         cov_mort_tot = cov_mort_tot_out,
         poll_bau_avg = poll_bau_avg,
         poll_cov_avg = poll_cov_avg)

#%% plot inputs, HIA results, and differences
# timeseries of pop-weighted pollutant differences
popw_bau = np.nansum(poll_bau*pop,axis=(1,2))/np.nansum(pop)
popw_cov = np.nansum(poll_cov*pop,axis=(1,2))/np.nansum(pop)
fig, ax = plt.subplots()
ax.plot(dates,popw_bau-popw_cov)
ax.set_xticks(['2020_04_01','2020_08_01','2020_12_01','2021_03_31'],
              ['1-April-20','1-Sept-20', '1-Dec-20','1-April-21'])
ax.set_ylabel('Difference, US pop-weighted\n '+poll+' [$\mu$g m$^{-3}$]')
ax.set_xlabel('Date')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.savefig(prj_folder + 'figures/preliminary/concentrations/pop_w_timerseries_'+poll+'.png',dpi=300)
fig.show()

# map of pollutant differences
diff = poll_bau_avg-poll_cov_avg
plt_map(lon,lat,
        [poll_bau_avg,poll_cov_avg,diff,100.0*diff/poll_bau_avg],
        1,
        ['jet','jet','seismic','seismic'],
        ['Surface '+poll,'Surface '+poll,poll+' difference',poll+' % difference [%]'],
        [base_name+' case',cov_name+' case',base_name+' case - '+cov_name+' case',
         base_name+' case - '+cov_name+' case'],
        clim = [[0,60],[0,60],[-5,5],[-50,50]],
        multi=[2,2],
        outname = prj_folder + 'figures/preliminary/concentrations/no2_concentrations/compare_'+poll+'_conc_bc.png')
 