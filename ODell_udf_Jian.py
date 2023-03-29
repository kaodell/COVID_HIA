#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ODell_udf_Jian.py
    python script of functions written by me or by others passed on to me
Created on Wed Sep  8 09:09:22 2021
@author: kodell
"""
#%% packages needed
import numpy as np
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib as mplt
from matplotlib import colors
mplt.rcParams['font.size'] = '14'
mplt.rcParams['font.family'] = 'sans-serif'
#mplt.rcParams['font.sans-serif'] = 'Veranda'

#%% make a basic map of the US using cartopy
# NOTE this assumes a PlateCaree projection
def plt_map(dlon,dlat,data,size,cmap,clabel,title,**kwargs):
    """Function creates a pcolormesh figure using cartopy mapping modules. 
    Can create a multi-panel figure and write figure to a file given user keyword args.
    
    Parameters
    ----------
    dlon : numpy.ndarray 
        2D array of longitude in degrees
    dlat : numpy.ndarray 
        2D array of latitude in degrees
    data : numpy.ndarry or tuple of numpy.ndarrys if multi is indicated
        2D array of the data to be plotted. must be on the same grid as other
        data if indicating multiple data grids to plot
    size: int
        size for scatter points (no longer needed with switch to pcolormesh, 
                                 but kept for ease of switch back to scatter)
    cmap: tuple
        string indicating cmap to use
    clabel: tuple
        string for colorbar label
    title: tuple
        string for figure title
    
    Keyword Arguments
    -------
    clim: tuple
        interger or float limits for the colorbar colors
        otherwise limits are chosen by matplotlib defaults
    outname: string
        filepath and name for file to save figure to
        otherwise figure is not written to a file
    cpts: tuple
        interter or float of min, middle, and max of colorbar
        otherwise limits are chosen by matplotlib defaults
    multi: tuple
        [n rows, n cols] for a multi panel figure. must match shape of data,
        cmap, clabel, and title. default is [1,1].
    bkcolor: string 
        color for plot background. default is white.
    norm: string
        normalization for colorbar if desired. cannot be used with cpts.
    
    Returns
    -------
    none. figure is shown and saved if outname keyword is used.
    ------
    written by Katelyn O'Dell
    """
    vlim = kwargs.get('clim', None)
    outpath = kwargs.get('outname',None)
    vpts = kwargs.get('cpts',None)
    multi = kwargs.get('multi',None)
    bkcolor = kwargs.get('bkcolor',None)
    norm = kwargs.get('norm',None)
    if multi:
        nd = len(data)
        if bkcolor:
            fig, axarr = plt.subplots(nrows=multi[0],ncols=multi[1],subplot_kw={'projection': ccrs.PlateCarree()},
                                      figsize=(11,8.5),facecolor=bkcolor)
        else:
            fig, axarr = plt.subplots(nrows=multi[0],ncols=multi[1],subplot_kw={'projection': ccrs.PlateCarree()},
                                      figsize=(11,8.5))
        axarr = axarr.flatten()
        
        for di in range(nd):
            ax = axarr[di]
            ax.patch.set_visible(False)
            ax.add_feature(cfeature.STATES.with_scale('50m'),edgecolor='gray')
            ax.set_extent([-125, -66, 25, 50], crs=ccrs.PlateCarree())

            if bkcolor:
                ax.outline_patch.set_edgecolor(bkcolor)
            else:
                ax.outline_patch.set_edgecolor('white')
            if vlim:
                cs = ax.pcolormesh(dlon,dlat,data[di],shading='nearest',
                            transform=ccrs.PlateCarree(),cmap=cmap[di],vmin=vlim[di][0],vmax=vlim[di][1])
            elif vpts:
                divnorm=colors.TwoSlopeNorm(vmin=vpts[di][0], vcenter=vpts[di][1], vmax=vpts[di][2])
                cs = ax.pcolormesh(dlon,dlat,data[di],shading='nearest',
                            transform=ccrs.PlateCarree(),cmap=cmap[di],norm=divnorm)
            elif norm:
                cs = ax.pcolormesh(dlon,dlat,data[di],shading='nearest',
                            transform=ccrs.PlateCarree(),cmap=cmap[di],norm=norm)
            else:
                cs = ax.pcolormesh(dlon,dlat,data[di],shading='nearest',
                            transform=ccrs.PlateCarree(),cmap=cmap[di])
            cbar = fig.colorbar(cs,ax=ax,orientation='horizontal',pad=0,shrink=0.6)
            
            if bkcolor:
                cbar.set_label(label=clabel[di],size=16,color='white')
                cbar.ax.xaxis.set_tick_params(color='white', labelcolor='white')
                ax.set_title(title[di],fontsize=18,color='white')
            else:
                cbar.set_label(label=clabel[di],size=16)
                ax.set_title(title[di],fontsize=18)
                
            plt.tight_layout()
    # now for case without a multi-panel figure (not used for the version for Jian's paper)
    else:
        if bkcolor:
            fig, ax = plt.subplots(nrows=1,ncols=1,
                                      subplot_kw={'projection': ccrs.PlateCarree()},
                                      figsize=(11,8.5),facecolor=bkcolor)
            ax.outline_patch.set_edgecolor(bkcolor)

        else:
            fig, ax = plt.subplots(nrows=1,ncols=1,
                                      subplot_kw={'projection': ccrs.PlateCarree()},
                                      figsize=(11,8.5))
            ax.outline_patch.set_edgecolor('white')

        ax.patch.set_visible(False)
        ax.add_feature(cfeature.STATES.with_scale('50m'),edgecolor='gray')
        ax.outline_patch.set_edgecolor('white')
        ax.set_extent([-125, -66, 24, 50], crs=ccrs.PlateCarree())

        if vlim:
            cs = ax.pcolormesh(dlon,dlat,data,shading='nearest',
                        transform=ccrs.PlateCarree(),cmap=cmap,vmin=vlim[0],vmax=vlim[1])
        elif vpts:
            divnorm=colors.TwoSlopeNorm(vmin=vpts[0], vcenter=vpts[1], vmax=vpts[2])
            cs = ax.pcolormesh(dlon,dlat,data,shading='nearest',
                        transform=ccrs.PlateCarree(),cmap=cmap,norm=divnorm)
        elif norm:
            cs = ax.pcolormesh(dlon,dlat,data,shading='nearest',
                        transform=ccrs.PlateCarree(),cmap=cmap,norm=norm)
        else:
            cs = ax.pcolormesh(dlon,dlat,data,shading='nearest',
                        transform=ccrs.PlateCarree(),cmap=cmap)
        cbar = fig.colorbar(cs,ax=ax,orientation='horizontal',pad=0,shrink=0.6)
        if bkcolor:
            cbar.set_label(label=clabel,size=16,color='white')
            ax.set_title(title,fontsize=18,color='white')
            cbar.ax.xaxis.set_tick_params(color='white', labelcolor='white')
        else:
            cbar.set_label(label=clabel,size=16)
            ax.set_title(title,fontsize=18)
        plt.tight_layout()

    if outpath:
        plt.savefig(outpath,dpi=400)
    plt.show()

#%% make map on a created axis
def mk_map(ax):
    """Function crates a map on cartopy axis.
    
    Parameters
    ----------
    ax: matplotlib figure axis
    
    Returns
    -------
    none.
    ------
    written by Katelyn O'Dell
    """

    ax.patch.set_visible(False)
    # plot shapfile with colors
    ax.add_feature(cfeature.LAND.with_scale('50m'),facecolor='gray',alpha=0.5)
    ax.add_feature(cfeature.OCEAN.with_scale('50m'))
    ax.add_feature(cfeature.STATES.with_scale('50m'),edgecolor='lightgray')
    ax.outline_patch.set_edgecolor('white')

#%% calculate wet pm2.5 at 35% RH and STP from wrf chem 
def calc_pmw(nc_fid,p25_str):
    """Function calculate wet pm2.5 at 35% RH and STP from wrf chem files from 
    Jian for GEO-XO project. Calculation based on email discussion between Jian He, 
    Brian McDonald,and Daven Henze and summarized in readme files 
    with the data from Jian.
    
    Parameters
    ----------
    nc_fid: loaded netCDF file ID
    
    Returns
    -------
    wet PM2.5 mass in ug/m3
    
    ------
    written by Katelyn O'Dell
    """
    so4a = nc_fid['so4a'][:]
    nh4a = nc_fid['nh4a'][:]
    no3a = nc_fid['no3a'][:]
    ec = nc_fid['ec'][:]
    orgpa = nc_fid['orgpa'][:]
    soa = nc_fid['soa'][:]
    p25 = nc_fid[p25_str][:]
    naa = nc_fid['naa'][:]
    cla = nc_fid['cla'][:]
    p = nc_fid['pres'][:]
    t = nc_fid['temp'][:]
    pm25w = (1.1*(so4a + nh4a + no3a) + ec + orgpa + 1.05*soa + p25 +1.86*(naa + cla))*(101325/p)*(t/298.0)
    return pm25w

#%% haversine function from Dr. Will Lassman
def haversine(lon0,lon1,lat0,lat1):
    """Function calculates distance between two lat/lon points assuming
    earth radius of 6,371 km
    
    Parameters
    ----------
    lon0: longitude of point 0 in degrees
    lon1: longitude of point 1 in degrees
    lat0: latitude of point 0 in degrees
    lat1: latitude of point 1 in degrees
    
    Returns
    -------
    Distance between point 0 and point 1 in meters.
    
    ------
    written by Will Lassman
    """

    r = 6371000. #m # mean earth radius                                                                                                                                                                                                                                              
    lon0 = lon0*np.pi/180

    lon1 = lon1*np.pi/180

    lat0 = lat0*np.pi/180

    lat1 = lat1*np.pi/180

    return 2*r*np.arcsin(np.sqrt(np.sin((lat1 - lat0)/2.)**2 +\
		 np.cos(lat0)*np.cos(lat1)*np.sin((lon1 - lon0)/2.)**2))
   
#%% chronic HIA function - generic
def chronic_HIA(conc, cf, pop, base_rate, betas, grid_area): 
    # beta is array of beta calculated from [rr,rr_lci,rr_uci]
    """Function performs a health impact assessment to calculates annual events 
    due to to long term exposure to a pollutant using the health impact function 
    from Anenberg et al., 2010 https://doi.org/10.1289/ehp.0901220
    
    Parameters
    ----------
    conc: numpy.ndarray 
        gridded annual mean pollutant concentration in the same units used for beta definition
    cf: float
        threshold concentration of minimal risk (conc level at which no excess risk is assumed)
    pop: numpy.ndarray
        gridded population counts on the same grid as conc
    base_rate: numpy.ndarray
        gridded baseline mortality rates per person on the same grid as conc
    betas: list
        list of betas for the health impact function
    grid_area: numpy.ndarray
        area of each grid cell in km for the conc grid
    
    Returns
    -------
    paf_avg_out2: numpy.ndarray
        population attributable fraction array with shape: nBetas, gridx, gridy
    events_tot_out2: numpy.ndarray
        total number of attributable events per year with shape: nBetas, gridx, gridy
    events_tot_pp_out2: numpy.ndarray
        total number of attributable events per year per 10^5 people with shape: nBetas, gridx, gridy
    events_tot_pk_out2: numpy.ndarray
        total number of attributable events per year per km with shape: nBetas, gridx, gridy
    
    ------
    written by Katelyn O'Dell
    """
    paf_avg_out = []
    events_tot_out = []
    events_tot_pp_out = []
    events_tot_pk_out = []
    z = np.where(conc<cf,0,conc-cf)
    for beta in betas:      
        paf_avg = 1.0 - np.exp(-beta*z)
        events = (paf_avg)*pop*(base_rate)
        events_tot = events
        events_tot_pk = events_tot/grid_area
        events_tot_pp = events_tot/(pop/100000)
        
        paf_avg_out.append(paf_avg)
        events_tot_out.append(events_tot)
        events_tot_pp_out.append(events_tot_pp)
        events_tot_pk_out.append(events_tot_pk)
    
    paf_avg_out2 = np.reshape(paf_avg_out,[len(betas),conc.shape[0],conc.shape[1]])
    events_tot_out2 = np.reshape(events_tot_out,[len(betas),conc.shape[0],conc.shape[1]])
    events_tot_pp_out2 = np.reshape(events_tot_pp_out,[len(betas),conc.shape[0],conc.shape[1]])
    events_tot_pk_out2 = np.reshape(events_tot_pk_out,[len(betas),conc.shape[0],conc.shape[1]])

            
    return paf_avg_out2, events_tot_out2, events_tot_pp_out2, events_tot_pk_out2


