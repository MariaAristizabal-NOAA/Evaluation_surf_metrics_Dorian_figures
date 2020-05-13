#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  4 14:59:17 2020

@author: aristizabal
"""

#%% User input

lon_lim = [-98.5,-60.0]
lat_lim = [10.0,45.0]

# Server erddap url IOOS glider dap
server = 'https://data.ioos.us/gliders/erddap'

#gliders sg666, sg665, sg668, silbo
gdata_sg665 = 'http://gliders.ioos.us/thredds/dodsC/deployments/aoml/SG665-20190718T1155/SG665-20190718T1155.nc3.nc'
gdata_sg666 = 'http://gliders.ioos.us/thredds/dodsC/deployments/aoml/SG666-20190718T1206/SG666-20190718T1206.nc3.nc'
gdata_sg668 = 'http://gliders.ioos.us/thredds/dodsC/deployments/aoml/SG668-20190819T1217/SG668-20190819T1217.nc3.nc'

# Time bounds
min_time = '2019/08/20/00'
max_time = '2019/09/07/00'

# Bathymetry file
bath_file = '/Users/aristizabal/Desktop/MARACOOS_project/Maria_scripts/nc_files/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'

# KMZ file
kmz_file = '/Users/aristizabal/Desktop/MARACOOS_project/Maria_scripts/KMZ_files/al052019_best_track-5.kmz'

cycle = '2019082800'

delta_lon = 0 # delta longitude around hurricane track to calculate
               # statistics
Nini = 0 # 0 is the start of forecating cycle (2019082800)
      # 1 is 6 hours of forecasting cycle   (2019082806)
      # 2 is 12 hours ...... 20 is 120 hours 

Nend = 22 # indicates how far in the hurricabe track you want
          # include in the analysis. This is helpful if for ex:
          # you onl want to analyse the portion of the track
          # where the storm intensifies
          # 22 corresponds to all the hurricane track forecasted in a cycle
#Nend = 13

# Bathymetry file
bath_file = '/home/aristizabal/bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'

# KMZ file best track Dorian
kmz_file_Dorian = '/home/aristizabal/KMZ_files/al052019_best_track-5.kmz'

# url for GOFS 3.1
url_GOFS = 'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_93.0/ts3z'

# figures
folder_fig = '/home/aristizabal/Figures/'

# folder nc files POM
folder_pom19 =  '/home/aristizabal/HWRF2019_POM_Dorian/'
folder_pom20 =  '/home/aristizabal/HWRF2020_POM_Dorian/'

# folde HWRF2020_HYCOM
folder_hycom20 = '/home/aristizabal/HWRF2020_HYCOM_Dorian/'

###################

# folder nc files POM
folder_pom_oper = folder_pom19 + 'HWRF2019_POM_dorian05l.' + cycle + '_pom_files_oper/'
folder_pom_exp = folder_pom20 + 'HWRF2020_POM_dorian05l.'  + cycle + '_pom_files_exp/'
prefix_pom = 'dorian05l.' + cycle + '.pom.00'

pom_grid_oper = folder_pom_oper + 'dorian05l.' + cycle + '.pom.grid.nc'
pom_grid_exp = folder_pom_exp + 'dorian05l.' + cycle + '.pom.grid.nc'

# Dorian track files
hwrf_pom_track_oper = folder_pom_oper + 'dorian05l.' + cycle + '.trak.hwrf.atcfunix'
hwrf_pom_track_exp = folder_pom_exp + 'dorian05l.' + cycle + '.trak.hwrf.atcfunix'

# folder nc files hwrf
folder_hwrf_pom19_oper = folder_pom19 + 'HWRF2019_POM_dorian05l.' + cycle + '_grb2_to_nc_oper/'
folder_hwrf_pom20_exp = folder_pom20 + 'HWRF2020_POM_dorian05l.' + cycle + '_grb2_to_nc_exp/'

##################
# folder ab files HYCOM
folder_hycom_exp = folder_hycom20 + 'HWRF2020_HYCOM_dorian05l.' + cycle + '_hycom_files_exp/'
prefix_hycom = 'dorian05l.' + cycle + '.hwrf_rtofs_hat10_3z'

#Dir_HMON_HYCOM = '/Volumes/aristizabal/ncep_model/HMON-HYCOM_Michael/'
Dir_HMON_HYCOM = '/home/aristizabal/ncep_model/HWRF-Hycom-WW3_exp_Michael/'
# RTOFS grid file name
hycom_grid_exp = Dir_HMON_HYCOM + 'hwrf_rtofs_hat10.basin.regional.grid'

# Dorian track files
hwrf_hycom_track_exp = folder_hycom_exp + 'dorian05l.' + cycle + '.trak.hwrf.atcfunix'

# folder nc files hwrf
folder_hwrf_hycom20_exp = folder_hycom20 + 'HWRF2020_HYCOM_dorian05l.' + cycle + '_grb2_to_nc_exp/'

#%%
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import os
import os.path
import glob
#import netCDF4
from bs4 import BeautifulSoup
from zipfile import ZipFile
#import matplotlib.patches as patches
import cmocean
from erddapy import ERDDAP
import pandas as pd
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)

import sys
sys.path.append('/home/aristizabal/glider_model_comparisons_Python')
from read_glider_data import read_glider_data_thredds_server

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

#%% Function Grid glider variables according to depth

def varsg_gridded(depth,time,temp,salt,dens,delta_z):
             
    depthg_gridded = np.arange(0,np.nanmax(depth),delta_z)
    tempg_gridded = np.empty((len(depthg_gridded),len(time)))
    tempg_gridded[:] = np.nan
    saltg_gridded = np.empty((len(depthg_gridded),len(time)))
    saltg_gridded[:] = np.nan
    densg_gridded = np.empty((len(depthg_gridded),len(time)))
    densg_gridded[:] = np.nan

    for t,tt in enumerate(time):
        depthu,oku = np.unique(depth[:,t],return_index=True)
        tempu = temp[oku,t]
        saltu = salt[oku,t]
        densu = dens[oku,t]
        okdd = np.isfinite(depthu)
        depthf = depthu[okdd]
        tempf = tempu[okdd]
        saltf = saltu[okdd]
        densf = densu[okdd]
 
        okt = np.isfinite(tempf)
        if np.sum(okt) < 3:
            temp[:,t] = np.nan
        else:
            okd = np.logical_and(depthg_gridded >= np.min(depthf[okt]),\
                                 depthg_gridded < np.max(depthf[okt]))
            tempg_gridded[okd,t] = np.interp(depthg_gridded[okd],depthf[okt],tempf[okt])
            
        oks = np.isfinite(saltf)
        if np.sum(oks) < 3:
            saltg_gridded[:,t] = np.nan
        else:
            okd = np.logical_and(depthg_gridded >= np.min(depthf[okt]),\
                        depthg_gridded < np.max(depthf[okt]))
            saltg_gridded[okd,t] = np.interp(depthg_gridded[okd],depthf[oks],saltf[oks])
    
        okdd = np.isfinite(densf)
        if np.sum(okdd) < 3:
            densg_gridded[:,t] = np.nan
        else:
            okd = np.logical_and(depthg_gridded >= np.min(depthf[okdd]),\
                        depthg_gridded < np.max(depthf[okdd]))
            densg_gridded[okd,t] = np.interp(depthg_gridded[okd],depthf[okdd],densf[okdd])
        
    return depthg_gridded, tempg_gridded, saltg_gridded, densg_gridded
           
#%% Get latent and sensible heat fluxes at glider position and time 
def get_sensible_latent_heat_fluxes_HWRF(ncfile_list,latg,long,timeg):

    shtfl_hwrf = []
    lhtfl_hwrf = []
    time_hwrf = []
    
    for N,file in enumerate(ncfile_list):
        print(N)
        HWRF = xr.open_dataset(file)
        lat_hwrf = np.asarray(HWRF.variables['latitude'][:])
        lon_hwrf = np.asarray(HWRF.variables['longitude'][:])
        t_hwrf = np.asarray(HWRF.variables['time'][:])
        SHTFL_hwrf = np.asarray(HWRF.variables['SHTFL_surface'][0,:,:])
        LHTFL_hwrf = np.asarray(HWRF.variables['LHTFL_surface'][0,:,:])
        
        # Changing times to timestamp
        tstamp_glider = [mdates.date2num(timeg[i]) for i in np.arange(len(timeg))]
        tstamp_hwrf = [mdates.date2num(t_hwrf[i]) for i in np.arange(len(t_hwrf))]
        
        # interpolating glider lon and lat to lat and lon on model time
        sublon_hwrf = np.interp(tstamp_hwrf,tstamp_glider,long)
        sublat_hwrf = np.interp(tstamp_hwrf,tstamp_glider,latg)
        
        # getting the model grid positions for sublonm and sublatm
        oklon_hwrf = np.round(np.interp(sublon_hwrf,lon_hwrf,np.arange(len(lon_hwrf)))).astype(int)
        oklat_hwrf = np.round(np.interp(sublat_hwrf,lat_hwrf,np.arange(len(lat_hwrf)))).astype(int)
            
        shtfl_hwrf.append(SHTFL_hwrf[oklat_hwrf,oklon_hwrf][0])
        lhtfl_hwrf.append(LHTFL_hwrf[oklat_hwrf,oklon_hwrf][0])
        time_hwrf.append(t_hwrf[0])
        
    shtfl_hwrf = np.asarray(shtfl_hwrf)
    lhtfl_hwrf = np.asarray(lhtfl_hwrf)
    time_hwrf = np.asarray(time_hwrf)
    
    return shtfl_hwrf, lhtfl_hwrf, time_hwrf

#%%
def MLD_temp_and_dens_criteria(dt,drho,time,depth,temp,salt,dens):

    MLD_temp_crit = np.empty(temp.shape[1]) 
    MLD_temp_crit[:] = np.nan
    Tmean_temp_crit = np.empty(temp.shape[1]) 
    Tmean_temp_crit[:] = np.nan
    Smean_temp_crit = np.empty(temp.shape[1]) 
    Smean_temp_crit[:] = np.nan
    Td_temp_crit = np.empty(temp.shape[1]) 
    Td_temp_crit[:] = np.nan
    MLD_dens_crit = np.empty(temp.shape[1])
    MLD_dens_crit[:] = np.nan
    Tmean_dens_crit = np.empty(temp.shape[1])
    Tmean_dens_crit[:] = np.nan
    Smean_dens_crit = np.empty(temp.shape[1]) 
    Smean_dens_crit[:] = np.nan
    Td_dens_crit = np.empty(temp.shape[1]) 
    Td_dens_crit[:] = np.nan
    for t in np.arange(temp.shape[1]):
        if depth.ndim == 1:
            d10 = np.where(depth >= 10)[0][0]
        if depth.ndim == 2:
            d10 = np.where(depth[:,t] >= -10)[0][-1]
        T10 = temp[d10,t]
        delta_T = T10 - temp[:,t] 
        ok_mld_temp = np.where(delta_T <= dt)[0]
        rho10 = dens[d10,t]
        delta_rho = -(rho10 - dens[:,t])
        ok_mld_rho = np.where(delta_rho <= drho)[0]
        
        if ok_mld_temp.size == 0:
            MLD_temp_crit[t] = np.nan
            Td_temp_crit[t] = np.nan
            Tmean_temp_crit[t] = np.nan
            Smean_temp_crit[t] = np.nan            
        else:                             
            if depth.ndim == 1:
                MLD_temp_crit[t] = depth[ok_mld_temp[-1]]
                ok_mld_plus1m = np.where(depth >= depth[ok_mld_temp[-1]] + 1)[0][0]                 
            if depth.ndim == 2:
                MLD_temp_crit[t] = depth[ok_mld_temp[-1],t]
                ok_mld_plus1m = np.where(depth >= depth[ok_mld_temp[-1],t] + 1)[0][0]
            Td_temp_crit[t] = temp[ok_mld_plus1m,t]
            Tmean_temp_crit[t] = np.nanmean(temp[ok_mld_temp,t])
            Smean_temp_crit[t] = np.nanmean(salt[ok_mld_temp,t])
                
        if ok_mld_rho.size == 0:
            MLD_dens_crit[t] = np.nan
            Td_dens_crit[t] = np.nan
            Tmean_dens_crit[t] = np.nan
            Smean_dens_crit[t] = np.nan           
        else:
            if depth.ndim == 1:
                MLD_dens_crit[t] = depth[ok_mld_rho[-1]]
                ok_mld_plus1m = np.where(depth >= depth[ok_mld_rho[-1]] + 1)[0][0] 
            if depth.ndim == 2:
                MLD_dens_crit[t] = depth[ok_mld_rho[-1],t]
                ok_mld_plus1m = np.where(depth >= depth[ok_mld_rho[-1],t] + 1)[0][0] 
            Td_dens_crit[t] = temp[ok_mld_plus1m,t]        
            Tmean_dens_crit[t] = np.nanmean(temp[ok_mld_rho,t])
            Smean_dens_crit[t] = np.nanmean(salt[ok_mld_rho,t]) 

    return MLD_temp_crit,Tmean_temp_crit,Smean_temp_crit,Td_temp_crit,\
           MLD_dens_crit,Tmean_dens_crit,Smean_dens_crit,Td_dens_crit

#%% Read best storm track from kmz file
    
def read_kmz_file_storm_best_track(kmz_file):
    
    os.system('cp ' + kmz_file + ' ' + kmz_file[:-3] + 'zip')
    os.system('unzip -o ' + kmz_file + ' -d ' + kmz_file[:-4])
    kmz = ZipFile(kmz_file[:-3]+'zip', 'r')
    kml_file = kmz_file.split('/')[-1].split('_')[0] + '.kml'
    kml_best_track = kmz.open(kml_file, 'r').read()
    
    # best track coordinates
    soup = BeautifulSoup(kml_best_track,'html.parser')
    
    lon_best_track = np.empty(len(soup.find_all("point")))
    lon_best_track[:] = np.nan
    lat_best_track = np.empty(len(soup.find_all("point")))
    lat_best_track[:] = np.nan
    for i,s in enumerate(soup.find_all("point")):
        lon_best_track[i] = float(s.get_text("coordinates").split('coordinates')[1].split(',')[0])
        lat_best_track[i] = float(s.get_text("coordinates").split('coordinates')[1].split(',')[1])
             
    #  get time stamp
    time_best_track = []
    for i,s in enumerate(soup.find_all("atcfdtg")):
        tt = datetime.strptime(s.get_text(' '),'%Y%m%d%H')
        time_best_track.append(tt)
    time_best_track = np.asarray(time_best_track)    
    
    # get type 
    wind_int_mph = []
    for i,s in enumerate(soup.find_all("intensitymph")):
        wind_int_mph.append(s.get_text(' ')) 
    wind_int_mph = np.asarray(wind_int_mph)
    wind_int_mph = wind_int_mph.astype(float)  
    
    wind_int_kt = []
    for i,s in enumerate(soup.find_all("intensity")):
        wind_int_kt.append(s.get_text(' ')) 
    wind_int_kt = np.asarray(wind_int_kt)
    wind_int_kt = wind_int_kt.astype(float)
      
    cat = []
    for i,s in enumerate(soup.find_all("styleurl")):
        cat.append(s.get_text('#').split('#')[-1]) 
    cat = np.asarray(cat)
    
    mslp = []
    for i,s in enumerate(soup.find_all("minsealevelpres")):
        mslp.append(s.get_text('#').split('#')[-1]) 
    mslp = np.asarray(mslp)
    
    return lon_best_track, lat_best_track, time_best_track, wind_int_mph, wind_int_kt, cat, mslp 

#%%
def get_max_winds_10m_HWRF(folder_nc_files): 
    
    files = sorted(glob.glob(os.path.join(folder_nc_files,'*.nc')))
    max_wind_10m_hwrf = []
    time_hwrf = []
    for i,fl in enumerate(files):
        print(i)
        HWRF = xr.open_dataset(fl)
        t_hwrf = np.asarray(HWRF.variables['time'][:])
        UGRD_hwrf = np.asarray(HWRF.variables['UGRD_10maboveground'][0,:,:])
        VGRD_hwrf = np.asarray(HWRF.variables['VGRD_10maboveground'][0,:,:])
        max_wind_10m_hwrf.append(np.max(np.sqrt(UGRD_hwrf**2 + VGRD_hwrf**2)))
        time_hwrf.append(t_hwrf)
        
    return max_wind_10m_hwrf, time_hwrf

#%% Reading bathymetry data
ncbath = xr.open_dataset(bath_file)
bath_lat = ncbath.variables['lat'][:]
bath_lon = ncbath.variables['lon'][:]
bath_elev = ncbath.variables['elevation'][:]

oklatbath = np.logical_and(bath_lat >= lat_lim[0],bath_lat <= lat_lim[-1])
oklonbath = np.logical_and(bath_lon >= lon_lim[0],bath_lon <= lon_lim[-1])

bath_latsub = bath_lat[oklatbath]
bath_lonsub = bath_lon[oklonbath]
bath_elevs = bath_elev[oklatbath,:]
bath_elevsub = bath_elevs[:,oklonbath]

#%% Get Dorian best track and winds

lon_best_track, lat_best_track, time_best_track, _, wind_int_kt, cat, _= \
read_kmz_file_storm_best_track(kmz_file_Dorian)

#%% Get list HWRF files

HWRF_POM_oper = sorted(glob.glob(os.path.join(folder_hwrf_pom19_oper,'*.nc')))
HWRF_POM_exp = sorted(glob.glob(os.path.join(folder_hwrf_pom20_exp,'*.nc')))
HWRF_HYCOM_exp = sorted(glob.glob(os.path.join(folder_hwrf_hycom20_exp,'*.nc')))

#%% Get winds in knots at 10 meters height
    
folder_nc_files = folder_hwrf_pom19_oper
max_wind_10m_hwrf_pom19_oper, time_hwrf = get_max_winds_10m_HWRF(folder_nc_files)

folder_nc_files = folder_hwrf_pom20_exp
max_wind_10m_hwrf_pom20_exp, _ = get_max_winds_10m_HWRF(folder_nc_files)

folder_nc_files = folder_hwrf_hycom20_exp
max_wind_10m_hwrf_hycom20_exp, _ = get_max_winds_10m_HWRF(folder_nc_files)

#%% wind speed in knots

max_wind_10m_hwrf_pom19_oper = 1.94384 * np.asarray(max_wind_10m_hwrf_pom19_oper)
max_wind_10m_hwrf_pom20_exp = 1.94384 * np.asarray(max_wind_10m_hwrf_pom20_exp)
max_wind_10m_hwrf_hycom20_exp = 1.94384 * np.asarray(max_wind_10m_hwrf_hycom20_exp)

#%% Reading glider data
#Time window
date_ini = cycle[0:4]+'/'+cycle[4:6]+'/'+cycle[6:8]+'/'+cycle[8:]+'/00/00'
tini = datetime.strptime(date_ini,'%Y/%m/%d/%H/%M/%S')
tend = tini + timedelta(hours=126)
date_end = tend.strftime('%Y/%m/%d/%H/%M/%S')

#%% SG665  
  
url_glider = gdata_sg665

var_name = 'temperature'
scatter_plot = 'no'
kwargs = dict(date_ini=date_ini[0:-6],date_end=date_end[0:-6])
             
varg, timeg, latg, long, depthg, dataset_id = \
             read_glider_data_thredds_server(url_glider,var_name,scatter_plot,**kwargs)
             
tempg = varg  

var_name = 'salinity'  
varg, timeg, latg, long, depthg, dataset_id = \
             read_glider_data_thredds_server(url_glider,var_name,scatter_plot,**kwargs)
            
saltg = varg
 
var_name = 'density'  
varg, timeg, latg, long, depthg, dataset_id = \
             read_glider_data_thredds_server(url_glider,var_name,scatter_plot,**kwargs)
             
densg = varg
depthg = depthg              
             
# Grid glider variables according to depth

delta_z = 0.5
depthg_gridded,tempg_gridded,saltg_gridded,densg_gridded = \
varsg_gridded(depthg,timeg,tempg,saltg,densg,delta_z)  

# Calculation of mixed layer depth based on temperature and density critria
# Tmean: mean temp within the mixed layer and 
# td: temp at 1 meter below the mixed layer            
dt = 0.2
drho = 0.125

# for glider data
MLD_temp_crit_glid_sg665, _, _, _, MLD_dens_crit_glid, Tmean_dens_crit_glid, Smean_dens_crit_glid, _ = \
MLD_temp_and_dens_criteria(dt,drho,timeg,depthg_gridded,tempg_gridded,saltg_gridded,densg_gridded)             

# Get latent and sensible heat fluxes at glider position and time 
  
shtfl_hwrf19_pom, lhtfl_hwrf19_pom, time_hwrf19_pom = \
           get_sensible_latent_heat_fluxes_HWRF(HWRF_POM_oper,latg,long,timeg)
           
shtfl_hwrf20_pom, lhtfl_hwrf20_pom, time_hwrf20_pom = \
           get_sensible_latent_heat_fluxes_HWRF(HWRF_POM_exp,latg,long,timeg)
         
shtfl_hwrf20_hycom, lhtfl_hwrf20_hycom,time_hwrf20_hycom = \
          get_sensible_latent_heat_fluxes_HWRF(HWRF_HYCOM_exp,latg,long,timeg)
          
enth_hwrf19_pom_sg665 = shtfl_hwrf19_pom + lhtfl_hwrf19_pom
enth_hwrf20_pom_sg665 = shtfl_hwrf20_pom + lhtfl_hwrf20_pom
enth_hwrf20_hycom_sg665 = shtfl_hwrf20_hycom + lhtfl_hwrf20_hycom
dataset_id_sg665 = dataset_id
          
#%% SG666  
  
url_glider = gdata_sg666

var_name = 'temperature'
scatter_plot = 'no'
kwargs = dict(date_ini=date_ini[0:-6],date_end=date_end[0:-6])
             
varg, timeg, latg, long, depthg, dataset_id = \
             read_glider_data_thredds_server(url_glider,var_name,scatter_plot,**kwargs)
             
tempg = varg  

var_name = 'salinity'  
varg, timeg, latg, long, depthg, dataset_id = \
             read_glider_data_thredds_server(url_glider,var_name,scatter_plot,**kwargs)
            
saltg = varg
 
var_name = 'density'  
varg, timeg, latg, long, depthg, dataset_id = \
             read_glider_data_thredds_server(url_glider,var_name,scatter_plot,**kwargs)
             
densg = varg
depthg = depthg              
             
# Grid glider variables according to depth

delta_z = 0.5
depthg_gridded,tempg_gridded,saltg_gridded,densg_gridded = \
varsg_gridded(depthg,timeg,tempg,saltg,densg,delta_z)  

# Calculation of mixed layer depth based on temperature and density critria
# Tmean: mean temp within the mixed layer and 
# td: temp at 1 meter below the mixed layer            
dt = 0.2
drho = 0.125

# for glider data
MLD_temp_crit_glid_sg665, _, _, _, MLD_dens_crit_glid, Tmean_dens_crit_glid, Smean_dens_crit_glid, _ = \
MLD_temp_and_dens_criteria(dt,drho,timeg,depthg_gridded,tempg_gridded,saltg_gridded,densg_gridded)             

# Get latent and sensible heat fluxes at glider position and time 
  
shtfl_hwrf19_pom, lhtfl_hwrf19_pom, time_hwrf19_pom = \
           get_sensible_latent_heat_fluxes_HWRF(HWRF_POM_oper,latg,long,timeg)
           
shtfl_hwrf20_pom, lhtfl_hwrf20_pom, time_hwrf20_pom = \
           get_sensible_latent_heat_fluxes_HWRF(HWRF_POM_exp,latg,long,timeg)
         
shtfl_hwrf20_hycom, lhtfl_hwrf20_hycom,time_hwrf20_hycom = \
          get_sensible_latent_heat_fluxes_HWRF(HWRF_HYCOM_exp,latg,long,timeg)

enth_hwrf19_pom_sg666 = shtfl_hwrf19_pom + lhtfl_hwrf19_pom
enth_hwrf20_pom_sg666 = shtfl_hwrf20_pom + lhtfl_hwrf20_pom
enth_hwrf20_hycom_sg666 = shtfl_hwrf20_hycom + lhtfl_hwrf20_hycom
dataset_id_sg666 = dataset_id

#%% SG668 
  
url_glider = gdata_sg668

var_name = 'temperature'
scatter_plot = 'no'
kwargs = dict(date_ini=date_ini[0:-6],date_end=date_end[0:-6])
             
varg, timeg, latg, long, depthg, dataset_id = \
             read_glider_data_thredds_server(url_glider,var_name,scatter_plot,**kwargs)
             
tempg = varg  

var_name = 'salinity'  
varg, timeg, latg, long, depthg, dataset_id = \
             read_glider_data_thredds_server(url_glider,var_name,scatter_plot,**kwargs)
            
saltg = varg
 
var_name = 'density'  
varg, timeg, latg, long, depthg, dataset_id = \
             read_glider_data_thredds_server(url_glider,var_name,scatter_plot,**kwargs)
             
densg = varg
depthg = depthg              
             
# Grid glider variables according to depth

delta_z = 0.5
depthg_gridded,tempg_gridded,saltg_gridded,densg_gridded = \
varsg_gridded(depthg,timeg,tempg,saltg,densg,delta_z)  

# Calculation of mixed layer depth based on temperature and density critria
# Tmean: mean temp within the mixed layer and 
# td: temp at 1 meter below the mixed layer            
dt = 0.2
drho = 0.125

# for glider data
MLD_temp_crit_glid_sg665, _, _, _, MLD_dens_crit_glid, Tmean_dens_crit_glid, Smean_dens_crit_glid, _ = \
MLD_temp_and_dens_criteria(dt,drho,timeg,depthg_gridded,tempg_gridded,saltg_gridded,densg_gridded)             

# Get latent and sensible heat fluxes at glider position and time 
  
shtfl_hwrf19_pom, lhtfl_hwrf19_pom, time_hwrf19_pom = \
           get_sensible_latent_heat_fluxes_HWRF(HWRF_POM_oper,latg,long,timeg)
           
shtfl_hwrf20_pom, lhtfl_hwrf20_pom, time_hwrf20_pom = \
           get_sensible_latent_heat_fluxes_HWRF(HWRF_POM_exp,latg,long,timeg)
         
shtfl_hwrf20_hycom, lhtfl_hwrf20_hycom,time_hwrf20_hycom = \
          get_sensible_latent_heat_fluxes_HWRF(HWRF_HYCOM_exp,latg,long,timeg)

enth_hwrf19_pom_sg668 = shtfl_hwrf19_pom + lhtfl_hwrf19_pom
enth_hwrf20_pom_sg668 = shtfl_hwrf20_pom + lhtfl_hwrf20_pom
enth_hwrf20_hycom_sg668 = shtfl_hwrf20_hycom + lhtfl_hwrf20_hycom
dataset_id_sg668 = dataset_id

#%% Map detail

# Time bounds
#min_time = '2019-08-24T00:00:00Z'
#max_time = '2019-09-07T00:00:00Z'

# Time bounds
min_time = '2019/08/28/00'
max_time = '2019/08/29/00'
tmin = datetime.strptime(min_time,'%Y/%m/%d/%H')
tmax = datetime.strptime(max_time,'%Y/%m/%d/%H')

lat_lim = [15,23]
lon_lim = [-70,-63.5]

# Reading bathymetry data
ncbath = xr.open_dataset(bath_file)
bath_lat = ncbath.variables['lat'][:]
bath_lon = ncbath.variables['lon'][:]
bath_elev = ncbath.variables['elevation'][:]

oklatbath = np.logical_and(bath_lat >= lat_lim[0],bath_lat <= lat_lim[-1])
oklonbath = np.logical_and(bath_lon >= lon_lim[0],bath_lon <= lon_lim[-1])

bath_latsub = bath_lat[oklatbath]
bath_lonsub = bath_lon[oklonbath]
bath_elevs = bath_elev[oklatbath,:]
bath_elevsub = bath_elevs[:,oklonbath]

print('Looking for glider data sets')
e = ERDDAP(server = server)

# Grab every dataset available
datasets = pd.read_csv(e.get_search_url(response='csv', search_for='all'))

# Search constraints
kw = {
    'min_lon': lon_lim[0],
    'max_lon': lon_lim[1],
    'min_lat': lat_lim[0],
    'max_lat': lat_lim[1],
    'min_time': min_time,
    'max_time': max_time,
}

search_url = e.get_search_url(response='csv', **kw)

# Grab the results
search = pd.read_csv(search_url)

# Extract the IDs
gliders = search['Dataset ID'].values

msg = 'Found {} Glider Datasets:\n\n{}'.format
print(msg(len(gliders), '\n'.join(gliders)))

#% get entire deployment (lat and lon) during hurricane season

# Time bounds
min_time2 = '2019/06/01/00'
max_time2 = '2019/11/30/00'

# Search constraints
kw = {
    'min_lon': lon_lim[0],
    'max_lon': lon_lim[1],
    'min_lat': lat_lim[0],
    'max_lat': lat_lim[1],
    'min_time': min_time2,
    'max_time': max_time2,
}

search_url = e.get_search_url(response='csv', **kw)

# Grab the results
search = pd.read_csv(search_url)

# Extract the IDs
gliders_all = search['Dataset ID'].values

msg = 'Found {} Glider Datasets:\n\n{}'.format
print(msg(len(gliders_all), '\n'.join(gliders_all)))

# Setting constraints
constraints = {
        'time>=': min_time2,
        'time<=': max_time2,
        'latitude>=': lat_lim[0],
        'latitude<=': lat_lim[1],
        'longitude>=': lon_lim[0],
        'longitude<=': lon_lim[1],
        }

variables = [
        'depth',
        'latitude',
        'longitude',
        'time'
        ]

e = ERDDAP(
        server=server,
        protocol='tabledap',
        response='nc'
        )

lev = np.arange(-9000,9100,100)
fig, ax2 = plt.subplots(figsize=(10, 5))
plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
ax2.contourf(bath_lonsub,bath_latsub,bath_elevsub,lev,cmap=cmocean.cm.topo)
ax2.set_yticks([])
ax2.set_xticks([])

# Include in gliders oly the SGs
gliders = gliders[0:-1]

for i,id_all in enumerate(gliders_all):
    id = [id for id in gliders if id == id_all]
    if len(id) != 0:
        print(id[0])     
        e.dataset_id = id[0]
        e.constraints = constraints
        e.variables = variables
        
        df = e.to_pandas(
                index_col='time (UTC)',
                parse_dates=True,
                skiprows=(1,)  # units information can be dropped.
                ).dropna()
        
        print(len(df))
               
        timeg, ind = np.unique(df.index.values,return_index=True)
        latg = df['latitude (degrees_north)'].values[ind]
        long = df['longitude (degrees_east)'].values[ind]
        ax2.plot(long,latg,'.-',color='darkorange',markersize=1)
        
        idd = id[0].split('-')[0]
        if np.logical_or.reduce((idd =='SG663',idd =='SG664',idd =='SG667',idd =='SG668')):
            okt = np.where(mdates.date2num(timeg) >= \
                           mdates.date2num(datetime(2019,8,28,15)))[0][0]
        if idd =='SG666':
            okt = np.where(mdates.date2num(timeg) >= \
                           mdates.date2num(datetime(2019,8,29,6)))[0][0]
        if idd =='SG665':
            okt = np.where(mdates.date2num(timeg) >= \
                           mdates.date2num(datetime(2019,8,29,0)))[0][0]
        ax2.plot(long[okt],latg[okt],'*',color='k',markersize=6)
            #ax.plot(df['longitude (degrees_east)'][len(df['longitude (degrees_east)'])-1],\
            #     df['latitude (degrees_north)'][len(df['longitude (degrees_east)'])-1],\
            #     '-',color=col[i],\
            #     marker = mark[i],markeredgecolor = 'k',markersize=6,\
            #     label=id.split('-')[0])
        
        if idd =='SG667':
            ax2.text(long[okt],latg[okt]+0.3,id[0].split('-')[0],weight='bold',
                bbox=dict(facecolor='white',alpha=0.4,edgecolor='none'))
        else:
            ax2.text(long[okt]-1.3,latg[okt],id[0].split('-')[0],weight='bold',
                bbox=dict(facecolor='white',alpha=0.4,edgecolor='none'))

ax2.plot(lon_best_track[15],lat_best_track[15],'o',markersize=6,color='mediumturquoise',\
         markeredgecolor='k',label='ts')
ax2.plot(lon_best_track[18],lat_best_track[18],'o',markersize=6,color='yellow',\
         markeredgecolor='k',label='cat1')
#ax2.legend(loc='upper right',bbox_to_anchor=[1.1,1.0])
ax2.legend(loc='lower right')  
for i in np.arange(16,17):
    if cat[i] == 'ts':
        col = 'mediumturquoise'
        print(col)
    if cat[i] == 'cat1':
        col = 'yellow'
        print(col)
    ax2.plot(lon_best_track[i],lat_best_track[i],'o',markersize=8,color=col,\
         markeredgecolor='k')
    ax2.text(lon_best_track[i]-1,lat_best_track[i]-0.4,str(time_best_track[i])[5:13],\
             bbox=dict(facecolor='white',alpha=0.4,edgecolor='none'))         
for i in np.arange(17,23):
    if cat[i] == 'ts':
        col = 'mediumturquoise'
        print(col)
    if cat[i] == 'cat1':
        col = 'yellow'
        print(col)
    ax2.plot(lon_best_track[i],lat_best_track[i],'o',markersize=8,color=col,\
         markeredgecolor='k')
    ax2.text(lon_best_track[i]+0.2,lat_best_track[i]-0.2,str(time_best_track[i])[5:13],\
             bbox=dict(facecolor='white',alpha=0.4,edgecolor='none'))    
ax2.axis('scaled')
ax2.axis([lon_lim[0],lon_lim[1],lat_lim[0],lat_lim[1]])

#ax2.text(-69.5,22.5,'(b)',fontsize=16,color='white',weight='bold')
plt.subplots_adjust(left=0.05,right=0.95, wspace=-0.1, top=0.9)
#file = folder + 'Daily_map_North_Atlantic_gliders_in_DAC_detail_' + min_time[0:10] + '_' + max_time[0:10] + '.png'
file = folder_fig + 'Daily_map_North_Atlantic_gliders_in_DAC_detail3_' + str(tmin)[0:10] + '_' + str(tmax)[0:10] 
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)  

#%% Figure forecasted intensity models vs best intensity

okt = np.logical_and(time_best_track >= tini,time_best_track <= tend)

lead_time_pom_oper = np.arange(0,129,3)

fig,ax1 = plt.subplots(figsize=(10, 5))
plt.ion()
plt.plot(lead_time_pom_oper[::2],wind_int_kt[okt],'o-k',label='Best')
plt.plot(lead_time_pom_oper,max_wind_10m_hwrf_pom19_oper,'X-',color='mediumorchid',label='HWRF2010-POM (IC Clim.)',markeredgecolor='k',markersize=7)
plt.plot(lead_time_pom_oper,max_wind_10m_hwrf_pom20_exp,'^-',color='teal',label='HWRF2020-POM (IC RTOFS)',markeredgecolor='k',markersize=7)
plt.plot(lead_time_pom_oper,max_wind_10m_hwrf_hycom20_exp,'H-',color='darkorange',label='HWRF2020-HYCOM (IC RTOFS)',markeredgecolor='k',markersize=7)
plt.plot(np.tile(15,len(np.arange(0,100))),np.arange(0,100),'--',color='silver')
plt.plot(np.tile(24,len(np.arange(0,100))),np.arange(0,100),'--',color='grey')
plt.plot(np.tile(30,len(np.arange(0,100))),np.arange(0,100),'--k')

ax1.tick_params(which='major', width=2)
ax1.tick_params(which='major', length=7)
ax1.tick_params(which='minor', length=4, color='k')

ax1.xaxis.set_major_locator(MultipleLocator(12))
ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax1.xaxis.set_minor_locator(MultipleLocator(3))
ax1.xaxis.set_ticks(np.arange(0,126,12))
ax1.xaxis.set_ticklabels(['28-Aug \n 0','\n 12','29-Aug \n 24','\n 36','30-Aug \n 48',\
                          '\n 60','21-Aug \n 72','\n 84','01-Sep \n 96','\n 108','02-Sep \n 120'])
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
plt.legend(loc='upper left',fontsize=14)
plt.ylim([20,165])
plt.xlim([0,126])
plt.xticks(np.arange(0,126,12))
plt.title('Intensity Forecast Dorian '+ cycle,fontsize=18)
plt.ylabel('Max 10m Wind (kt)',fontsize=14)

ax2 = ax1.twinx()
plt.ylim([20,165])
yticks = [64,83,96,113,137]
plt.yticks(yticks,['Cat 1','Cat 2','Cat 3','Cat 4','Cat 5'])
plt.grid(True)

file = folder_fig + 'best_intensity_vs_forec_intensity_' + cycle 
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 

#%% Heat fluxes figure SG668

fig,ax = plt.subplots(figsize=(10,5))
plt.plot(time_hwrf19_pom,enth_hwrf19_pom_sg668,'X-',color='mediumorchid',label='HWRF2019-POM (IC Clim.)',\
         linewidth=2,markeredgecolor='k',markersize=7)
plt.plot(time_hwrf20_pom,enth_hwrf20_pom_sg668,'^-',color='teal',label='HWRF2020-POM (IC RTOFS)',\
         linewidth=2,markeredgecolor='k',markersize=7)
plt.plot(time_hwrf20_hycom,enth_hwrf20_hycom_sg668,'H-',color='darkorange',label='HWRF2020-HYCOM (IC RTOFS)',\
         linewidth=2,markeredgecolor='k',markersize=7)
#plt.plot(timeg[oktg],Qs_glider2,'s',color='royalblue',label='Sensible '+inst_id.split('-')[0],markeredgecolor='k',markersize=7)
#plt.plot(timeg[oktg],Ql_glider2,'o',color='royalblue',label='Latent '+inst_id.split('-')[0],markeredgecolor='k',markersize=7)
plt.legend()
plt.title('Enthalpy Fluxes on '+ dataset_id_sg668.split('-')[0] + ' location. Cycle ' + cycle,fontsize = 16)
plt.ylabel('($W/m^2$)',fontsize = 14)

ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')

xfmt = mdates.DateFormatter('%d-%b')
ax.xaxis.set_major_formatter(xfmt)
dt =  mdates.date2num(datetime(2019,8,28,3)) - mdates.date2num(datetime(2019,8,28))
ax.xaxis.set_minor_locator(MultipleLocator(dt))
plt.plot(np.tile(datetime(2019,8,28,15),len(np.arange(1000))),np.arange(1000),'--',color='silver')
plt.xlim([time_hwrf[0],time_hwrf[-1]])
plt.ylim([0,1000])

file_name = folder_fig + 'Dorian_enthalpy_' + dataset_id_sg668 + '_' + cycle
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1) 


#%% Heat fluxes figure SG665

fig,ax = plt.subplots(figsize=(10,5))
plt.plot(time_hwrf19_pom,enth_hwrf19_pom_sg665,'X-',color='mediumorchid',label='HWRF2019-POM (IC Clim.)',\
         linewidth=2,markeredgecolor='k',markersize=7)
plt.plot(time_hwrf20_pom,enth_hwrf20_pom_sg665,'^-',color='teal',label='HWRF2020-POM (IC RTOFS)',\
         linewidth=2,markeredgecolor='k',markersize=7)
plt.plot(time_hwrf20_hycom,enth_hwrf20_hycom_sg665,'H-',color='darkorange',label='HWRF2020-HYCOM (IC RTOFS)',\
         linewidth=2,markeredgecolor='k',markersize=7)
#plt.plot(timeg[oktg],Qs_glider2,'s',color='royalblue',label='Sensible '+inst_id.split('-')[0],markeredgecolor='k',markersize=7)
#plt.plot(timeg[oktg],Ql_glider2,'o',color='royalblue',label='Latent '+inst_id.split('-')[0],markeredgecolor='k',markersize=7)
plt.legend()
plt.title('Enthalpy Fluxes on '+ dataset_id_sg665.split('-')[0] + ' location. Cycle ' + cycle,fontsize = 16)
plt.ylabel('($W/m^2$)',fontsize = 14)

ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')

xfmt = mdates.DateFormatter('%d-%b')
ax.xaxis.set_major_formatter(xfmt)
dt =  mdates.date2num(datetime(2019,8,28,3)) - mdates.date2num(datetime(2019,8,28))
ax.xaxis.set_minor_locator(MultipleLocator(dt))
plt.plot(np.tile(datetime(2019,8,29,0),len(np.arange(1000))),np.arange(1000),'--',color='grey')
plt.xlim([time_hwrf[0],time_hwrf[-1]])
plt.ylim([0,1000])

file_name = folder_fig + 'Dorian_enthalpy_' + dataset_id_sg665 + '_' + cycle
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1) 
          
#%% Heat fluxes figure SG666

fig,ax = plt.subplots(figsize=(10,5))
plt.plot(time_hwrf19_pom,enth_hwrf19_pom_sg666,'X-',color='mediumorchid',label='HWRF2019-POM (IC Clim.)',\
         linewidth=2,markeredgecolor='k',markersize=7)
plt.plot(time_hwrf20_pom,enth_hwrf20_pom_sg666,'^-',color='teal',label='HWRF2020-POM (IC RTOFS)',\
         linewidth=2,markeredgecolor='k',markersize=7)
plt.plot(time_hwrf20_hycom,enth_hwrf20_hycom_sg666,'H-',color='darkorange',label='HWRF2020-HYCOM (IC RTOFS)',\
         linewidth=2,markeredgecolor='k',markersize=7)
#plt.plot(timeg[oktg],Qs_glider2,'s',color='royalblue',label='Sensible '+inst_id.split('-')[0],markeredgecolor='k',markersize=7)
#plt.plot(timeg[oktg],Ql_glider2,'o',color='royalblue',label='Latent '+inst_id.split('-')[0],markeredgecolor='k',markersize=7)
plt.legend()
plt.title('Enthalpy Fluxes on '+ dataset_id_sg666.split('-')[0] + ' location. Cycle ' + cycle,fontsize = 16)
plt.ylabel('($W/m^2$)',fontsize = 14)

ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')

xfmt = mdates.DateFormatter('%d-%b')
ax.xaxis.set_major_formatter(xfmt)
dt =  mdates.date2num(datetime(2019,8,28,3)) - mdates.date2num(datetime(2019,8,28))
ax.xaxis.set_minor_locator(MultipleLocator(dt))
plt.plot(np.tile(datetime(2019,8,29,6),len(np.arange(1000))),np.arange(1000),'--',color='k')
plt.xlim([time_hwrf[0],time_hwrf[-1]])
plt.ylim([0,1000])

file_name = folder_fig + 'Dorian_enthalpy_' + dataset_id_sg666 + '_' + cycle
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1) 

#%% Estimate bulk heat fluxes during glider deployments
'''          
url_NDBC = 'https://dods.ndbc.noaa.gov/thredds/dodsC/data/stdmet/41043/41043h2019.nc'

wind_NDBC = xr.open_dataset(url_NDBC,decode_times=False)

tt = wind_NDBC['time']
time_NDBC = netCDF4.num2date(tt[:],tt.units)
okt = np.logical_and(time_NDBC >= tini, time_NDBC < tend)

t_NDBC = time_NDBC[okt]
wspd_NDBC = np.asarray(wind_NDBC['wind_spd'][:])[okt,0,0] # winds measured at 3.8 meters height
wdir_NDBC = np.asarray(wind_NDBC['wind_dir'][:])[okt,0,0] 
Tair_NDBC = np.asarray(wind_NDBC['air_temperature'][:])[okt,0,0]
Tsea_NDBC = np.asarray(wind_NDBC['sea_surface_temperature'][:])[okt,0,0]

# Simple estimate from S. A. Hsu et al. (1994) (https://www.ndbc.noaa.gov/adjust_wind.shtml)
# Uz2 = Uz1 (z2/z1)^(0.11)
wspd_NDBC_10m = wspd_NDBC *(10/3.8)**(0.11)

# Sensible heat flux NDBC
rho = 1.184 # kg/m^3
cp = 1004 # J/(kg K)
cs = 1 * 10**(-3)
Tair_NDBC_to_glid = np.interp(mdates.date2num(timeg),mdates.date2num(t_NDBC),Tair_NDBC)
del_t = Tmean_dens_crit_glid - Tair_NDBC_to_glid
oktg = np.logical_and(timeg >= datetime(2019,8,29,6),timeg <= datetime(2019,8,29,6,30))
oktt = time_best_track == datetime(2019,8,29,6)
#Qs_glider = rho * cp * cs * wind_int_ms[oktt] * del_t[oktg]
Qs_glider2 = 111.3 #SG665
#Qs_glider2 = 42.1  #SG666

# Latent heat flux NDBC
rho = 1.184 # kg/m^3
Le = 2.5 * 10**(6) # J/(kg K)
cl = 1.2 * 10**(-3)
qs = 1
rh = 0.97
#Ql_glider = rho * Le * cl * wind_int_ms[oktt] * qs *(1 - rh)   
Ql_glider2 = 419.7 #SG665
#Ql_glider2 = 184.9 #SG666
'''
