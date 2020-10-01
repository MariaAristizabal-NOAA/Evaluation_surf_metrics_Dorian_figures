#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 10:01:02 2020

@author: aristizabal
"""

#%% User input

dir_home = '/home'
#dir_home = '/Volumes'

lon_lim = [-98.5,-60.0]
lat_lim = [10.0,45.0]

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
#bath_file = '/Users/aristizabal/Desktop/MARACOOS_project/Maria_scripts/nc_files/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'
bath_file = dir_home +'/aristizabal/bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'

# KMZ file best track Dorian
#kmz_file_Dorian = '/Users/aristizabal/Desktop/MARACOOS_project/Maria_scripts/KMZ_files/al052019_best_track-5.kmz'
kmz_file_Dorian = dir_home + '/aristizabal/KMZ_files/al052019_best_track-5.kmz'

# url for GOFS 3.1
url_GOFS = 'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_93.0/ts3z'

# figures
#folder_fig = '/Users/aristizabal/Desktop/MARACOOS_project/Maria_scripts/Figures/Model_glider_comp2/'
folder_fig = dir_home + '/aristizabal/Figures/'

# folder nc files POM
folder_pom19 =  dir_home + '/aristizabal/HWRF2019_POM_Dorian/'
folder_pom20 =  dir_home + '/aristizabal/HWRF2020_POM_Dorian/'

# folde HWRF2020_HYCOM
folder_hycom20 = dir_home + '/aristizabal/HWRF2020_HYCOM_Dorian/'

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
Dir_HMON_HYCOM = '/Volumes/aristizabal/ncep_model/HWRF-Hycom-WW3_exp_Michael/'
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
import os
import os.path
import glob
import cmocean
import seawater

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

#%% Get storm track from HWRF/POM output

def get_storm_track_POM(file_track):

    ff = open(file_track,'r')
    f = ff.readlines()

    latt = []
    lont = []
    lead_time = []
    for l in f:
        lat = float(l.split(',')[6][0:4])/10
        if l.split(',')[6][4] == 'N':
            lat = lat
        else:
            lat = -lat
        lon = float(l.split(',')[7][0:5])/10
        if l.split(',')[7][4] == 'E':
            lon = lon
        else:
            lon = -lon
        latt.append(lat)
        lont.append(lon)
        lead_time.append(int(l.split(',')[5][1:4]))

    latt = np.asarray(latt)
    lont = np.asarray(lont)
    lead_time, ind = np.unique(lead_time,return_index=True)
    lat_track = latt[ind]
    lon_track = lont[ind]

    return lon_track, lat_track, lead_time

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

#%% Get Dorian track from models

lon_forec_track_pom_oper, lat_forec_track_pom_oper, lead_time_pom_oper = get_storm_track_POM(hwrf_pom_track_oper)

lon_forec_track_pom_exp, lat_forec_track_pom_exp, lead_time_pom_exp = get_storm_track_POM(hwrf_pom_track_exp)

lon_forec_track_hycom_exp, lat_forec_track_hycom_exp, lead_time_hycom_exp = get_storm_track_POM(hwrf_hycom_track_exp)

#%% Get list HWRF files

HWRF_POM_oper = sorted(glob.glob(os.path.join(folder_hwrf_pom19_oper,'*.nc')))
HWRF_POM_exp = sorted(glob.glob(os.path.join(folder_hwrf_pom20_exp,'*.nc')))
HWRF_HYCOM_exp = sorted(glob.glob(os.path.join(folder_hwrf_hycom20_exp,'*.nc')))

#%% Get Rmax and R values POM operational

lev = np.arange(-9000,9100,100)

fig,ax = plt.subplots()
plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,lev,cmap=cmocean.cm.topo)
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='seashell')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.plot(lon_forec_track_pom_oper, lat_forec_track_pom_oper,'X-',color='mediumorchid',\
         markeredgecolor='k',label='HWRF2019-POM Oper (IC Clim.)',markersize=7)
plt.axis('scaled')
plt.xlim([np.min(lon_forec_track_pom_oper)-2,np.max(lon_forec_track_pom_oper)+2])
plt.ylim([np.min(lat_forec_track_pom_oper)-2,np.max(lat_forec_track_pom_oper)+2])
plt.title('Forecastes Track and 8Rmax \n Dorian cycle='+ cycle,fontsize=18)
plt.legend()

for indx in np.arange(len(HWRF_POM_oper)):    
    print(HWRF_POM_oper[indx])
    HWRF = xr.open_dataset(HWRF_POM_oper[indx])
    lat_hwrf = np.asarray(HWRF.variables['latitude'][:])
    lon_hwrf = np.asarray(HWRF.variables['longitude'][:])
    time_hwrf = np.asarray(HWRF.variables['time'][:])
    UGRD_hwrf = np.asarray(HWRF.variables['UGRD_10maboveground'][0,:,:])
    VGRD_hwrf = np.asarray(HWRF.variables['VGRD_10maboveground'][0,:,:])
    
    wind_int = np.sqrt(UGRD_hwrf**2 + VGRD_hwrf**2)
    max_wind= np.max(wind_int)
    okwind = np.where(wind_int == max_wind)
    lat_maxwind = lat_hwrf[okwind[0][0]]
    lon_maxwind = lon_hwrf[okwind[1][0]]
    Rmax = seawater.dist([lat_forec_track_pom_oper[indx],lat_maxwind],\
                         [lon_forec_track_pom_oper[indx],lon_maxwind],'km')[0][0]
        
    xlim = [lon_forec_track_pom_oper[indx]-4,lon_forec_track_pom_oper[indx]+4]
    ylim = [lat_forec_track_pom_oper[indx]-4,lat_forec_track_pom_oper[indx]+4]
    
    oklon = np.where(np.logical_and(lon_hwrf>xlim[0],lon_hwrf<xlim[1]))[0]
    oklat = np.where(np.logical_and(lat_hwrf>ylim[0],lat_hwrf<ylim[1]))[0]
    
    meshlon_lat = np.meshgrid(lon_hwrf[oklon],lat_hwrf[oklat])
    
    lat_lon_matrix = np.stack((np.ravel(meshlon_lat[0]),np.ravel(meshlon_lat[1])),axis=1).T
        
    R = np.empty(lat_lon_matrix.shape[1])
    R[:] = np.nan
    for i in np.arange(lat_lon_matrix.shape[1]):
        R[i] = seawater.dist([lat_forec_track_pom_oper[indx],lat_lon_matrix[1,i]],\
                             [lon_forec_track_pom_oper[indx],lat_lon_matrix[0,i]],'km')[0][0] 
    
    R_norm = R/Rmax
    okR8 = np.logical_and(R_norm <= 8.05,R_norm >= 7.95)

    plt.plot(lat_lon_matrix[0][okR8],lat_lon_matrix[1][okR8],'.k',markersize=1)    

plt.plot(lon_forec_track_pom_oper, lat_forec_track_pom_oper,'X-',color='mediumorchid',\
         markeredgecolor='k',markersize=7)
file = folder_fig + 'best_track_and_8Rmax_HWRF2019_POM_oper_' + cycle
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)

#%% Get Rmax and R values POM experimental

lev = np.arange(-9000,9100,100)

fig,ax = plt.subplots()
plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,lev,cmap=cmocean.cm.topo)
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='seashell')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.plot(lon_forec_track_pom_exp, lat_forec_track_pom_exp,'^-',color='teal',\
         markeredgecolor='k',label='HWRF2020-POM Exp (IC RTOFS)',markersize=7)
plt.axis('scaled')
plt.xlim([np.min(lon_forec_track_pom_oper)-2,np.max(lon_forec_track_pom_oper)+2])
plt.ylim([np.min(lat_forec_track_pom_oper)-2,np.max(lat_forec_track_pom_oper)+2])
plt.title('Forecastes Track and 8Rmax \n Dorian cycle='+ cycle,fontsize=18)
plt.legend()

for indx in np.arange(len(HWRF_POM_exp)):
    print(HWRF_POM_exp[indx])
    HWRF = xr.open_dataset(HWRF_POM_exp[indx])
    lat_hwrf = np.asarray(HWRF.variables['latitude'][:])
    lon_hwrf = np.asarray(HWRF.variables['longitude'][:])
    time_hwrf = np.asarray(HWRF.variables['time'][:])
    UGRD_hwrf = np.asarray(HWRF.variables['UGRD_10maboveground'][0,:,:])
    VGRD_hwrf = np.asarray(HWRF.variables['VGRD_10maboveground'][0,:,:])
    
    wind_int = np.sqrt(UGRD_hwrf**2 + VGRD_hwrf**2)
    max_wind= np.max(wind_int)
    okwind = np.where(wind_int == max_wind)
    lat_maxwind = lat_hwrf[okwind[0][0]]
    lon_maxwind = lon_hwrf[okwind[1][0]]
    Rmax = seawater.dist([lat_forec_track_pom_exp[indx],lat_maxwind],\
                         [lon_forec_track_pom_exp[indx],lon_maxwind],'km')[0][0]
        
    xlim = [lon_forec_track_pom_exp[indx]-4,lon_forec_track_pom_exp[indx]+4]
    ylim = [lat_forec_track_pom_exp[indx]-4,lat_forec_track_pom_exp[indx]+4]
    
    oklon = np.where(np.logical_and(lon_hwrf>xlim[0],lon_hwrf<xlim[1]))[0]
    oklat = np.where(np.logical_and(lat_hwrf>ylim[0],lat_hwrf<ylim[1]))[0]
    
    meshlon_lat = np.meshgrid(lon_hwrf[oklon],lat_hwrf[oklat])
    
    lat_lon_matrix = np.stack((np.ravel(meshlon_lat[0]),np.ravel(meshlon_lat[1])),axis=1).T
    
    R = np.empty(lat_lon_matrix.shape[1])
    R[:] = np.nan
    for i in np.arange(lat_lon_matrix.shape[1]):
        R[i] = seawater.dist([lat_forec_track_pom_exp[indx],lat_lon_matrix[1,i]],\
                             [lon_forec_track_pom_exp[indx],lat_lon_matrix[0,i]],'km')[0][0] 
    
    R_norm = R/Rmax
    okR8 = np.logical_and(R_norm <= 8.05,R_norm >= 7.95)
    
    plt.plot(lat_lon_matrix[0][okR8],lat_lon_matrix[1][okR8],'.k',markersize=1)    

plt.plot(lon_forec_track_pom_exp, lat_forec_track_pom_exp,'^-',color='teal',\
         markeredgecolor='k',markersize=7)
file = folder_fig + 'best_track_and_8Rmax_HWRF2020_POM_exp_' + cycle
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)

#%% Get Rmax and R values HYCOM experimental

lev = np.arange(-9000,9100,100)

fig,ax = plt.subplots()
plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,lev,cmap=cmocean.cm.topo)
plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='seashell')
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.plot(lon_forec_track_hycom_exp, lat_forec_track_hycom_exp,'H-',color='orange',\
         markeredgecolor='k',label='HWRF2020-HYCOM Exp (IC RTOFS)',markersize=7)
plt.axis('scaled')
plt.xlim([np.min(lon_forec_track_pom_oper)-2,np.max(lon_forec_track_pom_oper)+2])
plt.ylim([np.min(lat_forec_track_pom_oper)-2,np.max(lat_forec_track_pom_oper)+2])
plt.title('Forecastes Track and 8Rmax \n Dorian cycle='+ cycle,fontsize=18)
plt.legend()

for indx in np.arange(len(HWRF_POM_exp)):
    print(HWRF_HYCOM_exp[indx])
    HWRF = xr.open_dataset(HWRF_HYCOM_exp[indx])
    lat_hwrf = np.asarray(HWRF.variables['latitude'][:])
    lon_hwrf = np.asarray(HWRF.variables['longitude'][:])
    time_hwrf = np.asarray(HWRF.variables['time'][:])
    UGRD_hwrf = np.asarray(HWRF.variables['UGRD_10maboveground'][0,:,:])
    VGRD_hwrf = np.asarray(HWRF.variables['VGRD_10maboveground'][0,:,:])
    
    wind_int = np.sqrt(UGRD_hwrf**2 + VGRD_hwrf**2)
    max_wind= np.max(wind_int)
    okwind = np.where(wind_int == max_wind)
    lat_maxwind = lat_hwrf[okwind[0][0]]
    lon_maxwind = lon_hwrf[okwind[1][0]]
    Rmax = seawater.dist([lat_forec_track_hycom_exp[indx],lat_maxwind],\
                         [lon_forec_track_hycom_exp[indx],lon_maxwind],'km')[0][0]
        
    xlim = [lon_forec_track_hycom_exp[indx]-4,lon_forec_track_hycom_exp[indx]+4]
    ylim = [lat_forec_track_hycom_exp[indx]-4,lat_forec_track_hycom_exp[indx]+4]
    
    oklon = np.where(np.logical_and(lon_hwrf>xlim[0],lon_hwrf<xlim[1]))[0]
    oklat = np.where(np.logical_and(lat_hwrf>ylim[0],lat_hwrf<ylim[1]))[0]
    
    meshlon_lat = np.meshgrid(lon_hwrf[oklon],lat_hwrf[oklat])
    
    lat_lon_matrix = np.stack((np.ravel(meshlon_lat[0]),np.ravel(meshlon_lat[1])),axis=1).T
    
    R = np.empty(lat_lon_matrix.shape[1])
    R[:] = np.nan
    for i in np.arange(lat_lon_matrix.shape[1]):
        R[i] = seawater.dist([lat_forec_track_hycom_exp[indx],lat_lon_matrix[1,i]],\
                             [lon_forec_track_hycom_exp[indx],lat_lon_matrix[0,i]],'km')[0][0] 
    
    R_norm = R/Rmax
    okR8 = np.logical_and(R_norm <= 8.05,R_norm >= 7.95)
    
    plt.plot(lat_lon_matrix[0][okR8],lat_lon_matrix[1][okR8],'.k',markersize=1)    

plt.plot(lon_forec_track_hycom_exp, lat_forec_track_hycom_exp,'H-',color='orange',\
         markeredgecolor='k',markersize=7)
file = folder_fig + 'best_track_and_8Rmax_HWRF2020_HYCOM_exp_' + cycle
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)