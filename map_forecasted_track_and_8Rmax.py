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
folder_fig = '/www/web/rucool/aristizabal/Figures/'

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

Rmax_POM_oper = []
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
    Rmax_POM_oper.append(seawater.dist([lat_forec_track_pom_oper[indx],lat_maxwind],\
                         [lon_forec_track_pom_oper[indx],lon_maxwind],'km')[0][0])

lev = np.arange(-9000,9100,100)

fig,ax = plt.subplots()
plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,lev,cmap=cmocean.cm.topo)
plt.plot(lon_forec_track_pom_oper, lat_forec_track_pom_oper,'X-',color='mediumorchid',\
         markeredgecolor='k',label='HWRF2019-POM Oper (IC Clim.)',markersize=7)
plt.axis('scaled')
plt.xlim([np.min(lon_forec_track_pom_oper)-2,np.max(lon_forec_track_pom_oper)+2])
plt.ylim([np.min(lat_forec_track_pom_oper)-2,np.max(lat_forec_track_pom_oper)+2])
plt.title('Forecastes Track and 8Rmax \n Dorian cycle='+ cycle,fontsize=18)
plt.legend()
for x,R in enumerate(Rmax_POM_oper):
    print(x)
    circle = plt.Circle((lon_forec_track_pom_oper[x],lat_forec_track_pom_oper[x]),8*R/110,\
             color='k',fill=False)
    ax.add_artist(circle)
plt.plot(lon_forec_track_pom_oper, lat_forec_track_pom_oper,'X-',color='mediumorchid',\
         markeredgecolor='k',markersize=7)
file = folder_fig + 'best_track_and_8Rmax_HWRF2019_POM_oper_' + cycle
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)

#%% Get Rmax and R values POM experimental

Rmax_POM_exp = []
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
    Rmax_POM_exp.append(seawater.dist([lat_forec_track_pom_exp[indx],lat_maxwind],\
                         [lon_forec_track_pom_exp[indx],lon_maxwind],'km')[0][0])

lev = np.arange(-9000,9100,100)

fig,ax = plt.subplots()
plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,lev,cmap=cmocean.cm.topo)
plt.plot(lon_forec_track_pom_exp, lat_forec_track_pom_exp,'^-',color='teal',\
         markeredgecolor='k',label='HWRF2020-POM Exp (IC RTOFS)',markersize=7)
plt.axis('scaled')
plt.xlim([np.min(lon_forec_track_pom_oper)-2,np.max(lon_forec_track_pom_oper)+2])
plt.ylim([np.min(lat_forec_track_pom_oper)-2,np.max(lat_forec_track_pom_oper)+2])
plt.title('Forecastes Track and 8Rmax \n Dorian cycle='+ cycle,fontsize=18)
plt.legend()
for x,R in enumerate(Rmax_POM_oper):
    print(x)
    circle = plt.Circle((lon_forec_track_pom_oper[x],lat_forec_track_pom_oper[x]),8*R/110,\
             color='k',fill=False)
    ax.add_artist(circle)
plt.plot(lon_forec_track_pom_exp, lat_forec_track_pom_exp,'^-',color='teal',\
         markeredgecolor='k',markersize=7)
file = folder_fig + 'best_track_and_8Rmax_HWRF2020_POM_exp_' + cycle
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)

#%% Get Rmax and R values HYCOM experimental

Rmax_HYCOM_exp = []
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
    Rmax_HYCOM_exp.append(seawater.dist([lat_forec_track_hycom_exp[indx],lat_maxwind],\
                         [lon_forec_track_hycom_exp[indx],lon_maxwind],'km')[0][0])

lev = np.arange(-9000,9100,100)

fig,ax = plt.subplots()
plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,lev,cmap=cmocean.cm.topo)
plt.plot(lon_forec_track_hycom_exp, lat_forec_track_hycom_exp,'H-',color='orange',\
         markeredgecolor='k',label='HWRF2020-HYCOM Exp (IC RTOFS)',markersize=7)
plt.axis('scaled')
plt.xlim([np.min(lon_forec_track_pom_oper)-2,np.max(lon_forec_track_pom_oper)+2])
plt.ylim([np.min(lat_forec_track_pom_oper)-2,np.max(lat_forec_track_pom_oper)+2])
plt.title('Forecastes Track and 8Rmax \n Dorian cycle='+ cycle,fontsize=18)
plt.legend()
for x,R in enumerate(Rmax_POM_oper):
    print(x)
    circle = plt.Circle((lon_forec_track_pom_oper[x],lat_forec_track_pom_oper[x]),8*R/110,\
             color='k',fill=False)
    ax.add_artist(circle)

plt.plot(lon_forec_track_hycom_exp, lat_forec_track_hycom_exp,'H-',color='orange',\
         markeredgecolor='k',markersize=7)
file = folder_fig + 'best_track_and_8Rmax_HWRF2020_HYCOM_exp_' + cycle
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)

#%% Rmax vs latitude

plt.figure(figsize=(10,5))
plt.plot(lat_forec_track_pom_oper,Rmax_POM_oper,'X-',color='mediumorchid',\
         markeredgecolor='k',markersize=7,label='HWRF2019-POM Oper (IC Clim.)')
plt.plot(lat_forec_track_pom_exp,Rmax_POM_exp,'^-',color='teal',\
          markeredgecolor='k',markersize=7,label='HWRF2020-POM Exp (IC RTOFS)')
plt.plot(lat_forec_track_hycom_exp,Rmax_HYCOM_exp, 'H-',color='orange',\
         markeredgecolor='k',markersize=7,label='HWRF2020-HYCOM Exp (IC RTOFS)')
plt.ylabel('Rmax (km)',fontsize=14)
plt.xlabel('Latitude ($^o$)',fontsize=14)

##############################################
#%% All figures together

fig, ax = plt.subplots(figsize=(10, 8))
grid = plt.GridSpec(2, 3, wspace=0.1, hspace=0.2,left=0.05,right=0.95)

#########
ax1 = plt.subplot(grid[0, 0])
plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,lev,cmap=cmocean.cm.topo)
plt.plot(lon_forec_track_pom_oper, lat_forec_track_pom_oper,'X-',color='mediumorchid',\
         markeredgecolor='k',label='HWRF2019-POM Oper (IC Clim.)',markersize=7)
plt.axis('scaled')
plt.xlim([np.min(lon_forec_track_pom_oper)-2,np.max(lon_forec_track_pom_oper)+2])
plt.ylim([np.min(lat_forec_track_pom_oper)-2,np.max(lat_forec_track_pom_oper)+2])
for x,R in enumerate(Rmax_POM_oper):
    print(x)
    circle = plt.Circle((lon_forec_track_pom_oper[x],lat_forec_track_pom_oper[x]),8*R/110,\
             color='k',fill=False)
    ax1.add_artist(circle)
plt.plot(lon_forec_track_pom_oper, lat_forec_track_pom_oper,'X-',color='mediumorchid',\
         markeredgecolor='k',markersize=7)
plt.text(-63.5,30,'(a)',fontsize=16)

#########

ax2 = plt.subplot(grid[0, 1])
plt.title('Hurricane Dorian Forecast Cycle '+ cycle + '\n Forecasted Track and 8Rmax',fontsize=16)
plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,lev,cmap=cmocean.cm.topo)
plt.plot(lon_forec_track_pom_exp, lat_forec_track_pom_exp,'^-',color='teal',\
         markeredgecolor='k',label='HWRF2020-POM Exp (IC RTOFS)',markersize=7)
plt.axis('scaled')
plt.xlim([np.min(lon_forec_track_pom_oper)-2,np.max(lon_forec_track_pom_oper)+2])
plt.ylim([np.min(lat_forec_track_pom_oper)-2,np.max(lat_forec_track_pom_oper)+2])
ax2.set_yticklabels(' ')
for x,R in enumerate(Rmax_POM_exp):
    print(x)
    circle = plt.Circle((lon_forec_track_pom_exp[x],lat_forec_track_pom_exp[x]),8*R/110,\
             color='k',fill=False)
    ax2.add_artist(circle)
plt.plot(lon_forec_track_pom_exp, lat_forec_track_pom_exp,'^-',color='teal',\
         markeredgecolor='k',markersize=7)
plt.text(-63.5,30,'(b)',fontsize=16)

##########

ax3 = plt.subplot(grid[0, 2])
plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,lev,cmap=cmocean.cm.topo)
plt.plot(lon_forec_track_hycom_exp, lat_forec_track_hycom_exp,'H-',color='orange',\
         markeredgecolor='k',label='HWRF2020-HYCOM Exp (IC RTOFS)',markersize=7)
plt.axis('scaled')
plt.xlim([np.min(lon_forec_track_pom_oper)-2,np.max(lon_forec_track_pom_oper)+2])
plt.ylim([np.min(lat_forec_track_pom_oper)-2,np.max(lat_forec_track_pom_oper)+2])
ax3.set_yticklabels(' ')
for x,R in enumerate(Rmax_HYCOM_exp):
    print(x)
    circle = plt.Circle((lon_forec_track_hycom_exp[x],lat_forec_track_hycom_exp[x]),8*R/110,\
             color='k',fill=False)
    ax3.add_artist(circle)
plt.plot(lon_forec_track_hycom_exp, lat_forec_track_hycom_exp,'H-',color='orange',\
         markeredgecolor='k',markersize=7)
plt.text(-63.5,30,'(c)',fontsize=16)

############

ax4 = plt.subplot(grid[1, :])
plt.plot(lat_forec_track_pom_oper,Rmax_POM_oper,'X-',color='mediumorchid',\
         markeredgecolor='k',markersize=7,label='HWRF2019-POM Oper (IC Clim.)')
plt.plot(lat_forec_track_pom_exp,Rmax_POM_exp,'^-',color='teal',\
          markeredgecolor='k',markersize=7,label='HWRF2020-POM Exp (IC RTOFS)')
plt.plot(lat_forec_track_hycom_exp,Rmax_HYCOM_exp, 'H-',color='orange',\
         markeredgecolor='k',markersize=7,label='HWRF2020-HYCOM Exp (IC RTOFS)')
plt.ylabel('Rmax (km)',fontsize=14)
plt.xlabel('Latitude ($^o$)',fontsize=14)
plt.legend()
plt.title('Radius of Maximum Wind',fontsize=16)
plt.ylim([0,42])
plt.text(15.5,38,'(d)',fontsize=16)

file = folder_fig + 'best_track_and_8Rmax_All_models_Dorian_' + cycle
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)
