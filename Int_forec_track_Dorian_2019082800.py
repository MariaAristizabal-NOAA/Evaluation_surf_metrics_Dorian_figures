#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 14:01:57 2020

@author: root
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
kmz_file_Dorian = dir_home + '/aristizabal/KMZ_files/2019/al052019_best_track.kmz'

# url for GOFS 3.1
url_GOFS = 'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_93.0/ts3z'

# figures
#folder_fig = '/Users/aristizabal/Desktop/MARACOOS_project/Maria_scripts/Figures/Model_glider_comp2/'
#folder_fig = dir_home + '/aristizabal/Figures/'
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
from datetime import datetime, timedelta
import os
import os.path
from bs4 import BeautifulSoup
from zipfile import ZipFile
import glob
import cmocean
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)

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

    return lon_best_track, lat_best_track, time_best_track, wind_int_mph, wind_int_kt, cat

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

#%% Get Dorian best track

lon_best_track, lat_best_track, time_best_track, _, wind_int_kt, _ = \
read_kmz_file_storm_best_track(kmz_file_Dorian)

#%% Obtain lat and lon band around forecated track operational

dlon = 0.1
nlevels = int(2*delta_lon /dlon) + 1

lon_bnd = np.linspace(lon_forec_track_pom_oper[2*Nini:2*Nend-1]-delta_lon,lon_forec_track_pom_oper[2*Nini:2*Nend-1]+delta_lon,nlevels)
lon_band = lon_bnd.ravel()
lat_bd = np.tile(lat_forec_track_pom_oper[2*Nini:2*Nend-1],lon_bnd.shape[0])
lat_bnd = lat_bd.reshape(lon_bnd.shape[0],lon_bnd.shape[1])
lat_band = lat_bnd.ravel()

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

#%% Figure forecasted track models vs best track

#Time window
date_ini = cycle[0:4]+'-'+cycle[4:6]+'-'+cycle[6:8]+' '+cycle[8:]+':00:00'
tini = datetime.strptime(date_ini,'%Y-%m-%d %H:%M:%S')
tend = tini + timedelta(hours=float(lead_time_pom_oper[-1]))
date_end = str(tend)

okt = np.logical_and(time_best_track >= tini,time_best_track <= tend)

# time forecasted track_exp
time_forec_track_oper = np.asarray([tini + timedelta(hours = float(t)) for t in lead_time_pom_oper])
oktt = [np.where(t == time_forec_track_oper)[0][0] for t in time_best_track[okt]]

str_time = [str(tt)[5:13] for tt in time_forec_track_oper[oktt]]

lev = np.arange(-9000,9100,100)

fig,ax = plt.subplots()
plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,lev,cmap=cmocean.cm.topo)
#plt.contourf(bath_lon,bath_lat,bath_elev,[0,10000],colors='seashell')
#plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')
plt.plot(lon_forec_track_pom_oper[oktt], lat_forec_track_pom_oper[oktt],'X-',color='mediumorchid',\
         markeredgecolor='k',label='POM Oper',markersize=7)
plt.plot(lon_forec_track_pom_exp[oktt], lat_forec_track_pom_exp[oktt],'^-',color='teal',\
         markeredgecolor='k',label='POM Exp',markersize=7)
plt.plot(lon_forec_track_hycom_exp[oktt], lat_forec_track_hycom_exp[oktt],'H-',color='orange',\
         markeredgecolor='k',label='HYCOM Exp',markersize=7)
plt.plot(lon_best_track[okt], lat_best_track[okt],'o-',color='k',label='Best Track')
#plt.legend()
plt.title('Track Forecast Dorian cycle '+ cycle,fontsize=18)
plt.axis('scaled')
plt.xlim([np.min(lon_forec_track_pom_oper[oktt])-0.5,np.max(lon_forec_track_pom_oper[oktt])+0.5])
plt.ylim([np.min(lat_forec_track_pom_oper[oktt])-0.5,np.max(lat_forec_track_pom_oper[oktt])+0.5])

for i,t in enumerate(lead_time_pom_oper[oktt][::2][0:6]):
    ax.text(lon_forec_track_pom_oper[oktt][::2][i]-1,lat_forec_track_pom_oper[oktt][::2][i]-0.3,\
            str(t),fontsize=10,bbox=dict(facecolor='white', alpha=0.5))

for i,t in enumerate(lead_time_pom_oper[oktt][::2][6:10]):
    ax.text(lon_forec_track_pom_oper[oktt][::2][6:10][i]+0.5,lat_forec_track_pom_oper[oktt][::2][6:10][i],\
            str(t),fontsize=10,bbox=dict(facecolor='white', alpha=0.3))

plt.text(-63.9,19.2,'T1',fontsize=14)
plt.arrow(-65.4,17.7,1.3,1.2,fc='k',head_width=0.3, head_length=0.3)
plt.text(-69.8,26.5,'T2',fontsize=14)
plt.arrow(-71.3,25.1,1.3,1.2,fc='k',head_width=0.3, head_length=0.3)
plt.text(-72.2,28.3,'T3',fontsize=14)
plt.arrow(-73.7,26.75,1.3,1.2,fc='k',head_width=0.3, head_length=0.3)

file = folder_fig + 'best_track_vs_forec_track_' + cycle
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)

#%% Figure forecasted intensity models vs best intensity

okt = np.logical_and(time_best_track >= tini,time_best_track <= tend)

lead_time_pom_oper = np.arange(0,129,3)

fig,ax1 = plt.subplots(figsize=(10, 4))
plt.plot(lead_time_pom_oper[::2],wind_int_kt[okt],'o-k',label='Best')
plt.plot(lead_time_pom_oper,max_wind_10m_hwrf_pom19_oper,'X-',color='mediumorchid',label='HWRF2010-POM (IC Clim.)',markeredgecolor='k',markersize=7)
plt.plot(lead_time_pom_oper,max_wind_10m_hwrf_pom20_exp,'^-',color='teal',label='HWRF2020-POM (IC RTOFS)',markeredgecolor='k',markersize=7)
plt.plot(lead_time_pom_oper,max_wind_10m_hwrf_hycom20_exp,'H-',color='darkorange',label='HWRF2020-HYCOM (IC RTOFS)',markeredgecolor='k',markersize=7)
plt.plot(np.tile(18,len(np.arange(0,110))),np.arange(0,110),'--',color='k')
plt.plot(np.tile(66,len(np.arange(0,110))),np.arange(0,110),'--',color='k')
plt.plot(np.tile(84,len(np.arange(0,110))),np.arange(0,110),'--k')
plt.text(14,25,'T1',fontsize=14)
plt.text(62,25,'T2',fontsize=14)
plt.text(80,25,'T3',fontsize=14)

ax1.tick_params(which='major', width=2)
ax1.tick_params(which='major', length=7)
ax1.tick_params(which='minor', length=4, color='k')

ax1.xaxis.set_major_locator(MultipleLocator(12))
ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax1.xaxis.set_minor_locator(MultipleLocator(3))
ax1.xaxis.set_ticks(np.arange(0,126,12))
ax1.xaxis.set_ticklabels(['28-Aug \n 0','\n 12','29-Aug \n 24','\n 36','30-Aug \n 48',\
                          '\n 60','31-Aug \n 72','\n 84','01-Sep \n 96','\n 108','02-Sep \n 120'])
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
plt.legend(loc='upper left',fontsize=14)
plt.ylim([20,165])
plt.xlim([0,126])
plt.xticks(np.arange(0,126,12))
plt.title('Intensity Forecast Dorian Cycle '+ cycle,fontsize=18)
plt.ylabel('Max 10m Wind (kt)',fontsize=14)

ax2 = ax1.twinx()
plt.ylim([20,165])
yticks = [64,83,96,113,137]
plt.yticks(yticks,['Cat 1','Cat 2','Cat 3','Cat 4','Cat 5'])
plt.grid(True)

file = folder_fig + 'best_intensity_vs_forec_intensity_' + cycle
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)

#%% Intensity error

int_err_hwrf19_pom_oper = (wind_int_kt[okt] - max_wind_10m_hwrf_pom19_oper[::2]) #*100/wind_int_kt[okt]
int_err_hwrf20_pom_exp = (wind_int_kt[okt] - max_wind_10m_hwrf_pom20_exp[::2]) #*100/wind_int_kt[okt]
int_err_hwrf20_hycom_exp = (wind_int_kt[okt] - max_wind_10m_hwrf_hycom20_exp[::2]) #*100/wind_int_kt[okt]

fig,ax1 = plt.subplots(figsize=(10, 4))
plt.plot(lead_time_pom_oper[::2],int_err_hwrf19_pom_oper,'X-',color='mediumorchid',label='HWRF2019-POM Oper',markeredgecolor='k',markersize=7)
plt.plot(lead_time_pom_exp[::2],int_err_hwrf20_pom_exp,'^-',color='teal',label='HRWF2020-POM Exp',markeredgecolor='k',markersize=7)
plt.plot(lead_time_hycom_exp[::2],int_err_hwrf20_hycom_exp,'H-',color='darkorange',label='HWRF2020-HYCOM Exp',markeredgecolor='k',markersize=7)
plt.plot(lead_time_pom_oper[::2],np.tile(0,len(lead_time_pom_oper[::2])),'--k')

plt.plot(np.tile(18,len(np.arange(-10,40))),np.arange(-10,40),'--',color='k')
plt.plot(np.tile(66,len(np.arange(-10,40))),np.arange(-10,40),'--',color='k')
plt.plot(np.tile(84,len(np.arange(-10,40))),np.arange(-10,40),'--k')
plt.text(14,25,'T1',fontsize=14)
plt.text(62,25,'T2',fontsize=14)
plt.text(80,25,'T3',fontsize=14)

ax1.tick_params(which='major', width=2)
ax1.tick_params(which='major', length=7)
ax1.tick_params(which='minor', length=4, color='k')

ax1.xaxis.set_major_locator(MultipleLocator(12))
ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax1.xaxis.set_minor_locator(MultipleLocator(3))
ax1.xaxis.set_ticks(np.arange(0,126,12))
ax1.xaxis.set_ticklabels(['28-Aug \n 0','\n 12','29-Aug \n 24','\n 36','30-Aug \n 48',\
                          '\n 60','31-Aug \n 72','\n 84','01-Sep \n 96','\n 108','02-Sep \n 120'])
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
plt.legend(loc='upper left',fontsize=14)
#plt.ylim([20,165])
plt.xlim([0,126])
plt.xticks(np.arange(0,126,12))

plt.title('Intensity Forecast Error Dorian Cycle '+ cycle,fontsize=18)
plt.ylabel('Forecast Error (Kt)',fontsize=14)
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14)
#plt.legend()

file = folder_fig + 'intensity_error_' + cycle
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)

#%% Three figures combined

fig, ax = plt.subplots(figsize=(15, 7))
grid = plt.GridSpec(2, 5, wspace=0.6, hspace=0.3,left=0.05,right=0.95)
plt.suptitle('Hurricane Dorian Forecast cycle '+ cycle,fontsize=18)

ax = plt.subplot(grid[0:2, 0:2])
plt.contourf(bath_lonsub,bath_latsub,bath_elevsub,lev,cmap=cmocean.cm.topo)
plt.plot(lon_forec_track_pom_oper[oktt], lat_forec_track_pom_oper[oktt],'X-',color='mediumorchid',\
         markeredgecolor='k',label='POM Oper',markersize=7)
plt.plot(lon_forec_track_pom_exp[oktt], lat_forec_track_pom_exp[oktt],'^-',color='teal',\
         markeredgecolor='k',label='POM Exp',markersize=7)
plt.plot(lon_forec_track_hycom_exp[oktt], lat_forec_track_hycom_exp[oktt],'H-',color='orange',\
         markeredgecolor='k',label='HYCOM Exp',markersize=7)
plt.plot(lon_best_track[okt], lat_best_track[okt],'o-',color='k',label='Best Track')
plt.title('Track Forecast',fontsize=16)
plt.axis('scaled')
plt.xlim([np.min(lon_forec_track_pom_oper[oktt])-0.5,np.max(lon_forec_track_pom_oper[oktt])+0.5])
plt.ylim([np.min(lat_forec_track_pom_oper[oktt])-0.5,np.max(lat_forec_track_pom_oper[oktt])+0.5])

for i,t in enumerate(lead_time_pom_oper[oktt][::2][0:6]):
    ax.text(lon_forec_track_pom_oper[oktt][::2][i]-1,lat_forec_track_pom_oper[oktt][::2][i]-0.3,\
            str(t),fontsize=10,bbox=dict(facecolor='white', alpha=0.5))

for i,t in enumerate(lead_time_pom_oper[oktt][::2][6:10]):
    ax.text(lon_forec_track_pom_oper[oktt][::2][6:10][i]+0.5,lat_forec_track_pom_oper[oktt][::2][6:10][i],\
            str(t),fontsize=10,bbox=dict(facecolor='white', alpha=0.3))

plt.text(-63.9,19.2,'T1',fontsize=14)
plt.arrow(-65.4,17.7,1.3,1.2,fc='k',head_width=0.3, head_length=0.3)
plt.text(-69.8,26.5,'T2',fontsize=14)
plt.arrow(-71.3,25.1,1.3,1.2,fc='k',head_width=0.3, head_length=0.3)
plt.text(-72.2,28.3,'T3',fontsize=14)
plt.arrow(-73.7,26.75,1.3,1.2,fc='k',head_width=0.3, head_length=0.3)
plt.text(-81.5,31,'(a)',fontsize=16)

ax1 = plt.subplot(grid[0, 2:])
plt.plot(lead_time_pom_oper[::2],wind_int_kt[okt],'o-k',label='Best')
plt.plot(lead_time_pom_oper,max_wind_10m_hwrf_pom19_oper,'X-',color='mediumorchid',label='HWRF2010-POM (IC Clim.)',markeredgecolor='k',markersize=7)
plt.plot(lead_time_pom_oper,max_wind_10m_hwrf_pom20_exp,'^-',color='teal',label='HWRF2020-POM (IC RTOFS)',markeredgecolor='k',markersize=7)
plt.plot(lead_time_pom_oper,max_wind_10m_hwrf_hycom20_exp,'H-',color='darkorange',label='HWRF2020-HYCOM (IC RTOFS)',markeredgecolor='k',markersize=7)
plt.plot(np.tile(18,len(np.arange(0,110))),np.arange(0,110),'--',color='k')
plt.plot(np.tile(66,len(np.arange(0,110))),np.arange(0,110),'--',color='k')
plt.plot(np.tile(84,len(np.arange(0,110))),np.arange(0,110),'--k')
plt.text(14,25,'T1',fontsize=14)
plt.text(62,25,'T2',fontsize=14)
plt.text(80,25,'T3',fontsize=14)

ax1.tick_params(which='major', width=2)
ax1.tick_params(which='major', length=7)
ax1.tick_params(which='minor', length=4, color='k')

ax1.xaxis.set_major_locator(MultipleLocator(12))
ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax1.xaxis.set_minor_locator(MultipleLocator(3))
ax1.xaxis.set_ticks(np.arange(0,126,12))
#ax1.xaxis.set_ticklabels(['28-Aug \n 0','\n 12','29-Aug \n 24','\n 36','30-Aug \n 48',\
#                          '\n 60','31-Aug \n 72','\n 84','01-Sep \n 96','\n 108','02-Sep \n 120'])
ax1.xaxis.set_ticklabels([])
#plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
legend = plt.legend(loc='upper left',fontsize=12)
legend.get_frame().set_facecolor('white')
plt.ylim([20,165])
plt.xlim([0,126])
plt.xticks(np.arange(0,126,12))
plt.title('Intensity Forecast',fontsize=16)
plt.ylabel('Max 10m Wind (kt)',fontsize=14)

ax2 = ax1.twinx()
plt.ylim([20,165])
yticks = [64,83,96,113,137]
plt.yticks(yticks,['Cat 1','Cat 2','Cat 3','Cat 4','Cat 5'])
plt.grid(True)
ax1.legend(loc='upper left',fontsize=12)
plt.text(0,170,'(b)',fontsize=16)

ax1 = plt.subplot(grid[1, 2:])
plt.plot(lead_time_pom_oper[::2],int_err_hwrf19_pom_oper,'X-',color='mediumorchid',label='HWRF2019-POM Oper',markeredgecolor='k',markersize=7)
plt.plot(lead_time_pom_exp[::2],int_err_hwrf20_pom_exp,'^-',color='teal',label='HRWF2020-POM Exp',markeredgecolor='k',markersize=7)
plt.plot(lead_time_hycom_exp[::2],int_err_hwrf20_hycom_exp,'H-',color='darkorange',label='HWRF2020-HYCOM Exp',markeredgecolor='k',markersize=7)
plt.plot(lead_time_pom_oper[::2],np.tile(0,len(lead_time_pom_oper[::2])),'--k')

plt.plot(np.tile(18,len(np.arange(-10,50))),np.arange(-10,50),'--',color='k')
plt.plot(np.tile(66,len(np.arange(-10,50))),np.arange(-10,50),'--',color='k')
plt.plot(np.tile(84,len(np.arange(-10,50))),np.arange(-10,50),'--k')
plt.text(14,40,'T1',fontsize=14)
plt.text(62,40,'T2',fontsize=14)
plt.text(80,40,'T3',fontsize=14)

ax1.tick_params(which='major', width=2)
ax1.tick_params(which='major', length=7)
ax1.tick_params(which='minor', length=4, color='k')

ax1.xaxis.set_major_locator(MultipleLocator(12))
ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax1.xaxis.set_minor_locator(MultipleLocator(3))
ax1.xaxis.set_ticks(np.arange(0,126,12))
ax1.xaxis.set_ticklabels(['28-Aug \n 0','\n 12','29-Aug \n 24','\n 36','30-Aug \n 48',\
                          '\n 60','31-Aug \n 72','\n 84','01-Sep \n 96','\n 108','02-Sep \n 120'])
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
#plt.legend(loc='upper left',fontsize=14)
#plt.ylim([20,165])
plt.xlim([0,126])
plt.xticks(np.arange(0,126,12))

plt.title('Intensity Forecast Error',fontsize=16)
plt.ylabel('Forecast Error (Kt)',fontsize=14)
plt.xlabel('Forecast Lead Time (Hr)',fontsize=14)
plt.text(0,60,'(c)',fontsize=16)

file = folder_fig + 'Int_forec_track_Dorian_' + cycle
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)
