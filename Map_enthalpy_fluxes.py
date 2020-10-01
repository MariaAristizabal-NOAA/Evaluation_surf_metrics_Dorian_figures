#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 09:37:22 2020

@author: aristizabal
"""

#%% User input

lon_lim = [-98.5,-60.0]
lat_lim = [10.0,45.0]

# Time bounds
min_time = '2019/08/20/00'
max_time = '2019/09/07/00'

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
kmz_file_Dorian = '/home/aristizabal/KMZ_files/2019/al052019_best_track.kmz'

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
from datetime import datetime, timedelta
import matplotlib.dates as mdates
import os
import os.path
import glob
#import netCDF4
from bs4 import BeautifulSoup
from zipfile import ZipFile
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
import seawater

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

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

#%% Get Dorian track from models

lon_forec_track_pom_oper, lat_forec_track_pom_oper, lead_time_pom_oper = get_storm_track_POM(hwrf_pom_track_oper)

lon_forec_track_pom_exp, lat_forec_track_pom_exp, lead_time_pom_exp = get_storm_track_POM(hwrf_pom_track_exp)

lon_forec_track_hycom_exp, lat_forec_track_hycom_exp, lead_time_hycom_exp = get_storm_track_POM(hwrf_hycom_track_exp)
   
#%% Time window
date_ini = cycle[0:4]+'/'+cycle[4:6]+'/'+cycle[6:8]+'/'+cycle[8:]+'/00/00'
tini = datetime.strptime(date_ini,'%Y/%m/%d/%H/%M/%S')
tend = tini + timedelta(hours=126)
date_end = tend.strftime('%Y/%m/%d/%H/%M/%S')

#%% Acumulated heat fluxes

#RR_norm_bins = np.arange(0,8,0.1)
RR_norm_bins = np.arange(0,8.1,0.5)
RR_norm_bins_mid = RR_norm_bins[0:-1] + np.diff(RR_norm_bins)/2 

sshla_vec_max_pom_oper = np.empty((3,len(RR_norm_bins)))
sshla_vec_max_pom_oper[:] = np.nan
sshla_vec_mean_pom_oper = np.empty((3,len(RR_norm_bins)))
sshla_vec_mean_pom_oper[:] = np.nan
sshla_vec_min_pom_oper = np.empty((3,len(RR_norm_bins)))
sshla_vec_min_pom_oper[:] = np.nan
area_int_flux_pom_oper = np.empty((3))
area_int_flux_pom_oper[:] = np.nan

time = []

for N,indx in enumerate(np.asarray([6,22,28])):    
    print(HWRF_POM_oper[indx])
    HWRF = xr.open_dataset(HWRF_POM_oper[indx])
    lat_hwrf = np.asarray(HWRF.variables['latitude'][:])
    lon_hwrf = np.asarray(HWRF.variables['longitude'][:])
    time_hwrf = np.asarray(HWRF.variables['time'][:])
    UGRD_hwrf = np.asarray(HWRF.variables['UGRD_10maboveground'][0,:,:])
    VGRD_hwrf = np.asarray(HWRF.variables['VGRD_10maboveground'][0,:,:])
    SHTFL_hwrf = np.asarray(HWRF.variables['SHTFL_surface'][0,:,:])
    LHTFL_hwrf = np.asarray(HWRF.variables['LHTFL_surface'][0,:,:])
    enth_hwrf = SHTFL_hwrf + LHTFL_hwrf
    
    time.append(np.asarray(HWRF.variables['time'][:]))
    
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
    #meshlat = np.meshgrid(lat_hwrf[oklat],lon_hwrf[oklon])
    
    eye_lon = np.tile(lon_forec_track_pom_oper[indx],meshlon_lat[0].shape[1])
    eye_lat = np.tile(lat_forec_track_pom_oper[indx],meshlon_lat[0].shape[0])
    eye_matrix = np.meshgrid(eye_lon,eye_lat)
    
    lat_lon_matrix = np.stack((np.ravel(meshlon_lat[0]),np.ravel(meshlon_lat[1])),axis=1).T
    eye_lat_lon_matrix = np.stack((np.ravel(eye_matrix[0]),np.ravel(eye_matrix[1])),axis=1).T
    
    delta_x = np.diff(lon_hwrf[oklon])[0] 
    delta_y = np.diff(lat_hwrf[oklat])[0]
    area_hwrf = delta_x * delta_y * np.cos(lat_hwrf[oklat]*np.pi/180) * (111111*10**2)**2
    _,area_matrix = np.meshgrid(np.arange(0,len(oklon)),area_hwrf)
        
    R = np.empty(lat_lon_matrix.shape[1])
    R[:] = np.nan
    for i in np.arange(lat_lon_matrix.shape[1]):
        R[i] = seawater.dist([lat_forec_track_pom_oper[indx],lat_lon_matrix[1,i]],\
                             [lon_forec_track_pom_oper[indx],lat_lon_matrix[0,i]],'km')[0][0] 
    
    R_norm = R/Rmax
    okR_4 = R_norm <= 4
    okR_8 = R_norm <= 8
    okR1 = np.logical_and(R_norm <= 1.05,R_norm >= 0.95)
    okR2 = np.logical_and(R_norm <= 2.05,R_norm >= 1.95)
    okR3 = np.logical_and(R_norm <= 3.05,R_norm >= 2.95)
    okR4 = np.logical_and(R_norm <= 4.05,R_norm >= 3.95)
    okR5 = np.logical_and(R_norm <= 5.05,R_norm >= 4.95)
    okR6 = np.logical_and(R_norm <= 6.05,R_norm >= 5.95)
    okR7 = np.logical_and(R_norm <= 7.05,R_norm >= 6.95)
    okR8 = np.logical_and(R_norm <= 8.05,R_norm >= 7.95)
    
    RR_norm = R_norm[okR_8]
    
    #units kj/cm^2
    sea_surf_heat_loss_area_hurr = enth_hwrf[oklat,:][:,oklon] * 3*3600*10**(-7)
    
    sshla_vec = np.ravel(sea_surf_heat_loss_area_hurr)[okR_8]
    sshla_vec[sshla_vec<=0.05] = np.nan
    area_matrix_vec = np.ravel(area_matrix)[okR_8]
    
    area_int_flux_pom_oper[N] = np.nansum(sshla_vec * area_matrix_vec)/np.nansum(area_matrix_vec)
    
    okRR = np.argsort(RR_norm)
    sshla_vec_interp = np.interp(RR_norm_bins,RR_norm[okRR],sshla_vec[okRR])
    
    for n in np.arange(1,len(RR_norm_bins)):
        print(n)
        okR = np.logical_and(RR_norm >= RR_norm_bins[n-1],RR_norm <= RR_norm_bins[n])
        '''
        sshla_vec_max_pom_oper[N,n] = np.nanmax(sshla_vec[okR])
        sshla_vec_mean_pom_oper[N,n] = np.nanmean(sshla_vec[okR])
        sshla_vec_min_pom_oper[N,n] = np.nanmin(sshla_vec[okR])
        '''
        if len(np.where(okR)[0])==0:
            sshla_vec_max_pom_oper[N,n] = np.nan
            sshla_vec_mean_pom_oper[N,n] = np.nan
            sshla_vec_min_pom_oper[N,n] = np.nan
        else:
            maxv = np.nanmean(sshla_vec[okR]) + 3*np.nanstd(sshla_vec[okR])
            minv = np.nanmean(sshla_vec[okR]) - 3*np.nanstd(sshla_vec[okR])
            okk = np.logical_and(sshla_vec[okR] >= minv,sshla_vec[okR] <= maxv)
            sshla_vec_max_pom_oper[N,n] = np.nanmax(sshla_vec[okR][okk])
            sshla_vec_mean_pom_oper[N,n] = np.nanmean(sshla_vec[okR][okk])
            sshla_vec_min_pom_oper[N,n] = np.nanmin(sshla_vec[okR][okk])
                
    plt.figure()
    plt.plot(RR_norm[okRR],sshla_vec[okRR],'.')
    plt.plot(RR_norm_bins_mid,sshla_vec_max_pom_oper[N,1:],'.-')
    plt.plot(RR_norm_bins_mid,sshla_vec_min_pom_oper[N,1:],'.-')
    plt.plot(RR_norm_bins_mid,sshla_vec_mean_pom_oper[N,1:],'.-')
    
    c = sshla_vec
    fig,ax = plt.subplots()
    plt.scatter(lat_lon_matrix[0][okR_8],lat_lon_matrix[1][okR_8],c=c,cmap='YlOrRd',vmin=0,vmax=0.8)
    c=plt.colorbar()
    c.set_label('$KJ/cm^2$',rotation=90, labelpad=15, fontsize=14)
    plt.title('Sea Surface Heat Loss per area on ' + str(time_hwrf[0])[0:13] + '\n HWRF2019-POM (IC Clim.)',fontsize=15)
    plt.plot(eye_lon[0],eye_lat[0],'*k',markersize=10)
    plt.plot(lon_forec_track_pom_oper,lat_forec_track_pom_oper,'.-k') 
    plt.plot(lat_lon_matrix[0][okR1],lat_lon_matrix[1][okR1],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR2],lat_lon_matrix[1][okR2],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR3],lat_lon_matrix[1][okR3],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR4],lat_lon_matrix[1][okR4],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR5],lat_lon_matrix[1][okR5],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR6],lat_lon_matrix[1][okR6],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR7],lat_lon_matrix[1][okR7],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR8],lat_lon_matrix[1][okR8],'.k',markersize=1)
    plt.plot(lon_maxwind,lat_maxwind,'*k',markersize=10)  
    plt.axis('scaled') 
    plt.xlim(np.min(lat_lon_matrix[0][okR_8])-0.1,np.max(lat_lon_matrix[0][okR_8])+0.1) 
    plt.ylim(np.min(lat_lon_matrix[1][okR_8])-0.1,np.max(lat_lon_matrix[1][okR_8])+0.1) 
    
    kw = dict(levels=np.arange(0,0.81,0.1))
    sea_surf_heat_loss_area_hurr[sea_surf_heat_loss_area_hurr>0.79] = 0.79
    fig,ax = plt.subplots()
    plt.contourf(lon_hwrf[oklon],lat_hwrf[oklat],sea_surf_heat_loss_area_hurr,cmap='YlOrRd',\
                 **kw,vmin=0,vmax=0.8)
    c=plt.colorbar()
    #plt.contour(lon_hwrf[oklon],lat_hwrf[oklat],sea_surf_heat_loss_area_hurr,\
    #            cmap='YlOrRd',linestyle='--',**kw)    
    c.set_label('$KJ/cm^2$',rotation=90, labelpad=15, fontsize=14)
    plt.title('Sea Surface Heat Loss per area on ' + str(time_hwrf[0])[0:13] + '\n HWRF2019-POM (IC Clim.)',fontsize=15)
    plt.plot(eye_lon[0],eye_lat[0],'*k',markersize=10)
    plt.plot(lon_forec_track_pom_oper,lat_forec_track_pom_oper,'.-k') 
    plt.plot(lat_lon_matrix[0][okR1],lat_lon_matrix[1][okR1],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR2],lat_lon_matrix[1][okR2],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR3],lat_lon_matrix[1][okR3],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR4],lat_lon_matrix[1][okR4],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR5],lat_lon_matrix[1][okR5],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR6],lat_lon_matrix[1][okR6],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR7],lat_lon_matrix[1][okR7],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR8],lat_lon_matrix[1][okR8],'.k',markersize=1)
    plt.plot(lon_maxwind,lat_maxwind,'*k',markersize=10)  
    plt.axis('scaled') 
    plt.xlim(np.min(lat_lon_matrix[0][okR_8])-0.1,np.max(lat_lon_matrix[0][okR_8])+0.1) 
    plt.ylim(np.min(lat_lon_matrix[1][okR_8])-0.1,np.max(lat_lon_matrix[1][okR_8])+0.1)  

    
#%% Acumulated heat fluxes
        
RR_norm_bins = np.arange(0,8,0.1)
RR_norm_bins_mid = RR_norm_bins[0:-1] + np.diff(RR_norm_bins)/2 

sshla_vec_max_pom_exp = np.empty((3,len(RR_norm_bins)))
sshla_vec_max_pom_exp[:] = np.nan
sshla_vec_mean_pom_exp = np.empty((3,len(RR_norm_bins)))
sshla_vec_mean_pom_exp[:] = np.nan
sshla_vec_min_pom_exp = np.empty((3,len(RR_norm_bins)))
sshla_vec_min_pom_exp[:] = np.nan 
area_int_flux_pom_exp = np.empty((3))
area_int_flux_pom_exp[:] = np.nan       

for N,indx in enumerate(np.asarray([6,22,28])):    
    print(HWRF_POM_exp[indx])
    HWRF = xr.open_dataset(HWRF_POM_exp[indx])
    lat_hwrf = np.asarray(HWRF.variables['latitude'][:])
    lon_hwrf = np.asarray(HWRF.variables['longitude'][:])
    time_hwrf = np.asarray(HWRF.variables['time'][:])
    UGRD_hwrf = np.asarray(HWRF.variables['UGRD_10maboveground'][0,:,:])
    VGRD_hwrf = np.asarray(HWRF.variables['VGRD_10maboveground'][0,:,:])
    SHTFL_hwrf = np.asarray(HWRF.variables['SHTFL_surface'][0,:,:])
    LHTFL_hwrf = np.asarray(HWRF.variables['LHTFL_surface'][0,:,:])
    enth_hwrf = SHTFL_hwrf + LHTFL_hwrf
    
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
    #meshlat = np.meshgrid(lat_hwrf[oklat],lon_hwrf[oklon])
    
    eye_lon = np.tile(lon_forec_track_pom_exp[indx],meshlon_lat[0].shape[1])
    eye_lat = np.tile(lat_forec_track_pom_exp[indx],meshlon_lat[0].shape[0])
    eye_matrix = np.meshgrid(eye_lon,eye_lat)
    
    delta_x = np.diff(lon_hwrf[oklon])[0] 
    delta_y = np.diff(lat_hwrf[oklat])[0]
    area_hwrf = delta_x * delta_y * np.cos(lat_hwrf[oklat]*np.pi/180) * (111111*10**2)**2
    _,area_matrix = np.meshgrid(np.arange(0,len(oklon)),area_hwrf)
    
    lat_lon_matrix = np.stack((np.ravel(meshlon_lat[0]),np.ravel(meshlon_lat[1])),axis=1).T
    eye_lat_lon_matrix = np.stack((np.ravel(eye_matrix[0]),np.ravel(eye_matrix[1])),axis=1).T
    
    R = np.empty(lat_lon_matrix.shape[1])
    R[:] = np.nan
    for i in np.arange(lat_lon_matrix.shape[1]):
        R[i] = seawater.dist([lat_forec_track_pom_exp[indx],lat_lon_matrix[1,i]],\
                             [lon_forec_track_pom_exp[indx],lat_lon_matrix[0,i]],'km')[0][0] 
    
    R_norm = R/Rmax
    okR_4 = R_norm <= 4
    okR_8 = R_norm <= 8
    okR1 = np.logical_and(R_norm <= 1.05,R_norm >= 0.95)
    okR2 = np.logical_and(R_norm <= 2.05,R_norm >= 1.95)
    okR3 = np.logical_and(R_norm <= 3.05,R_norm >= 2.95)
    okR4 = np.logical_and(R_norm <= 4.05,R_norm >= 3.95)
    okR5 = np.logical_and(R_norm <= 5.05,R_norm >= 4.95)
    okR6 = np.logical_and(R_norm <= 6.05,R_norm >= 5.95)
    okR7 = np.logical_and(R_norm <= 7.05,R_norm >= 6.95)
    okR8 = np.logical_and(R_norm <= 8.05,R_norm >= 7.95)
    
    RR_norm = R_norm[okR_8]
    
    #units kj/cm^2
    sea_surf_heat_loss_area_hurr = enth_hwrf[oklat,:][:,oklon] * 3*3600*10**(-7)
    
    sshla_vec = np.ravel(sea_surf_heat_loss_area_hurr)[okR_8]
    sshla_vec[sshla_vec<=0.05] = np.nan
    area_matrix_vec = np.ravel(area_matrix)[okR_8]
    
    area_int_flux_pom_exp[N] = np.nansum(sshla_vec * area_matrix_vec)/np.nansum(area_matrix_vec)
    
    okRR = np.argsort(RR_norm)
    sshla_vec_interp = np.interp(RR_norm_bins,RR_norm[okRR],sshla_vec[okRR])
    
    for n in np.arange(1,len(RR_norm_bins)):
        print(n)
        okR = np.logical_and(RR_norm >= RR_norm_bins[n-1],RR_norm <= RR_norm_bins[n])
        if len(np.where(okR)[0])==0:
            sshla_vec_max_pom_exp[N,n] = np.nan
            sshla_vec_mean_pom_exp[N,n] = np.nan
            sshla_vec_min_pom_exp[N,n] = np.nan
        else:
            
            sshla_vec_max_pom_exp[N,n] = np.nanmax(sshla_vec[okR])
            sshla_vec_mean_pom_exp[N,n] = np.nanmean(sshla_vec[okR])
            sshla_vec_min_pom_exp[N,n] = np.nanmin(sshla_vec[okR])
            '''
            maxv = np.nanmean(sshla_vec[okR]) + 2*np.nanstd(sshla_vec[okR])
            minv = np.nanmean(sshla_vec[okR]) - 2*np.nanstd(sshla_vec[okR])
            okk = np.logical_and(sshla_vec[okR] >= minv,sshla_vec[okR] <= maxv)
            sshla_vec_max_pom_exp[N,n] = np.nanmax(sshla_vec[okR][okk])
            sshla_vec_mean_pom_exp[N,n] = np.nanmean(sshla_vec[okR][okk])
            sshla_vec_min_pom_exp[N,n] = np.nanmin(sshla_vec[okR][okk])
            '''       
    plt.figure()
    plt.plot(RR_norm[okRR],sshla_vec[okRR],'.')
    plt.plot(RR_norm_bins_mid,sshla_vec_max_pom_exp[N,1:],'.-')
    plt.plot(RR_norm_bins_mid,sshla_vec_min_pom_exp[N,1:],'.-')
    plt.plot(RR_norm_bins_mid,sshla_vec_mean_pom_exp[N,1:],'.-')
    plt.ylim(0,1)
    
    sshla_vec = np.ravel(sea_surf_heat_loss_area_hurr)[okR_8]
    c = sshla_vec
    fig,ax = plt.subplots()
    plt.scatter(lat_lon_matrix[0][okR_8],lat_lon_matrix[1][okR_8],c=c,cmap='YlOrRd',vmin=0,vmax=0.8)
    c=plt.colorbar()
    c.set_label('$KJ/cm^2$',rotation=90, labelpad=15, fontsize=14)
    plt.title('Sea Surface Heat Loss per area on ' + str(time_hwrf[0])[0:13] + '\n HWRF2020-POM (IC RTOFS)',fontsize=15)
    plt.plot(eye_lon[0],eye_lat[0],'*k',markersize=10)
    plt.plot(lon_forec_track_pom_exp,lat_forec_track_pom_exp,'.-k') 
    plt.plot(lat_lon_matrix[0][okR1],lat_lon_matrix[1][okR1],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR2],lat_lon_matrix[1][okR2],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR3],lat_lon_matrix[1][okR3],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR4],lat_lon_matrix[1][okR4],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR5],lat_lon_matrix[1][okR5],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR6],lat_lon_matrix[1][okR6],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR7],lat_lon_matrix[1][okR7],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR8],lat_lon_matrix[1][okR8],'.k',markersize=1)
    plt.plot(lon_maxwind,lat_maxwind,'*k',markersize=10)  
    plt.axis('scaled') 
    plt.xlim(np.min(lat_lon_matrix[0][okR_8])-0.1,np.max(lat_lon_matrix[0][okR_8])+0.1) 
    plt.ylim(np.min(lat_lon_matrix[1][okR_8])-0.1,np.max(lat_lon_matrix[1][okR_8])+0.1)  
    #ax.text(np.min(lat_lon_matrix[0][okR8]),np.min(lat_lon_matrix[1][okR8])-0.6,\
    #        'Rmax = ' + str(np.round(Rmax,2)) + ' Km',fontsize=14)
        
#%% Acumulated heat fluxes
        
RR_norm_bins = np.arange(0,8,0.1)
RR_norm_bins_mid = RR_norm_bins[0:-1] + np.diff(RR_norm_bins)/2 

sshla_vec_max_hycom_exp = np.empty((3,len(RR_norm_bins)))
sshla_vec_max_hycom_exp[:] = np.nan
sshla_vec_mean_hycom_exp = np.empty((3,len(RR_norm_bins)))
sshla_vec_mean_hycom_exp[:] = np.nan
sshla_vec_min_hycom_exp = np.empty((3,len(RR_norm_bins)))
sshla_vec_min_hycom_exp[:] = np.nan
area_int_flux_hycom_exp = np.empty((3))
area_int_flux_hycom_exp[:] = np.nan     

for N,indx in enumerate(np.asarray([6,22,28])):    
    print(HWRF_HYCOM_exp[indx])
    HWRF = xr.open_dataset(HWRF_HYCOM_exp[indx])
    lat_hwrf = np.asarray(HWRF.variables['latitude'][:])
    lon_hwrf = np.asarray(HWRF.variables['longitude'][:])
    time_hwrf = np.asarray(HWRF.variables['time'][:])
    UGRD_hwrf = np.asarray(HWRF.variables['UGRD_10maboveground'][0,:,:])
    VGRD_hwrf = np.asarray(HWRF.variables['VGRD_10maboveground'][0,:,:])
    SHTFL_hwrf = np.asarray(HWRF.variables['SHTFL_surface'][0,:,:])
    LHTFL_hwrf = np.asarray(HWRF.variables['LHTFL_surface'][0,:,:])
    enth_hwrf = SHTFL_hwrf + LHTFL_hwrf
    
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
    #meshlat = np.meshgrid(lat_hwrf[oklat],lon_hwrf[oklon])
    
    eye_lon = np.tile(lon_forec_track_hycom_exp[indx],meshlon_lat[0].shape[1])
    eye_lat = np.tile(lat_forec_track_hycom_exp[indx],meshlon_lat[0].shape[0])
    eye_matrix = np.meshgrid(eye_lon,eye_lat)
    
    lat_lon_matrix = np.stack((np.ravel(meshlon_lat[0]),np.ravel(meshlon_lat[1])),axis=1).T
    eye_lat_lon_matrix = np.stack((np.ravel(eye_matrix[0]),np.ravel(eye_matrix[1])),axis=1).T
    
    delta_x = np.diff(lon_hwrf[oklon])[0] 
    delta_y = np.diff(lat_hwrf[oklat])[0]
    area_hwrf = delta_x * delta_y * np.cos(lat_hwrf[oklat]*np.pi/180) * (111111*10**2)**2
    _,area_matrix = np.meshgrid(np.arange(0,len(oklon)),area_hwrf)
    
    R = np.empty(lat_lon_matrix.shape[1])
    R[:] = np.nan
    for i in np.arange(lat_lon_matrix.shape[1]):
        R[i] = seawater.dist([lat_forec_track_hycom_exp[indx],lat_lon_matrix[1,i]],\
                             [lon_forec_track_hycom_exp[indx],lat_lon_matrix[0,i]],'km')[0][0] 
    
    R_norm = R/Rmax
    okR_4 = R_norm <= 4
    okR_8 = R_norm <= 8
    okR1 = np.logical_and(R_norm <= 1.05,R_norm >= 0.95)
    okR2 = np.logical_and(R_norm <= 2.05,R_norm >= 1.95)
    okR3 = np.logical_and(R_norm <= 3.05,R_norm >= 2.95)
    okR4 = np.logical_and(R_norm <= 4.05,R_norm >= 3.95)
    okR5 = np.logical_and(R_norm <= 5.05,R_norm >= 4.95)
    okR6 = np.logical_and(R_norm <= 6.05,R_norm >= 5.95)
    okR7 = np.logical_and(R_norm <= 7.05,R_norm >= 6.95)
    okR8 = np.logical_and(R_norm <= 8.05,R_norm >= 7.95)
    
    RR_norm = R_norm[okR_8]
    
    #units kj/cm^2
    sea_surf_heat_loss_area_hurr = enth_hwrf[oklat,:][:,oklon] * 3*3600*10**(-7)
    
    sshla_vec = np.ravel(sea_surf_heat_loss_area_hurr)[okR_8]
    sshla_vec[sshla_vec<=0.05] = np.nan
    area_matrix_vec = np.ravel(area_matrix)[okR_8]
    
    area_int_flux_hycom_exp[N] = np.nansum(sshla_vec * area_matrix_vec)/np.nansum(area_matrix_vec)
    
    okRR = np.argsort(RR_norm)
    sshla_vec_interp = np.interp(RR_norm_bins,RR_norm[okRR],sshla_vec[okRR])
    
    for n in np.arange(1,len(RR_norm_bins)):
        print(n)
        okR = np.logical_and(RR_norm >= RR_norm_bins[n-1],RR_norm <= RR_norm_bins[n])
        '''
        sshla_vec_max_hycom_exp[N,n] = np.nanmax(sshla_vec[okR])
        sshla_vec_mean_hycom_exp[N,n] = np.nanmean(sshla_vec[okR])
        sshla_vec_min_hycom_exp[N,n] = np.nanmin(sshla_vec[okR])
        '''
        maxv = np.nanmean(sshla_vec[okR]) + 2*np.nanstd(sshla_vec[okR])
        minv = np.nanmean(sshla_vec[okR]) - 2*np.nanstd(sshla_vec[okR])
        okk = np.logical_and(sshla_vec[okR] >= minv,sshla_vec[okR] <= maxv)
        sshla_vec_max_hycom_exp[N,n] = np.nanmax(sshla_vec[okR][okk])
        sshla_vec_mean_hycom_exp[N,n] = np.nanmean(sshla_vec[okR][okk])
        sshla_vec_min_hycom_exp[N,n] = np.nanmin(sshla_vec[okR][okk])
        
    plt.figure()
    plt.plot(RR_norm[okRR],sshla_vec[okRR],'.')
    plt.plot(RR_norm_bins_mid,sshla_vec_max_hycom_exp[N,1:],'.-')
    plt.plot(RR_norm_bins_mid,sshla_vec_min_hycom_exp[N,1:],'.-')
    plt.plot(RR_norm_bins_mid,sshla_vec_mean_hycom_exp[N,1:],'.-')
    
    sshla_vec = np.ravel(sea_surf_heat_loss_area_hurr)[okR_8]
    c = sshla_vec
    fig,ax = plt.subplots()
    plt.scatter(lat_lon_matrix[0][okR_8],lat_lon_matrix[1][okR_8],c=c,cmap='YlOrRd',vmin=0,vmax=0.8)
    c=plt.colorbar()
    c.set_label('$KJ/cm^2$',rotation=90, labelpad=15, fontsize=14)
    plt.title('Sea Surface Heat Loss per area on ' + str(time_hwrf[0])[0:13] + '\n HWRF2020-HYCOM (IC RTOFS)',fontsize=15)
    plt.plot(eye_lon[0],eye_lat[0],'*k',markersize=10)
    plt.plot(lon_forec_track_hycom_exp,lat_forec_track_hycom_exp,'.-k') 
    plt.plot(lat_lon_matrix[0][okR1],lat_lon_matrix[1][okR1],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR2],lat_lon_matrix[1][okR2],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR3],lat_lon_matrix[1][okR3],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR4],lat_lon_matrix[1][okR4],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR5],lat_lon_matrix[1][okR5],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR6],lat_lon_matrix[1][okR6],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR7],lat_lon_matrix[1][okR7],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR8],lat_lon_matrix[1][okR8],'.k',markersize=1)
    plt.plot(lon_maxwind,lat_maxwind,'*k',markersize=10)  
    plt.axis('scaled') 
    plt.xlim(np.min(lat_lon_matrix[0][okR_8])-0.1,np.max(lat_lon_matrix[0][okR_8])+0.1) 
    plt.ylim(np.min(lat_lon_matrix[1][okR_8])-0.1,np.max(lat_lon_matrix[1][okR_8])+0.1)  
    #ax.text(np.min(lat_lon_matrix[0][okR8]),np.min(lat_lon_matrix[1][okR8])-0.6,\
    #        'Rmax = ' + str(np.round(Rmax,2)) + ' Km',fontsize=14)
            
#%%

for N in np.arange(3):
    '''
    aif_pom_oper = str(format(area_int_flux_pom_oper[N],'0.2E'))
    aif_pom_exp = str(format(area_int_flux_pom_exp[N],'0.2E'))
    aif_hycom_exp = str(format(area_int_flux_hycom_exp[N],'0.2E'))
    '''
    aif_pom_oper = str(format(area_int_flux_pom_oper[N],'0.2')) + ' Kj/cm^2'
    aif_pom_exp = str(format(area_int_flux_pom_exp[N],'0.2')) + ' Kj/cm^2'
    aif_hycom_exp = str(format(area_int_flux_hycom_exp[N],'0.2')) + ' Kj/cm^2'
    
    fig,ax = plt.subplots() 
    plt.plot(RR_norm_bins_mid,sshla_vec_mean_pom_oper[N,1:],'X-',color='mediumorchid',\
             label='HWRF2019-POM (IC clim.) '+aif_pom_oper,markeredgecolor='k',markersize=7)
    plt.plot(RR_norm_bins_mid,sshla_vec_max_pom_oper[N,1:],'-',color='mediumorchid',alpha=0.5)
    plt.plot(RR_norm_bins_mid,sshla_vec_min_pom_oper[N,1:],'-',color='mediumorchid',alpha=0.5)
    ax.fill_between(RR_norm_bins_mid,sshla_vec_min_pom_oper[N,1:],sshla_vec_max_pom_oper[N,1:],color='mediumorchid',alpha=0.1)
    
    plt.plot(RR_norm_bins_mid,sshla_vec_mean_pom_exp[N,1:],'^-',color='teal',\
             label='HWRF2020-POM (IC RTOFS) '+aif_pom_exp,markeredgecolor='k',markersize=7)
    plt.plot(RR_norm_bins_mid,sshla_vec_max_pom_exp[N,1:],'-',color='teal',alpha=0.5)
    plt.plot(RR_norm_bins_mid,sshla_vec_min_pom_exp[N,1:],'-',color='teal',alpha=0.5)
    ax.fill_between(RR_norm_bins_mid,sshla_vec_min_pom_exp[N,1:],sshla_vec_max_pom_exp[N,1:],color='teal',alpha=0.1)
    
    plt.plot(RR_norm_bins_mid,sshla_vec_mean_hycom_exp[N,1:],'^-',color='darkorange',\
             label='HWRF2020-HYCOM (IC RTOFS) '+aif_hycom_exp,markeredgecolor='k',markersize=7)
    plt.plot(RR_norm_bins_mid,sshla_vec_max_hycom_exp[N,1:],'-',color='darkorange',alpha=0.5)
    plt.plot(RR_norm_bins_mid,sshla_vec_min_hycom_exp[N,1:],'-',color='darkorange',alpha=0.5)
    ax.fill_between(RR_norm_bins_mid,sshla_vec_min_hycom_exp[N,1:],sshla_vec_max_hycom_exp[N,1:],color='darkorange',alpha=0.1)
    
    plt.title('Sea Surface Heat Loss per area on ' + str(time[N])[2:15],fontsize=16)
    plt.ylim(0,1.8)
    plt.ylabel('$KJ/cm^2$',fontsize=14)
    plt.xlabel('r/Rmax',fontsize=14)
    plt.xticks(np.arange(9))
    plt.legend()
    
#%% Figure forecasted intensity models vs best intensity

okt = np.logical_and(time_best_track >= tini,time_best_track <= tend)

lead_time_pom_oper = np.arange(0,129,3)

fig,ax1 = plt.subplots(figsize=(10, 5))
plt.ion()
plt.plot(lead_time_pom_oper[::2],wind_int_kt[okt],'o-k',label='Best')
plt.plot(lead_time_pom_oper,max_wind_10m_hwrf_pom19_oper,'X-',color='mediumorchid',label='HWRF2010-POM (IC Clim.)',markeredgecolor='k',markersize=7)
plt.plot(lead_time_pom_oper,max_wind_10m_hwrf_pom20_exp,'^-',color='teal',label='HWRF2020-POM (IC RTOFS)',markeredgecolor='k',markersize=7)
plt.plot(lead_time_pom_oper,max_wind_10m_hwrf_hycom20_exp,'H-',color='darkorange',label='HWRF2020-HYCOM (IC RTOFS)',markeredgecolor='k',markersize=7)
plt.plot(np.tile(15,len(np.arange(0,110))),np.arange(0,110),'--',color='k')
#plt.plot(np.tile(24,len(np.arange(0,100))),np.arange(0,100),'--',color='grey')
plt.plot(np.tile(66,len(np.arange(0,110))),np.arange(0,110),'--k')

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
plt.title('Intensity Forecast Dorian '+ cycle,fontsize=18)
plt.ylabel('Max 10m Wind (kt)',fontsize=14)

ax2 = ax1.twinx()
plt.ylim([20,165])
yticks = [64,83,96,113,137]
plt.yticks(yticks,['Cat 1','Cat 2','Cat 3','Cat 4','Cat 5'])
plt.grid(True)

file = folder_fig + 'best_intensity_vs_forec_intensity2_' + cycle 
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 


############################################### 
#%% Acumulated heat fluxes

max_flux1 = []
max_flux2 = []
max_flux3 = []
max_flux4 = []
time_hwrf = []
for N,file in enumerate(HWRF_POM_oper[1:]):
    print(N)
    HWRF = xr.open_dataset(file)
    lat_hwrf = np.asarray(HWRF.variables['latitude'][:])
    lon_hwrf = np.asarray(HWRF.variables['longitude'][:])
    time_hwrf.append(np.asarray(HWRF.variables['time'][:]))
    UGRD_hwrf = np.asarray(HWRF.variables['UGRD_10maboveground'][0,:,:])
    VGRD_hwrf = np.asarray(HWRF.variables['VGRD_10maboveground'][0,:,:])
    SHTFL_hwrf = np.asarray(HWRF.variables['SHTFL_surface'][0,:,:])
    LHTFL_hwrf = np.asarray(HWRF.variables['LHTFL_surface'][0,:,:])
    enth_hwrf = SHTFL_hwrf + LHTFL_hwrf
    
    okflux = np.where(enth_hwrf == np.max(enth_hwrf))
    max_flux1.append(enth_hwrf[okflux[0][0],okflux[1][0]]) # units watts/m^2
    
    wind_int = np.sqrt(UGRD_hwrf**2 + VGRD_hwrf**2)
    max_wind= np.max(wind_int)
    okwind = np.where(wind_int == max_wind)
    max_flux2.append(enth_hwrf[okwind[0][0],okwind[1][0]]) # units watts/m^2
    
    xlim = [lon_forec_track_pom_oper[N]-2,lon_forec_track_pom_oper[N]+2]
    ylim = [lat_forec_track_pom_oper[N]-2,lat_forec_track_pom_oper[N]+2]
    
    oklon = np.where(np.logical_and(lon_hwrf>xlim[0],lon_hwrf<xlim[1]))[0]
    oklat = np.where(np.logical_and(lat_hwrf>ylim[0],lat_hwrf<ylim[1]))[0]
    max_flux3.append(np.mean(enth_hwrf[oklat,:][:,oklon])) # units watts/m^2
    
    max_flux4.append(np.sum(enth_hwrf[oklat,:][:,oklon])*(0.015*111*10**3)**2) # units watts
        
time_hwrf = np.asarray(time_hwrf)
timestamp_hwrf = mdates.date2num(time_hwrf)[:,0]
max_flux1_POM_oper = np.asarray(max_flux1)
max_flux2_POM_oper = np.asarray(max_flux2)
max_flux3_POM_oper = np.asarray(max_flux3)
max_flux4_POM_oper = np.asarray(max_flux4)

int_max_flux1_POM_oper_area1 = np.trapz(max_flux1_POM_oper[0:10],timestamp_hwrf[0:10])*3600*10**-7
int_max_flux2_POM_oper_area1 = np.trapz(max_flux2_POM_oper[0:10],timestamp_hwrf[0:10])*3600*10**-7
int_max_flux3_POM_oper_area1 = np.trapz(max_flux3_POM_oper[0:10],timestamp_hwrf[0:10])*3600*10**-7
int_max_flux4_POM_oper_area1 = np.trapz(max_flux4_POM_oper[0:10],timestamp_hwrf[0:10])*3600*10**-7

int_max_flux1_POM_oper_area2 = np.trapz(max_flux1_POM_oper[0:26],timestamp_hwrf[0:26])*3600*10**-7
int_max_flux2_POM_oper_area2 = np.trapz(max_flux2_POM_oper[0:26],timestamp_hwrf[0:26])*3600*10**-7
int_max_flux3_POM_oper_area2 = np.trapz(max_flux3_POM_oper[0:26],timestamp_hwrf[0:26])*3600*10**-7
int_max_flux4_POM_oper_area2 = np.trapz(max_flux4_POM_oper[0:26],timestamp_hwrf[0:26])*3600*10**-7

#%% Acumulated heat fluxes

max_flux1 = []
max_flux2 = []
max_flux3 = []
max_flux4 = []
time_hwrf = []
for N,file in enumerate(HWRF_POM_exp[1:]):
    print(N)
    HWRF = xr.open_dataset(file)
    lat_hwrf = np.asarray(HWRF.variables['latitude'][:])
    lon_hwrf = np.asarray(HWRF.variables['longitude'][:])
    time_hwrf.append(np.asarray(HWRF.variables['time'][:]))
    UGRD_hwrf = np.asarray(HWRF.variables['UGRD_10maboveground'][0,:,:])
    VGRD_hwrf = np.asarray(HWRF.variables['VGRD_10maboveground'][0,:,:])
    SHTFL_hwrf = np.asarray(HWRF.variables['SHTFL_surface'][0,:,:])
    LHTFL_hwrf = np.asarray(HWRF.variables['LHTFL_surface'][0,:,:])
    enth_hwrf = SHTFL_hwrf + LHTFL_hwrf
    
    okflux = np.where(enth_hwrf == np.max(enth_hwrf))
    max_flux1.append(enth_hwrf[okflux[0][0],okflux[1][0]])
    
    wind_int = np.sqrt(UGRD_hwrf**2 + VGRD_hwrf**2)
    max_wind= np.max(wind_int)
    okwind = np.where(wind_int == max_wind)
    max_flux2.append(enth_hwrf[okwind[0][0],okwind[1][0]])
    
    xlim = [lon_forec_track_pom_oper[N]-2,lon_forec_track_pom_oper[N]+2]
    ylim = [lat_forec_track_pom_oper[N]-2,lat_forec_track_pom_oper[N]+2]
    
    oklon = np.where(np.logical_and(lon_hwrf>xlim[0],lon_hwrf<xlim[1]))[0]
    oklat = np.where(np.logical_and(lat_hwrf>ylim[0],lat_hwrf<ylim[1]))[0]
    max_flux3.append(np.mean(enth_hwrf[oklat,:][:,oklon]))
    
    max_flux4.append(np.sum(enth_hwrf[oklat,:][:,oklon])*(0.015*111*10**3)**2) # units watts
    
time_hwrf = np.asarray(time_hwrf)
timestamp_hwrf = mdates.date2num(time_hwrf)[:,0]
max_flux1_POM_exp = np.asarray(max_flux1)
max_flux2_POM_exp = np.asarray(max_flux2)
max_flux3_POM_exp = np.asarray(max_flux3)
max_flux4_POM_exp = np.asarray(max_flux4)

int_max_flux1_POM_exp_area1 = np.trapz(max_flux1_POM_exp[0:10],timestamp_hwrf[0:10])*3600*10**-7
int_max_flux2_POM_exp_area1 = np.trapz(max_flux2_POM_exp[0:10],timestamp_hwrf[0:10])*3600*10**-7
int_max_flux3_POM_exp_area1 = np.trapz(max_flux3_POM_exp[0:10],timestamp_hwrf[0:10])*3600*10**-7
int_max_flux4_POM_exp_area1 = np.trapz(max_flux4_POM_exp[0:10],timestamp_hwrf[0:10])*360*10**-70

int_max_flux1_POM_exp_area2 = np.trapz(max_flux1_POM_exp[0:26],timestamp_hwrf[0:26])*3600*10**-7
int_max_flux2_POM_exp_area2 = np.trapz(max_flux2_POM_exp[0:26],timestamp_hwrf[0:26])*3600*10**-7
int_max_flux3_POM_exp_area2 = np.trapz(max_flux3_POM_exp[0:26],timestamp_hwrf[0:26])*3600*10**-7
int_max_flux4_POM_exp_area2 = np.trapz(max_flux4_POM_exp[0:26],timestamp_hwrf[0:26])*3600*10**-7

#%% Acumulated heat fluxes

max_flux1 = []
max_flux2 = []
max_flux3 = []
max_flux4 = []
time_hwrf = []
for N,file in enumerate(HWRF_HYCOM_exp[1:]):
    print(N)
    HWRF = xr.open_dataset(file)
    lat_hwrf = np.asarray(HWRF.variables['latitude'][:])
    lon_hwrf = np.asarray(HWRF.variables['longitude'][:])
    time_hwrf.append(np.asarray(HWRF.variables['time'][:]))
    UGRD_hwrf = np.asarray(HWRF.variables['UGRD_10maboveground'][0,:,:])
    VGRD_hwrf = np.asarray(HWRF.variables['VGRD_10maboveground'][0,:,:])
    SHTFL_hwrf = np.asarray(HWRF.variables['SHTFL_surface'][0,:,:])
    LHTFL_hwrf = np.asarray(HWRF.variables['LHTFL_surface'][0,:,:])
    enth_hwrf = SHTFL_hwrf + LHTFL_hwrf
    
    okflux = np.where(enth_hwrf == np.max(enth_hwrf))
    max_flux1.append(enth_hwrf[okflux[0][0],okflux[1][0]])
    
    wind_int = np.sqrt(UGRD_hwrf**2 + VGRD_hwrf**2)
    max_wind= np.max(wind_int)
    okwind = np.where(wind_int == max_wind)
    max_flux2.append(enth_hwrf[okwind[0][0],okwind[1][0]])
    
    xlim = [lon_forec_track_pom_oper[N]-2,lon_forec_track_pom_oper[N]+2]
    ylim = [lat_forec_track_pom_oper[N]-2,lat_forec_track_pom_oper[N]+2]
    
    oklon = np.where(np.logical_and(lon_hwrf>xlim[0],lon_hwrf<xlim[1]))[0]
    oklat = np.where(np.logical_and(lat_hwrf>ylim[0],lat_hwrf<ylim[1]))[0]
    max_flux3.append(np.mean(enth_hwrf[oklat,:][:,oklon]))
    
    max_flux4.append(np.sum(enth_hwrf[oklat,:][:,oklon])*(0.015*111*10**3)**2) # units watts
           
time_hwrf = np.asarray(time_hwrf)
timestamp_hwrf = mdates.date2num(time_hwrf)[:,0]
max_flux1_HYCOM_exp = np.asarray(max_flux1)
max_flux2_HYCOM_exp = np.asarray(max_flux2)
max_flux3_HYCOM_exp = np.asarray(max_flux3)
max_flux4_HYCOM_exp = np.asarray(max_flux4)

int_max_flux1_HYCOM_exp_area1 = np.trapz(max_flux1_HYCOM_exp[0:10],timestamp_hwrf[0:10])*3600*10**-7
int_max_flux2_HYCOM_exp_area1 = np.trapz(max_flux2_HYCOM_exp[0:10],timestamp_hwrf[0:10])*3600*10**-7
int_max_flux3_HYCOM_exp_area1 = np.trapz(max_flux3_HYCOM_exp[0:10],timestamp_hwrf[0:10])*3600*10**-7
int_max_flux4_HYCOM_exp_area1 = np.trapz(max_flux4_HYCOM_exp[0:10],timestamp_hwrf[0:10])*3600*10**-7

int_max_flux1_HYCOM_exp_area2 = np.trapz(max_flux1_HYCOM_exp[0:26],timestamp_hwrf[0:26])*3600*10**-7
int_max_flux2_HYCOM_exp_area2 = np.trapz(max_flux2_HYCOM_exp[0:26],timestamp_hwrf[0:26])*3600*10**-7
int_max_flux3_HYCOM_exp_area2 = np.trapz(max_flux3_HYCOM_exp[0:26],timestamp_hwrf[0:26])*3600*10**-7
int_max_flux4_HYCOM_exp_area2 = np.trapz(max_flux4_HYCOM_exp[0:26],timestamp_hwrf[0:26])*3600*10**-7

#%% Max fluxes
    
fig,ax = plt.subplots(figsize=(10, 5))
plt.plot(time_hwrf,max_flux1_POM_oper,'X-',color='mediumorchid',label='HWRF2010-POM (IC Clim.)',markeredgecolor='k',markersize=7)
plt.plot(time_hwrf,max_flux1_POM_exp,'^-',color='teal',label='HWRF2020-POM (IC RTOFS)',markeredgecolor='k',markersize=7)
plt.plot(time_hwrf,max_flux1_HYCOM_exp,'H-',color='darkorange',label='HWRF2020-HYCOM (IC RTOFS)',markeredgecolor='k',markersize=7)
ax.fill_between(time_hwrf[0:10][:,0],max_flux1_POM_oper[0:10],np.tile(0,10),color='mediumorchid',alpha=0.1)
ax.fill_between(time_hwrf[0:10][:,0],max_flux1_POM_exp[0:10],np.tile(0,10),color='teal',alpha=0.1)
ax.fill_between(time_hwrf[0:10][:,0],max_flux1_HYCOM_exp[0:10],np.tile(0,10),color='darkorange',alpha=0.1)
plt.legend(loc='upper left')
ax.text(time_hwrf[3],50,format(int_max_flux1_POM_oper_area1,'0.2e')+' $j/m^2$',color='mediumorchid',size=14)
ax.text(time_hwrf[3],200,format(int_max_flux1_POM_exp_area1,'0.2e')+' $j/m^2$',color='teal',size=14)
ax.text(time_hwrf[3],350,format(int_max_flux1_HYCOM_exp_area1,'0.2e')+' $j/m^2$',color='darkorange',size=14)
plt.title('Maximum Enthalpy Flux',size=16)

#%% Max fluxes
    
fig,ax = plt.subplots(figsize=(10, 5))
plt.plot(time_hwrf,max_flux1_POM_oper,'X-',color='mediumorchid',label='HWRF2010-POM (IC Clim.)',markeredgecolor='k',markersize=7)
plt.plot(time_hwrf,max_flux1_POM_exp,'^-',color='teal',label='HWRF2020-POM (IC RTOFS)',markeredgecolor='k',markersize=7)
plt.plot(time_hwrf,max_flux1_HYCOM_exp,'H-',color='darkorange',label='HWRF2020-HYCOM (IC RTOFS)',markeredgecolor='k',markersize=7)
ax.fill_between(time_hwrf[0:26][:,0],max_flux1_POM_oper[0:26],np.tile(0,26),color='mediumorchid',alpha=0.1)
ax.fill_between(time_hwrf[0:26][:,0],max_flux1_POM_exp[0:26],np.tile(0,26),color='teal',alpha=0.1)
ax.fill_between(time_hwrf[0:26][:,0],max_flux1_HYCOM_exp[0:26],np.tile(0,26),color='darkorange',alpha=0.1)
plt.legend(loc='upper left')
ax.text(time_hwrf[13],50,format(int_max_flux1_POM_oper_area2,'0.2e')+' $j/m^2$',color='mediumorchid',size=14)
ax.text(time_hwrf[13],200,format(int_max_flux1_POM_exp_area2,'0.2e')+' $j/m^2$',color='teal',size=14)
ax.text(time_hwrf[13],350,format(int_max_flux1_HYCOM_exp_area2,'0.2e')+' $j/m^2$',color='darkorange',size=14)
plt.title('Maximum Enthalpy Flux',size=16)

#%% fluxes at max wind
    
fig,ax = plt.subplots(figsize=(10, 5))
plt.plot(time_hwrf,max_flux2_POM_oper,'X-',color='mediumorchid',label='HWRF2010-POM (IC Clim.)',markeredgecolor='k',markersize=7)
plt.plot(time_hwrf,max_flux2_POM_exp,'^-',color='teal',label='HWRF2020-POM (IC RTOFS)',markeredgecolor='k',markersize=7)
plt.plot(time_hwrf,max_flux2_HYCOM_exp,'H-',color='darkorange',label='HWRF2020-HYCOM (IC RTOFS)',markeredgecolor='k',markersize=7)
ax.fill_between(time_hwrf[0:10][:,0],max_flux2_POM_oper[0:10],np.tile(0,10),color='mediumorchid',alpha=0.1)
ax.fill_between(time_hwrf[0:10][:,0],max_flux2_POM_exp[0:10],np.tile(0,10),color='teal',alpha=0.1)
ax.fill_between(time_hwrf[0:10][:,0],max_flux2_HYCOM_exp[0:10],np.tile(0,10),color='darkorange',alpha=0.1)
plt.legend(loc='upper left')
ax.text(time_hwrf[3],50,format(int_max_flux2_POM_oper_area1,'0.2e')+' $j/m^2$',color='mediumorchid',size=14)
ax.text(time_hwrf[3],200,format(int_max_flux2_POM_exp_area1,'0.2e')+' $j/m^2$',color='teal',size=14)
ax.text(time_hwrf[3],350,format(int_max_flux2_HYCOM_exp_area1,'0.2e')+' $j/m^2$',color='darkorange',size=14)
plt.title('Enthalpy Flux at location of Maximum Wind Intensity',size=16)

#%% fluxes at max wind
    
fig,ax = plt.subplots(figsize=(10, 5))
plt.plot(time_hwrf,max_flux2_POM_oper,'X-',color='mediumorchid',label='HWRF2010-POM (IC Clim.)',markeredgecolor='k',markersize=7)
plt.plot(time_hwrf,max_flux2_POM_exp,'^-',color='teal',label='HWRF2020-POM (IC RTOFS)',markeredgecolor='k',markersize=7)
plt.plot(time_hwrf,max_flux2_HYCOM_exp,'H-',color='darkorange',label='HWRF2020-HYCOM (IC RTOFS)',markeredgecolor='k',markersize=7)
ax.fill_between(time_hwrf[0:26][:,0],max_flux2_POM_oper[0:26],np.tile(0,26),color='mediumorchid',alpha=0.1)
ax.fill_between(time_hwrf[0:26][:,0],max_flux2_POM_exp[0:26],np.tile(0,26),color='teal',alpha=0.1)
ax.fill_between(time_hwrf[0:26][:,0],max_flux2_HYCOM_exp[0:26],np.tile(0,26),color='darkorange',alpha=0.1)
plt.legend(loc='upper left')
ax.text(time_hwrf[13],50,format(int_max_flux2_POM_oper_area2,'0.2e')+' $j/m^2$',color='mediumorchid',size=14)
ax.text(time_hwrf[13],200,format(int_max_flux2_POM_exp_area2,'0.2e')+' $j/m^2$',color='teal',size=14)
ax.text(time_hwrf[13],350,format(int_max_flux2_HYCOM_exp_area2,'0.2e')+' $j/m^2$',color='darkorange',size=14)
plt.title('Enthalpy Flux at location of Maximum Wind Intensity',size=16)

#%% Mean fluxes 
    
fig,ax = plt.subplots(figsize=(10, 5))
plt.plot(time_hwrf,max_flux3_POM_oper,'X-',color='mediumorchid',label='HWRF2010-POM (IC Clim.)',markeredgecolor='k',markersize=7)
plt.plot(time_hwrf,max_flux3_POM_exp,'^-',color='teal',label='HWRF2020-POM (IC RTOFS)',markeredgecolor='k',markersize=7)
plt.plot(time_hwrf,max_flux3_HYCOM_exp,'H-',color='darkorange',label='HWRF2020-HYCOM (IC RTOFS)',markeredgecolor='k',markersize=7)
ax.fill_between(time_hwrf[0:10][:,0],max_flux3_POM_oper[0:10],np.tile(0,10),color='mediumorchid',alpha=0.1)
ax.fill_between(time_hwrf[0:10][:,0],max_flux3_POM_exp[0:10],np.tile(0,10),color='teal',alpha=0.1)
ax.fill_between(time_hwrf[0:10][:,0],max_flux3_HYCOM_exp[0:10],np.tile(0,10),color='darkorange',alpha=0.1)
plt.legend(loc='lower right')
plt.ylabel('$W/m^2$',fontsize=14)
ax.text(time_hwrf[3],50,format(int_max_flux3_POM_oper_area1,'0.2e')+' $KJ/cm^2$',color='mediumorchid',size=14)
ax.text(time_hwrf[3],100,format(int_max_flux3_POM_exp_area1,'0.2e')+' $KJ/cm^2$',color='teal',size=14)
ax.text(time_hwrf[3],150,format(int_max_flux3_HYCOM_exp_area1,'0.2e')+' $KJ/cm^2$',color='darkorange',size=14)
plt.title('Mean Enthalpy Flux +- 2 Degrees around the Eye',size=16)

#%% Mean fluxes
    
fig,ax = plt.subplots(figsize=(10, 5))
plt.plot(time_hwrf,max_flux3_POM_oper,'X-',color='mediumorchid',label='HWRF2010-POM (IC Clim.)',markeredgecolor='k',markersize=7)
plt.plot(time_hwrf,max_flux3_POM_exp,'^-',color='teal',label='HWRF2020-POM (IC RTOFS)',markeredgecolor='k',markersize=7)
plt.plot(time_hwrf,max_flux3_HYCOM_exp,'H-',color='darkorange',label='HWRF2020-HYCOM (IC RTOFS)',markeredgecolor='k',markersize=7)
ax.fill_between(time_hwrf[0:26][:,0],max_flux3_POM_oper[0:26],np.tile(0,26),color='mediumorchid',alpha=0.1)
ax.fill_between(time_hwrf[0:26][:,0],max_flux3_POM_exp[0:26],np.tile(0,26),color='teal',alpha=0.1)
ax.fill_between(time_hwrf[0:26][:,0],max_flux3_HYCOM_exp[0:26],np.tile(0,26),color='darkorange',alpha=0.1)
plt.legend(loc='lower right')
plt.ylabel('$W/m^2$',fontsize=14)
ax.text(time_hwrf[11],50,format(int_max_flux3_POM_oper_area2,'0.2e')+' $KJ/cm^2$',color='mediumorchid',size=14)
ax.text(time_hwrf[11],100,format(int_max_flux3_POM_exp_area2,'0.2e')+' $KJ/cm^2$',color='teal',size=14)
ax.text(time_hwrf[11],150,format(int_max_flux3_HYCOM_exp_area2,'0.2e')+' $KJ/cm^2$',color='darkorange',size=14)
plt.title('Mean Enthalpy Flux +- 2 Degrees around the Eye',size=16)

#%% total power: j/s
    
fig,ax = plt.subplots(figsize=(10, 5))
plt.plot(time_hwrf,max_flux4_POM_oper,'X-',color='mediumorchid',label='HWRF2010-POM (IC Clim.)',markeredgecolor='k',markersize=7)
plt.plot(time_hwrf,max_flux4_POM_exp,'^-',color='teal',label='HWRF2020-POM (IC RTOFS)',markeredgecolor='k',markersize=7)
plt.plot(time_hwrf,max_flux4_HYCOM_exp,'H-',color='darkorange',label='HWRF2020-HYCOM (IC RTOFS)',markeredgecolor='k',markersize=7)
ax.fill_between(time_hwrf[0:10][:,0],max_flux4_POM_oper[0:10],np.tile(0,10),color='mediumorchid',alpha=0.1)
ax.fill_between(time_hwrf[0:10][:,0],max_flux4_POM_exp[0:10],np.tile(0,10),color='teal',alpha=0.1)
ax.fill_between(time_hwrf[0:10][:,0],max_flux4_HYCOM_exp[0:10],np.tile(0,10),color='darkorange',alpha=0.1)
plt.legend(loc='lower right')
ax.text(time_hwrf[3],0,format(int_max_flux4_POM_oper_area1,'0.2e')+' $J$',color='mediumorchid',size=14)
ax.text(time_hwrf[3],10**13,format(int_max_flux4_POM_exp_area1,'0.2e')+' $J$',color='teal',size=14)
ax.text(time_hwrf[3],2*10**13,format(int_max_flux4_HYCOM_exp_area1,'0.2e')+' $J$',color='darkorange',size=14)
plt.title('Total Power (j/s) +- 2 Degrees around the Eye',size=16)

#%% total power: j/s
    
fig,ax = plt.subplots(figsize=(10, 5))
plt.plot(time_hwrf,max_flux4_POM_oper,'X-',color='mediumorchid',label='HWRF2010-POM (IC Clim.)',markeredgecolor='k',markersize=7)
plt.plot(time_hwrf,max_flux4_POM_exp,'^-',color='teal',label='HWRF2020-POM (IC RTOFS)',markeredgecolor='k',markersize=7)
plt.plot(time_hwrf,max_flux4_HYCOM_exp,'H-',color='darkorange',label='HWRF2020-HYCOM (IC RTOFS)',markeredgecolor='k',markersize=7)
ax.fill_between(time_hwrf[0:26][:,0],max_flux4_POM_oper[0:26],np.tile(0,26),color='mediumorchid',alpha=0.1)
ax.fill_between(time_hwrf[0:26][:,0],max_flux4_POM_exp[0:26],np.tile(0,26),color='teal',alpha=0.1)
ax.fill_between(time_hwrf[0:26][:,0],max_flux4_HYCOM_exp[0:26],np.tile(0,26),color='darkorange',alpha=0.1)
plt.legend(loc='lower right')
ax.text(time_hwrf[11],0,format(int_max_flux4_POM_oper_area2,'0.2e')+' $J$',color='mediumorchid',size=14)
ax.text(time_hwrf[11],10**13,format(int_max_flux4_POM_exp_area2,'0.2e')+' $J$',color='teal',size=14)
ax.text(time_hwrf[11],2*10**13,format(int_max_flux4_HYCOM_exp_area2,'0.2e')+' $J$',color='darkorange',size=14)
plt.title('Total Power (j/s) +- 2 Degrees around the Eye',size=16)



#%% map wind vectors HWRP2019-POM operational

N = 5
kw = dict(levels=np.linspace(0,1200,13))
ylim = [17,18]
xlim = [-65.5,-64.5]

file = HWRF_POM_oper[N]
HWRF = xr.open_dataset(file)
lat_hwrf = np.asarray(HWRF.variables['latitude'][:])
lon_hwrf = np.asarray(HWRF.variables['longitude'][:])
t_hwrf = np.asarray(HWRF.variables['time'][:])
UGRD_hwrf = np.asarray(HWRF.variables['UGRD_10maboveground'][0,:,:])
VGRD_hwrf = np.asarray(HWRF.variables['VGRD_10maboveground'][0,:,:])
SHTFL_hwrf = np.asarray(HWRF.variables['SHTFL_surface'][0,:,:])
LHTFL_hwrf = np.asarray(HWRF.variables['LHTFL_surface'][0,:,:])
enth_hwrf = SHTFL_hwrf + LHTFL_hwrf

fig,ax = plt.subplots(figsize=(5,5)) 
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')   
plt.contourf(lon_hwrf,lat_hwrf,enth_hwrf,cmap='RdYlBu_r',**kw)
c = plt.colorbar(shrink=0.8)
#q = plt.quiver(lon_hwrf[::30], lat_hwrf[::30],UGRD_hwrf[::30,::30],VGRD_hwrf[::30,::30]) #,units='xy' ,scale=0.01)
#plt.quiverkey(q,-63.0,28.5,30,"30 m/s",coordinates='data',color='k',fontproperties={'size': 14})
q = plt.quiver(lon_hwrf[::5], lat_hwrf[::5],UGRD_hwrf[::5,::5],VGRD_hwrf[::5,::5],scale=400)#,units='xy' ,scale=0.01)
plt.quiverkey(q,-64.36,18.05,30,"30 m/s",coordinates='data',color='k',fontproperties={'size': 14})
plt.title('HWRF2019-POM Enthalpy Flux \n on '+str(t_hwrf[0])[0:13],fontsize=14)
c.set_label('$W/m^2$',rotation=90, labelpad=1, fontsize=14)
c.ax.tick_params(labelsize=14)
plt.axis('scaled')
plt.ylim([ylim[0],ylim[1]])
plt.xlim([xlim[0],xlim[1]])

file_name = folder_fig + 'Dorian_heat_fluxes_HWRF19_POM_' + str(t_hwrf[0])[0:13] + '_' + cycle
plt.savefig(file_name)#,bbox_inches = 'tight',pad_inches = 0.1) 

#%% map wind vectors HWRP2020-POM experimental

file = HWRF_POM_exp[N]
HWRF = xr.open_dataset(file)
lat_hwrf = np.asarray(HWRF.variables['latitude'][:])
lon_hwrf = np.asarray(HWRF.variables['longitude'][:])
t_hwrf = np.asarray(HWRF.variables['time'][:])
UGRD_hwrf = np.asarray(HWRF.variables['UGRD_10maboveground'][0,:,:])
VGRD_hwrf = np.asarray(HWRF.variables['VGRD_10maboveground'][0,:,:])
SHTFL_hwrf = np.asarray(HWRF.variables['SHTFL_surface'][0,:,:])
LHTFL_hwrf = np.asarray(HWRF.variables['LHTFL_surface'][0,:,:])
enth_hwrf = SHTFL_hwrf + LHTFL_hwrf

fig,ax = plt.subplots(figsize=(5,5)) 
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')   
plt.contourf(lon_hwrf,lat_hwrf,enth_hwrf,cmap='RdYlBu_r',**kw)
c = plt.colorbar(shrink=0.8)
#q = plt.quiver(lon_hwrf[::30], lat_hwrf[::30],UGRD_hwrf[::30,::30],VGRD_hwrf[::30,::30]) #,units='xy' ,scale=0.01)
q = plt.quiver(lon_hwrf[::5], lat_hwrf[::5],UGRD_hwrf[::5,::5],VGRD_hwrf[::5,::5],scale=400)#,units='xy' ,scale=0.01)
plt.quiverkey(q,-64.36,18.05,30,"30 m/s",coordinates='data',color='k',fontproperties={'size': 14})
plt.title('HWRF2020-POM Enthalpy Flux \n on '+str(t_hwrf[0])[0:13],fontsize=14)
c.set_label('$W/m^2$',rotation=90, labelpad=1, fontsize=16)
c.ax.tick_params(labelsize=14)
plt.axis('scaled')
plt.ylim([ylim[0],ylim[1]])
plt.xlim([xlim[0],xlim[1]])

file_name = folder_fig + 'Dorian_heat_fluxes_HWRF20_POM_' + str(t_hwrf[0])[0:13] + '_' + cycle
plt.savefig(file_name)#,bbox_inches = 'tight',pad_inches = 0.1) 

#%% map wind vectors HWRP2020-HYCOM experimental

file = HWRF_HYCOM_exp[N]
HWRF = xr.open_dataset(file)
lat_hwrf = np.asarray(HWRF.variables['latitude'][:])
lon_hwrf = np.asarray(HWRF.variables['longitude'][:])
t_hwrf = np.asarray(HWRF.variables['time'][:])
UGRD_hwrf = np.asarray(HWRF.variables['UGRD_10maboveground'][0,:,:])
VGRD_hwrf = np.asarray(HWRF.variables['VGRD_10maboveground'][0,:,:])
SHTFL_hwrf = np.asarray(HWRF.variables['SHTFL_surface'][0,:,:])
LHTFL_hwrf = np.asarray(HWRF.variables['LHTFL_surface'][0,:,:])
enth_hwrf = SHTFL_hwrf + LHTFL_hwrf

fig,ax = plt.subplots(figsize=(5,5)) 
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')   
plt.contourf(lon_hwrf,lat_hwrf,enth_hwrf,cmap='RdYlBu_r',**kw)
c = plt.colorbar(shrink=0.8)
q = plt.quiver(lon_hwrf[::5], lat_hwrf[::5],UGRD_hwrf[::5,::5],VGRD_hwrf[::5,::5],scale=400)#,units='xy' ,scale=0.01)
plt.quiverkey(q,-64.36,18.05,30,"30 m/s",coordinates='data',color='k',fontproperties={'size': 14})
plt.title('HWRF2020-HYCOM Enthalpy Flux \n on '+str(t_hwrf[0])[0:13],fontsize=14)
c.set_label('$W/m^2$',rotation=90, labelpad=1, fontsize=16)
c.ax.tick_params(labelsize=14)
plt.axis('scaled')
plt.ylim([ylim[0],ylim[1]])
plt.xlim([xlim[0],xlim[1]])

file_name = folder_fig + 'Dorian_heat_fluxes_HWRF20_HYCOM_' + str(t_hwrf[0])[0:13] + '_' + cycle
plt.savefig(file_name)#,bbox_inches = 'tight',pad_inches = 0.1) 

#%%

N = 22
ylim = [24.6,26.1]
xlim = [-72.5,-71]

file = HWRF_POM_oper[N]
HWRF = xr.open_dataset(file)
lat_hwrf = np.asarray(HWRF.variables['latitude'][:])
lon_hwrf = np.asarray(HWRF.variables['longitude'][:])
t_hwrf = np.asarray(HWRF.variables['time'][:])
UGRD_hwrf = np.asarray(HWRF.variables['UGRD_10maboveground'][0,:,:])
VGRD_hwrf = np.asarray(HWRF.variables['VGRD_10maboveground'][0,:,:])
SHTFL_hwrf = np.asarray(HWRF.variables['SHTFL_surface'][0,:,:])
LHTFL_hwrf = np.asarray(HWRF.variables['LHTFL_surface'][0,:,:])
enth_hwrf = SHTFL_hwrf + LHTFL_hwrf

fig,ax = plt.subplots(figsize=(5,5)) 
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')   
plt.contourf(lon_hwrf,lat_hwrf,enth_hwrf,cmap='RdYlBu_r',**kw)
c = plt.colorbar(shrink=0.8)
#q = plt.quiver(lon_hwrf[::30], lat_hwrf[::30],UGRD_hwrf[::30,::30],VGRD_hwrf[::30,::30]) #,units='xy' ,scale=0.01)
q = plt.quiver(lon_hwrf[::5], lat_hwrf[::5],UGRD_hwrf[::5,::5],VGRD_hwrf[::5,::5],scale=400)#,units='xy' ,scale=0.01)
plt.quiverkey(q,-70.8,26.2,30,"30 m/s",coordinates='data',color='k',fontproperties={'size': 14})
plt.title('HWRF2019-POM Enthalpy Flux \n on '+str(t_hwrf[0])[0:13],fontsize=14)
c.set_label('$W/m^2$',rotation=90, labelpad=1, fontsize=14)
c.ax.tick_params(labelsize=14)
plt.axis('scaled')
plt.ylim([ylim[0],ylim[1]])
plt.xlim([xlim[0],xlim[1]])
plt.yticks([24.8,25.2,25.6,26])

file_name = folder_fig + 'Dorian_heat_fluxes_HWRF19_POM_' + str(t_hwrf[0])[0:13] + '_' + cycle
plt.savefig(file_name)#,bbox_inches = 'tight',pad_inches = 0.1) 

#%% map wind vectors HWRP2020-POM experimental

file = HWRF_POM_exp[N]
HWRF = xr.open_dataset(file)
lat_hwrf = np.asarray(HWRF.variables['latitude'][:])
lon_hwrf = np.asarray(HWRF.variables['longitude'][:])
t_hwrf = np.asarray(HWRF.variables['time'][:])
UGRD_hwrf = np.asarray(HWRF.variables['UGRD_10maboveground'][0,:,:])
VGRD_hwrf = np.asarray(HWRF.variables['VGRD_10maboveground'][0,:,:])
SHTFL_hwrf = np.asarray(HWRF.variables['SHTFL_surface'][0,:,:])
LHTFL_hwrf = np.asarray(HWRF.variables['LHTFL_surface'][0,:,:])
enth_hwrf = SHTFL_hwrf + LHTFL_hwrf

fig,ax = plt.subplots(figsize=(5,5)) 
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')   
plt.contourf(lon_hwrf,lat_hwrf,enth_hwrf,cmap='RdYlBu_r',**kw)
c = plt.colorbar(shrink=0.8)
#q = plt.quiver(lon_hwrf[::30], lat_hwrf[::30],UGRD_hwrf[::30,::30],VGRD_hwrf[::30,::30]) #,units='xy' ,scale=0.01)
q = plt.quiver(lon_hwrf[::5], lat_hwrf[::5],UGRD_hwrf[::5,::5],VGRD_hwrf[::5,::5],scale=400)#,units='xy' ,scale=0.01)
plt.quiverkey(q,-70.8,26.2,30,"30 m/s",coordinates='data',color='k',fontproperties={'size': 14})
plt.title('HWRF2020-POM Enthalpy Flux \n on '+str(t_hwrf[0])[0:13],fontsize=14)
c.set_label('$W/m^2$',rotation=90, labelpad=1, fontsize=16)
c.ax.tick_params(labelsize=14)
plt.axis('scaled')
plt.ylim([ylim[0],ylim[1]])
plt.xlim([xlim[0],xlim[1]])
plt.yticks([24.8,25.2,25.6,26])

file_name = folder_fig + 'Dorian_heat_fluxes_HWRF20_POM_' + str(t_hwrf[0])[0:13] + '_' + cycle
plt.savefig(file_name)#,bbox_inches = 'tight',pad_inches = 0.1) 

#%% map wind vectors HWRP2020-HYCOM experimental

file = HWRF_HYCOM_exp[N]
HWRF = xr.open_dataset(file)
lat_hwrf = np.asarray(HWRF.variables['latitude'][:])
lon_hwrf = np.asarray(HWRF.variables['longitude'][:])
t_hwrf = np.asarray(HWRF.variables['time'][:])
UGRD_hwrf = np.asarray(HWRF.variables['UGRD_10maboveground'][0,:,:])
VGRD_hwrf = np.asarray(HWRF.variables['VGRD_10maboveground'][0,:,:])
SHTFL_hwrf = np.asarray(HWRF.variables['SHTFL_surface'][0,:,:])
LHTFL_hwrf = np.asarray(HWRF.variables['LHTFL_surface'][0,:,:])
enth_hwrf = SHTFL_hwrf + LHTFL_hwrf

fig,ax = plt.subplots(figsize=(5,5)) 
plt.contour(bath_lon,bath_lat,bath_elev,[0],colors='k')   
plt.contourf(lon_hwrf,lat_hwrf,enth_hwrf,cmap='RdYlBu_r',**kw)
c = plt.colorbar(shrink=0.8)
#q = plt.quiver(lon_hwrf[::30], lat_hwrf[::30],UGRD_hwrf[::30,::30],VGRD_hwrf[::30,::30]) #,units='xy' ,scale=0.01)
q = plt.quiver(lon_hwrf[::5], lat_hwrf[::5],UGRD_hwrf[::5,::5],VGRD_hwrf[::5,::5],scale=400)#,units='xy' ,scale=0.01)
plt.quiverkey(q,-70.8,26.2,30,"30 m/s",coordinates='data',color='k',fontproperties={'size': 14})
plt.title('HWRF2020-HYCOM Enthalpy Flux \n on '+str(t_hwrf[0])[0:13],fontsize=14)
c.set_label('$W/m^2$',rotation=90, labelpad=1, fontsize=16)
c.ax.tick_params(labelsize=14)
plt.axis('scaled')
plt.ylim([ylim[0],ylim[1]])
plt.xlim([xlim[0],xlim[1]])
plt.yticks([24.8,25.2,25.6,26])

file_name = folder_fig + 'Dorian_heat_fluxes_HWRF20_HYCOM_' + str(t_hwrf[0])[0:13] + '_' + cycle
plt.savefig(file_name)#,bbox_inches = 'tight',pad_inches = 0.1) 
