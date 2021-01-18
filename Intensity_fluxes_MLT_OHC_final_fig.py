#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 10:56:28 2021

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

Dir_HMON_HYCOM = '/home/aristizabal/HWRF_RTOFS_Michael_2018/HWRF-Hycom_exp_Michael/'
# RTOFS grid file name
hycom_grid_exp = Dir_HMON_HYCOM + 'hwrf_rtofs_hat10.basin.regional.grid'

# Dorian track files
hwrf_hycom_track_exp = folder_hycom_exp + 'dorian05l.' + cycle + '.trak.hwrf.atcfunix'

# folder nc files hwrf
folder_hwrf_hycom20_exp = folder_hycom20 + 'HWRF2020_HYCOM_dorian05l.' + cycle + '_grb2_to_nc_exp/'

#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import xarray as xr
from datetime import datetime, timedelta
import os
import os.path
from bs4 import BeautifulSoup
from zipfile import ZipFile
import glob
import cmocean
import seawater
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)

import sys
sys.path.append('/home/aristizabal/NCEP_scripts')
from utils4HYCOM import readgrids
from utils4HYCOM import readBinz

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)
plt.rc('legend',fontsize=12)

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

#%%  Calculation of mixed layer depth based the temperature and density criteria, Tmean: mean temp within the
# mixed layer

def MLD_temp_and_dens_criteria_from_3D_fields(dt,drho,depth,temp,salt,dens):

    MLD_temp_crit = np.empty((temp.shape[1],temp.shape[2]))
    MLD_temp_crit[:] = np.nan
    Tmean_temp_crit = np.empty((temp.shape[1],temp.shape[2]))
    Tmean_temp_crit[:] = np.nan
    Smean_temp_crit = np.empty((temp.shape[1],temp.shape[2]))
    Smean_temp_crit[:] = np.nan
    MLD_dens_crit = np.empty((temp.shape[1],temp.shape[2]))
    MLD_dens_crit[:] = np.nan
    Tmean_dens_crit = np.empty((temp.shape[1],temp.shape[2]))
    Tmean_dens_crit[:] = np.nan
    Smean_dens_crit = np.empty((temp.shape[1],temp.shape[2]))
    Smean_dens_crit[:] = np.nan
    for y in np.arange(temp.shape[1]):
        print(y)
        for x in np.arange(temp.shape[2]):
            if depth.ndim == 1:
                d10 = np.where(depth >= 10)[0][0]
            if depth.ndim == 3:
                d10 = np.where(depth[:,y,x] >= -10)[0][-1]
            T10 = temp[d10,y,x]
            delta_T = T10 - temp[:,y,x]
            ok_mld_temp = np.where(delta_T <= dt)[0]
            rho10 = dens[d10,y,x]
            delta_rho = -(rho10 - dens[:,y,x])
            ok_mld_rho = np.where(delta_rho <= drho)[0]

            if ok_mld_temp.size == 0:
                MLD_temp_crit[y,x] = np.nan
                Tmean_temp_crit[y,x] = np.nan
                Smean_temp_crit[y,x] = np.nan
            else:
                if depth.ndim == 1:
                    MLD_temp_crit[y,x] = depth[ok_mld_temp[-1]]
                if depth.ndim == 3:
                    MLD_temp_crit[y,x] = depth[ok_mld_temp[-1],y,x]
                Tmean_temp_crit[y,x] = np.nanmean(temp[ok_mld_temp,y,x])
                Smean_temp_crit[y,x] = np.nanmean(salt[ok_mld_temp,y,x])

            if ok_mld_rho.size == 0:
                MLD_dens_crit[y,x] = np.nan
                Tmean_dens_crit[y,x] = np.nan
                Smean_dens_crit[y,x] = np.nan
            else:
                if depth.ndim == 1:
                    MLD_dens_crit[y,x] = depth[ok_mld_rho[-1]]
                if depth.ndim == 3:
                    MLD_dens_crit[y,x] = depth[ok_mld_rho[-1],y,x]
                Tmean_dens_crit[y,x] = np.nanmean(temp[ok_mld_rho,y,x])
                Smean_dens_crit[y,x] = np.nanmean(salt[ok_mld_rho,y,x])

    return MLD_temp_crit,Tmean_temp_crit,Smean_temp_crit,\
               MLD_dens_crit,Tmean_dens_crit,Smean_dens_crit

#%% Function Ocean Heat Content

def OHC_from_3D_fields(temp,depth,dens):
    cp = 3985 #Heat capacity in J/(kg K)

    OHC = np.empty((temp.shape[1],temp.shape[2]))
    OHC[:] = np.nan
    for y in np.arange(temp.shape[1]):
        print(y)
        for x in np.arange(temp.shape[2]):
            ok26 = temp[:,y,x] >= 26
            if depth.ndim == 1:
                if len(depth[ok26]) != 0:
                    if np.nanmin(np.abs(depth[ok26]))>10:
                        OHC[y,x] = np.nan
                    else:
                        rho0 = np.nanmean(dens[ok26,y,x])
                        OHC[y,x] = np.abs(cp * rho0 * np.trapz(temp[ok26,y,x]-26,depth[ok26]))
                else:
                    OHC[y,x] = np.nan

            if depth.ndim == 3:
                if len(depth[ok26,y,x]) != 0:
                    if np.nanmin(np.abs(depth[ok26,y,x]))>10:
                        OHC[y,x] = np.nan
                    else:
                        rho0 = np.nanmean(dens[ok26,y,x])
                        OHC[y,x] = np.abs(cp * rho0 * np.trapz(temp[ok26,y,x]-26,depth[ok26,y,x]))
                else:
                    OHC[y,x] = np.nan

    return OHC

#%% Reading bathymetry data
'''
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
'''

#%% Get list HWRF files

HWRF_POM_oper = sorted(glob.glob(os.path.join(folder_hwrf_pom19_oper,'*.nc')))
HWRF_POM_exp = sorted(glob.glob(os.path.join(folder_hwrf_pom20_exp,'*.nc')))
HWRF_HYCOM_exp = sorted(glob.glob(os.path.join(folder_hwrf_hycom20_exp,'*.nc')))

#%% Get Dorian track from models

lon_forec_track_pom_oper, lat_forec_track_pom_oper, lead_time_pom_oper = get_storm_track_POM(hwrf_pom_track_oper)

lon_forec_track_pom_exp, lat_forec_track_pom_exp, lead_time_pom_exp = get_storm_track_POM(hwrf_pom_track_exp)

lon_forec_track_hycom_exp, lat_forec_track_hycom_exp, lead_time_hycom_exp = get_storm_track_POM(hwrf_hycom_track_exp)

#%% Get Dorian best track

lon_best_track, lat_best_track, time_best_track, _, wind_int_kt, _ = \
read_kmz_file_storm_best_track(kmz_file_Dorian)

#%% Obtain lat and lon band around forecated track operational
'''
dlon = 0.1
nlevels = int(2*delta_lon /dlon) + 1

lon_bnd = np.linspace(lon_forec_track_pom_oper[2*Nini:2*Nend-1]-delta_lon,lon_forec_track_pom_oper[2*Nini:2*Nend-1]+delta_lon,nlevels)
lon_band = lon_bnd.ravel()
lat_bd = np.tile(lat_forec_track_pom_oper[2*Nini:2*Nend-1],lon_bnd.shape[0])
lat_bnd = lat_bd.reshape(lon_bnd.shape[0],lon_bnd.shape[1])
lat_band = lat_bnd.ravel()
'''
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

#%%Time window
date_ini = cycle[0:4]+'-'+cycle[4:6]+'-'+cycle[6:8]+' '+cycle[8:]+':00:00'
tini = datetime.strptime(date_ini,'%Y-%m-%d %H:%M:%S')
tend = tini + timedelta(hours=float(lead_time_pom_oper[-1]))
date_end = str(tend)

okt = np.logical_and(time_best_track >= tini,time_best_track <= tend)

#%% Area integrated sea surface heat loss POM oper

area_int_sshl_pom_oper = np.empty((len(HWRF_POM_oper)))
area_int_sshl_pom_oper[:] = np.nan

sshla_R4_mean_pom_oper = np.empty((len(HWRF_POM_oper)))
sshla_R4_mean_pom_oper[:] = np.nan
sshla_R4_min_pom_oper = np.empty((len(HWRF_POM_oper)))
sshla_R4_min_pom_oper[:] = np.nan
sshla_R4_max_pom_oper = np.empty((len(HWRF_POM_oper)))
sshla_R4_max_pom_oper[:] = np.nan

area_int_sshl_pom_oper_R4 = np.empty((len(HWRF_POM_oper)))
area_int_sshl_pom_oper_R4[:] = np.nan

for N,indx in enumerate(np.arange(len(HWRF_POM_oper))):
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

    wind_int = np.sqrt(UGRD_hwrf**2 + VGRD_hwrf**2)
    max_wind= np.max(wind_int)
    okwind = np.where(wind_int == max_wind)
    lat_maxwind = lat_hwrf[okwind[0][0]]
    lon_maxwind = lon_hwrf[okwind[1][0]]
    Rmax = seawater.dist([lat_forec_track_pom_oper[indx],lat_maxwind],\
                         [lon_forec_track_pom_oper[indx],lon_maxwind],'km')[0][0]

    xlim = [lon_forec_track_pom_oper[N]-2,lon_forec_track_pom_oper[N]+2]
    ylim = [lat_forec_track_pom_oper[N]-2,lat_forec_track_pom_oper[N]+2]

    oklon = np.where(np.logical_and(lon_hwrf>xlim[0],lon_hwrf<xlim[1]))[0]
    oklat = np.where(np.logical_and(lat_hwrf>ylim[0],lat_hwrf<ylim[1]))[0]

    meshlon_lat = np.meshgrid(lon_hwrf[oklon],lat_hwrf[oklat])

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

    #units kj/cm^2
    sea_surf_heat_loss_area_hurr = enth_hwrf[oklat,:][:,oklon] * 3*3600*10**(-7)

    sshla_R4 = np.ravel(sea_surf_heat_loss_area_hurr)[okR_4]
    area_matrix_vec_R4 = np.ravel(area_matrix)[okR_4]
    area_int_sshl_pom_oper_R4[N] = np.nansum(sshla_R4 * area_matrix_vec_R4)

    sshla_R4_mean_pom_oper[N] = np.nanmean(sshla_R4)
    sshla_R4_min_pom_oper[N] = np.nanmin(sshla_R4)
    sshla_R4_max_pom_oper[N] = np.nanmax(sshla_R4)

    sshla_vec = np.ravel(sea_surf_heat_loss_area_hurr)
    sshla_vec[sshla_vec<=0.05] = np.nan
    area_matrix_vec = np.ravel(area_matrix)

    area_int_sshl_pom_oper[N] = np.nansum(sshla_vec * area_matrix_vec)
    #area_int_sshl_pom_oper2[N] = np.nansum(sshla_vec * area_matrix_vec)/np.nansum(area_matrix_vec)

#%% Area integrated sea surface heat loss POM exp

area_int_sshl_pom_exp = np.empty((len(HWRF_POM_exp)))
area_int_sshl_pom_exp[:] = np.nan

sshla_R4_mean_pom_exp = np.empty((len(HWRF_POM_exp)))
sshla_R4_mean_pom_exp[:] = np.nan
sshla_R4_min_pom_exp = np.empty((len(HWRF_POM_exp)))
sshla_R4_min_pom_exp[:] = np.nan
sshla_R4_max_pom_exp = np.empty((len(HWRF_POM_exp)))
sshla_R4_max_pom_exp[:] = np.nan

area_int_sshl_pom_exp_R4 = np.empty((len(HWRF_POM_exp)))
area_int_sshl_pom_exp_R4[:] = np.nan

for N,indx in enumerate(np.arange(len(HWRF_POM_exp))):
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

    xlim = [lon_forec_track_pom_oper[N]-2,lon_forec_track_pom_oper[N]+2]
    ylim = [lat_forec_track_pom_oper[N]-2,lat_forec_track_pom_oper[N]+2]

    oklon = np.where(np.logical_and(lon_hwrf>xlim[0],lon_hwrf<xlim[1]))[0]
    oklat = np.where(np.logical_and(lat_hwrf>ylim[0],lat_hwrf<ylim[1]))[0]

    meshlon_lat = np.meshgrid(lon_hwrf[oklon],lat_hwrf[oklat])

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

    #units kj/cm^2
    sea_surf_heat_loss_area_hurr = enth_hwrf[oklat,:][:,oklon] * 3*3600*10**(-7)

    sshla_R4 = np.ravel(sea_surf_heat_loss_area_hurr)[okR_4]
    area_matrix_vec_R4 = np.ravel(area_matrix)[okR_4]
    area_int_sshl_pom_exp_R4[N] = np.nansum(sshla_R4 * area_matrix_vec_R4)

    sshla_R4 = np.ravel(sea_surf_heat_loss_area_hurr)[okR_4]
    sshla_R4_mean_pom_exp[N] = np.nanmean(sshla_R4)
    sshla_R4_min_pom_exp[N] = np.nanmin(sshla_R4)
    sshla_R4_max_pom_exp[N] = np.nanmax(sshla_R4)

    sshla_vec = np.ravel(sea_surf_heat_loss_area_hurr)
    sshla_vec[sshla_vec<=0.05] = np.nan
    area_matrix_vec = np.ravel(area_matrix)

    area_int_sshl_pom_exp[N] = np.nansum(sshla_vec * area_matrix_vec)

#%% Area integrated sea surface heat loss HYCOM exp

area_int_sshl_hycom_exp = np.empty((len(HWRF_HYCOM_exp)))
area_int_sshl_hycom_exp[:] = np.nan

sshla_R4_mean_hycom_exp = np.empty((len(HWRF_HYCOM_exp)))
sshla_R4_mean_hycom_exp[:] = np.nan
sshla_R4_min_hycom_exp = np.empty((len(HWRF_HYCOM_exp)))
sshla_R4_min_hycom_exp[:] = np.nan
sshla_R4_max_hycom_exp = np.empty((len(HWRF_HYCOM_exp)))
sshla_R4_max_hycom_exp[:] = np.nan

area_int_sshl_hycom_exp_R4 = np.empty((len(HWRF_HYCOM_exp)))
area_int_sshl_hycom_exp_R4[:] = np.nan

for N,indx in enumerate(np.arange(len(HWRF_HYCOM_exp))):
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

    xlim = [lon_forec_track_pom_oper[N]-2,lon_forec_track_pom_oper[N]+2]
    ylim = [lat_forec_track_pom_oper[N]-2,lat_forec_track_pom_oper[N]+2]

    oklon = np.where(np.logical_and(lon_hwrf>xlim[0],lon_hwrf<xlim[1]))[0]
    oklat = np.where(np.logical_and(lat_hwrf>ylim[0],lat_hwrf<ylim[1]))[0]

    meshlon_lat = np.meshgrid(lon_hwrf[oklon],lat_hwrf[oklat])

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

    #units kj/cm^2
    sea_surf_heat_loss_area_hurr = enth_hwrf[oklat,:][:,oklon] * 3*3600*10**(-7)

    sshla_R4 = np.ravel(sea_surf_heat_loss_area_hurr)[okR_4]
    area_matrix_vec_R4 = np.ravel(area_matrix)[okR_4]
    area_int_sshl_hycom_exp_R4[N] = np.nansum(sshla_R4 * area_matrix_vec_R4)

    sshla_R4 = np.ravel(sea_surf_heat_loss_area_hurr)[okR_4]
    sshla_R4_mean_hycom_exp[N] = np.nanmean(sshla_R4)
    sshla_R4_min_hycom_exp[N] = np.nanmin(sshla_R4)
    sshla_R4_max_hycom_exp[N] = np.nanmax(sshla_R4)

    sshla_vec = np.ravel(sea_surf_heat_loss_area_hurr)
    sshla_vec[sshla_vec<=0.05] = np.nan
    area_matrix_vec = np.ravel(area_matrix)

    area_int_sshl_hycom_exp[N] = np.nansum(sshla_vec * area_matrix_vec)

#%% Calculate acumulated area integrated sea surface heat loss 84 hours lead time

acum_ai_sshl_pom_oper = np.sum(area_int_sshl_pom_oper[0:29])
acum_ai_sshl_pom_exp = np.sum(area_int_sshl_pom_exp[0:29])
acum_ai_sshl_hycom_exp = np.sum(area_int_sshl_hycom_exp[0:29])

#%%%%%%%%%%%%%%%%%%

#%% Get list HWRF files

HWRF_POM_oper = sorted(glob.glob(os.path.join(folder_hwrf_pom19_oper,'*.nc')))[::2]
HWRF_POM_exp = sorted(glob.glob(os.path.join(folder_hwrf_pom20_exp,'*.nc')))[::2]
HWRF_HYCOM_exp = sorted(glob.glob(os.path.join(folder_hwrf_hycom20_exp,'*.nc')))[::2]

#%% Get list ocean model files

POM_oper = sorted(glob.glob(os.path.join(folder_pom_oper,prefix_pom+'*.nc')))
POM_exp = sorted(glob.glob(os.path.join(folder_pom_exp,prefix_pom+'*.nc')))
HYCOM_exp = sorted(glob.glob(os.path.join(folder_hycom_exp,prefix_hycom+'*.a')))

#%% Read POM grid

print('Retrieving coordinates from POM')
POM_grid_oper = xr.open_dataset(pom_grid_oper,decode_times=False)
lon_pom_oper = np.asarray(POM_grid_oper['east_e'][:])
lat_pom_oper = np.asarray(POM_grid_oper['north_e'][:])
zlev_pom_oper = np.asarray(POM_grid_oper['zz'][:])
hpom_oper = np.asarray(POM_grid_oper['h'][:])
zmatrix = np.dot(hpom_oper.reshape(-1,1),zlev_pom_oper.reshape(1,-1)).T
zmatrix_pom_oper = zmatrix.reshape(zlev_pom_oper.shape[0],hpom_oper.shape[0],hpom_oper.shape[1])

POM_grid_exp = xr.open_dataset(pom_grid_exp,decode_times=False)
lon_pom_exp = np.asarray(POM_grid_exp['east_e'][:])
lat_pom_exp = np.asarray(POM_grid_exp['north_e'][:])
zlev_pom_exp = np.asarray(POM_grid_exp['zz'][:])
hpom_exp = np.asarray(POM_grid_exp['h'][:])
zmatrix = np.dot(hpom_exp.reshape(-1,1),zlev_pom_exp.reshape(1,-1)).T
zmatrix_pom_exp = zmatrix.reshape(zlev_pom_exp.shape[0],hpom_exp.shape[0],hpom_exp.shape[1])

#%% Reading HYCOM grid

# Reading lat and lon
lines_grid = [line.rstrip() for line in open(hycom_grid_exp+'.b')]
lon_hycom = np.array(readgrids(hycom_grid_exp,'plon:',[0]))
lat_hycom = np.array(readgrids(hycom_grid_exp,'plat:',[0]))

# Extracting the longitudinal and latitudinal size array
idm=int([line.split() for line in lines_grid if 'longitudinal' in line][0][0])
jdm=int([line.split() for line in lines_grid if 'latitudinal' in line][0][0])

afiles = sorted(glob.glob(os.path.join(folder_hycom_exp,prefix_hycom+'*.a')))

# Reading depths
lines=[line.rstrip() for line in open(afiles[0][:-2]+'.b')]
z = []
for line in lines[6:]:
    if line.split()[2]=='temp':
        #print(line.split()[1])
        z.append(float(line.split()[1]))
depth_HYCOM_exp = np.asarray(z)

nz = len(depth_HYCOM_exp)

#%% MLT and OHC under Dorian POM operational

RR_norm_bins = np.arange(0,8.1,0.5)
RR_norm_bins_mid = RR_norm_bins[0:-1] + np.diff(RR_norm_bins)/2

Tmean_max_pom_oper = np.empty(len(HWRF_POM_oper))
Tmean_max_pom_oper[:] = np.nan
Tmean_min_pom_oper = np.empty(len(HWRF_POM_oper))
Tmean_min_pom_oper[:] = np.nan
Tmean_mean_pom_oper = np.empty(len(HWRF_POM_oper))
Tmean_mean_pom_oper[:] = np.nan

OHC_max_pom_oper = np.empty(len(HWRF_POM_oper))
OHC_max_pom_oper[:] = np.nan
OHC_min_pom_oper = np.empty(len(HWRF_POM_oper))
OHC_min_pom_oper[:] = np.nan
OHC_mean_pom_oper = np.empty(len(HWRF_POM_oper))
OHC_mean_pom_oper[:] = np.nan

time_HWRF = []
time_POM = []

for n,file in enumerate(HWRF_POM_oper[0:-1]):
    print(file)
    HWRF = xr.open_dataset(file)
    lat_hwrf = np.asarray(HWRF.variables['latitude'][:])
    lon_hwrf = np.asarray(HWRF.variables['longitude'][:])
    time_hwrf = np.asarray(HWRF.variables['time'][:])
    UGRD_hwrf = np.asarray(HWRF.variables['UGRD_10maboveground'][0,:,:])
    VGRD_hwrf = np.asarray(HWRF.variables['VGRD_10maboveground'][0,:,:])

    time_HWRF.append(np.asarray(HWRF.variables['time'][:]))

    wind_int = np.sqrt(UGRD_hwrf**2 + VGRD_hwrf**2)
    max_wind= np.max(wind_int)
    okwind = np.where(wind_int == max_wind)
    lat_maxwind = lat_hwrf[okwind[0][0]]
    lon_maxwind = lon_hwrf[okwind[1][0]]
    Rmax = seawater.dist([lat_forec_track_pom_oper[::2][n],lat_maxwind],\
                         [lon_forec_track_pom_oper[::2][n],lon_maxwind],'km')[0][0]

    xlim = [lon_forec_track_pom_oper[::2][n]-4,lon_forec_track_pom_oper[::2][n]+4]
    ylim = [lat_forec_track_pom_oper[::2][n]-4,lat_forec_track_pom_oper[::2][n]+4]

    oklon = np.where(np.logical_and(lon_hwrf>xlim[0],lon_hwrf<xlim[1]))[0]
    oklat = np.where(np.logical_and(lat_hwrf>ylim[0],lat_hwrf<ylim[1]))[0]

    # POM oper
    print(POM_oper[n])
    pom = xr.open_dataset(POM_oper[n])
    tpom = np.asarray(pom['time'][:])
    temp_pom = np.asarray(pom['t'][0,:,:,:])
    salt_pom = np.asarray(pom['s'][0,:,:,:])
    dens_pom = np.asarray(pom['rho'][0,:,:,:]) * 1000 + 1000
    dens_pom[dens_pom==1000.0] = np.nan

    time_POM.append(tpom)

    oklon_pom = np.where(np.logical_and(lon_pom_oper[0,:]>xlim[0],lon_pom_oper[0,:]<xlim[1]))[0]
    oklat_pom = np.where(np.logical_and(lat_pom_oper[:,0]>ylim[0],lat_pom_oper[:,0]<ylim[1]))[0]

    meshlon_lat_pom = np.meshgrid(lon_pom_oper[0,oklon_pom],lat_pom_oper[oklat_pom,0])

    lat_lon_matrix_pom = np.stack((np.ravel(meshlon_lat_pom[0]),np.ravel(meshlon_lat_pom[1])),axis=1).T

    R_pom = np.empty(lat_lon_matrix_pom.shape[1])
    R_pom[:] = np.nan
    for i in np.arange(lat_lon_matrix_pom.shape[1]):
        R_pom[i] = seawater.dist([lat_forec_track_pom_oper[::2][n],lat_lon_matrix_pom[1,i]],\
                             [lon_forec_track_pom_oper[::2][n],lat_lon_matrix_pom[0,i]],'km')[0][0]

    R_pom_norm = R_pom/Rmax
    okR_pom_8 = R_pom_norm <= 8
    okR_pom_4 = R_pom_norm <= 4


    dt = 0.2
    drho = 0.125

    timestamp_POM_oper = mdates.date2num(tpom)
    target_temp_POM_oper = temp_pom[:,oklat_pom,:][:,:,oklon_pom]
    target_salt_POM_oper = salt_pom[:,oklat_pom,:][:,:,oklon_pom]
    target_dens_POM_oper = dens_pom[:,oklat_pom,:][:,:,oklon_pom]
    target_depth_POM_oper = zmatrix_pom_oper[:,oklat_pom,:][:,:,oklon_pom]

    # for POM operational
    MLD_temp_crit,Tmean_temp_crit,Smean_temp_crit,\
    MLD_dens_crit,Tmean_dens_crit,Smean_dens_crit = \
    MLD_temp_and_dens_criteria_from_3D_fields(dt,drho,target_depth_POM_oper,target_temp_POM_oper,target_salt_POM_oper,target_dens_POM_oper)

    Tmean_dens_crit_vec = np.ravel(Tmean_dens_crit)[okR_pom_4]

    # for POM operational
    OHC = OHC_from_3D_fields(target_temp_POM_oper,target_depth_POM_oper,\
                             target_dens_POM_oper)*10**(-7)

    OHC_vec = np.ravel(OHC)[okR_pom_4]

    Tmean_max_pom_oper[n] = np.nanmax(Tmean_dens_crit_vec)
    Tmean_mean_pom_oper[n] = np.nanmean(Tmean_dens_crit_vec)
    Tmean_min_pom_oper[n] = np.nanmin(Tmean_dens_crit_vec)

    OHC_max_pom_oper[n] = np.nanmax(OHC_vec)
    OHC_mean_pom_oper[n] = np.nanmean(OHC_vec)
    OHC_min_pom_oper[n] = np.nanmin(OHC_vec)

#%% MLT and OHC under Dorian POM experimental

Tmean_max_pom_exp = np.empty(len(HWRF_POM_exp))
Tmean_max_pom_exp[:] = np.nan
Tmean_min_pom_exp = np.empty(len(HWRF_POM_exp))
Tmean_min_pom_exp[:] = np.nan
Tmean_mean_pom_exp = np.empty(len(HWRF_POM_exp))
Tmean_mean_pom_exp[:] = np.nan

OHC_max_pom_exp = np.empty(len(HWRF_POM_exp))
OHC_max_pom_exp[:] = np.nan
OHC_min_pom_exp = np.empty(len(HWRF_POM_exp))
OHC_min_pom_exp[:] = np.nan
OHC_mean_pom_exp = np.empty(len(HWRF_POM_exp))
OHC_mean_pom_exp[:] = np.nan

time_HWRF_exp = []
time_POM_exp = []

for n,file in enumerate(HWRF_POM_exp[0:-1]):
    print(HWRF_POM_exp[n])
    HWRF = xr.open_dataset(HWRF_POM_exp[n])
    lat_hwrf = np.asarray(HWRF.variables['latitude'][:])
    lon_hwrf = np.asarray(HWRF.variables['longitude'][:])
    time_hwrf = np.asarray(HWRF.variables['time'][:])
    UGRD_hwrf = np.asarray(HWRF.variables['UGRD_10maboveground'][0,:,:])
    VGRD_hwrf = np.asarray(HWRF.variables['VGRD_10maboveground'][0,:,:])

    time_HWRF_exp.append(np.asarray(HWRF.variables['time'][:]))

    wind_int = np.sqrt(UGRD_hwrf**2 + VGRD_hwrf**2)
    max_wind= np.max(wind_int)
    okwind = np.where(wind_int == max_wind)
    lat_maxwind = lat_hwrf[okwind[0][0]]
    lon_maxwind = lon_hwrf[okwind[1][0]]
    Rmax = seawater.dist([lat_forec_track_pom_exp[::2][n],lat_maxwind],\
                         [lon_forec_track_pom_exp[::2][n],lon_maxwind],'km')[0][0]

    xlim = [lon_forec_track_pom_exp[::2][n]-4,lon_forec_track_pom_exp[::2][n]+4]
    ylim = [lat_forec_track_pom_exp[::2][n]-4,lat_forec_track_pom_exp[::2][n]+4]

    oklon = np.where(np.logical_and(lon_hwrf>xlim[0],lon_hwrf<xlim[1]))[0]
    oklat = np.where(np.logical_and(lat_hwrf>ylim[0],lat_hwrf<ylim[1]))[0]

    # POM exp
    print(POM_exp[n])
    pom = xr.open_dataset(POM_exp[n])
    tpom = np.asarray(pom['time'][:])
    temp_pom = np.asarray(pom['t'][0,:,:,:])
    salt_pom = np.asarray(pom['s'][0,:,:,:])
    dens_pom = np.asarray(pom['rho'][0,:,:,:]) * 1000 + 1000
    dens_pom[dens_pom==1000.0] = np.nan

    time_POM_exp.append(tpom)

    oklon_pom = np.where(np.logical_and(lon_pom_exp[0,:]>xlim[0],lon_pom_exp[0,:]<xlim[1]))[0]
    oklat_pom = np.where(np.logical_and(lat_pom_exp[:,0]>ylim[0],lat_pom_exp[:,0]<ylim[1]))[0]

    meshlon_lat_pom = np.meshgrid(lon_pom_exp[0,oklon_pom],lat_pom_exp[oklat_pom,0])

    lat_lon_matrix_pom = np.stack((np.ravel(meshlon_lat_pom[0]),np.ravel(meshlon_lat_pom[1])),axis=1).T

    R_pom = np.empty(lat_lon_matrix_pom.shape[1])
    R_pom[:] = np.nan
    for i in np.arange(lat_lon_matrix_pom.shape[1]):
        R_pom[i] = seawater.dist([lat_forec_track_pom_exp[::2][n],lat_lon_matrix_pom[1,i]],\
                             [lon_forec_track_pom_exp[::2][n],lat_lon_matrix_pom[0,i]],'km')[0][0]

    R_pom_norm = R_pom/Rmax
    okR_pom_8 = R_pom_norm <= 8
    okR_pom_4 = R_pom_norm <= 4

    dt = 0.2
    drho = 0.125

    timestamp_POM_exp = mdates.date2num(tpom)
    target_temp_POM_exp = temp_pom[:,oklat_pom,:][:,:,oklon_pom]
    target_salt_POM_exp = salt_pom[:,oklat_pom,:][:,:,oklon_pom]
    target_dens_POM_exp = dens_pom[:,oklat_pom,:][:,:,oklon_pom]
    target_depth_POM_exp = zmatrix_pom_exp[:,oklat_pom,:][:,:,oklon_pom]

    # for POM operational
    _,_,_,_,Tmean_dens_crit,_ = \
    MLD_temp_and_dens_criteria_from_3D_fields(dt,drho,target_depth_POM_exp,target_temp_POM_exp,\
                                              target_salt_POM_exp,target_dens_POM_exp)

    Tmean_dens_crit_vec = np.ravel(Tmean_dens_crit)[okR_pom_4]

    # for POM operational
    OHC = OHC_from_3D_fields(target_temp_POM_exp,target_depth_POM_exp,\
                             target_dens_POM_exp)*10**(-7)

    OHC_vec = np.ravel(OHC)[okR_pom_4]

    Tmean_max_pom_exp[n] = np.nanmax(Tmean_dens_crit_vec)
    Tmean_mean_pom_exp[n] = np.nanmean(Tmean_dens_crit_vec)
    Tmean_min_pom_exp[n] = np.nanmin(Tmean_dens_crit_vec)

    OHC_max_pom_exp[n] = np.nanmax(OHC_vec)
    OHC_mean_pom_exp[n] = np.nanmean(OHC_vec)
    OHC_min_pom_exp[n] = np.nanmin(OHC_vec)

#%% MLT and OHC under Dorian HYCOM experimental

Tmean_max_hycom_exp = np.empty(len(HWRF_HYCOM_exp))
Tmean_max_hycom_exp[:] = np.nan
Tmean_min_hycom_exp = np.empty(len(HWRF_HYCOM_exp))
Tmean_min_hycom_exp[:] = np.nan
Tmean_mean_hycom_exp = np.empty(len(HWRF_HYCOM_exp))
Tmean_mean_hycom_exp[:] = np.nan

OHC_max_hycom_exp = np.empty(len(HWRF_HYCOM_exp))
OHC_max_hycom_exp[:] = np.nan
OHC_min_hycom_exp = np.empty(len(HWRF_HYCOM_exp))
OHC_min_hycom_exp[:] = np.nan
OHC_mean_hycom_exp = np.empty(len(HWRF_HYCOM_exp))
OHC_mean_hycom_exp[:] = np.nan

time_hwrf_hycom_exp = []
time_hycom_exp = []

afiles = sorted(glob.glob(os.path.join(folder_hycom_exp,prefix_hycom+'*.a')))

for n,file in enumerate(HWRF_HYCOM_exp[0:-1]):
    print(HWRF_HYCOM_exp[n])
    HWRF = xr.open_dataset(HWRF_HYCOM_exp[n])
    lat_hwrf = np.asarray(HWRF.variables['latitude'][:])
    lon_hwrf = np.asarray(HWRF.variables['longitude'][:])
    time_hwrf = np.asarray(HWRF.variables['time'][:])
    UGRD_hwrf = np.asarray(HWRF.variables['UGRD_10maboveground'][0,:,:])
    VGRD_hwrf = np.asarray(HWRF.variables['VGRD_10maboveground'][0,:,:])
    SHTFL_hwrf = np.asarray(HWRF.variables['SHTFL_surface'][0,:,:])
    LHTFL_hwrf = np.asarray(HWRF.variables['LHTFL_surface'][0,:,:])
    enth_hwrf = SHTFL_hwrf + LHTFL_hwrf

    time_hwrf_hycom_exp.append(np.asarray(HWRF.variables['time'][:]))

    wind_int = np.sqrt(UGRD_hwrf**2 + VGRD_hwrf**2)
    max_wind= np.max(wind_int)
    okwind = np.where(wind_int == max_wind)
    lat_maxwind = lat_hwrf[okwind[0][0]]
    lon_maxwind = lon_hwrf[okwind[1][0]]
    Rmax = seawater.dist([lat_forec_track_hycom_exp[::2][n],lat_maxwind],\
                         [lon_forec_track_hycom_exp[::2][n],lon_maxwind],'km')[0][0]

    xlim = [lon_forec_track_hycom_exp[::2][n]-4,lon_forec_track_hycom_exp[::2][n]+4]
    ylim = [lat_forec_track_hycom_exp[::2][n]-4,lat_forec_track_hycom_exp[::2][n]+4]

    oklon = np.where(np.logical_and(lon_hwrf>xlim[0],lon_hwrf<xlim[1]))[0]
    oklat = np.where(np.logical_and(lat_hwrf>ylim[0],lat_hwrf<ylim[1]))[0]

    # Read HYCOM exp
    file = afiles[n]

    #Reading time stamp
    year = int(file.split('/')[-1].split('.')[1][0:4])
    month = int(file.split('/')[-1].split('.')[1][4:6])
    day = int(file.split('/')[-1].split('.')[1][6:8])
    hour = int(file.split('/')[-1].split('.')[1][8:10])
    dt = int(file.split('/')[-1].split('.')[3][1:])
    timestamp_hycom = mdates.date2num(datetime(year,month,day,hour)) + dt/24
    time_hycom_exp.append(datetime(year,month,day,hour))

    # Reading 3D variable from binary file
    temp_hycom = np.asarray(readBinz(file[:-2],'3z','temp'))
    temp_hycom[temp_hycom>100] = np.nan
    salt_hycom = np.asarray(readBinz(file[:-2],'3z','salinity'))
    salt_hycom[salt_hycom>100] = np.nan

    oklon_hycom = np.where(np.logical_and(lon_hycom[0,:]-360>xlim[0],lon_hycom[0,:]-360<xlim[1]))[0]
    oklat_hycom = np.where(np.logical_and(lat_hycom[:,0]>ylim[0],lat_hycom[:,0]<ylim[1]))[0]

    meshlon_lat_hycom = np.meshgrid(lon_hycom[0,oklon_hycom],lat_hycom[oklat_hycom,0])

    lat_lon_matrix_hycom = np.stack((np.ravel(meshlon_lat_hycom[0]),np.ravel(meshlon_lat_hycom[1])),axis=1).T

    R_hycom = np.empty(lat_lon_matrix_hycom.shape[1])
    R_hycom[:] = np.nan
    for i in np.arange(lat_lon_matrix_hycom.shape[1]):
        R_hycom[i] = seawater.dist([lat_forec_track_hycom_exp[::2][n],lat_lon_matrix_hycom[1,i]],\
                             [lon_forec_track_hycom_exp[::2][n],lat_lon_matrix_hycom[0,i]],'km')[0][0]

    R_hycom_norm = R_hycom/Rmax
    okR_hycom_8 = R_hycom_norm <= 8
    okR_hycom_4 = R_hycom_norm <= 4

    dt = 0.2
    drho = 0.125

    target_temp_hycom_exp = np.asarray(temp_hycom[oklat_hycom,:,:][:,oklon_hycom,:])
    target_salt_hycom_exp = np.asarray(salt_hycom[oklat_hycom,:,:][:,oklon_hycom,:])
    target_depth_hycom_exp = np.asarray(depth_HYCOM_exp)

    target_dens_hycom_exp = seawater.dens(target_salt_hycom_exp,target_temp_hycom_exp,target_depth_hycom_exp)

    # for POM operational
    _,_,_,_,Tmean_dens_crit,_ = \
    MLD_temp_and_dens_criteria_from_3D_fields(dt,drho,target_depth_hycom_exp,target_temp_hycom_exp.T,\
                                              target_salt_hycom_exp.T,target_dens_hycom_exp.T)

    Tmean_dens_crit = Tmean_dens_crit.T

    Tmean_dens_crit_vec = np.ravel(Tmean_dens_crit)[okR_hycom_4]

    # for HYCOM exp
    OHC = OHC_from_3D_fields(target_temp_hycom_exp.T,target_depth_hycom_exp,\
                             target_dens_hycom_exp.T)*10**(-7)
    OHC = OHC.T

    OHC_vec = np.ravel(OHC)[okR_hycom_4]

    Tmean_max_hycom_exp[n] = np.nanmax(Tmean_dens_crit_vec)
    Tmean_mean_hycom_exp[n] = np.nanmean(Tmean_dens_crit_vec)
    Tmean_min_hycom_exp[n] = np.nanmin(Tmean_dens_crit_vec)

    OHC_max_hycom_exp[n] = np.nanmax(OHC_vec)
    OHC_mean_hycom_exp[n] = np.nanmean(OHC_vec)
    OHC_min_hycom_exp[n] = np.nanmin(OHC_vec)

#%% Start figure all panels

fig, ax = plt.subplots(figsize=(5, 15))
grid = plt.GridSpec(4, 1, wspace=0.1, hspace=0.3,left=0.05,right=0.95)
fig.suptitle('Intensity Forecast Cycle '+cycle,fontsize=16)

#%% Figure forecasted intensity models vs best intensity

okt = np.logical_and(time_best_track >= tini,time_best_track <= tend)

lead_time_pom_oper = np.arange(0,129,3)

ax1 = plt.subplot(grid[0,0])
ax1.plot(lead_time_pom_oper,np.tile(64,len(lead_time_pom_oper)),'-',color='grey',alpha=0.5,linewidth=1)
ax1.plot(lead_time_pom_oper,np.tile(83,len(lead_time_pom_oper)),'-',color='grey',alpha=0.5,linewidth=1)
ax1.plot(lead_time_pom_oper,np.tile(96,len(lead_time_pom_oper)),'-',color='grey',alpha=0.5,linewidth=1)
ax1.plot(lead_time_pom_oper,np.tile(113,len(lead_time_pom_oper)),'-',color='grey',alpha=0.5,linewidth=1)
ax1.plot(lead_time_pom_oper,np.tile(137,len(lead_time_pom_oper)),'-',color='grey',alpha=0.5,linewidth=1)

ax1.plot(lead_time_pom_oper[::2],wind_int_kt[okt],'o-k',label='Best')
ax1.plot(lead_time_pom_oper,max_wind_10m_hwrf_pom19_oper,'X-',color='mediumorchid',label='HWRF2019-POM (IC Clim.)',markeredgecolor='k',markersize=7)
ax1.plot(lead_time_pom_oper,max_wind_10m_hwrf_pom20_exp,'^-',color='teal',label='HWRF2020-POM (IC RTOFS)',markeredgecolor='k',markersize=7)
ax1.plot(lead_time_pom_oper,max_wind_10m_hwrf_hycom20_exp,'H-',color='darkorange',label='HWRF2020-HYCOM (IC RTOFS)',markeredgecolor='k',markersize=7)

ax1.tick_params(which='major', width=2)
ax1.tick_params(which='major', length=7)
ax1.tick_params(which='minor', length=4, color='k')

ax1.xaxis.set_major_locator(MultipleLocator(12))
ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax1.xaxis.set_minor_locator(MultipleLocator(3))
ax1.xaxis.set_ticks(np.arange(0,86,12))
ax1.xaxis.set_ticklabels([''])
#ax1.xaxis.set_ticklabels(['28-Aug \n 0','\n 12','29-Aug \n 24','\n 36','30-Aug \n 48',\
#                          '\n 60','31-Aug \n 72','\n 84','01-Sep \n 96','\n 108','02-Sep \n 120'])
#plt.xlabel('Forecast Lead Time (Hr)',fontsize=14,labelpad=10)
legend = ax1.legend(loc='upper left',fontsize=10,bbox_to_anchor=(0,1.5),framealpha=1)
ax1.set_ylim([20,150])
ax1.set_xlim([0,86])
#plt.title('Intensity Forecast',fontsize=14,loc='right')
ax1.set_ylabel('Max 10m Wind (kt)',fontsize=12)

ax2 = ax1.twinx()
plt.ylim([20,150])
yticks = [64,83,96,113,137]
plt.yticks(yticks,['Cat 1','Cat 2','Cat 3','Cat 4','Cat 5'])
plt.text(82,155,'(a)',fontsize=14)

#%%  Figure Area Integrated Sea Surface Heat Loss vs forecasted lead time

aif_pom_oper = str(format(acum_ai_sshl_pom_oper,'0.2E'))
aif_pom_exp = str(format(acum_ai_sshl_pom_exp,'0.2E'))
aif_hycom_exp = str(format(acum_ai_sshl_hycom_exp,'0.2E'))

perc_sshl_pom_oper = (acum_ai_sshl_hycom_exp - acum_ai_sshl_pom_oper)*100/acum_ai_sshl_hycom_exp
perc_sshl_pom_exp = (acum_ai_sshl_hycom_exp - acum_ai_sshl_pom_exp)*100/acum_ai_sshl_hycom_exp

ax2 = plt.subplot(grid[1,0])
plt.plot(lead_time_pom_oper[1:29],sshla_R4_mean_pom_oper[1:29],'X-',color='mediumorchid',\
         label=aif_pom_oper+ ' KJ',markeredgecolor='k',markersize=7)
#ax.fill_between(lead_time_pom_oper[1:29],sshla_R4_min_pom_oper[1:29],sshla_R4_max_pom_oper[1:29],color='mediumorchid',alpha=0.1)

plt.plot(lead_time_pom_exp[1:29],sshla_R4_mean_pom_exp[1:29],'^-',color='teal',\
         label=aif_pom_exp+ ' KJ',markeredgecolor='k',markersize=7)
#ax2.fill_between(lead_time_pom_exp[1:29],sshla_R4_min_pom_exp[1:29],sshla_R4_max_pom_exp[1:29],color='teal',alpha=0.1)

plt.plot(lead_time_hycom_exp[1:29],sshla_R4_mean_hycom_exp[1:29],'H-',color='darkorange',\
         label=aif_hycom_exp+ ' KJ',markeredgecolor='k',markersize=7)
#ax2.fill_between(lead_time_hycom_exp[1:29],sshla_R4_min_hycom_exp[1:29],sshla_R4_max_hycom_exp[1:29],color='darkorange',alpha=0.1)

plt.title('Mean Sea Surface Heat Loss within 4Rmax',fontsize=12)
plt.ylabel('$KJ/cm^2$',fontsize=12)
#plt.xlabel('Forecasted Lead Time (Hours)',fontsize=14)
#plt.legend(loc='lower right',fontsize=9)
plt.xlim([0,86])
#plt.ylim([3*10**14,7*10**14])

ax2.tick_params(which='major', width=2)
ax2.tick_params(which='major', length=7)
ax2.tick_params(which='minor', length=4, color='k')

ax2.xaxis.set_major_locator(MultipleLocator(12))
ax2.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax2.xaxis.set_minor_locator(MultipleLocator(3))
ax2.xaxis.set_ticks(np.arange(0,86,12))
ax2.xaxis.set_ticklabels([''])
#ax.xaxis.set_ticklabels(['28-Aug \n 0','\n 12','29-Aug \n 24','\n 36','30-Aug \n 48',\
#                          '\n 60','31-Aug \n 72','\n 84']) #,'01-Sep \n 96','\n 108','02-Sep \n 120'])
plt.text(82,0.88,'(b)',fontsize=14)

#%%
ax1 = plt.subplot(grid[2, 0])
plt.title('Mean Mixed Layer Temparature within 4Rmax', fontsize=12)
plt.plot(lead_time_pom_oper[::2],Tmean_mean_pom_oper,'X-',color='mediumorchid',label='HWRF2019-POM (IC clim.)',markeredgecolor='k',markersize=7)
ax1.fill_between(lead_time_pom_oper[::2],Tmean_min_pom_oper,Tmean_max_pom_oper,color='mediumorchid',alpha=0.1)

plt.plot(lead_time_pom_exp[::2],Tmean_mean_pom_exp,'^-',color='teal',label='HWRF2020-POM (RTOFS)',markeredgecolor='k',markersize=7)
ax1.fill_between(lead_time_pom_exp[::2],Tmean_min_pom_exp,Tmean_max_pom_exp,color='teal',alpha=0.1)

plt.plot(lead_time_hycom_exp[::2],Tmean_mean_hycom_exp,'H-',color='darkorange',label='HWRF2020-HYCOM (RTOFS)',markeredgecolor='k',markersize=7)
ax1.fill_between(lead_time_hycom_exp[::2],Tmean_min_hycom_exp,Tmean_max_hycom_exp,color='darkorange',alpha=0.1)

ax1.xaxis.set_major_locator(MultipleLocator(12))
ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax1.xaxis.set_minor_locator(MultipleLocator(3))
ax1.xaxis.set_ticks(np.arange(0,126,12))
#ax.xaxis.set_ticklabels(['28-Aug \n 0','\n 12','29-Aug \n 24','\n 36','30-Aug \n 48',\
#                          '\n 60','31-Aug \n 72','\n 84']) #,'01-Sep \n 96','\n 108','02-Sep \n 120'])
ax1.xaxis.set_ticklabels(' ')

ax1.tick_params(which='major', width=2)
ax1.tick_params(which='major', length=7)
ax1.tick_params(which='minor', length=4, color='k')

#plt.ylim(27.0,29.5)
plt.ylim(27.5,29.5)
plt.xlim(0,86)
plt.ylabel('$^oC$',fontsize=12)
#plt.xlabel('Forecasted Lead Time (Hours)',fontsize=14)
#plt.legend(loc='lower left')
plt.text(82,29.6,'(c)',fontsize=14)


#%%
ax2 = plt.subplot(grid[3, 0])
plt.plot(lead_time_pom_oper[::2],OHC_mean_pom_oper,'X-',color='mediumorchid',label='HWRF2019-POM (IC clim.)',markeredgecolor='k',markersize=7)
ax2.fill_between(lead_time_pom_oper[::2],OHC_min_pom_oper,OHC_max_pom_oper,color='mediumorchid',alpha=0.1)

plt.plot(lead_time_pom_oper[::2],OHC_mean_pom_exp,'^-',color='teal',label='HWRF2020-POM (RTOFS)',markeredgecolor='k',markersize=7)
ax2.fill_between(lead_time_pom_oper[::2],OHC_min_pom_exp,OHC_max_pom_exp,color='teal',alpha=0.1)

plt.plot(lead_time_pom_oper[::2],OHC_mean_hycom_exp,'H-',color='darkorange',label='HWRF2020-HYCOM (RTOFS)',markeredgecolor='k',markersize=7)
ax2.fill_between(lead_time_pom_oper[::2],OHC_min_hycom_exp,OHC_max_hycom_exp,color='darkorange',alpha=0.1)

plt.plot(lead_time_pom_oper[::2],np.tile(60,len(lead_time_pom_oper[::2])),'--k')

ax2.xaxis.set_major_locator(MultipleLocator(12))
ax2.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax2.xaxis.set_minor_locator(MultipleLocator(3))
ax2.xaxis.set_ticks(np.arange(0,126,12))
ax2.xaxis.set_ticklabels(['28-Aug \n 0','\n 12','29-Aug \n 24','\n 36','30-Aug \n 48',\
                          '\n 60','31-Aug \n 72','\n 84']) #,'01-Sep \n 96','\n 108','02-Sep \n 120'])

ax2.tick_params(which='major', width=2)
ax2.tick_params(which='major', length=7)
ax2.tick_params(which='minor', length=4, color='k')

plt.title('Mean Ocean Heat Content within 4Rmax',fontsize=12)
#plt.ylim(0,100)
plt.ylim(25,100)
plt.xlim(0,86)
plt.ylabel('$Kj/cm^2$',fontsize=12)
plt.xlabel('Forecasted Lead Time (Hours)',fontsize=12)
#plt.legend(loc='lower left')
plt.text(82,103,'(d)',fontsize=14)

file = folder_fig + 'Mean_SSHL_MLT_OHC_final_fig' + '_' + cycle
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)

#%%%%%%%%%%%%%%%%%%%

#%%  Figure Area Integrated Sea Surface Heat Loss vs forecasted lead time

aif_pom_oper = str(format(acum_ai_sshl_pom_oper,'0.2E'))
aif_pom_exp = str(format(acum_ai_sshl_pom_exp,'0.2E'))
aif_hycom_exp = str(format(acum_ai_sshl_hycom_exp,'0.2E'))

perc_sshl_pom_oper = (acum_ai_sshl_hycom_exp - acum_ai_sshl_pom_oper)*100/acum_ai_sshl_hycom_exp
perc_sshl_pom_exp = (acum_ai_sshl_hycom_exp - acum_ai_sshl_pom_exp)*100/acum_ai_sshl_hycom_exp

fig, ax = plt.subplots(figsize=(7, 5))

plt.plot(lead_time_pom_oper[1:29],area_int_sshl_pom_oper[1:29],'X-',color='mediumorchid',\
         label=aif_pom_oper+ ' KJ',markeredgecolor='k',markersize=7)
ax.fill_between(lead_time_pom_oper[1:29],area_int_sshl_pom_oper[1:29],color='mediumorchid',alpha=0.1)

plt.plot(lead_time_pom_exp[1:29],area_int_sshl_pom_exp[1:29],'^-',color='teal',\
         label=aif_pom_exp+ ' KJ',markeredgecolor='k',markersize=7)
ax.fill_between(lead_time_pom_exp[1:29],area_int_sshl_pom_exp[1:29],color='teal',alpha=0.1)

plt.plot(lead_time_hycom_exp[1:29],area_int_sshl_hycom_exp[1:29],'H-',color='darkorange',\
         label=aif_hycom_exp+ ' KJ',markeredgecolor='k',markersize=7)
ax.fill_between(lead_time_hycom_exp[1:29],area_int_sshl_hycom_exp[1:29],color='darkorange',alpha=0.1)

plt.title('Area Integrated Sea Surface Heat Loss',fontsize=16)
plt.ylabel('$KJ$',fontsize=12)
plt.xlabel('Forecasted Lead Time (Hours)',fontsize=14)
plt.legend(loc='lower right',fontsize=14)
plt.xlim([0,86])
plt.ylim([3*10**14,7*10**14])

ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')

ax.xaxis.set_major_locator(MultipleLocator(12))
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,86,12))
ax.xaxis.set_ticklabels([''])
ax.xaxis.set_ticklabels(['28-Aug \n 0','\n 12','29-Aug \n 24','\n 36','30-Aug \n 48',\
                          '\n 60','31-Aug \n 72','\n 84']) #,'01-Sep \n 96','\n 108','02-Sep \n 120'])

file = folder_fig + 'Int_SSHL_final_fig' + '_' + cycle
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)

#%%%%%%%%%%%  Figure Accumulated Area Integrated Sea Surface Heat Loss vs forecasted lead time

aif_pom_oper = str(format(acum_ai_sshl_pom_oper,'0.2E'))
aif_pom_exp = str(format(acum_ai_sshl_pom_exp,'0.2E'))
aif_hycom_exp = str(format(acum_ai_sshl_hycom_exp,'0.2E'))

perc_sshl_pom_oper = (acum_ai_sshl_hycom_exp - acum_ai_sshl_pom_oper)*100/acum_ai_sshl_hycom_exp
perc_sshl_pom_exp = (acum_ai_sshl_hycom_exp - acum_ai_sshl_pom_exp)*100/acum_ai_sshl_hycom_exp

fig, ax = plt.subplots(figsize=(7, 5))

plt.plot(lead_time_pom_oper[1:29],np.cumsum(area_int_sshl_pom_oper[1:29]),'X-',color='mediumorchid',\
         label=aif_pom_oper+ ' KJ',markeredgecolor='k',markersize=7)

plt.plot(lead_time_pom_exp[1:29],np.cumsum(area_int_sshl_pom_exp[1:29]),'^-',color='teal',\
         label=aif_pom_exp+ ' KJ',markeredgecolor='k',markersize=7)


plt.plot(lead_time_hycom_exp[1:29],np.cumsum(area_int_sshl_hycom_exp[1:29]),'H-',color='darkorange',\
         label=aif_hycom_exp+ ' KJ',markeredgecolor='k',markersize=7)

plt.title('Area Integrated Sea Surface Heat Loss',fontsize=16)
plt.ylabel('$KJ$',fontsize=12)
plt.xlabel('Forecasted Lead Time (Hours)',fontsize=14)
#plt.legend(loc='lower right',fontsize=9)
plt.xlim([0,86])
#plt.ylim([3*10**14,7*10**14])

ax.tick_params(which='major', width=2)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')

ax.xaxis.set_major_locator(MultipleLocator(12))
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(3))
ax.xaxis.set_ticks(np.arange(0,86,12))
ax.xaxis.set_ticklabels([''])
ax.xaxis.set_ticklabels(['28-Aug \n 0','\n 12','29-Aug \n 24','\n 36','30-Aug \n 48',\
                          '\n 60','31-Aug \n 72','\n 84']) #,'01-Sep \n 96','\n 108','02-Sep \n 120'])

file = folder_fig + 'Accum_Int_SSHL_final_fig' + '_' + cycle
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)
