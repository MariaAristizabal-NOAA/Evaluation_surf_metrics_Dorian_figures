#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 15:08:51 2020

@author: aristizabal
"""

#%% User input

home_folder = '/home/'
#home_folder = '/Volumes'

#lon_lim = [-100.0,-55.0]
#lat_lim = [10.0,45.0]

lon_lim = [-85.0,-60.0]
lat_lim = [15.0,35.0]

# forecasting cycle to be used
cycle = '2019082800'

#Time window
#date_ini = '2019/08/28/00/00'
#date_end = '2019/09/02/06/00'

# Bathymetry file
bath_file = home_folder+'aristizabal/bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'

# KMZ file
kmz_file_Dorian = home_folder+'aristizabal/KMZ_files/al052019_best_track-5.kmz'

# Folder where to save figure
folder_fig = '/www/web/rucool/aristizabal/Figures/'

# folder nc files POM
folder_pom19 =  home_folder+'aristizabal/HWRF2019_POM_Dorian/'
folder_pom20 =  home_folder+'aristizabal/HWRF2020_POM_Dorian/'

# folde HWRF2020_HYCOM
folder_hycom20 = home_folder+'aristizabal/HWRF2020_HYCOM_Dorian/'

###################\
# folder nc files POM
folder_pom_oper = folder_pom19 + 'HWRF2019_POM_dorian05l.' + cycle + '_pom_files_oper/'
folder_pom_exp = folder_pom20 + 'HWRF2020_POM_dorian05l.'  + cycle + '_pom_files_exp/'
prefix_pom = 'dorian05l.' + cycle + '.pom.00'

pom_grid_oper = folder_pom_oper + 'dorian05l.' + cycle + '.pom.grid.nc'
pom_grid_exp = folder_pom_exp + 'dorian05l.' + cycle + '.pom.grid.nc'

# Dorian track files
hwrf_pom_track_oper = folder_pom_oper + 'dorian05l.' + cycle + '.trak.hwrf.atcfunix'
hwrf_pom_track_exp = folder_pom_exp + 'dorian05l.' + cycle + '.trak.hwrf.atcfunix'

##################
# folder ab files HYCOM
folder_hycom_exp = folder_hycom20 + 'HWRF2020_HYCOM_dorian05l.' + cycle + '_hycom_files_exp/'
prefix_hycom = 'dorian05l.' + cycle + '.hwrf_rtofs_hat10_3z'

Dir_HMON_HYCOM = home_folder+'aristizabal/HWRF_RTOFS_Michael_2018/HWRF-Hycom-WW3_exp_Michael/'
# RTOFS grid file name
hycom_grid_exp = Dir_HMON_HYCOM + 'hwrf_rtofs_hat10.basin.regional.grid'

# Dorian track files
hwrf_hycom_track_exp = folder_hycom_exp + 'dorian05l.' + cycle + '.trak.hwrf.atcfunix'

####################
# folder nc files hwrf
folder_hwrf_pom19_oper = folder_pom19 + 'HWRF2019_POM_dorian05l.' + cycle + '_grb2_to_nc_oper/'
folder_hwrf_pom20_exp = folder_pom20 + 'HWRF2020_POM_dorian05l.' + cycle + '_grb2_to_nc_exp/'

# folder nc files hwrf
folder_hwrf_hycom20_exp = folder_hycom20 + 'HWRF2020_HYCOM_dorian05l.' + cycle + '_grb2_to_nc_exp/'

#%%

from matplotlib import pyplot as plt
import numpy as np
import xarray as xr
from netCDF4 import Dataset
from datetime import datetime, timedelta
import matplotlib.dates as mdates
import sys
import os
import os.path
import glob
import seawater
import cmocean

#sys.path.append('/Users/aristizabal/Desktop/MARACOOS_project/NCEP_scripts')
sys.path.append('/home/aristizabal/NCEP_scripts')
from utils4HYCOM import readgrids
from utils4HYCOM import readBinz

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)
plt.rc('legend',fontsize=12)

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

ncbath = Dataset(bath_file)
bath_lat = ncbath.variables['lat'][:]
bath_lon = ncbath.variables['lon'][:]
bath_elev = ncbath.variables['elevation'][:]

oklatbath = np.logical_and(bath_lat >= lat_lim[0],bath_lat <= lat_lim[-1])
oklonbath = np.logical_and(bath_lon >= lon_lim[0],bath_lon <= lon_lim[-1])

bath_latsub = bath_lat[oklatbath]
bath_lonsub = bath_lon[oklonbath]
bath_elevs = bath_elev[oklatbath,:]
bath_elevsub = bath_elevs[:,oklonbath]

#%% Get list HWRF files

HWRF_POM_oper = sorted(glob.glob(os.path.join(folder_hwrf_pom19_oper,'*.nc')))[::2]
HWRF_POM_exp = sorted(glob.glob(os.path.join(folder_hwrf_pom20_exp,'*.nc')))[::2]
HWRF_HYCOM_exp = sorted(glob.glob(os.path.join(folder_hwrf_hycom20_exp,'*.nc')))[::2]

#%% Get list ocean model files

POM_oper = sorted(glob.glob(os.path.join(folder_pom_oper,prefix_pom+'*.nc')))
POM_exp = sorted(glob.glob(os.path.join(folder_pom_exp,prefix_pom+'*.nc')))
HYCOM_exp = sorted(glob.glob(os.path.join(folder_hycom_exp,prefix_hycom+'*.a')))

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

###################### Start figure all panels MLT ##############################

fig, ax = plt.subplots(figsize=(15, 15))
grid = plt.GridSpec(3, 9, wspace=0.1, hspace=0.3,left=0.05,right=0.95)
#rid = plt.GridSpec(3, 6, wspace=0.1, hspace=0.3,left=0.01,right=0.99)
fig.suptitle('Mixed Layer Temperature',fontsize=20)

RR_norm_bins = np.arange(0,8.1,0.5)
RR_norm_bins_mid = RR_norm_bins[0:-1] + np.diff(RR_norm_bins)/2

Tmean_max_pom_oper = np.empty((3,len(RR_norm_bins)))
Tmean_max_pom_oper[:] = np.nan
Tmean_min_pom_oper = np.empty((3,len(RR_norm_bins)))
Tmean_min_pom_oper[:] = np.nan
Tmean_mean_pom_oper = np.empty((3,len(RR_norm_bins)))
Tmean_mean_pom_oper[:] = np.nan

time = []

for N,indx in enumerate(np.asarray([3,11,14])):
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
    Rmax = seawater.dist([lat_forec_track_pom_oper[::2][indx],lat_maxwind],\
                         [lon_forec_track_pom_oper[::2][indx],lon_maxwind],'km')[0][0]

    xlim = [lon_forec_track_pom_oper[::2][indx]-4,lon_forec_track_pom_oper[::2][indx]+4]
    ylim = [lat_forec_track_pom_oper[::2][indx]-4,lat_forec_track_pom_oper[::2][indx]+4]

    oklon = np.where(np.logical_and(lon_hwrf>xlim[0],lon_hwrf<xlim[1]))[0]
    oklat = np.where(np.logical_and(lat_hwrf>ylim[0],lat_hwrf<ylim[1]))[0]

    meshlon_lat = np.meshgrid(lon_hwrf[oklon],lat_hwrf[oklat])
    #meshlat = np.meshgrid(lat_hwrf[oklat],lon_hwrf[oklon])

    eye_lon = np.tile(lon_forec_track_pom_oper[::2][indx],meshlon_lat[0].shape[1])
    eye_lat = np.tile(lat_forec_track_pom_oper[::2][indx],meshlon_lat[0].shape[0])
    eye_matrix = np.meshgrid(eye_lon,eye_lat)

    lat_lon_matrix = np.stack((np.ravel(meshlon_lat[0]),np.ravel(meshlon_lat[1])),axis=1).T
    eye_lat_lon_matrix = np.stack((np.ravel(eye_matrix[0]),np.ravel(eye_matrix[1])),axis=1).T

    R = np.empty(lat_lon_matrix.shape[1])
    R[:] = np.nan
    for i in np.arange(lat_lon_matrix.shape[1]):
        R[i] = seawater.dist([lat_forec_track_pom_oper[::2][indx],lat_lon_matrix[1,i]],\
                             [lon_forec_track_pom_oper[::2][indx],lat_lon_matrix[0,i]],'km')[0][0]

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

    # POM oper
    print(POM_oper[indx])
    pom = xr.open_dataset(POM_oper[indx])
    tpom = np.asarray(pom['time'][:])
    temp_pom = np.asarray(pom['t'][0,:,:,:])
    salt_pom = np.asarray(pom['s'][0,:,:,:])
    dens_pom = np.asarray(pom['rho'][0,:,:,:]) * 1000 + 1000
    dens_pom[dens_pom==1000.0] = np.nan

    oklon_pom = np.where(np.logical_and(lon_pom_oper[0,:]>xlim[0],lon_pom_oper[0,:]<xlim[1]))[0]
    oklat_pom = np.where(np.logical_and(lat_pom_oper[:,0]>ylim[0],lat_pom_oper[:,0]<ylim[1]))[0]

    meshlon_lat_pom = np.meshgrid(lon_pom_oper[0,oklon_pom],lat_pom_oper[oklat_pom,0])

    lat_lon_matrix_pom = np.stack((np.ravel(meshlon_lat_pom[0]),np.ravel(meshlon_lat_pom[1])),axis=1).T

    R_pom = np.empty(lat_lon_matrix_pom.shape[1])
    R_pom[:] = np.nan
    for i in np.arange(lat_lon_matrix_pom.shape[1]):
        R_pom[i] = seawater.dist([lat_forec_track_pom_oper[::2][indx],lat_lon_matrix_pom[1,i]],\
                             [lon_forec_track_pom_oper[::2][indx],lat_lon_matrix_pom[0,i]],'km')[0][0]

    R_pom_norm = R_pom/Rmax
    okR_pom_8 = R_pom_norm <= 8
    okR_pom_larger_8 = R_pom_norm > 8

    RR_pom_norm = R_pom_norm[okR_pom_8]

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

    Tmean_dens_crit_vec = np.ravel(Tmean_dens_crit)[okR_pom_8]

    okRR_pom = np.argsort(RR_pom_norm)

    for n in np.arange(1,len(RR_norm_bins)):
        print(n)
        okR = np.logical_and(RR_pom_norm >= RR_norm_bins[n-1],RR_pom_norm <= RR_norm_bins[n])
        if len(np.where(okR)[0])==0:
            Tmean_max_pom_oper[N,n] = np.nan
            Tmean_mean_pom_oper[N,n] = np.nan
            Tmean_min_pom_oper[N,n] = np.nan
        else:
            Tmean_max_pom_oper[N,n] = np.nanmax(Tmean_dens_crit_vec[okR])
            Tmean_mean_pom_oper[N,n] = np.nanmean(Tmean_dens_crit_vec[okR])
            Tmean_min_pom_oper[N,n] = np.nanmin(Tmean_dens_crit_vec[okR])

    kw = dict(levels=np.arange(27,29.6,0.2))
    #fig,ax = plt.subplots()
    ax = plt.subplot(grid[N, 0:2])
    plt.contourf(lon_pom_oper[0,oklon_pom],lat_pom_oper[oklat_pom,0],Tmean_dens_crit,cmap=cmocean.cm.thermal,\
                 **kw,vmin=27,vmax=29.5)
    #c=plt.colorbar()
    #c.set_label('$^oC$',rotation=90, labelpad=15, fontsize=14)
    plt.contour(lon_pom_oper[0,oklon_pom],lat_pom_oper[oklat_pom,0],Tmean_dens_crit,\
                colors='k',linestyle='-',**kw,alpha=0.3)
    plt.plot(eye_lon[0],eye_lat[0],'ok',markersize=7,markeredgecolor='lawngreen')
    plt.plot(lon_forec_track_pom_oper,lat_forec_track_pom_oper,'.-k')
    plt.plot(lat_lon_matrix[0][okR1],lat_lon_matrix[1][okR1],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR2],lat_lon_matrix[1][okR2],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR3],lat_lon_matrix[1][okR3],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR4],lat_lon_matrix[1][okR4],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR5],lat_lon_matrix[1][okR5],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR6],lat_lon_matrix[1][okR6],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR7],lat_lon_matrix[1][okR7],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR8],lat_lon_matrix[1][okR8],'.k',markersize=0.5)
    plt.plot(lon_maxwind,lat_maxwind,'*k',markersize=10,markeredgecolor='lawngreen')
    plt.axis('scaled')
    plt.xlim(np.min(lat_lon_matrix[0][okR_8])-0.1,np.max(lat_lon_matrix[0][okR_8])+0.1)
    plt.ylim(np.min(lat_lon_matrix[1][okR_8])-0.1,np.max(lat_lon_matrix[1][okR_8])+0.1)
    if N==0:
        plt.title('HWRF2019-POM \n'+ str(time_hwrf[0])[0:13],fontsize=14)
    else:
        plt.title(str(time_hwrf[0])[0:13],fontsize=14)
    if N==0:
        plt.text(-68.5, 17.5, 'T1',fontsize=16)
    if N==1:
        plt.text(-75.7, 25.5, 'T2',fontsize=16)
    if N==2:
        plt.text(-78.2, 27.5, 'T3',fontsize=16)


#%% MLT and OHC under Dorian POM experimental

Tmean_max_pom_exp = np.empty((3,len(RR_norm_bins)))
Tmean_max_pom_exp[:] = np.nan
Tmean_min_pom_exp = np.empty((3,len(RR_norm_bins)))
Tmean_min_pom_exp[:] = np.nan
Tmean_mean_pom_exp = np.empty((3,len(RR_norm_bins)))
Tmean_mean_pom_exp[:] = np.nan

time = []

for N,indx in enumerate(np.asarray([3,11,14])):
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

    time.append(np.asarray(HWRF.variables['time'][:]))

    wind_int = np.sqrt(UGRD_hwrf**2 + VGRD_hwrf**2)
    max_wind= np.max(wind_int)
    okwind = np.where(wind_int == max_wind)
    lat_maxwind = lat_hwrf[okwind[0][0]]
    lon_maxwind = lon_hwrf[okwind[1][0]]
    Rmax = seawater.dist([lat_forec_track_pom_exp[::2][indx],lat_maxwind],\
                         [lon_forec_track_pom_exp[::2][indx],lon_maxwind],'km')[0][0]

    xlim = [lon_forec_track_pom_exp[::2][indx]-4,lon_forec_track_pom_exp[::2][indx]+4]
    ylim = [lat_forec_track_pom_exp[::2][indx]-4,lat_forec_track_pom_exp[::2][indx]+4]

    oklon = np.where(np.logical_and(lon_hwrf>xlim[0],lon_hwrf<xlim[1]))[0]
    oklat = np.where(np.logical_and(lat_hwrf>ylim[0],lat_hwrf<ylim[1]))[0]

    meshlon_lat = np.meshgrid(lon_hwrf[oklon],lat_hwrf[oklat])
    #meshlat = np.meshgrid(lat_hwrf[oklat],lon_hwrf[oklon])

    eye_lon = np.tile(lon_forec_track_pom_exp[::2][indx],meshlon_lat[0].shape[1])
    eye_lat = np.tile(lat_forec_track_pom_exp[::2][indx],meshlon_lat[0].shape[0])
    eye_matrix = np.meshgrid(eye_lon,eye_lat)

    lat_lon_matrix = np.stack((np.ravel(meshlon_lat[0]),np.ravel(meshlon_lat[1])),axis=1).T
    eye_lat_lon_matrix = np.stack((np.ravel(eye_matrix[0]),np.ravel(eye_matrix[1])),axis=1).T

    R = np.empty(lat_lon_matrix.shape[1])
    R[:] = np.nan
    for i in np.arange(lat_lon_matrix.shape[1]):
        R[i] = seawater.dist([lat_forec_track_pom_exp[::2][indx],lat_lon_matrix[1,i]],\
                             [lon_forec_track_pom_exp[::2][indx],lat_lon_matrix[0,i]],'km')[0][0]

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

    # POM exp
    print(POM_exp[indx])
    pom = xr.open_dataset(POM_exp[indx])
    tpom = np.asarray(pom['time'][:])
    temp_pom = np.asarray(pom['t'][0,:,:,:])
    salt_pom = np.asarray(pom['s'][0,:,:,:])
    dens_pom = np.asarray(pom['rho'][0,:,:,:]) * 1000 + 1000
    dens_pom[dens_pom==1000.0] = np.nan

    oklon_pom = np.where(np.logical_and(lon_pom_exp[0,:]>xlim[0],lon_pom_exp[0,:]<xlim[1]))[0]
    oklat_pom = np.where(np.logical_and(lat_pom_exp[:,0]>ylim[0],lat_pom_exp[:,0]<ylim[1]))[0]

    meshlon_lat_pom = np.meshgrid(lon_pom_exp[0,oklon_pom],lat_pom_exp[oklat_pom,0])

    lat_lon_matrix_pom = np.stack((np.ravel(meshlon_lat_pom[0]),np.ravel(meshlon_lat_pom[1])),axis=1).T

    R_pom = np.empty(lat_lon_matrix_pom.shape[1])
    R_pom[:] = np.nan
    for i in np.arange(lat_lon_matrix_pom.shape[1]):
        R_pom[i] = seawater.dist([lat_forec_track_pom_exp[::2][indx],lat_lon_matrix_pom[1,i]],\
                             [lon_forec_track_pom_exp[::2][indx],lat_lon_matrix_pom[0,i]],'km')[0][0]

    R_pom_norm = R_pom/Rmax
    okR_pom_8 = R_pom_norm <= 8

    RR_pom_norm = R_pom_norm[okR_pom_8]

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

    Tmean_dens_crit_vec = np.ravel(Tmean_dens_crit)[okR_pom_8]

    okRR_pom = np.argsort(RR_pom_norm)

    for n in np.arange(1,len(RR_norm_bins)):
        print(n)
        okR = np.logical_and(RR_pom_norm >= RR_norm_bins[n-1],RR_pom_norm <= RR_norm_bins[n])
        if len(np.where(okR)[0])==0:
            Tmean_max_pom_exp[N,n] = np.nan
            Tmean_mean_pom_exp[N,n] = np.nan
            Tmean_min_pom_exp[N,n] = np.nan
        else:
            Tmean_max_pom_exp[N,n] = np.nanmax(Tmean_dens_crit_vec[okR])
            Tmean_mean_pom_exp[N,n] = np.nanmean(Tmean_dens_crit_vec[okR])
            Tmean_min_pom_exp[N,n] = np.nanmin(Tmean_dens_crit_vec[okR])

    kw = dict(levels=np.arange(27,29.6,0.2))
    ax = plt.subplot(grid[N, 2:4])
    cs=plt.contourf(lon_pom_exp[0,oklon_pom],lat_pom_exp[oklat_pom,0],Tmean_dens_crit,cmap=cmocean.cm.thermal,\
                 **kw,vmin=27,vmax=29.5)
    plt.contour(lon_pom_exp[0,oklon_pom],lat_pom_exp[oklat_pom,0],Tmean_dens_crit,\
                colors='k',linestyle='-',**kw,alpha=0.3)
    plt.plot(eye_lon[0],eye_lat[0],'ok',markersize=7,markeredgecolor='lawngreen')
    plt.plot(lon_forec_track_pom_exp,lat_forec_track_pom_exp,'.-k')
    plt.plot(lat_lon_matrix[0][okR1],lat_lon_matrix[1][okR1],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR2],lat_lon_matrix[1][okR2],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR3],lat_lon_matrix[1][okR3],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR4],lat_lon_matrix[1][okR4],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR5],lat_lon_matrix[1][okR5],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR6],lat_lon_matrix[1][okR6],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR7],lat_lon_matrix[1][okR7],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR8],lat_lon_matrix[1][okR8],'.k',markersize=0.5)
    plt.plot(lon_maxwind,lat_maxwind,'*k',markersize=10,markeredgecolor='lawngreen')
    plt.axis('scaled')
    plt.xlim(np.min(lat_lon_matrix[0][okR_8])-0.1,np.max(lat_lon_matrix[0][okR_8])+0.1)
    plt.ylim(np.min(lat_lon_matrix[1][okR_8])-0.1,np.max(lat_lon_matrix[1][okR_8])+0.1)
    if N==0:
        plt.title('HWRF2020-POM \n' + str(time_hwrf[0])[0:13],fontsize=14)
    else:
        plt.title(str(time_hwrf[0])[0:13],fontsize=14)

#%% MLT and OHC under Dorian HYCOM experimental

Tmean_max_hycom_exp = np.empty((3,len(RR_norm_bins)))
Tmean_max_hycom_exp[:] = np.nan
Tmean_min_hycom_exp = np.empty((3,len(RR_norm_bins)))
Tmean_min_hycom_exp[:] = np.nan
Tmean_mean_hycom_exp = np.empty((3,len(RR_norm_bins)))
Tmean_mean_hycom_exp[:] = np.nan

time = []

afiles = sorted(glob.glob(os.path.join(folder_hycom_exp,prefix_hycom+'*.a')))

for N,indx in enumerate(np.asarray([3,11,14])):
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

    time.append(np.asarray(HWRF.variables['time'][:]))

    wind_int = np.sqrt(UGRD_hwrf**2 + VGRD_hwrf**2)
    max_wind= np.max(wind_int)
    okwind = np.where(wind_int == max_wind)
    lat_maxwind = lat_hwrf[okwind[0][0]]
    lon_maxwind = lon_hwrf[okwind[1][0]]
    Rmax = seawater.dist([lat_forec_track_hycom_exp[::2][indx],lat_maxwind],\
                         [lon_forec_track_hycom_exp[::2][indx],lon_maxwind],'km')[0][0]

    xlim = [lon_forec_track_hycom_exp[::2][indx]-4,lon_forec_track_hycom_exp[::2][indx]+4]
    ylim = [lat_forec_track_hycom_exp[::2][indx]-4,lat_forec_track_hycom_exp[::2][indx]+4]

    oklon = np.where(np.logical_and(lon_hwrf>xlim[0],lon_hwrf<xlim[1]))[0]
    oklat = np.where(np.logical_and(lat_hwrf>ylim[0],lat_hwrf<ylim[1]))[0]

    meshlon_lat = np.meshgrid(lon_hwrf[oklon],lat_hwrf[oklat])

    eye_lon = np.tile(lon_forec_track_hycom_exp[::2][indx],meshlon_lat[0].shape[1])
    eye_lat = np.tile(lat_forec_track_hycom_exp[::2][indx],meshlon_lat[0].shape[0])
    eye_matrix = np.meshgrid(eye_lon,eye_lat)

    lat_lon_matrix = np.stack((np.ravel(meshlon_lat[0]),np.ravel(meshlon_lat[1])),axis=1).T
    eye_lat_lon_matrix = np.stack((np.ravel(eye_matrix[0]),np.ravel(eye_matrix[1])),axis=1).T

    R = np.empty(lat_lon_matrix.shape[1])
    R[:] = np.nan
    for i in np.arange(lat_lon_matrix.shape[1]):
        R[i] = seawater.dist([lat_forec_track_hycom_exp[::2][indx],lat_lon_matrix[1,i]],\
                             [lon_forec_track_hycom_exp[::2][indx],lat_lon_matrix[0,i]],'km')[0][0]

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

    # Read HYCOM exp
    file = afiles[indx]

    #Reading time stamp
    year = int(file.split('/')[-1].split('.')[1][0:4])
    month = int(file.split('/')[-1].split('.')[1][4:6])
    day = int(file.split('/')[-1].split('.')[1][6:8])
    hour = int(file.split('/')[-1].split('.')[1][8:10])
    dt = int(file.split('/')[-1].split('.')[3][1:])
    timestamp_hycom = mdates.date2num(datetime(year,month,day,hour)) + dt/24
    #time_hycom.append(mdates.num2date(timestamp_hycom))

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
        R_hycom[i] = seawater.dist([lat_forec_track_hycom_exp[::2][indx],lat_lon_matrix_hycom[1,i]],\
                             [lon_forec_track_hycom_exp[::2][indx],lat_lon_matrix_hycom[0,i]],'km')[0][0]

    R_hycom_norm = R_hycom/Rmax
    okR_hycom_8 = R_hycom_norm <= 8

    RR_hycom_norm = R_hycom_norm[okR_hycom_8]

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

    Tmean_dens_crit_vec = np.ravel(Tmean_dens_crit)[okR_hycom_8]

    okRR_hycom = np.argsort(RR_hycom_norm)

    for n in np.arange(1,len(RR_norm_bins)):
        print(n)
        okR = np.logical_and(RR_hycom_norm >= RR_norm_bins[n-1],RR_hycom_norm <= RR_norm_bins[n])
        if len(np.where(okR)[0])==0:
            Tmean_max_hycom_exp[N,n] = np.nan
            Tmean_mean_hycom_exp[N,n] = np.nan
            Tmean_min_hycom_exp[N,n] = np.nan
        else:
            Tmean_max_hycom_exp[N,n] = np.nanmax(Tmean_dens_crit_vec[okR])
            Tmean_mean_hycom_exp[N,n] = np.nanmean(Tmean_dens_crit_vec[okR])
            Tmean_min_hycom_exp[N,n] = np.nanmin(Tmean_dens_crit_vec[okR])

    Tmean_dens_crit[Tmean_dens_crit<27] = 27
    kw = dict(levels=np.arange(27,29.6,0.2))
    ax = plt.subplot(grid[N, 4:6])
    plt.contourf(lon_hycom[0,oklon_hycom]-360,lat_hycom[oklat_hycom,0],Tmean_dens_crit,cmap=cmocean.cm.thermal,\
                 **kw,vmin=27,vmax=29.5)
    c=plt.colorbar()
    c.set_label('$^oC$',rotation=90, labelpad=10, fontsize=14)
    plt.contour(lon_hycom[0,oklon_hycom]-360,lat_hycom[oklat_hycom,0],Tmean_dens_crit,\
                colors='k',linestyle='-',**kw,alpha=0.3)
    plt.plot(eye_lon[0],eye_lat[0],'ok',markersize=7,markeredgecolor='lawngreen')
    plt.plot(lon_forec_track_hycom_exp,lat_forec_track_hycom_exp,'.-k')
    plt.plot(lat_lon_matrix[0][okR1],lat_lon_matrix[1][okR1],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR2],lat_lon_matrix[1][okR2],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR3],lat_lon_matrix[1][okR3],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR4],lat_lon_matrix[1][okR4],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR5],lat_lon_matrix[1][okR5],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR6],lat_lon_matrix[1][okR6],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR7],lat_lon_matrix[1][okR7],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR8],lat_lon_matrix[1][okR8],'.k',markersize=0.5)
    plt.plot(lon_maxwind,lat_maxwind,'*k',markersize=10,markeredgecolor='lawngreen')
    plt.axis('scaled')
    plt.xlim(np.min(lat_lon_matrix[0][okR_8])-0.1,np.max(lat_lon_matrix[0][okR_8])+0.1)
    plt.ylim(np.min(lat_lon_matrix[1][okR_8])-0.1,np.max(lat_lon_matrix[1][okR_8])+0.1)
    if N==0:
        plt.title('HWRF2020-HYCOM \n' + str(time_hwrf[0])[0:13],fontsize=14)
    else:
        plt.title(str(time_hwrf[0])[0:13],fontsize=14)

for N in np.arange(3):
    ax = plt.subplot(grid[N, 7:])
    plt.plot(RR_norm_bins_mid,Tmean_mean_pom_oper[N,1:],'X-',color='mediumorchid',label='HWRF2019-POM (IC clim.)',markeredgecolor='k',markersize=7)
    ax.fill_between(RR_norm_bins_mid,Tmean_min_pom_oper[N,1:],Tmean_max_pom_oper[N,1:],color='mediumorchid',alpha=0.1)

    plt.plot(RR_norm_bins_mid,Tmean_mean_pom_exp[N,1:],'^-',color='teal',label='HWRF2020-POM (RTOFS)',markeredgecolor='k',markersize=7)
    ax.fill_between(RR_norm_bins_mid,Tmean_min_pom_exp[N,1:],Tmean_max_pom_exp[N,1:],color='teal',alpha=0.1)

    plt.plot(RR_norm_bins_mid,Tmean_mean_hycom_exp[N,1:],'H-',color='darkorange',label='HWRF2020-HYCOM (RTOFS)',markeredgecolor='k',markersize=7)
    ax.fill_between(RR_norm_bins_mid,Tmean_min_hycom_exp[N,1:],Tmean_max_hycom_exp[N,1:],color='darkorange',alpha=0.1)

    plt.title(str(time[N])[2:15],fontsize=14)
    plt.ylim(27.0,29.5)
    plt.ylabel('$^oC$',fontsize=14)
    plt.xticks(np.arange(9))
    if N==1:
        plt.legend(loc='lower left',fontsize=10)
    if N==2:
        plt.xlabel('r/Rmax',fontsize=14)

file_name = folder_fig + 'MLT_around_eye'
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

#%%
###################### Start figure all panels OHC ##############################

fig, ax = plt.subplots(figsize=(15, 15))
grid = plt.GridSpec(3, 9, wspace=0.1, hspace=0.3,left=0.05,right=0.95)
fig.suptitle('Ocean Heat Content',fontsize=20)

RR_norm_bins = np.arange(0,8.1,0.5)
RR_norm_bins_mid = RR_norm_bins[0:-1] + np.diff(RR_norm_bins)/2

RR_norm_bins = np.arange(0,8.1,0.5)
RR_norm_bins_mid = RR_norm_bins[0:-1] + np.diff(RR_norm_bins)/2

OHC_max_pom_oper = np.empty((3,len(RR_norm_bins)))
OHC_max_pom_oper[:] = np.nan
OHC_min_pom_oper = np.empty((3,len(RR_norm_bins)))
OHC_min_pom_oper[:] = np.nan
OHC_mean_pom_oper = np.empty((3,len(RR_norm_bins)))
OHC_mean_pom_oper[:] = np.nan

time = []

for N,indx in enumerate(np.asarray([3,11,14])):
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
    Rmax = seawater.dist([lat_forec_track_pom_oper[::2][indx],lat_maxwind],\
                         [lon_forec_track_pom_oper[::2][indx],lon_maxwind],'km')[0][0]

    xlim = [lon_forec_track_pom_oper[::2][indx]-4,lon_forec_track_pom_oper[::2][indx]+4]
    ylim = [lat_forec_track_pom_oper[::2][indx]-4,lat_forec_track_pom_oper[::2][indx]+4]

    oklon = np.where(np.logical_and(lon_hwrf>xlim[0],lon_hwrf<xlim[1]))[0]
    oklat = np.where(np.logical_and(lat_hwrf>ylim[0],lat_hwrf<ylim[1]))[0]

    meshlon_lat = np.meshgrid(lon_hwrf[oklon],lat_hwrf[oklat])

    eye_lon = np.tile(lon_forec_track_pom_oper[::2][indx],meshlon_lat[0].shape[1])
    eye_lat = np.tile(lat_forec_track_pom_oper[::2][indx],meshlon_lat[0].shape[0])
    eye_matrix = np.meshgrid(eye_lon,eye_lat)

    lat_lon_matrix = np.stack((np.ravel(meshlon_lat[0]),np.ravel(meshlon_lat[1])),axis=1).T
    eye_lat_lon_matrix = np.stack((np.ravel(eye_matrix[0]),np.ravel(eye_matrix[1])),axis=1).T

    R = np.empty(lat_lon_matrix.shape[1])
    R[:] = np.nan
    for i in np.arange(lat_lon_matrix.shape[1]):
        R[i] = seawater.dist([lat_forec_track_pom_oper[::2][indx],lat_lon_matrix[1,i]],\
                             [lon_forec_track_pom_oper[::2][indx],lat_lon_matrix[0,i]],'km')[0][0]

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
    okRR = np.argsort(RR_norm)

    # POM oper
    print(POM_oper[indx])
    pom = xr.open_dataset(POM_oper[indx])
    tpom = np.asarray(pom['time'][:])
    temp_pom = np.asarray(pom['t'][0,:,:,:])
    salt_pom = np.asarray(pom['s'][0,:,:,:])
    dens_pom = np.asarray(pom['rho'][0,:,:,:]) * 1000 + 1000
    dens_pom[dens_pom==1000.0] = np.nan

    oklon_pom = np.where(np.logical_and(lon_pom_oper[0,:]>xlim[0],lon_pom_oper[0,:]<xlim[1]))[0]
    oklat_pom = np.where(np.logical_and(lat_pom_oper[:,0]>ylim[0],lat_pom_oper[:,0]<ylim[1]))[0]

    meshlon_lat_pom = np.meshgrid(lon_pom_oper[0,oklon_pom],lat_pom_oper[oklat_pom,0])

    lat_lon_matrix_pom = np.stack((np.ravel(meshlon_lat_pom[0]),np.ravel(meshlon_lat_pom[1])),axis=1).T

    R_pom = np.empty(lat_lon_matrix_pom.shape[1])
    R_pom[:] = np.nan
    for i in np.arange(lat_lon_matrix_pom.shape[1]):
        R_pom[i] = seawater.dist([lat_forec_track_pom_oper[::2][indx],lat_lon_matrix_pom[1,i]],\
                             [lon_forec_track_pom_oper[::2][indx],lat_lon_matrix_pom[0,i]],'km')[0][0]

    R_pom_norm = R_pom/Rmax
    okR_pom_8 = R_pom_norm <= 8
    okR_pom_larger_8 = R_pom_norm > 8

    RR_pom_norm = R_pom_norm[okR_pom_8]

    dt = 0.2
    drho = 0.125

    timestamp_POM_oper = mdates.date2num(tpom)
    target_temp_POM_oper = temp_pom[:,oklat_pom,:][:,:,oklon_pom]
    target_salt_POM_oper = salt_pom[:,oklat_pom,:][:,:,oklon_pom]
    target_dens_POM_oper = dens_pom[:,oklat_pom,:][:,:,oklon_pom]
    target_depth_POM_oper = zmatrix_pom_oper[:,oklat_pom,:][:,:,oklon_pom]

    okRR_pom = np.argsort(RR_pom_norm)

    # for POM operational
    OHC = OHC_from_3D_fields(target_temp_POM_oper,target_depth_POM_oper,\
                             target_dens_POM_oper)*10**(-7)

    OHC_vec = np.ravel(OHC)[okR_pom_8]
    okRR_pom = np.argsort(RR_pom_norm)

    for n in np.arange(1,len(RR_norm_bins)):
        print(n)
        okR = np.logical_and(RR_pom_norm >= RR_norm_bins[n-1],RR_pom_norm <= RR_norm_bins[n])
        if len(np.where(okR)[0])==0:
            OHC_max_pom_oper[N,n] = np.nan
            OHC_mean_pom_oper[N,n] = np.nan
            OHC_min_pom_oper[N,n] = np.nan
        else:
            OHC_max_pom_oper[N,n] = np.nanmax(OHC_vec[okR])
            OHC_mean_pom_oper[N,n] = np.nanmean(OHC_vec[okR])
            OHC_min_pom_oper[N,n] = np.nanmin(OHC_vec[okR])

    kw = dict(levels=np.arange(0,101,10))

    ax = plt.subplot(grid[N, 0:2])
    plt.contourf(lon_pom_oper[0,oklon_pom],lat_pom_oper[oklat_pom,0],OHC,cmap='magma',\
                 **kw,vmin=0,vmax=100)
    #c=plt.colorbar()
    #c.set_label('($KJ/cm^2$)',rotation=90, labelpad=15, fontsize=14)
    plt.contour(lon_pom_oper[0,oklon_pom],lat_pom_oper[oklat_pom,0],OHC,\
                colors='k',linestyle='-',**kw,alpha=0.3)
    plt.plot(eye_lon[0],eye_lat[0],'ok',markersize=7,markeredgecolor='lawngreen')
    plt.plot(lon_forec_track_pom_oper,lat_forec_track_pom_oper,'.-k')
    plt.plot(lat_lon_matrix[0][okR1],lat_lon_matrix[1][okR1],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR2],lat_lon_matrix[1][okR2],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR3],lat_lon_matrix[1][okR3],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR4],lat_lon_matrix[1][okR4],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR5],lat_lon_matrix[1][okR5],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR6],lat_lon_matrix[1][okR6],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR7],lat_lon_matrix[1][okR7],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR8],lat_lon_matrix[1][okR8],'.k',markersize=0.5)
    plt.plot(lon_maxwind,lat_maxwind,'*k',markersize=10,markeredgecolor='lawngreen')
    plt.axis('scaled')
    plt.xlim(np.min(lat_lon_matrix[0][okR_8])-0.1,np.max(lat_lon_matrix[0][okR_8])+0.1)
    plt.ylim(np.min(lat_lon_matrix[1][okR_8])-0.1,np.max(lat_lon_matrix[1][okR_8])+0.1)
    if N==0:
        plt.title('HWRF2019-POM \n'+ str(time_hwrf[0])[0:13],fontsize=14)
    else:
        plt.title(str(time_hwrf[0])[0:13],fontsize=14)
    if N==0:
        plt.text(-68.5, 17.5, 'T1',fontsize=16)
    if N==1:
        plt.text(-75.7, 25.5, 'T2',fontsize=16)
    if N==2:
        plt.text(-78.2, 27.5, 'T3',fontsize=16)

#%% MLT and OHC under Dorian POM experimental
OHC_max_pom_exp = np.empty((3,len(RR_norm_bins)))
OHC_max_pom_exp[:] = np.nan
OHC_min_pom_exp = np.empty((3,len(RR_norm_bins)))
OHC_min_pom_exp[:] = np.nan
OHC_mean_pom_exp = np.empty((3,len(RR_norm_bins)))
OHC_mean_pom_exp[:] = np.nan

time = []

for N,indx in enumerate(np.asarray([3,11,14])):
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

    time.append(np.asarray(HWRF.variables['time'][:]))

    wind_int = np.sqrt(UGRD_hwrf**2 + VGRD_hwrf**2)
    max_wind= np.max(wind_int)
    okwind = np.where(wind_int == max_wind)
    lat_maxwind = lat_hwrf[okwind[0][0]]
    lon_maxwind = lon_hwrf[okwind[1][0]]
    Rmax = seawater.dist([lat_forec_track_pom_exp[::2][indx],lat_maxwind],\
                         [lon_forec_track_pom_exp[::2][indx],lon_maxwind],'km')[0][0]

    xlim = [lon_forec_track_pom_exp[::2][indx]-4,lon_forec_track_pom_exp[::2][indx]+4]
    ylim = [lat_forec_track_pom_exp[::2][indx]-4,lat_forec_track_pom_exp[::2][indx]+4]

    oklon = np.where(np.logical_and(lon_hwrf>xlim[0],lon_hwrf<xlim[1]))[0]
    oklat = np.where(np.logical_and(lat_hwrf>ylim[0],lat_hwrf<ylim[1]))[0]

    meshlon_lat = np.meshgrid(lon_hwrf[oklon],lat_hwrf[oklat])

    eye_lon = np.tile(lon_forec_track_pom_exp[::2][indx],meshlon_lat[0].shape[1])
    eye_lat = np.tile(lat_forec_track_pom_exp[::2][indx],meshlon_lat[0].shape[0])
    eye_matrix = np.meshgrid(eye_lon,eye_lat)

    lat_lon_matrix = np.stack((np.ravel(meshlon_lat[0]),np.ravel(meshlon_lat[1])),axis=1).T
    eye_lat_lon_matrix = np.stack((np.ravel(eye_matrix[0]),np.ravel(eye_matrix[1])),axis=1).T

    R = np.empty(lat_lon_matrix.shape[1])
    R[:] = np.nan
    for i in np.arange(lat_lon_matrix.shape[1]):
        R[i] = seawater.dist([lat_forec_track_pom_exp[::2][indx],lat_lon_matrix[1,i]],\
                             [lon_forec_track_pom_exp[::2][indx],lat_lon_matrix[0,i]],'km')[0][0]

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
    okRR = np.argsort(RR_norm)

    # POM exp
    print(POM_exp[indx])
    pom = xr.open_dataset(POM_exp[indx])
    tpom = np.asarray(pom['time'][:])
    temp_pom = np.asarray(pom['t'][0,:,:,:])
    salt_pom = np.asarray(pom['s'][0,:,:,:])
    dens_pom = np.asarray(pom['rho'][0,:,:,:]) * 1000 + 1000
    dens_pom[dens_pom==1000.0] = np.nan

    oklon_pom = np.where(np.logical_and(lon_pom_exp[0,:]>xlim[0],lon_pom_exp[0,:]<xlim[1]))[0]
    oklat_pom = np.where(np.logical_and(lat_pom_exp[:,0]>ylim[0],lat_pom_exp[:,0]<ylim[1]))[0]

    meshlon_lat_pom = np.meshgrid(lon_pom_exp[0,oklon_pom],lat_pom_exp[oklat_pom,0])

    lat_lon_matrix_pom = np.stack((np.ravel(meshlon_lat_pom[0]),np.ravel(meshlon_lat_pom[1])),axis=1).T

    R_pom = np.empty(lat_lon_matrix_pom.shape[1])
    R_pom[:] = np.nan
    for i in np.arange(lat_lon_matrix_pom.shape[1]):
        R_pom[i] = seawater.dist([lat_forec_track_pom_exp[::2][indx],lat_lon_matrix_pom[1,i]],\
                             [lon_forec_track_pom_exp[::2][indx],lat_lon_matrix_pom[0,i]],'km')[0][0]

    R_pom_norm = R_pom/Rmax
    okR_pom_8 = R_pom_norm <= 8

    RR_pom_norm = R_pom_norm[okR_pom_8]

    dt = 0.2
    drho = 0.125

    timestamp_POM_exp = mdates.date2num(tpom)
    target_temp_POM_exp = temp_pom[:,oklat_pom,:][:,:,oklon_pom]
    target_salt_POM_exp = salt_pom[:,oklat_pom,:][:,:,oklon_pom]
    target_dens_POM_exp = dens_pom[:,oklat_pom,:][:,:,oklon_pom]
    target_depth_POM_exp = zmatrix_pom_exp[:,oklat_pom,:][:,:,oklon_pom]

    # for POM operational
    OHC = OHC_from_3D_fields(target_temp_POM_exp,target_depth_POM_exp,\
                             target_dens_POM_exp)*10**(-7)

    OHC_vec = np.ravel(OHC)[okR_pom_8]

    okRR_pom = np.argsort(RR_pom_norm)

    for n in np.arange(1,len(RR_norm_bins)):
        print(n)
        okR = np.logical_and(RR_pom_norm >= RR_norm_bins[n-1],RR_pom_norm <= RR_norm_bins[n])
        if len(np.where(okR)[0])==0:
            OHC_max_pom_exp[N,n] = np.nan
            OHC_mean_pom_exp[N,n] = np.nan
            OHC_min_pom_exp[N,n] = np.nan
        else:
            OHC_max_pom_exp[N,n] = np.nanmax(OHC_vec[okR])
            OHC_mean_pom_exp[N,n] = np.nanmean(OHC_vec[okR])
            OHC_min_pom_exp[N,n] = np.nanmin(OHC_vec[okR])

    kw = dict(levels=np.arange(0,101,10))
    ax = plt.subplot(grid[N, 2:4])
    plt.contourf(lon_pom_exp[0,oklon_pom],lat_pom_exp[oklat_pom,0],OHC,cmap='magma',\
                 **kw,vmin=0,vmax=100)
    #=plt.colorbar()
    #c.set_label('($KJ/cm^2$)',rotation=90, labelpad=15, fontsize=14)
    plt.contour(lon_pom_exp[0,oklon_pom],lat_pom_exp[oklat_pom,0],OHC,\
                colors='k',linestyle='-',**kw,alpha=0.3)
    plt.plot(eye_lon[0],eye_lat[0],'ok',markersize=7,markeredgecolor='lawngreen')
    plt.plot(lon_forec_track_pom_exp,lat_forec_track_pom_exp,'.-k')
    plt.plot(lat_lon_matrix[0][okR1],lat_lon_matrix[1][okR1],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR2],lat_lon_matrix[1][okR2],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR3],lat_lon_matrix[1][okR3],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR4],lat_lon_matrix[1][okR4],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR5],lat_lon_matrix[1][okR5],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR6],lat_lon_matrix[1][okR6],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR7],lat_lon_matrix[1][okR7],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR8],lat_lon_matrix[1][okR8],'.k',markersize=0.5)
    plt.plot(lon_maxwind,lat_maxwind,'*k',markersize=10,markeredgecolor='lawngreen')
    plt.axis('scaled')
    plt.xlim(np.min(lat_lon_matrix[0][okR_8])-0.1,np.max(lat_lon_matrix[0][okR_8])+0.1)
    plt.ylim(np.min(lat_lon_matrix[1][okR_8])-0.1,np.max(lat_lon_matrix[1][okR_8])+0.1)
    if N==0:
        plt.title('HWRF2020-POM \n' + str(time_hwrf[0])[0:13],fontsize=14)
    else:
        plt.title(str(time_hwrf[0])[0:13],fontsize=14)

#%% MLT and OHC under Dorian HYCOM experimental
OHC_max_hycom_exp = np.empty((3,len(RR_norm_bins)))
OHC_max_hycom_exp[:] = np.nan
OHC_min_hycom_exp = np.empty((3,len(RR_norm_bins)))
OHC_min_hycom_exp[:] = np.nan
OHC_mean_hycom_exp = np.empty((3,len(RR_norm_bins)))
OHC_mean_hycom_exp[:] = np.nan

time = []

afiles = sorted(glob.glob(os.path.join(folder_hycom_exp,prefix_hycom+'*.a')))

for N,indx in enumerate(np.asarray([3,11,14])):
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

    time.append(np.asarray(HWRF.variables['time'][:]))

    wind_int = np.sqrt(UGRD_hwrf**2 + VGRD_hwrf**2)
    max_wind= np.max(wind_int)
    okwind = np.where(wind_int == max_wind)
    lat_maxwind = lat_hwrf[okwind[0][0]]
    lon_maxwind = lon_hwrf[okwind[1][0]]
    Rmax = seawater.dist([lat_forec_track_hycom_exp[::2][indx],lat_maxwind],\
                         [lon_forec_track_hycom_exp[::2][indx],lon_maxwind],'km')[0][0]

    xlim = [lon_forec_track_hycom_exp[::2][indx]-4,lon_forec_track_hycom_exp[::2][indx]+4]
    ylim = [lat_forec_track_hycom_exp[::2][indx]-4,lat_forec_track_hycom_exp[::2][indx]+4]

    oklon = np.where(np.logical_and(lon_hwrf>xlim[0],lon_hwrf<xlim[1]))[0]
    oklat = np.where(np.logical_and(lat_hwrf>ylim[0],lat_hwrf<ylim[1]))[0]

    meshlon_lat = np.meshgrid(lon_hwrf[oklon],lat_hwrf[oklat])

    eye_lon = np.tile(lon_forec_track_hycom_exp[::2][indx],meshlon_lat[0].shape[1])
    eye_lat = np.tile(lat_forec_track_hycom_exp[::2][indx],meshlon_lat[0].shape[0])
    eye_matrix = np.meshgrid(eye_lon,eye_lat)

    lat_lon_matrix = np.stack((np.ravel(meshlon_lat[0]),np.ravel(meshlon_lat[1])),axis=1).T
    eye_lat_lon_matrix = np.stack((np.ravel(eye_matrix[0]),np.ravel(eye_matrix[1])),axis=1).T

    R = np.empty(lat_lon_matrix.shape[1])
    R[:] = np.nan
    for i in np.arange(lat_lon_matrix.shape[1]):
        R[i] = seawater.dist([lat_forec_track_hycom_exp[::2][indx],lat_lon_matrix[1,i]],\
                             [lon_forec_track_hycom_exp[::2][indx],lat_lon_matrix[0,i]],'km')[0][0]

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
    okRR = np.argsort(RR_norm)

    # Read HYCOM exp
    file = afiles[indx]

    #Reading time stamp
    year = int(file.split('/')[-1].split('.')[1][0:4])
    month = int(file.split('/')[-1].split('.')[1][4:6])
    day = int(file.split('/')[-1].split('.')[1][6:8])
    hour = int(file.split('/')[-1].split('.')[1][8:10])
    dt = int(file.split('/')[-1].split('.')[3][1:])
    timestamp_hycom = mdates.date2num(datetime(year,month,day,hour)) + dt/24
    #time_hycom.append(mdates.num2date(timestamp_hycom))

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
        R_hycom[i] = seawater.dist([lat_forec_track_hycom_exp[::2][indx],lat_lon_matrix_hycom[1,i]],\
                             [lon_forec_track_hycom_exp[::2][indx],lat_lon_matrix_hycom[0,i]],'km')[0][0]

    R_hycom_norm = R_hycom/Rmax
    okR_hycom_8 = R_hycom_norm <= 8

    RR_hycom_norm = R_hycom_norm[okR_hycom_8]

    dt = 0.2
    drho = 0.125

    target_temp_hycom_exp = np.asarray(temp_hycom[oklat_hycom,:,:][:,oklon_hycom,:])
    target_salt_hycom_exp = np.asarray(salt_hycom[oklat_hycom,:,:][:,oklon_hycom,:])
    target_depth_hycom_exp = np.asarray(depth_HYCOM_exp)

    target_dens_hycom_exp = seawater.dens(target_salt_hycom_exp,target_temp_hycom_exp,target_depth_hycom_exp)

    okRR_hycom = np.argsort(RR_hycom_norm)

    # for HYCOM exp
    OHC = OHC_from_3D_fields(target_temp_hycom_exp.T,target_depth_hycom_exp,\
                             target_dens_hycom_exp.T)*10**(-7)
    OHC = OHC.T

    OHC_vec = np.ravel(OHC)[okR_hycom_8]

    okRR_hycom = np.argsort(RR_hycom_norm)

    for n in np.arange(1,len(RR_norm_bins)):
        print(n)
        okR = np.logical_and(RR_hycom_norm >= RR_norm_bins[n-1],RR_hycom_norm <= RR_norm_bins[n])
        if len(np.where(okR)[0])==0:
            OHC_max_hycom_exp[N,n] = np.nan
            OHC_mean_hycom_exp[N,n] = np.nan
            OHC_min_hycom_exp[N,n] = np.nan
        else:
            OHC_max_hycom_exp[N,n] = np.nanmax(OHC_vec[okR])
            OHC_mean_hycom_exp[N,n] = np.nanmean(OHC_vec[okR])
            OHC_min_hycom_exp[N,n] = np.nanmin(OHC_vec[okR])

    kw = dict(levels=np.arange(0,101,10))
    ax = plt.subplot(grid[N, 4:6])
    plt.contourf(lon_hycom[0,oklon_hycom]-360,lat_hycom[oklat_hycom,0],OHC,cmap='magma',\
                 **kw,vmin=0,vmax=100)
    c=plt.colorbar()
    c.set_label('$KJ/cm^2$',rotation=90, labelpad=10, fontsize=14)
    plt.contour(lon_hycom[0,oklon_hycom]-360,lat_hycom[oklat_hycom,0],OHC,\
                colors='k',linestyle='-',**kw,alpha=0.3)
    plt.plot(eye_lon[0],eye_lat[0],'ok',markersize=7,markeredgecolor='lawngreen')
    plt.plot(lon_forec_track_hycom_exp,lat_forec_track_hycom_exp,'.-k')
    plt.plot(lat_lon_matrix[0][okR1],lat_lon_matrix[1][okR1],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR2],lat_lon_matrix[1][okR2],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR3],lat_lon_matrix[1][okR3],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR4],lat_lon_matrix[1][okR4],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR5],lat_lon_matrix[1][okR5],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR6],lat_lon_matrix[1][okR6],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR7],lat_lon_matrix[1][okR7],'.k',markersize=0.5)
    plt.plot(lat_lon_matrix[0][okR8],lat_lon_matrix[1][okR8],'.k',markersize=0.5)
    plt.plot(lon_maxwind,lat_maxwind,'*k',markersize=10,markeredgecolor='lawngreen')
    plt.axis('scaled')
    plt.xlim(np.min(lat_lon_matrix[0][okR_8])-0.1,np.max(lat_lon_matrix[0][okR_8])+0.1)
    plt.ylim(np.min(lat_lon_matrix[1][okR_8])-0.1,np.max(lat_lon_matrix[1][okR_8])+0.1)
    if N==0:
        plt.title('HWRF2020-HYCOM \n' + str(time_hwrf[0])[0:13],fontsize=14)
    else:
        plt.title(str(time_hwrf[0])[0:13],fontsize=14)

for N in np.arange(3):
    ax = plt.subplot(grid[N, 7:])
    plt.plot(RR_norm_bins_mid,OHC_mean_pom_oper[N,1:],'X-',color='mediumorchid',label='HWRF2019-POM (IC clim.)',markeredgecolor='k',markersize=7)
    ax.fill_between(RR_norm_bins_mid,OHC_min_pom_oper[N,1:],OHC_max_pom_oper[N,1:],color='mediumorchid',alpha=0.1)

    plt.plot(RR_norm_bins_mid,OHC_mean_pom_exp[N,1:],'^-',color='teal',label='HWRF2020-POM (RTOFS)',markeredgecolor='k',markersize=7)
    ax.fill_between(RR_norm_bins_mid,OHC_min_pom_exp[N,1:],OHC_max_pom_exp[N,1:],color='teal',alpha=0.1)

    plt.plot(RR_norm_bins_mid,OHC_mean_hycom_exp[N,1:],'H-',color='darkorange',label='HWRF2020-HYCOM (RTOFS)',markeredgecolor='k',markersize=7)
    ax.fill_between(RR_norm_bins_mid,OHC_min_hycom_exp[N,1:],OHC_max_hycom_exp[N,1:],color='darkorange',alpha=0.1)

    plt.title(str(time[N])[2:15],fontsize=14)
    plt.ylim(0,100)
    plt.ylabel('$Kj/cm^2$',fontsize=14)
    plt.xticks(np.arange(9))
    if N==1:
        plt.legend(loc='upper left',fontsize=10)
    if N==2:
        plt.xlabel('r/Rmax',fontsize=14)

file_name = folder_fig + 'OHC_around_eye'
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

#%%
################################################################################
#%% MLT and OHC under Dorian POM operational

#RR_norm_bins = np.arange(0,8,0.1)
RR_norm_bins = np.arange(0,8.1,0.5)
RR_norm_bins_mid = RR_norm_bins[0:-1] + np.diff(RR_norm_bins)/2

Tmean_max_pom_oper = np.empty((3,len(RR_norm_bins)))
Tmean_max_pom_oper[:] = np.nan
Tmean_min_pom_oper = np.empty((3,len(RR_norm_bins)))
Tmean_min_pom_oper[:] = np.nan
Tmean_mean_pom_oper = np.empty((3,len(RR_norm_bins)))
Tmean_mean_pom_oper[:] = np.nan

OHC_max_pom_oper = np.empty((3,len(RR_norm_bins)))
OHC_max_pom_oper[:] = np.nan
OHC_min_pom_oper = np.empty((3,len(RR_norm_bins)))
OHC_min_pom_oper[:] = np.nan
OHC_mean_pom_oper = np.empty((3,len(RR_norm_bins)))
OHC_mean_pom_oper[:] = np.nan

time = []

for N,indx in enumerate(np.asarray([3,11,14])):
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
    Rmax = seawater.dist([lat_forec_track_pom_oper[::2][indx],lat_maxwind],\
                         [lon_forec_track_pom_oper[::2][indx],lon_maxwind],'km')[0][0]

    xlim = [lon_forec_track_pom_oper[::2][indx]-4,lon_forec_track_pom_oper[::2][indx]+4]
    ylim = [lat_forec_track_pom_oper[::2][indx]-4,lat_forec_track_pom_oper[::2][indx]+4]

    oklon = np.where(np.logical_and(lon_hwrf>xlim[0],lon_hwrf<xlim[1]))[0]
    oklat = np.where(np.logical_and(lat_hwrf>ylim[0],lat_hwrf<ylim[1]))[0]

    meshlon_lat = np.meshgrid(lon_hwrf[oklon],lat_hwrf[oklat])
    #meshlat = np.meshgrid(lat_hwrf[oklat],lon_hwrf[oklon])

    eye_lon = np.tile(lon_forec_track_pom_oper[::2][indx],meshlon_lat[0].shape[1])
    eye_lat = np.tile(lat_forec_track_pom_oper[::2][indx],meshlon_lat[0].shape[0])
    eye_matrix = np.meshgrid(eye_lon,eye_lat)

    lat_lon_matrix = np.stack((np.ravel(meshlon_lat[0]),np.ravel(meshlon_lat[1])),axis=1).T
    eye_lat_lon_matrix = np.stack((np.ravel(eye_matrix[0]),np.ravel(eye_matrix[1])),axis=1).T

    R = np.empty(lat_lon_matrix.shape[1])
    R[:] = np.nan
    for i in np.arange(lat_lon_matrix.shape[1]):
        R[i] = seawater.dist([lat_forec_track_pom_oper[::2][indx],lat_lon_matrix[1,i]],\
                             [lon_forec_track_pom_oper[::2][indx],lat_lon_matrix[0,i]],'km')[0][0]

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

    okRR = np.argsort(RR_norm)
    sshla_vec_interp = np.interp(RR_norm_bins,RR_norm[okRR],sshla_vec[okRR])

    # POM oper
    print(POM_oper[indx])
    pom = xr.open_dataset(POM_oper[indx])
    tpom = np.asarray(pom['time'][:])
    temp_pom = np.asarray(pom['t'][0,:,:,:])
    salt_pom = np.asarray(pom['s'][0,:,:,:])
    dens_pom = np.asarray(pom['rho'][0,:,:,:]) * 1000 + 1000
    dens_pom[dens_pom==1000.0] = np.nan

    oklon_pom = np.where(np.logical_and(lon_pom_oper[0,:]>xlim[0],lon_pom_oper[0,:]<xlim[1]))[0]
    oklat_pom = np.where(np.logical_and(lat_pom_oper[:,0]>ylim[0],lat_pom_oper[:,0]<ylim[1]))[0]

    meshlon_lat_pom = np.meshgrid(lon_pom_oper[0,oklon_pom],lat_pom_oper[oklat_pom,0])

    lat_lon_matrix_pom = np.stack((np.ravel(meshlon_lat_pom[0]),np.ravel(meshlon_lat_pom[1])),axis=1).T

    R_pom = np.empty(lat_lon_matrix_pom.shape[1])
    R_pom[:] = np.nan
    for i in np.arange(lat_lon_matrix_pom.shape[1]):
        R_pom[i] = seawater.dist([lat_forec_track_pom_oper[::2][indx],lat_lon_matrix_pom[1,i]],\
                             [lon_forec_track_pom_oper[::2][indx],lat_lon_matrix_pom[0,i]],'km')[0][0]

    R_pom_norm = R_pom/Rmax
    okR_pom_8 = R_pom_norm <= 8
    okR_pom_larger_8 = R_pom_norm > 8

    RR_pom_norm = R_pom_norm[okR_pom_8]

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

    Tmean_dens_crit_vec = np.ravel(Tmean_dens_crit)[okR_pom_8]

    okRR_pom = np.argsort(RR_pom_norm)
    '''
    c = Tmean_dens_crit_vec
    fig,ax = plt.subplots()
    plt.scatter(lat_lon_matrix_pom[0][okR_pom_8],lat_lon_matrix_pom[1][okR_pom_8],\
                c=c,s=100,cmap=cmocean.cm.thermal,vmin=27,vmax=29.5)
    c=plt.colorbar()
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
    plt.xlim(np.min(lat_lon_matrix[0][okR_8])-0.1,np.max(lat_lon_matrix[0][okR_8])+0.1)
    plt.ylim(np.min(lat_lon_matrix[1][okR_8])-0.1,np.max(lat_lon_matrix[1][okR_8])+0.1)
    c.set_label('$^oC$',rotation=90, labelpad=15, fontsize=14)
    plt.title('Mixed Layer Temperature ' + str(time_hwrf[0])[0:13] + '\n HWRF2019-POM (IC Clim.)',fontsize=15)
    '''
    kw = dict(levels=np.arange(27,29.6,0.2))
    fig,ax = plt.subplots()
    plt.contourf(lon_pom_oper[0,oklon_pom],lat_pom_oper[oklat_pom,0],Tmean_dens_crit,cmap=cmocean.cm.thermal,\
                 **kw,vmin=27,vmax=29.5)
    c=plt.colorbar()
    plt.contour(lon_pom_oper[0,oklon_pom],lat_pom_oper[oklat_pom,0],Tmean_dens_crit,\
                colors='k',linestyle='-',**kw,alpha=0.3)
    plt.plot(eye_lon[0],eye_lat[0],'ok',markersize=7,markeredgecolor='lawngreen')
    plt.plot(lon_forec_track_pom_oper,lat_forec_track_pom_oper,'.-k')
    plt.plot(lat_lon_matrix[0][okR1],lat_lon_matrix[1][okR1],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR2],lat_lon_matrix[1][okR2],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR3],lat_lon_matrix[1][okR3],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR4],lat_lon_matrix[1][okR4],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR5],lat_lon_matrix[1][okR5],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR6],lat_lon_matrix[1][okR6],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR7],lat_lon_matrix[1][okR7],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR8],lat_lon_matrix[1][okR8],'.k',markersize=1)
    plt.plot(lon_maxwind,lat_maxwind,'*k',markersize=10,markeredgecolor='lawngreen')
    plt.xlim(np.min(lat_lon_matrix[0][okR_8])-0.1,np.max(lat_lon_matrix[0][okR_8])+0.1)
    plt.ylim(np.min(lat_lon_matrix[1][okR_8])-0.1,np.max(lat_lon_matrix[1][okR_8])+0.1)
    c.set_label('$^oC$',rotation=90, labelpad=15, fontsize=14)
    plt.title('HWRF2019-POM (IC Clim.) \n'+ str(time_hwrf[0])[0:13],fontsize=15)
    plt.title('Mixed Layer Temperature ' + str(time_hwrf[0])[0:13] + '\n HWRF2019-POM (IC Clim.)',fontsize=15)
    file_name = folder_fig + 'HWRF2019_POM_Dorian_Temp_MLD_cycle_' + cycle + '_' + str(tpom[0])[0:13]
    plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1)

    '''
    Tmean_dens_crit_vec2 = np.ravel(Tmean_dens_crit)
    Tmean_dens_crit_vec2[okR_pom_larger_8] = np.nan
    Tmean_dens_crit2 = np.reshape(Tmean_dens_crit_vec2,(Tmean_dens_crit.shape[0],Tmean_dens_crit.shape[1]))
    kw = dict(levels=np.arange(27,29.6,0.2))
    fig,ax = plt.subplots()
    plt.contourf(lon_pom_oper[0,oklon_pom],lat_pom_oper[oklat_pom,0],Tmean_dens_crit2,cmap=cmocean.cm.thermal,\
                 **kw,vmin=27,vmax=29.5)
    c=plt.colorbar()
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
    plt.xlim(np.min(lat_lon_matrix[0][okR_8])-0.1,np.max(lat_lon_matrix[0][okR_8])+0.1)
    plt.ylim(np.min(lat_lon_matrix[1][okR_8])-0.1,np.max(lat_lon_matrix[1][okR_8])+0.1)
    c.set_label('$^oC$',rotation=90, labelpad=15, fontsize=14)
    plt.title('Mixed Layer Temperature ' + str(time_hwrf[0])[0:13] + '\n HWRF2019-POM (IC Clim.)',fontsize=15)
    '''
    # for POM operational
    OHC = OHC_from_3D_fields(target_temp_POM_oper,target_depth_POM_oper,\
                             target_dens_POM_oper)*10**(-7)

    OHC_vec = np.ravel(OHC)[okR_pom_8]

    okRR_pom = np.argsort(RR_pom_norm)
    '''
    c = OHC_vec
    fig,ax = plt.subplots()
    plt.scatter(lat_lon_matrix_pom[0][okR_pom_8],lat_lon_matrix_pom[1][okR_pom_8],c=c,cmap='magma',vmin=20,vmax=90)
    c=plt.colorbar()
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
    plt.xlim(np.min(lat_lon_matrix[0][okR_8])-0.1,np.max(lat_lon_matrix[0][okR_8])+0.1)
    plt.ylim(np.min(lat_lon_matrix[1][okR_8])-0.1,np.max(lat_lon_matrix[1][okR_8])+0.1)
    c.set_label('($KJ/cm^2$)',rotation=90, labelpad=15, fontsize=14)
    plt.title('Ocean Heat Content ' + str(time_hwrf[0])[0:13] + '\n HWRF2019-POM (IC Clim.)',fontsize=15)
    '''

    kw = dict(levels=np.arange(0,101,10))
    fig,ax = plt.subplots()
    plt.contourf(lon_pom_oper[0,oklon_pom],lat_pom_oper[oklat_pom,0],OHC,cmap='magma',\
                 **kw,vmin=0,vmax=100)
    c=plt.colorbar()
    plt.contour(lon_pom_oper[0,oklon_pom],lat_pom_oper[oklat_pom,0],OHC,\
                colors='k',linestyle='-',**kw,alpha=0.3)
    c.set_label('$KJ/cm^2$',rotation=90, labelpad=15, fontsize=14)
    plt.plot(eye_lon[0],eye_lat[0],'ok',markersize=7,markeredgecolor='lawngreen')
    plt.plot(lon_forec_track_pom_oper,lat_forec_track_pom_oper,'.-k')
    plt.plot(lat_lon_matrix[0][okR1],lat_lon_matrix[1][okR1],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR2],lat_lon_matrix[1][okR2],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR3],lat_lon_matrix[1][okR3],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR4],lat_lon_matrix[1][okR4],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR5],lat_lon_matrix[1][okR5],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR6],lat_lon_matrix[1][okR6],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR7],lat_lon_matrix[1][okR7],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR8],lat_lon_matrix[1][okR8],'.k',markersize=1)
    plt.plot(lon_maxwind,lat_maxwind,'*k',markersize=10,markeredgecolor='lawngreen')
    plt.xlim(np.min(lat_lon_matrix[0][okR_8])-0.1,np.max(lat_lon_matrix[0][okR_8])+0.1)
    plt.ylim(np.min(lat_lon_matrix[1][okR_8])-0.1,np.max(lat_lon_matrix[1][okR_8])+0.1)
    c.set_label('($KJ/cm^2$)',rotation=90, labelpad=15, fontsize=14)
    plt.title('Ocean Heat Content ' + str(time_hwrf[0])[0:13] + '\n HWRF2019-POM (IC Clim.)',fontsize=15)
    file_name = folder_fig + 'HWRF2019_POM_Dorian_OHC_cycle_' + cycle + '_' + str(tpom[0])[0:13]
    plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1)

    for n in np.arange(1,len(RR_norm_bins)):
        print(n)
        okR = np.logical_and(RR_pom_norm >= RR_norm_bins[n-1],RR_pom_norm <= RR_norm_bins[n])
        if len(np.where(okR)[0])==0:
            Tmean_max_pom_oper[N,n] = np.nan
            Tmean_mean_pom_oper[N,n] = np.nan
            Tmean_min_pom_oper[N,n] = np.nan
        else:
            Tmean_max_pom_oper[N,n] = np.nanmax(Tmean_dens_crit_vec[okR])
            Tmean_mean_pom_oper[N,n] = np.nanmean(Tmean_dens_crit_vec[okR])
            Tmean_min_pom_oper[N,n] = np.nanmin(Tmean_dens_crit_vec[okR])

    for n in np.arange(1,len(RR_norm_bins)):
        print(n)
        okR = np.logical_and(RR_pom_norm >= RR_norm_bins[n-1],RR_pom_norm <= RR_norm_bins[n])
        if len(np.where(okR)[0])==0:
            OHC_max_pom_oper[N,n] = np.nan
            OHC_mean_pom_oper[N,n] = np.nan
            OHC_min_pom_oper[N,n] = np.nan
        else:
            OHC_max_pom_oper[N,n] = np.nanmax(OHC_vec[okR])
            OHC_mean_pom_oper[N,n] = np.nanmean(OHC_vec[okR])
            OHC_min_pom_oper[N,n] = np.nanmin(OHC_vec[okR])
        '''
        maxv = np.nanmean(sshla_vec[okR]) + 2*np.nanstd(sshla_vec[okR])
        minv = np.nanmean(sshla_vec[okR]) - 2*np.nanstd(sshla_vec[okR])
        okk = np.logical_and(sshla_vec[okR] >= minv,sshla_vec[okR] <= maxv)
        sshla_vec_max_pom_exp[N,n] = np.nanmax(sshla_vec[okR][okk])
        sshla_vec_mean_pom_exp[N,n] = np.nanmean(sshla_vec[okR][okk])
        sshla_vec_min_pom_exp[N,n] = np.nanmin(sshla_vec[okR][okk])
        '''

#%% MLT and OHC under Dorian POM experimental

#RR_norm_bins = np.arange(0,8,0.1)
#RR_norm_bins_mid = RR_norm_bins[0:-1] + np.diff(RR_norm_bins)/2

Tmean_max_pom_exp = np.empty((3,len(RR_norm_bins)))
Tmean_max_pom_exp[:] = np.nan
Tmean_min_pom_exp = np.empty((3,len(RR_norm_bins)))
Tmean_min_pom_exp[:] = np.nan
Tmean_mean_pom_exp = np.empty((3,len(RR_norm_bins)))
Tmean_mean_pom_exp[:] = np.nan

OHC_max_pom_exp = np.empty((3,len(RR_norm_bins)))
OHC_max_pom_exp[:] = np.nan
OHC_min_pom_exp = np.empty((3,len(RR_norm_bins)))
OHC_min_pom_exp[:] = np.nan
OHC_mean_pom_exp = np.empty((3,len(RR_norm_bins)))
OHC_mean_pom_exp[:] = np.nan

time = []

for N,indx in enumerate(np.asarray([3,11,14])):
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

    time.append(np.asarray(HWRF.variables['time'][:]))

    wind_int = np.sqrt(UGRD_hwrf**2 + VGRD_hwrf**2)
    max_wind= np.max(wind_int)
    okwind = np.where(wind_int == max_wind)
    lat_maxwind = lat_hwrf[okwind[0][0]]
    lon_maxwind = lon_hwrf[okwind[1][0]]
    Rmax = seawater.dist([lat_forec_track_pom_exp[::2][indx],lat_maxwind],\
                         [lon_forec_track_pom_exp[::2][indx],lon_maxwind],'km')[0][0]

    xlim = [lon_forec_track_pom_exp[::2][indx]-4,lon_forec_track_pom_exp[::2][indx]+4]
    ylim = [lat_forec_track_pom_exp[::2][indx]-4,lat_forec_track_pom_exp[::2][indx]+4]

    oklon = np.where(np.logical_and(lon_hwrf>xlim[0],lon_hwrf<xlim[1]))[0]
    oklat = np.where(np.logical_and(lat_hwrf>ylim[0],lat_hwrf<ylim[1]))[0]

    meshlon_lat = np.meshgrid(lon_hwrf[oklon],lat_hwrf[oklat])
    #meshlat = np.meshgrid(lat_hwrf[oklat],lon_hwrf[oklon])

    eye_lon = np.tile(lon_forec_track_pom_exp[::2][indx],meshlon_lat[0].shape[1])
    eye_lat = np.tile(lat_forec_track_pom_exp[::2][indx],meshlon_lat[0].shape[0])
    eye_matrix = np.meshgrid(eye_lon,eye_lat)

    lat_lon_matrix = np.stack((np.ravel(meshlon_lat[0]),np.ravel(meshlon_lat[1])),axis=1).T
    eye_lat_lon_matrix = np.stack((np.ravel(eye_matrix[0]),np.ravel(eye_matrix[1])),axis=1).T

    R = np.empty(lat_lon_matrix.shape[1])
    R[:] = np.nan
    for i in np.arange(lat_lon_matrix.shape[1]):
        R[i] = seawater.dist([lat_forec_track_pom_exp[::2][indx],lat_lon_matrix[1,i]],\
                             [lon_forec_track_pom_exp[::2][indx],lat_lon_matrix[0,i]],'km')[0][0]

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

    okRR = np.argsort(RR_norm)
    sshla_vec_interp = np.interp(RR_norm_bins,RR_norm[okRR],sshla_vec[okRR])

    # POM exp
    print(POM_exp[indx])
    pom = xr.open_dataset(POM_exp[indx])
    tpom = np.asarray(pom['time'][:])
    temp_pom = np.asarray(pom['t'][0,:,:,:])
    salt_pom = np.asarray(pom['s'][0,:,:,:])
    dens_pom = np.asarray(pom['rho'][0,:,:,:]) * 1000 + 1000
    dens_pom[dens_pom==1000.0] = np.nan

    oklon_pom = np.where(np.logical_and(lon_pom_exp[0,:]>xlim[0],lon_pom_exp[0,:]<xlim[1]))[0]
    oklat_pom = np.where(np.logical_and(lat_pom_exp[:,0]>ylim[0],lat_pom_exp[:,0]<ylim[1]))[0]

    meshlon_lat_pom = np.meshgrid(lon_pom_exp[0,oklon_pom],lat_pom_exp[oklat_pom,0])

    lat_lon_matrix_pom = np.stack((np.ravel(meshlon_lat_pom[0]),np.ravel(meshlon_lat_pom[1])),axis=1).T

    R_pom = np.empty(lat_lon_matrix_pom.shape[1])
    R_pom[:] = np.nan
    for i in np.arange(lat_lon_matrix_pom.shape[1]):
        R_pom[i] = seawater.dist([lat_forec_track_pom_exp[::2][indx],lat_lon_matrix_pom[1,i]],\
                             [lon_forec_track_pom_exp[::2][indx],lat_lon_matrix_pom[0,i]],'km')[0][0]

    R_pom_norm = R_pom/Rmax
    okR_pom_8 = R_pom_norm <= 8

    RR_pom_norm = R_pom_norm[okR_pom_8]

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

    Tmean_dens_crit_vec = np.ravel(Tmean_dens_crit)[okR_pom_8]

    okRR_pom = np.argsort(RR_pom_norm)

    kw = dict(levels=np.arange(27,29.6,0.2))
    fig,ax = plt.subplots()
    plt.contourf(lon_pom_exp[0,oklon_pom],lat_pom_exp[oklat_pom,0],Tmean_dens_crit,cmap=cmocean.cm.thermal,\
                 **kw,vmin=27,vmax=29.5)
    c=plt.colorbar()
    plt.contour(lon_pom_exp[0,oklon_pom],lat_pom_exp[oklat_pom,0],Tmean_dens_crit,\
                colors='k',linestyle='-',**kw,alpha=0.3)
    plt.plot(eye_lon[0],eye_lat[0],'ok',markersize=7,markeredgecolor='lawngreen')
    plt.plot(lon_forec_track_pom_exp,lat_forec_track_pom_exp,'.-k')
    plt.plot(lat_lon_matrix[0][okR1],lat_lon_matrix[1][okR1],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR2],lat_lon_matrix[1][okR2],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR3],lat_lon_matrix[1][okR3],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR4],lat_lon_matrix[1][okR4],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR5],lat_lon_matrix[1][okR5],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR6],lat_lon_matrix[1][okR6],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR7],lat_lon_matrix[1][okR7],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR8],lat_lon_matrix[1][okR8],'.k',markersize=1)
    plt.plot(lon_maxwind,lat_maxwind,'*k',markersize=10,markeredgecolor='lawngreen')
    plt.xlim(np.min(lat_lon_matrix[0][okR_8])-0.1,np.max(lat_lon_matrix[0][okR_8])+0.1)
    plt.ylim(np.min(lat_lon_matrix[1][okR_8])-0.1,np.max(lat_lon_matrix[1][okR_8])+0.1)
    c.set_label('$^oC$',rotation=90, labelpad=15, fontsize=14)
    plt.title('Mixed Layer Temperature ' + str(time_hwrf[0])[0:13] + '\n HWRF2020-POM (IC RTOFS)',fontsize=15)
    file_name = folder_fig + 'HWRF2020_POM_Dorian_Temp_MLD_cycle_' + cycle + '_' + str(tpom[0])[0:13]
    plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1)

    '''
    c = Tmean_dens_crit_vec
    fig,ax = plt.subplots()
    plt.scatter(lat_lon_matrix_pom[0][okR_pom_8],lat_lon_matrix_pom[1][okR_pom_8],c=c,cmap=cmocean.cm.thermal,vmin=27,vmax=29.5)
    c=plt.colorbar()
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
    plt.xlim(np.min(lat_lon_matrix[0][okR_8])-0.1,np.max(lat_lon_matrix[0][okR_8])+0.1)
    plt.ylim(np.min(lat_lon_matrix[1][okR_8])-0.1,np.max(lat_lon_matrix[1][okR_8])+0.1)
    c.set_label('$^oC$',rotation=90, labelpad=15, fontsize=14)
    plt.title('Mixed Layer Temperature ' + str(time_hwrf[0])[0:13] + '\n HWRF2020-POM (IC RTOFS)',fontsize=15)
    '''
    # for POM operational
    OHC = OHC_from_3D_fields(target_temp_POM_exp,target_depth_POM_exp,\
                             target_dens_POM_exp)*10**(-7)

    OHC_vec = np.ravel(OHC)[okR_pom_8]

    okRR_pom = np.argsort(RR_pom_norm)

    kw = dict(levels=np.arange(0,101,10))
    fig,ax = plt.subplots()
    plt.contourf(lon_pom_exp[0,oklon_pom],lat_pom_exp[oklat_pom,0],OHC,cmap='magma',\
                 **kw,vmin=0,vmax=100)
    c=plt.colorbar()
    plt.contour(lon_pom_exp[0,oklon_pom],lat_pom_exp[oklat_pom,0],OHC,\
                colors='k',linestyle='-',**kw,alpha=0.3)
    c.set_label('$KJ/cm^2$',rotation=90, labelpad=15, fontsize=14)
    plt.plot(eye_lon[0],eye_lat[0],'ok',markersize=7,markeredgecolor='lawngreen')
    plt.plot(lon_forec_track_pom_exp,lat_forec_track_pom_exp,'.-k')
    plt.plot(lat_lon_matrix[0][okR1],lat_lon_matrix[1][okR1],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR2],lat_lon_matrix[1][okR2],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR3],lat_lon_matrix[1][okR3],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR4],lat_lon_matrix[1][okR4],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR5],lat_lon_matrix[1][okR5],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR6],lat_lon_matrix[1][okR6],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR7],lat_lon_matrix[1][okR7],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR8],lat_lon_matrix[1][okR8],'.k',markersize=1)
    plt.plot(lon_maxwind,lat_maxwind,'*k',markersize=10,markeredgecolor='lawngreen')
    plt.xlim(np.min(lat_lon_matrix[0][okR_8])-0.1,np.max(lat_lon_matrix[0][okR_8])+0.1)
    plt.ylim(np.min(lat_lon_matrix[1][okR_8])-0.1,np.max(lat_lon_matrix[1][okR_8])+0.1)
    c.set_label('($KJ/cm^2$)',rotation=90, labelpad=15, fontsize=14)
    plt.title('Ocean Heat Content ' + str(time_hwrf[0])[0:13] + '\n HWRF2020-POM (IC RTOFS)',fontsize=15)
    file_name = folder_fig + 'HWRF2020_POM_Dorian_OHC_cycle_' + cycle + '_' + str(tpom[0])[0:13]
    plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1)

    '''
    c = OHC_vec
    fig,ax = plt.subplots()
    plt.scatter(lat_lon_matrix_pom[0][okR_pom_8],lat_lon_matrix_pom[1][okR_pom_8],c=c,cmap='magma',vmin=20,vmax=90)
    c=plt.colorbar()
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
    plt.xlim(np.min(lat_lon_matrix[0][okR_8])-0.1,np.max(lat_lon_matrix[0][okR_8])+0.1)
    plt.ylim(np.min(lat_lon_matrix[1][okR_8])-0.1,np.max(lat_lon_matrix[1][okR_8])+0.1)
    c.set_label('($KJ/cm^2$)',rotation=90, labelpad=15, fontsize=14)
    plt.title('Ocean Heat Content ' + str(time_hwrf[0])[0:13] + '\n HWRF2020-POM (RTOFS)',fontsize=15)
    '''
    for n in np.arange(1,len(RR_norm_bins)):
        print(n)
        okR = np.logical_and(RR_pom_norm >= RR_norm_bins[n-1],RR_pom_norm <= RR_norm_bins[n])
        if len(np.where(okR)[0])==0:
            Tmean_max_pom_exp[N,n] = np.nan
            Tmean_mean_pom_exp[N,n] = np.nan
            Tmean_min_pom_exp[N,n] = np.nan
        else:
            Tmean_max_pom_exp[N,n] = np.nanmax(Tmean_dens_crit_vec[okR])
            Tmean_mean_pom_exp[N,n] = np.nanmean(Tmean_dens_crit_vec[okR])
            Tmean_min_pom_exp[N,n] = np.nanmin(Tmean_dens_crit_vec[okR])

    for n in np.arange(1,len(RR_norm_bins)):
        print(n)
        okR = np.logical_and(RR_pom_norm >= RR_norm_bins[n-1],RR_pom_norm <= RR_norm_bins[n])
        if len(np.where(okR)[0])==0:
            OHC_max_pom_exp[N,n] = np.nan
            OHC_mean_pom_exp[N,n] = np.nan
            OHC_min_pom_exp[N,n] = np.nan
        else:
            OHC_max_pom_exp[N,n] = np.nanmax(OHC_vec[okR])
            OHC_mean_pom_exp[N,n] = np.nanmean(OHC_vec[okR])
            OHC_min_pom_exp[N,n] = np.nanmin(OHC_vec[okR])
        '''
        maxv = np.nanmean(sshla_vec[okR]) + 2*np.nanstd(sshla_vec[okR])
        minv = np.nanmean(sshla_vec[okR]) - 2*np.nanstd(sshla_vec[okR])
        okk = np.logical_and(sshla_vec[okR] >= minv,sshla_vec[okR] <= maxv)
        sshla_vec_max_pom_exp[N,n] = np.nanmax(sshla_vec[okR][okk])
        sshla_vec_mean_pom_exp[N,n] = np.nanmean(sshla_vec[okR][okk])
        sshla_vec_min_pom_exp[N,n] = np.nanmin(sshla_vec[okR][okk])
        '''

#%% MLT and OHC under Dorian HYCOM experimental

#RR_norm_bins = np.arange(0,8,0.1)
#RR_norm_bins_mid = RR_norm_bins[0:-1] + np.diff(RR_norm_bins)/2

Tmean_max_hycom_exp = np.empty((3,len(RR_norm_bins)))
Tmean_max_hycom_exp[:] = np.nan
Tmean_min_hycom_exp = np.empty((3,len(RR_norm_bins)))
Tmean_min_hycom_exp[:] = np.nan
Tmean_mean_hycom_exp = np.empty((3,len(RR_norm_bins)))
Tmean_mean_hycom_exp[:] = np.nan

OHC_max_hycom_exp = np.empty((3,len(RR_norm_bins)))
OHC_max_hycom_exp[:] = np.nan
OHC_min_hycom_exp = np.empty((3,len(RR_norm_bins)))
OHC_min_hycom_exp[:] = np.nan
OHC_mean_hycom_exp = np.empty((3,len(RR_norm_bins)))
OHC_mean_hycom_exp[:] = np.nan

time = []

afiles = sorted(glob.glob(os.path.join(folder_hycom_exp,prefix_hycom+'*.a')))

for N,indx in enumerate(np.asarray([3,11,14])):
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

    time.append(np.asarray(HWRF.variables['time'][:]))

    wind_int = np.sqrt(UGRD_hwrf**2 + VGRD_hwrf**2)
    max_wind= np.max(wind_int)
    okwind = np.where(wind_int == max_wind)
    lat_maxwind = lat_hwrf[okwind[0][0]]
    lon_maxwind = lon_hwrf[okwind[1][0]]
    Rmax = seawater.dist([lat_forec_track_hycom_exp[::2][indx],lat_maxwind],\
                         [lon_forec_track_hycom_exp[::2][indx],lon_maxwind],'km')[0][0]

    xlim = [lon_forec_track_hycom_exp[::2][indx]-4,lon_forec_track_hycom_exp[::2][indx]+4]
    ylim = [lat_forec_track_hycom_exp[::2][indx]-4,lat_forec_track_hycom_exp[::2][indx]+4]

    oklon = np.where(np.logical_and(lon_hwrf>xlim[0],lon_hwrf<xlim[1]))[0]
    oklat = np.where(np.logical_and(lat_hwrf>ylim[0],lat_hwrf<ylim[1]))[0]

    meshlon_lat = np.meshgrid(lon_hwrf[oklon],lat_hwrf[oklat])

    eye_lon = np.tile(lon_forec_track_hycom_exp[::2][indx],meshlon_lat[0].shape[1])
    eye_lat = np.tile(lat_forec_track_hycom_exp[::2][indx],meshlon_lat[0].shape[0])
    eye_matrix = np.meshgrid(eye_lon,eye_lat)

    lat_lon_matrix = np.stack((np.ravel(meshlon_lat[0]),np.ravel(meshlon_lat[1])),axis=1).T
    eye_lat_lon_matrix = np.stack((np.ravel(eye_matrix[0]),np.ravel(eye_matrix[1])),axis=1).T

    R = np.empty(lat_lon_matrix.shape[1])
    R[:] = np.nan
    for i in np.arange(lat_lon_matrix.shape[1]):
        R[i] = seawater.dist([lat_forec_track_hycom_exp[::2][indx],lat_lon_matrix[1,i]],\
                             [lon_forec_track_hycom_exp[::2][indx],lat_lon_matrix[0,i]],'km')[0][0]

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

    okRR = np.argsort(RR_norm)
    sshla_vec_interp = np.interp(RR_norm_bins,RR_norm[okRR],sshla_vec[okRR])

    # Read HYCOM exp
    file = afiles[indx]

    #Reading time stamp
    year = int(file.split('/')[-1].split('.')[1][0:4])
    month = int(file.split('/')[-1].split('.')[1][4:6])
    day = int(file.split('/')[-1].split('.')[1][6:8])
    hour = int(file.split('/')[-1].split('.')[1][8:10])
    dt = int(file.split('/')[-1].split('.')[3][1:])
    timestamp_hycom = mdates.date2num(datetime(year,month,day,hour)) + dt/24
    #time_hycom.append(mdates.num2date(timestamp_hycom))

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
        R_hycom[i] = seawater.dist([lat_forec_track_hycom_exp[::2][indx],lat_lon_matrix_hycom[1,i]],\
                             [lon_forec_track_hycom_exp[::2][indx],lat_lon_matrix_hycom[0,i]],'km')[0][0]

    R_hycom_norm = R_hycom/Rmax
    okR_hycom_8 = R_hycom_norm <= 8

    RR_hycom_norm = R_hycom_norm[okR_hycom_8]

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

    Tmean_dens_crit_vec = np.ravel(Tmean_dens_crit)[okR_hycom_8]

    okRR_hycom = np.argsort(RR_hycom_norm)

    Tmean_dens_crit[Tmean_dens_crit<27] = 27
    kw = dict(levels=np.arange(27,29.6,0.2))
    fig,ax = plt.subplots()
    plt.contourf(lon_hycom[0,oklon_hycom]-360,lat_hycom[oklat_hycom,0],Tmean_dens_crit,cmap=cmocean.cm.thermal,\
                 **kw,vmin=27,vmax=29.5)
    c=plt.colorbar()
    plt.contour(lon_hycom[0,oklon_hycom]-360,lat_hycom[oklat_hycom,0],Tmean_dens_crit,\
                colors='k',linestyle='-',**kw,alpha=0.3)
    plt.plot(eye_lon[0],eye_lat[0],'ok',markersize=7,markeredgecolor='lawngreen')
    plt.plot(lon_forec_track_hycom_exp,lat_forec_track_hycom_exp,'.-k')
    plt.plot(lat_lon_matrix[0][okR1],lat_lon_matrix[1][okR1],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR2],lat_lon_matrix[1][okR2],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR3],lat_lon_matrix[1][okR3],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR4],lat_lon_matrix[1][okR4],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR5],lat_lon_matrix[1][okR5],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR6],lat_lon_matrix[1][okR6],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR7],lat_lon_matrix[1][okR7],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR8],lat_lon_matrix[1][okR8],'.k',markersize=1)
    plt.plot(lon_maxwind,lat_maxwind,'*k',markersize=10,markeredgecolor='lawngreen')
    plt.xlim(np.min(lat_lon_matrix[0][okR_8])-0.1,np.max(lat_lon_matrix[0][okR_8])+0.1)
    plt.ylim(np.min(lat_lon_matrix[1][okR_8])-0.1,np.max(lat_lon_matrix[1][okR_8])+0.1)
    c.set_label('$^oC$',rotation=90, labelpad=15, fontsize=14)
    plt.title('Mixed Layer Temperature ' + str(time_hwrf[0])[0:13] + '\n HWRF2020-HYCOM (IC RTOFS)',fontsize=15)
    file_name = folder_fig + 'HWRF2020_HYCOM_Dorian_Temp_MLD_cycle_' + cycle + '_' + str(time_hwrf[0])[0:13]
    plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1)

    '''
    c = Tmean_dens_crit_vec
    fig,ax = plt.subplots()
    plt.scatter(lat_lon_matrix_hycom[0][okR_hycom_8]-360,lat_lon_matrix_hycom[1][okR_hycom_8],c=c,cmap=cmocean.cm.thermal,vmin=27,vmax=29.5)
    c=plt.colorbar()
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
    plt.xlim(np.min(lat_lon_matrix[0][okR_8])-0.1,np.max(lat_lon_matrix[0][okR_8])+0.1)
    plt.ylim(np.min(lat_lon_matrix[1][okR_8])-0.1,np.max(lat_lon_matrix[1][okR_8])+0.1)
    c.set_label('$^oC$',rotation=90, labelpad=15, fontsize=14)
    plt.title('Mixed Layer Temperature ' + str(time_hwrf[0])[0:13] + '\n HWRF2020-HYCOM (IC RTOFS)',fontsize=15)
    '''
    # for HYCOM exp
    OHC = OHC_from_3D_fields(target_temp_hycom_exp.T,target_depth_hycom_exp,\
                             target_dens_hycom_exp.T)*10**(-7)
    OHC = OHC.T

    OHC_vec = np.ravel(OHC)[okR_hycom_8]

    okRR_hycom = np.argsort(RR_hycom_norm)

    kw = dict(levels=np.arange(0,101,10))
    fig,ax = plt.subplots()
    plt.contourf(lon_hycom[0,oklon_hycom]-360,lat_hycom[oklat_hycom,0],OHC,cmap='magma',\
                 **kw,vmin=0,vmax=100)
    c=plt.colorbar()
    plt.contour(lon_hycom[0,oklon_hycom]-360,lat_hycom[oklat_hycom,0],OHC,\
                colors='k',linestyle='-',**kw,alpha=0.3)
    c.set_label('$KJ/cm^2$',rotation=90, labelpad=15, fontsize=14)
    plt.plot(eye_lon[0],eye_lat[0],'ok',markersize=7,markeredgecolor='lawngreen')
    plt.plot(lon_forec_track_hycom_exp,lat_forec_track_hycom_exp,'.-k')
    plt.plot(lat_lon_matrix[0][okR1],lat_lon_matrix[1][okR1],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR2],lat_lon_matrix[1][okR2],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR3],lat_lon_matrix[1][okR3],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR4],lat_lon_matrix[1][okR4],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR5],lat_lon_matrix[1][okR5],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR6],lat_lon_matrix[1][okR6],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR7],lat_lon_matrix[1][okR7],'.k',markersize=1)
    plt.plot(lat_lon_matrix[0][okR8],lat_lon_matrix[1][okR8],'.k',markersize=1)
    plt.plot(lon_maxwind,lat_maxwind,'*k',markersize=10,markeredgecolor='lawngreen')
    plt.xlim(np.min(lat_lon_matrix[0][okR_8])-0.1,np.max(lat_lon_matrix[0][okR_8])+0.1)
    plt.ylim(np.min(lat_lon_matrix[1][okR_8])-0.1,np.max(lat_lon_matrix[1][okR_8])+0.1)
    c.set_label('($KJ/cm^2$)',rotation=90, labelpad=15, fontsize=14)
    plt.title('Ocean Heat Content ' + str(time_hwrf[0])[0:13] + '\n HWRF2020-HYCOM (IC RTOFS)',fontsize=15)
    file_name = folder_fig + 'HWRF2020_HYCOM_Dorian_OHC_cycle_' + cycle + '_' + str(time_hwrf[0])[0:13]
    plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1)

    '''
    c = OHC_vec
    fig,ax = plt.subplots()
    plt.scatter(lat_lon_matrix_hycom[0][okR_hycom_8]-360,lat_lon_matrix_hycom[1][okR_hycom_8],c=c,cmap='magma',vmin=20,vmax=90)
    c=plt.colorbar()
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
    plt.xlim(np.min(lat_lon_matrix[0][okR_8])-0.1,np.max(lat_lon_matrix[0][okR_8])+0.1)
    plt.ylim(np.min(lat_lon_matrix[1][okR_8])-0.1,np.max(lat_lon_matrix[1][okR_8])+0.1)
    c.set_label('($KJ/cm^2$)',rotation=90, labelpad=15, fontsize=14)
    plt.title('Ocean Heat Content ' + str(time_hwrf[0])[0:13] + '\n HWRF2020-HYCOM (RTOFS)',fontsize=15)
    '''
    for n in np.arange(1,len(RR_norm_bins)):
        print(n)
        okR = np.logical_and(RR_hycom_norm >= RR_norm_bins[n-1],RR_hycom_norm <= RR_norm_bins[n])
        if len(np.where(okR)[0])==0:
            Tmean_max_hycom_exp[N,n] = np.nan
            Tmean_mean_hycom_exp[N,n] = np.nan
            Tmean_min_hycom_exp[N,n] = np.nan
        else:
            Tmean_max_hycom_exp[N,n] = np.nanmax(Tmean_dens_crit_vec[okR])
            Tmean_mean_hycom_exp[N,n] = np.nanmean(Tmean_dens_crit_vec[okR])
            Tmean_min_hycom_exp[N,n] = np.nanmin(Tmean_dens_crit_vec[okR])

    for n in np.arange(1,len(RR_norm_bins)):
        print(n)
        okR = np.logical_and(RR_hycom_norm >= RR_norm_bins[n-1],RR_hycom_norm <= RR_norm_bins[n])
        if len(np.where(okR)[0])==0:
            OHC_max_hycom_exp[N,n] = np.nan
            OHC_mean_hycom_exp[N,n] = np.nan
            OHC_min_hycom_exp[N,n] = np.nan
        else:
            OHC_max_hycom_exp[N,n] = np.nanmax(OHC_vec[okR])
            OHC_mean_hycom_exp[N,n] = np.nanmean(OHC_vec[okR])
            OHC_min_hycom_exp[N,n] = np.nanmin(OHC_vec[okR])
        '''
        maxv = np.nanmean(sshla_vec[okR]) + 2*np.nanstd(sshla_vec[okR])
        minv = np.nanmean(sshla_vec[okR]) - 2*np.nanstd(sshla_vec[okR])
        okk = np.logical_and(sshla_vec[okR] >= minv,sshla_vec[okR] <= maxv)
        sshla_vec_max_pom_exp[N,n] = np.nanmax(sshla_vec[okR][okk])
        sshla_vec_mean_pom_exp[N,n] = np.nanmean(sshla_vec[okR][okk])
        sshla_vec_min_pom_exp[N,n] = np.nanmin(sshla_vec[okR][okk])
        '''

#%%
for N in np.arange(3):
    fig,ax = plt.subplots()
    plt.plot(RR_norm_bins_mid,Tmean_mean_pom_oper[N,1:],'X-',color='mediumorchid',label='HWRF2019-POM (IC clim.)',markeredgecolor='k',markersize=7)
    #plt.plot(RR_norm_bins_mid,Tmean_max_pom_oper[N,1:],'-',color='mediumorchid',alpha=0.5)
    #plt.plot(RR_norm_bins_mid,Tmean_min_pom_oper[N,1:],'-',color='mediumorchid',alpha=0.5)
    ax.fill_between(RR_norm_bins_mid,Tmean_min_pom_oper[N,1:],Tmean_max_pom_oper[N,1:],color='mediumorchid',alpha=0.1)

    plt.plot(RR_norm_bins_mid,Tmean_mean_pom_exp[N,1:],'^-',color='teal',label='HWRF2020-POM (RTOFS)',markeredgecolor='k',markersize=7)
    #plt.plot(RR_norm_bins_mid,Tmean_max_pom_exp[N,1:],'-',color='teal',alpha=0.5)
    #plt.plot(RR_norm_bins_mid,Tmean_min_pom_exp[N,1:],'-',color='teal',alpha=0.5)
    ax.fill_between(RR_norm_bins_mid,Tmean_min_pom_exp[N,1:],Tmean_max_pom_exp[N,1:],color='teal',alpha=0.1)

    plt.plot(RR_norm_bins_mid,Tmean_mean_hycom_exp[N,1:],'H-',color='darkorange',label='HWRF2020-HYCOM (RTOFS)',markeredgecolor='k',markersize=7)
    #plt.plot(RR_norm_bins_mid,Tmean_max_hycom_exp[N,1:],'-',color='darkorange',alpha=0.5)
    #plt.plot(RR_norm_bins_mid,Tmean_min_hycom_exp[N,1:],'-',color='darkorange',alpha=0.5)
    ax.fill_between(RR_norm_bins_mid,Tmean_min_hycom_exp[N,1:],Tmean_max_hycom_exp[N,1:],color='darkorange',alpha=0.1)

    plt.title('Mixed Layer Temperature ' + str(time[N])[2:15],fontsize=16)
    plt.ylim(27.0,29.5)
    plt.ylabel('$^oC$',fontsize=14)
    plt.xlabel('r/Rmax',fontsize=14)
    plt.xticks(np.arange(9))
    plt.legend(loc='lower left')
    file_name = folder_fig + 'Dorian_Temp_MLD_vs_Rnorm_cycle_' + cycle + '_' + str(time[N])[0:13]
    plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1)

#%%
for N in np.arange(3):
    fig,ax = plt.subplots()
    plt.plot(RR_norm_bins_mid,OHC_mean_pom_oper[N,1:],'X-',color='mediumorchid',label='HWRF2019-POM (IC clim.)',markeredgecolor='k',markersize=7)
    #plt.plot(RR_norm_bins_mid,OHC_max_pom_oper[N,1:],'-',color='mediumorchid',alpha=0.5)
    #plt.plot(RR_norm_bins_mid,OHC_min_pom_oper[N,1:],'-',color='mediumorchid',alpha=0.5)
    ax.fill_between(RR_norm_bins_mid,OHC_min_pom_oper[N,1:],OHC_max_pom_oper[N,1:],color='mediumorchid',alpha=0.1)

    plt.plot(RR_norm_bins_mid,OHC_mean_pom_exp[N,1:],'^-',color='teal',label='HWRF2020-POM (RTOFS)',markeredgecolor='k',markersize=7)
    #plt.plot(RR_norm_bins_mid,OHC_max_pom_exp[N,1:],'-',color='teal',alpha=0.5)
    #plt.plot(RR_norm_bins_mid,OHC_min_pom_exp[N,1:],'-',color='teal',alpha=0.5)
    ax.fill_between(RR_norm_bins_mid,OHC_min_pom_exp[N,1:],OHC_max_pom_exp[N,1:],color='teal',alpha=0.1)

    plt.plot(RR_norm_bins_mid,OHC_mean_hycom_exp[N,1:],'H-',color='darkorange',label='HWRF2020-HYCOM (RTOFS)',markeredgecolor='k',markersize=7)
    #plt.plot(RR_norm_bins_mid,OHC_max_hycom_exp[N,1:],'-',color='darkorange',alpha=0.5)
    #plt.plot(RR_norm_bins_mid,OHC_min_hycom_exp[N,1:],'-',color='darkorange',alpha=0.5)
    ax.fill_between(RR_norm_bins_mid,OHC_min_hycom_exp[N,1:],OHC_max_hycom_exp[N,1:],color='darkorange',alpha=0.1)

    plt.title('Ocean Heat Content ' + str(time[N])[2:15],fontsize=16)
    plt.ylim(0,100)
    plt.ylabel('$Kj/cm^2$',fontsize=14)
    plt.xlabel('r/Rmax',fontsize=14)
    plt.xticks(np.arange(9))
    plt.legend(loc='lower left')
    file_name = folder_fig + 'Dorian_OHC_vs_Rnorm_cycle_' + cycle + '_' + str(time[N])[0:13]
    plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1)
