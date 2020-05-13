#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 09:37:22 2020

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
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)

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

#%% map wind vectors HWRP2019-POM operational

N = 5
#kw = dict(levels=np.linspace(0,900,10))
kw = dict(levels=np.linspace(0,1100,12))
ylim = [13,25]
xlim = [-72,-63]

'''
N = 22
kw = dict(levels=np.linspace(0,1100,12))
ylim = [16,28]
xlim = [-74,-65]
'''

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
c = plt.colorbar()
q = plt.quiver(lon_hwrf[::30], lat_hwrf[::30],UGRD_hwrf[::30,::30],VGRD_hwrf[::30,::30]) #,units='xy' ,scale=0.01)
#plt.quiverkey(q,-61,25.5,30,"30 m/s",coordinates='data',color='k',fontproperties={'size': 14})
plt.quiverkey(q,-63.0,28.5,30,"30 m/s",coordinates='data',color='k',fontproperties={'size': 14})
plt.title('HWRF2019-POM Enthalpy Flux \n on '+str(t_hwrf[0])[0:13],fontsize=14)
c.set_label('$W/m^2$',rotation=90, labelpad=1, fontsize=16)
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
c = plt.colorbar()
q = plt.quiver(lon_hwrf[::30], lat_hwrf[::30],UGRD_hwrf[::30,::30],VGRD_hwrf[::30,::30]) #,units='xy' ,scale=0.01)
#plt.quiverkey(q,-61,25.5,30,"30 m/s",coordinates='data',color='k',fontproperties={'size': 14})
plt.quiverkey(q,-63.0,28.5,30,"30 m/s",coordinates='data',color='k',fontproperties={'size': 14})
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
c = plt.colorbar()
q = plt.quiver(lon_hwrf[::30], lat_hwrf[::30],UGRD_hwrf[::30,::30],VGRD_hwrf[::30,::30]) #,units='xy' ,scale=0.01)
#plt.quiverkey(q,-61,25.5,30,"30 m/s",coordinates='data',color='k',fontproperties={'size': 14})
plt.quiverkey(q,-63.0,28.5,30,"30 m/s",coordinates='data',color='k',fontproperties={'size': 14})
plt.title('HWRF2020-HYCOM Enthalpy Flux \n on '+str(t_hwrf[0])[0:13],fontsize=14)
c.set_label('$W/m^2$',rotation=90, labelpad=1, fontsize=16)
c.ax.tick_params(labelsize=14)
plt.axis('scaled')
plt.ylim([ylim[0],ylim[1]])
plt.xlim([xlim[0],xlim[1]])

file_name = folder_fig + 'Dorian_heat_fluxes_HWRF20_HYCOM_' + str(t_hwrf[0])[0:13] + '_' + cycle
plt.savefig(file_name)#,bbox_inches = 'tight',pad_inches = 0.1) 
