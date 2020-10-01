#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 14:55:54 2020

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
kmz_file_Dorian = '/home/aristizabal/KMZ_files/al052019_best_track-5.kmz'

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
from datetime import datetime
import os
import os.path
import glob
#import netCDF4
#from bs4 import BeautifulSoup
#from zipfile import ZipFile
#from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
import cmocean
import matplotlib.patches as patches   

import sys
sys.path.append('/home/aristizabal/NCEP_scripts/')
from utils4HYCOM import readBinz, readgrids

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

#%% Read POM grid

print('Retrieving coordinates from POM')
POM_grid_oper = xr.open_dataset(pom_grid_oper,decode_times=False)
lon_pom_oper = np.asarray(POM_grid_oper['east_e'][:])
lat_pom_oper = np.asarray(POM_grid_oper['north_e'][:])
zlev_pom_oper = np.asarray(POM_grid_oper['zz'][:])
hpom_oper = np.asarray(POM_grid_oper['h'][:])
zmatrix = np.dot(hpom_oper.reshape(-1,1),zlev_pom_oper.reshape(1,-1))
zmatrix_pom_oper = zmatrix.reshape(hpom_oper.shape[0],hpom_oper.shape[1],zlev_pom_oper.shape[0])

POM_grid_exp = xr.open_dataset(pom_grid_exp,decode_times=False)
lon_pom_exp = np.asarray(POM_grid_exp['east_e'][:])
lat_pom_exp = np.asarray(POM_grid_exp['north_e'][:])
zlev_pom_exp = np.asarray(POM_grid_exp['zz'][:])
hpom_exp = np.asarray(POM_grid_exp['h'][:])
zmatrix = np.dot(hpom_exp.reshape(-1,1),zlev_pom_exp.reshape(1,-1))
zmatrix_pom_exp = zmatrix.reshape(hpom_exp.shape[0],hpom_exp.shape[1],zlev_pom_exp.shape[0])
        
#%% Reading HYCOM grid

print('Retrieving coordinates from HYCOM')
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

#%% Get list ocean files

#POM_oper = sorted(glob.glob(os.path.join(folder_pom_oper,'*.nc')))
#POM_exp = sorted(glob.glob(os.path.join(folder_pom_exp,'*.nc')))
#HYCOM_exp = sorted(glob.glob(os.path.join(folder_hycom20,'*.nc')))

#%% Surface temperature POM operational

folder_pom =folder_pom_oper
prefix = prefix_pom
N = 0
ylim = [17,18]
xlim = [-65.5,-64.5]

pom_ncfiles = sorted(glob.glob(os.path.join(folder_pom,prefix+'*.nc')))   
file = pom_ncfiles[N]
pom = xr.open_dataset(file)
time_pom = np.asarray(pom['time'])
temp_POM = np.asarray(pom['t'][0,0,:,:])
temp_POM[temp_POM==0] = np.nan
#temp_POM[temp_POM<27] = 27
#temp_POM[temp_POM>32] = 32

kw = dict(levels = np.linspace(28,30,41))
#kw = dict(levels = np.linspace(26,32,17))
#fig,ax = plt.subplots(figsize=(5,5)) 
plt.figure()
plt.contourf(lon_pom_oper[0,:],lat_pom_oper[:,0],temp_POM,cmap=cmocean.cm.thermal,**kw)
#plt.xlim([np.min(lon_forec_track_pom_oper)-0.5,np.max(lon_forec_track_pom_oper)+0.5])
#plt.ylim([np.min(lat_forec_track_pom_oper)-0.5,np.max(lat_forec_track_pom_oper)+0.5])
plt.colorbar()
plt.title('HWRF2019-POM SST \n on '+str(time_pom[0])[0:13],fontsize=14)
plt.axis('scaled')
plt.ylim([ylim[0],ylim[1]])
plt.xlim([xlim[0],xlim[1]])

#%% Surface temperature POM experimental

folder_pom =folder_pom_exp
prefix = prefix_pom
N = 4
ylim = [17,18]
xlim = [-65.5,-64.5]

pom_ncfiles = sorted(glob.glob(os.path.join(folder_pom,prefix+'*.nc')))   
file = pom_ncfiles[N]
pom = xr.open_dataset(file)
time_pom = np.asarray(pom['time'])
temp_POM = np.asarray(pom['t'][0,0,:,:])
temp_POM[temp_POM==0] = np.nan
#temp_POM[temp_POM<27] = 27
#temp_POM[temp_POM>32] = 32

kw = dict(levels = np.linspace(28,30,41))
#kw = dict(levels = np.linspace(26,32,17))
#fig,ax = plt.subplots(figsize=(5,5)) 
plt.figure()
plt.contourf(lon_pom_exp[0,:],lat_pom_exp[:,0],temp_POM,cmap=cmocean.cm.thermal,**kw)
#plt.xlim([np.min(lon_forec_track_pom_oper)-0.5,np.max(lon_forec_track_pom_oper)+0.5])
#plt.ylim([np.min(lat_forec_track_pom_oper)-0.5,np.max(lat_forec_track_pom_oper)+0.5])
plt.colorbar()
plt.title('HWRF2020-POM SST \n on '+str(time_pom[0])[0:13],fontsize=14)
plt.axis('scaled')
plt.ylim([ylim[0],ylim[1]])
plt.xlim([xlim[0],xlim[1]])

#%% Surface temperature HYCOM experimental

folder_hycom =folder_hycom_exp
prefix = prefix_hycom
N = 0
var = 'temp'

afiles = sorted(glob.glob(os.path.join(folder_hycom,prefix+'*.a')))    
file = afiles[N]

# Reading 3D variable from binary file 
var_hyc = readBinz(file[:-2],'3z',var)
temp_HYCOM = var_hyc[:,:,0]

#Reading time stamp
year = int(file.split('/')[-1].split('.')[1][0:4])
month = int(file.split('/')[-1].split('.')[1][4:6])
day = int(file.split('/')[-1].split('.')[1][6:8])
hour = int(file.split('/')[-1].split('.')[1][8:10])
dt = int(file.split('/')[-1].split('.')[3][1:])
timestamp_hycom = mdates.date2num(datetime(year,month,day,hour)) + dt/24
time_hycom = mdates.num2date(timestamp_hycom)

kw = dict(levels = np.linspace(28,30,41))
#fig,ax = plt.subplots(figsize=(5,5)) 
plt.figure()
plt.contourf(lon_hycom[0,:]-360,lat_hycom[:,0],temp_HYCOM,cmap=cmocean.cm.thermal,**kw)
#plt.xlim([np.min(lon_forec_track_pom_oper[oktt])-0.5,np.max(lon_forec_track_pom_oper[oktt])+0.5])
#plt.ylim([np.min(lat_forec_track_pom_oper[oktt])-0.5,np.max(lat_forec_track_pom_oper[oktt])+0.5])
plt.colorbar()
plt.title('HWRF2020-HYCOM SST \n on '+str(time_hycom)[0:13],fontsize=14)
plt.axis('scaled')
plt.ylim([ylim[0],ylim[1]])
plt.xlim([xlim[0],xlim[1]])

#%% Delta Surface temperature POM operational

ylim = [17,18]
xlim = [-65.5,-64.5]

folder_pom =folder_pom_oper
prefix = prefix_pom

N = 0
pom_ncfiles = sorted(glob.glob(os.path.join(folder_pom,prefix+'*.nc')))   
file = pom_ncfiles[N]
pom = xr.open_dataset(file)
time_pom1 = np.asarray(pom['time'])
temp_POM1 = np.asarray(pom['t'][0,0,:,:])
temp_POM1[temp_POM1==0] = np.nan

N = 4
pom_ncfiles = sorted(glob.glob(os.path.join(folder_pom,prefix+'*.nc')))   
file = pom_ncfiles[N]
pom = xr.open_dataset(file)
time_pom4 = np.asarray(pom['time'])
temp_POM4 = np.asarray(pom['t'][0,0,:,:])
temp_POM4[temp_POM4==0] = np.nan

delta_temp = temp_POM4 - temp_POM1

kw = dict(levels = np.linspace(-2,2,21))
#fig,ax = plt.subplots(figsize=(5,5)) 
fig, ax = plt.subplots()
plt.contourf(lon_pom_oper[0,:],lat_pom_oper[:,0],delta_temp,cmap=cmocean.cm.curl,**kw)
plt.colorbar()
plt.title('HWRF2019-POM Delta SST \n'+\
          str(time_pom4[0])[0:13]+'-'+str(time_pom1[0])[0:13],fontsize=14)
plt.axis('scaled')
rect = patches.Rectangle((xlim[0],ylim[0]),\
                 xlim[1]-xlim[0],ylim[1]-ylim[0],\
                 linewidth=1,edgecolor='k',facecolor='none')
ax.add_patch(rect)

plt.xlim([-80,-60])
plt.ylim([15,35])

plt.ylim([ylim[0],ylim[1]])
plt.xlim([xlim[0],xlim[1]])

#%% Delta Surface temperature POM experimental

ylim = [17,18]
xlim = [-65.5,-64.5]

folder_pom =folder_pom_exp
prefix = prefix_pom

N = 0
pom_ncfiles = sorted(glob.glob(os.path.join(folder_pom,prefix+'*.nc')))   
file = pom_ncfiles[N]
pom = xr.open_dataset(file)
time_pom1 = np.asarray(pom['time'])
temp_POM1 = np.asarray(pom['t'][0,0,:,:])
temp_POM1[temp_POM1==0] = np.nan

N = 4
pom_ncfiles = sorted(glob.glob(os.path.join(folder_pom,prefix+'*.nc')))   
file = pom_ncfiles[N]
pom = xr.open_dataset(file)
time_pom4 = np.asarray(pom['time'])
temp_POM4 = np.asarray(pom['t'][0,0,:,:])
temp_POM4[temp_POM4==0] = np.nan

delta_temp = temp_POM4 - temp_POM1

kw = dict(levels = np.linspace(-2,2,21))
#fig,ax = plt.subplots(figsize=(5,5)) 
fig, ax = plt.subplots()
plt.contourf(lon_pom_exp[0,:],lat_pom_exp[:,0],delta_temp,cmap=cmocean.cm.curl,**kw)
plt.colorbar()
plt.title('HWRF2020-POM Delta SST \n'+\
          str(time_pom4[0])[0:13]+'-'+str(time_pom1[0])[0:13],fontsize=14)
rect = patches.Rectangle((xlim[0],ylim[0]),\
                 xlim[1]-xlim[0],ylim[1]-ylim[0],\
                 linewidth=1,edgecolor='k',facecolor='none')
ax.add_patch(rect)

plt.xlim([-80,-60])
plt.ylim([15,35])

plt.ylim([ylim[0],ylim[1]])
plt.xlim([xlim[0],xlim[1]])

#%% Delta Surface temperature HYCOM experimental

ylim = [17,18]
xlim = [-65.5,-64.5]

folder_hycom =folder_hycom_exp
prefix = prefix_hycom

N = 0
var = 'temp'
afiles = sorted(glob.glob(os.path.join(folder_hycom,prefix+'*.a')))    
file = afiles[N]
# Reading 3D variable from binary file 
var_hyc = readBinz(file[:-2],'3z',var)
temp_hycom1 = var_hyc[:,:,0]
#Reading time stamp
year = int(file.split('/')[-1].split('.')[1][0:4])
month = int(file.split('/')[-1].split('.')[1][4:6])
day = int(file.split('/')[-1].split('.')[1][6:8])
hour = int(file.split('/')[-1].split('.')[1][8:10])
dt = int(file.split('/')[-1].split('.')[3][1:])
timestamp_hycom = mdates.date2num(datetime(year,month,day,hour)) + dt/24
time_hycom1 = mdates.num2date(timestamp_hycom)

N = 4
var = 'temp'
afiles = sorted(glob.glob(os.path.join(folder_hycom,prefix+'*.a')))    
file = afiles[N]
# Reading 3D variable from binary file 
var_hyc = readBinz(file[:-2],'3z',var)
temp_hycom4 = var_hyc[:,:,0]
#Reading time stamp
year = int(file.split('/')[-1].split('.')[1][0:4])
month = int(file.split('/')[-1].split('.')[1][4:6])
day = int(file.split('/')[-1].split('.')[1][6:8])
hour = int(file.split('/')[-1].split('.')[1][8:10])
dt = int(file.split('/')[-1].split('.')[3][1:])
timestamp_hycom = mdates.date2num(datetime(year,month,day,hour)) + dt/24
time_hycom4 = mdates.num2date(timestamp_hycom)

delta_temp = temp_hycom4 - temp_hycom1
delta_temp[delta_temp<-0.5] = -0.5

kw = dict(levels = np.linspace(-2,2,21))
#fig,ax = plt.subplots(figsize=(5,5)) 
fig, ax = plt.subplots()
plt.contourf(lon_hycom[0,:]-360,lat_hycom[:,0],delta_temp,cmap=cmocean.cm.curl,**kw)
plt.colorbar()
plt.title('HWRF2020-HYCOM Delta SST \n'+\
          str(time_hycom4)[0:13]+'-'+ str(time_hycom1)[0:13],fontsize=14)
plt.axis('scaled')
rect = patches.Rectangle((xlim[0],ylim[0]),\
                 xlim[1]-xlim[0],ylim[1]-ylim[0],\
                 linewidth=1,edgecolor='k',facecolor='none')
ax.add_patch(rect)

plt.xlim([-80,-60])
plt.ylim([15,35])

plt.ylim([ylim[0],ylim[1]])
plt.xlim([xlim[0],xlim[1]])

#%% Delta Surface temperature POM operational

ylim = [24.6,26.1]
xlim = [-72.5,-71]

folder_pom =folder_pom_oper
prefix = prefix_pom

N = 9
pom_ncfiles = sorted(glob.glob(os.path.join(folder_pom,prefix+'*.nc')))   
file = pom_ncfiles[N]
pom = xr.open_dataset(file)
time_pomi = np.asarray(pom['time'])
temp_POMi = np.asarray(pom['t'][0,0,:,:])
temp_POMi[temp_POMi==0] = np.nan

N = 13
pom_ncfiles = sorted(glob.glob(os.path.join(folder_pom,prefix+'*.nc')))   
file = pom_ncfiles[N]
pom = xr.open_dataset(file)
time_pomf = np.asarray(pom['time'])
temp_POMf = np.asarray(pom['t'][0,0,:,:])
temp_POMf[temp_POMf==0] = np.nan

delta_temp = temp_POMf - temp_POMi

kw = dict(levels = np.linspace(-2,2,21))
#fig,ax = plt.subplots(figsize=(5,5)) 
fig, ax = plt.subplots()
plt.contourf(lon_pom_oper[0,:],lat_pom_oper[:,0],delta_temp,cmap=cmocean.cm.curl,**kw)
plt.colorbar()
plt.title('HWRF2019-POM Delta SST \n'+\
          str(time_pomf[0])[0:13]+'-'+str(time_pomi[0])[0:13],fontsize=14)
plt.axis('scaled')
rect = patches.Rectangle((xlim[0],ylim[0]),\
                 xlim[1]-xlim[0],ylim[1]-ylim[0],\
                 linewidth=1,edgecolor='k',facecolor='none')
ax.add_patch(rect)

plt.xlim([-80,-60])
plt.ylim([15,35])

plt.ylim([ylim[0],ylim[1]])
plt.xlim([xlim[0],xlim[1]])

#%% Delta Surface temperature POM experimental

ylim = [24.6,26.1]
xlim = [-72.5,-71]

folder_pom =folder_pom_exp
prefix = prefix_pom

N = 9
pom_ncfiles = sorted(glob.glob(os.path.join(folder_pom,prefix+'*.nc')))   
file = pom_ncfiles[N]
pom = xr.open_dataset(file)
time_pomi = np.asarray(pom['time'])
temp_POMi = np.asarray(pom['t'][0,0,:,:])
temp_POMi[temp_POMi==0] = np.nan

N = 13
pom_ncfiles = sorted(glob.glob(os.path.join(folder_pom,prefix+'*.nc')))   
file = pom_ncfiles[N]
pom = xr.open_dataset(file)
time_pomf = np.asarray(pom['time'])
temp_POMf = np.asarray(pom['t'][0,0,:,:])
temp_POMf[temp_POMf==0] = np.nan

delta_temp = temp_POMf - temp_POMi

kw = dict(levels = np.linspace(-2,2,21))
#fig,ax = plt.subplots(figsize=(5,5)) 
fig, ax = plt.subplots()
plt.contourf(lon_pom_exp[0,:],lat_pom_exp[:,0],delta_temp,cmap=cmocean.cm.curl,**kw)
plt.colorbar()
plt.title('HWRF2020-POM Delta SST \n'+\
          str(time_pomf[0])[0:13]+'-'+str(time_pomi[0])[0:13],fontsize=14)
plt.axis('scaled')
rect = patches.Rectangle((xlim[0],ylim[0]),\
                 xlim[1]-xlim[0],ylim[1]-ylim[0],\
                 linewidth=1,edgecolor='k',facecolor='none')
ax.add_patch(rect)

plt.xlim([-80,-60])
plt.ylim([15,35])

plt.ylim([ylim[0],ylim[1]])
plt.xlim([xlim[0],xlim[1]])

#%% Delta Surface temperature HYCOM experimental

ylim = [24.6,26.1]
xlim = [-72.5,-71]

folder_hycom =folder_hycom_exp
prefix = prefix_hycom

N = 9
var = 'temp'
afiles = sorted(glob.glob(os.path.join(folder_hycom,prefix+'*.a')))    
file = afiles[N]
# Reading 3D variable from binary file 
var_hyc = readBinz(file[:-2],'3z',var)
temp_hycomi = var_hyc[:,:,0]
#Reading time stamp
year = int(file.split('/')[-1].split('.')[1][0:4])
month = int(file.split('/')[-1].split('.')[1][4:6])
day = int(file.split('/')[-1].split('.')[1][6:8])
hour = int(file.split('/')[-1].split('.')[1][8:10])
dt = int(file.split('/')[-1].split('.')[3][1:])
timestamp_hycom = mdates.date2num(datetime(year,month,day,hour)) + dt/24
time_hycomi = mdates.num2date(timestamp_hycom)

N = 13
var = 'temp'
afiles = sorted(glob.glob(os.path.join(folder_hycom,prefix+'*.a')))    
file = afiles[N]
# Reading 3D variable from binary file 
var_hyc = readBinz(file[:-2],'3z',var)
temp_hycomf = var_hyc[:,:,0]
#Reading time stamp
year = int(file.split('/')[-1].split('.')[1][0:4])
month = int(file.split('/')[-1].split('.')[1][4:6])
day = int(file.split('/')[-1].split('.')[1][6:8])
hour = int(file.split('/')[-1].split('.')[1][8:10])
dt = int(file.split('/')[-1].split('.')[3][1:])
timestamp_hycom = mdates.date2num(datetime(year,month,day,hour)) + dt/24
time_hycomf = mdates.num2date(timestamp_hycom)

delta_temp = temp_hycomf - temp_hycomi
#delta_temp[delta_temp<-0.5] = -0.5

kw = dict(levels = np.linspace(-2,2,21))
#fig,ax = plt.subplots(figsize=(5,5)) 
fig, ax = plt.subplots()
plt.contourf(lon_hycom[0,:]-360,lat_hycom[:,0],delta_temp,cmap=cmocean.cm.curl,**kw)
plt.colorbar()
plt.title('HWRF2020-HYCOM Delta SST \n'+\
          str(time_hycomf)[0:13]+'-'+ str(time_hycomi)[0:13],fontsize=14)
plt.axis('scaled')
rect = patches.Rectangle((xlim[0],ylim[0]),\
                 xlim[1]-xlim[0],ylim[1]-ylim[0],\
                 linewidth=1,edgecolor='k',facecolor='none')
ax.add_patch(rect)

plt.xlim([-80,-60])
plt.ylim([15,35])

plt.ylim([ylim[0],ylim[1]])
plt.xlim([xlim[0],xlim[1]])

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%% Surface temperature POM operational

folder_pom =folder_pom_oper
prefix = prefix_pom
N = 11
ylim = [24.6,26.1]
xlim = [-72.5,-71]

pom_ncfiles = sorted(glob.glob(os.path.join(folder_pom,prefix+'*.nc')))   
file = pom_ncfiles[N]
pom = xr.open_dataset(file)
time_pom = np.asarray(pom['time'])
temp_POM = np.asarray(pom['t'][0,0,:,:])
temp_POM[temp_POM==0] = np.nan
#temp_POM[temp_POM<27] = 27
#temp_POM[temp_POM>32] = 32

kw = dict(levels = np.linspace(28,30,41))
#kw = dict(levels = np.linspace(26,32,17))
plt.figure()
plt.contourf(lon_pom_oper[0,:],lat_pom_oper[:,0],temp_POM,cmap=cmocean.cm.thermal,**kw)
#plt.xlim([np.min(lon_forec_track_pom_oper)-0.5,np.max(lon_forec_track_pom_oper)+0.5])
#plt.ylim([np.min(lat_forec_track_pom_oper)-0.5,np.max(lat_forec_track_pom_oper)+0.5])
plt.colorbar()
plt.title('HWRF2019-POM SST \n on '+str(time_pom[0])[0:13],fontsize=14)
plt.axis('scaled')
plt.ylim([ylim[0],ylim[1]])
plt.xlim([xlim[0],xlim[1]])

#%% Surface temperature POM experimental

folder_pom =folder_pom_exp
prefix = prefix_pom
N = 11

pom_ncfiles = sorted(glob.glob(os.path.join(folder_pom,prefix+'*.nc')))   
file = pom_ncfiles[N]
pom = xr.open_dataset(file)
time_pom = np.asarray(pom['time'])
temp_POM = np.asarray(pom['t'][0,0,:,:])
temp_POM[temp_POM==0] = np.nan
#temp_POM[temp_POM<27] = 27
#temp_POM[temp_POM>32] = 32

kw = dict(levels = np.linspace(28,30,41))
#kw = dict(levels = np.linspace(26,32,17))
plt.figure()
plt.contourf(lon_pom_oper[0,:],lat_pom_oper[:,0],temp_POM,cmap=cmocean.cm.thermal,**kw)
#plt.xlim([np.min(lon_forec_track_pom_oper)-0.5,np.max(lon_forec_track_pom_oper)+0.5])
#plt.ylim([np.min(lat_forec_track_pom_oper)-0.5,np.max(lat_forec_track_pom_oper)+0.5])
plt.colorbar()
plt.title('HWRF2020-POM SST \n on '+str(time_pom[0])[0:13],fontsize=14)
plt.axis('scaled')
plt.ylim([ylim[0],ylim[1]])
plt.xlim([xlim[0],xlim[1]])

#%% Surface temperature HYCOM experimental

folder_hycom =folder_hycom_exp
prefix = prefix_hycom
N = 11
var = 'temp'

afiles = sorted(glob.glob(os.path.join(folder_hycom,prefix+'*.a')))    
file = afiles[N]

# Reading 3D variable from binary file 
var_hyc = readBinz(file[:-2],'3z',var)
temp_HYCOM = var_hyc[:,:,0]

    #Reading time stamp
year = int(file.split('/')[-1].split('.')[1][0:4])
month = int(file.split('/')[-1].split('.')[1][4:6])
day = int(file.split('/')[-1].split('.')[1][6:8])
hour = int(file.split('/')[-1].split('.')[1][8:10])
dt = int(file.split('/')[-1].split('.')[3][1:])
timestamp_hycom = mdates.date2num(datetime(year,month,day,hour)) + dt/24
time_hycom = mdates.num2date(timestamp_hycom)

kw = dict(levels = np.linspace(28,30,41))
plt.figure()
plt.contourf(lon_hycom[0,:]-360,lat_hycom[:,0],temp_HYCOM,cmap=cmocean.cm.thermal,**kw)
#plt.xlim([np.min(lon_forec_track_pom_oper[oktt])-0.5,np.max(lon_forec_track_pom_oper[oktt])+0.5])
#plt.ylim([np.min(lat_forec_track_pom_oper[oktt])-0.5,np.max(lat_forec_track_pom_oper[oktt])+0.5])
plt.colorbar()
plt.title('HWRF2020-HYCOM SST \n on '+str(time_hycom)[0:13],fontsize=14)
plt.axis('scaled')
plt.ylim([ylim[0],ylim[1]])
plt.xlim([xlim[0],xlim[1]])

#%% Surface salinity POM operational

folder_pom =folder_pom_oper
prefix = prefix_pom
N = 2
ylim = [24.6,26.1]
xlim = [-72.5,-71]

pom_ncfiles = sorted(glob.glob(os.path.join(folder_pom,prefix+'*.nc')))   
file = pom_ncfiles[N]
pom = xr.open_dataset(file)
time_pom = np.asarray(pom['time'])
salt_POM = np.asarray(pom['s'][0,0,:,:])
salt_POM[salt_POM==0] = np.nan
#temp_POM[temp_POM<27] = 27
#temp_POM[temp_POM>32] = 32

kw = dict(levels = np.linspace(33,36,31))
#kw = dict(levels = np.linspace(26,32,17))
plt.figure()
plt.contourf(lon_pom_oper[0,:],lat_pom_oper[:,0],salt_POM,cmap=cmocean.cm.haline,**kw)
#plt.xlim([np.min(lon_forec_track_pom_oper)-0.5,np.max(lon_forec_track_pom_oper)+0.5])
#plt.ylim([np.min(lat_forec_track_pom_oper)-0.5,np.max(lat_forec_track_pom_oper)+0.5])
plt.colorbar()
plt.title('HWRF2019-POM SSS \n on '+str(time_pom[0])[0:13],fontsize=14)
plt.axis('scaled')
plt.ylim([ylim[0],ylim[1]])
plt.xlim([xlim[0],xlim[1]])

#%% Surface temperature POM experimental

folder_pom =folder_pom_exp
prefix = prefix_pom
N = 2

pom_ncfiles = sorted(glob.glob(os.path.join(folder_pom,prefix+'*.nc')))   
file = pom_ncfiles[N]
pom = xr.open_dataset(file)
time_pom = np.asarray(pom['time'])
salt_POM = np.asarray(pom['s'][0,0,:,:])
#salt_POM[salt_POM==0] = np.nan
#temp_POM[temp_POM<27] = 27
#temp_POM[temp_POM>32] = 32

#kw = dict(levels = np.linspace(28,30,41))
#kw = dict(levels = np.linspace(26,32,17))
plt.figure()
plt.contourf(lon_pom_oper[0,:],lat_pom_oper[:,0],salt_POM,cmap=cmocean.cm.haline,**kw)
#plt.xlim([np.min(lon_forec_track_pom_oper)-0.5,np.max(lon_forec_track_pom_oper)+0.5])
#plt.ylim([np.min(lat_forec_track_pom_oper)-0.5,np.max(lat_forec_track_pom_oper)+0.5])
plt.colorbar()
plt.title('HWRF2020-POM SSS \n on '+str(time_pom[0])[0:13],fontsize=14)
plt.axis('scaled')
plt.ylim([ylim[0],ylim[1]])
plt.xlim([xlim[0],xlim[1]])

#%% Surface temperature HYCOM experimental

folder_hycom =folder_hycom_exp
prefix = prefix_hycom
N = 2
var = 'salinity'

afiles = sorted(glob.glob(os.path.join(folder_hycom,prefix+'*.a')))    
file = afiles[N]

# Reading 3D variable from binary file 
var_hyc = readBinz(file[:-2],'3z',var)
salt_HYCOM = var_hyc[:,:,0]

    #Reading time stamp
year = int(file.split('/')[-1].split('.')[1][0:4])
month = int(file.split('/')[-1].split('.')[1][4:6])
day = int(file.split('/')[-1].split('.')[1][6:8])
hour = int(file.split('/')[-1].split('.')[1][8:10])
dt = int(file.split('/')[-1].split('.')[3][1:])
timestamp_hycom = mdates.date2num(datetime(year,month,day,hour)) + dt/24
time_hycom = mdates.num2date(timestamp_hycom)

#kw = dict(levels = np.linspace(28,30,41))
plt.figure()
plt.contourf(lon_hycom[0,:]-360,lat_hycom[:,0],salt_HYCOM,cmap=cmocean.cm.haline,**kw)
#plt.xlim([np.min(lon_forec_track_pom_oper[oktt])-0.5,np.max(lon_forec_track_pom_oper[oktt])+0.5])
#plt.ylim([np.min(lat_forec_track_pom_oper[oktt])-0.5,np.max(lat_forec_track_pom_oper[oktt])+0.5])
plt.colorbar()
plt.title('HWRF2020-HYCOM SST \n on '+str(time_hycom)[0:13],fontsize=14)
plt.axis('scaled')
plt.ylim([ylim[0],ylim[1]])
plt.xlim([xlim[0],xlim[1]])

#%% Surface salinity POM operational

folder_pom =folder_pom_oper
prefix = prefix_pom
N = 11
ylim = [17,18]
xlim = [-65.5,-64.5]

pom_ncfiles = sorted(glob.glob(os.path.join(folder_pom,prefix+'*.nc')))   
file = pom_ncfiles[N]
pom = xr.open_dataset(file)
time_pom = np.asarray(pom['time'])
salt_POM = np.asarray(pom['s'][0,0,:,:])
salt_POM[salt_POM==0] = np.nan
#temp_POM[temp_POM<27] = 27
#temp_POM[temp_POM>32] = 32

kw = dict(levels = np.linspace(33,36,31))
#kw = dict(levels = np.linspace(26,32,17))
plt.figure()
plt.contourf(lon_pom_oper[0,:],lat_pom_oper[:,0],salt_POM,cmap=cmocean.cm.haline,**kw)
#plt.xlim([np.min(lon_forec_track_pom_oper)-0.5,np.max(lon_forec_track_pom_oper)+0.5])
#plt.ylim([np.min(lat_forec_track_pom_oper)-0.5,np.max(lat_forec_track_pom_oper)+0.5])
plt.colorbar()
plt.title('HWRF2019-POM SSS \n on '+str(time_pom[0])[0:13],fontsize=14)
plt.axis('scaled')
plt.ylim([ylim[0],ylim[1]])
plt.xlim([xlim[0],xlim[1]])

#%% Surface temperature POM experimental

folder_pom =folder_pom_exp
prefix = prefix_pom
N = 11
ylim = [17,18]
xlim = [-65.5,-64.5]

pom_ncfiles = sorted(glob.glob(os.path.join(folder_pom,prefix+'*.nc')))   
file = pom_ncfiles[N]
pom = xr.open_dataset(file)
time_pom = np.asarray(pom['time'])
salt_POM = np.asarray(pom['s'][0,0,:,:])
#salt_POM[salt_POM==0] = np.nan
#temp_POM[temp_POM<27] = 27
#temp_POM[temp_POM>32] = 32

kw = dict(levels = np.linspace(32.6,36,31))

#kw = dict(levels = np.linspace(28,30,41))
#kw = dict(levels = np.linspace(26,32,17))
plt.figure()
plt.contourf(lon_pom_oper[0,:],lat_pom_oper[:,0],salt_POM,cmap=cmocean.cm.haline,**kw)
#plt.xlim([np.min(lon_forec_track_pom_oper)-0.5,np.max(lon_forec_track_pom_oper)+0.5])
#plt.ylim([np.min(lat_forec_track_pom_oper)-0.5,np.max(lat_forec_track_pom_oper)+0.5])
plt.colorbar()
plt.title('HWRF2020-POM SSS \n on '+str(time_pom[0])[0:13],fontsize=14)
plt.axis('scaled')
plt.ylim([ylim[0],ylim[1]])
plt.xlim([xlim[0],xlim[1]])

#%% Surface temperature HYCOM experimental

folder_hycom =folder_hycom_exp
prefix = prefix_hycom
N = 11
var = 'salinity'

afiles = sorted(glob.glob(os.path.join(folder_hycom,prefix+'*.a')))    
file = afiles[N]

# Reading 3D variable from binary file 
var_hyc = readBinz(file[:-2],'3z',var)
salt_HYCOM = var_hyc[:,:,0]

    #Reading time stamp
year = int(file.split('/')[-1].split('.')[1][0:4])
month = int(file.split('/')[-1].split('.')[1][4:6])
day = int(file.split('/')[-1].split('.')[1][6:8])
hour = int(file.split('/')[-1].split('.')[1][8:10])
dt = int(file.split('/')[-1].split('.')[3][1:])
timestamp_hycom = mdates.date2num(datetime(year,month,day,hour)) + dt/24
time_hycom = mdates.num2date(timestamp_hycom)

#kw = dict(levels = np.linspace(28,30,41))
plt.figure()
plt.contourf(lon_hycom[0,:]-360,lat_hycom[:,0],salt_HYCOM,cmap=cmocean.cm.haline,**kw)
#plt.xlim([np.min(lon_forec_track_pom_oper[oktt])-0.5,np.max(lon_forec_track_pom_oper[oktt])+0.5])
#plt.ylim([np.min(lat_forec_track_pom_oper[oktt])-0.5,np.max(lat_forec_track_pom_oper[oktt])+0.5])
plt.colorbar()
plt.title('HWRF2020-HYCOM SST \n on '+str(time_hycom)[0:13],fontsize=14)
plt.axis('scaled')
plt.ylim([ylim[0],ylim[1]])
plt.xlim([xlim[0],xlim[1]])

