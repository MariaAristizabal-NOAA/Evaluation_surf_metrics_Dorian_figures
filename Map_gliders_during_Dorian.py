#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 15:18:14 2020

@author: root
"""

#%% User input

#lon_lim = [-100.0,-55.0]
#lat_lim = [10.0,45.0]

#lon_lim = [-80.0,-60.0]
#lat_lim = [15.0,35.0]

lat_lim = [10,50]
lon_lim = [-100,-50]

# limits for rectangle on map
lat_lim2 = [15,24]
lon_lim2 = [-70.5,-63] 


# Server erddap url IOOS glider dap
server = 'https://data.ioos.us/gliders/erddap'

# Time bounds
min_time = '2019/08/20/00'
max_time = '2019/09/07/00'

# Bathymetry file
bath_file = '/Users/aristizabal/Desktop/MARACOOS_project/Maria_scripts/nc_files/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'

# KMZ file
kmz_file = '/Users/aristizabal/Desktop/MARACOOS_project/Maria_scripts/KMZ_files/al052019_best_track-5.kmz'

# Hurricane category figures
ts_fig = '/Users/aristizabal/Desktop/MARACOOS_project/Maria_scripts/KMZ_files/ts_nhemi.png'

# Folder where to save figure
folder_fig = '/Users/aristizabal/Desktop/my_papers/Evaluation_surf_metrics_Dorian/figures/'

#%%

from matplotlib import pyplot as plt
import matplotlib.patches as patches
import numpy as np
from netCDF4 import Dataset
import cmocean
import os
import matplotlib.dates as mdates
from bs4 import BeautifulSoup
from zipfile import ZipFile
from erddapy import ERDDAP
import pandas as pd
from datetime import datetime

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

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

#%% Reading KMZ file

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
    print(s.get_text("coordinates"))
    lon_best_track[i] = float(s.get_text("coordinates").split('coordinates')[1].split(',')[0])
    lat_best_track[i] = float(s.get_text("coordinates").split('coordinates')[1].split(',')[1])
           
#  get time stamp
time_best_track = []
for i,s in enumerate(soup.find_all("atcfdtg")):
    tt = datetime.strptime(s.get_text(' '),'%Y%m%d%H')
    time_best_track.append(tt)
time_best_track = np.asarray(time_best_track)    

# get type 
wind_int = []
for i,s in enumerate(soup.find_all("intensitymph")):
    wind_int.append(s.get_text(' ')) 
wind_int = np.asarray(wind_int)
wind_int = wind_int.astype(float)  

wind_int_kt = []
for i,s in enumerate(soup.find_all("intensity")):
    wind_int_kt.append(s.get_text(' ')) 
wind_int_kt = np.asarray(wind_int_kt)
wind_int_kt = wind_int_kt.astype(float)
  
cat = []
for i,s in enumerate(soup.find_all("styleurl")):
    cat.append(s.get_text('#').split('#')[-1]) 
cat = np.asarray(cat)  

#%% Map of North Atlantic with glider tracks

lat_lim = [10,50]
lon_lim = [-100,-50]  

# Reading bathymetry data

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
fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(10, 5))
fig.suptitle('Active Glider Deployments '+ min_time[0:4]+'/'+min_time[5:7]+'/'+min_time[8:10]+\
             '-'+ max_time[0:4]+'/'+max_time[5:7]+'/'+max_time[8:10],fontsize=16)

ax1.contourf(bath_lonsub,bath_latsub,bath_elevsub,lev,cmap=cmocean.cm.topo)
ax1.set_xticks([])
ax1.set_yticks([])
ax1.plot(lon_best_track[5:-2],lat_best_track[5:-2],'o',color='firebrick',markersize=3,\
         markeredgecolor='k')
#ax1.set_aspect(1)
ax1.axis('scaled')
rect = patches.Rectangle((lon_lim2[0],lat_lim2[0]),\
                             lon_lim2[1]-lon_lim2[0],lat_lim2[1]-lat_lim2[0],\
                             linewidth=2,edgecolor='k',facecolor='none')
ax1.add_patch(rect)

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
        ax1.plot(long,latg,'.-',color='darkorange',markersize=1)

ax1.text(-99,47,'(a)',fontsize=16,color='white',weight='bold')   
ax1.axis('scaled')
ax1.axis([lon_lim[0],lon_lim[1],lat_lim[0],lat_lim[1]])

'''
folder = '/Users/aristizabal/Desktop/MARACOOS_project/Maria_scripts/Figures/Model_glider_comp/'
file = folder + 'Daily_map_North_Atlantic_gliders_in_DAC_' + min_time[0:10] + '_' + max_time[0:10] + '.png'
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)
'''

############## Map detail

# Time bounds
#min_time = '2019-08-24T00:00:00Z'
#max_time = '2019-09-07T00:00:00Z'

# Time bounds
min_time = '2019/08/28/00'
max_time = '2019/08/29/00'
tmin = datetime.strptime(min_time,'%Y/%m/%d/%H')
tmax = datetime.strptime(max_time,'%Y/%m/%d/%H')

lat_lim = lat_lim2
lon_lim = lon_lim2

# Reading bathymetry data
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

#lev = np.arange(-9000,9100,100)
#fig, ax2 = plt.subplots(figsize=(10, 5))
#plt.contour(bath_lonsub,bath_latsub,bath_elevsub,[0],colors='k')
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
ax2.legend(loc='upper right')  
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
for i in np.arange(17,24):
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

ax2.text(-70.4,23.4,'(b)',fontsize=16,color='white',weight='bold')
plt.subplots_adjust(left=0.05,right=0.95, wspace=-0.1, top=0.9)
file = folder_fig + 'Daily_map_North_Atlantic_gliders_in_DAC_detail' + '.png' #+ min_time[0:10] + '_' + max_time[0:10] + '.png'
#file = 'Daily_map_North_Atlantic_gliders_in_DAC_detail_' + str(tmin)[0:10] + '_' + str(tmax)[0:10] 
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)  