#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 16:06:49 2020

@author: root
"""

#%% User input

home_folder = '/home/'
#home_folder = '/Volumes'

#lon_lim = [-100.0,-55.0]
#lat_lim = [10.0,45.0]

lon_lim = [-85.0,-60.0]
lat_lim = [15.0,35.0]

# Server erddap url IOOS glider dap
server = 'https://data.ioos.us/gliders/erddap'

#gliders sg666, sg665, sg668, silbo
gdata_sg665 = 'http://gliders.ioos.us/thredds/dodsC/deployments/aoml/SG665-20190718T1155/SG665-20190718T1155.nc3.nc'
gdata_sg666 = 'http://gliders.ioos.us/thredds/dodsC/deployments/aoml/SG666-20190718T1206/SG666-20190718T1206.nc3.nc'
gdata_sg668 = 'http://gliders.ioos.us/thredds/dodsC/deployments/aoml/SG668-20190819T1217/SG668-20190819T1217.nc3.nc'

#gliders sg666, sg665, sg668, silbo
url_aoml = 'http://gliders.ioos.us/thredds/dodsC/deployments/aoml/'

gdata_sg665 = url_aoml+'SG665-20190718T1155/SG665-20190718T1155.nc3.nc'
gdata_sg666 = url_aoml+'SG666-20190718T1206/SG666-20190718T1206.nc3.nc'
gdata_sg668 = url_aoml+'SG668-20190819T1217/SG668-20190819T1217.nc3.nc'
gdata_sg664 = url_aoml+'SG664-20190716T1218/SG664-20190716T1218.nc3.nc'
gdata_sg663 = url_aoml+'SG663-20190716T1159/SG663-20190716T1159.nc3.nc'
gdata_sg667 = url_aoml+'SG667-20190815T1247/SG667-20190815T1247.nc3.nc'

gdata = gdata_sg665

# forecasting cycle to be used
cycle = '2019082800'

#Time window
#date_ini = '2019/08/28/00/00'
#date_end = '2019/09/02/06/00'

# Bathymetry file
bath_file = home_folder+'aristizabal/bathymetry_files/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'

# KMZ file
kmz_file_Dorian = home_folder+'aristizabal/KMZ_files/al052019_best_track-5.kmz'

# url for GOFS
url_GOFS = 'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_93.0/ts3z'

# Folder where to save figure
#folder_fig = home_folder+'aristizabal/Figures/'
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

#Dir_HMON_HYCOM = '/Volumes/aristizabal/ncep_model/HMON-HYCOM_Michael/'
Dir_HMON_HYCOM = home_folder+'aristizabal/HWRF_RTOFS_Michael_2018/HWRF-Hycom-WW3_exp_Michael/'

# RTOFS grid file name
hycom_grid_exp = Dir_HMON_HYCOM + 'hwrf_rtofs_hat10.basin.regional.grid'

# Dorian track files
hwrf_hycom_track_exp = folder_hycom_exp + 'dorian05l.' + cycle + '.trak.hwrf.atcfunix'

#%%

from matplotlib import pyplot as plt
import numpy as np
import xarray as xr
import netCDF4
from netCDF4 import Dataset
import cmocean
from datetime import datetime, timedelta
import matplotlib.dates as mdates
import sys
import seawater as sw
import os
import os.path
import glob
from bs4 import BeautifulSoup
from zipfile import ZipFile

#sys.path.append('/Users/aristizabal/Desktop/MARACOOS_project/Maria_scripts/All_code/Remote_repos/glider_model_comparisons_Python')
sys.path.append('/home/aristizabal/glider_model_comparisons_Python')
from read_glider_data import read_glider_data_thredds_server
#from process_glider_data import grid_glider_data_thredd

import sys
#sys.path.append('/Users/aristizabal/Desktop/MARACOOS_project/NCEP_scripts')
sys.path.append('/home/aristizabal/NCEP_scripts')
from utils4HYCOM import readBinz, readgrids

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

#%% Function Conversion from glider longitude and latitude to GOFS convention

def glider_coor_to_GOFS_coord(long,latg):

    target_lon = np.empty((len(long),))
    target_lon[:] = np.nan
    for i,ii in enumerate(long):
        if ii < 0:
            target_lon[i] = 360 + ii
        else:
            target_lon[i] = ii
    target_lat = latg

    return target_lon, target_lat

#%%  Function Conversion from GOFS convention to glider longitude and latitude

def GOFS_coor_to_glider_coord(lon_GOFS,lat_GOFS):

    lon_GOFSg = np.empty((len(lon_GOFS),))
    lon_GOFSg[:] = np.nan
    for i in range(len(lon_GOFS)):
        if lon_GOFS[i] > 180:
            lon_GOFSg[i] = lon_GOFS[i] - 360
        else:
            lon_GOFSg[i] = lon_GOFS[i]
    lat_GOFSg = lat_GOFS

    return lon_GOFSg, lat_GOFSg

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

#%% Function Getting glider transect from GOFS

def get_glider_transect_from_GOFS(depth_GOFS,oktime_GOFS,oklon_GOFS,oklat_GOFS):

    print('Getting glider transect from GOFS')
    target_temp_GOFS = np.empty((len(depth_GOFS),len(oktime_GOFS[0])))
    target_temp_GOFS[:] = np.nan
    target_salt_GOFS = np.empty((len(depth_GOFS),len(oktime_GOFS[0])))
    target_salt_GOFS[:] = np.nan
    for i in range(len(oktime_GOFS[0])):
        print(len(oktime_GOFS[0]),' ',i)
        target_temp_GOFS[:,i] = GOFS.variables['water_temp'][oktime_GOFS[0][i],:,oklat_GOFS[i],oklon_GOFS[i]]
        target_salt_GOFS[:,i] = GOFS.variables['salinity'][oktime_GOFS[0][i],:,oklat_GOFS[i],oklon_GOFS[i]]

    return target_temp_GOFS,target_salt_GOFS

#%%  Function Getting glider transect from POM

def get_glider_transect_from_POM(folder,prefix,zlev,zmatrix_pom,lon_pom,lat_pom,tstamp_glider,long,latg):

    ncfiles = sorted(glob.glob(os.path.join(folder,prefix+'*.nc')))

    target_temp_POM = np.empty((len(zlev),len(ncfiles)))
    target_temp_POM[:] = np.nan
    target_salt_POM = np.empty((len(zlev),len(ncfiles)))
    target_salt_POM[:] = np.nan
    target_dens_POM = np.empty((len(zlev),len(ncfiles)))
    target_dens_POM[:] = np.nan
    target_depth_POM = np.empty((len(zlev),len(ncfiles)))
    target_depth_POM[:] = np.nan
    time_POM = []

    for x,file in enumerate(ncfiles):
        print(x)
        pom = xr.open_dataset(file)

        tpom = pom['time'][:]
        timestamp_pom = mdates.date2num(tpom)[0]
        time_POM.append(mdates.num2date(timestamp_pom))

        # Interpolating latg and longlider into RTOFS grid
        sublonpom = np.interp(timestamp_pom,tstamp_glider,long)
        sublatpom = np.interp(timestamp_pom,tstamp_glider,latg)
        oklonpom = np.int(np.round(np.interp(sublonpom,lon_pom[0,:],np.arange(len(lon_pom[0,:])))))
        oklatpom = np.int(np.round(np.interp(sublatpom,lat_pom[:,0],np.arange(len(lat_pom[:,0])))))

        target_temp_POM[:,x] = np.asarray(pom['t'][0,:,oklatpom,oklonpom])
        target_salt_POM[:,x] = np.asarray(pom['s'][0,:,oklatpom,oklonpom])
        target_rho_pom = np.asarray(pom['rho'][0,:,oklatpom,oklonpom])
        target_dens_POM[:,x] = target_rho_pom * 1000 + 1000
        target_depth_POM[:,x] = zmatrix_pom[oklatpom,oklonpom,:].T

    target_temp_POM[target_temp_POM==0] = np.nan
    target_salt_POM[target_salt_POM==0] = np.nan
    target_dens_POM[target_dens_POM==1000.0] = np.nan

    return time_POM, target_temp_POM, target_salt_POM, target_dens_POM, target_depth_POM

#%%

def get_glider_transect_from_HYCOM(folder_hycom,prefix,nz,lon_hycom,lat_hycom,var,timestamp_glider,lon_glider,lat_glider):

    afiles = sorted(glob.glob(os.path.join(folder_hycom,prefix+'*.a')))

    target_var_hycom = np.empty((nz,len(afiles)))
    target_var_hycom[:] = np.nan
    time_hycom = []
    for x, file in enumerate(afiles):
        print(x)
        #lines=[line.rstrip() for line in open(file[:-2]+'.b')]

        #Reading time stamp
        year = int(file.split('/')[-1].split('.')[1][0:4])
        month = int(file.split('/')[-1].split('.')[1][4:6])
        day = int(file.split('/')[-1].split('.')[1][6:8])
        hour = int(file.split('/')[-1].split('.')[1][8:10])
        dt = int(file.split('/')[-1].split('.')[3][1:])
        timestamp_hycom = mdates.date2num(datetime(year,month,day,hour)) + dt/24
        time_hycom.append(mdates.num2date(timestamp_hycom))

        # Interpolating latg and longlider into HYCOM grid
        sublon_hycom = np.interp(timestamp_hycom,timestamp_glider,lon_glider)
        sublat_hycom = np.interp(timestamp_hycom,timestamp_glider,lat_glider)
        oklon_hycom = np.int(np.round(np.interp(sublon_hycom,lon_hycom[0,:],np.arange(len(lon_hycom[0,:])))))
        oklat_hycom = np.int(np.round(np.interp(sublat_hycom,lat_hycom[:,0],np.arange(len(lat_hycom[:,0])))))

        # Reading 3D variable from binary file
        var_hycom = readBinz(file[:-2],'3z',var)
        #ts=readBin(afile,'archive','temp')
        target_var_hycom[:,x] = var_hycom[oklat_hycom,oklon_hycom,:]

    time_hycom = np.asarray(time_hycom)
    #timestamp_hycom = mdates.date2num(time_hycom)

    return target_var_hycom, time_hycom

#%%  Calculation of mixed layer depth based on dt, Tmean: mean temp within the
# mixed layer and td: temp at 1 meter below the mixed layer

def MLD_temp_and_dens_criteria(dt,drho,time,depth,temp,salt,dens):

    MLD_temp_crit = np.empty(len(time))
    MLD_temp_crit[:] = np.nan
    Tmean_temp_crit = np.empty(len(time))
    Tmean_temp_crit[:] = np.nan
    Smean_temp_crit = np.empty(len(time))
    Smean_temp_crit[:] = np.nan
    Td_temp_crit = np.empty(len(time))
    Td_temp_crit[:] = np.nan
    MLD_dens_crit = np.empty(len(time))
    MLD_dens_crit[:] = np.nan
    Tmean_dens_crit = np.empty(len(time))
    Tmean_dens_crit[:] = np.nan
    Smean_dens_crit = np.empty(len(time))
    Smean_dens_crit[:] = np.nan
    Td_dens_crit = np.empty(len(time))
    Td_dens_crit[:] = np.nan
    for t,tt in enumerate(time):
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

#%% Function Ocean Heat Content

def OHC_surface(time,temp,depth,dens):
    cp = 3985 #Heat capacity in J/(kg K)

    OHC = np.empty((len(time)))
    OHC[:] = np.nan
    for t,tt in enumerate(time):
        ok26 = temp[:,t] >= 26
        if len(depth[ok26]) != 0:
            if np.nanmin(np.abs(depth[ok26]))>10:
                OHC[t] = np.nan
            else:
                rho0 = np.nanmean(dens[ok26,t])
                if depth.ndim == 1:
                    OHC[t] = np.abs(cp * rho0 * np.trapz(temp[ok26,t]-26,depth[ok26]))
                if depth.ndim == 2:
                    OHC[t] = np.abs(cp * rho0 * np.trapz(temp[ok26,t]-26,depth[ok26,t]))
        else:
            OHC[t] = np.nan

    return OHC

#%%
def depth_aver_top_100(depth,var):

    varmean100 = np.empty(var.shape[1])
    varmean100[:] = np.nan
    if depth.ndim == 1:
        okd = np.abs(depth) <= 100
        if len(depth[okd]) != 0:
            for t in np.arange(var.shape[1]):
                if len(np.where(np.isnan(var[okd,t]))[0])>10:
                    varmean100[t] = np.nan
                else:
                    varmean100[t] = np.nanmean(var[okd,t],0)
    else:
        for t in np.arange(depth.shape[1]):
            okd = np.abs(depth[:,t]) <= 100
            if len(depth[okd,t]) != 0:
                if len(np.where(np.isnan(var[okd,t]))[0])>10:
                    varmean100[t] = np.nan
                else:
                    varmean100[t] = np.nanmean(var[okd,t])
            else:
                varmean100[t] = np.nan

    return varmean100

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

#%% Calculate time series of potential Energy Anomaly over the top 100 m

def Potential_Energy_Anomaly(time,depth,density):
    g = 9.8 #m/s
    PEA = np.empty((len(time)))
    PEA[:] = np.nan
    for t,tstamp in enumerate(time):
        print(t)
        if np.ndim(depth) == 2:
            dindex = np.fliplr(np.where(np.asarray(np.abs(depth[:,t])) <= 100))[0]
        else:
            dindex = np.fliplr(np.where(np.asarray(np.abs(depth)) <= 100))[0]
        if len(dindex) == 0:
            PEA[t] = np.nan
        else:
            if np.ndim(depth) == 2:
                zz = np.asarray(np.abs(depth[dindex,t]))
            else:
                zz = np.asarray(np.abs(depth[dindex]))
            denss = np.asarray(density[dindex,t])
            ok = np.isfinite(denss)
            z = zz[ok]
            dens = denss[ok]
            if len(z)==0 or len(dens)==0 or np.min(zz) > 10:
                PEA[t] = np.nan
            else:
                if z[-1] - z[0] > 0:
                    # So PEA is < 0
                    #sign = -1
                    # Adding 0 to sigma integral is normalized
                    z = np.append(0,z)
                else:
                    # So PEA is < 0
                    #sign = 1
                    # Adding 0 to sigma integral is normalized
                    z = np.flipud(z)
                    z = np.append(0,z)
                    dens = np.flipud(dens)

                # adding density at depth = 0
                densit = np.interp(z,z[1:],dens)
                densit = np.flipud(densit)

                # defining sigma
                max_depth = np.nanmax(zz[ok])
                sigma = -1*z/max_depth
                sigma = np.flipud(sigma)

                rhomean = np.trapz(densit,sigma,axis=0)
                drho = rhomean-densit
                torque = drho * sigma
                PEA[t] = g* max_depth * np.trapz(torque,sigma,axis=0)
                #print(max_depth, ' ',PEA[t])

    return PEA

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

#%% Read GOFS 3.1 grid

#Time window
date_ini = cycle[0:4]+'/'+cycle[4:6]+'/'+cycle[6:8]+'/'+cycle[8:]+'/00/00'
tini = datetime.strptime(date_ini,'%Y/%m/%d/%H/%M/%S')
#tend = tini + timedelta(hours=120)
tend = tini + timedelta(hours=126)
date_end = tend.strftime('%Y/%m/%d/%H/%M/%S')

print('Retrieving coordinates from GOFS')
GOFS = xr.open_dataset(url_GOFS,decode_times=False)

tt_G = GOFS.time
t_G = netCDF4.num2date(tt_G[:],tt_G.units)

tmin = datetime.strptime(date_ini[0:-6],'%Y/%m/%d/%H')
tmax = datetime.strptime(date_end[0:-6],'%Y/%m/%d/%H')
oktime_GOFS = np.where(np.logical_and(t_G >= tmin, t_G <= tmax))
time_GOFS = np.asarray(t_G[oktime_GOFS])
#timestamp_GOFS = mdates.date2num(time_GOFS)
timestamp_GOFS = [mdates.date2num(datetime.strptime(str(t),t.format)) for t in time_GOFS]

lat_G = np.asarray(GOFS.lat[:])
lon_G = np.asarray(GOFS.lon[:])

# Conversion from glider longitude and latitude to GOFS convention
lon_limG, lat_limG = glider_coor_to_GOFS_coord(lon_lim,lat_lim)

oklat_GOFS = np.where(np.logical_and(lat_G >= lat_limG[0], lat_G <= lat_limG[1]))
oklon_GOFS = np.where(np.logical_and(lon_G >= lon_limG[0], lon_G <= lon_limG[1]))

lat_GOFS = lat_G[oklat_GOFS]
lon_GOFS = lon_G[oklon_GOFS]

depth_GOFS = np.asarray(GOFS.depth[:])

# Conversion from GOFS longitude and latitude to glider convention
lon_GOFSg, lat_GOFSg = GOFS_coor_to_glider_coord(lon_GOFS,lat_GOFS)

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

#%% Reading glider data

#Time window
date_ini = cycle[0:4]+'/'+cycle[4:6]+'/'+cycle[6:8]+'/'+cycle[8:]+'/00/00'
tini = datetime.strptime(date_ini,'%Y/%m/%d/%H/%M/%S')
tend = tini + timedelta(hours=120)
date_end = tend.strftime('%Y/%m/%d/%H/%M/%S')

url_glider = gdata

var_name = 'temperature'
scatter_plot = 'no'
kwargs = dict(date_ini=date_ini[0:-6],date_end=date_end[0:-6])

varg, timeg, latg, long, depthg, dataset_id = \
             read_glider_data_thredds_server(url_glider,var_name,scatter_plot)

tempg = varg

var_name = 'salinity'
varg, timeg, latg, long, depthg, dataset_id = \
             read_glider_data_thredds_server(url_glider,var_name,scatter_plot)

saltg = varg

var_name = 'density'
varg, timeg, latg, long, depthg, dataset_id = \
             read_glider_data_thredds_server(url_glider,var_name,scatter_plot)

densg = varg
depthg = depthg

#%% Grid glider variables according to depth

delta_z = 0.5
depthg_gridded,tempg_gridded,saltg_gridded,densg_gridded = \
varsg_gridded(depthg,timeg,tempg,saltg,densg,delta_z)

#%% Getting glider transect from GOFS

# Conversion from glider longitude and latitude to GOFS convention
target_lon, target_lat = glider_coor_to_GOFS_coord(long,latg)

# Changing times to timestamp
tstamp_glider = [mdates.date2num(timeg[i]) for i in np.arange(len(timeg))]
#tstamp_model = [mdates.date2num(time_GOFS[i]) for i in np.arange(len(time_GOFS))]
tstamp_model = [mdates.date2num(datetime.strptime(str(t),t.format)) for t in time_GOFS]

# interpolating glider lon and lat to lat and lon on model time
sublon_GOFS = np.interp(tstamp_model,tstamp_glider,target_lon)
sublat_GOFS = np.interp(tstamp_model,tstamp_glider,target_lat)

# Conversion from GOFS convention to glider longitude and latitude
sublon_GOFSg,sublat_GOFSg = GOFS_coor_to_glider_coord(sublon_GOFS,sublat_GOFS)

# getting the model grid positions for sublonm and sublatm
oklon_GOFS = np.round(np.interp(sublon_GOFS,lon_G,np.arange(len(lon_G)))).astype(int)
oklat_GOFS = np.round(np.interp(sublat_GOFS,lat_G,np.arange(len(lat_G)))).astype(int)

# Getting glider transect from model
target_temp_GOFS, target_salt_GOFS = \
                          get_glider_transect_from_GOFS(depth_GOFS,oktime_GOFS,oklon_GOFS,oklat_GOFS)

#%% Calculate density for GOFS

target_dens_GOFS = sw.dens(target_salt_GOFS,target_temp_GOFS,np.tile(depth_GOFS,(len(time_GOFS),1)).T)

#%% Retrieve glider transect from POM operational

tstamp_glider = [mdates.date2num(timeg[i]) for i in np.arange(len(timeg))]

folder_pom = folder_pom_oper
prefix = prefix_pom
zlev = zlev_pom_oper
zmatrix_pom = zmatrix_pom_oper
lon_pom = lon_pom_oper
lat_pom = lat_pom_oper
tstamp_glider = tstamp_glider
long = long
latg = latg

time_POM_oper, target_temp_POM_oper, target_salt_POM_oper, \
    target_dens_POM_oper, target_depth_POM_oper = \
    get_glider_transect_from_POM(folder_pom,prefix,zlev,zmatrix_pom,lon_pom,lat_pom,\
                                 tstamp_glider,long,latg)

timestamp_POM_oper = mdates.date2num(time_POM_oper)

#%% Retrieve glider transect from POM experimental

tstamp_glider = [mdates.date2num(timeg[i]) for i in np.arange(len(timeg))]

folder_pom = folder_pom_exp
prefix = prefix_pom
zlev = zlev_pom_exp
zmatrix_pom = zmatrix_pom_exp
lon_pom = lon_pom_exp
lat_pom = lat_pom_exp
tstamp_glider = tstamp_glider
long = long
latg = latg

time_POM_exp, target_temp_POM_exp, target_salt_POM_exp,\
    target_dens_POM_exp, target_depth_POM_exp = \
    get_glider_transect_from_POM(folder_pom,prefix,zlev,zmatrix_pom,lon_pom,lat_pom,\
                                 tstamp_glider,long,latg)
timestamp_POM_exp = mdates.date2num(time_POM_exp)

#%% Get glider transect from HYCOM

folder_hycom = folder_hycom_exp
prefix = prefix_hycom

# Changing times to timestamp
tstamp_glider = [mdates.date2num(timeg[i]) for i in np.arange(len(timeg))]

# Conversion from glider longitude and latitude to GOFS convention
target_lonG, target_latG = glider_coor_to_GOFS_coord(long,latg)

lon_glider = target_lonG
lat_glider = target_latG

var = 'temp'
target_temp_HYCOM_exp, time_HYCOM_exp = \
    get_glider_transect_from_HYCOM(folder_hycom,prefix,nz,\
    lon_hycom,lat_hycom,var,tstamp_glider,lon_glider,lat_glider)

var = 'salinity'
target_salt_HYCOM_exp, _ = \
    get_glider_transect_from_HYCOM(folder_hycom,prefix,nz,\
      lon_hycom,lat_hycom,var,tstamp_glider,lon_glider,lat_glider)

#%% Calculate density for HYCOM

target_dens_HYCOM_exp = sw.dens(target_salt_HYCOM_exp,target_temp_HYCOM_exp,np.tile(depth_HYCOM_exp,(len(time_HYCOM_exp),1)).T)

#%% Calculation of mixed layer depth based on temperature and density critria
# Tmean: mean temp within the mixed layer and
# td: temp at 1 meter below the mixed layer

dt = 0.2
drho = 0.125

# for glider data
MLD_temp_crit_glid, _, _, _, MLD_dens_crit_glid, Tmean_dens_crit_glid, Smean_dens_crit_glid, _ = \
MLD_temp_and_dens_criteria(dt,drho,timeg,depthg_gridded,tempg_gridded,saltg_gridded,densg_gridded)

# for GOFS 3.1 output
MLD_temp_crit_GOFS, _, _, _, MLD_dens_crit_GOFS, Tmean_dens_crit_GOFS, Smean_dens_crit_GOFS, _ = \
MLD_temp_and_dens_criteria(dt,drho,time_GOFS,depth_GOFS,target_temp_GOFS,target_salt_GOFS,target_dens_GOFS)

# for POM operational
MLD_temp_crit_POM_oper, _, _, _, MLD_dens_crit_POM_oper, Tmean_dens_crit_POM_oper, Smean_dens_crit_POM_oper, _ = \
MLD_temp_and_dens_criteria(dt,drho,timestamp_POM_oper,target_depth_POM_oper,target_temp_POM_oper,target_salt_POM_oper,target_dens_POM_oper)

# for POM experimental
MLD_temp_crit_POM_exp, _, _, _, MLD_dens_crit_POM_exp, Tmean_dens_crit_POM_exp, Smean_dens_crit_POM_exp, _ = \
MLD_temp_and_dens_criteria(dt,drho,timestamp_POM_exp,target_depth_POM_exp,target_temp_POM_exp,target_salt_POM_exp,target_dens_POM_exp)

# for HYCOM experimental
timestamp_HYCOM_exp = mdates.date2num(time_HYCOM_exp)
MLD_temp_crit_HYCOM_exp, _, _, _, MLD_dens_crit_HYCOM_exp, Tmean_dens_crit_HYCOM_exp, Smean_dens_crit_HYCOM_exp, _ = \
MLD_temp_and_dens_criteria(dt,drho,timestamp_HYCOM_exp,depth_HYCOM_exp,target_temp_HYCOM_exp,target_salt_HYCOM_exp,target_dens_HYCOM_exp)

#%% Figure transets

def figure_transect(var_name,var1,min_var1,max_var1,nlevels,var2,var3,\
                    time,tini,tend,depth,max_depth,color_map):

    if depth.ndim == 1:
        time_matrix = time
    else:
        time_matrix = np.tile(time,(depth.shape[0],1))

    kw = dict(levels = np.linspace(min_var1,max_var1,nlevels))

    # plot
    fig, ax = plt.subplots(figsize=(12, 2))
    cs = plt.contourf(time_matrix,depth,var1,cmap=color_map,**kw)
    plt.contour(time_matrix,depth,var1,[26],colors='k')
    plt.plot(time,var2,'-',label='MLD dt',color='indianred',linewidth=2 )
    plt.plot(time,var3,'-',label='MLD drho',color='darkgreen',marker='.',\
             markeredgecolor='k',linewidth=2 )

    cb = fig.colorbar(cs, orientation='vertical',pad=0.02)

    if var_name == 'temp':
        cb.set_label('($^oC$)',fontsize=14)

    ax.set_ylim(-200, 0)
    ax.set_ylabel('Depth (m)',fontsize=14)
    xticks = [tini+nday*timedelta(1) for nday in np.arange(14)]
    xticks = np.asarray(xticks)
    plt.xticks(xticks)
    xfmt = mdates.DateFormatter('%d-%b')
    ax.xaxis.set_major_formatter(xfmt)
    plt.legend()
    ax.set_xlim(tini,tend)

def figure_transect_subplot(ax,var_name,var1,min_var1,max_var1,nlevels,var2,var3,\
                    time,tini,tend,depth,max_depth,color_map):

    if depth.ndim == 1:
        time_matrix = time
    else:
        time_matrix = np.tile(time,(depth.shape[0],1))

    kw = dict(levels = np.linspace(min_var1,max_var1,nlevels))

    # plot

    cs = ax.contourf(time_matrix,depth,var1,cmap=color_map,**kw)
    ax.contour(time_matrix,depth,var1,[26],colors='k')
    ax.plot(time,var2,'-',label='MLD dt',color='indianred',linewidth=2 )
    ax.plot(time,var3,'-',label='MLD drho',color='darkgreen',marker='.',\
             markeredgecolor='k',linewidth=2 )

    cb = fig.colorbar(cs, orientation='vertical',pad=0.02)

    if var_name == 'temp':
        cb.set_label('($^oC$)',fontsize=12)

    ax.set_ylim(-200, 0)
    #ax.set_ylabel('Depth (m)',fontsize=14)
    xticks = [tini+nday*timedelta(1) for nday in np.arange(14)]
    xticks = np.asarray(xticks)
    ax.set_xticks(xticks)
    xfmt = mdates.DateFormatter('%d-%b')
    ax.xaxis.set_major_formatter(xfmt)
    #ax.legend()
    ax.set_xlim(tini,tend)

#%% All figures together
fig, ax = plt.subplots(figsize=(25, 25))
grid = plt.GridSpec(5, 2, wspace=0.1, hspace=0.3,left=0.05,right=0.95)

# Top 200 m glider temperature from 2019/08/28/00
color_map = cmocean.cm.thermal
var_name = 'temp'
var1 = tempg_gridded
var2 = -MLD_temp_crit_glid
var3 = -MLD_dens_crit_glid
time = timeg
depth = -depthg_gridded
max_depth = -200
min_var1 = 19
max_var1 = 31
nlevels = 13
tini = tini
tend = tend
ax1 = plt.subplot(grid[0, 0])
figure_transect_subplot(ax1,var_name,var1,min_var1,max_var1,nlevels,var2,var3,\
                time,tini,tend,depth,max_depth,color_map)
ax1.legend()
ax1.set_xticklabels([])
ax1.set_ylabel('Depth (m)',fontsize=14)
tDorian = np.tile(datetime(2019,8,29,0),len(np.arange(max_depth,0))) #ng665, ng666
ax1.plot(tDorian,np.arange(max_depth,0),'--k')
ax1.set_title('Along Track ' + 'Temperature' + ' Profile \n' + dataset_id,fontsize=14)
ax1.text(datetime(2019,8,28,1),10,'(a)',fontsize=14)

ax2 = plt.subplot(grid[0, 1])
color_map = cmocean.cm.haline
var_name = 'salt'
var1 = saltg_gridded
var2 = -MLD_temp_crit_glid
var3 = -MLD_dens_crit_glid
time = timeg
depth = -depthg_gridded
max_depth = -200
min_var1 = 35.5
max_var1 = 37.6
nlevels = 22
tini = tini
tend = tend
figure_transect_subplot(ax2,var_name,var1,min_var1,max_var1,nlevels,var2,var3,\
                time,tini,tend,depth,max_depth,color_map)
ax2.legend()
ax2.set_xticklabels([])
tDorian = np.tile(datetime(2019,8,29,0),len(np.arange(max_depth,0))) #ng665, ng666
ax2.plot(tDorian,np.arange(max_depth,0),'--k')
ax2.set_title('Along Track ' + 'Salinity' + ' Profile \n' +dataset_id,fontsize=14)
ax2.text(datetime(2019,8,28,1),10,'(f)',fontsize=14)

ax3 = plt.subplot(grid[1, 0])
color_map = cmocean.cm.thermal
var_name = 'temp'
var1 = target_temp_GOFS
var2 = -MLD_temp_crit_GOFS
var3 = -MLD_dens_crit_GOFS
time = timestamp_GOFS
depth = -depth_GOFS
max_depth = -200
min_var1 = 19
max_var1 = 31
nlevels = 13
tini = tini
tend = tend
figure_transect_subplot(ax3,var_name,var1,min_var1,max_var1,nlevels,var2,var3,\
                time,tini,tend,depth,max_depth,color_map)
ax3.set_xticklabels([])
ax3.set_ylabel('Depth (m)',fontsize=14)
tDorian = np.tile(datetime(2019,8,29,0),len(np.arange(max_depth,0))) #ng665, ng666
ax3.plot(tDorian,np.arange(max_depth,0),'--k')
ax3.set_title('GOFS 3.1',fontsize=14)
ax3.text(datetime(2019,8,28,1),10,'(b)',fontsize=14)

ax4 = plt.subplot(grid[1, 1])
color_map = cmocean.cm.haline
var_name = 'salt'
var1 = target_salt_GOFS
var2 = -MLD_temp_crit_GOFS
var3 = -MLD_dens_crit_GOFS
time = timestamp_GOFS
depth = -depth_GOFS
max_depth = -200
min_var1 = 35.5
max_var1 = 37.6
nlevels = 22
tini = tini
tend = tend
figure_transect_subplot(ax4,var_name,var1,min_var1,max_var1,nlevels,var2,var3,\
                time,tini,tend,depth,max_depth,color_map)
ax4.set_xticklabels([])
tDorian = np.tile(datetime(2019,8,29,0),len(np.arange(max_depth,0))) #ng665, ng666
ax4.plot(tDorian,np.arange(max_depth,0),'--k')
ax4.set_title('GOFS 3.1',fontsize=14)
ax4.text(datetime(2019,8,28,1),10,'(g)',fontsize=14)

ax5 = plt.subplot(grid[2, 0])
color_map = cmocean.cm.thermal
var_name = 'temp'
var1 = target_temp_POM_oper
var2 = MLD_temp_crit_POM_oper
var3 = MLD_dens_crit_POM_oper
time = mdates.date2num(time_POM_oper)
depth = target_depth_POM_oper
max_depth = -200
min_var1 = 19
max_var1 = 31
nlevels = 13
tini = tini
tend = tend
figure_transect_subplot(ax5,var_name,var1,min_var1,max_var1,nlevels,var2,var3,\
                time,tini,tend,depth,max_depth,color_map)
ax5.set_xticklabels([])
ax5.set_ylabel('Depth (m)',fontsize=14)
tDorian = np.tile(datetime(2019,8,29,0),len(np.arange(max_depth,0))) #ng665, ng666
ax5.plot(tDorian,np.arange(max_depth,0),'--k')
ax5.set_title('HWRF2019-POM (IC clim.)',fontsize=14)
ax5.text(datetime(2019,8,28,1),10,'(c)',fontsize=14)

ax6 = plt.subplot(grid[2, 1])
color_map = cmocean.cm.haline
var_name = 'salt'
var1 = target_salt_POM_oper
var2 = MLD_temp_crit_POM_oper
var3 = MLD_dens_crit_POM_oper
time = mdates.date2num(time_POM_oper)
depth = target_depth_POM_oper
max_depth = -200
min_var1 = 35.5
max_var1 = 37.6
nlevels = 22
tini = tini
tend = tend
figure_transect_subplot(ax6,var_name,var1,min_var1,max_var1,nlevels,var2,var3,\
                time,tini,tend,depth,max_depth,color_map)
ax6.set_xticklabels([])
tDorian = np.tile(datetime(2019,8,29,0),len(np.arange(max_depth,0))) #ng665, ng666
ax6.plot(tDorian,np.arange(max_depth,0),'--k')
ax6.set_title('HWRF2019-POM (IC clim.)',fontsize=14)
ax6.text(datetime(2019,8,28,1),10,'(h)',fontsize=14)

ax7 = plt.subplot(grid[3, 0])
color_map = cmocean.cm.thermal
var_name = 'temp'
var1 = target_temp_POM_exp
var2 = MLD_temp_crit_POM_exp
var3 = MLD_dens_crit_POM_exp
time = mdates.date2num(time_POM_exp)
depth = target_depth_POM_exp
max_depth = -200
min_var1 = 19
max_var1 = 31
nlevels = 13
tini = tini
tend = tend
figure_transect_subplot(ax7,var_name,var1,min_var1,max_var1,nlevels,var2,var3,\
                time,tini,tend,depth,max_depth,color_map)
ax7.set_xticklabels([])
ax7.set_ylabel('Depth (m)',fontsize=14)
tDorian = np.tile(datetime(2019,8,29,0),len(np.arange(max_depth,0))) #ng665, ng666
ax7.plot(tDorian,np.arange(max_depth,0),'--k')
ax7.set_title('HWRF2020-POM (IC RTOFS)',fontsize=14)
ax7.text(datetime(2019,8,28,1),10,'(d)',fontsize=14)

ax8 = plt.subplot(grid[3, 1])
color_map = cmocean.cm.haline
var_name = 'salt'
var1 = target_salt_POM_exp
var2 = MLD_temp_crit_POM_exp
var3 = MLD_dens_crit_POM_exp
time = mdates.date2num(time_POM_exp)
depth = target_depth_POM_exp
max_depth = -200
min_var1 = 35.5
max_var1 = 37.6
nlevels = 22
tini = tini
tend = tend
figure_transect_subplot(ax8,var_name,var1,min_var1,max_var1,nlevels,var2,var3,\
                time,tini,tend,depth,max_depth,color_map)
ax8.set_xticklabels([])
tDorian = np.tile(datetime(2019,8,29,0),len(np.arange(max_depth,0))) #ng665, ng666
ax8.plot(tDorian,np.arange(max_depth,0),'--k')
ax8.set_title('HWRF2020-POM (IC RTOFS)',fontsize=14)
ax8.text(datetime(2019,8,28,1),10,'(i)',fontsize=14)

ax9 = plt.subplot(grid[4, 0])
color_map = cmocean.cm.thermal
var_name = 'temp'
var1 = target_temp_HYCOM_exp
var2 = -MLD_temp_crit_HYCOM_exp
var3 = -MLD_dens_crit_HYCOM_exp
time = mdates.date2num(time_HYCOM_exp)
depth = -depth_HYCOM_exp
max_depth = -200
min_var1 = 19
max_var1 = 31
nlevels = 13
tini = tini
tend = tend
figure_transect_subplot(ax9,var_name,var1,min_var1,max_var1,nlevels,var2,var3,\
                time,tini,tend,depth,max_depth,color_map)
ax9.set_ylabel('Depth (m)',fontsize=14)
tDorian = np.tile(datetime(2019,8,29,0),len(np.arange(max_depth,0))) #ng665, ng666
ax9.plot(tDorian,np.arange(max_depth,0),'--k')
ax9.set_title('HWRF2020-HYCOM (IC RTOFS)',fontsize=14)
ax9.text(datetime(2019,8,28,1),10,'(e)',fontsize=14)

ax10 = plt.subplot(grid[4, 1])
color_map = cmocean.cm.haline
var_name = 'salt'
var1 = target_salt_HYCOM_exp
var2 = -MLD_temp_crit_HYCOM_exp
var3 = -MLD_dens_crit_HYCOM_exp
time = mdates.date2num(time_HYCOM_exp)
depth = -depth_HYCOM_exp
max_depth = -200
min_var1 = 35.5
max_var1 = 37.6
nlevels = 22
tini = tini
tend = tend
figure_transect_subplot(ax10,var_name,var1,min_var1,max_var1,nlevels,var2,var3,\
                time,tini,tend,depth,max_depth,color_map)
tDorian = np.tile(datetime(2019,8,29,0),len(np.arange(max_depth,0))) #ng665, ng666
ax10.plot(tDorian,np.arange(max_depth,0),'--k')
ax10.set_title('HWRF2020-HYCOM (IC RTOFS)',fontsize=14)
ax10.text(datetime(2019,8,28,1),10,'(j)',fontsize=14)

file = folder_fig + 'temperature_salinity_transects3.png'
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1,dpi=350)

##################################################################################
#%% Top 200 m glider temperature from 2019/08/28/00
color_map = cmocean.cm.thermal
var_name = 'temp'
var1 = tempg_gridded
var2 = -MLD_temp_crit_glid
var3 = -MLD_dens_crit_glid
time = timeg
depth = -depthg_gridded
max_depth = -200
min_var1 = 19
max_var1 = 31
nlevels = 13
tini = tini
tend = tend

figure_transect(var_name,var1,min_var1,max_var1,nlevels,var2,var3,\
                time,tini,tend,depth,max_depth,color_map)

#tDorian = np.tile(datetime(2019,8,28,18),len(np.arange(max_depth,0))) #ng668
tDorian = np.tile(datetime(2019,8,29,0),len(np.arange(max_depth,0))) #ng665, ng666
plt.plot(tDorian,np.arange(max_depth,0),'--k')
plt.title('Along Track ' + 'Temperature' + ' Profile ' + dataset_id,fontsize=14)
file = folder_fig + 'along_track_temp_top200_' + dataset_id
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)

#%% Top 200 m glider salinity from 2019/08/28/00

color_map = cmocean.cm.haline

var_name = 'salt'
var1 = saltg_gridded
var2 = -MLD_temp_crit_glid
var3 = -MLD_dens_crit_glid
time = timeg
depth = -depthg_gridded
max_depth = -200
min_var1 = 35.5
max_var1 = 37.6
nlevels = 22
tini = tini
tend = tend

figure_transect(var_name,var1,min_var1,max_var1,nlevels,var2,var3,\
                time,tini,tend,depth,max_depth,color_map)

#tDorian = np.tile(datetime(2019,8,28,18),len(np.arange(max_depth,0))) #ng668
tDorian = np.tile(datetime(2019,8,29,0),len(np.arange(max_depth,0))) #ng665, ng666
plt.plot(tDorian,np.arange(max_depth,0),'--k')
plt.title('Along Track ' + 'Salinity' + ' Profile ' + dataset_id,fontsize=14)
file = folder_fig + 'along_track_salt_top200_' + dataset_id
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)

#%% Top 200 m GOFS 3.1 temperature from 2019/08/28/00

color_map = cmocean.cm.thermal

var_name = 'temp'
var1 = target_temp_GOFS
var2 = -MLD_temp_crit_GOFS
var3 = -MLD_dens_crit_GOFS
time = mdates.date2num(time_GOFS)
depth = -depth_GOFS
max_depth = -200
min_var1 = 19
max_var1 = 31
nlevels = 13
tini = tini
tend = tend

figure_transect(var_name,var1,min_var1,max_var1,nlevels,var2,var3,\
                time,tini,tend,depth,max_depth,color_map)

#tDorian = np.tile(datetime(2019,8,28,18),len(np.arange(max_depth,0))) #ng668
tDorian = np.tile(datetime(2019,8,29,0),len(np.arange(max_depth,0))) #ng665, ng666
plt.plot(tDorian,np.arange(max_depth,0),'--k')
plt.title('Along Track ' + 'Temperature' + ' Profile ' + 'GOFS 3.1',fontsize=14)
file = folder_fig + 'along_track_temp_top200_GOFS_' + dataset_id
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)

#%% Top 200 m GOFS 3.1 salinity from 2019/08/28/00

color_map = cmocean.cm.haline

var_name = 'salt'
var1 = target_salt_GOFS
var2 = -MLD_temp_crit_GOFS
var3 = -MLD_dens_crit_GOFS
time = mdates.date2num(time_GOFS)
depth = -depth_GOFS
max_depth = -200
min_var1 = 35.5
max_var1 = 37.6
nlevels = 22
tini = tini
tend = tend

figure_transect(var_name,var1,min_var1,max_var1,nlevels,var2,var3,\
                time,tini,tend,depth,max_depth,color_map)

#tDorian = np.tile(datetime(2019,8,28,18),len(np.arange(max_depth,0))) #ng668
tDorian = np.tile(datetime(2019,8,29,0),len(np.arange(max_depth,0))) #ng665, ng666
plt.plot(tDorian,np.arange(max_depth,0),'--k')
plt.title('Along Track ' + 'Salinity' + ' Profile ' + 'GOFS 3.1',fontsize=14)
file = folder_fig + 'along_track_salt_top200_GOFS_' + dataset_id
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)

#%% Top 200 m POM operational temperature from 2019/08/28/00

color_map = cmocean.cm.thermal

var_name = 'temp'
var1 = target_temp_POM_oper
var2 = MLD_temp_crit_POM_oper
var3 = MLD_dens_crit_POM_oper
time = mdates.date2num(time_POM_oper)
depth = target_depth_POM_oper
max_depth = -200
min_var1 = 19
max_var1 = 31
nlevels = 13
tini = tini
tend = tend

figure_transect(var_name,var1,min_var1,max_var1,nlevels,var2,var3,\
                time,tini,tend,depth,max_depth,color_map)

#tDorian = np.tile(datetime(2019,8,28,18),len(np.arange(max_depth,0))) #ng668
tDorian = np.tile(datetime(2019,8,29,0),len(np.arange(max_depth,0))) #ng665, ng666
plt.plot(tDorian,np.arange(max_depth,0),'--k')
plt.title('Along Track Temperature Profile HWRF2019-POM (IC clim.)',fontsize=14)
file = folder_fig + 'along_track_temp_top200_POM_oper_' + dataset_id
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)

#%% Top 200 m POM operational salinity from 2019/08/28/00

color_map = cmocean.cm.haline

var_name = 'salt'
var1 = target_salt_POM_oper
var2 = MLD_temp_crit_POM_oper
var3 = MLD_dens_crit_POM_oper
time = mdates.date2num(time_POM_oper)
depth = target_depth_POM_oper
max_depth = -200
min_var1 = 35.5
max_var1 = 37.6
nlevels = 22
tini = tini
tend = tend

figure_transect(var_name,var1,min_var1,max_var1,nlevels,var2,var3,\
                time,tini,tend,depth,max_depth,color_map)

#tDorian = np.tile(datetime(2019,8,28,18),len(np.arange(max_depth,0))) #ng668
tDorian = np.tile(datetime(2019,8,29,0),len(np.arange(max_depth,0))) #ng665, ng666
plt.plot(tDorian,np.arange(max_depth,0),'--k')
plt.title('Along Track Salinity Profile HWRF2019-POM (IC clim.)',fontsize=14)
file = folder_fig + 'along_track_salt_top200_POM_oper_' + dataset_id
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)

#%% Top 200 m POM experimental temperature from 2019/08/28/00

color_map = cmocean.cm.thermal

var_name = 'temp'
var1 = target_temp_POM_exp
var2 = MLD_temp_crit_POM_exp
var3 = MLD_dens_crit_POM_exp
time = mdates.date2num(time_POM_exp)
depth = target_depth_POM_exp
max_depth = -200
min_var1 = 19
max_var1 = 31
nlevels = 13
tini = tini
tend = tend

figure_transect(var_name,var1,min_var1,max_var1,nlevels,var2,var3,\
                time,tini,tend,depth,max_depth,color_map)

#tDorian = np.tile(datetime(2019,8,28,18),len(np.arange(max_depth,0))) #ng668
tDorian = np.tile(datetime(2019,8,29,0),len(np.arange(max_depth,0))) #ng665, ng666
plt.plot(tDorian,np.arange(max_depth,0),'--k')
plt.plot(tDorian,np.arange(max_depth,0),'--k')
plt.title('Along Track Salinity Profile HWRF2020-POM (IC RTOFS)',fontsize=14)
file = folder_fig + 'along_track_salt_top200_POM_exp_' + dataset_id
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)

#%% Top 200 m POM experimental salinity from 2019/08/28/00

color_map = cmocean.cm.haline

var_name = 'salt'
var1 = target_salt_POM_exp
var2 = MLD_temp_crit_POM_exp
var3 = MLD_dens_crit_POM_exp
time = mdates.date2num(time_POM_exp)
depth = target_depth_POM_exp
max_depth = -200
min_var1 = 35.5
max_var1 = 37.6
nlevels = 22
tini = tini
tend = tend

figure_transect(var_name,var1,min_var1,max_var1,nlevels,var2,var3,\
                time,tini,tend,depth,max_depth,color_map)

#tDorian = np.tile(datetime(2019,8,28,18),len(np.arange(max_depth,0))) #ng668
tDorian = np.tile(datetime(2019,8,29,0),len(np.arange(max_depth,0))) #ng665, ng666
plt.plot(tDorian,np.arange(max_depth,0),'--k')
plt.plot(tDorian,np.arange(max_depth,0),'--k')
plt.title('Along Track Salinity Profile HWRF2020-POM (IC RTOFS)',fontsize=14)
file = folder_fig + 'along_track_salt_top200_POM_exp_' + dataset_id
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)

#%% Top 200 m HYCOM experimental temperature from 2019/08/28/00

color_map = cmocean.cm.thermal

var_name = 'temp'
var1 = target_temp_HYCOM_exp
var2 = -MLD_temp_crit_HYCOM_exp
var3 = -MLD_dens_crit_HYCOM_exp
time = mdates.date2num(time_HYCOM_exp)
depth = -depth_HYCOM_exp
max_depth = -200
min_var1 = 19
max_var1 = 31
nlevels = 13
tini = tini
tend = tend

figure_transect(var_name,var1,min_var1,max_var1,nlevels,var2,var3,\
                time,tini,tend,depth,max_depth,color_map)

#tDorian = np.tile(datetime(2019,8,28,18),len(np.arange(max_depth,0))) #ng668
tDorian = np.tile(datetime(2019,8,29,0),len(np.arange(max_depth,0))) #ng665, ng666
plt.plot(tDorian,np.arange(max_depth,0),'--k')
plt.plot(tDorian,np.arange(max_depth,0),'--k')
plt.title('Along Track Temperature Profile HWRF2020-HYCOM (IC RTOFS)',fontsize=14)
file = folder_fig + 'along_track_temp_top200_HYCOM_exp_' + dataset_id
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)

#%% Top 200 m HYCOM experimental salinity from 2019/08/28/00

color_map = cmocean.cm.haline

var_name = 'salt'
var1 = target_salt_HYCOM_exp
var2 = -MLD_temp_crit_HYCOM_exp
var3 = -MLD_dens_crit_HYCOM_exp
time = mdates.date2num(time_HYCOM_exp)
depth = -depth_HYCOM_exp
max_depth = -200
min_var1 = 35.5
max_var1 = 37.6
nlevels = 22
tini = tini
tend = tend

figure_transect(var_name,var1,min_var1,max_var1,nlevels,var2,var3,\
                time,tini,tend,depth,max_depth,color_map)

#tDorian = np.tile(datetime(2019,8,28,18),len(np.arange(max_depth,0))) #ng668
tDorian = np.tile(datetime(2019,8,29,0),len(np.arange(max_depth,0))) #ng665, ng666
plt.plot(tDorian,np.arange(max_depth,0),'--k')
plt.plot(tDorian,np.arange(max_depth,0),'--k')
plt.title('Along Track Salinity Profile HWRF2020-HYCOM (IC RTOFS)',fontsize=14)
file = folder_fig + 'along_track_salt_top200_HYCOM_exp_' + dataset_id
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1)
