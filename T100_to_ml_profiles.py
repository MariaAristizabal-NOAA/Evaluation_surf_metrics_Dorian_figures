#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 10:25:25 2020

@author: root
"""

#%% User input

home_folder = '/home/'
#home_folder = 'Volumes/'

lon_lim = [-80.0,-60.0]
lat_lim = [15.0,35.0]

# Server erddap url IOOS glider dap
server = 'https://data.ioos.us/gliders/erddap'

#gliders sg666, sg665, sg668, silbo
url_aoml = 'http://gliders.ioos.us/thredds/dodsC/deployments/aoml/'
#url_RU = 'http://gliders.ioos.us/thredds/dodsC/deployments/rutgers/'
gdata = [url_aoml+'SG665-20190718T1155/SG665-20190718T1155.nc3.nc',\
         url_aoml+'SG666-20190718T1206/SG666-20190718T1206.nc3.nc',\
         url_aoml+'SG668-20190819T1217/SG668-20190819T1217.nc3.nc',\
         url_aoml+'SG664-20190716T1218/SG664-20190716T1218.nc3.nc',\
         url_aoml+'SG663-20190716T1159/SG663-20190716T1159.nc3.nc',\
         url_aoml+'SG667-20190815T1247/SG667-20190815T1247.nc3.nc']
         
#url_RU+'silbo-20190717T1917/silbo-20190717T1917.nc3.nc']

cycle = '2019082800'    
    
#Time window
#date_ini = '2019/08/28/00/00'
#date_end = '2019/09/02/00/06'

# url for GOFS 3.1
url_GOFS = 'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_93.0/ts3z'

# figures
#folder_fig = '/Users/aristizabal/Desktop/my_papers/Evaluation_surf_metrics_Dorian/figures/'
folder_fig = '/home/aristizabal/Figures/'

# folder nc files POM
folder_pom19 =  home_folder+'aristizabal/HWRF2019_POM_Dorian/'
folder_pom20 =  home_folder+'aristizabal/HWRF2020_POM_Dorian/'

# folde HWRF2020_HYCOM
folder_hycom20 = home_folder+'aristizabal/HWRF2020_HYCOM_Dorian/'

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


##################
# folder ab files HYCOM
folder_hycom_exp = folder_hycom20 + 'HWRF2020_HYCOM_dorian05l.' + cycle + '_hycom_files_exp/'
prefix_hycom = 'dorian05l.' + cycle + '.hwrf_rtofs_hat10_3z'

#Dir_HMON_HYCOM = '/Volumes/aristizabal/ncep_model/HMON-HYCOM_Michael/'
Dir_HMON_HYCOM = home_folder+'aristizabal/ncep_model/HWRF-Hycom-WW3_exp_Michael/'
# RTOFS grid file name
hycom_grid_exp = Dir_HMON_HYCOM + 'hwrf_rtofs_hat10.basin.regional.grid'

# Dorian track files
hwrf_hycom_track_exp = folder_hycom_exp + 'dorian05l.' + cycle + '.trak.hwrf.atcfunix'

#%%
#import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import xarray as xr
import netCDF4
from datetime import datetime
from datetime import timedelta
import matplotlib.dates as mdates
import sys
import seawater as sw
import os
import os.path
import glob 

sys.path.append('/home/aristizabal/glider_model_comparisons_Python')
#sys.path.append('/Users/aristizabal/Desktop/MARACOOS_project/Maria_scripts/All_code/Remote_repos/glider_model_comparisons_Python')
from read_glider_data import read_glider_data_thredds_server
#from process_glider_data import grid_glider_data_thredd

import sys
#sys.path.append('/Users/aristizabal/Desktop/MARACOOS_project/NCEP_scripts')
sys.path.append('/home/aristizabal/NCEP_scripts')
from utils4HYCOM import readBinz
from utils4HYCOM import readgrids

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

#%% Function Getting glider transect from GOFS
    
def get_glider_transect_from_GOFS(GOFS,depth_GOFS,oktime_GOFS,oklat_GOFS,oklon_GOFS):
    
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
    
def get_glider_transect_from_POM(folder_pom,prefix,zlev,zmatrix_pom,lon_pom,lat_pom,tstamp_glider,long,latg):

    ncfiles = sorted(glob.glob(os.path.join(folder_pom,prefix+'*.nc')))

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
        year = int(file.split('.')[1][0:4])
        month = int(file.split('.')[1][4:6])
        day = int(file.split('.')[1][6:8])
        hour = int(file.split('.')[1][8:10])
        dt = int(file.split('.')[-2][1:])
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
            if np.nanmin(depth[ok26])>10:
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
def interp_datasets_in_z(temp_orig,salt_orig,depth_orig,depth_target,indexes_time,dim_vars):

    temp_interp = np.empty((dim_vars[0],dim_vars[1]))
    temp_interp[:] = np.nan
    salt_interp = np.empty((dim_vars[0],dim_vars[1]))
    salt_interp[:] = np.nan
    for i in np.arange(len(indexes_time)):
        pos = np.argsort(depth_orig[:,indexes_time[i]])
        if depth_target.ndim == 1:
            temp_interp[:,i] = np.interp(depth_target,depth_orig[pos,indexes_time[i]],temp_orig[pos,indexes_time[i]])
            salt_interp[:,i] = np.interp(depth_target,depth_orig[pos,indexes_time[i]],salt_orig[pos,indexes_time[i]])
        if depth_target.ndim == 2:
            temp_interp[:,i] = np.interp(depth_target[:,i],depth_orig[pos,indexes_time[i]],temp_orig[pos,indexes_time[i]])
            salt_interp[:,i] = np.interp(depth_target[:,i],depth_orig[pos,indexes_time[i]],salt_orig[pos,indexes_time[i]])

    return temp_interp, salt_interp

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

#%%    
def depth_aver_temp_from_100_to_base_mixed_layer(dt,drho,time,depth,temp,salt,dens):

    _, _, _, _, MLD_dens_crit_glid,_,_,_ = \
        MLD_temp_and_dens_criteria(dt,drho,time,depth,temp,salt,dens)    
        
    tempmean100_to_ml = np.empty(temp.shape[1])
    tempmean100_to_ml[:] = np.nan
    if depth.ndim == 1:
       for t in np.arange(temp.shape[1]):
           okd = np.logical_and(np.abs(depth) <= 100,np.abs(depth) >= np.abs(MLD_dens_crit_glid[t]))
           if len(depth[okd]) != 0:
              tempmean100_to_ml[t] = np.nanmean(temp[okd,t],0)
           else:
               tempmean100_to_ml[t] = np.nan
        
    else:
        for t in np.arange(temp.shape[1]):
            okd = np.logical_and(np.abs(depth[:,t]) <= 100,np.abs(depth[:,t]) >= np.abs(MLD_dens_crit_glid[t]))
            if len(depth[okd,t]) != 0:
                tempmean100_to_ml[t] = np.nanmean(temp[okd,t])
            else:
                tempmean100_to_ml[t] = np.nan
    
    return tempmean100_to_ml 

#%% Taylor Diagram

def taylor_template(angle_lim,std_lim):
    
    import mpl_toolkits.axisartist.floating_axes as floating_axes
    from matplotlib.projections import PolarAxes
    from mpl_toolkits.axisartist.grid_finder import (FixedLocator,
                                                 DictFormatter)

    fig = plt.figure()
    tr = PolarAxes.PolarTransform()
    
    min_corr = np.round(np.cos(angle_lim),1)
    CCgrid= np.concatenate((np.arange(min_corr*10,10,2.0)/10.,[0.9,0.95,0.99]))
    CCpolar=np.arccos(CCgrid)
    gf=FixedLocator(CCpolar)
    tf=DictFormatter(dict(zip(CCpolar, map(str,CCgrid))))
    
    STDgrid=np.arange(0,std_lim,.5)
    gfs=FixedLocator(STDgrid)
    tfs=DictFormatter(dict(zip(STDgrid, map(str,STDgrid))))
    
    ra0, ra1 =0, angle_lim
    cz0, cz1 = 0, std_lim
    grid_helper = floating_axes.GridHelperCurveLinear(
        tr, extremes=(ra0, ra1, cz0, cz1),
        grid_locator1=gf,
        tick_formatter1=tf,
        grid_locator2=gfs,
        tick_formatter2=tfs)
    
    ax1 = floating_axes.FloatingSubplot(fig, 111, grid_helper=grid_helper)
    fig.add_subplot(ax1)
    
    ax1.axis["top"].set_axis_direction("bottom")  
    ax1.axis["top"].toggle(ticklabels=True, label=True)
    ax1.axis["top"].major_ticklabels.set_axis_direction("top")
    ax1.axis["top"].label.set_axis_direction("top")
    ax1.axis["top"].label.set_text("Correlation")
    ax1.axis['top'].label.set_size(14)
       
    ax1.axis["left"].set_axis_direction("bottom") 
    ax1.axis["left"].label.set_text("Normalized Standard Deviation")
    ax1.axis['left'].label.set_size(14)

    ax1.axis["right"].set_axis_direction("top")  
    ax1.axis["right"].toggle(ticklabels=True)
    ax1.axis["right"].major_ticklabels.set_axis_direction("left")
    
    ax1.axis["bottom"].set_visible(False) 
    ax1 = ax1.get_aux_axes(tr)
    
    plt.grid(linestyle=':',alpha=0.5)
    
    return fig,ax1

#%% Create a plotting function for Taylor diagrams.

def taylor(scores,colors,units,angle_lim):
    
    import mpl_toolkits.axisartist.floating_axes as floating_axes
    from matplotlib.projections import PolarAxes
    from mpl_toolkits.axisartist.grid_finder import (FixedLocator,
                                                 DictFormatter)

    fig = plt.figure()
    tr = PolarAxes.PolarTransform()
    
    min_corr = np.round(np.cos(angle_lim),1)
    CCgrid= np.concatenate((np.arange(min_corr*10,10,2.0)/10.,[0.9,0.95,0.99]))
    CCpolar=np.arccos(CCgrid)
    gf=FixedLocator(CCpolar)
    tf=DictFormatter(dict(zip(CCpolar, map(str,CCgrid))))
    
    max_std = np.nanmax([scores.OSTD,scores.MSTD])
    
    if np.round(scores.OSTD[0],2)<=0.2:
        #STDgrid=np.linspace(0,np.round(scores.OSTD[0]+0.01,2),3)
        STDgrid=np.linspace(0,np.round(scores.OSTD[0]+max_std+0.02,2),3)
    if np.logical_and(np.round(scores.OSTD[0],2)>0.2,np.round(scores.OSTD[0],2)<=1):
        STDgrid=np.linspace(0,np.round(scores.OSTD[0]+0.1,2),3)
    if np.logical_and(np.round(scores.OSTD[0],2)>1,np.round(scores.OSTD[0],2)<=5):
        STDgrid=np.arange(0,np.round(scores.OSTD[0]+2,2),1)
    if np.round(scores.OSTD[0],2)>5:
        STDgrid=np.arange(0,np.round(scores.OSTD[0]+5,1),2)
    
    gfs=FixedLocator(STDgrid)
    tfs=DictFormatter(dict(zip(STDgrid, map(str,STDgrid))))
    
    ra0, ra1 =0, angle_lim
    if np.round(scores.OSTD[0],2)<=0.2:
        cz0, cz1 = 0, np.round(max_std+0.1,2)        
    else:
        #cz0, cz1 = 0, np.round(scores.OSTD[0]+0.1,2)
        cz0, cz1 = 0, np.round(max_std+0.1,2)
    grid_helper = floating_axes.GridHelperCurveLinear(
        tr, extremes=(ra0, ra1, cz0, cz1),
        grid_locator1=gf,
        tick_formatter1=tf,
        grid_locator2=gfs,
        tick_formatter2=tfs)
    
    ax1 = floating_axes.FloatingSubplot(fig, 111, grid_helper=grid_helper)
    fig.add_subplot(ax1)
    
    ax1.axis["top"].set_axis_direction("bottom")  
    ax1.axis["top"].toggle(ticklabels=True, label=True)
    ax1.axis["top"].major_ticklabels.set_axis_direction("top")
    ax1.axis["top"].label.set_axis_direction("top")
    ax1.axis["top"].label.set_text("Correlation")
    ax1.axis['top'].label.set_size(14)
   
    
    ax1.axis["left"].set_axis_direction("bottom") 
    ax1.axis["left"].label.set_text("Standard Deviation "+ '(' + units +')' )
    ax1.axis['left'].label.set_size(14)

    ax1.axis["right"].set_axis_direction("top")  
    ax1.axis["right"].toggle(ticklabels=True)
    ax1.axis["right"].major_ticklabels.set_axis_direction("left")
    
    ax1.axis["bottom"].set_visible(False) 
    ax1 = ax1.get_aux_axes(tr)
    
    plt.grid(linestyle=':',alpha=0.5)
    
    for i,r in enumerate(scores.iterrows()):
        theta=np.arccos(r[1].CORRELATION)
        rr=r[1].MSTD
        
        ax1.plot(theta,rr,'o',label=r[0],color = colors[i])
    
    ax1.plot(0,scores.OSTD[0],'o',label='Obs')    
    plt.legend(loc='upper right',bbox_to_anchor=[1.3,1.15])    
    plt.show()
    
    rs,ts = np.meshgrid(np.linspace(0,np.round(max_std+0.1,2)),np.linspace(0,angle_lim))
    
    rms = np.sqrt(scores.OSTD[0]**2 + rs**2 - 2*rs*scores.OSTD[0]*np.cos(ts))
    
    ax1.contour(ts, rs, rms,5,colors='0.5')
    #contours = ax1.contour(ts, rs, rms,5,colors='0.5')
    #plt.clabel(contours, inline=1, fontsize=10)
    plt.grid(linestyle=':',alpha=0.5)
    
    for i,r in enumerate(scores.iterrows()):
        #crmse = np.sqrt(r[1].OSTD**2 + r[1].MSTD**2 \
        #           - 2*r[1].OSTD*r[1].MSTD*r[1].CORRELATION) 
        crmse = np.sqrt(scores.OSTD[0]**2 + r[1].MSTD**2 \
                   - 2*scores.OSTD[0]*r[1].MSTD*r[1].CORRELATION) 
        print(crmse)
        c1 = ax1.contour(ts, rs, rms,[crmse],colors=colors[i])
        plt.clabel(c1, inline=1, fontsize=10,fmt='%1.2f')
        
    return fig,ax1
        
#%% Create a plotting function for normalized Taylor diagrams.

def taylor_normalized(scores,colors,angle_lim):
    
    import mpl_toolkits.axisartist.floating_axes as floating_axes
    from matplotlib.projections import PolarAxes
    from mpl_toolkits.axisartist.grid_finder import (FixedLocator,
                                                 DictFormatter)

    fig = plt.figure()
    tr = PolarAxes.PolarTransform()
    
    min_corr = np.round(np.cos(angle_lim),1)
    CCgrid= np.concatenate((np.arange(min_corr*10,10,2.0)/10.,[0.9,0.95,0.99]))
    CCpolar=np.arccos(CCgrid)
    gf=FixedLocator(CCpolar)
    tf=DictFormatter(dict(zip(CCpolar, map(str,CCgrid))))
    
    STDgrid=np.arange(0,2.0,.5)
    gfs=FixedLocator(STDgrid)
    tfs=DictFormatter(dict(zip(STDgrid, map(str,STDgrid))))
    
    ra0, ra1 =0, angle_lim
    cz0, cz1 = 0, 2
    grid_helper = floating_axes.GridHelperCurveLinear(
        tr, extremes=(ra0, ra1, cz0, cz1),
        grid_locator1=gf,
        tick_formatter1=tf,
        grid_locator2=gfs,
        tick_formatter2=tfs)
    
    ax1 = floating_axes.FloatingSubplot(fig, 111, grid_helper=grid_helper)
    fig.add_subplot(ax1)
    
    ax1.axis["top"].set_axis_direction("bottom")  
    ax1.axis["top"].toggle(ticklabels=True, label=True)
    ax1.axis["top"].major_ticklabels.set_axis_direction("top")
    ax1.axis["top"].label.set_axis_direction("top")
    ax1.axis["top"].label.set_text("Correlation")
    ax1.axis['top'].label.set_size(14)
       
    ax1.axis["left"].set_axis_direction("bottom") 
    ax1.axis["left"].label.set_text("Normalized Standard Deviation")
    ax1.axis['left'].label.set_size(14)

    ax1.axis["right"].set_axis_direction("top")  
    ax1.axis["right"].toggle(ticklabels=True)
    ax1.axis["right"].major_ticklabels.set_axis_direction("left")
    
    ax1.axis["bottom"].set_visible(False) 
    ax1 = ax1.get_aux_axes(tr)
    
    plt.grid(linestyle=':',alpha=0.5)
   
    for i,r in enumerate(scores.iterrows()):
        theta=np.arccos(r[1].CORRELATION)            
        rr=r[1].MSTD/r[1].OSTD
        print(rr)
        print(theta)
        
        ax1.plot(theta,rr,'o',label=r[0],color = colors[i])
    
    ax1.plot(0,1,'o',label='Obs')    
    plt.legend(loc='upper right',bbox_to_anchor=[1.3,1.15])    
    plt.show()
   
    rs,ts = np.meshgrid(np.linspace(0,2),np.linspace(0,angle_lim))
    rms = np.sqrt(1 + rs**2 - 2*rs*np.cos(ts))
    
    ax1.contour(ts, rs, rms,3,colors='0.5')
    #contours = ax1.contour(ts, rs, rms,3,colors='0.5')
    #plt.clabel(contours, inline=1, fontsize=10)
    plt.grid(linestyle=':',alpha=0.5)
    
    for i,r in enumerate(scores.iterrows()):
        crmse = np.sqrt(1 + (r[1].MSTD/scores.OSTD[i])**2 \
                   - 2*(r[1].MSTD/scores.OSTD[i])*r[1].CORRELATION) 
        print(crmse)
        c1 = ax1.contour(ts, rs, rms,[crmse],colors=colors[i])
        plt.clabel(c1, inline=1, fontsize=10,fmt='%1.2f')

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
date_ini = cycle[0:4]+'-'+cycle[4:6]+'-'+cycle[6:8]+' '+cycle[8:]+':00:00'
tini = datetime.strptime(date_ini,'%Y-%m-%d %H:%M:%S')
tend = tini + timedelta(hours=126)
date_end = str(tend)

print('Retrieving coordinates from GOFS')
GOFS = xr.open_dataset(url_GOFS,decode_times=False) 

tt_G = GOFS.time
t_G = netCDF4.num2date(tt_G[:],tt_G.units) 

tmin = datetime.strptime(date_ini[0:-6],'%Y-%m-%d %H')
tmax = datetime.strptime(date_end[0:-6],'%Y-%m-%d %H')
oktime_GOFS = np.where(np.logical_and(t_G >= tmin, t_G <= tmax)) 
time_GOFS = np.asarray(t_G[oktime_GOFS])
timestamp_GOFS = mdates.date2num(time_GOFS)

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

#%% Reading glider data

url_glider = gdata[0] #SG665

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

#%% Get glider transect from GOFS 3.1

# Conversion from glider longitude and latitude to GOFS convention
target_lon, target_lat = glider_coor_to_GOFS_coord(long,latg)

# Changing times to timestamp
tstamp_glider = [mdates.date2num(timeg[i]) for i in np.arange(len(timeg))]
tstamp_model = [mdates.date2num(time_GOFS[i]) for i in np.arange(len(time_GOFS))]

# interpolating glider lon and lat to lat and lon on model time
sublon_GOFS = np.interp(tstamp_model,tstamp_glider,target_lon)
sublat_GOFS = np.interp(tstamp_model,tstamp_glider,target_lat)

# Conversion from GOFS convention to glider longitude and latitude
#sublon_GOFSg,sublat_GOFSg = GOFS_coor_to_glider_coord(sublon_GOFS,sublat_GOFS,oklat_GOFS,oklon_GOFS)
sublon_GOFSg,sublat_GOFSg = GOFS_coor_to_glider_coord(sublon_GOFS,sublat_GOFS)

# getting the model grid positions for sublonm and sublatm
oklon_GOFS = np.round(np.interp(sublon_GOFS,lon_G,np.arange(len(lon_G)))).astype(int)
oklat_GOFS = np.round(np.interp(sublat_GOFS,lat_G,np.arange(len(lat_G)))).astype(int)
    
# Getting glider transect from model
target_temp_GOFS, target_salt_GOFS = \
                          get_glider_transect_from_GOFS(GOFS,depth_GOFS,oktime_GOFS,oklat_GOFS,oklon_GOFS)

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

#%% Surface Ocean Heat Content

# glider
OHC_glid = OHC_surface(timeg,tempg_gridded,depthg_gridded,densg_gridded)

# GOFS
OHC_GOFS = OHC_surface(time_GOFS,target_temp_GOFS,depth_GOFS,target_dens_GOFS)

# POM operational    
OHC_POM_oper = OHC_surface(timestamp_POM_oper,target_temp_POM_oper,target_depth_POM_oper,target_dens_POM_oper)

# POM experimental
OHC_POM_exp = OHC_surface(timestamp_POM_exp,target_temp_POM_exp,target_depth_POM_exp,target_dens_POM_exp)

# HYCOM experimental
OHC_HYCOM_exp = OHC_surface(timestamp_HYCOM_exp,target_temp_HYCOM_exp,depth_HYCOM_exp,target_dens_HYCOM_exp)

#%% Calculate T100

# glider
T100_glid = depth_aver_top_100(depthg_gridded,tempg_gridded)

# GOFS
T100_GOFS = depth_aver_top_100(depth_GOFS,target_temp_GOFS)

# POM operational
T100_POM_oper = depth_aver_top_100(target_depth_POM_oper,target_temp_POM_oper) 

# POM experimental
T100_POM_exp = depth_aver_top_100(target_depth_POM_exp,target_temp_POM_exp)

# HYCOM experimental
T100_HYCOM_exp = depth_aver_top_100(depth_HYCOM_exp,target_temp_HYCOM_exp)

#%% Calculate T100_to_ml

dt = 0.2
drho = 0.125

# glider
T100_to_ml_glid = depth_aver_temp_from_100_to_base_mixed_layer\
                  (dt,drho,timeg,depthg_gridded,tempg_gridded,saltg_gridded,densg_gridded)

# GOFS
T100_to_ml_GOFS = depth_aver_temp_from_100_to_base_mixed_layer\
                 (dt,drho,time_GOFS,depth_GOFS,target_temp_GOFS,target_salt_GOFS,target_dens_GOFS)

# POM operational
T100_to_ml_POM_oper = depth_aver_temp_from_100_to_base_mixed_layer\
                      (dt,drho,timestamp_POM_oper,target_depth_POM_oper,target_temp_POM_oper,target_salt_POM_oper,target_dens_POM_oper) 

# POM experimental
T100_to_ml_POM_exp = depth_aver_temp_from_100_to_base_mixed_layer\
                     (dt,drho,timestamp_POM_exp,target_depth_POM_exp,target_temp_POM_exp,target_salt_POM_exp,target_dens_POM_exp) 
                     
# POM experimental
T100_to_ml_HYCOM_exp = depth_aver_temp_from_100_to_base_mixed_layer\
                     (dt,drho,timestamp_HYCOM_exp,depth_HYCOM_exp,target_temp_HYCOM_exp,target_salt_HYCOM_exp,target_dens_HYCOM_exp) 

#%% Vertican profile glider, GOFS 3.1, POM , HYCOM During Dorian
'''
tDorian = datetime(2019,8,29,3)
okg = np.where(timeg < tDorian)[0][-1]
okgofs = np.where(time_GOFS >= tDorian )[0][0] #tdorian
okpom_oper = np.where(timestamp_POM_oper >= mdates.date2num(tDorian))[0][0] #tdorian
okpom_exp = np.where(timestamp_POM_exp >= mdates.date2num(tDorian))[0][0]
okhycom_exp = np.where(timestamp_HYCOM_exp >= mdates.date2num(tDorian))[0][0]

fig,ax = plt.subplots(figsize=(4,7))    
#plt.figure(figsize=(4,7))
plt.plot(tempg[:,okg],-depthg[:,okg],'-',color='royalblue',linewidth=4,label=dataset_id.split('-')[0]) 
plt.plot(target_temp_GOFS[:,okgofs],-depth_GOFS,'s-',color='indianred',linewidth=2,label='GOFS 3.1',markeredgecolor='k',markersize=7)
plt.plot(target_temp_POM_oper[:,okpom_oper],target_depth_POM_oper[:,okpom_oper],'X-',color='mediumorchid',label='HWRF2019-POM Oper',markeredgecolor='k',markersize=7)
plt.plot(target_temp_POM_exp[:,okpom_exp],target_depth_POM_exp[:,okpom_exp],'^-',color='teal',label='HRWF2020-POM Exp',markeredgecolor='k',markersize=7)
plt.plot(target_temp_HYCOM_exp[:,okhycom_exp],-depth_HYCOM_exp,'H-',color='darkorange',label='HWRF2020-HYCOM Exp',markeredgecolor='k',markersize=7)
plt.plot(np.arange(20,31),np.tile(-MLD_dens_crit_glid[okg],len(np.arange(20,31))),'--k')
plt.plot(np.arange(20,31),np.tile(-100,len(np.arange(20,31))),'--k')
ax.fill_between(np.arange(20,31), 0,np.tile(-MLD_dens_crit_glid[okg],len(np.arange(20,31))),\
                color='peachpuff',alpha=0.2)
ax.fill_between(np.arange(20,31),np.tile(-100,len(np.arange(20,31))),\
                                 np.tile(-MLD_dens_crit_glid[okg],len(np.arange(20,31))),\
                color='powderblue',alpha=0.2)  
ax.text(20.5,-MLD_dens_crit_glid[okg]+2,'MLD '+dataset_id.split('-')[0],fontsize=14)
plt.ylim([-200,0])
plt.xlim([20,30])
plt.title('Temperature on '+str(timeg[okg])[0:13],fontsize=16)
plt.ylabel('Depth (m)',fontsize=14)
plt.xlabel('($^oC$)',fontsize=14)
plt.legend(loc='lower right',bbox_to_anchor=(1.6,0),fontsize=14)

file = folder_fig + 'temp_profile_glider_GOFS_POM_' + dataset_id.split('-')[0]
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
'''

#%% Vertican profile glider, GOFS 3.1, POM , HYCOM During Dorian

tDorian = datetime(2019,8,29,3)
okg = np.where(timeg < tDorian)[0][-1]
okgofs = np.where(time_GOFS >= tDorian )[0][0] #tdorian
okpom_oper = np.where(timestamp_POM_oper >= mdates.date2num(tDorian))[0][0] #tdorian
okpom_exp = np.where(timestamp_POM_exp >= mdates.date2num(tDorian))[0][0]
okhycom_exp = np.where(timestamp_HYCOM_exp >= mdates.date2num(tDorian))[0][0]    

fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(8, 7))

ax1.plot(tempg[:,okg],-depthg[:,okg],'-',color='royalblue',linewidth=4,label=dataset_id.split('-')[0]) 
ax1.plot(target_temp_GOFS[:,okgofs],-depth_GOFS,'s-',color='indianred',linewidth=2,label='GOFS 3.1',markeredgecolor='k',markersize=7)
ax1.plot(target_temp_POM_oper[:,okpom_oper],target_depth_POM_oper[:,okpom_oper],'X-',color='mediumorchid',label='HWRF2019-POM (IC clim.)',markeredgecolor='k',markersize=7)
ax1.plot(target_temp_POM_exp[:,okpom_exp],target_depth_POM_exp[:,okpom_exp],'^-',color='teal',linewidth=2,label='HWRF2020-POM (IC RTOFS)',markeredgecolor='k',markersize=7)
ax1.plot(target_temp_HYCOM_exp[:,okhycom_exp],-depth_HYCOM_exp,'H-',color='darkorange',linewidth=2,label='HWRF2020-HYCOM (IC RTOFS)',markeredgecolor='k',markersize=7)
ax1.plot(tempg[:,okg],-depthg[:,okg],'-',color='royalblue',linewidth=4) 
ax1.plot(np.arange(20,31),np.tile(-MLD_dens_crit_glid[okg],len(np.arange(20,31))),'--k')
ax1.plot(np.arange(20,31),np.tile(-100,len(np.arange(20,31))),'--k')
ax1.fill_between(np.arange(20,31), 0,np.tile(-MLD_dens_crit_glid[okg],len(np.arange(20,31))),\
               color='peachpuff',alpha=0.2)
ax1.fill_between(np.arange(20,31),np.tile(-100,len(np.arange(20,31))),\
                                 np.tile(-MLD_dens_crit_glid[okg],len(np.arange(20,31))),\
               color='powderblue',alpha=0.2)  
ax1.text(20.5,-MLD_dens_crit_glid[okg]+2,'MLD '+dataset_id.split('-')[0],fontsize=14)
ax1.set_ylim([-200,0])
ax1.set_xlim([20,30])
ax1.set_title('Temperature on \n'+str(timeg[okg])[0:13],fontsize=16)
ax1.set_ylabel('Depth (m)',fontsize=14)
ax1.set_xlabel('($^oC$)',fontsize=14)

ax2.plot(densg[:,okg],-depthg[:,okg],'-',color='royalblue',linewidth=4,label=dataset_id.split('-')[0]) 
ax2.plot(target_dens_GOFS[:,okgofs],-depth_GOFS,'s-',color='indianred',linewidth=2,label='GOFS 3.1',markeredgecolor='k',markersize=7)
ax2.plot(target_dens_POM_oper[:,okpom_oper],target_depth_POM_oper[:,okpom_oper],'X-',color='mediumorchid',label='HWRF2019-POM (IC clim.)',markeredgecolor='k',markersize=7)
ax2.plot(target_dens_POM_exp[:,okpom_exp],target_depth_POM_exp[:,okpom_exp],'^-',color='teal',label='HRWF2020-POM (IC RTOFS)',markeredgecolor='k',markersize=7)
ax2.plot(target_dens_HYCOM_exp[:,okhycom_exp],-depth_HYCOM_exp,'H-',color='darkorange',linewidth=2,label='HWRF2020-HYCOM (IC RTOFS)',markeredgecolor='k',markersize=7)
ax2.plot(densg[:,okg],-depthg[:,okg],'-',color='royalblue',linewidth=4)
ax2.plot(np.arange(1022,1028),np.tile(-MLD_dens_crit_glid[okg],len(np.arange(1022,1028))),'--k')
ax2.plot(np.arange(1022,1028),np.tile(-100,len(np.arange(1022,1028))),'--k')
ax2.fill_between(np.arange(1022,1028), 0,np.tile(-MLD_dens_crit_glid[okg],len(np.arange(1022,1028))),\
               color='peachpuff',alpha=0.2)
ax2.fill_between(np.arange(1022,1028),np.tile(-100,len(np.arange(1022,1028))),\
                                 np.tile(-MLD_dens_crit_glid[okg],len(np.arange(1022,1028))),\
               color='powderblue',alpha=0.2) 
ax2.set_ylim([-200,0])
ax2.set_xlim([1022,1027])
ax2.set_yticks([])
ax2.set_xlabel('($kg/m^3$)',fontsize=14)
ax2.set_title('Density on \n'+str(timeg[okg])[0:13],fontsize=16)
ax2.legend(loc='lower left',bbox_to_anchor=(-0.75,0),fontsize=14)

file = folder_fig + 'temp_dens_profile_glider_GOFS_POM_' + dataset_id.split('-')[0]
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
