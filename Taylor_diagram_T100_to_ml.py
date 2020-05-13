#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 12 17:38:31 2020

@author: aristizabal
"""

#%% User input
folder_fig = '/Users/aristizabal/Desktop/my_papers/Evaluation_surf_metrics_Dorian/figures/'
#folder_fig = '/home/aristizabal/Figures/'

#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#import xarray as xr
#import netCDF4
#from datetime import datetime
#from datetime import timedelta
#import matplotlib.dates as mdates
import sys
#import seawater as sw
#import os
#import os.path
#import glob 

sys.path.append('/home/aristizabal/glider_model_comparisons_Python')
#sys.path.append('/Users/aristizabal/Desktop/MARACOOS_project/Maria_scripts/All_code/Remote_repos/glider_model_comparisons_Python')
#from read_glider_data import read_glider_data_thredds_server
#from process_glider_data import grid_glider_data_thredd

import sys
#sys.path.append('/Users/aristizabal/Desktop/MARACOOS_project/NCEP_scripts')
sys.path.append('/home/aristizabal/NCEP_scripts')
#from utils4HYCOM import readBinz
#from utils4HYCOM import readgrids

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

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
        
#%% Load all data frames 

import feather

DF_GOFS_temp_salt = feather.read_dataframe('DF_GOFS_temp_salt.feather')
DF_GOFS_MLD = feather.read_dataframe('DF_GOFS_MLD.feather')
DF_GOFS_OHC = feather.read_dataframe('DF_GOFS_OHC.feather')
DF_GOFS_T100 = feather.read_dataframe('DF_GOFS_T100.feather')
DF_GOFS_T100_to_ml = feather.read_dataframe('DF_GOFS_T100_to_ml.feather')

DF_POM_temp_salt = feather.read_dataframe('DF_POM_temp_salt.feather')
DF_POM_MLD = feather.read_dataframe('DF_POM_MLD.feather')
DF_POM_OHC = feather.read_dataframe('DF_POM_OHC.feather')
DF_POM_T100 = feather.read_dataframe('DF_POM_T100.feather')
DF_POM_T100_to_ml = feather.read_dataframe('DF_POM_T100_to_ml.feather')

DF_HYCOM_temp_salt = feather.read_dataframe('DF_HYCOM_temp_salt.feather')
DF_HYCOM_MLD = feather.read_dataframe('DF_HYCOM_MLD.feather')
DF_HYCOM_OHC = feather.read_dataframe('DF_HYCOM_OHC.feather')
DF_HYCOM_T100 = feather.read_dataframe('DF_HYCOM_T100.feather')
DF_HYCOM_T100_to_ml = feather.read_dataframe('DF_HYCOM_T100_to_ml.feather')

#%% Temperature statistics.

DF_GOFS = DF_GOFS_temp_salt.dropna()
DF_POM = DF_POM_temp_salt.dropna()
DF_HYCOM = DF_HYCOM_temp_salt.dropna()
       
NGOFS = len(DF_GOFS)-1  #For Unbiased estimmator.
NPOM = len(DF_POM)-1
NHYCOM = len(DF_HYCOM)-1

cols = ['CORRELATION','OSTD','MSTD','CRMSE','BIAS']
    
tskill = np.empty((4,5))
tskill[:] = np.nan

#CORR
tskill[0,0] = DF_GOFS.corr()['temp_obs']['temp_GOFS']
tskill[1,0] = DF_POM.corr()['temp_obs_to_oper']['temp_POM_oper']
tskill[2,0] = DF_POM.corr()['temp_obs_to_exp']['temp_POM_exp']
tskill[3,0] = DF_HYCOM.corr()['temp_obs_to_exp']['temp_HYCOM_exp']

#OSTD
tskill[0,1] = DF_GOFS.std().temp_obs
tskill[1,1] = DF_POM.std().temp_obs_to_oper
tskill[2,1] = DF_POM.std().temp_obs_to_exp
tskill[3,1] = DF_HYCOM.std().temp_obs_to_exp

#MSTD
tskill[0,2] = DF_GOFS.std().temp_GOFS
tskill[1,2] = DF_POM.std().temp_POM_oper
tskill[2,2] = DF_POM.std().temp_POM_exp
tskill[3,2] = DF_HYCOM.std().temp_HYCOM_exp

#CRMSE
tskill[0,3] = np.sqrt(np.nansum(((DF_GOFS.temp_obs-DF_GOFS.mean().temp_obs)-\
                                 (DF_GOFS.temp_GOFS-DF_GOFS.mean().temp_GOFS))**2)/NGOFS)
tskill[1,3] = np.sqrt(np.nansum(((DF_POM.temp_obs_to_exp-DF_POM.mean().temp_obs_to_oper)-\
                                 (DF_POM.temp_POM_oper-DF_POM.mean().temp_POM_oper))**2)/NPOM)
tskill[2,3] = np.sqrt(np.nansum(((DF_POM.temp_obs_to_exp-DF_POM.mean().temp_obs_to_exp)-\
                                 (DF_POM.temp_POM_exp-DF_POM.mean().temp_POM_exp))**2)/NPOM)
tskill[3,3] = np.sqrt(np.nansum(((DF_HYCOM.temp_obs_to_exp-DF_HYCOM.mean().temp_obs_to_exp)-\
                                 (DF_HYCOM.temp_HYCOM_exp-DF_HYCOM.mean().temp_HYCOM_exp))**2)/NHYCOM)

#BIAS
tskill[0,4] = DF_GOFS.mean().temp_obs - DF_GOFS.mean().temp_GOFS
tskill[1,4] = DF_POM.mean().temp_obs_to_oper - DF_POM.mean().temp_POM_oper
tskill[2,4] = DF_POM.mean().temp_obs_to_exp - DF_POM.mean().temp_POM_exp
tskill[3,4] = DF_HYCOM.mean().temp_obs_to_exp - DF_HYCOM.mean().temp_HYCOM_exp

#color
colors = ['indianred','seagreen','darkorchid','darkorange']
    
temp_skillscores = pd.DataFrame(tskill,
                        index=['GOFS','POM_oper','POM_exp','HYCOM exp'],
                        columns=cols)
print(temp_skillscores)

#%% Salinity statistics.

DF_GOFS = DF_GOFS_temp_salt.dropna()
DF_POM = DF_POM_temp_salt.dropna()
DF_HYCOM = DF_HYCOM_temp_salt.dropna()

NGOFS = len(DF_GOFS)-1  #For Unbiased estimmator.
NPOM = len(DF_POM)-1
NHYCOM = len(DF_HYCOM)-1

cols = ['CORRELATION','OSTD','MSTD','CRMSE','BIAS']
    
tskill = np.empty((4,5))
tskill[:] = np.nan

#CORR
tskill[0,0] = DF_GOFS.corr()['salt_obs']['salt_GOFS']
tskill[1,0] = DF_POM.corr()['salt_obs_to_oper']['salt_POM_oper']
tskill[2,0] = DF_POM.corr()['salt_obs_to_exp']['salt_POM_exp']
tskill[3,0] = DF_HYCOM.corr()['salt_obs_to_exp']['salt_HYCOM_exp']

#OSTD
tskill[0,1] = DF_GOFS.std().salt_obs
tskill[1,1] = DF_POM.std().salt_obs_to_oper
tskill[2,1] = DF_POM.std().salt_obs_to_exp
tskill[3,1] = DF_HYCOM.std().salt_obs_to_exp

#MSTD
tskill[0,2] = DF_GOFS.std().salt_GOFS
tskill[1,2] = DF_POM.std().salt_POM_oper
tskill[2,2] = DF_POM.std().salt_POM_exp
tskill[3,2] = DF_HYCOM.std().salt_HYCOM_exp

#CRMSE
tskill[0,3] = np.sqrt(np.nansum(((DF_GOFS.salt_obs-DF_GOFS.mean().salt_obs)-\
                                 (DF_GOFS.salt_GOFS-DF_GOFS.mean().salt_GOFS))**2)/NGOFS)
tskill[1,3] = np.sqrt(np.nansum(((DF_POM.salt_obs_to_exp-DF_POM.mean().salt_obs_to_oper)-\
                                 (DF_POM.salt_POM_oper-DF_POM.mean().salt_POM_oper))**2)/NPOM)
tskill[2,3] = np.sqrt(np.nansum(((DF_POM.salt_obs_to_exp-DF_POM.mean().salt_obs_to_exp)-\
                                 (DF_POM.salt_POM_exp-DF_POM.mean().salt_POM_exp))**2)/NPOM)
tskill[3,3] = np.sqrt(np.nansum(((DF_HYCOM.salt_obs_to_exp-DF_HYCOM.mean().salt_obs_to_exp)-\
                                 (DF_HYCOM.salt_HYCOM_exp-DF_HYCOM.mean().salt_HYCOM_exp))**2)/NHYCOM)    

#BIAS
tskill[0,4] = DF_GOFS.mean().salt_obs - DF_GOFS.mean().salt_GOFS
tskill[1,4] = DF_POM.mean().salt_obs_to_oper - DF_POM.mean().salt_POM_oper
tskill[2,4] = DF_POM.mean().salt_obs_to_exp - DF_POM.mean().salt_POM_exp
tskill[3,4] = DF_HYCOM.mean().salt_obs_to_exp - DF_HYCOM.mean().salt_HYCOM_exp

#color
#color
colors = ['indianred','seagreen','darkorchid','darkorange']
    
salt_skillscores = pd.DataFrame(tskill,
                        index=['GOFS','POM_oper','POM_exp','HYCOM_exp'],
                        columns=cols)
print(salt_skillscores)

#%% Mixed layer statistics Temperature.

DF_GOFS = DF_GOFS_MLD.dropna()
DF_POM = DF_POM_MLD.dropna()
DF_HYCOM = DF_HYCOM_MLD.dropna()

NGOFS = len(DF_GOFS)-1  #For Unbiased estimmator.
NPOM = len(DF_POM)-1
NHYCOM = len(DF_HYCOM)-1

cols = ['CORRELATION','OSTD','MSTD','CRMSE','BIAS']
    
tskill = np.empty((4,5))
tskill[:] = np.nan

#CORR
tskill[0,0] = DF_GOFS.corr()['Tmean_obs']['Tmean_GOFS']
tskill[1,0] = DF_POM.corr()['Tmean_obs_to_oper']['Tmean_POM_oper']
tskill[2,0] = DF_POM.corr()['Tmean_obs_to_exp']['Tmean_POM_exp']
tskill[3,0] = DF_HYCOM.corr()['Tmean_obs_to_exp']['Tmean_HYCOM_exp']

#OSTD
tskill[0,1] = DF_GOFS.std().Tmean_obs
tskill[1,1] = DF_POM.std().Tmean_obs_to_oper
tskill[2,1] = DF_POM.std().Tmean_obs_to_exp
tskill[3,1] = DF_HYCOM.std().Tmean_obs_to_exp

#MSTD
tskill[0,2] = DF_GOFS.std().Tmean_GOFS
tskill[1,2] = DF_POM.std().Tmean_POM_oper
tskill[2,2] = DF_POM.std().Tmean_POM_exp
tskill[3,2] = DF_HYCOM.std().Tmean_HYCOM_exp

#CRMSE
tskill[0,3] = np.sqrt(np.nansum(((DF_GOFS.Tmean_obs-DF_GOFS.mean().Tmean_obs)-\
                                 (DF_GOFS.Tmean_GOFS-DF_GOFS.mean().Tmean_GOFS))**2)/NGOFS)
tskill[1,3] = np.sqrt(np.nansum(((DF_POM.Tmean_obs_to_exp-DF_POM.mean().Tmean_obs_to_oper)-\
                                 (DF_POM.Tmean_POM_oper-DF_POM.mean().Tmean_POM_oper))**2)/NPOM)
tskill[2,3] = np.sqrt(np.nansum(((DF_POM.Tmean_obs_to_exp-DF_POM.mean().Tmean_obs_to_exp)-\
                                 (DF_POM.Tmean_POM_exp-DF_POM.mean().Tmean_POM_exp))**2)/NPOM)
tskill[3,3] = np.sqrt(np.nansum(((DF_HYCOM.Tmean_obs_to_exp-DF_HYCOM.mean().Tmean_obs_to_exp)-\
                                 (DF_HYCOM.Tmean_HYCOM_exp-DF_HYCOM.mean().Tmean_HYCOM_exp))**2)/NHYCOM)

#BIAS
tskill[0,4] = DF_GOFS.mean().Tmean_obs - DF_GOFS.mean().Tmean_GOFS
tskill[1,4] = DF_POM.mean().Tmean_obs_to_oper - DF_POM.mean().Tmean_POM_oper
tskill[2,4] = DF_POM.mean().Tmean_obs_to_exp - DF_POM.mean().Tmean_POM_exp
tskill[3,4] = DF_HYCOM.mean().Tmean_obs_to_exp - DF_HYCOM.mean().Tmean_HYCOM_exp

# colors
colors = ['indianred','seagreen','darkorchid','darkorange']
    
Tmean_mld_skillscores = pd.DataFrame(tskill,
                        index=['GOFS','POM_oper','POM_exp', 'HYCOM_exp'],
                        columns=cols)
print(Tmean_mld_skillscores)

#%% Mixed layer statistics Salinity.

DF_GOFS = DF_GOFS_MLD.dropna()
DF_POM = DF_POM_MLD.dropna()
DF_HYCOM = DF_HYCOM_MLD.dropna()

NGOFS = len(DF_GOFS)-1  #For Unbiased estimmator.
NPOM = len(DF_POM)-1
NHYCOM = len(DF_HYCOM)-1

cols = ['CORRELATION','OSTD','MSTD','CRMSE','BIAS']
    
tskill = np.empty((4,5))
tskill[:] = np.nan

#CORR
tskill[0,0] = DF_GOFS.corr()['Smean_obs']['Smean_GOFS']
tskill[1,0] = DF_POM.corr()['Smean_obs_to_oper']['Smean_POM_oper']
tskill[2,0] = DF_POM.corr()['Smean_obs_to_exp']['Smean_POM_exp']
tskill[3,0] = DF_HYCOM.corr()['Smean_obs_to_exp']['Smean_HYCOM_exp']

#OSTD
tskill[0,1] = DF_GOFS.std().Smean_obs
tskill[1,1] = DF_POM.std().Smean_obs_to_oper
tskill[2,1] = DF_POM.std().Smean_obs_to_exp
tskill[3,1] = DF_HYCOM.std().Smean_obs_to_exp

#MSTD
tskill[0,2] = DF_GOFS.std().Smean_GOFS
tskill[1,2] = DF_POM.std().Smean_POM_oper
tskill[2,2] = DF_POM.std().Smean_POM_exp
tskill[3,2] = DF_HYCOM.std().Smean_HYCOM_exp

#CRMSE
tskill[0,3] = np.sqrt(np.nansum(((DF_GOFS.Smean_obs-DF_GOFS.mean().Smean_obs)-\
                                 (DF_GOFS.Smean_GOFS-DF_GOFS.mean().Smean_GOFS))**2)/NGOFS)
tskill[1,3] = np.sqrt(np.nansum(((DF_POM.Smean_obs_to_exp-DF_POM.mean().Smean_obs_to_oper)-\
                                 (DF_POM.Smean_POM_oper-DF_POM.mean().Smean_POM_oper))**2)/NPOM)
tskill[2,3] = np.sqrt(np.nansum(((DF_POM.Smean_obs_to_exp-DF_POM.mean().Smean_obs_to_exp)-\
                                 (DF_POM.Smean_POM_exp-DF_POM.mean().Smean_POM_exp))**2)/NPOM)
tskill[3,3] = np.sqrt(np.nansum(((DF_HYCOM.Smean_obs_to_exp-DF_HYCOM.mean().Smean_obs_to_exp)-\
                                 (DF_HYCOM.Smean_HYCOM_exp-DF_HYCOM.mean().Smean_HYCOM_exp))**2)/NHYCOM)

#BIAS
tskill[0,4] = DF_GOFS.mean().Smean_obs - DF_GOFS.mean().Smean_GOFS
tskill[1,4] = DF_POM.mean().Smean_obs_to_oper - DF_POM.mean().Smean_POM_oper
tskill[2,4] = DF_POM.mean().Smean_obs_to_exp - DF_POM.mean().Smean_POM_exp
tskill[3,4] = DF_HYCOM.mean().Smean_obs_to_exp - DF_HYCOM.mean().Smean_HYCOM_exp

# colors
colors = ['indianred','seagreen','darkorchid','darkorange']
    
Smean_mld_skillscores = pd.DataFrame(tskill,
                        index=['GOFS','POM_oper','POM_exp', 'HYCOM_exp'],
                        columns=cols)
print(Smean_mld_skillscores)

#%% OHC statistics 

DF_GOFS = DF_GOFS_OHC.dropna()
DF_POM = DF_POM_OHC.dropna()
DF_HYCOM = DF_HYCOM_OHC.dropna()

NGOFS = len(DF_GOFS)-1  #For Unbiased estimmator.
NPOM = len(DF_POM)-1
NHYCOM = len(DF_HYCOM)-1

cols = ['CORRELATION','OSTD','MSTD','CRMSE','BIAS']
    
tskill = np.empty((4,5))
tskill[:] = np.nan

#CORR
tskill[0,0] = DF_GOFS.corr()['OHC_obs']['OHC_GOFS']
tskill[1,0] = DF_POM.corr()['OHC_obs_to_oper']['OHC_POM_oper']
tskill[2,0] = DF_POM.corr()['OHC_obs_to_exp']['OHC_POM_exp']
tskill[3,0] = DF_HYCOM.corr()['OHC_obs_to_exp']['OHC_HYCOM_exp']

#OSTD
tskill[0,1] = DF_GOFS.std().OHC_obs
tskill[1,1] = DF_POM.std().OHC_obs_to_oper
tskill[2,1] = DF_POM.std().OHC_obs_to_exp
tskill[3,1] = DF_HYCOM.std().OHC_obs_to_exp

#MSTD
tskill[0,2] = DF_GOFS.std().OHC_GOFS
tskill[1,2] = DF_POM.std().OHC_POM_oper
tskill[2,2] = DF_POM.std().OHC_POM_exp
tskill[3,2] = DF_HYCOM.std().OHC_HYCOM_exp

#CRMSE
tskill[0,3] = np.sqrt(np.nansum(((DF_GOFS.OHC_obs-DF_GOFS.mean().OHC_obs)-\
                                 (DF_GOFS.OHC_GOFS-DF_GOFS.mean().OHC_GOFS))**2)/NGOFS)
tskill[1,3] = np.sqrt(np.nansum(((DF_POM.OHC_obs_to_exp-DF_POM.mean().OHC_obs_to_oper)-\
                                 (DF_POM.OHC_POM_oper-DF_POM.mean().OHC_POM_oper))**2)/NPOM)
tskill[2,3] = np.sqrt(np.nansum(((DF_POM.OHC_obs_to_exp-DF_POM.mean().OHC_obs_to_exp)-\
                                 (DF_POM.OHC_POM_exp-DF_POM.mean().OHC_POM_exp))**2)/NPOM)
tskill[3,3] = np.sqrt(np.nansum(((DF_HYCOM.OHC_obs_to_exp-DF_HYCOM.mean().OHC_obs_to_exp)-\
                                 (DF_HYCOM.OHC_HYCOM_exp-DF_HYCOM.mean().OHC_HYCOM_exp))**2)/NHYCOM)

#BIAS
tskill[0,4] = DF_GOFS.mean().OHC_obs - DF_GOFS.mean().OHC_GOFS
tskill[1,4] = DF_POM.mean().OHC_obs_to_oper - DF_POM.mean().OHC_POM_oper
tskill[2,4] = DF_POM.mean().OHC_obs_to_exp - DF_POM.mean().OHC_POM_exp
tskill[3,4] = DF_HYCOM.mean().OHC_obs_to_exp - DF_HYCOM.mean().OHC_HYCOM_exp

#color
colors = ['indianred','seagreen','darkorchid','darkorange']
    
OHC_skillscores = pd.DataFrame(tskill,
                        index=['GOFS','POM_oper','POM_exp','HYCOM_exp'],
                        columns=cols)
print(OHC_skillscores)

#%% T100 statistics 

DF_GOFS = DF_GOFS_T100.dropna()
DF_POM = DF_POM_T100.dropna()
DF_HYCOM = DF_HYCOM_T100.dropna()

NGOFS = len(DF_GOFS)-1  #For Unbiased estimmator.
NPOM = len(DF_POM)-1
NHYCOM = len(DF_HYCOM)-1

cols = ['CORRELATION','OSTD','MSTD','CRMSE','BIAS']
    
tskill = np.empty((4,5))
tskill[:] = np.nan

#CORR
tskill[0,0] = DF_GOFS.corr()['T100_obs']['T100_GOFS']
tskill[1,0] = DF_POM.corr()['T100_obs_to_oper']['T100_POM_oper']
tskill[2,0] = DF_POM.corr()['T100_obs_to_exp']['T100_POM_exp']
tskill[3,0] = DF_HYCOM.corr()['T100_obs_to_exp']['T100_HYCOM_exp']

#OSTD
tskill[0,1] = DF_GOFS.std().T100_obs
tskill[1,1] = DF_POM.std().T100_obs_to_oper
tskill[2,1] = DF_POM.std().T100_obs_to_exp
tskill[3,1] = DF_HYCOM.std().T100_obs_to_exp

#MSTD
tskill[0,2] = DF_GOFS.std().T100_GOFS
tskill[1,2] = DF_POM.std().T100_POM_oper
tskill[2,2] = DF_POM.std().T100_POM_exp
tskill[3,2] = DF_HYCOM.std().T100_HYCOM_exp

#CRMSE
tskill[0,3] = np.sqrt(np.nansum(((DF_GOFS.T100_obs-DF_GOFS.mean().T100_obs)-\
                                 (DF_GOFS.T100_GOFS-DF_GOFS.mean().T100_GOFS))**2)/NGOFS)
tskill[1,3] = np.sqrt(np.nansum(((DF_POM.T100_obs_to_exp-DF_POM.mean().T100_obs_to_oper)-\
                                 (DF_POM.T100_POM_oper-DF_POM.mean().T100_POM_oper))**2)/NPOM)
tskill[2,3] = np.sqrt(np.nansum(((DF_POM.T100_obs_to_exp-DF_POM.mean().T100_obs_to_exp)-\
                                 (DF_POM.T100_POM_exp-DF_POM.mean().T100_POM_exp))**2)/NPOM)
tskill[3,3] = np.sqrt(np.nansum(((DF_HYCOM.T100_obs_to_exp-DF_HYCOM.mean().T100_obs_to_exp)-\
                                 (DF_HYCOM.T100_HYCOM_exp-DF_HYCOM.mean().T100_HYCOM_exp))**2)/NHYCOM)

#BIAS
tskill[0,4] = DF_GOFS.mean().T100_obs - DF_GOFS.mean().T100_GOFS
tskill[1,4] = DF_POM.mean().T100_obs_to_oper - DF_POM.mean().T100_POM_oper
tskill[2,4] = DF_POM.mean().T100_obs_to_exp - DF_POM.mean().T100_POM_exp
tskill[3,4] = DF_HYCOM.mean().T100_obs_to_exp - DF_HYCOM.mean().T100_HYCOM_exp

#color
colors = ['indianred','seagreen','darkorchid','darkorange']
    
T100_skillscores = pd.DataFrame(tskill,
                        index=['GOFS','POM_oper','POM_exp','HYCOM_exp'],
                        columns=cols)
print(T100_skillscores)

#%% T100_to_ml statistics 

DF_GOFS = DF_GOFS_T100_to_ml.dropna()
DF_POM = DF_POM_T100_to_ml.dropna()
DF_HYCOM = DF_HYCOM_T100_to_ml.dropna()

NGOFS = len(DF_GOFS)-1  #For Unbiased estimmator.
NPOM = len(DF_POM)-1
NHYCOM = len(DF_HYCOM)-1

cols = ['CORRELATION','OSTD','MSTD','CRMSE','BIAS']
    
tskill = np.empty((4,5))
tskill[:] = np.nan

#CORR
tskill[0,0] = DF_GOFS.corr()['T100_to_ml_obs']['T100_to_ml_GOFS']
tskill[1,0] = DF_POM.corr()['T100_to_ml_obs_to_oper']['T100_to_ml_POM_oper']
tskill[2,0] = DF_POM.corr()['T100_to_ml_obs_to_exp']['T100_to_ml_POM_exp']
tskill[3,0] = DF_HYCOM.corr()['T100_to_ml_obs_to_exp']['T100_to_ml_HYCOM_exp']

#OSTD
tskill[0,1] = DF_GOFS.std().T100_to_ml_obs
tskill[1,1] = DF_POM.std().T100_to_ml_obs_to_oper
tskill[2,1] = DF_POM.std().T100_to_ml_obs_to_exp
tskill[3,1] = DF_HYCOM.std().T100_to_ml_obs_to_exp

#MSTD
tskill[0,2] = DF_GOFS.std().T100_to_ml_GOFS
tskill[1,2] = DF_POM.std().T100_to_ml_POM_oper
tskill[2,2] = DF_POM.std().T100_to_ml_POM_exp
tskill[3,2] = DF_HYCOM.std().T100_to_ml_HYCOM_exp

#CRMSE
tskill[0,3] = np.sqrt(np.nansum(((DF_GOFS.T100_to_ml_obs-DF_GOFS.mean().T100_to_ml_obs)-\
                                 (DF_GOFS.T100_to_ml_GOFS-DF_GOFS.mean().T100_to_ml_GOFS))**2)/NGOFS)
tskill[1,3] = np.sqrt(np.nansum(((DF_POM.T100_to_ml_obs_to_exp-DF_POM.mean().T100_to_ml_obs_to_oper)-\
                                 (DF_POM.T100_to_ml_POM_oper-DF_POM.mean().T100_to_ml_POM_oper))**2)/NPOM)
tskill[2,3] = np.sqrt(np.nansum(((DF_POM.T100_to_ml_obs_to_exp-DF_POM.mean().T100_to_ml_obs_to_exp)-\
                                 (DF_POM.T100_to_ml_POM_exp-DF_POM.mean().T100_to_ml_POM_exp))**2)/NPOM)
tskill[3,3] = np.sqrt(np.nansum(((DF_HYCOM.T100_to_ml_obs_to_exp-DF_HYCOM.mean().T100_to_ml_obs_to_exp)-\
                                 (DF_HYCOM.T100_to_ml_HYCOM_exp-DF_HYCOM.mean().T100_to_ml_HYCOM_exp))**2)/NPOM)    

#BIAS
tskill[0,4] = DF_GOFS.mean().T100_to_ml_obs - DF_GOFS.mean().T100_to_ml_GOFS
tskill[1,4] = DF_POM.mean().T100_to_ml_obs_to_oper - DF_POM.mean().T100_to_ml_POM_oper
tskill[2,4] = DF_POM.mean().T100_to_ml_obs_to_exp - DF_POM.mean().T100_to_ml_POM_exp
tskill[3,4] = DF_HYCOM.mean().T100_to_ml_obs_to_exp - DF_HYCOM.mean().T100_to_ml_HYCOM_exp

#color
colors = ['indianred','seagreen','darkorchid','darkorange']
    
T100_to_ml_skillscores = pd.DataFrame(tskill,
                        index=['GOFS','POM_oper','POM_exp','HYCOM_exp'],
                        columns=cols)
print(T100_to_ml_skillscores)

#%% Combine all metrics into one normalized Taylor diagram 

angle_lim = np.pi/2
std_lim = 1.5
fig,ax1 = taylor_template(angle_lim,std_lim)
markers = ['s','X','^','H']
  
scores = temp_skillscores  
for i,r in enumerate(scores.iterrows()):
    theta=np.arccos(r[1].CORRELATION)
    rr=r[1].MSTD/r[1].OSTD
    ax1.plot(theta,rr,markers[i],color = 'darkorange',markersize=8)
ax1.plot(theta,rr,markers[i],label='Temp',color = 'darkorange',markersize=8)
        
scores = Tmean_mld_skillscores  
for i,r in enumerate(scores.iterrows()):
    theta=np.arccos(r[1].CORRELATION)            
    rr=r[1].MSTD/r[1].OSTD
    ax1.plot(theta,rr,markers[i],color = 'darkorchid',markersize=8)
ax1.plot(theta,rr,markers[i],label='MLT',color = 'darkorchid',markersize=8)
      

scores = T100_skillscores  
for i,r in enumerate(scores.iterrows()):
    theta=np.arccos(r[1].CORRELATION)            
    rr=r[1].MSTD/r[1].OSTD
    ax1.plot(theta,rr,markers[i],color = 'royalblue',markersize=8) 
ax1.plot(theta,rr,markers[i],label='T100',color = 'royalblue',markersize=8) 

scores = T100_to_ml_skillscores  
for i,r in enumerate(scores.iterrows()):
    theta=np.arccos(r[1].CORRELATION)            
    rr=r[1].MSTD/r[1].OSTD
    ax1.plot(theta,rr,markers[i],color = 'darkblue',markersize=8) 
ax1.plot(theta,rr,markers[i],label='T100_ml',color = 'darkblue',markersize=8) 
        
ax1.plot(0,1,'o',label='Obs',markersize=8) 
ax1.plot(0,0,'sk',label='GOFS 3.1',markersize=8)
ax1.plot(0,0,'Xk',label='HWRF2019-POM (IC clim.)',markersize=8)
ax1.plot(0,0,'^k',label='HWRF2020-POM (IC RTOFS)',markersize=8)
ax1.plot(0,0,'Hk',label='HWRF2020-HYCOM (IC RTOFS)',markersize=8)
   
plt.legend(loc='upper left',bbox_to_anchor=[0,1.85])    

rs,ts = np.meshgrid(np.linspace(0,std_lim),np.linspace(0,angle_lim))
rms = np.sqrt(1 + rs**2 - 2*rs*np.cos(ts))
    
contours = ax1.contour(ts, rs, rms,3,colors='0.5')
plt.clabel(contours, inline=1, fontsize=10)
plt.grid(linestyle=':',alpha=0.5)

file = folder_fig + 'Taylor_norm_2019082800_v3'
plt.savefig(file,bbox_inches = 'tight',pad_inches = 0.1) 
