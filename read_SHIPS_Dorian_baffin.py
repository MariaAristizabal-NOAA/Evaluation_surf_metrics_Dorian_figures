#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 11:40:58 2020

@author: aristizabal
"""

#%%

year = '2019'
basin = 'NATL'
name_storm = 'DORIAN05L'
cycle = '2019082800'

folder_SHIPS_Dorian = '/home/aristizabal/SHIPS_shear_Dorian_2019/'

# figures
folder_fig = '/www/web/rucool/aristizabal/Figures/'

#%%

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
import os
import glob

# Increase fontsize of labels globally
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.rc('legend',fontsize=14)

#%%

SHIPS_Dorian1 = sorted(glob.glob(os.path.join(folder_SHIPS_Dorian,'*POM*.txt')))
SHIPS_Dorian2 = sorted(glob.glob(os.path.join(folder_SHIPS_Dorian,'*HYCOM*.txt')))
SHIPS_Dorian = np.hstack([SHIPS_Dorian1,SHIPS_Dorian2])

markers = ['X-','^-','H-']
colors = ['mediumorchid','teal','darkorange']
labels = ['HWRF2019-POM (IC clim.) ','HWRF2020-POM (IC RTOFS) ','HWRF2020-HYCOM (IC RTOFS) ']

#%%
fig,ax1 = plt.subplots(figsize=(10, 5))
plt.title('Shear Magnitude Forecast Dorian Cycle ' + cycle, fontsize=16)
plt.ylabel('Shear Magnitude (KT)',fontsize=14)
plt.xlabel('Forecast Lead Time (Hours)',fontsize=14)
ax1.xaxis.set_major_locator(MultipleLocator(12))
ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax1.xaxis.set_minor_locator(MultipleLocator(3))

for f,file in enumerate(SHIPS_Dorian):
    print(file)
    r = open(file)
    data = r.read()

    for s in data.split('\n'):
        if s[0:4] == 'TIME':
            lead_h = s.split()[2:]
        if s[0:7] == 'SHR_MAG':
            SHR_M = s.split()[2:]

    lead_hours = np.asarray([int(hh) for hh in lead_h])
    SHR_MAG_HWRF = np.asarray([float(ss) for ss in SHR_M])

    plt.plot(lead_hours,SHR_MAG_HWRF,markers[f],color=colors[f],label=labels[f],markeredgecolor='k',markersize=7)
    plt.xlim(lead_hours[0],lead_hours[-1])
    plt.ylim([0,30])
    plt.legend()

file_name = folder_fig + 'Dorian_shear_magnitude_cycle_' + cycle
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1)

#%%
fig,ax1 = plt.subplots(figsize=(10, 5))
plt.title('Translation Speed Forecast Dorian Cycle ' + cycle, fontsize=16)
plt.ylabel('Translation Speed (m/s)',fontsize=14)
plt.xlabel('Forecast Lead Time (Hours)',fontsize=14)
ax1.xaxis.set_major_locator(MultipleLocator(12))
ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax1.xaxis.set_minor_locator(MultipleLocator(3))

for f,file in enumerate(SHIPS_Dorian):
    print(file)
    r = open(file)
    data = r.read()

    for s in data.split('\n'):
        if s[0:4] == 'TIME':
            lead_h = s.split()[2:]
        if s[0:7] == 'STM_SPD':
            STM_S = s.split()[2:]

    lead_hours = np.asarray([int(hh) for hh in lead_h])
    STM_SPD_HWRF = np.asarray([float(ss) for ss in STM_S])* 0.5144

    plt.plot(lead_hours,STM_SPD_HWRF,markers[f],color=colors[f],label=labels[f],markeredgecolor='k',markersize=7)
    plt.xlim(lead_hours[0],lead_hours[-1])
    plt.ylim([0,10])
    plt.legend()

file_name = folder_fig + 'Dorian_trans_speed_cycle_' + cycle
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1)

#%%

fig, ax = plt.subplots(figsize=(8, 8))
grid = plt.GridSpec(2, 1, wspace=0.3, hspace=0.2,left=0.05,right=0.95)

ax1 = plt.subplot(grid[0, 0])
plt.title('Hurricane Dorian Forecast Cycle '+ cycle + '\n Shear Magnitude', fontsize=16)
plt.ylabel('Shear Magnitude (KT)',fontsize=14)
#plt.xlabel('Forecast Lead Time (Hours)',fontsize=14)
ax1.xaxis.set_major_locator(MultipleLocator(12))
ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax1.xaxis.set_minor_locator(MultipleLocator(3))
ax1.xaxis.set_ticklabels([])
ax1.set_xlim([0,84])
plt.plot(np.arange(84),np.tile(10,len(np.arange(84))),'--k')
plt.plot(np.arange(84),np.tile(20,len(np.arange(84))),'--k')
plt.text(6,4,'Low',fontsize=18)
plt.text(6,14,'Moderate',fontsize=18)
plt.text(6,24,'High',fontsize=18)
plt.text(0,31,'(a)',fontsize=16)

for f,file in enumerate(SHIPS_Dorian):
    print(file)
    r = open(file)
    data = r.read()

    for s in data.split('\n'):
        if s[0:4] == 'TIME':
            lead_h = s.split()[2:]
        if s[0:7] == 'SHR_MAG':
            SHR_M = s.split()[2:]

    lead_hours = np.asarray([int(hh) for hh in lead_h])
    SHR_MAG_HWRF = np.asarray([float(ss) for ss in SHR_M])

    plt.plot(lead_hours,SHR_MAG_HWRF,markers[f],color=colors[f],label=labels[f],markeredgecolor='k',markersize=7)
    plt.ylim([0,30])
    plt.legend(loc='upper right',fontsize=14)

ax2 = plt.subplot(grid[1, 0])
plt.title('Translation Speed', fontsize=16)
plt.ylabel('Translation Speed (m/s)',fontsize=14)
plt.xlabel('Forecast Lead Time (Hours)',fontsize=14)
ax2.xaxis.set_major_locator(MultipleLocator(12))
ax2.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax2.xaxis.set_minor_locator(MultipleLocator(3))
ax2.xaxis.set_ticks(np.arange(0,86,12))
ax2.xaxis.set_ticklabels(['28-Aug \n 0','\n 12','29-Aug \n 24','\n 36','30-Aug \n 48',\
                          '\n 60','31-Aug \n 72','\n 84']) #,'01-Sep' \n 96','\n 108','02-Sep \n 120'])
ax1.set_xlim([0,84])
plt.text(0,10.3,'(b)',fontsize=16)

for f,file in enumerate(SHIPS_Dorian):
    print(file)
    r = open(file)
    data = r.read()

    for s in data.split('\n'):
        if s[0:4] == 'TIME':
            lead_h = s.split()[2:]
        if s[0:7] == 'STM_SPD':
            STM_S = s.split()[2:]

    lead_hours = np.asarray([int(hh) for hh in lead_h])
    STM_SPD_HWRF = np.asarray([float(ss) for ss in STM_S])* 0.5144

    plt.plot(lead_hours,STM_SPD_HWRF,markers[f],color=colors[f],label=labels[f],markeredgecolor='k',markersize=7)
    plt.ylim([0,10])
    plt.xlim([0,84])

file_name = folder_fig + 'Shear_mag_trans_speed_cycle_' + cycle
plt.savefig(file_name,bbox_inches = 'tight',pad_inches = 0.1)
