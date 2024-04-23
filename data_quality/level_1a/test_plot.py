#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 14:01:54 2019

@author: idariash
"""

import pyart
import netCDF4
import datetime
from matplotlib import pyplot as plt


#filename_RHI = '/net/denali/storage2/radar2/tmp/RELAMPAGO/RAW/2018/11/29/COL20181129_170604'
filename_RHI = '/net/denali/storage2/radar2/tmp/Ivan/CSU/CHIVO/src/python/data_quality/level_1a/output/chivo.1a.20181129_170604.REL_RHI30.nc'

radar = pyart.io.read(filename_RHI)
angles = radar.fixed_angle['data'];
gatefilter = pyart.filters.GateFilter(radar)
gatefilter.exclude_transition()
gatefilter.exclude_invalid('velocity')
gatefilter.exclude_invalid('reflectivity')
gatefilter.exclude_outside('reflectivity', 0, 100)
gatefilter.exclude_outside('normalized_coherent_power', 0.2, 1)
gatefilter.exclude_outside('cross_correlation_ratio', 0.6, 1)
nyq = radar.instrument_parameters['nyquist_velocity']['data'][0]
time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
time_text = ' ' + time_start.isoformat() + ' UTC '
display = pyart.graph.RadarDisplay(radar)


fig = plt.figure(figsize = [15,10])
x = 1
#file_name = filename[len(filename)-15:len(filename)]
#fig_name_prefix = PNG_Path + 'RHI/Many_Fields/'
    
#  Reflectivity                
ax = fig.add_subplot(321)
display.plot_rhi('reflectivity', sweep = x , ax = ax, #axislabels = ['', 'Distance Above radar (km)'],
                title = '', colorbar_label='Reflectivity (dBZ)',
                vmin = 0, vmax = 60, mask_outside = False,
                cmap = pyart.graph.cm.NWSRef,
                gatefilter = gatefilter)
  
display.plot_grid_lines(ax=None, col='k', ls=':')
display.set_limits(xlim=[0,150], ylim=[0, 21])

#   Differential Reflectivity
ax = fig.add_subplot(322)                
display.plot_rhi('differential_reflectivity', sweep = x , ax = ax, #axislabels = ['', ''],
                vmin = -2, vmax = 6, mask_outside = False,
                cmap = pyart.graph.cm.RefDiff, 
                title='', colorbar_label='Diff. Reflectivity (dB)',
                gatefilter = gatefilter)
  
display.plot_grid_lines(ax=None, col='k', ls=':')
display.set_limits(xlim=[0,150], ylim=[0, 21])

#   RHOHV
ax = fig.add_subplot(323)
display.plot_rhi('cross_correlation_ratio', sweep = x , ax =  ax, #axislabels = ['', 'Distance Above radar (km)'],
                vmin = 0.7, vmax = 1, mask_outside = False,
                title='', colorbar_label='Cross Polar Correlation Ratio', 
                cmap = 'viridis', gatefilter = gatefilter)
  
display.plot_grid_lines(ax=None, col='k', ls=':')
display.set_limits(xlim=[0,150], ylim=[0, 21])

#   Specific Differential Phase
ax = fig.add_subplot(324)
display.plot_rhi('specific_differential_phase', sweep = x , ax =  ax, #axislabels = ['', ''],
                vmin = -1, vmax = 6, mask_outside = False,
                title='', colorbar_label= 'Specific Diff. Phase(deg/km)',
                gatefilter = gatefilter)
  
display.plot_grid_lines(ax=None, col='k', ls=':')
display.set_limits(xlim=[0,150], ylim=[0, 21])


#   Velocity
ax = fig.add_subplot(325)
display.plot_rhi('velocity', sweep = x , ax =  ax,  #axislabels = ['', ''],
                vmin = -nyq, vmax = nyq, mask_outside = False,
                title='', colorbar_label= 'Radial Velocity (m/s)',
                cmap = pyart.graph.cm.NWSVel, gatefilter = gatefilter)
  
display.plot_grid_lines(ax=None, col='k', ls=':')
display.set_limits(xlim=[0,150], ylim=[0, 21])

  
#   Differential Phase
ax = fig.add_subplot(326)
display.plot_rhi('differential_phase', sweep = x , ax =  ax, #axislabels = ['Distance from radar (km)', ''],
                vmin = -180, vmax = 180, #mask_outside = False,
                title='', colorbar_label='Differential Phase (deg.)',
                gatefilter = gatefilter)
  
display.plot_grid_lines(ax=None, col='k', ls=':')
display.set_limits(xlim=[0,150], ylim=[0, 21])


    
plt.suptitle('\n CHIVO| Azimuth: ' + str( round(angles[x],2)) + ' Deg.| ' + time_text, fontsize=16)

#%%

fig = plt.figure(figsize = [7,5])
display.plot_rhi('differential_phase', sweep = x , #ax =  ax, #axislabels = ['Distance from radar (km)', ''],
                vmin = -180, vmax = 180, mask_outside = False)
                #title='', colorbar_label='Differential Phase (deg.)', cmap = 'jet',
                #gatefilter = gatefilter)
  
#display.plot_grid_lines(ax=None, col='k', ls=':')
#display.set_limits(xlim=[0,150], ylim=[0, 21])
