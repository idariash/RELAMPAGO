#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 14:04:57 2019

@author: idariash
"""
import matplotlib.pyplot as plt
import pyart
import os
import numpy as np
import sys
print (sys.getdefaultencoding())

def Zdr_correction(radar, Zdr_bias):
    radar.fields['differential_reflectivity']['data'] = radar.fields['differential_reflectivity']['data'] - Zdr_bias
    
def azimuth_correction(radar, azimuth_error):
    radar.azimuth['data'] = radar.azimuth['data'] - azimuth_error

#%%

filename = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Velocity/Data/20190125/RAW/COL20190125_193048'
#'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Velocity/Data/20190125/RAW/COL20190125_210046'
radar = pyart.io.read(filename)
azimuth_correction(radar, 13.5)
angles = radar.fixed_angle['data'];
nyq = radar.instrument_parameters['nyquist_velocity']['data'][0]
display = pyart.graph.RadarDisplay(radar)

#%% Filter



gatefilter = pyart.filters.GateFilter(radar)
gatefilter.exclude_transition()
gatefilter.exclude_invalid('velocity')
gatefilter.exclude_invalid('reflectivity')
gatefilter.exclude_outside('reflectivity', 0, 100)
#gatefilter.exclude_outside('normalized_coherent_power', 0.2, 1)
#gatefilter.exclude_outside('cross_correlation_ratio', 0.6, 1)

#%%
# (xlim=[-50,-30], ylim=[-30, -10])
x = 7
fig = plt.figure(figsize = [12,10])
display.plot_ppi('reflectivity', sweep = x , #axislabels = ['', 'North South distance from radar (km)'],
                #title = 'Reflectivity', colorbar_label='Z (dBZ)',
                vmin = 0, vmax = 60, mask_outside = False,
                cmap = pyart.graph.cm.NWSRef,
                gatefilter = gatefilter)
display.plot_range_rings([50, 100, 150])
display.plot_grid_lines(ax=None, col='k', ls=':')
display.set_limits (xlim=[-50,-30], ylim=[-30, -10])

#%show velocity, sqi and spectral width

#   Spectum Width
#ax = fig.add_subplot(232)
fig = plt.figure(figsize = [12,10])               
display.plot_ppi('spectrum_width', sweep = x  , #axislabels = ['', ''],
                vmin = 0, vmax = 6, mask_outside = False,
                #cmap = pyart.graph.cm.RefDiff, 
                #title='spectrum_width', colorbar_label='Zdr (dB)',
                gatefilter = gatefilter)
display.plot_range_rings([50, 100, 150])
display.plot_grid_lines(ax=None, col='k', ls=':')
display.set_limits (xlim=[-50,-30], ylim=[-30, -10])

#%

#   Velocity
#ax = fig.add_subplot(233)
fig = plt.figure(figsize = [12,10])
display.plot_ppi('velocity', sweep = x,  #axislabels = ['', ''],
                vmin = -nyq, vmax = nyq, mask_outside = False,
                #title='Radial Velocity', colorbar_label= 'Vel. (m/s)',
                cmap = pyart.graph.cm.NWSVel, gatefilter = gatefilter)
display.plot_range_rings([50, 100, 150])
display.plot_grid_lines(ax=None, col='k', ls=':')
display.set_limits (xlim=[-50,-30], ylim=[-30, -10])
  
#%
#   SQI
#ax = fig.add_subplot(234)
fig = plt.figure(figsize = [12,10])
display.plot_ppi('normalized_coherent_power', sweep = x, #axislabels = ['km', 'km'],
                vmin = 0, vmax = 1, mask_outside = False, gatefilter = gatefilter)
                #title='Cross Polar Correlation Ratio', colorbar_label='RhoHV', 
                #cmap = 'viridis', gatefilter = gatefilter)
display.plot_range_rings([50, 100, 150])
display.plot_grid_lines(ax=None, col='k', ls=':')
display.set_limits (xlim=[-50,-30], ylim=[-30, -10])