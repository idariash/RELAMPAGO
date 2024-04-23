#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 15:46:17 2019

@author: idariash

This file contains fuctions for ploting CHIVO data during RELAMPAGO
"""
import matplotlib.pyplot as plt
import pyart
import matplotlib.colors as colors
import numpy as np
import netCDF4
#import importlib
#importlib.reload(pyart)

def adjust_fhc_colorbar_for_pyart(cb):
    cb.set_ticks(np.arange(5.5, 16.5, 0.92))
    cb.ax.set_yticklabels(['Not Classified', 'Large Drops', 'Drizzle', 'Rain',
                           'Heavy Rain', 'Rain and Hail' , 'Hail',
                           'Graupel', 'Wet Snow', 'Dry Snow', 'Crystals', 'Dendrites'])
    cb.ax.set_ylabel('')
    cb.ax.tick_params(length=0)
    return cb

def chivo_rhi_drops(filename, sweep):
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarMapDisplay(radar)
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('corrected_reflectivity')
    gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 0, 100)
    #gatefilter.exclude_outside('corrected_cross_correlation_ratio', 0.7, 1)
    
    angle = str( round(angles[x] - 7,5)) # azimuth correction of 13 deg
    fig = plt.figure(figsize = [15,10])
    
    #  Reflectivity                
    ax = fig.add_subplot(321)
    display.plot_rhi('corrected_reflectivity', sweep = x , ax = ax, title = '',
                     axislabels = ['', 'Distance Above radar (km)'],
                    vmin = -15, vmax = 80, mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef, colorbar_label =  'Corrected Reflectivity (dBZ)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, 80])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   Differential Reflectivity
    ax = fig.add_subplot(322)                
    display.plot_rhi('corrected_differential_reflectivity', sweep = x , ax = ax,
                     axislabels = ['', ''],
                    vmin = -2, vmax = 8, mask_outside = False,
                    cmap = pyart.graph.cm.RefDiff, title='', colorbar_label='Diff. Reflectivity (dB)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, 80])
    display.plot_grid_lines(ax=None, col='k', ls=':')

    #   KDP
    ax = fig.add_subplot(323)
    display.plot_rhi('corrected_specific_differential_phase', sweep = x , ax =  ax, 
                    axislabels = ['', 'Distance Above radar (km)'],
                    vmin = 0, vmax = 6, mask_outside = False,
                    title='', colorbar_label= 'Specific Diff. Phase (deg/km)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, 80])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   PHIDP
    ax = fig.add_subplot(324)
    display.plot_rhi('corrected_differential_phase', sweep = x , ax =  ax,
                    axislabels = ['', ' '],
                    vmin = -180, vmax = 180, mask_outside = False,
                     title='', colorbar_label='Corrected Diff. Phase (°)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, 80])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   RhoHV
    ax = fig.add_subplot(325)
    display.plot_rhi('corrected_copol_correlation_ratio', sweep = x , ax =  ax,
                     vmin = 0.7, vmax = 1, mask_outside = False, cmap = 'jet',
                     title='', colorbar_label='Copolar Correlcation',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, 80])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   HyC
    
    hid_colors = ['White', 'Black','LightBlue', 'MediumBlue', 'DarkOrange', 'LightPink',
                  'Cyan', 'DarkGray', 'Lime', 'Yellow', 'Red', 'Fuchsia']
    cmaphid = colors.ListedColormap(hid_colors)
    ax = fig.add_subplot(326)
    display.plot_rhi('HydroClass', sweep = x ,ax = ax, 
                     axislabels = ['Distance from radar (km)', ' '],
                     cmap = cmaphid, colorbar_label='',
                     vmin = 5, vmax = 16, mask_outside = False, gatefilter = gatefilter,
                     title = ' ')
    display.set_limits(ylim=[0, 21])
    display.cbs[5] = adjust_fhc_colorbar_for_pyart(display.cbs[5])
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, 80])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    angle = str( round(angles[x],2))    
    plt.suptitle('\n CSU-CHIVO| Azimuth: ' + angle + ' Deg.| ' + time_text, fontsize=16)
    print('check line 385 at ' + '/top/students/GRAD/ECE/idariash/home/anaconda3/lib/python3.6/site-packages/pyart/core/radar.py')
    return fig

def chivo_rhi_drops_ref(filename, sweep, azimuth_offset):
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarMapDisplay(radar)
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('corrected_reflectivity')
    gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 0, 100)
    #gatefilter.exclude_outside('corrected_cross_correlation_ratio', 0.7, 1)
    
    angle = str( round(angles[x] - azimuth_offset)) # azimuth correction of 13 deg
    fig = plt.figure(figsize = [10,4])
    
    #  Reflectivity                
    display.plot_rhi('corrected_reflectivity', sweep = x, title = '',
                    # axislabels = ['', 'Distance Above radar (km)'],
                    vmin = -15, vmax = 80, mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef, colorbar_label =  'Corrected Reflectivity (dBZ)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, 80])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    angle = str( round(angles[x],2))    
    plt.suptitle('CSU-CHIVO| Azimuth: ' + angle + ' Deg.| ' + time_text, fontsize=14)
    print('check line 385 at ' + '/top/students/GRAD/ECE/idariash/home/anaconda3/lib/python3.6/site-packages/pyart/core/radar.py')
    return fig

def chivo_rhi_1a_ref(filename, sweep, azimuth_offset):
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarMapDisplay(radar)
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('reflectivity')
    gatefilter.exclude_invalid('differential_phase')
    gatefilter.exclude_outside('reflectivity', 0, 100)
    #gatefilter.exclude_outside('cross_correlation_ratio', 0.7, 1)
    gatefilter.exclude_outside('normalized_coherent_power', 0.1, 1)
    
    angle = str( round(angles[x] - azimuth_offset)) # azimuth correction of 13 deg
    fig = plt.figure(figsize = [10,4])
    
    #  Reflectivity                
    display.plot_rhi('reflectivity', sweep = x, title = '',
                    # axislabels = ['', 'Distance Above radar (km)'],
                    vmin = 0, vmax = 70, mask_outside = False,
                    colorbar_label =  'Reflectivity (dBZ)', #cmap = 'jet', 
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, 80])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    angle = str( round(angles[x],2))    
    plt.suptitle('CSU-CHIVO| Azimuth: ' + angle + ' Deg.| ' + time_text, fontsize=14)
    print('check line 385 at ' + '/top/students/GRAD/ECE/idariash/home/anaconda3/lib/python3.6/site-packages/pyart/core/radar.py')
    return fig


def chivo_rhi_1b_HyC(filename, sweep, azimuth_offset):
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarMapDisplay(radar)
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('corrected_reflectivity')
    gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 10, 100)
    #gatefilter.exclude_outside('cross_correlation_ratio', 0.7, 1)
    #gatefilter.exclude_outside('normalized_coherent_power', 0.2, 1)
    
    angle = str( round(angles[x] - azimuth_offset)) # azimuth correction of 13 deg
    fig = plt.figure(figsize = [8, 4])
    
    #  HyC
    hid_colors = ['White', 'Black','LightBlue', 'MediumBlue', 'DarkOrange', 'LightPink',
                  'Cyan', 'DarkGray', 'Lime', 'Yellow', 'Red', 'Fuchsia']
    cmaphid = colors.ListedColormap(hid_colors)
    display.plot_rhi('HydroClass', sweep = x , cmap = cmaphid, colorbar_label='',
                     vmin = 5, vmax = 16, mask_outside = False, gatefilter = gatefilter,
                     title = ' ', axislabels = ['', ''],)
    display.set_limits(ylim=[0, 14])
    display.set_limits(xlim=[20, 50])
    display.cbs[0] = adjust_fhc_colorbar_for_pyart(display.cbs[0])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    angle = str( round(angles[x],2))    
    #plt.suptitle('CSU-CHIVO| Azimuth: ' + angle + ' Deg.| ' + time_text, fontsize=14)
    print('check line 385 at ' + '/top/students/GRAD/ECE/idariash/home/anaconda3/lib/python3.6/site-packages/pyart/core/radar.py')
    return fig


def chivo_rhi_1b_ref(filename, sweep, azimuth_offset):
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarMapDisplay(radar)
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('corrected_reflectivity')
    gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 0, 100)
    #gatefilter.exclude_outside('cross_correlation_ratio', 0.7, 1)
    gatefilter.exclude_outside('normalized_coherent_power', 0.1, 1)
    
    angle = str( round(angles[x] - azimuth_offset)) # azimuth correction of 13 deg
    fig = plt.figure(figsize = [10,4])
    
    #  Reflectivity                
    display.plot_rhi('corrected_reflectivity', sweep = x, title = '',
                    # axislabels = ['', 'Distance Above radar (km)'],
                    vmin = 0, vmax = 70, mask_outside = False,
                    colorbar_label =  'Reflectivity (dBZ)', #cmap = 'jet', 
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, 80])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    angle = str( round(angles[x],2))    
    plt.suptitle('CSU-CHIVO| Azimuth: ' + angle + ' Deg.| ' + time_text, fontsize=14)
    print('check line 385 at ' + '/top/students/GRAD/ECE/idariash/home/anaconda3/lib/python3.6/site-packages/pyart/core/radar.py')
    return fig

def chivo_rhi(filename, sweep, xlim):
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarMapDisplay(radar)
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('reflectivity')
    gatefilter.exclude_invalid('differential_phase')
    gatefilter.exclude_outside('reflectivity', 0, 100)
    #gatefilter.exclude_outside('cross_correlation_ratio', 0.3, 1)
    
    angle = str( round(angles[x],2))
    fig = plt.figure(figsize = [15,10])
    
    #  Reflectivity                
    ax = fig.add_subplot(321)
    display.plot_rhi('reflectivity', sweep = x , ax = ax, title = '',
                     axislabels = ['', 'Distance Above radar (km)'],
                    vmin = -15, vmax = 80, mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef, colorbar_label =  'Reflectivity (dBZ)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, xlim])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   Differential Reflectivity
    ax = fig.add_subplot(322)                
    display.plot_rhi('differential_reflectivity', sweep = x , ax = ax,
                     axislabels = ['', ''],
                    vmin = -2, vmax = 8, mask_outside = False,
                    cmap = pyart.graph.cm.RefDiff, title='', colorbar_label='Diff. Reflectivity (dB)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, xlim])
    display.plot_grid_lines(ax=None, col='k', ls=':')

    #   KDP
    ax = fig.add_subplot(323)
    display.plot_rhi('specific_differential_phase', sweep = x , ax =  ax, 
                    axislabels = ['', 'Distance Above radar (km)'],
                    vmin = 0, vmax = 6, mask_outside = False,
                    title='', colorbar_label= 'Specific Diff. Phase (deg/km)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, xlim])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   PHIDP
    ax = fig.add_subplot(324)
    display.plot_rhi('differential_phase', sweep = x , ax =  ax,
                    axislabels = ['', ' '],
                    vmin = 0, vmax = 180, mask_outside = False,
                     title='', colorbar_label='Diff. Phase (°)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, xlim])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   RhoHV
    ax = fig.add_subplot(325)
    display.plot_rhi('cross_correlation_ratio', sweep = x , ax =  ax,
                     vmin = 0.0, vmax = 1, mask_outside = False, cmap = 'jet',
                     title='', colorbar_label='Copolar Correlcation',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, xlim])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   RhoHV
    ax = fig.add_subplot(326)
    display.plot_rhi('spectrum_width', sweep = x , ax =  ax,
                     vmin = 0, vmax = 6, mask_outside = False, 
                     title='', colorbar_label='Spectrum Width (m/s)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, xlim])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
  
        
    plt.suptitle('\n CSU-CHIVO| Azimuth: ' + angle + ' Deg.| ' + time_text, fontsize=16)
    return fig

def chivo_rhi_drops_4dd(filename_drops, filename_sigmet ,sweep, xlim, azimuth_offset):
    radar = pyart.io.read(filename_drops)
    x = sweep
    angles = radar.fixed_angle['data'] # - 13
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
    display = pyart.graph.RadarMapDisplay(radar)
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('corrected_reflectivity')
    gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 0, 100)
    #gatefilter.exclude_outside('corrected_cross_correlation_ratio', 0.6, 1)
    
    #angle = str( round(angles[x],2))
    angle = str( round(angles[x] - azimuth_offset))
    fig = plt.figure(figsize = [15,10])

    #  Reflectivity                
    ax = fig.add_subplot(321)
    display.plot_rhi('corrected_reflectivity', sweep = x , ax = ax, title = '',
                     axislabels = ['', 'Distance Above radar (km)'],
                    vmin = 0, vmax = 80, mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef, colorbar_label =  'Corrected Reflectivity (dbZ)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, xlim])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   Differential Reflectivity
    ax = fig.add_subplot(322)                
    display.plot_rhi('corrected_differential_reflectivity', sweep = x , ax = ax,
                     axislabels = ['', ''],
                    vmin = -2, vmax = 8, mask_outside = False,
                    cmap = pyart.graph.cm.RefDiff, title='', colorbar_label='Diff. Reflectivity (dB)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, xlim])
    display.plot_grid_lines(ax=None, col='k', ls=':')

    #   KDP
    ax = fig.add_subplot(323)
    display.plot_rhi('corrected_specific_differential_phase', sweep = x , ax =  ax, 
                    axislabels = ['', 'Distance Above radar (km)'],
                    vmin = -1, vmax = 6, mask_outside = False,
#                    vmin = 0, vmax = 6, mask_outside = False,
                    title='', colorbar_label= 'Specific Diff. Phase (deg/km)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, xlim])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   Velocity
    last_radar = pyart.io.read(filename_sigmet)
    radar_sigmet = pyart.io.read(filename_sigmet)
    display_sigmet = pyart.graph.RadarMapDisplay(last_radar) # just region base
    #display_sigmet = pyart.graph.RadarMapDisplay(radar_sigmet) # for 4dd
    
    gatefilter_sigmet = pyart.filters.GateFilter(last_radar)
    gatefilter_sigmet.exclude_transition()
    gatefilter_sigmet.exclude_invalid('velocity')
    gatefilter_sigmet.exclude_invalid('reflectivity')
    gatefilter_sigmet.exclude_outside('reflectivity', -20, 80)
    gatefilter_sigmet.exclude_outside('normalized_coherent_power', 0.2, 1)
    gatefilter_sigmet.exclude_outside('cross_correlation_ratio', 0.7, 1)
    
    nyq = radar_sigmet.instrument_parameters['nyquist_velocity']['data'][0]
    corr_vel = pyart.correct.dealias_region_based(
            last_radar, vel_field='velocity', keep_original=False, 
            gatefilter = gatefilter_sigmet, nyquist_vel=nyq, centered = True)
    last_radar.add_field('corrected_velocity', corr_vel, replace_existing = True)
    
    dealias_data = pyart.correct.dealias_fourdd(
            radar_sigmet, last_radar = last_radar, gatefilter = gatefilter_sigmet)
    radar_sigmet.add_field('corrected_velocity', dealias_data, replace_existing = True)
    
    ax = fig.add_subplot(324)
    display_sigmet.plot_rhi('corrected_velocity', sweep = x , ax =  ax,
                    axislabels = ['', ' '],
                    vmin = -25, vmax = 25, mask_outside = False,
                    title='', colorbar_label= 'Unlfolded Velocity (m/s)',
                    gatefilter = gatefilter)
    display_sigmet.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, xlim])
    display_sigmet.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   RhoHV
    ax = fig.add_subplot(325)
    display.plot_rhi('corrected_cross_correlation_ratio', sweep = x , ax =  ax,
                     vmin = 0.7, vmax = 1, mask_outside = False, cmap = 'jet',
                     title='', colorbar_label='Copolar Correlcation',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, xlim])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   HyC
    
    hid_colors = ['White', 'Black','LightBlue', 'MediumBlue', 'DarkOrange', 'LightPink',
                  'Cyan', 'DarkGray', 'Lime', 'Yellow', 'Red', 'Fuchsia']
    cmaphid = colors.ListedColormap(hid_colors)
    ax = fig.add_subplot(326)
    display.plot_rhi('HydroClass', sweep = x ,ax = ax, 
                     axislabels = ['Distance from radar (km)', ' '],
                     cmap = cmaphid, colorbar_label='',
                     vmin = 5, vmax = 16, mask_outside = False, gatefilter = gatefilter,
                     title = ' ')
    display.set_limits(ylim=[0, 21])
    display.cbs[4] = adjust_fhc_colorbar_for_pyart(display.cbs[4])
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, xlim])
    display.plot_grid_lines(ax=None, col='k', ls=':')
        
    plt.suptitle('\n CSU-CHIVO| Azimuth: ' + angle + ' Deg.| ' + time_text, fontsize=16)
    return fig
    
def ppi_drops(filename, sweep):
    
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarDisplay(radar)
    radar.azimuth['data'] = radar.azimuth['data'] - 13
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('corrected_reflectivity')
    gatefilter.exclude_invalid('corrected_differential_phase')
    #gatefilter.exclude_outside('corrected_reflectivity', 0, 100)
    #gatefilter.exclude_outside('corrected_cross_correlation_ratio', 0.7, 1)
    
    
    angle = str( round(angles[x],2))
    fig = plt.figure(figsize = [19,10])
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
           
    #  Reflectivity                
    ax = fig.add_subplot(231)
    display.plot_ppi('corrected_reflectivity', sweep = x , ax = ax, axislabels = ['', 'North South distance from radar (km)'],
                    title = 'Corrected Reflectivity', colorbar_label='Z (dBZ)',
                    vmin = 0, vmax = 80, mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef,
                    gatefilter = gatefilter)
    display.plot_range_rings([50, 100, 150])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-100,100], ylim=[-100, 100])
    
    #   Differential Reflectivity
    ax = fig.add_subplot(232)                
    display.plot_ppi('corrected_differential_reflectivity', sweep = x , ax = ax, axislabels = ['', ''],
                    vmin = -2, vmax = 8, mask_outside = False,
                    cmap = pyart.graph.cm.RefDiff, 
                    title='Diff. Reflectivity', colorbar_label='Zdr (dB)',
                    gatefilter = gatefilter)
    display.plot_range_rings([50, 100, 150])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-100,100], ylim=[-100, 100])

    
    #   HyC
    ax = fig.add_subplot(233)
    hid_colors = ['White', 'Black','LightBlue', 'MediumBlue', 'DarkOrange', 'LightPink',
                  'Cyan', 'DarkGray', 'Lime', 'Yellow', 'Red', 'Fuchsia']
    cmaphid = colors.ListedColormap(hid_colors)
    display.plot_ppi('HydroClass', sweep = x ,ax = ax, 
                     axislabels = ['', ''],
                     cmap = cmaphid, colorbar_label='',
                     vmin = 5, vmax = 16, mask_outside = False, gatefilter = gatefilter,
                     title = 'Hydrometeor Classification')
    display.cbs[2] = adjust_fhc_colorbar_for_pyart(display.cbs[2])
    display.plot_range_rings([50, 100, 150])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-100,100], ylim=[-100, 100])
      
    #   RHOHV
    ax = fig.add_subplot(234)
    display.plot_ppi('corrected_copol_correlation_ratio', sweep = x , ax =  ax, #axislabels = ['km', 'km'],
                    vmin = 0.7, vmax = 1, mask_outside = False,
                    title='Copolar Correlation Ratio', colorbar_label='RhoHV', 
                    cmap = 'jet', gatefilter = gatefilter)
    display.plot_range_rings([50, 100, 150])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-100,100], ylim=[-100, 100])
    
    #   Specific Differential Phase
    ax = fig.add_subplot(235)
    display.plot_ppi('corrected_specific_differential_phase', sweep = x , ax =  ax, axislabels = ['East West distance from radar (km)', ''],
                    vmin = 0, vmax = 6, mask_outside = False,
                    title='Specific Diff. Phase', colorbar_label= 'Kdp (deg/km)',
                    gatefilter = gatefilter)
    display.plot_range_rings([50, 100, 150])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-100,100], ylim=[-100, 100])
      
    #   Differential Phase
    ax = fig.add_subplot(236)
    display.plot_ppi('corrected_differential_phase', sweep = x , ax =  ax, axislabels = ['East West distance from radar (km)', ''],
                    vmin = -180, vmax = 180, mask_outside = False,
                    title='Differential Phase', colorbar_label='Phidp (deg.)', 
                    gatefilter = gatefilter)
    display.plot_range_rings([50, 100, 150])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-100,100], ylim=[-100, 100])
    plt.suptitle('\n CSU-CHIVO| Elevation: ' + angle + ' Deg.| ' + time_text, fontsize=16)
    
    return fig


def chivo_ppi_drops_region(filename, sweep, azimuth_offset):
    
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarDisplay(radar)
    radar.azimuth['data'] = radar.azimuth['data'] - azimuth_offset#13
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('corrected_reflectivity')
    gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 0, 100)
    gatefilter.exclude_outside('corrected_cross_correlation_ratio', 0.7, 1)
    
    
    angle = str( round(angles[x],2))
    fig = plt.figure(figsize = [19,10])
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
           
    #  Reflectivity                
    ax = fig.add_subplot(231)
    display.plot_ppi('corrected_reflectivity', sweep = x , ax = ax, axislabels = ['', 'North South distance from radar (km)'],
                    title = 'Corrected Reflectivity', colorbar_label='Z (dBZ)',
                    vmin = 0, vmax = 80, mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef,
                    gatefilter = gatefilter)
    display.plot_range_rings([25, 50, 75])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-80,80], ylim=[-80,80])
    
    #   Differential Reflectivity
    ax = fig.add_subplot(232)                
    display.plot_ppi('corrected_differential_reflectivity', sweep = x , ax = ax, axislabels = ['', ''],
                    vmin = -2, vmax = 8, mask_outside = False,
                    cmap = pyart.graph.cm.RefDiff, 
                    title='Diff. Reflectivity', colorbar_label='Zdr (dB)',
                    gatefilter = gatefilter)
    display.plot_range_rings([25, 50, 75])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-80,80], ylim=[-80,80])

    
    #   HyC
    ax = fig.add_subplot(233)
    hid_colors = ['White', 'Black','LightBlue', 'MediumBlue', 'DarkOrange', 'LightPink',
                  'Cyan', 'DarkGray', 'Lime', 'Yellow', 'Red', 'Fuchsia']
    cmaphid = colors.ListedColormap(hid_colors)
    display.plot_ppi('HydroClass', sweep = x ,ax = ax, 
                     axislabels = ['', ''],
                     cmap = cmaphid, colorbar_label='',
                     vmin = 5, vmax = 16, mask_outside = False, gatefilter = gatefilter,
                     title = 'Hydrometeor Classification')
    display.cbs[2] = adjust_fhc_colorbar_for_pyart(display.cbs[2])
    display.plot_range_rings([25, 50, 75])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-80,80], ylim=[-80,80])
      
    #   RHOHV
    ax = fig.add_subplot(234)
    display.plot_ppi('corrected_cross_correlation_ratio', sweep = x , ax =  ax, #axislabels = ['km', 'km'],
                    vmin = 0.7, vmax = 1, mask_outside = False,
                    title='Copolar Correlation Ratio', colorbar_label='RhoHV', 
                    cmap = 'jet', gatefilter = gatefilter)
    display.plot_range_rings([25, 50, 75])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-80,80], ylim=[-80,80])
    
    #   Specific Differential Phase
    ax = fig.add_subplot(235)
    display.plot_ppi('corrected_specific_differential_phase', sweep = x , ax =  ax, axislabels = ['East West distance from radar (km)', ''],
                    vmin = 0, vmax = 6, mask_outside = False,
                    title='Specific Diff. Phase', colorbar_label= 'Kdp (deg/km)',
                    gatefilter = gatefilter)
    display.plot_range_rings([25, 50, 75])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-80,80], ylim=[-80,80])
      
    #   Differential Phase
    ax = fig.add_subplot(236)
    display.plot_ppi('corrected_differential_phase', sweep = x , ax =  ax, axislabels = ['East West distance from radar (km)', ''],
                    vmin = -180, vmax = 180, mask_outside = False,
                    title='Differential Phase', colorbar_label='Phidp (deg.)', 
                    gatefilter = gatefilter)
    display.plot_range_rings([25, 50, 75])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-80,80], ylim=[-80,80])
    plt.suptitle('\n CSU-CHIVO| Elevation: ' + angle + ' Deg.| ' + time_text, fontsize=16)
    
    return fig

def chivo_ppi_drops_region_field(filename, sweep, azimuth_offset, field_name):
    
    corrected_field_name = {
            'reflectivity': 'corrected_reflectivity',
            'differential_reflectivity': 'corrected_differential_reflectivity',
            'specific_differential_phase': 'corrected_specific_differential_phase'
            }
    
    colorMap = {
            'reflectivity': pyart.graph.cm.NWSRef,
            'differential_reflectivity': pyart.graph.cm.RefDiff,
            'rhohv': 'jet',
            'specific_differential_phase': pyart.graph.cm.Theodore16            
            
            }
    minVal = {
            'reflectivity': 0,
            'differential_reflectivity': -2,
            'specific_differential_phase': -1
            }
    maxVal = {
            'reflectivity': 80,
            'differential_reflectivity': 8,
            'specific_differential_phase': 6
            }
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarDisplay(radar)
    radar.azimuth['data'] = radar.azimuth['data'] - azimuth_offset
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('corrected_reflectivity')
    gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 18, 100)
    #gatefilter.exclude_outside('corrected_copol_correlation_ratio', 0.7, 1)
    gatefilter.exclude_inside('normalized_coherent_power', 0, 0.2)
    
    
    angle = str( round(angles[x],2))
    fig = plt.figure(figsize = [5,4])
    
    # time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    # time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
    
    # time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    # time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S') #time_start.strftime('%Y-%m-%dT%H:%M:%S')
           

    display.plot_ppi(corrected_field_name[field_name], sweep = x , title = '',
                    vmin = minVal[field_name], vmax = maxVal[field_name], mask_outside = False,
                    cmap = colorMap[field_name],  gatefilter = gatefilter)
   
    display.plot_range_rings([10, 20, 30, 40, 50])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-55, 20], ylim=[-55, 20])
    angle = str( round(angles[x],2))    
    plt.suptitle('CSU-CHIVO| Azimuth: ' + angle + ' Deg.| ' + time_text, fontsize=14)
    
    return fig


 
def chivo_ppi_drops_region_ref(filename, sweep, azimuth_offset):
    
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarDisplay(radar)
    radar.azimuth['data'] = radar.azimuth['data'] - azimuth_offset#13
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('corrected_reflectivity')
    gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 13, 100)
    #gatefilter.exclude_outside('corrected_cross_correlation_ratio', 0.7, 1)
    
    
    angle = str( round(angles[x],2))
    fig = plt.figure(figsize = [10,8])
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
           
    #  Reflectivity                
    display.plot_ppi('corrected_reflectivity', sweep = x , # axislabels = ['', 'North South distance from radar (km)'],
                    title = 'Reflectivity', colorbar_label='Z (dBZ)',
                    vmin = 0, vmax = 80, mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef,
                    gatefilter = gatefilter)
    display.plot_range_rings([25, 50, 75])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-80,80], ylim=[-80,80])
    
    
    plt.suptitle('\n CSU-CHIVO| Elevation: ' + angle + ' Deg.| ' + time_text, fontsize=16)
    
    return fig

def chivo_ppi_1a_vel(filename, sweep):
    
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarDisplay(radar)
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('reflectivity')
    #gatefilter.exclude_invalid('differential_phase')
    gatefilter.exclude_outside('reflectivity', 10, 100)
    #gatefilter.exclude_outside('normalized_coherent_power', 0.3, 1)
    #gatefilter.exclude_outside('corrected_cross_correlation_ratio', 0.7, 1)
    
    nyq = radar.instrument_parameters['nyquist_velocity']['data'][0]
    angle = str( round(angles[x],2))
    fig = plt.figure(figsize = [10,8])
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
           
    #  Velocity               
    display.plot_ppi('velocity', sweep = x , # axislabels = ['', 'North South distance from radar (km)'],
                    title = 'Velocity', colorbar_label='Vel. (m/s)',
                    vmin = -nyq, vmax = nyq, mask_outside = False,
                    cmap = pyart.graph.cm.NWSVel, gatefilter = gatefilter)
    display.plot_range_rings([25, 50, 75])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-120,120], ylim=[-120,120])
    
    
    plt.suptitle('\n CSU-CHIVO| Elevation: ' + angle + ' Deg.| ' + time_text, fontsize=16)
    
    return fig


def chivo_ppi_1a_vel_anomalies(filename, sweep):
    
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarDisplay(radar)
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('reflectivity')
    #gatefilter.exclude_invalid('differential_phase')
    gatefilter.exclude_outside('reflectivity', 0, 100)
    #gatefilter.exclude_outside('corrected_cross_correlation_ratio', 0.7, 1)
    
    nyq = radar.instrument_parameters['nyquist_velocity']['data'][0]
    angle = str( round(angles[x],2))
    fig = plt.figure(figsize = [10,8])
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
           
    #  Reflectivity                
    display.plot_ppi('velocity', sweep = x , # axislabels = ['', 'North South distance from radar (km)'],
                    title = 'Velocity', colorbar_label='Vel. (m/s)',
                    vmin = -nyq, vmax = nyq, mask_outside = False,
                    gatefilter = gatefilter)
    display.plot_range_rings([40, 60])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-60,-30], ylim=[-35,-5])
    
    
    plt.suptitle('\n CSU-CHIVO| Elevation: ' + angle + ' Deg.| ' + time_text, fontsize=16)
    
    return fig

def chivo_ray(filename, field, ray):
    
    radar = pyart.io.read(filename)
    display = pyart.graph.RadarDisplay(radar)
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('reflectivity')
    #gatefilter.exclude_invalid('differential_phase')
    gatefilter.exclude_outside('reflectivity', 10, 100)
    #gatefilter.exclude_outside('corrected_cross_correlation_ratio', 0.7, 1)
    
    nyq = radar.instrument_parameters['nyquist_velocity']['data'][0]
#    angle = str( round(angles[x],2))
    fig = plt.figure(figsize = [10,5])
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
           
    #  Reflectivity                
    display.plot_ray(field, ray = ray , axislabels = ['Distance from radar (km)', 'Radial Velocity (m/s)'],
                    #title = 'Velocity', colorbar_label='Vel. (m/s)',
                    #vmin = -nyq, vmax = nyq, mask_outside = False,
                    gatefilter = gatefilter)
    #display.plot_range_rings([40, 60])
    #display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[0,80], ylim=[-nyq,nyq])
    
    
    return fig

def chivo_ppi_1a_vel_unfold(filename, sweep):
    
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarDisplay(radar)
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('reflectivity')
    #gatefilter.exclude_invalid('differential_phase')
    gatefilter.exclude_outside('reflectivity', 8, 100)
    gatefilter.exclude_outside('normalized_coherent_power', 0.5, 1)
    #gatefilter.exclude_outside('cross_correlation_ratio', 0.7, 1)
    
    nyq = radar.instrument_parameters['nyquist_velocity']['data'][0]
    corr_vel = pyart.correct.dealias_region_based(
            radar, vel_field='velocity', keep_original=False, 
            gatefilter = gatefilter, nyquist_vel=nyq, centered = True)
    radar.add_field('corrected_velocity', corr_vel, replace_existing = True)
    
    angle = str( round(angles[x],2))
    fig = plt.figure(figsize = [10,8])
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
                          
    display.plot_ppi('corrected_velocity', sweep = x , # axislabels = ['', 'North South distance from radar (km)'],
                    title = 'Unfoled Velocity', colorbar_label='Vel. (m/s)',
                    vmin = -2*nyq, vmax = 2*nyq, mask_outside = False,
                    cmap = pyart.graph.cm.NWSVel)#, gatefilter = gatefilter)
    #display.plot_range_rings([25, 50, 75])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-60,60], ylim=[-60,60])
    
    
    plt.suptitle('\n CSU-CHIVO| Elevation: ' + angle + ' Deg.| ' + time_text, fontsize=16)
    
    return fig

def ppi_drops_tmp(filename, sweep):
    
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarDisplay(radar)
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('corrected_reflectivity')
    gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 0, 100)
    gatefilter.exclude_outside('corrected_cross_correlation_ratio', 0.7, 1)
    
    radar.azimuth['data'] = radar.azimuth['data'] - 13
    
    angle = str( round(angles[x],2))
    fig = plt.figure(figsize = [19,10])
    file_name = filename[len(filename)-48:len(filename)-33]
           
    #  Reflectivity                
    ax = fig.add_subplot(231)
    display.plot_ppi('corrected_reflectivity', sweep = x , ax = ax, axislabels = ['', 'North South distance from radar (km)'],
                    title = 'Corrected Reflectivity', colorbar_label='Z (dBZ)',
                    vmin = 0, vmax = 80, mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef,
                    gatefilter = gatefilter)
    display.plot_range_rings([25, 50, 75])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-80,80], ylim=[-80, 80])
    
    #   Differential Reflectivity
    ax = fig.add_subplot(232)                
    display.plot_ppi('corrected_differential_reflectivity', sweep = x , ax = ax, axislabels = ['', ''],
                    vmin = -2, vmax = 8, mask_outside = False,
                    cmap = pyart.graph.cm.RefDiff, 
                    title='Diff. Reflectivity', colorbar_label='Zdr (dB)',
                    gatefilter = gatefilter)
    display.plot_range_rings([25, 50, 75])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-80,80], ylim=[-80, 80])

    
    #   HyC
    ax = fig.add_subplot(233)
    hid_colors = ['White', 'Black','LightBlue', 'MediumBlue', 'DarkOrange', 'LightPink',
                  'Cyan', 'DarkGray', 'Lime', 'Yellow', 'Red', 'Fuchsia']
    cmaphid = colors.ListedColormap(hid_colors)
    display.plot_ppi('HydroClass', sweep = x ,ax = ax, 
                     axislabels = ['', ''],
                     cmap = cmaphid, colorbar_label='',
                     vmin = 5, vmax = 16, mask_outside = False, gatefilter = gatefilter,
                     title = 'Hydrometeor Classification')
    display.cbs[2] = adjust_fhc_colorbar_for_pyart(display.cbs[2])
    display.plot_range_rings([25, 50, 75])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-80,80], ylim=[-80, 80])
      
    #   RHOHV
    ax = fig.add_subplot(234)
    display.plot_ppi('corrected_cross_correlation_ratio', sweep = x , ax =  ax, #axislabels = ['km', 'km'],
                    vmin = 0.7, vmax = 1, mask_outside = False,
                    title='Cross Polar Correlation Ratio', colorbar_label='RhoHV', 
                    cmap = 'jet', gatefilter = gatefilter)
    display.plot_range_rings([25, 50, 75])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-80,80], ylim=[-80, 80])
    
    #   Specific Differential Phase
    ax = fig.add_subplot(235)
    display.plot_ppi('corrected_specific_differential_phase', sweep = x , ax =  ax, axislabels = ['East West distance from radar (km)', ''],
                    vmin = 0, vmax = 6, mask_outside = False,
                    title='Specific Diff. Phase', colorbar_label= 'Kdp (deg/km)',
                    gatefilter = gatefilter)
    display.plot_range_rings([25, 50, 75])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-80,80], ylim=[-80, 80])
      
    #   Differential Phase
    ax = fig.add_subplot(236)
    display.plot_ppi('corrected_differential_phase', sweep = x , ax =  ax, axislabels = ['East West distance from radar (km)', ''],
                    vmin = -180, vmax = 180, mask_outside = False,
                    title='Differential Phase', colorbar_label='Phidp (deg.)', 
                    gatefilter = gatefilter)
    display.plot_range_rings([25, 50, 75])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-80,80], ylim=[-80, 80])
    plt.suptitle('\n CSU-CHIVO| Elevation: ' + angle + ' Deg.| ' + file_name, fontsize=16)
    
    return fig



#----------------------------- CSAPR -------------------------------

def csapr_ppi_1a_ref(filename, sweep):
    
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarDisplay(radar)
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('reflectivity')
    #gatefilter.exclude_invalid('differential_phase')
    gatefilter.exclude_outside('reflectivity', -20, 100)
    gatefilter.exclude_outside('normalized_coherent_power', 0.3, 1)
    #gatefilter.exclude_outside('corrected_cross_correlation_ratio', 0.7, 1)
    
    angle = str( round(angles[x],2))
    fig = plt.figure(figsize = [10,8])
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
           
    #  Reflectivity               
    display.plot_ppi('reflectivity', sweep = x , # axislabels = ['', 'North South distance from radar (km)'],
                    title = 'Reflectivity', colorbar_label='Ref. (dBZ)',
                    vmin = 0, vmax = 60, mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef, gatefilter = gatefilter)
    #display.plot_range_rings([25, 50, 75])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-100,100], ylim=[-100,100])
    
    
    plt.suptitle('\n CSAPR | Elevation: ' + angle + ' Deg.| ' + time_text, fontsize=16)
    
    return fig

def csapr_ppi_1a_vel_unfold(filename, sweep):
    
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarDisplay(radar)
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('reflectivity')
    #gatefilter.exclude_invalid('differential_phase')
    gatefilter.exclude_outside('reflectivity', 1.5, 100)
    gatefilter.exclude_outside('normalized_coherent_power', 0.2, 1)
    #gatefilter.exclude_outside('copol_correlation_coeff', 0.7, 1)
    
    nyq = radar.instrument_parameters['nyquist_velocity']['data'][0]
    corr_vel = pyart.correct.dealias_region_based(
            radar, vel_field='mean_doppler_velocity', keep_original=False, 
            gatefilter = gatefilter, nyquist_vel=nyq, centered = True)
    radar.add_field('corrected_velocity', corr_vel, replace_existing = True)
    
    angle = str( round(angles[x],2))
    fig = plt.figure(figsize = [10,8])
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
                          
    display.plot_ppi('corrected_velocity', sweep = x , # axislabels = ['', 'North South distance from radar (km)'],
                    title = 'Unfoled Velocity', colorbar_label='Vel. (m/s)',
                    vmin = -2*nyq, vmax = 2*nyq, mask_outside = False,
                    cmap = pyart.graph.cm.NWSVel, gatefilter = gatefilter)
    #display.plot_range_rings([25, 50, 75])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-100,100], ylim=[-100,100])
    
    
    plt.suptitle('\n CSAPR| Elevation: ' + angle + ' Deg.| ' + time_text, fontsize=16)
    
    return fig



def csapr_ppi_drops(filename, sweep):
    
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarDisplay(radar)
    radar.azimuth['data'] = radar.azimuth['data']
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('corrected_reflectivity')
    gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 10, 100)
    #gatefilter.exclude_outside('corrected_cross_correlation_ratio', 0.7, 1)
    
    
    #angle = str( round(sweep*30,2))
    angle = str( round(angles[x],2))
    fig = plt.figure(figsize = [19,10])
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
           
    #  Reflectivity                
    ax = fig.add_subplot(231)
    display.plot_ppi('corrected_reflectivity', sweep = x , ax = ax, axislabels = ['', 'North South distance from radar (km)'],
                    title = 'Corrected Reflectivity', colorbar_label='Z (dBZ)',
                    vmin = 0, vmax = 80, mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef,
                    gatefilter = gatefilter)
    display.plot_range_rings([50, 100, 150])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-100,100], ylim=[-100, 100])
    
    #   Differential Reflectivity
    ax = fig.add_subplot(232)                
    display.plot_ppi('corrected_differential_reflectivity', sweep = x , ax = ax, axislabels = ['', ''],
                    vmin = -2, vmax = 8, mask_outside = False,
                    cmap = pyart.graph.cm.RefDiff, 
                    title='Diff. Reflectivity', colorbar_label='Zdr (dB)',
                    gatefilter = gatefilter)
    display.plot_range_rings([50, 100, 150])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-100,100], ylim=[-100, 100])

    
    #   HyC
    ax = fig.add_subplot(233)
    hid_colors = ['White', 'Black','LightBlue', 'MediumBlue', 'DarkOrange', 'LightPink',
                  'Cyan', 'DarkGray', 'Lime', 'Yellow', 'Red', 'Fuchsia']
    cmaphid = colors.ListedColormap(hid_colors)
    display.plot_ppi('HydroClass', sweep = x ,ax = ax, 
                     axislabels = ['', ''],
                     cmap = cmaphid, colorbar_label='',
                     vmin = 5, vmax = 16, mask_outside = False, gatefilter = gatefilter,
                     title = 'Hydrometeor Classification')
    display.cbs[2] = adjust_fhc_colorbar_for_pyart(display.cbs[2])
    display.plot_range_rings([50, 100, 150])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-100,100], ylim=[-100, 100])
      
    #   RHOHV
    ax = fig.add_subplot(234)
    display.plot_ppi('corrected_cross_correlation_ratio', sweep = x , ax =  ax, #axislabels = ['km', 'km'],
                    vmin = 0.7, vmax = 1, mask_outside = False,
                    title='Copolar Correlation Ratio', colorbar_label='RhoHV', 
                    cmap = 'jet', gatefilter = gatefilter)
    display.plot_range_rings([50, 100, 150])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-100,100], ylim=[-100, 100])
    
    #   Specific Differential Phase
    ax = fig.add_subplot(235)
    display.plot_ppi('corrected_specific_differential_phase', sweep = x , ax =  ax, axislabels = ['East West distance from radar (km)', ''],
                    vmin = 0, vmax = 6, mask_outside = False,
                    title='Specific Diff. Phase', colorbar_label= 'Kdp (deg/km)',
                    gatefilter = gatefilter)
    display.plot_range_rings([50, 100, 150])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-100,100], ylim=[-100, 100])
      
    #   Differential Phase
    ax = fig.add_subplot(236)
    display.plot_ppi('corrected_differential_phase', sweep = x , ax =  ax, axislabels = ['East West distance from radar (km)', ''],
                    vmin = -180, vmax = 180, mask_outside = False,
                    title='Differential Phase', colorbar_label='Phidp (deg.)', 
                    gatefilter = gatefilter)
    display.plot_range_rings([50, 100, 150])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-100,100], ylim=[-100, 100])
    plt.suptitle('\n CSAPR| Elevation: ' + angle + ' Deg.| ' + time_text, fontsize=16)
    
    return fig

def csapr_ppi_drops_region(filename, sweep):
    
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarDisplay(radar)
    radar.azimuth['data'] = radar.azimuth['data']
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('corrected_reflectivity')
    gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 0, 100)
    gatefilter.exclude_outside('corrected_cross_correlation_ratio', 0.7, 1)
    
    
    #angle = str( round(sweep*30,2))
    angle = str( round(angles[x],2))
    fig = plt.figure(figsize = [19,10])
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S') #time_start.strftime('%Y-%m-%dT%H:%M:%S')
           
    #  Reflectivity                
    ax = fig.add_subplot(231)
    display.plot_ppi('corrected_reflectivity', sweep = x , ax = ax, axislabels = ['', 'North South distance from radar (km)'],
                    title = 'Corrected Reflectivity', colorbar_label='Z (dBZ)',
                    vmin = 0, vmax = 80, mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef,
                    gatefilter = gatefilter)
    display.plot_range_rings([25, 50, 75])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-80,80], ylim=[-80,80])
    
    #   Differential Reflectivity
    ax = fig.add_subplot(232)                
    display.plot_ppi('corrected_differential_reflectivity', sweep = x , ax = ax, axislabels = ['', ''],
                    vmin = -2, vmax = 8, mask_outside = False,
                    cmap = pyart.graph.cm.RefDiff, 
                    title='Diff. Reflectivity', colorbar_label='Zdr (dB)',
                    gatefilter = gatefilter)
    display.plot_range_rings([25, 50, 75])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-80,80], ylim=[-80,80])

    
    #   HyC
    ax = fig.add_subplot(233)
    hid_colors = ['White', 'Black','LightBlue', 'MediumBlue', 'DarkOrange', 'LightPink',
                  'Cyan', 'DarkGray', 'Lime', 'Yellow', 'Red', 'Fuchsia']
    cmaphid = colors.ListedColormap(hid_colors)
    display.plot_ppi('HydroClass', sweep = x ,ax = ax, 
                     axislabels = ['', ''],
                     cmap = cmaphid, colorbar_label='',
                     vmin = 5, vmax = 16, mask_outside = False, gatefilter = gatefilter,
                     title = 'Hydrometeor Classification')
    display.cbs[2] = adjust_fhc_colorbar_for_pyart(display.cbs[2])
    display.plot_range_rings([25, 50, 75])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-80,80], ylim=[-80,80])
      
    #   RHOHV
    ax = fig.add_subplot(234)
    display.plot_ppi('corrected_cross_correlation_ratio', sweep = x , ax =  ax, #axislabels = ['km', 'km'],
                    vmin = 0.7, vmax = 1, mask_outside = False,
                    title='Copolar Correlation Ratio', colorbar_label='RhoHV', 
                    cmap = 'jet', gatefilter = gatefilter)
    display.plot_range_rings([25, 50, 75])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-80,80], ylim=[-80,80])
    
    #   Specific Differential Phase
    ax = fig.add_subplot(235)
    display.plot_ppi('corrected_specific_differential_phase', sweep = x , ax =  ax, axislabels = ['East West distance from radar (km)', ''],
                    vmin = 0, vmax = 6, mask_outside = False,
                    title='Specific Diff. Phase', colorbar_label= 'Kdp (deg/km)',
                    gatefilter = gatefilter)
    display.plot_range_rings([25, 50, 75])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-80,80], ylim=[-80,80])
      
    #   Differential Phase
    ax = fig.add_subplot(236)
    display.plot_ppi('corrected_differential_phase', sweep = x , ax =  ax, axislabels = ['East West distance from radar (km)', ''],
                    vmin = -180, vmax = 180, mask_outside = False,
                    title='Differential Phase', colorbar_label='Phidp (deg.)', 
                    gatefilter = gatefilter)
    display.plot_range_rings([25, 50, 75])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-80,80], ylim=[-80,80])
    plt.suptitle('\n CSAPR| Elevation: ' + angle + ' Deg.| ' + time_text, fontsize=16)
    
    return fig


def csapr_ppi_drops_region_field(filename, sweep, field_name):
    
    corrected_field_name = {
            'reflectivity': 'corrected_reflectivity',
            'differential_reflectivity': 'corrected_differential_reflectivity',
            'specific_differential_phase': 'corrected_specific_differential_phase'
            }
    
    colorMap = {
            'reflectivity': pyart.graph.cm.NWSRef,
            'differential_reflectivity': pyart.graph.cm.RefDiff,
            'rhohv': 'jet',
            'specific_differential_phase': pyart.graph.cm.Theodore16            
            
            }
    minVal = {
            'reflectivity': 0,
            'differential_reflectivity': -2,
            'specific_differential_phase': -1
            }
    maxVal = {
            'reflectivity': 80,
            'differential_reflectivity': 8,
            'specific_differential_phase': 6
            }
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarDisplay(radar)#, shift = (-52000, -52000))
    radar.fields['corrected_differential_reflectivity']['data'] = radar.fields['corrected_differential_reflectivity']['data'] #- 0.7
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('corrected_reflectivity')
    gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 18, 100)
    #gatefilter.exclude_outside('corrected_copol_correlation_ratio', 0.7, 1)
    gatefilter.exclude_inside('normalized_coherent_power', 0, 0.2)
    
    
    angle = str( round(angles[x],2))
    fig = plt.figure(figsize = [5,4])
    
    # time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    # time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
    
    # time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    # time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S') #time_start.strftime('%Y-%m-%dT%H:%M:%S')
           

    display.plot_ppi(corrected_field_name[field_name], sweep = x , title = '',
                    vmin = minVal[field_name], vmax = maxVal[field_name], mask_outside = False,
                    cmap = colorMap[field_name],  gatefilter = gatefilter)
   
    display.plot_range_rings([10, 20, 30, 40, 50])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    #display.set_limits(xlim=[-55, 20], ylim=[-55, 20])
    display.set_limits(xlim=[-20, 80], ylim=[-20, 80])
    angle = str( round(angles[x],2))    
    #plt.suptitle('CSARP-2| Azimuth: ' + angle + ' Deg.| ' + time_text, fontsize=14)
    
    return fig


def csapr_ppi_drops_region_ref(filename, sweep, azimuth_offset):
    
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarDisplay(radar)
    radar.azimuth['data'] = radar.azimuth['data'] - azimuth_offset#13
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    #gatefilter.exclude_invalid('corrected_reflectivity')
    #gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 13, 100)
    #gatefilter.exclude_outside('corrected_cross_correlation_ratio', 0.7, 1)
    
    
    angle = str( round(angles[x],2))
    fig = plt.figure(figsize = [10,8])
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
           
    #  Reflectivity                
    display.plot_ppi('corrected_reflectivity', sweep = x , # axislabels = ['', 'North South distance from radar (km)'],
                    title = 'Reflectivity', colorbar_label='Z (dBZ)',
                    vmin = 0, vmax = 80, mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef,
                    gatefilter = gatefilter)
    display.plot_range_rings([25, 50, 75])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-80,80], ylim=[-80,80])
    
    
    plt.suptitle('\n CSAPR | Elevation: ' + angle + ' Deg.| ' + time_text, fontsize=16)
    
    return fig


def csapr_ppi_1a_vel(filename, sweep):
    
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarDisplay(radar)
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('mean_doppler_velocity')
    gatefilter.exclude_invalid('reflectivity')
    gatefilter.exclude_outside('reflectivity', 10, 100)
    gatefilter.exclude_outside('normalized_coherent_power', 0.2, 1)
    gatefilter.exclude_outside('copol_correlation_coeff', 0.7, 1)
    
    nyq = radar.instrument_parameters['nyquist_velocity']['data'][0]
    angle = str( round(angles[x],2))
    fig = plt.figure(figsize = [10,8])
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
                          
    display.plot_ppi('mean_doppler_velocity', sweep = x , # axislabels = ['', 'North South distance from radar (km)'],
                    title = 'Velocity', colorbar_label='Vel. (m/s)',
                    vmin = -nyq, vmax = nyq, mask_outside = False,
                    cmap = pyart.graph.cm.NWSVel,#cmap = pyart.graph.cm.BuDRd18,
                    gatefilter = gatefilter)
    display.plot_range_rings([25, 50, 75])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-100,100], ylim=[-100,100])
    
    
    plt.suptitle('\n CSAPR| Elevation: ' + angle + ' Deg.| ' + time_text, fontsize=16)
    
    return fig


def csapr_ppi(filename, sweep):
    
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarDisplay(radar)
    radar.azimuth['data'] = radar.azimuth['data']
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('reflectivity')
    #gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('reflectivity', 10, 100)
    #gatefilter.exclude_outside('corrected_cross_correlation_ratio', 0.7, 1)
    
    
    #angle = str( round(sweep*30,2))
    angle = str( round(angles[x],2))
    fig = plt.figure(figsize = [19,10])
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
           
    #  Reflectivity                
    ax = fig.add_subplot(231)
    display.plot_ppi('reflectivity', sweep = x , ax = ax, axislabels = ['', 'North South distance from radar (km)'],
                    title = 'Reflectivity', colorbar_label='Z (dBZ)',
                    vmin = 0, vmax = 70, mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef,
                    gatefilter = gatefilter)
    display.plot_range_rings([50, 100, 150])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-100,100], ylim=[-100, 100])
    
    #   Differential Reflectivity
    ax = fig.add_subplot(232)                
    display.plot_ppi('differential_reflectivity', sweep = x , ax = ax, axislabels = ['', ''],
                    vmin = -2, vmax = 8, mask_outside = False,
                    cmap = pyart.graph.cm.RefDiff, 
                    title='Diff. Reflectivity', colorbar_label='Zdr (dB)',
                    gatefilter = gatefilter)
    display.plot_range_rings([50, 100, 150])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-100,100], ylim=[-100, 100])

     
    #   RHOHV
    ax = fig.add_subplot(234)
    display.plot_ppi('copol_correlation_coeff', sweep = x , ax =  ax, #axislabels = ['km', 'km'],
                    vmin = 0.7, vmax = 1, mask_outside = False,
                    title='Copolar Correlation Ratio', colorbar_label='RhoHV', 
                    cmap = 'jet', gatefilter = gatefilter)
    display.plot_range_rings([50, 100, 150])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-100,100], ylim=[-100, 100])
    
    #   Specific Differential Phase
    ax = fig.add_subplot(235)
    display.plot_ppi('specific_differential_phase', sweep = x , ax =  ax, axislabels = ['East West distance from radar (km)', ''],
                    vmin = 0, vmax = 6, mask_outside = False,
                    title='Specific Diff. Phase', colorbar_label= 'Kdp (deg/km)',
                    gatefilter = gatefilter)
    display.plot_range_rings([50, 100, 150])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-100,100], ylim=[-100, 100])
      
    #   Differential Phase
    ax = fig.add_subplot(236)
    display.plot_ppi('differential_phase', sweep = x , ax =  ax, axislabels = ['East West distance from radar (km)', ''],
                    vmin = -180, vmax = 180, mask_outside = False,
                    title='Differential Phase', colorbar_label='Phidp (deg.)', 
                    gatefilter = gatefilter)
    display.plot_range_rings([50, 100, 150])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-100,100], ylim=[-100, 100])
    plt.suptitle('\n CSAPR| Elevation: ' + angle + ' Deg.| ' + time_text, fontsize=16)
    
    return fig

def csapr_rhi_drops(filename, sweep, xlim):
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarMapDisplay(radar)
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
    
    gatefilter = pyart.filters.GateFilter(radar)
    #gatefilter.exclude_transition()
    #gatefilter.exclude_invalid('corrected_reflectivity')
    #gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', -15, 100)
    #gatefilter.exclude_outside('corrected_cross_correlation_ratio', 0.7, 1)
    
    angle = str( round(sweep*30,2)) # azimuth correction of 13 deg
    fig = plt.figure(figsize = [15,10])
    file_name = filename[len(filename)-46:len(filename)-31]
    
    #  Reflectivity                
    ax = fig.add_subplot(321)
    display.plot_rhi('corrected_reflectivity', sweep = x , ax = ax, title = '',
                     axislabels = ['', 'Distance Above radar (km)'],
                    vmin = -15, vmax = 80, mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef, colorbar_label =  'Corrected Reflectivity (dBZ)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, xlim])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   Differential Reflectivity
    ax = fig.add_subplot(322)                
    display.plot_rhi('corrected_differential_reflectivity', sweep = x , ax = ax,
                     axislabels = ['', ''],
                    vmin = -2, vmax = 8, mask_outside = False,
                    cmap = pyart.graph.cm.RefDiff, title='', colorbar_label='Diff. Reflectivity (dB)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, xlim])
    display.plot_grid_lines(ax=None, col='k', ls=':')

    #   KDP
    ax = fig.add_subplot(323)
    display.plot_rhi('corrected_specific_differential_phase', sweep = x , ax =  ax, 
                    axislabels = ['', 'Distance Above radar (km)'],
                    vmin = 0, vmax = 6, mask_outside = False,
                    title='', colorbar_label= 'Specific Diff. Phase (deg/km)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, xlim])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   PHIDP
    ax = fig.add_subplot(324)
    display.plot_rhi('corrected_differential_phase', sweep = x , ax =  ax,
                    axislabels = ['', ' '],
                    vmin = -180, vmax = 180, mask_outside = False,
                     title='', colorbar_label='Corrected Diff. Phase (°)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, xlim])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   RhoHV
    ax = fig.add_subplot(325)
    display.plot_rhi('corrected_cross_correlation_ratio', sweep = x , ax =  ax,
                     vmin = 0.7, vmax = 1, mask_outside = False, cmap = 'jet',
                     title='', colorbar_label='Copolar Correlcation',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, xlim])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   HyC
    
    hid_colors = ['White', 'Black','LightBlue', 'MediumBlue', 'DarkOrange', 'LightPink',
                  'Cyan', 'DarkGray', 'Lime', 'Yellow', 'Red', 'Fuchsia']
    cmaphid = colors.ListedColormap(hid_colors)
    ax = fig.add_subplot(326)
    display.plot_rhi('HydroClass', sweep = x ,ax = ax, 
                     axislabels = ['Distance from radar (km)', ' '],
                     cmap = cmaphid, colorbar_label='',
                     vmin = 5, vmax = 16, mask_outside = False, gatefilter = gatefilter,
                     title = ' ')
    display.set_limits(ylim=[0, 21])
    display.cbs[5] = adjust_fhc_colorbar_for_pyart(display.cbs[5])
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, xlim])
    display.plot_grid_lines(ax=None, col='k', ls=':')
        
    plt.suptitle('\n CSAPR| Azimuth: ' + angle + ' Deg.| ' + time_text, fontsize=16)
    print('check line 385 at ' + '/top/students/GRAD/ECE/idariash/home/anaconda3/lib/python3.6/site-packages/pyart/core/radar.py')
    return fig

def csapr_rhi_drops_ref(filename, sweep, xlim):
    radar = pyart.io.read(filename)
    x = sweep
    display = pyart.graph.RadarMapDisplay(radar)
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
    
    gatefilter = pyart.filters.GateFilter(radar)
    #gatefilter.exclude_transition()
    #gatefilter.exclude_invalid('corrected_reflectivity')
    #gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', -15, 100)
    #gatefilter.exclude_outside('corrected_cross_correlation_ratio', 0.7, 1)
    
    angle = str( round(sweep*30,2)) # azimuth correction of 13 deg
    fig = plt.figure(figsize = [10,4])
    
    #  Reflectivity                
    display.plot_rhi('corrected_reflectivity', sweep = x, title = '',
                     #axislabels = ['', 'Distance Above radar (km)'],
                    vmin = -15, vmax = 80, mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef, colorbar_label =  'Corrected Reflectivity (dBZ)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, xlim])
    display.plot_grid_lines(ax=None, col='k', ls=':')
        
    plt.suptitle('CSAPR| Azimuth: ' + angle + ' Deg.| ' + time_text, fontsize=14)
    print('check line 385 at ' + '/top/students/GRAD/ECE/idariash/home/anaconda3/lib/python3.6/site-packages/pyart/core/radar.py')
    return fig


def csapr_rhi_drops_4dd(filename_drops, filename_sigmet ,sweep, xlim):
    radar = pyart.io.read(filename_drops)
    x = sweep
    angles = radar.fixed_angle['data'] # - 13
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
    display = pyart.graph.RadarMapDisplay(radar)
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    #gatefilter.exclude_invalid('corrected_reflectivity')
    #gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 0, 100)
    #gatefilter.exclude_outside('corrected_cross_correlation_ratio', 0.6, 1)
    
    angle = str( round(x*30,2))
    fig = plt.figure(figsize = [15,10])
    file_name = filename_drops[len(filename_drops)-46:len(filename_drops)-31]
    
    #  Reflectivity                
    ax = fig.add_subplot(321)
    display.plot_rhi('corrected_reflectivity', sweep = x , ax = ax, title = '',
                     axislabels = ['', 'Distance Above radar (km)'],
                    vmin = 0, vmax = 80, mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef, colorbar_label =  'Corrected Reflectivity (dbZ)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, xlim])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   Differential Reflectivity
    ax = fig.add_subplot(322)                
    display.plot_rhi('corrected_differential_reflectivity', sweep = x , ax = ax,
                     axislabels = ['', ''],
                    vmin = -2, vmax = 8, mask_outside = False,
                    cmap = pyart.graph.cm.RefDiff, title='', colorbar_label='Diff. Reflectivity (dB)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, xlim])
    display.plot_grid_lines(ax=None, col='k', ls=':')

    #   KDP
    ax = fig.add_subplot(323)
    display.plot_rhi('corrected_specific_differential_phase', sweep = x , ax =  ax, 
                    axislabels = ['', 'Distance Above radar (km)'],
                    vmin = -1, vmax = 6, mask_outside = False,
#                    vmin = 0, vmax = 6, mask_outside = False,
                    title='', colorbar_label= 'Specific Diff. Phase (deg/km)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, xlim])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   Velocity
    last_radar = pyart.io.read(filename_sigmet)
    radar_sigmet = pyart.io.read(filename_sigmet)
    display_sigmet = pyart.graph.RadarMapDisplay(last_radar)
    
    gatefilter_sigmet = pyart.filters.GateFilter(last_radar)
    gatefilter_sigmet.exclude_transition()
    gatefilter_sigmet.exclude_invalid('mean_doppler_velocity')
    gatefilter_sigmet.exclude_invalid('reflectivity')
    gatefilter_sigmet.exclude_outside('reflectivity', 20, 100)
    #gatefilter_sigmet.exclude_outside('normalized_coherent_power', 0.2, 1)
    #gatefilter_sigmet.exclude_outside('copol_correlation_coeff', 0.7, 1)
    
    nyq = radar_sigmet.instrument_parameters['nyquist_velocity']['data'][0]
    corr_vel = pyart.correct.dealias_region_based(
            last_radar, vel_field='mean_doppler_velocity', keep_original=False, 
            gatefilter = gatefilter_sigmet, nyquist_vel=nyq, centered = True)
    last_radar.add_field('corrected_velocity', corr_vel, replace_existing = True)
#    
#    dealias_data = pyart.correct.dealias_fourdd(vel_field='mean_doppler_velocity',
#            radar_sigmet, last_radar = last_radar, gatefilter = gatefilter_sigmet)
#    radar_sigmet.add_field('corrected_velocity', dealias_data, replace_existing = True)
    
    ax = fig.add_subplot(324)
    display_sigmet.plot_rhi('corrected_velocity', sweep = x , ax =  ax,
                    axislabels = ['', ' '],
                    vmin = -2*nyq, vmax = nyq, mask_outside = False,
                    title='', colorbar_label= 'Unlfolded Velocity (m/s)',
                    cmap = pyart.graph.cm.BuDRd18,
                    gatefilter = gatefilter)
    display_sigmet.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, xlim])
    display_sigmet.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   RhoHV
    ax = fig.add_subplot(325)
    display.plot_rhi('corrected_cross_correlation_ratio', sweep = x , ax =  ax,
                     vmin = 0.7, vmax = 1, mask_outside = False, cmap = 'jet',
                     title='', colorbar_label='Copolar Correlcation',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, xlim])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   HyC
    
    hid_colors = ['White', 'Black','LightBlue', 'MediumBlue', 'DarkOrange', 'LightPink',
                  'Cyan', 'DarkGray', 'Lime', 'Yellow', 'Red', 'Fuchsia']
    cmaphid = colors.ListedColormap(hid_colors)
    ax = fig.add_subplot(326)
    display.plot_rhi('HydroClass', sweep = x ,ax = ax, 
                     axislabels = ['Distance from radar (km)', ' '],
                     cmap = cmaphid, colorbar_label='',
                     vmin = 5, vmax = 16, mask_outside = False, gatefilter = gatefilter,
                     title = ' ')
    display.set_limits(ylim=[0, 21])
    display.cbs[4] = adjust_fhc_colorbar_for_pyart(display.cbs[4])
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, xlim])
    display.plot_grid_lines(ax=None, col='k', ls=':')
        
    plt.suptitle('\n CSAPR | Azimuth: ' + angle + ' Deg.| ' + time_text, fontsize=16)
    return fig

def csapr_rhi(filename, sweep, xlim):
    radar = pyart.io.read(filename)
    x = sweep
    display = pyart.graph.RadarMapDisplay(radar)
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('reflectivity')
    #gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('reflectivity', 0, 100)
    #gatefilter.exclude_outside('corrected_cross_correlation_ratio', 0.7, 1)
    
    radar.fields['differential_reflectivity']['data'] = radar.fields['differential_reflectivity']['data'] - 0.51 #(-3.5)
    
    angle = str( round(sweep*30,2))
    fig = plt.figure(figsize = [15,10])
    
    #  Reflectivity                
    ax = fig.add_subplot(321)
    display.plot_rhi('reflectivity', sweep = x , ax = ax, title = '',
                     axislabels = ['', 'Distance Above radar (km)'],
                    vmin = -15, vmax = 80, mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef, colorbar_label =  'Reflectivity (dBZ)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, xlim])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   Differential Reflectivity
    ax = fig.add_subplot(322)                
    display.plot_rhi('differential_reflectivity', sweep = x , ax = ax,
                     axislabels = ['', ''],
                    vmin = -2, vmax = 8, mask_outside = False,
                    cmap = pyart.graph.cm.RefDiff, title='', colorbar_label='Diff. Reflectivity (dB)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, xlim])
    display.plot_grid_lines(ax=None, col='k', ls=':')

    #   KDP
    ax = fig.add_subplot(323)
    display.plot_rhi('specific_differential_phase', sweep = x , ax =  ax, 
                    axislabels = ['', 'Distance Above radar (km)'],
                    vmin = 0, vmax = 6, mask_outside = False,
                    title='', colorbar_label= 'Specific Diff. Phase (deg/km)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, xlim])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   PHIDP
    ax = fig.add_subplot(324)
    display.plot_rhi('differential_phase', sweep = x , ax =  ax,
                    axislabels = ['', ' '],
                    vmin = -180, vmax = 180, mask_outside = False,
                     title='', colorbar_label='Corrected Diff. Phase (°)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, xlim])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   RhoHV
    ax = fig.add_subplot(325)
    display.plot_rhi('copol_correlation_coeff', sweep = x , ax =  ax,
                     vmin = 0.7, vmax = 1, mask_outside = False, cmap = 'jet',
                     title='', colorbar_label='Copolar Correlcation',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, xlim])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
        
    plt.suptitle('\n CSU-CHIVO| Azimuth: ' + angle + ' Deg.| ' + time_text, fontsize=16)
    return fig

def csapr_rhi_from_ppi(filename, azimuth):
    radar = pyart.io.read(filename)
    x = 0
    angles = radar.fixed_angle['data'];
    xsect = pyart.util.cross_section_ppi(radar, [azimuth])
    display = pyart.graph.RadarMapDisplay(xsect)
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    #gatefilter.exclude_invalid('corrected_reflectivity')
    #gatefilter.exclude_invalid('corrected_differential_phase')
    #gatefilter.exclude_outside('corrected_reflectivity', -15, 100)
    #gatefilter.exclude_outside('corrected_cross_correlation_ratio', 0.7, 1)
    
    angle = str( round(angles[x],2))
    fig = plt.figure(figsize = [15,10])
    file_name = filename[len(filename)-46:len(filename)-31]
    
    #  Reflectivity                
    ax = fig.add_subplot(321)
    display.plot_rhi('corrected_reflectivity', sweep = x , ax = ax, title = '',
                     axislabels = ['', 'Distance Above radar (km)'],
                    vmin = -15, vmax = 80, mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef, colorbar_label =  'Corrected Reflectivity (dBZ)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, 100])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   Differential Reflectivity
    ax = fig.add_subplot(322)                
    display.plot_rhi('corrected_differential_reflectivity', sweep = x , ax = ax,
                     axislabels = ['', ''],
                    vmin = -2, vmax = 8, mask_outside = False,
                    cmap = pyart.graph.cm.RefDiff, title='', colorbar_label='Diff. Reflectivity (dB)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, 100])
    display.plot_grid_lines(ax=None, col='k', ls=':')

    #   KDP
    ax = fig.add_subplot(323)
    display.plot_rhi('corrected_specific_differential_phase', sweep = x , ax =  ax, 
                    axislabels = ['', 'Distance Above radar (km)'],
                    vmin = 0, vmax = 6, mask_outside = False,
                    title='', colorbar_label= 'Specific Diff. Phase (deg/km)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, 100])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   PHIDP
    ax = fig.add_subplot(324)
    display.plot_rhi('corrected_differential_phase', sweep = x , ax =  ax,
                    axislabels = ['', ' '],
                    vmin = -180, vmax = 180, mask_outside = False,
                     title='', colorbar_label='Corrected Diff. Phase (°)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, 100])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   RhoHV
    ax = fig.add_subplot(325)
    display.plot_rhi('corrected_cross_correlation_ratio', sweep = x , ax =  ax,
                     vmin = 0.7, vmax = 1, mask_outside = False, cmap = 'jet',
                     title='', colorbar_label='Copolar Correlcation',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, 100])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   HyC
    
    hid_colors = ['White', 'Black','LightBlue', 'MediumBlue', 'DarkOrange', 'LightPink',
                  'Cyan', 'DarkGray', 'Lime', 'Yellow', 'Red', 'Fuchsia']
    cmaphid = colors.ListedColormap(hid_colors)
    ax = fig.add_subplot(326)
    display.plot_rhi('HydroClass', sweep = x ,ax = ax, 
                     axislabels = ['Distance from radar (km)', ' '],
                     cmap = cmaphid, colorbar_label='',
                     vmin = 5, vmax = 16, mask_outside = False, gatefilter = gatefilter,
                     title = ' ')
    display.set_limits(ylim=[0, 21])
    display.cbs[5] = adjust_fhc_colorbar_for_pyart(display.cbs[5])
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[0, 100])
    display.plot_grid_lines(ax=None, col='k', ls=':')
        
    plt.suptitle('\n CSAPR| Azimuth: ' + str(azimuth) + ' Deg.| ' + time_text, fontsize=16)
    return fig

#------------RMA1------------------
    
def rma_ppi_1a_vel(filename, sweep):
    
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarDisplay(radar)
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('DBZH')
    #gatefilter.exclude_invalid('differential_phase')
    gatefilter.exclude_outside('DBZH', 8, 100)
    #gatefilter.exclude_outside('normalized_coherent_power', 0.3, 1)
    #gatefilter.exclude_outside('corrected_cross_correlation_ratio', 0.7, 1)
    
    nyq = np.amax(radar.fields['VRAD']['data'])

    angle = str( round(angles[x],2))
    fig = plt.figure(figsize = [10,8])
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
           
    #  Velocity               
    display.plot_ppi('VRAD', sweep = x , # axislabels = ['', 'North South distance from radar (km)'],
                    title = 'Velocity', colorbar_label='Vel. (m/s)',
                    vmin = -nyq, vmax = nyq, mask_outside = False,
                    cmap = pyart.graph.cm.NWSVel, gatefilter = gatefilter)
    #display.plot_range_rings([25, 50, 75])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-100,100], ylim=[-100,100])
    
    
    plt.suptitle('\n RMA-1 Elevation: ' + angle + ' Deg.| ' + time_text, fontsize=16)
    
    return fig

def rma_ppi_1a_vel_unfold(filename, sweep):
    
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarDisplay(radar)
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('DBZH')
    #gatefilter.exclude_invalid('differential_phase')
    gatefilter.exclude_outside('DBZH', 8, 100)
    #gatefilter.exclude_outside('normalized_coherent_power', 0.5, 1)
    #gatefilter.exclude_outside('cross_correlation_ratio', 0.7, 1)
    
    nyq = np.amax(radar.fields['VRAD']['data'])
    corr_vel = pyart.correct.dealias_region_based(
            radar, vel_field='VRAD', keep_original=False, 
            gatefilter = gatefilter, nyquist_vel=nyq, centered = True)
    radar.add_field('corrected_velocity', corr_vel, replace_existing = True)
    
    angle = str( round(angles[x],2))
    fig = plt.figure(figsize = [10,8])
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
                          
    display.plot_ppi('corrected_velocity', sweep = x , # axislabels = ['', 'North South distance from radar (km)'],
                    title = 'Unfoled Velocity', colorbar_label='Vel. (m/s)',
                    vmin = -3*nyq, vmax = 3*nyq, mask_outside = False,
                    cmap = pyart.graph.cm.NWSVel)#, gatefilter = gatefilter)
    #display.plot_range_rings([25, 50, 75])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-100,100], ylim=[-100,100])
    
    
    plt.suptitle('\n CSU-CHIVO| Elevation: ' + angle + ' Deg.| ' + time_text, fontsize=16)
    
    return fig


#----------------- Functions for data quality
    
def chivo_ppi_1a_ref(filename, sweep):
    
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarDisplay(radar)
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('reflectivity')
    #gatefilter.exclude_invalid('differential_phase')
    gatefilter.exclude_outside('reflectivity', 0, 100)
    gatefilter.exclude_outside('normalized_coherent_power', 0.5, 1)
    gatefilter.exclude_outside('cross_correlation_ratio', 0.95, 1)
    
    angle = str( round(angles[x],2))
    fig = plt.figure(figsize = [10,8])
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
           
    #  Reflectivity               
    display.plot_ppi('reflectivity', sweep = x , # axislabels = ['', 'North South distance from radar (km)'],
                    title = 'Reflectivity', colorbar_label='Ref. (dBZ)',
                    vmin = 0, vmax = 70, mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef, gatefilter = gatefilter) #pyart.graph.cm.NWSRef
    #display.plot_range_rings([25, 50, 75])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-100,100], ylim=[-100,100])
    
    
    plt.suptitle('\n CSU-CHIVO| Elevation: ' + angle + ' Deg.| ' + time_text, fontsize=16)
    
    return fig





def chivo_ppi_1b_ref(filename, sweep):
    
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarDisplay(radar)
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('corrected_reflectivity')
    #gatefilter.exclude_invalid('differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 0, 100)
    #gatefilter.exclude_outside('normalized_coherent_power', 0.3, 1)
    gatefilter.exclude_outside('corrected_copol_correlation_ratio', 0.7, 1)
    
    angle = str( round(angles[x],2))
    fig = plt.figure(figsize = [10,8])
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
           
    #  Reflectivity               
    display.plot_ppi('corrected_reflectivity', sweep = x , # axislabels = ['', 'North South distance from radar (km)'],
                    title = 'Reflectivity', colorbar_label='Ref. (dBZ)',
                    vmin = 0, vmax = 70, mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef, gatefilter = gatefilter)
    #display.plot_range_rings([25, 50, 75])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-100,100], ylim=[-100,100])
    
    
    plt.suptitle('\n CSU-CHIVO| Elevation: ' + angle + ' Deg.| ' + time_text, fontsize=16)
    
    return fig



def chivo_ppi_1a(filename, sweep, radar_name):
    
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarDisplay(radar)
    radar.azimuth['data'] = radar.azimuth['data']
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('reflectivity')
    #gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('reflectivity', 10, 100)
    #gatefilter.exclude_outside('corrected_copol_correlation_ratio', 0.7, 1)
    
    
    #angle = str( round(sweep*30,2))
    angle = str( round(angles[x],2))
    fig = plt.figure(figsize = [19,10])
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S') #time_start.strftime('%Y-%m-%dT%H:%M:%S')
           
    #  Reflectivity                
    ax = fig.add_subplot(231)
    display.plot_ppi('reflectivity', sweep = x , ax = ax, axislabels = ['', 'North South distance from radar (km)'],
                    title = 'Reflectivity', colorbar_label='Z (dBZ)',
                    vmin = 0, vmax = 70, mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef,
                    gatefilter = gatefilter)
    display.plot_range_rings([25, 50, 75, 100])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-100,100], ylim=[-100,100])
    
    #   Differential Reflectivity
    ax = fig.add_subplot(232)                
    display.plot_ppi('differential_reflectivity', sweep = x , ax = ax, axislabels = ['', ''],
                    vmin = -2, vmax = 8, mask_outside = False,
                    cmap = pyart.graph.cm.RefDiff, 
                    title='Diff. Reflectivity', colorbar_label='Zdr (dB)',
                    gatefilter = gatefilter)
    display.plot_range_rings([25, 50, 75, 100])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-100,100], ylim=[-100,100])

    
    #   HyC
    # ax = fig.add_subplot(233)
    # hid_colors = ['White', 'Black','LightBlue', 'MediumBlue', 'DarkOrange', 'LightPink',
    #               'Cyan', 'DarkGray', 'Lime', 'Yellow', 'Red', 'Fuchsia']
    # cmaphid = colors.ListedColormap(hid_colors)
    # display.plot_ppi('HydroClass', sweep = x ,ax = ax, 
    #                  axislabels = ['', ''],
    #                  cmap = cmaphid, colorbar_label='',
    #                  vmin = 5, vmax = 16, mask_outside = False, gatefilter = gatefilter,
    #                  title = 'Hydrometeor Classification')
    # display.cbs[2] = adjust_fhc_colorbar_for_pyart(display.cbs[2])
    # display.plot_range_rings([25, 50, 75, 100])
    # display.plot_grid_lines(ax=None, col='k', ls=':')
    # display.set_limits(xlim=[-100,100], ylim=[-100,100])
      
    #   RHOHV
    ax = fig.add_subplot(234)
    display.plot_ppi('cross_correlation_ratio', sweep = x , ax =  ax, #axislabels = ['km', 'km'],
                    vmin = 0.7, vmax = 1, mask_outside = False,
                    title='Copolar Correlation Ratio', colorbar_label='RhoHV', 
                    cmap = 'jet', gatefilter = gatefilter)
    display.plot_range_rings([25, 50, 75, 100])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-100,100], ylim=[-100,100])
    
    #   Specific Differential Phase
    ax = fig.add_subplot(235)
    display.plot_ppi('specific_differential_phase', sweep = x , ax =  ax, axislabels = ['East West distance from radar (km)', ''],
                    vmin = 0, vmax = 6, mask_outside = False,
                    title='Specific Diff. Phase', colorbar_label= 'Kdp (deg/km)',
                    gatefilter = gatefilter)
    display.plot_range_rings([25, 50, 75, 100])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-100,100], ylim=[-100,100])
      
    #   Differential Phase
    ax = fig.add_subplot(236)
    display.plot_ppi('differential_phase', sweep = x , ax =  ax, axislabels = ['East West distance from radar (km)', ''],
                    vmin = -180, vmax = 180, mask_outside = False,
                    title='Differential Phase', colorbar_label='Phidp (deg.)', )
                    #gatefilter = gatefilter)
    display.plot_range_rings([25, 50, 75, 100])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-100,100], ylim=[-100,100])
    plt.suptitle('\n ' + radar_name +' | Elevation: ' + angle + ' Deg.| ' + time_text, fontsize=16)
    
    return fig

  
def ppi_drops_1b(filename, sweep, radar_name):
    
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarDisplay(radar)
    radar.azimuth['data'] = radar.azimuth['data']
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    #gatefilter.exclude_invalid('corrected_reflectivity')
    #gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 20, 100)
    gatefilter.exclude_outside('normalized_coherent_power', 0.3, 1)
    gatefilter.exclude_outside('corrected_copol_correlation_ratio', 0.5, 1)
    
    
    #angle = str( round(sweep*30,2))
    angle = str( round(angles[x],2))
    fig = plt.figure(figsize = [19,10])
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S') #time_start.strftime('%Y-%m-%dT%H:%M:%S')
           
    #  Reflectivity                
    ax = fig.add_subplot(231)
    display.plot_ppi('corrected_reflectivity', sweep = x , ax = ax, axislabels = ['', 'North South distance from radar (km)'],
                    title = 'Corrected Reflectivity', colorbar_label='Z (dBZ)',
                    vmin = 0, vmax = 70, mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef,
                    gatefilter = gatefilter)
    display.plot_range_rings([10, 20, 30, 40, 50])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-50,50], ylim=[-50,50])
    
    #   Differential Reflectivity
    ax = fig.add_subplot(232)                
    display.plot_ppi('corrected_differential_reflectivity', sweep = x , ax = ax, axislabels = ['', ''],
                    vmin = -2, vmax = 8, mask_outside = False,
                    cmap = pyart.graph.cm.RefDiff, 
                    title='Diff. Reflectivity', colorbar_label='Zdr (dB)',
                    gatefilter = gatefilter)
    display.plot_range_rings([10, 20, 30, 40, 50])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-50,50], ylim=[-50,50])

    
    #   HyC
    ax = fig.add_subplot(233)
    hid_colors = ['White', 'Black','LightBlue', 'MediumBlue', 'DarkOrange', 'LightPink',
                  'Cyan', 'DarkGray', 'Lime', 'Yellow', 'Red', 'Fuchsia']
    cmaphid = colors.ListedColormap(hid_colors)
    display.plot_ppi('HydroClass', sweep = x ,ax = ax, 
                     axislabels = ['', ''],
                     cmap = cmaphid, colorbar_label='',
                     vmin = 5, vmax = 16, mask_outside = False, gatefilter = gatefilter,
                     title = 'Hydrometeor Classification')
    display.cbs[2] = adjust_fhc_colorbar_for_pyart(display.cbs[2])
    display.plot_range_rings([10, 20, 30, 40, 50])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-50,50], ylim=[-50,50])
      
    #   RHOHV
    ax = fig.add_subplot(234)
    display.plot_ppi('corrected_copol_correlation_ratio', sweep = x , ax =  ax, #axislabels = ['km', 'km'],
                    vmin = 0.7, vmax = 1, mask_outside = False,
                    title='Copolar Correlation Ratio', colorbar_label='RhoHV', 
                    cmap = 'jet', gatefilter = gatefilter)
    display.plot_range_rings([10, 20, 30, 40, 50])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-50,50], ylim=[-50,50])
    
    #   Specific Differential Phase
    ax = fig.add_subplot(235)
    display.plot_ppi('corrected_specific_differential_phase', sweep = x , ax =  ax, axislabels = ['East West distance from radar (km)', ''],
                    vmin = 0, vmax = 4, mask_outside = False,
                    title='Specific Diff. Phase', colorbar_label= 'Kdp (deg/km)',
                    gatefilter = gatefilter)
    display.plot_range_rings([10, 20, 30, 40, 50])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-50,50], ylim=[-50,50])
      
    #   Differential Phase
    ax = fig.add_subplot(236)
    display.plot_ppi('corrected_differential_phase', sweep = x , ax =  ax, axislabels = ['East West distance from radar (km)', ''],
                    vmin = -180, vmax = 180, mask_outside = False,
                    title='Differential Phase', colorbar_label='Phidp (deg.)', 
                    gatefilter = gatefilter)
    display.plot_range_rings([10, 20, 30, 40, 50])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-50,50], ylim=[-50,50])
    plt.suptitle('\n ' + radar_name +' | Elevation: ' + angle + ' Deg.| ' + time_text, fontsize=16)
    
    return fig



def rhi_drops_1b(filename, sweep):
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarMapDisplay(radar)
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('corrected_reflectivity')
    #gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 0, 100)
    #gatefilter.exclude_outside('corrected_copol_correlation_ratio', 0.7, 1)
    
    angle = str( round(angles[x] - 7,5)) # azimuth correction of 13 deg
    fig = plt.figure(figsize = [15,10])
    
    #  Reflectivity                
    ax = fig.add_subplot(321)
    display.plot_rhi('corrected_reflectivity', sweep = x , ax = ax, title = '',
                     axislabels = ['', 'Distance Above radar (km)'],
                    vmin = 0, vmax = 70, mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef, colorbar_label =  'Corrected Reflectivity (dBZ)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[30, 120])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   Differential Reflectivity
    ax = fig.add_subplot(322)                
    display.plot_rhi('corrected_differential_reflectivity', sweep = x , ax = ax,
                     axislabels = ['', ''],
                    vmin = -1, vmax = 6, mask_outside = False,
                    cmap = pyart.graph.cm.RefDiff, title='', colorbar_label='Diff. Reflectivity (dB)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[30, 120])
    display.plot_grid_lines(ax=None, col='k', ls=':')

    #   KDP
    ax = fig.add_subplot(323)
    display.plot_rhi('corrected_specific_differential_phase', sweep = x , ax =  ax, 
                    axislabels = ['', 'Distance Above radar (km)'],
                    vmin = 0, vmax = 6, mask_outside = False,
                    title='', colorbar_label= 'Specific Diff. Phase (deg/km)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[30, 120])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   Velocity
    ax = fig.add_subplot(324)
    display.plot_rhi('velocity', sweep = x , ax =  ax,
                    axislabels = ['', ' '], 
                    vmin = -15, vmax = 15, mask_outside = False, cmap = pyart.graph.cm.NWSVel,
                     title='', colorbar_label='Radial Velocity (m/s)',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[30, 120])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   RhoHV
    ax = fig.add_subplot(325)
    display.plot_rhi('corrected_copol_correlation_ratio', sweep = x , ax =  ax,
                     vmin = 0.8, vmax = 1, mask_outside = False, cmap = 'jet',
                     title='', colorbar_label='Copolar Correlcation',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[30, 120])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   HyC
    
    hid_colors = ['White', 'Black','LightBlue', 'MediumBlue', 'DarkOrange', 'LightPink',
                  'Cyan', 'DarkGray', 'Lime', 'Yellow', 'Red', 'Fuchsia']
    cmaphid = colors.ListedColormap(hid_colors)
    ax = fig.add_subplot(326)
    display.plot_rhi('HydroClass', sweep = x ,ax = ax, 
                     axislabels = ['Distance from radar (km)', ' '],
                     cmap = cmaphid, colorbar_label='',
                     vmin = 5, vmax = 16, mask_outside = False, gatefilter = gatefilter,
                     title = ' ')
    display.set_limits(ylim=[0, 21])
    display.cbs[5] = adjust_fhc_colorbar_for_pyart(display.cbs[5])
    display.set_limits(ylim=[0, 21])
    display.set_limits(xlim=[30, 120])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    angle = str( round(angles[x],2))    
    plt.suptitle('\n CSU-CHIVO| Azimuth: ' + angle + ' Deg.| ' + time_text, fontsize=16)
    print('check line 385 at ' + '/top/students/GRAD/ECE/idariash/home/anaconda3/lib/python3.6/site-packages/pyart/core/radar.py')
    return fig

def rhi_drops_1b_4plots(filename, sweep):
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarMapDisplay(radar)
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('corrected_reflectivity')
    #gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 10, 100)
    #gatefilter.exclude_outside('corrected_copol_correlation_ratio', 0.7, 1)
    
    fig = plt.figure(figsize = [12,6])
    
    #  Reflectivity                
    ax = fig.add_subplot(221)
    display.plot_rhi('corrected_reflectivity', sweep = x , ax = ax, title = '',
                     axislabels = ['', ''],
                    vmin = 0, vmax = 70, mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef, colorbar_label =  '',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 14])
    display.set_limits(xlim=[20, 50])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   Differential Reflectivity
    ax = fig.add_subplot(222)                
    display.plot_rhi('corrected_differential_reflectivity', sweep = x , ax = ax,
                    axislabels = ['', ''],
                    vmin = -3, vmax = 6, mask_outside = False,
                    cmap = pyart.graph.cm.RefDiff, title='', colorbar_label='',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 14])
    display.set_limits(xlim=[20, 50])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   KDP
    ax = fig.add_subplot(223)
    display.plot_rhi('corrected_specific_differential_phase', sweep = x , ax =  ax, 
                    axislabels = ['', ''],
                    vmin = 0, vmax = 6, mask_outside = False,
                    title='', colorbar_label= '',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 14])
    display.set_limits(xlim=[20, 50])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    
   
    #   rhohv
    
    ax = fig.add_subplot(224)
    display.plot_rhi('corrected_copol_correlation_ratio', sweep = x , ax =  ax,
                     axislabels = ['',''],
                     vmin = 0.5, vmax = 1, mask_outside = False, cmap = 'jet',
                     title='', colorbar_label='',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 14])
    display.set_limits(xlim=[20, 50])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    
    angle = str( round(angles[x],2))    
    #plt.suptitle('CSU-CHIVO| Azimuth: ' + angle + ' Deg.| ' + time_text , fontsize=12)
    print('check line 385 at ' + '/top/students/GRAD/ECE/idariash/home/anaconda3/lib/python3.6/site-packages/pyart/core/radar.py')
    return fig

def rhi_drops_1b_4plots_phidp(filename, sweep):
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarMapDisplay(radar)
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('corrected_reflectivity')
    #gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 10, 100)
    #gatefilter.exclude_outside('corrected_copol_correlation_ratio', 0.7, 1)
    
    fig = plt.figure(figsize = [12,6])
    
    #  Reflectivity                
    ax = fig.add_subplot(221)
    display.plot_rhi('corrected_reflectivity', sweep = x , ax = ax, title = '',
                     axislabels = ['', ''],
                    vmin = 0, vmax = 70, mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef, colorbar_label =  '',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 14])
    display.set_limits(xlim=[20, 50])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   Differential Reflectivity
    ax = fig.add_subplot(222)                
    display.plot_rhi('corrected_differential_reflectivity', sweep = x , ax = ax,
                    axislabels = ['', ''],
                    vmin = -3, vmax = 6, mask_outside = False,
                    cmap = pyart.graph.cm.RefDiff, title='', colorbar_label='',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 14])
    display.set_limits(xlim=[20, 50])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   PHIDP
    ax = fig.add_subplot(223)
    
    
    display.plot_rhi('corrected_differential_phase', sweep = x , ax =  ax,
                    axislabels = ['', ''],
                    vmin = 0, vmax = 60, mask_outside = False,
                    title='', colorbar_label='',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 14])
    display.set_limits(xlim=[20, 50])
    display.plot_grid_lines(ax=None, col='k', ls=':') 
   
    #   rhohv
    
    ax = fig.add_subplot(224)
    display.plot_rhi('corrected_copol_correlation_ratio', sweep = x , ax =  ax,
                     axislabels = ['',''],
                     vmin = 0.5, vmax = 1, mask_outside = False, cmap = 'jet',
                     title='', colorbar_label='',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 14])
    display.set_limits(xlim=[20, 50])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    
    angle = str( round(angles[x],2))    
    #plt.suptitle('CSU-CHIVO| Azimuth: ' + angle + ' Deg.| ' + time_text , fontsize=12)
    print('check line 385 at ' + '/top/students/GRAD/ECE/idariash/home/anaconda3/lib/python3.6/site-packages/pyart/core/radar.py')
    return fig


def rhi_drops_1b_4plots_vel(filename_1a, filename_1b, sweep):
    
    radar = pyart.io.read(filename_1b)
    radar_1a = pyart.io.read(filename_1a)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarMapDisplay(radar)
    display_1a = pyart.graph.RadarMapDisplay(radar_1a)
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('corrected_reflectivity')
    #gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 0, 100)
    #gatefilter.exclude_outside('corrected_copol_correlation_ratio', 0.7, 1)
    
    
    gatefilter_1a = pyart.filters.GateFilter(radar_1a)
    gatefilter_1a.exclude_transition()
    gatefilter_1a.exclude_invalid('reflectivity')
    #gatefilter.exclude_invalid('differential_phase')
    gatefilter_1a.exclude_outside('reflectivity', 0, 100)
    #gatefilter_1a.exclude_outside('normalized_coherent_power', 0.2, 1)
    #gatefilter_1a.exclude_outside('cross_correlation_ratio', 0.5, 1)
    
    fig = plt.figure(figsize = [12,6])
    
    #  Reflectivity                
    ax = fig.add_subplot(221)
    display.plot_rhi('corrected_reflectivity', sweep = x , ax = ax, title = '',
                     axislabels = ['', ''],
                    vmin = 0, vmax = 70, mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef, colorbar_label =  '',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 18])
    display.set_limits(xlim=[20, 60])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   Differential Reflectivity
    ax = fig.add_subplot(222)                
    display.plot_rhi('corrected_differential_reflectivity', sweep = x , ax = ax,
                    axislabels = ['', ''],
                    vmin = -3, vmax = 6, mask_outside = False,
                    cmap = pyart.graph.cm.RefDiff, title='', colorbar_label='',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 18])
    display.set_limits(xlim=[20, 60])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    #   VEL
    ax = fig.add_subplot(223)
    last_radar = pyart.io.read(filename_1a)
    
    nyq = radar_1a.instrument_parameters['nyquist_velocity']['data'][0]
    corr_vel = pyart.correct.dealias_region_based(
            last_radar, vel_field='velocity', keep_original=False, 
            gatefilter = gatefilter_1a, nyquist_vel=nyq, centered = True)
    radar_1a.add_field('corrected_velocity', corr_vel, replace_existing = True)
    
    radar_1a.fields['corrected_velocity']['data'] = (radar_1a.fields['corrected_velocity']['data'] + 1*nyq)
   
    # corr_vel_2nd_time = pyart.correct.dealias_fourdd(
    #         radar_1a, last_radar = last_radar, gatefilter = gatefilter_1a)
    # radar_1a.add_field('corrected_velocity', corr_vel_2nd_time, replace_existing = True)
    
    
    display_1a.plot_rhi('corrected_velocity', sweep = x , ax =  ax, 
                    axislabels = ['', ''],
                    #vmin = -3*nyq, vmax = 3*nyq, mask_outside = False,
                    vmin = -40, vmax = 40, mask_outside = False,
                    title='', colorbar_label= '',
                    gatefilter = gatefilter, cmap = pyart.graph.cm.NWSVel)
    display.set_limits(ylim=[0, 18])
    display.set_limits(xlim=[20, 60])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    
   
    #   rhohv
    
    ax = fig.add_subplot(224)
    display.plot_rhi('corrected_copol_correlation_ratio', sweep = x , ax =  ax,
                     axislabels = ['',''],
                     vmin = 0.5, vmax = 1, mask_outside = False, cmap = 'jet',
                     title='', colorbar_label='',
                    gatefilter = gatefilter)
    display.set_limits(ylim=[0, 18])
    display.set_limits(xlim=[20, 60])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    
    
    angle = str( round(angles[x],2))    
    #plt.suptitle('CSU-CHIVO| Azimuth: ' + angle + ' Deg.| ' + time_text , fontsize=12)
    print('check line 385 at ' + '/top/students/GRAD/ECE/idariash/home/anaconda3/lib/python3.6/site-packages/pyart/core/radar.py')
    return fig
    


#-------------------- For Test -----------------------

#filename = '/net/denali/storage/radar/RELAMPAGO/DROPS/2018/12/14/cfrad.20181214_020602.120_to_20181214_020923.124_col-radar_REL_RHI45_RHI.nc'
#filename_sigmet = '/net/denali/storage/radar/RELAMPAGO/RAW/2018/12/14/COL20181214_020602'
#fig = chivo_rhi_drops_4dd(filename, filename_sigmet, 11, 80, 1.7)
