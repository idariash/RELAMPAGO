#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 14:35:01 2021

@author: idariash
"""

import pyart
import netCDF4
import matplotlib.pyplot as plt

def ppi_reflectivity(filename, sweep):
    
    radar = pyart.io.read(filename)
    x = sweep
    angles = radar.fixed_angle['data'];
    display = pyart.graph.RadarDisplay(radar)
    
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('DBZ')
    #gatefilter.exclude_invalid('differential_phase')
    gatefilter.exclude_outside('DBZ', 0, 100)
    gatefilter.exclude_outside('RHOHV', 0.90, 1)
    
    angle = str( round(angles[x],2))
    fig = plt.figure(figsize = [10,8])
    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M:%S')
           
    #  Reflectivity               
    display.plot_ppi('DBZ', sweep = x , # axislabels = ['', 'North South distance from radar (km)'],
                    title = 'Reflectivity', colorbar_label='Ref. (dBZ)',
                    vmin = 0, vmax = 70, mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef, gatefilter = gatefilter) #pyart.graph.cm.NWSRef
    #display.plot_range_rings([25, 50, 75])
    display.plot_grid_lines(ax=None, col='k', ls=':')
    display.set_limits(xlim=[-100,100], ylim=[-100,100])
    
    
    plt.suptitle('\n CSU-CHIVO| Elevation: ' + angle + ' Deg.| ' + time_text, fontsize=16)
    
    return fig