#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 16:13:52 2019

@author: idariash
"""

import pyart

filename = '/net/denali/storage/radar/RELAMPAGO/quality_controlled_data/level_1a/2019/01/30/chivo.1a.20190130_083700.REL_RHI45.nc'
radar = pyart.io.read(filename)
display = pyart.graph.RadarMapDisplay(radar)
gatefilter = pyart.filters.GateFilter(radar)
x = 4
display.plot_rhi('reflectivity', sweep = x, #title = '',
                    axislabels = ['', 'Distance Above radar (km)'],
                    vmin = -15, vmax = 80, mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef, colorbar_label =  'Corrected Reflectivity (dBZ)',
                    gatefilter = gatefilter)
display.set_limits(ylim=[0, 21])
display.set_limits(xlim=[0, 80])
display.plot_grid_lines(ax=None, col='k', ls=':')