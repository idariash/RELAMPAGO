#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 12:05:00 2020

@author: idariash
"""

import pyart
import time

def grid_radar(radar, grid_shape=(40, 601, 601), xlim=(-150000, 150000), # xlim=(-70000, 60000), #grid_shape=(20, 401, 401), xlim=(-100000, 100000)
               ylim=(-150000, 150000), zlim=(1000, 20000),
               fields=['corrected_reflectivity', 'corrected_differential_reflectivity'], origin=None):
        
    if origin is None:
        origin = (radar.latitude['data'][0],
                  radar.longitude['data'][0])

    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('corrected_reflectivity')
    gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 0, 100)
    gatefilter.exclude_outside('normalized_coherent_power', 0.6, 1)
    #gatefilter.exclude_outside('corrected_copol_correlation_ratio', 0.3, 1)    

    grid = pyart.map.grid_from_radars(
        radar, grid_shape=grid_shape,
        grid_limits=(zlim, ylim, xlim),
        grid_origin=origin, fields=fields, gatefilters = gatefilter,
        gridding_algo='map_gates_to_grid', grid_origin_alt=0.0)

    
    return grid
