#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 15:44:57 2020

@author: idariash
"""

import pyart
import glob
import netCDF4
import numpy as np
from operator import and_

def grid_radar(radar, grid_shape=(40, 401, 401), xlim=(-100000, 100000), 
               ylim=(-100000, 100000), zlim=(1000, 20000),
               fields=['corrected_reflectivity', 
                       'corrected_differential_reflectivity', 'RainRate'], origin=None):
        
    if origin is None:
        origin = (radar.latitude['data'][0],
                  radar.longitude['data'][0])
        
    grid = pyart.map.grid_from_radars(
        radar, grid_shape=grid_shape,
        grid_limits=(zlim, ylim, xlim),
        grid_origin=origin, fields=fields,
        gridding_algo='map_gates_to_grid', grid_origin_alt=0.0)

    
    return grid

DataPath = '/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/test/data_1b_whitney'
'/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/WRF_Radar/rsc/in/'
files = [f for f in glob.glob(DataPath +"/*P*.nc", recursive=True)]
files.sort()
radar = pyart.io.read(files[0])
grid = grid_radar(radar)
accumulated_total_rain  = (grid.fields['RainRate']['data'].copy())*0 # every 10 min per file
accumulated_total_rain = np.ma.masked_where(accumulated_total_rain < 0, accumulated_total_rain)
accumulated_total_rain.filled(0)
accumulated_total_rain.mask = False

for file in files:
    radar = pyart.io.read(file)
    grid = grid_radar(radar)
    total_rain_file = (grid.fields['RainRate']['data'].copy())/6 
    total_rain_file = np.ma.masked_where(total_rain_file < 0, total_rain_file)
    total_rain_file.filled(0)
    total_rain_file.mask = False
    accumulated_total_rain  = accumulated_total_rain + total_rain_file
    accumulated_total_rain_dict = {'data': accumulated_total_rain, 'units': 'mm',
                                   'missing_value': grid.fields['RainRate']['missing_value'],
                                   }
    grid.add_field('accumulated_total_rain', accumulated_total_rain_dict, True)
    time_start = netCDF4.num2date(grid.time['data'][0], grid.time['units'])
    time_text = time_start.strftime('%Y%m%dT%H%M')
    grid_fname = ('../rsc/out/whitney_1b/chivo_grid_'+ time_text + '.nc')
    grid.write(grid_fname)
    