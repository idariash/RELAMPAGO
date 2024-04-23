#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 10:45:02 2021

@author: idariash
"""

import sys
sys.path.append('/net/denali/storage2/radar2/tmp/Ivan/CSU/RELAMPAGO/plot')

import pyart
import glob
import netCDF4
import numpy as np

import chivo_grid_utl as utl
import scipy.ndimage as spyi
import scipy.io


DataPath = '/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/test/data_1b_whitney'

files = [f for f in glob.glob(DataPath +"/*P*.nc", recursive=True)]
files.sort()

for file in files:
    radar = pyart.io.read(file)
    grid = utl.grid_radar(radar)
    y, z, x = np.meshgrid(0.001*grid.y['data'],  0.001*grid.z['data'],  0.001*grid.x['data'])
    reflectivity = grid.fields['corrected_reflectivity']['data']
    RainRate = grid.fields['RainRate']['data']
    HydroClass = grid.fields['HydroClass']['data'][2]
    #reflectivity = np.swapaxes(reflectivity, 0, 2)
    reflectivity = spyi.gaussian_filter(reflectivity, sigma = 4)
    
    z = np.ma.masked_where( reflectivity < 18, z)
    echo_top = np.amax(z, axis=0)
    echo_top = np.ma.masked_where( echo_top < 4, echo_top)
    
    
    #reflectivity = np.swapaxes(reflectivity, 0, 2)
    max_ref = np.amax(reflectivity, axis=0)
    RainRate = np.amax(RainRate, axis=0)
    x, y = np.meshgrid(0.001*grid.x['data'],  0.001*grid.y['data'])
    
    time_start = netCDF4.num2date(grid.time['data'][0], grid.time['units'])
    time_text = time_start.strftime('%Y%m%dT%H%M')
    
    scipy.io.savemat('./mat_files/' + 'chivo-' + time_text + '.mat',
                     dict(x=x, y=y, time = time_text, echo_top = echo_top,
                         RainRate = RainRate, max_ref = max_ref, HydroClass = HydroClass))
    