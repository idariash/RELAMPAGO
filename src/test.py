#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 17:57:52 2019

@author: idariash
"""

# test for echotops in 2D

import pyart
import numpy as np
import scipy.ndimage as spyi
import netCDF4
import gc


filename = '/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/level_1b.2/2018/12/14/chivo.1b.2.20181214_024046.REL_PNL360A.nc'
   

radar = pyart.io.read(filename)
echo_top_height = 0
angles = radar.fixed_angle['data'];
X, Y, z = radar.get_gate_x_y_z(len(angles) - 1, edges=False)


for sweep in range(len(angles) - 1, -1, -1): #start from the highest sweep
    reflectivity = radar.get_field(sweep, 'corrected_reflectivity')
    rhohv = radar.get_field(sweep, 'corrected_copol_correlation_ratio')
    ncp = radar.get_field(sweep, 'normalized_coherent_power')
    x, y, z = radar.get_gate_x_y_z(sweep, edges=False)
    # x /= 1000.0
    # y /= 1000.0
    # z /= 1000.0
    # calculate (R)ange
    R = np.sqrt(x ** 2 + y ** 2)
    
    x_diff = X - x
    y_diff = Y- y
    
    print(np.amax(X - x))
    print(np.amax(Y - y))
    X = x
    Y = y
    continue

    reflectivity = spyi.gaussian_filter(reflectivity, sigma = 2)
    rhohv = spyi.gaussian_filter(rhohv, sigma = 2)
    
    # Masking
    z = np.ma.masked_where( z < echo_top_height, z)
    z = np.ma.masked_where( R > 100, z)
    z = np.ma.masked_where( reflectivity < 18, z)
    z = np.ma.masked_where( rhohv < 0.9 , z)
    z = np.ma.masked_where( ncp < 0.3 , z)
    echo_top_sweep = np.amax(np.amax(z))
    if echo_top_sweep > echo_top_height:
        echo_top_height = echo_top_sweep
        ind_max_z = np.unravel_index(np.argmax(z, axis=None), z.shape)
        x_echo_top_height = x[ind_max_z]
        y_echo_top_height = y[ind_max_z]

