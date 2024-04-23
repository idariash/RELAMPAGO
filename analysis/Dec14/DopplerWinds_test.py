#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 18:41:54 2021

@author: idariash
"""

import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
import multidop
import pyart
import pydda
import tempfile
import os
import glob
import time

filename_chivo = '/net/denali/storage2/radar2/people/idariash/home/CSU/RELAMPAGO/analysis/DualRadar/rsc/chivo.1a.20181214_020047.REL_PFAR360.nc'
filename_csapr = '/net/denali/storage2/radar2/people/idariash/home/CSU/RELAMPAGO/analysis/DualRadar/rsc/corcsapr2cfrppiM1.a1.20181214.020004.nc'

r2 = pyart.io.read(filename_chivo)
r1 = pyart.io.read(filename_csapr)

#%% Dealias velocity

# CHIVO r2
gatefilter = pyart.filters.GateFilter(r2)
gatefilter.exclude_transition()
gatefilter.exclude_invalid('reflectivity')
#gatefilter.exclude_invalid('differential_phase')
gatefilter.exclude_outside('reflectivity', 8, 100)
#gatefilter.exclude_outside('cross_correlation_ratio', 0.7, 1)

nyq = r2.instrument_parameters['nyquist_velocity']['data'][0]

corr_vel = pyart.correct.dealias_region_based(
            r2, vel_field='velocity', keep_original=False, 
            gatefilter = gatefilter, nyquist_vel=nyq, centered = True)

r2.add_field('corrected_velocity', corr_vel, replace_existing = True)

# CSAPR r1
gatefilter = pyart.filters.GateFilter(r1)
gatefilter.exclude_transition()
gatefilter.exclude_invalid('mean_doppler_velocity')
gatefilter.exclude_invalid('reflectivity')
gatefilter.exclude_outside('reflectivity', 1.5, 100)
gatefilter.exclude_outside('normalized_coherent_power', 0.2, 1)
#gatefilter.exclude_outside('copol_correlation_coeff', 0.7, 1)

nyq = r1.instrument_parameters['nyquist_velocity']['data'][0]

corr_vel = pyart.correct.dealias_region_based(
            r1, vel_field='mean_doppler_velocity', keep_original=False, 
            gatefilter = gatefilter, nyquist_vel=nyq, centered = True)

r1.add_field('corrected_velocity', corr_vel, replace_existing = True)

cp = deepcopy(r1.fields['reflectivity']['data'])
r1.add_field_like('reflectivity', 'DT', cp, replace_existing=True)
cp = deepcopy(r1.fields['corrected_velocity']['data'])
r1.add_field_like('corrected_velocity', 'VT', cp, replace_existing=True)

cp = deepcopy(r2.fields['reflectivity']['data'])
r2.add_field_like('reflectivity', 'DT', cp, replace_existing=True)
cp = deepcopy(r2.fields['corrected_velocity']['data'])
r2.add_field_like('corrected_velocity', 'VT', cp, replace_existing=True)

#%%
# The analysis engine currently expects the "missing_value" attribute
r1.fields['DT']['missing_value'] = 1.0 * r1.fields['DT']['_FillValue']
r2.fields['DT']['missing_value'] = 1.0 * r2.fields['DT']['_FillValue']
r1.fields['VT']['missing_value'] = 1.0 * r1.fields['VT']['_FillValue']
r2.fields['VT']['missing_value'] = 1.0 * r2.fields['VT']['_FillValue']

#%% Now grid the volumes and add azimuths and elevations afterward
def grid_radar(radar, grid_shape=(20, 301, 301), xlim=(-150000, 150000),
               ylim=(-150000, 150000), zlim=(1000, 20000),
               fields=['DT', 'VT'], origin=None):
    bt = time.time()
    radar_list = [radar]
    if origin is None:
        origin = (radar.latitude['data'][0],
                  radar.longitude['data'][0])
    grid = pyart.map.grid_from_radars(
        radar_list, grid_shape=grid_shape,
        grid_limits=(zlim, ylim, xlim),
        grid_origin=origin, fields=fields,
        gridding_algo='map_gates_to_grid', grid_origin_alt=0.0)
    print(time.time()-bt, 'seconds to grid radar')
    return grid

#%%
# Fix for Py-ART problem with r2's data structure
# This can be skipped if you don't have TypeErrors while gridding
r2.longitude['data'] = np.array([r2.longitude['data'][0]])
r2.altitude['data'] = np.array([r2.altitude['data'][0]])
r2.latitude['data'] = np.array([r2.latitude['data'][0]])

#%%
g1 = grid_radar(r1, origin=(r2.latitude['data'][0], r2.longitude['data'][0]),
                xlim=(-100000, 50000), ylim=(-100000, 50000), grid_shape=(20, 151, 151))
g2 = grid_radar(r2, origin=(r2.latitude['data'][0], r2.longitude['data'][0]),
                xlim=(-100000, 50000), ylim=(-100000, 50000), grid_shape=(20, 151, 151))

# Set initialization and do retrieval
u_init, v_init, w_init = pydda.initialization.make_constant_wind_field(g1, vel_field='VT')
new_grids = pydda.retrieval.get_dd_wind_field([g1, g2],
                                              u_init, v_init, w_init,
                                              vel_name='VT', refl_field='DT',
                                              mask_outside_opt=True)
#%%
# Make a neat plot
fig = plt.figure(figsize=(10,7))
ax = pydda.vis.plot_horiz_xsection_quiver_map(new_grids, background_field='DT', level=3,
                                              show_lobes=False, bg_grid_no=3, vmin=0, vmax=60,
                                              quiverkey_len=40.0,
                                              quiver_spacing_x_km=2.0, quiver_spacing_y_km=2.0,
                                              quiverkey_loc='top', colorbar_contour_flag=True,
                                              cmap='pyart_HomeyerRainbow')