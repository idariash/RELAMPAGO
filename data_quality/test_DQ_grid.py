#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 15:45:55 2019

@author: idariash
"""

import pyart
import matplotlib.pyplot as plt

chivo = '/net/denali/storage2/radar2/tmp/Ivan/CSU/CHIVO/res/data/test_data_quality/chivo.1a.20190125_193049.REL_PNL135A.nc'
csapr = '/net/denali/storage2/radar2/tmp/Ivan/CSU/CHIVO/res/data/test_data_quality/corcsapr2cfrppiM1.a1.20190125.193003.nc'

CHIVO = pyart.io.read(chivo)
CSAPR =  pyart.io.read(csapr)

# filter out gates with differential_reflectivity > 100 from both radars
gatefilter_CHIVO = pyart.filters.GateFilter(CHIVO)
gatefilter_CHIVO.exclude_transition()
gatefilter_CHIVO.exclude_outside('reflectivity', 10, 100)
gatefilter_CHIVO.exclude_outside('normalized_coherent_power', 0.3, 1)

gatefilter_CSAPR = pyart.filters.GateFilter(CSAPR)
gatefilter_CSAPR.exclude_transition()
gatefilter_CSAPR.exclude_outside('reflectivity', 10, 100)
gatefilter_CSAPR.exclude_outside('normalized_coherent_power', 0.3, 1)
nyq_chivo = CSAPR.instrument_parameters['nyquist_velocity']['data'][0]
nyq_chivo = CHIVO.instrument_parameters['nyquist_velocity']['data'][0]

csapr_fields = CSAPR.fields


lat_mid = (CHIVO.latitude['data'][0] + CSAPR.latitude['data'][0])/2

lon_mid = (CHIVO.longitude['data'][0] + CSAPR.longitude['data'][0])/2

alt_mid = (CHIVO.altitude['data'][0] + CSAPR.altitude['data'][0])/2

#-----------------CSAPR Grid ----------------------------

grids = pyart.map.grid_from_radars(
        (CSAPR,), gatefilters=(gatefilter_CSAPR,), grid_shape=(41, 251, 251),
        grid_limits=((alt_mid, 20000.0),(-150000, 150000), (-150000, 150000)),
        fields=['mean_doppler_velocity'], gridding_algo= 'map_gates_to_grid', # 'map_to_grid', # 
        grid_origin = (lat_mid, lon_mid), grid_origin_alt = alt_mid,
        weighting_function='BARNES')

display = pyart.graph.GridMapDisplay(grids)
#fig = plt.figure(figsize=[15, 7])
fig = plt.figure(figsize=[15, 15])

# panel sizes
map_panel_axes = [0.05, 0.05, .4, .80]
x_cut_panel_axes = [0.55, 0.10, .4, .25]
y_cut_panel_axes = [0.55, 0.50, .4, .25]

# parameters
level = 10
level_1 = 20 
vmin = 0
vmax = 6
lat = -31.6342 #36.5
lon = -64.1586 #-98.0

# panel 1, basemap, radar differential_reflectivity and NARR overlay
#ax1 = fig.add_axes(map_panel_axes)
#display.plot_basemap(lon_lines = np.arange(-104, -93, 1) )
display.plot_grid('mean_doppler_velocity', level=level,  
                  vmin= -nyq_chivo, vmax= nyq_chivo, 
                  cmap = pyart.graph.cm.NWSVel)
display.plot_crosshairs(lon=lon, lat=lat)
display.plot_crosshairs(CSAPR.longitude['data'][0],CSAPR.latitude['data'][0])

display = pyart.graph.GridMapDisplay(grids)
#fig = plt.figure(figsize=[15, 7])
fig = plt.figure(figsize=[15, 15])

# panel sizes
map_panel_axes = [0.05, 0.05, .4, .80]
x_cut_panel_axes = [0.55, 0.10, .4, .25]
y_cut_panel_axes = [0.55, 0.50, .4, .25]

# parameters
level = 5
vmin = 0
vmax = 6
lat = -31.6342 #36.5
lon = -64.1586 #-98.0

# panel 1, basemap, radar differential_reflectivity and NARR overlay
#ax1 = fig.add_axes(map_panel_axes)
#display.plot_basemap(lon_lines = np.arange(-104, -93, 1) )
display.plot_grid('mean_doppler_velocity', level=level, 
                  vmin= -nyq_chivo, vmax= nyq_chivo, 
                  cmap = pyart.graph.cm.NWSVel)
display.plot_crosshairs(lon=lon, lat=lat)
display.plot_crosshairs(CSAPR.longitude['data'][0],CSAPR.latitude['data'][0])

#----------------CSAPR Corrected Velocity
grids = pyart.map.grid_from_radars(
        (CSAPR,), gatefilters=(gatefilter_CSAPR,), grid_shape=(41, 251, 251),
        grid_limits=((alt_mid, 20000.0),(-150000, 150000), (-150000, 150000)),
        fields=['corrected_velocity'], gridding_algo= 'map_gates_to_grid', # 'map_to_grid', # 
        grid_origin = (lat_mid, lon_mid), grid_origin_alt = alt_mid,
        weighting_function='BARNES')

display = pyart.graph.GridMapDisplay(grids)
#fig = plt.figure(figsize=[15, 7])
fig = plt.figure(figsize=[15, 15])

# panel sizes
map_panel_axes = [0.05, 0.05, .4, .80]
x_cut_panel_axes = [0.55, 0.10, .4, .25]
y_cut_panel_axes = [0.55, 0.50, .4, .25]

# parameters
level = 10
level_1 = 20 
vmin = 0
vmax = 6
lat = -31.6342 #36.5
lon = -64.1586 #-98.0

# panel 1, basemap, radar differential_reflectivity and NARR overlay
#ax1 = fig.add_axes(map_panel_axes)
#display.plot_basemap(lon_lines = np.arange(-104, -93, 1) )
display.plot_grid('corrected_velocity', level=level,  
                  vmin= -nyq_chivo, vmax= nyq_chivo, 
                  cmap = pyart.graph.cm.NWSVel)
display.plot_crosshairs(lon=lon, lat=lat)
display.plot_crosshairs(CSAPR.longitude['data'][0],CSAPR.latitude['data'][0])

display = pyart.graph.GridMapDisplay(grids)
#fig = plt.figure(figsize=[15, 7])
fig = plt.figure(figsize=[15, 15])

# panel sizes
map_panel_axes = [0.05, 0.05, .4, .80]
x_cut_panel_axes = [0.55, 0.10, .4, .25]
y_cut_panel_axes = [0.55, 0.50, .4, .25]

# parameters
level = 5
vmin = 0
vmax = 6
lat = -31.6342 #36.5
lon = -64.1586 #-98.0

# panel 1, basemap, radar differential_reflectivity and NARR overlay
#ax1 = fig.add_axes(map_panel_axes)
#display.plot_basemap(lon_lines = np.arange(-104, -93, 1) )
display.plot_grid('corrected_velocity', level=level, 
                  vmin= -nyq_chivo, vmax= nyq_chivo, 
                  cmap = pyart.graph.cm.NWSVel)
display.plot_crosshairs(lon=lon, lat=lat)
display.plot_crosshairs(CSAPR.longitude['data'][0],CSAPR.latitude['data'][0])

#-----------------CHIVO Grid ----------------------------

grids = pyart.map.grid_from_radars(
        (CHIVO,), gatefilters=(gatefilter_CHIVO,), grid_shape=(41, 251, 251),
        grid_limits=((alt_mid, 20000.0),(-150000, 150000), (-150000, 150000)),
        fields=['velocity'], gridding_algo= 'map_gates_to_grid', # 'map_to_grid', # 
        grid_origin = (lat_mid, lon_mid), grid_origin_alt = alt_mid,
        weighting_function='BARNES')

display = pyart.graph.GridMapDisplay(grids)
#fig = plt.figure(figsize=[15, 7])
fig = plt.figure(figsize=[15, 15])

# panel sizes
map_panel_axes = [0.05, 0.05, .4, .80]
x_cut_panel_axes = [0.55, 0.10, .4, .25]
y_cut_panel_axes = [0.55, 0.50, .4, .25]

# parameters
level = 10
level_1 = 20 
vmin = 0
vmax = 6
lat = -31.6342 #36.5
lon = -64.1586 #-98.0

# panel 1, basemap, radar differential_reflectivity and NARR overlay
#ax1 = fig.add_axes(map_panel_axes)
#display.plot_basemap(lon_lines = np.arange(-104, -93, 1) )
display.plot_grid('velocity', level=level,  
                  vmin= -nyq_chivo, vmax= nyq_chivo, 
                  cmap = pyart.graph.cm.NWSVel)
display.plot_crosshairs(lon=lon, lat=lat)
display.plot_crosshairs(CSAPR.longitude['data'][0],CSAPR.latitude['data'][0])

display = pyart.graph.GridMapDisplay(grids)
#fig = plt.figure(figsize=[15, 7])
fig = plt.figure(figsize=[15, 15])

# panel sizes
map_panel_axes = [0.05, 0.05, .4, .80]
x_cut_panel_axes = [0.55, 0.10, .4, .25]
y_cut_panel_axes = [0.55, 0.50, .4, .25]

# parameters
level = 5
vmin = 0
vmax = 6
lat = -31.6342 #36.5
lon = -64.1586 #-98.0

# panel 1, basemap, radar reflectivity and NARR overlay
#ax1 = fig.add_axes(map_panel_axes)
#display.plot_basemap(lon_lines = np.arange(-104, -93, 1) )
display.plot_grid('velocity', level=level, 
                  vmin= -nyq_chivo, vmax= nyq_chivo, 
                  cmap = pyart.graph.cm.NWSVel)
display.plot_crosshairs(lon=lon, lat=lat)
display.plot_crosshairs(CSAPR.longitude['data'][0],CSAPR.latitude['data'][0])

#------------CHIVO Corrected Velocity ----------------
grids = pyart.map.grid_from_radars(
        (CHIVO,), gatefilters=(gatefilter_CHIVO,), grid_shape=(41, 251, 251),
        grid_limits=((alt_mid, 20000.0),(-150000, 150000), (-150000, 150000)),
        fields=['corrected_velocity'], gridding_algo= 'map_gates_to_grid', # 'map_to_grid', # 
        grid_origin = (lat_mid, lon_mid), grid_origin_alt = alt_mid,
        weighting_function='BARNES')

display = pyart.graph.GridMapDisplay(grids)
#fig = plt.figure(figsize=[15, 7])
fig = plt.figure(figsize=[15, 15])

# panel sizes
map_panel_axes = [0.05, 0.05, .4, .80]
x_cut_panel_axes = [0.55, 0.10, .4, .25]
y_cut_panel_axes = [0.55, 0.50, .4, .25]

# parameters
level = 10
level_1 = 20 
vmin = 0
vmax = 6
lat = -31.6342 #36.5
lon = -64.1586 #-98.0

# panel 1, basemap, radar differential_reflectivity and NARR overlay
#ax1 = fig.add_axes(map_panel_axes)
#display.plot_basemap(lon_lines = np.arange(-104, -93, 1) )
display.plot_grid('corrected_velocity', level=level,  
                  vmin= -nyq_chivo, vmax= nyq_chivo, 
                  cmap = pyart.graph.cm.NWSVel)
display.plot_crosshairs(lon=lon, lat=lat)
display.plot_crosshairs(CSAPR.longitude['data'][0],CSAPR.latitude['data'][0])

display = pyart.graph.GridMapDisplay(grids)
#fig = plt.figure(figsize=[15, 7])
fig = plt.figure(figsize=[15, 15])

# panel sizes
map_panel_axes = [0.05, 0.05, .4, .80]
x_cut_panel_axes = [0.55, 0.10, .4, .25]
y_cut_panel_axes = [0.55, 0.50, .4, .25]

# parameters
level = 5
vmin = 0
vmax = 6
lat = -31.6342 #36.5
lon = -64.1586 #-98.0

# panel 1, basemap, radar reflectivity and NARR overlay
#ax1 = fig.add_axes(map_panel_axes)
#display.plot_basemap(lon_lines = np.arange(-104, -93, 1) )
display.plot_grid('corrected_velocity', level=level, 
                  vmin= -nyq_chivo, vmax= nyq_chivo, 
                  cmap = pyart.graph.cm.NWSVel)
display.plot_crosshairs(lon=lon, lat=lat)
display.plot_crosshairs(CSAPR.longitude['data'][0],CSAPR.latitude['data'][0])