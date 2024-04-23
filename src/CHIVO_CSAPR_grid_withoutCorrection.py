#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  7 11:30:21 2019

@author: idariash
"""

import pyart
from matplotlib import pyplot as plt
import numpy as np
from time import time
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pyart



chivo = '/net/denali/storage/radar/RELAMPAGO/quality_controlled_data/level_1a/2019/01/25/chivo.1a.20190125_193049.REL_PNL135A.nc'; 
#'/net/denali/storage/radar/RELAMPAGO/quality_controlled_data/level_1a/2019/01/25/chivo.1a.20190125_210047.REL_PNL135A.nc' 
#'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Attenuation/NetCDF/20190125_2100/cfrad.20190125_210047.138_to_20190125_210408.475_col-radar_REL_PNL135A_PPI.nc' 
#'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Attenuation/NetCDF/20190125_2100/DROPS/cfrad.20190125_210047.138_to_20190125_210408.475_col-radar_REL_PNL135A_PPI.nc' 
#'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/CHIVO_CSAPR/data/DROPS/NetCDF_drops/cfrad.20190126_053009.397_to_20190126_053642.486_col-radar_PPINEARLT360_SUR.nc'
#'/net/denali/storage/radar/RELAMPAGO/quality_controlled_data/level_1a/2019/01/26/chivo.1a.20190126_053023.PPINEARLT360.nc'
csapr = '/net/denali/storage2/radar2/tmp/Ivan/DOE/RELAMPAGO/CSAPR_PPI/corcsapr2cfrppiM1.a1.20190125.193003.nc' 
'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Attenuation/NetCDF/20190125_2100/corcsapr2cfrppiM1.a1.20190125.210003.nc'
#'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Attenuation/NetCDF/20190125_2100/DROPS/corcsapr2cfrppiM1.a1.20190125.210003.nc'
#'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/CHIVO_CSAPR/data/DROPS/NetCDF_drops/corcsapr2cfrppiM1.a1.20190126.053003.nc'
#'/net/denali/storage2/radar2/tmp/Ivan/DOE/RELAMPAGO/20190126/CSAPR/NetCDF/211907/corcsapr2cfrppiM1.a1.20190126.053003.nc'

CHIVO = pyart.io.read(chivo)
CSAPR =  pyart.io.read(csapr)

# Zdr bias correction
CSAPR.fields['differential_reflectivity']['data'] = CSAPR.fields['differential_reflectivity']['data'] - (-3.54)
#CHIVO.fields['differential_reflectivity']['data'] = CHIVO.fields['differential_reflectivity']['data'] - 0.75
#CHIVO.azimuth['data'] = CHIVO.azimuth['data'] - 13

gatefilter_CHIVO = pyart.filters.GateFilter(CHIVO)
gatefilter_CHIVO.exclude_transition()
gatefilter_CHIVO.exclude_outside('reflectivity', 10, 100)
gatefilter_CHIVO.exclude_outside('cross_correlation_ratio', 0.7, 1)

gatefilter_CSAPR = pyart.filters.GateFilter(CSAPR)
gatefilter_CSAPR.exclude_transition()
gatefilter_CSAPR.exclude_outside('reflectivity', 10, 100)
gatefilter_CSAPR.exclude_outside('copol_correlation_coeff', 0.7, 1)

a = CSAPR.fields


#Grid origin
lat_mid = (CHIVO.latitude['data'][0] + CSAPR.latitude['data'][0])/2
lon_mid = (CHIVO.longitude['data'][0] + CSAPR.longitude['data'][0])/2
alt_mid = (CHIVO.altitude['data'][0] + CSAPR.altitude['data'][0])/2

#Creating the grids, resolution 100m
CSAPR_grid = pyart.map.grid_from_radars(
        (CSAPR,), gatefilters=(gatefilter_CSAPR,), grid_shape=(21, 351, 351),
        grid_limits=((1200, 11200.0),(-35000, 35000), (-35000, 35000)),
        fields=['reflectivity', 'differential_reflectivity', 'copol_correlation_coeff'], 
        gridding_algo= 'map_gates_to_grid', # 'map_to_grid', # 
        grid_origin = (lat_mid, lon_mid), grid_origin_alt = alt_mid,
        weighting_function='BARNES')

CHIVO_grid = pyart.map.grid_from_radars(
        (CHIVO,), gatefilters=(gatefilter_CHIVO,), grid_shape=(21, 351, 351),
        grid_limits=((1200, 11200.0),(-35000, 35000), (-35000, 35000)),
        fields=['reflectivity', 'differential_reflectivity', 'cross_correlation_ratio'],
        gridding_algo= 'map_gates_to_grid', # 'map_to_grid', # 
        grid_origin = (lat_mid, lon_mid), grid_origin_alt = alt_mid,
        weighting_function='BARNES')



#%%Display Grid
display = pyart.graph.GridMapDisplay(CSAPR_grid)
#fig = plt.figure(figsize=[15, 7])
fig = plt.figure(figsize=[8, 8])

# panel sizes
map_panel_axes = [1, 0.5, 0.5, 0.5]
#map_panel_axes = [0.05, 0.05, .4, .80]
x_cut_panel_axes = [0.55, 0.10, .4, .25]
y_cut_panel_axes = [0.55, 0.50, .4, .25]

# parameters
level = 3
vmin = 0
vmax = 3

#ax1 = fig.add_axes(map_panel_axes)
display.plot_basemap(lon_lines = np.arange(-100, 0, 0.1))
display.plot_basemap(lat_lines = np.arange(-100, 0, 0.1)) 

display.plot_grid('differential_reflectivity', level=level,  
                  vmin= -2, vmax= 6, #cmap = 'viridis') 
                  cmap = pyart.graph.cm.RefDiff)

display.plot_crosshairs(CSAPR.longitude['data'][0],CSAPR.latitude['data'][0])
display.plot_crosshairs(CHIVO.longitude['data'][0],CHIVO.latitude['data'][0])

a = CHIVO.longitude['data'][0]
b = CSAPR.longitude['data'][0]

display = pyart.graph.GridMapDisplay(CHIVO_grid)
#fig = plt.figure(figsize=[15, 7])
fig = plt.figure(figsize=[8, 8])

# panel sizes
map_panel_axes = [0.05, 0.05, .4, .80]
x_cut_panel_axes = [0.55, 0.10, .4, .25]
y_cut_panel_axes = [0.55, 0.50, .4, .25]

# parameters
vmin = 0
vmax = 6

display.plot_basemap(lon_lines = np.arange(-100, 0, 0.1))
display.plot_basemap(lat_lines = np.arange(-100, 0, 0.1))

display.plot_grid('differential_reflectivity', level=level,  
                  vmin= -2, vmax= 6, #cmap = 'viridis') 
                  cmap = pyart.graph.cm.RefDiff)
display.plot_crosshairs(CSAPR.longitude['data'][0],CSAPR.latitude['data'][0])
display.plot_crosshairs(CHIVO.longitude['data'][0],CHIVO.latitude['data'][0])

#%%Export grids

#csapr_gridFilename = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/CHIVO_CSAPR/data/csapr_grid_20190126-0530.nc'
#chivo_gridFilename = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/CHIVO_CSAPR/data/chivo_grid_20190126-0530.nc'
#
#pyart.io.write_grid(csapr_gridFilename, CSAPR_grid, format='NETCDF4') 
#                    #write_proj_coord_sys=True)
#
#pyart.io.write_grid(chivo_gridFilename, CHIVO_grid, format='NETCDF4') 
#                    #write_proj_coord_sys=True)




