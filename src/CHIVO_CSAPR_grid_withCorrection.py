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



chivo = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Attenuation/NetCDF/20190125_2100/DROPS/cfrad.20190125_210047.138_to_20190125_210408.475_col-radar_REL_PNL135A_PPI.nc'
#'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Attenuation/NetCDF/20190125_1930/DROPS/cfrad.20190125_193049.064_to_20190125_193410.629_col-radar_REL_PNL135A_PPI.nc' 
#'/net/denali/storage/radar/RELAMPAGO/quality_controlled_data/level_1a/2019/01/25/chivo.1a.20190125_193049.REL_PNL135A.nc'; 
#'/net/denali/storage/radar/RELAMPAGO/quality_controlled_data/level_1a/2019/01/25/chivo.1a.20190125_210047.REL_PNL135A.nc' 
#'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Attenuation/NetCDF/20190125_2100/cfrad.20190125_210047.138_to_20190125_210408.475_col-radar_REL_PNL135A_PPI.nc' 
#'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Attenuation/NetCDF/20190125_2100/DROPS/cfrad.20190125_210047.138_to_20190125_210408.475_col-radar_REL_PNL135A_PPI.nc' 
#'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/CHIVO_CSAPR/data/DROPS/NetCDF_drops/cfrad.20190126_053009.397_to_20190126_053642.486_col-radar_PPINEARLT360_SUR.nc'
#'/net/denali/storage/radar/RELAMPAGO/quality_controlled_data/level_1a/2019/01/26/chivo.1a.20190126_053023.PPINEARLT360.nc'
csapr = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Attenuation/NetCDF/20190125_2100/DROPS/corcsapr2cfrppiM1.a1.20190125.210003.nc'
#'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Attenuation/NetCDF/20190125_1930/DROPS/corcsapr2cfrppiM1.a1.20190125.193003.nc'
#'/net/denali/storage2/radar2/tmp/Ivan/DOE/RELAMPAGO/CSAPR_PPI/corcsapr2cfrppiM1.a1.20190125.193003.nc' 
#'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Attenuation/NetCDF/20190125_2100/corcsapr2cfrppiM1.a1.20190125.210003.nc'
#'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Attenuation/NetCDF/20190125_2100/DROPS/corcsapr2cfrppiM1.a1.20190125.210003.nc'
#'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/CHIVO_CSAPR/data/DROPS/NetCDF_drops/corcsapr2cfrppiM1.a1.20190126.053003.nc'
#'/net/denali/storage2/radar2/tmp/Ivan/DOE/RELAMPAGO/20190126/CSAPR/NetCDF/211907/corcsapr2cfrppiM1.a1.20190126.053003.nc'

CHIVO = pyart.io.read(chivo)
CSAPR =  pyart.io.read(csapr)

# Zdr bias correction
#CSAPR.fields['corrected_differential_reflectivity']['data'] = CSAPR.fields['corrected_differential_reflectivity']['data'] - (-3.54)
#CHIVO.fields['corrected_differential_reflectivity']['data'] = CHIVO.fields['corrected_differential_reflectivity']['data'] - 0.75
CHIVO.azimuth['data'] = CHIVO.azimuth['data'] - 13

gatefilter_CHIVO = pyart.filters.GateFilter(CHIVO)
gatefilter_CHIVO.exclude_transition()
gatefilter_CHIVO.exclude_outside('corrected_reflectivity', 13, 100)
gatefilter_CHIVO.exclude_outside('corrected_cross_correlation_ratio', 0.7, 1)
gatefilter_CHIVO.exclude_outside('corrected_differential_reflectivity', -2, 6)

gatefilter_CSAPR = pyart.filters.GateFilter(CSAPR)
gatefilter_CSAPR.exclude_transition()
gatefilter_CSAPR.exclude_outside('corrected_reflectivity', 13, 100)
gatefilter_CSAPR.exclude_outside('corrected_cross_correlation_ratio', 0.7, 1)


a = CSAPR.fields


#Grid origin
lat_mid = (CHIVO.latitude['data'][0] + CSAPR.latitude['data'][0])/2
lon_mid = (CHIVO.longitude['data'][0] + CSAPR.longitude['data'][0])/2
alt_mid = (CHIVO.altitude['data'][0] + CSAPR.altitude['data'][0])/2

#Creating the grids, resolution 100m
CSAPR_grid = pyart.map.grid_from_radars(
        (CSAPR,), gatefilters=(gatefilter_CSAPR,), grid_shape=(21, 351, 351),
        grid_limits=((1200, 11200.0),(-35000, 35000), (-35000, 35000)),
        fields=['corrected_reflectivity', 'corrected_differential_reflectivity', 'corrected_cross_correlation_ratio', 'corrected_specific_differential_phase'], 
        gridding_algo= 'map_gates_to_grid', # 'map_to_grid', # 
        grid_origin = (lat_mid, lon_mid), grid_origin_alt = alt_mid,
        weighting_function='BARNES')

CHIVO_grid = pyart.map.grid_from_radars(
        (CHIVO,), gatefilters=(gatefilter_CHIVO,), grid_shape=(21, 351, 351),
        grid_limits=((1200, 11200.0),(-35000, 35000), (-35000, 35000)),
        fields=['corrected_reflectivity', 'corrected_differential_reflectivity', 'corrected_cross_correlation_ratio', 'corrected_specific_differential_phase'],
        gridding_algo= 'map_gates_to_grid', # 'map_to_grid', # 
        grid_origin = (lat_mid, lon_mid), grid_origin_alt = alt_mid,
        weighting_function='BARNES')



#%%Display Grid

# parameters
sufix = '2.7 km 2019-01-25 19:30 UTC \n '
level = 3

vmin_Z =  0
vmax_Z =  70
cmap_Z = pyart.graph.cm.NWSRef

vmin_Zrd =  0
vmax_Zdr =  10
cmap_Zdr = pyart.graph.cm.RefDiff

vmin_rhohv =  0.7
vmax_rhohv=  1.0
cmap = 'viridis' 

vmin_Kdp =  0
vmax_Kdp =  10

def plot_grid(grid, field, title, vmax, vmin, level, cmap):
    display = pyart.graph.GridMapDisplay(grid)
    fig = plt.figure(figsize=[8, 8])
    display.plot_basemap(lon_lines = np.arange(-100, 0, 0.1))
    display.plot_basemap(lat_lines = np.arange(-100, 0, 0.1)) 
    display.plot_grid(field, level=level,  
                      vmin= vmin, vmax= vmax, title = title,
                      cmap = cmap)
    display.plot_crosshairs(CSAPR.longitude['data'][0],CSAPR.latitude['data'][0])
    display.plot_crosshairs(CHIVO.longitude['data'][0],CHIVO.latitude['data'][0])
    

# CSAPR -----------------------------------

#   Reflectivity
title = 'CSAPR ' + sufix + 'Corrected Reflectivity'
plot_grid(CSAPR, 'corrected_reflectivity', title, vmax_Z, vmin_Z, level, cmap)
#display = pyart.graph.GridMapDisplay(CSAPR_grid)
#fig = plt.figure(figsize=[8, 8])
#csapr_title = 'CSAPR-2 ' + sufix + field
#display.plot_basemap(lon_lines = np.arange(-100, 0, 0.1))
#display.plot_basemap(lat_lines = np.arange(-100, 0, 0.1)) 
#display.plot_grid('corrected_reflectivity', level=level,  
#                  vmin= vmin, vmax= vmax, title = csapr_title)#,
#                  #cmap = cmap)
#display.plot_crosshairs(CSAPR.longitude['data'][0],CSAPR.latitude['data'][0])
#display.plot_crosshairs(CHIVO.longitude['data'][0],CHIVO.latitude['data'][0])


#   Specific Differential Phase
display = pyart.graph.GridMapDisplay(CSAPR_grid)
fig = plt.figure(figsize=[8, 8])
csapr_title = 'CSAPR-2 ' + sufix
display.plot_basemap(lon_lines = np.arange(-100, 0, 0.1))
display.plot_basemap(lat_lines = np.arange(-100, 0, 0.1)) 
display.plot_grid('corrected_specific_differential_phase', level=level,  
                  vmin= vmin, vmax= vmax, title = csapr_title)#,
                  #cmap = cmap)
display.plot_crosshairs(CSAPR.longitude['data'][0],CSAPR.latitude['data'][0])
display.plot_crosshairs(CHIVO.longitude['data'][0],CHIVO.latitude['data'][0])


# CHIVO  ------------------------------

display = pyart.graph.GridMapDisplay(CHIVO_grid)
fig = plt.figure(figsize=[8, 8])

display.plot_basemap(lon_lines = np.arange(-100, 0, 0.1))
display.plot_basemap(lat_lines = np.arange(-100, 0, 0.1))
chivo_title  = 'CSU-CHIVO ' + sufix

display.plot_grid('corrected_specific_differential_phase', level=level,  
                  vmin= vmin, vmax= vmax, title = chivo_title)#, 
                  #cmap = cmap)

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




