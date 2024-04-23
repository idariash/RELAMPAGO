#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 15:45:55 2019

@author: idariash
"""

#%% Import, define and read
import pyart
import matplotlib.pyplot as plt

def Zdr_correction(radar, Zdr_bias):
    radar.fields['differential_reflectivity']['data'] = radar.fields['differential_reflectivity']['data'] - Zdr_bias
    
def azimuth_correction(radar, azimuth_error):
    radar.azimuth['data'] = radar.azimuth['data'] - azimuth_error

chivo = '/net/denali/storage2/radar2/tmp/RELAMPAGO/RAW/2019/01/26/COL20190126_053008'
#'/net/denali/storage2/radar2/tmp/RELAMPAGO/RAW/2019/01/25/COL20190125_193048' 
#'/net/denali/storage2/radar2/tmp/RELAMPAGO/RAW/2019/01/26/COL20190126_053008'
#'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Attenuation/NetCDF/CHIVO/cfrad.20190126_053009.397_to_20190126_053642.486_col-radar_PPINEARLT360_SUR.nc'
csapr = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Attenuation/NetCDF/CSAPR/corcsapr2cfrppiM1.a1.20190126.053003.nc'
#'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Attenuation/NetCDF/CSAPR/corcsapr2cfrppiM1.a1.20190126.053003.nc'

CHIVO = pyart.io.read(chivo)
CSAPR =  pyart.io.read(csapr)

#%%Correct Azimuth and Zdr
azimuth_correction(CHIVO, 13.5)
Zdr_correction(CHIVO, 0.75)
Zdr_correction(CSAPR, -3.65)
a = CHIVO.fields 



#%% filter out gates with differential_reflectivity > 100 from both radars
gatefilter_CHIVO = pyart.filters.GateFilter(CHIVO)
gatefilter_CHIVO.exclude_transition()
gatefilter_CHIVO.exclude_outside('reflectivity', 10, 100)
gatefilter_CHIVO.exclude_outside('normalized_coherent_power', 0.6, 1)
gatefilter_CHIVO.exclude_outside('cross_correlation_ratio', 0.9, 1)

gatefilter_CSAPR = pyart.filters.GateFilter(CSAPR)
gatefilter_CSAPR.exclude_transition()
gatefilter_CSAPR.exclude_outside('reflectivity', 10, 100)
gatefilter_CSAPR.exclude_outside('normalized_coherent_power', 0.3, 1)
#60 = CSAPR.instrument_parameters['nyquist_reflectivity']['data'][0]
#60 = CHIVO.instrument_parameters['nyquist_reflectivity']['data'][0]

csapr_fields = CSAPR.fields


lat_mid = (CHIVO.latitude['data'][0] + CSAPR.latitude['data'][0])/2

lon_mid = (CHIVO.longitude['data'][0] + CSAPR.longitude['data'][0])/2

alt_mid = (CHIVO.altitude['data'][0] + CSAPR.altitude['data'][0])/2

#%%-----------------CSAPR Grid ----------------------------

grids = pyart.map.grid_from_radars(
        (CSAPR,), gatefilters=(gatefilter_CSAPR,), grid_shape=(41, 201, 201),
        grid_limits=((alt_mid, 20000.0),(-50000, 50000), (-50000, 50000)),
        fields=['reflectivity', 'differential_reflectivity'], gridding_algo= 'map_gates_to_grid', # 'map_to_grid', # 
        grid_origin = (lat_mid, lon_mid), grid_origin_alt = alt_mid,
        weighting_function='BARNES')
#%%

display = pyart.graph.GridMapDisplay(grids)
#fig = plt.figure(figsize=[15, 7])
fig = plt.figure(figsize=[15, 15])

# panel sizes
map_panel_axes = [0.05, 0.05, .4, .80]
x_cut_panel_axes = [0.55, 0.10, .4, .25]
y_cut_panel_axes = [0.55, 0.50, .4, .25]

# parameters
level = 3
level_1 = 20 
vmin = 0
vmax = 6
lat = -31.6342 #36.5
lon = -64.1586 #-98.0

# panel 1, basemap, radar differential_reflectivity and NARR overlay
#ax1 = fig.add_axes(map_panel_axes)
#display.plot_basemap(lon_lines = np.arange(-104, -93, 1) )
display.plot_grid('reflectivity', level=level,  
                  vmin= 0, vmax= 60, 
                  cmap = pyart.graph.cm.NWSRef)
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
#level = 3
vmin = 0
vmax = 6
lat = -31.6342 #36.5
lon = -64.1586 #-98.0

# panel 1, basemap, radar differential_reflectivity and NARR overlay
#ax1 = fig.add_axes(map_panel_axes)
#display.plot_basemap(lon_lines = np.arange(-104, -93, 1) )
display.plot_grid('differential_reflectivity', level=level, 
                  vmin= -2, vmax= 6, 
                  cmap = pyart.graph.cm.RefDiff)
display.plot_crosshairs(lon=lon, lat=lat)
display.plot_crosshairs(CSAPR.longitude['data'][0],CSAPR.latitude['data'][0])


#%%-----------------CHIVO Grid ----------------------------

grids = pyart.map.grid_from_radars(
        (CHIVO,), gatefilters=(gatefilter_CHIVO,), grid_shape=(41, 201, 201),
        grid_limits=((alt_mid, 20000.0),(-50000, 50000), (-50000, 50000)),
        fields=['reflectivity', 'differential_reflectivity'], gridding_algo= 'map_gates_to_grid', # 'map_to_grid', # 
        grid_origin = (lat_mid, lon_mid), grid_origin_alt = alt_mid,
        weighting_function='BARNES')
#%%

display = pyart.graph.GridMapDisplay(grids)
#fig = plt.figure(figsize=[15, 7])
fig = plt.figure(figsize=[15, 15])

# panel sizes
map_panel_axes = [0.05, 0.05, .4, .80]
x_cut_panel_axes = [0.55, 0.10, .4, .25]
y_cut_panel_axes = [0.55, 0.50, .4, .25]

# parameters
level = 3
level_1 = 20 
vmin = 0
vmax = 6
lat = -31.6342 #36.5
lon = -64.1586 #-98.0

# panel 1, basemap, radar differential_reflectivity and NARR overlay
#ax1 = fig.add_axes(map_panel_axes)
#display.plot_basemap(lon_lines = np.arange(-104, -93, 1) )
display.plot_grid('reflectivity', level=level,  
                  vmin= 0, vmax= 60, 
                  cmap = pyart.graph.cm.NWSRef)
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
#level = 3
vmin = 0
vmax = 6
lat = -31.6342 #36.5
lon = -64.1586 #-98.0

# panel 1, basemap, radar reflectivity and NARR overlay
#ax1 = fig.add_axes(map_panel_axes)
#display.plot_basemap(lon_lines = np.arange(-104, -93, 1) )
display.plot_grid('differential_reflectivity', level=level, 
                  vmin= -2, vmax= 6, 
                  cmap = pyart.graph.cm.RefDiff)
display.plot_crosshairs(lon=lon, lat=lat)
display.plot_crosshairs(CSAPR.longitude['data'][0],CSAPR.latitude['data'][0])
