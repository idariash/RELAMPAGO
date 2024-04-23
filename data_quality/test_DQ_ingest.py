#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 19:32:34 2019

@author: idariash
"""
import pyart
import netCDF4
import datetime
from matplotlib import pyplot as plt

def Zdr_correction(radar, Zdr_bias):
    radar.fields['differential_reflectivity']['data'] = radar.fields['differential_reflectivity']['data'] - Zdr_bias
    
def azimuth_correction(radar, azimuth_error):
    radar.azimuth['data'] = radar.azimuth['data'] - azimuth_error
    

chivo = '/net/denali/storage2/radar2/tmp/RELAMPAGO/RAW/2019/01/25/COL20190125_193048'
csapr = '/net/denali/storage2/radar2/tmp/Ivan/DOE/RELAMPAGO/CSAPR_PPI/corcsapr2cfrppiM1.a1.20190125.193003.nc'

CHIVO = pyart.io.read(chivo)
CSAPR =  pyart.io.read(csapr)

#Correct Azimuth and Zdr
azimuth_correction(CHIVO, 13.5)
Zdr_correction(CHIVO, 0.75)
Zdr_correction(CSAPR, -3.65)

#Correct Velocity CSAPR
csapr_fields = CSAPR.fields
gatefilter = pyart.filters.GateFilter(CSAPR)
gatefilter.exclude_transition()
gatefilter.exclude_invalid('mean_doppler_velocity')
gatefilter.exclude_invalid('reflectivity')
gatefilter.exclude_outside('reflectivity', 0, 80)
gatefilter.exclude_outside('normalized_coherent_power', 0.2, 1)
gatefilter.exclude_outside('copol_correlation_coeff', 0.8, 1)

nyq = CSAPR.instrument_parameters['nyquist_velocity']['data'][0]
last_radar = CSAPR
corr_vel = pyart.correct.dealias_region_based(
    last_radar, vel_field='mean_doppler_velocity', keep_original=False, 
    gatefilter = gatefilter, nyquist_vel=nyq, centered = True)
last_radar.add_field('corrected_velocity', corr_vel, replace_existing = True)

# perform dealiasing 4DD
dealias_data = pyart.correct.dealias_fourdd(
    CSAPR, vel_field = 'mean_doppler_velocity', last_radar = last_radar, gatefilter=gatefilter)
CSAPR.add_field('corrected_velocity', dealias_data, replace_existing = True)

#Create new NetCDF file
output_path = '/net/denali/storage2/radar2/tmp/Ivan/CSU/CHIVO/res/data/test_data_quality'

filename = output_path + '/corcsapr2cfrppiM1.a1.20190125.193003.nc'
pyart.io.write_cfradial(filename, CSAPR, format='NETCDF4', time_reference=True, arm_time_variables=False)

print('end here')

#Correct Velocity CHIVO,I can maka a function for this

gatefilter = pyart.filters.GateFilter(CHIVO)
gatefilter.exclude_transition()
gatefilter.exclude_invalid('velocity')
gatefilter.exclude_invalid('reflectivity')
gatefilter.exclude_outside('reflectivity', 0, 80)
gatefilter.exclude_outside('normalized_coherent_power', 0.2, 1)
gatefilter.exclude_outside('cross_correlation_ratio', 0.8, 1)

nyq = CHIVO.instrument_parameters['nyquist_velocity']['data'][0]
last_radar = CHIVO
corr_vel = pyart.correct.dealias_region_based(
    last_radar, vel_field='velocity', keep_original=False, 
    gatefilter = gatefilter, nyquist_vel=nyq, centered = True)
last_radar.add_field('corrected_velocity', corr_vel, replace_existing = True)

# perform dealiasing 4DD
dealias_data = pyart.correct.dealias_fourdd(
    CHIVO, last_radar = last_radar, gatefilter=gatefilter)
CHIVO.add_field('corrected_velocity', dealias_data, replace_existing = True)

#Create new NetCDF file
output_path = '/net/denali/storage2/radar2/tmp/Ivan/CSU/CHIVO/res/data/test_data_quality'
CHIVO.metadata['title'] = 'Level 1 CHIVO Radar RELAMPAGO'
CHIVO.metadata['institution'] = 'Colorado State University'
CHIVO.metadata['reference'] = 'RELAMPAGO CHIVO report and Data quality control of CHIVO report'
CHIVO.metadata['comment'] = 'Data collected during RELAMPAGO field campaign at Lozada site in Cordoba, Argentina'
CHIVO.metadata['instrument_name'] = 'CSU-CHIVO'
CHIVO.metadata['history'] = 'Created by Ivan Arias from RAW sigmet data, July 2019'

scan_strategy = str(CHIVO.metadata['sigmet_task_name'] )
scan_strategy = scan_strategy[2 : len(scan_strategy) - 2]
time_start = netCDF4.num2date(CHIVO.time['data'][0], CHIVO.time['units'])
time_text = time_start.strftime('%Y%m%d_%H%M%S')
filename = output_path + '/chivo.1a.' + time_text + '.' + scan_strategy + '.nc'
pyart.io.write_cfradial(filename, CHIVO, format='NETCDF4', time_reference=True, arm_time_variables=False)

print(scan_strategy)
print('end of test')

# filter out gates with differential_reflectivity > 100 from both radars
gatefilter_CHIVO = pyart.filters.GateFilter(CHIVO)
gatefilter_CHIVO.exclude_transition()
gatefilter_CHIVO.exclude_outside('reflectivity', 10, 100)
gatefilter_CHIVO.exclude_outside('normalized_coherent_power', 0.3, 1)

gatefilter_CSAPR = pyart.filters.GateFilter(CSAPR)
gatefilter_CSAPR.exclude_transition()
gatefilter_CSAPR.exclude_outside('reflectivity', 10, 100)
gatefilter_CSAPR.exclude_outside('normalized_coherent_power', 0.3, 1)
nyq_csapr = CSAPR.instrument_parameters['nyquist_velocity']['data'][0]
nyq_chivo = CHIVO.instrument_parameters['nyquist_velocity']['data'][0]

csapr_fields = CSAPR.fields


lat_mid = (CHIVO.latitude['data'][0] + CSAPR.latitude['data'][0])/2

lon_mid = (CHIVO.longitude['data'][0] + CSAPR.longitude['data'][0])/2

alt_mid = (CHIVO.altitude['data'][0] + CSAPR.altitude['data'][0])/2

grids = pyart.map.grid_from_radars(
        (CSAPR,), gatefilters=(gatefilter_CSAPR,), grid_shape=(41, 251, 251),
        grid_limits=((alt_mid, 20000.0),(-150000, 150000), (-150000, 150000)),
        fields=['differential_reflectivity'], gridding_algo= 'map_gates_to_grid', # 'map_to_grid', # 
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
display.plot_grid('differential_reflectivity', level=level,  
                  vmin= -2, vmax= 6, 
                  cmap = pyart.graph.cm.RefDiff)
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
display.plot_grid('differential_reflectivity', level=level, 
                  vmin= -2, vmax= 6, 
                  cmap = pyart.graph.cm.RefDiff)
display.plot_crosshairs(lon=lon, lat=lat)
display.plot_crosshairs(CSAPR.longitude['data'][0],CSAPR.latitude['data'][0])


grids = pyart.map.grid_from_radars(
        (CHIVO,), gatefilters=(gatefilter_CHIVO,), grid_shape=(41, 251, 251),
        grid_limits=((alt_mid, 20000.0),(-150000, 150000), (-150000, 150000)),
        fields=['differential_reflectivity'], gridding_algo= 'map_gates_to_grid', # 'map_to_grid', # 
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
display.plot_grid('differential_reflectivity', level=level,  
                  vmin= -2, vmax= 6, 
                  cmap = pyart.graph.cm.RefDiff)
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
display.plot_grid('differential_reflectivity', level=level, 
                  vmin= -2, vmax= 6, 
                  cmap = pyart.graph.cm.RefDiff)
display.plot_crosshairs(lon=lon, lat=lat)
display.plot_crosshairs(CSAPR.longitude['data'][0],CSAPR.latitude['data'][0])

## panel 2, longitude slice.
#ax2 = fig.add_axes(x_cut_panel_axes)
#display.plot_longitude_slice('differential_reflectivity', lon=lon, lat=lat, vmin=vmin, vmax=vmax,
#                             )
#ax2.set_ylim([0,15])
#ax2.set_xlim([-150,150])
##ax2.set_xlabel('Distance from SGP CF (km)')
#
## panel 3, latitude slice
#ax3 = fig.add_axes(y_cut_panel_axes)
#ax3.set_ylim([0,15])
#ax3.set_xlim([-150,150])
#display.plot_latitude_slice('differential_reflectivity', lon=lon, lat=lat, vmin=vmin, vmax=vmax,
#                            )
#
##
### panel 1, basemap, radar differential_reflectivity and NARR overlay
##ax1 = fig.add_axes(map_panel_axes)
##display.plot_basemap(lon_lines = np.arange(-104, -93, 2) )
##display.plot_grid('differential_reflectivity', level=level, vmin=vmin, vmax=vmax,
##                  )
##display.plot_crosshairs(lon=lon, lat=lat)
##
##
### panel 2, longitude slice.
##ax2 = fig.add_axes(x_cut_panel_axes)
##display.plot_longitude_slice('differential_reflectivity', lon=lon, lat=lat, vmin=vmin, vmax=vmax,
##                             )
##ax2.set_ylim([0,15])
##ax2.set_xlim([-400,400])
##ax2.set_xlabel('Distance from SGP CF (km)')
##
### panel 3, latitude slice
##ax3 = fig.add_axes(y_cut_panel_axes)
##ax3.set_ylim([0,15])
##ax3.set_xlim([-200,200])
##display.plot_latitude_slice('differential_reflectivity', lon=lon, lat=lat, vmin=vmin, vmax=vmax,
##                            )
