#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 14:56:45 2022

@author: idariash
"""

import sys
sys.path.append('/net/nasstore/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/plot')
import chivo_grid_utl as utl
import pyart
import glob
import netCDF4
import numpy as np
from operator import and_
import multiprocessing as mp
import gc
import os

def grid_EchoTop(radar, grid_shape=(20, 801, 801), xlim=(-150000, 150000), 
               ylim=(-150000, 150000), zlim=(1000, 20000),
               fields=['corrected_reflectivity', 'corrected_differential_reflectivity', 
                       'RainRate', 'HydroClass'], origin=None):
  
    
    reflectivity = (radar.fields['corrected_reflectivity']['data'].copy())
    reflectivity.mask = False
    reflectivity = np.ma.masked_where(reflectivity < 100, reflectivity)
    reflectivity.filled(20)
    reflectivity.mask = False
    for i in range(0, len(reflectivity)):
        for j in range(0, len(reflectivity[i])):
            reflectivity[i][j] = 20
    
    
    
    radar.add_field_like('corrected_reflectivity', 'corrected_reflectivity', 
                         reflectivity, replace_existing=True)
    
    gatefilter = pyart.filters.GateFilter(radar)
    # gatefilter.exclude_transition()
    
    # gatefilter.exclude_invalid('corrected_reflectivity')
    # gatefilter.exclude_invalid('corrected_differential_phase')
    # gatefilter.exclude_outside('corrected_reflectivity', 10, 100)
    # try:
    #     gatefilter.exclude_outside('normalized_coherent_power', 0.3, 1)
    # except:
    #     print('no NCP')
    # gatefilter.exclude_outside('corrected_copol_correlation_ratio', 0.9, 1)

    
    if origin is None:
        origin = (radar.latitude['data'][0],
                  radar.longitude['data'][0])
        
    radar_lower_sweep = radar#radar.extract_sweeps([0,1])
        
    grid = pyart.map.grid_from_radars(
        radar_lower_sweep, grid_shape=grid_shape,
        grid_limits=(zlim, ylim, xlim),
        grid_origin=origin, fields=fields,gatefilters = gatefilter,
        gridding_algo='map_gates_to_grid', grid_origin_alt=0, 
        min_radius = 500, nb = 1, bsp = 1, )#weighting_function = 'Nearest')

    
    return grid

file = '/net/k2/storage/projects/CSAPR/level_1b.2/corcsapr2cfrppiqcM1.1b.20181214.020004.nc'

'/net/k2/storage/projects/RELAMPAGO/RMA1/quality_controlled/level_1b/cfrad.20181214_020509.000_to_20181214_021157.928_1_SUR.nc'

'/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/level_1b.2/2018/12/14/chivo.1b.2.20181214_020047.REL_PFAR360.nc'

radar = pyart.io.read(file)
grid = grid_EchoTop(radar = radar,  origin = [-31.6342, -64.1686])
ref_EchoTop = grid.fields['corrected_reflectivity']['data']

# Max. 18 dBZ reflectivity
y, z, x = np.meshgrid(0.001*grid.y['data'], 0.001*grid.z['data'],  0.001*grid.x['data'])

z_EchoTop = np.ma.masked_where(ref_EchoTop < 18, z)

EchoTop = np.amax(z_EchoTop, axis=0)

#%% Create figures
    
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from netCDF4 import Dataset
import netCDF4
import scipy.ndimage as spyi
import cmocean
from mpl_toolkits.basemap import Basemap


# Read topography

topo_filename = '/net/k2/storage/people/idariash/home/CSU/RELAMPAGO/analysis/Mountain/terr.nc'
topography = Dataset(topo_filename)

lons = topography.variables['X'][:]
lats = topography.variables['Y'][:]
terr = 0.001*topography.variables['topo'][:]


height_toCut = 2
reflectivity_toCut = 5
range_toCut = 200
 

lon_0 = grid.origin_longitude['data']
lat_0 = grid.origin_latitude['data']

delta = 6#deg
lats_mask = np.absolute(lats - lat_0) < delta

lons_mask = np.absolute(lons - lon_0) < delta

Sierras_terr = terr[lats_mask, :]
Sierras_terr = Sierras_terr[:, lons_mask]
Sierras_terr = np.ma.masked_where( Sierras_terr < 0.5, Sierras_terr)

Sierras_lons = lons[lons_mask]
Sierras_lats = lats[lats_mask]

Sierras_lons_grid, Sierras_lats_grid = np.meshgrid(Sierras_lons, Sierras_lats)

Sierras_terr = spyi.gaussian_filter(Sierras_terr, sigma = 2)

#%%

m = Basemap(projection='lcc', resolution=None,
            width=300E3, height=300E3, 
            lat_0=lat_0, lon_0=lon_0)

#m.etopo(scale=2, alpha=0.5)

parallels = np.arange(-40,-20,1)
# labels = [left,right,top,bottom]
m.drawparallels(parallels,labels=[True,False,False, False])

meridians = np.arange(-70,-60.,1)
m.drawmeridians(meridians,labels=[True,False,False,True])

x_Cordoba, y_Cordoba = m(-64.19, -31.42)
x_CHIVO, y_CHIVO = m(-64.1686, -31.6342)
x_ARM, y_ARM = m(-64.7284, -32.1264)
x_Rio4to, y_Rio4to = m(-64.3493, -33.1232)

plt.plot(x_CHIVO, y_CHIVO, 'ok', markersize=5)
plt.text(x_CHIVO, y_CHIVO, ' C', fontsize=12);

plt.plot(x_ARM, y_ARM, 'ok', markersize=5)
plt.text(x_ARM, y_ARM, ' Y', fontsize=12);

plt.plot(x_Cordoba, y_Cordoba, 'ok', markersize=5)
plt.text(x_Cordoba, y_Cordoba, ' R', fontsize=12);

levels_terr = [1, 2]

m.contour(Sierras_lons_grid, Sierras_lats_grid, Sierras_terr, levels_terr, 
               latlon = True, linewidths = 1, colors = 'k',
                      linestyles='solid', antialiased=True, animated=True)



lon_0 = grid.origin_longitude['data']
lat_0 = grid.origin_latitude['data']
chivo_lons, chivo_lats = grid.get_point_longitude_latitude()

x, y = np.meshgrid(0.001*grid.x['data'],  0.001*grid.y['data'])
R = np.sqrt(x ** 2 + y ** 2)

# RainRate = np.ma.masked_where( R > 100, RainRate)
# RainRate = np.ma.masked_where( RainRate < 1, RainRate)
# HydroClass = np.ma.masked_where(HydroClass < 10, HydroClass)
# HydroClass = np.ma.masked_where(HydroClass > 11, HydroClass)
# HydroClass = np.ma.masked_where(R > 80, HydroClass)
# HydroClass = np.ma.masked_where(max_ref <  40, HydroClass)

cs = m.pcolormesh(chivo_lons, chivo_lats, EchoTop, vmin=0, vmax=20, latlon=True,
        cmap='jet')#, norm=colors.LogNorm()) #  cmocean.cm.rain #pyart.graph.cm.LangRainbow12

#m.colorbar(cs, label='Total Rain (mm)')
m.colorbar(cs, label='Echo Top Heigh (km)')

# time_start = netCDF4.num2date(grid.time['data'][0], grid.time['units'])
# time_text = time_start.strftime('%Y-%m-%dT%H:%M')

radar_name = 'CSU-CHIVO'


#plt.title('CSU-CHIVO | '+ '2018/11/10 to 2019/01/31' + ' \n Accumulated Rain Map')
plt.title(radar_name + ' | '+ '2018/12/14 02:30 UTC')


#accumulated_total_rain = sum(results)

#%%
# import chivo_plot as chivo_plot

# file = '/net/k2/storage/projects/RELAMPAGO/RMA1/quality_controlled/level_1b/cfrad.20181214_032757.000_to_20181214_032919.922_1_SUR.nc'

# chivo_plot.ppi_drops(file, 0)

