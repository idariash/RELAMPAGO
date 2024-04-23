#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 08:54:11 2022

@author: idariash

Plot the rain from the radars and the gauges
"""

import pandas as pd
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
from datetime import datetime

rainData_file = '/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Rain/Rain Dec 14 gauges.xlsx'

df = pd.read_excel(rainData_file, sheet_name = 1)

# lats = df['Lat'][1]

# print(df)

# Create figures
    
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

lon_0 = -64.1686
lat_0 = -31.6342

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

rain_X = np.zeros([len(df)]); 
rain_Y = np.zeros([len(df)]);
rain_values = np.zeros([len(df)]);

for i in range(0, len(df)):
    x, y = m(df['Lon'][i]/1E3, df['Lat'][i]/1E3)
    rain_X[i] = x
    rain_Y[i] = y
    rain_values[i] = df['Rain_UTC'][i]
    
cs = plt.scatter(rain_X, rain_Y, c = rain_values, vmin=1, vmax=200, 
                 cmap=pyart.graph.cm.LangRainbow12, norm=colors.LogNorm())
m.colorbar(cs, label='Rain (mm)')

# lon_0 = grid.origin_longitude['data']
# lat_0 = grid.origin_latitude['data']
# chivo_lons, chivo_lats = grid.get_point_longitude_latitude()

# x, y = np.meshgrid(0.001*grid.x['data'],  0.001*grid.y['data'])
# R = np.sqrt(x ** 2 + y ** 2)

# # RainRate = np.ma.masked_where( R > 100, RainRate)
# # RainRate = np.ma.masked_where( RainRate < 1, RainRate)
# # HydroClass = np.ma.masked_where(HydroClass < 10, HydroClass)
# # HydroClass = np.ma.masked_where(HydroClass > 11, HydroClass)
# # HydroClass = np.ma.masked_where(R > 80, HydroClass)
# # HydroClass = np.ma.masked_where(max_ref <  40, HydroClass)

# hour = 2
# accumulated_total_rain = rma_rain_per_hour[hour]/(rma_number_of_files_per_hour[hour]*2)

# accumulated_total_rain = np.ma.masked_where(accumulated_total_rain <  2, accumulated_total_rain)
# accumulated_total_rain = np.ma.masked_where(R > 130, accumulated_total_rain)
# accumulated_total_rain = rain #rain_plot#rain_std#rain

# cs = m.pcolormesh(chivo_lons, chivo_lats, accumulated_total_rain, vmin=1, vmax=200, latlon=True,
#         cmap=pyart.graph.cm.LangRainbow12, norm=colors.LogNorm()) #  cmocean.cm.rain #pyart.graph.cm.LangRainbow12

# #m.colorbar(cs, label='Total Rain (mm)')
# m.colorbar(cs, label='Rain Rate (mm/h)')

# # time_start = netCDF4.num2date(grid.time['data'][0], grid.time['units'])
# # time_text = time_start.strftime('%Y-%m-%dT%H:%M')

# radar_name = 'CSU-CHIVO'

# if 'csapr' in DataPath:
#     radar_name = 'CSAPR-2'
    
# if 'rma' in DataPath:
#     radar_name = 'RMA-1'

# #plt.title('CSU-CHIVO | '+ '2018/11/10 to 2019/01/31' + ' \n Accumulated Rain Map')
# plt.title(radar_name + ' | '+ '2018/12/14 02:30 UTC')


#accumulated_total_rain = sum(results)

#%%
# import chivo_plot as chivo_plot

# file = '/net/k2/storage/projects/RELAMPAGO/RMA1/quality_controlled/level_1b/cfrad.20181214_032757.000_to_20181214_032919.922_1_SUR.nc'

# chivo_plot.ppi_drops(file, 0)

file = '/net/nasstore/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Climatology/RadarsRain_Dec14.npz'
radarsRain = np.load(file)['radarsRain']
radar_lons = np.load(file)['chivo_lons']
radar_lats = np.load(file)['chivo_lats']

from scipy import interpolate

Lons = np.mean(radar_lons, 0)
Lats = np.mean(radar_lats, 1)

f = interpolate.interp2d(Lons, Lats, radarsRain, kind='cubic')

radarsRain_atGauge = np.zeros(len(df))

for i in range(0, len(df)):
    radarsRain_atGauge[i] = f(df['Lon'][i]/1E3, df['Lat'][i]/1E3)
  
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
    
cs = plt.scatter(rain_X, rain_Y, c = radarsRain_atGauge, vmin=1, vmax=200, 
                 cmap=pyart.graph.cm.LangRainbow12, norm=colors.LogNorm())
m.colorbar(cs, label='Rain (mm)')

#%%
rain_values = np.ma.masked_where(rain_values < 10, rain_values)

plt.scatter(rain_values, radarsRain_atGauge)
plt.xlim(xmin=-5,xmax=100)
plt.ylim(ymin=-5,ymax=100)



