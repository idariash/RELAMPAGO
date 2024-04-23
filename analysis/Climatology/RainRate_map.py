#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 15:44:57 2020

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



DataPath = '/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/level_1b.2/2019/01'
files = [f for f in glob.glob(DataPath +"/**/*P*.nc", recursive=True)]
files.sort()
radar = pyart.io.read(files[0])
grid = utl.grid_radar(radar)
accumulated_total_rain  = (grid.fields['RainRate']['data'].copy())*0 # every 10 min per file
accumulated_total_rain = np.ma.masked_where(accumulated_total_rain < 0, accumulated_total_rain)
accumulated_total_rain.filled(0)
accumulated_total_rain.mask = False
accumulated_total_rain = np.amax(accumulated_total_rain, axis=0)


for file in files:
    try:
        print (file)
        radar = pyart.io.read(file)
        grid = utl.grid_radar(radar)
        total_rain_file = (grid.fields['RainRate']['data'].copy())/6 
        total_rain_file = np.ma.masked_where(total_rain_file < 0, total_rain_file)
        total_rain_file.filled(0)
        total_rain_file.mask = False
        total_rain_file = np.amax(total_rain_file, axis=0)
        
        accumulated_total_rain  = accumulated_total_rain + total_rain_file
        
    except:
        print('error: ' + file)
    
    
    # accumulated_total_rain_dict = {'data': accumulated_total_rain, 'units': 'mm',
    #                                'missing_value': grid.fields['RainRate']['missing_value'],
    #                                }
    # grid.add_field('accumulated_total_rain', accumulated_total_rain_dict, True)
    # time_start = netCDF4.num2date(grid.time['data'][0], grid.time['units'])
    # time_text = time_start.strftime('%Y%m%dT%H%M')
    # grid_fname = ('../rsc/out/chivo_grid_'+ time_text + '.nc')
    # grid.write(grid_fname)
    
#%%

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

m = Basemap(projection='lcc', resolution=None,
            width=150E3, height=150E3, 
            lat_0=lat_0, lon_0=lon_0)

#m.etopo(scale=2, alpha=0.5)

parallels = np.arange(-40,-20,0.5)
# labels = [left,right,top,bottom]
m.drawparallels(parallels,labels=[True,False,False, False])

meridians = np.arange(-70,-60.,0.5)
m.drawmeridians(meridians,labels=[True,False,False,True])

x_Cordoba, y_Cordoba = m(-64.19, -31.42)
x_CHIVO, y_CHIVO = m(-64.1686, -31.6342)
x_ARM, y_ARM = m(-64.7284, -32.1264)
x_Rio4to, y_Rio4to = m(-64.3493, -33.1232)

plt.plot(x_CHIVO, y_CHIVO, 'ok', markersize=5)
plt.text(x_CHIVO, y_CHIVO, ' C', fontsize=12);

# plt.plot(x_ARM, y_ARM, 'ok', markersize=5)
# plt.text(x_ARM, y_ARM, ' Y', fontsize=12);

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

cs = m.pcolormesh(chivo_lons, chivo_lats, accumulated_total_rain, vmin=1, vmax=10000, latlon=True,
        cmap=pyart.graph.cm.LangRainbow12, norm=colors.LogNorm())

m.colorbar(cs, label='Total Rain (mm)')

# time_start = netCDF4.num2date(grid.time['data'][0], grid.time['units'])
# time_text = time_start.strftime('%Y-%m-%dT%H:%M')


plt.title('CSU-CHIVO | '+ '2018/12/14' + ' UTC'+ ' \n Accumulated Rain')
