#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 11:57:45 2021

@author: idariash
"""

import pyart
import chivo_grid_utl as utl
import numpy as np


import pyart
import time
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import netCDF4

import scipy.ndimage as spyi
import matplotlib.colors as mcolors
import cmocean
import pandas as pd

file = '/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/level_1b.2/2018/12/14/chivo.1b.2.20181214_013054.REL_PFAR360.nc'
'/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/level_1b.2/2018/12/14/chivo.1b.2.20181214_020047.REL_PFAR360.nc'
'/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/level_1b.2/2018/12/14/chivo.1b.2.20181214_023049.REL_PNL360A.nc'
'/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/level_1b.2/2018/12/14/chivo.1b.2.20181214_024046.REL_PNL360A.nc'
radar = pyart.io.read(file)
grid = utl.grid_radar(radar)
#%%
y, z, x = np.meshgrid(0.001*grid.y['data'],  0.001*grid.z['data'],  0.001*grid.x['data'])
reflectivity = grid.fields['corrected_reflectivity']['data']
RainRate = grid.fields['RainRate']['data']
HydroClass = grid.fields['HydroClass']['data'][2]
#reflectivity = np.swapaxes(reflectivity, 0, 2)
reflectivity = spyi.gaussian_filter(reflectivity, sigma = 4)

z = np.ma.masked_where( reflectivity < 18, z)
echo_top = np.amax(z, axis=0)
echo_top = np.ma.masked_where( echo_top < 4, echo_top)


#reflectivity = np.swapaxes(reflectivity, 0, 2)
max_ref = np.amax(reflectivity, axis=0)
RainRate = np.amax(RainRate, axis=0)

#%%





# table_filename = '/net/nasstore/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Climatology/relampago_echoTop_distribution_noBirdBath.xlsx'
# radar_grid_filename = '/net/nasstore/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Mountain/chivo_grid_20181214T0130.nc'
topo_filename = '/net/k2/storage/people/idariash/home/CSU/RELAMPAGO/analysis/Mountain/terr.nc'

# fig = utl.plot_grid_echoTop(table_filename, radar_grid_filename, topo_filename)

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

#(lat_0 - delta < lats) & (lats < lat_0 + 1)
lons_mask = np.absolute(lons - lon_0) < delta
#(lon_0 - delta < lons) & (lons < lon_0 + 1)  

Sierras_terr = terr[lats_mask, :]
Sierras_terr = Sierras_terr[:, lons_mask]
Sierras_terr = np.ma.masked_where( Sierras_terr < 0.5, Sierras_terr)

Sierras_lons = lons[lons_mask]
Sierras_lats = lats[lats_mask]

Sierras_lons_grid, Sierras_lats_grid = np.meshgrid(Sierras_lons, Sierras_lats)

Sierras_terr = spyi.gaussian_filter(Sierras_terr, sigma = 2)

#%% 


fig = plt.figure(figsize=(6, 5))



m = Basemap(projection='lcc', resolution=None,
            width=300E3, height=300E3, 
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
echo_top = np.ma.masked_where( R > 80, echo_top)
RainRate = np.ma.masked_where( R > 80, RainRate)
RainRate = np.ma.masked_where( RainRate < 1, RainRate)
HydroClass = np.ma.masked_where(HydroClass < 10, HydroClass)
HydroClass = np.ma.masked_where(HydroClass > 11, HydroClass)
HydroClass = np.ma.masked_where(R > 80, HydroClass)
HydroClass = np.ma.masked_where(max_ref <  40, HydroClass)
  
cs = m.pcolormesh(chivo_lons, chivo_lats, HydroClass, vmin=0, vmax=100, latlon=True,
        cmap=cmocean.cm.rain)#cmap = pyart.graph.cm.NWSRef)#cmap = pyart.graph.cm.LangRainbow12)#cmap=cmocean.cm.rain)#, norm=norm_precip)
   
m.colorbar(cs, label='Echo top height (km)')

time_start = netCDF4.num2date(grid.time['data'][0], grid.time['units'])
time_text = time_start.strftime('%Y-%m-%dT%H:%M UTC')


plt.title('CSU-CHIVO | '+ time_text+ ' \n Echo top height')


