#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 17:18:56 2020

@author: idariash
"""
import pyart
import time

import AMS_Mountain_utl as utl
import importlib

from mpl_toolkits.basemap import Basemap, shiftgrid, cm
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import netCDF4

import scipy.ndimage as spyi

filename = 'chivo_grid_20181214T0130.nc'
#filename = 'csapr_grid_20181214T01.nc'
height_toCut = 1.0
reflectivity_toCut = 20
range_toCut = 40
chivo_grid = pyart.io.read_grid(filename)

lon_0 = chivo_grid.origin_longitude['data']
lat_0 = chivo_grid.origin_latitude['data']

Z = chivo_grid.fields['corrected_reflectivity']['data']

Z = np.ma.masked_where( Z < reflectivity_toCut, Z)

x, y, z = np.meshgrid(0.001*chivo_grid.x['data'], 0.001*chivo_grid.y['data'],
                     0.001*chivo_grid.z['data'])

x = np.swapaxes(x,1,2)
x = np.swapaxes(x,0,1)

y = np.swapaxes(y,1,2)
y = np.swapaxes(y,0,1)

z = np.swapaxes(z,1,2)
z = np.swapaxes(z,0,1)

R = np.sqrt(x ** 2 + y ** 2 + z ** 2)

Z = np.ma.masked_where( z < height_toCut, Z)
Z = np.ma.masked_where( R > range_toCut, Z)


#Zdr_max = np.ma.masked_where( Z_max < 20, Zdr_max)

lon_0 = chivo_grid.origin_longitude['data']
lat_0 = chivo_grid.origin_latitude['data']
chivo_lons, chivo_lats = chivo_grid.get_point_longitude_latitude()

x, z = np.meshgrid(0.001*chivo_grid.x['data'], 0.001*chivo_grid.z['data'])

Z = np.swapaxes(Z,0,1)
Z = np.swapaxes(Z,0,2)

Z_max = np.amax(Z, axis=0)
lats_vector = np.mean(chivo_lats, axis=1)
lats_mat, z = np.meshgrid(lats_vector, 0.001*chivo_grid.z['data'])

fig = plt.figure(figsize=(7, 3))
cs = plt.pcolormesh(lats_mat, z, Z_max, vmin=0, vmax=70,# latlon=True,
             cmap='pyart_NWSRef')

plt.colorbar(cs, label='Reflectivity [dBZ]')

time_start = netCDF4.num2date(chivo_grid.time['data'][0], chivo_grid.time['units'])
time_text = time_start.strftime('%Y-%m-%dT%H:%M UTC')

plt.title('CSU-CHIVO | '+ time_text+ ' \n Max. Reflectivity')

plt.xlim(-32, -31.5)
#plt.xlim(-65, -63.5)
plt.ylim(0, 20)

plt.xlabel('Longitude (deg.)')
plt.ylabel('Hieght (km)')

time_text = time_start.strftime('%Y%m%dT%H%M')
fig_name = ('./fig/lon_maxRef_' + time_text + '.png')
fig.savefig(fig_name)