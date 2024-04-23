#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 17:25:59 2020

@author: idariash
"""

import pyart
import netCDF4
import numpy as np
import matplotlib.pyplot as plt

filename = 'chivo_grid_20181214T0140.nc'
chivo_grid = pyart.io.read_grid(filename)

Z = chivo_grid.fields['corrected_reflectivity']['data']

Z = np.ma.masked_where( Z < 20, Z)

x, y, z = np.meshgrid(0.001*chivo_grid.x['data'], 0.001*chivo_grid.y['data'],
                     0.001*chivo_grid.z['data'])

x = np.swapaxes(x,1,2)
x = np.swapaxes(x,0,1)

y = np.swapaxes(y,1,2)
y = np.swapaxes(y,0,1)

z = np.swapaxes(z,1,2)
z = np.swapaxes(z,0,1)

R = np.sqrt(x ** 2 + y ** 2 + z ** 2)

Z = np.ma.masked_where( R > 60, Z)

lon_0 = chivo_grid.origin_longitude['data']
lat_0 = chivo_grid.origin_latitude['data']
chivo_lons, chivo_lats = chivo_grid.get_point_longitude_latitude()

#Z = np.swapaxes(Z,1,2)
Z = np.swapaxes(Z,0,1)

Z_max = np.amax(Z, axis=0)
lons_vector = np.mean(chivo_lons, axis=0)
lons_mat, z = np.meshgrid(lons_vector, 0.001*chivo_grid.z['data'])

fig = plt.figure(figsize=(7, 3))
cs = plt.pcolormesh(lons_mat, z, Z_max, vmin=0, vmax=70,# latlon=True,
             cmap='pyart_NWSRef')

plt.colorbar(cs, label='Reflectivity [dBZ]')

time_start = netCDF4.num2date(chivo_grid.time['data'][0], chivo_grid.time['units'])
time_text = time_start.strftime('%Y-%m-%dT%H:%M UTC')

plt.title('CSU-CHIVO | '+ time_text+ ' \n Max. Reflectivity')

plt.xlim(-65, -64)
plt.ylim(0, 20)

plt.xlabel('Longitude [°]')
plt.ylabel('Hieght [km]')
