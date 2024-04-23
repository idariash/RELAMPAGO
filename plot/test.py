#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 17:39:57 2021

@author: idariash
"""

import scipy
import pyart

import numpy as np
from scipy.interpolate import interpn

from numpy import array
from scipy.interpolate import RegularGridInterpolator as rgi
import matplotlib.pyplot as plt
import netCDF4

filename = '/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Dec14/chivo_1b_20181214_0230.nc'
filename_dd = '/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Dec14/dda_3radars_20181214_0230.nc'

pyart_grid_1b = pyart.io.read_grid(filename)
pyart_grid_dd = pyart.io.read_grid(filename_dd)

x = 0.001*pyart_grid_1b.y['data']
y = 0.001*pyart_grid_1b.y['data']
z = 0.001*pyart_grid_1b.z['data']

d = 0

x0 = -90
y0 = -10

x1 = -60 + 20
y1 = 22  + 20

x0 = x0 + d
x1 = x1 + d

print(np.arctan((y1-y0)/(x1 - x0))*180/np.pi)

distance = np.sqrt((x1- x0)**2 + (y1 - y0)**2)


xi = np.linspace(x0, x1, int(distance))
yi = np.linspace(y0, y1, int(distance))


U = pyart_grid_dd.fields['u']['data']
V = pyart_grid_dd.fields['v']['data']
W = pyart_grid_dd.fields['w']['data']
Z = pyart_grid_1b.fields['corrected_reflectivity']['data']

U = np.swapaxes(U, 0, 2)
V = np.swapaxes(V, 0, 2)
W = np.swapaxes(W, 0, 2)
Z = np.swapaxes(Z, 0, 2)

fn_u =rgi((x, y, z), U)
fn_v =rgi((x, y, z), V)
fn_w =rgi((x, y, z), W)
fn_z =rgi((x, y, z), Z)

Ui_pj = np.zeros((len(yi), len(z))) # Projection over V1-V0
Vi_pj = np.zeros((len(yi), len(z)))
Wi = np.zeros((len(yi), len(z)))
Zi = np.zeros((len(yi), len(z)))

for i in range(len(yi),):
    for j in range(len(z)):
        u = fn_u(np.array([xi[i], yi[i], z[j]]))
        v = fn_v(np.array([xi[i], yi[i], z[j]]))
        Ui_pj[i, j] = ((x1-x0)*u + (y1-y0)*v)/distance
        Wi[i, j] = fn_w(np.array([xi[i], yi[i], z[j]]))
        Zi[i, j] = fn_z(np.array([xi[i], yi[i], z[j]]))
        
Ui_pj = np.ma.masked_where( Ui_pj < -100, Ui_pj)
Zi = np.ma.masked_where( Zi < 0, Zi)

fig = plt.figure(figsize=(6, 5))
cs = plt.pcolormesh(yi, z, Zi.T, vmin=0, vmax=70, cmap='pyart_NWSRef')

thin = 1

q = plt.quiver(yi, z, Ui_pj.T, Wi.T, scale = 400, color='black',)

plt.quiverkey(q, X=0.3, Y=1.1, U=10, label='Vertical Motion 10(m/s)', labelpos='E')

plt.xlim(0, 60)
plt.ylim(0, 20)

field = 'Reflectivity'
units = 'dBZ'

time_start = netCDF4.num2date(pyart_grid_1b.time['data'][0], pyart_grid_1b.time['units'])
time_text = time_start.strftime('%Y-%m-%dT%H:%M UTC')
plt.colorbar(cs, label= field + ' (' + units + ')')
plt.title('CSU-CHIVO '+ field + time_text + '\n \n \n')
plt.xlabel('North South distance from CHIVO (km)')
plt.ylabel('Altitude (km)')

