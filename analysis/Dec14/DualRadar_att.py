#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 12 13:13:25 2021

@author: idariash
"""

# Dual radar attenuation retrieval

import sys
import dual_radar_plot as plot
import importlib
import pyart
import numpy as np
import matplotlib.pyplot as plt
import netCDF4

ind = 1

filename_chivo_1a = 'chivo_1a_20181214_0230.nc'
filename_chivo_1b = 'chivo_1b_test_20181214_0230.nc'

filename_csapr_1a = 'csapr_b1_20181214_0230.nc'
filename_csapr_1b = 'csapr_1b_test_20181214_0230.nc'

# Create grids
grid_chivo_1a = pyart.io.read_grid(filename_chivo_1a)
grid_chivo_1b = pyart.io.read_grid(filename_chivo_1b)

grid_csapr_1a = pyart.io.read_grid(filename_csapr_1a)
grid_csapr_1b = pyart.io.read_grid(filename_csapr_1b)

# Read reflectivity
ref_chivo_1a = grid_chivo_1a.fields['reflectivity']['data']
ref_chivo_1b = grid_chivo_1b.fields['corrected_reflectivity']['data']

ref_csapr_1a = grid_csapr_1a.fields['reflectivity']['data']
ref_csapr_1b = grid_csapr_1b.fields['corrected_reflectivity']['data']


# Mask low values
threshold_ref = 0
ref_chivo_1a = np.ma.masked_where( ref_chivo_1a < threshold_ref, ref_chivo_1a)
ref_chivo_1b = np.ma.masked_where( ref_chivo_1b < threshold_ref, ref_chivo_1b)

ref_csapr_1a = np.ma.masked_where( ref_csapr_1a < threshold_ref, ref_csapr_1a)
ref_csapr_1b = np.ma.masked_where( ref_csapr_1b < threshold_ref, ref_csapr_1b)

# Compute attenuation

att_chivo = ref_chivo_1b - ref_chivo_1a
att_csapr = ref_csapr_1b - ref_csapr_1a
min_att = np.minimum(att_chivo, att_csapr)
min_att = np.ma.masked_where( min_att > 3, min_att)


# Combined reflectivity 

combined_field_var = np.divide((np.multiply(att_csapr, ref_chivo_1b) + 
                                    np.multiply(att_chivo, ref_csapr_1b)), 
                                   att_chivo + att_csapr)
    
combined_field_var = np.ma.masked_where( combined_field_var < 18, combined_field_var)    

for layer in range(1,6):
    
    cappi = combined_field_var[layer]
    
    for row in cappi:
        for element in row:
            



fig = plt.figure(figsize=(6, 5))
x, y = np.meshgrid(0.001*grid_chivo_1a.x['data'], 0.001*grid_chivo_1a.y['data'])

field = 'Horizontal att.'
units = 'dB' 
 
cs = plt.pcolormesh(x, y, min_att[ind], vmin=0, vmax=15, cmap='pyart_NWSRef')



 
  
plt.contour(x, y, combined_field_var[ind], levels=[30], colors=['k', 'k', 'k'], antialised = True)

 #CHIVO mark
plt.plot(0, 0, 'ok', markersize=5)
plt.text(0, 0, ' C', fontsize=12)
# CSAPR mark
plt.plot(-54, -54, 'ok', markersize=5)
plt.text(-54, -54, ' Y', fontsize=12)
# # RMA mark
# plt.plot(-10, 20, 'ok', markersize=5)
# plt.text(-10, 20, ' R', fontsize=12)

plt.xlim(-90, 40)
plt.ylim(-90, 40)

time_start = netCDF4.num2date(grid_chivo_1a.time['data'][0], grid_chivo_1a.time['units'])
time_text = time_start.strftime('%Y-%m-%dT%H:%M UTC')
plt.colorbar(cs, label= 'Attenuation '+ ' (' + units + ')')
plt.title('CHIVO '+ field + ' Attenuation at '+ str(ind + 1) +' km \n'+ time_text)
plt.xlabel('Distance east of CHIVO (km)')
plt.ylabel('Distance north of CHIVO (km)')