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
import multiprocessing as mp
import gc
import os



#Read data

chivo_1a_filename = '/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/level_1a/2018/12/14/chivo.1a.20181214_023049.REL_PNL360A.nc'
chivo_1b_filename = '/net/k2/storage/people/idariash/home/Field_Campaigns/Relampago/Analisis/EOL/data/chivo_ppi/chivo.1b.2.20181214_023049.REL_PNL360A.nc'

csapr_b1_filename = '/net/k2/storage/projects/CSAPR/DOE_b1/corcsapr2cfrppiqcM1.b1.20181214.023003.nc'
csapr_1b_filename = '/net/k2/storage/people/idariash/home/Field_Campaigns/Relampago/Analisis/EOL/data/csapr/corcsapr2cfrppiqcM1.1b.20181214.023003.nc'

rma_cfradial_filename = '/net/k2/storage/people/idariash/home/Field_Campaigns/Relampago/Analisis/EOL/data/rma/original/cfrad.20181214_023044.000_to_20181214_023724.950_1_SUR.nc'
rma_1b_filename = '/net/k2/storage/people/idariash/home/Field_Campaigns/Relampago/Analisis/EOL/data/rma/corrected/cfrad.20181214_023044.000_to_20181214_023724.950_1_SUR.nc'

# Create radar objects

chivo_1a_radar = pyart.io.read(chivo_1a_filename)
chivo_1b_radar = pyart.io.read(chivo_1b_filename)

csapr_b1_radar = pyart.io.read(csapr_b1_filename)
csapr_1b_radar = pyart.io.read(csapr_1b_filename)

rma_cfradial_radar = pyart.io.read(rma_cfradial_filename)
rma_1b_radar = pyart.io.read(rma_1b_filename)

# Create grids

chivo_1a_grid = utl.grid_rainRate_oneSweep_original(radar = chivo_1a_radar, origin = [-31.6342, -64.1686])
chivo_1b_grid = utl.grid_rainRate_oneSweep(radar = chivo_1b_radar, origin = [-31.6342, -64.1686])

csapr_b1_grid = utl.grid_rainRate_oneSweep_original(radar = csapr_b1_radar, origin = [-31.6342, -64.1686])
csapr_1b_grid = utl.grid_rainRate_oneSweep(radar = csapr_1b_radar, origin = [-31.6342, -64.1686])

rma_cfradial_grid = utl.grid_rainRate_oneSweep_original(radar = rma_cfradial_radar, origin = [-31.6342, -64.1686])
rma_1b_grid = utl.grid_rainRate_oneSweep(radar = rma_1b_radar, origin = [-31.6342, -64.1686])


# Read reflectivity from grids

chivo_1a_ref = chivo_1a_grid.fields['reflectivity']['data']
chivo_1b_ref = chivo_1b_grid.fields['corrected_reflectivity']['data']

csapr_b1_ref = csapr_b1_grid.fields['reflectivity']['data']
csapr_1b_ref = csapr_1b_grid.fields['corrected_reflectivity']['data']

rma_cfradial_ref = rma_cfradial_grid.fields['reflectivity']['data']
rma_1b_ref = rma_1b_grid.fields['corrected_reflectivity']['data']

# Compute attenuation
chivo_att = chivo_1b_ref[0] - chivo_1a_ref[0]
csapr_att = csapr_1b_ref[0] - csapr_b1_ref[0]
rma_att = rma_1b_ref[0] - rma_cfradial_ref[0]

# Read RainRate

chivo_rainRate = chivo_1b_grid.fields['RainRate']['data'][0]
csapr_rainRate = csapr_1b_grid.fields['RainRate']['data'][0]
rma_rainRate = rma_1b_grid.fields['RainRate']['data'][0]/1.5

# Masking RainRate
chivo_rainRate = np.ma.masked_where(chivo_rainRate < 1, chivo_rainRate)
#chivo_rainRate.filled(0)
#chivo_rainRate.mask = False

csapr_rainRate = np.ma.masked_where(csapr_rainRate < 1, csapr_rainRate)
#csapr_rainRate.filled(0)
#csapr_rainRate.mask = False

rma_rainRate = np.ma.masked_where(rma_rainRate < 1, rma_rainRate)
#rma_rainRate.filled(0)
#rma_rainRate.mask = False

# Masking attenuation

#chivo_att = np.ma.masked_where(chivo_att < -2, chivo_att)
#chivo_att.filled(0)
chivo_att.mask = False

#csapr_att = np.ma.masked_where(csapr_att < -2, csapr_att)
#csapr_att.filled(0)
csapr_att.mask = False

#rma_att = np.ma.masked_where(rma_att < -2, rma_att)
#rma_att.filled(0)
rma_att.mask = False

array_test = np.ma.array([chivo_rainRate, csapr_rainRate, rma_rainRate])
test_average = np.ma.mean(array_test, 0)






#         total_rain_file = (grid.fields['RainRate']['data'].copy())/1.5
#         total_rain_file = np.ma.masked_where(total_rain_file < 0, total_rain_file)
#         total_rain_file.filled(0)
#         total_rain_file.mask = False
# 

# Combine rain rate
RainRate_combined = np.divide(np.multiply(np.multiply(csapr_att, rma_att), chivo_rainRate) +
                            np.multiply(np.multiply(chivo_att, rma_att), csapr_rainRate) + 
                            np.multiply(np.multiply(chivo_att, csapr_att), rma_rainRate),
                            np.multiply(csapr_att, rma_att) + np.multiply(chivo_att, rma_att) +
                            np.multiply(chivo_att, csapr_att))

RainRate_combined = test_average


# Plot

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
 

lon_0 = chivo_1b_grid.origin_longitude['data']
lat_0 = chivo_1b_grid.origin_latitude['data']

delta = 6#deg
lats_mask = np.absolute(lats - lat_0) < delta

lons_mask = np.absolute(lons - lon_0) < delta

Sierras_terr = terr[lats_mask, :]
Sierras_terr = Sierras_terr[:, lons_mask]
Sierras_terr = np.ma.masked_where( Sierras_terr < 0.5, Sierras_terr)

Sierras_lons = lons[lons_mask]
Sierras_lats = lats[lats_mask]

Sierras_lons_grid, Sierras_lats_grid = np.meshgrid(Sierras_lons, Sierras_lats)

#%%

Sierras_terr = spyi.gaussian_filter(Sierras_terr, sigma = 2)

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



lon_0 = chivo_1b_grid.origin_longitude['data']
lat_0 = chivo_1b_grid.origin_latitude['data']
chivo_lons, chivo_lats = chivo_1b_grid.get_point_longitude_latitude()

x, y = np.meshgrid(0.001*chivo_1b_grid.x['data'],  0.001*chivo_1b_grid.y['data'])
R = np.sqrt(x ** 2 + y ** 2)

# RainRate = np.ma.masked_where( R > 100, RainRate)
# RainRate = np.ma.masked_where( RainRate < 1, RainRate)
# HydroClass = np.ma.masked_where(HydroClass < 10, HydroClass)
# HydroClass = np.ma.masked_where(HydroClass > 11, HydroClass)
# HydroClass = np.ma.masked_where(R > 80, HydroClass)
# HydroClass = np.ma.masked_where(max_ref <  40, HydroClass)

RainRate_combined = np.ma.masked_where(RainRate_combined <  2, RainRate_combined)
RainRate_combined = np.ma.masked_where(R > 130, RainRate_combined)

cs = m.pcolormesh(chivo_lons, chivo_lats, RainRate_combined, vmin=1, vmax=100, latlon=True,
        cmap=pyart.graph.cm.LangRainbow12, norm=colors.LogNorm()) #  cmocean.cm.rain

#m.colorbar(cs, label='Total Rain (mm)')
m.colorbar(cs, label='Rain Rate (mm/h)')

# time_start = netCDF4.num2date(grid.time['data'][0], grid.time['units'])
# time_text = time_start.strftime('%Y-%m-%dT%H:%M')


#plt.title('CSU-CHIVO | '+ '2018/11/10 to 2019/01/31' + ' \n Accumulated Rain Map')
plt.title('Combined Rain Rate | '+ '2018/12/14 02:30 UTC')


# #files = [f for f in glob.glob(DataPath +"/**/*20181214*P*.nc", recursive=True)]
# #files = [f for f in glob.glob(DataPath +"/**/*20181214*HYDRO*.nc", recursive=True)]
# files = [f for f in glob.glob(DataPath +"/**/*20181214*.nc", recursive=True)]
# files.sort()
# radar = pyart.io.read(files[0])
# grid = utl.grid_rainRate_oneSweep(radar = radar, origin = [-31.6342, -64.1686])
# # accumulated_total_rain  = (grid.fields['RainRate']['data'].copy())*0 # every 10 min per file
# # accumulated_total_rain = np.ma.masked_where(accumulated_total_rain < 0, accumulated_total_rain)
# # accumulated_total_rain.filled(0)
# # accumulated_total_rain.mask = False
# # accumulated_total_rain = np.amax(accumulated_total_rain, axis=0)

# def get_rain(file):
#     try:
#         filename = os.path.basename(file)
#         # if '20181214' in file or '2018111' in file or '20190125' in file:
#         #     print('not processing ' + filename)
#         #     return 0
#         radar = pyart.io.read(file)
#         grid = utl.grid_rainRate_oneSweep(radar = radar, origin = [-31.6342, -64.1686])
#         total_rain_file = (grid.fields['RainRate']['data'].copy())/1.5
#         total_rain_file = np.ma.masked_where(total_rain_file < 0, total_rain_file)
#         total_rain_file.filled(0)
#         total_rain_file.mask = False
#         total_rain_file = np.ma.amax(total_rain_file, axis=0)
#         del(radar)
#         del(grid)
#         gc.collect()        
#         return total_rain_file
#         print(filename)
#     except:
#         print('error ' + filename)

# pool = mp.Pool(40)
# results = pool.map_async(get_rain, files).get()
# pool.close()
# pool.join()

# accumulated_total_rain = sum(results)
# #total = np.sum((accumulated_total_rain, 0)
        
        
# #%%

# import matplotlib.pyplot as plt
# import matplotlib.colors as colors
# from netCDF4 import Dataset
# import netCDF4
# import scipy.ndimage as spyi
# import cmocean
# from mpl_toolkits.basemap import Basemap

# # Read topography

# topo_filename = '/net/k2/storage/people/idariash/home/CSU/RELAMPAGO/analysis/Mountain/terr.nc'
# topography = Dataset(topo_filename)

# lons = topography.variables['X'][:]
# lats = topography.variables['Y'][:]
# terr = 0.001*topography.variables['topo'][:]


# height_toCut = 2
# reflectivity_toCut = 5
# range_toCut = 200
 

# lon_0 = grid.origin_longitude['data']
# lat_0 = grid.origin_latitude['data']

# delta = 6#deg
# lats_mask = np.absolute(lats - lat_0) < delta

# lons_mask = np.absolute(lons - lon_0) < delta

# Sierras_terr = terr[lats_mask, :]
# Sierras_terr = Sierras_terr[:, lons_mask]
# Sierras_terr = np.ma.masked_where( Sierras_terr < 0.5, Sierras_terr)

# Sierras_lons = lons[lons_mask]
# Sierras_lats = lats[lats_mask]

# Sierras_lons_grid, Sierras_lats_grid = np.meshgrid(Sierras_lons, Sierras_lats)

# #%%

# Sierras_terr = spyi.gaussian_filter(Sierras_terr, sigma = 2)

# m = Basemap(projection='lcc', resolution=None,
#             width=300E3, height=300E3, 
#             lat_0=lat_0, lon_0=lon_0)

# #m.etopo(scale=2, alpha=0.5)

# parallels = np.arange(-40,-20,1)
# # labels = [left,right,top,bottom]
# m.drawparallels(parallels,labels=[True,False,False, False])

# meridians = np.arange(-70,-60.,1)
# m.drawmeridians(meridians,labels=[True,False,False,True])

# x_Cordoba, y_Cordoba = m(-64.19, -31.42)
# x_CHIVO, y_CHIVO = m(-64.1686, -31.6342)
# x_ARM, y_ARM = m(-64.7284, -32.1264)
# x_Rio4to, y_Rio4to = m(-64.3493, -33.1232)

# plt.plot(x_CHIVO, y_CHIVO, 'ok', markersize=5)
# plt.text(x_CHIVO, y_CHIVO, ' C', fontsize=12);

# plt.plot(x_ARM, y_ARM, 'ok', markersize=5)
# plt.text(x_ARM, y_ARM, ' Y', fontsize=12);

# plt.plot(x_Cordoba, y_Cordoba, 'ok', markersize=5)
# plt.text(x_Cordoba, y_Cordoba, ' R', fontsize=12);

# levels_terr = [1, 2]

# m.contour(Sierras_lons_grid, Sierras_lats_grid, Sierras_terr, levels_terr, 
#                latlon = True, linewidths = 1, colors = 'k',
#                       linestyles='solid', antialiased=True, animated=True)



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

# accumulated_total_rain = np.ma.masked_where(accumulated_total_rain <  2, accumulated_total_rain)
# accumulated_total_rain = np.ma.masked_where(R > 130, accumulated_total_rain)

# cs = m.pcolormesh(chivo_lons, chivo_lats, accumulated_total_rain, vmin=1, vmax=100, latlon=True,
#         cmap=pyart.graph.cm.LangRainbow12, norm=colors.LogNorm()) #  cmocean.cm.rain

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
