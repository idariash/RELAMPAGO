#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 16:31:22 2022

@author: idariash
"""

# This code compine the rain rate from many radars

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


def get_rain(file):
    try:
        filename = os.path.basename(file)
        radar = pyart.io.read(file)
        time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
        time_UTC = time_start.strftime('%Y-%m-%dT%H:%M:%S') 
        if 'HYDRO' in file or 'RHI' in file:
            print('not processing ' + filename)
            return (time_UTC, 0)
        if 'RMA' in file:
            #file_size = os.path.getsize(file)/1E6
            nsweeps = radar.nsweeps
            if nsweeps > 4:#file_size > 18:
                return (time_UTC, 0)
        
        grid = utl.grid_rainRate_oneSweep(radar = radar, origin = [-31.6342, -64.1686])
        total_rain_file = (grid.fields['RainRate']['data'].copy())
        total_rain_file = np.ma.masked_where(total_rain_file < 0, total_rain_file)
        total_rain_file.filled(0)
        total_rain_file.mask = False
        total_rain_file = np.ma.amax(total_rain_file, axis=0)
        del(radar)
        del(grid)
        gc.collect()        
        return (time_UTC, total_rain_file)
        print(filename)
    except:
        print('error ' + filename)

# For CHIVO
DataPath = '/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/level_1b.2'
files = [f for f in glob.glob(DataPath +"/**/*20181214*.nc", recursive=True)]
files.sort()

pool = mp.Pool(40)
results = pool.map_async(get_rain, files).get()
pool.close()
pool.join()


chivo_rain_per_hour = np.zeros((24,801,801))
chivo_number_of_files_per_hour = np.zeros((24,1))

for i in range(0, len(results)):
    if len(str(results[i][1])) == 1:
        continue
    result_datetime = datetime.strptime(results[i][0], '%Y-%m-%dT%H:%M:%S')
    hour = result_datetime.hour
    chivo_rain_per_hour[hour] = chivo_rain_per_hour[hour] + results[i][1]
    chivo_number_of_files_per_hour[hour] = chivo_number_of_files_per_hour[hour] + 1    

# For CSAPR
DataPath = '/net/k2/storage/projects/CSAPR/level_1b.2'
files = [f for f in glob.glob(DataPath +"/**/*20181214*.nc", recursive=True)]
files.sort()

pool = mp.Pool(40)
results = pool.map_async(get_rain, files).get()
pool.close()
pool.join()


csapr_rain_per_hour = np.zeros((24,801,801))
csapr_number_of_files_per_hour = np.zeros((24,1))

for i in range(0, len(results)):
    if len(str(results[i][1])) == 1:
        continue
    result_datetime = datetime.strptime(results[i][0], '%Y-%m-%dT%H:%M:%S')
    hour = result_datetime.hour
    csapr_rain_per_hour[hour] = csapr_rain_per_hour[hour] + results[i][1]
    csapr_number_of_files_per_hour[hour] = csapr_number_of_files_per_hour[hour] + 1
    
    
# For RMA
DataPath = '/net/k2/storage/projects/RELAMPAGO/RMA1/quality_controlled/level_1b'
files = [f for f in glob.glob(DataPath +"/**/*20181214*.nc", recursive=True)]
files.sort()

pool = mp.Pool(40)
results = pool.map_async(get_rain, files).get()
pool.close()
pool.join()


rma_rain_per_hour = np.zeros((24,801,801))
rma_number_of_files_per_hour = np.zeros((24,1))

for i in range(0, len(results)):
    if len(str(results[i][1])) == 1:
        continue
    result_datetime = datetime.strptime(results[i][0], '%Y-%m-%dT%H:%M:%S')
    hour = result_datetime.hour
    rma_rain_per_hour[hour] = rma_rain_per_hour[hour] + results[i][1]
    rma_number_of_files_per_hour[hour] = rma_number_of_files_per_hour[hour] + 1
    
#%% Combine rain rate

# file_long = '/net/k2/storage/projects/RELAMPAGO/RMA1/quality_controlled/level_1b/cfrad.20181213_235847.000_to_20181214_000527.950_1_SUR.nc'
# file_short = '/net/k2/storage/projects/RELAMPAGO/RMA1/quality_controlled/level_1b/cfrad.20181214_000539.000_to_20181214_000701.922_1_SUR.nc'
# radar_long = pyart.io.read(file_long)
# radar_short = pyart.io.read(file_short)

radar = pyart.io.read(files[0])
radar_2 = pyart.io.read(files[1])
grid = utl.grid_rainRate_oneSweep(radar = radar, origin = [-31.6342, -64.1686])

#%%
Xo_chivo, Yo_chivo = 0,0
projparams = grid.get_projparams()
Xo_csapr, Yo_csapr = pyart.core.geographic_to_cartesian(-64.7284, -32.1264, projparams)
Xo_rma, Yo_rma = pyart.core.geographic_to_cartesian(-64.19, -31.42, projparams)

x, y = np.meshgrid(grid.x['data'],  grid.y['data'])
R_chivo = np.sqrt((x - Xo_chivo)** 2 + (y - Yo_chivo)** 2)
R_csapr = np.sqrt((x - Xo_csapr)** 2 + (y - Yo_csapr)** 2)
R_rma = np.sqrt((x - Xo_rma)** 2 + (y - Yo_rma)** 2)

chivo_csapr_distance = np.sqrt((Xo_chivo - Xo_csapr)**2+(Yo_chivo - Yo_csapr)**2)
chivo_rma_distance = np.sqrt((Xo_chivo - Xo_rma)**2+(Yo_chivo - Yo_rma)**2)
csapr_rma_distance = np.sqrt((Xo_csapr - Xo_rma)**2+(Yo_csapr - Yo_rma)**2)

chivo_azi = 90 - 180*np.arctan2(y,x)/np.pi
rma_azi = 90 - 180*np.arctan2(y - Yo_rma, x - Xo_rma)/np.pi

rain_per_hour = np.zeros((24,801,801))

for hour in range(0, 24):

    chivo_rain = chivo_rain_per_hour[hour]
    csapr_rain = csapr_rain_per_hour[hour]
    rma_rain = rma_rain_per_hour[hour]
    
    # Masking
    chivo_rain = np.ma.masked_where( chivo_rain < 1, chivo_rain)
    csapr_rain = np.ma.masked_where( csapr_rain < 1, csapr_rain)
    rma_rain = np.ma.masked_where( rma_rain < 1, rma_rain)
    
    chivo_rain = np.ma.masked_where( R_chivo > 100E3, chivo_rain)
    csapr_rain = np.ma.masked_where( R_csapr > 100E3, csapr_rain)
    rma_rain = np.ma.masked_where( R_rma > 100E3, rma_rain)
    
    chivo_rain = np.ma.masked_where(abs(chivo_azi - 90) < 45, chivo_rain)
    chivo_rain = np.ma.masked_where(abs(chivo_azi - -32) < 5, chivo_rain)
    
    rma_rain = np.ma.masked_where(abs(rma_azi - -20) < 2, rma_rain)
    rma_rain = np.ma.masked_where(abs(rma_azi - 247) < 2, rma_rain)
    
    chivo_rain = chivo_rain/chivo_number_of_files_per_hour[hour]
    csapr_rain = csapr_rain/(1.5*csapr_number_of_files_per_hour[hour])
    rma_rain = rma_rain/(2*rma_number_of_files_per_hour[hour])
    
    chivo_rain = np.ma.masked_where( chivo_rain < 5, chivo_rain)
    csapr_rain = np.ma.masked_where( csapr_rain < 5, csapr_rain)
    rma_rain = np.ma.masked_where( rma_rain < 5, rma_rain)
    combined_rain = [chivo_rain, csapr_rain, rma_rain]
    rain = np.ma.mean(combined_rain, 0)
    rain_per_hour[hour] = rain
    
    if hour == 5:
        rain_std = np.ma.std(combined_rain, 0)
        rain_plot = rma_rain
    
rain = np.ma.sum(rain_per_hour, 0)

#%%
rain_chivo = np.zeros((24,801,801))
for hour in range(0, 24):
    chivo_rain = csapr_rain_per_hour[hour]
    chivo_rain = np.ma.masked_where( chivo_rain < 1, chivo_rain)
    chivo_rain = chivo_rain/csapr_number_of_files_per_hour[hour]
    rain_chivo[hour] = chivo_rain
    
chivo_rain = np.ma.sum(rain_chivo, 0)
    
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

hour = 2
accumulated_total_rain = rma_rain_per_hour[hour]/(rma_number_of_files_per_hour[hour]*2)

accumulated_total_rain = np.ma.masked_where(accumulated_total_rain <  2, accumulated_total_rain)
accumulated_total_rain = np.ma.masked_where(R > 130, accumulated_total_rain)
accumulated_total_rain = rain #rain_plot#rain_std#rain

levels = np.arange(10, 110, 10)

cs = m.contourf(chivo_lons, chivo_lats, accumulated_total_rain, levels,
                  cmap= pyart.graph.cm.LangRainbow12, latlon=True)

# cs = m.pcolormesh(chivo_lons, chivo_lats, accumulated_total_rain, vmin=10, vmax=100, latlon=True,
#         cmap=pyart.graph.cm.LangRainbow12,) # norm=colors.LogNorm()) #  cmocean.cm.rain #pyart.graph.cm.LangRainbow12

#m.colorbar(cs, label='Total Rain (mm)')
m.colorbar(cs, label='Rain Rate (mm/h)')

# time_start = netCDF4.num2date(grid.time['data'][0], grid.time['units'])
# time_text = time_start.strftime('%Y-%m-%dT%H:%M')

radar_name = 'CSU-CHIVO'

if 'csapr' in DataPath:
    radar_name = 'CSAPR-2'
    
if 'rma' in DataPath:
    radar_name = 'RMA-1'

#plt.title('CSU-CHIVO | '+ '2018/11/10 to 2019/01/31' + ' \n Accumulated Rain Map')
plt.title(radar_name + ' | '+ '2018/12/14 02:30 UTC')


#accumulated_total_rain = sum(results)

#%%
# import chivo_plot as chivo_plot

# file = '/net/k2/storage/projects/RELAMPAGO/RMA1/quality_controlled/level_1b/cfrad.20181214_032757.000_to_20181214_032919.922_1_SUR.nc'

# chivo_plot.ppi_drops(file, 0)

arr0= accumulated_total_rain
arr1 = chivo_lons
arr2 = chivo_lats

explort_file = '/net/nasstore/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Climatology/RadarsRain_Dec14.npz'
np.savez(explort_file, radarsRain = arr0, chivo_lons = arr1, chivo_lats = arr2)

