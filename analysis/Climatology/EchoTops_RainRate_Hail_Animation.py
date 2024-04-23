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

#%% Read data paths

chivo_dataPath = '/net/k2/storage/people/idariash/home/Field_Campaigns/Relampago/Analisis/EOL/data/20181214/chivo'
csapr_dataPath = '/net/k2/storage/people/idariash/home/Field_Campaigns/Relampago/Analisis/EOL/data/20181214/csapr'

chivo_files = [f for f in glob.glob(chivo_dataPath +"/**/*.nc", recursive=True)]
chivo_files.sort()

csapr_files = [f for f in glob.glob(csapr_dataPath +"/**/*.nc", recursive=True)]
csapr_files.sort()

for file_index in range(len(chivo_files)):
    
    #%% Combine RainRate
    
    #Read data
    
    chivo_1b_filename = chivo_files[file_index]
    
    csapr_1b_filename = csapr_files[file_index]
    
    # Create radar objects
    
    chivo_1b_radar = pyart.io.read(chivo_1b_filename)
    
    csapr_1b_radar = pyart.io.read(csapr_1b_filename)
    
    # Create grids
    
    chivo_1b_grid = utl.grid_rainRate_oneSweep(radar = chivo_1b_radar, origin = [-31.6342, -64.1686])
    
    csapr_1b_grid = utl.grid_rainRate_oneSweep(radar = csapr_1b_radar, origin = [-31.6342, -64.1686])
    
    # Read reflectivity from grids
    
    chivo_1b_ref = chivo_1b_grid.fields['corrected_reflectivity']['data']
    
    csapr_1b_ref = csapr_1b_grid.fields['corrected_reflectivity']['data']
       
    # Read RainRate
    
    chivo_rainRate = chivo_1b_grid.fields['RainRate']['data'][0]
    csapr_rainRate = csapr_1b_grid.fields['RainRate']['data'][0]
    
    # Masking RainRate
    chivo_rainRate = np.ma.masked_where(chivo_rainRate < 1, chivo_rainRate)
    #chivo_rainRate.filled(0)
    #chivo_rainRate.mask = False
    
    csapr_rainRate = np.ma.masked_where(csapr_rainRate < 1, csapr_rainRate)
    #csapr_rainRate.filled(0)
    #csapr_rainRate.mask = False
    
    RainRate_Radars = []
    RainRate_Radars = np.ma.array([chivo_rainRate, csapr_rainRate])#, rma_rainRate])
    RainRate_average = np.ma.mean(RainRate_Radars, 0)
    
    RainRate_combined = RainRate_average
    print('done with rainrate')
    
    #%% Joint EchoTop
    
    import scipy.ndimage as spyi
    
    #Data read in RainRate cell
    
    # Radar objects created in the RainRate cell 
    
    # Create grids
    
    chivo_grid_EchoTop = utl.grid_EchoTop(radar = chivo_1b_radar, origin = [-31.6342, -64.1686])
    
    csapr_grid_EchoTop = utl.grid_EchoTop(radar = csapr_1b_radar, origin = [-31.6342, -64.1686])
    
    # Read reflectivity from grids
    
    chivo_ref_EchoTop = chivo_grid_EchoTop.fields['corrected_reflectivity']['data']
    
    csapr_ref_EchoTop = csapr_grid_EchoTop.fields['corrected_reflectivity']['data']
    
    # Max. 18 dBZ reflectivity
    y, z, x = np.meshgrid(0.001*chivo_grid_EchoTop.y['data'],  
                          0.001*chivo_grid_EchoTop.z['data'],  0.001*chivo_grid_EchoTop.x['data'])
    
    # chivo_ref_EchoTop = spyi.gaussian_filter(chivo_ref_EchoTop, sigma = 4)
    # csapr_ref_EchoTop = spyi.gaussian_filter(csapr_ref_EchoTop, sigma = 4)
    # rma_ref_EchoTop = spyi.gaussian_filter(rma_ref_EchoTop, sigma = 4)
    
    # Combine the reflectivities
    chivo_ref_EchoTop = np.ma.masked_where(chivo_ref_EchoTop < 10, chivo_ref_EchoTop)
    csapr_ref_EchoTop = np.ma.masked_where(csapr_ref_EchoTop < 10, csapr_ref_EchoTop)
    
    Ref_Radars = []
    Ref_Radars = np.ma.array([chivo_ref_EchoTop, csapr_ref_EchoTop])#, rma_ref_EchoTop])
    Ref_average = np.ma.mean(Ref_Radars, 0)
    Ref_average = spyi.gaussian_filter(Ref_average, sigma = 4)
    
    z_EchoTop = np.ma.masked_where(Ref_average < 18, z)
    
    # z_chivo = np.ma.masked_where( chivo_ref_EchoTop < 18, z)
    # z_csapr = np.ma.masked_where( csapr_ref_EchoTop < 18, z)
    # z_rma = np.ma.masked_where( csapr_ref_EchoTop < 18, z)
    
    EchoTop_combined = np.amax(z_EchoTop, axis=0)
    EchoTop_combined = np.ma.masked_where( EchoTop_combined < 4, EchoTop_combined)
    
    print('done with echotops')
    
    #%% Combined hailmap
    
    # Data read in RainRate cell
    
    # Radar objects created in the RainRate cell 
    
    # Create grids
    
    chivo_grid_hail = utl.grid_hail_oneSweep(radar = chivo_1b_radar, origin = [-31.6342, -64.1686])
    csapr_grid_hail = utl.grid_hail_oneSweep(radar = csapr_1b_radar, origin = [-31.6342, -64.1686])
    
    # Read HydroClass from grids
    
    chivo_HydroClass = chivo_grid_hail.fields['HydroClass']['data'][0]
    csapr_HydroClass = csapr_grid_hail.fields['HydroClass']['data'][0]
    
    chivo_HydroClass.filled(0)
    chivo_HydroClass.mask = False
    
    csapr_HydroClass.filled(0)
    csapr_HydroClass.mask = False
    
    def hail_flag(HID):
        if HID < 10 or 11 < HID:
            return 0
        else:
            return 1
    
    for i in  range(len(chivo_HydroClass)):
        for j in range (len(chivo_HydroClass[i])):
            chivo_HydroClass[i][j] = hail_flag(chivo_HydroClass[i][j])
            csapr_HydroClass[i][j] = hail_flag(csapr_HydroClass[i][j])
            
    # Combine hail detection
    
    hail_detection = []       
    hail_detection = np.ma.array([chivo_HydroClass, csapr_HydroClass,])# rma_HydroClass])
    hail_detection = np.ma.amax(hail_detection, 0)
    
    print('done with hail')
    
    #%% Read the Sierras
    
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
    Sierras_terr = spyi.gaussian_filter(Sierras_terr, sigma = 2)
    #%%
    
    fig = plt.figure(figsize=(6, 5))
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
    
    EchoTop_combined = np.ma.masked_where( R > 130, EchoTop_combined)
    EchoTop_combined = np.ma.masked_where(EchoTop_combined < 11, EchoTop_combined)
    
    hail_detection = np.ma.masked_where( R > 130, hail_detection)
    hail_detection = np.ma.masked_where( hail_detection < 1, hail_detection)

    # cs = m.pcolormesh(chivo_lons, chivo_lats, RainRate_combined, vmin=1, vmax=100, latlon=True,
    #         cmap=pyart.graph.cm.LangRainbow12, norm=colors.LogNorm()) #  cmocean.cm.rain
    # m.colorbar(cs, label='Rain Rate (mm/h)')
    
    m.scatter(chivo_lons, chivo_lats, s = 1 ,c = hail_detection, marker = 'o',
              vmin=0, vmax=1.3, latlon=True, cmap=cmocean.cm.rain, alpha = 0.2)
    
    m.contour(chivo_lons, chivo_lats, EchoTop_combined, [12, 14, 18], 
                   latlon = True, linewidths = 1, colors = 'k',
                          linestyles='dashed',antialiased=True,  animated=True)
    
    
    time_start = netCDF4.num2date(chivo_1b_grid.time['data'][0], chivo_1b_grid.time['units'])
    time_text = time_start.strftime('%Y/%m/%d %H:%M')
    time_text = time_start.strftime('%Y/%m/%dT%H:%M')
    
    
    #plt.title('CSU-CHIVO | '+ '2018/11/10 to 2019/01/31' + ' \n Accumulated Rain Map')
    plt.title('CHIVO & CSAPR \n'+ time_text+' UTC')
    
    time_text = time_start.strftime('%Y%m%dT%H%M')
    fig_name = './fig/multiple_radars/echoTop_hail/radar_retrieval_' + time_text + 'png'
    plt.savefig(fig_name)
    del (m)
    gc.collect()
