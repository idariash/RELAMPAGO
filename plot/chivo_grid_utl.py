#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 14:47:53 2020

@author: idariash
"""

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

clevs = [0,1,5,10,20,30,40,50,75,100,125,150,200,250,300,350,400,500]
#[0,1,2.5,5,10,20,30,40,50,60,80,100,140,160,180,200,300,400]
cmap_data = [(1.000000, 1.000000, 1.000000),
(0.498039, 1.000000, 0.000000),
(0.000000, 0.803922, 0.000000),
(0.000000, 0.545098, 0.000000),
(0.062745, 0.305882, 0.545098),
(0.117647, 0.564706, 1.000000),
(0.000000, 0.698039, 0.933333),
(0.000000, 0.933333, 0.933333),
(0.537255, 0.407843, 0.803922),
(0.568627, 0.172549, 0.933333),
(0.545098, 0.000000, 0.545098),
(0.545098, 0.000000, 0.000000),
(0.803922, 0.000000, 0.000000),
(0.933333, 0.250980, 0.000000),
(1.000000, 0.498039, 0.000000),
(0.803922, 0.521569, 0.000000),
(1.000000, 0.843137, 0.000000)]
cmap_precip = mcolors.ListedColormap(cmap_data, 'precip_rss')
norm_precip = mcolors.BoundaryNorm(clevs, cmap_precip.N)

def plot_grid_echoTop(table_filename, radar_grid_filename, topo_filename):
    
    topography = Dataset(topo_filename)

    lons = topography.variables['X'][:]
    lats = topography.variables['Y'][:]
    terr = 0.001*topography.variables['topo'][:]
    
    
    height_toCut = 2
    reflectivity_toCut = 5
    range_toCut = 200
    radar_grid = pyart.io.read_grid(radar_grid_filename)    
    
    lon_0 = radar_grid.origin_longitude['data']
    lat_0 = radar_grid.origin_latitude['data']
    
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
    
    # Create scatter plot
    
    df = pd.read_excel(table_filename)
    X = df.x_echo_top.to_numpy()*1000
    Y = df.y_echo_top.to_numpy()*1000
    echoTops = df.echo_top.to_numpy()
    echoTops = np.ma.masked_where( echoTops < 10, echoTops)
    [lons_echoTop, lats_echoTops] = pyart.core.cartesian_to_geographic(X, Y, radar_grid.get_projparams())
    
    cs = m.scatter(lons_echoTop, lats_echoTops, s = 10, c = echoTops, vmin=4, vmax=20, latlon=True,
                   cmap=cmocean.cm.rain)#, norm=norm_precip)

        
    

    
    # accumulated_total_rain = chivo_grid.fields['accumulated_total_rain']['data']
    
    # x, y, z = np.meshgrid(0.001*chivo_grid.x['data'], 0.001*chivo_grid.y['data'],
    #                      0.001*chivo_grid.z['data'])
    
    # x = np.swapaxes(x,1,2)
    # x = np.swapaxes(x,0,1)
    
    # y = np.swapaxes(y,1,2)
    # y = np.swapaxes(y,0,1)
    
    # z = np.swapaxes(z,1,2)
    # z = np.swapaxes(z,0,1)
    
    # R = np.sqrt(x ** 2 + y ** 2)
    
    # accumulated_total_rain = np.ma.masked_where( z < 1.9, accumulated_total_rain)
    # accumulated_total_rain = np.ma.masked_where( z > 2.1, accumulated_total_rain)
    # accumulated_total_rain = np.ma.masked_where( R > range_toCut, accumulated_total_rain)
    
    # accumulated_total_rain_max = np.ma.mean(accumulated_total_rain, axis=0)
    # #accumulated_total_rain_max = np.ma.masked_where( accumulated_total_rain_max < reflectivity_toCut, accumulated_total_rain_max)
    
    # lon_0 = chivo_grid.origin_longitude['data']
    # lat_0 = chivo_grid.origin_latitude['data']
    # chivo_lons, chivo_lats = chivo_grid.get_point_longitude_latitude()
    
   
    # cs = m.pcolormesh(chivo_lons, chivo_lats, accumulated_total_rain_max, vmin=0, vmax=150, latlon=True,
    #              cmap=cmocean.cm.rain)#, norm=norm_precip)
    
    m.colorbar(cs, label='Height (km)')
    
    time_start = netCDF4.num2date(radar_grid.time['data'][0], radar_grid.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M UTC')
    
    
    plt.title('CSU-CHIVO Echo Tops')
    
    return fig




def plot_grid_AccumulateTotalPrecip(grid_filename, topo_filename):
    
    topography = Dataset(topo_filename)

    lons = topography.variables['X'][:]
    lats = topography.variables['Y'][:]
    terr = 0.001*topography.variables['topo'][:]
    
    
    height_toCut = 2
    reflectivity_toCut = 5
    range_toCut = 200
    chivo_grid = pyart.io.read_grid(grid_filename)
    
    lon_0 = chivo_grid.origin_longitude['data']
    lat_0 = chivo_grid.origin_latitude['data']
    
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
                                           
    plt.plot(x_ARM, y_ARM, 'ok', markersize=5)
    plt.text(x_ARM, y_ARM, ' Y', fontsize=12);
                                           
    levels_terr = [1, 2]
    
    m.contour(Sierras_lons_grid, Sierras_lats_grid, Sierras_terr, levels_terr, 
                   latlon = True, linewidths = 1, colors = 'k',
                          linestyles='solid', antialiased=True, animated=True)
    
    accumulated_total_rain = chivo_grid.fields['accumulated_total_rain']['data']
    
    x, y, z = np.meshgrid(0.001*chivo_grid.x['data'], 0.001*chivo_grid.y['data'],
                         0.001*chivo_grid.z['data'])
    
    x = np.swapaxes(x,1,2)
    x = np.swapaxes(x,0,1)
    
    y = np.swapaxes(y,1,2)
    y = np.swapaxes(y,0,1)
    
    z = np.swapaxes(z,1,2)
    z = np.swapaxes(z,0,1)
    
    R = np.sqrt(x ** 2 + y ** 2)
    
    accumulated_total_rain = np.ma.masked_where( z < 1.9, accumulated_total_rain)
    accumulated_total_rain = np.ma.masked_where( z > 2.1, accumulated_total_rain)
    accumulated_total_rain = np.ma.masked_where( R > range_toCut, accumulated_total_rain)
    
    accumulated_total_rain_max = np.ma.mean(accumulated_total_rain, axis=0)
    #accumulated_total_rain_max = np.ma.masked_where( accumulated_total_rain_max < reflectivity_toCut, accumulated_total_rain_max)
    
    lon_0 = chivo_grid.origin_longitude['data']
    lat_0 = chivo_grid.origin_latitude['data']
    chivo_lons, chivo_lats = chivo_grid.get_point_longitude_latitude()
    
   
    cs = m.pcolormesh(chivo_lons, chivo_lats, accumulated_total_rain_max, vmin=0, vmax=150, latlon=True,
                 cmap=cmocean.cm.rain)#, norm=norm_precip)
    
    m.colorbar(cs, label='precipitation (mm)')
    
    time_start = netCDF4.num2date(chivo_grid.time['data'][0], chivo_grid.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M UTC')
    
    
    plt.title('CSU-CHIVO | '+ time_text+ ' \n Event total Precip.')
    
    return fig

def plot_grid_RainRate(grid_filename, topo_filename):
    
    topography = Dataset(topo_filename)

    lons = topography.variables['X'][:]
    lats = topography.variables['Y'][:]
    terr = 0.001*topography.variables['topo'][:]
    
    
    height_toCut = 2
    reflectivity_toCut = 5
    range_toCut = 200
    chivo_grid = pyart.io.read_grid(grid_filename)
    
    lon_0 = chivo_grid.origin_longitude['data']
    lat_0 = chivo_grid.origin_latitude['data']
    
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
                                           
    plt.plot(x_ARM, y_ARM, 'ok', markersize=5)
    plt.text(x_ARM, y_ARM, ' Y', fontsize=12);
                                           
    levels_terr = [0.5, 1.5]
    
    m.contour(Sierras_lons_grid, Sierras_lats_grid, Sierras_terr, levels_terr, 
                   latlon = True, linewidths = 1, colors = 'k',
                          linestyles='solid', antialiased=True, animated=True)
    
    RainRate = chivo_grid.fields['RainRate']['data']
    
    x, y, z = np.meshgrid(0.001*chivo_grid.x['data'], 0.001*chivo_grid.y['data'],
                         0.001*chivo_grid.z['data'])
    
    x = np.swapaxes(x,1,2)
    x = np.swapaxes(x,0,1)
    
    y = np.swapaxes(y,1,2)
    y = np.swapaxes(y,0,1)
    
    z = np.swapaxes(z,1,2)
    z = np.swapaxes(z,0,1)
    
    R = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    
    RainRate = np.ma.masked_where( z < 1.5, RainRate)
    RainRate = np.ma.masked_where( z > 2.5, RainRate)
    RainRate = np.ma.masked_where( R > range_toCut, RainRate)
    
    RainRate_max = np.mean(RainRate, axis=0)
    #RainRate_max = np.ma.masked_where( RainRate_max < reflectivity_toCut, RainRate_max)
    
    lon_0 = chivo_grid.origin_longitude['data']
    lat_0 = chivo_grid.origin_latitude['data']
    chivo_lons, chivo_lats = chivo_grid.get_point_longitude_latitude()
    
   
    cs = m.pcolormesh(chivo_lons, chivo_lats, RainRate_max, vmin=0, vmax=100, latlon=True,
                 cmap=cmap_precip, norm=norm_precip)
    
    m.colorbar(cs, label='Rain Rate (mm/h)')
    
    time_start = netCDF4.num2date(chivo_grid.time['data'][0], chivo_grid.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M UTC')
    
    
    plt.title('CSU-CHIVO | '+ time_text+ ' \n Rain Rate')
    
    return fig

def plot_grid_maxRef(grid_filename, topo_filename):
    
    topography = Dataset(topo_filename)

    lons = topography.variables['X'][:]
    lats = topography.variables['Y'][:]
    terr = 0.001*topography.variables['topo'][:]
    
    
    height_toCut = 2
    reflectivity_toCut = 5
    range_toCut = 200
    chivo_grid = pyart.io.read_grid(grid_filename)
    
    lon_0 = chivo_grid.origin_longitude['data']
    lat_0 = chivo_grid.origin_latitude['data']
    
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
    
    
    
    fig = plt.figure(figsize=(6, 5))
    
    
    
    m = Basemap(projection='lcc', resolution=None,
                width=300E3, height=300E3, 
                lat_0=lat_0, lon_0=lon_0)
    
    m.etopo(scale=2, alpha=0.5)
    
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
                                           
    plt.plot(x_ARM, y_ARM, 'ok', markersize=5)
    plt.text(x_ARM, y_ARM, ' Y', fontsize=12);
                                           
    levels_terr = [0.5, 1.5]
    
    m.contour(Sierras_lons_grid, Sierras_lats_grid, Sierras_terr, levels_terr, 
                   latlon = True, linewidths = 1, colors = 'k',
                          linestyles='solid', antialiased=True, animated=True)
    
    Z = chivo_grid.fields['corrected_reflectivity']['data']
    
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
    
    Z_max = np.amax(Z, axis=0)
    Z_max = np.ma.masked_where( Z_max < reflectivity_toCut, Z_max)
    
    lon_0 = chivo_grid.origin_longitude['data']
    lat_0 = chivo_grid.origin_latitude['data']
    chivo_lons, chivo_lats = chivo_grid.get_point_longitude_latitude()
    
    
    cs = m.pcolormesh(chivo_lons, chivo_lats, Z_max, vmin=0, vmax=70, latlon=True,
                 cmap='pyart_NWSRef')
    
    m.colorbar(cs, label='Reflectivity (dBZ)')
    
    time_start = netCDF4.num2date(chivo_grid.time['data'][0], chivo_grid.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M UTC')
    
    
    plt.title('CSU-CHIVO | '+ time_text+ ' \n Max. Reflectivity')
    
    return fig


def grid_radar(radar, grid_shape=(5, 401, 401), xlim=(-100000, 100000), 
               ylim=(-100000, 100000), zlim=(500, 4500),
               fields=['corrected_reflectivity', 'corrected_differential_reflectivity', 
                       'RainRate', 'HydroClass'], origin=None):
  
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('corrected_reflectivity')
    gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 10, 100)
    gatefilter.exclude_outside('normalized_coherent_power', 0.3, 1)
    gatefilter.exclude_outside('corrected_copol_correlation_ratio', 0.9, 1)

    
    if origin is None:
        origin = (radar.latitude['data'][0],
                  radar.longitude['data'][0])
        
    grid = pyart.map.grid_from_radars(
        radar, grid_shape=grid_shape,
        grid_limits=(zlim, ylim, xlim),
        grid_origin=origin, fields=fields,gatefilter = gatefilter,
        gridding_algo='map_gates_to_grid', grid_origin_alt=0.0)#, weighting_function = 'Nearest')

    
    return grid

def grid_hail(radar, grid_shape=(8, 801, 801), xlim=(-150000, 150000), 
               ylim=(-100000, 100000), zlim=(0, 2000),
               fields=['corrected_reflectivity', 'corrected_differential_reflectivity', 
                       'RainRate', 'HydroClass'], origin=None):
  
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('corrected_reflectivity')
    gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 0, 100)
    #gatefilter.exclude_outside('normalized_coherent_power', 0.3, 1)
    #gatefilter.exclude_outside('corrected_copol_correlation_ratio', 0.9, 1)

    
    if origin is None:
        origin = (radar.latitude['data'][0],
                  radar.longitude['data'][0])
        
    grid = pyart.map.grid_from_radars(
        radar, grid_shape=grid_shape,
        grid_limits=(zlim, ylim, xlim),
        grid_origin=origin, fields=fields,gatefilter = gatefilter,
        gridding_algo='map_gates_to_grid', grid_origin_alt=0, 
        min_radius = 100, nb = 1, bsp = 1, weighting_function = 'Nearest')

    
    return grid

def grid_rainRate(radar, grid_shape=(8, 801, 801), xlim=(-150000, 150000), 
               ylim=(-100000, 100000), zlim=(0, 2000),
               fields=['corrected_reflectivity', 'corrected_differential_reflectivity', 
                       'RainRate', 'HydroClass'], origin=None):
  
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('corrected_reflectivity')
    gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 10, 100)
    gatefilter.exclude_outside('normalized_coherent_power', 0.3, 1)
    gatefilter.exclude_outside('corrected_copol_correlation_ratio', 0.9, 1)

    
    if origin is None:
        origin = (radar.latitude['data'][0],
                  radar.longitude['data'][0])
        
    grid = pyart.map.grid_from_radars(
        radar, grid_shape=grid_shape,
        grid_limits=(zlim, ylim, xlim),
        grid_origin=origin, fields=fields,gatefilters = gatefilter,
        gridding_algo='map_gates_to_grid', grid_origin_alt=0, 
        min_radius = 100, nb = 1, bsp = 1, )#weighting_function = 'Nearest')

    
    return grid



def grid_EchoTop_max(radar, grid_shape=(80, 801, 801), xlim=(-150000, 150000), 
               ylim=(-150000, 150000), zlim=(1000, 20000),
               fields=['corrected_reflectivity', 'corrected_differential_reflectivity', 
                       'RainRate', 'HydroClass'], origin=None):
  
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('corrected_reflectivity')
    gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 10, 100)
    try:
        gatefilter.exclude_outside('normalized_coherent_power', 0.3, 1)
    except:
        print('no NCP')
    gatefilter.exclude_outside('corrected_copol_correlation_ratio', 0.9, 1)

    
    if origin is None:
        origin = (radar.latitude['data'][0],
                  radar.longitude['data'][0])
        
    radar_lower_sweep = radar#radar.extract_sweeps([0,1])
        
    grid = pyart.map.grid_from_radars(
        radar_lower_sweep, grid_shape=grid_shape,
        grid_limits=(zlim, ylim, xlim),
        grid_origin=origin, fields=fields,gatefilters = gatefilter,
        gridding_algo='map_gates_to_grid', grid_origin_alt=0, 
        min_radius = 250, nb = 0.5, bsp = 1, weighting_function = 'Nearest')

    
    return grid

def grid_EchoTop(radar, grid_shape=(20, 801, 801), xlim=(-150000, 150000), 
               ylim=(-150000, 150000), zlim=(1000, 20000),
               fields=['corrected_reflectivity', 'corrected_differential_reflectivity', 
                       'RainRate', 'HydroClass'], origin=None):
  
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('corrected_reflectivity')
    gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 10, 100)
    try:
        gatefilter.exclude_outside('normalized_coherent_power', 0.3, 1)
    except:
        print('no NCP')
    gatefilter.exclude_outside('corrected_copol_correlation_ratio', 0.9, 1)

    
    if origin is None:
        origin = (radar.latitude['data'][0],
                  radar.longitude['data'][0])
        
    radar_lower_sweep = radar#radar.extract_sweeps([0,1])
        
    grid = pyart.map.grid_from_radars(
        radar_lower_sweep, grid_shape=grid_shape,
        grid_limits=(zlim, ylim, xlim),
        grid_origin=origin, fields=fields,gatefilters = gatefilter,
        gridding_algo='map_gates_to_grid', grid_origin_alt=0, 
        min_radius = 500, nb = 1, bsp = 1, )#weighting_function = 'Nearest')

    
    return grid

def grid_rainRate_oneSweep(radar, grid_shape=(1, 801, 801), xlim=(-150000, 150000), 
               ylim=(-150000, 150000), zlim=(2500, 2500),
               fields=['corrected_reflectivity', 'corrected_differential_reflectivity', 
                       'RainRate', 'HydroClass'], origin=None):
  
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('corrected_reflectivity')
    gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 10, 100)
    gatefilter.exclude_outside('corrected_copol_correlation_ratio', 0.95, 1)
    try:
        gatefilter.exclude_outside('normalized_coherent_power', 0.3, 1)
    except:
        print('no NCP')
        gatefilter.exclude_outside('corrected_copol_correlation_ratio', 0.95, 1)
        gatefilter.exclude_outside('corrected_reflectivity', 20, 100)
    

    
    if origin is None:
        origin = (radar.latitude['data'][0],
                  radar.longitude['data'][0])
        
    radar_lower_sweep = radar #radar.extract_sweeps([1])
        
    grid = pyart.map.grid_from_radars(
        radar_lower_sweep, grid_shape=grid_shape,
        grid_limits=(zlim, ylim, xlim),
        grid_origin=origin, fields=fields,gatefilters = gatefilter,
        gridding_algo='map_gates_to_grid', grid_origin_alt=0, 
        min_radius = 2000, nb = 1, bsp = 1, )#weighting_function = 'Nearest')

    
    return grid

def grid_rainRate_oneSweep_original(radar, grid_shape=(1, 801, 801), xlim=(-150000, 150000), 
               ylim=(-150000, 150000), zlim=(2000, 2000),
               fields=['reflectivity',], origin=None):
    

    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    try:    
        gatefilter.exclude_invalid('reflectivity')
        #gatefilter.exclude_invalid('corrected_differential_phase')
        gatefilter.exclude_outside('reflectivity', 10, 100)
        #gatefilter.exclude_outside('normalized_coherent_power', 0.3, 1)
        #gatefilter.exclude_outside('copol_correlation_ratio', 0.9, 1)
    except:
        reflectivity = radar.fields['DBZH'].copy()
        radar.add_field('reflectivity', reflectivity, replace_existing=True)  
    
    if origin is None:
        origin = (radar.latitude['data'][0],
                  radar.longitude['data'][0])
        
    radar_lower_sweep = radar#radar.extract_sweeps([0,1])
        
    grid = pyart.map.grid_from_radars(
        radar_lower_sweep, grid_shape=grid_shape,
        grid_limits=(zlim, ylim, xlim),
        grid_origin=origin, fields=fields,gatefilters = gatefilter,
        gridding_algo='map_gates_to_grid', grid_origin_alt=0, 
        min_radius = 500, nb = 1, bsp = 1, weighting_function = 'Nearest')

    
    return grid

def grid_hail_oneSweep(radar, grid_shape=(1, 801, 801), xlim=(-150000, 150000), 
               ylim=(-150000, 150000), zlim=(1000, 1000),
               fields=['corrected_reflectivity', 'corrected_differential_reflectivity', 
                       'RainRate', 'HydroClass'], origin=None):
  
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('corrected_reflectivity')
    gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 40, 100)
    try:
        gatefilter.exclude_outside('normalized_coherent_power', 0.3, 1)
    except:
        print('no NCP')
        gatefilter.exclude_outside('corrected_reflectivity', 50, 100)
    #gatefilter.exclude_outside('normalized_coherent_power', 0.3, 1)
    #gatefilter.exclude_outside('corrected_copol_correlation_ratio', 0.9, 1)

    
    if origin is None:
        origin = (radar.latitude['data'][0],
                  radar.longitude['data'][0])
        
    radar_lower_sweep = radar.extract_sweeps([1])
    
    grid = pyart.map.grid_from_radars(
        radar_lower_sweep, grid_shape=grid_shape,
        grid_limits=(zlim, ylim, xlim),
        grid_origin=origin, fields=fields,gatefilters = gatefilter,
        gridding_algo='map_gates_to_grid', grid_origin_alt=0, 
        min_radius = 1000, nb = 1, bsp = 1, weighting_function = 'Nearest')

    
    return grid


def grid_dualRadar_original(chivo_filename, csapr_filename, grid_shape=(20, 301, 301), 
                   xlim=(-150000, 150000), ylim=(-150000, 150000), zlim=(1000, 20000),
               fields=['reflectivity', 'differential_reflectivity']):
    
    chivo_radar = pyart.io.read(chivo_filename)
    csapr_radar = pyart.io.read(csapr_filename)
  
    chivo_gatefilter = pyart.filters.GateFilter(chivo_radar)
    chivo_gatefilter.exclude_transition()
    chivo_gatefilter.exclude_invalid('reflectivity')
    chivo_gatefilter.exclude_invalid('differential_phase')
    chivo_gatefilter.exclude_outside('reflectivity', 0, 100)
    chivo_gatefilter.exclude_outside('normalized_coherent_power', 0.3, 1)
    #chivo_gatefilter.exclude_outside('corrected_copol_correlation_ratio', 0.9, 1)
    
    csapr_gatefilter = pyart.filters.GateFilter(csapr_radar)
    csapr_gatefilter.exclude_transition()
    csapr_gatefilter.exclude_invalid('reflectivity')
    csapr_gatefilter.exclude_invalid('differential_phase')
    csapr_gatefilter.exclude_outside('reflectivity', 0, 100)
    csapr_gatefilter.exclude_outside('normalized_coherent_power', 0.3, 1)
    #csapr_gatefilter.exclude_outside('corrected_copol_correlation_ratio', 0.9, 1)

    
    origin = (chivo_radar.latitude['data'][0], chivo_radar.longitude['data'][0])
        
    chivo_grid = pyart.map.grid_from_radars(chivo_radar, grid_shape=grid_shape,
        grid_limits=(zlim, ylim, xlim), grid_origin=origin, fields=fields,
        gatefilters = chivo_gatefilter, gridding_algo='map_gates_to_grid', grid_origin_alt=0, 
        min_radius = 100, nb = 1, bsp = 1)#, weighting_function = 'Nearest')
    
    csapr_grid = pyart.map.grid_from_radars(csapr_radar, grid_shape=grid_shape,
        grid_limits=(zlim, ylim, xlim), grid_origin=origin, fields=fields,
        gatefilters = csapr_gatefilter, gridding_algo='map_gates_to_grid', grid_origin_alt=0, 
        min_radius = 100, nb = 1, bsp = 1)#, weighting_function = 'Nearest')
    
    output_path = '/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Dec14/grided_data/'
    
    time_start = netCDF4.num2date(chivo_grid.time['data'][0], chivo_grid.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H%M')
    chivo_grid_filename = output_path + 'chivo_original_' + time_text + '.nc'
    
    time_start = netCDF4.num2date(csapr_grid.time['data'][0], csapr_grid.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H_%M')
    csapr_grid_filename = output_path + 'csapr_original_' + time_text + '.nc'
    
    pyart.io.write_grid(chivo_grid_filename, chivo_grid)
    pyart.io.write_grid(csapr_grid_filename, csapr_grid)



def grid_dualRadar_att_corrected(chivo_filename, csapr_filename, grid_shape=(20, 301, 301), 
                   xlim=(-150000, 150000), ylim=(-150000, 150000), zlim=(1000, 20000),
               fields=['corrected_reflectivity', 'corrected_differential_reflectivity', 
                       'RainRate', 'HydroClass']):
    
    chivo_radar = pyart.io.read(chivo_filename)
    csapr_radar = pyart.io.read(csapr_filename)
  
    chivo_gatefilter = pyart.filters.GateFilter(chivo_radar)
    chivo_gatefilter.exclude_transition()
    chivo_gatefilter.exclude_invalid('corrected_reflectivity')
    chivo_gatefilter.exclude_invalid('corrected_differential_phase')
    chivo_gatefilter.exclude_outside('corrected_reflectivity', 0, 100)
    chivo_gatefilter.exclude_outside('normalized_coherent_power', 0.3, 1)
    #chivo_gatefilter.exclude_outside('corrected_copol_correlation_ratio', 0.9, 1)
    
    csapr_gatefilter = pyart.filters.GateFilter(csapr_radar)
    csapr_gatefilter.exclude_transition()
    csapr_gatefilter.exclude_invalid('corrected_reflectivity')
    csapr_gatefilter.exclude_invalid('corrected_differential_phase')
    csapr_gatefilter.exclude_outside('corrected_reflectivity', 0, 100)
    csapr_gatefilter.exclude_outside('normalized_coherent_power', 0.3, 1)
    #csapr_gatefilter.exclude_outside('corrected_copol_correlation_ratio', 0.9, 1)

    
    origin = (chivo_radar.latitude['data'][0], chivo_radar.longitude['data'][0])
        
    chivo_grid = pyart.map.grid_from_radars(chivo_radar, grid_shape=grid_shape,
        grid_limits=(zlim, ylim, xlim), grid_origin=origin, fields=fields,
        gatefilters = chivo_gatefilter, gridding_algo='map_gates_to_grid', grid_origin_alt=0, 
        min_radius = 100, nb = 1, bsp = 1)#, weighting_function = 'Nearest')
    
    csapr_grid = pyart.map.grid_from_radars(csapr_radar, grid_shape=grid_shape,
        grid_limits=(zlim, ylim, xlim), grid_origin=origin, fields=fields,
        gatefilters = csapr_gatefilter, gridding_algo='map_gates_to_grid', grid_origin_alt=0, 
        min_radius = 100, nb = 1, bsp = 1)#, weighting_function = 'Nearest')
    
    output_path = '/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Dec14/grided_data/'
    
    time_start = netCDF4.num2date(chivo_grid.time['data'][0], chivo_grid.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H_%M')
    chivo_grid_filename = output_path + 'chivo_att_corrected_' + time_text + '.nc'
    
    time_start = netCDF4.num2date(csapr_grid.time['data'][0], csapr_grid.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H_%M')
    csapr_grid_filename = output_path + 'csapr_att_corrected_' + time_text + '.nc'
    
    pyart.io.write_grid(chivo_grid_filename, chivo_grid)
    pyart.io.write_grid(csapr_grid_filename, csapr_grid)

    

def grid_hail_csapr(radar, grid_shape=(4, 801, 801), xlim=(-120000, 120000), 
               ylim=(-120000, 120000), zlim=(0, 1000),
               fields=['corrected_reflectivity', 'corrected_differential_reflectivity', 
                       'RainRate', 'HydroClass'], origin=None):
  
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('corrected_reflectivity')
    gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 0, 100)
    #gatefilter.exclude_outside('normalized_coherent_power', 0.3, 1)
    #gatefilter.exclude_outside('corrected_copol_correlation_ratio', 0.9, 1)

    
    if origin is None:
        origin = (radar.latitude['data'][0],
                  radar.longitude['data'][0])
        
    grid = pyart.map.grid_from_radars(
        radar, grid_shape=grid_shape,
        grid_limits=(zlim, ylim, xlim),
        grid_origin=origin, fields=fields,gatefilter = gatefilter,
        gridding_algo='map_gates_to_grid', grid_origin_alt=0.0)

    
    return grid

def grid_hail_PPI(radar, grid_shape=(3, 801, 801), xlim=(-100000, 100000), 
               ylim=(-100000, 100000), zlim=(0, 3),
               fields=['corrected_reflectivity', 'corrected_differential_reflectivity', 
                       'RainRate', 'HydroClass'], origin=None):
  
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('corrected_reflectivity')
    gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 0, 100)
    #gatefilter.exclude_outside('normalized_coherent_power', 0.3, 1)
    #gatefilter.exclude_outside('corrected_copol_correlation_ratio', 0.9, 1)

    
    if origin is None:
        origin = (radar.latitude['data'][0],
                  radar.longitude['data'][0])
        
    grid = pyart.map.grid_from_radars(
        radar, grid_shape=grid_shape,
        grid_limits=(zlim, ylim, xlim),
        grid_origin=origin, fields=fields,gatefilter = gatefilter,
        gridding_algo='map_gates_to_grid', grid_origin_alt=0.0, weighting_function = 'Nearest')

    
    return grid

def grid_hail_hydro(radar, grid_shape=(1, 801, 801), xlim=(-100000, 100000), 
               ylim=(-100000, 100000), zlim=(0, 2000),
               fields=['corrected_reflectivity', 'corrected_differential_reflectivity', 
                       'RainRate', 'HydroClass'], origin=None):
  
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_invalid('corrected_reflectivity')
    gatefilter.exclude_invalid('corrected_differential_phase')
    gatefilter.exclude_outside('corrected_reflectivity', 0, 100)
    #gatefilter.exclude_outside('normalized_coherent_power', 0.3, 1)
    #gatefilter.exclude_outside('corrected_copol_correlation_ratio', 0.9, 1)

    
    if origin is None:
        origin = (radar.latitude['data'][0],
                  radar.longitude['data'][0])
        
    grid = pyart.map.grid_from_radars(
        radar, grid_shape=grid_shape,
        grid_limits=(zlim, ylim, xlim),
        grid_origin=origin, fields=fields,gatefilter = gatefilter,
        gridding_algo='map_gates_to_grid', grid_origin_alt=0.0, weighting_function = 'Nearest')

    
    return grid

