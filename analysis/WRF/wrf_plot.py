#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 15:29:40 2023

@author: idariash
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature

from wrf import (getvar, interplevel, to_np, latlon_coords, get_cartopy,
                 cartopy_xlim, cartopy_ylim, GeoBounds, CoordPair)

import pyart
import wrf
import os
import cmocean
import matplotlib.colors as mcolors


def max_wind_updraft(filename, timeidx):
    
  
    # Open the NetCDF file
    ncfile = Dataset(filename)  
    
    # Get the WRF variables
    terr = getvar(ncfile, "ter", timeidx = timeidx)
    max_dbz = getvar(ncfile, 'W', timeidx = timeidx) #"QRAIN"
        
    # Get the lat/lon coordinates
    lats, lons = latlon_coords(max_dbz)
    
    # Get the map projection information
    cart_proj = get_cartopy(max_dbz)
    
    # Create the figure
    fig = plt.figure(figsize=(12,9))
    ax = plt.axes(projection=cart_proj)
    
    # Download and add the states and coastlines
    states = NaturalEarthFeature(category="cultural", scale="50m",
                                 facecolor="none",
                                 name="admin_1_states_provinces_shp")
    ax.add_feature(states, linewidth=0.5, edgecolor="black")
    ax.coastlines('50m', linewidth=0.8)
    
    
    # Add reflectivity contours
    levels = np.arange(5., 70., 5.)
    maxRef_contours = plt.contourf(to_np(lons), to_np(lats), np.amax(to_np(max_dbz), 0),
                                 levels=levels,
                                 cmap=pyart.graph.cm.NWSRef,
                                 transform=crs.PlateCarree())
    plt.colorbar(maxRef_contours, ax=ax, orientation="vertical", pad=0.05)
    
    # Add terrain contours
    levels_terr = [1000,  2000]
    plt.contour(to_np(lons), to_np(lats), to_np(terr),
                                 levels=levels_terr,
                                linewidths = 1, colors = 'k',
                                  linestyles='solid', antialiased=True,
                                 #cmap=pyart.graph.cm.NWSRef,
                                 transform=crs.PlateCarree())
    
    Cordoba_lon, Cordoba_lat = -64.1686, -31.6342#Chivo -64.19, -31.42 # Cordoba
    plt.plot(Cordoba_lon, Cordoba_lat, 'ok', markersize = 5, transform=crs.PlateCarree())
    plt.text(Cordoba_lon, Cordoba_lat, ' C', fontsize=16, transform=crs.PlateCarree())
    
    ARM_lon, ARM_lat = -64.7284, -32.1264
    plt.plot(ARM_lon, ARM_lat, 'ok', markersize = 5, transform=crs.PlateCarree())
    plt.text(ARM_lon, ARM_lat, ' Y', fontsize=16, transform=crs.PlateCarree())
    
    BuenosAires_lon, BuenosAires_lat = -58.3816, -34.6037
    plt.plot(BuenosAires_lon, BuenosAires_lat, 'ok', markersize = 5, transform=crs.PlateCarree())
    plt.text(BuenosAires_lon, BuenosAires_lat, ' B', fontsize=16, transform=crs.PlateCarree())
    
    
    # Set the map bounds
    delta = 4.5
    REL_domain = GeoBounds(CoordPair(lat=Cordoba_lat - delta*(1 + 0.0), lon= Cordoba_lon - delta*(1 - 0.0)),
                       CoordPair(lat=Cordoba_lat + delta*(1 - 0.0), lon= Cordoba_lon + delta*(1 + 0.0)))
    # REL_domain = GeoBounds(CoordPair(lat=Cordoba_lat - delta*(1 + 0.5), lon= Cordoba_lon - delta*(1 - 0.5)),
    #                    CoordPair(lat=Cordoba_lat + delta*(1 - 0.5), lon= Cordoba_lon + delta*(1 + 0.5)))
    ax.set_xlim(cartopy_xlim(max_dbz, geobounds=REL_domain))
    ax.set_ylim(cartopy_ylim(max_dbz, geobounds=REL_domain))
    
    ax.gridlines()
    
    time = wrf.extract_times(ncfile, timeidx = timeidx)
    time_text = str(time)[:16]
    
    plt.title('WRF output  |  '+ time_text+ ' \n Max. Wind Updraft (m/s)')
       
    return fig

# Test
    
# filename = '/glade/scratch/ivanarias/WRF_test/domains/RELAMPAGO_IOP17_Aerosol/wrfout_d03_2018-12-14_03:00:00'
# #filename = '/glade/scratch/ivanarias/WRF_test/domains/RELAMPAGO_IOP17_Aerosol/ndown_wrf/wrfout_d01_2018-12-14_03:31:00'

# fig = max_wind_updraft(filename, 3)

def accumulated_total_precipitation(filename, timeidx):
    
    initial_filename = os.path.join(os.path.dirname(filename), 'wrfout_d03_2018-12-13_19:00:00')
    
    # Open the NetCDF file
    ncfile_initial_precip = Dataset(initial_filename) 
    ncfile = Dataset(filename)  
    
    # Get the WRF variables
    terr = getvar(ncfile, "ter", timeidx = timeidx)
    total_initial_precip = getvar(ncfile_initial_precip, 'RAINNC', timeidx = 5)
    total_precip = getvar(ncfile, 'RAINNC', timeidx = timeidx) #"QRAIN"
    
    # substract preveous preciptation
    total_precip.values = total_precip.values -  total_initial_precip.values
        
    # Get the lat/lon coordinates
    lats, lons = latlon_coords(total_precip)
    
    # Get the map projection information
    cart_proj = get_cartopy(total_precip)
    
    # Create the figure
    fig = plt.figure(figsize=(12,9))
    ax = plt.axes(projection=cart_proj)
    
    # Download and add the states and coastlines
    states = NaturalEarthFeature(category="cultural", scale="50m",
                                 facecolor="none",
                                 name="admin_1_states_provinces_shp")
    ax.add_feature(states, linewidth=0.5, edgecolor="black")
    ax.coastlines('50m', linewidth=0.8)
    
    
    # Add reflectivity contours
    
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
    cmap = mcolors.ListedColormap(cmap_data, 'precip_rss')
    norm = mcolors.BoundaryNorm(clevs, cmap.N)
    
    levels = np.arange(0, 150, 5.)
    RainC_contours = plt.contourf(to_np(lons), to_np(lats), to_np(total_precip),
                                 levels, cmap=cmocean.cm.rain,# norm=norm,
                                 transform=crs.PlateCarree())
    
    cb = plt.colorbar(RainC_contours, ax=ax, orientation="vertical", pad=0.05)
    cb.set_label('precipitation (mm)', fontsize=13)
        
    # Add terrain contours
    levels_terr = [1000,  2000]
    plt.contour(to_np(lons), to_np(lats), to_np(terr),
                                 levels=levels_terr,
                                linewidths = 1, colors = 'k',
                                  linestyles='solid', antialiased=True,
                                 #cmap=pyart.graph.cm.NWSRef,
                                 transform=crs.PlateCarree())
    
    Cordoba_lon, Cordoba_lat = -64.1686, -31.6342#Chivo -64.19, -31.42 # Cordoba
    plt.plot(Cordoba_lon, Cordoba_lat, 'ok', markersize = 5, transform=crs.PlateCarree())
    plt.text(Cordoba_lon, Cordoba_lat, ' C', fontsize=16, transform=crs.PlateCarree())
    
    ARM_lon, ARM_lat = -64.7284, -32.1264
    plt.plot(ARM_lon, ARM_lat, 'ok', markersize = 5, transform=crs.PlateCarree())
    plt.text(ARM_lon, ARM_lat, ' Y', fontsize=16, transform=crs.PlateCarree())
    
    BuenosAires_lon, BuenosAires_lat = -58.3816, -34.6037
    plt.plot(BuenosAires_lon, BuenosAires_lat, 'ok', markersize = 5, transform=crs.PlateCarree())
    plt.text(BuenosAires_lon, BuenosAires_lat, ' B', fontsize=16, transform=crs.PlateCarree())
    
    
    # Set the map bounds
    delta = 4.5
    REL_domain = GeoBounds(CoordPair(lat=Cordoba_lat - delta*(1 + 0.5), lon= Cordoba_lon - delta*(1 - 0.5)),
                        CoordPair(lat=Cordoba_lat + delta*(1 - 0.5), lon= Cordoba_lon + delta*(1 + 0.5)))
    
    # REL_domain = GeoBounds(CoordPair(lat=Cordoba_lat - delta, lon= Cordoba_lon - delta),
    #                     CoordPair(lat=Cordoba_lat + delta, lon= Cordoba_lon + delta))
    
    ax.set_xlim(cartopy_xlim(total_precip, geobounds=REL_domain))
    ax.set_ylim(cartopy_ylim(total_precip, geobounds=REL_domain))
       
    ax.gridlines()
    
    time = wrf.extract_times(ncfile, timeidx = timeidx)
    time_text = str(time)[:16]
    
    plt.title('WRF output  |  '+ time_text+ ' \n Event accumulated total precip.')
       
    return fig


def horizontal_max_ref(filename, timeidx):
    
  
    # Open the NetCDF file
    ncfile = Dataset(filename)  
    
    # Get the WRF variables
    terr = getvar(ncfile, "ter", timeidx = timeidx)
    max_dbz = getvar(ncfile, 'mdbz', timeidx = timeidx) #"QRAIN"
        
    # Get the lat/lon coordinates
    lats, lons = latlon_coords(max_dbz)
    
    # Get the map projection information
    cart_proj = get_cartopy(max_dbz)
    
    # Create the figure
    fig = plt.figure(figsize=(12,9))
    ax = plt.axes(projection=cart_proj)
    
    # Download and add the states and coastlines
    states = NaturalEarthFeature(category="cultural", scale="50m",
                                 facecolor="none",
                                 name="admin_1_states_provinces_shp")
    ax.add_feature(states, linewidth=0.5, edgecolor="black")
    ax.coastlines('50m', linewidth=0.8)
    
    
    # Add reflectivity contours
    levels = np.arange(5., 75., 5.)
    maxRef_contours = plt.contourf(to_np(lons), to_np(lats), to_np(max_dbz),
                                 levels=levels,
                                 cmap=pyart.graph.cm.NWSRef,
                                 transform=crs.PlateCarree())
    plt.colorbar(maxRef_contours, ax=ax, orientation="vertical", pad=0.05)
    
    # Add terrain contours
    levels_terr = [1000,  2000]
    plt.contour(to_np(lons), to_np(lats), to_np(terr),
                                 levels=levels_terr,
                                linewidths = 1, colors = 'k',
                                  linestyles='solid', antialiased=True,
                                 #cmap=pyart.graph.cm.NWSRef,
                                 transform=crs.PlateCarree())
    
    Cordoba_lon, Cordoba_lat = -64.1686, -31.6342#Chivo -64.19, -31.42 # Cordoba
    plt.plot(Cordoba_lon, Cordoba_lat, 'ok', markersize = 5, transform=crs.PlateCarree())
    plt.text(Cordoba_lon, Cordoba_lat, ' C', fontsize=16, transform=crs.PlateCarree())
    
    ARM_lon, ARM_lat = -64.7284, -32.1264
    plt.plot(ARM_lon, ARM_lat, 'ok', markersize = 5, transform=crs.PlateCarree())
    plt.text(ARM_lon, ARM_lat, ' Y', fontsize=16, transform=crs.PlateCarree())
    
    BuenosAires_lon, BuenosAires_lat = -58.3816, -34.6037
    plt.plot(BuenosAires_lon, BuenosAires_lat, 'ok', markersize = 5, transform=crs.PlateCarree())
    plt.text(BuenosAires_lon, BuenosAires_lat, ' B', fontsize=16, transform=crs.PlateCarree())
    
    
    # Set the map bounds
    delta = 4.5
    REL_domain = GeoBounds(CoordPair(lat=Cordoba_lat - delta*(1 + 0.0), lon= Cordoba_lon - delta*(1 - 0.0)),
                       CoordPair(lat=Cordoba_lat + delta*(1 - 0.0), lon= Cordoba_lon + delta*(1 + 0.0)))
    # REL_domain = GeoBounds(CoordPair(lat=Cordoba_lat - delta*(1 + 0.5), lon= Cordoba_lon - delta*(1 - 0.5)),
    #                    CoordPair(lat=Cordoba_lat + delta*(1 - 0.5), lon= Cordoba_lon + delta*(1 + 0.5)))
    ax.set_xlim(cartopy_xlim(max_dbz, geobounds=REL_domain))
    ax.set_ylim(cartopy_ylim(max_dbz, geobounds=REL_domain))
    
    ax.gridlines()
    
    time = wrf.extract_times(ncfile, timeidx = timeidx)
    time_text = str(time)[:16]
    
    plt.title('WRF output  |  '+ time_text+ ' \n Max. Reflectivity (dBZ)')
       
    return fig

# filename = '/glade/scratch/ivanarias/WRF_test/domains/RELAMPAGO_IOP17_NSSL_2Moments/wrfout_d03_2018-12-14_06:00:00'

# fig = accumulated_total_precipitation(filename, 3)