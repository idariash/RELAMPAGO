#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 12:28:39 2019

@author: idariash

This file contians functions for radar climatology
"""

import sys
sys.path.append('/net/nasstore/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/plot')
import chivo_grid_utl as utl
import pyart
import numpy as np
import scipy.ndimage as spyi
import netCDF4
import gc



#-------------------Following function are for DROPS files--------------

def echo_top_2d(filename):
# Reference: "An Improved Method for Estimating Radar Echo-Top Height" 2012 for reference
    
    if ('HYDRO' in filename) or ('BIRDBATH' in filename) or ('RHI' in filename):
        return (-100, -1000, -1000)
    
    radar = pyart.io.read(filename)
    echo_top_height = 0
    angles = radar.fixed_angle['data'];
    for sweep in range(len(angles) - 1, -1, -1): #start from the highest sweep
        reflectivity = radar.get_field(sweep, 'corrected_reflectivity')
        rhohv = radar.get_field(sweep, 'corrected_copol_correlation_ratio')
        ncp = radar.get_field(sweep, 'normalized_coherent_power')
        x, y, z = radar.get_gate_x_y_z(sweep, edges=False)
        x /= 1000.0
        y /= 1000.0
        z /= 1000.0
        # calculate (R)ange
        R = np.sqrt(x ** 2 + y ** 2)
        
        reflectivity = spyi.gaussian_filter(reflectivity, sigma = 2)
        rhohv = spyi.gaussian_filter(rhohv, sigma = 2)
        
        # Masking
        z = np.ma.masked_where( z < echo_top_height, z)
        z = np.ma.masked_where( R > 100, z)
        z = np.ma.masked_where( reflectivity < 18, z)
        z = np.ma.masked_where( rhohv < 0.9 , z)
        z = np.ma.masked_where( ncp < 0.3 , z)
        echo_top_sweep = np.amax(np.amax(z))
        if echo_top_sweep > echo_top_height:
            echo_top_height = echo_top_sweep
            ind_max_z = np.unravel_index(np.argmax(z, axis=None), z.shape)
            x_echo_top_height = x[ind_max_z]
            y_echo_top_height = y[ind_max_z]

    del(radar)
    gc.collect()
    return (echo_top_height, x_echo_top_height, y_echo_top_height)


def max_echoTop_grid(filename):
# Reference: "An Improved Method for Estimating Radar Echo-Top Height" 2012 for reference
    
    if ('HYDRO' in filename) or ('BIRDBATH' in filename):
        return -100

    if ('RHI' in filename): #not
        return -100
    
    radar = pyart.io.read(filename)
    # Create grids
    grid_echoTop = utl.grid_EchoTop_max(radar = radar, origin = [-31.6342, -64.1686])
    ref_echoTop = grid_echoTop.fields['corrected_reflectivity']['data']
    y, z, x = np.meshgrid(0.001*grid_echoTop.y['data'],  
                          0.001*grid_echoTop.z['data'],  0.001*grid_echoTop.x['data'])
    z_echoTop = np.ma.masked_where(ref_echoTop < 18, z)
    echoTop = np.amax(z_echoTop, axis=0)
    echoTop= np.ma.masked_where( echoTop < 4, echoTop )
    
    # Mask for region of computation
    lons, lats = grid_echoTop.get_point_longitude_latitude()
    echoTop = np.ma.masked_where( lats < -32, echoTop)
    echoTop = np.ma.masked_where( lats > -31, echoTop)
    echoTop = np.ma.masked_where( lons < -65, echoTop)
    echoTop = np.ma.masked_where( lons > -64, echoTop)
    
    echoTop_max = np.amax(np.amax(echoTop))
    
    del(radar, grid_echoTop)
    gc.collect()
    return (echoTop_max)


def accumulated_RR_Hail(filename):
# Reference: "An Improved Method for Estimating Radar Echo-Top Height" 2012 for reference
    
    if ('HYDRO' in filename) or ('BIRDBATH' in filename):
        return (-100, -1000, -1000)

    if ('RHI' in filename): #not
        return (-100, -1000, -1000)
    
    radar = pyart.io.read(filename)
    # Create grids
    grid_rainRate = utl.grid_rainRate_oneSweep(radar = radar, origin = [-31.6342, -64.1686])
    grid_hail = utl.grid_hail_oneSweep(radar = radar, origin = [-31.6342, -64.1686])
    # Reads fields
    radar_rainRate = grid_rainRate.fields['RainRate']['data'][0]
    radar_HydroClass = grid_hail.fields['HydroClass']['data'][0]
    lons, lats = grid_rainRate.get_point_longitude_latitude()
    #lons, lats = np.meshgrid(radar_lons, radar_lats)
    # Mask for region of computation
    radar_rainRate = np.ma.masked_where( lats < -32, radar_rainRate)
    radar_rainRate = np.ma.masked_where( lats > -31, radar_rainRate)
    radar_rainRate = np.ma.masked_where( lons < -65, radar_rainRate)
    radar_rainRate = np.ma.masked_where( lons > -64, radar_rainRate)
    #radar_rainRate = np.ma.masked_where( radar_rainRate < 1, radar_rainRate)
    
    radar_HydroClass = np.ma.masked_where( lats < -32, radar_HydroClass)
    radar_HydroClass = np.ma.masked_where( lats > -31, radar_HydroClass)
    radar_HydroClass = np.ma.masked_where( lons < -65, radar_HydroClass)
    radar_HydroClass = np.ma.masked_where( lons > -64, radar_HydroClass)
    # Filled the gabs
    # radar_rainRate.filled(0)
    # radar_rainRate.mask = False
    radar_HydroClass.filled(0)
    radar_HydroClass.mask = False
    
    def hail_flag(HID):
        if HID < 10 or 11 < HID:
            return 0
        else:
            return 1
    
    for i in  range(len(radar_HydroClass)):
        for j in range (len(radar_HydroClass[i])):
            radar_HydroClass[i][j] = hail_flag(radar_HydroClass[i][j])
            
    radar_HydroClass = np.ma.masked_where( lats < -32, radar_HydroClass)
    radar_HydroClass = np.ma.masked_where( lats > -31, radar_HydroClass)
    radar_HydroClass = np.ma.masked_where( lons < -65, radar_HydroClass)
    radar_HydroClass = np.ma.masked_where( lons > -64, radar_HydroClass)        
    
    rainRate_accumulation = np.ma.sum(np.sum(radar_rainRate))
    hail_accumulation = np.ma.sum(np.sum(radar_HydroClass))
    
    del(radar, grid_rainRate, grid_hail)
    gc.collect()
    return (rainRate_accumulation, hail_accumulation)
    
    

def echo_top_distribution(filename):
# Reference: "An Improved Method for Estimating Radar Echo-Top Height" 2012 for reference
    
    if ('HYDRO' in filename) or ('BIRDBATH' in filename):
        return (-100, -1000, -1000)

    if ('RHI' in filename): #not
        return (-100, -1000, -1000)
    
    radar = pyart.io.read(filename)
    echo_top_height = 0
    angles = radar.fixed_angle['data'];
    for sweep in range(len(angles) - 1, -1, -1): #start from the highest sweep
        reflectivity = radar.get_field(sweep, 'corrected_reflectivity')
        rhohv = radar.get_field(sweep, 'corrected_copol_correlation_ratio')
        x, y, z = radar.get_gate_x_y_z(sweep, edges=False)
        x /= 1000.0
        y /= 1000.0
        z /= 1000.0
        # calculate (R)ange
        R = np.sqrt(x ** 2 + y ** 2)
        
        reflectivity = spyi.gaussian_filter(reflectivity, sigma = 2)
        rhohv = spyi.gaussian_filter(rhohv, sigma = 2)
        
        # Masking
        z = np.ma.masked_where( z < echo_top_height, z)
        z = np.ma.masked_where( R > 100, z)
        z = np.ma.masked_where( reflectivity < 18, z)
        z = np.ma.masked_where( rhohv < 0.9 , z)
        try:
             ncp = radar.get_field(sweep, 'normalized_coherent_power')
             z = np.ma.masked_where( ncp < 0.3 , z)
        except:
            print('No NCP ')

        echo_top_sweep = np.amax(np.amax(z))
        if echo_top_sweep > echo_top_height:
            echo_top_height = echo_top_sweep
            ind_max_z = np.unravel_index(np.argmax(z, axis=None), z.shape)
            x_echo_top_height = x[ind_max_z]
            y_echo_top_height = y[ind_max_z]

    del(radar)
    gc.collect()
    return (echo_top_height, x_echo_top_height, y_echo_top_height)


def echo_top_drops(filename):
# Check "An Improved Method for Estimating Radar Echo-Top Height" 2012 for reference
    radar = pyart.io.read(filename)
    echo_top_height = 0
    angles = radar.fixed_angle['data'];
    for sweep in range(len(angles) - 1, -1, -1): #start from the highest sweep
        reflectivity = radar.get_field(sweep, 'corrected_reflectivity')
        rhohv = radar.get_field(sweep, 'corrected_copol_correlation_ratio')
        x, y, z = radar.get_gate_x_y_z(sweep, edges=False)
        x /= 1000.0
        y /= 1000.0
        z /= 1000.0
        # calculate (R)ange
        R = np.sqrt(x ** 2 + y ** 2)
        
        reflectivity = spyi.gaussian_filter(reflectivity, sigma = 2)
        rhohv = spyi.gaussian_filter(rhohv, sigma = 2)
        
        # Masking
        z = np.ma.masked_where( z < echo_top_height, z)
        z = np.ma.masked_where( R > 100, z)
        z = np.ma.masked_where( reflectivity < 18, z)
        z = np.ma.masked_where( rhohv < 0.9 , z)
        echo_top_sweep = np.amax(np.amax(z))
        if echo_top_sweep > echo_top_height:
            echo_top_height = echo_top_sweep

    del(radar)
    gc.collect()
    return echo_top_height

def max_reflectivity_drops(filename):
    radar = pyart.io.read(filename)
    max_reflectivity = -100
    angles = radar.fixed_angle['data'];
    for sweep in range(0, len(angles)):
        reflectivity = radar.get_field(sweep, 'corrected_reflectivity')
        rhohv = radar.get_field(sweep, 'corrected_copol_correlation_ratio')
        x, y, z = radar.get_gate_x_y_z(sweep, edges=False)
        x /= 1000.0
        y /= 1000.0
        z /= 1000.0
        # calculate (R)ange
        R = np.sqrt(x ** 2 + y ** 2)
        
        reflectivity = spyi.gaussian_filter(reflectivity, sigma = 2)
        rhohv = spyi.gaussian_filter(rhohv, sigma = 2)
        
        # Masking
        reflectivity = np.ma.masked_where( reflectivity < max_reflectivity, reflectivity)
        reflectivity = np.ma.masked_where( R > 100, reflectivity)
        reflectivity = np.ma.masked_where( z < 2, reflectivity)
        reflectivity = np.ma.masked_where( rhohv < 0.8 , reflectivity)
        max_reflectivity_sweep = np.amax(np.amax(reflectivity))
        if max_reflectivity_sweep > max_reflectivity:
            max_reflectivity = max_reflectivity_sweep
    
    del(radar)
    gc.collect()
    return max_reflectivity

def max_Kdp_drops(filename):
    radar = pyart.io.read(filename)
    max_Kdp = 0
    angles = radar.fixed_angle['data'];
    for sweep in range(0, len(angles)):
        Kdp = radar.get_field(sweep, 'corrected_specific_differential_phase')
        rhohv = radar.get_field(sweep, 'corrected_copol_correlation_ratio')
        x, y, z = radar.get_gate_x_y_z(sweep, edges=False)
        x /= 1000.0
        y /= 1000.0
        z /= 1000.0
        # calculate (R)ange
        R = np.sqrt(x ** 2 + y ** 2)
        
        Kdp = spyi.gaussian_filter(Kdp, sigma = 2)
        rhohv = spyi.gaussian_filter(rhohv, sigma = 2)
        
        # Masking
        Kdp = np.ma.masked_where( Kdp < max_Kdp, Kdp)
        Kdp = np.ma.masked_where( R > 100, Kdp)
        Kdp = np.ma.masked_where( z < 2, Kdp)
        Kdp = np.ma.masked_where( rhohv < 0.8 , Kdp)
        max_Kdp_sweep = np.amax(np.amax(Kdp))
        if max_Kdp_sweep > max_Kdp:
            max_Kdp = max_Kdp_sweep
    
    return max_Kdp

def max_Kdp_drops_filterCloudBoundary(filename):
    radar = pyart.io.read(filename)
    max_Kdp = 0
    angles = radar.fixed_angle['data'];
    for sweep in range(0, len(angles)):
        Kdp = radar.get_field(sweep, 'corrected_specific_differential_phase')
        rhohv = radar.get_field(sweep, 'corrected_copol_correlation_ratio')
        reflectivity = radar.get_field(sweep, 'corrected_reflectivity')
        x, y, z = radar.get_gate_x_y_z(sweep, edges=False)
        x /= 1000.0
        y /= 1000.0
        z /= 1000.0
        # calculate (R)ange
        R = np.sqrt(x ** 2 + y ** 2)
        
        reflectivity = spyi.gaussian_filter(reflectivity, sigma = 2)
        Kdp = spyi.gaussian_filter(Kdp, sigma = 2)
        rhohv = spyi.gaussian_filter(rhohv, sigma = 2)
        
        # Masking
        Kdp = np.ma.masked_where( Kdp < max_Kdp, Kdp)
        Kdp = np.ma.masked_where( R > 100, Kdp)
        Kdp = np.ma.masked_where( z < 2, Kdp)
        Kdp = np.ma.masked_where( z > 8, Kdp)
        Kdp = np.ma.masked_where( rhohv < 0.8 , Kdp)
        Kdp = np.ma.masked_where( reflectivity < 25 , Kdp)
        max_Kdp_sweep = np.amax(np.amax(Kdp))
        if max_Kdp_sweep > max_Kdp:
            max_Kdp = max_Kdp_sweep
    
    return max_Kdp

def max_Kdp_drops_filterCloudBoundary_manyReturns(filename):
    radar = pyart.io.read(filename)
    max_Kdp = 0
    angles = radar.fixed_angle['data'];
    for sweep in range(0, len(angles)):
        Kdp = radar.get_field(sweep, 'corrected_specific_differential_phase')
        rhohv = radar.get_field(sweep, 'corrected_copol_correlation_ratio')
        reflectivity = radar.get_field(sweep, 'corrected_reflectivity')
        x, y, z = radar.get_gate_x_y_z(sweep, edges=False)
        x /= 1000.0
        y /= 1000.0
        z /= 1000.0
        # calculate (R)ange
        R = np.sqrt(x ** 2 + y ** 2)
        
        reflectivity = spyi.gaussian_filter(reflectivity, sigma = 2)
        Kdp = spyi.gaussian_filter(Kdp, sigma = 2)
        rhohv = spyi.gaussian_filter(rhohv, sigma = 2)
        
        # Masking
        Kdp = np.ma.masked_where( Kdp < max_Kdp, Kdp)
        Kdp = np.ma.masked_where( R > 100, Kdp)
        Kdp = np.ma.masked_where( z < 2, Kdp)
        Kdp = np.ma.masked_where( z > 8, Kdp)
        Kdp = np.ma.masked_where( rhohv < 0.8 , Kdp)
        Kdp = np.ma.masked_where( reflectivity < 25 , Kdp)
        max_Kdp_sweep = np.amax(np.amax(Kdp))
        if max_Kdp_sweep > max_Kdp:
            max_Kdp = max_Kdp_sweep
            max_sweepAngle = angles[sweep]
            max_index = np.unravel_index(Kdp.argmax(), Kdp.shape)
            max_height = z[max_index]
            max_range = R[max_index]
            
    return (max_Kdp, max_sweepAngle, max_height, max_range)

def max_Kdp_info(filename):
    radar = pyart.io.read(filename)
    max_Kdp = 0
    angles = radar.fixed_angle['data']
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.isoformat()
    for sweep in range(0, len(angles)):
        Kdp = radar.get_field(sweep, 'corrected_specific_differential_phase')
        rhohv = radar.get_field(sweep, 'corrected_copol_correlation_ratio')
        reflectivity = radar.get_field(sweep, 'corrected_reflectivity')
        x, y, z = radar.get_gate_x_y_z(sweep, edges=False)
        x /= 1000.0
        y /= 1000.0
        z /= 1000.0
        # calculate (R)ange
        R = np.sqrt(x ** 2 + y ** 2)
        
        reflectivity = spyi.gaussian_filter(reflectivity, sigma = 2)
        Kdp = spyi.gaussian_filter(Kdp, sigma = 2)
        rhohv = spyi.gaussian_filter(rhohv, sigma = 2)
        
        # Masking
        Kdp = np.ma.masked_where( Kdp < max_Kdp, Kdp)
        Kdp = np.ma.masked_where( R > 100, Kdp)
        Kdp = np.ma.masked_where( z < 2, Kdp)
        Kdp = np.ma.masked_where( z > 8, Kdp)
        #Kdp = np.ma.masked_where( rhohv < 0.8 , Kdp)
        Kdp = np.ma.masked_where( reflectivity < 25 , Kdp)
        max_Kdp_sweep = np.amax(np.amax(Kdp))
        if max_Kdp_sweep > max_Kdp:
            max_Kdp = max_Kdp_sweep
            max_sweepAngle = angles[sweep]
            max_index = np.unravel_index(Kdp.argmax(), Kdp.shape)
            max_height = z[max_index]
            max_range = R[max_index]
            max_ref = reflectivity[max_index]
            max_rhohv = rhohv[max_index]
            sweep_number = sweep
            
    return {'time_UTC': time_text, 'Kdp_max': max_Kdp, 'sweep_number': sweep_number,
            'sweep_angle': max_sweepAngle, 'height': max_height, 'range': max_range,
            'reflectivity': max_ref, 'rhohv': max_rhohv
            }

def max_ref_info(filename, sweepNumber):
    radar = pyart.io.read(filename)
    max_ref = 0
    angles = radar.fixed_angle['data']
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.isoformat()
    for sweep in range(sweepNumber - 1, sweepNumber):
        Kdp = radar.get_field(sweep, 'corrected_specific_differential_phase')
        rhohv = radar.get_field(sweep, 'corrected_copol_correlation_ratio')
        reflectivity = radar.get_field(sweep, 'corrected_reflectivity')
        x, y, z = radar.get_gate_x_y_z(sweep, edges=False)
        x /= 1000.0
        y /= 1000.0
        z /= 1000.0
        # calculate (R)ange
        R = np.sqrt(x ** 2 + y ** 2)
        
        reflectivity = spyi.gaussian_filter(reflectivity, sigma = 2)
        Kdp = spyi.gaussian_filter(Kdp, sigma = 2)
        rhohv = spyi.gaussian_filter(rhohv, sigma = 2)
        
        # Masking
        #Kdp = np.ma.masked_where( Kdp < max_Kdp, Kdp)
        reflectivity = np.ma.masked_where( R > 30, reflectivity)
        reflectivity = np.ma.masked_where( z < 3, reflectivity)
        Kdp = np.ma.masked_where( z > 8, Kdp)
        #Kdp = np.ma.masked_where( rhohv < 0.8 , Kdp)
        #reflectivity = np.ma.masked_where( reflectivity < 25 , Kdp)
        max_ref_sweep = np.amax(np.amax(reflectivity))
        if max_ref_sweep > max_ref:
            max_sweepAngle = angles[sweep]
            max_index = np.unravel_index(reflectivity.argmax(), reflectivity.shape)
            max_height = z[max_index]
            max_range = R[max_index]
            max_ref = reflectivity[max_index]
            #max_rhohv = rhohv[max_index]
            sweep_number = sweep
            
    return {'time_UTC': time_text,  'sweep_number': sweep_number,
            'sweep_angle': max_sweepAngle, 'height': max_height, 'range': max_range,
            'reflectivity': max_ref
            }
    
def max_ref_info_noDrops(filename):
    radar = pyart.io.read(filename)
    max_ref = 0
    angles = radar.fixed_angle['data']
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.isoformat()
    scan_type = radar.scan_type
    for sweep in range(0, len(angles)):
        #rhohv = radar.get_field(sweep, 'cross_correlation_ratio')
        reflectivity = radar.get_field(sweep, 'reflectivity')
        x, y, z = radar.get_gate_x_y_z(sweep, edges=False)
        x /= 1000.0
        y /= 1000.0
        z /= 1000.0
        # calculate (R)ange
        R = np.sqrt(x ** 2 + y ** 2)
        
        reflectivity = spyi.gaussian_filter(reflectivity, sigma = 2)
        #rhohv = spyi.gaussian_filter(rhohv, sigma = 2)
        
        # Masking
        reflectivity = np.ma.masked_where( R > 80, reflectivity)
        reflectivity = np.ma.masked_where( z < 3, reflectivity)
        max_ref_sweep = np.amax(np.amax(reflectivity))
        if max_ref_sweep > max_ref:
            max_sweepAngle = angles[sweep]
            max_index = np.unravel_index(reflectivity.argmax(), reflectivity.shape)
            max_height = z[max_index]
            max_range = R[max_index]
            max_ref = reflectivity[max_index]
            #max_rhohv = rhohv[max_index]
            sweep_number = sweep
    try:
        return {'time_UTC': time_text,  'scan_type': scan_type, 'reflectivity': max_ref , 
                'height': max_height, 'range': max_range, 
                'sweep_number': sweep_number, 'sweep_angle': max_sweepAngle
                }
    except:
        return {'time_UTC': time_text, 'scan_type': scan_type, 'reflectivity': -100 , 'height': -1,
                'range': -1, 'sweep_number': -1, 'sweep_angle': -1
                }

def max_ref_info_drops(filename):
    radar = pyart.io.read(filename)
    max_ref = 0
    angles = radar.fixed_angle['data']
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.isoformat()
    for sweep in range(0, len(angles)):
        #Kdp = radar.get_field(sweep, 'corrected_specific_differential_phase')
        #rhohv = radar.get_field(sweep, 'corrected_copol_correlation_ratio')
        reflectivity = radar.get_field(sweep, 'corrected_reflectivity')
        x, y, z = radar.get_gate_x_y_z(sweep, edges=False)
        x /= 1000.0
        y /= 1000.0
        z /= 1000.0
        # calculate (R)ange
        R = np.sqrt(x ** 2 + y ** 2)
        
        reflectivity = spyi.gaussian_filter(reflectivity, sigma = 2)
        #Kdp = spyi.gaussian_filter(Kdp, sigma = 2)
        #rhohv = spyi.gaussian_filter(rhohv, sigma = 2)
        
        # Masking
        #Kdp = np.ma.masked_where( reflectivity < 25 , Kdp)
        reflectivity = np.ma.masked_where( R > 70, reflectivity)
        reflectivity = np.ma.masked_where( z < 5, reflectivity)
        #Kdp = np.ma.masked_where( z > 8, Kdp)
        #Kdp = np.ma.masked_where( rhohv < 0.8 , Kdp)
        max_ref_sweep = np.amax(np.amax(reflectivity))
        if max_ref_sweep > max_ref:
            max_sweepAngle = angles[sweep]
            max_index = np.unravel_index(reflectivity.argmax(), reflectivity.shape)
            max_height = z[max_index]
            max_range = R[max_index]
            max_ref = reflectivity[max_index]
            #max_rhohv = rhohv[max_index]
            sweep_number = sweep
    try:            
        return {'time_UTC': time_text,  'reflectivity': max_ref , 'height': max_height,
                'range': max_range, 'sweep_number': sweep_number, 'sweep_angle': max_sweepAngle
                }
    except:
        return {'time_UTC': time_text,  'reflectivity': -100 , 'height': -1,
                'range': -1, 'sweep_number': -1, 'sweep_angle': -1
                }

def max_Zdr_drops(filename):
    radar = pyart.io.read(filename)
    max_Zdr = -100
    angles = radar.fixed_angle['data'];
    for sweep in range(0, len(angles)):
        Zdr = radar.get_field(sweep, 'corrected_differential_reflectivity')
        reflectivity = radar.get_field(sweep, 'corrected_reflectivity')
        rhohv = radar.get_field(sweep, 'corrected_copol_correlation_ratio')
        x, y, z = radar.get_gate_x_y_z(sweep, edges=False)
        x /= 1000.0
        y /= 1000.0
        z /= 1000.0
        # calculate (R)ange
        R = np.sqrt(x ** 2 + y ** 2)
        
        Zdr = spyi.gaussian_filter(Zdr, sigma = 2)
        rhohv = spyi.gaussian_filter(rhohv, sigma = 2)
        
        # Masking
        Zdr = np.ma.masked_where( Zdr < max_Zdr, Zdr)
        Zdr = np.ma.masked_where( R > 100, Zdr)
        Zdr = np.ma.masked_where( z < 2, Zdr)
        Zdr = np.ma.masked_where( rhohv < 0.8 , Zdr)
        Zdr = np.ma.masked_where( reflectivity < 20 , Zdr)
        max_Zdr_sweep = np.amax(np.amax(Zdr))
        if max_Zdr_sweep > max_Zdr:
            max_Zdr = max_Zdr_sweep
    
    return max_Zdr

def min_Zdr_drops(filename):
    radar = pyart.io.read(filename)
    min_Zdr = 100
    angles = radar.fixed_angle['data'];
    for sweep in range(0, len(angles)):
        Zdr = radar.get_field(sweep, 'corrected_differential_reflectivity')
        reflectivity = radar.get_field(sweep, 'corrected_reflectivity')
        rhohv = radar.get_field(sweep, 'corrected_copol_correlation_ratio')
        x, y, z = radar.get_gate_x_y_z(sweep, edges=False)
        x /= 1000.0
        y /= 1000.0
        z /= 1000.0
        # calculate (R)ange
        R = np.sqrt(x ** 2 + y ** 2)
        
        Zdr = spyi.gaussian_filter(Zdr, sigma = 2)
        rhohv = spyi.gaussian_filter(rhohv, sigma = 2)
        
        # Masking
        Zdr = np.ma.masked_where( Zdr > min_Zdr, Zdr)
        Zdr = np.ma.masked_where( R > 100, Zdr)
        Zdr = np.ma.masked_where( z < 2, Zdr)
        Zdr = np.ma.masked_where( rhohv < 0.9 , Zdr)
        Zdr = np.ma.masked_where( reflectivity < 25 , Zdr)
        min_Zdr_sweep = np.amin(np.amin(Zdr))
        if min_Zdr_sweep < min_Zdr:
            min_Zdr = min_Zdr_sweep
    
    return min_Zdr




def echo_top(filename):
# Check "An Improved Method for Estimating Radar Echo-Top Height" 2012 for reference
    radar = pyart.io.read(filename)
    echo_top_height = 0
    angles = radar.fixed_angle['data'];
    for sweep in range(len(angles) - 1, -1, -1): #start from the highest sweep
        reflectivity = radar.get_field(sweep, 'reflectivity')
        rhohv = radar.get_field(sweep, 'cross_correlation_ratio')
        x, y, z = radar.get_gate_x_y_z(sweep, edges=False)
        x /= 1000.0
        y /= 1000.0
        z /= 1000.0
        # calculate (R)ange
        R = np.sqrt(x ** 2 + y ** 2)
        
        reflectivity = spyi.gaussian_filter(reflectivity, sigma = 2)
        rhohv = spyi.gaussian_filter(rhohv, sigma = 2)
        
        # Masking
        z = np.ma.masked_where( z < echo_top_height, z)
        z = np.ma.masked_where( R > 100, z)
        z = np.ma.masked_where( reflectivity < 18, z)
        z = np.ma.masked_where( rhohv < 0.9 , z)
        echo_top_sweep = np.amax(np.amax(z))
        if echo_top_sweep > echo_top_height:
            echo_top_height = echo_top_sweep
    
    return echo_top_height

def max_reflectivity(filename):
    radar = pyart.io.read(filename)
    max_reflectivity = -100
    angles = radar.fixed_angle['data'];
    for sweep in range(0, len(angles)):
        reflectivity = radar.get_field(sweep, 'reflectivity')
        rhohv = radar.get_field(sweep, 'cross_correlation_ratio')
        x, y, z = radar.get_gate_x_y_z(sweep, edges=False)
        x /= 1000.0
        y /= 1000.0
        z /= 1000.0
        # calculate (R)ange
        R = np.sqrt(x ** 2 + y ** 2)
        
        reflectivity = spyi.gaussian_filter(reflectivity, sigma = 2)
        rhohv = spyi.gaussian_filter(rhohv, sigma = 2)
        
        # Masking
        reflectivity = np.ma.masked_where( reflectivity < max_reflectivity, reflectivity)
        reflectivity = np.ma.masked_where( R > 100, reflectivity)
        reflectivity = np.ma.masked_where( z < 2, reflectivity)
        reflectivity = np.ma.masked_where( rhohv < 0.8 , reflectivity)
        max_reflectivity_sweep = np.amax(np.amax(reflectivity))
        if max_reflectivity_sweep > max_reflectivity:
            max_reflectivity = max_reflectivity_sweep
    
    return max_reflectivity

def max_Kdp(filename):
    radar = pyart.io.read(filename)
    max_Kdp = 0
    angles = radar.fixed_angle['data'];
    for sweep in range(0, len(angles)):
        Kdp = radar.get_field(sweep, 'specific_differential_phase')
        rhohv = radar.get_field(sweep, 'cross_correlation_ratio')
        x, y, z = radar.get_gate_x_y_z(sweep, edges=False)
        x /= 1000.0
        y /= 1000.0
        z /= 1000.0
        # calculate (R)ange
        R = np.sqrt(x ** 2 + y ** 2)
        
        Kdp = spyi.gaussian_filter(Kdp, sigma = 2)
        rhohv = spyi.gaussian_filter(rhohv, sigma = 2)
        
        # Masking
        Kdp = np.ma.masked_where( Kdp < max_Kdp, Kdp)
        Kdp = np.ma.masked_where( R > 100, Kdp)
        Kdp = np.ma.masked_where( z < 2, Kdp)
        Kdp = np.ma.masked_where( rhohv < 0.8 , Kdp)
        max_Kdp_sweep = np.amax(np.amax(Kdp))
        if max_Kdp_sweep > max_Kdp:
            max_Kdp = max_Kdp_sweep
    
    return max_Kdp

def max_Zdr(filename):
    radar = pyart.io.read(filename)
    max_Zdr = -100
    angles = radar.fixed_angle['data'];
    for sweep in range(0, len(angles)):
        Zdr = radar.get_field(sweep, 'differential_reflectivity')
        reflectivity = radar.get_field(sweep, 'reflectivity')
        rhohv = radar.get_field(sweep, 'cross_correlation_ratio')
        x, y, z = radar.get_gate_x_y_z(sweep, edges=False)
        x /= 1000.0
        y /= 1000.0
        z /= 1000.0
        # calculate (R)ange
        R = np.sqrt(x ** 2 + y ** 2)
        
        Zdr = spyi.gaussian_filter(Zdr, sigma = 2)
        rhohv = spyi.gaussian_filter(rhohv, sigma = 2)
        
        # Masking
        Zdr = np.ma.masked_where( Zdr < max_Zdr, Zdr)
        Zdr = np.ma.masked_where( R > 100, Zdr)
        Zdr = np.ma.masked_where( z < 2, Zdr)
        Zdr = np.ma.masked_where( rhohv < 0.8 , Zdr)
        Zdr = np.ma.masked_where( reflectivity < 20 , Zdr)
        max_Zdr_sweep = np.amax(np.amax(Zdr))
        if max_Zdr_sweep > max_Zdr:
            max_Zdr = max_Zdr_sweep
    
    return max_Zdr




#--------------- CSAPR Functions -----------------------
def csapr_get_maxVal(filename):
    radar = pyart.io.read(filename)
    max_Kdp = 0.1
    max_reflectivity = -100
    echo_top_height = 0
    angles = radar.fixed_angle['data']
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.isoformat()
    for sweep in range(0, len(angles)):
        Kdp = radar.get_field(sweep, 'corrected_specific_differential_phase')
        rhohv = radar.get_field(sweep, 'corrected_copol_correlation_ratio')
        reflectivity = radar.get_field(sweep, 'corrected_reflectivity')
        x, y, z = radar.get_gate_x_y_z(sweep, edges=False)
        x /= 1000.0
        y /= 1000.0
        z /= 1000.0
        # calculate (R)ange
        R = np.sqrt(x ** 2 + y ** 2)
        
        reflectivity = spyi.gaussian_filter(reflectivity, sigma = 2)
        Kdp = spyi.gaussian_filter(Kdp, sigma = 2)
        rhohv = spyi.gaussian_filter(rhohv, sigma = 2)
        
        # Masking Kdp
        Kdp = np.ma.masked_where( Kdp < max_Kdp, Kdp)
        Kdp = np.ma.masked_where( R > 100, Kdp)
        Kdp = np.ma.masked_where( z < 2, Kdp)
        Kdp = np.ma.masked_where( z > 8, Kdp)
        #Kdp = np.ma.masked_where( rhohv < 0.8 , Kdp)
        Kdp = np.ma.masked_where( reflectivity < 25 , Kdp)
        
        # Masking height for echo top
        z = np.ma.masked_where( z < echo_top_height, z)
        z = np.ma.masked_where( R > 100, z)
        z = np.ma.masked_where( reflectivity < 18, z)
        z = np.ma.masked_where( rhohv < 0.9 , z)
        
        # Masking for reflectivity
        reflectivity = np.ma.masked_where( reflectivity < max_reflectivity, reflectivity)
        reflectivity = np.ma.masked_where( R > 100, reflectivity)
        #reflectivity = np.ma.masked_where( z < 2, reflectivity)
        reflectivity = np.ma.masked_where( rhohv < 0.8 , reflectivity)
        
        max_Kdp_sweep = np.amax(np.amax(Kdp))
        max_reflectivity_sweep = np.amax(np.amax(reflectivity))
        echo_top_sweep = np.amax(np.amax(z))
        
        if max_Kdp_sweep > max_Kdp:
            max_Kdp = max_Kdp_sweep
            
        if echo_top_sweep > echo_top_height:
            echo_top_height = echo_top_sweep
            
        if max_reflectivity_sweep > max_reflectivity:
            max_reflectivity = max_reflectivity_sweep
            
    return {'time_utc': time_text, 'echo_top': echo_top_height, 
            'max_reflectivity': max_reflectivity, 'max_Kdp': max_Kdp
            }



def echo_top_csapr(filename):
# Check "An Improved Method for Estimating Radar Echo-Top Height" 2012 for reference
    radar = pyart.io.read(filename)
    echo_top_height = 0
    angles = radar.fixed_angle['data'];
    for sweep in range(len(angles) - 1, -1, -1): #start from the highest sweep
        reflectivity = radar.get_field(sweep, 'reflectivity')
        rhohv = radar.get_field(sweep, 'copol_correlation_coeff')
        x, y, z = radar.get_gate_x_y_z(sweep, edges=False)
        x /= 1000.0
        y /= 1000.0
        z /= 1000.0
        # calculate (R)ange
        R = np.sqrt(x ** 2 + y ** 2)
        
        reflectivity = spyi.gaussian_filter(reflectivity, sigma = 2)
        rhohv = spyi.gaussian_filter(rhohv, sigma = 2)
        
        # Masking
        z = np.ma.masked_where( z < echo_top_height, z)
        z = np.ma.masked_where( R > 100, z)
        z = np.ma.masked_where( reflectivity < 18, z)
        z = np.ma.masked_where( rhohv < 0.9 , z)
        echo_top_sweep = np.amax(np.amax(z))
        if echo_top_sweep > echo_top_height:
            echo_top_height = echo_top_sweep

    del(radar)
    gc.collect()
    return echo_top_height

def max_reflectivity_csapr(filename):
    radar = pyart.io.read(filename)
    max_reflectivity = -100
    angles = radar.fixed_angle['data'];
    for sweep in range(0, len(angles)):
        reflectivity = radar.get_field(sweep, 'reflectivity')
        rhohv = radar.get_field(sweep, 'copol_correlation_coeff')
        x, y, z = radar.get_gate_x_y_z(sweep, edges=False)
        x /= 1000.0
        y /= 1000.0
        z /= 1000.0
        # calculate (R)ange
        R = np.sqrt(x ** 2 + y ** 2)
        
        reflectivity = spyi.gaussian_filter(reflectivity, sigma = 2)
        rhohv = spyi.gaussian_filter(rhohv, sigma = 2)
        
        # Masking
        reflectivity = np.ma.masked_where( reflectivity < max_reflectivity, reflectivity)
        reflectivity = np.ma.masked_where( R > 100, reflectivity)
        reflectivity = np.ma.masked_where( z < 2, reflectivity)
        reflectivity = np.ma.masked_where( rhohv < 0.8 , reflectivity)
        max_reflectivity_sweep = np.amax(np.amax(reflectivity))
        if max_reflectivity_sweep > max_reflectivity:
            max_reflectivity = max_reflectivity_sweep
    
    del(radar)
    gc.collect()
    return max_reflectivity







#------- Functions for level 1a data ---------------------
    
def max_ref_info_level1a(filename):
    radar = pyart.io.read(filename)
    max_ref = -100
    angles = radar.fixed_angle['data']
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.isoformat()
    print(time_text)
    # Following is to make sure is a volume scan
    if len(angles) < 3 or radar.scan_type == 'rhi' :
        return {'time_UTC': time_text,  'reflectivity': -100, 
                'NorthSouth_distance': 1, 'EastWest_distance': 1, 'height': -1, 
                'range': -100, 'sweep_number': -1, 'sweep_angle': -1, 
                'IsVolumeScan': -1, 'fileLoccation': filename
                }
    try:    
        for sweep in range(0, len(angles)):
            reflectivity = radar.get_field(sweep, 'reflectivity')
            x, y, z = radar.get_gate_x_y_z(sweep, edges=False)
            x /= 1000.0
            y /= 1000.0
            z /= 1000.0
            # calculate (R)ange
            R = np.sqrt(x ** 2 + y ** 2)
            
            reflectivity = spyi.gaussian_filter(reflectivity, sigma = 2)
            
            # Masking
            reflectivity = np.ma.masked_where( reflectivity < 20, reflectivity)
            reflectivity = np.ma.masked_where( R > 80, reflectivity)
            reflectivity = np.ma.masked_where( z < 3, reflectivity)
            # folloring is to just take azimuth from 180 to 270
            reflectivity = np.ma.masked_where( x > 0, reflectivity)
            reflectivity = np.ma.masked_where( y > 0, reflectivity)
            max_ref_sweep = np.amax(np.amax(reflectivity))
            if max_ref_sweep > max_ref:
                max_sweepAngle = angles[sweep]
                max_index = np.unravel_index(reflectivity.argmax(), reflectivity.shape)
                max_height = z[max_index]
                max_range = R[max_index]
                max_ref = reflectivity[max_index]
                #max_rhohv = rhohv[max_index]
                sweep_number = sweep
                NorthSouth_distance = y[max_index]
                EastWest_distance = x[max_index]
                
                   
        return {'time_UTC': time_text,  'reflectivity': max_ref , 
                'NorthSouth_distance': NorthSouth_distance,
                'EastWest_distance': EastWest_distance, 
                'height': max_height, 'range': max_range, 
                'sweep_number': sweep_number, 'sweep_angle': max_sweepAngle, 
                'IsVolumeScan': 1, 'fileLoccation': filename
                }
    except:
        return {'time_UTC': time_text,  'reflectivity': -100, 
                'NorthSouth_distance': 1, 'EastWest_distance': 1, 'height': -1, 
                'range': -100, 'sweep_number': -1, 'sweep_angle': -1, 
                'IsVolumeScan': 1, 'fileLoccation': filename
                }

def max_inf_region(filename):
    radar = pyart.io.read(filename)
    max_ref = -100
    echo_top_height = 0
    angles = radar.fixed_angle['data']
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.isoformat()
    try:    
        for sweep in range(0, len(angles)):
            reflectivity = radar.get_field(sweep, 'reflectivity')
            reflectivity_echoTop = radar.get_field(sweep, 'reflectivity')
            x, y, z = radar.get_gate_x_y_z(sweep, edges=False)
            x /= 1000.0
            y /= 1000.0
            z /= 1000.0
            # calculate (R)ange
            R = np.sqrt(x ** 2 + y ** 2)
            
            reflectivity = spyi.gaussian_filter(reflectivity, sigma = 2)
            
            # Masking
            reflectivity = np.ma.masked_where( reflectivity < 20, reflectivity)
            reflectivity = np.ma.masked_where( R > 80, reflectivity)
            reflectivity = np.ma.masked_where( z < 3, reflectivity)
            
            
            # Region Masking 25/sqrt(2) grid 3x3 distance between CHIVO and CSAPR 75 km
            reflectivity = np.ma.masked_where(x < -35.35, reflectivity)
            reflectivity = np.ma.masked_where(x > -17.68, reflectivity)
            reflectivity = np.ma.masked_where(y < -35.35, reflectivity)
            reflectivity = np.ma.masked_where(y > -17.68, reflectivity)
            
            z = np.ma.masked_where( reflectivity_echoTop < 18, z)
            z = np.ma.masked_where(x < -35.35, z)
            z = np.ma.masked_where(x > -17.68, z)
            z = np.ma.masked_where(y < -35.35, z)
            z = np.ma.masked_where(y > -17.68, z)
            
            
            max_ref_sweep = np.amax(np.amax(reflectivity))
            echo_top_sweep = np.amax(np.amax(z))
            
            if max_ref_sweep > max_ref:
                max_sweepAngle = angles[sweep]
                max_index = np.unravel_index(reflectivity.argmax(), reflectivity.shape)
                max_height = z[max_index]
                max_range = R[max_index]
                max_ref = reflectivity[max_index]
                #max_rhohv = rhohv[max_index]
                sweep_number = sweep
            if echo_top_sweep > echo_top_height:
                echo_top_height = echo_top_sweep
        
        print(time_text + ' OK')
                       
        return {'time_UTC': time_text,  'reflectivity': max_ref , 
                'height_maxRef': max_height, 'echo_top': echo_top_height, 
                'range': max_range, 'sweep_number': sweep_number, 
                'sweep_angle': max_sweepAngle, 'fileLoccation': filename
                }
    except:
        print(time_text + ' Error')
        return {'time_UTC': time_text,  'reflectivity': -100, 
                'height_maxRef': -1, 'echo_top': -1, 
                'range': -1, 'sweep_number': -1, 
                'sweep_angle': -1, 'fileLoccation': filename
                }
#-------------------------------------------------
#This is for testing

#filename = '/net/denali/storage/radar/RELAMPAGO/quality_controlled_data/level_1a/2018/11/30/chivo.1a.20181130_033058.REL_PNL360A.nc'
#results  = max_inf_region(filename)    
    
    
    
