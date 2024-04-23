#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 16:15:31 2021

@author: idariash
"""

import scipy.ndimage as spyi
import pyart
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
import matplotlib.colors as colors
from scipy.interpolate import RegularGridInterpolator


def cappi_strom_contour(filename_chivo_1a, filename_chivo_1b, 
                        filename_csapr_1a, filename_csapr_1b, ind):
    
        
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
    #min_att = np.ma.masked_where( min_att > 3, min_att)    
    
    fig = plt.figure(figsize=(6, 5))
    x, y = np.meshgrid(0.001*grid_chivo_1a.x['data'], 0.001*grid_chivo_1a.y['data'])
    
    field = 'Minimun attenuation'
    units = 'dB' 
    
    cs = plt.pcolormesh(x, y, min_att[ind], vmin=0, vmax=15, cmap='pyart_NWSRef')
    
    
    combined_field_var = np.divide((np.multiply(att_csapr, ref_chivo_1b) + 
                                        np.multiply(att_chivo, ref_csapr_1b)), 
                                       att_chivo + att_csapr)
        
    combined_field_var = np.ma.masked_where( combined_field_var < 0, combined_field_var)
     
      
    plt.contour(x, y, combined_field_var[ind], levels=[30], colors=['k', 'k', 'k'], antialised = True)
    
     #CHIVO mark
    plt.plot(0, 0, 'ok', markersize=5)
    plt.text(0, 0, ' CHIVO', fontsize=12)
    # CSAPR mark
    plt.plot(-54, -54, 'ok', markersize=5)
    plt.text(-54, -54, ' CSAPR', fontsize=12)
    # # RMA mark
    # plt.plot(-10, 20, 'ok', markersize=5)
    # plt.text(-10, 20, ' R', fontsize=12)
    
    plt.xlim(-100, 50)
    plt.ylim(-100, 50)
    
    time_start = netCDF4.num2date(grid_chivo_1a.time['data'][0], grid_chivo_1a.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M UTC')
    plt.colorbar(cs, label= 'Attenuation '+ ' (' + units + ')')
    plt.title('CHIVO '+ field + ' at '+ str(ind + 1) +' km \n'+ time_text)
    plt.xlabel('Distance east of CHIVO (km)')
    plt.ylabel('Distance north of CHIVO (km)')
    
    return fig


def dual_radar_cappi(radar1_filename_1a, radar1_filename_1b, 
                     radar2_filename_1a, radar2_filename_1b, field, ind):
    
    radar1_grid_1a = pyart.io.read_grid(radar1_filename_1a)
    radar1_grid_1b = pyart.io.read_grid(radar1_filename_1b)
    radar2_grid_1a = pyart.io.read_grid(radar2_filename_1a)
    radar2_grid_1b = pyart.io.read_grid(radar2_filename_1b)
    
    fig = plt.figure(figsize=(6, 5))
    x, y = np.meshgrid(0.001*radar1_grid_1a.x['data'], 0.001*radar1_grid_1a.y['data'])
    
    Z = radar1_grid_1a.fields['reflectivity']['data']
    
    if ('REF' in field.upper()) and (not 'DIF' in field.upper()):
        field = 'Reflectivity'
        units = 'dBZ' 
        radar1_field_var_1a = radar1_grid_1a.fields['reflectivity']['data']
        radar1_field_var_1a = np.ma.masked_where( Z < 0, radar1_field_var_1a)
        radar1_field_var_1b = radar1_grid_1b.fields['corrected_reflectivity']['data']
        radar1_att = radar1_field_var_1b[ind] - radar1_field_var_1a 
        
        radar2_field_var_1a = radar2_grid_1a.fields['reflectivity']['data']
        radar2_field_var_1b = radar2_grid_1b.fields['corrected_reflectivity']['data']
        radar2_att = radar2_field_var_1b[ind] - radar2_field_var_1a        
        
        combined_field_var = np.divide((np.multiply(1./radar1_att, radar1_field_var_1b) + 
                                        np.multiply(1./radar2_att, radar2_field_var_1b)), 
                                       1./radar1_att + 1./radar2_att)
        
        combined_field_var = np.ma.masked_where( combined_field_var < 0, combined_field_var)
         
        cs = plt.pcolormesh(x, y, combined_field_var[ind] , vmin=0, vmax=70, cmap='pyart_NWSRef')
          
        plt.contour(x, y, combined_field_var[ind], levels=[30], colors=['k', 'k', 'k'], antialised = True)
        
        
    elif 'REF' in field.upper() and 'DIF' in field.upper():
        field = 'Diff. Reflectivity'
        units = 'dB'        
        radar1_field_var_1a = radar1_grid_1a.fields['differential_reflectivity']['data']
        radar1_field_var_1a = np.ma.masked_where( Z < 5, radar1_field_var_1a)
        radar1_field_var_1b = radar1_grid_1b.fields['corrected_differential_reflectivity']['data']
        radar1_att = radar1_field_var_1b[ind] - radar1_field_var_1a 
        
        radar2_field_var_1a = radar2_grid_1a.fields['differential_reflectivity']['data']
        radar2_field_var_1b = radar2_grid_1b.fields['corrected_differential_reflectivity']['data']
        radar2_att = radar2_field_var_1b[ind] - radar2_field_var_1a        
        
        combined_field_var = np.divide((np.multiply(1./radar1_att, radar1_field_var_1b) + 
                                        np.multiply(1./radar2_att, radar2_field_var_1b)), 
                                       1./radar1_att + 1./radar2_att)
         
        cs = plt.pcolormesh(x, y, combined_field_var[ind], vmin=-2, vmax=6, cmap=pyart.graph.cm.RefDiff)
    
    
     #CHIVO mark
    plt.plot(0, 0, 'ok', markersize=5)
    plt.text(0, 0, ' CHIVO', fontsize=12)
    # CSAPR mark
    plt.plot(-54, -54, 'ok', markersize=5)
    plt.text(-54, -54, ' CSAPR', fontsize=12)
    # # RMA mark
    # plt.plot(-10, 20, 'ok', markersize=5)
    # plt.text(-10, 20, ' R', fontsize=12)
    
    plt.xlim(-100, 50)
    plt.ylim(-100, 50)
    
    time_start = netCDF4.num2date(radar1_grid_1a.time['data'][0], radar1_grid_1a.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M UTC')
    plt.colorbar(cs, label= field + ' (' + units + ')')
    plt.title('Combined '+ field + ' at '+ str(ind + 1) +' km \n'+ time_text)
    plt.xlabel('Distance east of CHIVO (km)')
    plt.ylabel('Distance north of CHIVO (km)')
    
    return fig

def cappi_att(filename_1a, filename_1b, field, ind):
    
    grid_1a = pyart.io.read_grid(filename_1a)
    grid_1b = pyart.io.read_grid(filename_1b)
    fig = plt.figure(figsize=(6, 5))
    x, y = np.meshgrid(0.001*grid_1a.x['data'], 0.001*grid_1a.y['data'])
    
    Z = grid_1a.fields['reflectivity']['data']
    
    if ('REF' in field.upper()) and (not 'DIF' in field.upper()):
        field = 'Reflectivity'
        units = 'dB' 
        field_var_1a = grid_1a.fields['reflectivity']['data']
        field_var_1a = np.ma.masked_where( Z < 0, field_var_1a)
        field_var_1b = grid_1b.fields['corrected_reflectivity']['data']
        field_var_att = field_var_1b[ind] - field_var_1a[ind]
         
        cs = plt.pcolormesh(x, y, field_var_att, vmin=0, vmax=30, cmap='pyart_NWSRef')
          
    elif 'REF' in field.upper() and 'DIF' in field.upper():
        field = 'Diff. Reflectivity'
        units = 'dB'        
        field_var_1a = grid_1a.fields['differential_reflectivity']['data']
        field_var_1a = np.ma.masked_where( Z < 5, field_var_1a)
        field_var_1b = grid_1b.fields['corrected_differential_reflectivity']['data']
        field_var_att = field_var_1b[ind] - field_var_1a[ind]
         
        cs = plt.pcolormesh(x, y, field_var_att, vmin=0, vmax=10, cmap='pyart_NWSRef')
    
    
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
    
    time_start = netCDF4.num2date(grid_1a.time['data'][0], grid_1a.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M UTC')
    plt.colorbar(cs, label= 'Attenuation '+ ' (' + units + ')')
    plt.title('CHIVO '+ field + ' Attenuation at '+ str(ind + 1) +' km \n'+ time_text)
    plt.xlabel('Distance east of CHIVO (km)')
    plt.ylabel('Distance north of CHIVO (km)')
    
    return fig


def cappi_1a(filename_1b, field, ind, radar_name):

    pyart_grid_1b = pyart.io.read_grid(filename_1b)
    fig = plt.figure(figsize=(6, 5))
    Z = pyart_grid_1b.fields['reflectivity']['data']
    x, y = np.meshgrid(0.001*pyart_grid_1b.x['data'], 0.001*pyart_grid_1b.y['data']) 
    
    if ('REF' in field.upper()) and (not 'DIF' in field.upper()):
        field = 'Reflectivity'
        units = 'dBZ'        
        field_var = pyart_grid_1b.fields['reflectivity']['data']
        field_var = np.ma.masked_where( Z < 0, field_var)        
        cs = plt.pcolormesh(x, y, field_var[ind], vmin=0, vmax=70, cmap='pyart_NWSRef')
          
    elif 'REF' in field.upper() and 'DIF' in field.upper():
        field = 'Diff. Reflectivity'
        units = 'dB'        
        field_var = pyart_grid_1b.fields['differential_reflectivity']['data']       
        field_var = np.ma.masked_where( Z < 10, field_var)        
        cs = plt.pcolormesh(x, y, field_var[ind], vmin=-2, vmax=6, cmap=pyart.graph.cm.RefDiff)
        
    elif 'SPECIFIC' in field.upper():
        field = 'Specific Diff. Phase'
        units = 'deg/km'        
        field_var = pyart_grid_1b.fields['specific_differential_phase']['data']      
        field_var = np.ma.masked_where( Z < 10, field_var)        
        cs = plt.pcolormesh(x, y, field_var[ind], vmin=0, vmax=6, cmap=pyart.graph.cm.Theodore16)
        
    elif 'CORRELATION' in field.upper():
        field = 'Copolar Corr. Ratio'
        units = ''        
        field_var = pyart_grid_1b.fields['cross_correlation_ratio']['data']       
        field_var = np.ma.masked_where( Z < 10, field_var)        
        cs = plt.pcolormesh(x, y, field_var[ind], vmin=0.5, vmax=1, cmap='jet')
        
    #CHIVO mark
    plt.plot(0, 0, 'ok', markersize=5)
    plt.text(0, 0, ' C', fontsize=12)
    # CSAPR mark
    plt.plot(-54, -54, 'ok', markersize=5)
    plt.text(-54, -54, ' Y', fontsize=12)
    # RMA mark
    plt.plot(-10, 20, 'ok', markersize=5)
    plt.text(-10, 20, ' R', fontsize=12)
    
    plt.xlim(-90, 40)
    plt.ylim(-90, 40)
    
    
    time_start = netCDF4.num2date(pyart_grid_1b.time['data'][0], pyart_grid_1b.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M UTC')
    plt.colorbar(cs, label= field + ' (' + units + ')')
    plt.title(radar_name + ' ' + field + ' at '+ str(ind + 1) +' km \n'+ time_text)
    plt.xlabel('Distance east of CHIVO (km)')
    plt.ylabel('Distance north of CHIVO (km)')
    return fig



def cappi(filename_1b, field, ind):

    pyart_grid_1b = pyart.io.read_grid(filename_1b)
    fig = plt.figure(figsize=(6, 5))
    Z = pyart_grid_1b.fields['corrected_reflectivity']['data']
    x, y = np.meshgrid(0.001*pyart_grid_1b.x['data'], 0.001*pyart_grid_1b.y['data']) 
    
    if ('REF' in field.upper()) and (not 'DIF' in field.upper()):
        field = 'Reflectivity'
        units = 'dBZ'        
        field_var = pyart_grid_1b.fields['corrected_reflectivity']['data']
        field_var = np.ma.masked_where( Z < 10, field_var)        
        cs = plt.pcolormesh(x, y, field_var[ind], vmin=0, vmax=70, cmap='pyart_NWSRef')
          
    elif 'REF' in field.upper() and 'DIF' in field.upper():
        field = 'Diff. Reflectivity'
        units = 'dB'        
        field_var = pyart_grid_1b.fields['corrected_differential_reflectivity']['data']       
        field_var = np.ma.masked_where( Z < 10, field_var)        
        cs = plt.pcolormesh(x, y, field_var[ind], vmin=-2, vmax=6, cmap=pyart.graph.cm.RefDiff)
        
    elif 'SPECIFIC' in field.upper():
        field = 'Specific Diff. Phase'
        units = 'deg/km'        
        field_var = pyart_grid_1b.fields['corrected_specific_differential_phase']['data']      
        field_var = np.ma.masked_where( Z < 10, field_var)        
        cs = plt.pcolormesh(x, y, field_var[ind], vmin=0, vmax=6, cmap=pyart.graph.cm.Theodore16)
        
    elif 'CORRELATION' in field.upper():
        field = 'Copolar Corr. Ratio'
        units = ''        
        field_var = pyart_grid_1b.fields['corrected_copol_correlation_ratio']['data']       
        field_var = np.ma.masked_where( Z < 10, field_var)        
        cs = plt.pcolormesh(x, y, field_var[ind], vmin=0.5, vmax=1, cmap='jet')
        
    #CHIVO mark
    plt.plot(0, 0, 'ok', markersize=5)
    plt.text(0, 0, ' C', fontsize=12)
    # CSAPR mark
    plt.plot(-54, -54, 'ok', markersize=5)
    plt.text(-54, -54, ' Y', fontsize=12)
    # RMA mark
    plt.plot(-10, 20, 'ok', markersize=5)
    plt.text(-10, 20, ' R', fontsize=12)
    
    plt.xlim(-90, 40)
    plt.ylim(-90, 40)
    
    
    time_start = netCDF4.num2date(pyart_grid_1b.time['data'][0], pyart_grid_1b.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M UTC')
    plt.colorbar(cs, label= field + ' (' + units + ')')
    plt.title('CSAPR '+ field + ' at '+ str(ind + 1) +' km \n'+ time_text)
    plt.xlabel('Distance east of CHIVO (km)')
    plt.ylabel('Distance north of CHIVO (km)')
    return fig

def cappi_grid(grid, field, ind):

    pyart_grid_1b = grid
    fig = plt.figure(figsize=(6, 5))
    Z = pyart_grid_1b.fields['corrected_reflectivity']['data']
    x, y = np.meshgrid(0.001*pyart_grid_1b.x['data'], 0.001*pyart_grid_1b.y['data']) 
    
    if ('REF' in field.upper()) and (not 'DIF' in field.upper()):
        field = 'Reflectivity'
        units = 'dBZ'        
        field_var = pyart_grid_1b.fields['corrected_reflectivity']['data']
        field_var = np.ma.masked_where( Z < 10, field_var)        
        cs = plt.pcolormesh(x, y, field_var[ind], vmin=0, vmax=70, cmap='pyart_NWSRef')
          
    elif 'REF' in field.upper() and 'DIF' in field.upper():
        field = 'Diff. Reflectivity'
        units = 'dB'        
        field_var = pyart_grid_1b.fields['corrected_differential_reflectivity']['data']       
        field_var = np.ma.masked_where( Z < 10, field_var)        
        cs = plt.pcolormesh(x, y, field_var[ind], vmin=-2, vmax=6, cmap=pyart.graph.cm.RefDiff)
        
    elif 'SPECIFIC' in field.upper():
        field = 'Specific Diff. Phase'
        units = 'deg/km'        
        field_var = pyart_grid_1b.fields['corrected_specific_differential_phase']['data']      
        field_var = np.ma.masked_where( Z < 10, field_var)        
        cs = plt.pcolormesh(x, y, field_var[ind], vmin=0, vmax=6, cmap=pyart.graph.cm.Theodore16)
        
    elif 'CORRELATION' in field.upper():
        field = 'Copolar Corr. Ratio'
        units = ''        
        field_var = pyart_grid_1b.fields['corrected_copol_correlation_ratio']['data']       
        field_var = np.ma.masked_where( Z < 10, field_var)        
        cs = plt.pcolormesh(x, y, field_var[ind], vmin=0.5, vmax=1, cmap='jet')
        
    elif 'RAIN' in field.upper():
        field = 'Rain Rate'
        units = 'mm/h'        
        field_var = pyart_grid_1b.fields['RainRate']['data']       
        field_var = np.ma.masked_where( Z < 10, field_var)        
        cs = plt.pcolormesh(x, y, field_var[ind], vmin=1, vmax=100,
                           cmap=pyart.graph.cm.LangRainbow12, norm=colors.LogNorm())
        
    elif 'CLASS' in field.upper():
        field = 'Hydrometoer Classification'
        units = 'mm/h'        
        field_var = pyart_grid_1b.fields['HydroClass']['data']       
        field_var = np.ma.masked_where( Z < 10, field_var)
        # field_var = np.ma.masked_where( field_var < 10, field_var)
        # field_var = np.ma.masked_where( field_var > 12, field_var)
        cs = plt.pcolormesh(x, y, field_var[ind], vmin=9, vmax=13, cmap='jet')
        
    #CHIVO mark
    plt.plot(0, 0, 'ok', markersize=5)
    plt.text(0, 0, ' C', fontsize=12)
    # CSAPR mark
    plt.plot(-54, -54, 'ok', markersize=5)
    plt.text(-54, -54, ' Y', fontsize=12)
    # RMA mark
    plt.plot(-10, 20, 'ok', markersize=5)
    plt.text(-10, 20, ' R', fontsize=12)
    
    plt.xlim(-90, 40)
    plt.ylim(-90, 40)
    
    
    time_start = netCDF4.num2date(pyart_grid_1b.time['data'][0], pyart_grid_1b.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M UTC')
    plt.colorbar(cs, label= field + ' (' + units + ')')
    plt.title('CSAPR '+ field + ' at '+ str(ind + 1) +' km \n'+ time_text)
    plt.xlabel('Distance east of CHIVO (km)')
    plt.ylabel('Distance north of CHIVO (km)')
    return fig


def chivo_csapr_cappi(filename_chivo, filename_csapr, field, ind):
    
    grid_chivo = pyart.io.read_grid(filename_chivo)
    grid_csapr = pyart.io.read_grid(filename_csapr)
    fig = plt.figure(figsize=(6, 5))
    x, y = np.meshgrid(0.001*grid_chivo.x['data'], 0.001*grid_chivo.y['data'])
    
    
    if ('REF' in field.upper()) and (not 'DIF' in field.upper()):
        field = 'Reflectivity'
        units = 'dBZ' 
        field_var_chivo = grid_chivo.fields['corrected_reflectivity']['data']
        field_var_csapr = grid_csapr.fields['corrected_reflectivity']['data']
        field_var_diff = field_var_chivo[ind] - field_var_csapr[ind]
         
        cs = plt.pcolormesh(x, y, field_var_diff, vmin=-10, vmax=10, cmap='pyart_NWSRef')
          
    elif 'REF' in field.upper() and 'DIF' in field.upper():
        field = 'Diff. Reflectivity'
        units = 'dB'        
        field_var_chivo = grid_chivo.fields['corrected_differential_reflectivity']['data']
        field_var_csapr = grid_csapr.fields['corrected_differential_reflectivity']['data']
        field_var_diff = field_var_chivo[ind] - field_var_csapr[ind]
         
        cs = plt.pcolormesh(x, y, field_var_diff, vmin=-2, vmax=2, cmap='pyart_NWSRef')
    
    
    
    #CHIVO mark
    plt.plot(0, 0, 'ok', markersize=5)
    plt.text(0, 0, ' C', fontsize=12)
    # CSAPR mark
    plt.plot(-54, -54, 'ok', markersize=5)
    plt.text(-54, -54, ' Y', fontsize=12)
    # RMA mark
    plt.plot(-10, 20, 'ok', markersize=5)
    plt.text(-10, 20, ' R', fontsize=12)
    
    plt.xlim(-90, 40)
    plt.ylim(-90, 40)
    time_start = netCDF4.num2date(grid_chivo.time['data'][0], grid_chivo.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M UTC')
    plt.colorbar(cs, label= field + ' (' + units + ')')
    plt.title('CSU-CHIVO '+ field + ' at '+ str(ind + 1) +' km \n'+ time_text)
    plt.xlabel('Distance east of CHIVO (km)')
    plt.ylabel('Distance north of CHIVO (km)')
    
    return fig

def chivo_csapr_cappi_1a(filename_chivo, filename_csapr, field, ind):
    
    grid_chivo = pyart.io.read_grid(filename_chivo)
    grid_csapr = pyart.io.read_grid(filename_csapr)
    fig = plt.figure(figsize=(6, 5))
    x, y = np.meshgrid(0.001*grid_chivo.x['data'], 0.001*grid_chivo.y['data'])
    
    
    if ('REF' in field.upper()) and (not 'DIF' in field.upper()):
        field = 'Reflectivity'
        units = 'dBZ' 
        field_var_chivo = grid_chivo.fields['reflectivity']['data']
        field_var_csapr = grid_csapr.fields['reflectivity']['data']
        field_var_diff = field_var_chivo[ind] - field_var_csapr[ind]
         
        cs = plt.pcolormesh(x, y, field_var_diff, vmin=-10, vmax=10, cmap='pyart_NWSRef')
          
    elif 'REF' in field.upper() and 'DIF' in field.upper():
        field = 'Diff. Reflectivity'
        units = 'dB'        
        field_var_chivo = grid_chivo.fields['differential_reflectivity']['data']
        field_var_csapr = grid_csapr.fields['differential_reflectivity']['data']
        field_var_diff = field_var_chivo[ind] - field_var_csapr[ind]
         
        cs = plt.pcolormesh(x, y, field_var_diff, vmin=-2, vmax=2, cmap='pyart_NWSRef')
    
    
    
    #CHIVO mark
    plt.plot(0, 0, 'ok', markersize=5)
    plt.text(0, 0, ' C', fontsize=12)
    # CSAPR mark
    plt.plot(-54, -54, 'ok', markersize=5)
    plt.text(-54, -54, ' Y', fontsize=12)
    # RMA mark
    plt.plot(-10, 20, 'ok', markersize=5)
    plt.text(-10, 20, ' R', fontsize=12)
    
    plt.xlim(-90, 40)
    plt.ylim(-90, 40)
    time_start = netCDF4.num2date(grid_chivo.time['data'][0], grid_chivo.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M UTC')
    plt.colorbar(cs, label= field + ' (' + units + ')')
    plt.title('CHIVO - CSAPR Original '+ field + ' at '+ str(ind + 1) +' km \n'+ time_text)
    plt.xlabel('Distance east of CHIVO (km)')
    plt.ylabel('Distance north of CHIVO (km)')
    
    return fig

def chivo_csapr_cappi_1b(filename_chivo, filename_csapr, field, ind):
    
    grid_chivo = pyart.io.read_grid(filename_chivo)
    grid_csapr = pyart.io.read_grid(filename_csapr)
    fig = plt.figure(figsize=(6, 5))
    x, y = np.meshgrid(0.001*grid_chivo.x['data'], 0.001*grid_chivo.y['data'])
    
    
    if ('REF' in field.upper()) and (not 'DIF' in field.upper()):
        field = 'Reflectivity'
        units = 'dBZ' 
        field_var_chivo = grid_chivo.fields['corrected_reflectivity']['data']
        field_var_csapr = grid_csapr.fields['corrected_reflectivity']['data']
        field_var_diff = field_var_chivo[ind] - field_var_csapr[ind]
         
        cs = plt.pcolormesh(x, y, field_var_diff, vmin=-10, vmax=10, cmap='pyart_NWSRef')
          
    elif 'REF' in field.upper() and 'DIF' in field.upper():
        field = 'Diff. Reflectivity'
        units = 'dB'        
        field_var_chivo = grid_chivo.fields['corrected_differential_reflectivity']['data']
        field_var_csapr = grid_csapr.fields['corrected_differential_reflectivity']['data']
        field_var_diff = field_var_chivo[ind] - field_var_csapr[ind]
         
        cs = plt.pcolormesh(x, y, field_var_diff, vmin=-2, vmax=2, cmap='pyart_NWSRef')
    
    
    
    #CHIVO mark
    plt.plot(0, 0, 'ok', markersize=5)
    plt.text(0, 0, ' C', fontsize=12)
    # CSAPR mark
    plt.plot(-54, -54, 'ok', markersize=5)
    plt.text(-54, -54, ' Y', fontsize=12)
    # RMA mark
    plt.plot(-10, 20, 'ok', markersize=5)
    plt.text(-10, 20, ' R', fontsize=12)
    
    plt.xlim(-90, 40)
    plt.ylim(-90, 40)
    time_start = netCDF4.num2date(grid_chivo.time['data'][0], grid_chivo.time['units'])
    time_text = time_start.strftime('%Y-%m-%dT%H:%M UTC')
    plt.colorbar(cs, label= field + ' (' + units + ')')
    plt.title('CHIVO - CSAPR Corrected '+ field + ' at '+ str(ind + 1) +' km \n'+ time_text)
    plt.xlabel('Distance east of CHIVO (km)')
    plt.ylabel('Distance north of CHIVO (km)')
    
    return fig
