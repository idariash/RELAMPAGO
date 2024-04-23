#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 10:21:42 2019

@author: idariash

This file contains fuctions to perform Level 1a of correction for data quality 
of CHIVO during RELAMPAGO campaign
"""

import pyart
import netCDF4
import pandas
from datetime import datetime
from datetime import timedelta

def Zdr_correction(radar, Zdr_bias):
    radar.fields['differential_reflectivity']['data'] = radar.fields['differential_reflectivity']['data'] - Zdr_bias
    
def azimuth_correction(radar, azimuth_error):
    radar.azimuth['data'] = radar.azimuth['data'] - azimuth_error
    if radar.scan_type == 'rhi':
        radar.fixed_angle['data'] = radar.fixed_angle['data'] - azimuth_error
    
def exportFile(radar, output_path):
    radar.metadata['title'] = 'RELAMPAGO CSU-CHIVO Level 1b data'
    radar.metadata['institution'] = 'Colorado State University'
    radar.metadata['references'] = 'Data quality control of CHIVO report'
    radar.metadata['comment'] = 'Data corrected using disdrometer coefficient from the field'
    radar.metadata['instrument_name'] = 'CSU-CHIVO'
    radar.metadata['history'] = 'Created by idariash from RAW sigmet data, June 2020'
    radar.metadata['source'] = 'Retrieved from sigmet original file'
    
    try:
        del radar.metadata['_NCProperties']
    except:
        print(' ')
        
    del radar.metadata['history_of_appended_files']
    del radar.metadata['NCO']
    
    scan_strategy = str(radar.metadata['scan_name'] )
    scan_strategy = scan_strategy.replace(' ', '');
    scan_strategy = scan_strategy.replace('b', '');
    scan_strategy = scan_strategy.replace("'", "");    
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.strftime('%Y%m%d_%H%M%S')
    #date_directory = time_start.strftime('%Y/%m/%d')
    #filename = output_path + date_directory + '/chivo.1b.' + time_text + '.' + scan_strategy + '.nc'
    filename = output_path  + '/chivo.1b.' + time_text + '.' + scan_strategy + '.nc'
    pyart.io.write_cfradial(filename, radar, time_reference=True, arm_time_variables=False)
    #print(filename)
    
def remove_field(radar, field_name):
    if field_name in radar.fields:
        del radar.fields[field_name]
    
def flip_phidp (radar, degreesToFlip):
    radar.fields['differential_phase']['data'] = radar.fields['differential_phase']['data'] - degreesToFlip
    
def find_correctionAzimuth(filename, azimuthTable_path):
    
    azimuthTable = pandas.read_excel(azimuthTable_path)
    azimuthTable_date = [datetime.strptime(x, "%d-%b-%Y %H:%M:%S'") for x in azimuthTable['Date']]
    radar = pyart.io.read(filename)
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    closest_correctionTime = 999999999; # Closest time to the scan

    for i in range(0,len(azimuthTable_date)):
        deltaTime = abs(timedelta.total_seconds(azimuthTable_date[i] - time_start))
        if deltaTime < closest_correctionTime:
            index_forAzimuthError = i
            closest_correctionTime = deltaTime
            
    return azimuthTable['AzimuthError'][index_forAzimuthError]
    
    
def correctData(filename, output_path, azimuthTable_path):
    #try:
        azimuth_error = find_correctionAzimuth(filename, azimuthTable_path)
        radar = pyart.io.read(filename)
        remove_field(radar, 'RainRateCode')
        #azimuth_correction(radar, azimuth_error)
        exportFile(radar, output_path)
        print(azimuth_error)
        return 0
    #except:
        #return 1       

