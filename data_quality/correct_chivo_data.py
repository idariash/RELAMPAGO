#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 19:32:34 2019

@author: idariash
"""
import pyart
import netCDF4
import datetime
from matplotlib import pyplot as plt

def Zdr_correction(radar, Zdr_bias):
    radar.fields['differential_reflectivity']['data'] = radar.fields['differential_reflectivity']['data'] - Zdr_bias
    
def azimuth_correction(radar, azimuth_error):
    radar.azimuth['data'] = radar.azimuth['data'] - azimuth_error
    

chivo = '/net/denali/storage2/radar2/tmp/RELAMPAGO/RAW/2019/01/25/COL20190125_192419'
#'/net/denali/storage2/radar2/tmp/RELAMPAGO/RAW/2019/01/25/COL20190125_193421'
#'/net/denali/storage2/radar2/tmp/RELAMPAGO/RAW/2019/01/25/COL20190125_193048'

CHIVO = pyart.io.read(chivo)

#Correct Azimuth and Zdr
azimuth_correction(CHIVO, 13.5)
Zdr_correction(CHIVO, 0.75)

#Correct Velocity CHIVO,I can maka a function for this

gatefilter = pyart.filters.GateFilter(CHIVO)
gatefilter.exclude_transition()
gatefilter.exclude_invalid('velocity')
gatefilter.exclude_invalid('reflectivity')
gatefilter.exclude_outside('reflectivity', 0, 80)
gatefilter.exclude_outside('normalized_coherent_power', 0.2, 1)
gatefilter.exclude_outside('cross_correlation_ratio', 0.8, 1)

nyq = CHIVO.instrument_parameters['nyquist_velocity']['data'][0]
last_radar = CHIVO
corr_vel = pyart.correct.dealias_region_based(
    last_radar, vel_field='velocity', keep_original=False, 
    gatefilter = gatefilter, nyquist_vel=nyq, centered = True)
last_radar.add_field('corrected_velocity', corr_vel, replace_existing = True)

# perform dealiasing 4DD
dealias_data = pyart.correct.dealias_fourdd(
    CHIVO, last_radar = last_radar, gatefilter=gatefilter)
CHIVO.add_field('corrected_velocity', dealias_data, replace_existing = True)

#Create new NetCDF file
output_path = '/net/denali/storage2/radar2/tmp/Ivan/CSU/CHIVO/res/data/test_data_quality'
CHIVO.metadata['title'] = 'Level 1 CHIVO Radar RELAMPAGO'
CHIVO.metadata['institution'] = 'Colorado State University'
CHIVO.metadata['reference'] = 'RELAMPAGO CHIVO report and Data quality control of CHIVO report'
CHIVO.metadata['comment'] = 'Data collected during RELAMPAGO field campaign at Lozada site in Cordoba, Argentina'
CHIVO.metadata['instrument_name'] = 'CSU-CHIVO'
CHIVO.metadata['history'] = 'Created by Ivan Arias from RAW sigmet data, July 2019'

scan_strategy = str(CHIVO.metadata['sigmet_task_name'] )
scan_strategy = scan_strategy[2 : len(scan_strategy) - 2]
time_start = netCDF4.num2date(CHIVO.time['data'][0], CHIVO.time['units'])
time_text = time_start.strftime('%Y%m%d_%H%M%S')
filename = output_path + '/chivo.1a.' + time_text + '.' + scan_strategy + '.nc'
pyart.io.write_cfradial(filename, CHIVO, format='NETCDF4', time_reference=True, arm_time_variables=False)
