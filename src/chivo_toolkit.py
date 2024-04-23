#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 15:51:40 2019

@author: idariash
"""
import pandas as pd
import netCDF4
import pyart

def get_sigmetFile(filename):
    radar = pyart.io.read(filename)
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = time_start.isoformat()
    time_info = pd.read_excel (r'/net/denali/storage2/radar2/tmp/Ivan/CSU/RELAMPAGO/analysis/Climatology/sigmet_time_info.xlsx')
    closest_timeIndex = time_info.iloc[(time_info['time']- time_text).abs().argsort()[:2]]
    sigmetFile = time_info.loc[time_info['time'] == time_text]['filename']
    return sigmetFile

#--------------- for test ------------------

filename = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Attenuation/NetCDF/Tallest_Storms/DROPS/cfrad.20190125_211012.947_to_20190125_211704.658_col-radar_PPINEARLT360_SUR.nc'
sigmetFile = get_sigmetFile(filename)

df.iloc[(df['num']-input).abs().argsort()[:2]]

a = time_start + 1