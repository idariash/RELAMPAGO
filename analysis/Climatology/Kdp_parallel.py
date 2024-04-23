#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 18:46:55 2019

@author: idariash
"""

import sys
sys.path.append('/net/denali/storage2/radar2/tmp/Ivan/CSU/RELAMPAGO/src')

import radar_climatology as climatology
import os
import datetime
import netCDF4

import importlib
importlib.reload(climatology)

import glob
from pandas import DataFrame
import pyart

import multiprocessing as mp

def radar_statistics(file):
    radar = pyart.io.read(file)
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_UTC = time_start.isoformat()
    try:
        max_Kdp = climatology.max_Kdp_drops_filterCloudBoundary(file)
    except:
        max_Kdp = -100
        print('error: ' + time_UTC)
    print(time_UTC)
    return (time_UTC, max_Kdp)
    
pool = mp.Pool(10)
path = '/net/denali/storage/radar/RELAMPAGO/DROPS/'
files = [f for f in glob.glob(path + "**/*RHI*.nc", recursive=True)]
files.sort()

results = pool.map_async(radar_statistics, files).get()
pool.close()
pool.join()

resutls = sorted(results, key = lambda x:x[0])

relampago_statistics = {'time_UTC':[], 'max_Kdp': []}

results = sorted(results, key = lambda x:x[0])

for i in range(0, len(results)):
    relampago_statistics['time_UTC'].append(results[i][0])
    relampago_statistics['max_Kdp'].append(results[i][3])


df = DataFrame(relampago_statistics, 
               columns= ['time_UTC', 'max_Kdp'])

export_excel = df.to_excel (r'./relampago_statistics_Drops_Kdp.xlsx', index = None, header=True) #Don't forget to add '.xlsx' at the end of the path
