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
    try:
        radar = pyart.io.read(file)
        time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
        time_UTC = time_start.strftime('%Y-%m-%dT%H:%M:%S')
        #max_Kdp = climatology.max_Kdp_drops_filterCloudBoundary(file)
        echo_top = climatology.echo_top_csapr(file)
        max_reflectivity = climatology.max_reflectivity_csapr(file)
        del(radar)
        print(time_UTC)
    except:
        time_UTC = '2000-01-01T00:00:00'
        echo_top = - 100
        max_reflectivity = -100
        print('error: ' + time_UTC)
    return (time_UTC, echo_top, max_reflectivity)
    
pool = mp.Pool(4)
path = '/net/denali/storage/radar/CSAPR/level_1a/'
files = [f for f in glob.glob(path + "**/*201811*.nc", recursive=True)]
files.sort()

results = pool.map_async(radar_statistics, files).get()
pool.close()
pool.join()

resutls = sorted(results, key = lambda x:x[0])

relampago_statistics = {'time_UTC':[], 'echo_top': [], 'max_ref': []}

results = sorted(results, key = lambda x:x[0])

for i in range(0, len(results)):
    relampago_statistics['time_UTC'].append(results[i][0])
    relampago_statistics['echo_top'].append(results[i][1])
    relampago_statistics['max_ref'].append(results[i][2])


df = DataFrame(relampago_statistics, 
               columns= ['time_UTC', 'echo_top','max_ref'])

export_excel = df.to_excel (r'./relampago_statistics_csapr_nov.xlsx', index = None, header=True) #Don't forget to add '.xlsx' at the end of the path


#export_excel = df.to_excel (r'./relampago_statistics_csapr_ppi.xlsx', index = None, header=True) #Don't forget to add '.xlsx' at the end of the path
