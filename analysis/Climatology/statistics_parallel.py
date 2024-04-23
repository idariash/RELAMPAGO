#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 18:46:55 2019

@author: idariash
"""

import sys
sys.path.append('/net/nasstore/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/src')

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
import gc

def radar_statistics(file):
    radar = pyart.io.read(file)
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_UTC = time_start.strftime('%Y-%m-%dT%H:%M:%S')
    try:
        #max_Kdp = climatology.max_Kdp_drops_filterCloudBoundary(file)
        [echo_top, x_echo_top, y_echo_top] = climatology.echo_top_distribution(file)
        #max_reflectivity = climatology.max_reflectivity_drops(file)
    except:
        echo_top = - 100
        x_echo_top = - 1000
        y_echo_top = - 1000
        print('error: ' + time_UTC)
    del(radar)
    gc.collect()
    #print(time_UTC)
    return (time_UTC, echo_top, x_echo_top, y_echo_top)
    
pool = mp.Pool(40)

path = '/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/level_1b.2/'
'/net/k2/storage/projects/RELAMPAGO/RMA1/quality_controlled/level_1b'
'/net/k2/storage/projects/CSAPR/level_1b.2'
'/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/level_1b.2/'
files = [f for f in glob.glob(path + "**/*.nc", recursive=True)]
files.sort()

#climatology.echo_top_distribution(files[5])


results = pool.map_async(radar_statistics, files).get()
pool.close()
pool.join()

resutls = sorted(results, key = lambda x:x[0])

relampago_statistics = {'time_UTC':[], 'echo_top': [], 'x_echo_top': [], 'y_echo_top': []}

results = sorted(results, key = lambda x:x[0])

for i in range(0, len(results)):
    if results[i][1] < 3:
        continue
    
    relampago_statistics['time_UTC'].append(results[i][0])
    relampago_statistics['echo_top'].append(results[i][1])
    relampago_statistics['x_echo_top'].append(results[i][2])
    relampago_statistics['y_echo_top'].append(results[i][3])


df = DataFrame(relampago_statistics, 
                columns= ['time_UTC', 'echo_top','x_echo_top', 'y_echo_top'])

export_excel = df.to_excel (r'/net/nasstore/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Climatology/chivo_echoTop_distribution_unbias.xlsx', 
    index = None, header=True) #Don't forget to add '.xlsx' at the end of the path
