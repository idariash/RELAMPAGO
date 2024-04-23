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
import numpy as np

def radar_statistics(file):
    radar = pyart.io.read(file)
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_UTC = time_start.strftime('%Y-%m-%dT%H:%M:%S')
    try:
        echo_top = climatology.max_echoTop_grid(file)
    except:
        echo_top = - 100
        print('error: ' + time_UTC)
    del(radar)
    gc.collect()
    #print(time_UTC)
    return (time_UTC, echo_top)
    


path = '/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/level_1b.2'
'/net/k2/storage/projects/RELAMPAGO/RMA1/quality_controlled/level_1b'
'/net/k2/storage/projects/CSAPR/level_1b.2'
'/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/level_1b.2/'
files = [f for f in glob.glob(path + "/**/*.nc", recursive=True)]
files.sort()

# file = '/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/level_1b.2/2018/12/14/chivo.1b.2.20181214_023049.REL_PNL360A.nc'

# resutls = radar_statistics(file)

pool = mp.Pool(20)
results = pool.map_async(radar_statistics, files).get()
pool.close()
pool.join()
results = sorted(results, key = lambda x:x[0])
relampago_statistics = {'time_UTC':[], 'echo_top': []}



for i in range(0, len(results)):
    if results[i][1] < 6 or np.ma.is_masked(results[i][1]):
        continue
    
    relampago_statistics['time_UTC'].append(results[i][0])
    relampago_statistics['echo_top'].append(results[i][1])


df = DataFrame(relampago_statistics, 
                columns= ['time_UTC', 'echo_top',])

export_excel = df.to_excel (r'/net/nasstore/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Climatology/chivo_echoTop_grid.xlsx', 
    index = None, header=True) #Don't forget to add '.xlsx' at the end of the path
