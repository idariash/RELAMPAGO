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
import numpy as np

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
        [rain_accumulation, hail_counts] = climatology.accumulated_RR_Hail(file)
    except:
        rain_accumulation = - 100
        hail_counts = -100
        print('error: ' + time_UTC)
    del(radar)
    gc.collect()
    #print(time_UTC)
    return (time_UTC, rain_accumulation, hail_counts)
    

# filename = '/net/k2/storage/people/idariash/home/Field_Campaigns/Relampago/Analisis/EOL/data/chivo_ppi/chivo.1b.2.20181214_023049.REL_PNL360A.nc'
# [rain_accumulation, hail_counts] = climatology.accumulated_RR_Hail(filename)
# [time_UTC, rain_accumulation, hail_counts] = radar_statistics(filename)

path = '/net/k2/storage/projects/RELAMPAGO/RMA1/quality_controlled/level_1b'
'/net/k2/storage/projects/RELAMPAGO/RMA1/quality_controlled/level_1b'
'/net/k2/storage/projects/CSAPR/level_1b.2'
'/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/level_1b.2/'
files = [f for f in glob.glob(path + "**/*.nc", recursive=True)]
files.sort()

#climatology.echo_top_distribution(files[5])



pool = mp.Pool(40)
results = pool.map_async(radar_statistics, files).get()
pool.close()
pool.join()

resutls = sorted(results, key = lambda x:x[0])

relampago_statistics = {'time_UTC':[], 'rain_accumulation': [], 'hail_counts': [],}

results = sorted(results, key = lambda x:x[0])

for i in range(0, len(results)):
    if results[i][1] < 10 or np.ma.is_masked(results[i][1]):
        continue
    
    relampago_statistics['time_UTC'].append(results[i][0])
    relampago_statistics['rain_accumulation'].append(results[i][1])
    relampago_statistics['hail_counts'].append(results[i][2])


df = DataFrame(relampago_statistics, 
                columns= ['time_UTC', 'rain_accumulation','hail_counts',])

export_excel = df.to_excel (r'/net/nasstore/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Climatology/rma_rainRate_hail.xlsx', 
    index = None, header=True) #Don't forget to add '.xlsx' at the end of the path
