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
        echo_top = climatology.echo_top_drops(file)
        max_ref = climatology.max_reflectivity_drops(file)
        max_Kdp = climatology.max_Kdp_drops(file)
        max_Zdr = climatology.max_Zdr_drops(file)   
        min_Zdr = climatology.min_Zdr_drops(file)
    except:
        echo_top = -1
        max_ref = -100
        max_Kdp = -100
        max_Zdr = -100   
        min_Zdr = 100
    print(time_UTC)
    return (time_UTC, echo_top, max_ref, max_Kdp, max_Zdr, min_Zdr)
    
pool = mp.Pool(30)
path = '/net/denali/storage/radar/RELAMPAGO/DROPS/'
files = [f for f in glob.glob(path + "**/*RHI*.nc", recursive=True)]
files.sort()

results = pool.map_async(radar_statistics, files).get()
pool.close()
pool.join()

resutls = sorted(results, key = lambda x:x[0])

relampago_statistics = {'time_UTC':[],'echo_top':[], 
                        'max_ref': [], 'max_Kdp': [], 
                        'max_Zdr': [], 'min_Zdr':[]}

results = sorted(results, key = lambda x:x[0])

for i in range(0, len(results)):
    relampago_statistics['time_UTC'].append(results[i][0])
    relampago_statistics['echo_top'].append(results[i][1])
    relampago_statistics['max_ref'].append(results[i][2])
    relampago_statistics['max_Kdp'].append(results[i][3])
    relampago_statistics['max_Zdr'].append(results[i][4])
    relampago_statistics['min_Zdr'].append(results[i][5])


df = DataFrame(relampago_statistics, 
               columns= ['time_UTC', 'echo_top', 'max_ref',
                         'max_Kdp', 'max_Zdr', 'min_Zdr'])

export_excel = df.to_excel (r'./relampago_statistics_all_Drops.xlsx', index = None, header=True) #Don't forget to add '.xlsx' at the end of the path
