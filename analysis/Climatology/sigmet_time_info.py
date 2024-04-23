#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 22:19:41 2019

@author: idariash
"""

import os
import datetime
import netCDF4


import glob
from pandas import DataFrame
import pyart
import multiprocessing as mp

def get_time_info(file):
    try:
        radar = pyart.io.read(file)
        time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
        time_text = time_start.isoformat()
        print(time_text)
        return (time_text, file)
    except:
        return ('error', file)

pool = mp.Pool(30)
path = '/net/denali/storage/radar/RELAMPAGO/RAW/'
files = [f for f in glob.glob(path + "**/*COL*", recursive=True)]
files.sort()

results = pool.map_async(get_time_info, files).get()
pool.close()
pool.join()

resutls = sorted(results, key = lambda x:x[0])

time_info = {'time':[], 'filename':[]}
for i in range(0, len(results)):
    time_info['time'].append(results[i][0])
    time_info['filename'].append(results[i][1])
    
df = DataFrame(time_info,columns= ['time', 'filename'])

export_excel = df.to_excel (r'./sigmet_time_info.xlsx', index = None, header=True) 