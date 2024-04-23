#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 16:11:19 2019

@author: idariash
"""

import pandas
import pyart
from datetime import datetime
from datetime import timedelta
import netCDF4

azimuthTable = pandas.read_excel('/net/denali/storage2/radar2/tmp/Ivan/CSU/CHIVO/src/python/data_quality/level_1a/CHIVO_AzimuthCorrectionTable.xlsx')
filename = '/net/denali/storage2/radar2/tmp/RELAMPAGO/RAW/2018/11/29/COL20181129_170604'

a = azimuthTable['Date'][2]
b = azimuthTable['AzimuthError'][2]

aa = a.replace('Nov','12')
c = datetime.strptime(a, "%d-%b-%Y %H:%M:%S'")

azimuthTable_date = [datetime.strptime(x, "%d-%b-%Y %H:%M:%S'") for x in azimuthTable['Date']]

radar = pyart.io.read(filename)
time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])

d = azimuthTable_date[4] - time_start

ee = timedelta.total_seconds(d)

closest_correctionTime = 999999999; # Closest time to the scan

for i in range(0,len(azimuthTable_date)):
    deltaTime = abs(timedelta.total_seconds(azimuthTable_date[i] - time_start))
    if deltaTime < closest_correctionTime:
        index_forAzimuthError = i
        closest_correctionTime = deltaTime
        
    
print(azimuthTable['AzimuthError'][index_forAzimuthError])

#[datetime.strptime(x, '%m/%d/%Y') for x in attack_dates]
