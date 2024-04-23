#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 17:50:55 2019

@author: idariash
This script creates a new version of CHIVO data with azimuth and Zdr correction 
"""


import chivo_dataQuality_1a as dq_1a
import pandas
import pyart
from datetime import datetime
from datetime import timedelta
import netCDF4
import os

def main(day_to_process, start_file):   #YYYYMMDD
    day_to_process = str(day_to_process)
    YYYY = day_to_process [0:4]
    MM = day_to_process [4:6]
    DD = day_to_process [6:8]
    DataPath = '/net/denali/storage2/radar2/tmp/RELAMPAGO/RAW/' + YYYY + '/' + MM + '/' + DD + '/'
    output_path = '/net/denali/storage2/radar2/tmp/RELAMPAGO/quality_controlled_data/level_1a/' + YYYY + '/' + MM + '/' + DD + '/'
    azimuthTable_path = '/net/denali/storage2/radar2/tmp/Ivan/CSU/CHIVO/src/python/data_quality/level_1a/CHIVO_AzimuthCorrectionTable_Extended.xlsx'
    Zdr_bias = 0.75
    log_name = output_path + YYYY + MM + DD + '.log'
    Directory = os.listdir(DataPath)
    log = open(log_name,'a')
    log.write('This log is for ' + YYYY+MM+DD + '\n')
    log.write(str(datetime.now()) + '\n')
    log.close()
    for x in range(start_file,len(Directory)):
        filename = DataPath + Directory[x]
        print (str(x) + '/' + str(len(Directory)) + ' ' + Directory[x])
        error = dq_1a.correctData(filename, output_path, azimuthTable_path, Zdr_bias)
        log = open(log_name,'a')
        
        if error == 0:
            log.write(str(x) + '/' + str(len(Directory)) + ' ' + Directory[x] + ' done \n')
        else:
            log.write(str(x) + '/' + str(len(Directory)) + ' ' + Directory[x] + ' error \n')
       
        log.close()
        
        

main(20190101, 0)
main(20190102, 0)
main(20190103, 0)
main(20190104, 0)
main(20190105, 0)
main(20190106, 0)
main(20190107, 0)
main(20190108, 0)
main(20190109, 0)
main(20190110, 0)
main(20190111, 0)


