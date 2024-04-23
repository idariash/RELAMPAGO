#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 17:50:55 2019

@author: idariash
This script creates a new version of CHIVO data with azimuth and Zdr correction 
"""


import chivo_dataQuality_1b as dq_1b
import pandas
import pyart
from datetime import datetime
from datetime import timedelta
import netCDF4
import os
import multiprocessing as mp
import glob


azimuthTable_path = '/net/denali/storage2/radar2/people/idariash/home/CSU/CHIVO/data_quality/level_1a/CHIVO_AzimuthCorrectionTable_Extended.xlsx'

output_path_drops = '/net/nasstore/students/GRAD/ECE/idariash/home/tmp/';
HydroClass_exe = '/net/denali/storage2/radar2/people/idariash/home/Utiles/DROPS2/Whitney/hydroclass'

'/net/denali/storage2/radar2/people/idariash/home/Utiles/DROPS2/DROPSX2/hydroclassX2'

'/net/denali/storage2/radar2/people/idariash/home/Utiles/DROPS2/DROPSX2/hydroclassX2'

'/net/denali/storage2/radar2/people/idariash/home/Utiles/DROPS2/Whitney/hydroclass'

#'/net/denali/storage2/radar2/people/idariash/home/Utiles/DROPS2/DROPSX2/hydroclassX2'

#'/net/denali/storage2/radar2/people/idariash/home/Utiles/DROPS2/Whitney/hydroclass'

#'/net/denali/storage2/radar2/people/idariash/home/Utiles/DROPS2/hydroclass'
#'/net/denali/storage2/radar2/people/idariash/home/Utiles/DROPS2/Whitney/hydroclass'

#'/net/denali/storage2/radar2/people/idariash/home/Utiles/DROPS2/DROPS_20200603/hydroclass'
#'/net/denali/storage2/radar2/people/idariash/home/Utiles/DROPS2/whitney_202008/hydroclass'

 
 
#'/net/denali/storage2/radar2/people/idariash/home/Utiles/DROPS2/DROPS_20200603/hydroclass'
#'/net/denali/storage2/radar2/people/idariash/home/Utiles/DROPS2/DROPS_20200603/hydroclass'
 

sounding_file = ' ';

membership_functions = '/net/denali/storage2/radar2/people/idariash/home/Utiles/DROPS2/DROPSX2/membership_functions/' 

#'/net/denali/storage2/radar2/people/idariash/home/Utiles/DROPS2/membership_functions/';

def drops_chivo(file, output_path):
    
    radar = pyart.io.read(file)
    radar_name = radar.metadata['instrument_name']
    
    if 'COL' in radar_name.upper():
        
        '/net/denali/storage2/radar2/people/idariash/home/Utiles/DROPS2/DROPSX2/RE_CHIVO_RR_update.v5.ini'
        config_file = '/net/denali/storage2/radar2/people/idariash/home/Utiles/DROPS2/RE_CHIVO_RR_update.v5.ini'
        #'/net/denali/storage2/radar2/people/idariash/home/Utiles/DROPS2/DROPSX2/RE_CHIVO_RR_update.v5_Sounak.ini' 
        #'/net/denali/storage2/radar2/people/idariash/home/Utiles/DROPS2/RE_CHIVO_RR_update.v5.ini';
        
        [input_path, filename] = os.path.split(file)
        file_drops_tmp = os.path.join(output_path_drops, filename)
        
        cmd = (HydroClass_exe + ' ' + file + ' -o ' + file_drops_tmp + 
        ' -c ' + config_file + ' -m VHS -d ' + membership_functions)
        
        os.system(cmd)
        
        cmd = 'ncks -A -v fixed_angle ' + file + ' ' + file_drops_tmp 
        os.system(cmd)
        
        cmd = 'ncks -A -v sweep_mode ' + file + ' ' + file_drops_tmp
        os.system(cmd)
        
        # cmd = 'ncks -A -v time ' + file + ' ' + file_drops_tmp
        # os.system(cmd)
        
        error = dq_1b.correctData(file_drops_tmp, output_path, azimuthTable_path)
        cmd = 'rm ' + file_drops_tmp 
        os.system(cmd)
        return error
    
    elif 'CSAPR' in radar_name.upper():
        config_file = '/net/denali/storage2/radar2/people/idariash/home/Utiles/DROPS2/RE_CSAPR_RR_update.v6_b1.ini'
        '/net/denali/storage2/radar2/people/idariash/home/Utiles/DROPS2/RE_CSAPR_RR_update.v5.ini';
        [input_path, filename] = os.path.split(file)
        file_drops_tmp = os.path.join(output_path, filename)
        
        cmd = (HydroClass_exe + ' ' + file + ' -o ' + file_drops_tmp + 
        ' -c ' + config_file + ' -s ' + sounding_file + ' -m VHS -d ' + membership_functions)
        
        os.system(cmd)
        
        cmd = 'ncks -A -v fixed_angle ' + file + ' ' + file_drops_tmp 
        os.system(cmd)
        
        cmd = 'ncks -A -v sweep_mode ' + file + ' ' + file_drops_tmp
        os.system(cmd)
        
        cmd = 'ncks -A -v time ' + file + ' ' + file_drops_tmp
        os.system(cmd)
             
        return 0
     
    elif '1' in radar_name.upper():
        config_file = '/net/denali/storage2/radar2/people/idariash/home/Utiles/DROPS2/RE_RMA_RR_update.v5.ini';
        [input_path, filename] = os.path.split(file)
        file_drops_tmp = os.path.join(output_path, filename)
        
        cmd = (HydroClass_exe + ' ' + file + ' -o ' + file_drops_tmp + 
        ' -c ' + config_file + ' -s ' + sounding_file + ' -m VHS -d ' + membership_functions)
        
        os.system(cmd)
        
        cmd = 'ncks -A -v fixed_angle ' + file + ' ' + file_drops_tmp 
        os.system(cmd)
        
        cmd = 'ncks -A -v sweep_mode ' + file + ' ' + file_drops_tmp
        os.system(cmd)
        
        cmd = 'ncks -A -v time ' + file + ' ' + file_drops_tmp
        os.system(cmd)
             
        return 0
        
    else:
        print('Radar not found ')
        return 0

    
    
def run_drops(input_path, output_path):   #YYYYMMDD
    DataPath = input_path
        
    files = [f for f in glob.glob(DataPath +"/*.nc", recursive=True)]
    files.sort()
    
    
    for file in files:
        drops_chivo(file, output_path)
        #break 
    
    # pool = mp.Pool(20)
    # pool.map_async(drops_chivo, files)
    # pool.close()
    # pool.join()
    
    # cmd = ('rm ' + output_path_drops + '/*.nc')
    # os.system(cmd)
    

