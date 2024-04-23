#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 12:55:11 2020

@author: idariash
"""

import sys
import os

import glob

import multiprocessing as mp


output_path_drops = '/net/denali/storage2/radar2/people/idariash/home/Field_Campaigns/Relampago/Analisis/GPM/20190113/azimuth_corrected_data/'

HydroClass_exe = '/net/denali/storage2/radar2/tmp/Ivan/Utiles/DROPS2/hydroclass';
sounding_file = ' ';

config_file = '/net/denali/storage2/radar2/tmp/Ivan/Utiles/DROPS2/ARDEC_CHIVO.v1.ini';
membership_functions = '/net/denali/storage2/radar2/tmp/Ivan/Utiles/DROPS2/membership_functions/';

def drops_csapr(file):
     
    [input_path, filename] = os.path.split(file)
    file_drops_tmp = os.path.join(output_path_drops, filename)
    
    cmd = (HydroClass_exe + ' ' + file + ' -o ' + file_drops_tmp + 
    ' -c ' + config_file + ' -s ' + sounding_file + ' -m VHS -d ' + membership_functions)
    
    os.system(cmd)
    
    cmd = 'ncks -A -v fixed_angle ' + file + ' ' + file_drops_tmp 
    os.system(cmd)
    
    cmd = 'ncks -A -v sweep_mode ' + file + ' ' + file_drops_tmp
    os.system(cmd)
    
    cmd = 'ncks -A -v time ' + file + ' ' + file_drops_tmp
    os.system(cmd)
    print(filename)
    

file = '/net/denali/storage2/radar2/tmp/Ivan/CHIVO/20200520/20200521/cfrad.20200521_024736.208_to_20200521_025157.999_col-radar_ARDEC_PPI_SUR.nc'

drops_csapr(file)

