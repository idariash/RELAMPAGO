#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ivan Arias
2021/02/15

This funtion are to run_drops from 1a data

"""

import os
import glob
import pyart
import netCDF4

config_file = '/net/k2/storage/people/idariash/home/Utiles/DROPS2/Sail_Xband_ConfigFile.ini'


HydroClass_exe = '/net/denali/storage2/radar2/people/idariash/home/Utiles/DROPS2/hydroclass' # This config file is for elbert

membership_functions = '/net/k2/storage/people/idariash/home/Utiles/DROPS2/membership_functions/';


def run_drops_commands_only(input_path, output_path): 
    # run DROPS for the -disag files
    files = [f for f in glob.glob(input_path +"/**/*.nc", recursive=True)]
    
    for file in files:
        
        [input_path, filename] = os.path.split(file)
        file_drops_tmp = os.path.join(output_path, filename)
        
        print(HydroClass_exe + ' ' + file + ' -o ' + file_drops_tmp + 
        ' -c ' + config_file + ' -m VHS -d ' + membership_functions) # this is the comand that runs drops
        
        cmd = 'ncks -A -v fixed_angle ' + file + ' ' + file_drops_tmp # comands to make drops output compatible with pyart
        print(cmd)
        
        cmd = 'ncks -A -v sweep_mode ' + file + ' ' + file_drops_tmp
        print(cmd)
        
    return 0

input_path = '/net/k2/storage/people/idariash/home/tmp/drops_pat/original_nc'
output_path = '/net/k2/storage/people/idariash/home/tmp/drops_pat/drops'

run_drops_commands_only(input_path, output_path)
