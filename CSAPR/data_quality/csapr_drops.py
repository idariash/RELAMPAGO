#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 12:55:11 2020

@author: idariash
"""

import os



output_path_drops = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/GPM/20190113/data/NetCDF_drops/level_1b_test'
#output_path_drops = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/GPM_CSAPR/GPM/DROPS/';
HydroClass_exe = '/net/denali/storage2/radar2/tmp/Ivan/Utiles/DROPS2/hydroclass';
sounding_file = ' ';

config_file = '/net/denali/storage2/radar2/tmp/Ivan/Utiles/DROPS2/RE_CHIVO_RR_update.v5.ini' 
#config_file = '/net/denali/storage2/radar2/tmp/Ivan/Utiles/DROPS2/RE_CSAPR_RR_update.v5.ini';
#config_file = '/net/denali/storage2/radar2/tmp/Ivan/Utiles/DROPS2/RE_CSAPR_Cband.v1.ini';
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
    
   
    
file = '/net/denali/storage/radar/RELAMPAGO/NetCDF/20190113/cfrad.20190113_040044.721_to_20190113_040600.628_col-radar_REL_PFAR360_SUR.nc'
#file = '/net/denali/storage/radar/CSAPR/DOE_b1/corcsapr2cfrppiqcM1.b1.20190131.223003.nc'
#file = '/net/denali/storage/radar/CSAPR/DOE_a1/corcsapr2cfrppiM1.a1.20190131.223003.nc'
drops_csapr(file)
