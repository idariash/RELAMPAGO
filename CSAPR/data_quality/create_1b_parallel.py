#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 12:55:11 2020

@author: idariash
"""

import os

import glob

import multiprocessing as mp

output_path_drops = '/net/denali/storage/radar/CSAPR/DROPS/';
HydroClass_exe = '/net/denali/storage2/radar2/tmp/Ivan/Utiles/DROPS2/hydroclass';
sounding_file = ' ';

config_file = '/net/denali/storage2/radar2/tmp/Ivan/Utiles/DROPS2/RE_CSAPR_RR_update.v5.ini';
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
    
    
path = '/net/denali/storage/radar/CSAPR/DOE_b1/'
files = [f for f in glob.glob(path + "**/*2018*.nc", recursive=True)]
files.sort()


pool = mp.Pool(8)
pool.map_async(drops_csapr, files)
pool.close()
pool.join()
    
    # filename_drops = file_drops_tmp       
    
    # print (str(x) + '/' + str(len(Directory)) + ' ' + Directory[x])
    
    # error = dq_1b.correctData(file_drops_tmp, output_path, azimuthTable_path)
    # log = open(log_name,'a')
    
    # if error == 0:
    #     log.write(str(x) + '/' + str(len(Directory)) + ' ' + Directory[x] + ' done \n')
    # else:
    #     log.write(str(x) + '/' + str(len(Directory)) + ' ' + Directory[x] + ' error \n')
   
    # log.close()
    # cmd = 'rm ' + file_drops_tmp
    # os.system(cmd)
    
    
