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

'/net/k2/storage/people/idariash/home/Utiles/DROPS2/RE_CHIVO_RR_update.v8_1b_rain.ini'

'/net/k2/storage/people/idariash/home/Utiles/DROPS2/RE_CSAPR_RR_update.v8_b1_enhancedCoeff.ini'

'/net/k2/storage/people/idariash/home/Utiles/DROPS2/RE_CHIVO_RR_update.v8_1a_water_coated_ice.ini'

'/net/k2/storage/people/idariash/home/Utiles/DROPS2/RE_CHIVO_RR_update.v8_1a_water_coated_ice.ini'

'/net/k2/storage/people/idariash/home/Utiles/DROPS2/RE_RMA_RR_update.v5.ini'

'/net/k2/storage/people/idariash/home/Utiles/DROPS2/RE_CSAPR_RR_update.v7_b1.ini'

'/net/denali/storage2/radar2/people/idariash/home/Utiles/DROPS2/RE_CHIVO_RR_update.v7_1a.ini'

'/net/denali/storage2/radar2/people/idariash/home/Utiles/DROPS2/RE_CHIVO_RR_update.v6_1a.ini'

HydroClass_exe = '/net/k2/storage/people/idariash/home/Utiles/DROPS2/Whitney/hydroclass'

'/net/denali/storage2/radar2/people/idariash/home/Utiles/DROPS2/hydroclass'





 

membership_functions = '/net/k2/storage/people/idariash/home/Utiles/DROPS2/membership_functions/';

def add_fields_from_1a(radar, radar_1a):
    ncp_1a = radar_1a.fields['normalized_coherent_power'].copy()
    vel_1a = radar_1a.fields['velocity'].copy()
    width_1a = radar_1a.fields['spectrum_width'].copy()
    
    radar.add_field('normalized_coherent_power', ncp_1a, replace_existing=True)
    radar.add_field('velocity', vel_1a, replace_existing=True)
    radar.add_field('spectrum_width', width_1a, replace_existing=True)    
    
    return radar

def add_fields_from_b1(radar, radar_b1):
    ncp_b1 = radar_b1.fields['normalized_coherent_power'].copy()
    vel_b1 = radar_b1.fields['mean_doppler_velocity'].copy()
    width_b1 = radar_b1.fields['spectral_width'].copy()
    
    radar.add_field('normalized_coherent_power', ncp_b1, replace_existing=True)
    radar.add_field('velocity', vel_b1, replace_existing=True)
    radar.add_field('spectrum_width', width_b1, replace_existing=True)    
    
    return radar

def add_fields_from_rma_cfradial(radar, radar_cfradial):
    #ncp_b1 = radar_cfradial.fields['normalized_coherent_power'].copy()
    vel_cfradial = radar_cfradial.fields['VRAD'].copy()
    width_cfradial = radar_cfradial.fields['CM'].copy()
    
    #radar.add_field('normalized_coherent_power', ncp_cfradial, replace_existing=True)
    radar.add_field('velocity', vel_cfradial, replace_existing=True)
    radar.add_field('spectrum_width', width_cfradial, replace_existing=True)    
    
    return radar

def exportFile(radar, file, output_path):
    
    # radar.metadata['title'] = 'RELAMPAGO CSU-CHIVO Level 1b data'
    # radar.metadata['institution'] = 'Colorado State University'
    # radar.metadata['references'] = 'Data quality control of CHIVO'
    # radar.metadata['comment'] = 'Attenuation correction using disdrometer coefficient from the field'
    # radar.metadata['instrument_name'] = 'CSU-CHIVO'
    # radar.metadata['history'] = 'Created by idariash from 1a data, March 2021'
    # radar.metadata['source'] = 'Retrieved from sigmet original file'
    
    [path, filename_1b] = os.path.split(file)
    
    filename_1b = filename_1b.replace('b1', '1b')
    
    #filename_1b = filename_1b.replace('1a', '1b')

    
    filename = os.path.join(output_path, filename_1b)
    pyart.io.write_cfradial(filename, radar, time_reference=True, arm_time_variables=False)
    
    return 0


def join_files(input_path):
    files = [f for f in glob.glob(input_path +"/**/*.nc", recursive=True)]
    files.sort()
    
    radar = pyart.io.read(files[0])
    files.pop(0)
    for file in files:
        sweep = pyart.io.read(file)
        radar = pyart.util.join_radar(radar, sweep)
        
    return radar

def RadxConvert_1a(file, output_path):
    # Radx convert to make compatible the file with DROPS
    cmd = 'rm -r ' + os.path.join(output_path, '*')
    os.system(cmd)
    cmd = ('/usr/local/lrose/bin/RadxConvert -f ' + file + ' -disag ' + 
           ' -outdir ' +  output_path)
    os.system(cmd)
    
    return 0
    
    
def run_drops(input_path, output_path): 
    # run DROPS for the -disag files
    cmd = 'rm -r ' + os.path.join(output_path, '*')
    os.system(cmd)
    files = [f for f in glob.glob(input_path +"/**/*.nc", recursive=True)]
    
    for file in files:
        
        [input_path, filename] = os.path.split(file)
        file_drops_tmp = os.path.join(output_path, filename)
        
        cmd = (HydroClass_exe + ' ' + file + ' -o ' + file_drops_tmp + 
        ' -c ' + config_file + ' -m VHS -d ' + membership_functions) # this is the comand that runs drops
        
        os.system(cmd)
        
        cmd = 'ncks -A -v fixed_angle ' + file + ' ' + file_drops_tmp # comands to make drops output compatible with pyart
        os.system(cmd)
        
        cmd = 'ncks -A -v sweep_mode ' + file + ' ' + file_drops_tmp
        os.system(cmd)
        
    return 0

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
