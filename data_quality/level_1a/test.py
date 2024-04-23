#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 10:52:51 2019

@author: idariash
"""

import pyart
import chivo_dataQuality_1a as dq_1a

filename_HYDRO = '/net/denali/storage2/radar2/tmp/RELAMPAGO/RAW/2018/11/29/COL20181129_170002' 
#'/net/denali/storage2/radar2/tmp/RELAMPAGO/RAW/2018/11/29/COL20181129_164002'
filename_PPI = '/net/denali/storage2/radar2/tmp/RELAMPAGO/RAW/2018/11/29/COL20181129_170042'
#'/net/denali/storage2/radar2/tmp/RELAMPAGO/RAW/2018/11/29/COL20181129_164043'
filename_RHI = '/net/denali/storage2/radar2/tmp/RELAMPAGO/RAW/2018/11/29/COL20181129_170604' 
#'/net/denali/storage2/radar2/tmp/RELAMPAGO/RAW/2018/11/29/COL20181129_170512' 
#'/net/denali/storage2/radar2/tmp/RELAMPAGO/RAW/2018/11/29/COL20181129_164510'
output_path = '/net/denali/storage2/radar2/tmp/Ivan/CSU/CHIVO/src/python/data_quality/level_1a/output'
azimuth_error = 1.97
Zdr_bias = 0.75

HYDRO = pyart.io.read(filename_HYDRO)
PPI = pyart.io.read(filename_PPI)
RHI = pyart.io.read(filename_RHI)

dq_1a.Zdr_correction(HYDRO, Zdr_bias)
dq_1a.azimuth_correction(HYDRO, azimuth_error)
dq_1a.remove_field(HYDRO, 'radar_echo_classification')
dq_1a.flip_phidp(HYDRO, 180)
dq_1a.exportFile(HYDRO, output_path)



dq_1a.Zdr_correction(PPI, Zdr_bias)
dq_1a.azimuth_correction(PPI, azimuth_error)
dq_1a.remove_field(PPI, 'radar_echo_classification')
dq_1a.flip_phidp(PPI, 180)
dq_1a.exportFile(PPI, output_path)



dq_1a.Zdr_correction(RHI, Zdr_bias)
dq_1a.azimuth_correction(RHI, azimuth_error)
dq_1a.remove_field(RHI, 'radar_echo_classification')
dq_1a.flip_phidp(RHI, 180)
dq_1a.exportFile(RHI, output_path)


#%%
a = RHI.fields





