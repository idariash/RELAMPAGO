#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 12:50:54 2020

@author: idariash
"""

import pyart

filename = '/net/denali/storage/radar/CSAPR/DOE_b1/corcsapr2cfrppiqcM1.b1.20190131.223003.nc'
#'/net/denali/storage/radar/RELAMPAGO/quality_controlled_data/level_1b_3/2019/01/13/chivo.1b.20190113_040044.REL_PFAR360.nc'
#filename = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/GPM/20190113/data/NetCDF_drops/level_1b_test/cfrad.20190113_040044.721_to_20190113_040600.628_col-radar_REL_PFAR360_SUR.nc'
#filename = '/net/denali/storage/radar/RELAMPAGO/quality_controlled_data/level_1b/2019/01/13/chivo.1b.20190113_040044.REL_PFAR360.nc'
filename_out = '/net/nasstore/students/GRAD/ECE/idariash/home/IDL/ground_radar/CSAPR/1CUF/corcsapr2cfrppiqcM1.b1.20190131.223003.uf'
uf_field_names = {'attenuation_corrected_reflectivity_h': 'CZ'}
radar = pyart.io.read(filename)
#radar.azimuth['data'] = radar.azimuth['data']
pyart.io.write_uf(filename_out, radar, uf_field_names = uf_field_names)
