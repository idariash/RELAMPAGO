#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 12:50:54 2020

@author: idariash
"""

import pyart

filename = '/net/denali/storage2/radar2/people/idariash/home/Field_Campaigns/Relampago/Analisis/GPM/20181206/CHIVO/NetCDF_drops/chivo.1b.20181206_052042.REL_PFAR360.nc'
#'/net/denali/storage2/radar2/people/idariash/home/Field_Campaigns/Relampago/Analisis/GPM/20181206/NetCDF_drops/cfrad.20181206_052008.000_to_20181206_052650.950_1_SUR.nc'
#'/net/denali/storage2/radar2/people/idariash/home/Field_Campaigns/Relampago/Analisis/GPM/20190131/drops/corcsapr2cfrppiM1.a1.20190131.223003_disdrometer.nc'
#'/net/denali/storage2/radar2/people/idariash/home/Field_Campaigns/Relampago/Analisis/GPM/20190113/azimuth_corrected_data/drops_book/chivo.1b.20190113_040044.REL_PFAR360_book.nc'
#'/net/denali/storage2/radar2/people/idariash/home/Field_Campaigns/Relampago/Analisis/GPM/20190113/azimuth_corrected_data/drops_book/chivo.1b.20190113_040044.REL_PFAR360_disdrometer.nc' 
#'/net/denali/storage/radar/RELAMPAGO/quality_controlled_data/level_1b/2018/12/06/chivo.1b.20181206_052042.REL_PFAR360.nc'
#'/net/denali/storage/radar/RELAMPAGO/quality_controlled_data/level_1b_2/2019/01/13/chivo.1b.20190113_040044.REL_PFAR360.nc'
#'/net/denali/storage/radar/RELAMPAGO/quality_controlled_data/level_1b/2019/01/13/chivo.1b.20190113_040044.REL_PFAR360.nc' 
#'/net/denali/storage/radar/CSAPR/DROPS/corcsapr2cfrppiqcM1.b1.20190131.223003.nc'
#'/net/denali/storage/radar/RELAMPAGO/quality_controlled_data/level_1b_3/2019/01/13/chivo.1b.20190113_040044.REL_PFAR360.nc'
#filename = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/GPM/20190113/data/NetCDF_drops/level_1b_test/cfrad.20190113_040044.721_to_20190113_040600.628_col-radar_REL_PFAR360_SUR.nc'

filename_out = '/net/denali/storage2/radar2/people/idariash/home/Field_Campaigns/Relampago/Analisis/GPM/20181206/CHIVO/NetCDF_drops/chivo.1b.20181206_052042.REL_PFAR360.uf' 
#'/net/denali/storage2/radar2/people/idariash/home/Field_Campaigns/Relampago/Analisis/GPM/20181206/NetCDF_drops/cfrad.20181206_052008.000_to_20181206_052650.950_1_SUR.uf'
#'/net/denali/storage2/radar2/people/idariash/home/Field_Campaigns/Relampago/Analisis/GPM/20190131/drops/corcsapr2cfrppiM1.a1.20190131.223003_disdrometer.uf' 
#'/net/denali/storage2/radar2/people/idariash/home/Field_Campaigns/Relampago/Analisis/GPM/20190113/azimuth_corrected_data/drops_book/chivo.1b_2.20190113_040044.REL_PFAR360_book.uf' 

#uf_field_names = {'attenuation_corrected_reflectivity_h': 'CZ'}
radar = pyart.io.read(filename)
#radar.azimuth['data'] = radar.azimuth['data']


pyart.io.write_uf(filename_out, radar)
#pyart.io.write_uf(filename_out, radar, uf_field_names = uf_field_names)
