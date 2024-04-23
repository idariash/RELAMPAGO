#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 15:34:49 2019

@author: idariash
"""

import pyart
import os
import sys

filename = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/GPM/20190113/data/NetCDF_drops/cfrad.20190113_040044.721_to_20190113_040600.628_col-radar_REL_PFAR360_SUR_book.nc'
filename_out = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/GPM/20190113/data/NetCDF_drops/cfrad.20190113_040044.721_to_20190113_040600.628_col-radar_REL_PFAR360_SUR_book.uf'

radar =  pyart.io.read(filename)
radar.azimuth['data'] = radar.azimuth['data'] - 13
pyart.io.write_uf(filename_out, radar)

