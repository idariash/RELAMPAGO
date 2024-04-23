#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 16:24:52 2021

@author: idariash
"""

import sys
sys.path.append('/net/nasstore/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/plot')
import chivo_grid_utl as utl

chivo_att_corrected_filename = '/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/level_1b.2/2018/12/14/chivo.1b.2.20181214_020047.REL_PFAR360.nc'
'/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/level_1b.2/2018/12/14/chivo.1b.2.20181214_023049.REL_PNL360A.nc'

csapr_att_corrected_filename = '/net/k2/storage/projects/CSAPR/level_1b.2/corcsapr2cfrppiqcM1.1b.20181214.020004.nc'
'/net/k2/storage/projects/CSAPR/level_1b.2/corcsapr2cfrppiqcM1.1b.20181214.023003.nc'

chivo_original_filename = '/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/level_1a/2018/12/14/chivo.1a.20181214_020047.REL_PFAR360.nc'
'/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/level_1a/2018/12/14/chivo.1a.20181214_023049.REL_PNL360A.nc'

csapr_original_filename = '/net/k2/storage/projects/CSAPR/DOE_b1/corcsapr2cfrppiqcM1.b1.20181214.020004.nc'
'/net/k2/storage/projects/CSAPR/DOE_b1/corcsapr2cfrppiqcM1.b1.20181214.023003.nc'


utl.grid_dualRadar_att_corrected(chivo_att_corrected_filename, 
                                 csapr_att_corrected_filename)

utl.grid_dualRadar_original(chivo_original_filename, 
                                 csapr_original_filename)