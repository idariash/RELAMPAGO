#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 14:15:08 2021

@author: idariash
"""

import sys
sys.path.append('./') # path where your ploting function are
import radar_plot as plot
import importlib

importlib.reload(plot)

filename = './cfrad.20181130_024042.935_to_20181130_024558.468_col-radar_REL_PFAR360_SUR.nc'

fig = plot.ppi_reflectivity(filename, 2)
