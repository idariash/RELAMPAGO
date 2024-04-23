#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 15:44:57 2020

@author: idariash
"""

import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
import multidop
import pyart
import pydda
import tempfile
import os
import glob
import time
import math
from wrf import  vertcross


filename_chivo = '/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/level_1a/2018/12/14/chivo.1a.20181214_030048.REL_PNL360A.nc'
filename_csapr = '/net/k2/storage/projects/CSAPR/DOE_b1/corcsapr2cfrppiqcM1.b1.20181214.030003.nc'


r2 = pyart.io.read(filename_chivo)
r1 = pyart.io.read(filename_csapr)


# Rotate axis

r2.azimuth['data'] = (r2.azimuth['data'] + 45) % 360
r1.azimuth['data'] = (r1.azimuth['data'] + 45) % 360

chivo_lon = r2.longitude['data']
chivo_lat = r2.latitude['data']

csapr_lon = r1.longitude['data']
csapr_lat = r1.latitude['data']

factor = 1.08

chivo_csapr_angle = math.pi + math.atan((chivo_lat - csapr_lat)/(chivo_lon - csapr_lon))

chivo_csapr_r = math.sqrt( (chivo_lat - csapr_lat)**2 + (chivo_lon - csapr_lon)**2)

chivo_csapr_new_angle = chivo_csapr_angle - math.pi/4

csapr_lon_r = chivo_lon + factor*chivo_csapr_r*math.cos(chivo_csapr_new_angle)

csapr_lat_r = chivo_lat + factor*chivo_csapr_r*math.sin(chivo_csapr_new_angle)