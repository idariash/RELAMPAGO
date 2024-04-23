#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 12:55:10 2019

@author: idariash
"""

import pyart
import netCDF4
import datetime
from matplotlib import pyplot as plt
import numpy

filename = '/net/denali/storage2/radar2/tmp/Ivan/CSU/CHIVO/res/data/test_data_quality/corcsapr2cfrppiM1.a1.20190125.193003.nc'
CSAPR =  pyart.io.read(filename)

#%%

a = CSAPR.fields
latitude_csapr = CSAPR.latitude

#%%
chivo = '/net/denali/storage2/radar2/tmp/Ivan/CSU/CHIVO/res/data/test_data_quality/chivo.1a.20190125_193049.REL_PNL135A.nc'
CHIVO = pyart.io.read(chivo)
c = CHIVO.fields
latitude_chivo = CHIVO.latitude

CSAPR.latitude['_FillValue'] = numpy.float64(CSAPR.latitude['_FillValue'])
CSAPR.latitude['data'] = numpy.float64(CSAPR.latitude['data'])
CSAPR.latitude['valid_max'] = numpy.float64(CSAPR.latitude['valid_max'])
CSAPR.latitude['valid_min'] = numpy.float64(CSAPR.latitude['valid_min'])

CSAPR.longitude['_FillValue'] = numpy.float64(CSAPR.longitude['_FillValue'])
CSAPR.longitude['data'] = numpy.float64(CSAPR.longitude['data'])
CSAPR.longitude['valid_max'] = numpy.float64(CSAPR.longitude['valid_max'])
CSAPR.longitude['valid_min'] = numpy.float64(CSAPR.longitude['valid_min'])

CSAPR.altitude['_FillValue'] = numpy.float64(CSAPR.altitude['_FillValue'])
CSAPR.altitude['data'] = numpy.float64(CSAPR.altitude['data'])
#CSAPR.altitude['valid_max'] = numpy.float64(CSAPR.altitude['valid_max'])
#CSAPR.altitude['valid_min'] = numpy.float64(CSAPR.altitude['valid_min'])

latitude_csapr_new = CSAPR.latitude
altitude =  CSAPR.altitude

#Create new NetCDF file
output_path = '/net/denali/storage2/radar2/tmp/Ivan/CSU/CHIVO/res/data/test_data_quality'

filename = output_path + '/corcsapr2cfrppiM1.a1.20190125.193003_reformat_test.nc'
pyart.io.write_cfradial(filename, CSAPR, format='NETCDF4', time_reference=True, arm_time_variables=False)

#%%
#cp = deepcopy(CHIVO.fields['latitude']['data'])
#r1.add_field_like('REF', 'DT', cp, replace_existing=True)

