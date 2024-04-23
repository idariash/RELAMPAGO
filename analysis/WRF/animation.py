#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 15:34:53 2023

@author: idariash
"""

import glob
import wrf_plot as plot
import wrf
from netCDF4 import Dataset

#'/glade/scratch/ivanarias/WRF_test/domains/RELAMPAGO_IOP17_Aerosol/'

path = '/net/k2/storage/people/idariash/home/WRF/WRF_test/domains/RELAMPAGO_IOP17_Morrison/ndown_wrf'

'/net/k2/storage/people/idariash/home/WRF/WRF_test/domains/RELAMPAGO_IOP17_Aerosol/ndown_wrf/'

#files = [f for f in glob.glob(path + "**/wrfout_d03_2018-12-14_0*", recursive=True)]

files = [f for f in glob.glob(path + "/wrfout_d01*", recursive=True)]

files.sort()

# for file in files:
#     for i in range(6):
#         fig = plot.horizontal_max_ref(file, i)
#         ncfile = Dataset(file)
#         time = wrf.extract_times(ncfile, timeidx = i)
#         time_text = str(time)[:16]
#         time_text = time_text.replace(':', '')
#         fig_name = ('/glade/u/home/ivanarias/CSU/RELAMPAGO/analysis/fig/ndwon/domain_3/wrfout_d3_aerosol_' + 
#                     time_text + '.png')
        
#         fig.savefig(fig_name)
        
for file in files:
        for i in range(10):
            fig = plot.horizontal_max_ref(file, i)
            ncfile = Dataset(file)
            time = wrf.extract_times(ncfile, timeidx = i)
            time_text = str(time)[:15]
            time_text = time_text.replace(':', '')
            fig_name = ('/net/k2/storage/people/idariash/home/WRF/analysis/Morrison/figures/wrfout_d4_aerosol_' + 
                        time_text  + str(i) + '.png')
            
            fig.savefig(fig_name)