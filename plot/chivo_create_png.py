#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 21:04:55 2019

@author: idariash
"""

import chivo_plot as plot
import os
import matplotlib.pyplot as plt


DataPath = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Attenuation/NetCDF/20190125_1900to2130/DROPS'
PNG_Path = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Attenuation/NetCDF/20190125_1900to2130/PNG'

Directory = os.listdir(DataPath)

for file in Directory:
    if 'RHI' in file:
        continue
        print (file)
        filename = os.path.join(DataPath, file)
        file_name = filename[len(filename)-46:len(filename)-31]
        i = 0
        while True:
            try: 
                fig = plot.rhi_drops(filename, i)
                fig_name = os.path.join(PNG_Path, 'rhi/rhi_' + file_name + '_' + str(i) + '.png')
                fig.savefig(fig_name)
                fig.clear()
                plt.close(fig)
                i = i + 1
            except:
                break
    if 'PNL' in file:
        filename = os.path.join(DataPath, file)
        file_name = file_name = filename[len(filename)-49:len(filename)-34]
        i = 0
        while True:
            try: 
                fig = plot.ppi_drops(filename, i)
                fig_name = os.path.join(PNG_Path, 'ppi_100km/ppi_' + file_name + '_' + str(i) + '.png')
                fig.savefig(fig_name)
                fig.clear()
                plt.close(fig)
                i = i + 1
            except:
                break