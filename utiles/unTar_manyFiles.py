#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 16:31:15 2022

@author: idariash
"""
import os
import glob

path = '/net/k2/storage/projects/RELAMPAGO/RMA1/download_03'
'/net/k2/storage/projects/RELAMPAGO/RMA1/new_download_01'
'/net/k2/storage/projects/RELAMPAGO/RMA1'
files = [f for f in glob.glob(path + "/**/*.tar", recursive=True)]
files.sort()
for file in files:
    cmd = ('tar -xvf ' + file + 
           ' -C /net/k2/storage/projects/RELAMPAGO/RMA1/Original_data')
    os.system(cmd)