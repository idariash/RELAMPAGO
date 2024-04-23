#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 12:20:37 2021

@author: idariash
"""

import os
import glob


input_path = '/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/level_1b.2'

files = [f for f in glob.glob(input_path +"/**/*.nc", recursive=True)]
files.sort()

for file in files:
    [parent_path, filename] = os.path.split(file)
    AAAA = filename[9:13]
    MM = filename[13:15]
    DD = filename[15:17]
    file_organized = os.path.join(parent_path, AAAA, MM, DD, filename)
    
    try:    
        os.replace(file, file_organized)
    except:
        os.makedirs(os.path.join(parent_path, AAAA, MM, DD))
        os.replace(file, file_organized)