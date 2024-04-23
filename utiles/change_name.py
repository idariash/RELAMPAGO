#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  3 16:26:17 2021

@author: idariash
"""

import os
import glob


input_path = '/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/level_1b.2'

files = [f for f in glob.glob(input_path +"/**/*.nc", recursive=True)]
files.sort()

for file in files:
    [parent_path, filename] = os.path.split(file)

    
    filename = filename.replace('1b', '1b.2')
    
    renamed_file = os.path.join(parent_path, filename)
    
    cmd = 'mv  ' + file + ' ' + renamed_file 
    
    os.system(cmd)
    