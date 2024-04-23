#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 14:48:21 2022

Read rain from files

@author: idariash
"""

import pandas as pd
from datetime import datetime
import glob
import numpy as np

DataPath = '/net/k2/storage/people/idariash/home/RELAMPAGO/Analisis/Rain_Validation/SMN_RainGauges'

files = [f for f in glob.glob(DataPath +"/**/*.xls", recursive=True)]

files.sort()

rainGauge_IDs = np.zeros(len(files))
rainGauges_accumulation = np.zeros(len(files))
I = 0

for file in files:
    
    xl = pd.ExcelFile(file)
    
    sheet_names = xl.sheet_names
    
    Dec14_UTC_start = datetime.strptime('13/12/2018 21:00:00', '%d/%m/%Y %H:%M:%S')
    
    Dec14_UTC_end = datetime.strptime('14/12/2018 20:59:59', '%d/%m/%Y %H:%M:%S')
    
    rainGauge_accumulation = 0
    
    
    for i in range(0, len(sheet_names)):
        
        df = pd.read_excel(file, sheet_name = i, header= 2)
        
        
        if 'bla' in sheet_names[i]:
            continue
        
        rainGauge_ID = df.iloc[4,2]
        
        for j in range(0, len(df)):
            
            try:
                rainGauge_time = datetime.strptime(df.iloc[j,1], '%d/%m/%Y %H:%M:%S')
            except:
                continue
            
            if (Dec14_UTC_start < rainGauge_time) & (rainGauge_time < Dec14_UTC_end):
                
                try:
                    rainGauge_accumulation = rainGauge_accumulation + df.iloc[j,3]
                except:
                    rainGauge_accumulation = rainGauge_accumulation + int(df.iloc[j,3])
                
                print(rainGauge_time)
                
    print(rainGauge_ID)
    rainGauge_IDs[I] = rainGauge_ID
    rainGauges_accumulation[I] = rainGauge_accumulation
    I = I + 1
