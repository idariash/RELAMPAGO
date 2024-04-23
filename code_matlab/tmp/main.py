# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 13:59:23 2018

@author: idariash
"""

import datetime
myFile = "/home/idariash/Desktop/Link to Ivan/Field_Campaigns/Relampago/Data/Cordoba_Radar/Time_Series/20171129T191331_524/RMA01_0123_02_20171129T191331_524Z_17.IQ"
test = invapIQ( myFile )
#test.dump()

test.getScans()

print("\nEnd")