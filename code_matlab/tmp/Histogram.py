# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 13:32:55 2018

@author: idariash
"""

#!/usr/bin/env python
 
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
 
 
# example data
mu = 0 # mean of distribution
sigma = 15 # standard deviation of distribution
#x = mu + sigma * np.random.randn(10000)
V_I = test.IQdata[100].V_I
V_Q = test.IQdata[100].V_Q
V_I = np.asarray(V_I)
V_Q = np.asarray(V_Q)
#x = test.IQdata[7].V_I
x = V_Q**2 + V_I**2
 
num_bins = 32
# the histogram of the data
n, bins, patches = plt.hist(x, num_bins, [0,0.000001], facecolor='blue', alpha=0.5)
 
# add a 'best fit' line
#y = mlab.exppdf(bins, mu, sigma)
#plt.plot(bins, y, 'r--')
plt.xlabel('$I^2+Q^2$')
plt.ylabel('Counts')
plt.title(r'Horizontal Channel')
 
# Tweak spacing to prevent clipping of ylabel
plt.subplots_adjust(left=0.15)
plt.show()