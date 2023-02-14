#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 15:23:33 2022

@author: eesopee
"""

from astropy.io import fits
import sys
import matplotlib.pyplot as plt
import numpy as np


scan = fits.open("output_test.fits")
data_table = scan["DATA TABLE INTERP"].data
iq_table = scan["IQ TABLE"].data

feed_table = scan["FEED TABLE"].data

px_x = feed_table["xOffset"]
px_y = feed_table["yOffset"]


#plt.plot(data_table["raj2000"], data_table["decj2000"])
fig = plt.figure(figsize=(7,7))

start=100
stop=102
skip=1

from scipy.stats import binned_statistic_2d

feed_maps = []

for id in feed_table["id"]:
    
    xoffset = feed_table["xOffset"][id]*2/3600
    yoffset = feed_table["yOffset"][id]*2/3600
    gain = feed_table["relativePower"][id]
    
    ra = data_table["raj2000_interpolated"][start:stop:skip]+xoffset
    dec = data_table["decj2000_interpolated"][start:stop:skip]+yoffset
    s = iq_table["chI_{:s}".format(str(id).zfill(3))][start:stop:skip]/gain
    
    
    
    fluxmap, a,b,c = binned_statistic_2d(x=np.array(ra)[~np.isnan(ra)], 
                                 y=np.array(dec)[~np.isnan(ra)], 
                                 values=np.array(s)[~np.isnan(ra)],
                                 bins=32,
                                 statistic=np.nanmean)
    
    feed_maps.append(fluxmap)
    
plt.imshow(np.nanmean(np.array(feed_maps),axis=0))
    
delta=0.03

'''
plt.xlim(min(data_table["raj2000_interpolated"]-delta),
         max(data_table["raj2000_interpolated"]+delta))
         
plt.ylim(min(data_table["decj2000_interpolated"]-delta),
         max(data_table["decj2000_interpolated"]+delta))
'''
plt.grid()
plt.xlabel("raj2000")
plt.ylabel("decj2000")
