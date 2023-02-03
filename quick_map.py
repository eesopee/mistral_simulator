#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 15:23:33 2022

@author: eesopee
"""

from astropy.io import fits
import sys
import matplotlib.pyplot as plt

scan = fits.open(sys.argv[1])

data_table = scan["DATA TABLE"].data

#plt.plot(data_table["raj2000"], data_table["decj2000"])
plt.scatter(data_table["az"], data_table["el"])

plt.show()