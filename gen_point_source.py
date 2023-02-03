#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 11:22:45 2022

@author: eesopee
"""

'''
This script generates the image of a point source
'''

import numpy as np
from astropy.convolution import convolve
import matplotlib.pyplot as plt
import mistral
import astropy.units as u

shape = 200
flux=1

img = np.zeros(shape=(shape,shape)) 

psf = mistral.psf(pixel_scale = 1 * u.arcsec/u.pix)

img[shape//2,shape//2] = flux #arbitrary flux

img = convolve(img, psf)

plt.imshow(img)

print(np.max(img))


