#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 14:33:06 2022

@author: eesopee
"""

'''
This file contains useful info about MISTRAL
'''

import numpy as np
import os
import sys
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.convolution import Gaussian2DKernel


mistral_path = "/home/eesopee/mistral/mistral_simulator/config_files"

sampling_frequency = 244 * u.Hz #Hz

focal_plane_scale = 2.5 * u.arcsec / u.mm #arcsec/mm

fwhm = 12 * u.arcsec

focal_plane = np.genfromtxt(os.path.join(mistral_path, "pixel_positions_415.csv"), unpack=True, delimiter=",", skip_header=1)

feed_id = focal_plane[0]
feed_x = focal_plane[1]
feed_y = focal_plane[2]

def psf(pixel_scale, fwhm=fwhm):
    
    fwhm_pixel = fwhm / pixel_scale
    print(fwhm_pixel)
    sigma = fwhm_pixel.value/2.633
    print(sigma)
    return Gaussian2DKernel(x_stddev = sigma, y_stddev = sigma)

'''
La feed table vuole:
    Index xOffset yOffset relativePower

Forse dovremmo metterci pure i parametri dei cerchi per la pipeline.

'''
