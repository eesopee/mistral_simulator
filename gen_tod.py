#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 11:06:43 2022

@author: eesopee
"""

'''

This script generates "clean" TODs from a synthetic fits map.

Inputs:
    
    DISCOS fits file, to provide the telescope pointings
    synthetic fits map

Outputs:
    
    DISCOS fits file with attached data table

'''

from astropy.io import fits
import astropy
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
import mistral
import argparse
from astropy.coordinates import EarthLocation
import astropy.time
import time
import sys

import mistral



if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--discos_fits")
    parser.add_arguments("--map")
    
    args = parser.parse_args()
    
    
    
    
    
    