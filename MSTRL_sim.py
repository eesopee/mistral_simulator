# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 16:05:56 2023

file precedente a gen_tod

@author: eleobar
"""

import scanlib
import mistralib as ml
import sys
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import scipy
from scipy import stats
import time
import numpy as np
import pandas as pd
import astropy.units.si as u
import astropy.wcs
import astropy.io.fits as fits
import astropy.coordinates as apcoord
from tqdm import tqdm

sys.path.append('C:/Users/Eleonora/Desktop/mistralib')


"""
DEFINE A DETECTOR ARRAY AND ITS TELESCOPE
"""
T_env = 275

#This line initializes SRT with default parameters. It also contains information about optics such as efficiency and temperature

SRT = ml.Antenna(antenna = "SRT")

#defining the optical system of MISTRAL: lenses and filter chain

elements = [ml.OpticalElement("window", eff = lambda x : 0.8+ 0*x, emi=3e-2, T=T_env),
            ml.OpticalElement("L1", eff = lambda x : 0.9 + 0*x, emi=0, T=0),
            ml.OpticalElement("L2", eff = lambda x : 0.9 + 0*x, emi=0, T=0),
            ml.OpticalElement("bandpass20", eff = lambda x : ml.Filter(filename="C:/Users/Eleonora/Desktop/mistralib/mistral_bandpass_cm-1.txt", unit=1/u.cm).efficiency(x), emi=0, T=0)
            #,ml.OpticalElement("bandpass30", eff = lambda x : ml.Filter(filename="mistral_bandpass_cm-1.txt", unit=1/u.cm).efficiency(x), emi=0, T=0)
            ]

#defining the focal plane of MISTRAL: xy position of each pixel and its response.

fp_filename = "C:/Users/Eleonora/Desktop/mistralib/MISTRAL_smol.txt" #fp_filename = "/content/drive/My Drive/university/magistrale/Tesi/MISTRAL/scanning/mistralib/pixels_paiella.txt"
eff_filename = "C:/Users/Eleonora/Desktop/mistralib/kids_efficiency.txt"

fp = ml.FocalPlane(fp_filename, eff_filename)

#Now we define the receiver object containing everything we need:

mistral = ml.Receiver(antenna=SRT, optical_elements = elements, focal_plane = fp)



'''
Opening the FITS file we want to scan. It should have a valid header with the info needed to build a World Coordinate System (WCS)
'''

fits_name = "C:/Users/Eleonora/Desktop/mistralib/point_source_f=1e-3.fits" 

hdul = fits.open(fits_name)
map_hdu = hdul["PRIMARY"]

hdul.info()
map_hdu.header

map_wcs = astropy.wcs.WCS(map_hdu.header) #this function extracts the WCS from the FITS header
sky_map = map_hdu.data



'''
Defining the scan curve and then generate the timestream (without noise)
'''

T_scan = 60*u.s #540 * u.s
scan_center = apcoord.SkyCoord(map_wcs.wcs.crval[0], map_wcs.wcs.crval[1],unit=(u.deg, u.deg))

tau_radial = T_scan/22
scan_radius = 1 * u.arcmin
phi_0=0
phi_1=1

daisy_args = (scan_radius, phi_0, phi_1, tau_radial)

daisy = scanlib.DaisyScan(scanlib.daisy, T_scan, mistral.sampling_frequency, scan_center, *daisy_args)

timestreams = np.array(scanlib.gen_timestreams(daisy,mistral,map_hdu)) #The output is a list of arrays, each array is the time-stream of the i-th feed.
 
#_=[plt.plot(daisy.t,timestreams[i][2]) for i in range(0,mistral.focal_plane.number_of_feeds)]


#copertura telescopio
totx = np.array([])
toty = np.array([])
timestream = np.array([])
for feed in range(0,mistral.focal_plane.number_of_feeds):
  timestream = np.concatenate((timestream, timestreams[feed][2]))
  totx = np.concatenate((totx, timestreams[feed][0]))
  toty = np.concatenate((toty, timestreams[feed][1]))

image = np.ones(len(daisy.x))
cop1, x_edge, y_edge, binnumber = scipy.stats.binned_statistic_2d(x=np.array(totx,dtype="float"),y=np.array(toty,dtype="float"), values=np.array(timestream,dtype="float"), bins=64, statistic="count")

cop2=[plt.plot(timestreams[i][0],timestreams[i][1], color="k", alpha=0.01) for i in range(0,mistral.focal_plane.number_of_feeds)]
'''
fig, (ax1,ax2) = plt.subplots(ncols=2,nrows=1, figsize=(15,7))
ax1.imshow(cop1, norm=LogNorm())
plt.title("Copertura telescopio")

ax2.imshow(cop2, norm=LogNorm())
plt.title("Daisy scan path")
'''

daisy_path = plt.scatter(totx, toty, c=timestream, marker=".")
plt.title("Daisy scan")


'''
METTERE INSIEME QUESTE TRE FIGURE
'''
fig, (ax1,ax2,ax3) = plt.subplots(ncols=3,nrows=1, figsize=(23,7))

ax1.imshow(daisy_path)
plt.title("Daisy scan")
 
ax2.imshow(cop1, norm=LogNorm())
plt.title("Copertura telescopio")

ax3.imshow(cop2, norm=LogNorm())
plt.title("Copertura telescopio pixel")




'''
"real" map-making: we will use the scipy.binned_statistic_2d module -> see Airy disk

#timestream += np.random.normal(0.5,1, size=(timestream.shape)) #adding a gaussian noise without any particular physical meaning

map, x_edge, y_edge, binnumber = scipy.stats.binned_statistic_2d(x=np.array(totx, dtype=float),y=np.array(toty, dtype=float), values=np.array(timestream, dtype=float), bins=64, statistic="median")
std, x_edge, y_edge, binnumber = scipy.stats.binned_statistic_2d(x=np.array(totx, dtype=float),y=np.array(toty, dtype=float), values=np.array(timestream, dtype=float), bins=64, statistic="std")
count, x_edge, y_edge, binnumber = scipy.stats.binned_statistic_2d(x=np.array(totx, dtype=float),y=np.array(toty, dtype=float), values=np.array(timestream, dtype=float), bins=64, statistic="count")

fig, (ax1,ax2) = plt.subplots(ncols=2,nrows=1,figsize=(15,7))
ax1.imshow(map, norm=LogNorm())
plt.title("Map")

snr_map = map/(std/np.sqrt(count))

ax2.imshow(snr_map,norm=LogNorm())
plt.title("S/N map")
'''














