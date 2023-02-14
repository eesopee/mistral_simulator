#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 15:26:00 2023

@author: eesopee
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.time import Time
from scipy.interpolate import interp1d
import astropy
import glob

sys.path.append("/home/eesopee/libgetdata2/getdata/bindings/python") #this is needed to import pygetdata

import pygetdata

def count_channels(df):

    '''
    This function counts the number of channels in the dirfile.
    '''

    i_files = glob.glob(df.name+"/chI_*")
    q_files = glob.glob(df.name+"/chQ_*")
    
    channel_number = len(i_files)
    channel_number2 = len(q_files)
    
    if channel_number != channel_number2:
        print("Inconsistent number of I and Q channels!")
    else:
        return channel_number


dirfile = "./dirfile_20230202_154055/"
fitsfile = "./M31_daisy5_1676027931.071732.fits"

discos_hdul = fits.open(fitsfile)
df = pygetdata.dirfile(dirfile) 

nchannels = count_channels(df)
                                                #Here I'm starting the dirfile 5s before the scan.

'''
The acquisition starts slightly before the scan. We search for the first timestamp in the scan to be able to interpolate
'''

'''
Here we define the time axes for the dirfile and the discos fits.
The discos fits gives a timestamp in mjd with a frequency of 5-10 Hz
The dirfile timestamp is in ctime at the sampling frequency (244 Hz)
'''

times = df.getdata("time")
t = Time(times, format='unix', scale='utc')
t.format = 'mjd'

dirfile_time = t.value - 5/86400 - 19389.9681

# *almost* syncronizing the two datasets

############################################################################################################################
# SECTION TABLE
############################################################################################################################



############################################################################################################################
# RF INPUTS
############################################################################################################################



############################################################################################################################
# FEED TABLE
############################################################################################################################

index, px_x, px_y = np.loadtxt("./config_files/pixel_positions_415.csv", unpack=True, delimiter=",", skiprows=1)

'''
Since this is a simulation for now, we randomly extract a number of pixels equal to the channels in the dirfile. 
In general we will have:
- KIDs with scientific data
- Out of resonance channels to monitor the read-out noise 

OoR channels will have no position in the focal plane.

'''

np.random.seed(69)

pxindex = np.random.choice(index, size=nchannels, replace=False)

px_x = [px_x[int(i-1)] for i in pxindex]
px_y = [px_y[int(i-1)] for i in pxindex]

pxindex = np.arange(nchannels)

gains = np.random.normal(loc=1, scale=0.5, size=nchannels)

centers = (np.ones(nchannels), np.ones(nchannels))
radii = np.ones(nchannels)
rotations = np.ones(nchannels)

#fig = plt.figure(figsize=(7,7))
#plt.scatter(px_x, px_y, c=index, marker=",", s=100)


cols = []

idCol = fits.Column(name="id",
                    format="J",
                    array=pxindex,
                    unit="")

xOffsetCol = fits.Column(name="xOffset",
                         format="D", 
                         array=px_x,
                         unit="mm")

yOffsetCol = fits.Column(name="yOffset",
                         format="D", 
                         array=px_y,
                         unit="mm")

relPowCol = fits.Column(name="relativePower",
                        format="D",
                        array=gains,
                        unit="")

cols.append(idCol)
cols.append(xOffsetCol)
cols.append(yOffsetCol)
cols.append(relPowCol)

#centers
#radii
#rotations

#va parsato format_complex. Contiene tutto.

format_complex_extra_filename = df.name + "/format_complex_extra"
format_complex_extra_file = open(format_complex_extra_filename)
format_complex_extra = format_complex_extra_file.readlines()

# centered = [::8] I_center[i] Q_center[i]
# centered and rotated [1::8] cosi[i] + sini[i]
# centered rotated and scaled [2::8]  1/radii[i]


    
    
def parse_format_complex_extra(f, nchannels):
    
    centers = []
    radii = []
    angles = []
    
    for i in range(0, nchannels):
        
        center_I, center_Q = format_complex_extra[0::8][i].split(" ")[7].split(";")
        center_I = float(center_I)
        center_Q = float(center_Q)
       
        center = (center_I, center_Q) 
       
        centers.append(center)
       
        radius = 1/float(format_complex_extra[2::8][i].split(" ")[3])
        
        radii.append(radius)
        
        cosi, sini = format_complex_extra[1::8][i].split(" ")[4].split(";")
        cosi = float(cosi)
        sini = float(sini)
        angle = np.arcsin(sini)
        
        angles.append(angle)
        
    return np.array(centers), np.array(radii), np.array(rotations)

centers, radii, rotations = parse_format_complex_extra(format_complex_extra, nchannels)

centersCol = fits.Column(name="center",
                        format="2D",
                        array=centers,
                        unit="")

radiiCol = fits.Column(name="radius",
                        format="D",
                        array=radii,
                        unit="")

rotationsCol = fits.Column(name="rotation",
                        format="D",
                        array=rotations,
                        unit="")

cols.append(centersCol)
cols.append(radiiCol)
cols.append(rotationsCol)

'''
Qui scriviamo le frequenze di risonanza e se è in o out of resonance. Da prendere dal target_freqs_new o simili
'''

freqs0 = np.linspace(435-256, 435+256, nchannels)
freqs = freqs0 + np.random.uniform(low=-1.5, high=+1.5, size=nchannels) #scatteriamo le frequenze di ~1 MHz

res_flag = np.ones(nchannels)
res_flag[0:10] = 0 #metto a 0 il flag per i toni fuori risonanza


freqsCol = fits.Column("resFreq",
                        format="D",
                        array=freqs,
                        unit="MHz")

resFlagCol = fits.Column("resFlag",
                        format="D",
                        array=res_flag,
                        unit="Bool")

cols.append(freqsCol)
cols.append(resFlagCol)

hdu_feedtable = fits.BinTableHDU.from_columns(cols)

index_new = discos_hdul.index_of("FEED TABLE")
discos_hdul.insert(index=index_new, hdu=hdu_feedtable)
discos_hdul[index_new].header.insert(-1, ("EXTNAME","FEED TABLE"))

############################################################################################################################
# DATA TABLE
############################################################################################################################

discos_time = discos_hdul["DATA TABLE"].data["TIME"]

dt = discos_hdul["DATA TABLE"]
number_of_fields = len(dt.header["TTYPE*"])

cols = []

for index in range(1, number_of_fields+1):
    print(index)
    field_name = dt.header["TTYPE"+str(index)]
    format = dt.header["TFORM"+str(index)]
    
    try:
        unit = dt.header["TUNIT"+str(index)]
    except:
        unit=None
        pass
    
    if field_name != "weather" and field_name != "time":
        print(field_name)
        data = dt.data[field_name]
        data_interpolated = np.interp(x = dirfile_time,
                                      xp = discos_time,
                                      fp = data,
                                      left = np.nan,
                                      right = np.nan)
        
        null_mask = [np.isnan(data_interpolated) == False]
        null_mask = null_mask[0]
        
        field_name_new = field_name + "_interpolated"
        
        col = fits.Column(name=field_name_new,
                          format=format,
                          array=data_interpolated, 
                          unit=unit)
        cols.append(col)

    elif(field_name == "weather"):
        #print("interpolating weather")

        weather_data_new = []
        weather_data = dt.data[field_name]
        weather_data = np.transpose(weather_data)


        for i in range(0,3):
            weather_interp = np.interp(x = dirfile_time,
                          xp = discos_time,
                          fp =weather_data[i],
                          left=np.nan,
                          right=np.nan)
            
            weather_data_new.append(weather_interp)
        
        weather_data_new = np.transpose(weather_data_new)
        
        field_name_new = field_name + "_interpolated"

        col = fits.Column(name=field_name_new, 
                          format=format, 
                          array=weather_data_new, 
                          unit=unit)
        cols.append(col)

    elif(field_name == "time"):

        time_new = dirfile_time
        field_name_new = field_name + "_interpolated"

        col = fits.Column(name=field_name_new, 
                          format=format, 
                          array=time_new, 
                          unit=unit)
        cols.append(col)

hdu2 = fits.BinTableHDU.from_columns(cols)

index_new = discos_hdul.index_of("DATA TABLE")+1
discos_hdul.insert(index=index_new, hdu=hdu2)
discos_hdul[index_new].header.insert(-1, ("EXTNAME","DATA TABLE INTERP"))



##################################################################################
#
#       WRITING THE NEW FITS FILE
#
##################################################################################



############################################################################################################################
# IQ TABLE
############################################################################################################################

def count_channels(df):

    '''
    This function counts the number of channels in the dirfile.
    '''

    i_files = glob.glob(df.name+"/chI_*")
    q_files = glob.glob(df.name+"/chQ_*")
    
    channel_number = len(i_files)
    channel_number2 = len(q_files)
    
    if channel_number != channel_number2:
        print("Inconsistent number of I and Q channels!")
    else:
        return channel_number

def get_tod(df):
    
    number_of_channels = count_channels(df)
    
    entries = [["ch{:s}_{:03d}".format(I_or_Q,ch) for ch in range(0,number_of_channels)] for I_or_Q in ("I", "Q")]
    
    tod_I = [df.getdata(entry) for entry in entries[0]]
    tod_Q = [df.getdata(entry) for entry in entries[1]]

    entries = np.concatenate((entries[0], entries[1]))

    return tod_I, tod_Q, entries

entry_list = df.entry_list()

print("Retrieving TODs from dirfile")
I, Q, entries = get_tod(df)

I = np.transpose(np.array(I, dtype='int32'))
Q = np.transpose(np.array(Q, dtype='int32'))

print("I table shape =",I.shape)
print("Q table shape =",Q.shape) 

print("Writing TOD hdu to ANTENNA TEMP TABLE")

iq_table = np.array(np.concatenate((I,Q),axis=1), dtype=np.int32) #joining I and Q tables into a single table
table = astropy.table.QTable(iq_table, names=entries)
TOD_hdu = fits.BinTableHDU(data=table)

discos_hdul.pop(index=5) #removes antenna temp table
discos_hdul.insert(index=5, hdu=TOD_hdu)
discos_hdul[5].header.insert(-1, ("EXTNAME","IQ TABLE"))


############################################################################################################################
# SERVO TABLE
############################################################################################################################

############################################################################################################################
# HOUSEKEEPING TABLE
############################################################################################################################



discos_hdul.writeto('output_test.fits', overwrite=True)
