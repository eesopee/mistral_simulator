#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 16 13:05:55 2021

@author: eesopee
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("/home/eesopee/libgetdata2/getdata/bindings/python") #this is needed to import pygetdata

import pygetdata

def count_channels(df):
    
    channels_splitted = [entry.decode().split("chp_")for entry in entry_list[2:]] #searching for chp entry
    #only entries starting with chp will have size > 1
    len_mask = np.array([len(chitem) for chitem in channels_splitted])
    mask = np.array(len_mask) > 1
    #here I select the arrays with len > 1    
    channels = [int(channels_splitted[i][-1]) for i in range(len(channels_splitted)) if mask[i]]
    #then I retrieve the maximum number
    channel_number = np.max(channels)

    return channel_number


df = pygetdata.dirfile("./dirfile_20230202_154055") #this generates a dirfile object

entry_list = df.entry_list() #this prints the list of ENTRIES in the dirfile

#entry = entry_list[5] #taking a random entry
#data = df.getdata(entry) #this retrieves the data for a given entry

'''
First we determine how many channels we have
'''

time = df.getdata(entry_list[0])
print(time)
index = df.getdata(entry_list[1])
    
count_channels(df)

def plot_channels(df):
    
    number_of_channels = count_channels(df)
    print(number_of_channels)
    
    entries = ["chp_{:03d}".format(ch) for ch in range(0,number_of_channels)]
    print(entries)    
    
    [plt.plot(time, df.getdata(entry)) for entry in entries]
    
def retrieve_phases(df):
    
    number_of_channels = count_channels(df)
    entries = ["chp_{:03d}".format(ch) for ch in range(0,number_of_channels)]
    phases = [df.getdata(entry) for entry in entries]

    phases = np.transpose(phases)
    return phases

phases_table = np.transpose(np.array(retrieve_phases(df))) #each row is a sample. Columns are channels.

'''

this converts olimpo pc time to mjd
'''
from astropy.time import Time
times = df.getdata("time")
t = Time(times, format='unix', scale='utc')
t.format = 'mjd'
t = t.value