#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 10:24:40 2023

@author: eesopee
"""

'''
This script simulates a scan with SRT.

It generates a FITS file with information about pointing of the telescope and other ancilliary data.
It then starts an acquisition with the read-out system of MISTRAL or a simulated version of it.
Finally, it merges the dirfile with the fits file in a usable fits.

'''

class RoachCommands:
    
    '''
    This class runs the read-out commands that should be sent by discos via ethernet.
    Here is just a python wrapper.
    '''
    
    def __init__(self, clientpath):
        self.clientpath = clientpath
    
    def initialize(self):
        
        return
    
    def target_sweep(self):
        
        return
    
    def start(self):
        
        return
    
    def stop(self):
        
        return

import gen_discos_fits as discosfits
import time

ro = RoachCommands("/home/mistral/src/mistral_simulator/mistral_client_simulator.py")

ro.start()

duration=60

discosfits.gen_fits(target="M31",
                    t0="now",
                    scan="cross_radec",
                    radius=1,
                    duration=duration,
                    outname="simulator_test_aaa")

#time.wait(duration)

ro.stop()