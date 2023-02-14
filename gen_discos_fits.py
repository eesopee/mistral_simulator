#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 16:48:29 2022

@author: eesopee
"""

'''
This script generates a fits file in the discos format with an user-defined scan curve. It also updates it with the correct
feed table with the MISTRAL array.
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
import astroplan

time_start = time.time()

def gen_fits(target, 
             scan, 
             site = "SRT", 
             t0 = "now",
             duration=60, 
             radius=1, 
             outname="output"):
    
    loc = EarthLocation.of_site(site) #retrieving lon-lat of SRT to calculate AZ EL and vice versa
    discos_sampling_frequency = 10
     
    observer = astroplan.Observer(loc)
    
    script_start = time.time()

    if t0 == "now":
        start_time = time.time()
        t0 = astropy.time.Time(float(start_time), format="unix").mjd #59200.84973159712
    else:
        start_time = t0 #t0 should be in MJD
        
    target = SkyCoord.from_name(str(target), frame='icrs')
 
    radius = float(radius)
    radius /= 60
 
    duration = float(duration) / 86400 #scan duration in days
    
    if scan == "daisy":    
        print("Daisy")
        t,ra,dec,alt,az = daisy_scan(target=target, 
                                     t0=t0, 
                                     duration=duration, 
                                     tau_radial=duration/22, 
                                     radius=radius, 
                                     loc=loc,
                                     discos_sampling_frequency=discos_sampling_frequency)
    
    if scan == "cross_radec":
        print("radec")
        t,ra,dec,alt,az = radec_cross_scan(target, 
                                           span=radius, 
                                           t0=t0, 
                                           duration=duration, 
                                           loc=loc,
                                           discos_sampling_frequency=discos_sampling_frequency)
    
    if scan == "ra_subscan":
        print("ra_subscan")
        t,ra,dec,alt,az = ra_subscan(target, 
                                           span=radius, 
                                           t0=t0, 
                                           duration=duration, 
                                           loc=loc,
                                           discos_sampling_frequency=discos_sampling_frequency)

    
    if np.any(alt < 0):
        print("WARNING: target below horizon")
    if np.any(alt < 30) :
        print("WARNING: target at low elevation")
        print("min EL=", np.min(alt))
    if np.any(alt > 80):
        print("WARNING: target at high elevation")
    else:
        pass
    
    time_axis = astropy.time.Time(t, format="mjd")
    
    par_angle = observer.parallactic_angle(time_axis, target=target) #np.ones(len(t))*(-9999)
    derot_angle = np.ones(len(t))*(-9999)
    flag_cal = np.ones(len(t))
    flag_track=np.ones(len(t))
    weather = np.ones(shape=(len(t),3))

    discos_hdul = fits.open("./discos_template.fits")
    discos_dt_hdu = discos_hdul["DATA TABLE"]
    discos_dt = discos_dt_hdu.data
 
        
    '''
    data table format:
        
    t [MJD], raJ2000 [rad], decJ2000 [rad], az [rad], el [rad], par_angle [rad], derot_angle [rad], flag_cal [0/1], flag_track [0/1], weather (a,b,c)
    
    We don't need all entries. t, ra, dec, alt, ax are sufficient. We don't need par_angle too.
    '''
     
    datas = (t, ra, dec, az, alt, par_angle, derot_angle, flag_cal, flag_track, weather)
 
    cols = []
 
    for field in discos_dt_hdu.header["TTYPE*"]:
        index = field.split("TTYPE")[-1]
        label = discos_dt_hdu.header[field]
        
        fmt = discos_dt_hdu.header["TFORM"+str(index)]
        
        index = int(index)
        
        try:
            unit = discos_dt_hdu.header["TUNIT"+str(index)]
            col = fits.Column(name=label, format=fmt, array=datas[index-1], unit=unit)
        except:
            col = fits.Column(name=label, format=fmt, array=datas[index-1])
            
        cols.append(col)
 
    hdu_new = fits.BinTableHDU.from_columns(cols)
 
    discos_hdul["DATA TABLE"] = hdu_new
    index = discos_hdul.index_of("FEED TABLE")+1
    discos_hdul[index].header.insert(-1, ("EXTNAME","DATA TABLE"))
    discos_hdul.writeto(outname+'_{:f}.fits'.format(start_time), overwrite=True)
 
    print("Done in ", time.time()-script_start,"s")


def radec_to_altaz(t, coords, loc):
    obstime = astropy.time.Time(t, format="mjd")
    aa = astropy.coordinates.AltAz(location=loc, obstime=obstime)
    coords_new = coords.transform_to(aa)
    
    return coords_new


def radec_cross_scan(target, span, t0, duration, include_overhead=False, loc=None,discos_sampling_frequency=10):

    overhead = 5 / 86400  #time to repoint the antenna after the first scan. 
    
    th = np.linspace(t0, t0+duration/2, int(duration/2*86400*discos_sampling_frequency)) #time array for the horizontal branch
    to = np.linspace(th[-1], th[-1]+overhead,int(overhead*86400*discos_sampling_frequency)) #time for the overhead
    tv = np.linspace(to[-1], to[-1]+duration/2, int(duration/2*86400*discos_sampling_frequency)) #time for the vertical branch
    
    ra_h = np.linspace(target.ra.deg-span/2, target.ra.deg+span/2, len(th))
    dec_h = np.ones(len(th))*target.dec.deg 
    ra_v = np.ones(len(tv))*target.ra.deg
    dec_v = np.linspace(target.dec.deg-span/2, target.dec.deg+span/2, len(tv))

    t = np.concatenate((th,tv))
    ra = np.concatenate((ra_h, ra_v))
    dec = np.concatenate((dec_h, dec_v))

    coords = SkyCoord(ra = ra, dec = dec, unit="deg")
    
    coords_new = radec_to_altaz(t, coords, loc)
    
    alt = coords_new.alt
    az = coords_new.az
    
    return t, ra,dec,alt.value,az.value

def ra_subscan(target, span, t0, duration, include_overhead=False, loc=None, discos_sampling_frequency=10):

    overhead = 5 / 86400  #time to repoint the antenna after the first scan. 
    
    t = np.linspace(t0, t0+duration, int(duration*86400*discos_sampling_frequency)) #time array for the horizontal branch
    
    ra = np.linspace(target.ra.deg-span/2, target.ra.deg+span/2, len(t))
    dec = np.ones(len(t))*target.dec.deg 

    coords = SkyCoord(ra = ra, dec = dec, unit="deg")
    
    coords_new = radec_to_altaz(t, coords, loc)
    
    alt = coords_new.alt
    az = coords_new.az
    
    return t, ra,dec,alt.value,az.value


def daisy_scan(target, radius, t0, duration, tau_radial, loc, discos_sampling_frequency=10):
    
    t = np.linspace(t0, t0+duration, int(duration*86400*discos_sampling_frequency))
    
    puls = 2*np.pi/tau_radial
    alpha = 1
    beta= 1
    phi_1=0
    phi_2=0
    r_0 = radius

    x_0 = target.ra.deg
    y_0 = target.dec.deg 
    
    x = x_0 + r_0*np.sin(puls*alpha*t+phi_1)*np.cos(beta*puls*t/np.pi+phi_2)/np.cos(y_0)
    y = y_0 + r_0*np.sin(puls*alpha*t+phi_1)* np.sin(beta*puls*t/np.pi+phi_2)
    
    coords = SkyCoord(ra=x, dec=y, unit="deg")
    
    coords_new = radec_to_altaz(t,coords,loc)
    
    return t, x,y,coords_new.alt.value,coords_new.az.value


'''
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", help ="Name of the target")
    parser.add_argument("--scan", help="Scan curve. Options: cross_radec, daisy")
    parser.add_argument("--site", default = "SRT", help="Observatory. Default: SRT")
    parser.add_argument("--t0", default = time.time(), help="Scan start time. Default: now")
    parser.add_argument("--duration", default = 60, help="Scan duration in seconds. Default: 60s")
    parser.add_argument("--radius", default = 1, help="Scan radius in arcmin. Default: 1arcmin")
    parser.add_argument("--outname", default="output", help="Name of the output filename. Default: Output")
    
    args = parser.parse_args()
    
    
    gen_fits(args.target,
             args.scan,
             args.site,
             args.t0,
             args.duration,
             args.radius,
             args.outname)
'''  