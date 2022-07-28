#!/usr/bin/env python
# coding: utf-8

import numpy as np
import h5netcdf

import xarray as xr
import xarray.ufuncs as xu
import xrft
import pandas as pd

# Constante
rlat = 0
pi    = np.pi
radius = 6.37122e06    # [m]   average radius of earth
g     = 9.80665        # [m/s] gravity at 45 deg lat used by the WMO
omega = 7.292e-05      # [1/s] earth's angular vel
ll    = 2.*pi*radius*np.cos(np.abs(rlat))
Beta  = 2.*omega*np.cos(np.abs(rlat))/radius
fillval = 1e20

#############################################
#############################################
def deletebug(da, ds):
    _da = da.real.values + da.conj().real.values
    da = xr.DataArray(_da,
                       dims=("lat","time","lon"), 
                       coords={"lat":ds.lat,
                               "time":ds.time,
                               "lon":ds.lon

                               })
    # da = da.sel(time = str(f))
    # da['time'] = date_Re
      
    return da

#############################################
#############################################


def kelvinfilter(da, ds, fmin=None, fmax=None, kmin=None, kmax=None, hmin=None, hmax=None, rlat = 0, freq_lon_Save = None, freq_time_Save = None):
    """kelvin wave filter
    Arguments:
       'fmin/fmax' -- unit is cycle per day
       'kmin/kmax' -- zonal wave number
       'hmin/hmax' --equivalent depth
    """
    knum = da.freq_lon
    freq = da.freq_time
    
    # filtering ############################################################
    mask = da.copy()
    #wavenumber cut-off
    if kmin != None:
        mask = mask.where((da.freq_lon > kmin) & (da.freq_lon < kmax), drop = False)
    
    #frequency cutoff
    if fmin != None:
        mask = mask.where((da.freq_time > fmin) & (da.freq_time < fmax), drop = False)

    #dispersion filter
    if hmin != None:
        c = np.sqrt(g*hmin)
        omega = -2.*pi*freq/24./3600. / np.sqrt(Beta*c) #adusting day^-1 to s^-1
        k     = 1/((radius)/(knum-1)) * np.sqrt(c/Beta)         #adusting wavenumber to m
        mask = mask.where(lambda da: -(omega - k) <= 0., drop = False)
    if hmax != None:
        c = np.sqrt(g*hmax)
        omega = -2.*pi*freq/24./3600. / np.sqrt(Beta*c) #adusting day^-1 to s^-1
        k     = 1/((radius)/(knum+1)) * np.sqrt(c/Beta)         #adusting wavenumber to m
        mask = mask.where(lambda da: -(omega - k) >= 0., drop = False)

    _varhat = mask.persist()
    _varhat['freq_lon'] = freq_lon_Save
    _varhat['freq_time'] = freq_time_Save
    filterd = xrft.ifft(_varhat.fillna(0.), dim = ['freq_time','freq_lon'], 
                                  true_phase=False, true_amplitude=True) # Signal in direct space
    filterd['time'] = ds['time']
    return mask, filterd

def mjofilter(da, ds, fmin=None, fmax=None, kmin=None, kmax=None, freq_lon_Save = None, freq_time_Save = None) :
    """kelvin wave filter
    Arguments:
       'fmin/fmax' -- unit is cycle per day
       'kmin/kmax' -- zonal wave number
       'hmin/hmax' --equivalent depth
    """
    
    knum = da.freq_lon
    freq = da.freq_time

    # filtering ############################################################
    mask = da.copy()
    #wavenumber cut-off
    if kmin != None:
        mask = mask.where((da.freq_lon > kmin) & (da.freq_lon < kmax), drop = False)
    
    #frequency cutoff
    if fmin != None:
        mask = mask.where((da.freq_time > fmin) & (da.freq_time < fmax), drop = False)

    _varhat = mask.persist()
    _varhat['freq_lon'] = freq_lon_Save
    _varhat['freq_time'] = freq_time_Save
    filterd = xrft.ifft(_varhat.fillna(0.), dim = ['freq_time','freq_lon'], 
                                  true_phase=False, true_amplitude=True)  # Signal in direct space
    filterd['time'] = ds['time']
    return mask, filterd

def erfilter(da, ds, fmin=None, fmax=None, kmin=None, kmax=None, hmin=None, hmax=None, n=1, rlat = 0, freq_lon_Save = None, freq_time_Save = None):
    """equatorial wave filter
    Arguments:
       'fmin/fmax' -- unit is cycle per day
       'kmin/kmax' -- zonal wave number
       'hmin/hmax' -- equivalent depth
       'n'         -- meridional mode number
    """

    if n <=0 or n%1 !=0:
        print ("n must be n>=1 integer")
        sys.exit()

    knum = da.freq_lon
    freq = da.freq_time

    
    # filtering ############################################################
    mask = da.copy()
    #wavenumber cut-off
    if kmin != None:
        mask = mask.where((da.freq_lon > kmin) & (da.freq_lon < kmax), drop = False)
    
    #frequency cutoff
    if fmin != None:
        mask = mask.where((da.freq_time > fmin) & (da.freq_time < fmax), drop = False)
        
        
    # filtering ############################################################
    if hmin != None:
        c = np.sqrt(g*hmin)
        omega = -2.*pi*freq/24./3600. / np.sqrt(Beta*c)           #adusting day^-1 to s^-1
        k     = 1/((radius)/(knum+1)) * np.sqrt(c/Beta)                          #adusting ^2pia to ^m               
        mask = mask.where(lambda da: -(omega*(k**2 + (2*n+1)) + k) <= 0., drop = False)   
    if hmax != None:
        c = np.sqrt(g*hmax)
        omega = -2.*pi*freq/24./3600. / np.sqrt(Beta*c)           #adusting day^-1 to s^-1
        k     = 1/((radius)/(knum-1)) * np.sqrt(c/Beta)   
        mask = mask.where(lambda da: -(omega*(k**2 + (2*n+1)) + k) >= 0., drop = False)             

    _varhat = mask.persist()
    _varhat['freq_lon'] = freq_lon_Save
    _varhat['freq_time'] = freq_time_Save
    filterd = xrft.ifft(_varhat.fillna(0.), dim = ['freq_time','freq_lon'], 
                                  true_phase=False, true_amplitude=True) # Signal in direct space
    filterd['time'] = ds['time']
    return mask, filterd

def mrgfilter(da, ds, fmin=None, fmax=None, kmin=None, kmax=None, hmin=None, hmax=None, freq_lon_Save = None, freq_time_Save = None):
    """mixed Rossby gravity wave
    Arguments:
       'fmin/fmax' -- unit is cycle per day
       'kmin/kmax' -- zonal wave number. negative is westward, positive is
                      eastward
       'hmin/hmax' -- equivalent depth
    """
    
    knum = da.freq_lon
    freq = da.freq_time
    
        
    # filtering ############################################################
    mask = da.copy()
    #wavenumber cut-off
    if kmin != None:
        mask = mask.where((da.freq_lon > kmin) & (da.freq_lon < kmax), drop = False)
    
    #frequency cutoff
    if fmin != None:
        mask = mask.where((da.freq_time > fmin) & (da.freq_time < fmax), drop = False)
        
        
    #dispersion filter
    mask = mask.where((da.freq_lon <= -1.), drop = False)
    if hmin != None:
        c = np.sqrt(g*hmin)
        omega = -2.*pi*freq/24./3600. / np.sqrt(Beta*c)           #adusting day^-1 to s^-1
        k     = 1/((radius)/(knum-1)) * np.sqrt(c/Beta)                         #adusting ^2pia to ^m               
        mask = mask.where(lambda da: -(omega**2 - k*omega - 1) <= 0., drop = False) 
    if hmax != None:
        c = np.sqrt(g*hmax)
        omega = -2.*pi*freq/24./3600. / np.sqrt(Beta*c)           #adusting day^-1 to s^-1
        k     = 1/((radius)/(knum+1)) * np.sqrt(c/Beta)                 #adusting ^2pia to ^m               
        mask = mask.where(lambda da: -(omega**2 - k*omega - 1) >= 0., drop = False) 
  
    _varhat = mask.persist()
    _varhat['freq_lon'] = freq_lon_Save
    _varhat['freq_time'] = freq_time_Save
    filterd = xrft.ifft(_varhat.fillna(0.), dim = ['freq_time','freq_lon'], 
                                  true_phase=False, true_amplitude=True) # Signal in direct space
    filterd['time'] = ds['time']
    return mask, filterd

def igfilter(da, ds, fmin=None, fmax=None, kmin=None, kmax=None, hmin=None, hmax=None, n=1,freq_lon_Save = None, freq_time_Save = None):
    """n>=1 inertio gravirt wave filter. default is n=1 WIG.
    Arguments:
       'fmin/fmax' -- unit is cycle per day
       'kmin/kmax' -- zonal wave number. negative is westward, positive is
                      eastward
       'hmin/hmax' -- equivalent depth
       'n'         -- meridional mode number
    """

    # filtering ############################################################
    mask = da.copy()
    #wavenumber cut-off
    if kmin != None:
        mask = mask.where((da.freq_lon > kmin) & (da.freq_lon < kmax), drop = False)
    
    #frequency cutoff
    if fmin != None:
        mask = mask.where((da.freq_time > fmin) & (da.freq_time < fmax), drop = False)
        
    knum = da.freq_lon
    freq = da.freq_time
    #dispersion filter
    if hmin != None:
        c = np.sqrt(g*hmin)
        omega = -2.*pi*freq/24./3600. / np.sqrt(Beta*c)           #adusting day^-1 to s^-1
        k     = 1/((radius)/(knum-1)) * np.sqrt(c/Beta)                   #adusting ^2pia to ^m               
        mask = mask.where(lambda da: -(omega**2 - k**2 - (2*n+1)) <= 0., drop = False) 
    if hmax != None:
        c = np.sqrt(g*hmax)
        omega = -2.*pi*freq/24./3600. / np.sqrt(Beta*c)           #adusting day^-1 to s^-1
        k     = 1/((radius)/(knum-1)) * np.sqrt(c/Beta)                      #adusting ^2pia to ^m               
        mask = mask.where(lambda da: -(omega**2 - k**2 - (2*n+1)) >= 0., drop = False) 
        
        
    _varhat = mask.persist()
    _varhat['freq_lon'] = freq_lon_Save
    _varhat['freq_time'] = freq_time_Save
    filterd = xrft.ifft(_varhat.fillna(0.), dim = ['freq_time','freq_lon'], 
                                  true_phase=False, true_amplitude=True) # Signal in direct space
    filterd['time'] = ds['time']
    return mask, filterd

