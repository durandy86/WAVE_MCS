#!/usr/bin/env python
# coding: utf-8

import os, sys

import numpy as np
import h5netcdf

import xarray as xr
import xarray.ufuncs as xu
import xrft
import pandas as pd

from waveFiltering import *
###########################################################################
###########################################################################

units = "m.s^{-1}"

def filterWave(ds, tcwvhat, var, freq_lon_Save, freq_time_Save):
    # for j in range(0,ds.lat.size,50):
    #    _ds = ds.isel(lat = slice(j,j+50))
    _ds = ds
    # _tcwvhat = tcwvhat.isel(lat = slice(j,j+50)).persist()
    _tcwvhat = tcwvhat

    ### Delete LF
    _dahat_LF, da_LFfilter = mjofilter(_tcwvhat, ds, fmin = -1/96, fmax= -1/999, kmax = 10, kmin = -10, freq_lon_Save = freq_lon_Save, freq_time_Save = freq_time_Save)
    _tcwvhat = _tcwvhat - _dahat_LF.fillna(0.)
    ###

    
    #####################################
    # Kelvin Filtering
    # kelvin: (Straub Kiladis, Kiladis 2009)
    tMin = 2.5
    tMax = 30.
    kMin = 1
    kMax = 14
    hMin = 8
    hMax = 90


    _dahat_Kev, da_Kevfilter = kelvinfilter(_tcwvhat, ds, fmin = -1/tMin, fmax = -1/tMax, kmax = kMax, kmin = 1, hmin=8, hmax=90, freq_lon_Save = freq_lon_Save, freq_time_Save = freq_time_Save)
    da_Kevfilter = deletebug(da_Kevfilter, _ds)
    
    da_Kevfilter.attrs['units'] = units
    da_Kevfilter.attrs['tMin'] = str(tMin)
    da_Kevfilter.attrs['tMax'] = str(tMax)
    da_Kevfilter.attrs['kMin'] = str(kMin)
    da_Kevfilter.attrs['kMax'] = str(kMax)
    da_Kevfilter.attrs['hMin'] = str(hMin)
    da_Kevfilter.attrs['hMax'] = str(hMax)
    
    ######################################
    # Rossby Filtering
    # n=1 Eq. Rossby  Kiladis et al 2009
    tMin = 9.7
    tMax = 72.
    kMin = -10
    kMax = -1
    hMin = 8 
    hMax = 90

    _dahat_Ross, da_Rossfilter = erfilter(_tcwvhat, ds, fmin = -1/tMin, fmax = -1/tMax, kmax = kMax, kmin = kMin, hmin=8, hmax=90,freq_lon_Save = freq_lon_Save, freq_time_Save = freq_time_Save)
    da_Rossfilter = deletebug(da_Rossfilter, _ds)
    
    da_Rossfilter.attrs['units'] = units
    da_Rossfilter.attrs['tMin'] = str(tMin)
    da_Rossfilter.attrs['tMax'] = str(tMax)
    da_Rossfilter.attrs['kMin'] = str(kMin)
    da_Rossfilter.attrs['kMax'] = str(kMax)
    da_Rossfilter.attrs['hMin'] = str(hMin)
    da_Rossfilter.attrs['hMax'] = str(hMax)
    
    #########################################
    # MJO Filtering
    # MJO
    # waveName = "MJO"   ; kiladis et al. 2009
    tMin = 30.
    tMax = 96.
    kMin = 1
    kMax = 5  # 5 is also possible

    _dahat_MJO, da_MJOfilter = mjofilter(_tcwvhat, ds, fmin = -1/tMin, fmax= -1/tMax, kmax = kMax, kmin = kMin,freq_lon_Save = freq_lon_Save, freq_time_Save = freq_time_Save)
    da_MJOfilter = deletebug(da_MJOfilter, _ds)

    da_MJOfilter.attrs['units'] = units
    da_MJOfilter.attrs['tMin'] = str(tMin)
    da_MJOfilter.attrs['tMax'] = str(tMax)
    da_MJOfilter.attrs['kMin'] = str(kMin)
    da_MJOfilter.attrs['kMax'] = str(kMax)

    #############################################
    # MRG Filtering
    # MRG
    tMin = 3.
    tMax = 8.
    kMin = -10
    kMax = -1
    hMin = 8
    hMax = 90

    _dahat_MRG, da_MRGfilter = mrgfilter(_tcwvhat, ds, fmin = -1/3, fmax = -1/8, kmax = -1, kmin = -10, hmin = 8, hmax = 90, freq_lon_Save = freq_lon_Save, freq_time_Save = freq_time_Save)
    da_MRGfilter = deletebug(da_MRGfilter, _ds)

    da_MRGfilter.attrs['units'] = units
    da_MRGfilter.attrs['tMin'] = str(tMin)
    da_MRGfilter.attrs['tMax'] = str(tMax)
    da_MRGfilter.attrs['kMin'] = str(kMin)
    da_MRGfilter.attrs['kMax'] = str(kMax)
    da_MRGfilter.attrs['hMin'] = str(hMin)
    da_MRGfilter.attrs['hMax'] = str(hMax)

    ################################################
    # IGW Filtering
    # EIG  (Wheeler Kiladis 1998)    
    tMin = 1.
    tMax = 5.
    kMin = 14
    kMax = 1
    hMin = 12  
    hMax = 50
   
    _dahat_EIG, da_EIGfilter = igfilter(_tcwvhat, ds, fmin = -1, fmax = -1/5, kmin=1, kmax=14, hmin=12, freq_lon_Save = freq_lon_Save, freq_time_Save = freq_time_Save)
    da_EIGfilter = deletebug(da_EIGfilter, _ds)

    da_EIGfilter.attrs['units'] = units
    da_EIGfilter.attrs['tMin'] = str(tMin)
    da_EIGfilter.attrs['tMax'] = str(tMax)
    da_EIGfilter.attrs['kMin'] = str(kMin)
    da_EIGfilter.attrs['kMax'] = str(kMax)
    da_EIGfilter.attrs['hMin'] = str(hMin)
    da_EIGfilter.attrs['hMax'] = str(hMax)

    # WIG  (Wheeler Kiladis 1998)
    tMin = 1.
    tMax = 5.
    kMin = -14
    kMax = -1
    hMin = 12 
    hMax = 50
    
    _dahat_WIG, da_WIGfilter = igfilter(_tcwvhat, ds, fmin = -1/tMin, fmax = -1/tMax, kmin=kMin, kmax=kMax, hmin=12, hmax=50, freq_lon_Save = freq_lon_Save, freq_time_Save = freq_time_Save)
    da_WIGfilter = deletebug(da_WIGfilter, _ds)
    
    da_WIGfilter.attrs['units'] = units
    da_WIGfilter.attrs['tMin'] = str(tMin)
    da_WIGfilter.attrs['tMax'] = str(tMax)
    da_WIGfilter.attrs['kMin'] = str(kMin)
    da_WIGfilter.attrs['kMax'] = str(kMax)
    da_WIGfilter.attrs['hMin'] = str(hMin)
    da_WIGfilter.attrs['hMax'] = str(hMax)
    ######################################################
    # TD filtering
    
    # TD (Easterly wave) : 
    tMin = 2.5
    tMax = 7.
    kMin = -20
    kMax = -6
    # hMin = mis  ==> Min curve = up curve of MRG Domain
    # hMax = mis
    _dahat_TD, da_TDfilter = mrgfilter(_tcwvhat, ds, fmin = -1/2.5, fmax = -1/7, kmax = -6, kmin = -20, hmin = 100, freq_lon_Save = freq_lon_Save, freq_time_Save = freq_time_Save)
    da_TDfilter = deletebug(da_TDfilter, _ds)
    
    da_TDfilter.attrs['units'] = units
    da_TDfilter.attrs['tMin'] = str(tMin)
    da_TDfilter.attrs['tMax'] = str(tMax)
    da_TDfilter.attrs['kMin'] = str(kMin)
    da_TDfilter.attrs['kMax'] = str(kMax)
    """
    ########################################################
    
    # Low Frequency Filtering
    # Low frequency 
    tMin = 120.
    tMax = 999
    kMin = -10
    kMax = 10
    
    _dahat_LF, da_LFfilter = mjofilter(_tcwvhat, ds, fmin = -1/120, fmax= -1/999, kmax = 10, kmin = -10, freq_lon_Save = freq_lon_Save, freq_time_Save = freq_time_Save)
    da_LFfilter = deletebug(da_LFfilter, _ds)

    da_LFfilter.attrs['units'] = units
    da_LFfilter.attrs['tMin'] = str(tMin)
    da_LFfilter.attrs['tMax'] = str(tMax)
    da_LFfilter.attrs['kMin'] = str(kMin)
    da_LFfilter.attrs['kMax'] = str(kMax)
    """


    #######################################################
    #######################################################
    # NetCDF save
    ds_NetCdf = da_Kevfilter.to_dataset(name = var + '_Kelvin')
    ds_NetCdf[var + '_Rossby'] = da_Rossfilter
    ds_NetCdf[var + '_MJO'] = da_MJOfilter
    ds_NetCdf[var + '_MRG'] = da_MRGfilter
    ds_NetCdf[var + '_EIG'] = da_EIGfilter
    ds_NetCdf[var + '_WIG'] = da_WIGfilter
    ds_NetCdf[var + '_TD'] = da_TDfilter
    # ds_NetCdf[var + '_LF'] = da_LFfilter 

    ds_NetCdf.time.encoding['units'] = "hours since 1900-01-01"
    ds_NetCdf.time.encoding['calendar'] = "proleptic_gregorian"
    ds_NetCdf.time.encoding['standard_name'] = "time"
    return ds_NetCdf

