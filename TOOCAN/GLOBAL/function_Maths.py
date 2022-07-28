#!/usr/bin/env python
# coding: utf-8

import os, sys

import numpy as np
import xarray as xr
import pandas as pd

import math
####################
####################

def compute_distribution_init(da, ds_T, data_MCSselection):
    nb_MCS = 0
    distri_init = np.zeros(32, dtype = np.int32)
    surface = np.zeros(32, dtype = np.float32)
    bMin = np.zeros(32, dtype = np.float32)
    for iMCS in np.arange(0,len(data_MCSselection),1):
        nb_phase = da.sel(time = ds_T.time.isel(time = iMCS),    
                 lat = data_MCSselection[iMCS].latInit,
                 lon = data_MCSselection[iMCS].lonInit,
                 method = 'nearest').values

        if math.isnan(nb_phase) != True:
            nb_phase = int(nb_phase)
            distri_init[nb_phase] =  distri_init[nb_phase] + 1
            surface[nb_phase] = surface[nb_phase] + data_MCSselection[iMCS].surfmaxkm2_235K
            bMin[nb_phase] = bMin[nb_phase] + data_MCSselection[iMCS].tbmin
            nb_MCS = nb_MCS + 1

    return distri_init, nb_MCS,  surface/distri_init, bMin/distri_init


####################
####################
def compute_distribution_end(da, ds_T, data_MCSselection):
    nb_MCS = 0
    distri_init = np.zeros(32, dtype = np.int)
    for iMCS in np.arange(0,ds_T.time.size,1):
        nb_phase = da.sel(time = (ds_T.time.isel(time = iMCS)), 
                 lat = data_MCSselection[iMCS].latEnd,
                 lon = data_MCSselection[iMCS].lonEnd,
                 method = 'nearest').values

        if math.isnan(nb_phase) != True:
            nb_phase = int(nb_phase)
            distri_init[nb_phase] =  distri_init[nb_phase] + 1
            nb_MCS = nb_MCS + 1

    return distri_init, nb_MCS

####################
def compute_distribution_life(da, ds_T, data_MCSselection): 
    nb_MCS = 0
    distri_init = np.zeros(32, dtype = np.int)
    for iMCS in np.arange(0,len(data_MCSselection),1):
        for ilife in range( (int(data_MCSselection[iMCS].nbslot))) :
            _time = ds_T.isel(time = iMCS).time.values + np.timedelta64((30*ilife), 'm')
            nb_phase = da.sel(time = _time, 
                     lat = data_MCSselection[iMCS].clusters.lat[ilife],
                     lon = data_MCSselection[iMCS].clusters.lon[ilife],
                     method = 'nearest').values

        if math.isnan(nb_phase) != True:
            nb_phase = int(nb_phase)
            distri_init[nb_phase] =  distri_init[nb_phase] + 1
            nb_MCS = nb_MCS + 1

    return distri_init, nb_MCS

####################
####################
####################


