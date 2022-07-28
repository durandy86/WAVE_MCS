#!/usr/bin/env python
# coding: utf-8

import os, sys
os.environ["MALLOC_TRIM_THRESHOLD_"] = '0'

import numpy as np
import xarray as xr
import pandas as pd
import glob

from Load_TOOCAN import *
####################
####################
path = '/cnrm/tropics/commun/DATASAT/NO_SAVE/TRACKING_LMD/TOOCAN_v2.07/PostProcessing/'

year = 2016
month = 5
day = 15

def load_Data_TOOCAN(iregion):
    path = '/cnrm/tropics/commun/DATASAT/NO_SAVE/TRACKING_LMD/TOOCAN_v2.07/PostProcessing/'
    if(iregion == 0):
            REGION      = 'AFRICA'
    if(iregion == 1):
            REGION      = 'INDIA'
    #
    # if Region == WESTERNPACIFIC and the date is before 20150601
    # then the satellite used is MTSAT-2
    # For technical reasons inherent to the Geostationary satellite, the SouthScan is not processed
    ###############################################################################################
    if( (iregion == 2) & (year*100+month < 201506)) :
            REGION      = 'WESTERNPACIFIC'
    #
    # if Region == WESTERNPACIFIC and the date is after 20150601
    # then the satellite used is HIMAWARI-8
    ###############################################################
    if( (iregion == 2) & (year*100+month >= 201506)) :
            REGION      = 'WESTERNPACIFIC'
    if(iregion == 3):
            REGION      = 'EASTERNPACIFIC'
    if(iregion == 4):
            REGION      = 'AMERICA'
    print(REGION)
    ##########################################################
    # Load the MCS parameters into a class object
    # ##########################################################
    All_FileTOOCAN = sorted(glob.glob(path+REGION+ '/*.dat'))

    for ifileTOOCAN in np.arange(0,len(All_FileTOOCAN) - (12*3) ,1):
        FileTOOCAN = All_FileTOOCAN[ifileTOOCAN]
        print(FileTOOCAN)
        if(ifileTOOCAN == 0) : 
            data     = load_TOOCAN(FileTOOCAN)
        else:
            data     = data + load_TOOCAN(FileTOOCAN)

    # selection of the MCS correctly identified and tracked 
    #for example : qc_MCS =  11108
    #First digit =1        MCS initiation OK 
    #Second digit = 1     MCS dissipation OK 
    #Third digit = 1          MCS not impacted by the image boundaries
    #Two last digit = 08     less than 8 interpolated GEO images occuring during the MCS tracking
    data = [data[iMCS] for iMCS in np.arange(0,len(data),1) if(int(data[iMCS].qc_MCS) <= 11108)]
    return data 


########################################
####################
def load_phase(wave, lat_min, lon_min, lat_max, lon_max, time_init, time_end, coeff, name_phase, iregion):
    indir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/PHASE/' + name_phase + '/'
    ds = xr.open_mfdataset(indir_data + '*' + wave + '*.nc', chunks = {'time' : 1}, parallel = True)
    if ((iregion == 0) | (iregion == 4)) :
        ds = ds.assign_coords(lon=(((ds.lon + 180) % 360) - 180)).sortby('lon')
    
    ds = ds.sel(time = slice(time_init, time_end), lat = slice(lat_min, lat_max), lon = slice(lon_min, lon_max))
    ds = xr.where(ds['magnitude'] < coeff, np.nan, ds)
    return ds



########################################
####################
def separate_data_MCS(data, iregion, region, latmin, latmax, lonmin, lonmax):
    ##########################################################################
    # Selection of the MCSs which occured within the region of interest
    ##########################################################################
    # And we separate in two groups with the classif, one < 5 hours and the anothers > 5 hours
    # data_MCSselection_G1 is for < 5h
    # data_MCSselection_G2 is for > 5h
    ######################################
    #### Select the data for Groupe 1  ###
    ######################################
    data_MCSselection_G1 = []
    for iMCS in np.arange(0,len(data),1):
            _id = [ilife for ilife in np.arange(0,data[iMCS].nbslot,1,dtype=int)
                  if( (data[iMCS].clusters.lon[ilife] > lonmin)  &  (data[iMCS].clusters.lon[ilife] < lonmax) 
                     &  (data[iMCS].clusters.lat[ilife] > latmin)  &  (data[iMCS].clusters.lat[ilife] < latmax)
                     &  (data[iMCS].classif == 1) ) ]
    #         for ilife in np.arange(0,data[iMCS].nbslot,1,dtype=int): 
    #             print(data[iMCS].clusters.lon[ilife])
            if(np.size(_id) > 0):
                data_MCSselection_G1.append(data[iMCS])
    print(len(data),len(data_MCSselection_G1))

    data_birth_G1 = []
    # Load the time of birth hours since 1970 into xarrayDataset
    for iMCS in np.arange(0,len(data_MCSselection_G1),1):
        birth = data_MCSselection_G1[iMCS].Utime_Init
        data_birth_G1.append(birth)
    attrs = {'units': 'days since 1970-01-01'}
    ds_T_I_G1 = xr.Dataset({'time': ('time', data_birth_G1, attrs)})
    ds_T_I_G1 = xr.decode_cf(ds_T_I_G1)

    data_birth_G1 = []
    # Load the time of end hours since 1970 into xarrayDataset
    for iMCS in np.arange(0,len(data_MCSselection_G1),1):
        birth = data_MCSselection_G1[iMCS].Utime_End
        data_birth_G1.append(birth)
    attrs = {'units': 'days since 1970-01-01'}
    ds_T_E_G1 = xr.Dataset({'time': ('time', data_birth_G1, attrs)})
    ds_T_E_G1 = xr.decode_cf(ds_T_E_G1)

    ######################################
    #### Select the data for Groupe 2  ###
    ######################################
    data_MCSselection_G2 = []
    for iMCS in np.arange(0,len(data),1):
            _id = [ilife for ilife in np.arange(0,data[iMCS].nbslot,1,dtype=int)
                  if( (data[iMCS].clusters.lon[ilife] > lonmin)  &  (data[iMCS].clusters.lon[ilife] < lonmax) 
                     &  (data[iMCS].clusters.lat[ilife] > latmin)  &  (data[iMCS].clusters.lat[ilife] < latmax)
                     &  (data[iMCS].classif > 1) ) ]
    #         for ilife in np.arange(0,data[iMCS].nbslot,1,dtype=int): 
    #             print(data[iMCS].clusters.lon[ilife])
            if(np.size(_id) > 0):
                data_MCSselection_G2.append(data[iMCS])
    print(len(data),len(data_MCSselection_G2))

    data_birth_G2 = []
    # Load the time of birth hours since 1970 into xarrayDataset
    for iMCS in np.arange(0,len(data_MCSselection_G2),1):
        birth = data_MCSselection_G2[iMCS].Utime_Init
        data_birth_G2.append(birth)
    attrs = {'units': 'days since 1970-01-01'}
    ds_T_I_G2 = xr.Dataset({'time': ('time', data_birth_G2, attrs)})
    ds_T_I_G2 = xr.decode_cf(ds_T_I_G2)

    data_birth_G2 = []
    for iMCS in np.arange(0,len(data_MCSselection_G1),1):
        birth = data_MCSselection_G1[iMCS].Utime_End
        data_birth_G2.append(birth)
    attrs = {'units': 'days since 1970-01-01'}
    ds_T_E_G2 = xr.Dataset({'time': ('time', data_birth_G2, attrs)})
    ds_T_E_G2 = xr.decode_cf(ds_T_E_G2)
    
    return data_MCSselection_G1, data_MCSselection_G2, ds_T_I_G1, ds_T_I_G2, ds_T_E_G1, ds_T_E_G2


        
