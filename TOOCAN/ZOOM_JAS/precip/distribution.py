#!/usr/bin/env python
# coding: utf-8

import os, sys
# os.environ["MALLOC_TRIM_THRESHOLD_"] = '0'
import numpy as np
import xarray as xr
import pandas as pd

from matplotlib import pyplot as plt

from Load_TOOCAN import *
from Load_DATA import *
from function_Maths import *
from calcul import *

plt.rc("figure", figsize=(12,10))
plt.rc("font", size=14)
tempDir = os.environ['TMPDIR']

from dask.distributed import Client, LocalCluster
# Initialisation d'un cluster de 32 coeurs
cluster = LocalCluster(processes=False, n_workers=1, threads_per_worker=32, protocol = 'tcp', local_directory = tempDir)#, host = '137.129.155.67')
client = Client(cluster)
client

##################################################
##################################################
##################################################

outfig = '/cnrm/tropics/commun/DATACOMMUN/WAVE/RAPPORT/FIGURES/TOOCAN/Distribution/'
positions = (np.arange(0,32*2,8)) # Pour Kelvin
labels = (np.arange(0,32*2,8))%32

time_init = '2012'
time_end = '2018'
coeff = 1
iregion = int(sys.argv[1])
bins = np.arange(0,32*2, dtype = np.int)

##################################################
region, latmin, latmax, lonmin, lonmax = parameter_TOOCAN(iregion)
data = load_Data_TOOCAN(iregion)
################################
#--- Create different box ---###
################################
zone_geo = ['Congo','Atlantique','Afrique_O']
for zone in zone_geo :
    if iregion == 0 :
        if zone == 'Congo':
            print('COUCOU \n')
            latmin, latmax =  -10., 10.
            lonmin, lonmax = 10., 35.
        elif zone == 'Atlantique':
            latmin, latmax =  0., 20.
            lonmin, lonmax = -40., -20.
        elif zone == 'Afrique_O':
            latmin, latmax =  0., 20.
            lonmin, lonmax = -20., 10.

    data_MCSselection_G1, data_MCSselection_G2, ds_T_I_G1, ds_T_I_G2, ds_T_E_G1, ds_T_E_G2 = separate_data_MCS(data, iregion, region, latmin, latmax, lonmin, lonmax)

    #####################################
    ###--- Wave and date to select ---###
    #####################################
    # Be carefull that the date match TOOCAN
    wave_SEL = ['Kelvin', 'Rossby', 'TD', 'MRG', 'MJO', 'WIG','EIG']
    for wave in wave_SEL :	
        #########################################################################
        ### data phase OLR ###
        #########################################################################
        ds = load_phase(wave, latmin, lonmin, latmax, lonmax, time_init, time_end, coeff, 'OLR', iregion)
        ds = ds.phase.load()
        print(ds.isel(time = 1).values)

        nb_occu = np.zeros(32, dtype = np.int32)
        for p in range(32):
            __ = xr.where(ds == p, 1, np.nan)
            nb_occu[p] = __.sum()
        nb_occu = nb_occu/np.sum(nb_occu)
        nb_occu[-1] = nb_occu[0]

        #### Load Parameter for MCS and TOOCAN
        p_distri_G1, p_distri_G2, nb_MCS_1, nb_MCS_2, surface_G1, TbMin_G1, surface_G2, TbMin_G2 = \
                            calcul_Proba_wave(ds, ds_T_I_G1, ds_T_E_G1, data_MCSselection_G1, ds_T_I_G2, ds_T_E_G2, data_MCSselection_G2, 'init', 'OLR')


        ### Compute nb of MCS in an active wave
        nb_occu_with_MCS = np.zeros(32, dtype = np.int32)
        for iMCS in np.arange(0,len(data_MCSselection_G2),1):
            for ilife in np.arange(0,data_MCSselection_G2[iMCS].nbslot,1,dtype=int):
                _time = ds_T_I_G2.isel(time = iMCS).time.values + np.timedelta64((30*ilife), 'm')
                nb_phase = ds.sel(time = _time, 
                         lat = data_MCSselection_G2[iMCS].clusters.lat[ilife],
                         lon = data_MCSselection_G2[iMCS].clusters.lon[ilife],
                         method = 'nearest').values

                if np.isnan(nb_phase) == False :
                    nb_phase = int(nb_phase)
                    nb_occu_with_MCS[nb_phase] = nb_occu_with_MCS[nb_phase] + 1

        nb_occu_with_MCS = nb_occu_with_MCS/np.sum(nb_occu_with_MCS)
        nb_occu_with_MCS[-1] = nb_occu_with_MCS[0]
        varFiltre = 'OLR'
        if ((varFiltre != 'OLR') | (varFiltre != 'divergence/200')):
            nb_occu = np.roll(nb_occu, 16, axis=0) 
            nb_occu_with_MCS = np.roll(nb_occu_with_MCS, 16, axis=0) 

        #########################################################################
        ############################### ploting #################################
        #########################################################################
        __ = np.arange(1,33,4)
        theta = np.linspace(0, 2*np.pi, 32)
        theta[-1] = theta[0]
        r = p_distri_G2
        r[-1] = r[0]
        nb_occu_with_MCS[-1] = nb_occu_with_MCS[0]

        plt.figure(figsize=(12,10))
        ax = plt.subplot(111, projection='polar')
        ax.plot(theta, r, label='Pdf MCS initiation')
        ax.plot(theta, nb_occu, label='Pdf active wave')
        ax.set_theta_direction(-1)  # theta increasing clockwise
        ax.set_theta_offset(pi)
        ax.set_xticklabels(__)
        ax.set_rmax(0.06)
        plt.legend()
        plt.savefig(outfig + zone + '_initiation_' + wave + '_OLR_sup5H.png')
        plt.close()

        plt.figure(figsize=(12,10))
        ax = plt.subplot(111, projection='polar')
        ax.plot(theta, nb_occu_with_MCS, label='Pdf with MCS')
        ax.plot(theta, nb_occu, label='Pdf active wave')
        ax.set_theta_direction(-1)  # theta increasing clockwise
        ax.set_theta_offset(pi)
        ax.set_xticklabels(__)
        ax.set_rmax(0.06)
        plt.legend()
        plt.savefig(outfig + zone + '_nbMCS_' + wave + '_OLR_sup5H.png')
        plt.close()


        #########################################################################
        ### data phase TCWV ###
        #########################################################################
        ds = load_phase(wave, latmin, lonmin, latmax, lonmax, time_init, time_end, coeff, 'TCWV', iregion)
        ds['phase'] = ds.phase
        ds = ds['phase'].load()

        nb_occu = np.zeros(32, dtype = np.int32)
        for p in range(32):
            __ = xr.where(ds == p, 1, np.nan)
            nb_occu[p] = __.sum()
        nb_occu = nb_occu/np.sum(nb_occu)
        nb_occu[-1] = nb_occu[0]


        #### Load Parameter for MCS and TOOCAN
        p_distri_G1, p_distri_G2, nb_MCS_1, nb_MCS_2, surface_G1, TbMin_G1, surface_G2, TbMin_G2 = \
                            calcul_Proba_wave(ds, ds_T_I_G1, ds_T_E_G1, data_MCSselection_G1, ds_T_I_G2, ds_T_E_G2, data_MCSselection_G2, 'init', 'TCWV')

        ### Compute nb of MCS in an active wave
        for iMCS in np.arange(0,len(data_MCSselection_G2),1):
            for ilife in np.arange(0,data_MCSselection_G2[iMCS].nbslot,1,dtype=int):
                _time = ds_T_I_G2.isel(time = iMCS).time.values + np.timedelta64((30*ilife), 'm')
                nb_phase = ds.sel(time = _time, 
                         lat = data_MCSselection_G2[iMCS].clusters.lat[ilife],
                         lon = data_MCSselection_G2[iMCS].clusters.lon[ilife],
                         method = 'nearest').values

                if np.isnan(nb_phase) == False :
                    nb_phase = int(nb_phase)
                    nb_occu_with_MCS[nb_phase] = nb_occu_with_MCS[nb_phase] + 1

        nb_occu_with_MCS = nb_occu_with_MCS/np.sum(nb_occu_with_MCS)
        nb_occu_with_MCS[-1] = nb_occu_with_MCS[0]

        #########################################################################
        ############################### ploting #################################
        #########################################################################
        __ = np.arange(1,33,4)
        theta = np.linspace(0, 2*np.pi, 32)
        theta[-1] = theta[0]
        r = p_distri_G2
        r[-1] = r[0]
        nb_occu_with_MCS[-1] = nb_occu_with_MCS[0]

        plt.figure(figsize=(12,10))
        ax = plt.subplot(111, projection='polar')
        ax.plot(theta, r, label='Pdf MCS initiation')
        ax.plot(theta, nb_occu, label='Pdf active wave')
        ax.set_theta_direction(-1)  # theta increasing clockwise
        ax.set_theta_offset(pi)
        ax.set_xticklabels(__)
        ax.set_rmax(0.06)
        plt.legend()
        plt.savefig(outfig + zone + '_initiation_' + wave + '_TCWV_sup5H.png')
        plt.close()

        plt.figure(figsize=(12,10))
        ax = plt.subplot(111, projection='polar')
        ax.plot(theta, nb_occu_with_MCS, label='Pdf with MCS')
        ax.plot(theta, nb_occu, label='Pdf active wave')
        ax.set_theta_direction(-1)  # theta increasing clockwise
        ax.set_theta_offset(pi)
        ax.set_xticklabels(__)
        ax.set_rmax(0.06)
        plt.legend()
        plt.savefig(outfig + zone + '_nbMCS_' + wave + '_TCWV_sup5H.png')
        plt.close()

        #########################################################################
        ### data phase divergence 850 hPa ###
        #########################################################################
        ds = load_phase(wave, latmin, lonmin, latmax, lonmax, time_init, time_end, coeff, 'divergence/850', iregion)
        ds['phase'] = ds.phase
        ds = ds['phase'].load()

        nb_occu = np.zeros(32, dtype = np.int32)
        for p in range(32):
            __ = xr.where(ds == p, 1, np.nan)
            nb_occu[p] = __.sum()
        nb_occu = nb_occu/np.sum(nb_occu)
        nb_occu[-1] = nb_occu[0]

        #### Load Parameter for MCS and TOOCAN
        p_distri_G1, p_distri_G2, nb_MCS_1, nb_MCS_2, surface_G1, TbMin_G1, surface_G2, TbMin_G2 = \
                            calcul_Proba_wave(ds, ds_T_I_G1, ds_T_E_G1, data_MCSselection_G1, ds_T_I_G2, ds_T_E_G2, data_MCSselection_G2, 'init', 'divergence/850')


        ### Compute nb of MCS in an active wave
        for iMCS in np.arange(0,len(data_MCSselection_G2),1):
            for ilife in np.arange(0,data_MCSselection_G2[iMCS].nbslot,1,dtype=int):
                _time = ds_T_I_G2.isel(time = iMCS).time.values + np.timedelta64((30*ilife), 'm')
                nb_phase = ds.sel(time = _time, 
                         lat = data_MCSselection_G2[iMCS].clusters.lat[ilife],
                         lon = data_MCSselection_G2[iMCS].clusters.lon[ilife],
                         method = 'nearest').values

                if np.isnan(nb_phase) == False :
                    nb_phase = int(nb_phase)
                    nb_occu_with_MCS[nb_phase] = nb_occu_with_MCS[nb_phase] + 1

        nb_occu_with_MCS = nb_occu_with_MCS/np.sum(nb_occu_with_MCS)
        nb_occu_with_MCS[-1] = nb_occu_with_MCS[0]

        #########################################################################
        ############################### ploting #################################
        #########################################################################
        __ = np.arange(1,33,4)
        theta = np.linspace(0, 2*np.pi, 32)
        theta[-1] = theta[0]
        r = p_distri_G2
        r[-1] = r[0]
        nb_occu_with_MCS[-1] = nb_occu_with_MCS[0]

        plt.figure(figsize=(12,10))
        ax = plt.subplot(111, projection='polar')
        ax.plot(theta, r, label='Pdf MCS initiation')
        ax.plot(theta, nb_occu, label='Pdf active wave')
        ax.set_theta_direction(-1)  # theta increasing clockwise
        ax.set_theta_offset(pi)
        ax.set_xticklabels(__)
        ax.set_rmax(0.06)
        plt.legend()
        plt.savefig(outfig + zone + '_initiation_' + wave + '_div850_sup5H.png')
        plt.close()

        plt.figure(figsize=(12,10))
        ax = plt.subplot(111, projection='polar')
        ax.plot(theta, nb_occu_with_MCS, label='Pdf with MCS')
        ax.plot(theta, nb_occu, label='Pdf active wave')
        ax.set_theta_direction(-1)  # theta increasing clockwise
        ax.set_theta_offset(pi)
        ax.set_xticklabels(__)
        ax.set_rmax(0.06)
        plt.legend()
        plt.savefig(outfig + zone + '_nbMCS_' + wave + '_div850_sup5H.png')
        plt.close()

        #########################################################################
        ### data phase divergence 200 hPa ###
        #########################################################################
        ds = load_phase(wave, latmin, lonmin, latmax, lonmax, time_init, time_end, coeff, 'divergence/200', iregion)
        ds['phase'] = ds.phase
        ds = ds['phase'].load()

        nb_occu = np.zeros(32, dtype = np.int32)
        for p in range(32):
            __ = xr.where(ds == p, 1, np.nan)
            nb_occu[p] = __.sum()
        nb_occu = nb_occu/np.sum(nb_occu)
        nb_occu[-1] = nb_occu[0]

        #### Load Parameter for MCS and TOOCAN
        p_distri_G1, p_distri_G2, nb_MCS_1, nb_MCS_2, surface_G1, TbMin_G1, surface_G2, TbMin_G2 = \
                            calcul_Proba_wave(ds, ds_T_I_G1, ds_T_E_G1, data_MCSselection_G1, ds_T_I_G2, ds_T_E_G2, data_MCSselection_G2, 'init', 'divergence/200')


        ### Compute nb of MCS in an active wave
        for iMCS in np.arange(0,len(data_MCSselection_G2),1):
            for ilife in np.arange(0,data_MCSselection_G2[iMCS].nbslot,1,dtype=int):
                _time = ds_T_I_G2.isel(time = iMCS).time.values + np.timedelta64((30*ilife), 'm')
                nb_phase = ds.sel(time = _time, 
                         lat = data_MCSselection_G2[iMCS].clusters.lat[ilife],
                         lon = data_MCSselection_G2[iMCS].clusters.lon[ilife],
                         method = 'nearest').values

                if np.isnan(nb_phase) == False :
                    nb_phase = int(nb_phase)
                    nb_occu_with_MCS[nb_phase] = nb_occu_with_MCS[nb_phase] + 1

        nb_occu_with_MCS = nb_occu_with_MCS/np.sum(nb_occu_with_MCS)
        nb_occu_with_MCS[-1] = nb_occu_with_MCS[0]

        varFiltre = 'divergence/200'
        if ((varFiltre != 'OLR') | (varFiltre != 'divergence/200')):
            nb_occu = np.roll(nb_occu, 16, axis=0) 
            nb_occu_with_MCS = np.roll(nb_occu_with_MCS, 16, axis=0) 

        #########################################################################
        ############################### ploting #################################
        #########################################################################
        __ = np.arange(1,33,4)
        theta = np.linspace(0, 2*np.pi, 32)
        theta[-1] = theta[0]
        r = p_distri_G2
        r[-1] = r[0]
        nb_occu_with_MCS[-1] = nb_occu_with_MCS[0]

        plt.figure(figsize=(12,10))
        ax = plt.subplot(111, projection='polar')
        ax.plot(theta, r, label='Pdf MCS initiation')
        ax.plot(theta, nb_occu, label='Pdf active wave')
        ax.set_theta_direction(-1)  # theta increasing clockwise
        ax.set_theta_offset(pi)
        ax.set_xticklabels(__)
        ax.set_rmax(0.06)
        plt.legend()
        plt.savefig(outfig + zone + '_initiation_' + wave + '_div200_sup5H.png')
        plt.close()

        plt.figure(figsize=(12,10))
        ax = plt.subplot(111, projection='polar')
        ax.plot(theta, nb_occu_with_MCS, label='Pdf with MCS')
        ax.plot(theta, nb_occu, label='Pdf active wave')
        ax.set_theta_direction(-1)  # theta increasing clockwise
        ax.set_theta_offset(pi)
        ax.set_xticklabels(__)
        ax.set_rmax(0.06)
        plt.legend()
        plt.savefig(outfig + zone + '_nbMCS_' + wave + '_div200_sup5H.png')
        plt.close()
  
