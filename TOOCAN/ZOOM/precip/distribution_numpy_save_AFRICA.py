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

folder_out = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/ANALYSIS_TOOCAN/Distribution/'


time_init = '2012'
time_end = '2017'
coeff = 1
iregion = 0
phase_sel = ['divergence/200'] #'OLR', 'TCWV', 'divergence/850', 
phase_p = ['divergence_200'] #'OLR', 'TCWV', 'divergence_850', 
wave_SEL = ['TD', 'MRG', 'MJO'] #'Kelvin', 'Rossby', 

##################################################
region, latmin, latmax, lonmin, lonmax = parameter_TOOCAN(iregion)
data = load_Data_TOOCAN(iregion)
################################
#--- Create different box ---###
################################
geos = int(sys.argv[1])
k = 0 
if geos == 0 :
    zone = 'Congo'
    latmin, latmax =  -8., 8.
    lonmin, lonmax = 15., 35.
elif geos == 1 :
    zone = 'Atlantique'
    latmin, latmax =  0., 10.
    lonmin, lonmax = -40., -15.
elif geos == 2 :
    zone = 'Afrique_O'
    latmin, latmax =  5., 15.
    lonmin, lonmax = -15., 10.
"""
name_text = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/ANALYSIS_TOOCAN/Distribution/life_data_MCS_' + zone + '.txt'
file_text = open(name_text, 'w')
file_text.write('START \n')
file_text.close
"""
data_MCSselection_G1, data_MCSselection_G2, ds_T_I_G1, ds_T_I_G2, ds_T_E_G1, ds_T_E_G2 = separate_data_MCS(data, iregion, region, latmin, latmax, lonmin, lonmax)

for name_phase in phase_sel :
    #####################################
    ###--- Wave and date to select ---###
    #####################################
    # Be carefull that the date match TOOCAN
    for wave in wave_SEL :	
        #########################################################################
        ### data phase OLR ###
        #########################################################################
        ds = load_phase(wave, latmin, lonmin, latmax, lonmax, time_init, time_end, coeff, name_phase, iregion)
        ds = ds.phase.load()

        nb_occu = np.zeros(32, dtype = np.int32)
        for p in range(32):
            __ = xr.where(ds == p, 1, np.nan)
            nb_occu[p] = __.sum()
        nb_occu = nb_occu/np.sum(nb_occu)
        """
        #### Load Parameter for MCS and TOOCAN
        ### p_distri_G1, p_distri_G2, nb_MCS_1, nb_MCS_2, surface_G1, TbMin_G1, surface_G2, TbMin_G2
        p_distri_G1, p_distri_G2, nb_MCS_1, nb_MCS_2, surface_G1, TbMin_G1, surface_G2, TbMin_G2 = \
                            calcul_Proba_wave(ds, ds_T_I_G1, ds_T_E_G1, data_MCSselection_G1, ds_T_I_G2, ds_T_E_G2, data_MCSselection_G2, 'init', name_phase)

        ### save numpy array
        save_file = np.stack((nb_occu, p_distri_G1, p_distri_G2, surface_G1, TbMin_G1, surface_G2, TbMin_G2))
        outfile = folder_out + zone + '_nbMCS_' + wave + '_' + phase_p[k] + '_init.npy'
        np.save(outfile, save_file)

        file_text = open(name_text, 'a')
        str_file = 'Pour ' + name_phase + ' ' + zone + ' et ' + wave + ' Nb MCS <5h = ' + str(nb_MCS_1) + ' Nb MCS > 5H = ' + str(nb_MCS_2) + '\n'
        file_text.write(str_file)
        file_text.close() 

        """
        #######################################
        ### Compute nb of MCS in an active wave
        ### p_distri_G1, p_distri_G2, nb_MCS_1, nb_MCS_2, nb_total_1, 0, nb_total_G2, 0
        p_distri_G1, p_distri_G2, nb_MCS_G1, nb_MCS_G2, surface_G1, TbMin_G1, surface_G2, TbMin_G2 = \
                            calcul_Proba_wave(ds, ds_T_I_G1, ds_T_E_G1, data_MCSselection_G1, ds_T_I_G2, ds_T_E_G2, data_MCSselection_G2, 'life', name_phase, iregion)

        ### save numpy array
        save_file = np.stack((nb_occu, p_distri_G1, p_distri_G2, surface_G1, TbMin_G1, surface_G2, TbMin_G2))
        outfile = folder_out + zone + '_nbMCS_' + wave + '_' + phase_p[k] + '_life.npy'
        np.save(outfile, save_file)
    k = k + 1
