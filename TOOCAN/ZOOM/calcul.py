#!/usr/bin/env python
# coding: utf-8

import os, sys

import numpy as np
import xarray as xr
import pandas as pd

from function_Maths import *

def calcul_Proba_wave(ds, ds_T_I_G1, ds_T_E_G1, data_MCSselection_G1, ds_T_I_G2, ds_T_E_G2, data_MCSselection_G2, condition, varFiltre):

    if condition == 'init':
        # # Initialisation des MCS
        distri_G1, nb_MCS_G1, surface_G1, TbMin_G1 = compute_distribution_init(ds, ds_T_I_G1, data_MCSselection_G1)
        ## Pour les MCS supèrieur à 5h ##
        distri_G2, nb_MCS_G2, surface_G2, TbMin_G2 = compute_distribution_init(ds, ds_T_I_G2, data_MCSselection_G2)
        p_distri_G1, p_distri_G2 = distri_G1/(nb_MCS_G1 + nb_MCS_G2), distri_G2/(nb_MCS_G1 + nb_MCS_G2)

    elif condition == 'end' :
        # # Initialisation des MCS
        ## Pour les MCS inférieur à 5h ##
        distri_G1, nb_MCS_G1 = compute_distribution_end(ds, ds_T_E_G1, data_MCSselection_G1)
        ## Pour les MCS supèrieur à 5h ##
        distri_G2, nb_MCS_G2 = compute_distribution_end(ds, ds_T_I_G2, data_MCSselection_G2)

        p_distri_G1, p_distri_G2 = distri_G1/(nb_MCS_G1 + nb_MCS_G2), distri_G2/(nb_MCS_G1 + nb_MCS_G2)

        surface_G1, TbMin_G1 = [0],[0]
        surface_G2, TbMin_G2 = [0],[0]

    elif condition == 'life' : 
        ### surface == nb de MCS total par pas de temps
        # # Initialisation des MCS
        ## Pour les MCS inférieur à 5h ##
        distri_G1, nb_MCS_G1, surface_G1 = compute_distribution_life(ds, ds_T_I_G1, data_MCSselection_G1)
        ## Pour les MCS supèrieur à 5h ##
        distri_G2, nb_MCS_G2, surface_G2 = compute_distribution_life(ds, ds_T_I_G2, data_MCSselection_G2)

        p_distri_G1, p_distri_G2 = distri_G1/(nb_MCS_G1 + nb_MCS_G2), distri_G2/(nb_MCS_G1 + nb_MCS_G2)

        TbMin_G1 = [0]
        TbMin_G2 = [0]
    """
    if condition == 'init':
        if ((varFiltre != 'OLR') | (varFiltre != 'divergence/200')):
            p_distri_G1 = np.roll(p_distri_G1, 16, axis=0) 
            p_distri_G2 = np.roll(p_distri_G2, 16, axis=0) 
            surface_G1 = np.roll(surface_G1, 16, axis=0) 
            TbMin_G1 = np.roll(TbMin_G1, 16, axis=0) 
            surface_G2 = np.roll(surface_G2, 16, axis=0) 
            TbMin_G2 = np.roll(TbMin_G2, 16, axis=0) 
    """
 
    return p_distri_G1, p_distri_G2, nb_MCS_G1, nb_MCS_G2, surface_G1, TbMin_G1, surface_G2, TbMin_G2
