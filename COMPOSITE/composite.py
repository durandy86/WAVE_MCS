#!/usr/bin/env python
# coding: utf-8

import os, sys
os.environ["MALLOC_TRIM_THRESHOLD_"] = '0'

# tempDir = os.environ['TMPDIR']
import numpy as np
import xarray as xr
import pandas as pd

from dask.distributed import Client, LocalCluster
#
# Initialisation d'un cluster de 32 coeurs
cluster = LocalCluster(processes=False, n_workers=1, threads_per_worker=32, protocol = 'tcp://')#, dashboard_address= None, local_directory = tempDir)
client = Client(cluster)
client
#####################################################
#####################################################
#####################################################

name_wave = ['Kelvin','Rossby','MRG','TD','MJO']
level = [850, 200] 
outdir = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/COMPOSITE/'
indir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/PHASE/OLR/'
year = np.arange(2001,2019)
coeff = 1 # Sueil Sigma

########################################################################
########################################################################
########################################################################
for l in level :
    for w in name_wave :
        for y in year :
            indir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/PHASE/OLR/'
            ds_phase = xr.open_mfdataset(indir_data + '*' + w + '*' + str(y) + '.nc', chunks = {'time' : 1 })
            ds_phase = ds_phase.sel( lat = slice(-25,25))

            ### When level is l == 850 we make the composite for the 2D variable 
            if l == 850 :
                #######################################################################################
                indir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/RAW_ANOMALY/OLR/'
                ds_OLR = xr.open_mfdataset(indir_data + '*.nc', chunks = {'time' : 1}, parallel = True)
                ds_OLR = ds_OLR.sel(time = ds_phase.time, lat = ds_phase.lat).astype(np.float32)

                #######################################################################################
                indir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/RAW_ANOMALY/TCWV/reduce/'
                ds_tw = xr.open_mfdataset(indir_data + '*.nc', chunks = {'time' : 1}, parallel = True)
                ds_tw = ds_tw.sel(time = ds_phase.time, lat = ds_phase.lat).astype(np.float32)
            
            #######################################################################################
            indir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/RAW_ANOMALY/u/reduce/'
            ds_u = xr.open_mfdataset(indir_data + '*z' + str(l) + '*.nc', chunks = {'time' : 1}, parallel = True)
            ds_u = ds_u.sel(time = ds_phase.time, lat = ds_phase.lat).astype(np.float32)
          
            #######################################################################################
            indir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/RAW_ANOMALY/v/reduce/'
            ds_v = xr.open_mfdataset(indir_data + '*z' + str(l) + '*.nc', chunks = {'time' : 1}, parallel = True)
            ds_v = ds_v.sel(time = ds_phase.time, lat = ds_phase.lat).astype(np.float32)

            #######################################################################################
            indir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/RAW_ANOMALY/d/reduce/'
            ds_d = xr.open_mfdataset(indir_data + '*z' + str(l) + '*.nc', chunks = {'time' : 1}, parallel = True)
            ds_d = ds_d.sel(time = ds_phase.time, lat = ds_phase.lat).astype(np.float32)

            #######################################################################################
            indir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/RAW_ANOMALY/vo/reduce/'
            ds_vo = xr.open_mfdataset(indir_data + '*z' + str(l) + '*.nc', chunks = {'time' : 1}, parallel = True)
            ds_vo = ds_vo.sel(time = ds_phase.time, lat = ds_phase.lat).astype(np.float32)
            
            #######################################################################################
            indir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/RAW_ANOMALY/t/reduce/'
            ds_t = xr.open_mfdataset(indir_data + '*z' + str(l) + '*.nc', chunks = {'time' : 1}, parallel = True)
            ds_t = ds_t.sel(time = ds_phase.time, lat = ds_phase.lat).astype(np.float32)
            ds_t = ds_t.isel(level = 0) 

            if l == 850 :
                ds = xr.merge([ds_OLR, ds_u, ds_v, ds_d, ds_t, ds_vo, ds_tw])
            else :
                ds = xr.merge([ds_u, ds_v, ds_d, ds_t, ds_vo])


            ds['magnitude'] = ds_phase['magnitude'].astype(np.float32)

            ### Crate a new dimension for the phase
            ds = ds.expand_dims(phase = np.arange(0,32,1, dtype = 'int'), axis = 3).copy() 
            ds = ds.load()
            _ds_phase = ds_phase.phase.astype(np.int)
            _ds_phase = _ds_phase.compute()
            ### When level is l == 850 we make the composite for the 2D variable 
            if l == 850 :
                var = ['OLR_ano', 'tcwv_ano', 'u_ano', 'v_ano', 'd_ano', 't_ano', 'vo_ano', 'magnitude']
            else :
                var = ['u_ano', 'v_ano', 'd_ano', 't_ano', 'vo_ano', 'magnitude']

            for v in var :
                for i in range(32):
                    ds[v][...,i] = xr.where(_ds_phase == i, ds[v][...,i], np.nan)
 
            ds = xr.where(ds.magnitude < 1, np.nan, ds) ### 
            ds = ds.mean(dim = ['time'])
            ds = ds.compute()
            ds.to_netcdf(outdir + 'mean_anomalie_z' + str(l) + '_phase_OLR_'+ w + '_' + str(y) + '_JAS.nc')

