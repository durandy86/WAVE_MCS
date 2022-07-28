#!/usr/bin/env python
# coding: utf-8

import os, sys
os.environ["MALLOC_TRIM_THRESHOLD_"] = '0'
tempDir = os.environ['TMPDIR']

import numpy as np
import xarray as xr
import pandas as pd

from functionMaths import *

from dask.distributed import Client, LocalCluster
#
# Initialisation d'un cluster de 32 coeurs
cluster = LocalCluster(processes=False, n_workers=1, threads_per_worker=32, local_directory = tempDir, protocol = 'tcp://')
client = Client(cluster)
client
################################################
################################################
v = sys.argv[1]
nb_phase = [int(sys.argv[2])]
wave = [v]
lat_Ref = 0
outdir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/SCRIPTS/NETcdf/plotPHASE/2DMAP/phase_OLR/shear/'
coeff = 1

################################################
################################################
v = sys.argv[1]
nb_phase = [int(sys.argv[2])]
wave = [v]
lat_Ref = 0
coeff = 1
level = '850'
outdir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/SCRIPTS/NETcdf/plotPHASE/2DMAP/phase_OLR/z_' + level + '/'
# nb_phase = np.arange(0,32, dtype = 'int')

########################################################################################
indir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/RAW_ANOMALY/OLR/'
ds = xr.open_mfdataset(indir_data + 'anom_OLR_brut_ERA5_*.nc')
ds = ds.sel(time = slice('2001','2018')).astype(np.float32)

########################################################################################
indir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/RAW_ANOMALY/u/reduce/'
ds_u = xr.open_mfdataset(indir_data + '*z' + level + '*.nc', chunks = {'time' : 100}, parallel = True)
ds_u = ds_u.sel(time = slice('2001','2018')).astype(np.float32)

########################################################################################
indir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/RAW_ANOMALY/v/reduce/'
ds_v = xr.open_mfdataset(indir_data + '*z' + level + '*.nc', chunks = {'time' : 100}, parallel = True)
ds_v = ds_v.sel(time = slice('2001','2018')).astype(np.float32)
"""
########################################################################################
indir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/RAW_ANOMALY/d/reduce/'
ds_d = xr.open_mfdataset(indir_data + '*z' + level + '*.nc', chunks = {'time' : 100}, parallel = True)
ds_d = ds_d.sel(time = slice('2001','2018'))

########################################################################################
indir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/RAW_ANOMALY/vo/reduce/'
ds_vo = xr.open_mfdataset(indir_data + '*z' + level + '*.nc', chunks = {'time' : 100}, parallel = True)
ds_vo = ds_vo.sel(time = slice('2001','2018')).astype(np.float32)

########################################################################################
indir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/RAW_ANOMALY/t/reduce/'
ds_t = xr.open_mfdataset(indir_data + '*z' + level + '*.nc', chunks = {'time' : 100}, parallel = True)
ds_t = ds_t.sel(time = slice('2001','2018')).astype(np.float32)


ds = xr.merge([ds, ds_u, ds_v, ds_d, ds_vo, ds_t])
"""
### Shear
outdir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/SCRIPTS/NETcdf/plotPHASE/2DMAP/phase_OLR/shear/'
indir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/RAW_ANOMALY/u/reduce/'
_ds_u = xr.open_mfdataset(indir_data + '*z200*.nc', chunks = {'time' : 100}, parallel = True)
_ds_u = _ds_u.sel(time = slice('2001','2018'))
ds_u = _ds_u - ds_u 

indir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/RAW_ANOMALY/v/reduce/'
_ds_v = xr.open_mfdataset(indir_data + '*z200*.nc', chunks = {'time' : 100}, parallel = True)
_ds_v = _ds_v.sel(time = slice('2001','2018'))
ds_v = _ds_v - ds_v 

### TCWV
indir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/RAW_ANOMALY/TCWV/reduce/'
ds_t = xr.open_mfdataset(indir_data + '*.nc', chunks = {'time' : 100}, parallel = True)
ds_t = ds_t.sel(time = slice('2001','2018')).astype(np.float32)

### Merge all
ds = xr.merge([ds, ds_u, ds_v, ds_t])

########################################################################################
indir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/PHASE/OLR/'
for w in wave:
    if (w == 'OLR_Kelvin_SYM' or w == 'OLR_Rossby_SYM') :
        indir_data = indir_data + 'SYM/'
    if (w == 'OLR_Rossby' or w == 'OLR_Rossby_SYM'):
        lat_Ref = 10
    else :
        lat_Ref = 0

    for p in nb_phase:
        ds_phase = xr.open_mfdataset(indir_data + '*' + w + '*.nc', 
                                 chunks = {'time':100})

        ### Select only active wave
        ds_phase = xr.where(ds_phase['magnitude'] < 1, np.nan, ds_phase)
        ### Select the phase
        ds_phase = xr.where(ds_phase['phase'] != p, np.nan, ds_phase)

        
        _ds_phase = ds_phase['magnitude']
        ### Stock datasets ds into tempory to avoid conflict        
        _ds = ds
        _ds = xr.where(_ds_phase.isnull() != True,
                                 _ds, np.nan)

        _ds = _ds.mean(dim = ['time']).compute()
        _ds.to_netcdf(outdir_data + w + '_without_ref_phase_nb' + str(p) + '.nc')
        """
        _ds_phase = ds_phase['magnitude'].sel(lat = lat_Ref, method = 'nearest')
        del _ds_phase['lat']
        _ds_phase = _ds_phase.expand_dims(lat = ds['lat'], axis = 1)
        
        _ds = ds
        _ds = xr.where(_ds_phase.isnull() != True,
                                 _ds, np.nan)
        _ds = _ds.mean(dim = ['time']).compute()
        _ds.to_netcdf(outdir_data + w + '_with_ref_' + str(lat_Ref) + '_phase_nb' + str(p) + '.nc')
        """



