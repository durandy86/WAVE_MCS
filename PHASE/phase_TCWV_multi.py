#!/usr/bin/env python
# coding: utf-8

import numpy as np
import xarray as xr
import xarray.ufuncs as xu
import xrft
import pandas as pd

from functionMaths import *


import os, sys
tempDir = os.environ['TMPDIR']

from dask.distributed import Client, LocalCluster
# Initialisation d'un cluster de 32 coeurs
cluster = LocalCluster(processes=False, n_workers=1, threads_per_worker=32, local_directory = tempDir, protocol = 'tcp://')
client = Client(cluster)
client

############################################
############################################
############################################

wave = ['TCWV_Kelvin','TCWV_Rossby','TCWV_TD','TCWV_MJO','TCWV_MRG','TCWV_WIG','TCWV_EIG']
outdir = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/ANALYSIS/PHASE/TCWV/'
indir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/FILTERED_ANOMALY/TCWV/reduce/'
coeff = 1

ds = xr.open_mfdataset(indir_data + '*.nc', parallel = True)
ds = ds.sel(lat = slice(-30,30))

### Compute temporal derivative
ds_DDt = ds.differentiate('time', datetime_unit = "s")

ds_std = ds.std(dim = ['time'])
ds_DDt_std = ds_DDt.std(dim = ['time'])

# Normalizing by standard deviation
ds_norm = ds/ds_std
ds_DDt = ds_DDt/ds_DDt_std

### Stock in one datasets
space_phase = np.linspace(np.pi, -np.pi, 33, endpoint = True)
d_space_phase = (space_phase[1] - space_phase[2])/2
space_phase = space_phase[:-1] - d_space_phase
nb_phase = np.linspace(0, 31, 32, dtype = 'int')


for w in wave:  
    ds_norm[w + '_DT'] = ds_DDt[w]

year = [sys.argv[1]]
for y in year:
    for w in wave:
        v1 = ds_norm[w].sel(time = str(y))
        v2 = ds_norm[w+'_DT'].sel(time = str(y))
        ampl, magn = amplitudeMagnitude(v1, v2)
        phase = magn.copy()
        phase = (phase * 0).astype(np.int8)
        phase = phase.rename('phase')
        for i in range(nb_phase.size-1) :
            print(i)
            phase = xr.where( ((ampl <=  space_phase[i]) & ( ampl > space_phase[i+1]) == True),
                                               i+1, phase)


        phase = xr.where( ((ampl <=  space_phase[-1]) | ( ampl > space_phase[0]) == True),
                                               0, phase)
        # __ds['magnitude'] = xr.where(_ds['amplitude'].isnull() != True, _ds['magnitude'], np.nan)
        __ds = magn.to_dataset(name = 'magnitude')
        __ds['phase'] = phase
        ds_netcdf = __ds.compute()
        ds_netcdf.to_netcdf(outdir + 'phase_' + w + '_' + str(y) + '.nc')
