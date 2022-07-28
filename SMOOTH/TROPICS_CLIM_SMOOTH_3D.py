#!/usr/bin/env python
# coding: utf-8

import os, sys
tempDir = os.environ['TMPDIR']
os.environ["MALLOC_TRIM_THRESHOLD_"] = '0'
tempdir = tempDir + '/level_CLIM/'

import numpy as np
import h5netcdf

import xarray as xr
import xrft
import pandas as pd

from dask.distributed import Client, LocalCluster
cluster = LocalCluster(processes=False, n_workers=1, threads_per_worker=32, protocol = 'tcp', local_directory = tempDir)
client = Client(cluster)
#client


############################################################################
indir = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/RAW_CLIM/'
outdir = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/SMOTHED_CLIM/'
minYear = '1990'
maxYear = '2020'
zLevel = [1000, 950, 800, 700, 600, 400, 300] #, 900, 700, 500, 200]

#number of harmonics to keep
nbSampl = 8
nbHarm = 3  # 0 for the mean, 1 for the the annual cycle, etc...
nbHarmKeep = nbHarm
#Variable to smooth
var = ['t']

############################################################################
nbHarmKeep = nbHarm + 1 #+1 car on garde la moyenne
for v in var :
    # 
    for z in zLevel :
        # i_ds = ds.sel(level = z)
        # _ds = _ds.expand_dims(level = 1)
        # _ds['level'] = z
        exp = 'clim_' + v + '_' + 'z' + str(z) + '_brut_ERA5_3H_1990_2020'
        ds = xr.open_dataset(indir + exp + '.nc', chunks = {'latitude' : 1})
        prefix = 'SMOTHEDCLIM_' + v + '_'
        # ## Selection
        # Dans la cellule au-dessus on a sortie dans un dataset la climatologie brute de ERA5, et on souhaite "smooth" le signal. On ne garde que les premières harmonique du signal. On utilise le package pour faire une transformée de Fourier sur des tableaux Xarray. https://xrft.readthedocs.io/en/latest/index.html

        tcwvhat = xrft.fft(ds[v], dim="time", true_phase=False, true_amplitude=True)
        tcwvhat = xr.where(tcwvhat.freq_time <= -nbHarmKeep/(86400*366),  0., tcwvhat)
        tcwvhat = xr.where(tcwvhat.freq_time >= nbHarmKeep/(86400*366),  0., tcwvhat)
        tcwv_Sm = xrft.ifft(tcwvhat, dim = 'freq_time', true_phase=False, true_amplitude=True) # Signal in direct space
        tcwv_Sm['time'] = ds.time


        __ds = tcwv_Sm.real
        __ds = __ds.astype(np.float32)
        __ds = __ds.to_dataset(name = v+'_smooth')
        __ds = __ds.compute()
        """
        datasets = list(split_by_chunks(ds))
        paths_l = [create_filepath_l(ds, 
                          prefix = prefix,
                          root_path = tempdir) for ds in datasets]
        """
        # print(tempdir + prefix + str(z) + '.nc')
        __ds.to_netcdf(outdir + 'clim_' + v + '_' + 'z' + str(z) + '_smooth_ERA5_3H_' + str(minYear) + '_' + str(maxYear) + '.nc', unlimited_dims={'time':True})

print('script end')
