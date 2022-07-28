#!/usr/bin/env python
# coding: utf-8

import os, sys
import numpy as np
import itertools

import xarray as xr
import xrft
import pandas as pd

from dask.diagnostics import ProgressBar
from dask.distributed import Client, LocalCluster
#
# Initialisation d'un cluster de 32 coeurs
tempDir = os.environ['TMPDIR']
cluster = LocalCluster(processes=False, n_workers=1, threads_per_worker=32, local_directory = tempDir)
client = Client(cluster)
print(client)

############################################################################
indir = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/RAW_CLIM/'
outdir = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/SMOTHED_CLIM/'
minYear = '1990'
maxYear = '2020'
#number of harmonics to keep
nbSampl = 4
nbHarm = 3  # 0 for the mean, 1 for the the annual cycle, etc...
nbHarmKeep = nbSampl*nbHarm
#Variable to smooth
var = ['u','v','vo']
############################################################################
def split_by_chunks(dataset):
    chunk_slices = {}
    for dim, chunks in dataset.chunks.items():
        slices = []
        start = 0
        for chunk in chunks:
            if start >= dataset.sizes[dim]:
                break
            stop = start + chunk
            slices.append(slice(start, stop))
            start = stop
        chunk_slices[dim] = slices
    for slices in itertools.product(*chunk_slices.values()):
        selection = dict(zip(chunk_slices.keys(), slices))
        yield dataset[selection]

def create_filepath_l(ds, prefix='filename', root_path="/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/"):
    """
    Generate a filepath when given an xarray dataset
    """
    start = str(ds.level.values[0])
    start1 = str(int(ds.latitude.data[0]))
    end = str(int(ds.latitude.data[-1]))
    filepath = root_path + prefix + '_' + start + '_' + start1 + '_' + end + '_ERA5_2.nc'
    return filepath

############################################################################
tempdir = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/SMOTHED_CLIM/temp/'
v = var[0]

# ds.to_netcdf(outdir + 'clim_'+v+'_smooth_ERA5_3H_'+str(minYear)+'_'+str(maxYear)+'.nc')
new_ds = xr.open_mfdataset(tempdir+'*.nc', combine='by_coords', parallel=True)
new_ds.to_netcdf(outdir + 'clim_'+v+'_smooth_ERA5_3H_'+str(minYear)+'_'+str(maxYear)+'.nc')
print('script end')
