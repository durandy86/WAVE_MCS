#!/usr/bin/env python
# coding: utf-8

import os, sys
tempDir = os.environ['TMPDIR']
os.environ["MALLOC_TRIM_THRESHOLD_"] = '0'
import numpy as np

import xarray as xr
import xrft
import pandas as pd

from mpi4py import MPI

from dask_mpi import initialize
initialize(dashboard=False, local_directory = tempDir, nthreads=2, protocol = 'tcp://')

from distributed import Client
client = Client()

################################################################################

def isLeapYear (yearN):
    if ((yearN % 4 == 0) and (yearN % 100 != 0)) or (yearN % 400 == 0):
        reponse = True
    else:
        reponse = False
    print(reponse, '\n')
    return reponse
    
#################################################################################

indir_clim = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/SMOTHED_CLIM/'

## Pour TCWV
indir_data = '/cnrm/tropics/commun/DATACOMMUN/ERA5/0.25/netcdf/sfc_3h/'
outdir = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/RAW_ANOMALY/TCWV/'
var = 'tcwv'

ds_data = xr.open_mfdataset(indir_data+'*'+var+'*.nc', chunks = {'time' : 150}, parallel=True)
ds_data = ds_data.sel(latitude = slice(40.1,-40.1))
ds_clim = xr.open_dataset(indir_clim + 'clim_tcwv_smooth_ERA5_3H_1990_2020.nc', chunks = {'time' : 150})

v = var
year = np.arange(1990,2021)

for y in year :
    print(y)
    ds = ds_data.sel(time = str(y))
    if isLeapYear(ds.time.dt.year.values[0]) == False :
        _ds_clim = ds_clim.isel(time = slice(0,365*int(24/3)))
        _ds_clim['time'] = ds['time']
    else :
        _ds_clim = ds_clim
        _ds_clim['time'] = ds['time']

    da = ds[v] - _ds_clim[v + '_smooth']
    ### add on 10 november
    # da = da - (ds[v].mean('time') - _ds_clim[v + '_smooth'].mean('time'))
    ###
    da = da.to_dataset(name = "tcwv_ano")
    da = da.compute()
    da.to_netcdf(outdir + 'anom_tcwv_brut_ERA5_3H_'+str(y) + '.nc', unlimited_dims={'time':True}, mode = 'w')
