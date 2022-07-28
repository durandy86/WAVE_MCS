#!/usr/bin/env python
# coding: utf-8

import os, sys
import numpy as np
tempDir = os.environ['TMPDIR']
os.environ["MALLOC_TRIM_THRESHOLD_"] = '0'

import xarray as xr
import pandas as pd


#
from dask_mpi import initialize
initialize(dashboard=False, local_directory = tempDir, nthreads=8, protocol = 'tcp://')

from dask.distributed import Client
client = Client() 
###################################################################################################################

def isLeapYear (yearN):
    if ((yearN % 4 == 0) and (yearN % 100 != 0)) or (yearN % 400 == 0):
        reponse = True
    else:
        reponse = False
    print(reponse, '\n')
    return reponse

def hour_mean(x):
     return x.groupby('time.hour').mean('time')
    
def hour_sum(x):
     return x.groupby('time.hour').sum('time')
    
def hour_std(x):
     return x.groupby('time.hour').std('time')

#######################################################################################################################
var = ['tcwv']
# indir = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/RAW_CLIM/'
indir = '/cnrm/tropics/commun/DATACOMMUN/ERA5/0.25/netcdf/sfc_3h/'
outdir = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/RAW_CLIM/'
minYear = '1990'
maxYear = '2020'

attrs = {"units": "hours since 1900-01-01", "calendar":"gregorian"}
ds_T = xr.Dataset({"time": ("time", np.arange(4*365*24, 4*365*24+366*24, 3), attrs)}) # Reconstruction du calendrier

for v in var :
    ds = xr.open_mfdataset(indir+'*.nc', chunks = {'time' : 10}, parallel=True) #!!!!
    # ds = ds.isel(time = slice(0,None,3)) #pour OLR 1H
    ds = ds.sel(time = slice(minYear,maxYear), latitude = slice(30.1,-30.1)) #Pour période TOUCAN 
 
    # group by dayofyear, then apply the function:
    ds = ds.groupby('time.dayofyear').apply(hour_mean)
    # ds_clim.to_netcdf(outdir + 'clim_'+v+'_brut_ERA5_2000_2020.nc')
    #### Reconstruction of the data
    da = xr.DataArray(
        data = np.reshape(ds[v].values,(np.shape(ds.dayofyear)[0]*np.shape(ds.hour)[0],np.shape(ds.latitude)[0],np.shape(ds.longitude)[0])),
            dims = ["time","latitude","longitude"],
            coords=dict(
            time = ds_T.time,
            latitude = ds.latitude,
            longitude = ds.longitude
            ),
            attrs=dict(
            description="tcwv",
            units="kg m-2"
            ),
        )
    ds = da.to_dataset(name = 'tcwv')
    ds = ds.astype(np.float32)
    ds = ds.compute()
    ds.to_netcdf(outdir + 'clim_'+v+'_brut_ERA5_3H_'+str(minYear)+'_'+str(maxYear)+'.nc', unlimited_dims={'time':True}, mode = 'w')

print('end of script')


