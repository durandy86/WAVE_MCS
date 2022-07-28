#!/usr/bin/env python
# coding: utf-8

import os, sys
tempDir = os.environ['TMPDIR']
os.environ["MALLOC_TRIM_THRESHOLD_"] = '1'
import numpy as np

import xarray as xr
import pandas as pd

#
from dask_mpi import initialize
initialize(dashboard=False, local_directory = tempDir, nthreads=1, protocol = 'tcp://')

from dask.distributed import Client
client = Client() 

ds = xr.open_mfdataset('/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/SMOTHED_CLIM/' + 'clim_tcwv_smooth_ERA5_3H_1990_2020_3harm.nc', chunks = {'time' : 1}, parallel=True)
ds = ds.sel(latitude = slice(40.1,-40.1))
ds.to_netcdf('/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/SMOTHED_CLIM/' + 'clim_tcwv_smooth_ERA5_3H_1990_2020_3harm_2.nc', unlimited_dims={'time':True})

