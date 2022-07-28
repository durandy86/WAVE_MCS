#!/usr/bin/env python
# coding: utf-8

import os, sys
tempDir = os.environ['TMPDIR']
os.environ["MALLOC_TRIM_THRESHOLD_"] = '5'
tempdir_z = tempDir + '/level_CLIM/'

import numpy as np
import h5netcdf

import xarray as xr
import xarray.ufuncs as xu
import xrft
import pandas as pd
from scipy.signal import convolve2d, detrend

from definition import *
from waveFiltering import *

#from dask.distributed import Client, LocalCluster
#
# Initialisation d'un cluster de 32 coeurs
#cluster = LocalCluster(processes=False, n_workers=1, threads_per_worker=32, silence_logs='error', local_directory = tempDir)
#client = Client(cluster)
#client

from dask_mpi import initialize
initialize(dashboard=False, local_directory = tempDir, nthreads=8, protocol = 'tcp://')

from dask.distributed import Client
client = Client()


############################################################################################
############################################################################################
############################################################################################
var = 'vo'
units = "m/s^1"
filenames = np.arange(2009,2010)
varF = 'vo'
prefix = 'FILTEREDANOM_' + varF + '_'
zLevel = [850]

addDay = 180
spd = 8

indir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/RAW_ANOMALY/' + var +'/'
outdir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/FILTERED_ANOMALY/' + var + '/'

var_file = 'anom_u_z850_brut_ERA5_3H_'


# Constante
rlat = 0
pi    = np.pi
radius = 6.37122e06    # [m]   average radius of earth
g     = 9.80665        # [m/s] gravity at 45 deg lat used by the WMO
omega = 7.292e-05      # [1/s] earth's angular vel
ll    = 2.*pi*radius*np.cos(np.abs(rlat))
Beta  = 2.*omega*np.cos(np.abs(rlat))/radius
fillval = 1e20

###################################################################################
def createArray(year) :
    _ds_m1 = xr.open_mfdataset(indir_data+'*'+var_file+'*'+str(year-1)+'.nc', chunks = {'lat' : 1}, parallel=True)
    _ds_m1 = _ds_m1.isel(time = slice(-addDay*spd,None))
    _ds = xr.open_mfdataset(indir_data+'*'+var_file+'*'+str(year)+'.nc', chunks = {'lat' : 1}, parallel=True)
    _ds1 = xr.open_mfdataset(indir_data+'*'+var_file+'*'+str(year+1)+'.nc', chunks = {'lat' : 1}, parallel=True)
    _ds1 = _ds1.isel(time = slice(None,addDay*spd))

    ds = xr.concat([_ds_m1,_ds,_ds1], dim='time', coords='minimal', compat='override')
    
    return ds

#####################################
#####################################
for f in filenames:
    date_Re = pd.date_range(start = (str(f) + '-01-01'), end = (str(f + 1) + '-01-01'), freq = '3H', closed = 'left')
    ds = createArray(f)
    ds = ds.rename({'latitude':'lat','longitude':'lon'})
    # Part for the test
    # date_Re = pd.date_range(start = (str(f) + '-01-01'), end = (str(f + 1) + '-01-01'), freq = '24H', closed = 'left')
    # ds = ds.isel(time = slice(0,None,spd), lon = slice(0,None,4),lat = slice(0,None,4))
    # End part for the test
###########################################################################
###########################################################################
    for z in zLevel :
        print("Level = ", z,"\n")
        _ds = ds.sel(level = z)
        _ds = _ds.expand_dims(level = 1)
        x_wintap = _ds[var + '_ano'].chunk({"time" : -1, "lat": 1})  #!!!!!!!!!!!!!!!!!!!!!!! Variable to change

        tcwvhat = xrft.fft(x_wintap, detrend='linear',
                    dim=['time','lon'], true_phase=False, true_amplitude=True)


        #####################################
        ### We save freq_lon in ° and freq_time in second, because we use wavenumber and cycle by day for filtering
        freq_lon_Save = tcwvhat['freq_lon']
        freq_time_Save = tcwvhat['freq_time']
        #####################################
        wavenumber = np.zeros(tcwvhat.freq_lon.size)
        for i in range( tcwvhat.freq_lon.size):
            j= - int(360/2) + i
            wavenumber[i] = tcwvhat.freq_lon[int(360/2)+j]*360 + 1
     
        tcwvhat['freq_lon'] = wavenumber
        tcwvhat['freq_time'] = tcwvhat.freq_time*86400
        print(_ds['time'])
        __ds = filterWave(_ds, tcwvhat, var, freq_lon_Save, freq_time_Save)
        
        __ds = __ds.sel(time = str(f)).astype(np.float32)
        __ds.to_netcdf(tempdir_z + prefix + str(z) + '.nc') #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!, engine = 'h5netcdf'
        del __ds
        #xr.save_mfdataset(datasets=datasets, paths=paths_l)
    # ds.to_netcdf(outdir + 'clim_'+v+'_smooth_ERA5_3H_'+str(minYear)+'_'+str(maxYear)+'.nc')
    ds = xr.open_mfdataset(tempdir_z + prefix + '*.nc' , combine='by_coords', parallel=True)
    ds.to_netcdf(outdir_data + varF + '_z850_FILTER_ERA5_3H_' + str(f) + '.nc', unlimited_dims='time', mode = 'w')

print('script end')

