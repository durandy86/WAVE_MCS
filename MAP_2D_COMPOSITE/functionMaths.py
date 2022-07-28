#!/usr/bin/env python
# coding: utf-8

import xarray as xr
import numpy as np

def amplitudeMagnitude(var1, var2):
    ampl = np.arctan2(var2,var1)  ### arctan2 variable on y first, on x second
    magn = np.sqrt(var1**2 + var2**2)
    
    ampl = ampl.rename('amplitude')
    magn = magn.rename('magnitude')
    # Create new dataset to sctock all the variable
    _ds = xr.merge([ampl, magn]) 

    return _ds

outdir = '/cnrm/tropics/commun/DATACOMMUN/WAVE/SCRIPTS/NETcdf/composite/'
def composite_200(name_wave, name_v, ds_phase):
    if name_v != 'tcwv' :
        indir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/RAW_ANOMALY/' + name_v + '/reduce/'
        ds_v = xr.open_mfdataset(indir_data + '*z200*.nc', chunks = {'time' : 1}, parallel = True)
    else :
        indir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/RAW_ANOMALY/TCWV/reduce/'
        ds_v = xr.open_mfdataset(indir_data + 'TCWV_1990_2020.nc', chunks = {'time' : 1}, parallel = True)
    ds_v = ds_v.sel(time = slice('2001','2006'))    
    #########################################
    if (name_v == 'd' or name_v == 't'):
        ds_v = ds_v[name_v + '_ano'].isel(level = 0).expand_dims(phase = ds_phase['phase'], axis = 3)
    else : 
        ds_v = ds_v[name_v + '_ano'].expand_dims(phase = ds_phase['phase'], axis = 3)
    ds_v = xr.where(ds_phase.isnull() != True,
                             ds_v, np.nan)

    ds_v = ds_v.mean(dim = ['time','lon']).compute()

    ds_v.name = name_v
    ds_dyn = ds_v.to_dataset()
    ds_dyn.to_netcdf(outdir + name_v + '_ano_z200_' + name_wave + '_NAKA.nc')

def composite_850(name_wave, name_v, ds_phase):
    indir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/RAW_ANOMALY/' + name_v + '/reduce/'
    ds_v = xr.open_mfdataset(indir_data + '*z850*.nc', chunks = {'time' : 1}, parallel = True)
    ds_v = ds_v.sel(time = slice('2001','2006'))
    #########################################
    ds_v = ds_v[name_v + '_ano'].isel(level = 0).expand_dims(phase = ds_phase['phase'], axis = 3)
    ds_v = xr.where(ds_phase.isnull() != True,
                             ds_v, np.nan)

    ds_v = ds_v.mean(dim = ['time','lon']).compute()

    ds_v.name = name_v
    ds_dyn = ds_v.to_dataset()

    ds_dyn.to_netcdf(outdir + name_v + '_ano_z850_' + name_wave + '_NAKA.nc')

