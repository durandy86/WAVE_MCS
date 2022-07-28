#!/usr/bin/env python
# coding: utf-8

import os, sys
tempDir = os.environ['TMPDIR']
os.environ["MALLOC_TRIM_THRESHOLD_"] = '5'

import numpy as np
import h5netcdf

import xarray as xr
import xarray.ufuncs as xu
import xrft
import pandas as pd


from dask.distributed import Client, LocalCluster
#
# Initialisation d'un cluster de 32 coeurs
cluster = LocalCluster(processes=False, n_workers=1, threads_per_worker=32, silence_logs='error', local_directory = tempDir)
client = Client(cluster)
client


############################################################################################
############################################################################################
############################################################################################
var = 'TCWV'
units = "kg/m^2"
filenames = np.arange(2019,2020)

addDay = 180
spd = 8

indir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/RAW_ANOMALY/TCWV/'
outdir_data = '/cnrm/tropics/commun/DATACOMMUN/WAVE/NO_SAVE/DATA/FILTERED_ANOMALY/TCWV/'

var_file = 'anom_tcwv_brut_ERA5_3H_'


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

def kelvinfilter(da, fmin=None, fmax=None, kmin=None, kmax=None, hmin=None, hmax=None, rlat = 0):
    """kelvin wave filter
    Arguments:
       'fmin/fmax' -- unit is cycle per day
       'kmin/kmax' -- zonal wave number
       'hmin/hmax' --equivalent depth
    """
    knum = da.freq_lon
    freq = da.freq_time
    
    # filtering ############################################################
    mask = da.copy()
    #wavenumber cut-off
    if kmin != None:
        mask = mask.where((da.freq_lon > kmin) & (da.freq_lon < kmax), drop = False)
    
    #frequency cutoff
    if fmin != None:
        mask = mask.where((da.freq_time > fmin) & (da.freq_time < fmax), drop = False)

    #dispersion filter
    if hmin != None:
        c = np.sqrt(g*hmin)
        omega = -2.*pi*freq/24./3600. / np.sqrt(Beta*c) #adusting day^-1 to s^-1
        k     = 1/((radius)/(knum-1)) * np.sqrt(c/Beta)         #adusting wavenumber to m
        mask = mask.where(lambda da: -(omega - k) <= 0., drop = False)
    if hmax != None:
        c = np.sqrt(g*hmax)
        omega = -2.*pi*freq/24./3600. / np.sqrt(Beta*c) #adusting day^-1 to s^-1
        k     = 1/((radius)/(knum+1)) * np.sqrt(c/Beta)         #adusting wavenumber to m
        mask = mask.where(lambda da: -(omega - k) >= 0., drop = False)

    _varhat = mask
    _varhat['freq_lon'] = freq_lon_Save
    _varhat['freq_time'] = freq_time_Save
    filterd = xrft.ifft(_varhat.fillna(0.), dim = ['freq_time','freq_lon'], 
                                  true_phase=False, true_amplitude=True) # Signal in direct space
    filterd['time'] = ds['time']
    return mask.persist(), filterd

def mjofilter(da, fmin=None, fmax=None, kmin=None, kmax=None) :
    """kelvin wave filter
    Arguments:
       'fmin/fmax' -- unit is cycle per day
       'kmin/kmax' -- zonal wave number
       'hmin/hmax' --equivalent depth
    """
    
    knum = da.freq_lon
    freq = da.freq_time

    # filtering ############################################################
    mask = da.copy()
    #wavenumber cut-off
    if kmin != None:
        mask = mask.where((da.freq_lon > kmin) & (da.freq_lon < kmax), drop = False)
    
    #frequency cutoff
    if fmin != None:
        mask = mask.where((da.freq_time > fmin) & (da.freq_time < fmax), drop = False)

    _varhat = mask.compute()
    _varhat['freq_lon'] = freq_lon_Save
    _varhat['freq_time'] = freq_time_Save
    filterd = xrft.ifft(_varhat.fillna(0.), dim = ['freq_time','freq_lon'], 
                                  true_phase=False, true_amplitude=True)  # Signal in direct space
    filterd['time'] = ds['time']
    return mask, filterd

def erfilter(da, fmin=None, fmax=None, kmin=None, kmax=None, hmin=None, hmax=None, n=1, rlat = 0):
    """equatorial wave filter
    Arguments:
       'fmin/fmax' -- unit is cycle per day
       'kmin/kmax' -- zonal wave number
       'hmin/hmax' -- equivalent depth
       'n'         -- meridional mode number
    """

    if n <=0 or n%1 !=0:
        print ("n must be n>=1 integer")
        sys.exit()

    knum = da.freq_lon
    freq = da.freq_time

    
    # filtering ############################################################
    mask = da.copy()
    #wavenumber cut-off
    if kmin != None:
        mask = mask.where((da.freq_lon > kmin) & (da.freq_lon < kmax), drop = False)
    
    #frequency cutoff
    if fmin != None:
        mask = mask.where((da.freq_time > fmin) & (da.freq_time < fmax), drop = False)
        
        
    # filtering ############################################################
    if hmin != None:
        c = np.sqrt(g*hmin)
        omega = -2.*pi*freq/24./3600. / np.sqrt(Beta*c)           #adusting day^-1 to s^-1
        k     = 1/((radius)/(knum+1)) * np.sqrt(c/Beta)                          #adusting ^2pia to ^m               
        mask = mask.where(lambda da: -(omega*(k**2 + (2*n+1)) + k) <= 0., drop = False)   
    if hmax != None:
        c = np.sqrt(g*hmax)
        omega = -2.*pi*freq/24./3600. / np.sqrt(Beta*c)           #adusting day^-1 to s^-1
        k     = 1/((radius)/(knum-1)) * np.sqrt(c/Beta)   
        mask = mask.where(lambda da: -(omega*(k**2 + (2*n+1)) + k) >= 0., drop = False)             

    _varhat = mask
    _varhat['freq_lon'] = freq_lon_Save
    _varhat['freq_time'] = freq_time_Save
    filterd = xrft.ifft(_varhat.fillna(0.), dim = ['freq_time','freq_lon'], 
                                  true_phase=False, true_amplitude=True) # Signal in direct space
    filterd['time'] = ds['time']
    return mask.persist(), filterd

def mrgfilter(da, fmin=None, fmax=None, kmin=None, kmax=None, hmin=None, hmax=None,):
    """mixed Rossby gravity wave
    Arguments:
       'fmin/fmax' -- unit is cycle per day
       'kmin/kmax' -- zonal wave number. negative is westward, positive is
                      eastward
       'hmin/hmax' -- equivalent depth
    """
    
    knum = da.freq_lon
    freq = da.freq_time
    
        
    # filtering ############################################################
    mask = da.copy()
    #wavenumber cut-off
    if kmin != None:
        mask = mask.where((da.freq_lon > kmin) & (da.freq_lon < kmax), drop = False)
    
    #frequency cutoff
    if fmin != None:
        mask = mask.where((da.freq_time > fmin) & (da.freq_time < fmax), drop = False)
        
        
    #dispersion filter
    mask = mask.where((da.freq_lon <= -1.), drop = False)
    if hmin != None:
        c = np.sqrt(g*hmin)
        omega = -2.*pi*freq/24./3600. / np.sqrt(Beta*c)           #adusting day^-1 to s^-1
        k     = 1/((radius)/(knum-1)) * np.sqrt(c/Beta)                         #adusting ^2pia to ^m               
        mask = mask.where(lambda da: -(omega**2 - k*omega - 1) <= 0., drop = False) 
    if hmax != None:
        c = np.sqrt(g*hmax)
        omega = -2.*pi*freq/24./3600. / np.sqrt(Beta*c)           #adusting day^-1 to s^-1
        k     = 1/((radius)/(knum+1)) * np.sqrt(c/Beta)                 #adusting ^2pia to ^m               
        mask = mask.where(lambda da: -(omega**2 - k*omega - 1) >= 0., drop = False) 
  
    _varhat = mask
    _varhat['freq_lon'] = freq_lon_Save
    _varhat['freq_time'] = freq_time_Save
    filterd = xrft.ifft(_varhat.fillna(0.), dim = ['freq_time','freq_lon'], 
                                  true_phase=False, true_amplitude=True) # Signal in direct space
    filterd['time'] = ds['time']
    return mask.persist(), filterd

def igfilter(da, fmin=None, fmax=None, kmin=None, kmax=None, hmin=None, hmax=None, n=1):
    """n>=1 inertio gravirt wave filter. default is n=1 WIG.
    Arguments:
       'fmin/fmax' -- unit is cycle per day
       'kmin/kmax' -- zonal wave number. negative is westward, positive is
                      eastward
       'hmin/hmax' -- equivalent depth
       'n'         -- meridional mode number
    """

    # filtering ############################################################
    mask = da.copy()
    #wavenumber cut-off
    if kmin != None:
        mask = mask.where((da.freq_lon > kmin) & (da.freq_lon < kmax), drop = False)
    
    #frequency cutoff
    if fmin != None:
        mask = mask.where((da.freq_time > fmin) & (da.freq_time < fmax), drop = False)
        
    knum = da.freq_lon
    freq = da.freq_time
    #dispersion filter
    if hmin != None:
        c = np.sqrt(g*hmin)
        omega = -2.*pi*freq/24./3600. / np.sqrt(Beta*c)           #adusting day^-1 to s^-1
        k     = 1/((radius)/(knum-1)) * np.sqrt(c/Beta)                   #adusting ^2pia to ^m               
        mask = mask.where(lambda da: -(omega**2 - k**2 - (2*n+1)) <= 0., drop = False) 
    if hmax != None:
        c = np.sqrt(g*hmax)
        omega = -2.*pi*freq/24./3600. / np.sqrt(Beta*c)           #adusting day^-1 to s^-1
        k     = 1/((radius)/(knum-1)) * np.sqrt(c/Beta)                      #adusting ^2pia to ^m               
        mask = mask.where(lambda da: -(omega**2 - k**2 - (2*n+1)) >= 0., drop = False) 
        
        
    _varhat = mask
    _varhat['freq_lon'] = freq_lon_Save
    _varhat['freq_time'] = freq_time_Save
    filterd = xrft.ifft(_varhat.fillna(0.), dim = ['freq_time','freq_lon'], 
                                  true_phase=False, true_amplitude=True) # Signal in direct space
    filterd['time'] = ds['time']
    return mask.persist(), filterd


def deletebug(da, ds):
    _da = da.real.values + xu.conj(da).real.values
    da = xr.DataArray(_da,
                       dims=("time","lat","lon"), 
                       coords={"time":ds.time,
                               "lat":ds.lat,
                               "lon":ds.lon})
      
    return da.persist()

#####################################
#####################################

for f in filenames:
    date_Re = pd.date_range(start = (str(f) + '-01-01'), end = (str(f + 1) + '-01-01'), freq = '3H', closed = 'left')

    ds = createArray(f)
    # Part for the test
    # date_Re = pd.date_range(start = (str(f) + '-01-01'), end = (str(f + 1) + '-01-01'), freq = '24H', closed = 'left')
    # ds = ds.isel(time = slice(0,None,spd)).load()
    # End part for the test

    ds = ds.rename({'latitude':'lat','longitude':'lon'})
    x_wintap = ds['tcwv_ano'].chunk({"time" : -1, "lat": 1})

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

    for j in range(0,ds.lat.size,50):
        _ds = ds.isel(lat = slice(j,j+50))
        _tcwvhat = tcwvhat.isel(lat = slice(j,j+50))

        ### Delete LF
        _dahat_LF, da_LFfilter = mjofilter(_tcwvhat, fmin = -1/96, fmax= -1/999, kmax = 10, kmin = -10)
        _tcwvhat = _tcwvhat - _dahat_LF.fillna(0.)
        ###

        #####################################
        # Kelvin Filtering
        # kelvin: (Straub Kiladis, Kiladis 2009)
        tMin = 2.5
        tMax = 30.
        kMin = 1
        kMax = 14
        hMin = 8
        hMax = 90


        _dahat_Kev, da_Kevfilter = kelvinfilter(_tcwvhat, fmin = -1/tMin, fmax = -1/tMax, kmax = kMax, kmin = 1, hmin=8, hmax=90)
        da_Kevfilter = deletebug(da_Kevfilter, _ds)
        
        da_Kevfilter.attrs['units'] = units
        da_Kevfilter.attrs['tMin'] = str(tMin)
        da_Kevfilter.attrs['tMax'] = str(tMax)
        da_Kevfilter.attrs['kMin'] = str(kMin)
        da_Kevfilter.attrs['kMax'] = str(kMax)
        da_Kevfilter.attrs['hMin'] = str(hMin)
        da_Kevfilter.attrs['hMax'] = str(hMax)

        ######################################
        # Rossby Filtering
        # n=1 Eq. Rossby  Kiladis et al 2009
        tMin = 9.7
        tMax = 72.
        kMin = -10
        kMax = -1
        hMin = 8 
        hMax = 90

        _dahat_Ross, da_Rossfilter = erfilter(_tcwvhat, fmin = -1/tMin, fmax = -1/tMax, kmax = kMax, kmin = kMin, hmin=hMin, hmax=90)
        da_Rossfilter = deletebug(da_Rossfilter, _ds)
        
        da_Rossfilter.attrs['units'] = units
        da_Rossfilter.attrs['tMin'] = str(tMin)
        da_Rossfilter.attrs['tMax'] = str(tMax)
        da_Rossfilter.attrs['kMin'] = str(kMin)
        da_Rossfilter.attrs['kMax'] = str(kMax)
        da_Rossfilter.attrs['hMin'] = str(hMin)
        da_Rossfilter.attrs['hMax'] = str(hMax)
        
        #########################################
        # MJO Filtering
        # MJO
        # waveName = "MJO"   ; kiladis et al. 2009
        tMin = 30.
        tMax = 96.
        kMin = 1
        kMax = 5  # 5 is also possible

        _dahat_MJO, da_MJOfilter = mjofilter(_tcwvhat, fmin = -1/tMin, fmax= -1/tMax, kmax = kMax, kmin = kMin)
        da_MJOfilter = deletebug(da_MJOfilter, _ds)

        da_MJOfilter.attrs['units'] = units
        da_MJOfilter.attrs['tMin'] = str(tMin)
        da_MJOfilter.attrs['tMax'] = str(tMax)
        da_MJOfilter.attrs['kMin'] = str(kMin)
        da_MJOfilter.attrs['kMax'] = str(kMax)

        #############################################
        # MRG Filtering
        # MRG
        tMin = 3.
        tMax = 8.
        kMin = -10
        kMax = -1
        hMin = 8
        hMax = 90

        _dahat_MRG, da_MRGfilter = mrgfilter(_tcwvhat, fmin = -1/3, fmax = -1/8, kmax = -1, kmin = -10, hmin = 8, hmax = 90)
        da_MRGfilter = deletebug(da_MRGfilter, _ds)

        da_MRGfilter.attrs['units'] = units
        da_MRGfilter.attrs['tMin'] = str(tMin)
        da_MRGfilter.attrs['tMax'] = str(tMax)
        da_MRGfilter.attrs['kMin'] = str(kMin)
        da_MRGfilter.attrs['kMax'] = str(kMax)
        da_MRGfilter.attrs['hMin'] = str(hMin)
        da_MRGfilter.attrs['hMax'] = str(hMax)

        ################################################
        # IGW Filtering
        # EIG  (Wheeler Kiladis 1998)    
        tMin = 1.
        tMax = 5.
        kMin = 14
        kMax = 1
        hMin = 12  
        hMax = 50
       
        _dahat_EIG, da_EIGfilter = igfilter(_tcwvhat, fmin = -1, fmax = -1/5, kmin=1, kmax=14, hmin=12, hmax = hMax)
        da_EIGfilter = deletebug(da_EIGfilter, _ds)

        da_EIGfilter.attrs['units'] = units
        da_EIGfilter.attrs['tMin'] = str(tMin)
        da_EIGfilter.attrs['tMax'] = str(tMax)
        da_EIGfilter.attrs['kMin'] = str(kMin)
        da_EIGfilter.attrs['kMax'] = str(kMax)
        da_EIGfilter.attrs['hMin'] = str(hMin)
        da_EIGfilter.attrs['hMax'] = str(hMax)

        # WIG  (Wheeler Kiladis 1998)
        tMin = 1.
        tMax = 5.
        kMin = -14
        kMax = -1
        hMin = 12 
        hMax = 50
        
        _dahat_WIG, da_WIGfilter = igfilter(_tcwvhat, fmin = -1/tMin, fmax = -1/tMax, kmin=-14, kmax=-1, hmin=12, hmax=50)
        da_WIGfilter = deletebug(da_WIGfilter, _ds)
        
        da_WIGfilter.attrs['units'] = units
        da_WIGfilter.attrs['tMin'] = str(tMin)
        da_WIGfilter.attrs['tMax'] = str(tMax)
        da_WIGfilter.attrs['kMin'] = str(kMin)
        da_WIGfilter.attrs['kMax'] = str(kMax)
        da_WIGfilter.attrs['hMin'] = str(hMin)
        da_WIGfilter.attrs['hMax'] = str(hMax)
        ######################################################
        # TD filtering
        
        # TD (Easterly wave) : 
        tMin = 2.5
        tMax = 7.
        kMin = -20
        kMax = -6
        # hMin = mis  ==> Min curve = up curve of MRG Domain
        # hMax = mis
        _dahat_TD, da_TDfilter = mrgfilter(_tcwvhat, fmin = -1/2.5, fmax = -1/7, kmax = -6, kmin = -20, hmin = 100)
        da_TDfilter = deletebug(da_TDfilter, _ds)
        
        da_TDfilter.attrs['units'] = units
        da_TDfilter.attrs['tMin'] = str(tMin)
        da_TDfilter.attrs['tMax'] = str(tMax)
        da_TDfilter.attrs['kMin'] = str(kMin)
        da_TDfilter.attrs['kMax'] = str(kMax)
        """
        ########################################################
        # Low Frequency Filtering
        # Low frequency 
        tMin = 120.
        tMax = 999
        kMin = -10
        kMax = 10
        
        _dahat_LF, da_LFfilter = mjofilter(_tcwvhat, fmin = -1/120, fmax= -1/999, kmax = 10, kmin = -10)
        da_LFfilter = deletebug(da_LFfilter, _ds)

        da_LFfilter.attrs['units'] = units
        da_LFfilter.attrs['tMin'] = str(tMin)
        da_LFfilter.attrs['tMax'] = str(tMax)
        da_LFfilter.attrs['kMin'] = str(kMin)
        da_LFfilter.attrs['kMax'] = str(kMax)
        """
  
        #######################################################
        #######################################################
        # NetCDF save
        ds_NetCdf = da_Kevfilter.to_dataset(name = var + '_Kelvin')
        ds_NetCdf[var + '_Rossby'] = da_Rossfilter
        ds_NetCdf[var + '_MJO'] = da_MJOfilter
        ds_NetCdf[var + '_MRG'] = da_MRGfilter
        ds_NetCdf[var + '_EIG'] = da_EIGfilter
        ds_NetCdf[var + '_WIG'] = da_WIGfilter
        ds_NetCdf[var + '_TD'] = da_TDfilter
        # ds_NetCdf[var + '_LF'] = da_LFfilter 

        ds_NetCdf.time.encoding['units'] = "hours since 1900-01-01"
        ds_NetCdf.time.encoding['calendar'] = "proleptic_gregorian"
        ds_NetCdf.time.encoding['standard_name'] = "time"
        ds_NetCdf = ds_NetCdf.sel(time = str(f))
        ds_NetCdf.to_netcdf(outdir_data + 'LAT/TCWV_ERA5_3H_' + str(j) + '_F.nc', unlimited_dims='time', mode = 'w')
#     ds_NetCdf = ds_NetCdf.sel(time = str(f))
#     ds_NetCdf.to_netcdf(outdir_data + 'TCWV_ERA5_3H_' + str(f) +'.nc', unlimited_dims='time')
    ds_NetCdf = ds_NetCdf.compute()
    ds_NetCdf = xr.open_mfdataset(outdir_data + 'LAT/TCWV_ERA5_3H_*_F.nc', parallel=True)
    ds_NetCdf.to_netcdf(outdir_data + 'TCWV_FILTER_ERA5_3H_' + str(f) + '.nc', unlimited_dims='time', mode = 'w')
print('script end')

