#!/usr/bin/env python
# coding: utf-8

import xarray as xr
import numpy as np

def amplitudeMagnitude(var1, var2):
    ampl = np.arctan2(var2,var1)  ### arctan2 variable on y first, on x second
    magn = np.sqrt(var1**2 + var2**2)
    
    ampl = ampl.rename('amplitude')
    magn = magn.rename('magnitude').astype('float32')
    # Create new dataset to sctock all the variable
    # _ds = xr.merge([ampl, magn]) 

    return ampl.persist(), magn.persist()

