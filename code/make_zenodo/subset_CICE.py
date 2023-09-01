#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 10:15:29 2023

@author: maxthomas_UO
"""

import xarray as xr
import sys

if __name__ == '__main__':
    lats = [-90.,-60]
    depths = [0,1000]
    variables = ['sidmassth']
    
    infile, outfile = sys.argv[1], sys.argv[2]
        
    raw = xr.open_dataset(infile)[variables]
    print(raw)
    
    raw = raw.where( (raw.TLAT>lats[0]) & (raw.TLAT<lats[1]), drop=True )
    
    print(raw)
    
    raw.to_netcdf(outfile)
    
    
    
    
    
