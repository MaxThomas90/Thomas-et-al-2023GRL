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
    variables = ['thetao','so','soicecov','sowflisf','opottempadvect','friver']
    
    infile, outfile = sys.argv[1], sys.argv[2]
        
    raw = xr.open_dataset(infile)[variables]
    
    raw = raw.where( (raw.nav_lat>lats[0]) & (raw.nav_lat<lats[1]), drop=True )
    raw = raw.sel(deptht=slice(0,1000))
    
    
    raw.to_netcdf(outfile)
    
    
    
    
    
