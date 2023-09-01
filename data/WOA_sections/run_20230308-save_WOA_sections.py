#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 14:16:15 2023

@author: maxthomas_UO

This script 
1) loads in the files:
    woa18_decav81B0_t00_04.nc
    woa18_decav81B0_s00_04.nc
2) exctracts data for shelf ocean near the Thwaites and Amery glaciers, and 
    for the western Ross (251 E, 71 E, and 177 E, respectively)
3) saves those sections as netcdfs

The original data files are 1981-2010 time means from the 2018 World Ocean 
Atlas on 1/4 degree grids https://www.ncei.noaa.gov/access/world-ocean-atlas-2018/

The profiles are extracted and saved (rather than calculated from raw) as we 
need a very small fraction of the data and moving the whole dataset to the 
supercomputer would be time consuming
"""

import xarray as xr
import os
import shutil
import datetime


if __name__ == '__main__':
    sections = ['Thwaites','Amery','W.Ross']
    depths = [1000, 0]
    extract_properties = {
        'Thwaites': {'lon': 251, 'lats':(-76,-68)},
        'Amery': {'lon': 71, 'lats':(-73,-65)},
        'W.Ross': {'lon': 177, 'lats':(-78.2,-70)}
        }
    
    processed_data_to = '../../data/WOA_sections/'
    raw_data_at = '../../data/WOA_raw/'
    
    os.makedirs(processed_data_to, exist_ok=True)
    for section in sections:
        lats = extract_properties[section]['lats']
        lon = extract_properties[section]['lon']
        lats = extract_properties[section]['lats']
        for vn in ['t','s']:
            loadfile = raw_data_at + 'woa18_decav81B0_' + vn + '00_04.nc'
            savefile = processed_data_to + 'woa-' + section + '-' + vn + '.nc'
            raw_data = xr.open_dataset(loadfile, decode_times=False)[ vn + '_an' ]
            raw_data['lon'] = raw_data['lon'].where(raw_data['lon'] > 0, raw_data['lon'] + 360)
            raw_data = raw_data.where(raw_data.lat<lats[1],drop=True)
            raw_data = raw_data.where(raw_data.lat>lats[0],drop=True)
            raw_data = raw_data.where(raw_data.lon<(lon+1),drop=True)
            raw_data = raw_data.where(raw_data.lon>(lon-1),drop=True)
            raw_data = raw_data.mean(dim='lon')
            
            raw_data.to_netcdf(savefile)
            
            #print(raw_data)
     
    dandt = datetime.datetime.now()
    dandt = str(dandt.year) + str(dandt.month).zfill(2) + str(dandt.day).zfill(2)
    shutil.copy('save_WOA_sections.py', processed_data_to + 'run_'+dandt+'-save_WOA_sections.py')
            
