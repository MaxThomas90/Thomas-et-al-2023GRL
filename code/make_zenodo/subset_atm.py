#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import xarray as xr
import iris
import glob
import os
import sys

def vars2exclude(example_file, vars2keep):
    if type(vars2keep) is str:
        vars2keep = [vars2keep]
    all_vars = xr.open_dataset(example_file)
    drop_vars = []
    for vb in all_vars:
        if not all_vars[vb].name in vars2keep:
            drop_vars.append(all_vars[vb].name)
    return drop_vars

if __name__ == '__main__':
    suite = sys.argv[1]
    data_at = '/nesi/nobackup/nesi00442/thoma97p/cylc-run/' + suite + '/share/data/History_Data/'
    stash_codes = {'m01s03i236':'tas',
                   'm01s05i216':'pr',
                   'm01s03i209':'uwind'}
    
    files = glob.glob(data_at + '*p5*.pp')
    files.sort()
    example_file = files[0]
    #test=iris.load(example_file)#, iris.AttributeConstraint(STASH=stash_codes))
    #print(test)
    for fn in files:
        print(fn[-10:])
        datacubes = iris.load(fn)
        for cube in datacubes:
            # get the STASH code
            cubeSTASH = cube.attributes['STASH']
            if cubeSTASH in list(stash_codes.keys()):
                # Make the output file name
                outfile = '../../data/Zenodo_archive/'+ suite + '/atm/'  + stash_codes[str(cubeSTASH)] + '-' + fn[-10:-3] + '.nc'
                print(outfile)
                # Save the file
                iris.save(cube, outfile)
                print('Saved {}'.format(outfile))


