#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 16:10:17 2023

@author: maxthomas_UO
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.interpolate import griddata
import xarray as xr
import cmocean

if __name__ == '__main__':
    
    matplotlib.rcParams.update({'font.size': 14}) 
    
    def regrid_cmip_to_woa(cmip, woa, lats):
        shape1 = cmip.shape
        shape2 = woa.shape
        xcoords1 = lats
        xcoords2 = woa.lat
        ycoords1 = cmip.lev
        ycoords2 = woa.depth
        X1, Y1 = np.meshgrid(xcoords1, ycoords1)
        X2, Y2 = np.meshgrid(xcoords2, ycoords2)
        points_1 = np.column_stack((X1.ravel(), Y1.ravel()))
        points_2 = np.column_stack((X2.ravel(), Y2.ravel()))
        data1 = cmip.values
        data2 = woa.values
        regridded_data = griddata(points_1, data1.ravel(), points_2, method='linear')        
        regridded_data = regridded_data.reshape(shape2)
        mask = np.isnan(data2)
        difference = data2 - regridded_data
        difference[mask] = np.nan
        to_plot = {'cmip': [X1,Y1,data1],
                   'woa': [X2,Y2,data2],
                   'woa-cmip': [X2,Y2,difference]}
        return to_plot
    
    def fmt(x):
        s = f"{x:.1f}"
        if s.endswith("0"):
            s = f"{x:.0f}"
        return rf"{s} " if plt.rcParams["text.usetex"] else f"{s} "
        
    
    WOA = {'to': xr.open_dataset('ISMIPfw-obs_comparison-woa_to.nc', decode_times=False).s_an,
           'so': xr.open_dataset('ISMIPfw-obs_comparison-woa_so.nc', decode_times=False).s_an,
           'sigma0': xr.open_dataset('ISMIPfw-obs_comparison-woa_sigma0.nc', decode_times=False).s_an}
    
    CMIP = {'to': xr.open_dataset('ISMIPfw-obs_comparison-cmip_to.nc', decode_times=False).thetao,
            'so': xr.open_dataset('ISMIPfw-obs_comparison-cmip_so.nc', decode_times=False).so,
            'sigma0': xr.open_dataset('ISMIPfw-obs_comparison-cmip_sigma0.nc', decode_times=False).so}
    
    lats = xr.open_dataset('ISMIPfw-obs_comparison-lats.nc').lat.values
    
    
    fig = plt.figure(dpi=300,figsize=(12,18))
    levels = [27.4,27.6,27.8,28]    
    gs = fig.add_gridspec(3,2)
    
    to_plot = {
            'to': regrid_cmip_to_woa(CMIP['to'], WOA['to'], lats),
            'so': regrid_cmip_to_woa(CMIP['so'], WOA['so'], lats),
            'sigma0': regrid_cmip_to_woa(CMIP['sigma0'], WOA['sigma0'], lats)
              }
    
    axes_so = {}
    for iax, ax in enumerate(['WOA', 'cmip', 'WOA-cmip']):
        axes_so[ax] = fig.add_subplot(gs[iax,0])
        axes_so[ax].set_ylim([1000,0])
        axes_so[ax].set_xlim([-78,-70])
        axes_so[ax].set_facecolor((0.8,0.8,0.8))
        axes_so[ax].set_ylabel(r'depth, $z$ / m')
        if iax<2:
            axes_so[ax].set_xticklabels([])
        else:
            axes_so[ax].set_xlabel(r'latitude, $\phi$ / $^\circ$N')

    axes_to = {}
    for iax, ax in enumerate(['WOA', 'cmip', 'WOA-cmip']):
        axes_to[ax] = fig.add_subplot(gs[iax,1]) 
        axes_to[ax].set_ylim([1000,0])
        axes_to[ax].set_xlim([-78,-70])
        axes_to[ax].set_facecolor((0.8,0.8,0.8))
        axes_to[ax].set_yticklabels([])
        if iax<2:
            axes_to[ax].set_xticklabels([])
        else:
            axes_to[ax].set_xlabel(r'latitude, $\phi$ / $^\circ$N')
        
    cm = axes_so['WOA'].pcolormesh(to_plot['so']['woa'][0], to_plot['so']['woa'][1], to_plot['so']['woa'][2], vmin=33., vmax=35., cmap=cmocean.cm.haline)
    sig = axes_so['WOA'].contour(to_plot['sigma0']['woa'][0], to_plot['sigma0']['woa'][1], to_plot['sigma0']['woa'][2], levels=levels, colors='k')
    axes_so['WOA'].clabel(sig, sig.levels, inline=True, fmt=fmt, fontsize=9) 
    fig.colorbar(cm, ax=axes_so['WOA'], location='right', label=r'$S$ / g/kg')
    
    cm = axes_so['cmip'].pcolormesh(to_plot['so']['cmip'][0], to_plot['so']['cmip'][1], to_plot['so']['cmip'][2], vmin=33., vmax=35., cmap=cmocean.cm.haline)
    sig = axes_so['cmip'].contour(to_plot['sigma0']['cmip'][0], to_plot['sigma0']['cmip'][1], to_plot['sigma0']['cmip'][2], levels=levels, colors='k', linestyles=[':'])
    axes_so['cmip'].clabel(sig, sig.levels, inline=True, fmt=fmt, fontsize=9) 
    fig.colorbar(cm, ax=axes_so['cmip'], location='right', label=r'$S$ / g/kg')
    
    
    cm = axes_so['WOA-cmip'].pcolormesh(to_plot['so']['woa-cmip'][0], to_plot['so']['woa-cmip'][1], to_plot['so']['woa-cmip'][2], vmin=-1, vmax=1, cmap=cmocean.cm.balance)
    fig.colorbar(cm, ax=axes_so['WOA-cmip'], location='right', label=r'$\Delta S$ / g/kg')
        
    cm = axes_to['WOA'].pcolormesh(to_plot['to']['woa'][0], to_plot['to']['woa'][1], to_plot['to']['woa'][2], vmin=-2.5, vmax=2.5, cmap=cmocean.cm.delta)
    sig = axes_to['WOA'].contour(to_plot['sigma0']['woa'][0], to_plot['sigma0']['woa'][1], to_plot['sigma0']['woa'][2], levels=levels, colors='k')
    axes_to['WOA'].clabel(sig, sig.levels, inline=True, fmt=fmt, fontsize=9) 
    fig.colorbar(cm, ax=axes_to['WOA'], location='right', label=r'$\theta$ / $^\circ$C')
    
    cm = axes_to['cmip'].pcolormesh(to_plot['to']['cmip'][0], to_plot['to']['cmip'][1], to_plot['to']['cmip'][2], vmin=-2.5, vmax=2.5, cmap=cmocean.cm.delta)
    sig = axes_to['cmip'].contour(to_plot['sigma0']['cmip'][0], to_plot['sigma0']['cmip'][1], to_plot['sigma0']['cmip'][2], levels=levels, colors='k')
    axes_so['cmip'].clabel(sig, sig.levels, inline=True, fmt=fmt, fontsize=9) 
    fig.colorbar(cm, ax=axes_to['cmip'], location='right', label=r'$\theta$ / $^\circ$C')
    
    cm = axes_to['WOA-cmip'].pcolormesh(to_plot['to']['woa-cmip'][0], to_plot['to']['woa-cmip'][1], to_plot['to']['woa-cmip'][2], vmin=-1, vmax=1, cmap=cmocean.cm.balance)
    fig.colorbar(cm, ax=axes_to['WOA-cmip'], location='right', label=r'$\Delta\theta$ / $^\circ$C')
        
    fig.tight_layout()






    
    