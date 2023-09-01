#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 10:53:38 2022

@author: maxthomas_UO
"""

import parse_arguments
import lazy_load_cmip
from functools import partial
import xarray as xr
import iris
import iris.analysis as irisa
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import cartopy.crs        as ccrs
import cartopy.feature    as cfeature
import cftime
import dask
import glob
from matplotlib import colors
import cmip6_preprocessing.preprocessing       as cmip6
import gsw
from scipy.stats import ttest_ind
import warnings
import cmocean
warnings.filterwarnings("ignore")

dask.config.set({"array.slicing.split_large_chunks": True})

def qpm(ds):
    ds.attrs['source_id'] = 'HadGEM3-GC31-LL'
    ds = cmip6.rename_cmip6(ds)
    ds = cmip6.correct_coordinates(ds)
    ds = cmip6.correct_lon(ds)
    ds = cmip6.correct_units(ds)
    ds = cmip6.promote_empty_dims(ds)
    return ds

def process_maui(maui, model='oce'):  
    if model=='oce':
        maui = maui.rename({'nav_lat': 'lat', 'nav_lon': 'lon'})
        if 'deptht' in maui.coords:
            maui = maui.rename({'deptht': 'lev'})
        elif 'depthv' in maui.coords:
            maui = maui.rename({'depthv': 'lev'})
        elif 'depthu' in maui.coords:
            maui = maui.rename({'depthu': 'lev'})
        elif 'depthw' in maui.coords:
            maui = maui.rename({'depthw': 'lev'})
        maui = maui.swap_dims({'time_counter' : 'time_centered'}).rename({'time_centered':'time'}).drop_isel(y=[0,-1],x=[0,-1]).drop_vars('time_counter')
        # maui = maui.assign_coords(x=maui['x'].values-1)
        # maui = maui.assign_coords(y=maui['y'].values-1)
        maui.attrs['source_id'] = 'HadGEM3-GC31-LL'
        lon = maui['lon'].where(maui['lon'] > 0, maui['lon'] + 360)
        maui = maui.assign_coords(lon=lon)
        #maui = maui.assign_coords({'member_id': suite_pair(maui.suite_id.values)})
        #maui = maui.swap_dims({'suite_id': 'member_id'})
        # maui = lazy_load_cmip.combined_preprocessing.rename_cmip6(maui)
        # maui = lazy_load_cmip.combined_preprocessing.correct_coordinates(maui)
        # maui = lazy_load_cmip.combined_preprocessing.correct_lon(maui)
        # maui = lazy_load_cmip.combined_preprocessing.correct_units(maui)
        # maui = lazy_load_cmip.combined_preprocessing.promote_empty_dims(maui)  
    elif model=='atm':
        maui = maui.rename({'latitude': 'lat', 'longitude': 'lon'})
        maui = maui.assign_coords({'y': maui.lat, 'x': maui.lon})
        maui = maui.swap_dims({'lat': 'y', 'lon': 'x'})
    elif model=='ice':
        maui = maui.rename({'TLAT': 'lat', 'TLON': 'lon'})
        lon = maui['lon'].where(maui['lon'] > 0, maui['lon'] + 360)
        maui = maui.assign_coords(lon=lon)
        maui = maui.rename({'nj': 'y', 'ni': 'x'})
    return maui

def lon_select(ds, lon_low=-78, lon_high=-82):
    # lons = ds.lon.values
    # gt = lons > lon_low
    # lt = lons < lon_high
    # selector = np.multiply(gt,lt)
    # return ds.where(selector)

    # mean_lon = (lon_low + lon_high)/2
    # target_lon = (np.abs(ds.lon.values - mean_lon)).argmin()
    # target_lon = ds.lon.values.flatten()[target_lon]
    # return ds.where(ds.lon==target_lon, drop=True)
    return ds.where( (ds.lon>lon_low) & (ds.lon<lon_high), drop=True )

def lat_select(ds, lat_low=-90, lat_high=-60):
    return ds.where( (ds.lat>lat_low) & (ds.lat<lat_high), drop=True )

def depth_slice(ds, deep=1400, shallow=0):
    return ds.sel(lev=slice(shallow,deep))

def extract_depths(ds, depths):
    if 'lev' in ds.dims:
        if len(depths)==2:
            z0, z1, zx = depths[0], depths[1], ds.lev.values
            i0, i1 = np.argmin(np.abs(zx - z0)), np.argmin(np.abs(zx - z1))
            ds = ds.isel(lev=i0) - ds.isel(lev=i1)
        else:
            ds = ds.isel(lev=np.argmin(np.abs(ds.lev.values - depths[0])))
    return ds

def KtoC(ds):
    return ds - 273.15

def extract_season(ds, season='winter'):
    def is_season(month, season):
        if season == 'autumn':
            return (month >= 3) & (month <= 5)
        elif season == 'spring':
            return (month >= 9) & (month <= 11)
        elif season == 'winter':
            return (month >= 6) & (month <= 8)
        elif season == 'summer':
            return (month >= 12) | (month <= 2)
        elif season == 'year':
            return month > -1
        elif season == 'september':
            return month==9
        elif season == 'march':
            return month==3
    return ds.sel(time=is_season(ds['time.month'], season=season))

def make_mean(ds, dim='time'):
    return ds.mean(dim=dim)

def suite_pair(suite):
    if suite == 'ci501':
        return 'r1i1p1f3'
    elif suite == 'cm483':
        return 'r2i1p1f3'
    elif suite == 'cn043':
        return 'r3i1p1f3'
    elif suite == 'cn077':
        return 'r4i1p1f3'
    
def datestring2cftime(datestring='18500101'):
    return cftime.Datetime360Day(int(datestring[0:4]), int(datestring[4:6]), int(datestring[6:8]))


def read_netcdfs(files, dim='time', drop_vars=None, transform_func=None, density=False, model='oce'):  
    if model=='oce':
        if density:
            if 'so' in drop_vars:
                drop_vars.remove('so')
            if 'thetao' in drop_vars:
                drop_vars.remove('thetao')    
        def process_one_path(path):
            # use a context manager, to ensure the file gets closed after use
            with xr.open_dataset(path, drop_variables=drop_vars) as ds:
                # transform_func should do some sort of selection or
                # aggregation
                if transform_func is not None:
                    #if type(transform_func) is not list:
                    #    transform_func = [transform_func]
                    for tfunc in transform_func:
                        ds = tfunc(ds)
                # load all data from the transformed dataset, to ensure we can
                # use it after closing each original file
                ds.load()
                print(path)
                return ds
    #paths = sorted(am.glob(files))
    paths = sorted(files)
    datasets = [process_one_path(p) for p in paths]
    combined = xr.concat(datasets, dim)
    return combined

## Get directories
def get_path(get_this_dir, computer=None):
    """
    lookup and return path to a given directory using the table in:
        ../../documents/notebooks/paths/paths_maui

    Parameters
    ----------
    get_this_dir : string
        descriptive name of desired directory, currently:
            proj_home          /nesi/project/nesi00442/thoma97p/UO_postdoc/
            raw_data            /nesi/nobackup/nesi00442/thoma97p/cylc-run/
            processed_data                                            data/
            notes                                                documents/
            suite2atm                              share/data/History_Data/
            suite2oce       share/data/History_Data/NEMOhist/archive_ready/
            suite2ice       share/data/History_Data/CICEhist/archive_ready/
            variables                   documents/notebooks/variable_lists/
    computer : string
        name of computer running script
            maui (default), desktop, or laptop
    Returns
    -------
    string
        path to directory, ending with a slash.
    """
    if computer is None:
        paths = pd.read_table('paths', sep='\s+', index_col=0)
    else:
        paths = pd.read_table('../../documents/notebooks/paths/paths_'+computer, sep='\s+', index_col=0)
    return paths['path'][get_this_dir]

def list_maui_files(suite='ci501', model='oce', resolution='1m', grid='T'):
    data_at = get_path('raw_data')
    data_at = data_at + 'u-' + suite + '/'
    data_at = data_at + get_path('suite2' + model)
    if model=='oce':
        file_list = glob.glob(data_at+'/*1m*'+grid+'*')
    elif model=='atm':
        file_list = glob.glob(data_at+'/*p5*')
    elif model=='ice':
        file_list = glob.glob(data_at+'/*1m*')
    file_list = sorted(file_list)
    return file_list

def get_um_stashcode(variable_name):
    um_variables = pd.read_table(get_path('proj_home')+get_path('variables')+'um', sep='\s+', index_col=0)
    return um_variables.stash.loc[variable_name]
    
def vars2exclude(example_file, vars2keep):
    if type(vars2keep) is str:
        vars2keep = [vars2keep]
    all_vars = xr.open_dataset(example_file)
    drop_vars = []
    for vb in all_vars:
        if not all_vars[vb].name in vars2keep:
            drop_vars.append(all_vars[vb].name)
    return drop_vars

def make_cmap(ds, cmap='PiYG', abs_lims=None, diff_lims=None):
    if (abs_lims is None) and (not diff_lims==-999):
        # define via diff_lims
        vmin, vmax = diff_lims[0], diff_lims[1]
        diff_map = True
    elif (diff_lims is None) and (not abs_lims==-999):
        # define via abs lims
        vmin, vmax = abs_lims[0], abs_lims[1]
        diff_map=False
    else:
        # get from ds
        vmin, vmax = np.nanmin(ds), np.nanmax(ds)
        diff_map = False
        
    if cmap=='PiYG':
        cmap = 'divergent' # just for consistency with older plots while editing. Use continuous/divergent in future
        
    
    if (cmap=='divergent') | (diff_map): # if we want a divergent cmap centered on 0 or a difference map (i.e. temperature or difference in anything)
        return colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    else: # if we want a continuous map (i.e. salinity)
        return colors.TwoSlopeNorm(vmin=vmin, vcenter=(vmin+vmax)/2, vmax=vmax)

def sigmax(maui=None, so=None, thetao=None, i_sigma=0):
    if not maui is None:
        so = maui['so']
        thetao = maui['thetao']
    if i_sigma == 0:
        maui = gsw.sigma0(so, thetao)
    elif i_sigma == 1:
        maui = gsw.sigma1(so, thetao)
    elif i_sigma == 2:
        maui = gsw.sigma2(so, thetao)
    elif i_sigma == 3:
        maui = gsw.sigma3(so, thetao)
    elif i_sigma == 4:
        maui = gsw.sigma4(so, thetao)
    return maui        


def fmt(x):
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    return rf"{s} " if plt.rcParams["text.usetex"] else f"{s} "

def pcm(ds, norm, cmap_name, ax):
    try:
        cmap_out = xr.plot.pcolormesh(ds, 'lon', 'lat', transform=ccrs.PlateCarree(), norm=norm, cmap=cmap_name, ax=ax, add_colorbar=False, add_labels=False) 
    except: # weirdly u wind drops lon coord during subtraction. this is just a quick fix as lon and x are the same for um data
        print('Using x and x in pcolormesh rather than lat lon')
        cmap_out = xr.plot.pcolormesh(ds, 'x', 'y', transform=ccrs.PlateCarree(), norm=norm, cmap=cmap_name, ax=ax, add_colorbar=False, add_labels=False)        
    return cmap_out

def calculate_ekman(ds_u, ds_v, e1v, e2u):
    tau_y = ds_v.values
    tau_x = ds_u.values
    Omega = np.pi / 43082.
    f_u = 2 * Omega * np.sin(np.deg2rad(ds_u.lat.values))
    f_v = 2 * Omega * np.sin(np.deg2rad(ds_v.lat.values))
    rho_0 = 1026.
    
    tau_y = tau_y / f_v
    tau_x = tau_x / f_u

    d_tau_y = np.zeros_like(e2u)
    d_tau_x = np.zeros_like(e2u)
    d_y      = np.zeros_like(e2u)
    d_x      = np.zeros_like(e2u)
    
    # calculate terms for differentials
    d_tau_y[:,:-1] = tau_y[:,1:] - tau_y[:,:-1]
    d_x[:,:-1] = (e1v[:,1:] + e1v[:,:-1]) / 2
    d_tau_x[:-1,:] = tau_x[1:,:] - tau_x[:-1,:]
    d_y[:-1,:] = (e2u[1:,:] + e2u[:-1,:])/2
    
    d_tau_y[:,-1] = tau_y[:,0] - tau_y[:,-1]
    d_x[:,-1] = (e1v[:,0] + e1v[:,-1])/2
    
    ### dirty patch as these dont matter
    d_tau_x[-1,:] = d_tau_x[-2,:]
    d_y[-1,:] = d_y[-2,:]
    
    #d_tau_y = d_tau_y / f_v
    #d_tau_x = d_tau_x / f_u
    
    ekm = ( ( d_tau_y / d_x ) -  ( d_tau_x / d_y ) ) / rho_0
    
    return ekm  

def shelf_mask(ds):
    mm = xr.open_dataset('../../data/grid/mesh_mask_eORCA1_v2.2x.nc').isel(t=0).isel(y=slice(1,-1)).isel(x=slice(1,-1))    
    z_bot = (mm['e3t_0'] * mm['tmask']).sum(dim='z')
    z_bot = z_bot.assign_coords({'lon': ds.lon, 'lat': ds.lat})
    shelf = z_bot.where( (z_bot>0) & (z_bot<1000))
    return (shelf>0) & (mm.nav_lat<-60)

def sea_ice_contour(soicecov, levels=[0.15], threshold=0.15):
    
    sia_map = soicecov.copy()
        
    SIE = np.ma.where((sia_map[:,:])>=threshold,1,0)
    
    plt.figure()
    cs=plt.contour(SIE, levels=levels)
    Xs, Ys = [], []
    for item in cs.collections:
       for i in item.get_paths():
          v = i.vertices
          x = v[:, 0]
          y = v[:, 1]
          Xs.append(x.astype(int))
          Ys.append(y.astype(int))

    lats, lons = [], []
    for X,Y in zip(Xs,Ys):
        for ip in range(len(X)):
            lats.append(float(sia_map.lat[Y[ip],X[ip]]))
            lons.append(float(sia_map.lon[Y[ip],X[ip]]) )
    
    return lons, lats
    
def plot_sea_ice_contour(soicecov, ax, levels=[0.15], threshold=0.15, color='k'):
    lons, lats = sea_ice_contour(soicecov, levels=levels, threshold=threshold)
    ax.scatter(lons, lats, c=color, s=0.1, transform=ccrs.PlateCarree())
    
def get_sea_ice_area(args, cfdates=None, transform_functions=None):
    if cfdates is None:
        cfdates = make_cfdates(args['run_type'])[0]
    if transform_functions is None:
        # transform functions
        transform_functions = [partial(process_maui, model='oce'),
                               partial(lat_select, lat_low=args['lats'][0], lat_high=args['lats'][1]),
                               partial(extract_season, season=args['season'])]
    print('starting soicecov map')
    variable_id = 'soicecov'
    variable_info = pd.read_csv('../../documents/notebooks/variable_lists/variables_for_analysis.csv', index_col='name')
    label = variable_info.at[variable_id,'label1']
    table_id = model[0].upper() + 'mon'
    cmap_name  = variable_info.at[variable_id, 'cmap']
    cmip_units = 100. 
    
    
    cmip = {'hist': lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='historical', variable_id=variable_id) / cmip_units,
            'ssp': lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='ssp585', variable_id=variable_id) / cmip_units}
    
    for expid in ['hist', 'ssp']:
        cmip[expid] = cmip[expid].sel(time=slice(cfdates[expid][0], cfdates[expid][1])) 
        for tf in transform_functions[1:]:
            cmip[expid] = tf(cmip[expid])
                
    maui = {}
    for suite in args['suite']:
        file_list = list_maui_files(suite=suite, model='oce')[-maui_months:]
        maui[suite] = read_netcdfs(file_list, drop_vars=vars2exclude(file_list[0], variable_id), transform_func=transform_functions)[variable_id] 
    maui = xr.concat([maui[sn] for sn in args['suite']],
           pd.Index([sn for sn in args['suite']], name="suite_id"))
    maui = maui.sel(time=slice(cfdates['ssp'][0],cfdates['ssp'][1]))
    return cmip, maui
    
def make_cfdates(run_type='main'):#hist_times=('19941101','20141116'), ssp_times=('20791101','20991116')):
    if run_type=='main':
        hist_times=('19941101','20141116')
        ssp_times=('20791101','20991116')
        maui_months = 240
    elif run_type=='test': # just do three months for speed
        hist_times=('20131101','20141116')
        ssp_times=('20981101','20991116')
        maui_months = 12
    cfdates_out =  {'hist': [datestring2cftime(hist_times[0]), datestring2cftime(hist_times[1])],
                    'ssp':  [datestring2cftime(ssp_times[0]), datestring2cftime(ssp_times[1])]}
    return cfdates_out, maui_months

def regrid_mm():
    mm = xr.open_dataset('../../data/grid/mesh_mask_eORCA1_v2.2x.nc')
    mm = mm.rename({'nav_lon': 'longitude', 'nav_lat': 'latitude'})
    
    e3t_0 = mm['e3t_0'].assign_coords({'longitude': mm['longitude'],
                                       'latitude': mm['latitude']})
    tmask = mm['tmask'].assign_coords({'longitude': mm['longitude'],
                                       'latitude': mm['latitude']})
    for l in ['latitude','longitude']:
        e3t_0[l] = e3t_0[l].assign_attrs({'standard_name':l,
                                          'units':'degrees'})
        tmask[l] = tmask[l].assign_attrs({'standard_name':l,
                                          'units':'degrees'})
        
    e3t_0 = iris.analysis.cartography.project(
                e3t_0.to_iris(),
                ccrs.PlateCarree(), nx=360, ny=330)[0]
    
    tmask = iris.analysis.cartography.project(
                    tmask.to_iris(),
                    ccrs.PlateCarree(), nx=360, ny=330)[0]
    
    e3t_0 = xr.DataArray.from_iris(e3t_0).squeeze().rename({'dim_1': 'z'})
    tmask = xr.DataArray.from_iris(tmask).squeeze().rename({'dim_1': 'z'})
    
    
    return e3t_0, tmask 

def make_shelf_mask():
    e3t_0, tmask = regrid_mm()
    z_bot = (e3t_0 * tmask).where(e3t_0.latitude<-60).sum(dim='z')      
    shelf = z_bot.where( (z_bot>0) & (z_bot<1000))
    shelf2plot = np.where((shelf<=1000), 1, 0)
    return shelf, shelf2plot

def read_and_regrid(files, dim='time_counter', transform_func=None, density=False, model='oce', regrid=True, variable_id=None):  
    var_selector = variable_id
    def open_ds(path2load, var_selector):                  
        opened_ds = iris.analysis.cartography.project(
                    iris.load_cube(path2load, var_selector),
                    ccrs.PlateCarree(), nx=360, ny=330)
        return xr.DataArray.from_iris(opened_ds[0])
    
        
    def process_one_path(path, var_selector=var_selector):
        # use a context manager, to ensure the file gets closed after use
        with open_ds(path, var_selector) as ds:
            # transform_func should do some sort of selection or
            # aggregation
            
            if transform_func is not None:
                #if type(transform_func) is not list:
                #    transform_func = [transform_func]
                for tfunc in transform_func:
                    ds = tfunc(ds)
            # load all data from the transformed dataset, to ensure we can
            # use it after closing each original file
            ds.load()
            ds = ds.rename({'longitude': 'lon', 'latitude': 'lat'})
            print(path)
            return ds
        
    #paths = sorted(am.glob(files))
    paths = sorted(files)
    datasets = [process_one_path(p) for p in paths]
    combined = xr.concat(datasets, dim)
    return combined

def regridded_seaice_contour(files, threshold=0.15):
    sic = read_and_regrid(files, variable_id='soicecov')
    sic = sic.mean(dim='time_counter')    
    SIE = np.where((sic>threshold), 1, 0)
    return sic, SIE

def get_plot_colours(colour_set='main'):
    colours = {'cmaps': {}, 'lines': {}}
    if colour_set=='main':
        colours['cmaps']['continuous'] = cmocean.cm.haline
        colours['cmaps']['divergent'] = cmocean.cm.delta
        colours['cmaps']['difference'] = cmocean.cm.balance
        colours['lines']['coast'] = 'k'
        colours['lines']['shelf'] = 'k'
        colours['lines']['hist'] = 'c'
        colours['lines']['ssp'] = 'grey'
        colours['lines']['maui'] = 'k'
        return colours

def make_map_axes(panels=['hist','ssp_hist','ssp','maui_hist','maui','maui_ssp']):
    projection = ccrs.SouthPolarStereo()
    #panels = ['hist','ssp','maui','ssp_hist','maui_hist','maui_ssp']
    fig = plt.figure(dpi=300, figsize=(12,18))
    gs  = fig.add_gridspec(3,2)
    gs.update(wspace=0.1, hspace=0.1)
    axes = {}
    for ip, panel in enumerate(panels):
        axes[panel] = fig.add_subplot(gs[ip], projection=projection)
        axes[panel].set_facecolor((0.8,0.8,0.8))
        axes[panel].gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
        #axes[panel].add_feature(cfeature.COASTLINE) 
        axes[panel].set_extent([-180, 180, -90, -65], ccrs.PlateCarree()) 
    return fig, axes

if __name__ == '__main__':
    args = parse_arguments.inputs2arguments()
    variable_id = args['variable'][0]
    variable_info = pd.read_csv('../../documents/notebooks/variable_lists/variables_for_analysis.csv', index_col='name')
    model = variable_info.at[variable_id, 'model']
    grid  = variable_info.at[variable_id, 'grid']
    label = variable_info.at[variable_id,'label1']
    table_id = model[0].upper() + 'mon'
    cmap_name  = variable_info.at[variable_id, 'cmap']
    colours = get_plot_colours()
    cfdates, maui_months = make_cfdates(args['run_type'])
    # cfdates = {'hist': [datestring2cftime('19941101'), datestring2cftime('20141116')],
    #            'ssp':  [datestring2cftime('20791101'), datestring2cftime('20991116')]}
    #cfdates = {'hist': [datestring2cftime('19950101'), datestring2cftime('19951116')],
    #           'ssp':  [datestring2cftime('20990101'), datestring2cftime('20991116')]}
    #maui_months = 240
    lon_pairs = [  [300,360], [240,300], [180,240], [120,180], [60,120], [0,60], [170,180], [210,290] ]
    levels_fontsize = 9
    section_linewidth = 2
    # number of monthly means to load in for maui. 240 gives the last 20 years
    #cmip_scaler = variable_info.at[variable_id, 'cmip_scaler'] 
    if variable_id=='soicecov':
        cmip_units = 100.
    else:
        cmip_units = 1.
        
    meshmask = xr.open_dataset('../../data/grid/mesh_mask_eORCA1_v2.2x.nc')
    land_mask = xr.DataArray.from_iris(iris.load('../../data/grid/qrparm.mask')[0])    
    
    
    ###########################################################################
    # lon, lat map
    # end of historical, end of maui, maui - hist, maui - cmip
    ###########################################################################
    if 'map' in args['analysis']:
        
        # transform functions
        transform_functions = [partial(process_maui, model=model),
                               partial(lat_select, lat_low=args['lats'][0], lat_high=args['lats'][1]),
                               partial(extract_depths, depths=args['depths']),
                               partial(extract_season, season=args['season'])]#,
        
        def cmip2plot(cmip_data, variable_id=variable_id, tfs=transform_functions.copy(), cfdates=cfdates, cmip_units=cmip_units):
            # making data
            data = {'hist': cmip_data['historical'].copy(),
                    'ssp':  cmip_data['ssp585'].copy()}
            to_plot = {}
            if variable_id=='soicecov':
                tfs[-1] = partial(extract_season, season='winter')
            for eid in ['hist', 'ssp']:
                for tf in transform_functions[1:]:
                    data[eid] = tf(data[eid])
                data[eid] = data[eid].sel(time=slice(cfdates[eid][0],cfdates[eid][1])).mean(dim='time') / cmip_units
                if variable_id == '1p5m_air_temperature':
                    data[eid] = data[eid] - 273.15
                to_plot[eid] = data[eid].mean(dim='member_id').compute()
                print('Done ' + eid)
            return to_plot
        
        def maui2plot(args, variable_id=variable_id, tfs=transform_functions.copy(), grid=grid, model=model, maui_months=maui_months):
            if variable_id=='soicecov':
                tfs[-1] = partial(extract_season, season='winter')
            print('starting maui load')
            maui = {}
            for suite in args['suite']:
                if model=='atm':
                    print(list_maui_files(suite=suite, model=model)[-maui_months:])
                    maui[suite] = xr.DataArray.from_iris(iris.load(list_maui_files(suite=suite, model=model)[-maui_months:], \
                                    iris.AttributeConstraint(STASH=get_um_stashcode(variable_id)))[0])
                    for tf in tfs[:2]:
                        maui[suite] = tf(maui[suite])  
                    if variable_id == '1p5m_air_temperature':
                        maui[suite] = maui[suite] - 273.15
                else:
                    maui[suite] = read_netcdfs(list_maui_files(suite=suite, model=model, grid=grid)[-maui_months:], drop_vars=vars2exclude(list_maui_files(suite=suite, model=model, grid=grid)[0], variable_id), transform_func=transform_functions)[variable_id] 
            maui = xr.concat([maui[sn] for sn in args['suite']],
                   pd.Index([sn for sn in args['suite']], name="suite_id"))
            maui = maui.sel(time=slice(cfdates['ssp'][0],cfdates['ssp'][1])).mean(dim='time')
            return maui
        
        def sic_contour(sic_data, ax, levels=[0.15], colour='k', linestyles='solid'):
            sic = sic_data.where(sic_data<0.15)*0
            sic_mask = sic.isin(0)
            sic_mask.plot.contour(x='lon', y='lat', transform=ccrs.PlateCarree(), colors=colour, ax=ax)
            #xr.plot.contour(sic_data, 'lon', 'lat', transform=ccrs.PlateCarree(), levels=levels, ax=ax, colors=colour, linestyles=linestyles) 
                    
        cmip = {'historical': lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='historical', variable_id=variable_id),
                'ssp585': lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='ssp585', variable_id=variable_id)}
        
        #cmip_sic = {'historical': lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='historical', variable_id='soicecov'),
        #            'ssp585': lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='ssp585', variable_id='soicecov')}
        
        #land_mask_map = transform_functions[1](land_mask)
        land_mask_map = land_mask.where( (land_mask.latitude>args['lats'][0]) & (land_mask.latitude<args['lats'][1]), drop=True )
        land_mask_map = np.where(land_mask_map > 0.5,True,False) 
        
            
        to_plot = cmip2plot(cmip)
        
        maui = maui2plot(args)
  
        
        # map axes
        #projection = ccrs.NearsidePerspective(central_latitude=-90, satellite_height=3.6e6)
        # projection = ccrs.SouthPolarStereo()
        # #panels = ['hist','ssp','maui','ssp_hist','maui_hist','maui_ssp']
        # panels = ['hist','ssp_hist','ssp','maui_hist','maui','maui_ssp']
        # fig = plt.figure(dpi=300, figsize=(12,18))
        # gs  = fig.add_gridspec(3,2)
        # gs.update(wspace=0.1, hspace=0.1)
        # axes = {}
        # for ip, panel in enumerate(panels):
        #     axes[panel] = fig.add_subplot(gs[ip], projection=projection)
        #     axes[panel].set_facecolor((0.8,0.8,0.8))
        #     axes[panel].gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
        #     #axes[panel].add_feature(cfeature.COASTLINE) 
        #     axes[panel].set_extent([-180, 180, -90, -65], ccrs.PlateCarree()) 
        
        panels=['hist','ssp_hist','ssp','maui_hist','maui','maui_ssp']
        fig, axes = make_map_axes(panels=['hist','ssp_hist','ssp','maui_hist','maui','maui_ssp'])
 
        to_plot['maui'] = maui.mean(dim=('suite_id')).compute()
        #to_plot_sic['maui'] = maui_sic.mean(dim=('suite_id')).compute()
        print(to_plot['hist'])
        print(to_plot['ssp'])
        print(to_plot['maui'])        
        
        to_plot['ssp_hist'] = to_plot['ssp'] - to_plot['hist']
        to_plot['maui_ssp'] = to_plot['maui'] - to_plot['ssp']
        to_plot['maui_hist'] = to_plot['maui'] - to_plot['hist']
        
        if args['land_mask']:
            for panel in panels:
                to_plot[panel].values = np.ma.MaskedArray(to_plot[panel], mask=land_mask_map)
        
        print(cmap_name)
        print(colours['cmaps'][cmap_name])
        
        if not 'lon' in to_plot['maui_ssp'].coords:
            to_plot['maui_ssp'] = to_plot['maui_ssp'].assign_coords({'lon': to_plot['hist'].lon}) # for some reason uo and vo are dropping lon coord during subtraction. This is just an easy fix
        cmaps = [make_cmap(to_plot['maui'], cmap_name, abs_lims=args['abs_lims']),
                 make_cmap(to_plot['maui'], cmap_name, abs_lims=args['abs_lims']),
                 make_cmap(to_plot['maui'], cmap_name, abs_lims=args['abs_lims']),
                 make_cmap(to_plot['ssp_hist'], 'difference', diff_lims=args['diff_lims']),
                 make_cmap(to_plot['maui_ssp'], 'difference', diff_lims=args['diff_lims']),#]
                 make_cmap(to_plot['maui_hist'], 'difference', diff_lims=args['diff_lims'])]

        cmabs = pcm(to_plot['hist'], cmaps[0], colours['cmaps'][cmap_name], axes['hist'])
        cmabs = pcm(to_plot['ssp'], cmaps[1], colours['cmaps'][cmap_name], axes['ssp'])
        cmabs = pcm(to_plot['maui'], cmaps[2], colours['cmaps'][cmap_name], axes['maui'])
        cmhdif = pcm(to_plot['ssp_hist'], cmaps[3], colours['cmaps']['difference'], axes['ssp_hist'])
        cmpdif = pcm(to_plot['maui_hist'], cmaps[4], colours['cmaps']['difference'], axes['maui_hist'])
        cmpdif = pcm(to_plot['maui_ssp'], cmaps[5], colours['cmaps']['difference'], axes['maui_ssp'])
        
        axes['hist'].plot([252,252],[-76,-68],'m-',transform=ccrs.PlateCarree())
        axes['hist'].plot([177,177],[-78.2,-70],'m-',transform=ccrs.PlateCarree())
        axes['hist'].plot([71,71],[-69,-65],'m-',transform=ccrs.PlateCarree())

        axes['hist'].plot([210,210], [-75,-60], 'm:',transform=ccrs.PlateCarree())
        axes['hist'].plot([290,290], [-75,-60], 'm:',transform=ccrs.PlateCarree())
        axes['hist'].plot([170,170], [-75,-60], 'm:',transform=ccrs.PlateCarree())
        axes['hist'].plot([180,180], [-75,-60], 'm:',transform=ccrs.PlateCarree())
        axes['hist'].plot([120,120], [-75,-60], 'm:',transform=ccrs.PlateCarree())
        axes['hist'].plot([60,60], [-75,-60], 'm:',transform=ccrs.PlateCarree())
        #sic_contour(to_plot_sic['hist'], axes['hist'], colour='c')
        #sic_contour(to_plot_sic['ssp'], axes['ssp'])
        #sic_contour(to_plot_sic['maui'], axes['maui'], linestyles='dashed')
        # ssp_hist
        #sic_contour(to_plot_sic['hist'], axes['ssp_hist'], colour='c')
        #sic_contour(to_plot_sic['ssp'], axes['ssp_hist'])
        # maui_hist
        #sic_contour(to_plot_sic['hist'], axes['maui_hist'], colour='c')
        #sic_contour(to_plot_sic['maui'], axes['maui_hist'], linestyles='dashed')
        # maui_ssp
        #sic_contour(to_plot_sic['maui'], axes['maui_ssp'], linestyles='dashed')
        #sic_contour(to_plot_sic['ssp'], axes['maui_ssp'])
        
        
        # fig.colorbar(cmabs, ax=axes['hist'], shrink=0.5, location='right', label=variable_info.at[variable_id,'label1'])
        # fig.colorbar(cmabs, ax=axes['maui'], shrink=0.5, location='right', label=variable_info.at[variable_id,'label1'])
        # fig.colorbar(cmabs, ax=axes['ssp'], shrink=0.5, location='right', label=variable_info.at[variable_id,'label1'])
        # fig.colorbar(cmhdif, ax=axes['ssp_hist'], shrink=0.5, location='right',label='$\Delta$'+variable_info.at[variable_id,'label1'])
        # fig.colorbar(cmpdif, ax=axes['maui_hist'], shrink=0.5, location='right',label='$\Delta$'+variable_info.at[variable_id,'label1'])
        # fig.colorbar(cmpdif, ax=axes['maui_ssp'], shrink=0.5, location='right',label='$\Delta$'+variable_info.at[variable_id,'label1'])
        
        if args['plot_sia_contour']:
            cmip_sia, maui_sia = get_sea_ice_area(args)
            for expid in ['hist', 'ssp']:
                cmip_sia[expid] = cmip_sia[expid].mean(dim=('member_id','time')).compute()
            maui_sia = maui_sia.mean(dim=('suite_id','time')).compute()
            
            to_plot_sia = {'hist': cmip_sia['hist'],
                           'ssp': cmip_sia['ssp'],
                           'maui': maui_sia}
            to_plot_sia['ssp_hist'] = to_plot_sia['ssp'] - to_plot_sia['hist']
            to_plot_sia['maui_hist'] = to_plot_sia['maui'] - to_plot_sia['hist']
            to_plot_sia['maui_ssp'] = to_plot_sia['maui'] - to_plot_sia['ssp']
            
            plot_sea_ice_contour(to_plot_sia['hist'], axes['hist'], color=colours['lines']['hist'])
            plot_sea_ice_contour(to_plot_sia['hist'], axes['ssp_hist'], color=colours['lines']['hist'])
            plot_sea_ice_contour(to_plot_sia['hist'], axes['maui_hist'], color=colours['lines']['hist'])
            
            plot_sea_ice_contour(to_plot_sia['ssp'], axes['ssp'], color=colours['lines']['ssp'])
            plot_sea_ice_contour(to_plot_sia['ssp'], axes['ssp_hist'], color=colours['lines']['ssp'])
            plot_sea_ice_contour(to_plot_sia['ssp'], axes['maui_ssp'], color=colours['lines']['ssp'])
            
            plot_sea_ice_contour(to_plot_sia['maui'], axes['maui'], color=colours['lines']['maui'])
            plot_sea_ice_contour(to_plot_sia['maui'], axes['maui_hist'], color=colours['lines']['maui'])
            plot_sea_ice_contour(to_plot_sia['maui'], axes['maui_ssp'], color=colours['lines']['maui'])
            
        if args['plot_shelf']:
            shelf, shelf2plot = make_shelf_mask()
            for ax in panels:
                axes[ax].contour(shelf.longitude, shelf.latitude, shelf2plot, transform=ccrs.PlateCarree(), levels=[1], linestyles=['-'], colors=[colours['lines']['shelf']])
                
        if args['plot_0_wind']: # plot u wind 0 line. This is not an ideal way to do this but I was getting wrapping problems with xarrays
            hist_cube = to_plot['hist'].isel(x=slice(1,-1)).to_iris()
            ssp_cube = to_plot['ssp'].isel(x=slice(1,-1)).to_iris()
            maui_cube = to_plot['maui'].isel(x=slice(1,-1)).to_iris()

            
            print(hist_cube)
            print(ssp_cube)
            print(maui_cube)

            axes['hist'].contour(hist_cube.coord('lon').points, hist_cube.coord('lat').points, np.transpose(hist_cube.data),transform=ccrs.PlateCarree(),levels=[0],colors=colours['lines']['hist'])
            axes['ssp'].contour(ssp_cube.coord('lon').points, ssp_cube.coord('lat').points, np.transpose(ssp_cube.data),transform=ccrs.PlateCarree(),levels=[0],colors=colours['lines']['ssp'])
            axes['maui'].contour(ssp_cube.coord('lon').points, ssp_cube.coord('lat').points, np.transpose(maui_cube.data),transform=ccrs.PlateCarree(),levels=[0],colors=colours['lines']['maui'],linestyles='dashed')

            axes['ssp_hist'].contour(hist_cube.coord('lon').points, hist_cube.coord('lat').points, np.transpose(hist_cube.data),transform=ccrs.PlateCarree(),levels=[0],colors=colours['lines']['hist'])
            axes['ssp_hist'].contour(ssp_cube.coord('lon').points, ssp_cube.coord('lat').points, np.transpose(ssp_cube.data),transform=ccrs.PlateCarree(),levels=[0],colors=colours['lines']['ssp'])

            axes['maui_hist'].contour(hist_cube.coord('lon').points, hist_cube.coord('lat').points, np.transpose(hist_cube.data),transform=ccrs.PlateCarree(),levels=[0],colors=colours['lines']['hist'])
            axes['maui_hist'].contour(ssp_cube.coord('lon').points, ssp_cube.coord('lat').points, np.transpose(maui_cube.data),transform=ccrs.PlateCarree(),levels=[0],colors=colours['lines']['maui'],linestyles='dashed')

            axes['maui_ssp'].contour(ssp_cube.coord('lon').points, ssp_cube.coord('lat').points, np.transpose(ssp_cube.data),transform=ccrs.PlateCarree(),levels=[0],colors=colours['lines']['ssp'])
            axes['maui_ssp'].contour(ssp_cube.coord('lon').points, ssp_cube.coord('lat').points, np.transpose(maui_cube.data),transform=ccrs.PlateCarree(),levels=[0],colors=colours['lines']['maui'],linestyles='dashed')
            #xr.plot.contour(to_plot['hist'], 'lon', 'lat', transform=ccrs.PlateCarree(), color=colours['lines']['hist'], ax=axes['hist'], levels=[0]) 
            #xr.plot.contour(to_plot['ssp'], 'lon', 'lat', transform=ccrs.PlateCarree(), color=colours['lines']['ssp'], ax=axes['ssp'], levels=[0]) 
            #xr.plot.contour(to_plot['maui'], 'lon', 'lat', transform=ccrs.PlateCarree(), color=colours['lines']['maui'], ax=axes['maui'], levels=[0]) 
           
            #xr.plot.contour(to_plot['hist'], 'lon', 'lat', transform=ccrs.PlateCarree(), color=colours['lines']['hist'], ax=axes['ssp_hist'], levels=[0]) 
            #xr.plot.contour(to_plot['ssp'], 'lon', 'lat', transform=ccrs.PlateCarree(), color=colours['lines']['ssp'], ax=axes['ssp_hist'], levels=[0]) 
            
            #xr.plot.contour(to_plot['hist'], 'lon', 'lat', transform=ccrs.PlateCarree(), color=colours['lines']['hist'], ax=axes['maui_hist'], levels=[0]) 
            #xr.plot.contour(to_plot['maui'], 'lon', 'lat', transform=ccrs.PlateCarree(), color=colours['lines']['maui'], ax=axes['maui_hist'], levels=[0]) 
            
            #xr.plot.contour(to_plot['maui'], 'lon', 'lat', transform=ccrs.PlateCarree(), color=colours['lines']['maui'], ax=axes['maui_ssp'], levels=[0]) 
            #xr.plot.contour(to_plot['ssp'], 'lon', 'lat', transform=ccrs.PlateCarree(), color=colours['lines']['ssp'], ax=axes['maui_ssp'], levels=[0]) 
            

                
        if args['save_data']:
            data_file_name = 'ISMIPfw-data_for_map-'+variable_id+'-'
            for expid in panels:
                to_plot[expid].to_netcdf(data_file_name+expid+'.nc')
        
        fig.tight_layout()
        fig.savefig('ISMIPfw-map-'+variable_id+'.png')
       
        fig.colorbar(cmabs, ax=axes['maui'], shrink=0.75, location='bottom', label=variable_info.at[variable_id,'label1'])
        fig.colorbar(cmpdif, ax=axes['maui_ssp'], shrink=0.75, location='bottom',label='$\Delta$'+variable_info.at[variable_id,'label1']) 
        fig.savefig('ISMIPfw-map-'+variable_id+'-colorbars.png') 

    if 'section' in args['analysis']:
        
        
        if 'mackie' in args['suite']:
            Mackie = True
            args['suite'] = args['suite'][1:] # so if you want to plot ay187 (FW), do --suites=mackie,ay187  can be: ['ay187', 'az576', 'bb819']
            cfdates = {'hist': [datestring2cftime('19291101'), datestring2cftime('19491116')],
                       'ssp':  [datestring2cftime('19291101'), datestring2cftime('19491116')]}
        else:
            Mackie = False 
        # transform functions
        transform_functions = [partial(process_maui, model='oce'),
                               partial(lat_select, lat_low=args['lats'][0], lat_high=args['lats'][1]),
                               partial(lon_select, lon_low=args['lons'][0], lon_high=args['lons'][1]),
                               partial(depth_slice, deep=args['depths'][1], shallow=args['depths'][0]),
                               partial(extract_season, season=args['season'])]
        
        if not Mackie:
            cmip_so = {'historical': lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='historical', variable_id='so'),
                       'ssp585': lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='ssp585', variable_id='so')}
            cmip_thetao = {'historical': lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='historical', variable_id='thetao'),
                       'ssp585': lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='ssp585', variable_id='thetao')}
        else: # this is a bit of a fudge to get mackie data plotting with the same code. the historical control becomes pi, and the ssp585 control becomes 1%co2
            print('loading PI and 1pctCO2')
            cmip_so = {'historical': lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='piControl', variable_id='so'),
                       'ssp585': lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='1pctCO2', variable_id='so')}
            cmip_thetao = {'historical': lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='piControl', variable_id='thetao'),
                       'ssp585': lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='1pctCO2', variable_id='thetao')}
            
        
        # making data
        # cmip
        # we need so and thetao to calculate densities
        # data = {'hist': cmip['historical'].copy(),
        #         'ssp':  cmip['ssp585'].copy()}
        data_so = {'hist': cmip_so['historical'].copy(),
                'ssp':  cmip_so['ssp585'].copy()}
        data_thetao = {'hist': cmip_thetao['historical'].copy(),
                'ssp':  cmip_thetao['ssp585'].copy()}
        #to_plot = {}
        for eid in ['hist', 'ssp']:
            for tf in transform_functions[1:]:
                #data[eid] = tf(data[eid])
                data_so[eid] = tf(data_so[eid])
                data_thetao[eid] = tf(data_thetao[eid])
            #data[eid] = data[eid].sel(time=slice(cfdates[eid][0],cfdates[eid][1])).mean(dim=('member_id','time')) / cmip_units
            data_so[eid] = data_so[eid].sel(time=slice(cfdates[eid][0],cfdates[eid][1])).mean(dim=('member_id','time')) 
            data_thetao[eid] = data_thetao[eid].sel(time=slice(cfdates[eid][0],cfdates[eid][1])).mean(dim=('member_id','time'))
            #to_plot[eid] = data[eid].compute()
            print('Done ' + eid)
        
        print('starting maui load')
        #maui = {}
        maui_so = {}
        maui_thetao = {}
        for suite in args['suite']:
            if not Mackie:
                #maui[suite] = read_netcdfs(list_maui_files(suite=suite)[-maui_months:], drop_vars=vars2exclude(list_maui_files()[0], variable_id), transform_func=transform_functions)
                maui_thetao[suite] = read_netcdfs(list_maui_files(suite=suite)[-maui_months:], drop_vars=vars2exclude(list_maui_files()[0], 'thetao'), transform_func=transform_functions)
                maui_so[suite] = read_netcdfs(list_maui_files(suite=suite)[-maui_months:], drop_vars=vars2exclude(list_maui_files()[0], 'so'), transform_func=transform_functions)
            else:
                file_list = sorted(glob.glob('../../data/'+suite+'/nemo_'+suite+'o_1m*'))[-maui_months:]
                maui_thetao[suite] = read_netcdfs(file_list, drop_vars=vars2exclude(list_maui_files()[0], 'thetao'), transform_func=transform_functions)
                maui_so[suite] = read_netcdfs(file_list, drop_vars=vars2exclude(list_maui_files()[0], 'so'), transform_func=transform_functions)
        
        maui_so = xr.concat([maui_so[sn] for sn in args['suite']],
               pd.Index([sn for sn in args['suite']], name="suite_id"))
        maui_so = maui_so['so'] 
        
        maui_thetao = xr.concat([maui_thetao[sn] for sn in args['suite']],
               pd.Index([sn for sn in args['suite']], name="suite_id"))
        maui_thetao = maui_thetao['thetao'] 
        
        # sigma0 = {}
        # sigma0['hist'] = sigmax(so=data_so['hist'], thetao=data_thetao['hist']).mean(dim='x')
        # sigma0['ssp'] = sigmax(so=data_so['ssp'], thetao=data_thetao['ssp']).mean(dim='x')
        # sigma0['maui'] = sigmax(so=maui_so.mean(dim='suite_id'), thetao=maui_thetao).mean(dim='x')
        
        maui_so = maui_so.sel(time=slice(cfdates[eid][0],cfdates[eid][1])).mean(dim=('suite_id','time')).compute()
        maui_thetao = maui_thetao.sel(time=slice(cfdates[eid][0],cfdates[eid][1])).mean(dim=('suite_id','time')).compute()
        

        data_so['hist']['lat'] = data_so['hist']['lat'].mean(dim='x')
        data_so['ssp']['lat'] = data_so['ssp']['lat'].mean(dim='x')
        maui_so['lat'] = maui_so['lat'].mean(dim='x')
        
        data_thetao['hist']['lat'] = data_thetao['hist']['lat'].mean(dim='x')
        data_thetao['ssp']['lat'] = data_thetao['ssp']['lat'].mean(dim='x')
        maui_thetao['lat'] = maui_thetao['lat'].mean(dim='x')
                
        to_plot = {}
        to_plot['so'] = {
                         'hist': data_so['hist'].mean(dim='x').compute(),
                         'ssp': data_so['ssp'].mean(dim='x').compute(),
                         'maui': maui_so.mean(dim='x').compute()
                        }
        to_plot['so']['ssp_hist'] = to_plot['so']['ssp'] - to_plot['so']['hist']
        to_plot['so']['maui_ssp'] = to_plot['so']['maui'] - to_plot['so']['ssp']
        to_plot['so']['maui_hist'] = to_plot['so']['maui'] - to_plot['so']['hist']
        
        to_plot['thetao'] = {
                         'hist': data_thetao['hist'].mean(dim='x').compute(),
                         'ssp': data_thetao['ssp'].mean(dim='x').compute(),
                         'maui': maui_thetao.mean(dim='x').compute()
                        }
        to_plot['thetao']['ssp_hist'] = to_plot['thetao']['ssp'] - to_plot['thetao']['hist']
        to_plot['thetao']['maui_ssp'] = to_plot['thetao']['maui'] - to_plot['thetao']['ssp']
        to_plot['thetao']['maui_hist'] = to_plot['thetao']['maui'] - to_plot['thetao']['hist']
        
        to_plot['sigma0'] = {
                             'hist': sigmax(so=to_plot['so']['hist'], thetao=to_plot['thetao']['hist']),
                             'ssp': sigmax(so=to_plot['so']['ssp'], thetao=to_plot['thetao']['ssp']),
                             'maui': sigmax(so=to_plot['so']['maui'], thetao=to_plot['thetao']['maui'])
                             }
        args['abs_lims'] = [33,35]
        args['diff_lims'] = [-1.5,1.5]
        cmap_name = 'continuous'       
        cmaps = [make_cmap(to_plot['so']['maui'], colours['cmaps'][cmap_name], abs_lims=args['abs_lims']),
                 make_cmap(to_plot['so']['maui'], colours['cmaps'][cmap_name], abs_lims=args['abs_lims']),
                 make_cmap(to_plot['so']['maui'], colours['cmaps'][cmap_name], abs_lims=args['abs_lims']),
                 make_cmap(to_plot['so']['ssp_hist'], colours['cmaps']['difference'], diff_lims=args['diff_lims']),
                 make_cmap(to_plot['so']['maui_ssp'], colours['cmaps']['difference'], diff_lims=args['diff_lims']),
                 make_cmap(to_plot['so']['maui_hist'], colours['cmaps']['difference'], diff_lims=args['diff_lims'])]
    
        cols = [cmap_name, cmap_name, cmap_name, 'difference', 'difference', 'difference']
        
        vmm = [args['abs_lims'], args['abs_lims'], args['abs_lims'], args['diff_lims'], args['diff_lims'], args['diff_lims']]
    
        levels = [26.4, 26.6, 26.8, 27.0, 27.2, 27.4, 27.6, 27.8]

        panels = ['hist', 'ssp', 'maui', 'ssp_hist', 'maui_hist', 'maui_ssp']
        Xs, Ys = np.meshgrid(to_plot['so']['maui'].lat, to_plot['so']['maui'].lev)
        axes = {}
        fig = plt.figure(dpi=300, figsize=(5,15))
        gs = fig.add_gridspec(len(panels),1)
        #ax = fig.add_subplot(gs[0])
        #ax.set_ylabel('test')
        for igs, pn in enumerate(panels):
            axes[pn] = fig.add_subplot(gs[igs])
            axes[pn].set_ylabel('depth, $z$ / m', fontsize=14)
            if igs==len(panels)-1:
                axes[pn].set_xlabel('latitude, $\phi$ / $^\mathrm{o}$N', fontsize=14)
            else:
                axes[pn].set_xticklabels([])
            axes[pn].grid('major')

            axes[pn].invert_yaxis()
            axes[pn].set_facecolor((0.8,0.8,0.8))
            axes[pn].set_xlim((args['lats'][0], args['lats'][1]))            

            print('start panel '+pn)

            cmh = axes[pn].pcolormesh(Xs,Ys,to_plot['so'][pn], norm=cmaps[igs], cmap=colours['cmaps'][cols[igs]])
            if cols[igs]=='difference':
                fig.colorbar(cmh, ax=axes[pn], label=r'$\Delta S$ / g.kg$^{-1}$', shrink=0.75, location='right')
            else:
                fig.colorbar(cmh, ax=axes[pn], label=r'$S$ / g.kg$^{-1}$', shrink=0.75, location='right')
              
                
            
            print('doing sigma0')
            if pn=='hist':
                sig = axes['hist'].contour(Xs, Ys, to_plot['sigma0']['hist'], colors=colours['lines']['hist'], linewidths=section_linewidth, levels=levels)
                axes['hist'].clabel(sig, sig.levels, inline=True, fmt=fmt, fontsize=levels_fontsize)
            elif pn=='ssp':
                sig = axes['ssp'].contour(Xs, Ys, to_plot['sigma0']['ssp'], colors=colours['lines']['ssp'], linewidths=section_linewidth, levels=levels)
                axes['ssp'].clabel(sig, sig.levels, inline=True, fmt=fmt, fontsize=levels_fontsize)
            elif pn=='maui':
                sig = axes['maui'].contour(Xs, Ys, to_plot['sigma0']['maui'], colors=colours['lines']['maui'], linewidths=section_linewidth, levels=levels)
                axes['maui'].clabel(sig, sig.levels, inline=True, fmt=fmt, fontsize=levels_fontsize)
            elif pn=='ssp_hist':
                sig = axes['ssp_hist'].contour(Xs, Ys, to_plot['sigma0']['ssp'], colors=colours['lines']['ssp'], linewidths=section_linewidth, levels=levels)
                axes['ssp_hist'].clabel(sig, sig.levels, inline=True, fmt=fmt, fontsize=levels_fontsize)
                sig = axes['ssp_hist'].contour(Xs, Ys, to_plot['sigma0']['hist'], colors=colours['lines']['hist'], linewidths=section_linewidth, levels=levels)
                axes['ssp_hist'].clabel(sig, sig.levels, inline=True, fmt=fmt, fontsize=levels_fontsize)
            elif pn=='maui_hist':
                sig = axes['maui_hist'].contour(Xs, Ys, to_plot['sigma0']['maui'], colors=colours['lines']['maui'], linewidths=section_linewidth, levels=levels)
                axes['maui_hist'].clabel(sig, sig.levels, inline=True, fmt=fmt, fontsize=levels_fontsize)
                sig = axes['maui_hist'].contour(Xs, Ys, to_plot['sigma0']['hist'], colors=colours['lines']['hist'], linewidths=section_linewidth, levels=levels)
                axes['maui_hist'].clabel(sig, sig.levels, inline=True, fmt=fmt, fontsize=levels_fontsize)
            elif pn=='maui_ssp':
                sig = axes['maui_ssp'].contour(Xs, Ys, to_plot['sigma0']['ssp'], colors=colours['lines']['ssp'], linewidths=section_linewidth, levels=levels)
                axes['maui_ssp'].clabel(sig, sig.levels, inline=True, fmt=fmt, fontsize=levels_fontsize)
                sig = axes['maui_ssp'].contour(Xs, Ys, to_plot['sigma0']['maui'], colors=colours['lines']['maui'], linewidths=section_linewidth, levels=levels)
                axes['maui_ssp'].clabel(sig, sig.levels, inline=True, fmt=fmt, fontsize=levels_fontsize)
               
            
        fig.tight_layout()
        fig.savefig('ISMIPfw-section-so.png')

        args['abs_lims'] = [-2.5,2.5]
        args['diff_lims'] = [-1.5,1.5]
        cmap_name = 'divergent'
        cmaps = [make_cmap(to_plot['thetao']['maui'], colours['cmaps'][cmap_name], abs_lims=args['abs_lims']),
                 make_cmap(to_plot['thetao']['maui'], colours['cmaps'][cmap_name], abs_lims=args['abs_lims']),
                 make_cmap(to_plot['thetao']['maui'], colours['cmaps'][cmap_name], abs_lims=args['abs_lims']),
                 make_cmap(to_plot['thetao']['ssp_hist'], colours['cmaps']['difference'], diff_lims=args['diff_lims']),
                 make_cmap(to_plot['thetao']['maui_hist'], colours['cmaps']['difference'], diff_lims=args['diff_lims']),
                 make_cmap(to_plot['thetao']['maui_ssp'], colours['cmaps']['difference'], diff_lims=args['diff_lims'])]
    
        cols = [cmap_name, cmap_name, cmap_name,'difference', 'difference', 'difference']
        
        vmm = [args['abs_lims'], args['abs_lims'], args['abs_lims'], args['diff_lims'], args['diff_lims'], args['diff_lims']]
        axes = {}
        fig = plt.figure(dpi=300, figsize=(5,15))
        gs = fig.add_gridspec(len(panels),1)
        #ax = fig.add_subplot(gs[0])
        #ax.set_ylabel('test')
        for igs, pn in enumerate(panels):
            axes[pn] = fig.add_subplot(gs[igs])
            axes[pn].set_ylabel('depth, $z$ / m', fontsize=14)
            if igs==len(panels)-1:
                axes[pn].set_xlabel('latitude, $\phi$ / $^\mathrm{o}$N', fontsize=14)
            else:
                axes[pn].set_xticklabels([])
            axes[pn].grid('major')

            axes[pn].invert_yaxis()
            axes[pn].set_facecolor((0.8,0.8,0.8))
            axes[pn].set_xlim((args['lats'][0], args['lats'][1]))            


            cmh = axes[pn].pcolormesh(Xs,Ys,to_plot['thetao'][pn], norm=cmaps[igs], cmap=colours['cmaps'][cols[igs]])
            if cols[igs]=='difference':
                fig.colorbar(cmh, ax=axes[pn], label=r'$\Delta \theta_\mathrm{o} / ^\mathrm{o}\mathrm{C}$', shrink=0.75, location='right')
            else:
                fig.colorbar(cmh, ax=axes[pn], label=r'$\theta_\mathrm{o} / ^\mathrm{o}\mathrm{C}$', shrink=0.75, location='right')
              
                
            
            print('doing sigma0')
            if pn=='hist':
                sig = axes['hist'].contour(Xs, Ys, to_plot['sigma0']['hist'], colors=colours['lines']['hist'], linewidths=section_linewidth, levels=levels)
                axes['hist'].clabel(sig, sig.levels, inline=True, fmt=fmt, fontsize=levels_fontsize)
            elif pn=='ssp':
                sig = axes['ssp'].contour(Xs, Ys, to_plot['sigma0']['ssp'], colors=colours['lines']['ssp'], linewidths=section_linewidth, levels=levels)
                axes['ssp'].clabel(sig, sig.levels, inline=True, fmt=fmt, fontsize=levels_fontsize)
            elif pn=='maui':
                sig = axes['maui'].contour(Xs, Ys, to_plot['sigma0']['maui'], colors=colours['lines']['maui'], linewidths=section_linewidth, levels=levels)
                axes['maui'].clabel(sig, sig.levels, inline=True, fmt=fmt, fontsize=levels_fontsize)
            elif pn=='ssp_hist':
                sig = axes['ssp_hist'].contour(Xs, Ys, to_plot['sigma0']['ssp'], colors=colours['lines']['ssp'], linewidths=section_linewidth, levels=levels)
                axes['ssp_hist'].clabel(sig, sig.levels, inline=True, fmt=fmt, fontsize=levels_fontsize)
                sig = axes['ssp_hist'].contour(Xs, Ys, to_plot['sigma0']['hist'], colors=colours['lines']['hist'], linewidths=section_linewidth, levels=levels)
                axes['ssp_hist'].clabel(sig, sig.levels, inline=True, fmt=fmt, fontsize=levels_fontsize)
            elif pn=='maui_hist':
                sig = axes['maui_hist'].contour(Xs, Ys, to_plot['sigma0']['maui'], colors=colours['lines']['maui'], linewidths=section_linewidth, levels=levels)
                axes['maui_hist'].clabel(sig, sig.levels, inline=True, fmt=fmt, fontsize=levels_fontsize)
                sig = axes['maui_hist'].contour(Xs, Ys, to_plot['sigma0']['hist'], colors=colours['lines']['hist'], linewidths=section_linewidth, levels=levels)
                axes['maui_hist'].clabel(sig, sig.levels, inline=True, fmt=fmt, fontsize=levels_fontsize)
            elif pn=='maui_ssp':
                sig = axes['maui_ssp'].contour(Xs, Ys, to_plot['sigma0']['ssp'], colors=colours['lines']['ssp'], linewidths=section_linewidth, levels=levels)
                axes['maui_ssp'].clabel(sig, sig.levels, inline=True, fmt=fmt, fontsize=levels_fontsize)
                sig = axes['maui_ssp'].contour(Xs, Ys, to_plot['sigma0']['maui'], colors=colours['lines']['maui'], linewidths=section_linewidth, levels=levels)
                axes['maui_ssp'].clabel(sig, sig.levels, inline=True, fmt=fmt, fontsize=levels_fontsize)
               
            
        fig.tight_layout()
        fig.savefig('ISMIPfw-section-thetao.png')
        
        if args['save_data']:
            for vn in ['thetao','so','sigma0']:
                data_file_name = 'ISMIPfw-data_for_section-'+vn+'-'
                for expid in panels:
                    if not vn=='sigma0':
                        to_plot[vn][expid].to_netcdf(data_file_name+expid+'.nc')
                    elif (expid=='hist') or (expid=='ssp') or (expid=='maui'):
                        to_plot[vn][expid].to_netcdf(data_file_name+expid+'.nc')
        
        
    if 'ts_profile' in args['analysis']:
        seasons = ['winter']
        # transform functions
        transform_functions = [partial(process_maui, model='oce'),
                               partial(lat_select, lat_low=args['lats'][0], lat_high=args['lats'][1]),
                               partial(lon_select, lon_low=args['lons'][0], lon_high=args['lons'][1]),
                               partial(depth_slice, deep=args['depths'][1], shallow=args['depths'][0])]#,
                               #partial(extract_season, season=args['season'])]

        cmip_so = {'historical': lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='historical', variable_id='so'),
                   'ssp585': lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='ssp585', variable_id='so')}
        cmip_thetao = {'historical': lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='historical', variable_id='thetao'),
                   'ssp585': lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='ssp585', variable_id='thetao')}
        
        
        # making data
        # cmip
        # we need so and thetao to calculate densities
       
        data_so = {'hist': cmip_so['historical'].copy(),
                'ssp':  cmip_so['ssp585'].copy()}
        data_thetao = {'hist': cmip_thetao['historical'].copy(),
                'ssp':  cmip_thetao['ssp585'].copy()}
        #to_plot = {}
        for eid in ['hist', 'ssp']:
            for tf in transform_functions[1:]:
                #data[eid] = tf(data[eid])
                data_so[eid] = tf(data_so[eid])
                data_thetao[eid] = tf(data_thetao[eid])
            #data[eid] = data[eid].sel(time=slice(cfdates[eid][0],cfdates[eid][1])).mean(dim=('member_id','time')) / cmip_units
            data_so[eid] = data_so[eid].sel(time=slice(cfdates[eid][0],cfdates[eid][1])).mean(dim=('member_id','x','y')) 
            data_thetao[eid] = data_thetao[eid].sel(time=slice(cfdates[eid][0],cfdates[eid][1])).mean(dim=('member_id','x','y'))
            #to_plot[eid] = data[eid].compute()
            print('Done ' + eid)
        
        print('starting maui load')
        #maui = {}
        maui_so = {}
        maui_thetao = {}
        for suite in args['suite']:
            if args['computer'] == 'laptop':
                maui_thetao[suite] = read_netcdfs(list_maui_files(suite=suite)[:13], drop_vars=vars2exclude(list_maui_files()[0], 'thetao'), transform_func=transform_functions)
                maui_so[suite] = read_netcdfs(list_maui_files(suite=suite)[:13], drop_vars=vars2exclude(list_maui_files()[0], 'so'), transform_func=transform_functions)
            else:
                maui_thetao[suite] = read_netcdfs(list_maui_files(suite=suite)[-maui_months:], drop_vars=vars2exclude(list_maui_files()[0], 'thetao'), transform_func=transform_functions)
                maui_so[suite] = read_netcdfs(list_maui_files(suite=suite)[-maui_months:], drop_vars=vars2exclude(list_maui_files()[0], 'so'), transform_func=transform_functions)
                
                
        
        data_so['maui'] = xr.concat([maui_so[sn] for sn in args['suite']],
               pd.Index([sn for sn in args['suite']], name="suite_id"))['so']
        
        data_thetao['maui'] = xr.concat([maui_thetao[sn] for sn in args['suite']],
               pd.Index([sn for sn in args['suite']], name="suite_id"))['thetao']
        
        data_so['maui'] = data_so['maui'].sel(time=slice(cfdates[eid][0],cfdates[eid][1])).mean(dim=('suite_id','x','y')) 
        data_thetao['maui'] = data_thetao['maui'].sel(time=slice(cfdates[eid][0],cfdates[eid][1])).mean(dim=('suite_id','x','y')) 



        to_plot = {'hist': {}, 'ssp': {}, 'maui': {}}
        for expid in ['hist', 'ssp', 'maui']:
            for season in seasons:
                to_plot[expid][season] = {'thetao': extract_season(data_thetao[expid], season).mean(dim='time').compute(),
                                          'so': extract_season(data_so[expid], season).mean(dim='time').compute(),
                                          'maui': extract_season(data_so[expid], season).mean(dim='time').compute()}
        
        fig = plt.figure(dpi=300, figsize=(4,8))
        gs = fig.add_gridspec(2,len(seasons))
        axes = {}
        for ip, season in enumerate(seasons):
            axes[season+'_so'] = fig.add_subplot(gs[0,ip])
            axes[season+'_so'].invert_yaxis()
            axes[season+'_thetao'] = axes[season+'_so'].twiny()
            axes[season+'_so'].set_xlim((32,35))
            axes[season+'_so'].set_xlabel(r'$S$ / g.kg$^{-1}$')
            #axes[season+'_so'].grid('major')
            axes[season+'_thetao'].set_xlim((-2.5,4))
            axes[season+'_thetao'].set_xlabel(r'$\theta_\mathrm{o}$ / $^\mathrm{o}$C')
            #axes[season+'_thetao'].grid('major')
            axes[season+'_so'].set_ylabel(r'$z$ / m$')
                        
            axes[season+'_thetao'].plot(to_plot['hist'][season]['thetao'], to_plot['hist'][season]['thetao'].lev, 'c')
            axes[season+'_thetao'].plot(to_plot['ssp'][season]['thetao'], to_plot['ssp'][season]['thetao'].lev, 'k')
            axes[season+'_thetao'].plot(to_plot['maui'][season]['thetao'], to_plot['ssp'][season]['thetao'].lev, 'm')

            axes[season+'_so'].plot([-999,-1000], [0,0], 'c', label=r'hist, $\theta_\mathrm{o}$')
            axes[season+'_so'].plot([-999,-1000], [0,0], 'k', label=r'ssp585, $\theta_\mathrm{o}$')
            axes[season+'_so'].plot([-999,-1000], [0,0], 'm', label=r'ssp585FW, $\theta_\mathrm{o}$')
        
            axes[season+'_so'].plot(to_plot['hist'][season]['so'], to_plot['hist'][season]['so'].lev, 'c:', label=r'hist, $S$')
            axes[season+'_so'].plot(to_plot['ssp'][season]['so'], to_plot['ssp'][season]['so'].lev, 'k:', label=r'ssp585, $S$')
            axes[season+'_so'].plot(to_plot['maui'][season]['so'], to_plot['maui'][season]['so'].lev, 'm:', label=r'ssp585FW, $S$')
            
            axes[season+'_sigma0'] = fig.add_subplot(gs[1,ip])
            axes[season+'_sigma0'].invert_yaxis()
            axes[season+'_sigma0'].set_xlim((26,28))
            axes[season+'_sigma0'].set_ylabel(r'$z$ / m$')
            axes[season+'_sigma0'].set_xlabel(r'$\sigma_\mathrm{0}$ / kg.m$^\mathrm{-3}$')
            axes[season+'_sigma0'].grid('both')
            
            axes[season+'_sigma0'].plot(sigmax(so=to_plot['hist'][season]['so'], thetao=to_plot['hist'][season]['thetao']), to_plot['hist'][season]['so'].lev, 'c', label=r'hist, $\sigma_\mathrm{0}$')
            axes[season+'_sigma0'].plot(sigmax(so=to_plot['ssp'][season]['so'], thetao=to_plot['ssp'][season]['thetao']), to_plot['ssp'][season]['so'].lev, 'k', label=r'ssp585, $\sigma_\mathrm{0}$')
            axes[season+'_sigma0'].plot(sigmax(so=to_plot['maui'][season]['so'], thetao=to_plot['maui'][season]['thetao']), to_plot['maui'][season]['so'].lev, 'm', label=r'ssp585FW, $\sigma_\mathrm{0}$')
            
            if ip==0:
                axes[season+'_so'].legend(fontsize=7)
                axes[season+'_sigma0'].legend(fontsize=7)
            
            
        
        fig.tight_layout()
        fig.savefig('ISMIPfw-ts_profile.png')
                
        
    if 'fw_timeseries' in args['analysis']:
        
        fig = plt.figure(dpi=300)
        gs = fig.add_gridspec(1,1)
        ax1 = fig.add_subplot(gs[0])
        cfdates = [datestring2cftime('20150101'), datestring2cftime('20991116')]
        
        Gta_to_kgs = 10**3 * 10**9 / (360 * 24 * 60 * 60)
        Gta_to_Sv  = 10**-6 * 10**9 / (360 * 24 * 60 * 60) # approx
        
        transform_functions = [partial(process_maui, model='oce'),
                               partial(lat_select, lat_low=-90, lat_high=-60)]
        
        areacello = lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='piControl', variable_id='areacello', table_id='Ofx')        
        areacello = transform_functions[1](areacello)
        data = {} 
        # cmip friver only available by wget so loading those in
        friver_at = '../../data/cmip_wgets/friver/data/friver_Omon_HadGEM3-GC31-LL_ssp585_'
        cmip_friver = {}
        ems = ['r1','r2','r3','r4']
        for em in ems:
            cmip1 = xr.open_dataset(friver_at+em+'i1p1f3_gn_201501-204912.nc')
            cmip2 = xr.open_dataset(friver_at+em+'i1p1f3_gn_205001-210012.nc')
            cmip_friver[em] = qpm(xr.concat((cmip1,cmip2), dim='time'))
        data['cmip_friver'] = xr.concat([cmip_friver[sn] for sn in ems],\
                       pd.Index([sn for sn in ems], name="ensemble_member"))
        print('before averaging')
        print(data['cmip_friver'])
        data['cmip_friver'] = transform_functions[1](data['cmip_friver']).sel(time=slice(cfdates[0],cfdates[1])).mean(dim='ensemble_member')['friver'] * areacello
        print('extensive and time selected')
        print(data['cmip_friver'])
        
        #data['cmip_friver'] = cmip_friver * areacello
        
        print('starting maui load friver')
        #maui = {}
        maui_friver = {}
        #maui_sowflisf = {}
        for suite in args['suite']:
            if args['computer'] == 'laptop':
                maui_friver[suite] = read_netcdfs(list_maui_files(suite=suite)[:13], drop_vars=vars2exclude(list_maui_files()[0], 'friver'), transform_func=transform_functions)
                #maui_sowflisf[suite] = read_netcdfs(list_maui_files(suite=suite)[:13], drop_vars=vars2exclude(list_maui_files()[0], 'so'), transform_func=transform_functions)
            else:
                maui_friver[suite] = read_netcdfs(list_maui_files(suite=suite), drop_vars=vars2exclude(list_maui_files()[0], 'friver'), transform_func=transform_functions)
                #maui_so[suite] = read_netcdfs(list_maui_files(suite=suite)[-maui_months:], drop_vars=vars2exclude(list_maui_files()[0], 'so'), transform_func=transform_functions)
                
                
        
        data['maui_friver'] = xr.concat([maui_friver[sn] for sn in args['suite']],
                       pd.Index([sn for sn in args['suite']], name="suite_id"))['friver'].sel(time=slice(cfdates[0],cfdates[1])).mean(dim='suite_id') * areacello
        
        if args['computer'] == 'laptop':
            data['maui_sowflisf'] = read_netcdfs(list_maui_files(suite='ci501')[:13], drop_vars=vars2exclude(list_maui_files()[0], 'sowflisf'), transform_func=transform_functions)['sowflisf'] * areacello
            #maui_sowflisf[suite] = read_netcdfs(list_maui_files(suite=suite)[:13], drop_vars=vars2exclude(list_maui_files()[0], 'so'), transform_func=transform_functions)
        else:
            data['maui_sowflisf'] = read_netcdfs(list_maui_files(suite='ci501'), drop_vars=vars2exclude(list_maui_files()[0], 'sowflisf'), transform_func=transform_functions)['sowflisf'] * areacello
        ax1.plot(data['cmip_friver'].time, data['cmip_friver'].sum(dim=('x','y')).rolling(time=12).mean() / Gta_to_kgs, colours['lines']['ssp'],label='runoff in SSP585')
        ax1.plot(data['maui_friver'].time, data['maui_friver'].sum(dim=('x','y')).rolling(time=12).mean() / Gta_to_kgs, colours['lines']['maui']+'--',label='runoff in SSP585FW')
        
        to_plot = {}
        #to_plot['friver_maui_cmip'] = (data['maui_friver'] - data['cmip_friver']).sum(dim=('x','y')) # difference in runoff as circumpolar sum (kg/s)
        to_plot['sowflisf_maui'] = (data['maui_sowflisf'].sum(dim=('x','y'))).compute()*-1  # sum of southern ocean basal melt and calved bergs (kg/s)
        
        PItotal = 5.693e7 # pre-industrial Antarctic fw (kg/s)
        basal_proportion = 0.55
        PIbasal = PItotal * basal_proportion 
        #plt.plot(to_plot['sowflisf_maui'].time, to_plot['sowflisf_maui'] / Gta_to_kgs, 'k:', label='Total basal melt')
        #plt.plot(to_plot['sowflisf_maui'].time, (to_plot['sowflisf_maui'] / basal_proportion) / Gta_to_kgs, 'k', label='Total freshwater forcing')
        ax1.plot(to_plot['sowflisf_maui'].time, (to_plot['sowflisf_maui'] - PIbasal) / Gta_to_kgs, 'k:', label='Extra basal melt in SSP585FW')
        ax1.plot(to_plot['sowflisf_maui'].time, ((to_plot['sowflisf_maui'] / basal_proportion) - PItotal) / Gta_to_kgs, 'k', label='Extra basal melt and calving in SSP585FW')



        
        
        
        experiments = ['exp05', 'exp07', 'expA3', 'expA4', 'exp01', 'exp03']
        descriptions = ['NorESM1, RCP8.5 (standard)', 'NorESM1, RCP2.6 (standard)', 'ISPL-CM5A-MR, RCP8.5', 'ISPL-CM5A-MR, RCP2.6', 'NorESM1, RCP8.5 (open)', 'NorESM, RCP2.6 (open)']
        linestyles = ['r','b','r--','b--','r:','b:']
        
        data_home = '../../data/seroussiEA20/ComputedScalarsPaper/'
        model_centers = glob.glob(data_home + '*')
        
        def get_experiment(exp_id='exp05', variable='bmbfl_minus'):
            basal_data = {}
            for moc in model_centers:
                models = glob.glob(moc + '/*')
                for mo in models:
                    runs = glob.glob(mo + '/*')
                    for ru in runs:
                        if exp_id in ru:
                            data_at = glob.glob(ru + '/*' + variable + '*')[0]
                            basal_data[ru] = xr.open_dataset(data_at, decode_times=False)
            
            models = []
            for key in basal_data.keys():
                models.append(key)
                
            times = basal_data[models[0]].time.values
                
            return basal_data, models, times
        
        exps = {}
        for exp in experiments:
            exps[exp] = {}
            exps[exp]['data'], exps[exp]['models'], exps[exp]['time'] = \
                get_experiment(exp_id = exp)
        
        
        def syr2dtsr(seroussi_years):
            return [datestring2cftime(str(xv)[:4]+'0601') for xv in seroussi_years]
        ## Replicate Seroussi et al Figure 12 for standard, medium ocean melt
        r85 = exps['exp05']
        xvar = r85['time']
        mean_basal = np.zeros_like(r85['time'])
        for mn in r85['models']:
            #plt.plot(r85['data'][mn].time.values, - r85['data'][mn].bmbfl.values / Gta_to_kgs, 'r-', alpha=0.1)
            ax1.plot(syr2dtsr(r85['data'][mn].time.values), - r85['data'][mn].bmbfl.values / Gta_to_kgs, 'r-', alpha=0.1)
            for ts, year in enumerate(r85['data'][mn].time.values):
                if year<2100: # we need this if statement because some of the models include 2100.5
                    mean_basal[ts] += r85['data'][mn].bmbfl.values[ts]
        
        mean_basal_85 = - mean_basal / Gta_to_kgs / len(r85['models'])
        ax1.plot(syr2dtsr(xvar), mean_basal_85, 'r', label='ISMIP6 basal melt projection')


        ax1.set_ylabel(r'freshwater input / Gt.a$^{-1}$')
        ax1.set_xlabel('Date')
        ax1.legend(fontsize=8)
        ax2 = ax1.twinx()
        mn, mx = ax1.get_ylim()
        ax2.set_ylim(mn*Gta_to_Sv, mx*Gta_to_Sv)
        ax2.set_ylabel(r'freshwater input / Sv')
        fig.savefig('fw_timeseries.png')
        
        # plt.figure(dpi=300)
        # ax1.plot(data['cmip_friver'].time, data['cmip_friver'].sum(dim=('x','y')).rolling(time=12).mean() / Gta_to_kgs, colours['lines']['ssp'], label='SSP585')
        # ax1.plot(data['maui_friver'].time, data['maui_friver'].sum(dim=('x','y')).rolling(time=12).mean() / Gta_to_kgs, colours['lines']['maui'], label='SSP585FW')
        # ax1.xlabel('Date')
        # ax1.ylabel('Runoff / Gt.a$^{-1}$')
        # ax1.legend()
        # ax1.savefig('runoff_timeseries.png')
        
        if args['save_data']:
            data_file_name = 'ISMIPfw-data_for_fw_timeseries-'
            data['cmip_friver'].sum(dim=('x','y')).rolling(time=12).mean().to_netcdf(data_file_name+'ssp585_runoff.nc')
            data['maui_friver'].sum(dim=('x','y')).rolling(time=12).mean().to_netcdf(data_file_name+'ssp585FW_runoff.nc')
            ((to_plot['sowflisf_maui'] - PIbasal) / Gta_to_kgs).to_netcdf(data_file_name+'SSP585FW_extra_basal.nc')
            

    if 'mackie' in args['analysis']:
        
        suites = ['ay187', 'az576', 'bb819']
        descriptions = ['FW', 'FWberg', 'FWCO2']
        mackie_months = 958
        mackie_months_end = 1200 
        cfdates = [datestring2cftime('19291101'), datestring2cftime('19491116')]
        map_depth = 400
        #projection = ccrs.NearsidePerspective(central_latitude=-90, satellite_height=3.6e6)
        transform_functions0 = [partial(process_maui, model='oce'),
                                partial(lat_select, lat_low=-90, lat_high=0),
                                partial(depth_slice, deep=1400, shallow=0)]
        
        projection = ccrs.SouthPolarStereo()
        # #panels = ['hist','ssp','maui','ssp_hist','maui_hist','maui_ssp']
        # fig = plt.figure(dpi=300, figsize=(12,18))
        # gs  = fig.add_gridspec(3,2)
        # gs.update(wspace=0.1, hspace=0.1)
        # axes = {}
        # for ip, panel in enumerate(panels):
        #     axes[panel] = fig.add_subplot(gs[ip], projection=projection)
        #     axes[panel].set_facecolor((0.8,0.8,0.8))
        #     axes[panel].gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
        #     #axes[panel].add_feature(cfeature.COASTLINE) 
        #     axes[panel].set_extent([-180, 180, -90, -65], ccrs.PlateCarree()) 
        
        # load pi control from cmip
        data = {'cmip': {} }
        for variable in ['thetao', 'so']:
            data['cmip'][variable] = lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='piControl', variable_id=variable).sel(time=slice(cfdates[0],cfdates[1]))
            for tf in transform_functions0[1:]:
                data['cmip'][variable] = tf(data['cmip'][variable])
                
        data['cmip_1pctCO2'] = {}
        for variable in ['thetao', 'so']:
            data['cmip_1pctCO2'][variable] = lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='1pctCO2', variable_id=variable).sel(time=slice(cfdates[0],cfdates[1]))
            for tf in transform_functions0[1:]:
                data['cmip_1pctCO2'][variable] = tf(data['cmip_1pctCO2'][variable])
        
        # load mackie data
        data['mackie'] = {}
        for suite in suites:
            data_at = '/nesi/project/nesi00442/thoma97p/UO_postdoc/data/' + suite + '/*'
            files = sorted(glob.glob(data_at))[mackie_months:mackie_months_end] # only load the last 2 years of monthly files as this is slow step and the data will get quite large
            drop_vars = vars2exclude(files[0], ['thetao','so'])
            data['mackie'][suite] = read_netcdfs(files, drop_vars=drop_vars, transform_func=transform_functions0).sel(time=slice(cfdates[0],cfdates[1]))
            
        # make map data
        maps = {'cmip': extract_depths(data['cmip']['thetao'], [map_depth]).mean(dim=('time','member_id')).compute()}
        for suite in suites:
            maps[suite] = extract_depths(data['mackie'][suite]['thetao'], [map_depth]).mean(dim=('time')).compute() - maps['cmip']
            maps[suite] = maps[suite].compute()
            #print(maps[suite]) 
        maps['cmip_1pctCO2'] = extract_depths(data['mackie']['bb819']['thetao'], [map_depth]).mean(dim=('time')).compute() \
                                - extract_depths(data['cmip_1pctCO2']['thetao'], [map_depth]).mean(dim=('time','member_id')).compute()
        # plot map data
        fig = plt.figure(dpi=300, figsize=(5,15))
        gs = fig.add_gridspec(len(suites) + 2,1)
        axes = {}
        panels = ['cmip'] + suites + ['cmip_1pctCO2']
        for ip, suite in enumerate(panels):
            # axes[suite] = fig.add_subplot(gs[igs], projection=projection)
            # axes[suite].set_facecolor((0.8,0.8,0.8))
            # axes[suite].gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
            # axes[suite].add_feature(cfeature.COASTLINE)
            axes[suite] = fig.add_subplot(gs[ip], projection=projection)
            axes[suite].set_facecolor((0.8,0.8,0.8))
            axes[suite].gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
            #axes[panel].add_feature(cfeature.COASTLINE) 
            axes[suite].set_extent([-180, 180, -90, -65], ccrs.PlateCarree()) 
            if ip == 0:
                vmin, vmax = -2.5, 7.5
                cmap = colours['cmaps']['divergent']
                label = r'$\theta$ / $^\mathrm{o}$C'
            else:
                vmin, vmax = -1.5, 1.5
                cmap = colours['cmaps']['difference']
                label = r'$\Delta\theta$ / $^\mathrm{o}$C'
            norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
            #axes[suite].scatter(maps[suite].lon, maps[suite].lat, c=maps[suite])
            cmh = pcm(maps[suite], norm=norm, cmap_name=cmap, ax=axes[suite])
            fig.colorbar(cmh, ax=axes[suite], shrink=0.75, location='right', label=label)
            
        if args['plot_shelf']:
            shelf, shelf2plot = make_shelf_mask()
            for ax in panels:
                axes[ax].contour(shelf.longitude, shelf.latitude, shelf2plot, transform=ccrs.PlateCarree(), levels=[1], linestyles=['-'], colors=[colours['lines']['shelf']])
         
            
        fig.tight_layout()
        fig.savefig('mackie_maps.png')
        
    if 'mackie_uo' in args['analysis']:
        
        suites = ['ay187', 'az576', 'bb819']
        descriptions = ['FW', 'FWberg', 'FWCO2']
        mackie_months = 958
        mackie_months_end = 1200 
        cfdates = [datestring2cftime('19291101'), datestring2cftime('19491116')]
        #print('TESTING DATES ONLY')
        map_depth = args['depths']
        projection = ccrs.NearsidePerspective(central_latitude=-90, satellite_height=3.6e6)
        transform_functions0 = [partial(process_maui, model='oce'),
                                partial(lat_select, lat_low=-90, lat_high=0),
                                partial(depth_slice, deep=1400, shallow=0)]
        
        # load pi control from cmip
        data = {'cmip': {} }
        for variable in ['uo']:
            data['cmip'][variable] = lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='piControl', variable_id=variable).sel(time=slice(cfdates[0],cfdates[1]))
            for tf in transform_functions0[1:]:
                data['cmip'][variable] = tf(data['cmip'][variable])
                
        data['cmip_1pctCO2'] = {}
        for variable in ['uo']:
            data['cmip_1pctCO2'][variable] = lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='1pctCO2', variable_id=variable).sel(time=slice(cfdates[0],cfdates[1]))
            for tf in transform_functions0[1:]:
                data['cmip_1pctCO2'][variable] = tf(data['cmip_1pctCO2'][variable])
        
        # load mackie data
        data['mackie'] = {}
        for suite in suites:
            data_at = '/nesi/project/nesi00442/thoma97p/UO_postdoc/data/mackie_Ugrid/' + suite + '/*'
            files = sorted(glob.glob(data_at))[mackie_months:mackie_months_end] # only load the last 2 years of monthly files as this is slow step and the data will get quite large
            drop_vars = vars2exclude(files[0], ['uo'])
            print(xr.open_dataset(files[0]))
            data['mackie'][suite] = read_netcdfs(files, drop_vars=drop_vars, transform_func=transform_functions0).sel(time=slice(cfdates[0],cfdates[1]))
        #print(data['mackie']['ay187'])    
        # make map data
        maps = {'cmip': extract_depths(data['cmip']['uo'], [map_depth]).mean(dim=('time','member_id')).compute()}
        for suite in suites:
            maps[suite] = extract_depths(data['mackie'][suite]['uo'], [map_depth]).mean(dim=('time')).compute() - maps['cmip']
            maps[suite] = maps[suite].compute()
            #print(maps[suite]) 
        maps['cmip_1pctCO2'] = extract_depths(data['mackie']['bb819']['uo'], [map_depth]).mean(dim=('time')).compute() \
                                - extract_depths(data['cmip_1pctCO2']['uo'], [map_depth]).mean(dim=('time','member_id')).compute()
        # plot map data
        fig = plt.figure(dpi=300, figsize=(5,15))
        gs = fig.add_gridspec(len(suites) + 2,1)
        axes = {}
        panels = ['cmip'] + suites + ['cmip_1pctCO2']
        for igs, suite in enumerate(panels):
            axes[suite] = fig.add_subplot(gs[igs], projection=projection)
            axes[suite].set_facecolor((0.8,0.8,0.8))
            axes[suite].gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
            axes[suite].add_feature(cfeature.COASTLINE)
            if igs == 0:
                vmin, vmax = -0.1, 0.1
                cmap = colours['cmaps']['divergent']
                label = r'$\theta$ / $^\mathrm{o}$C'
            else:
                vmin, vmax = -0.1, 0.1
                cmap = colours['cmaps']['difference']
                label = r'$\Delta\theta$ / $^\mathrm{o}$C'
            norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
            #axes[suite].scatter(maps[suite].lon, maps[suite].lat, c=maps[suite])
            cmh = pcm(maps[suite], norm=norm, cmap_name=cmap, ax=axes[suite])
            fig.colorbar(cmh, ax=axes[suite], shrink=0.75, location='right', label=label)
            
        if args['plot_shelf']:
            shelf, shelf2plot = make_shelf_mask()
            for ax in panels:
                axes[ax].contour(shelf.longitude, shelf.latitude, shelf2plot, transform=ccrs.PlateCarree(), levels=[1], linestyles=[':'], colors=[colours['lines']['shelf']])
         
            
        fig.tight_layout()
        fig.savefig('mackie_maps_uo.png')
        
        #if 'archive_test' in 
            
        
    if 'nearline_test' in args['analysis']:
        
        # transform functions
        transform_functions = [partial(process_maui, model=model),
                               partial(lat_select, lat_low=args['lats'][0], lat_high=args['lats'][1]),
                               partial(extract_depths, depths=args['depths']),
                               partial(extract_season, season=args['season'])]#,
        
        
        # map axes
        projection = ccrs.NearsidePerspective(central_latitude=-90, satellite_height=3.6e6)
        panels = ['atm','si','o_T','o_U','o_V','o_W', 'diaptr', 'scalar']
        fig = plt.figure(dpi=300, figsize=(12,7))
        gs  = fig.add_gridspec(3,len(panels))
        axes = {}
        for ip, panel in enumerate(panels):
            axes[panel] = fig.add_subplot(gs[ip], projection=projection)
            axes[panel].set_facecolor((0.3,0.3,0.3))
            axes[panel].gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
            axes[panel].add_feature(cfeature.COASTLINE)        
        # making data
        data = {'hist': cmip['historical'].copy(),
                'ssp':  cmip['ssp585'].copy()}
        to_plot = {}
        for eid in ['hist', 'ssp']:
            for tf in transform_functions[1:]:
                data[eid] = tf(data[eid])
            data[eid] = data[eid].sel(time=slice(cfdates[eid][0],cfdates[eid][1])).mean(dim='time') / cmip_units
            if variable_id == '1p5m_air_temperature':
                data[eid] = data[eid] - 273.15
            to_plot[eid] = data[eid].mean(dim='member_id').compute()
            print('Done ' + eid)
            
        print('starting maui load')
        maui = {}
        for suite in args['suite']:
            if model=='atm':
                maui[suite] = xr.DataArray.from_iris(iris.load(list_maui_files(suite=suite, model=model)[-maui_months:], \
                                iris.AttributeConstraint(STASH=get_um_stashcode(variable_id)))[0])
                for tf in transform_functions[:2]:
                    maui[suite] = tf(maui[suite])  
                if variable_id == '1p5m_air_temperature':
                    maui[suite] = maui[suite] - 273.15
            else:
                maui[suite] = read_netcdfs(list_maui_files(suite=suite, model=model, grid=grid)[-maui_months:], drop_vars=vars2exclude(list_maui_files(suite=suite, model=model, grid=grid)[0], variable_id), transform_func=transform_functions)[variable_id] 
        maui = xr.concat([maui[sn] for sn in args['suite']],
               pd.Index([sn for sn in args['suite']], name="suite_id"))
        maui = maui.sel(time=slice(cfdates[eid][0],cfdates[eid][1])).mean(dim='time')
        
        
        pvals = {}
        tstat, pvals['ssp_hist'] = ttest_ind(data['ssp'], data['hist'], axis=0, equal_var=False)
        tstat, pvals['maui_ssp'] = ttest_ind(maui, data['ssp'], axis=0, equal_var=False)
        
        to_plot['confident-ssp_hist'] = np.where(pvals['ssp_hist']<0.05, 1, 0)
        to_plot['confident-maui_ssp'] = np.where(pvals['maui_ssp']<0.05, 1, 0)
 
        to_plot['maui'] = maui.mean(dim=('suite_id')).compute()
        
        
        to_plot['ssp_hist'] = to_plot['ssp'] - to_plot['hist']
        to_plot['maui_ssp'] = to_plot['maui'] - to_plot['ssp']
        to_plot['maui_hist'] = to_plot['maui'] - to_plot['hist']
        if not 'lon' in to_plot['maui_ssp'].coords:
            to_plot['maui_ssp'] = to_plot['maui_ssp'].assign_coords({'lon': to_plot['hist'].lon}) # for some reason uo and vo are dropping lon coord during subtraction. This is just an easy fix
        cmaps = [make_cmap(to_plot['maui'], cmap_name, abs_lims=args['abs_lims']),
                 make_cmap(to_plot['maui'], cmap_name, abs_lims=args['abs_lims']),
                 make_cmap(to_plot['maui'], cmap_name, abs_lims=args['abs_lims']),
                 make_cmap(to_plot['ssp_hist'], 'seismic', diff_lims=args['diff_lims']),
                 make_cmap(to_plot['maui_ssp'], 'seismic', diff_lims=args['diff_lims']),#]
                 make_cmap(to_plot['maui_hist'], 'seismic', diff_lims=args['diff_lims'])]

        cmabs = pcm(to_plot['hist'], cmaps[0], cmap_name, axes['hist'])
        cmabs = pcm(to_plot['ssp'], cmaps[1], cmap_name, axes['ssp'])
        cmabs = pcm(to_plot['maui'], cmaps[2], cmap_name, axes['maui'])
        cmhdif = pcm(to_plot['ssp_hist'], cmaps[3], 'seismic', axes['ssp_hist'])
        cmpdif = pcm(to_plot['maui_hist'], cmaps[4], 'seismic', axes['maui_ssp'])
        cmpdif = pcm(to_plot['maui_ssp'], cmaps[5], 'seismic', axes['maui_hist'])
        
        fig.colorbar(cmabs, ax=axes['hist'], shrink=0.5, location='top', label=variable_info.at[variable_id,'label1'])
        fig.colorbar(cmabs, ax=axes['maui'], shrink=0.5, location='top', label=variable_info.at[variable_id,'label1'])
        fig.colorbar(cmabs, ax=axes['ssp'], shrink=0.5, location='top', label=variable_info.at[variable_id,'label1'])
        fig.colorbar(cmhdif, ax=axes['ssp_hist'], shrink=0.5, location='top',label='$\Delta$'+variable_info.at[variable_id,'label1'])
        fig.colorbar(cmpdif, ax=axes['maui_ssp'], shrink=0.5, location='top',label='$\Delta$'+variable_info.at[variable_id,'label1'])
        fig.colorbar(cmpdif, ax=axes['maui_hist'], shrink=0.5, location='top',label='$\Delta$'+variable_info.at[variable_id,'label1'])
        
        fig.tight_layout()
        fig.savefig('ISMIPfw-map-'+variable_id+'.png')
        

    if 'ekman' in args['analysis']:
        #variable_id = 'ekman'    
        print('starting wind stress curl map')
        # transform functions
        transform_functions = [partial(process_maui, model=model),
                               partial(lat_select, lat_low=args['lats'][0], lat_high=args['lats'][1]),
                               partial(extract_season, season=args['season'])]#,
        
        
        cmip_tauuo = {'hist': lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='historical', variable_id='tauuo'),
                             'ssp': lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='ssp585', variable_id='tauuo')}
        cmip_tauvo = {'hist': lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='historical', variable_id='tauvo'),
                             'ssp': lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='ssp585', variable_id='tauvo')}
        for eid in ['hist', 'ssp']:
            for tf in transform_functions[1:]:
                cmip_tauuo[eid] = tf(cmip_tauuo[eid])
                cmip_tauvo[eid] = tf(cmip_tauvo[eid])
            cmip_tauuo[eid] = cmip_tauuo[eid].sel(time=slice(cfdates[eid][0],cfdates[eid][1])).mean(dim=('time','member_id')).compute()
            cmip_tauvo[eid] = cmip_tauvo[eid].sel(time=slice(cfdates[eid][0],cfdates[eid][1])).mean(dim=('time','member_id')).compute()
            print('Done ' + eid)
        
        print('starting maui load')
        maui_tauuo = {}
        maui_tauvo = {}
        for suite in args['suite']:
            maui_tauuo[suite] = read_netcdfs(list_maui_files(suite=suite, model='oce', grid='U')[-maui_months:], drop_vars=vars2exclude(list_maui_files(suite=suite, model='oce', grid='U')[0], 'tauuo'), transform_func=transform_functions)['tauuo'] 
            maui_tauvo[suite] = read_netcdfs(list_maui_files(suite=suite, model='oce', grid='V')[-maui_months:], drop_vars=vars2exclude(list_maui_files(suite=suite, model='oce', grid='V')[0], 'tauvo'), transform_func=transform_functions)['tauvo'] 
            
        maui_tauuo = xr.concat([maui_tauuo[sn] for sn in args['suite']],
               pd.Index([sn for sn in args['suite']], name="suite_id"))
        maui_tauuo = maui_tauuo.sel(time=slice(cfdates[eid][0],cfdates[eid][1])).mean(dim=('time','suite_id')).compute()
        
        maui_tauvo = xr.concat([maui_tauvo[sn] for sn in args['suite']],
               pd.Index([sn for sn in args['suite']], name="suite_id"))
        maui_tauvo = maui_tauvo.sel(time=slice(cfdates[eid][0],cfdates[eid][1])).mean(dim=('time','suite_id')).compute()
        
        
        print('starting wind stress curl calculation')
        e1v = meshmask['e1v'][0,1:-1,1:-1]
        e2u = meshmask['e2u'][0,1:-1,1:-1]
        
        ekm = {'hist': cmip_tauuo['hist'].copy(),
               'ssp':     cmip_tauuo['ssp'].copy(),
               'maui':       maui_tauuo.copy()}
        ekm_historical = calculate_ekman(cmip_tauuo['hist'], cmip_tauvo['hist'], e1v, e2u)
        ekm_ssp585     = calculate_ekman(cmip_tauuo['ssp'], cmip_tauvo['ssp'], e1v, e2u)
        ekm_maui       = calculate_ekman(maui_tauuo, maui_tauvo, e1v, e2u)
        
        ekm['hist'].values = ekm_historical
        ekm['ssp'].values     = ekm_ssp585
        ekm['maui'].values       = ekm_maui
        
        to_plot = {} 
        to_plot['hist'] = ekm['hist'].copy()
        to_plot['ssp'] = ekm['ssp'].copy()
        to_plot['maui'] = ekm['maui'].copy()
        to_plot['ssp_hist'] = to_plot['ssp'] - to_plot['hist']
        to_plot['maui_hist'] = to_plot['maui'] - to_plot['hist']
        to_plot['maui_ssp'] = to_plot['maui'] - to_plot['ssp']
        
        
        # map axes
        # projection = ccrs.NearsidePerspective(central_latitude=-90, satellite_height=3.6e6)
        projection = ccrs.SouthPolarStereo()
        # panels = ['hist','ssp','maui','ssp_hist','maui_hist','maui_ssp']
        # fig = plt.figure(dpi=300, figsize=(12,7))
        # gs  = fig.add_gridspec(2,3)
        # axes = {}
        # for ip, panel in enumerate(panels):
        #     axes[panel] = fig.add_subplot(gs[ip], projection=projection)
        #     axes[panel].set_facecolor((0.8,0.8,0.8))
        #     axes[panel].gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
        #     axes[panel].add_feature(cfeature.COASTLINE)       
        panels = ['hist','ssp_hist','ssp','maui_hist','maui','maui_ssp']
        fig, axes = make_map_axes(panels)
        

        if not 'lon' in to_plot['maui_ssp'].coords:
            to_plot['maui_ssp'] = to_plot['maui_ssp'].assign_coords({'lon': to_plot['hist'].lon}) # for some reason uo and vo are dropping lon coord during subtraction. This is just an easy fix
        cmaps = [make_cmap(to_plot['maui'], 'divergent', abs_lims=args['abs_lims']),
                 make_cmap(to_plot['maui'], 'divergent', abs_lims=args['abs_lims']),
                 make_cmap(to_plot['maui'], 'divergent', abs_lims=args['abs_lims']),
                 make_cmap(to_plot['ssp_hist'], 'difference', diff_lims=args['diff_lims']),
                 make_cmap(to_plot['maui_ssp'], 'difference', diff_lims=args['diff_lims']),#]
                 make_cmap(to_plot['maui_hist'], 'difference', diff_lims=args['diff_lims'])]

        cmabs = pcm(to_plot['hist'], cmaps[0], colours['cmaps']['divergent'], axes['hist'])
        cmabs = pcm(to_plot['ssp'], cmaps[1], colours['cmaps']['divergent'], axes['ssp'])
        cmabs = pcm(to_plot['maui'], cmaps[2], colours['cmaps']['divergent'], axes['maui'])
        cmhdif = pcm(to_plot['ssp_hist'], cmaps[3], colours['cmaps']['difference'], axes['ssp_hist'])
        cmpdif = pcm(to_plot['maui_hist'], cmaps[4], colours['cmaps']['difference'], axes['maui_hist'])
        cmpdif = pcm(to_plot['maui_ssp'], cmaps[5], colours['cmaps']['difference'], axes['maui_ssp'])
 
        axes['hist'].plot([252,252],[-76,-68],'m-',transform=ccrs.PlateCarree())
        axes['hist'].plot([177,177],[-78.2,-70],'m-',transform=ccrs.PlateCarree())
        axes['hist'].plot([71,71],[-69,-65],'m-',transform=ccrs.PlateCarree())
        axes['hist'].plot([210,210], [-75,-60], 'm:',transform=ccrs.PlateCarree())
        axes['hist'].plot([290,290], [-75,-60], 'm:',transform=ccrs.PlateCarree())
        axes['hist'].plot([170,170], [-75,-60], 'm:',transform=ccrs.PlateCarree())
        axes['hist'].plot([180,180], [-75,-60], 'm:',transform=ccrs.PlateCarree())
        axes['hist'].plot([120,120], [-75,-60], 'm:',transform=ccrs.PlateCarree())
        axes['hist'].plot([60,60], [-75,-60], 'm:',transform=ccrs.PlateCarree())    
        # fig.colorbar(cmabs, ax=axes['hist'], shrink=0.5, location='right', label=r'$w_\mathrm{Ek} /$ m.s$^{-1}$')
        # fig.colorbar(cmabs, ax=axes['maui'], shrink=0.5, location='right', label=r'$w_\mathrm{Ek} /$ m.s$^{-1}$')
        # fig.colorbar(cmabs, ax=axes['ssp'], shrink=0.5, location='right', label=r'$w_\mathrm{Ek} /$ m.s$^{-1}$')
        # fig.colorbar(cmhdif, ax=axes['ssp_hist'], shrink=0.5, location='right',label=r'$\Delta$'+'$w_\mathrm{Ek} /$ m.s$^{-1}$')
        # fig.colorbar(cmpdif, ax=axes['maui_ssp'], shrink=0.5, location='right',label=r'$\Delta$'+'$w_\mathrm{Ek} /$ m.s$^{-1}$')
        # fig.colorbar(cmpdif, ax=axes['maui_hist'], shrink=0.5, location='right',label=r'$\Delta$'+'$w_\mathrm{Ek} /$ m.s$^{-1}$')
        
        if args['plot_sia_contour']:
            cmip_sia, maui_sia = get_sea_ice_area(args)
            for expid in ['hist', 'ssp']:
                cmip_sia[expid] = cmip_sia[expid].mean(dim=('member_id','time')).compute()
            maui_sia = maui_sia.mean(dim=('suite_id','time')).compute()
            
            to_plot_sia = {'hist': cmip_sia['hist'],
                           'ssp': cmip_sia['ssp'],
                           'maui': maui_sia}
            to_plot_sia['ssp_hist'] = to_plot_sia['ssp'] - to_plot_sia['hist']
            to_plot_sia['maui_hist'] = to_plot_sia['maui'] - to_plot_sia['hist']
            to_plot_sia['maui_ssp'] = to_plot_sia['maui'] - to_plot_sia['ssp']
            
            plot_sea_ice_contour(to_plot_sia['hist'], axes['hist'], color='c')
            plot_sea_ice_contour(to_plot_sia['hist'], axes['ssp_hist'], color='c')
            plot_sea_ice_contour(to_plot_sia['hist'], axes['maui_hist'], color='c')
            
            plot_sea_ice_contour(to_plot_sia['ssp'], axes['ssp'], color='m')
            plot_sea_ice_contour(to_plot_sia['ssp'], axes['ssp_hist'], color='m')
            plot_sea_ice_contour(to_plot_sia['ssp'], axes['maui_ssp'], color='m')
            
            plot_sea_ice_contour(to_plot_sia['maui'], axes['maui'])
            plot_sea_ice_contour(to_plot_sia['maui'], axes['maui_hist'])
            plot_sea_ice_contour(to_plot_sia['maui'], axes['maui_ssp'])
        
        if args['plot_shelf']:
            shelf, shelf2plot = make_shelf_mask()
            for ax in panels:
                axes[ax].contour(shelf.longitude, shelf.latitude, shelf2plot, transform=ccrs.PlateCarree(), levels=[1], linestyles=['-'], colors=[colours['lines']['shelf']])
        
        
        fig.tight_layout()
        fig.savefig('ISMIPfw-map-ekman.png')
    
        fig.colorbar(cmabs, ax=axes['maui'], shrink=0.75, location='bottom', label=r'$w_\mathrm{Ek} /$ m.s$^{-1}$')
        fig.colorbar(cmpdif, ax=axes['maui_ssp'], shrink=0.75, location='bottom',label=r'$\Delta$'+'$w_\mathrm{Ek} /$ m.s$^{-1}$')
        fig.savefig('ISMIPfw-map-ekman-colorbars.png')
    
        if args['save_data']:
            data_file_name = 'ISMIPfw-data_for_map-ekman-'
            for expid in panels:
                to_plot[expid].to_netcdf(data_file_name+expid+'.nc')
        
        areacello = lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='piControl', variable_id='areacello', table_id='Ofx')    
        
        weights = areacello.fillna(0)
        mask_1000 = shelf_mask(areacello)
        
        print('printing shelf masked weighted mean ekman')
        
        for lp in lon_pairs:
            fig = plt.figure(dpi=300)
            gs = fig.add_gridspec(1,3)
            mask_lon = (areacello.lon>=lp[0]) & (areacello.lon<lp[1])
            mask = (mask_1000) & (mask_lon)
            to_print = {}
            print('for %s>=lon>%s:' % (lp[0], lp[1]))
            for sp, expid in enumerate(['hist','ssp','maui']):
                to_print[expid] = to_plot[expid].where(mask)
                to_print[expid] = to_print[expid].assign_coords({'lon': areacello['lon'], 'lat': areacello['lat']})
                print('%s: %s' % (expid, to_print[expid].weighted(weights).mean(('x','y')).values))
                ax = fig.add_subplot(gs[sp], projection=projection)
                ax.set_facecolor((0.8,0.8,0.8))
                ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
                ax.add_feature(cfeature.COASTLINE)  
                pcm(to_print[expid], cmaps[sp], 'PiYG', ax)
            fig.savefig('ISMIPfw-map-ekman-'+args['season']+str(lp[0])+'-'+str(lp[1])+'.png')
            
                
                
        
        
    if 'sidmassth' in args['analysis']:
        print('starting sidmassth map')
        variable_id = 'sidmassth'
        variable_info = pd.read_csv('../../documents/notebooks/variable_lists/variables_for_analysis.csv', index_col='name')
        model = variable_info.at[variable_id, 'model']
        grid  = variable_info.at[variable_id, 'grid']
        label = variable_info.at[variable_id,'label1']
        table_id = model[0].upper() + 'mon'
        cmap_name  = variable_info.at[variable_id, 'cmap']
        maui_months = 240
        
        # transform functions
        transform_functions = [partial(process_maui, model=model),
                               partial(lat_select, lat_low=args['lats'][0], lat_high=args['lats'][1]),
                               partial(extract_season, season=args['season'])]#,
        
        cmip_ssp585_at = '../../data/cmip_wgets/sidmassth/data/sidmassth_SImon_HadGEM3-GC31-LL_ssp585_r'
        cmip_hist_at   = '../../data/cmip_wgets/sidmassth/data/sidmassth_SImon_HadGEM3-GC31-LL_historical_r'
        
        cmip = {}
        for isuite in range(4):
            cmip_str = cmip_ssp585_at+str(isuite+1)+'i1p1f3_gn_205001-210012.nc'
            cmip[str(isuite)] = xr.open_dataset(cmip_str).sel(time=slice(cfdates['ssp'][0],cfdates['ssp'][1]))['sidmassth']
            cmip[str(isuite)] = cmip[str(isuite)].rename({'latitude':'lat', 'longitude':'lon', 'i':'x', 'j':'y'})
            
        cmip_ssp = xr.concat([cmip[str(sn)] for sn in range(4)],
                   pd.Index([str(sn) for sn in range(4)], name="member_id"))
            
        cmip = {}
        for isuite in range(5):
            cmip_str = cmip_hist_at+str(isuite+1)+'i1p1f3_gn_195001-201412.nc'
            cmip[str(isuite)] = xr.open_dataset(cmip_str).sel(time=slice(cfdates['hist'][0],cfdates['hist'][1]))['sidmassth']
            cmip[str(isuite)] = cmip[str(isuite)].rename({'latitude':'lat', 'longitude':'lon', 'i':'x', 'j':'y'})
            
        cmip_hist = xr.concat([cmip[str(sn)] for sn in range(5)],
                   pd.Index([str(sn) for sn in range(5)], name="member_id"))
        
        for tf in transform_functions[1:]:
            cmip_ssp = tf(cmip_ssp)
            cmip_hist = tf(cmip_hist)
        maui = {}
        for suite in args['suite']:
            file_list = list_maui_files(suite=suite, model=model)[-maui_months:]
            maui[suite] = read_netcdfs(file_list, drop_vars=vars2exclude(file_list[0], variable_id), transform_func=transform_functions)[variable_id] 
        maui = xr.concat([maui[sn] for sn in args['suite']],
               pd.Index([sn for sn in args['suite']], name="suite_id"))
        maui = maui.sel(time=slice(cfdates['ssp'][0],cfdates['ssp'][1]))
        
        to_plot = {'hist': cmip_hist.mean(dim=('time','member_id')).compute(),
                   'ssp': cmip_ssp.mean(dim=('time','member_id')).compute(),
                   'maui': maui.mean(dim=('time','suite_id')).compute()}
        to_plot['ssp_hist'] = to_plot['ssp'] - to_plot['hist']
        #to_plot['maui_hist'] = to_plot['maui'] - to_plot['hist']
        #to_plot['maui_ssp'] = to_plot['maui'] - to_plot['ssp']
        to_plot['maui_hist'] = to_plot['hist'].copy()
        to_plot['maui_hist'].values = (to_plot['maui'] - to_plot['hist']).values
        to_plot['maui_ssp'] = to_plot['hist'].copy()
        to_plot['maui_ssp'].values = (to_plot['maui'] - to_plot['ssp']).values
   
        # map axes
        #projection = ccrs.NearsidePerspective(central_latitude=-90, satellite_height=3.6e6)
        projection = ccrs.SouthPolarStereo()
        panels = ['hist','ssp_hist','ssp','maui_hist','maui','maui_ssp']
        fig, axes = make_map_axes(panels)
        # fig = plt.figure(dpi=300, figsize=(12,7))
        # gs  = fig.add_gridspec(2,3)
        # axes = {}
        # for ip, panel in enumerate(panels):
        #     axes[panel] = fig.add_subplot(gs[ip], projection=projection)
        #     axes[panel].set_facecolor((0.8,0.8,0.8))
        #     axes[panel].gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
        #     axes[panel].add_feature(cfeature.COASTLINE)    
        
        

        

        if not 'lon' in to_plot['maui_ssp'].coords:
            to_plot['maui_ssp'] = to_plot['maui_ssp'].assign_coords({'lon': to_plot['hist'].lon}) # for some reason uo and vo are dropping lon coord during subtraction. This is just an easy fix
        cmaps = [make_cmap(to_plot['maui'], 'divergent', abs_lims=args['abs_lims']),
                 make_cmap(to_plot['maui'], 'divergent', abs_lims=args['abs_lims']),
                 make_cmap(to_plot['maui'], 'divergent', abs_lims=args['abs_lims']),
                 make_cmap(to_plot['ssp_hist'], 'difference', diff_lims=args['diff_lims']),
                 make_cmap(to_plot['maui_ssp'], 'difference', diff_lims=args['diff_lims']),#]
                 make_cmap(to_plot['maui_hist'], 'difference', diff_lims=args['diff_lims'])]

        cmabs = pcm(to_plot['hist'], cmaps[0], colours['cmaps']['divergent'], axes['hist'])
        cmabs = pcm(to_plot['ssp'], cmaps[1], colours['cmaps']['divergent'], axes['ssp'])
        cmabs = pcm(to_plot['maui'], cmaps[2], colours['cmaps']['divergent'], axes['maui'])
        cmhdif = pcm(to_plot['ssp_hist'], cmaps[3], colours['cmaps']['difference'], axes['ssp_hist'])
        cmpdif = pcm(to_plot['maui_hist'], cmaps[4], colours['cmaps']['difference'], axes['maui_hist'])
        cmpdif = pcm(to_plot['maui_ssp'], cmaps[5], colours['cmaps']['difference'], axes['maui_ssp'])
        
        axes['hist'].plot([252,252],[-76,-68],'m-',transform=ccrs.PlateCarree())
        axes['hist'].plot([177,177],[-78.2,-70],'m-',transform=ccrs.PlateCarree())
        axes['hist'].plot([71,71],[-69,-65],'m-',transform=ccrs.PlateCarree()) 
        axes['hist'].plot([210,210], [-75,-60], 'm:',transform=ccrs.PlateCarree())
        axes['hist'].plot([290,290], [-75,-60], 'm:',transform=ccrs.PlateCarree())
        axes['hist'].plot([170,170], [-75,-60], 'm:',transform=ccrs.PlateCarree())
        axes['hist'].plot([180,180], [-75,-60], 'm:',transform=ccrs.PlateCarree())
        axes['hist'].plot([120,120], [-75,-60], 'm:',transform=ccrs.PlateCarree())
        axes['hist'].plot([60,60], [-75,-60], 'm:',transform=ccrs.PlateCarree())
        # fig.colorbar(cmabs, ax=axes['hist'], shrink=0.5, location='right', label='$\delta m_\mathrm{si}$ / kg.m$^{-2}$s$^{-1}$')
        # fig.colorbar(cmabs, ax=axes['maui'], shrink=0.5, location='right', label='$\delta m_\mathrm{si}$ / kg.m$^{-2}$s$^{-1}$')
        # fig.colorbar(cmabs, ax=axes['ssp'], shrink=0.5, location='right', label='$\delta m_\mathrm{si}$ / kg.m$^{-2}$s$^{-1}$')
        # fig.colorbar(cmhdif, ax=axes['ssp_hist'], shrink=0.5, location='right',label='$\Delta$'+'$\Delta m_\mathrm{si}$ / kg.m$^{-2}$s$^{-1}$')
        # fig.colorbar(cmpdif, ax=axes['maui_ssp'], shrink=0.5, location='right',label='$\Delta$'+'$\Delta m_\mathrm{si}$ / kg.m$^{-2}$s$^{-1}$')
        # fig.colorbar(cmpdif, ax=axes['maui_hist'], shrink=0.5, location='right',label='$\Delta$'+'$\Delta m_\mathrm{si}$ / kg.m$^{-2}$s$^{-1}$')
        
        if args['plot_sia_contour']:
            cmip_sia, maui_sia = get_sea_ice_area(args)
            for expid in ['hist', 'ssp']:
                cmip_sia[expid] = cmip_sia[expid].mean(dim=('member_id','time')).compute()
            maui_sia = maui_sia.mean(dim=('suite_id','time')).compute()
            
            to_plot_sia = {'hist': cmip_sia['hist'],
                           'ssp': cmip_sia['ssp'],
                           'maui': maui_sia}
            to_plot_sia['ssp_hist'] = to_plot_sia['ssp'] - to_plot_sia['hist']
            to_plot_sia['maui_hist'] = to_plot_sia['maui'] - to_plot_sia['hist']
            to_plot_sia['maui_ssp'] = to_plot_sia['maui'] - to_plot_sia['ssp']
            
            plot_sea_ice_contour(to_plot_sia['hist'], axes['hist'], color='c')
            plot_sea_ice_contour(to_plot_sia['hist'], axes['ssp_hist'], color='c')
            plot_sea_ice_contour(to_plot_sia['hist'], axes['maui_hist'], color='c')
            
            plot_sea_ice_contour(to_plot_sia['ssp'], axes['ssp'], color='m')
            plot_sea_ice_contour(to_plot_sia['ssp'], axes['ssp_hist'], color='m')
            plot_sea_ice_contour(to_plot_sia['ssp'], axes['maui_ssp'], color='m')
            
            plot_sea_ice_contour(to_plot_sia['maui'], axes['maui'])
            plot_sea_ice_contour(to_plot_sia['maui'], axes['maui_hist'])
            plot_sea_ice_contour(to_plot_sia['maui'], axes['maui_ssp'])
 
        if args['plot_shelf']:
            shelf, shelf2plot = make_shelf_mask()
            for ax in panels:
                axes[ax].contour(shelf.longitude, shelf.latitude, shelf2plot, transform=ccrs.PlateCarree(), levels=[1], linestyles=['-'], colors=[colours['lines']['shelf']])

        fig.tight_layout()
        fig.savefig('ISMIPfw-map-sidmassth.png')
    
        cbar = fig.colorbar(cmabs, ax=axes['maui'], shrink=0.75, location='bottom', label=variable_info.at[variable_id,'label1'])
        cbar.ax.tick_params(rotation=45)
        cbar = fig.colorbar(cmpdif, ax=axes['maui_ssp'], shrink=0.75, location='bottom',label='$\Delta$'+variable_info.at[variable_id,'label1']) 
        cbar.ax.tick_params(rotation=45)
        fig.savefig('ISMIPfw-map-sidmassth-colorbars.png')
    
        if args['save_data']:
            data_file_name = 'ISMIPfw-data_for_map-sidmassth-'
            for expid in panels:
                to_plot[expid].to_netcdf(data_file_name+expid+'.nc')
        
    
        areacello = lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='piControl', variable_id='areacello', table_id='Ofx').isel(member_id=0)
        to_plot_Xarea = {'hist': (to_plot['hist'] * areacello).compute(),
                         'ssp': (to_plot['ssp'] * areacello).compute(),
                         'maui': (to_plot['maui'] * areacello).compute()}
        weights = areacello.fillna(0)
        mask_1000 = shelf_mask(areacello)

        print('printing shelf masked sum sidmassth')
        print(areacello)
        for lp in lon_pairs:
            fig = plt.figure(dpi=300)
            gs = fig.add_gridspec(1,3)
            mask_lon = (areacello.lon>=lp[0]) & (areacello.lon<lp[1])
            mask = (mask_1000) & (mask_lon)
            to_print = {}
            print('for %s>=lon>%s:' % (lp[0], lp[1]))
            for sp, expid in enumerate(['hist','ssp','maui']):
                to_print[expid] = to_plot_Xarea[expid].where(mask)
                to_print[expid] = to_print[expid].assign_coords({'lon': areacello['lon'], 'lat': areacello['lat']})
                print('%s: %s' % (expid, to_print[expid].sum(('x','y')).values))
                ax = fig.add_subplot(gs[sp], projection=projection)
                ax.set_facecolor((0.8,0.8,0.8))
                ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
                ax.add_feature(cfeature.COASTLINE)
                pcm(to_print[expid], cmaps[sp], 'PiYG', ax)
                #xr.plot.pcolormesh(to_print[expid], 'lon', 'lat', transform=ccrs.PlateCarree(), norm=cmaps[sp], cmap='YlGn', ax=ax, add_colorbar=False, add_labels=False)
            fig.savefig('ISMIPfw-map-sidmassth-'+args['season'][0]+str(lp[0])+'-'+str(lp[1])+'.png')
   

    if 'sea_ice_area' in args['analysis']:
        
        variable_id = 'soicecov'
        
        cmip, maui = get_sea_ice_area(args)
        
        to_plot = {'hist': cmip['hist'].mean(dim=('time','member_id')).compute(),
                   'ssp': cmip['ssp'].mean(dim=('time','member_id')).compute(),
                   'maui': maui.mean(dim=('time','suite_id')).compute()}
        
        to_plot['ssp_hist'] = to_plot['ssp'] - to_plot['hist']
        #to_plot['maui_hist'] = to_plot['maui'] - to_plot['hist']
        #to_plot['maui_ssp'] = to_plot['maui'] - to_plot['ssp']
        to_plot['maui_hist'] = to_plot['hist'].copy()
        to_plot['maui_hist'].values = (to_plot['maui'] - to_plot['hist']).values
        to_plot['maui_ssp'] = to_plot['hist'].copy()
        to_plot['maui_ssp'].values = (to_plot['maui'] - to_plot['ssp']).values
   
        # map axes
        # projection = ccrs.NearsidePerspective(central_latitude=-90, satellite_height=3.6e6)
        projection = ccrs.SouthPolarStereo()
        panels = ['hist','ssp_hist','ssp','maui_hist','maui','maui_ssp']
        fig, axes = make_map_axes(panels)
        # fig = plt.figure(dpi=300, figsize=(12,7))
        # gs  = fig.add_gridspec(2,3)
        # axes = {}
        # for ip, panel in enumerate(panels):
        #     axes[panel] = fig.add_subplot(gs[ip], projection=projection)
        #     axes[panel].set_facecolor((0.8,0.8,0.8))
        #     axes[panel].gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
        #     axes[panel].add_feature(cfeature.COASTLINE)       

        

        if not 'lon' in to_plot['maui_ssp'].coords:
            to_plot['maui_ssp'] = to_plot['maui_ssp'].assign_coords({'lon': to_plot['hist'].lon}) # for some reason uo and vo are dropping lon coord during subtraction. This is just an easy fix
        cmaps = [make_cmap(to_plot['maui'], 'YlGn', abs_lims=args['abs_lims']),
                 make_cmap(to_plot['maui'], 'YlGn', abs_lims=args['abs_lims']),
                 make_cmap(to_plot['maui'], 'YlGn', abs_lims=args['abs_lims']),
                 make_cmap(to_plot['ssp_hist'], 'seismic', diff_lims=args['diff_lims']),
                 make_cmap(to_plot['maui_ssp'], 'seismic', diff_lims=args['diff_lims']),#]
                 make_cmap(to_plot['maui_hist'], 'seismic', diff_lims=args['diff_lims'])]

        cmabs = pcm(to_plot['hist'], cmaps[0], colours['cmaps']['divergent'], axes['hist'])
        cmabs = pcm(to_plot['ssp'], cmaps[1], colours['cmaps']['divergent'], axes['ssp'])
        cmabs = pcm(to_plot['maui'], cmaps[2], colours['cmaps']['divergent'], axes['maui'])
        cmhdif = pcm(to_plot['ssp_hist'], cmaps[3], colours['cmaps']['difference'], axes['ssp_hist'])
        cmpdif = pcm(to_plot['maui_hist'], cmaps[4], colours['cmaps']['difference'], axes['maui_hist'])
        cmpdif = pcm(to_plot['maui_ssp'], cmaps[5], colours['cmaps']['difference'], axes['maui_ssp'])
        
        axes['hist'].plot([252,252],[-76,-68],'m-',transform=ccrs.PlateCarree())
        axes['hist'].plot([177,177],[-78.2,-70],'m-',transform=ccrs.PlateCarree())
        axes['hist'].plot([71,71],[-69,-65],'m-',transform=ccrs.PlateCarree())
        axes['hist'].plot([210,210], [-75,-60], 'm:',transform=ccrs.PlateCarree())
        axes['hist'].plot([290,290], [-75,-60], 'm:',transform=ccrs.PlateCarree())
        axes['hist'].plot([170,170], [-75,-60], 'm:',transform=ccrs.PlateCarree())
        axes['hist'].plot([180,180], [-75,-60], 'm:',transform=ccrs.PlateCarree())
        axes['hist'].plot([120,120], [-75,-60], 'm:',transform=ccrs.PlateCarree())
        axes['hist'].plot([60,60], [-75,-60], 'm:',transform=ccrs.PlateCarree())
        
        #cbar_fig = plt.figure(dpi=300)
        #fig.colorbar(cmabs, ax=axes['hist'], shrink=0.5, location='right', label=variable_info.at[variable_id,'label1'])
        #fig.colorbar(cmabs, ax=axes['maui'], shrink=0.5, location='right', label=variable_info.at[variable_id,'label1'])
        #fig.colorbar(cmabs, ax=axes['ssp'], shrink=0.5, location='right', label=variable_info.at[variable_id,'label1'])
        #fig.colorbar(cmhdif, ax=axes['ssp_hist'], shrink=0.5, location='right',label='$\Delta$'+variable_info.at[variable_id,'label1'])
        #fig.colorbar(cmpdif, ax=axes['maui_ssp'], shrink=0.5, location='right',label='$\Delta$'+variable_info.at[variable_id,'label1'])
        #fig.colorbar(cmpdif, ax=axes['maui_hist'], shrink=0.5, location='right',label='$\Delta$'+variable_info.at[variable_id,'label1'])
        
        if args['plot_sia_contour']:
            cmip_sia, maui_sia = get_sea_ice_area(args)
            for expid in ['hist', 'ssp']:
                cmip_sia[expid] = cmip_sia[expid].mean(dim=('member_id','time')).compute()
            maui_sia = maui_sia.mean(dim=('suite_id','time')).compute()
            
            to_plot_sia = {'hist': cmip_sia['hist'],
                           'ssp': cmip_sia['ssp'],
                           'maui': maui_sia}
            to_plot_sia['ssp_hist'] = to_plot_sia['ssp'] - to_plot_sia['hist']
            to_plot_sia['maui_hist'] = to_plot_sia['maui'] - to_plot_sia['hist']
            to_plot_sia['maui_ssp'] = to_plot_sia['maui'] - to_plot_sia['ssp']
            
            plot_sea_ice_contour(to_plot_sia['hist'], axes['hist'], color='c')
            plot_sea_ice_contour(to_plot_sia['hist'], axes['ssp_hist'], color='c')
            plot_sea_ice_contour(to_plot_sia['hist'], axes['maui_hist'], color='c')
            
            plot_sea_ice_contour(to_plot_sia['ssp'], axes['ssp'], color='m')
            plot_sea_ice_contour(to_plot_sia['ssp'], axes['ssp_hist'], color='m')
            plot_sea_ice_contour(to_plot_sia['ssp'], axes['maui_ssp'], color='m')
            
            plot_sea_ice_contour(to_plot_sia['maui'], axes['maui'])
            plot_sea_ice_contour(to_plot_sia['maui'], axes['maui_hist'])
            plot_sea_ice_contour(to_plot_sia['maui'], axes['maui_ssp'])
 
        if args['plot_shelf']:
            shelf, shelf2plot = make_shelf_mask()
            for ax in panels:
                axes[ax].contour(shelf.longitude, shelf.latitude, shelf2plot, transform=ccrs.PlateCarree(), levels=[1], linestyles=['-'], colors=[colours['lines']['shelf']])

        
        
        
        
        fig.tight_layout()
        fig.savefig('ISMIPfw-map-soicecov.png')
        fig.colorbar(cmabs, ax=axes['maui'], shrink=0.75, location='bottom', label=variable_info.at[variable_id,'label1'])
        fig.colorbar(cmpdif, ax=axes['maui_ssp'], shrink=0.75, location='bottom',label='$\Delta$'+variable_info.at[variable_id,'label1'])
        fig.savefig('ISMIPfw-map-soicecov-colorbars.png')
 
        if args['save_data']:
            data_file_name = 'ISMIPfw-data_for_map-sia-'
            for expid in panels:
                to_plot[expid].to_netcdf(data_file_name+expid+'.nc')
        
    
        areacello = lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='piControl', variable_id='areacello', table_id='Ofx').isel(member_id=0)
        to_plot_Xarea = {'hist': (to_plot['hist'] * areacello).compute(),
                         'ssp': (to_plot['ssp'] * areacello).compute(),
                         'maui': (to_plot['maui'] * areacello).compute()}
        weights = areacello.fillna(0)
        mask_1000 = shelf_mask(areacello)

        print('printing shelf masked sum sea-ice area')
        print(areacello)
        for lp in lon_pairs:
            fig = plt.figure(dpi=300)
            gs = fig.add_gridspec(1,3)
            mask_lon = (areacello.lon>=lp[0]) & (areacello.lon<lp[1])
            mask = (mask_1000) & (mask_lon)
            to_print = {}
            print('for %s>=lon>%s:' % (lp[0], lp[1]))
            for sp, expid in enumerate(['hist','ssp','maui']):
                to_print[expid] = to_plot_Xarea[expid].where(mask)
                to_print[expid] = to_print[expid].assign_coords({'lon': areacello['lon'], 'lat': areacello['lat']})
                print('%s: %s' % (expid, to_print[expid].sum(('x','y')).values))
                ax = fig.add_subplot(gs[sp], projection=projection)
                ax.set_facecolor((0.8,0.8,0.8))
                ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
                ax.add_feature(cfeature.COASTLINE)
                pcm(to_print[expid], cmaps[sp], 'YlGn', ax)
                #xr.plot.pcolormesh(to_print[expid], 'lon', 'lat', transform=ccrs.PlateCarree(), norm=cmaps[sp], cmap='YlGn', ax=ax, add_colorbar=False, add_labels=False)
            fig.savefig('ISMIPfw-map-soicecov-'+args['season'][0]+str(lp[0])+'-'+str(lp[1])+'.png')
   

 
    if 'test_contour' in args['analysis']:
        ff = '../../data/nemo_ci501o_1m_20150101-20150201_grid-T.nc'
        fh = process_maui(xr.open_dataset(ff))
        aice = fh['soicecov'][0]
              
        lons, lats = sea_ice_contour(aice)
        
        fig = plt.figure()
        gs = fig.add_gridspec(1,1)
        ax = fig.add_subplot(gs[0], projection=ccrs.NearsidePerspective(central_latitude=-90, satellite_height=3.6e6))
        ax.set_facecolor((0.8,0.8,0.8))
        ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
        ax.add_feature(cfeature.COASTLINE)
        ax.scatter(lons,lats,c='k',s=0.01, transform=ccrs.PlateCarree())
        
        
    if 'test_regrid' in args['analysis']:
        
        
        
        
        
        ff = ['../../data/nemo_ci501o_1m_20150101-20150201_grid-T.nc',
              '../../data/nemo_ci501o_1m_20150801-20150901_grid-T.nc']
        # test = read_and_regrid(ff, variable_id='soicecov')
        # test=test.mean(dim='time_counter')
        
        # SIE = np.where((test>0.15), 1, 0)
        #SIE = np.where(test,np.isclose(test.values, np.ones_like(test.values)*0.85, atol=0.01))
        
        sic, SIE = regridded_seaice_contour(ff)
        
        fig = plt.figure()
        gs = fig.add_gridspec(1,1)
        ax = fig.add_subplot(gs[0], projection=ccrs.NearsidePerspective(central_latitude=-90, satellite_height=3.6e6))
        ax.set_facecolor((0.3,0.3,0.3))
        ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
        ax.add_feature(cfeature.COASTLINE)
        
        
        
        shelf, shelf2plot = make_shelf_mask()
        ax.contour(shelf.longitude, shelf.latitude, shelf2plot, transform=ccrs.PlateCarree(), levels=[1], linestyles=[':'], colors=['y'])
        ax.contour(sic.lon, sic.lat, SIE, transform=ccrs.PlateCarree(), levels=[1], colors=['k'])
        
    if 'test_heat_budget' in args['analysis']:
        # load data and make time means
        data_at = '../../data/'
        #Tdata = xr.open_dataset('~/Desktop/nemo_ci501o_1m_20981201-20990101_grid-T.nc')
        Tdata = xr.open_dataset(data_at + 'nemo_ci501o_1m_20150101-20150201_grid-T.nc')
        #Vdata = xr.open_dataset(data_at + 'nemo_ci501o_1m_20981201-20990101_grid-V.nc')
        #Udata = xr.open_dataset(data_at + 'nemo_ci501o_1m_20981201-20990101_grid-U.nc')
        
        # get opottemptend and opottempadvect
        tend = Tdata.opottemptend 
        adve = Tdata.opottempadvect 
        area = Tdata.area
        
       
       
        
        # apply transform functions
        transform_functions = [partial(process_maui, model=model),
                               partial(lat_select, lat_low=args['lats'][0], lat_high=args['lats'][1]),
                               partial(extract_season, season=args['season'])]#,
        
        area = area.assign_coords({'lat': tend.nav_lat, 'lon': tend.nav_lon})
        area = area[1:-1,1:-1]

        mask1000 = shelf_mask(area)
        mask1000 = transform_functions[1](mask1000).astype(bool)

        
        area = transform_functions[1](area)
        
        for tf in transform_functions:
            tend = tf(tend)
            adve = tf(adve)
        
        
    
        #DELETE LATER, just gets rid of time and depth dimensions for testing
        tend=tend[0,0]
        adve=adve[0,0]
            
        # mask shelf
        shelf, shelf2plot = make_shelf_mask()
        
        # make plot data
        to_plot = {'tend': tend,
                   'adve': adve}
        
        cmaps = [make_cmap(to_plot['tend'], 'divergent', abs_lims=args['abs_lims']),
                 make_cmap(to_plot['adve'], 'divergent', abs_lims=args['abs_lims'])]
        

        # map axes
        projection = ccrs.NearsidePerspective(central_latitude=-90, satellite_height=3.6e6)
        #projection = ccrs.SouthPolarStereo()
        panels = ['tend','adve']
        fig = plt.figure(dpi=300, figsize=(12,7))
        gs  = fig.add_gridspec(1,2)
        axes = {}
        masklon = (tend.lon>=200) & (tend.lon<300)
        
        mask = (masklon * mask1000).astype(bool)
    
        for ip, panel in enumerate(panels):
            axes[panel] = fig.add_subplot(gs[ip], projection=projection)
            axes[panel].set_facecolor((0.8,0.8,0.8))
            axes[panel].gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
            axes[panel].add_feature(cfeature.COASTLINE) 
            cmap = pcm(to_plot[panel].where(mask), cmaps[ip], colours['cmaps']['divergent'], axes[panel])
            fig.colorbar(cmap, ax=axes[panel], shrink=0.5, location='right', label=panel+' / W/m2')


    if 'woa' in args['analysis']:
        

        def load_woa(section, variable_name):
            woa = xr.open_dataset('../../data/WOA_sections/woa-'+section+'-'+variable_name+'.nc', decode_times=False)
            try:
                return woa[variable_name+'_an'].isel(time=0)
            except:
                return woa[variable_name+'_an']
            
        cfdates = [datestring2cftime('19810101'), datestring2cftime('20100101')]
        cmip = {'to': lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='historical', variable_id='thetao'),
                'so': lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='historical', variable_id='so')}
        
        for vn in ['to','so']:
            cmip[vn] = cmip[vn].sel(time=slice(cfdates[0],cfdates[1]))
        
        
        # make figure
        
        cmap_name = 'continuous' 
        sections = ['W.Ross']
        plot_features = {'Thwaites': {'lon': 251, 'lats':(-76,-68), 'title':'Thwaites -- 252 E'},
                         'Amery': {'lon': 71, 'lats':(-73,-65), 'title':'Amery -- 71 E'},
                         'W.Ross': {'lon': 177, 'lats':(-78.2,-70), 'title':'W. Ross -- 177 E'}}
        levels = [26.4,26.6,26.8,27.0,27.2,27.4,27.6,27.8,28.0,28.2,28.4,28.6]
       
        
        args['abs_lims'] = [33,35]
        fig_so = plt.figure(dpi=300,figsize=(18,12))
        gs = fig_so.add_gridspec(2,3)
        axes_so = {'WOA': {}, 'cmip': {}}
        for iax, section in enumerate(sections):
            axes_so['WOA'][section] = fig_so.add_subplot(gs[0,iax])
            axes_so['WOA'][section].set_ylim([1000,0])
            axes_so['WOA'][section].set_xlim(plot_features[section]['lats'])
            axes_so['WOA'][section].set_title(plot_features[section]['title'])
            axes_so['WOA'][section].set_facecolor((0.8,0.8,0.8))
            
        for iax, section in enumerate(sections):
            axes_so['cmip'][section] = fig_so.add_subplot(gs[1,iax])
            axes_so['cmip'][section].set_ylim([1000,0])
            axes_so['cmip'][section].set_xlim(plot_features[section]['lats'])
            axes_so['cmip'][section].set_facecolor((0.8,0.8,0.8))
            #axes['maui'][section].set_title(plot_features[section]['title'])
            

        print('woa plotting') 
        # plot WOA
        woa_sections_so = {}
        for section in sections:
            woa_sections_so[section] = load_woa(section, 's')
            woa_sections_so[section] = gsw.SA_from_SP(SP=woa_sections_so[section], p=woa_sections_so[section]['depth'], lon=plot_features[section]['lon'], lat=woa_sections_so[section]['lat'])
            cmap = make_cmap(woa_sections_so[section], colours['cmaps']['continuous'], abs_lims=args['abs_lims'])
            cbr = axes_so['WOA'][section].pcolormesh(woa_sections_so[section]['lat'],woa_sections_so[section]['depth'],woa_sections_so[section], norm=cmap, cmap=colours['cmaps']['continuous'])
            fig_so.colorbar(cbr, ax=axes_so['WOA'][section], label=r'$S_\mathrm{o} /$ g/kg', shrink=0.75, location='right')
        print('cmip plotting')
        # plot cmip historical
        cmip_sections_so = {}
        for section in sections:
            lon = plot_features[section]['lon']
            lats = plot_features[section]['lats']
            
            cmip_sections_so[section] = cmip['so'].copy()
            cmip_sections_so[section] = lat_select(cmip_sections_so[section],lats[0], lats[1])
            cmip_sections_so[section] = lon_select(cmip_sections_so[section],lon-1, lon+1)
            lat2ass = cmip_sections_so[section].lat
            
            cmip_sections_so[section] = cmip_sections_so[section].mean(dim=('time','x','member_id'))           
             
            cbr = axes_so['cmip'][section].pcolormesh(lat2ass.mean(dim='x'),cmip_sections_so[section]['lev'],cmip_sections_so[section], norm=cmap, cmap=colours['cmaps']['continuous'])
            fig_so.colorbar(cbr, ax=axes_so['cmip'][section], label=r'$S_\mathrm{o} / $g/kg', shrink=0.75, location='right')
            
        
        
        # temperature
        args['abs_lims'] = [-2.5,2.5]
        fig_to = plt.figure(dpi=300,figsize=(18,12))
        gs = fig_to.add_gridspec(2,3)
        axes_to = {'WOA': {}, 'cmip': {}}
        for iax, section in enumerate(sections):
            axes_to['WOA'][section] = fig_to.add_subplot(gs[0,iax])
            axes_to['WOA'][section].set_ylim([1000,0])
            axes_to['WOA'][section].set_xlim(plot_features[section]['lats'])
            axes_to['WOA'][section].set_title(plot_features[section]['title'])
            axes_to['WOA'][section].set_facecolor((0.8,0.8,0.8))
            
        for iax, section in enumerate(sections):
            axes_to['cmip'][section] = fig_to.add_subplot(gs[1,iax])
            axes_to['cmip'][section].set_ylim([1000,0])
            axes_to['cmip'][section].set_xlim(plot_features[section]['lats'])
            axes_to['cmip'][section].set_facecolor((0.8,0.8,0.8))
            #axes['maui'][section].set_title(plot_features[section]['title'])
        print('starting woa woa fig')
        # plot WOA
        woa_sections_to = {}
        for section in sections:
            woa_sections_to[section] = load_woa(section, 't')
            woa_sections_to[section] = gsw.CT_from_t(SA=woa_sections_so[section], t=woa_sections_to[section], p=woa_sections_to[section]['depth'])
            cmap = make_cmap(woa_sections_to[section], colours['cmaps'][cmap_name], abs_lims=args['abs_lims'])
            cbr = axes_to['WOA'][section].pcolormesh(woa_sections_to[section]['lat'],woa_sections_to[section]['depth'],woa_sections_to[section], norm=cmap, cmap=colours['cmaps']['divergent'])
            fig_to.colorbar(cbr, ax=axes_to['WOA'][section], label=r'$\theta_\mathrm{o} / ^\mathrm{o}\mathrm{C}$', shrink=0.75, location='right')
            
        print('starting cmip woa fig')
        # plot cmip historical
        cmip_sections_to = {}
        for section in sections:
            lon = plot_features[section]['lon']
            lats = plot_features[section]['lats']
            
            cmip_sections_to[section] = cmip['to'].copy()
            cmip_sections_to[section] = lat_select(cmip_sections_to[section],lats[0], lats[1])
            cmip_sections_to[section] = lon_select(cmip_sections_to[section],lon-1, lon+1)
            lat2ass = cmip_sections_to[section].lat

            cmip_sections_to[section] = cmip_sections_to[section].mean(dim=('time','x','member_id'))
            
            cbr = axes_to['cmip'][section].pcolormesh(lat2ass.mean(dim='x'),cmip_sections_to[section]['lev'],cmip_sections_to[section], norm=cmap, cmap=colours['cmaps']['divergent'])
            fig_to.colorbar(cbr, ax=axes_to['cmip'][section], label=r'$\theta_\mathrm{o} / ^\mathrm{o}\mathrm{C}$', shrink=0.75, location='right')
        print('starting sigmas woa fig')
        woa_sections_sigma0={}
        cmip_sections_sigma0={}
        for section in sections:
            woa_sections_sigma0[section] = gsw.sigma0(woa_sections_so[section], woa_sections_to[section])
            cmip_sections_sigma0[section] = gsw.sigma0(cmip_sections_so[section], cmip_sections_to[section])
            
            #lat2ass = cmip_sections_so[section].lat.mean(dim='x')
            #lat2ass_sig = lat2ass.mean(dim='x')

            sig = axes_so['WOA'][section].contour(woa_sections_sigma0[section].lat,woa_sections_sigma0[section].depth,woa_sections_sigma0[section], levels=levels, colors='k')
            axes_so['WOA'][section].clabel(sig, sig.levels, inline=True, fmt=fmt, fontsize=9)
           
            sig =  axes_so['cmip'][section].contour(lat2ass.mean(dim='x'),cmip_sections_sigma0[section].lev,cmip_sections_sigma0[section], levels=levels, colors='k')
            axes_so['cmip'][section].clabel(sig, sig.levels, inline=True, fmt=fmt, fontsize=9)            
            
            sig = axes_to['WOA'][section].contour(woa_sections_sigma0[section].lat,woa_sections_sigma0[section].depth,woa_sections_sigma0[section], levels=levels, colors='k')
            axes_to['WOA'][section].clabel(sig, sig.levels, inline=True, fmt=fmt, fontsize=9)   
            
            sig = axes_to['cmip'][section].contour(lat2ass.mean(dim='x'),cmip_sections_sigma0[section].lev,cmip_sections_sigma0[section], levels=levels, colors='k')
            axes_to['cmip'][section].clabel(sig, sig.levels, inline=True, fmt=fmt, fontsize=9)   
            
        print('starting save woa fig')    
        data_file_name = 'ISMIPfw-obs_comparison-'
        lat2ass.mean(dim='x').to_netcdf(data_file_name+'lats.nc')
        woa_sections_so['W.Ross'].to_netcdf(data_file_name+'woa_so.nc')
        woa_sections_to['W.Ross'].to_netcdf(data_file_name+'woa_to.nc')  
        woa_sections_sigma0['W.Ross'].to_netcdf(data_file_name+'woa_sigma0.nc')  
        cmip_sections_so['W.Ross'].to_netcdf(data_file_name+'cmip_so.nc')  
        cmip_sections_to['W.Ross'].to_netcdf(data_file_name+'cmip_to.nc')  
        cmip_sections_sigma0['W.Ross'].to_netcdf(data_file_name+'cmip_sigma0.nc')      
        fig_so.savefig('S_woa.png')
        fig_to.savefig('T_woa.png')

    if 'section_uo' in args['analysis']:        
        
        if 'mackie' in args['suite']:
            Mackie = True
            args['suite'] = args['suite'][1:] # so if you want to plot ay187 (FW), do --suites=mackie,ay187  can be: ['ay187', 'az576', 'bb819']
            cfdates = {'hist': [datestring2cftime('19291101'), datestring2cftime('19491116')],
                       'ssp':  [datestring2cftime('19291101'), datestring2cftime('19491116')]}
        else:
            Mackie = False 
        # transform functions
        transform_functions = [partial(process_maui, model='oce'),
                               partial(lat_select, lat_low=args['lats'][0], lat_high=args['lats'][1]),
                               partial(lon_select, lon_low=args['lons'][0], lon_high=args['lons'][1]),
                               partial(depth_slice, deep=args['depths'][1], shallow=args['depths'][0]),
                               partial(extract_season, season=args['season'])]
        
        if not Mackie:
            cmip_uo = {'historical': lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='historical', variable_id='uo'),
                       'ssp585': lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='ssp585', variable_id='uo')}
        else: # this is a bit of a fudge to get mackie data plotting with the same code. the historical control becomes pi, and the ssp585 control becomes 1%co2
            print('loading PI and 1pctCO2')
            cmip_uo = {'historical': lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='piControl', variable_id='uo'),
                       'ssp585': lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='1pctCO2', variable_id='uo')}
            cmip_thetao = {'historical': lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='piControl', variable_id='thetao'),
                       'ssp585': lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='1pctCO2', variable_id='thetao')}
            
        
        # making data
        # cmip
        # we need so and thetao to calculate densities
        # data = {'hist': cmip['historical'].copy(),
        #         'ssp':  cmip['ssp585'].copy()}
        data_uo = {'hist': cmip_uo['historical'].copy(),
                'ssp':  cmip_uo['ssp585'].copy()}

        #to_plot = {}
        for eid in ['hist', 'ssp']:
            for tf in transform_functions[1:]:
                #data[eid] = tf(data[eid])
                data_uo[eid] = tf(data_uo[eid])
            #data[eid] = data[eid].sel(time=slice(cfdates[eid][0],cfdates[eid][1])).mean(dim=('member_id','time')) / cmip_units
            data_uo[eid] = data_uo[eid].sel(time=slice(cfdates[eid][0],cfdates[eid][1])).mean(dim=('member_id','time')) 
            #to_plot[eid] = data[eid].compute()
            print('Done ' + eid)
        
        print('starting maui load')
        #maui = {}
        maui_uo = {}
        for suite in args['suite']:
            if not Mackie:
                #maui[suite] = read_netcdfs(list_maui_files(suite=suite)[-maui_months:], drop_vars=vars2exclude(list_maui_files()[0], variable_id), transform_func=transform_functions)
                maui_uo[suite] = read_netcdfs(list_maui_files(suite=suite,grid='U')[-maui_months:], drop_vars=vars2exclude(list_maui_files(grid='U')[0], 'uo'), transform_func=transform_functions)
            else:
                file_list = sorted(glob.glob('../../data/'+suite+'/nemo_'+suite+'o_1m*'))[-maui_months:]
                maui_uo[suite] = read_netcdfs(file_list, drop_vars=vars2exclude(list_maui_files()[0], 'uo'), transform_func=transform_functions)
        print(maui_uo)
        maui_uo = xr.concat([maui_uo[sn] for sn in args['suite']],
               pd.Index([sn for sn in args['suite']], name="suite_id"))
       
        maui_uo = maui_uo['uo'] 
        
        maui_uo = maui_uo.sel(time=slice(cfdates[eid][0],cfdates[eid][1])).mean(dim=('suite_id','time')).compute()        

        data_uo['hist']['lat'] = data_uo['hist']['lat'].mean(dim='x')
        data_uo['ssp']['lat'] = data_uo['ssp']['lat'].mean(dim='x')
        maui_uo['lat'] = maui_uo['lat'].mean(dim='x')
       
        to_plot = {}
        to_plot['uo'] = {
                         'hist': data_uo['hist'].mean(dim='x').compute(),
                         'ssp': data_uo['ssp'].mean(dim='x').compute(),
                         'maui': maui_uo.mean(dim='x').compute()
                        }
        to_plot['uo']['ssp_hist'] = to_plot['uo']['ssp'] - to_plot['uo']['hist']
        to_plot['uo']['maui_ssp'] = to_plot['uo']['maui'] - to_plot['uo']['ssp']
        to_plot['uo']['maui_hist'] = to_plot['uo']['maui'] - to_plot['uo']['hist']
        
    
        #args['abs_lims'] = [33,35]
        #args['diff_lims'] = [-1.5,1.5]
        cmap_name = 'continuous'       
        cmaps = [make_cmap(to_plot['uo']['maui'], colours['cmaps'][cmap_name], abs_lims=args['abs_lims']),
                 make_cmap(to_plot['uo']['maui'], colours['cmaps'][cmap_name], abs_lims=args['abs_lims']),
                 make_cmap(to_plot['uo']['maui'], colours['cmaps'][cmap_name], abs_lims=args['abs_lims']),
                 make_cmap(to_plot['uo']['ssp_hist'], colours['cmaps']['difference'], diff_lims=args['diff_lims']),
                 make_cmap(to_plot['uo']['maui_ssp'], colours['cmaps']['difference'], diff_lims=args['diff_lims']),
                 make_cmap(to_plot['uo']['maui_hist'], colours['cmaps']['difference'], diff_lims=args['diff_lims'])]
    
        cols = [cmap_name, cmap_name, cmap_name, 'difference', 'difference', 'difference']
        
        vmm = [args['abs_lims'], args['abs_lims'], args['abs_lims'], args['diff_lims'], args['diff_lims'], args['diff_lims']]
    
        panels = ['hist', 'ssp', 'maui', 'ssp_hist', 'maui_hist', 'maui_ssp']
        Xs, Ys = np.meshgrid(to_plot['uo']['maui'].lat, to_plot['uo']['maui'].lev)
        axes = {}
        fig = plt.figure(dpi=300, figsize=(5,15))
        gs = fig.add_gridspec(len(panels),1)
        #ax = fig.add_subplot(gs[0])
        #ax.set_ylabel('test')
        for igs, pn in enumerate(panels):
            axes[pn] = fig.add_subplot(gs[igs])
            axes[pn].set_ylabel('depth, $z$ / m', fontsize=14)
            if igs==len(panels)-1:
                axes[pn].set_xlabel('latitude, $\phi$ / $^\mathrm{o}$N', fontsize=14)
            else:
                axes[pn].set_xticklabels([])
            axes[pn].grid('major')

            axes[pn].invert_yaxis()
            axes[pn].set_facecolor((0.8,0.8,0.8))
            axes[pn].set_xlim((args['lats'][0], args['lats'][1]))            

            print('start panel '+pn)

            cmh = axes[pn].pcolormesh(Xs,Ys,to_plot['uo'][pn], norm=cmaps[igs], cmap=colours['cmaps'][cols[igs]])
            if cols[igs]=='seismic':
                fig.colorbar(cmh, ax=axes[pn], label=r'$\Delta u$ / m.s$^{-1}$', shrink=0.75, location='right')
            else:
                fig.colorbar(cmh, ax=axes[pn], label=r'$u$ / m.s$^{-1}$', shrink=0.75, location='right')
              
            
        fig.tight_layout()
        fig.savefig('ISMIPfw-section-uo.png')        
        
        
    if 'latent_heat' in args['analysis']:
        vn = 'vowflisf'
        LH = 333.55 * 10**3 # latent heat of fusion J/kg
        
        dummy_file = process_maui(xr.open_dataset('../../data/nemo_ci501o_1m_20150101-20150201_grid-T.nc'))
        mask = shelf_mask(dummy_file)
        
        def apply_shelf_mask(ds,mask=None):
            if mask is None:
                ds = ds.where(shelf_mask(ds), drop=True)
            else:
                ds = ds.where(mask, drop=True)
            return ds
        
        transform_functions = [partial(process_maui, model='oce'),
                              partial(apply_shelf_mask,mask=mask),
                              partial(lat_select, lat_low=args['lats'][0], lat_high=args['lats'][1]),
                              partial(lon_select, lon_low=args['lons'][0], lon_high=args['lons'][1]),
                              partial(depth_slice, deep=args['depths'][1], shallow=args['depths'][0])]  
 
        # load maui data (only one suite, all the same as prescribed forcing)
        cfdates = [datestring2cftime('20150101'), datestring2cftime('20991116')]
        
        areacello = lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='piControl', variable_id='areacello', table_id='Ofx')  
        for tf in transform_functions[1:4]:
            areacello = tf(areacello)
            
        
        # maui vowflisf only available by wget so loading those in (files like nemo_ci501o_1m_20981201-20990101_grid-T.nc)
        maui = read_netcdfs(list_maui_files(suite='ci501'), drop_vars=vars2exclude(list_maui_files()[0], 'vowflisf'), transform_func=transform_functions)['vowflisf']
              
        # mulitply by area
        maui = areacello * maui # kg/m2/s * m2 -> kg/s
        
        
        #fig = plt.figure(dpi=300)
        #gs = fig.add_gridspec(1,1)
        #ax = fig.add_subplot(gs[0], projection=ccrs.NearsidePerspective(central_latitude=-90, satellite_height=3.6e6))
        #ax.add_feature(cfeature.COASTLINE)  
        #ax.pcolormesh(maui.lon,maui.lat,maui.isel(time=0).sum(dim='lev'),transform=ccrs.PlateCarree())
        #fig.savefig('mask_test-lh.png')
        maui.to_netcdf('lh_spatial.nc')        

        # sum
        summed = maui.sum(('x','y','lev')) # kg/s
        
        summed.to_netcdf('lh_spatial-sum'+str(args['depths'][1])+'_lons'+str(args['lons'][0])+'to'+str(args['lons'][1])+'_shelf.nc') 
        # convert to latent heat
        latent_heat_timeseries = summed * LH # kg/s * J/kg -> J/s -> W
        latent_heat_timeseries = latent_heat_timeseries * 10**-12 # W * 10**-12 -> TW
        
        print(latent_heat_timeseries)    
        
            
            
        # plot timeseries
        #fig = plt.figure(dpi=300)
        #plt.plot(latent_heat_timeseries.time, latent_heat_timeseries)
        #plt.title('latent heat for '+str(args['lons'][0])+' to '+str(args['lons'][1])+' E')
        #plt.ylabel('latent heat / TW')
        #plt.xlabel('Date')
        #fig.savefig('lh.png')
        
    #    if args['save_data']:
    #        data_file_name = 'ISMIPfw-data_for_latent_heat'
    #        latent_heat_timeseries.to_netcdf(data_file_name+'-SSP585FW.nc')
        
        # save results
        
    if 'input_regions' in  args['analysis']:
        example_tgrid = process_maui(xr.open_dataset('../../data/nemo_ci501o_1m_20150101-20150201_grid-T.nc')['sowflisf']).isel(time=0)
        areas = process_maui(xr.open_dataset('../../data/nemo_ci501o_1m_20150101-20150201_grid-T.nc'))['area']
        fig = plt.figure(dpi=300,figsize=(4,4))
        
        
        #projection = ccrs.NearsidePerspective(central_latitude=-90, satellite_height=3.6e6)
        #projection = ccrs.SouthPolarStereo()
        #panels = ['hist','ssp','maui','ssp_hist','maui_hist','maui_ssp']
        #panels = ['hist']
        #gs  = fig.add_gridspec(1,1)
        #axes = fig.add_subplot(gs[0], projection=projection)
        #axes.set_facecolor((0.8,0.8,0.8))
        #axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
        #axes.set_extent([-180, 180, -90, -65], ccrs.PlateCarree())
        #axes.add_feature(cfeature.COASTLINE)  
        
        example_tgrid = example_tgrid * areas
        basals = example_tgrid.where(example_tgrid<0, drop=True) * (360*24*60*60) * 10**-12 *-1
        print(np.nansum(basals))
        #axes.scatter(example_tgrid.lon, example_tgrid.lat, s=20*example_tgrid/np.nanmax(example_tgrid), transform=ccrs.PlateCarree())
        #axes.scatter(basals.lon, basals.lat, s=basals/np.nanmax(basals), marker='.', transform=ccrs.PlateCarree())
        #axes.scatter(basals.lon.where(~np.isnan(basals),drop=True), basals.lat.where(~np.isnan(basals),drop=True), s=0.5, c='r', marker='.', transform=ccrs.PlateCarree())
        plt.scatter(basals.lon, basals, c='k', marker='x')
        plt.ylabel(r'Basal melt / Gt.a$^{-1}$')
        plt.xlabel('Longitude / $^\circ$E')
        fig.savefig('input_regions.png')
    


    if 'opottemprmadvect' in args['analysis']:
        #months = 240
        
        # shelf mask
        dummy_file = process_maui(xr.open_dataset('../../data/nemo_ci501o_1m_20150101-20150201_grid-T.nc'))
        mask = shelf_mask(dummy_file)
        
        def apply_shelf_mask(ds,mask=None):
            if mask is None:
                ds = ds.where(shelf_mask(ds), drop=True)
            else:
                ds = ds.where(mask, drop=True)
            return ds
        transform_functions = [partial(process_maui, model='oce'),
                              partial(apply_shelf_mask,mask=mask),
                              partial(lat_select, lat_low=args['lats'][0], lat_high=args['lats'][1]),
                              partial(lon_select, lon_low=args['lons'][0], lon_high=args['lons'][1]),
                              partial(depth_slice, deep=args['depths'][1], shallow=args['depths'][0])]  
        
        #cmip = xr.open_dataset('../../data/cmip_wgets/opottemprmadvect/processed/cmip6-opottemprmadvect-1000m_AAcoast.nc')
        #cmip = xr.open_dataset('test.nc')
        #cmip = cmip.mean(dim='member_id').isel(time=slice(0,-1))
        cmip = xr.open_dataset('../../data/cmip_wgets/opottemprmadvect/processed/cmip6-opottemprmadvect-sum1000m_shelf.nc')
        for tf in transform_functions[2:-1]:
            cmip = tf(cmip)
 
        # var_at = '../../data/cmip_wgets/opottemprmadvect/data/opottemprmadvect_Emon_HadGEM3-GC31-LL_ssp585_r'
        # cmip = {}
        # for isuite in range(4):
        #     cmip_str = var_at+str(isuite+1)+'i1p1f3_gn_205001-209912.nc'
        #     cmip[str(isuite)] = xr.open_dataset(cmip_str)['opottemprmadvect'].isel(time=slice(0,599))
        #     cmip[str(isuite)] = cmip[str(isuite)].rename({'latitude':'lat', 'longitude':'lon', 'i':'x', 'j':'y'})
        #     lon = cmip[str(isuite)]['lon'].where(cmip[str(isuite)]['lon'] > 0, cmip[str(isuite)]['lon'] + 360)
        #     cmip[str(isuite)] = cmip[str(isuite)].assign_coords(lon=lon)
        #     for tf in transform_functions[1:]:
        #         cmip[str(isuite)] = tf(cmip[str(isuite)])
        # cmip_late = xr.concat([cmip[str(sn)] for sn in range(4)],
        #            pd.Index([str(sn) for sn in range(4)], name="member_id"))

        # cmip = {}
        # for isuite in range(4):
        #     cmip_str = var_at+str(isuite+1)+'i1p1f3_gn_201501-204912.nc'
        #     cmip[str(isuite)] = xr.open_dataset(cmip_str)['opottemprmadvect']
        #     cmip[str(isuite)] = cmip[str(isuite)].rename({'latitude':'lat', 'longitude':'lon', 'i':'x', 'j':'y'})
        #     lon = cmip[str(isuite)]['lon'].where(cmip[str(isuite)]['lon'] > 0, cmip[str(isuite)]['lon'] + 360)
        #     cmip[str(isuite)] = cmip[str(isuite)].assign_coords(lon=lon)
        #     for tf in transform_functions[1:]:
        #         cmip[str(isuite)] = tf(cmip[str(isuite)])
        # cmip_early = xr.concat([cmip[str(sn)] for sn in range(4)],
        #            pd.Index([str(sn) for sn in range(4)], name="member_id"))
        # cmip_ssp = xr.concat([cmip_early,cmip_late], dim='time')

        # get maui opottemprmadvec
        #maui = read_netcdfs(list_maui_files(suite='ci501'), drop_vars=vars2exclude(list_maui_files(suite='ci501')[0],'opottempadvect'), transform_func=transform_functions)['opottempadvect'] 
        suites = ['ci501','cm483','cn043','cn077']
        maui = {}
        for suite in suites:
            maui[suite] = read_netcdfs(list_maui_files(suite=suite), drop_vars=vars2exclude(list_maui_files(suite=suite)[0],'opottempadvect'), transform_func=transform_functions)['opottempadvect'] 
        maui = xr.concat([maui[sn] for sn in suites],
                   pd.Index([sn for sn in suites], name="suite_id"))
        maui = maui.mean(dim='suite_id').sum(dim='lev')
        areacello = lazy_load_cmip.get_hadgem_cmip6_experiment(experiment_id='piControl', variable_id='areacello', table_id='Ofx')
        for tf in transform_functions[1:4]:
            areacello = tf(areacello)
        areacello = areacello.isel(member_id=0)
       
        print(areacello) 
        print(maui)
        print(cmip)

        cmip = cmip*areacello
        maui=maui*areacello 
        
        cmip.to_netcdf('cmip-opot-spatial.nc')
        maui.to_netcdf('maui-opot-spatial.nc')
                
        
        
        # fig = plt.figure(dpi=300)
        # gs = fig.add_gridspec(1,1)
        # ax = fig.add_subplot(gs[0], projection=ccrs.NearsidePerspective(central_latitude=-90, satellite_height=3.6e6))
        # ax.add_feature(cfeature.COASTLINE)  
        # ax.pcolormesh(maui.lon,maui.lat,maui,transform=ccrs.PlateCarree())
        # fig.savefig('mask_test-opt.png')
        
        
        
        
       #print(areacello)
        
    #    SSP585 = cmip.sum(dim=('x','y')).rolling(time=12).mean() * 10**-12
    #    SSP585FW = maui.sum(dim=('x','y')).rolling(time=12).mean() * 10**-12 
        #print( (SSP585FW - SSP585).values )
        
        
     #   if args['save_data']:
      #      data_file_name = 'ISMIPfw-data_for_opottemprmadvect'
      #      SSP585.to_netcdf(data_file_name+'-SSP585.nc')
      #      SSP585FW.to_netcdf(data_file_name+'-SSP585FW.nc')
        
        
      #  fig = plt.figure()
      #  plt.plot(SSP585FW.time, SSP585FW, label='SSP585FW')
      #  plt.plot(SSP585.time, SSP585, label='SSP585')
        #plt.plot(SSP585.time, (SSP585FW - SSP585).values, label='difference')
       # plt.legend()
        #plt.plot(SSP585.time, (SSP585FW - SSP585).values, 'r:')
        
        #fig.savefig('opottemprmadvect.png')
        
    if 'make_opottemprmadvect' in args['analysis']:
        
        var_at = '../../data/cmip_wgets/opottemprmadvect/data/opottemprmadvect_Emon_HadGEM3-GC31-LL_ssp585_r'
        early_suffix = 'i1p1f3_gn_201501-204912.nc'
        late_suffix = 'i1p1f3_gn_205001-209912.nc'
        vn = 'opottemprmadvect'
        lons = (210,290) # longitudes for selecting AmBe
        
        def subset_3D(ds, lats=args['lats'], levs=args['depths']):
            ds = ds.rename({'i':'x','j':'y','latitude':'lat','longitude':'lon'})
            lon = ds['lon'].where(ds.lon>0, ds.lon+360)
            ds = ds.assign_coords(lon=lon)
            ds = depth_slice(ds, deep=args['depths'][1], shallow=args['depths'][0])
            ds = ds.sum(dim='lev')
            return ds
        
        cmip_members = {}
        for imember, member in enumerate(['1','2','3','4']):
            early = subset_3D(xr.open_dataset(var_at + member + early_suffix)[vn])
            late = subset_3D(xr.open_dataset(var_at + member + late_suffix)[vn])
            cmip_members[str(imember)] = xr.concat((early,late), dim='time')
            print('done member '+member)
            
            
        cmip = xr.concat( [cmip_members[str(sn)] for sn in range(4)],
                   pd.Index([str(sn) for sn in range(4)], name="member_id"))
        
        cmip = cmip.mean(dim='member_id')
        print(cmip)
        cmip.to_netcdf('cmip6-opottemprmadvect-sum1000m_global.nc')
      
        cmip_shelf = cmip.where(shelf_mask(cmip), drop=True)
        print(cmip_shelf)
        cmip_shelf.to_netcdf('cmip6-opottemprmadvect-sum1000m_shlef.nc')
        
    if 'opt_lh_plot' in args['analysis']:
        lons = (210, 290)
        #maui_lh = xr.open_dataarray('../../data/lh_spatial_sum1000m.nc').isel(member_id=0)
        #maui_ad = xr.open_dataarray('../../data/maui-opot-spatial.nc')
        #cmip_ad = xr.open_dataarray('../../data/cmip-opot-spatial.nc').isel(time=slice(0,-1))
        maui_lh = xr.open_dataarray('lh_spatial-sum1000.0_lons210.0to290.0_shelf.nc').isel(member_id=0)
        maui_ad = xr.open_dataarray('maui-opot-spatial.nc')
        cmip_ad = xr.open_dataarray('cmip-opot-spatial.nc').isel(time=slice(0,-1))

        lh_control = xr.open_dataset('../../data/cmip_wgets/ficeberg/ficeberg_Omon_HadGEM3-GC31-LL_ssp585_r3i1p1f3_gn_210001-210012.nc')['ficeberg']
        #lh_control = xr.open_dataset('../../data/ficeberg_Omon_HadGEM3-GC31-LL_ssp585_r3i1p1f3_gn_210001-210012.nc')['ficeberg'] # vowflisf is constant in ssp585 so we just need one value. ficeberg contains vowflisf in cells deeper than the surface
        lh_control = lh_control.rename({'latitude':'lat','longitude':'lon','j':'y','i':'x'})
        #print(lh_control)
        lh_control = lh_control.isel(time=0).where(shelf_mask(lh_control), drop=True)
        lon = lh_control['lon'].where(lh_control['lon'] > 0, lh_control['lon'] + 360)
        lh_control = lh_control.assign_coords(lon=lon)
        lh_control = depth_slice(lh_control, 1000, 2) # 2 here cuts off berg melt
        lh_control = lh_control.sum(dim='lev')
        lh_control = lon_select(lh_control, lons[0], lons[1])
        
        #area = process_maui(xr.open_dataset('../../data/nemo_ci501o_1m_20981201-20990101_grid-T.nc'))['area']
        area = process_maui(xr.open_dataset('../../data/nemo_ci501o_1m_20150101-20150201_grid-T.nc'))['area']
        area = area.where(shelf_mask(area), drop=True)
        area = lon_select(area, lons[0], lons[1])
        
        lh_control = lh_control*area
        lh_control = lh_control *333.55 * 10**3* 10**-12
        lh_control = lh_control.sum(dim=('x','y')).values
       
        print(lh_control)
        print(maui_lh)
        print(maui_ad)
        print(cmip_ad) 
        AmBe_lh = maui_lh
        #AmBe_lh = lon_select(maui_lh, lons[0], lons[1])
        AmBe_ad_maui = maui_ad
        AmBe_ad_cmip = cmip_ad
        
        plot_lh = AmBe_lh * 333.55 * 10**3 * 10**-12
        #plot_lh = AmBe_lh.sum(dim=('x','y'))*333.55 * 10**3 * 10**-12
        plot_ad_maui = AmBe_ad_maui.sum(dim=('x','y')).rolling(time=12).mean()*10**-12
        plot_ad_cmip = AmBe_ad_cmip.sum(dim=('x','y')).rolling(time=12).mean()*10**-12
        
        print(plot_lh.time[0])
        print(plot_lh.time[-1]) 
        plt.figure(dpi=300)
        #plt.title('latent heat and advective temperature tendencies for 220 to 290 E \n on continental shelf')
        plt.plot(plot_lh.time, -plot_lh,'k:', label='SSP585FW latent heat')
        plt.plot(plot_lh.time, np.ones_like(plot_lh)*lh_control, 'b:', label=r'SSP585 latent heat')
        
        plt.plot(plot_ad_maui.time,plot_ad_maui,'k', label=r'SSP585FW advective $\theta_\mathrm{o}$ tendency')
        plt.plot(plot_ad_cmip.time,plot_ad_cmip,'b', label=r'SSP585 advective $\theta_\mathrm{o}$ tendency')
        
        plt.plot(plot_ad_maui.time, plot_ad_maui - plot_ad_cmip,'r', label=r'SSP585FW-SSP585, advective $\theta_\mathrm{o}$ tendency')
        plt.plot(plot_lh.time, -plot_lh-lh_control, 'r:', label='SSP585FW-SSP585, latent heat')
        plt.xlabel('Date')
        plt.ylabel(r'latent heat or advcetive $\theta_\mathrm{o}$ tendency / TW')
        plt.legend(fontsize=8)
        plt.grid(True)
        
        plt.savefig('opottemprmadvect_latentHeat.png')
        
        
        
        
        
        
