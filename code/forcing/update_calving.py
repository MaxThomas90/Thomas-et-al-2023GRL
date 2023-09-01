#' % Update NEMO Antarctic calving distribution
#' % Max Thomas

#' # Summary
#' Here we change the values in the forcing file:
#'      ecalving_v2.2x.nc
#' taken from /nesi/project/niwa00013/tids/OCEAN/hadgem3/forcing/ocean/eORCA1v2.2x/
#' This script is copied from update_FW.py, which works on the runoff ancil



#' # Dependencies

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature
import os

home_dir = '../..'
#data_dir = home_dir + '/data/forcing/'
data_dir = home_dir + '/configurations/forcing/calving_ancils/'
docs_dir = home_dir + '/documents/notebooks/suite_info/'

#' ## Region definition
#' The first thing we do is load the region definitions to a dataframe
regions = pd.read_csv(docs_dir + 'region_selectors.csv', delimiter='\s+', index_col=0)
print('----------------------------------------------------')
print('Updating calving for ' + os.sys.argv[1])
#' Next we find our region using the argument passed to the script
region = regions.loc[os.sys.argv[1]]
print(region)

#' ## Fuctions
#' ### Antarctica centered plotting

# sets up plots centered on Antarctica 
def AA_axes(subplot_index, title, lonlon_latlat=[-180,180,-90,-60]):
    ax_out = plt.subplot(2,3,subplot_index, projection=ccrs.SouthPolarStereo())
    ax_out.set_extent([-180,180,-90,-60], ccrs.PlateCarree())
    ax_out.add_feature(cartopy.feature.LAND, edgecolor='k', facecolor='w')
    #ax1.add_feature(cartopy.feature.OCEAN)
    ax_out.gridlines()
    ax_out.set_title(title)
    return ax_out

#' ### rescale intensive
#' we first convert intensive to extensive by multiplying through by area
def intensive2extensive(intensive_values, areas, rescaler = 1.):
    extensive_values = intensive_values.copy()
    extensive_values = np.multiply(extensive_values, areas)
    extensive_values = np.multiply(extensive_values, rescaler)
    if not( type(rescaler) is float):
        print('Rescaling extensive values')
        print('Sum of rescaled extensive = ' + str(km3yr2Gtyr(extensive_values.sum())) + ' Gt/yr')
    return extensive_values

#' ### get back intensive
#' recover intensive values, needed as forcing, from extensive values
def extensive2intensive(extensive_values, areas):
    intensive_values = np.divide(extensive_values, areas)
    return intensive_values

#' ### convert kg/s to Gt/yr
def km3yr2Gtyr(kgs):
   # gtyr = ( 1/( 1000 * 10**9 ) ) * 365*24*60*60 * kgs # Converts kg/s to Gt/yr
    gtyr = 1 * kgs
    return gtyr

#' ### calculate value of rescaler
#' first we sum the extensive fluxes for areas inside logical selection - F_changed
#' next we sum the extensive fluxes for areas outside the logical selection - F_unchanged
#' next we calculate the rescaler as: rescaler = ( F_scaled - F_unchanged ) / F_changed
def calculate_scaler(intensive_values, areas, logical_selector, F_scaled, selection_name='this selection'):
    vex = intensive2extensive(intensive_values, areas, 1.)
    F_changed = km3yr2Gtyr( vex[logical_selector].sum() )
    F_unchanged = km3yr2Gtyr( vex[np.invert(logical_selector)].sum() )
    scaler = (F_scaled - F_unchanged) / F_changed
    print('scaler = ' + str(round(scaler,3)) + ' for ' + selection_name + ' and F_scaled = ' + str(F_scaled) + ' Gt/yr')
    return scaler

#' ### Main function
#' Here, we update the freshwater inputs, save a new .nc, and plot the results as a sanity check.

# updates freshwater inputs, saves new .nc, plots results
print('----------------------------------------------------')
def update_FWF(filename='ecalving_v2.2x.nc',
               longitude_selection_limits=[-150,160], # less than, more than
               longitude_selection_method='or', # or, and, not
               selection_name='Ross_TEST',
               month=0, 
               F_scaled=3150.):
    
    data_in = xr.open_dataset(data_dir + filename)
    zvars_in = data_in.calvingmask[month,:,:]
    #FW_data = data_in[month,:,:]
    lsl = longitude_selection_limits
    lsm = longitude_selection_method

    areas = xr.open_dataset('areas.nc').area.values # This is just some NEMO output where we can get the grid cell areas
    # dirty patch as values already extensive...
    areas = 1.
    ls_lat = data_in.nav_lat.values<-50 # cut off northern hemisphere
    
    # Here, we make a boolean matrix that picks out the required longitude
    if lsm=='or': # useful when we span 0 degrees (Ross)
        ls_lon = np.add(data_in.nav_lon.values<lsl[0], 
                        data_in.nav_lon.values>lsl[1])
    elif lsm=='and': # useful otherwise
        ls_lon = np.multiply(data_in.nav_lon.values>lsl[0], 
                             data_in.nav_lon.values<lsl[1])
    elif lsm=='not': # Untested
        ls_lon = np.invert(np.add(data_in.nav_lon.values<lsl[0], 
                        data_in.nav_lon.values>lsl[1]))
    
    ls_lat_lon = np.multiply(ls_lat,ls_lon) # this selector picks out our desired data
    
    # calculate required value to rescale extensive fluxes
    # this number will preserve F_scaled when the extensive fluxes for the given region are multiplied by it
    scaler = calculate_scaler(zvars_in.values,
                              areas,
                              ls_lat_lon,
                              F_scaled,
                              selection_name)
    
    #rescaler = ls_lat_lon.astype(float) * scaling_factor # this is a matrix. Where we want to change the value, the value is the rescaler. Else, the value is 1
    rescaler = ls_lat_lon.astype(float) * scaler # this is a matrix that will recover F_scaled when multiplying through the extensive fluxes
    rescaler[rescaler==0] = 1.
    
    # this variable will be updated and written to the .nc
    zvars_out = zvars_in.copy()
    
    # rescale extensive, then recover intensive
    zvars_out.values = extensive2intensive(
                            intensive2extensive(zvars_in.values, areas, rescaler)
                            , areas)
    
    # sanity check
    print('Sum of prescribed extensive = ' + str(F_scaled) + ' Gt/yr')
    
    # plotting
    x = data_in.nav_lon.values
    y = data_in.nav_lat.values
    
    # Intensive fluxes for plotting (kg/s/m2)
    intensive_0 = zvars_in.values       
    intensive_new = zvars_out.values   
    intensive_diff = np.subtract(zvars_out.values, zvars_in.values)
    
    # extensive values for plotting (Gt/yr)
    extensive_0 = km3yr2Gtyr(intensive2extensive(zvars_in.values,areas,rescaler=1.))
    extensive_new = km3yr2Gtyr(intensive2extensive(zvars_out.values,areas,rescaler=1.))
    extensive_diff = np.subtract(extensive_new, extensive_0)
    
    fig = plt.figure(figsize=(20,10))
    
    ax1 = AA_axes(1, 'Freshwater input, default, extensive km3/yr?')
    im1 = ax1.scatter(x,y,c=intensive_0,s=intensive_0,cmap='Blues',transform=ccrs.PlateCarree())

    ax2 = AA_axes(2, 'Freshwater input, '+selection_name+', F_scaled = '+str(F_scaled) + ' Gt/yr, extensive km3/yr?')
    im2 = ax2.scatter(x,y,c=intensive_new,s=intensive_new,cmap='Blues',transform=ccrs.PlateCarree())
    
    ax3 = AA_axes(3, '$F_\mathrm{FW, 0} - F_\mathrm{FW, new}$, extensive km3/yr?')
    im3 = ax3.scatter(x,y,c=intensive_diff,s=intensive_diff,cmap='Blues',transform=ccrs.PlateCarree())
    
    ax4 = AA_axes(4, 'Freshwater input, default, extensive')
    im4 = ax4.scatter(x, y, c=extensive_0/extensive_0.max(), s=100*extensive_0/extensive_0.max(), cmap='Blues', transform=ccrs.PlateCarree())
    print('Default $F_\mathrm{FW} =$ ' + str(extensive_0.sum().round(0)) + 'Gt/yr')
    
    ax5 = AA_axes(5, 'Freshwater input, '+selection_name+', F_scaled = '+str(F_scaled) + ' Gt/yr extensive')
    im5 = ax5.scatter(x, y, c=extensive_new/extensive_0.max(), s=100*extensive_new/extensive_0.max(), cmap='Blues', transform=ccrs.PlateCarree())
    print('New $F_\mathrm{FW} =$ ' + str(extensive_new.sum().round(0)) + 'Gt/yr')
    
    ax6 = AA_axes(6, '$F_\mathrm{FW, 0} - F_\mathrm{FW, new}$, extensive')
    im6 = ax6.scatter(x, y, c=extensive_diff/extensive_0.max(), s=100*extensive_diff/extensive_0.max(), cmap='Blues', transform=ccrs.PlateCarree())
    print('$F_\mathrm{FW, 0} - F_\mathrm{FW, new}=$ ' + str(extensive_diff.sum().round(0)) + 'Gt/yr')
    
   # if F_scaled==2494:
   #     fn_out = 'eORCA1_runoff_v2.2x_' + selection_name + 'double'
   # elif F_scaled==4988:
   #     fn_out = 'eORCA1_runoff_v2.2x_' + selection_name + 'quadruple'
   # else:
   #     fn_out = 'eORCA1_runoff_v2.2x_' + selection_name + str(int(F_scaled))
    fn_out = 'ecalving_v2.2x_' + selection_name    

    data_out = data_in.copy()
    for mn in range(12):
        data_out.calvingmask.values[mn,:,:] = zvars_out.values
     
    data_out.to_netcdf(data_dir + fn_out + '.nc')
    plt.show()
    fig.savefig(data_dir + fn_out + '.png')
    
update_FWF(longitude_selection_limits=[region.longitude1,region.longitude2],
           longitude_selection_method=region.longitude_method,
           selection_name=os.sys.argv[1],
           F_scaled=region.F_scaled_calv)
print('----------------------------------------------------')




