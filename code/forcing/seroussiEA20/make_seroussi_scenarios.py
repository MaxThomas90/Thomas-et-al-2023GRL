#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 10:04:01 2021

@author: mthomas
"""

import xarray as xr
import numpy as np
import glob
import matplotlib.pyplot as plt


experiments = ['exp05', 'exp07', 'expA3', 'expA4', 'exp01', 'exp03']
descriptions = ['NorESM1, RCP8.5 (standard)', 'NorESM1, RCP2.6 (standard)', 'ISPL-CM5A-MR, RCP8.5', 'ISPL-CM5A-MR, RCP2.6', 'NorESM1, RCP8.5 (open)', 'NorESM, RCP2.6 (open)']
linestyles = ['r','b','r--','b--','r:','b:']
kgs2gty = 10**-12 * 365 * 24 * 60**2

data_home = '../../../data/seroussiEA20/ComputedScalarsPaper/'
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


## First replicate Seroussi et al Figure 12 for standard, medium ocean melt
fig = plt.figure(dpi=600)
r85 = exps['exp05']
xvar = r85['time']
mean_basal = np.zeros_like(r85['time'])
for mn in r85['models']:
    plt.plot(r85['data'][mn].time.values, - kgs2gty * r85['data'][mn].bmbfl.values, 'r-', alpha=0.1)
    for ts, year in enumerate(r85['data'][mn].time.values):
        if year<2100: # we need this if statement because some of the models include 2100.5
            mean_basal[ts] += r85['data'][mn].bmbfl.values[ts]

mean_basal_85 = - kgs2gty * mean_basal / len(r85['models'])
plt.plot(xvar, mean_basal_85, 'r')
plt.title('Basal melt, standard ocean sensitivity, no ice shelf fracture, NorESM1, RCP8.5 \n from Seroussi et al 2020 (Figure 12)')
plt.xlabel('Date')
plt.ylabel('Basal melt rate / Gt/a')

plt.savefig('seroussi_et_al_2020-Figure12_replica.png')
    

## similar to above but for open, medium ocean melt
# fig = plt.figure(dpi=600)
# r85 = exps['exp05']
# xvar = r85['time']
# mean_basal = np.zeros_like(r85['time'])
# for mn in r85['models']:
#     plt.plot(r85['data'][mn].time.values, - kgs2gty * r85['data'][mn].bmbfl.values, 'r-', alpha=0.1)
#     for ts, year in enumerate(r85['data'][mn].time.values):
#         if year<2100: # we need this if statement because some of the models include 2100.5
#             mean_basal[ts] += r85['data'][mn].bmbfl.values[ts]

# mean_basal_85 = - kgs2gty * mean_basal / len(r85['models'])
# plt.plot(xvar, mean_basal_85, 'r')
# plt.title('Basal melt, standard ocean sensitivity, no ice shelf fracture, NorESM1, RCP8.5 \n from Seroussi et al 2020 (Figure 12)')
# plt.xlabel('Date')
# plt.ylabel('Basal melt rate / Gt/a')
    
## Then do the same for RCP2.6
fig = plt.figure(dpi=600)
r26 = exps['exp07']
xvar = r26['time']
mean_basal = np.zeros_like(r26['time'])
for mn in r26['models']:
    plt.plot(r26['data'][mn].time.values, - kgs2gty * r26['data'][mn].bmbfl.values, 'b-', alpha=0.1)
    for ts, year in enumerate(r26['data'][mn].time.values):
        if year<2100:
            mean_basal[ts] += r26['data'][mn].bmbfl.values[ts]

mean_basal_26 = - kgs2gty * mean_basal / len(r26['models'])
plt.plot(xvar, mean_basal_26, 'b')
plt.title('Basal melt, standard ocean sensitivity, no ice shelf fracture, NorESM1, RCP2.6 \n from Seroussi et al 2020')
plt.xlabel('Date')
plt.ylabel('Basal melt rate / Gt/a')
plt.savefig('seroussi_et_al_2020-Figure12_with_RCP26.png')
    

## For context, we can look at the mean basal for all the runs with standard medium melt
fig = plt.figure(dpi=600)
for iexp, exp in enumerate(experiments):
    mean_basal = np.zeros_like(exps['exp05']['time'])
    print(exp)
    #min_basal = np.zeros_like(exps['exp05']['time'])
    #max_basal = np.zeros_like(exps['exp05']['time'])
    for mn in exps[exp]['models']:
        for ts, year in enumerate(exps['exp05']['time']):
            if year<2100:
                mean_basal[ts] += exps[exp]['data'][mn].bmbfl.values[ts]
    mean_basal = - kgs2gty * mean_basal / len(exps[exp]['models'])
    plt.plot(exps['exp05']['time'], mean_basal, linestyles[iexp], label=descriptions[iexp])
        
plt.legend()
plt.title('Basal melt, no ice shelf fracture \n from Seroussi et al 2020')
plt.xlabel('Date')
plt.ylabel('Basal melt rate / Gt/a')
plt.savefig('seroussi_et_al_2020-all_runs_with_standard_melt.png')


## Now we define functions to capture the melt
trel = xvar - xvar[0]
# total freshwater input
# integral_85 = sum(mean_basal_85) # each timestep is one year, so the integral is just the sum
# integral_26 = sum(mean_basal_26)

# try a few degrees of polynomials
plt.figure(dpi=600)
plt.plot(trel, mean_basal_85, 'r-', label='NorESM1, RCP8.5')#', total = ' + str(int(integral_85)) + ' Gt')
max_degree = 10
coeffs_85 = {}
coeffs_26 = {}
for coeffs in range(2, max_degree):
    coeffsp1 = str(coeffs+1) # this is the number of coefficients in the polynomial
    coeffs_85['coeffs'+coeffsp1] = np.polyfit(trel, mean_basal_85, coeffs)
    plt.plot(trel, np.poly1d(np.polyfit(trel,mean_basal_85,coeffs))(trel), alpha=0.7, label=coeffsp1 + ' coefficients')# + ', total = ' + str(int(sum(np.poly1d(np.polyfit(trel,mean_basal_85,deg))(trel)))) + ' Gt')
plt.legend()
plt.title('Polynomial representations of NorESM1 RCP8.5 basal melt')
plt.xlabel('Years elapsed')
plt.ylabel('Basal melt rate / Gt/a')
plt.savefig('polynomials_for_RCP85.png')


# try a few degrees of polynomials
plt.figure(dpi=600)
plt.plot(trel, mean_basal_26, 'b-', label='NorESM1, RCP2.6')#', total = ' + str(int(integral_26)) + ' Gt')
for coeffs in range(2, max_degree):
    coeffsp1 = str(coeffs+1) # this is the true degree of the polynomial
    coeffs_26['coeffs'+coeffsp1] = np.polyfit(trel, mean_basal_26, coeffs)
    plt.plot(trel, np.poly1d(np.polyfit(trel,mean_basal_26,coeffs))(trel), alpha=0.7, label=coeffsp1 + ' coefficients')# + ', total = ' + str(int(sum(np.poly1d(np.polyfit(trel,mean_basal_26,deg))(trel)))) + ' Gt')
plt.legend()
plt.title('Polynomial representations of NorESM1 RCP2.6 basal melt')
plt.xlabel('Years elapsed')
plt.ylabel('Basal melt rate / Gt/a')
plt.ylim((0,1000))
plt.savefig('polynomials_for_RCP26.png')

# Noting that all polynomials capture the total melt well, I pick the polynomial with 8 terms (including offset)
chosen_coeffs = {}
chosen_coeffs['RCP85'] = coeffs_85['coeffs8']
chosen_coeffs['RCP26'] = coeffs_26['coeffs8']

plt.figure(dpi=600)
plt.title('Proposed basal melt parameterisations for HadGEM')
plt.plot(xvar, mean_basal_85, 'r-', label='NorESM1, RCP8.5, ISMIP')
#plt.plot(xvar, np.poly1d(np.polyfit(trel,mean_basal_85,7))(trel), 'r-.', label='Parameterised, degree '+str(7))
plt.plot(xvar, np.poly1d(chosen_coeffs['RCP85'])(trel), 'r-.', label='Parameterised, degree '+str(8))

plt.plot(xvar, mean_basal_26, 'b-', label='NorESM1, RCP2.6, ISMIP')
plt.plot(xvar, np.poly1d(chosen_coeffs['RCP26'])(trel), 'b-.', label='Parameterised, degree '+str(8))
plt.xlabel('Date')
plt.ylabel('Basal melt rate / Gt/a')
plt.legend()
plt.savefig('chosen_basal_parameterisation-RCP85.png')

# And converting to a basal + calving flux
plt.figure(dpi=600)
plt.title('Proposed basal + calving flux parameterisations for HadGEM')
plt.plot(xvar, mean_basal_85 / 0.55, 'r-', label='NorESM1, RCP8.5, ISMIP, /0.55')
plt.plot(xvar, np.poly1d(chosen_coeffs['RCP85'])(trel) / 0.55, 'r-.', label='Parameterised, degree '+str(7))

plt.plot(xvar, mean_basal_26 / 0.55, 'b-', label='NorESM1, RCP2.6, ISMIP, /0.55')
plt.plot(xvar, np.poly1d(chosen_coeffs['RCP26'])(trel) / 0.55, 'b-.', label='Parameterised, degree '+str(7))
plt.xlabel('Date')
plt.ylabel('Basal melt + calving flux / Gt/a')
plt.legend()
plt.savefig('chosen_basal_plus_calving_parameterisation-Gta.png')

# For ease of model input, deriving coefficients for kg/s and already scaled
targets = mean_basal_85 / (0.55*kgs2gty) # mean basal melt, scaled to calving flux (HadGEM standard), units to kg/s (as in HadGEM)
years = trel
xs_85 = np.polyfit(trel, targets, 7) # coefficients for polynomial for RCP8.5
plt.figure(dpi=600)
plt.plot(years + 2015, targets, 'r-', label='NorESM1, RCP8.5, ISMIP, /0.55')
plt.plot(years + 2015, np.poly1d(xs_85)(years), 'r-.', label='Parameterised, degree '+str(7))
plt.ylabel('Basal + calving / kg/s')
plt.xlabel('Date')
plt.title('Proposed basal + calving flux parameterisations for HadGEM \n')

# For ease of model input, deriving coefficients for kg/s and already scaled
targets = mean_basal_26 / (0.55*kgs2gty)
xs_26 = np.polyfit(trel, targets, 7)
#plt.figure(dpi=600)
plt.plot(years + 2015, targets, 'b-', label='NorESM1, RCP2.6, ISMIP, /0.55')
plt.plot(years + 2015, np.poly1d(xs_26)(years), 'b-.', label='Parameterised, degree '+str(7))
plt.legend()
plt.grid('major')
plt.tight_layout()
plt.savefig('chosen_basal_plus_calving_parameterisation-kgs.png')

with open('coefficients.txt', "w") as coefficients:
    coefficients.write('Coefficients for kg/s total HadGEM freshwater input for RCP8.5 are:\n')
    coefficients.write(str(xs_85))
    coefficients.write('\n\n')
    coefficients.write('Coefficients for kg/s total HadGEM freshwater input for RCP2.6 are:\n')
    coefficients.write(str(xs_26))
















