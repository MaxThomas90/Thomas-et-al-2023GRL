#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 15:27:12 2021

@author: maxthomas
"""

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly

print('Extracting timeseries of Antarctic SLR contribution from Edwards et al. (2021)')
if 'risk' in sys.argv[1].lower():
    data_dir = 'risk_averse/'
    print('Risk averse (pessemistic) projections')
elif 'main' in sys.argv[1].lower():
    data_dir = 'main/'
    print('Main projections')

if sys.argv[2].lower()=='slr':
    units, ylab1, ylab2 = 1., ' / cmSLR', ' / cmSLR/a' 
elif sys.argv[2].lower()=='gt':
    units, ylab1, ylab2 = 3600., ' / Gt', ' / Gt/a'
elif sys.argv[2].lower()=='sv':
    units, ylab1, ylab2 = 0.114, ' / Sv.a', ' / Sv'
#cmSLR_to_Gt = 3600
#cmSLR_to_Sva = 0.114
SSP = ['119', '126', '245', '370', '585', 'NDC']
files = []
for ssp in SSP:
    files.append('projections_FAIR_SSP' + ssp + '.csv')

fig1 = plt.figure(dpi=300, figsize=(7,7))
gs = fig1.add_gridspec(2,1)    
ax1 = fig1.add_subplot(gs[0])
#ax2 = ax1.twinx()
AA_ALL = {'mean': []}
for ssp in SSP:
    AA_ALL[ssp] = []
for ifn, fn in enumerate(files):
    data = pd.read_csv(data_dir + fn)
    WAIS = data.loc[lambda data: data['region'] == 'WAIS', :]
    EAIS = data.loc[lambda data: data['region'] == 'EAIS', :]
    PENi = data.loc[lambda data: data['region'] == 'PEN',  :]
    years = []
    for yr in np.unique(data.year):
        aa_year = np.mean(WAIS.loc[lambda WAIS: WAIS['year'] == yr, 'SLE'])
        aa_year += np.mean(EAIS.loc[lambda EAIS: EAIS['year'] == yr, 'SLE'])
        aa_year += np.mean(PENi.loc[lambda PENi: PENi['year'] == yr, 'SLE'])
        AA_ALL[SSP[ifn]] = np.append(AA_ALL[SSP[ifn]], aa_year)
        years = np.append(years, yr)
    ax1.plot(years, AA_ALL[SSP[ifn]] * units, label='ssp'+SSP[ifn])
aa_mean = 0
for ssp in SSP:
    aa_mean += AA_ALL[ssp]
AA_ALL['mean'] = aa_mean / len(SSP)
ax1.plot(years, AA_ALL['mean'] * units, 'r', label='mean of scenarios')

ps = poly.polyfit(years, AA_ALL['mean'], 5)
fit_AA = poly.polyval(years, ps)*units
ax1.plot(years, fit_AA, 'k:', label='fit through mean')

ax1.legend()
#ax1.set_xlabel('Year')
ax1.set_ylabel('Total AA contribution'+ylab1)
#plt.show()
#fig1.savefig('cumulative_AA_mass_loss_'+sys.argv[1]+'.png')


ax2 = fig1.add_subplot(gs[1])
#fig2 = plt.figure()
AA_rate = {'mean': []}
dfit_dt = np.gradient(fit_AA)
for ssp in SSP:
    AA_rate[ssp] = np.diff(AA_ALL[ssp])
    ax2.plot(years[1:], AA_rate[ssp] * units, label='ssp'+ssp)
ax2.plot(years[1:], np.diff(AA_ALL['mean'])*units, 'r', label='mean')
ax2.plot(years, dfit_dt, 'k:', label='d/dt fit through mean')
#ax2.legend()
ax2.set_xlabel('Year')
ax2.set_ylabel('Rate of AA contribution' + ylab2)

gs.tight_layout(fig1)
plt.show()

fig1.savefig('AA_mass_loss_'+sys.argv[1]+'_'+sys.argv[2].lower()+'.png')




