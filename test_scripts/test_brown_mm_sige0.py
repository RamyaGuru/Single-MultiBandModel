#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 18:03:19 2021

@author: ramyagurunathan

Test brown family
"""

import sys
sys.path.append('../electronic_models/')
import mixing_models_sigmae0 as mms
import sigmae0_alloy as sig
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import ternary
import csv


'''
constants
'''
e = 1.602e-19
m_e = 9.109e-31
hbar = 1.045e-34
kb = 1.38e-23

data_path = '/Users/ramyagurunathan/Documents/PhDProjects/Argonne_TECCA/SigmaE0_values/'
datafile = data_path + 'XCoSb_gauss_data_brown.csv'
sigmae0_df = sig.fetch_sigmae0_dataframe(datafile)
x_val, y_val, z_val = sig.fetch_endmember_values(sigmae0_df)
'''
Initalize Endmember Properties
Order of compounds: TaCoSn, NbCoSn, VCoSn

Average Effective Mass values: Obtained from the MPContribs database
Dividing the carr. conc., eff. mass, 
'''
#Bulk Modulus sources Ta: Materials Project
mat_props = {'formula': ['TiCoSb', 'ZrCoSb', 'HfCoSb'], 'BulkMod': [147E9, 139E9, 144E9],\
             'AtmV' : [(5.885E-10)**3 / 12, (6.0697E-10)**3 / 12, (6.0390E-10)**3 / 12]}

mat_props['sigma_endpts'] = [x_val, y_val, z_val]


#Add the coefficients
#mat_props = sig.fit_sigma_e0_coefficient([x_val, y_val, z_val], mat_props, 'p', T = 300)

'''
Split into binary data
'''
bin_list = mms.split_binary_data(sigmae0_df)

'''
Fit separate U values
'''
U_list = mms.fit_binary_U_list(bin_list, mat_props, d_key_list =  ['ZrCoSb', 'TiCoSb', 'TiCoSb'], d_h_pairs = [[1,2],[0,2],[0,1]])

'''
Print binaries
'''
first = 1e-10
last = 9.9999999999999e-1
j = 0
XY_bin = []
XZ_bin = []
YZ_bin = []
for c in np.arange(first,1, (last - first) / 100):
    XY_bin.append(sig.binary_sigma_e0_endpts(c, U_list[2], 0, 1, mat_props, 'n', 300))
    XZ_bin.append(sig.binary_sigma_e0_endpts(c, U_list[1], 0, 2, mat_props, 'n', 300))
    YZ_bin.append(sig.binary_sigma_e0_endpts(c, U_list[0], 1, 2, mat_props, 'n', 300)) 
c_list = np.arange(first,1, (last - first) / 100)
plt.figure()
plt.plot(c_list, XY_bin)
plt.scatter(bin_list[2]['TiCoSb'], bin_list[2]['Sigma_E0'])
plt.ylabel(r'$\sigma_{e0}$')
plt.xlabel(r'x in Ti$_x$Zr$_{1-x}$CoSb')
plt.figure()
plt.plot(c_list, XZ_bin)
plt.scatter(bin_list[1]['TiCoSb'], bin_list[1]['Sigma_E0'])
plt.ylabel(r'$\sigma_{e0}$')
plt.xlabel(r'x in Ti$_x$Hf$_{1-x}$CoSb')
plt.figure()
plt.scatter(bin_list[0]['ZrCoSb'], bin_list[0]['Sigma_E0'])
plt.ylabel(r'$\sigma_{e0}$')
plt.xlabel(r'x in Zr$_x$Hf$_{1-x}$CoSb')
plt.plot(c_list, YZ_bin)

'''
Print ternary
'''
#sig_data_dict = mms.convert_data_df_to_dict(sigmae0_df, 'TaCoSn', 'VCoSn')
rho_excess_dict, sig_mu_dict = mms.run_sigmae0_muggianu_dict(sigmae0_df, bin_list, U_list, mat_props, 'p', T = 300, n = 200)

fig, ax = plt.subplots()
ax.axis("off")
figure, tax = ternary.figure(ax=ax, scale = 100)

tax.heatmap(sig_mu_dict, style = 'h', cmap=plt.cm.get_cmap('Spectral_r', 20),\
     cbarlabel=r'$\sigma_{e0}$',vmin = 20000, vmax = 90000,\
     scientific = False)

Tt = [(sigmae0_df['TiCoSb'][i] * 100, sigmae0_df['ZrCoSb'][i] * 100, sigmae0_df['HfCoSb'][i] * 100) for i in  list(sigmae0_df.index)]
At = sigmae0_df['Sigma_E0']
tax.scatter(Tt, c = list(At), colormap=plt.cm.get_cmap('Spectral_r', 20),\
     cmap=plt.cm.get_cmap('Spectral_r', 20),vmin = 20000, vmax = 90000,\
     scientific = False, s = 30, edgecolors = 'k', zorder = 10, clip_on = False)

tax.boundary(linewidth=2.0)

tax.top_corner_label('ZrCoSb')
tax.left_corner_label('HfCoSb', position = (0,0.04, 0))
tax.right_corner_label('TiCoSb', position = (0.95,0.04, 0))

tax.savefig('XCoSb_mm_sige0_ternary.pdf', bbox_inches = 'tight')


fig, ax = plt.subplots()
ax.axis("off")
figure, tax = ternary.figure(ax=ax, scale = 100)

tax.heatmap(rho_excess_dict, style = 'h', cmap=plt.cm.get_cmap('Spectral_r', 20),\
     cbarlabel=r'excess $\rho_{e0}$',\
     scientific = False)

tax.boundary(linewidth=2.0)

tax.top_corner_label('ZrCoSb')
tax.left_corner_label('HfCoSb', position = (0,0.04, 0))
tax.right_corner_label('TiCoSb', position = (0.95,0.04, 0))


tax.savefig('XCoSb_mm_exc_rho_ternary.pdf', bbox_inches = 'tight')


