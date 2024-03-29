#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 15:28:06 2021

@author: ramyagurunathan

TaCoSn Black Family
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
datafile = data_path + 'XCoSn_gauss_data_black.csv'
sigmae0_df = sig.fetch_sigmae0_dataframe(datafile)
x_val, y_val, z_val = sig.fetch_endmember_values(sigmae0_df)
'''
Initalize Endmember Properties
Order of compounds: TaCoSn, NbCoSn, VCoSn

Average Effective Mass values: Obtained from the MPContribs database
Dividing the carr. conc., eff. mass, 
'''
#Bulk Modulus sources Ta: Materials Project
mat_props = {'formula': ['TaCoSn', 'NbCoSn', 'VCoSn'], 'BulkMod': [163E9, 102E9, 201E9],\
             'AtmV' : [( 5.948E-10)**3 / 12, (5.946E-10)**3 / 12, 1.745E-29]}

mat_props['sigma_endpts'] = [x_val, y_val, z_val]

#Add effective masses. Sources: Nb, V: A. Page J. of MAteriomics 2016, Ta: 
mat_props['eff_mass_p'] = [1.65, 5.1, 1.6] #Check the VFeSb value

mat_props['Nv'] = [{'p' : 10, 'n' : 3}, {'p' : 10, 'n' : 3}, {'p' : 10, 'n' : 3}] #P may be different

#Add the coefficients
#mat_props = sig.fit_sigma_e0_coefficient([x_val, y_val, z_val], mat_props, 'p', T = 300)

'''
Split into binary data
'''
bin_list = mms.split_binary_data(sigmae0_df)

'''
Fit separate U values
'''
U_list = mms.fit_binary_U_list(bin_list, mat_props, d_key_list =  ['NbCoSn', 'TaCoSn', 'TaCoSn'], d_h_pairs = [[1,2],[0,2],[0,1]])

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
plt.scatter(bin_list[2]['TaCoSn'], bin_list[2]['Sigma_E0'])
plt.ylabel(r'$\sigma_{e0}$')
plt.xlabel(r'x in Ta$_x$V$_{1-x}$CoSn')
plt.figure()
plt.plot(c_list, XZ_bin)
plt.scatter(bin_list[1]['TaCoSn'], bin_list[1]['Sigma_E0'])
plt.ylabel(r'$\sigma_{e0}$')
plt.xlabel(r'x in Ta$_x$Nb$_{1-x}$CoSn')
plt.figure()
plt.scatter(bin_list[0]['NbCoSn'], bin_list[0]['Sigma_E0'])
plt.ylabel(r'$\sigma_{e0}$')
plt.xlabel(r'x in Nb$_x$V$_{1-x}$CoSn')
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

Tt = [(sigmae0_df['TaCoSn'][i] * 100, sigmae0_df['NbCoSn'][i] * 100, sigmae0_df['VCoSn'][i] * 100) for i in  list(sigmae0_df.index)]
At = sigmae0_df['Sigma_E0']
tax.scatter(Tt, c = list(At), colormap=plt.cm.get_cmap('Spectral_r', 20),\
     cmap=plt.cm.get_cmap('Spectral_r', 20),vmin = 20000, vmax = 90000,\
     scientific = False, s = 30, edgecolors = 'k', zorder = 10, clip_on = False)

tax.boundary(linewidth=2.0)

tax.top_corner_label('NbCoSn')
tax.left_corner_label('VCoSn', position = (0,0.04, 0))
tax.right_corner_label('TaCoSn', position = (0.95,0.04, 0))

tax.savefig('XCoSm_mm_sige0_ternary.pdf', bbox_inches = 'tight')


fig, ax = plt.subplots()
ax.axis("off")
figure, tax = ternary.figure(ax=ax, scale = 100)

tax.heatmap(rho_excess_dict, style = 'h', cmap=plt.cm.get_cmap('Spectral_r', 20),\
     cbarlabel=r'excess $\rho_{e0}$',\
     scientific = False)

tax.boundary(linewidth=2.0)

tax.top_corner_label('NbCoSn')
tax.left_corner_label('VCoSn', position = (0,0.04, 0))
tax.right_corner_label('TaCoSn', position = (0.95,0.04, 0))


tax.savefig('XCoSm_mm_exc_rho_ternary.pdf', bbox_inches = 'tight')
