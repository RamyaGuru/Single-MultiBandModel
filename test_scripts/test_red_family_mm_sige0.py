#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 21:17:52 2021

@author: ramyagurunathan

Test sigma_e0 mixing model Red Family
"""

import mixing_models_sigmae0 as mms
import sigmae0_alloy as sig
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import ternary


mat_props = {'formula': ['TaFeSb', 'VFeSb', 'NbFeSb'], 'BulkMod': [175E9, 167E9, 146E9],\
             'mpid': ['mp-631267', 'mp-567636', 'mp-9437'],\
             'AtmV' : [1.743E-29, 1.624E-29, 1.772E-29]}

#Add effective masses. Sources: Nb, V: A. Page J. of MAteriomics 2016, Ta: 
mat_props['eff_mass_p'] = [1.65, 2.5, 1.6] #Check the VFeSb value

mat_props['Nv'] = [{'p' : 4, 'n' : 3}, {'p' : 4, 'n' : 3}, {'p' : 4, 'n' : 3}]

#Source: PCCP 2016, 18 14107
mat_props['eff_mass_n'] = [0.36, 0.56, 0.38]
mat_props['eff_mass_p'] = [0.22, 0.39, 0.25]

'''
Get the sigmae0 data
'''

data_path = '/Users/ramyagurunathan/Documents/PhDProjects/Argonne_TECCA/SigmaE0_values/'
datafile = data_path + 'XFeSb_gauss_data_red_2nd_Gn.csv'
sigmae0_df = sig.fetch_sigmae0_dataframe(datafile)
#Drop multiphase samples
sigmae0_df = sigmae0_df.drop([8, 29, 36])
x_val, y_val, z_val = sig.fetch_endmember_values(sigmae0_df)

'''
Fit coefficients
'''
mat_props = sig.fit_sigma_e0_coefficient([x_val, y_val, z_val], mat_props, 'n', T = 300)

mat_props['sigma_endpts'] = [x_val, y_val, z_val]

'''
Split into binary data
'''
bin_list = mms.split_binary_data(sigmae0_df)

'''
Fit separate U values
'''
U_list = mms.fit_binary_U_list(bin_list, mat_props, d_key_list =  ['VFeSb', 'TaFeSb', 'TaFeSb'], d_h_pairs = [[1,2],[0,2],[0,1]])


'''
Print binary curves
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
plt.scatter(bin_list[2]['TaFeSb'], bin_list[2]['Sigma_E0'])
plt.ylabel(r'$\sigma_{e0}$')
plt.xlabel(r'x in Ta$_x$V$_{1-x}$FeSb')
plt.figure()
plt.plot(c_list, XZ_bin)
plt.scatter(bin_list[1]['TaFeSb'], bin_list[1]['Sigma_E0'])
plt.ylabel(r'$\sigma_{e0}$')
plt.xlabel(r'x in Ta$_x$Nb$_{1-x}$FeSb')
plt.figure()
plt.scatter(bin_list[0]['VFeSb'], bin_list[0]['Sigma_E0'])
plt.ylabel(r'$\sigma_{e0}$')
plt.xlabel(r'x in V$_x$Nb$_{1-x}$FeSb')
plt.plot(c_list, YZ_bin)


'''
Muggianu Result
'''
#Look at the endpoints
c_Ti = [0.01, 0.98]
sig_mu_Hf = mms.muggianu_model(sigmae0_df, bin_list, U_list, mat_props, 'n', c_Ti, T = 300)

'''
Plot ternary for muggianu result
'''
sig_data_dict = mms.convert_data_df_to_dict(sigmae0_df, 'TaFeSb', 'VFeSb')
sig_mu_dict = mms.run_sigmae0_muggianu_dict(sigmae0_df, bin_list, U_list, mat_props, 'n', T = 300, n = 100)
fig, ax = plt.subplots()
ax.axis("off")
figure, tax = ternary.figure(ax=ax, scale = 100)

tax.heatmap(sig_mu_dict, style = 'h', cmap=plt.cm.get_cmap('Spectral_r', 20),\
     cbarlabel=r'$\sigma_{e0}$', cb_kwargs = {'shrink' : 0.8, 'pad' : 0.01}, vmin = 25000, vmax = 125000,\
     scientific = False)

Tt = [(sigmae0_df['TaFeSb'][i] * 100, sigmae0_df['VFeSb'][i] * 100, sigmae0_df['NbFeSb'][i] * 100) for i in  list(sigmae0_df.index)]
At = sigmae0_df['Sigma_E0']
tax.scatter(Tt, c = list(At), colormap=plt.cm.get_cmap('Spectral_r', 20),\
     cmap=plt.cm.get_cmap('Spectral_r', 20), vmin = 25000, vmax = 125000,\
     scientific = False, s = 30, edgecolors = 'k', zorder = 10, clip_on = False)

tax.boundary(linewidth=2.0)

tax.top_corner_label('VFeSb')
tax.left_corner_label('NbFeSb', position = (0,0., 0))
tax.right_corner_label('TaFeSb', position = (0.95,0, 0))

tax.savefig('XFeSB_mm_sige0_ternary.pdf', bbox_inches = 'tight')



