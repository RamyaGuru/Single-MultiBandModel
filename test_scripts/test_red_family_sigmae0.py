#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 13:47:18 2021

@author: ramyagurunathan


Red Family Sigma_e0 Test
"""
import sys
sys.path.append('../electronic_models/')
import sigmae0_alloy as sig
import ternary
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv

'''
constants
'''
e = 1.602e-19
m_e = 9.109e-31
hbar = 1.045e-34
kb = 1.38e-23


data_path = '/Users/ramyagurunathan/Documents/PhDProjects/Argonne_TECCA/SigmaE0_values/'
datafile = data_path + 'XFeSb_gauss_data_red_2nd_Gn.csv'
sigmae0_df = sig.fetch_sigmae0_dataframe(datafile)
#Drop multiphase samples
sigmae0_df = sigmae0_df.drop([8, 29, 36])
x_val, y_val, z_val = sig.fetch_endmember_values(sigmae0_df)
'''
Initalize Endmember Properties
Order of compounds: TaFeSb, VFeSb, NbFeSb

Average Effective Mass values: Obtained from the MPContribs database
Dividing the carr. conc., eff. mass, 
'''
#Bulk Modulus sources Ta: Naydenov J Phys Mater, V: Bo Kong Physica B, Nb: Silpawilawan 
mat_props = {'formula': ['TaFeSb', 'VFeSb', 'NbFeSb'], 'BulkMod': [175E9, 167E9, 146E9],\
             'mpid': ['mp-631267', 'mp-567636', 'mp-9437'],\
             'AtmV' : [1.743E-29, 1.624E-29, 1.772E-29]}



#Add effective masses. Sources: Nb, V: A. Page J. of MAteriomics 2016, Ta: 
mat_props['eff_mass_p'] = [1.65, 2.5, 1.6] #Check the VFeSb value

mat_props['Nv'] = [{'p' : 4, 'n' : 3}, {'p' : 4, 'n' : 3}, {'p' : 4, 'n' : 3}]

#Add the coefficients
mat_props = sig.fit_sigma_e0_coefficient([x_val, y_val, z_val], mat_props, 'p', T = 300)

'''
Print ternary
'''
u, cov = sig.fit_U(sigmae0_df, 'TaFeSb', 'VFeSb', mat_props, 'p', T = 300, p0 = 0)
sig_dict = sig.run_sigmae0_tern_data_dict(mat_props, 'p', T = 300, n = 100, U = u)
fig, ax = plt.subplots()
ax.axis("off")
figure, tax = ternary.figure(ax=ax, scale = 100)

tax.heatmap(sig_dict, style = 'h', cmap=plt.cm.get_cmap('Spectral_r', 20),\
     cbarlabel=r'$\sigma_{e0}$',\
     scientific = False)

Tt = [(sigmae0_df['TaFeSb'][i] * 100, sigmae0_df['VFeSb'][i] * 100, sigmae0_df['NbFeSb'][i] * 100) for i in  list(sigmae0_df.index)]
At = sigmae0_df['Sigma_E0']
tax.scatter(Tt, c = list(At), colormap=plt.cm.get_cmap('Spectral_r', 20),\
     cmap=plt.cm.get_cmap('Spectral_r', 20), vmin = 90000, vmax = 150000,\
     scientific = False, s = 30, edgecolors = 'k', zorder = 10, clip_on = False)

tax.boundary(linewidth=2.0)

tax.top_corner_label('VFeSb')
tax.left_corner_label('NbFeSb', position = (0,0.04, 0))
tax.right_corner_label('TaFeSb', position = (0.95,0.04, 0))

with open('XFeSb_sige0_data.csv', 'w') as csvfile:
    field_names = ['% (Ta)', '% (V)', '% (Nb)', 'Sigma_e0']
    writer = csv.DictWriter(csvfile, fieldnames  = field_names)
    writer.writeheader()
    for k,v in sig_dict.items():
        writer.writerow({'% (Ta)': k[0], '% (V)' : k[1], '% (Nb)' : 100 - (k[0] + k[1]), 'Sigma_e0' : v})
        
'''
Divide by Kappa_L and get quality factor
'''

fig, ax = plt.subplots()
ax.axis("off")
figure, tax = ternary.figure(ax=ax, scale = 100)

kL_df = pd.read_csv('FeXSb_data.csv')

sig_values = list(sig_dict.values())

T = 300

B = (kb / e)**2 * 300 * (np.array(sig_values) / np.array(kL_df['kappa_lattice']))


B_dict = sig_dict

j = 0
for k in B_dict.keys():
    B_dict[k] = B[j]
    j = j+1
    
tax.heatmap(B_dict, style = 'h', cmap=plt.cm.get_cmap('Spectral_r', 20),\
     cbarlabel=r'B',\
     scientific = False)

tax.boundary(linewidth=2.0)

tax.top_corner_label('VFeSb')
tax.left_corner_label('NbFeSb', position = (0,0.04, 0))
tax.right_corner_label('TaFeSb', position = (0.95,0.04, 0))

with open('XFeSb_Bfactor_data.csv', 'w') as csvfile:
    field_names = ['% (Ta)', '% (V)', '% (Nb)', 'BFactor']
    writer = csv.DictWriter(csvfile, fieldnames  = field_names)
    writer.writeheader()
    for k,v in B_dict.items():
        writer.writerow({'% (Ta)': k[0], '% (V)' : k[1], '% (Nb)' : 100 - (k[0] + k[1]), 'BFactor' : v})