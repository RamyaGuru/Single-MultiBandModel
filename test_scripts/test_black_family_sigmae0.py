#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 15:28:06 2021

@author: ramyagurunathan

TaCoSn Black Family
"""

data_path = '/Users/ramyagurunathan/Documents/PhDProjects/Argonne_TECCA/SigmaE0_values/'
datafile = data_path + 'XFeSb_gauss_data_red_2nd_Gn.csv'
sigmae0_df = sig.fetch_sigmae0_dataframe(datafile)
x_val, y_val, z_val = sig.fetch_endmember_values(sigmae0_df)
'''
Initalize Endmember Properties
Order of compounds: TaFeSb, VFeSb, NbFeSb

Average Effective Mass values: Obtained from the MPContribs database
Dividing the carr. conc., eff. mass, 
'''
#Bulk Modulus sources Ta: Materials Project
mat_props = {'formula': ['TaCoSn', 'NbCoSn', 'VCoSn'], 'BulkMod': [163E9, 102E9, 201E9],\
             'AtmV' : [( 5.948E-10)**3 / 12, (5.946E-10)**3 / 12, 1.745E-29]}



#Add effective masses. Sources: Nb, V: A. Page J. of MAteriomics 2016, Ta: 
mat_props['eff_mass_p'] = [1.65, 5.1, 1.6] #Check the VFeSb value

mat_props['Nv'] = [{'p' : 10, 'n' : 3}, {'p' : 10, 'n' : 3}, {'p' : 10, 'n' : 3}] #P may be different

#Add the coefficients
mat_props = sig.fit_sigma_e0_coefficient([x_val, y_val, z_val], mat_props, 'p', T = 300)

'''
Print ternary
'''
u, cov = fit_U(sigmae0_df, 'TaFeSb', 'VFeSb', mat_props, 'p', T = 300, p0 = 0)
sig_dict = run_sigmae0_tern_data_dict(mat_props, 'p', T = 300, n = 100, U = u)
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

