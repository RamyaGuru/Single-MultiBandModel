#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 10:09:17 2020

@author: ramyagurunathan

XNiSn Family

Need to be aware of the huge miscibility gap
"""
import ternary_tcond as tt
import numpy as np
import matplotlib.pyplot as plt
import ternary
import csv

'''
Constants
'''
kB = 1.38E-23
hbar = 1.054E-34
Na = 6.02E23

#Source: 6.95

ti = {
      'atmV' : (5.95E-10)**3 / 12,
      'natoms' : 3,
      'stoich': [1,1,1],
      'atmMass' : [47.87, 58.69, 118.71],
      'atmRadius' : [1, 0.69, 0.69],
      'debyeT' : 407,
      'k0': 6.95,
      'gruneisen' : 1.6,
      }

ti['vs'] = tt.vs_from_debyeT(ti['atmV'], ti['debyeT'])

#Source: Sekimoto IEEE Xplore
hf = {
      'atmV' : (6.12E-10)**3 / 12, 
      'natoms' : 3,
      'stoich': [1,1,1],
      'atmMass' : [178.49, 58.69, 118.71],
      'atmRadius' : [0.85, 0.69, 0.69],
      'debyeT' : 307,
      'k0': 8.25,
      }

hf['vs'] = tt.vs_from_debyeT(hf['atmV'], hf['debyeT'])


zr = {
      'atmV' : (6.15E-10)**3 / 12,
      'natoms' : 3,
      'stoich': [1,1,1],
      'atmMass' : [91.22, 58.69, 118.71],
      'atmRadius' : [0.86, 0.69, 0.69],
      'debyeT' : 323,
      'k0': 8.75,
      }
zr['vs'] = tt.vs_from_debyeT(zr['atmV'], zr['debyeT'])

#Binary Thermal Conductivity Data
#Source:
ti_hf = [[],[]]
ti_hf[0] = [1e-10, 0.05, 9.999999e-1] #Hf site fraction
ti_hf[1] = [6.95, 5.21, 8.25]
ti_hf = np.array(ti_hf)

eps_ti_hf, kL_ti_hf, kL_full_ti_hf = tt.fit_eps_kL(tt.gammaM_vac, tt.gammaV, ti_hf, ti, hf)
#Source: Zou J.Appl.Phys. 2013
zr_hf = [[],[]]
zr_hf[0] = [1e-10, 0.5, 9.999999e-1] #Hf site fraction
zr_hf[1] = [8.75, 5, 8.25] #Hf site fraction
zr_hf = np.array(zr_hf)

eps_zr_hf, kL_zr_hf, kL_full_zr_hf = tt.fit_eps_kL(tt.gammaM_vac, tt.gammaV, zr_hf, zr, hf)
#Source: Rogl, Acta Materialia 113
zr_ti = [[],[]]
zr_ti[0] = [1e-10, 0.5, 9.999999e-1]
zr_ti[1] = [8.72, 5.25, 6.95]
zr_ti = np.array(zr_ti)

eps_zr_ti, kL_zr_ti, kL_full_zr_ti = tt.fit_eps_kL(tt.gammaM_vac, tt.gammaV, zr_ti, zr, ti)

'''
Calculate the full ternary
'''
kL_tern = tt.run_kL_tern_data_dict(tt.gamma_tern, tt.gamma_tern, eps_ti_hf[0], 200, ti, [hf, zr])


'''
Generate data for a ternary diagram
'''
data_scatter = [(0, 0, 100), (5, 0, 99.5),\
                (100, 0, 0), (0, 100, 0),\
                (50, 50, 0), (0, 50, 50)]
kL_scatter = [6.95, 5.21, 8.25, 8.75, 5, 5.25]

fig, ax = plt.subplots()
ax.axis("off")
figure, tax = ternary.figure(ax=ax, scale = 100)

tax.heatmap(kL_tern, style = 'h', cmap=plt.cm.get_cmap('Spectral_r', 20),
            vmin = 2, vmax = 9,
             cbarlabel=r'$\kappa$ (W/m/K)',
             scientific = False)

tax.scatter(data_scatter, c = kL_scatter, colormap=plt.cm.get_cmap('Spectral_r', 20),\
     cmap=plt.cm.get_cmap('Spectral_r', 20), vmin = 2, vmax = 9,\
     scientific = False, s = 30, edgecolors = 'k', zorder = 10, clip_on = False)    
    
tax.boundary(linewidth=2.0)

tax.top_corner_label('ZrNiSn')
tax.left_corner_label('TiNiSn', position = (0,0.04, 0))
tax.right_corner_label('HfNiSn', position = (0.95,0.04, 0))

tax.savefig('klemens_model_xnisn_old.pdf', bbox_inches = 'tight')

with open('XNiSn_data.csv', 'w') as csvfile:
    field_names = ['% (Hf)', '% (Ti)', '% (Zr)', 'kappa_lattice']
    writer = csv.DictWriter(csvfile, fieldnames  = field_names)
    writer.writeheader()
    for k,v in kL_tern.items():
        writer.writerow({'% (Hf)': k[0], '% (Ti)' : k[1], '% (Zr)' : 100 - (k[0] + k[1]), 'kappa_lattice' : v})
    