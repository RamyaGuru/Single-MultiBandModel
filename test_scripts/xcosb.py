#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 12:13:55 2020

@author: ramyagurunathan

Brown Compound Family:(Ti,Hf,Zr)CoSb
"""

import sys
sys.path.append('/Users/ramyagurunathan/Documents/PhDProjects/Single_Multiband_Models')
from math import pi
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib as mpl
import ternary
import plotly.figure_factory as ff
from plotly.offline import plot
import xlwt
import csv
from pymatgen.ext.matproj import MPRester
import ternary_tcond as tt

mpr = MPRester('pfSJBa1OwitR5uNL')

mpl.rcdefaults()
mpl.rcParams['font.sans-serif'] = 'Apple Symbols'
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['lines.markersize'] = 4
mpl.rcParams['axes.xmargin'] = 0.1
mpl.rcParams['axes.formatter.useoffset'] = False
mpl.rcParams['font.weight'] = 'bold'
mpl.rcParams['font.size'] = '18'

mpl.rcParams['mathtext.fontset'] = 'custom'

mpl.rcParams['mathtext.bf'] = 'Apple Symbols'

'''
Constants
'''
kB = 1.38E-23
hbar = 1.054E-34
Na = 6.02E23


ti = {
      'atmV' : (5.885E-10)**3 / 12,
      'natoms' : 3,
      'stoich': [1,1,1],
      'atmMass' : [47.87, 58.9, 121.76],
      'atmRadius' : [1, 0.54, 0.9],
      'debyeT' : 416,
      'k0': 14.691,
      }

ti['vs'] = tt.vs_from_debyeT(ti['atmV'], ti['debyeT'])

#Source: Sekimoto IEEE Xplore
hf = {
      'atmV' : (6.0390E-10)**3 / 12,
      'natoms' : 3,
      'stoich': [1,1,1],
      'atmMass' : [178.49, 58.9, 121.76],
      'atmRadius' : [0.85, 0.54, 0.9],
      'debyeT' : 340,
      'k0': 11.88613116,
      }

hf['vs'] = tt.vs_from_debyeT(hf['atmV'], hf['debyeT'])

#ZrCoSn, Source: Silpawilawan JMC C
zr = {
      'atmV' : (6.0697E-10)**3 / 12,
      'natoms' : 3,
      'stoich': [1,1,1],
      'atmMass' : [91.224, 58.9, 121.76],
      'atmRadius' : [0.86, 0.54, 0.9],
      'debyeT' : 392,
      'k0': 15.58614565,
      'gruneisen' : 1.44,
      }
zr['vs'] = tt.vs_from_debyeT(zr['atmV'], zr['debyeT'])


'''
TCond Data: Use the study that is likely mulitple phases..?
'''

zr_ti = [[],[]]
zr_ti[0] = [1e-10, 0.5, 0.6, 0.8 ,0.9, 9.999999e-1]
zr_ti[1] = [15.586, 5.96, 5.65, 5.98, 7.79 ,14.691]
zr_ti = np.array(zr_ti)

eps_zr_ti, kL_zr_ti, kL_full_zr_ti = tt.fit_eps_kL(tt.gammaM_vac, tt.gammaV, zr_ti, zr, ti)

'''
Calculate the full ternary
'''
kL_tern = tt.run_kL_tern_data_dict(tt.gamma_tern, tt.gamma_tern, eps_zr_ti[0], zr, [hf, ti])

fig, ax = plt.subplots()
ax.axis("off")
figure, tax = ternary.figure(ax=ax, scale = 100)

tax.heatmap(kL_tern, style = 'h', cmap=plt.cm.get_cmap('viridis_r', 15),
             cbarlabel=r'$\kappa$ (W/m/K)',
             vmax=10, vmin=3.4, scientific = False)
    
    
tax.boundary(linewidth=2.0)

tax.top_corner_label('TiCoSb')
tax.left_corner_label('ZrCoSb', position = (0,0.04, 0))
tax.right_corner_label('HfCoSb', position = (0.95,0.04, 0))

fig.savefig('klemens_model_xcosb.pdf', bbox_inches = 'tight')

with open('XCoSb_data.csv', 'w') as csvfile:
    field_names = ['% (Hf)', '% (Ti)', '% (Zr)', 'kappa_lattice']
    writer = csv.DictWriter(csvfile, fieldnames  = field_names)
    writer.writeheader()
    for k,v in kL_tern.items():
        writer.writerow({'% (Hf)': k[0], '% (Ti)' : k[1], '% (Zr)' : 100 - (k[0] + k[1]), 'kappa_lattice' : v})