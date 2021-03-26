#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 11:39:13 2021

@author: ramyagurunathan

Kazu Mg3X2 ternary
"""
import sys
sys.path.append('../thermal_models/')
import ternary
import ternary_tcond as tt
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from math import pi as pi


'''
Constants
'''
kB = 1.38E-23
hbar = 1.054E-34
Na = 6.02E23

#Source: 6.95

'''
Endmemeber Data: Temperature of 400K

Gruneisen: http://www.rsc.org/suppdata/c8/ta/c8ta07285j/c8ta07285j1.pdf
https://spj.sciencemag.org/journals/research/2020/1934848/

DebyeT: https://www.sciencedirect.com/science/article/abs/pii/S2542529318301147

'''

sb_ = {
      'atmV' : 26.1E-30,
      'natoms' : 5,
      'stoich': [3,2],
      'atmMass' : [24.31, 121.76],
      'atmRadius' : [.57, 0.76],
      'debyeT' : 230,
      'k0': 1.181,
      'gruneisen' : 1.83,
      }
sb_['vs'] = tt.vs_from_debyeT(sb_['atmV'], sb_['debyeT'])

bi_ = {
      'atmV' : 27.7E-30,
      'natoms' : 5,
      'stoich': [3,2],
      'atmMass' : [24.31, 208.98],
      'atmRadius' : [.57, 1.03],
      'debyeT' : 177,
      'k0': 1.08, #1.549,
      'gruneisen' : 1.94,
      }
bi_['vs'] = tt.vs_from_debyeT(bi_['atmV'], bi_['debyeT'])


#Velocity for Mg3Bi2
vl = 3220
vt = 1690

avg_vs_Bi = ((1/3) * (2 / (vt**3)) + (1 / (vl**3)))**(-1/3)

#Velocity for Mg3Sb2
vl = 4130
vt = 2190

avg_vs_Sb = ((1/3) * (2 / (vt**3)) + (1 / (vl**3)))**(-1/3)

#Velocity for Mg3As2
vl = 5120
vt = 2800

avg_vs_As = ((1/3) * (2 / (vt**3)) + (1 / (vl**3)))**(-1/3)


as_ =  {
      'atmV' : 21.69E-30,
      'natoms' : 5,
      'stoich': [3,2],
      'atmMass' : [24.31, 74.92],
      'atmRadius' : [.57, .58],
      'vs' : avg_vs_As,
      'k0': 4.76, #fit this value? scaled from Mg3Bii2
      'gruneisen' : 1.3,
      }

def guess_kappa_endpt(known, guess):
    avgM1 = sum(known['atmMass']) / len(known['atmMass'])
    coeff1 = avgM1 * known['vs']**3 / (known['atmV']**(2/3) * known['gruneisen']**2)
    avgM2 = sum(guess['atmMass']) / len(guess['atmMass'])
    coeff2 = avgM2 * guess['vs']**3 / (guess['atmV']**(2/3) * guess['gruneisen']**2)
    ratio = coeff2 / coeff1
    return ratio * known['k0']

def toberer_kappa_U_model(props, T):
    coeff = (6 * pi**2)**(2/3) / ((4 * pi**2) * T)
    avgM = sum(props['atmMass']) / len(props['atmMass'])
    return coeff * (props['vs']**3 * (avgM / Na) * 1e-3) / (props['atmV']**(2/3) * props['gruneisen']**2) *\
    (1 / (props['natoms']**(1/3)))
    
T = 400
    
k_tob = toberer_kappa_U_model(as_, T)    

ko1 = guess_kappa_endpt(sb_,as_)
ko2 = guess_kappa_endpt(bi_,as_)

datafile = '/Users/ramyagurunathan/Documents/PhDProjects/Argonne_TECCA/Mg3X2_data/400K_ternary'

data_df = pd.read_csv(datafile, sep = '\t')

sb_bi = np.array([[data_df['Mg2Bi3'][i] for i in list(data_df.index)[:-3]], [data_df['kL (W/m/K)'][i] for i in list(data_df.index)[:-3]]])
eps_bin,kL,kL_full = tt.fit_eps_kL(tt.gammaM_vac, tt.gammaV, sb_bi, sb_, bi_)

xyz_scatter = [(data_df['Mg2Sb3'][i] * 100, data_df['Mg2Bi3'][i] * 100, data_df['Mg2As3'][i] * 100)\
               for i in list(data_df.index)]
kL_scatter = data_df['kL (W/m/K)']

sb_conc = [data_df['Mg2Sb3'][i] for i in list(data_df.index)]
bi_conc = [data_df['Mg2Bi3'][i] for i in list(data_df.index)]

sb_bi_array = np.array([sb_conc, bi_conc])

'''
Plot binary
'''

plt.figure()

plt.scatter(sb_bi[0] * 100, sb_bi[1])
plt.plot(kL_full)
plt.xlabel(r'Mg$_3$(Sb$_{1-x)}$Bi$_x$)$_2$')
plt.ylabel(r'$\kappa_L$ (W/m/K)')
plt.savefig('Mg3Sb_Bi2_binary_kL.pdf', bbox_inches = 'tight')

plt.figure()
conc = np.linspace(1e-10,9.9999999e-1,100)
as_bi_binary = tt.run_kL(tt.gammaM_vac, tt.gammaV, eps_bin, bi_, as_)
plt.plot(conc, as_bi_binary)
plt.xlabel(r'Mg$_3$(Bi$_{1-x}$As$_x$)$_2$')
plt.ylabel(r'$\kappa_L$ (W/m/K)')

#%%

'''
Calculate the kL
'''
eps_tern = tt.fit_eps_kL_tern(tt.gamma_tern, tt.gamma_tern, sb_bi_array, kL_scatter,\
                              as_, [bi_, sb_]) 

eps_zero = 0
kL_tern = tt.run_kL_tern_data_dict(tt.gamma_tern, tt.gamma_tern, eps_bin, 200, as_, [bi_, sb_])
'''
Plot the ternary with data
'''
mpl.rcParams['figure.figsize'] = [10,6] 
mpl.rcParams['font.size'] = 18

fig, ax = plt.subplots()
ax.axis("off")
figure, tax = ternary.figure(ax=ax, scale = 100)

tax.boundary(linewidth=2.0)

tax.heatmap(kL_tern, style = 'h', cmap=plt.cm.get_cmap('Spectral_r', 20),\
     cbarlabel=r'$\kappa_L$ (W/m/K)',vmin = 0.5, vmax = 1.5,\
     scientific = False)

tax.scatter(xyz_scatter, c = kL_scatter, colormap=plt.cm.get_cmap('Spectral_r', 20),\
     cmap=plt.cm.get_cmap('Spectral_r', 20), cbarlabel = r'$\kappa_L$ (W/m/K)',\
     scientific = False, s = 50, edgecolors = 'k', zorder = 10,\
     vmin = 0.5, vmax = 1.5, clip_on = False)  

tax.boundary(linewidth=2.0)

tax.top_corner_label(r'Mg$_3$Bi$_2$', position = (-0.08, 1.15, 0))
tax.left_corner_label(r'Mg$_3$As$_2$', position = (0,0.05, 0))
tax.right_corner_label(r'Mg$_3$Sb$_2$', position = (0.95,0.05, 0))

tax.savefig('Mg3X2_ternary_4_kL.pdf', bbox_inches = 'tight')



