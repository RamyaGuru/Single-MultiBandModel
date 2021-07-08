#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 11 10:52:52 2021

@author: ramyagurunathan

Ternayr Klemens Toy Model
"""

import ternary_tcond as tt
import numpy as np
import sys
sys.path.append('../thermal_models/')
from math import pi as pi
import matplotlib.pyplot as plt
import ternary

'''
Constants
'''
kB = 1.38E-23
hbar = 1.054E-34
Na = 6.02E23

def debyeT(atmV, vs):
    return (hbar /kB) * (6*pi**2 / atmV)**(1/3) * vs

def vs_from_debyeT(atmV, debyeT):
    return (kB / hbar) * (6*pi**2 / atmV)**(-1/3) * debyeT

'''
FeNbSb Material Properties : Silpawilawan JMC C
'''
B = {
    'vs' : 3052,
    'atmMass' : [55.845, 75, 121.76],
    'atmRadius': [.75, .5, 0.9],
    'atmV' : 1.8e-29,
    'natoms' : 3,
    'stoich': [1,1,1],
#    'k0' : 17.9,
    'k0' : 10,
    'gruneisen' : 1.49
    }

#nb['dT'] = debyeT(nb['atmV'], nb['vs'])

'''
TaFeSb Material Properties: Grytsiv Intermetallics 111
'''

C = {
      'vs' : 3052,
      'atmMass' :[55.845, 100, 121.76],
      'atmRadius': [.75, 1.0, 0.9],
      'atmV': 1.8e-29,
      'natoms': 3,
      'stoich': [1,1,1],
#      'k0' : 12.9, #W/o the Titanium doping
      'k0' : 10, #W/ Ti doping
      'gruneisen': 1.49
      }
'''
VFeSb Material Properties : C. Fu, et al. JAP 112, 124915 (2012)
'''

A = {
     'vs' : 3052,
     'atmMass' : [55.845, 50 ,121.76],
     'atmRadius': [.75, 1.5, 0.9],
     'atmV' : 1.8e-29,
     'natoms' : 3,
     'stoich': [1,1,1],
     'k0': 10, #W/o the Titanium doping
     'gruneisen' : 1.49
     }

D = {
     'vs' : 3052,
     'atmMass' : [55.845, 100 ,121.76],
     'atmRadius': [.75, 0.5, 0.9],
     'atmV' : 1.8e-29,
     'natoms' : 3,
     'stoich': [1,1,1],
     'k0': 10, #W/o the Titanium doping
     'gruneisen' : 1.49
     }

epsilon = .25

kL_tern, gamma_tern, gammaM_tern, gammaV_tern = tt.run_kL_gamma_tern_data_dict(tt.gamma_tern, tt.gamma_tern, epsilon, 100, B, [C,A,D])

    
fig, ax = plt.subplots()
ax.axis("off")
figure, tax = ternary.figure(ax=ax, scale = 100)

tax.heatmap(kL_tern, style = 'h', cmap=plt.cm.get_cmap('Spectral_r', 20),\
     cbarlabel=r'$\kappa$ (W/m/K)',vmin = 0,\
     scientific = False)


tax.boundary(linewidth=2.0)

tax.top_corner_label(r'A')
tax.left_corner_label(r'B', position = (0,0.04, 0))
tax.right_corner_label(r'C', position = (0.90,0.04, 0))
    
fig, ax = plt.subplots()
ax.axis("off")
figure, tax = ternary.figure(ax=ax, scale = 100)

tax.heatmap(gamma_tern, style = 'h', cmap=plt.cm.get_cmap('Spectral_r', 20),\
     cbarlabel=r'$\Gamma$' ,vmin = 0,\
     scientific = False)


tax.boundary(linewidth=2.0)

tax.top_corner_label(r'A')
tax.left_corner_label(r'B', position = (0,0.04, 0))
tax.right_corner_label(r'C', position = (0.90,0.04, 0))
    

fig, ax = plt.subplots()
ax.axis("off")
figure, tax = ternary.figure(ax=ax, scale = 100)

tax.heatmap(gammaM_tern, style = 'h', cmap=plt.cm.get_cmap('Spectral_r', 20),\
     cbarlabel=r'$\Gamma$' ,vmin = 0,\
     scientific = False)


tax.boundary(linewidth=2.0)

tax.top_corner_label(r'A')
tax.left_corner_label(r'B', position = (0,0.04, 0))
tax.right_corner_label(r'C', position = (0.90,0.04, 0))



fig, ax = plt.subplots()
ax.axis("off")
figure, tax = ternary.figure(ax=ax, scale = 100)

tax.heatmap(gammaV_tern, style = 'h', cmap=plt.cm.get_cmap('Spectral_r', 20),\
     cbarlabel=r'$\Gamma$' ,vmin = 0,\
     scientific = False)


tax.boundary(linewidth=2.0)

tax.top_corner_label(r'A')
tax.left_corner_label(r'B', position = (0,0.04, 0))
tax.right_corner_label(r'C', position = (0.90,0.04, 0))


'''
A-D binary
'''
plt.figure()
kL, gamma = tt.run_kL_gamma(tt.gammaM_vac, tt.gammaV, epsilon, A, D)

plt.plot(np.linspace(1e-10,9.9999999e-1,100), gamma)


