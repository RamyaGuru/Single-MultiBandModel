#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 21 11:15:03 2021

@author: ramyagurunathan

Thermal Transport Model Quarternary

Need to try excluding the bad corner around SnSe?

And may need to fit an epsilon value. 

Also should fit epsilon to each binary and make sure they are similar? If not may need
to itnerpolate binary values in some way.


Epsilon value for PbTe-SnTe is much higher: decide how to interpolate epsilon? Is there a way to determine it experimentally for comparison?
"""
import sys
sys.path.append('../thermal_models/')
import ternary_tcond as tt

import numpy as np
import matplotlib.pyplot as plt

import math

'''
Constants
'''
kB = 1.38E-23
hbar = 1.054E-34
Na = 6.02E23


'''
endmember properties
'''

PbTe = {
      'vs' : 1850,
      'atmMass' :[207.2, 127.60],
      'atmRadius': [0.775, 0.97],
      'atmV': 35.379e-30,
      'natoms': 2,
      'stoich': [1,1],
      'k0' : 1.73, #W/ Ti doping
      }

PbSe = {
      'vs' : 1960,
      'atmMass' :[207.2, 78.971],
      'atmRadius': [0.775, 0.5],
      'atmV': 30.127e-30,
      'natoms': 2,
      'stoich': [1,1],
      'k0' : 1.43, #W/ Ti doping
      }

SnTe = {
      'vs' : 1800,
      'atmMass' :[118.71, 127.60], #changed from 118 to 80
      'atmRadius': [.69, 0.97],
      'atmV': 33.0335e-30,
      'natoms': 2,
      'stoich': [1,1],
#      'k0' : 12.9, #W/o the Titanium doping
      'k0' : 2.02, #W/ Ti doping
      }

#Just include this endmember even though it's a different  structure type?
SnSe = {
      'vs' : 1420,
      'atmMass' :[118.71, 78.971], #changed 78.971 60 127.60
      'atmRadius': [.69, 0.5],
      'atmV': 28.13e-30,
      'natoms': 2,
      'stoich': [1,1],
#      'k0' : 12.9, #W/o the Titanium doping
      'k0' : 0.99, #W/ Ti doping
      }

'''
Binary Data Lists
'''


#Could combine epsilons using the Muggianu model?
epsilon = 10


pred_kL, gamma, gammaM, gammaV = tt.run_kL_gamma_quat_data_dict(tt.gamma_tern, tt.gamma_tern, epsilon, 11, PbTe, [SnSe, PbSe, SnTe])

#Implement the quartenrary model for the htemral conductivity prediction

exp_kL = np.loadtxt('../datafiles/quat_tcond_chalc.csv',delimiter = ',')

print(exp_kL)

'''
Plot the Experimental Thermal Conductivity
'''
plt.matshow(exp_kL, cmap = plt.cm.get_cmap('rainbow'), vmin = 0.5, vmax = 2)
cbar = plt.colorbar()    
cbar.set_label(r'Lattice Thermal Conductivity (W/m/K)')

'''
Plot the Prediicted Themral Conductivity
'''
plt.matshow(pred_kL, cmap = plt.cm.get_cmap('rainbow'), vmin = 0.5, vmax = 2)
cbar = plt.colorbar()    
cbar.set_label(r'Lattice Thermal Conductivity (W/m/K)')


'''
Plot the Predicted Mass Variance
'''
print(gammaM)
plt.matshow(gammaM, cmap = plt.cm.get_cmap('rainbow'), vmin = 0, vmax = 0.08)
cbar = plt.colorbar()    
cbar.set_label(r'Mass Variance Term')


'''
doesn;t seem to be showing the endmember values for PbSe.. basically going to 0 for some reason?
'''
print(pred_kL)


'''

Get binaries and fit epislon to each one 
'''
PbTe_PbSe = np.array([np.linspace(1e-10,1 - 1e-10, 11), exp_kL[0,:]])

PbTe_SnTe = np.array([np.linspace(1e-10,1 - 1e-10, 11),exp_kL[:,0]])

SnTe_SnSe = np.array([np.linspace(1e-10,1 - 1e-10, 11), exp_kL[-1, :]])

PbSe_SnSe = np.array([np.linspace(1e-10,1 - 1e-10, 11), exp_kL[:,-1]])

PbTe_SnSe = np.array([np.linspace(1e-10,1 - 1e-10, 11), [exp_kL[j,j] for j in range(len(exp_kL))]])

SnTe_PbSe = np.array([np.linspace(1e-10,1 - 1e-10, 11), [exp_kL[len(exp_kL)-1-j,j] for j in range(len(exp_kL)-1, -1, -1)]])


eps_pbte_pbse, kL_pbte_pbse, kL_full_pbte_pbse = tt.fit_eps_kL(tt.gammaM_vac, tt.gammaV, PbTe_PbSe, PbTe, PbSe)

eps_pbte_snte, kL_pbte_snte, kL_full_pbte_snte = tt.fit_eps_kL(tt.gammaM_vac, tt.gammaV, PbTe_SnTe, PbTe, SnTe)

eps_snte_snse, kL_snte_snse, kL_full_snte_snse = tt.fit_eps_kL(tt.gammaM_vac, tt.gammaV, SnTe_SnSe, SnTe, SnSe)

eps_pbse_snse, kL_pbse_snse, kL_full_pbse_snse = tt.fit_eps_kL(tt.gammaM_vac, tt.gammaV, PbSe_SnSe, PbSe, SnSe)

eps_pbte_snse, kL_pbte_snse, kL_full_pbte_snse = tt.fit_eps_kL(tt.gammaM_vac, tt.gammaV, PbTe_SnSe, PbTe, SnSe)

eps_snte_pbse, kL_snte_pbse, kL_full_snte_pbse = tt.fit_eps_kL(tt.gammaM_vac, tt.gammaV, SnTe_PbSe, SnTe, PbSe)

'''
Print
'''

print(eps_pbte_pbse)
print(eps_pbte_snte) #This strain scattering term is much higher? combine using muggianu type model? weighted sum of epsilon values?
print(eps_snte_snse)
print(eps_pbse_snse)
print(eps_pbte_snse)
print(eps_snte_pbse)

'''
Make sure this is in the correct order: SnSe-PbSe, SnSe-SnTe, SnSe-PbTe, PbSe-SnTe, PbSe-PbTe, SnTe-PbTe
'''
eps_list = [eps_pbse_snse[0], eps_snte_snse[0], eps_pbte_snse[0], eps_snte_pbse[0], eps_pbte_pbse[0], eps_pbte_snte[0]]


'''
Updated quaternary plots: Muggianu model --> need to write a general muggianu model
'''

quat_kL, quat_gamma, quat_gammaM = tt.run_kL_gamma_quat_data_dict_mugg(tt.gammaM_vac, tt.gammaV, eps_list, 11, PbTe, [SnSe, PbSe, SnTe])

'''
Plot the Prediicted Themral Conductivity
'''
plt.matshow(quat_kL, cmap = plt.cm.get_cmap('rainbow'))
cbar = plt.colorbar()    
cbar.set_label(r'Lattice Thermal Conductivity (W/m/K)')


'''
Plot the Predicted Gamma
'''
plt.matshow(quat_gamma, cmap = plt.cm.get_cmap('rainbow'))
cbar = plt.colorbar()    
cbar.set_label(r'Total Gamma')

'''
Plot the Predicted Gamma
'''
plt.matshow(quat_gammaM, cmap = plt.cm.get_cmap('rainbow'), vmin = 0, vmax = 0.08)
cbar = plt.colorbar()    
cbar.set_label(r'Mass Gamma')
