#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 12:55:34 2021

@author: ramyagurunathan
XCoSn Using the ternary_tcond class

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
import csv
from pymatgen.ext.matproj import MPRester
import ternary_tcond as tt

mpr = MPRester('pfSJBa1OwitR5uNL')

mpl.rcParams['figure.figsize'] = [8,6]
mpl.rcParams['font.size'] = 18

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

def average_phase_velocity(vp):
    avg_vp = (1/3) * ( vp[0]**(-3) + vp[1]**(-3) + vp[2]**(-3))**(-1/3)
    return avg_vp

'''
NbCoSn Material Properties : Li ACS Appl. Mater. Interfaces
'''
nb = {
    'vs' : 3368,
    'atmMass' : [92.906, 58.933, 118.71],
    'atmRadius': [.86, .54, 0.69],
    'atmV' : (5.946E-10)**3 / 12,
    'natoms' : 3,
    'stoich': [1,1,1],
    'k0' : 7.72,
    'gruneisen' : (1.97 + 2.16 + 1.95) / 3
    }

#nb['dT'] = debyeT(nb['atmV'], nb['vs'])

'''
TaCoSn Material Properties: Li ACS Appl. Mater. Interfaces
'''

ta = {
      'atmMass' :[180.95, 58.933 ,118.71],
      'atmRadius': [.86, .54, 0.69],
      'atmV': (5.948E-10)**3 / 12,
      'natoms': 3,
      'stoich': [1,1,1],
      'k0' : 5.94, #W/o the Titanium doping
      'gruneisen': (2.15 + 2.35 + 2.23) / 3
      }
ta['vs'] = 3077

'''
VCoSn Material Properties : Zaferani (definitely mulitple phases in the sample)
'''
struct = mpr.get_structure_by_material_id('mp-1018119')
vol = struct.volume
natoms = struct.composition.num_atoms
nsites = struct.num_sites
elem = struct.species
speciesMass = struct.composition.weight
atmMass = speciesMass/natoms
        
mass_density = 1.6605E3 * nsites * speciesMass /\
(natoms * vol)

elast = mpr.get_data('mp-1018119') 
bulkMod = elast[0]['elasticity']['K_VRH'] #Bulk modulus
shearMod = elast[0]['elasticity']['G_VRH'] #shear modulus
trans_v = (1E9 * abs(shearMod) / mass_density)**0.5
long_v = (1E9 * (abs(bulkMod) + 4./3. * abs(shearMod)) /\
                mass_density)**0.5
avg_v = 3**(1/3)*(1/long_v**3 + 2/trans_v**3)**(-1/3)

v = {
     'vs' : avg_v,
     'atmMass' : [50.942, 58.933 ,118.71],
     'atmRadius': [.78, .54, 0.69],
     'atmV' : vol * 1E-30 / nsites,
     'natoms' : 3,
     'stoich': [1,1,1],
     'k0': 14.45, #W/o the Titanium doping
     }

'''
Data for the (Nb,Ta)CoSn 
'''
x_Ta = [1e-10, 0.6, 9.9999999e-01]
kappa = [7.72, 4.8, 5.94]
tcond_data=[[],[]]
tcond_data[0] = x_Ta
tcond_data[1] = kappa
tcond_data = np.array(tcond_data)



eps, kL, kL_full = tt.fit_eps_kL(tt.gammaM_vac, tt.gammaV, tcond_data, nb, ta)
eps_cov = 0
'''
Plot Results ((Nb, Ta)CoSn)
'''
plt.scatter(tcond_data[0,:], tcond_data[1,:])
plt.plot(np.linspace(1e-10,9.9999999e-1,100), kL_full, color = 'xkcd:light indigo')
plt.ylabel(r'$\kappa_L$ (W/m/K)')
plt.xlabel(r'Nb$_{1-x}$Ta$_x$CoSn')
plt.savefig('FeV_NbSb_kL.pdf', bbox_inches = 'tight')
'''
Plot Results (V,Ta)CoSn
'''
plt.figure()
kL_VTa = tt.run_kL(tt.gammaM_vac, tt.gammaV,eps, v, ta)
plt.plot(np.linspace(1e-10,9.9999999e-1,100), kL_VTa, color = 'xkcd:tree green')
plt.ylabel(r'$\kappa_L$ (W/m/K)')
plt.xlabel(r'V$_{1-x}$Ta$_x$CoSn')
plt.savefig('FeV_TaSb_kL.pdf', bbox_inches = 'tight')
'''
Plot Results (V,Nb)CoSn
'''
plt.figure()
kL_NbTa = tt.run_kL(tt.gammaM_vac, tt.gammaV, eps, v, nb)
plt.plot(np.linspace(1e-10,9.9999999e-1,100), kL_NbTa, color = 'xkcd:blood red')
plt.ylabel(r'$\kappa_L$ (W/m/K)')
plt.xlabel(r'V$_{1-x}$Nb$_x$CoSn')
plt.savefig('FeNb_TaSb_kL.pdf', bbox_inches = 'tight')


kL_tern, gamma_full = tt.run_kL_cov_tern_data_dict(tt.gamma_tern, tt.gamma_tern, eps[0], eps_cov, 200, nb, [ta, v])

eps_list = [eps[0], eps[0], eps[0]]
kL_mugg, gamma_mugg = tt.run_kL_tern_data_dict_muggianu(tt.gammaM_vac, tt.gammaV, eps_list, 200, nb, [ta, v])

data_scatter = [(0, 100, 0), (0, 0, 100), (100, 0, 0), (60, 0, 40)]
kL_scatter = [14.45, 7.72, 5.94, 4.8]

'''
OG Plots
'''
fig, ax = plt.subplots()
ax.axis("off")
figure, tax = ternary.figure(ax=ax, scale = 100)

tax.heatmap(kL_tern, style = 'h', cmap=plt.cm.get_cmap('Spectral_r', 15),
             cbarlabel=r'$\kappa$ (W/m/K)',
             vmax=10, vmin=3.4, scientific = False)

tax.scatter(data_scatter, c = kL_scatter, colormap=plt.cm.get_cmap('Spectral_r', 20),\
     cmap=plt.cm.get_cmap('Spectral_r', 20), vmin = 3.5, vmax = 10,\
     scientific = False, s = 30, edgecolors = 'k', zorder = 10, clip_on = False)   

tax.boundary(linewidth=2.0)

tax.top_corner_label('VCoSn')
tax.left_corner_label('NbCoSn', position = (0,0.04, 0))
tax.right_corner_label('TaCoSn', position = (0.95,0.04, 0))


tax.savefig('xcosn_model_experiment.pdf', bbox_inches = 'tight')

fig, ax = plt.subplots()
ax.axis("off")
figure, tax = ternary.figure(ax=ax, scale = 100)

tax.heatmap(gamma_full, style = 'h', cmap=plt.cm.get_cmap('Spectral_r', 15),
             scientific = False)


'''
Muggianu Plots
'''

fig, ax = plt.subplots()
ax.axis("off")
figure, tax = ternary.figure(ax=ax, scale = 100)

tax.heatmap(kL_mugg, style = 'h', cmap=plt.cm.get_cmap('Spectral_r', 15),
             cbarlabel=r'$\kappa$ (W/m/K)',
             vmax=10, vmin=3.4, scientific = False)

tax.scatter(data_scatter, c = kL_scatter, colormap=plt.cm.get_cmap('Spectral_r', 20),\
     cmap=plt.cm.get_cmap('Spectral_r', 20), vmin = 3.5, vmax = 10,\
     scientific = False, s = 30, edgecolors = 'k', zorder = 10, clip_on = False)   

tax.boundary(linewidth=2.0)

tax.top_corner_label('NbCoSn')
tax.left_corner_label('VCoSn', position = (0,0.04, 0))
tax.right_corner_label('TaCoSn', position = (0.95,0.04, 0))

fig, ax = plt.subplots()
ax.axis("off")
figure, tax = ternary.figure(ax=ax, scale = 100)

tax.heatmap(gamma_mugg, style = 'h', cmap=plt.cm.get_cmap('Spectral_r', 15),
             scientific = False)
