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

mpr = MPRester('pfSJBa1OwitR5uNL')

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

zr_ti = [[],[]]
zr_ti[0] = [1e-10, 0.5, 0.6, 0.8 ,0.9, 9.999999e-1]
zr_ti[1] = [15.586, 5.96, 5.65, 5.98, 7.79 ,14.691]
zr_ti = np.array(zr_ti)

eps_zr_ti, kL_zr_ti, kL_full_zr_ti = tt.fit_eps_kL(tt.gammaM_vac, tt.gammaV, zr_ti, zr, ti)

plt.plot(np.linspace(0, 1, 100), kL_full_zr_ti)

