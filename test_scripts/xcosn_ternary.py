#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 08:53:08 2020

@author: ramyagurunathan

XCoSn thermal conductivity prediction
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
      'atmV': ( 5.948E-10)**3 / 12,
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

'''
Binary PDScattering Model
'''

"""
Function: Mass difference equation from the current paper-- with vacancies
"""
def gammaM_vac(stoich, mass, subst, c):
    natoms = sum(stoich)
    delM2 = 0
    denom = 0
    for n in range(len(mass)):
        msite = subst[n]*c + mass[n]*(1-c)
        #delM2 = delM2 + stoich[n]*(c*(subst[n] - msite)**2 + (1-c)*(mass[n] - msite)**2)
        if subst[n] == 0 or mass[n] == 0:
            delM2 = delM2 + stoich[n]*c*(1-c)*(3*(subst[n] - mass[n]))**2
            denom = denom + stoich[n]*msite
        else:
            delM2 = delM2 + stoich[n]*c*(1-c)*(subst[n] - mass[n])**2
            denom = denom + stoich[n]*msite               
    gamma = (delM2/natoms)/((denom/natoms)**2)
    return gamma

'''
Function for radius difference scattering
'''
def gammaV(stoich, rad, subst, c):
    natoms = sum(stoich)
    delR2 = 0
    denom = 0
    for n in range(len(rad)):
        rsite = subst[n]*c + rad[n]*(1-c)        
        delR2 = delR2 + stoich[n]*c*(1-c)*(subst[n] - rad[n])**2
        denom = denom + stoich[n]*rsite               
    gamma = (delR2/natoms)/((denom/natoms)**2)
    return gamma    

'''
Fit Eps and Calculate Thermal Conductivity
'''
def kL_from_gamma(gamma, propA, propB, c):
    atmV = (1-c) * propA['atmV'] + c * propB['atmV']
    vs = (1-c) * propA['vs'] + c * propB['vs']
    k0 = (1-c) * propA['k0'] + c * propB['k0']
    prefix = (6**(1/3)/2)*(pi**(5/3)/kB)*(atmV**(2/3)/vs)
    u = (prefix * gamma * k0)**(1/2)
    kL = k0*np.arctan(u)/u
    return kL

def kL_tot(c, eps, Mfunc, Rfunc, propA, propB):
    gamma = Mfunc(propA['stoich'], propA['atmMass'], propB['atmMass'], c) +\
    eps * Rfunc(propA['stoich'], propA['atmRadius'], propB['atmRadius'], c)
    kL = kL_from_gamma(gamma, propA, propB, c)
    return kL

def fit_eps_kL(Mfunc, Rfunc, data, propA, propB):
    data = data[data[:,0].argsort()]
    eps, cov = curve_fit(lambda c, eps:\
                         kL_tot(c, eps, Mfunc, Rfunc, propA, propB), data[0,:],\
                                data[1,:], bounds = (0,np.inf))
    kL = np.zeros(100)
    i = 0
    for c in np.linspace(data[0,0], data[0, -1], 100):
        kL[i] = kL_tot(c, eps, Mfunc, Rfunc, propA, propB)
        i = i+1
    j = 0
    kL_full = np.zeros(100)
    for d in np.linspace(1e-10,9.9999999e-1,100):
        kL_full[j] = kL_tot(d, eps, Mfunc, Rfunc, propA, propB)
        j = j+1
    return eps, kL, kL_full

def run_kL(Mfunc, Rfunc, eps, propA, propB):
    kL_full = np.zeros(100)
    j = 0
    for d in np.linspace(1e-10,9.9999999e-1,100):
        kL_full[j] = kL_tot(d, eps, Mfunc, Rfunc, propA, propB)
        j = j+1
    return kL_full

#First ,fit epsilon to the Fe(V,Nb)Sb Data

eps, kL, kL_full = fit_eps_kL(gammaM_vac, gammaV, tcond_data, nb, ta)


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
kL_VTa = run_kL(gammaM_vac, gammaV,0, v, ta)
plt.plot(np.linspace(1e-10,9.9999999e-1,100), kL_VTa, color = 'xkcd:tree green')
plt.ylabel(r'$\kappa_L$ (W/m/K)')
plt.xlabel(r'V$_{1-x}$Ta$_x$CoSn')
plt.savefig('FeV_TaSb_kL.pdf', bbox_inches = 'tight')
'''
Plot Results (V,Nb)CoSn
'''
plt.figure()
kL_NbTa = run_kL(gammaM_vac, gammaV, 0, v, nb)
plt.plot(np.linspace(1e-10,9.9999999e-1,100), kL_NbTa, color = 'xkcd:blood red')
plt.ylabel(r'$\kappa_L$ (W/m/K)')
plt.xlabel(r'V$_{1-x}$Nb$_x$CoSn')
plt.savefig('FeNb_TaSb_kL.pdf', bbox_inches = 'tight')

'''
Ternary Methods: One method for gamma, just put radius in for mass when doing the volume term
'''  
#Need to Add Vacancy-handling back in 
def gamma_tern(stoich, mass, subst : list, c : list):
    n_sub = len(subst)
    defect_conc = sum(c)
    natoms = sum(stoich)
    delM2 = 0
    denom = 0
    for n in range(len(mass)):
        msite = mass[n] * (1-defect_conc)
        for s1 in range(n_sub):
            msite = msite  + c[s1] * subst[s1][n]
        delM2 = delM2 + stoich[n]* (1- defect_conc) * (mass[n] - msite)**2
        for s2 in range(n_sub):
            delM2 = delM2 + stoich[n]*c[s2]*(subst[s2][n] - msite)**2
        denom = denom + stoich[n]*msite  
#    print(delM2)             
    gamma = (delM2/natoms)/((denom/natoms)**2)
    return gamma  

def kL_from_gamma_tern(gamma, propA, propB : list, c : list):
    defect_conc = sum(c)
    atmV = (1 - defect_conc) * propA['atmV'] +  sum(c[i] * propB[i]['atmV'] for i in range(len(c)))
    vs = (1 - defect_conc) * propA['vs'] +  sum(c[i] * propB[i]['vs'] for i in range(len(c)))
    k0 = (1 - defect_conc) * propA['k0'] +  sum(c[i] * propB[i]['k0'] for i in range(len(c)))
    prefix = (6**(1/3)/2)*(pi**(5/3)/kB)*(atmV**(2/3)/vs)
    u = (prefix * gamma * k0)**(1/2)
    kL = k0*np.arctan(u)/u
    return kL   

def kL_tot_tern(c : list, eps, Mfunc, Rfunc, propA, propB : list): 
    gamma = Mfunc(propA['stoich'], propA['atmMass'], [p['atmMass'] for p in propB], c) +\
    eps * Rfunc(propA['stoich'], propA['atmRadius'], [p['atmRadius'] for p in propB], c)
    kL = kL_from_gamma_tern(gamma, propA, propB, c)
    return kL  

def fit_eps_kL_tern(Mfunc, Rfunc, data, propA, propB):
    data = data[data[:,0].argsort()]
    eps, cov = curve_fit(lambda c, eps:\
                         kL_tot_tern(c, eps, Mfunc, Rfunc, propA, propB), data[0,:],\
                                data[1,:], bounds = (0,np.inf))
    kL = np.zeros(100)
    i = 0
    for c in np.linspace(data[0,0], data[0, -1], 100):
        kL[i] = kL_tot_tern(c, eps, Mfunc, Rfunc, propA, propB)
        i = i+1
    j = 0
    kL_full = np.zeros(100)
    for d in np.linspace(1e-10,9.9999999e-1,100):
        kL_full[j] = kL_tot_tern(d, eps, Mfunc, Rfunc, propA, propB)
        j = j+1
    return eps, kL, kL_full  


#Currently only works for tenrary
def run_kL_tern(Mfunc, Rfunc, eps, propA, propB : list, n):
    kL_full = np.zeros(n**2)
#    n_sub = len(propB)
    k = 0
    t = []
    u = []
    v = []
    for c in np.linspace(1e-10,9.9999999e-1,n):
        for d in np.linspace(1e-10, 9.9999999e-1 - c, n):
            t.append(c)
            u.append(d)
            v.append(1-c-d)
            kL_full[k] = kL_tot_tern([c,d], eps, Mfunc, Rfunc, propA, propB)
            k = k+1
    return kL_full, [t,u,v]

first = 1e-10
last = 9.99999999e-1

def run_kL_tern_data_dict(Mfunc, Rfunc, eps, propA, propB : list, n = 10):
    kL_full = dict()
#    n_sub = len(propB)
    j = 0
    for c in np.arange(first,1, (last - first) / n):
        k = 0
        for d in np.arange(first, 1 - c, (last - first) / n):
            kL_full[(c*100,d*100)] = kL_tot_tern([c,d], eps, Mfunc, Rfunc, propA, propB)
            k = k+1
        j = j+1
    return kL_full

'''
Plot Results (V,Nb,Ta)CoSn
'''
#plt.figure()
kL_tern = run_kL_tern_data_dict(gamma_tern, gamma_tern, 0, v, [ta, nb], 200)
#plt.plot(np.linspace(1e-10,9.9999999e-1,100), kL_VTa_tern, color = 'xkcd:tree green')
#plt.ylabel(r'$\kappa_L$ (W/m/K)')
#plt.xlabel(r'FeV$_{1-x}$Ta$_x$Sb')
#plt.savefig('FeV_TaSb_kL.pdf', bbox_inches = 'tight')

'''
Generate data for a ternary diagram
'''
data_scatter = [(0, 0, 100), (0, 100, 0), (100, 0, 0), (60, 40, 0)]
kL_scatter = [14.45, 7.72, 5.94, 4.8]

fig, ax = plt.subplots()
ax.axis("off")
figure, tax = ternary.figure(ax=ax, scale = 100)

tax.heatmap(kL_tern, style = 'h', cmap=plt.cm.get_cmap('Spectral_r', 20),
             cbarlabel=r'$\kappa$ (W/m/K)',
             vmax=10.0, vmin=3.5, scientific = False)

tax.scatter(data_scatter, c = kL_scatter, colormap=plt.cm.get_cmap('Spectral_r', 20),\
     cmap=plt.cm.get_cmap('Spectral_r', 20), vmin = 3.5, vmax = 10,\
     scientific = False, s = 30, edgecolors = 'k', zorder = 10, clip_on = False)   
    
    
tax.boundary(linewidth=2.0)

tax.top_corner_label('NbCoSn')
tax.left_corner_label('VCoSn', position = (0,0.04, 0))
tax.right_corner_label('TaCoSn', position = (0.95,0.04, 0))

tax.savefig('klemens_model_xcosn.pdf', bbox_inches = 'tight')



with open('XCoSn_data.csv', 'w') as csvfile:
    field_names = ['% (Ta)', '% (Nb)', '% (V)', 'kappa_lattice']
    writer = csv.DictWriter(csvfile, fieldnames  = field_names)
    writer.writeheader()
    for k,v in kL_tern.items():
        writer.writerow({'% (Ta)': k[0], '% (Nb)' : k[1], '% (V)' : 100 - (k[0] + k[1]), 'kappa_lattice' : v})
    