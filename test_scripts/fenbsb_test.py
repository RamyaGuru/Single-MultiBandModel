#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 23:21:48 2020

@author: ramyagurunathan

(Nb, V, Ta)FeSb
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

'''
FeNbSb Material Properties : Silpawilawan JMC C
'''
nb = {
    'vs' : 3052,
    'atmMass' : [55.845, 92.906, 121.76],
    'atmRadius': [.75, .86, 0.9],
    'atmV' : (53.167E-30) / 3,
    'natoms' : 3,
    'stoich': [1,1,1],
    'k0' : 17.9,
    'gruneisen' : 1.49
    }

#nb['dT'] = debyeT(nb['atmV'], nb['vs'])

'''
TaFeSb Material Properties: Grytsiv Intermetallics 111
'''

ta = {
      'atmMass' :[55.845, 180.95, 121.76],
      'atmRadius': [.75, .86, 0.9],
      'atmV': (5.9355E-10)**3 / 12,
      'natoms': 3,
      'stoich': [1,1,1],
      'k0' : 12.9, #W/o the Titanium doping
      'gruneisen': 1.49
      }
ta['vs'] = vs_from_debyeT(ta['atmV'], 351)
'''
VFeSb Material Properties : C. Fu, et al. JAP 112, 124915 (2012)
'''

v = {
     'vs' : 2374,
     'atmMass' : [55.845, 50.942 ,121.76],
     'atmRadius': [.75, .78, 0.9],
     'atmV' : (5.7976E-10)**3 / 12,
     'natoms' : 3,
     'stoich': [1,1,1],
     'k0': 12.5, #W/o the Titanium doping
     'gruneisen' : 1.49
     }

'''
(Nb, V)FeSb Binary Data
'''
x_Nb = [1e-10, 0.1, 0.2 ,0.3 ,0.4, 9.999999E-1]
kappa = [12.5, 7.46, 6.39, 6.07, 5.65, 17.9]
tcond_data=[[],[]]
tcond_data[0] = x_Nb
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
    print(vs)
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

def fit_eps_kL(Mfunc, Rfunc, data, propA, propB, n):
    data = data[data[:,0].argsort()]
    eps, cov = curve_fit(lambda c, eps:\
                         kL_tot(c, eps, Mfunc, Rfunc, propA, propB), data[0,:],\
                                data[1,:], bounds = (0,np.inf))
    kL = np.zeros(n)
    i = 0
    for c in np.linspace(data[0,0], data[0, -1], n):
        kL[i] = kL_tot(c, eps, Mfunc, Rfunc, propA, propB)
        i = i+1
    j = 0
    kL_full = np.zeros(n)
    for d in np.linspace(1e-10,9.9999999e-1,n):
        kL_full[j] = kL_tot(d, eps, Mfunc, Rfunc, propA, propB)
        j = j+1
    return eps, kL, kL_full

def run_kL(Mfunc, Rfunc, eps, propA, propB, n):
    kL_full = np.zeros(n)
    j = 0
    for d in np.linspace(1e-10,9.9999999e-1,n):
        kL_full[j] = kL_tot(d, eps, Mfunc, Rfunc, propA, propB)
        j = j+1
    return kL_full

#First ,fit epsilon to the Fe(V,Nb)Sb Data

eps, kL, kL_full = fit_eps_kL(gammaM_vac, gammaV, tcond_data, v, nb, 4)


'''
Plot Results (Fe(V,Nb)Sb)
'''
plt.scatter(tcond_data[0,:], tcond_data[1,:])
plt.plot(np.linspace(1e-10,9.9999999e-1,4), kL_full, color = 'xkcd:light indigo')
plt.ylabel(r'$\kappa_L$ (W/m/K)')
plt.xlabel(r'FeV$_{1-x}$Nb$_x$Sb')
plt.savefig('FeV_NbSb_kL.pdf', bbox_inches = 'tight')
'''
Plot Results Fe(V,Ta)Sb
'''
plt.figure()
kL_VTa = run_kL(gammaM_vac, gammaV, eps, v, ta, 4)
plt.plot(np.linspace(1e-10,9.9999999e-1,4), kL_VTa, color = 'xkcd:tree green')
plt.ylabel(r'$\kappa_L$ (W/m/K)')
plt.xlabel(r'FeV$_{1-x}$Ta$_x$Sb')
plt.savefig('FeV_TaSb_kL.pdf', bbox_inches = 'tight')
'''
Plot Results Fe(Nb,Ta)Sb
'''
plt.figure()
kL_NbTa = run_kL(gammaM_vac, gammaV, eps, nb, ta, 4)
plt.plot(np.linspace(1e-10,9.9999999e-1,4), kL_NbTa, color = 'xkcd:blood red')
plt.ylabel(r'$\kappa_L$ (W/m/K)')
plt.xlabel(r'FeNb$_{1-x}$Ta$_x$Sb')
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
#    trum(delM2)             
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
def run_kL_tern(Mfunc, Rfunc, eps, propA, propB : list):
    kL_full = np.zeros(10000)
#    n_sub = len(propB)
    k = 0
    t = []
    u = []
    v = []
    for c in np.linspace(1e-10,9.9999999e-1,100):
        for d in np.linspace(1e-10, 9.9999999e-1 - c, 100):
            t.append(c)
            u.append(d)
            v.append(1-c-d)
            kL_full[k] = kL_tot_tern([c,d], eps, Mfunc, Rfunc, propA, propB)
            k = k+1
    return kL_full, [t,u,v]

first = 1e-10
last = 9.99999999e-1
def run_kL_tern_data_dict(Mfunc, Rfunc, eps, propA, propB : list):
    kL_full = dict()
#    n_sub = len(propB)
    j = 0
    for c in np.arange(first,1,(last - first) / 10):
        k = 0
        for d in np.arange(first, 1 - c, (last - first) / 10):
            kL_full[(c*100,d*100)] = kL_tot_tern([c,d], eps, Mfunc, Rfunc, propA, propB)
            k = k+1
        j = j+1
    return kL_full


'''
Plot Results Fe(V,Ta)Sb
'''
#plt.figure()
kL_tern = run_kL_tern_data_dict(gamma_tern, gamma_tern, eps[0], v, [ta, nb])
#plt.plot(np.linspace(1e-10,9.9999999e-1,100), kL_VTa_tern, color = 'xkcd:tree green')
#plt.ylabel(r'$\kappa_L$ (W/m/K)')
#plt.xlabel(r'FeV$_{1-x}$Ta$_x$Sb')
#plt.savefig('FeV_TaSb_kL.pdf', bbox_inches = 'tight')

'''
Generate data for a ternary diagram
'''
fig, ax = plt.subplots()
ax.axis("off")
figure, tax = ternary.figure(ax=ax, scale = 100)

tax.heatmap(kL_tern, style = 'h', cmap=plt.cm.get_cmap('viridis_r', 12),
             cbarlabel=r'$\kappa$ (W/m/K)',
             vmax=12.0, vmin=3.0, scientific = False)
    
    
tax.boundary(linewidth=2.0)

tax.top_corner_label('FeNbSb')
tax.left_corner_label('FeVSb', position = (0,0.04, 0))
tax.right_corner_label('FeTaSb', position = (0.95,0.04, 0))


fig.savefig('klemens_model.pdf', bbox_inches = 'tight')

'''
Trying plotly instead to get the contour lines to show up
'''

with open('FeXSb_data.csv', 'w') as csvfile:
    field_names = ['% (Ta)', '% (Nb)', '% (V)', 'kappa_lattice']
    writer = csv.DictWriter(csvfile, fieldnames  = field_names)
    writer.writeheader()
    for k,v in kL_tern.items():
        writer.writerow({'% (Ta)': k[0], '% (Nb)' : k[1], '% (V)' : 100 - (k[0] + k[1]), 'kappa_lattice' : v})
    
#kL, coords = run_kL_tern(gamma_tern, gamma_tern, eps[0], v, [ta, nb])
#fig = ff.create_ternary_contour(coords, kL, ncontours = 20, colorscale='Viridis', showscale=True)
#plot(fig, auto_open = True)