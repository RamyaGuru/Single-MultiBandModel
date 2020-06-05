#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 10:25:01 2020

@author: ramyagurunathan

Fe2VAl

Thermal Model for Fe2VAl
"""

import numpy as np
from math import pi
from scipy.optimize import curve_fit


'''
Constants
'''
kB = 1.38e-23 # in V / K
hbar = 1.054e-34 # in J * s
Na = 6.02e23

'''
Material Properties
'''
atmMass = 47.403 #Average atomic mass in atomic mass units
atmV = 11.62e-30 #Avg. atomic volume in cubic meters
natoms = 4

'''
Sound velocities
'''
v_l = 7750
v_t = 4530

def v_sound(v_l, v_t):
    v_s = ((1 / (3 * v_l**3)) + (2 / (3 * v_t**3)))**(-1/3)
    return v_s

def debyeT(atmV, vs):
    return (hbar /kB) * (6*pi / atmV)**(1/3) * vs

vs = v_sound(v_l, v_t)

dT = debyeT(atmV, vs)

debyef = (dT * kB)/hbar

'''
Original Masses
'''
stoich = [2, 1, 1]
og_mass = [55.845, 50.942, 26.982]
og_rad = [.75, .86, .675]
'''
New Masses
'''
new_mass = {
        'Germanium': [55.845, 50.942, 72.630],
        'Silicon': [55.845, 50.942, 28.085],
        'Cobalt':[58.933, 50.942, 26.982],
        'Rhenium' : [186.21, 50.942, 26.982],
        'Titanium' : [55.845, 47.867, 26.982],
        'Tantalum_Al': [55.845, 180.95, 26.982],
                }
new_radius = {
        'Germanium' : [.75, .86, .87],
        'Silicon' : [.75, .86, .54],
        'Cobalt' : [.79, .86, .675],
        'Rhenium': [.77, .86, .675],
        'Titanium': [.75, .88, .675],
        'Tantalum_Al': [.75, .86, .88],
        }

'''
Elastic constants
'''
nu = .24
G = 4



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


def kL_from_gamma(gamma, kap_pure):
    prefix = (6**(1/3)/2)*(pi**(5/3)/kB)*(atmV**(2/3)/vs)
    u = (prefix*gamma*kap_pure)**(1/2)
    kL = kap_pure*np.arctan(u)/u
    return kL
    
     
def kLcalc(data, stoich, mass, subst):
    data = data[data[:,0].argsort()]
    kL = np.zeros(len(data))
    i = 0
    for c, kappa in zip(data[:,0], data[:,1]):
        kap_pure = data[0,1]
#        kap_pure = (1-c)*data[0,1] + c*14
        gamma = gammaM_vac(stoich, mass, subst, c)
        print(gamma)
        prefix = (6**(1/3)/2)*(pi**(5/3)/kB)*(atmV**(2/3)/vs)
        u = (prefix*gamma*kap_pure)**(1/2)
        kL[i] = kap_pure*np.arctan(u)/u
        i = i+1
    return kL

def kL_tot(c, eps, Mfunc, Rfunc, stoich, mass, msubst, rad, rsubst, data):
    gamma = Mfunc(stoich, mass, msubst, c) + eps * Rfunc(stoich, rad, rsubst, c)
    kap_pure = data[0,1]        
#        kap_pure = (1-c)*data[0,1] + c*14
    kL = kL_from_gamma(gamma, kap_pure)
    return kL

def fit_eps_kL(Mfunc, Rfunc, stoich, mass, msubst, rad, rsubst, data):
    data = data[data[:,0].argsort()]
    eps, cov = curve_fit(lambda c, eps:\
                         kL_tot(c, eps, Mfunc, Rfunc, stoich, mass,\
                                msubst, rad, rsubst, data), data[:,0],\
                                data[:,1], bounds = (0,np.inf))
    kL = np.zeros(100)
    i = 0
    for c in np.linspace(data[0,0], data[-1, 0], 100):
        kL[i] = kL_tot(c, eps, Mfunc, Rfunc, stoich, mass, msubst, rad, rsubst, data)
        i = i+1
    j = 0
    kL_full = np.zeros(100)
    for d in np.linspace(1e-10,9.9999999e-1,100):
        kL_full[j] = kL_tot(d, eps, Mfunc, Rfunc, stoich, mass, msubst, rad, rsubst, data)
        j = j+1
    return eps, kL, kL_full



      
'''
kappa vs. temperature: Cobalt system?

Use two-band model for the Lorentz number when calculating the electronic portion
Conductivity has come directly from experimental data
'''  
#keep Gruneisen as a fitting parameter
def umklapp_tau(T, grun, freq):
    return (6 * pi**2)**(1/3)/2 * ((atmMass / Na) * 1e-3 * vs**3)\
 / (kB * atmV**(1/3) * float(grun)**2 * freq**2 * T)

#PD scattering term
def pd_tau(gamma, freq):
    return 4 * pi * vs**3/ (atmV * freq**4 * gamma)

def spectral_C(freq, T):
    x = hbar * freq / (kB * T)
    C = (3 / (2 * pi**2)) * kB * (freq**2/ vs**3) * (x**2 * np.exp(x))/ (np.exp(x) - 1)**2
    return C

def kL_spectral(T, grun, gamma, freq, i):
    tau = 1/(umklapp_tau(T, grun, freq)**(-1) + pd_tau(gamma, freq)**(-1))
    return spectral_C(freq, T) * vs**2 * tau


n_f = 1000
dfreq = debyef/n_f   


def kL_T_integral(T, grun, gamma):
    kL_int = 0
    i = 0
    for freq in np.arange(dfreq, debyef, dfreq):
        kL_int = kL_int + kL_spectral(T, grun, gamma, freq, i) * dfreq
        i = i+1
    return (1/3) * kL_int



def fit_grun(gamma, data_df, c):
    grun = curve_fit(lambda t, grun: kL_T_integral(t, grun, gamma), data_df['T'],\
                     data_df[str(c)], bounds = (0, np.inf))
    return grun
    

        
    