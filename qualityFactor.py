#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 21:57:15 2020

@author: ramyagurunathan

Quality Factor
"""

'''
Constants
'''

from math import pi as pi
from fdint import fdk
import numpy as np


'''
Constants
'''

kB = 1.38e-23 # V / K
h = 6.02e-34 # J * s
hbar = 1.054e-34

e = 1.602e-19 # C
me = 9.11e-31 # kg
Na = 6.02e23 # /mol

'''
Material Properties/ inputs
'''

mol_mass = 2 * 55.845 + 50.942 + 26.98
density = 6770 # kg/m**3

T = 300 #K
kL = 1 #W/m/K


'''
n-type and p-type specific properties
'''

dpg_lvl = np.linspace(-8,4,1000)

eff_mass = {'n': 12.75, 'p': 4.7}

weighted_mobility = {'n':700e-4, 'p': 300e-4} # m**2/V/s

weighted_mobility_p = {'n': 700e-4, 'p': 300e-4} #m**2/V/s

band_gap = np.array([0.04, -0.1])

red_bgap = (band_gap * e) / (kB * T)

'''
Single Parabolic Band Model:
    Relationships between transport coefficients and doping level (Fermi level)
'''

def carr_conc_from_eta(eta, mstar, T=300):
    return 4 * pi * (2 * mstar * me * kB * T / h**2)**(3/2) * fdk(1/2, eta)


def seebeck(eta):
    return (kB / e) * (2 * fdk(1, eta) / fdk(0, eta) - eta)

def lorentz(eta):
    return (kB**2 / e**2) * (3 * fdk(0, eta) * fdk(2, eta) - 4 * fdk(1, eta)**2)\
/(fdk(0,eta)**2)

def conductivity(eta, sigmae0, s=1):
    if s == 0:  # s=0 requires analytic simplification
        return sigmae0 / (1. + np.exp(-eta))
    else:
        return sigmae0 * s * fdk(s - 1, eta)
     

def sigmae0_from_muW(muW, T):
    return (8 * pi * e * (2 * me * kB * T)**(3/2) / (3 * h**3)) * muW

def zT(eta, B):
    S = seebeck(eta)
    L = lorentz(eta)
    return S**2 / (((kB / e)**2 / (B * np.log(1 + np.exp(eta)))) + L)

'''
Two Band Model (EXTEND TO ARBITRARY NUMBER OF BANDS)
'''


def twoband_carrconc(eta_n, eta_p, mstar1, mstar2):
    '''
    Carrier concentration in m^{-3}
    '''
    return carr_conc_from_eta(eta_n, mstar1) + -1 * carr_conc_from_eta(eta_p, mstar2)
    
    
def twoband_conductivity(eta_n, eta_p, sigmae01, sigmae02, s=1):
    '''
    Note: sigmae01 should be the n-type value
    Conductivity in S/m
    '''
    return conductivity(eta_n, sigmae01, s=1) + conductivity(eta_p, sigmae02, s=1)

def cond_n_and_p(eta_n, eta_p, sigmae01, sigmae02, s=1):
    return conductivity(eta_n, sigmae01, s=1), conductivity(eta_p, sigmae02, s=1)

def twoband_seebeck(eta_n, eta_p, sigmae01, sigmae02):
    '''
    Note: sigmae01 should be the n-type value
    Seebeck in V/K
    '''
    return (-1 * seebeck(eta_n) * conductivity(eta_n, sigmae01, s=1) + seebeck(eta_p) *\
 conductivity(eta_p, sigmae02, s=1)) / (twoband_conductivity(eta_n, eta_p, sigmae01, sigmae02, s=1))

def twoband_kappa_e(eta_n, eta_p, sigmae01, sigmae02, T):
    '''
    may need to double check term 1
    '''
    sigma1 = conductivity(eta_n, sigmae01, 1)
    sigma2 = conductivity(eta_p, sigmae02, 1)
    term1 = T * (lorentz(eta_n) * sigma1 + lorentz(eta_p) * sigma2)
    term2 = T * ((seebeck(eta_n)**2 * sigma1 + seebeck(eta_p)**2 * sigma2) -\
                 (-1 * seebeck(eta_n) * sigma1 + seebeck(eta_p) * sigma2)**2 / (sigma1 + sigma2))
    return term1 +term2

def twoband_kappa(eta_n, eta_p, sigmae01, sigmae02, kL, T):
    '''
    Thermal conductivity in W/m/K
    '''
    kappa_e = twoband_kappa_e(eta_n, eta_p, sigmae01, sigmae02, T)
    return kappa_e + kL

def twobandzT(eta_n, eta_p, sigmae01, sigmae02, kL, T):  
    cond = twoband_conductivity(eta_n, eta_p, sigmae01, sigmae02, s=1)
    seebeck = twoband_seebeck(eta_n, eta_p, sigmae01, sigmae02)
    kappa = twoband_kappa(eta_n, eta_p, sigmae01, sigmae02, kL, T)
    return seebeck**2 * cond * T / kappa

'''
Convert between carrier concentration and dopant concentration assuming 100%
doping efficiency
'''
    

def carrconc_from_dopantconc(x):
    '''
    x is a dopant; will assume it contributes 1 carrier and scatters
    more or less like Ge
    '''
    return x * Na * (1/mol_mass) * (density * 1e3)
    
def dopant_conc_from_carrconc(n):
    return n * (1 / Na) * (mol_mass) * (1/(density * 1e3))


carrconc = np.zeros(len(dpg_lvl))
seeb = np.zeros(len(dpg_lvl))
kappa = np.zeros(len(dpg_lvl))
cond = np.zeros(len(dpg_lvl))
if __name__ == "__main__":
    for rbg in red_bgap:
        i = 0
        for eta_n in dpg_lvl:
            eta_p = -(rbg + eta_n)
            sigmae0_n = sigmae0_from_muW(weighted_mobility['n'], T)
            sigmae0_p = sigmae0_from_muW(weighted_mobility['p'], T)
            seeb[i] = twoband_seebeck(eta_n, eta_p, sigmae0_n, sigmae0_p)
            cond[i] = twoband_conductivity(eta_n, eta_p, sigmae0_n, sigmae0_p)
            kappa[i] = twoband_kappa(eta_n, eta_p, sigmae0_n, sigmae0_p, kL, T)
            carrconc[i] = twoband_carrconc(eta_n, eta_p, eff_mass['n'], eff_mass['p'])
            i = i+1


        
    