#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 09:56:23 2021

@author: ramyagurunathan

Agne diffson model
"""

#import eigenvector_overlap as eo
from phonopy import Phonopy
from phonopy.interface.vasp import read_vasp
from phonopy.file_IO import parse_FORCE_SETS, parse_BORN
import numpy as np
import matplotlib.pyplot as plt
from math import pi as pi
import ternary_tcond as tt

'''
Si : Material properties
'''
V = (2.70E-10)**3
avgM = 28
debyeT = 498
vs = tt.vs_from_debyeT(V, debyeT)
grun = 1.0

'''
Constants
'''
kB = 1.38e-23
hbar = 1.054e-34
Na = 6.022e23

T = 300

debyef = (kB / hbar) * debyeT

#def get_phonon_dos():
#    '''
#    Phonon density-of-states
#    '''


def analytic_tau_U(omega, T):
    '''
    omega : frequency
    '''
    return ((6 * pi**2)**(1 / 3) / 2) * ((avgM / Na) * vs**3) / (kB * V ** (1/3) *\
             grun ** 2 * omega ** 2 * T)

def diffusion_coefficient_agne(omega):
    return step_size**2 * (2 * omega / (2 * pi))

def diffusion_coefficient_tau(tau, step_size):
    '''
    Diffusion coeffcient calculation from the DFT phonon properties
    step_size : center of mass of vibration?
    '''
    return step_size**2 * (1/tau)

def heat_capacity(omega):
    return 3 * kB * omega**2 / (2 * pi**2 * vs**3)
    
def integrate_kappaL(Dcoeff_array, Cp_array, step_size, T, n_freq):
    d_omega = debyef / n_freq
    kL_array = Cp_array * Dcoeff_array * d_omega
    return sum(kL_array)

if __name__ == "__main__":
    step_size = V**(1/3)
    n_freq = 100
    omega = np.linspace(1e-10, debyef, n_freq)
    tau_array = analytic_tau_U(omega, T)
    Dcoeff_tau_array = diffusion_coefficient_tau(tau_array, step_size)
    Dcoeff_agne_array = diffusion_coefficient_agne(omega)
    Cp_array = heat_capacity(omega)
    kL_agne = integrate_kappaL(Cp_array, Dcoeff_agne_array, step_size, T, n_freq)
    kL_tau = integrate_kappaL(Cp_array, Dcoeff_tau_array, step_size, T, n_freq)
    '''
    kL_tau is significantly lower than kL_agne because of the lower energy transfer rate in this model
    '''
    
    
