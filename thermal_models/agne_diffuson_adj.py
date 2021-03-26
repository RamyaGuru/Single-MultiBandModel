#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 13:04:04 2021

@author: ramyagurunathan


Implementation of the random walk model

Using DFT inputs

Avoid usage of the 

jump distance ** 2 * Gamma (scattering rate)
"""

import numpy as np
import ternary_tcond as tt
from math import pi as pi

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

T = 300

def tau_U(omega, T):
    '''
    omega : frequency
    '''
    return ((6 * pi**2)**(1 / 3) / 2) * (avgM * vs**3) / (kB * V ** (1/3) *\
             grun ** 2 * omega ** 2 * T)
    

#def gamma_PD(mat_props):
#    '''
#    mat_props : dict for the 
#    '''
#    return
#    
#def tau_PD(omega):
#    '''
#    omega : frequency
#    '''
#    return 

def diffusivity():
    '''
    
    '''