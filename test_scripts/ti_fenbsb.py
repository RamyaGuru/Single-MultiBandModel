#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 11 23:32:35 2020

@author: ramyagurunathan

Thermal Model for the Ti-doped FeNbSb
"""

import thermal_model as tm
import sys
sys.path.append('/Users/ramyagurunathan/Documents/PhDProjects/UQ_SPB')
import klat_temperature as kt

'''
Constants
'''
kB = 1.38e-23 # in V / K
hbar = 1.054e-34 # in J * s
Na = 6.02e23

D = {}
D['og'] = {'mass': [55.845, 92.906, 121.76],'rad': [.75, .86, 0.9]}
D['subst'] = {'mass': [55.845, 92.906, 121.76] ,'rad': [.75, .67, 0.9]}
D['vs'] = 3052
D['avgM'] = sum([55.845, 92.906, 121.76]) / 3
D['avgV'] =(53.167E-30) / 3 #in cubic meters
D['N'] = 3
D['stoich'] = [1,1,1]
D['nfreq'] = 100
D['c'] = 0.2
D['grun'] = 1.5

D['debyef'] = kt.debyeT(D['vs'], D['avgV']) * (hbar / kB)
'''
Data
'''
D['Tt'], D['At'], D['Et'], D['It'] = \
kt.read_data_pd('/Users/ramyagurunathan/Documents/PhDProjects/Argonne_TECCA/UQData/FeNbSb/Titanium','20') 


'''
Fit an epsilon value
'''

n_f = 1000
dfreq = D['debyef']/n_f
eps = 100

gammaM = kt.gamma(D['stoich'], D['og']['mass'], D['subst']['mass'], D['c'])

gammaV = kt.gamma(D['stoich'], D['og']['rad'], D['subst']['rad'], D['c'])

gamma = gammaM + eps * gammaV

kappa = tm.kL_T_integral(300, D['grun'], gamma)