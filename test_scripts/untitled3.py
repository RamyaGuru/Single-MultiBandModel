#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 11:39:13 2021

@author: ramyagurunathan

Kazu Mg3X2 ternary
"""
import ternary
import ternary_tcond as tt
import pandas as pd


'''
Constants
'''
kB = 1.38E-23
hbar = 1.054E-34
Na = 6.02E23

#Source: 6.95

'''
Endmemeber Data: Temperature of 400K

Gruneisen: http://www.rsc.org/suppdata/c8/ta/c8ta07285j/c8ta07285j1.pdf
https://spj.sciencemag.org/journals/research/2020/1934848/

DebyeT: https://www.sciencedirect.com/science/article/abs/pii/S2542529318301147

'''

sb_ = {
      'atmV' : 26.1E-30,
      'natoms' : 5,
      'stoich': [3,2],
      'atmMass' : [24.31, 121.76],
      'atmRadius' : [.57, 0.76],
      'debyeT' : 230,
      'k0': 1.181,
      'gruneisen' : 1.83,
      }

bi_ = {
      'atmV' : 27.7E-30,
      'natoms' : 5,
      'stoich': [3,2],
      'atmMass' : [24.31, 208.98],
      'atmRadius' : [.57, 1.03],
      'debyeT' : 177,
      'k0': 1.549,
      'gruneisen' : 1.94,
      }

as_ =  {
      'atmV' : 21.69E-30,
      'natoms' : 5,
      'stoich': [3,2],
      'atmMass' : [24.31, 74.92],
      'atmRadius' : [.57, .58],
      'debyeT' : 250,
      'k0': 1.7, #fit this value?
      'gruneisen' : 1.94,
      }

datafile = ''

data_df = pd.read_csv()
