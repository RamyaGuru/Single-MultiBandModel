#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 09:14:04 2021

@author: ramyagurunathan

Silicon spectral thermal conductivity
"""

'''
constants
'''

from math import pi as pi

kB = 1.38E-23
hbar = 1.054E-34
Na = 6.02E23

'''
Silicon values
'''
vs = 6084
T = 300
a = 2.7e-10
grun = 0.56
atmM = (28.05 / Na) / 1e3

k_full = (6 * pi**2)**(2 / 3) * (1 / (4 * pi**2)) * (vs**3 / (T * a**2 * grun**2)) * atmM


k_spec = k_full / (vs * (2 * pi) / a)