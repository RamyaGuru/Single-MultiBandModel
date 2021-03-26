#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 13:15:25 2021

@author: ramyagurunathan


Toy Plot of the Mixing Models
"""

import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib as mpl

mpl.rcParams['figure.figsize'] = [5,5]
mpl.rcParams['font.size'] = 14

x = np.linspace(1e-5, 9.9999E-1, 40)

L = 20

R = 8.3145

T = 3

#def V(x):
#    return

def S(x):
    return -R * (x * math.log(x) + (1-x) * math.log(1-x))

def H(x):
    return x * (1-x) * L

def G(x, G1, G2):
    g_veg = G1 * (1-x) + G2 * x
    h = H(x)
    s = S(x)
    return g_veg + h - T * s


plt.figure()

plt.plot(x, [S(i) for i in x])
plt.plot(x, [H(i) for i in x])

plt.plot(x, [G(i, 0, 3) for i in x])

plt.plot(x, [3 * i for i in x], ':', color = 'xkcd:black')

plt.xlim([0, 1])
plt.yticks([])

plt.savefig('toy_model_mixing.pdf', bbox_inches = 'tight')
