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

x = np.linspace(1e-5, 9.9999E-1, 40)

L = 20

R = 8.3145

T = 3

def S(x):
    return -R * (x * math.log(x) + (1-x) * math.log(1-x))

def H(x):
    return x * (1-x) * L

def G(x):
    h = H(x)
    s = S(x)
    return h - T * s


plt.figure()

plt.plot(x, [S(i) for i in x])
plt.plot(x, [H(i) for i in x])

plt.plot(x, [G(i) for i in x])

plt.plot(x, [0 for i in x], ':', color = 'xkcd:black')

plt.xlim([0, 1])
plt.yticks([])


