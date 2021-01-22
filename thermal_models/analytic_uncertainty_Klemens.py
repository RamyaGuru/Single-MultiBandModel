#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  9 20:12:54 2021

@author: ramyagurunathan

Propagate uncertainity through Klemens model

If only variable is the epsilon parameter, then can just solve for Gamma and do
a linear fit

If fitting kappa0 and epsilon.. then should use curvefit 

--> curve_fit gives variance of the parameters (epsilon)
--> Gamma_tot variance should be that times the Gamma_strain parameter
--> kappa variance should be more complicated... propagate uncertainty (std. dev.) through the trig function?
take the derivative and then multiply by the uncertainty

--> plot kappa versus gamma to see validity of Newton's law.. 
"""

import sys
sys.path.append('/Users/ramyagurunathan/Documents/PhDProjects/UQ_SPB')

import ternary_tcond as tt
import thermal_model as tm
import ternary_UQ_XNiSn as tuq
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from math import pi

mpl.rcParams['font.size'] = 14

'''
Constants
'''
kB = 1.38E-23
hbar = 1.054E-34
Na = 6.02E23

'''
Endmember Properties
'''
ti = {
      'atmV' : (5.95E-10)**3 / 12,
      'natoms' : 3,
      'stoich': [1,1,1],
      'atmMass' : [47.87, 58.69, 118.71],
      'atmRadius' : [1, 0.69, 0.69],
      'debyeT' : 407,
      'k0': 3.176,#6.95,
      'gruneisen' : 1.6,
      }

ti['vs'] = tt.vs_from_debyeT(ti['atmV'], ti['debyeT'])

#Source: Sekimoto IEEE Xplore
hf = {
      'atmV' : (6.12E-10)**3 / 12, #are these correct??
      'natoms' : 3,
      'stoich': [1,1,1],
      'atmMass' : [178.49, 58.69, 118.71],
      'atmRadius' : [0.85, 0.69, 0.69],
      'debyeT' : 307,
      'k0': 3.4626, #8.25,
      }

hf['vs'] = tt.vs_from_debyeT(hf['atmV'], hf['debyeT'])

#ZrCoSn, Source: Silpawilawan JMC C
zr = {
      'atmV' : (6.15E-10)**3 / 12,
      'natoms' : 3,
      'stoich': [1,1,1],
      'atmMass' : [91.22, 58.69, 118.71],
      'atmRadius' : [0.86, 0.69, 0.69],
      'debyeT' : 323,
      'k0': 3.486, #8.75,
      }
zr['vs'] = tt.vs_from_debyeT(zr['atmV'], zr['debyeT'])


data_dir = '/Users/ramyagurunathan/Documents/PhDProjects/Argonne_TECCA/UQData/TiNiSn/PurpleFamily/'
allk_name_list = ['Akram2016', 'Appel2014' ,'Appel2015', 'Bhardwaj2012', 'Fan2014', 'Gurth2016',\
             'Kim2006','Kurosaki2004','Muta2005', 'Populoh2013', 'Sakurada2005', 'Schwall2013',\
              'Shutoh2004', 'Tang2009', 'Yu2012']
name_list = ['Gurth2016']
data = tuq.read_data_pd(data_dir, name_list)
Tt, At,  Et, It = tuq.get_data(data)
eps,cov = tt.fit_eps_cov_kL_tern(tt.gamma_tern, tt.gamma_tern, list(Tt), At, ti, [hf, zr], p0 = 0)

'''
Plot of Kappa versus Gamma
'''
n = 10

first = 1e-10
last = 9.999999e-1

i = 1

kL_list = []
gamma_list = []

c = 0.1

for d in np.arange(first, 1 - c, (last - first) / n):
    kL_tot, gamma = tt.kL_tot_tern_gamma((c, d), eps, tt.gamma_tern, tt.gamma_tern, ti, [hf, zr])
    kL_list.append(kL_tot)
    gamma_list.append(gamma[0])


plt.plot(gamma_list, kL_list)
plt.ylabel(r'$\kappa_L$')

plt.xlabel(r'$\Gamma$')

'''
Trivially propagate through Gamma
'''

def Gamma_cov(C, propA, propB, Rfunc, eps_cov):
    Gamma_R = Rfunc(propA['stoich'], propA['atmRadius'], [p['atmRadius'] for p in propB], list(C))
    return eps_cov * Gamma_R

'''
Test
'''
Gam_cov = []
for d in np.arange(first, 1 - c, (last - first) / n):
    Gam_cov.append(Gamma_cov((c,d), ti, [hf, zr], tt.gamma_tern, cov))

'''
Use Newton's laws to propagate to kappa
'''

def d_kL_d_Gamma(gamma, propA, propB : list, c : list):
    defect_conc = sum(c)
    host_conc = (1 + 1e-10) - defect_conc
    atmV = (host_conc) * propA['atmV'] +  sum(c[i] * propB[i]['atmV'] for i in range(len(c)))
    vs = (host_conc) * propA['vs'] +  sum(c[i] * propB[i]['vs'] for i in range(len(c)))
    k0 = (host_conc) * propA['k0'] +  sum(c[i] * propB[i]['k0'] for i in range(len(c)))
    prefix = (6**(1/3)/2)*(pi**(5/3)/kB)*(atmV**(2/3)/vs)
    prefix_u = ((6**(1/3)/2)*(pi**(5/3)/kB)*(atmV**(2/3)/vs) * k0)**(1/2)
    d_kappa_d_gamma = (k0 / (2 * gamma * (prefix_u**2*gamma + 1))) -\
    (k0 * np.arctan(prefix_u * gamma**(1/2)) / (2 * prefix_u * gamma**(3/2)))
    return d_kappa_d_gamma

kappa_cov = []

c = 0.1
for d in np.arange(first, 1 - c, (last - first) / n):
    gam_cov = Gamma_cov((c,d), ti, [hf, zr], tt.gamma_tern, cov)
    kL_tot, gamma = tt.kL_tot_tern_gamma((c, d), eps, tt.gamma_tern, tt.gamma_tern, ti, [hf, zr])
    slope = d_kL_d_Gamma(gamma, ti, [hf, zr], [c,d])
    kap_cov = gam_cov * slope**2
    kappa_cov.append(kap_cov[0][0])

kappa_std_dev = [k**(1/2) for k in kappa_cov]