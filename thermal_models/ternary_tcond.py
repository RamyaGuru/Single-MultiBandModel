#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 09:31:33 2020

@author: ramyagurunathan

Thermal Conductivity Ternary General Methods
"""
from math import pi
import numpy as np
from scipy.optimize import curve_fit 
from itertools import combinations

'''
Constants
'''
kB = 1.38E-23
hbar = 1.054E-34
Na = 6.02E23

'''
Methods for Property Conversions
'''

def debyeT(atmV, vs):
    return (hbar /kB) * (6*pi**2 / atmV)**(1/3) * vs

def vs_from_debyeT(atmV, debyeT):
    return (kB / hbar) * (6*pi**2 / atmV)**(-1/3) * debyeT

def average_phase_velocity(vp):
    avg_vp = (1/3) * ( vp[0]**(-3) + vp[1]**(-3) + vp[2]**(-3))**(-1/3)
    return avg_vp

'''
Methods to fit the Epsilon Values of Binary Compounds
'''

'''
Binary PDScattering Model
'''

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

'''
Binary Thermal Conductivity
'''
def kL_from_gamma(gamma, propA, propB, c):
    atmV = (1-c) * propA['atmV'] + c * propB['atmV']
    vs = (1-c) * propA['vs'] + c * propB['vs']
    k0 = (1-c) * propA['k0'] + c * propB['k0']
    prefix = (6**(1/3)/2)*(pi**(5/3)/kB)*(atmV**(2/3)/vs)
    u = (prefix * gamma * k0)**(1/2)
    kL = k0*np.arctan(u)/u
    return kL

def kL_tot(c, eps, Mfunc, Rfunc, propA, propB):
    gamma = Mfunc(propA['stoich'], propA['atmMass'], propB['atmMass'], c) +\
    eps * Rfunc(propA['stoich'], propA['atmRadius'], propB['atmRadius'], c)
    kL = kL_from_gamma(gamma, propA, propB, c)
    return kL

def kL_tot_gamma(c, eps, Mfunc, Rfunc, propA, propB):
    gamma = Mfunc(propA['stoich'], propA['atmMass'], propB['atmMass'], c) +\
    eps * Rfunc(propA['stoich'], propA['atmRadius'], propB['atmRadius'], c)
    kL = kL_from_gamma(gamma, propA, propB, c)
    return kL, gamma

def fit_eps_kL(Mfunc, Rfunc, data, propA, propB):
    data = data[data[:,0].argsort()]
    eps, cov = curve_fit(lambda c, eps:\
                         kL_tot(c, eps, Mfunc, Rfunc, propA, propB), data[0,:],\
                                data[1,:], bounds = (0,np.inf))
    kL = np.zeros(100)
    i = 0
    for c in np.linspace(data[0,0], data[0, -1], 100):
        kL[i] = kL_tot(c, eps, Mfunc, Rfunc, propA, propB)
        i = i+1
    j = 0
    kL_full = np.zeros(100)
    for d in np.linspace(1e-10,9.9999999e-1,100):
        kL_full[j] = kL_tot(d, eps, Mfunc, Rfunc, propA, propB)
        j = j+1
    return eps, kL, kL_full

def run_kL(Mfunc, Rfunc, eps, propA, propB):
    kL_full = np.zeros(100)
    j = 0
    for d in np.linspace(1e-10,9.9999999e-1,100):
        kL_full[j] = kL_tot(d, eps, Mfunc, Rfunc, propA, propB)
        j = j+1
    return kL_full

def run_kL_gamma(Mfunc, Rfunc, eps, propA, propB):
    kL_full = np.zeros(100)
    gamma_full = np.zeros(100)
    j = 0
    for d in np.linspace(1e-10,9.9999999e-1,100):
        kL_full[j], gamma_full[j] = kL_tot_gamma(d, eps, Mfunc, Rfunc, propA, propB)
        j = j+1
    return kL_full, gamma_full
'''
Ternary Thermal Conductivity
'''
def gamma_tern(stoich, mass, subst : list, c : list):
    n_sub = len(subst)
    defect_conc = sum(c)
    host_conc = (1 + 1e-10) - defect_conc
    natoms = sum(stoich)
    delM2 = 0
    denom = 0
    for n in range(len(mass)):
        msite = mass[n] * (host_conc)
        for s1 in range(n_sub):
            msite = msite  + c[s1] * subst[s1][n]
        delM2 = delM2 + stoich[n]* (host_conc) * (mass[n] - msite)**2
        for s2 in range(n_sub):
            delM2 = delM2 + stoich[n]*c[s2]*(subst[s2][n] - msite)**2
        denom = denom + stoich[n]*msite  
#    print(delM2)             
    gamma = (delM2/natoms)/((denom/natoms)**2)
    return gamma  

def kL_from_gamma_tern(gamma, propA, propB : list, c : list):
    defect_conc = sum(c)
    atmV = (1 - defect_conc) * propA['atmV'] +  sum([c[i] * propB[i]['atmV'] for i in range(len(c))]) #hard-coded for now... need to change this!
    vs = (1 - defect_conc) * propA['vs'] +  sum([c[i] * propB[i]['vs'] for i in range(len(c))])
    k0 = (1 - defect_conc) * propA['k0'] +  sum([c[i] * propB[i]['k0'] for i in range(len(c))])
    prefix = (6**(1/3)/2)*(pi**(5/3)/kB)*(atmV**(2/3)/vs)
    u = (prefix * gamma * k0)**(1/2)
    kL = k0*np.arctan(u)/u
    return kL   

def kL_tot_tern(C, eps, Mfunc, Rfunc, propA, propB : list): 
#    c1, c2 = C
    gamma = Mfunc(propA['stoich'], propA['atmMass'], [p['atmMass'] for p in propB], list(C)) +\
    eps * Rfunc(propA['stoich'], propA['atmRadius'], [p['atmRadius'] for p in propB], list(C))
    kL = kL_from_gamma_tern(gamma, propA, propB, C)
    try:
        kL = kL[0]
    except:
        pass
    return kL  

def kL_tot_tern_gamma(C, eps, Mfunc, Rfunc, propA, propB : list): 
#    c1, c2 = C
    gammaM =  Mfunc(propA['stoich'], propA['atmMass'], [p['atmMass'] for p in propB], list(C))
    gammaS = eps * Rfunc(propA['stoich'], propA['atmRadius'], [p['atmRadius'] for p in propB], list(C))
    gamma = gammaM + gammaS
    kL = kL_from_gamma_tern(gamma, propA, propB, C)
    try:
        kL = kL[0]
    except:
        pass
    return kL, gamma, gammaM, gammaS


'''
Ternary Thermal Conductivity: Muggianu Model
Think I would need to apply this to the Gamma calculation

--> only working for ternary at the moment
'''

def muggianu_model_gamma(comp, eps_list, Mfunc, Rfunc, propA, propB : list):
    '''
    propA: zr, c[0]
    
    propB[0] : hf, c[1]
    
    propB[1] : ti, c[2]
    '''
    c = comp.copy()
    c.insert(0, (1 + 1e-10) - sum(c))
    bin1 = [((1 + c[0] - c[1]) / 2) , ((1 + c[1] - c[0]) / 2)] # don't really need the second value in the pair, just 1 - the other
    bin2 = [((1 + c[0] - c[2]) / 2) , ((1 + c[2] - c[0]) / 2)]
    bin3 = [((1 + c[1] - c[2]) / 2) , ((1 + c[2] - c[1]) / 2)]
    #should probably apply this to Gamma?
    gamma1 = Mfunc(propA['stoich'], propA['atmMass'], propB[0]['atmMass'], bin1[1]) +\
    eps_list[0] * Rfunc(propA['stoich'], propA['atmRadius'], propB[0]['atmRadius'], bin1[1]) #defect conc  would be the second value in the binary??
    gamma2 = Mfunc(propA['stoich'], propA['atmMass'], propB[1]['atmMass'], bin2[1]) +\
    eps_list[1] * Rfunc(propA['stoich'], propA['atmRadius'], propB[1]['atmRadius'], bin2[1])
    gamma3 = Mfunc(propA['stoich'], propB[0]['atmMass'], propB[1]['atmMass'], bin3[1]) +\
    eps_list[2] * Rfunc(propA['stoich'], propB[0]['atmRadius'], propB[1]['atmRadius'], bin3[1])
    comp_coeff = []
    for j,k in zip([0,0,1],[1,2,2]):
        comp_coeff.append((4 * c[j] * c[k]) / ((1 + c[j] - c[k]) * (1 + c[k] - c[j])))
    gamma_tern = (gamma3 * comp_coeff[2] + gamma2 * comp_coeff[1] + gamma1 * comp_coeff[0])
    return gamma_tern

def quat_muggianu_model_gamma(comp, eps_list, Mfunc, Rfunc, propA, propB : list):
    '''
    c : Quaternary composition (%SnSe, %PbSe, %SnTe, %PbTe)
    '''
    c = list(comp.copy())
    c.insert(0, (1 + 1e-10) - sum(c))
    k = 0
    quat_gamma = 0
    quat_gammaM = 0
    prop_list = [propA, propB[0], propB[1], propB[2]]
    for i, j in zip([0,0,0,1,1,2], [1,2,3,2,3,3]): #should be easy to genneralize with a permutation generator
        coeffX = (4 * c[i] * c[j]) / ((1 + c[i] - c[j]) * (1 + c[j] - c[i]))
        binX = [((1 + c[i] - c[j]) / 2) , ((1 + c[j] - c[i]) / 2)]
        gammaM = Mfunc(prop_list[i]['stoich'], prop_list[i]['atmMass'], prop_list[j]['atmMass'], binX[1])
        gammaX = Mfunc(prop_list[i]['stoich'], prop_list[i]['atmMass'], prop_list[j]['atmMass'], binX[1]) +\
        eps_list[k] * Rfunc(prop_list[i]['stoich'], prop_list[i]['atmRadius'], prop_list[j]['atmRadius'], binX[1])
        quat_gamma = quat_gamma + coeffX * gammaX
        quat_gammaM = quat_gammaM + coeffX * gammaM
        k = k+1
    return quat_gamma, quat_gammaM
    

#def muggianu_model_gamma_general(comp, eps_list, Mfunc, Rfunc, propA, propB : list):
#    c.insert(0, (1 + 1e-10) - sum(c))
#    bins = []
#    for pair in combinations(range(len(comp+1)), 2):
#        bins.append((1 + pair[0] - pair[1]) / 2)
#    '''
#    NEED TO FINISH THIS FUNCTION
#    '''

def kL_tot_tern_muggianu(c : list, eps_list, Mfunc, Rfunc, propA, propB : list):
    gamma_tern = muggianu_model_gamma(c, eps_list, Mfunc, Rfunc, propA, propB)
    kL = kL_from_gamma_tern(gamma_tern, propA, propB, c)
    return kL, gamma_tern

def kL_tot_quat_muggianu(c : list, eps_list, Mfunc, Rfunc, propA, propB : list):
    gamma_quat, gammaM_quat = quat_muggianu_model_gamma(c, eps_list, Mfunc, Rfunc, propA, propB)
    kL = kL_from_gamma_tern(gamma_quat, propA, propB, c)
    return kL, gamma_quat, gammaM_quat


def fit_eps_kL_tern(Mfunc, Rfunc, Tt, At, propA, propB, p0 = 0):
    eps, cov = curve_fit(lambda c, eps:\
                         kL_tot_tern(c, eps, Mfunc, Rfunc, propA, propB), Tt,\
                                At, p0 = p0, bounds = (0,np.inf), method = 'dogbox')
    return eps

def fit_eps_cov_kL_tern(Mfunc, Rfunc, Tt, At, propA, propB, p0 = 0):
    eps, cov = curve_fit(lambda c, eps:\
                         kL_tot_tern(c, eps, Mfunc, Rfunc, propA, propB), Tt,\
                                At, p0 = p0, bounds = (0,np.inf))
    return eps, cov

def d_kL_d_Gamma(gamma, propA, propB : list, c : list):
    defect_conc = sum(c)
    host_conc = (1 + 1e-10) - defect_conc
    atmV = (host_conc) * propA['atmV'] +  sum(c[i] * propB[i]['atmV'] for i in range(len(c)))
    vs = (host_conc) * propA['vs'] +  sum(c[i] * propB[i]['vs'] for i in range(len(c)))
    k0 = (host_conc) * propA['k0'] +  sum(c[i] * propB[i]['k0'] for i in range(len(c)))
    prefix_u = ((6**(1/3)/2)*(pi**(5/3)/kB)*(atmV**(2/3)/vs) * k0)**(1/2)
    d_kappa_d_gamma = (k0 / (2 * gamma * (prefix_u**2*gamma + 1))) -\
    (k0 * np.arctan(prefix_u * gamma**(1/2)) / (2 * prefix_u * gamma**(3/2)))
    return d_kappa_d_gamma

def Gamma_cov(C, propA, propB, Rfunc, eps_cov):
    Gamma_R = Rfunc(propA['stoich'], propA['atmRadius'], [p['atmRadius'] for p in propB], list(C))
    return eps_cov * Gamma_R

def kL_cov_tern(C, propA, propB : list, Rfunc, gamma, eps_cov):
    gam_cov = Gamma_cov(C, propA, propB, Rfunc, eps_cov)
    slope = d_kL_d_Gamma(gamma, propA, propB, C)
    kappa_cov = gam_cov * slope**2 
    return kappa_cov[0][0]  
       

#Currently only works for tenrary
def run_kL_tern(Mfunc, Rfunc, eps, propA, propB : list):
    kL_full = np.zeros(10000)
#    n_sub = len(propB)
    k = 0
    t = []
    u = []
    v = []
    for c in np.linspace(1e-10,9.9999999e-1,100):
        for d in np.linspace(1e-10, 9.9999999e-1 - c, 100):
            t.append(c)
            u.append(d)
            v.append((1 + 1e-10)-c-d)
            kL_full[k] = kL_tot_tern([c,d], eps, Mfunc, Rfunc, propA, propB)
            k = k+1
    return kL_full, [t,u,v]


def run_kL_tern_data_dict(Mfunc, Rfunc, eps, n, propA, propB : list):
    kL_full = dict()
    first = 1e-10
    last = 9.9999999999999e-1
    j = 0
    for c in np.arange(first,1, (last - first) / n):
        k = 0
        for d in np.arange(first, 1 - c, (last - first) / n):
            kL_full[(c*100,d*100)] = kL_tot_tern((c,d), eps, Mfunc, Rfunc, propA, propB)
            k = k+1
        j = j+1
    return kL_full

def run_kL_gamma_tern_data_dict(Mfunc, Rfunc, eps, n, propA, propB : list):
    kL_full = dict()
    gamma_full = dict()
    gammaM = dict()
    gammaV = dict()
    first = 1e-10
    last = 9.9999999999999e-1
    j = 0
    for c in np.arange(first,1, (last - first) / n):
        k = 0
        for d in np.arange(first, 1 - c, (last - first) / n):
            kL_full[(c*100,d*100)], gamma_full[(c*100, d*100)],gammaM[(c*100, d*100)] ,gammaV[(c*100, d*100)]  = kL_tot_tern_gamma((c,d), eps, Mfunc, Rfunc, propA, propB)
            k = k+1
        j = j+1
    return kL_full, gamma_full, gammaM, gammaV

'''
QUATERNARY METHOD
'''

def run_kL_gamma_quat_data_dict(Mfunc, Rfunc, eps_quat, n, propA, propB : list):
    '''
    x : concentration of Se
    y : concentration of Sn
    '''
    kL_full = np.zeros([n+1, n+1])
    gamma_full = np.zeros([n+1, n+1])
    gammaM = np.zeros([n+1, n+1])
    gammaV = np.zeros([n+1, n+1])
    first = 1e-10
    last = 9.999999999e-01
    i = 0
    for x in np.arange(first, 1, (last-first)/n):
        j = 0
        for y in np.arange(first, 1, (last- first)/n):
            c_list = np.array([x * y, x * (1-y) , (1-x) * y])  #(SnSe, PbSe, SnTe)
#            eps_quat = sum(eps * c_list)
            kL_full[j,i], gamma_full[j,i],gammaM[j,i]\
            ,gammaV[j,i]  = kL_tot_tern_gamma(c_list, eps_quat, Mfunc, Rfunc, propA, propB)
            j = j+1
        i = i+1
    return kL_full, gamma_full, gammaM, gammaV

def run_kL_gamma_quat_data_dict_mugg(Mfunc, Rfunc, eps_list, n, propA, propB : list):
    '''
    x : concentration of Sn
    y : concentration of Se?
    '''
    kL_full = np.zeros([n+1, n+1])
    gamma_full = np.zeros([n+1, n+1])
    gammaM_full = np.zeros([n+1, n+1])
    first = 1e-10
    last = 9.999999999e-01
    i = 0
    for x in np.arange(first, 1, (last-first)/n):
        j = 0
        for y in np.arange(first, 1, (last- first)/n):
            c_list = np.array([x * y, x * (1-y) , (1-x) * y]) 
            kL_full[j,i], gamma_full[j,i], gammaM_full[j,i]  = kL_tot_quat_muggianu(c_list, eps_list, Mfunc, Rfunc, propA, propB)
            j = j+1
        i = i+1
    return kL_full, gamma_full, gammaM_full


def run_kL_gamma_tern_data_dict_dop(Mfunc, Rfunc, eps, n, propA, propB : list, dopant_conc = 0):
    kL_full = dict()
    gamma_full = dict()
    gammaM = dict()
    gammaV = dict()
    first = 1e-10
    last = 9.999999e-1
    weight = 1- dopant_conc
    last = last * weight
    j = 0
    for c in np.arange(first,1 * weight, (last - first) / n):
        k = 0
        for d in np.arange(first, 1 * weight - c, (last - first) / n):
            kL_full[(c*100,d*100)], gamma_full[(c*100, d*100)],  gammaM[(c*100, d*100)],  gammaV[(c*100, d*100)]  = kL_tot_tern_gamma([c,d,dopant_conc - 1e-10], eps, Mfunc, Rfunc, propA, propB)
            k = k+1
        j = j+1
    return kL_full, gamma_full, gammaM, gammaV

def run_kL_cov_tern_data_dict(Mfunc, Rfunc, eps, eps_cov, n, propA, propB : list):
    kL_full = dict()
    kL_std_dev = dict()
    gamma = dict()
    gammaM = 0
    gammaS = 0
    first = 1e-10
    last = 9.9999999999999e-1
    j = 0
    for c in np.arange(first,1, (last - first) / n):
        k = 0
        for d in np.arange(first, 1 - c, (last - first) / n):
            kL_full[(c*100,d*100)], gamma[(c*100, d*100)], gammaM,  gammaS = kL_tot_tern_gamma((c,d), eps, Mfunc, Rfunc, propA, propB)
#            kL_std_dev[(c*100,d*100)] = (kL_cov_tern((c,d), propA, propB, Rfunc, gamma, eps_cov))**(1/2)
            k = k+1
        j = j+1
    return kL_full, gamma


def run_kL_tern_data_dict_muggianu(Mfunc, Rfunc, eps, n, propA, propB : list):
    kL_full = dict()
    gamma_full = dict()
    first = 1e-10
    last = 9.9999999999999e-1
    j = 0
    for c in np.arange(first,1,(last - first) / n):
        k = 0
        for d in np.arange(first, 1 - c, (last - first) / n):
            kL_full[(c*100,d*100)],  gamma_full[(c*100, d*100)] = kL_tot_tern_muggianu([c,d], eps, Mfunc, Rfunc, propA, propB)
            k = k+1
        j = j+1
    return kL_full, gamma_full


