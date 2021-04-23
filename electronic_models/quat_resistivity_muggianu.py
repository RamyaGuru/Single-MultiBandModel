#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 11:24:19 2021

@author: ramyagurunathan

Quaternary Pb(Sn)-Te(Se) System
Get Resistivity using Muggianu model



Should do this on mobility instead?
"""

import numpy as np
import mixing_models_sigmae0 as mms
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib as mpl

'''
Data from Ortiz 2019 (PbTe-PbSe-SnTe-SnSe Alloys)

indexes = [SnSe, PbSe, SnTe, PbTe]
'''
c_list = list(np.arange(0,1.1,0.1))

# c is the %PbSe
PbTe_PbSe = np.array([c_list, [5.40, 13.68, 20.10, 24.67, 27.29, 22.44, 23.79,\
                               18.25, 15.27, 11.30, 3.36]])
# c is the %SnTe
PbTe_SnTe = np.array([c_list, [5.40, 2.72, 1.08, 0.71, 0.47, 0.38, 0.33, 0.25,\
                              0.21, 0.17, 0.13]])
# c is the %SnSe
SnTe_SnSe = np.array([c_list, [0.13, 0.17, 0.20, 0.26, 0.50, 1.17, 8.52, 210.23,\
                              232.24, 269.08, 319.44]])

# c is the %SnSe
PbSe_SnSe = np.array([c_list, [3.36, 9.23, 2.08, 0.96, 0.70, 1.00, 4.13, 244.84,\
                              1458.12, 639.35, 319.44]])

# c is the %PbSe
SnTe_PbSe = np.array([c_list, [0.13, 0.20, 0.28, 0.38, 0.48, 0.61, 0.80, 1.20,\
                              2.40, 7.91, 3.36]])

# c is the %SnSe
PbTe_SnSe = np.array([c_list, [5.40, 3.77, 1.54, 1.02, 0.73, 0.61, 0.51, 2.06,\
                              28.65, 430.61, 319.44]])


binary_data_list = [PbSe_SnSe, SnTe_SnSe, PbTe_SnSe, SnTe_PbSe, PbTe_PbSe, PbTe_SnTe]


def Nordheim_resistivity_model(x, C, rho_data):
    rho0 = (1-x) * rho_data[1][0] + x * rho_data[1][-1]
    return rho0 + x * (1-x) * C * rho0

def Nordheim_exc_resistivity(x, C, rho_data):
    rho0 = (1-x) * rho_data[1][0] + x * rho_data[1][-1]
    return x * (1-x) * C * rho0

def fit_Nord_coeff(x, binary_data_list):
    C_list = []
    for rho_data in binary_data_list:
#        rho0 = (1-x) * rho_data[0][0] + x * rho_data[1][-1]
        C, eps = curve_fit(lambda x, C: Nordheim_resistivity_model(x, C, rho_data), rho_data[0],\
                      rho_data[1], bounds = (0, np.inf))
        C_list.append(C[0])
    return C_list

def Nordheim_vegard_resistivity(c: list, binary_data_list):
    return sum([c[0] * binary_data_list[0][1][-1], c[1] * binary_data_list[0][1][0],\
               c[2] * binary_data_list[-1][1][-1], c[3] * binary_data_list[-1][1][0]])

def quat_muggianu_model(c: list, binary_exc_fxn, binary_vegard_fxn, binary_coeff_list, binary_data_list):
    '''
    c : Quaternary composition (%SnSe, %PbSe, %SnTe, %PbTe)
    '''
    k = 0
    quat_exc_rho = 0
    quat_vegard_rho = binary_vegard_fxn(c, binary_data_list)
    for i, j in zip([0,0,0,1,1,2], [1,2,3,2,3,3]):
        coeffX = (4 * c[i] * c[j]) / ((1 + c[i] - c[j]) * (1 + c[j] - c[i]))
        binX = [((1 + c[i] - c[j]) / 2) , ((1 + c[j] - c[i]) / 2)]
#        rho0 = c[i] * binary_data_list[k][-1] + c[j] * binary_data_list[k][0] # defect + host
        rhoX = binary_exc_fxn(binX[0], binary_coeff_list[k], binary_data_list[k])
        quat_exc_rho = quat_exc_rho + coeffX * rhoX
        k = k+1
    quat_rho = quat_vegard_rho + quat_exc_rho
    return quat_rho

def run_quat_rho(binary_exc_fxn, binary_vegard_fxn, binary_coeff_list, binary_data_list, n):
    '''
    Generate matrix that follows the heatmap from Brenden's paper
    '''
    quat_rho = np.zeros([n+1,n+1])
    first = 1e-10
    last = 9.999999999e-01
    i = 0
    for x in np.arange(first, 1, (last-first)/n):
        j = 0
        for y in np.arange(first, 1, (last- first)/n):
            c_list = [x * y, x * (1-y) , (1-x) * y, (1-x) * (1-y)] #Replace any negative values with 0? 
#            print(c_list)
            quat_rho[j,i] = quat_muggianu_model(c_list, binary_exc_fxn, binary_vegard_fxn, binary_coeff_list, binary_data_list)
#            print(quat_rho[i,j])
            j = j+1
        i = i+1
    return quat_rho


if __name__ =='__main__':
    C_list = fit_Nord_coeff(c_list, binary_data_list)
    quat_rho = run_quat_rho(Nordheim_exc_resistivity, Nordheim_vegard_resistivity,\
                            C_list, binary_data_list, 10)
    
    fig,ax = plt.subplots()
    plt.matshow(np.log10(quat_rho), cmap = plt.cm.get_cmap('rainbow'), vmin = -1, vmax = 3)
    cbar = plt.colorbar()    
    cbar.set_label(r'Log Resistivity (m$\Omega$cm)')
    plt.savefig('quat_res.pdf', bbox_inches = 'tight')      
   
    mpl.rcParams['figure.figsize'] = [5,3]
    
    '''
    Print Binaries
    '''
    plt.figure()
    plt.scatter(PbTe_PbSe[0], PbTe_PbSe[1])
    bin4_model = Nordheim_resistivity_model(PbTe_PbSe[0], C_list[4], binary_data_list[4])
    plt.plot(PbTe_PbSe[0], bin4_model)
    plt.ylabel(r'Resistivity m$\Omega$cm')
    plt.xlabel('%Se') 
    
    
    plt.figure()
    plt.scatter(SnTe_SnSe[0], SnTe_SnSe[1])
    bin2_model = Nordheim_resistivity_model(SnTe_SnSe[0], C_list[1], binary_data_list[1])
    plt.plot(PbTe_PbSe[0], bin2_model) 
    plt.ylabel(r'Resistivity m$\Omega$cm')
    plt.xlabel('%Se')
