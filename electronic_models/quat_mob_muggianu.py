#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 14:48:53 2021

@author: ramyagurunathan

Quaternary Intrinsic Mobility Model!!
"""


import numpy as np
import quat_resistivity_muggianu as qrm
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

'''
Data from Ortiz 2019 (PbTe-PbSe-SnTe-SnSe Alloys)

indexes = [SnSe, PbSe, SnTe, PbTe]
'''

c_list = list(np.arange(0,1.1,0.1))

# c is the %PbSe
PbTe_PbSe = np.array([c_list, [616.50, 524.08, 586.37, 535.90, 532.11, 567.62,\
                               558.68, 676.73, 643.78, 661.31, 1224.12]])
PbTe_PbSe_i = np.array([PbTe_PbSe[0], 1/PbTe_PbSe[1]])
# c is the %SnTe
PbTe_SnTe = np.array([c_list, [616.50, 531.85, 760.61, 802.79, 845.53, 544.83,\
                               518.12, 250.52, 282.46, 300.42, 283.64]])
PbTe_SnTe_i = np.array([PbTe_SnTe[0], 1/PbTe_SnTe[1]])
# c is the %SnSe
SnTe_SnSe = np.array([c_list, [283.64, 267.70, 250.82, 159.69, 73.05 ,26.72,\
                               3.71, 14.83, 28.50, 24.71, 31.74]])
SnTe_SnSe_i = np.array([SnTe_SnSe[0], 1/SnTe_SnSe[1]])
# c is the %SnSe
PbSe_SnSe = np.array([[ c for c in c_list if c not in [0.1, 0.2]], [1224.12, 812.22, 529.17, 346.07,\
                               97.10, 3.06, 13.50, 12.79, 31.74]])
PbSe_SnSe_i = np.array([PbSe_SnSe[0], 1/PbSe_SnSe[1]])

# c is the %PbSe
SnTe_PbSe = np.array([c_list, [283.64, 238.09, 215.05, 238.47, 352.83, 470.67,\
                               575.34, 661.49, 608.29, 537.69, 1224.12]])
SnTe_PbSe_i = np.array([SnTe_PbSe[0], 1/SnTe_PbSe[1]])

# c is the %SnSe
PbTe_SnSe = np.array([c_list, [616.50, 441.94, 585.51, 630.75, 584.67, 470.67,\
                               268.52, 124.39, 12.87, 15.26, 31.74]])
PbTe_SnSe_i = np.array([PbTe_SnSe[0], 1/PbTe_SnSe[1]])


binary_mob_list = [PbSe_SnSe, SnTe_SnSe, PbTe_SnSe, SnTe_PbSe, PbTe_PbSe, PbTe_SnTe]

binary_inv_mob_list = [PbSe_SnSe_i, SnTe_SnSe_i, PbTe_SnSe_i, SnTe_PbSe_i, PbTe_PbSe_i, PbTe_SnTe_i]


if __name__ == '__main__':
    C_list = qrm.fit_Nord_coeff(c_list, binary_inv_mob_list)
    C_list[0] = 0
    quat_inv_mob = qrm.run_quat_rho(qrm.Nordheim_exc_resistivity,\
                                    qrm.Nordheim_vegard_resistivity, C_list, binary_inv_mob_list, 10)
    
    fig,ax = plt.subplots()
    plt.matshow(1/quat_inv_mob, cmap = plt.cm.get_cmap('rainbow'), vmin = 0, vmax = 800)
    cbar = plt.colorbar()        
    
    
    '''
    Print Binaries
    '''
    plt.figure()
    plt.scatter(PbTe_PbSe_i[0],1/ PbTe_PbSe_i[1])
    bin4_model = qrm.Nordheim_resistivity_model(PbTe_PbSe_i[0], C_list[4], binary_inv_mob_list[4])
    plt.plot(PbTe_PbSe[0], 1/bin4_model)
    
    
    plt.figure()
    plt.scatter(SnTe_SnSe_i[0], 1/ SnTe_SnSe_i[1])
    bin2_model = qrm.Nordheim_resistivity_model(SnTe_SnSe_i[0], C_list[1], binary_inv_mob_list[1])
    plt.plot(PbTe_PbSe[0], 1/ bin2_model) 
    
    plt.figure()
    plt.scatter(PbSe_SnSe_i[0], 1/PbSe_SnSe_i[1])
    bin0_model = qrm.Nordheim_resistivity_model(PbSe_SnSe_i[0], C_list[0], binary_inv_mob_list[0])
    plt.plot(PbSe_SnSe[0], 1/bin0_model)    
    
    
    plt.figure()
    plt.scatter(PbTe_SnTe_i[0], 1/PbTe_SnTe_i[1])
    bin5_model = qrm.Nordheim_resistivity_model(PbTe_SnTe_i[0], C_list[5], binary_inv_mob_list[5])
    plt.plot(PbTe_SnTe[0], 1/bin5_model)