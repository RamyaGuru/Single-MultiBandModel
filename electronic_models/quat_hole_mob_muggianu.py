#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 16:53:32 2021

@author: ramyagurunathan

Quaternary Hole Mobility Model
"""

import numpy as np
import quat_resistivity_muggianu as qrm
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
import matplotlib as mpl

'''
Data from Ortiz 2019 (PbTe-PbSe-SnTe-SnSe Alloys)

indexes = [SnSe, PbSe, SnTe, PbTe]
'''

c_list = list(np.arange(0,1.1,0.1))

# c is the %PbSe
PbTe_PbSe = np.array([c_list, [558.93, 495.13, 563.45, 495.98, 504.28, 529.51,\
                               522.19, 631.78, 594.66, 604.31, 1098.34 ]])
PbTe_PbSe_i = np.array([PbTe_PbSe[0], 1/PbTe_PbSe[1]])
# c is the %SnTe
PbTe_SnTe = np.array([c_list, [558.93, 438.49, 516.30, 445.63, 327.65, 118.70,\
                               166.06, 92.37, 106.53, 117.10, 116.02]])
PbTe_SnTe_i = np.array([PbTe_SnTe[0], 1/PbTe_SnTe[1]])
# c is the %SnSe
SnTe_SnSe = np.array([c_list[:3], [116.02, 97.90, 92.03]])#, 65.41, 32.05, 12.37, 2.22,\
                               #14.37, 28.08, 24.54, 31.53]])
SnTe_SnSe_i = np.array([SnTe_SnSe[0], 1/SnTe_SnSe[1]])
# c is the %SnSe
PbSe_SnSe = np.array([c_list, [1098.34, 262.35, 413.88, 451.53, 250.26, 162.10,\
                               42.14, 1.86, 13.45, 12.72, 31.53]])
PbSe_SnSe_i = np.array([PbSe_SnSe[0], 1/PbSe_SnSe[1]])

# c is the %PbSe
SnTe_PbSe = np.array([c_list, [116.02, 85.22, 73.35, 75.29, 104.72, 156.19,\
                               251.65, 379.70, 449.06, 464.30, 1098.34]])
SnTe_PbSe_i = np.array([SnTe_PbSe[0], 1/SnTe_PbSe[1]])

# c is the %SnSe
PbTe_SnSe = np.array([c_list, [558.93, 377.00, 413.00, 357.28, 246.85, 156.19,\
                               82.52, 42.14, 5.17, 15.18, 31.53]])
PbTe_SnSe_i = np.array([PbTe_SnSe[0], 1/PbTe_SnSe[1]])

#Full quaternary data

binary_mob_list = [PbSe_SnSe, SnTe_SnSe, PbTe_SnSe, SnTe_PbSe, PbTe_PbSe, PbTe_SnTe]

binary_inv_mob_list = [PbSe_SnSe_i, SnTe_SnSe_i, PbTe_SnSe_i, SnTe_PbSe_i, PbTe_PbSe_i, PbTe_SnTe_i]



if __name__ == '__main__':
    C_list = qrm.fit_Nord_coeff(c_list, binary_inv_mob_list)
    C_list[0] = 0
    quat_inv_mob = qrm.run_quat_rho(qrm.Nordheim_exc_resistivity,\
                                    qrm.Nordheim_vegard_resistivity, C_list, binary_inv_mob_list, 10)

    '''
    Replace multiphase or Pnma samples with NaNs
    '''
    for r in range(6):
        quat_inv_mob[5+r][-(3+r):] = math.nan
    fig,ax = plt.subplots()
    plt.matshow(1/quat_inv_mob, cmap = plt.cm.get_cmap('rainbow'), vmin = 0, vmax = 800)
    cbar = plt.colorbar()    
    cbar.set_label(r'Hall Hole Mobility (cm$^2$/V/s)')
    plt.savefig('quat_hole_mob.pdf', bbox_inches = 'tight')    
    
    mpl.rcParams['figure.figsize'] = [5,3]
    
    
    '''
    Print Binaries
    '''
    plt.figure()
    plt.scatter(PbTe_PbSe_i[0],PbTe_PbSe[1])
    bin4_model = qrm.Nordheim_resistivity_model(PbTe_PbSe_i[0], C_list[4], binary_inv_mob_list[4])
    plt.plot(PbTe_PbSe[0], 1/bin4_model)
    plt.ylabel(r'Hall Hole Mobility (cm$^2$/V/s)')
    plt.xlabel('%Se')
    
    
    plt.figure()
    plt.scatter(SnTe_SnSe_i[0], SnTe_SnSe[1])
    bin2_model = qrm.Nordheim_resistivity_model(SnTe_SnSe_i[0], C_list[1], binary_inv_mob_list[1])
    plt.plot(PbTe_PbSe[0], 1/ bin2_model) 
    plt.ylabel(r'Hall Hole Mobility (cm$^2$/V/s)')
    plt.xlabel('%Se')
    
    plt.figure()
    plt.scatter(PbSe_SnSe_i[0], PbSe_SnSe[1])
    bin0_model = qrm.Nordheim_resistivity_model(PbSe_SnSe_i[0], C_list[0], binary_inv_mob_list[0])
    plt.plot(PbSe_SnSe[0], 1/bin0_model)  
    plt.ylabel(r'Hall Hole Mobility (cm$^2$/V/s)')
    plt.xlabel('%Sn')
    
    
    plt.figure()
    plt.scatter(PbTe_SnTe_i[0], PbTe_SnTe[1])
    bin5_model = qrm.Nordheim_resistivity_model(PbTe_SnTe_i[0], C_list[5], binary_inv_mob_list[5])
    plt.plot(PbTe_SnTe[0], 1/bin5_model)
    plt.ylabel(r'Hall Hole Mobility (cm$^2$/V/s)')
    plt.xlabel('%Sn')
