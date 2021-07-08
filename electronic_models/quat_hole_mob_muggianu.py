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
import numpy as np

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
PbSe_SnSe = np.array([c_list[:5], [1098.34, 262.35, 413.88, 451.53, 250.26]])
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
full_quat = np.loadtxt('../datafiles/quat_chalcogenides.txt', delimiter = ',')

binary_mob_list = [PbSe_SnSe, SnTe_SnSe, PbTe_SnSe, SnTe_PbSe, PbTe_PbSe, PbTe_SnTe]

binary_inv_mob_list = [PbSe_SnSe_i, SnTe_SnSe_i, PbTe_SnSe_i, SnTe_PbSe_i, PbTe_PbSe_i, PbTe_SnTe_i]


def write_data_file(quat_rho, quat_comp, out_file):
    '''
    Write out easy to read data file for either the experimental data or the alloy scattering model predictions
    
    quat_rho: 2D array of experimental or modelled values
    '''
    rho_new = quat_rho.flatten()
    rho_new = np.array([rho_new]).T
    comp_new = quat_comp.reshape(121, 4)
    full_array = np.append(comp_new, rho_new, axis = 1)
    np.savetxt(out_file, full_array, delimiter = ',', header = 'x_SnSe, x_PbSe, x_SnTe, x_PbTe, Hole_Mobility (cm2/Vs)')


if __name__ == '__main__':
    C_list = qrm.fit_Nord_coeff(c_list, binary_inv_mob_list)
    C_list[0] = 0
    quat_inv_mob, comp_array = qrm.run_quat_rho(qrm.Nordheim_exc_resistivity,\
                                    qrm.Nordheim_vegard_resistivity, C_list, binary_inv_mob_list, 10)
    
    quat_dev_alloy = full_quat-(1/quat_inv_mob)
    
    dev_dos_effect = qrm.deviation_from_alloy_model(full_quat, (1/quat_inv_mob))
    
    '''
    Replace multiphase or Pnma samples with NaNs
    '''
    for r in range(6):
        quat_inv_mob[5+r][-(3+r):] = math.nan
    fig,ax = plt.subplots()
    plt.matshow(1/quat_inv_mob, cmap = plt.cm.get_cmap('rainbow'), vmin = 0, vmax = 800)
    plt.title('Alloy Model')
    cbar = plt.colorbar()    
    cbar.set_label(r'Hall Hole Mobility (cm$^2$/V/s)')
    plt.savefig('quat_hole_mob.pdf', bbox_inches = 'tight')    
    

    '''
    Plot the experimental data
    '''
    for r in range(6):
        full_quat[5+r][-(3+r):] = math.nan
    fig,ax = plt.subplots()
    plt.matshow(full_quat, cmap = plt.cm.get_cmap('rainbow'), vmin = 0, vmax = 800)
    plt.title('Experimental')
    cbar = plt.colorbar()    
    cbar.set_label(r'Hall Hole Mobility (cm$^2$/V/s)')
    plt.savefig('quat_hole_mob.pdf', bbox_inches = 'tight')    
       
    
    '''
    Plot the deviation from alloy scattering: Experimental - Prediction
    '''
    for r in range(6):
        quat_dev_alloy[5+r][-(3+r):] = math.nan
    fig,ax = plt.subplots()
    plt.matshow(quat_dev_alloy, cmap = plt.cm.get_cmap('rainbow'))
    plt.title('Deviations from Alloy Model')
    cbar = plt.colorbar()     
    
    '''
    Plot the deviation from alloy scattering: Experimental - Prediction
    '''
    for r in range(6):
        dev_dos_effect[5+r][-(3+r):] = math.nan
    fig,ax = plt.subplots()
    plt.matshow(dev_dos_effect, cmap = plt.cm.get_cmap('rainbow'))
    plt.title('Deviations from Alloy Model')
    cbar = plt.colorbar()  
    
    '''
    Print Binaries
    '''
    mpl.rcParams['figure.figsize'] = [5,3]
    
    plt.figure()
    plt.scatter(PbTe_PbSe_i[0],PbTe_PbSe[1])
    bin4_model = qrm.Nordheim_resistivity_model(PbTe_PbSe_i[0], C_list[4], binary_inv_mob_list[4])
    plt.plot(PbTe_PbSe[0], 1/bin4_model)
    plt.ylabel(r'Hall Hole Mobility (cm$^2$/V/s)')
    plt.xlabel('%Se')
    
    
    plt.figure()
    plt.scatter(SnTe_SnSe_i[0], SnTe_SnSe[1])
    bin2_model = qrm.Nordheim_resistivity_model(SnTe_SnSe_i[0], C_list[1], binary_inv_mob_list[1])
    plt.plot(SnTe_SnSe[0], 1/ bin2_model) 
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


    '''
    Write Datafile
    '''
    exp_outfile = '../datafiles/quat_exp_mobility_full.csv'
    write_data_file(full_quat, comp_array, exp_outfile)
    
    model_outfile = '../datafiles/quat_model_mobility_full.csv'
    write_data_file(1/quat_inv_mob, comp_array, model_outfile)
