#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 17:52:40 2021

@author: ramyagurunathan

sigmae0 modle CALPHAD mixing strategies


To Do Today:
    --> implement a verson of the Redich Kister polynomial
"""

'''
constants
'''
e = 1.602e-19
m_e = 9.109e-31
hbar = 1.045e-34
kb = 1.38e-23


import sigmae0_alloy as sig
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib as mpl
import math

mpl.rcParams['figure.figsize'] = [5,3] 


def split_binary_data(sig_df):
    '''
    Inputs:
        c : array of molar fractions
    Note: fits separate U parameters for each binary. 
    They don't change with composition
    
    bin_list = [B-C, A-C, A-B]
    '''
    #Split up data for each binary system
    bin_list = []
    for comp in sigmae0_df.columns[1:4]:
        binary = sigmae0_df.loc[sigmae0_df[comp] == 0.0]
        bin_list.append(binary)                
    return bin_list


def fit_binary_U(bin_list, A_endpt, B_endpt):
    '''
    Inputs:
        sigmae0_df: dataframe with ternary compositions and sigma_e0 values
    Output:
        U_list: []
    '''
    U_list = []
    for b in bin_list:
        u, cov = sig.fit_U(b, A_endpt, B_endpt, mat_props, 'n', T = 300, p0 = 0)
        U_list.append(u)    
    return U_list


def trivial_U_mixing(sig_df, A_endpt, B_endpt, mat_props, dpg, c : list, T = 300):
    '''
    Calculate alloy scattering parameters using the correct average of binary U parameters?
    
    U_list = [B-C, A-C, A-B]
    '''
    bin_list = split_binary_data(sig_df)
    U_list = fit_binary_U(bin_list, A_endpt, B_endpt)
    defect_conc = sum(c)
    c_2 = (1 + 1e-10) - defect_conc
    m_eff = mat_props['eff_mass_{}'.format(dpg)][0] * c[0] +\
    mat_props['eff_mass_{}'.format(dpg)][1] * c[1] +\
    mat_props['eff_mass_{}'.format(dpg)][2] * c_2
    
    Cl = mat_props['BulkMod'][0] * c[0] + mat_props['BulkMod'][1] * c[1] +\
    mat_props['BulkMod'][2] * c_2
    
    coeff = mat_props['coeff_{}'.format(dpg)][0] * c[0] + mat_props['coeff_{}'.format(dpg)][1] * c[1] +\
    mat_props['coeff_{}'.format(dpg)][2] * c_2
    
    AtmV = mat_props['AtmV'][0] * c[0] + mat_props['AtmV'][1] * c[1] +\
    mat_props['AtmV'][2] * c_2
    
    Nv = mat_props['Nv'][2][dpg]
    
    U_tern = U_list[0] * c[1] * c_2 + U_list[1] * c[0] * c_2 + U_list[2] * c[0] * c[1]
    #Does weighting coefficient for U change depending on the 
    sigmae0 = (2 * hbar * Cl * Nv * e**2 * coeff) / (3 * math.pi * m_eff * m_e * e**2 *\
               (1 + U_tern * (3 * math.pi**2 * Cl * AtmV)/\
                (8 * kb * T)))
    return sigmae0


    

def redlich_kister_polynomial(sig_df):
    '''
    Inputs:
        c : array of molar fractions
        Lcoeff: array of coefficients in U polynomial --> how to fit these?
    NOTE: calculates based on the U parameters fit to each of the binaries
    '''
    
    return 
#    

def muggianu_model(sigmae0_df, A_endpt, B_endpt, mat_props, dpg, c : list, T = 300):
    '''
    Inputs:
        c : array of molar fractions (i,j) k = 1-i-j
        U : array of pairwise alloy scattering U coefficients
        Will these be symmetric such that U_ij = U_ji?
        Kind of depends on what the "host" bandstructure is?
    NOTE: calculates based on the binary sigma_e0 value. Chooses binary composition
    based on the geometry of the Muggianu diagram
    U_list = [B-C, A-C, A-B]
    '''
    bin_list = split_binary_data(sig_df)
    U_list = fit_binary_U(bin_list, A_endpt, B_endpt)
    c.append((1 + 1e-10) - sum(c))
    bin3 = [((1 + c[0] - c[1]) / 2) , ((1 + c[1] - c[0]) / 2)]
    bin2 = [((1 + c[0] - c[2]) / 2) , ((1 + c[2] - c[0]) / 2)]
    bin1 = [((1 + c[1] - c[2]) / 2) , ((1 + c[2] - c[1]) / 2)]
    sig3 = sig.sigma_e0_model(bin3, U_list[2], mat_props, dpg, T)
    sig2 = sig.sigma_e0_model(bin2, U_list[1], mat_props, dpg, T)
    sig1 = sig.sigma_e0_model(bin1, U_list[0], mat_props, dpg, T)  
    comp_coeff = []
    for j,k in zip([1,0,0],[2,2,1]):
        comp_coeff.append((4 * c[j] * c[k]) / ((1 + c[j] - c[k]) * (1 + c[k] - c[j])))
    sigma_tern = sig3 * comp_coeff[2] + sig2 * comp_coeff[1] + sig1 * comp_coeff[0]
    return sigma_tern


if __name__ == '__main__':
    mat_props = {'formula': ['HfNiSn', 'TiNiSn', 'ZrNiSn'], 'BulkMod': [105E9, 122E9, 152E9],\
         'mpid': ['mp-20523', 'mp-924130', 'mp-30806'],\
         'AtmV' : [1.91E-29, 1.755E-29, 1.93E-29]}

#Add the valley degeneracy
    X_Nv = {'p' : 1, 'n' : 3}
    Y_Nv = {'p' : 1, 'n' : 3} #mtminer_valley_degeneracy('mp-924130')
    Z_Nv = {'p' : 1, 'n' : 3}
    

    mat_props['Nv'] = [X_Nv, Y_Nv, Z_Nv]

    #Source: PCCP 2016, 18 14107
    mat_props['eff_mass_n'] = [0.36, 0.56, 0.38]
    mat_props['eff_mass_p'] = [0.22, 0.39, 0.25]


    '''
    Get sigmae0_df
    '''
    data_path = '/Users/ramyagurunathan/Documents/PhDProjects/Argonne_TECCA/SigmaE0_values/'
    datafile = data_path + 'XNiSn_gauss_data_purple.csv'
    sigmae0_df = sig.fetch_sigmae0_dataframe(datafile)
    #Drop multiphase samples
#    sigmae0_df = sigmae0_df.drop([8, 29, 36])
    x_val, y_val, z_val = sig.fetch_endmember_values(sigmae0_df)
    
#    datafile2 = data_path + 'XNiSn_S_and_sigma.csv'
#    sige0_list = sig.fit_sigmae0(datafile2)
#    for item in sige0_list:
#        sigmae0_df = sigmae0_df.append(item.iloc[0])
#        
    
    '''
    Fit coefficients
    '''
    mat_props = sig.fit_sigma_e0_coefficient([x_val, y_val, z_val], mat_props, 'n', T = 300)
    
    '''
    Split into binary data
    '''
    bin_list = split_binary_data(sigmae0_df)
    
    '''
    Fit separate U values
    '''
    U_list = []
    for b in bin_list:
        u, cov = sig.fit_U(b, 'HfNiSn', 'TiNiSn', mat_props, 'n', T = 300, p0 = 0)
        U_list.append(u)
    #Why does the Hf-Ti binary have such a huge U parameter?

    

    