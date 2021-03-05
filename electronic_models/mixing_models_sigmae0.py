#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 17:52:40 2021

@author: ramyagurunathan

sigmae0 modle CALPHAD mixing strategies


To Do Today:
    --> implement a verson of the Redich Kister polynomial
    
Is this potentially backwards?
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
import ternary

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
    for comp in sig_df.columns[1:4]:
        binary = sig_df.loc[sig_df[comp] == 0.0]
        bin_list.append(binary)                
    return bin_list


def convert_data_df_to_dict(sig_df, A_endpt, B_endpt):
    sig_dict = {}
    for row in sig_df.index:
        sig_dict[(sig_df[A_endpt][row] * 100, sig_df[B_endpt][row] * 100)] =  sig_df['Sigma_E0'][row]
    return sig_dict

def fit_binary_U_list(bin_list, mat_props, d_key_list =  ['TiNiSn', 'HfNiSn', 'HfNiSn'],\
                      d_h_pairs = [[1,2],[0,2],[0,1]], T = 300, p0 = 0):
    '''
    Inputs:
        sigmae0_df: dataframe with ternary compositions and sigma_e0 values
    Output:
        U_list: []
    '''
    U_list = []
    for b in [0,1,2]:
        u, cov = sig.fit_binary_U_endpts(bin_list[b], d_key_list[b],\
                                         d_h_pairs[b][0], d_h_pairs[b][1],\
                                         mat_props, 'n', T = 300, p0 = 0)
        U_list.append(u)    
    return U_list


    
#def trivial_U_mixing(sig_df, A_endpt, B_endpt, mat_props, dpg, c : list, T = 300):
#    '''
#    Calculate alloy scattering parameters using the correct average of binary U parameters?
#    
#    U_list = [B-C, A-C, A-B]
#    '''
#    bin_list = split_binary_data(sig_df)
#    U_list = fit_binary_U_list(bin_list, A_endpt, B_endpt)
#    defect_conc = sum(c)
#    c_2 = (1 + 1e-10) - defect_conc
#    m_eff = mat_props['eff_mass_{}'.format(dpg)][0] * c[0] +\
#    mat_props['eff_mass_{}'.format(dpg)][1] * c[1] +\
#    mat_props['eff_mass_{}'.format(dpg)][2] * c_2
#    
#    Cl = mat_props['BulkMod'][0] * c[0] + mat_props['BulkMod'][1] * c[1] +\
#    mat_props['BulkMod'][2] * c_2
#    
#    coeff = mat_props['coeff_{}'.format(dpg)][0] * c[0] + mat_props['coeff_{}'.format(dpg)][1] * c[1] +\
#    mat_props['coeff_{}'.format(dpg)][2] * c_2
#    
#    AtmV = mat_props['AtmV'][0] * c[0] + mat_props['AtmV'][1] * c[1] +\
#    mat_props['AtmV'][2] * c_2
#    
#    Nv = mat_props['Nv'][2][dpg]
#    
##    U_tern = 
#    #Does weighting coefficient for U change depending on the 
#    sigmae0 = (2 * hbar * Cl * Nv * e**2 * coeff) / (3 * math.pi * m_eff * m_e * e**2 *\
#               (1 + U_tern * (3 * math.pi**2 * Cl * AtmV)/\
#                (8 * kb * T)))
#    return sigmae0


    

#def redlich_kister_polynomial(sig_df):
#    '''
#    Inputs:
#        c : array of molar fractions
#        Lcoeff: array of coefficients in U polynomial --> how to fit these?
#    NOTE: calculates based on the U parameters fit to each of the binaries
#    '''
#    
#    return 
#    

def muggianu_model(sigmae0_df, bin_list, U_list, mat_props, dpg, c : list, T = 300):
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
    c.append((1 + 1e-10) - sum(c))
    bin3 = [((1 + c[0] - c[1]) / 2) , ((1 + c[1] - c[0]) / 2)]
    bin2 = [((1 + c[0] - c[2]) / 2) , ((1 + c[2] - c[0]) / 2)]
    bin1 = [((1 + c[1] - c[2]) / 2) , ((1 + c[2] - c[1]) / 2)]
    sig3 = sig.binary_sigma_e0_endpts(bin3[0], U_list[2], 0, 1, mat_props, dpg, T)
    sig2 = sig.binary_sigma_e0_endpts(bin2[0], U_list[1], 0, 2, mat_props, dpg, T)
    sig1 = sig.binary_sigma_e0_endpts(bin1[0], U_list[0], 1, 2, mat_props, dpg, T) 
    comp_coeff = []
    for j,k in zip([1,0,0],[2,2,1]):
        comp_coeff.append((4 * c[j] * c[k]) / ((1 + c[j] - c[k]) * (1 + c[k] - c[j])))
#    print(comp_coeff)
#    print([bin3, bin2, bin1])
#    print([sig3, sig2, sig1])
    sigma_tern = (sig3 * comp_coeff[2] + sig2 * comp_coeff[1] + sig1 * comp_coeff[0])
    return sigma_tern

def muggianu_model_sig_diff(sigmae0_df, bin_list, U_list, mat_props, dpg, c : list, T = 300):
    '''
    Calculate excess sigma_e0 term by substracting off the Vegard's law value

    Inputs:
        c : array of molar fractions (i,j) k = 1-i-j
        U : array of pairwise alloy scattering U coefficients
        Will these be symmetric such that U_ij = U_ji?
        Kind of depends on what the "host" bandstructure is?
    NOTE: calculates based on the binary sigma_e0 value. Chooses binary composition
    based on the geometry of the Muggianu diagram
    U_list = [B-C, A-C, A-B]
    '''
    sigma_endpts = sig.fetch_endmember_values(sigmae0_df)
    c.append((1 + 1e-10) - sum(c))
    bin3 = [((1 + c[0] - c[1]) / 2) , ((1 + c[1] - c[0]) / 2)]
    bin2 = [((1 + c[0] - c[2]) / 2) , ((1 + c[2] - c[0]) / 2)]
    bin1 = [((1 + c[1] - c[2]) / 2) , ((1 + c[2] - c[1]) / 2)]
    excess3 = -1 * (bin3[0] * (1/ sigma_endpts[0]) + bin3[1] * (1/sigma_endpts[1])) +\
    (1/ sig.binary_sigma_e0_endpts(bin3[0], U_list[2], 0, 1, mat_props, dpg, T))
    excess2 = -1 * (bin2[0] * (1/ sigma_endpts[0]) + bin2[1] * (1/ sigma_endpts[2])) +\
        (1/ sig.binary_sigma_e0_endpts(bin2[0], U_list[1], 0, 2, mat_props, dpg, T))
    excess1 = -1 * (bin1[0] * (1/ sigma_endpts[1]) + bin1[1] * (1/ sigma_endpts[2])) +\
        (1/ sig.binary_sigma_e0_endpts(bin1[0], U_list[0], 1, 2, mat_props, dpg, T))

    comp_coeff = []
    for j,k in zip([1,0,0],[2,2,1]):
        comp_coeff.append((4 * c[j] * c[k]) / ((1 + c[j] - c[k]) * (1 + c[k] - c[j])))
    rho_excess = (excess3 * comp_coeff[2] + excess2 * comp_coeff[1] + excess1 * comp_coeff[0])
    rho_vegard = sum([conc * 1/sig for conc, sig in zip(c, sigma_endpts)])
    rho_final = rho_vegard + rho_excess

    return (1 / rho_final)

def muggianu_model_redo(sigmae0_df, bin_list, U_list, mat_props, dpg, c : list, T = 300):
    c.append((1 + 1e-10) - sum(c))
#    bin3 = [((1 + c[0] - c[1]) / 2) , ((1 + c[1] - c[0]) / 2)]
#    bin2 = [((1 + c[0] - c[2]) / 2) , ((1 + c[2] - c[0]) / 2)]
#    bin1 = [((1 + c[1] - c[2]) / 2) , ((1 + c[2] - c[1]) / 2)]    
    U_coeff = []
    for j,k in zip([1,0,0],[2,2,1]):
        U_coeff.append((4 * c[j] * c[k]) / ((1 + c[j] - c[k]) * (1 + c[k] - c[j])))  
    U_tern = (U_list[2] * U_coeff[2] + U_list[1] * U_coeff[1] + U_list[0] * U_coeff[0])
    sigma_tern = sig.sigma_e0_model(c[0:2], U_tern, mat_props, dpg, T)
    return sigma_tern

def run_sigmae0_muggianu_dict(sigmae0_df, bin_list, U_list, mat_props, dpg, T = 300, n = 100):
    sig_tern = dict()
    first = 1e-10
    last = 9.9999999999999e-1
    j = 0
    for c in np.arange(first,1, (last - first) / n):
        k = 0
        for d in np.arange(first, 1 - c, (last - first) / n):
            sig_tern[(c*100,d*100)] = muggianu_model_sig_diff(sigmae0_df, bin_list, U_list,\
                     mat_props, dpg, [c,d], T)
            k = k+1
        j = j+1
    return sig_tern  


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
    #Drop band convergence
    sigmae0_df = sigmae0_df.drop([1,66])
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
    
    mat_props['sigma_endpts'] = [x_val, y_val, z_val]
    
    '''
    Split into binary data
    '''
    bin_list = split_binary_data(sigmae0_df)
    
    '''
    Fit separate U values
    '''
    U_list = fit_binary_U_list(bin_list, mat_props, d_key_list =  ['TiNiSn', 'HfNiSn', 'HfNiSn'], d_h_pairs = [[1,2],[0,2],[0,1]])
    #Why does the Hf-Ti binary have such a huge U parameter?
    
    
    '''
    Print binaries
    '''
    first = 1e-10
    last = 9.9999999999999e-1
    j = 0
    XY_bin = []
    XZ_bin = []
    YZ_bin = []
    for c in np.arange(first,1, (last - first) / 100):
        XY_bin.append(sig.binary_sigma_e0_endpts(c, U_list[2], 0, 1, mat_props, 'n', 300))
        XZ_bin.append(sig.binary_sigma_e0_endpts(c, U_list[1], 0, 2, mat_props, 'n', 300))
        YZ_bin.append(sig.binary_sigma_e0_endpts(c, U_list[0], 1, 2, mat_props, 'n', 300)) 
    c_list = np.arange(first,1, (last - first) / 100)
    plt.figure()
    plt.plot(c_list, XY_bin)
    plt.scatter(bin_list[2]['HfNiSn'], bin_list[2]['Sigma_E0'])
    plt.ylabel(r'$\sigma_{e0}$')
    plt.xlabel(r'x in Hf$_x$Ti$_{1-x}$NiSn')
    plt.figure()
    plt.plot(c_list, XZ_bin)
    plt.scatter(bin_list[1]['HfNiSn'], bin_list[1]['Sigma_E0'])
    plt.ylabel(r'$\sigma_{e0}$')
    plt.xlabel(r'x in Hf$_x$Zr$_{1-x}$NiSn')
    plt.figure()
    plt.scatter(bin_list[0]['TiNiSn'], bin_list[0]['Sigma_E0'])
    plt.ylabel(r'$\sigma_{e0}$')
    plt.xlabel(r'x in Ti$_x$Zr$_{1-x}$NiSn')
    plt.plot(c_list, YZ_bin)
    
    '''
    Muggianu Result
    '''
    #Look at the endpoints
    c_Ti = [0.01, 0.98]
    sig_mu_Hf = muggianu_model(sigmae0_df, bin_list, U_list, mat_props, 'n', c_Ti, T = 300)
    
    '''
    Plot ternary for muggianu result
    '''
    sig_data_dict = convert_data_df_to_dict(sigmae0_df, 'HfNiSn', 'TiNiSn')
    sig_mu_dict = run_sigmae0_muggianu_dict(sigmae0_df, bin_list, U_list, mat_props, 'n', T = 300, n = 100)
    fig, ax = plt.subplots()
    ax.axis("off")
    figure, tax = ternary.figure(ax=ax, scale = 100)
    
    tax.heatmap(sig_mu_dict, style = 'h', cmap=plt.cm.get_cmap('Spectral_r', 20),\
         cbarlabel=r'$\sigma_{e0}$', cb_kwargs = {'shrink' : 0.8, 'pad' : 0.01}, vmin = 70000, vmax = 110000,\
         scientific = False)
    
    Tt = [(sigmae0_df['HfNiSn'][i] * 100, sigmae0_df['TiNiSn'][i] * 100, sigmae0_df['ZrNiSn'][i] * 100) for i in  list(sigmae0_df.index)]
    At = sigmae0_df['Sigma_E0']
    tax.scatter(Tt, c = list(At), colormap=plt.cm.get_cmap('Spectral_r', 20),\
         cmap=plt.cm.get_cmap('Spectral_r', 20), vmin = 70000, vmax = 110000,\
         scientific = False, s = 30, edgecolors = 'k', zorder = 10, clip_on = False)
    
    tax.boundary(linewidth=2.0)
    
    tax.top_corner_label('TiNiSn')
    tax.left_corner_label('ZrNiSn', position = (0,0, 0))
    tax.right_corner_label('HfNiSn', position = (0.95,0, 0))
    
    tax.savefig('XNiSn_mm_sige0_ternary.pdf', bbox_inches = 'tight')
    
    '''
    trivial U mixing result
    '''
#    trivial_u = trivial_U_mixing(sigmae0_df, 'HfNiSn', 'ZrNiSn', mat_props, 'n', c = [0.1, 0.1, 0.8], T = 300)
    

    