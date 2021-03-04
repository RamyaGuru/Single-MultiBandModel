#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 13:12:53 2021

@author: ramyagurunathan


Alloy Model for Sigma_e0


Plan:
    --> Fit the sigma_e0 values to end members
    --> Try to 
"""

import numpy as np
import sys
sys.path.append('/Users/ramyagurunathan/Documents/PhDProjects/Argonne_TECCA')
import pandas as pd
from mpcontribs.client import Client
api_key = 'pb59iL8clkMEjpJ6WV8xYWLDp084BJ8m'
client = Client(api_key)
import math
import re
import matplotlib.pyplot as plt
import matplotlib as mpl
import matminer as mm
from matminer.data_retrieval.retrieve_MP import MPDataRetrieval
from bandstructure_mm_mod import BandFeaturizer
import ternary
from scipy.optimize import curve_fit
import csv
from fdint import fdk
import qualityFactor as qf
from scipy.optimize import minimize

mpl.rcParams['font.size'] = 12

mpl.rcParams['figure.figsize'] = [6.4,4.8] 

project = 'carrier_transport'

mpd = MPDataRetrieval('pfSJBa1OwitR5uNL')

'''
constants
'''
e = 1.602e-19
m_e = 9.109e-31
hbar = 1.045e-34
kb = 1.38e-23

def conductivity_from_sige0(S, sigmae0):
    A = np.abs(S) / (kb / e)
    cond = sigmae0 / ((np.exp(A - 2) / (1 + np.exp(-5 * (A - 1)))) +\
                     ((3 * A / math.pi**2) / (1 + np.exp(5 * (A - 1)))))
    return cond

def fit_sigmae0(datafile, p0 = 1e6):
    '''
    Instead of using Jeff's formula.. fit using the Fermi itnegrals
    '''
    data_df = pd.read_csv(datafile)
    grps = data_df.groupby(by = [data_df.columns[1], data_df.columns[2]])
    split_comp = [grps.get_group(x) for x in grps.groups]
    sig_df_list = []
    for data in split_comp:
        sig = data['Sigma (S/m)']
        S = data['S (V/K)']
        sigmaE0, cov = curve_fit(lambda s, sigmaE0 : conductivity_from_sige0(s, sigmaE0), S, sig,\
                            p0 = p0, bounds = (0.1, np.inf), method = 'dogbox')
        
#    sigmae0 = sig * ((np.exp(A - 2) / (1 + np.exp(-5 * (A - 1)))) +\
#                     ((3 * A / math.pi**2) / (1 + np.exp(5 * (A - 1)))))
        data['Sigma_E0'] = np.ones(len(sig)) * sigmaE0
        sig_df_list.append(data)
    return sig_df_list


def fetch_sigmae0_dataframe(datafile):
    sigmae0_df = pd.read_csv(datafile)
    sigmae0_df = pd.DataFrame.dropna(sigmae0_df)
    return sigmae0_df

def fetch_endmember_values(sig_df):
    col_list = list(sig_df.columns)
    col_list.pop(0)
    col_list.remove('Sigma_E0')
#    print(sig_df.loc[sig_df[col_list[0]] > 9.99e-01]['Sigma_E0'])
    x_value = float(sig_df.loc[sig_df[col_list[0]] > 9.99e-01]['Sigma_E0'])
#    print(sig_df.loc[sig_df[col_list[1]] > 9.99e-01]['Sigma_E0'])
    y_value = float(sig_df.loc[sig_df[col_list[1]] > 9.99e-01]['Sigma_E0'])
    z_value = float(sig_df.loc[sig_df[col_list[2]] > 9.99e-01]['Sigma_E0'])
    return [x_value, y_value, z_value]

def get_mpcontribs_entry(table_id):
    df = client.get_table(table_id)
    return df

def get_eff_mass(df, doping_type):
    m_star = df['data'][0]['mₑᶜ'][doping_type]
    avg_m = re.sub(' mₑ', '', m_star['ε̄']['display'])
    avg_m = float(avg_m)
    return avg_m
    

def eff_mass_from_cond(df, doping_type):
    '''
    df : dataframe for the compound form MPContribs
    doping_type : 'n' or 'p'
    '''
    scinot_int = re.compile('([0-9]+\.?[0-9]+)(×10)([⁰¹²³⁴⁵⁶⁷⁸⁹]+)')  
    for t in df['tables'][0]:
        if t['name'] == 'σ(' + doping_type + ')':
            data_clean = np.zeros((t['total_data_rows'], len(t['columns'])))
            # Data table is rows: Temperature, doping conc: columns 
            #Need to clean strings more? better way?
            for r in range(len(t['data'])):
                for c in range(len(t['data'][0])):   
                    data_str = t['data'][r][c]
                    data_cln = re.sub(',', '', data_str)
                    str_list= re.findall(scinot_int, data_cln)
                    if str_list:
                        expnt = str_list[0][2]
                        #Fix this exponent replacement...
                        expnt = re.sub(u'\u00B9\u2079', '19', expnt)
                        expnt = re.sub(u'\u00B9\u2078', '18', expnt)
                        expnt = re.sub(u'\u00B9\u2077', '17', expnt)
                        expnt = re.sub(u'\u00B9\u2076', '16', expnt)
                        data_cln = float(str_list[0][0]) * 10 ** (float(expnt))
                    data_cln = float(data_cln)
                    data_clean[r][c] = (data_cln)**(-1)  
            n_array = np.array([float(n) * 1e6 for n in t['columns']])
            T_list = [float(T) for T in t['index']]
            n_rows = t['total_data_rows']
            cols = t['columns']
        else:
            continue
    n_matrix = np.tile(n_array, (n_rows, 1))
    eff_mass = data_clean * n_matrix * (e**2 / m_e)
    df['tables'][0].append({'name' : 'Conductivity effective mass', 'data': eff_mass,\
      'index' : T_list, 'columns' : cols, 'total_data_rows' : n_rows})
    return df

def print_eff_mass(df):
    plt.figure()
    eff_mass = df['tables'][0][-1]['data']
    for c in range(len(eff_mass[0])):
        plt.plot(df['tables'][0][-1]['index'] ,  eff_mass[:, c])
        plt.ylabel(r'm*$_c$')
        plt.xlabel('T (C)')
        plt.title(str(df['formula'][0]))
        
def mtminer_valley_degeneracy(mpid):
    '''
    Return effective mass as list: [n,p]
    '''
    try:
        bstruct = mpd.get_dataframe(mpid, ['bandstructure'])
        bfeat = BandFeaturizer()
        degen = bfeat.feature_nv(bstruct['bandstructure'][mpid])   
    except:
        degen = math.nan
    return degen  

def fit_sigma_e0_coefficient(sigma_e0_data, mat_props, dpg_type, T = 300): 
    coeff = []  
    mat_props['coeff_{}'.format(dpg_type)] = []    
    for n in range(len(mat_props['formula'])):
        s_e0_exp = sigma_e0_data[n]
        #e**2 in denominator is signifying deformation potential of 1eV
        A = s_e0_exp / ((2 * hbar * mat_props['BulkMod'][n] * mat_props['Nv'][n][dpg_type] * e**2)/\
        (3 * math.pi * mat_props['eff_mass_{}'.format(dpg_type)][n] * m_e * e**2))
        coeff.append(A)
    mat_props['coeff_{}'.format(dpg_type)] = coeff 
    return mat_props

def binary_sigma_e0_model(c, U, defect, host, mat_props, dpg, T = 300):
    '''
    Use Vegard's law for the coefficients, mat_props.
    Otherwise just use the sigma_e0 model?
    
    We'll call ZrNiSn the host
    And HfNiSn + TiNiSn the defects
    '''
    m_eff = mat_props['eff_mass_{}'.format(dpg)][defect] * c +\
    mat_props['eff_mass_{}'.format(dpg)][defect] * (1-c)
    
    Cl = mat_props['BulkMod'][defect] * c + mat_props['BulkMod'][host] * (1-c)
    
    coeff = mat_props['coeff_{}'.format(dpg)][defect] * c + mat_props['coeff_{}'.format(dpg)][host] * (1-c)
    
    AtmV = mat_props['AtmV'][defect] * c + mat_props['AtmV'][host] * (1-c)
    
    Nv = mat_props['Nv'][2][dpg]
    #Need to figure out what's going on with A factor
    sigmae0 = (2 * hbar * Cl * Nv * e**2 * coeff) / (3 * math.pi * m_eff * m_e * e**2 *\
               (1 + U * (3 * math.pi**2 * c * (1-c) * Cl * AtmV)/\
                (8 * kb * T)))
    return sigmae0   

def binary_sigma_e0_endpts(c, U, defect, host, mat_props, dpg, T = 300):
    sigmae0_vegard = mat_props['sigma_endpts'][defect] * c +\
    mat_props['sigma_endpts'][host] * (1-c)
    
    Cl = mat_props['BulkMod'][defect] * c + mat_props['BulkMod'][host] * (1-c)
    
    AtmV = mat_props['AtmV'][defect] * c + mat_props['AtmV'][host] * (1-c)
    
    A = (3 * math.pi**2 * c * (1-c) * Cl * AtmV) / (8 * kb * T)
    sigma_e0 = sigmae0_vegard * (1 / (1 + A * U))
    return sigma_e0
    
def sigma_e0_model(c: list, U, mat_props, dpg, T = 300):
    '''
    Use Vegard's law for the coefficients, mat_props.
    Otherwise just use the sigma_e0 model?
    
    We'll call ZrNiSn the host
    And HfNiSn + TiNiSn the defects
    '''
    defect_conc = sum(c)
    host_conc = (1 + 1e-10) - defect_conc
    m_eff = mat_props['eff_mass_{}'.format(dpg)][0] * c[0] +\
    mat_props['eff_mass_{}'.format(dpg)][1] * c[1] +\
    mat_props['eff_mass_{}'.format(dpg)][2] * host_conc
    
    Cl = mat_props['BulkMod'][0] * c[0] + mat_props['BulkMod'][1] * c[1] +\
    mat_props['BulkMod'][2] * host_conc
    
    coeff = mat_props['coeff_{}'.format(dpg)][0] * c[0] + mat_props['coeff_{}'.format(dpg)][1] * c[1] +\
    mat_props['coeff_{}'.format(dpg)][2] * host_conc
    
    AtmV = mat_props['AtmV'][0] * c[0] + mat_props['AtmV'][1] * c[1] +\
    mat_props['AtmV'][2] * host_conc
    
    Nv = mat_props['Nv'][2][dpg]
    #Need to figure out what's going on with A factor
    sigmae0 = (2 * hbar * Cl * Nv * e**2 * coeff) / (3 * math.pi * m_eff * m_e * e**2 *\
               (1 + U * (3 * math.pi**2 * host_conc * defect_conc * Cl * AtmV)/\
                (8 * kb * T)))
    return sigmae0

def fit_binary_U(data_df, D_key, defect, host, mat_props, dpg, T = 300, p0 = 0):
    Tt = np.array(data_df[D_key])
    At = list(data_df['Sigma_E0'])
    U, cov = curve_fit(lambda C, U: binary_sigma_e0_model(C, U, defect, host, mat_props, dpg, T), Tt, At,\
                       p0 = p0, bounds = (0,np.inf), method = 'dogbox')
    return U[0], cov 

def fit_binary_U_endpts(data_df, D_key, defect, host, mat_props, dpg, T = 300, p0 = 0):
    Tt = np.array(data_df[D_key])
    At = list(data_df['Sigma_E0'])
    U, cov = curve_fit(lambda C, U: binary_sigma_e0_model(C, U, defect, host, mat_props, dpg, T), Tt, At,\
                       p0 = p0, bounds = (0,np.inf), method = 'dogbox')
    return U[0], cov

def fit_U(data_df, X_key, Y_key, mat_props, dpg, T = 300, p0 = 0):
#    Tt_X = np.array([data_df[X_key][i]  for i in  list(data_df.index)])
#    Tt_Y = np.array([data_df[Y_key][i]  for i in  list(data_df.index)])
    Tt_X = np.array(data_df[X_key])
    Tt_Y = np.array(data_df[Y_key])
    Tt = (Tt_X, Tt_Y)
    At = list(data_df['Sigma_E0'])
    U, cov = curve_fit(lambda C, U: sigma_e0_model(C, U, mat_props, dpg, T), Tt, At,\
                       p0 = p0, bounds = (0,np.inf), method = 'dogbox')
    return U[0], cov
    

def run_sigmae0_tern(mat_props, dpg, T = 300, n = 100, U = 1):
    length = n * n
    sig_tern = np.zeros(length)
    k = 0
    t = []
    u = []
    v = []
    for c in np.linspace(1e-10,9.9999999e-1,100):
        for d in np.linspace(1e-10, 9.9999999e-1 - c, 100):
            t.append(c)
            u.append(d)
            v.append((1 + 1e-10)-c-d)
            sig_tern[k] = sigma_e0_model([c,d], U, mat_props, dpg, T)
            k = k+1
    return sig_tern, [t,u,v]

def run_sigmae0_tern_data_dict(mat_props, dpg, T = 300, n = 100, U = 1):
    sig_tern = dict()
    first = 1e-10
    last = 9.9999999999999e-1
    j = 0
    for c in np.arange(first,1, (last - first) / n):
        k = 0
        for d in np.arange(first, 1 - c, (last - first) / n):
            sig_tern[(c*100,d*100)] = sigma_e0_model([c,d], U, mat_props, dpg, T)
            k = k+1
        j = j+1
    return sig_tern

#Add naive model that is just combination of Vegard's law for sigma_e0 values and x(1-x) mixing rule
def naive_sigmae0_data_dict(sig_endpts, n = 100):
    sig_tern = dict()
    first = 1e-10
    last = 9.9999999999999e-1
    j = 0
    for c in np.arange(first,1, (last - first) / n):
        k = 0
        for d in np.arange(first, 1 - c, (last - first) / n):
            sig_tern[(c*100,d*100)] = c * sig_endpts[0] + d * sig_endpts[1] +\
            (1-c-d) * sig_endpts[2]# / (1 + (c + d) * (1 - (c + d)))
            k = k+1
        j = j+1
    return sig_tern    

#def fit_coeff_sigmae0()

if __name__ == '__main__':
    data_path = '/Users/ramyagurunathan/Documents/PhDProjects/Argonne_TECCA/SigmaE0_values/'
    datafile = data_path + 'XNiSn_gauss_data_purple.csv'
    sigmae0_df = fetch_sigmae0_dataframe(datafile)
    #Drop the band convergence sample
    sigmae0_df = sigmae0_df.drop(66)
    x_val, y_val, z_val = fetch_endmember_values(sigmae0_df)
    '''
    Initalize Endmember Properties
    Order of compounds: HfNiSn, TiNiSn, ZrNiSn
    
    Average Effective Mass values: Obtained from the MPContribs database
    Dividing the carr. conc., eff. mass, 
    '''
    mat_props = {'formula': ['HfNiSn', 'TiNiSn', 'ZrNiSn'], 'BulkMod': [105E9, 122E9, 152E9],\
                 'mpid': ['mp-20523', 'mp-924130', 'mp-30806'],\
                 'AtmV' : [1.91E-29, 1.755E-29, 1.93E-29]}
    mat_props['sigma_endpts'] = [x_val, y_val, z_val]
    property_path = '/Users/ramyagurunathan/Documents/PhDProjects/Argonne_TECCA/'
    
#    mp_entry = get_mpcontribs_entry('5f8a62d6f02214a1d07e2d2a')
    
    X_df = pd.read_json(property_path + mat_props['formula'][0])
    Y_df = pd.read_json(property_path + mat_props['formula'][1])
    Z_df = pd.read_json(property_path + mat_props['formula'][2])    
    
#    ex_avg_m = get_eff_mass(ex_df, 'n')
    X_m_n = eff_mass_from_cond(X_df, 'n')
    Y_m_n = eff_mass_from_cond(Y_df, 'n')
    Z_m_n = eff_mass_from_cond(Z_df, 'n')

#    ex_avg_m = get_eff_mass(ex_df, 'n')
    X_m_p = eff_mass_from_cond(X_df, 'p')
    Y_m_p = eff_mass_from_cond(Y_df, 'p')
    Z_m_p = eff_mass_from_cond(Z_df, 'p')  

#Add the effectivs mass to the dictionary
    #Source: PCCP 2016, 18 14107
    mat_props['eff_mass_n'] = [0.36, 0.56, 0.38]
    mat_props['eff_mass_p'] = [0.22, 0.39, 0.25]
    
    print_eff_mass(X_m_n)
    print_eff_mass(Y_m_n)
    
    print_eff_mass(Y_m_p)

#Add the valley degeneracy
    X_Nv = {'p' : 1, 'n' : 3}
    Y_Nv = {'p' : 1, 'n' : 3} #mtminer_valley_degeneracy('mp-924130')
    Z_Nv = {'p' : 1, 'n' : 3}
    

    mat_props['Nv'] = [X_Nv, Y_Nv, Z_Nv]

#Add the coefficients
    mat_props = fit_sigma_e0_coefficient([x_val, y_val, z_val], mat_props, 'n', T = 300)

#Sigmae0
    sigmae0 = sigma_e0_model([0.1, 0.2], 1, mat_props, 'n')

#%%

    '''
    Print ternary
    '''
    u, cov = fit_U(sigmae0_df, 'HfNiSn', 'TiNiSn', mat_props, 'n', T = 300, p0 = 0)
    sig_dict = run_sigmae0_tern_data_dict(mat_props, 'n', T = 300, n = 100, U = u)
    fig, ax = plt.subplots()
    ax.axis("off")
    figure, tax = ternary.figure(ax=ax, scale = 100)
    
    tax.heatmap(sig_dict, style = 'h', cmap=plt.cm.get_cmap('Spectral_r', 20),\
         cbarlabel=r'$\sigma_{e0}$',\
         scientific = False)
    
    Tt = [(sigmae0_df['HfNiSn'][i] * 100, sigmae0_df['TiNiSn'][i] * 100, sigmae0_df['ZrNiSn'][i] * 100) for i in  list(sigmae0_df.index)]
    At = sigmae0_df['Sigma_E0']
    tax.scatter(Tt, c = list(At), colormap=plt.cm.get_cmap('Spectral_r', 20),\
         cmap=plt.cm.get_cmap('Spectral_r', 20), vmin = 80000, vmax = 110000,\
         scientific = False, s = 30, edgecolors = 'k', zorder = 10, clip_on = False)
    
    tax.boundary(linewidth=2.0)
    
    tax.top_corner_label('TiNiSn')
    tax.left_corner_label('ZrNiSn', position = (0,0.04, 0))
    tax.right_corner_label('HfNiSn', position = (0.95,0.04, 0))
    
#Naive Ternary Plot
    naive_dict = naive_sigmae0_data_dict([x_val, y_val, z_val])
    fig1, ax1 = plt.subplots()
    ax1.axis("off")
    figure1, tax1 = ternary.figure(ax=ax1, scale = 100)    
    tax1.heatmap(naive_dict, style = 'h', cmap=plt.cm.get_cmap('Spectral_r', 20),\
         cbarlabel=r'$\sigma_{e0}$',\
         scientific = False) 
    tax1.scatter(Tt, c = list(At), colormap=plt.cm.get_cmap('Spectral_r', 20),\
         cmap=plt.cm.get_cmap('Spectral_r', 20), vmin = 80000, vmax = 110000,\
         scientific = False, s = 30, edgecolors = 'k', zorder = 10, clip_on = False)

    tax1.boundary(linewidth=2.0)
    
    tax1.top_corner_label('TiNiSn')
    tax1.left_corner_label('ZrNiSn', position = (0,0.04, 0))
    tax1.right_corner_label('HfNiSn', position = (0.95,0.04, 0))
#    tax.savefig('klemens_model_xnisn_lowk.pdf', bbox_inches = 'tight')   
    
    with open('XNiSn_sige0_data.csv', 'w') as csvfile:
        field_names = ['% (Hf)', '% (Ti)', '% (Zr)', 'Sigma_e0']
        writer = csv.DictWriter(csvfile, fieldnames  = field_names)
        writer.writeheader()
        for k,v in sig_dict.items():
            writer.writerow({'% (Hf)': k[0], '% (Ti)' : k[1], '% (Zr)' : 100 - (k[0] + k[1]), 'Sigma_e0' : v})
            
    '''
    Plot Quality Factor
    '''
    fig, ax = plt.subplots()
    ax.axis("off")
    figure, tax = ternary.figure(ax=ax, scale = 100)
    
    kL_df = pd.read_csv('../XNiSn_data.csv')
    
    sig_values = list(sig_dict.values())
    
    T = 300
    
    B = (kb / e)**2 * 300 * (np.array(sig_values) / np.array(kL_df['kappa_lattice']))
    
    
    B_dict = sig_dict
    
    j = 0
    for k in B_dict.keys():
        B_dict[k] = B[j]
        j = j+1
        
    tax.heatmap(B_dict, style = 'h', cmap=plt.cm.get_cmap('Spectral_r', 20),\
         cbarlabel=r'B',\
         scientific = False)
    
    tax.boundary(linewidth=2.0)
    
    tax.top_corner_label('TiNiSn')
    tax.left_corner_label('ZrNiSn', position = (0,0.04, 0))
    tax.right_corner_label('HfNiSn', position = (0.95,0.04, 0))

    with open('XNiSn_Bfactor_data.csv', 'w') as csvfile:
        field_names = ['% (Hf)', '% (Ti)', '% (Zr)', 'BFactor']
        writer = csv.DictWriter(csvfile, fieldnames  = field_names)
        writer.writeheader()
        for k,v in B_dict.items():
            writer.writerow({'% (Hf)': k[0], '% (Ti)' : k[1], '% (Zr)' : 100 - (k[0] + k[1]), 'BFactor' : v})
    
    mpl.rcParams['figure.figsize'] = [5,5]
    fig, ax = plt.subplots()
    ax.axis("off")
    figure, tax = ternary.figure(ax=ax)
    
    tax.boundary(linewidth=2.0)
    
    tax.top_corner_label('TiNiSn')
    tax.left_corner_label('ZrNiSn', position = (0,0.04, 0))
    tax.right_corner_label('HfNiSn', position = (0.95,0.04, 0))
    #%%
            
    '''
    Sigma_e0 from Seebeck and conductivity
    '''
    
    new_datafile = data_path + 'XNiSn_S_and_sigma.csv'
    sige0 = fit_sigmae0(new_datafile)
    
    '''
    Jonker Plots
    '''
    for data in sige0:
        plt.figure()
        plt.scatter(data['Sigma (S/m)'], data['S (V/K)'] * 1e6)
        cond = conductivity_from_sige0(data['S (V/K)'], data['Sigma_E0'])
        plt.scatter(cond, data['S (V/K)'] * 1e6)
    
    
    