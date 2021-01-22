#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 11:31:31 2018

@author: guru

script to play aorund with pymatgen phonon tools and matminer
"""

from __future__ import division, unicode_literals, print_function

import numpy as np
import ternary_tcond as tt

#from pymatgen.matproj.rest import MPRester
#mpr = MPRester('pfSJBa1OwitR5uNL')

from matminer.featurizers.base import BaseFeaturizer
from pymatgen.phonon.bandstructure import PhononBandStructure, PhononBandStructureSymmLine
from pymatgen.phonon.dos import PhononDos, CompletePhononDos
from pymatgen.phonon.plotter import PhononDosPlotter, PhononBSPlotter
from pymatgen import MPRester
from pymatgen.ext.matproj import MPRestError
from pymatgen.util.plotting import pretty_plot
import matplotlib as mpl
import matplotlib.pyplot as plt
from math import pi, sin

from matminer.data_retrieval.retrieve_MP import MPDataRetrieval 
#import pandas as pd

mpl.rcParams['xtick.major.pad'] = 0

#Plotting DOS
def plot_dos(bs_dos):
    dosplt = PhononDosPlotter()
    dosplt.add_dos("Total", bs_dos)
    dosplt.get_plot().show()


#Bandstructure plotting
def plot_bands(bs_symm, segment = 0, figname = None):
    bspltr = PhononBSPlotter(bs_symm)
    bsplt_data = bspltr.bs_plot_data()
    bsplt = bspltr.get_plot()
    bsplt.xlim(0, bsplt_data['distances'][segment][-1])
    bsplt.ylim(top = 9)
    bsplt.xlabel('Wavevector (k)', fontsize = 18)
    bsplt.ylabel('Frequencies (THz)', fontsize = 18)
    bsplt.title('DFT', fontsize = 16)
    fig = bsplt.gcf()
    fig.set_size_inches(2,7) 
    bsplt.axes().tick_params(labelsize = 16)
    if figname:
        bsplt.savefig(figname + '.pdf', bbox_inches = 'tight')

#Plot the partial density of states from each of the sites
def plot_partial_dos(bs_symm, bs_dos, dosplt):
    i=0
    for site in bs_symm.structure.sites:
        site_dos = bs_dos.get_site_dos(site)
        sitename = 'site' + str(i)
        element_dos = bs_dos.get_element_dos()
        dosplt.add_dos(sitename, site_dos)
        i=i+1
    dosplt.get_plot().show()


if __name__ =='__main__':
    mpd = MPDataRetrieval('pfSJBa1OwitR5uNL')

    #this fetches the phonon bandstructure along a specific path in k-space
    bstruct = mpd.get_dataframe("mp-924128", ["phonon_bandstructure", "phonon_dos"])
    
    bs_symm = bstruct["phonon_bandstructure"][0]
    
    bs_dos = bstruct["phonon_dos"][0]
    
    dos_zeros = np.where(bs_dos.densities > 1e-02)
    
    branch = bs_symm.get_branch(0)
    
    #Plot the bandstructure
    plot_bands(bs_symm, figname = 'HfNiSn_bnd_diagram')
    
    #Output the bandstructure data
    bspltr = PhononBSPlotter(bs_symm)
    bsplt_data = bspltr.bs_plot_data()
    
    '''
    1D: Debye Model  + Einstein Plot-- might have to do some mass ratio thing?
    '''
    atmV = (6.12E-10)**3 / 12
    atmV_MP = (4.324E-10)**3 * 2**(1/2) / 3
    debyeT = 307
    vs = tt.vs_from_debyeT(atmV, 307)
    omegaA = np.zeros(50)
    omegaO1 = np.zeros(50)
    omegaO2 = np.zeros(50)
    kmax = 1.0383633591718584
    k_debye = (6 * pi**2 / atmV_MP)**(1/3)
    k_array = np.linspace(0, kmax)
    omegaD1 = vs * (6 * pi**2 / (2 * atmV))**(1/3)
    omegaD2 = vs * (6 * pi**2 / (atmV))**(1/3)
    i = 0
    for k in k_array:
        omegaA[i] = vs * k * 1E-3
        omegaO1[i] = vs * kmax * 2 * 1E-3
        omegaO2[i] = vs * kmax * 3 * 1E-3
        i = i+1
    pplt = pretty_plot() 
    fig = pplt.gcf()
    fig.set_size_inches(2,6) 
    pplt.axes().tick_params(labelsize = 16, pad = 4)
    pplt.axes().set_xlim(0,1)
    pplt.axes().set_ylim(-0.4,9)
    plt.axes().set_xticks([0,1])
    pplt.axes().set_xticklabels([r'$\Gamma$', 'X'])
    pplt.plot(k_array, omegaA, linewidth = 1.5, color = 'blue')
    pplt.plot(k_array, omegaO1, linewidth = 1.5, color = 'blue')
    pplt.plot(k_array, omegaO2, linewidth = 1.5, color = 'blue')
    fig.tight_layout(pad = 0)
            # Main X and Y Labels
#    pplt.xlabel(r'$\mathrm{Wave\ Vector\ (k)}$', fontsize=18)
    ylabel = r'$\mathrm{{Frequencies\ (THz)}}$'
#    pplt.ylabel(ylabel, fontsize=18)
    pplt.title('Debye+Einstein', fontsize = 16)
    pplt.savefig('1D_debye_ein_hfnisn.pdf', bbox_inches = 'tight')
    
    '''
    1D: Debye Model Plot
    '''
    omegaA0 = np.zeros(50)
    omegaA1 = np.zeros(50)
    omegaA2 = np.zeros(50)
    j = 0
    for k in k_array:
        omegaA0[j] = vs * k * 1E-3
        omegaA1[j] = vs * kmax * 1E-3 + (kmax - k) * vs * 1E-3
        omegaA2[j] = vs * kmax * 2 * 1E-3 + vs * k * 1E-3
        j = j+1
    pplt = pretty_plot() 
    fig = pplt.gcf()
    fig.set_size_inches(2,6) 
    pplt.axes().tick_params(labelsize = 16, pad = 4)
    pplt.axes().set_xlim(0,1)
    pplt.axes().set_ylim(-0.4,9)
    plt.axes().set_xticks([0,1])
    pplt.axes().set_xticklabels([r'$\Gamma$', 'X'])
    pplt.plot(k_array, omegaA0, linewidth = 1.5, color = 'blue')
    pplt.plot(k_array, omegaA1, linewidth = 1.5, color = 'blue')
    pplt.plot(k_array, omegaA2, linewidth = 1.5, color = 'blue')
    fig.tight_layout(pad = 0)
            # Main X and Y Labels
#    pplt.xlabel(r'$\mathrm{Wave\ Vector\ (k)}$', fontsize=18)
    ylabel = r'$\mathrm{{Frequencies\ (THz)}}$'
#    pplt.ylabel(ylabel, fontsize=18)  
    pplt.title('Debye', fontsize = 16)
    pplt.savefig('1D_debye_hfnisn.pdf', bbox_inches = 'tight')
    
    '''
    1D: Born von Karmann Plot
    '''
    omegaA0 = np.zeros(50)
    omegaA1 = np.zeros(50)
    omegaA2 = np.zeros(50)
    j = 0
    for k in k_array:
        omegaA0[j] =(2 / pi) * vs * 3 * kmax * sin(pi * k / (2 * 3 * kmax)) * 1E-3
        omegaA1[j] = (2 / pi) * vs * 3 * kmax * sin(pi * (kmax + (kmax - k)) / (2 * 3 * kmax)) * 1E-3
        omegaA2[j] = (2 / pi) * vs * 3 * kmax * sin(pi * (2 * kmax + k) / (2 * 3 * kmax)) * 1E-3
        j = j+1
    pplt = pretty_plot() 
    fig = pplt.gcf()
    fig.set_size_inches(2,6) 
    pplt.axes().tick_params(labelsize = 16, pad = 4)
    pplt.axes().set_xlim(0,1)
    pplt.axes().set_ylim(-0.4,9)
    plt.axes().set_xticks([0,1])
    pplt.axes().set_xticklabels([r'$\Gamma$', 'X'])
    pplt.plot(k_array, omegaA0, linewidth = 1.5, color = 'blue')
    pplt.plot(k_array, omegaA1, linewidth = 1.5, color = 'blue')
    pplt.plot(k_array, omegaA2, linewidth = 1.5, color = 'blue')
    fig.tight_layout(pad = 0)
            # Main X and Y Labels
#    pplt.xlabel(r'$\mathrm{Wave\ Vector\ (k)}$', fontsize=18)
    ylabel = r'$\mathrm{{Frequencies\ (THz)}}$'
#    pplt.ylabel(ylabel, fontsize=18)
    pplt.title('Born von Karmann', fontsize = 16)
    pplt.savefig('1D_bvk_hfnisn.pdf', bbox_inches = 'tight')
      

    

