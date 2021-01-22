#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 18:06:39 2020

@author: ramyagurunathan

Dispersion comparison: fixed issue and with subplot layout
"""
from __future__ import division, unicode_literals, print_function

import numpy as np
import ternary_tcond as tt

import disp_comp as dp
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

mpl.rcParams['font.sans-serif'] = 'Apple Symbols'
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['lines.linewidth'] = 3

if __name__ == '__main__':
    mpd = MPDataRetrieval('pfSJBa1OwitR5uNL')

    #this fetches the phonon bandstructure along a specific path in k-space
    bstruct = mpd.get_dataframe("mp-924128", ["phonon_bandstructure", "phonon_dos"])
    
    bs_symm = bstruct["phonon_bandstructure"][0]
    
    bs_dos = bstruct["phonon_dos"][0]
    
    dos_zeros = np.where(bs_dos.densities > 1e-02)
    
    branch = bs_symm.get_branch(0)
    #Generate the subplots
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharey = True)
    fig.tight_layout(pad=1)
    #Output the bandstructure data
    bspltr = PhononBSPlotter(bs_symm)
    bsplt_data = bspltr.bs_plot_data()
    segment = 0    
    for d,f in zip(bsplt_data['distances'], bsplt_data['frequency']):
        if not any(d > bsplt_data['distances'][segment][-1]):
            for freq in f:
                ax1.plot(d, freq, color = 'xkcd:blurple')
                ax1.tick_params(labelsize = 16, pad = 4)
                ax1.set_ylim(0, 8)
                ax1.set_xlim(0, bsplt_data['distances'][segment][-1])
                ax1.set_xticks([0, bsplt_data['distances'][segment][-1]])
                ax1.set_xticklabels([r'$\Gamma$', '$X$'])
                ax1.title.set_text('DFT')
    #Debye dispersion
    Vat = (4.324E-10)**3 * 2**(1/2) / 3 #atomic volume form Materials Project
    kDebye = (6 * pi**2 / Vat)**(1/3)
    debyeT = 307
    vs = tt.vs_from_debyeT(Vat, 307)
    k_array = np.linspace(0, kDebye)
    omegaD = (vs / (2 * pi)) * k_array * 1E-12
    ax2.tick_params(labelsize = 16, pad = 4)
    ax2.plot(k_array, omegaD, color = 'xkcd:blurple')
    ax2.set_xlim(0, kDebye)
    ax2.set_xticks([0, kDebye])
    ax2.set_xticklabels([r'$\Gamma$', r'k$_{max}$'])
    ax2.title.set_text('Debye')
    
    #BvK dispersion
    k_array = np.linspace(0, kDebye)
    omegaB = (2 / pi) * vs * kDebye * np.sin(pi * k_array / (2 * kDebye)) * 1E-12 * (1 / (2 * pi))
    ax3.tick_params(labelsize = 16, pad = 4)
    ax3.plot(k_array, omegaB, color = 'xkcd:blurple')
    ax3.set_xlim(0, kDebye)
    ax3.set_xticks([0, kDebye])
    ax3.set_xticklabels([r'$\Gamma$', r'k$_{max}$'])  
    ax3.title.set_text('Born von Karmann')
    
    #Debye + Einstein dispersion
    kD1 = (6 * pi**2 / (3 * Vat))**(1/3)
    kD2 = (6 * pi**2 / (2 * Vat))**(1/3)
    k_array1 = np.linspace(0, kD1)
    k_array2 = np.linspace(kD1, kD2)
    k_array3 = np.linspace(kD2, kDebye)
    omegaE1 = k_array1 * (vs / (2 * pi)) * 1E-12
    omegaE2 = np.ones(50) * kD2 * vs * 1E-12 * (1/(2 * pi))
    omegaE3 = np.ones(50) * kDebye * vs * 1E-12 * (1/(2 * pi))
    ax4.tick_params(labelsize = 16, pad = 4)
    ax4.plot(k_array1, omegaE1, color = 'xkcd:blurple')
    ax4.plot(k_array2, omegaE2, color = 'xkcd:blurple')
    ax4.plot(k_array3, omegaE3, color = 'xkcd:blurple')
    ax4.set_xlim(0, kDebye)
    ax4.set_xticks([0, kDebye])
    ax4.set_xticklabels([r'$\Gamma$', r'k$_{max}$'])
    ax4.title.set_text('Debye + Einstein')
    
    #Add common axis labels
    fig.text(-0.01, 0.35, 'Frequency (THz)', rotation = 'vertical', fontsize = 16)
    plt.savefig('disp_comp_new.pdf', bbox_inches = 'tight')
