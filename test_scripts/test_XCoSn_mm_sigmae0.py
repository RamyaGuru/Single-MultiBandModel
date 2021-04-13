#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 16:47:13 2021

@author: ramyagurunathan


Test script for XCoSn family
"""

import sys
sys.path.append('../electronic_models/')
import mixing_models_sigmae0 as mms
import sigmae0_alloy as sig
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import ternary
import csv


'''
constants
'''
e = 1.602e-19
m_e = 9.109e-31
hbar = 1.045e-34
kb = 1.38e-23

mat_props = {'formula': ['TaCoSn', 'VCoSn', 'NbCoSn'], 'BulkMod': [226.2E9, 201E9, 134E9],\
             'mpid': ['n/a', 'mp-1018119', 'mp-1094088'],\
             'AtmV' : [(5.948E-10)**3 / 12, 1.744861001976855e-29, (5.946E-10)**3 / 12]}



'''
Get the sigmae0 data
'''

data_path = '/Users/ramyagurunathan/Documents/PhDProjects/Argonne_TECCA/SigmaE0_values/'
datafile = data_path + 'XCoSn_gauss_data_black.csv'
sigmae0_df = sig.fetch_sigmae0_dataframe(datafile)