#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# =============================================================================
# Python script to optimize the values of k_ext and k_tu to collectively fit
# the first set of dendrimer data with tumour PK of released active AZD4320
#
# Fits data to Equation 3 in the manuscript
#
# iain.grant@astrazeneca.com
# =============================================================================

import numpy as np
from scipy.optimize import fmin


def conc_tu(t_tu, t_hy, t_res, dose, V_tu, k_ext, t):

    # Returns the concentration of released API in tumour tissue
    # with a fixed release rate in both tumour and plasma
    k_res = np.log(2.0) / t_res
    k_tu = np.log(2.0) / t_tu
    k_h = np.log(2.0) / t_hy  
    nm = (k_res * np.exp(-k_tu * t) + (k_h - k_tu) * np.exp(-(k_h + k_res) * t)
    + (k_tu - k_h - k_res) * np.exp(-k_h * t))
    dnm = (k_h - k_tu + k_res) * (k_h - k_tu)

    return((k_h * nm / dnm) * (dose * k_ext / (k_res * V_tu)))



def objective_fn(x, t_res, t_hyd, dose, V_tu, t_exp, c_tu_rel, sd_c_tu_rel):

    # Objective function to fit k_res and k_ext to first dataset
    t_tu = x[0]
    k_ext = x[1]
    tssq = 0.0

    for j in range(0,3):
        for i in range(0,len(c_tu_rel[j])):
            y_hat = conc_tu(t_tu, t_hyd[j], t_res, dose, V_tu, 
                            k_ext, t_exp[j][i])
            tssq += ((c_tu_rel[j][i]-y_hat) ** 2 / y_hat) \
            * (1.0 / sd_c_tu_rel[j][i])

    return(tssq)
    

# =============================================================================
# Experimental Data in order SPL-8933, SPL-8932 and SPL-8931
# =============================================================================

# Experimental time points for tumour PK
t_exp = [[0.5, 1.0, 2.0, 5.0, 25.0],
         [1.0, 2.5, 5.0, 24.0],
         [5.0, 24.0, 31.0, 48.0]]
         
# Concentrations of Releasd AZD4320 in tumour tissue
c_tu_rel = [[0.03775, 0.0414, 0.167, 0.3355, 0.308],
            [0.01193, 0.066, 0.123, 0.5835],
            [0.010165, 0.0233, 0.0739, 0.11645]]

# Standard deviation on 3 replication concentration measurements            
sd_c_tu_rel = [[0.004596, 0.0280, 0.0396, 0.0572, 0.09],
               [0.00603869, 0.03040559, 0.01414214, 0.23829499],
               [0.00528209, 0.00282843, 0.01103087, 0.02906209]]

# Experimental hydrolysis half-lives of the first set of dendrimers
t_hyd = [1.7, 4.4,  201.0]


# =============================================================================
# Fixed Inputs
# =============================================================================

# Half-life of dendrimer in systemic circulation in hours
t_res = 3.301

# Total dose in mg
dose = 0.25

# Tumour volume in cubic litres
V_tu = 3.0e-4

# =============================================================================
# Minimise 'objective function' to find optimised values for t_res and k_ext
# =============================================================================

init_val = [12.0, 1.0e-4]
res = fmin(objective_fn, init_val, args=(t_res, t_hyd, dose, V_tu, t_exp, 
                                         c_tu_rel, sd_c_tu_rel), disp=False)
t_tu = res[0]
k_ext = res[1]
k_tu = np.log(2.0)/t_tu

print('\n' + 'Fitted Parameters' + '\n')
print ('k_tu = ', int(1000.0 * k_tu + 0.5) / 1000.0, ' per hr')
print ('k_ext = ', int(1.0e6 * k_ext + 0.5) / 1.0e6, ' per hr')