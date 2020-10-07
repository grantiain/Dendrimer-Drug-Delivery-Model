#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from scipy.optimize import curve_fit

# Python3 code used to extract release half-lives from raw release data
# from in vitro release data


def release_dend(t, k_h, rmax):
    # Release model (assume negligible free drug at initial time)
    return(rmax * (1.0  - np.exp(-k_h * t)))


def fit_release(t_data, rel_data):
    # Returns fitted parameters release rate constant, max extent
    # and release half-life
    popt, pcov = curve_fit(release_dend, t_data, rel_data, 
                           bounds = ([0.0, 99.9], [1.0, 100.1]))    
    k_h = popt[0]
    rel_max = popt[1]
    t_h = np.log(2) / k_h
    return(k_h, rel_max, t_h)

# Data Mean Values

# SPL-8932
t_exp = [0.34, 1.10, 2.04, 2.80, 3.75, 4.51, 5.46, 6.22, 7.16, 7.92, 8.87, 
         9.63, 10.58, 11.33, 12.28, 13.04, 13.99, 14.75, 15.70, 16.45]
r_exp = [3.37, 11.00, 18.76, 25.73, 35.22, 43.64, 52.39, 58.59, 63.16, 67.43, 
         71.59, 72.38, 74.87, 76.53, 79.43, 80.07, 82.34, 82.38, 83.65, 84.70]

# SPL-8931
t_exp2 = [1.70, 11.68, 21.66, 31.64, 41.61, 67.81, 77.79, 87.77, 97.74, 107.72, 
          117.70, 127.68, 137.65, 147.63]
r_exp2 = [0.00, 3.55, 7.02, 9.61, 11.66, 18.44, 20.76, 23.77, 26.44, 28.96, 
          31.25, 34.17, 36.12, 38.96]

# SPL-8933                 
t_exp3 = [0.12, 0.31, 0.49, 0.68, 0.86, 1.18, 1.36, 1.55, 1.74]
r_exp3 = [1.45, 9.66, 17.21, 25.03, 32.41, 36.11, 42.14, 46.32, 49.89]

# Get optimised values for each dataset
# SPL-8932
k_h, rel_max, t_h = fit_release(t_exp, r_exp)
# SPL-8931
k_h2, rel_max2, t_h2 = fit_release(t_exp2, r_exp2)
# SPL-8933
k_h3, rel_max3, t_h3 = fit_release(t_exp3, r_exp3)

# Show Results  
              
print('\n' + 'Fitted Dendrimer Release Half-Lives' + '\n')
print('SPL-3931 t_half = ', int(10.0 * t_h2 + 0.5) / 10.0, ' hr')
print('SPL-3932 t_half = ', int(10.0 * t_h + 0.5) / 10.0, ' hr')
print('SPL-3933 t_half = ', int(10.0 * t_h3 + 0.5) / 10.0, ' hr')