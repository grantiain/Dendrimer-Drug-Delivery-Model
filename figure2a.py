#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from scipy.optimize import curve_fit

# Python3 code used to extract release half-lives from raw release data
# from in vitro release data

# Questions to iain.grant@astrazeneca.com

def release_dend(t, k_h, rmax):
    # Release model (assume negligible free drug at initial time)
    return(rmax * (1.0  - np.exp(-k_h * t)))


def fit_release(t_data, rel_data):
    # Returns fitted parameters release rate constant, max extent
    # and release half-life
    popt, pcov = curve_fit(release_dend, t_data, rel_data, 
                           bounds = ([0.0, 99.0], [1.0, 101]))    
    k_h = popt[0]
    rel_max = popt[1]
    t_h = np.log(2) / k_h
    return(k_h, rel_max, t_h)

# Experimental times and % released

# SPL-8932
t_exp = [0.23, 0.44, 0.99, 1.2, 1.94, 2.15, 2.7, 2.91, 3.64, 3.86, 4.4, 
              4.62, 5.35, 5.56, 6.11, 6.32, 7.06, 7.27, 7.82, 8.03, 8.76, 8.98, 
              9.52, 9.73, 10.47, 10.68, 11.23, 11.44, 12.18, 12.39, 12.93, 
              13.15, 13.88, 14.1, 14.64, 14.85, 15.59, 15.8, 16.35, 16.56]
r_exp = [1.63, 5.11, 10.36, 14.17, 19.17, 22.73, 26.62, 30.5, 37.96, 
                40.54, 48.37, 49.36, 57.57, 59.68, 66.07, 65.37, 67.49, 73.02, 
                72.7, 77.31, 77.2, 82.81, 78, 83.86, 80.54, 85.99, 82.09, 
                88.69, 86.09, 91.17, 86.56, 92.11, 88.57, 94.62, 88.4, 94.56, 
                90.06, 95.82, 90.47, 97.32]

# SPL-8931
t_exp2 = [1.59, 11.57, 21.55, 31.52, 41.5, 67.7, 77.68, 87.65, 97.63, 107.61, 
          117.59, 127.56, 137.54, 147.52, 1.7, 11.68, 21.66, 31.64, 41.61, 
          67.81, 77.79, 87.77, 97.74, 107.72, 117.7, 127.68, 137.65, 147.63]
r_exp2 = [0.11, 3.84, 7.65, 10.45, 12.68, 20.4, 22.74, 25.88, 28.4, 31.04, 
          33.67, 37.36, 39.07, 42.18, 0.02, 3.83, 7.49, 10.31, 12.62, 19.37, 
          22.54, 25.28, 28.62, 31.23, 33.52, 36.26, 38.81, 41.65]

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