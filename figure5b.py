# -*- coding: utf-8 -*-
"""

@author: GRANT_I
AstraZeneca, Macclesfield, UK

Questions to iain.grant@astrazeneca.com

"""

# Python 3 script to fit equation 1 in manuscript to experimental data
# Total in Plasma over time

import numpy as np
from scipy.optimize import curve_fit

def plasma_total(t, k_res, cmax):
    # Total concentration in plasma (assume negligible contribution from the
    # released drug)
    # Equation 1 in manuscript
    cpl_total = cmax * np.exp(-(k_hy + k_res) * t)
    return cpl_total


def fit_release(t_data, rel_data):
    # Returns fitted parameters release rate constant, max extent
    # and release half-life
    popt, pcov = curve_fit(plasma_total, t_8932, c_pl8932)    
    k_res = popt[0]
    cmax = popt[1]

    return(k_res, cmax)


# From figure 2a we have a release half-life for SPL-8932 of 4.4 hours
    
# convert to hydrolysis rate constant
k_hy = np.log(2) / 4.4

# Experimental Data for SPL-8932 (time in hours and total concentration in
# plasma in micrograms per ml)
t_8932 = [1, 3, 6, 9, 16, 28]
c_pl8932 = [116.333, 59.2, 17.467, 2.977, 0.189, 0.03]

# Fit model to experimental data for SPL-8932 total in plasma
k_res, c_max = fit_release(t_8932, c_pl8932)
dose = 250.0 # micrograms (10 mg/kg to a 25 g mouse)

# Effecitve volume of distribution for dendrimer 'initially'
Vc = dose / (c_max * 0.025)

# Show fitted results
print('\n' + 'Fitted Value for k_res RES Uptake Constant' + '\n')
print ('k_res = ', int(100.0 * k_res + 0.5) / 100.0, ' per hr')
print ('Vc = ', int(10.0 * Vc + 0.5) / 10.0, ' ml per kg')