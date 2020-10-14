#!/usr/bin/env python
# coding: utf-8

# =============================================================================
# Python script to optimize the PK parameters of the Naked API
#
# iain.grant@astrazeneca.com
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import minimize


def c_total(t, k_res, k_h, dose, V_c):
    #c_max = dose / V_c
    c_max = 170.0
    c_tot = c_max * np.exp(-(k_res+k_h) * t)
    return c_tot
    


# Definition of Model
def dxdt(X, t, CL, Vc, k12, k21, k13, k31):
    
    # Iain Grant, AstraZeneca, Macclesfield
    # iain.grant@astrazeneca.com
    dx0 = -((CL / Vc) + k12 + k13) * X[0] + k21 * X[1] + k31 * X[2]
    dx1 = k12 * X[0] - k21 * X[1]
    dx2 = k13 * X[0] - k31 * X[2]
    return dx0, dx1, dx2



# Objective function returning the sum of the residuals squared
def objective_function(parms, t_data, c_data, dose):
    
    # Iain Grant, AstraZeneca, Macclesfield
    # iain.grant@astrazeneca.com
    # Parameters to be fitted
    CL = parms[0]
    Vc = parms[1]
    k12 = parms[2]
    k21 = parms[3]
    k13 = parms[4]
    k31 = parms[5]
    
    # initial conditions
    X_init = [dose, 0.0, 0.0]

    soln = odeint(dxdt, X_init, t_data, args=(CL, Vc, k12, k21, k13, k31))
    
    # Convert to concentration
    
    c_model = soln[:, 0] / Vc
    # calculate sum of the residuals squared (with 1/yhat ** 2 weighting)
    res_sq = 0.0
    for i in range(0,len(c_data)):
        res_sq += ((c_data[i] - c_model[i]) ** 2 / c_data[i] ** 2)
    return res_sq



# Plot model and experimental data
def plot_data(soln, t_soln, t, c_plasma, Vc, gcol):
    
    # Iain Grant, AstraZeneca, Macclesfield
    # iain.grant@astrazeneca.com
    plt.plot(t_soln, soln[:,0] / Vc, c = gcol)
    plt.plot(t, c_plasma, ls = 'none', marker = 'o', c = gcol, 
             markeredgecolor = 'k', markersize = '8')
    plt.xlabel('Time [hr]', fontsize = 16)
    plt.ylabel('Concentration [$\mu$g/mL]', fontsize = 16)
    plt.grid(b=True)
    plt.show()



def fit_model(init_vals, dose, t_data, c_data, bnds):
    
    # Iain Grant, AstraZeneca, Macclesfield
    # iain.grant@astrazeneca.com
    res = minimize(objective_function, init_vals, args=(t_data, c_data, dose), 
                   bounds = bnds)    
    return res.x


# Mouse PK at 5 mg/kg
t_data = [0.0, 0.05, 0.25, 0.5, 1, 2.0, 4.0, 6.0, 8.0]
c_data = [30.0, 15.3, 2.753, 0.977, 0.436, 0.0753, 0.0863, 0.0601, 0.0413]

t_data_all = [0.05, 0.25, 0.5, 1, 2, 4, 6, 8, 0.05, 0.25, 0.5, 1, 2, 4, 6, 
              8, 0.05, 0.25, 0.5, 1, 2, 4, 6, 8]
c_data_all = [13.9, 2.89, 1.08, 0.579, 0.0797, 0.131, 0.0921, 0.0529, 16, 
              2.74, 1.03, 0.402, 0.0678, 0.0544, 0.0395, 0.0534, 16, 2.63, 
              0.822, 0.328, 0.0784, 0.0736, 0.0488, 0.0176]


parms_init = [0.03, 0.004, 2.0, 0.02, 1.0, 1.0]
pmin = [0.01, 0.001, 0.001, 0.001, 0.001, 0.001]
pmax = [0.1, 0.02, 5.0, 5.0, 5.0, 5.0]
bnds = np.c_[pmin, pmax]


ts = np.linspace(0, 12, 300)

# 5 mg/kg in 25 g mouse
dose = 0.125

# Fit model find optimised values for PK parameters
CL, Vc, k12, k21, k13, k31 = fit_model(parms_init, dose, t_data, c_data, bnds)

# Solve model with optimised parameters
X_init = [dose, 0.0, 0.0]
soln = odeint(dxdt, X_init, ts, args=(CL, Vc, k12, k21, k13, k31))


# Plot the experimental data and the fitted model
plt.plot(ts, soln[:, 0] / Vc, color = 'k')
plt.plot(t_data_all, c_data_all, marker = '^', color = 'k', ls = 'None',
         markersize = 8, label = 'AZD4320')

plt.ylabel(r'$C_{pl}$ AZD4320 ($\mu$g/ml)', fontsize=16, labelpad=12)
plt.xlabel(r'Time (h)', fontsize=16, labelpad=12)

plt.grid(linestyle='--', linewidth=0.5)
plt.xscale('linear')
plt.yscale('log')
plt.legend(loc = 0 , numpoints=1, prop={'size': 13.5})
plt.xlim(-0.4,12.4)
plt.ylim(0.008, 1500)
plt.tick_params(labelsize=15)
plt.xticks(np.arange(0.0, 14.0, step=2.0))
plt.tight_layout()
plt.savefig('figure2b_chart.svg',  format='svg')
plt.show()   


# Optimised Parameters (5 mg/kg in mouse)
print('CL = ', CL, ' L/hr')
print('Vc = ', Vc, ' L')
print('k12 = ', k12, ' per hour')
print('k21 = ', k21, ' per hour')
print('k13 = ', k13, ' per hour')
print('k31 = ', k31, ' per hour')