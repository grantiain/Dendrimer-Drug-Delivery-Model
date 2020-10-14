#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


def conc_tu(t_tu, t_hy, t):
    # Concentration of Released API in tumour
    k_res = np.log(2) / t_res
    k_tu = np.log(2) / t_tu
    k_h = np.log(2) / t_hy   
    nm = k_res * np.exp(-k_tu * t)+(k_h - k_tu) * np.exp(-(k_h + k_res) * t) \
         + (k_tu - k_h - k_res) * np.exp(-k_h * t)
    dnm = (k_h - k_tu + k_res) * (k_h - k_tu)
    return ((k_h * nm/dnm) * (dose * k_ext / (k_res * V_tu)))
        
        
def dcdt(c, t, t_hy, t_res, CL, Vc, k12, k21, k13, k31):
    # PK Model for Plasma Concentrations
    k_res = np.log(2)/t_res
    k_h = np.log(2)/t_hy
    k_el = CL/Vc
    dc0 = -(k_res+k_h)*c[0]
    dc1 = k_h*c[0] - (k_el+k12+k13)*c[1] + k21*c[2] + k31*c[3]
    dc2 = k12*c[1]-k21*c[2]
    dc3 = k13*c[1]-k31*c[3]
    return [dc0,dc1,dc2,dc3]


# PK Parameter for API
CL = 0.03
Vc = 0.0038925
k12 = 1.90724
k21 = 0.0262561
k13 = 1.00656
k31 = 1.04364

CL =  0.029618110766202254  
Vc =  0.004228781280327043  
k12 =  2.306847594088607  
k21 =  0.10833051046041449 
k13 =  3.9181302136597678  
k31 =  4.35611563799507  

t_res = 3.301
k_res = np.log(2) / t_res
t_tu = 13.1

# 1 mg/kg in a 25 g mouse
dose = 0.025

V_tu = 3.0e-4
k_ext = 0.00038

x0 = [dose, 0.0, 0.0]
c0 = [dose, 0.0, 0.0, 0.0]
ts = np.linspace(0,200,600)     
t = np.linspace(0,200,600)
t_hyd = np.linspace(1,49,600)
fn_val_max1 = []
fn_val_max2 = []
c_pl_val_max = []

# Find maximum concentrations of released API in plasma over a range
# of hydrolysis half-lives from 1 to 49 hours

for i in range(0,600):
    fn_val = conc_tu(t_tu,t_hyd[i], t)
    fn_val_max1.append(max(fn_val))
    cs = odeint(dcdt, c0, ts, args=(t_hyd[i], t_res, CL, Vc, k12, k21, k13, k31)) 
    c_rel_pl = cs[:,1]/Vc
    c_pl_val_max.append(max(c_rel_pl))

# Find maximum concentrations of released API in tumour tissue over a range
# of hydrolysis half-lives from 1 to 49 hours
    
for i in range(0,600):
    fn_val = conc_tu(t_tu,t_hyd[i], t)
    fn_val_max2.append(max(fn_val))

       
plt.plot(t_hyd, fn_val_max1, label = r'$k_{res}$ = 0.23 $h^{-1}$')
plt.ylabel(r'$\bar{C}_{max,tu}$ ($\mu$g/ml/(mg/kg))', fontsize=16, labelpad=12)
plt.xlabel(r'Half-Life (h)', fontsize=16, labelpad=12)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.grid(b=True)
plt.ylim(ymin=0, ymax=0.061)
plt.grid(linestyle='--', linewidth=0.5)
plt.tight_layout()
# Save graphics file
plt.savefig('Figure_5D.svg',  format='svg')
plt.show()

plt.plot(t_hyd, c_pl_val_max, color = 'red')
plt.ylabel(r'$\bar{C}_{max,pl}$ ($\mu$g/ml/(mg/kg))', fontsize=16, labelpad=12)
plt.xlabel(r'Half-Life (h)', fontsize=16, labelpad=12)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.grid(b=True)
plt.grid(linestyle='--', linewidth=0.5)
plt.ylim(ymin=0)
plt.tight_layout()
plt.savefig('Figure_5E.svg',  format='svg')
plt.show()

# Calculate Optimisation Index Values
ti_fn = np.array(fn_val_max1)**2/(np.array(c_pl_val_max))
plt.plot(t_hyd, ti_fn, color = 'green')
plt.ylabel(r'$C_{max,tu}^2\, /\, C_{max,pl}$ ($\mu$g/ml/(mg/kg))', 
           fontsize=13, labelpad=10)
plt.xlabel(r'Half-Life (h)', fontsize=16, labelpad=12)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.grid(b=True)
plt.ylim(ymin=0, ymax=0.1)
plt.grid(linestyle='--', linewidth=0.5)
plt.title('  ')
plt.tight_layout()
# Save graphics file
plt.savefig('Figure_5F.svg',  format='svg')
plt.show()