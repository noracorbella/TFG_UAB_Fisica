#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 23 19:30:20 2025

@author: nora
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)


cell_count_file = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/CC3D_Tutorials/Mitosis/Mitosis_3_1/cell_count.txt"
cell_count_png = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/CC3D_Tutorials/Mitosis/Mitosis_3_1/cell_count_fit.png"

MCS, cell_count_lst = [], []

with open(cell_count_file,'r') as f1:
    next(f1)
        
    for l1 in f1:
        line1 = l1.split()
        step1 = int(line1[0])
        cell_count = float(line1[1])
        MCS.append(step1)
        cell_count_lst.append(cell_count)
        
MCS = np.array(MCS, dtype=float)
Cell_Count_lst = np.array(cell_count_lst, dtype=float)


def sigmoid(x, a, b, c):
    return a / (1 + np.exp(-b * (x - c)))


p0=(max(Cell_Count_lst), 1, np.median(MCS))
bounds = ([0, 0, -np.inf], [np.inf, 1, np.inf])

popt, pcov = curve_fit(sigmoid, MCS, Cell_Count_lst, p0=p0)

a_fit, b_fit, c_fit = popt
print(f"Fitted parameters: a = {a_fit}, b = {b_fit}, c = {c_fit}")

residuals = Cell_Count_lst - sigmoid(MCS, *popt)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((Cell_Count_lst - np.mean(Cell_Count_lst))**2)
R2 = 1 - (ss_res / ss_tot)
print(f"$R**2$ value: {R2}")


y_fit = sigmoid(MCS, *popt)

plt.plot(MCS, Cell_Count_lst, label='Simulation Data', alpha=0.5)
plt.plot(MCS, y_fit, 'r-', label=f'Fit: {a_fit:.2f} / (1 + exp(-{b_fit:.2f} * (x - {c_fit:.2f})))')
plt.xlabel('MCS')
plt.ylabel('Cell Count')
plt.legend()
plt.grid(True)
plt.savefig(cell_count_png, format='png', dpi=300)
plt.show()

'''
Fitted parameters: a = 1384.8419624513772, b = 0.002973253394504467, c = 1636.917299577077
$R**2$ value: 0.9980982809767125
'''