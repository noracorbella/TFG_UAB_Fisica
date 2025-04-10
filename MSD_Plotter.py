#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 14:50:21 2025

@author: nora
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

file = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/TcellsDynabeadSystem_DBTracker/TcellsDynabeadSystem_DBTracker_1/MDAnalysis_MSD_DB.dat"
plot = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/TcellsDynabeadSystem_DBTracker/TcellsDynabeadSystem_DBTracker_1/MDAnalysis_MSD_DB.png"


data = np.loadtxt(file)

time_window = data[:, 0]
MSD = data[:, 1]

def linear_function(x, m, n):
    return m * x + n

params, _ = curve_fit(linear_function, time_window, MSD)
m, n = params

MSD_fit = linear_function(time_window, m, n)

plt.figure(figsize=(9, 6))
plt.grid()
plt.scatter(time_window, MSD, color='blue', s=10, label='Simulation Data')
plt.plot(time_window, MSD_fit, color='red', linewidth=2, label='Linear fit')
plt.xlabel('Time Window')
plt.ylabel('MSD')
plt.legend()
plt.savefig(plot, format='png', dpi=300)
plt.show()

