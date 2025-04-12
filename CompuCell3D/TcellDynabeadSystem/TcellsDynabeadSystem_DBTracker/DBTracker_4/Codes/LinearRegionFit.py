#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 15:32:23 2025

@author: nora
"""

import numpy as np
from scipy.stats import linregress


MSD_data = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/DBTracker/DBTracker_5/MDAnalysis_MSD_DB5.dat"

debug_file = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/DBTracker/DBTracker_5/DB5_R2_data.txt"




plot = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/DBTracker/DBTracker_5/MSD_DB5_LinearFit.png"

output_data = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/DBTracker/DBTracker_5/MSD_DB5_LinearFit_data.txt"


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def linear_model(x, a, b):
    return a * x + b

def linear_fit(filename, debug_file, col_x=0, col_y=1, delimiter=None, min_points=10000):
    # Load data from file
    data = np.loadtxt(filename, delimiter=delimiter)
    x = data[:, col_x]
    y = data[:, col_y]
    
    best_start, best_end, best_r2 = 0, min_points, -1
    
    with open(debug_file, "w") as debug:
        debug.write("Start Index \t End Index \t rsquared \n")
        for start in range(len(x) - min_points + 1):
            for end in range(start + min_points, len(x)):
                x_subset, y_subset = x[start:end], y[start:end]
                popt, _ = curve_fit(linear_model, x_subset, y_subset)
                slope, intercept = popt
                
                y_fit = linear_model(x_subset, slope, intercept)
                residuals = y_subset - y_fit
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((y_subset - np.mean(y_subset))**2)
                r_squared = 1 - (ss_res / ss_tot)
                
                print(f"Range {start} to {end}: R² = {r_squared:.4f}")
                debug.write(f"{start} \t {end} \t {r_squared:.4f}\n")
                
                if r_squared > best_r2:
                    best_start, best_end, best_r2 = start, end, r_squared
                

    
    # Compute final fit
    x_final, y_final = x[best_start:best_end], y[best_start:best_end]
    popt, _ = curve_fit(linear_model, x_final, y_final)
    slope, intercept = popt
    y_fit = linear_model(x_final, slope, intercept)
    
    # Plot data and fit
    plt.scatter(x, y, label='Data')
    plt.plot(x_final, y_fit, color='red', label=f'Fit: y={slope:.2f}x+{intercept:.2f}')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend()
    plt.title(f'Linear Fit (R²={best_r2:.4f})')
    plt.savefig(plot, format='png', dpi=300)
    plt.show()
    
    with open(output_data, "w") as summaryfile:
        summaryfile.write(f'Slope: {slope} \n')
        summaryfile.write(f'Slope: {slope}\n')
        summaryfile.write(f'Intercept: {intercept}\n')
        summaryfile.write(f'R²: {best_r2}\n')
        summaryfile.write(f'Linear region from index {best_start} to {best_end}\n')
    
    return best_start, best_end, slope, intercept, best_r2

linear_fit(MSD_data, debug_file)