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

cell_count_file = r"TC_count.txt"
cell_count_png = r"TC_count_plot.png"

MCS, cell_count_lst = [], []
with open(cell_count_file,'r') as f1:
    next(f1)
    
    for l1 in f1:
        line1 = l1.split()
        step1 = int(line1[0])
        cell_count = float(line1[3])
        MCS.append(step1)
        cell_count_lst.append(cell_count)

MCS = np.array(MCS, dtype=float)
Cell_Count_lst = np.array(cell_count_lst, dtype=float)
print(len(Cell_Count_lst))

Cell_Count_cropped = Cell_Count_lst[0:110000]
MCS_cropped = MCS[0:110000]
print(len(Cell_Count_cropped))
print(type(MCS), type(Cell_Count_cropped))

# For step-like mitosis data, try multiple models
def sigmoid(x, a, b, c, d):
    """Modified sigmoid with baseline offset"""
    return d + (a - d) / (1 + np.exp(-b * (x - c)))

def exponential_saturation(x, a, b, c):
    """Exponential approach to saturation"""
    return a * (1 - np.exp(-b * x)) + c

def logistic_growth(x, K, r, t0, N0):
    """Standard logistic growth equation"""
    return K / (1 + ((K - N0) / N0) * np.exp(-r * (x - t0)))

# Analyze the data to get better initial parameters
y_min = np.min(Cell_Count_cropped)
y_max = np.max(Cell_Count_cropped)
x_min = np.min(MCS_cropped)
x_max = np.max(MCS_cropped)

# Find approximate inflection point
midpoint_y = (y_min + y_max) / 2
midpoint_idx = np.argmin(np.abs(Cell_Count_cropped - midpoint_y))
inflection_x = MCS_cropped[midpoint_idx] if midpoint_idx < len(MCS_cropped) else np.median(MCS_cropped)

print(f"Data analysis:")
print(f"Y range: {y_min:.2f} to {y_max:.2f}")
print(f"X range: {x_min:.0f} to {x_max:.0f}")
print(f"Estimated inflection point: {inflection_x:.0f}")

# Try model and choose the best one
models = {
    'Modified Sigmoid': {
        'func': sigmoid,
        'p0': [y_max - y_min, 1e-5, inflection_x, y_min],
        'bounds': ([0, 1e-8, x_min, 0], [2*(y_max-y_min), 1e-3, x_max, y_max])
    },
    'Exponential Saturation': {
        'func': exponential_saturation,
        'p0': [y_max - y_min, 1e-5, y_min],
        'bounds': ([0, 1e-8, 0], [2*(y_max-y_min), 1e-3, y_max])
    },
    'Logistic Growth': {
        'func': logistic_growth,
        'p0': [y_max, 1e-5, inflection_x, y_min],
        'bounds': ([y_max*0.5, 1e-8, x_min, 0], [y_max*2, 1e-3, x_max, y_max*0.5])
    }
}

best_model = None
best_r2 = -np.inf
best_params = None
best_name = ""

for name, model_info in models.items():
    print(f"\nTrying {name}...")
    
    try:
        # Try with TRF method first
        popt, pcov = curve_fit(
            model_info['func'], 
            MCS_cropped, 
            Cell_Count_cropped, 
            p0=model_info['p0'],
            bounds=model_info['bounds'],
            maxfev=15000,
            method='trf'
        )
        
        # Calculate R²
        y_pred = model_info['func'](MCS_cropped, *popt)
        ss_res = np.sum((Cell_Count_cropped - y_pred) ** 2)
        ss_tot = np.sum((Cell_Count_cropped - np.mean(Cell_Count_cropped)) ** 2)
        r2 = 1 - (ss_res / ss_tot)
        
        print(f"{name} - R² = {r2:.4f}")
        print(f"Parameters: {popt}")
        
        if r2 > best_r2:
            best_r2 = r2
            best_model = model_info['func']
            best_params = popt
            best_name = name
            
    except Exception as e:
        print(f"{name} failed: {e}")
        continue

if best_model is not None:
    print(f"\nBest model: {best_name}")
    print(f"Best R² = {best_r2:.3f}")
    
    # Generate curve for plotting
    x_smooth = np.linspace(np.min(MCS_cropped), np.max(MCS_cropped), 1000)
    y_smooth = best_model(x_smooth, *best_params)
    y_fit = best_model(MCS_cropped, *best_params)
    
    # Plot results
    plt.figure(figsize=(10, 6))
    
    # Original data 
    plt.plot(MCS_cropped, Cell_Count_cropped, 'b-', alpha=0.8, linewidth=0.8, label='Simulation Data')
    
    # Fitted curve
    plt.plot(x_smooth, y_smooth, 'r-', linewidth=2, label=f'Sigmoid Fit ($R²={best_r2:.3f}$)')
    
    plt.xlabel('MCS', fontsize=14)
    plt.ylabel('Cell Count', fontsize=14)
    #plt.title('Cell Growth Curve Fitting', fontsize=16)
    plt.legend(loc='upper left', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    
    #plt.text(0.6, 0.15, f'Model: {best_name}\nR² = {best_r2:.4f}', 
    #         transform=plt.gca().transAxes, fontsize=10,
    #         bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(cell_count_png, format='png', dpi=300, bbox_inches='tight')
    plt.show()
    
           
else:
    print("No model found.")
    
    # Plot the data for visualization
    plt.figure(figsize=(10, 6))
    plt.plot(MCS_cropped, Cell_Count_cropped, 'b-', alpha=0.8, linewidth=1, label='Simulation 1')
    plt.xlabel('MCS', fontsize=14)
    plt.ylabel('Cell Count', fontsize=14)
    #plt.title('Cell Growth Data (No Suitable Model Found)', fontsize=16)
    plt.legend(loc='upper left', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.savefig(cell_count_png, format='png', dpi=300, bbox_inches='tight')
    plt.show()
