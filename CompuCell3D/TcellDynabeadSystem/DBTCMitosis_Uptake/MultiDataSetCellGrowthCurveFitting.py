#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 31 14:52:25 2025

@author: nora
"""


import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

datasets = [
    {
        'file': r"TC_count_4.txt",
        'label': 'Simulation 4',
        'color': 'tab:blue',
        'line_style': '-'
    },
    {
        'file': r"TC_count_5.txt",
        'label': 'Simulation 5',
        'color': 'tab:orange',
        'line_style': '-'
    },
    {
        'file': r"TC_count_6.txt",
        'label': 'Simulation 6',
        'color': 'tab:green',
        'line_style': '-'
    },
    {
        'file': r"TC_count_7.txt",        
        'label': 'Simulation 7',
        'color': 'tab:red',
        'line_style': '-'
    }]

# Output plot
output_png = r"/TC_count_plot_4_5_6_7_combined.png"

# Growth Functions
def sigmoid(x, a, b, c, d):
    """Modified sigmoid with baseline offset"""
    return d + (a - d) / (1 + np.exp(-b * (x - c)))

def exponential_saturation(x, a, b, c):
    """Exponential approach to saturation"""
    return a * (1 - np.exp(-b * x)) + c

def logistic_growth(x, K, r, t0, N0):
    """Standard logistic growth equation"""
    return K / (1 + ((K - N0) / N0) * np.exp(-r * (x - t0)))

def load_and_process_data(file_path, crop_length=110000):
    MCS, cell_count_lst = [], []
    
    try:
        with open(file_path, 'r') as f1:
            next(f1)  # Skip header
            
            for l1 in f1:
                line1 = l1.split()
                step1 = int(line1[0])
                cell_count = float(line1[3])
                MCS.append(step1)
                cell_count_lst.append(cell_count)
        
        MCS = np.array(MCS, dtype=float)
        Cell_Count_lst = np.array(cell_count_lst, dtype=float)
        
        # Crop data if specified
        if crop_length and len(Cell_Count_lst) > crop_length:
            Cell_Count_cropped = Cell_Count_lst[0:crop_length]
            MCS_cropped = MCS[0:crop_length]
        else:
            Cell_Count_cropped = Cell_Count_lst
            MCS_cropped = MCS
            
        return MCS_cropped, Cell_Count_cropped
        
    except FileNotFoundError:
        print(f"Warning: File {file_path} not found. Skipping this dataset.")
        return None, None
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return None, None

def fit_best_model(MCS_cropped, Cell_Count_cropped):
    """Find the best fitting model for the data"""
    
    # Analyze the data to get better initial parameters
    y_min = np.min(Cell_Count_cropped)
    y_max = np.max(Cell_Count_cropped)
    x_min = np.min(MCS_cropped)
    x_max = np.max(MCS_cropped)
    
    # Find approximate inflection point
    midpoint_y = (y_min + y_max) / 2
    midpoint_idx = np.argmin(np.abs(Cell_Count_cropped - midpoint_y))
    inflection_x = MCS_cropped[midpoint_idx] if midpoint_idx < len(MCS_cropped) else np.median(MCS_cropped)
    
    # Define models to try
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
            
            if r2 > best_r2:
                best_r2 = r2
                best_model = model_info['func']
                best_params = popt
                best_name = name
                
        except Exception as e:
            print(f"{name} failed: {e}")
            continue
    
    return best_model, best_params, best_r2, best_name

plt.figure(figsize=(10, 6))

results = []
valid_datasets = []

for i, dataset in enumerate(datasets):
    print(f"\nProcessing {dataset['label']}...")
    
    MCS_cropped, Cell_Count_cropped = load_and_process_data(dataset['file'])
    
    if MCS_cropped is None or Cell_Count_cropped is None:
        continue
    
    print(f"Data points: {len(Cell_Count_cropped)}")
    
    # Fit best model
    best_model, best_params, best_r2, best_name = fit_best_model(MCS_cropped, Cell_Count_cropped)
    
    if best_model is not None:
        print(f"Best model: {best_name}")
        print(f"R² = {best_r2:.4f}")
        
        results.append({
            'MCS': MCS_cropped,
            'Cell_Count': Cell_Count_cropped,
            'model': best_model,
            'params': best_params,
            'r2': best_r2,
            'name': best_name,
            'dataset': dataset
        })
        valid_datasets.append(dataset)
        
        # Generate smooth curve for plotting
        x_smooth = np.linspace(np.min(MCS_cropped), np.max(MCS_cropped), 1000)
        y_smooth = best_model(x_smooth, *best_params)
        
        # Priginal data
        plt.plot(MCS_cropped, Cell_Count_cropped, 
                color=dataset['color'], linestyle='--', 
                alpha=0.7, linewidth=1.5, 
                label=f'{dataset["label"]}')
        
        # Fitted curve
        plt.plot(x_smooth, y_smooth, 
                color=dataset['color'], linestyle='-', 
                linewidth=2, alpha=0.9,
                label=f'Sigmoid Fit (R²={best_r2:.3f})')
    else:
        print(f"No suitable model found for {dataset['label']}")
        # Plot the data
        plt.plot(MCS_cropped, Cell_Count_cropped, 
                color=dataset['color'], linestyle=dataset['line_style'], 
                alpha=0.7, linewidth=1.5, 
                label=f'{dataset["label"]} (No Fit)')

plt.xlabel('MCS', fontsize=14)
plt.ylabel('Cell Count', fontsize=14)
plt.legend(loc='upper left', fontsize=12)
plt.grid(True, alpha=0.3)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)

plt.tight_layout()
plt.savefig(output_png, format='png', dpi=300, bbox_inches='tight')
plt.show()

# Print results
print(f"\n{'='*60}")
print("SUMMARY OF RESULTS")
print(f"{'='*60}")

for i, result in enumerate(results):
    print(f"\n{result['dataset']['label']}:")
    print(f"  Best model: {result['name']}")
    print(f"  R²: {result['r2']:.4f}")
    print(f"  Parameters: {result['params']}")
    
    if 'Logistic' in result['name']:
        print(f"  Carrying capacity (K) = {result['params'][0]:.2f}")
        print(f"  Growth rate (r) = {result['params'][1]:.2e}")
        print(f"  Initial population (N0) = {result['params'][3]:.2f}")
