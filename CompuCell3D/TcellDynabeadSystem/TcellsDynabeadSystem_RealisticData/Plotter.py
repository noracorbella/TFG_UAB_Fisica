# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


diffusion_file = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/TcellsDynabeadSystem_ODE/Diffusion_Data.txt"


diffusion_plot = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/TcellsDynabeadSystem_ODE/Diffusion_plot_ODE_RealisticData.png"



mcs, D_DB_lst, D_TC_lst = [], [], []

with open(diffusion_file,'r', encoding='utf-8') as f:
    next(f)
    for l in f:
        line = l.strip().split()
        if not line:  # Skip empty lines
            continue
        try:
            s, D_DB, D_TC = float(line[0]), float(line[1]), float(line[2])
            D_DB_lst.append(D_DB)
            D_TC_lst.append(D_TC)
            mcs.append(s)
        except (ValueError, IndexError) as e:
            print(f"Error processing line: {l}. Error: {e}")
  
    
MCS = np.array(mcs, dtype = int)
D_Dynabead = np.array(D_DB_lst, dtype=float)
D_TCell = np.array(D_TC_lst, dtype=float)




plt.figure(figsize=(8,6))
plt.plot(MCS, D_Dynabead, linestyle = '-', color = 'g', label = 'Dynabead')
plt.plot(MCS, D_TCell, linestyle = '-', color = 'r', label = 'TCell')
plt.xlabel("Step", fontsize=12)
plt.ylabel("Diffusion", fontsize=12)
plt.grid(True)
plt.legend()
plt.savefig(diffusion_plot, format='png', dpi=300)
plt.show()


