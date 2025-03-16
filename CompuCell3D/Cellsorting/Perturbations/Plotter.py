# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import numpy as np
import matplotlib.pyplot as plt

E_log = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/CC3D_Tutorials/CellSorting/Perturbations/Perturbation_0_1/Perturbation_0_1_cc3d_02_04_2025_11_46_14_e220f1/simulation.log"

out_file = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/CC3D_Tutorials/CellSorting/erturbation_0_1/Perturbation_0_1_cc3d_02_04_2025_11_46_14_e220f1/energy_plot_not_isolated.png"


MCS, E = [], []


with open(E_log,'r') as f1:
    for l1 in f1:
        if l1.startswith('INFO: Step'):
            line1 = l1.split()
            step1 = int(line1[2])
            energy1 = float(line1[6])
            MCS.append(step1)
            E.append(energy1)
            
MSC = np.array(MCS)
Energy1 = np.array(E)



plt.figure(figsize=(8,6))
plt.plot(MSC, Energy1, linestyle = '-', color = 'r')
plt.xlabel("MSC", fontsize=12)
plt.ylabel("Energy", fontsize=12)
plt.grid(True)
plt.savefig(out_file, format='png', dpi=300)
plt.show()