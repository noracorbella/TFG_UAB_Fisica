# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


cell_count_file = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/CC3D_Tutorials/Mitosis/Mitosis_3_1/cell_count.txt"
pressure_file = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/CC3D_Tutorials/Mitosis/Mitosis_3_1/pressure_file.txt"
volume_file = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/CC3D_Tutorials/Mitosis/Mitosis_3_1/volume_file.txt"


cell_count_plot = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/CC3D_Tutorials/Mitosis/Mitosis_3_1/cell_count_plot.png"
pressure_plot = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/CC3D_Tutorials/Mitosis/Mitosis_3_1/pressure_plot.png"
volume_plot = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/CC3D_Tutorials/Mitosis/Mitosis_3_1/volume_plot.png"
combined_plot = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/CC3D_Tutorials/Mitosis/Mitosis_3_1/combined_plot.png"


MCS, cell_count_lst, pressure_lst, volume_lst = [], [], [], []

with open(cell_count_file,'r') as f1, open(pressure_file, 'r') as f2, open(volume_file, 'r') as f3:
    next(f1)
    next(f2)
    next(f3)
    
    for l1, l2, l3 in zip(f1, f2, f3):
        line1, line2, line3 = l1.split(), l2.split(), l3.split()
        step1 = int(line2[0])
        cell_count, pressure, volume = float(line1[1]), float(line2[1]), float(line3[1])
        MCS.append(step1)
        cell_count_lst.append(cell_count)
        pressure_lst.append(pressure)
        volume_lst.append(volume)


MCS = np.array(MCS, dtype=float)
Cell_Count_lst = np.array(cell_count_lst, dtype=float)
Pressure_lst = np.array(pressure_lst)
Volume_lst = np.array(volume_lst)


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




fig, axs = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

axs[0].plot(MCS, Cell_Count_lst, linestyle='-', color='r', alpha=0.5, label = 'Simulation Data')
axs[0].plot(MCS, y_fit, 'b', label = f'Fitted Data, R2 = {R2: .3f}')
axs[0].set_ylabel("Cell Count", fontsize=12)
axs[0].legend()
axs[0].grid(True)

axs[1].plot(MCS, Pressure_lst, linestyle='-', color='r')
axs[1].set_ylabel("Pressure", fontsize=12)
axs[1].grid(True)

axs[2].plot(MCS, Volume_lst, linestyle='-', color='g')
axs[2].set_xlabel("MSC", fontsize=12)
axs[2].set_ylabel("Volume", fontsize=12)
axs[2].grid(True)

plt.tight_layout()
plt.savefig(combined_plot, format='png', dpi=300)
plt.show()


'''
plt.figure(figsize=(8,6))
plt.plot(MSC, Cell_Count_lst, linestyle = '-', color = 'b')
plt.xlabel("MSC", fontsize=12)
plt.ylabel("Cell Count", fontsize=12)
plt.grid(True)
plt.savefig(cell_count_plot, format='png', dpi=300)
plt.show()

plt.figure(figsize=(8,6))
plt.plot(MSC, Pressure_lst, linestyle = '-', color = 'b')
plt.xlabel("MSC", fontsize=12)
plt.ylabel("Volume", fontsize=12)
plt.grid(True)
plt.savefig(pressure_plot, format='png', dpi=300)
plt.show()

plt.figure(figsize=(8,6))
plt.plot(MSC, Volume_lst, linestyle = '-', color = 'b')
plt.xlabel("MSC", fontsize=12)
plt.ylabel("Volume", fontsize=12)
plt.grid(True)
plt.savefig(volume_plot, format='png', dpi=300)
plt.show()
'''
