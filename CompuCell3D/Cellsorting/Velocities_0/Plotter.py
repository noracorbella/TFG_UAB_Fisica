# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import numpy as np
import matplotlib.pyplot as plt

E_contactareabytype_log = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/CC3D_Tutorials/CellSorting/Velocities_0/contact_areabytype_data.txt"
E_contactareawithmedium_log = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/CC3D_Tutorials/CellSorting/Velocities_0/contact_areabymedium_data.txt"
velocities_log = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/CC3D_Tutorials/CellSorting/Velocities_0/velocity_file.txt"


E_contactareabytype_plot = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/CC3D_Tutorials/CellSorting/Velocities_0/contact_areabytype_data.png"
E_contactareawithmedium_plot = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/CC3D_Tutorials/CellSorting/Velocities_0/contact_areawithmedium_data.png"
velocities_plot = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/CC3D_Tutorials/CellSorting/Velocities_0/velocities_plot.png"

MCS, ACC_lst, ACN_lst, ANN_lst = [], [], [], []
ACM_lst, ANM_lst = [], []
vel_lst = []

with open(E_contactareabytype_log,'r') as f1, open(E_contactareawithmedium_log, 'r') as f2, open(velocities_log) as f3:
    next(f1)
    next(f2)
    next(f3)
    
    for l1, l2, l3 in zip(f1, f2, f3):
        line1 = l1.split()
        step1 = int(line1[0])
        acc, acn, ann = float(line1[1]), float(line1[2]), float(line1[3])
        line2 = l2.split()
        acm, anm = float(line2[1]), float(line2[2])
        line3 = l3.split()
        vel = float(line3[1])
        
        
        if not step1%10:
        
            MCS.append(step1)
            ACC_lst.append(acc)
            ACN_lst.append(acn)
            ANN_lst.append(ann)
            ACM_lst.append(acm)
            ANM_lst.append(anm)
            vel_lst.append(vel)



MSC = np.array(MCS)
ACC = np.array(ACC_lst)
ACN = np.array(ACN_lst)
ANN = np.array(ANN_lst)
ACM = np.array(ACM_lst)
ANM = np.array(ANM_lst)
Velocities = np.array(vel_lst)



plt.figure(figsize=(8,6))
plt.plot(MSC, ACC, linestyle = '-', color = 'b', label='ACC')
plt.plot(MSC, ACN, linestyle = '-', color = 'r', label='ACN')
plt.plot(MSC, ANN, linestyle = '-', color = 'g', label='ANN')
plt.xlabel("MSC", fontsize=12)
plt.ylabel("Contact Area by Cell Type", fontsize=12)
plt.legend()
plt.grid(True)
plt.savefig(E_contactareabytype_plot, format='png', dpi=300)
plt.show()

plt.figure(figsize=(8,6))
plt.plot(MSC, ACM, linestyle = '-', color = 'b', label='ACM')
plt.plot(MSC, ACC, linestyle = '-', color = 'g', label='ANM')
plt.xlabel("MSC", fontsize=12)
plt.ylabel("Contact Area with Medium", fontsize=12)
plt.legend()
plt.grid(True)
plt.savefig(E_contactareawithmedium_plot, format='png', dpi=300)
plt.show()

plt.figure(figsize=(8,6))
plt.plot(MSC, Velocities, linestyle = '-', color = 'b')
plt.xlabel("MSC", fontsize=12)
plt.ylabel("Velocity", fontsize=12)
plt.grid(True)
plt.savefig(velocities_plot, format='png', dpi=300)
plt.show()