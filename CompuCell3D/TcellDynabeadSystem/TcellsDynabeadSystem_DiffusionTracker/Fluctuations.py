# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import matplotlib.pyplot as plt
import random
import numpy as np

positions_file = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/TcellsDynabeadSystem_DiffusionTracker/cell_positions.txt"
TCell_fluctuations_hist = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/TcellsDynabeadSystem_DiffusionTracker/Tcell_fluctuations_hist.png"
Dynabead_fluctuations_hist = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/TcellsDynabeadSystem_DiffusionTracker/Dynabead_fluctuations_hist.png"

TCells = {}
Dynabeads = {}

with open(positions_file, 'r') as f:
    next(f)
    for l in f:
        mcs, cell_id, cell_type, x, y = map(float, l.split())
        
        if cell_type == 1:

            if cell_id not in TCells:
                TCells[cell_id] = ([], [], [])
            TCells[cell_id][0].append(mcs)
            TCells[cell_id][1].append(x)
            TCells[cell_id][2].append(y)

                
        elif cell_type == 2:
            if cell_id not in Dynabeads:
                Dynabeads[cell_id] = ([], [], [])
            Dynabeads[cell_id][0].append(mcs)
            Dynabeads[cell_id][1].append(x)
            Dynabeads[cell_id][2].append(y)
  
    
DB_x_fluct = []
DB_y_fluct = []
DB_tot_fluct = []
TC_x_fluct = []
TC_y_fluct = []
TC_tot_fluct = []

for cell_id, (mcs, x_coords, y_coords) in Dynabeads.items():
    for i in range(1, len(x_coords)):
        x_fluct = x_coords[i] - x_coords[i-1]
        y_fluct = y_coords[i] - y_coords[i-1]
        tot_fluct = np.sqrt( x_fluct**2 + y_fluct**2 )
        DB_x_fluct.append(x_fluct)
        DB_y_fluct.append(y_fluct)
        DB_tot_fluct.append(tot_fluct)
        
for cell_id, (mcs, x_coords, y_coords) in TCells.items():
    for i in range(1, len(x_coords)):
        x_fluct = x_coords[i] - x_coords[i-1]
        y_fluct = y_coords[i] - y_coords[i-1]
        tot_fluct = np.sqrt( x_fluct**2 + y_fluct**2 )
        TC_x_fluct.append(x_fluct)
        TC_y_fluct.append(y_fluct)
        TC_tot_fluct.append(tot_fluct)        

    

plt.figure(figsize=(15, 5))

plt.subplot(1, 2, 1)
plt.hist(DB_x_fluct, bins=30, density=True, alpha=0.7, color='blue')
plt.title("Dynabead X Fluctuations")
plt.xlabel("X Fluctuation")
plt.ylabel("Frequency")
plt.xlim(-100, 100)

plt.subplot(1, 2, 2)
plt.hist(DB_y_fluct, bins=30, density=True, alpha=0.7, color='orange')
plt.title("Dynabead Y Fluctuations")
plt.xlabel("Y Fluctuation")
plt.ylabel("Frequency")
plt.xlim(-100, 100)

#plt.subplot(1, 3, 3)
#plt.hist(DB_tot_fluct, bins=30, density=True, alpha=0.7, color='purple')
#plt.title("Dynabead Total Fluctuations")
#plt.xlabel("Total Fluctuation")
#plt.ylabel("Frequency")

plt.tight_layout()
plt.savefig(Dynabead_fluctuations_hist, format='png', dpi=300)
plt.show()



plt.figure(figsize=(15, 5))

plt.subplot(1, 2, 1)
plt.hist(TC_x_fluct, bins=30, density=True, alpha=0.7, color='green')
plt.title("TCell X Fluctuations")
plt.xlabel("X Fluctuation")
plt.ylabel("Frequency")

plt.subplot(1, 2, 2)
plt.hist(TC_y_fluct, bins=30, density=True, alpha=0.7, color='red')
plt.title("TCell Y Fluctuations")
plt.xlabel("Y Fluctuation")
plt.ylabel("Frequency")

#plt.subplot(1, 3, 3)
#plt.hist(TC_tot_fluct, bins=30, density=True, alpha=0.7, color='cyan')
#plt.title("TCell Total Fluctuations")
#plt.xlabel("Total Fluctuation")
#plt.ylabel("Frequency")

plt.tight_layout()
plt.savefig(TCell_fluctuations_hist, format='png', dpi=300)
plt.show()        
  
    
  
    
  
    
            
'''
random_tcell = random.choice(list(TCells.keys()))
random_dynabead = random.choice(list(Dynabeads.keys()))

plt.figure(figsize=(10,5))
#for cell_id, (mcs_list, x_list, y_list) in TCells.items():
#    plt.plot(mcs_list, x_list, marker='o', linewidth=1, label=f'TCell {int(cell_id)}')
    
plt.plot(TCells[random_tcell][0], TCells[random_tcell][1], marker='o', linewidth=1, label=f'TCell {int(random_tcell)}')
plt.xlabel("MCS")
plt.ylabel("X Position")
plt.title("TCell")
plt.legend(loc='best')
plt.grid(True)
plt.savefig(TCell_positions_plot, format='png', dpi=300)
plt.show()

plt.figure(figsize=(10, 5))
#for cell_id, (mcs_list, x_list, y_list) in Dynabeads.items():
#    plt.plot(mcs_list, x_list, marker='o', linewidth=1, label=f'TCell {int(cell_id)}')
plt.plot(Dynabeads[random_dynabead][0], Dynabeads[random_dynabead][1], marker='o', linewidth=1, label=f'Dynabead {int(random_dynabead)}')
plt.xlabel("MCS")
plt.ylabel("X Position")
plt.title("Dynabead")
plt.legend(loc='best')
plt.grid(True)
plt.savefig(Dynabead_positions_plot, format='png', dpi=300)
plt.show()   
    
'''

                    

                
