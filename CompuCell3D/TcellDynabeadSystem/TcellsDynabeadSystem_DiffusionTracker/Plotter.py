# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import matplotlib.pyplot as plt
import random
positions_file = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/TcellsDynabeadSystem_DiffusionTracker/cell_positions.txt"
TCell_positions_plot = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/TcellsDynabeadSystem_DiffusionTracker/Tcell_x_vs_mcs.png"
Dynabead_positions_plot = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/TcellsDynabeadSystem_DiffusionTracker/Dynabead_x_vs_mcs.png"

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
    


                    

                
