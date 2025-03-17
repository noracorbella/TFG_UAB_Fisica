# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import matplotlib.pyplot as plt
import random

positions_file = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/TcellsDynabeadSystem_DiffusionTracker/cell_positions.txt"
positions_plot = r"/mnt/c/Users/norac/OneDrive - UAB/Escritorio/uab/5/TFGJordi/ExperimentalData/TcellsDynabeadSystem/TcellsDynabeadSystem_DiffusionTracker/cell_positions.png"



'''
tcell_ids = set()
dynabead_ids = set()

with open(positions_file, "r") as file:
    next(file)
    for line in file:
        _, cell_id, cell_type, _, _ = map(float, line.split())
        if cell_type == 1:
            tcell_ids.add(cell_id)
        elif cell_type == 2:
            dynabead_ids.add(cell_id)

random_tcells = random.sample(list(tcell_ids), 24)
random_dynabeads = random.sample(list(dynabead_ids), 24)

selected_cells = {cell_id: ([], []) for cell_id in random_tcells + random_dynabeads}

with open(positions_file, "r") as file:
    next(file) 
    for line in file:
        mcs, cell_id, cell_type, x, y = map(float, line.split())
        if cell_id in selected_cells:
            selected_cells[cell_id][0].append(x)
            selected_cells[cell_id][1].append(y)

# Step 4: Plot the selected cells' trajectories
plt.figure(figsize=(10, 8))

# Plot the Tcell trajectories (solid lines)
for cell_id in random_tcells:
    plt.plot(selected_cells[cell_id][0], selected_cells[cell_id][1],
             label=f"Tcell {int(cell_id)}", linewidth=1.5, color="blue")

# Plot the Dynabead trajectories (dashed lines)
for cell_id in random_dynabeads:
    plt.plot(selected_cells[cell_id][0], selected_cells[cell_id][1],
             label=f"Dynabead {int(cell_id)}", linewidth=1.5, color="green")

# Add titles, labels, legend, and grid
plt.title("Trajectories of 2 Random Tcells and 2 Random Dynabeads")
plt.xlabel("X Coordinate")
plt.ylabel("Y Coordinate")
#plt.legend()
plt.grid(True)
plt.savefig(positions_plot, format='png', dpi=300)

plt.show()

'''



TCells = {}
Dynabeads = {}



with open(positions_file,'r') as f:
    next(f)
    for l in f:
        mcs, cell_id, cell_type, x, y = map(float, l.split())

        if cell_type == 1:

            if cell_id not in TCells:
                TCells[cell_id] = ([], [])
            TCells[cell_id][0].append(x)
            TCells[cell_id][1].append(y)
                
        elif cell_type == 2:
            if cell_id not in Dynabeads:
                Dynabeads[cell_id] = ([], [])
            Dynabeads[cell_id][0].append(x)
            Dynabeads[cell_id][1].append(y)
            

            
plt.figure(figsize=(10, 8))
for cell_id, (x_coord, y_coord) in TCells.items():
    plt.plot(x_coord, y_coord, color = 'blue', marker = 'o', linewidth = 1, label='TCells' if cell_id == list(TCells.keys())[0] else "")
    
for cell_id, (x_coord, y_coord) in Dynabeads.items():
    plt.plot(x_coord, y_coord, color = 'green', linewidth = 1, label='Dynabeads' if cell_id == list(Dynabeads.keys())[0] else "")

plt.legend()
plt.grid(True)
plt.savefig(positions_plot, format='png', dpi=300)
plt.show()

    


                    

                
